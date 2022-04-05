!
! module m_likelihood
! likelihood function here
!
module m_body_likelihood

    use m_exception, only : exception_raiseWarning, exception_raiseError
    use m_logger, only    : log_msg
    use m_utils, only     : ii10, itoa, rtoa
    use m_settings, only  : T_GRID, T_MOD, out_bnd, mod_setup
    use like_settings, only: T_LIKE_BASE, T_LIKE_SET, T_DATA
    use run_info, only    : T_RUN_INFO
    use fastmarching, only: fastmarching3d, T_RAY, gradientdescent
    use cgal_delaunay,only: d3
    use kdtree2_precision_module, only : kdkind
    use kdtree2_module, only    : kdtree2, kdtree2_result, kdtree2_create,&
                                  kdtree2_n_nearest, kdtree2_destroy

    use omp_lib
    use iso_c_binding

    implicit none
    private

    public :: body_likelihood, body_noise_likelihood, body_location_likelihood
    public :: body_likelihood_grads
    public :: point2idx

    ! debug
    real( kind=ii10 ) :: t1, t2
    ! static value
    real( kind=ii10 ), parameter :: ERAD = 6371
    real( kind=ii10 ), parameter :: PI2 = 6.283185
    real( kind=ii10 ), parameter :: EPS = tiny(0d0)*10


contains

    subroutine body_noise_likelihood(dat,RTI,settings,like)
        implicit none
        type(T_DATA), intent(in)                 :: dat
        type(T_RUN_INFO), intent(in)             :: RTI
        type(T_LIKE_SET), intent(in)             :: settings
        type(T_LIKE_BASE), intent(inout)         :: like

        integer nrr
        integer nsrc, nrev
        integer i, j, k

        nsrc = dat%nsrc
        nrev = dat%nrev

        do i = 1, dat%np
            like%sigma(:,i) = RTI%bnoise0(i)*like%srdist(:,i) + RTI%bnoise1(i)
        enddo

        like%like = 0
        like%misfit = 0
        like%unweighted_misfit = 0
        do i = 1, dat%np
            nrr  = 0
            do j = 1, nsrc
                   do k = 1, nrev
                    nrr = nrr + 1
                    if(dat%raystat(nrr,1,i) == 1) then
                        like%like = like%like + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
                            -dat%ttime(nrr,1,i))**2/( 2*(like%sigma(nrr,i))**2 )
                        like%misfit = like%misfit + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
                            -dat%ttime(nrr,1,i))**2/(like%sigma(nrr,i)**2)
                        like%unweighted_misfit = like%unweighted_misfit + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
                            -dat%ttime(nrr,1,i))**2
                    else
                        like%sigma(nrr,i) = 1.0
                    endif
                enddo
            enddo
        enddo

        like%like = like%like + sum(log(like%sigma)) + dat%nrays/2.0 * log(PI2)
        return
    end subroutine body_noise_likelihood

    subroutine body_location_likelihood(dat,model,RTI,perturbed_box,settings,like)
        implicit none
        type(T_DATA), intent(in)                    :: dat
        type(T_MOD), intent(in)                     :: model
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(d3), dimension(2), intent(in)          :: perturbed_box
        type(T_LIKE_SET), intent(in)                :: settings
        type(T_LIKE_BASE), intent(inout)            :: like

        integer nrr
        integer nsrc, nrev
        integer i, j, k

        if(settings%dynamic == 0)then
            call body_likelihood(dat,model,RTI,perturbed_box,settings,like)
            return
        endif

        nsrc = dat%nsrc
        nrev = dat%nrev

        ! interpolate the source locations
        do k = 1, dat%np
            do i = 1, nsrc
                do j = 1, nrev
                    like%phaseTime(j,i,k)=Interpolate3d(like%field5d(:,:,:,j,k),&
                        settings%grid,RTI%locations(1:3,i))
                enddo
            enddo
        enddo

        ! using travel time instead of source-receiver distance
        do j = 1, nsrc
            do k = 1, nrev
                nrr = k + (j-1)*nrev
                like%srdist(nrr,:) = like%phaseTime(k,j,:)
            enddo
        enddo

        do i = 1, dat%np
            like%sigma(:,i) = RTI%bnoise0(i)*like%srdist(:,i) + RTI%bnoise1(i)
        enddo

        like%like = 0
        like%misfit = 0
        like%unweighted_misfit = 0
        do i = 1, dat%np
            nrr  = 0
            do j = 1, nsrc
                   do k = 1, nrev
                    nrr = nrr + 1
                    if(dat%raystat(nrr,1,i) == 1) then
                        like%like = like%like + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
                            -dat%ttime(nrr,1,i))**2/( 2*(like%sigma(nrr,i))**2 )
                        like%misfit = like%misfit + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
                            -dat%ttime(nrr,1,i))**2/(like%sigma(nrr,i)**2)
                        like%unweighted_misfit = like%unweighted_misfit + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
                            -dat%ttime(nrr,1,i))**2
                    else
                        like%sigma(nrr,i) = 1.0
                    endif
                enddo
            enddo
        enddo

        like%like = like%like + sum(log(like%sigma)) + dat%nrays/2.0 * log(PI2)
        return
    end subroutine body_location_likelihood

    subroutine body_likelihood(dat,model,RTI,perturbed_box,settings,like)
    
        implicit none
        type(T_DATA), intent(in)                    :: dat
        type(T_MOD), intent(in)                     :: model
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(d3), dimension(2), intent(in)          :: perturbed_box  ! TODO: for dynamic modelling
        type(T_LIKE_SET), intent(in)                :: settings
        type(T_LIKE_BASE), intent(inout)            :: like

        ! > local variable for travel times and rays
        integer(c_int) nsrc, nrev
        !integer, dimension(:), allocatable :: crazyray

        ! > local grid
        type(T_GRID) :: grid

        integer nrr
        integer i, j, k

        nsrc = dat%nsrc
        nrev = dat%nrev
        grid = settings%grid

        ! noise level
        if(settings%sigdep /= 0)then
            do i = 1, dat%np
                like%sigma(:,i) = RTI%bnoise0(i)*like%srdist(:,i) + RTI%bnoise1(i)
            enddo
        else
            like%sigma = dat%ttime(:,2,:)
        endif
        ! set no data sigma to 1, this is important since sigma is log added to
        ! lglikelihood
        do i = 1, dat%np
            nrr  = 0
            do j = 1, nsrc
                   do k = 1, nrev
                    nrr = nrr + 1
                    if(dat%raystat(nrr,1,i) /= 1) then
                        like%sigma(nrr,i) = 1.0
                    endif
                enddo
            enddo
        enddo

        if(settings%dynamic==0)then

#ifdef _OPENMP
            !t1 = omp_get_wtime ( )
            !call omp_set_num_threads(settings%nthreads)
#endif
            if(nsrc<=nrev)then
                !$omp parallel
                !$omp do private(i,j)
                do i = 1, nsrc
                    if (settings%datatype==0)then
                        call fastmarching3d(model%vp,grid%nx,grid%ny,grid%nz,& 
                            grid%dx, grid%dy, grid%dz/grid%scaling,&
                            grid%xmin,grid%ymin,grid%zmin/grid%scaling,rti%locations(1:3,i),&
                            nrev, dat%rev,settings%order,like%phaseTime(:,i,1))
                    else
                        call fastmarching3d(model%vp,grid%nx,grid%ny,grid%nz,& 
                            grid%dx, grid%dy, grid%dz/grid%scaling,&
                            grid%xmin,grid%ymin,grid%zmin/grid%scaling,rti%locations(1:3,i),&
                            nrev, dat%rev,settings%order,like%phaseTime(:,i,1))
                        call fastmarching3d(model%vs,grid%nx,grid%ny,grid%nz,& 
                            grid%dx, grid%dy, grid%dz/grid%scaling,&
                            grid%xmin,grid%ymin,grid%zmin/grid%scaling,rti%locations(1:3,i),&
                            nrev, dat%rev,settings%order,like%phaseTime(:,i,2))
                    endif
                enddo
                !$omp end do
                !$omp end parallel
                !t2 = omp_get_wtime ( )
            else
                ! receivers is less than sources, tracing from receivers
                !$omp parallel
                !$omp do private(i,j)
                do i = 1, nrev
                    if (settings%datatype==0)then
                        call fastmarching3d(model%vp,grid%nx,grid%ny,grid%nz,& 
                            grid%dx, grid%dy, grid%dz/grid%scaling,&
                            grid%xmin,grid%ymin,grid%zmin/grid%scaling,dat%rev(:,i),&
                            nsrc,rti%locations(1:3,:),settings%order,like%phaseTime(i,:,1))
                    else
                        call fastmarching3d(model%vp,grid%nx,grid%ny,grid%nz,& 
                            grid%dx, grid%dy, grid%dz/grid%scaling,&
                            grid%xmin,grid%ymin,grid%zmin/grid%scaling,dat%rev(:,i),&
                            nsrc,rti%locations(1:3,:),settings%order,like%phaseTime(i,:,1))
                        call fastmarching3d(model%vs,grid%nx,grid%ny,grid%nz,& 
                            grid%dx, grid%dy, grid%dz/grid%scaling,&
                            grid%xmin,grid%ymin,grid%zmin/grid%scaling,dat%rev(:,i),&
                            nsrc,rti%locations(1:3,:),settings%order,like%phaseTime(i,:,2))
                    endif
                enddo
                !$omp end do
                !$omp end parallel
                !t2 = omp_get_wtime ( )
            endif
        else
#ifdef _OPENMP
            !t1 = omp_get_wtime ( )
            !call omp_set_num_threads(settings%nthreads)
#endif
            !$omp parallel
            !$omp do private(i,j)
            do i = 1, nrev
                if (settings%datatype==0)then
                    call fastmarching3d(model%vp,grid%nx,grid%ny,grid%nz,& 
                        grid%dx, grid%dy, grid%dz/grid%scaling,&
                        grid%xmin,grid%ymin,grid%zmin/grid%scaling,dat%rev(:,i),&
                        nsrc,rti%locations(1:3,:),settings%order,like%phaseTime(i,:,1),&
                        like%field5d(:,:,:,i,1))
                else
                    call fastmarching3d(model%vp,grid%nx,grid%ny,grid%nz,& 
                        grid%dx, grid%dy, grid%dz/grid%scaling,&
                        grid%xmin,grid%ymin,grid%zmin/grid%scaling,dat%rev(:,i),&
                        nsrc,rti%locations(1:3,:),settings%order,like%phaseTime(i,:,1),&
                        like%field5d(:,:,:,i,1))
                    call fastmarching3d(model%vs,grid%nx,grid%ny,grid%nz,& 
                        grid%dx, grid%dy, grid%dz/grid%scaling,&
                        grid%xmin,grid%ymin,grid%zmin/grid%scaling,dat%rev(:,i),&
                        nsrc,rti%locations(1:3,:),settings%order,like%phaseTime(i,:,2),&
                        like%field5d(:,:,:,i,2))
                endif
            enddo
            !$omp end do
            !$omp end parallel
            !t2 = omp_get_wtime ( )
        endif
        !call log_msg(itoa(omp_get_num_threads() ) )
        !call log_msg('parallelized dispersion curve code: ' //rtoa(t2-t1) )
        !call write_vel(pvel,'phaseVel.dat')
        ! discard models
        !if(any(ierr==1))then
        !    like%like = -1
        !    RTI%num_bad_model = RTI%num_bad_model + 1
        !    return
        !endif
        
        ! using travel time instead of source-receiver distance
        do j = 1, nsrc
            do k = 1, nrev
                nrr = k + (j-1)*nrev
                like%srdist(nrr,:) = like%phaseTime(k,j,:)
            enddo
        enddo

        !! likelihood
        !! noise level
        if(settings%sigdep /= 0)then
            do i = 1, dat%np
                like%sigma(:,i) = RTI%bnoise0(i)*like%srdist(:,i) + RTI%bnoise1(i)
            enddo
        else
            like%sigma = dat%ttime(:,2,:)
        endif

        like%like = 0
        like%misfit = 0
        like%unweighted_misfit = 0
        do i = 1, dat%np
            nrr  = 0
            do j = 1, nsrc
                   do k = 1, nrev
                    nrr = nrr + 1
                    if(dat%raystat(nrr,1,i) == 1) then
                        if(like%sigma(nrr,i)<EPS) call exception_raiseError('The noise level is 0!')
                        like%like = like%like + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
                            -dat%ttime(nrr,1,i))**2/( 2*(like%sigma(nrr,i))**2 )
                        like%misfit = like%misfit + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
                            -dat%ttime(nrr,1,i))**2/(like%sigma(nrr,i)**2)
                        like%unweighted_misfit = like%unweighted_misfit + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
                            -dat%ttime(nrr,1,i))**2
                    else
                        like%sigma(nrr,i) = 1.0
                    endif
                enddo
            enddo
        enddo

        like%like = like%like + sum(log(like%sigma)) + dat%nrays/2.0 * log(PI2)

        return

    end subroutine body_likelihood

    subroutine body_likelihood_grads(dat,model,RTI,settings,like)
    
        implicit none
        type(T_DATA), intent(in)                    :: dat
        type(T_MOD), intent(in)                     :: model
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(T_LIKE_SET), intent(in)                :: settings
        type(T_LIKE_BASE), intent(inout)            :: like

        ! > local variable for travel times and rays
        integer(c_int) nsrc, nrev
        real(kind=ii10), dimension(4*dat%nsrc) :: grads_src, grads_src2

        ! > local grid
        type(T_GRID) :: grid

        integer nrr, nsites
        integer i, j, k

        nsrc = dat%nsrc
        nrev = dat%nrev
        grid = settings%grid
        nsites = RTI%ncells
        grads_src = 0
        grads_src2 = 0

        ! noise level
        if(settings%sigdep /= 0)then
            do i = 1, dat%np
                like%sigma(:,i) = RTI%bnoise0(i)*like%srdist(:,i) + RTI%bnoise1(i)
            enddo
        else
            like%sigma = dat%ttime(:,2,:)
        endif
        ! set no data sigma to 1, this is important since sigma is log added to
        ! lglikelihood
        do i = 1, dat%np
            nrr  = 0
            do j = 1, nsrc
                   do k = 1, nrev
                    nrr = nrr + 1
                    if(dat%raystat(nrr,1,i) /= 1) then
                        like%sigma(nrr,i) = 1.0
                    endif
                enddo
            enddo
        enddo


        if(nsrc<=nrev)then

            call trace_from_source(settings,model%vp,RTI,dat%rev,reshape(dat%ttime(:,1,1),[nrev,nsrc]),&
                reshape(like%sigma(:,1),[nrev,nsrc]),like%phaseTime(:,:,1),like%grads(1:nsites),grads_src)
            like%grads(1:nsites) = -like%grads(1:nsites)/RTI%parameters(1,1:nsites)**2
            grads_src2 = grads_src2 + grads_src
            if(settings%datatype>0)then
                call trace_from_source(settings,model%vs,RTI,dat%rev,reshape(dat%ttime(:,1,2),[nrev,nsrc]),&
                    reshape(like%sigma(:,2),[nrev,nsrc]),like%phaseTime(:,:,2),like%grads(nsites+1:2*nsites),grads_src)
                like%grads(nsites+1:2*nsites) = -like%grads(nsites+1:2*nsites)/RTI%parameters(2,1:nsites)**2
                grads_src2 = grads_src2 + grads_src
            endif
        else
            call trace_from_receiver(settings,model%vp,RTI,dat%rev,reshape(dat%ttime(:,1,1),[nrev, nsrc]),&
                reshape(like%sigma(:,1),[nrev,nsrc]),like%phaseTime(:,:,1),like%grads(1:nsites),grads_src)
            like%grads(1:nsites) = -like%grads(1:nsites)/RTI%parameters(1,1:nsites)**2
            grads_src2 = grads_src2 + grads_src
            if(settings%datatype>0)then
                call trace_from_receiver(settings,model%vs,RTI,dat%rev,reshape(dat%ttime(:,1,2),[nrev, nsrc]),&
                    reshape(like%sigma(:,2),[nrev,nsrc]),like%phaseTime(:,:,2),like%grads(nsites+1:2*nsites),grads_src)
                like%grads(nsites+1:2*nsites) = -like%grads(nsites+1:2*nsites)/RTI%parameters(2,1:nsites)**2
                grads_src2 = grads_src2 + grads_src
            endif
        endif
        like%grads(2*nsites+1:ubound(like%grads,1)) = grads_src2

        like%like = 0
        like%misfit = 0
        like%unweighted_misfit = 0
        do i = 1, dat%np
            nrr  = 0
            do j = 1, nsrc
                   do k = 1, nrev
                    nrr = nrr + 1
                    if(dat%raystat(nrr,1,i) == 1) then
                        if(like%sigma(nrr,i)<EPS) call exception_raiseError('The noise level is 0!')
                        like%like = like%like + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
                            -dat%ttime(nrr,1,i))**2/( 2*(like%sigma(nrr,i))**2 )
                        like%misfit = like%misfit + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
                            -dat%ttime(nrr,1,i))**2/(like%sigma(nrr,i)**2)
                        like%unweighted_misfit = like%unweighted_misfit + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
                            -dat%ttime(nrr,1,i))**2
                    endif
                enddo
            enddo
        enddo

        like%like = like%like + sum(log(like%sigma)) + dat%nrays/2.0 * log(PI2)

        return

    end subroutine body_likelihood_grads

    ! likelihood gradient given we know the ray path already
    !subroutine body_likelihood_grads_approx(dat,RTI,settings,like)
    !
    !    implicit none
    !    type(T_DATA), intent(in)                    :: dat
    !    type(T_RUN_INFO), intent(inout)             :: RTI
    !    type(T_LIKE_SET), intent(in)                :: settings
    !    type(T_LIKE_BASE), intent(inout)            :: like

    !    ! > local variable for travel times and rays
    !    integer(c_int) nsrc, nrev
    !    real(kind=ii10), dimension(size(RTI%points,2)) :: grads_src, grads_src2

    !    ! > local grid
    !    type(T_GRID) :: grid
    !    type(T_MOD)  :: model

    !    integer nrr, nsites
    !    integer i, j, k

    !    nsrc = dat%nsrc
    !    nrev = dat%nrev
    !    grid = settings%grid
    !    nsites = RTI%ncells
    !    grads_src = 0
    !    grads_src2 = 0

    !    call mod_setup(model,settings%grid)
    !    do i = 1, settings%grid%nx
    !        do j = 1, settings%grid%ny
    !            do k = 1, settings%grid%nz
    !                model%vp(k,j,i) = RTI%parameters(1,RTI%sites_id(k,j,i))
    !                model%vs(k,j,i) = RTI%parameters(2,RTI%sites_id(k,j,i))
    !                model%rho(k,j,i) = RTI%parameters(3,RTI%sites_id(k,j,i))
    !            enddo
    !        enddo
    !    enddo

    !    ! noise level
    !    if(settings%sigdep /= 0)then
    !        do i = 1, dat%np
    !            like%sigma(:,i) = RTI%bnoise0(i)*like%srdist(:,i) + RTI%bnoise1(i)
    !        enddo
    !    else
    !        like%sigma = dat%ttime(:,2,:)
    !    endif
    !    ! set no data sigma to 1, this is important since sigma is log added to
    !    ! lglikelihood
    !    do i = 1, dat%np
    !        nrr  = 0
    !        do j = 1, nsrc
    !               do k = 1, nrev
    !                nrr = nrr + 1
    !                if(dat%raystat(nrr,1,i) /= 1) then
    !                    like%sigma(nrr,i) = 1.0
    !                endif
    !            enddo
    !        enddo
    !    enddo

    !    do i = 1, nsrc
    !        do j = 1, nrev
    !            call integerate_time(like%ray(j,i,1), time_out(j,i,1))
    !            call derivative(rti%sites_id,model%vp,grid,like%ray(j,i),one_grads,one_grad_src,1)
    !            if(settings%datatype>0)then
    !                call derivative(rti%sites_id,model%vs,grid,like%ray(j,i),one_grads,one_grad_src,1)
    !            endif
    !            ! test finite difference derivatives of source locations
    !            !call derivative_src(rti%locations(1:3,j),timefield, grid, one_grad_src)

    !            residual_factor =(time_out(i,j)+RTI%locations(4,j)-time_data(i,j))/time_noise(i,j)**2  
    !            grads = grads + one_grads * residual_factor

    !            !src_idx = (j-1)*4 + 1
    !            !grads_src(src_idx:src_idx+3) = grads_src(src_idx:src_idx+3) + one_grad_src*residual_factor

    !        enddo
    !    enddo
    !        
    !    like%like = 0
    !    like%misfit = 0
    !    like%unweighted_misfit = 0
    !    do i = 1, dat%np
    !        nrr  = 0
    !        do j = 1, nsrc
    !               do k = 1, nrev
    !                nrr = nrr + 1
    !                if(dat%raystat(nrr,1,i) == 1) then
    !                    if(like%sigma(nrr,i)<EPS) call exception_raiseError('The noise level is 0!')
    !                    like%like = like%like + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
    !                        -dat%ttime(nrr,1,i))**2/( 2*(like%sigma(nrr,i))**2 )
    !                    like%misfit = like%misfit + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
    !                        -dat%ttime(nrr,1,i))**2/(like%sigma(nrr,i)**2)
    !                    like%unweighted_misfit = like%unweighted_misfit + ( like%phaseTime(k,j,i)+RTI%locations(4,j)&
    !                        -dat%ttime(nrr,1,i))**2
    !                endif
    !            enddo
    !        enddo
    !    enddo

    !    like%like = like%like + sum(log(like%sigma)) + dat%nrays/2.0 * log(PI2)

    !    return

    !end subroutine body_likelihood_grads_approx

    ! data is needed for gradient of -lglikelihood
    subroutine trace_from_source(settings,velocity,rti,rev,time_data,time_noise,time_out,grads,grads_src)
        implicit none
        type(T_LIKE_SET), intent(in) :: settings
        real(kind=ii10),dimension(:,:,:), intent(in) :: velocity
        type(T_RUN_INFO), intent(in) :: rti
        real(kind=ii10), dimension(:,:), intent(in) :: rev
        real(kind=ii10), dimension(:,:), intent(in) :: time_data
        real(kind=ii10), dimension(:,:), intent(in) :: time_noise
        real(kind=ii10), dimension(:,:), intent(out) :: time_out
        real(kind=ii10), dimension(:), intent(out) :: grads
        real(kind=ii10), dimension(:), intent(out) :: grads_src

        real(kind=ii10), dimension(size(velocity,1),size(velocity,2),size(velocity,3)) :: timefield
        real(kind=ii10), dimension(:), allocatable :: one_grads
        real(kind=ii10), dimension(4) :: one_grad_src
        integer i, j
        integer src_idx
        integer nsrc, nrev
        type(T_GRID) grid
        type(T_RAY) ray
        real(kind=ii10) residual_factor

        nsrc = size(time_data,2)
        nrev = size(time_data,1)
        grid = settings%grid

        ! for efficiency, allocate ray first
        allocate(ray%points(3,grid%nx*grid%ny*grid%nz))
        !allocate(ray%point_velocity(grid%nx*grid%ny*grid%nz))
        ray%points = 0
        !ray%point_velocity = 0

        time_out = 0.0
        grads = 0.0
        grads_src = 0.0
        one_grads = grads
        one_grad_src = 0.0
        !$omp parallel private(i,j,ray,timefield,one_grads,one_grad_src,residual_factor,src_idx) reduction(+:grads, grads_src)
        if(.not.allocated(ray%points)) allocate(ray%points(3,grid%nx*grid%ny*grid%nz))
        ray%points = 0
        !$omp do 
        do i = 1, nsrc
            call fastmarching3d(velocity,grid%nx,grid%ny,grid%nz,grid%dx, grid%dy, grid%dz/grid%scaling,&
                grid%xmin,grid%ymin,grid%zmin/grid%scaling,rti%locations(1:3,i),nrev, rev,settings%order,time_out(:,i),timefield)
            do j = 1, nrev
                ray%npoints = 0
                ray%points = 0
                !ray%point_velocity = 0
                if(abs(time_data(j,i)+1)>eps)then ! if data is valid
                    ray%srcid = i
                    ray%revid = j
                    call gradientdescent(grid%nx,grid%ny,grid%nz,grid%dx, grid%dy, grid%dz/grid%scaling,&
                        grid%xmin,grid%ymin,grid%zmin/grid%scaling,velocity,rti%locations(1:3,i),rev(:,j),timefield,&
                        ray%npoints, ray%points)
                    call derivative(rti%sites_id,rti%points(:,1:rti%ncells),grid,ray,one_grads,one_grad_src,0)

                    residual_factor =(time_out(j,i)+RTI%locations(4,i)-time_data(j,i))/time_noise(j,i)**2  
                    grads = grads + one_grads * residual_factor

                    src_idx = (i-1)*4 + 1
                    grads_src(src_idx:src_idx+3) = grads_src(src_idx:src_idx+3) + one_grad_src*residual_factor
                endif
            enddo
        enddo
        !$omp end do
        deallocate(ray%points)
        !$omp end parallel
    endsubroutine trace_from_source

    subroutine trace_from_receiver(settings,velocity,rti,rev,time_data,time_noise,time_out,grads,grads_src)
        implicit none
        type(T_LIKE_SET), intent(in) :: settings
        real(kind=ii10),dimension(:,:,:), intent(in) :: velocity
        type(T_RUN_INFO), intent(in) :: rti
        real(kind=ii10), dimension(:,:), intent(in) :: rev
        real(kind=ii10), dimension(:,:), intent(in) :: time_data
        real(kind=ii10), dimension(:,:), intent(in) :: time_noise
        real(kind=ii10), dimension(:,:), intent(out) :: time_out
        real(kind=ii10), dimension(:), intent(out) :: grads
        real(kind=ii10), dimension(:), intent(out) :: grads_src


        !real(kind=ii10), dimension(size(velocity,1),size(velocity,2),size(velocity,3)) :: timefield
        real(kind=ii10), dimension(:,:,:), allocatable :: timefield
        real(kind=ii10), dimension(:), allocatable :: one_grads
        real(kind=ii10), dimension(4) :: one_grad_src
        integer i, j
        integer src_idx
        integer nsrc, nrev
        type(T_GRID) grid
        type(T_RAY) ray
        real(kind=ii10) residual_factor

        nsrc = size(time_data,2)
        nrev = size(time_data,1)
        grid = settings%grid

        ! for efficiency, allocate ray first
        !allocate(ray%points(3,grid%nx*grid%ny*grid%nz))
        !allocate(ray%point_velocity(grid%nx*grid%ny*grid%nz))
        !ray%points = 0
        !ray%point_velocity = 0
        !allocate(timefield(grid%nz,grid%ny,grid%nx,nrev))
        !timefield = 0

        time_out = 0.0
        grads = 0.0
        grads_src = 0.0
        one_grads = grads
        one_grad_src = 0.0
        !$omp parallel private(i,j,ray,timefield,one_grads,one_grad_src,residual_factor,src_idx) reduction(+:grads, grads_src)
        if(.not.allocated(ray%points)) allocate(ray%points(3,grid%nx*grid%ny*grid%nz))
        ray%points = 0
        if(.not.allocated(timefield)) allocate(timefield(grid%nz,grid%ny,grid%nx))
        timefield = 0
        !$omp do 
        !if(.not.allocated(ray%point_velocity)) allocate(ray%point_velocity(grid%nx*grid%ny*grid%nz))
        do i = 1, nrev
            call fastmarching3d(velocity,grid%nx,grid%ny,grid%nz,grid%dx, grid%dy, grid%dz/grid%scaling,&
                grid%xmin,grid%ymin,grid%zmin/grid%scaling,rev(:,i),nsrc,rti%locations(1:3,:),settings%order,time_out(i,:),timefield)
            do j = 1, nsrc
                ray%npoints = 0
                ray%points = 0
                !ray%point_velocity = 0
                if(abs(time_data(i,j)+1)>eps)then ! if data is valid
                    ray%srcid = j
                    ray%revid = i
                    call gradientdescent(grid%nx,grid%ny,grid%nz,grid%dx, grid%dy, grid%dz/grid%scaling,&
                        grid%xmin,grid%ymin,grid%zmin/grid%scaling,velocity,rev(:,i),rti%locations(1:3,j),timefield,&
                        ray%npoints, ray%points)
                    !ray%npoints = size(ray%points,2)
                    !write(*,*) 'path: ', rev(:,i), rti%locations(1:3,j)
                    !write(*,*) 'path points: ', ray%points
                    call derivative(rti%sites_id,rti%points(:,1:rti%ncells),grid,ray,one_grads,one_grad_src,1)

                    ! test finite difference derivatives of source locations
                    call derivative_src(rti%locations(1:3,j),timefield, grid, one_grad_src)

                    residual_factor =(time_out(i,j)+RTI%locations(4,j)-time_data(i,j))/time_noise(i,j)**2  
                    grads = grads + one_grads * residual_factor

                    src_idx = (j-1)*4 + 1
                    grads_src(src_idx:src_idx+3) = grads_src(src_idx:src_idx+3) + one_grad_src*residual_factor
                endif
            enddo
        enddo
        !$omp end do
        deallocate(ray%points)
        deallocate(timefield)
        !$omp end parallel
        !write(*,*) 'Velocity grads:', grads
        !write(*,*) 'Source grads:', grads_src
    endsubroutine trace_from_receiver

    ! derivative of travel time to Voronoi velocities and source parameters
    ! for each ray
    subroutine derivative(sites_id,points,grid,ray,grads,grads_src,istart)
        implicit none
        integer, dimension(:,:,:), intent(in) :: sites_id
        real(kind=ii10), dimension(:,:), intent(in) :: points
        !real(kind=ii10), dimension(:), intent(in) :: parameters
        type(T_GRID), intent(in) :: grid
        type(T_RAY), intent(in) :: ray
        real(kind=ii10), dimension(:), intent(inout) :: grads
        real(kind=ii10), dimension(4), intent(out) :: grads_src
        integer, intent(in) :: istart

        type(kdtree2), pointer :: tree
        integer i, j, k
        integer ix, iy, iz, prev_id, curr_id
        integer src_idx
        real(kind=ii10) step, dt, ds
        real(kind=ii10), dimension(3) :: dpoint, pt

        grads = 0
        grads_src = 0
        step = minval([grid%dx,grid%dy,grid%dz/grid%scaling])*1.0

        ! if ray is not valid, return
        if(ray%npoints==1) return

        ! create kd-tree
        tree => kdtree2_create(points,sort=.false.,rearrange=.true.)

        ! currently, simply use the minimum corner id as the path point id
        !call point2idx(grid,ray%points(:,1),ix,iy,iz)
        !prev_id = sites_id(iz,iy,ix)
        prev_id =  sites_locate(tree,sites_id,grid,ray%points(:,1))
        do j = 2, ray%npoints
            !call point2idx(grid,ray%points(:,j),ix,iy,iz)
            !curr_id = sites_id(iz,iy,ix)
            curr_id = sites_locate(tree,sites_id,grid,ray%points(:,j))
            
            step = norm2(ray%points(:,j)-ray%points(:,j-1))
            if(prev_id==curr_id)then
                ! quick, not general, note id not checked for efficiency
                grads(curr_id) = grads(curr_id) + step
            else
                ds = step/10
                dpoint = ray%points(:,j) - ray%points(:,j-1)
                pt = ray%points(:,j-1)
                do k = 1, 10 
                    pt = pt + ds*dpoint/step
                    curr_id = sites_locate(tree,sites_id,grid,pt) 
                    grads(prev_id) = grads(prev_id) + ds/2.0
                    grads(curr_id) = grads(curr_id) + ds/2.0
                    prev_id =  curr_id
                enddo
            endif

            prev_id = curr_id
        enddo

        ! for source derivative
        if(istart==1)then
            !call point2idx(grid,ray%points(:,1),ix,iy,iz)
            curr_id = sites_locate(tree,sites_id,grid,ray%points(:,1))
            dt = step/1.0
            ! x, y, z, t
            grads_src(1) = -dt/(ray%points(1,2)-ray%points(1,1))
            grads_src(2) = -dt/(ray%points(2,2)-ray%points(2,1))
            grads_src(3) = -dt/(ray%points(3,2)-ray%points(3,1))
            grads_src(4) = 1.0
        else
            !call point2idx(grid,ray%points(:,ray%npoints),ix,iy,iz)
            curr_id = sites_locate(tree,sites_id,grid,ray%points(:,ray%npoints))
            dt = step/1.0
            src_idx = ray%npoints
            ! x, y, z, t
            grads_src(1) = -dt/(ray%points(1,src_idx-1)-ray%points(1,src_idx))
            grads_src(2) = -dt/(ray%points(2,src_idx-1)-ray%points(2,src_idx))
            grads_src(3) = -dt/(ray%points(3,src_idx-1)-ray%points(3,src_idx))
            grads_src(4) = 1.0
        endif
        !write(*,*) 'grad_src from ray', grads_src

        call kdtree2_destroy(tree)

    endsubroutine derivative

    integer function sites_locate(tree,sites_id,grid,point)
        implicit none
        type(kdtree2), pointer, intent(in) :: tree
        integer, dimension(:,:,:), intent(in) :: sites_id
        type(T_GRID), intent(in) :: grid
        real(kind=ii10), dimension(:), intent(in) :: point

        type(kdtree2_result), dimension(1) :: results
        integer ix, iy, iz, idx
        integer i, j, k, num 

        call point2idx(grid,point,ix,iy,iz)
        idx =  sites_id(iz,iy,ix)

        num = 0
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    num = num + abs(sites_id(iz+k-1,iy+j-1,ix+i-1)-idx)
                enddo
            enddo
        enddo

        if(num==0)then
            sites_locate = idx
        else
            call kdtree2_n_nearest(tp=tree,qv=point,nn=1,results=results)
            sites_locate = results(1)%idx
        endif

    end function

    ! derivative of travel time to Voronoi velocities and source parameters
    ! for each ray
    subroutine derivative2(sites_id,vel,grid,ray,grads,grads_src,istart)
        implicit none
        integer, dimension(:,:,:), intent(in) :: sites_id
        real(kind=ii10), dimension(:,:,:) :: vel
        type(T_GRID), intent(in) :: grid
        type(T_RAY), intent(in) :: ray
        real(kind=ii10), dimension(:), intent(inout) :: grads
        real(kind=ii10), dimension(4), intent(out) :: grads_src
        integer, intent(in) :: istart

        integer i, j, k, l, m, n
        integer ix, iy, iz, prev_id, curr_id
        integer pix, piy, piz
        integer src_idx
        real(kind=ii10) step, dt, ds
        real(kind=ii10), dimension(3) :: dpoint, pt
        real(kind=ii10), dimension(2,2,2) :: w
        integer, parameter :: nsteps = 5

        grads = 0
        grads_src = 0
        step = minval([grid%dx,grid%dy,grid%dz/grid%scaling])*1.0
        ds = step/nsteps

        ! if ray is not valid, return
        if(ray%npoints==1) return

        ! currently, simply use the minimum corner id as the path point id
        call point2idx(grid,ray%points(:,1),pix,piy,piz)
        prev_id = sites_id(piz,piy,pix)
        do j = 2, ray%npoints
            step = norm2(ray%points(:,j)-ray%points(:,j-1))
            ds = step/nsteps
            dpoint = ray%points(:,j) - ray%points(:,j-1)
            pt = ray%points(:,j-1)
            do k = 1, nsteps
                call GetWeight(grid,pt,w)
                call point2idx(grid,pt,ix,iy,iz)
                do l = 1, 2
                  do m = 1, 2
                    do n = 1, 2
                        piz = iz + (n-1)
                        piy = iy + (m-1)
                        pix = ix + (l-1)
                        grads(sites_id(piz,piy,pix)) = grads(sites_id(piz,piy,pix)) + ds*w(n,m,l)
                    enddo
                  enddo
                enddo
                pt = pt + ds*dpoint/step
            enddo
            if(abs(pt(1)-ray%points(1,j))>1E-5)then
                write(*,*) 'error: ', pt, ray%points(:,j)
            endif
        enddo

        ! for source derivative
        if(istart==1)then
            call point2idx(grid,ray%points(:,1),ix,iy,iz)
            dt = step/vel(iz,iy,ix)
            ! x, y, z, t
            grads_src(1) = -dt/(ray%points(1,2)-ray%points(1,1))
            grads_src(2) = -dt/(ray%points(2,2)-ray%points(2,1))
            grads_src(3) = -dt/(ray%points(3,2)-ray%points(3,1))
            grads_src(4) = 1.0
        else
            call point2idx(grid,ray%points(:,ray%npoints),ix,iy,iz)
            dt = step/vel(iz,iy,ix)
            src_idx = ray%npoints
            ! x, y, z, t
            grads_src(1) = -dt/(ray%points(1,src_idx-1)-ray%points(1,src_idx))
            grads_src(2) = -dt/(ray%points(2,src_idx-1)-ray%points(2,src_idx))
            grads_src(3) = -dt/(ray%points(3,src_idx-1)-ray%points(3,src_idx))
            grads_src(4) = 1.0
        endif

    endsubroutine derivative2

    ! using finite difference to calculate source locations derivatives
    subroutine derivative_src(src, timefield, grid, grads_src)
        implicit none
        real(kind=ii10), dimension(:), intent(in) :: src
        real(kind=ii10), dimension(:,:,:), intent(in) :: timefield
        type(T_GRID), intent(in) :: grid
        real(kind=ii10), dimension(4), intent(out) :: grads_src

        integer i, j
        real(kind=ii10) step, t1, t2
        real(kind=ii10), dimension(3) :: dsrc

        step = minval([grid%dx,grid%dy,grid%dz/grid%scaling])/2

        ! x direction
        dsrc = [step, 0.0_ii10, 0.0_ii10]
        t1 = Interpolate3d(timefield,grid,src-dsrc)
        t2 = Interpolate3d(timefield,grid,src+dsrc)
        grads_src(1) = (t2-t1)/(2*step)

        ! y direction
        dsrc = [0.0_ii10, step, 0.0_ii10]
        t1 = Interpolate3d(timefield,grid,src-dsrc)
        t2 = Interpolate3d(timefield,grid,src+dsrc)
        grads_src(2) = (t2-t1)/(2*step)

        ! z direction
        dsrc = [0.0_ii10, 0.0_ii10, step]
        t1 = Interpolate3d(timefield,grid,src-dsrc)
        t2 = Interpolate3d(timefield,grid,src+dsrc)
        grads_src(3) = (t2-t1)/(2*step)

        grads_src(4) = 1.0

        !write(*,*) 'grad_src from timefield', grads_src
    endsubroutine derivative_src

    function Interpolate2d(array,grid,point) result(qv)
        implicit none
        real(kind=ii10), dimension(:,:), intent(in) :: array
        type(T_GRID), intent(in) :: grid
        real(kind=ii10), dimension(:), intent(in) :: point
        real(kind=ii10) :: qv

        integer ix, iy, i, j
        real(kind=ii10) dsx, dsy
        real(kind=ii10) weight

        ix = floor((point(1)-grid%xmin)/grid%dx) + 1
        iy = floor((point(2)-grid%ymin)/grid%dy) + 1
        if(ix < 1) ix = 1
        if(iy < 1) iy = 1
        if(ix >= grid%nx) ix = grid%nx - 1
        if(iy >= grid%ny) iy = grid%ny - 1
        dsx = point(1) - (grid%xmin + (ix-1)*grid%dx)
        dsy = point(2) - (grid%ymin + (iy-1)*grid%dy)

        qv = 0
        do i = 1, 2
            do j = 1, 2
                weight = (1.0-abs((i-1)*grid%dx-dsx)/grid%dx)*(1.0-abs((j-1)*grid%dy-dsy)/grid%dy)
                qv = qv + weight*array(iy+j-1,ix+i-1)
            enddo
        enddo
    end function

    function Interpolate3d(array,grid,point) result(qv)
        implicit none
        real(kind=ii10), dimension(:,:,:), intent(in) :: array
        type(T_GRID), intent(in) :: grid
        real(kind=ii10), dimension(:), intent(in) :: point
        real(kind=ii10) :: qv

        integer ix, iy, iz, i, j, k
        real(kind=ii10) dsx, dsy, dsz
        real(kind=ii10) weight

        call point2idx(grid,point,ix,iy,iz)
        dsx = point(1) - (grid%xmin + (ix-1)*grid%dx)
        dsy = point(2) - (grid%ymin + (iy-1)*grid%dy)
        dsz = point(3) - (grid%zmin + (iz-1)*grid%dz)/grid%scaling

        qv = 0
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    weight = (1.0-abs((i-1)*grid%dx-dsx)/grid%dx)*&
                             (1.0-abs((j-1)*grid%dy-dsy)/grid%dy)*&
                             (1.0-abs((k-1)*grid%dz/grid%scaling-dsz)/(grid%dz/grid%scaling))
                    qv = qv + weight*array(iz+k-1,iy+j-1,ix+i-1)
                enddo
            enddo
        enddo

        ! debug
        !print *, array(iz:iz+1,iy:iy+1,ix:ix+1)
        !print *, qv
    end function

    subroutine GetWeight(grid,point,w)
        implicit none
        type(T_GRID), intent(in) :: grid
        real(kind=ii10), dimension(:), intent(in) :: point
        real(kind=ii10), dimension(2,2,2), intent(out) :: w

        integer ix, iy, iz, i, j, k
        real(kind=ii10) dsx, dsy, dsz

        call point2idx(grid,point,ix,iy,iz)
        dsx = point(1) - (grid%xmin + (ix-1)*grid%dx)
        dsy = point(2) - (grid%ymin + (iy-1)*grid%dy)
        dsz = point(3) - (grid%zmin + (iz-1)*grid%dz)/grid%scaling

        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    w(k,j,i) = (1.0-abs((i-1)*grid%dx-dsx)/grid%dx)*&
                             (1.0-abs((j-1)*grid%dy-dsy)/grid%dy)*&
                             (1.0-abs((k-1)*grid%dz/grid%scaling-dsz)/(grid%dz/grid%scaling))
                enddo
            enddo
        enddo

        ! debug
        !print *, array(iz:iz+1,iy:iy+1,ix:ix+1)
        !print *, qv
    end subroutine

    pure subroutine point2idx(grid,point,ix,iy,iz)
        implicit none
        type(T_GRID), intent(in) :: grid
        real(kind=c_double), dimension(3), intent(in) :: point
        integer(c_int), intent(out) :: ix, iy, iz

        ix = floor((point(1)-grid%xmin)/grid%dx) + 1
        iy = floor((point(2)-grid%ymin)/grid%dy) + 1
        iz = floor((point(3)-grid%zmin/grid%scaling)/(grid%dz/grid%scaling)) + 1
        if(ix < 1) ix = 1
        if(iy < 1) iy = 1
        if(iz < 1) iz = 1
        if(ix >= grid%nx) ix = grid%nx - 1
        if(iy >= grid%ny) iy = grid%ny - 1
        if(iz >= grid%nz) iz = grid%nz - 1

    endsubroutine point2idx

end module m_body_likelihood
