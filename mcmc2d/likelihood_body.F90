!
! module m_likelihood
! likelihood function here
!
module m_likelihood

    use m_exception, only : exception_raiseWarning, exception_raiseError
    use m_logger, only    : log_msg
    use m_utils, only     : ii10, EPS, write_resume_unit, write_doubles, itoa, rtoa, vs2vp, vp2rho
    use m_settings, only  : T_GRID, T_MOD, mod_setup, out_bnd
    use run_info, only    : T_RUN_INFO
    use m_fm2d, only      : modrays, T_RAY
    use like_settings,only: T_LIKE_BASE, T_LIKE_SET, T_DATA

    use omp_lib

    implicit none
    private

    public :: body_likelihood, body_noise_likelihood
    public :: body_like_grads, update_lgP_grads_approx

    ! debug
    real( kind=ii10 ) :: t1, t2
    ! static value
    real( kind=ii10 ), parameter :: ERAD = 6371
    real( kind=ii10 ), parameter :: PI2 = 6.283185

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
            like%sigma(:,i) = RTI%snoise0(i)*like%srdist(:,i) + RTI%snoise1(i)
        enddo

        like%like = 0
        like%misfit = 0
        nrr  = 0
        do j = 1, nsrc
               do k = 1, nrev
                nrr = nrr + 1
                if(dat%raystat(nrr,1,1) == 1) then
                    if(like%sigma(nrr,1)==0) call exception_raiseError('The noise level is 0!')
                    like%like = like%like + ( like%phaseTime(k,j,1)-dat%ttime(nrr,1,1)&
                        )**2/( 2*(like%sigma(nrr,1))**2 )
                    like%misfit = like%misfit + ( like%phaseTime(k,j,1)-dat%ttime(nrr,1,1)&
                        )**2/(like%sigma(nrr,1)**2)
                else
                    like%sigma(nrr,1) = 1.0
                endif
            enddo
        enddo

        if(any(like%sigma<EPS))then
            call exception_raiseError('The noise level is 0!')
        endif

        like%like = like%like + sum(log(like%sigma)) + dat%nrays/2.0 * log(PI2)
        return
    end subroutine body_noise_likelihood

    subroutine body_likelihood_grads(dat,RTI,settings,like)
        implicit none
        type(T_DATA), intent(in)                    :: dat
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(T_LIKE_SET), intent(in)                :: settings
        type(T_LIKE_BASE), intent(inout)            :: like

        ! > local variable for fast marching
        integer gridx, gridy
        integer sgref
        integer sgdic, sgext
        integer order
        integer uar
        real( kind=ii10 ) band
    
        ! > local variable for travel times and rays
        integer nsrc, nrev
        integer :: crazyray
        type(T_RAY), dimension(:,:), allocatable           :: rays
        real(kind=ii10), dimension(settings%grid%ny+2,settings%grid%nx+2) :: ext_vel

        ! > local grid
        type(T_GRID) :: grid
        type(T_MOD)  :: model

        integer nrr
        integer i, j, k

        nsrc = dat%nsrc
        nrev = dat%nrev
        grid = settings%grid
        
        call mod_setup(model,settings%grid)
        do i = 1, settings%grid%nx
            do j = 1, settings%grid%ny
                model%vp(j,i) = RTI%parameters(1,RTI%sites_id(j,i))
                model%vs(j,i) = RTI%parameters(2,RTI%sites_id(j,i))
                model%rho(j,i) = RTI%parameters(3,RTI%sites_id(j,i))
            enddo
        enddo

        ! noise level
        if(settings%sigdep /= 0)then
            do i = 1, dat%np
                like%sigma(:,i) = RTI%snoise0(i)*like%srdist(:,i) + RTI%snoise1(i)
            enddo
        else
            like%sigma = dat%ttime(:,2,:)
        endif

        ! if using straight rays
        if(settings%isStraight == 1)then
            if(.not.like%straightRaySet)then
                allocate( like%rays(dat%nrev*dat%nsrc,1) )
                like%rays%srcid = 0
                like%rays%revid = 0
                like%rays%npoints = 0
                call setup_straightRays(dat, grid, like%rays,like%srdist)
                like%straightRaySet = .true.
            endif
            call IntegrateRays(model%vp,grid,like%rays(:,1),like%phaseTime(:,:,1))
        else
            ! calculate travel time of rayleigh/love wave using fast marching code
            ! settings
            gridx = settings%gridx
            gridy = settings%gridy
            sgref = settings%sgref
            sgdic = settings%sgdic
            sgext = settings%sgext
            order = settings%order
            band  = settings%band
            uar   = 0

            ! prepare velocity model for the fast marching code
            ! assign boundary value
            ext_vel(2:grid%ny+1,2:grid%nx+1) = model%vp
            ext_vel(:,1) = ext_vel(:,2)
            ext_vel(:,grid%nx+2) = ext_vel(:,grid%nx+1)
            ext_vel(1,:) = ext_vel(2,:)
            ext_vel(grid%ny+2,:) = ext_vel(grid%ny+1,:)

            ! prepare the beginning time
            ! call modrays to calculate travel times
            allocate( rays(dat%nrev*dat%nsrc,1) )
            rays%srcid = 0
            rays%revid = 0
            rays%npoints = 0
            crazyray = 0
#ifdef _OPENMP
            t1 = omp_get_wtime ( )
            !call omp_set_num_threads(settings%nthreads)
#endif
            !!$omp parallel
            !!$omp do private(nrr,i,j,k)
            call modrays(nsrc,dat%src(1,:),dat%src(2,:), &
                       nrev,dat%rev(1,:),dat%rev(2,:), &
                   dat%raystat(:,:,1),0, &
                   grid%nx,grid%ny,grid%xmin,grid%ymin/grid%scaling,&
                   grid%dx,grid%dy/grid%scaling,ext_vel, &
                   gridx,gridy,sgref, &
                   sgdic,sgext,ERAD, &
                   order,band,like%phaseTime(:,:,1), &
                   rays(:,1),crazyray,uar)
            !!$omp end do
            !call log_msg(itoa(omp_get_num_threads() ) )
            !!$omp end parallel
            t2 = omp_get_wtime ( )
            !call log_msg('parallelized fast marching code: ' //rtoa(t2-t1) )
    
        endif
    
        !call derivative2(RTI%points(:,1:RTI%ncells),RTI%parameters(1,1:RTI%ncells),grid,rays(:,1),&
        !    reshape(dat%ttime(:,1,1),[nrev,nsrc]),reshape(like%sigma(:,1),[nrev,nsrc]),&
        !    like%phaseTime(:,:,1),like%grads)
        call derivative(RTI%sites_id,RTI%points(:,1:RTI%ncells),RTI%parameters(1,1:RTI%ncells),grid,rays(:,1),&
            reshape(dat%ttime(:,1,1),[nrev,nsrc]),reshape(like%sigma(:,1),[nrev,nsrc]),&
            like%phaseTime(:,:,1),like%grads)
        ! velocity to m/s testing
        !like%grads = like%grads/1000


        like%like = 0
        like%misfit = 0
        like%unweighted_misfit = 0
        nrr  = 0
        do j = 1, nsrc
               do k = 1, nrev
                nrr = nrr + 1
                if(dat%raystat(nrr,1,1) == 1) then
                    like%like = like%like + ( like%phaseTime(k,j,1)-dat%ttime(nrr,1,1)&
                        )**2/( 2*(like%sigma(nrr,1))**2 )
                    like%misfit = like%misfit + ( like%phaseTime(k,j,1)-dat%ttime(nrr,1,1)&
                        )**2/(like%sigma(nrr,1)**2)
                    like%unweighted_misfit = like%unweighted_misfit + ( like%phaseTime(k,j,1)-dat%ttime(nrr,1,1)&
                        )**2
                else
                    like%sigma(nrr,1) = 1.0
                endif
            enddo
        enddo

        like%like = like%like + sum(log(like%sigma)) + dat%nrays/2.0 * log(PI2)
        if(like%like/=like%like)then
            print *, like%sigma
        endif
        RTI%like_count = RTI%like_count + 1
        return

    endsubroutine body_likelihood_grads

    subroutine update_lgP_grads_approx(dat,RTI,settings,like)
        implicit none
        type(T_DATA), intent(in)                    :: dat
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(T_LIKE_SET), intent(in)                :: settings
        type(T_LIKE_BASE), intent(inout)            :: like

        ! > local variable for fast marching
        integer gridx, gridy
        integer sgref
        integer sgdic, sgext
        integer order
        integer uar
        real( kind=ii10 ) band
    
        ! > local variable for travel times and rays
        integer nsrc, nrev
        integer :: crazyray
        type(T_RAY), dimension(:,:), allocatable           :: rays
        real(kind=ii10), dimension(settings%grid%ny+2,settings%grid%nx+2) :: ext_vel

        ! > local grid
        type(T_GRID) :: grid
        type(T_MOD)  :: model

        integer nrr
        integer i, j, k

        nsrc = dat%nsrc
        nrev = dat%nrev
        grid = settings%grid
        
        call mod_setup(model,settings%grid)
        do i = 1, settings%grid%nx
            do j = 1, settings%grid%ny
                model%vp(j,i) = RTI%parameters(1,RTI%sites_id(j,i))
                model%vs(j,i) = RTI%parameters(2,RTI%sites_id(j,i))
                model%rho(j,i) = RTI%parameters(3,RTI%sites_id(j,i))
            enddo
        enddo

        ! noise level
        if(settings%sigdep /= 0)then
            do i = 1, dat%np
                like%sigma(:,i) = RTI%snoise0(i)*like%srdist(:,i) + RTI%snoise1(i)
            enddo
        else
            like%sigma = dat%ttime(:,2,:)
        endif

        ! if using straight rays
        call IntegrateRays(model%vp,grid,like%rays(:,1),like%phaseTime(:,:,1))
        !call derivative2(RTI%points(:,1:RTI%ncells),RTI%parameters(1,1:RTI%ncells),grid,rays(:,1),&
        !    reshape(dat%ttime(:,1,1),[nrev,nsrc]),reshape(like%sigma(:,1),[nrev,nsrc]),&
        !    like%phaseTime(:,:,1),like%grads)
        call derivative(RTI%sites_id,RTI%points(:,1:RTI%ncells),RTI%parameters(1,1:RTI%ncells),grid,like%rays(:,1),&
            reshape(dat%ttime(:,1,1),[nrev,nsrc]),reshape(like%sigma(:,1),[nrev,nsrc]),&
            like%phaseTime(:,:,1),like%grads)
        ! velocity to m/s testing
        !like%grads = like%grads/1000


        like%like = 0
        like%misfit = 0
        like%unweighted_misfit = 0
        nrr  = 0
        do j = 1, nsrc
               do k = 1, nrev
                nrr = nrr + 1
                if(dat%raystat(nrr,1,1) == 1) then
                    like%like = like%like + ( like%phaseTime(k,j,1)-dat%ttime(nrr,1,1)&
                        )**2/( 2*(like%sigma(nrr,1))**2 )
                    like%misfit = like%misfit + ( like%phaseTime(k,j,1)-dat%ttime(nrr,1,1)&
                        )**2/(like%sigma(nrr,1)**2)
                    like%unweighted_misfit = like%unweighted_misfit + ( like%phaseTime(k,j,1)-dat%ttime(nrr,1,1)&
                        )**2
                else
                    like%sigma(nrr,1) = 1.0
                endif
            enddo
        enddo

        like%like = like%like + sum(log(like%sigma)) + dat%nrays/2.0 * log(PI2)
        if(like%like/=like%like)then
            print *, like%sigma
        endif
        RTI%like_count = RTI%like_count + 1
        return

    endsubroutine update_lgP_grads_approx

    subroutine body_likelihood(dat,RTI,settings,like)
    
        implicit none
        type(T_DATA), intent(in)                    :: dat
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(T_LIKE_SET), intent(in)                :: settings
        type(T_LIKE_BASE), intent(inout)            :: like

        ! > local variable for fast marching
        integer gridx, gridy
        integer sgref
        integer sgdic, sgext
        integer order
        integer uar
        real( kind=ii10 ) band
    
        ! > local variable for travel times and rays
        integer nsrc, nrev
        integer :: crazyray
        type(T_RAY), dimension(:,:), allocatable           :: rays
        real(kind=ii10), dimension(settings%grid%ny,settings%grid%nx) :: vel
        real(kind=ii10), dimension(settings%grid%ny+2,settings%grid%nx+2) :: ext_vel

        ! > local grid
        type(T_GRID) :: grid

        integer nrr
        integer i, j, k

        nsrc = dat%nsrc
        nrev = dat%nrev
        grid = settings%grid

        if(settings%sigdep /= 0)then
            do i = 1, dat%np
                like%sigma(:,i) = RTI%snoise0(i)*like%srdist(:,i) + RTI%snoise1(i)
            enddo
        else
            like%sigma = dat%ttime(:,2,:)
        endif

        do i = 1, settings%grid%nx
            do j = 1, settings%grid%ny
                vel(j,i) = RTI%parameters(1,RTI%sites_id(j,i))
                !model%vs(j,i) = RTI%parameters(2,RTI%sites_id(j,i))
                !model%rho(j,i) = RTI%parameters(3,RTI%sites_id(j,i))
            enddo
        enddo

        ! if using straight rays
        if(settings%isStraight == 1)then
            if(.not.like%straightRaySet)then
                allocate( like%rays(dat%nrev*dat%nsrc,1) )
                like%rays%srcid = 0
                like%rays%revid = 0
                like%rays%npoints = 0
                call setup_straightRays(dat, grid, like%rays,like%srdist)
                like%straightRaySet = .true.
            endif
            call IntegrateRays(vel,grid,like%rays(:,1),like%phaseTime(:,:,1))
        else
            ! calculate travel time of rayleigh/love wave using fast marching code
            ! settings
            gridx = settings%gridx
            gridy = settings%gridy
            sgref = settings%sgref
            sgdic = settings%sgdic
            sgext = settings%sgext
            order = settings%order
            band  = settings%band
            uar   = 1

            ! prepare velocity model for the fast marching code
            ! assign boundary value
            ext_vel(2:grid%ny+1,2:grid%nx+1) = vel
            ext_vel(:,1) = ext_vel(:,2)
            ext_vel(:,grid%nx+2) = ext_vel(:,grid%nx+1)
            ext_vel(1,:) = ext_vel(2,:)
            ext_vel(grid%ny+2,:) = ext_vel(grid%ny+1,:)

            ! prepare the beginning time
            ! call modrays to calculate travel times
            allocate( rays(dat%nrev*dat%nsrc,1) )
            rays%srcid = 0
            rays%revid = 0
            rays%npoints = 0
            crazyray = 0
#ifdef _OPENMP
            t1 = omp_get_wtime ( )
            !call omp_set_num_threads(settings%nthreads)
#endif
            !!$omp parallel
            !!$omp do private(nrr,i,j,k)
            call modrays(nsrc,dat%src(1,:),dat%src(2,:), &
                       nrev,dat%rev(1,:),dat%rev(2,:), &
                   dat%raystat(:,:,1),0, &
                   grid%nx,grid%ny,grid%xmin,grid%ymin/grid%scaling,&
                   grid%dx,grid%dy/grid%scaling,ext_vel, &
                   gridx,gridy,sgref, &
                   sgdic,sgext,ERAD, &
                   order,band,like%phaseTime(:,:,1), &
                   rays(:,1),crazyray,uar)
            !!$omp end do
            !call log_msg(itoa(omp_get_num_threads() ) )
            !!$omp end parallel
            t2 = omp_get_wtime ( )
            !call log_msg('parallelized fast marching code: ' //rtoa(t2-t1) )
    
        endif


        like%like = 0
        like%misfit = 0
        like%unweighted_misfit = 0
        nrr  = 0
        do j = 1, nsrc
               do k = 1, nrev
                nrr = nrr + 1
                if(dat%raystat(nrr,1,1) == 1) then
                    like%like = like%like + ( like%phaseTime(k,j,1)-dat%ttime(nrr,1,1)&
                        )**2/( 2*(like%sigma(nrr,1))**2 )
                    like%misfit = like%misfit + ( like%phaseTime(k,j,1)-dat%ttime(nrr,1,1)&
                        )**2/(like%sigma(nrr,1)**2)
                    like%unweighted_misfit = like%unweighted_misfit + ( like%phaseTime(k,j,1)-dat%ttime(nrr,1,1)&
                        )**2
                else
                    like%sigma(nrr,1) = 1.0
                endif
            enddo
        enddo

        like%like = like%like + sum(log(like%sigma)) + dat%nrays/2.0 * log(PI2)
        if(like%like/=like%like)then
            print *, like%sigma
        endif
        RTI%like_count = RTI%like_count + 1
        return

    end subroutine body_likelihood

    ! derivative of travel time to Voronoi velocities
    ! for each ray
    subroutine derivative2(points,parameters,grid,rays,time_data,time_noise,time_out,grads)
        use kdtree2_precision_module, only : kdkind
        use kdtree2_module, only    : kdtree2, kdtree2_result, kdtree2_create,&
                                  kdtree2_n_nearest
        implicit none
        real(kind=ii10), dimension(:,:), intent(in) :: points
        real(kind=ii10), dimension(:), intent(in) :: parameters
        type(T_GRID), intent(in) :: grid
        type(T_RAY), dimension(:), intent(in) :: rays
        real(kind=ii10), dimension(:,:), intent(in) :: time_data
        real(kind=ii10), dimension(:,:), intent(in) :: time_noise
        real(kind=ii10), dimension(:,:), intent(in) :: time_out
        real(kind=ii10), dimension(:), intent(inout) :: grads

        ! kd-tree
        type(kdtree2), pointer :: tree
        type(kdtree2_result), dimension(1) :: results
        real(kdkind), dimension(3) :: qv

        real(kind=ii10), dimension(:), allocatable :: one_grads
        real(kind=ii10) residual_factor
        integer i, j
        integer prev_id, curr_id
        integer src_idx, rev_idx
        real(kind=ii10) step, dt

        grads = 0
        one_grads = grads

        ! create kd-tree
        tree => kdtree2_create(points(:,:), sort=.false., rearrange=.true.)
        
        !$omp parallel private(i,j,one_grads,residual_factor,src_idx,rev_idx,step,results,prev_id,curr_id) reduction(+:grads)
        !$omp do 
        do i = 1, size(rays)
            ! if ray is not valid, return
            if(rays(i)%npoints<=1) cycle

            one_grads = 0
            call kdtree2_n_nearest(tp=tree,qv=rays(i)%points(:,1),nn=1,results=results)
            prev_id = results(1)%idx
            do j = 2, rays(i)%npoints
                call kdtree2_n_nearest(tp=tree,qv=rays(i)%points(:,j),nn=1,results=results)
                curr_id = results(1)%idx

                step = norm2(rays(i)%points(:,j)-rays(i)%points(:,j-1))
                if(prev_id==curr_id)then
                    ! quick, not general, note id not checked for efficiency
                    one_grads(curr_id) = one_grads(curr_id) + step
                else
                    one_grads(prev_id) = one_grads(prev_id) + step/2.0
                    one_grads(curr_id) = one_grads(curr_id) + step/2.0
                endif

                prev_id = curr_id
            enddo

            src_idx = rays(i)%srcid
            rev_idx = rays(i)%revid
            residual_factor =(time_out(rev_idx,src_idx)-time_data(rev_idx,src_idx))/time_noise(rev_idx,src_idx)**2  
            grads = grads + one_grads * residual_factor
        enddo
        !$omp enddo
        !$omp end parallel
        grads(1:size(parameters)) = -grads(1:size(parameters))/parameters**2

    endsubroutine derivative2

    ! derivative of travel time to Voronoi velocities and source parameters
    ! for each ray
    subroutine derivative(sites_id,points,parameters,grid,rays,time_data,time_noise,time_out,grads)
        use kdtree2_precision_module, only : kdkind
        use kdtree2_module, only    : kdtree2, kdtree2_result, kdtree2_create,&
                                  kdtree2_n_nearest, kdtree2_destroy
        implicit none
        integer, dimension(:,:), intent(in) :: sites_id
        real(kind=ii10), dimension(:,:), intent(in) :: points
        real(kind=ii10), dimension(:), intent(in) :: parameters
        type(T_GRID), intent(in) :: grid
        type(T_RAY), dimension(:), intent(in) :: rays
        real(kind=ii10), dimension(:,:), intent(in) :: time_data
        real(kind=ii10), dimension(:,:), intent(in) :: time_noise
        real(kind=ii10), dimension(:,:), intent(in) :: time_out
        real(kind=ii10), dimension(:), intent(inout) :: grads

        type(kdtree2), pointer :: tree

        real(kind=ii10), dimension(:), allocatable :: one_grads
        real(kind=ii10) residual_factor
        integer i, j, k
        integer ix, iy, prev_id, curr_id
        integer src_idx, rev_idx
        real(kind=ii10) step, dt, ds
        real(kind=ii10), dimension(2) :: dpoint, pt

        grads = 0
        one_grads = grads
        step = minval([grid%dx,grid%dy/grid%scaling])*1.0

        ! create kd-tree
        tree => kdtree2_create(points(:,:), sort=.false., rearrange=.true.)
        
        ! if ray is not valid, return
        !$omp parallel private(i,j,ix,iy,one_grads,residual_factor,src_idx,rev_idx,step,prev_id,curr_id) reduction(+:grads)
        !$omp do 
        do i = 1, size(rays)
            ! if ray is not valid, return
            if(rays(i)%npoints<=1) cycle

            one_grads = 0
            ! currently, simply use the minimum corner id as the path point id
            !call point2idx(grid,rays(i)%points(:,1),ix,iy)
            !prev_id = sites_id(iy,ix)
            prev_id = sites_locate(tree,sites_id,grid,rays(i)%points(:,1))
            do j = 2, rays(i)%npoints
                !call point2idx(grid,rays(i)%points(:,j),ix,iy)
                !curr_id = sites_id(iy,ix)
                curr_id = sites_locate(tree,sites_id,grid,rays(i)%points(:,j))
                
                step = norm2(rays(i)%points(:,j)-rays(i)%points(:,j-1))
                if(prev_id==curr_id)then
                    ! quick, not general, note id not checked for efficiency
                    one_grads(curr_id) = one_grads(curr_id) + step
                else
                    ds = step/10
                    dpoint = rays(i)%points(:,j) - rays(i)%points(:,j-1)
                    pt = rays(i)%points(:,j-1)
                    do k = 1, 10 
                        pt = pt + ds*dpoint/step
                        curr_id = sites_locate(tree,sites_id,grid,pt) 
                        one_grads(prev_id) = one_grads(prev_id) + ds/2.0
                        one_grads(curr_id) = one_grads(curr_id) + ds/2.0
                        prev_id =  curr_id
                    enddo
                endif

                prev_id = curr_id
            enddo

            src_idx = rays(i)%srcid
            rev_idx = rays(i)%revid
            residual_factor =(time_out(rev_idx,src_idx)-time_data(rev_idx,src_idx))/time_noise(rev_idx,src_idx)**2  
            grads = grads + one_grads * residual_factor
        enddo
        !$omp enddo
        !$omp end parallel
        grads(1:size(parameters)) = -grads(1:size(parameters))/parameters**2

        call kdtree2_destroy(tree)

    endsubroutine derivative

    integer function sites_locate(tree,sites_id,grid,point)
        use kdtree2_module, only    : kdtree2, kdtree2_result, kdtree2_create,&
                                  kdtree2_n_nearest
        implicit none
        type(kdtree2), pointer, intent(in) :: tree
        integer, dimension(:,:), intent(in) :: sites_id
        type(T_GRID), intent(in) :: grid
        real(kind=c_double), dimension(:), intent(in) :: point

        type(kdtree2_result), dimension(1) :: results
        integer ix, iy

        call point2idx(grid,point,ix,iy)

        if(sites_id(iy,ix)==sites_id(iy+1,ix) .and. sites_id(iy,ix)==sites_id(iy,ix+1)&
            .and. sites_id(iy,ix)==sites_id(iy+1,ix+1))then
            sites_locate = sites_id(iy,ix)
        else
            call kdtree2_n_nearest(tp=tree,qv=point,nn=1,results=results)
            sites_locate = results(1)%idx
        endif
    end function

    pure subroutine point2idx(grid,point,ix,iy)
        implicit none
        type(T_GRID), intent(in) :: grid
        real(kind=c_double), dimension(:), intent(in) :: point
        integer(c_int), intent(out) :: ix, iy

        ix = floor((point(1)-grid%xmin)/grid%dx) + 1
        iy = floor((point(2)-grid%ymin/grid%scaling)/(grid%dy/grid%scaling)) + 1
        if(ix < 1) ix = 1
        if(iy < 1) iy = 1
        if(ix >= grid%nx) ix = grid%nx - 1
        if(iy >= grid%ny) iy = grid%ny - 1

    endsubroutine point2idx

    subroutine setup_straightRays(dat, grid, rays, dist)
        implicit none
        type(T_DATA), intent(in) :: dat
        type(T_GRID), intent(in) :: grid
        type(T_RAY), dimension(:,:), intent(inout) :: rays
        real(kind=ii10), dimension(:,:), intent(inout) :: dist

        integer i, j, k, l, n
        real(kind=ii10) :: dl, dx, dy, ds
        real(kind=ii10), dimension(4) :: bnd

        bnd = [grid%xmin,grid%xmax,grid%ymin/grid%scaling,grid%ymax/grid%scaling]
        
        dl = min(grid%dx,grid%dy/grid%scaling)/2
        do i = 1, dat%np
            do j = 1, dat%nsrc
                do k = 1, dat%nrev
                    n = k + (j-1)*dat%nrev
                    dx = dat%rev(1,k)-dat%src(1,j)
                    dy = dat%rev(2,k)-dat%src(2,j)
                    ds = sqrt(dx**2+dy**2)
                    dist(n,i) = ds
                    ! set up straight rays
                    rays(n,i)%srcid = j
                    rays(n,i)%revid = k
                    if( dat%raystat(n,1,i)==1 )then
                        rays(n,i)%npoints = floor(dist(n,i)/dl) + 1
                        allocate( rays(n,i)%points(2,rays(n,i)%npoints) )
                        do l = 1, rays(n,i)%npoints-1
                            rays(n,i)%points(1,l) = dat%src(1,j) + (l-1)*dl*dx/ds
                            rays(n,i)%points(2,l) = dat%src(2,j) + (l-1)*dl*dy/ds
                            if(out_bnd(rays(n,i)%points(:,l),bnd)) &
                                call exception_raiseError('Straight ray out of area') 
                        enddo
                        rays(n,i)%points(:,rays(n,i)%npoints) = dat%rev(:,k)
                    endif
                enddo
            enddo
        enddo

    end subroutine

    subroutine IntegrateRays(vel,grid,rays,time)
        implicit none
        real(kind=ii10), dimension(:,:), intent(in) :: vel
        type(T_GRID), intent(in) :: grid
        type(T_RAY), dimension(:), intent(in) :: rays
        real(kind=ii10), dimension(:,:), intent(inout) :: time

        ! local 
        integer nrays
        integer i, j, k, n
        real(kind=ii10) vhead, vtail, dist

        time = 0
        nrays = 0
        !$omp parallel
        !$omp do private(nrays,vhead,dist,vtail)
        do j = 1, size(time,2)
            do k = 1, size(time,1)
                nrays = nrays + 1
                if(rays(nrays)%npoints>=2)then
                    ! calculate velocity for the first point
                    vhead=GetVelocity(vel,grid,rays(nrays)%points(:,1))
                    do n = 2, rays(nrays)%npoints
                        dist = (rays(nrays)%points(1,n)-rays(nrays)%points(1,n-1))**2 + & 
                               (rays(nrays)%points(2,n)-rays(nrays)%points(2,n-1))**2
                        dist = sqrt(dist)
                        vtail = GetVelocity(vel,grid,rays(nrays)%points(:,n))
                        time(k,j) = time(k,j) + dist*2/(vhead+vtail)
                        vhead = vtail
                    enddo
                endif
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine

    function GetVelocity(vel,grid,point) result(qv)
        implicit none
        real(kind=ii10), dimension(:,:), intent(in) :: vel
        type(T_GRID), intent(in) :: grid
        real(kind=ii10), dimension(:), intent(in) :: point
        real(kind=ii10) :: qv

        integer ix, iy, i, j
        real(kind=ii10) dsx, dsy
        real(kind=ii10) weight

        ix = floor((point(1)-grid%xmin)/grid%dx) + 1
        iy = floor((point(2)-grid%ymin/grid%scaling)/(grid%dy/grid%scaling)) + 1
        if(ix < 1) ix = 1
        if(iy < 1) iy = 1
        if(ix >= grid%nx) ix = grid%nx - 1
        if(iy >= grid%ny) iy = grid%ny - 1
        dsx = point(1) - (grid%xmin + (ix-1)*grid%dx)
        dsy = point(2) - (grid%ymin/grid%scaling + (iy-1)*grid%dy/grid%scaling)

        qv = 0
        do i = 1, 2
            do j = 1, 2
                weight = (1.0-abs((i-1)*grid%dx-dsx)/grid%dx)*(1.0-abs((j-1)*grid%dy/grid%scaling-dsy)/(grid%dy/grid%scaling))
                qv = qv + weight*vel(iy+j-1,ix+i-1)
            enddo
        enddo
    end function

end module m_body_likelihood
