!
! module m_likelihood
! likelihood function here
!
module m_surf_likelihood

    use m_exception, only : exception_raiseWarning, exception_raiseError
    use m_logger, only    : log_msg
    use m_utils, only     : ii10, write_resume_unit, write_doubles, itoa, rtoa, vs2vp, vp2rho
    use m_settings, only  : T_GRID, T_MOD, out_bnd
    use like_settings, only: T_LIKE_BASE, T_LIKE_SET, T_DATA
    use run_info, only    : T_RUN_INFO
    use m_fm2d, only      : modrays, T_RAY
    use m_surfmodes, only : surfmodes, T_MODES_PARA

    use omp_lib

    implicit none
    private

    public :: surf_likelihood, surf_noise_likelihood
    public :: surf_likelihood_grads

    ! debug
    real( kind=ii10 ) :: t1, t2
    ! static value
    real( kind=ii10 ), parameter :: ERAD = 6371
    ! maximum layers
    integer, parameter :: NMAX = 100
    ! water P wave velocity and density
    real( kind=ii10 ), parameter :: waterVel = 1.5
    real( kind=ii10 ), parameter :: waterDensity = 1
    ! vp/vs ratio, currently, used fixed value 1.73
    real( kind=ii10 ), parameter :: POISSON = 1.730
    real( kind=ii10 ), parameter :: PI2 = 6.283185
    real( kind=ii10 ), parameter :: EPS = 1.0E-10
    ! temporarily use
    integer, parameter :: expand = 1
    integer, parameter :: ext = 2

    type T_LAY_MOD
        integer nlayers
        real(kind=ii10), dimension(:), allocatable :: alpha
        real(kind=ii10), dimension(:), allocatable :: beta
        real(kind=ii10), dimension(:), allocatable :: vpvs
        real(kind=ii10), dimension(:), allocatable :: rho
        real(kind=ii10), dimension(:), allocatable :: thick
        real(kind=ii10), dimension(:), allocatable :: ax
        real(kind=ii10), dimension(:), allocatable :: ap 
        real(kind=ii10), dimension(:), allocatable :: ae
    end type T_LAY_MOD

contains

    subroutine surf_noise_likelihood(dat,RTI,like)
        implicit none
        type(T_DATA), intent(in)                 :: dat
        type(T_RUN_INFO), intent(in)             :: RTI
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
        do i = 1, dat%np
            nrr  = 0
            do j = 1, nsrc
                   do k = 1, nrev
                    nrr = nrr + 1
                    if(dat%raystat(nrr,1,i) == 1) then
                        if(like%sigma(nrr,i)==0) call exception_raiseError('The noise level is 0!')
                        like%like = like%like + ( like%phaseTime(k,j,i)-dat%ttime(nrr,1,i)&
                            )**2/( 2*(like%sigma(nrr,i))**2 )
                        like%misfit = like%misfit + ( like%phaseTime(k,j,i)-dat%ttime(nrr,1,i)&
                            )**2/(like%sigma(nrr,i)**2)
                    else
                        like%sigma(nrr,i) = 1.0
                    endif
                enddo
            enddo
        enddo

        if(any(like%sigma<EPS))then
            call exception_raiseError('The noise level is 0!')
        endif

        like%like = like%like + sum(log(like%sigma)) + dat%nrays/2.0 * log(PI2)
        return
    end subroutine surf_noise_likelihood

    subroutine surf_likelihood_grads(dat,RTI,settings,like)
        implicit none
        type(T_DATA), intent(in)                    :: dat
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(T_LIKE_SET), intent(in)                :: settings
        type(T_LIKE_BASE), intent(inout)            :: like

    endsubroutine surf_likelihood_grads

    subroutine surf_likelihood(dat,RTI,settings,like)
    
        implicit none
        type(T_DATA), intent(in)                    :: dat
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(T_LIKE_SET), intent(in)                :: settings
        type(T_LIKE_BASE), intent(inout)            :: like

        ! > local variable for normal mode modeling
        type(T_MODES_PARA) :: paras
        integer         :: ix0, ix1
        integer         :: idx0, idx1
        !integer         :: iper, ier, ia
        !real(kind=ii10) :: w, ekd, y0l(6), y0r(3), yij(15)
        integer, dimension(:), allocatable            :: ierr
        real(kind=ii10), dimension(:,:), allocatable :: pvel, gvel
        type(T_LAY_MOD), dimension(:), allocatable    :: layer

        ! > local variable for travel times and rays
        integer nsrc, nrev

        ! > local grid
        type(T_GRID) :: grid
        type(T_MOD) :: model

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
                model%vp(j,i) = RTI%parameters(1,RTI%sites_id(j,i))
                model%vs(j,i) = RTI%parameters(2,RTI%sites_id(j,i))
                model%rho(j,i) = RTI%parameters(3,RTI%sites_id(j,i))
            enddo
        enddo

        ! first, calculate dispersion curve grid by grid
        ! convert 3d grid vp, vs and rho model to layered model
        ix0 = 1
        ix1 = grid%nx
        call convert_to_layer( model, grid, ix0, ix1, layer )

        ! first set up the parameters for surface modes code
        paras%modetype = settings%raylov
        paras%phaseGroup = settings%phaseGroup
        paras%tolmin = settings%tol
        paras%tolmax = 10*settings%tol
        paras%smin_min = 1E-3
        paras%smin_max = 5E-3
        paras%dc = settings%dPhaseVel
        paras%dcm = settings%dPhaseVel
        paras%dc1 = settings%dPhaseVel
        paras%dc2 = settings%dPhaseVel
        ! calculate dispersion curve
        ! TODO: currently, discard any velocity model which could not produce
        ! surface waves in one or more frequencies
        allocate( pvel(dat%np, ix0:ix1) )
        allocate( gvel(dat%np, ix0:ix1) )
        pvel = 100.0 ! safe
        gvel = 100.0 ! safe
        allocate( ierr(ix0:ix1) )
        ierr = 0
#ifdef _OPENMP
        t1 = omp_get_wtime ( )
        !call log_msg('Begin normal modes...')
        !call omp_set_num_threads(settings%nthreads)
#endif
        !$omp parallel
        !$omp do private(i,j)
        do i = ix0, ix1
           call surfmodes(layer(i)%thick,layer(i)%alpha,layer(i)%beta, &
               layer(i)%rho,dat%freqs,paras,pvel(:,i),gvel(:,i),ierr(i))
        enddo
        !$omp end do
        !$omp end parallel
        t2 = omp_get_wtime ( )
        !call log_msg(itoa(omp_get_num_threads() ) )
        !call log_msg('parallelized dispersion curve code: ' //rtoa(t2-t1) )
        !call write_vel(pvel,'phaseVel.dat')
        ! discard models
        if(any(ierr==1))then
            like%like = huge(like%like)
            RTI%num_bad_model = RTI%num_bad_model + 1
            return
        endif

        ! if using straight rays
        if(.not.like%straightRaySet)then
            allocate( like%rays(dat%nsrc*dat%nrev,1) )
            like%rays%srcid = 0
            like%rays%revid = 0
            like%rays%npoints = 0
            call setup_straightRays(dat, grid, like%rays,like%srdist)
            like%straightRaySet = .true.
        endif
        if(settings%phaseGroup==0)then
            call CalGroupTime(pvel,grid,like%rays,like%phaseTime)
        else
            call CalGroupTime(gvel,grid,like%rays,like%phaseTime)
        endif

        ! likelihood
        ! noise level
        if(settings%sigdep /= 0)then
            do i = 1, dat%np
                nrr  = 0
                do j = 1, nsrc
                    do k = 1, nrev
                        nrr = nrr + 1
                        if(dat%raystat(nrr,1,i) == 1) then
                            like%sigma(nrr,i) = RTI%snoise0(i)*like%srdist(nrr,i) + RTI%snoise1(i) 
                        else
                            like%sigma(nrr,i) = 1.0
                        endif
                    enddo
                enddo
            enddo
        else
            like%sigma = dat%ttime(:,2,:)
        endif

        !if(any(like%sigma==0))then
        !    call exception_raiseError('The noise level is 0!')
        !endif

        like%like = 0
        like%misfit = 0
        like%unweighted_misfit = 0
        do i = 1, dat%np
            nrr  = 0
            do j = 1, nsrc
                   do k = 1, nrev
                    nrr = nrr + 1
                    if(dat%raystat(nrr,1,i) == 1) then
                        if(like%sigma(nrr,i)<EPS)then
                            !print * , k, j, i, like%sigma(nrr,i), like%srdist(nrr,i)
                            call exception_raiseError('The noise level is 0!')
                        endif
                        like%like = like%like + ( like%phaseTime(k,j,i)-dat%ttime(nrr,1,i)&
                            )**2/( 2*(like%sigma(nrr,i))**2 )
                        like%misfit = like%misfit + ( like%phaseTime(k,j,i)-dat%ttime(nrr,1,i)&
                            )**2/(like%sigma(nrr,i)**2)
                        like%unweighted_misfit = like%unweighted_misfit + ( like%phaseTime(k,j,i)-dat%ttime(nrr,1,i)&
                            )**2
                    else
                        like%sigma(nrr,i) = 1.0
                    endif
                enddo
            enddo
        enddo

        like%like = like%like + sum(log(like%sigma)) + dat%nrays/2.0 * log(PI2)
        if(like%like/=like%like)then
            print *, like%sigma
        endif
        return

    end subroutine surf_likelihood

    subroutine setup_straightRays(dat, grid, rays, dist)
        implicit none
        type(T_DATA), intent(in) :: dat
        type(T_GRID), intent(in) :: grid
        type(T_RAY), dimension(:,:), intent(inout) :: rays
        real(kind=ii10), dimension(:,:), intent(inout) :: dist

        integer i, j, k, l, n
        real(kind=ii10) :: dl, dx, dy, ds
        real(kind=ii10), dimension(4) :: bnd

        bnd = [grid%xmin,grid%xmax,0.0_ii10,0.0_ii10]
        
        dl = grid%dx/2
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
                        allocate( rays(n,i)%points(rays(n,i)%npoints,2) )
                        do l = 1, rays(n,i)%npoints-1
                            rays(n,i)%points(l,1) = dat%src(1,j) + (l-1)*dl*dx/ds
                            rays(n,i)%points(l,2) = dat%src(2,j) + (l-1)*dl*dy/ds
                            if(out_bnd(rays(n,i)%points(l,:),bnd)) &
                                call exception_raiseError('Straight ray out of area') 
                        enddo
                        rays(n,i)%points(rays(n,i)%npoints,:) = dat%rev(:,k)
                    endif
                enddo
            enddo
        enddo

    end subroutine

    subroutine CalGroupTime(vel,grid,rays,time)
        implicit none
        real(kind=ii10), dimension(:,:), intent(in) :: vel
        type(T_GRID), intent(in) :: grid
        type(T_RAY), dimension(:,:), intent(in) :: rays
        real(kind=ii10), dimension(:,:,:), intent(inout) :: time

        ! local 
        integer nrays
        integer i, j, k, n
        real(kind=ii10) vhead, vtail, dist

        nrays = 0
        time = 0
        !$omp parallel
        !$omp do private(nrays,vhead,dist,vtail)
        do i = 1, size(time,3)
            nrays = 0
            do j = 1, size(time,2)
                do k = 1, size(time,1)
                    nrays = nrays + 1
                    if(rays(nrays,i)%npoints>=2)then
                        ! calculate velocity for the first point
                        vhead=GetVelocity(vel(i,:),grid,rays(nrays,i)%points(1,:))
                        do n = 2, rays(nrays,i)%npoints
                            dist = (rays(nrays,i)%points(n,1)-rays(nrays,i)%points(n-1,1))**2 + & 
                                   (rays(nrays,i)%points(n,2)-rays(nrays,i)%points(n-1,2))**2
                            dist = sqrt(dist)
                            vtail = GetVelocity(vel(i,:),grid,rays(nrays,i)%points(n,:))
                            time(k,j,i) = time(k,j,i) + dist*2/(vhead+vtail)
                            vhead = vtail
                        enddo
                    endif
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine

    function GetVelocity(vel,grid,point) result(qv)
        implicit none
        real(kind=ii10), dimension(:), intent(in) :: vel
        type(T_GRID), intent(in) :: grid
        real(kind=ii10), dimension(:), intent(in) :: point
        real(kind=ii10) :: qv

        integer ix, iy, i, j
        real(kind=ii10) dsx, dsy
        real(kind=ii10) weight

        ix = floor((point(1)-grid%xmin)/grid%dx) + 1
        if(ix < 1) ix = 1
        if(ix >= grid%nx) ix = grid%nx - 1
        iy = 1
        dsx = point(1) - (grid%xmin + (ix-1)*grid%dx)

        qv = 0
        do i = 1, 2
            weight = (1.0-abs((i-1)*grid%dx-dsx)/grid%dx)
            qv = qv + weight*vel(ix+i-1)
        enddo
    end function

    subroutine convert_to_layer(model, grid, ix0, ix1, layer)
        implicit none

        type(T_MOD), intent(in) :: model
        type(T_GRID), intent(in) :: grid
        integer, intent(in) :: ix0, ix1
        type(T_LAY_MOD), dimension(:),allocatable, intent(out) :: layer

        ! > local variables
        integer i, k
        integer nlayers
        integer last_k
        real(kind=ii10) :: last_vp, last_vs, last_rho
        real(kind=ii10), dimension(NMAX) :: alpha, beta, rho_k, thick

        allocate( layer(ix0:ix1) )
        !$omp parallel
        !$omp do private(i,k,last_vp,last_vs,last_rho,last_k,nlayers,alpha)&
        !$omp& private(beta,rho_k,thick)
        do i = ix0, ix1
            if( grid%waterDepth>EPS )then
                nlayers = 1
                alpha(1) = waterVel
                beta(1) = 0
                rho_k(1) = waterDensity
                thick(1) = grid%waterDepth
            else
                nlayers = 0
            endif
            last_vp = model%vp(1,i)
            last_vs = model%vs(1,i)
            last_rho = model%rho(1,i)
            last_k = 1
            do k = 2, grid%ny
                if( abs(model%vs(k,i)-last_vs) > EPS ) then
                    nlayers = nlayers + 1
                    ! store the parameters value for this layer
                    alpha(nlayers) = last_vp
                    beta(nlayers) = last_vs
                    rho_k(nlayers) = last_rho
                    thick(nlayers) = (k-last_k) * grid%dz
                    ! update last_* values
                    last_vp = model%vp(k,i)
                    last_vs = model%vs(k,i)
                    last_rho = model%rho(k,i)
                    last_k = k
                endif
            enddo
            ! last 2 layers, last thick is 0
            nlayers = nlayers + 1
            if( last_k /= grid%ny) then
                alpha(nlayers) = model%vp(grid%ny,i)
                beta(nlayers) = model%vs(grid%ny,i)
                rho_k(nlayers) = model%rho(grid%ny,i)
                thick(nlayers) = (grid%ny-last_k) * grid%dz
                thick(nlayers) = 0
            else
                alpha(nlayers) = model%vp(grid%ny,i)
                beta(nlayers) = model%vs(grid%ny,i)
                rho_k(nlayers) = model%rho(grid%ny,i)
                thick(nlayers) = 0
            endif
            ! allocate and give value for layer model
            allocate( layer(i)%alpha(nlayers), layer(i)%beta(nlayers),&
            layer(i)%vpvs(nlayers), layer(i)%rho(nlayers),&
            layer(i)%thick(nlayers) )
            !allocate( layer(j,i)%ax(nlayers), layer(j,i)%ap(nlayers),&
            !layer(j,i)%ae(nlayers) )
            layer(i)%nlayers = nlayers
            layer(i)%alpha(1:nlayers) = alpha(1:nlayers)
            layer(i)%beta(1:nlayers) = beta(1:nlayers)
            layer(i)%vpvs = POISSON
            layer(i)%rho(1:nlayers) = rho_k(1:nlayers)
            layer(i)%thick(1:nlayers) = thick(1:nlayers)/grid%scaling
            ! waterDepth deos not scale, need to be recovered
            if(grid%waterDepth>0) layer(i)%thick(1)=grid%waterDepth

            !layer(j,i)%ax = 1
            !layer(j,i)%ap = 1
            !layer(j,i)%ae = 1
        enddo
        !$omp end do
        !$omp end parallel

        return

    end subroutine


    subroutine write_layer(layer,filename)
        use m_utils, only : write_resume_unit
        implicit none
        character(len=*), intent(in) :: filename
        type(T_LAY_MOD), intent(in) :: layer

        integer i

        open(unit=write_resume_unit, file=filename, status='unknown')
        write(write_resume_unit,*) layer%nlayers
        do i = 1, layer%nlayers
            write(write_resume_unit,*) i,layer%thick(i)*1000,layer%beta(i)*1000,&
            layer%vpvs(i), layer%rho(i)
        enddo
        close(write_resume_unit)
        return
    end subroutine

    subroutine write_layers(layer,filename)
        use m_utils, only : write_resume_unit
        implicit none
        character(len=*), intent(in) :: filename
        type(T_LAY_MOD), dimension(:,:), intent(in) :: layer

        integer i, j, k

        open(unit=write_resume_unit, file=filename, status='unknown')
        do j = lbound(layer,2), ubound(layer,2)
        do k = lbound(layer,1), ubound(layer,1)
            write(write_resume_unit,*) layer(k,j)%nlayers
            do i = 1, layer(k,j)%nlayers
                write(write_resume_unit,*) i,layer(k,j)%thick(i)*1000,layer(k,j)%beta(i)*1000,&
                layer(k,j)%vpvs(i), layer(k,j)%rho(i)
            enddo
        enddo
        enddo
        close(write_resume_unit)
        return
    end subroutine

    subroutine write_vel(pvel,filename)
        use m_utils, only : write_resume_unit, write_doubles
        implicit none
        character(len=*), intent(in) :: filename
        real(kind=ii10), dimension(:,:,:), intent(in) :: pvel

        open(unit=write_resume_unit, file=filename, status='unknown',&
        position='append')
        call write_doubles(pvel)
        close(write_resume_unit)
        return
    end subroutine
end module m_surf_likelihood
