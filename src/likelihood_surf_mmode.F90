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
    use m_surfmodes, only : surfmmodes, T_MODES_PARA
    use cgal_delaunay,only: d3

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

        integer nrr, nmodes
        integer nsrc, nrev
        integer i, j, k

        nsrc = dat%nsrc
        nrev = dat%nrev
        nmodes = dat%nmodes

        do i = 1, dat%np*nmodes
            like%sigma(:,i) = RTI%snoise0(i)*like%srdist(:,i) + RTI%snoise1(i)
        enddo

        like%like = 0
        like%misfit = 0
        do i = 1, dat%np*nmodes
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

    subroutine surf_likelihood(dat,model,RTI,perturbed_box,settings,like)
    
        implicit none
        type(T_DATA), intent(in)                    :: dat
        type(T_MOD), intent(inout)                  :: model
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(d3), dimension(2), intent(in)          :: perturbed_box
        type(T_LIKE_SET), intent(in)                :: settings
        type(T_LIKE_BASE), intent(inout)            :: like

        ! > local variable for normal mode modeling
        type(T_MODES_PARA) :: paras
        integer         :: ix0, ix1, iy0, iy1
        integer         :: idx0, idx1, idy0, idy1
        !integer         :: iper, ier, ia
        !real(kind=ii10) :: w, ekd, y0l(6), y0r(3), yij(15)
        integer, dimension(:,:), allocatable            :: ierr
        real(kind=ii10), dimension(:,:,:), allocatable :: pvel, gvel
        type(T_LAY_MOD), dimension(:,:), allocatable    :: layer

        ! > local variable for fast marching
        integer gridx, gridy
        integer sgref
        integer sgdic, sgext
        integer order
        integer uar
        real( kind=ii10 ) band
    
        ! > local variable for travel times and rays
        integer nsrc, nrev
        integer, dimension(:), allocatable :: crazyray
        type(T_RAY), dimension(:,:), allocatable           :: phaseRays

        ! > local grid
        type(T_GRID) :: grid

        integer nrr
        integer i, j, k

        nsrc = dat%nsrc
        nrev = dat%nrev
        grid = settings%grid

        ! first, calculate dispersion curve grid by grid
        ! convert 3d grid vp, vs and rho model to layered model
        ix0 = floor( (perturbed_box(1)%x-grid%xmin)/grid%dx ) + 1 - expand
        ix1 = floor( (perturbed_box(2)%x-grid%xmin)/grid%dx ) + 1 + expand
        iy0 = floor( (perturbed_box(1)%y-grid%ymin)/grid%dy ) + 1 - expand
        iy1 = floor( (perturbed_box(2)%y-grid%ymin)/grid%dy ) + 1 + expand

        ! check model, discard those the top layer is not smallest
        if(check_model(model%vs))then
            like%like = huge(like%like)
            return
        endif

        if(ix0<1) ix0 = 1
        if(iy0<1) iy0 = 1
        if(ix1>grid%nx) ix1 = grid%nx
        if(iy1>grid%ny) iy1 = grid%ny
        call convert_to_layer( model, grid, ix0, ix1, iy0, iy1, layer )

        ! first set up the parameters for surface modes code
        paras%modetype = settings%raylov
        paras%nmodes = dat%nmodes
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
        allocate( pvel(dat%np*dat%nmodes, iy0:iy1, ix0:ix1) )
        allocate( gvel(dat%np*dat%nmodes, iy0:iy1, ix0:ix1) )
        pvel = 100.0 ! safe
        gvel = 100.0 ! safe
        allocate( ierr(iy0:iy1,ix0:ix1) )
        ierr = 0
#ifdef _OPENMP
        t1 = omp_get_wtime ( )
        !call log_msg('Begin normal modes...')
        !call omp_set_num_threads(settings%nthreads)
#endif
        !$omp parallel
        !$omp do private(i,j)
        do i = ix0, ix1
            do j = iy0, iy1
                call surfmodes(layer(j,i)%thick,layer(j,i)%alpha,layer(j,i)%beta, &
                    layer(j,i)%rho,dat%freqs,paras,pvel(:,j,i),gvel(:,j,i),ierr(j,i))
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        t2 = omp_get_wtime ( )
        !call log_msg(itoa(omp_get_num_threads() ) )
        !call log_msg('parallelized dispersion curve code: ' //rtoa(t2-t1) )
        !call write_vel(pvel,'phaseVel.dat')
        ! discard models
        if(any(ierr==1))then
            ! prepare velocity model for the fast marching code
            like%vel(:, iy0+1:iy1+1, ix0+1:ix1+1) = pvel
            like%gvel(:,iy0:iy1, ix0:ix1) = gvel
            ! assign boundary value
            if(ix0==1) like%vel(:,:,1) = like%vel(:,:,2)
            if(ix1==grid%nx) like%vel(:,:,grid%nx+2) = like%vel(:,:,grid%nx+1)
            if(iy0==1) like%vel(:,1,:) = like%vel(:,2,:)
            if(iy1==grid%ny) like%vel(:,grid%ny+2,:) = like%vel(:,grid%ny+1,:)
            like%like = huge(like%like)
            RTI%num_bad_model = RTI%num_bad_model + 1
            return
        endif

        ! phase/group velocity
        if(settings%phaseGroup == 1)then
            like%gvel(:,iy0:iy1, ix0:ix1) = gvel
        else
            like%gvel(:,iy0:iy1, ix0:ix1) = pvel
        endif

        ! if using straight rays
        if(settings%isStraight == 1)then
            if(.not.like%straightRaySet)then
                allocate( like%rays(dat%nsrc*dat%nrev, dat%np) )
                like%rays%srcid = 0
                like%rays%revid = 0
                like%rays%npoints = 0
                call setup_straightRays(dat, grid, like%rays,like%srdist)
                like%straightRaySet = .true.
            endif
            call CalGroupTime(like%gvel,grid,like%rays,like%phaseTime)
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
            !uar   = settings%uar
            uar   = abs(settings%phaseGroup-1)


            ! prepare velocity model for the fast marching code
            like%vel(:, iy0+1:iy1+1, ix0+1:ix1+1) = pvel
            ! assign boundary value
            if(ix0==1) like%vel(:,:,1) = like%vel(:,:,2)
            if(ix1==grid%nx) like%vel(:,:,grid%nx+2) = like%vel(:,:,grid%nx+1)
            if(iy0==1) like%vel(:,1,:) = like%vel(:,2,:)
            if(iy1==grid%ny) like%vel(:,grid%ny+2,:) = like%vel(:,grid%ny+1,:)

            ! prepare the beginning time
            idx0 = ix0 - ext
            idx1 = ix1 + ext
            idy0 = iy0 - ext
            idy1 = iy1 + ext
            if(idx0<1) idx0 = 1
            if(idy0<1) idy0 = 1
            if(idx1>grid%nx) idx1 = grid%nx
            if(idy1>grid%ny) idy1 = grid%ny
            if(settings%dynamic == 2)then
                do i = 1, dat%np*dat%nmodes
                    do j = 1, nsrc
                        like%btime(j,i) = &
                        minval(like%field4d(idy0:idy1,idx0:idx1,j,i))
                    enddo
                enddo
            endif

            ! call modrays to calculate travel times
            allocate( phaseRays(dat%nrev*dat%nsrc, dat%np) )
            phaseRays%srcid = 0
            phaseRays%revid = 0
            phaseRays%npoints = 0
            allocate( crazyray(dat%np) )
            crazyray = 0
#ifdef _OPENMP
            t1 = omp_get_wtime ( )
            !call omp_set_num_threads(settings%nthreads)
#endif
            !$omp parallel
            !$omp do private(nrr,i,j,k)
            do i = 1, dat%np*dat%nmodes 
                !if( .not.allocated(ttime) ) allocate( ttime(nrev, nsrc) )
                if(settings%dynamic == 2)then
                    call modrays(nsrc,dat%src(1,:),dat%src(2,:), &
                            nrev,dat%rev(1,:),dat%rev(2,:), &
                        dat%raystat(:,:,i),0, &
                        grid%nx,grid%ny,grid%xmin,grid%ymin,&
                        grid%dx,grid%dy,like%vel(i,:,:), &
                        gridx,gridy,sgref, &
                        sgdic,sgext,ERAD, &
                        order,band,like%phaseTime(:,:,i), &
                        phaseRays(:,i),crazyray(i),uar,&
                        like%field4d(:,:,:,i),like%btime(:,i))
                else
                    call modrays(nsrc,dat%src(1,:),dat%src(2,:), &
                            nrev,dat%rev(1,:),dat%rev(2,:), &
                        dat%raystat(:,:,i),0, &
                        grid%nx,grid%ny,grid%xmin,grid%ymin,&
                        grid%dx,grid%dy,like%vel(i,:,:), &
                        gridx,gridy,sgref, &
                        sgdic,sgext,ERAD, &
                        order,band,like%phaseTime(:,:,i), &
                        phaseRays(:,i),crazyray(i),uar)
                endif

                ! allocate and storage the ray info
                do j = 1, nsrc
                    do k = 1, nrev
                        nrr = k + (j-1)*nrev
                        if(settings%phaseGroup ==1)then
                            like%srdist(nrr,i) = phaseRays(nrr,i)%length()
                        else
                            like%srdist(nrr,i) = like%phaseTime(k,j,i)
                        endif
                    enddo
                enddo
            enddo
            !$omp end do
            !call log_msg(itoa(omp_get_num_threads() ) )
            !$omp end parallel
            t2 = omp_get_wtime ( )
            !call log_msg('parallelized fast marching code: ' //rtoa(t2-t1) )
    
            ! if crazyray, return
            if( any(crazyray > 0) )then
                like%like = huge(like%like)
                RTI%num_bad_ray = RTI%num_bad_ray + 1
                return
            endif

            ! if group
            if(settings%phaseGroup == 1)then
                ! calculate group delay along rays
                call CalGroupTime(like%gvel,grid,phaseRays,like%phaseTime)
            endif
        endif

        ! likelihood
        ! noise level
        if(settings%sigdep /= 0)then
            do i = 1, dat%np*dat%nmodes
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
        do i = 1, dat%np*dat%nmodes
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

        bnd = [grid%xmin,grid%xmax,grid%ymin,grid%ymax]
        
        dl = min(grid%dx,grid%dy)/2
        do i = 1, dat%np*dat%nmodes
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
        real(kind=ii10), dimension(:,:,:), intent(in) :: vel
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
                        vhead=GetVelocity(vel(i,:,:),grid,rays(nrays,i)%points(:,1))
                        do n = 2, rays(nrays,i)%npoints
                            dist = (rays(nrays,i)%points(1,n)-rays(nrays,i)%points(1,n-1))**2 + & 
                                   (rays(nrays,i)%points(2,n)-rays(nrays,i)%points(2,n-1))**2
                            dist = sqrt(dist)
                            vtail = GetVelocity(vel(i,:,:),grid,rays(nrays,i)%points(:,n))
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
        real(kind=ii10), dimension(:,:), intent(in) :: vel
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
                qv = qv + weight*vel(iy+j-1,ix+i-1)
            enddo
        enddo
    end function

    subroutine convert_to_layer(model, grid, ix0, ix1, iy0,&
                                iy1, layer, vbnd)
        implicit none

        type(T_MOD), intent(in) :: model
        type(T_GRID), intent(in) :: grid
        integer, intent(in) :: ix0, ix1, iy0, iy1
        type(T_LAY_MOD), dimension(:,:),allocatable, intent(out) :: layer
        real( kind=ii10), intent(in), optional :: vbnd

        ! > local variables
        integer i, j, k
        integer nlayers
        integer last_k
        real(kind=ii10) :: last_vp, last_vs, last_rho
        real(kind=ii10), dimension(NMAX) :: alpha, beta, rho_k, thick

        allocate( layer(iy0:iy1, ix0:ix1) )
        !$omp parallel
        !$omp do private(i,j,k,last_vp,last_vs,last_rho,last_k,nlayers,alpha)&
        !$omp& private(beta,rho_k,thick)
        do i = ix0, ix1
            do j = iy0, iy1
                if( grid%waterDepth>EPS )then
                    nlayers = 1
                    alpha(1) = waterVel
                    beta(1) = 0
                    rho_k(1) = waterDensity
                    thick(1) = grid%waterDepth
                else
                    nlayers = 0
                endif
                last_vp = model%vp(1,j,i)
                last_vs = model%vs(1,j,i)
                last_rho = model%rho(1,j,i)
                last_k = 1
                do k = 2, grid%nz
                    if( abs(model%vs(k,j,i)-last_vs) > EPS ) then
                        nlayers = nlayers + 1
                        ! store the parameters value for this layer
                        alpha(nlayers) = last_vp
                        beta(nlayers) = last_vs
                        rho_k(nlayers) = last_rho
                        thick(nlayers) = (k-last_k) * grid%dz
                        ! update last_* values
                        last_vp = model%vp(k,j,i)
                        last_vs = model%vs(k,j,i)
                        last_rho = model%rho(k,j,i)
                        last_k = k
                    endif
                enddo
                ! last 2 layers, last thick is 0
                nlayers = nlayers + 1
                if( last_k /= grid%nz) then
                    alpha(nlayers) = model%vp(grid%nz,j,i)
                    beta(nlayers) = model%vs(grid%nz,j,i)
                    rho_k(nlayers) = model%rho(grid%nz,j,i)
                    thick(nlayers) = (grid%nz-last_k) * grid%dz
                    if(present(vbnd)) then
                        nlayers = nlayers + 1
                        !alpha(nlayers) = model%vp(grid%nz,j,i)
                        !beta(nlayers) = model%vs(grid%nz,j,i)
                        !rho_k(nlayers) = model%rho(grid%nz,j,i)
                    !else
                        alpha(nlayers) = vs2vp(vbnd)
                        beta(nlayers) = vbnd
                        rho_k(nlayers) = vp2rho(alpha(nlayers))
                    endif
                    thick(nlayers) = 0
                else
                    if(present(vbnd)) then
                        alpha(nlayers) = vs2vp(vbnd)
                        beta(nlayers) = vbnd
                        rho_k(nlayers) = vp2rho(alpha(nlayers))
                    else
                        alpha(nlayers) = model%vp(grid%nz,j,i)
                        beta(nlayers) = model%vs(grid%nz,j,i)
                        rho_k(nlayers) = model%rho(grid%nz,j,i)
                    endif
                    thick(nlayers) = 0
                endif
                ! allocate and give value for layer model
                allocate( layer(j,i)%alpha(nlayers), layer(j,i)%beta(nlayers),&
                layer(j,i)%vpvs(nlayers), layer(j,i)%rho(nlayers),&
                layer(j,i)%thick(nlayers) )
                !allocate( layer(j,i)%ax(nlayers), layer(j,i)%ap(nlayers),&
                !layer(j,i)%ae(nlayers) )
                layer(j,i)%nlayers = nlayers
                layer(j,i)%alpha(1:nlayers) = alpha(1:nlayers)
                layer(j,i)%beta(1:nlayers) = beta(1:nlayers)
                layer(j,i)%vpvs = POISSON
                layer(j,i)%rho(1:nlayers) = rho_k(1:nlayers)
                layer(j,i)%thick(1:nlayers) = thick(1:nlayers)/grid%scaling
                ! waterDepth deos not scale, need to be recovered
                if(grid%waterDepth>0) layer(j,i)%thick(1)=grid%waterDepth

                !layer(j,i)%ax = 1
                !layer(j,i)%ap = 1
                !layer(j,i)%ae = 1
            enddo
        enddo
        !$omp end do
        !$omp end parallel

        return

    end subroutine

    function check_model(vs) result(valid)
        real(kind=ii10), dimension(:,:,:), intent(in) :: vs
        logical valid
        integer i, j
        
        valid = .false.
        do i = 1, size(vs,3)
            do j = 1, size(vs,2)
                if(any(vs(2:size(vs,1),j,i)<vs(1,j,i)))then
                    valid = .true.
                    return
                endif
            enddo
        enddo
    
    endfunction


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
