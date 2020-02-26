
! module for modelling travel time data

module m_modelling 

    use iso_c_binding
    use m_settings, only : T_MOD, T_GRID
    use m_utils, only : ii10, vs2vp, vp2rho
    use m_fm2d, only  : T_RAY

    implicit none
    private


    public :: convert_to_layer, T_LAY_MOD, CalGroupTime

    ! static value
    real( kind=ii10 ), parameter, public :: ERAD = 6371
    ! maximum layers
    integer, parameter :: NMAX = 100
    ! water P wave velocity and density
    real( kind=ii10 ), parameter :: waterVel = 1.5
    real( kind=ii10 ), parameter :: waterDensity = 1
    ! vp/vs ratio, currently, used fixed value 1.73
    real( kind=ii10 ), parameter :: POISSON = 1.730
    real( kind=ii10 ), parameter :: PI2 = 6.283185
    real( kind=ii10 ), parameter :: EPS = 1.0E-5

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
                if( grid%waterDepth>0 )then
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
                layer(j,i)%thick(1:nlayers) = thick(1:nlayers)
                !layer(j,i)%ax = 1
                !layer(j,i)%ap = 1
                !layer(j,i)%ae = 1
            enddo
        enddo
        !$omp end do
        !$omp end parallel

        return

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

end module m_modelling

! modelling body and surface wave travel time data for a 3D model

program modelling

    use iso_c_binding
    use fastmarching, only : fastmarching3d
    use m_utils, only : integer_format, double_format, vp2rho_3d, ii10
    use m_settings, only : T_GRID, T_MOD
    use m_modelling, only: T_LAY_MOD, convert_to_layer, ERAD, CalGroupTime
    use m_fm2d, only      : modrays, T_RAY
    use m_surfmodes, only : surfmodes, T_MODES_PARA
    use mt19937, only : gasdev

    use omp_lib

    implicit none

    ! static value
    integer, parameter :: file_unit = 8

    ! grid infomation
    integer(c_int) :: nx, ny, nz
    real(c_double) :: dx, dy, dz
    real(c_double) :: xmin, ymin, zmin
    type(T_GRID) :: grid

    ! data info
    real(c_double), dimension(:,:,:), allocatable :: vel, velg, vp, vs
    type(T_MOD) :: model
    real(c_double), dimension(:,:), allocatable :: src
    real(c_double), dimension(:,:), allocatable :: rev
    real(c_double), dimension(:,:,:), allocatable :: ttime
    real(c_double), dimension(:,:,:), allocatable :: tfield
    integer(c_int) :: nsrc, nrev, ndim

    ! >surface wave related
    ! periods/frequencies
    integer np
    real(c_double), dimension(:), allocatable :: freqs

    ! phase velocities
    real(c_double), dimension(:,:,:), allocatable :: pvel, gvel
    type(T_LAY_MOD), dimension(:,:), allocatable    :: layer
    integer, dimension(:,:), allocatable            :: ierr
    integer ix0, ix1, iy0, iy1
    type(T_MODES_PARA) :: paras

    ! > local variable for fast marching
    integer gridx, gridy
    integer sgref
    integer sgdic, sgext
    integer order
    integer uar
    real( kind=ii10 ) band
    
    ! > local variable for travel times and rays
    integer :: npairs
    integer, dimension(:,:,:), allocatable :: raystat
    integer, dimension(:), allocatable :: crazyray
    type(T_RAY), dimension(:,:), allocatable           :: phaseRays

    ! time
    real(c_double) :: t1, t2
    integer i, j, k

    ! noise level
    real(c_double), parameter :: noise=0.01

    nx = 59
    ny = 71
    nz = 61
    
    allocate(vel(nz,ny,nx*2))
    !vel = 1.0
    open(unit=file_unit, file='model.dat',status='old',action='read')
    read(file_unit,*) vel
    close(file_unit)
    vp = vel(:,:,1:nx)
    vs = vel(:,:,nx+1:2*nx)

    xmin = -3.0
    ymin = -3.5
    zmin = 0
    nsrc = 10
    nrev = 10
    open(unit=file_unit, file='bsources.dat',status='old',action='read')
    read(file_unit,*) nsrc, ndim
    allocate(src(ndim,nsrc))
    read(file_unit,*) src
    close(file_unit)
    open(unit=file_unit, file='breceivers.dat',status='old',action='read')
    read(file_unit,*) nrev, ndim
    allocate(rev(ndim,nrev))
    read(file_unit,*) rev
    close(file_unit)
    write(*,*) 'Number of sources and receivers: ', nsrc, nrev
    write(*,*) 'sources: ', src(:,2)
    write(*,*) 'receivers: ', rev(:,2)
    !call random_number(src)
    !call random_number(rev)
    
    ! body wave first
    dx = 0.1
    dy = 0.1
    dz = 0.05

    ! print some information
    write(*,*) 'xmin: , ymin: , zmin:', xmin, ymin, zmin
    write(*,*) 'nx: , ny: , nz:', nx, ny, nz
    write(*,*) 'dx: , dy: , dz:', dx, dy, dz

    allocate(ttime(nrev,nsrc,2))
    ttime = 0.0
    allocate(tfield(nz,ny,nx))
    tfield =  0.0
    call cpu_time(t1)
    t1 = omp_get_wtime()
    !$omp parallel
    !$omp do 
    do i = 1, nsrc
        call fastmarching3d(vp,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,src(1:3,i),nrev,&
            rev, 1, ttime(:,i,1))
        call fastmarching3d(vs,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,src(1:3,i),nrev,&
            rev, 1, ttime(:,i,2))
        !call fastmarching3d(vp,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,rev(:,i),nsrc,&
        !    src, 1, ttime(i,:,1))
        !call fastmarching3d(vs,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,rev(:,i),nsrc,&
        !    src, 1, ttime(i,:,2))
    enddo
    !$omp enddo
    !$omp end parallel
    call cpu_time(t2)
    t2 = omp_get_wtime()

    write(*,*) 'Time consuming: ', t2-t1

    open(unit=file_unit,file='btimes.dat',action='write')
    write(file_unit,integer_format(1)) 2
    write(file_unit,double_format(2)) 0.0, 1.0
    do j = 1, nsrc
        do k = 1, nrev
            write(file_unit,integer_format(1)) 1
            do i = 1, 2
                write(file_unit,double_format(2)) ttime(k,j,i)+src(4,j)+noise*ttime(k,j,i)*gasdev(),&
                noise*ttime(k,j,i)
            enddo
        enddo
    enddo
    close(file_unit)

    ! Next, surface wave

    ! read periods/frequencies
    open(unit=file_unit,file='freqs.dat',action='read')
    read(file_unit,*) np
    allocate(freqs(np))
    read(file_unit,*) freqs
    freqs = 1/freqs
    ! first set up the parameters for surface modes code
    grid%nx = nx
    grid%ny = ny
    grid%nz = nz
    grid%xmin = xmin
    grid%ymin = ymin
    grid%zmin = zmin
    grid%dx = dx
    grid%dy = dy
    grid%dz = dz
    allocate( model%vp(grid%nz, grid%ny, grid%nx) )
    allocate( model%vs(grid%nz, grid%ny, grid%nx) )
    allocate( model%rho(grid%nz, grid%ny, grid%nx) )
    model%vp = vp
    model%vs = vs
    call vp2rho_3d(vp,model%rho)
    call convert_to_layer( model, grid, 1, nx, 1, ny, layer )

    paras%modetype = 1
    paras%phaseGroup = 1
    paras%tolmin = 1E-6
    paras%tolmax = 1e-7
    paras%smin_min = 1E-3
    paras%smin_max = 5E-3
    paras%dc = 1E-3
    paras%dcm = 1E-3
    paras%dc1 = 1E-3
    paras%dc2 = 1E-3
    iy0 = 1
    iy1 = ny
    ix0 = 1
    ix1 = nx
    ! calculate dispersion curve
    ! TODO: currently, discard any velocity model which could not produce
    ! surface waves in one or more frequencies
    allocate( pvel(np, iy0:iy1, ix0:ix1) )
    allocate( gvel(np, iy0:iy1, ix0:ix1) )
    pvel = 1000.0 ! safe
    gvel = 1000.0 ! safe
    allocate( ierr(iy0:iy1,ix0:ix1) )
    ierr = 0
    t1 = omp_get_wtime ( )
    !call omp_set_num_threads(settings%nthreads)
    !$omp parallel
    !$omp do private(i,j)
    do i = ix0, ix1
        do j = iy0, iy1
            call surfmodes(layer(j,i)%thick,layer(j,i)%alpha,layer(j,i)%beta, &
                layer(j,i)%rho,freqs,paras,pvel(:,j,i),gvel(:,j,i),ierr(j,i))
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    t2 = omp_get_wtime ( )
    write(*,*) 'parallelized dispersion curve code: ', t2-t1
    !call write_vel(pvel,'phaseVel.dat')
    if(any(ierr==1))then
        write(*,*) 'Disperison curve calculation error!'
        stop
    endif


    ! calculate travel time of rayleigh/love wave using fast marching code
    ! settings
    gridx = 1
    gridy = 1
    sgref = 1
    sgdic = 4
    sgext = 8
    order = 1
    band  = 1
    uar   = 0


    ! prepare velocity model for the fast marching code
    if(allocated(vel)) deallocate(vel)
    if(allocated(velg)) deallocate(velg)
    allocate(vel(np, ny+2,nx+2))
    allocate(velg(np, ny,nx))
    vel(:, iy0+1:iy1+1, ix0+1:ix1+1) = pvel
    ! assign boundary value
    if(ix0==1) vel(:,:,1) = vel(:,:,2)
    if(ix1==grid%nx) vel(:,:,grid%nx+2) = vel(:,:,grid%nx+1)
    if(iy0==1) vel(:,1,:) = vel(:,2,:)
    if(iy1==grid%ny) vel(:,grid%ny+2,:) = vel(:,grid%ny+1,:)
    velg(:,iy0:iy1, ix0:ix1) = gvel

    ! prepare raystat
    nsrc = nrev
    allocate(raystat(nrev*nsrc,2,np))
    raystat = 0
    npairs = 0
    do i = 1, nsrc
        do j = 1, nrev
            npairs = npairs + 1
            if(j>i)then
                raystat(npairs,1,:) = 1
                raystat(npairs,2,:) = npairs
            endif
        enddo
    enddo

    ! prepare the beginning time
    if(allocated(ttime)) deallocate(ttime)
    allocate(ttime(nrev,nsrc,np))
    allocate( crazyray(np) )
    allocate( phaseRays(nrev*nsrc,np) )
    crazyray = 0
    t1 = omp_get_wtime ( )
    !$omp parallel
    !$omp do private(i)
    do i = 1, np 
        call modrays(nrev,rev(1,:),rev(2,:), &
                nrev,rev(1,:),rev(2,:), &
            raystat(:,:,i),0, &
            grid%nx,grid%ny,grid%xmin,grid%ymin,&
            grid%dx,grid%dy,vel(i,:,:), &
            gridx,gridy,sgref, &
            sgdic,sgext,ERAD, &
            order,band,ttime(:,:,i), &
            phaseRays(:,i),crazyray(i),uar)
    enddo
    !$omp end do
    !$omp end parallel
    t2 = omp_get_wtime ( )
    write(*,*) 'parallelized fast marching code: ', t2-t1

    ! write the travel times
    open(unit=file_unit,file='stimes.dat',action='write')
    write(file_unit,integer_format(1)) np
    write(file_unit,double_format(np)) freqs
    do j = 1, nsrc
        do k = 1, nrev
            if(any(raystat(k+(j-1)*nrev,1,:)==1))then
                write(file_unit,integer_format(1)) 1
                do i = 1, np
                    write(file_unit,double_format(2)) ttime(k,j,i)+noise*ttime(k,j,i)*gasdev(),&
                    noise*ttime(k,j,i)
                enddo
            else
                write(file_unit,integer_format(1)) 0
            endif
        enddo
    enddo
    close(file_unit)

    ! calculate group travel time for surface waves
    call CalGroupTime(velg,grid,phaseRays,ttime)
    ! write the travel times
    open(unit=file_unit,file='gtimes.dat',action='write')
    write(file_unit,integer_format(1)) np
    write(file_unit,double_format(np)) freqs
    do j = 1, nsrc
        do k = 1, nrev
            if(any(raystat(k+(j-1)*nrev,1,:)==1))then
                write(file_unit,integer_format(1)) 1
                do i = 1, np
                    write(file_unit,double_format(2)) ttime(k,j,i)+noise*ttime(k,j,i)*gasdev(),&
                    noise*ttime(k,j,i)
                enddo
            else
                write(file_unit,integer_format(1)) 0
            endif
        enddo
    enddo
    close(file_unit)

end program
