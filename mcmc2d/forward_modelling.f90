! modelling 2d surface wave travel time data for a 2D model

program modelling

    use iso_c_binding
    use m_utils, only : integer_format, double_format, vp2rho_3d, ii10
    use m_fm2d, only      : modrays, T_RAY
    use mt19937, only : gasdev

    use omp_lib

    implicit none

    ! static value
    integer, parameter :: file_unit = 8
    real(c_double), parameter :: EARD = 6371.0

    ! grid infomation
    integer(c_int) :: nx, ny, nz
    real(c_double) :: dx, dy, dz
    real(c_double) :: xmin, ymin, zmin

    ! data info
    real(c_double), dimension(:,:), allocatable :: vel, vp, vs
    real(c_double), dimension(:,:), allocatable :: src
    real(c_double), dimension(:,:), allocatable :: rev
    real(c_double), dimension(:,:), allocatable :: ttime
    real(c_double), dimension(:,:,:), allocatable :: tfield
    integer(c_int) :: nsrc, nrev, ndim


    ! > local variable for fast marching
    integer gridx, gridy
    integer sgref
    integer sgdic, sgext
    integer order
    integer uar
    real( kind=ii10 ) band
    
    ! > local variable for travel times and rays
    integer :: ix0, ix1, iy0, iy1
    integer :: npairs
    integer, dimension(:,:), allocatable :: raystat
    integer :: crazyray
    type(T_RAY), dimension(:), allocatable           :: phaseRays
    logical ex
    integer nt

    ! time
    real(c_double) :: t1, t2
    integer i, j, k

    ! noise level
    real(c_double), parameter :: noise=0.02

    nx = 101
    ny = 281
    nz = 51
    
    !vel = 1.0
    open(unit=file_unit, file='model.dat',status='old',action='read')
    read(file_unit,*) ny, nx
    allocate(vel(ny,nx))
    read(file_unit,*) vel
    close(file_unit)
    vp = vel

    xmin = 0
    ymin = 0 
    zmin = 0
    dx = 0.05
    dy = 0.05
    nsrc = 10
    nrev = 20
    open(unit=file_unit, file='sources.dat',status='old',action='read')
    read(file_unit,*) nsrc
    allocate(src(2,nsrc))
    read(file_unit,*) src
    close(file_unit)
    open(unit=file_unit, file='receivers.dat',status='old',action='read')
    read(file_unit,*) nrev
    allocate(rev(2,nrev))
    read(file_unit,*) rev
    close(file_unit)
    write(*,*) 'Number of sources and receivers: ', nsrc, nrev
    write(*,*) 'sources: ', src(:,2)
    write(*,*) 'receivers: ', rev(:,2)
    !call random_number(src)
    !call random_number(rev)
    ! read ray status if applicable
    allocate(raystat(nrev*nsrc,2))
    raystat = 0
    npairs = 0
    inquire(file='otimes_stat.dat',exist=ex)
    if(ex)then
        open(unit=file_unit,file='otimes_stat.dat',action='read')
        read(file_unit,integer_format(1)) nt
        read(file_unit,*)
        do j = 1, nsrc
            do k = 1, nrev
                npairs = npairs + 1
                read(file_unit,integer_format(1)) raystat(npairs,1)
                raystat(npairs,2) = npairs
                if(raystat(npairs,1)==0) cycle
                do i =  1, nt
                    read(file_unit,double_format(2))
                enddo
            enddo
        enddo
        close(file_unit)
    else
        ! prepare raystat
        do i = 1, nsrc
            do j = 1, nrev
                npairs = npairs + 1
                if(any(rev(:,j)/=src(:,i)))then
                    raystat(npairs,1) = 1
                    raystat(npairs,2) = npairs
                endif
            enddo
        enddo
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
    uar   = 1


    ! prepare velocity model for the fast marching code
    if(allocated(vel)) deallocate(vel)
    allocate(vel(ny+2,nx+2))
    iy0 = 1
    iy1 = ny
    ix0 = 1
    ix1 = nx
    vel(iy0+1:iy1+1, ix0+1:ix1+1) = vp
    ! assign boundary value
    if(ix0==1) vel(:,1) = vel(:,2)
    if(ix1==nx) vel(:,nx+2) = vel(:,nx+1)
    if(iy0==1) vel(1,:) = vel(2,:)
    if(iy1==ny) vel(ny+2,:) = vel(ny+1,:)

    ! prepare the beginning time
    if(allocated(ttime)) deallocate(ttime)
    allocate(ttime(nrev,nsrc))
    crazyray = 0
    t1 = omp_get_wtime ( )
    call modrays(nsrc,src(1,:),src(2,:), &
            nrev,rev(1,:),rev(2,:), &
        raystat(:,:),0, &
        nx,ny,xmin,ymin,&
        dx,dy,vel(:,:), &
        gridx,gridy,sgref, &
        sgdic,sgext,EARD, &
        order,band,ttime(:,:), &
        phaseRays,crazyray,uar)
    t2 = omp_get_wtime ( )
    write(*,*) 'parallelized fast marching code: ', t2-t1

    ! write the travel times
    open(unit=file_unit,file='otimes.dat',action='write')
    write(file_unit,integer_format(1)) 1
    write(file_unit,double_format(1)) 0.9
    do j = 1, nsrc
        do k = 1, nrev
            if(raystat(k+(j-1)*nrev,1)==1)then
                write(file_unit,integer_format(1)) 1
                write(file_unit,double_format(2)) ttime(k,j)+noise*ttime(k,j)*gasdev(),noise*ttime(k,j)
            else
                write(file_unit,integer_format(1)) 0
            endif
        enddo
    enddo
    close(file_unit)

end program
