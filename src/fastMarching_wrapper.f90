! wrapper for fast marching code
module fastmarching
    use iso_c_binding

    use m_utils, only : write_doubles, ii10
    
    implicit none

    private

    ! the fortran interface to fastmarching
    public :: fastmarching3d
    public :: T_RAY
    public :: gradientdescent

    interface

    subroutine fastmarching3d_(nx,ny,nz,dx,dy,dz,xmin,ymin,&
            zmin, vel, src, rev, nrev, order, ttime, timeField) bind(C,name="fastMarching3D")
        import
        integer(c_int), intent(in), value :: nx, ny, nz
        real(c_double), intent(in), value :: dx, dy, dz
        real(c_double), intent(in), value :: xmin, ymin, zmin
        type(c_ptr), value :: vel, src, rev
        integer(c_int), intent(in), value :: nrev
        integer(c_int), intent(in), value :: order
        type(c_ptr), value :: ttime
        type(c_ptr), value :: timeField
    end subroutine

    end interface

    ! ray type
    type T_RAY
        integer srcid, revid
        integer npoints
        real(kind=c_double), dimension(:,:), allocatable :: points
        real(kind=c_double), dimension(:), allocatable :: point_velocity
    endtype

    ! shared variables (also shared across openmp)
    !integer(c_int), intent(in) :: nx, ny, nz
    !real(c_double), intent(in) :: dx, dy, dz
    !real(c_double), intent(in) :: xmin, ymin, zmin
    !logical initialised

contains

    ! The argument list is generally speaking too long, this is for
    ! compatibility, it seems there is a bug related to c_loc in gfortran
    subroutine fastmarching3d(vel,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,src,nrev,rev,order,ttime,timeField)
        implicit none
        integer(c_int), intent(in) :: nx, ny, nz
        real(c_double), intent(in) :: dx, dy, dz
        real(c_double), dimension(nz,ny,nx), intent(in), target :: vel
        real(c_double), intent(in) :: xmin, ymin, zmin
        real(c_double), dimension(3), intent(in), target :: src
        integer(c_int), intent(in) :: nrev
        real(c_double), dimension(3,nrev), intent(in), target :: rev
        integer(c_int), intent(in) :: order
        real(c_double), dimension(nrev), target, intent(inout) :: ttime
        real(c_double), dimension(nz,ny,nx), target, intent(inout), optional :: timeField

        if(present(timeField))then
            call fastmarching3d_(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,&
                c_loc(vel),c_loc(src),c_loc(rev),nrev,order,c_loc(ttime),&
                c_loc(timeField) )
        else
            call fastmarching3d_(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,&
                c_loc(vel),c_loc(src),c_loc(rev),nrev,order,c_loc(ttime),&
                c_null_ptr)
        endif

    end subroutine

    subroutine gradientdescent(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,vel,src,rev,timefield,npoints, points,point_velocity)
        implicit none
        integer(c_int), intent(in) :: nx, ny, nz
        real(c_double), intent(in) :: dx, dy, dz
        real(c_double), intent(in) :: xmin, ymin, zmin
        real(c_double), dimension(nz,ny,nx), intent(in) :: vel
        real(kind=c_double), dimension(3), intent(in) :: src
        real(kind=c_double), dimension(3), intent(in) :: rev
        real(c_double), dimension(nz,ny,nx), intent(in) :: timefield
        integer, intent(out) :: npoints
        real(c_double), dimension(:,:), intent(out) :: points
        real(c_double), dimension(:), optional, intent(out) :: point_velocity

        !real(c_double), dimension(:,:,:), allocatable :: timefield
        integer i, ix, iy, iz
        integer max_npoints
        real(c_double), dimension(3) :: grad
        real(c_double) step

        max_npoints = nx*ny*nz ! probably not larger than this
        points = 0.0
        npoints = 0
        if(present(point_velocity)) point_velocity =  0.0

        !allocate(timefield(nz+2,ny+2,nx+2))
        !timefield(2:nz+1,2:ny+1,2:nx+1) = timefield_in
        !timefield(1,:,:) = timefield(2,:,:)
        !timefield(nz+2,:,:) = timefield(nz+1,:,:)
        !timefield(:,1,:) = timefield(:,2,:)
        !timefield(:,ny+2,:) = timefield(:,ny+1,:)
        !timefield(:,:,1) = timefield(:,:,2)
        !timefield(:,:,nx+2) = timefield(:,:,nx+1)

        step = minval([dx,dy,dz])*1.0

        ! first point is receiver
        points(:,1) = rev
        call point2idx(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,points(:,1),ix,iy,iz)
        if(present(point_velocity)) point_velocity(1) = vel(iz,iy,ix)

        ! if rev point reach the src point, end it
        if(reached_end(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,rev,src,step))then
            !allocate(points(3,2))
            points(:,1) = rev
            points(:,2) = src
            !allocate(point_velocity(2))
            if(present(point_velocity)) point_velocity(2) = point_velocity(1)
            npoints = 2
            return
        endif

        i = 1
        do while(i<max_npoints-1)
            call gradient(dx,dy,dz,ix,iy,iz,timefield,grad)
            !write(*,*) 'Grad: ', grad
            !call gradient2(points(:,i),dx,dy,dz,xmin,ymin,zmin,timefield,grad)
            !write(*,*) 'Grad2: ', grad

            ! update next point
            points(:,i+1) = points(:,i) - step*grad/norm2(grad)
            i = i + 1

            ! check boundary
            if(points(1,i)<xmin)then
                points(1,i)=xmin
            elseif(points(1,i)>xmin+(nx-1)*dx)then
                points(1,i)=xmin+(nx-1)*dx
            endif
            if(points(2,i)<ymin)then
                points(2,i)=ymin
            elseif(points(2,i)>ymin+(ny-1)*dy)then
                points(2,i)=ymin+(ny-1)*dy
            endif
            if(points(3,i)<zmin)then
                points(3,i)=zmin
            elseif(points(3,i)>zmin+(nz-1)*dz)then
                points(3,i)=zmin+(nz-1)*dz
            endif
           
            call point2idx(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,points(:,i),ix,iy,iz)
            if(present(point_velocity))then
                point_velocity(i) = vel(iz,iy,ix)
            endif


            ! test if we have reached the end point
            if(reached_end(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,points(:,i),src,step)) exit

        enddo

        ! add src point finally
        npoints = i + 1
        points(:,npoints) = src

        if(present(point_velocity)) point_velocity(npoints) = point_velocity(npoints-1)

        ! using straight-ray to approximate it
        if(i==max_npoints-1)then 
            npoints = 1
            write(*,*) 'src: ', src
            write(*,*) 'receiver: ', rev
            write(*,*) 'ray cannot be found'
            !call write_doubles(timefield)
            !stop
            i = 1
            grad = - (src - rev)
            do while(i<=floor(norm2(grad)/step))

                ! update next point
                i = i + 1
                points(:,i) = points(:,i-1) - step*grad/norm2(grad)

                if(present(point_velocity))then
                    call point2idx(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,points(:,i),ix,iy,iz)
                    point_velocity(i) = vel(iz,iy,ix)
                endif

                ! test if we have reached the end point
                !if(reached_end(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,points(:,i),src,step)) exit

            enddo
            npoints = i+1
            points(:,npoints) = src
        endif

        ! debug, write out ray paths
        !open(unit=12,file='raypath.txt',position='append')
        !write(12,*) npoints
        !do i = 1, npoints
        !    write(12,*) points(:,i)
        !enddo
        !close(12)

        !allocate(rpoints(3,npoints))
        !rpoints = points(:,1:npoints)
        !allocate(rpoint_velocity(npoints))
        !rpoint_velocity = point_velocity(1:npoints)

    endsubroutine gradientdescent

    pure subroutine point2idx(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,point,ix,iy,iz)
        implicit none
        integer(c_int), intent(in) :: nx, ny, nz
        real(c_double), intent(in) :: dx, dy, dz
        real(c_double), intent(in) :: xmin, ymin, zmin
        real(kind=c_double), dimension(3), intent(in) :: point
        integer(c_int), intent(out) :: ix, iy, iz

        ix = floor((point(1)-xmin)/dx) + 1
        iy = floor((point(2)-ymin)/dy) + 1
        iz = floor((point(3)-zmin)/dz) + 1
        if(ix < 1) ix = 1
        if(iy < 1) iy = 1
        if(iz < 1) iz = 1
        if(ix >= nx) ix = nx - 1
        if(iy >= ny) iy = ny - 1
        if(iz >= nz) iz = nz - 1

    endsubroutine point2idx

    ! boundary not checked, efficiency is more important
    pure subroutine gradient(dx,dy,dz,ix,iy,iz,timefield,grad)
        implicit none
        real(c_double), intent(in) :: dx, dy, dz
        integer(c_int), intent(in) :: ix, iy, iz
        real(c_double), dimension(:,:,:), intent(in) :: timefield
        real(c_double), dimension(3), intent(out) :: grad

        integer nz, ny, nx
        nz = size(timefield,1)
        ny = size(timefield,2)
        nx = size(timefield,3)
        ! simple first order, high order could be used
        if(ix+1<=nx .and. ix-1>=1)then
            grad(1) = (timefield(iz,iy,ix+1) - timefield(iz,iy,ix-1))/(2.0*dx)
        elseif(ix+1>nx)then
            grad(1) = (timefield(iz,iy,ix) - timefield(iz,iy,ix-1))/dx
        elseif(ix-1<1)then
            grad(1) = (timefield(iz,iy,ix+1) - timefield(iz,iy,ix))/dx
        endif
        if(iy+1<=ny .and. iy-1>=1)then
            grad(2) = (timefield(iz,iy+1,ix) - timefield(iz,iy-1,ix))/(2.0*dy)
        elseif(iy+1>ny)then
            grad(2) = (timefield(iz,iy,ix) - timefield(iz,iy-1,ix))/dy
        elseif(iy-1<1)then
            grad(2) = (timefield(iz,iy+1,ix) - timefield(iz,iy,ix))/dy
        endif
        if(iz+1<=nz .and. iz-1>=1)then
            grad(3) = (timefield(iz+1,iy,ix) - timefield(iz-1,iy,ix))/(2.0*dz)
        elseif(iz+1>nz)then
            grad(3) = (timefield(iz,iy,ix) - timefield(iz-1,iy,ix))/(dz)
        elseif(iz-1<1)then
            grad(3) = (timefield(iz+1,iy,ix) - timefield(iz,iy,ix))/(dz)
        endif

    endsubroutine gradient

    ! using interpolated timefield to get gradient (test)
    subroutine gradient2(point,dx,dy,dz,xmin,ymin,zmin,timefield,grad)
        implicit none
        real(c_double), dimension(:), intent(in) :: point
        real(c_double), intent(in) :: dx, dy, dz
        real(c_double), intent(in) :: xmin, ymin, zmin
        real(c_double), dimension(:,:,:), intent(in) :: timefield
        real(c_double), dimension(3), intent(out) :: grad

        integer nz, ny, nx
        integer ix, iy, iz
        real(c_double) step, stepl 
        real(c_double) t1, t2
        real(c_double), dimension(3) :: dp, pt1, pt2

        nz = size(timefield,1)
        ny = size(timefield,2)
        nx = size(timefield,3)

        step = minval([dx,dy,dz])

        ! x direction
        dp = [step, 0.0_ii10, 0.0_ii10]
        pt1 = point - dp
        pt2 = point + dp
        stepl = 2*step
        if(pt1(1)<xmin)then
            pt1 = point
            stepl = step
        elseif(pt2(1)>xmin+(nx-1)*dx)then
            pt2 = point
            stepl = step
        endif
        t1 = Interpolate3d(timefield,dx,dy,dz,xmin,ymin,zmin,point-dp)
        t2 = Interpolate3d(timefield,dx,dy,dz,xmin,ymin,zmin,point+dp)
        grad(1) = (t2-t1)/stepl
        write(*,*) 'point: ', point, ' t2: ', t2, ' t1: ', t1, ' grad: ', grad(1)

        ! y direction
        dp = [0.0_ii10, step, 0.0_ii10]
        pt1 = point - dp
        pt2 = point + dp
        stepl = 2*step
        if(pt1(2)<ymin)then
            pt1 = point
            stepl = step
        elseif(pt2(2)>ymin+(ny-1)*dy)then
            pt2 = point
            stepl = step
        endif
        t1 = Interpolate3d(timefield,dx,dy,dz,xmin,ymin,zmin,point-dp)
        t2 = Interpolate3d(timefield,dx,dy,dz,xmin,ymin,zmin,point+dp)
        grad(2) = (t2-t1)/stepl

        ! z direction
        dp = [0.0_ii10, 0.0_ii10, step]
        pt1 = point - dp
        pt2 = point + dp
        stepl = 2*step
        if(pt1(3)<zmin)then
            pt1 = point
            stepl = step
        elseif(pt2(3)>zmin+(nz-1)*dz)then
            pt2 = point
            stepl = step
        endif
        t1 = Interpolate3d(timefield,dx,dy,dz,xmin,ymin,zmin,point-dp)
        t2 = Interpolate3d(timefield,dx,dy,dz,xmin,ymin,zmin,point+dp)
        grad(3) = (t2-t1)/stepl

    endsubroutine gradient2

    pure function reached_end(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,point,end_point,step)
        implicit none
        integer(c_int), intent(in) :: nx, ny, nz
        real(c_double), intent(in) :: dx, dy, dz
        real(c_double), intent(in) :: xmin, ymin, zmin
        real(c_double), dimension(:), intent(in) :: point
        real(c_double), dimension(:), intent(in) :: end_point
        real(c_double), intent(in) :: step
        logical reached_end

        integer ix, iy, iz, ix_end, iy_end, iz_end

        reached_end = .false.

        if(norm2(point-end_point)<step)then
            reached_end = .true.
            return
        endif

        ! if we are in end point cell, end it
        call point2idx(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,point,ix,iy,iz)
        call point2idx(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,end_point,ix_end,iy_end,iz_end)
        if(ix==ix_end .and. iy==iy_end .and. iz==iz_end)then
            reached_end = .true.
        endif

    endfunction reached_end

    pure function Interpolate3d(array,dx,dy,dz,xmin,ymin,zmin,point) result(qv)
        implicit none
        real(kind=ii10), dimension(:,:,:), intent(in) :: array
        real(c_double), intent(in) :: dx, dy, dz
        real(c_double), intent(in) :: xmin, ymin, zmin
        real(kind=ii10), dimension(:), intent(in) :: point
        real(kind=ii10) :: qv

        integer nx,ny,nz
        integer ix, iy, iz, i, j, k
        real(kind=ii10) dsx, dsy, dsz
        real(kind=ii10) weight

        nz = size(array,1)
        ny = size(array,2)
        nx = size(array,3)

        call point2idx(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,point,ix,iy,iz)
        dsx = point(1) - (xmin + (ix-1)*dx)
        dsy = point(2) - (ymin + (iy-1)*dy)
        dsz = point(3) - (zmin + (iz-1)*dz)

        qv = 0
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    weight = (1.0-abs((i-1)*dx-dsx)/dx)*&
                             (1.0-abs((j-1)*dy-dsy)/dy)*&
                             (1.0-abs((k-1)*dz-dsz)/dz)
                    qv = qv + weight*array(iz+k-1,iy+j-1,ix+i-1)
                enddo
            enddo
        enddo

        ! debug
        !print *, array(iz:iz+1,iy:iy+1,ix:ix+1)
        !print *, qv
    end function


end module fastmarching
