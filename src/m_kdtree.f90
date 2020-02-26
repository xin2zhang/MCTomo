! module to convert the voronoi diagram to grid based model
! using kd-tree, which is faster than walker in the delaunay
! triangulation
!
module m_kdtree

    use kdtree2_precision_module, only : kdkind
    use kdtree2_module, only : kdtree2, kdtree2_result, kdtree2_create,&
    kdtree2_n_nearest
    
    use iso_c_binding

    implicit none

contains

    subroutine kdtree_to_grid(points,parameters,np, x0,y0,z0, x1,y1,z1, bnd_box,&
            nx,ny,nz, vp,vs,rho)
        use cgal_delaunay, only : d3
        use kdtree2_precision_module, only : kdkind
        use kdtree2_module, only : kdtree2, kdtree2_result, kdtree2_create,&
            kdtree2_n_nearest
        implicit none
        real(kdkind), dimension(:,:), intent(in) :: points, parameters
        integer, intent(in) :: np
	real(kdkind), intent(in) :: x0,y0,z0
	real(kdkind), intent(in) :: x1,y1,z1
	integer(c_int), intent(in) :: nx,ny,nz
	type(d3), dimension(2), intent(in) :: bnd_box
	real(kind=c_double), dimension(:,:,:), target, intent(inout) :: vp, vs, rho

        ! kd-tree
        type(kdtree2), pointer :: tree
        type(kdtree2_results), dimension(1) :: results
        real(kdkind), dimension(3) :: qv

        ! grid
        real(c_double) :: dx, dy, dz
        integer        :: ix0, ix1, iy0, iy1, iz0, iz1

        integer i,j,k

        ! create kd-tree
        tree => kdtree2_create(points(:,1:np), sort=.false., rearrange=.true.)

        ! convert to grid based model using kd-tree nearest search
        dx = (x1-x0)/nx
        dy = (y1-y0)/ny
        dz = (z1-z0)/nz

        ! the lowest and highest indices for grid model vp, vs and rho
        ix0 = floor( (bnd(1)%x-x0)/dx ) + 1
        ix1 = floor( (bnd(2)%x-x0)/dx ) + 1
        iy0 = floor( (bnd(1)%y-y0)/dy ) + 1
        iy1 = floor( (bnd(2)%y-y0)/dy ) + 1
        iz0 = floor( (bnd(1)%z-z0)/dz ) + 1
        iz1 = floor( (bnd(2)%z-z0)/dz ) + 1
        if(ix1>nx) ix1 = nx
        if(iy1>ny) iy1 = ny
        if(iz1>nz) iz1 = nz

        ! call nearest search of kd-tree
        do i = ix0, ix1
            do j = iy0, iy1
                do k = iz0, iz1
                    qv = [x0+(i-1)*dx, y0+(j-1)*dy, z0+(k-1)*dz]
                    call kdtree2_n_nearest(tp=tree,qv=qv,nn=1,results=results)
                    vp(k,j,i) = parameters(1,results(1)%idx)
                    vs(k,j,i) = parameters(2,results(1)%idx)
                    rho(k,j,i) = parameters(3,results(1)%idx)
                enddo
            enddo
        enddo

        return

    end subroutine
end module

