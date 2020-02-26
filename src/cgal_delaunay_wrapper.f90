!GNU GPL license v3
module cgal_delaunay
    use iso_c_binding
    use m_exception, only : exception_raiseError

    implicit none

    private

    ! the fortran interface to cgal cpp library
    public :: cgal_delaunay_build
    public :: cgal_delaunay_copy
    public :: delaunay_size
    public :: cgal_get_vertices
    public :: cgal_reset_values
    public :: delaunay_vertex
    public :: vertex_point
    public :: vertex_value
    public :: cgal_delaunay_insert
    public :: cgal_delaunay_remove
    public :: cgal_delaunay_move
    public :: cgal_delaunay_value
    public :: cgal_delaunay_box
    public :: cgal_delaunay_locate
    public :: cgal_delaunay_to_grid
    public :: cgal_delaunay_read
    public :: cgal_delaunay_write
    public :: cgal_vertex_write
    public :: cgal_delaunay_delete

    ! derived type for points and associated parameters and its transformation
    ! with arrays
    public :: p3, d3
    public :: points_array
    public :: array_points
    public :: array_parameters

    interface points_array
        module procedure points_array_1, points_array_2
    end interface

    interface array_points
        module procedure array_points_1, array_points_2
    end interface

    interface array_parameters
        module procedure array_parameters_1, array_parameters_2
    end interface

    ! derived type for parameter binding with c struct
    type, bind(C) :: p3
        real(c_double) :: vp, vs, rho
    end type

    type, bind(C) :: d3
	    real(c_double) :: x, y, z
    end type

    interface
	
	function delaunay_build(ppoint, pparameter, n) result(pt)&
        bind(C, name="delaunay_build")
	    import
	    type(c_ptr) :: pt
	    type(c_ptr), value :: ppoint
	    type(c_ptr), value :: pparameter
	    integer(c_size_t), intent(in), value :: n
	end function delaunay_build

    function delaunay_size(pt) result(n) bind(C, name="delaunay_size")
        import
        type(c_ptr), value :: pt
        integer(c_size_t) :: n
    end function

    function delaunay_copy(pt) result(pt_copy) bind(C, name="delaunay_copy")
	    import
        type(c_ptr), value, intent(in) :: pt
        type(c_ptr) :: pt_copy
    end function

    subroutine get_vertices(pt,points,parameters) bind(C,name="get_vertices")
        import
        type(c_ptr), value, intent(in) :: pt
        type(c_ptr), value :: points
        type(c_ptr), value :: parameters
    end subroutine

    subroutine reset_values(pt,parameters) bind(C,name="reset_values")
        import
        type(c_ptr), value, intent(in) :: pt
        type(c_ptr), value :: parameters
    end subroutine

    function vertex_point(vt) result(point) bind(C, name="vertex_point")
        import
        type(c_ptr), value, intent(in) :: vt
        type(d3) :: point
    end function

    function vertex_value(vt) result(pm) bind(C, name="vertex_value")
        import
        type(c_ptr), value, intent(in) :: vt
        type(p3) :: pm
    end function

    function delaunay_vertex(pt, ind) result(vt) bind(C, name="delaunay_vertex")
        import
        type(c_ptr), value, intent(in) :: pt
        integer(c_int), intent(in), value :: ind
        type(c_ptr) :: vt
    end function

	subroutine delaunay_insert(pt, pinsert, pparameter, point0, point1, bnd_box)&
        bind(C, name="delaunay_insert")
	    import
	    type(c_ptr), value :: pt
	    type(d3), intent(in) :: pinsert
	    type(p3), intent(in) :: pparameter
	    type(d3), intent(in) :: point0
	    type(d3), intent(in) :: point1
	    type(c_ptr), value :: bnd_box
	end subroutine

	subroutine delaunay_remove(pt, removal, point0, point1, removed_pm, bnd_box)&
        bind(C, name="delaunay_remove")
	    import
	    type(c_ptr), value :: pt
            type(d3), intent(in) :: removal
	    type(d3), intent(in) :: point0
	    type(d3), intent(in) :: point1
            type(p3), intent(out) :: removed_pm
	    type(c_ptr), value :: bnd_box
	end subroutine delaunay_remove

	subroutine delaunay_move(pt, p_src, p_dst, point0, point1, bnd_box, verbose)&
        bind(C, name="delaunay_move")
	    import
	    type(c_ptr), value :: pt
	    type(d3), intent(in) :: p_src
	    type(d3), intent(in) :: p_dst
	    type(d3), intent(in) :: point0
	    type(d3), intent(in) :: point1
	    type(c_ptr), value :: bnd_box
            integer(c_int), intent(out) :: verbose
	end subroutine delaunay_move

    subroutine delaunay_value(pt, p_src, pm_dst, point0, point1, &
        bnd_box, verbose) bind(C, name="delaunay_value")
        import
        type(c_ptr), value :: pt
        type(d3), intent(in) :: p_src
        type(p3), intent(in) :: pm_dst
	    type(d3), intent(in) :: point0
	    type(d3), intent(in) :: point1
	    type(c_ptr), value :: bnd_box
        integer(c_int) :: verbose
    end subroutine

    subroutine delaunay_box(pt, p_src, point0, point1, &
        bnd_box) bind(C, name="delaunay_box")
        import
        type(c_ptr), value :: pt
        type(d3), intent(in) :: p_src
	    type(d3), intent(in) :: point0
	    type(d3), intent(in) :: point1
	    type(c_ptr), value :: bnd_box
    end subroutine

	subroutine cgal_delaunay_locate(pt, ppoint, pparameter) bind(C, name="delaunay_locate")
	    import
	    type(c_ptr), value, intent(in) :: pt
	    type(d3), intent(in) :: ppoint
	    type(p3) :: pparameter
	end subroutine cgal_delaunay_locate

    subroutine delaunay_to_grid(pt, ppoint0, bnd0, bnd1, dx, dy, dz,&
                nx, ny, nz, vp, vs,rho) bind(C, name="delaunay_fast_to_grid")
        import
        type(c_ptr), value, intent(in) :: pt
        type(d3), intent(in) :: ppoint0
        type(d3), intent(in) :: bnd0
        type(d3), intent(in) :: bnd1
        real(c_double), value, intent(in) :: dx, dy, dz
        integer(c_int), intent(in), value :: nx, ny, nz
        type(c_ptr), value :: vp, vs, rho
    end subroutine 

    function delaunay_read(filename) result(pt) bind(C, name="delaunay_read")
        import
        type(c_ptr) :: pt
        character(len=1,kind=c_char), dimension(*), intent(in) :: filename
    end function
	
    subroutine delaunay_write(pt,filename) bind(C, name="delaunay_write")
        import
        type(c_ptr), value :: pt
        character(len=1,kind=c_char), dimension(*), intent(in) :: filename
    end subroutine
	
    subroutine vertex_write(pt,filename) bind(C, name="vertex_write")
        import
        type(c_ptr), value :: pt
        character(len=1,kind=c_char), dimension(*), intent(in) :: filename
    end subroutine
	
	subroutine delaunay_delete(pt) bind(C, name="delaunay_delete")
	    import
	    type(c_ptr), value :: pt
	end subroutine delaunay_delete
    end interface

contains

    subroutine cgal_delaunay_build(points, parameters, nsites, triangulation)
	    integer(c_size_t) nsites
	    type(d3), dimension(nsites), target, intent(in) :: points
	    type(p3), dimension(nsites), target, intent(in) :: parameters
	    type(c_ptr), intent(out) :: triangulation

        integer(c_size_t) n
	    n = nsites
	    if( n < 3) call exception_raiseError('Error: the number of points smaller than 3')

	    triangulation = delaunay_build(c_loc(points(1)), c_loc(parameters(1)), n)

	    return
    end subroutine

    subroutine cgal_delaunay_copy(pt1,pt2)
        type(c_ptr),intent(in) :: pt1
        type(c_ptr),intent(inout) :: pt2

        if(c_associated(pt2) ) call cgal_delaunay_delete(pt2)
        pt2 = delaunay_copy(pt1)
    end subroutine

    subroutine cgal_get_vertices(triangulation, nsites, points, parameters)
        type(c_ptr), intent(in) :: triangulation
        ! the explictly shaped array is needed by gfortran < 4.9
        integer, intent(in) :: nsites
        real(c_double), dimension(3,nsites), target, intent(inout) :: points
        real(c_double), dimension(3,nsites), target, intent(inout) :: parameters

        call get_vertices(triangulation,c_loc(points),c_loc(parameters))

        return

    end subroutine
	
    subroutine cgal_reset_values(triangulation, n, parameters)
        type(c_ptr), intent(in) :: triangulation
        ! the explictly shaped array is needed by gfortran < 4.9
        integer, intent(in) :: n
        real(c_double), dimension(3,n), target, intent(inout) :: parameters

        call reset_values(triangulation,c_loc(parameters))

        return

    end subroutine
	
    subroutine cgal_delaunay_insert(triangulation, insert_pt, insert_pm, point0, point1, bnd_box)
	    type(c_ptr), intent(inout) :: triangulation
	    type(d3), intent(in) :: insert_pt
	    type(p3), intent(in) :: insert_pm
        type(d3), intent(in) :: point0
	    type(d3), intent(in) :: point1
	    type(d3), dimension(2), target, intent(inout) :: bnd_box

	    if( .not.c_associated(triangulation) ) &
	        call exception_raiseError('Error: triangulation is empty!')

	    call delaunay_insert(triangulation, insert_pt, insert_pm, point0, point1, c_loc(bnd_box))

	    return
    end subroutine

    subroutine cgal_delaunay_remove(triangulation,  point0, point1,&
    removed_pt, removed_pm, bnd_box)
	    type(c_ptr) :: triangulation
        type(d3), intent(in) :: point0
	    type(d3), intent(in) :: point1
        type(d3), intent(in) :: removed_pt
        type(p3), intent(out) :: removed_pm
	    type(d3), dimension(2), target, intent(inout) :: bnd_box

	    if( .not.c_associated(triangulation) ) &
	        call exception_raiseError('Error: triangulation is empty!')

	    call delaunay_remove(triangulation, removed_pt, point0, point1,&
                removed_pm, c_loc(bnd_box) )

	    return
    end subroutine

    subroutine cgal_delaunay_move(triangulation, p_src, p_dst, point0, point1,&
            bnd_box, verbose)
	    type(c_ptr), intent(inout) :: triangulation
	    type(d3), intent(in) :: p_src
	    type(d3), intent(in) :: p_dst
        type(d3), intent(in) :: point0
	    type(d3), intent(in) :: point1
	    type(d3), dimension(2), intent(inout), target :: bnd_box
        integer(c_int), intent(out) :: verbose
	
	    if( .not.c_associated(triangulation) ) &
	        call exception_raiseError('Error: triangulation is empty!')

	    call delaunay_move(triangulation, p_src, p_dst, point0, point1, &
                c_loc(bnd_box), verbose )
	    !if( .not.c_associated(vertex) ) &
	    !    call exception_raiseError('Error: move point conflicts!')

	    return
    end subroutine

    subroutine cgal_delaunay_value(triangulation, p_src, pm_dst, point0, point1,&
            bnd_box, verbose)
	    type(c_ptr), intent(inout) :: triangulation
	    type(d3), intent(in) :: p_src
	    type(p3), intent(in) :: pm_dst
        type(d3), intent(in) :: point0
	    type(d3), intent(in) :: point1
	    type(d3), dimension(2), intent(inout), target :: bnd_box
        integer(c_int), intent(out) :: verbose
	    
	    if( .not.c_associated(triangulation) ) &
	        call exception_raiseError('Error: triangulation is empty!')

	    call delaunay_value(triangulation, p_src, pm_dst, point0, point1, &
            c_loc(bnd_box), verbose )

	    return
    end subroutine

    subroutine cgal_delaunay_box(triangulation, p_src, point0, point1,bnd_box)
	    type(c_ptr), intent(inout) :: triangulation
	    type(d3), intent(in) :: p_src
        type(d3), intent(in) :: point0
	    type(d3), intent(in) :: point1
	    type(d3), dimension(2), intent(inout), target :: bnd_box
	    
	    if( .not.c_associated(triangulation) ) &
	        call exception_raiseError('Error: triangulation is empty!')

	    call delaunay_box(triangulation, p_src, point0, point1, c_loc(bnd_box))

    end subroutine

    subroutine cgal_delaunay_read(dt, filename)
        type(c_ptr), intent(out) :: dt
        character(len=*), intent(in) :: filename
        character(len=1,kind=c_char), dimension(len(filename)+1) :: c_filename

        c_filename = f_c_string_func(filename)
        dt = delaunay_read(c_filename)
        
        return
    end subroutine

    subroutine cgal_delaunay_write(dt, filename)
        type(c_ptr), intent(in) :: dt
        character(len=*), intent(in) :: filename
        character(len=1,kind=c_char), dimension(len(filename)+1) :: c_filename

        c_filename = f_c_string_func(filename)
        call delaunay_write(dt,c_filename)
        
        return
    end subroutine

    subroutine cgal_vertex_write(dt, filename)
        type(c_ptr), intent(in) :: dt
        character(len=*), intent(in) :: filename
        character(len=1,kind=c_char), dimension(len(filename)+1) :: c_filename

        c_filename = f_c_string_func(filename)
        call vertex_write(dt,c_filename)
        
        return
    end subroutine

    subroutine cgal_delaunay_delete(dt1, dt2)
	    type(c_ptr) :: dt1
        type(c_ptr),optional :: dt2

	    call delaunay_delete(dt1)
	    if( present(dt2) ) call delaunay_delete(dt2)

	    return
    end subroutine

    subroutine cgal_delaunay_to_grid(triangulation, x0,y0,z0, x1,y1,z1, bnd_box, &
	    				nx,ny,nz, vp,vs,rho)
	    type(c_ptr) :: triangulation
	    real(c_double), intent(in) :: x0,y0,z0
	    real(c_double), intent(in) :: x1,y1,z1
	    integer(c_int), intent(in) :: nx,ny,nz
	    type(d3), dimension(2), intent(in) :: bnd_box
	    real(kind=c_double), dimension(nz,ny,nx), target, intent(inout) :: vp, vs, rho

        real(c_double) :: dx, dy, dz
            
	    if( .not.c_associated(triangulation) ) &
	        call exception_raiseError('Error: triangulation is empty!')

        dx = (x1-x0)/nx
        dy = (y1-y0)/ny
        dz = (z1-z0)/nz

	    call delaunay_to_grid( triangulation, d3(x0,y0,z0), bnd_box(1), bnd_box(2),&
			      dx,dy,dz,nx,ny,nz,c_loc(vp),c_loc(vs),c_loc(rho) )

        return   

    end subroutine

    pure function f_c_string_func (f_string) result (c_string)
        use, intrinsic :: iso_c_binding, only: c_char, c_null_char
        implicit none
        character(len=*), intent(in) :: f_string
        character(len=1,kind=c_char) :: c_string(len_trim(f_string)+1)
        integer                      :: n, i
    
        n = len_trim(f_string)
        do i = 1, n
          c_string(i) = f_string(i:i)
        end do
        c_string(n + 1) = c_null_char
    
    end function f_c_string_func

    subroutine points_array_1(d3_point,point)
        implicit none
        type(d3), intent(in) :: d3_point
        real(kind=c_double), dimension(3), intent(inout) :: point

        point(1) = d3_point%x
        point(2) = d3_point%y
        point(3) = d3_point%z
        return
    end subroutine

    subroutine points_array_2(d3_points,points)
        implicit none
        type(d3), dimension(:), intent(in) :: d3_points
        real(kind=c_double), dimension(:,:), intent(inout) :: points
        integer i

        do i = 1, size(d3_points)
            call points_array_1(d3_points(i),points(:,i))
        enddo
        return
    end subroutine

    subroutine array_points_1(d3_point,point)
        implicit none
        type(d3), intent(inout) :: d3_point
        real(kind=c_double), dimension(3), intent(in) :: point

        d3_point%x = point(1)
        d3_point%y = point(2)
        d3_point%z = point(3)
        return
    end subroutine

    subroutine array_points_2(d3_points,points)
        implicit none
        type(d3), dimension(:), intent(inout) :: d3_points
        real(kind=c_double), dimension(:,:), intent(in) :: points
        integer i

        do i = 1, size(points,2)
            call array_points_1(d3_points(i),points(:,i))
        enddo
        return
    end subroutine

    subroutine array_parameters_1(d3_point,point)
        implicit none
        type(p3), intent(inout) :: d3_point
        real(kind=c_double), dimension(3), intent(in) :: point

        d3_point%vp = point(1)
        d3_point%vs = point(2)
        d3_point%rho = point(3)
        return
    end subroutine

    subroutine array_parameters_2(d3_points,points)
        implicit none
        type(p3), dimension(:), intent(inout) :: d3_points
        real(kind=c_double), dimension(:,:), intent(in) :: points
        integer i

        do i = 1, size(points,2)
            call array_parameters_1(d3_points(i),points(:,i))
        enddo
        return
    end subroutine

end module cgal_delaunay
