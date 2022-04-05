! settings for McMC running
!
module m_settings

    use iso_c_binding, only : c_size_t, c_int, c_double
    use m_utils, only : ii10, STRLEN, itoa, eps
    use m_exception,only : exception_raiseError
    use m_logger, only          : log_msg
    implicit none

    private

    public :: T_GRID, grid_setup, point2idx
    public :: T_MOD, mod_setup
    public :: T_MCMC_SET
    public :: out_bnd
    public :: settings_check
    public :: prior_check

    ! grid settings
    type T_GRID
        integer(c_int)    :: nx, ny, nz
        real( kind=ii10 ) :: xmin, xmax
        real( kind=ii10 ) :: ymin, ymax
        real( kind=ii10 ) :: zmin, zmax
        real( kind=ii10 ) :: dx, dy, dz
        real( kind=ii10 ) :: waterDepth
        character(len=STRLEN) :: waterDepthFile
        real( kind=ii10 ) :: scaling
    endtype

    ! model
    type T_MOD
        real( kind=ii10 ), dimension(:,:,:), allocatable :: vp, vs, rho
        real( kind=ii10 ), dimension(:,:), allocatable :: waterdepth
        real(kind=ii10), dimension(:,:,:), allocatable :: vpmin_array, vpmax_array 
        real(kind=ii10), dimension(:,:,:), allocatable :: vsmin_array, vsmax_array 
        real(kind=ii10), dimension(:,:,:), allocatable :: rhomin_array, rhomax_array 
    end type T_MOD

    ! mcmc settings
    type T_MCMC_SET
        ! initial model
        integer         :: initialise
        real(kind=ii10) :: init_vsmin, init_vsmax
        character(len=STRLEN) :: initial_model
        ! number of samples
        integer         :: processor
        integer(8)      :: burn_in
        integer         :: thin
        integer         :: nsamples
        integer         :: resume
        integer         :: display
        integer         :: runtime_step
        integer         :: sigdep 
        integer         :: ncell_max
        integer         :: ncell_min
        integer :: kernel
        real(kind=ii10) :: vpmin, vpmax
        real(kind=ii10) :: vsmin, vsmax
        real(kind=ii10) :: rhomin, rhomax
        character(len=STRLEN) :: prior_model
        !real(kind=ii10), dimension(3) :: pm1, pm2
        real(kind=ii10) :: bn0_min, bn0_max
        real(kind=ii10) :: bn1_min, bn1_max
        real(kind=ii10) :: sn0_min, sn0_max
        real(kind=ii10) :: sn1_min, sn1_max
        real(kind=ii10) :: sigma_vp, sigma_vp2
        real(kind=ii10) :: sigma_vs, sigma_vs2
        real(kind=ii10) :: sigma_rho, sigma_rho2
        real(kind=ii10) :: sigma_bn0, sigma_bn1
        real(kind=ii10) :: sigma_sn0, sigma_sn1
        real(kind=ii10) :: pd, pd2

        ! datatype
        integer         :: datatype

        ! located
        integer         :: locate, nloc
        real(kind=ii10) :: sigma_x, sigma_y, sigma_z, sigma_t
        real(kind=ii10) :: xwidth, ywidth, zwidth, twidth
        ! grid
        type(T_GRID)    :: grid

        ! tempering
        integer         :: tempering
        character(STRLEN):: temperfile
        integer         :: tempering_start
        integer         :: tempering_step
        integer         :: jump_type
        integer         :: number_of_temperatures, number_of_1s

        ! hmc
        integer         :: priortype
        integer         :: hmc
        integer         :: n_leapfrogs
        real(kind=ii10) :: vp_step, vs_step, src_step, src_tstep
        real(kind=ii10) :: vp_mass, vs_mass, src_mass

        ! slice sampling
        integer         :: slicesample

    endtype


contains
    subroutine grid_setup(grid)
        implicit none
        type(T_GRID), intent(inout) :: grid

        ! if appropriate, scaling first
        grid%zmin = grid%zmin*grid%scaling
        grid%zmax = grid%zmax*grid%scaling

        if( grid%nx/=0 .and. grid%ny/=0 .and. grid%nz/=0)then
            grid%dx = (grid%xmax - grid%xmin)/(grid%nx-1)
            grid%dy = (grid%ymax - grid%ymin)/(grid%ny-1)
            grid%dz = (grid%zmax - grid%zmin)/(grid%nz-1)
        else
            grid%xmin = 0
            grid%xmax = 0
            grid%ymin = 0
            grid%ymax = 0
            grid%zmin = 0
            grid%zmax = 0
            grid%dx = 0
            grid%dy = 0
            grid%dz = 0
        endif
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

    subroutine mod_setup(model,mcmc_set)
        implicit none
        type(T_MOD), intent(out) :: model
        type(T_MCMC_SET), intent(in) :: mcmc_set

        type(T_GRID) grid
        integer ierr
        integer nz, ny, nx, i, j
	    logical lexist
        real(kind=ii10), dimension(:,:), allocatable :: array

        grid = mcmc_set%grid
        nz = mcmc_set%grid%nz
        ny = mcmc_set%grid%ny
        nx = mcmc_set%grid%nx
        allocate( model%vp(grid%nz, grid%ny, grid%nx) )
        allocate( model%vs(grid%nz, grid%ny, grid%nx) )
        allocate( model%rho(grid%nz, grid%ny, grid%nx) )

        model%vp = 1.0 ! safe
        model%vs = 1.0
        model%rho = 1.0

        allocate(model%waterdepth(grid%ny, grid%nx))
        model%waterdepth = 0
        if(grid%waterDepth>eps)then
            call read_waterdepth(grid%waterdepthfile,model%waterdepth,ierr)
        endif
        if(ierr /= 0)then
            model%waterdepth = grid%waterDepth
        endif

        allocate(model%vpmin_array(nz,ny,nx))
        allocate(model%vpmax_array(nz,ny,nx))
        allocate(model%vsmin_array(nz,ny,nx))
        allocate(model%vsmax_array(nz,ny,nx))

        inquire(file=mcmc_set%prior_model,exist=lexist)
        if(.not.lexist)then
            call log_msg('No prior file found, using scalar prior')
            model%vpmin_array = mcmc_set%vpmin
            model%vpmax_array = mcmc_set%vpmax
            model%vsmin_array = mcmc_set%vsmin
            model%vsmax_array = mcmc_set%vsmax
        else
            call readtxt(array,mcmc_set%prior_model)
            do i = 1, nx
                do j = 1, ny
                    model%vpmin_array(:,j,i) = array(0,:)
                    model%vpmax_array(:,j,i) = array(1,:)
                    model%vsmin_array(:,j,i) = array(2,:)
                    model%vsmax_array(:,j,i) = array(3,:)
                enddo
            enddo
        endif

    end subroutine mod_setup

    subroutine read_waterdepth(filename,waterdepth,ierr)
        use m_utils, only : read_resume_unit
        implicit none
        character( len=* ), intent(in) :: filename
        real(kind=ii10), dimension(:,:), intent(out) :: waterdepth
        integer, intent(out) :: ierr

        integer i
        
        ierr = 0
        open(unit=read_resume_unit,file=filename,status='old',action='read')
        do i = 1, size(waterdepth,2)
            read(read_resume_unit,*,iostat=ierr) waterdepth(:,i) 
            if(ierr /= 0)then
                call log_msg('Reading water depth error! Set waterdepth to a scalar value provided in waterDepth!')
                waterdepth = 0
                exit
            endif
        enddo
        close(read_resume_unit)

    end subroutine


    logical function out_bnd(coord, bnd)
        implicit none
        real( kind=ii10 ), dimension(:), intent(in) :: coord
        real( kind=ii10 ), dimension(:), intent(in) :: bnd
    
        if( any(coord < bnd(1:size(bnd):2)) .or. any(coord > bnd(2:size(bnd):2)) )then
            out_bnd = .true.
        else
            out_bnd = .false.
        endif

    end function

    ! check prior 
    logical function prior_check(model)
        implicit none
        type(T_MOD), intent(in) :: model

        prior_check = .false.
        if(any(model%vp < model%vpmin_array) .or. any(model%vp > model%vpmax_array) &
            .or. any(model%vs < model%vsmin_array) .or. any(model%vs > model%vsmax_array) )then
            prior_check = .true.
        endif
        
    end function

    ! default values for mcmc settings and check validity of settings
    subroutine settings_check(mcmc_set)
        implicit none
        type(T_MCMC_SET), intent(inout) :: mcmc_set

        ! check the validity of settings
	    logical lexist

        if(mcmc_set%initialise==1)then
            inquire(file=mcmc_set%initial_model,exist=lexist)
            if(.not.lexist)&
                call exception_raiseError('Initial model file does not exist!')
        endif

        if(mcmc_set%thin == 0)then
            call exception_raiseError('Thinning of the chain cannot be zero.')
        endif

        if(mcmc_set%nsamples == 0)then
            call exception_raiseError('The number of samples cannot be zero.')
        endif

        if(mcmc_set%display == 0 .or. mcmc_set%runtime_step==0)then
            call exception_raiseError('Display and run time step cannot be zero.')
        endif

        if(mcmc_set%ncell_min<4 .or. mcmc_set%ncell_max<mcmc_set%ncell_min)then
            call exception_raiseError('The minimun number of cells cannot be &
            &smaller than 4; the max number of cells cannot be smaller than &
            &the minimum number!')
        endif

        if(mcmc_set%datatype > 3)then
            call exception_raiseError('Invalid data type.')
        endif

        if(mcmc_set%datatype==2 .and. mcmc_set%locate == 1)then
            call exception_raiseError('The code does not support to change &
            &source locations for surface wave inversion!')
        endif

        if(mcmc_set%tempering==1)then
            inquire(file=mcmc_set%temperfile,exist=lexist)
            if(.not.lexist)&
                call exception_raiseError('Temperature file does not exist!')
            if(mcmc_set%tempering_step==0)&
                call exception_raiseError('Tmepering step cannot be zero, can be&
                &one.')
            if(mcmc_set%number_of_temperatures<=mcmc_set%number_of_1s)&
                call exception_raiseError('The number of temperatures must be &
                &larger than the number of chains with temperature of one.')
        endif

        if(mcmc_set%locate==1 .and. mcmc_set%nloc==0)then
            mcmc_set%nloc=1
        endif
        call log_msg('Number of location change per iter:'//itoa(mcmc_set%nloc))

        !mcmc_set%pm1(1) = mcmc_set%vpmin
        !mcmc_set%pm1(2) = mcmc_set%vsmin
        !mcmc_set%pm1(3) = mcmc_set%rhomin

        !mcmc_set%pm2(1) = mcmc_set%vpmax
        !mcmc_set%pm2(2) = mcmc_set%vsmax
        !mcmc_set%pm2(3) = mcmc_set%rhomax
        
        ! disable slice sampling and HMC sampling for now
        mcmc_set%slicesample = 0
        mcmc_set%hmc = 0
        mcmc_set%n_leapfrogs = 50
        mcmc_set%vp_step = 0.01
        mcmc_set%vs_step = 0.01
        mcmc_set%src_step = 0.01
        mcmc_set%vp_mass = 1
        mcmc_set%vs_mass = 1
        mcmc_set%src_mass = 1

    endsubroutine

    subroutine readtxt(array, filename)
        use m_utils, only : read_resume_unit, read_doubles
        implicit none
        character(len=*), intent(in) :: filename
        real(kind=ii10), dimension(:,:), allocatable, intent(out) :: array
        integer n1, n2, i

        open(unit=read_resume_unit,file=filename,status='old',action='read')
        read(read_resume_unit,*) n1, n2
        allocate(array(n2,n1))
        do i = 1, n1
            read(read_resume_unit,*) array(:,i)
        enddo
        close(read_resume_unit)
    endsubroutine


end module
