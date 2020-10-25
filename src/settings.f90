! settings for McMC running
!
module m_settings

    use iso_c_binding, only : c_size_t, c_int
    use m_utils, only : ii10, STRLEN, itoa
    use m_exception,only : exception_raiseError
    use m_logger, only          : log_msg
    implicit none

    private

    public :: T_GRID, grid_setup
    public :: T_MOD, mod_setup
    public :: T_MCMC_SET
    public :: out_bnd
    public :: settings_check

    ! grid settings
    type T_GRID
        integer(c_int)    :: nx, ny, nz
        real( kind=ii10 ) :: xmin, xmax
        real( kind=ii10 ) :: ymin, ymax
        real( kind=ii10 ) :: zmin, zmax
        real( kind=ii10 ) :: dx, dy, dz
        real( kind=ii10 ) :: waterDepth
        real( kind=ii10 ) :: scaling
    endtype

    ! model
    type T_MOD
        real( kind=ii10 ), dimension(:,:,:), allocatable :: vp, vs, rho
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

    subroutine mod_setup(model,grid)
        implicit none
        type(T_MOD), intent(out) :: model
        type(T_GRID), intent(in) :: grid

        !type(T_GRID) grid
        !grid = mcmc_set%grid

        allocate( model%vp(grid%nz, grid%ny, grid%nx) )
        allocate( model%vs(grid%nz, grid%ny, grid%nx) )
        allocate( model%rho(grid%nz, grid%ny, grid%nx) )

        model%vp = 1.0 ! safe
        model%vs = 1.0
        model%rho = 1.0

    end subroutine mod_setup

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

end module
