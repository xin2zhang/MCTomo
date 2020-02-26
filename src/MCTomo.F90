! main program
! MCTomo: implement 3d body wave and surface wave tomography using reversible jump
! mcmc algorithm and 3d Voronoi parameterization.
! GNU GPL license v3

program MCTomo

    use iso_c_binding
    use omp_lib

    use m_logger,       only : log_startup, log_shutdown, log_msg
    use m_exception,    only : exception_raiseError
    use m_utils,        only : newunit, ii10, itoa, rtoa, results_dir, &
                               FILEDIR_SEPARATOR, log_dir, create_dir
    use m_settings,     only : T_MCMC_SET, T_GRID, grid_setup
    use like_settings,  only : T_LIKE_SET, T_DATA
    use run_info,       only : T_SAMPLE, init_samples, T_RUN_INFO, write_samples
    use m_mcmc,         only : mcmc 
    use m_hmcmc,        only : hmcmc, T_HMC_SAMPLES
    use read_write,     only : write_likelihood, write_number_of_cells, write_temperatures
    use m_initialise,   only : initialise
    

#ifdef MPI
    use mpi
    use m_mcmc,         only : mcmc_pt 
#endif

#ifdef NETCDF
    use netcdf_read_write, only : netcdf_write, T_NC
#endif
            
    implicit none

    ! variable relevant to data
    type(T_DATA), dimension(:), allocatable    :: dat

    ! likelihood settings
    type(T_LIKE_SET)    :: like_set

    ! mcmc settings
    type(T_MCMC_SET)    :: mcmc_set

    ! grid tesselation
    type(T_GRID) :: grid

    ! mcmc
    type(T_SAMPLE), dimension(:), allocatable :: samples
    type(T_HMC_SAMPLES) :: hmcsamples
    type(T_RUN_INFO)  :: RTI
    real(kind=ii10), dimension(:,:), allocatable :: temperatures

    ! file unit
    integer :: iunit

    ! cpu time
    real(kind=ii10) :: t1, t2
    integer(8)      :: count1, count2
    integer(8)      :: count_rate, count_max

    ! mpi
#ifdef MPI
    integer rank
    integer nbproc, ierr
#endif

    ! netcdf specific variables
#ifdef NETCDF
    type(T_NC) :: nc_model
#endif

    character(8) :: date
    character(10) :: time

    namelist /grid_settings/ grid
    namelist /likelihood_settings/ like_set
    namelist /mcmc_settings/ mcmc_set

    !
    ! Parallelization using mpi
    !
#ifdef MPI
    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, nbproc, ierr)
    call mpi_comm_rank(mpi_comm_world, rank, ierr)
#endif

    !
    ! Set up log file
    !
#ifdef MPI
    if(rank==0)then
        if( create_dir(trim(log_dir)) == -1)&
            call exception_raiseError('Error creating the log folder &
            &'//trim(log_dir))
    endif
#else
    if( create_dir(trim(log_dir)) == -1)&
        call exception_raiseError('Error creating the log folder &
        &'//trim(log_dir))
#endif

#ifdef MPI
    call mpi_barrier(mpi_comm_world,ierr)
    call log_startup(trim(log_dir)//FILEDIR_SEPARATOR//'MCTomo_'//itoa(rank+1)//'.log')
#else
    call log_startup(trim(log_dir)//FILEDIR_SEPARATOR//'MCTomo.log')
#endif

    call date_and_time(date,time)
    call log_msg('Start MCTomo at '//date//' '//time//'...')

    ! initialize system clock
    call system_clock(count1,count_rate,count_max)

    ! 
    ! Read settings for mcmc and likelihood
    !
    call log_msg('Reading input file...')
    open(unit=newunit(iunit),file='MCTomo.inp',status='old',delim='apostrophe')
    read(iunit, nml=grid_settings)
    read(iunit, nml=likelihood_settings)
    read(iunit, nml=mcmc_settings)
    close(iunit)
    ! initial some variable of settings
    call grid_setup(grid)
    mcmc_set%grid = grid
    like_set%grid = grid
    like_set%sigdep = mcmc_set%sigdep
    mcmc_set%datatype = like_set%datatype
#ifdef MPI
    if(nbproc>1) mcmc_set%processor = rank + 1
#endif

    !
    ! initialise each chain
    !
    call initialise(RTI,dat,mcmc_set,like_set)


    !
    ! start rj-mcmc sampling
    !
    call log_msg('Start rjmcmc sampling...')
    call system_clock(count2)
    t1 = (count2-count1)*1.0/count_rate
    if(mcmc_set%tempering == 0 .and. mcmc_set%hmc == 0)then
        allocate( samples(RTI%nsamples) )
        call init_samples(samples)
        call mcmc( samples, dat, RTI, mcmc_set, like_set )
    elseif(mcmc_set%tempering == 1)then
        allocate( samples(RTI%nsamples) )
        call init_samples(samples)
        allocate( temperatures(4,RTI%nsamples) )
        temperatures = 0
#ifdef MPI
        call mcmc_pt( samples, temperatures, dat, RTI, mcmc_set, like_set )
#endif
    elseif(mcmc_set%hmc == 1)then
        call hmcmc(hmcsamples, dat, RTI, mcmc_set, like_set)
    endif
    call system_clock(count2)
    t2 = (count2-count1)*1.0/count_rate
    ! end rj-mcmc
    call log_msg('Accepted samples: '//itoa(sum(RTI%acceptedcount)) //&
            '/' // itoa(sum(RTI%samplecount)) )
    call log_msg('Sampling of ' // itoa(mcmc_set%nsamples) // ' samples: '// rtoa(t2-t1) )

    !
    ! write out result
    !
    ! write out the samples
    call write_samples(trim(results_dir)//FILEDIR_SEPARATOR//'samples_'//itoa(mcmc_set%processor)//'.out',samples)
    call write_likelihood(trim(results_dir)//FILEDIR_SEPARATOR//'likelihood_'//itoa(mcmc_set%processor)//'.dat',samples)
    call write_number_of_cells(trim(results_dir)//FILEDIR_SEPARATOR//'ncells_'//itoa(mcmc_set%processor)//'.dat',samples)


#ifdef NETCDF
    call netcdf_write(trim(results_dir)//FILEDIR_SEPARATOR//'samples_'//itoa(mcmc_set%processor)//'.nc',samples)
#endif

    if(mcmc_set%tempering == 1)then
        call write_temperatures(trim(results_dir)//FILEDIR_SEPARATOR//'temperatures_'//itoa(mcmc_set%processor)//'.dat',temperatures)
    endif

    call system_clock(count2)
    t2 = (count2-count1)*1.0/count_rate
    call log_msg('Time taken by the code was ' // rtoa(t2) //&
        ' seconds')

    call log_shutdown()

#ifdef MPI
    call mpi_finalize(ierr)
#endif

end program
