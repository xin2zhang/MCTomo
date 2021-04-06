!
! modules m_mcmc
! reversible jump Markov Chain Monte Carlo algorithm
! 
module m_mcmc

    use m_exception, only       : exception_raiseError
    use m_logger, only          : log_msg
    use m_utils, only           : ii10, itoa, rtoa, delete_file, vs2vp, vp2rho, last_dir, resume_dir, results_dir, FILEDIR_SEPARATOR
    use m_settings, only        : T_MCMC_SET, T_GRID, T_MOD, mod_setup
    use like_settings, only     : T_LIKE_SET, T_DATA, T_LIKE_BASE, likeBaseSetup, write_times
    use m_likelihood, only      : T_LIKE, likelihood, noise_likelihood,like_setup
    use run_info, only          : T_RUN_INFO, T_SAMPLE, write_info, write_samples
    use mt19937, only           : gasdev, unirand
    use kdtree2_precision_module, only : kdkind
    use kdtree2_module, only    : kdtree2, kdtree2_result, kdtree2_create,&
                                  kdtree2_n_nearest
    use read_write, only        : write_mean, write_var, write_vertices, write_temperatures, write_array_1d, write_data

    use sliceSample, only       : slice_sample, slice_sample_pos
    ! use all the available procedure, types and variable in the following
    ! module
    use iso_c_binding
#ifdef MPI
    use mpi
#endif

    implicit none
    private

    public :: mcmc
#ifdef MPI
    public :: mcmc_pt
#endif

    public :: kdtree_to_grid
    public :: ImportanceSample
    public :: proposeSigmaType
    public :: cell_birth, cell_death
    public :: cell_move
    public :: ssigma_change
    public :: stat_rti
   
    ! static value
    real(kind=ii10), parameter :: one = 1.0
    real(kind=ii10), parameter :: PII = 3.1415926
    real(kind=ii10), parameter :: VPVS = 1.732
    real(kind=ii10), parameter :: eps = 1e-8
    integer, parameter :: maxTryNum = 1000

    ! debug
    real t1, t2
    logical debug_mode


contains

    !
    ! subroutine mcmc
    ! implement the reversible jump Markov Chain Monte Carlo mtheod
    ! input  :
    !   data
    !   RTI: run time information
    !   likelihood settings
    !   settings
    ! output :
    !   samples
    !
    subroutine mcmc(samples,dat,RTI,mcmc_set,like_set)

        implicit none

        type(T_SAMPLE), dimension(:), intent(inout)     :: samples
        type(T_DATA),dimension(:),  intent(in)          :: dat
        type(T_RUN_INFO), intent(inout)                 :: RTI
        type(T_MCMC_SET), intent(in)                    :: mcmc_set
        type(T_LIKE_SET), intent(in)                    :: like_set

        ! local variables
        ! voronoi tessellation
        integer                                         :: ncells_copy
        real(ii10), dimension(2,mcmc_set%ncell_max)     :: points_copy
        real(ii10), dimension(3,mcmc_set%ncell_max)     :: parameters_copy
        real(c_double), dimension(:), allocatable       :: snoise0_copy, snoise1_copy
        type(T_GRID)                                    :: grid
        type(T_MOD)                                     :: model, model_copy
        integer, dimension(:,:), allocatable            :: sites_id_copy

        ! propose related variable
        type( T_LIKE )    :: like, like_copy
        real( kind=8 )    :: propose
        integer           :: ptype, icycle
        real( kind=ii10 ) :: prob
        real( kind=ii10 ) :: alpha, random
        integer(c_int)    :: iremove, imove, ivalue, iperiod, iloc
        logical           :: accepted, lerr
        real(kind=ii10),dimension(2)   :: point
        real(kind=ii10),dimension(3)   :: pm_src, pm
        integer           :: iter, countFM, nthin
        integer           :: i, j
        integer           :: accDepMv, totalDepMv, accDepV, totalDepV
        integer           :: accUpMv, totalUpMv, accUpV, totalUpV

        ! debug
        integer, dimension(RTI%nsamples) :: likecounts
        !character(100) :: debug_info

        ! delelte run-time sample file
        call delete_file(trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_sample_'//itoa(mcmc_set%processor)//'.dat')
        ! initialize and allocate several variables
        grid = mcmc_set%grid
        call mod_setup(model, grid)
        call mod_setup(model_copy, grid)

        ! initialize points and parameters corresponding to the current delaunay

        ! calculate likelihood of the initial sample
        call kdtree_to_grid(RTI, grid, model)

        call like_setup(like, dat, like_set, mcmc_set%ncell_max)
        call like_setup(like_copy, dat, like_set, mcmc_set%ncell_max)
        call likelihood(dat,RTI,like_set,like)

        call log_msg('Initial number of cells: '// itoa(RTI%ncells) )
        call log_msg('Initial -loglikelihood: ' // &
            &rtoa(like%like) )
        call log_msg('Initial misfits: ' // &
            &rtoa(like%misfit) )
        call log_msg('Initial unweighted misfits: ' // &
            &rtoa(like%unweighted_misfit) )
        
        !
        ! Markov Chain Monte Carlo
        !

        ! > initialize some setting values
        accDepMv = 0
        totalDepMv = 0
        accDepV = 0
        totalDepV = 0
        accUpMv = 0
        totalUpMv = 0
        accUpV = 0
        totalUpV = 0
        iter = 0
        icycle = 0
        !allocate( samples(mcmc_set%nsamples) )
        !call init_sample(samples) ! initialize sample
        do while(iter < RTI%nsamples)

            iter = iter + 1
            RTI%sampletotal = RTI%sampletotal + 1

            ! back up some values, like delaunay triangualtion, vp, vs, rho,
            ! vel, likelihood, these will possibley be recovered if new sample
            ! not accepteded
            ncells_copy = RTI%ncells
            points_copy = RTI%points
            parameters_copy = RTI%parameters
            snoise0_copy = RTI%snoise0
            snoise1_copy = RTI%snoise1
            sites_id_copy = RTI%sites_id
            model_copy = model
            like_copy = like

            ! initial some values, eg. logical accepted, propose type, random number
            accepted = .false.
            if(icycle ==0 .or. countFM>=maxTryNum)then
                random = unirand(RTI%randcount)
                propose = unirand(RTI%randcount)
                ptype = proposeType2(propose,mcmc_set,iter)
                countFM = 0
            endif

            icycle = 0
            select case ( ptype )
            case(1) 
                ! count
                RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
                samples(iter)%step = ptype

                call cell_birth(RTI,mcmc_set,point,pm,prob,lerr)
                if(.not.lerr) goto 100
                call kdtree_to_grid(RTI,grid,model)
                if(all(RTI%sites_id==sites_id_copy)) goto 100

                ! calculate travel time and likelihood value
                !call cpu_time(t1)
                call likelihood(dat,RTI,like_set,like)
                !call cpu_time(t2)
                !write(*,*) 'likelihood: ', t2-t1

                ! calculate acceptance ratio
                !if(like%like<0) goto 100
                if(abs(like%like-huge(like%like))<eps)then
                    iter = iter - 1
                    RTI%sampletotal = RTI%sampletotal - 1
                    RTI%samplecount(ptype) = RTI%samplecount(ptype) - 1
                    icycle = 1
                    countFM = countFM + 1
                    goto 100
                endif

                if(mcmc_set%kernel==0) then
                    alpha = minval([log(one),prob-&
                        log(mcmc_set%vpmax-mcmc_set%vpmin)-like%like+like_copy%like])
                else
                    alpha = minval([log(one), -like%like+like_copy%like])
                endif
                if(log(random)<alpha)then
                    accepted = .true.
                    ! update accepted ratio and samples
                    samples(iter)%step = ptype
                    samples(iter)%accepted = accepted
                    samples(iter)%vindex = 0
                    samples(iter)%ncells = RTI%ncells
                    samples(iter)%like = like%like
                    samples(iter)%misfit = like%misfit
                    samples(iter)%unweighted_misfit = like%unweighted_misfit
                    samples(iter)%noise0 = 0
                    samples(iter)%noise1 = 0
                    samples(iter)%coord = point
                    samples(iter)%values = pm
                    RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                endif
                ! debug

            case(2)
                ! count
                RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
                samples(iter)%step = ptype
                ! remove a cell
                call cell_death(RTI,mcmc_set,iremove,pm,prob,lerr)
                if(.not.lerr) goto 100
                call kdtree_to_grid(RTI,grid,model)

                ! calculate travel time and likelihood value
                call likelihood(dat,RTI,like_set,like)

                ! calculate acceptance ratio
                !if(like%like<0) goto 100
                if(abs(like%like-huge(like%like))<eps)then
                    iter = iter - 1
                    RTI%sampletotal = RTI%sampletotal - 1
                    RTI%samplecount(ptype) = RTI%samplecount(ptype) - 1
                    icycle = 1
                    countFM = countFM + 1
                    goto 100
                endif

                if(mcmc_set%kernel==0) then
                    alpha = minval([log(one),prob + &
                        log(mcmc_set%vpmax-mcmc_set%vpmin)-like%like+like_copy%like])
                else
                    alpha = minval([log(one), -like%like+like_copy%like])
                endif
                if(log(random)<alpha)then
                    accepted = .true.
                    ! update accepted ratio and samples
                    samples(iter)%step = ptype
                    samples(iter)%accepted = accepted
                    samples(iter)%vindex = iremove
                    samples(iter)%ncells = RTI%ncells
                    samples(iter)%like = like%like
                    samples(iter)%misfit = like%misfit
                    samples(iter)%unweighted_misfit = like%unweighted_misfit
                    samples(iter)%noise0 = 0
                    samples(iter)%noise1 = 0
                    samples(iter)%coord = 0
                    samples(iter)%values = 0
                    RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                endif
                ! debug

            case(3)
                ! count
                RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
                samples(iter)%step = ptype
                if(mcmc_set%slicesample==2)then
                    imove = ceiling(unirand(RTI%randcount)*RTI%ncells)
                    if( imove==0 ) imove = 1
                    call slice_sample_pos(dat,RTI,like_set,mcmc_set,1,imove,like,mcmc_set%pd/100*(grid%xmax-grid%xmin))
                    call slice_sample_pos(dat,RTI,like_set,mcmc_set,2,imove,like,mcmc_set%pd/100*(grid%ymax-grid%ymin))
                    do i = 1, grid%nx
                        do j = 1, grid%ny
                            model%vp(j,i) = RTI%parameters(1,RTI%sites_id(j,i))
                            model%vs(j,i) = RTI%parameters(2,RTI%sites_id(j,i))
                            model%rho(j,i) = RTI%parameters(3,RTI%sites_id(j,i))
                        enddo
                    enddo
                    accepted=.true.
                else
                    ! move a cell
                    call cell_move(RTI,mcmc_set,imove,point,lerr)
                    ! convert voronoi to grid
                    if(lerr) then
                        call kdtree_to_grid(RTI,grid,model)
                        ! calculate travel time and likelihood value
                        call likelihood(dat,RTI,like_set,like)
                        ! calculate acceptance ratio
                        !if(like%like<0) goto 100
                        if(abs(like%like-huge(like%like))<eps)then
                            iter = iter - 1
                            RTI%sampletotal = RTI%sampletotal - 1
                            RTI%samplecount(ptype) = RTI%samplecount(ptype) - 1
                            icycle = 1
                            countFM = countFM + 1
                            goto 100
                        endif
                        ! calculate acceptance ratio
                        alpha = minval([log(one), -like%like+like_copy%like])
                        if( log(random)<alpha )then
                            accepted = .true.
                            samples(iter)%step = ptype
                        endif
                    endif
                endif

                if(mcmc_set%datatype==2 .and. points_copy(2,imove)>(grid%ymin+grid%ymax)/2)then
                    if(accepted) accDepMv = accDepMv + 1
                    totalDepMv = totalDepMv + 1
                else
                    if(accepted) accUpMv = accUpMv + 1
                    totalUpMv = totalUpMv + 1
                endif
                if( accepted )then
                    samples(iter)%accepted = accepted
                    samples(iter)%vindex = imove
                    samples(iter)%ncells = RTI%ncells
                    samples(iter)%like = like%like
                    samples(iter)%misfit = like%misfit
                    samples(iter)%unweighted_misfit = like%unweighted_misfit
                    samples(iter)%noise0 = 0
                    samples(iter)%noise1 = 0
                    samples(iter)%coord = point
                    samples(iter)%values = 0
                    RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                endif

                ! debug
                ! if not accepted firstly, start delayed rejection
            case(4)
                ! count
                RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
                samples(iter)%step = ptype
                if(mcmc_set%slicesample>=1)then
                    ivalue = ceiling(unirand(RTI%randcount)*RTI%ncells)
                    if( ivalue==0 ) ivalue = 1
                    call slice_sample(dat,RTI,like_set,mcmc_set,1,ivalue,like,mcmc_set%sigma_vp)
                    do i = 1, grid%nx
                        do j = 1, grid%ny
                            model%vp(j,i) = RTI%parameters(1,RTI%sites_id(j,i))
                            model%vs(j,i) = RTI%parameters(2,RTI%sites_id(j,i))
                            model%rho(j,i) = RTI%parameters(3,RTI%sites_id(j,i))
                        enddo
                    enddo
                    accepted=.true.
                else
                    ! change a velocity
                    call vel_change(RTI,mcmc_set,ivalue,pm_src,pm,lerr)
                    ! accept or not based on acceptance ratio
                    if(lerr) then
                        ! update model
                        do i = 1, grid%nx
                            do j = 1, grid%ny
                                 model%vp(j,i) = RTI%parameters(1,RTI%sites_id(j,i))
                                 model%vs(j,i) = RTI%parameters(2,RTI%sites_id(j,i))
                                 model%rho(j,i) = RTI%parameters(3,RTI%sites_id(j,i))
                            enddo
                        enddo
                        ! calculate travel time and likelihood value
                        call likelihood(dat,RTI,like_set,like)
                        ! calculate acceptance ratio
                        !if(like%like<0) goto 100
                        if(abs(like%like-huge(like%like))<eps)then
                            iter = iter - 1
                            RTI%sampletotal = RTI%sampletotal - 1
                            RTI%samplecount(ptype) = RTI%samplecount(ptype) - 1
                            icycle = 1
                            countFM = countFM + 1
                            goto 100
                        endif
                        ! calculate acceptance ratio
                        alpha = minval([log(one), -like%like+like_copy%like])
                        if( log(random)<alpha )then
                            accepted = .true.
                            samples(iter)%step = ptype
                        endif
                    endif
                endif
                if(mcmc_set%datatype==2 .and. points_copy(2,ivalue)>(grid%ymin+grid%ymax)/2)then
                    if(accepted) accDepV = accDepV + 1
                    totalDepV = totalDepV + 1
                else
                    if(accepted) accUpV = accUpV + 1
                    totalUpV = totalUpV + 1
                endif
                ! accept
                if( accepted )then
                    accepted = .true.
                    samples(iter)%step = ptype
                    samples(iter)%accepted = accepted
                    samples(iter)%vindex = ivalue
                    samples(iter)%ncells = RTI%ncells
                    samples(iter)%like = like%like
                    samples(iter)%misfit = like%misfit
                    samples(iter)%unweighted_misfit = like%unweighted_misfit
                    samples(iter)%noise0 = 0
                    samples(iter)%noise1 = 0
                    !samples(iter)%coord = [point%x,point%y,point%z]
                    samples(iter)%coord = 0
                    samples(iter)%values = RTI%parameters(:,ivalue)
                    RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                endif
                ! if not accepted, delay reject
                ! debug
            case(5)
                ! count
                RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
                samples(iter)%step = ptype
                ! change the data noise sigma
                call ssigma_change(RTI,mcmc_set,iperiod,lerr)
                if(.not.lerr) goto 100
                ! calculate likelihood
                call noise_likelihood(dat,RTI,like_set,like)
                ! calculate acceptance ratio
                alpha = minval([log(one), -like%like+like_copy%like])
                ! if accepted
                if(log(random)<alpha)then
                    accepted = .true.
                    ! update accepted ratio and samples
                    samples(iter)%step = ptype
                    samples(iter)%accepted = accepted
                    samples(iter)%vindex = iperiod
                    samples(iter)%ncells = RTI%ncells
                    samples(iter)%like = like%like
                    samples(iter)%misfit = like%misfit
                    samples(iter)%unweighted_misfit = like%unweighted_misfit
                    samples(iter)%noise0 = RTI%snoise0(iperiod)
                    samples(iter)%noise1 = RTI%snoise1(iperiod)
                    samples(iter)%coord = 0
                    samples(iter)%values = 0
                    RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                endif

            case default 
                call exception_raiseError("wrong proposal type (1 for birth, 2&
                    &for death, 3 for move, 4 for velocity change, 5 for surface wave noise change, 6 for body wave noise change&
                    &and 7 for source locations change)")
              end select
            ! exit of select structure, might be changed to exit instead of goto
            ! in fortran 2008
            100 continue     
            ! if accepteded, update delaunay triangulation, vp, vs, rho, vel,
            ! like
            if(.not.accepted) then
                model = model_copy
                like = like_copy
                RTI%ncells = ncells_copy
                RTI%points = points_copy
                RTI%parameters = parameters_copy
                RTI%snoise0 = snoise0_copy
                RTI%snoise1 = snoise1_copy
                RTI%sites_id = sites_id_copy
            endif

            ! if icylce = 1, redoing the step thus skip next section
            if(icycle==1) cycle

            samples(iter)%ncells = RTI%ncells
            samples(iter)%like = like%like
            samples(iter)%misfit = like%misfit
            samples(iter)%unweighted_misfit = like%unweighted_misfit

            likecounts(iter) = RTI%like_count

            ! if after burn-in, begin sampling
            if(RTI%sampletotal > mcmc_set%burn_in .and.&
            mod(RTI%sampletotal-mcmc_set%burn_in,mcmc_set%thin)==0)then
                if(RTI%sampletotal-mcmc_set%burn_in==mcmc_set%thin)then
                    call log_msg('###########################')
                    call log_msg('Start sampling here: '//itoa(RTI%sampletotal))
                endif
                call stat_rti(RTI,model,mcmc_set)
            endif
            ! write run-time info into files
            if(mcmc_set%runtime_step > 0 .and. mod(iter,mcmc_set%runtime_step) == 0)then
                call write_samples(trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_sample_'//itoa(mcmc_set%processor)&
                    &//'.dat',samples(RTI%nsampled+1:iter))
                call write_array_1d(likecounts(RTI%nsampled+1:iter),trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_likecounts_'//itoa(mcmc_set%processor)//'.dat')
                RTI%nsampled = iter
                call write_info(RTI,trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_info_'//itoa(mcmc_set%processor)//'.dat')
                call write_vertices(RTI,trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_vertices_'//itoa(mcmc_set%processor))
                call write_mean(RTI,trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_average_'//itoa(mcmc_set%processor)//'.dat')
                call write_var(RTI,trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_var_'//itoa(mcmc_set%processor)//'.dat')
                call write_data(trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_predicated_data_'//itoa(mcmc_set%processor),dat,RTI,like%likelihoods,like_set)
            endif
            ! write info during run
            ! display some information of the markov chain
            if(mod(iter,mcmc_set%display) == 0)then
                call log_msg('')
                call log_msg( 'Processor number: '//itoa(mcmc_set%processor) )
                call log_msg( 'Sample: '//itoa(iter)// '/'//itoa(RTI%nsamples) )
                call log_msg( 'Number of cells: '//itoa(int(RTI%ncells)))
                call log_msg( '-loglikelihood: '//rtoa(like%like) )
                call log_msg( 'Weighted misfit: '//rtoa(like%misfit) )
                call log_msg( 'Unweighted misfit: '//rtoa(like%unweighted_misfit) )
                call log_msg('')
                call log_msg( 'Accpetance rate for velocity' )
                call log_msg( 'ARV '//rtoa(accUpV/real(totalUpV,kind=ii10)) )
                call log_msg( 'Accpetance rate for position' )
                call log_msg( 'ARP '//rtoa(accUpMv/real(totalUpMv,kind=ii10)) )
                !call log_msg( 'ARV '//rtoa(RTI%acceptedcount(4)/real(RTI%samplecount(4),kind=ii10)) )
                !call log_msg( 'Accpetance rate for position' )
                !call log_msg( 'ARP '//rtoa(RTI%acceptedcount(3)/real(RTI%samplecount(3),kind=ii10)) )
                call log_msg( 'Accpetance rate for deep velocity' )
                call log_msg( 'ARV '//rtoa(accDepV/real(totalDepV,kind=ii10)) )
                call log_msg( 'Accpetance rate for deep position' )
                call log_msg( 'ARP '//rtoa(accDepMv/real(totalDepMv,kind=ii10)) )
                call log_msg( 'Accpetance rate for birth' )
                call log_msg( 'ARB '//rtoa(RTI%acceptedcount(1)/real(RTI%samplecount(1),kind=ii10)) )
                call log_msg( 'Accpetance rate for death' )
                call log_msg( 'ARD '//rtoa(RTI%acceptedcount(2)/real(RTI%samplecount(2),kind=ii10)) )
                call log_msg( 'Accpetance rate for sigma change of body waves' )
                call log_msg( 'ARS '//rtoa(RTI%acceptedcount(5)/real(RTI%samplecount(5),kind=ii10)) )
                call log_msg( 'Likelihood calculation count' )
                call log_msg( 'Number of likelihood for slice sampling '//itoa(RTI%like_count))
                call log_msg( 'Rate of bad models' )
                call log_msg( 'Bad dispersion curves '//rtoa(RTI%num_bad_model/real(RTI%sampletotal,kind=ii10)) )
                call log_msg( 'Bad rays '//rtoa(RTI%num_bad_ray/real(RTI%sampletotal,kind=ii10)) )
                call log_msg( '---------------------------------------------------------' )
            endif
        enddo

        ! write out the last delaunay triangulation
        ! write both a text and binary file, so no suffix specified
        call write_vertices(RTI,trim(last_dir)//FILEDIR_SEPARATOR//'last_vertices_'//itoa(mcmc_set%processor) )
        ! write last parameters to file
        call write_info(RTI,trim(last_dir)//FILEDIR_SEPARATOR//'last_info_'//itoa(mcmc_set%processor)//'.dat')
        call write_mean(RTI,trim(last_dir)//FILEDIR_SEPARATOR//'last_average_'//itoa(mcmc_set%processor)//'.dat')
        call write_var(RTI,trim(last_dir)//FILEDIR_SEPARATOR//'last_var_'//itoa(mcmc_set%processor)//'.dat')
        call write_array_1d(likecounts,trim(results_dir)//FILEDIR_SEPARATOR//'likecounts_'//itoa(mcmc_set%processor)//'.dat')
        return

    end subroutine mcmc

#ifdef MPI
    subroutine mcmc_pt(samples,temperatures,dat,RTI,mcmc_set,like_set)

        implicit none

        type(T_SAMPLE), dimension(:), intent(inout)     :: samples
        real(kind=ii10), dimension(:,:), intent(inout)  :: temperatures
        type(T_DATA), dimension(:), intent(in)          :: dat
        type(T_RUN_INFO), intent(inout)                 :: RTI
        type(T_MCMC_SET), intent(in)                    :: mcmc_set
        type(T_LIKE_SET), intent(in)                    :: like_set

        ! local variables
        ! voronoi tessellation
        integer                                         :: ncells_copy
        real(ii10), dimension(2,mcmc_set%ncell_max)     :: points_copy
        real(ii10), dimension(3,mcmc_set%ncell_max)     :: parameters_copy
        real(c_double), dimension(:), allocatable       :: snoise0_copy, snoise1_copy
        type(T_GRID)                                    :: grid
        type(T_MOD)                                     :: model, model_copy
        integer, dimension(:,:), allocatable            :: sites_id_copy

        ! propose related variable
        type( T_LIKE )    :: like, like_copy
        real( kind=8 )    :: propose
        integer           :: ptype, icycle
        real( kind=ii10 ) :: prob
        real( kind=ii10 ) :: alpha, random
        integer(c_int)    :: iremove, imove, ivalue, iperiod, iloc
        logical           :: accepted, lerr
        real(kind=ii10),dimension(2)   :: point
        real(kind=ii10),dimension(3)   :: pm_src, pm
        integer           :: iter, countFM, nthin
        integer           :: i, j

        ! parallel tempering
        real(kind=ii10) :: temperature
        integer         :: t1_chain, t2_chain
        real(kind=ii10) :: likelihood_t1, likelihood_t2

        ! mpi
        integer ierror, astatus
        integer status(MPI_STATUS_SIZE)

        ! debug
        integer, dimension(RTI%nsamples) :: likecounts

        ! debug
        debug_mode = .false.
        !character(100) :: debug_info
        temperature = 1

        ! delelte run-time sample file
        call delete_file(trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_sample_'//itoa(mcmc_set%processor)//'.dat')
        call delete_file(trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_temperatures_'//itoa(mcmc_set%processor)//'.dat')
        ! initialize and allocate several variables
        grid = mcmc_set%grid
        call mod_setup(model, grid)
        call mod_setup(model_copy, grid)

        ! initialize tempering variable
        likelihood_t1 = 0.0
        likelihood_t2 = 0.0

        ! initialize points and parameters corresponding to the current delaunay

        ! calculate likelihood of the initial sample
        call kdtree_to_grid(RTI, grid, model)

        call like_setup(like, dat, like_set, mcmc_set%ncell_max)
        call like_setup(like_copy, dat, like_set, mcmc_set%ncell_max)
        call likelihood(dat,RTI,like_set,like)

        call log_msg('Initial number of cells: '// itoa(RTI%ncells) )
        call log_msg('Initial -loglikelihood: ' // &
            &rtoa(like%like) )
        call log_msg('Initial misfits: ' // &
            &rtoa(like%misfit) )
        call log_msg('Initial unweighted misfits: ' // &
            &rtoa(like%unweighted_misfit) )
        
        !
        ! Markov Chain Monte Carlo
        !

        ! > initialize some setting values
        iter = 0
        icycle = 0
        !allocate( samples(mcmc_set%nsamples) )
        !call init_sample(samples) ! initialize sample
        do while(iter < RTI%nsamples)

            iter = iter + 1
            RTI%sampletotal = RTI%sampletotal + 1
            if(RTI%sampletotal>=mcmc_set%tempering_start)then
                temperature = RTI%temperature_values(mcmc_set%processor)
            else
                temperature = 1
            endif


            ! back up some values, like delaunay triangualtion, vp, vs, rho,
            ! vel, likelihood, these will possibley be recovered if new sample
            ! not accepteded
            ncells_copy = RTI%ncells
            points_copy = RTI%points
            parameters_copy = RTI%parameters
            snoise0_copy = RTI%snoise0
            snoise1_copy = RTI%snoise1
            sites_id_copy = RTI%sites_id
            model_copy = model
            like_copy = like

            ! initial some values, eg. logical accepted, propose type, random number
            accepted = .false.
            if(icycle ==0 .or. countFM>=maxTryNum)then
                random = unirand(RTI%randcount)
                propose = unirand(RTI%randcount)
                ptype = proposeType2(propose,mcmc_set,iter)
                countFM = 0
            endif

            icycle = 0
            select case ( ptype )
            case(1) 
                ! count
                RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
                samples(iter)%step = ptype

                call cell_birth(RTI,mcmc_set,point,pm,prob,lerr)
                if(.not.lerr) goto 100
                call kdtree_to_grid(RTI,grid,model)
                if(all(RTI%sites_id==sites_id_copy)) goto 100

                ! calculate travel time and likelihood value
                !call cpu_time(t1)
                call likelihood(dat,RTI,like_set,like)
                !call cpu_time(t2)
                !write(*,*) 'likelihood: ', t2-t1

                ! calculate acceptance ratio
                !if(like%like<0) goto 100
                if(abs(like%like-huge(like%like))<eps)then
                    iter = iter - 1
                    RTI%sampletotal = RTI%sampletotal - 1
                    RTI%samplecount(ptype) = RTI%samplecount(ptype) - 1
                    countFM = countFM + 1
                    icycle = 1
                    goto 100
                endif

                if(mcmc_set%kernel==0) then
                    alpha = minval([log(one),prob-&
                        log(mcmc_set%vsmax-mcmc_set%vsmin)+(-like%like+like_copy%like)/temperature])
                else
                    alpha = minval([log(one),&
                    (-like%like+like_copy%like)/temperature])
                endif
                if(log(random)<alpha)then
                    accepted = .true.
                    ! update accepted ratio and samples
                    samples(iter)%step = ptype
                    samples(iter)%accepted = accepted
                    samples(iter)%vindex = 0
                    samples(iter)%ncells = RTI%ncells
                    samples(iter)%like = like%like
                    samples(iter)%misfit = like%misfit
                    samples(iter)%unweighted_misfit = like%unweighted_misfit
                    samples(iter)%noise0 = 0
                    samples(iter)%noise1 = 0
                    samples(iter)%coord = point
                    samples(iter)%values = pm
                    RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                endif
                ! debug

            case(2)
              ! count
              RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
              samples(iter)%step = ptype

              ! remove a cell
              call cell_death(RTI,mcmc_set,iremove,pm,prob,lerr)
              if(.not.lerr) goto 100
              call kdtree_to_grid(RTI,grid,model)

              ! calculate travel time and likelihood value
              call likelihood(dat,RTI,like_set,like)

              ! calculate acceptance ratio
              !if(like%like<0) goto 100
              if(abs(like%like-huge(like%like))<eps)then
                  iter = iter - 1
                  RTI%sampletotal = RTI%sampletotal - 1
                  RTI%samplecount(ptype) = RTI%samplecount(ptype) - 1
                  icycle = 1
                  countFM = countFM + 1
                  goto 100
              endif

              if(mcmc_set%kernel==0) then
                  alpha = minval([log(one),prob + &
                      log(mcmc_set%vsmax-mcmc_set%vsmin)+(-like%like+like_copy%like)/temperature])
              else
                  alpha = minval([log(one),&
                  (-like%like+like_copy%like)/temperature])
              endif
              if(log(random)<alpha)then
                  accepted = .true.
                  ! update accepted ratio and samples
                  samples(iter)%step = ptype
                  samples(iter)%accepted = accepted
                  samples(iter)%vindex = iremove
                  samples(iter)%ncells = RTI%ncells
                  samples(iter)%like = like%like
                  samples(iter)%misfit = like%misfit
                  samples(iter)%unweighted_misfit = like%unweighted_misfit
                  samples(iter)%noise0 = 0
                  samples(iter)%noise1 = 0
                  samples(iter)%coord = 0
                  samples(iter)%values = 0
                  RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
              endif
              ! debug

            case(3)
              ! count
              RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
              samples(iter)%step = ptype
              ! move a cell
              call cell_move(RTI,mcmc_set,imove,point,lerr)
              ! convert voronoi to grid
              if(lerr) then
                  call kdtree_to_grid(RTI,grid,model)
                  ! calculate travel time and likelihood value
                  call likelihood(dat,RTI,like_set,like)
                  ! calculate acceptance ratio
                  !if(like%like<0) goto 100
                  if(abs(like%like-huge(like%like))<eps)then
                      iter = iter - 1
                      RTI%sampletotal = RTI%sampletotal - 1
                      RTI%samplecount(ptype) = RTI%samplecount(ptype) - 1
                      icycle = 1
                      countFM = countFM + 1
                      goto 100
                  endif
                  ! calculate acceptance ratio
                  alpha = minval([log(one),&
                  (-like%like+like_copy%like)/temperature])
                  if( log(random)<alpha )then
                      accepted = .true.
                      samples(iter)%step = ptype
                      samples(iter)%accepted = accepted
                      samples(iter)%vindex = imove
                      samples(iter)%ncells = RTI%ncells
                      samples(iter)%like = like%like
                      samples(iter)%misfit = like%misfit
                      samples(iter)%unweighted_misfit = like%unweighted_misfit
                      samples(iter)%noise0 = 0
                      samples(iter)%noise1 = 0
                      samples(iter)%coord = point
                      samples(iter)%values = 0
                      RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                  endif
              endif

              ! debug
              ! if not accepted firstly, start delayed rejection
            case(4)
              ! count
              RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
              samples(iter)%step = ptype
              if(mcmc_set%slicesample>=1)then
                  ivalue = ceiling(unirand(RTI%randcount)*RTI%ncells)
                  if( ivalue==0 ) ivalue = 1
                  call slice_sample(dat,RTI,like_set,mcmc_set,1,ivalue,like,mcmc_set%sigma_vp)
                  do i = 1, grid%nx
                      do j = 1, grid%ny
                          model%vp(j,i) = RTI%parameters(1,RTI%sites_id(j,i))
                          model%vs(j,i) = RTI%parameters(2,RTI%sites_id(j,i))
                          model%rho(j,i) = RTI%parameters(3,RTI%sites_id(j,i))
                      enddo
                  enddo
                  accepted=.true.
              else
                  ! change a velocity
                  call vel_change(RTI,mcmc_set,ivalue,pm_src,pm,lerr)
                  ! accept or not based on acceptance ratio
                  if(lerr) then
                      ! update model
                      do i = 1, grid%nx
                          do j = 1, grid%ny
                               model%vp(j,i) = RTI%parameters(1,RTI%sites_id(j,i))
                               model%vs(j,i) = RTI%parameters(2,RTI%sites_id(j,i))
                               model%rho(j,i) = RTI%parameters(3,RTI%sites_id(j,i))
                          enddo
                      enddo
                      ! calculate travel time and likelihood value
                      call likelihood(dat,RTI,like_set,like)
                      ! calculate acceptance ratio
                      !if(like%like<0) goto 100
                      if(abs(like%like-huge(like%like))<eps)then
                          iter = iter - 1
                          RTI%sampletotal = RTI%sampletotal - 1
                          RTI%samplecount(ptype) = RTI%samplecount(ptype) - 1
                          icycle = 1
                          countFM = countFM + 1
                          goto 100
                      endif
                      ! calculate acceptance ratio
                      alpha = minval([log(one),&
                      (-like%like+like_copy%like)/temperature])
                      if( log(random)<alpha )then
                          accepted = .true.
                          samples(iter)%step = ptype
                      endif
                  endif
              endif
              ! accept
              if( accepted )then
                  accepted = .true.
                  samples(iter)%step = ptype
                  samples(iter)%accepted = accepted
                  samples(iter)%vindex = ivalue
                  samples(iter)%ncells = RTI%ncells
                  samples(iter)%like = like%like
                  samples(iter)%misfit = like%misfit
                  samples(iter)%unweighted_misfit = like%unweighted_misfit
                  samples(iter)%noise0 = 0
                  samples(iter)%noise1 = 0
                  !samples(iter)%coord = [point%x,point%y,point%z]
                  samples(iter)%coord = 0
                  samples(iter)%values = RTI%parameters(:,ivalue)
                  RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
              endif
              ! if not accepted, delay reject
              ! debug
            case(5)
                ! count
                RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
                samples(iter)%step = ptype
                ! change the data noise sigma
                call ssigma_change(RTI,mcmc_set,iperiod,lerr)
                if(.not.lerr) goto 100
                ! calculate likelihood
                call noise_likelihood(dat,RTI,like_set, like)
                ! calculate acceptance ratio
                alpha = minval([log(one),&
                (-like%like+like_copy%like)/temperature])
                ! if accepted
                if(log(random)<alpha)then
                    accepted = .true.
                    ! update accepted ratio and samples
                    samples(iter)%step = ptype
                    samples(iter)%accepted = accepted
                    samples(iter)%vindex = iperiod
                    samples(iter)%ncells = RTI%ncells
                    samples(iter)%like = like%like
                    samples(iter)%misfit = like%misfit
                    samples(iter)%unweighted_misfit = like%unweighted_misfit
                    samples(iter)%noise0 = RTI%snoise0(iperiod)
                    samples(iter)%noise1 = RTI%snoise1(iperiod)
                    samples(iter)%coord = 0
                    samples(iter)%values = 0
                    RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                endif


              case default 
                call exception_raiseError("wrong proposal type (1 for birth, 2&
                    &for death, 3 for move, 4 for velocity change, 5 for surface wave noise change, 6 for body wave noise change&
                    &and 7 for source locations change)")
              end select
            ! exit of select structure, might be changed to exit instead of goto
            ! in fortran 2008
            100 continue     
            ! if not accepteded, resume delaunay triangulation, vp, vs, rho, vel,
            ! like, noise
            if(.not.accepted) then
                model = model_copy
                like = like_copy
                RTI%ncells = ncells_copy
                RTI%points = points_copy
                RTI%parameters = parameters_copy
                RTI%snoise0 = snoise0_copy
                RTI%snoise1 = snoise1_copy
                RTI%sites_id = sites_id_copy
            endif

            ! if icylce = 1, redoing the step thus skip next section
            if(icycle==1) cycle

            samples(iter)%ncells = RTI%ncells
            RTI%ncells = samples(iter)%ncells
            samples(iter)%like = like%like
            samples(iter)%misfit = like%misfit
            samples(iter)%unweighted_misfit = like%unweighted_misfit

            ! at this point, apply parallel tempering
            if(RTI%sampletotal>=mcmc_set%tempering_start .and. &
               mod(RTI%sampletotal-mcmc_set%tempering_start,mcmc_set%tempering_step)==0)then
                if(mcmc_set%processor == 1)then
                    call choose_t1(RTI,t1_chain)
                    call choose_t2(RTI,mcmc_set,t1_chain,t2_chain)
                endif
                CALL MPI_BCAST(T1_chain,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
                CALL MPI_BCAST(T2_chain,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
                
                ! Get likelihood
                if(mcmc_set%processor == t1_chain)then
                    likelihood_t1 = like%like
                endif
                if(mcmc_set%processor == t2_chain)then
                    likelihood_t2 = like%like
                endif

                if (t1_chain.ne.1) then
                    if (mcmc_set%processor.eq.1) then
                            call mpi_recv(likelihood_t1,1,mpi_double_precision,&
                                          t1_chain-1,1,mpi_comm_world,status,ierror)
                    else if (mcmc_set%processor.eq.t1_chain) then
                            call mpi_ssend(likelihood_t1,1,mpi_double_precision,&
                                           0,1,mpi_comm_world,ierror)
                    end if
                end if
                if (t2_chain.ne.1) then
                    if (mcmc_set%processor.eq.1) then
                            call mpi_recv(likelihood_t2,1,mpi_double_precision,&
                                          t2_chain-1,1,mpi_comm_world,status,ierror)
                    else if (mcmc_set%processor.eq.t2_chain) then
                            call mpi_ssend(likelihood_t2,1,mpi_double_precision,&
                                           0,1,mpi_comm_world,ierror)
                    end if
                end if
                
                ! check whether the temperatures are swapping
                if (mcmc_set%processor.eq.1) then
                    call do_tempering(RTI,t1_chain,t2_chain,likelihood_t1,likelihood_t2,astatus)
                    if(astatus==0)then
                        temperatures(1:3,iter) = (/1,t1_chain,t2_chain/)
                    else
                        temperatures(1:3,iter) = (/0,t1_chain,t2_chain/)
                    endif
                end if
                call mpi_bcast(rti%temperature_values,rti%number_of_temperatures,mpi_double_precision,&
                               0,mpi_comm_world,ierror)
                call mpi_bcast(rti%temperature_indices,rti%number_of_temperatures,mpi_int,&
                               0,mpi_comm_world,ierror)
                temperature=rti%temperature_values(mcmc_set%processor)
                temperatures(4,iter) = temperature
                call mpi_barrier(mpi_comm_world,ierror)
            endif
            temperature=rti%temperature_values(mcmc_set%processor)
            temperatures(4,iter) = temperature
            
            ! if after burn-in, begin sampling
            if(RTI%sampletotal > mcmc_set%burn_in .and. &
            mod(RTI%sampletotal-mcmc_set%burn_in,mcmc_set%thin)==0 .and. &
            abs(temperature-1)<eps)then
                if(RTI%sampletotal-mcmc_set%burn_in==mcmc_set%thin)then
                    call log_msg('###########################')
                    call log_msg('Start sampling here: '//itoa(RTI%sampletotal))
                endif
                call stat_rti(RTI,model,mcmc_set)
            endif
            ! write run-time info into files
            if(mcmc_set%runtime_step > 0 .and. mod(iter,mcmc_set%runtime_step) == 0)then
                call write_samples(trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_sample_'//itoa(mcmc_set%processor)&
                    &//'.dat',samples(RTI%nsampled+1:iter))
                call write_array_1d(likecounts(RTI%nsampled+1:iter),trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_likecounts_'//itoa(mcmc_set%processor)//'.dat')
                call write_temperatures(trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_temperatures_'//itoa(mcmc_set%processor)&
                    &//'.dat',temperatures(:,RTI%nsampled+1:iter))
                RTI%nsampled = iter
                call write_info(RTI,trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_info_'//itoa(mcmc_set%processor)//'.dat')
                call write_vertices(RTI,trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_vertices_'//itoa(mcmc_set%processor))
                call write_mean(RTI,trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_average_'//itoa(mcmc_set%processor)//'.dat')
                call write_var(RTI,trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_var_'//itoa(mcmc_set%processor)//'.dat')
                call write_data(trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_predicated_data_'//itoa(mcmc_set%processor),dat,RTI,like%likelihoods,like_set)
            endif
            ! write info during run
            ! display some information of the markov chain
            if(mod(iter,mcmc_set%display) == 0)then
                call log_msg('')
                call log_msg( 'Processor number: '//itoa(mcmc_set%processor) )
                call log_msg( 'Sample: '//itoa(iter)// '/'//itoa(RTI%nsamples) )
                call log_msg( 'Number of cells: '//itoa(RTI%ncells) )
                call log_msg( '-loglikelihood: '//rtoa(like%like) )
                call log_msg( 'Weighted misfit: '//rtoa(like%misfit) )
                call log_msg( 'Unweighted misfit: '//rtoa(like%unweighted_misfit) )
                call log_msg( 'Temperature: '//rtoa(temperature) )
                call log_msg('')
                call log_msg( 'Accpetance rate for velocity' )
                call log_msg( 'ARV '//rtoa(RTI%acceptedcount(4)/real(RTI%samplecount(4),kind=ii10)) )
                call log_msg( 'Accpetance rate for position' )
                call log_msg( 'ARP '//rtoa(RTI%acceptedcount(3)/real(RTI%samplecount(3),kind=ii10)) )
                call log_msg( 'Accpetance rate for birth' )
                call log_msg( 'ARB '//rtoa(RTI%acceptedcount(1)/real(RTI%samplecount(1),kind=ii10)) )
                call log_msg( 'Accpetance rate for death' )
                call log_msg( 'ARD '//rtoa(RTI%acceptedcount(2)/real(RTI%samplecount(2),kind=ii10)) )
                call log_msg( 'Accpetance rate for sigma change ' )
                call log_msg( 'ARS '//rtoa(RTI%acceptedcount(5)/real(RTI%samplecount(5),kind=ii10)) )
                call log_msg( 'Accpetance rate for sigma change of surface waves' )
                call log_msg( 'ARS '//rtoa(RTI%acceptedcount(7)/real(RTI%samplecount(7),kind=ii10)) )
                call log_msg( 'Accpetance rate for tempering' )
                call log_msg( 'ART '//rtoa(RTI%accepted_tempering/real(RTI%total_tempering,kind=ii10)) )
                call log_msg( 'Likelihood calculation count' )
                call log_msg( 'Number of likelihood for slice sampling '//itoa(RTI%like_count))
                call log_msg( 'Rate of bad models' )
                call log_msg( 'Bad dispersion curves '//rtoa(RTI%num_bad_model/real(RTI%sampletotal,kind=ii10)) )
                call log_msg( 'Bad rays '//rtoa(RTI%num_bad_ray/real(RTI%sampletotal,kind=ii10)) )
                call log_msg( '---------------------------------------------------------' )
            endif
        enddo

        ! write out the last delaunay triangulation
        ! write both a text and binary file, so no suffix specified
        call write_vertices(RTI,trim(last_dir)//FILEDIR_SEPARATOR//'last_vertices_'//itoa(mcmc_set%processor) )
        ! write last parameters to file
        call write_info(RTI,trim(last_dir)//FILEDIR_SEPARATOR//'last_info_'//itoa(mcmc_set%processor)//'.dat')
        call write_mean(RTI,trim(last_dir)//FILEDIR_SEPARATOR//'last_average_'//itoa(mcmc_set%processor)//'.dat')
        call write_var(RTI,trim(last_dir)//FILEDIR_SEPARATOR//'last_var_'//itoa(mcmc_set%processor)//'.dat')
        call write_array_1d(likecounts,trim(results_dir)//FILEDIR_SEPARATOR//'likecounts_'//itoa(mcmc_set%processor)//'.dat')

        return

    end subroutine mcmc_pt
#endif

    integer function proposeType2(propose,set,iter)
    
        implicit none
        real(kind=8), intent(in) :: propose
        type(T_MCMC_SET), intent(in) :: set
        integer, intent(in) :: iter

        real(kind=ii10) random
        call random_number(random)

        proposeType2 = 4
        select case (mod(iter,2))
        case ( 0 )
            if (propose < 0.333) then
                proposeType2 = 1
            elseif(propose < 0.666) then
                proposeType2 = 2
            else
                proposeType2 = 3
            endif
        case ( 1 )
            if(set%sigdep /= 0) then
                if (propose < 0.950)then
                    proposeType2 = 4
                else
                    proposeType2 = 5
                endif
            else
                proposeType2 = 4
            endif
        end select

        return

    end function

    integer function proposeType(propose,set,iter)
    
        implicit none
        real(kind=8), intent(in) :: propose
        type(T_MCMC_SET), intent(in) :: set
        integer, intent(in) :: iter

        real(kind=ii10) random
        call random_number(random)

        proposeType = 4
        if(set%sigdep /= 0)then
            if (propose < 0.2) then
                proposeType = 1
            elseif(propose < 0.4) then
                proposeType = 2
            elseif(propose < 0.6) then
                proposeType = 3
            elseif(propose < 0.95) then
                proposeType = 4
            else
                proposeType = proposeSigmaType(set%datatype)
            endif
        else
            if (propose < 0.2) then
                proposeType = 1
            elseif(propose < 0.4) then
                proposeType = 2
            elseif(propose < 0.6) then
                proposeType = 3
            else
                proposeType = 4
            endif
        endif

        return

    end function

    integer function proposeSigmaType(datatype)
        implicit none
        integer, intent(in) :: datatype

        real(kind=ii10) random
        call random_number(random)

        select case(datatype)
        case(0,1)
            proposeSigmaType = 5
        case(2)
            proposeSigmaType = 5
        case(3)
            if(random<0.5)then
                proposeSigmaType = 5
            else
                proposeSigmaType = 6
            endif
        endselect

    endfunction

    subroutine cell_birth(RTI,mcmc_set,point,pm,prob,lerr)
        implicit none
        type(T_MCMC_SET), intent(in)            :: mcmc_set
        real(kind=ii10),dimension(2), intent(out):: point
        real(kind=ii10),dimension(3), intent(out):: pm
        real(kind=ii10), intent(out)            :: prob
        logical, intent(out)                    :: lerr
        type(T_RUN_INFO), intent(inout)         :: RTI

        ! local variable
        real(c_double) :: x, y, z
        real(c_double) :: vp, vs, rho
        integer(c_size_t)  :: ncells
        integer idx
        ! debug    

        ncells = RTI%ncells
        if( ncells == mcmc_set%ncell_max)then
            lerr = .false.
            prob = 0
            return
        endif

        x = mcmc_set%grid%xmin + unirand(RTI%randcount)*(mcmc_set%grid%xmax-mcmc_set%grid%xmin)
        y = mcmc_set%grid%ymin + unirand(RTI%randcount)*(mcmc_set%grid%ymax-mcmc_set%grid%ymin)
        ! gaussian
        if(mcmc_set%kernel == 0) then
            call kdtree_locate(RTI%points(:,1:ncells),[x,y],idx)
            pm = RTI%parameters(:,idx)
            vs = RTI%parameters(2,idx) + gasdev(RTI%randcount)*mcmc_set%sigma_vs2
            vp = RTI%parameters(1,idx) + gasdev(RTI%randcount)*mcmc_set%sigma_vp2
            !rho = pm%rho + gasdev(RTI%randcount)*mcmc_set%sigma_rho
            !vp = vs2vp(vs)
            rho = vp2rho(vp)
            if(mcmc_set%datatype == 0)then
                prob = log(mcmc_set%sigma_vp2*sqrt(2*PII)) + (vp-pm(1))**2/(2*mcmc_set%sigma_vp2**2)
            else
                prob = log(mcmc_set%sigma_vp2*sqrt(2*PII)) + (vp-pm(1))**2/(2*mcmc_set%sigma_vp2**2)&
                     + log(mcmc_set%sigma_vs2*sqrt(2*PII)) + (vs-pm(2))**2/(2*mcmc_set%sigma_vs2**2)
            endif

        else
            ! prior
            vs = mcmc_set%vsmin + unirand(RTI%randcount)*(mcmc_set%vsmax-mcmc_set%vsmin)
            vp = mcmc_set%vpmin + unirand(RTI%randcount)*(mcmc_set%vpmax-mcmc_set%vpmin)
            !rho = mcmc_set%rhomin + unirand(RTI%randcount)*(mcmc_set%rhomax-mcmc_set%rhomin)
            rho = vp2rho(vp)
            prob = 0
        endif

        if(mcmc_set%datatype==2)then
            vp = vs2vp(vs)
        endif
        rho = vp2rho(vp)

        if( vp<mcmc_set%vpmin .or. vp>mcmc_set%vpmax .or. &
            vs<mcmc_set%vsmin .or. vs>mcmc_set%vsmax) then
            lerr = .false.
            prob = 0
            return
        endif

        ! update the points and parameters array
        !call cgal_get_vertices(delaunay_ptr,mcmc_set%ncell_max,RTI%points,RTI%parameters)
        RTI%points(1,ncells+1) = x
        RTI%points(2,ncells+1) = y

        RTI%parameters(1,ncells+1) = vp
        RTI%parameters(2,ncells+1) = vs
        RTI%parameters(3,ncells+1) = rho
        pm = RTI%parameters(:,ncells+1)
        point = RTI%points(:,ncells+1)

        RTI%ncells = ncells+1
        lerr = .true.

        return

    end subroutine cell_birth

    subroutine cell_death(RTI,mcmc_set,iremove,pm,prob,lerr)

        implicit none

        type(T_MCMC_SET), intent(in)            :: mcmc_set
        integer(c_int), intent(out)             :: iremove
        real(kind=ii10),dimension(3), intent(out):: pm
        real(kind=ii10), intent(out)            :: prob
        logical, intent(out)                    :: lerr
        type(T_RUN_INFO), intent(inout)         :: RTI

        ! local variable
        real(c_double), dimension(2,mcmc_set%ncell_max):: points_copy
        real(c_double), dimension(3,mcmc_set%ncell_max):: parameters_copy
        real(c_double), dimension(2):: pt
        real(c_double), dimension(3):: pm1
        integer(c_size_t)   :: ncells
        integer idx
        !debug
        !integer i

        ncells = RTI%ncells
        if( ncells == mcmc_set%ncell_min)then
            lerr = .false.
            prob = 0
            return
        endif

        iremove = ceiling(unirand(RTI%randcount)*ncells)
        if( iremove==0 ) iremove = 1

        pt = RTI%points(:,iremove)
        ! update the points and parameters array
        points_copy = RTI%points
        parameters_copy = RTI%parameters
        RTI%points(:,iremove:ncells-1) = points_copy(:,iremove+1:ncells)
        RTI%parameters(:,iremove:ncells-1) = parameters_copy(:,iremove+1:ncells)
        RTI%ncells = ncells-1
        pm = RTI%parameters(:,iremove)

        if(mcmc_set%kernel == 0) then !gaussian kernel
            call kdtree_locate(RTI%points(:,1:RTI%ncells),pt,idx)
            pm1 = RTI%parameters(:,idx)
            if(mcmc_set%datatype == 0)then
                prob = log( 1/(mcmc_set%sigma_vp2*sqrt(2*PII)) ) - &
                    (pm1(1)-pm(1))**2/(2*mcmc_set%sigma_vp2**2)
            else
                prob = log( 1/(mcmc_set%sigma_vp2*sqrt(2*PII)) ) - &
                    (pm1(1)-pm(1))**2/(2*mcmc_set%sigma_vp2**2) + &
                    log( 1/(mcmc_set%sigma_vs2*sqrt(2*PII)) )    - &
                    (pm1(2)-pm(2))**2/(2*mcmc_set%sigma_vs2**2)
            endif
        else ! prior kernel
            prob = 0
        endif
        lerr = .true.


        return

    end subroutine cell_death
    
    subroutine cell_move(RTI,mcmc_set,imove,moved_dst,lerr)
        implicit none

        type(T_MCMC_SET), intent(in)            :: mcmc_set
        integer(c_int), intent(out)             :: imove
        real(kind=ii10),dimension(2), intent(out):: moved_dst
        logical, intent(out)                    :: lerr
        type(T_RUN_INFO), intent(inout)         :: RTI

        ! local variable
        real(kind=ii10),dimension(2) ::pt_dst
        integer(c_int) :: verbose
        real(kind=8) :: propose
        ! debug
        
        imove = ceiling(unirand(RTI%randcount)*RTI%ncells)
        if( imove==0 ) imove = 1
        
        ! choose one coordinate to perturb the position
        ! initiate first
        pt_dst = RTI%points(:,imove)
        !pt_dst%x = RTI%points(1,imove) + gasdev(RTI%randcount)*mcmc_set%pd*(mcmc_set%grid%xmax-mcmc_set%grid%xmin)/100
        !pt_dst%y = RTI%points(2,imove) + gasdev(RTI%randcount)*mcmc_set%pd*(mcmc_set%grid%ymax-mcmc_set%grid%ymin)/100
        !pt_dst%z = RTI%points(3,imove) + gasdev(RTI%randcount)*mcmc_set%pd*(mcmc_set%grid%zmax-mcmc_set%grid%zmin)/100
        propose = unirand(RTI%randcount)
        if (propose < 0.5) then
            pt_dst(1) = RTI%points(1,imove) + gasdev(RTI%randcount)*mcmc_set%pd*(mcmc_set%grid%xmax-mcmc_set%grid%xmin)/100
            if(mcmc_set%datatype==2.and. RTI%points(2,imove)>(mcmc_set%grid%ymin+mcmc_set%grid%ymax)/2)&
                pt_dst(1) = RTI%points(1,imove) + gasdev(RTI%randcount)*mcmc_set%pd2*(mcmc_set%grid%xmax-mcmc_set%grid%xmin)/100
        else
            pt_dst(2) = RTI%points(2,imove) + gasdev(RTI%randcount)*mcmc_set%pd*(mcmc_set%grid%ymax-mcmc_set%grid%ymin)/100
            if(mcmc_set%datatype==2.and. RTI%points(2,imove)>(mcmc_set%grid%ymin+mcmc_set%grid%ymax)/2)&
                pt_dst(2) = RTI%points(2,imove) + gasdev(RTI%randcount)*mcmc_set%pd2*(mcmc_set%grid%ymax-mcmc_set%grid%ymin)/100
        endif

        if( pt_dst(1)<mcmc_set%grid%xmin .or. pt_dst(1)>mcmc_set%grid%xmax &
            .or. pt_dst(2)<mcmc_set%grid%ymin .or. pt_dst(2)>mcmc_set%grid%ymax)then
            lerr = .false.
            return
        endif

        ! update the points and parameters array
        !call cgal_get_vertices(delaunay_ptr,mcmc_set%ncell_max,RTI%points,RTI%parameters)
        RTI%points(:,imove) = pt_dst
        moved_dst = pt_dst
        lerr = .true.

        ! debug
        !if(any(RTI%points(:,1:RTI%ncells) /= points(:,1:RTI%ncells)))then
        !    call exception_raiseError('Local vertices not equal to the one&
        !        & from delaunay in move step!')
        !endif
        return

    end subroutine cell_move

    subroutine vel_change(RTI,mcmc_set,ivalue,pm_src,pm,lerr)
        implicit none

        type(T_MCMC_SET), intent(in)            :: mcmc_set
        real(kind=ii10),dimension(:), intent(out):: pm, pm_src
        integer(c_int), intent(out)             :: ivalue
        logical, intent(out)                    :: lerr
        type(T_RUN_INFO), intent(inout)         :: RTI

        ! local
        real(kind=ii10),dimension(2):: pt_src
        real(kind=ii10),dimension(3):: pm_dst
        integer(c_int) :: verbose
        ! debug

        ivalue = ceiling(unirand(RTI%randcount)*RTI%ncells)
        if( ivalue==0 ) ivalue = 1
        
        pt_src = RTI%points(:,ivalue)

        pm_dst(2) = RTI%parameters(2,ivalue) + gasdev(RTI%randcount)*mcmc_set%sigma_vs
        pm_dst(1) = RTI%parameters(1,ivalue) + gasdev(RTI%randcount)*mcmc_set%sigma_vp
        if(mcmc_set%datatype==2)then
            if(pt_src(2)>(mcmc_set%grid%ymin+mcmc_set%grid%ymax)/2) &
                 pm_dst(2) = RTI%parameters(2,ivalue) + gasdev(RTI%randcount)*mcmc_set%sigma_vs2
            pm_dst(1) = vs2vp(pm_dst(2))
        endif
        pm_dst(3) = vp2rho(pm_dst(1))

        !dev%vp = gasdev(RTI%randcount)*mcmc_set%sigma_vp
        !dev%rho = gasdev(RTI%randcount)*mcmc_set%sigma_rho

        ! out of prior or vp<1.1*vs (vp is generally bigger than vs)
        if( pm_dst(1)<mcmc_set%vpmin .or. pm_dst(1)>mcmc_set%vpmax .or. &
            pm_dst(2)<mcmc_set%vsmin .or. pm_dst(2)>mcmc_set%vsmax) then
            lerr = .false.
            return
        endif

        ! update the points and parameters array
        !call cgal_get_vertices(delaunay_ptr,mcmc_set%ncell_max,points,parameters)
        pm_src = RTI%parameters(:,ivalue)
        RTI%parameters(:,ivalue) = pm_dst
        pm = pm_dst
        lerr = .true.

        return

    end subroutine vel_change

    subroutine ssigma_change(RTI,mcmc_set,iperiod,lerr)
        implicit none
        type(T_RUN_INFO), intent(inout) :: RTI
        type(T_MCMC_SET), intent(in) :: mcmc_set
        integer(c_int), intent(out) :: iperiod
        logical :: lerr

        ! local
        real( kind=ii10 ) random

        iperiod = ceiling(unirand(RTI%randcount)*size(RTI%snoise0))
        if( iperiod == 0 ) iperiod = 1

        random = unirand(RTI%randcount)
        if( random < 0.5 .or. mcmc_set%datatype==2) then
            RTI%snoise0(iperiod) = RTI%snoise0(iperiod) + &
            GASDEV(RTI%randcount)*RTI%sigma_n0(iperiod)
        else
            RTI%snoise1(iperiod) = RTI%snoise1(iperiod) + &
            GASDEV(RTI%randcount)*RTI%sigma_n1(iperiod)
        endif

        if( RTI%snoise0(iperiod) < RTI%n0_min(iperiod) .or. RTI%snoise0(iperiod) > RTI%n0_max(iperiod) .or.&
            RTI%snoise1(iperiod) < RTI%n1_min(iperiod) .or. RTI%snoise1(iperiod) > RTI%n1_max(iperiod)) then
            lerr = .false.
            return
        else
            lerr = .true.
        endif

        return

    end subroutine

    subroutine stat_rti(RTI,model,mcmc_set)
        implicit none
        type(T_RUN_INFO), intent(inout) :: RTI
        type(T_MOD), intent(in) :: model
        type(T_MCMC_SET), intent(in) :: mcmc_set

        RTI%aveS = RTI%aveS + model%vs
        RTI%stdS = RTI%stdS + model%vs**2
        RTI%aveP = RTI%aveP + model%vp
        RTI%stdP = RTI%stdP + model%vp**2

        RTI%nthin = RTI%nthin + 1

    endsubroutine

    subroutine kdtree_to_grid(RTI, grid, model)
        use kdtree2_precision_module, only : kdkind
        use kdtree2_module, only : kdtree2, kdtree2_result, kdtree2_create,&
            kdtree2_destroy, kdtree2_n_nearest
        use omp_lib
        implicit none
        type(T_RUN_INFO), intent(inout) :: RTI
        type(T_GRID), intent(in) :: grid
        type(T_MOD), intent(inout) :: model

        ! kd-tree
        type(kdtree2), pointer :: tree
        type(kdtree2_result), dimension(1) :: results
        real(kdkind), dimension(2) :: qv

        ! grid
        integer        :: ix0, ix1, iy0, iy1, iz0, iz1

        integer i,j,k

        !real(kind=8) :: t

        ! create kd-tree
        tree => kdtree2_create(RTI%points(:,1:RTI%ncells), sort=.false., rearrange=.true.)

        ! convert to grid based model using kd-tree nearest search

        ! the lowest and highest indices for grid model vp, vs and rho
        ix0 = 1 
        ix1 = grid%nx
        iy0 = 1
        iy1 = grid%ny

        ! call nearest search of kd-tree
        !t =  omp_get_wtime()
        !$omp parallel
        !$omp do private(qv,i,j,results)
        do i = ix0, ix1
            do j = iy0, iy1
                    qv = [grid%xmin+(i-1)*grid%dx, grid%ymin+(j-1)*grid%dy]
                    call kdtree2_n_nearest(tp=tree,qv=qv,nn=1,results=results)
                    RTI%sites_id(j,i) = results(1)%idx
                    model%vp(j,i) = RTI%parameters(1,results(1)%idx)
                    model%vs(j,i) = RTI%parameters(2,results(1)%idx)
                    model%rho(j,i) = RTI%parameters(3,results(1)%idx)
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        !t =  omp_get_wtime() -  t
        !call log_msg('parallelized kdtree code: '//rtoa(t) )

        call kdtree2_destroy(tree)

        return

    end subroutine


    subroutine kdtree_locate(points, pt, idx)
        use kdtree2_precision_module, only : kdkind
        use kdtree2_module, only : kdtree2, kdtree2_result, kdtree2_create,&
            kdtree2_destroy, kdtree2_n_nearest
        use omp_lib
        implicit none
        real(kdkind), dimension(:,:), intent(in) :: points
        real(kdkind), dimension(:), intent(in) :: pt
        integer, intent(out) :: idx

        ! kd-tree
        type(kdtree2), pointer :: tree
        type(kdtree2_result), dimension(1) :: results

        ! create kd-tree
        tree => kdtree2_create(points, sort=.false., rearrange=.true.)

        call kdtree2_n_nearest(tp=tree,qv=pt,nn=1,results=results)
        idx = results(1)%idx

        call kdtree2_destroy(tree)

        return

    end subroutine

    !
    ! tempering related subroutines
    !

    subroutine choose_t1(RTI,t1_chain)
    ! the first temperature is chosen randomly
    implicit none
    type(T_RUN_INFO), intent(inout) :: RTI
    integer, intent(out)    :: t1_chain
    integer         :: ii, t1_index
    real(kind=ii10) :: t1_value
    
    t1_chain=0
    do while (.not.(any((/(ii,ii=1,RTI%number_of_chains)/).eq.t1_chain)))
            !t1_chain=nint(grnd()*(RTI%number_of_chains+1))
            t1_chain=ceiling(unirand(RTI%randcount)*RTI%number_of_chains)
            if(t1_chain==0) t1_chain = 1
    end do
    t1_value=RTI%temperature_values(t1_chain)
    t1_index=RTI%temperature_indices(t1_chain)
    
    if (debug_mode) call log_msg('----- t1: chain='//trim(itoa(t1_chain))//&
                              ' index='//trim(itoa(t1_index))//&
                              ', value='//trim(rtoa(t1_value))//' -----')
    
    end subroutine
    
    subroutine choose_t2(RTI,mcmc_set,t1_chain,t2_chain)
    ! choose the second temperature depending on the choice made in mksamples.in
    implicit none
    type(T_RUN_INFO), intent(inout) :: RTI
    type(T_MCMC_SET), intent(in) :: mcmc_set
    integer, intent(in)    :: t1_chain
    integer, intent(out)    :: t2_chain

    integer         :: t1_index, t2_index
    real(kind=ii10) :: t1_value, t2_value
    integer         :: ii
    
    t1_value=RTI%temperature_values(t1_chain)
    t1_index=RTI%temperature_indices(t1_chain)
    t2_chain=0
    if (mcmc_set%jump_type.eq.0) then ! choose a random temperature different from t1
       do while (.not.(any((/(ii,ii=1,RTI%number_of_chains)/).eq.t2_chain)).or.abs(t1_value-t2_value)<eps)
           !t2_chain=nint(grnd()*(RTI%number_of_chains+1))
           t2_chain=ceiling(unirand(RTI%randcount)*RTI%number_of_chains)
           if(t2_chain==0) t2_chain = 1
           t2_value=RTI%temperature_values(t2_chain)
           t2_index=RTI%temperature_indices(t2_chain)
       end do
    else if (mcmc_set%jump_type.eq.1) then 
        ! choose a temperature randomly between the two neighbours of t1
        if (abs(t1_value-1.0)<eps) then ! if t1=1, choose the next temperature above it
            t2_value=RTI%temperature_values0(mcmc_set%number_of_1s+1)
            t2_index=RTI%temperature_indices0(mcmc_set%number_of_1s+1)
        else if (abs(t1_value-RTI%temperature_values0(RTI%number_of_temperatures))<eps) then 
            ! if t1 is the maximum temperature, choose the first below it
            t2_value=RTI%temperature_values0(RTI%number_of_temperatures-1)
            t2_index=RTI%temperature_indices0(RTI%number_of_temperatures-1)
        else ! in all other cases choose randomly between the two neighbours of t1
            if (unirand(RTI%randcount).lt.0.5) then
                ii=-1 ! one below
            else
                ii=1 ! one above
            end if
            t2_value=RTI%temperature_values0(t1_index+ii)
            t2_index=RTI%temperature_indices0(t1_index+ii)
        end if
        t2_chain=minloc(abs(RTI%temperature_indices-t2_index),1)
    else if (mcmc_set%jump_type.eq.2) then ! choose the nearest temperature to t1
        if (abs(t1_value-1.0)<eps) then ! if t1=1, choose the next temperature above it
            t2_value=RTI%temperature_values0(mcmc_set%number_of_1s+1)
            t2_index=RTI%temperature_indices0(mcmc_set%number_of_1s+1)
        else if (abs(t1_value-RTI%temperature_values0(RTI%number_of_temperatures))<eps) then 
            ! if t1 is the maximum temperature, choose the first below it
            t2_value=RTI%temperature_values0(RTI%number_of_temperatures-1)
            t2_index=RTI%temperature_indices0(RTI%number_of_temperatures-1)
        else ! in all other cases choose the nearest temperature between the two neighbours of t1
            ii=minloc(abs(RTI%temperature_values0-t1_value),1,abs(RTI%temperature_values0-t1_value)>eps)
            if(ii==1)then
                ii = ceiling(unirand(RTI%randcount)*mcmc_set%number_of_1s)
                if(ii==0) ii=1
            endif
            t2_value=RTI%temperature_values0(ii)
            t2_index=RTI%temperature_indices0(ii)
        end if
        t2_chain=minloc(abs(RTI%temperature_indices-t2_index),1)
    end if
    
    if (debug_mode) call log_msg('----- t2: chain='//trim(itoa(t2_chain))//&
                              ' index='//trim(itoa(t2_index))//&
                              ', value='//trim(rtoa(t2_value))//' -----')
    
    end subroutine
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine&
        do_tempering(RTI,t1_chain,t2_chain,likelihood_t1,likelihood_t2,astatus)
    !TODO: arguments list are too long, using derived type
    
    implicit none
    type(T_RUN_INFO), intent(inout) :: RTI
    integer, intent(in)         :: t1_chain, t2_chain
    real(kind=ii10), intent(in) :: likelihood_t1, likelihood_t2
    integer, intent(out) :: astatus

    real(kind=ii10)        :: alpha_swap
    real(kind=ii10)        :: random_deviate
    real(kind=ii10)        :: t1_value, t2_value
    logical  swap_temperatures
    
    astatus = -1
    swap_temperatures = .false.
    
    t1_value=RTI%temperature_values(t1_chain)
    t2_value=RTI%temperature_values(t2_chain)

    alpha_swap=minval([dble(0),-likelihood_t1/t2_value&
                               +likelihood_t2/t2_value&
                               -likelihood_t2/t1_value&
                               +likelihood_t1/t1_value])
    
    random_deviate=unirand(RTI%randcount)
    if (log(random_deviate).le.alpha_swap)then
        RTI%temperature_values((/t1_chain,t2_chain/))=RTI%temperature_values((/t2_chain,t1_chain/))
        RTI%temperature_indices((/t1_chain,t2_chain/))=RTI%temperature_indices((/t2_chain,t1_chain/))
        RTI%accepted_tempering = RTI%accepted_tempering + 1
        swap_temperatures = .true.
        astatus = 0
    end if
    RTI%total_tempering = RTI%total_tempering + 1
    
    if (debug_mode) then
        if (swap_temperatures) call log_msg('----- swapping: yes ------')
        if (.not.swap_temperatures) call log_msg('----- swapping: no------')
    end if
    
    end subroutine

    ! importance sampling, i.e. weight average through all temperatures
    ! generally, a solution of reverse logistic regression would be better. But
    ! it is more obvious to implement after sampling
    function ImportanceSample(temperature,loglike,norm) result(w)
        implicit none
        real(kind=ii10), intent(in) :: temperature
        real(kind=ii10), intent(in) :: loglike
        real(kind=ii10), intent(in) :: norm
        real(kind=16) :: w

        w = exp(loglike*(1.0-1.0/temperature)-norm*(1.0-1.0/temperature))

        !print *, w
        if(w /= w)then
            call log_msg('loglike: '//rtoa(loglike))
            call log_msg('temperature: '//rtoa(temperature))
            call exception_raiseError('Error when importance sampling')
        endif

    end function
    
end module m_mcmc
