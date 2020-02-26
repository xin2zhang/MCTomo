!
! modules m_mcmc
! reversible jump Markov Chain Monte Carlo algorithm
! 
module m_hmcmc

    use m_exception, only       : exception_raiseError
    use m_logger, only          : log_msg
    use m_utils, only           : ii10, itoa, rtoa, delete_file, vs2vp, vp2rho, results_dir, last_dir, resume_dir, FILEDIR_SEPARATOR
    use m_settings, only        : T_MCMC_SET, T_GRID, T_MOD, mod_setup
    use like_settings, only     : T_LIKE_SET, T_DATA, write_sources
    use m_likelihood, only      : likelihood, T_LIKE, like_setup, noise_likelihood
    use run_info, only          : T_RUN_INFO, T_SAMPLE, write_info, write_samples
    use read_write, only        : write_mean, write_var, write_vertices, write_data
    use mt19937, only           : unirand

    use m_hmc, only             : hmc_step, T_HMC_SET


    ! use all the available procedure, types and variable in the following
    ! module
    use m_mcmc
    use cgal_delaunay
    use iso_c_binding
#ifdef MPI
    use mpi
#endif

    implicit none
    private

    public :: hmcmc
    public :: T_HMC_SAMPLES
    public :: write_likelihood_hmc, write_number_of_cells_hmc
    public :: create_hmcsamples, read_hmcsamples, write_hmcsamples
   
    ! static value
    real(kind=ii10), parameter :: one = 1.0
    real(kind=ii10), parameter :: PII = 3.1415926
    real(kind=ii10), parameter :: VPVS = 1.732
    real(kind=ii10), parameter :: eps = 1e-8
    integer, parameter :: maxTryNum = 1000

    ! debug
    real t1, t2
    logical debug_mode

    ! for hmc, usually the model is changed everywhere
    ! need to save the whole model every time (large space!)
    type T_HMC_SAMPLES
        integer(c_int), dimension(:), allocatable     :: ncells
        real(c_double), dimension(:,:,:), allocatable :: points
        real(c_double), dimension(:,:,:), allocatable :: parameters
        real(c_double), dimension(:,:), allocatable   :: bnoise0, bnoise1
        real(c_double), dimension(:,:), allocatable   :: snoise0, snoise1
        real(c_double), dimension(:,:,:), allocatable :: locations
        real(c_double), dimension(:), allocatable     :: negative_lglike
        real(c_double), dimension(:), allocatable     :: weighted_misfits
        real(c_double), dimension(:), allocatable     :: unweighted_misfits
    endtype T_HMC_SAMPLES


contains

    subroutine setup_hmc(mcmc_set, nsrc, hmc_set)
        implicit none
        type(T_MCMC_SET), intent(in) :: mcmc_set
        integer, intent(in) :: nsrc
        type(T_HMC_SET), intent(out) :: hmc_set

        integer ndims

        hmc_set%n_leapfrogs = mcmc_set%n_leapfrogs
        hmc_set%priortype = mcmc_set%priortype

        ndims = mcmc_set%ncell_max*2
        if(mcmc_set%locate==1)then
            ndims = ndims + 4*nsrc
        endif

        allocate(hmc_set%stepsize(ndims))
        allocate(hmc_set%diagonal_matrix(ndims))
        allocate(hmc_set%pm1(ndims))
        allocate(hmc_set%pm2(ndims))

        hmc_set%stepsize = 0.01
        hmc_set%diagonal_matrix = 1.0
        hmc_set%pm1 = -huge(0.0_ii10)
        hmc_set%pm2 = huge(0.0_ii10)
    endsubroutine setup_hmc

    subroutine create_hmcsamples(hmcsamples,nsamples,dat,mcmc_set)
        implicit none
        type(T_HMC_SAMPLES), intent(out) :: hmcsamples
        integer, intent(in) :: nsamples
        type(T_DATA), dimension(:), intent(in) :: dat
        type(T_MCMC_SET), intent(in) :: mcmc_set

        integer nd, np

        !nsamples = mcmc_set%nsamples

        allocate(hmcsamples%ncells(nsamples))
        hmcsamples%ncells = mcmc_set%ncell_min
        allocate(hmcsamples%negative_lglike(nsamples))
        hmcsamples%negative_lglike = 0.0
        allocate(hmcsamples%weighted_misfits(nsamples))
        hmcsamples%weighted_misfits = 0.0
        allocate(hmcsamples%unweighted_misfits(nsamples))
        hmcsamples%unweighted_misfits = 0.0
        allocate(hmcsamples%points(3,mcmc_set%ncell_max,nsamples))
        hmcsamples%points = 0.0
        allocate(hmcsamples%parameters(3,mcmc_set%ncell_max,nsamples))
        hmcsamples%parameters = 0.0

        select case (mcmc_set%datatype)
        case (0,1)
            ! allocate one element for surface wave, just for safe
            nd = dat(1)%np
            np = 1
        case (2)
            ! allocate one element for body wave
            nd = 1
            np = dat(1)%np
        case (3)
            ! for body waves
            nd = dat(1)%np
            np = dat(2)%np
        end select

        allocate(hmcsamples%bnoise0(nd,nsamples))
        allocate(hmcsamples%bnoise1(nd,nsamples))
        allocate(hmcsamples%snoise0(np,nsamples))
        allocate(hmcsamples%snoise1(np,nsamples))
        hmcsamples%bnoise0 = 0.0
        hmcsamples%bnoise1 = 0.0
        hmcsamples%snoise0 = 0.0
        hmcsamples%snoise1 = 0.0

        if(mcmc_set%locate==1)then
            allocate(hmcsamples%locations(4,dat(1)%nsrc,nsamples))
            hmcsamples%locations = 0.0
        endif

    endsubroutine create_hmcsamples

    integer function proposeType(propose,set,iter)
    
        implicit none
        real(kind=8), intent(in) :: propose
        type(T_MCMC_SET), intent(in) :: set
        integer, intent(in) :: iter

        real(kind=ii10) random
        call random_number(random)

        proposeType = 4
        if (propose < 0.3)then
            proposeType = 1
        elseif(propose < 0.6)then
            proposeType = 2
        elseif(propose < 0.9)then
            proposeType = 3
        elseif(propose < 0.95)then
            proposeType = 4
        else
            proposeType = proposeSigmaType(set%datatype)
        endif
        !select case (mod(iter,2))
        !case ( 0 )
        !    if (propose < 0.333) then
        !        proposeType = 1
        !    elseif(propose < 0.666) then
        !        proposeType = 2
        !    else
        !        proposeType = 3
        !    endif
        !case ( 1 )
        !    if(set%sigdep /= 0) then
        !        if (propose < 0.500)then
        !            proposeType = 4
        !        else
        !            proposeType = proposeSigmaType(set%datatype)
        !        endif
        !    else
        !        proposeType = 4
        !    endif
        !end select

        return

    end function

    !
    ! subroutine mcmc
    ! implement the reversible jump Markov Chain Monte Carlo mtheod with hmc
    ! input  :
    !   data
    !   RTI: run time information
    !   likelihood settings
    !   settings
    ! output :
    !   samples
    !
    subroutine hmcmc(hmcsamples,dat,RTI,mcmc_set,like_set)

        implicit none

        type(T_HMC_SAMPLES), intent(out)                :: hmcsamples
        type(T_DATA), dimension(:), intent(in)          :: dat
        type(T_RUN_INFO), intent(inout)                 :: RTI
        type(T_MCMC_SET), intent(in)                    :: mcmc_set
        type(T_LIKE_SET), intent(in)                    :: like_set

        ! local variables
        ! voronoi tessellation
        type(c_ptr)                                     :: delaunay_ptr
        type(d3), dimension(:), allocatable, target     :: points
        type(p3), dimension(:), allocatable, target     :: parameters
        integer(c_size_t)                               :: ncells_copy
        real(ii10), dimension(3,mcmc_set%ncell_max)     :: points_copy, parameters_copy
        real(c_double), dimension(:), allocatable       :: bnoise0_copy, bnoise1_copy
        real(c_double), dimension(:), allocatable       :: snoise0_copy, snoise1_copy
        real(c_double), dimension(:,:), allocatable     :: sources_copy
        type(c_ptr)                                     :: delaunay_ptr_copy
        type(T_GRID)                                    :: grid
        type(d3), dimension(2)                          :: bnd_box
        type(T_MOD)                                     :: model, model_copy
        integer, dimension(:,:,:), allocatable          :: sites_id_copy

        ! hmc step
        type(T_HMC_SET)   :: hmc_set
        real(kind=ii10), dimension(2*mcmc_set%ncell_max+dat(1)%nsrc*4) :: qvals_in
        integer           :: ndim
        ! propose related variable
        type( T_LIKE )    :: like, like_copy
        real( kind=8 )    :: propose
        integer           :: ptype, icycle
        real( kind=ii10 ) :: prob
        real( kind=ii10 ) :: alpha, random
        integer(c_int)    :: iremove, imove, iperiod
        logical           :: accepted, lerr
        type(d3)          :: point
        type(p3)          :: pm
        integer           :: iter, countFM, nthin

        integer i, j, k
        ! debug
        !character(100) :: debug_info

        ! if zero sample size, return
        if(RTI%nsamples==0) return

        ! delelte run-time sample file
        call delete_file(trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_sample_'//itoa(mcmc_set%processor)//'.dat')
        ! initialize and allocate several variables
        call create_hmcsamples(hmcsamples,RTI%nsamples, dat, mcmc_set)
        call setup_hmc(mcmc_set,dat(1)%nsrc, hmc_set)
        grid = mcmc_set%grid
        call mod_setup(model, grid)
        call mod_setup(model_copy, grid)

        ! initialize points and parameters corresponding to the current delaunay

        ! calculate likelihood of the initial sample
        bnd_box(1) = d3(grid%xmin,grid%ymin,grid%zmin)
        bnd_box(2) = d3(grid%xmax,grid%ymax,grid%zmax)
        call kdtree_to_grid(RTI, grid, bnd_box, model)

        call like_setup(like, dat, like_set, mcmc_set%ncell_max)
        call like_setup(like_copy, dat, like_set, mcmc_set%ncell_max)
        call likelihood(dat,model,RTI,bnd_box,like_set,like)

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

        ! > create cgal delaunay triangulation first
        allocate( points(mcmc_set%ncell_max) )
        allocate( parameters(mcmc_set%ncell_max) )
        call array_points(points,RTI%points)
        call array_parameters(parameters,RTI%parameters)
        call cgal_delaunay_build(points, parameters, RTI%ncells, delaunay_ptr )

        ! > initialize some setting values
        iter = 0
        icycle = 0
        delaunay_ptr_copy = C_NULL_PTR
        !allocate( samples(mcmc_set%nsamples) )
        !call init_sample(samples) ! initialize sample
        do while(iter < RTI%nsamples)

            iter = iter + 1
            RTI%sampletotal = RTI%sampletotal + 1

            ! back up some values, like delaunay triangualtion, vp, vs, rho,
            ! vel, likelihood, these will possibley be recovered if new sample
            ! not accepteded
            call cgal_delaunay_copy(delaunay_ptr,delaunay_ptr_copy)
            ncells_copy = RTI%ncells
            points_copy = RTI%points
            parameters_copy = RTI%parameters
            bnoise0_copy = RTI%bnoise0
            bnoise1_copy = RTI%bnoise1
            snoise0_copy = RTI%snoise0
            snoise1_copy = RTI%snoise1
            sources_copy = RTI%locations
            sites_id_copy = RTI%sites_id
            model_copy = model
            like_copy = like

            ! initial some values, eg. logical accepted, propose type, random number
            accepted = .false.
            if(icycle ==0 .or. countFM>=maxTryNum)then
                random = unirand(RTI%randcount)
                propose = unirand(RTI%randcount)
                ptype = proposeType(propose,mcmc_set,iter)
                countFM = 0
            endif

            icycle = 0
            select case ( ptype )
            case(1) 
                ! count
                RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
                !samples(iter)%step = ptype

                call cell_birth(delaunay_ptr,RTI,mcmc_set,&
                                point,pm,bnd_box,prob,lerr)
                if(.not.lerr) goto 100
                call kdtree_to_grid(RTI,grid,bnd_box,model)
                !call log_msg('birth')
                !call check_grid(RTI,grid,model)

                ! calculate travel time and likelihood value
                !call cpu_time(t1)
                call likelihood(dat,model,RTI,bnd_box,like_set,like)
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
                        log(mcmc_set%vsmax-mcmc_set%vsmin)-like%like+like_copy%like])
                else
                    alpha = minval([log(one), -like%like+like_copy%like])
                endif
                if(log(random)<alpha)then
                    accepted = .true.
                    ! update accepted ratio and samples
                    !samples(iter)%step = ptype
                    !samples(iter)%accepted = accepted
                    !samples(iter)%vindex = 0
                    !samples(iter)%ncells = delaunay_size(delaunay_ptr)
                    !samples(iter)%like = like%like
                    !samples(iter)%misfit = like%misfit
                    !samples(iter)%unweighted_misfit = like%unweighted_misfit
                    !samples(iter)%noise0 = 0
                    !samples(iter)%noise1 = 0
                    !samples(iter)%coord = [point%x,point%y,point%z]
                    !samples(iter)%values = [pm%vp,pm%vs,pm%rho]
                    RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                endif
                ! debug

            case(2)
                ! count
                RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
                !samples(iter)%step = ptype
                ! remove a cell
                call cell_death(delaunay_ptr,RTI,mcmc_set,iremove,pm,bnd_box,prob,lerr)
                if(.not.lerr) goto 100
                ! update sites id, TODO: using a hash table might be faster
                where(RTI%sites_id > iremove)
                    RTI%sites_id = RTI%sites_id - 1
                endwhere
                call kdtree_to_grid(RTI,grid,bnd_box,model)
                !call log_msg('death: '//itoa(iremove)//' vp:'//rtoa(pm%vp))
                !write(*,*) 'vp: ', RTI%parameters(1,1:RTI%ncells)
                !call check_grid(RTI,grid,model)

                ! calculate travel time and likelihood value
                call likelihood(dat,model,RTI,bnd_box,like_set,like)

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
                        log(mcmc_set%vsmax-mcmc_set%vsmin)-like%like+like_copy%like])
                else
                    alpha = minval([log(one), -like%like+like_copy%like])
                endif
                if(log(random)<alpha)then
                    accepted = .true.
                    ! update accepted ratio and samples
                    !samples(iter)%step = ptype
                    !samples(iter)%accepted = accepted
                    !samples(iter)%vindex = iremove
                    !samples(iter)%ncells = delaunay_size(delaunay_ptr)
                    !samples(iter)%like = like%like
                    !samples(iter)%misfit = like%misfit
                    !samples(iter)%unweighted_misfit = like%unweighted_misfit
                    !samples(iter)%noise0 = 0
                    !samples(iter)%noise1 = 0
                    !samples(iter)%coord = 0
                    !samples(iter)%values = 0
                    RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                endif
                ! debug

            case(3)
                ! count
                RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
                !samples(iter)%step = ptype
                ! move a cell
                call cell_move(delaunay_ptr,RTI,mcmc_set,imove,point,bnd_box,lerr)
                ! convert voronoi to grid
                if(lerr) then
                    call kdtree_to_grid(RTI,grid,bnd_box,model)
                    !call log_msg('move')
                    !call check_grid(RTI,grid,model)
                    ! calculate travel time and likelihood value
                    call likelihood(dat,model,RTI,bnd_box,like_set,like)
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
                        !samples(iter)%step = ptype
                        !samples(iter)%accepted = accepted
                        !samples(iter)%vindex = imove
                        !samples(iter)%ncells = delaunay_size(delaunay_ptr)
                        !samples(iter)%like = like%like
                        !samples(iter)%misfit = like%misfit
                        !samples(iter)%unweighted_misfit = like%unweighted_misfit
                        !samples(iter)%noise0 = 0
                        !samples(iter)%noise1 = 0
                        !samples(iter)%coord = [point%x,point%y,point%z]
                        !samples(iter)%values = 0
                        RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                    endif
                endif

                ! debug
                ! if not accepted firstly, start delayed rejection
            case(4)
                ! count
                RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1

                ! using hmc step
                call prepare_hmc(mcmc_set,RTI,dat(1)%src, hmc_set,qvals_in,ndim)
                call hmc_step(dat,RTI,like_set,hmc_set,qvals_in(1:ndim), like, accepted)
                ! accept or not based on acceptance ratio
                if(accepted) then
                    ! rebuild delaunay triangulation
                    !call array_points(points,RTI%points)
                    !call array_parameters(parameters,RTI%parameters)
                    !call cgal_delaunay_build(points, parameters, RTI%ncells, delaunay_ptr )
                    ! update model
                    do i = 1, grid%nx
                        do j = 1, grid%ny
                            do k = 1, grid%nz
                                model%vp(k,j,i) = RTI%parameters(1,RTI%sites_id(k,j,i))
                                model%vs(k,j,i) = RTI%parameters(2,RTI%sites_id(k,j,i))
                                model%rho(k,j,i) = RTI%parameters(3,RTI%sites_id(k,j,i))
                            enddo
                        enddo
                    enddo
                    !call log_msg('hmc')
                    !call check_grid(RTI,grid,model)
                    RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                endif
                ! if not accepted, delay reject
                ! debug
            case(5)
                ! count
                RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
                !samples(iter)%step = ptype
                ! change the data noise sigma
                call bsigma_change(RTI,mcmc_set,iperiod,lerr)
                if(.not.lerr) goto 100
                ! calculate likelihood
                call noise_likelihood(dat,RTI,like_set,like)
                ! calculate acceptance ratio
                alpha = minval([log(one), -like%like+like_copy%like])
                ! if accepted
                if(log(random)<alpha)then
                    accepted = .true.
                    ! update accepted ratio and samples
                    !samples(iter)%step = ptype
                    !samples(iter)%accepted = accepted
                    !samples(iter)%vindex = iperiod
                    !samples(iter)%ncells = delaunay_size(delaunay_ptr)
                    !samples(iter)%like = like%like
                    !samples(iter)%misfit = like%misfit
                    !samples(iter)%unweighted_misfit = like%unweighted_misfit
                    !samples(iter)%noise0 = RTI%bnoise0(iperiod)
                    !samples(iter)%noise1 = RTI%bnoise1(iperiod)
                    !samples(iter)%coord = 0
                    !samples(iter)%values = 0
                    RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                endif

            case(6)
                ! count
                RTI%samplecount(ptype) = RTI%samplecount(ptype) + 1
                !samples(iter)%step = ptype
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
                    !samples(iter)%step = ptype
                    !samples(iter)%accepted = accepted
                    !samples(iter)%vindex = iperiod
                    !samples(iter)%ncells = delaunay_size(delaunay_ptr)
                    !samples(iter)%like = like%like
                    !samples(iter)%misfit = like%misfit
                    !samples(iter)%unweighted_misfit = like%unweighted_misfit
                    !samples(iter)%noise0 = RTI%snoise0(iperiod)
                    !samples(iter)%noise1 = RTI%snoise1(iperiod)
                    !samples(iter)%coord = 0
                    !samples(iter)%values = 0
                    RTI%acceptedcount(ptype) = RTI%acceptedcount(ptype) + 1
                endif

            case default 
                call exception_raiseError("wrong proposal type (1 for birth, 2&
                    &for death, 3 for move, 4 for velocity change, 5 for surface wave noise change, 6 for body wave noise change&
                    &)")
              end select
            ! exit of select structure, might be changed to exit instead of goto
            ! in fortran 2008
            100 continue     
            ! if accepteded, update delaunay triangulation, vp, vs, rho, vel,
            ! like
            if(.not.accepted) then
                call cgal_delaunay_copy(delaunay_ptr_copy, delaunay_ptr)
                model = model_copy
                like = like_copy
                RTI%ncells = ncells_copy
                RTI%points = points_copy
                RTI%parameters = parameters_copy
                RTI%bnoise0 = bnoise0_copy
                RTI%bnoise1 = bnoise1_copy
                RTI%snoise0 = snoise0_copy
                RTI%snoise1 = snoise1_copy
                RTI%locations = sources_copy
                RTI%sites_id = sites_id_copy
            endif

            ! if icylce = 1, redoing the step thus skip next section
            if(icycle==1) cycle

            hmcsamples%ncells(iter) = RTI%ncells
            hmcsamples%points(:,:,iter) = RTI%points
            hmcsamples%parameters(:,:,iter) = RTI%parameters
            hmcsamples%negative_lglike(iter) = like%like
            hmcsamples%weighted_misfits(iter) = like%misfit
            hmcsamples%unweighted_misfits(iter) = like%unweighted_misfit
            hmcsamples%bnoise0(:,iter) = RTI%bnoise0
            hmcsamples%bnoise1(:,iter) = RTI%bnoise1
            hmcsamples%snoise0(:,iter) = RTI%snoise0
            hmcsamples%snoise1(:,iter) = RTI%snoise1

            if(mcmc_set%locate==1)then
                hmcsamples%locations(:,:,iter) = RTI%locations
            endif
            !samples(iter)%ncells = delaunay_size(delaunay_ptr)
            !RTI%ncells = samples(iter)%ncells
            !samples(iter)%like = like%like
            !samples(iter)%misfit = like%misfit
            !samples(iter)%unweighted_misfit = like%unweighted_misfit

            ! if after burn-in, begin sampling
            if(RTI%sampletotal > mcmc_set%burn_in .and.&
            mod(RTI%sampletotal-mcmc_set%burn_in,mcmc_set%thin)==0)then
                if(RTI%sampletotal-mcmc_set%burn_in==mcmc_set%thin)then
                    call log_msg('###########################')
                    call log_msg('Start sampling here: '//itoa(RTI%sampletotal))
                endif
                call stat_rti(RTI,model,dat,mcmc_set)
            endif
            ! write run-time info into files
            if(mcmc_set%runtime_step > 0 .and. mod(iter,mcmc_set%runtime_step) == 0)then
                call write_hmcsamples(hmcsamples,trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_sample_'//itoa(mcmc_set%processor)& 
                    &//'.dat',RTI%nsampled+1,iter)
                RTI%nsampled = iter
                call write_info(RTI,trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_info_'//itoa(mcmc_set%processor)//'.dat')
                call write_sources(RTI%locations,trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_sources_'//itoa(mcmc_set%processor)//'.dat')
                if(mcmc_set%locate == 1)then
                    nthin = RTI%nthin
                    if(RTI%nthin ==0) nthin = 1
                    call write_sources(RTI%aveL/nthin,trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_mean_sources_'//itoa(mcmc_set%processor)//'.dat')
                    call write_sources(RTI%stdL/nthin,trim(resume_dir)//FILEDIR_SEPARATOR//'run_time_var_sources_'//itoa(mcmc_set%processor)//'.dat')
                endif
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
                call log_msg( 'Number of cells: '//itoa(int(delaunay_size(delaunay_ptr))) )
                call log_msg( '-loglikelihood: '//rtoa(like%like) )
                call log_msg( 'Weighted misfit: '//rtoa(like%misfit) )
                call log_msg( 'Unweighted misfit: '//rtoa(like%unweighted_misfit) )
                call log_msg('')
                call log_msg( 'Accpetance rate for velocity' )
                call log_msg( 'ARV '//rtoa(RTI%acceptedcount(4)/real(RTI%samplecount(4),kind=ii10)) )
                call log_msg( 'Accpetance rate for position' )
                call log_msg( 'ARP '//rtoa(RTI%acceptedcount(3)/real(RTI%samplecount(3),kind=ii10)) )
                call log_msg( 'Accpetance rate for birth' )
                call log_msg( 'ARB '//rtoa(RTI%acceptedcount(1)/real(RTI%samplecount(1),kind=ii10)) )
                call log_msg( 'Accpetance rate for death' )
                call log_msg( 'ARD '//rtoa(RTI%acceptedcount(2)/real(RTI%samplecount(2),kind=ii10)) )
                call log_msg( 'Accpetance rate for sigma change of body waves' )
                call log_msg( 'ARS '//rtoa(RTI%acceptedcount(5)/real(RTI%samplecount(5),kind=ii10)) )
                call log_msg( 'Accpetance rate for sigma change of surface waves' )
                call log_msg( 'ARS '//rtoa(RTI%acceptedcount(6)/real(RTI%samplecount(6),kind=ii10)) )
                call log_msg( 'Accpetance rate for source locations change' )
                call log_msg( 'ARS '//rtoa(RTI%acceptedcount(7)/real(RTI%samplecount(7),kind=ii10)) )
                call log_msg( 'Rate of bad models' )
                call log_msg( 'Bad dispersion curves '//rtoa(RTI%num_bad_model/real(RTI%sampletotal,kind=ii10)) )
                call log_msg( 'Bad rays '//rtoa(RTI%num_bad_ray/real(RTI%sampletotal,kind=ii10)) )
                call log_msg( '---------------------------------------------------------' )
            endif
        enddo

        !call cgal_delaunay_write(delaunay_ptr,'last_delaunay_'//itoa(mcmc_set%processor)//'.dat')
        call cgal_delaunay_delete(delaunay_ptr)
        call cgal_delaunay_delete(delaunay_ptr_copy)

        ! write out the last delaunay triangulation
        ! write both a text and binary file, so no suffix specified
        call write_vertices(RTI,trim(last_dir)//FILEDIR_SEPARATOR//'last_vertices_'//itoa(mcmc_set%processor) )
        ! write last parameters to file
        call write_info(RTI,trim(last_dir)//FILEDIR_SEPARATOR//'last_info_'//itoa(mcmc_set%processor)//'.dat')
        call write_mean(RTI,trim(last_dir)//FILEDIR_SEPARATOR//'last_average_'//itoa(mcmc_set%processor)//'.dat')
        call write_var(RTI,trim(last_dir)//FILEDIR_SEPARATOR//'last_var_'//itoa(mcmc_set%processor)//'.dat')
        call write_sources(RTI%locations,trim(last_dir)//FILEDIR_SEPARATOR//'last_sources_'//itoa(mcmc_set%processor)//'.dat')
        if(mcmc_set%locate == 1)then
            nthin = RTI%nthin
            if(RTI%nthin ==0) nthin = 1
            call write_sources(RTI%aveL/nthin,trim(last_dir)//FILEDIR_SEPARATOR//'last_mean_sources_'//itoa(mcmc_set%processor)//'.dat')
            call write_sources(RTI%stdL/nthin,trim(last_dir)//FILEDIR_SEPARATOR//'last_var_sources_'//itoa(mcmc_set%processor)//'.dat')
        endif

        call write_hmcsamples(hmcsamples,trim(results_dir)//FILEDIR_SEPARATOR//'hmcsamples_'//itoa(mcmc_set%processor)//'.dat')
        call write_likelihood_hmc(trim(results_dir)//FILEDIR_SEPARATOR//'likelihood_'//itoa(mcmc_set%processor)//'.dat',hmcsamples)
        call write_number_of_cells_hmc(trim(results_dir)//FILEDIR_SEPARATOR//'ncells_'//itoa(mcmc_set%processor)//'.dat',hmcsamples)

        return

    end subroutine hmcmc

    subroutine prepare_hmc(mcmc_set, RTI, src, hmc_set, qvals_in, ndim)
        implicit none
        type(T_MCMC_SET), intent(in) :: mcmc_set
        type(T_RUN_INFO), intent(in) :: RTI
        real(kind=ii10), dimension(:,:), intent(in) :: src
        type(T_HMC_SET), intent(inout) :: hmc_set
        real(kind=ii10), dimension(:), intent(out) :: qvals_in
        integer, intent(out) :: ndim

        integer ncells, nsrc

        ncells = RTI%ncells
        nsrc = size(RTI%locations,2)
        hmc_set%priortype = mcmc_set%priortype
        hmc_set%n_leapfrogs = mcmc_set%n_leapfrogs

        ! prepare qvalues for hmc
        select case(mcmc_set%datatype)
        case (0)
            qvals_in(1:ncells) = RTI%parameters(1,1:ncells)
            hmc_set%stepsize(1:ncells) = mcmc_set%vp_step
            hmc_set%diagonal_matrix(1:ncells) = mcmc_set%vp_mass
            hmc_set%pm1(1:ncells) = mcmc_set%vpmin
            hmc_set%pm2(1:ncells) = mcmc_set%vpmax
            
            ndim = ncells
        case (1)
            qvals_in(1:ncells) = RTI%parameters(1,1:ncells)
            hmc_set%stepsize(1:ncells) = mcmc_set%vp_step
            hmc_set%diagonal_matrix(1:ncells) = mcmc_set%vp_mass
            hmc_set%pm1(1:ncells) = mcmc_set%vpmin
            hmc_set%pm2(1:ncells) = mcmc_set%vpmax
            qvals_in(ncells+1:2*ncells) = RTI%parameters(2,1:ncells)
            hmc_set%stepsize(ncells+1:2*ncells) = mcmc_set%vs_step
            hmc_set%diagonal_matrix(ncells+1:2*ncells) = mcmc_set%vs_mass
            hmc_set%pm1(ncells+1:2*ncells) = mcmc_set%vsmin
            hmc_set%pm2(ncells+1:2*ncells) = mcmc_set%vsmax
            ndim = ncells*2
        case (2)
            qvals_in(1:ncells) = RTI%parameters(2,1:ncells)
            hmc_set%stepsize(1:ncells) = mcmc_set%vs_step
            hmc_set%diagonal_matrix(1:ncells) = mcmc_set%vs_mass
            hmc_set%pm1(1:ncells) = mcmc_set%vsmin
            hmc_set%pm2(1:ncells) = mcmc_set%vsmax
            ndim = ncells
        case (3)
            qvals_in(1:ncells) = RTI%parameters(1,1:ncells)
            hmc_set%stepsize(1:ncells) = mcmc_set%vp_step
            hmc_set%diagonal_matrix(1:ncells) = mcmc_set%vp_mass
            hmc_set%pm1(1:ncells) = mcmc_set%vpmin
            hmc_set%pm2(1:ncells) = mcmc_set%vpmax
            qvals_in(ncells+1:2*ncells) = RTI%parameters(2,1:ncells)
            hmc_set%stepsize(ncells+1:2*ncells) = mcmc_set%vs_step
            hmc_set%diagonal_matrix(ncells+1:2*ncells) = mcmc_set%vs_mass
            hmc_set%pm1(ncells+1:2*ncells) = mcmc_set%vsmin
            hmc_set%pm2(ncells+1:2*ncells) = mcmc_set%vsmax
            ndim = ncells*2
        case default
            ! default
        endselect

        ! add source locations
        if(mcmc_set%locate==1)then
            qvals_in(ndim+1:4*nsrc+ndim) =&
            reshape(RTI%locations,[4*nsrc])
            hmc_set%stepsize(ndim+1:4*nsrc+ndim) = mcmc_set%src_step
            hmc_set%stepsize(ndim+4:4*nsrc+ndim:4) = mcmc_set%src_tstep
            hmc_set%diagonal_matrix(ndim+1:4*nsrc+ndim) = mcmc_set%src_mass
            hmc_set%pm1(ndim+1:4*nsrc+ndim:4) = src(1,:) - mcmc_set%xwidth
            hmc_set%pm1(ndim+2:4*nsrc+ndim:4) = src(2,:) - mcmc_set%ywidth
            hmc_set%pm1(ndim+3:4*nsrc+ndim:4) = src(3,:) - mcmc_set%zwidth
            hmc_set%pm1(ndim+4:4*nsrc+ndim:4) = src(4,:) - mcmc_set%twidth
            where(hmc_set%pm1(ndim+1:4*nsrc+ndim:4)<mcmc_set%grid%xmin)
                hmc_set%pm1(ndim+1:4*nsrc+ndim:4) = mcmc_set%grid%xmin
            endwhere
            where(hmc_set%pm1(ndim+2:4*nsrc+ndim:4)<mcmc_set%grid%ymin)
                hmc_set%pm1(ndim+2:4*nsrc+ndim:4)= mcmc_set%grid%ymin
            endwhere
            where(hmc_set%pm1(ndim+3:4*nsrc+ndim:4)<mcmc_set%grid%zmin)
                hmc_set%pm1(ndim+3:4*nsrc+ndim:4) = mcmc_set%grid%zmin
            endwhere
            hmc_set%pm2(ndim+1:4*nsrc+ndim:4) = src(1,:) + mcmc_set%xwidth
            hmc_set%pm2(ndim+2:4*nsrc+ndim:4) = src(2,:) + mcmc_set%ywidth
            hmc_set%pm2(ndim+3:4*nsrc+ndim:4) = src(3,:) + mcmc_set%zwidth
            hmc_set%pm2(ndim+4:4*nsrc+ndim:4) = src(4,:) + mcmc_set%twidth
            where(hmc_set%pm2(ndim+1:4*nsrc+ndim:4)>mcmc_set%grid%xmax)
                hmc_set%pm2(ndim+1:4*nsrc+ndim:4) = mcmc_set%grid%xmax
            endwhere
            where(hmc_set%pm2(ndim+2:4*nsrc+ndim:4)>mcmc_set%grid%ymax)
                hmc_set%pm2(ndim+2:4*nsrc+ndim:4) = mcmc_set%grid%ymax
            endwhere
            where(hmc_set%pm2(ndim+3:4*nsrc+ndim:4)>mcmc_set%grid%zmax)
                hmc_set%pm2(ndim+3:4*nsrc+ndim:4) = mcmc_set%grid%zmax
            endwhere
            ndim = ndim + 4*nsrc
        endif

    endsubroutine prepare_hmc

    subroutine write_hmcsamples(hmcsamples,filename, ibegin, iend)
        use m_utils, only : write_resume_unit
        implicit none
	    type(T_HMC_SAMPLES), intent(in) :: hmcsamples
    	character( len=* ), intent(in) :: filename
        integer, intent(in), optional :: ibegin, iend

	    integer i, begin_idx
	    integer nsamples, ncells, nsrc, ndim, nd, np
	    logical lexist

	    ! > open file for writing
	    inquire(file=filename,exist=lexist)
	    if(.not.lexist)then
	    	open(unit=write_resume_unit,file=filename,status='new',access='sequential')
	    else
	    	open(unit=write_resume_unit,file=filename,status='old',access='sequential',position='append')
	    endif

	    ! > write samples to the file sample by sample
        begin_idx = 1
        if(present(ibegin)) begin_idx = ibegin
	    nsamples = size(hmcsamples%ncells)
        if(present(iend)) nsamples = iend
        if(allocated(hmcsamples%locations))then
            nsrc = size(hmcsamples%locations,2)
            ndim = size(hmcsamples%locations,1)
        else
            nsrc = 0
        endif
        nd = size(hmcsamples%bnoise0,1)
        np = size(hmcsamples%snoise0,1)
        if(.not.lexist) write(write_resume_unit,*) nsrc, ndim, nd, np
        if(nsrc /= 0)then
	        do i = begin_idx, nsamples
                ncells = hmcsamples%ncells(i)
                write(write_resume_unit,*) hmcsamples%negative_lglike(i),&
                hmcsamples%weighted_misfits(i), hmcsamples%unweighted_misfits(i),&
                hmcsamples%ncells(i) 
                write(write_resume_unit,*) hmcsamples%points(:,1:ncells,i),&
                hmcsamples%parameters(:,1:ncells,i), hmcsamples%locations(:,:,i),&
                hmcsamples%bnoise0(:,i), hmcsamples%bnoise1(:,i),&
                hmcsamples%snoise0(:,i), hmcsamples%snoise1(:,i) 
	        enddo
        else
	        do i = begin_idx, nsamples
                ncells = hmcsamples%ncells(i)
                write(write_resume_unit,*) hmcsamples%negative_lglike(i),&
                hmcsamples%weighted_misfits(i), hmcsamples%unweighted_misfits(i),&
                hmcsamples%ncells(i) 
                write(write_resume_unit,*) hmcsamples%points(:,1:ncells,i),&
                hmcsamples%parameters(:,1:ncells,i), &
                hmcsamples%bnoise0(:,i), hmcsamples%bnoise1(:,i),&
                hmcsamples%snoise0(:,i), hmcsamples%snoise1(:,i) 
	        enddo
        endif
	    close(write_resume_unit)

    end subroutine

    subroutine read_hmcsamples(hmcsamples,filename)
        use m_utils, only : write_resume_unit
        implicit none
	    type(T_HMC_SAMPLES), intent(inout) :: hmcsamples
    	character( len=* ), intent(in) :: filename

	    integer i
	    integer ncells, nsrc, ndim, nd, np
	    logical lexist

	    ! > open file for writing
	    !inquire(file=filename,exist=lexist)
	    open(unit=write_resume_unit,file=filename,status='old',access='sequential',action='read')

	    ! > read samples from the file sample by sample
        if(allocated(hmcsamples%locations))then
            nsrc = size(hmcsamples%locations,2)
            ndim = size(hmcsamples%locations,1)
        else
            nsrc = 0
        endif
        nd = size(hmcsamples%bnoise0,1)
        np = size(hmcsamples%snoise0,1)
        read(write_resume_unit,*) nsrc, ndim, nd, np
        if(nsrc /= 0)then
	        do i = 1, size(hmcsamples%ncells)
                read(write_resume_unit,*) hmcsamples%negative_lglike(i),&
                hmcsamples%weighted_misfits(i), hmcsamples%unweighted_misfits(i),&
                hmcsamples%ncells(i) 
                ncells = hmcsamples%ncells(i)
                read(write_resume_unit,*) hmcsamples%points(:,1:ncells,i),&
                hmcsamples%parameters(:,1:ncells,i), hmcsamples%locations(:,:,i),&
                hmcsamples%bnoise0(:,i), hmcsamples%bnoise1(:,i),&
                hmcsamples%snoise0(:,i), hmcsamples%snoise1(:,i) 
	        enddo
        else
	        do i = 1, size(hmcsamples%ncells)
                ncells = hmcsamples%ncells(i)
                read(write_resume_unit,*) hmcsamples%negative_lglike(i),&
                hmcsamples%weighted_misfits(i), hmcsamples%unweighted_misfits(i),&
                hmcsamples%ncells(i) 
                ncells = hmcsamples%ncells(i)
                read(write_resume_unit,*) hmcsamples%points(:,1:ncells,i),&
                hmcsamples%parameters(:,1:ncells,i), &
                hmcsamples%bnoise0(:,i), hmcsamples%bnoise1(:,i),&
                hmcsamples%snoise0(:,i), hmcsamples%snoise1(:,i) 
	        enddo
        endif
	    close(write_resume_unit)

    end subroutine

    subroutine write_likelihood_hmc(filename,samples)
        use m_utils, only : write_resume_unit, write_doubles
        implicit none
        character( len=* ), intent(in) :: filename
        type(T_HMC_SAMPLES), intent(in) :: samples

        integer i
        integer nsamples
        logical lexist

        ! > open file for writing
        inquire(file=filename,exist=lexist)
        if(.not.lexist)then
            open(unit=write_resume_unit,file=filename,status='new',access='sequential')
        else
            open(unit=write_resume_unit,file=filename,status='old',access='sequential',position='append')
        endif

        ! > write samples to the file sample by sample
        nsamples = size(samples%ncells)
        do i = 1, nsamples
                call write_doubles([samples%negative_lglike(i), samples%weighted_misfits(i), samples%unweighted_misfits(i)])
        enddo
        close(write_resume_unit)

    end subroutine

    subroutine write_number_of_cells_hmc(filename,samples)
        use m_utils, only : write_resume_unit, integer_format
        implicit none
        character( len=* ), intent(in) :: filename
        type(T_HMC_SAMPLES), intent(in) :: samples

        integer i
        integer nsamples
        logical lexist

        ! > open file for writing
        inquire(file=filename,exist=lexist)
        if(.not.lexist)then
            open(unit=write_resume_unit,file=filename,status='new',access='sequential')
        else
            open(unit=write_resume_unit,file=filename,status='old',access='sequential',position='append')
        endif

        ! > write samples to the file sample by sample
        nsamples = size(samples%ncells)
        do i = 1, nsamples
            write(write_resume_unit,integer_format(1)) samples%ncells(i)
        enddo
        close(write_resume_unit)

    end subroutine

end module m_hmcmc
