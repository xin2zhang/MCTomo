!
! Sample: post process these samples produced by
! reversible jump Markov Chain Monte Carlo algorithm
! 
module m_sample

    use m_exception, only       : exception_raiseError
    use m_logger, only          : log_msg
    use m_utils, only           : ii10, itoa, rtoa, newunit
    use m_settings, only        : T_GRID
    use run_info, only          : T_SAMPLE, T_RUN_INFO
    ! use all the available procedure, types and variable in the following
    ! module
    use cgal_delaunay
    use iso_c_binding
    implicit none

    private
    public :: march_sample
    public :: T_OUT, output_setup
    public :: write_output
    public :: read_samples
    public :: read_array2d


    ! static values
    integer, parameter, public :: temp_unit = 7
    integer, parameter, public :: like_unit = 8 
    integer, parameter, public :: sample_unit = 9

    type T_OUT
        integer :: nthin
        real(kind=16), dimension(:,:,:), allocatable :: aveS
        real(kind=16), dimension(:,:,:), allocatable :: stdS
        real(kind=16), dimension(:,:,:), allocatable :: aveP
        real(kind=16), dimension(:,:,:), allocatable :: stdP
        real(kind=ii10), dimension(:,:,:,:), allocatable :: postS
        !TODO : other output, like evidence of number of cells
    endtype

contains

    subroutine output_setup(output,grid,nv)
        implicit none
        type(T_OUT), intent(out) :: output
        type(T_GRID), intent(in) :: grid
        integer, intent(in) :: nv

        allocate( output%aveS(grid%nz,grid%ny,grid%nx) )
        allocate( output%stdS(grid%nz,grid%ny,grid%nx) )
        allocate( output%aveP(grid%nz,grid%ny,grid%nx) )
        allocate( output%stdP(grid%nz,grid%ny,grid%nx) )
        allocate( output%postS(nv,grid%nz,grid%ny,grid%nx) )

        output%nthin = 0
        output%aveS = 0.0
        output%stdS = 0.0
        output%aveP = 0.0
        output%stdP = 0.0
        output%postS = 0.0

        return
    end subroutine

    subroutine march_sample(sample,RTI)
        implicit none
        type(T_SAMPLE), intent(in)      :: sample
        type(T_RUN_INFO), intent(inout) :: RTI

        real(kind=ii10), dimension(3,size(RTI%points,2)) :: points_copy
        real(kind=ii10), dimension(3,size(RTI%parameters,2)) :: parameters_copy

        ! not accepted, return
        if(.not.sample%accepted) return

        select case (sample%step)
        case(1)
            RTI%ncells = RTI%ncells + 1
            RTI%points(:,RTI%ncells) = sample%coord
            RTI%parameters(:,RTI%ncells) = sample%values
        
        case(2)
            points_copy = RTI%points
            parameters_copy = RTI%parameters

            RTI%points(:,sample%vindex:RTI%ncells-1) = points_copy(:,sample%vindex+1:RTI%ncells)
            RTI%parameters(:,sample%vindex:RTI%ncells-1) = parameters_copy(:,sample%vindex+1:RTI%ncells)
            RTI%ncells = RTI%ncells - 1

        case(3)
            RTI%points(:,sample%vindex) = sample%coord

        case(4)
            RTI%parameters(:,sample%vindex) = sample%values

        case(5)
            RTI%bnoise0(sample%vindex) = sample%noise0
            RTI%bnoise1(sample%vindex) = sample%noise1

        case(6)
            RTI%snoise0(sample%vindex) = sample%noise0
            RTI%snoise1(sample%vindex) = sample%noise1

        case(7)
            RTI%locations(1:3,sample%vindex) = sample%coord
            RTI%locations(4,sample%vindex) = sample%values(1)

        case default
            call exception_raiseError("wrong proposal type (1 for birth, 2&
                 &for death, 3 for move and 4 for velocity change)")
        end select
        
        if(RTI%ncells/=sample%ncells)then
            call log_msg('Warning: number of cells not equal, force equaled')
            RTI%ncells = sample%ncells
        endif

        return
    end subroutine

    subroutine read_array2d(iunit,array2d, ntemps)
        use m_utils, only : double_format
        implicit none
        integer, intent(in) :: iunit
	    real(kind=ii10), dimension(:,:), intent(out) :: array2d
        integer, intent(out) :: ntemps

	    integer i, nd
        integer stat
	    integer nsamples

        ntemps = 0
        nd = size(array2d,1)
	    nsamples = size(array2d,2)
	    do i = 1, nsamples
            read(iunit,double_format(nd),iostat=stat) array2d(:,i)
            if(stat/=0) exit
            ntemps = ntemps + 1
	    enddo

    end subroutine

    subroutine read_samples(samples, nsamples)
        type(T_SAMPLE), dimension(:), intent(out) :: samples
        integer, intent(out) :: nsamples

        integer i
        integer stat

        nsamples = 0
        do i = 1, size(samples)
            read(sample_unit,iostat=stat) samples(i)%step, samples(i)%accepted,&
            samples(i)%vindex, samples(i)%ncells, samples(i)%misfit,&
            samples(i)%unweighted_misfit,samples(i)%like,  &
            samples(i)%coord, samples(i)%values, samples(i)%noise0, &
            samples(i)%noise1
            !read(iunit,iostat=stat) samples(i)
            if(stat /= 0) exit
            nsamples = nsamples + 1
        enddo

    end subroutine

    subroutine write_doubles(arr,str,n)
        use m_utils, only : write_resume_unit, double_format
        implicit none
        real(kind=16),dimension(:,:,:), intent(in) :: arr
        integer, intent(in),dimension(size(arr,2)),optional :: n
        character(len=*), intent(in),optional :: str
        integer :: i, j, k

        if(present(str)) write(write_resume_unit,'("'//trim(str)//'")')
        do i=1,size(arr,3)
            do j = 1, size(arr,2)
                write(write_resume_unit,double_format(size(arr,1))) arr(:,j,i)
            enddo
        end do

    end subroutine write_doubles
    
    subroutine write_output(output,prefix)
        use m_utils, only : write_resume_unit, results_dir, FILEDIR_SEPARATOR
        implicit none
        type(T_OUT), intent(in) :: output
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: str
        real(kind=16), dimension(:,:,:), allocatable :: stdS, stdP

        allocate( stdS(size(output%aveS,1),size(output%aveS,2),size(output%aveS,3)) )
        allocate( stdP(size(output%aveP,1),size(output%aveP,2),size(output%aveP,3)) )

        stdS = output%stdS - output%aveS**2
        stdP = output%stdP - output%aveP**2

        where(stdS<0)
            stdS = 0
        endwhere
        where(stdP<0)
            stdP=0
        endwhere

        stdP = sqrt(stdP)
        stdS = sqrt(stdS)

        ! average
        open(unit=write_resume_unit,file=trim(results_dir)//FILEDIR_SEPARATOR//prefix//'_average.dat',status='replace',action='write')
        str = itoa(size(output%aveS,1))//' '//itoa(size(output%aveS,2))&
            &//' '//itoa(size(output%aveS,3)*2)
        call write_doubles(output%aveP,str)
        call write_doubles(output%aveS)
        close(write_resume_unit)
        ! standard deviation
        open(unit=write_resume_unit,file=trim(results_dir)//FILEDIR_SEPARATOR//prefix//'_std.dat',status='replace',action='write')
        str = itoa(size(output%aveS,1))//' '//itoa(size(output%aveS,2))&
            &//' '//itoa(size(output%aveS,3)*2)
        call write_doubles(stdP,str)
        call write_doubles(stdS)
        close(write_resume_unit)
        return
    end subroutine

end module m_sample

! main program
program sample

    use iso_c_binding

    use m_exception,                    only : exception_raiseError
    use m_logger,                       only : log_msg, log_startup,log_shutdown
    use m_utils,                        only : itoa, rtoa, STRLEN, newunit,&
                                                FILEDIR_SEPARATOR, results_dir, eps, ini_dir, log_dir
    use kdtree2_module,                 only : kdtree2, kdtree2_result, kdtree2_create,&
                                                kdtree2_n_nearest
    use m_settings,                     only : T_MCMC_SET, T_GRID, grid_setup, T_MOD, mod_setup
    use like_settings,                  only : T_DATA, T_LIKE_SET, read_data
    use run_info,                       only : T_SAMPLE, T_RUN_INFO,&
                                                init_samples, init_run_info
    use m_mcmc,                         only : kdtree_to_grid, ImportanceSample
    use read_write,                     only : read_vertices, write_vertices, read_temperatures0
    use cgal_delaunay,                  only : p3, d3

    ! use all the available procedure, types and variable in the following
    ! module
    use m_sample

    implicit none

    ! static value
    integer, parameter :: N = 1000000
    integer, parameter :: NCOLLECT = 3000000

    ! files
    character(len=STRLEN), parameter :: fsamples = 'samples'
    character(len=:), allocatable :: dirname

    ! likelihood settings
    type(T_LIKE_SET)    :: like_set

    ! data
    type(T_DATA), dimension(:), allocatable :: dat
    ! Voronoi tesselation
    type(T_MCMC_SET)    :: mcmc_set

    ! grid tesselation
    type(T_GRID) :: grid
    type(T_MOD)  :: model
    type(d3) :: bnd_box(2)

    ! samples
    type(T_SAMPLE), dimension(:), allocatable  :: samples
    type(T_RUN_INFO)  :: RTI
    !type(T_SAMPLE), dimension(:,:,:), allocatable  :: thinnedSamples
    real(c_double), dimension(:,:), allocatable :: lglikelihoods
    real(c_double), dimension(:,:), allocatable :: tlglikelihoods

    ! tempering
    real(c_double), dimension(:,:), allocatable :: temperatures
    real(c_double), dimension(:), allocatable :: temperatures0
    real(kind=16), dimension(:,:,:,:), allocatable :: tvp, tvs
    real(kind=16), dimension(:,:,:,:), allocatable :: tvarp, tvars
    real(kind=16), dimension(:), allocatable :: tnormalise
    integer, dimension(:), allocatable :: tnthin
    integer idx, ntemperatures, nthin, totalSamples ! num of temperatures not equal to 1
    real(kind=c_double), dimension(:), allocatable :: norm
    real(kind=c_double) :: tnorm
    real(kind=16) :: weight
    real(kind=16), dimension(:,:), allocatable :: importance

    ! sampling output
    integer         :: nv
    real(c_double), parameter  :: dvel = 0.01
    type(T_OUT)     :: output

    ! chains
    integer, dimension(:), allocatable :: chains
    integer, dimension(:), allocatable :: burn_in

    ! unit
    integer iunit
    ! iterations
    integer     nsamples, nsamples2, nlikelihoods
    integer     ncount
    integer     nchains
    integer     ichain
    integer     iter, i 

    ! debug
    character(100) :: debug_info

    !namelist /datafiles/ fbase, fsamples, np, dvel
    namelist /grid_settings/ grid
    namelist /likelihood_settings/ like_set
    namelist /mcmc_settings/ mcmc_set

    ! set up log file
    call log_startup(trim(log_dir)//FILEDIR_SEPARATOR//'Sample.log')

    ! read input file
    call log_msg('Reading input file...')
    open(unit=newunit(iunit),file='MCTomo.inp',status='old',delim='apostrophe')
    !read(iunit, nml=datafiles)
    read(iunit, nml=grid_settings)
    read(iunit, nml=likelihood_settings)
    read(iunit, nml=mcmc_settings)
    close(iunit)

    ! initial some variable of settings
    call grid_setup(grid)
    mcmc_set%grid = grid
    mcmc_set%datatype = like_set%datatype
    bnd_box(1) = d3(grid%xmin,grid%ymin,grid%zmin)
    bnd_box(2) = d3(grid%xmax,grid%ymax,grid%zmax)
    
    ! read data
    allocate(dat(2))
    call read_data(dat,like_set)

    call init_run_info(RTI,dat,mcmc_set)
    ! initialize and allocate several variables related to grid
    allocate( model%vp(grid%nz,grid%ny,grid%nx) )
    allocate( model%vs(grid%nz,grid%ny,grid%nx) )
    allocate( model%rho(grid%nz,grid%ny,grid%nx) )
    model%vp = 0
    model%vs = 0
    model%rho = 0

    ! if parallel tempering, allocate and initialise vp, vs for importance
    ! sampling TODO: using a solution of reverse logistic regression would be
    ! better
    ntemperatures = 0
    if(mcmc_set%tempering == 1)then
        ntemperatures = mcmc_set%number_of_temperatures - mcmc_set%number_of_1s
        allocate(tvp(grid%nz,grid%ny,grid%nx,ntemperatures))
        allocate(tvs(grid%nz,grid%ny,grid%nx,ntemperatures))
        allocate(tvarp(grid%nz,grid%ny,grid%nx,ntemperatures))
        allocate(tvars(grid%nz,grid%ny,grid%nx,ntemperatures))
        tvp = 0
        tvs = 0
        tvarp = 0
        tvars = 0
        allocate(temperatures0(ntemperatures))
        call read_temperatures0(trim(mcmc_set%temperfile),&
             temperatures0,1)
        allocate(tnormalise(ntemperatures))
        tnormalise = 0.0
        allocate(tnthin(ntemperatures))
        tnthin = 0
        allocate(norm(ntemperatures))
        norm = 0
    endif
    

    ! read chains number
    open(unit=newunit(iunit),file='chains.inp',status='old')
    read(iunit, *) nchains
    allocate( chains(nchains) )
    allocate( burn_in(nchains) )
    do i = 1, nchains
        read(iunit, *) chains(i), burn_in(i)
    enddo
    close(iunit)

    allocate(tlglikelihoods(NCOLLECT*nchains/mcmc_set%thin,ntemperatures))
    allocate(importance(NCOLLECT*nchains/mcmc_set%thin,ntemperatures))
    allocate(temperatures(4,N))
    allocate(lglikelihoods(3,N))
    dirname = adjustl(trim(results_dir))//FILEDIR_SEPARATOR
    ! read negtive log likelihood and tempratures to determine the weight for
    ! importance sampling
    if(mcmc_set%tempering == 1)then
        do ichain = 1, nchains
            ! open temperature file and sample file
            call log_msg('Processing lglikelihoods and temperatures for chain '//itoa(chains(ichain))//' ...')
            open(unit=temp_unit,file=dirname//'temperatures_'//itoa(chains(ichain))//'.dat',status='old',action='read')
            open(unit=like_unit,file=dirname//'likelihood_'//itoa(chains(ichain))//'.dat',status='old',action='read')
            read(like_unit,*)
            ncount = 0
            do
                call read_array2d(temp_unit,temperatures,nsamples2)
                call read_array2d(like_unit,lglikelihoods,nlikelihoods)
                iter = 0
                do while(iter<nsamples2)
                    iter = iter + 1
                    ncount = ncount + 1
                    if(ncount >= burn_in(ichain) .and. mod(ncount-burn_in(ichain), mcmc_set%thin) == 0 &
                        )then
                        if(abs(temperatures(4,iter) - 1.0)>eps)then
                            idx=minloc(abs(temperatures0-temperatures(4,iter)),1)
                            tnthin(idx) = tnthin(idx) + 1
                            tlglikelihoods(tnthin(idx),idx) = -lglikelihoods(1,iter)
                        endif
                    endif
                enddo

                ! if EOF has reached
                if(nsamples2 < N) exit

            enddo
            close(temp_unit)
            close(like_unit)
        enddo
    endif

    ! get the importance and norm
    do i = 1, ntemperatures
        norm(i) = maxval(tlglikelihoods(1:tnthin(i),i))
        importance(1:tnthin(i),i) = tlglikelihoods(1:tnthin(i),i)*(1-1.0/temperatures0(i)) - norm(i)*(1-1.0/temperatures0(i)) 
        tnorm = log(sum(exp((1-1.0/temperatures0(i))*(tlglikelihoods(1:tnthin(i),i)-norm(i)))))
        importance(1:tnthin(i),i) = importance(1:tnthin(i),i) - tnorm 
        importance(1:tnthin(i),i) = exp(importance(1:tnthin(i),i))
    enddo


    ! loop over all chains
    allocate(samples(N))
    nv = ceiling((mcmc_set%vsmax - mcmc_set%vsmin)/dvel)
    call output_setup(output,grid,nv)
    call init_samples(samples)
    temperatures = 1.0
    tnthin = 0
    do ichain = 1, nchains
        call log_msg('Processing chain '//itoa(chains(ichain))//' ...')
        ! read the initial vertices from file
        call read_vertices(RTI,trim(ini_dir)//FILEDIR_SEPARATOR//'InitialSample_'//itoa(chains(ichain))//'.dat' )

        ! open temperature file and sample file
        if(mcmc_set%tempering == 1)then
            open(unit=temp_unit,file=dirname//'temperatures_'//itoa(chains(ichain))//'.dat',status='old',action='read')
        endif
        open(unit=sample_unit,file=dirname//adjustl(trim(fsamples))//'_'//itoa(chains(ichain))//'.out',&
            status='old',access='stream',action='read')

        ncount = 0
        do
            ! read N tempratures and N samples if possible
            call read_samples(samples,nsamples)
            if(mcmc_set%tempering==1)then
                call read_array2d(temp_unit,temperatures,nsamples2)
                if(nsamples /= nsamples2)&
                    call exception_raiseError('The number of temperatures not equal to the number of samples!')
            endif

            iter = 0
            do while(iter<nsamples)
                iter = iter + 1
                ncount = ncount + 1
                call march_sample(samples(iter),RTI)
                !if(samples%accepted)then
                !    call write_vertices(RTI,'sampleVertices_'//itoa(iter)//'.dat')
                !endif
                if(ncount >= burn_in(ichain) .and. mod(ncount-burn_in(ichain), mcmc_set%thin) == 0 &
                    )then
                    if(ncount-burn_in(ichain) == 0)then
                        call log_msg('Start sampling '//itoa(ncount))
                        !if(ichain==1) norm = -1.05*samples(iter)%like
                        !cnorm = -0.95*samples(iter)%like

                    endif
                    call kdtree_to_grid(RTI,grid,bnd_box,model)
                    if(abs(temperatures(4,iter) - 1.0)<eps)then
                        output%aveS = output%aveS + model%vs
                        output%stdS = output%stdS + model%vs**2
                        output%aveP = output%aveP + model%vp
                        output%stdP = output%stdP + model%vp**2
                        output%nthin = output%nthin + 1
                    else
                        idx=minloc(abs(temperatures0-temperatures(4,iter)),1)
                        weight = ImportanceSample(temperatures(4,iter),-samples(iter)%like,norm(idx))
                        tvp(:,:,:,idx) = tvp(:,:,:,idx) + weight * model%vp
                        tvs(:,:,:,idx) = tvs(:,:,:,idx) + weight * model%vs
                        tvarp(:,:,:,idx) = tvarp(:,:,:,idx) + weight * model%vp**2
                        tvars(:,:,:,idx) = tvars(:,:,:,idx) + weight * model%vs**2
                        tnormalise(idx) = tnormalise(idx) + weight
                        tnthin(idx) = tnthin(idx) + 1
                        !print *, 'like: ', -samples(iter)%like, 'norm: ', norm(idx), 'weight1: ', weight, 'weight2: ', importance(tnthin(idx),idx)
                    endif
                endif
            enddo

            ! if EOF has reached
            if(nsamples < N) exit

        enddo

        close(sample_unit)
        if(mcmc_set%tempering==1) close(temp_unit)

        call write_vertices(RTI,trim(results_dir)//FILEDIR_SEPARATOR//'last_vertices_'//itoa(chains(ichain)) )
        call log_msg('Total samples: '//itoa(ncount) // ' on chain &
        &'//itoa(chains(ichain)) )
    enddo

    ! first, only sample those chains with temperatures = 1
    nthin = output%nthin
    if(output%nthin==0) nthin=1
    output%aveS = output%aveS/nthin
    output%stdS = output%stdS/nthin
    output%aveP = output%aveP/nthin
    output%stdP = output%stdP/nthin
    ! write out the information
    call log_msg('Output posterior information...')
    call write_output(output,'Unweighted')

    ! Then, sample all those chains with importance sampling
    if(mcmc_set%tempering==1)then
        totalSamples = sum(tnthin) + output%nthin
        weight = output%nthin*1.0/totalSamples
        output%aveP = weight * output%aveP 
        output%aveS = weight * output%aveS 
        output%stdP = weight * output%stdP 
        output%stdS = weight * output%stdS 
        do i = 1, ntemperatures
            weight = tnthin(i)*1.0/totalSamples
            if(tnthin(i)>0 .and. tnormalise(i)>0)then
                output%aveP = output%aveP + tvp(:,:,:,i)/tnormalise(i) * weight
                output%aveS = output%aveS + tvs(:,:,:,i)/tnormalise(i) * weight
                output%stdP = output%stdP + tvarp(:,:,:,i)/tnormalise(i) * weight
                output%stdS = output%stdS + tvars(:,:,:,i)/tnormalise(i) * weight
            endif
        enddo
        call log_msg('Output posterior information from importance sampling...')
        call write_output(output,'Weighted')
    endif

    call log_shutdown()

end program

