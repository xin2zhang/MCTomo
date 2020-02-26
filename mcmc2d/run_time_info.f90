! run time information
!
module run_info

    use m_exception,    only : exception_raiseError
    use m_utils, only    : newunit, itoa
    use m_settings, only : T_MCMC_SET
    use like_settings, only: T_DATA
    use mt19937, only    : unirand
    use iso_c_binding

    implicit none

    private
    public :: T_RUN_INFO
    public :: init_run_info, read_info, write_info
    public :: T_SAMPLE
    public :: init_sample, init_samples
    public :: read_samples, write_samples
    ! run time info in case to recover later
    ! It is better to write info during run in case to recover run after
    ! iterruption

    ! > number of proposal types
    integer, parameter :: NTYPES = 7
    type T_RUN_INFO
        integer(8)      randcount
        integer         nsampled, nsamples
        integer         like_count
        integer(8)      sampletotal
        integer         nthin
        integer(8),dimension(NTYPES) :: samplecount, acceptedcount
        
        ! bad model due to dispersion curve modelling error and
        ! fast marching ray tracing error
        integer(c_int) num_bad_model, num_bad_ray

        ! tempering
        integer number_of_chains
        integer number_of_temperatures
        integer accepted_tempering, total_tempering
        real(c_double), dimension(:), allocatable :: temperature_values0,temperature_values
        integer, dimension(:), allocatable :: temperature_indices0,temperature_indices

        ! model information
        integer(c_size_t) ncells
        real(c_double), dimension(:,:), allocatable :: points
        real(c_double), dimension(:,:), allocatable :: parameters
        integer, dimension(:,:), allocatable :: sites_id
        real(c_double), dimension(:), allocatable   :: sigma_n0, sigma_n1
        real(c_double), dimension(:), allocatable   :: n0_min, n0_max, n1_min, n1_max
        real(c_double), dimension(:), allocatable   :: snoise0, snoise1
        real(c_double), dimension(:,:), allocatable   :: locations
        real(c_double), dimension(:,:), allocatable :: aveP, aveS
        real(c_double), dimension(:,:), allocatable :: stdP, stdS
        real(c_double), dimension(:,:), allocatable :: aveL, stdL
    end type

    !> derived type of each sample
    !  use c_double real type to keep consistence with other language( eg.
    !  c, cpp, matlab, python) to post process samples
    type T_SAMPLE
        integer(c_int)  step    ! propose type
	    logical(c_bool) accepted  ! accepted or not
	    integer(c_int)  vindex  ! index of vertex
        integer(c_size_t)  ncells  ! number of cells
	    real( kind=c_double ) misfit, unweighted_misfit
	    real( kind=c_double ) like
	    real( kind=c_double ), dimension(2) :: coord
	    real( kind=c_double ), dimension(3) :: values
        real( kind=c_double )  :: noise0
        real( kind=c_double )  :: noise1
        !TODO: noise update
	    !real( kind=ii10 ), dimension(3) :: sigma
	    !integer(c_long), dimension(4) :: counts
    endtype

contains
    subroutine init_run_info(RTI,dat,set)
        use iso_c_binding, only : c_ptr
        implicit none
        type(T_RUN_INFO), intent(inout) :: RTI
        type(T_DATA), dimension(:), intent(in) :: dat
        type(T_MCMC_SET), intent(in) :: set

        integer np, nd
        integer i

        RTI%randcount = 0
        RTI%nsampled = 0
        RTI%nsamples = set%nsamples
        RTI%like_count = 0
        RTI%sampletotal = 0
        RTI%samplecount = 0
        RTI%acceptedcount = 0
        RTI%ncells = 0
        RTI%nthin = 0
        RTI%num_bad_model = 0
        RTI%num_bad_ray = 0
    
        allocate(RTI%points(2,set%ncell_max))
        allocate(RTI%parameters(3,set%ncell_max))
        RTI%points = 0
        RTI%parameters = 0


        ! allocate and initialise the noise
        if(set%datatype==1)then
            allocate( RTI%snoise0(dat(1)%np*2) )
            allocate( RTI%snoise1(dat(1)%np*2) )
            allocate( RTI%sigma_n0(dat(1)%np*2) )
            allocate( RTI%sigma_n1(dat(1)%np*2) )
            allocate( RTI%n0_min(dat(1)%np*2) )
            allocate( RTI%n0_max(dat(1)%np*2) )
            allocate( RTI%n1_min(dat(1)%np*2) )
            allocate( RTI%n1_max(dat(1)%np*2) )
        elseif(set%datatype==3)then
            allocate( RTI%snoise0(dat(1)%np*2+dat(2)%np) )
            allocate( RTI%snoise1(dat(1)%np*2+dat(2)%np) )
            allocate( RTI%sigma_n0(dat(1)%np*2+dat(2)%np) ) 
            allocate( RTI%sigma_n1(dat(1)%np*2+dat(2)%np) )
            allocate( RTI%n0_min(dat(1)%np*2+dat(2)%np) )
            allocate( RTI%n0_max(dat(1)%np*2+dat(2)%np) )
            allocate( RTI%n1_min(dat(1)%np*2+dat(2)%np) )
            allocate( RTI%n1_max(dat(1)%np*2+dat(2)%np) )
        else
            allocate( RTI%snoise0(dat(1)%np) )
            allocate( RTI%snoise1(dat(1)%np) )
            allocate( RTI%sigma_n0(dat(1)%np) )
            allocate( RTI%sigma_n1(dat(1)%np) )
            allocate( RTI%n0_min(dat(1)%np) )
            allocate( RTI%n0_max(dat(1)%np) )
            allocate( RTI%n1_min(dat(1)%np) )
            allocate( RTI%n1_max(dat(1)%np) )
        endif

        do i = 1, dat(1)%np
            RTI%snoise0(i) = set%bn0_min + unirand(RTI%randcount)*(set%bn0_max-set%bn0_min)
            RTI%snoise1(i) = set%bn1_min + unirand(RTI%randcount)*(set%bn1_max-set%bn1_min)
            RTI%sigma_n0(i) = set%sigma_bn0
            RTI%sigma_n1(i) = set%sigma_bn1
            RTI%n0_min(i) = set%bn0_min
            RTI%n0_max(i) = set%bn0_max
            RTI%n1_min(i) = set%bn1_min
            RTI%n1_max(i) = set%bn1_max
            !RTI%snoise0(i) = (set%sn0_min + set%sn0_max)/2
            !RTI%snoise1(i) = (set%sn1_min + set%sn1_max)/2
        enddo

        ! if using phase velocity line for surface wave
        if(set%datatype==2)then
            do i = 1, dat(1)%np
                RTI%snoise0(i) = set%sn0_min + unirand(RTI%randcount)*(set%sn0_max-set%sn0_min)
                RTI%snoise1(i) = set%sn1_min + unirand(RTI%randcount)*(set%sn1_max-set%sn1_min)
                RTI%snoise1(i) = 0
                RTI%sigma_n0(i) = set%sigma_sn0
                RTI%sigma_n1(i) = set%sigma_sn1
                RTI%n0_min(i) = set%sn0_min
                RTI%n0_max(i) = set%sn0_max
                RTI%n1_min(i) = set%sn1_min
                RTI%n1_max(i) = set%sn1_max
                !RTI%snoise0(i) = (set%sn0_min + set%sn0_max)/2
                !RTI%snoise1(i) = (set%sn1_min + set%sn1_max)/2
            enddo
        endif


        ! if use differential time
        if(set%datatype==1 .or. set%datatype==3)then
            do i = dat(1)%np+1, 2*dat(1)%np
                RTI%snoise0(i) = set%dn0_min + unirand(RTI%randcount)*(set%dn0_max-set%dn0_min)
                RTI%snoise1(i) = set%dn1_min + unirand(RTI%randcount)*(set%dn1_max-set%dn1_min)
                RTI%sigma_n0(i) = set%sigma_dn0
                RTI%sigma_n1(i) = set%sigma_dn1
                RTI%n0_min(i) = set%dn0_min
                RTI%n0_max(i) = set%dn0_max
                RTI%n1_min(i) = set%dn1_min
                RTI%n1_max(i) = set%dn1_max
                !RTI%snoise0(i) = (set%sn0_min + set%sn0_max)/2
                !RTI%snoise1(i) = (set%sn1_min + set%sn1_max)/2
            enddo
        endif

        ! if jointly inversion
        if(set%datatype == 3)then
            do i = 2*dat(1)%np+1,dat(1)%np*2+dat(2)%np 
                RTI%snoise0(i) = set%sn0_min + unirand(RTI%randcount)*(set%sn0_max-set%sn0_min)
                RTI%snoise1(i) = set%sn1_min + unirand(RTI%randcount)*(set%sn1_max-set%sn1_min)
                RTI%sigma_n0(i) = set%sigma_sn0
                RTI%sigma_n1(i) = set%sigma_sn1
                RTI%n0_min(i) = set%sn0_min
                RTI%n0_max(i) = set%sn0_max
                RTI%n1_min(i) = set%sn1_min
                RTI%n1_max(i) = set%sn1_max
                !RTI%snoise0(i) = (set%sn0_min + set%sn0_max)/2
                !RTI%snoise1(i) = (set%sn1_min + set%sn1_max)/2
            enddo
        endif


        ! allocate and initialise the mean and std model
        allocate(RTI%aveP(set%grid%ny,set%grid%nx))
        allocate(RTI%stdP(set%grid%ny,set%grid%nx))
        allocate(RTI%aveS(set%grid%ny,set%grid%nx))
        allocate(RTI%stdS(set%grid%ny,set%grid%nx))
        RTI%aveP = 0
        RTI%stdP = 0
        RTI%aveS = 0
        RTI%stdS = 0

        ! sites id for each grid point
        allocate(RTI%sites_id(set%grid%ny,set%grid%nx))
        RTI%sites_id = 0

        ! tempering
        RTI%number_of_chains = set%number_of_temperatures
        RTI%number_of_temperatures = set%number_of_temperatures
        RTI%accepted_tempering = 0
        RTI%total_tempering = 0
        if(set%tempering==1)then
            allocate(RTI%temperature_values0(set%number_of_temperatures))
            allocate(RTI%temperature_indices0(set%number_of_temperatures))
            allocate(RTI%temperature_values(set%number_of_temperatures))
            allocate(RTI%temperature_indices(set%number_of_temperatures))
            RTI%temperature_indices0 = (/(i,i=1, set%number_of_temperatures)/)
            RTI%temperature_indices = RTI%temperature_indices0
            RTI%temperature_values0 = 1.0
            RTI%temperature_values = 1.0
        else
            ! for safe, allocate 1 temperature
            allocate(RTI%temperature_values(1))
            allocate(RTI%temperature_indices(1))
            RTI%temperature_indices = 1
            RTI%temperature_values = 1.0
        endif

        return
    end subroutine
        
    subroutine init_sample(sample)
        implicit none
        type(T_SAMPLE), intent(inout) :: sample

        sample%step = 0
        sample%accepted = .false.
        sample%vindex = 0
        sample%ncells = 0
        sample%misfit = 0
        sample%unweighted_misfit = 0
        sample%like = 0
        sample%coord = 0
        sample%values = 0
        sample%noise0 = 0
        sample%noise1 = 0
    end subroutine init_sample

    subroutine init_samples(samples)
        implicit none
        type(T_SAMPLE), dimension(:), intent(inout) :: samples
        integer  i

        do i = 1, size(samples)
            call init_sample(samples(i))
        enddo
    end subroutine init_samples

    subroutine read_info(RTI,filename)
        implicit none
        type(T_RUN_INFO), intent(inout) :: RTI
        character(len=*), intent(in) :: filename
        integer iunit
        open(unit=newunit(iunit),file=filename,status='old',action='read')
        read(iunit,*) RTI%randcount, RTI%like_count
        read(iunit,*) RTI%nsampled, RTI%nsamples
        read(iunit,*) RTI%nthin, RTI%sampletotal
        read(iunit,*) RTI%samplecount
        read(iunit,*) RTI%acceptedcount
        read(iunit,*) RTI%snoise0
        read(iunit,*) RTI%snoise1
        read(iunit,*) RTI%accepted_tempering, RTI%total_tempering
        read(iunit,*) RTI%temperature_indices
        read(iunit,*) RTI%temperature_values
        read(iunit,*) RTI%num_bad_model, RTI%num_bad_ray
        close(iunit)
        return
    end subroutine

    subroutine write_info(RTI,filename)
        implicit none
        type(T_RUN_INFO), intent(in) :: RTI
        character(len=*), intent(in) :: filename
        integer iunit
        open(unit=newunit(iunit),file=filename,status='replace',action='write')
        write(iunit,*) RTI%randcount, RTI%like_count
        write(iunit,*) RTI%nsampled, RTI%nsamples
        write(iunit,*) RTI%nthin, RTI%sampletotal
        write(iunit,*) RTI%samplecount
        write(iunit,*) RTI%acceptedcount
        write(iunit,*) RTI%snoise0
        write(iunit,*) RTI%snoise1
        write(iunit,*) RTI%accepted_tempering, RTI%total_tempering
        write(iunit,*) RTI%temperature_indices
        write(iunit,*) RTI%temperature_values
        write(iunit,*) RTI%num_bad_model, RTI%num_bad_ray
        close(iunit)
        return
    end subroutine

    subroutine read_samples(samples,filename)
    	character( len=* ), intent(in) :: filename
	    type(T_SAMPLE), dimension(:), intent(inout) :: samples

	    integer iunit, i
        integer stat
	    logical lexist

	    inquire(file=filename,exist=lexist)
	    if(.not.lexist) call exception_raiseError('sample file does not exist!')
	    open(unit=newunit(iunit),file=filename,status='old',access='stream',action='read')
	    !read(iunit) samples
            do i = 1, size(samples)
                read(iunit,iostat=stat) samples(i)%step, samples(i)%accepted,&
                samples(i)%vindex, samples(i)%ncells, samples(i)%misfit,&
                samples(i)%unweighted_misfit,samples(i)%like,  &
                samples(i)%coord, samples(i)%values, samples(i)%noise0, &
                samples(i)%noise1
                !read(iunit) samples(i)%noise0, samples(i)%noise1
                !read(iunit,iostat=stat) samples(i)
                if(stat /= 0)then
                    call exception_raiseError('Read total samples:&
                    & '//itoa(i-1))
                endif
            enddo
	    close(iunit)

    end subroutine

    subroutine write_samples(filename,samples)
        use m_utils, only : write_resume_unit
        use m_exception, only : exception_raiseWarning
        implicit none
    	character( len=* ), intent(in) :: filename
	    type(T_SAMPLE), dimension(:), intent(in) :: samples

	    integer i
	    integer nsamples
	    logical lexist

	    ! > open file for writing
	    inquire(file=filename,exist=lexist)
	    if(.not.lexist)then
	    	!call exception_raiseWarning('Sample file does not exist!')
	    	!call exception_raiseWarning('Create a new sample file!')
	    	open(unit=write_resume_unit,file=filename,status='new',access='stream')
	    else
	    	open(unit=write_resume_unit,file=filename,status='old',access='stream',position='append')
	    endif

	    ! > write samples to the file sample by sample
	    nsamples = size(samples)
	    !write(write_resume_unit) nsamples
	    do i = 1, nsamples
            write(write_resume_unit) samples(i)%step, samples(i)%accepted,&
            samples(i)%vindex, samples(i)%ncells, samples(i)%misfit,&
            samples(i)%unweighted_misfit,samples(i)%like,  &
            samples(i)%coord, samples(i)%values,&
            samples(i)%noise0, samples(i)%noise1
            !write(write_resume_unit) samples(i)%noise0, samples(i)%noise1
            !write(write_resume_unit) samples(i)
	    enddo
	    close(write_resume_unit)

    end subroutine

end module
