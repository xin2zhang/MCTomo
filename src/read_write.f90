module read_write

    use m_utils, only       : newunit, ii10, itoa,  write_doubles, &
                              integer_format, double_format, EPS, ini_dir, FILEDIR_SEPARATOR
    use m_exception, only   : exception_raiseError, exception_raiseWarning
    use run_info, only      : T_RUN_INFO, T_SAMPLE
    use m_likelihood, only  : T_LIKE
    use like_settings, only : T_LIKE_SET, T_LIKE_BASE, T_DATA

    use iso_c_binding
    implicit none

    private
    public :: read_freqs
    public :: read_vertices, write_vertices, read_vertex, read_txt_vertices
    public :: read_init_sigmas
    public :: read_mean, read_var
    public :: write_mean, write_var
    public :: write_initial_sample
    public :: write_likelihood, write_sigma
    public :: write_number_of_cells
    public :: write_temperatures
    public :: read_temperatures
    public :: read_temperatures0
    public :: write_data, write_atimes
    public :: write_src_samples


contains

    subroutine read_freqs(filename, freqs, ftype)
        implicit none
        character(len=*), intent(in) :: filename
        real(kind=ii10), dimension(:), intent(inout) :: freqs
        integer, intent(in), optional :: ftype

        !integer i
        integer iunit, iostatus
        integer np
        
        open(unit=newunit(iunit), file = filename, status = 'old', iostat = iostatus )
        if( iostatus /= 0 )then
            call exception_raiseError( 'error when open the file: ' // filename )
        endif

        read(iunit,*) np
        read(iunit,*) freqs(1:np)
        if(present(ftype)) then
            select case(ftype)
            case(0)
                freqs = freqs
            case(1)
                freqs = 1/freqs
            case default
                call exception_raiseError('Wrong file type for&
                &frequency/period')
            end select
        endif

        return
    end subroutine

    subroutine read_vertex(points,parameters,filename)
        use cgal_delaunay, only : p3, d3
        use m_utils, only : newunit
        implicit none
        type(d3), dimension(:), allocatable, intent(out) :: points
        type(p3), dimension(:), allocatable, intent(out) :: parameters
        character(len=*), intent(in) :: filename

        integer iunit
        integer(c_size_t) ncells, i

        open(unit=newunit(iunit), file=filename, access='stream',&
            status='old',action='read')
        read(iunit) ncells
        allocate( points(ncells) )
        allocate( parameters(ncells) )
        do i = 1, ncells
            read(iunit) points(i),parameters(i)
        enddo
        close(iunit)

        return
    end subroutine

    subroutine write_initial_sample(RTI,processor)
        !use m_utils, only: write_resume_unit,double_format
        implicit none
        type(T_RUN_INFO), intent(inout) :: RTI
        integer, intent(in) :: processor

        integer(c_size_t) i
        integer iunit

        !call cgal_delaunay_write(delaunay_ptr,'InitialDT_'//itoa(processor)//'.out')

        !call write_vertices(RTI,'InitialSample_'//itoa(processor))
        open(unit=newunit(iunit),file=trim(ini_dir)//FILEDIR_SEPARATOR//'InitialSample_'//itoa(processor)&
            //'.dat',status='new',access='stream')
        write(iunit) RTI%ncells
        do i = 1, RTI%ncells
            write(iunit) RTI%points(:,i), RTI%parameters(:,i)
        enddo
        close(iunit)
        
        ! write a txt version for man view
        open(unit=newunit(iunit),file=trim(ini_dir)//FILEDIR_SEPARATOR//'InitialSample_'//itoa(processor)&
            //'.txt',status='new', action='write')
        write(iunit,integer_format(1)) RTI%ncells
        do i = 1, RTI%ncells
            write(iunit,double_format(6)) RTI%points(:,i), RTI%parameters(:,i)
        enddo
        close(iunit)

        ! write initial sigma
        open(unit=newunit(iunit),file=trim(ini_dir)//FILEDIR_SEPARATOR//'InitialSigma_'//itoa(processor)&
            //'.dat',status='new',access='sequential')
        write(iunit,*) size(RTI%bnoise0), size(RTI%snoise0)
        do i = 1, size(RTI%bnoise0)
            write(iunit,*) RTI%bnoise0(i), RTI%bnoise1(i)
        enddo
        do i = 1, size(RTI%snoise0)
            write(iunit,*) RTI%snoise0(i), RTI%snoise1(i)
        enddo
        close(iunit)
        return

    end subroutine

    subroutine read_init_sigmas(RTI,filename)
        use m_utils, only : double_format, integer_format
        implicit none
        type(T_RUN_INFO), intent(inout) :: RTI
        character(len=*), intent(in) :: filename
        integer i, nsigma1, nsigma2, iunit
        open(unit=newunit(iunit),file=filename,status='old',access='sequential',action='read')
        read(iunit,*) nsigma1, nsigma2
        do i = 1, nsigma1
            read(iunit,*) RTI%bnoise0(i), RTI%bnoise1(i)
        enddo
        do i = 1, nsigma2
            read(iunit,*) RTI%snoise0(i), RTI%snoise1(i)
        enddo
        close(iunit)
        return
    end subroutine


    subroutine read_txt_vertices(RTI,filename)
        use m_utils, only : double_format, integer_format
        implicit none
        type(T_RUN_INFO), intent(inout) :: RTI
        character(len=*), intent(in) :: filename
        integer iunit
        integer(c_size_t) i
        open(unit=newunit(iunit),file=filename,status='old',access='sequential',action='read')
        read(iunit) RTI%ncells
        do i = 1, RTI%ncells
            read(iunit) RTI%points(:,i), RTI%parameters(:,i)
        enddo
        close(iunit)
        return
    end subroutine

    subroutine read_vertices(RTI,filename)
        use m_utils, only : double_format, integer_format
        implicit none
        type(T_RUN_INFO), intent(inout) :: RTI
        character(len=*), intent(in) :: filename
        integer iunit
        integer(c_size_t) i
        open(unit=newunit(iunit),file=filename,status='old',access='stream',action='read')
        read(iunit) RTI%ncells
        do i = 1, RTI%ncells
            read(iunit) RTI%points(:,i), RTI%parameters(:,i)
        enddo
        close(iunit)
        return
    end subroutine

    subroutine write_vertices(RTI,filename)
        use m_utils, only : double_format, integer_format
        implicit none
        type(T_RUN_INFO), intent(in) :: RTI
        character(len=*), intent(in) :: filename
        integer iunit
        integer(c_size_t) i
        open(unit=newunit(iunit),file=filename//'.dat',status='replace',&
            access='stream',action='write')
        write(iunit) RTI%ncells
        do i = 1, RTI%ncells
            write(iunit) RTI%points(:,i), RTI%parameters(:,i)
        enddo
        close(iunit)

        open(unit=newunit(iunit),file=filename//'.txt',status='replace',&
            action='write')
        write(iunit,integer_format(1)) RTI%ncells
        do i = 1, RTI%ncells
            write(iunit,double_format(6)) RTI%points(:,i), RTI%parameters(:,i)
        enddo
        close(iunit)
        return
    end subroutine

    subroutine read_mean(RTI,filename)
        use m_utils, only : read_resume_unit, read_doubles
        implicit none
        type(T_RUN_INFO), intent(inout) :: RTI
        character(len=*), intent(in) :: filename
        !character(len=:), allocatable :: str
        !real(kind=ii10), dimension(size(RTI%aveS,1),size(RTI%aveS,2),size(RTI%aveS,3)) :: aveS

        open(unit=read_resume_unit,file=filename,status='old',action='read')
        call read_doubles(RTI%aveP,' ')
        call read_doubles(RTI%aveS)
        close(read_resume_unit)
        RTI%aveS = RTI%aveS*RTI%nthin
        RTI%aveP = RTI%aveP*RTI%nthin
    end subroutine

    subroutine write_mean(RTI,filename)
        use m_utils, only : write_resume_unit, write_doubles
        implicit none
        type(T_RUN_INFO), intent(in) :: RTI
        character(len=*), intent(in) :: filename
        character(len=:), allocatable :: str
        real(kind=ii10), dimension(:,:,:), allocatable :: aveS, aveP

        allocate( aveP(size(RTI%aveP,1),size(RTI%aveP,2),size(RTI%aveP,3)) )
        allocate( aveS(size(RTI%aveS,1),size(RTI%aveS,2),size(RTI%aveS,3)) )
        if(RTI%nthin == 0)then
            aveS = 0.0
            aveP = 0.0
        else
            aveS = RTI%aveS/RTI%nthin
            aveP = RTI%aveP/RTI%nthin
        endif
        open(unit=write_resume_unit,file=filename,status='replace',action='write')
        str = itoa(size(RTI%aveS,1))//' '//itoa(size(RTI%aveS,2))&
            &//' '//itoa(2*size(RTI%aveS,3))
        call write_doubles(aveP,str)
        call write_doubles(aveS)
        close(write_resume_unit)
        return
    end subroutine

    subroutine read_var(RTI,filename)
        use m_utils, only : read_resume_unit, read_doubles
        implicit none
        type(T_RUN_INFO), intent(inout) :: RTI
        character(len=*), intent(in) :: filename
        !character(len=:), allocatable :: str
        real(kind=ii10), dimension(:,:,:), allocatable :: stdS

        allocate( stdS(size(RTI%stdS,1),size(RTI%stdS,2),size(RTI%stdS,3)) )

        open(unit=read_resume_unit,file=filename,status='old',action='read')
        call read_doubles(RTI%stdP,' ')
        call read_doubles(stdS)
        close(read_resume_unit)
        RTI%stdP =  RTI%stdP*RTI%nthin
        RTI%stdS = stdS*RTI%nthin
    end subroutine

    subroutine write_var(RTI,filename)
        use m_utils, only : write_resume_unit, write_doubles
        implicit none
        type(T_RUN_INFO), intent(in) :: RTI
        character(len=*), intent(in) :: filename
        character(len=:), allocatable :: str
        real(kind=ii10), dimension(:,:,:), allocatable :: stdP, stdS

        allocate( stdS(size(RTI%stdS,1),size(RTI%stdS,2),size(RTI%stdS,3)) )
        allocate( stdP(size(RTI%stdP,1),size(RTI%stdP,2),size(RTI%stdP,3)) )

        if( RTI%nthin == 0)then
            stdS = 0.0
            stdP = 0.0
        else
            stdS = RTI%stdS/RTI%nthin
            stdP = RTI%stdP/RTI%nthin
        endif
        open(unit=write_resume_unit,file=filename,status='replace',action='write')
        str = itoa(size(RTI%stdS,1))//' '//itoa(size(RTI%stdS,2))&
            &//' '//itoa(2*size(RTI%stdS,3))
        call write_doubles(stdP,str)
        call write_doubles(stdS)
        close(write_resume_unit)
        return
    end subroutine


    subroutine write_likelihood(filename,samples,initial_like)
        use m_utils, only : write_resume_unit, write_double
        implicit none
        character( len=* ), intent(in) :: filename
        type(T_SAMPLE), dimension(:), intent(in) :: samples
        type(T_LIKE), intent(in), optional :: initial_like

        integer i
        integer nsamples
        logical lexist

        if(size(samples)==0) return

        ! > open file for writing
        inquire(file=filename,exist=lexist)
        if(.not.lexist)then
            open(unit=write_resume_unit,file=filename,status='new',access='sequential')
        else
            open(unit=write_resume_unit,file=filename,status='old',access='sequential',position='append')
        endif

        ! > write samples to the file sample by sample
        if(present(initial_like))&
            call write_doubles([initial_like%like,initial_like%misfit,initial_like%unweighted_misfit])
        nsamples = size(samples)
        do i = 1, nsamples
                call write_doubles([samples(i)%like, samples(i)%misfit, samples(i)%unweighted_misfit])
        enddo
        close(write_resume_unit)

    end subroutine

    subroutine write_number_of_cells(filename,samples,initial_number)
        use m_utils, only : write_resume_unit, integer_format
        implicit none
        character( len=* ), intent(in) :: filename
        type(T_SAMPLE), dimension(:), intent(in) :: samples
        integer(c_size_t), intent(in), optional  :: initial_number

        integer i
        integer nsamples
        logical lexist

        if(size(samples)==0) return

        ! > open file for writing
        inquire(file=filename,exist=lexist)
        if(.not.lexist)then
            open(unit=write_resume_unit,file=filename,status='new',access='sequential')
        else
            open(unit=write_resume_unit,file=filename,status='old',access='sequential',position='append')
        endif

        ! > write samples to the file sample by sample
            if(present(initial_number))&
                write(write_resume_unit,integer_format(1)) initial_number
        nsamples = size(samples)
        do i = 1, nsamples
                write(write_resume_unit,integer_format(1)) samples(i)%ncells
        enddo
        close(write_resume_unit)

    end subroutine

    subroutine write_sigma(filename,samples)
        use m_utils, only : write_resume_unit, write_double
        implicit none
        character( len=* ), intent(in) :: filename
        type(T_SAMPLE), dimension(:), intent(in) :: samples

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
            !if(present(initial_like))&
            !    call write_doubles([initial_like%like,initial_like%misfit])
        nsamples = size(samples)
        do i = 1, nsamples
                !call write_doubles(samples(i)%noise0)
                !call write_doubles(samples(i)%noise1)
                write(write_resume_unit,double_format(1))&
                sqrt(samples(i)%misfit/samples(i)%like/2)
        enddo
        close(write_resume_unit)

    end subroutine

    subroutine read_temperatures0(filename,values,begin)
        implicit none
        character(len=*), intent(in) :: filename
        real(kind=ii10), dimension(:), intent(inout) :: values
        integer, intent(in) :: begin

        !integer i
        integer iunit, iostatus
        
        open(unit=newunit(iunit), file = filename, status = 'old', iostat = iostatus )
        if( iostatus /= 0 )then
            call exception_raiseError( 'error when open the file: ' // filename )
        endif

        read(iunit,*) values(begin:ubound(values,1))
        !do i = begin, size(indices)
        !    read(iunit,*) indices(i), values(i)
        !enddo

    end subroutine

    subroutine write_temperatures(filename,temperatures)
        use m_utils, only : write_resume_unit, write_doubles
        implicit none
        character( len=* ), intent(in) :: filename
        real(kind=ii10), dimension(:,:), intent(in) :: temperatures

        integer i, nd
        integer nsamples
        logical lexist

        if(size(temperatures,2)==0) return

        ! > open file for writing
        inquire(file=filename,exist=lexist)
        if(.not.lexist)then
            open(unit=write_resume_unit,file=filename,status='new',access='sequential')
        else
            open(unit=write_resume_unit,file=filename,status='old',access='sequential',position='append')
        endif

        nd = size(temperatures,1)
        nsamples = size(temperatures,2)
        do i = 1, nsamples
            write(write_resume_unit,double_format(nd)) temperatures(:,i)
        enddo
        close(write_resume_unit)

    end subroutine

    subroutine read_temperatures(filename,temperatures)
        use m_utils, only : write_resume_unit
        implicit none
        character( len=* ), intent(in) :: filename
        real(kind=ii10), dimension(:,:), intent(out) :: temperatures

        integer i, nd
        integer nsamples

        ! > open file for writing
        open(unit=write_resume_unit,file=filename,status='old',action='read')

        nd = size(temperatures,1)
        nsamples = size(temperatures,2)
        do i = 1, nsamples
            read(write_resume_unit,double_format(nd)) temperatures(:,i)
        enddo
        close(write_resume_unit)

    end subroutine

    subroutine write_data(filename,dat,RTI,likelihoods,like_set)
        implicit none
        character( len=* ), intent(in) :: filename
        type(T_DATA), dimension(:), intent(in) :: dat
        type(T_RUN_INFO), intent(in) :: RTI
        type(T_LIKE_BASE), dimension(:), intent(in) :: likelihoods
        type(T_LIKE_SET), intent(in) :: like_set

        select case (like_set%datatype)
        case(0,1)
            ! write data for body waves
            call write_atimes(filename//'_body.dat',dat(1),RTI,likelihoods(1))
        case(2)
            ! write data for surface waves
            call write_atimes(filename//'_surf.dat',dat(1),RTI,likelihoods(1))
        case(3)
            ! write data for body waves
            call write_atimes(filename//'_body.dat',dat(1),RTI,likelihoods(1))
            ! write data for surface waves
            call write_atimes(filename//'_surf.dat',dat(2),RTI,likelihoods(2))
        end select
    end subroutine write_data

    subroutine write_atimes(filename,dat,RTI,like)
        implicit none
        character( len=* ), intent(in)  :: filename
        type(T_DATA), intent(in)        :: dat
        type(T_RUN_INFO), intent(in)    :: RTI
        type(T_LIKE_BASE), intent(in)   :: like

        ! local variable
        integer, parameter :: yesrays = 1, norays =  0
        integer i, j, k !iterator
        integer npairs
        integer nsrc, nrev
        integer iunit1, iostatus

        open(unit=newunit(iunit1), file = filename, status = 'unknown', action='write', iostat = iostatus )
        if( iostatus /= 0 )then
            call exception_raiseError( 'error when open the file: ' // filename )
        endif

        write(iunit1,*) dat%np
        write(iunit1,*) dat%freqs

        nsrc = dat%nsrc
        nrev = dat%nrev

        npairs = 0
        do i = 1, nsrc
           do j = 1, nrev
              npairs = npairs + 1
              if( any(dat%raystat(npairs,1,:)==1) ) then
                  write(iunit1,*) yesrays
              else
                  write(iunit1,*) norays
                  cycle
              endif

              do k = 1, dat%np
                  if( dat%raystat(npairs,1,k)==1 ) then
                      write(iunit1,*) like%phaseTime(j,i,k)+RTI%locations(4,i), like%sigma(npairs,k)
                  else
                      write(iunit1,*) like%phaseTime(j,i,k), like%sigma(npairs,k)
                      !write(iunit1,*) dat%ttime(npairs,1,k),dat%ttime(npairs,2,k)
                  endif
              enddo
           enddo    
        enddo

        close(iunit1)

    end subroutine

    subroutine write_src_samples(samples,filename)
        use m_utils, only : write_resume_unit, write_double
        implicit none
        character( len=* ), intent(in) :: filename
        real(kind=ii10), dimension(:,:,:), intent(in) :: samples

        integer i
        integer nsamples
        logical lexist

        if(size(samples,3)==0) return

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
        write(write_resume_unit) samples
        close(write_resume_unit)

    end subroutine

end module
