module like_settings

    use iso_c_binding, only : c_int
    use m_exception,    only : exception_raiseError
    use m_utils, only : ii10, STRLEN, itoa, rtoa, EPS, newunit, float_format
    use m_fm2d, only : T_RAY
    use m_settings, only : T_GRID, out_bnd

    implicit none

    private

    public :: T_DATA
    public :: read_data, read_sources, write_sources, write_times
    public :: T_LIKE_BASE, likeBaseSetup
    public :: T_LIKE_SET
    public :: write_phase

    ! data 
    type T_DATA
        integer  :: datatype
        integer  :: nrays, np, ndtime
        integer  :: nsrc, nrev
        real( kind=ii10 ), dimension(:,:), allocatable   :: src
        real( kind=ii10 ), dimension(:,:), allocatable   :: rev
        integer, dimension(:,:,:), allocatable           :: raystat
        real( kind=ii10 ), dimension(:,:,:), allocatable :: ttime
        real( kind=ii10 ), dimension(:,:,:), allocatable :: dtime
        real( kind=ii10 ), dimension(:,:,:), allocatable :: dcc
        ! if body waves, for different types; if surface waves, for freqs
        real( kind=ii10 ), dimension(:), allocatable     :: freqs
    end type T_DATA

    ! likelihood
    type T_LIKE_BASE
        real( kind=ii10 )       :: like, misfit, unweighted_misfit
        ! vel and gvel for surface wave phase and group velocities
        !real( kind=ii10 ), dimension(:,:,:), allocatable :: vel, gvel
        real(kind=ii10), dimension(:,:,:), allocatable :: phaseTime
        real(kind=ii10), dimension(:,:,:), allocatable :: dtime, dsigma
        real(kind=ii10), dimension(:), allocatable :: grads
        real(kind=ii10), dimension(:,:,:,:), allocatable :: field4d
        real(kind=ii10), dimension(:,:,:,:,:), allocatable :: field5d
        real(kind=ii10), dimension(:,:), allocatable :: btime
        type(T_RAY), dimension(:,:), allocatable :: rays ! currently, only for straight rays
        real( kind=ii10 ), dimension(:,:), allocatable :: srdist, sigma
        logical :: straightRaySet
    end type T_LIKE_BASE

    ! likelihood settings
    type T_LIKE_SET
        ! data type, 0 for P wave, 1 for P&S wave, 2 surface wave, 3 for P, S and surface waves
        integer :: datatype
        ! source files
        character(len=STRLEN) :: bsources_file, ssources_file
        character(len=STRLEN) :: breceivers_file, sreceivers_file
        character(len=STRLEN) :: bdata_file, sdata_file
        character(len=STRLEN) :: dt_file
        ! nd: number of body wave data types; np: number of periods
        ! integer :: nd, np
        ! sigma dependence
        integer :: sigdep
        ! parallelized threads number
        integer :: nthreads
        ! grid info
        type(T_GRID) :: grid
        ! straight ray or not
        integer :: isStraight
        ! input for dispersion code
        ! minimu and maxmum phase velocity
        real(kind=ii10) :: dPhaseVel
        integer :: raylov, phaseGroup
        real(kind=ii10) :: tol
        ! input for fast marching code
        integer :: gridx
        integer :: gridy
        integer :: sgref
        integer :: sgdic
        integer :: sgext
        integer :: order
        real( kind=ii10 ) :: band
        integer :: uar
        integer :: dynamic
    endtype

contains

    subroutine likeBaseSetup(like,dat,set,ncell_max)
        implicit none
        type(T_LIKE_BASE), intent(out) :: like
        type(T_DATA), intent(in)  :: dat
        type(T_LIKE_SET), intent(in) :: set
        integer, intent(in) :: ncell_max

        integer nsrc, nrev
        integer dim_src, dim_rev
        integer ierr
        integer i, j, k, nrr, ndtime

        select case(set%datatype)
        case (0,1,3)
            nsrc = dat%nsrc
            nrev = dat%nrev
            dim_src=size(dat%src,1)
            dim_rev=size(dat%rev,1)
            if(dim_rev>dim_src) call exception_raiseError('dimension of source is less than receiver')

            allocate( like%phaseTime(nrev,nsrc,dat%np) )
            allocate( like%srdist(nrev*nsrc,dat%np) )
            allocate( like%sigma(nrev*nsrc,dat%np) )
            like%srdist = 1.0 ! safe
            like%sigma = 1.0 ! safe
            ! initialise source-receiver ray lenth using source-receiver distance
	        do i = 1, dat%np
	            nrr  = 0
	            do j = 1, nsrc
	   	            do k = 1, nrev
	                    nrr = nrr + 1
                        like%srdist(nrr,:) = norm2(dat%src(1:dim_rev,j)-dat%rev(:,k))
	                enddo
	            enddo
	        enddo
        case (2)
            allocate( like%phaseTime(1,dat%nrays,dat%np) )
            allocate( like%srdist(dat%nrays,dat%np) )
            allocate( like%sigma(dat%nrays,dat%np) )
            ! initialise source-receiver ray lenth using noise level from data
            like%srdist = 1.0 ! safe
            like%sigma = 1.0 ! safe
            like%srdist = dat%ttime(:,2,:)
        case default
            call exception_raiseError('Datatype not supported yet!')
        endselect

        like%like = huge(like%like)
        like%misfit = huge(like%misfit)
        like%unweighted_misfit = huge(like%unweighted_misfit)
        like%phaseTime = 0.0
        like%straightRaySet = .false.


        ! using differential time
        if( set%datatype == 1)then
            ndtime = (nsrc-1)*nsrc/2
            allocate( like%dtime(nrev,ndtime,dat%np) )
            allocate( like%dsigma(nrev,ndtime,dat%np) )
            like%dtime = 0
            like%dsigma = 1.0
        endif

        ! gradient
        if( dat%datatype == 0)then
            allocate( like%grads(2*ncell_max+4*dat%nsrc) )
            like%grads = 0
        else
            allocate( like%grads(2*ncell_max) )
            like%grads = 0
        endif

        if( set%dynamic == 1 )then
            select case (dat%datatype)
            case (0)
                allocate( like%field5d(set%grid%nz,set%grid%ny,set%grid%nx,dat%nrev,dat%np),stat=ierr )
                if(ierr /= 0)then
                    call exception_raiseError('Cannot allocate the memory for&
                        & travel time field! Switch off the dynamic option!')
                endif
                like%field5d = 0
                allocate( like%btime(nsrc,dat%np))
                like%btime = 0
            case(1)
                allocate(like%field4d(set%grid%ny,set%grid%nx,nsrc,dat%np),stat=ierr)
                if(ierr /= 0)then
                    call exception_raiseError('Cannot allocate the memory for&
                        & travel time field! Switch off the dynamic option!')
                endif
                like%field4d = 0
                allocate(like%btime(nsrc,dat%np))
                like%btime = 0
            endselect
        endif

        return
    end subroutine likeBaseSetup

    subroutine read_data(dat,like_set,bnd)
        implicit none
        type(T_DATA), dimension(:), intent(out) :: dat
        type(T_LIKE_SET), intent(in) :: like_set
        !> bnd boudary of latitude/longtitude, longmin longmax latmin latmax
        real( kind=ii10 ), dimension(:), intent(in), optional :: bnd

        select case (like_set%datatype)
        case(0,1)
            ! read data for body waves
            call read_sources(like_set%bsources_file,dat(1)%src,bnd)
            call read_sources(like_set%breceivers_file,dat(1)%rev,bnd)
            call read_times_3(like_set%bdata_file,dat(1),like_set%dt_file)
            dat(1)%datatype = 0
        case(2)
            ! read data for surface waves
            !call read_sources(like_set%ssources_file,dat(1)%src,bnd)
            !call read_sources(like_set%sreceivers_file,dat(1)%rev,bnd)
            call read_phase(like_set%sdata_file,dat(1))
            dat(1)%datatype = 1
        case(3)
            ! read data for body waves
            call read_sources(like_set%bsources_file,dat(1)%src,bnd)
            call read_sources(like_set%breceivers_file,dat(1)%rev,bnd)
            call read_times_3(like_set%bdata_file,dat(1),like_set%dt_file)
            dat(1)%datatype = 0
            ! read data for surface waves
            call read_sources(like_set%ssources_file,dat(2)%src,bnd)
            call read_sources(like_set%sreceivers_file,dat(2)%rev,bnd)
            call read_times_3(like_set%sdata_file,dat(2))
            dat(2)%datatype = 1
        end select
    end subroutine read_data

    subroutine read_sources(filename,src,bnd)
        implicit none
        character( len=* ), intent(in) :: filename
        real( kind=ii10 ), dimension(:,:), allocatable, intent(out) :: src
        !> bnd boudary of latitude/longtitude, longmin longmax latmin latmax
        real( kind=ii10 ), dimension(:), intent(in), optional :: bnd

        !> local variable
        integer iunit
        integer iostatus
        integer nsrc, ndim
        integer i ! iterator

        if(allocated(src)) deallocate(src)
        open(unit=newunit(iunit), file = filename, status = 'old', iostat = iostatus )
        if( iostatus /= 0 )then
            call exception_raiseError( 'error when open the file: ' // filename )
        else
            read(iunit,*) nsrc, ndim
            allocate( src(ndim,nsrc) )
            do i = 1, nsrc
                read(iunit,*) src(:,i)
                if( present(bnd) ) then
                    if( out_bnd(src(:,i),bnd) ) call exception_raiseError( 'Source ' // itoa(i) // ' is out of bounds!' )
                endif
            enddo
        endif
        close(iunit)

    end subroutine

    subroutine write_sources(src,filename)
        implicit none
        character( len=* ), intent(in) :: filename
        real( kind=ii10 ), dimension(:,:), intent(in) :: src

        !> local variable
        integer iunit
        integer iostatus
        integer nsrc, ndim
        integer i ! iterator

        ndim =  size(src,1)
        nsrc =  size(src,2)


        open(unit=newunit(iunit), file = filename, status = 'replace', iostat = iostatus, action='write' )
        if( iostatus /= 0 )then
            call exception_raiseError( 'error when open the file: ' // filename )
        else
            write(iunit,*) nsrc, ndim
            do i = 1, nsrc
                write(iunit,*) src(:,i)
	        enddo
	    endif
	    close(iunit)
	
    end subroutine

    subroutine read_times_2(filename,time,raystat,src,rev,dsfile)
    	implicit none
	    character( len=* ), intent(in) :: filename
	    real( kind=ii10 ), dimension(:,:), allocatable, intent(out) :: time
	    integer, dimension(:,:), allocatable, intent(out) :: raystat
	    real( kind=ii10 ), dimension(:,:), intent(in) :: src
	    real( kind=ii10 ), dimension(:,:), intent(in) :: rev
	    character( len=* ), intent(in), optional :: dsfile

	    ! local variable
	    integer i, j !iterator
	    integer nrays, npairs
	    integer nsrc, nrev
	    integer validity
	    integer iunit1, iunit2
	    integer iostatus
	    real( kind=ii10 ) t, n
	    real( kind=ii10 ) dssval
	    real( kind=ii10 ), dimension(:,:), allocatable :: time_src_rev

        open(unit=newunit(iunit1), file = filename, status = 'old', iostat = iostatus )
        if( iostatus /= 0 )then
            call exception_raiseError( 'error when open the file: ' // filename )
	    endif
	    if( present(dsfile) ) then
            open(unit=newunit(iunit2), file = filename, status = 'old', iostat = iostatus )
            if( iostatus /= 0 )then
                 call exception_raiseError( 'error when open the file: ' // dsfile )
	        endif
	    endif

	    nsrc = size(src,2)
	    nrev = size(rev,2)

	    allocate( raystat(nsrc*nrev,2) )
	    allocate( time_src_rev(nsrc*nrev,7) )

	    nrays = 0
	    npairs = 0
	    do i = 1, nsrc
	       do j = 1, nrev
	          npairs = npairs + 1
	          read(iunit1, *) validity, t, n
	          if( present(dsfile) ) then
	              read(iunit2,*) dssval
	          else
	          	  dssval = 1
	          endif

	          if(validity /= 0) then
	              nrays = nrays + 1
	    	      time_src_rev(nrays,1:2) = src(1:2,i)
	    	      time_src_rev(nrays,3:4) = rev(1:2,i)
	    	      time_src_rev(nrays,5) = t
	    	      time_src_rev(nrays,6) = 1
	    	      time_src_rev(nrays,7) = dssval
	    	      raystat(npairs,1) = 1
	    	      raystat(npairs,2) = nrays
	          else
	              raystat(npairs,:) = 0
	          endif
           enddo	
        enddo

	    close(iunit1)
	    if(present(dsfile)) close(iunit2)

	    allocate( time(nrays,7) )
	    time(:,:) = time_src_rev(1:nrays,:)
	    deallocate( time_src_rev )

    end subroutine

    subroutine read_times_3(filename,dat,dtfile,dsfile)
    	implicit none
	    character( len=* ), intent(in)  :: filename
	    character( len=* ), intent(in), optional  :: dtfile
        type(T_DATA), intent(inout)       :: dat
	    character( len=* ), intent(in), optional :: dsfile

	    ! local variable
	    integer i, j, k !iterator
        integer isrc1, isrc2, idt, irec, dtime, dcc, iphase, iswap
	    integer npairs, nrays, ndtime
	    integer nsrc, nrev
	    integer validity
        integer nfreqs
	    integer iunit1, iunit2
	    integer iostatus
	    real( kind=ii10 ) t, n
	    real( kind=ii10 ) dssval
        logical lexist
        integer stat
	    !real( kind=ii10 ), dimension(:,:,:), allocatable :: time_src_rev

        open(unit=newunit(iunit1), file = filename, status = 'old', iostat = iostatus )
        if( iostatus /= 0 )then
            call exception_raiseError( 'error when open the file: ' // filename )
	    endif
	    if( present(dsfile) ) then
            open(unit=newunit(iunit2), file = filename, status = 'old', iostat = iostatus )
            if( iostatus /= 0 )then
                 call exception_raiseError( 'error when open the file: ' // dsfile )
	        endif
	    endif

        read(iunit1,*) nfreqs
        allocate( dat%freqs(nfreqs) )
        read(iunit1,*) dat%freqs
        dat%np = nfreqs

	    nsrc = size(dat%src,2)
	    nrev = size(dat%rev,2)
        dat%nsrc = nsrc
        dat%nrev = nrev

	    allocate( dat%raystat(nsrc*nrev,2,nfreqs) )
	    allocate( dat%ttime(nsrc*nrev,3,nfreqs) )
	    !allocate( time_src_rev(nsrc*nrev,3,nfreqs) )
        dat%raystat =  0
        dat%ttime = 0
        !time_src_rev =  0

	    npairs = 0
        nrays = 0
	    do i = 1, nsrc
	       do j = 1, nrev
	          npairs = npairs + 1
	          read(iunit1, *) validity
	          if( present(dsfile) ) then
	              read(iunit2,*) dssval
	          else
	          	  dssval = 1
	          endif

	          if(validity /= 0) then
                  do k = 1, nfreqs
                      read(iunit1,*) t, n
	                  dat%ttime(npairs,1,k) = t
                      dat%ttime(npairs,2,k) = n
	         	      dat%ttime(npairs,3,k) = dssval
                      if( abs(t+1) > EPS ) then
	    	              dat%raystat(npairs,1,k) = 1
	    	              dat%raystat(npairs,2,k) = npairs
                          nrays = nrays + 1
                      endif
                  enddo
	          else
	              dat%raystat(npairs,:,:) = 0
                  dat%ttime(npairs,1,:) =  -1
                  dat%ttime(npairs,2,:) = 1
	          endif
           enddo	
        enddo

	    close(iunit1)
	    if(present(dsfile)) close(iunit2)

        dat%nrays = nrays
	    !allocate( dat%ttime(nrays,3,nfreqs) )
	    !dat%ttime(:,:,:) = time_src_rev(1:nrays,:,:)
	    !deallocate( time_src_rev )
        
        if(.not. present(dtfile)) return
        ! if dt file does not exist, just return
        inquire(file=dtfile,exist=lexist)
        if(.not.lexist)then
            ! initialize to zero
            ndtime = (nsrc-1)*nsrc/2
            allocate( dat%dtime(nrev,ndtime,dat%np) )
            allocate( dat%dcc(nrev,ndtime,dat%np) )
            dat%dtime = 0
            dat%dcc = 0
            dat%ndtime = 0
            return
        endif

        ! read dtimes
        open(unit=newunit(iunit1), file = dtfile, status = 'old', iostat = iostatus )
        if( iostatus /= 0 )then
            call exception_raiseError( 'error when open the file: ' // dtfile )
	    endif

        ndtime = (nsrc-1)*nsrc/2
        allocate( dat%dtime(nrev,ndtime,dat%np) )
        allocate( dat%dcc(nrev,ndtime,dat%np) )
        dat%dtime = 0
        dat%dcc = 0
        ndtime = 0
        do
            read(iunit1,*,iostat=stat) isrc1, isrc2, irec, dtime, dcc, iphase
            if(stat/=0) exit
            if(isrc1>isrc2)then
                iswap = isrc1
                isrc1 = isrc2
                isrc2 = iswap
                dtime = -dtime
            endif
            idt = (nsrc-1+nsrc-isrc1+1)*(isrc1-1)/2 + isrc2 - isrc1
            dat%dtime(irec,idt,iphase) = dtime
            dat%dcc(irec,idt,iphase) = dcc
            ndtime = ndtime + 1
        enddo

        dat%ndtime = ndtime

        close(iunit1)

    end subroutine

    subroutine write_times(filename,dat,like)
    	implicit none
	    character( len=* ), intent(in)  :: filename
        type(T_DATA), intent(in)        :: dat
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
                      write(iunit1,*) like%phaseTime(j,i,k), like%sigma(npairs,k)
                  else
                      write(iunit1,*) like%phaseTime(j,i,k), like%sigma(npairs,k)
                      !write(iunit1,*) dat%ttime(npairs,1,k),dat%ttime(npairs,2,k)
                  endif
              enddo
           enddo	
        enddo

	    close(iunit1)

    end subroutine

    subroutine read_phase(filename,dat)
    	implicit none
	    character( len=* ), intent(in)  :: filename
        type(T_DATA), intent(inout)       :: dat
        !character( len=* ), intent(in), optional :: dsfile

	    ! local variable
	    integer i, j !iterator
	    integer ndata
        integer nfreqs
	    integer iunit1
	    integer iostatus
        logical lexist
        integer stat

        open(unit=newunit(iunit1), file = filename, status = 'old', iostat = iostatus )
        if( iostatus /= 0 )then
            call exception_raiseError( 'error when open the file: ' // filename )
	    endif

        read(iunit1,*) ndata, nfreqs
        allocate( dat%freqs(nfreqs) )
        read(iunit1,*) dat%freqs
        dat%np = nfreqs

	    allocate( dat%raystat(ndata,2,nfreqs) )
	    allocate( dat%ttime(ndata,3,nfreqs) )
        dat%raystat =  0
        dat%ttime = 0

        ! read phase velocity
	    do i = 1, ndata
	        read(iunit1, *) dat%ttime(i,1,:)
            where( abs(dat%ttime(i,1,:)+1) > EPS )
	    	    dat%raystat(i,1,:) = 1
	    	    dat%raystat(i,2,:) = i
            endwhere
        enddo

        ! read phase velocity noise
	    do i = 1, ndata
	        read(iunit1, *) dat%ttime(i,2,:)
        enddo
	    close(iunit1)

        dat%nrays = ndata 
        
    end subroutine

    subroutine write_phase(filename,dat,like)
    	implicit none
	    character( len=* ), intent(in)  :: filename
        type(T_DATA), intent(in)        :: dat
        type(T_LIKE_BASE), intent(in)   :: like
        !character( len=* ), intent(in), optional :: dsfile

	    ! local variable
	    integer i, j !iterator
	    integer ndata
        integer nfreqs
	    integer iunit1
	    integer iostatus
        logical lexist
        integer stat

        open(unit=newunit(iunit1), file = filename, iostat = iostatus )
        if( iostatus /= 0 )then
            call exception_raiseError( 'error when open the file: ' // filename )
	    endif

        ndata = dat%nrays
        nfreqs = dat%np
        write(iunit1,*) ndata, nfreqs
        write(iunit1,*) dat%freqs

        ! read phase velocity
	    do i = 1, ndata
	        write(iunit1, float_format(nfreqs)) like%phaseTime(1,i,:)
        enddo

        ! read phase velocity noise
	    do i = 1, ndata
	        write(iunit1, float_format(nfreqs)) like%sigma(i,:)
        enddo
	    close(iunit1)
        
    end subroutine

end module like_settings
