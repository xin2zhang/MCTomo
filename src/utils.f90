!> modules for some frequently used functions and subroutines
module m_utils
    use iso_c_binding
    implicit none

    private

    public :: itoa
    public :: rtoa
    public :: ishuge
    public :: delete_file
    public :: newunit
    public :: bisearch
    public :: integer_format
    public :: double_format
    public :: write_double
    public :: write_doubles
    public :: read_doubles
    public :: brute_search
    public :: vs2vp, vp2rho, vs2vp_3d, vp2rho_3d
    public :: create_dir
    public :: init_random_seed
    public :: init_random_seed_internal
    
    interface read_doubles
        module procedure read_doubles_1, read_doubles_2, read_doubles_3
    end interface read_doubles

    interface write_doubles
        module procedure write_doubles_1, write_doubles_2, write_doubles_3
    end interface write_doubles

    interface itoa
        module procedure itoa4, itoa8
    end interface itoa

    interface ishuge
        module procedure ishuge_4, ishuge_8, ishuge_16
    end interface ishuge

    interface
        function makedir(path) bind(C,name="makedir")
            import
            integer(c_int) :: makedir
            character(kind=c_char,len=1), dimension(*), intent(in) :: path
        endfunction
    endinterface

    !> real type
    integer, parameter, public :: ii5 = selected_real_kind(5,10)
    !integer, parameter, public :: ii10 = selected_real_kind(10,100)
    integer, parameter, public :: ii10 = c_double

    ! DEBUG
    logical, parameter, public :: DEBUG = .true.
    !> mathmatical const
    real(kind=ii10), parameter, public :: eps = 0.000000001
    real(kind=ii10), parameter, public :: pii = 3.1415926535898

    !> files/dir const
    character, parameter, public :: FILEDIR_SEPARATOR = '/'
    character, parameter, public :: NAME_SEPARATOR = '.'

    !> string length, directory length
    integer, parameter, public :: STRLEN = 100
    integer, parameter, public :: DIRLEN = 100
    character(len=STRLEN), parameter, public :: results_dir = 'Results'
    character(len=DIRLEN), parameter, public :: last_dir = 'Results'//FILEDIR_SEPARATOR//'last_run'
    character(len=DIRLEN), parameter, public :: resume_dir = 'Results'//FILEDIR_SEPARATOR//'resume'
    character(len=DIRLEN), parameter, public :: ini_dir = 'Results'//FILEDIR_SEPARATOR//'ini'
    character(len=STRLEN), parameter, public :: log_dir = 'Log'

    !> exit failure const
    integer, parameter, public :: exit_failure = -1

    integer, parameter :: delete_unit = 99

    !> The default writing formats
    integer, parameter :: fmt_len = 200
    character(7) :: DB_FMT='E17.8E3'
    character(5) :: FLT_FMT='F10.4'
    character(3) :: INT_FMT='I12'
    integer, parameter,public :: read_resume_unit = 11
    integer, parameter,public :: write_resume_unit = 12

contains

    !> experience equation for vs to vp and vp to density
    !  For different rock types, please use different laws
    function vs2vp(vs) result(vp)
        implicit none
        real(kind=c_double), intent(in) :: vs
        real(kind=c_double) :: vp

        real( kind=c_double ), parameter :: POISSON = 1.730
        ! crust
        vp = vs * POISSON
        ! sedimentary
        !vp = 1.16*vs + 1.36
    end function

    subroutine vs2vp_3d(vs,vp)
        implicit none
        real(kind=c_double), dimension(:,:,:), intent(in)  :: vs
        real(kind=c_double), dimension(:,:,:), intent(out) :: vp

        real( kind=c_double ), parameter :: POISSON = 1.730
        ! crust
        vp = vs * POISSON
        ! sedimentary
        !vp = 1.16*vs + 1.36
    end subroutine

    function vp2rho(vp) result(rho)
        implicit none
        real(kind=c_double), intent(in) :: vp
        real(kind=c_double) :: rho

        ! crust
        rho = 2.35 +  0.036*(vp-3)**2
        ! sedimentary
        !rho = 1.74*vp**0.25
    endfunction

    subroutine vp2rho_3d(vp, rho)
        implicit none
        real(kind=c_double), dimension(:,:,:), intent(in) :: vp
        real(kind=c_double), dimension(:,:,:), intent(out) :: rho

        ! crust
        rho = 2.35 +  0.036*(vp-3)**2
        ! sedimentary
        !rho = 1.74*vp**0.25
    end subroutine

    function itoa4(n) result(newstring)
        implicit none
        integer(kind=4), intent(in) :: n
        character(:), allocatable :: newstring ! fortran 2003
        character(len=range(n)+2) :: tmp

        write(tmp,'(I0)') n
        newstring = trim(tmp)
    end function

    function itoa8(n) result(newstring)
        implicit none
        integer(kind=8), intent(in) :: n
        character(:), allocatable :: newstring ! fortran 2003
        character(len=range(n)+2) :: tmp

        write(tmp,'(I0)') n
        newstring = trim(tmp)
    end function

    function rtoa(n) result(newstring)
        implicit none
        integer, parameter :: str_len = 30
        integer, parameter :: ii10 = selected_real_kind(10,100)
        real(kind=8), parameter :: max_flt = 1E4
        !character(5)  :: db_fmt = 'f15.9'

        real(kind=ii10), intent(in) :: n
        character(:), allocatable :: newstring 
        character(len=str_len) :: tmp

        if( abs(n) < max_flt)then
            write(tmp,"(" // FLT_FMT // ")") n
        else
            write(tmp,"(" // DB_FMT // ")") n
        endif
        newstring = trim(adjustl(tmp))
    end function

    logical function ishuge_4(r) result(ishuge)
        real(kind=4), intent(in) :: r

        ishuge = .false.
        if(abs(huge(r)-r) < 1.0E10*tiny(r)) ishuge=.true.
    end function

    logical function ishuge_8(r) result(ishuge)
        real(kind=8), intent(in) :: r

        ishuge = .false.
        if(abs(huge(r)-r) < 1.0E10*tiny(r)) ishuge=.true.
    end function

    logical function ishuge_16(r) result(ishuge)
        real(kind=16), intent(in) :: r

        ishuge = .false.
        if(abs(huge(r)-r) < 1.0E10*tiny(r)) ishuge=.true.
    end function

    subroutine delete_file(file_name,feedback)
        implicit none
        character(*),intent(in) :: file_name
        logical, optional, intent(in) :: feedback

        logical :: deleted ! whether or not there was a file to be deleted

        ! Check that file exists:
        inquire( file=trim(file_name), exist=deleted)

        if(deleted) then
            if(present(feedback)) then
                if(feedback) write(*,'("Deleting file: ", A)') trim(file_name)
            end if
            ! open the file
            open(delete_unit,file=trim(file_name)) 
            ! Delete it if it exists
            close(delete_unit,status='delete')
        end if

    end subroutine delete_file

    integer function newunit(unit)
        integer, intent(out), optional :: unit
        ! local
        integer, parameter :: LUN_MIN=10, LUN_MAX=1000
        logical :: opened
        integer :: lun
        ! begin
        newunit=-1
        do lun=LUN_MIN,LUN_MAX
           inquire(unit=lun,opened=opened)
           if (.not. opened) then
              newunit=lun
              exit
           end if
        end do
        if (present(unit)) unit=newunit
    end function newunit

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

    ! Write format for n integers
    function integer_format(n)
        implicit none
        integer, intent(in) :: n
        character(len=fmt_len) :: integer_format

        write(integer_format,'("(",I0,A,")")') n,INT_FMT   ! define the integer format

    end function integer_format

    ! Write format for n doubles
    function double_format(n)
        implicit none
        integer, intent(in) :: n
        character(len=fmt_len) :: double_format

        write(double_format,'("(",I0,A,")")') n,DB_FMT   ! define the integer format

    end function double_format

    subroutine read_doubles_1(arr,str)
        implicit none
        real(kind=ii10),dimension(:), intent(out) :: arr
        character(len=*), intent(in),optional :: str

        if(present(str)) read(read_resume_unit,*)
        read(read_resume_unit,double_format(size(arr))) arr

    end subroutine read_doubles_1

    subroutine read_doubles_2(arr,str)
        implicit none
        real(kind=ii10),dimension(:,:), intent(out) :: arr
        !integer,intent(in) :: n1,n2
        character(len=*), intent(in),optional :: str
        integer :: i2

        if(present(str)) read(read_resume_unit,*)
        do i2=1,size(arr,2)
            read(read_resume_unit,double_format(size(arr,1))) arr(:,i2)
        end do

    end subroutine read_doubles_2

    subroutine read_doubles_3(arr,str,n)
        implicit none
        real(kind=ii10),dimension(:,:,:), intent(out) :: arr
        !integer,intent(in) :: n1,n2,n3
        integer,optional,intent(in),dimension(size(arr,3)) :: n
        character(len=*), intent(in),optional :: str
        integer :: i2,i3,m

        if(present(str)) read(read_resume_unit,*)
        do i3=1,size(arr,3)
            !read(read_resume_unit,*)
            if(present(n)) then
                m=n(i3)
            else
                m=size(arr,2)
            end if
            do i2=1,m
                read(read_resume_unit,double_format(size(arr,1))) arr(:,i2,i3)
            end do
        end do

    end subroutine read_doubles_3

    subroutine write_double(a,str)
        implicit none
        real(kind=ii10), intent(in) :: a
        character(len=*), intent(in),optional :: str

        if(present(str)) write(write_resume_unit,'("'//trim(str)//'")')
        write(write_resume_unit,double_format(1)) a

    end subroutine write_double

    subroutine write_doubles_1(arr,str)
        implicit none
        real(kind=ii10),dimension(:), intent(in) :: arr
        character(len=*), intent(in),optional :: str

        if(present(str)) write(write_resume_unit,'("'//trim(str)//'")')
        if(size(arr)>0) write(write_resume_unit,double_format(size(arr))) arr

    end subroutine write_doubles_1

    subroutine write_doubles_2(arr,str,n)
        implicit none
        real(kind=ii10),dimension(:,:), intent(in) :: arr
        integer, intent(in), optional :: n
        character(len=*), intent(in),optional :: str
        integer :: i

        if(present(str)) write(write_resume_unit,'("'//trim(str)//'")')
        if(present(n)) then
            do i=1,n
                call write_doubles( arr(:,i) )
            end do
        else
            do i=1,size(arr,2)
                call write_doubles( arr(:,i) )
            end do
        end if

    end subroutine write_doubles_2

    subroutine write_doubles_3(arr,str,n)
        implicit none
        real(kind=ii10),dimension(:,:,:), intent(in) :: arr
        integer, intent(in),dimension(size(arr,2)),optional :: n
        character(len=*), intent(in),optional :: str
        integer :: i

        if(present(str)) write(write_resume_unit,'("'//trim(str)//'")')
        do i=1,size(arr,3)
            if(present(n)) then
                call write_doubles( arr(:,:,i), '---------------------------------------',n(i) )
            else
                call write_doubles( arr(:,:,i) )
            end if
        end do

    end subroutine write_doubles_3

    integer function bisearch(array,istart,iend,element)

        implicit none
        real(kind=ii10),dimension(:),intent(in) :: array
        real(kind=ii10),intent(in) :: element
        integer,intent(in) :: istart, iend

        !local variable
        integer low, middle, up

        low = istart
        middle = istart
        up = iend

        if(element < array(low))then 
            bisearch = exit_failure  
            return 
        endif

        if(element > array(up))then
            bisearch = up
            return
        endif

        do while(low <= up)

           middle = low + (up-low) / 2

           if( abs(element- array(middle)) < EPS) then
               bisearch = middle
               return
           elseif(element < array(middle))then
               up = middle - 1
           elseif(element > array(middle))then
               low = middle + 1
           endif

        enddo

        bisearch = up 

        return

    end function bisearch

    integer function bisearch2(array,istart,iend,element)

        implicit none
        real(kind=ii10),dimension(:),intent(in) :: array
        real(kind=ii10),intent(in) :: element
        integer,intent(in) :: istart, iend

        !local variable
        integer low, middle, up

        low = istart
        middle = istart
        up = iend

        if(element > array(low))then 
            bisearch2 = exit_failure  
            return 
        endif

        if(element < array(up))then
            bisearch2 = exit_failure
            return
        endif

        do while(low <= up)

           middle = low + (up-low) / 2

           if( abs(element- array(middle)) < EPS) then
               bisearch2 = middle
               return
           elseif(element > array(middle))then
               up = middle - 1
           elseif(element < array(middle))then
               low = middle + 1
           endif

        enddo

        bisearch2 = up

        return

    end function bisearch2

    integer function brute_search(array,element)
        implicit none
        real(kind=ii10), dimension(:,:), intent(in) :: array
        real(kind=ii10), dimension(:)  :: element

        real(kind=ii10), parameter :: eps = 1e-5
        integer i

        do i = 1, size(array,2)
            if(all(abs(array(:,i) - element) < eps))then
                brute_search = i
                return
            endif
        enddo

        brute_search = -1

        return
    end function

    integer function create_dir(path)
        character(len=*), intent(in) :: path
        character(len=1,kind=c_char), dimension(len_trim(path)+1) :: c_path

        c_path =  f_c_string_func(path)
        create_dir = makedir(c_path)

    endfunction

    ! generate an random seed, for parallel processes, let users to provide
    ! an pid since there is no standard way to get pid
    integer function init_random_seed(pid_in) result(seed)
        implicit none
        integer, intent(in), optional :: pid_in
        integer :: un, istat, dt(8), pid, t(2), s
        integer(8) :: count, tms

        ! First try if the OS provides a random number generator
        open(newunit=un, file="/dev/urandom", access="stream", &
          form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
          read(un) seed
          close(un)
        else
          ! Fallback to XOR:ing the current time and pid. The PID is
          ! useful in case one launches multiple instances of the same
          ! program in parallel.
          call system_clock(count)
          if (count /= 0) then
            t = transfer(count, t)
          else
            call date_and_time(values=dt)
            tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
              + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
              + dt(3) * 24 * 60 * 60 * 60 * 1000 &
              + dt(5) * 60 * 60 * 1000 &
              + dt(6) * 60 * 1000 + dt(7) * 1000 &
              + dt(8)
            t = transfer(tms, t)
          end if
          s = ieor(t(1), t(2))
          if(present(pid_in)) pid = pid_in + 1099279 ! Add a prime
          s = ieor(s, pid)
          seed = s
        endif

    end function init_random_seed

    ! See: http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
    subroutine init_random_seed_internal (pid_in)

      implicit none

      integer, intent(in), optional :: pid_in
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid, t(2), s
      integer(8) :: count, tms

      call random_seed(size = n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream", &
        form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
        read(un) seed
        close(un)
      else
        ! Fallback to XOR:ing the current time and pid. The PID is
        ! useful in case one launches multiple instances of the same
        ! program in parallel.
        call system_clock(count)
        if (count /= 0) then
          t = transfer(count, t)
        else
          call date_and_time(values=dt)
          tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
            + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
            + dt(3) * 24 * 60 * 60 * 60 * 1000 &
            + dt(5) * 60 * 60 * 1000 &
            + dt(6) * 60 * 1000 + dt(7) * 1000 &
            + dt(8)
          t = transfer(tms, t)
        end if
        s = ieor(t(1), t(2))
        if(present(pid_in)) pid = pid_in + 1099279 ! Add a prime
        s = ieor(s, pid)
        if (n >= 3) then
          seed(1) = t(1) + 36269
          seed(2) = t(2) + 72551
          seed(3) = pid
          if (n > 3) then
            seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          end if
        else
          seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        end if
      end if
      call random_seed(put=seed)
    end subroutine init_random_seed_internal
end module
