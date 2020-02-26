module m_GRT

    use iso_c_binding
      
    implicit none

    private

    public :: T_GRT, init_grt
    public :: csq
    public :: dp, nmode
    public :: pi, eps, expo
    public :: ai, IC

    ! static values
    integer, parameter :: dp=c_double
    integer, parameter :: nmode = 100
    complex*16,parameter::ai=(0d0,1d0)
    complex*16,parameter::IC=(1d0,0d0)
    real(kind=dp),parameter::expo=46d0,eps=1d-10
    real(kind=dp),parameter :: edNN=.5d0,pi=3.1415926535897932d0

    ! parameters used in generalized R/T method
    type T_GRT
        integer nlayers
        integer index0, index_a
        real(kind=dp)  w
        real(kind=dp)  smin, tol
        real(kind=dp)  dc, dc2, dcm

        ! the first layer for calculation
        integer ll
        ! the first low velocity layer except water layers
        integer L1
        ! the last low velocity layer
        integer lvlast
        ! the number of low velocity layers in fluid
        integer no_lvl_fl
        ! the number of low velocity layers
        integer no_lvl
        ! the number of layers whose velocity < first vs
        integer nlvl1, nlvls1
        ! the number of fluid layers
        integer ifs
        ! the index of velocity layers in ascending order
        integer, dimension(:), allocatable ::  lvls

        ! velocity and mu
        integer :: ilastvs
        real(kind=dp), dimension(:), allocatable :: d
        real(kind=dp), dimension(:), allocatable :: vs
        real(kind=dp), dimension(:), allocatable :: vp
        real(kind=dp), dimension(:), allocatable :: rho
        real(kind=dp), dimension(:), allocatable :: v
        real(kind=dp), dimension(:), allocatable :: mu
        real(kind=dp) :: mu0
        real(kind=dp) :: vsy, vs1, vsm, vss1

        ! storeage for the root
        real(kind=dp) :: root1(nmode)
        integer ncr1
    endtype T_GRT

contains

    subroutine init_grt(GRT, nlayers)
        type(T_GRT), intent(out) :: GRT
        integer, intent(in) :: nlayers

        GRT%nlayers = nlayers
        GRT%w = 0
        GRT%index0 = 0
        GRT%index_a = 0
        GRT%smin = 1E-4
        GRT%tol = 1E-5
        GRT%dc = 1E-4
        GRT%dc2 = 1E-4
        GRT%dcm = 1E-4

        GRT%ll = 0
        GRT%L1 = 0
        GRT%lvlast = 0
        GRT%no_lvl_fl = 0
        GRT%no_lvl = 0
        GRT%nlvl1 = 0
        GRT%nlvls1 = 0
        GRT%ifs = 0
        GRT%ilastvs = 0

        allocate( GRT%d(nlayers) )
        allocate( GRT%vp(nlayers) )
        allocate( GRT%vs(nlayers) )
        allocate( GRT%rho(nlayers) )
        allocate( GRT%v(2*nlayers) )
        allocate( GRT%mu(nlayers) )
        allocate( GRT%lvls(nlayers/2+1) )
        GRT%d = 0
        GRT%vp = 0
        GRT%vs = 0
        GRT%rho = 0
        GRT%v = 0
        GRT%mu = 0
        GRT%mu0 = 0
        GRT%lvls = 0
        GRT%vsy = huge(grt%vsy)
        GRT%vs1 = 0
        GRT%vss1 = 0
        GRT%vsm = 0
    end subroutine init_grt

    complex*16 function csq(c,vel)
        implicit none
        real(kind=dp) c,vel
        csq=sqrt(dcmplx(1-(c/vel)**2))
       !csq=sqrt(dcmplx(1./(c*c)-1./(vel*vel)))
    end function csq
    
end module m_GRT
