module Love
    use m_GRT, only: dp, T_GRT, csq, IC
    
    implicit none
    
    private
    
    public :: SecFuns_L 
    public :: init_love, delete_love
    
    
    complex*16,allocatable::RduL(:),RudL(:),TdL(:),TuL(:)
    complex*16,private::r
    complex*16::a22(2,2),b22(2,2)
    complex*16::cs(0:1),la(2)
    !$omp threadprivate(RduL, RudL, TdL, TuL, a22, b22, cs, la)
    
CONTAINS
    
    subroutine init_love(nly)
        implicit none
        integer, intent(in) :: nly

        allocate(RduL(1:nly-1),TdL(1:nly-1),&
            RudL(0:nly-2),TuL(1:nly-2))
        RduL = 0
        TdL = 0
        RudL = 0
        TuL = 0
    end subroutine

    subroutine delete_love
        implicit none

        deallocate(RduL,TdL,RudL,TuL)

    end subroutine

    function EinvE_L(j,c,iq,GRT)
        implicit none
        complex*16 EinvE_L(2,2)
        real*8,intent(in):: c
        integer,intent(in):: iq,j  ! interface number
        type(T_GRT), intent(in) :: GRT

        integer k
        real(kind=dp) am, as
        
        do k=1,iq,-1
            as=GRT%vs(j+k)
            cs(k)=csq(c,as)
            am=GRT%mu(j+k)
            select case(k)
            case(1) ! below the interface
                a22(:,2)=[IC,am*cs(k)]
                a22(:,1)=[IC,-a22(2,2)]
            case(0) ! above the interface
                b22(:,1)=am*cs(k)
                b22(:,2)=[-IC,IC]
                b22=b22/(2.*b22(1,1))
            end select
        enddo
        if(iq==0) then
            EinvE_L=matmul(b22,a22)
        else
            EinvE_L=a22
        endif
    end function EinvE_L
    
    subroutine propdn_L(c,j1,j2,GRT)
        implicit none
        real*8,intent(in)::c
        integer,intent(in)::j1,j2
        type(T_GRT), intent(in) :: GRT

        integer j
        real(kind=dp) vk
        
        vk = GRT%w/c
        a22=EinvE_L(GRT%ifs,c,1,GRT)
        la(2)=exp(-GRT%d(1+GRT%ifs)*vk*cs(1))
        RudL(GRT%ifs)=-a22(2,2)*la(2)/a22(2,1)
        
        ! now loop downward over other interfaces 
        do j=j1,j2-1
            a22=EinvE_L(j,c,0,GRT)
            la=[exp(-GRT%d(j)*vk*cs(0)),exp(-GRT%d(j+1)*vk*cs(1))]
            la(1)=la(1)*RudL(j-1)
            RudL(j)=(a22(1,2)-la(1)*a22(2,2))*la(2) &
                /(la(1)*a22(2,1)-a22(1,1))
            TuL(j)=a22(2,1)*RudL(j)+a22(2,2)*la(2)
        enddo
    end subroutine propdn_L
    
    subroutine propup_L(c,j2,j1,GRT)
        implicit none
        real*8,intent(in)::c
        integer,intent(in)::j1,j2
        type(T_GRT), intent(in) :: GRT

        integer j
        real(kind=dp) vk
        
        vk = GRT%w/c
        a22=EinvE_L(j2,c,0,GRT)
        TdL(j2)=exp(-GRT%d(j2)*vk*cs(0))/a22(1,1)
        RduL(j2)=a22(2,1)*TdL(j2)
        ! now loop upward over other interfaces 
        do j=j2-1,j1,-1
            a22=EinvE_L(j,c,0,GRT)
            r=exp(-GRT%d(j+1)*vk*cs(1))*RduL(j+1)
            TdL(j)=exp(-GRT%d(j)*vk*cs(0))/(a22(1,1)+a22(1,2)*r)
            RduL(j)=(a22(2,1)+a22(2,2)*r)*TdL(j)
        enddo
    end subroutine propup_L
    
    real*8 function SecFuns_L(lay,c,GRT,imf)
        implicit none
        integer,intent(in)::lay
        real*8,intent(in)::c
        type(T_GRT), intent(in) :: GRT
        real*8,intent(out)::Imf

        complex*16 dsp
        real(kind=dp) vk
        vk=GRT%w/c
        call propdn_L(c,1+GRT%ifs,lay,GRT)
        call propup_L(c,GRT%ll-1,lay,GRT)
        dsp=IC-RudL(lay-1)*RduL(lay)
        SecFuns_L=aimag(dsp)
        Imf=dble(dsp)
    end function SecFuns_L

end module Love
