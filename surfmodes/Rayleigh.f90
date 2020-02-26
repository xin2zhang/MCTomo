module Rayleigh
    use m_GRT, only : T_GRT, dp, IC, csq, expo

    implicit none

    private

    public :: init_rayleigh, delete_rayleigh
    public :: startl
    public :: SecFunSurf
    public :: SecFunSt, SecFunSt2
    public :: SecFuns

    ! some global variable
    complex*16,target,allocatable::Rdu(:,:,:),Rud(:,:,:),Td(:,:,:),Tu(:,:,:)
    complex*16::a22(2,2),b22(2,2),a33(3,3),b33(3,3),a32(3,2),b3(3)
    complex*16,target::a44(4,4),b44(4,4)
    complex*16,pointer::pp(:,:),ps(:,:)
    complex*16::cp(0:1),cs(0:1),la(2)
    !$omp threadprivate (Rdu,Rud,Td,Tu,a22,b22,a33,b33,a32,b3,a44,b44,pp,ps,cp,cs,la)
    
CONTAINS

    function inv2(a)
        implicit none
        complex*16 inv2(2,2)
        complex*16,intent(in)::a(2,2)
        inv2=reshape([a(2,2),-a(2,1),-a(1,2),a(1,1)],[2,2])/&
            (a(1,1)*a(2,2)-a(2,1)*a(1,2))
    end function inv2

    subroutine init_rayleigh(nly)
        implicit none
        integer, intent(in) :: nly

        if(allocated(Rdu)) deallocate(Rdu)
        if(allocated(Td)) deallocate(Td)
        if(allocated(Rud)) deallocate(Rud)
        if(allocated(Tu)) deallocate(Tu)
        allocate(Rdu(2,2,1:nly-1),Td(2,2,1:nly-1),&
            Rud(2,2,0:nly-2),Tu(2,2,1:nly-2))
        Rdu = 0
        Td = 0
        Rud = 0
        Tu = 0
    end subroutine

    subroutine delete_rayleigh
        implicit none

        deallocate(Rdu,Td,Rud,Tu)

    end subroutine
    
    subroutine startl(c,GRT)
        implicit none
        real(kind=dp),intent(in)::c
        type(T_GRT), intent(inout) :: GRT

        real(kind=dp) su
        real(kind=dp) vk
        integer sl

        integer j

        sl=GRT%nlayers
        ! if(ifs>0) return
        vk=GRT%w/c
        !c=w/vk
        ! if calculating from a water layer, propup_f needs to be modified
        ! if calculating from the first solid layer, cannot use propup, cause
        ! cannot get the R/T coefficients at the solid/fluid boundary
        su=0
        do j=GRT%lvlast,GRT%nlayers-1
            if(c<GRT%vs(j)) then
                su=su+vk*dble(csq(c,GRT%vs(j)))*GRT%d(j)
                if(su>expo) then
                    sl=j
                    exit
                endif
            else
                su=0
            endif
        enddo
        GRT%ll = sl
        ! do j=lvlast,nly-1
        !     if(c<vs(j)) then    ! 没考虑液体层
        !         ipl=j
        !         su=0
        !         do k=ipl,nly-1
        !             su=su+vk*dble(csq(c,vs(k)))*d(k)
        !             if(su>expo) then
        !                 ll=k
        !                 exit
        !             endif
        !         enddo
        !         exit
        !     endif
        ! enddo
    end subroutine startl
    
    ! used only to find fundamental mode
    real(kind=dp) function SecFunSurf(isurf,c,GRT,Imf)
        implicit none
        integer,intent(in)::isurf   ! isurf must be zero
        real(kind=dp),intent(in)::c
        type(T_GRT), intent(in) :: GRT
        real(kind=dp),intent(out)::Imf

        complex*16 dsp
        real(kind=dp) vk
        integer k

        !c=w/vk
        vk=GRT%w/c
        call propup(c,GRT%ll-1,1,GRT)
        a44=EinvE(isurf,c,1,GRT)
        b22=a44(3:4,3:4)
        do k=1,2
            b22(:,k)=b22(:,k)*la(k)
        enddo
        a22=a44(3:4,1:2)+matmul(b22,Rdu(:,:,1))
        dsp=a22(1,1)*a22(2,2)-a22(1,2)*a22(2,1)
        SecFunSurf=dble(dsp)
        Imf=aimag(dsp)
        !print*,'c:',c,w
        !do l=1,ll
        !    print '(I2,":")',l
        !    print '(1(8g11.3))',(FMRay(k,:,l),k=1,4)
        !enddo
        !print*,'a44:'
        !print '(1(8g11.3))',(a44(k,:),k=1,4)
    end function SecFunSurf
    
    ! used only to find Stoneley mode
    real(kind=dp) function SecFunSt(ifs,c,GRT,Imf)
        implicit none
        integer,intent(in)::ifs
        real(kind=dp),intent(in)::c
        type(T_GRT), intent(in) :: GRT
        real(kind=dp),intent(out)::Imf
        complex*16 dsp
        a33=Stoneley(ifs,c,GRT)
        call det3(a33,dsp)
        SecFunSt=dble(dsp)
        Imf=aimag(dsp)
    end function SecFunSt
    
    real(kind=dp) function SecFunSt2(ifs,c,GRT,Imf)
        implicit none
        integer,intent(in)::ifs
        real(kind=dp),intent(in)::c
        type(T_GRT), intent(in) :: GRT
        real(kind=dp),intent(out)::Imf
        complex*16 dsp
        a33=Stoneley(ifs,c,GRT)
        call det3(a33,dsp)
        SecFunSt2=dimag(dsp)
        Imf=dble(dsp)
    end function SecFunSt2
    
    function Stoneley(ifs,c,GRT)
        implicit none
        integer,intent(in)::ifs
        real(kind=dp),intent(in)::c
        type(T_GRT), intent(in) :: GRT

        complex*16 Stoneley(3,3)
        real(kind=dp) vk
        integer k

        !c=w/vk
        vk=GRT%w/c
        call propdn_f(c,1,ifs,GRT)     ! get Rud^(ifs-1)
        call propup(c,GRT%ll-1,ifs+1,GRT)  ! get Rdu^(ifs+1)
        ! obtain E matrix below f/s interface
        a44=EinvE(ifs,c,1,GRT) 
        la=exp(-GRT%d(ifs+1)*vk*[cp(1),cs(1)])
        b33=0.
        b33(1:2,1:2)=Rdu(:,:,ifs+1)
        do k=1,2
            b33(k,1:2)=b33(k,1:2)*la(k)
        enddo
        b22=EinvE_f(ifs-1,c,1,GRT)
        b33(3,3)=Rud(1,1,ifs-1)*exp(-GRT%d(ifs)*vk*cp(1))
        
        a33(:,1:2)=a44(2:4,3:4)
        a33(:,3)=-[b22(1,1),(0d0,0d0),b22(2,1)]
        b33=matmul(a33,b33)
        a33(:,1:2)=a44(2:4,1:2)
        a33(:,3)=-[b22(1,2),(0d0,0d0),b22(2,2)]
        Stoneley=a33+b33
    end function Stoneley
    
    real(kind=dp) function SecFuns(lay,c,GRT,Imf)
        implicit none
        integer,intent(in)::lay
        real(kind=dp),intent(in)::c
        type(T_GRT), intent(in) :: GRT
        real(kind=dp),intent(out)::Imf
        complex*16 dsp

        integer k
        real(kind=dp) vk

        vk=GRT%w/c
        !c=w/vk
        if(GRT%ifs>=lay) then
            ! in fluid: from layer 1 downward to lay
            call propdn_f(c,1,lay,GRT)
            ! now from ll-1 upward to the 1st solid layer
            call propup(c,GRT%ll-1,GRT%ifs+1,GRT)
            call propup_f(c,GRT%ifs,lay,GRT)
            dsp=Rud(1,1,lay-1)*Rdu(1,1,lay)-1
        else
            ! from layer 1 downward to lay
            if(GRT%ifs>0) call propdn_f(c,1,GRT%ifs,GRT)
            call propdn(c,1+GRT%ifs,lay,GRT)
            ! now from ll-1 upward to lay; ll=nly
            call propup(c,GRT%ll-1,lay,GRT)
        
            a22=matmul(Rud(:,:,lay-1),Rdu(:,:,lay))
            do k=1,2
                a22(k,k)=a22(k,k)-IC
            enddo
            dsp=a22(1,1)*a22(2,2)-a22(1,2)*a22(2,1)
        endif
        SecFuns=aimag(dsp)
        Imf=dble(dsp)
    end function SecFuns
    
    function EinvE_f(j,c,iq,GRT)
        implicit none
        complex*16 EinvE_f(2,2)
        real(kind=dp),intent(in):: c
        type(T_GRT), intent(in) :: GRT
        integer,intent(in):: iq,j  ! interface number
        complex*16 xi

        integer k
        real(kind=dp)  ap

        do k=1,iq,-1
            ap=GRT%vp(j+k); cp(k)=csq(c,ap)
            xi=dcmplx(GRT%rho(j+k)*c**2/(2.*GRT%mu0))
            select case(k)
            case(1) ! below the interface
                a22(:,1)=[IC, xi/cp(k)]
                a22(:,2)=[IC, -a22(2,1)]
            case(0) ! above the interface
                b22(:,1)=xi/cp(k)
                b22(:,2)=[IC,-IC]
                b22=b22/(2.*b22(1,1))
            end select
        enddo
        if(iq==0) then
            EinvE_f=matmul(b22,a22)
        else
            EinvE_f=a22
        endif
    end function EinvE_f
    
    subroutine propup_f(c,j2,j1,GRT)
    ! j2层的下界面到j1层的下界面
    ! j2=ll-1
        implicit none
        real(kind=dp),intent(in)::c
        integer,intent(in)::j1,j2
        type(T_GRT), intent(in) :: GRT

        integer j
        real(kind=dp)  vk

        ! j2=ifs f/s界面号
        vk = GRT%w/c
        call up_fs(c,j2,GRT)
        do j=j2-1,j1,-1
            a22=EinvE_f(j,c,0,GRT)
            b22(1,1)=Rdu(1,1,j+1)*exp(-GRT%d(j+1)*vk*cp(1))
            Td(1,1,j)=exp(-GRT%d(j)*vk*cp(0))/ &
                (a22(1,1)+a22(1,2)*b22(1,1))
            Rdu(1,1,j)=(a22(2,1)+a22(2,2)*b22(1,1))*Td(1,1,j)
        enddo
    end subroutine propup_f
        
    subroutine up_fs(c,j,GRT)
        real(kind=dp),intent(in)::c
        integer,intent(in)::j
        type(T_GRT), intent(in) :: GRT
        integer k
        real(kind=dp)  vk

        vk = GRT%w/c
        a44=EinvE(j,c,1,GRT) ! obtain E matrix below f/s interface
        la=exp(-GRT%d(j+1)*vk*[cp(1),cs(1)])
        do k=1,2
            a44(:,k+2)=a44(:,k+2)*la(k)
        enddo
        b22=EinvE_f(j-1,c,1,GRT)
        b22(:,1)=b22(:,1)*exp(-GRT%d(j)*vk*cp(1))
        
        a33(:,1:2)=matmul(a44(2:4,3:4),Rdu(:,:,j+1))+a44(2:4,1:2)
        a33(:,3)=-[b22(1,2),(0d0,0d0),b22(2,2)]
        b3=[b22(1,1),(0d0,0d0),b22(2,1)]
        call lucc(a33,3,b3,1)
        Td(:,1,j)=b3(1:2)
        Rdu(1,1,j)=b3(3)
    end subroutine up_fs
        
    subroutine propdn_f(c,j1,j2,GRT)
        ! j1层的上界面到j2层的上界面
        ! j1=1
        implicit none
        real(kind=dp),intent(in)::c
        integer,intent(in)::j1,j2
        type(T_GRT), intent(in) :: GRT

        integer j
        real(kind=dp) vk

        vk = GRT%w/c
        Rud(1,1,0)=exp(-GRT%d(1)*vk*csq(c,GRT%vp(1)))
        do j=j1,j2-1    ! j is interface number
            a22=EinvE_f(j,c,0,GRT)
            b22(1,1)=Rud(1,1,j-1)*exp(-GRT%d(j)*vk*cp(0))
            la(1)=exp(-GRT%d(j+1)*vk*cp(1))
            Rud(1,1,j)=(a22(1,2)-b22(1,1)*a22(2,2))*la(1) &
                /(b22(1,1)*a22(2,1)-a22(1,1))
            Tu(1,1,j)=a22(2,1)*Rud(1,1,j)+a22(2,2)*la(1)
        enddo
    end subroutine propdn_f
        
    subroutine dn_fs(c,j,GRT)
        real(kind=dp),intent(in)::c
        integer,intent(in)::j
        type(T_GRT), intent(in) :: GRT

        integer k
        real(kind=dp) vk

        vk = GRT%w/c
        a44=EinvE(j,c,1,GRT) ! obtain E matrix below f/s interface
        la=exp(-GRT%d(j+1)*vk*[cp(1),cs(1)])
        do k=1,2
            a44(:,k+2)=a44(:,k+2)*la(k)
        enddo
        b22=EinvE_f(j-1,c,1,GRT)
        b22(:,1)=b22(:,1)*exp(-GRT%d(j)*vk*cp(1))
        
        a33(:,1:2)=a44(2:4,1:2)
        b3=[b22(1,1),(0d0,0d0),b22(2,1)]*Rud(1,1,j-1)
        a33(:,3)=-[b22(1,2),(0d0,0d0),b22(2,2)]-b3
        a32=-a44(2:4,3:4)
        call lucc(a33,3,a32,2)
        Rud(:,:,j)=a32(1:2,:)
        Tu(1,:,j)=a32(3,:)
    end subroutine dn_fs
        
    function EinvE(j,c,iq,GRT)
        implicit none
        complex*16 EinvE(4,4)
        real(kind=dp),intent(in):: c
        integer,intent(in):: iq,j  ! interface number
        type(T_GRT), intent(in) :: GRT
        complex*16 xi

        integer k
        real(kind=dp) ap, as, am

        do k=1,iq,-1
            ap=GRT%vp(j+k); cp(k)=csq(c,ap) ! cp(1)即界面下的量; cp(0)即界面上的量
            as=GRT%vs(j+k); cs(k)=csq(c,as)
            am=GRT%mu(j+k)
            xi=dcmplx(1-(c/as)**2/2.)
            select case(k)
            case(1) ! below the interface
                a44(:,1)=[IC,cp(k),-am*cp(k),-am*xi]
                a44(:,2)=[cs(k),IC,-am*xi,-am*cs(k)]
                a44(:,3:4)=a44(:,1:2)
                pp=>a44(2:3,3:4)
                pp=-pp
            case(0) ! above the interface
                b44(1,:)=[IC,-xi/cp(k),-1/(cp(k)*am),IC/am]
                b44(2,:)=[-xi/cs(k),IC,IC/am,-1/(cs(k)*am)]
                b44(3:4,:)=b44(1:2,:)
                pp=>b44(3:4,2:3)
                pp=-pp
                b44=b44/(2*(1-xi))
            end select
        enddo
        if(iq==0) then
            EinvE=matmul(b44,a44)
        else
            EinvE=a44
        endif
    end function EinvE
        
    subroutine propdn(c,j1,j2,GRT)
        ! j1层的上界面到j2层的上界面
        ! j1=1
        real(kind=dp),intent(in)::c
        integer,intent(in)::j1,j2
        type(T_GRT), intent(in) :: GRT

        real(kind=dp) vk
        integer j, k

        vk = GRT%w/c
        if(j1==1) then
        a44=EinvE(0,c,1,GRT)
        a22=a44(3:4,3:4)
        la=exp(-GRT%d(1)*vk*[cp(1),cs(1)])
        do k=1,2
            a22(:,k)=a22(:,k)*la(k)
        enddo
        b22=-a44(3:4,1:2)
        call LUCC(b22,2,a22,2)
        Rud(:,:,0)=a22
        !Rud(:,:,j1-1)=-matmul(inv2(a44(3:4,1:2)),a22)
        else
            ! get Rud(:,:,ifs), Tu(1,:,ifs)
            call dn_fs(c,j1-1,GRT)
        endif
        
        ! now loop downward over other interfaces 
        do j=j1,j2-1
            a44=EinvE(j,c,0,GRT)
            a22=Rud(:,:,j-1)
            la=exp(-GRT%d(j)*vk*[cp(0),cs(0)])
            do k=1,2
                a22(k,:)=la(k)*a22(k,:)
            enddo
        
            pp=>Tu(:,:,j)
            pp=a44(3:4,3:4)
            b22=a44(1:2,3:4)
            la=exp(-GRT%d(j+1)*vk*[cp(1),cs(1)])
            do k=1,2
                b22(:,k)=b22(:,k)*la(k)
                pp(:,k)=pp(:,k)*la(k)
            enddo
            b22=b22-matmul(a22,pp)
        
            ps=>a44(3:4,1:2)
            a22=matmul(a22,ps)-a44(1:2,1:2)
            call lucc(a22,2,b22,2)
            Rud(:,:,j)=b22
            pp=matmul(ps,b22)+pp
           !Rud(:,:,j)=matmul(inv2(a22),b22)
        enddo
    end subroutine propdn
        
    subroutine propup(c,j2,j1,GRT)
        ! j2层的下界面到j1层的下界面
        ! j2=ll-1
        real(kind=dp),intent(in)::c
        integer,intent(in)::j1,j2
        type(T_GRT), intent(in) :: GRT
        
        integer j, k
        real(kind=dp) vk

        vk = GRT%w/c
        a44=EinvE(j2,c,0,GRT)
        a22=inv2(a44(1:2,1:2))
        la=exp(-GRT%d(j2)*vk*[cp(0),cs(0)])
        do k=1,2
            a22(:,k)=a22(:,k)*la(k)
        enddo
        Td(:,:,j2)=a22
        Rdu(:,:,j2)=matmul(a44(3:4,1:2),a22)
        ! now loop upward over other interfaces 
        do j=j2-1,j1,-1
            a44=EinvE(j,c,0,GRT)
            b22=Rdu(:,:,j+1)
            la=exp(-GRT%d(j+1)*vk*[cp(1),cs(1)])
            do k=1,2
                b22(k,:)=b22(k,:)*la(k)
            enddo
            ! have obtained Λ_u^{j+1}Rdu^{j+1}
        
            a22=inv2(a44(1:2,1:2)+matmul(a44(1:2,3:4),b22))
            la=exp(-GRT%d(j)*vk*[cp(0),cs(0)])
            do k=1,2
                a22(:,k)=a22(:,k)*la(k)
            enddo
            Td(:,:,j)=a22
            b22=a44(3:4,1:2)+matmul(a44(3:4,3:4),b22)
            Rdu(:,:,j)=matmul(b22,a22)
        enddo
    end subroutine propup

end module Rayleigh
