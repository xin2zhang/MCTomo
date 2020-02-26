module eigenfunctions
use Rayleigh
implicit none
real*8 c,absym,absym2
complex*16 y(4),ym,ym2
real*8,allocatable::eigfun(:,:,:),zz(:)
complex*16,allocatable:: coef(:,:),VUSR(:,:)
integer,allocatable:: jz(:,:)
integer::ns=0,km  ! ns<1000

contains
subroutine eig
implicit none
integer k

vk=w/c
! the following if-structure are refered to that in SecFuns
if(ilay>0) then
    if(ifs>=ilay) then
        call propdn_f(c,1,ilay)
        call propup(c,ll-1,ifs+1)
        call propup_f(c,ifs,ilay)

        coef(1,ilay)=1.
        coef(2,ilay)=Rdu(1,1,ilay)
        call layerCup_f(ilay-1,1)
        call layerCdn_f(ilay+1,ifs)
        ! 转换
        k=ifs+1
        coef(1:2,k)=Td(:,1,ifs)*coef(1,ifs)
        pp=>Rdu(:,:,k)
        coef(3:4,k)=coef(1,k)*pp(:,1)+coef(2,k)*pp(:,2)
        call layerCdn(k+1,ll)
    else
        if(ifs>0) call propdn_f(c,1,ifs)
        call propdn(c,1+ifs,ilay)  ! Rud,Tu
        call propup(c,ll-1,ilay)   ! Rdu,Td
        ! 得到层系数
        pp=>Rdu(:,:,ilay)
        a22=matmul(Rud(:,:,ilay-1),pp)
        do k=1,2
            a22(k,k)=a22(k,k)-IC
        enddo

        coef(1:2,ilay)=[-a22(1,2),a22(1,1)]
        coef(3:4,ilay)=coef(1,ilay)*pp(:,1)+coef(2,ilay)*pp(:,2)
        ! 沿相反方向传递
        call layerCup(ilay-1,1+ifs)
        if(ifs>0) then
            ! 转换
            k=ifs
            coef(2,k)=dot_product(Tu(1,:,k),coef(3:4,k+1))
            coef(1,k)=coef(2,k)*Rud(1,1,k-1)
            call layerCup_f(k-1,1)
        endif
        call layerCdn(ilay+1,ll)
    endif

! for fundamental mode; solid-layer model
else if(ilay==0) then
    call propup(c,ll-1,1)  ! 传递系数，得到指数函数值
    ! 自由表面边界条件
    a44=EinvE(0,c,1)    ! 得到vertical wavenumber
    b22=a44(3:4,3:4)
    do k=1,2
        b22(:,k)=b22(:,k)*la(k)
    enddo
    pp=>Rdu(:,:,1)
    a22=a44(3:4,1:2)+matmul(b22,pp)
    ! 得到层系数
    coef(1:2,1)=[-a22(1,2),a22(1,1)]
    coef(3:4,1)=coef(1,1)*pp(:,1)+coef(2,1)*pp(:,2)
    call layerCdn(2,ll)

! for St mode or a few stubborn modes; ilay = -ifs
else            
    ilay=-ilay
    a33=Stoneley(ilay,c)
    k=ilay+1; pp=>Rdu(:,:,k)
    coef(1:2,k)=[-det2(a33(1:2,2:3)),det2(a33(1:2,[1,3]))]
    coef(3:4,k)=coef(1,k)*pp(:,1)+coef(2,k)*pp(:,2)
    coef(2,ilay)=-sum(a33(1,1:2)*coef(1:2,ilay+1))/a33(1,3)
    coef(1,ilay)=Rud(1,1,ilay-1)*coef(2,ilay)
    call layerCup_f(ilay-1,1)
    call layerCdn(k+1,ll)
endif
! 计算最终的波场
call wavefield
end subroutine eig

function det2(a)
implicit none
complex*16 det2
complex*16 a(2,2)
det2=a(1,1)*a(2,2)-a(1,2)*a(2,1)
end function det2

subroutine layerCup_f(j2,j1)
implicit none
integer,intent(in)::j2,j1
integer k

! j2<ifs
do k=j2,j1,-1
coef(2,k)=coef(2,k+1)*Tu(1,1,k)
coef(1,k)=coef(2,k)*Rud(1,1,k-1)
enddo
end subroutine layerCup_f

subroutine layerCdn_f(j1,j2)
implicit none
integer,intent(in)::j1,j2
integer k

! j1>=2
do k=j1,j2
coef(1,k)=coef(1,k-1)*Td(1,1,k-1)
coef(2,k)=coef(1,k)*Rdu(1,1,k)
enddo
end subroutine layerCdn_f

subroutine layerCup(j2,j1)
implicit none
integer,intent(in)::j2,j1
integer k

do k=j2,j1,-1
!do k=ilay-1,1,-1
ps=>Tu(:,:,k); pp=Rud(:,:,k-1)
coef(3:4,k)=coef(3,k+1)*ps(:,1)+coef(4,k+1)*ps(:,2)
coef(1:2,k)=coef(3,k)*pp(:,1)+coef(4,k)*pp(:,2)
enddo
end subroutine layerCup

subroutine layerCdn(j1,j2)
    implicit none
integer,intent(in)::j1,j2
integer k

do k=j1,j2-1
!do k=2,ll-1
ps=>Td(:,:,k-1); pp=>Rdu(:,:,k)
coef(1:2,k)=coef(1,k-1)*ps(:,1)+coef(2,k-1)*ps(:,2)
coef(3:4,k)=coef(1,k)*pp(:,1)+coef(2,k)*pp(:,2)
enddo
! for the last layer ll
ps=>Td(:,:,k-1)
coef(1:2,k)=coef(1,k-1)*ps(:,1)+coef(2,k-1)*ps(:,2)
coef(3:4,k)=0.
end subroutine layerCdn

subroutine wavefield
    implicit none
    integer i, j, k
    integer nnode
! for fluid layers
do k=1,ifs
    a22=EinvE_f(k-1,c,1)
    do j=jz(1,k),jz(2,k)
        y(1)=exp(-vk*(zz(j)-z(k-1))*cp(1))*coef(1,k)
        y(2)=exp( vk*(zz(j)-z(k))  *cp(1))*coef(2,k)
        y(1:2)=y(1)*a22(:,1)+y(2)*a22(:,2)
        y([2,4])=y(1:2)
        y(3)=0
        y(1)=2*mu0*y(4)/(rho(k)*c*c)
       !VU(j,:)=y(1:2)
        VUSR(j,:)=y
    enddo
enddo

! for solid layers
do k=1+ifs,ll
    a44=EinvE(k-1,c,1)  ! layer matrix
    do j=jz(1,k),jz(2,k)
        y=Laz(zz(j),k)*coef(:,k)
        y=y(1)*a44(:,1)+y(2)*a44(:,2)+&
            y(3)*a44(:,3)+y(4)*a44(:,4)
       !VU(j,:)=y(1:2)
        VUSR(j,:)=y
    enddo
enddo

do k=ll+1,nly
    VUSR(jz(1,k):jz(2,k),:)=0.
enddo

! to ensure the maximum of U and S is one
nnode=jz(2,ll)
ym=VUSR(0,2); absym=abs(ym)
ym2=VUSR(0,km); absym2=abs(ym2)
do i=1,nnode
    if(abs(VUSR(i,2))>absym) then
        ym=VUSR(i,2)
        absym=abs(ym)
    endif
    if(abs(VUSR(i,km))>absym2) then
        ym2=VUSR(i,km)
        absym2=abs(ym2)
    endif
enddo
VUSR(0:nnode,1:2)=VUSR(0:nnode,1:2)/ym
VUSR(0:nnode,3:4)=VUSR(0:nnode,3:4)/ym2

ns=ns+1
eigfun(:,:,ns)=dble(VUSR)

CONTAINS
function Laz(zx,il)
    implicit none
complex*16 Laz(4)
real*8,intent(in):: zx
integer,intent(in):: il
Laz(1:2)=exp(-vk*(zx-z(il-1))*[cp(1),cs(1)])
Laz(3:4)=exp( vk*(zx-z(il))*[cp(1),cs(1)])
end function Laz
end subroutine wavefield

end module eigenfunctions
