module eigenfunctions_L
use Love
use eigenfunctions,only:c,absym,absym2,y,ym,ym2,&
    eigfun,zz,coef,VUSR,jz,ns

contains
subroutine eig_L
    implicit none
vk=w/c
! the following if-structure are refered to that in SecFuns
call propdn_L(c,1+ifs,ilay)  ! Rud,Tu
call propup_L(c,ll-1,ilay)   ! Rdu,Td
! 得到层系数
coef(1:2,ilay)=[IC,RduL(ilay)]
! 沿相反方向传递
call layerCup_L(ilay-1,1+ifs)
call layerCdn_L(ilay+1,ll)
! 计算最终的波场
call wavefield_L
end subroutine eig_L

subroutine layerCup_L(j2,j1)
    implicit none
integer,intent(in)::j2,j1
integer k

do k=j2,j1,-1
coef(2,k)=coef(2,k+1)*TuL(k)
coef(1,k)=coef(2,k)*RudL(k-1)
enddo
end subroutine layerCup_L

subroutine layerCdn_L(j1,j2)
    implicit none
integer,intent(in)::j1,j2
integer k

do k=j1,j2-1
coef(1,k)=coef(1,k-1)*TdL(k-1)
coef(2,k)=RduL(k)*coef(1,k)
enddo
! for the last layer ll
coef(1,k)=coef(1,k-1)*TdL(k-1)
coef(2,k)=0.
end subroutine layerCdn_L

subroutine wavefield_L
    implicit none
    integer i, j, k
    integer nnode
! for fluid layers
do k=1,ifs
    do j=jz(1,k),jz(2,k)
        y(1)=0.
        y(2)=0.
        VUSR(j,:)=y(1:2)
    enddo
enddo

! for solid layers
do k=1+ifs,ll
    a22=EinvE_L(k-1,c,1)  ! layer matrix
    do j=jz(1,k),jz(2,k)
        y(1:2)=Laz_L(zz(j),k)*coef(1:2,k)
        y(1:2)=y(1)*a22(:,1)+y(2)*a22(:,2)
       !VU(j,:)=y(1:2)
        VUSR(j,:)=y(1:2)
    enddo
enddo

do k=ll+1,nly
    VUSR(jz(1,k):jz(2,k),:)=0.
enddo

! to ensure the maximum of U and S is one
nnode=jz(2,ll)
ym=VUSR(0,1); absym=abs(ym)
ym2=VUSR(0,2); absym2=abs(ym2)
do i=1,nnode
    if(abs(VUSR(i,1))>absym) then
        ym=VUSR(i,1)
        absym=abs(ym)
    endif
    if(abs(VUSR(i,2))>absym2) then
        ym2=VUSR(i,2)
        absym2=abs(ym2)
    endif
enddo
VUSR(0:nnode,1)=VUSR(0:nnode,1)/ym
VUSR(0:nnode,2)=VUSR(0:nnode,2)/ym2

ns=ns+1
eigfun(:,:,ns)=dble(VUSR)

CONTAINS
function Laz_L(zx,il)
    implicit none
complex*16 Laz_L(2)
real*8,intent(in):: zx
integer,intent(in):: il
Laz_L(1)=exp(-vk*(zx-z(il-1))*cs(1))
Laz_L(2)=exp( vk*(zx-z(il))*cs(1))
end function Laz_L
end subroutine wavefield_L

end module eigenfunctions_L
