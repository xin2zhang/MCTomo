real*8 function zbrent(func,ilay,x1,x2,f1,f2,imf,GRT,iq)
    ! internal procedures passed as actual argument, 
    ! only allowed in CVF and ifort
    !use hash,only:smin,imf
    use m_GRT, only : T_GRT
    implicit none
    real*8,external::func
    integer,intent(in) :: ilay
    real*8, intent(in) :: x1,x2,f1,f2
    real*8, intent(inout) :: imf
    type(T_GRT), intent(inout) :: GRT
    integer,intent(out) :: iq
    real*8, parameter :: eps=epsilon(x1)
    integer, parameter :: itmax=120,DP=kind(9d0)
    integer :: iter
    real*8 :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    iq=-1
    a=x1; b=x2
    !fa=func(ilay,a,imf); fb=func(ilay,b,imf)
    fa=f1; fb=f2
    !if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
    !    zbrent=0.
    !    return
    !endif
    
    c=b
    fc=fb
    do iter=1,itmax
        if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
            c=a ! rename a, b, c and adjust bounding interval d.
            fc=fa
            d=b-a
            e=d
        end if
    
        if (abs(fc) < abs(fb)) then  ! swap
            a=b; b=c; c=a
            fa=fb; fb=fc; fc=fa
        end if
    
        tol1=2.0_dp*eps*abs(b)+0.5_dp*GRT%tol  ! convergence check.
    !   xm=(c-b)  ! 不可行，不然会导致迭代次数过大而终止
        xm=0.5_dp*(c-b)
    !   if (abs(xm) <= tol1 .or. fb == 0.0) then
        if(abs(xm) <= tol1 .and. abs(imf)<GRT%smin) then
            if(abs(fb)<GRT%smin+GRT%smin) then
                iq=0
                zbrent=0.5_dp*(c+b)
            else
                zbrent=0.
            endif
            return  ! root found
        end if
    
        if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
            s=fb/fa      ! attempt inverse quadratic interpolation.
            if (a == c) then
                p=2.0_dp*xm*s
                q=1.0_dp-s
            else
                q=fa/fc
                r=fb/fc
                p=s*(2.0_dp*xm*q*(q-r)-(b-a)*(r-1.0_dp))
                q=(q-1.0_dp)*(r-1.0_dp)*(s-1.0_dp)
            end if
            if (p > 0.0) q=-q  ! check whether in bounds.
            p=abs(p)
            if (2.0_dp*p < min(3.0_dp*xm*q-abs(tol1*q),abs(e*q))) then
                e=d             ! accept interpolation.
                d=p/q
            else
                d=xm            ! interpolation failed; use bisection.
                e=d
            end if
        else                    ! bounds decreasing too slowly; use bisection
            d=xm
            e=d
        end if
    
        a=b                     ! move last best guess to a.
        fa=fb
        b=b+merge(d,sign(tol1,xm), abs(d) > tol1 ) ! evaluate new trial root.
        fb=func(ilay,b,imf)
    end do
    
    !zbrent=b
    !stop 'zbrent exceeded maximum iterations'
end function zbrent

real*8 function bisecim(f,ilay,x1,x2,f1,f2,GRT,iq)
    ! half method and reverse interpolation method
    ! for y = x^2, inverse function is x=sqrt(y)
    !use hash,only:smin,imf
    use m_GRT, only : T_GRT
    implicit none
    real*8,external::f
    real*8,intent(in)::x1,x2,f1,f2
    type(T_GRT), intent(in) :: GRT
    integer,intent(in) ::ilay
    integer,intent(out)::iq
    real*8 imf
    real*8 x(5),y(5),xt1,xt2
    real*8 fa,dx,dxt,u1,u2,imf2,smin2

    integer, parameter :: nmax =  1000
    integer::nc
    !real*8, parameter :: eps = 1E-32

    smin2=GRT%smin*GRT%smin*2.
    !smin2=8d-6  ! 2d-3**2*2
    nc=0
    dx=abs(x1-x2)
    fa=f1   !fa=f(x1)
    x(1:2)=[x1,x2]
    y(1:2)=[fa,f2]
    xt1=(x(1)+x(2))/2.
    do
        x(3)=(x(1)+x(2))/2.
        y(3)=f(ilay,x(3),GRT,imf)   !f(x(3))
        ! 三点牛顿反插值
        u1=(x(2)-x(1))/(y(2)-y(1))
        u2=(x(2)-x(3))/(y(2)-y(3))
        !if(abs(y(3)-y(1))<eps) y(3)=y(1)+eps
        xt2=x(1)-y(1)*(u1-y(2)*((u2-u1)/(y(3)-y(1))))
    !   xt2=Ninterp(y,x,3,0d0) ! 抛物反插值
    !   yt=f(xt2)
        dxt=abs(xt2-xt1)
        dx=dx/2.
        if(dx<dxt) dxt=dx
        if(dxt<GRT%tol) then
           !u1=abs(f(ilay,xt2,imf2))
            u1=f(ilay,xt2,GRT,imf2)
            u2=y(3)
            if(u1*u1+imf2*imf2<smin2 .or. u2*u2+imf*imf<smin2) then
                iq=0
                bisecim=merge(xt2,x(3),abs(u1)<abs(u2))
            else
                iq=-1
                bisecim=0.
            endif
            return
           !if(abs(imf2)<smin.or.abs(imf)<smin) then
           !    iq=0
           !    u2=abs(y(3))
           !    bisecim=merge(xt2,x(3),u1<u2)
           !else
           !    iq=-1
           !    bisecim=0.
           !endif
           !return
        endif
        xt1=xt2
        if(fa*y(3)<0) then
            x(2)=x(3)
            y(2)=y(3)
        else
            x(1)=x(3)
            y(1)=y(3)
        endif
       nc=nc+1
       if(nc>=nmax)then
           iq=-1
           bisecim=0
           exit
       endif
    enddo
end function bisecim

subroutine TwoRoots(f,ilay,xa,xb,GRT,x,nr,imf, iq)
    ! find one or two roots in interval [xa,xb]
    ! xa<xb. tol is the spacing between adjacent roots
    !use hash,only:imf
    use m_GRT, only : T_GRT
    implicit none
    real*8,external::f,bisecim
    integer,intent(in) :: ilay
    real*8,intent(in)::xa,xb
    type(T_GRT), intent(inout) :: GRT
    real*8, intent(inout) :: imf
    real*8,intent(out)::x(2) ! store two roots
    integer,intent(out)::nr,iq
    real*8 a,b,tol0,xt(8),yt(8),yt1,xt1,xn,yn
    integer,parameter::nd=7
    integer iq0,i
    integer, dimension(nd-1) :: temp=(/(i,i=1,nd-1)/)
    a=xa+.5d-7
    b=xb-.5d-7
    if(a>b) then
        !pause 'a>b'
        return
    endif
    xt([1,8])=[a,b]
    yt([1,8])=[f(ilay,a,GRT,imf),f(ilay,b,GRT,imf)]
    yt1=yt(1); xt1=xt(1)
    nr=0
    if(yt1*yt(8)<0.) then
        x(1)=bisecim(f,ilay,xt1,b,yt1,yt(8),GRT,iq0)
        nr=1
        iq=1
    else
        xt(2)=.5d0*(a+b)
        yt(2)=f(ilay,xt(2),GRT,imf)
        if(yt(2)*yt1<0.) then
            x(1)=bisecim(f,ilay,xt1,xt(2),yt1,yt(2),GRT,iq0)
            x(2)=bisecim(f,ilay,xt(2),b,yt(2),yt(8),GRT,iq0)
            nr=2
            iq=2
        else
            xt(2:nd)=a+(b-a)/nd*dble(temp)
            do i=2,nd
                yt(i)=f(ilay,xt(i),GRT,imf)
                if(yt(i)*yt1<0.) then
                    x(1)=bisecim(f,ilay,xt1,xt(i),yt1,yt(i),GRT,iq0)
                    x(2)=bisecim(f,ilay,xt(i),b,yt(i),yt(8),GRT,iq0)
                    nr=2
                    iq=nd
                    return
                endif
            enddo
            xn=xt1
            tol0=3*GRT%tol
    66      xn=xn+tol0
            yn=f(ilay,xn,GRT,imf)
            if(yn*yt1<0.) then
               !x(1)=bisecim(f,xn-tol0,xn,eps,iq0)
                x(1)=bisecim(f,ilay,xt1,xn,yt1,yn,GRT,iq0)
                x(2)=bisecim(f,ilay,xn,b,yn,yt(8),GRT,iq0)
                nr=2
                iq=4
                return
            endif
            if(xn<b) goto 66
        endif
    endif
end subroutine TwoRoots

subroutine det3(a,v)
    implicit none
    complex*16,intent(in)::a(3,3)
    complex*16,intent(out)::v
    complex*16 q1,q2,q3
    q1 = a(2,2)*a(3,3)-a(2,3)*a(3,2)
    q2 = a(2,1)*a(3,3)-a(2,3)*a(3,1) 
    q3 = a(2,1)*a(3,2)-a(2,2)*a(3,1)
    v  = a(1,1)*q1-a(1,2)*q2+a(1,3)*q3
end subroutine det3

subroutine det4(a,value)
    implicit none
    complex*16 a(4,4),value
    complex*16 q1,q2,q3,q4,q5,q6,q7
    complex*16 q8,q9,q10,q11,q12
    complex*16 qq1,qq2,qq3,qq4
    q1=a(3,3)*a(4,4)-a(3,4)*a(4,3)
    q2=a(3,2)*a(4,4)-a(3,4)*a(4,2)
    q3=a(3,2)*a(4,3)-a(3,3)*a(4,2)
    q4=a(3,3)*a(4,4)-a(3,4)*a(4,3)
    q5=a(3,1)*a(4,4)-a(3,4)*a(4,1)
    q6=a(3,1)*a(4,3)-a(3,3)*a(4,1)
    q7=a(3,2)*a(4,4)-a(3,4)*a(4,2)
    q8=a(3,1)*a(4,4)-a(3,4)*a(4,1)
    q9=a(3,1)*a(4,2)-a(3,2)*a(4,1)
    q10=a(3,2)*a(4,3)-a(3,3)*a(4,2)
    q11=a(3,1)*a(4,3)-a(3,3)*a(4,1)
    q12=a(3,1)*a(4,2)-a(3,2)*a(4,1)
    qq1=a(2,2)*q1-a(2,3)*q2+a(2,4)*q3
    qq2=a(2,1)*q4-a(2,3)*q5+a(2,4)*q6
    qq3=a(2,1)*q7-a(2,2)*q8+a(2,4)*q9
    qq4=a(2,1)*q10-a(2,2)*q11+a(2,3)*q12
    value=a(1,1)*qq1-a(1,2)*qq2+a(1,3)*qq3-a(1,4)*qq4
end subroutine det4

subroutine inv4(a0)
    implicit none
    complex*16      a0(4,4)
    complex*16      d, a(4, 4), y(4, 4)
    complex*16      a12a23, a12a24, a13a24, a11a23, a11a24, a11a22
    complex*16      a33a44, a32a44, a32a43, a31a44, a31a43, a31a42
    real*8          epsi04
    
    epsi04=1.e-8
    a=a0
    
    a12a23 = a(1,2)*a(2,3)-a(1,3)*a(2,2)
    a12a24 = a(1,2)*a(2,4)-a(1,4)*a(2,2)
    a13a24 = a(1,3)*a(2,4)-a(1,4)*a(2,3)
    a11a23 = a(1,1)*a(2,3)-a(1,3)*a(2,1)
    a11a24 = a(1,1)*a(2,4)-a(1,4)*a(2,1)
    a11a22 = a(1,1)*a(2,2)-a(1,2)*a(2,1)
    a33a44 = a(3,3)*a(4,4)-a(3,4)*a(4,3)
    a32a44 = a(3,2)*a(4,4)-a(3,4)*a(4,2)
    a32a43 = a(3,2)*a(4,3)-a(3,3)*a(4,2)
    a31a44 = a(3,1)*a(4,4)-a(3,4)*a(4,1)
    a31a43 = a(3,1)*a(4,3)-a(3,3)*a(4,1)
    a31a42 = a(3,1)*a(4,2)-a(3,2)*a(4,1)
    
    y(1,1) = a(2,2)*a33a44-a(2,3)*a32a44+a(2,4)*a32a43
    y(2,1) =-a(2,1)*a33a44+a(2,3)*a31a44-a(2,4)*a31a43
    y(3,1) = a(2,1)*a32a44-a(2,2)*a31a44+a(2,4)*a31a42
    y(4,1) =-a(2,1)*a32a43+a(2,2)*a31a43-a(2,3)*a31a42
    
    y(1,2) =-a(1,2)*a33a44+a(1,3)*a32a44-a(1,4)*a32a43
    y(2,2) = a(1,1)*a33a44-a(1,3)*a31a44+a(1,4)*a31a43
    y(3,2) =-a(1,1)*a32a44+a(1,2)*a31a44-a(1,4)*a31a42
    y(4,2) = a(1,1)*a32a43-a(1,2)*a31a43+a(1,3)*a31a42
    
    y(1,3) = a(4,4)*a12a23-a(4,3)*a12a24+a(4,2)*a13a24
    y(2,3) =-a(4,4)*a11a23+a(4,3)*a11a24-a(4,1)*a13a24
    y(3,3) = a(4,4)*a11a22-a(4,2)*a11a24+a(4,1)*a12a24
    y(4,3) =-a(4,3)*a11a22+a(4,2)*a11a23-a(4,1)*a12a23
    
    y(1,4) =-a(3,4)*a12a23+a(3,3)*a12a24-a(3,2)*a13a24
    y(2,4) = a(3,4)*a11a23-a(3,3)*a11a24+a(3,1)*a13a24
    y(3,4) =-a(3,4)*a11a22+a(3,2)*a11a24-a(3,1)*a12a24
    y(4,4) = a(3,3)*a11a22-a(3,2)*a11a23+a(3,1)*a12a23	
    
    d = a(1,1)*y(1,1)+a(1,2)*y(2,1)+a(1,3)*y(3,1)+a(1,4)*y(4,1)
    if ( cdabs(d).lt.epsi04 ) then
        !pause ' Singular Matrix_4'
        d = epsi04
    endif
    
    a0=y/d
end subroutine inv4
!
!subroutine inv44(a)
!real*8 a(4,4),b(4,4)
!b=0
!forall(k=1:4) b(k,k)=1d0
!call LU(a,4,b,4)
!a=b
!end subroutine inv44

!  Compact and Complex LU
subroutine LUCC(a,n,b,k)
    ! ax=b; the value of x is stored in b
    implicit none
    integer,intent(in)::n,k
    complex(8),intent(inout)::a(n,n),b(n,k)
    integer::i,j,jj,maxl
    complex(8)::rep(n),rex(k)
    real(8)::maxv,t
    do j=1,n-1
        maxv=real(a(j,j))*real(a(j,j))+imag(a(j,j))*imag(a(j,j)) !set the diagnal element to be the max
        maxl=j           ! keep the row index
        do i=j+1,n
            t=real(a(i,j))*real(a(i,j))+imag(a(i,j))*imag(a(i,j))
            if(t>maxv) then
                maxv=t
                maxl=i
            endif
        enddo
    
        ! if max1/=j, switch row and vector b
        if(maxl/=j) then
            rep=a(j,:)
            a(j,:)=a(maxl,:)
            a(maxl,:)=rep
            rex=b(j,:)
            b(j,:)=b(maxl,:)
            b(maxl,:)=rex
        endif
    
        jj=j+1     ! calculate L, U, saving storage
        do i=j+1,n 
            a(i,j)=a(i,j)/a(j,j)
            a(i,jj:n)=a(i,jj:n)-a(i,j)*a(j,jj:n)
        enddo
    enddo
    
    !x=0  !(forward substitution)
    do i=2,n
        b(i,:)=b(i,:)-matmul(a(i,1:i-1),b(1:i-1,:))
        !do j=1,i-1
        !    b(i,:)=b(i,:)-a(i,j)*b(j,:)
        !enddo
    enddo
    
    b(n,:)=b(n,:)/a(n,n) !(back substitution)
    do i=n-1,1,-1
        b(i,:)=b(i,:)-matmul(a(i,i+1:n),b(i+1:n,:))
        !do j=i+1,n
        !    b(i,:)=b(i,:)-a(i,j)*b(j,:)
        !enddo
        b(i,:)=b(i,:)/a(i,i)
    enddo
end subroutine LUCC

SUBROUTINE sort(arr,n,iq)
    IMPLICIT NONE
    integer,intent(in)::iq,n
    real*8,intent(inout)::arr(n)
    !Sorts an array arr into ascending numerical order 
    !by Shell's method (diminishing increment sort).
    !arr is replaced on output by its sorted rearrangement.
    INTEGER :: i,j,inc
    REAL*8 :: v
    inc=1
    do                  ! Determine the starting increment.
        inc=3*inc+1
        if (inc > n) exit
    end do
    
    if(iq>0) then
        do                  ! Loop over the partial sorts.
            inc=inc/3
            do i=inc+1,n    ! Outer loop of straight insertion.
                v=arr(i)
                j=i
                do          ! Inner loop of straight insertion.
                    if (arr(j-inc) <= v) exit
                    arr(j)=arr(j-inc)
                    j=j-inc
                    if (j <= inc) exit
                end do
                arr(j)=v
            end do
            if (inc <= 1) exit
        end do
    else    ! decreasing order
        do                  ! Loop over the partial sorts.
            inc=inc/3
            do i=inc+1,n    ! Outer loop of straight insertion.
                v=arr(i)
                j=i
                do          ! Inner loop of straight insertion.
                    if (arr(j-inc) >= v) exit
                    arr(j)=arr(j-inc)
                    j=j-inc
                    if (j <= inc) exit
                end do
                arr(j)=v
            end do
            if (inc <= 1) exit
        end do
    endif
END SUBROUTINE sort

SUBROUTINE sort2(arr,iq)
    IMPLICIT NONE
    integer,intent(in)::iq
    real*8,intent(inout) :: arr(:,:)
    ! sort according to the first column of arr
    !Sorts an array arr into ascending numerical order 
    !by Shell's method (diminishing increment sort).
    !arr is replaced on output by its sorted rearrangement.
    INTEGER :: i,j,inc,n
    REAL*8::v(size(arr,2))
    n=size(arr,1)
    inc=1
    do                  ! Determine the starting increment.
        inc=3*inc+1
        if (inc > n) exit
    end do
    
    if(iq>0) then
        do                  ! Loop over the partial sorts.
            inc=inc/3
            do i=inc+1,n    ! Outer loop of straight insertion.
                v=arr(i,:)
                j=i
                do          ! Inner loop of straight insertion.
                    if (arr(j-inc,1) <= v(1)) exit
                    arr(j,:)=arr(j-inc,:)
                    j=j-inc
                    if (j <= inc) exit
                end do
                arr(j,:)=v
            end do
            if (inc <= 1) exit
        end do
    else    ! decreasing order
        do                  ! Loop over the partial sorts.
            inc=inc/3
            do i=inc+1,n    ! Outer loop of straight insertion.
                v=arr(i,:)
                j=i
                do          ! Inner loop of straight insertion.
                    if (arr(j-inc,1) >= v(1)) exit
                    arr(j,:)=arr(j-inc,:)
                    j=j-inc
                    if (j <= inc) exit
                end do
                arr(j,:)=v
            end do
            if (inc <= 1) exit
        end do
    endif
END SUBROUTINE sort2

subroutine timecost(difftime)
    implicit none
    real difftime
    integer::ncall=0
    ncall=ncall+1
    write(*,'("time ",i2,": ")') ncall
    if(difftime>24*3600) then
        write(*,'(i2,"d")') int(difftime/24./3600)
        difftime=mod(difftime,24*3600.)
    endif
    if(difftime>3600) then
        write(*,'(i2,"h")') int(difftime/3600)
        difftime=mod(difftime,3600.)
    endif
    if(difftime>60) then
        write(*,'(i2,"m")') int(difftime/60)
        difftime=mod(difftime,60.)
    endif
    write(*,'(f5.2,"s")') difftime
end subroutine timecost

