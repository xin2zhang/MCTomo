subroutine C_Interval(GRT,ccc,ncc,im1)
    use m_GRT, only : T_GRT, dp, eps
    implicit none
    type(T_GRT), intent(in) :: GRT
    real(kind=dp), dimension(20000), intent(out) :: ccc
    integer, intent(out) :: ncc, im1
    real*8,external::N_cf
    
    real(kind=dp), parameter :: pi = 3.1415926
    real(kind=dp), parameter :: edNN=.5d0
    integer index0
    integer i,j,ij,intemp,NN, NNc,ii
    real*8 vvv(20000),c0,c1,c2,dc,dN,c00
    real*8 hi,Ni,c01,tol0
    real(kind=dp) :: freq
    
    freq  = GRT%w/(2.0*pi)
    NN = nint(N_cf(GRT%vsy,GRT)-N_cf(GRT%vsm,GRT)) ! anint: 四舍五入到最近的整数
    ! If the frequency is below the critical frequency that brings up the
    ! inverse dispersion, or the number of the roots is too small
    ! not necessary to use FYHs method
    !     !!!! FOR THIS CASE !!!!                 ! byme
    if (freq .LT. 0.12 .OR. NN .LT. 2) then
        index0 = (NN+1)*1024   ! byme
        !index0 = (NN+1)*256    ! byme
    
        do i = 1,index0 ! 从v(1)到grt%vsy
            vvv(i) = GRT%v(1)-0.1d0+dble(i)*(grt%vsy-GRT%v(1)+0.1d0)/dble(index0) ! byme
            !vvv(i) = 0.95*GRT%v(1)+dble(i)*(grt%vsy-GRT%v(1)*0.95)/dble(index0) ! byme
        end do
    
        !     --insert 100 points near grt%vsy
        do i=1,100      ! 从0.2*grt%vsy到~grt%vsy
            index0 = index0+1  ! index0又增加了100 ! byme
           !c0 = grt%vsy-i*grt%vsy*0.8d0/100.
            c0=grt%vsy*(1.-.008d0*i)
            vvv(index0) = c0
        end do
       !do i=1,10
       !    c0 = grt%vsy-(i*10)*tol
       !    index0 = index0+1  ! index0又增加了10 ! byme
       !    vvv(index0) = c0
       !end do 
       !index0 = index0+1     ! 再加1个
       !vvv(index0) = grt%vsy-3.0*tol
       !index0 = index0+1     ! 又加1个
       !vvv(index0) = grt%vsy
       !call sort(vvv,index0,1)
    !   open(88,file='ccc') ! byme
    !   write(88,*) ccc(:index0), index0
    else
    
        !     FOR THE OTHER CASE
        !	--get the c sequence according to the N_cf----	   
        index0=0
        c1=GRT%v(1)
        index0=index0+1   ! index is the variable playing the role similar to index0
        vvv(index0)=c1   ! vvv vs. ccc
        c2=grt%vsy
        c0=c2
        NNc  = NN
        
        do while(c0.GT.c1) 
            if(NNc.GT.0)  then
                dc = (c2-c1)/(NNc)   ! interval
            endif 
            c0 = c2-dc
            do
                dN = N_cf(c2,GRT)-N_cf(c0,GRT)
                if (dN.LT.edNN)  then    ! dN<.5
                    index0 = index0+1  ! !!
                    vvv(index0) = c0
                    NNc = NNc-1
                    c2=c0     
                    exit  
                else
                    c0 = (c2+c0)/2.0   ! 使得c0点的根的数目与c2点的根的数目相差小于.5
                endif
            end do
        end do
        
        !	--insert the point according the Ni_cf
        ij=1
        do while(grt%v(ij).LE.grt%vsy)
            ij=ij+1 
        end do
        
        do i=1,ij-2 ! print*, ij
            c01 = GRT%v(i) 
            do j=1,GRT%nlayers
                if (abs(GRT%vs(j)-c01)<eps .OR. abs(GRT%vp(j)-c01)<eps) then    
                    ii = j   !
                endif
            end do
            !hi = z(ii)-z(ii-1)
            hi = GRT%d(ii)
            Ni = 2.0*freq*hi/c01 + eps   ! Ni means?
            do j=1,floor(Ni)
                c00 = c01/sqrt(1.0d0-(dble(j)/Ni)**2)
                if(c00.LE.grt%vsy) then
                    index0 = index0+1    !继续增加index
                    vvv(index0) = c00
                end if
            end do
        end do   
        
        !	--sort vvv(i)
        call sort(vvv,index0,1)
        
        ! insert the middle point of every sequential point
        do j=1,2
            intemp=index0
            do i=1,intemp-1
                index0=index0+1    ! 继续增加index
                vvv(index0)=(vvv(i)+vvv(i+1))/2.0
            end do
        enddo
        ! --insert 100 points near grt%vsy
        do i=1,100      ! 在grt%vsy附近处加密100个点
            c0 = grt%vsy-dble(i)/100.d0*0.1d0
            index0 = index0+1
            vvv(index0) = c0
        end do
    
    endif
    
    !print*,'tol:',tol
    !! finally,
    do i=1,10       ! 又加密10个点
        c0 = grt%vsy-(i*10)*GRT%tol
        index0 = index0+1
        vvv(index0) = c0
    !   print*,c0
    end do 
    index0 = index0+1
    vvv(index0) = grt%vsy-3.0*GRT%tol
    !   print*,c0
    ! add one point a little larger than grt%vs1
    index0 = index0+1
    vvv(index0) = grt%vs1+4*GRT%tol
    index0 = index0+1
    vvv(index0) = grt%vs1-4*GRT%tol
    
    ! sort vvv(i)
    call sort(vvv,index0,1)
    ! and ensure no points are less than vsm
    do i=1,index0
        if(vvv(i)>grt%vsm) then
            ii=index0
            index0=index0-(i-1)+1 ! final points
            vvv(2:index0)=vvv(i:ii)
            vvv(1)=grt%vsm+GRT%tol
            exit
        endif
    enddo
    
    ! ensure no points are greater than grt%vsy
    ii=index0
    do while(vvv(ii)>=grt%vsy) 
        ii=ii-1
    enddo
    index0=ii
    
    ! exclude the same point 
    tol0 = 10*GRT%tol   
    !index0 = 0
    ccc(1)=vvv(1)
    index0=1
    do i=2,ii
        if(vvv(i)-vvv(i-1)<tol0) cycle
        index0=index0+1
        ccc(index0)=vvv(i)
    enddo
    
    ! determine where points begin to be larger than grt%vs1
    im1 = 0
    ncc = index0
    do i=1,index0
        if(ccc(i)>=grt%vs1) then
            im1=i
            exit
        endif
    enddo
    if(im1==0) im1 = index0

end subroutine C_Interval
    
real*8 function  N_cf(c,GRT)
    use m_GRT, only : T_GRT, dp
    implicit none
    real(kind=dp), intent(in) :: c
    type(T_GRT), intent(in) :: GRT
    
    real*8, parameter :: pi = 3.1415926
    complex*16  yp,ys,sum
    integer i
    sum=0.
    do i=1,GRT%nlayers-1
        yp=sqrt(dcmplx((c/GRT%vp(i))**2-1))
        if(GRT%vs(i).GT.0.0) then
            ys=sqrt(dcmplx((c/GRT%vs(i))**2-1))
        else
            ys=0.0
        end if 
        sum=sum+2.0*(yp+ys)*GRT%d(i)/dcmplx(c) 
    end do
    N_cf=(GRT%w/(2.0*pi)*dble(sum))
end function N_cf
