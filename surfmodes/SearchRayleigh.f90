subroutine SearchRayleigh(GRT, c0, roots, ncr, allmodes, ierr)
    use m_GRT, only : pi, dp, T_GRT, nmode
    use Rayleigh, only : startl, SecFunSurf, SecFunSt, SecFunSt2, SecFuns, &
        init_rayleigh, delete_rayleigh
    implicit none
    type(T_GRT), intent(inout) :: GRT
    real(kind=dp), intent(in) :: c0
    real(kind=dp), dimension(nmode), intent(out) :: roots
    integer, intent(out) :: ncr
    integer, intent(in) :: allmodes
    integer, intent(out) :: ierr
    !interface
    !    subroutine sort2(arr,iq)
    !    integer,intent(in)::iq
    !    real*8,intent(inout)::arr(:,:)
    !    end subroutine sort2
    !end interface
    real*8,external::bisecim, N_cf
    real*8 c,k1,k2,kt,f1,f2,f1t,f2t
    real*8 c3(3)
    real(kind=dp) freq
    real(kind=dp), dimension(20000) :: ccc
    
    integer i, ip, j, i2
    integer ib, iq, NN
    integer index0, im1
    integer ilay
    real(kind=dp) cr1, cray, imf
    real(kind=dp) cr(3,nmode) 
    real(kind=dp), dimension(nmode) :: root2
    !real t1, t2
    
    ierr = 0
    freq = grt%w/(2*pi)
    
    NN=nint(N_cf(grt%vsy,GRT)-N_cf(grt%vsm,GRT))
    !print*,NN,vsy,vsm,w; stop
    ccc = 0
    im1 = 0
    index0 = 0
    !call cpu_time(t1)
    call c_interval(GRT,ccc,index0,im1)
    !call cpu_time(t2)
    !write(*,*) 'Time of interval: ', t2-t1
    
    !!!! test C_interval
    !open(33,file='interval_1p36.txt')
    !do ip=1,index0
    !write(33,*) ccc(ip)
    !enddo
    !stop
    !!!! test C_interval
    
    call init_rayleigh(grt%nlayers)
    ncr = 0 ! 对应该频率下的根的数目
    ! used only once to exclude roots very near to layer velocities
    !  tol0=1d-7   
    !tol0=10*tol
    !tol0=5d-6
    
    !! search fundamental mode or Stoneley mode first
    !cr(1,1)=0.
    !cr1=0.
    !call cpu_time(t1)
    cray = 0.0
    if(grt%ifs==0) then
        cray = c0
        !if(grt%nlvl1>0) cray = 0.0
        call FundaMode(grt,ccc,index0,im1,cray,ierr)
        cr1 = cray
        !cray=cr1
    else
        cray = c0
        if(cray<=0) call St_Finder(grt%ifs,grt,cray)
        call StMode(grt,ccc,index0,im1,cray,ierr)
    endif
    !call cpu_time(t2)
    !write(*,*) 'Time of search: ', t2-t1
    ncr = ncr + 1
    roots(ncr) = cray
    cr(:,ncr) = [cray,0.0_dp,dble(grt%ll)] 
    !write(80,*) index_a, freq, cray
    
    !!!!!!!!!!!!TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   print*,w/(2*pi)
    !  !open(90,file='testLVLfull.txt')
    !   open(90,file='testLay5full.txt')
    !   open(90,file='testLVL5.796.txt')
    !   open(90,file='testLVL5.820.txt')
    !   open(90,file='testLay2_5.808.txt')
    !   open(90,file='testLVL5.808b.txt')
      ! open(90,file='testL1_156p45.txt')
      ! open(90,file='testL3_5p8086.txt')
       !open(90,file='testL0_0p023.txt')
     !  open(90,file='testL0_1Hz.txt')
      ! open(90,file='testSt_6p324.txt')
      !!k1=0.11d0; k2=1.80d0
      ! k1=3.85d0; k2=3.94d0
      ! open(90,file='testL1_0p216Hz.txt')
    !   ll=nly
    !   open(90,file='testL0_1p0129Hz_media05.txt')
    !   print*,'test:',freq
    !   k1=2.3d0; k2=5.8d0
    !  !k1=4.55d0; k2=4.56d0
    !  !k1=5.63d0; k2=5.80d0
    !   do while(k1<k2)
    !      !f1=SecFuns(4,k1,imf)
    !       f1=SecFunSurf(0,k1,imf)
    !      !f1=SecFunSt(ifs,k1,imf)
    !       write(90,*) k1,f1,imf
    !       k1=k1+dc
    !      !k1=k1+1d-5
    !   enddo
    !   close(90)
    !   stop
      ! k1=4.03800476719322d0
      ! f1=SecFuns(1,k1,imf)
      ! print*,k1,f1,imf
      ! k1=k1+tol0
      ! f1=SecFuns(1,k1,imf)
      ! print*,k1,f1,imf
      ! k1=4.03800482602158d0
      ! f1=SecFuns(1,k1,imf)
      ! print*,k1,f1,imf
      ! stop
    ! k1=5.09602519847487d0
    ! print*,SecFuns(1,k1,imf)
    ! k2=5.548d0
    ! print*,SecFuns(3,k2,imf)
    !!!!!!!!!!!!TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !! first apply continuity condition on the low velocity layer boundary, not the first layer
    if(allmodes == 1) then
        k2=ccc(1)
        if(im1>2) then
            ilay=grt%lvls(1)
            call startl(k2,grt)
            f2=SecFuns(ilay,k2,GRT,imf)
            do ip=2,im1-1
                k1=k2
                f1=f2
                k2=ccc(ip)
                ilay=grt%lvls(1)
                call startl(k2,grt)
                f2=SecFuns(ilay,k2,grt,imf)
                ! here assume that there do not exist two different
                ! secular functions that have opposite sign in 
                ! certain interval and that the mode can be found
                ! in one function but not in another.
                if(f1*f2<0.) then
                    kt=bisecim(SecFuns,ilay,k1,k2,f1,f2, GRT, iq)
                    if(iq==0) then
                        c=kt
                        if(abs(c-cr1)>grt%tol*4) call check(ierr)
                        if(ierr==0)then
                            ncr = ncr + 1
                            roots(ncr) = c
                            cr(:,ncr) = [c,dble(ilay),dble(grt%ll)]
                        endif
                    endif
                else
                    do j=grt%nlvl1,2,-1
                    ilay=grt%lvls(j)
                    if(k1>grt%vs(ilay)) then
                        f1t=SecFuns(ilay,k1,grt,imf)
                        f2t=SecFuns(ilay,k2,grt,imf)
                        if(f1t*f2t<0) then
                            kt=bisecim(SecFuns,ilay,k1,k2,f1t,f2t,grt,iq)
                            if(iq==0) then
                                c=kt
                                if(abs(c-cr1)>grt%tol*4) call check(ierr)
                                if(ierr==0)then
                                    ncr = ncr + 1
                                    roots(ncr) = c
                                    cr(:,ncr) = [c,dble(ilay),dble(grt%ll)]
                                endif
                            endif
                            exit
                        endif
                    endif
                    enddo
                endif
            enddo
        endif
        
        !!!! then search modes of phase velocities > vs1
        ilay=1
        k2=ccc(im1)
        call startl(k2,grt)
        f2=SecFuns(ilay,k2,grt,imf)
        do ip=im1+1,index0
            k1=k2
            f1=f2
            k2=ccc(ip)
            ilay=1
            call startl(k2,grt)
            f2=SecFuns(ilay,k2,grt,imf)
            iq=-1
            if(f1*f2<0.) then
                kt=bisecim(SecFuns,ilay,k1,k2,f1,f2,grt,iq)
                if(iq==0) then
                    c=kt
                    ! the 1st SecFun may also find St
                    if(abs(c-cr1)>grt%tol*4) call check(ierr)
                    if(ierr==0)then
                        ncr = ncr + 1
                        roots(ncr) = c
                        cr(:,ncr) = [c,dble(ilay),dble(grt%ll)]
                    endif
                endif
            endif
            if(iq/=0) then
                do j=grt%no_lvl,1,-1
                ilay=grt%lvls(j)
                if(k1>grt%vs(ilay)) then
                    f1t=SecFuns(ilay,k1,grt,imf)
                    f2t=SecFuns(ilay,k2,grt,imf)
                    if(f1t*f2t<0) then
                        kt=bisecim(SecFuns,ilay,k1,k2,f1t,f2t,grt,iq)
                        if(iq==0) then
                            c=kt
                            if(abs(c-cr1)>grt%tol*4) call check(ierr)
                            if(ierr==0)then
                                ncr = ncr + 1
                                roots(ncr) = c
                                cr(:,ncr) = [c,dble(ilay),dble(grt%ll)]
                            endif
                            exit
                        endif
                    endif
                endif
                enddo
            endif
        enddo
        
        !! insert the fundamental or St mode into the proper place
        !call sort(cr,ncr,1)
        do i=1,ncr-1
            if(cr(1,i)>cr(1,i+1)) then
                c3=cr(:,i)
                cr(:,i)=cr(:,i+1)
                cr(:,i+1)=c3
            else
                exit
            endif
        enddo
        
        !! to determine if roots are lost
        if(freq>.15d0 .and. grt%index_a>1) then
           !i1=allroots(index_a-1,3)
           !NN0=allroots(index_a-1,2)
           root2=roots
           call SearchLost
        endif
        !allroots(index_a,1:3) = [freq,dble(NN),dble(ncr)]
        !allroots(index_a,4:ncr+3) = cr(1,1:ncr)
    endif
    
    !! to be used in the next iteration
    !i1=ncr; NN0=NN
    !root1(1:i1)=cr(1,1:ncr)
    ! free memory
    call delete_rayleigh

contains
    subroutine check(ierr)
        implicit none
        integer, intent(out), optional :: ierr
        integer i
    
        ierr = 0
        ! exclude layer velocity
        do i=1,grt%ilastvs
           !if(abs(c-v(i))<tol0) exit
            if(abs(c-grt%v(i))<2d-5) exit
        enddo
        if(i==grt%ilastvs+1) then
            ierr = 0
            !ncr=ncr+1
            !cr(:,ncr)=[c,dble(ilay),dble(grt%ll)]
            !if(index_a>593) print*,ncr,c,ilay
            ! print*,ncr,c,ilay
        else
            ierr = 1
        endif
    end subroutine check
    
    subroutine StMode(grt,ccc,index0,im1,cSt,ierr)
        implicit none
        type(T_GRT), intent(inout) :: grt
        real(kind=dp), dimension(:), intent(in) :: ccc
        integer, intent(in) :: index0, im1
        real(kind=dp), intent(inout) :: cSt
        integer, intent(out) :: ierr
        real*8,parameter::dc=1d-3
        real*8 cmn,cmx !,kk(nk)
        integer im2
        !real*8::cSt=0.
        !if(.not.(cSt>0.)) call St_Finder(ifs,cSt)

        ierr = 1
        !write(*,*) 'Estimated stoneley mode: ', cSt
        cmn=cSt*.75d0
        !write(*,*) 'Estimated stoneley mode: ', cmn, cSt
        cmx=min(grt%vs(grt%ifs+1),grt%vp(grt%ifs))-4*grt%tol
        call startl(cmx,grt)
        !ll=nly
        iq = -1
        k2=cmx
        f2=SecFunSt(grt%ifs,k2,grt,imf)
        do 
            k1=k2
            f1=f2
            k2=k2-dc
            if(k2<cmn) exit
            f2=SecFunSt(grt%ifs,k2,grt,imf)
            if(f1*f2<0.) then
                kt=bisecim(SecFunSt,grt%ifs,k1,k2,f1,f2,grt,iq)
                if(iq==0) then
                    !c=kt
                    cSt = kt
                    ilay=-grt%ifs
                    ierr = 0
                    !call check(ierr)
                    !if(ierr == 0)then
                    !    cSt = c
                    !    exit
                    !endif
                   !cr(:,ncr)=[c,dble(0),dble(ll)]
                    exit
                endif
            endif
        enddo
        
        if(iq /= 0 )then
            !print *, 'Search larger than vs1: ', freq, im1
            ! search larger than vs1
            im2 = 1
            do while(ccc(im2)<cmx)
                im2 = im2 + 1
                if(im2==index0)then
                    !print *, index0, ccc(index0), cmx, grt%vsy, grt%vs(grt%ifs+1)
                    exit
                endif
            enddo
            k2=ccc(im2)
            call startl(k2,grt)
            f2=SecFunSt(grt%ifs,k2,grt,imf)
            do ip=im2+1,index0
                k1=k2
                f1=f2
                k2=ccc(ip)
                ilay=1
                call startl(k2,grt)
                f2=SecFunSt(grt%ifs,k2,grt,imf)
                iq=-1
                if(f1*f2<0.) then
                    kt=bisecim(SecFunSt,grt%ifs,k1,k2,f1,f2,grt,iq)
                    if(iq==0) then
                        cSt = kt
                        ierr = 0
                        exit
                        !c=kt
                        !call check(ierr)
                        !if(ierr == 0)then
                        !    cSt = c
                        !    exit
                        !endif
                        !print *, c
                    endif
                endif
            enddo
        
        endif
        if(iq /= 0) cSt = 0
    !print*,cst,cmn,cmx
    ! sorted in decreasing order
    ! kk(11:nk)=cmn+(cmx-cmn)/dble(90)*[90:1:-1]
    ! kk(1:11)=kk(11)+([1:11]-1)*((kk(12)-kk(11))/11.)
    ! !kk(1:nk)=w/kk(1:nk) ! from small to great
    ! call startl(kk(1))
    
    ! TEST
    ! print*,cSt,cmn,cmx
    !open(97,file='testSt_1p03.txt')
    !open(97,file='testL1_0p68.txt')
    !!open(97,file='testL5_0p68.txt')
    !open(97,file='testSt_1p226.txt')
    !open(97,file='testL2_1p226.txt')
    !open(97,file='testL1_1p275.txt')
    !open(97,file='testL2_0p68b.txt')
    !open(97,file='testSt_0p68b.txt')
    !open(97,file='testL1_4p61.txt')
    !!open(97,file='testSt_1p275.txt')
    !open(97,file='testSt_6p868.txt')
    !ll=nly
    !open(97,file='testL1_0p534.txt')
    !!open(97,file='testSt_7p856.txt')
    !!k2=cmn;k1=cmx
    !k2=4.20d0; k1=4.38d0
    !do while(k2<k1)
    !    f2=SecFuns(1,k2,imf)
    !   !f2=SecFunSt(ifs,k2,imf)
    !    write(97,*) k2,f2,imf
    !    k2=k2+dc
    !enddo
    !close(97)
    !stop 'St'
    !k2=1.54911657405712d0
    !k2=3.90128755145500d0
    !f2=SecFunSt(ifs,k2,imf)
    !print*,k2,f2,imf
    !stop 'St'
    ! TEST
    end subroutine StMode
    
    subroutine FundaMode(grt,ccc,index0,im1,cray,ierr)
    implicit none
    type(T_GRT), intent(inout) :: grt
    real(kind=dp), dimension(:), intent(in) :: ccc
    integer, intent(in) :: index0, im1
    real(kind=dp), intent(inout) :: cray
    integer, intent(out) :: ierr
    real*8 cmn,cmx,kk(100)
    
    integer nk
    integer i
    real(kind=dp) imf
    integer im2
    
    ierr = 1
    if(cray>0.) then
        !cmn=cray-.1d0
        cmn = 0.90*cray
        nk=10
    else
        call CR0_Finder(grt%vs1,grt%vp(1),cmn)
        cmn=cmn-.1d0
        nk=100
    endif
    
    cmx=grt%vs1-4*grt%tol
    ! sorted in increasing order
    !kk(1:nk)=cmn+(cmx-cmn)/dble(nk)*[1:nk]
    iq = -1
    if(cmx>cmn)then
        do i = 1, nk
            kk(i) = cmn + (cmx-cmn)/dble(nk)*i
        enddo
        
        call startl(kk(nk),grt)
        
        j=1
        iq = 1
        k2=kk(j)
        f2=SecFunSurf(0,k2,grt,imf)
        do
            k1=k2
            f1=f2
            j=j+1
            k2=kk(j)
            f2=SecFunSurf(0,k2,grt,imf)
            if(f1*f2<0.) then
                kt=bisecim(SecFunSurf,0,k1,k2,f1,f2,grt,iq)
                if(iq==0) then
                    c=kt
                    ilay=0
                    ierr = 0
                    !call check(ierr)
                    !if(ierr==1)then
                    !    iq = 1
                    !else
                    cray = c
                        !print *, freq, c
                    exit
                    !endif
                   !cr(:,ncr)=[c,dble(0),dble(ll)]
                endif
            endif
            if(j==nk) exit
        enddo
    endif
    
    if( iq /= 0 .and. grt%nlvl1>0)then
        ! exist of low velocity layers (vs < first layers)
        ! search lower velocity
        !print *, 'Search lower than vs1'
        !im2 = im1
        !do while(ccc(im2)>cmn)
        !    im2 = im2 - 1
        !    if(im2==0) exit
        !enddo

        !ierr = 1
        !if(im2>1) then
        !    k2=ccc(im2)
        !    call startl(k2,grt)
        !    f2=SecFunSurf(0,k2,grt,imf)
        !    do ip=im2-1, -1, 1
        !        k1=k2
        !        f1=f2
        !        k2=ccc(ip)
        !        ilay=grt%lvls(1)
        !        call startl(k2,grt)
        !        f2=SecFunSurf(0,k2,grt,imf)
        !        ! here assume that there do not exist two different
        !        ! secular functions that have opposite sign in 
        !        ! certain interval and that the mode can be found
        !        ! in one function but not in another.
        !        if(f1*f2<0.) then
        !            kt=bisecim(SecFunSurf,0,k1,k2,f1,f2,grt,iq)
        !            if(iq==0) then
        !                c=kt
        !                call check(ierr)
        !                if(ierr == 0)then
        !                    cray = c
        !                    exit
        !                endif
        !            endif
        !        endif
        !    enddo
        !endif
    
        !if(ierr == 1)then
            ! search c larger than vs1
            k2=ccc(im1)
            call startl(k2,grt)
            f2=SecFunSurf(0,k2,grt,imf)
            do ip=im1+1,index0
                k1=k2
                f1=f2
                k2=ccc(ip)
                ilay=1
                call startl(k2,grt)
                f2=SecFunSurf(0,k2,grt,imf)
                iq=-1
                if(f1*f2<0.) then
                    kt=bisecim(SecFunSurf,0,k1,k2,f1,f2,grt,iq)
                    if(iq==0) then
                        c=kt
                        ierr = 0
                        !call check(ierr)
                        !if(ierr == 0)then
                        cray = c
                        exit
                        !endif
                    endif
                endif
            enddo
        !endif
    
    elseif(iq /= 0 .and.  grt%nlvl1==0)then
        !print *, 'Search larger than vs1: ', freq
        ! search larger than vs1
        im2 = im1
        !if(cmn>cmx)then
            do while(ccc(im2)<cray-2*grt%dc)
                im2 = im2 + 1
                if(im2==index0) exit
            enddo
        !endif

        k2=ccc(im2)
        call startl(k2,grt)
        f2=SecFunSurf(0,k2,grt,imf)
        do ip=im2+1,index0
            k1=k2
            f1=f2
            k2=ccc(ip)
            ilay=1
            call startl(k2,grt)
            f2=SecFunSurf(0,k2,grt,imf)
            iq=-1
            if(f1*f2<0.) then
                kt=bisecim(SecFunSurf,0,k1,k2,f1,f2,grt,iq)
                if(iq==0) then
                    c=kt
                    ierr = 0
                    !call check(ierr)
                    !if(ierr == 0)then
                    cray = c
                    exit
                    !endif
                    !print *, c
                    !exit
                endif
            endif
        enddo
    
    endif
    ! TEST
    !open(97,file='testL1_4p61.txt')
    !open(97,file='testSurf_0p0478.txt')
    !open(97,file='testSurf_0p0358.txt')
    !open(97,file='testSurf_0p0453.txt')
    !open(97,file='testSurf_0p0679.txt')
    !k1=kk(1)
    !k2=kk(nk)
    !do while(k1<k2)
    !    f1=SecFunSurf(0,k1,imf)
    !    write(97,*) k1,f1,imf
    !    k1=k1+dc
    !enddo
    !close(97)
    !stop 'FundaMode'
    ! TEST
    end subroutine FundaMode
    
    subroutine SearchLost
    !implicit none
    !real*8,external::SecFunSt2
    !real*8,parameter::dcm=8d-7
    real*8 a,b
    integer ji
    integer i1
    
    !integer ib,iq,i2
    !common /trivials/ ib,iq,i2
    i1 = grt%ncr1
    i2=ncr
    do i=1,grt%ncr1
        if(root2(i)-grt%root1(i)>grt%dc) then
            iq=1
            do ji=1,-3,-2
            select case(iq)
            case(1)
                if(i>1) then
                    a=root2(i-1)+grt%dcm
                    b=grt%root1(i)-grt%dcm
                else ! i=1
                    a=ccc(1)
                    b=root2(1)-grt%dcm
                endif
                ib=i
                call BeginSearch(a,b)
            case(-1)
                if(i>2) then
                    a=root2(i-2)+grt%dcm
                    b=root2(i-1)-grt%dcm
                else ! i=2
                    a=ccc(1)
                    b=root2(1)-grt%dcm
                endif
                ib=i-1
                call BeginSearch(a,b)
                if(iq==-1) iq=-3
            case(-3)
                a=root2(i-3)+grt%dcm
                b=root2(i-2)-grt%dcm
                ib=i-2
                call BeginSearch(a,b)
                if(iq/=0) then
                    if(grt%ifs/=0) then
                        a=root2(i-1)+grt%dcm
                        b=grt%root1(i)-grt%dcm
                        ib=i
                        grt%ll=grt%nlayers ! 在有水的情况下可能造成极个别mode不显现
                       !iq=1
                       !do jii=1,-1,-2
                       call FinerSearch(SecFunSt2,a,b,c,grt%ifs,grt,iq)
                        if(iq==0) then
                            ncr=ncr+1
                            cr(:,ib+1:ncr)=cr(:,ib:i2)
                            cr(:,ib)=[c,dble(-grt%ifs),dble(grt%ll)]
                            i2=ncr
                            print '(g22.15,2(i3,1x),g22.15,2i3)', &
                                freq,ib,i2,c,1,-grt%ifs
                            !pause
                           ! test how many roots are lost
                            exit
                        else
                            call FinerSearch(SecFunSt,a,b,c,grt%ifs,grt,iq)
                           !print*,a,b
                           !print*,SecFunSt(ifs,a,imf),SecFunSt(ifs,b,imf)
                            if(iq==0) then
                                ncr=ncr+1
                                cr(:,ib+1:ncr)=cr(:,ib:i2)
                                cr(:,ib)=[c,dble(-grt%ifs),dble(grt%ll)]
                                i2=ncr
                                print '(g22.15,2(i3,1x),g22.15,2i3,a)', &
                                    freq,ib,i2,c,1,-grt%ifs,"2"
                                !pause
                               ! test how many roots are lost
                                exit
                            endif
                        endif
                    endif
                    print '(I3,"th mode lost at freq",&
                        &f9.5)',i,freq
                    print*,'root2(i:i-5:-1)'
                    print*,root2(i:i-5:-1)
                    print*,'root1(i:i-5:-1):'
                    print*,grt%root1(i:i-5:-1)
                    print*,'tol0:',grt%tol*4
                    print*,'smin:',grt%smin
                    !pause
                    return
                endif
            end select
            enddo
        endif
    enddo
    do while(i1>i2)
        iq=1
        do ji=1,-1,-2
            select case(iq)
            case(1)
                a=root2(i2)+grt%dcm
                b=grt%root1(i1)-grt%dcm
                ib=ncr+1
                call BeginSearch(a,b)
            case(-1)
                a=root2(i2-1)+grt%dcm
                b=root2(i2)-grt%dcm
                ib=ncr
                call BeginSearch(a,b)
                if(iq/=0) then
                    print '("at least ",I1," modes lost at freq",&
                        &f9.5)',i1-i2,freq
                    !pause
                endif
            end select
        enddo
    end do
    ! now i2>=i1
    
    ! 当频率间隔比较稀疏的时候，有必要做如下判断
    ! then test if a root is lost in this case:
    ! i2==i1 .and. i2==NN .and. i1==NN0
    
    !if(i2==NN0) then
    !    iq=1
    !    do ji=1,-1,-2
    !        select case(iq)
    !        case(1)
    !            a=root2(i2)+dcm
    !            b=vsy-dcm
    !            ib=ncr+1
    !            call BeginSearch(a,b)
    !        case(-1)
    !            a=root2(i2-1)+dcm
    !            b=root2(i2)-dcm
    !            ib=ncr
    !            call BeginSearch(a,b)
    !        end select
    !    enddo
    !endif
    end subroutine SearchLost
    
    subroutine BeginSearch(a,b)
    !real*8,external::SecFuns
    !real*8,parameter::dcm=5d-7
    real*8 a,b,vt
    integer iqt, jl
    !common /trivials/ ib,iq,i2
    iqt=iq
    !if(a>b) then
    !    a=a-dcm+5d-8
    !    b=b+dcm-5d-8
    !endif
    !integer ib,iq
    call startl(b, grt)
    do j=grt%no_lvl+1,1,-1
        jl=grt%lvls(j)
        vt=merge(grt%vs(jl),grt%vp(jl),jl>grt%ifs)
        if(b>vt) then
           !call FinerSearch(SecFuns,a,b,c,jl,tol,dc1,iq)
            call FinerSearch(SecFuns,max(a,vt),b,c,jl,grt,iq)
            ! iq returned as -1 or 0
            if(iq==0) then
                ncr=ncr+1
                cr(:,ib+1:ncr)=cr(:,ib:i2)
                cr(:,ib)=[c,dble(jl),dble(grt%ll)]
                i2=ncr
                print '(g22.15,2(i3,1x),g22.15,2i2)', freq,ib,i2,c,iqt,jl
               ! test how many roots are lost
                exit
            endif
        endif
    enddo
    end subroutine BeginSearch
end subroutine SearchRayleigh

subroutine FinerSearch(f,a,b,xc,jl,grt,iq)
    use m_GRT, only : dp, T_GRT
    implicit none
    real*8,external::bisecim,f
    real*8,intent(in)::a,b
    real*8,intent(out)::xc
    integer,intent(in)::jl
    type(T_GRT), intent(in) :: grt
    integer,intent(out)::iq
    integer,parameter::nd=8,nd1=nd+1
    real*8 xt(nd1),yt(nd1),yt1,xt1,xn,yn,dx
    real(kind=dp) imf
    integer i
    iq=-1
    xt([1,nd1])=[a,b]
    yt([1,nd1])=[f(jl,a,grt,imf),f(jl,b,grt,imf)]
    yt1=yt(1); xt1=xt(1)
    if(yt1*yt(nd1)<0.) then
        xc=bisecim(f,jl,xt1,b,yt1,yt(nd1),grt,iq)
        if(iq==0) return
    endif
    ! into two intervals
    xt(2)=.5d0*(a+b)
    yt(2)=f(jl,xt(2),grt,imf)
    if(yt(2)*yt1<0.) then
        xc=bisecim(f,jl,xt1,xt(2),yt1,yt(2),grt,iq)
        if(iq==0) return
        xc=bisecim(f,jl,b,xt(2),yt(nd1),yt(2),grt,iq)
        if(iq==0) return
    endif
    ! into nd intervals
    xn=(b-a)/nd
    do i=2,nd1
        if(i/=nd1) then
            xt(i)=xt(i-1)+xn
            yt(i)=f(jl,xt(i),grt,imf)
        endif
        if(yt(i)*yt(i-1)<0.) then
            xc=bisecim(f,jl,xt(i-1),xt(i),&
                yt(i-1),yt(i),grt,iq)
            if(iq==0) return
        endif
    enddo
    ! into finer steps
    dx=xn/50
    if(dx>grt%dc2) dx=grt%dc2
    !xn=xt1
    xn=b
    yn=yt(nd1)
    do
       !xn=xn+dx
       !if(xn>b) exit
       !yn=f(jl,xn,imf)
        xt1=xn-dx
        if(xt1<a) exit
        yt1=f(jl,xt1,imf)
        if(yn*yt1<0.) then
            xc=bisecim(f,jl,xt1,xn,yt1,yn,grt,iq)
            if(iq==0) return
        endif
       !yt1=yn; xt1=xn
        yn=yt1; xn=xt1
    enddo
end subroutine FinerSearch

! estimate Stoneley mode
subroutine St_Finder(n,grt,cSt)
    use m_GRT, only : dp, T_GRT
    ! n=ifs
    !use hash,only:vp,vs,rho
    implicit none
    type(T_GRT), intent(in) :: grt
    integer,intent(in)::n
    real(kind=dp),intent(out)::cSt
    real(kind=dp),parameter::tol_St=1d-7
    real(kind=dp) c1,c2,c0,a1,a2,a0,dc
    c2=min(grt%vp(n),grt%vs(n+1))
    c1=c2*.8d0
    c2=c2-(c2-c1)/1d4
    
    a2=getSt(n,c2)
    a1=getSt(n,c1)
    dc=c2-c1
    c0=(c2+c1)/2d0
    
    if(a1*a2<0.) then
        do while(dc>=tol_St)
            a0=getSt(n,c0)
            if(a0*a1<0) then
                a2=a0
                c2=c0
            else
                a1=a0
                c1=c0
            endif
            dc=c2-c1
            c0=(c2+c1)/2d0
        enddo
    else
        stop 'Wrong input for Stoneley mode!'
    endif
    a0=getSt(n,c0)
    if(abs(a0)<.5d0) cSt=c0
    
CONTAINS
    real*8 function getSt(n,x)
        implicit none
        integer, intent(in) :: n
        real*8 x,a,b,c,c1,c2
        a=1-(x/grt%vp(n))**2
        b=1-(x/grt%vp(n+1))**2
        c1=(x/grt%vs(n+1))**2*(grt%rho(n)/grt%rho(n+1))
        c2=(grt%vs(n+1)/x)**2
        c=1-1d0/c2
        getSt=c1*sqrt(b/a)+c2*((1+c)**2-4*sqrt(b*c))
    end function getSt
end subroutine St_Finder

! estimate fundamental mode
subroutine CR0_Finder(v1, v2, CRo)
    implicit none
    real*8,parameter::tol_cr0=1d-7
    real*8     v1, v2, CRo, c, R, DR
    c = 0.8d0*v1  ! initial value
    do
        call Rayhomo
        c = c - R/DR
        if (v1-c<tol_cr0 .or. v2-c<tol_cr0 .or. c/=c) exit
        if ( dabs(R/DR).LT.tol_cr0 ) exit
    end do 
    CRo = c
    if(v1-c<tol_cr0 .or. v2-c<tol_cr0 .or. c/=c ) c = 0.8d0*v1
CONTAINS
    subroutine Rayhomo
        implicit none
        real*8    p, ps, pp, sps, spp, p2, ps2, pp2
        ps  = 1.0/v1
        pp  = 1.0/v2
        p   = 1.0/c
        p2  = p*p
        ps2 = ps*ps
        pp2 = pp*pp
        sps = dsqrt( p2 - ps2 )
        spp = dsqrt( p2 - pp2 )
        R   = ( ps2 - 2.0*p2 )**2 - 4.0*p2*sps*spp  ! Rayleigh function
        DR  = p2*( 8.0*p*( ps2 - 2.0*p2 ) + 8.0*p*sps*spp+&
              4.0*(p2*p)*( spp/sps + sps/spp ) )
    end subroutine Rayhomo
end subroutine CR0_Finder

