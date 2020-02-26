subroutine SearchLove(GRT, c0, roots, ncr, allmodes, ierr)
    use m_GRT, only : pi, dp, T_GRT, nmode
    use Rayleigh,only:startl
    use Love, only : init_love, delete_love, SecFuns_L
    implicit none
    type(T_GRT), intent(inout) :: GRT
    real(kind=dp), intent(in) :: c0
    real(kind=dp), dimension(nmode), intent(out) :: roots
    integer, intent(out) :: ncr
    integer, intent(in) :: allmodes
    integer, intent(out) :: ierr


    real*8,external::bisecim, N_cf_L
    real*8 c,k1,k2,kt,f1,f2,f1t,f2t
    real(kind=dp) freq
    real(kind=dp), dimension(20000) :: ccc
    
    integer ip, j, i2
    integer ib, iq, NN
    integer index0, im1
    integer ilay
    real(kind=dp) cray, imf
    real(kind=dp) cr(3,nmode) 
    real(kind=dp), dimension(nmode) :: root2
    
    
    ierr = 0
    freq = grt%w/(2*pi)
    !private check,SearchLost,BeginSearch
    
    NN=nint(N_cf_L(GRT%vsy,GRT)-N_cf_L(GRT%vsm,GRT))
    ccc = 0
    index0 = 0
    im1 = 0
    call c_interval_L(GRT,ccc,index0,im1)
    
    !!!! test C_interval
    !open(33,file='interval_.txt')
    !do ip=1,index0
    !write(33,*) ccc(ip)
    !enddo
    !stop
    !!!! test C_interval
    
    call init_love(grt%nlayers)
    ncr = 0 
    ! used only once to exclude roots very near to layer velocities
    !  tol0=1d-7   
    !tol0=10*tol
    !tol0=5d-6
    if(allmodes == 0)then
        cray = c0
        call FundaMode(grt,ccc,index0,cray,ierr)
        ncr = ncr + 1
        roots(ncr) = cray
        cr(:,ncr) = [cray,0.0_dp,dble(grt%ll)] 
        
    else
        ! search all modes  
        k2=ccc(1)
        if(im1>2) then
            ilay=grt%lvls(grt%L1)
            call startl(k2,grt)
            f2=SecFuns_L(ilay,k2,grt,imf)
            do ip=2,im1-1
                k1=k2
                f1=f2
                k2=ccc(ip)
                ilay=grt%lvls(grt%L1)
                call startl(k2,grt)
                f2=SecFuns_L(ilay,k2,grt,imf)
                ! here assume that there do not exist two different
                ! secular functions that have opposite sign in 
                ! certain interval and that the mode can be found
                ! in one function but not in another.
                if(f1*f2<0.) then
                    kt=bisecim(SecFuns_L,ilay,k1,k2,f1,f2,grt,iq)
                    if(iq==0) then
                        c=kt
                        call check(ierr)
                        if(ierr==0)then
                            ncr = ncr + 1
                            roots(ncr) = c
                            cr(:,ncr) = [c,dble(ilay),dble(grt%ll)]
                        endif
                    endif
                else
                    do j=grt%L1+grt%nlvl1-1,grt%L1+1,-1
                        ilay=grt%lvls(j)
                        if(k1>grt%vs(ilay)) then
                            f1t=SecFuns_L(ilay,k1,grt,imf)
                            f2t=SecFuns_L(ilay,k2,grt,imf)
                            if(f1t*f2t<0) then
                                kt=bisecim(SecFuns_L,ilay,k1,k2,f1t,f2t,grt,iq)
                                if(iq==0) then
                                    c=kt
                                    call check(ierr)
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
        ilay=1+grt%ifs
        k2=ccc(im1)
        call startl(k2,grt)
        f2=SecFuns_L(ilay,k2,grt,imf)
        do ip=im1+1,index0
            k1=k2
            f1=f2
            k2=ccc(ip)
            ilay=1+grt%ifs
            call startl(k2,grt)
            f2=SecFuns_L(ilay,k2,grt,imf)
            iq=-1
            if(f1*f2<0.) then
                kt=bisecim(SecFuns_L,ilay,k1,k2,f1,f2,grt,iq)
                if(iq==0) then
                    c=kt
                    call check
                    if(ierr==0)then
                        ncr = ncr + 1
                        roots(ncr) = c
                        cr(:,ncr) = [c,dble(ilay),dble(grt%ll)]
                    endif
                endif
            endif
            if(iq/=0) then
                do j=grt%no_lvl,grt%no_lvl_fl+1,-1
                ilay=grt%lvls(j)
                if(k1>grt%vs(ilay)) then
                    f1t=SecFuns_L(ilay,k1,grt,imf)
                    f2t=SecFuns_L(ilay,k2,grt,imf)
                    if(f1t*f2t<0) then
                        kt=bisecim(SecFuns_L,ilay,k1,k2,f1t,f2t,grt,iq)
                        if(iq==0) then
                            c=kt
                            call check
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
        
        !! to determine if roots are lost
        if(freq>.15d0 .and. grt%index_a>1) then
           !i1=allroots(index_a-1,3)
           !NN0=allroots(index_a-1,2)
           !root1(1:i1)=allroots(index_a-1,4:i1+3)
            root2=roots
            call SearchLost
        endif
        !allroots(index_a,1:3) = [freq,dble(NN),dble(ncr)]
        !allroots(index_a,3+[1:ncr]) = cr(1,1:ncr)
        !allroots(index_a,4:ncr+3) = cr(1,1:ncr)
        roots(1:ncr) = cr(1,1:ncr)
    endif
    
    !! to be used in the next iteration
    !i1=ncr; NN0=NN
    !root1(1:i1)=cr(1,1:ncr)
    call delete_love
    
contains
    subroutine check(ierr)
        implicit none
        integer, intent(out), optional :: ierr
        integer i

        ierr =  0
        ! exclude layer velocity
        do i=1,grt%ilastvs
           !if(abs(c-v(i))<tol0) exit
            if(abs(c-grt%v(i))<2d-5) exit
        enddo
        if(i==grt%ilastvs+1) then
            ierr = 0
            !ncr=ncr+1
            !cr(:,ncr)=[c,dble(ilay),dble(ll)]
            !if(index_a>593) print*,ncr,c,ilay
            !print*,ncr,c,ilay
        else
            ierr =  1
        endif
    end subroutine check
    
    subroutine FundaMode(grt,ccc,index0,cray,ierr)
        implicit none
        type(T_GRT), intent(inout) :: grt
        real(kind=dp), dimension(:), intent(in) :: ccc
        integer, intent(in) :: index0
        real(kind=dp), intent(inout) :: cray
        integer, intent(out) :: ierr

        integer ilay, ip, iq, verbose
        real(kind=dp) k1, k2, f1, f2
        real(kind=dp) imf, c0

        ierr = 1 
        c0 = cray

        ilay=1+grt%ifs
        k2=ccc(1)
        call startl(k2,grt)
        f2=SecFuns_L(ilay,k2,grt,imf)
        do ip=2,index0
            k1=k2
            f1=f2
            k2=ccc(ip)
            call startl(k2,grt)
            f2=SecFuns_L(ilay,k2,grt,imf)
            iq=-1
            if(f1*f2<0.) then
                kt=bisecim(SecFuns_L,ilay,k1,k2,f1,f2,grt,iq)
                if(iq==0) then
                    call check(verbose)
                    if(verbose==0)then
                        cray = kt
                        ierr = 0
                        exit
                    endif
                endif
            endif
        enddo
        
        if(ierr == 1)then
            ! Finer search
            k2 = 1.1*c0
            call startl(k2,grt)
            f2 = SecFuns_L(ilay,k2,grt,imf)
            do
                k1=k2-grt%dc
                if(k1<grt%vsm) exit
                call startl(k1,grt)
                f1 = SecFuns_L(ilay,k1,grt,imf)
                iq=-1
                if(f1*f2<0.) then
                    kt=bisecim(SecFuns_L,ilay,k1,k2,f1,f2,grt,iq)
                    if(iq==0) then
                        cray = kt
                        ierr = 0
                        !call check
                        exit
                    endif
                endif
                k2 = k1
                f2 = f1
            enddo
        endif

        if(ierr==1) cray=0

    endsubroutine

    subroutine SearchLost
        implicit none
        !real*8,external::SecFunSt2
        !real*8,parameter::dcm=8d-7
        real*8 a,b
        integer i, ji
        integer i1
        !integer ib,iq,i2
        !common /trivials/ ib,iq,i2
        i1=grt%ncr1
        i2=ncr
        do i=1,i1
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
                        print '(I3,"th mode lost at freq",&
                            &f9.5)',i,freq
                        print*,'root2(i:i-5:-1)'
                        print*,root2(i:i-5:-1)
                        print*,'root1(i:i-5:-1):'
                        print*,grt%root1(i:i-5:-1)
                        print*,'tol0:',grt%tol*3
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
                   !a=merge(root2(i2)+grt%dcm,vsm+grt%dcm,i2>0)
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
        real*8 a,b
        integer iqt
        integer j, jl
        !common /trivials/ ib,iq,i2
        iqt=iq
        !if(a>b) then
        !    a=a-dcm+5d-8
        !    b=b+dcm-5d-8
        !endif
        !integer ib,iq
        call startl(b,grt)
        ! lvls(L1)对应的S波速在固体层中最小
        do j=grt%no_lvl+1,grt%no_lvl_fl+1,-1
            jl=grt%lvls(j)
           !vt=merge(vs(jl),vp(jl),jl>ifs)
            if(b>grt%vs(jl)) then
                call FinerSearch(SecFuns_L,a,b,c,jl,grt,iq)
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
    
end subroutine SearchLove

