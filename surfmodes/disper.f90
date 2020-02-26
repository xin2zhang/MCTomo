! This program is used to calculate eigensolutions of normal
! modes in horizontal layered model without fluid layers.
! 
! INPUT FILES: 
! input: 3 lines; containing two another filenames for input,
!       e.g.,
!           Media06.dat,
!           SynSurfW.dat
!       and the 3rd line is a flag indicating whether to compute
!       eigenfunctions (1) or not (0);
!
!       if the flag is 1, then another file is needed:
! eigfuninput: the first two lines are filenames, e.g.,
!              'all', 'eigf';
!              3rd line: the number of frequencies at which 
!              eigenfunctions are to be computed, and the
!              maximum frequency of these modes;
!              the following lines: No. of each frequency and 
!              the starting and ending values of n
! all (default): the output when flag is 0 will be as input
! eigf (default): the output when flag is 0 will be as input
!
!
! OUTPUT FILES:
! cr0.txt: df; fundamental mode at each frequency point
! eigf: means the flag is 0; it is a binary file containing
!       each line: angular frequency, phase velocity, ilay, ll
! all: freq, expected number of roots, practical number of roots,
!      and the difference between them
! allroots: freq, the root found, root No.
!
! the following two files are produced if flag is 1:
! info4eigfun.txt: 
!   line1: dz
!   line2: number of nodes
!   line3-N: freq No., freq, number of modes
!            w, c, ilay, ll for each mode
! eigfun: binay file. Containing zz, and V, U, S, R for each mode
          
                
program disper
use eigenfunctions
use eigenfunctions_L
implicit none
character*80 file1,file2
integer expect,pract,ntype
real timeo,timee
real*8 dz
character*10::real_clock(2)=' ', Ctype*2
integer iw, iwp
integer n1, n2
integer nnode, nt, nw
integer i, j, k, ii

call readinput
print '(2(/,a))', ' Check the printed messages above;',&
    ' press Enter if OK, or press Ctrl-C to stop running.'
!pause

select case(modetype)
case(1)
Ctype='_R'
allocate(Rdu(2,2,1:nly-1),Td(2,2,1:nly-1),&
        Rud(2,2,0:nly-2),Tu(2,2,1:nly-2))

! if ifun=1, then calculate eigenfunctions
! read from file 'input'
if(ifun==1) goto 111
call cpu_time(timeo)
open(80,file='cr0'//Ctype//suf)
call date_and_time(date=real_clock(1),time=real_clock(2))
write(80,*) 'executed at: ',real_clock(2),' ',real_clock(1)
write(80,*) 'df=',df
!open(88,file='eigf'//Ctype//suf,form='binary',access='direct',recl=2*8+4+4)
open(88,file='eigf'//Ctype//suf,access='stream')
! index_a is used to record the number of the frequency
index_a = 0
nt=0
!ni=486  ! byme at  5.80859238901878Hz
root2=>cr(1,:)
cr1=>cr(1,1)
if(ifs==0) then
    call CR0_Finder(vs1,vp(1),cray)
    write(80,*) 'estimated cray:',cray
    cray=0.
else
    call St_Finder(ifs,cSt)
    write(80,*) 'estimated cSt:',cSt
endif
do k = 1, 100 ! test cpu time
do i = nf1, nf4    ! nf1th - nf4th frequency points
    if(nlvl1>0) cray = 0.0
    index_a = index_a + 1
    freq=(i*df)
    w = 2.0*pi*freq
    tol = tolmin+(real(nf4t-i+nf1t))*(tolmax-tolmin)/real(nf4t)
    tol0=4*tol
    smin= smin_min+(real(i-nf1t))*(smin_max-smin_min)/real(nf4t)
    call SearchRayleigh
    print '(i3,1x,f9.5,3(1x,i3))',index_a,freq,NN,ncr,NN-ncr
    do ii=1,ncr
        ! angular freq, root, ilay, ll
        write(88) w,cr(1,ii),int(cr(2:3,ii))
    enddo
    nt=nt+ncr
enddo
enddo
call cpu_time(timee)
deallocate(Rdu,Td,Rud,Tu)


! For Love waves
case(2)
Ctype='_L'
allocate(RduL(1:nly-1),TdL(1:nly-1),&
        RudL(0:nly-2),TuL(1:nly-2))

! L1: 固体层中小于第一层固体层S波速的层号，且对应波速最小
! L1: 数组lvls的index number
L1=0
! lvls has been sorted
do i=1,no_lvl
    if(lvls(i)>ifs) then
        if(vs(lvls(i))<vs1) then
            L1=i
            exit
        endif
    endif
enddo

! if ifun=1, then calculate eigenfunctions
! read from file 'input'
if(ifun==1) goto 111
call cpu_time(timeo)
open(80,file='cr0'//Ctype//suf)
call date_and_time(date=real_clock(1),time=real_clock(2))
write(80,*) 'executed at: ',real_clock(2),' ',real_clock(1)
write(80,*) 'df=',df
!open(88,file='eigf'//Ctype//suf,form='binary',access='direct',recl=2*8+4+4)
open(88,file='eigf'//Ctype//suf,access='stream')
! index_a is used to record the number of the frequency
index_a = 0
nt=0
!ni=486  ! byme at  5.80859238901878Hz
root2=>cr(1,:)
cr1=>cr(1,1)
do i = nf1, nf4    ! nf1th - nf4th frequency points
    index_a = index_a + 1
    freq=(i*df)
    w = 2.0*pi*freq
    tol = tolmin+(real(nf4t-i+nf1t))*(tolmax-tolmin)/real(nf4t)
    tol0=4*tol
    smin= smin_min+(real(i-nf1t))*(smin_max-smin_min)/real(nf4t)
    call SearchLove
    print '(i3,1x,f9.5,3(1x,i3))',index_a,freq,NN,ncr,NN-ncr
    do ii=1,ncr
        ! angular freq, root, ilay, ll
        write(88) w,cr(1,ii),int(cr(2:3,ii))
    enddo
    nt=nt+ncr
enddo
call cpu_time(timee)
deallocate(RduL,TdL,RudL,TuL)
end select

!!! END FOR MODETYPE

call timecost(timee-timeo)
write(80,*) 'time consumed (s):',timee-timeo
close(80)

open(30,file='all'//Ctype//suf)
open(40,file='allroots'//Ctype//suf)
do i=1,index_a
    freq=allroots(i,1)
    expect=allroots(i,2)
    pract=allroots(i,3)
    write(30,*) freq,expect,pract,expect-pract
    do j=4,pract+3
        write(40,*) freq,allroots(i,j),j-3
    enddo
enddo
close(30)
close(40)
close(88)
goto 110

!!! CALCULATE EIGENFUNCTIONS !!!

111 continue
open(77,file='eigfuninput')
read(77,'(a80)') file1
read(77,'(a80)') file2
read(77,*) nw,ff4
select case(modetype)
case(1)
    ntype=4
case(2)
    ntype=2
end select
allocate(jz(2,nly),coef(ntype,nly))
dz=vsm/ff4/8.   ! lamda/8
!d(nly)=50
do i=1,nly
    jz(2,i)=ceiling(d(i)/dz)
enddo
nnode=sum(jz(2,:))
allocate(zz(0:nnode),VUSR(0:nnode,ntype),&
    eigfun(0:nnode,ntype,999),stat=i)
if(i/=0) stop 'allocate eigfun failure!'

open(99,file='info4eigfun'//Ctype//suf)
write(99,*) 'dz:',dz
write(99,*) 'the number of nodes:',nnode+1
write(99,*) 'the number of frequencies:',nw
print*, 'dz:',dz
print*, 'the number of nodes:',nnode+1
print*, 'the number of frequencies:',nw

! jz记录每一层的节点序号（从0开始）
jz(1,:)=0
do i=2,nly
    jz(:,i)=jz(:,i)+jz(2,i-1)
enddo
jz(2,1:nly-1)=jz(2,1:nly-1)-1
do i=1,nly
    nnode=(jz(2,i)-jz(1,i))
    do k = 0, nnode
        !zz(jz(1,i):jz(2,i))=z(i-1)+(d(i)/(nnode+1))*[0:nnode]
        zz(jz(1,i)+k)=z(i-1)+(d(i)/(nnode+1))*k
    enddo
enddo

if(file1(1:2)=='de') then
    file1='all'//Ctype//suf
    file2='eigf'//Ctype//suf
endif
!print*,nf1,nf4,trim(file1),trim(file2); stop
open(30,file=file1)
open(88,file=file2,access='stream')
call cpu_time(timeo)
iwp=0
nt=0
select case(modetype)
case(1)
km=merge(4,3,ifs>0)
do i=1,nw
    read(77,*) iw,n1,n2
    ! 从文件'all'中读取第几个频率
    do ii=1,iw-iwp
        read(30,*) freq,j,ncr
        nt=nt+ncr
    enddo

    if(n1==0) then  ! 意味着该频率下的所有modes
        n1=1; n2=ncr
    endif
    print*,'freq:',iw,freq,ncr
    write(99,'(a,i3,f10.6,3(1x,i3))')'freq:',iw,freq,ncr,n1,n2

    do j=nt-ncr+1+(n1-1), nt-(ncr-n2)
        read(88) w,c,ilay,ll
        write(99,'(2g17.8,2i3)') w/(2*pi),c,ilay,ll
        call eig
    enddo
    iwp=iw
enddo

! For Love waves
case(2)
do i=1,nw
    read(77,*) iw,n1,n2
    ! 从文件'all'中读取第几个频率
    do ii=1,iw-iwp
        read(30,*) freq,j,ncr
        nt=nt+ncr
    enddo

    if(n1==0) then  ! 意味着该频率下的所有modes
        n1=1; n2=ncr
    endif
    print*,'freq:',iw,freq,ncr
    write(99,'(a,i3,f10.6,3(1x,i3))')'freq:',iw,freq,ncr,n1,n2

    do j=nt-ncr+1+(n1-1), nt-(ncr-n2)
        read(88) w,c,ilay,ll
        write(99,'(2g17.8,2i3)') w/(2*pi),c,ilay,ll
        call eig_L
    enddo
    iwp=iw
enddo
end select
! end for modetype

print*,'total modes:',ns
write(99,*) 'total modes:',ns
close(99)
call cpu_time(timee)
call timecost(timee-timeo)
close(77)
close(30)
close(88)
! output eigenfunctions
! sequential binary file
open(81,file='eigfun'//Ctype//suf,access='stream')
write(81) zz   !, eigfun(:,:,1:ns)
! 若不使用循环输出，则会出现stack overflow
! 需要increase stack size
do i=1,ns
    write(81) eigfun(:,:,i) !可能要关闭杀毒软件
enddo
close(81)
deallocate(jz,coef,zz,VUSR,eigfun)

110 continue
deallocate(z,d,rho,vs,vp,mu,v,lvls)
end program disper

subroutine ReadInput
use hash
implicit none
real*8,external::N_cf,N_cf_L
character*80 list,file1,file2
!real*8 mu0
integer::indx=0
integer i, j, k
ifs=0
open(33,file='input')
read(33,*) modetype,ifun
! modetype=1 (for Rayleigh), or 2 (for Love)
read(33,'(a80)') file1
read(33,'(a80)') file2
read(33,'(a20)') suf
close(33)
!  Media_data:
open (31,file=file1,status='old',form='formatted')
read (31, '(a80)' ) list
read (31, '(a80)' ) list
read (31,    *    ) nly
write(*,*) 'number of layers: ', nly
read (31, '(a80)' ) list
allocate(z(0:nly),d(nly),rho(nly),vs(nly),&
        vp(nly),mu(nly),v(2*nly),lvls(nly/2+1))
do i = 1, nly
    read (31, * ) j, z(j-1), rho(j), vs(j), vp(j)
    if(vs(j)/=0) then
        indx=indx+1
        v(indx)=vs(j)
    else
        ifs=ifs+1   ! record the number of fluid layers
    endif
    indx=indx+1
    v(indx)=vp(j)
end do
do i=1,indx
    if(v(i)==vs(nly)) then
        jj=i  !Record the nly-th S wave velocity
        exit
    endif
enddo
!z(nly)=z(nly-1)+50
!z(nly)=1.6d0*z(nly-1)
!z(nly)=3d0*z(nly-1)
read(31,*) z(nly)
! 如果所给的数是0~2之间的数，则视为倍数
if(abs(z(nly)-1)<1.) z(nly)=z(nly)*z(nly-1)

print*,'rho, ','vs, ','vp'
do i=1,nly
    write(*,*) real(rho(i)),real(vs(i)),real(vp(i))
end do

! Sort the wave velocity in increasing order
call sort(v,indx,1)
print*,'vs, vp in solid layers, after sorted:'
do i=1,indx
    write(*,*) real(v(i))
end do
print*,'jj=',jj

d=z(1:nly)-z(0:nly-1)
close (31)

!  Calculate the "mu(i)":
j=0
do i = 1, nly
    mu(i) = rho(i)*vs(i)**2
    if(mu(i)/=0.) then
        j=j+1
        mu0=mu0+mu(i)
    endif
end do
mu0=mu0/j   ! normalized constant
mu=mu/mu0

! check the presence of low velocity layers
nlvl1=0; no_lvl=0
ll=nly
!vsy=vs(nly)
vsy=maxval(vs)
! vs1
select case(modetype)
case(1)
    if(ifs>0) then
        vs1=vp(1)
    else
        vs1=vs(1)
    endif
    vsm=v(1) ! 模型中最小的速度值
case(2)
    vs1=vs(1+ifs)
    vsm=minval(vs(ifs+1:nly))! the lowest S wave velocity
end select
!vsm=vs1 
! number of lvls in fluid
no_lvl_fl=0
do i=2,nly-1
    if(vp(i)<vp(i+1) .and. vp(i)<vp(i-1)) then
        print*,'LVL:',i
        no_lvl=no_lvl+1
        lvls(no_lvl)=i
!       if(i>ifs+1) no_lvl_L=no_lvl_L+1
!       if(vs(i)<vsm) vsm=vs(i)
        if(ifs==0) then
            if(vs(i)<vs1) nlvl1=nlvl1+1
        else
            select case(modetype)
            case(1)
                if(vp(i)<vs1) nlvl1=nlvl1+1
            case(2)
                if(vs(i)>0.) then
                    if(vs(i)<vs1) nlvl1=nlvl1+1
                else
                    no_lvl_fl=no_lvl_fl+1
                endif
            end select
        endif
    endif
enddo
!
select case(modetype)
case(1)
    lvls(no_lvl+1)=1
case(2)
    lvls(no_lvl+1)=1+ifs
end select
! lvlast将在startl中用到
if(no_lvl==0) then
    lvlast=1
    print*,'LVLs are absent.'
else
    ! 将低速层根据层速大小进行排序，使得lvls的第一个
    ! 元素对应的层的层速最小
    lvlast=lvls(no_lvl)
    do i=1,no_lvl-1
        do j=1,no_lvl-i
            if(vp(lvls(j))>vp(lvls(j+1))) then
                k=lvls(j)
                lvls(j)=lvls(j+1)
                lvls(j+1)=k
            endif
        enddo
    enddo
    print '(2(a,i2))','nlvl1=',nlvl1,', no_lvl=',no_lvl
    print '(a,(/),1(i4,f7.3,1x,f7.3))', 'LVLs:', &
        (lvls(i),vs(lvls(i)),vp(lvls(i)),i=1,no_lvl)
    print '(a,f7.3)', 'vs1:',vs1
endif
lvlast=max(ifs+2,lvlast)

if(ifun==0) call ReadParas
! file2 (paras.dat) is not used in computing eigenfunctions

CONTAINS
subroutine ReadParas
    implicit none
! parameters used for computing eigenfrequencies
open (12,file=file2,status='old',form='formatted')
!read (12, '(a80)' ) list
read (12, '(a80)' ) list
read (12,    *    ) ff1
read (12,    *    ) ff4
read (12,    *    ) df

select case(modetype)
case(1)
    NNN_max=N_cf(vs(nly),2*pi*ff4)
case(2)
    NNN_max=N_cf_L(vs(nly),2*pi*ff4)
end select
if(df==0.) df = ff4/((NNN_max+1)*2) !~ 1.0Hz/double maximum number of roots
nf1t = (ff1/df)+1 
nf4t = (ff4/df)+1
read (12,    *    ) nf1,nf4
if(nf1==0) then
    nf1 = nf1t
    nf4 = nf4t
endif
print*,'nf1=',nf1, ' nf4=',nf4,'  df=',df
print*, "the maximum number of the roots",NNN_max

read (12, '(a80)' ) list
read (12, '(a80)' ) list
read (12,    *    ) tolmin
read (12,    *    ) tolmax
read (12,    *    ) Smin_min
read (12,    *    ) Smin_max
read (12, '(a80)' ) list
read (12,    *    ) dc
read (12,    *    ) dcm
read (12,    *    ) dc1
read (12,    *    ) dc2
!read (12,    *    ) dc,dcm  ! for my code
!read (12,    *    ) dc1,dc2  ! for my code
!read (12, '(a80)' ) list
!read (12,    *    ) ileaky
close (12)
print*,'dc ','dcm ','dc1 ','dc2'
print*, dc, dcm, dc1, dc2
end subroutine ReadParas
end subroutine ReadInput

real*8 function  N_cf(c,w)
use hash,only:nly,pi,vp,vs,d
implicit none
real*8      c,w
complex*16  yp,ys,sum
integer i
sum=0.
do i=1,nly-1
    yp=sqrt(dcmplx((c/vp(i))**2-1))
    if(vs(i).NE.0.0) then
        ys=sqrt(dcmplx((c/vs(i))**2-1))
    else
        ys=0.0
    end if 
    sum=sum+2.0*(yp+ys)*d(i)/dcmplx(c) 
end do
N_cf=(w/(2.0*pi)*dble(sum))
end function N_cf

real*8 function  N_cf_L(c,w)
use hash,only:nly,pi,vs,d
implicit none
real*8      c,w
complex*16  ys,sum
integer i
sum=0.
do i=1,nly-1
    if(vs(i)>0.0) then
        ys=sqrt(dcmplx((c/vs(i))**2-1))
    else
        ys=0.0
    end if 
    sum=sum+2.0*(ys)*d(i)/dcmplx(c) 
end do 
N_cf_L=(w/(2.0*pi)*dble(sum))
end function N_cf_L
