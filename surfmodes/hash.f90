module hash
    implicit none
    !real*8,allocatable::FMRay(:,:,:)
    !real*8,target::a44(4,4),b44(4,4)
    !complex*16,target,allocatable::Rdu(:,:,:),Rud(:,:,:),Td(:,:,:),Tu(:,:,:)
    complex*16,parameter::ai=(0d0,1d0)
    complex*16,parameter::IC=(1d0,0d0)
    real*8,allocatable::d(:),z(:),vs(:),vp(:),v(:),mu(:),rho(:)
    real*8,parameter::expo=46d0,eps=1d-10,edNN=.5d0,pi=3.1415926535897932d0
    integer,parameter::nmode=1000,nmax=3000
    integer nf1,nf4,nf1t,nf4t,NN,ncr,index_a,index0
    integer::i1,i2,NN0
    integer modetype
    integer,allocatable::lvls(:)
    
    ! ll:		the first calculation layer
    ! lvlast:	the last low velocity layer
    ! no_lvl:	the number of low velocity layers
    ! nlvl1:	the number of layers whose velocity < first S wave velocity
    ! lvls:		# of velocity layers in ascending order
    integer ifun,ilay,ll,jj,nly,lvl,lvlast,&
            nlvl1,no_lvl,no_lvl_fl,ileaky,im1,ifs,L1
    logical Lend
    character*20 suf
    real*8 ap,as,am,vs1,vsm,vsmin,vsy,imf, mu0
    real*8 vk,w,freq,df,ff1,ff4
    real*8 tol0,tol,smin, dc,dcm, dc1,dc2
    real*8 tolmin,tolmax,smin_min,smin_max
    !real*8 dc1min,dc1max,dc2min,dc2max
    real*8,target::cr(3,nmode)
    real*8,pointer::root2(:),cr1
    real*8 allroots(3000,nmode),root1(nmode),ccc(20000),NNN_max
    !equivalence(cr1,cr(1,1))
    !target cr
end module hash
