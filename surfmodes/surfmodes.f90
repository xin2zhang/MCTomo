! Subroutine to calculate surface modes for layered model which can
! contain fluid layers above the solid layers (e.g. sea floor).
!
 
module m_surfmodes
      
      use m_GRT, only : T_GRT, init_grt, dp
      use omp_lib

      implicit none

      private

      public :: surfmodes, surfmmodes
      public :: T_MODES_PARA

      ! static values
      real(kind=dp), parameter :: eps = 1E-6
      real(kind=dp), parameter :: pi = 3.1415926
      integer, parameter :: nmodes = 100
      integer, parameter :: nfreqs = 100
      ! parameters for surface wave modes calculation
      type T_MODES_PARA
          integer :: modetype, nmodes ! 1 for rayleigh wave, 0 for love wave
          integer :: phaseGroup ! 0 for phase velocity, 1 for group velocity
          real(kind=dp) :: tolmin, tolmax
          real(kind=dp) :: Smin_min, Smin_max
          real(kind=dp) :: dc, dcm
          real(kind=dp) :: dc1, dc2
      endtype

      ! debug
      real(kind=dp) :: t1, t2
      !$omp threadprivate(t1,t2)


contains

      subroutine surfmodes(thick,vp,vs,rho,freqs,paras,phase,group,ierr)
          implicit none
          real(kind=dp), dimension(:), intent(in) :: thick
          real(kind=dp), dimension(:), intent(in) :: vp
          real(kind=dp), dimension(:), intent(in) :: vs
          real(kind=dp), dimension(:), intent(in) :: rho
          real(kind=dp), dimension(:), intent(in) :: freqs
          type(T_MODES_PARA), intent(in) :: paras
          real(kind=dp), dimension(:), intent(out) :: phase
          real(kind=dp), dimension(:), intent(out) :: group
          integer, intent(out) :: ierr

          type(T_GRT) GRT
          integer nlayers

          nlayers = size(thick)

          ! prepare for the GRT method
          call init_grt(GRT, nlayers)
          GRT%d = thick
          GRT%vp = vp
          GRT%vs = vs
          GRT%rho = rho
          !t1=omp_get_wtime()
          call setup_grt(GRT,paras)
          !t2=omp_get_wtime()
          !write(*,*) 'Preparing time: ', t2-t1
          !open(22,file='1dmodel_tmp.dat',access='append')
          !write(22,*) nlayers
          !do i =  1, nlayers
          !    write(22,*) thick(i), vp(i), vs(i), rho(i)
          !enddo
          !close(22)

          !t1=omp_get_wtime()
          ierr = 0
          select case(paras%modetype)

          case(1)
              ! rayleigh waves
              if(grt%nlvls1==0)then
                  !write(*,*) 'No low velocity layer, using surfdisp96:'
                  call surfdisp96(real(thick,4),real(vp,4),real(vs,4),real(rho,4),&
                      size(thick),0,2,1,paras%phaseGroup,size(freqs),dble(1/freqs),&
                      paras%dc,phase,group,ierr)
              else
                  !write(*,*) 'low velocity layer, using generalized R/T:'
                  call RayleighModes(GRT, freqs, paras,&
                  phase,group,ierr)
              endif
          case(0)
              ! love waves
              if(grt%nlvls1==0)then
                  !write(*,*) 'No low velocity layer, using surfdisp96:'
                  call surfdisp96(real(thick,4),real(vp,4),real(vs,4),real(rho,4),&
                      size(thick),0,1,1,paras%phaseGroup,size(freqs),dble(1/freqs),&
                      paras%dc,phase,group,ierr)
              else
                  !write(*,*) 'low velocity layer, using generalized R/T:'
                  call&
                      LoveModes(GRT,freqs,paras,phase,group,ierr)
              endif
          case default
              ! 
          end select
          !t2=omp_get_wtime()
          !write(*,*) 'dispersion curve calculation: ', t2-t1


      end subroutine

      subroutine surfmmodes(thick,vp,vs,rho,freqs,paras,phase,group,ierr)
          implicit none
          real(kind=dp), dimension(:), intent(in) :: thick
          real(kind=dp), dimension(:), intent(in) :: vp
          real(kind=dp), dimension(:), intent(in) :: vs
          real(kind=dp), dimension(:), intent(in) :: rho
          real(kind=dp), dimension(:), intent(in) :: freqs
          type(T_MODES_PARA), intent(in) :: paras
          real(kind=dp), dimension(:), intent(out) :: phase
          real(kind=dp), dimension(:), intent(out) :: group
          integer, intent(out) :: ierr

          type(T_GRT) GRT
          integer nlayers, nfreqs
          real(kind=dp), dimension(size(freqs),paras%nmodes) :: phase2d, group2d

          nlayers = size(thick)
          nfreqs =  size(freqs)

          ! prepare for the GRT method
          call init_grt(GRT, nlayers)
          GRT%d = thick
          GRT%vp = vp
          GRT%vs = vs
          GRT%rho = rho
          !t1=omp_get_wtime()
          call setup_grt(GRT,paras)
          !t2=omp_get_wtime()
          !write(*,*) 'Preparing time: ', t2-t1
          !open(22,file='1dmodel_tmp.dat',access='append')
          !write(22,*) nlayers
          !do i =  1, nlayers
          !    write(22,*) thick(i), vp(i), vs(i), rho(i)
          !enddo
          !close(22)

          !t1=omp_get_wtime()
          ierr = 0
          select case(paras%modetype)

          case(1)
              ! rayleigh waves
              if(grt%nlvls1==0)then
                  !write(*,*) 'No low velocity layer, using surfdisp96:'
                  call surfdisp_mmodes(real(thick,4),real(vp,4),real(vs,4),real(rho,4),&
                      size(thick),0,2,paras%nmodes,paras%phaseGroup,size(freqs),dble(1/freqs),&
                      paras%dc,phase2d,group2d,ierr)
              else
                  write(*,*) 'low velocity layer, not supported yet:'
                  !call RayleighModes(GRT, freqs, paras,&
                  !phase,group,ierr)
              endif
          case(0)
              ! love waves
              if(grt%nlvls1==0)then
                  !write(*,*) 'No low velocity layer, using surfdisp96:'
                  call surfdisp_mmodes(real(thick,4),real(vp,4),real(vs,4),real(rho,4),&
                      size(thick),0,1,paras%nmodes,paras%phaseGroup,size(freqs),dble(1/freqs),&
                      paras%dc,phase2d,group2d,ierr)
              else
                  write(*,*) 'low velocity layer, not supported yet:'
                  !call&
                  !    LoveModes(GRT,freqs,paras,phase,group,ierr)
              endif
          case default
              ! 
          end select
          !t2=omp_get_wtime()
          !write(*,*) 'dispersion curve calculation: ', t2-t1
          phase = reshape(phase2d,(/nfreqs*paras%nmodes/))
          group = reshape(group2d,(/nfreqs*paras%nmodes/))


      end subroutine

      subroutine RayleighModes(GRT,freqs,paras,&
              phase,group,ierr,allroots)
          implicit none
          type(T_GRT), intent(inout) :: GRT
          real(kind=dp), dimension(:), intent(in) :: freqs
          type(T_MODES_PARA), intent(in) :: paras
          real(kind=dp), dimension(:), intent(out) :: phase
          real(kind=dp), dimension(:), intent(out) :: group
          integer, intent(out) :: ierr
          real(kind=dp), dimension(:,:), intent(out), optional :: allroots

          real(kind=dp), parameter :: dh = 0.005
          real(kind=dp) freq0, cp0
          real(kind=dp), dimension(nmodes,nfreqs) :: roots, roots0
          real(kind=dp) c0
          integer allmodes, nfreqs
          integer nroots, ierr1
          integer i

          allmodes = 0
          if(present(allroots))then
              allmodes = 1
          endif
          
          ierr = 0
          roots = 0
          c0 = 0
          nfreqs = size(freqs)
          do i = 1, size(freqs)
              grt%w = freqs(i)*2*pi
              grt%tol = paras%tolmin + (nfreqs+1-i)*(paras%tolmax-paras%tolmin)/nfreqs
              grt%smin = paras%smin_min + (i-1)*(paras%smin_max-paras%smin_min)/nfreqs
              grt%index_a = i
              call SearchRayleigh(GRT, c0, &
                  roots(:,i),nroots, allmodes,ierr1)
              if(ierr1 == 1)then
                  ierr = 1
                  exit
              endif
              phase(i) = roots(1,i)
              c0 = phase(i)
              grt%ncr1 = nroots
              grt%root1 = roots(:,i)
              if(paras%phaseGroup==1)then
                  freq0 = freqs(i) + dh
                  grt%w = freq0*2*pi
                  call SearchRayleigh(GRT,c0,&
                      roots0(:,i),nroots, allmodes,ierr)
                  if(ierr==1) exit
                  cp0 = roots0(1,i)
                  call CalGroup(phase(i),cp0,freqs(i),dh,group(i))
              endif
          enddo

          if(present(allroots))then
              allroots = roots(1:size(allroots,1),1:size(allroots,2))
          endif

      end subroutine

      subroutine LoveModes(GRT,freqs,paras,&
              phase,group,ierr,allroots)
          implicit none
          type(T_GRT), intent(inout) :: GRT
          real(kind=dp), dimension(:), intent(in) :: freqs
          type(T_MODES_PARA), intent(in) :: paras
          real(kind=dp), dimension(:), intent(out) :: phase
          real(kind=dp), dimension(:), intent(out) :: group
          integer, intent(out) :: ierr
          real(kind=dp), dimension(:,:), intent(out), optional :: allroots

          real(kind=dp), parameter :: dh = 0.005
          real(kind=dp) freq0, cp0
          real(kind=dp), dimension(nmodes,nfreqs) :: roots, roots0
          real(kind=dp) c0
          integer allmodes, nfreqs
          integer nroots, ierr1
          integer i

          allmodes = 0
          if(present(allroots))then
              allmodes = 1
          endif
          
          nfreqs = size(freqs)
          ierr = 0
          roots = 0
          c0 = 0
          do i = 1, size(freqs)
              grt%w = freqs(i)*2*pi
              grt%tol = paras%tolmin + (nfreqs+1-i)*(paras%tolmax-paras%tolmin)/nfreqs
              grt%smin = paras%smin_min + (i-1)*(paras%smin_max-paras%smin_min)/nfreqs
              grt%index_a = i
              call SearchLove(GRT,c0,&
                  roots(:,i),nroots,allmodes,ierr1)
              if(ierr1 == 1)then
                  ierr = 1
                  exit
              endif
              phase(i) = roots(1,i)
              c0 = phase(i)
              grt%ncr1 = nroots
              grt%root1 = roots(:,i)
              if(paras%phaseGroup==1)then
                  freq0 = freqs(i) + dh
                  grt%w = freq0*2*pi
                  call SearchLove(GRT,c0,&
                      roots0(:,i),nroots,allmodes,ierr)
                  if(ierr==1) exit
                  cp0 = roots0(1,i)
                  call CalGroup(phase(i),cp0,freqs(i),dh,group(i))
              endif
          enddo

          if(present(allroots))then
              allroots = roots(1:size(allroots,1),1:size(allroots,2))
          endif

      end subroutine

      subroutine CalGroup(c1,c2,freq,dh,group)
          implicit none
          real(kind=dp), intent(in) :: c1, c2
          real(kind=dp), intent(in) :: freq
          real(kind=dp), intent(in) :: dh
          real(kind=dp), intent(out) :: group

          group = (freq+dh)/c2 - freq/c1
          if(group>0)then
              group = dh/group
          else
              group = 0
          endif
      endsubroutine

      subroutine setup_grt(GRT, para)
          implicit none
          type(T_GRT), intent(inout) :: GRT
          type(T_MODES_PARA), intent(in) :: para

          integer idx, i, j, k
          integer nlayers
          real(kind=dp) mu0, mu

          GRT%dc = para%dc
          GRT%dc2 = para%dc2
          GRT%dcm = para%dcm
          GRT%ifs = 0

          idx = 0
          nlayers = GRT%nlayers
          do i = 1, nlayers
              if(abs(GRT%vs(i))>eps) then
                  idx = idx + 1
                  GRT%v(idx) = GRT%vs(i)
                  if(i==nlayers) GRT%ilastvs = idx
              else
                  if(i>1)then
                      write(*,*) 'Error: currently only allow first layer to be water'
                      stop
                  endif
                  GRT%ifs = GRT%ifs + 1
              endif
              idx = idx + 1
              GRT%v(idx) = GRT%vp(i)
          enddo
          ! sort the wave velocity
          call sort(GRT%v,idx,1)

          ! prepare mu
          idx = 0
          mu0 = 0
          do i = 1, nlayers
              mu = GRT%rho(i)*GRT%vs(i)**2
              if(abs(mu)>eps) then
                  idx = idx + 1
                  mu0 = mu0 + mu
              endif
              GRT%mu(i) = mu
          end do
          mu0 = mu0/idx
          GRT%mu =GRT%mu/mu0
          GRT%mu0 = mu0

          GRT%vsy = maxval(GRT%vs)
          select case(para%modetype)
          case(1)
              if(GRT%ifs>0) then
                  GRT%vs1=GRT%vp(1)
                  GRT%vss1=GRT%vs(grt%ifs+1)
              else
                  GRT%vs1=GRT%vs(1)
              endif
              GRT%vsm=GRT%v(1)
          case(0)
              GRT%vs1=GRT%vs(1+GRT%ifs)
              GRT%vsm=minval(GRT%vs(GRT%ifs+1:nlayers))! the lowest S wave velocity
          end select
          !vsm=vs1 
          ! number of lvls in fluid
          GRT%no_lvl_fl=0
          do i=2,nlayers-1
              if(i>grt%ifs .and. grt%vs(i)<grt%vss1) grt%nlvls1 = grt%nlvls1 + 1
              if(GRT%vp(i)<GRT%vp(i+1) .and. GRT%vp(i)<GRT%vp(i-1)) then
                  !print*,'LVL:',i
                  GRT%no_lvl=GRT%no_lvl+1
                  GRT%lvls(GRT%no_lvl)=i
          !       if(i>ifs+1) no_lvl_L=no_lvl_L+1
          !       if(vs(i)<vsm) vsm=vs(i)
                  if(GRT%ifs==0) then
                      if(GRT%vs(i)<GRT%vs1) GRT%nlvl1=GRT%nlvl1+1
                  else
                      select case(para%modetype)
                      case(1)
                          if(GRT%vp(i)<GRT%vs1) grt%nlvl1=grt%nlvl1+1
                      case(0)
                          if(grt%vs(i)>0.) then
                              if(grt%vs(i)<grt%vs1) grt%nlvl1=grt%nlvl1+1
                          else
                              grt%no_lvl_fl=grt%no_lvl_fl+1
                          endif
                      end select
                  endif
              endif
          enddo
          if(grt%ifs==0 .or. para%modetype==0) grt%nlvls1 = grt%nlvl1
          !
          select case(para%modetype)
          case(1)
              grt%lvls(grt%no_lvl+1)=1
          case(0)
              grt%lvls(grt%no_lvl+1)=1+grt%ifs
          end select

          ! prepare lvlast
          if(grt%no_lvl==0) then
              grt%lvlast=1
              !print*,'LVLs are absent.'
          else
              grt%lvlast=grt%lvls(grt%no_lvl)
              do i=1,grt%no_lvl-1
                  do j=1,grt%no_lvl-i
                      if(grt%vp(grt%lvls(j))>grt%vp(grt%lvls(j+1))) then
                          k=grt%lvls(j)
                          grt%lvls(j)=grt%lvls(j+1)
                          grt%lvls(j+1)=k
                      endif
                  enddo
              enddo
              !print '(2(a,i2))','nlvl1=',grt%nlvl1,', no_lvl=',grt%no_lvl
              !print '(a,(/),1(i4,f7.3,1x,f7.3))', 'LVLs:', &
              !    (grt%lvls(i),grt%vs(grt%lvls(i)),grt%vp(grt%lvls(i)),i=1,grt%no_lvl)
              !print '(a,f7.3)', 'vs1:',grt%vs1
          endif
          grt%lvlast=max(grt%ifs+2,grt%lvlast)

          ! lvls has been sorted
          do i=1,grt%no_lvl
              if(grt%lvls(i)>grt%ifs) then
                  if(grt%vs(grt%lvls(i))<grt%vs1) then
                      grt%L1=i
                      exit
                  endif
              endif
          enddo
      end subroutine

end module
