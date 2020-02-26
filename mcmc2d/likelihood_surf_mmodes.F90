!
! module m_likelihood
! likelihood function here
!
module m_surf_likelihood

    use m_exception, only : exception_raiseWarning, exception_raiseError
    use m_logger, only    : log_msg
    use m_utils, only     : ii10, write_resume_unit, write_doubles, itoa, rtoa, vs2vp, vp2rho
    use m_settings, only  : T_GRID, T_MOD, mod_setup
    use like_settings, only: T_LIKE_BASE, T_LIKE_SET, T_DATA
    use run_info, only    : T_RUN_INFO
    use m_fm2d, only      : modrays, T_RAY
    use m_surfmodes, only : surfmmodes, T_MODES_PARA

    use omp_lib

    implicit none
    private

    public :: surf_likelihood, surf_noise_likelihood
    public :: surf_likelihood_grads

    ! debug
    real( kind=ii10 ) :: t1, t2
    ! static value
    real( kind=ii10 ), parameter :: ERAD = 6371
    ! maximum layers
    integer, parameter :: NMAX = 100
    ! water P wave velocity and density
    real( kind=ii10 ), parameter :: waterVel = 1.5
    real( kind=ii10 ), parameter :: waterDensity = 1
    ! vp/vs ratio, currently, used fixed value 1.73
    real( kind=ii10 ), parameter :: POISSON = 1.730
    real( kind=ii10 ), parameter :: PI2 = 6.283185
    real( kind=ii10 ), parameter :: EPS = 1.0E-10
    ! temporarily use
    integer, parameter :: expand = 1
    integer, parameter :: ext = 2

    type T_LAY_MOD
        integer nlayers
        real(kind=ii10), dimension(:), allocatable :: alpha
        real(kind=ii10), dimension(:), allocatable :: beta
        real(kind=ii10), dimension(:), allocatable :: vpvs
        real(kind=ii10), dimension(:), allocatable :: rho
        real(kind=ii10), dimension(:), allocatable :: thick
        real(kind=ii10), dimension(:), allocatable :: ax
        real(kind=ii10), dimension(:), allocatable :: ap 
        real(kind=ii10), dimension(:), allocatable :: ae
    end type T_LAY_MOD

contains

    subroutine surf_noise_likelihood(dat,RTI,like)
        implicit none
        type(T_DATA), intent(in)                 :: dat
        type(T_RUN_INFO), intent(in)             :: RTI
        type(T_LIKE_BASE), intent(inout)         :: like

        integer nrr, nmodes
        integer nsrc, nrev
        integer i, j, k

        nsrc = dat%nsrc
        nrev = dat%nrev
        nmodes = size(dat%raystat,3)

        do j = 1, nmodes
            do i = 1, dat%np
                k = i + (j-1)*dat%np
                like%sigma(:,i,j) = RTI%snoise0(k)*like%srdist(:,i,j) + RTI%snoise1(k)
            enddo
        enddo

        like%like = 0
        like%misfit = 0
        do k = 1, dat%nmodes
            do i = 1, dat%nrays
                do j = 1, dat%np
                    if(dat%raystat(j,i,k) == 1) then
                        like%like = like%like + ( like%phaseTime(j,i,k)-dat%ttime(j,i,2*k-1)&
                            )**2/( 2*(like%sigma(j,i,k))**2 )
                        like%misfit = like%misfit + ( like%phaseTime(j,i,k)-dat%ttime(j,i,2*k-1)&
                            )**2/(like%sigma(j,i,k)**2)
                    else
                        like%sigma(j,i,k) = 1.0
                    endif
                enddo
            enddo
        enddo

        if(any(like%sigma<EPS))then
            call exception_raiseError('The noise level is 0!')
        endif

        like%like = like%like + sum(log(like%sigma)) + sum(dat%raystat)/2.0 * log(PI2)
        return
    end subroutine surf_noise_likelihood

    subroutine surf_likelihood_grads(dat,RTI,settings,like)
        implicit none
        type(T_DATA), intent(in)                    :: dat
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(T_LIKE_SET), intent(in)                :: settings
        type(T_LIKE_BASE), intent(inout)            :: like

    endsubroutine surf_likelihood_grads

    subroutine surf_likelihood(dat,RTI,settings,like)
    
        implicit none
        type(T_DATA), intent(in)                    :: dat
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(T_LIKE_SET), intent(in)                :: settings
        type(T_LIKE_BASE), intent(inout)            :: like

        ! > local variable for normal mode modeling
        type(T_MODES_PARA) :: paras
        integer         :: ix0, ix1
        integer         :: idx0, idx1
        !integer         :: iper, ier, ia
        !real(kind=ii10) :: w, ekd, y0l(6), y0r(3), yij(15)
        integer, dimension(:), allocatable            :: ierr
        real(kind=ii10), dimension(:,:,:), allocatable :: pvel, gvel
        type(T_LAY_MOD), dimension(:), allocatable    :: layer

        ! > local variable for travel times and rays
        integer nsrc, nrev

        ! > local grid
        type(T_GRID) :: grid
        type(T_MOD) :: model

        integer nrr, iflat
        integer i, j, k

        nsrc = dat%nsrc
        nrev = dat%nrev
        grid = settings%grid

        if(settings%sigdep /= 0)then
            do j = 1, dat%nmodes
                do i = 1, dat%np
                    k = i + (j-1)*dat%np
                    like%sigma(:,i,j) = RTI%snoise0(k)*like%srdist(:,i,j) + RTI%snoise1(k)
                enddo
            enddo
        endif
        
        call mod_setup(model,settings%grid)
        do i = 1, settings%grid%nx
            do j = 1, settings%grid%ny
                model%vp(j,i) = RTI%parameters(1,RTI%sites_id(j,i))
                model%vs(j,i) = RTI%parameters(2,RTI%sites_id(j,i))
                model%rho(j,i) = RTI%parameters(3,RTI%sites_id(j,i))
            enddo
        enddo

        if(check_model(model%vs))then
            like%like = huge(like%like)
            return
        endif

        ! first, calculate dispersion curve grid by grid
        ! convert 3d grid vp, vs and rho model to layered model
        ix0 = 1
        ix1 = grid%nx
        call convert_to_layer( model, grid, ix0, ix1, layer )
        !call write_layers(layer,'layer.dat')

        ! first set up the parameters for surface modes code
        paras%modetype = settings%raylov
        paras%phaseGroup = settings%phaseGroup
        paras%tolmin = settings%tol
        paras%tolmax = 10*settings%tol
        paras%smin_min = 1E-3
        paras%smin_max = 5E-3
        paras%dc = settings%dPhaseVel
        paras%dcm = settings%dPhaseVel
        paras%dc1 = settings%dPhaseVel
        paras%dc2 = settings%dPhaseVel
        paras%nmode = dat%nmodes
        ! calculate dispersion curve
        ! TODO: currently, discard any velocity model which could not produce
        ! surface waves in one or more frequencies
        allocate( pvel(dat%np, ix0:ix1, dat%nmodes) )
        allocate( gvel(dat%np, ix0:ix1, dat%nmodes) )
        pvel = 100.0 ! safe
        gvel = 100.0 ! safe
        allocate( ierr(ix0:ix1) )
        ierr = 0
#ifdef _OPENMP
        t1 = omp_get_wtime ( )
        !call log_msg('Begin normal modes...')
        !call omp_set_num_threads(settings%nthreads)
#endif
        !$omp parallel
        !$omp do private(i,j)
        do i = ix0, ix1
           call surfmmodes(layer(i)%thick,layer(i)%alpha,layer(i)%beta, &
               layer(i)%rho,dat%freqs,paras,pvel(:,i,:),gvel(:,i,:),ierr(i))
        enddo
        !$omp end do
        !$omp end parallel
        t2 = omp_get_wtime ( )
        !call log_msg(itoa(omp_get_num_threads() ) )
        !call log_msg('parallelized dispersion curve code: ' //rtoa(t2-t1) )
        !call write_vel(pvel,'phaseVel.dat')
        ! discard models
        if(any(ierr==1))then
            like%like = huge(like%like)
            RTI%num_bad_model = RTI%num_bad_model + 1
            return
        endif

        if(settings%phaseGroup==0)then
            like%phaseTime=pvel
        else
            like%phaseTime=gvel
        endif

        ! likelihood
        ! noise level
        if(settings%sigdep /= 0)then
            do k = 1, dat%nmodes
                do i = 1, dat%nrays
                    do j = 1, dat%np
                        iflat = j + (k-1)*dat%np
                        if(dat%raystat(j,i,k) == 1) then
                            like%sigma(j,i,k) = RTI%snoise0(iflat)*like%srdist(j,i,k) + RTI%snoise1(iflat) 
                        else
                            like%sigma(j,i,k) = 1.0
                        endif
                    enddo
                enddo
            enddo
        else
            like%sigma = dat%ttime(:,:,2:2*dat%nmodes:2)
        endif

        !if(any(like%sigma==0))then
        !    call exception_raiseError('The noise level is 0!')
        !endif

        like%like = 0
        like%misfit = 0
        like%unweighted_misfit = 0
        do k = 1, dat%nmodes
            do i = 1, dat%nrays
                do j = 1, dat%np 
                    if(dat%raystat(j,i,k) == 1) then
                        if(like%sigma(j,i,k)<EPS)then
                            !print * , k, j, i, like%sigma(nrr,i), like%srdist(nrr,i)
                            call exception_raiseError('The noise level is 0!')
                        endif
                        like%like = like%like + ( like%phaseTime(j,i,k)-dat%ttime(j,i,2*k-1)&
                            )**2/( 2*(like%sigma(j,i,k))**2 )
                        like%misfit = like%misfit + ( like%phaseTime(j,i,k)-dat%ttime(j,i,2*k-1)&
                            )**2/(like%sigma(j,i,k)**2)
                        like%unweighted_misfit = like%unweighted_misfit + ( like%phaseTime(j,i,k)-dat%ttime(j,i,2*k-1)&
                            )**2
                    else
                        like%sigma(j,i,k) = 1.0
                    endif
                enddo
            enddo
        enddo

        like%like = like%like + sum(log(like%sigma)) + sum(dat%raystat)/2.0 * log(PI2)
        if(like%like/=like%like)then
            print *, like%sigma
        endif
        return

    end subroutine surf_likelihood

    subroutine convert_to_layer(model, grid, ix0, ix1, layer)
        implicit none

        type(T_MOD), intent(in) :: model
        type(T_GRID), intent(in) :: grid
        integer, intent(in) :: ix0, ix1
        type(T_LAY_MOD), dimension(:),allocatable, intent(out) :: layer

        ! > local variables
        integer i, k
        integer nlayers
        integer last_k
        real(kind=ii10) :: last_vp, last_vs, last_rho
        real(kind=ii10), dimension(NMAX) :: alpha, beta, rho_k, thick

        allocate( layer(ix0:ix1) )
        !$omp parallel
        !$omp do private(i,k,last_vp,last_vs,last_rho,last_k,nlayers,alpha)&
        !$omp& private(beta,rho_k,thick)
        do i = ix0, ix1
            if( grid%waterDepth>EPS )then
                nlayers = 1
                alpha(1) = waterVel
                beta(1) = 0
                rho_k(1) = waterDensity
                thick(1) = grid%waterDepth
            else
                nlayers = 0
            endif
            last_vp = model%vp(1,i)
            last_vs = model%vs(1,i)
            last_rho = model%rho(1,i)
            last_k = 1
            do k = 2, grid%ny
                if( abs(model%vs(k,i)-last_vs) > EPS ) then
                    nlayers = nlayers + 1
                    ! store the parameters value for this layer
                    alpha(nlayers) = last_vp
                    beta(nlayers) = last_vs
                    rho_k(nlayers) = last_rho
                    thick(nlayers) = (k-last_k) * grid%dz
                    ! update last_* values
                    last_vp = model%vp(k,i)
                    last_vs = model%vs(k,i)
                    last_rho = model%rho(k,i)
                    last_k = k
                endif
            enddo
            ! last 2 layers, last thick is 0
            nlayers = nlayers + 1
            if( last_k /= grid%ny) then
                alpha(nlayers) = model%vp(grid%ny,i)
                beta(nlayers) = model%vs(grid%ny,i)
                rho_k(nlayers) = model%rho(grid%ny,i)
                thick(nlayers) = (grid%ny-last_k) * grid%dz
                thick(nlayers) = 0
            else
                alpha(nlayers) = model%vp(grid%ny,i)
                beta(nlayers) = model%vs(grid%ny,i)
                rho_k(nlayers) = model%rho(grid%ny,i)
                thick(nlayers) = 0
            endif
            ! allocate and give value for layer model
            allocate( layer(i)%alpha(nlayers), layer(i)%beta(nlayers),&
            layer(i)%vpvs(nlayers), layer(i)%rho(nlayers),&
            layer(i)%thick(nlayers) )
            !allocate( layer(j,i)%ax(nlayers), layer(j,i)%ap(nlayers),&
            !layer(j,i)%ae(nlayers) )
            layer(i)%nlayers = nlayers
            layer(i)%alpha(1:nlayers) = alpha(1:nlayers)
            layer(i)%beta(1:nlayers) = beta(1:nlayers)
            layer(i)%vpvs = POISSON
            layer(i)%rho(1:nlayers) = rho_k(1:nlayers)
            layer(i)%thick(1:nlayers) = thick(1:nlayers)/grid%scaling
            ! waterDepth deos not scale, need to be recovered
            if(grid%waterDepth>0) layer(i)%thick(1)=grid%waterDepth

            !layer(j,i)%ax = 1
            !layer(j,i)%ap = 1
            !layer(j,i)%ae = 1
        enddo
        !$omp end do
        !$omp end parallel

        return

    end subroutine

    function check_model(vs) result(valid)
        real(kind=ii10), dimension(:,:), intent(in) :: vs
        logical valid
        integer i
        
        valid = .false.
        do i = 1, size(vs,2)
            if(any(vs(2:size(vs,1),i)<vs(1,i)))then
                valid = .true.
                return
            endif
        enddo
    
    endfunction

    subroutine write_layer(layer,filename)
        use m_utils, only : write_resume_unit
        implicit none
        character(len=*), intent(in) :: filename
        type(T_LAY_MOD), intent(in) :: layer

        integer i

        open(unit=write_resume_unit, file=filename, status='unknown')
        write(write_resume_unit,*) layer%nlayers
        do i = 1, layer%nlayers
            write(write_resume_unit,*) i,layer%thick(i)*1000,layer%beta(i)*1000,&
            layer%vpvs(i), layer%rho(i)
        enddo
        close(write_resume_unit)
        return
    end subroutine

    subroutine write_layers(layer,filename)
        use m_utils, only : write_resume_unit
        implicit none
        character(len=*), intent(in) :: filename
        type(T_LAY_MOD), dimension(:), intent(in) :: layer

        integer i, j, k

        open(unit=write_resume_unit, file=filename, status='unknown')
        do k = lbound(layer,1), ubound(layer,1)
            if(layer(k)%beta(2)>minval(layer(k)%beta(2:layer(k)%nlayers)))then
            write(write_resume_unit,*) layer(k)%nlayers
            do i = 1, layer(k)%nlayers
                write(write_resume_unit,*) i,layer(k)%thick(i)*1000,layer(k)%alpha(i)*1000,layer(k)%beta(i)*1000,&
                layer(k)%rho(i)
            enddo
            endif
        enddo
        close(write_resume_unit)
        return
    end subroutine

    subroutine write_vel(pvel,filename)
        use m_utils, only : write_resume_unit, write_doubles
        implicit none
        character(len=*), intent(in) :: filename
        real(kind=ii10), dimension(:,:), intent(in) :: pvel

        open(unit=write_resume_unit, file=filename, status='unknown',&
        position='append')
        call write_doubles(pvel)
        close(write_resume_unit)
        return
    end subroutine
end module m_surf_likelihood
