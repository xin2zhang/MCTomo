! module of likelihood
! 
module m_likelihood

    use m_exception,       only : exception_raiseError
    use m_utils,           only : ii10, vs2vp_3d, vp2rho_3d
    use m_settings,        only : T_MOD
    use run_info,          only : T_RUN_INFO
    use like_settings,     only : T_DATA, T_LIKE_BASE, T_LIKE_SET, likeBaseSetup
    use cgal_delaunay,     only : d3
    use m_body_likelihood, only : body_likelihood, body_likelihood_grads, body_noise_likelihood, body_location_likelihood
    use m_surf_likelihood, only : surf_likelihood, surf_noise_likelihood, surf_likelihood_grads

    implicit none

    private

    public :: T_LIKE, like_setup
    public :: likelihood, noise_likelihood, location_likelihood
    public :: update_lgP_grads
    
    ! define a likelihood type based on the base type
    type T_LIKE
        real(kind=ii10) :: like, misfit, unweighted_misfit
        real(kind=ii10), dimension(:), allocatable :: grads
        type(T_LIKE_BASE), dimension(:), allocatable :: likelihoods
    end type

contains

    subroutine like_setup(like, dat, like_set, ncell_max)
        implicit none
        type(T_LIKE), intent(out) :: like
        type(T_DATA), dimension(:), intent(in) :: dat
        type(T_LIKE_SET), intent(in) :: like_set
        integer, intent(in) :: ncell_max

        integer i
        
        like%like = huge(like%like)
        like%misfit = huge(like%misfit)
        like%unweighted_misfit = huge(like%unweighted_misfit)
        
        allocate(like%grads(2*ncell_max+4*dat(1)%nsrc))
        if(like_set%datatype==3)then
            allocate(like%likelihoods(2))
        else
            allocate(like%likelihoods(1))
        endif

        do i = 1, size(dat)
            call likeBaseSetup(like%likelihoods(i),dat(i),like_set,ncell_max)
        enddo

    end subroutine like_setup

    subroutine likelihood(dat,model,RTI,perturbed_box,like_set,like)
        implicit none
        type(T_DATA), dimension(:), intent(in)      :: dat
        type(T_MOD), intent(inout)                  :: model
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(d3), dimension(2), intent(in)          :: perturbed_box
        type(T_LIKE_SET), intent(in)                :: like_set
        type(T_LIKE), intent(inout)                 :: like

        select case (like_set%datatype)
        case (0,1)
            call body_likelihood(dat(1),model,RTI,perturbed_box,like_set,like%likelihoods(1))
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
        case (2)
            ! only surface waves are used, set the vp according to vs and
            ! density according to vp
            call vs2vp_3d(model%vs,model%vp)
            call vp2rho_3d(model%vp,model%rho)
            call surf_likelihood(dat(1),model,RTI,perturbed_box,like_set,like%likelihoods(1))
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
        case (3)
            call body_likelihood(dat(1),model,RTI,perturbed_box,like_set,like%likelihoods(1))
            call surf_likelihood(dat(2),model,RTI,perturbed_box,like_set,like%likelihoods(2))
            like%like = like%likelihoods(1)%like + like%likelihoods(2)%like
            like%misfit = like%likelihoods(1)%misfit + like%likelihoods(2)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit + like%likelihoods(2)%unweighted_misfit
        end select

        !write(*,*) '-loglikelihood: ', like%like
        !write(*,*) 'misfits: ', like%misfit
        !write(*,*) 'unweighted misfits: ', like%unweighted_misfit

    end subroutine likelihood

    subroutine update_lgP_grads(dat,model,RTI,like_set,like)
        implicit none
        type(T_DATA), dimension(:), intent(in) :: dat
        type(T_MOD), intent(in)                :: model
        type(T_RUN_INFO), intent(inout) :: RTI
        type(T_LIKE_SET), intent(in) :: like_set
        type(T_LIKE), intent(inout) :: like

        select case (like_set%datatype)
        case (0,1)
            call body_likelihood_grads(dat(1),model,RTI,like_set,like%likelihoods(1))
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
            like%grads = like%likelihoods(1)%grads
        case (2)
            ! only surface waves are used, set the vp according to vs and
            ! density according to vp
            !call vs2vp_3d(model%vs,model%vp)
            !call vp2rho_3d(model%vp,model%rho)
            call surf_likelihood_grads(dat(1),RTI,like_set,like%likelihoods(1))
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
            like%grads(1:RTI%ncells) = like%likelihoods(1)%grads(1:RTI%ncells)
        case (3)
            call body_likelihood_grads(dat(1),model,RTI,like_set,like%likelihoods(1))
            call surf_likelihood_grads(dat(2),RTI,like_set,like%likelihoods(2))
            like%like = like%likelihoods(1)%like + like%likelihoods(2)%like
            like%misfit = like%likelihoods(1)%misfit + like%likelihoods(2)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit + like%likelihoods(2)%unweighted_misfit
            like%grads = like%likelihoods(1)%grads
            like%grads(RTI%ncells+1:2*RTI%ncells) = like%likelihoods(1)%grads(RTI%ncells+1:2*RTI%ncells) + like%likelihoods(2)%grads(1:RTI%ncells)
        end select

    endsubroutine update_lgP_grads

    subroutine noise_likelihood(dat,RTI,like_set,like)
        implicit none
        type(T_DATA), dimension(:), intent(in)      :: dat
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(T_LIKE_SET), intent(in)                :: like_set
        type(T_LIKE), intent(inout)                 :: like

        select case (like_set%datatype)
        case (0,1)
            call body_noise_likelihood(dat(1),RTI,like_set,like%likelihoods(1))
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
        case (2)
            call surf_noise_likelihood(dat(1),RTI,like%likelihoods(1))
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
        case (3)
            call body_noise_likelihood(dat(1),RTI,like_set,like%likelihoods(1))
            call surf_noise_likelihood(dat(2),RTI,like%likelihoods(2))
            like%like = like%likelihoods(1)%like + like%likelihoods(2)%like
            like%misfit = like%likelihoods(1)%misfit + like%likelihoods(2)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit + like%likelihoods(2)%unweighted_misfit
        end select

        !write(*,*) '-loglikelihood: ', like%like
        !write(*,*) 'misfits: ', like%misfit
        !write(*,*) 'unweighted misfits: ', like%unweighted_misfit
    end subroutine noise_likelihood

    
    subroutine location_likelihood(dat,model,RTI,perturbed_box,like_set,like)
        implicit none
        type(T_DATA), dimension(:), intent(in)      :: dat
        type(T_MOD), intent(inout)                  :: model
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(d3), dimension(2), intent(in)          :: perturbed_box
        type(T_LIKE_SET), intent(in)                :: like_set
        type(T_LIKE), intent(inout)                 :: like

        if(like_set%datatype == 2) &
            call exception_raiseError('Error. Surface wave tomography cannot change source locations!')

        call body_location_likelihood(dat(1),model,RTI,perturbed_box,like_set,like%likelihoods(1))

        select case (like_set%datatype)
        case(0,1)
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
        case(3)
            like%like = like%likelihoods(1)%like + like%likelihoods(2)%like
            like%misfit = like%likelihoods(1)%misfit + like%likelihoods(2)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit + like%likelihoods(2)%unweighted_misfit
        end select

        !write(*,*) '-loglikelihood: ', like%like
        !write(*,*) 'misfits: ', like%misfit
        !write(*,*) 'unweighted misfits: ', like%unweighted_misfit
    end subroutine location_likelihood

end module m_likelihood
