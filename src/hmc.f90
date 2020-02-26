! Hamiltonial monte carlo
!
module m_hmc
    use m_logger,      only : log_msg
    use m_utils,       only : ii10, rtoa
    use like_settings, only : T_DATA, T_LIKE_SET
    use m_likelihood,  only : T_LIKE, update_lgP_grads
    use run_info,      only : T_RUN_INFO
    use mt19937,       only : grnd, gasdev

    implicit none

    public :: hmc_step

    type T_HMC_SET
        ! number of leapforg steps
        integer n_leapfrogs 
        ! step size of leapfrog
        real(kind=ii10), dimension(:), allocatable :: stepsize 
        ! preconditioned matrix for momentum
        real(kind=ii10), dimension(:), allocatable :: diagonal_matrix

        ! prior
        integer priortype
        real(kind=ii10), dimension(:), allocatable :: pm1, pm2

    endtype

    real(kind=ii10), parameter :: threshold = -0.5


contains

    subroutine hmc_step(target_data, RTI, like_set, hmc_set, qvals_in, lgP,accept) 
        implicit none
        type(T_DATA), dimension(:), intent(in)      :: target_data
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(T_LIKE_SET), intent(in)                :: like_set
        type(T_HMC_SET), intent(in)                 :: hmc_set
        real(kind=ii10), dimension(:), intent(in)   :: qvals_in
        !real(kind=ii10), dimension(:), intent(out)  :: qvals_out
        type(T_LIKE), intent(inout)                 :: lgP
        !type(T_LIKE), intent(in)                    :: lgP_in
        logical, intent(out)                        :: accept

        integer ndim, i
        real(kind=ii10), dimension(size(qvals_in)) :: pvals, qvals, qvals_start, grad
        real(kind=ii10), dimension(size(qvals_in)) :: diagonal_matrix, sqrt_diagonal_matrix, inv_diagonal_matrix
        real(kind=ii10) prev_K, prev_U, prop_K, prop_U
        real(kind=ii10) alpha, random
        !type(T_LIKE) :: lgP

        accept = .false.

        ! preconditioned diagonal matrix
        ndim = size(qvals_in)
        diagonal_matrix = hmc_set%diagonal_matrix(1:ndim)
        sqrt_diagonal_matrix = sqrt(diagonal_matrix)
        inv_diagonal_matrix = 1.0/diagonal_matrix

        ! regenerate momentum(p) values
        do i = 1, ndim
            pvals(i) = gasdev()
            !pvals(i) = random_normal()
        enddo
        pvals = pvals*sqrt_diagonal_matrix
        qvals = qvals_in
        ! if constrained, transform first
        if(hmc_set%priortype==1) call transform(qvals_in, qvals, hmc_set%pm1(1:ndim), hmc_set%pm2(1:ndim))
        prev_K = dot_product(pvals,inv_diagonal_matrix*pvals) / 2.0
        qvals_start =  qvals

        call log_msg('Beginning likelihood: '// rtoa(lgP%like))
        call update_lgP_grads(target_data, RTI, like_set, lgP)
        prev_U = posterior(qvals, lgP%like, hmc_set)
        call gradient(qvals, grad, lgP%grads(1:ndim), hmc_set)

        call log_msg('Begin leapfrog steps: ')
        !call random_number(random)
        random = grnd() 
        !call log_msg('Begining U: '//rtoa(prev_U) // 'Beginning K: '//rtoa(prev_K))
        !write(*,*) 'Beginning likelihood', lgP%like
        !write(13,*) 'Beginning velocity', RTI%parameters(1,1:rti%ncells)
        !write(13,*) 'Beginning sources', RTI%locations
        ! leapfrog steps
        do i = 1, hmc_set%n_leapfrogs

            ! debug
            !write(*,*) 'Velocity grads', lgP%grads(1:2*RTI%ncells)
            !write(*,*) 'Velocity', qvals(1:RTI%ncells)
            !write(*,*) 'Velocity', RTI%parameters(1,1:RTI%ncells)
            !write(*,*) 'Velocity grads', grad(1:2*RTI%ncells)
            !write(*,*) 'location grads', grad(2*RTI%ncells:ubound(grad,1))
            
            pvals = pvals - hmc_set%stepsize(1:ndim) * grad / 2.0

            qvals = qvals + hmc_set%stepsize(1:ndim)*inv_diagonal_matrix*pvals

            ! update potential (-loglikelihood) and gradient
            call update_model(qvals, hmc_set, RTI, like_set)
            call update_lgP_grads(target_data, RTI, like_set, lgP)
            call gradient(qvals, grad, lgp%grads(1:ndim), hmc_set)

            pvals = pvals - hmc_set%stepsize(1:ndim) * grad / 2.0

            !call log_msg('likelihood: '// rtoa(lgP%like))
            !call log_msg('distance: '// rtoa(norm2(qvals-qvals_start)))
            !call log_msg('distance: '// rtoa(dot_product(qvals-qvals_start,pvals)))
            ! reject it if it is larger than threshold
            prop_U = posterior(qvals, lgP%like, hmc_set)
            prop_K = dot_product(pvals,inv_diagonal_matrix*pvals)/2.0
            if(-prop_U-prop_K+prev_U+prev_K<log(random)+threshold) exit
        enddo

        !prop_U = posterior(qvals, lgP%like, hmc_set)
        !prop_K = dot_product(pvals,inv_diagonal_matrix*pvals)/2.0
        
        !call log_msg('End U: '//rtoa(prop_U) // 'End K: '//rtoa(prop_K))
        call log_msg('End likelihood: '// rtoa(lgP%like))
        !write(13,*) 'End sources', RTI%locations
        !write(13,*) 'End velocity', RTI%parameters(1,1:rti%ncells)
        !if(.not.(prop_U)<huge(prop_U)) prop_U = huge(prop_U)
        alpha = minval([0.0_ii10, -prop_U - prop_K + prev_U + prev_K])
        call log_msg('alpha: '//rtoa(alpha)//' random:'//rtoa(random))
        if(log(random)<alpha)then
            !qvals_out = qvals
            accept = .true.
        endif

    end subroutine

    ! negtive log posterior 
    ! only needed up to a constant, this is for Gaussian prior 
    ! for uniform prior, transformation is needed
    function posterior(qvals_trans, lglike, hmc_set)
        implicit none
        real(kind=ii10), dimension(:), intent(in) :: qvals_trans
        real(kind=ii10), intent(in) :: lglike
        type(T_HMC_SET), intent(in) :: hmc_set
        real(kind=ii10) :: posterior

        ! note the lglike is alreadly negative log likelihood
        posterior = lglike
        select case(hmc_set%priortype)
        case(0)
            ! if gaussian prior
            posterior = posterior + sum((qvals_trans-hmc_set%pm1)**2/(2*hmc_set%pm2**2))
        case(1)
            ! if uniform prior, transformation needed
            posterior = posterior - log_jacobian(qvals_trans, hmc_set%pm1, hmc_set%pm2)
        case default
            ! not recognized prior type
        endselect

    endfunction posterior

    subroutine gradient(qvals_trans, grad, grad_in, hmc_set)
        implicit none
        real(kind=ii10), dimension(:), intent(in) :: qvals_trans
        real(kind=ii10), dimension(:), intent(out) :: grad
        real(kind=ii10), dimension(:), intent(in) :: grad_in
        type(T_HMC_SET), intent(in) :: hmc_set

        grad = grad_in
        select case(hmc_set%priortype)
        case(0)
            ! if gaussian prior
            grad = grad + (qvals_trans-hmc_set%pm1)/hmc_set%pm2**2
        case(1)
            ! if uniform prior, transformation needed
            call jacobian_trans_grad(grad, qvals_trans, hmc_set%pm1, hmc_set%pm2)
        case default
            ! not recognized prior type
        endselect

    endsubroutine gradient

    subroutine update_p(pvals, gradient, stepsize)
        implicit none
        real(kind=ii10), dimension(:), intent(inout) :: pvals
        real(kind=ii10), dimension(:), intent(inout) ::  gradient
        real(kind=ii10), dimension(:), intent(in) :: stepsize

        ! update p
        pvals = pvals - stepsize * gradient / 2.0

    endsubroutine

    subroutine update_model(qvals, hmc_set, rti, like_set)
        implicit none
        real(kind=ii10), dimension(:), intent(in) :: qvals
        type(T_HMC_SET), intent(in) :: hmc_set
        type(T_RUN_INFO), intent(inout) :: rti
        type(T_LIKE_SET), intent(in) :: like_set

        real(kind=ii10), dimension(size(qvals)) :: inv_qvals
        integer ndim, ncells, nsrc

        ndim = size(qvals)
        ncells = RTI%ncells
        inv_qvals = qvals
        nsrc = size(rti%locations,2)

        ! if constrained, transform back to real values first
        if(hmc_set%priortype==1)then
            call inv_transform(qvals, inv_qvals, hmc_set%pm1(1:ndim), hmc_set%pm2(1:ndim))
        endif

        select case(like_set%datatype)
        case (0)
            rti%parameters(1,1:ncells) = inv_qvals(1:ncells)
            ndim = ncells
        case (1,3)
            rti%parameters(1,1:ncells) = inv_qvals(1:ncells)
            rti%parameters(2,1:ncells) = inv_qvals(ncells+1:2*ncells)
            ndim = 2*ncells
        case (2)
            rti%parameters(2,1:ncells) = inv_qvals(1:ncells)
            ndim = ncells
        case default
            ! default
        endselect

        ! if source included
        if(size(qvals)>ndim)then
            rti%locations(1,:) = inv_qvals(ndim+1:ndim+4*nsrc:4)
            rti%locations(2,:) = inv_qvals(ndim+2:ndim+4*nsrc:4)
            rti%locations(3,:) = inv_qvals(ndim+3:ndim+4*nsrc:4)
            rti%locations(4,:) = inv_qvals(ndim+4:ndim+4*nsrc:4)
        endif
        
    endsubroutine update_model

    ! transform uniform prior into infinite space for hmc
    ! Aternative is to use "billiards" (Neal, 2010)
    subroutine transform(qvals, qvals_out, lower_bounds, upper_bounds)
        !use ieee_arithmetic
        implicit none
        real(kind=ii10), dimension(:), intent(in) :: qvals
        real(kind=ii10), dimension(:), intent(out) :: qvals_out
        real(kind=ii10), dimension(:), intent(in) :: lower_bounds
        real(kind=ii10), dimension(:), intent(in) :: upper_bounds

        integer i

        qvals_out = 0.
        do i = 1, size(qvals)
            if(ieee_is_finite(lower_bounds(i)) .and. .not.ieee_is_finite(upper_bounds(i)))then
                qvals_out(i) = log(qvals(i)-lower_bounds(i))
            elseif(.not.ieee_is_finite(lower_bounds(i)) .and. ieee_is_finite(upper_bounds(i)))then
                qvals_out(i) = log(upper_bounds(i)-qvals(i))
            else
                qvals_out(i) = log(qvals(i)-lower_bounds(i)) - log(upper_bounds(i)-qvals(i))
            endif
        enddo
           
    endsubroutine transform

    subroutine inv_transform(qvals, qvals_out, lower_bounds, upper_bounds)
        !use ieee_arithmetic
        implicit none
        real(kind=ii10), dimension(:), intent(in) :: qvals
        real(kind=ii10), dimension(:), intent(out) :: qvals_out
        real(kind=ii10), dimension(:), intent(in) :: lower_bounds
        real(kind=ii10), dimension(:), intent(in) :: upper_bounds

        integer i
        real(kind=ii10) exp_term

        qvals_out = 0.
        do i = 1, size(qvals)
            if(ieee_is_finite(lower_bounds(i)) .and. .not.ieee_is_finite(upper_bounds(i)))then
                qvals_out(i) = lower_bounds(i) + exp(qvals(i))
            elseif(.not.ieee_is_finite(lower_bounds(i)) .and. ieee_is_finite(upper_bounds(i)))then
                qvals_out(i) = upper_bounds(i) - exp(qvals(i))
            else
                if(qvals(i)>0)then
                    exp_term = exp(-qvals(i))
                    qvals_out(i) = lower_bounds(i) + (upper_bounds(i)-lower_bounds(i))&
                        / (1+exp_term)
                else
                    exp_term = exp(qvals(i))
                    qvals_out(i) = upper_bounds(i) + (lower_bounds(i) - upper_bounds(i))&
                        / (1+exp_term)
                endif
            endif
        enddo
           
    endsubroutine inv_transform

    function log_jacobian(qvals_trans, lower_bounds, upper_bounds)
        !use ieee_arithmetic
        implicit none
        real(kind=ii10), dimension(:), intent(in) :: qvals_trans
        real(kind=ii10), dimension(:), intent(in) :: lower_bounds
        real(kind=ii10), dimension(:), intent(in) :: upper_bounds
        real(kind=ii10) log_jacobian

        integer i

        log_jacobian = 0.
        do i = 1, size(qvals_trans)
            if(ieee_is_finite(lower_bounds(i)) .and. .not.ieee_is_finite(upper_bounds(i)))then
                log_jacobian = log_jacobian + qvals_trans(i)
            elseif(.not.ieee_is_finite(lower_bounds(i)) .and. ieee_is_finite(upper_bounds(i)))then
                log_jacobian = log_jacobian + qvals_trans(i)
            else
                log_jacobian = log_jacobian + log(upper_bounds(i)-lower_bounds(i)) + qvals_trans(i) - 2*log(1+exp(qvals_trans(i)))
            endif
        enddo

    endfunction log_jacobian

    ! calculate new gradient after transformation, note it is actually a diagonal matrix
    subroutine jacobian_trans_grad(grad, qvals_trans, lower_bounds, upper_bounds)
        !use ieee_arithmetic
        implicit none
        real(kind=ii10), dimension(:), intent(inout) :: grad
        real(kind=ii10), dimension(:), intent(in) :: qvals_trans
        real(kind=ii10), dimension(:), intent(in) :: lower_bounds
        real(kind=ii10), dimension(:), intent(in) :: upper_bounds

        integer i
        real(kind=ii10) exp_term

        do i = 1, size(qvals_trans)
            exp_term = exp(qvals_trans(i))
            if(ieee_is_finite(lower_bounds(i)) .and. .not.ieee_is_finite(upper_bounds(i)))then
                grad(i) = exp_term*grad(i) - 1
            elseif(.not.ieee_is_finite(lower_bounds(i)) .and. ieee_is_finite(upper_bounds(i)))then
                grad(i) = -exp_term*grad(i) - 1
            else
                grad(i) = exp_term*(upper_bounds(i)-lower_bounds(i))/(1+exp_term)**2 * grad(i) - &
                    1 + 2*exp_term/(1+exp_term)
            endif
        enddo

    endsubroutine jacobian_trans_grad

    function ieee_is_finite(num) result(isfinite)
        real(kind=ii10), intent(in) :: num
        logical :: isfinite

        isfinite=.false.
        if(num<huge(num)) isfinite=.true.
    end function


    function random_normal() result(k)

      implicit none
      real(kind=ii10) :: k
      double precision :: s, u, v

      do 
         call random_number(u)
         call random_number(v)
         u = u*2.0d0 - 1.0d0
         s = u * u + v * v
         if ((s < 1).and.(s > 0)) exit
      end do
      k =  u * sqrt(-2 * log(s) / s)
    end function

end module
