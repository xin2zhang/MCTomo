! Slice sampling
!
module sliceSample
    use m_logger,      only : log_msg
    use m_utils,       only : ii10, rtoa
    use m_settings,    only : T_MCMC_SET, T_MOD
    use like_settings, only : T_DATA, T_LIKE_SET
    use m_likelihood,  only : T_LIKE, likelihood
    use run_info,      only : T_RUN_INFO
    use mt19937,       only : grnd, gasdev
    use cgal_delaunay, only : d3

    implicit none

    public :: slice_sample

contains

    subroutine slice_sample(likelihoodfunc,target_data,model,RTI,perturbed_box,like_set,qvals,w,bnd,idx,lgP,n)
        implicit none

        ! These are needed by likelihood TODO: reconfigure the code to make
        ! likelihood calculation clean, not through the whole code
        interface
            subroutine likelihoodfunc(target_data,model,RTI,perturbed_box,like_set,lgP)
                use m_settings,    only : T_MOD
                use like_settings, only : T_DATA, T_LIKE_SET
                use m_likelihood,  only : T_LIKE
                use run_info,      only : T_RUN_INFO
                use cgal_delaunay, only : d3
                type(T_DATA), dimension(:), intent(in)      :: target_data
                type(T_MOD), intent(inout)                  :: model
                type(T_RUN_INFO), intent(inout)             :: RTI
                type(d3), dimension(2), intent(in)          :: perturbed_box
                type(T_LIKE_SET), intent(in)                :: like_set
                type(T_LIKE), intent(inout)                 :: lgP
            endsubroutine
        endinterface

        type(T_DATA), dimension(:), intent(in)      :: target_data
        type(T_MOD), intent(inout)                  :: model
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(d3), dimension(2), intent(in)          :: perturbed_box
        type(T_LIKE_SET), intent(in)                :: like_set
        real(kind=ii10), dimension(:), intent(inout):: qvals
        !> The initial width
        real(kind=ii10), dimension(:), intent(in)   :: w
        !> Constrained boundary
        real(kind=ii10), dimension(:,:), intent(in)   :: bnd
        ! index of row of parameters (i.e. vp, vs or rho)
        integer, intent(in)                         :: idx
        !> The Loglikelihood bound
        type(T_LIKE), intent(inout)            :: lgP
        !> The number of likelihood calls
        integer, intent(inout) :: n

        ! The upper bound
        real(kind=ii10)    :: R, lgR
        ! The lower bound
        real(kind=ii10)    :: L, lgL
        real(kind=ii10)    :: x0, x
        real(kind=ii10)    :: lgP0, slice_lgp, lgX


        real(kind=ii10) :: temp_random

        integer :: i_step

        !call log_msg('Begin likelihood: '//rtoa(lgP%like))
        !call log_msg('Begin vp: '//rtoa(qvals(idx)))
        !call update_model(qvals, RTI, model)
        !call likelihoodfunc(target_data,model,RTI,perturbed_box,like_set,lgP)
        !call log_msg('Begin likelihood: '//rtoa(lgP%like))
        ! define a slice
        lgP0 = lgP%like
        slice_lgp = lgP0 -  log(grnd())

        !call log_msg('slice like: '//rtoa(slice_lgp))
        ! Select initial start and end points
        x0 = qvals(idx)
        temp_random = grnd()
        L = x0 -   temp_random   * w(idx) 
        R = x0 + (1-temp_random) * w(idx)

        ! temporaly use
        RTI%randcount =  RTI%randcount + 2
        ! Calculate initial likelihoods
        qvals(idx) = L
        call update_model(qvals, RTI, model)
        call likelihoodfunc(target_data,model,RTI,perturbed_box,like_set,lgP)
        lgL = lgP%like
        qvals(idx) = R
        call update_model(qvals, RTI, model)
        call likelihoodfunc(target_data,model,RTI,perturbed_box,like_set,lgP)
        lgR = lgP%like
        !call log_msg('likeL: '//rtoa(lgL))
        !call log_msg('likeR: '//rtoa(lgR))

        ! expand R until it's outside the likelihoodfunc region
        n = n + 2
        i_step=0
        do while( lgR <= slice_lgp .and. lgR /= huge(lgR) )
            i_step=i_step+1
            R = x0 + w(idx) * i_step
            if(R>bnd(idx,2))then
                lgR = huge(lgR)
                R = bnd(idx,2)
                exit
            endif
            qvals(idx) = R
            call update_model(qvals, RTI, model)
            call likelihoodfunc(target_data,model,RTI,perturbed_box,like_set,lgP)
            lgR = lgP%like
            n = n + 1
            !call log_msg('vp: '//rtoa(R))
            !call log_msg('likelihoodfunc: '//rtoa(lgR))
        end do
        if(i_step>100) write(*,'(" too many R steps (",I10,")")') i_step

        ! expand L until it's outside the likelihoodfunc region
        i_step=0
        do while(lgL <= slice_lgP .and. lgL /= huge(lgL) )
            i_step=i_step+1
            L = x0 - w(idx) * i_step
            if(L<bnd(idx,1))then
                lgL = huge(lgL)
                L = bnd(idx,1)
                exit
            endif
            qvals(idx) = L
            call update_model(qvals, RTI, model)
            call likelihoodfunc(target_data,model,RTI,perturbed_box,like_set,lgP)
            lgL = lgP%like
            n = n + 1
            !call log_msg('vp: '//rtoa(L))
            !call log_msg('likelihoodfunc: '//rtoa(lgL))
        end do
        if(i_step>100) write(*,'(" too many L steps (",I10,")")') i_step

        ! Sample within this bound
        i_step=0
        x = find_positive_within(L,R)

        contains

        recursive function find_positive_within(L,R) result(x1)
            implicit none
            !> The upper bound
            real(kind=ii10), intent(inout)   :: R
            !> The lower bound
            real(kind=ii10), intent(inout)   :: L

            ! The output finish point
            real(kind=ii10)     :: x1

            real(kind=ii10) :: x0Rd
            real(kind=ii10) :: x0Ld

            i_step=i_step+1
            if (i_step>100) then
                write(*,'("Polychord Warning: Non deterministic loglikelihoodfunc")')
                lgX = huge(lgX)
                x1 = (L+R)/2
                return
            end if

            ! Find the distance between x0 and L 
            x0Ld= abs(x0-L)
            ! Find the distance between x0 and R 
            x0Rd= abs(x0-R)

            ! Draw a random point within L and R
            x1 = x0+ (grnd() * (x0Rd+x0Ld) - x0Ld)

            ! temporaly use
            RTI%randcount =  RTI%randcount + 2
            ! calculate the likelihoodfunc 
            qvals(idx) = x1
            call update_model(qvals, RTI, model)
            call likelihoodfunc(target_data,model,RTI,perturbed_box,like_set,lgP)
            lgX = lgP%like
            n =  n + 1

            ! If we're not within the likelihoodfunc bound then we need to sample further
            if( lgX > slice_lgP .or. lgX == huge(lgX) ) then

                if ( x1-x0 > 0d0 ) then
                    ! If x1 is on the R side of x0, then
                    ! contract R
                    R = x1
                else
                    ! If x1 is on the L side of x0, then
                    ! contract L
                    L = x1
                end if

                ! Call the function again
                x1 = find_positive_within(L,R)

            end if
            ! otherwise x1 is returned

        end function find_positive_within


    end subroutine slice_sample

    subroutine update_model(qvals, RTI,model)
        implicit none
        real(kind=ii10), dimension(:), intent(in) :: qvals
        type(T_RUN_INFO), intent(inout) :: RTI
        type(T_MOD), intent(inout) :: model

        integer i, j, k
        integer ndim, ncells, nsrc
        
        ndim = size(qvals)
        ncells = RTI%ncells
        nsrc = size(rti%locations,2)
        ! update model
        rti%parameters(1,1:ncells) = qvals(1:ncells)
        rti%parameters(2,1:ncells) = qvals(ncells+1:2*ncells)
        ndim = 2*ncells

        ! if source included
        if(size(qvals)>ndim)then
            rti%locations(1,:) = qvals(ndim+1:ndim+4*nsrc:4)
            rti%locations(2,:) = qvals(ndim+2:ndim+4*nsrc:4)
            rti%locations(3,:) = qvals(ndim+3:ndim+4*nsrc:4)
            rti%locations(4,:) = qvals(ndim+4:ndim+4*nsrc:4)
        endif
        
        do i = 1, size(model%vp,3)
            do j = 1, size(model%vp,2)
                do k = 1, size(model%vp,1)
                    model%vp(k,j,i) = RTI%parameters(1,RTI%sites_id(k,j,i))
                    model%vs(k,j,i) = RTI%parameters(2,RTI%sites_id(k,j,i))
                    model%rho(k,j,i) = RTI%parameters(3,RTI%sites_id(k,j,i))
                enddo
            enddo
        enddo

    endsubroutine update_model

end module
