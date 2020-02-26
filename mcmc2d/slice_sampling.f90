! Slice sampling
!
module sliceSample
    use m_logger,      only : log_msg
    use m_utils,       only : ii10, rtoa
    use m_settings,    only : T_MCMC_SET
    use like_settings, only : T_DATA, T_LIKE_SET, T_LIKE_BASE
    use m_likelihood,  only : T_LIKE, likelihood
    use run_info,      only : T_RUN_INFO
    use mt19937,       only : grnd, gasdev

    implicit none

    public :: slice_sample
    public :: slice_sample_pos

contains

    subroutine slice_sample(target_data,RTI,like_set,mcmc_set,idxp,idx,lgP,w,n) 
        implicit none

        type(T_DATA), dimension(:), intent(in)      :: target_data
        type(T_RUN_INFO), intent(inout)             :: RTI
	    type(T_LIKE_SET), intent(in)                :: like_set
	    type(T_MCMC_SET), intent(in)                :: mcmc_set
        ! index of row of parameters (i.e. vp, vs or rho)
        integer, intent(in)                         :: idxp
        integer, intent(in)                         :: idx
        !> The Loglikelihood bound
        type(T_LIKE), intent(inout)            :: lgP
        !> The initial width
        real(kind=ii10), intent(in) :: w
        !> The number of likelihood calls
        integer, intent(inout), optional :: n

        ! The upper bound
        real(kind=ii10)    :: R, lgR
        ! The lower bound
        real(kind=ii10)    :: L, lgL
        real(kind=ii10)    :: x0, x
        real(kind=ii10)    :: lgP0, slice_lgp, lgX


        real(kind=ii10) :: temp_random

        integer :: i_step

        !call log_msg('Begin likelihood: '//rtoa(lgP%like))
        !call log_msg('Begin vp: '//rtoa(RTI%parameters(idxp,idx)))
        ! define a slice
        lgP0 = lgP%like
        slice_lgp = lgP0 -  log(grnd())

        !call log_msg('slice like: '//rtoa(slice_lgp))
        ! Select initial start and end points
        x0 = RTI%parameters(idxp,idx)
        temp_random = grnd()
        L = x0 -   temp_random   * w 
        R = x0 + (1-temp_random) * w 

        ! temporaly use
        RTI%randcount =  RTI%randcount + 2
        ! Calculate initial likelihoods
        RTI%parameters(idxp,idx) = L
        call likelihood(target_data,RTI,like_set, lgP)
        lgL = lgP%like
        RTI%parameters(idxp,idx) = R
        call likelihood(target_data,RTI,like_set, lgP)
        lgR = lgP%like

        ! expand R until it's outside the likelihood region
        if(present(n)) n = n + 2
        i_step=0
        do while( lgR <= slice_lgp .and. lgR /= huge(lgR) )
            i_step=i_step+1
            R = x0 + w * i_step
            if(R>mcmc_set%vpmax)then
                lgR = huge(lgR)
                R = mcmc_set%vpmax
                exit
            endif
            RTI%parameters(idxp,idx) = R
            call likelihood(target_data,RTI,like_set, lgP)
            lgR = lgP%like
            if(present(n))    n = n + 1
            !call log_msg('vp: '//rtoa(R))
            !call log_msg('likelihood: '//rtoa(lgR))
        end do
        if(i_step>100) write(*,'(" too many R steps (",I10,")")') i_step

        ! expand L until it's outside the likelihood region
        i_step=0
        do while(lgL <= slice_lgP .and. lgL /= huge(lgL) )
            i_step=i_step+1
            L = x0 - w * i_step
            if(L<mcmc_set%vpmin)then
                lgL = huge(lgL)
                L = mcmc_set%vpmin
                exit
            endif
            RTI%parameters(idxp,idx) = L
            call likelihood(target_data,RTI,like_set, lgP)
            lgL = lgP%like
            if(present(n)) n = n + 1
            !call log_msg('vp: '//rtoa(L))
            !call log_msg('likelihood: '//rtoa(lgL))
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
                write(*,'("Polychord Warning: Non deterministic loglikelihood")')
                lgX = huge(lgX)
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
            ! calculate the likelihood 
            RTI%parameters(idxp,idx) = x1
            call likelihood(target_data,RTI,like_set, lgP)
            lgX = lgP%like
            if(present(n)) n =  n + 1

            ! If we're not within the likelihood bound then we need to sample further
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

    subroutine slice_sample_pos(target_data,RTI,like_set,mcmc_set,idxp,idx,lgP,w,n) 
        implicit none

        type(T_DATA), dimension(:), intent(in)      :: target_data
        type(T_RUN_INFO), intent(inout)             :: RTI
	    type(T_LIKE_SET), intent(in)                :: like_set
	    type(T_MCMC_SET), intent(in)                :: mcmc_set
        ! index of row of parameters (i.e. vp, vs or rho)
        integer, intent(in)                         :: idxp
        integer, intent(in)                         :: idx
        !> The Loglikelihood bound
        type(T_LIKE), intent(inout)            :: lgP
        !> The initial width
        real(kind=ii10), intent(in) :: w
        !> The number of likelihood calls
        integer, intent(inout), optional :: n
        
        ! boundary
        real(kind=ii10), dimension(2) :: bmin, bmax

        ! The upper bound
        real(kind=ii10)    :: R, lgR
        ! The lower bound
        real(kind=ii10)    :: L, lgL
        real(kind=ii10)    :: x0, x
        real(kind=ii10)    :: lgP0, slice_lgp, lgX


        real(kind=ii10) :: temp_random

        integer :: i_step

        ! set up bounday
        bmin(1) =  mcmc_set%grid%xmin
        bmin(2) =  mcmc_set%grid%ymin
        bmax(1) =  mcmc_set%grid%xmax
        bmax(2) =  mcmc_set%grid%ymax
        !call log_msg('Begin likelihood: '//rtoa(lgP%like))
        !call log_msg('Begin x: '//rtoa(RTI%points(idxp,idx)))
        ! define a slice
        lgP0 = lgP%like
        slice_lgp = lgP0 -  log(grnd())

        !call log_msg('slice like: '//rtoa(slice_lgp))
        ! Select initial start and end points
        x0 = RTI%points(idxp,idx)
        temp_random = grnd()
        L = x0 -   temp_random   * w 
        R = x0 + (1-temp_random) * w 

        ! temporaly use
        RTI%randcount =  RTI%randcount + 2
        ! Calculate initial likelihoods
        RTI%points(idxp,idx) = L
        call update_model(RTI, mcmc_set%grid)
        call likelihood(target_data,RTI,like_set, lgP)
        lgL = lgP%like
        RTI%points(idxp,idx) = R
        call update_model(RTI, mcmc_set%grid)
        call likelihood(target_data,RTI,like_set, lgP)
        lgR = lgP%like
        !call log_msg('L '//rtoa(L))
        !call log_msg('L likelihood: '//rtoa(lgL))
        !call log_msg('R '//rtoa(R))
        !call log_msg('R likelihood: '//rtoa(lgR))

        ! expand R until it's outside the likelihood region
        if(present(n)) n = n + 2
        i_step=0
        do while( lgR <= slice_lgp .and. lgR /= huge(lgR) )
            i_step=i_step+1
            R = x0 + w * i_step
            if(R>bmax(idxp))then
                lgR = huge(lgR)
                R = bmax(idxp)
                exit
            endif
            RTI%points(idxp,idx) = R
            call update_model(RTI, mcmc_set%grid)
            call likelihood(target_data,RTI,like_set, lgP)
            lgR = lgP%like
            if(present(n)) n = n + 1
            !call log_msg('x: '//rtoa(R))
            !call log_msg('likelihood: '//rtoa(lgR))
        end do
        if(i_step>100) write(*,'(" too many R steps (",I10,")")') i_step

        ! expand L until it's outside the likelihood region
        i_step=0
        do while(lgL <= slice_lgP .and. lgL /= huge(lgL) )
            i_step=i_step+1
            L = x0 - w * i_step
            if(L<bmin(idxp))then
                lgL = huge(lgL)
                L = bmin(idxp)
                exit
            endif
            RTI%points(idxp,idx) = L
            call update_model(RTI, mcmc_set%grid)
            call likelihood(target_data,RTI,like_set, lgP)
            lgL = lgP%like
            if(present(n)) n = n + 1
            !call log_msg('x: '//rtoa(L))
            !call log_msg('likelihood: '//rtoa(lgL))
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
                write(*,'("Polychord Warning: Non deterministic loglikelihood")')
                lgX = huge(lgX)
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
            ! calculate the likelihood 
            RTI%points(idxp,idx) = x1
            call update_model(RTI, mcmc_set%grid)
            call likelihood(target_data,RTI,like_set, lgP)
            lgX = lgP%like
            if(present(n)) n =  n + 1

            ! If we're not within the likelihood bound then we need to sample further
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


    end subroutine slice_sample_pos

    subroutine update_model(RTI, grid)
        use kdtree2_precision_module, only : kdkind
        use kdtree2_module, only : kdtree2, kdtree2_result, kdtree2_create,&
            kdtree2_destroy, kdtree2_n_nearest
        use m_settings, only : T_GRID
        use omp_lib
        implicit none
        type(T_RUN_INFO), intent(inout) :: RTI
        type(T_GRID), intent(in) :: grid

        ! kd-tree
        type(kdtree2), pointer :: tree
        type(kdtree2_result), dimension(1) :: results
        real(kdkind), dimension(2) :: qv

        ! grid
        integer        :: ix0, ix1, iy0, iy1, iz0, iz1

        integer i,j,k

        !real(kind=8) :: t

        ! create kd-tree
        tree => kdtree2_create(RTI%points(:,1:RTI%ncells), sort=.false., rearrange=.true.)

        ! convert to grid based model using kd-tree nearest search

        ! the lowest and highest indices for grid model vp, vs and rho
        ix0 = 1 
        ix1 = grid%nx
        iy0 = 1
        iy1 = grid%ny

        ! call nearest search of kd-tree
        !t =  omp_get_wtime()
        !$omp parallel
        !$omp do private(qv,i,j,results)
        do i = ix0, ix1
            do j = iy0, iy1
                    qv = [grid%xmin+(i-1)*grid%dx, grid%ymin+(j-1)*grid%dy]
                    call kdtree2_n_nearest(tp=tree,qv=qv,nn=1,results=results)
                    RTI%sites_id(j,i) = results(1)%idx
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        !t =  omp_get_wtime() -  t
        !call log_msg('parallelized kdtree code: '//rtoa(t) )

        call kdtree2_destroy(tree)

        return

    end subroutine

end module
