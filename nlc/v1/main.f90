      PROGRAM main
        IMPLICIT NONE
        real, parameter :: PI_math = 4 * atan(1.0)
        integer, parameter :: L = 10
        integer, dimension(L) :: NI, PI
        integer, dimension(L,L) :: fi
        real, dimension(2,2) :: Q
        integer :: k = 230000
        real P2(-180:180)
        integer :: i, counter, j
        real :: T, eig
        real, dimension(100) :: TA

        call init_ta()

        do i=1,L
            NI(I) = I+1
            PI(I) = I-1
        enddo
        NI(L) = 1
        PI(1) = L

        do i = -180,180
            P2(i) = 0.5 * (3 * cos(i*PI_math/180.0)**2 - 1.0)
        enddo

        open(unit=2, file='result.csv')
        write(2, *) "T, e"
        ! start loop for every temp
        do j=1,100

        T = TA(j)

        Q(1,1) = 0
        Q(1,2) = 0
        Q(2,1) = 0
        Q(2,2) = 0

        call init_fi()

        counter = 0
        ! start loop for k MCS steps
        do i=1,k
            call mcs(T)
            if ((i >= 30000) .and. (mod(i,1000) == 0)) then
                call calc_q()
                counter = counter + 1
            endif
        enddo

        Q(2,1) = Q(1,2)
        Q(2,2) = - Q(1,1)

        ! print *, fi
        ! print *, Q(1,1), Q(1,2)
        ! print *, Q(2,1), Q(2,2)
        eig = calc_eigen()/float(counter)
        write(2,*) T, "," , eig
        print *, j, "%"
        enddo
        ! end of loop for temp
        close(2)

        contains

            subroutine init_ta
                integer :: i
                do i=1,100
                    TA(I) = 0.1 + I*((3.7-0.1)/100.)
                enddo
            end subroutine init_ta

            function calc_eigen() result(eigen)
                !     .. Parameters ..
                INTEGER          N
                PARAMETER        ( N = 2 )
                INTEGER          LDA, LDVL, LDVR
                PARAMETER        ( LDA = N, LDVL = N, LDVR = N )
                INTEGER          LWMAX
                PARAMETER        ( LWMAX = 1000 )
                ! 
                INTEGER          INFO, LWORK
                !
                REAL             VL( LDVL, N ), VR( LDVR, N )
                REAL             WR( N ), WI( N ), WORK( LWMAX )
                REAL eigen
                LWORK = -1
                CALL SGEEV( 'Vectors', 'Vectors', N, Q, LDA, WR, WI, VL, LDVL, &
                    VR, LDVR, WORK, LWORK, INFO )
                LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
                CALL SGEEV( 'Vectors', 'Vectors', N, Q, LDA, WR, WI, VL, LDVL, &
                    VR, LDVR, WORK, LWORK, INFO )

                ! print *, INFO
                ! print *, WR
                ! print *, WI
                ! print *, MAX(WR(1), WR(2))
                eigen = MAX(WR(1), WR(2))
            end function calc_eigen

            subroutine calc_q()
                integer :: i,j
                DO i=1,L
                    DO j=1,l
                        Q(1,1) = Q(1,1) + &
            ((2.0 * cos(fi(i,j)*PI_math/180.0)**2 - 1.0)/ (L*L))
                        Q(1,2) = Q(1,2) + &
            ((2.0 * cos(fi(i,j)*PI_math/180.0) &
            * sin(fi(i,j)*PI_math/180.0)) &
            / (L*L))
                    enddo
                enddo
            end subroutine calc_q

            subroutine init_fi()
                integer :: i, j
                do i = 1,L
                    do j = 1,L
                        fi(i, j) = 45
                    enddo
                enddo
            end subroutine init_fi

            subroutine mcs(T)
                real :: T
                integer :: i, j
                real :: pi, pj
                do i=1,L
                    do j=1,L
                        call random_number(pi)
                        call random_number(pj)
                        call trial(int(pi * L) + 1, int(pj * L) + 1, T)
                    enddo
                enddo
            end subroutine mcs

            subroutine trial(i, j, T)
                integer, intent(in) :: i,j
                real, intent(in) :: T
                integer :: fi_new, diff
                real :: p, p2
                call random_number(p)
                p = p - 0.5
                fi_new = fi(i,j) + int(p * 180)
                if (fi_new .gt. 90) then
                    fi_new = fi_new - 180
                else if (fi_new .lt. (-90)) then
                    fi_new = fi_new + 180
                endif
                ! print *, "TRIAL", i, j
                ! print *, fi_new
                call random_number(p2)
                diff = calc_diff_ene(i,j,fi_new)
                if (diff .gt. 0) then
                    fi(i,j) = fi_new
                    ! print *, "changing due to diff energy"
                else if (p .lt. exp(diff/T)) then
                    fi(i,j) = fi_new
                    ! print *, "changing due to exp"
                else
                    ! print *, "Not changing"
                endif

            end subroutine trial

            function calc_diff_ene(i, j, new_fi) result(diff)
                integer, intent(in) :: i, j, new_fi
                integer :: old, new
                integer :: diff
                old = P2(fi(i,j) - fi(NI(i), j)) + &
                    P2(fi(i,j) - fi(PI(i), j)) +  &
                    P2(fi(i,j) - fi(i, NI(j))) +  &
                    P2(fi(i,j) - fi(i, PI(j)))
                old = old * -1

                new = P2(new_fi - fi(NI(i), j)) + &
                    P2(new_fi - fi(PI(i), j)) + &
                    P2(new_fi - fi(i, NI(j))) + &
                    P2(new_fi - fi(i, PI(j)))
                new = new * -1

                diff = old - new

            end function calc_diff_ene

      END PROGRAM