      PROGRAM main
        IMPLICIT NONE
        real, parameter :: PI_math = 4 * atan(1.0)
        integer, parameter :: L = 3
        integer, dimension(L) :: NI, PI
        integer, dimension(L,L) :: fi
        real, dimension(2,2) :: Q
        integer :: k = 1230000
        real P2(-180:180)
        integer :: i, counter
        real :: T
        REAL :: A( 2, 2 ), W( 2 ), WORK( 1000 )
        INTEGER :: INFO, LWORK, LIWORK
        INTEGER :: IWORK( 1000 )

        Q(1,1) = 0
        Q(1,2) = 0
        Q(2,1) = 0
        Q(2,2) = 0

        do i=1,L
            NI(I) = I+1
            PI(I) = I-1
        enddo
        NI(L) = 1
        PI(1) = L

        do i = -180,180
            P2(i) = 0.5 * (3 * cos(i*PI_math/180.0)**2 - 1.0)
        enddo

        call init_fi()
        call mcs()

        T = 0.1

        counter = 0
        do i=1,k
            call mcs()
            if ((i >= 30000) .and. (mod(i,100) == 0)) then
                call calc_q()
                counter = counter + 1
            endif
        enddo

        Q(2,1) = Q(1,2)
        Q(2,2) = - Q(1,1)

        print *, fi
        print *, "DUpa"

        print *,Q

        CALL SSYEVD( 'Vectors', 'Upper', 2, Q, 2, W, WORK, -1, &
               IWORK, -1, INFO )
        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
        print *, INFO


        contains

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

            subroutine mcs()
                integer :: i, j
                do i=1,L
                    do j=1,L
                        call trial(i,j, 0.5)
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
                else if (p .lt. exp(diff/T)) then
                    fi(i,j) = fi_new
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