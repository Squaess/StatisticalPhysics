      PROGRAM main
        IMPLICIT NONE
        real, parameter :: PI = 4 * atan(1.0)
        integer, parameter :: L = 3
        integer, dimension(L) :: NI, PI
        integer, dimension(L,L) :: fi
        real P2(-180:180)
        integer :: i

        do i=1,L
            NI(I) = I+1
            PI(I) = I-1
        enddo
        NI(L) = 1
        PI(1) = L

        do i = -180,180
            P2(i) = 0.5 * (3 * cos(i*PI/180.0)**2 - 1.0)
        enddo

        call init_fi()
        call mcs()

        print *, PI
        print *, fi
        print *, "DUpa"

        contains
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
                        call trial(i,j)
                    enddo
                enddo
            end subroutine mcs

            subroutine trial(i, j)
                integer, intent(in) :: i,j
                print *, "TRIAL", i, j
            end subroutine trial
      END PROGRAM