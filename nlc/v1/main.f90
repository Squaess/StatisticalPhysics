      PROGRAM main
        IMPLICIT NONE
        real, parameter :: PI_math = 4 * atan(1.0)
        integer, parameter :: L = 20
        integer, dimension(L) :: NI, PI
        integer, dimension(L,L) :: fi
        real, dimension(2,2) :: Q
        integer :: k = 230000
        real P2(-180:180)
        integer :: i, counter, j
        real :: T, eig, q11, q12
        real, dimension(100) :: TA

        call init_ta()

        ! periodic boundary condition
        do i=1,L
            NI(I) = I+1
            PI(I) = I-1
        enddo
        NI(L) = 1
        PI(1) = L

        do i = -180,180
            P2(i) = 0.5 * (3 * (cos(i*PI_math/180.0)**2) - 1.0)
        enddo

        open(unit=2, file='result.csv')
        write(2, *) "T, e"
        ! start loop for every temp
        do j=1,SIZE(TA)

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
            if ((i >= 30000) .and. (MOD(i,100) == 0)) then
                call calc_q(q11, q12)
                Q(1,1) = Q(1,1) + q11
                Q(1,2) = Q(1,2) + q12
                counter = counter + 1
            endif
        enddo

        Q(1,1) = Q(1,1)/float(counter)
        Q(1,2) = Q(1,2)/float(counter)
        Q(2,1) = Q(1,2)
        Q(2,2) = - Q(1,1)

        ! eig = calc_eigen()/float(counter)
        eig = calc_eigen()
        write(2,*) T, "," , eig
        print *, int(j/float(size(ta)) * 100), "%"

        enddo
        ! end of loop for temp
        close(2)

        contains

            subroutine init_ta
                integer :: i
                do i=1,100
                    TA(I) = 0.0 + I*((1.5-0.0)/100.)
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

            subroutine calc_q(q11, q12)
                real, intent(inout) :: q11, q12
                integer :: i,j
                q11 = 0
                q12 = 0
                DO i=1,L
                    DO j=1,l
                        ! moze na koncu prze L**2?
                        ! Q(1,1) = Q(1,1) + &
                        !     (2.0 * cos(fi(i,j)*PI_math/180.0)**2 - 1.0)
                        ! Q(1,2) = Q(1,2) + &
                        !     (2.0 * cos(fi(i,j)*PI_math/180.0) &
                        !     * sin(fi(i,j)*PI_math/180.0))
                        q11 = q11 + &
                            (2.0 * cos(fi(i,j)*PI_math/180.0)**2 - 1.0)
                        q12 = q12 + &
                            (2.0 * cos(fi(i,j)*PI_math/180.0) &
                            * sin(fi(i,j)*PI_math/180.0))
                    enddo
                enddo
                ! Q(1,1) = Q(1,1)/float(L*L)
                ! Q(1,2) = Q(1,2)/float(L*L)
                q11 = q11 / float(L*L)
                q12 = q12 / float(L*L)
            end subroutine calc_q

            subroutine init_fi()
                ! how about rand value for angle
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
                ! real :: pi, pj
                do i=1,L
                    do j=1,L
                        ! call random_number(pi)
                        ! call random_number(pj)
                        ! call trial(int(pi * L) + 1, int(pj * L) + 1, T)
                        call trial(i, j, T)
                    enddo
                enddo
            end subroutine mcs

            subroutine trial(i, j, T)
                integer, intent(in) :: i,j
                real, intent(in) :: T
                integer :: fi_new
                real :: p, diff
                fi_new = get_new_fi(fi(i,j))
                call random_number(p)
                diff = calc_diff_ene(i,j,fi_new)
                if (diff .lt. 0) then
                    fi(i,j) = fi_new
                else if (p .le. exp(-diff/T)) then
                    fi(i,j) = fi_new
                endif
            end subroutine trial

            function get_new_fi(old_fi) result(new_fi)
                integer, intent(in) :: old_fi
                integer :: new_fi
                real :: p
                call random_number(p)
                new_fi = old_fi + nint((p-0.5)*10)
                if (new_fi > 90) then
                    new_fi = new_fi - 180
                else if (new_fi < -90) then
                    new_fi = new_fi + 180
                endif
            end function get_new_fi

            function calc_diff_ene(i, j, new_fi) result(diff)
                integer, intent(in) :: i, j, new_fi
                real :: old, new
                real :: diff
                old = calc_energy(i,j,fi(i,j))
                new = calc_energy(i,j,new_fi)
                diff = new - old
            end function calc_diff_ene

            function calc_energy(i, j, angle) result(res)
                integer, intent(in) :: i,j,angle
                real :: res
                res = -P2(angle - fi(NI(i), j)) &
                      -P2(angle - fi(PI(i), j)) &
                      -P2(angle - fi(i, NI(j))) &
                      -P2(angle - fi(i, PI(j)))

            end function calc_energy


      END PROGRAM