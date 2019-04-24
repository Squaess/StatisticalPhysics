      PROGRAM ising
        IMPLICIT NONE
        integer, parameter :: L = 10
        integer, dimension(L,L) :: S
        integer, dimension(L) :: NI, PI
        integer :: k = 230 000
        real :: T = 1.7

        integer :: i, j

        DO I=1,L
            NI(I)=I+1
            PI(I)=I-1
        ENDDO 
        NI(L)=1 
        PI(1)=L

        do i=1,L
            do j=1,L
                S(i,j) = 1
            enddo
        enddo

        open(unit=2, file="flips.csv")
        write(2, *) "k,m"
        do i=1,k
            call mcs(t)
            if ((k .ge. 30 000) .and. (mod(k,100) .eq. 0)) then 
                write(2, *) i, ",", calc_m()
            end if
        enddo

        contains
            function calc_m() result(m)
                real :: m
                integer :: i,j
                m = 0
                do i=1,L
                    do j=1,L
                        m = m + S(i,j)
                    enddo
                enddo
                m = m/(L*L)
            end function calc_m

            subroutine trial(i, j, t)
                integer, intent(in) :: i, j
                real, intent(in) :: t
                integer :: su
                real :: du, paccept, ptrial

                su = S(NI(I),J) + S(PI(I), J) + S(I,NI(J)) + S(I,PI(J))
                du = 2 * S(i,j) * su
                if (du .lt. 0) then
                    S(i,j) = -S(i,j)
                else
                    call random_number(paccept)
                    ptrial = EXP(-(du)/t)
                    if (paccept .le. ptrial) then
                        s(i,j) = -s(i,j)
                    end if
                endif

            end subroutine trial

            subroutine mcs(t)
                real, intent(in) :: t
                integer i, j

                do i=1,L
                    do j=1,L
                        call trial(i,j,t)
                    enddo
                enddo
            end subroutine mcs

            subroutine print_S()
                integer :: i, j
                open(unit=2, file="lattice5-100.csv")
                write(2, *) "x,y,m"
                do i=1,L
                    do j=1,L
                        write(2, *) i,",", j, "," ,S(i,j)
                    enddo
                enddo
                close(2)
! 100             format (I2," ")
            end subroutine print_S

      END PROGRAM


