      PROGRAM ising
        IMPLICIT NONE
        integer, parameter :: L = 50
        integer, dimension(L,L) :: S
        integer, dimension(L) :: NI, PI
        integer :: k = 530000
        real :: T = 0
        real, dimension(-8:8) :: Tex
        real, dimension(50) :: Ta
        real :: ma, ea, avgm, avgm2
        real :: sus, avge, avge2, heatc
        real :: bind_m4, bind_m2
        real :: bind

        integer :: i, j, counter

        DO I=1,L
            NI(I)=I+1
            PI(I)=I-1
        ENDDO 
        NI(L)=1 
        PI(1)=L

        call init_Ta_bind(Ta)

        call initS(S, L)

        ! open(unit=2, file="sus100.csv")
        ! write(2, *) "T,m,L"
        ! open(unit=3, file="heat100.csv")
        ! write(3, *) "T,e,L"
        open(unit=4, file="binder50.csv")
        write(4, *) "T,b,L"

        do j=1,50
        print *, 2*j, "%"
        T = ta(j)
        ! initialize exp(-de/T) array
        call initTex(tex,T)
        ma = 0
        ea = 0
        avge = 0
        avge2 = 0
        avgm = 0
        avgm2 = 0
        counter = 0
        bind_m4 = 0
        bind_m2 = 0
        do i=1,k
            call mcs(t)
            if ((i >= 30000) .and. (mod(i,100) == 0)) then 
                ! write(2, *) i, ",", calc_m()
                avgm = avgm + abs(calc_m())
                avgm2 = avgm2 + (calc_m()**2)
                ma = ma + abs(calc_m())
                ea = ea + calc_ener()
                avge = avge + calc_ener2()
                avge2 = avge2 + (calc_ener2()**2)
                bind_m4 = bind_m4 + (calc_m2()**4)
                bind_m2 = bind_m2 + (calc_m2()**2)
                counter = counter + 1
            end if
        enddo
        ma = ma/float(counter)
        ! write(2, *) T, ",", ma, ",", L
        ea = ea/float(counter)
        avge = avge/float(counter)
        avge2 = avge2/float(counter)
        avgm = avgm/float(counter)
        avgm2 = avgm2/float(counter)
        bind_m4 = bind_m4/float(counter)
        bind_m2 = bind_m2/float(counter)
        sus = ((L*L)/T)*(avgm2 - (avgm**2))
        heatc = (1/(float(L**2)*(T**2)))*(avge2 - (avge**2))
        bind = 1 - (bind_m4/(3* (bind_m2**2)))
        ! write(3, *) T, ",", ea, ",", L
        ! write(2, *) T, ",", sus, ",", L
        ! write(3, *) T, ",", heatc, ",", L
        write(4, *) T, ",", bind, ",", L
        enddo
        ! close(2)
        close(3)
        close(4)

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

            function calc_m2() result(m)
                real :: m
                integer :: i,j
                m = 0
                do i=1,L
                    do j=1,L
                        m = m + S(i,j)
                    enddo
                enddo
            end function calc_m2

            function calc_ener() result(en)
                real :: en
                integer :: i, j, su
                en = 0
                do i=1,L
                    do j=1,L
                        su = S(NI(i), j) + &
                            S(PI(i), j) + &
                            S(i, NI(j))+ &
                            S(i, PI(j))
                        en = en + ((-su) * S(i,j))
                    enddo
                enddo
                en = en/((L*L)*2.0)
            end function calc_ener

            function calc_ener2() result(en)
                real :: en
                integer :: i, j, su
                en = 0
                do i=1,L
                    do j=1,L
                        su = S(NI(i), j) + &
                            S(PI(i), j) + &
                            S(i, NI(j))+ &
                            S(i, PI(j))
                        en = en + ((-su) * S(i,j))
                    enddo
                enddo
                en = en/(2.0)
            end function calc_ener2

            function calc_du(x, y) result(du)
                integer :: du
                integer, intent(in) :: x,y
                integer :: su ! suma

                su = S(NI(x), y) + &
                    S(PI(x), y) + &
                    S(x, NI(y))+ &
                    S(x, PI(y))
                du = 2 * S(x,y) * su
            end function calc_du

            subroutine trial(i, j, t)
                integer, intent(in) :: i, j
                real, intent(in) :: t
                integer :: du
                real :: paccept, ptrial

                du = calc_du(i, j)
                if (du < 0) then
                    S(i,j) = -S(i,j)
                else
                    call random_number(paccept)
                    ptrial = tex(du)
                    if (paccept <= ptrial) then
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

            subroutine print_S_out()
                integer :: i, j
                do i=1,L
                    print *, S(i,:)
                enddo
            end subroutine print_S_out

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
            end subroutine print_S

      END PROGRAM

      subroutine initTex(tex,tmp)
          real, dimension(-8:8), intent(inout) :: tex
          real, intent(inout) :: tmp
          integer :: i
          do i=-8,8,4
            tex(i) = exp((-i)/tmp)
          enddo
      end subroutine initTex

      subroutine initS(s, l)
          integer, intent(in) :: L
          integer, dimension(l,l) :: S
          integer :: i,j
          do i=1,L
            do j=1,L
                S(i,j) = 1
            enddo
          enddo
      end subroutine initS

      subroutine init_Ta(ta)
          real, intent(inout), dimension(50) :: ta
          integer :: i, j = 0
          do i=1,10
            ta(i) = 1.5 + i * ((2.0-1.5)/10)
          enddo

          j=1
          do i=11,40
            ta(i) = 2.0 + j * ((2.7-2.0)/30)
            j = j+1
          enddo

          j = 0
          do i=41,50
            ta(i) = 2.7 + j * ((3.5-2.7)/10)
            j = j+1
          enddo
      end subroutine init_Ta

      subroutine init_Ta_100(ta)
          real, intent(inout), dimension(50) :: ta
          integer :: i, j = 0
          do i=1,10
            ta(i) = 1.5 + i * ((2.2-1.5)/10)
          enddo

          j=1
          do i=11,40
            ta(i) = 2.2 + j * ((2.5-2.2)/30)
            j = j+1
          enddo

          j = 1
          do i=41,50
            ta(i) = 2.5 + j * ((3.5-2.5)/10)
            j = j+1
          enddo
      end subroutine init_Ta_100

      subroutine init_Ta_reg(ta)
        real, intent(inout), dimension(50) :: ta
        integer :: i
        do i=1,50
            ta(i) = 1.5 + i * ((3.5-1.5)/50)
        enddo
      end subroutine init_Ta_reg

      subroutine init_Ta_bind(ta)
        real, intent(inout), dimension(50) :: ta
        integer :: i
        do i=1,50
            ta(i) = 1.8 + i * ((2.5-1.8)/50)
        enddo
      end subroutine init_Ta_bind
