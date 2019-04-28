      implicit none

      integer, parameter :: L = 10
      ! table of nearest neighbour
      integer, dimension(L) :: NI, PI
      ! table of spins
      integer, dimension(L,L) :: S
      integer :: MCS = 530 000
      ! iterator
      integer :: i, j, k
      real :: T = 1.7
      ! number of temperature points
      integer, parameter :: nt = 100
       real, dimension(nt) :: TA
      ! average magnetization
      real :: ma = 0
      ! real :: m = 0
      integer :: maCounter = 0

      ! initialize PBC
      do i=1,L
        NI(i) = i+1
        PI(i) = i-1
      enddo
      NI(L) = 1
      PI(1) = L

      do i=1,L
        do j=1,L
            S(i,j) = 1
        enddo
      enddo

      ! initialize temp values
      do i=1,nt
        TA(i) = 1.2 + i*(3.0/float(nt))
      enddo

      open(unit=2, file='mag10.csv')
      write(2,*) "T, m, L"
      ! write(2,*) "MCS, m"

      do i=1,nt
        print *, i/float(nt)*100, "%"
        T = TA(i)
        ! T = TA(1)
        ma = 0.
        maCounter = 0
        do k=1,MCS
            call mcstep(T)
            if (k .gt. 30 000 .AND. mod(k,100) .eq. 0) then
                call calc_ma(ma)
                maCounter = maCounter + 1
                ! print *, k, ",", calc_m()
                ! write(2,*) k, ",", calc_m()
            end if
        enddo
        ma = ma/float(maCounter)
        write(2,*) T, ",", ma,",",L
      enddo

      close(2)

      contains
        subroutine calc_ma(ma)
            real, intent(inout) :: ma
            real :: m = 0
            integer :: i, j
            m = 0
            do i=1,L
                do j=1,L
                    m = m + S(i,j)
                enddo
            enddo
            m = m/float(L*L)
            ma = ma + abs(m)
        end subroutine calc_ma

        function calc_m() result(m)
            real :: m
            integer :: i, j
            m = 0
            do i=1,L
                do j=1,L
                    m = m + S(i,j)
                enddo
            enddo
            m = m/(L*L)
        end function calc_m

        subroutine trial(i,j,t)
            integer, intent(in) :: i, j
            real, intent(in) :: t
            integer :: dU, suma
            real :: r
            suma = S(NI(i),j) + S(PI(i), j) + S(i, NI(j)) + S(i, PI(j))
            du = 2 * S(i,j) * suma
            if (du .lt. 0) then
                S(i,j) = S(i,j) * (-1)
            else
                call random_number(r)
                if (r .lt. exp(-du/t)) then
                    S(i,j) = S(i,j) * (-1)
                endif
            end if
        end subroutine trial

        subroutine mcstep(t)
            real, intent(in) :: t
            integer :: i, j
            i = 0
            j = 0
            do i=1,L
                do j=1,L
                    call trial(i,j,t)
                enddo
            enddo
        end subroutine mcstep
      end

      subroutine printS(S, L)
        integer, intent(in) :: L
        integer, dimension(L,L), intent(in) :: S
        integer :: i = 0
        do i=1,L
            print *, S(i,:)
        enddo
        print *, "**************************"
        read (*,*)
      end subroutine printS
