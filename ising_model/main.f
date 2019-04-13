      implicit none

      integer, parameter :: L = 25
      ! table of nearest neighbour
      integer, dimension(L) :: IN, IP
      ! table of spins
      integer, dimension(L,L) :: S
      integer :: MCS = 530 000
      ! iterator
      integer :: i, j, k
      real :: T = 0
      ! number of temperature points
      integer, parameter :: nt = 300
      real, dimension(nt) :: TA
      ! average magnetization
      real :: ma = 0
      integer :: maCounter = 0

      ! initialize PBC
      do i=1,L
        IN(i) = i+1
        IP(i) = i-1
      enddo
      IN(L) = 1
      IP(1) = L

      do i=1,L
        do j=1,L
            S(i,j) = 1
        enddo
      enddo

      ! initialize temp values
      do i=1,nt
        TA(i) = 0. + i*(5./float(nt))
      enddo

      open(unit=2, file='magnetization25.csv')
      write(2,*) "T, m, L"

      do i=1,nt
        print *, i/float(nt)*100, "%"
        T = TA(i)
        ma = 0.
        maCounter = 0
        do k=1,MCS
            call mcstep(T)
            if (k .gt. 30 000 .AND. mod(k,100) .eq. 0) then
                call calc_ma(ma)
                maCounter = maCounter + 1
            end if
        enddo
        ma = ma/float(maCounter)
        write(2,*) T, ",", ma,",",L
      enddo

      close(2)

      contains
        subroutine calc_ma(ma)
            real, intent(inout) :: ma
            real :: m
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

        subroutine trial(i,j,t)
            integer, intent(in) :: i, j
            real, intent(in) :: t
            integer :: dU, suma
            real :: r
            suma = S(IN(i),j) + S(IP(i), j) + S(i, IN(j)) + S(i, IP(j))
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
