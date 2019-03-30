      implicit none
      REAL :: random
      integer, parameter :: L=5            ! size of lattice
      integer, parameter :: S=1             ! number of realization
      integer :: K=2                        ! steps
      integer, parameter ::N=10            ! particles
      integer, dimension(L) :: IN, IP       ! next and previous index 
      logical, dimension(L,L) :: A          ! if occupies?
      integer :: i,j,ii                   ! loops iterator
      integer, dimension(N) :: x, y         ! point location
      integer :: xt, yt, dxt, dyt
      integer, dimension(N) :: dx, dy

      ! initialization of index arrays
      do i=1,L
        IN(i) = i+1
        IP(i) = i-1
      enddo
      IN(L) = 1
      IP(1) = L
     
      print *, IN
      print *,"**********************"
      print *, IP
      print *,(N/float(L*L))


      ! this should be for every realization 
      ! initialize displacement array
      do i=1,N
        dx(i) = 0
        dy(i) = 0
      enddo

      ! initialize A array
      do i=1,L
        do j=1,L
            A(i,j) = .TRUE.
        enddo
      enddo

      ! initialize x,y for all particles
      i = N
      do while ( i .ge. 1)
        call random_number(random)
        xt = int(L*random) + 1
        call random_number(random)
        yt = int(L*random) + 1
        if (A(xt, yt)) then
            A(xt, yt) = .FALSE.
            x(i) = xt
            y(i) = yt
            i = i-1
        endif
      enddo

      call print2d(A, L)

      do i=1,K                              ! for every step
          do j=1,N                          ! for every particle
              xt = x(j)
              yt = y(j)
              dxt = 0
              dyt = 0
              call random_number(random)
              SELECT CASE (int(4*random))
                CASE (0)
                    xt = IN(x(j))
                    dxt = 1 
                CASE (1)
                    xt = IP(x(j))
                    dxt = -1
                CASE (2)
                    yt = IN(y(j))
                    dyt = 1
                CASE (3)
                    yt = IP(y(j))
                    dyt = -1
              ENDSELECT
              if (A(xt,yt)) then
                  A(xt, yt) = .FALSE.
                  A(x(j), y(j)) = .TRUE.
                  x(j) = xt
                  y(j) = yt
                  dx(j) = dx(j) + dxt
                  dy(j) = dy(j) + dyt
              endif
          enddo
      enddo
      end

      subroutine positions(x, y, n)
      implicit none
      integer N
      integer, dimension(N) :: x, y
      integer :: i
      do i=1,N
        print *,i,',',x(i),',',y(i)
      enddo
      end subroutine positions

      subroutine print2d(x, n)
      implicit none
      integer N
      integer, dimension(N,N) :: x
      integer :: i
      do i=1,N
        print *,x(i, :)
      enddo
      end subroutine print2d