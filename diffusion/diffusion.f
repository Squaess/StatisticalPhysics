      implicit none
      REAL :: random
      integer, parameter :: L=20            ! size of lattice
      integer, parameter :: S=10            ! number of realization
      integer :: K=20000                    ! steps
      integer, parameter ::N=200            ! particles
      integer, dimension(L) :: IN, IP       ! next and previous index 
      logical, dimension(L,L) :: A          ! if occupies?
      integer :: i,j,ii,t                   ! loops iterator
      integer, dimension(N) :: x, y         ! point location
      integer :: xt, yt, dxt, dyt
      integer, dimension(N) :: dx, dy
      integer :: r2
      real :: r2s, r2ii
      integer :: jj

      print *,(N/float(L*L))

      r2ii = 0
      do t=1,50
        do jj=1,t
          r2s = 0
          do ii=1,S 
        
              ! initialize displacement array
              do i=1,N
                dx(i) = 0
                dy(i) = 0
              enddo

              ! initialization of index arrays
              do i=1,L
                IN(i) = i+1
                IP(i) = i-1
              enddo
              IN(L) = 1
              IP(1) = L

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
                            dyt = -1
                            yt = IP(y(j))
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

              r2 = 0
              do i=1,N
                r2 = r2 + (dx(i)**2 + dy(i)**2)
              enddo
              r2s = r2s + (float(r2)/float(N))
          enddo
          r2s = r2s/float(S)
          r2ii = r2ii + r2s
          ! print *,r2s/float(4*t)

        enddo
        print *, r2ii/float(4*t)
      !print *,r2ii/(4*float(t))

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
