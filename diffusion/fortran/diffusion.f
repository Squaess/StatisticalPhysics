      implicit none
      REAL :: random
      integer, parameter :: L=20            ! size of lattice
      integer, parameter :: S=10            ! number of realization
      integer :: K, MCS=50                ! steps
      integer, parameter ::N=400            ! particles
      integer, dimension(L) :: IN, IP       ! next and previous index 
      logical, dimension(L,L) :: A          ! if occupies?
      integer :: i,j,ii                   ! loops iterator
      integer, dimension(N) :: x, y         ! point location
      integer :: xt, yt, dxt, dyt
      integer, dimension(N) :: dx, dy
      real :: dr=0, drS = 0

      ! initialization of index arrays
      do i=1,L
        IN(i) = i+1
        IP(i) = i-1
      enddo
      IN(L) = 1
      IP(1) = L
     
      ! print *, IN
      ! print *,"**********************"
      ! print *, IP
      print *,(N/float(L*L))

      open(unit=2, file="max1.csv")
      write(2,*) "C ,", (N/float(L*L))

      do k=1,MCS
      ! this should be for every realization 
      drS = 0
      do ii=1,S
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

      ! call print2d(A, L)
      ! print *, "Prticles positions beggining"
      ! call positions(x, y, n)

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

      dr = 0
      do i=1,n
        dr = dr + dx(i)**2 + dy(i)**2
      enddo
      ! average for 1 simulation for N particles
      dr = dr/float(n)
      drS = drS + dr
      ! print *, dr

      ! end loop for many simluations
      enddo
      !average drS for simulations
      drS = drS/float(S)
      ! print *, K,drS/(4 * K)
      ! write(2,*) K,",",drS/(4 * K)
      

      enddo
      print *, K,drS/(4 * K)
      close(2)
      ! print *, "Perticle positions after K steps"
      ! call positions(x, y, n)
      ! ! displacement after K steps for N particles
      ! print *, "Displaycment"
      ! call positions(dx, dy, N)
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