      program sum
      implicit none
      real xNa, xN2a
      integer c, N, R, d, J, i, K, x
      K = 10
      N = 3
      do J=1,K
        x=0
        do i=1,N
            if (ran1(1) < 0.5) then
                x = x + 1
            else
                x = x -1
            endif
        enddo
        write(*,*) x
      enddo
      end program

