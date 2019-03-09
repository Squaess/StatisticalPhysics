      integer d,N,step,x,K,i,j
      real xNa, xNa2, u
c ilosc krokow dla jednego pijaka
      integer, dimension(10) :: nu
      nu(1) = 10000
      nu(1) = 10
      nu(2) = 21.5443469
      nu(3) = 46.41588834
      nu(4) = 100.
      nu(5) = 215.443469
      nu(6) = 464.15888336
      nu(7) = 1000.
      nu(8) = 2154.43469003
      nu(9) = 4641.58883361
      nu(10) = 10000.
c ilosc pijakow
      K = 30000
      open(1, file='data.txt')
      d=-1
      do j=1,size(nu)
          N = nu(j)
          xNa = 0
          xNa2 = 0

          do i=1,K
              x=0
              do step=1,N
                if (ran1(d)<0.5) then
                    x = x + 1
                else
                    x = x - 1
                endif
              enddo
              xNa = xNa + x
              xN2a = xN2a + x*x
          enddo
c          write(*,*) "N", log(float(N)) 
c          write(*,*) xNa, xN2a
          u = sqrt((xN2a/float(K)) - ((xNa/float(K)) * (xNa/float(K))))
          write(1,*) N, "," ,u
c          write(*,*) log(u)
      enddo
      close(1)
      end