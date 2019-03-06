      integer d,N,step,x,K,i,j
      real xNa, xNa2, u
c ilosc krokow dla jednego pijaka
      integer, dimension(3) :: nu
      nu(1) = 10
      nu(2) = 100
      nu(3) = 10000
c ilosc pijakow
      K = 30 000

      xNa = 0
      xNa2 = 0

      d=-1
      do j=1,size(nu)
          N = nu(j)
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
          write(*,*) "N", log(float(N)) 
          write(*,*) xNa, XN2a
          u = sqrt((xN2a/float(K)) - ((xNa/float(K)) * (xNa/float(K))))

          write(*,*) log(u)
      enddo

      end

