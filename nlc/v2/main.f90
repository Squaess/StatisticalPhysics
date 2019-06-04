program name
implicit none
real, parameter :: PI_MATH = 4 * atan(1.0)
! lettice size
integer, parameter :: L = 100
integer, dimension(L,L) :: fi
! number of monte carlo steps
integer :: k = 50000
! PBC
integer, dimension(L) :: NI, PI
! array of temperatures
real, dimension(50) :: TA
real, dimension(50) :: eigens
real :: T, eig
real P2(-180:180)
real, dimension(2,2) :: Q

integer :: i, j, c

call init_PBC()
call init_TA()
call init_P2()
call init_eigens()

! start loop for different temperatures
open(unit=2, file="result3.csv")
write(2, *) "T,S, S2"

do i=1,size(TA)
    Q(1,1) = 0
    Q(1,2) = 0
    Q(2,1) = 0
    Q(2,2) = 0

    call init_fi()
    c = 0
    T = TA(i)
    do j=1,k
        call mcs(T)
        if ((j > 5000) .and. (MOD(j,100) .eq. 0)) then
            call update_q()
            Q(2,1) = Q(1,2)
            Q(2,2) = -Q(1,1)
            call update_eigens(i)
            c = c + 1
        endif
    enddo
    Q(2,1) = Q(1,2)
    Q(2,2) = -Q(1,1)
    print *, "Completed: ", int(i/float(size(TA))*100), "%"
    eig = calc_eigen() / float(c)
    write(2, *) T, ",", eig, ",", eigens(i)/float(c)

enddo
! end of the loop for different temperatures
close(2)
contains

subroutine init_P2()
    integer :: i
    do i=-180,180
        P2(i) = 3.0/2.0 * cos(I * PI_MATH / 180.)**2 - 0.5
    enddo
end subroutine init_p2

subroutine init_PBC()
    integer :: i
    do i=1,L
        NI(i) = i+1
        pI(i) = i-1
    enddo
    NI(L) = 1
    PI(1) = L
end subroutine init_PBC

subroutine init_TA()
    integer :: i
    ! do i=1,40
    !     TA(i) = 0.0 + i*(0.7)/float(40)
    ! enddo
    ! do i=1,20
    !     TA(i+40) = 0.7 + i*(0.9-0.7)/float(20)
    ! enddo
    ! do i=1,40
    !     TA(i+60) = 0.9 + i*(1.5-0.9)/float(40)
    ! enddo
    do i=1,15
        TA(i) = 0.0 + i*(0.75)/float(15)
    enddo
    do i=1,20
        TA(i+15) = 0.75 + i*(1.-0.75)/float(20)
    enddo
    do i=1,15
        TA(i+30) = 1 + i*(1.5-1.0)/float(15)
    enddo
end subroutine init_TA

subroutine init_fi()
    integer :: i, j
    do i=1,L
        do j=1,L
            fi(i,j) = 0
        enddo
    enddo
end subroutine init_fi

subroutine mcs(T)
    real, intent(in) :: T
    integer :: i,j
    real :: pi, pj
    do i=1,L
        do j=1,L
            call random_number(pi)
            call random_number(pj)
            ! call trial(i,j,T)
            call trial(int(pi*L)+1,int(pj*L)+1,T)
        enddo
    enddo
end subroutine mcs

function calc_energy(i,j, angle) result(energy)
    real :: energy
    integer, intent(in) :: i, j, angle
    energy = -P2(angle - fi(NI(i), j)) &
             -P2(angle - fi(PI(i), j)) &
             -P2(angle - fi(i, NI(j))) &
             -P2(angle - fi(i, PI(j)))
    
end function calc_energy

function get_new_fi(old_fi) result(new_fi)
    integer, intent(in) :: old_fi
    integer :: new_fi
    real p
    call random_number(p)
    new_fi = old_fi + nint((p - 0.5) * 10)
    if (new_fi > 90) then
        new_fi = new_fi - 180
    else if (new_fi < -90) then
        new_fi = new_fi + 180
    endif
end function get_new_fi
    
subroutine trial(i,j,T)
    integer, intent(in) :: i,j
    real, intent(in) :: T
    real :: Uo, Ut, Ud, p, w
    integer :: new_fi
    new_fi = get_new_fi(fi(i,j))
    Uo = calc_energy(i,j, fi(i, j))
    Ut = calc_energy(i,j, new_fi)
    Ud = Ut - Uo
    if (Ud < 0) then
        fi(i,j) = new_fi
    else
        call random_number(p)
        w = exp(-Ud/T)
        if (p <= w) then
            fi(i,j) = new_fi
        endif
    endif
end subroutine trial

subroutine update_q()
    integer :: i, j
    real :: q11, q12
    q11 = 0
    q12 = 0
    do i=1,L
        do j=1,L
            q11 = q11 + 2.*(cos(fi(i,j) * PI_MATH / 180.)**2) - 1.
            q12 = q12 + 2.*cos(fi(i,j) * PI_MATH / 180.) * sin(fi(i,j)* PI_MATH / 180.)
        enddo
    enddo
    q11 = q11/float(L*L)
    q12 = q12/float(L*L)
    Q(1,1) = Q(1,1) + q11
    Q(1,2) = Q(1,2) + q12
end subroutine update_q

subroutine update_eigens(idx)
    integer, intent(in) :: idx
    integer :: i, j
    real :: q11, q12
    real :: t_q11, t_q12, t_q21, t_q22
    q11 = 0
    q12 = 0
    do i=1,L
        do j=1,L
            q11 = q11 + 2.*(cos(fi(i,j) * PI_MATH / 180.)**2) - 1.
            q12 = q12 + 2.*cos(fi(i,j) * PI_MATH / 180.) * sin(fi(i,j)* PI_MATH / 180.)
        enddo
    enddo
    q11 = q11/float(L*L)
    q12 = q12/float(L*L)
    t_q11 = Q(1,1)
    t_q12 = Q(1,2)
    t_q21 = Q(2,1)
    t_q22 = Q(2,2)
    Q(1,1) = q11
    Q(1,2) = q12
    Q(2,1) = Q(1,2)
    Q(2,2) = -Q(1,1)
    eigens(idx) = eigens(idx) + calc_eigen()
    Q(1,1) = t_q11
    Q(1,2) = t_q12
    Q(2,1) = t_q21
    Q(2,2) = t_q22
end subroutine update_eigens

subroutine init_eigens()
    integer :: i
    do i=1, size(eigens)    
        eigens(i) = 0 
    enddo
end subroutine init_eigens

function calc_eigen() result(eigen)
    !     .. Parameters ..
    INTEGER          N
    PARAMETER        ( N = 2 )
    INTEGER          LDA, LDVL, LDVR
    PARAMETER        ( LDA = N, LDVL = N, LDVR = N )
    INTEGER          LWMAX
    PARAMETER        ( LWMAX = 1000 )
    ! 
    INTEGER          INFO, LWORK
    !
    REAL             VL( LDVL, N ), VR( LDVR, N )
    REAL             WR( N ), WI( N ), WORK( LWMAX )
    REAL eigen
    LWORK = -1
    CALL SGEEV( 'Vectors', 'Vectors', N, Q, LDA, WR, WI, VL, LDVL, &
        VR, LDVR, WORK, LWORK, INFO )
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    CALL SGEEV( 'Vectors', 'Vectors', N, Q, LDA, WR, WI, VL, LDVL, &
        VR, LDVR, WORK, LWORK, INFO )

    ! print *, INFO
    ! print *, WR
    ! print *, WI
    ! print *, MAX(WR(1), WR(2))
    eigen = MAX(WR(1), WR(2))
end function calc_eigen

end program name
