	    PROGRAM ising

	    IMPLICIT NONE
	    INTEGER ::  MCS
	    INTEGER ::  monte
	    INTEGER ::  L ! size
	    REAL ::     T_start ! reduced temp
	    REAL ::     actual_T
	    REAL ::     T_step
	    INTEGER ::  K_zero ! nb of MCS to ensure equilibrium
	    INTEGER ::  observation_step
	    INTEGER :: temp_iterator
	    INTEGER :: really_all_q_iterator
	    INTEGER, DIMENSION(:,:), ALLOCATABLE :: spins
	    INTEGER, DIMENSION(:), ALLOCATABLE :: NEXT, PREV
	    REAL, DIMENSION(41) :: temps_tab
	    REAL, DIMENSION(-180:180) :: P2
	    REAL, PARAMETER :: Pi = 3.1415927
	    INTEGER :: new_fi
	    REAL, DIMENSION(2,2) :: Q
	    REAL, DIMENSION(164) :: allQ
	    REAL, DIMENSION(:), ALLOCATABLE :: really_all_q
	    
	    integer :: values(1:8), k
        integer, dimension(:), allocatable :: seed
        real rand_num
        integer :: time
	    
	    MCS = 2*10**5
	    L = 20
	    K_zero = 30000
	    observation_step = 100
	    temps_tab = (/ 0.0000001, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,&
	    & 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, &
	    & 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0 /)




	    ALLOCATE(spins(L,L))
	    ALLOCATE(NEXT(L))
	    ALLOCATE(PREV(L))
	    ALLOCATE(really_all_q(SIZE(temps_tab)*(MCS/observation_step)*2))

	    
	    call date_and_time(values=values)
        call random_seed(size=k)
        allocate(seed(1:k))
        seed(:) = values(8)
        call random_seed(put=seed)
	    
	    print *, "fasten seatbelts, time started"
	    print *, "L: ", L, " MCS: ", MCS, "  K_zero: ", K_zero
	    print *, "temps_tab ", temps_tab
	    call tick(time)
	    
	    call generate_TNN()
	    !print *, "generated tnn"
	    
	    call generate_P2()
	    print *, "P2", P2
	    
	    really_all_q_iterator = 1
	    
	    DO temp_iterator=1, SIZE(temps_tab)
	    
	        Q(1,1) = 0
	        Q(1,2) = 0
	        Q(2,1) = 0
	        Q(2,2) = 0
	    
	        actual_T = temps_tab(temp_iterator)
	        print *, "actual temp: ", actual_T, " done in ", tock(time), "s"
	        
	        call generate_spins()
	        !print *, "generated spins"
	        
	        call execute_K_zero_MCS()
	        !print *, "executed K_zero"
	        
	        call execute_all_MCS() ! current values stored there
	        
	        !SGEEV()
	        
	        allQ((temp_iterator-1)*4 + 1) = Q(1,1)
	        allQ((temp_iterator-1)*4 + 2) = Q(1,2)
	        allQ((temp_iterator-1)*4 + 3) = Q(2,1)
	        allQ((temp_iterator-1)*4 + 4) = Q(2,2)
	        
	        print *, spins
	        print *, Q
	        print *, really_all_q_iterator
	        
	    END DO
        
        ! save results
        print *, "writing results to files"
        open(unit=3,file='nlc1_AAAAAAAAAAAAAAAAAAAAAAaa.dat',form='unformatted')
        write(3) really_all_q
        close(3)
        
	    DEALLOCATE(spins)
	    DEALLOCATE(NEXT)
	    DEALLOCATE(PREV)
	    DEALLOCATE(really_all_q)
	    
	    CONTAINS
	    
	    SUBROUTINE execute_all_MCS()
	    REAL :: Q11, Q12
	    integer :: i,j
	        DO monte=1, MCS
!	            IF (MOD(monte,10**5) == 0) THEN
!	                print*, monte, " MCS done in ", tock(time), "s"
!	            END IF
	            
	            call monte_carlo_step()
	            IF (MOD(monte,observation_step) == 0) THEN
	                Q11 = 0
	                Q12 = 0
!	                Q(1,1) = Q(1,1) + SUM(2*cos(REAL(spins))**2 - 1)/(real(L*L))
!	                Q(1,2) = Q(1,2) + SUM(2*cos(REAL(spins)) * sin(REAL(spins)))/(real(L*L))
                    DO i=1, L
	                    DO j=1, L
	                        Q11 = Q11 + real(2)*cos(real(spins(i,j))*Pi/real(180))**2 - real(1)
	                        Q12 = Q12 + real(2)*cos(real(spins(i,j))*Pi/real(180)) * sin(real(spins(i,j))*Pi/real(180))
	                    END DO
	                END DO
!	                print *, Q11, " i dwa ", Q12
	                really_all_q(really_all_q_iterator) = Q11/real(L*L)
	                really_all_q_iterator = really_all_q_iterator + 1
	                really_all_q(really_all_q_iterator) =  Q12/real(L*L)
	                really_all_q_iterator = really_all_q_iterator + 1
	            END IF
	        END DO
	    
	    END SUBROUTINE execute_all_MCS
	    
	    FUNCTION get_energy_delta_after_trial_conf(x,y) result(energy_delta)
	        real :: energy_delta
	        real :: r_p2, l_p2, u_p2, d_p2
	        real :: U0, Ut
	        integer :: x,y
	        
	        r_p2 = P2(spins(x,y) - spins(NEXT(x), y))
	        l_p2 = P2(spins(x,y) - spins(PREV(x), y))
	        u_p2 = P2(spins(x,y) - spins(x, NEXT(y)))
	        d_p2 = P2(spins(x,y) - spins(x, PREV(y)))
	        
	        U0 = -r_p2 - l_p2 - u_p2 - d_p2
	        
	        call random_number(rand_num)
	        new_fi = spins(x,y) + nint((rand_num - 0.5)*10)
	        
	        IF (new_fi > 90) THEN
	            new_fi = new_fi - 180
            ELSE IF (new_fi < -90) THEN
                new_fi = new_fi + 180
            END IF
                

	        r_p2 = P2(new_fi - spins(NEXT(x), y))
	        l_p2 = P2(new_fi - spins(PREV(x), y))
	        u_p2 = P2(new_fi - spins(x, NEXT(y)))
	        d_p2 = P2(new_fi - spins(x, PREV(y)))
	        
	        Ut = -r_p2 - l_p2 - u_p2 - d_p2
	        
	        energy_delta = Ut - U0
	    END FUNCTION get_energy_delta_after_trial_conf
	    
	    FUNCTION should_i_accept_trial_conf(energy_del) result(should_accept)
	        real energy_del, w
	        logical should_accept
	        
	        IF (energy_del < 0) THEN
	            should_accept = .TRUE.
	        ELSE
	            w = EXP(-energy_del/actual_T)
	            call random_number(rand_num)
	            
	            IF (rand_num <= w) THEN
	                should_accept = .TRUE.
	            ELSE
	                should_accept = .FALSE.
	            END IF
            END IF
	    END FUNCTION should_i_accept_trial_conf
	    
	    SUBROUTINE change_spin(x_ch, y_ch)
	        integer :: x_ch, y_ch
	        
	        spins(x_ch,y_ch) = new_fi
	    END SUBROUTINE
	    
	    SUBROUTINE monte_carlo_step()
	        integer :: i,j
	        
	        DO i=1, L
	            DO j=1, L
	                IF (should_i_accept_trial_conf(get_energy_delta_after_trial_conf(i,j))) THEN
	                    call change_spin(i,j)
	                ENDIF
	            END DO
	        END DO
	    END SUBROUTINE
	    
	    SUBROUTINE execute_K_zero_MCS()
	        integer :: z
	        
	        DO z=1, K_zero
	            call monte_carlo_step()
	        END DO
	    END SUBROUTINE
	    
	    SUBROUTINE generate_spins()
	        integer :: i,j,rand_fi
	        
	        DO i=1, L
	            DO j=1, L
	                call random_number(rand_num)
	                rand_fi = MOD(nint(rand_num*1000),180)-90
	                spins(i,j) = rand_fi
!                    spins(i,j) = 0
	            END DO
	        END DO
	    END SUBROUTINE
	    
	    SUBROUTINE generate_TNN()
	        integer :: i
	        
	        DO i=1, L
	            NEXT(i) = i+1
	            PREV(i) = i-1
	        END DO
	        
	        PREV(1) = L
	        NEXT(L) = 1
	    END SUBROUTINE
	    
	    SUBROUTINE generate_P2()
	        integer :: i
	        
	        DO i=-180, 180
	            P2(i) = 3/real(2) * cos(real(i)*Pi/real(180))**2 - 0.5
	        END DO
	    END SUBROUTINE
	    
	    subroutine tick(t)
            integer, intent(OUT) :: t
            call system_clock(t)
        end subroutine tick

        ! returns time in seconds from now to time described by t 
        real function tock(t)
            integer, intent(in) :: t
            integer :: now, clock_rate
            call system_clock(now,clock_rate)
            tock = real(now - t)/real(clock_rate)
        end function tock
        
	    END PROGRAM ising
