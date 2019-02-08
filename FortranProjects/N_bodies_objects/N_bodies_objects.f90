!  N_bodies_objects.f90 
!
!

!****************************************************************************
!
!  PROGRAM: N_bodies
!
!  PURPOSE:  With Object oriented programming. Only 3 bodies
!
!****************************************************************************

    program N_bodies_objects

    use dislin
    use types
    
    implicit none

    integer,parameter :: Steps=100000
    real:: pos(0:steps,3,2), Time(0:Steps), t0, tf
    integer::i, k
    type (planet) :: Earth, Moon, Sat
    type (system) :: Earth_SOI !Sphere of Influence
    
    t0=0.
    tf=50000000.
    print*,"1"
    Time = [ (t0 + (tf -t0 ) * i / (1d0 * Steps), i=0, Steps ) ]
    print*,"2"
    !INITIAL CONDITIONS
    Earth % r = [ 0., 0. ]
    Moon % r = [ 400e6, 0. ]
    Sat % r = [ 200e6, 346e6 ]
    Earth % v = [ 0.,-4.08]
    Moon % v = [ 0., 408.3]
    Sat % v = [ -355.474, 205.233 ]
    Earth % m = 1e24
    Moon % m = 1e22
    Sat % m = 0.
    
    Earth_SOI % bodies = [ Earth, Moon, Sat ]
    
    pos(0,:,1) = Earth_SOI % bodies % r(1) !pos(time, body, coordinate)
    pos(0,:,2) = Earth_SOI % bodies % r(2)
    
    do k=1,steps
        
        call update_positions(Earth_SOI, time(k)-time(k-1))
        
        pos(k,:,1) = Earth_SOI % bodies % r(1)
        pos(k,:,2) = Earth_SOI % bodies % r(2)
    enddo
    
    
    call plots
    
    
    
    contains
    
    
    
    !Function of the positions  of the bodies of a certain system, it returns an object that contains the accelerations and velocities
    function gravity3 (system_SOI) result (F) 
    
        real, parameter:: G = 6.67e-11
        type(system) :: system_SOI
        type(gravity_planet):: F(size(system_SOI % bodies(:)))
        integer:: N , j
        real:: dij(2)
        
        N = size(system_SOI % bodies(:))
    
    
    do i=1,N
        F(i) % ai = 0.
        do j=1,N

            dij = system_SOI % bodies(j) % r- system_SOI % bodies(i) % r
		
            if ((dij(1)**2+dij(2)**2)<=100) then
                !write(*,*) "Body ", j, " is too close to ", i, ". Cycling..."
                cycle
            endif
    
            F(i) % ai= F(i) % ai + G * system_SOI % bodies(j) % m * dij/sqrt(dij(1)**2+dij(2)**2)**3
        enddo
    
    
        F(i) % drdt = system_SOI % bodies(i) % v

        F(i) % dvdt = F(i) % ai
        

    enddo
    
    
    end function
    
    
    subroutine update_positions(system_SOI, time_step) !Euler propagator
        type(system), intent(inout) :: system_SOI
        real,intent(in) :: time_step
        type(gravity_planet) :: instant_gravity(3)
        integer::j
        
        instant_gravity = gravity3(system_SOI)
        
        do j=1,3
            
            system_SOI % bodies(j) % v = system_SOI % bodies(j) % v + instant_gravity(j) % dvdt * time_step
            system_SOI % bodies(j) % r = system_SOI % bodies(j) % r + instant_gravity(j) % drdt * time_step
        enddo
           
    end subroutine
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    subroutine plots
    
        CHARACTER*40 :: CBUF
    
    call qplot(pos(:,1,1), pos(:,1,2), steps+1)
    call qplot(pos(:,2,1), pos(:,2,2), steps+1)
    call qplot(pos(:,3,1), pos(:,3,2), steps+1)
    
    
    
call METAFL('PDF')
CALL LABELS ('EXP','FLOAT')
CALL DISINI
CALL PAGERA
CALL COMPLX
CALL AXSPOS(450,1800)
CALL AXSLEN(2200,1200)
CALL NAME('X','X')
CALL NAME('Y','Y')
CALL LABDIG(-1,'X')
CALL TICKS(9,'XY')
CALL TITLIN('Solucion para 3 cuerpos',1)
CALL TITLIN('Tierra',3)
CALL GRAF(-5e5,85e5,0.,1e6,-45e5,45e5,-4e6,1e6)
CALL TITLE
CALL CHNCRV('LINE')
CALL CURVE( pos(:,1,1), pos(:,1,2), Steps+1 )

CALL LEGINI(CBUF,2,20) ! Legend statements
CALL LEGPOS(NXPOSN(1.),NYPOSN(1.75))
CALL LEGLIN(CBUF,'Numerical solution',1)
CALL LEGTIT('Legend')
CALL LEGEND(CBUF,3)

CALL DASH
CALL AXGIT
CALL DISFIN




call METAFL('PDF')
CALL LABELS ('EXP','FLOAT')
CALL DISINI
CALL PAGERA
CALL COMPLX
CALL AXSPOS(450,1800)
CALL AXSLEN(2200,1200)
CALL NAME('X','X')
CALL NAME('Y','Y')
CALL LABDIG(-1,'X')
CALL TICKS(9,'XY')
CALL TITLIN('Solucion para 3 cuerpos, ',1)
CALL TITLIN('Luna',3)
CALL GRAF(-450e6,450e6,-450e6,300e6,-450e6,450e6,-450e6,300e6)
CALL TITLE
CALL CHNCRV('LINE')
CALL CURVE( pos(:,2,1), pos(:,2,2), Steps+1 )

CALL LEGINI(CBUF,2,20) ! Legend statements
CALL LEGPOS(NXPOSN(4e6),NYPOSN(1.75))
CALL LEGLIN(CBUF,'Numerical solution',1)
CALL LEGTIT('Legend')
CALL LEGEND(CBUF,3)

CALL DASH
CALL AXGIT
CALL DISFIN



call METAFL('PDF')
CALL LABELS ('EXP','FLOAT')
CALL DISINI
CALL PAGERA
CALL COMPLX
CALL AXSPOS(450,1800)
CALL AXSLEN(2200,1200)
CALL NAME('X','X')
CALL NAME('Y','Y')
CALL LABDIG(-1,'X')
CALL TICKS(9,'XY')
CALL TITLIN('Solucion para 3 cuerpos',1)
CALL TITLIN('Sat',3)
CALL GRAF(-450e6,450e6,-450e6,300e6,-450e6,450e6,-450e6,300e6)
CALL TITLE
CALL CHNCRV('LINE')
CALL CURVE( pos(:,3,1), pos(:,3,2), Steps+1 )

CALL LEGINI(CBUF,2,20) ! Legend statements
CALL LEGPOS(NXPOSN(1.),NYPOSN(1.75))
CALL LEGLIN(CBUF,'Numerical solution',1)
CALL LEGTIT('Legend')
CALL LEGEND(CBUF,3)

CALL DASH
CALL AXGIT
CALL DISFIN



call METAFL('PDF')
CALL LABELS ('EXP','FLOAT')
CALL DISINI
CALL PAGERA
CALL COMPLX
CALL AXSPOS(450,1800)
CALL AXSLEN(2200,1200)
CALL NAME('X','X')
CALL NAME('Y','Y')
CALL LABDIG(-1,'X')
CALL TICKS(9,'XY')
CALL TITLIN('Solucion para 3 cuerpos',1)
CALL TITLIN('Sistema Tierra-Luna',3)
CALL GRAF(-450e6,450e6,-450e6,300e6,-450e6,450e6,-450e6,300e6)
CALL TITLE
CALL CHNCRV('LINE')
CALL CURVE( pos(:,1,1), pos(:,1,2), Steps+1 )
CALL CURVE( pos(:,2,1), pos(:,2,2), Steps+1 )

CALL LEGINI(CBUF,2,20) ! Legend statements
CALL LEGPOS(NXPOSN(150e6),NYPOSN(1.75))
CALL LEGLIN(CBUF,'Tierra',1)
CALL LEGLIN(CBUF,'Luna',2)
CALL LEGTIT('Legend')
CALL LEGEND(CBUF,3)

CALL DASH
CALL AXGIT
CALL DISFIN



    
    
    end subroutine
    
    end program 


