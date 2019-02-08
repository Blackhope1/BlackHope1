!  N_bodies.f90 
!

!****************************************************************************
!
!  PROGRAM: N_bodies
!
!  PURPOSE:  Propagator of a three body system
!
!****************************************************************************

    program N_bodies

use Cauchy_problem
use dislin

    implicit none

    ! Variables
    integer,parameter::N=3, Steps=1000000 !N: Number of bodies
    integer:: i
    real:: t0, tf, Time(0:Steps), U(0:Steps,N*4) ! U = [ x1, y1, x2, y2, ... , u1, v1, u2, v2, ... ]
    
    
    t0=0.
    tf=50000000.
    
    !Initial Conditions U = [ x1, y1, x2, y2, ... , u1, v1, u2, v2, ... ]
    U(0,:) = [0.,0. , 400e6,0. , 200e6,346e6 , 0.,-4.083 , 0.,408.3 , -355.474,205.233] !346e6!-355.474
    
    Time = [ (t0 + (tf -t0 ) * i / (1d0 * Steps), i=0, Steps ) ]
    
    
    

call Cauchy_ProblemS(Time_Domain = Time , Differential_operator = Gravity3, Scheme = Euler , Solution = U )

!PLOTTING

call plots





contains


function Gravity3( U, t ) result(F)
    real :: U(:), t                         ! U = [ x1, y1, x2, y2, ... , u1, v1, u2, v2, ... ]
    real :: F(size(U))
    real :: r(2,size(U)/4) , v(2,size(U)/4) , m(size(U)/4) , drdt(2,size(U)/4) , dvdt(2,size(U)/4) , G
    real :: dij(2),ai(2)
    integer :: N, j, i

    N=size(U)/4

    G = 6.674e-11
    m(1) = 1e24 !Earth 
    m(2) = 1e22 !Moon
    m(3) = 1e20   !Sat


    !State vector U --> Variables
    r(:,1) = U(1:2)
    r(:,2) = U(3:4)
    r(:,3) = U(5:6)
    v(:,1) = U(7:8)
    v(:,2) = U(9:10)
    v(:,3) = U(11:12)



    !Acceleration on i caused by j
    do i=1,N
        ai = 0.
        do j=1,N

            dij = r(:,j)-r(:,i)
		
            if ((dij(1)**2+dij(2)**2)<=1) then
                !write(*,*) "Body ", j, " is too close to ", i, ". Cycling..."
                cycle
            endif
    
            ai= ai + G* m(j)*dij/sqrt(dij(1)**2+dij(2)**2)**3  !Gravity Law
        enddo

        drdt(:,i) = v(:,i)

        dvdt(:,i) = ai(:)



    enddo

    !Accelerations --> derivative of state vector
    F(1:2) = drdt(:,1)
    F(3:4) = drdt(:,2)
    F(5:6) = drdt(:,3)
    F(7:8) = dvdt(:,1)
    F(9:10) = dvdt(:,2)
    F(11:12) = dvdt(:,3)

end function







subroutine plots

real:: m(0:Steps), d(0:Steps,2)
CHARACTER*40 :: CBUF



call qplot( U(0:Steps,1) , U(0:Steps,2) , Steps+1 )   !Movement of Earth

call qplot( U(0:Steps,3) , U(0:Steps,4) , Steps+1 )   !Movement of Moon

call qplot( U(0:Steps,5) , U(0:Steps,6) , Steps+1 )   !Movement of Sat


m = ( U(0:Steps,4) - U(0:Steps,2) )/( U(0:Steps,3) - U(0:Steps,1) )

d(:,2) = abs(U(0:Steps,5) * m - U(0:Steps,6) - m * U(0:Steps,1) + U(0:Steps,2))/(sqrt(m**2+1))

d(:,1) = sqrt(((U(0:Steps,6)-U(0:Steps,2))**2 + (U(0:Steps,5)-U(0:Steps,1))**2)-d(:,2)**2)

call qplot( d(0:Steps,1) , d(0:Steps,2) , Steps+1 )   !Movement of Sat with non inertial reference






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
CALL GRAF(-100e6,8e6,0.,1e6,-4e6,7e6,-4e6,1e6)
CALL TITLE
CALL CHNCRV('LINE')
CALL CURVE( U(0:Steps,1) , U(0:Steps,2) , Steps+1 )

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
!CALL TITLIN('Luna',3)
CALL GRAF(-700e6,500e6,-500e6,250e6,-1000e6,450e6,-1000e6,200e6)
CALL TITLE
CALL CHNCRV('LINE')
CALL CURVE( U(0:Steps,3) , U(0:Steps,4) , Steps+1 )
CALL CURVE( U(0:Steps,1) , U(0:Steps,2) , Steps+1 )
CALL CURVE( U(0:Steps,5) , U(0:Steps,6) , Steps+1 )

CALL LEGINI(CBUF,2,20) ! Legend statements
CALL LEGPOS(NXPOSN(1.),NYPOSN(1.75))
CALL LEGLIN(CBUF,'Luna',1)
CALL LEGLIN(CBUF,'Sat',2)
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
CALL GRAF(-450e6,450e6,-450e6,200e6,-450e6,450e6,-450e6,200e6)
CALL TITLE
CALL CHNCRV('LINE')
CALL CURVE( U(0:Steps,5) , U(0:Steps,6) , Steps+1 )

CALL LEGINI(CBUF,2,20) ! Legend statements
CALL LEGPOS(NXPOSN(1.),NYPOSN(1.75))
CALL LEGLIN(CBUF,'Numerical solution',1)
CALL LEGTIT('Legend')
CALL LEGEND(CBUF,3)

CALL DASH
CALL AXGIT
CALL DISFIN




call METAFL('PDF')
CALL DISINI
CALL LABELS ('EXP','FLOAT')
CALL PAGERA
CALL COMPLX
CALL AXSPOS(450,1800)
CALL AXSLEN(2200,1200)
CALL NAME('X','X')
CALL NAME('Y','Y')
CALL LABDIG(-1,'X')
CALL TICKS(9,'XY')
CALL TITLIN('Solucion para 3 cuerpos',1)
CALL TITLIN('Sat en sistema no inercial',3)
CALL GRAF(185e6,220e6,190e6,10e6,330e6,370e6,330e6,10e6)
CALL TITLE
CALL CHNCRV('LINE')
CALL CURVE( d(0:Steps,1) , d(0:Steps,2) , Steps+1 )

CALL LEGINI(CBUF,2,20) ! Legend statements
CALL LEGPOS(NXPOSN(1.),NYPOSN(220e6))
CALL LEGLIN(CBUF,'Numerical solution',1)
CALL LEGTIT('Legend')
CALL LEGEND(CBUF,3)

CALL DASH
CALL AXGIT
CALL DISFIN


end subroutine

    end program N_bodies

