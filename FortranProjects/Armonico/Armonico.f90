!  Armonico.f90 
!
!

!****************************************************************************
!
!  PROGRAM: Armonico
!
!  PURPOSE:  Numerical solutions for tha harmonic movement.
!
!****************************************************************************

    program Armonico

use cauchy_problem
use dislin

    implicit none

    ! Variables
    integer,parameter::N=100
    integer:: i
    real:: t0, tf, Time(0:N), U(0:N,2), error (0:N)
    CHARACTER*40 :: CBUF
    
	t0=0.
	tf=10.
	
    U(0,:)=[1.,0.]
    
    
    Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
    
    
    

call Cauchy_ProblemS(Time_Domain = Time , Differential_operator = Harmonic_Equation, Scheme = Inverse_Euler , Solution = U )

error = abs(U(:,1)-[(cos(t0 + (tf -t0 ) * i / (1d0 * N)), i=0, N )])







!GENERATING PLOTS IN PDF
call plots




contains

function Harmonic_Equation( U, t ) result(F)
real :: U(:), t
real :: F(size(U))
F(1) = U(2)
F(2) = -U(1)
end function


















subroutine plots


call METAFL('PDF')
CALL DISINI
CALL PAGERA
CALL COMPLX
CALL AXSPOS(450,1800)
CALL AXSLEN(2200,1200)
CALL NAME('Time','X')
CALL NAME('X','Y')
CALL LABDIG(-1,'X')
CALL TICKS(9,'XY')
CALL TITLIN('Solution to the',1)
CALL TITLIN('Harmonic Equation',3)
CALL GRAF(0.,10.,0.,1.,-2.,2.,-2.,0.5)
CALL TITLE
CALL CHNCRV('LINE')
CALL CURVE(Time(0:N),U(0:N,1),N+1)
CALL CURVE(Time(0:N),[(cos(t0 + (tf -t0 ) * i / (1d0 * N)), i=0, N )],N+1)

CALL LEGINI(CBUF,2,20) ! Legend statements
CALL LEGPOS(NXPOSN(1.),NYPOSN(1.75))
CALL LEGLIN(CBUF,'Numerical solution',1)
CALL LEGLIN(CBUF,'Analytical solution',2)
CALL LEGTIT('Legend')
CALL LEGEND(CBUF,3)

CALL DASH
CALL AXGIT
CALL DISFIN


CALL DISINI
CALL PAGERA
CALL COMPLX
CALL AXSPOS(450,1800)
CALL AXSLEN(2200,1200)
CALL NAME('Time','X')
CALL NAME('X','Y')
CALL LABDIG(-1,'X')
CALL LABDIG(6,'Y')
CALL TICKS(9,'XY')
CALL TITLIN('Numerical error of the',1)
CALL TITLIN('Harmonic Equation',3)
CALL GRAF(0.,10.,0.,1.,minval(error),maxval(error),0.,0.1)
CALL TITLE
CALL CHNCRV('LINE')
CALL CURVE(Time(0:N),error,N+1)

CALL LEGINI(CBUF,2,20) ! Legend statements
CALL LEGPOS(NXPOSN(1.),NYPOSN(1.75))
CALL LEGLIN(CBUF,'Numerical solution',1)
CALL LEGLIN(CBUF,'Analytical solution',2)
CALL LEGTIT('Legend')
CALL LEGEND(CBUF,3)

CALL DASH
CALL AXGIT
CALL DISFIN




end subroutine



    end program Armonico

