
!***********************************************************************************************
! Example program to test finite differences of a given Order based on  a set of points : x, y
!***********************************************************************************************
program Test_FD

 
  use Finite_differences
  use dislin 
  implicit none 

  integer, parameter :: N= 50; 
  real ::  x(0:N), y(0:N)

  real :: U(0:N), Ux(0:N), Uxx(0:N), dop(0:N) 
  real :: W(0:N, 0:N), Wx(0:N, 0:N), Wxx(0:N, 0:N), Wyy(0:N, 0:N), doperator(0:N, 0:N) 

  integer :: Order = 20, i  
  real :: pi = 4*atan(1.) 

  
  ! *** Define grid nodes 
        call Grid_Initialization( "uniform", "x", Order, x )
        call Grid_Initialization( "nonuniform", "y", Order, y )
 
  ! *** 1D function to be interpolated 
        U =  cos(pi*x); 
  
  ! *** 1D derivatives of U and difference operator: dop 
        call Derivative( "x", 1, U, Ux) 
        call Derivative( "x", 2, U, Uxx) 
  
        dop = Ux * ( U + Uxx )  ! 1D spatial operator 

        call scrmod('revers') 
        call qplot(x, U, N+1)
        call qplot(x, dop, N+1)
        
        
  ! *** 2D function to be interpolated  
        forall(i=0:N) W(i,:) = sin(pi*x(i)) * sin(pi*y)
  
  ! *** partial derivatives of W and difference operator: doperator 
        call Derivative( "x", 1, W, Wx)
        call Derivative( "x", 2, W, Wxx)
        call Derivative( "y", 2, W, Wyy)
  
        doperator = W * Wxx + Wx * Wyy  ! 2D spatial operator 
  
        call qplcon( W, N+1, N+1, 20);
        call qplcon( doperator, N+1, N+1, 20);
 
end program 
