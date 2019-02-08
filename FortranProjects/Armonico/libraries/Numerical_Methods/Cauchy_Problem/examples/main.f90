


!***********************************************************************************************
! Example program to test the temporal integration of the Cauchy Problem 
!***********************************************************************************************
program Test_Cauchy_Problems 

    use Cauchy_Problem

    implicit none
 
  
  call Test_mass_spring_damper
  
  call Test_Lorenz_Attractor
  
  call Test_Lorenz_Attractor_implicit
  
 

contains 
!-----------------------------------------------------------------------
function  F_spring( U, t )  result(F) 

    real :: U(:), t 
    real :: F(size(U)) 
            
    real, parameter :: M = 1.0 , D = .1 , P = 0, K = 10.0;
           
        F(1) = U(2)
        F(2) = ( P - D * U(2) -K * U(1) )/ M

 end function 
!-----------------------------------------------------------------------
subroutine Test_mass_spring_damper
        
    real :: t0 = 0, tf = 5
    integer :: i 
    integer, parameter :: N = 100;   !Time steps
    real :: Time (0:N), U(0:N, 2); 


    Time = [ (t0 + (tf -t0 )*i/(1d0*N), i=0, N ) ]
        
    U(0,:) = [ 1, 0 ] 
        
   

    call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = F_spring, Solution = U )

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TITLIN (" Mass-spring-damper system  ", 2)
    call qplot( time(0:N), U(0:N, 1), N+1 );
             


end subroutine
 !----------------------------------------------------------------------
function  Lorenz_attractor(U, t)  result(F) 

     real  ::  U(:), t  
     real   :: F( size(U) )

    real :: a=10., b=28., c=2.6666666666
    real :: x, y , z
      
     x = U(1); y = U(2); z = U(3)   

     F(1) =   a * ( y - x ) 
     F(2) =   x * ( b - z ) - y 
     F(3) =   x * y - c * z   

end function  
!-----------------------------------------------------------------------
subroutine Test_Lorenz_Attractor


    real :: t0 = 0, tf = 30
    integer :: i 
    integer, parameter :: N = 2000;   !Time steps
    real :: Time (0:N), U(0:N, 3); 


    Time = [ (t0 + (tf -t0 )*i/(1d0*N), i=0, N ) ]
        
    U(0,:) = [ 12, 15, 30 ] 

    call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = Lorenz_attractor, Solution = U )
    
    
    call disini() 
    CALL TITLIN (" Lorenz atractor  ", 2)
    call qplot( U(:, 1), U(:, 2), N+1 ); 

end subroutine 
!-----------------------------------------------------------------------
subroutine Test_Lorenz_Attractor_implicit


    real :: t0 = 0, tf = 30
    integer :: i 
    integer, parameter :: N = 2000;   !Time steps
    real :: Time (0:N), U(0:N, 3); 


    Time = [ (t0 + (tf -t0 )*i/(1d0*N), i=0, N ) ]
        
    U(0,:) = [ 12, 15, 30 ] 

    call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = Lorenz_attractor, Solution = U,Scheme = Crank_Nicolson )
    
    
    call disini() 
    CALL TITLIN (" Lorenz atractor  ", 2)
    CALL TITLIN (" Crank Nicolson temporal scheme  ", 4)
    call qplot( U(:, 1), U(:, 2), N+1 ); 

end subroutine         
       
       

end program 
  

   
