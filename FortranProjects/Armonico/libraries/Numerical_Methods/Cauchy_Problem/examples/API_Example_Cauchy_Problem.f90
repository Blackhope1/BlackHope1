module API_Example_Cauchy_Problem

    use Cauchy_Problem
    implicit none
    
    real, parameter :: PI = 4 *atan(1d0) 
    
    
    contains  
  
subroutine Cauchy_problem_examples
   
   call First_order_ODE
   
   call Linear_Spring
   
   call Lorenz_Attractor
   
end subroutine   

! line 19 
subroutine Linear_Spring 
        
    real :: t0 = 0, tf = 5
    integer :: i 
    integer, parameter :: N = 100   !Time steps
    real :: Time (0:N), U(0:N, 2)

    Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
        
    U(0,:) = [ 5, 0] 

    call Cauchy_ProblemS( Time_Domain = Time ,                     & 
                          Differential_operator = F_spring,        & 
                          Scheme = Crank_Nicolson , Solution = U )
    
    call qplot( time(0:N) , U(0:N,1) , N+1 )
      
contains

function F_spring( U, t )  result(F) 

    real :: U(:), t 
    real :: F(size(U)) 
            
    real, parameter :: a = 3.0,b = 0.0 
    
    F(1) = U(2)
    F(2) = -a * t * U(1) + b
     
end function

end subroutine





! line 59
subroutine First_order_ODE
        
    real :: t0 = 0, tf = 4
    integer :: i 
    integer, parameter :: N = 1000  !Time steps
    real :: Time (0:N), U(0:N,1)

    Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]

    U(0,1) =  1
   
    call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = F1, & 
                          Scheme = Runge_Kutta4, Solution = U )
    
    call qplot( Time(0:N) , U(0:N,1) , N+1 )
       
contains

function F1( U, t ) result(F) 
    real :: U(:), t 
    real :: F(size(U)) 
    
    F(1) = -2*U(1)
  
end function 


end subroutine



! line 94
subroutine Lorenz_Attractor
   integer, parameter :: N = 10000
   real :: Time(0:N), U(0:N,3) 
   real :: a=10., b=28., c=2.6666666666
   real :: t0 =0, tf=25 
   integer :: i
 
   Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
   
   U(0,1) = 12 
   U(0,2) = 15 
   U(0,3) = 30  

   
    call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = F_L, & 
                          Scheme = Runge_Kutta4, Solution = U )
    
    call qplot(  U(0:N,1) , U(0:N,2) , N+1 ) 
    
   
contains

 function F_L(U, t) result(F)           
     real :: U(:),t
     real :: F(size(U))
     
     real :: x, y , z
      
     x = U(1); y = U(2); z = U(3)   

     F(1) =   a * ( y - x ) 
     F(2) =   x * ( b - z ) - y 
     F(3) =   x * y - c * z   

end function
  


end subroutine

subroutine Lorenz_Attractor_Cash_Karp
   integer, parameter :: N = 10000
   real :: Time(0:N), U(0:N,3) 
   real :: a=10., b=28., c=2.6666666666
   real :: t0 =0, tf=25 
   integer :: i
 
   Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
   
   U(0,1) = 12 
   U(0,2) = 15 
   U(0,3) = 30  

   
    call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = F, & 
                          Scheme = Cash_Karp, Solution = U )
    
    call qplot(  U(0:N,1) , U(0:N,2) , N+1 ) 
    
   
contains

 function  F(U, t)            
     real :: U(:),t
     real :: F(size(U))
     
     real :: x, y , z
      
     x = U(1); y = U(2); z = U(3)   

     F(1) =   a * ( y - x ) 
     F(2) =   x * ( b - z ) - y 
     F(3) =   x * y - c * z   

 end function
 

end subroutine

!******************************************************
!* Test Floquet_multipliers 
!******************************************************
subroutine Test_floquet_multipliers 

      integer, parameter :: N = 4, M=10000 
      complex :: rho(N)
      real :: Period, Phi(N,N), PM(N,N), Up(N)    
      integer :: i, j 
      
     Period =   2 * PI 
     
     Up = [ 1, 0, 0, 1 ] 

     call Floquet_multipliers( F, U0, Period, rho, Phi) 

     do i=1, N 
         write(*,*)  " abs(rho)=", abs(rho(i))  
     end do 
     
     write(*,'(a20, <N>f13.2)')  " U(0) = ",  Up  
     PM = transpose(Phi) 
     !do i=1, N 
     !      write(*,'(a20, <N>f13.2)') " Fundamental matrix = ", PM(i,:)   
     !end do
     !write(*,'(a20, <N>f13.2)') " U(nT) = ", matmul(PM, Up ) 
     do j=1, M 
         PM = matmul( PM, transpose(Phi))  
     end do 
     write(*,'(a20, <N>f13.2)') " U(nT) = ", matmul(PM, Up ) 
     write(*,*) 
     do i=1, N 
           write(*,'(a20, <N>f13.2)') " Fundamental matrix = ", PM(i,:)   
     end do 
          
     
contains
!---------------
! Kepler orbits 
!---------------
function F(U, t)  
  real :: U(:), t
  real :: F ( size(U) ) 
    
    real :: a = 0.15
     
    F = [ U(3), U(4), - U(1)/ norm2(U(1:2))**3, - U(2)/ norm2(U(1:2))**3 ] 
    
    
end function  
!-------------------
! Periodic solution 
!-------------------
function U0(N, t)  
  integer, intent(in) :: N 
  real, intent(in) ::   t
  real :: U0 ( N ) 
    
     real :: x, y, dx, dy 
     
     x = cos(t); dx = - sin(t) 
     y = sin(t); dy =   cos(t) 
     
     U0 = [ x, y, dx, dy  ] 
     
    
end function  


end subroutine 

!******************************************************
!* Test Floquet_multipliers 
!******************************************************
subroutine Test_floquet_multipliers2

      integer, parameter :: N = 3 
      complex :: rho(N)
      real :: Period, Phi(N,N)  
      integer :: i 
      
     Period =   2 * PI 

     call Floquet_multipliers( F, U0, Period, rho, Phi ) 

     do i=1, N 
         write(*,'(a20, 2f8.3, a10, f8.3)') " Floquet multiplier = ", rho(i), " abs(rho)=", abs(rho(i))  
     end do 
     
     
contains
!---------------
! Chapter 9 Stability II: maps and periodic orbits 
!---------------
function F(U, t)  
  real :: U(:), t
  real :: F ( size(U) ) 
    
    real :: a = -0.15, x, y, z
    
    x = U(1); y = U(2); z= U(3) 
     
    F = [ -y + x * (1 - x**2 - y**2), x + y * (1 - x**2 - y**2), a*z ] 
    
    
end function  
!-------------------
! Periodic solution 
!-------------------
function U0(N, t)  
  integer, intent(in) :: N 
  real, intent(in) ::   t
  real :: U0 ( N ) 
    
     real :: x, y
     
     x = cos(t); 
     y = sin(t); 
     
     U0 = [ x, y, 0.  ] 
     
    
end function  


end subroutine 




!***********************************************************************
!*  If gives the linearized operator F(U,t) in some point U0 
!***********************************************************************
subroutine test_System_matrix
          
          integer, parameter :: N = 4 
          real :: A(N,N),  U0(N)   
          integer :: j  
          real :: x, y, dx ,dy, t 
                  
          t = 0
          x = cos(t) 
          y = sin(t) 
          dx = -sin(t) 
          dy = cos(t) 
          
          U0  = [ x, y, dx, dy ]
          call System_matrix( F, U0, t, A) 
          
          do j=1, N 
             write(*,'(10f8.3)') A(j, :) 
          end do 
          
        
contains 
!----------------------
! Kepler force 
!----------------------
function F(U, t) 
  real :: U(:), t 
  real :: F( size(U) ) 

  F = [ U(3), U(4), -U(1) / norm2( U(1:2) )**3 ,  -U(2) / norm2( U(1:2) )**3 ] 
   
end function           
             
end subroutine  
 
 





end module
