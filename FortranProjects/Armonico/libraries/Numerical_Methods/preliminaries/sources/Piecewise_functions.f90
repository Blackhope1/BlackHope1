module Piecewise_functions 
    
    use dislin                 ! scientific graphics library
    implicit none              ! implicit specification of variables is not allowed 
    
    real :: PI = 4 * atan(1d0) ! PI definition 
   
    
    contains
!*****************************************************************************
! An example of  piecewise functions 
!***************************************************************************** 
elemental real function f(x) 
real, intent(in) :: x 
 

 if (x<= - PI ) then 
     
        f =  - (x+PI) / (x**2 + 1) 
        
 elseif ( x < PI ) then 
     
        f = x * sin(  x )
 else
        f =   (x-PI) / (x**2 + 1) 
 endif 

end function 
!*****************************************************************************
! It plots  piecewise functions 
!*****************************************************************************   
subroutine Function_examples 

  real, allocatable :: x(:)  ! domain point 
  real, allocatable :: y(:)  ! image point 
  integer :: i               ! index of domain point 
  integer :: N               ! number of points 
  real :: a=-5, b=5          ! domain segment  x in [a, b]  
  real ::  dx                ! step between domin points 
    
  N = 10 
  dx = (b-a)/N 
  
  ! automatic allocation 
  x = [ (a + i * dx , i=0, N) ] 
  
  ! vectorial operation (note: elemntal real function f(x)) 
  y = f(x)
  write(*,*) "sum from 0..10 = ", dx *sum(y) 
  
  N = 100 
  dx = (b-a)/N 
  x = [ (a + i * dx , i=0, N) ] 
  y = f(x)
 
  write(*,*) "sum from 0..100 = ", dx *sum(y)
 
  
  N = 1000 
  dx = (b-a)/N 
  x = [ (a + i * dx , i=0, N) ] 
  y = f(x) 
  
  write(*,*) "sum from 0..1000 = ",  dx *sum(y)  
 
  ! plot y= f(x)  x in [a,b]
  call scrmod("reverse")
  call qplot(x, y, N+1) 

end subroutine 

 
end module   