module Roots
    
    implicit none  ! no implicit specifications of variables
                   ! real, integer of any variable type must be given 
    
contains     
!*****************************************************************************
!* Roots a second grade equation:
!       a x**2 + b x + c = 0,    with complex numbers a,b,c 
!*****************************************************************************
 subroutine Roots_2th 
 
    complex :: a, b, c, d  
    complex :: x1, x2 
    
    a = (1., 1.) 
    b = (0., 1.) 
    c = (0., 0.) 
    
    if (abs(a)==0) then 
         if (abs(b)==0) then 
             write(*,*) "There is no solution " 
         else
             write(*,*) "There is only one solution x1 =", -c/b  
         end if 
    else
            d = sqrt( b**4 - 4*a*c ) 
            x1 = ( -b + d )/( 2*a ) 
            x2 = ( -b - d )/( 2*a ) 
            write(*,*) " x1 = ", x1 
            write(*,*) " x2 = ", x2 
    end if 
     
 end subroutine 
 
 end module 