module Series_expansion 
    
    implicit none 
    real :: PI = 4 * atan(1d0) ! PI definition 
    
    abstract interface 
      real function ak(x, k) 
         real, intent(in) :: x 
         integer, intent(in) :: k 
      end function
    end interface 
    
    
    contains

!*****************************************************************************
!* Examples of truncated and infinite expansions 
!***************************************************************************** 
 subroutine Expansion_examples 
  
   write(*,*) "truncated exp. sin(x)=",  expansion( a_k = ak_sin,  M = 8,   x0 = 0., x = PI/2 ) 
   write(*,*) "truncated exp. 1/(1+x)=", expansion( a_k = ak_frac, M = 10,  x0 = 0., x = 0.9 ) 
   
   write(*,*) "infinity exp. sin(x)=",  expansion( a_k = ak_sin,  x0 = 0., x = PI/2 ) 
   write(*,*) "infinity exp. 1/(1+x)=", expansion( a_k = ak_frac, x0 = 0., x = 0.9 ) 
   
 
 end subroutine    
!*************************************************************************** 
!                    Taylor expansion 
! expansion = sum _{k=0} ^{M} fk(x0) (x-x0)**k if M is given  
!***************************************************************************
real function expansion( a_k, M, x0, x)  
       procedure (ak) :: a_k               ! generic coefficent of the expasion a_k 
       integer, optional, intent(in) :: M  ! optional truncation order. 
                                           ! if it is not present, it assumes infinity  
       real, intent(in) :: x0              ! origin of Taylor expansion 
       real, intent(in) :: x               ! independent variable where the expansion is evaluated 
    
    integer :: k       ! index of the expansion 
    integer ::  Mmax   ! max value of terms if no convergence is reached to avoid an infinite loop
    real :: b_k        ! value of  a_k(x0, k) * ( x-x0)**k
    real ::  S         ! sum 
    
    if (present(M)) then 
                            Mmax = M  
    else
                            Mmax = 10000  
    end if 
    
    S = 0; k = 0;  b_k = 1  
       
    do while( ( (abs(b_k) > 1d-14) .or. (b_k == 0)) .and. (k <= Mmax) )  
           
        b_k = a_k(x0, k) * (x-x0)**k  
        S = S + b_k  
        k = k + 1 
    end do 
    
    expansion = S 
 
 end function 
   
!  k-th derivative of sin(x)  
 real function ak_sin(x, k) 
      real, intent(in) :: x 
      integer, intent(in) :: k 
      
      ak_sin = sin( x + PI*k/2 ) 
          
 end function   
 
! k-th derivative of 1/(1+x)   
 real function ak_frac(x, k) 
      real, intent(in) :: x 
      integer, intent(in) :: k 
    
      ak_frac = 1 / ( 1 - x )**(k+1)
  
 end function
 
! factorial(n) 
 real function factorial(n) 
    integer, intent(in) :: n 
    
     integer :: i ! index of the product 
     real ::  p   ! partial product of factorial 
    
     p = 1 
     
     do i = 1, n 
                  p = p * i 
     end do 
     factorial = p 
  
 end function 
 

end module   