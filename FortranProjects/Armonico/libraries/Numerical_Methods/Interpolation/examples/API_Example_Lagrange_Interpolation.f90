module API_Example_Lagrange_Interpolation

    use Interpolation
    
    implicit none
    real :: PI = 4 * atan(1d0)
    
    
    contains
    
  !line 11  
subroutine Lagrange_Interpolation_examples

     call Interpolated_value_example 
     call Interpolant_example 
     call Integral_example 
     
     call Lagrange_polynomial_example 
     call Ill_posed_interpolation_example
     
     call Lebesgue_and_PI_functions
     
     call Chebyshev_polynomials
     call Chebyshev_expansion( t_kind = 1 ) 
     call Chebyshev_expansion( t_kind = 2 )
   
end subroutine  



!*************************************************************************************
! It interpolates a function at a given abscisa x= xp for a given set of values:
!       (xj, fj) j=0, ..., N  
!line 34******************************************************************************
subroutine Interpolated_value_example 

    real :: x(4) = [ 0.0, 0.1, 0.2, 0.5 ] 
    real :: f(4) = [ 0.3, 0.5, 0.8, 0.2 ] 
    
    real :: xp = 0.15 
    real :: yp1, yp2  
      
   
    yp1 = Interpolated_value( x , f , xp )
    yp2 = Interpolated_value( x , f , xp, 3)  
    
    write (*,*) 'The interpolated value (order 2)  at xp is = ', yp1
    write (*,*) 'The interpolated value (order 3)  at xp is = ', yp2
   
end subroutine
 
  
 



!********************************************************************
!*  It builds an interpolant for a given set of values 
!line 59 ************************************************************
subroutine Interpolant_example

 integer, parameter :: N=3, M=400 
 real ::  xp(0:M) 
 real :: x(0:N) = [ 0.0, 0.1, 0.2, 0.5 ] 
 real :: f(0:N) = [ 0.3, 0.5, 0.8, 0.2 ]
 real :: I_N(0:N, 0:M)  ! first index: derivative
                        ! second index: point where the interpolant is evaluated 
 real :: a, b 
 integer :: i 
 
   a = x(0); b = x(N) 
   xp = [ (a + (b-a)*i/M, i=0, M) ] 

   I_N = Interpolant(x, f, N, xp) 
 
   call scrmod("reverse")
   call qplot( xp, I_N(0, :), M+1) ! plot the Interplant 
   call qplot( xp, I_N(1, :), M+1) ! first derivative of the interpolant 
 
 
end subroutine    
    
 



!***************************************************************************************
! Integral_example
!line 89************************************************************************************************************************************
subroutine Integral_example
 
 integer, parameter :: N=6
 real :: x(0:N), f(0:N) 
 real :: a = 0, b = 1, I0    
 integer :: i 
 
   x = [ (a + (b-a)*i/N, i=0, N) ] 
   f = sin ( PI * x ) 

   I0 = Integral( x , f, 4 )
   write (*,*) 'The integral [0,1] of sin( PI x ) is: ', I0
   write (*,*) 'Error : ',  ( 1 -cos(PI) )/PI - I0
   
end subroutine






!**************************************************************************************************
! Lagrange_polynomial_example
!* line 114****************************************************************************************
subroutine Lagrange_polynomial_example
 
 integer, parameter :: N=4, M=400 
 real :: x(0:N), f(0:N), xp(0:M) 
 real :: Lg(-1:N, 0:N, 0:M)   ! Lagrange polynomials 
                              ! -1:N (integral, lagrange, derivatives) 
                              !  0:N ( L_0(x), L_1(x),... L_N(x)     ) 
                              !  0:M ( points where L_j(x) is evaluated  ) 
 real :: Lebesgue_N(-1:N, 0:M) 
                            
 integer :: i 
 real :: a=-1, b=1 
 
 x  = [ (a + (b-a)*i/N, i=0, N) ] 
 f = 1 
 xp = [ (a + (b-a)*i/M, i=0, M) ] 
 
 do i=0, M   
   Lg(:, :, i) = Lagrange_polynomials( x, xp(i) )
 end do 
 
 Lebesgue_N = Lebesgue_functions( x, xp ) 
   
end subroutine

 !call scrmod("reverse")
 !call metafl("xwin")
 !call disini  
 !call graf(-1., 1., -1., 0.2, -1., 2., -1., 0.2); 
 !do i=0, N 
 !          call curve(xp, Lg(0,i, :), M+1) 
 !end do 
 !
 !call color("blue")
 !call marker(21)
 !call incmrk(-1) 
 !call curve(x, f, N+1 )
 !call disfin 


!********************************************************************
!* Ill_posed_interpolation_example 
!**line 159 *********************************************************
subroutine Ill_posed_interpolation_example 
 
 integer, parameter :: N=64, M=300
 real :: x(0:N), f(0:N)
 real :: I_N(0:N, 0:M)          
 real :: Lebesgue_N(-1:N, 0:M) 
 real :: xp(0:M)   
 real :: a=-1, b=1
 integer :: i  
 
 x  = [ (a + (b-a)*i/N, i=0, N) ] 
 xp = [ (a + (b-a)*i/M, i=0, M) ] 
 f = sin ( PI * x ) 
 
 I_N = Interpolant(x, f, N, xp) 
 Lebesgue_N = Lebesgue_functions( x, xp ) 
 
 call scrmod("reverse")
 call qplot( xp, I_N(0, :), M+1) 
 write(*,*) " maxval Lebesgue =", maxval( Lebesgue_N(0,:) )  
 
end subroutine 





!********************************************************************
!* Lebesgue and PI functions
!**line 189 *********************************************************
subroutine Lebesgue_and_PI_functions 
 
 integer, parameter :: N=10, M=700
 real :: x(0:N), xp(0:M)
 real :: Lebesgue_N(-1:N, 0:M),  PI_N(0:N, 0:M) 
 real :: a=-1, b=1
 integer :: i, k  
 
  x  = [ (a + (b-a)*i/N, i=0, N) ] 
  xp = [ (a + (b-a)*i/M, i=0, M) ] 
 
  Lebesgue_N = Lebesgue_functions( x, xp ) 
  PI_N = PI_error_polynomial( x, xp ) 
 
  call scrmod("reverse")
  
  do k=0, 2      ! k=0 function, k=1 first derivative, k=2 second derivative
    call qplot( xp, Lebesgue_N(k, :), M+1) 
    call qplot( xp, PI_N(k, :), M+1) 
  end do  
 
end subroutine 









 
!**************************************************************************************************
!
!***line 224***************************************************************************************
subroutine Chebyshev_polynomials
 
    integer, parameter :: N = 100 
    integer, parameter :: M = 5 
    real :: x(0:N), y1(0:N), y2(0:N)
    real :: x0=-0.99, xf= 0.99 
    integer :: i, k
    
    call scrmod("reverse") 
    
    x = [ ( x0 + (xf-x0)*i/N, i=0, N ) ] 
    
    do k=1,  M 
        y1 = Chebyshev( 1, k, x)
        y2 = Chebyshev( 2, k, x); 
        call qplot(x, y1, N+1 ) 
        call qplot(x, y2, N+1 ) 
    end do 
   
end subroutine


!**************************************************************************************************
!
!***line 249 **************************************************************************************
elemental real function Chebyshev( t_kind,  k, x ) 
    integer, intent(in) :: t_kind, k  
    real, intent(in) :: x 
    
       real :: theta 
       
    theta = acos ( x ) 
    
    if (t_kind==1) then 
                                Chebyshev = cos ( k * theta ) 
    else if (t_kind==2) then 
                                Chebyshev = sin ( k * theta ) / sin ( theta )
    end if  
    
    
end function

!**************************************************************************************************
!
!****line 269 *************************************************************************************
elemental real function w_Chebyshev( t_kind, x ) 
    integer, intent(in) :: t_kind  
    real, intent(in) :: x 
    
       
    if (t_kind==1) then 
                                w_Chebyshev = 1 / sqrt( 1 -x**2 ) 
    else if (t_kind==2) then 
                                w_Chebyshev =     sqrt( 1 -x**2 )
    end if  
    
    
end function




!**************************************************************************************************
!
!***line 289 **************************************************************************************
elemental real function f_Chebyshev( x ) 
    real, intent(in) :: x 
    
       
 !   f_Chebyshev = sin ( PI * x ) 
 !    f_Chebyshev = 1 / ( 1 + 25 * x**2 ) 
 !    f_Chebyshev = cos ( PI * x )
     f_Chebyshev = x**20 ! cos ( PI * x ) 
    
end function



!**************************************************************************************************
!
!****line 304 *************************************************************************************
subroutine Chebyshev_expansion( t_kind ) 
    integer, intent(in ) :: t_kind
 
    integer, parameter :: N = 6, M = 500
    real :: x(0:N), f(0:N)
    real :: I_N(0:N, 0:M)            
    real :: xp(0:M), yp(0:M), fp(0:M), Integrand(0:M)  
    
    integer :: i, k
    real :: c_k
    real :: a=-0.999999, b = 0.999999, ymax, gamma
    
    xp = [ (cos( (PI/2 + PI*i)/(M+1)  ), i=M, 0, -1) ]
    fp = f_Chebyshev( xp ) 
  
    yp = 0
    do k=0, N 
       
        Integrand = fp * Chebyshev( t_kind, k, xp ) * w_Chebyshev( t_kind, xp )
        
        if (k==0) then;    gamma = 1 / PI; 
                  else;    gamma = 2 / PI; 
        end if 
        
        c_k = Integral( xp , Integrand )  * gamma
        
        yp = yp + c_k * Chebyshev( t_kind, k, xp)  
                                   ! Truncated Chebyshev series
    end do 
    
    if (t_kind == 1) then 
                         x = [ (cos( (PI/2 + PI*i)/(N+1)  ), i=N, 0, -1) ]
    else if ( t_kind == 2) then 
                         x = [ (cos(PI*i/N), i=N, 0, -1) ] 
    end if 
    
    f = f_Chebyshev( x )
    
    I_N = Interpolant(x, f, N, xp) ! Equal to Discrete Chebyshev series  
   
end subroutine

    !ymax = 1.2 
    !
    !call scrmod("reverse") 
    !call metafl("xwin")
    !call disini
    !call graf(-1., 1., -1., 0.2, -ymax, ymax, -ymax, ymax/10); 
    !
    !  call curve(xp, fp, M+1 )
    !
    !call color("red")
    !call curve(xp, I_N(0, :), M+1 )
    !call color("green")
    !call curve(xp, yp, M+1) 
    !
    !call color("blue")
    !call marker(21)
    !call incmrk(1) 
    !call curve(x, f, N+1 )
    !call disfin 

end module 
    
    
    
    