module API_Example_Chebyshev_Interpolation

    use Interpolation
    
    implicit none
    real :: PI = 4 * atan(1d0)
    
    
    contains
    
    
subroutine Chebyshev_Interpolation_examples

    call test_dct2
    call test_fft2
   
   
end subroutine  



!*****************************************
subroutine  test_dct2
!*****************************************

  integer, parameter :: N = 64
  complex :: U(0:N) 
  

  integer :: i
  real :: k(0:N)
  complex ::  ck(0:N), C1(0:N)   
  real :: x(0:N), y(0:N)  
  real :: a= -1, b= 1
  
  
  x = [ ( cos(PI*i/N), i=N, 0, -1) ]  
  U = 1  + 2 * x**2 -1 ! sin( PI * X ) ! 2*x**2 !1/(1+25*x**2) !   + 4* x**3 - 3*x + 8*x**4 - 8*x**2 + 1 !  + 5* x**3
 ! u= [ 1, 0, 0, 1, 2 ] 
!  write(*,*) U 
  
  call scrmod("reverse")
  call qplot( x, real(U), N+1 )
 
  call DCT2(N, U )
 ! write(*,*) " u = ", U 
  
   ck(0) =  U(0) / N 
   ck(1:N) = 2 * U(1:N)  / N
   call qplot([(real(i),i=0,N)], real(ck), N+1)
  
   U = 2 * U / N 
  
   call DCT2(N, U)  
 
   call qplot( x, real(U), N+1 )
 !  call qplot( x, imag(U), N+1 )
 
   
  call Chebyshev_Derivative(Ck, C1)
  
  call DCT2(N, C1) 
  call qplot( x, real(C1), N+1 )
   
  
end subroutine      
    
    
  

!*****************************************
subroutine  test_fft2
!*****************************************

  integer, parameter :: N = 256 
  complex :: u(0:N-1) 
  

  integer :: i
  real :: k(0:N/2-1), ck(0:N/2-1)  
  real :: x(0:N-1), y(0:N-1)  
  real :: a, b 
  
  x = [ ( 2*PI*i/N, i=0, N-1 ) ]  
  u =  cos( x)  + sin (4 * x ) 
  y = u  
  
  call scrmod("reverse")
  call qplot( x, real(u), N ) 
 
  call FFT2(N, u )
  
   do i=1, N/2-1 
    k(i) = i 
    ck(i) = 2 * abs(u(i)) / N 
   end do 
   call qplot( k, ck, N/2 ) 
  
  
  u = conjg( u ) / N    ! (Duhamel et al., 1988)
 
  
  
  call FFT2(N, u)  
  u = conjg( u) 
   call qplot( x, real(u), N )
   call qplot( x, imag(u), N )
 
  
  
 
   
  
end subroutine  



end module 
    
    
    
    