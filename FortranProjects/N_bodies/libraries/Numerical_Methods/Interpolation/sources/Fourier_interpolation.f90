module Fourier_interpolation

 implicit none 
 real :: PI = 4 * atan(1d0)
 complex   :: II = (0, 1)
 
 private 
 public :: FFT2
 
    contains 
    
    
    
    
!*************************************************************************
! Cooley-Tukey FFT
!            Fast Fourier Transform
!
!    Inputs:
!              - N        : number of collocation points(must be N=2^r)
!              - U(0:N-1) : collocation points 
!                           u_0, u_1,.. u_N-1
!
!    Outputs:
!              - U(0:N-1) : discrete fourier coefficients 
!                           N *(c_0, c_1, ... c_N/2, c_N/2-1, ..., c_-1)
!
!
!  Author: Juan A. Hernandez, Feb, 2018    
!*************************************************************************
  recursive subroutine FFT2(N, U)
     integer, intent(in) :: N 
     complex, intent(inout)  :: U(0:N-1)
    
    complex   :: w, Even(0:N/2-1), Odd( 0:N/2-1 )
    integer   :: k 
 
    if(N <= 1) return
    
    Even = U(0:N-1:2)
    Odd  = U(1:N-1:2)
    
    call FFT2( N/2, Even)
    call FFT2( N/2, Odd)
   
    do k = 0, N/2-1
        
       w = exp( -2*PI*II * k/N )
       U(k)     = Even(k) + w * Odd(k)
       U(k+N/2) = Even(k) - w * Odd(k)
       
    end do
 
  end subroutine  



end module 
