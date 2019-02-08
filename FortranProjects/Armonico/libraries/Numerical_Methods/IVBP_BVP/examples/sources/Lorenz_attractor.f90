

module Lorenz

 use Cauchy_problem
 use dislin 

 implicit none

 private 

 real :: a=10., b=28., c=2.6666666666, d=50. 
 real :: t0 =0, tf=20;
 integer, save :: N 

 real, save, allocatable :: tv(:), V(:,:)  ! grahs 
 integer, save :: it =0; 

 public :: Test_Lorenz_Attractor

contains 

!******************************************************************************
!*
!******************************************************************************
subroutine Time_Domain( Time ) 
 real, pointer :: Time(:) 
 
    integer :: i


     N = 2000
     allocate ( Time(0:N) ) 

     Time = [ (t0 + (tf-t0)*i/N, i=0, N ) ]  
     

end subroutine 

!******************************************************************
!*
!*****************************************************************
subroutine  Lorenz_IC( t0, U)
   real, intent(in) :: t0 
   real, pointer :: U(:)

   
     allocate( U(3) ) 

	 U(1) = 12 ! sqrt( c *(b-1) ) + 1 
     U(2) = 15 ! U(1) + 1  
	 U(3) = 30 ! b - 1  + 0.5  

     allocate( tv(0:N), V(3, 0:N) ); 

end subroutine 

!*******************************************************************
!*
!*******************************************************************
subroutine  Lorenz_attractor(t, U,  F)            
              real, intent(in)  :: t
              real, intent(inout)  ::  U(:)
              real, intent(out) :: F(:)

      
      real :: x, y , z
	  
	 x = U(1); y = U(2); z = U(3)   

     F(1) =   a * ( y - x ) 
	 F(2) =   x * ( b - z ) - y 
	 F(3) =   x * y - c * z   

end subroutine 

!******************************************************************************
!*
!******************************************************************************
subroutine Lorenz_graphs(t, U, S ) 
          real, intent(in) :: t
          real, intent(inout) :: U(:) 
          procedure (ODES_Outputs), optional :: S 
  
    
  
  real :: x0, xf, y0, yf, z0, zf, dh


 !  write(*,*) U(:) 
   tv(it) = t 
   V(:, it) = U(:)
    
   it = it + 1; 

   if (t==tf) then 

         call scrmod('revers') 
         call qplot( V(1, :), V(2,:), N+1 ); 

   endif 


  
end subroutine  



!*******************************************************************
!*
!*******************************************************************
subroutine Test_Lorenz_Attractor


    call Cauchy_Problem_Solution( Domain = Time_Domain,  Initial_C = Lorenz_IC, System = Lorenz_Attractor, Outputs = Lorenz_graphs ) 


end subroutine
 



end module 
