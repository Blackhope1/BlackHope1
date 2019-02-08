module Waves1D
 
 use Finite_differences
 use Cauchy_Problem 
 

 !use nr
 use dislin

 implicit none 
 private 

 public :: Test_Wave1D
 real, save :: t0 = 0., tf = 0.6 
 integer, save :: Nx = 100, Order = 50; 
 real, save, allocatable :: x(:)

 


contains 


!******************************************************************************
!*
!******************************************************************************
subroutine Time_Domain( Time ) 
 real, pointer :: Time(:) 
 
    integer :: i, N 
    
     N = 6000
     allocate ( Time(0:N) ) 

     Time = [ (t0 + (tf-t0)*i/N, i=0, N ) ]  

end subroutine 

!******************************************************************
!*
!*****************************************************************
subroutine  Wave1D_IC( t0, U)
   real, intent(in) :: t0 
   real, pointer :: U(:)

     integer :: M 
     real :: a=0, b=1;  
     
      M = (Nx+1); 

     allocate( x(0:Nx),  U(2*M) )   
   
  

    call Grid_Initialization( "uniform", "x", Order, x )
   

    call IC(Nx,  U(1:M), U(M+1:2*M) ) 


contains 

subroutine IC( Nx, P, Pt  ) 
      integer, intent(in) :: Nx
      real, intent(out) ::   P(0:Nx), Pt(0:Nx) 

    integer :: i
    real :: a

    a = x(Nx/2); 

    do i=0, Nx 
          P(i) =  exp( -500*(x(i)-a)**2 ) 
    enddo 
    Pt = 0; 


end subroutine 



end subroutine 

!*******************************************************************
!*
!*******************************************************************
subroutine  Wave1D_equation(t, U,  F)            
              real, intent(in)  :: t
              real, intent(inout)  ::  U(:)
              real, intent(out) :: F(:)

  integer :: M 

  M = (Nx+1)
 
  call Wave1D_equations( Nx, U(1:M), U(1+M:2*M), F(1:M), F(1+M:2*M) ) 

contains 

subroutine Wave1D_equations( Nx, P, W, Fp, Fw  ) 
      integer, intent(in) :: Nx 
      real, intent(inout) ::   P(0:Nx),  W(0:Nx) 
      real, intent(out) :: Fp(0:Nx), Fw(0:Nx) 

      real :: Pxx(0:Nx), Pyy(0:Nx) 

 !** Boundary conditions 
    call Newmann( 0, P, 0.0 );   call Newmann( Nx, P, 0.0 );
 !   P(0) = 0;  P(Nx) = 0; 


  !** Wave equation
      call Derivative( "x", 2, P, Pxx )
   
      Fp =  W 
      Fw =  Pxx 

end subroutine 

end subroutine 

!******************************************************************************
!*
!******************************************************************************
subroutine Wave1D_graphs(t, U, S ) 
          real, intent(in) :: t
          real, intent(inout) :: U(:) 
          procedure (ODES_Outputs), optional :: S 
  
    integer :: i, M 
    real :: pi = 4 * atan(1d0); 

    real index(Nx), dx(Nx) 

    M = (Nx+1)

  ! write(*,*) " t = ", t, maxval(U), minval(U) 
  ! read(*,*) 
 
  if ( (t==tf) .or. (t==0) ) then 
     call scrmod('revers') 
     call metafl('xwin') 
     call Qplot( x, U, Nx+1);
  endif    

  
end subroutine  



!*******************************************************************
!*
!*******************************************************************
subroutine Test_Wave1D


   call Cauchy_Problem_Solution( Domain = Time_Domain,  Initial_C = Wave1D_IC, System = Wave1D_equation, Outputs = Wave1D_graphs ) 


end subroutine






end module 

