module Waves2D
 
 use Finite_differences
 use Cauchy_Problem 
 
 use dislin

 implicit none 
 private 

 public :: Test_Wave2D
 real, save :: t0 = 0., tf = 7.5 
 integer, save :: Nx = 20, Ny = 20, Order = 8; 
 real, save, allocatable :: x(:), y(:) 

 
integer :: i_graph = 0 

contains 


!******************************************************************************
!*
!******************************************************************************
subroutine Time_Domain( Time ) 
 real, pointer :: Time(:) 
 
    integer :: i, N 
    
     N = 75000
     allocate ( Time(0:N) ) 

     Time = [ (t0 + (tf-t0)*i/N, i=0, N ) ]  

end subroutine 

!******************************************************************
!*
!*****************************************************************
subroutine  Wave2D_IC( t0, U)
   real, intent(in) :: t0 
   real, pointer :: U(:)

     integer :: M 
     real :: a=0, b=1;  
     
    M = (Nx+1)*(Ny+1); 

    allocate( x(0:Nx),  y(0:Ny), U(2*M) ) 

    call Grid_Initialization( "uniform", "x", Order, x )
    call Grid_Initialization( "nonuniform", "y", Order, y )


    call IC(Nx, Ny, U(1:M), U(M+1:2*M) ) 


contains 

subroutine IC( Nx, Ny, P, Pt  ) 
      integer, intent(in) :: Nx, Ny 
      real, intent(out) ::   P(0:Nx, 0:Ny), Pt(0:Nx, 0:Ny) 

    integer :: i, j 
    real :: x0, y0  

    x0 = x(Nx/2); 
    y0 = y(Ny/2); 

    do i=0, Nx 
       do j=0, Ny 
          P(i,j) =  exp( -25*( (x(i)-x0)**2 + (y(j)-y0)**2 ) ) 
       enddo 
    enddo 

    Pt = 0; 


end subroutine 



end subroutine 

!*******************************************************************
!*
!*******************************************************************
subroutine  Wave2D_equation(t, U,  F)            
              real, intent(in)  :: t
              real, intent(inout)  ::  U(:)
              real, intent(out) :: F(:)

  integer :: M 

  M = (Nx+1)*(Ny+1) 
 
  call Wave2D_equations( Nx, Ny, U(1:M), U(1+M:2*M), F(1:M), F(1+M:2*M) ) 

contains 

subroutine Wave2D_equations( Nx, Ny, P, W, Fp, Fw  ) 
      integer, intent(in) :: Nx, Ny 
      real, intent(inout) ::   P(0:Nx, 0:Ny),  W(0:Nx, 0:Ny) 
      real, intent(out) :: Fp(0:Nx, 0:Ny), Fw(0:Nx, 0:Ny) 

      real :: Pxx(0:Nx, 0:Ny), Pyy(0:Nx, 0:Ny) 

 !** Boundary conditions 
      call Newmann("x", 0, P, 0.0 );   call Newmann("x", Nx, P, 0.0 );
      call Newmann("y", 0, P, 0.0 );   call Newmann("y", Ny, P, 0.0 ); 

  !** Wave equation
      call Derivative( "x", 2, P, Pxx ) 
      call Derivative( "y", 2, P, Pyy ) 

      Fp =  W 
      Fw =  4 * (Pxx + Pyy) 

end subroutine 

end subroutine 

!******************************************************************************
!*
!******************************************************************************
subroutine Wave2D_graphs(t, U, S ) 
          real, intent(in) :: t
          real, intent(inout) :: U(:) 
          procedure (ODES_Outputs), optional :: S 
  
 !   integer, parameter :: NNx = 100, NNY = 100 
  !  real :: W(0:NNx, 0:NNY) 
    integer :: M 

  M = (Nx+1)*(Ny+1) 
    
  i_graph = i_graph + 1 
  
  if (mod(i_graph, 500) == 0) then    
     write(*,*) " t = ", t, maxval(U), minval(U) 
  endif   
  ! read(*,*) 
 
  if ( (t==tf) .or. (t==0) ) then 
     call scrmod('revers') 
  !   call metafl('xwin') 
   !  call QPLCON( U(1:M), Nx+1, Ny+1, 10);
  !   call metafl('PDF') 
     call QPLClr( U(1:M), Nx+1, Ny+1);
   !  call Generate_seeds( Order, Seeds ) 
    ! do i=0, NNx
     !  do j=0, NNy  
      !   jp = minloc( abs(x_nodes-xp(i)) ) 
       !  j = jp(1) - 1
       !  s = Seeds(j) - 1
       !  q = Order 
       !  Error(0:2, i) =  Lagrange_Error_Polynomial( x_nodes(s:s+q), xp(i) )
       !enddo 





  endif    

  
end subroutine  



!*******************************************************************
!*
!*******************************************************************
subroutine Test_Wave2D


    call Cauchy_Problem_Solution( Domain = Time_Domain,  Initial_C = Wave2D_IC, System = Wave2D_equation, Outputs = Wave2D_graphs ) 


end subroutine



end module 

