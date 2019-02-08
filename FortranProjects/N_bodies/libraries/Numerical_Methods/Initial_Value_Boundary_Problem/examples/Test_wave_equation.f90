
!*******************************************************************
!*
!*******************************************************************
module Test_Wave_equation 



use Finite_differences
use Cauchy_Problem  
implicit none 




contains 

!*******************************************************************
!*
!*******************************************************************
subroutine Test_Wave1D

       integer, parameter :: Nx = 100, Nt = 200, Nv = 2 
       real ::  x(0:Nx), Time(0:Nt), U(0:Nt, 2*(Nx+1))
       
       real ::  x0 = -1, xf = 1
       real ::  t0 = 0, tf = 2.5 
       integer ::  i, M, Order = 8 
       real :: a 
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     a = x(Nx/2);  
     M = Nx+1 
     
     call Grid_Initialization( "nonuniform", "x", Order, x )
     
     U = 0 
     U(0, 1:Nx+1) =   exp( -100*(x-a)**2 )       
    
     
     call scrmod('revers')   
     call metafl('xwin')
     call disini() 
     CALL TITLIN (" Initial condition ", 2)
     call qplot(x, U(0, :), Nx+1) 
   
     call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = Wave1D, Solution = U )
     
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $  Utt = Uxx  $ ", 2)
     CALL TITLIN (" BC: Neumann   ", 4)
     call qplot(x, U(Nt, :), Nx+1)                                 

   
contains 
!-----------------------------------------------------------------------
function  Wave1D( V, t ) result(F)            
              real  ::  V(:), t 
              real  :: F(size(V))
    
    call Wave1D_equations( V(1:M), V(1+M:2*M), F(1:M), F(1+M:2*M) ) 
    
 
end function
!-----------------------------------------------------------------------
subroutine Wave1D_equations(  P, W, Fp, Fw  ) 
      real, intent(inout) ::   P(0:Nx),  W(0:Nx) 
      real, intent(out) :: Fp(0:Nx), Fw(0:Nx) 

      real :: Pxx(0:Nx)
   
 !** Boundary conditions 
      call Neumann("x", 0, P, 0. );   call Neumann("x", Nx, P, 0. );

  !** Wave equation
      call Derivative( "x", 2, P, Pxx ) 

      Fp =  W 
      Fw =  Pxx  

end subroutine 

end subroutine 



!*******************************************************************
!*
!*******************************************************************
subroutine Test_Wave2D

       integer, parameter :: Nx = 50, Ny = 50, Nt = 200, Nv = 2 
       real ::  x(0:Nx), y(0:Ny), Time(0:Nt), U(0:Nt, 2*(Nx+1)*(Ny+1) )
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
       real ::  t0 = 0, tf = 2.5 
       integer ::  i, M, Order = 8 
           
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*i/Nx, i=0, Ny ) ]
     
     call Grid_Initialization( "nonuniform", "x", Order, x )
     call Grid_Initialization( "nonuniform", "y", Order, y )
    
     M = (Nx+1)*(Ny+1) 
     call  Wave2D_IC( U(0,1:M), U(0,M+1:2*M) )
     
     
     call scrmod('revers')   
     call metafl('xwin')
     call disini() 
     CALL TITLIN (" Initial condition ", 2)
     call QPLClr( U(0, 1:M), Nx+1, Ny+1);
     
   
     call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = Wave2D, Solution = U )
     
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $  Utt = Uxx + Uyy $ ", 2)
     CALL TITLIN (" BC: Neumann   ", 4)
     call QPLClr( U(Nt,1:M), Nx+1, Ny+1);

contains  
!-----------------------------------------------------------------------
subroutine Wave2D_IC( P, Pt  ) 
      real, intent(out) ::   P(0:Nx, 0:Ny), Pt(0:Nx, 0:Ny) 

    integer :: i, j 
    
    do i=0, Nx 
       do j=0, Ny 
          P(i,j) =  exp( -25*( x(i)**2 + y(j)**2 ) ) 
       enddo 
    enddo 

    Pt = 0; 

end subroutine 
!-----------------------------------------------------------------------
function  Wave2D( U, t ) result(F)            
              real  :: U(:), t
              real :: F( size(U) )
  
  call Wave2D_equations( U(1:M), U(1+M:2*M), F(1:M), F(1+M:2*M) ) 


end function
!----------------------------------------------------------------------
subroutine Wave2D_equations(  P, W, Fp, Fw  ) 
      real, intent(inout) ::   P(0:Nx, 0:Ny),  W(0:Nx, 0:Ny) 
      real, intent(out) :: Fp(0:Nx, 0:Ny), Fw(0:Nx, 0:Ny) 

      real :: Pxx(0:Nx, 0:Ny), Pyy(0:Nx, 0:Ny) 
      
      real :: BC_neu(0:Nx)

 !** Boundary conditions 
      BC_neu(:) = 0.
      call Neumann(["x","y"], 1, 0, P, BC_neu );   call Neumann(["x","y"], 1, Nx, P, BC_neu );
      call Neumann(["x","y"], 2, 0, P, BC_neu );   call Neumann(["x","y"], 2, Ny, P, BC_neu ); 

  !** Wave equation
      call Derivative( ["x","y"], 1, 2, P, Pxx ) 
      call Derivative( ["x","y"], 2, 2, P, Pyy ) 

      Fp =  W 
      Fw =  4 * (Pxx + Pyy) 

end subroutine 

end subroutine 


end module 


