module NSB_Convection_cell

    use Finite_differences
    use Cauchy_Problem
    use Non_Linear_Systems
    implicit none

contains



!***********************************************************************
!*
!***********************************************************************
subroutine Convection_cell2D

       integer, parameter :: Nx = 20, Ny = 20, Nt = 5, Nv = 2
       real ::  x(0:Nx), y(0:Ny), Time(0:Nt), Us(0:Nt, Nv*(Nx+1)*(Ny+1) )
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
       real ::  t0 = 0, tf = 0.01
       integer ::  i, M, Order = 8 
    
           
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*i/Nx, i=0, Ny ) ]
     
     call Grid_Initialization( "nonuniform", "x", Order, x )
     call Grid_Initialization( "nonuniform", "y", Order, y )
    
     M = (Nx+1)*(Ny+1) 
     call  Vorticity_IC( Us(0,1:M), Us(0,M+1:2*M) )
     
     
     call scrmod('revers')   
     call metafl('xwin')
     call disini() 
     CALL TITLIN (" Initial condition: Temperature ", 2)
     call QPLClr( Us(0, M+1:2*M), Nx+1, Ny+1);
     
     call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = Navier_StokesB, Solution = Us )
     
     call disini() 
     CALL TITLIN (" Convection cell, Temperature ", 2)
     call QPLClr( Us(Nt,M+1:2*M), Nx+1, Ny+1);

    
contains
!------------------------------------------------------------------------
subroutine Vorticity_IC(  W, T )
   real :: W(0:Nx, 0:Ny), T(0:Nx, 0:Ny)

   integer :: i, j

    W = 0 
    do i=0, Nx 
        T(i,:) = ( 1 + x(i)**3 ) / 2
    enddo

      

end subroutine
!-----------------------------------------------------------------------
function  Navier_StokesB( V, t ) result(F) 
          
              real :: V(:), t 
              real :: F(size(V))


  call Navier_StokesB_equations( Nx, Ny, V(1:M), V(M+1:2*M), F(1:M), F(M+1:2*M) ) 
           
   
end function  
end subroutine 

!-----------------------------------------------------------------------
subroutine Navier_StokesB_equations( Nx, Ny, W, T, Wt, Tt )
        integer, intent(in) :: Nx, Ny 
        real :: W(0:Nx, 0:Ny), T(0:Nx, 0:Ny), Wt(0:Nx, 0:Ny), Tt(0:Nx, 0:Ny)   
           
        
       real :: Wx(0:Nx, 0:Ny), Wy(0:Nx, 0:Ny), Wxx(0:Nx, 0:Ny), Wyy(0:Nx, 0:Ny)
       real :: Tx(0:Nx, 0:Ny), Ty(0:Nx, 0:Ny), Txx(0:Nx, 0:Ny), Tyy(0:Nx, 0:Ny)
       real :: Psi(0:Nx, 0:Ny), Psixx(0:Nx, 0:Ny), Psiyy(0:Nx, 0:Ny)  

       real :: u(0:Nx,0:Ny), v(0:Nx,0:Ny)   ! velocities  
       integer :: i, j 
       real :: Re = 1, Pr = 1 , Ra  = 1 
       
       real :: BC_neu(0:Nx)
   
   
      call Stream_function( Nx, Ny, Psi, W ) 
      
     call disini() 
     CALL TITLIN (" Stream function ", 2)
     call qplcon( Psi, Nx+1, Ny+1, 20);
     
     call disini() 
     CALL TITLIN (" Temperature ", 2)
     call QPLClr( T, Nx+1, Ny+1);

      call  Derivative( ["x","y"], 2, 1, Psi, u )
      call  Derivative( ["x","y"], 1, 1, Psi, v ); v = - v 

      call  Derivative( ["x","y"], 2, 2, Psi, Psiyy )  
      call  Derivative( ["x","y"], 1, 2, Psi, Psixx )
    

! *** boundary_conditions
      W(0, :)  = - Psixx(0,   :) - Psiyy(0, :)
      W(Nx, :) = - Psixx(Nx,  :) - Psiyy(Nx, :)
      W(:, 0)  = - Psixx(:,  0)  - Psiyy(:, 0)
      W(:, Ny) = - Psixx(:,  Ny) - Psiyy(:, Ny)

      BC_neu(:) = 0.
  !    call Neumann( ["x","y"], 2, 0, T, BC_neu ); call Neumann( ["x","y"], 2, Ny, T, BC_neu )
      T(0,:)= 0. ; T(Nx,:)= 1
      
      call  Derivative( ["x","y"], 2, 1, W, Wy );  call  Derivative( ["x","y"], 1, 1, W, Wx )
      call  Derivative( ["x","y"], 2, 2, W, Wyy ); call  Derivative( ["x","y"], 1, 2, W, Wxx )
      call  Derivative( ["x","y"], 2, 1, T, Ty );  call  Derivative( ["x","y"], 1, 1, T, Tx )
      call  Derivative( ["x","y"], 2, 2, T, Tyy ); call  Derivative( ["x","y"], 1, 2, T, Txx )
  

! *** inner grid points
      Wt = -u * Wx - v * Wy + (1/Re) * ( Wxx + Wyy ) - (Ra/(Re*Re*Pr)) * Tx
     
      Tt = -u * Tx - v * Ty + (1/(Re*Pr)) * ( Txx + Tyy )
     
end subroutine


!********************************************************************************************
!*
!********************************************************************************************
subroutine Stream_function( Nx, Ny, Psi, W ) 
          integer, intent(in) :: Nx, Ny 
          real, intent(out) :: Psi(0:Nx, 0:Ny)
          real, intent(in)  ::   W(0:Nx, 0:Ny) 
     
     
    integer :: i,j, ij, M 
    real, allocatable, save :: Difference_operator(:,:), b(:), V(:)  
    real :: F(0:Nx, 0:Ny) 
                    
        
       M = (Nx+1)*(Ny+1)
       allocate( Difference_operator(M, M), b(M), V(M)  )
     
!  *** Delta kronecker to calculate the difference operator
       do i=0, Nx
           do j=0, Ny
              Psi = 0
              Psi(i,j) = 1.0
              
              F =  Laplacian_equation( Nx, Ny, Psi)
              
              ij = i + j * (Nx+1) + 1
              Difference_operator(1:M,ij) = reshape( F, [ M ])

          enddo
       enddo
     

!  *** solve the linear system of equations
       call LU_Factorization(Difference_operator)

       b = -reshape( W, [ M ])

       V = Solve_LU( Difference_operator, b )

       Psi = reshape( V, [Nx+1, Ny+1] )


      deallocate( Difference_operator, b, V  )

end subroutine 

!-----------------------------------------------------------------
!
!-----------------------------------------------------------------
function Laplacian_equation(Nx, Ny, Psi) result (F) 
     integer, intent(in) :: Nx, Ny 
     real, intent(in) :: Psi(0:Nx, 0:Ny) 
     real :: F(0:Nx, 0:Ny) 
         
    real :: Psixx(0:Nx, 0:Ny),  Psix(0:Nx, 0:Ny) 
    real :: Psiyy(0:Nx, 0:Ny),  Psiy(0:Nx, 0:Ny) 

    integer ::  i, j 

        call Derivative( ["x","y"], 1, 1, Psi, Psix )
        call Derivative( ["x","y"], 2, 1, Psi, Psiy )
       
        call Derivative( ["x","y"], 1, 2, Psi, Psixx )
        call Derivative( ["x","y"], 2, 2, Psi, Psiyy )

!   *** inner grid points
        F = Psixx + Psiyy
        
!  ***  boundary conditions
        do j = 0, Ny
           F(0,  j)   = Psi (0, j) 
           F(1,  j)   = Psix(0, j) 
           F(Nx, j)   = Psi (Nx, j) 
           F(Nx-1, j) = Psix(Nx, j) 
        enddo


        do i=1, Nx-1
           F(i,  0)   = Psi (i, 0) 
           F(i,  1)   = Psiy(i, 0) 
           F(i,  Ny)  = Psi (i, Ny) 
           F(i, Ny-1) = Psiy(i, Ny)
        end do 



end function 



end module


