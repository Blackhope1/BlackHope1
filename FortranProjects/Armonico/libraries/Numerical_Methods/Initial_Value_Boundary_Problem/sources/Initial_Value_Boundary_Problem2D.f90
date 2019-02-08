
!***********************************************************************
! It integrates in time the Initial value boundary problem 2D.
! Given the differential operator and the boundary equations,
! the discrete solution is calculated.
! author: Juan A Hernandez, juanantonio.hernandez@upm.es 
!***********************************************************************
module Initial_Value_Boundary_Problem2D

use Cauchy_Problem
use Finite_differences
use Non_Linear_Systems
use Utilities
use Dependencies

implicit none   

private

public :: IVBP2D
public :: IVBP2D_system

   
abstract interface  

             
       real function DifferentialOperator2D(x, y, t, u, ux, uy, uxx, uyy, uxy) 
                        real, intent(in) :: x, y, t, u, ux, uy, uxx, uyy, uxy
       end function  

       real function      BC2D(x, y, t, u, ux, uy) 
           real, intent(in) :: x, y, t, u, ux, uy 
       end function 
       
      function DifferentialOperator2D_system(x, y, t, u, ux, uy, uxx, uyy, uxy) 
                        real, intent(in) :: x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
                        real :: DifferentialOperator2D_system (size(u)) 
       end function  

       function      BC2D_system(x, y, t, u, ux, uy) 
           real, intent(in) :: x, y, t, u(:), ux(:), uy(:) 
           real :: BC2D_system(size(u)) ! maximum number of boundary conditions at each point
       end function  


       function DifferentialOperatorODE(U, t) result(F) 
                        real, intent(in) :: U(:), t 
                        real :: F(size(U)) 
       end function 

end interface


 type Boundary2D 
     
        integer :: i, j 
        real :: value 
        real :: equation 
        
 end type
 
  type Boundary2D_system 
     
        integer :: i, j 
        real, allocatable :: value(:) 
        real, allocatable :: equation(:)  
        
end type
 

contains



!*******************************************************************
!*
!*******************************************************************
subroutine IVBP2D( Time_Domain, x_nodes, y_nodes, Order, Differential_operator,  Boundary_conditions, Scheme, Solution) 

     real, intent(in) :: Time_Domain(:) 
     real, intent(inout) :: x_nodes(0:), y_nodes(0:)
     integer, intent(in) :: Order
     procedure (DifferentialOperator2D) :: Differential_operator 
     procedure (BC2D) ::  Boundary_conditions
     procedure (Temporal_Scheme), optional :: Scheme
     real, intent(out) :: Solution(0:, 0:, 0:) 

     integer :: Nt, Nx, Ny, M 
     integer :: Nt_u, Nx_u, Ny_u
     real, pointer :: pSolution(:, :) 
     real :: t_BC 
     integer ::  it  
     integer :: N1, N2, N3, N4
     type (Boundary2D), allocatable :: B(:)
     integer :: Nb 
     logical :: dU(5) ! matrix of dependencies(variables, derivative)  
     
     Nt = size(Time_Domain) -1 
     Nx = size(x_nodes) - 1
     Ny = size(y_nodes) - 1
     
     Nx_u = size( Solution(0, :, 0) ) - 1 
     Ny_u = size( Solution(0, 0, :) ) - 1 
     Nt_u = size( Solution(:, 0, 0) ) - 1 
     
     if (Nx /= Nx_u .or. Ny /= Ny_u .or. Nt /= Nt_u ) then 
         
             write(*,*) " ERROR : Non conformal dimensions Time_Domain, x_nodes, y_nodes, Solution"
             write(*,*) " Nx = ", Nx,  " Ny = ", Ny,  " Nt = ", Nt
             write(*,*) " Nx_u = ", Nx_u,  " Ny_u = ", Ny_u,  " Nt_u = ", Nt_u
             stop 
             
     end if         
     
     N1 = Ny + 1; N2 = Ny + 1; N3 = Nx -1; N4 = Nx - 1;  
     M = (Nx+1) * (Ny+1) 
   
     dU = Dependencies_IVBP_2D( Differential_operator ) 
        
     call Data_pointer( N1 = Nt+1, N2 =  (Nx+1) * (Ny+1),  U = Solution(0:Nt, 0:Nx, 0:Ny) , pU = pSolution ) 
    
     call Cauchy_ProblemS(  Time_Domain = Time_Domain, Differential_operator = Space_discretization, Scheme = Scheme, Solution = pSolution )
      
   
contains
!----------------------------------------------------------------------
function  Space_discretization( U, t ) result(F) 
          real ::  U(:), t         
          real :: F(size(U))    

         call Space_discretization2D( U, t, F) 
         
end function 
!---------------------------------------------------------------------------------------------
subroutine Space_discretization2D( U, t, F )
                           real :: U(0:Nx, 0:Ny), t 
                           real :: F(0:Nx, 0:Ny)  
                 

    integer :: i, j, k, Nb
    real :: Ux(0:Nx, 0:Ny), Uxx(0:Nx, 0:Ny), Uxy(0:Nx, 0:Ny), Uy(0:Nx, 0:Ny), Uyy(0:Nx, 0:Ny)
   
    real, allocatable :: Ub(:)
    logical :: l1, l2, NO_BC 
    
        t_BC = t 
        call Binary_search(t_BC, Time_Domain, it)
        Solution(it, :, :)  = U(:, :)

        call Boundary_unknowns( Ub ) 
       
        call Newton( BCs, Ub )
        
        if ( dU(1) ) call Derivative( [ "x", "y" ], 1, 1, U, Ux  )
        if ( dU(2) ) call Derivative( [ "x", "y" ], 2, 1, U, Uy  )
        if ( dU(3) ) call Derivative( [ "x", "y" ], 1, 2, U, Uxx )
        if ( dU(4) ) call Derivative( [ "x", "y" ], 2, 2, U, Uyy )
        if ( dU(5) ) call Derivative( [ "x", "y" ], 2, 1, Ux, Uxy)

!  *** inner grid points
        F = 0
        do i=0, Nx
             do j=0, Ny
                 
                 k = boundary_index(i,j) 
                 if (k>0) then 
                    l1 = B(k) % equation == FREE_BOUNDARY_CONDITION
                    l2 = B(k) % equation == PERIODIC_BOUNDARY_CONDITION
                    NO_BC = l1 .or. l2 
                 end if 
                 
                 if ( k<0 .or. NO_BC) then 
                       F(i,j) = Differential_operator( x_nodes(i), y_nodes(j), t, U(i,j), Ux(i,j), Uy(i,j), Uxx(i,j), Uyy(i,j), Uxy(i,j) )
                 end if   
             enddo
        enddo
        
        deallocate( Ub, B ) 


end subroutine 


!-----------------------------------------------------------------------
subroutine Boundary_unknowns( Ub ) 
real, allocatable, intent(out) :: Ub(:)
   
       integer :: i, j, k
       real :: u0 = 1., ux0= 2., uy0 = 3. 
  
  allocate( B( N1 + N2 + N3 + N4 ) )    
             
  do k=1, size(B)
      
    call  ij_index( k, i, j) 
    
    B(k) % equation   =   Boundary_conditions (x_nodes(i),  y_nodes(j),  t_BC, u0,  ux0,  uy0 ) 
    B(k) % i = i 
    B(k) % j = j 
    B(k) % value = Solution(it, i, j) 
   
  end do
  
  Nb = count( B % equation /= FREE_BOUNDARY_CONDITION .and. B % equation /= PERIODIC_BOUNDARY_CONDITION ) 
  
  allocate ( Ub(Nb) )
  Ub = pack( B(:) % value, B % equation /= FREE_BOUNDARY_CONDITION .and.  B % equation /= PERIODIC_BOUNDARY_CONDITION ) 
  
end subroutine 

!-----------------------------------------------------------------------
function BCs(Y) result(G) 
    real, intent(in) :: Y(:)
    real :: G(size(Y))

 integer :: i, j, k, m  
 real :: Wx(0:Nx,0:Ny), Wy(0:Nx,0:Ny)
  
  m = 1           
  do k=1, size(B) 
      
      if( B(k)% equation /= FREE_BOUNDARY_CONDITION) then 
            i = B(k)% i
            j = B(k)% j
            Solution(it,i,j) = Y(m)
            m = m + 1 
      end if 
      
  end do 
            
  call Derivative( [ "x", "y" ], 1, 1, Solution(it,:,:), Wx(:,:))
  call Derivative( [ "x", "y" ], 2, 1, Solution(it,:,:), Wy(:,:))

  m = 1 
  do k=1, size(B) 
       
       if (B(k)% equation /= FREE_BOUNDARY_CONDITION) then 
           
           i = B(k) % i 
           j = B(k) % j 
           G(m)   = Boundary_conditions (x_nodes(i),  y_nodes(j),  t_BC, Solution(it, i, j),  Wx(i, j),  Wy(i, j) )  
           m = m + 1 
           
       end if 
       
  end do

end function 

!-----------------------------------------------------------------------
subroutine ij_index(k, i, j) 
   integer, intent(in) :: k
   integer, intent(out) :: i, j 
   
    if (k<=N1) then 
                       i = 0; j = k-1 
                       
    else if (k<=N1+N2) then 
        
                       i = Nx; j = k - N1 - 1
                       
    else if (k<=N1+N2+N3) then 
        
                       j = 0; i = k - N1 - N2 - 1
    else 
                       j = Ny; i = k - N1 - N2 - N3 - 1
    endif 

end subroutine 

!------------------------------------------------------------------------
integer function boundary_index(i,j) result(k) 
         integer, intent(in) :: i, j 

    if (i==0) then 
        
               k = j + 1  
               
    else if (i==Nx) then 
        
               k = 1 + N1 + j 
               
    else if (j==0) then 
        
               k = 1 + N1 + N2 + i 
               
    else if (j==Ny) then  
        
               k = 1 + N1 + N2 + N3 + i 
    else 
               k = -1 
    end if 
    

end function 


end subroutine




























!********************************************************************************************************************************
!*
!********************************************************************************************************************************
subroutine IVBP2D_system( Time_Domain, x_nodes, y_nodes, Order, N_variables, Differential_operator,  Boundary_conditions, Scheme, Solution) 

     real, intent(in) :: Time_Domain(:) 
     real, intent(in) :: x_nodes(0:), y_nodes(0:)
     integer, intent(in) :: Order, N_variables 
     procedure (DifferentialOperator2D_system) :: Differential_operator 
     procedure (BC2D_system) ::  Boundary_conditions
     procedure (Temporal_Scheme), optional :: Scheme
     real, intent(out) :: Solution(0:, 0:, 0:, :) 
   
     
     real, pointer :: pSolution(:, :) 
     real :: t_BC 
     integer ::  it  
     type (Boundary2D_system), allocatable :: B(:)
     integer :: Nb 
     logical :: dU(N_variables, 5) ! matrix of dependencies(variables, derivative)
     
     integer :: Nx, Ny, Nt, Nv, Nt_u, Nx_u, Ny_u, Nv_u
     
    
     Nx = size(x_nodes) - 1
     Ny = size(y_nodes) - 1
     Nt = size(Time_Domain) - 1 
     Nv = N_variables 
     
     Nt_u = size( Solution(:, 0, 0, 1) ) - 1
     Nx_u = size( Solution(0, :, 0, 1) ) - 1 
     Ny_u = size( Solution(0, 0, :, 1) ) - 1 
     Nv_u = size( Solution(0, 0, 0, :) ) 
     
     if (Nx /= Nx_u .or. Ny /= Ny_u .or. Nt /= Nt_u .or. Nv /= Nv_u ) then 
         
             write(*,*) " ERROR : Non conformal dimensions Time_Domain, x_nodes, y_nodes, Solution"
             write(*,*) " Nx = ", Nx," Ny = ", Ny, " Nt = ", Nt, " Nv = ", Nv 
             write(*,*) " Nx_u = ", Nx_u, " Ny_u = ", Ny_u, " Nt_u = ", Nt_u, " Nv_u =", Nv_u 
             stop 
     end if  
     
     dU = .true. 
 
    call Data_pointer( N1 = Nt+1, N2 = Nv*  (Nx+1) * (Ny+1),  U = Solution(0:Nt, 0:Nx, 0:Ny, 1:Nv) , pU = pSolution ) 
    
    call Cauchy_ProblemS(  Time_Domain = Time_Domain, Differential_operator = Space_discretization, Scheme = Scheme, Solution = pSolution )
  
    
contains

!----------------------------------------------------------------------
function  Space_discretization( U, t ) result(F) 
          real ::  U(:), t         
          real :: F(size(U))   

         call Space_discretization_2D_system( U, t, F ) 
         
         
end function    
!-----------------------------------------------------------------------
subroutine  Space_discretization_2D_system( U, t, F )
          real :: U(0:Nx,0:Ny, Nv), t 
          real :: F(0:Nx,0:Ny, Nv)
              
    integer :: i, j, k, iv
    real :: Ux(0:Nx,0:Ny, Nv), Uxx(0:Nx,0:Ny, Nv), Uxy(0:Nx,0:Ny, Nv)
    real :: Uy(0:Nx,0:Ny, Nv), Uyy(0:Nx,0:Ny, Nv)
    
    real, allocatable ::  Ub(:)    
    real :: Fv(Nv), eq 
    logical :: l1, l2, Boundary 
 
       
        t_BC = t 
        call Binary_search(t_BC, Time_Domain, it)
        Solution(it, :, :, :)  = U(:, :, :)
        
!  ***  It identifies the boundary value unknowns       
        call Boundary_unknowns( Ub ) 
   
 !  *** Solve boundary equations  
        call Newton( BCs, Ub )
  
!  *** inner grid points
        do k=1, Nv 
                    
           if ( dU(k,1) ) call Derivative( [ "x", "y" ], 1, 1, U(:,:,k), Ux(:,:,k)  )
           if ( dU(k,2) ) call Derivative( [ "x", "y" ], 2, 1, U(:,:,k), Uy(:,:,k)  )
           if ( dU(k,3) ) call Derivative( [ "x", "y" ], 1, 2, U(:,:,k), Uxx(:,:,k) )
           if ( dU(k,4) ) call Derivative( [ "x", "y" ], 2, 2, U(:,:,k), Uyy(:,:,k) )
           if ( dU(k,5) ) call Derivative( [ "x", "y" ], 2, 1, Ux(:,:,k),Uxy(:,:,k) )
           
        end do 
            
        F = 0
        do i=0, Nx
             do j=0, Ny
                 Fv = Differential_operator( x_nodes(i), y_nodes(j), t, U(i,j,:), Ux(i,j,:), Uy(i,j,:), Uxx(i,j,:), Uyy(i,j,:), Uxy(i,j,:) ) 
                 k = boundary_index(i,j)  
                 
            ! ** Boundary point when k>0
                 if (k>0) then  
                    do iv=1, Nv
                      eq =  B(k) % equation(iv)
                      if ( eq == FREE_BOUNDARY_CONDITION ) then
                              F(i, j, iv) = Fv(iv) 
                      end if   
                    end do
                    
            ! ** inner point    
                 else 
                              F(i, j, :) = Fv(:)
                 end if   
                 
             enddo
        enddo
        
        deallocate( Ub, B ) 
        
end subroutine


!-----------------------------------------------------------------------
subroutine Boundary_unknowns( Ub ) 
real, allocatable, intent(out) :: Ub(:) 
   
       integer :: i, j, k, l, m 
       real :: u0(Nv), ux0(Nv), uy0(Nv) 
       real :: eq 
  
     
  allocate( B( 2*Nx+2 + 2*(Ny-1) ) )
  do k=1, size(B) 
         allocate( B(k) % value(Nv), B(k) % equation(Nv) ) 
  end do 
  
  u0 = 1.; ux0 = 2.; uy0 = 3  
  Nb = 0 
  do k=1, size(B)
      
    call  ij_index( k, i, j) 
    
    B(k) % equation(:)   =   Boundary_conditions (x_nodes(i),  y_nodes(j),  t_BC, u0,  ux0,  uy0 ) 
    B(k) % i = i 
    B(k) % j = j 
    B(k) % value(:) = Solution(it, i, j, :)
    Nb = Nb + count( B(k)% equation /= FREE_BOUNDARY_CONDITION  .and. B(k) % equation /= PERIODIC_BOUNDARY_CONDITION )
  end do
 
  
  allocate ( Ub(Nb) ) 
  m = 1 
  do k=1, size(B)
     do l=1, Nv 
         eq = B(k)% equation(l)
         if ( eq == PERIODIC_BOUNDARY_CONDITION ) then
             i = B(k) % i
             j = B(k) % j
             if (i==0)     Solution(it, 0, :, l) = Solution(it, Nx, :, l)
             if (j==0)     Solution(it, :, 0, l) = Solution(it, :, Ny, l)
             
         else if (eq /= FREE_BOUNDARY_CONDITION ) then 
             Ub(m) = B(k)% value(l)
             m = m + 1 
             
         end if 
     end do  
  end do
  
  
end subroutine 

!-----------------------------------------------------------------------
!subroutine ij_index(k, i, j) 
!   integer, intent(in) :: k
!   integer, intent(out) :: i, j 
!   
!   
!    if (k<=Nx+1) then 
!                       i = 0; j = k-1 
!                       
!    else if (k<=2*nx+2) then 
!        
!                       i = Nx; j = k - Nx - 2
!                       
!    else if (k<=2*Nx+Ny+1) then 
!        
!                       j = 0; i = k - 2*Nx - 3 + 1 
!    else 
!                       j = Ny; i = k - 2 *Nx - 2 * Ny - 2 + 1 
!    endif 
!
!end subroutine 
!-----------------------------------------------------------------------
subroutine ij_index(k, i, j) 
   integer, intent(in) :: k
   integer, intent(out) :: i, j 
   
   
    if (k<=Nx+1) then 
                       i = 0; j = k-1 
                       
                       !write(*,*) "i = ", i, "j = ", j, "k =", k
                       
    else if (k<=2*Nx+2) then ! .and.(k>=Nx+2)) then 
        
                       i = Nx; j = k - Nx - 2
                       !write(*,*) "i = ", i, "j = ", j, "k =", k
    else if (k<=2*Nx+Ny+1) then !.and. (k>=2*Nx+3-1) )then 
        
                       j = 0; i = k - 2*Nx - 3 + 1 
                       !write(*,*) "i = ", i, "j = ", j, "k =", k
    else 
                       j = Ny; i = -(k - 2 *Nx - 2 * Ny - 2 + 1) 
                       !write(*,*) "i = ", i, "j = ", j, "k =", k
    endif 

end subroutine 

!------------------------------------------------------------------------
integer function boundary_index(i,j) result(k) 
         integer, intent(in) :: i, j 

    if (i==0) then 
        
               k = j + 1  
               
    else if (i==Nx) then 
        
               k = Nx+2 + j 
               
    else if (j==0) then 
        
               k =   2*Nx + 2 + i
               
    else if (j==Ny) then  
        
               k =  2*Nx + 2 +  Ny-1 + i 
    else 
               k = -1 
    end if 
    
 end function 


!-----------------------------------------------------------------------
function BCs(Y) result(G) 
        real, intent (in) :: Y(:)
        real :: G(size(Y))
    
       real :: Wx(0:Nx, 0:Ny, Nv), Wy(0:Nx, 0:Ny, Nv)
       integer :: i, j, k, m, iv  
       real :: Gv(Nv), eq  
      
  m = 1             
  do k=1, size(B) 
     do iv=1, Nv  
                 if( B(k)% equation(iv) /= FREE_BOUNDARY_CONDITION) then 
                    i = B(k)% i
                    j = B(k)% j
                    Solution(it, i, j, iv) = Y(m)
                    m = m + 1 
                 end if 
     end do 
  end do 
  
  do iv=1, Nv           
     call Derivative( [ "x", "y" ], 1, 1, Solution(it, :, :, iv), Wx(:, :, iv) )
     call Derivative( [ "x", "y" ], 2, 1, Solution(it, :, :, iv), Wy(:, :, iv) )
  end do 
 
 
  m = 1 
  do k=1, size(B) 
       i = B(k) % i 
       j = B(k) % j 
       Gv = Boundary_conditions (x_nodes(i),  y_nodes(j),  t_BC, Solution(it, i, j, :),  Wx(i, j, :),  Wy(i, j, :) )  
       
       do iv=1, Nv 
                   eq = B(k)% equation(iv)
                   if ( eq /= FREE_BOUNDARY_CONDITION .and. eq /= PERIODIC_BOUNDARY_CONDITION) then 
                      G(m) = Gv(iv) 
                      m = m + 1 
                   end if 
       end do 
  end do
      
end function

end subroutine 




end module 
    
    
 
  