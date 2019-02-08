
!***********************************************************************
! It integrates in time the Initial value boundary problem 1D.
! Given the differential operator and the boundary equations,
! the discrete solution is calculated.
! author: Juan A Hernandez, juanantonio.hernandez@upm.es 
!***********************************************************************
module Initial_Value_Boundary_Problem1D

use Cauchy_Problem
use Finite_differences
use Non_Linear_Systems
use Utilities
use Dependencies

implicit none   

private

public :: IVBP1D
public :: IVBP1D_system

   
abstract interface  

       real function DifferentialOperator1D(x, t, u, ux, uxx) 
                        real, intent(in) :: x, t, u, ux, uxx 
       end function  

       real function BC1D(x, t, u, ux) 
           real, intent(in) :: x, t, u, ux 
       end function  

       function DifferentialOperator1D_system(x, t, u, ux, uxx) 
                        real, intent(in) :: x, t, u(:), ux(:), uxx(:) 
                        real ::  DifferentialOperator1D_system( size(u) ) 
       end function  

       function BC1D_system(x, t, u, ux) 
           real, intent(in) :: x, t, u(:), ux(:) 
           real :: BC1D_system( size(u) )  
       end function  


       function DifferentialOperatorODE(U, t) result(F) 
                        real, intent(in) :: U(:), t 
                        real :: F(size(U)) 
       end function 

end interface


    contains

!**********************************************************************************************************************
!
!***********************************************************************************************************************
subroutine IVBP1D( Time_Domain, x_nodes, Order, Differential_operator,  Boundary_conditions, Scheme, Solution) 

     real, intent(in) :: Time_Domain(:) 
     real, intent(inout) :: x_nodes(0:)
     integer, intent(in) :: Order
     procedure (DifferentialOperator1D) :: Differential_operator 
     procedure (BC1D) ::  Boundary_conditions
     procedure (Temporal_Scheme), optional :: Scheme
     real, intent(out) :: Solution(0:,0:) 
     real :: t_BC
     integer ::  it 
     logical :: dU(2) ! matrix of dependencies( direction, order ) 

     integer :: Nx, Nt, Nt_u, Nx_u
     Nx = size(x_nodes) - 1
     Nt = size(Time_Domain) - 1 
     
     Nx_u = size( Solution(0, :) ) - 1 
     Nt_u = size( Solution(:, 0) ) - 1 
     
     dU = Dependencies_IVBP_1D( Differential_operator ) 
     
     if (Nx /= Nx_u .or. Nt /= Nt_u ) then 
         
             write(*,*) " ERROR : Non conformal dimensions Time_Domain, x_nodes, Solution"
             write(*,*) " Nx = ", Nx, " Nt = ", Nt
             write(*,*) " Nx_u = ", Nx_u, " Nt_u = ", Nt_u
             stop 
             
     end if         
          
    call Cauchy_ProblemS(  Time_Domain = Time_Domain,                    & 
                           Differential_operator = Space_discretization, & 
                           Scheme = Scheme, Solution = Solution )
      

contains
!----------------------------------------------------------------------
function Space_discretization( U, t ) result(F) 
          real ::  U(:), t         
          real :: F(size(U))    

         call Space_discretization1D( U, t, F ) 
         
      
end function     
!-----------------------------------------------------------------------
subroutine Space_discretization1D( U, t, F )
          real ::  U(0:Nx), t 
          real ::  F(0:Nx)
              
    integer :: k, N_int 
    real :: Ux(0:Nx), Uxx(0:Nx)
    real, allocatable :: Ub(:)
  
        t_BC = t 
        call Binary_search(t_BC, Time_Domain, it)
        Solution(it, :)  = U(:) 

!  ***  solve a system of two equations to yield the values at the boundaries  
        if (Boundary_conditions( x_nodes(0),  t_BC, 1.,  2.) ==       & 
                                PERIODIC_BOUNDARY_CONDITION) then 
            
            N_int = Nx 
            U(0)  = U(Nx)
           
            
       else if (Boundary_conditions( x_nodes(Nx),  t_BC, 1.,  2.) ==  & 
                               FREE_BOUNDARY_CONDITION) then 
            
            N_int = Nx 
            allocate( Ub(1) ) 
            Ub = [ U(0) ]  
            call Newton( BCs1, Ub )
            U(0)  = Ub(1)
      
        else 
            
            N_int = Nx-1 
            allocate( Ub(2) ) 
            Ub = [ U(0), U(Nx) ] 
            call Newton( BCs2, Ub )
            U(0:Nx:Nx)  = Ub(1:2)
           
        end if 
        
           
!  *** inner grid points  
        F = 0   
        if (dU(1)) call Derivative( "x", 1, U, Ux)
        if (dU(2)) call Derivative( "x", 2, U, Uxx)  
          
        do k=1, N_int
            F(k) = Differential_operator( x_nodes(k), t, U(k), Ux(k), Uxx(k) )
        enddo
        
        if (allocated(Ub)) deallocate( Ub ) 


end subroutine

!-----------------------------------------------------------------------
function BCs2(Y) result(G) 
        real, intent (in) :: Y(:)
        real :: G(size(Y))

    real :: Wx(0:Nx)


    Solution(it, 0)  = Y(1) 
    Solution(it, Nx) = Y(2) 

    call Derivative( "x", 1, Solution(it,:), Wx)

    G(1)= Boundary_conditions( x_nodes(0),  t_BC, Solution(it,0),  Wx(0)  )
    G(2)= Boundary_conditions( x_nodes(Nx), t_BC, Solution(it, Nx), Wx(Nx) )

end function 

!-----------------------------------------------------------------------
function BCs1(Y) result(G) 
        real, intent (in) :: Y(:)
        real :: G(size(Y))

    real :: Wx(0:Nx)

    Solution(it, 0) = Y(1) 

    call Derivative( "x", 1, Solution(it,:), Wx)

    G(1)= Boundary_conditions( x_nodes(0),  t_BC, Solution(it,0),  Wx(0)  )

end function 


end subroutine
!**********************************************************************************************************************
!
!***********************************************************************************************************************
subroutine IVBP1D_system( Time_Domain, x_nodes, Order, N_variables, Differential_operator,  Boundary_conditions, Scheme, Solution) 

     real, intent(in) :: Time_Domain(:) 
     real, intent(inout) :: x_nodes(0:)
     integer, intent(in) :: Order, N_variables 
     procedure (DifferentialOperator1D_system) :: Differential_operator 
     procedure (BC1D_system) ::  Boundary_conditions
     procedure (Temporal_Scheme), optional :: Scheme
     real, intent(out) :: Solution(0:,0:, :) 

     real, pointer :: pSolution(:, :) 

     integer,  allocatable :: iU0(:), iUN(:)
     integer :: N0, N1 
     real :: t_BC 
     integer ::  it 
     logical :: dU(N_variables, 2) 
     
     integer :: Nx, Nt, Nv, Nt_u, Nx_u, Nv_u
     Nx = size(x_nodes) - 1
     Nt = size(Time_Domain) - 1 
     Nv = N_variables 
     
     Nx_u = size( Solution(0, :, 1) ) - 1 
     Nt_u = size( Solution(:, 0, 1) ) - 1 
     Nv_u = size( Solution(0, 0, :) ) 
     
     if (Nx /= Nx_u .or. Nt /= Nt_u .or. Nv /= Nv_u ) then 
         
             write(*,*) " ERROR : Non conformal dimensions Time_Domain, x_nodes, Solution"
             write(*,*) " Nx = ", Nx, " Nt = ", Nt, " Nv = ", Nv 
             write(*,*) " Nx_u = ", Nx_u, " Nt_u = ", Nt_u, " Nv_u =", Nv_u 
             stop 
             
     end if         
     
     if (present(scheme) ) then 
     endif 
  
    dU = .true. 
    
    call Data_pointer( N1 = Nt+1, N2 = Nv*  (Nx+1) ,  U = Solution(0:Nt, 0:Nx, 1:Nv) , pU = pSolution ) 
    call Cauchy_ProblemS(  Time_Domain = Time_Domain, Differential_operator = Space_discretization, Scheme = Scheme, Solution = pSolution )
      
 
    
contains
!----------------------------------------------------------------------
function  Space_discretization( U, t ) result(F) 
          real ::  U(:), t         
          real :: F(size(U))    


         call Space_discretization1D( U, t, F ) 
         
         if (t>0) then 
         end if 
         
end function     
!-----------------------------------------------------------------------
subroutine  Space_discretization1D( U, t, F )
          real ::  U(0:Nx, Nv), t 
          real :: F(0:Nx, Nv)
              
    integer :: i,k
    real :: Ux(0:Nx, Nv), Uxx(0:Nx, Nv)
    real, allocatable ::  Uc(:)  
   
    

!  ***  solve a system of two equations to yield the values at the boundaries

        t_BC = t 


        call Binary_search(t_BC, Time_Domain, it)
        Solution(it, :, :)  = U(:, :)

     
        call Equations_at_boundary( x = x_nodes(0),  W = U(0,:),  Wx = Ux(0,:),  iU = iU0 )
        call Equations_at_boundary( X = x_nodes(Nx), W = U(Nx,:), Wx = Ux(Nx,:), iU =  iUN )
        
        N0 = size(iU0); N1 = size(iUN) 
     
        allocate( Uc(N0+N1) ) 
        do k=1, N0;     Uc(k)    = U(0,  iU0(k) ); end do 
        do k=1, N1;     Uc(k+N0) = U(Nx, iUN(k) ); end do 
      
        call Newton( BCs, Uc )
       
        do k=1, N0;   U(0,  iU0(k) ) = Uc(k);    end do 
        do k=1, N1;   U(Nx, iUN(k) ) = Uc(k+N0); end do 
      
        
!  *** inner grid points
        do i=1, Nv 
           if (dU(i,1)) call Derivative( "x", 1, U(0:,i), Ux(0:,i) )
           if (dU(i,2)) call Derivative( "x", 2, U(0:,i), Uxx(0:,i) )
        end do    
          
        do k=0, Nx
            F(k, :) = Differential_operator( x_nodes(k), t, U(k,:), Ux(k,:), Uxx(k, :) )
        enddo
    
        deallocate(iU0, iUN, Uc) 

end subroutine


!-----------------------------------------------------------------------
function BCs(Y) result(G) 
     real, intent (in) :: Y(:)
     real :: G(size(Y))
     
     real, allocatable :: G0(:), GN(:)   
     integer :: i,k    
     
     real :: Wx(0:Nx, Nv)

      
     do k=1, N0;       Solution(it,  0, iU0(k) ) = Y(k);    end do
     do k=1, N1;       Solution(it, Nx, iUN(k) ) = Y(k+N0); end do 
     
     do i=1, Nv 
           call Derivative( "x", 1, Solution(it, 0:, i), Wx(0:,i) )
     end do    
    
     call Equations_at_boundary( x_nodes(0),  Solution(it,  0, :), Wx( 0,:), G0 )
     call Equations_at_boundary( x_nodes(Nx), Solution(it, Nx, :), Wx(Nx,:), GN )
      
     G = [ G0, GN ] 

     deallocate(G0, GN) 

end function 

!*******************************************************************************************
!*
!*******************************************************************************************
subroutine Equations_at_boundary( x, W,  Wx, F, iU )
                real, intent(in) ::  x, W(:),  Wx(:)
                real, allocatable, optional, intent(out) ::  F(:)
                integer, allocatable, optional, intent(out) ::  iU(:) 

        integer :: i, k,  Nc           
        real :: F1( size(W) )   
         
        
        F1 = Boundary_conditions( x, t_BC, W(:),  Wx(:) ) 
       
        Nc = count(F1 /= FREE_BOUNDARY_CONDITION ) 
    
        if (present(F) ) allocate( F(Nc) ) 
        if (present(iU)) allocate( iU(Nc) ) 
        
        k = 0 
        do i=1, Nv 
            if ( F1(i) /= FREE_BOUNDARY_CONDITION ) then 
              k = k + 1 
              if (present(F))   F(k) = F1(i) 
              if (present(iU))  iU(k) = i 
           end if
        end do  

end subroutine 



end subroutine



end module 
    
    