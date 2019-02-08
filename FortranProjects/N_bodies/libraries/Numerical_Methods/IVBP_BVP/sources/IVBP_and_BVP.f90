!***********************************************************************
module IVBP_and_BVP
!***********************************************************************

use Cauchy_problem
use Finite_differences
use Non_Linear_Systems
use Boundary_value_problems

implicit none   

        
   
abstract interface  

      
 function DifferentialOperator2D_system_mixed(x,  y,   t,               &
                                              u, ux, uy, uxx, uyy, uxy, &
                                              v, vx, vy, vxx, vyy, vxy  ) 
       real, intent(in) :: x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
       real, intent(in) ::          v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
       real :: DifferentialOperator2D_system_mixed (size(u)) 
 end function  

 function BC2D_system(x, y, t, u, ux, uy) 
       real, intent(in) :: x, y, t, u(:), ux(:), uy(:) 
       real :: BC2D_system(size(u)) ! maximum number of BCs at each point
 end function  


end interface




interface  IVBP_and_BVP_problem
                                     module procedure IVBP2D_system_mixed
end interface




contains
    
 
!*****************************************************************
!*
!*****************************************************************
subroutine Data_pointer( N1, N2, U, pU )
     integer, intent(in) :: N1, N2
     real, target, intent(in) :: U(N1, N2)  
     real, pointer, intent(out) :: pU(:,:) 
          
     pU => U  
     
end subroutine 
 
 
!*****************************************************************************
!* It searches the element x in an ordered list V ( x = V(ix) ) 
!* If the element x is not in the list it returns -1 
!*****************************************************************************
subroutine Binary_search(x, V, ix)
    real, intent (in) :: x, V(0:)
    integer, intent (out) :: ix

  integer:: m, upper_limit, lower_limit, N
  N = size(V)-1
  
  
  lower_limit = 0; upper_limit = N;
 
  do while( upper_limit > lower_limit ) 
     
     m = ( upper_limit + lower_limit ) / 2
     
     if ( x > V(m) ) then 
                          lower_limit = m + 1 
     else 
                          upper_limit = m 
     endif 
     
  end do 
 
  
  ix = lower_limit  
  
 
end subroutine
   
    
   

!********************************************************************************************************************************
!*
!********************************************************************************************************************************
subroutine IVBP2D_system_mixed(Time_Domain, x_nodes, y_nodes,          &
                                     Order,     N_u,     N_v,          &
                                     Differential_operator_u,          &
                                     Differential_operator_v,          &
                                       Boundary_conditions_u,          &
                                       Boundary_conditions_v,          &
                                              Scheme, Ut, Vt           ) 

   real, intent(in) :: Time_Domain(:) 
   real, intent(inout) :: x_nodes(0:), y_nodes(0:)
   integer, intent(in) :: Order, N_u, N_v
   procedure (DifferentialOperator2D_system_mixed) ::                  &
                        Differential_operator_u, Differential_operator_v
   procedure (BC2D_system) :: Boundary_conditions_u, Boundary_conditions_v
   procedure (Temporal_Scheme), optional :: Scheme
   real, intent(out) :: Ut(0:, 0:, 0:, :),  Vt(0:, 0:, 0:, :)  
  
   real, pointer :: pUt(:, :) 
   real :: t_BC 
   integer ::  it, Nx, Ny, Nt, Nu, Nv, Nt_u, Nx_u, Ny_u, Nu_u, Nv_u  
   integer ::  M1, M2, M3, M4
   
   real, allocatable :: Ux(:,:,:), Uxx(:,:,:), Uxy(:,:,:)
   real, allocatable :: Uy(:,:,:), Uyy(:,:,:) 
   
   Nx = size(x_nodes) - 1 ; Ny = size(y_nodes) - 1
   Nt = size(Time_Domain) - 1 
   Nu = N_u   ;   Nv = N_v
   
   Nt_u = size( Ut(:, 0, 0, 1) ) - 1
   Nx_u = size( Ut(0, :, 0, 1) ) - 1 
   Ny_u = size( Ut(0, 0, :, 1) ) - 1 
   Nu_u = size( Ut(0, 0, 0, :) ) 
   Nv_u = size( Vt(0, 0, 0, :) ) 
   
   if (Nx /= Nx_u .or. Ny /= Ny_u .or. Nt /= Nt_u .or.                 &
       Nu /= Nu_u .or. Nv /= Nv_u ) then 
       
     write(*,*) " ERROR : Non conformal dimensions" ,                  & 
                " Time_Domain, x_nodes, y_nodes, Solution"
     write(*,*) " Nx = ", Nx, " Ny = ", Ny, " Nt = ", Nt,              &
                " Nu = ", Nu, " Nv = ", Nv 
     write(*,*) " Nx_u = ", Nx_u, " Ny_u = ", Ny_u, " Nt_u = ", Nt_u,  &
                " Nu_u = ", Nu_u, " Nv_u =", Nv_u 
     stop 
           
   end if  
          
   M1 = Nu*(Ny - 1); M2 = Nu*(Ny - 1); M3 = Nu*(Nx - 1); M4 = Nu*(Nx - 1) 
  
   allocate( Ux(0:Nx,0:Ny, Nu), Uxx(0:Nx,0:Ny, Nu), Uxy(0:Nx,0:Ny, Nu) )
   allocate( Uy(0:Nx,0:Ny, Nu), Uyy(0:Nx,0:Ny, Nu) ) 
  
  call Grid_Initialization( "nonuniform", "x", Order, x_nodes ) 
  call Grid_Initialization( "nonuniform", "y", Order, y_nodes ) 
   
!  call Check_grid( "x", x_nodes, Order, Nx+1 ) 
!  call Check_grid( "y", y_nodes, Order, Ny+1 ) 
  
  call Data_pointer( N1 = Nt+1, N2 = Nu* (Nx+1) * (Ny+1),              &
                      U = Ut(0:Nt, 0:Nx, 0:Ny, 1:Nu) , pU = pUt ) 
  
  call Cauchy_ProblemS( Time_Domain = Time_Domain, &
          Differential_operator = Space_discretization, Solution = pUt )
  
   deallocate( Ux, Uxx, Uxy, Uy, Uyy  )
  
contains

!----------------------------------------------------------------------
function Space_discretization( U, t ) result(F) 
          real ::  U(:), t         
          real :: F(size(U))   

         call Space_discretization_2D_system( U, t, F )                        
         
         if (t>0) then 
         end if 

         
end function    

!-----------------------------------------------------------------------
subroutine Space_discretization_2D_system( U, t, F_u )
          real :: U(0:Nx,0:Ny, Nu),  t 
          real :: F_u(0:Nx,0:Ny, Nu)
              
    integer :: i, j, k 
    real :: Vx(0:Nx,0:Ny, Nv), Vxx(0:Nx,0:Ny, Nv), Vxy(0:Nx,0:Ny, Nv)
    real :: Vy(0:Nx,0:Ny, Nv), Vyy(0:Nx,0:Ny, Nv) 
    real ::  Uc(M1 + M2 + M3 + M4)    
  
        t_BC = t 

        call Binary_search(t_BC, Time_Domain, it)
        write(*,*) " it = ", it  
        
 !  *** initial boundary value :  Uc 
        call Asign_BV2s( U( 0, 1:Ny-1, 1:Nu ), U( Nx, 1:Ny-1, 1:Nu ),  & 
                         U( 1:Nx-1, 0, 1:Nu ), U( 1:Nx-1, Ny, 1:Nu ), Uc ) 
        
 !  *** Boundary points Uc from inner points U 
        call Newton( BCs, Uc )
        F_u=0
        
!   *** asign boundary points Uc  to U        
        call Asign_BVs( U( 0, 1:Ny-1, 1:Nu ), U( Nx, 1:Ny-1, 1:Nu  ),  &
                        U( 1:Nx-1, 0, 1:Nu ), U( 1:Nx-1, Ny, 1:Nu ), Uc )         
        
!  *** Derivatives of U for inner grid points
       do k=1, Nu 
          call Derivative( ["x","y"], 1, 1, U(0:,0:, k),  Ux (0:,0:,k) )
          call Derivative( ["x","y"], 1, 2, U(0:,0:, k),  Uxx(0:,0:,k) )
          call Derivative( ["x","y"], 2, 1, U(0:,0:, k),  Uy (0:,0:,k) )
          call Derivative( ["x","y"], 2, 2, U(0:,0:, k),  Uyy(0:,0:,k) )
          call Derivative( ["x","y"], 2, 1, Ux(0:,0:,k),  Uxy(0:,0:,k) )
       end do  
       
       call Boundary_Value_Problem(                                 &
                x_nodes = x_nodes,                                  &
                y_nodes = y_nodes,                                  &
                Order = Order, N_variables = Nv,                    & 
                Differential_operator = Differential_operator_v_R,  &
                Boundary_conditions = Boundary_conditions_v_R,      &
                Solution = Vt(it, 0:, 0:, :)    )
        
!  *** Derivatives for V
    do k=1, Nv 
      call Derivative( ["x","y"], 1, 1, Vt(it, 0:,0:, k), Vx (0:,0:,k) )
      call Derivative( ["x","y"], 1, 2, Vt(it, 0:,0:, k), Vxx(0:,0:,k) )
      call Derivative( ["x","y"], 2, 1, Vt(it, 0:,0:, k), Vy (0:,0:,k) )
      call Derivative( ["x","y"], 2, 2, Vt(it, 0:,0:, k), Vyy(0:,0:,k) )
      call Derivative( ["x","y"], 2, 1, Vx(0:    ,0:, k), Vxy(0:,0:,k) )
    end do      
      
!  *** Differential operator L_u(U,V) at inner grid points

 do i = 1, Nx - 1
   do j = 1, Ny - 1
         
   F_u(i,j,:) = Differential_operator_u( x_nodes(i), y_nodes(j),       &
                                                 t ,   U(i, j,:),      &
                                        Ux(i, j, :),  Uy(i, j, :),     &
                                       Uxx(i, j, :), Uyy(i, j, :),     &
                                       Uxy(i, j, :), Vt(it, i, j, :),  &
                                        Vx(i, j, :), Vy(i, j, :),      &
                                       Vxx(i, j, :), Vyy(i, j, :),     &
                                                     Vxy(i, j, :)      )
   end  do
 end do 
               
end subroutine


!-----------------------------------------------------------------------

function BCs(Y) result(G) 
        real, intent (in) :: Y(:)
        real :: G(size(Y))
    
     real :: G1(M1), G2(M2), G3(M3), G4(M4) 

! ** Asign Newton's iteration Y to Solution     
     call Asign_BVs(Ut(it, 0 ,1:Ny-1, 1:Nu), Ut(it, Nx, 1:Ny-1, 1:Nu), &
                    Ut(it, 1:Nx-1, 0, 1:Nu), Ut(it, 1:Nx-1, Ny, 1:Nu),Y ) 
 
! ** Calculate boundary conditions G      
     call Asign_BCs(  G1,  G2, G3,  G4 ) 
  
     G = [ G1, G2, G3, G4 ] 

end function

!------------------------------------------------------------------------------------
subroutine Asign_BV2s(   U1,  U2,  U3,  U4, Uc  )
   real, intent(in) :: U1(M1), U2(M2), U3(M3), U4(M4) 
   real, intent(out) :: Uc(M1+M2+M3+M4)
   
     Uc = [ U1, U2, U3, U4 ] 
   
end subroutine

!---------------------------------------------------------------------------------
subroutine Asign_BVs( U1,  U2,  U3, U4, Y) 
   real, intent(in) :: Y(M1+M2+M3+M4)
   real, intent(out) :: U1(M1), U2(M2), U3(M3), U4(M4)  
   
    integer :: i1, i2, i3, i4 
  
    i1 = 1  + M1 
    i2 = i1 + M2 
    i3 = i2 + M3 
    i4 = i3 + M4 

     U1 = Y(1  : i1-1)
     U2 = Y(i1 : i2-1) 
     U3 = Y(i2 : i3-1) 
     U4 = Y(i3 : i4-1) 
   
end subroutine


!***********************************************************************************
!*
!***********************************************************************************
subroutine Asign_BCs(  G1, G2,  G3,  G4  ) 
   real, intent(out) :: G1(1:Ny-1,Nu), G2(1:Ny-1,Nu),                  &
                        G3(1:Nx-1,Nu), G4(1:Nx-1,Nu)  
  
   real :: Wx(0:Nx, 0:Ny, Nu), Wy(0:Nx, 0:Ny, Nu)
   integer :: i, j, k 
    
   do k=1, Nu 
    call Derivative( ["x","y"], 1, 1, Ut(it, 0:, 0:, k), Wx(0:, 0:, k) )
    call Derivative( ["x","y"], 2, 1, Ut(it, 0:, 0:, k), Wy(0:, 0:, k) )
   end do 
        
   do j = 1, Ny-1
     G1(j,:) = Boundary_conditions_u ( x_nodes(0), y_nodes(j), t_BC,   &
                               Ut(it, 0, j, : ), Wx(0, j,:), Wy(0, j,:)) 
                                      
     G2(j,:) = Boundary_conditions_u (x_nodes(Nx), y_nodes(j), t_BC,   &
                          Ut(it, Nx, j, : ), Wx(Nx, j, :), Wy(Nx, j, :)) 
   end do
     
   do i = 1, Nx-1
     G3(i,:) = Boundary_conditions_u ( x_nodes(i), y_nodes(0), t_BC,   &
                                Ut(it, i, 0,:), Wx(i, 0, :), Wy(i, 0,:)) 
                                      
     G4(i,:) = Boundary_conditions_u ( x_nodes(i), y_nodes(Ny), t_BC,  &
                          Ut(it, i, Ny, : ), Wx(i, Ny, :), Wy(i, Ny, :)) 
                                  
   end do
   
end subroutine




!-----------------------------------------------------------------------

function Differential_operator_v_R(x, y, V, Vx, Vy, Vxx, Vyy, Vxy) result (Fv)
     real, intent(in) :: x, y, V(:), Vx(:), Vy(:), Vxx(:), Vyy(:), Vxy(:)
     real :: Fv(size(V))      
 
     integer :: ix, iy 
     
     call Binary_search(x, x_nodes, ix)  
     call Binary_search(y, y_nodes, iy)
     
     Fv = Differential_operator_v( x, y, t_BC, V(:),  Vx(:),  Vy(:),   &
                                             Vxx(:), Vyy(:), Vxy(:),   & 
                                  Ut(it, ix, iy, :),  Ux(ix, iy, :),   &
                                      Uy(ix, iy, :), Uxx(ix, iy, :),   &
                                     Uyy(ix, iy, :), Uxy(ix, iy, :)    ) 
      
end function


!-----------------------------------------------------------------------

function Boundary_conditions_v_R(x, y, V, Vx, Vy) result (G_v)
     real, intent(in) :: x, y, V(:), Vx(:), Vy(:)
     real :: G_v(size(V))      
 
    G_v = Boundary_conditions_v (x, y, t_BC, V, Vx, Vy )
    
        
end function






end subroutine   




end module 
    
    
  
