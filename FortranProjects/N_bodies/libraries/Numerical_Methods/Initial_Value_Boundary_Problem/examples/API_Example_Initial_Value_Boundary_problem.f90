module  API_Example_Initial_Value_Boundary_Problem 


use Linear_systems
use Initial_Value_Boundary_Problem
implicit none 


    contains 
    
  !line 10   
subroutine IVBP_examples 
  
  
  implicit none 
  
   call Heat_equation_1D
   call Advection_Diffusion_1D
   
   call Heat_equation_2D
   call Advection_Diffusion_2D
 
   call Waves_equation_1D
   call Waves_equation_2D
   call Plate_vibration 
   
end subroutine 
   














































 
!*************************************************************************************************
!*
!line 77********************************************************************************************
subroutine Waves_equation_1D

       integer, parameter :: Nx = 40, Nt = 1000, Nv = 2 
       real :: Time(0:Nt), U(0:Nt,0:Nx, Nv), x(0:Nx)  
       real ::  x0 = -1, xf = 1, t0 = 0, tf = 2. 
       integer :: i, Order = 4;  
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x(0) = x0; x(Nx) = xf;  
     call Grid_Initialization( "nonuniform", "x", Order, x )
       
     U(0, :,1) = exp( - 15 * x**2 );   U(0, :,2) = 0    
     
     call Initial_Value_Boundary_ProblemS(                                  &  
          Time_Domain = Time, x_nodes = x, Order = Order, N_variables = Nv, & 
          Differential_operator =  Wave_equation1D,                         &
          Boundary_conditions   =  Waves_BC1D,  Solution = U ) 
     
    
contains 


function Wave_equation1D( x, t, u, ux, uxx) result(F)
        real, intent(in) ::  x, t, u(:), ux(:), uxx(:)
        real :: F(size(u))
            
            real :: v, vxx, w 
            
            v = u(1);   w = u(2) 
            vxx = uxx(1) 
        
            F(1)  = w 
            F(2)  = vxx 
                 
end function 
!-------------------------------------------------------
function Waves_BC1D(x, t, u, ux) result(BC) 
  real, intent(in) :: x, t, u(:), ux(:) 
  real :: BC( size(u) ) 
  
       real :: v
             
       v = u(1)
  
        if (x==x0) then
                            BC(1) = v
                            BC(2) = FREE_BOUNDARY_CONDITION 

        else if (x==xf) then
                            BC(1) = v 
                            BC(2) =  FREE_BOUNDARY_CONDITION

        else
             write(*,*)  "Error in BC_Burgers"
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif

end function



end subroutine 
 


!*************************************************************************************************
!*
!******************************************************************************************
subroutine Waves_equation_2D

       integer, parameter :: Nx = 20, Ny = 20, Nt = 100, Nv = 2 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv)  
       
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
       real ::  t0 = 0, tf = 2 
       integer :: i, j, Order = 8 
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x(0) = x0;  x(Nx) = xf; y(0) = y0;  y(Ny) = yf;
     
     call Grid_Initialization( "nonuniform", "x", Order, x )
     call Grid_Initialization( "nonuniform", "y", Order, y )

  
     U(0, :, :, 1) = Tensor_product( exp(-10*x**2) , exp(-10*y**2) )  
     U(0, :, :, 2) = 0
  
     call Initial_Value_Boundary_ProblemS(                                 & 
     
                                Time_Domain = Time,                        &
                                x_nodes = x, y_nodes = y,                  & 
                                Order = Order, N_variables = Nv,           &    
                                Differential_operator = Wave_equation2D,  & 
                                Boundary_conditions = Wave_BC2D,           &
                                Solution = U ) 
     
     call scrmod('revers') 
     call qplclr( U(0, 0:Nx, 0:Ny, 1) ,  Nx+1, Ny+1)
     call qplclr( U(Nt, 0:Nx, 0:Ny, 1) , Nx+1, Ny+1) 
     

contains 



function Wave_equation2D( x, y, t, u, ux, uy, uxx, uyy, uxy ) result(L)
      real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
      real :: L(size(u))
            
            real :: v, vxx, vyy, w 
            
            v = u(1);   w = u(2) 
            vxx = uxx(1) ; vyy = uyy(1)
        
            L(1)  = w 
            L(2)  = vxx +vyy
                 
end function 



!-------------------------------------------------------
function Wave_BC2D(x,y, t, u, ux, uy) result(BC) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BC( size(u) ) 
  
       real :: v, w 
            
        v = u(1)
        w = u(2) 
  
        if (x==x0) then
                            BC(1) = v
                            BC(2) = w 
        else if (x==xf) then
                            BC(1) = v 
                            BC(2) = w
        else if (y==y0) then
                            BC(1) = v
                            BC(2) = w 

        else if (y==yf) then
                            BC(1) = v 
                            BC(2) = w 

        else
             write(*,*)  "Error in BC2D_waves"
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif
      
end function

end subroutine 






















































































































!line 753------------------------------------------------------
subroutine Heat_equation_1D

       integer, parameter :: Nx = 20, Nt = 1000
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx)  
       
       real ::  x0 = -1, xf = 1
       real :: t0 = 0, tf = 4  
       integer :: i, j, Order = 6 
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x(0) = x0; x(Nx) = xf 
     call Grid_Initialization( "nonuniform", "x", Order, x )
     
      
     U(0, :)  =  exp(-25*x**2 )
  
       
     call Initial_Value_Boundary_ProblemS(                                 & 
         
                       Time_Domain = Time, x_nodes = x,  Order = Order,    & 
                       Differential_operator =  Heat_equation1D,           & 
                       Boundary_conditions   =  Heat_BC1D,                 & 
                       Solution = U ) 
     
     call scrmod("reverse")
     call qplot(x, U(0,:),  Nx+1)  
     call qplot(x, U(Nt,:), Nx+1)    
contains 
!line 783----------------------------------------------------------
real function Heat_equation1D( x, t, u, ux, uxx) result(F)
        real, intent(in) ::  x, t, u, ux, uxx
            
            real :: nu = 0.02
            
            F =  nu * uxx
                 
end function 
!-------------------------------------------------------
real function Heat_BC1D(x, t, u, ux) result(BC) 
  real, intent(in) :: x, t, u, ux 
  
        if (x==x0) then
                            BC = u 

        else if (x==xf) then
                            BC = u

        else
             write(*,*)  "Error in Heat_BC"
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif

end function


end subroutine 
 

!*************************************************************************************************
!*
!line 817********************************************************************************************
subroutine Advection_diffusion_1D

       integer, parameter :: Nx = 30, Nt = 1000
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx)  
       
       real ::  x0 = -1, xf = 1
       real :: t0 = 0, tf = 0.5  
       integer :: i, Order = 4 
        
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]
     x(0) = x0; x(Nx) = xf 
     call Grid_Initialization( "nonuniform", "x", Order, x )
    
  
     
     U(0, :)  =  exp( -25*x**2 )
       
     call Initial_Value_Boundary_ProblemS(                                 & 
         
                 Time_Domain = Time, x_nodes = x, Order = Order,           &
                 Differential_operator =  Advection_equation1D,            &
                 Boundary_conditions   =  Advection_BC1D,  Solution = U ) 
      
     call scrmod("reverse")  
     call qplot(x, U(0,:),  Nx+1)
     call qplot(x, U(Nt,:), Nx+1)    
                             
contains 
!line 849 -------------------------------------------------------
real function Advection_equation1D( x, t, u, ux, uxx) result(F) 
        real, intent(in) ::  x, t, u, ux, uxx
            
            real :: nu = 0.02
            
            F = - ux + nu * uxx
end function 

!-------------------------------------------------------
real function Advection_BC1D(x, t, u, ux) result(BC) 
  real, intent(in) :: x, t, u, ux 
  
        if (x==x0) then
                            BC = u 

        else if (x==xf) then
                            BC = u
        else
             write(*,*)  "Error in Advection_BC"
             write(*,*) " x = ", x 
             stop 
        endif

end function



end subroutine 
 


!*************************************************************************************************
!*
!line 883*********************************************************************************************
subroutine Advection_diffusion_2D

      integer, parameter :: Nx = 20, Ny = 20, Nt = 200
      real :: x(0:Nx), y(0:Ny), Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny) 
       
       real :: x0 = -1, xf = 1
       real :: y0 = -1, yf = 1
       real :: t0 = 0, tf = 1.0
       integer :: i, j, Order = 8
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x(0)= x0; x(Nx) = xf; y(0)=y0; y(Ny) = yf 
     call Grid_Initialization( "nonuniform", "x", Order, x )
      call Grid_Initialization( "nonuniform", "y", Order, y )
  
     U(0, :, :) = Tensor_product( exp(-25*x**2), exp(-25*y**2) ) 
       
     call Initial_Value_Boundary_ProblemS(                                  & 
               Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, &
               Differential_operator =  Advection_equation2D,               & 
               Boundary_conditions   =  Advection_BC2D,  Solution = U ) 
     
     call scrmod("reverse")  
     call qplcon( U(0,:, :),   Nx+1, Ny+1, 20)
     call qplcon( U(Nt, :, :), Nx+1, Ny+1, 20)

contains





!line 916-----------------------------------------------------
function Advection_equation2D( x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy ) result(F) 
        real,intent(in) :: x, y, t 
        real, intent(in) ::   U, Ux, Uy, Uxx, Uyy, Uxy
        real :: F 
        
        real :: nu = 0.02

        F = - Ux + nu * ( Uxx + Uyy )

end function


!-------------------------------------------------------
function Advection_BC2D( x, y, t, U, Ux, Uy ) result (BC) 

      real, intent(in) :: x, y, t 
      real, intent(in) :: U, Ux, Uy
      real :: BC 

        if (x==x0) then
                               BC = U  
        else if (x==xf) then
                               BC = U
        else if (y==y0) then
                               BC = U
        else if (y==yf) then
                               BC = U
        else
            Write(*,*)  "Error in Advection_BC "
        end if

end function

end subroutine  

!*************************************************************************************************
!* line 953
!*************************************************************************************************  
subroutine Heat_equation_2D

      integer, parameter :: Nx = 20, Ny = 20, Nt = 200
      real :: x(0:Nx), y(0:Ny), Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny) 
       
       real :: x0 = -1, xf = 1, y0 = -1, yf = 1  
       real :: t0 = 0, tf = 4
       integer :: i, j, Order = 4       
      
    Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
    x(0) = x0;  x(Nx) = xf; y(0) = y0;  y(Ny) = yf;
    
    call Grid_Initialization( "nonuniform", "x", Order, x )
    call Grid_Initialization( "nonuniform", "y", Order, y )
   
    U(0, :, :) = Tensor_product( exp(-25*x**2), exp(-25*y**2) ) 
    
    call Initial_Value_Boundary_ProblemS(                                    & 
                Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, &
                Differential_operator =  Heat_equation2D,                    & 
                Boundary_conditions   =  Heat_BC2D,  Solution = U ) 
    
    call scrmod("reverse")  
    call qplcon( U(0,:, :),   Nx+1, Ny+1, 20)
    call qplcon( U(Nt, :, :), Nx+1, Ny+1, 20)
     
contains


function Heat_equation2D( x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy ) result(F) 
        real,intent(in) :: x, y, t 
        real, intent(in) ::   U, Ux, Uy, Uxx, Uyy, Uxy
        real :: F 
             real :: nu = 0.02
        
             F =   nu * ( Uxx + Uyy )
end function

!-------------------------------------------------------
function Heat_BC2D( x, y, t, U, Ux, Uy ) result (BC) 

      real, intent(in) :: x, y, t 
      real, intent(in) :: U, Ux, Uy
      real :: BC 

        if (x==x0) then
                               BC = U  
        else if (x==xf) then
                               BC = U
        else if (y==y0) then
                               BC = U
        else if (y==yf) then
                               BC = U
        else
            write(*,*)  "Error in Advection_BC "
        end if

end function

end subroutine








subroutine Plate_Vibration

       integer, parameter :: Nx = 10, Ny = 10, Nt = 200, Nv = 3 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv)   
       real ::  x0 = -1 , xf = 1 , y0 = -1 , yf = 1
       real ::  t0 = 0 , tf  , pi
       integer :: i, j,m, Order = 6  
       real :: thickness = 16d-3, rho_g = 2500
       real :: E = 72d6 , poisson = 0.22 
       real :: mu, lambda, D  
       
       pi = 4 * atan(1d0)
       tf = 1
       Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]   
               
       x(0) = x0;  x(Nx) = xf 
       y(0) = y0;  y(Ny) = yf  
   
     call Grid_Initialization( "nonuniform", "x", Order, x )
     call Grid_Initialization( "nonuniform", "y", Order, y )

      D = E * thickness**3 / (12 * (1- poisson**2)) 
      
      U(0, :, :, 1) =  0; U(0, :, :, 3) =  0
      
      U(0, :, :, 2) = Tensor_product( exp(-10*x**2), exp(-10*y**2) ) 
     
      call Initial_Value_Boundary_ProblemS( Time_Domain = Time,             & 
                                   x_nodes = x, y_nodes = y,                & 
                                   Order = Order, N_variables = Nv,         & 
                                   Differential_operator = Plate_equation,  & 
                                   Boundary_conditions = Plate_BC,          &
                                   Solution = U ) 
      
     call scrmod("reverse")  
     call qplcon( U(0,:, :, 2),   Nx+1, Ny+1, 20)
     call qplcon( U(Nt, :, :, 2), Nx+1, Ny+1, 20)
     
contains 

function Plate_equation( x, y, t, u, ux, uy, uxx, uyy, uxy ) result(L)
    real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
    real :: L(size(u))
        
            L(1)  =   u(2) 
            L(2)  = -( uxx(3) + uyy(3) ) * D / ( rho_g*thickness )
            L(3)  =  ( uxx(2) + uyy(2) )
                 
end function 

function Plate_BC( x, y, t, u, ux, uy ) result(BC) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BC( size(u) ) 
         
        if (x==x0) then
                           BC(1) = u(1)
                           BC(2) = u(2)
                           BC(3) = u(3)

        else if (x==xf) then
                           BC(1) = u(1)
                           BC(2) = u(2)
                           BC(3) = u(3)

        else if (y==y0) then
                           BC(1) = u(1)
                           BC(2) = u(2)
                           BC(3) = u(3)

        else if (y==yf) then
                           BC(1) = u(1)
                           BC(2) = u(2)
                           BC(3) = u(3)
        else
             write(*,*)  "Error in BC1 "
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif

end function


end subroutine


end module 
