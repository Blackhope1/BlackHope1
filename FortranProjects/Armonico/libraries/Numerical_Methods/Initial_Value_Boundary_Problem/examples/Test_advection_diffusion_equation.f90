module Test_advection_diffusion_equation 


use Initial_Value_Boundary_Problem
implicit none 


contains 

!*************************************************************************************************
!*
!*************************************************************************************************
subroutine Test_IVBP1D

       integer, parameter :: Nx = 60, Nt = 1000
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx)  
       
       real ::  x0 = -1, xf = 1
       real :: t0 = 0, tf = 3  
       integer :: i, Order = 4 
        
     U = 0. 
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
       
       
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, Order = Order, & 
                                           Differential_operator =  Burgers_equation, & 
                                           Boundary_conditions   =  Burgers_BC,  Solution = U ) 
      
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $  \frac{\partial{u}}{\partial{t}} + u \frac{\partial{u}}{\partial{x}}  = \nu \frac{\partial{^2 u}}{\partial{x^2}}  $ ", 2)
     CALL TITLIN (" $ y(0) = 1, \quad y_{x}(1) = - \pi $  ", 4)
     call qplot(x, U(Nt,:), Nx+1)                                 


contains 
!----------------------------------------------------------
real function Burgers_equation( x, u, ux, uxx) 
        real, intent(in) ::  x, u, ux, uxx
            
            real :: nu = 0.01
            
            Burgers_equation = - u * ux + nu * uxx
                 
end function 
!-------------------------------------------------------
real function Burgers_BC(x, u, ux)
  real, intent(in) :: x, u, ux 
  
        if (x==x0) then
                            Burgers_BC = u - 1

        else if (x==xf) then
                            Burgers_BC = ux

        else
             write(*,*)  "Error in BC_Burgers"
        endif

end function



end subroutine 
 
 

    
!*************************************************************************************************
!*
!*************************************************************************************************  
subroutine Test_IVBP2D

      integer, parameter :: Nx = 10, Ny = 10, Nt = 200
      real :: x(0:Nx), y(0:Ny), Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny) 
       
       real :: x0 = -1, xf = 1
       real :: y0 = -1, yf = 1
       real :: t0 = 0, tf = 2
       integer :: i, j, Order = 8
       character(len=500) :: title 
       
       title = " $  Ut  + U  Ux   = \nu ( Uxx + Uyy ) $ "
       
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*i/Ny, i=0, Ny ) ]
       
     do i=0,Nx 
       do j=0, Ny 
            U(0, i, j)  =  exp( -25*x(i)**2  -25*y(j)**2 )
       end do 
     end do
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TITLIN (" Initial condition for the advection diffusion equation ", 2)
     call qplcon( U(0, :, :), Nx+1, Ny+1, 20);
       
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, & 
                                           Differential_operator =  Advection_equation, & 
                                           Boundary_conditions   =  Advection_BC,  Solution = U ) 
      
     
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN ( title, 2) 
     CALL TITLIN (" $ U = 0 $ at the boundary of $[-1,1] \times [-1,1] $  ", 4)
     call qplcon( U(Nt, :, :), Nx+1, Ny+1, 20);
                                   
 

contains

!----------------------------------------------------------
function Advection_equation( x, y, U, Ux, Uy, Uxx, Uyy, Uxy ) result(F) 
        real,intent(in) :: x, y
        real, intent(in) ::   U, Ux, Uy, Uxx, Uyy, Uxy
        real :: F 
        
        real :: nu = 0.02

        F = -U * Ux + nu * ( Uxx + Uyy )
        
   

end function
!-------------------------------------------------------
function Advection_BC( x, y, U, Ux, Uy ) result (BC) 

      real, intent(in) :: x, y
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
 


end module 
