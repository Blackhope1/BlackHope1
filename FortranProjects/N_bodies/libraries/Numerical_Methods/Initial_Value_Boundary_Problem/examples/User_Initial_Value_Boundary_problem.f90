module  API_Example_Initial_Value_Boundary_problem 


use Initial_Value_Boundary_Problem
use Finite_differences

implicit none 


contains 

!*************************************************************************************************
!*
!*************************************************************************************************
subroutine Test_IVBP1D

       integer, parameter :: Nx = 80, Nt = 1000
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx)  
       
       real ::  x0 = 0, xf = 1
       real :: t0 = 0, tf = 1  
       integer :: i, Order = 6 
        
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
     CALL TITLIN (" $ u(0) = 1, \quad u_{x}(1) = 0 $  ", 4)
     call qplot(x, U(Nt,:), Nx+1)                                 


contains 
!----------------------------------------------------------
real function Burgers_equation( x, t, u, ux, uxx) 
        real, intent(in) ::  x, t, u, ux, uxx
            
            real :: nu = 0.005
            
            Burgers_equation = - u * ux + nu * uxx
                 
end function 
!-------------------------------------------------------
real function Burgers_BC(x, t, u, ux)
  real, intent(in) :: x, t, u, ux 
  
        if (x==x0) then
                            Burgers_BC = u - 1

        else if (x==xf) then
                            Burgers_BC = ux

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
!*************************************************************************************************
subroutine Test_IVBP1D_system

       integer, parameter :: Nx = 20, Nt = 1000, Nv = 2 
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx, Nv)  
       
       real ::  x0 = -1, xf = 1
       real :: t0 = 0, tf = 4.5 
       integer :: i, Order = 4 
        
     
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
       
     U(0, :,1) = exp( - 15 * x**2 )    
     U(0, :,2) = 0 
     
     
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, Order = Order, N_variables = Nv, & 
                                           Differential_operator =  Wave_equation, & 
                                           Boundary_conditions   =  Waves_BC,  Solution = U ) 
      
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $  \frac{\partial{^2 u}}{\partial{t^2}}   = \frac{\partial{^2 u}}{\partial{x^2}}  $ ", 2)
     CALL TITLIN (" $ u(-1) = 0, \quad u(1) = 0 $  ", 4)
     
     call qplot(x, U(0, :, 1), Nx+1)                                  
     call qplot(x, U(Nt,:, 1), Nx+1)

contains 
!----------------------------------------------------------
function Wave_equation( x, t, u, ux, uxx) 
        real, intent(in) ::  x, t, u(:), ux(:), uxx(:)
        real :: Wave_equation(size(u))
            
            real :: v, vxx, w 
            
            
            v = u(1);   w = u(2) 
            vxx = uxx(1) 
        
        
            Wave_equation(1)  = w 
            Wave_equation(2)  = vxx 
       
                 
end function 
!-------------------------------------------------------
function Waves_BC(x, t, u, ux)
  real, intent(in) :: x, t, u(:), ux(:) 
  real :: Waves_BC( size(u) ) 
  
       real :: v
            
       v = u(1)
       
  
        if (x==x0) then
                            Waves_BC(1) = v
                            Waves_BC(2) = 1d100 

        else if (x==xf) then
                            Waves_BC(1) = v 
                            Waves_BC(2) =  1d100

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
!*************************************************************************************************  
subroutine Test_IVBP2D

      integer, parameter :: Nx = 10, Ny = 10, Nt = 200
      real :: x(0:Nx), y(0:Ny), Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny) 
       
       real :: x0 = 0, xf = 1
       real :: y0 = 0, yf = 1
       real :: t0 = 0, tf = 4
       integer :: i, j, Order = 8
       character(len=500) :: title 
       
       title = " $  Ut  + U  Ux   = \nu ( Uxx + Uyy ) $ "
       
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*i/Ny, i=0, Ny ) ]
       
     do i=0,Nx 
       do j=0, Ny 
            U(0, i, j)  =  exp( -25*(x(i)- (x0+xf)/2 ) **2  -25*(y(j) - (x0+xf)/2 )  **2 )
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
     CALL TITLIN (" $ U = 0 $ at the boundary of $[0,1] \times [0,1] $  ", 4)
     call qplcon( U(Nt, :, :), Nx+1, Ny+1, 20);
                                   
 

contains

!----------------------------------------------------------
function Advection_equation( x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy ) result(F) 
        real,intent(in) :: x, y, t 
        real, intent(in) ::   U, Ux, Uy, Uxx, Uyy, Uxy
        real :: F 
        
        real :: nu = 0.02

        F = -U * Ux + nu * ( Uxx + Uyy )
        
   

end function
!-------------------------------------------------------
function Advection_BC( x, y, t, U, Ux, Uy ) result (BC) 

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
!*
!*************************************************************************************************  
subroutine Test_Heat2D

      integer, parameter :: Nx = 15, Ny = 15, Nt = 200
      real :: x(0:Nx), y(0:Ny), Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny) 
       
       real :: x0 = 0, xf = 1
       real :: y0 = 0, yf = 1
       real :: t0 = 0, tf = 2
       integer :: i, j, Order = 4
       character(len=500) :: title 
       real :: zmax = 2., zmin = 0. 
       
       title = " $  Ut  = \nu ( Uxx + Uyy ) $ "
       
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*i/Ny, i=0, Ny ) ]
     
           
     !do i=0,Nx 
     !  do j=0, Ny 
     !       U(0, i, j)  =  0; !exp( -25*(x(i)- (x0+xf)/2 ) **2  -25*(y(j) - (x0+xf)/2 )  **2 )
     !  end do 
     !end do
     U(0,:,:) = 0 
     
    
           
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, & 
                                           Differential_operator =  Heat_equation, & 
                                           Boundary_conditions   =  Heat_BC,  Solution = U ) 
     ! 
     !call scrmod('revers')                          
     !call metafl('xwin')
     !call disini() 
     !CALL TEXMOD('ON')
     !CALL TITLIN ( title, 2) 
     !CALL TITLIN (" $ U = 0 $ at the boundary of $[0,1] \times [0,1] $  ", 4)
     !CALL GRAF(x0, xf, x0, (xf-x0)/5, y0, yf, y0, (yf-y0)/5 )
     !CALL LABELS ('float', 'CONTUR') 
     !write(*,*) maxval( U(Nt,:,:) ) 
   !  do i=0, 10       
    ! CALL CONTUR (x, Nx+1,  y, Ny+1, U(Nt,:,:), 1.5*i)
     
     !write(*,*) " xmax, xmin  =", maxval(x), minval(x) 
     !write(*,*) " ymax, ymin  =", maxval(y), minval(y) 
     CALL qplclr (U(Nt,:,:), Nx+1, Ny+1 )
   !  end do 
   !  call disfin
                                   
 

contains

!----------------------------------------------------------
function Heat_equation( x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy ) result(F) 

     real,intent(in) :: x, y, t 
     real, intent(in) ::   U, Ux, Uy, Uxx, Uyy, Uxy  
     real :: F 

     
        real :: nu = 0.02
        logical :: sx, sy, inner_region 
        
   !  write(*,*) x, y 
        
        sx = x >= 1 .and. x <= 3
        sy = y >= 0.2 .and. y <= 0.8
        inner_region = sx .and. sy ! .false. ! 
        
       

        !if (inner_region) then 
        !    F =  0; ! ( Ux  ) 
        !else     
            F =  nu * ( Uxx + Uyy )
        !end if     
        
   

end function
!-------------------------------------------------------
function Heat_BC( x, y, t, U, Ux, Uy ) result (BC) 

      real, intent(in) :: x, y, t 
      real, intent(in) :: U, Ux, Uy
      real :: BC 

        if (x==x0) then
                               BC = U - 15   
        else if (x==xf) then
                               BC = U - 15 
        else if (y==y0) then
                               BC = U - 15 
        else if (y==yf) then
                               BC = U - 15 
        else
            Write(*,*)  "Error in Heat_BC "
        end if

end function

end subroutine
 
!*************************************************************************************************
!*
!*************************************************************************************************
subroutine Test_pointers2 


   integer, parameter :: N = 10 
   real, target :: U(0:N, 0:N) 
    
  
   type real_pointer
      real, pointer :: p 
   end type 
   
   type (real_pointer) :: Uc((N-1)*4)

   
   real, pointer :: pU(:), p1(:),  p2(:), p3(:), p4(:)   
   real :: s 
   
   integer :: i 
   
 
   
   U = 1000 
   U(0,:) = 100
   U(:, N) = 200
   U(N, :) = 300 
   U(:, 0) = 400
   
   do i=1, N-1 
      Uc(i) % p         => U(0, i) 
      Uc(N-1+i) % p     => U(i, N) 
      Uc(2*(N-1)+i) % p => U(N, i) 
      Uc(3*(N-1)+i) % p => U(i, 0) 
   end do 
   
   write(*,*) Uc(1) % p, Uc(N) % p, Uc( 2*N ) % p, Uc(4*(N-1)) % p 
   
   do i=1, 4*(N-1) 
      Uc(i) % p = Uc(i) % p + 1
   end do 
     
   do i=0, N
       write(*, '(100f8.0)' ) U(i, :) 
   end do 
   
  ! pU =>
   p1 =>  U(0, 1:N-1) 
   p2 =>  U(1:N-1, N)
   p3 =>  U(N, 1:N-1)
   p4 =>  U(1:N-1, 0) 
   
  
   
!   call Newtonv(  p1, p2, p3, p4  ) 
 
   
   write(*,*)
    do i=0, N
       write(*, '(100f8.0)' ) U(i, :) 
    end do  
 
  

  
   
end subroutine 



!*************************************************************************************************
!*
!*************************************************************************************************
subroutine Test_pointers 


   integer, parameter :: N = 10 
   real, target :: U(0:N, 0:N) 
    
   
   real, pointer :: p1(:), p2(:), p3(:), p4(:)  
   
   real :: f1(N-1), f2(N-1), f3(N-1), f4(N-1)
   real :: s 
   
   integer :: i 
   
   U = 1000 
   U(0,:) = 100
   U(:, N) = 200
   U(N, :) = 300 
   U(:, 0) = 400
    
   do i=0, N
       write(*, '(100f8.0)' ) U(i, :) 
   end do 

   p1 =>  U(0, 1:N-1) 
   p2 =>  U(1:N-1, N)
   p3 =>  U(N, 1:N-1)
   p4 =>  U(1:N-1, 0) 
   
   f1 = 100 
   f2 = 100 
   f3 = 100 
   f4 = 100 
   
   
   
   call Newtonv( p1, f1, p2, f2, p3, f3, p4, f4 ) 
 
   
   write(*,*)
    do i=0, N
       write(*, '(100f8.0)' ) U(i, :) 
    end do  
 
  

  
   
end subroutine 

 
!*************************************************************************************************
!*
!*************************************************************************************************
subroutine Newtonv(x1, F1, x2, F2, x3, F3, x4, F4) 
          real, intent(inout) :: x1(:), x2(:), x3(:), x4(:) 
          real, intent(in) :: F1(:), F2(:), F3(:), F4(:) 
   
   real :: suma 
   integer :: n
   
   real :: x(4*size(x1)), F(4*size(x1))  
   
   x = [ x1, x2, x3, x4 ] 
   F = [F1, F2, F3, F4 ] 
   
   call Newtonr(x, F)  
   
   n = size(x1) 
   
   x1 = x(1:n) 
   x2 = x(n+1:2*n) 
   x3 = x(2*n+1:3*n) 
   x4 = x(3*n+1:4*n) 
   
  
   
end subroutine 

!*************************************************************************************************
!*
!*************************************************************************************************
subroutine Newtonr(x, F) 
          real, intent(inout) :: x(:) 
          real, intent(in) :: F(:) 
   
   write(*,*) " F = ", F 
   write(*,*) " X = ", X 
             
   X = X - F 
  
   
end subroutine 



!*************************************************************************************************
!*
!*************************************************************************************************
subroutine Test_IVBP2D_system

       integer, parameter :: Nx = 30, Ny = 30, Nt = 200, Nv = 2 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv)  
       
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
       real ::  t0 = 0, tf = 10 
       integer :: i, j, Order = 4  
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*j/Ny, j=0, Ny ) ]

    
   
     
  do i = 0, Nx
    do j = 0, Nx
     U(0, i, j, 1) = exp( - 4 *( x(i)**2 + y(j)**2 ) )    
     U(0, i, j, 2) = 0   
    end do  
  end do
    
     
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, N_variables = Nv, & 
                                           Differential_operator =  Wave_equation_2D, & 
                                           Boundary_conditions   =  Wave_BC2D,  Solution = U ) 
     
     
     
   
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $  \frac{\partial{^2 u}}{\partial{t^2}}   = \frac{\partial{^2 u}}{\partial{x^2}} + \frac{\partial{^2 u}}{\partial{y^2}}  $ ", 2)
     CALL TITLIN (" $ u(-1) = 0, \quad u(1) = 0 $  ", 4)

    call qplclr( U(0 ,:,:, 1), Nx+1, Ny+1)   
 
     call qplclr( U(0.25*Nt ,:,:, 1), Nx+1, Ny+1)    
     
     call qplclr( U(0.5*Nt ,:,:, 1), Nx+1, Ny+1)    
    
     call qplclr( U(0.75*Nt ,:,:, 1), Nx+1, Ny+1)    
                                    
     call qplclr( U(Nt,:,:, 1), Nx+1, Ny+1) 

contains 
!----------------------------------------------------------
function Wave_equation_2D( x, y, t, u, ux, uy, uxx, uyy, uxy ) result(L)
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

subroutine Linear_Plate_Vibrations_Validation

       integer, parameter :: Nx = 20, Ny = 20, Nt = 800, Nv = 3 
       real ::  x(0:Nx), y(0:Ny),Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny,Nv)  
       real ::  x0 , xf , y0 , yf
       real ::  t0 = 0 , tf  , pi
       integer :: i, j,m, Order = 4  

       pi = 4 * atan(1d0)
       tf = 1

       x0 = -1 ;  y0 = -1
       xf = 1 ;  yf =  1  
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*j/Ny, j=0, Ny ) ]

      D = E * (thick**3) / (12 * (1- poisson**2)) 

  do i = 0, Nx
   do j= 0, Ny
  U(0,i,j,1) =  1d-3*sin(pi*x(i))*sin(pi*y(j))
  U(0,i,j,2) =  0
  U(0,i,j,3) = 1d-3* (-2)*(pi**2) *sin(pi*x(i))*sin(pi*y(j))
   end do
  end do


     
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, N_variables = Nv, & 
                                           Differential_operator =  L1, Boundary_conditions   =  BC1, Solution = U ) 
   

contains 

!----------------------------------------------------------
function L1( x, y, t, u, ux, uy, uxx, uyy, uxy   ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real :: L1(size(u))
            
           
        
            L1(1)  =   u(2) 
            L1(2)  = -(uxx(3) + uyy(3)) + 1d-3*( (4*pi**4)   - (4*pi**2)  ) *sin(pi*x) * sin(pi*y)* cos(2*pi*t)
            L1(3)  =  ( uxx(2) + uyy(2) )
       
                 
end function 
!-------------------------------------------------------
function BC1(x,y, t, u, ux, uy) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BC1( size(u) ) 
 
  
        if (x==x0) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (x==xf) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (y==y0) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (y==yf) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)
        else
             write(*,*)  "Error in BC1 "
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif
      

end function


!----------------------------------------------------------
end subroutine


subroutine Linear_Vibrations_Glazings


    integer, parameter :: Nx = 20, Ny = 20, Nt = 500, Nv = 3 
    real ::  x(0:Nx), y(0:Ny)
    real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv)
    real :: Ux(0:Nt, 0:Nx, 0:Ny, Nv), Uy(0:Nt, 0:Nx, 0:Ny, Nv)
    real :: Uxx(0:Nt, 0:Nx, 0:Ny, Nv), Uyy(0:Nt, 0:Nx, 0:Ny, Nv), Uxy(0:Nt, 0:Nx, 0:Ny, Nv)
    real :: Sxx(0:Nt, 0:Nx, 0:Ny), Syy(0:Nt, 0:Nx, 0:Ny), Sxy(0:Nt, 0:Nx, 0:Ny)
    real :: thick = 8d-3, Lx = 1.3, Ly = 3, rho_g = 2500, d_0 = 16d-3
    real :: x0 , xf , y0 , yf , tf, t0 = 0
    integer :: i, j,l
    integer, parameter :: Order = 6
    real :: E = 72d9 , poisson = 0.22 
    real :: mu, lambda, D, rho_eq, D_eq, pi
    real :: rho_water = 1000, p_atm = 101325 , g = 9.817
    real :: dmin 
    
    
    x0 = 0 ; xf =  Lx
    y0 = 0 ; yf =  Ly

    pi = acos(-1d0)
    
    
    !Malla equiespaciada
       tf = 0.00001
      
    Time = [ (t0 + (tf-t0)*l/Nt, l=0, Nt ) ]     
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]

    
    
    mu= 0.5 * E / (1+poisson) 
    lambda = E * poisson / ((1 - 2*poisson)*(1 + poisson))
    D = E * (thick**3) / (12 * (1- poisson**2))
    rho_eq = 2* rho_g * thick + rho_water* d_0
  D_eq = 2*D
! do i = 0, Nx
!  do j= 0, Ny
! U(0,i,j,1) = 0
! U(0,i,j,2) =  exp(-15*((x(i)-0.5*Lx)**2 + (y(j)-0.5*Ly)**2 ))
! U(0,i,j,3) =  0
!  end do
! end do
! do i = 0, Nx
!  do j= 0, Ny
! U(0,i,j,1) = 0
! U(0,i,j,2) =  0
! U(0,i,j,3) =  0
!  end do
! end do  
  
  do i = 0, Nx
   do j= 0, Ny
  U(0,i,j,1) =  sin(pi*x(i)/(Lx*0.5)) * sin(pi*y(j)/(Ly*0.5))
  U(0,i,j,2) =  0
  U(0,i,j,3) =  -(pi**2) *sin(pi*x(i)/(Lx*0.5))*sin(pi*y(j)/(Ly*0.5))* (4/(Lx**2) + 4/(Ly**2))
   end do
  end do
   
     write(*,*) " example definition OK "
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini()     
     CALL TEXMOD('ON')
     CALL TITLIN (" Introduced Initial value for U(1)  ", 2)
      call qplclr( U(0 ,:,:, 1), Nx+1, Ny+1)   

     
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, N_variables = Nv, & 
                                           Differential_operator =  L1, Boundary_conditions   =  BC1, Solution = U ) 
     

     
   call Check_grid( "x", x, Order, size( U(0,:,0,1) ) ) 
   call Check_grid( "y", y, Order, size( U(0,0,:,1) ) )
  
  do l = 0, Nt
     do i=1, Nv 
        call Derivative( ["x","y"], 1, 1, U (l, 0:,0:,i), Ux (l, 0:,0:,i) )
        call Derivative( ["x","y"], 1, 2, U (l, 0:,0:,i), Uxx(l, 0:,0:,i))
                                              
        call Derivative( ["x","y"], 2, 1, U (l, 0:,0:,i), Uy (l, 0:,0:,i)  )
        call Derivative( ["x","y"], 2, 2, U (l, 0:,0:,i), Uyy(l, 0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, Ux(l, 0:,0:,i), Uxy(l, 0:,0:,i))
     end do   
  end do
 Sxx =  -1d-6 * D * (Uxx(:,:,:,1) + poisson * Uyy(:,:,:,1) ) * 6 / (thick)**2 !+ 1d-6 * ( Uyy(:,:,:,3))
 Syy =  -1d-6 * D * (Uyy(:,:,:,1) + poisson * Uxx(:,:,:,1) ) * 6 / (thick)**2 !+ 1d-6 * ( Uxx(:,:,:,3))
 Sxy =  -1d-6 * D * (1 - poisson) * Uxy(:,:,:,1) * 6 / (thick)**2 !- 1d-6 * ( Uxy(:,:,:,3))
    
     

     write(*,*) " problem resolution OK "
     
     
   
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $ Desplazamiento U(x,y, t_0)  $ ", 2)

   

    call qplclr( U(0,:,:, 1), Nx+1, Ny+1)  
     
    call qplclr( U(0.25*Nt ,:,:, 1), Nx+1, Ny+1)    
    
    call qplclr( U(0.5*Nt ,:,:, 1), Nx+1, Ny+1)    
    
    call qplclr( U(0.75*Nt ,:,:, 1), Nx+1, Ny+1)    
                                   
    call qplclr( U(Nt,:,:, 1), Nx+1, Ny+1) 
     


  open(101, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Linear_Vibrations_Glazings\Linear_Vibrations_Glazings.plt')
write (101,*) 'VARIABLES = "X", "Y", "U0", "U1", "U2", "U3", "U4" '
write (101,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
  do j = 0, Ny 
  do  i=0 , Nx
 write (101,*) x(i), y(j), U(0,i,j, 1), U(0.25*Nt ,i,j, 1), U(0.5*Nt ,i,j, 1), U(0.75*Nt ,i,j, 1), U(Nt,i,j, 1)
  end do
 end do
close (101) 
     
  open(102, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Linear_Vibrations_Glazings\Sxx.plt')
write (102,*) 'VARIABLES = "X", "Y", "Sxx0", "Sxx1", "Sxx2", "Sxx3", "Sxx4" '
write (102,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
  do j = 0, Ny 
  do  i=0 , Nx
 write (102,*) x(i), y(j), Sxx(0,i,j), Sxx(0.25*Nt ,i,j), Sxx(0.5*Nt ,i,j), Sxx(0.75*Nt ,i,j), Sxx(Nt,i,j)
  end do
 end do
close (102) 

  open(104, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Linear_Vibrations_Glazings\Syy.plt')
write (104,*) 'VARIABLES = "X", "Y", "Syy0", "Syy1", "Syy2", "Syy3", "Syy4" '
write (104,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
  do j = 0, Ny 
  do  i=0 , Nx
 write (104,*) x(i), y(j), Syy(0,i,j), Syy(0.25*Nt ,i,j), Syy(0.5*Nt ,i,j), Syy(0.75*Nt ,i,j), Syy(Nt,i,j)
  end do
 end do
close (104) 

  open(103, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Linear_Vibrations_Glazings\Sxy.plt')
write (103,*) 'VARIABLES = "X", "Y", "Sxy0", "Sxy1", "Sxy2", "Sxy3", "Sxy4" '
write (103,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
  do j = 0, Ny 
  do  i=0 , Nx
 write (103,*) x(i), y(j), Sxy(0,i,j), Sxy(0.25*Nt ,i,j), Sxy(0.5*Nt ,i,j), Sxy(0.75*Nt ,i,j), Sxy(Nt,i,j)
  end do
 end do
close (103) 

contains 
!----------------------------------------------------------
function L1( x, y, t, u, ux, uy, uxx, uyy, uxy   ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real :: L1(size(u))
            
           
        
            L1(1)  =   u(2) 
            L1(2)  = -(uxx(3) + uyy(3))*(D_eq  /rho_eq )+ ((pi**4)* ((4/(Lx**2) + 4/(Ly**2))**2 )*D_eq  -rho_eq*2000*2000)*sin(pi*x/(Lx*0.5)) * sin(pi*y/(Ly*0.5))* cos(2000*t)
            L1(3)  =   uxx(2) + uyy(2)  
       
                 
end function 
!-------------------------------------------------------
function BC1(x,y, t, u, ux, uy) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BC1( size(u) ) 
  
         
  
       
  
        if (x==x0) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (x==xf) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (y==y0) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (y==yf) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)
        else
             write(*,*)  "Error in BC1 "
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif
      

end function


!----------------------------------------------------------
end subroutine


subroutine Linear_Plate_Pulse

       integer, parameter :: Nx = 20, Ny = 20, Nt = 550, Nv = 3 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv)   
       real ::  x0 , xf , y0 , yf
       real ::  t0 = 0 , tf  , pi
       integer :: i, j,m, Order = 4  
       real :: thick = 16d-3, rho_g = 2500
       real :: E = 72d6 , poisson = 0.22 
       real :: mu, lambda, D  
       
       pi = 4 * atan(1d0)
       tf = 1
       Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]   
               
       x0 = -1 ;  y0 = -1
       xf = 1  ;  yf =  1  
   
       x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
       y    = [ (y0 + (yf-y0)*j/Ny, j=0, Ny ) ]

      D = E * (thick**3) / (12 * (1- poisson**2)) 

  do i = 0, Nx
   do j= 0, Ny
   U(0,i,j,1) =  0
   U(0,i,j,2) = exp(-10*(x(i)**2+ y(j)**2 )   ) 
   U(0,i,j,3) =  0
   end do
  end do
   


     
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, N_variables = Nv, & 
                                           Differential_operator =  L1, Boundary_conditions   =  BC1, Solution = U ) 
     


contains 
!----------------------------------------------------------
function L1( x, y, t, u, ux, uy, uxx, uyy, uxy   ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real :: L1(size(u))
            
           
        
            L1(1)  =   u(2) 
            L1(2)  = -(uxx(3) + uyy(3))*D/(rho_g*thick)
            L1(3)  =  ( uxx(2) + uyy(2) )
       
                 
end function 
!-------------------------------------------------------
function BC1(x,y, t, u, ux, uy) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BC1( size(u) ) 
  
         
        if (x==x0) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (x==xf) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (y==y0) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (y==yf) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)
        else
             write(*,*)  "Error in BC1 "
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif
      

end function


!----------------------------------------------------------
end subroutine






subroutine Linear_Plate_Small_Pulse


       integer, parameter :: Nx = 20, Ny = 20, Nt = 550, Nv = 3 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv), Ux(0:Nt, 0:Nx, 0:Ny)
       real ::  Uxx(0:Nt, 0:Nx, 0:Ny), Uyy(0:Nt, 0:Nx, 0:Ny), Uxy(0:Nt, 0:Nx, 0:Ny)
       real ::  Mxx(0:Nt, 0:Nx, 0:Ny), Myy(0:Nt, 0:Nx, 0:Ny), Mxy(0:Nt, 0:Nx, 0:Ny)
       real ::  Sxx(0:Nt, 0:Nx, 0:Ny), Syy(0:Nt, 0:Nx, 0:Ny), Sxy(0:Nt, 0:Nx, 0:Ny)
       real :: Lx ,Ly       
       real ::  x0 , xf , y0 , yf
       real ::  t0 = 0 , tf  , pi
       integer :: i, j,m, Order = 4  
       real :: thick = 16d-3, rho_g = 2500
       real :: E = 72d6 , poisson = 0.22 
       real :: mu, lambda, D, eps  
       
       pi = 4 * atan(1d0)
       tf = 1
        


     ! x0 = -Lx*0.5 ;  y0 = -Ly*0.5
     ! xf =  Lx*0.5 ;  yf =  Ly*0.5
       x0 = -1 ;  y0 = -1
       xf = 1 ;  yf =  1  
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*j/Ny, j=0, Ny ) ]

      D = E * (thick**3) / (12 * (1- poisson**2)) 

      lambda = E/(rho_g * (xf-x0)*(yf-y0))
      eps =( thick**2/((xf-x0)*(yf-y0)))**0.5
      
  do i = 0, Nx
   do j= 0, Ny
  U(0,i,j,1) =  0
  U(0,i,j,2) = exp(-10*(x(i)**2+ y(j)**2 )   ) 
  U(0,i,j,3) =  0
   end do
  end do
   

   
     write(*,*) " example definition OK "
     write(*,*) " lambda = ", lambda, " eps = ", eps
     call scrmod('revers')                          
     call metafl('xwin')
     call disini()     
     CALL TEXMOD('ON')
     CALL TITLIN (" Introduced Initial value for U(1)  ", 2)
      call qplclr( U(0 ,:,:, 1), Nx+1, Ny+1)   

     
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, N_variables = Nv, & 
                                           Differential_operator =  L1, Boundary_conditions   =  BC1, Solution = U ) 
     


     write(*,*) " problem resolution OK "
     
  
  call Check_grid( "x", x, Order, size( U(0,:,0,1) ) ) 
  call Check_grid( "y", y, Order, size( U(0,0,:,1) ) )
     
    do i = 0, Nt
     
      call Derivative( ["x","y"], 1, 1,  U(i, 0:,0:,1),  Ux(i,0:,0:)  )
      call Derivative( ["x","y"], 1, 2,  U(i, 0:,0:,1), Uxx(i,0:,0:) )
                                                            

      call Derivative( ["x","y"], 2, 2,  U(i, 0:,0:,1), Uyy(i,0:,0:) )
                                                          
      call Derivative( ["x","y"], 2, 1, Ux(i, 0:,0:),Uxy(i,0:,0:))
      
    end do

 Mxx = - D * (Uxx(:,:,:) + poisson * Uyy(:,:,:) )
 Myy = - D * (Uyy(:,:,:) + poisson * Uxx(:,:,:) )
 Mxy = - D * (1 - poisson) * Uxy(:,:,:)

do m = 0,Nt
  do i = 0,Nx
     do j = 0,Ny
 Sxx(m,i,j) = ( Mxx(m,i,j) )* 6 / (thick)**2    !  1d-3 *
 Syy(m,i,j) = ( Myy(m,i,j) )* 6 / (thick)**2    !  1d-3 *
 Sxy(m,i,j) = ( Mxy(m,i,j) )* 6 / (thick)**2    !  1d-3 *
    end do
  end do
end do
     
     
     
     
   
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $ Desplazamiento U(x,y, t_0)  $ ", 2)

   

    call qplclr( U(0,:,:, 1), Nx+1, Ny+1)  
     
    call qplclr( U(0.25*Nt ,:,:, 1), Nx+1, Ny+1)    
    
    call qplclr( U(0.5*Nt ,:,:, 1), Nx+1, Ny+1)    
    
    call qplclr( U(0.75*Nt ,:,:, 1), Nx+1, Ny+1)    
                                   
    call qplclr( U(Nt,:,:, 1), Nx+1, Ny+1) 
     


open(87, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Linear_Vibrations_Plates_Small_Pulse\Desplazamiento.plt')
write (87,*) 'VARIABLES = "X", "Y", "U0", "U1", "U2", "U3", "U4" '
write (87,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (87,*) x(i), y(j), U(0,i,j, 1), U(0.25*Nt ,i,j, 1), U(0.5*Nt ,i,j, 1), U(0.75*Nt ,i,j, 1), U(Nt,i,j, 1)
  end do
  end do
close (87)  

open(87, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Linear_Vibrations_Plates_Small_Pulse\EsfuerzosFlexion.plt')
write (87,*) 'VARIABLES = "X", "Y", "Sxx0", "Sxx1", "Sxx2", "Sxx3", "Sxx4" '
write (87,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (87,*) x(i), y(j), Sxx(0,i,j), Sxx(0.25*Nt ,i,j), Sxx(0.5*Nt ,i,j), Sxx(0.75*Nt ,i,j), Sxx(Nt,i,j)
  end do
  end do
close (87)  

contains 
!----------------------------------------------------------
function L1( x, y, t, u, ux, uy, uxx, uyy, uxy   ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real :: L1(size(u))
            
           
        
            L1(1)  =   u(2) 
            L1(2)  = -(uxx(3) + uyy(3))*(D/(rho_g*thick))
            L1(3)  =  ( uxx(2) + uyy(2) )
       
                 
end function 
!-------------------------------------------------------
function BC1(x,y, t, u, ux, uy) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BC1( size(u) ) 
  
         
  
       
  
        if (x==x0) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (x==xf) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (y==y0) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (y==yf) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)
        else
             write(*,*)  "Error in BC1 "
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif
      

end function








!----------------------------------------------------------
end subroutine


subroutine Linear_Plate_Step

       integer, parameter :: Nx = 20, Ny = 20, Nt = 550, Nv = 3 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv), Ux(0:Nt, 0:Nx, 0:Ny)
       real ::  Uxx(0:Nt, 0:Nx, 0:Ny), Uyy(0:Nt, 0:Nx, 0:Ny), Uxy(0:Nt, 0:Nx, 0:Ny)
       real ::  Mxx(0:Nt, 0:Nx, 0:Ny), Myy(0:Nt, 0:Nx, 0:Ny), Mxy(0:Nt, 0:Nx, 0:Ny)
       real ::  Sxx(0:Nt, 0:Nx, 0:Ny), Syy(0:Nt, 0:Nx, 0:Ny), Sxy(0:Nt, 0:Nx, 0:Ny)
       real :: Lx ,Ly       
       real ::  x0 , xf , y0 , yf
       real ::  t0 = 0 , tf  , pi
       integer :: i, j,m, Order = 4  
       real :: thick = 16d-3, rho_g = 2500
       real :: E = 72d6 , poisson = 0.22 
       real :: mu, lambda, D  
       
       pi = 4 * atan(1d0)
       tf = 1
        


     ! x0 = -Lx*0.5 ;  y0 = -Ly*0.5
     ! xf =  Lx*0.5 ;  yf =  Ly*0.5
       x0 = -1 ;  y0 = -1
       xf = 1 ;  yf =  1  
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*j/Ny, j=0, Ny ) ]

      D = E * (thick**3) / (12 * (1- poisson**2)) 

  do i = 0, Nx
   do j= 0, Ny
  U(0,i,j,1) = exp(-10*(x(i)**2+ y(j)**2 )   )
  U(0,i,j,2) = 0 
  U(0,i,j,3) = -40*  U(0,i,j,1) + 40* (x(i)**2+ y(j)**2 )*  U(0,i,j,1)
   end do
  end do
   

   
     write(*,*) " example definition OK "
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini()     
     CALL TEXMOD('ON')
     CALL TITLIN (" Introduced Initial value for U(1)  ", 2)
      call qplclr( U(0 ,:,:, 1), Nx+1, Ny+1)   

     
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, N_variables = Nv, & 
                                           Differential_operator =  L1, Boundary_conditions   =  BC1, Solution = U ) 
     


     write(*,*) " problem resolution OK "
     
  
  call Check_grid( "x", x, Order, size( U(0,:,0,1) ) ) 
  call Check_grid( "y", y, Order, size( U(0,0,:,1) ) )
     
    do i = 0, Nt
     
      call Derivative( ["x","y"], 1, 1,  U(i, 0:,0:,1),  Ux(i,0:,0:)  )
      call Derivative( ["x","y"], 1, 2,  U(i, 0:,0:,1), Uxx(i,0:,0:) )
                                                            

      call Derivative( ["x","y"], 2, 2,  U(i, 0:,0:,1), Uyy(i,0:,0:) )
                                                          
      call Derivative( ["x","y"], 2, 1, Ux(i, 0:,0:),Uxy(i,0:,0:))
      
    end do

 Mxx = - D * (Uxx(:,:,:) + poisson * Uyy(:,:,:) )
 Myy = - D * (Uyy(:,:,:) + poisson * Uxx(:,:,:) )
 Mxy = - D * (1 - poisson) * Uxy(:,:,:)

do m = 0,Nt
  do i = 0,Nx
     do j = 0,Ny
 Sxx(m,i,j) = ( Mxx(m,i,j) )* 6 / (thick)**2    !  1d-3 *
 Syy(m,i,j) = ( Myy(m,i,j) )* 6 / (thick)**2    !  1d-3 *
 Sxy(m,i,j) = ( Mxy(m,i,j) )* 6 / (thick)**2    !  1d-3 *
    end do
  end do
end do
     
     
     
     
   
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $ Desplazamiento U(x,y, t_0)  $ ", 2)

   

    call qplclr( U(0,:,:, 1), Nx+1, Ny+1)  
     
    call qplclr( U(0.25*Nt ,:,:, 1), Nx+1, Ny+1)    
    
    call qplclr( U(0.5*Nt ,:,:, 1), Nx+1, Ny+1)    
    
    call qplclr( U(0.75*Nt ,:,:, 1), Nx+1, Ny+1)    
                                   
    call qplclr( U(Nt,:,:, 1), Nx+1, Ny+1) 
     


open(87, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Linear_Plates_Step\Desplazamiento.plt')
write (87,*) 'VARIABLES = "X", "Y", "U0", "U1", "U2", "U3", "U4" '
write (87,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (87,*) x(i), y(j), U(0,i,j, 1), U(0.25*Nt ,i,j, 1), U(0.5*Nt ,i,j, 1), U(0.75*Nt ,i,j, 1), U(Nt,i,j, 1)
  end do
  end do
close (87)  

open(87, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Linear_Plates_Step\EsfuerzosFlexion.plt')
write (87,*) 'VARIABLES = "X", "Y", "Sxx0", "Sxx1", "Sxx2", "Sxx3", "Sxx4" '
write (87,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (87,*) x(i), y(j), Sxx(0,i,j), Sxx(0.25*Nt ,i,j), Sxx(0.5*Nt ,i,j), Sxx(0.75*Nt ,i,j), Sxx(Nt,i,j)
  end do
  end do
close (87)  

open(86, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Linear_Plates_Step\DesplazamientoyEsfuerzosijk.plt')
write (86,*) 'VARIABLES = "X","Y", "T", "U", "Sxx"'
write (86,*) 'ZONE I=', Nx+1,', J=', Ny+1,', K=', Nt+1,',DATAPACKING=POINT'
 do  m= 0 , Nt
  do  i=0 , Nx
   do j = 0, Ny
 write (86,*)  x(i), y(j), Time(m),U(m,i,j,1), Sxx(m ,i,j)
   end do
  end do
end do
close (86)  

contains 
!----------------------------------------------------------
function L1( x, y, t, u, ux, uy, uxx, uyy, uxy   ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real :: L1(size(u))
            
           
        
            L1(1)  =   u(2) 
            L1(2)  = -(uxx(3) + uyy(3))*D/(rho_g*thick)
            L1(3)  =  ( uxx(2) + uyy(2) )
       
                 
end function 
!-------------------------------------------------------
function BC1(x,y, t, u, ux, uy) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BC1( size(u) ) 
  
         
  
       
  
        if (x==x0) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (x==xf) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (y==y0) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)

        else if (y==yf) then
                           BC1(1) = u(1)
                           BC1(2) = u(2)
                           BC1(3) = u(3)
        else
             write(*,*)  "Error in BC1 "
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif
      

end function


!----------------------------------------------------------
end subroutine
!-----------------------------------------------------------
subroutine Heat_equation_1D

       integer, parameter :: Nx = 30, Nt = 1000
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx)  
       
       real ::  x0 = 0, xf = 1
       real :: t0 = 0, tf = 1  
       integer :: i, Order = 6 
        
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x    = [ (x0 + (xf-x0)*j/Nx, j=0, Nx ) ]
     
     forall (i=0:Nx) U(0, i)  =  exp(-(5*x(i))**2 )
  
       
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, Order = Order, & 
                                           Differential_operator =  Heat_equation, & 
                                           Boundary_conditions   =  Heat_BC,  Solution = U ) 


     call qplot(x, U(4,:), Nx+1)           
     call qplot(x, U(0.25*Nt,:), Nx+1)    
     call qplot(x, U(0.5*Nt,:), Nx+1)       
     call qplot(x, U(Nt,:), Nx+1)                             


contains 
!----------------------------------------------------------
real function Heat_equation( x, t, u, ux, uxx) 
        real, intent(in) ::  x, t, u, ux, uxx
            
            real :: nu = 0.005
            
            Heat_equation =  nu * uxx
                 
end function 
!-------------------------------------------------------
real function Heat_BC(x, t, u, ux)
  real, intent(in) :: x, t, u, ux 
  
        if (x==x0) then
                            Burgers_BC = u - 1

        else if (x==xf) then
                            Burgers_BC = ux

        else
             write(*,*)  "Error in BC_Burgers"
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif

end function


end subroutine 
 


subroutine Heat_equation_2D

      integer, parameter :: Nx = 10, Ny = 10, Nt = 200
      real :: x(0:Nx), y(0:Ny), Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny) 
       
       real :: x0 = -1, xf = 1    
       real :: y0 = -1, yf = 1
       real :: t0 = 0, tf = 4
       integer :: i, j, Order = 8       
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*i/Ny, i=0, Ny ) ]
       
     !Condicion inicial
     forall (i=0:Nx,j=0:Ny) U(0, i, j)  =  exp(-(5*x(i))**2  -(5*y(j))**2)
    
     
     !Llamada a la función 
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, & 
                                           Differential_operator =  Heat_equation, & 
                                           Boundary_conditions   =  Heat_BC,  Solution = U ) 
     
  
     ! Graficas en cuatro instantes
     call qplcon( U(0, :, :), Nx+1, Ny+1, 20)

     call qplcon( U(Nt/4, :, :), Nx+1, Ny+1, 20)

     call qplcon( U(Nt/2, :, :), Nx+1, Ny+1, 20)

     call qplcon( U(Nt, :, :), Nx+1, Ny+1, 20)
                                   
 

contains

!----------------------------------------------------------
function Heat_equation( x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy ) result(F) 
        real,intent(in) :: x, y, t 
        real, intent(in) ::   U, Ux, Uy, Uxx, Uyy, Uxy
        real :: F 
        
        real :: nu = 0.02

        F =   nu * ( Uxx + Uyy )
        
   

end function
!-------------------------------------------------------
function Heat_BC( x, y, t, U, Ux, Uy ) result (BC) 

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



end module 
