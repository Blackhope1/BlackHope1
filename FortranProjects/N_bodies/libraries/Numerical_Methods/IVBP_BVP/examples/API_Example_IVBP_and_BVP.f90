module  API_Example_IVBP_and_BVP 


use IVBP_and_BVP 

implicit none 


contains 







!*************************************************************************************************
!*
!*************************************************************************************************
subroutine Test_IVBP_and_BVP

       integer, parameter :: Nx = 15, Ny = 15, Nt = 1, Nv1 = 1 , Nv2 = 2 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv1) , V (0:Nt, 0:Nx, 0:Ny, Nv2)
       
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
       real ::  t0 = 0, tf = 0.5 , pi
       integer :: i, j, Order = 4  
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*j/Ny, j=0, Ny ) ]

    pi = acos(-1d0)
   
    do i = 0, Nx
     do j= 0, Ny
         U(0,i,j,1) = exp( - x(i)**2 - y(j)**2 )  
     end do
    end do
    
     write(*,*) " example definition OK "
     
     !call scrmod('revers')                          
     !call metafl('xwin')
     !call disini()     

   !  call qplclr( U(0 ,:,:, 1), Nx+1, Ny+1)   
    
     
     call IVBP_and_BVP_Problem( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, N_u = Nv1, N_v = Nv2, & 
                                Differential_operator_u =  L1, Differential_operator_v =  L2, & 
                                Boundary_conditions_u   =  BC1, Boundary_conditions_v =  BC2,  Ut = U,  Vt = V ) 
     
     write(*,*) " problem resolution OK "
     
     
   
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $  \frac{\partial{^2 u}}{\partial{t^2}}   = \nabla ^4 u  $ ", 2)
     CALL TITLIN (" $ u(-1) = 0, \quad u(1) = 0 $  ", 4)

     call qplclr( U(0 ,:,:, 1), Nx+1, Ny+1)   
 
     call qplclr( U(0.25*Nt ,:,:, 1), Nx+1, Ny+1)    
     
     call qplclr( U(0.5*Nt ,:,:, 1), Nx+1, Ny+1)    
    
     call qplclr( U(0.75*Nt ,:,:, 1), Nx+1, Ny+1)    
                                    
     call qplclr( U(Nt,:,:, 1), Nx+1, Ny+1) 
   

contains 
!----------------------------------------------------------
function L1( x, y, t, u, ux, uy, uxx, uyy, uxy , v, vx, vy, vxx, vyy, vxy  ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
        real :: L1(size(u))
            
            real :: w, wx, wy, wxx, wyy, v1, v2  
        
            w = u(1); wx = ux(1); wy = uy(1); wxx = uxx(1); wyy = uyy(1);
            v1 = v(1); v2 = v(2)
            
            
            L1(1) = - dot_product( [v1, v2 ] , [ wx, wy ] )  + wxx + wyy 
        
                   
end function 
!-------------------------------------------------------
function BC1(x,y, t, u, ux, uy) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BC1( size(u) ) 
  
       real ::  w 
            
       w = u(1)
  
        if (x==x0) then
                           BC1(1) = w

        else if (x==xf) then
                            BC1(1) = w
                          
        else if (y==y0) then
                           
                            BC1(1) = w 

        else if (y==yf) then
                        
                            BC1(1) = w 

        else
             write(*,*)  "Error in BC1 "
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif
      

end function

!----------------------------------------------------------
function L2( x, y, t, v, vx, vy, vxx, vyy, vxy , u, ux, uy, uxx, uyy, uxy  ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
        real :: L2(size(v))
            
            
          real :: w, v1x, v1y, v2x, v2y 
          
          w = u(1); 
          v1x = vx(1); v1y = vy(1) 
          v2x = vx(2); v2y = vy(2) 
            
          L2(1)  = v1x + v2y 
        
          L2(2)  = v2x - v1y + w  
            
       
                 
end function 
!-------------------------------------------------------
function BC2(x,y, t, v, vx, vy)  
  real, intent(in) :: x, y,  t, v(:), vx(:), vy(:) 
  real :: BC2( size(v) ) 
  
        
       
  
        if (x==x0) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)

        else if (x==xf) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)

        else if (y==y0) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)


        else if (y==yf) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)

        else
             write(*,*)  "Error in BC2 "
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif
      

end function


end subroutine 
subroutine Linear_Plate_Vibrations_IVBP


       integer, parameter :: Nx = 24, Ny = 24, Nt  = 800, Nv1 = 3 , Nv2 = 1 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv1) , V (0:Nt, 0:Nx, 0:Ny, Nv2)
       
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
       real ::  t0 = 0 , tf  , pi
       integer :: i, j, m, Order = 4 
       
       pi = 4 * atan(1d0)
       tf = 1
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*j/Ny, j=0, Ny ) ]

   !  D = E * (thick**3) / (12 * (1- poisson**2)) 

  do i = 0, Nx
   do j= 0, Ny
  U(0,i,j,1) =  sin(pi*x(i))*sin(pi*y(j))
  U(0,i,j,2) =  0
  U(0,i,j,3) =  -2*(pi**2)*sin(pi*x(i))*sin(pi*y(j))
   end do
  end do
   
   

   
     write(*,*) " example definition OK "
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini()     
     CALL TEXMOD('ON')
     CALL TITLIN (" Introduced Initial value for U(1)  ", 2)
      call qplclr( U(0 ,:,:, 1), Nx+1, Ny+1)   

     
     call IVBP_and_BVP_Problem( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, N_u = Nv1, N_v = Nv2, & 
                                           Differential_operator_u =  L1, Differential_operator_v =  L2, & 
                                           Boundary_conditions_u   =  BC1, Boundary_conditions_v =  BC2, Ut = U,  Vt = V ) 
     


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
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $ Desplazamiento V(x,y, t_0)  $ ", 2)

    call qplclr( V(0,:,:, 1), Nx+1, Ny+1)  
     
    call qplclr( V(0.25*Nt ,:,:, 1), Nx+1, Ny+1)    
    
    call qplclr( V(0.5*Nt ,:,:, 1), Nx+1, Ny+1)    
    
    call qplclr( V(0.75*Nt ,:,:, 1), Nx+1, Ny+1)    
                                   
    call qplclr( V(Nt,:,:, 1), Nx+1, Ny+1) 
    

  open(101, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto Lopez\figuras\Graficas\Dynamic_Behaviour\Linear_Vibrations_Plates_IVBP_BVP\Linear_Vibrations_Plates.plt')
write (101,*) 'VARIABLES = "X", "Y", "U0", "U1", "U2", "U3", "U4" '
write (101,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
  do j = 0, Ny 
  do  i=0 , Nx
 write (101,*) x(i), y(j), 1d3*U(0,i,j, 1), 1d3*U(0.25*Nt ,i,j, 1), 1d3*U(0.5*Nt ,i,j, 1), 1d3*U(0.75*Nt ,i,j, 1), 1d3*U(Nt,i,j, 1)
  end do
 end do
close (101) 
  
   open(86, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Linear_Vibrations_Plates_IVBP_BVP\DesplazamientoyEsfuerzosijk.plt')
write (86,*) 'VARIABLES = "X","Y", "T", "U" '
write (86,*) 'ZONE I=', Nx+1,', J=', Ny+1,', K=', Nt+1,',DATAPACKING=POINT'
 do  m= 0 , Nt
  do  i=0 , Nx
   do j = 0, Ny
 write (86,*)  x(i), y(j), Time(m),U(m,i,j,1)
   end do
  end do
end do
close (86)    

contains 
!----------------------------------------------------------
function L1( x, y, t, u, ux, uy, uxx, uyy, uxy , v, vx, vy, vxx, vyy, vxy  ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
        real :: L1(size(u))
            
           
        
            L1(1)  =   u(2) 
            L1(2)  = -(uxx(3)+uyy(3))  + (4*pi**4-1) *sin(pi*x) * sin(pi*y)* cos(t)
            L1(3)  = uxx(2) + uyy(2)
       
                 
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
function L2( x, y, t, v, vx, vy, vxx, vyy, vxy , u, ux, uy, uxx, uyy, uxy  ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
        real :: L2(size(v))
            
               
            L2(1)  = v(1)    - u(1)

            L2(2)  = v(2)    - u(2)
       
                 
end function 
!-------------------------------------------------------
function BC2(x,y, t, v, vx, vy)  
  real, intent(in) :: x, y,  t, v(:), vx(:), vy(:) 
  real :: BC2( size(v) ) 
  
        
       
  
        if (x==x0) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)


        else if (x==xf) then
                           BC2(1) = v(1)
                          BC2(2) = v(2)
        else if (y==y0) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)

        else if (y==yf) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)
        else
             write(*,*)  "Error in BC2 "
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif
      

end function
end subroutine













subroutine Non_Linear_Plate_Vibrations

     integer, parameter :: Nx = 10, Ny = 10, Nt = 200, Nv1 = 3 , Nv2 = 2 
     real ::  x(0:Nx), y(0:Ny)
     real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv1) , V (0:Nt, 0:Nx, 0:Ny, Nv2)
     real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
     real ::  t0 = 0, tf = 1, pi = 4*atan(1d0), thickness = 16d-3, D 
     real ::  poisson = 0.22, E = 72d3, rho_g=2.500
     integer :: i, j, m, Order = 4 
      
       
    Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
    x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
    y    = [ (y0 + (yf-y0)*j/Ny, j=0, Ny ) ]

    D = E * thickness**3 / ( 12*(1-poisson**2)  )
       
 
   U(0,:,:,1) =  1d-3 * Tensor_product( sin(pi*x), sin(pi*y) )  
   U(0,:,:,2) =  0
   U(0,:,:,3) =  -2*1d-3 * pi**2 * Tensor_product( sin(pi*x), sin(pi*y) ) 
 
     
   call IVBP_and_BVP_Problem( Time_Domain = Time,                  & 
                              x_nodes = x, y_nodes = y,            & 
                              Order = Order, N_u = Nv1, N_v = Nv2, & 
                              Differential_operator_u = Lu,        & 
                              Differential_operator_v = Lv,        & 
                              Boundary_conditions_u   = BCu,       & 
                              Boundary_conditions_v   = BCv,       &    
                              Ut = U,  Vt = V ) 
     
   call scrmod("reverse")  
   call qplcon( U(0,:, :, 1),   Nx+1, Ny+1, 20)
   call qplcon( U(Nt, :, :, 1), Nx+1, Ny+1, 20)
   
contains 

function Lu( x, y, t, u, ux, uy, uxx, uyy, uxy , v, vx, vy, vxx, vyy, vxy  ) 
 real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
 real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
 real :: Lu(size(u))
            
  Lu(1)  =   u(2) 
  
  Lu(2)  =    (-(uxx(3) + uyy(3)) * D / (rho_g*thickness) +               & 
      
           +  ( vyy(1)*uxx(1) + vxx(1)*uyy(1) - 2* vxy(1)*uxy(1))/(rho_g) & 
      
            + 1d-3*(4*(D/(rho_g*thickness))*pi**4-4*pi**2)                & 
                    * cos(pi*x/2) * cos(pi*y)* cos(2*pi*t))
  
  Lu(3)  =   uxx(2) + uyy(2) 
       
                 
end function 

function BCu(x,y, t, u, ux, uy) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BCu( size(u) ) 
       
        if (x==x0) then
                           BCu(1) = u(1)
                           BCu(2) = u(2)
                           BCu(3) = u(3)

        else if (x==xf) then
                           BCu(1) = u(1)
                           BCu(2) = u(2)
                           BCu(3) = u(3)

        else if (y==y0) then
                           BCu(1) = u(1)
                           BCu(2) = u(2)
                           BCu(3) = u(3)

        else if (y==yf) then
                           BCu(1) = u(1)
                           BCu(2) = u(2)
                           BCu(3) = u(3)
        else
             write(*,*)  "Error in BC1 "
             stop 
        endif

end function

function Lv( x, y, t, v, vx, vy, vxx, vyy, vxy , u, ux, uy, uxx, uyy, uxy  ) 
 real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
 real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
 real :: Lv(size(v))
            
               
 Lv(1)  = vxx(1) + vyy(1) - v(2)  
 
 Lv(2)  =  vxx(2) + vyy(2)           & 
     
         + 0.5*E*( uyy(1)*uxx(1) + uxx(1)*uyy(1) - 2* uxy(1)*uxy(1) )
       
                 
end function 


function BCv(x,y, t, v, vx, vy)  
  real, intent(in) :: x, y,  t, v(:), vx(:), vy(:) 
  real :: BCv( size(v) ) 
    
  
        if (x==x0) then
                           BCv(1) = v(1)
                           BCv(2) = v(2)

        else if (x==xf) then
                           BCv(1) = v(1)
                           BCv(2) = v(2)

        else if (y==y0) then
                           BCv(1) = v(1)
                           BCv(2) = v(2)


        else if (y==yf) then
                           BCv(1) = v(1)
                           BCv(2) = v(2)


        else
             write(*,*)  "Error in BC2 "
             stop 
        endif

end function
end subroutine

subroutine Non_Linear_Plate_Pulse

       integer, parameter :: Nx = 20, Ny = 20, Nt = 550, Nv1 = 3 , Nv2 = 2 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv1) , V (0:Nt, 0:Nx, 0:Ny, Nv2), Ux(0:Nt, 0:Nx, 0:Ny)
       real ::  Uxx(0:Nt, 0:Nx, 0:Ny), Uyy(0:Nt, 0:Nx, 0:Ny), Uxy(0:Nt, 0:Nx, 0:Ny)
       real ::  Mxx(0:Nt, 0:Nx, 0:Ny), Myy(0:Nt, 0:Nx, 0:Ny), Mxy(0:Nt, 0:Nx, 0:Ny)
       real ::  Sxx(0:Nt, 0:Nx, 0:Ny), Syy(0:Nt, 0:Nx, 0:Ny), Sxy(0:Nt, 0:Nx, 0:Ny)
       
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
       real ::  t0 = 0, tf  , pi , thick = 16d-3, D, poisson = 0.22, E = 72d6, rho_g=2500
       integer :: i, j,m, Order = 4 
       
       pi = 4*atan(1d0)
       tf = 1
       
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*j/Ny, j=0, Ny ) ]

   
       D = E* thick**3 / ( 12*(1-poisson**2)  )

  do i = 0, Nx
   do j= 0, Ny
  U(0,i,j,1) =  0
  U(0,i,j,2) = exp(-10*(x(i)**2+ y(j)**2 )   ) 
  U(0,i,j,3) =  0
   end do
  end do

   
     write(*,*) " example definition OK "
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini()     
     CALL TEXMOD('ON')
     CALL TITLIN (" Introduced Initial value for U(1)  ", 2)
      call qplclr( U(0 ,:,:, 1), Nx+1, Ny+1)   

     
     call IVBP_and_BVP_Problem( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, N_u = Nv1, N_v = Nv2, & 
                                           Differential_operator_u =  L1, Differential_operator_v =  L2, & 
                                           Boundary_conditions_u   =  BC1, Boundary_conditions_v =  BC2, Ut = U,  Vt = V ) 
     


     write(*,*) " problem resolution OK "
     
     
   
  call Grid_Initialization( "nonuniform", "x", Order, x ) 
  call Grid_Initialization( "nonuniform", "y", Order, y )
 ! call Check_grid( "x", x, Order, size( U(0,:,0,1) ) ) 
 ! call Check_grid( "y", y, Order, size( U(0,0,:,1) ) )
     
    do i = 0, Nt
     
      call Derivative( ["x","y"], 1, 1,  U(i, 0:,0:,1),  Ux(i,0:,0:)  )
      call Derivative( ["x","y"], 1, 2,  U(i, 0:,0:,1), Uxx(i,0:,0:) )
                                                            

      call Derivative( ["x","y"], 2, 2,  U(i, 0:,0:,1), Uyy(i,0:,0:) )
                                                          
      call Derivative( ["x","y"], 2, 1, Ux(i, 0:,0:),Uxy(i,0:,0:))
      
    end do

 Mxx = -1d-3 *D * (Uxx(:,:,:) + poisson * Uyy(:,:,:) )
 Myy = -1d-3 *D * (Uyy(:,:,:) + poisson * Uxx(:,:,:) )
 Mxy = -1d-3 *D * (1 - poisson) * Uxy(:,:,:)

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
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $ Desplazamiento V(x,y, t_0)  $ ", 2)

    call qplclr( V(0,:,:, 2), Nx+1, Ny+1)  
     
    call qplclr( V(0.25*Nt ,:,:, 2), Nx+1, Ny+1)    
    
    call qplclr( V(0.5*Nt ,:,:, 2), Nx+1, Ny+1)    
    
    call qplclr( V(0.75*Nt ,:,:, 2), Nx+1, Ny+1)    
                              
    call qplclr( V(Nt,:,:, 2), Nx+1, Ny+1) 


open(10, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Non_Linear_Vibrations_Plates_Pulse\Desplazamiento.plt')
write (10,*) 'VARIABLES = "X", "Y", "U0", "U1", "U2", "U3", "U4" '
write (10,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (10,*) x(i), y(j), U(0,i,j, 1), U(0.25*Nt ,i,j, 1), U(0.5*Nt ,i,j, 1), U(0.75*Nt ,i,j, 1), U(Nt,i,j, 1)
  end do
  end do
close (10)  

open(12, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Non_Linear_Vibrations_Plates_Pulse\Esfuerzos.plt')
write (12,*) 'VARIABLES = "X", "Y", "V0", "V1", "V2", "V3", "V4" '
write (12,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (12,*) x(i), y(j), 1d-3*V(0,i,j, 2), 1d-3*V(0.25*Nt ,i,j, 2), 1d-3*V(0.5*Nt ,i,j, 2), 1d-3*V(0.75*Nt ,i,j, 2), 1d-3*V(Nt,i,j, 2)
  end do
 end do
close (12)  

open(17, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Non_Linear_Vibrations_Plates_Pulse\EsfuerzosFlexion.plt')
write (17,*) 'VARIABLES = "X", "Y", "V0", "V1", "V2", "V3", "V4" '
write (17,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (17,*) x(i), y(j), Sxx(0,i,j), Sxx(0.25*Nt ,i,j), Sxx(0.5*Nt ,i,j), Sxx(0.75*Nt ,i,j), Sxx(Nt,i,j)
  end do
 end do
close (17) 

open(86, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Non_Linear_Vibrations_Plates_Pulse\DesplazamientoyEsfuerzosijk.plt')
write (86,*) 'VARIABLES = "X","Y", "T", "U", "Sxx", "Phixx"'
write (86,*) 'ZONE I=', Nx+1,', J=', Ny+1,', K=', Nt+1,',DATAPACKING=POINT'
 do  m= 0 , Nt
  do  i=0 , Nx
   do j = 0, Ny
 write (86,*)  x(i), y(j), Time(m),U(m,i,j,1), Sxx(m ,i,j),V(m,i,j,2)
   end do
  end do
end do
close (86)  
contains 
!----------------------------------------------------------
function L1( x, y, t, u, ux, uy, uxx, uyy, uxy , v, vx, vy, vxx, vyy, vxy  ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
        real :: L1(size(u))
            
           
        
            L1(1)  =   u(2) 
            L1(2)  = -(uxx(3) + uyy(3))*D/(rho_g*thick) + (vyy(1)*uxx(1) + vxx(1)*uyy(1) - 2* vxy(1)*uxy(1)) /(rho_g)
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
function L2( x, y, t, v, vx, vy, vxx, vyy, vxy , u, ux, uy, uxx, uyy, uxy  ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
        real :: L2(size(v))
            
               
            L2(1)  = vxx(1) + vyy(1) - v(2)  
            L2(2)  = vxx(2) + vyy(2) + E*0.5*1d-6*( uyy(1)*uxx(1) + uxx(1)*uyy(1) - 2* uxy(1)*uxy(1) )
            
       
                 
end function 
!-------------------------------------------------------
function BC2(x,y, t, v, vx, vy)  
  real, intent(in) :: x, y,  t, v(:), vx(:), vy(:) 
  real :: BC2( size(v) ) 
  
        
       
  
        if (x==x0) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)

        else if (x==xf) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)

        else if (y==y0) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)


        else if (y==yf) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)


        else
             write(*,*)  "Error in BC2 "
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif
      

end function
end subroutine



subroutine Non_Linear_Plate_Small_Pulse

       integer, parameter :: Nx = 20, Ny = 20, Nt = 550, Nv1 = 3 , Nv2 = 2 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv1) , V (0:Nt, 0:Nx, 0:Ny, Nv2), Ux(0:Nt, 0:Nx, 0:Ny)
       real ::  Uxx(0:Nt, 0:Nx, 0:Ny), Uyy(0:Nt, 0:Nx, 0:Ny), Uxy(0:Nt, 0:Nx, 0:Ny)
       real ::  Mxx(0:Nt, 0:Nx, 0:Ny), Myy(0:Nt, 0:Nx, 0:Ny), Mxy(0:Nt, 0:Nx, 0:Ny)
       real ::  Sxx(0:Nt, 0:Nx, 0:Ny), Syy(0:Nt, 0:Nx, 0:Ny), Sxy(0:Nt, 0:Nx, 0:Ny)
       
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
       real ::  t0 = 0, tf  , pi , thick = 16d-3, D, poisson = 0.22, E = 72d6, rho_g = 2500
       integer :: i, j,m, Order = 4 
       
       pi = 4*atan(1d0)
       tf = 1
       
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*j/Ny, j=0, Ny ) ]

   
       D = E* thick**3 / ( 12*(1-poisson**2)  )

  do i = 0, Nx
   do j= 0, Ny
  U(0,i,j,1) =  0
  U(0,i,j,2) = 10*exp(-10*(x(i)**2+ y(j)**2 )   )    
  U(0,i,j,3) =  0
   end do
  end do

   
     write(*,*) " example definition OK "
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini()     
     CALL TEXMOD('ON')
     CALL TITLIN (" Introduced Initial value for U(1)  ", 2)
      call qplclr( U(0 ,:,:, 1), Nx+1, Ny+1)   

     
     call IVBP_and_BVP_Problem( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, N_u = Nv1, N_v = Nv2, & 
                                           Differential_operator_u =  L1, Differential_operator_v =  L2, & 
                                           Boundary_conditions_u   =  BC1, Boundary_conditions_v =  BC2, Ut = U,  Vt = V ) 
     


     write(*,*) " problem resolution OK "
     
     
   
  call Grid_Initialization( "nonuniform", "x", Order, x ) 
  call Grid_Initialization( "nonuniform", "y", Order, y ) 
  
 ! call Check_grid( "x", x, Order, size( U(0,:,0,1) ) ) 
 ! call Check_grid( "y", y, Order, size( U(0,0,:,1) ) )
     
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
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $ Desplazamiento V(x,y, t_0)  $ ", 2)

    call qplclr( V(0,:,:, 2), Nx+1, Ny+1)  
     
    call qplclr( V(0.25*Nt ,:,:, 2), Nx+1, Ny+1)    
    
    call qplclr( V(0.5*Nt ,:,:, 2), Nx+1, Ny+1)    
    
    call qplclr( V(0.75*Nt ,:,:, 2), Nx+1, Ny+1)    
                                   
    call qplclr( V(Nt,:,:, 2), Nx+1, Ny+1) 


open(10, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Non_Linear_Plates_Small_Pulse\Desplazamiento.plt')
write (10,*) 'VARIABLES = "X", "Y", "U0", "U1", "U2", "U3", "U4" '
write (10,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (10,*) x(i), y(j), U(0,i,j, 1), U(0.25*Nt ,i,j, 1), U(0.5*Nt ,i,j, 1), U(0.75*Nt ,i,j, 1), U(Nt,i,j, 1)
  end do
  end do
close (10)  

open(12, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Non_Linear_Plates_Small_Pulse\Esfuerzos.plt')
write (12,*) 'VARIABLES = "X", "Y", "V0", "V1", "V2", "V3", "V4" '
write (12,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (12,*) x(i), y(j), 1d-3*V(0,i,j, 2), 1d-3*V(0.25*Nt ,i,j, 2), 1d-3*V(0.5*Nt ,i,j, 2), 1d-3*V(0.75*Nt ,i,j, 2), 1d-3*V(Nt,i,j, 2)
  end do
 end do
close (12)  

open(17, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Non_Linear_Plates_Small_Pulse\EsfuerzosFlexion.plt')
write (17,*) 'VARIABLES = "X", "Y", "V0", "V1", "V2", "V3", "V4" '
write (17,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (17,*) x(i), y(j), Sxx(0,i,j), Sxx(0.25*Nt ,i,j), Sxx(0.5*Nt ,i,j), Sxx(0.75*Nt ,i,j), Sxx(Nt,i,j)
  end do
 end do
close (17) 

contains 
!----------------------------------------------------------
function L1( x, y, t, u, ux, uy, uxx, uyy, uxy , v, vx, vy, vxx, vyy, vxy  ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
        real :: L1(size(u))
            
           
        
            L1(1)  =   u(2) 
            L1(2)  = -(uxx(3) + uyy(3))*(D/(rho_g*thick)) + (vyy(1)*uxx(1) + vxx(1)*uyy(1) - 2* vxy(1)*uxy(1)) /(rho_g)
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
function L2( x, y, t, v, vx, vy, vxx, vyy, vxy , u, ux, uy, uxx, uyy, uxy  ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
        real :: L2(size(v))
            
               
            L2(1)  = vxx(1) + vyy(1) - v(2)  
            L2(2)  = vxx(2) + vyy(2) + 0.5*E*1d-3*( uyy(1)*uxx(1) + uxx(1)*uyy(1) - 2* uxy(1)*uxy(1) )
            
       
                 
end function 
!-------------------------------------------------------
function BC2(x,y, t, v, vx, vy)  
  real, intent(in) :: x, y,  t, v(:), vx(:), vy(:) 
  real :: BC2( size(v) ) 
  
        
       
  
        if (x==x0) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)

        else if (x==xf) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)

        else if (y==y0) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)


        else if (y==yf) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)


        else
             write(*,*)  "Error in BC2 "
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif
      

end function
end subroutine

subroutine Non_Linear_Plate_Step

       integer, parameter :: Nx = 20, Ny = 20, Nt = 550, Nv1 = 3 , Nv2 = 2 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv1) , V (0:Nt, 0:Nx, 0:Ny, Nv2), Ux(0:Nt, 0:Nx, 0:Ny)
       real ::  Uxx(0:Nt, 0:Nx, 0:Ny), Uyy(0:Nt, 0:Nx, 0:Ny), Uxy(0:Nt, 0:Nx, 0:Ny)
       real ::  Mxx(0:Nt, 0:Nx, 0:Ny), Myy(0:Nt, 0:Nx, 0:Ny), Mxy(0:Nt, 0:Nx, 0:Ny)
       real ::  Sxx(0:Nt, 0:Nx, 0:Ny), Syy(0:Nt, 0:Nx, 0:Ny), Sxy(0:Nt, 0:Nx, 0:Ny)
       
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
       real ::  t0 = 0, tf  , pi , thick = 16d-3, D, poisson = 0.22, E = 72d6, rho_g = 2500
       integer :: i, j,m, Order = 4 
       
       pi = 4*atan(1d0)
       tf = 1
       
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*j/Ny, j=0, Ny ) ]

   
       D = E* thick**3 / ( 12*(1-poisson**2)  )

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

     
     call IVBP_and_BVP_Problem( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, N_u = Nv1, N_v = Nv2, & 
                                           Differential_operator_u =  L1, Differential_operator_v =  L2, & 
                                           Boundary_conditions_u   =  BC1, Boundary_conditions_v =  BC2, Ut = U,  Vt = V ) 
     


     write(*,*) " problem resolution OK "
     
     
  call Grid_Initialization( "nonuniform", "x", Order, x ) 
  call Grid_Initialization( "nonuniform", "y", Order, y )  
       
!  call Check_grid( "x", x, Order, size( U(0,:,0,1) ) ) 
!  call Check_grid( "y", y, Order, size( U(0,0,:,1) ) )
     
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
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $ Desplazamiento V(x,y, t_0)  $ ", 2)

    call qplclr( V(0,:,:, 2), Nx+1, Ny+1)  
     
    call qplclr( V(0.25*Nt ,:,:, 2), Nx+1, Ny+1)    
    
    call qplclr( V(0.5*Nt ,:,:, 2), Nx+1, Ny+1)    
    
    call qplclr( V(0.75*Nt ,:,:, 2), Nx+1, Ny+1)    
                                   
    call qplclr( V(Nt,:,:, 2), Nx+1, Ny+1) 


open(10, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Non_Linear_Plates_Step\Desplazamiento.plt')
write (10,*) 'VARIABLES = "X", "Y", "U0", "U1", "U2", "U3", "U4" '
write (10,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (10,*) x(i), y(j), U(0,i,j, 1), U(0.25*Nt ,i,j, 1), U(0.5*Nt ,i,j, 1), U(0.75*Nt ,i,j, 1), U(Nt,i,j, 1)
  end do
  end do
close (10)  

open(12, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Non_Linear_Plates_Step\Esfuerzos.plt')
write (12,*) 'VARIABLES = "X", "Y", "V0", "V1", "V2", "V3", "V4" '
write (12,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (12,*) x(i), y(j), 1d-3*V(0,i,j, 2), 1d-3*V(0.25*Nt ,i,j, 2), 1d-3*V(0.5*Nt ,i,j, 2), 1d-3*V(0.75*Nt ,i,j, 2), 1d-3*V(Nt,i,j, 2)
  end do
 end do
close (12)  

open(17, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Non_Linear_Plates_Step\EsfuerzosFlexion.plt')
write (17,*) 'VARIABLES = "X", "Y", "V0", "V1", "V2", "V3", "V4" '
write (17,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (17,*) x(i), y(j), Sxx(0,i,j), Sxx(0.25*Nt ,i,j), Sxx(0.5*Nt ,i,j), Sxx(0.75*Nt ,i,j), Sxx(Nt,i,j)
  end do
 end do
close (17) 

open(86, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Dynamic_Behaviour\Non_Linear_Plates_Step\DesplazamientoyEsfuerzosijk.plt')
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
function L1( x, y, t, u, ux, uy, uxx, uyy, uxy , v, vx, vy, vxx, vyy, vxy  ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
        real :: L1(size(u))
            
           
        
            L1(1)  =   u(2) 
            L1(2)  = -(uxx(3) + uyy(3))*(D/(rho_g*thick)) + (vyy(1)*uxx(1) + vxx(1)*uyy(1) - 2* vxy(1)*uxy(1)) /(rho_g)
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
function L2( x, y, t, v, vx, vy, vxx, vyy, vxy , u, ux, uy, uxx, uyy, uxy  ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
        real :: L2(size(v))
            
               
            L2(1)  = vxx(1) + vyy(1) - v(2)  
            L2(2)  = vxx(2) + vyy(2) + 0.5*E*1d-3*( uyy(1)*uxx(1) + uxx(1)*uyy(1) - 2* uxy(1)*uxy(1) )
            
       
                 
end function 
!-------------------------------------------------------
function BC2(x,y, t, v, vx, vy)  
  real, intent(in) :: x, y,  t, v(:), vx(:), vy(:) 
  real :: BC2( size(v) ) 
  
        
       
  
        if (x==x0) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)

        else if (x==xf) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)

        else if (y==y0) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)


        else if (y==yf) then
                           BC2(1) = v(1)
                           BC2(2) = v(2)


        else
             write(*,*)  "Error in BC2 "
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif
      

end function
end subroutine


end module 
