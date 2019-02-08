module  API_Example_IVBP_and_BVP 


use IVBP_and_BVP 

implicit none 


contains 







!*************************************************************************************************
!*
!*************************************************************************************************
subroutine Test_IVBP_and_BVP

       integer, parameter :: Nx = 15, Ny = 15, Nt = 100, Nv1 = 3 , Nv2 = 1 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv1) , V(0:Nt, 0:Nx, 0:Ny, Nv2)
       
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
       real ::  t0 = 0, tf = 0.2 , pi
       integer :: i, j, Order = 4  
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
     x    = [ (x0 + (xf-x0)*i/Nx, i=0, Nx ) ]
     y    = [ (y0 + (yf-y0)*j/Ny, j=0, Ny ) ]

    pi = acos(-1d0)
   
    U(0, :,:,1) = 0 
    U(0, :,:,3) = 0 
    do i = 0, Nx
     do j= 0, Ny
         U(0,i,j,2) = exp( - 5 * x(i)**2 - 5 * y(j)**2 )  
     end do
    end do
    
     write(*,*) " example definition OK "
     write(*,*) " maxval U = ", maxval(U) 
     read(*,*)
     
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

     do i=0, Nt, Nt/10 
      call qplclr( U(i ,:,:, 2), Nx+1, Ny+1)   
     end do  
 
     

contains 
!----------------------------------------------------------
function L1( x, y, t, u, ux, uy, uxx, uyy, uxy , v, vx, vy, vxx, vyy, vxy  ) 
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
        real :: L1(size(u))
            
                 
         real :: w, dwdt, d2wdx2
            
            w = u(1);   dwdt = u(2); d2wdx2 = u(3) 
            
        
            L1(1)  = dwdt 
            L1(2)  = - uxx(3)  - uyy(3)
            L1(3)  =   uxx(2)  + uyy(2)
                   
end function 
!-------------------------------------------------------
function BC1(x,y, t, u, ux, uy) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BC1( size(u) ) 
  
      real ::  w, dwdt, d2wdx2 
            
       w = u(1)
       dwdt = u(2) 
       d2wdx2 = u(3) 
       
  
        if (x==x0) then
                           BC1(1) = w
                           BC1(2) = dwdt
                           BC1(3) = d2wdx2 

        else if (x==xf) then
                            BC1(1) = w
                            BC1(2) = dwdt
                            BC1(3) = d2wdx2
        else if (y==y0) then
                            BC1(1) = w
                            BC1(2) = dwdt
                            BC1(3) = d2wdx2

        else if (y==yf) then
                             BC1(1) = w
                             BC1(2) = dwdt
                             BC1(3) = d2wdx2
        else
             write(*,*)  "Error in BC2D_waves"
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
            
            
         
          L2(1)  = vxx(1) + vyy(1) - uxx(1) - uyy(1)
        
               
                 
end function 
!-------------------------------------------------------
function BC2(x,y, t, v, vx, vy)  
  real, intent(in) :: x, y,  t, v(:), vx(:), vy(:) 
  real :: BC2( size(v) ) 
  
        
       
  
        if (x==x0) then
                           BC2(1) = v(1)
                         

        else if (x==xf) then
                           BC2(1) = v(1)
                        

        else if (y==y0) then
                           BC2(1) = v(1)
                        

        else if (y==yf) then
                           BC2(1) = v(1)
                       

        else
             write(*,*)  "Error in BC2 "
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif
      

end function


end subroutine 



end module 
