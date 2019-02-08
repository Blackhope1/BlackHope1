module test_2D 
    
    use Linear_systems
    use Initial_Value_Boundary_Problem
    implicit none 
    
 contains    
!*************************************************************************************************
!*
!*************************************************************************************************
subroutine Test_IVBP1D

       integer, parameter :: Nx = 20, Nt = 1000
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx)  
       
       real ::  x0 = 0, xf = 1
       real :: t0 = 0, tf = 0.2  
       integer :: i, Order = 20 
       
        
       
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x(0) = x0; x(Nx) = xf 
     call Grid_Initialization( "nonuniform", "x", Order, x )
       
     U(0, :) = exp(-200*(x-0.5)**2 ) 
     
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, Order = Order, & 
                                           Differential_operator =  Burgers_equation, & 
                                           Boundary_conditions   =  Burgers_BC,  Solution = U ) 
      
     call scrmod('revers')   
     call qplot(x, U(0,:), Nx+1)  
     call qplot(x, U(Nt,:), Nx+1)                                 


contains 
!----------------------------------------------------------
real function Burgers_equation( x, t, u, ux, uxx) 
        real, intent(in) ::  x, t, u, ux, uxx
            
            real :: nu = 0.0
            
          !  Burgers_equation = - u * ux + nu * uxx
            Burgers_equation = -  ux + nu * uxx
                 
end function 
!-------------------------------------------------------
real function Burgers_BC(x, t, u, ux)
  real, intent(in) :: x, t, u, ux 
  
        if (x==x0) then
                            Burgers_BC = u !- 1

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
subroutine Waves_equation_2D_old

       integer, parameter :: Nx = 20, Ny = 20, Nt = 1000, Nv = 2 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv)  
       
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
       real ::  t0 = 0, tf = 4.5 
       integer :: i, j, Order = 6 
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
     call Grid_Initialization( "nonuniform", "x", Order, x )
      call Grid_Initialization( "nonuniform", "y", Order, y )

  
     U(0, :, :, 1) = Tensor_product( exp(-15*x**2), exp(-15*y**2) ) 
     U(0, :, :, 2) = 0
  
  
     
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, N_variables = Nv, & 
                                           Differential_operator =  Wave_equation_2D, & 
                                           Boundary_conditions   =  Wave_BC2D,  Solution = U ) 
     
   
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $  \frac{\partial{^2 u}}{\partial{t^2}}   = \frac{\partial{^2 u}}{\partial{x^2}} + \frac{\partial{^2 u}}{\partial{y^2}}  $ ", 2)
     CALL TITLIN (" $ u(-1) = 0, \quad u(1) = 0 $  ", 4)
     
     do i=0, Nt, 20 
      call qplclr( U(i, 0:Nx, 0:Ny, 1) , Nx+1, Ny+1)                                  
     end do 

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
    
    
!*************************************************************************************************
!*
!*************************************************************************************************
subroutine Test_IVBP2D_plate

       integer, parameter :: Nx = 15, Ny = 15, Nt = 40, Nv = 3
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv)  
       
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
       real ::  t0 = 0, tf = 0.1 
       integer :: i, j, Order = 4 
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
     call Grid_Initialization( "nonuniform", "x", Order, x )
      call Grid_Initialization( "nonuniform", "y", Order, y )
  
       
     U(0, :, :, 1) = 0
     U(0, :, :, 2) = Tensor_product( exp(-5*x**2), exp(-5*y**2) ) 
     U(0, :, :, 3) = 0
     
  
     !write(*,*) maxval( U(0,:,:,1) ) , minval( U(0,:,:,1) ) 
     ! do j=0, Ny 
     !     write(*,'(100f8.3)' ) U(0, :, j, 1) 
     ! enddo 
    !  call qplclr( U(i, 0:Nx, 0:Ny, 1) , Nx+1, Ny+1)     
     
     call Initial_Value_Boundary_ProblemS( Time_Domain = Time, x_nodes = x, y_nodes = y, Order = Order, N_variables = Nv, & 
                                           Differential_operator =  Wave_equation_2D, & 
                                           Boundary_conditions   =  Wave_BC2D,  Solution = U ) 
     
     
      !write(*,*) maxval( U(0,:,:,1) ) , minval( U(0,:,:,1) ) 
      !do j=0, Ny 
      !    write(*,'(100f8.3)' ) U(0, :, j, 1) 
      !enddo 
   
     
     call scrmod('revers')                          
     call metafl('xwin')
     call disini() 
     CALL TEXMOD('ON')
     CALL TITLIN (" $  \frac{\partial{^2 u}}{\partial{t^2}}   = \frac{\partial{^2 u}}{\partial{x^2}} + \frac{\partial{^2 u}}{\partial{y^2}}  $ ", 2)
     CALL TITLIN (" $ u(-1) = 0, \quad u(1) = 0 $  ", 4)
     
     do i=0, Nt, 5 
      call qplclr( U(i, 0:Nx, 0:Ny, 1) , Nx+1, Ny+1)                                  
   !  call qplclr( U(Nt, :, :, 1),  Nx+1, Ny+1) 
     end do 

contains 
!----------------------------------------------------------
function Wave_equation_2D( x, y, t, u, ux, uy, uxx, uyy, uxy ) result(L)
        real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
        real :: L(size(u))
            
            real :: w, dwdt, d2wdx2
            
            
            w = u(1);   dwdt = u(2); d2wdx2 = u(3) 
            
           
        
        
            L(1)  = dwdt 
            L(2)  = - uxx(3)  - uyy(3)
            L(3)  =   uxx(2)  + uyy(2)
       
                 
end function 
!-------------------------------------------------------
function Wave_BC2D(x,y, t, u, ux, uy) result(BC) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BC( size(u) ) 
  
       real ::  w, dwdt, d2wdx2 
            
       w = u(1)
       dwdt = u(2) 
       d2wdx2 = u(3) 
       
  
        if (x==x0) then
                           BC(1) = w
                           BC(2) = dwdt
                           BC(3) = d2wdx2 

        else if (x==xf) then
                            BC(1) = w
                            BC(2) = dwdt
                            BC(3) = d2wdx2
        else if (y==y0) then
                            BC(1) = w
                            BC(2) = dwdt
                            BC(3) = d2wdx2

        else if (y==yf) then
                             BC(1) = w
                             BC(2) = dwdt
                             BC(3) = d2wdx2
        else
             write(*,*)  "Error in BC2D_waves"
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
     x(0) = x0;  x(Nx) = xf; y(0) = y0;  y(Ny) = yf;
     call Grid_Initialization( "nonuniform", "x", Order, x )
      call Grid_Initialization( "nonuniform", "y", Order, y )
       
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
     x(0) = x0;  x(Nx) = xf; y(0) = y0;  y(Ny) = yf;
     call Grid_Initialization( "nonuniform", "x", Order, x )
      call Grid_Initialization( "nonuniform", "y", Order, y )
     
           
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

    
    end module 