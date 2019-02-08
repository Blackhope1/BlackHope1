module API_Example_Finite_Differences

    use Finite_differences
    use Non_Linear_Systems
    
    implicit none
    
contains  
  
!line 9 
subroutine Finite_difference_examples

   
   call Derivative_function_x
   call Derivative_function_xy
   call Derivative_error   
   call BVP_FD
   

end subroutine

subroutine Derivative_function_x

    integer, parameter :: Nx = 20, Order = 4
    real :: x(0:Nx)
    real :: x0 = -1, xf = 1
    integer :: i
    real :: pi = 4 * atan(1.) 
    real :: u(0:Nx), ux(0:Nx), uxx(0:Nx), ErrorUx(0:Nx), ErrorUxx(0:Nx)
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]   
    
    call Grid_Initialization( "nonuniform", "x", Order, x )
         
    u = sin(pi * x) 
       
    call Derivative( 'x' , 1 , u , ux ) 
    call Derivative( 'x' , 2 , u , uxx )  


    ErrorUx =  ux - pi* cos(pi * x) 
    ErrorUxx = uxx + pi**2 * u
    
    call scrmod("reverse")
    call qplot(x, ErrorUx , Nx+1)
    call qplot(x, ErrorUxx, Nx+1)



end subroutine



!***********************************************************************
!
!***********************************************************************

subroutine Derivative_function_xy

    integer, parameter :: Nx = 20, Ny = 20, Order = 6
    real :: x(0:Nx), y(0:Ny)
    real :: x0 = -1, xf = 1, y0 = -1, yf = 1
    integer :: i, j
    real :: pi = 4 * atan(1.0) 
    real :: u(0:Nx,0:Ny), uxx(0:Nx,0:Ny), uy(0:Nx,0:Ny), uxy(0:Nx,0:Ny)
    real :: Erroruxx(0:Nx,0:Ny), Erroruxy(0:Nx,0:Ny)
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]    
    
    call Grid_Initialization( "nonuniform", "x", Order, x )
    call Grid_Initialization( "nonuniform", "y", Order, y )
    
    
    do i=0, Nx; do j=0, Ny     
        u(i,j) = sin(pi * x(i)) * sin(pi * y(j))
    end do; end do    
   
    call Derivative( ["x", "y"], 1, 2, u, uxx )  
    call Derivative( ["x", "y"], 2, 1, u, uy ) 
    call Derivative( ["x", "y"], 1, 1, uy, uxy ) 

    Erroruxx = uxx + pi**2 * u
    
    do i=0, Nx; do j=0, Ny 
       Erroruxy(i,j) = uxy(i,j) - pi**2 * cos(pi*x(i))  * cos(pi*y(j))
    end do; end do    
    
end subroutine


    !call scrmod("reverse")
    !call qplcon( Uxx, Nx+1, Ny+1, 20)
    !call qplcon( Erroruxx, Nx+1, Ny+1, 20)







!********************************************************************
!* 
!line 99*****************************************************************
subroutine Derivative_error
 
 integer :: q = 2                  ! interpolant order 
 integer :: N                      ! # of nodes (piecewise pol. interpolation) 
 integer :: k = 2                  ! derivative order 
 integer :: p = 0                  !  where error is evaluated p=0, 1,...N
 integer, parameter :: M = 100     ! number of grids ( N = 10, .... N = 10**4) 
 real :: log_Error(M),  log_dx(M)  ! Error versus Dx 
 real :: epsilon = 1d-12           ! order of the random perturbation 
 
 real :: PI = 4 * atan(1d0), logN  
 integer ::  j
 
 real, allocatable :: x(:), f(:), dfdx(:)  ! function to be interpolated 
 real, allocatable :: dIdx(:)              ! derivative of the interpolant 
 
 do j=1, M 
     
   logN = 1 + 3.*(j-1)/(M-1)   
   N = 2*int(10**logN)
   
   allocate( x(0:N), f(0:N),  dfdx(0:N), dIdx(0:N) ) 
   x(0) = -1; x(N) = 1 
 
   call Grid_Initialization( "uniform", "x", q, x ) 
   
   call random_number(f)
   f = cos ( PI * x ) + epsilon * f
   dfdx = - PI**2 * cos ( PI * x ) 
   
   call Derivative( "x", k, f, dIdx )
   
   log_dx(j) = log( x(1)-x(0) ) 
   log_Error(j) = log( abs(dIdx(p)-dfdx(p)) ) 
 
   deallocate( x, f, dIdx, dfdx ) 
 
 end do 
 
 call scrmod("reverse")
 call qplot( log_dx, log_Error, M )
 
end subroutine 





!***********************************************************************
!
!***********************************************************************

subroutine Laplace_Equation_FD

    integer, parameter :: Nx = 20, Ny = 20, Order = 6
    real :: x_nodes(0:Nx), y_nodes(0:Ny)
    real :: x0 = -1, xf = 1, y0 = -1, yf = 1
    integer :: i, j, M
    real :: pi = 4 * atan(1.) 
    real :: u(0:Nx,0:Ny) 
    real, allocatable :: U_1(:)
    
    M = (Nx+1)*(Ny+1)
    
    allocate ( U_1(M) )
    
    x_nodes = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y_nodes = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]    
    
    call Grid_Initialization( "nonuniform", "x", Order, x_nodes )
    call Grid_Initialization( "nonuniform", "y", Order, y_nodes ) 

  
    U_1 = reshape(u, [ M ]);
  !  call Newton(System_BVP, U_1)
    
    u = reshape( U_1, [Nx+1, Ny+1] )
    
    
    call qplcon( u(:,:), Nx+1, Ny+1, 20)
    
    deallocate( U_1 )
    
contains

!-----------------------------------------------------------------------
 function Laplace_operator( x, y, u, uxx, uyy ) 
                real, intent (in) :: x, y, u, uxx, uyy
                real :: Laplace_operator
                
    
               Laplace_operator = uxx + uyy - u**4

 end function 

!-----------------------------------------------------------------------
  function Boundary_conditions( x, y, u,  ux, uy )
           real, intent (in) :: x, y, u, ux, uy
           real :: Boundary_conditions
            

            if (x==x0) then
            Boundary_conditions = u - ( cos(pi*y) + 1 ) 
        elseif (x==xf) then
            Boundary_conditions = u - ( cos(pi*y**2) + 1 ) 
        elseif (y==y0) then
            Boundary_conditions = uy
        elseif (y==yf) then
            Boundary_conditions = uy
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, xf=", x0, xf
            write(*,*) " y0, yf=", y0, yf

            stop 
        endif

   end function


!-----------------------------------------------------------------------
 function System_BVP(U) result(F)
                real, intent (in) :: U(:)
                real :: F(size(U))
  
                real :: UU(0:Nx,0:Ny), FF(0:Nx,0:Ny)
                
                UU = reshape( U, [Nx+1, Ny+1] )

                call Difference_equation(x_nodes, y_nodes, UU, FF)

                F=reshape(FF, [ M ]);

 end function 

 

 
 
!-------------------------------------------------------------------
 subroutine Difference_equation(x, y, W, Fxy)
            real, intent(in) :: x(0:), y(0:), W(0:,0:)
            real, intent(out) :: Fxy(0:, 0:)
            
  integer ::  k, l
  
           real :: Wx(0:Nx,0:Nx), Wy(0:Nx,0:Ny), Wxx(0:Nx,0:Ny), Wyy(0:Nx,0:Ny)
    
          
            call Derivative( ["x","y"], 1, 1, W, Wx  )
            call Derivative( ["x","y"], 1, 2, W, Wxx )

            call Derivative( ["x","y"], 2, 1, W, Wy  )
            call Derivative( ["x","y"], 2, 2, W, Wyy )



    !  ***  boundary conditions
            do k=0, Nx
                Fxy(k,  0) = Boundary_conditions( x(k), y(0),  W(k,  0),  Wx(k,  0),  Wy(k,  0) )
                Fxy(k, Ny) = Boundary_conditions( x(k), y(Ny), W(k, Ny),  Wx(k, Ny),  Wy(k, Ny) )
            enddo

            do l=0, Ny
                Fxy(0,  l)  = Boundary_conditions( x(0),  y(l),  W(0,  l),  Wx(0,  l),  Wy(0,  l) )
                Fxy(Nx, l)  = Boundary_conditions( x(Nx), y(l),  W(Nx, l),  Wx(Nx, l),  Wy(Nx, l) )
            enddo


     !  *** inner grid points
            do k=1, Nx-1
                do l=1, Ny-1
                    Fxy(k,l) = Laplace_operator( x(k), y(l),  W(k, l), Wxx(k,l), Wyy(k,l) )
                enddo
            enddo
            
         
 end subroutine
 

    
end subroutine

















!***********************************************************************
!
!***********************************************************************

!line 309-------------------------------------------
subroutine BVP_FD

    integer, parameter :: Nx = 20,  Order = 6 
    real :: x_nodes(0:Nx)        ! Grid definition
    real :: x0 = -1, xf = 1      ! Spatial domain
    integer :: i                 
    real :: pi = 4 * atan(1.)    
    real :: u(0:Nx)              ! Solution u(x)
          
    x_nodes = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ] 

    call Grid_Initialization( "nonuniform", "x", Order, x_nodes )

    u = 1                        ! Initial guess
    
    call Newton(System_BVP, u)   ! Difference equation 
                                 ! resolution

    call qplot( x_nodes, u, Nx+1)! Solution plot
    
contains 
 
!line 332-------------------------------------------------
 function System_BVP(U) result(F)
    real, intent (in) :: U(:)
    real :: F(size(U)) 
    
    call Difference_equation(x_nodes, U, F)

 end function  

!line 341-------------------------------------------------
 subroutine Difference_equation ( x, W, F)
    real, intent(in) :: x(0:Nx), W(0:Nx)
    real, intent(out):: F(0:Nx)
    
    real :: Wxx(0:Nx)
    
    call Derivative( "x" , 2 , W , Wxx )   ! Derivative computation
    
    F = Wxx + W   ! D.O. on inner points
    F(0) = W(0) - 1   ! B.C. on x = -1
    F(Nx) = W(Nx)     ! B.C. on x =  1

 end subroutine
 
    
end subroutine

    



end module 
