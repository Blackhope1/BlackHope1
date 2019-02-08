
!***********************************************************************************************
! Example program to test finite differences of a given Order based on  a set of points : x, y
!***********************************************************************************************
program Example_FD

 
  use Finite_differences
  use Lagrange_interpolation
  use Non_uniform_grids
 
  
  implicit none 
  
  
   write(*,*) " Test Non Uniform Grids " 
   call Test_Non_uniform_grids
 
   
   write(*,*) " Test Runge Phenomenon " 
   call Test_Runge_phenomenon 
   
   
   write(*,*) " Test Derivatives 2D " 
   call  Test_Derivatives2D 
   
 
  
contains 

!******************************************************************************
!*
!******************************************************************************
subroutine Test_Non_uniform_grids  
  
    integer, parameter :: N = 30
    integer :: Max_Order = N 
    real :: x_nodes(0:N), Dx(1:N), index(1:N)
    integer :: i
    integer :: Order
   
    

    do Order = 4, Max_Order, 2 
    
      
       call Grid_Initialization( "nonuniform", "x", Order, x_nodes  ) 
       write(*,*) "Order = ", Order, " N = ", N  
       
       do i=1, N 
          index(i) = i 
          Dx(i) = x_nodes(i) - x_nodes(i-1)
       end do 
       
       call scrmod('revers') 
       call metafl('xwin')
       call disini() 
       CALL TITLIN ("Dx_nodes", 2)
       call qplot( index, Dx, N )
       write(*,*) x_nodes
       
   enddo
     
 
end subroutine


!******************************************************************************
!*
!******************************************************************************
subroutine Test_Runge_phenomenon

    integer Order1
    integer N
    real, allocatable :: x_grid(:) 
    integer ::  j 

    
    real, allocatable :: W(:), Wx(:), Wxx(:)
    real :: Pi = 4 * atan(1.)  
   
    write(*,*) " Number of grid points 0..N "
    N = 30 
    write(*,*) "  N =  ", N 
  !  read(*,*) N 
    
    write(*,*) " Order of approximation "
    Order1 = 25
    write(*,*) " Order =  ", Order1 
 !   read(*,*) Order1 
    
    allocate( x_grid(0:N), W(0:N), Wx(0:N), Wxx(0:N) )
   
  
    call Grid_Initialization( "uniform", "x", Order1, x_grid  ) 

    
    write(*,*) " Uniform grid points  " 
    W =  1./ (1 +25*x_grid**2); 
  
    call Derivative("x", 1, W, Wx) 
    call Derivative("x", 2, W, Wxx)
    
    call scrmod('revers') 
    call disini() 
    CALL TITLIN (" W = 1/(1+25*x**2) interpolated in an uniform grid ", 2)
    call qplot( x_grid,  W,   N+1 )
    call disini() 
    CALL TITLIN (" dW/dx interpolated in an uniform grid ", 2)
    call qplot( x_grid,  Wx,  N+1 ) 
    call disini() 
    CALL TITLIN (" d2W/dxx interpolated in an uniform grid ", 2)
    call qplot( x_grid,  Wxx, N+1 )

    call Grid_Initialization( "nonuniform", "x", Order1, x_grid  ) 
 
    write(*,*) " Nonuniform grid points  " 
    W =  1./ (1 +25*x_grid**2); 
   
 
    call Derivative("x", 1, W, Wx) 
    call Derivative("x", 2, W, Wxx)
    
    call scrmod('revers') 
    call disini() 
    CALL TITLIN (" W = 1/(1+25*x**2) interpolated in an nonuniform grid ", 2)
    call qplot( x_grid,  W,   N+1 )
    call disini() 
    CALL TITLIN (" dW/dx interpolated in an nonuniform grid ", 2)
    call qplot( x_grid,  Wx,  N+1 )
    call disini() 
    CALL TITLIN (" d2W/dxx interpolated in an nonuniform grid ", 2)
    call qplot( x_grid,  Wxx, N+1 )

    



end subroutine 



 

!******************************************************************************
!*
!******************************************************************************

subroutine Test_Derivatives2D 

  integer, parameter :: N= 50; 
  real ::  x(0:N), y(0:N)

  real :: U(0:N), Ux(0:N), Uxx(0:N), dop(0:N) 
  real :: W(0:N, 0:N), Wx(0:N, 0:N), Wxx(0:N, 0:N), Wyy(0:N, 0:N), doperator(0:N, 0:N) 

  integer :: Order = 20, i  
  real :: pi = 4*atan(1.) 

  
  ! *** Define grid nodes 
        call Grid_Initialization( "uniform", "x", Order, x )
        call Grid_Initialization( "nonuniform", "y", Order, y )
 
  ! *** 1D function to be interpolated 
        U =  cos(pi*x); 
  
  ! *** 1D derivatives of U and difference operator: dop 
        call Derivative( "x", 1, U, Ux) 
        call Derivative( "x", 2, U, Uxx) 
  
        dop = Ux * ( U + Uxx )  ! 1D spatial operator 

        call scrmod('revers') 
        call disini() 
        CALL TITLIN (" U = cos( pi x)  ", 2)
        call qplot(x, U, N+1)
        call disini() 
        CALL TITLIN (" W = Ux ( U + Uxx ) ", 2)
        call qplot(x, dop, N+1)
        
        
  ! *** 2D function to be interpolated  
        forall(i=0:N) W(i,:) = sin(pi*x(i)) * sin(pi*y)
  
  ! *** partial derivatives of W and difference operator: doperator 
        call Derivative( (/"x","y"/), 1, 1, W, Wx)
        call Derivative( (/"x","y"/), 1, 2, W, Wxx)
        call Derivative( (/"x","y"/), 2, 2, W, Wyy)
  
        doperator = W * Wxx + Wx * Wyy  ! 2D spatial operator 
        
        call disini() 
        CALL TITLIN (" W = sin(pi*x) * sin(pi*y) ", 2)
        call qplcon( W, N+1, N+1, 20);
        
        call disini() 
        CALL TITLIN (" V = W * Wxx + Wx * Wyy ", 2)
        call qplcon( doperator, N+1, N+1, 20);
 
 end subroutine 
 


 
 
 
 
end program 






