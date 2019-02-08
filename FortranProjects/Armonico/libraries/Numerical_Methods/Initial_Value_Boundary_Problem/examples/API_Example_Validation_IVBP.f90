module API_Example_Validation_IVBP

    use Linear_systems
    use Initial_Value_Boundary_Problem
    use IVBP_Validation
    
implicit none

    contains
    
subroutine IVBP_validation_examples

  !  call optimization_state


   call Heat_equation_1D
   
   call Heat_equation_2D

   !call Waves_equation_2D
   

end subroutine

subroutine Heat_equation_1D

       integer, parameter :: Nx=2
       real ::  x(0:Nx)
       real ::  time(0:Nx)
       
       real, allocatable :: error(:), N_vec(:)
       real ::  x0 = -1, xf = 1
       real :: t0 = 0, tf = 4  
       integer :: i, j, k, Order
      
      ! write(*,*)  'Heat_equation_1D'
      ! write(*,*)
      ! write(*,*) 'Do you want to save times? [Yes=1/No=0]'
      ! read(*,*) k
      ! write(*,*)
        
     x(0)    = x0; x(Nx)    = xf
     time(0) = t0; time(Nx) = tf

     do order=2, 8, 2
         
     !    if (k==1) call Open_TimeSheet_1D_IVBP(order)
         
         
         call IBVP_1D_Validation(                                          &
                           Differential_operator =  Heat_equation,           &
                           Boundary_conditions   =  Heat_BC,                 &
                           Initial_condition     =  Heat_IC,                 &
                           Order = Order, Spatial_Domain = x,                &
                           Time_Domain= time, error = error, N_vec = N_vec) 
                       
          call Save_validation_IVBP_1D(N_vec, error, size(N_vec), order)
     
      end do   
contains 
!line 783----------------------------------------------------------
real function Heat_equation( x, t, u, ux, uxx) 
        real, intent(in) ::  x, t, u, ux, uxx
            
            real :: nu = 0.02
            
            Heat_equation =  nu * uxx
                 
end function 
!-------------------------------------------------------
real function Heat_BC(x, t, u, ux)
  real, intent(in) :: x, t, u, ux 
  
        if (x==x0) then
                            Heat_BC = u 

        else if (x==xf) then
                            Heat_BC = u

        else
             write(*,*)  "Error in Heat_BC"
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif

end function
!-------------------------------------------------------
function Heat_IC(x)
    real, intent(in) :: x(:)
    real :: Heat_IC(size(x))
    
    integer:: k
    
    do k=1, size(x)
    Heat_IC(k)  =  exp(-25*x(k)**2 )
    end do
        
    
end function



end subroutine 
    
subroutine Heat_equation_2D

      integer, parameter :: Ni = 2
      real :: x(0:Ni), y(0:Ni), t(0:Ni)
       
       real, allocatable :: error(:), N_vec(:)
       real :: x0 = -1, xf = 1, y0 = -1, yf = 1  
       real :: t0 = 0, tf = 4
       integer :: i, j, k, Order       
      
     t(0) = t0;  t(Ni) = tf
     x(0) = x0;  x(Ni) = xf; y(0) = y0;  y(Ni) = yf;
     
     !write(*,*)  'Heat_equation_2D'
     !write(*,*)
     !write(*,*) 'Do you want to save times? [Yes=1/No=0]'
     !read(*,*) k
     !write(*,*)
     
    do order=2, 6, 2
     
 !    if (k==1) call Open_TimeSheet_2D_IVBP(order)
     
     call IBVP_2D_Validation(                                                 &
                 Differential_operator =  Heat_equation,                      & 
                 Boundary_conditions   =  Heat_BC,                            &
                 Initial_condition     =  Heat_IC,                            &
                 Order = Order, Spatial_Domain_x = x, Spatial_Domain_y = y,   &
                 Time_Domain = t, error = error, N_vec = N_vec )
                 
      call Save_validation_IVBP_2D(N_vec, error, size(N_vec), order)
    end do
     
contains



!-line 995-------------------------------------------------------
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
!-------------------------------------------------------
function Heat_IC(x,y)
    real, intent(in) :: x(:), y(:)
    real :: Heat_IC(size(x), size(y))
    
    integer:: i,j
    
    Heat_IC  =  Tensor_product( exp(-25*x**2),  exp(-25*y**2) )
    
end function

end subroutine 
    end module
    