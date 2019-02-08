module IVBP_Validation
    
    use Initial_Value_Boundary_Problem
    use Finite_differences
    
implicit none

abstract interface  

       real function DifferentialOperator1D(x, u, ux, uxx) 
                        real, intent(in) :: x, u, ux, uxx 
       end function  
       
       real function DifferentialOperator2D(x, y, u, ux, uy, uxx, uyy, uxy) 
                        real, intent(in) :: x, y, u, ux, uy, uxx, uyy, uxy
       end function  

        function DifferentialOperator2D_system(x, y, u, ux, uy, uxx, uyy, uxy) 
                        real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
                        real :: DifferentialOperator2D_system(size(u))
        end function
        
       real function      BC1D(x, u, ux) 
           real, intent(in) :: x, u, ux 
       end function
       
       real function      BC2D(x, y, u, ux, uy) 
           real, intent(in) :: x, y, u, ux, uy 
       end function  
       
        function BC2D_system(x, y, u, ux, uy) 
           real, intent(in) :: x, y, u(:), ux(:), uy(:) 
           real :: BC2D_system(size(u)) 
       end function  
       
        function Analytical_Solution1D(x)
            real, intent(in) :: x(:)
            real :: Analytical_Solution1D( size(x) )
        end function
        
        function Analytical_Solution2D(x, y)
            real, intent(in) :: x(:), y(:)
            real :: Analytical_Solution2D( size(x), size(y) )
        end function
        
        real function DifferentialOperator1DIVBP( x, t, u, ux, uxx) 
            real, intent(in) ::  x, t, u, ux, uxx
        end function
        
        real function DifferentialOperator2DIVBP(  x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy ) 
            real, intent(in) :: x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy 
        end function
        
        function IC1D(x)
            real, intent(in) :: x(:)
            real :: IC1D(size(x))
        end function
        
        function IC2D(x,y)
            real, intent(in) :: x(:), y(:)
            real ::  IC2D(size(x), size(y))
        end function
        
            
        
        real function BC1D_IVBP(x, t, u, ux)
            real, intent(in) :: x, t, u, ux
        end function
        real function BC2D_IVBP( x, y, t, U, Ux, Uy )
            real, intent(in) ::  x, y, t, U, Ux, Uy 
        end function
        
end interface

    contains

    
       
    
 



subroutine IBVP_1D_Validation( Differential_operator, Boundary_conditions, Initial_condition, Order, Spatial_Domain,Time_Domain,  error, N_vec)

     procedure (DifferentialOperator1DIVBP) :: Differential_operator
     procedure (BC1D_IVBP) ::  Boundary_conditions
     procedure (IC1D) ::  Initial_condition
     integer, intent(in) :: Order
     real, intent(in) :: Spatial_Domain(:), Time_Domain(:)
     real, allocatable, intent(out) :: error(:), N_vec(:)
     
     real :: x0 , xf, t0, tf
     real, allocatable :: x_nodes(:),  U(:,:), U_pre(:,:), U_1(:,:), U_2(:,:), error_i(:,:), Time(:)
     
     integer :: i, j, k, N, Nt, Nt0
    
    integer :: N0 = 20, ni=5
     
     Nt0 = N0*1000/20
     allocate( error(2:ni), N_vec(2:ni), error_i(0:Nt0,0:N0) )
     
     
     x0 = Spatial_Domain(1)
     xf = Spatial_Domain(size(Spatial_Domain))
     t0 = Time_Domain(1)
     tf = Time_Domain(size(Time_Domain) )

         

     do i = 1, ni
         
         N  = N0*i
         Nt = N*1000/20

        
        allocate ( x_nodes(0:N), U(0:Nt,0:N), U_1(0:Nt0,0:N) )
        
        x_nodes(0) = x0; x_nodes(N) = xf
    
         call Grid_Initialization( grid_spacing = "uniform", direction = "x", q = Order, nodes = x_nodes )
         Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]
         
         U(0, :)  =  Initial_condition(x_nodes)

         call Initial_Value_Boundary_ProblemS(                          & 
         
                        Time_Domain = Time, x_nodes = x_nodes,          &
                        Order = Order,                                  & 
                        Differential_operator =  Differential_operator, & 
                        Boundary_conditions   =  Boundary_conditions,   & 
                        Solution = U ) 
        

        if (i==1) then
            
            allocate( U_pre(0:Nt,0:N), U_2(0:Nt0,0:N) )
            U_pre = U
            
        else
            do k=0, Nt0
                U_1(k,:) = U(i*k, :)
                U_2(k,:) = U_pre((i-1)*k,:)
            end do    
            
            do j=0, N0
                error_i(:,j) =  U_1(:, i*j) - U_2(:, (i-1)*j)
            end do
            
            
            error(i) = abs( norm2(error_i) )
            N_vec(i) = real(N)
           
            deallocate(U_pre, U_2)
            allocate(U_pre(0:Nt, 0:N), U_2(0:Nt0, 0:N))
            U_pre=U
        end if
        
        
    write(*,*) 'N = ', N, 'error = ', abs( norm2(error_i) ) 
    deallocate( x_nodes, U, U_1)
    
    end do
   
end subroutine




subroutine IBVP_2D_Validation( Differential_operator, Boundary_conditions, Initial_condition, Order, Spatial_Domain_x, Spatial_Domain_y, Time_Domain,  error, N_vec)

     procedure (DifferentialOperator2DIVBP) :: Differential_operator
     procedure (BC2D_IVBP) ::  Boundary_conditions
     procedure (IC2D) ::  Initial_condition
     integer, intent(in) :: Order
     real, intent(in) :: Spatial_Domain_x(:), Spatial_Domain_y(:), Time_Domain(:)
     real, allocatable, intent(out) :: error(:), N_vec(:)
     
     real :: x0 , xf, y0, yf, t0, tf
     real, allocatable :: x_nodes(:), y_nodes(:), Time(:),  U(:,:,:), U_pre(:,:,:), U_1(:,:,:), U_2(:,:,:), error_i(:,:,:), U_1x(:,:,:), U_1xy(:,:,:) , U_2x(:,:,:), U_2xy(:,:,:)
     
     integer :: i, j, k, N, Nt, Nt0, N_pre , i2
    
    integer :: N0 = 8, ni=4
     
    
    Nt0 = N0 * 200/20
    
     allocate( error(2:ni), N_vec(2:ni), error_i(0:Nt0,0:N0,0:N0) )
     
      
     x0 = Spatial_Domain_x(1)
     xf = Spatial_Domain_x(size(Spatial_Domain_x))
     y0 = Spatial_Domain_y(1)
     yf = Spatial_Domain_y(size(Spatial_Domain_y))
     t0 = Time_Domain(1)
     tf = Time_Domain(size(Time_Domain) )

     
     do i = 1, ni
         
         N  = N0*i
         Nt = N * 200/20
         !Nt = N**2/2
         !Nt = Nt0*(i**2)
         
         write(*,*) 'Nx = ', N, 'Nt = ', Nt
        allocate ( x_nodes(0:N), y_nodes(0:N), Time(0:Nt), U(0:Nt,0:N,0:N), U_1(0:Nt0,0:N,0:N) )
        
        
        x_nodes(0) = x0; x_nodes(N) = xf
        y_nodes(0) = y0; y_nodes(N) = yf
    
         call Grid_Initialization( "uniform", "x", Order, x_nodes )
         call Grid_Initialization( "uniform", "y", Order, y_nodes )
         
        Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]
        

         U(0, :,:)  =  Initial_condition(x_nodes,y_nodes)

         call Initial_Value_Boundary_ProblemS(                                          & 
                     Time_Domain = Time, x_nodes = x_nodes, y_nodes = y_nodes, Order = Order, &
                     Differential_operator =  Differential_operator,                    & 
                     Boundary_conditions   =  Boundary_conditions,  Solution = U ) 
        
        write(*,*) 'obtained', i

        if (i==1) then
            
            allocate( U_pre(0:Nt,0:N,0:N), U_2(0:Nt,0:N,0:N))
            U_pre = U
            N_pre=N
            
        else
            allocate(U_1x(0:Nt0, 0:N0, 0:N), U_2x(0:Nt0, 0:N0, 0:N_pre), U_1xy(0:Nt0, 0:N0, 0:N0),  U_2xy(0:Nt0, 0:N0, 0:N0) ) 
            U_1x=0; U_1xy=0; U_2x=0; U_2xy=0
            
            do k=0, Nt0
                U_1(k,:,:) = U(i*k, :,:)
                U_2(k,:,:) = U_pre((i-1)*k, :,:)
            end do
            do k=0, N0
                U_1x(:,k,:) = U_1(:, i*k, :)
                U_2x(:,k,:) = U_2(:, (i-1)*k, :)
            end do
            do k=0, N0
                U_1xy(:,:,k) = U_1x(:, :, i*k)
                U_2xy(:,:,k) = U_2x(:, :, (i-1)*k)
            end do
             error_i = 0
             error_i = U_1xy - U_2xy
            
            error(i) = abs( norm2(error_i) )
            N_vec(i) = real(N)
            write(*,*) 'N = ', N, 'error = ', error(i)
           
            
            
            deallocate(U_pre, U_2)
            allocate(U_pre(0:Nt, 0:N, 0:N), U_2(0:Nt0, 0:N, 0:N))
            
            U_pre=U
            N_pre=N

            deallocate(U_1x, U_1xy, U_2x, U_2xy)
        end if
        
        

    deallocate( x_nodes, y_nodes, Time, U, U_1)
     end do
     
    deallocate( U_pre, error_i)
   
end subroutine








subroutine crea_name(name_fichero_, numero, extension, name_creado)

     character,intent(in)    :: name_fichero_*(*)
     character,intent(in)    :: extension*(*)
     integer, intent(in)     :: numero
     character, intent(out)  :: name_creado*(*)

     integer     :: numero_caracteres
     character   :: indice_texto*2

     if (numero == 0) then
         name_creado = trim(name_fichero_)//extension
     else
         call intcha(numero, numero_caracteres, indice_texto)
         indice_texto = adjustl(indice_texto)
         name_creado = trim(name_fichero_)//&
                         indice_texto(1:numero_caracteres)//extension
     endif


endsubroutine





subroutine Save_validation_IVBP_1D(x, y, N, q)
    real, intent(in) :: x(:), y(:)
    integer, intent(in) :: N, q

    integer :: i
    character(len=100) :: name, names
    character(len=4) :: extension = '.txt'
    
    name = 'IVBP_1D_Validation_'
    
    call crea_name(trim(name), q, extension, names)
    
    open(62, file = names)
    
        do i = 1, N
            write (62,*) x(i), y(i)
        end do
        
    close (62)
    
    write(*,*) 'Saved as: ', names

end subroutine

subroutine Save_validation_IVBP_2D(x, y, N, q)
    real, intent(in) :: x(:), y(:)
    integer, intent(in) :: N, q

    integer :: i
    character(len=100) ::  name, names
    character(len=4) :: extension = '.txt'
    
    name = 'IVBP_2D_Validation_'
    
    
    call crea_name(trim(name), q, extension, names)
    
    open(62, file = names)
    
        do i = 1, N
            write (62,*) x(i), y(i)
        end do
        
    close (62)
    
    write(*,*) 'Saved as: ', names

end subroutine


 end module
