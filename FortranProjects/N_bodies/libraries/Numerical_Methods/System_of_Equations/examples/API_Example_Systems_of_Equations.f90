module API_Example_Systems_of_Equations

    use Linear_systems
    use Non_Linear_Systems 
    use Numerical_recipes

    implicit none
    
    contains    
 !line 9    
subroutine Systems_of_Equations_examples
    
    call LU_Solution
    call Newton_Solution
    call Test_Power_Method 
    call Vandermonde_condition_number
    
       
end subroutine   
  

!****************************************************************
! 
!****************************************************************
subroutine LU_Solution

    real :: A(4,4), b(4), x(4)
    
   
    A(1,:) = [ 4, 3, 6, 9]
    A(2,:) = [ 2, 5,  4, 2]
    A(3,:) = [ 1, 3, 2, 7]
    A(4,:) = [ 2, 4, 3, 8]
   
    b = [ 3, 1, 5, 2]
   
    call LU_factorization( A )    
    x = Solve_LU( A , b )
    
    write (*,*) 'The solution is = ', x
   
end subroutine 


!****************************************************************
! 
!****************************************************************
subroutine Newton_Solution

    real :: x0(3) = [1., 1., 1.  ]; 

    call Newton( F, x0 ) 
    
    write(*,*)  'Zeroes of F(x) are x = ', x0
   
end subroutine 


function F(xv)

    real, intent(in) :: xv(:)
    real:: F(size(xv))
 
    real :: x, y, z
 
    x = xv(1)
    y = xv(2)
    z = xv(3)
    
    F(1) = x**2 - y**3 - 2
    F(2) = 3 * x * y - z
    F(3) = z**2 - x
   
end function


!****************************************************************
! 
!****************************************************************
subroutine   test_power_method1 

 integer :: i, j, k  
 integer, parameter :: PI = 4 * atan(1d0) 
 integer, parameter :: N = 20
 real :: x(0:N), Vandermonde(0:N, 0:N), sigma   
 real :: a=-1, b=1 
 real V(0:N), V0(0:N) 
 
 
 x = [ ( a + (b-a)*i/N, i=0, N) ] 
 
 do i=0, N; do j=0, N; 
    Vandermonde(i,j) = x(i)**j
 end do; end do  
    
 V = 1 
 V0 = 0
 do while( abs(norm2(V)-norm2(V0)) > 1d-5 )  
     V0 = V
     V = matmul( Vandermonde, V ) / norm2(V) 
     write(*,*) maxval(V) 
 end do 
 sigma = dot_product( V, V ) 
 write(*,*) "sigma = ", sigma 

 
end subroutine 
!****************************************************************
! Milestone 2 
!****************************************************************
subroutine Test_Power_Method 
     
      integer :: i, j, k  
      integer, parameter :: N = 3 
      real :: A(N, N), B(N, N )   
      real :: lambda, U(N) 
      
      
      A(1,:) = [ 7, 4, 1 ] 
      A(2,:) = [ 4, 4, 4 ] 
      A(3,:) = [ 1, 4, 7 ] 
      
      U = 1
      
      call power_method(A, lambda, U) 
      
      A = A - lambda * Tensor_product( U, U ) 
      
      B = matmul( A, transpose(A) ) 
      
      B = matmul( transpose(A), A ) 
     
      U = [-1, 0,  1 ]
        
      call power_method(A, lambda, U) 
      
      A = A - lambda * Tensor_product( U, U ) 
      
      U = [1, 2,  1 ] 
       
      call power_method(A, lambda, U) 

end subroutine 

!****************************************************************
! Test eigenvalues 
!****************************************************************
subroutine   test_SVD2 

 integer :: i, j, k  
 integer, parameter :: N = 3 
 real :: A(N, N)
 real :: sigma(N), U(N, N), V(N,N)  
 
 
 A(1,:) = [ 7, 4, 1 ] 
 A(2,:) = [ 4, 4, 4 ] 
 A(3,:) = [ 1, 4, 7 ] 
 
 do i=1, N 
    write(*,'(A8, 3f8.3)') "A = ", A(i,:) 
 end do  
 
 call SVD(A, sigma, U, V) 

 do i=1, N 
    write(*,'(A8, f8.3, A15, 3f8.3, A15, 3f8.3)') "sigma = ", sigma(i), "U = ", U(i, :), " V = ",  V(i, :)
 end do   
  
end subroutine
!****************************************************************
! Test SVD
!****************************************************************
subroutine   test_Vandermonde

 integer :: i, j, k  
 integer, parameter :: N = 20
 real :: A(N, N), sigma(N), U(N, N), V(N,N)      
 
 
 
 do i=1, N; do j=1, N
    A(i,j) = (i/real(N))**j 
 end do; end do 

 call SVD(A, sigma, U, V) 
  
 do i=1, N 
    write(*,'(A8, e15.7)') "sigma = ", sigma(i)
 end do 
 
 
end subroutine

!****************************************************************
! Vandermonde_condition_number
!****************************************************************
subroutine Vandermonde_condition_number
   
    integer :: i, j, k  
    integer, parameter :: N = 6
    real :: A(0:N, 0:N), sigma(0:N), U(0:N, 0:N), V(0:N, 0:N)
    real :: A_SVD(0:N, 0:N), D(0:N, 0:N) 
    real :: kappa 
    
    

    do i=0, N; do j=0, N     
        A(i,j) = (i/real(N))**j 
    end do; end do 
    
    
    call SVD(A, sigma, U, V) 
    D = 0 
    do i=0, N
        D(i,i) = sigma(i) 
    end do 
    
    
    A_SVD = matmul(U , matmul( D ,  transpose(V) ) )
   
    call print_matrix("A", A) 
    call print_matrix("A_SVD", A_SVD) 
    call print_matrix("D", D) 
    
    kappa = Condition_number(A) 
    
    write(*,*) " Condition_ number =", sigma(0)/sigma(N) 
    write(*,*) " Condition_ number =", kappa
    
end subroutine

!*********************************************************************************
!*
!*********************************************************************************
subroutine print_matrix( legend, A) 
character(len=*), intent(in) :: legend 
real, intent(in) :: A(:,:) 
  
integer :: N, i 

N = size( A, dim=1) 

write(*,*) 
do i=1, N 
     write(*,'(A8, 100f8.3)') trim(legend)//" =    ", A(i,:) 
end do



end subroutine 

!****************************************************************
! Test SVD
!****************************************************************
subroutine   test_SVD 

 integer :: i, j, k  
 integer, parameter :: N = 3
 real :: A(N, N)
 real :: sigma(N), U(N, N), V(N,N), A_SVD(N,N), D(N,N), C(N,N), B(N,N)       
 
 
 A(1,:) = [ 7, 4, 1 ] 
 A(2,:) = [ 4, 4, 5 ] 
 A(3,:) = [ 1, 4, 7 ]
  
 
 call SVD(A, sigma, U, V) 
 
 D = 0 
 do i=1, N 
     D(i,i) = sigma(i) 
 end do 
 
 do i=1, N 
    write(*,'(A8, 100f8.3)') "A = ", A(i,:) 
 end do  
 
 do i=1, N 
    write(*,'(A8, f8.3)') "sigma = ", sigma(i)
 end do 
 
 !write(*,*) 
 ! do i=1, N 
 !    do j=1, N 
 !        write(*,'(A5, 2i3, 100f8.3)') "V dot V  = ", i, j, dot_product(V(:,i), V(:,j) )  
 !    end do 
 !end do
 
 C = matmul( V , matmul( matmul(D,D) , transpose(V) ) )
 do i=1, N 
    write(*,'(A20, 100f8.3)') "transpose(A) * A  = ", C(i,:) 
 end do 
 
 write(*,*) 
 B = matmul( transpose(A), A )
  do i=1, N 
    write(*,'(A20, 100f8.3)') "transpose(A) * A  = ", B(i,:) 
  end do 
 
 !write(*,*) 
 ! do i=1, N 
 !    do j=1, N 
 !        write(*,'(A5, 2i3, 100f8.3)') "U dot U  = ", i, j, dot_product(U(:,i), U(:,j) )  
 !    end do 
 ! end do
  
  !write(*,*) " B ="
  !B = matmul( transpose(U), U ) 
  !do i=1, N 
  !       write(*,'(A15, 100f8.3)') "transpose(U)  U  = ",  B(i,:)  
  !end do
  !
  !write(*,*) " B ="
  !B = matmul( U, transpose(U) ) 
  !do i=1, N 
  !       write(*,'(A15, 100f8.3)') " U transpose(U)   = ",  B(i,:)  
  !end do
  
  
 write(*,*) 
 B = matmul( A, transpose(A) )
  do i=1, N 
    write(*,'(A20, 100f8.3)') "A * transpose(A)  = ", B(i,:) 
  end do 
  
 write(*,*)  
 C = matmul( U , matmul( matmul(D,D) , transpose(U) ) )
 do i=1, N 
    write(*,'(A20, 100f8.3)') "A * transpose(A)  = ", C(i,:) 
 end do 
 
 
 
 !A_SVD = matmul( U , matmul( D ,  transpose(V) ) )
 
 A_SVD = matmul( transpose(U) , matmul( A ,  V ) )

 write(*,*) 
 do i=1, N 
    write(*,'(A8, 100f10.5)') "Sigma = ", A_SVD(i, :)
 end do 
 

 A_SVD = matmul(U , matmul( D ,  transpose(V) ) )

 write(*,*) 
 do i=1, N 
    write(*,'(A8, 100f10.5)') "A_SVD = ", A_SVD(i, :)
 end do 
 
 
 
end subroutine

!****************************************************************
! Test SVD
!****************************************************************
subroutine   test_SVD3 

 integer :: i, j, k  
 integer, parameter :: N = 2
 real :: A(N, N)
 real :: sigma(N), U(N, N), V(N,N), A_SVD(N,N), D(N,N), ATA(N,N), AAT(N,N)      
 
 
 A(1,:) = [ 0,  1 ] 
 A(2,:) = [ -1, 0 ] 
 
 
 do i=1, N 
    write(*,'(A8, 2f8.3)') "A = ", A(i,:) 
 end do  
 
 
 call SVD(A, sigma, U, V) 
 
 do i=1, N 
    write(*,'(A8, f8.3, A15, 2f8.3, A15, 2f8.3, A5, 2f8.3 )') "sigma = ", sigma(i), "U = ", U(:,i), " V = ",  V(:,i), " norm =", norm2(U(:,i)), norm2(V(:,i)) 
 end do   
 
 D = 0 
 do i=1, N 
     D(i,i) = sigma(i)
 end do 
 
  
 A_SVD = matmul( U , matmul( D , transpose(V) ) )
 ATA = matmul( V , matmul( D , transpose(V) ) )
 AAT = matmul( U , matmul( D , transpose(U) ) )
 
 do i=1, N 
    write(*,'(A5, 100f8.3)') "U = ", U(i,:) 
 end do
 do i=1, N 
    write(*,'(A5, 100f8.3)') "V = ", V(i,:) 
 end do
 
 
 do i=1, N 
    write(*,'(A20, 100f10.5)') "transpose(A) * A  = ", ATA(i,:) 
 end do 
 
 
  
 do i=1, N 
    write(*,'(A20, 100f10.5)') "A * transpose(A)   = ", AAT(i,:) 
 end do 
 
 do i=1, N 
    write(*,'(A8, 100f10.5)') "A_SVD = ", A_SVD(i, :)
 end do 
 
 
end subroutine

!****************************************************************
! Test eigenvalues 
!****************************************************************
subroutine   test_eigenvalues 

 integer :: i, j, k  
 integer, parameter :: N = 3 
 real :: A(N, N)
 real :: lambda(N), U(N, N) 
 
 
 A(1,:) = [ 7, 4, 1 ] 
 A(2,:) = [ 4, 4, 4 ] 
 A(3,:) = [ 1, 4, 7 ] 
 
 do i=1, N 
    write(*,'(A8, 3f8.3)') "A = ", A(i,:) 
 end do  
 
 
 call Eigenvalues(A, lambda, U) 
 
 do i=1, N 
    write(*,'(A8, f8.3, A15, 3f8.3, A10, f8.3)') "lambda = ", lambda(i), "eigenvector = ", U(:,i), " norm =", norm2(U(:,i)) 
 end do   
  
 
 
end subroutine 








real function Sphere(x,y,z)

    real, intent(in):: x, y, z
    
    Sphere = x**2 + y**2 + z**2 -1 
     
end function


subroutine Restricted_sphere
   
    real:: xi(1)
   
    xi = 0.5
       
    call Newton (f = g, x0 = xi)   
    
    write (*,*) "Zeroes of f(x) are x=", xi
    
 contains
 
   function g(x)
   real, intent(in):: x(:)
   real:: g(size(x))
   
   g = Sphere(x = x(1),y = 0. ,z = 0.)
   end function
    
    
      
 end subroutine
 
 
 
 
 
 
 
end module 
