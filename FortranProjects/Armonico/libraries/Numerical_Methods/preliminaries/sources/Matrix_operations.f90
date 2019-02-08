module Matrix_operations 
    
    implicit none 
    
    contains
    
!*****************************************************************************
!* Examples of vectorial and matrix operations
!*****************************************************************************   
subroutine Examples_matrix_operations 

   integer, parameter :: N = 20 
   real :: W(N), V(N), A(N,N)   
   integer :: i, j  
   logical :: diagonal(N,N) = .false. 

  V = [ ( 1./i**2, i=1, N ) ] 
  W = [ ( (-1)**(i+1)/(2*i+1.), i=1, N ) ] 

  do i=1, N; do j=1, N; 
    A(i,j) = (i / real(N) ) ** j; 
  end do; end do  

  do i=1, N  
      diagonal(i,i) = .true. 
  end do 
  

 
  write(*,*) "Sum of componentes of W >0 = ",  sum(W, W>0)
  write(*,*) "Dot product of V, W = ",  dot_product(V, W)
  write(*,*) "Dot product of A(1:N, 3), W = ",  dot_product( A(1:N, 3), W)  
  write(*,*) "sum ( matmul(A,W)  = ",  sum( matmul(A, W) ) 
  write(*,*) "trace of A   = ", sum( A, diagonal ) 
  

end subroutine    
!*****************************************************************************
! Dot product of two vectors of dimension N 
!    sum from i=1 to i=N ( u_i * v_i ) 
!*****************************************************************************   
real function my_dot_product( u, v ) 
     real, intent(in) :: u(:), v(:) 
     
     
    integer :: i       ! summation index 
    integer :: N       ! dimension of vectors u, v 
    real    :: S = 0   ! summation 
    
    N = size(u) 
    
    
    do i=1, N 
        S = S + u(i) * v(i)   
    end do 
    
    my_dot_product = S 
    
end function    

!*****************************************************************************
! Matrix multiplication  of two matrices A, B  of dimension NxM and MxL  
!    C_ij = A_ik * B_kj
!*****************************************************************************   
function my_matmult( A, B ) result(C) 
     real, intent(in) :: A(:,:), B(:,:) 
     real C( size(A, dim=1), size(B, dim=2) ) 
     
    integer :: i        ! row of matrix C 
    integer :: j        ! column of matrix C 
    integer :: k        ! summation index 
    integer :: N,M, L   ! dimension of A and B  
    real    :: S = 0    ! summation 
    
    N = size(A, dim=1) 
    M = size(A, dim=2) 
    L = size(B, dim=2) 
    
    do i = 1, N 
        do j = 1, L 
            
            S = 0 
            do k = 1, M 
                S = S + A(i,k) * B(k,j) 
            end do  
            C(i,j) = S 
            
        end do 
    end do 
    
end function  


end module   