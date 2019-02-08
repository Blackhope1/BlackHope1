module Dynamic_allocation 
  
    implicit none 
  
    contains

!*****************************************************************************
!* Examples of dynamic allocation 
!*****************************************************************************   
subroutine Matrices_allocation   

  real, allocatable :: A(:, :), Ak(:, :), B(:, :), Identity(:,:)   ! matrices
  
  integer ::  M       !  matrices are of variable dimension M x M 
  real :: S           ! it stores the trace of some matrix 
  integer :: i        ! index of row or column 
  integer :: k        ! index of sum 
  
! sum from M=0 to M=10 trace( A_M )   
  S = 0 
  do M=0, 10 
    call initialialization(M, A) 
    S = S + trace(A) 
  end do
  write(*,*) "sum traces ( A_M ) = ", S 
 
! sum from M=0 to M=5 trace( A**2_M ) 
  S = 0 
  do M=0, 5
    call initialialization(M, A) 
    S = S + trace( matmul(A,A) ) 
  end do 
  write(*,*) "sum traces ( A_M **2 ) = ", S 
  
! trace ( sum from k=0 to =10 A_5 **k  ) 
  M = 5  
  allocate( B(0:M, 0:M), Ak(0:M, 0:M), Identity(0:M, 0:M) ) 
  call initialialization(M, A)
  Identity = 0;   do i=0, M;  Identity(i,i) = 1; end do 
  
  Ak = Identity 
  B = Identity 
  
  do k=1,  10
              Ak = matmul( Ak, A )  
             B  = B +  Ak  
  end do 
  
  write(*,*) "trace ( sum  A_5**k ) = ", trace(B) 
 
end subroutine 
!*****************************************************************************
! Initialization of Vandermonde matrix A of dimension MxM  
!***************************************************************************** 
subroutine initialialization(M, A) 
  integer, intent(in) :: M 
  real, allocatable, intent(out) :: A(:,:) 

  integer :: i, j 
 
  allocate ( A(0:M, 0:M) ) 
  
  do i=0, M; do j=0, M
      A(i,j) =  ( i / real(M) ) **j 
  end do; end do 
  
  

end subroutine 
!*****************************************************************************
! Trace of a matrix A of dimension NxN  
!*****************************************************************************   
real function trace( A ) 
     real, intent(in) :: A(:, :)
     
    integer :: i           ! row of matrix A 
    integer :: N           ! dimension of matrix A 
    real    :: S = 0       ! trace  
    
    N = size(A, dim=1) 
    
    do i=1, N 
        S  = S + A(i,i) 
    end do 
    
    trace = S 
    
end function     

end module   