program TestSystemsofEquations

   
   use Linear_systems
   use Jacobian_module
   use Non_Linear_Systems 
   
   
   implicit none 

   !  linear system of equations 
      call Test_LU_factorization 
      call Test_Gauss  
      call Test_Gauss_indio
      call Test_Gauss_Sanz_Serna
     
    ! Jacobin of a vectorial function 
      call Test_Jacobian 
      
    ! Non linear system of equations 
      call Test_Newton


      Write(*,*) " Press any key to continue " 
      read(*,*) 


contains


!***********************************************************************
!* An example of a vector function to test the Jacobian
!***********************************************************************
function F( xv ) 
  real, intent(in) :: xv(:)
  real :: F(size(xv))  

  real :: x, y

  x = xv(1) 
  y = xv(2) 
  

  F(1) = 2 * x - 3 * y 
  F(2) = x + 2 * y 
 



end function 


!******************************************************************************
!* Test to validate the Jacobian 
!******************************************************************************
subroutine Test_Jacobian
  
integer, parameter :: N = 2 
real :: x(N) = [ 0, 1 ], J(N,N) 
integer :: i 


J = Jacobian( F, x) 

write(*,*)
write(*,*) ' ********************** TEST JACOBIAN ****************************** ' 
do i=1, N 

 write(*,'(a20, 3f15.7)')  ' Jacobian =', J(i,:) 

enddo 
write(*,*)


end subroutine 




!********************************************************************************   
! Dada la matriz A calculo su factorizacion LU 
!  A(1,:) = [ 8, 2, 9 ]
!  A(2,:) = [ 4, 9, 4 ]
!  A(3,:) = [ 6, 7, 9 ]

!****** L 
! 1 0 0
! 1/2 1 0
! 3/4 11/16 1

! ***** U 
! 8 2 9
! 0 8 −1/2
! 0 0 83/32
!*******************************************************************************
subroutine Test_LU_factorization  

 
  real :: A(3,3), b(3), x(3), Ac(3,3)     
  integer :: i 


  A(1,:) = [ 8, 2, 9 ]
  A(2,:) = [ 4, 9, 4 ]
  A(3,:) = [ 6, 7, 9 ]

  Ac = A 
  b = [1, 2, 3 ] 


 call LU_Factorization(A)
 
write(*,*) 
write(*,*) ' ********************** TEST LU factorization *********************** ' 
  
  write(*,*) ' Factorizacion L U ' 
  do i=1, 3 
     write(*,'(3f8.3)')  A(i, :) 
  end do 

  x = Solve_LU( A, b ) 
  write(*,'(a20,3f13.4)') 'The solution is: x = ',  x
  write(*,*) 
  
  write(*,'(a40,3f8.4)') 'Check : A x - b = ',  matmul(Ac,x) - b  




end subroutine 

!********************************************************************************   
! Dada la matriz A calculo su factorizacion LU 
!  A(1,:) = [ 8, 2, 9 ]
!  A(2,:) = [ 4, 9, 4 ]
!  A(3,:) = [ 6, 7, 9 ]

!****** L 
! 1 0 0
! 1/2 1 0
! 3/4 11/16 1

! ***** U 
! 8 2 9
! 0 8 −1/2
! 0 0 83/32
!*******************************************************************************
subroutine Test_Gauss  

 
  real :: A(3,3), b(3), x(3), Ac(3,3), bc(3)      
 


  A(1,:) = [ 8, 2, 9 ]
  A(2,:) = [ 4, 9, 4 ]
  A(3,:) = [ 6, 7, 9 ]

  Ac = A   
  b = [1, 2, 3 ] 
  bc = b 


  x = Gauss( A, b ) 

  write(*,*)
  write(*,*) ' ********************** GAUSS TEST  ****************************** ' 
  write(*,*) ' Matrix ' 
  write(*,*)  A(1,:) 
  write(*,*)  A(2,:) 
  write(*,*)  A(3,:)
  write(*,*) 
  write(*,*) ' The independent term  ' 
  write(*,*)  b(:)  
  write(*,*) 
  

  write(*,'(a20,3f13.4)') ' The solution is : x = ',  x
  write(*,*) 
  
  write(*,'(a40,3f8.4)') 'Check: A x - b = ',  matmul(Ac,x) - bc  
  write(*,*)

end subroutine 



!*******************************************************************************
subroutine Test_Gauss_indio  

 
  real :: A(3,3), b(3), x(3), Ac(3,3), bc(3)      
 

  A(1,:) = [ 20, 15, 10]
  A(2,:) = [ 0., 0.001, 8.5 ]
  A(3,:) = [0., 0. , 23375.]

  Ac = A   
  b = [45., 8.501, 23374.] 
  bc = b 


  x = Gauss( A, b ) 
  
  write(*,*)
  write(*,*) ' ********************** TEST GAUSS ****************************** ' 
  write(*,*) ' Matrix  ' 
  write(*,*)  A(1,:) 
  write(*,*)  A(2,:) 
  write(*,*)  A(3,:) 

  write(*,'(a20,3f13.4)') ' The solution is : x = ',  x
  write(*,*) 
  
  write(*,'(a40,3f8.4)') 'Check: A x - b = ',  matmul(Ac,x) - bc  
  write(*,*)



end subroutine 



!*******************************************************************************
subroutine Test_Gauss_Sanz_Serna

 
  real :: A(2,2), b(2), x(2), Ac(2,2), bc(2)      
 


  A(1,:) = [ 0.0001, 1.]
  A(2,:) = [ 1, 1 ]
  

  Ac = A   
  b = [1, 2] 
  bc = b 


  x = Gauss( A, b ) 

  write(*,*)
  write(*,*) ' ********************** GAUSS TEST  ****************************** ' 

  write(*,*) ' Matrix ' 
  write(*,*)  A(1,:) 
  write(*,*)  A(2,:) 
 
  write(*,'(a20,2f13.4)') ' The solution is : x = ',  x
  write(*,*) 
  
  write(*,'(a40,3f8.4)') 'Check: A x - b = ',  matmul(Ac,x) - bc  
  write(*,*)



end subroutine 



!*****************************************************************************
!* Test 
!*****************************************************************************
subroutine Test_Newton

  real :: x0(3) = [1., 1., 1.  ]; 

  write(*,*)
  write(*,*)  '********************** Test_Newton ************************************'
  write(*,*)  '             G(1) =  x**2 - y**2 + 1 '
  write(*,*)  '             G(2) = 2 * x *y           '
  write(*,*)  '             G(3) = z**2  - 2           '
  
  !write(*,'(a25, 3f13.6)')  'Un cero  de G(x) es x='  , Newton( G, x0 ) ! Error compiler 
  call Newton( G, x0 ) 
  write(*,'(a25, 3f13.6)')  'Zeroes of G(x)  x='  , x0
   
  write(*,*)  '********************** Press Enter ************************************'
  write(*,*)
 


end subroutine 
!***************************************************************************
!* Vector function for the test 
!***************************************************************************
function G(xv)  
 real, intent(in) :: xv(:) 
 real :: G(size(xv)) 
 
   real :: x, y, z 
   
   x =  xv(1) 
   y =  xv(2) 
   z =  xv(3) 
 
  G(1) = x**2 - y**2 + 1 
  G(2) = 2 * x *y 
  G(3) = z**2  - 2

end function


end program 

