module Sum_series 
    
    implicit none 
    
    
    contains 
    
!!*****************************************************************************
!* Summation of three numerical series
!*****************************************************************************
subroutine  Summation_examples 

   write(*,*) "Summation from n=1 to n= infinity "
   
   write(*,*) "an= 1/n**2,               Sn= ", Summation_n2()
   write(*,*) "an= (-1.)**(n+1) / 2**n,  Sn= ", Summation_2n()
   write(*,*) "an= 1/factorial(n),       Sn= ", Summation_factorialn()

    
end subroutine    

!*****************************************************************************
!* Summation of  Sn = sum from {n=1} to {infinity} 1/ n**2 
!*****************************************************************************
real function Summation_n2()

  integer :: n = 1    ! index of the sum
  real :: an = 1      ! general term 
  real :: Sn = 0      ! summation of first n terms 
  real :: eps = 1d-7  ! smallest term to sum 

  do while( abs(an) >    Sn * eps ) 
     
      an =  1. / n**2
      Sn = Sn + an 
      n = n + 1 
 
  end do 
 
  Summation_n2 = Sn 

end function  

!*****************************************************************************
!* Summation of  Sn = sum from {n=1} to {infinity} (-1.)**(n+1) / 2**n  
!*****************************************************************************
real function Summation_2n( )

  integer :: n = 1    ! index of the sum
  real :: an = 1      ! general term 
  real :: Sn = 0      ! summation of first n terms 
  real :: eps = 1d-7  ! smallest term to sum 

  do while( abs(an) >    Sn * eps ) 
     
      an =  (-1.)**(n+1) / 2**n 
      Sn = Sn + an 
      n = n + 1 
 
  end do 
 
  Summation_2n = Sn 

end function  

!*****************************************************************************
!* Summation of  Sn = sum from {n=1} to {infinity} 1/ factorial(n) 
!*****************************************************************************
real function Summation_factorialn() 

  integer :: n = 1    ! index of the sum 
  real :: an = 1      ! general term 
  real :: Sn = 0      ! summation of first n terms 
  real :: eps = 1d-7  ! smallest term to sum 

  do while( abs(an) >    Sn * eps ) 
     
      an =  an / n 
      Sn = Sn + an 
      n = n + 1 
 
  end do 
 
  Summation_factorialn = Sn 

end function  
  
end module   
