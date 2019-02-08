program TestLagrangeInterpolation 

   use Lagrange_interpolation
   
   
   implicit none 




   !  Lagrange interpolation  
      call Test_Lagrange_interpolation 
      
      call Test_Chebyshev_nodes 
      
   
      Write(*,*) " Press any key to continue " 
      read(*,*) 
      
contains



!******************************************************************************
!*
!******************************************************************************
subroutine Test_Lagrange_interpolation 
  
integer, parameter :: N = 5
real :: x(0:N) = [ 0, 1, 2, 3, 4, 5 ]
real :: y(0:N) = [ 0, 1, 4, 9, 16, 25 ]
real :: xp = 2
integer :: j 
integer :: s
real :: suma=0, Int 

integer, parameter :: Nk = 4 ! maximum order of derivative and width of the stencil 
 
real :: Weights(-1:Nk, 0:Nk, 0:N) ! Lagrange coefficients and their derivatives at xp 

write(*,*) '*********Test_Lagrange_interpolation***********'  
do j=0, N
  if (mod(Nk,2)==0) then 
                          s = max( 0, min(j-Nk/2, N-Nk) )   
  else 
                          s = max( 0, min(j-(Nk-1)/2, N-Nk) )
  endif 
!  write(*,'(a6,i3, a15, i3)' ) '** j = ', j, ' Stencil =', s   
  xp = x(j) ! x(s+Nk/2)
  Weights(-1:Nk, 0:Nk, j)  = Lagrange_polynomials( x = x(s:s+Nk), xp = xp ) 

 
  write(*,*) "Lagrange coefficients and their derivatives at xp "
  write(*,*) "Stencil of length :", Nk+1 

  write(*,*) 
  write(*,'(a10, 5f10.2, a6, f10.2)') 'Nodes =', x(s:s+Nk), 'xp =', xp 
  write(*,'(a10,f10.2)') 
  write(*,'(a25, 10f10.1)')  ' (l_j) xp =',              Weights(0, 0:Nk, j)
  write(*,'(a25, 10f10.1)')  ' (dl_j/dx) xp  =',         Weights(1, 0:Nk, j)
  write(*,'(a25, 10f10.1)')  ' (d2l_j/dx2) xp  =',       Weights(2, 0:Nk, j) 
  write(*,'(a25, 10f10.1)')  ' integral(l_j) 0..xp =',  Weights(-1, 0:Nk, j) 
  write(*,*) 


  write(*,'(a25, 2f10.3)') 'Interpolate value at xp =', xp, sum ( Weights(0, 0:Nk, j) * y(s:s+Nk) ) 
  write(*,'(a25, 2f10.3)') 'First derivative at xp =',  xp, sum ( Weights(1, 0:Nk, j) * y(s:s+Nk) )
  write(*,'(a25, 2f10.3)') 'Second derivative at xp =',  xp, sum ( Weights(2, 0:Nk, j) * y(s:s+Nk) )
  Int = sum ( Weights(-1, 0:Nk, j) * y(s:s+Nk) )
  write(*,'(a25, 2f10.3)') 'Integral from x(p-1) to xp =',  xp, Int 
  suma   = suma  + Int 
  write(*,*) " Press enter "
  read(*,*) 

enddo 

write(*,'(a15, f8.3, a8, f8.3, a10, f8.3 )' ) 'Integral from  ',  x(0), 'to xn =',  x(N), 'Int  = ',  suma 


write(*,*) " ************* Press enter ****************** "
read(*,*) 
 

end subroutine 

!******************************************************************************************
!*
!******************************************************************************************
 subroutine Test_Chebyshev_nodes 


 integer, parameter :: N = 200, Ng = 4, Ne = 91   
 integer :: i, k, l  
 integer ::   Np = 10

 real, allocatable :: x_nodes(:)
 real :: Error(0:N,0:1), x(0:N), E(0:2), maxError(Ne) 
 character(len=40) ::  tittle(Ng) 
 real :: pi 
 real ::  N_error(Ne), logNe(Ne)  
 character(len=10) :: graphs(Ng) = [ 'first', 'next ', 'next ', 'last ' ]
 
 pi = 4 * atan(1d0)     

 N_error = [ (10. + (Ne-1) * i/(Ne-1.), i=0, Ne-1) ] 
 logNe= log10( N_error ) 

 x =  [ (-1.0 + 2.0 * i/N, i=0, N )  ] 


 do k=1, Ng 

 do l=1, Ne
  if (allocated(x_nodes)) deallocate( x_nodes ) 
  Np = N_error(l)
  allocate( x_nodes(0:Np) )
  if (k==4) then      
       tittle(k) = ' Equally spaced nodes  ' ! ( canuto pag 67)
      x_nodes = [ (-1.0 + 2.0 * i/Np, i=0, Np )  ] 

  else if (k==2)  then 
       tittle(k) = ' Chebyshev-Gauss-Lobato nodes '  
       x_nodes = [ (cos(pi*i/Np), i=Np, 0, -1 )  ] 

  else if (k==3) then 
       tittle(k) = ' Chebyshev-Gauss ' 
       x_nodes = [ (cos(pi*(0.5+i)/(Np+1)), i=Np, 0, -1 )  ] 

  else if (k==1) then  
       tittle(k) = '  Chebyshev-Gauss-Radua '  
       x_nodes = [ (cos(pi*i/(Np+0.5)), i=Np, 0, -1 )  ] 

  endif 

  !do i=0, Np 
  !   write(*,*) 'x_nodes =', x_nodes(i) 
  !enddo 
  !read(*,*) 

 do i=0, N 
   E(0:2) = Lagrange_error_polynomial( x = x_nodes, xp = x(i) )
   Error(i,0:1)  = abs(E(0:1)) 

 enddo  
 
 maxError(l) = log10 ( maxval ( Error(:,0) ) )
  
! call qplot( x, Error(0:N,1), N+1 )



enddo  

 
enddo 


end subroutine 


end program 
