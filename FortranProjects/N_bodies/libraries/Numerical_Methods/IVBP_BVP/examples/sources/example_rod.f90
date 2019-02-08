module example_rod

use Boundary_value_problems

implicit none

contains

subroutine test_rod

   use dislin


    implicit none

    integer, parameter :: N=10
    real :: x(0:N), U(0:N)

    real :: x0=0d0 , xf=1d0
    integer :: i

        x = [ (x0 + (xf-x0)*i/N, i=0, N) ]


   call Linear_Boundary_Value_Problem( x_nodes=x, Order=4, Differential_operator= Equation,&
                                            Boundary_conditions= BCs, solution= U )

    write(*,*) " maxval, minval U =", maxval(U), minval(U)


       call scrmod('revers')
       call qplot(x, U, N+1)

    contains


            subroutine area(x, A, dA)

                real, intent(in) :: x
               real, intent(out) :: A, dA

              real, parameter :: A0=1., L=1.0

                A=A0*exp(-x/(1d0*L))

               dA= -A0/(1d0*L)*exp(-x/(1d0*L))

            end subroutine

    real function Equation(x, u, ux, uxx)

           real, intent(in) :: x, u, ux, uxx

           real :: A, dA

           call area(x, A, dA)

           Equation = uxx + ux*dA/(1d0*A)

    end function

    real function BCs(x, y, yx)

           real, intent(in) :: x, y, yx

           real, parameter :: E = 1. , P=1.
           real :: A, dA


         if (x==x0) then
                      BCs = y

         elseif (x==xf) then

                     call area(x, A, dA)
                      BCs = yx - P/(1d0*E*A)
         else
             write(*,*) " Error BCs x=", x
             write(*,*) " a, b=", x0, xf
             read(*,*)
         endif


    end function





end subroutine




end module example_rod
