module mass_spring_damper

    use Cauchy_Problem
    use dislin

    implicit none

       private

            real, save :: t0=0., tf=10.
            integer, save :: order= 20;

            integer,save :: N=2000;   !Time steps

            real, save, allocatable :: tv(:), V(:,:)  ! graphs
            integer, save :: it =0;

       public:: test_mass_spring_damper




    contains

        subroutine Time_Domain ( Time )

            real , pointer :: Time (:)
            integer :: i

                allocate ( Time (0:N) )

                Time = [ (t0 + (tf -t0 )*i/(1d0*N), i=0, N ) ]

                allocate( tv(0:N) )
                allocate( V(2, 0:N) );  !graphs

        end subroutine

       subroutine problem_IC (t0 ,U)

          real , intent (in) :: t0
          real , pointer :: U(:)

           allocate( U(2) )

              U(1)=1
              U(2)=0

       end subroutine

        subroutine system_mass_spring_damper (t,U,F)

            real , intent (in) :: t
            real , intent ( inout ) :: U (:)
            real , intent ( out) :: F (:)
            real , parameter :: M=1.0 , D=.0 , K =10.0
            real :: P
            real :: pi=4*atan(1.0);

                 P=0.
                !P= 1.*sin(pi*t)

                F (1)= U(2)
                F (2)=(P-D*U(2) -K*U (1))/ M

        end subroutine

        subroutine results(t,U,S)

             real, intent(in) :: t
             real, intent(inout) :: U(:)
             procedure (ODES_Outputs), optional :: S


                tv(it) = t
                V(:, it) = U(:)

                it=it+1

           if (t==tf) then

                 call scrmod('revers')
                 call qplot( tv(1:N), V(1,2:N), N );
                 call qplot( tv(1:N), V(2,2:N), N  );

           endif

        end subroutine

        subroutine test_mass_spring_damper

            call Cauchy_Problem_Solution( Domain = Time_Domain,  Initial_C = problem_IC, System = system_mass_spring_damper, &
                                              Scheme=Runge_Kutta4, Outputs=results)

        end subroutine


end module
