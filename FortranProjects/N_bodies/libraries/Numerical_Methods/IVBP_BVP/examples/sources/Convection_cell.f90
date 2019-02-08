module convection_cell

    use Boundary_value_problems
    use dislin
    use Finite_differences
	use Cauchy_problem

    implicit none

contains

subroutine test_convection_cell(Nx,Ny,Nt,TT)

       use dislin

       integer, intent(in) :: Nx, Ny, Nt
       real, intent(out) :: TT(0:,0:)


       integer :: M
       real :: x_nodes(0:Nx),y_nodes(0:Ny), a=-1, b=1, c=-1, d=1
       real, parameter :: Re=1, Pr=1, Ra=1E5
       integer :: i, Order = 10
       real :: t0 = 0., tf = 1.
       !integer, parameter :: Nt = 1000

       real,allocatable, save :: WW(:,:)!, TT(0:Nx,0:Ny) !vorticidad y temperatura
       real :: WWx(0:Nx,0:Ny), WWy(0:Nx,0:Ny), WWxx(0:Nx,0:Ny), WWyy(0:Nx,0:Ny)
       real :: FF(0:Nx,0:Ny), FFxx(0:Nx,0:Ny), FFyy(0:Nx,0:Ny)  !Función de corriente

       real :: U(0:Nx,0:Ny)  !velocidad X
       real :: V(0:Nx,0:Ny)  !velocidad Y

        real :: E(0:Nt)
        real :: time_1(0:Nt)
        integer:: ii


         !!   OPEN(UNIT=60, FILE="Results.txt", ACTION="write")


        M=(Nx+1)*(Ny+1)

        print*, Nx,Ny,Nt

        ii=0

        allocate(WW(0:Nx,0:Ny))

       call Grid_Initialization( "nonuniform", "x", Order, x_nodes )
       call Grid_Initialization( "nonuniform", "y", Order, y_nodes )

 !!!      write (60,*) "x="
 !!!             write (60,*) x_nodes
 !!!
 !!!      write (60,*) "y="
 !!!             write (60,*) y_nodes




       print*, "working..."

       call Cauchy_Problem_Solution( Domain = T_Domain, Initial_C = vorticity_IC, &
                                    System = System_1, Scheme= Runge_Kutta4, Outputs = vorticity_graphs)

                                 !call  Derivative( "y", 1, TT, Ty )
                                 !call  Derivative( "x", 1, TT, Tx )

                                 call High_order_Grid_deallocate
                                 deallocate(WW)


       !! close (60)
contains

!------------------------------------------------------------------------------
subroutine T_Domain( Time )
 real, pointer :: Time(:)

    integer :: i
    allocate ( Time(0:Nt) )

     Time = [ (t0 + (tf-t0)*i/(Nt*1d0), i=0, Nt ) ]

end subroutine
!------------------------------------------------------------------------
subroutine vorticity_IC( t0, X)
   real, intent(in) :: t0
   real, pointer :: X(:)

   integer :: i, j

     allocate(  X(2*M) )

  WW=0.
    do i=0, Nx
        TT(i,:)=(1+x_nodes(i)**3)/(b-a)
    enddo

      X(1:M)=reshape(WW, [ M ]);

      X(M+1:2*M)=reshape(TT, [ M ]);

end subroutine


!*******************************************************************
!*
!*******************************************************************
subroutine  System_1(t, X, F)
              real, intent(in)  :: t
              real, intent(inout)  :: X(:)
              real, intent(out) :: F(:)

           integer :: k,l

    real :: Wx(0:Nx,0:Ny), Wxx(0:Nx,0:Ny), Wxy(0:Nx,0:Ny), Wy(0:Nx,0:Ny), Wyy(0:Nx,0:Ny)
    real :: Tx(0:Nx,0:Ny), Ty(0:Nx,0:Ny), Txx(0:Nx,0:Ny), Tyy(0:Nx,0:Ny)
    real :: Fxy(0:Nx,0:Ny), Txy(0:Nx,0:Ny)

     WW = reshape( X(1:M), [Nx+1, Ny+1] )

     TT= reshape ( X(M+1:2*M), [Nx+1, Ny+1] )

      call Boundary_Value_Problem( x_nodes = x_nodes, y_nodes = y_nodes, Order = Order, &
                          Differential_operator = Equation,  &
                          Boundary_conditions1 = BCs1, &
                          Boundary_conditions2 = BCs2, Solution = FF )


      call  Derivative( "y", 1, FF, U )
      call  Derivative( "x", 1, FF, V )
          V=-V

      call  Derivative( "y", 2, FF, FFyy )
      call  Derivative( "x", 2, FF, FFxx )

    !boundary_conditions
      WW(0,:)= -(FFxx(0,:)+FFyy(0,:))
      WW(Nx,:)= -(FFxx(Nx,:)+FFyy(Nx,:))
      WW(:,0)= -(FFxx(:,0)+FFyy(:,0))
      WW(:,Ny)= -(FFxx(:,Ny)+FFyy(:,Ny))

      call  Derivative( "y", 1, WW, WWy )
      call  Derivative( "x", 1, WW, WWx )

      call  Derivative( "y", 2, WW, WWyy )
      call  Derivative( "x", 2, WW, WWxx )

      Fxy(0,:)= 0.
      Fxy(Nx,:)= 0.
      Fxy(:,0)= 0.
      Fxy(:,Ny)= 0.




      call Newmann( 'y', 0, TT, 0. )
      call Newmann( 'y', Ny, TT, 0. )

          TT(0,:)= 0.
      TT(Nx,:)= 1.

      call  Derivative( "y", 1, TT, Ty )
      call  Derivative( "x", 1, TT, Tx )
      call  Derivative( "y", 2, TT, Tyy )
      call  Derivative( "x", 2, TT, Txx )

!  *** inner grid points
        do k=1, Nx-1
             do l=1, Ny-1
                 Fxy(k,l) = -U(k,l)*WWx(k,l)-V(k,l)*WWy(k,l)+(1/Re)*(WWxx(k,l)+WWyy(k,l))-(Ra/(Re*Re*Pr))*Tx(k,l)
             enddo
         enddo

      Txy(0,:)= 0.
      Txy(Nx,:)= 0.
      Txy(:,0)= 0.
      Txy(:,Ny)= 0.

         do k=1, Nx-1
             do l=1, Ny-1
                 Txy(k,l) = -U(k,l)*Tx(k,l)-V(k,l)*Ty(k,l)+(1/(Re*Pr))*(Txx(k,l)+Tyy(k,l))
             enddo
         enddo

         F(1:M)=reshape(Fxy, [ M ] );
         F(1+M:2*M)=reshape(Txy, [ M ] );

end subroutine

!-------------------------------------------------------------------------------
subroutine vorticity_graphs(t, X, S )
          real, intent(in) :: t
          real, intent(inout) :: X(:)
          procedure (ODES_Outputs), optional :: S
          integer :: i, j
          real :: Ty(0:Nx,0:Ny)

          integer, save:: cont



  !
  if  (t==t0) then

  cont=0

  WW=reshape(X(1:M),[Nx+1,Ny+1])
        TT=reshape(X(M+1:2*M),[Nx+1,Ny+1])
            call Newmann( 'y', 0, TT, 0. )
           call Newmann( 'y', Ny, TT, 0. )
              call  Derivative( "y", 1, TT, Ty )


    call scrmod('revers')
    call qplcon( FF, Nx+1, Ny+1, 20);
    call qplcon( WW, Nx+1, Ny+1, 20);
    call qplcon( TT, Nx+1, Ny+1, 20);


 !!  write (60,*) "time=", t
 !!   write (60,*) "T="
 !!   do j=0,Ny
 !!      write (60,*) TT(:, j)
 !!   end do
 !!
 !!   write (60,*) "F="
 !!   do j=0,Ny
 !!      write (60,*) FF(:, j)
 !!   end do



  else

  cont=cont+1

     E(ii)=0

          do i=1, Nx-1
        do j=1, Ny-1
            E(ii)=E(ii)+(U(i,j)*U(i,j))+(V(i,j)*V(i,j))
        end do
        end do
        time_1(ii)= t


  end if



  if ((cont==100)) then


        cont=0
        WW=reshape(X(1:M),[Nx+1,Ny+1])
        TT=reshape(X(M+1:2*M),[Nx+1,Ny+1])
            call Newmann( 'y', 0, TT, 0. )
           call Newmann( 'y', Ny, TT, 0. )
              call  Derivative( "y", 1, TT, Ty )


    call scrmod('revers')
    call qplcon( FF, Nx+1, Ny+1, 20);
    call qplcon( WW, Nx+1, Ny+1, 20);
    call qplcon( TT, Nx+1, Ny+1, 20);


!!   write (60,*) "time=", t
!!    write (60,*) "T="
!!    do j=0,Ny
!!       write (60,*) TT(:, j)
!!    end do
!!
!!    write (60,*) "F="
!!    do j=0,Ny
!!       write (60,*) FF(:, j)
!!    end do



            call qplot(TT(Nx/2,:),y_nodes, Ny+1)

            call qplot(x_nodes, TT(:,Ny/2), Nx+1)

            call qplot(time_1, E, Nt)

    ! call qplcon( U, Nx+1, Ny+1, 20);
    ! call qplcon( V, Nx+1, Ny+1, 20);

  endif
ii=ii+1



end subroutine


!-----------------------------------------------------------------
!
!-----------------------------------------------------------------
real function Equation(x, y, U, Ux, Uy, Uxx, Uyy, Uxy)
           real, intent(in) :: x, y, u, ux, uy, uxx, uyy, uxy
           integer :: ix, iy

           call search_index(x, y, ix, iy)

            Equation = Uxx+Uyy+WW(ix,iy)

end function

real function BCs1(x, y, u, ux, uy)
           real, intent(in) :: x, y, u, ux, uy


   if (x==x_nodes(0)) then
                      BCs1 = uy
   elseif (y==y_nodes(0)) then
                      BCs1 = ux
   elseif (x==x_nodes(Nx)) then
                      BCs1 = uy
   elseif (y==y_nodes(Ny)) then
                      BCs1 = ux
   else
        write(*,*) " Error BCs x=", x
        write(*,*) " a, b=", a, b
        read(*,*)
   endif

end function

real function BCs2(x, y, u, ux, uy)
           real, intent(in) :: x, y, u, ux, uy


   if (x==x_nodes(0)) then
                      BCs2 = ux
   elseif (y==y_nodes(0)) then
                      BCs2 = uy
   elseif (x==x_nodes(Nx)) then
                      BCs2 = ux
   elseif (y==y_nodes(Ny)) then
                      BCs2 = uy
   else
        write(*,*) " Error BCs x=", x
        write(*,*) " a, b=", a, b
        read(*,*)
   endif

end function

subroutine search_index(x, y, ix, iy)

    real, intent(in):: x, y
    integer, intent(out):: ix, iy

    ix=0
    iy=0

    do while (abs(x-x_nodes(ix))>1E-5)
        ix=ix+1
    end do

    do while (abs(y-y_nodes(iy))>1E-5)
        iy=iy+1
    end do


end subroutine

end subroutine


end module


program main

    use convection_cell
    use dislin

    implicit none

    integer:: Nx, Ny,Nt
    integer:: i, j, k
    integer:: kk=4 !RungeKutta
    real, allocatable:: TT(:,:),T(:,:,:),R(:,:)



        Nx=10
        Ny=10

    allocate(T(0:Nx,0:Ny,2), TT(0:Nx,0:Ny), R(0:Nx, 0:Ny))

    do i=1,1

        Nt=10000*i

        print*, Nt

    call Test_convection_cell(Nx,Ny,Nt,TT)

    do j=0,Nx
    do k=0,Ny
    T(j,k,i)=TT(j,k)
    end do
    end do

    end do


    do i=0,Nx
    do j=0,Ny
    R(i,j)=(T(i,j,2)-T(i,j,1))/(2**kk-1)
    end do
    end do

     print*, "R="
     do j=Ny,0,-1
        print*, R(:, j)
     end do


 call scrmod('revers')
          call qplcon( R, Nx+1, Ny+1, 40);


end program main
