module Graphs

!USE IFPORT  
use dislin 
implicit none 

type contour 
     real, allocatable :: x(:)
     real, allocatable :: y(:) 
     real, allocatable :: z(:,:) 
     real, allocatable :: levels(:) 
     character(len=50) :: label = ' z label ' 
end type 

type parametric 
     real, allocatable :: x(:)
     real, allocatable :: y(:) 
     character(len=50) :: label = ' parametric k ' 
     logical :: marker =.false. 
     character(len=10) :: color='white'
end type 

type graph  ! polymorphic graph  ( parametric   contour ) 

   character(len=50) :: cabecera(2) =[  '  ', '  ' ]
   character(len=50) :: label_x = ' leyenda x ' 
   character(len=50) :: label_y = ' leyenda y ' 
   logical :: log = .false.  
   character(len=10) :: device= 'xwin' 
   logical :: XWIN_active = .true.
   logical :: PDF_active = .true. 
   logical :: PNG_active = .false.
   character(len=100) :: file_name = ' '
   character(len=10) :: max_min = 'auto'  
   real :: x_max
   real :: x_min
   real :: y_max
   real :: y_min 
   
   
   type (parametric), allocatable :: p(:) 
   type (contour) :: c 
   
end type 

 
character (len=10) :: colores(9) = [ character(len=10) :: 'RED', 'GREEN', 'BLUE', 'CYAN', 'YELLOW', 'ORANGE', 'MAGENTA', 'WHITE', 'BLACK' ] 
  

contains 

 
!******************************************************************************************************
subroutine Pinta( id_graf,  g ) 
!******************************************************************************************************
integer, intent(in) :: id_graf
type (graph) :: g 

   write(*,*) " BEGIN: Pinta " 
   
   if (g%XWIN_active) then  
      
      write(*,*) " Plot xwin " 
      g%device= 'xwin' 
      call Pinta_graficas( id_graf, g) 
      
   endif 
    
   if (g%PDF_active) then  
   
      write(*,*) " Plot pdf " 
      g%device= 'pdf'
      call Pinta_graficas( id_graf, g) 
      
   endif    

   if (g%PNG_active) then  
   
      write(*,*) " Plot png " 
      g%device= 'png' 
      call Pinta_graficas( id_graf, g) 
      
   endif

      
end subroutine 

!******************************************************************************************************
subroutine Pinta_graficas( id_graf,  g ) 
!******************************************************************************************************
integer, intent(in) :: id_graf
type (graph), intent(in) :: g 

   character(len=2000) :: car 
   integer :: i, j 
   integer :: Np  ! numero de parametrics 
   integer :: Ncx, Ncy 
   character(len=20) :: graficas 
!   integer, parameter :: nmax=1000
!   integer :: NTRI, I1(nmax), I2(nmax), I3(nmax)  
   
   real :: zmax, zmin 
   integer :: n_levels 
   integer :: nx, ny, nxl, nyl  
   real :: xymax 
   
   
   
   real :: x0, xf, y0, yf 
   real :: x1, x2, dx, y1, y2, dy, xor, yor 

   if (allocated(g%p)) then 
                                   graficas = 'parametrics'
                                   Np = size ( g%p ) 
   elseif (allocated(g%c%x)) then 
                                   graficas = 'contour'
                                   Ncx = size( g%c%x )
                                   Ncy = size( g%c%y ) 
                                    
   else
         write(*,*) ' Error: No se ha cargado la grafica'
         stop  
   endif      
   write(*,*) ' Tipo de graficas =',   graficas     
                       
    
  if ( g % max_min == 'auto' ) then 
    if (graficas=='parametrics') then 
 
      x0 = 1d50; y0 = 1d50; xf = -1d50; yf = -1d50        
      do j=1, Np 
        do i=1, size( g%p(j)%x ) 
          if ( x0 > g%p(j)%x(i) ) then; x0 = g%p(j)%x(i) ;  endif  
          if ( xf < g%p(j)%x(i) ) then; xf = g%p(j)%x(i) ;  endif
          if ( y0 > g%p(j)%y(i) ) then; y0 = g%p(j)%y(i) ;  endif  
          if ( yf < g%p(j)%y(i) ) then; yf = g%p(j)%y(i) ;  endif 
        enddo 
      enddo
      call Nice_value(x0, xf) ; call Nice_value(y0, yf)
    
      
    elseif (graficas=='contour') then 
    
       
       x0 = g % c % x(1);  xf = g % c % x(Ncx) 
       y0 = g % c % y(1);  yf = g % c % y(Ncy) 
    
    endif  
 else if (g % max_min == 'dados' ) then 
 
      x0 = g % x_min 
      xf = g % x_max
      y0 = g % y_min 
      yf = g % y_max
     
 else 
      write(*,*) ' Error en pinta_graficas: max_min debe ser auto o dados ' 
      stop 
 endif              
      
   if (id_graf /= 0 ) then   
    call setxid(id_graf, 'widget')
    call WINSIZ ( 885 , 470 )
    
   end if  
      
     
 !   call csruni('plot') 
    nx = 4000 !2969 
    ny = 2099 ! 2099
    call page(nx, ny)
    
   
    
    call scrmod('reverse')
  
    
    if (g%cabecera(1) == ' ' .AND. g%file_name == ' ')  then 
          call Setfil ('.\doc\'//'graf.'//trim(g%device))
    elseif ( g%file_name == ' ') then        
          call Setfil ('.\doc\'//trim(g%cabecera(1))//'.'//trim(g%device))
    else      
          call Setfil ('.\doc\'//trim(g%file_name)//'.'//trim(g%device))
    endif 
  
                                      
    call METAFL(trim(g%device))
    CALL FILMOD ('delete')
     
       
    call disini 
    call erase
    call height(50)
    call hname(50) 
    call hwfont  
    call chacod('ISO1')
    call texmod('ON')  
    call errmod('all', 'off')
   
    
 !  call paghdr('pagina', 'p2', 2, 0) 
      
   ! call axspos( 650, 2099-200)
   
    if ( graficas=='parametrics') then 
    
       nxl = 3*nx/4 - 300 !2*nx/4 - 300
       nyl = 2*ny/3 
       call axslen(nxl, nyl)   
       call axspos( 400, ny-300)   !(1000, ny-300) 
       
       
        call titlin(trim(adjustl(g%Cabecera(1))), 2)
        call titlin(trim(adjustl(g%Cabecera(2))), 4)
       
    elseif( graficas=='contour') then 
     
       xymax = max( xf-x0, yf-y0) 
       nxl = int ( 3*ny/4.*(xf-x0)/xymax ) 
       nyl = int ( 3*ny/4.*(yf-y0)/xymax ) 
       call axslen(nxl, nyl)   
       call axspos(1000, ny-250) 
       
       call titlin(trim(adjustl(g%Cabecera(1))), 4)
    
    
    endif 
    
    call setscl( [x0, xf], 2, 'x' ) 
    call setscl( [y0, yf], 2, 'y' )
  
  
    

    if (g%log) CALL axsscl('LOG','xy')
      
    call name(trim(g%label_x),'X')
    call name(trim(g%label_y),'Y') 
     
    call graf( x1, x2, xor, dx,   y1, y2, yor, dy )
    call frame(10)
    call box2d 
    call grid(1,1 ) 

    
   ! call chncrv('color') 
    call thkcrv( 5 ) 
      
      
    if  (graficas=='parametrics') then    
      do j=1, Np 
        if (g%p(j)%marker) then 
                      call incmrk(1)
                      call marker(21)
        else  
                      call incmrk(0)
                      call marker(0)   
        endif               
                 
        call color(g % p(j) % color )                        
        call curve( g%p(j)%x(:), g%p(j)%y(:),  size( g%p(j)%x ) )
      enddo   

     call color('white') 
     call legini(car, Np, 20)
 

     if (Np > 1) then 
     
      do j=1, Np 
         call leglin(car, trim(g%p(j)%label), j)
      enddo  
      
      call frame(3)
      call legtit(' ')
      call legopt(2.0, 0.5, 1.0) ! tamanio de las legends de las parametrics 
      call legpos(nx-800, 400) !(nx-1000, 400)   ! posicion de las legends 
      call legend(car, 3)
     
    endif  
   
  elseif (graficas=='contour') then   
    
         !   call zscale(-4.0, 4.0 ) 
            n_levels = size(g%c%levels)
            zmax = maxval(g%c%levels) 
            zmin = minval(g%c%levels)  
            !CALL LABDIS (500, 'xyz') ! everything is shifted (nmumbers and labels) 
            CALL NAMDIS (100, 'z')
            CALL ZAXIS (zmin, zmax, zmin , (zmax-zmin)/10., nyl, g%c%label, 0, 0, nx-1000, ny-250)
            CALL shdmod('cell','CONTUR')
            write(*,*) Ncx, Ncy 
           ! read(*,*)
            CALL conshd(g%c%x, Ncx, g%c%y, Ncy, g%c%z, g%c%levels,  n_levels )
            
  endif           
    
  call title
  call disfin

      
end subroutine 

!*****************************************************************************
!*
!*****************************************************************************
subroutine Nice_value(x0, xf)  
  real, intent(inout) :: x0, xf 
  
   integer :: f 
   real :: dx 
   
   f = 1 
   dx = xf - x0 
   if (dx==0) then 
      dx = 0.1 
   endif 
   
   do 
   !   write(*,*) 'dx=', dx, f 
   !   read(*,*) 
      if (int(f*dx) > 0 ) then 
           
           call busca_proximos(f, 'x0', x0)
           call busca_proximos(f, 'xf', xf)
           exit 
      else 
           if (f>10) then 
              call busca_proximos(f, 'x0', x0)
              call busca_proximos(f, 'xf', xf)
              exit
           endif 
           f = f * 10 
      endif 
   enddo               
   
  
  
  
end subroutine    

!*****************************************************************************
!*
!*****************************************************************************
subroutine busca_proximos(f, c, x)  
  integer, intent(in) :: f 
  real, intent(inout) :: x
  character(len=*), intent(in) :: c
  
 !  integer :: signo  
 !  real :: distancia
 !  integer :: m 
   
   
  ! write(*,*) ' x  la entrada=', x, f 
    
  ! if (x>=0) then 
  !               signo = 1 
  !else 
  !              signo = -1 
  !endif    
   
    if (c=='x0') then  
                           
                           if (x>=0) then 
                                x = int(f*x)/real(f)  
                           else 
                                x = (int(f*x)-1)/real(f) 
                           endif        
    elseif (c=='xf') then 
                           if (x>0) then 
                                x = (int(f*x)+1)/real(f)  
                           else 
                                x = int(f*x) / real(f) 
                           endif       
    endif  
                        
   
 !  distancia = -1d30 
 !  m = 0 
  ! do  
   !   if (distancia > f*abs(x)- 5**m  ) then 
      
    !       if (c=='x0') then 
     !          x = signo * 5**(m-1) / f 
      !     elseif (c=='xf') then 
       !        x = signo * 5**m / f 
       !    endif     
       !    exit 
    !  else 
     !      distancia = f * abs(x) - 5**m  
     ! endif 
    !  m = m + 1 
  ! enddo               
   
 !  write(*,*) ' x  la salida=', x 
  
  
end subroutine 





end module 
