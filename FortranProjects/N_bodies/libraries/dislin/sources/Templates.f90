module Templates

!USE IFPORT  
use dislin 
use Fonts_and_formats
implicit none 


  type coordinates 
      integer :: x 
      integer :: y 
  end type 

  type template 
        
        character(len=100) :: template_name 
        character(len=100) :: file_name = ' '
        
        integer :: N 
        type (coordinates) :: upper_left(100) 
        type (coordinates) :: lower_right(100) 
        integer :: x(50)
        integer :: y(50) 
        character(len=256) :: description(100)
        character(len=20) :: color(100) = 'white' 
        real :: content(100) = 0 
        logical :: editable(100) = .true.  
        integer :: font_size = 25 
        logical :: fixed_space = .false. 
        integer :: width  = int(3.16*885)
        integer :: N_layers = int(3.16*940)
       
        character(len=10) :: device= 'xwin'
        
        logical :: XWIN_active = .true. 
        logical :: PDF_active = .true. 
        logical :: PNG_active = .false.
        
        
  end type 
  
    type (template) :: p 

    contains 
    
!*******************************************************************************************************
!*
!******************************************************************************************************
subroutine Load_template( file_name, ulx, uly, lrx, lry, message, content  )
         character(len=*), intent(in) :: file_name, message(:) 
         real, intent(in) ::  content(:) 
         integer, intent(in) ::  ulx(:), uly(:), lrx(:), lry(:)
            
        
   integer ::  N
   real :: f 

   N = size(content) 
  
   p % template_name = file_name
   p % file_name = file_name 


 write(*,*) " file name = ", p % file_name,  p % template_name
  !
 
  f = 1.0
  p % N =  N 
  p % width = int(885*f) 
 ! p % N_layers = 523*f 
   p % N_layers = int(533*f)  
  
  p % upper_left(1:N) % x  = ulx  
  p % upper_left(1:N) % y  = uly
  
  p % lower_right(1:N) % x = lrx
  p % lower_right(1:N) % y = lry
  
      
  p % description(1:N) = message 
  p % content(1:N) = content 
  
  p % font_size =  45 ! 100 !
  p % fixed_space = .false. 
  
        
end subroutine    
    
    

!******************************************************************************************************
subroutine Print_templates( id_template ) 
!******************************************************************************************************
  integer, intent(in) :: id_template
 
   ! write(*,*) " Print template device XWIN = ", p%XWIN_active
  
 !  if (p%XWIN_active) then  
   
      p%device= 'xwin' 
      call Plot_template( id_template, p )
      
 !  endif
  
 !  if (p%PDF_active) then  
   
      p%device= 'pdf' 
      call Plot_template( id_template, p )
      
 !  endif    

 !  if (p%PNG_active) then  
   
 !     p%device= 'png' 
 !     call Plot_template( id_template, p ) 
      
 !  endif

      
end subroutine 



!************************************************************************************************************
!*
!************************************************************************************************************ 
 subroutine  Plot_template( id, p ) 
  integer, intent(in) :: id 
  type (template) :: p 
  
     integer :: i 
     integer :: ix, iy 
     integer :: font_N_layers 
     
     real :: fx, fy 
     integer :: nx, ny
   !  real :: x1 = 0., x2 = 1. , dx = 0.1,   y1 = 0. , y2 = 1. , yor = 0. , dy = 0.1 ! xor = 0., 
   !  integer :: ox 
     real :: xver 
  !   integer :: n1, n2 
  !   integer :: nx0, ny0 
   !  real :: f 
     character(len=400) :: content
     integer :: length
  
   write(*,*) 'BEGIN: template =', trim(p%template_name), p%device 
  
   call metafl(trim(p%device))
   CALL FILMOD ('delete')
  
   if (id /= 0 ) then   
    call setxid(id, 'widget')
    
   end if  
           
    call scrmod ('revers')
    
    if ( p%file_name == ' ')  then 
          call Setfil ('.\doc\template.'//trim(p%device))
    else        
          call Setfil ('.\doc\template_'//trim(p%file_name)//'.'//trim(p%device)) 
    endif 
    
    call IMGFMT ('RGB')
  
   
    !write(*,*) 'page = ', p%width, p%N_layers
    nx = 4000 
    ny = 4000*533/885
   ! nx = 885 
  !  ny = 533 
    call page(nx, ny)
 !   write(*,*) "Templates set page nx, ny =", nx, ny 
    
    
    
    
    
    
    fx =  nx / p%width 
    fy =  ny / p%N_layers
    
    
   ! call page( p%width, p%N_layers )
  !  CALL HWPAGE ( p%width, 2*p%N_layers)
    
    call disini 
     
      
       call DUPLX
       
       select case (trim(p%device))
           case ('xwin')
             call winfnt('Arial Bold')
           case ('pdf') 
             call psfont('Helvetica-Bold')
           !  call winfnt('Arial Bold')
            case ('png')
             call winfnt('Arial Bold')    
           
           case default
           
           end select 
   !   call winfnt('Arial Bold') 
    call getver(xver); 
 !   write(*,*) ' Dislin version = ', xver 
    
            
       call errmod('all', 'off')
       call erase;  !  call frame(1);  call filbox(0, 0, 1, 1 )
 !      write(*,*) 'template =', p%template_name 
       call frame(0)
     !  nx = p % width
      ! ny = p % N_layers 
       call filbox(0,0, 4000*nx, 2000*ny)
       call incfil( '.\templates\'//trim(p%template_name)//'.png' ) 

       
          
       ! call N_layers(p%font_size)  
        if (p%fixed_space) then 
                                call fixspc(1.0) 
        endif                         
       call texmod('ON')  
    
       call getpag(nx, ny) 
       write(*,*) "Templates get page nx, ny =", nx, ny 
       
    ! call messag("HELLO", nx/2, ny/2) 
   !  call messag("HELLO", 1, 1)  
       
     !  call axslen(nx, ny)
   !    call axspos(1, ny) 
       
       
    !   CALL SETGRF('none','none','none', 'none')
  !     call graf( x1, x2, xor, dx,   y1, y2, yor, dy )
       
     !  call getpos(nx0, ny0)
       
    !   write(*,*) nx0, ny0 
       
   !   call rlmess("Hello", 0.5, 0.5) 
       
       
      ! call frmess(1)
       do i=1, p%N 
            call color( p%color(i) ) 
            call gethgt(font_N_layers) 
            ix =  int ( ( 1.5* p%upper_left(i)%x +  0.5 * p%lower_right(i)%x )/2 ) 
            iy =  int ( p%upper_left(i)%y  ) ! +  p%lower_right(i)%y )/2 ! - font_N_layers / 2 
         
         !   write(*,*) p % description(i)
            
        !   if ( p%description(i) /=" " ) then 
               content = p % description(i)
               call messag(trim(content), nx*ix/885, ny*iy/533) 
                
         !   endif 
             length = len(trim(content)) * 50
            
           if ( p%content(i) /=-123456789 ) then 
       !       call rlnumb(  p%content(i), 3,  fx,  fy)  
               !write(*,*) p % content(i) 
               !read(*,*) 
               write(content, '(f7.3)') p%content(i)
               call messag(trim(content), nx*ix/885 + length, ny*iy/533) 
                
            endif 
            
          
            
            ! write(*,*) " ix, iy =",  ix, iy 
          !   fx = real(ix) / p%width   
           !  fy = 1 - real(iy) / p%N_layers
          !  call mensaje( p%description(i), p%content(i),  ix,  iy) 
             
        !    call rlmess(  p%description(i),  fx,  fy)
       !     ox = nlmess( p%description(i) ) 
            
            !write(*,*) 
            
        !    fx = fx  + real(ox) /  nx ! p%N_layers
         !   if ( p%content(i) /=-123456789 ) then 
       !       call rlnumb(  p%content(i), 3,  fx,  fy)  
                
       !     endif 
           
       enddo   
    
               
   call disfin
   write(*,*) ' END:  template'
 
 
 end subroutine 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

!*******************************************************************************************************
!*
!******************************************************************************************************
subroutine Edita_toda_la_template( id_g, p )
integer, intent(in) :: id_g  
type (template) :: p
 


! integer :: ix, iy 
 character(len=10) ::  car 
 logical :: encontrado
 integer :: i  

!10  call disini 

   encontrado = .false. 
   do i=1, p%N 
   
       if (  p % editable(i) ) then 
                  p%color(i) = 'red' 
                  call Plot_template( id_g, p )  
                  call dwgtxt( 'salida=s', car) 
                   p%color(i) = 'white' 
                  if (car=='s') then 
                         exit 
                  endif
                  
                  p % content(i) = character_to_real( car ) 
                  encontrado  = .true. 
                 
                
       endif
         
    enddo   
    call Plot_template( id_g, p )  
   
  !  call disfin 
   
   
 end subroutine 

!*******************************************************************************************************
!*
!******************************************************************************************************
subroutine Edita_template( id_g, p )
integer, intent(in) :: id_g  
type (template) :: p
 


 integer :: ix, iy 
! character(len=10) ::  car 
 logical :: encontrado 

write(*,*) 'entra en edita template' 
    call disini 
   
    do 
        call csrpt1(ix, iy) 
        call widget_txt_en_template(ix, iy, id_g, p, encontrado ) 
        if (encontrado) then 
          call Plot_template( id_g, p ) 
          exit 
        endif           
    enddo 
   
    call disfin 
    
  write(*,*) 'sale de edita_template'  
    
   ! goto 10 
    
   
 end subroutine 
       
!*******************************************************************************************************
!*
!******************************************************************************************************
subroutine widget_txt_en_template( ix, iy, id_g, p, encontrado ) 
integer :: ix, iy, id_g  
type (template) :: p
logical :: encontrado 

character(len=10) ::  car 
!real :: d 
integer :: i 
logical :: dentro

   encontrado = .false. 
   do i=1, p%N 
   
       if (  p % editable(i) ) then 
              dentro = inside ([ix, iy], [p%upper_left(i)%x, p%upper_left(i)%y], [p%lower_right(i)%x, p%lower_right(i)%y] ) 
              !distancia( [ix, iy], [p % x(i), p % y(i) ] )
              
            !  if (d > 0 .and. d < 100 ) then 
              if  (dentro) then 
              
                   p%color(i) = 'red' 
                   call Plot_template( id_g, p )  
               !   call Plot_template( id_g, p )  
                  call dwgtxt( 'salida=s', car) 
                  p%color(i) = 'white' 
                  if (car=='s') then 
                         encontrado  = .true. 
                         exit 
                  endif    
        
                 ! call dwgtxt( p % description(i), car) 
                  p % content(i) = character_to_real( car ) 
                  encontrado  = .true. 
                  
                  
                  exit 
             endif      
       endif 
  enddo 
   
 ! call disfin 
 ! call Plot_template( id_g3, esquema ) 
   
 end subroutine    
        
 




end module 
