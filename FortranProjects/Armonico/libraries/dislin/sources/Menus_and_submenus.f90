module Menus_and_submenus 


use dislin 
implicit none 

 integer, parameter :: N_menusmax=10, N_submenusmax = 10 
 character(len=30) :: menu_labels(N_menusmax),  sub_menu_labels(N_menusmax, N_submenusmax) 
 
 
 integer :: id_menu(N_menusmax), id_sub_menu(N_menusmax, N_submenusmax) 
 
 integer :: indice(10000) = 0
 character(len=30) :: etiquetas(10000) = ' '   
 integer :: i_menu, i_submenu  

character(len=30) :: actual_screen

character(len=300) :: ini_filename='sin inicializar' 

logical :: quit = .false.  
logical :: save = .false.  
logical :: open = .false. 

   
contains


!*************************************************************************************
!*
!*************************************************************************************
subroutine set_menus_y_submenus(N, i_screen) 

integer, intent(in) :: N(:), i_screen

     
      integer ::  id
      integer :: id_file,  id_new, id_open,  id_exit ! id_save, id_save_pdf,id_save_txt, 
      integer :: N_menus 
      
      integer :: i ! , j 
     
     call wgpop(i_screen, 'File', id_file) 
    ! call wgapp(id_file, 'New', id_new)   ;          call swgcbk(id_new,    new_file)
     call wgapp(id_file, 'Open', id_open) ;          call swgcbk(id_open,   new_file)
     call wgapp(id_file, 'Save ', id_open) ;          call swgcbk(id_open,   save_file)
     call wgapp(id_file, 'Exit', id_exit) ;          call swgcbk(id_exit,   exit)
    
    
     N_menus = size( N ) 
     do i_menu=1, N_menus
      
        call wgpop(i_screen, menu_labels(i_menu), id_menu(i_menu))
         
        do i=1, N(i_menu) 
     
           call wgapp(id_menu(i_menu),  sub_menu_labels(i_menu, i), id_sub_menu(i_menu,i) )
           etiquetas (id_sub_menu(i_menu,i)) = sub_menu_labels(i_menu, i)         
           call swgcbk( id_sub_menu(i_menu,i), subroutine_sub_menu ) 
                  
        enddo
    enddo        
     
 end subroutine
   
!**************************************************************************************
!
!**************************************************************************************

subroutine subroutine_sub_menu(id)  
 integer, intent(in) :: id  
  
          
   actual_screen = etiquetas( id ) 

   call sendok 
  

end subroutine 
!**************************************************************************************
!
!**************************************************************************************

subroutine new_file(id) 
   integer, intent(in) :: id  
 
  write(*,*) ' Enter new file '
  open = .true. 
 ! call dwgfil('Open initialization file', ini_filename, '*.ini') 

  write(*,*) ' filename = ', ini_filename 

  call sendok
     

end subroutine 
!**************************************************************************************
!
!!**************************************************************************************
!
!subroutine open_file(id) 
!    integer, intent(in) :: id  
!    
! 
!  call dwgfil('Open initialization file', ini_filename, '*.ini') 
!
!  write(*,*) ' filename = ', ini_filename 
!
!  call sendok
!
!
!end subroutine 
!**************************************************************************************
!
!**************************************************************************************

subroutine exit(id) 
    integer, intent(in) :: id  
 
   quit = .true.  
   call sendok

end subroutine 

 !**************************************************************************************
!
!**************************************************************************************

subroutine save_file(id) 
    integer, intent(in) :: id  
 
  ! call dwgfil('Save initialization file', ini_filename, '*.ini') 
  ! call gwgfil(id_filename, ini_file_name)
   save = .true.  
   call sendok

end subroutine 

end module 

