!********************************************************************************
! 
!********************************************************************************
module  GUI_dislin
     
    use Menus_and_submenus 
    use Fonts_and_formats
    use type_menu 
    implicit none 
    
      
 contains 
!**************************************************************************************
!
!**************************************************************************************
 subroutine GUI_Aplication( Titulo, Menu, imenu, isubmenu ) ! 
  character(len=*), intent(in) :: Titulo
  type (menus) :: Menu(:) 
  integer, intent(in) :: imenu, isubmenu 
   
    
  
   integer :: id_screen !, id_bottoms, id_graf  
   integer :: nx, ny   
   integer :: id_control, id_datos, id_graph1, id_graph2, id_graph3
   
! ip, 
   
   integer ::N_submenus(size(Menu)) 
   integer :: N_menus  
   character(len=300) :: ini 
   
   integer :: i ! , ipa
  
   N_menus = size(Menu) 
   do i=1, N_menus 
      N_submenus(i) = size( Menu(i)%sub_menu ) 
   enddo    
      
   actual_screen = Menu(imenu) % sub_menu(isubmenu)% description 
 
 
 10  call Set_fonts(nx, ny) 

     call swgtit (Titulo)
      
     menu_labels =  ' ' 
     
     do i_menu=1, N_menus
      
       menu_labels(i_menu) = Menu(i_menu) % description 
       
       do i_submenu= 1, N_submenus(i_menu)  
       
         sub_menu_labels(i_menu, i_submenu) =  menu(i_menu) % sub_menu(i_submenu)%description 
         
       enddo  
     enddo  
       
    call wgini('form',id_screen)  
         
       call set_menus_y_submenus( N_submenus, id_screen )
      	
      do i_menu = 1, N_menus; do i_submenu=1, N_submenus(i_menu) 
      
         if  (actual_screen == Menu(i_menu) % sub_menu(i_submenu) % description ) then 
                
                 if (trim(ini_filename) == 'sin inicializar') then 
                      
                 else
                      Menu(i_menu) % sub_menu(i_submenu) % ini_filename = ini_filename
                 endif   
                 ini = Menu(i_menu) % sub_menu(i_submenu) % ini_filename
              !   write(*,*) ' desde interfaz para entrar en sub submenu =', trim(ini) 
                                 
                 call Divide_screen(id_screen, id_control, id_datos, id_graph1, id_graph2, id_graph3) 
                 call Menu(i_menu) % sub_menu(i_submenu) % ps( ini, quit, save, open, id_control, id_datos, id_graph1, id_graph2, id_graph3) 
                    
         endif 
      enddo;  enddo 
    
    call wgfin 
    goto 10 

end  subroutine 

!********************************************************************************************************************
!* Divide screen  1280x1024 
!********************************************************************************************************************
subroutine Divide_screen(id_screen, id_control, id_datos, id_graph1, id_graph2, id_graph3)              
      integer, intent(in) :: id_screen 
      integer, intent(out) :: id_control, id_datos, id_graph1, id_graph2, id_graph3 
           
   integer :: sx, sy, sx1, sy1, sx2, sy2,  ox1, oy1, ox2, oy2, ox3, oy3, sx3, sy3  ! ox, oy,
   real :: fc, f 
   
   f = 1.0
   fc = 1.0
    
  
   sx = int(380*fc)
   sy = int(370*fc) ! 470
   call swgsiz(sx, sy) 
   call swgpos(0,0) 
   call wgbas(id_screen, 'form', id_control ) 
    
   
   call swgsiz(sx, sy) 
   call swgpos(0,sy) 
   call wgbas(id_screen, 'form', id_datos ) 
    
   
      
   sx1 = 885 ! 885*f 
   !sy1 = 800 ! 427*f
   sy1 = 700 ! 427*f 
   ox1 = sx; oy1 = 0 
   
   sx2 = 885 ! 885*f 
   sy2 = 700 ! 427*f 
   ox2 = sx; oy2 = 0 
   
   sx3 = 885 ! 885*f 
   sy3 = 200 ! 427*f 
   ox3 = sx; oy3 = 700 
   
   call swgsiz(sx1, sy1) 
   call swgpos(ox1, oy1) 
   call wgdraw(id_screen, id_graph1)
   
   call swgsiz(sx2, sy2) 
   call swgpos(ox2, oy2) 
   call wgdraw(id_screen, id_graph2)
   
   call swgsiz(sx3, sy3) 
   call swgpos(ox3, oy3) 
   call wgdraw(id_screen, id_graph3)
   
   
  
   
   !sx2 = 890
   !sy2 = 950 ! *f !- sy1
   !ox2 = sx; !oy = sy1 ! 472*f
   !oy2 = 0 
   
   !!f = 0.7 
   !!sx = 885 * f 
   !!sy1 = 427 * f 
   !!sy2 = 950 * f - sy1
   !!ox = 380;
   
  ! sx = 0.8 * sx 
   
   !call swgsiz(sx1, sy1+sy2) 
   !call swgpos(ox,0) 
   !call wgdraw(id_screen, id_graph3)    
   !
   !call swgsiz(sx1, sy1) 
   !call swgpos(ox,oy) 
   !call wgdraw(id_screen, id_graph1)    
        
  
        
  
        
end subroutine 

 
end module 
































































































