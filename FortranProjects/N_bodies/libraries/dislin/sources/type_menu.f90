
module type_menu 


abstract interface 
          subroutine screen(ini, quit, save, open, id1, id2, id3, id4, id5 )
            character(len=*), intent(in) :: ini 
            logical :: quit, save, open  
            integer, intent(in) :: id1, id2, id3, id4, id5 
           
          end subroutine  
 end interface 

type pointer_procedure_interface 
    procedure(screen), pointer, nopass :: ps
    character(len=30) description 
    character(len=300) :: ini_filename = ' '
end type 

type menus 
  type (pointer_procedure_interface), allocatable :: sub_menu(:) 
  character(len=30) :: description 
   
end type 


end module 