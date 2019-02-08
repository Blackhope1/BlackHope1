
module Fonts_and_formats 

use dislin 
!USE IFPORT
use widgets 

implicit none 

contains 

!***************************************************************************************
!*
!***************************************************************************************
subroutine Set_fonts(nx, ny) 
      integer, intent(out) :: nx, ny 
 
     
   integer :: nw, nh 
   real :: R, G, B 
   character(len=20) :: resolucion_pantalla 

  call getscr(nw, nh)  
  write(*,*)  "nw = ", nw, " nh = ", nh 

  
  
  if (nw==1024 .and. nh==768) then 

                                  resolucion_pantalla = '1024x768'
                                  call swgfnt('Arial Bold', 18) 
                                !  call swgfnt('Times New Roman Bold', 16)

  elseif (nw==800 .and. nh==600) then 

                                  resolucion_pantalla = '800x600'
                                  call swgfnt('Arial', 10) 
                                !  call swgfnt('Times New Roman Bold', 12)
  elseif (nw==1024 .and. nh==600) then 

                                  resolucion_pantalla = '1024x600'
                                  call swgfnt('Arial', 10) 
                                                             

  elseif (nw==640 .and. nh==480) then 

                                  resolucion_pantalla = '640x480'
                                  call swgfnt('Arial', 8)   
                                 ! call swgfnt('Times New Roman Bold', 10)
  elseif (nw==1280 .and. nh==1024) then 

                                  resolucion_pantalla = '1280x1024'
                                  call swgfnt('Arial',  14) ! 16) ! 

  elseif (nw==1280 .and. nh==800) then ! WARNING 

                                  resolucion_pantalla = '1280x800'
                                  call swgfnt('Arial', 14) 

   elseif (nw==1680 .and. nh==1050) then ! WARNING 

                                  resolucion_pantalla = '1680x1050'
                                  call swgfnt('Arial', 14) 

  elseif (nw==1366 .and. nh==768) then ! WARNING 

                                  resolucion_pantalla = '1366x768'
                                  call swgfnt('Arial', 14) 

  elseif (nw==1600 .and. nh==1200) then ! WARNING 

                                  resolucion_pantalla = '1600x1200'
                                  call swgfnt('Arial', 14) 
                                  
  elseif (nw==1920 .and. nh==967) then ! WARNING 

                                  resolucion_pantalla = '1920x967'
                                  call swgfnt('Arial', 14)
                                  
  elseif (nw==1882 .and. nh==916) then ! WARNING 

                                  resolucion_pantalla = '1882x916'
                                  call swgfnt('Arial', 14)
  elseif (nw==1920 .and. nh==1080) then ! WARNING 

                                  resolucion_pantalla = '920x1080'
                                  call swgfnt('Arial', 14) 
                                  
  elseif (nw==3840 .and. nh==2160) then ! WARNING 

                                  resolucion_pantalla = '3840x2160'
                                  call swgfnt('Arial', 12)   
                                  
  elseif (nw==2560 .and. nh==1440) then ! WARNING 

                                  resolucion_pantalla = '2560x1440'
                                  call swgfnt('Arial', 14)     
                                  

  else
   
  
         !write(debug_unit,*) ' Resolucion pantalla no conocida  =', nw, nh !bbbbbbbb
         !write(*,*) ' Resolucion pantalla no conocida  =', nw, nh !bbbbbbbb
         resolucion_pantalla = 'wwwwxhhhh'
         call swgfnt('Arial', 12)    
        ! stop
      !   nw = 3840
       !  nh = 2160
          

  endif 


  nx = nw 
  ny = nh 

  ! *** colores 
      R = 212./256; G = 208./256; B = 200./256  
      call swgclr(R, G, B, 'back') 
      call swgclr(R, G, B, 'scroll')
      call swgclr(R, G, B, 'ltext')
   !   call swgclr(R, G, B, 'fore')
      
       R = 212./256; G = 208./256; B = 200./256  
      call swgclr(R, G, B, 'back') 
      
      
       write(*,*) ' Pantalla =', resolucion_pantalla 
      
     call swgpop('NOHELP')
     call swgpop('NOOK')
     call swgpop('NOQUIT')
     
     

     call swgopt ('center', 'position')
 !    call winapp('windows')
     

    call swgjus('center', 'button')
! call swgopt ('change', 'callback') ! es mi modificacion helmut dislin para diccionarios. 
                                     ! sin embargo, cuando se utiliza wgtxt  para meter txt con swgtxt da problemas porque llama a la 
                                     ! subrutina que tiene vinculada mediante call back sin nosotros quererlo. 
                                     ! dejarlo por defecto en 'return' salvo necesidad justificada. 
                                     
                                     
    call swgopt ('standard', 'position') 
    


   ! call grey 
    
end subroutine 


!*****************************************************************************
!*
!*****************************************************************************
subroutine grey 

   integer :: ic 
   real :: R, G, B 
   integer :: i_grey  

    ic = 123
  !  call setind(ic, 112./256., 108/256., 100./256.)
  !  call setclr(ic)


     R = 212./256; G = 208./256; B = 200./256  
     R = 1; G = 0; B = 0 
     i_grey = indrgb( R, G, B )


end subroutine 

!**************************************************************************************************************
!*
!**************************************************************************************************************  
 function  Inside( x , x1, x2 ) 
 integer, intent(in) :: x(:), x1(:), x2(:)   
 logical :: Inside 
  
    if ( all( (x-x1)>0 ) .and. all( (x-x2)<0 ) ) then 
                Inside = .true. 
    else                            
                Inside = .false. 
    endif                  
        
  
  
  end function   
!**************************************************************************************************************
!*
!**************************************************************************************************************  
 function  Distancia( x , y ) 
 integer, intent(in) :: x(:), y(:)  
 real :: Distancia 
  
    if (any((x-y)<0) ) then 
                 Distancia = -1 
    else                            
                 Distancia =  maxval( x-y )
    endif                  
        
  
  
  end function   
!**************************************************************************************************************
!*
!**************************************************************************************************************  
 subroutine  Decrementa( N ) 
 integer, intent(inout) :: N 
  
   if (N>1) then 
      N = N - 1 
   else 
     N = 0 
  endif    
        
  
  
  end subroutine   
 
!**************************************************************************************************************
!*
!**************************************************************************************************************  
 subroutine  Incrementa( N, N_max ) 
 integer, intent(in) :: N_max 
 integer, intent(inout) :: N 
  
   if (N<N_max) then 
      N = N + 1 
   else 
     N = N_max 
  endif    
        
  
  
  end subroutine   




!********************************************************************************************
function integer_to_character(i)
!********************************************************************************************
 integer, intent(in) :: i 
 character(len=4) :: integer_to_character
 
       
    if (i<10) then  
         write(integer_to_character, '(i1)') i 
    elseif (i<100) then    
         write(integer_to_character, '(i2)') i
    elseif (i<1000) then    
         write(integer_to_character, '(i3)') i   
    elseif (i<10000) then    
         write(integer_to_character, '(i4)') i 
    endif           
   

end function


!********************************************************************************************
function real_to_character(r)
!********************************************************************************************
 real, intent(in) :: r 
 character(len=5) :: real_to_character
 
       
    write(real_to_character, '(f5.1)') r
   

end function


!********************************************************************************************
function character_to_integer(c)
!********************************************************************************************
 character(len=2), intent(in) :: c 
 integer :: character_to_integer
 
    open(8,  file='conv') 
    write(8,'(a2)') c 
    rewind(8) 
    read(8, '(i2)') character_to_integer
    close(8)    
   
   

end function

!********************************************************************************************
function character_to_real(c)
!********************************************************************************************
 character(len=*), intent(in) :: c 
 real :: character_to_real
 
    open(8,  file='conv') 
    write(8,* ) c 
    rewind(8) 
    read(8, *) character_to_real
    close(8)    
   
   

end function

!******************************************************************************************
subroutine  Busca_subcadena( Vector, car, ip) 
!******************************************************************************************

character(len=*), intent(in) :: Vector(:), car
integer, intent(out) :: ip   


   integer :: i, j(1), is(size(Vector)), N  
   character(len=10) carU  

N = size(Vector) 

!write(*,*) 'car = ', car 
!do i=1, N
!write(*,*) Vector(i) 
!enddo 

carU = car 
call upstr (carU)

if (car==' ' ) then

  do i=1, N
  
  if (trim(Vector(i)) == '' ) then 
      j(1) = i 
      exit 
  endif
  enddo  
 
else 

do i=1, N

! write(*,*) index(Vector(i), car)
 
 is(i) = max( index(trim(Vector(i)), trim(car)), index(trim(Vector(i)), trim(carU)) ) 

enddo 
j = minloc( is(1:N), is>0 ) 
endif 

if (j(1) == 0 ) then 
    write(*,*) ' Busca subcadena no encontrada '
    ip = -1 
     
else 
    ip = j(1) 
    
endif     
 



end subroutine 

!******************************************************************************************
subroutine  dislin_list( Vector, list ) 
!******************************************************************************************

character(len=*), intent(in) :: Vector(:)
character(len=*), intent(out) :: list  


   integer :: i, N, M, i2   
   
   N = size( Vector ) 
   list = trim( Vector(1) )
  
     
   i2 = 0 
   
   do i=2, N 
        M = len( trim(Vector(i-1)) ) 
       i2 = i2 + M + 1
       if (trim(Vector(i)) /= " ") then 
         list = list(1:i2) // "|" // trim(Vector(i))
       end if 
   end do 
   
 !  write(*,*) Vector
 !  write(*,*) trim(list) 
   
!   read(*,*) 


end subroutine 




!****************************************************************************************
!*
!***************************************************************************************
subroutine mensaje(  c, r, ix, iy ) 
character(len=*), intent(in) :: c 
real, intent(in) :: r 
integer, intent(in) :: ix, iy 

    integer ::  ox, ndigits ! n,
 
      call messag( c, ix, iy )
      ox = nlmess(c) 
      if (r/=-123456789) then 
      
         if (r == int(r)) then 
                                          ndigits = 0 
         elseif (10*r == int(10*r)) then  
                                          ndigits = 1 
         else 
                                          ndigits = 2                       
         endif            
                                         
         call number ( r, ndigits, ox + ix, iy ) 
      endif  
 

end subroutine  


!***********************************************************************************************************************
!*
!**********************************************************************************************************************
subroutine send_email( recipient, sender, password, file, subject ) 
character(len=*), intent(in) :: recipient, sender, password, file, subject 


character(len=400) :: email, t1, t2 
!integer :: ierr 

t1 = 'sendemail -t '// trim(recipient) // ' -f '// trim(sender) // ' -s smtp.gmail.com:587 '//' -u '// trim(subject)
t2 = ' -xu '// trim(sender)//' -xp '// trim(password)// ' -o  message-file='//trim(file) //' -o message-charset="LATIN1" '
email = trim(t1) // trim (t2)   
 
!ierr = system( email )
 
end subroutine 
 

!***********************************************************************
!*
!***********************************************************************
function email_sender(Instalacion) 
      character(len=*), intent(in) :: Instalacion 
      character(len=150) :: email_sender 
      
      
      email_sender = 'IntelliGlass.'//trim(Instalacion)//'@gmail.com'


end function 

!***********************************************************************
!*
!***********************************************************************
function filtro_primer_orden ( valor_anterior_filtrado, valor_actual )
    real, intent(in) :: valor_anterior_filtrado, valor_actual
    real :: filtro_primer_orden
    
    real :: A, t_muestreo = 60 , t_filtro = 1200
   
    A = exp ( - t_muestreo / t_filtro )
    
    filtro_primer_orden = A * valor_anterior_filtrado + ( 1 - A ) * valor_actual 

end function 

 

end module 
