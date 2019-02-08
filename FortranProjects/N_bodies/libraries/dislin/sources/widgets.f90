module widgets


use dislin 
implicit none 
private 
public :: Place_buttons_with_scales, set_of_widgets_scale, set_widgets_scale 
public :: k_set,  i_button

integer :: k_set,  i_button

type widget_button 

   character(len=50) :: label
   integer :: id 

end type 

type widget_scale 
    
   character(len=50) :: label
   real :: value 
   real :: x0 
   real :: xf 
   real :: step 
   integer :: digits 
   integer :: id 

end type 
       
type set_of_widgets_scale

  type (widget_button) :: b 
  type (widget_scale), allocatable :: s(:) 

end type 


type widgets_bs
   type (set_of_widgets_scale), pointer :: Buttons_con_deslizadores(:) => null()
   procedure(), nopass, pointer :: pVisualiza
end type 

type (widgets_bs) :: wg(300)
integer, save :: N_bs = 0 

  

 ! , pVisualizaD


   
contains

  
 !***************************************************************************************************
 !*
 !***************************************************************************************************
 subroutine  Place_buttons_with_scales(ox, oy, sx, sy, ie, w, Visualiza ) 
 
    integer, intent(in) :: ox, oy, sx, sy, ie 
    type (set_of_widgets_scale), target :: w(:) 
    
    interface
      subroutine Visualiza 
      end subroutine 
    end interface  
    

    integer :: oxi, oyj  
    integer :: i,  j, N, Ns   
    
    N = size( w )
    !do i=1, N 
    !  write(*,'(a45, a20, i2)')  " Button label of a set of widget scales :",  w(i) % b % label, i 
    !end do 
      
  
   
   do i=1, N
       
       call swgpos(ox+(i-1)*sx, oy+sy)
       call swgsiz(sx, sy)
       call wgpbut(ie, w(i) % b % label, w(i) % b % id )
       call swgcbk( w(i) % b %id, get_i_button )
       
     
       oxi = ox-110; oyj = oy + 50 
       Ns = size ( w(i) % s ) 
      ! write(*,*) ' button=', i, ' number of wgscl =', Ns 
       do j=1, Ns  
        call swgsiz(100, 50)
        if (oxi>220) then;   oyj = oyj + 50;   oxi = ox-110;  endif  
      
        if (w(i)%s(j)%label=='\')  then 
                                 oyj = oyj + 50 
                                 oxi = ox-110
                                 cycle 
        endif                             
          
        oxi = oxi + 116
        call swgpos(oxi, oyj)
        call swgstp( w(i) % s(j) % step )
        call wgscl(ie, w(i) % s(j) % label, w(i) % s(j) % x0, w(i) % s(j) % xf, w(i) % s(j) % value, & 
                                            w(i) % s(j) % digits, w(i) % s(j) % id )
        call swgatt( w(i) % s(j) % id, 'invisible', 'status')
        call swgcbk( w(i) % s(j) % id , get_widgets_scale)
       enddo  
     
       enddo 
    
  
!    write(*,*) ' ox =', ox, ' oy =', oy, ' sx = ', sx, ' sy =', sy 
     N_bs = N_bs + 1   
!     write(*,*)  ' N_bs = ', N_bs         
          
     !  Buttons_con_Deslizadores => w
      wg(N_bs) % Buttons_con_Deslizadores => w  
      wg(N_bs) % pVisualiza => Visualiza 
   
 
 end subroutine 
 
  
!***********************************************************************************************
!*
!**********************************************************************************************
 subroutine  get_i_button(id) 
 integer, intent(in) :: id 
 
     integer :: i, j, k ! , iloc(1) 
     integer :: N, Ns 
 

!write(*,*) " N_bs = ", N_bs 
!write(*,*) " get_i_button id  = ", id

do k=1, N_bs  
  
  N = size( wg(k) % Buttons_con_deslizadores ) 
  
 ! iloc = maxloc( wg(k) % Buttons_con_deslizadores(:)%b%id, wg(k) % Buttons_con_deslizadores(:)%b%id >= id ) 
 !  i_button = iloc(1) 
 !  write(*,*) '+++', i_button
   i_button = 0   
   do i=1, N
     if (wg(k) % Buttons_con_deslizadores(i) % b % id == id ) then 
            i_button = i 
            exit
    endif          
   enddo
   
   if (i_button /= 0 ) then 
    k_set = k   
 !   write(*,*) 'k_set=', k_set, 'i_button = ', i_button 
    
   
    do i = 1, N 
          Ns = size ( wg(k) % Buttons_con_deslizadores(i) % s )
             
          if (i==i_button) then 
              do j=1, Ns
                if ( wg(k) % Buttons_con_deslizadores(i)%s(j)%label /= '\' )  then 
                     call swgatt( wg(k) % Buttons_con_deslizadores(i) % s(j) % id, 'active', 'status')
                
                endif       
             enddo  
               
         else
             do j=1, Ns  
              if ( wg(k) % Buttons_con_deslizadores(i)%s(j)%label /= '\' )  then 
                call swgatt( wg(k) % Buttons_con_deslizadores(i)%s(j)%id, 'invisible', 'status')
            
              endif   
             enddo    
         endif     
    enddo 
    call wg(k_set) % pVisualiza  
    end if 
   
   
enddo  
 
 
 end subroutine 

!*********************************************************************************************************
!*
!*********************************************************************************************************
subroutine get_widgets_scale(id) 
integer, intent(in) :: id 
    
    integer :: N, Ns 
    integer :: i, j, k 
    
    logical :: found = .false.
  
!write(*,*) " get_widgets_scale id  = ", id
  
 do k=1, N_bs 
 
    N = size( wg(k) % Buttons_con_deslizadores ) 
    
    do i=1, N 
      Ns = size(wg(k) % Buttons_con_deslizadores(i)%s ) 
      
      do j=1, Ns 
          
 !       write(*,*) " widget_scale id  = ", wg(k) % Buttons_con_deslizadores(i) % s(j) % id  
        if  ( wg(k) % Buttons_con_deslizadores(i) % s(j) % id  == id ) then 
        if (wg(k) % Buttons_con_deslizadores(i) % s(j) % label  /= '\')  then
        
          call gwgscl(wg(k) % Buttons_con_deslizadores(i) % s(j)% id, wg(k) % Buttons_con_deslizadores(i) % s(j) % value )
          found = .true. 
          goto 10  
       endif    
       end if 
      enddo 
    enddo  
    
10  if (found)  then 
        call wg(k) % pVisualiza 
        found = .false. 
    end if 
    
    
 enddo    
    
 !call pVisualiza  
      
end subroutine    

!*********************************************************************************************************
!*
!*********************************************************************************************************
subroutine set_widgets_scale(ip) 
integer, intent(in) :: ip 
    
    integer :: N, Ns 
    integer :: i, j, k  
    integer :: id
    real :: value 
  
 ! write(*,*) " get_widgets_scale id  = ", id
  
 do k=1, N_bs 
 
    N = size( wg(k) % Buttons_con_deslizadores ) 
    
    do i=1, N 
      Ns = size(wg(k) % Buttons_con_deslizadores(i)%s ) 
      
      do j=1, Ns 
         id    =  wg(k) % Buttons_con_deslizadores(i) % s(j) % id  
         value =  wg(k) % Buttons_con_deslizadores(i) % s(j) % value
          
     !   if  ( wg(k) % Buttons_con_deslizadores(i) % s(j) % id  == id ) then 
       if (id  /= 0)  then
          write(*,*) " widget_scale id  = ", id, value 
          call swgscl( id, value)
            
       endif    
     !  end if 
      enddo 
    enddo  
  !  call wg(k) % pVisualiza  
    
 enddo    
    
 !call pVisualiza  
      
end subroutine    






 
 end module 