subroutine Test_preliminaries 

   use Roots  
   use Sum_series 
   use Matrix_operations
   use Dynamic_allocation
   use Piecewise_functions 
   use Series_expansion 
   use dislin
   
   
   implicit none 

   call Roots_2th 
   
   call Summation_examples
   
   call Examples_matrix_operations
   
   call Matrices_allocation
   
   call Function_examples 
   
   call Expansion_examples
   
   

end subroutine 

