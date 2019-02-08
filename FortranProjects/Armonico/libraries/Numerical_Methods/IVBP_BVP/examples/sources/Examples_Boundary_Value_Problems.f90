program High_order

    use Boundary_value_problems
    use example_rod

  implicit none 
  integer:: option

  	
	print*, "Please select the problem you want to see"
	print*, "1.- Example: rod with axial load"
	print*, "2.- Test linear_BVP_1D"
	print*, "3.- Test non_linear_BVP_1D_1"
	print*, "4.- Test non_linear_BVP_1D_2"
	print*, "5.- Test non_linear_BVP_1D_3"
	print*, "6.- Test linear_BVP_2D"
	print*, "7.- Test non_linear_BVP_2D"
	print*, "==========================="
	read*, option
	
	if (option==1) then
		call test_rod
	else if (option==2) then  
		call Test_BVP1D
	else if (option==3) then 
		call Test_non_linear_BVP_1
	else if (option==4) then
        call Test_non_linear_BVP_2
	else if (option==5) then
        call Test_non_linear_BVP_3
	else if (option==6) then
        call Test_BVP2D
	else if (option==7) then
        call Test_non_linear_BVP_4
	else
		print*, "Try again!"
	end if
		

end 










 

