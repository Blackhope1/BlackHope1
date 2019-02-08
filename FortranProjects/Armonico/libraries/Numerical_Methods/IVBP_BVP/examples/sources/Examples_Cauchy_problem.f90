program main

    use mass_spring_damper
    use Waves1d
	use Waves2d
	use Lorenz

    implicit none

	integer:: option
	
	print*, "Please select the problem you want to see"
	print*, "1.- Mass-spring-damper"
	print*, "2.- Lorenz attractor"
	print*, "3.- Waves 1D"
	print*, "4.- Waves 2D"
	print*, "==========================="
	read*, option
	
	if (option==1) then
		call test_mass_spring_damper
	else if (option==2) then  
		call Test_Lorenz_Attractor
	else if (option==3) then 
		call Test_Wave1D
	else if (option==4) then
        call Test_Wave2D
	else
		print*, "Try again!"
	end if
		


end program
