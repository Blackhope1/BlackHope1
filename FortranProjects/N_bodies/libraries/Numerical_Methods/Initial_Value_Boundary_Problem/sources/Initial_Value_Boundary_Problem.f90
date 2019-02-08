
!***********************************************************************
! It integrates in time the Initial value boundary problem.
! Given the differential operator and the boundary equations,
! the discrete solution is calculated.
! author: Juan A Hernandez, juanantonio.hernandez@upm.es 
!***********************************************************************
module Initial_Value_Boundary_Problem

use Initial_Value_Boundary_Problem1D
use Initial_Value_Boundary_Problem2D
use Utilities
use Temporal_Schemes
use Finite_differences


implicit none   


interface  Initial_Value_Boundary_ProblemS
                                         module procedure IVBP1D, IVBP1D_system, IVBP2D,  IVBP2D_system
end interface



end module 