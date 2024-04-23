!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*  Module to set globally the compiler-dependent kind (and thereby numerical           *
!*  precision) of floating point numbers and potentially other variable types used      *
!*  througout this program.                                                             *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2023-09-05                                                      *
!*   -> latest changes: 2023-09-05                                                      *
!*                                                                                      *
!****************************************************************************************
module Mod_kind_param

implicit none
integer,parameter,public :: wp = selected_real_kind(15)     !wp = working precision kind in terms of real kind; 
                                                            !set to 15 digits of floating point precision ("double precision")

!Note: for integer and logical kinds, we use the default kinds unless specified otherwise locally in a scoping unit.

end module Mod_kind_param