!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module containing numerical transformation functions and subroutines, such as for  *
!*   the treatment of exponential function calls while avoiding numerical overflow.     *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2022-11-02                                                      *
!*   -> latest changes: 2023-03-17                                                      *
!*                                                                                      *
!*   :: License ::                                                                      *
!*   This program is free software: you can redistribute it and/or modify it under the  *
!*   terms of the GNU General Public License as published by the Free Software          *
!*   Foundation, either version 3 of the License, or (at your option) any later         *
!*   version.                                                                           *
!*   The AIOMFAC model code is distributed in the hope that it will be useful, but      *
!*   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or      *
!*   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      *
!*   details.                                                                           *
!*   You should have received a copy of the GNU General Public License along with this  *
!*   program. If not, see <http://www.gnu.org/licenses/>.                               *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  pure elemental function safe_exp                                                *
!*   -  pure elemental function soft_bounds                                             *
!*   -  pure elemental function sigmoidal_map                                           *
!*   -  pure elemental function inverse_sigmoidal_map                                   *  
!*   -  pure elemental function extended_sigmoidal_map                                  *
!*                                                                                      *
!****************************************************************************************
module ModNumericalTransformations

use ModSystemProp, only : Mmass

implicit none
private                     !set default as private for subroutines and module variables
integer,parameter :: wp = kind(Mmass(1))                              !set real kind

!public procedures from this module:
public  ::  safe_exp
public  ::  soft_bounds
public  ::  sigmoidal_map
public  ::  inverse_sigmoidal_map
public  ::  extended_sigmoidal_map

contains

    !********************************************************************************************
    !*   :: Purpose ::                                                                          *
    !*  Elemental function to exponentiate an argument while considering potential need for     *
    !*  truncating the log-scale (x) values in a responsive way (not simply hard bounds), such  *
    !*  that their use in the exp(x) function will not lead to a numerical overflow or          *
    !*  underflow exception.                                                                    *
    !*                                                                                          *
    !*   :: Author & Copyright ::                                                               *
    !*   Andi Zuend,                                                                            *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                              *
    !*                                                                                          *
    !*   -> created:        2022-11-02                                                          *
    !*   -> latest changes: 2022-11-04                                                          *
    !*                                                                                          *
    !********************************************************************************************
    pure elemental function safe_exp(x, ln_limit) result(y)
    
    implicit none
    !interface arguments:
    real(wp),intent(in) :: x            !input; considered a natural log-scale value    
    real(wp),intent(in) :: ln_limit     !imposed hard upper limit of natural log-scale values; e.g. set as 0.48_wp*log(huge(1.0_wp))          
    real(wp)            :: y            !output; y = exp(x) usually, but potentially truncated in x to meet set bounds
    !.............................
    
    !apply soft bounds if necessary:
    y = soft_bounds(x, -ln_limit, ln_limit) 
    !now apply natural exponential function:    
    y = exp(y)                          
    
    end function safe_exp
    !------------------------------------------------------------------------------------------------------------
    
    
    !********************************************************************************************
    !*   :: Purpose ::                                                                          *
    !*  Elemental function to constrain an input value x to fall within provided lower and      *
    !*  upper limits. If x exceeds one of the bounds, the value of x will be adjusted in a      *
    !*  responsive manner, such that subsequent exceeding, yet different values of x would lead *
    !*  to a different output (not simply setting hard bounds). Hard limits are only enforced   *
    !*  in cases in which the "softness" of the applied bounds is exceeded by extreme x values. *
    !*                                                                                          *
    !*   :: Author & Copyright ::                                                               *
    !*   Andi Zuend,                                                                            *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                              *
    !*                                                                                          *
    !*   -> created:        2022-11-04                                                          *
    !*   -> latest changes: 2022-11-04                                                          *
    !*                                                                                          *
    !********************************************************************************************
    pure elemental function soft_bounds(x, l_lim, u_lim) result(y)
    
    implicit none
    !interface arguments:
    real(wp),intent(in) :: x            !input; considered a natural log-scale value    
    real(wp),intent(in) :: l_lim        !imposed hard lower limit
    real(wp),intent(in) :: u_lim        !imposed hard upper limit
    real(wp)            :: y            !output; x value, but potentially truncated to meet set bound constraints
    !local variables:
    real(wp),parameter  :: scaler = 1.0E-2_wp 
    real(wp)            :: local_low_lim, local_upper_lim, margin, trunc
    !.............................
    
    !(1) determine margin of x values in which soft bounds apply:
    margin = 1.0E-2_wp*abs(u_lim - l_lim)
    local_low_lim = l_lim + margin
    local_upper_lim = u_lim - margin
    
    !(2) apply soft bounds:
    if (x > local_upper_lim) then
        trunc = scaler*log(abs(x))
        y = local_upper_lim + trunc     !apply a scaled limit
        if (y > u_lim) then             !very rare case; apply hard bounds
            y = u_lim    
        endif
    else if (x < local_low_lim) then
        trunc = scaler*log(abs(x))
        y = local_low_lim - trunc
        if (y < l_lim) then
            y = l_lim    
        endif 
    else !no bound constraint necessary
        y = x
    endif                         
    
    end function soft_bounds
    !------------------------------------------------------------------------------------------------------------
    
    
    !********************************************************************************************
    !*   :: Purpose ::                                                                          *
    !*  A simple sigmoidal mapping function according to a slightly modified version of         *
    !*  Eq. (2.1) from Yun (2008), doi:10.1017/S1446181108000060                                *
    !*  This sigmoidal function is only valid for inputs x in interval [0, 1] and any exponent  *
    !*  r > 0; the modifications include the addition of parameters to set lower and upper      *
    !*  bounds of the output value space (y-axis range).                                        *
    !*                                                                                          *
    !*   :: Author & Copyright ::                                                               *
    !*   Andi Zuend,                                                                            *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                              *
    !*                                                                                          *
    !*   -> created:        2022-10-25                                                          *
    !*   -> latest changes: 2022-11-26                                                          *
    !*                                                                                          *
    !********************************************************************************************
    pure elemental function sigmoidal_map(x, r, l, u) result(y)
    
    implicit none
    !interface arguments:
    real(wp),intent(in) :: x            !(x-axis) value > 0.0; i.e. on positive real line domain              
    real(wp),intent(in) :: r            !the sigmoidal function order (exponent of x); r > 0; e.g. r = 3.0          
    real(wp),intent(in) :: l, u         !lower and upper asymptotic bounds on the y-value range (default: 0.0, 1.0)   
    real(wp)            :: y            !output; mapped value to [l, u] interval
    !.............................
    
    y = l + (u - l) * x**r / ( x**r + (abs(1.0_wp - x))**r )
    
    end function sigmoidal_map
    !------------------------------------------------------------------------------------------------------------
    
    
    !------------------------------------------------------------------------------------------------------------
    !* inverse of the modified simple sigmoidal function according to Yun (2008)                *
    !* here the function is only valid for inputs y in interval [l, u] and any exponent r > 0.
    pure elemental function inverse_sigmoidal_map(y, r, l, u) result(x)
    
    implicit none
    !interface arguments:
    real(wp),intent(in) :: y            !(x-axis) value > 0.0; i.e. on positive real line domain              
    real(wp),intent(in) :: r            !a sigmoidal function order (exponent of x); r > 0;
    real(wp),intent(in) :: l, u         !lower and upper asymptotic bounds on the y-value range
    real(wp)            :: x            !output; mapped value to [l, u] interval
    !local variables:
    real(wp)            :: ymod
    !.............................
    
    ymod = (y - l) / (u - l)
    ymod = max( min(ymod, 1.0_wp), 0.0_wp )     !to make sure that numerical rounding does not lead to ymod violating the unit interval bounds
    !for inverse, run the modified y value over the unit interval but with 1.0/r as parameter:
    x = sigmoidal_map(ymod, (1.0_wp/r), 0.0_wp, 1.0_wp)     
    
    end function inverse_sigmoidal_map
    !------------------------------------------------------------------------------------------------------------
    
    
    !------------------------------------------------------------------------------------------------------------
    !* extended sigmoidal function according to a modified version of Eq. (4.1) from Yun (2008), doi:10.1017/S1446181108000060 *
    pure elemental function extended_sigmoidal_map(x, r, l, u) result(y)
    
    implicit none
    !interface arguments:
    real(wp),intent(in) :: x            !(x-axis) value > 0.0; i.e. on positive real line domain (and < 9E18 because of floor function)           
    real(wp),intent(in) :: r            !a sigmoidal function order (exponent of x); r > 0;
    real(wp),intent(in) :: l, u         !lower and upper asymptotic bounds on the y-value range
    real(wp)            :: y            !output; mapped value to [l, u] interval
    !local variables:
    integer,parameter   :: i_long = selected_int_kind(18)
    real(wp)            :: xmod
    !.............................

    xmod = abs( x - 2.0_wp*floor(0.5_wp*(x + 1.0_wp), kind = i_long) )
    y = sigmoidal_map(xmod, r, l, u)
    
    end function extended_sigmoidal_map
    !------------------------------------------------------------------------------------------------------------
    
    
end module ModNumericalTransformations