!****************************************************************************
!*                                                                            *
!*  Subroutine to calculate the interval limits of an univariate function,  *
!*  such that the interval bounds bracket a root (zero) for input in a      *
!*  root bracketing solver (e.g. Brent's method).                           *
!*  Here we look inward into a given interval and try to find one or        *
!*  several limits that bracket roots.                                      *
!*  adapted and modified from Numerical Recipes by Andi Zuend, 2011.        *
!*                                                                          *
!****************************************************************************

! "fx" is here a user-provided external object (an external function fx(x) with input arguments x and return function value fx)

subroutine zerobracket_inwards(fx, xlower, xupper, ntry, geomscal, nb, xblow, xbup, success) 

use Mod_kind_param, only : wp

implicit none
!..
!interface variables:
real(wp),external :: fx                                     !user-provided univariate external function fx(x)
real(wp),intent(inout) :: xlower, xupper                    !user supplied lower and upper bounds for the bracketing search (adjust to huge (neg. / pos.) values if no limits are needed)
integer,intent(in) :: ntry                                  !number of interval subdivisions (an input value) (= max. number of roots found this way)
logical,intent(in) :: geomscal                              !set to true to use a geometric mean for scaling for chopping the interval between xlower and xupper.
integer,intent(out) :: nb                                   !number of root brackets found
real(wp),dimension(1:ntry+1),intent(out) :: xblow, xbup     !arrays storing the bounds for the nb roots
logical,intent(out) :: success
!..
!local variables:
integer :: i
real(wp),dimension(1:ntry+1) :: fval, x
real(wp),parameter :: deps = 2.0_wp*epsilon(1.0_wp), sqrt_huge = sqrt(huge(1.0_wp))
real(wp) :: dx, xtemp
!...............................

if (ntry > 0) then
    !check the overall intervall limits:
    if (abs(xlower-xupper) < deps*min(abs(xupper), 1.0_wp)) then
        !!$OMP CRITICAL
        !write(*,*) "Interval in zerobracket_inwards is chosen too small (xlower == xupper)."
        !!$OMP end CRITICAL
        success = .false.
        return
    else if (xlower > xupper) then !switch bounds
        xtemp = xlower
        xlower = xupper
        xupper = xtemp
    endif
    !initialize:
    xblow = xlower-abs(xlower)-99999.9_wp
    xbup =  xlower-abs(xlower)-99999.9_wp
    fval = 0.0_wp
    if (geomscal) then !use logarithmic (base-10) interval chopping
        if (xlower < 0.0_wp .AND. xupper > 0.0_wp) then !cannot use geometric scaling
            dx = (xupper-xlower)/real(ntry, kind=wp) !step in interval search
        endif
    else
        dx = (xupper-xlower)/real(ntry, kind=wp) !step in interval search
    endif
    x(1) = xlower
    do i = 2,ntry
        if (geomscal) then !use geometric mean for scaling
            x(i) = sqrt(x(i-1)*xupper)
        else !use linear step
            x(i) = x(i-1)+dx
        endif
    enddo
    x(ntry+1) = xupper
    fval(1) = fx(x(1))
    nb = 0
    do i = 2,ntry+1
        if (x(i) > xupper+deps) then
            exit
        endif
        fval(i) = fx(x(i)) !evaluate external function fx at point x(i)
        !check for numerical issues (overflow potential):
        if (abs(fval(i)) > sqrt_huge ) then
            fval(i) = sign(sqrt_huge , fval(i))
        endif
        if (fval(i)*fval(i-1) < 0.0_wp) then  ! => zero or pole is bracketed
            nb = nb+1 !count zeros
            xblow(nb) = x(i-1)
            xbup(nb) = x(i)
        endif
    enddo !i
    if (nb > 0) then
        success = .true.
    else
        success = .false.
    endif
else
    success = .false.
endif

end subroutine zerobracket_inwards
