!****************************************************************************
!*																		    *
!*  Subroutine to calculate the interval limits of an univariate function,  *
!*  such that the interval bounds bracket a root (zero) for input in a      *
!*  root bracketing solver (e.g. Brent's method).                           *
!*  Here we look inward into a given interval and try to find one or        *
!*  several limits that bracket roots.                                      *
!*  adapted and modified from Numerical Recipes by Andi Zuend, 2011.        *
!*                                                                          *
!****************************************************************************

SUBROUTINE zerobracket_inwards(fx, xlower, xupper, ntry, geomscal, nb, xblow, xbup, success) 

IMPLICIT NONE
!..
!interface variables:
REAL(8),EXTERNAL :: fx                   !user-provided univariate external function fx(x)
REAL(8),INTENT(INOUT) :: xlower, xupper  !user supplied lower and upper bounds for the bracketing search (adjust to huge (neg. / pos.) values if no limits are needed)
INTEGER(4),INTENT(IN) :: ntry            !number of interval subdivisions (an input value) (= max. number of roots found this way)
LOGICAL(4),INTENT(IN) :: geomscal        !set to true to use a geometric mean for scaling for chopping the interval between xlower and xupper.
INTEGER(4),INTENT(OUT) :: nb             !number of root brackets found
REAL(8),DIMENSION(1:ntry+1),INTENT(OUT) :: xblow, xbup !arrays storing the bounds for the nb roots
LOGICAL(4),INTENT(OUT) :: success
!..
!local variables:
INTEGER(4) :: i
REAL(8),DIMENSION(1:ntry+1) :: fval, x
REAL(8),PARAMETER :: DEPS = 2.0D0*EPSILON(1.0D0)
REAL(8) :: dx, xtemp
!...............................

IF (ntry > 0) THEN
    !check the overall intervall limits:
    IF (ABS(xlower-xupper) < DEPS*MIN(ABS(xupper),1.0D0)) THEN
        !!$OMP CRITICAL
        !WRITE(*,*) "Interval in zerobracket_inwards is chosen too small (xlower == xupper)."
        !!$OMP END CRITICAL
        success = .false.
        RETURN
    ELSE IF (xlower > xupper) THEN !switch bounds
        xtemp = xlower
        xlower = xupper
        xupper = xtemp
    ENDIF
    !initialize:
    xblow = xlower-ABS(xlower)-99999.9D0
    xbup =  xlower-ABS(xlower)-99999.9D0
    fval = 0.0D0
    IF (geomscal) THEN !use logarithmic (base-10) interval chopping
        IF (xlower < 0.0D0 .AND. xupper > 0.0D0) THEN !cannot use geometric scaling
            dx = (xupper-xlower)/REAL(ntry, KIND=8) !step in interval search
        ENDIF
    ELSE
        dx = (xupper-xlower)/REAL(ntry, KIND=8) !step in interval search
    ENDIF
    x(1) = xlower
    DO i = 2,ntry
        IF (geomscal) THEN !use geometric mean for scaling
            x(i) = SQRT(x(i-1)*xupper)
        ELSE !use linear step
            x(i) = x(i-1)+dx
        ENDIF
    ENDDO
    x(ntry+1) = xupper
    fval(1) = fx(x(1))
    nb = 0
    DO i = 2,ntry+1
        IF (x(i) > xupper+DEPS) THEN
            EXIT
        ENDIF
        fval(i) = fx(x(i)) !evaluate external function fx at point x(i)
        !check for numerical issues (overflow potential):
        IF (ABS(fval(i)) > 1.0D100) THEN
            fval(i) = SIGN(1.0D100, fval(i))
        ENDIF
        IF (fval(i)*fval(i-1) < 0.0D0) THEN  ! => zero or pole is bracketed
            nb = nb+1 !count zeros
            xblow(nb) = x(i-1)
            xbup(nb) = x(i)
        ENDIF
    ENDDO !i
    IF (nb > 0) THEN
        success = .true.
    ELSE
        success = .false.
    ENDIF
ELSE
    success = .false.
ENDIF

END SUBROUTINE zerobracket_inwards
