function glomin ( a, b, c, m, machep, e, t, f, x )

!*****************************************************************************80
!
!! GLOMIN seeks a global minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    This function assumes that F(X) is twice continuously differentiable
!    over [A,B] and that F''(X) <= M for all X in [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(wp) A, B, the endpoints of the interval.
!    It must be the case that A < B.
!
!    Input, real(wp) C, an initial guess for the global
!    minimizer.  If no good guess is known, C = A or B is acceptable.
!
!    Input, real(wp) M, the bound on the second derivative.
!
!    Input, real(wp) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real(wp) E, a positive tolerance, a bound for the
!    absolute error in the evaluation of F(X) for any X in [A,B].
!
!    Input, real(wp) T, a positive error tolerance.
!
!    Input, external real(wp) F, the name of a user-supplied
!    function, of the form "function F ( X )", which evaluates the
!    function whose global minimum is being sought.
!
!    Output, real(wp) X, the estimated value of the abscissa
!    for which F attains its global minimum value in [A,B].
!
!    Output, real(wp) GLOMIN, the value F(X).
!
use Mod_kind_param, only : wp
implicit none

real(wp) a
real(wp) a0
real(wp) a2
real(wp) a3
real(wp) b
real(wp) c
real(wp) d0
real(wp) d1
real(wp) d2
real(wp) e
real(wp) f
real(wp) glomin
real(wp) h
integer k
real(wp) m
real(wp) m2
real(wp) machep
real(wp) p
real(wp) q
real(wp) qs
real(wp) r
real(wp) s
real(wp) sc
real(wp) t
real(wp) x
real(wp) y
real(wp) y0
real(wp) y1
real(wp) y2
real(wp) y3
real(wp) yb
real(wp) z0
real(wp) z1
real(wp) z2

a0 = b
x = a0
a2 = a
y0 = f ( b )
yb = y0
y2 = f ( a )
y = y2

if ( y0 < y ) then
    y = y0
else
    x = a
end if

if ( m <= 0.0E0_wp .or. b <= a ) then
    glomin = y
    return
end if

m2 = 0.5E0_wp * ( 1.0E0_wp + 16.0E0_wp * machep ) * m

if ( c <= a .or. b <= c ) then
    sc = 0.5E0_wp * ( a + b )
else
    sc = c
end if

y1 = f ( sc )
k = 3
d0 = a2 - sc
h = 9.0E0_wp / 11.0E0_wp

if ( y1 < y ) then
    x = sc
    y = y1
end if

do

    d1 = a2 - a0
    d2 = sc - a0
    z2 = b - a2
    z0 = y2 - y1
    z1 = y2 - y0
    r = d1 * d1 * z0 - d0 * d0 * z1
    p = r
    qs = 2.0E0_wp * ( d0 * z1 - d1 * z0 )
    q = qs

    if ( k < 1000000 .or. y2 <= y ) then

        do

            if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) < &
                z2 * m2 * r * ( z2 * q - r ) ) then
                a3 = a2 + r / q
                y3 = f ( a3 )

                if ( y3 < y ) then
                    x = a3
                    y = y3
                end if
            end if

            k = mod ( 1611 * k, 1048576 )
            q = 1.0E0_wp
            r = ( b - a ) * 0.00001E0_wp * real ( k, kind=wp )

            if ( z2 <= r ) then
                exit
            end if

        end do

    else

        k = mod ( 1611 * k, 1048576 )
        q = 1.0E0_wp
        r = ( b - a ) * 0.00001E0_wp * k

        do while ( r < z2 )

            if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) < &
                z2 * m2 * r * ( z2 * q - r ) ) then
                a3 = a2 + r / q
                y3 = f ( a3 )

                if ( y3 < y ) then
                    x = a3
                    y = y3
                end if
            end if

            k = mod ( 1611 * k, 1048576 )
            q = 1.0E0_wp
            r = ( b - a ) * 0.00001E0_wp *k

        end do

    end if

    r = m2 * d0 * d1 * d2
    s = sqrt ( ( ( y2 - y ) + t ) / m2 )
    h = 0.5E0_wp * ( 1.0E0_wp + h )
    p = h * ( p + 2.0E0_wp * r * s )
    q = q + 0.5E0_wp * qs
    r = - 0.5E0_wp * ( d0 + ( z0 + 2.01E0_wp * e ) / ( d0 * m2 ) )

    if ( r < s .or. d0 < 0.0E0_wp ) then
        r = a2 + s
    else
        r = a2 + r
    end if

    if ( 0.0E0_wp < p * q ) then
        a3 = a2 + p / q
    else
        a3 = r
    end if

    do

        a3 = max ( a3, r )

        if ( b <= a3 ) then
            a3 = b
            y3 = yb
        else
            y3 = f ( a3 )
        end if

        if ( y3 < y ) then
            x = a3
            y = y3
        end if

        d0 = a3 - a2

        if ( a3 <= r ) then
            exit
        end if

        p = 2.0E0_wp * ( y2 - y3 ) / ( m * d0 )

        if ( ( 1.0E0_wp + 9.0E0_wp * machep ) * d0 <= abs ( p ) ) then
            exit
        end if

        if ( 0.5E0_wp * m2 * ( d0 * d0 + p * p ) <= &
            ( y2 - y ) + ( y3 - y ) + 2.0E0_wp * t ) then
            exit
        end if

        a3 = 0.5E0_wp * ( a2 + a3 )
        h = 0.9E0_wp * h

    end do

    if ( b <= a3 ) then
        exit
    end if

    a0 = sc
    sc = a2
    a2 = a3
    y0 = y1
    y1 = y2
    y2 = y3

end do

glomin = y

return
end
    
    
function local_min ( a, b, eps, t, f, x )

!*****************************************************************************80
!
!! LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much slower
!    than that for a Fibonacci search.  If F has a continuous second
!    derivative which is positive at the minimum (which is not at A or
!    B), then convergence is superlinear, and usually of the order of
!    about 1.324....
!
!    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
!    F is never evaluated at two points closer than TOL.
!
!    If F is a unimodal function and the computed values of F are always
!    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
!    LOCAL_MIN approximates the abscissa of the global minimum of F on the
!    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.
!
!    If F is not unimodal, then LOCAL_MIN may approximate a local, but
!    perhaps non-global, minimum to the same accuracy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(wp) A, B, the endpoints of the interval.
!
!    Input, real(wp) EPS, a positive relative error tolerance.
!    EPS should be no smaller than twice the relative machine precision,
!    and preferably not much less than the square root of the relative
!    machine precision.
!
!    Input, real(wp) T, a positive absolute error tolerance.
!
!    Input, external real(wp) F, the name of a user-supplied
!    function, of the form "function F ( X )", which evaluates the
!    function whose local minimum is being sought.
!
!    Output, real(wp) X, the estimated value of an abscissa
!    for which F attains a local minimum value in [A,B].
!
!    Output, real(wp) LOCAL_MIN, the value F(X).
!
use Mod_kind_param, only : wp

implicit none

real(wp) a
real(wp) b
real(wp) c
real(wp) d
real(wp) e
real(wp) eps
real(wp) f
real(wp) fu
real(wp) fv
real(wp) fw
real(wp) fx
real(wp) local_min
real(wp) m
real(wp) p
real(wp) q
real(wp) r
real(wp) sa
real(wp) sb
real(wp) t
real(wp) t2
real(wp) tol
real(wp) u
real(wp) v
real(wp) w
real(wp) x
!
!  C is the square of the inverse of the golden ratio.
!
c = 0.5E0_wp * ( 3.0E0_wp - sqrt ( 5.0E0_wp ) )

sa = a
sb = b
x = sa + c * ( b - a )
w = x
v = w
e = 0.0E0_wp
fx = f ( x )
fw = fx
fv = fw

do

    m = 0.5E0_wp * ( sa + sb )
    tol = eps * abs ( x ) + t
    t2 = 2.0E0_wp * tol
    !
    !  Check the stopping criterion.
    !
    if ( abs ( x - m ) <= t2 - 0.5E0_wp * ( sb - sa ) ) then
        exit
    end if
    !
    !  Fit a parabola.
    !
    r = 0.0E0_wp
    q = r
    p = q

    if ( tol < abs ( e ) ) then

        r = ( x - w ) * ( fx - fv )
        q = ( x - v ) * ( fx - fw )
        p = ( x - v ) * q - ( x - w ) * r
        q = 2.0E0_wp * ( q - r )

        if ( 0.0E0_wp < q ) then
            p = - p
        end if

        q = abs ( q )

        r = e
        e = d

    end if

    if ( abs ( p ) < abs ( 0.5E0_wp * q * r ) .and. &
        q * ( sa - x ) < p .and. &
        p < q * ( sb - x ) ) then
        !
        !  Take the parabolic interpolation step.
        !
        d = p / q
        u = x + d
        !
        !  F must not be evaluated too close to A or B.
        !
        if ( ( u - sa ) < t2 .or. ( sb - u ) < t2 ) then

            if ( x < m ) then
                d = tol
            else
                d = - tol
            end if

        end if
        !
        !  A golden-section step.
        !
    else

        if ( x < m ) then
            e = sb - x
        else
            e = a - x
        end if

        d = c * e

    end if
    !
    !  F must not be evaluated too close to X.
    !
    if ( tol <= abs ( d ) ) then
        u = x + d
    else if ( 0.0E0_wp < d ) then
        u = x + tol
    else
        u = x - tol
    end if

    fu = f ( u )
    !
    !  Update A, B, V, W, and X.
    !
    if ( fu <= fx ) then

        if ( u < x ) then
            sb = x
        else
            sa = x
        end if

        v = w
        fv = fw
        w = x
        fw = fx
        x = u
        fx = fu

    else

        if ( u < x ) then
            sa = u
        else
            sb = u
        end if

        if ( fu <= fw .or. w == x ) then
            v = w
            fv = fw
            w = u
            fw = fu
        else if ( fu <= fv .or. v == x .or. v== w ) then
            v = u
            fv = fu
        end if

    end if

end do

local_min = fx

return
end
    

subroutine local_min_rc ( a, b, arg, status, value )

!*****************************************************************************80
!
!! LOCAL_MIN_RC seeks a minimizer of a scalar function of a scalar variable.
!
!  Discussion:
!
!    This routine seeks an approximation to the point where a function
!    F attains a minimum on the interval (A,B).
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much
!    slower than that for a Fibonacci search.  If F has a continuous
!    second derivative which is positive at the minimum (which is not
!    at A or B), then convergence is superlinear, and usually of the
!    order of about 1.324...
!
!    The routine is a revised version of the Brent local minimization
!    algorithm, using reverse communication.
!
!    It is worth stating explicitly that this routine will NOT be
!    able to detect a minimizer that occurs at either initial endpoint
!    A or B.  If this is a concern to the user, then the user must
!    either ensure that the initial interval is larger, or to check
!    the function value at the returned minimizer against the values
!    at either endpoint.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters
!
!    Input/output, real(wp) A, B.  On input, the left and right
!    endpoints of the initial interval.  On output, the lower and upper
!    bounds for an interval containing the minimizer.  It is required
!    that A < B.
!
!    Output, real(wp) ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS zero, ARG is the routine's
!    estimate for the function minimizer.
!
!    Input/output, integer STATUS, used to communicate between
!    the user and the routine.  The user only sets STATUS to zero on the first
!    call, to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated minimizer.
!
!    Input, real(wp) VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
!  Local parameters:
!
!    C is the squared inverse of the golden ratio.
!
!    EPS is the square root of the relative machine precision.
!
use Mod_kind_param, only : wp

implicit none

real(wp) a
real(wp) arg
real(wp) b
real(wp), save :: c
real(wp), save :: d
real(wp), save :: e
real(wp), save :: eps
real(wp), save :: fu
real(wp), save :: fv
real(wp), save :: fw
real(wp), save :: fx
real(wp), save :: midpoint
real(wp), save :: p
real(wp), save :: q
real(wp), save :: r
integer status
real(wp), save :: tol
real(wp), save :: tol1
real(wp), save :: tol2
real(wp), save :: u
real(wp), save :: v
real(wp) value
real(wp), save :: w
real(wp), save :: x
!
!  STATUS (INPUT) = 0, startup.
!
if ( status == 0 ) then

    if ( b <= a ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LOCAL_MIN_RC - Fatal error!'
        write ( *, '(a)' ) '  A < B is required, but'
        write ( *, '(a,g14.6)' ) '  A = ', a
        write ( *, '(a,g14.6)' ) '  B = ', b
        status = -1
        stop
    end if

    c = 0.5E0_wp * ( 3.0E0_wp - sqrt ( 5.0E0_wp ) )

    eps = sqrt ( epsilon ( eps ) )
    tol = epsilon ( tol )

    v = a + c * ( b - a )
    w = v
    x = v
    e = 0.0E0_wp

    status = 1
    arg = x

    return
    !
    !  STATUS (INPUT) = 1, return with initial function value of FX.
    !
else if ( status == 1 ) then

    fx = value
    fv = fx
    fw = fx
    !
    !  STATUS (INPUT) = 2 or more, update the data.
    !
else if ( 2 <= status ) then

    fu = value

    if ( fu <= fx ) then

        if ( x <= u ) then
            a = x
        else
            b = x
        end if

        v = w
        fv = fw
        w = x
        fw = fx
        x = u
        fx = fu

    else

        if ( u < x ) then
            a = u
        else
            b = u
        end if

        if ( fu <= fw .or. w == x ) then
            v = w
            fv = fw
            w = u
            fw = fu
        else if ( fu <= fv .or. v == x .or. v == w ) then
            v = u
            fv = fu
        end if

    end if

end if
!
!  Take the next step.
!
midpoint = 0.5E0_wp * ( a + b )
tol1 = eps * abs ( x ) + tol / 3.0E0_wp
tol2 = 2.0E0_wp * tol1
!
!  If the stopping criterion is satisfied, we can exit.
!
if ( abs ( x - midpoint ) <= ( tol2 - 0.5E0_wp * ( b - a ) ) ) then
    status = 0
    return
end if
!
!  Is golden-section necessary?
!
if ( abs ( e ) <= tol1 ) then
    if ( midpoint <= x ) then
        e = a - x
    else
        e = b - x
    end if

    d = c * e
    !
    !  Consider fitting a parabola.
    !
else

    r = ( x - w ) * ( fx - fv )
    q = ( x - v ) * ( fx - fw )
    p = ( x - v ) * q - ( x - w ) * r
    q = 2.0E0_wp * ( q - r )
    if ( 0.0E0_wp < q ) then
        p = - p
    end if
    q = abs ( q )
    r = e
    e = d
    !
    !  Choose a golden-section step if the parabola is not advised.
    !
    if ( &
        ( abs ( 0.5E0_wp * q * r ) <= abs ( p ) ) .or. &
        ( p <= q * ( a - x ) ) .or. &
        ( q * ( b - x ) <= p ) ) then

        if ( midpoint <= x ) then
            e = a - x
        else
            e = b - x
        end if

        d = c * e
        !
        !  Choose a parabolic interpolation step.
        !
    else

        d = p / q
        u = x + d

        if ( ( u - a ) < tol2 ) then
            d = sign ( tol1, midpoint - x )
        end if

        if ( ( b - u ) < tol2 ) then
            d = sign ( tol1, midpoint - x )
        end if

    end if

end if
!
!  F must not be evaluated too close to X.
!
if ( tol1 <= abs ( d ) ) then
    u = x + d
end if

if ( abs ( d ) < tol1 ) then
    u = x + sign ( tol1, d )
end if
!
!  Request value of F(U).
!
arg = u
status = status + 1

return
end
    
    
subroutine timestampbrent ( )

!*****************************************************************************80
!
!! timestampbrent prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
implicit none

character (len=8) :: ampm
integer d
integer h
integer m
integer mm
character (len=9), parameter, dimension(12) :: month = [ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' ]
integer n
integer s
integer values(8)
integer y

call date_and_time ( values = values )

y = values(1)
m = values(2)
d = values(3)
h = values(5)
n = values(6)
s = values(7)
mm = values(8)

if ( h < 12 ) then
    ampm = 'AM'
else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
        ampm = 'Noon'
    else
        ampm = 'PM'
    end if
else
    h = h - 12
    if ( h < 12 ) then
        ampm = 'PM'
    else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
            ampm = 'Midnight'
        else
            ampm = 'AM'
        end if
    end if
end if

write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

return
end


function brentzero ( a, b, machep, t, f )

!*****************************************************************************80
!
!! brentzero seeks the root of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the zero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(wp) A, B, the endpoints of the change of sign interval.
!
!    Input, real(wp) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real(wp) T, a positive error tolerance.
!
!    Input, external real(wp) F, the name of a user-supplied
!    function, of the form "function F ( X )", which evaluates the
!    function whose zero is being sought.
!
!    Output, real(wp) zero, the estimated value of a zero of
!    the function F.
!
use Mod_kind_param, only : wp

implicit none

real(wp) a
real(wp) b
real(wp) c
real(wp) d
real(wp) e
real(wp) f
real(wp) fa
real(wp) fb
real(wp) fc
real(wp) m
real(wp) machep
real(wp) p
real(wp) q
real(wp) r
real(wp) s
real(wp) sa
real(wp) sb
real(wp) t
real(wp) tol
real(wp) brentzero
!
!  Make local copies of A and B.
!
sa = a
sb = b
fa = f ( sa )
fb = f ( sb )

c = sa
fc = fa
e = sb - sa
d = e

do
    if ( abs ( fc ) < abs ( fb ) ) then
        sa = sb
        sb = c
        c = sa
        fa = fb
        fb = fc
        fc = fa
    end if

    tol = 2.0E0_wp * machep * abs ( sb ) + t
    m = 0.5E0_wp * ( c - sb )

    if ( abs ( m ) <= tol .or. fb == 0.0E0_wp ) then
        exit
    end if

    if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then
        e = m
        d = e
    else
        s = fb / fa
        if ( sa == c ) then
            p = 2.0E0_wp * m * s
            q = 1.0E0_wp - s
        else
            q = fa / fc
            r = fb / fc
            p = s * ( 2.0E0_wp * m * a * ( q - r ) - ( sb - sa ) * ( r - 1.0E0_wp ) )
            q = ( q - 1.0E0_wp ) * ( r - 1.0E0_wp ) * ( s - 1.0E0_wp )
        end if

        if ( 0.0E0_wp < p ) then
            q = - q
        else
            p = - p
        end if

        s = e
        e = d

        if ( 2.0E0_wp * p < 3.0E0_wp * m * q - abs ( tol * q ) .and. &
            p < abs ( 0.5E0_wp * s * q ) ) then
            d = p / q
        else
            e = m
            d = e
        end if

    end if

    sa = sb
    fa = fb

    if ( tol < abs ( d ) ) then
        sb = sb + d
    else if ( 0.0E0_wp < m ) then
        sb = sb + tol
    else
        sb = sb - tol
    end if

    fb = f ( sb )

    if ( ( 0.0E0_wp < fb .and. 0.0E0_wp < fc ) .or. &
        ( fb <= 0.0E0_wp .and. fc <= 0.0E0_wp ) ) then
        c = sa
        fc = fa
        e = sb - sa
        d = e
    end if
enddo

brentzero = sb

return
end


subroutine zero_rc ( a, b, t, arg, status, value )

!*****************************************************************************80
!
!! ZERO_RC seeks the root of a function F(X) using reverse communication.
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the zero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    The routine is a revised version of the Brent zero finder
!    algorithm, using reverse communication.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(wp) A, B, the endpoints of the change of sign interval.
!
!    Input, real(wp) T, a positive error tolerance.
!
!    Output, real(wp) ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS zero, ARG is the routine's
!    estimate for the function's zero.
!
!    Input/output, integer STATUS, used to communicate between
!    the user and the routine.  The user only sets STATUS to zero on the first
!    call, to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated zero
!
!    Input, real(wp) VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
use Mod_kind_param, only : wp

implicit none

real(wp) a
real(wp) arg
real(wp) b
real(wp), save :: c
real(wp), save :: d
real(wp), save :: e
real(wp), save :: fa
real(wp), save :: fb
real(wp), save :: fc
real(wp) m
real(wp), save :: machep
real(wp) p
real(wp) q
real(wp) r
real(wp) s
real(wp), save :: sa
real(wp), save :: sb
integer status
real(wp) t
real(wp) tol
real(wp) value
!
!  Input STATUS = 0.
!  Initialize, request F(A).
!
if ( status == 0 ) then

    machep = epsilon ( a )

    sa = a
    sb = b
    e = sb - sa
    d = e

    status = 1
    arg = a
    return
    !
    !  Input STATUS = 1.
    !  Receive F(A), request F(B).
    !
else if ( status == 1 ) then

    fa = value

    status = 2
    arg = sb
    return
    !
    !  Input STATUS = 2
    !  Receive F(B).
    !
else if ( status == 2 ) then

    fb = value

    if ( 0.0E0_wp < fa * fb ) then
        status = -1
        return
    end if

    c = sa
    fc = fa

else

    fb = value

    if ( ( 0.0E0_wp < fb .and. 0.0E0_wp < fc ) .or. &
        ( fb <= 0.0E0_wp .and. fc <= 0.0E0_wp ) ) then
        c = sa
        fc = fa
        e = sb - sa
        d = e
    end if

end if
!
!  Compute the next point at which a function value is requested.
!
if ( abs ( fc ) < abs ( fb ) ) then

    sa = sb
    sb = c
    c = sa
    fa = fb
    fb = fc
    fc = fa

end if

tol = 2.0E0_wp * machep * abs ( sb ) + t
m = 0.5E0_wp * ( c - sb )

if ( abs ( m ) <= tol .or. fb == 0.0E0_wp ) then
    status = 0
    arg = sb
    return
end if

if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

    e = m
    d = e

else

    s = fb / fa

    if ( sa == c ) then

        p = 2.0E0_wp * m * s
        q = 1.0E0_wp - s

    else

        q = fa / fc
        r = fb / fc
        p = s * ( 2.0E0_wp * m * a * ( q - r ) - ( sb - sa ) * ( r - 1.0E0_wp ) )
        q = ( q - 1.0E0_wp ) * ( r - 1.0E0_wp ) * ( s - 1.0E0_wp )

    end if

    if ( 0.0E0_wp < p ) then
        q = - q
    else
        p = - p
    end if

    s = e
    e = d

    if ( 2.0E0_wp * p < 3.0E0_wp * m * q - abs ( tol * q ) .and. &
        p < abs ( 0.5E0_wp * s * q ) ) then
        d = p / q
    else
        e = m
        d = e
    end if

end if

sa = sb
fa = fb

if ( tol < abs ( d ) ) then
    sb = sb + d
else if ( 0.0E0_wp < m ) then
    sb = sb + tol
else
    sb = sb - tol
end if

arg = sb
status = status + 1

return
end
    