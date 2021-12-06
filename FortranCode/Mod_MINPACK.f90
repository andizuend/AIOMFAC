!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module containing the MINPACK methods for solving non-linear equations by Powell's *
!*   hybrid method, the Levenberg-Marquardt method and related procedures for finite-   *
!*   difference gradient calculations and other MINPACK utility functions .             *
!*                                                                                      *
!*   f90 code from https://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.html,     *
!*   See licence information and original netlib source in the comments further below.  *
!*                                                                                      *
!*   :: Authors & Copyright ::                                                          *
!*   Jorge More, Danny Sorenson, Burton Garbow, Kenneth Hillstrom                       *
!*   (User Guide for MINPACK-1, Argonne National Laboratory, Argonne, Illinois.)        *
!*                                                                                      *
!*   Minpack Copyright Notice (1999) University of Chicago.  All rights reserved.       *
!*   "This product includes software developed by the University of Chicago, as         *
!*   Operator of Argonne National Laboratory."                                          *
!*                                                                                      *
!*   Fortran 90 code and some modifications by John Burkardt,                           *
!*   Scientific Computing, Florida State Univ. (https://people.sc.fsu.edu/~jburkardt)   *
!*                                                                                      *
!*   Some modifications of Fortran 90 code version by Andi Zuend,                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*   Those modifications include use of set real(kind=wp) kind instead of real(8) and   *
!*   passing of subroutine/function optional parameters to hybrd1/hybrd via a derived   *
!*   type.                                                                              *
!*                                                                                      *
!*                                                                                      *
!**************************************************************************************** 
module Mod_MINPACK

implicit none

private  !default
integer,parameter :: wp = KIND(1.0D0)

!public module procedures (add others if needed):
public :: hybrd1, lmdif1

    contains
    
    
    !-----------------------------------------------------------

    !"This product includes software developed by the
    !   University of Chicago, as Operator of Argonne National
    !   Laboratory."

    
    !==========================================================================
    !DISCLAIMER as retrieved from http://www.netlib.org/minpack/  
    !-----------------------------------------------------------
    !
    !Minpack Copyright Notice (1999) University of Chicago.  All rights reserved
    !
    !Redistribution and use in source and binary forms, with or
    !without modification, are permitted provided that the
    !following conditions are met:
    !
    !1. Redistributions of source code must retain the above
    !copyright notice, this list of conditions and the following
    !disclaimer.
    !
    !2. Redistributions in binary form must reproduce the above
    !copyright notice, this list of conditions and the following
    !disclaimer in the documentation and/or other materials
    !provided with the distribution.
    !
    !3. The end-user documentation included with the
    !redistribution, if any, must include the following
    !acknowledgment:
    !
    !   "This product includes software developed by the
    !   University of Chicago, as Operator of Argonne National
    !   Laboratory."
    !
    !Alternately, this acknowledgment may appear in the software
    !itself, if and wherever such third-party acknowledgments
    !normally appear.
    !
    !4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
    !WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
    !UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
    !THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
    !IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
    !OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
    !OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
    !OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
    !USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
    !THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
    !DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
    !UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
    !BE CORRECTED.
    !
    !5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
    !HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
    !ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
    !INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
    !ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
    !PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
    !SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
    !(INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
    !EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
    !POSSIBILITY OF SUCH LOSS OR DAMAGES.
    
    !==========================================================================


    subroutine chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err )
    !*****************************************************************************80
    !
    !! chkder checks the gradients of M functions of N variables.
    !
    !  Discussion:
    !
    !    chkder checks the gradients of M nonlinear functions in N variables,
    !    evaluated at a point X, for consistency with the functions themselves.
    !
    !    The user calls chkder twice, first with MODE = 1 and then with MODE = 2.
    !
    !    MODE = 1.
    !      On input,
    !        X contains the point of evaluation.
    !      On output,
    !        XP is set to a neighboring point.
    !
    !    Now the user must evaluate the function and gradients at X, and the
    !    function at XP.  Then the subroutine is called again:
    !
    !    MODE = 2.
    !      On input,
    !        FVEC contains the function values at X,
    !        FJAC contains the function gradients at X.
    !        FVECP contains the functions evaluated at XP.
    !      On output,
    !        ERR contains measures of correctness of the respective gradients.
    !
    !    The subroutine does not perform reliably if cancellation or
    !    rounding errors cause a severe loss of significance in the
    !    evaluation of a function.  Therefore, none of the components
    !    of X should be unusually small (in particular, zero) or any
    !    other value which may cause loss of significance.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, integer :: M, is the number of functions.
    !
    !    Input, integer :: N, is the number of variables.
    !
    !    Input, real(wp) :: X(N), the point at which the jacobian is to be
    !    evaluated.
    !
    !    Input, real(wp) :: FVEC(M), is used only when MODE = 2.
    !    In that case, it should contain the function values at X.
    !
    !    Input, real(wp) :: FJAC(LDFJAC,N), an M by N array.  When MODE = 2,
    !    FJAC(I,J) should contain the value of dF(I)/dX(J).
    !
    !    Input, integer :: LDFJAC, the leading dimension of the array FJAC.
    !    LDFJAC must be at least M.
    !
    !    Output, real(wp) :: XP(N), on output with MODE = 1, is a neighboring
    !    point of X, at which the function is to be evaluated.
    !
    !    Input, real(wp) :: FVECP(M), on input with MODE = 2, is the function
    !    value at XP.
    !
    !    Input, integer :: MODE, should be set to 1 on the first call, and
    !    2 on the second.
    !
    !    Output, real(wp) :: ERR(M).  On output when MODE = 2, ERR contains
    !    measures of correctness of the respective gradients.  If there is no
    !    severe loss of significance, then if ERR(I):
    !      = 1.0E0_wp, the I-th gradient is correct,
    !      = 0.0E0_wp, the I-th gradient is incorrect.
    !      > 0.5E0_wp, the I-th gradient is probably correct.
    !      < 0.5E0_wp, the I-th gradient is probably incorrect.
    !
    implicit none

    integer :: ldfjac
    integer :: m
    integer :: n

    real(wp) :: eps
    real(wp) :: epsf
    real(wp) :: epslog
    real(wp) :: epsmch
    real(wp) :: err(m)
    real(wp) :: fjac(ldfjac,n)
    real(wp) :: fvec(m)
    real(wp) :: fvecp(m)
    integer :: i
    integer :: j
    integer :: mode
    real(wp) :: temp
    real(wp) :: x(n)
    real(wp) :: xp(n)

    epsmch = epsilon ( epsmch )
    eps = sqrt ( epsmch )
    !
    !  MODE = 1.
    !
    if ( mode == 1 ) then

        do j = 1, n
            temp = eps * abs( x(j) )
            if ( abs(temp) < epsmch) then
                temp = eps
            end if
            xp(j) = x(j) + temp
        end do
        !
        !  MODE = 2.
        !
    else if ( mode == 2 ) then

        epsf = 100.0E0_wp * epsmch
        epslog = log10 ( eps )

        err = 0.0E0_wp

        do j = 1, n
            temp = abs( x(j) )
            if (abs(temp) < epsmch) then
                temp = 1.0E0_wp
            end if
            err(1:m) = err(1:m) + temp * fjac(1:m,j)
        end do

        do i = 1, m

            temp = 1.0E0_wp

            if ( fvec(i) /= 0.0E0_wp .and. fvecp(i) /= 0.0E0_wp .and. &
                abs( fvecp(i)-fvec(i)) >= epsf * abs( fvec(i) ) ) then
            temp = eps * abs( (fvecp(i)-fvec(i)) / eps - err(i) ) &
                / ( abs( fvec(i) ) + abs( fvecp(i) ) )
            end if

            err(i) = 1.0E0_wp

            if ( epsmch < temp .and. temp < eps ) then
                err(i) = ( log10 ( temp ) - epslog ) / epslog
            end if

            if ( eps <= temp ) then
                err(i) = 0.0E0_wp
            end if

        end do

    end if

    return
    end
    !-----------------------------------------------------------

    
    pure subroutine d_swap( x, y )
    !*****************************************************************************80
    !
    !! d_swap switches two floating point values.
    !
    !
    !  Modified:
    !
    !    01 May 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, real(wp) :: X, Y.  On output, the values of X and
    !    Y have been interchanged.
    !
    implicit none

    real(wp),intent(inout) :: x
    real(wp),intent(inout) :: y
    real(wp) :: z

    z = x
    x = y
    y = z

    end subroutine d_swap
    !-----------------------------------------------------------


    subroutine dogleg( n, r, lr, diag, qtb, delta, x )
    !*****************************************************************************80
    !
    !! dogleg finds the minimizing combination of Gauss-Newton and gradient steps.
    !
    !  Discussion:
    !
    !    Given an M by N matrix A, an N by N nonsingular diagonal
    !    matrix D, an M-vector B, and a positive number DELTA, the
    !    problem is to determine the convex combination X of the
    !    Gauss-Newton and scaled gradient directions that minimizes
    !    (A*X - B) in the least squares sense, subject to the
    !    restriction that the euclidean norm of D*X be at most DELTA.
    !
    !    This subroutine completes the solution of the problem
    !    if it is provided with the necessary information from the
    !    QR factorization of A.  That is, if A = Q*R, where Q has
    !    orthogonal columns and R is an upper triangular matrix,
    !    then dogleg expects the full upper triangle of R and
    !    the first N components of Q'*B.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, integer :: N, the order of the matrix R.
    !
    !    Input, real(wp) :: R(LR), the upper triangular matrix R stored
    !    by rows.
    !
    !    Input, integer :: LR, the size of the R array, which must be no less
    !    than (N*(N+1))/2.
    !
    !    Input, real(wp) :: DIAG(N), the diagonal elements of the matrix D.
    !
    !    Input, real(wp) :: QTB(N), the first N elements of the vector Q'* B.
    !
    !    Input, real(wp) :: DELTA, is a positive upper bound on the
    !    euclidean norm of D*X(1:N).
    !
    !    Output, real(wp) :: X(N), the desired convex combination of the
    !    Gauss-Newton direction and the scaled gradient direction.
    !
    implicit none

    integer :: lr
    integer :: n

    real(wp) :: alpha
    real(wp) :: bnorm
    real(wp) :: delta
    real(wp) :: diag(n)
    real(wp) :: epsmch
    real(wp) :: gnorm
    integer :: i
    integer :: j
    integer :: jj
    integer :: k
    integer :: l
    real(wp) :: qnorm
    real(wp) :: qtb(n)
    real(wp) :: r(lr)
    real(wp) :: sgnorm
    real(wp) :: sum2
    real(wp) :: temp
    real(wp) :: wa1(n)
    real(wp) :: wa2(n)
    real(wp) :: x(n)

    epsmch = epsilon( epsmch )
    !
    !  Calculate the Gauss-Newton direction.
    !
    jj = ( n * ( n + 1 ) ) / 2 + 1

    do k = 1, n
        j = n - k + 1
        jj = jj - k
        l = jj + 1
        sum2 = 0.0E0_wp

        do i = j + 1, n
            sum2 = sum2 + r(l) * x(i)
            l = l + 1
        end do
        !correct if sum2 is too large or too small, which could cause numerical issues below; block introduced by A. Zuend, Oct. 2012.
        if (sum2 > 1.0E100_wp) then
            sum2 = 1.0E100_wp
        else if (sum2 < -1.0E100_wp) then
            sum2 = -1.0E100_wp
        endif

        temp = r(jj)
        if (abs(temp) <= epsmch) then
            l = j
            do i = 1, j
                temp = max ( temp, abs( r(l)) )
                l = l + n - i
            end do
            if ( temp <= epsmch  ) then
                temp = epsmch
            else
                temp = epsmch * temp
            end if
        end if

        if (qtb(j) > 1.0E100_wp*epsmch) then
            qtb(j) = 1.0E100_wp*epsmch
        else if (qtb(j) < -1.0E100_wp*epsmch) then
            qtb(j) = 1.0E100_wp*epsmch
        endif
        x(j) = ( qtb(j) - sum2 ) / temp
    end do
    !
    !  Test whether the Gauss-Newton direction is acceptable.
    !
    wa1(1:n) = 0.0E0_wp
    wa2(1:n) = diag(1:n) * x(1:n)
    qnorm = enorm2 ( n, wa2 ) !enorm ( n, wa2 )

    if ( qnorm <= delta ) then
        return
    end if
    !
    !  The Gauss-Newton direction is not acceptable.
    !  Calculate the scaled gradient direction.
    !
    l = 1
    do j = 1, n
        temp = qtb(j)
        do i = j, n
            wa1(i) = wa1(i) + r(l) * temp
            l = l + 1
        end do
        wa1(j) = wa1(j) / diag(j)
    end do
    !
    !  Calculate the norm of the scaled gradient.
    !  Test for the special case in which the scaled gradient is zero.
    !
    gnorm = enorm2( n, wa1 ) !enorm ( n, wa1 )
    sgnorm = 0.0E0_wp
    alpha = delta / qnorm

    if ( gnorm /= 0.0E0_wp ) then
        !
        !  Calculate the point along the scaled gradient which minimizes the quadratic.
        !
        wa1(1:n) = ( wa1(1:n) / gnorm ) / diag(1:n)

        l = 1
        do j = 1, n
            sum2 = 0.0E0_wp
            do i = j, n
                sum2 = sum2 + r(l) * wa1(i)
                l = l + 1
            end do
            wa2(j) = sum2
        end do

        temp = enorm2( n, wa2 ) !enorm ( n, wa2 )
        sgnorm = ( gnorm / temp ) / temp
        !
        !  Test whether the scaled gradient direction is acceptable.
        !
        alpha = 0.0E0_wp
        !
        !  The scaled gradient direction is not acceptable.
        !  Calculate the point along the dogleg at which the quadratic is minimized.
        !
        if ( sgnorm < delta ) then
            bnorm = enorm2( n, qtb ) !enorm ( n, qtb )
            temp = ( bnorm / gnorm ) * ( bnorm / qnorm ) * ( sgnorm / delta )
            temp = temp - ( delta / qnorm ) * ( sgnorm / delta)**2 &
                + sqrt ( ( temp - ( delta / qnorm ) )**2 &
                + ( 1.0E0_wp - ( delta / qnorm )**2 ) &
                * ( 1.0E0_wp - ( sgnorm / delta )**2 ) )

            alpha = ( ( delta / qnorm ) * ( 1.0E0_wp - ( sgnorm / delta )**2 ) ) / temp
        end if
    end if
    !
    !  Form appropriate convex combination of the Gauss-Newton
    !  direction and the scaled gradient direction.
    !
    temp = ( 1.0E0_wp - alpha ) * min( sgnorm, delta )

    x(1:n) = temp * wa1(1:n) + alpha * x(1:n)

    return
    end
    !-----------------------------------------------------------

    
    pure function enorm(n, x)
    !*****************************************************************************80
    !
    !! enorm computes the Euclidean norm of a vector.
    !
    !  Discussion:
    !
    !    This is an extremely simplified version of the original enorm
    !    routine, which has been renamed to "enorm2".
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, integer :: N, is the length of the vector.
    !
    !    Input, real(wp) :: X(N), the vector whose norm is desired.
    !
    !    Output, real(wp) :: enorm, the Euclidean norm of the vector.
    !
    implicit none

    integer,intent(in) :: n
    real(wp),intent(in) :: x(n)
    real(wp) :: enorm

    if (ANY(abs(x) > 1.0E100_wp)) then  !overflow exception
        enorm = sum(abs(x))
    else
        enorm = sqrt(sum(x**2))
    endif

    end function enorm
    !-----------------------------------------------------------

    
    pure function enorm2(n, x)
    !*****************************************************************************80
    !
    !! enorm2 computes the Euclidean norm of a vector.
    !
    !  Discussion:
    !
    !    This routine was named enorm.  It has been renamed "enorm2",
    !    and a simplified routine has been substituted.
    !
    !    The Euclidean norm is computed by accumulating the sum of
    !    squares in three different sums.  The sums of squares for the
    !    small and large components are scaled so that no overflows
    !    occur.  Non-destructive underflows are permitted.  Underflows
    !    and overflows do not occur in the computation of the unscaled
    !    sum of squares for the intermediate components.
    !
    !    The definitions of small, intermediate and large components
    !    depend on two constants, RDWARF and RGIANT.  The main
    !    restrictions on these constants are that RDWARF**2 not
    !    underflow and RGIANT**2 not overflow.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, integer :: n, is the length of the vector.
    !
    !    Input, real(wp) :: x(n), the vector whose norm is desired.
    !
    !    Output, real(wp) :: enorm2, the Euclidean norm of the vector.
    !
    implicit none

    integer,intent(in) :: n
    real(wp),intent(in) :: x(n)
    real(wp) :: enorm2  !return value

    integer :: i
    real(wp) :: agiant
    real(wp) :: rdwarf
    real(wp) :: rgiant
    real(wp) :: s1
    real(wp) :: s2
    real(wp) :: s3
    real(wp) :: xabs
    real(wp) :: x1max
    real(wp) :: x3max

    rdwarf = sqrt(tiny( rdwarf ))
    rgiant = sqrt(huge( rgiant ))

    s1 = 0.0E0_wp
    s2 = 0.0E0_wp
    s3 = 0.0E0_wp
    x1max = 0.0E0_wp
    x3max = 0.0E0_wp
    agiant = rgiant / real(n, kind=wp)

    do i = 1, n
        xabs = abs( x(i) )
        if ( xabs <= rdwarf ) then
            if ( x3max < xabs ) then
                s3 = 1.0E0_wp + s3 * ( x3max / xabs )**2
                x3max = xabs
            else if ( xabs /= 0.0E0_wp ) then
                s3 = s3 + ( xabs / x3max )**2
            end if
        else if ( agiant <= xabs ) then
            if ( x1max < xabs ) then
                s1 = 1.0E0_wp + s1 * ( x1max / xabs )**2
                x1max = xabs
            else
                s1 = s1 + ( xabs / x1max )**2
            end if
        else
            s2 = s2 + xabs**2
        end if
    end do
    !
    !  Calculation of norm.
    !
    if ( s1 /= 0.0E0_wp ) then
        enorm2 = x1max * sqrt( s1 + ( s2 / x1max ) / x1max )
    else if ( s2 /= 0.0E0_wp ) then
        if ( x3max <= s2 ) then
            enorm2 = sqrt( s2 * ( 1.0E0_wp + ( x3max / s2 ) * ( x3max * s3 ) ) )
        else
            enorm2 = sqrt( x3max * ( ( s2 / x3max ) + ( x3max * s3 ) ) )
        end if
    else
        enorm2 = x3max * sqrt(s3)
    end if

    end function enorm2
    !-----------------------------------------------------------

    
    subroutine fdjac1 ( fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn )
    !*****************************************************************************80
    !
    !! fdjac1 estimates an N by N jacobian matrix using forward differences.
    !
    !  Discussion:
    !
    !    This subroutine computes a forward-difference approximation
    !    to the N by N jacobian matrix associated with a specified
    !    problem of N functions in N variables. If the jacobian has
    !    a banded form, then function evaluations are saved by only
    !    approximating the nonzero terms.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, external fcn, the name of the user-supplied subroutine which
    !    calculates the functions.  The routine should have the form:
    !
    !      subroutine fcn( n, x, fvec, iflag )
    !
    !      integer :: n
    !
    !      real(wp) fvec(n)
    !      integer :: iflag
    !      real(wp) x(n)
    !
    !    The value of IFLAG should not be changed by fcn unless
    !    the user wants to terminate execution of the routine.
    !    In this case set IFLAG to a negative integer.
    !
    !    Input, integer :: N, the number of functions and variables.
    !
    !    Input, real(wp) :: X(N), the point where the jacobian is evaluated.
    !
    !    Input, real(wp) :: FVEC(N), the functions evaluated at X.
    !
    !    Output, real(wp) :: FJAC(LDFJAC,N), the N by N approximate
    !    jacobian matrix.
    !
    !    Input, integer :: LDFJAC, the leading dimension of FJAC, which must
    !    not be less than N.
    !
    !    Output, integer :: IFLAG, is an error flag returned by fcn.  If fcn
    !    returns a nonzero value of IFLAG, then this routine returns immediately
    !    to the calling program, with the value of IFLAG.
    !
    !    Input, integer :: ML, MU, specify the number of subdiagonals and
    !    superdiagonals within the band of the jacobian matrix.  If the
    !    jacobian is not banded, set ML and MU to N-1.
    !
    !    Input, real(wp) :: EPSFCN, is used in determining a suitable step
    !    length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    implicit none

    integer :: ldfjac
    integer :: n

    real(wp) :: eps
    real(wp) :: epsfcn
    real(wp) :: epsmch
    external fcn
    real(wp) :: fjac(ldfjac,n)
    real(wp) :: fvec(n)
    real(wp) :: h
    integer :: i
    integer :: iflag
    integer :: j
    integer :: k
    integer :: ml
    integer :: msum
    integer :: mu
    real(wp) :: temp
    real(wp) :: wa1(n)
    real(wp) :: wa2(n)
    real(wp) :: x(n)

    epsmch = epsilon ( epsmch )

    eps = sqrt ( max ( epsfcn, epsmch ) )
    msum = ml + mu + 1

    !
    !  Computation of dense approximate jacobian.
    !
    if ( n <= msum ) then
        do j = 1, n
            temp = x(j)
            h = eps * abs( temp )
            if ( abs(h) < epsmch) then
                IF (abs(temp) > epsmch) THEN
                    h = sqrt(abs(temp) * epsmch)    !changed by AZ
                ELSE
                    h = sqrt(eps * epsmch) !h = eps !changed by AZ
                ENDIF
            end if

            x(j) = temp + h
            call fcn ( n, x, wa1, iflag )

            if ( iflag < 0 ) then
                exit
            end if

            x(j) = temp
            IF (ANY(fvec(1:n) < -1.0E100_wp)) THEN  !floating point exception handling (by AZ)
                WHERE (fvec(1:n) < -1.0E100_wp)
                    fvec(1:n) = -1.0E100_wp
                ENDWHERE
            ENDIF
            IF (ANY(wa1(1:n) < -1.0E100_wp)) THEN   !floating point exception handling (by AZ)
                WHERE (wa1(1:n) < -1.0E100_wp)
                    wa1(1:n) = -1.0E100_wp
                ENDWHERE
            ENDIF
            fjac(1:n,j) = ( wa1(1:n) - fvec(1:n) ) / h
        end do
    else
        !
        !  Computation of banded approximate jacobian.
        !
        do k = 1, msum
            do j = k, n, msum
                wa2(j) = x(j)
                temp = x(j)
                h = eps * abs( wa2(j) )
                if ( abs(h) < epsmch ) then
                    IF (abs(temp) > epsmch) THEN
                        h = sqrt(abs(temp) * epsmch) !changed by AZ
                    ELSE
                        h = sqrt(eps * epsmch) !h = eps !changed by AZ
                    ENDIF
                end if
                x(j) = wa2(j) + h
            end do

            call fcn ( n, x, wa1, iflag )

            if ( iflag < 0 ) then
                exit
            end if

            do j = k, n, msum
                x(j) = wa2(j)
                temp = x(j)
                h = eps * abs( wa2(j) )
                if ( abs(h) < epsmch ) then
                    IF (abs(temp) > epsmch) THEN
                        h = sqrt(abs(temp) * epsmch) !changed by AZ
                    ELSE
                        h = sqrt(eps * epsmch) !h = eps !changed by AZ
                    ENDIF
                end if
                fjac(1:n,j) = 0.0E0_wp
                do i = 1, n
                    if ( j - mu <= i .and. i <= j + ml ) then
                        fjac(i,j) = ( wa1(i) - fvec(i) ) / h
                    end if
                end do
            end do
        end do
    end if

    return
    end
    !-----------------------------------------------------------

    
    subroutine fdjac2 ( fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn )
    !*****************************************************************************80
    !
    !! fdjac2 estimates an M by N jacobian matrix using forward differences.
    !
    !  Discussion:
    !
    !    This subroutine computes a forward-difference approximation
    !    to the M by N jacobian matrix associated with a specified
    !    problem of M functions in N variables.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, external fcn, the name of the user-supplied subroutine which
    !    calculates the functions.  The routine should have the form:
    !
    !      subroutine fcn ( m, n, x, fvec, iflag )
    !      integer :: n
    !      real(wp) fvec(m)
    !      integer :: iflag
    !      real(wp) x(n)
    !
    !    The value of IFLAG should not be changed by fcn unless
    !    the user wants to terminate execution of the routine.
    !    In this case set IFLAG to a negative integer.
    !
    !    Input, integer :: M, is the number of functions.
    !
    !    Input, integer :: N, is the number of variables.  N must not exceed M.
    !
    !    Input, real(wp) :: X(N), the point where the jacobian is evaluated.
    !
    !    Input, real(wp) :: FVEC(M), the functions evaluated at X.
    !
    !    Output, real(wp) :: FJAC(LDFJAC,N), the M by N approximate
    !    jacobian matrix.
    !
    !    Input, integer :: LDFJAC, the leading dimension of FJAC, which must
    !    not be less than M.
    !
    !    Output, integer :: IFLAG, is an error flag returned by fcn.  If fcn
    !    returns a nonzero value of IFLAG, then this routine returns immediately
    !    to the calling program, with the value of IFLAG.
    !
    !    Input, real(wp) :: EPSFCN, is used in determining a suitable
    !    step length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    implicit none

    integer :: ldfjac
    integer :: m
    integer :: n

    real(wp) :: eps
    real(wp) :: epsfcn
    real(wp) :: epsmch
    external fcn
    real(wp) :: fjac(ldfjac,n)
    real(wp) :: fvec(m)
    real(wp) :: h
    !integer :: i
    integer :: iflag
    integer :: j
    real(wp) :: temp
    real(wp) :: wa(m)
    real(wp) :: x(n)

    epsmch = epsilon ( epsmch )

    eps = sqrt ( max ( epsfcn, epsmch ) )

    do j = 1, n

        temp = x(j)
        h = eps * abs( temp )
        if ( abs(h) < epsmch ) then
            h = eps
        end if

        x(j) = temp + h
        call fcn ( m, n, x, wa, iflag )

        if ( iflag < 0 ) then
            exit
        end if

        x(j) = temp
        fjac(1:m,j) = ( wa(1:m) - fvec(1:m) ) / h

    end do

    return
    end
    !-----------------------------------------------------------


    subroutine hybrd ( fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
        factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf )

    !*****************************************************************************80
    !
    !! hybrd seeks a zero of N nonlinear equations in N variables.
    !
    !  Discussion:
    !
    !    hybrd finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  The user must provide a
    !    subroutine which calculates the functions.  The jacobian is
    !    then calculated by a forward-difference approximation.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, external fcn, the name of the user-supplied subroutine which
    !    calculates the functions.  The routine should have the form:
    !
    !      subroutine fcn(n, x, fvec, iflag)
    !      integer,intent(in)       :: n
    !      real(wp),intent(inout)   :: x(n)
    !      real(wp),intent(out)     :: fvec(n)
    !      integer,intent(out)      :: iflag
    !
    !    The value of IFLAG should not be changed by fcn unless
    !    the user wants to terminate execution of the routine.
    !    In this case set IFLAG to a negative integer.
    !
    !    Input, integer :: N, the number of functions and variables.
    !
    !    Input/output, real(wp) :: X(N).  On input, X must contain an initial
    !    estimate of the solution vector.  On output X contains the final
    !    estimate of the solution vector.
    !
    !    Output, real(wp) :: FVEC(N), the functions evaluated at the output X.
    !
    !    Input, real(wp) :: XTOL.  Termination occurs when the relative error
    !    between two consecutive iterates is at most XTOL.  XTOL should be
    !    nonnegative.
    !
    !    Input, integer :: MAXFEV.  Termination occurs when the number of calls
    !    to fcn is at least MAXFEV by the end of an iteration.
    !
    !    Input, integer :: ML, MU, specify the number of subdiagonals and
    !    superdiagonals within the band of the jacobian matrix.  If the jacobian
    !    is not banded, set ML and MU to at least n - 1.
    !
    !    Input, real(wp) :: EPSFCN, is used in determining a suitable step
    !    length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    !    Input/output, real(wp) :: DIAG(N).  If MODE = 1, then DIAG is set
    !    internally.  If MODE = 2, then DIAG must contain positive entries that
    !    serve as multiplicative scale factors for the variables.
    !
    !    Input, integer :: MODE, scaling option.
    !    1, variables will be scaled internally.
    !    2, scaling is specified by the input DIAG vector.
    !
    !    Input, real(wp) :: FACTOR, determines the initial step bound.  This
    !    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
    !    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
    !    in the interval (0.1, 100) with 100 the recommended value.
    !
    !    Input, integer :: NPRINT, enables controlled printing of iterates if it
    !    is positive.  In this case, fcn is called with IFLAG = 0 at the
    !    beginning of the first iteration and every NPRINT iterations thereafter
    !    and immediately prior to return, with X and FVEC available
    !    for printing.  If NPRINT is not positive, no special calls
    !    of fcn with IFLAG = 0 are made.
    !
    !    Output, integer :: INFO, error flag.  If the user has terminated execution,
    !    INFO is set to the (negative) value of IFLAG. See the description of fcn.
    !    Otherwise, INFO is set as follows:
    !    0, improper input parameters.
    !    1, relative error between two consecutive iterates is at most XTOL.
    !    2, number of calls to fcn has reached or exceeded MAXFEV.
    !    3, XTOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    4, iteration is not making good progress, as measured by the improvement
    !       from the last five jacobian evaluations.
    !    5, iteration is not making good progress, as measured by the improvement
    !       from the last ten iterations.
    !
    !    Output, integer :: NFEV, the number of calls to fcn.
    !
    !    Output, real(wp) :: FJAC(LDFJAC,N), an N by N array which contains
    !    the orthogonal matrix Q produced by the QR factorization of the final
    !    approximate jacobian.
    !
    !    Input, integer :: LDFJAC, the leading dimension of the array FJAC.
    !    LDFJAC must be at least N.
    !
    !    Output, real(wp) :: R(LR), the upper triangular matrix produced by
    !    the QR factorization of the final approximate jacobian, stored rowwise.
    !
    !    Input, integer :: LR, the size of the R array, which must be no less
    !    than (N*(N+1))/2.
    !
    !    Output, real(wp) :: QTF(N), contains the vector Q'*FVEC.
    !
    implicit none

    integer :: ldfjac
    integer :: lr
    integer :: n

    real(wp) :: actred
    real(wp) :: delta
    real(wp) :: diag(n)
    real(wp) :: epsfcn
    real(wp) :: epsmch
    real(wp) :: factor
    external fcn
    real(wp) :: fjac(ldfjac,n)
    real(wp) :: fnorm
    real(wp) :: fnorm1
    real(wp) :: fvec(n)
    integer :: i
    integer :: iflag
    integer :: info
    integer :: iter
    integer :: iwa(1)
    integer :: j
    logical :: jeval
    integer :: l
    integer :: maxfev
    integer :: ml
    integer :: mode
    integer :: msum
    integer :: mu
    integer :: ncfail
    integer :: nslow1
    integer :: nslow2
    integer :: ncsuc
    integer :: nfev
    integer :: nprint
    real(wp) :: pnorm
    real(wp) :: prered
    real(wp) :: qtf(n)
    real(wp) :: r(lr)
    real(wp) :: ratio
    logical :: sing
    real(wp) :: sum2
    real(wp) :: temp
    real(wp) :: wa1(n)
    real(wp) :: wa2(n)
    real(wp) :: wa3(n)
    real(wp) :: wa4(n)
    real(wp) :: x(n)
    real(wp) :: xnorm
    real(wp) :: xtol

    epsmch = epsilon( epsmch )

    info = 0
    iflag = 0
    nfev = 0
    !
    !  Check the input parameters for errors.
    !
    if ( n <= 0 ) then
        return
    else if ( xtol < 0.0E0_wp ) then
        return
    else if ( maxfev <= 0 ) then
        return
    else if ( ml < 0 ) then
        return
    else if ( mu < 0 ) then
        return
    else if ( factor <= 0.0E0_wp ) then
        return
    else if ( ldfjac < n ) then
        return
    else if ( lr < ( n * (n + 1) ) / 2 ) then
        return
    end if

    if (mode == 2) then
        do j = 1, n
            if (diag(j) <= 0.0E0_wp) then
                go to 300
            end if
        end do
    end if
    !
    !  Evaluate the function at the starting point
    !  and calculate its norm.
    !
    iflag = 1
    call fcn(n, x, fvec, iflag)
    nfev = 1

    if ( iflag < 0 ) then
        go to 300
    end if

    fnorm = enorm2(n, fvec) 
    !
    !  Determine the number of calls to fcn needed to compute the jacobian matrix.
    !
    msum = min(ml+mu+1, n)
    !
    !  Initialize iteration counter and monitors.
    !
    iter = 1
    ncsuc = 0
    ncfail = 0
    nslow1 = 0
    nslow2 = 0
    !
    !  Beginning of the outer loop.
    !
    30  continue

    jeval = .true.
    !
    !  Calculate the jacobian matrix.
    !
    iflag = 2
    call fdjac1(fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn)

    nfev = nfev + msum
    if ( iflag < 0 ) go to 300
    !
    !  Compute the QR factorization of the jacobian.
    !
    call qrfac(n, n, fjac, ldfjac, .false., iwa, 1, wa1, wa2)
    !
    !  On the first iteration, if MODE is 1, scale according
    !  to the norms of the columns of the initial jacobian.
    !
    if ( iter == 1 ) then
        if ( mode /= 2 ) then
            diag(1:n) = wa2(1:n)
            do j = 1, n
                if ( abs(wa2(j)) < epsmch ) then
                    diag(j) = 1.0E0_wp
                end if
            end do
        end if
        !
        !  On the first iteration, calculate the norm of the scaled X
        !  and initialize the step bound DELTA.
        !
        wa3(1:n) = diag(1:n) * x(1:n)
        xnorm = enorm2(n, wa3) !enorm ( n, wa3 )
        delta = factor * xnorm
        if ( abs(delta) < epsmch ) then
            delta = factor
        end if
    end if
    !
    !  Form Q' * FVEC and store in QTF.
    !
    qtf(1:n) = fvec(1:n)

    do j = 1, n
        if ( fjac(j,j) /= 0.0E0_wp ) then
            temp = - dot_product ( qtf(j:n), fjac(j:n,j) ) / fjac(j,j)
            qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp
        end if
    end do
    !
    !  Copy the triangular factor of the QR factorization into R.
    !
    sing = .false.

    do j = 1, n
        l = j
        do i = 1, j-1
            r(l) = fjac(i,j)
            l = l + n - i
        end do
        r(l) = wa1(j)
        if (abs( wa1(j) ) < epsmch) sing = .true.
    end do
    !
    !  Accumulate the orthogonal factor in FJAC.
    !
    call qform ( n, n, fjac, ldfjac )
    !
    !  Rescale if necessary.
    !
    if ( mode /= 2 ) then
        do j = 1, n
            diag(j) = max ( diag(j), wa2(j) )
        end do
    end if
    !
    !  Beginning of the inner loop.
    !
    180 continue
    !
    !  If requested, call fcn to enable printing of iterates.
    !
    if ( 0 < nprint ) then
        iflag = 0
        if (mod( iter-1, nprint ) == 0) then
            call fcn( n, x, fvec, iflag )
        end if
        if (iflag < 0) go to 300
    end if
    !
    !  Determine the direction P.
    !
    call dogleg( n, r, lr, diag, qtf, delta, wa1 )
    !
    !  Store the direction P and X + P.
    !  Calculate the norm of P.
    !
    wa1(1:n) = - wa1(1:n)
    wa2(1:n) = x(1:n) + wa1(1:n)
    wa3(1:n) = diag(1:n) * wa1(1:n)

    pnorm = enorm2(n, wa3) !enorm ( n, wa3 )
    !
    !  On the first iteration, adjust the initial step bound.
    !
    if ( iter == 1 ) then
        delta = min( delta, pnorm )
    end if
    !
    !  Evaluate the function at X + P and calculate its norm.
    !
    iflag = 1
    call fcn(n, wa2, wa4, iflag)
    nfev = nfev + 1

    if ( iflag < 0 ) then
        go to 300
    end if

    fnorm1 = enorm2 ( n, wa4 ) !enorm ( n, wa4 )
    !
    !  Compute the scaled actual reduction.
    !
    actred = -1.0E0_wp
    if (fnorm1 < fnorm) then
        actred = 1.0E0_wp - (fnorm1 / fnorm)**2
    endif
    !
    !  Compute the scaled predicted reduction.
    !
    l = 1
    do i = 1, n
        sum2 = 0.0E0_wp
        do j = i, n
            sum2 = sum2 + r(l) * wa1(j)
            l = l + 1
        end do
        wa3(i) = qtf(i) + sum2
    end do

    temp = enorm2 ( n, wa3 ) !enorm ( n, wa3 )
    prered = 0.0E0_wp
    if ( temp < fnorm ) then
        prered = 1.0E0_wp - ( temp / fnorm )**2
    end if
    !
    !  Compute the ratio of the actual to the predicted reduction.
    !
    ratio = 0.0E0_wp
    if ( 0.0E0_wp < prered ) then
        ratio = actred / prered
    end if
    !
    !  Update the step bound.
    !
    if ( ratio < 0.1E0_wp ) then
        ncsuc = 0
        ncfail = ncfail + 1
        delta = 0.5E0_wp * delta
    else
        ncfail = 0
        ncsuc = ncsuc + 1
        if ( 0.5E0_wp <= ratio .or. 1 < ncsuc ) then
            delta = max ( delta, pnorm / 0.5E0_wp )
        end if
        if ( abs( ratio - 1.0E0_wp ) <= 0.1E0_wp ) then
            delta = pnorm / 0.5E0_wp
        end if
    end if
    !
    !  Test for successful iteration.
    !
    !  Successful iteration.
    !  Update x, fvec, and their norms.
    !
    if ( 0.0001E0_wp <= ratio) then
        x(1:n) = wa2(1:n)
        wa2(1:n) = diag(1:n) * x(1:n)
        fvec(1:n) = wa4(1:n)
        xnorm = enorm2 ( n, wa2 ) !enorm ( n, wa2 )
        fnorm = fnorm1
        iter = iter + 1
    end if
    !
    !  Determine the progress of the iteration.
    !
    nslow1 = nslow1 + 1
    if ( 0.001E0_wp <= actred ) then
        nslow1 = 0
    end if

    if ( jeval ) then
        nslow2 = nslow2 + 1
    end if

    if ( 0.1E0_wp <= actred ) then
        nslow2 = 0
    end if
    !
    !  Test for convergence.
    !
    if ( delta <= xtol * xnorm .or. abs(fnorm) <= epsmch ) then
        info = 1
    end if

    if ( info /= 0 ) then
        go to 300
    end if
    !
    !  Tests for termination and stringent tolerances.
    !
    if ( maxfev <= nfev ) then
        info = 2
    end if

    if ( 0.1E0_wp * max ( 0.1E0_wp * delta, pnorm ) <= epsmch * xnorm ) then
        info = 3
    end if

    if ( nslow2 == 5 ) then
        info = 4
    end if

    if ( nslow1 == 10 ) info = 5
    if ( info /= 0 ) go to 300
    !
    !  Criterion for recalculating jacobian approximation
    !  by forward differences.
    !
    if ( ncfail == 2 ) go to 290
    !
    !  Calculate the rank one modification to the jacobian
    !  and update QTF if necessary.
    !
    do j = 1, n
        sum2 = dot_product ( wa4(1:n), fjac(1:n,j) )
        wa2(j) = ( sum2 - wa3(j) ) / pnorm
        wa1(j) = diag(j) * ( ( diag(j) * wa1(j) ) / pnorm )
        if ( 0.0001E0_wp <= ratio ) then
            qtf(j) = sum2
        end if
    end do
    !
    !  Compute the QR factorization of the updated jacobian.
    !
    call r1updt ( n, n, r, lr, wa1, wa2, wa3, sing )
    call r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 )
    call r1mpyq ( 1, n, qtf, 1, wa2, wa3 )
    !
    !  End of the inner loop.
    !
    jeval = .false.
    go to 180
    290 continue
    !
    !  End of the outer loop.
    !
    go to 30
    300 continue
    !
    !  Termination, either normal or user imposed.
    !
    if ( iflag < 0 ) then
        info = iflag
    end if

    iflag = 0

    if ( 0 < nprint ) then
        call fcn( n, x, fvec, iflag )
    end if

    end subroutine hybrd
    !-----------------------------------------------------------

        
    subroutine hybrd1(fcn, n, x, fvec, tol, info)
    !*****************************************************************************
    !
    !! hybrd1 seeks a zero of N nonlinear equations in N variables.
    !
    !  Discussion:
    !
    !    hybrd1 finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  This is done by using the
    !    more general nonlinear equation solver hybrd.  The user must provide a
    !    subroutine which calculates the functions.  The jacobian is then
    !    calculated by a forward-difference approximation.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, external fcn, the name of the user-supplied subroutine which
    !    calculates the functions.  The routine should have the form:
    !
    !      subroutine fcn(n, x, fvec, iflag)
    !      integer,intent(in)       :: n
    !      real(wp),intent(inout)   :: x(n)
    !      real(wp),intent(out)     :: fvec(n)
    !      integer,intent(out)      :: iflag
    !
    !    The value of IFLAG should not be changed by fcn unless
    !    the user wants to terminate execution of the routine.
    !    In this case set IFLAG to a negative integer.
    !
    !    Input, integer :: N, the number of functions and variables.
    !
    !    Input/output, real(wp) :: X(N).  On input, X must contain an initial
    !    estimate of the solution vector.  On output X contains the final
    !    estimate of the solution vector.
    !
    !    Output, real(wp) :: FVEC(N), the functions evaluated at the output X.
    !
    !    Input, real(wp) :: TOL.  Termination occurs when the algorithm
    !    estimates that the relative error between X and the solution is at
    !    most TOL.  TOL should be nonnegative.
    !
    !    Output, integer :: INFO, error flag.  If the user has terminated execution,
    !    INFO is set to the (negative) value of IFLAG. See the description of fcn.
    !    Otherwise, INFO is set as follows:
    !    0, improper input parameters.
    !    1, algorithm estimates that the relative error between X and the
    !       solution is at most TOL.
    !    2, number of calls to fcn has reached or exceeded 500*(N+1).
    !    3, TOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    4, the iteration is not making good progress.
    !
    implicit none
    
    !interface arguments:    
    external fcn
    integer,intent(in) :: n
    real(wp),intent(inout) :: x(n)
    real(wp),intent(out) :: fvec(n)
    real(wp),intent(in) :: tol
    integer,intent(out) :: info
    
    !local variables:
    real(wp) :: diag(n)
    real(wp) :: epsfcn
    real(wp) :: factor
    real(wp) :: fjac(n,n)
    integer :: ldfjac
    integer :: lr
    integer :: maxfev
    integer :: ml
    integer :: mode
    integer :: mu
    integer :: nfev
    integer :: nprint
    real(wp) :: qtf(n)
    real(wp) :: r((n*(n+1))/2)
    real(wp) :: xtol
    !..............................
    
    info = 0

    if (n <= 0) then
        return
    else if (tol < 0.0E0_wp) then
        return
    end if

    maxfev = 500*(n + 1)
    xtol = tol
    ml = n - 1
    mu = n - 1
    epsfcn = 0.0E0_wp
    mode = 2
    diag(1:n) = 1.0E0_wp
    nprint = 0
    lr = (n*(n + 1)) / 2
    factor = 100.0E0_wp
    ldfjac = n

    call hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
        factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf)

    if ( info == 5 ) then
        info = 4
    end if
    
    end subroutine hybrd1
    !-----------------------------------------------------------

    
    subroutine hybrj(fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, &
        factor, nprint, info, nfev, njev, r, lr, qtf)
    !*****************************************************************************80
    !
    !! hybrj seeks a zero of N nonlinear equations in N variables.
    !
    !  Discussion:
    !
    !    hybrj finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  The user must provide a
    !    subroutine which calculates the functions and the jacobian.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, external fcn, the name of the user-supplied subroutine which
    !    calculates the functions and the jacobian.  fcn should have the form:
    !
    !      subroutine fcn( n, x, fvec, fjac, ldfjac, iflag )
    !      integer :: ldfjac
    !      integer :: n
    !      real fjac(ldfjac,n)
    !      real(wp) fvec(n)
    !      integer :: iflag
    !      real(wp) x(n)
    !
    !    If IFLAG = 1 on intput, fcn should calculate the functions at X and
    !    return this vector in FVEC.
    !    If IFLAG = 2 on input, fcn should calculate the jacobian at X and
    !    return this matrix in FJAC.
    !    To terminate the algorithm, the user may set IFLAG negative.
    !
    !    Input, integer :: N, the number of functions and variables.
    !
    !    Input/output, real(wp) :: X(N).  On input, X must contain an initial
    !    estimate of the solution vector.  On output X contains the final
    !    estimate of the solution vector.
    !
    !    Output, real(wp) :: FVEC(N), the functions evaluated at the output X.
    !
    !    Output, real(wp) :: FJAC(LDFJAC,N), an N by N matrix, containing
    !    the orthogonal matrix Q produced by the QR factorization
    !    of the final approximate jacobian.
    !
    !    Input, integer :: LDFJAC, the leading dimension of the array FJAC.
    !    LDFJAC must be at least N.
    !
    !    Input, real(wp) :: XTOL.  Termination occurs when the relative error
    !    between two consecutive iterates is at most XTOL.  XTOL should be
    !    nonnegative.
    !
    !    Input, integer :: MAXFEV.  Termination occurs when the number of calls
    !    to fcn is at least MAXFEV by the end of an iteration.
    !
    !    Input/output, real(wp) :: DIAG(N).  If MODE = 1, then DIAG is set
    !    internally.  If MODE = 2, then DIAG must contain positive entries that
    !    serve as multiplicative scale factors for the variables.
    !
    !    Input, integer :: MODE, scaling option.
    !    1, variables will be scaled internally.
    !    2, scaling is specified by the input DIAG vector.
    !
    !    Input, real(wp) :: FACTOR, determines the initial step bound.  This
    !    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
    !    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
    !    in the interval (0.1, 100) with 100 the recommended value.
    !
    !    Input, integer :: NPRINT, enables controlled printing of iterates if it
    !    is positive.  In this case, fcn is called with IFLAG = 0 at the
    !    beginning of the first iteration and every NPRINT iterations thereafter
    !    and immediately prior to return, with X and FVEC available
    !    for printing.  If NPRINT is not positive, no special calls
    !    of fcn with IFLAG = 0 are made.
    !
    !    Output, integer :: INFO, error flag.  If the user has terminated execution,
    !    INFO is set to the (negative) value of IFLAG.  See the description of fcn.
    !    Otherwise, INFO is set as follows:
    !    0, improper input parameters.
    !    1, relative error between two consecutive iterates is at most XTOL.
    !    2, number of calls to fcn with IFLAG = 1 has reached MAXFEV.
    !    3, XTOL is too small.  No further improvement in
    !       the approximate solution X is possible.
    !    4, iteration is not making good progress, as measured by the
    !       improvement from the last five jacobian evaluations.
    !    5, iteration is not making good progress, as measured by the
    !       improvement from the last ten iterations.
    !
    !    Output, integer :: NFEV, the number of calls to fcn with IFLAG = 1.
    !
    !    Output, integer :: NJEV, the number of calls to fcn with IFLAG = 2.
    !
    !    Output, real(wp) :: R(LR), the upper triangular matrix produced
    !    by the QR factorization of the final approximate jacobian, stored rowwise.
    !
    !    Input, integer :: LR, the size of the R array, which must be no less
    !    than (N*(N+1))/2.
    !
    !    Output, real(wp) :: QTF(N), contains the vector Q'*FVEC.
    !
    implicit none

    integer :: ldfjac
    integer :: lr
    integer :: n

    real(wp) :: actred
    real(wp) :: delta
    real(wp) :: diag(n)
    real(wp) :: epsmch
    real(wp) :: factor
    external fcn
    real(wp) :: fjac(ldfjac,n)
    real(wp) :: fnorm
    real(wp) :: fnorm1
    real(wp) :: fvec(n)
    integer :: i
    integer :: iflag
    integer :: info
    integer :: iter
    integer :: iwa(1)
    integer :: j
    logical :: jeval
    integer :: l
    integer :: maxfev
    integer :: mode
    integer :: ncfail
    integer :: nslow1
    integer :: nslow2
    integer :: ncsuc
    integer :: nfev
    integer :: njev
    integer :: nprint
    real(wp) :: pnorm
    real(wp) :: prered
    real(wp) :: qtf(n)
    real(wp) :: r(lr)
    real(wp) :: ratio
    logical :: sing
    real(wp) :: sum2
    real(wp) :: temp
    real(wp) :: wa1(n)
    real(wp) :: wa2(n)
    real(wp) :: wa3(n)
    real(wp) :: wa4(n)
    real(wp) :: x(n)
    real(wp) :: xnorm
    real(wp) :: xtol

    epsmch = epsilon ( epsmch )

    info = 0
    iflag = 0
    nfev = 0
    njev = 0
    !
    !  Check the input parameters for errors.
    !
    if ( n <= 0 ) then
        go to 300
    end if

    if ( ldfjac < n .or. &
        xtol < 0.0E0_wp .or. &
        maxfev <= 0 .or. &
        factor <= 0.0E0_wp .or. &
        lr < (n*(n + 1))/2 ) then
    go to 300
    end if

    if ( mode == 2 ) then
        do j = 1, n
            if ( diag(j) <= 0.0E0_wp ) go to 300
        end do
    end if

    20  continue
    !
    !  Evaluate the function at the starting point
    !  and calculate its norm.
    !
    iflag = 1
    call fcn( n, x, fvec, fjac, ldfjac, iflag )
    nfev = 1
    if ( iflag < 0 ) go to 300
    fnorm = enorm2 ( n, fvec ) !enorm ( n, fvec )
    !
    !  Initialize iteration counter and monitors.
    !
    iter = 1
    ncsuc = 0
    ncfail = 0
    nslow1 = 0
    nslow2 = 0
    !
    !  Beginning of the outer loop.
    !
    30  continue

    jeval = .true.
    !
    !  Calculate the jacobian matrix.
    !
    iflag = 2
    call fcn( n, x, fvec, fjac, ldfjac, iflag )
    njev = njev + 1
    if ( iflag < 0 ) go to 300
    !
    !  Compute the QR factorization of the jacobian.
    !
    call qrfac ( n, n, fjac, ldfjac, .false., iwa, 1, wa1, wa2 )
    !
    !  On the first iteration, if MODE is 1, scale according
    !  to the norms of the columns of the initial jacobian.
    !
    if ( iter == 1 ) then

        if ( mode /= 2 ) then
            diag(1:n) = wa2(1:n)
            do j = 1, n
                if ( abs(wa2(j)) < epsmch ) diag(j) = 1.0E0_wp
            end do
        end if
        !
        !  On the first iteration, calculate the norm of the scaled X
        !  and initialize the step bound DELTA.
        !
        wa3(1:n) = diag(1:n) * x(1:n)
        xnorm = enorm2 ( n, wa3 ) !enorm ( n, wa3 )
        delta = factor * xnorm
        if ( abs(delta) < epsmch ) then
            delta = factor
        end if
    end if
    !
    !  Form Q'*FVEC and store in QTF.
    !
    qtf(1:n) = fvec(1:n)

    do j = 1, n
        if ( fjac(j,j) /= 0.0E0_wp ) then
            sum2 = 0.0E0_wp
            do i = j, n
                sum2 = sum2 + fjac(i,j) * qtf(i)
            end do
            temp = - sum2 / fjac(j,j)
            do i = j, n
                qtf(i) = qtf(i) + fjac(i,j) * temp
            end do
        end if
    end do
    !
    !  Copy the triangular factor of the QR factorization into R.
    !
    sing = .false.
    do j = 1, n
        l = j
        do i = 1, j-1
            r(l) = fjac(i,j)
            l = l + n - i
        end do
        r(l) = wa1(j)
        if ( abs(wa1(j)) < epsmch) then
            sing = .true.
        end if
    end do
    !
    !  Accumulate the orthogonal factor in FJAC.
    !
    call qform ( n, n, fjac, ldfjac )
    !
    !  Rescale if necessary.
    !
    if ( mode /= 2 ) then
        do j = 1, n
            diag(j) = max ( diag(j), wa2(j) )
        end do
    end if
    !
    !  Beginning of the inner loop.
    !
    180 continue
    !
    !  If requested, call fcn to enable printing of iterates.
    !
    if ( 0 < nprint ) then
        iflag = 0
        if ( mod ( iter-1, nprint ) == 0 ) then
            call fcn( n, x, fvec, fjac, ldfjac, iflag )
        end if
        if ( iflag < 0 ) then
            go to 300
        end if
    end if
    !
    !  Determine the direction P.
    !
    call dogleg ( n, r, lr, diag, qtf, delta, wa1 )
    !
    !  Store the direction P and X + P.
    !  Calculate the norm of P.
    !
    wa1(1:n) = - wa1(1:n)
    wa2(1:n) = x(1:n) + wa1(1:n)
    wa3(1:n) = diag(1:n) * wa1(1:n)
    pnorm = enorm2 ( n, wa3 ) !enorm ( n, wa3 )
    !
    !  On the first iteration, adjust the initial step bound.
    !
    if ( iter == 1 ) then
        delta = min ( delta, pnorm )
    end if
    !
    !  Evaluate the function at X + P and calculate its norm.
    !
    iflag = 1
    call fcn ( n, wa2, wa4, fjac, ldfjac, iflag )
    nfev = nfev + 1
    if ( iflag < 0 ) go to 300
    fnorm1 = enorm2 ( n, wa4 ) !enorm ( n, wa4 )
    !
    !  Compute the scaled actual reduction.
    !
    actred = -1.0E0_wp
    if ( fnorm1 < fnorm ) then
        actred = 1.0E0_wp - ( fnorm1 / fnorm )**2
    end if
    !
    !  Compute the scaled predicted reduction.
    !
    l = 1
    do i = 1, n
        sum2 = 0.0E0_wp
        do j = i, n
            sum2 = sum2 + r(l) * wa1(j)
            l = l + 1
        end do
        wa3(i) = qtf(i) + sum2
    end do

    temp = enorm2 ( n, wa3 ) !enorm ( n, wa3 )
    prered = 0.0E0_wp
    if ( temp < fnorm ) then
        prered = 1.0E0_wp - ( temp / fnorm )**2
    end if
    !
    !  Compute the ratio of the actual to the predicted reduction.
    !
    if ( 0.0E0_wp < prered ) then
        ratio = actred / prered
    else
        ratio = 0.0E0_wp
    end if
    !
    !  Update the step bound.
    !
    if ( ratio < 0.1E0_wp ) then
        ncsuc = 0
        ncfail = ncfail + 1
        delta = 0.5E0_wp * delta
    else
        ncfail = 0
        ncsuc = ncsuc + 1

        if ( 0.5E0_wp <= ratio .or. 1 < ncsuc ) then
            delta = max ( delta, pnorm / 0.5E0_wp )
        end if

        if ( abs( ratio - 1.0E0_wp ) <= 0.1E0_wp ) then
            delta = pnorm / 0.5E0_wp
        end if

    end if
    !
    !  Test for successful iteration.
    !

    !
    !  Successful iteration.
    !  Update X, FVEC, and their norms.
    !
    if ( 0.0001E0_wp <= ratio ) then
        x(1:n) = wa2(1:n)
        wa2(1:n) = diag(1:n) * x(1:n)
        fvec(1:n) = wa4(1:n)
        xnorm = enorm2( n, wa2 ) !enorm ( n, wa2 )
        fnorm = fnorm1
        iter = iter + 1
    end if
    !
    !  Determine the progress of the iteration.
    !
    nslow1 = nslow1 + 1
    if ( 0.001E0_wp <= actred ) then
        nslow1 = 0
    end if

    if ( jeval ) then
        nslow2 = nslow2 + 1
    end if

    if ( 0.1E0_wp <= actred ) then
        nslow2 = 0
    end if
    !
    !  Test for convergence.
    !
    if ( delta <= xtol * xnorm .or. abs(fnorm) < epsmch) then
        info = 1
    end if

    if ( info /= 0 ) go to 300
    !
    !  Tests for termination and stringent tolerances.
    !
    if ( maxfev <= nfev ) then
        info = 2
    end if

    if ( 0.1E0_wp *max ( 0.1E0_wp * delta, pnorm ) <= epsmch * xnorm ) then
        info = 3
    end if

    if ( nslow2 == 5 ) then
        info = 4
    end if

    if ( nslow1 == 10 ) then
        info = 5
    end if

    if ( info /= 0 ) then
        go to 300
    end if
    !
    !  Criterion for recalculating jacobian.
    !
    if ( ncfail == 2 ) then
        go to 290
    end if
    !
    !  Calculate the rank one modification to the jacobian
    !  and update QTF if necessary.
    !
    do j = 1, n
        sum2 = dot_product ( wa4(1:n), fjac(1:n,j) )
        wa2(j) = ( sum2 - wa3(j) ) / pnorm
        wa1(j) = diag(j) * ( ( diag(j) * wa1(j) ) / pnorm )
        if ( ratio >= 0.0001E0_wp ) then
            qtf(j) = sum2
        end if
    end do
    !
    !  Compute the QR factorization of the updated jacobian.
    !
    call r1updt ( n, n, r, lr, wa1, wa2, wa3, sing )
    call r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 )
    call r1mpyq ( 1, n, qtf, 1, wa2, wa3 )
    !
    !  End of the inner loop.
    !
    jeval = .false.
    go to 180

    290 continue
    !
    !  End of the outer loop.
    !
    go to 30

    300 continue
    !
    !  Termination, either normal or user imposed.
    !
    if ( iflag < 0 ) then
        info = iflag
    end if

    iflag = 0

    if ( nprint > 0 ) then
        call fcn( n, x, fvec, fjac, ldfjac, iflag )
    end if

    return
    end
    !-----------------------------------------------------------

        
    subroutine hybrj1( fcn, n, x, fvec, fjac, ldfjac, tol, info )
    !*****************************************************************************80
    !
    !! hybrj1 seeks a zero of N nonlinear equations in N variables by Powell's method.
    !
    !  Discussion:
    !
    !    hybrj1 finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  This is done by using the
    !    more general nonlinear equation solver hybrj.  The user
    !    must provide a subroutine which calculates the functions
    !    and the jacobian.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, external fcn, the name of the user-supplied subroutine which
    !    calculates the functions and the jacobian.  fcn should have the form:
    !
    !      subroutine fcn( n, x, fvec, fjac, ldfjac, iflag )
    !      integer :: ldfjac
    !      integer :: n
    !      real fjac(ldfjac,n)
    !      real(wp) fvec(n)
    !      integer :: iflag
    !      real(wp) x(n)
    !
    !    If IFLAG = 1 on intput, fcn should calculate the functions at X and
    !    return this vector in FVEC.
    !    If IFLAG = 2 on input, fcn should calculate the jacobian at X and
    !    return this matrix in FJAC.
    !    To terminate the algorithm, the user may set IFLAG negative.
    !
    !    Input, integer :: N, the number of functions and variables.
    !
    !    Input/output, real(wp) :: X(N).  On input, X must contain an initial
    !    estimate of the solution vector.  On output X contains the final
    !    estimate of the solution vector.
    !
    !    Output, real(wp) :: FVEC(N), the functions evaluated at the output X.
    !
    !    Output, real(wp) :: FJAC(LDFJAC,N), an N by N array which contains
    !    the orthogonal matrix Q produced by the QR factorization of the final
    !    approximate jacobian.
    !
    !    Input, integer :: LDFJAC, the leading dimension of the array FJAC.
    !    LDFJAC must be at least N.
    !
    !    Input, real(wp) :: TOL.  Termination occurs when the algorithm
    !    estimates that the relative error between X and the solution is at most
    !    TOL.  TOL should be nonnegative.
    !
    !    Output, integer :: INFO, error flag.  If the user has terminated execution,
    !    INFO is set to the (negative) value of IFLAG. See the description of fcn.
    !    Otherwise, INFO is set as follows:
    !    0, improper input parameters.
    !    1, algorithm estimates that the relative error between X and the
    !       solution is at most TOL.
    !    2, number of calls to fcn with IFLAG = 1 has reached 100*(N+1).
    !    3, TOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    4, iteration is not making good progress.
    !
    implicit none

    integer :: ldfjac
    integer :: n

    real(wp) :: diag(n)
    real(wp) :: factor
    external fcn
    real(wp) :: fjac(ldfjac,n)
    real(wp) :: fvec(n)
    integer :: info
    integer :: lr
    integer :: maxfev
    integer :: mode
    integer :: nfev
    integer :: njev
    integer :: nprint
    real(wp) :: qtf(n)
    real(wp) :: r((n*(n+1))/2)
    real(wp) :: tol
    real(wp) :: x(n)
    real(wp) :: xtol

    info = 0

    if ( n <= 0 ) then
        return
    else if ( ldfjac < n ) then
        return
    else if ( tol < 0.0E0_wp ) then
        return
    end if

    maxfev = 100 * ( n + 1 )
    xtol = tol
    mode = 2
    diag(1:n) = 1.0E0_wp
    factor = 100.0E0_wp
    nprint = 0
    lr = ( n * ( n + 1 ) ) / 2

    call hybrj( fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, &
        factor, nprint, info, nfev, njev, r, lr, qtf )

    if ( info == 5 ) then
        info = 4
    end if

    end subroutine hybrj1
    !-----------------------------------------------------------


    pure subroutine i_swap( i, j )
    !*****************************************************************************80
    !
    !! i_swap swaps two integer :: values.
    !
    !  Modified:
    !
    !    30 November 1998
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, integer :: I, J.  On output, the values of I and
    !    J have been interchanged.
    !
    implicit none

    integer,intent(inout) :: i
    integer,intent(inout) :: j
    integer :: k

    k = i
    i = j
    j = k

    end subroutine i_swap
    !-----------------------------------------------------------

    
    subroutine lmder ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
        diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )
    !*****************************************************************************80
    !
    !! lmder minimizes M functions in N variables using the Levenberg-Marquardt method.
    !
    !  Discussion:
    !
    !    lmder minimizes the sum of the squares of M nonlinear functions in
    !    N variables by a modification of the Levenberg-Marquardt algorithm.
    !    The user must provide a subroutine which calculates the functions
    !    and the jacobian.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, external fcn, the name of the user-supplied subroutine which
    !    calculates the functions and the jacobian.  fcn should have the form:
    !
    !      subroutine fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
    !      integer :: ldfjac
    !      integer :: n
    !      real fjac(ldfjac,n)
    !      real(wp) fvec(m)
    !      integer :: iflag
    !      real(wp) x(n)
    !
    !    If IFLAG = 1 on intput, fcn should calculate the functions at X and
    !    return this vector in FVEC.
    !    If IFLAG = 2 on input, fcn should calculate the jacobian at X and
    !    return this matrix in FJAC.
    !    To terminate the algorithm, the user may set IFLAG negative.
    !
    !    Input, integer :: M, is the number of functions.
    !
    !    Input, integer :: N, is the number of variables.  N must not exceed M.
    !
    !    Input/output, real(wp) :: X(N).  On input, X must contain an initial
    !    estimate of the solution vector.  On output X contains the final
    !    estimate of the solution vector.
    !
    !    Output, real(wp) :: FVEC(M), the functions evaluated at the output X.
    !
    !    Output, real(wp) :: FJAC(LDFJAC,N), an M by N array.  The upper
    !    N by N submatrix of FJAC contains an upper triangular matrix R with
    !    diagonal elements of nonincreasing magnitude such that
    !      P' * ( JAC' * JAC ) * P = R' * R,
    !    where P is a permutation matrix and JAC is the final calculated jacobian.
    !    Column J of P is column IPVT(J) of the identity matrix.  The lower
    !    trapezoidal part of FJAC contains information generated during
    !    the computation of R.
    !
    !    Input, integer :: LDFJAC, the leading dimension of the array FJAC.
    !    LDFJAC must be at least M.
    !
    !    Input, real(wp) :: FTOL.  Termination occurs when both the actual
    !    and predicted relative reductions in the sum of squares are at most FTOL.
    !    Therefore, FTOL measures the relative error desired in the sum of
    !    squares.  FTOL should be nonnegative.
    !
    !    Input, real(wp) :: XTOL.  Termination occurs when the relative error
    !    between two consecutive iterates is at most XTOL.  XTOL should be
    !    nonnegative.
    !
    !    Input, real(wp) :: GTOL.  Termination occurs when the cosine of the
    !    angle between FVEC and any column of the jacobian is at most GTOL in
    !    absolute value.  Therefore, GTOL measures the orthogonality desired
    !    between the function vector and the columns of the jacobian.  GTOL should
    !    be nonnegative.
    !
    !    Input, integer :: MAXFEV.  Termination occurs when the number of calls
    !    to fcn with IFLAG = 1 is at least MAXFEV by the end of an iteration.
    !
    !    Input/output, real(wp) :: DIAG(N).  If MODE = 1, then DIAG is set
    !    internally.  If MODE = 2, then DIAG must contain positive entries that
    !    serve as multiplicative scale factors for the variables.
    !
    !    Input, integer :: MODE, scaling option.
    !    1, variables will be scaled internally.
    !    2, scaling is specified by the input DIAG vector.
    !
    !    Input, real(wp) :: FACTOR, determines the initial step bound.  This
    !    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
    !    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
    !    in the interval (0.1, 100) with 100 the recommended value.
    !
    !    Input, integer :: NPRINT, enables controlled printing of iterates if it
    !    is positive.  In this case, fcn is called with IFLAG = 0 at the
    !    beginning of the first iteration and every NPRINT iterations thereafter
    !    and immediately prior to return, with X and FVEC available
    !    for printing.  If NPRINT is not positive, no special calls
    !    of fcn with IFLAG = 0 are made.
    !
    !    Output, integer :: INFO, error flag.  If the user has terminated execution,
    !    INFO is set to the (negative) value of IFLAG. See the description of fcn.
    !    Otherwise, INFO is set as follows:
    !    0, improper input parameters.
    !    1, both actual and predicted relative reductions in the sum of
    !       squares are at most FTOL.
    !    2, relative error between two consecutive iterates is at most XTOL.
    !    3, conditions for INFO = 1 and INFO = 2 both hold.
    !    4, the cosine of the angle between FVEC and any column of the jacobian
    !       is at most GTOL in absolute value.
    !    5, number of calls to fcn with IFLAG = 1 has reached MAXFEV.
    !    6, FTOL is too small.  No further reduction in the sum of squares
    !       is possible.
    !    7, XTOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    8, GTOL is too small.  FVEC is orthogonal to the columns of the
    !       jacobian to machine precision.
    !
    !    Output, integer :: NFEV, the number of calls to fcn with IFLAG = 1.
    !
    !    Output, integer :: NJEV, the number of calls to fcn with IFLAG = 2.
    !
    !    Output, integer :: IPVT(N), defines a permutation matrix P such that
    !    JAC*P = Q*R, where JAC is the final calculated jacobian, Q is
    !    orthogonal (not stored), and R is upper triangular with diagonal
    !    elements of nonincreasing magnitude.  Column J of P is column
    !    IPVT(J) of the identity matrix.
    !
    !    Output, real(wp) :: QTF(N), contains the first N elements of Q'*FVEC.
    !
    implicit none

    integer :: ldfjac
    integer :: m
    integer :: n

    real(wp) :: actred
    real(wp) :: delta
    real(wp) :: diag(n)
    real(wp) :: dirder
    real(wp) :: epsmch
    real(wp) :: factor
    external fcn
    real(wp) :: fjac(ldfjac,n)
    real(wp) :: fnorm
    real(wp) :: fnorm1
    real(wp) :: ftol
    real(wp) :: fvec(m)
    real(wp) :: gnorm
    real(wp) :: gtol
    !integer :: i
    integer :: iflag
    integer :: info
    integer :: ipvt(n)
    integer :: iter
    integer :: j
    integer :: l
    integer :: maxfev
    integer :: mode
    integer :: nfev
    integer :: njev
    integer :: nprint
    real(wp) :: par
    real(wp) :: pnorm
    real(wp) :: prered
    real(wp) :: qtf(n)
    real(wp) :: ratio
    real(wp) :: sum2
    real(wp) :: temp
    real(wp) :: temp1
    real(wp) :: temp2
    real(wp) :: wa1(n)
    real(wp) :: wa2(n)
    real(wp) :: wa3(n)
    real(wp) :: wa4(m)
    real(wp) :: xnorm
    real(wp) :: x(n)
    real(wp) :: xtol

    epsmch = epsilon ( epsmch )

    info = 0
    iflag = 0
    nfev = 0
    njev = 0
    !
    !  Check the input parameters for errors.
    !
    if ( n <= 0 ) then
        go to 300
    end if

    if ( m < n ) then
        go to 300
    end if

    if ( ldfjac < m &
        .or. ftol < 0.0E0_wp .or. xtol < 0.0E0_wp .or. gtol < 0.0E0_wp &
        .or. maxfev <= 0 .or. factor <= 0.0E0_wp ) then
    go to 300
    end if

    if ( mode == 2 ) then
        do j = 1, n
            if ( diag(j) <= 0.0E0_wp ) go to 300
        end do
    end if
    !
    !  Evaluate the function at the starting point and calculate its norm.
    !
    iflag = 1
    call fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
    nfev = 1
    if ( iflag < 0 ) then
        go to 300
    end if

    fnorm = enorm2 ( m, fvec ) !enorm ( m, fvec )
    !
    !  Initialize Levenberg-Marquardt parameter and iteration counter.
    !
    par = 0.0E0_wp
    iter = 1
    !
    !  Beginning of the outer loop.
    !
    30  continue
    !
    !  Calculate the jacobian matrix.
    !
    iflag = 2
    call fcn ( m, n, x, fvec, fjac, ldfjac, iflag )

    njev = njev + 1

    if ( iflag < 0 ) then
        go to 300
    end if
    !
    !  If requested, call fcn to enable printing of iterates.
    !
    if ( 0 < nprint ) then
        iflag = 0
        if ( mod ( iter-1, nprint ) == 0 ) then
            call fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
        end if
        if ( iflag < 0 ) then
            go to 300
        end if
    end if
    !
    !  Compute the QR factorization of the jacobian.
    !
    call qrfac ( m, n, fjac, ldfjac, .true., ipvt, n, wa1, wa2 )
    !
    !  On the first iteration and if mode is 1, scale according
    !  to the norms of the columns of the initial jacobian.
    !
    if ( iter == 1 ) then
        if ( mode /= 2 ) then
            diag(1:n) = wa2(1:n)
            do j = 1, n
                if ( abs(wa2(j)) < epsmch ) then
                    diag(j) = 1.0E0_wp
                end if
            end do
        end if
        !
        !  On the first iteration, calculate the norm of the scaled X
        !  and initialize the step bound DELTA.
        !
        wa3(1:n) = diag(1:n) * x(1:n)

        xnorm = enorm2 ( n, wa3 ) !enorm ( n, wa3 )
        delta = factor * xnorm
        if ( abs(delta) < epsmch ) then
            delta = factor
        end if
    end if
    !
    !  Form Q'*FVEC and store the first N components in QTF.
    !
    wa4(1:m) = fvec(1:m)

    do j = 1, n
        if ( fjac(j,j) /= 0.0E0_wp ) then
            sum2 = dot_product ( wa4(j:m), fjac(j:m,j) )
            temp = - sum2 / fjac(j,j)
            wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
        end if
        fjac(j,j) = wa1(j)
        qtf(j) = wa4(j)
    end do
    !
    !  Compute the norm of the scaled gradient.
    !
    gnorm = 0.0E0_wp

    if ( fnorm /= 0.0E0_wp ) then
        do j = 1, n
            l = ipvt(j)
            if ( wa2(l) /= 0.0E0_wp ) then
                sum2 = dot_product ( qtf(1:j), fjac(1:j,j) ) / fnorm
                gnorm = max ( gnorm, abs( sum2 / wa2(l) ) )
            end if
        end do
    end if
    !
    !  Test for convergence of the gradient norm.
    !
    if ( gnorm <= gtol ) then
        info = 4
        go to 300
    end if
    !
    !  Rescale if necessary.
    !
    if ( mode /= 2 ) then
        do j = 1, n
            diag(j) = max ( diag(j), wa2(j) )
        end do
    end if
    !
    !  Beginning of the inner loop.
    !
    200 continue
    !
    !  Determine the Levenberg-Marquardt parameter.
    !
    call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
    !
    !  Store the direction p and x + p. calculate the norm of p.
    !
    wa1(1:n) = - wa1(1:n)
    wa2(1:n) = x(1:n) + wa1(1:n)
    wa3(1:n) = diag(1:n) * wa1(1:n)

    pnorm = enorm2 ( n, wa3 ) !enorm ( n, wa3 )
    !
    !  On the first iteration, adjust the initial step bound.
    !
    if ( iter == 1 ) then
        delta = min ( delta, pnorm )
    end if
    !
    !  Evaluate the function at x + p and calculate its norm.
    !
    iflag = 1
    call fcn ( m, n, wa2, wa4, fjac, ldfjac, iflag )

    nfev = nfev + 1

    if ( iflag < 0 ) then
        go to 300
    end if

    fnorm1 = enorm2 ( m, wa4 ) !enorm ( m, wa4 )
    !
    !  Compute the scaled actual reduction.
    !
    actred = -1.0E0_wp
    if ( 0.1E0_wp * fnorm1 < fnorm ) then
        actred = 1.0E0_wp - ( fnorm1 / fnorm )**2
    end if
    !
    !  Compute the scaled predicted reduction and
    !  the scaled directional derivative.
    !
    do j = 1, n
        wa3(j) = 0.0E0_wp
        l = ipvt(j)
        temp = wa1(l)
        wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
    end do

    temp1 = enorm2 ( n, wa3 ) / fnorm !enorm ( n, wa3 ) / fnorm
    temp2 = ( sqrt ( par ) * pnorm ) / fnorm
    prered = temp1**2 + temp2**2 / 0.5E0_wp
    dirder = - ( temp1**2 + temp2**2 )
    !
    !  Compute the ratio of the actual to the predicted reduction.
    !
    if ( prered /= 0.0E0_wp ) then
        ratio = actred / prered
    else
        ratio = 0.0E0_wp
    end if
    !
    !  Update the step bound.
    !
    if ( ratio <= 0.25E0_wp ) then
        if ( actred >= 0.0E0_wp ) then
            temp = 0.5E0_wp
        end if
        if ( actred < 0.0E0_wp ) then
            temp = 0.5E0_wp * dirder / ( dirder + 0.5E0_wp * actred )
        end if
        if ( 0.1E0_wp * fnorm1 >= fnorm .or. temp < 0.1E0_wp ) then
            temp = 0.1E0_wp
        end if
        delta = temp * min ( delta, pnorm / 0.1E0_wp )
        par = par / temp

    else
        if (  abs(par) < epsmch .OR. ratio >= 0.75E0_wp ) then
            delta = 2.0E0_wp * pnorm
            par = 0.5E0_wp * par
        end if
    end if
    !
    !  Successful iteration.
    !
    !  Update x, fvec, and their norms.
    !
    if ( ratio >= 0.0001E0_wp ) then
        x(1:n) = wa2(1:n)
        wa2(1:n) = diag(1:n) * x(1:n)
        fvec(1:m) = wa4(1:m)
        xnorm = enorm2( n, wa2 ) !enorm ( n, wa2 )
        fnorm = fnorm1
        iter = iter + 1
    end if
    !
    !  Tests for convergence.
    !
    if ( abs( actred) <= ftol .and. &
        prered <= ftol .and. &
        0.5E0_wp * ratio <= 1.0E0_wp ) then
    info = 1
    end if

    if ( delta <= xtol * xnorm ) then
        info = 2
    end if

    if ( abs( actred) <= ftol .and. prered <= ftol &
        .and. 0.5E0_wp * ratio <= 1.0E0_wp .and. info == 2 ) then
    info = 3
    end if

    if ( info /= 0 ) go to 300
    !
    !  Tests for termination and stringent tolerances.
    !
    if ( nfev >= maxfev ) then
        info = 5
    end if

    if ( abs( actred ) <= epsmch .and. prered <= epsmch &
        .and. 0.5E0_wp * ratio <= 1.0E0_wp ) info = 6
    if ( delta <= epsmch * xnorm ) info = 7
    if ( gnorm <= epsmch ) info = 8
    if ( info /= 0 ) go to 300
    !
    !  End of the inner loop. repeat if iteration unsuccessful.
    !
    if ( ratio < 0.0001E0_wp ) then
        go to 200
    end if
    !
    !  End of the outer loop.
    !
    go to 30
    300 continue
    !
    !  Termination, either normal or user imposed.
    !
    if ( iflag < 0 ) then
        info = iflag
    end if

    iflag = 0

    if ( 0 < nprint ) then
        call fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
    end if

    return
    end
    !-----------------------------------------------------------

        
    subroutine lmder1 ( fcn, m, n, x, fvec, fjac, ldfjac, tol, info )
    !*****************************************************************************80
    !
    !! lmder1 minimizes M functions in N variables using the Levenberg-Marquardt method.
    !
    !  Discussion:
    !
    !    lmder1 minimizes the sum of the squares of M nonlinear functions in
    !    N variables by a modification of the Levenberg-Marquardt algorithm.
    !    This is done by using the more general least-squares solver lmder.
    !    The user must provide a subroutine which calculates the functions
    !    and the jacobian.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, external fcn, the name of the user-supplied subroutine which
    !    calculates the functions and the jacobian.  fcn should have the form:
    !
    !      subroutine fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
    !      integer :: ldfjac
    !      integer :: n
    !      real fjac(ldfjac,n)
    !      real(wp) fvec(m)
    !      integer :: iflag
    !      real(wp) x(n)
    !
    !    If IFLAG = 1 on intput, fcn should calculate the functions at X and
    !    return this vector in FVEC.
    !    If IFLAG = 2 on input, fcn should calculate the jacobian at X and
    !    return this matrix in FJAC.
    !    To terminate the algorithm, the user may set IFLAG negative.
    !
    !    Input, integer :: M, the number of functions.
    !
    !    Input, integer :: N, is the number of variables.  N must not exceed M.
    !
    !    Input/output, real(wp) :: X(N).  On input, X must contain an initial
    !    estimate of the solution vector.  On output X contains the final
    !    estimate of the solution vector.
    !
    !    Output, real(wp) :: FVEC(M), the functions evaluated at the output X.
    !
    !    Output, real(wp) :: FJAC(LDFJAC,N), an M by N array.  The upper
    !    N by N submatrix contains an upper triangular matrix R with
    !    diagonal elements of nonincreasing magnitude such that
    !      P' * ( JAC' * JAC ) * P = R' * R,
    !    where P is a permutation matrix and JAC is the final calculated
    !    jacobian.  Column J of P is column IPVT(J) of the identity matrix.
    !    The lower trapezoidal part of FJAC contains information generated during
    !    the computation of R.
    !
    !    Input, integer :: LDFJAC, is the leading dimension of FJAC,
    !    which must be no less than M.
    !
    !    Input, real(wp) :: TOL.  Termination occurs when the algorithm
    !    estimates either that the relative error in the sum of squares is at
    !    most TOL or that the relative error between X and the solution is at
    !    most TOL.
    !
    !    Output, integer :: INFO, error flag.  If the user has terminated execution,
    !    INFO is set to the (negative) value of IFLAG. See the description of fcn.
    !    Otherwise, INFO is set as follows:
    !    0, improper input parameters.
    !    1, algorithm estimates that the relative error in the sum of squares
    !       is at most TOL.
    !    2, algorithm estimates that the relative error between X and the
    !       solution is at most TOL.
    !    3, conditions for INFO = 1 and INFO = 2 both hold.
    !    4, FVEC is orthogonal to the columns of the jacobian to machine precision.
    !    5, number of calls to fcn with IFLAG = 1 has reached 100*(N+1).
    !    6, TOL is too small.  No further reduction in the sum of squares is
    !       possible.
    !    7, TOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !
    implicit none

    integer :: ldfjac
    integer :: m
    integer :: n

    real(wp) :: diag(n)
    real(wp) :: factor
    external fcn
    real(wp) :: fjac(ldfjac,n)
    real(wp) :: ftol
    real(wp) :: fvec(m)
    real(wp) :: gtol
    integer :: info
    integer :: ipvt(n)
    integer :: maxfev
    integer :: mode
    integer :: nfev
    integer :: njev
    integer :: nprint
    real(wp) :: qtf(n)
    real(wp) :: tol
    real(wp) :: x(n)
    real(wp) :: xtol

    info = 0

    if ( n <= 0 ) then
        return
    else if ( m < n ) then
        return
    else if ( ldfjac < m ) then
        return
    else if ( tol < 0.0E0_wp ) then
        return
    end if

    factor = 100.0E0_wp
    maxfev = 100*(n + 1)
    ftol = tol
    xtol = tol
    gtol = 1.0E-14_wp
    mode = 1
    nprint = 0

    call lmder ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
        diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )

    if ( info == 8 ) then
        info = 4
    end if

    return
    end
    !-----------------------------------------------------------

    
    subroutine lmdif ( fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
        diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf )

    !*****************************************************************************80
    !
    !! lmdif minimizes M functions in N variables using the Levenberg-Marquardt method.
    !
    !  Discussion:
    !
    !    lmdif minimizes the sum of the squares of M nonlinear functions in
    !    N variables by a modification of the Levenberg-Marquardt algorithm.
    !    The user must provide a subroutine which calculates the functions.
    !    The jacobian is then calculated by a forward-difference approximation.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, external fcn, the name of the user-supplied subroutine which
    !    calculates the functions.  The routine should have the form:
    !
    !      subroutine fcn ( m, n, x, fvec, iflag )
    !      integer :: m
    !      integer :: n
    !
    !      real(wp) fvec(m)
    !      integer :: iflag
    !      real(wp) x(n)
    !
    !    The value of IFLAG should not be changed by fcn unless
    !    the user wants to terminate execution of the routine.
    !    In this case set IFLAG to a negative integer.
    !
    !    Input, integer :: M, the number of functions.
    !
    !    Input, integer :: N, the number of variables.  N must not exceed M.
    !
    !    Input/output, real(wp) :: X(N).  On input, X must contain an initial
    !    estimate of the solution vector.  On output X contains the final
    !    estimate of the solution vector.
    !
    !    Output, real(wp) :: FVEC(M), the functions evaluated at the output X.
    !
    !    Input, real(wp) :: FTOL.  Termination occurs when both the actual
    !    and predicted relative reductions in the sum of squares are at most FTOL.
    !    Therefore, FTOL measures the relative error desired in the sum of
    !    squares.  FTOL should be nonnegative.
    !
    !    Input, real(wp) :: XTOL.  Termination occurs when the relative error
    !    between two consecutive iterates is at most XTOL.  Therefore, XTOL
    !    measures the relative error desired in the approximate solution.  XTOL
    !    should be nonnegative.
    !
    !    Input, real(wp) :: GTOL. termination occurs when the cosine of the
    !    angle between FVEC and any column of the jacobian is at most GTOL in
    !    absolute value.  Therefore, GTOL measures the orthogonality desired
    !    between the function vector and the columns of the jacobian.  GTOL should
    !    be nonnegative.
    !
    !    Input, integer :: MAXFEV.  Termination occurs when the number of calls
    !    to fcn is at least MAXFEV by the end of an iteration.
    !
    !    Input, real(wp) :: EPSFCN, is used in determining a suitable step length for
    !    the forward-difference approximation.  This approximation assumes that
    !    the relative errors in the functions are of the order of EPSFCN.
    !    If EPSFCN is less than the machine precision, it is assumed that the
    !    relative errors in the functions are of the order of the machine
    !    precision.
    !
    !    Input/output, real(wp) :: DIAG(N).  If MODE = 1, then DIAG is set
    !    internally.  If MODE = 2, then DIAG must contain positive entries that
    !    serve as multiplicative scale factors for the variables.
    !
    !    Input, integer :: MODE, scaling option.
    !    1, variables will be scaled internally.
    !    2, scaling is specified by the input DIAG vector.
    !
    !    Input, real(wp) :: FACTOR, determines the initial step bound.  This bound is
    !    set to the product of FACTOR and the euclidean norm of DIAG*X if
    !    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
    !    in the interval (0.1, 100) with 100 the recommended value.
    !
    !    Input, integer :: NPRINT, enables controlled printing of iterates if it
    !    is positive.  In this case, fcn is called with IFLAG = 0 at the
    !    beginning of the first iteration and every NPRINT iterations thereafter
    !    and immediately prior to return, with X and FVEC available
    !    for printing.  If NPRINT is not positive, no special calls
    !    of fcn with IFLAG = 0 are made.
    !
    !    Output, integer :: INFO, error flag.  If the user has terminated execution,
    !    INFO is set to the (negative) value of IFLAG. See the description of fcn.
    !    Otherwise, INFO is set as follows:
    !    0, improper input parameters.
    !    1, both actual and predicted relative reductions in the sum of squares
    !       are at most FTOL.
    !    2, relative error between two consecutive iterates is at most XTOL.
    !    3, conditions for INFO = 1 and INFO = 2 both hold.
    !    4, the cosine of the angle between FVEC and any column of the jacobian
    !       is at most GTOL in absolute value.
    !    5, number of calls to fcn has reached or exceeded MAXFEV.
    !    6, FTOL is too small.  No further reduction in the sum of squares
    !       is possible.
    !    7, XTOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    8, GTOL is too small.  FVEC is orthogonal to the columns of the
    !       jacobian to machine precision.
    !
    !    Output, integer :: NFEV, the number of calls to fcn.
    !
    !    Output, real(wp) :: FJAC(LDFJAC,N), an M by N array.  The upper N by N submatrix
    !    of FJAC contains an upper triangular matrix R with diagonal elements
    !    of nonincreasing magnitude such that
    !
    !      P' * ( JAC' * JAC ) * P = R' * R,
    !
    !    where P is a permutation matrix and JAC is the final calculated jacobian.
    !    Column J of P is column IPVT(J) of the identity matrix.  The lower
    !    trapezoidal part of FJAC contains information generated during
    !    the computation of R.
    !
    !    Input, integer :: LDFJAC, the leading dimension of the array FJAC.
    !    LDFJAC must be at least M.
    !
    !    Output, integer :: IPVT(N), defines a permutation matrix P such that
    !    JAC * P = Q * R, where JAC is the final calculated jacobian, Q is
    !    orthogonal (not stored), and R is upper triangular with diagonal
    !    elements of nonincreasing magnitude.  Column J of P is column IPVT(J)
    !    of the identity matrix.
    !
    !    Output, real(wp) :: QTF(N), the first N elements of Q'*FVEC.
    !
    implicit none

    integer :: ldfjac
    integer :: m
    integer :: n

    real(wp) :: actred
    real(wp) :: delta
    real(wp) :: diag(n)
    real(wp) :: dirder
    real(wp) :: epsfcn
    real(wp) :: epsmch
    real(wp) :: factor
    external fcn
    real(wp) :: fjac(ldfjac,n)
    real(wp) :: fnorm
    real(wp) :: fnorm1
    real(wp) :: ftol
    real(wp) :: fvec(m)
    real(wp) :: gnorm
    real(wp) :: gtol
    integer :: i
    integer :: iflag
    integer :: iter
    integer :: info
    integer :: ipvt(n)
    integer :: j
    integer :: l
    integer :: maxfev
    integer :: mode
    integer :: nfev
    integer :: nprint
    real(wp) :: par
    real(wp) :: pnorm
    real(wp) :: prered
    real(wp) :: qtf(n)
    real(wp) :: ratio
    real(wp) :: sum2
    real(wp) :: temp
    real(wp) :: temp1
    real(wp) :: temp2
    real(wp) :: wa1(n)
    real(wp) :: wa2(n)
    real(wp) :: wa3(n)
    real(wp) :: wa4(m)
    real(wp) :: x(n)
    real(wp) :: xnorm
    real(wp) :: xtol

    epsmch = epsilon ( epsmch )

    info = 0
    iflag = 0
    nfev = 0

    if ( n <= 0 ) then
        go to 300
    else if ( m < n ) then
        go to 300
    else if ( ldfjac < m ) then
        go to 300
    else if ( ftol < 0.0E0_wp ) then
        go to 300
    else if ( xtol < 0.0E0_wp ) then
        go to 300
    else if ( gtol < 0.0E0_wp ) then
        go to 300
    else if ( maxfev <= 0 ) then
        go to 300
    else if ( factor <= 0.0E0_wp ) then
        go to 300
    end if

    if ( mode == 2 ) then
        do j = 1, n
            if ( diag(j) <= 0.0E0_wp ) then
                go to 300
            end if
        end do
    end if
    !
    !  Evaluate the function at the starting point and calculate its norm.
    !
    iflag = 1
    call fcn ( m, n, x, fvec, iflag )
    nfev = 1

    if ( iflag < 0 ) then
        go to 300
    end if

    fnorm = enorm2 ( m, fvec ) !enorm ( m, fvec )
    !
    !  Initialize Levenberg-Marquardt parameter and iteration counter.
    !
    par = 0.0E0_wp
    iter = 1
    !
    !  Beginning of the outer loop.
    !
    30  continue
    !
    !  Calculate the jacobian matrix.
    !
    iflag = 2
    call fdjac2 ( fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn )
    nfev = nfev + n

    if ( iflag < 0 ) then
        go to 300
    end if
    !
    !  If requested, call fcn to enable printing of iterates.
    !
    if ( 0 < nprint ) then
        iflag = 0
        if ( mod ( iter-1, nprint ) == 0 ) then
            call fcn ( m, n, x, fvec, iflag )
        end if
        if ( iflag < 0 ) then
            go to 300
        end if
    end if
    !
    !  Compute the QR factorization of the jacobian.
    !
    call qrfac ( m, n, fjac, ldfjac, .true., ipvt, n, wa1, wa2 )
    !
    !  On the first iteration and if MODE is 1, scale according
    !  to the norms of the columns of the initial jacobian.
    !
    if ( iter == 1 ) then
        if ( mode /= 2 ) then
            diag(1:n) = wa2(1:n)
            do j = 1, n
                if ( abs(wa2(j)) < epsmch ) then
                    diag(j) = 1.0E0_wp
                end if
            end do
        end if
        !
        !  On the first iteration, calculate the norm of the scaled X
        !  and initialize the step bound DELTA.
        !
        wa3(1:n) = diag(1:n) * x(1:n)
        xnorm = enorm2( n, wa3 ) !enorm ( n, wa3 )
        delta = factor * xnorm
        if ( abs(delta) < epsmch) then
            delta = factor
        end if
    end if
    !
    !  Form Q' * FVEC and store the first N components in QTF.
    !
    wa4(1:m) = fvec(1:m)

    do j = 1, n
        if ( fjac(j,j) /= 0.0E0_wp ) then
            sum2 = dot_product ( wa4(j:m), fjac(j:m,j) )
            temp = - sum2 / fjac(j,j)
            wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
        end if
        fjac(j,j) = wa1(j)
        qtf(j) = wa4(j)
    end do
    !
    !  Compute the norm of the scaled gradient.
    !
    gnorm = 0.0E0_wp

    if ( fnorm /= 0.0E0_wp ) then
        do j = 1, n
            l = ipvt(j)
            if ( wa2(l) /= 0.0E0_wp ) then
                sum2 = 0.0E0_wp
                do i = 1, j
                    sum2 = sum2 + fjac(i,j) * ( qtf(i) / fnorm )
                end do
                gnorm = max ( gnorm, abs( sum2 / wa2(l) ) )
            end if
        end do
    end if
    !
    !  Test for convergence of the gradient norm.
    !
    if ( gnorm <= gtol ) then
        info = 4
        go to 300
    end if
    !
    !  Rescale if necessary.
    !
    if ( mode /= 2 ) then
        do j = 1, n
            diag(j) = max ( diag(j), wa2(j) )
        end do
    end if
    !
    !  Beginning of the inner loop.
    !
    200 continue
    !
    !  Determine the Levenberg-Marquardt parameter.
    !
    call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
    !
    !  Store the direction P and X + P.
    !  Calculate the norm of P.
    !
    wa1(1:n) = -wa1(1:n)
    wa2(1:n) = x(1:n) + wa1(1:n)
    wa3(1:n) = diag(1:n) * wa1(1:n)

    pnorm = enorm2 ( n, wa3 ) !enorm ( n, wa3 )
    !
    !  On the first iteration, adjust the initial step bound.
    !
    if ( iter == 1 ) then
        delta = min ( delta, pnorm )
    end if
    !
    !  Evaluate the function at X + P and calculate its norm.
    !
    iflag = 1
    call fcn ( m, n, wa2, wa4, iflag )
    nfev = nfev + 1
    if ( iflag < 0 ) then
        go to 300
    end if
    fnorm1 = enorm2 ( m, wa4 ) !enorm ( m, wa4 )
    !
    !  Compute the scaled actual reduction.
    !
    if ( 0.1E0_wp * fnorm1 < fnorm ) then
        actred = 1.0E0_wp - ( fnorm1 / fnorm )**2
    else
        actred = -1.0E0_wp
    end if
    !
    !  Compute the scaled predicted reduction and the scaled directional derivative.
    !
    do j = 1, n
        wa3(j) = 0.0E0_wp
        l = ipvt(j)
        temp = wa1(l)
        wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
    end do

    temp1 = enorm2 ( n, wa3 ) / fnorm !enorm ( n, wa3 ) / fnorm
    temp2 = ( sqrt ( par ) * pnorm ) / fnorm
    prered = temp1**2 + temp2**2 / 0.5E0_wp
    dirder = - ( temp1**2 + temp2**2 )
    !
    !  Compute the ratio of the actual to the predicted reduction.
    !
    ratio = 0.0E0_wp
    if ( prered /= 0.0E0_wp ) then
        ratio = actred / prered
    end if
    !
    !  Update the step bound.
    !
    if ( ratio <= 0.25E0_wp ) then

        if ( actred >= 0.0E0_wp ) then
            temp = 0.5E0_wp
        endif

        if ( actred < 0.0E0_wp ) then
            temp = 0.5E0_wp * dirder / ( dirder + 0.5E0_wp * actred )
        end if

        if ( 0.1E0_wp * fnorm1 >= fnorm .or. temp < 0.1E0_wp ) then
            temp = 0.1E0_wp
        end if

        delta = temp * min ( delta, pnorm / 0.1E0_wp  )
        par = par / temp

    else

        if (abs(par) < epsmch .or. ratio >= 0.75E0_wp ) then
            delta = 2.0E0_wp * pnorm
            par = 0.5E0_wp * par
        end if

    end if
    !
    !  Test for successful iteration.
    !

    !
    !  Successful iteration. update X, FVEC, and their norms.
    !
    if ( 0.0001E0_wp <= ratio ) then
        x(1:n) = wa2(1:n)
        wa2(1:n) = diag(1:n) * x(1:n)
        fvec(1:m) = wa4(1:m)
        xnorm = enorm2 ( n, wa2 ) !enorm ( n, wa2 )
        fnorm = fnorm1
        iter = iter + 1
    end if
    !
    !  Tests for convergence.
    !
    if ( abs( actred) <= ftol .and. prered <= ftol &
        .and. 0.5E0_wp * ratio <= 1.0E0_wp ) then
    info = 1
    end if

    if ( delta <= xtol * xnorm ) then
        info = 2
    end if

    if ( abs( actred) <= ftol .and. prered <= ftol &
        .and. 0.5E0_wp * ratio <= 1.0E0_wp .and. info == 2 ) info = 3
    if ( info /= 0 ) go to 300
    !
    !  Tests for termination and stringent tolerances.
    !
    if ( nfev >= maxfev ) then
        info = 5
    end if

    if ( abs( actred) <= epsmch .and. prered <= epsmch &
        .and. 0.5E0_wp * ratio <= 1.0E0_wp ) info = 6
    if ( delta <= epsmch * xnorm ) info = 7
    if ( gnorm <= epsmch ) info = 8

    if ( info /= 0 ) then
        go to 300
    end if
    !
    !  End of the inner loop.  Repeat if iteration unsuccessful.
    !
    if ( ratio < 0.0001E0_wp ) then
        go to 200
    end if
    !
    !  End of the outer loop.
    !
    go to 30

    300 continue
    !
    !  Termination, either normal or user imposed.
    !
    if ( iflag < 0 ) then
        info = iflag
    end if

    iflag = 0

    if ( nprint > 0 ) then
        call fcn ( m, n, x, fvec, iflag )
    end if

    end subroutine lmdif
    !-----------------------------------------------------------
       
        
    subroutine lmdif1 ( fcn, m, n, x, fvec, tol, info )
    !*****************************************************************************80
    !
    !! lmdif1 minimizes M functions in N variables using the Levenberg-Marquardt method.
    !
    !  Discussion:
    !
    !    lmdif1 minimizes the sum of the squares of M nonlinear functions in
    !    N variables by a modification of the Levenberg-Marquardt algorithm.
    !    This is done by using the more general least-squares solver lmdif.
    !    The user must provide a subroutine which calculates the functions.
    !    The jacobian is then calculated by a forward-difference approximation.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, external fcn, the name of the user-supplied subroutine which
    !    calculates the functions.  The routine should have the form:
    !
    !      subroutine fcn ( m, n, x, fvec, iflag )
    !      integer :: n
    !      real(wp) fvec(m)
    !      integer :: iflag
    !      real(wp) x(n)
    !
    !    The value of IFLAG should not be changed by fcn unless
    !    the user wants to terminate execution of the routine.
    !    In this case set IFLAG to a negative integer.
    !
    !    Input, integer :: M, the number of functions.
    !
    !    Input, integer :: N, the number of variables.  N must not exceed M.
    !
    !    Input/output, real(wp) :: X(N).  On input, X must contain an initial
    !    estimate of the solution vector.  On output X contains the final
    !    estimate of the solution vector.
    !
    !    Output, real(wp) :: FVEC(M), the functions evaluated at the output X.
    !
    !    Input, real(wp) :: TOL.  Termination occurs when the algorithm
    !    estimates either that the relative error in the sum of squares is at
    !    most TOL or that the relative error between X and the solution is at
    !    most TOL.  TOL should be nonnegative.
    !
    !    Output, integer :: INFO, error flag.  If the user has terminated execution,
    !    INFO is set to the (negative) value of IFLAG. See the description of fcn.
    !    Otherwise, INFO is set as follows:
    !    0, improper input parameters.
    !    1, algorithm estimates that the relative error in the sum of squares
    !       is at most TOL.
    !    2, algorithm estimates that the relative error between X and the
    !       solution is at most TOL.
    !    3, conditions for INFO = 1 and INFO = 2 both hold.
    !    4, FVEC is orthogonal to the columns of the jacobian to machine precision.
    !    5, number of calls to fcn has reached or exceeded 300*(N+1).
    !    6, TOL is too small.  No further reduction in the sum of squares
    !       is possible.
    !    7, TOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !
    implicit none

    integer :: m
    integer :: n

    real(wp) :: diag(n)
    real(wp) :: epsfcn
    real(wp) :: factor
    external fcn
    real(wp) :: fjac(m,n)
    real(wp) :: ftol
    real(wp) :: fvec(m)
    real(wp) :: gtol
    integer :: info
    integer :: ipvt(n)
    integer :: ldfjac
    integer :: maxfev
    integer :: mode
    integer :: nfev
    integer :: nprint
    real(wp) :: qtf(n)
    real(wp) :: tol
    real(wp) :: x(n)
    real(wp) :: xtol

    info = 0

    if ( n <= 0 ) then
        return
    else if ( m < n ) then
        return
    else if ( tol < 0.0E0_wp ) then
        return
    end if

    factor = 100.0E0_wp
    maxfev = 300*(n + 1)
    ftol = tol
    xtol = tol
    gtol = 5.0E0_wp*epsilon(gtol)
    epsfcn = 0.0E0_wp
    mode = 1
    nprint = 0
    ldfjac = m

    call lmdif ( fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
        diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf )

    if ( info == 8 ) then
        info = 4
    end if

    return
    end
    !-----------------------------------------------------------

    
    subroutine lmpar ( n, r, ldr, ipvt, diag, qtb, delta, par, x, sdiag )
    !*****************************************************************************80
    !
    !! lmpar computes a parameter for the Levenberg-Marquardt method.
    !
    !  Discussion:
    !
    !    Given an M by N matrix A, an N by N nonsingular diagonal
    !    matrix D, an M-vector B, and a positive number DELTA,
    !    the problem is to determine a value for the parameter
    !    PAR such that if X solves the system
    !
    !      A*X = B,
    !      sqrt ( PAR ) * D * X = 0,
    !
    !    in the least squares sense, and DXNORM is the euclidean
    !    norm of D*X, then either PAR is zero and
    !
    !      ( DXNORM - DELTA ) <= 0.1 * DELTA,
    !
    !    or PAR is positive and
    !
    !      abs( DXNORM - DELTA) <= 0.1 * DELTA.
    !
    !    This subroutine completes the solution of the problem
    !    if it is provided with the necessary information from the
    !    QR factorization, with column pivoting, of A.  That is, if
    !    A*P = Q*R, where P is a permutation matrix, Q has orthogonal
    !    columns, and R is an upper triangular matrix with diagonal
    !    elements of nonincreasing magnitude, then lmpar expects
    !    the full upper triangle of R, the permutation matrix P,
    !    and the first N components of Q'*B.  On output
    !    lmpar also provides an upper triangular matrix S such that
    !
    !      P' * ( A' * A + PAR * D * D ) * P = S'* S.
    !
    !    S is employed within lmpar and may be of separate interest.
    !
    !    Only a few iterations are generally needed for convergence
    !    of the algorithm.  If, however, the limit of 10 iterations
    !    is reached, then the output PAR will contain the best
    !    value obtained so far.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, integer :: N, the order of R.
    !
    !    Input/output, real(wp) :: R(LDR,N),the N by N matrix.  The full
    !    upper triangle must contain the full upper triangle of the matrix R.
    !    On output the full upper triangle is unaltered, and the strict lower
    !    triangle contains the strict upper triangle (transposed) of the upper
    !    triangular matrix S.
    !
    !    Input, integer :: LDR, the leading dimension of R.  LDR must be
    !    no less than N.
    !
    !    Input, integer :: IPVT(N), defines the permutation matrix P such that
    !    A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
    !
    !    Input, real(wp) :: DIAG(N), the diagonal elements of the matrix D.
    !
    !    Input, real(wp) :: QTB(N), the first N elements of the vector Q'*B.
    !
    !    Input, real(wp) :: DELTA, an upper bound on the euclidean norm of D*X.
    !    DELTA should be positive.
    !
    !    Input/output, real(wp) :: PAR.  On input an initial estimate of the
    !    Levenberg-Marquardt parameter.  On output the final estimate.
    !    PAR should be nonnegative.
    !
    !    Output, real(wp) :: X(N), the least squares solution of the system A*X = B,
    !    sqrt(PAR)*D*X = 0, for the output value of PAR.
    !
    !    Output, real(wp) :: SDIAG(N), the diagonal elements of the upper triangular
    !    matrix S.
    !
    implicit none

    integer :: ldr
    integer :: n

    real(wp) :: delta
    real(wp) :: diag(n)
    real(wp) :: dwarf
    real(wp) :: dxnorm
    real(wp) :: gnorm
    real(wp) :: fp
    integer :: ipvt(n)
    integer :: iter
    integer :: j
    integer :: k
    integer :: l
    integer :: nsing
    real(wp), PARAMETER :: epsmch = 5.0E0_wp*epsilon(epsmch)
    real(wp) :: par
    real(wp) :: parc
    real(wp) :: parl
    real(wp) :: paru
    real(wp) :: qtb(n)
    real(wp) :: r(ldr,n)
    real(wp) :: sdiag(n)
    real(wp) :: sum2
    real(wp) :: temp
    real(wp) :: wa1(n)
    real(wp) :: wa2(n)
    real(wp) :: x(n)
    !
    !  DWARF is the smallest positive magnitude.
    !
    dwarf = tiny ( dwarf )
    !
    !  Compute and store in X the Gauss-Newton direction.
    !
    !  If the jacobian is rank-deficient, obtain a least squares solution.
    !
    nsing = n

    do j = 1, n
        wa1(j) = qtb(j)
        if ( abs(r(j,j)) < epsmch .and. nsing == n ) then
            nsing = j - 1
        end if
        if ( nsing < n ) then
            wa1(j) = 0.0E0_wp
        end if
    end do

    do k = 1, nsing
        j = nsing - k + 1
        wa1(j) = wa1(j) / r(j,j)
        temp = wa1(j)
        wa1(1:j-1) = wa1(1:j-1) - r(1:j-1,j) * temp
    end do

    do j = 1, n
        l = ipvt(j)
        x(l) = wa1(j)
    end do
    !
    !  Initialize the iteration counter.
    !  Evaluate the function at the origin, and test
    !  for acceptance of the Gauss-Newton direction.
    !
    iter = 0
    wa2(1:n) = diag(1:n) * x(1:n)
    dxnorm = enorm2 ( n, wa2 ) !enorm ( n, wa2 )
    fp = dxnorm - delta

    if ( fp <= 0.1E0_wp * delta ) then
        go to 220
    end if
    !
    !  If the jacobian is not rank deficient, the Newton
    !  step provides a lower bound, PARL, for the zero of
    !  the function.
    !
    !  Otherwise set this bound to zero.
    !
    parl = 0.0E0_wp

    if ( nsing >= n ) then
        do j = 1, n
            l = ipvt(j)
            wa1(j) = diag(l) * ( wa2(l) / dxnorm )
        end do
        do j = 1, n
            sum2 = dot_product ( wa1(1:j-1), r(1:j-1,j) )
            wa1(j) = ( wa1(j) - sum2 ) / r(j,j)
        end do
        temp = enorm2 ( n, wa1 ) !enorm ( n, wa1 )
        parl = ( ( fp / delta ) / temp ) / temp
    end if
    !
    !  Calculate an upper bound, PARU, for the zero of the function.
    !
    do j = 1, n
        sum2 = dot_product ( qtb(1:j), r(1:j,j) )
        l = ipvt(j)
        wa1(j) = sum2 / diag(l)
    end do

    gnorm = enorm2 ( n, wa1 ) !enorm ( n, wa1 )
    paru = gnorm / delta
    if (abs(paru) < epsmch) then
        paru = dwarf / min ( delta, 0.1E0_wp )
    end if
    !
    !  If the input PAR lies outside of the interval (PARL, PARU),
    !  set PAR to the closer endpoint.
    !
    par = max ( par, parl )
    par = min ( par, paru )
    if (abs(par) < epsmch) then
        par = gnorm / dxnorm
    end if
    !
    !  Beginning of an iteration.
    !
    150 continue

    iter = iter + 1
    !
    !  Evaluate the function at the current value of PAR.
    !
    if (abs(paru) < epsmch) then
        par = max ( dwarf, 0.001E0_wp * paru )
    end if

    wa1(1:n) = sqrt ( par ) * diag(1:n)

    call qrsolv ( n, r, ldr, ipvt, wa1, qtb, x, sdiag )

    wa2(1:n) = diag(1:n) * x(1:n)
    dxnorm = enorm2 ( n, wa2 ) !enorm ( n, wa2 )
    temp = fp
    fp = dxnorm - delta
    !
    !  If the function is small enough, accept the current value of PAR.
    !
    if ( abs( fp ) <= 0.1E0_wp * delta ) then
        go to 220
    end if
    !
    !  Test for the exceptional cases where PARL
    !  is zero or the number of iterations has reached 10.
    !
    if (abs(parl) < epsmch .and. fp <= temp .and. temp < 0.0E0_wp ) then
        go to 220
    else if ( iter == 10 ) then
        go to 220
    end if
    !
    !  Compute the Newton correction.
    !
    do j = 1, n
        l = ipvt(j)
        wa1(j) = diag(l) * ( wa2(l) / dxnorm )
    end do

    do j = 1, n
        IF (sdiag(j) /= 0.0E0_wp) THEN
            wa1(j) = wa1(j) / sdiag(j)
        ELSE
            sdiag(j) = 1.0E-15_wp   !correction added by A. Zuend to prevent divison by zero here
            wa1(j) = wa1(j) / sdiag(j)
        ENDIF
        temp = wa1(j)
        wa1(j+1:n) = wa1(j+1:n) - r(j+1:n,j) * temp
    end do

    temp = enorm2 ( n, wa1 ) !enorm ( n, wa1 )
    parc = ( ( fp / delta ) / temp ) / temp
    !
    !  Depending on the sign of the function, update PARL or PARU.
    !
    if ( 0.0E0_wp < fp ) then
        parl = max ( parl, par )
    else if ( fp < 0.0E0_wp ) then
        paru = min ( paru, par )
    end if
    !
    !  Compute an improved estimate for PAR.
    !
    par = max ( parl, par + parc )
    !
    !  End of an iteration.
    !
    go to 150

    220 continue
    !
    !  Termination.
    !
    if ( iter == 0 ) then
        par = 0.0E0_wp
    end if

    return
    end
    !-----------------------------------------------------------

    
    subroutine lmstr ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
        diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )

    !*****************************************************************************80
    !
    !! lmstr minimizes M functions in N variables using the Levenberg-Marquardt method.
    !
    !  Discussion:
    !
    !    lmstr minimizes the sum of the squares of M nonlinear functions in
    !    N variables by a modification of the Levenberg-Marquardt algorithm
    !    which uses minimal storage.
    !
    !    The user must provide a subroutine which calculates the functions and
    !    the rows of the jacobian.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, external fcn, the name of the user-supplied subroutine which
    !    calculates the functions and the rows of the jacobian.
    !    fcn should have the form:
    !
    !      subroutine fcn ( m, n, x, fvec, fjrow, iflag )
    !
    !      integer :: m
    !      integer :: n
    !
    !      real fjrow(n)
    !      real(wp) fvec(m)
    !      integer :: iflag
    !      real(wp) x(n)
    !
    !    If the input value of IFLAG is 1, calculate the functions at X and
    !    return this vector in FVEC.
    !    If the input value of IFLAG is I > 1, calculate the (I-1)-st row of
    !    the jacobian at X, and return this vector in FJROW.
    !    To terminate the algorithm, set the output value of IFLAG negative.
    !
    !    Input, integer :: M, the number of functions.
    !
    !    Input, integer :: N, the number of variables.  N must not exceed M.
    !
    !    Input/output, real(wp) :: X(N).  On input, X must contain an initial
    !    estimate of the solution vector.  On output X contains the final
    !    estimate of the solution vector.
    !
    !    Output, real(wp) :: FVEC(M), the functions evaluated at the output X.
    !
    !    Output, real(wp) :: FJAC(LDFJAC,N), an N by N array.  The upper
    !    triangle of FJAC contains an upper triangular matrix R such that
    !
    !      P' * ( JAC' * JAC ) * P = R' * R,
    !
    !    where P is a permutation matrix and JAC is the final calculated jacobian.
    !    Column J of P is column IPVT(J) of the identity matrix.  The lower
    !    triangular part of FJAC contains information generated during
    !    the computation of R.
    !
    !    Input, integer :: LDFJAC, the leading dimension of the array FJAC.
    !    LDFJAC must be at least N.
    !
    !    Input, real(wp) :: FTOL.  Termination occurs when both the actual and
    !    predicted relative reductions in the sum of squares are at most FTOL.
    !    Therefore, FTOL measures the relative error desired in the sum of
    !    squares.  FTOL should be nonnegative.
    !
    !    Input, real(wp) :: XTOL.  Termination occurs when the relative error between
    !    two consecutive iterates is at most XTOL.  XTOL should be nonnegative.
    !
    !    Input, real(wp) :: GTOL. termination occurs when the cosine of the angle
    !    between FVEC and any column of the jacobian is at most GTOL in absolute
    !    value.  Therefore, GTOL measures the orthogonality desired between the
    !    function vector and the columns of the jacobian.  GTOL should
    !    be nonnegative.
    !
    !    Input, integer :: MAXFEV.  Termination occurs when the number of calls
    !    to fcn with IFLAG = 1 is at least MAXFEV by the end of an iteration.
    !
    !    Input/output, real(wp) :: DIAG(N).  If MODE = 1, then DIAG is set internally.
    !    If MODE = 2, then DIAG must contain positive entries that serve as
    !    multiplicative scale factors for the variables.
    !
    !    Input, integer :: MODE, scaling option.
    !    1, variables will be scaled internally.
    !    2, scaling is specified by the input DIAG vector.
    !
    !    Input, real(wp) :: FACTOR, determines the initial step bound.  This bound is
    !    set to the product of FACTOR and the euclidean norm of DIAG*X if
    !    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
    !    in the interval (0.1, 100) with 100 the recommended value.
    !
    !    Input, integer :: NPRINT, enables controlled printing of iterates if it
    !    is positive.  In this case, fcn is called with IFLAG = 0 at the
    !    beginning of the first iteration and every NPRINT iterations thereafter
    !    and immediately prior to return, with X and FVEC available
    !    for printing.  If NPRINT is not positive, no special calls
    !    of fcn with IFLAG = 0 are made.
    !
    !    Output, integer :: INFO, error flag.  If the user has terminated execution,
    !    INFO is set to the (negative) value of IFLAG. See the description of fcn.
    !    Otherwise, INFO is set as follows:
    !    0, improper input parameters.
    !    1, both actual and predicted relative reductions in the sum of squares
    !       are at most FTOL.
    !    2, relative error between two consecutive iterates is at most XTOL.
    !    3, conditions for INFO = 1 and INFO = 2 both hold.
    !    4, the cosine of the angle between FVEC and any column of the jacobian
    !       is at most GTOL in absolute value.
    !    5, number of calls to fcn with IFLAG = 1 has reached MAXFEV.
    !    6, FTOL is too small.  No further reduction in the sum of squares is
    !       possible.
    !    7, XTOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    8, GTOL is too small.  FVEC is orthogonal to the columns of the
    !       jacobian to machine precision.
    !
    !    Output, integer :: NFEV, the number of calls to fcn with IFLAG = 1.
    !
    !    Output, integer :: NJEV, the number of calls to fcn with IFLAG = 2.
    !
    !    Output, integer :: IPVT(N), defines a permutation matrix P such that
    !    JAC * P = Q * R, where JAC is the final calculated jacobian, Q is
    !    orthogonal (not stored), and R is upper triangular.
    !    Column J of P is column IPVT(J) of the identity matrix.
    !
    !    Output, real(wp) :: QTF(N), contains the first N elements of Q'*FVEC.
    !
    implicit none

    integer :: ldfjac
    integer :: m
    integer :: n

    real(wp) :: actred
    real(wp) :: delta
    real(wp) :: diag(n)
    real(wp) :: dirder
    real(wp) :: epsmch
    real(wp) :: factor
    external fcn
    real(wp) :: fjac(ldfjac,n)
    real(wp) :: fnorm
    real(wp) :: fnorm1
    real(wp) :: ftol
    real(wp) :: fvec(m)
    real(wp) :: gnorm
    real(wp) :: gtol
    integer :: i
    integer :: iflag
    integer :: info
    integer :: ipvt(n)
    integer :: iter
    integer :: j
    integer :: l
    integer :: maxfev
    integer :: mode
    integer :: nfev
    integer :: njev
    integer :: nprint
    real(wp) :: par
    real(wp) :: pnorm
    real(wp) :: prered
    real(wp) :: qtf(n)
    real(wp) :: ratio
    logical :: sing
    real(wp) :: sum2
    real(wp) :: temp
    real(wp) :: temp1
    real(wp) :: temp2
    real(wp) :: wa1(n)
    real(wp) :: wa2(n)
    real(wp) :: wa3(n)
    real(wp) :: wa4(m)
    real(wp) :: x(n)
    real(wp) :: xnorm
    real(wp) :: xtol

    epsmch = epsilon ( epsmch )

    info = 0
    iflag = 0
    nfev = 0
    njev = 0
    !
    !  Check the input parameters for errors.
    !
    if ( n <= 0 ) then
        go to 340
    else if ( m < n ) then
        go to 340
    else if ( ldfjac < n ) then
        go to 340
    else if ( ftol < 0.0E0_wp ) then
        go to 340
    else if ( xtol < 0.0E0_wp ) then
        go to 340
    else if ( gtol < 0.0E0_wp ) then
        go to 340
    else if ( maxfev <= 0 ) then
        go to 340
    else if ( factor <= 0.0E0_wp ) then
        go to 340
    end if

    if ( mode == 2 ) then
        do j = 1, n
            if ( diag(j) <= 0.0E0_wp ) then
                go to 340
            end if
        end do
    end if
    !
    !  Evaluate the function at the starting point and calculate its norm.
    !
    iflag = 1
    call fcn ( m, n, x, fvec, wa3, iflag )
    nfev = 1
    if ( iflag < 0 ) go to 340
    fnorm = enorm2 ( m, fvec ) !enorm ( m, fvec )
    !
    !  Initialize Levenberg-Marquardt parameter and iteration counter.
    !
    par = 0.0E0_wp
    iter = 1
    !
    !  Beginning of the outer loop.
    !
    30  continue
    !
    !  If requested, call fcn to enable printing of iterates.
    !
    if ( 0 < nprint ) then
        iflag = 0
        if ( mod ( iter-1, nprint ) == 0 ) then
            call fcn ( m, n, x, fvec, wa3, iflag )
        end if
        if ( iflag < 0 ) go to 340
    end if
    !
    !  Compute the QR factorization of the jacobian matrix calculated one row
    !  at a time, while simultaneously forming Q'* FVEC and storing
    !  the first N components in QTF.
    !
    qtf(1:n) = 0.0E0_wp
    fjac(1:n,1:n) = 0.0E0_wp
    iflag = 2

    do i = 1, m
        call fcn ( m, n, x, fvec, wa3, iflag )
        if ( iflag < 0 ) go to 340
        temp = fvec(i)
        call rwupdt ( n, fjac, ldfjac, wa3, qtf, temp, wa1, wa2 )
        iflag = iflag + 1
    end do

    njev = njev + 1
    !
    !  If the jacobian is rank deficient, call qrfac to
    !  reorder its columns and update the components of QTF.
    !
    sing = .false.
    do j = 1, n
        if ( abs(fjac(j,j)) < epsmch ) sing = .true.
        ipvt(j) = j
        wa2(j) = enorm2 ( j, fjac(1,j) ) !enorm ( j, fjac(1,j) )
    end do

    if (sing) then
        call qrfac ( n, n, fjac, ldfjac, .true., ipvt, n, wa1, wa2 )
        do j = 1, n
            if ( fjac(j,j) /= 0.0E0_wp ) then
                sum2 = dot_product ( qtf(j:n), fjac(j:n,j) )
                temp = - sum2 / fjac(j,j)
                qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp
            end if
            fjac(j,j) = wa1(j)
        end do
    end if
    !
    !  On the first iteration
    !    if mode is 1,
    !      scale according to the norms of the columns of the initial jacobian.
    !    calculate the norm of the scaled X,
    !    initialize the step bound delta.
    !
    if ( iter == 1 ) then
        if ( mode /= 2 ) then
            diag(1:n) = wa2(1:n)
            do j = 1, n
                if ( abs(wa2(j)) < epsmch ) then
                    diag(j) = 1.0E0_wp
                end if
            end do
        end if
        wa3(1:n) = diag(1:n) * x(1:n)
        xnorm = enorm2 ( n, wa3 ) !enorm ( n, wa3 )
        delta = factor * xnorm
        if ( abs(delta) < epsmch) then
            delta = factor
        end if
    end if
    !
    !  Compute the norm of the scaled gradient.
    !
    gnorm = 0.0E0_wp

    if ( fnorm /= 0.0E0_wp ) then
        do j = 1, n
            l = ipvt(j)
            if ( wa2(l) /= 0.0E0_wp ) then
                sum2 = dot_product ( qtf(1:j), fjac(1:j,j) ) / fnorm
                gnorm = max ( gnorm, abs( sum2 / wa2(l) ) )
            end if
        end do
    end if
    !
    !  Test for convergence of the gradient norm.
    !
    if ( gnorm <= gtol ) then
        info = 4
        go to 340
    end if
    !
    !  Rescale if necessary.
    !
    if ( mode /= 2 ) then
        do j = 1, n
            diag(j) = max ( diag(j), wa2(j) )
        end do
    end if
    !
    !  Beginning of the inner loop.
    !
    240 continue
    !
    !  Determine the Levenberg-Marquardt parameter.
    !
    call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
    !
    !  Store the direction P and X + P.
    !  Calculate the norm of P.
    !
    wa1(1:n) = -wa1(1:n)
    wa2(1:n) = x(1:n) + wa1(1:n)
    wa3(1:n) = diag(1:n) * wa1(1:n)
    pnorm = enorm2 ( n, wa3 ) !enorm ( n, wa3 )
    !
    !  On the first iteration, adjust the initial step bound.
    !
    if ( iter == 1 ) then
        delta = min ( delta, pnorm )
    end if
    !
    !  Evaluate the function at X + P and calculate its norm.
    !
    iflag = 1
    call fcn ( m, n, wa2, wa4, wa3, iflag )
    nfev = nfev + 1
    if ( iflag < 0 ) go to 340
    fnorm1 = enorm2 ( m, wa4 ) !enorm ( m, wa4 )
    !
    !  Compute the scaled actual reduction.
    !
    if ( 0.1E0_wp * fnorm1 < fnorm ) then
        actred = 1.0E0_wp - ( fnorm1 / fnorm )**2
    else
        actred = -1.0E0_wp
    end if
    !
    !  Compute the scaled predicted reduction and
    !  the scaled directional derivative.
    !
    do j = 1, n
        wa3(j) = 0.0E0_wp
        l = ipvt(j)
        temp = wa1(l)
        wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
    end do

    temp1 = enorm2 ( n, wa3 ) / fnorm !enorm ( n, wa3 ) / fnorm
    temp2 = ( sqrt(par) * pnorm ) / fnorm
    prered = temp1**2 + temp2**2 / 0.5E0_wp
    dirder = - ( temp1**2 + temp2**2 )
    !
    !  Compute the ratio of the actual to the predicted reduction.
    !
    if ( prered /= 0.0E0_wp ) then
        ratio = actred / prered
    else
        ratio = 0.0E0_wp
    end if
    !
    !  Update the step bound.
    !
    if ( ratio <= 0.25E0_wp ) then
        if ( actred >= 0.0E0_wp ) then
            temp = 0.5E0_wp
        else
            temp = 0.5E0_wp * dirder / ( dirder + 0.5E0_wp * actred )
        end if
        if ( 0.1E0_wp * fnorm1 >= fnorm .or. temp < 0.1E0_wp ) then
            temp = 0.1E0_wp
        end if
        delta = temp * min ( delta, pnorm / 0.1E0_wp )
        par = par / temp
    else
        if (abs(par) < epsmch .or. ratio >= 0.75E0_wp ) then
            delta = pnorm / 0.5E0_wp
            par = 0.5E0_wp * par
        end if
    end if
    !
    !  Test for successful iteration.
    !
    if ( ratio >= 0.0001E0_wp ) then
        x(1:n) = wa2(1:n)
        wa2(1:n) = diag(1:n) * x(1:n)
        fvec(1:m) = wa4(1:m)
        xnorm = enorm2 ( n, wa2 ) !enorm ( n, wa2 )
        fnorm = fnorm1
        iter = iter + 1
    end if
    !
    !  Tests for convergence, termination and stringent tolerances.
    !
    if ( abs( actred ) <= ftol .and. prered <= ftol &
        .and. 0.5E0_wp * ratio <= 1.0E0_wp ) then
    info = 1
    end if

    if ( delta <= xtol * xnorm ) then
        info = 2
    end if

    if ( abs( actred ) <= ftol .and. prered <= ftol &
        .and. 0.5E0_wp * ratio <= 1.0E0_wp .and. info == 2 ) then
    info = 3
    end if

    if ( info /= 0 ) then
        go to 340
    end if

    if ( nfev >= maxfev) then
        info = 5
    end if

    if ( abs( actred ) <= epsmch .and. prered <= epsmch &
        .and. 0.5E0_wp * ratio <= 1.0E0_wp ) info = 6
    if ( delta <= epsmch * xnorm ) info = 7
    if ( gnorm <= epsmch ) info = 8

    if ( info /= 0 ) then
        go to 340
    end if
    !
    !  End of the inner loop.  Repeat if iteration unsuccessful.
    !
    if ( ratio < 0.0001E0_wp ) then
        go to 240
    end if
    !
    !  End of the outer loop.
    !
    go to 30

    340 continue
    !
    !  Termination, either normal or user imposed.
    !
    if ( iflag < 0 ) then
        info = iflag
    end if

    iflag = 0

    if ( nprint > 0 ) then
        call fcn ( m, n, x, fvec, wa3, iflag )
    end if

    return
    end
    !-----------------------------------------------------------

        
    subroutine lmstr1 ( fcn, m, n, x, fvec, fjac, ldfjac, tol, info )
    !*****************************************************************************80
    !
    !! lmstr1 minimizes M functions in N variables using the Levenberg-Marquardt method.
    !
    !  Discussion:
    !
    !    lmstr1 minimizes the sum of the squares of M nonlinear functions in
    !    N variables by a modification of the Levenberg-Marquardt algorithm
    !    which uses minimal storage.
    !
    !    This is done by using the more general least-squares solver
    !    lmstr.  The user must provide a subroutine which calculates
    !    the functions and the rows of the jacobian.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, external fcn, the name of the user-supplied subroutine which
    !    calculates the functions and the rows of the jacobian.
    !    fcn should have the form:
    !
    !      subroutine fcn ( m, n, x, fvec, fjrow, iflag )
    !      integer :: m,n,iflag
    !      integer :: n
    !      real fjrow(n)
    !      real(wp) fvec(m)
    !      integer :: iflag
    !      real(wp) x(n)
    !
    !    If the input value of IFLAG is 1, calculate the functions at X and
    !    return this vector in FVEC.
    !    If the input value of IFLAG is I > 1, calculate the (I-1)-st row of
    !    the jacobian at X, and return this vector in FJROW.
    !    To terminate the algorithm, set the output value of IFLAG negative.
    !
    !    Input, integer :: M, the number of functions.
    !
    !    Input, integer :: N, the number of variables.  N must not exceed M.
    !
    !    Input/output, real(wp) :: X(N).  On input, X must contain an initial
    !    estimate of the solution vector.  On output X contains the final
    !    estimate of the solution vector.
    !
    !    Output, real(wp) :: FVEC(M), the functions evaluated at the output X.
    !
    !    Output, real(wp) :: FJAC(LDFJAC,N), an N by N array.  The upper
    !    triangle contains an upper triangular matrix R such that
    !
    !      P' * ( JAC' * JAC ) * P = R' * R,
    !
    !    where P is a permutation matrix and JAC is the final calculated
    !    jacobian.  Column J of P is column IPVT(J) of the identity matrix.
    !    The lower triangular part of FJAC contains information generated
    !    during the computation of R.
    !
    !    Input, integer :: LDFJAC, the leading dimension of the array FJAC.
    !    LDFJAC must be at least N.
    !
    !    Input, real(wp) :: TOL. Termination occurs when the algorithm estimates
    !    either that the relative error in the sum of squares is at most TOL
    !    or that the relative error between X and the solution is at most TOL.
    !    TOL should be nonnegative.
    !
    !    Output, integer :: INFO, error flag.  If the user has terminated execution,
    !    INFO is set to the (negative) value of IFLAG. See the description of fcn.
    !    Otherwise, INFO is set as follows:
    !    0, improper input parameters.
    !    1, algorithm estimates that the relative error in the sum of squares
    !       is at most TOL.
    !    2, algorithm estimates that the relative error between X and the
    !       solution is at most TOL.
    !    3, conditions for INFO = 1 and INFO = 2 both hold.
    !    4, FVEC is orthogonal to the columns of the jacobian to machine precision.
    !    5, number of calls to fcn with IFLAG = 1 has reached 100*(N+1).
    !    6, TOL is too small.  No further reduction in the sum of squares
    !       is possible.
    !    7, TOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !
    implicit none

    integer :: ldfjac
    integer :: m
    integer :: n

    real(wp) :: diag(n)
    real(wp) :: factor
    external fcn
    real(wp) :: fjac(ldfjac,n)
    real(wp) :: ftol
    real(wp) :: fvec(m)
    real(wp) :: gtol
    integer :: info
    integer :: ipvt(n)
    integer :: maxfev
    integer :: mode
    integer :: nfev
    integer :: njev
    integer :: nprint
    real(wp) :: qtf(n)
    real(wp) :: tol
    real(wp) :: x(n)
    real(wp) :: xtol

    info = 0

    if ( n <= 0 ) then
        return
    else if ( m < n ) then
        return
    else if ( ldfjac < n ) then
        return
    else if ( tol < 0.0E0_wp ) then
        return
    end if

    factor = 100.0E0_wp
    maxfev = 100 * ( n + 1 )
    ftol = tol
    xtol = tol
    gtol = 1.0E-14_wp
    mode = 1
    nprint = 0

    call lmstr ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
        diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )

    if ( info == 8 ) then
        info = 4
    end if

    return
    end
    !-----------------------------------------------------------

    
    subroutine qform ( m, n, q, ldq )
    !*****************************************************************************80
    !
    !! qform produces the explicit QR factorization of a matrix.
    !
    !  Discussion:
    !
    !    The QR factorization of a matrix is usually accumulated in implicit
    !    form, that is, as a series of orthogonal transformations of the
    !    original matrix.  This routine carries out those transformations,
    !    to explicitly exhibit the factorization construced by qrfac.
    !
    !  Parameters:
    !
    !    Input, integer :: M, is a positive integer :: input variable set to the number
    !    of rows of A and the order of Q.
    !
    !    Input, integer :: N, is a positive integer :: input variable set to the number
    !    of columns of A.
    !
    !    Input/output, real(wp) :: Q(LDQ,M).  Q is an M by M array.
    !    On input the full lower trapezoid in the first min(M,N) columns of Q
    !    contains the factored form.
    !    On output, Q has been accumulated into a square matrix.
    !
    !    Input, integer :: LDQ, is a positive integer :: input variable not less
    !    than M which specifies the leading dimension of the array Q.
    !
    implicit none

    integer :: ldq
    integer :: m
    integer :: n

    integer :: j
    integer :: k
    integer :: l
    integer :: minmn
    real(wp) :: q(ldq,m)
    real(wp) :: temp
    real(wp) :: wa(m)

    minmn = min ( m, n )

    do j = 2, minmn
        q(1:j-1,j) = 0.0E0_wp
    end do
    !
    !  Initialize remaining columns to those of the identity matrix.
    !
    q(1:m,n+1:m) = 0.0E0_wp

    do j = n+1, m
        q(j,j) = 1.0E0_wp
    end do
    !
    !  Accumulate Q from its factored form.
    !
    do l = 1, minmn

        k = minmn - l + 1

        wa(k:m) = q(k:m,k)

        q(k:m,k) = 0.0E0_wp
        q(k,k) = 1.0E0_wp

        if ( wa(k) /= 0.0E0_wp ) then

            do j = k, m
                temp = dot_product ( wa(k:m), q(k:m,j) ) / wa(k)
                q(k:m,j) = q(k:m,j) - temp * wa(k:m)
            end do

        end if

    end do

    return
    end
    !-----------------------------------------------------------

    
    subroutine qrfac ( m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm )
    !*****************************************************************************80
    !
    !! qrfac computes a QR factorization using Householder transformations.
    !
    !  Discussion:
    !
    !    This subroutine uses Householder transformations with column
    !    pivoting (optional) to compute a QR factorization of the
    !    M by N matrix A.  That is, qrfac determines an orthogonal
    !    matrix Q, a permutation matrix P, and an upper trapezoidal
    !    matrix R with diagonal elements of nonincreasing magnitude,
    !    such that A*P = Q*R.  The Householder transformation for
    !    column K, K = 1,2,...,min(M,N), is of the form
    !
    !      I - ( 1 / U(K) ) * U * U'
    !
    !    where U has zeros in the first K-1 positions.  The form of
    !    this transformation and the method of pivoting first
    !    appeared in the corresponding LINPACK subroutine.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, integer :: M, the number of rows of A.
    !
    !    Input, integer :: N, the number of columns of A.
    !
    !    Input/output, real(wp) :: A(LDA,N), the M by N array.
    !    On input, A contains the matrix for which the QR factorization is to
    !    be computed.  On output, the strict upper trapezoidal part of A contains
    !    the strict upper trapezoidal part of R, and the lower trapezoidal
    !    part of A contains a factored form of Q (the non-trivial elements of
    !    the U vectors described above).
    !
    !    Input, integer :: LDA, the leading dimension of A, which must
    !    be no less than M.
    !
    !    Input, logical :: PIVOT, is TRUE if column pivoting is to be carried out.
    !
    !    Output, integer :: IPVT(LIPVT), defines the permutation matrix P such
    !    that A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
    !    If PIVOT is false, IPVT is not referenced.
    !
    !    Input, integer :: LIPVT, the dimension of IPVT, which should be N if
    !    pivoting is used.
    !
    !    Output, real(wp) :: RDIAG(N), contains the diagonal elements of R.
    !
    !    Output, real(wp) :: ACNORM(N), the norms of the corresponding
    !    columns of the input matrix A.  If this information is not needed,
    !    then ACNORM can coincide with RDIAG.
    !
    implicit none

    integer :: lda
    integer :: lipvt
    integer :: n

    real(wp) :: a(lda,n)
    real(wp) :: acnorm(n)
    real(wp) :: ajnorm
    real(wp) :: epsmch
    integer :: i
    integer :: ipvt(lipvt)
    integer :: j
    integer :: k
    integer :: kmax
    integer :: m
    integer :: minmn
    logical :: pivot
    real(wp) :: rdiag(n)
    real(wp) :: temp
    real(wp) :: wa(n)

    epsmch = epsilon ( epsmch )
    !
    !  Compute the initial column norms and initialize several arrays.
    !
    IF (m < 1) THEN
        WRITE(*,*) "WARNING: m is < 1 in minpack2.f90, line 4486!! m: ",m
        READ(*,*)
        STOP
    ENDIF
    do j = 1, n
        WHERE (a(1:m,j) > 1.0E100_wp) !exception handling (by AZ)
            a(1:m,j) = 1.0E100_wp
        ELSEWHERE (a(1:m,j) < -1.0E100_wp)
            a(1:m,j) = -1.0E100_wp
        ENDWHERE
        acnorm(j) = enorm2(m, a(1:m,j)) !enorm2(m, a(1:m,j)) !enorm ( m, a(1:m,j) )
    end do

    rdiag(1:n) = acnorm(1:n)
    wa(1:n) = acnorm(1:n)

    if ( pivot ) then
        do j = 1, n
            ipvt(j) = j
        end do
    end if
    !
    !  Reduce A to R with Householder transformations.
    !
    minmn = min ( m, n )

    do j = 1, minmn
        !
        !  Bring the column of largest norm into the pivot position.
        !
        if ( pivot ) then

            kmax = j

            do k = j, n
                if ( rdiag(k) > rdiag(kmax) ) then
                    kmax = k
                end if
            end do

            if ( kmax /= j ) then

                do i = 1, m
                    call d_swap( a(i,j), a(i,kmax) )
                end do

                rdiag(kmax) = rdiag(j)
                wa(kmax) = wa(j)

                call i_swap( ipvt(j), ipvt(kmax) )

            end if

        end if
        !
        !  Compute the Householder transformation to reduce the
        !  J-th column of A to a multiple of the J-th unit vector.
        !
        ajnorm = enorm2(m-j+1, a(j,j)) !enorm2 ( m-j+1, a(j,j) ) !enorm ( m-j+1, a(j,j) )

        if ( ajnorm /= 0.0E0_wp ) then

            if ( a(j,j) < 0.0E0_wp ) then
                ajnorm = -ajnorm
            end if

            a(j:m,j) = a(j:m,j) / ajnorm
            a(j,j) = a(j,j) + 1.0E0_wp
            !
            !  Apply the transformation to the remaining columns and update the norms.
            !
            do k = j+1, n

                temp = dot_product ( a(j:m,j), a(j:m,k) ) / a(j,j)

                a(j:m,k) = a(j:m,k) - temp * a(j:m,j)

                if ( pivot .and. rdiag(k) /= 0.0E0_wp ) then

                    temp = a(j,k) / rdiag(k)
                    rdiag(k) = rdiag(k) * sqrt ( max ( 0.0E0_wp, 1.0E0_wp-temp**2 ) )

                    if ( 0.05E0_wp * ( rdiag(k) / wa(k) )**2 <= epsmch ) then
                        rdiag(k) = enorm2(m-j, a(j+1,k))
                        wa(k) = rdiag(k)
                    end if

                end if

            end do

        end if

        rdiag(j) = -ajnorm

    end do

    return
    end
    !-----------------------------------------------------------

    
    subroutine qrsolv ( n, r, ldr, ipvt, diag, qtb, x, sdiag )
    !*****************************************************************************80
    !
    !! qrsolv solves a rectangular linear system A*x=b in the least squares sense.
    !
    !  Discussion:
    !
    !    Given an M by N matrix A, an N by N diagonal matrix D,
    !    and an M-vector B, the problem is to determine an X which
    !    solves the system
    !
    !      A*X = B
    !      D*X = 0
    !
    !    in the least squares sense.
    !
    !    This subroutine completes the solution of the problem
    !    if it is provided with the necessary information from the
    !    QR factorization, with column pivoting, of A.  That is, if
    !    Q*P = Q*R, where P is a permutation matrix, Q has orthogonal
    !    columns, and R is an upper triangular matrix with diagonal
    !    elements of nonincreasing magnitude, then qrsolv expects
    !    the full upper triangle of R, the permutation matrix p,
    !    and the first N components of Q'*B.
    !
    !    The system is then equivalent to
    !
    !      R*Z = Q'*B
    !      P'*D*P*Z = 0
    !
    !    where X = P*Z.  If this system does not have full rank,
    !    then a least squares solution is obtained.  On output qrsolv
    !    also provides an upper triangular matrix S such that
    !
    !      P'*(A'*A + D*D)*P = S'*S.
    !
    !    S is computed within qrsolv and may be of separate interest.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, integer :: N, the order of R.
    !
    !    Input/output, real(wp) :: R(LDR,N), the N by N matrix.
    !    On input the full upper triangle must contain the full upper triangle
    !    of the matrix R.  On output the full upper triangle is unaltered, and
    !    the strict lower triangle contains the strict upper triangle
    !    (transposed) of the upper triangular matrix S.
    !
    !    Input, integer :: LDR, the leading dimension of R, which must be
    !    at least N.
    !
    !    Input, integer :: IPVT(N), defines the permutation matrix P such that
    !    A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
    !
    !    Input, real(wp) :: DIAG(N), the diagonal elements of the matrix D.
    !
    !    Input, real(wp) :: QTB(N), the first N elements of the vector Q'*B.
    !
    !    Output, real(wp) :: X(N), the least squares solution.
    !
    !    Output, real(wp) :: SDIAG(N), the diagonal elements of the upper
    !    triangular matrix S.
    !
    implicit none

    integer :: ldr
    integer :: n

    real(wp) :: c
    real(wp) :: cotan
    real(wp) :: diag(n)
    integer :: i
    integer :: ipvt(n)
    integer :: j
    integer :: k
    integer :: l
    integer :: nsing
    real(wp), parameter :: epsmch = 5.0_wp*epsilon(epsmch)
    real(wp) :: qtb(n)
    real(wp) :: qtbpj
    real(wp) :: r(ldr,n)
    real(wp) :: s
    real(wp) :: sdiag(n)
    real(wp) :: sum2
    real(wp) :: t
    real(wp) :: temp
    real(wp) :: wa(n)
    real(wp) :: x(n)
    !
    !  Copy R and Q'*B to preserve input and initialize S.
    !
    !  In particular, save the diagonal elements of R in X.
    !
    do j = 1, n
        r(j:n,j) = r(j,j:n)
        x(j) = r(j,j)
    end do

    wa(1:n) = qtb(1:n)
    !
    !  Eliminate the diagonal matrix D using a Givens rotation.
    !
    do j = 1, n
        !
        !  Prepare the row of D to be eliminated, locating the
        !  diagonal element using P from the QR factorization.
        !
        l = ipvt(j)

        if ( diag(l) /= 0.0E0_wp ) then

            sdiag(j:n) = 0.0E0_wp
            sdiag(j) = diag(l)
            !
            !  The transformations to eliminate the row of D
            !  modify only a single element of Q'*B
            !  beyond the first N, which is initially zero.
            !
            qtbpj = 0.0E0_wp

            do k = j, n
                !
                !  Determine a Givens rotation which eliminates the
                !  appropriate element in the current row of D.
                !
                if ( sdiag(k) /= 0.0E0_wp ) then

                    if ( abs( r(k,k) ) < abs( sdiag(k) ) ) then
                        cotan = r(k,k) / sdiag(k)
                        s = 0.5E0_wp / sqrt ( 0.25E0_wp + 0.25E0_wp * cotan**2 )
                        c = s * cotan
                    else
                        t = sdiag(k) / r(k,k)
                        c = 0.5E0_wp / sqrt ( 0.25E0_wp + 0.25E0_wp * t**2 )
                        s = c * t
                    end if
                    !
                    !  Compute the modified diagonal element of R and
                    !  the modified element of (Q'*B,0).
                    !
                    r(k,k) = c * r(k,k) + s * sdiag(k)
                    temp = c * wa(k) + s * qtbpj
                    qtbpj = - s * wa(k) + c * qtbpj
                    wa(k) = temp
                    !
                    !  Accumulate the tranformation in the row of S.
                    !
                    do i = k+1, n
                        temp = c * r(i,k) + s * sdiag(i)
                        sdiag(i) = - s * r(i,k) + c * sdiag(i)
                        r(i,k) = temp
                    end do

                end if

            end do

        end if
        !
        !  Store the diagonal element of S and restore
        !  the corresponding diagonal element of R.
        !
        sdiag(j) = r(j,j)
        r(j,j) = x(j)

    end do
    !
    !  Solve the triangular system for Z.  If the system is
    !  singular, then obtain a least squares solution.
    !
    nsing = n

    do j = 1, n

        if (abs(sdiag(j)) < epsmch .and. nsing == n ) then
            nsing = j - 1
        end if

        if ( nsing < n ) then
            wa(j) = 0.0E0_wp
        end if

    end do

    do j = nsing, 1, -1
        sum2 = dot_product ( wa(j+1:nsing), r(j+1:nsing,j) )
        wa(j) = ( wa(j) - sum2 ) / sdiag(j)
    end do
    !
    !  Permute the components of Z back to components of X.
    !
    do j = 1, n
        l = ipvt(j)
        x(l) = wa(j)
    end do

    return
    end
    !-----------------------------------------------------------

    
    subroutine rwupdt ( n, r, ldr, w, b, alpha, c, s )
    !*****************************************************************************80
    !
    !! rwupdt computes the decomposition of a triangular matrix augmented by one row.
    !
    !  Discussion:
    !
    !    Given an N by N upper triangular matrix R, this subroutine
    !    computes the QR decomposition of the matrix formed when a row
    !    is added to R.  If the row is specified by the vector W, then
    !    rwupdt determines an orthogonal matrix Q such that when the
    !    N+1 by N matrix composed of R augmented by W is premultiplied
    !    by Q', the resulting matrix is upper trapezoidal.
    !    The matrix Q' is the product of N transformations
    !
    !      G(N)*G(N-1)* ... *G(1)
    !
    !    where G(I) is a Givens rotation in the (I,N+1) plane which eliminates
    !    elements in the (N+1)-st plane.  rwupdt also computes the product
    !    Q'*C where C is the (N+1)-vector (B,ALPHA).  Q itself is not
    !    accumulated, rather the information to recover the G rotations is
    !    supplied.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, integer :: N, the order of R.
    !
    !    Input/output, real(wp) :: R(LDR,N).  On input the upper triangular
    !    part of R must contain the matrix to be updated.  On output R contains the
    !    updated triangular matrix.
    !
    !    Input, integer :: LDR, the leading dimension of the array R.
    !    LDR must not be less than N.
    !
    !    Input, real(wp) :: W(N), the row vector to be added to R.
    !
    !    Input/output, real(wp) :: B(N).  On input, the first N elements
    !    of the vector C.  On output the first N elements of the vector Q'*C.
    !
    !    Input/output, real(wp) :: ALPHA.  On input, the (N+1)-st element
    !    of the vector C.  On output the (N+1)-st element of the vector Q'*C.
    !
    !    Output, real(wp) :: C(N), S(N), the cosines and sines of the
    !    transforming Givens rotations.
    !
    implicit none

    integer :: ldr
    integer :: n

    real(wp) :: alpha
    real(wp) :: b(n)
    real(wp) :: c(n)
    real(wp) :: cotan
    integer :: i
    integer :: j
    real(wp) :: r(ldr,n)
    real(wp) :: rowj
    real(wp) :: s(n)
    real(wp) :: tan
    real(wp) :: temp
    real(wp) :: w(n)

    do j = 1, n

        rowj = w(j)
        !
        !  Apply the previous transformations to R(I,J), I=1,2,...,J-1, and to W(J).
        !
        do i = 1, j-1
            temp =   c(i) * r(i,j) + s(i) * rowj
            rowj = - s(i) * r(i,j) + c(i) * rowj
            r(i,j) = temp
        end do
        !
        !  Determine a Givens rotation which eliminates W(J).
        !
        c(j) = 1.0E0_wp
        s(j) = 0.0E0_wp

        if ( rowj /= 0.0E0_wp ) then

            if ( abs( r(j,j) ) < abs( rowj ) ) then
                cotan = r(j,j) / rowj
                s(j) = 0.5E0_wp / sqrt ( 0.25E0_wp + 0.25E0_wp * cotan**2 )
                c(j) = s(j) * cotan
            else
                tan = rowj / r(j,j)
                c(j) = 0.5E0_wp / sqrt ( 0.25E0_wp + 0.25E0_wp * tan**2 )
                s(j) = c(j) * tan
            end if
            !
            !  Apply the current transformation to R(J,J), B(J), and ALPHA.
            !
            r(j,j) =  c(j) * r(j,j) + s(j) * rowj
            temp =    c(j) * b(j)   + s(j) * alpha
            alpha = - s(j) * b(j)   + c(j) * alpha
            b(j) = temp

        end if

    end do

    return
    end
    !-----------------------------------------------------------

    
    subroutine r1mpyq ( m, n, a, lda, v, w )
    !*****************************************************************************80
    !
    !! r1mpyq computes A*Q, where Q is the product of Householder transformations.
    !
    !  Discussion:
    !
    !    Given an M by N matrix A, this subroutine computes A*Q where
    !    Q is the product of 2*(N - 1) transformations
    !
    !      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
    !
    !    and GV(I), GW(I) are Givens rotations in the (I,N) plane which
    !    eliminate elements in the I-th and N-th planes, respectively.
    !    Q itself is not given, rather the information to recover the
    !    GV, GW rotations is supplied.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, integer :: M, the number of rows of A.
    !
    !    Input, integer :: N, the number of columns of A.
    !
    !    Input/output, real(wp) :: A(LDA,N), the M by N array.
    !    On input, the matrix A to be postmultiplied by the orthogonal matrix Q.
    !    On output, the value of A*Q.
    !
    !    Input, integer :: LDA, the leading dimension of A, which must not
    !    be less than M.
    !
    !    Input, real(wp) :: V(N), W(N), contain the information necessary
    !    to recover the Givens rotations GV and GW.
    !
    implicit none

    integer :: lda
    integer :: m
    integer :: n

    real(wp) :: a(lda,n)
    real(wp) :: c
    integer :: i
    integer :: j
    real(wp) :: s
    real(wp) :: temp
    real(wp) :: v(n)
    real(wp) :: w(n)
    !
    !  Apply the first set of Givens rotations to A.
    !
    do j = n-1, 1, -1

        if ( 1.0E0_wp < abs( v(j) ) ) then
            c = 1.0E0_wp / v(j)
            s = sqrt ( 1.0E0_wp - c**2 )
        else
            s = v(j)
            c = sqrt ( 1.0E0_wp - s**2 )
        end if

        do i = 1, m
            temp =   c * a(i,j) - s * a(i,n)
            a(i,n) = s * a(i,j) + c * a(i,n)
            a(i,j) = temp
        end do

    end do
    !
    !  Apply the second set of Givens rotations to A.
    !
    do j = 1, n-1

        if ( abs( w(j) ) > 1.0E0_wp ) then
            c = 1.0E0_wp / w(j)
            s = sqrt ( 1.0E0_wp - c**2 )
        else
            s = w(j)
            c = sqrt ( 1.0E0_wp - s**2 )
        end if

        do i = 1, m
            temp =     c * a(i,j) + s * a(i,n)
            a(i,n) = - s * a(i,j) + c * a(i,n)
            a(i,j) = temp
        end do

    end do

    return
    end
    !-----------------------------------------------------------

    
    subroutine r1updt ( m, n, s, ls, u, v, w, sing )
    !*****************************************************************************80
    !
    !! r1updt re-triangularizes a matrix after a rank one update.
    !
    !  Discussion:
    !
    !    Given an M by N lower trapezoidal matrix S, an M-vector U, and an
    !    N-vector V, the problem is to determine an orthogonal matrix Q such that
    !
    !      (S + U * V' ) * Q
    !
    !    is again lower trapezoidal.
    !
    !    This subroutine determines Q as the product of 2 * (N - 1)
    !    transformations
    !
    !      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
    !
    !    where GV(I), GW(I) are Givens rotations in the (I,N) plane
    !    which eliminate elements in the I-th and N-th planes,
    !    respectively.  Q itself is not accumulated, rather the
    !    information to recover the GV and GW rotations is returned.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow and Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Parameters:
    !
    !    Input, integer :: M, the number of rows of S.
    !
    !    Input, integer :: N, the number of columns of S.  N must not exceed M.
    !
    !    Input/output, real(wp) :: S(LS).  On input, the lower trapezoidal
    !    matrix S stored by columns.  On output S contains the lower trapezoidal
    !    matrix produced as described above.
    !
    !    Input, integer :: LS, the length of the S array.  LS must be at least
    !    (N*(2*M-N+1))/2.
    !
    !    Input, real(wp) :: U(M), the U vector.
    !
    !    Input/output, real(wp) :: V(N).  On input, V must contain the vector V.
    !    On output V contains the information necessary to recover the Givens
    !    rotations GV described above.
    !
    !    Output, real(wp) :: W(M), contains information necessary to
    !    recover the Givens rotations GW described above.
    !
    !    Output, logical :: SING, is set to TRUE if any of the diagonal elements
    !    of the output S are zero.  Otherwise SING is set FALSE.
    !
    implicit none

    integer :: ls
    integer :: m
    integer :: n

    real(wp) :: cos
    real(wp) :: cotan
    real(wp) :: giant
    integer :: i
    integer :: j
    integer :: jj
    integer :: l
    real(wp) :: s(ls)
    real(wp) :: sin
    logical :: sing
    real(wp),parameter :: epsmch = 5.0_wp*epsilon(epsmch)
    real(wp) :: tan
    real(wp) :: tau
    real(wp) :: temp
    real(wp) :: u(m)
    real(wp) :: v(n)
    real(wp) :: w(m)
    !
    !  GIANT is the largest magnitude.
    !
    giant = huge( giant )
    !
    !  Initialize the diagonal element pointer.
    !
    jj = ( n * ( 2 * m - n + 1 ) ) / 2 - ( m - n )
    !
    !  Move the nontrivial part of the last column of S into W.
    !
    l = jj
    do i = n, m
        w(i) = s(l)
        l = l + 1
    end do
    !
    !  Rotate the vector V into a multiple of the N-th unit vector
    !  in such a way that a spike is introduced into W.
    !
    do j = n-1, 1, -1
        jj = jj - ( m - j + 1 )
        w(j) = 0.0E0_wp
        if ( v(j) /= 0.0E0_wp ) then
            !
            !  Determine a Givens rotation which eliminates the
            !  J-th element of V.
            !
            if ( abs( v(n) ) < abs( v(j) ) ) then
                cotan = v(n) / v(j)
                sin = 0.5E0_wp / sqrt ( 0.25E0_wp + 0.25E0_wp * cotan**2 )
                cos = sin * cotan
                tau = 1.0E0_wp
                if ( abs( cos ) * giant > 1.0E0_wp ) then
                    tau = 1.0E0_wp / cos
                end if
            else
                tan = v(j) / v(n)
                cos = 0.5E0_wp / sqrt ( 0.25E0_wp + 0.25E0_wp * tan**2 )
                sin = cos * tan
                tau = sin
            end if
            !
            !  Apply the transformation to V and store the information
            !  necessary to recover the Givens rotation.
            !
            v(n) = sin * v(j) + cos * v(n)
            v(j) = tau
            !
            !  Apply the transformation to S and extend the spike in W.
            !
            l = jj
            do i = j, m
                temp = cos * s(l) - sin * w(i)
                w(i) = sin * s(l) + cos * w(i)
                s(l) = temp
                l = l + 1
            end do
        end if
    end do
    !
    !  Add the spike from the rank 1 update to W.
    !
    w(1:m) = w(1:m) + v(n) * u(1:m)
    !
    !  Eliminate the spike.
    !
    sing = .false.

    do j = 1, n-1
        if ( w(j) /= 0.0E0_wp ) then
            !
            !  Determine a Givens rotation which eliminates the
            !  J-th element of the spike.
            !
            if ( abs( s(jj) ) < abs( w(j) ) ) then

                cotan = s(jj) / w(j)
                sin = 0.5E0_wp / sqrt ( 0.25E0_wp + 0.25E0_wp * cotan**2 )
                cos = sin * cotan

                if ( 1.0E0_wp < abs( cos ) * giant ) then
                    tau = 1.0E0_wp / cos
                else
                    tau = 1.0E0_wp
                end if
            else
                tan = w(j) / s(jj)
                cos = 0.5E0_wp / sqrt ( 0.25E0_wp + 0.25E0_wp * tan**2 )
                sin = cos * tan
                tau = sin
            end if
            !
            !  Apply the transformation to S and reduce the spike in W.
            !
            l = jj
            do i = j, m
                temp = cos * s(l) + sin * w(i)
                w(i) = - sin * s(l) + cos * w(i)
                s(l) = temp
                l = l + 1
            end do
            !
            !  Store the information necessary to recover the Givens rotation.
            !
            w(j) = tau
        end if
        !
        !  Test for zero diagonal elements in the output S.
        !
        if (abs(s(jj)) < epsmch ) then
            sing = .true.
        end if
        jj = jj + ( m - j + 1 )
    end do
    !
    !  Move W back into the last column of the output S.
    !
    l = jj
    do i = n, m
        s(l) = w(i)
        l = l + 1
    end do

    if ( abs(s(jj)) < epsmch) then
        sing = .true.
    end if

    return
    end
    !-----------------------------------------------------------

    
    subroutine timestamp ( )
    !*****************************************************************************80
    !
    !! timestamp prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    May 31 2001   9:45:54.872 AM
    !
    !  Modified:
    !
    !    31 May 2001
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

    character(len = 8) :: ampm
    integer :: d
    character(len = 8) :: date
    integer :: h
    integer :: m
    integer :: mm
    character(len = 9), parameter, dimension(12) :: month = [ &
        'January  ', 'February ', 'March    ', 'April    ', &
        'May      ', 'June     ', 'July     ', 'August   ', &
        'September', 'October  ', 'November ', 'December ' ]
    integer :: n
    integer :: s
    character(len = 10) :: time
    integer :: values(8)
    integer :: y
    character(len = 5) :: zone

    call date_and_time( date, time, zone, values )

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
    write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
        trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

    return
    end
    !-----------------------------------------------------------

end module Mod_MINPACK