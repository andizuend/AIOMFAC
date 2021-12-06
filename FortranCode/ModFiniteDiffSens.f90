!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module for finite difference computations of AIOMFAC activities, activity coeff.,  *
!*   and mixture viscosity at given temperature and composition. This serves the        *
!*   estimation of model sensitivities with respect to small changes in molar amounts.  *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2018                                                            *
!*   -> latest changes: 2021/11/24                                                      *
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
!*   -  SUBROUTINE DeltaActivities                                                      *
!*   -  SUBROUTINE partialdactcoeff                                                     *
!*                                                                                      *
!****************************************************************************************
MODULE ModFiniteDiffSens

IMPLICIT NONE

!public procedures
PUBLIC :: DeltaActivities, partialdactcoeff
PRIVATE

!============================================================================================
    CONTAINS
!============================================================================================
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to calculate total derivatives (finite differences) of the AIOMFAC      *
    !*   activities and activity coefficients at a given composition (xin) of mixture nd.   *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   Dept. Chem. Engineering, California Institute of Technology (2009 - 2012),         *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2010                                                            *
    !*   -> latest changes: 2021-11-24                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    SUBROUTINE DeltaActivities(xin, TKelvin, onlyDeltaVisc, dact, dactcoeff)

    USE ModSystemProp, ONLY : nindcomp
    USE ModAIOMFACvar, ONLY : partial_log10_etamix

    IMPLICIT NONE

    !interface variables:
    REAL(8),DIMENSION(nindcomp),INTENT(IN) :: xin
    REAL(8),INTENT(IN) :: TKelvin
    LOGICAL(4),INTENT(IN) :: onlyDeltaVisc                  !this is set .true. when only the first component is varied; use only for delta(etamix) calculation.
    REAL(8),DIMENSION(nindcomp),INTENT(OUT) :: dact, dactcoeff
    !local variables:
    INTEGER(4) :: ni, j
    REAL(8),DIMENSION(nindcomp) :: xinit
    REAL(8),DIMENSION(:,:),ALLOCATABLE :: partdact, partdactcoeff
    REAL(8) :: partdact_ji, partdactcoeff_ji
    LOGICAL(4),PARAMETER :: xinputtype = .true.             !set to mole fraction as input concentration
    !...........................................................
    !parameters:
    xinit = xin     !xinit is locally stored and unaffected by changes when reentrant code will be called that might feed back on xin!

    !calculate the finite differences (numerical "partial derivatives") of the component ni's activity at composition xinit while holding the moles of the other componentes fixed.
    !compositions for the forward / backward differences at the given point xinit (component nnvar is the one to be enhanced/diminished):
    IF (onlyDeltaVisc) THEN !e.g. for web-model version, only perform calculation with variation in component 1 (often water) and effect on component 1 only for partial_log10_etamix etc. (much faster).
        ALLOCATE( partdact(1,1), partdactcoeff(1,1) )
        partdact = 0.0D0
        partdactcoeff = 0.0D0
        CALL partialdactcoeff(xinit, TKelvin, 1, 1, partdact_ji, partdactcoeff_ji)
        partdact(1,1) = partdact_ji
        partdactcoeff(1,1) = partdactcoeff_ji    
    ELSE
        ALLOCATE( partdact(nindcomp,nindcomp), partdactcoeff(nindcomp,nindcomp) )
        partdact = 0.0D0
        partdactcoeff = 0.0D0
        DO ni = 1,nindcomp  !loop
            !calling the subroutine calculating the partial derivative of the activity and act. coeff. of component j with respect to moles of component ni:
            DO j = 1,nindcomp
                CALL partialdactcoeff(xinit, TKelvin, j, ni, partdact_ji, partdactcoeff_ji)
                partdact(j,ni) = partdact_ji
                partdactcoeff(j,ni) = partdactcoeff_ji
            ENDDO
        ENDDO !ni
    ENDIF
    !---
    !calculate the total derivative for each component from the absolute values of the partial derivatives:
    dact = 0.0D0
    dactcoeff = 0.0D0
    IF (.NOT. onlyDeltaVisc) THEN
        DO ni = 1,nindcomp
            dact(ni) = SUM(ABS(partdact(ni,1:nindcomp)))
            dactcoeff(ni) = SUM(ABS(partdactcoeff(ni,1:nindcomp)))
        ENDDO
    ENDIF
    !for viscosity sensitivity estimate as a function of a change in mole fraction of component 1 (water):
    partial_log10_etamix = ABS(partial_log10_etamix)
    
    DEALLOCATE(partdact, partdactcoeff)

    END SUBROUTINE DeltaActivities
    !==========================================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to calculate partial derivatives (finite differences) of the AIOMFAC    *
    !*   activities and activity coefficients of a certain mixture component "j" with       *
    !*   respect to a tiny change in moles of component "i" at a composition vector (xin).  *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   Dept. Chem. Engineering, California Institute of Technology (2009 - 2012),         *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2010                                                            *
    !*   -> latest changes: 2021-11-24                                                      *
    !*                                                                                      *
    !****************************************************************************************
    SUBROUTINE partialdactcoeff(xin, TKelvin, j, i, partdact_ji, partdactcoeff_ji)

    USE ModSystemProp, ONLY : nindcomp, nneutral, nd, Mmass, errorflagcalc, calcviscosity
    USE ModCompScaleConversion
    USE ModAIOMFACvar, ONLY : activity, meanmolalactcoeff, actcoeff_n, partial_log10_etamix, ln_etamix
    USE ModCalcActCoeff, ONLY : AIOMFAC_calc
    
    IMPLICIT NONE

    !interface variables:
    INTEGER(4),INTENT(IN) :: j, i
    REAL(8),DIMENSION(nindcomp),INTENT(IN) :: xin
    REAL(8),INTENT(IN) :: TKelvin
    REAL(8),INTENT(OUT) :: partdact_ji, partdactcoeff_ji
    !local variables:
    REAL(8),DIMENSION(nindcomp) :: xinit, xplus, wtf
    REAL(8) :: ntplus, actplus, actinit, actcoeffplus, actcoeffinit, ln_etamixplus, ln_etamixinit
    REAL(8),PARAMETER :: Dtiny = EPSILON(1.0D0)
    REAL(8),PARAMETER :: dn = SQRT(EPSILON(1.0D0))              ![mol] a small molar (change) for the numerical differentiation
    REAL(8),PARAMETER :: ln10 = LOG(10.0D0)
    LOGICAL(4),PARAMETER :: xinputtype = .true.                 !set mole fraction as input concentration
    !...........................................................
    
    !parameters:
    xinit = xin     !xinit is locally stored and unaffected by changes when reentrant code will be called that might feed back on xin!
    
    !calculate the finite differences (numerical "partial derivatives") of component j's activity and activity coeff. with respect to component i,
    !at composition xinit, while holding the moles of the other components fixed.
    !compositions for the forward / backward differences at the given point xinit (component i is the one to be enhanced/diminished):
    ntplus = 1.0D0 + dn     !as the molar sum of all components can be set arbitrarily to 1.0
    !calculate the mole fractions at the slightly changed compositions:
    xplus(1:i-1) = xinit(1:i-1)/ntplus
    xplus(i) = (xinit(i)+dn)/ntplus
    xplus(i+1:nindcomp) = xinit(i+1:nindcomp)/ntplus
    !..
    !Call AIOMFAC_calc to calculate the activity at composition xplus:
    CALL MoleFrac2MassFrac(xplus, Mmass, wtf) 
    IF (j == 1 .AND. i == 1) THEN
        CALL AIOMFAC_calc(wtf, TKelvin)
        ln_etamixplus = ln_etamix
    ELSE
        CALL AIOMFAC_calc(wtf, TKelvin)
    ENDIF
    actplus = activity(j)
    IF (j <= nneutral) THEN !neutral component
        actcoeffplus = actcoeff_n(j)
    ELSE !electrolyte component
        actcoeffplus = meanmolalactcoeff(j-nneutral)
    ENDIF
    !..
    !calculate the activity & activity coeff. at the initial composition xinit:
    CALL MoleFrac2MassFrac(xinit, Mmass, wtf) 
    IF (j == 1 .AND. i == 1) THEN
        CALL AIOMFAC_calc(wtf, TKelvin)
        ln_etamixinit = ln_etamix
        !partial forward difference for mixture viscosity:
        IF (calcviscosity) THEN
            partial_log10_etamix = (ln_etamixplus - ln_etamixinit) / (dn * ln10)
        ELSE
            partial_log10_etamix = 0.0D0
        ENDIF
    ELSE
        CALL AIOMFAC_calc(wtf, TKelvin)
    ENDIF
    actinit = activity(j)
    IF (j <= nneutral) THEN
        actcoeffinit = actcoeff_n(j)
    ELSE
        actcoeffinit = meanmolalactcoeff(j-nneutral)
    ENDIF
    !..
    !calculate the numerical forward differences with respect to activity at point xinit 
    !(partial derivative of component j with resp. to a small molar change, dn, in component i)
    IF (actplus > 0.0D0 .AND. actinit > 0.0D0) THEN
        IF (ABS(actplus -actinit) < 1.0D-292) THEN                  !floating point underflow risk
            IF (actplus -actinit < 0.0D0) THEN                      !"negative zero"
                partdact_ji = -1.0D-292                             !almost zero
            ELSE                                                    !"positive zero"
                partdact_ji = 1.0D-292                  
            ENDIF
        ELSE !no underflow
            IF (actplus -actinit > 1.0D292) THEN                    !floating overflow risk
                partdact_ji = 1.0D292 +(LOG(actplus) -LOG(actinit)) !set to huge positive number to allow following computations but prevent overflow
            ELSE IF (actplus -actinit < -1.0D292) THEN              !floating overflow risk
                partdact_ji = -1.0D292 -(LOG(actplus) -LOG(actinit))!set to huge negative number to allow following computations but prevent overflow
            ELSE
                partdact_ji = (actplus -actinit)/dn                 !the finite difference derivative with respect to component i (while all other component's moles are kept const.) 
            ENDIF
        ENDIF
        IF (ABS(actcoeffplus -actcoeffinit) < 1.0D-292) THEN        !floating underflow risk
            IF (actcoeffplus -actcoeffinit < 0.0D0) THEN 
                partdactcoeff_ji = -1.0D-292 
            ELSE !"positive zero"
                partdactcoeff_ji = 1.0D-292 
            ENDIF
        ELSE !ok, no underflow
            IF (actcoeffplus -actcoeffinit > 1.0D292) THEN 
                partdactcoeff_ji = 1.0D292 +(LOG(actcoeffplus) -LOG(actcoeffinit)) 
            ELSE IF (actcoeffplus -actcoeffinit < -1.0D292) THEN
                partdactcoeff_ji = -1.0D292 -(LOG(actcoeffplus) -LOG(actcoeffinit))
            ELSE
                partdactcoeff_ji = (actcoeffplus -actcoeffinit)/dn
            ENDIF
        ENDIF
    ELSE IF (xinit(j) > Dtiny .AND. (actplus < Dtiny .OR. actinit < Dtiny)) THEN 
        !there is some amount of a component in the mixture, but the model activity coefficient is too large (thus, very steep derivative)
        partdact_ji = 1.11111111D5  !value indicating a floating point overflow in AIOMFAC_calc, while not suppressing the feedback of such a value.
        partdactcoeff_ji = 1.11111111D5
        errorflagcalc = 2
    ELSE !exceptions where problem occurred
        partdact_ji = 0.0D0
        partdactcoeff_ji = 0.0D0
        errorflagcalc = 3
    ENDIF
   
    END SUBROUTINE partialdactcoeff
    !========================================================================================================================== 
    
END MODULE ModFiniteDiffSens