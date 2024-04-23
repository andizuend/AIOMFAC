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
!*   -> latest changes: 2021-11-24                                                      *
!*                                                                                      *
!*   :: License ::                                                                      *
!*   This program is free software: you can redistribute it and/or modify it under the  *
!*   terms of the GNU General Public License as published by the Free Software          *
!*   Foundation, either version 3 of the License, or (at your option) any later         *
!*   version.                                                                           *
!*   The AIOMFAC model code is distributed in the hope that it will be useful, but      *
!*   WITHOUT any WARRANTY; without even the implied warranty of MERCHANTABILITY or      *
!*   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      *
!*   details.                                                                           *
!*   You should have received a copy of the GNU General Public License along with this  *
!*   program. If not, see <http://www.gnu.org/licenses/>.                               *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  subroutine DeltaActivities                                                      *
!*   -  subroutine partialdactcoeff                                                     *
!*                                                                                      *
!****************************************************************************************
module ModFiniteDiffSens

use Mod_kind_param, only : wp

implicit none

!public procedures
public :: DeltaActivities, partialdactcoeff
private

!============================================================================================
    contains
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
    subroutine DeltaActivities(xin, TKelvin, onlyDeltaVisc, dact, dactcoeff)

    use ModSystemProp, only : nindcomp
    use ModAIOMFACvar, only : partial_log10_etamix

    implicit none

    !interface variables:
    real(wp),dimension(nindcomp),intent(in) :: xin
    real(wp),intent(in) :: TKelvin
    logical,intent(in) :: onlyDeltaVisc                  !this is set .true. when only the first component is varied; use only for delta(etamix) calculation.
    real(wp),dimension(nindcomp),intent(out) :: dact, dactcoeff
    !local variables:
    integer :: ni, j
    real(wp),dimension(nindcomp) :: xinit
    real(wp),dimension(:,:),allocatable :: partdact, partdactcoeff
    real(wp) :: partdact_ji, partdactcoeff_ji
    logical,parameter :: xinputtype = .true.             !set to mole fraction as input concentration
    !...........................................................
    !parameters:
    xinit = xin     !xinit is locally stored and unaffected by changes when reentrant code will be called that might feed back on xin!

    !calculate the finite differences (numerical "partial derivatives") of the component ni's activity at composition xinit while holding the moles of the other componentes fixed.
    !compositions for the forward / backward differences at the given point xinit (component nnvar is the one to be enhanced/diminished):
    if (onlyDeltaVisc) then !e.g. for web-model version, only perform calculation with variation in component 1 (often water) and effect on component 1 only for partial_log10_etamix etc. (much faster).
        allocate( partdact(1,1), partdactcoeff(1,1) )
        partdact = 0.0_wp
        partdactcoeff = 0.0_wp
        call partialdactcoeff(xinit, TKelvin, 1, 1, partdact_ji, partdactcoeff_ji)
        partdact(1,1) = partdact_ji
        partdactcoeff(1,1) = partdactcoeff_ji    
    else
        allocate( partdact(nindcomp,nindcomp), partdactcoeff(nindcomp,nindcomp) )
        partdact = 0.0_wp
        partdactcoeff = 0.0_wp
        do ni = 1,nindcomp  !loop
            !calling the subroutine calculating the partial derivative of the activity and act. coeff. of component j with respect to moles of component ni:
            do j = 1,nindcomp
                call partialdactcoeff(xinit, TKelvin, j, ni, partdact_ji, partdactcoeff_ji)
                partdact(j,ni) = partdact_ji
                partdactcoeff(j,ni) = partdactcoeff_ji
            enddo
        enddo !ni
    endif
    !---
    !calculate the total derivative for each component from the absolute values of the partial derivatives:
    dact = 0.0_wp
    dactcoeff = 0.0_wp
    if (.NOT. onlyDeltaVisc) then
        do ni = 1,nindcomp
            dact(ni) = sum(abs(partdact(ni,1:nindcomp)))
            dactcoeff(ni) = sum(abs(partdactcoeff(ni,1:nindcomp)))
        enddo
    endif
    !for viscosity sensitivity estimate as a function of a change in mole fraction of component 1 (water):
    partial_log10_etamix = abs(partial_log10_etamix)
    
    deallocate(partdact, partdactcoeff)

    end subroutine DeltaActivities
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
    subroutine partialdactcoeff(xin, TKelvin, j, i, partdact_ji, partdactcoeff_ji)

    use ModSystemProp, only : nindcomp, nneutral, nd, Mmass, errorflag_clist, calcviscosity
    use ModCompScaleConversion
    use ModAIOMFACvar, only : activity, meanmolalactcoeff, actcoeff_n, partial_log10_etamix, ln_etamix
    use ModCalcActCoeff, only : AIOMFAC_calc
    
    implicit none

    !interface variables:
    integer,intent(in) :: j, i
    real(wp),dimension(nindcomp),intent(in) :: xin
    real(wp),intent(in) :: TKelvin
    real(wp),intent(out) :: partdact_ji, partdactcoeff_ji
    !local variables:
    real(wp),dimension(nindcomp) :: xinit, xplus, wtf
    real(wp) :: ntplus, actplus, actinit, actcoeffplus, actcoeffinit, ln_etamixplus, ln_etamixinit
    real(wp),parameter :: deps = epsilon(1.0_wp), tinyn = tiny(1.0_wp), hugen = 1.0_wp/tinyn
    real(wp),parameter :: dn = sqrt(epsilon(1.0_wp))              ![mol] a small molar (change) for the numerical differentiation
    real(wp),parameter :: ln10 = log(10.0_wp)
    logical,parameter :: xinputtype = .true.                 !set mole fraction as input concentration
    !...........................................................
    
    !parameters:
    xinit = xin     !xinit is locally stored and unaffected by changes when reentrant code will be called that might feed back on xin!
    
    !calculate the finite differences (numerical "partial derivatives") of component j's activity and activity coeff. with respect to component i,
    !at composition xinit, while holding the moles of the other components fixed.
    !compositions for the forward / backward differences at the given point xinit (component i is the one to be enhanced/diminished):
    ntplus = 1.0_wp + dn     !as the molar sum of all components can be set arbitrarily to 1.0_wp
    !calculate the mole fractions at the slightly changed compositions:
    xplus(1:i-1) = xinit(1:i-1)/ntplus
    xplus(i) = (xinit(i)+dn)/ntplus
    xplus(i+1:nindcomp) = xinit(i+1:nindcomp)/ntplus
    !..
    !Call AIOMFAC_calc to calculate the activity at composition xplus:
    call MoleFrac2MassFrac(xplus, Mmass, wtf) 
    if (j == 1 .AND. i == 1) then
        call AIOMFAC_calc(wtf, TKelvin)
        ln_etamixplus = ln_etamix
    else
        call AIOMFAC_calc(wtf, TKelvin)
    endif
    actplus = activity(j)
    if (j <= nneutral) then !neutral component
        actcoeffplus = actcoeff_n(j)
    else !electrolyte component
        actcoeffplus = meanmolalactcoeff(j-nneutral)
    endif
    !..
    !calculate the activity & activity coeff. at the initial composition xinit:
    call MoleFrac2MassFrac(xinit, Mmass, wtf) 
    if (j == 1 .AND. i == 1) then
        call AIOMFAC_calc(wtf, TKelvin)
        ln_etamixinit = ln_etamix
        !partial forward difference for mixture viscosity:
        if (calcviscosity) then
            partial_log10_etamix = (ln_etamixplus - ln_etamixinit) / (dn * ln10)
        else
            partial_log10_etamix = 0.0_wp
        endif
    else
        call AIOMFAC_calc(wtf, TKelvin)
    endif
    actinit = activity(j)
    if (j <= nneutral) then
        actcoeffinit = actcoeff_n(j)
    else
        actcoeffinit = meanmolalactcoeff(j-nneutral)
    endif
    !..
    !calculate the numerical forward differences with respect to activity at point xinit 
    !(partial derivative of component j with resp. to a small molar change, dn, in component i)
    if (actplus > 0.0_wp .AND. actinit > 0.0_wp) then
        if (abs(actplus -actinit) < tinyn) then                  !floating point underflow risk
            if (actplus -actinit < 0.0_wp) then                      !"negative zero"
                partdact_ji = -tinyn                             !almost zero
            else                                                    !"positive zero"
                partdact_ji = tinyn                  
            endif
        else !no underflow
            if (actplus -actinit > hugen) then                    !floating overflow risk
                partdact_ji = hugen +(log(actplus) -log(actinit)) !set to huge positive number to allow following computations but prevent overflow
            else if (actplus -actinit < -hugen) then              !floating overflow risk
                partdact_ji = -hugen -(log(actplus) -log(actinit))!set to huge negative number to allow following computations but prevent overflow
            else
                partdact_ji = (actplus -actinit)/dn                 !the finite difference derivative with respect to component i (while all other component's moles are kept const.) 
            endif
        endif
        if (abs(actcoeffplus -actcoeffinit) < tinyn) then        !floating underflow risk
            if (actcoeffplus -actcoeffinit < 0.0_wp) then 
                partdactcoeff_ji = -tinyn 
            else !"positive zero"
                partdactcoeff_ji = tinyn 
            endif
        else !ok, no underflow
            if (actcoeffplus -actcoeffinit > hugen) then 
                partdactcoeff_ji = hugen +(log(actcoeffplus) -log(actcoeffinit)) 
            else if (actcoeffplus -actcoeffinit < -hugen) then
                partdactcoeff_ji = -hugen -(log(actcoeffplus) -log(actcoeffinit))
            else
                partdactcoeff_ji = (actcoeffplus -actcoeffinit)/dn
            endif
        endif
    else if (xinit(j) > deps .AND. (actplus < deps .OR. actinit < deps)) then 
        !there is some amount of a component in the mixture, but the model activity coefficient is too large (thus, very steep derivative)
        partdact_ji = 1.11111111E5_wp  !value indicating a floating point overflow in AIOMFAC_calc, while not suppressing the feedback of such a value.
        partdactcoeff_ji = 1.11111111E5_wp
        errorflag_clist(2) = .true.
    else !exceptions where problem occurred
        partdact_ji = 0.0_wp
        partdactcoeff_ji = 0.0_wp
        errorflag_clist(3) = .true.
    endif
   
    end subroutine partialdactcoeff
    !========================================================================================================================== 
    
end module ModFiniteDiffSens