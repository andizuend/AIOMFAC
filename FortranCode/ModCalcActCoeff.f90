!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module providing core AIOMFAC subroutines for the calculation of activity coeff.,  *
!*   activities and ion activity products for the present mixture composition and T.    *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2018                                                            *
!*   -> latest changes: 2021-12-02                                                      *
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
!*   -  subroutine AIOMFAC_calc                                                         *
!*   -  subroutine Gammas                                                               *
!*   -  submodule SubModDissociationEquil                                               *
!*                                                                                      *
!****************************************************************************************
module ModCalcActCoeff

use Mod_kind_param, only : wp
use ModAIOMFACvar           !provide access to all composition-dependent AIOMFAC variables

implicit none

real(wp),parameter,private :: lntiny = 0.49_wp*log(tiny(1.0_wp))
real(wp),parameter,private :: lnhuge = 0.49_wp*log(huge(1.0_wp))
real(wp),parameter,private :: logval_threshold = 0.4_wp*lnhuge

!interfaces to procedures in submodules:
interface
    module subroutine HSO4_dissociation()
    end subroutine HSO4_dissociation
    !--
    module subroutine DiffKsulfuricDissoc(mHSO4inp, diffK)
        real(wp),intent(inout) :: mHSO4inp
        real(wp),intent(out) :: diffK
    end subroutine DiffKsulfuricDissoc
    !--
    module function fHSO4dissoc(mHSO4in)
        real(wp), intent(inout) :: mHSO4in
        real(wp) :: fHSO4dissoc    
    end function fHSO4dissoc
    !--
    module subroutine HSO4_and_HCO3_dissociation()
    end subroutine HSO4_and_HCO3_dissociation
    !--
    pure elemental module subroutine rboundsCheck(rset, lowerbound, ntiny, maxlim)
        real(wp),intent(inout) :: rset
        real(wp),intent(in) :: lowerbound, ntiny, maxlim
    end subroutine rboundsCheck
    !--
    pure module function sum_sorted(list) result(summed)
        real(wp),dimension(:),intent(in) :: list
        real(wp) :: summed
    end function sum_sorted
    !-- 
end interface
    
!============================================================================================
    contains
!============================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Determine activity coefficients and activities from given input wtf composition    * 
    !*   for the current system. The LR, MR, & SR activity coeff. contributions are         *
    !*   calculated by the call to 'Gammas' and/or HSO4_dissociation.                   * 
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
    !*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA                *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2004                                                            *
    !*   -> latest changes: 2021-09-23                                                      *
    !*                                                                                      *
    !****************************************************************************************
    subroutine AIOMFAC_calc(WTFin, TKelvin) 
    
    !Public Variables:
    use ModSystemProp, only : anNr, catNr, ElectComps, ElectNues, Ianion, Ication, &
        Nanion, Ncation, nelectrol, nindcomp, nneutral, SolvMixRefnd, bisulfsyst, &
        errorflag_clist, bicarbsyst !, noCO2input
    use ModCompScaleConversion, only : MassFrac2IonMolalities
    use ModNumericalTransformations, only : safe_exp

    implicit none
    !interface variable declarations:
    real(wp),dimension(nindcomp),intent(in) :: WTFin
    real(wp),intent(in) :: TKelvin
    !local variable declarations:
    integer :: I, ii, k, ic, ia
    real(wp),parameter :: dtiny = sqrt(tiny(1.0_wp))
    real(wp) :: t1, ma1, nuelngc, nuelnga, nuec, nuea, trunc
    !.................................................................

    T_K = TKelvin
    wtf = WTFin     !input concentration is in mass fractions; other concentration scales are derived (when needed) from this.
    !if (noCO2input) then   !not needed in the normal call, may be needed with the gas--particle partitioning code
    !    wtf(nneutral+1:nindcomp) = wtf(nneutral:nindcomp-1)
    !    wtf(nneutral) = 0.0_wp
    !endif
    
    !initialize the variables/arrays:
    errorflag_clist = .false.
    gnmrln = 0.0_wp
    gnlrln = 0.0_wp
    gnsrln = 0.0_wp
    galrln = 0.0_wp
    gamrln = 0.0_wp
    gasrln = 0.0_wp
    gclrln = 0.0_wp
    gcmrln = 0.0_wp
    gcsrln = 0.0_wp
    actcoeff_ion = 0.0_wp
    molality_ion = 0.0_wp

    !convert mass fraction input (potentially involving mass frac. of electrolytes) to molalities 
    !of specific ions for use in MR and LR parts:
    call MassFrac2IonMolalities(wtf, SMC, SMA)
    SumIonMolalities = sum(SMA(1:Nanion)) +sum(SMC(1:Ncation))
    
    if (bicarbsyst) then            !includes case for 'bisulfsyst .AND. bicarbsyst'
        call HSO4_and_HCO3_dissociation()
    else if (bisulfsyst) then
        !Perform the dissociation check for bisulfate-system ions and, via this 
        !subroutine, the activity coefficient calculations:
        call HSO4_dissociation()    !many variables are referenced through the module ModAIOMFACvar
    else
        alphaHSO4 = -9.999999_wp
        alphaHCO3 = -9.999999_wp
        call Gammas()               !no partial dissociation of specific ions/electrolytes, so call Gammas for activity 
                                    !coefficient calculation separately.
    endif
    
    !initialize the output arrays:
    lnactcoeff_n = -9999.9_wp
    lnactcoeff_c = -9999.9_wp
    lnactcoeff_a = -9999.9_wp
    activity = 0.0_wp
    actcoeff_n = 0.0_wp
    actcoeff_c = 0.0_wp
    actcoeff_a = 0.0_wp
    actcoeff_ion = 0.0_wp        !to save activity coefficients of ions at their ID number
    molality_ion = 0.0_wp
    
    !ln of the activity coefficient for the neutrals:
    do I = 1,nneutral
        if (wtf(I) > dtiny) then
            lnactcoeff_n(I) = gnmrln(I) +gnsrln(I) +gnlrln(I)
            actcoeff_n(I) = safe_exp(lnactcoeff_n(I), logval_threshold)     !apply safe exponentiation function
        endif
    enddo

    !ln of the activity coefficient for the cations:
    do I = 1,Ncation
        if (SMC(I) > dtiny) then
            lnactcoeff_c(I) = gcmrln(I) +gcsrln(I) +gclrln(I) -Tmolal   !this term converts to the molality scale (basis/scale conversion)
            if (solvmixrefnd) then                                      !correction terms for MR and SR part, if reference solution is the solvent mixture
                lnactcoeff_c(I) = lnactcoeff_c(I) +solvmixcorrMRc(I) +Tmolal -TmolalSolvmix    
            endif
            actcoeff_c(I) = safe_exp(lnactcoeff_c(I), logval_threshold)
        endif
        !save activity coefficient and molality of this ion in an array by actual AIOMFAC ion index:
        ii = Ication(I)
        actcoeff_ion(ii) = actcoeff_c(I)
        molality_ion(ii) = SMC(I)
    enddo
    
    !ln of the activity coefficient for the anions:
    do I = 1,Nanion 
        if (SMA(I) > dtiny) then
            lnactcoeff_a(I) = gamrln(I) +gasrln(I) +galrln(I) -Tmolal
            if (solvmixrefnd) then
                lnactcoeff_a(I) = lnactcoeff_a(I) +solvmixcorrMRa(I) +Tmolal -TmolalSolvmix 
            endif
            actcoeff_a(I) = safe_exp(lnactcoeff_a(I), logval_threshold)
        endif
        !save activity coefficient and molality of this ion in an array by actual AIOMFAC ion index:
        ii = Ianion(I)
        actcoeff_ion(ii) = actcoeff_a(I)
        molality_ion(ii) = SMA(I)
        !note: single-ion activities are:  activity_a(I) = actcoeff_a(I)*SMA(I)  or  activity_ion(ii) = actcoeff_a(I)*SMA(I)
    enddo

    Xwdissoc = x(1)
    activity(1:nneutral) = actcoeff_n(1:nneutral)*x(1:nneutral)

    !if (dicarbsyst) then
    !    pH_calc = -log10(SMC(idH)*actcoeff_c(idH))
    !endif

    !loop over all identified electrolyte components and calculate the corresponding mean molal activity coefficient and molal ion activity product:
    ii = 0
    ionactivityprod = 0.0_wp
    meanmolalactcoeff = 1.0_wp
    lnmeanmactcoeff = 0.0_wp
    do ii = 1,nelectrol
        ic = ElectComps(ii,1)       !get cation identifier 
        ia = ElectComps(ii,2)       !get anion identifier
        i = CatNr(ic)               !array index number of cation ic in, e.g., SMC array
        k = AnNr(ia) 
        if (SMC(i) > 0.0_wp .AND. SMA(k) > 0.0_wp) then
            if (actcoeff_c(i) > 0.0_wp .AND. actcoeff_a(k) > 0.0_wp) then
                nuec = real(ElectNues(ii,1), kind=wp)
                nuea = real(ElectNues(ii,2), kind=wp)
                nuelngc = nuec*lnactcoeff_c(i)
                nuelnga = nuea*lnactcoeff_a(k)
                ma1 = (nuelngc + nuelnga)/(nuec+nuea)
                lnmeanmactcoeff(ii) = ma1
                if (ma1 > lntiny) then                          !no floating point underflow problem expected
                    if (ma1 < lnhuge) then                      !no floating point overflow problem expected
                        meanmolalactcoeff(ii) = exp(ma1)        !the mean molal activity coefficient of the ions of "electrolyte unit" ii
                    else                                        !numerical issue, so output a large number (smaller than overflow risk)
                        trunc = 0.1_wp*log( abs(ma1) )
                        meanmolalactcoeff(ii) = exp(lnhuge + trunc)
                        errorflag_clist(7) = .true.
                    endif
                else !underflow risk
                    meanmolalactcoeff(ii) = exp(lntiny)         !i.e. tiny number, almost zero
                    errorflag_clist(7) = .true.
                endif
                t1 = nuelngc + nuelnga + nuec*log(SMC(i)) + nuea*log(SMA(k))
                if (t1 > lntiny) then                           !no floating point underflow problem expected
                    if (t1 < lnhuge) then                       !no floating point overflow problem expected
                        ionactivityprod(ii) = exp(t1)           !the mean molal activity coefficient of the ions of "electrolyte unit" ii
                    else                                        !numerical issue, so output a large number (smaller than overflow risk)
                        trunc = 0.1_wp*log( abs(t1) )
                        ionactivityprod(ii) = exp(lnhuge + trunc)
                        errorflag_clist(7) = .true.
                    endif
                else                                            !underflow risk
                    trunc = 0.1_wp*log( abs(t1) )
                    ionactivityprod(ii) = exp(lntiny - trunc)   !i.e. tiny number, almost zero
                    errorflag_clist(7) = .true.
                endif
            endif
        endif
    enddo !ii
    activity(nneutral+1:nindcomp) = ionactivityprod(1:nelectrol)

    end subroutine AIOMFAC_calc 
    !==========================================================================================
    

    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to calculate the activity coefficient values ("gammas") by calling the  *
    !*   LR, MR and SR part routines with given mixture composition and temperature data.   *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich (2004 - 2009),                                                  *
    !*   Dept. Chem. Engineering, California Institute of Technology (2009 - 2012),         *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2004                                                            *
    !*   -> latest changes: 2023-03-17                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    subroutine Gammas()
    
    use ModSystemProp, only : nneutral, NKNpNGS, Ncation, Nanion, Mmass, bisulfsyst, bicarbsyst
    use ModSRunifac, only : SRsystm, SRunifac
    use ModMRpart, only : LR_MR_activity, GammaCO2
    use ModAIOMFACvar, only : meanSolventMW
    use ModCompScaleConversion, only : MassFrac2SolvMolalities

    implicit none
    !local variables:
    integer :: NKNpNcat, nnp1
    real(wp) :: sum_molalities, sumXN
    real(wp),dimension(NKNpNGS) :: lnGaSR
    real(wp),dimension(nneutral) :: wtfbycompMW, mNeutral
    logical :: refreshgref
    !........................................................................................................
    
    nnp1 = nneutral +1
    !initialize the mole fraction composition arrays:
    X = 0.0_wp
    XN = 0.0_wp
    !Calculation of the electrolyte-free mole fractions of the neutral components, XN:
    wtfbycompMW(1:nneutral) = wtf(1:nneutral)/Mmass(1:nneutral)
    sumXN = sum(wtfbycompMW)
    XN(1:nneutral) = wtfbycompMW/sumXN
    
    meanSolventMW = sum(Mmass(1:nneutral)*XN(1:nneutral))               !mean molar mass of the solvent mixture (only non-electrolytes)
    mNeutral = XN(1:nneutral)/meanSolventMW                             ![mol/kg], the molalities of the neutrals

    !Addition of the moles of substance (neutral and ionic) per 1 kg of electrolyte-free solvent mixture
    if (bisulfsyst .OR. bicarbsyst) then    !update since it changes in DiffKsulfuricDissoc and DiffKcarbonateDissoc
        SumIonMolalities = sum(SMA(1:Nanion)) +sum(SMC(1:Ncation))
    endif
    sum_molalities = sum(mNeutral) + SumIonMolalities                   !sum of all molalities

    !Calculation of the mole fraction (X) of the neutral components and the ions with respect to dissociated electrolytes/ions.
    !==> the structure of the mole fraction array X is: 
    !1) neutral components in component order,
    !2) ions: first the cations, then the anions
    X(1:nneutral) = mNeutral/sum_molalities                             !mole fraction of the neutral components (on the basis of dissociated electrolytes)    
    NKNpNcat = nneutral + Ncation
    X(nnp1:NKNpNcat) = SMC(1:Ncation)/sum_molalities                    !mole fractions of the cations
    X(NKNpNcat+1:NKNpNcat+Nanion) = SMA(1:Nanion)/sum_molalities        !mole fractios of the anions

    !check whether temperature-dependent parameters need to be updated in SR and LR parts:
    if (abs(T_K - lastTK) > 1.0E-2_wp) then !detected a change in temperature --> PsiT and other coeff. need to be updated
        refreshgref = .true.
        DebyeHrefresh = .true.
        lastTK = T_K
    else
        refreshgref = .false.
        DebyeHrefresh = .false.
    endif
    !calculate the LR and MR activity coefficient contributions (gammas):  
    call LR_MR_activity()
    
    !call the UNIFAC model part for the short-range interaction contributions.
    call SRunifac(NKNpNGS, T_K, X, XN, refreshgref, lnGaSR) 
    gnsrln(1:nneutral) = lnGaSR(1:nneutral)
    gcsrln(1:Ncation)  = lnGaSR(nnp1:NKNpNcat)
    gasrln(1:Nanion)   = lnGaSR(NKNpNcat+1:NKNpNcat+Nanion)

    !update activity coefficient of CO2(aq):
    if (bicarbsyst) then
        call GammaCO2()
    endif

    end subroutine Gammas 
    !==========================================================================================
    
end module ModCalcActCoeff