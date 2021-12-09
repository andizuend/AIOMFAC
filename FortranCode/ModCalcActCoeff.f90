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
!*   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or      *
!*   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      *
!*   details.                                                                           *
!*   You should have received a copy of the GNU General Public License along with this  *
!*   program. If not, see <http://www.gnu.org/licenses/>.                               *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  SUBROUTINE AIOMFAC_calc                                                         *
!*   -  SUBROUTINE Gammas                                                               *
!*   -  SUBMODULE SubModDissociationEquil                                               *
!*                                                                                      *
!****************************************************************************************
MODULE ModCalcActCoeff

USE ModAIOMFACvar           !provide access to all composition-dependent AIOMFAC variables

IMPLICIT NONE

REAL(8),PARAMETER,PRIVATE :: lntiny = 0.49D0*LOG(TINY(1.0D0))
REAL(8),PARAMETER,PRIVATE :: lnhuge = 0.49D0*LOG(HUGE(1.0D0))
REAL(8),PARAMETER,PRIVATE :: logval_threshold = 0.4D0*lnhuge

!interfaces to procedures in submodules:
INTERFACE
    MODULE SUBROUTINE HSO4_dissociation()
    END SUBROUTINE HSO4_dissociation
    !--
    MODULE SUBROUTINE DiffKsulfuricDissoc(mHSO4inp, diffK)
        REAL(8),INTENT(INOUT) :: mHSO4inp
        REAL(8),INTENT(OUT) :: diffK
    END SUBROUTINE DiffKsulfuricDissoc
    !--
    MODULE FUNCTION fHSO4dissoc(mHSO4in)
        REAL(8), INTENT(INOUT) :: mHSO4in
        REAL(8) :: fHSO4dissoc    
    END FUNCTION fHSO4dissoc
    !--
    MODULE SUBROUTINE HSO4_and_HCO3_dissociation()
    END SUBROUTINE HSO4_and_HCO3_dissociation
    !--
    PURE ELEMENTAL MODULE SUBROUTINE rboundsCheck(rset, lowerbound, ntiny, maxlim)
        REAL(8),INTENT(INOUT) :: rset
        REAL(8),INTENT(IN) :: lowerbound, ntiny, maxlim
    END SUBROUTINE rboundsCheck
    !--
    PURE MODULE FUNCTION sum_sorted(list) RESULT(summed)
        REAL(8),DIMENSION(:),INTENT(IN) :: list
        REAL(8) :: summed
    END FUNCTION sum_sorted
    !--
END INTERFACE
    
!============================================================================================
    CONTAINS
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
    SUBROUTINE AIOMFAC_calc(WTFin, TKelvin) 
    
    !Public Variables:
    USE ModSystemProp, ONLY : anNr, catNr, ElectComps, ElectNues, Ianion, Ication, &
        Nanion, Ncation, nelectrol, nindcomp, nneutral, SolvMixRefnd, bisulfsyst, &
        errorflagcalc, bicarbsyst, noCO2input
    USE ModCompScaleConversion, ONLY : MassFrac2IonMolalities

    IMPLICIT NONE
    !interface variables declarations
    REAL(8),DIMENSION(nindcomp),INTENT(IN) :: WTFin
    REAL(8),INTENT(IN) :: TKelvin
    !local variable declarations:
    INTEGER(4) :: I, ii, k, ic, ia
    REAL(8) :: t1, ma1, nuelngc, nuelnga, nuec, nuea, trunc
    REAL(8),PARAMETER :: dtiny = SQRT(TINY(1.0D0))
    !.................................................................

    T_K = TKelvin
    wtf = WTFin     !input concentration is in mass fractions; other concentration scales are derived (when needed) from this.
    !IF (noCO2input) THEN   !not needed in the normal call, may be needed with the gas--particle partitioning code
    !    wtf(nneutral+1:nindcomp) = wtf(nneutral:nindcomp-1)
    !    wtf(nneutral) = 0.0D0
    !ENDIF
    
    !initialize the variables/arrays:
    errorflagcalc = 0
    gnmrln = 0.0D0
    gnlrln = 0.0D0
    gnsrln = 0.0D0
    galrln = 0.0D0
    gamrln = 0.0D0
    gasrln = 0.0D0
    gclrln = 0.0D0
    gcmrln = 0.0D0
    gcsrln = 0.0D0
    actcoeff_ion = 0.0D0
    molality_ion = 0.0D0

    !convert mass fraction input to ion molalities for use in MR and LR parts:
    CALL MassFrac2IonMolalities(wtf, SMC, SMA)
    SumIonMolalities = SUM(SMA(1:Nanion)) +SUM(SMC(1:Ncation))
    
    IF (bicarbsyst) THEN            !includes case for 'bisulfsyst .AND. bicarbsyst'
        CALL HSO4_and_HCO3_dissociation()
    ELSE IF (bisulfsyst) THEN
        !Perform the dissociation check for bisulfate-system ions and, via this subroutine, the activity coefficient calculations:
        CALL HSO4_dissociation()    !many variables are referenced through the module ModAIOMFACvar
    ELSE
        alphaHSO4 = -9.999999D0
        alphaHCO3 = -9.999999D0
        CALL Gammas()   !No dissociation of specific ions/electrolytes, so call Gammas for activity coefficient calculation separately.
    ENDIF
    
    !initialize the output arrays:
    lnactcoeff_n = -9999.9D0
    lnactcoeff_c = -9999.9D0
    lnactcoeff_a = -9999.9D0
    activity = 0.0D0
    actcoeff_n = 0.0D0
    actcoeff_c = 0.0D0
    actcoeff_a = 0.0D0
    actcoeff_ion = 0.0D0 !to save activity coefficients of ions at their ID number
    molality_ion = 0.0D0
    
    !ln of the activity coefficient for the neutrals:
    DO I = 1,nneutral
        IF (wtf(I) > dtiny) THEN
            lnactcoeff_n(I) = gnmrln(I) +gnsrln(I) +gnlrln(I)
            IF (lnactcoeff_n(I) < -logval_threshold .OR. lnactcoeff_n(I) > logval_threshold) THEN 
                trunc = 1.0D-2*LOG( ABS(lnactcoeff_n(I)) )     !reduce magnitude via log and scaling
                IF (lnactcoeff_n(I) < -logval_threshold) THEN
                    actcoeff_n(I) = EXP(-logval_threshold - trunc)     
                ELSE
                    actcoeff_n(I) = EXP(logval_threshold + trunc)
                ENDIF
            ELSE
                actcoeff_n(I) = EXP(lnactcoeff_n(I))
            ENDIF
        ENDIF
    ENDDO

    !ln of the activity coefficient for the cations:
    DO I = 1,Ncation
        IF (SMC(I) > dtiny) THEN
            lnactcoeff_c(I) = gcmrln(I) +gcsrln(I) +gclrln(I) -Tmolal  !this term converts to the molality scale (basis/scale conversion)
            IF (solvmixrefnd) THEN !correction terms for MR and SR part, because reference solution is the solvent mixture
                lnactcoeff_c(I) = lnactcoeff_c(I) +solvmixcorrMRc(I) +Tmolal -TmolalSolvmix    
            ENDIF
            IF (lnactcoeff_c(I) < -logval_threshold .OR. lnactcoeff_c(I) > logval_threshold) THEN 
                trunc = 1.0D-2*LOG( ABS(lnactcoeff_c(I)) )
                IF (lnactcoeff_c(I) < -logval_threshold) THEN
                    actcoeff_c(I) = EXP(-logval_threshold - trunc)    
                ELSE
                    actcoeff_c(I) = EXP(logval_threshold + trunc)
                ENDIF
            ELSE
                actcoeff_c(I) = EXP(lnactcoeff_c(I))
            ENDIF
        ENDIF
        !save activity coefficient and molality of this ion in an array by actual AIOMFAC ion index:
        ii = Ication(I)
        actcoeff_ion(ii) = actcoeff_c(I)
        molality_ion(ii) = SMC(I)
    ENDDO
    
    !ln of the activity coefficient for the anions:
    DO I = 1,Nanion 
        IF (SMA(I) > dtiny) THEN
            lnactcoeff_a(I) = gamrln(I) +gasrln(I) +galrln(I) -Tmolal    !this term converts to the molality scale (basis/scale conversion)
            IF (solvmixrefnd) THEN  !correction terms for MR and SR because reference solution is the solvent mixture
                lnactcoeff_a(I) = lnactcoeff_a(I) +solvmixcorrMRa(I) +Tmolal -TmolalSolvmix 
            ENDIF
            IF (lnactcoeff_a(I) < -logval_threshold .OR. lnactcoeff_a(I) > logval_threshold) THEN 
                trunc = 1.0D-2*LOG( ABS(lnactcoeff_a(I)) )
                IF (lnactcoeff_a(I) < -logval_threshold) THEN
                    actcoeff_a(I) = EXP(-logval_threshold - trunc)  
                ELSE
                    actcoeff_a(I) = EXP(logval_threshold + trunc) 
                ENDIF
            ELSE
                actcoeff_a(I) = EXP(lnactcoeff_a(I))
            ENDIF
        ENDIF
        !save activity coefficient and molality of this ion in an array by actual AIOMFAC ion index:
        ii = Ianion(I)
        actcoeff_ion(ii) = actcoeff_a(I)
        molality_ion(ii) = SMA(I)
        !note: single-ion activities are:  activity_a(I) = actcoeff_a(I)*SMA(I)  or  activity_ion(ii) = actcoeff_a(I)*SMA(I)
    ENDDO

    Xwdissoc = x(1)
    activity(1:nneutral) = actcoeff_n(1:nneutral)*x(1:nneutral)

    !!!notify exception in case of floating point overflow problems and return:
    !!IF (floatingproblem) THEN
    !!    errorflagcalc = 6
    !!    activity(1:nindcomp) = -9999.9D0
    !!    actcoeff_n(1:nneutral) = -9999.9D0
    !!    actcoeff_c(1:nelectrol) = -9999.9D0
    !!    actcoeff_a(1:nelectrol) = -9999.9D0
    !!    meanmolalactcoeff(1:nelectrol) = -9999.9D0
    !!    ionactivityprod(1:nelectrol) = -9999.9D0
    !!    lnmeanmactcoeff = -9999.9D0
    !!    RETURN
    !!ENDIF

    !loop over all identified electrolyte components and calculate the corresponding mean molal activity coefficient and molal ion activity product:
    ii = 0
    ionactivityprod = 0.0D0
    meanmolalactcoeff = 1.0D0
    lnmeanmactcoeff = 0.0D0
    DO ii = 1,nelectrol
        ic = ElectComps(ii,1)       !get cation identifier 
        ia = ElectComps(ii,2)       !get anion identifier
        i = CatNr(ic)               !array index number of cation ic in, e.g., SMC array
        k = AnNr(ia) 
        IF (SMC(i) > 0.0D0 .AND. SMA(k) > 0.0D0) THEN
            IF (actcoeff_c(i) > 0.0D0 .AND. actcoeff_a(k) > 0.0D0) THEN
                nuec = REAL(ElectNues(ii,1), KIND=8)
                nuea = REAL(ElectNues(ii,2), KIND=8)
                nuelngc = nuec*lnactcoeff_c(i)
                nuelnga = nuea*lnactcoeff_a(k)
                ma1 = (nuelngc + nuelnga)/(nuec+nuea)
                lnmeanmactcoeff(ii) = ma1
                IF (ma1 > lntiny) THEN                          !no floating point underflow problem expected
                    IF (ma1 < lnhuge) THEN                      !no floating point overflow problem expected
                        meanmolalactcoeff(ii) = EXP(ma1)        !the mean molal activity coefficient of the ions of "electrolyte unit" ii
                    ELSE                                        !numerical issue, so output a large number (smaller than overflow risk)
                        trunc = 0.1D0*LOG( ABS(ma1) )
                        meanmolalactcoeff(ii) = EXP(lnhuge + trunc)
                        errorflagcalc = 7
                    ENDIF
                ELSE !underflow risk
                    meanmolalactcoeff(ii) = EXP(lntiny)         !i.e. tiny number, almost zero
                    errorflagcalc = 7
                ENDIF
                t1 = nuelngc + nuelnga + nuec*LOG(SMC(i)) + nuea*LOG(SMA(k))
                IF (t1 > lntiny) THEN                           !no floating point underflow problem expected
                    IF (t1 < lnhuge) THEN                       !no floating point overflow problem expected
                        ionactivityprod(ii) = EXP(t1)           !the mean molal activity coefficient of the ions of "electrolyte unit" ii
                    ELSE                                        !numerical issue, so output a large number (smaller than overflow risk)
                        trunc = 0.1D0*LOG( ABS(t1) )
                        ionactivityprod(ii) = EXP(lnhuge + trunc)
                        errorflagcalc = 7
                    ENDIF
                ELSE                                            !underflow risk
                    trunc = 0.1D0*LOG( ABS(t1) )
                    ionactivityprod(ii) = EXP(lntiny - trunc)   !i.e. tiny number, almost zero
                    errorflagcalc = 7
                ENDIF
            ENDIF
        ENDIF
    ENDDO !ii
    activity(nneutral+1:nindcomp) = ionactivityprod(1:nelectrol)

    END SUBROUTINE AIOMFAC_calc 
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
    !*   -> latest changes: 2018/05/27                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    SUBROUTINE Gammas()
    
    USE ModSystemProp, ONLY : nneutral, NKNpNGS, Ncation, Nanion, Mmass, bisulfsyst, bicarbsyst
    USE ModSRunifac, ONLY : SRsystm, SRunifac
    USE ModMRpart, ONLY : LR_MR_activity, GammaCO2

    IMPLICIT NONE
    !local variables:
    INTEGER(4) :: NKNpNcat, nnp1
    REAL(8) :: summolal, sumXN
    REAL(8),DIMENSION(NKNpNGS) :: lnGaSR
    REAL(8),DIMENSION(nneutral) :: wtfbycompMW, mNeutral
    LOGICAL(4) :: refreshgref
    !........................................................................................................
    
    nnp1 = nneutral +1
    !initialize the mole fraction composition arrays
    X = 0.0D0
    XN = 0.0D0
    !Calculation of the electrolyte-free mole fractions of the neutral components, XN:
    wtfbycompMW(1:nneutral) = wtf(1:nneutral)/Mmass(1:nneutral)
    sumXN = SUM(wtfbycompMW)
    XN(1:nneutral) = wtfbycompMW/sumXN
    
    meanSolventMW = SUM(Mmass(1:nneutral)*XN(1:nneutral)) !mean molecular mass of the solvent mixture (non-electrolytes)
    mNeutral = XN(1:nneutral)/meanSolventMW  ![mol/kg], the molality of the neutrals

    !Addition of the moles of substance (neutral and ionic) per 1 kg of electrolyte-free solvent mixture
    IF (bisulfsyst .OR. bicarbsyst) THEN    !update since it changes in DiffKsulfuricDissoc and DiffKcarbonateDissoc
        SumIonMolalities = SUM(SMA(1:Nanion)) +SUM(SMC(1:Ncation))
    ENDIF
    summolal = SUM(mNeutral) +SumIonMolalities !sum of all molalities

    !Calculation of the mole fraction (X) of the neutral components and the ions with respect to dissociated electrolytes/ions.
    !==> the structure of the mole fraction array X is: 
    !1) neutral components in component order,
    !2) ions: first the cations, then the anions
    X(1:nneutral) = mNeutral/summolal                          !Mole fraction of the neutral components (on the basis of dissociated electrolytes)    
    NKNpNcat = nneutral + Ncation
    X(nnp1:NKNpNcat) = SMC(1:Ncation)/summolal                 !Mole fractions of the cations
    X(NKNpNcat+1:NKNpNcat+Nanion) = SMA(1:Nanion)/summolal     !Mole fractios of the anions

    !check whether temperature-dependent parameters need to be updated in SR and LR parts:
    IF (ABS(T_K - lastTK) > 1.0D-2) THEN !detected a change in temperature --> PsiT and other coeff. need to be updated
        refreshgref = .true.
        DebyeHrefresh = .true.
        lastTK = T_K
    ELSE
        refreshgref = .false.
        DebyeHrefresh = .false.
    ENDIF
    !calculate the LR and MR activity coefficient contributions (gammas):  
    CALL LR_MR_activity()
    
    !call the UNIFAC model part for the short-range interaction contributions.
    CALL SRunifac(NKNpNGS, T_K, X, XN, refreshgref, lnGaSR) 
    gnsrln(1:nneutral) = lnGaSR(1:nneutral)
    gcsrln(1:Ncation)  = lnGaSR(nnp1:NKNpNcat)
    gasrln(1:Nanion)   = lnGaSR(NKNpNcat+1:NKNpNcat+Nanion)

    !update activity coefficient of CO2(aq):
    IF (bicarbsyst) THEN
        CALL GammaCO2()
    ENDIF

    END SUBROUTINE Gammas 
    !==========================================================================================

END MODULE ModCalcActCoeff