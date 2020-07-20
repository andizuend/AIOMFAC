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
!*   -> latest changes: 2018/05/28                                                      *
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
!*   -  SUBROUTINE HSO4DissociationCheck                                                *
!*                 CONTAINS:  -  SUBROUTINE DiffKsulfuricDissoc                         *
!*                            -  FUNCTION   fHSO4dissoc                                 *
!*   -  SUBROUTINE DeltaActivities                                                      *
!*   -  SUBROUTINE partialdactcoeff                                                     *
!*                                                                                      *
!****************************************************************************************
MODULE ModCalcActCoeff

USE ModAIOMFACvar !provide access to all composition-dependent AIOMFAC variables

IMPLICIT NONE

!==========================================================================================================================
    CONTAINS
!==========================================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Determine activity coefficients and activities from given input wtf composition    * 
    !*   for the current system. The LR, MR, & SR activity coeff. contributions are         *
    !*   calculated by the call to 'Gammas' and/or HSO4DissociationCheck.                   * 
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
    !*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA                *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2004                                                            *
    !*   -> latest changes: 2018/05/23                                                      *
    !*                                                                                      *
    !****************************************************************************************
    SUBROUTINE AIOMFAC_calc(WTFin, TKelvin) 
    
    !Public Variables:
    USE ModSystemProp, ONLY : anNr, catNr, ElectComps, ElectNues, Ianion, Ication, &
        Nanion, Ncation, nelectrol, nindcomp, nneutral, SolvMixRefnd, bisulfsyst, errorflagcalc
    USE ModCompScaleConversion, ONLY : MassFrac2IonMolalities

    IMPLICIT NONE
    !interface variables declarations
    REAL(8),DIMENSION(nindcomp),INTENT(IN) :: WTFin
    REAL(8),INTENT(IN) :: TKelvin
    !local variable declarations:
    INTEGER(4) :: I, ii, k, ic, ia
    !..
    REAL(8) :: t1, ma1, nuelngc, nuelnga, nuec, nuea, trunc
    REAL(8),PARAMETER :: dtiny = SQRT(TINY(1.0D0))
    REAL(8),PARAMETER :: lntiny = 0.49D0*LOG(TINY(1.0D0))
    REAL(8),PARAMETER :: lnhuge = 0.49D0*LOG(HUGE(1.0D0))
    !........................................................................................................

    wtf = WTFin  !input concentration is in mass fractions; other concentration scales are derived (when needed) from this.
    T_K = TKelvin
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

    !convert mass fraction input to ion molalities for use in MR and LR parts:
    CALL MassFrac2IonMolalities(wtf, SMC, SMA)
    SumIonMolalities = SUM(SMA(1:Nanion)) +SUM(SMC(1:Ncation))
    IF (bisulfsyst) THEN
        !Perform the dissociation check for bisulfate-system ions and, via this subroutine, the activity coefficient calculations:
        CALL HSO4DissociationCheck() !many variables are referenced through the module ModAIOMFACvar
    ELSE
        alphaHSO4 = -9.999999D0
        CALL Gammas() !No dissociation of bisulfate in mixture, so call Gammas for activity coefficient calculation separately.
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
            IF (lnactcoeff_n(I) < -340.0D0 .OR. lnactcoeff_n(I) > 340.0D0) THEN 
                trunc = AINT(lnactcoeff_n(I)) !get the truncated real value of lnactcoeff_n(I)
                IF (lnactcoeff_n(I) < -340.0D0) THEN
                    actcoeff_n(I) = EXP(-340.0D0 + 3.0D0*(lnactcoeff_n(I) - trunc))    
                ELSE
                    actcoeff_n(I) = EXP(340.0D0 + 3.0D0*(lnactcoeff_n(I) - trunc))
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
            IF (lnactcoeff_c(I) < -340.0D0 .OR. lnactcoeff_c(I) > 340.0D0) THEN 
                trunc = AINT(lnactcoeff_c(I)) !get the truncated real value of lnactcoeff_c(I)
                IF (lnactcoeff_c(I) < -340.0D0) THEN
                    actcoeff_c(I) = EXP(-340.0D0 + 3.0D0*(lnactcoeff_c(I) -trunc))    
                ELSE
                    actcoeff_c(I) = EXP(340.0D0 + 3.0D0*(lnactcoeff_c(I) -trunc))
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
                lnactcoeff_a(I) = lnactcoeff_a(I) +solvmixcorrMRa(I)+ Tmolal -TmolalSolvmix 
            ENDIF
            IF (lnactcoeff_a(I) < -340.0D0 .OR. lnactcoeff_a(I) > 340.0D0) THEN 
                trunc = AINT(lnactcoeff_a(I)) !get the truncated real value of lnactcoeff_a(I)
                IF (lnactcoeff_a(I) < -340.0D0) THEN
                    actcoeff_a(I) = EXP(-340.0D0 + 3.0D0*(lnactcoeff_a(I) -trunc))  
                ELSE
                    actcoeff_a(I) = EXP(340.0D0 + 3.0D0*(lnactcoeff_a(I) -trunc)) 
                ENDIF
            ELSE
                actcoeff_a(I) = EXP(lnactcoeff_a(I))
            ENDIF
        ENDIF
        !save activity coefficient and molality of this ion in an array by actual AIOMFAC ion index:
        ii = Ianion(I)
        actcoeff_ion(ii) = actcoeff_a(I)
        molality_ion(ii) = SMA(I)
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
    DO  ii = 1,nelectrol
        ic = ElectComps(ii,1)  !get cation identifier 
        ia = ElectComps(ii,2)  !get anion identifier
        i = CatNr(ic)          !array index number of cation ic in, e.g., SMC array
        k = AnNr(ia) 
        IF (SMC(i) > 0.0D0 .AND. SMA(k) > 0.0D0) THEN
            IF (actcoeff_c(i) > 0.0D0 .AND. actcoeff_a(k) > 0.0D0) THEN
                nuec = REAL(ElectNues(ii,1), KIND=8)
                nuea = REAL(ElectNues(ii,2), KIND=8)
                nuelngc = nuec*lnactcoeff_c(i)
                nuelnga = nuea*lnactcoeff_a(k)
                ma1 = (nuelngc + nuelnga)/(nuec+nuea)
                lnmeanmactcoeff(ii) = ma1
                IF (ma1 > lntiny) THEN !no floating point underflow problem expected
                    IF (ma1 < lnhuge) THEN !no floating point overflow problem expected
                        meanmolalactcoeff(ii) = EXP(ma1) !the mean molal activity coefficient of the ions of "electrolyte unit" ii
                    ELSE !numerical issue, so output a large number (smaller than overflow risk)
                        trunc = AINT(ma1)
                        meanmolalactcoeff(ii) = EXP(lnhuge + 3.0D0*(ma1 - trunc))
                        errorflagcalc = 7
                    ENDIF
                ELSE !underflow risk
                    meanmolalactcoeff(ii) = EXP(lntiny) !i.e. tiny number, almost zero
                    errorflagcalc = 7
                ENDIF
                t1 = nuelngc + nuelnga + nuec*LOG(SMC(i)) + nuea*LOG(SMA(k))
                IF (t1 > lntiny) THEN !no floating point underflow problem expected
                    IF (t1 < lnhuge) THEN !no floating point overflow problem expected
                        ionactivityprod(ii) = EXP(t1) !the mean molal activity coefficient of the ions of "electrolyte unit" ii
                    ELSE !numerical issue, so output a large number (smaller than overflow risk)
                        trunc = AINT(t1)
                        ionactivityprod(ii) = EXP(lnhuge + 3.0D0*(t1 - trunc))
                        errorflagcalc = 7
                    ENDIF
                ELSE !underflow risk
                    trunc = AINT(t1)
                    ionactivityprod(ii) = EXP(lntiny + 3.0D0*(t1 - trunc)) !i.e. tiny number, almost zero
                    errorflagcalc = 7
                ENDIF
            ENDIF
        ENDIF
    ENDDO !ii
    activity(nneutral+1:nindcomp) = ionactivityprod(1:nelectrol)

    END SUBROUTINE AIOMFAC_calc 
    !==========================================================================================================================
    

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
    
    !Public module variables used in Gammas():
    USE ModSystemProp, ONLY : nneutral, NKNpNGS, Ncation, Nanion, Mmass, bisulfsyst
    USE ModSRunifac, ONLY : SRsystm, SRunifac
    USE ModMRpart, ONLY : LR_MR_activity

    IMPLICIT NONE
    !Local Variables:
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
    IF (bisulfsyst) THEN !update since it changes in DiffKsulfuricDissoc
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

    END SUBROUTINE Gammas 
    !==========================================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   The degree of dissociation of bisulfate ions (HSO4- <--> H+ + SO4--) will be       *
    !*   determined iteratively for the given mixture and temperature.                      *
    !*                                                                                      *
    !*   This subroutine makes use of Brent's method for finding a root of an equation.     *
    !*   This third-party code is distributed under the GNU LGPL license (see brent.f90).   *
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
    SUBROUTINE HSO4DissociationCheck()

    USE ModSystemProp, ONLY : idH, idHSO4, idSO4, nneutral

    IMPLICIT NONE
    !Local variables and parameters:
    INTEGER(4) :: NTRY, nb
    REAL(8) :: mH, mHSO4, mHSO4min, mHSO4max, mSulfmax, mHmax
    !--
    REAL(8),PARAMETER :: dtiny = 1.0D-14, DEPS = EPSILON(1.0D0)
    REAL(8),PARAMETER :: T0 = 298.15D0                                                                                   ![K] reference temperature for dissociation constant temperature dependency calculation
    REAL(8),PARAMETER :: MolarMassHSO4 = 0.097071D0                                                                      !the molar mass of the bisulfate ion in [kg/mol]
    REAL(8),PARAMETER :: c1 = 2.63468D-01, c2 = 2.63204D-01, c3 = 1.72067D+00, c4 = 7.30039D+01, c5 = 7.45718D-01        !fitted parameters for the estimation of alphaHSO4pred based on scaled mass fraction of HSO4max
    REAL(8),PARAMETER :: DH0 = -1.8554748598D4, cp0 = -2.05594443486D2, dcpdT = -9.94992864240D-1                        !parameters for the temperature dependent bisulfate dissociation constant calculation.
    REAL(8),PARAMETER :: fc1 = 1.0D0/8.314472D0, fc2 = DH0-cp0*T0+0.5D0*dcpdT*T0*T0, fc3 = 1.0D0/T0, fc4 = cp0-dcpdT*T0  !terms/factors for the temperature dependent bisulfate dissociation constant calculation.
    REAL(8) :: mscale, diffK, mSulf, dm, dm2, mr, T, tol, lnKHSO4atT, lngKbisulf, wfHSO4max, alphaHSO4pred, &
        & mfHSO4maxinions, mHSO4Guess1, mHSO4Guess2, b, c
    REAL(8),DIMENSION(:),ALLOCATABLE :: mSb1, mSb2
    REAL(8),EXTERNAL :: brentzero
    !--
    LOGICAL(4) :: success, geomscal
    !................................................................................................

    !total [H+] molality: mH = 2*[H2SO4]+[HCl]+[HNO3]+[HBr] + ...; here this is the sum: [H+]+[HSO4-], 
    !because [H+] is already the sum of the molalities from all H+ containing electrolytes.
    mHmax = 0.0D0
    mHSO4max = 0.0D0
    mSulfmax = 0.0D0  !highest possible molality of SO4-- 
    mSulf = 0.0D0
    mHSO4 = 0.0D0
    mH = 0.0D0
    
    !calculate maximum possible molality amounts:
    !H+
    mHmax = mHmax + SMC(idH)
    mH = SMC(idH)
    !HSO4-
    mHmax = mHmax + SMA(idHSO4)
    mSulfmax = mSulfmax + SMA(idHSO4) ![SO4--] from [HSO4-]
    mHSO4 = SMA(idHSO4)  !initial molality of HSO4-
    !SO4--
    mSulfmax = mSulfmax + SMA(idSO4) ![SO4--]
    mSulf = SMA(idSO4)
    mHSO4max = mHSO4 +MIN(mH, mSulf)
    mHSO4min = 0.0D0 !this is always the case
    mscale = 1.0D-3*DEPS !minimal acceptable molality of any dissociation reagent (set to a small scaling value).

    !If no proton or no sulfate or no sulfuric acid are in the solution at the same time, then step over, otherwise calculate dissociation:
    IF (mHSO4max > mscale) THEN !@@**## there is some dissociation, so proceed...
        !calculate the reference value of the 2nd dissociation constant of H2SO4 as a function of temperature using parameterization by Knopf et al. (2003) (corrected equation 16). 
        !The parameterization is valid in the T-range [179 K < T_K < 474 K], for lower temperatures, replace T with 180 K temporarily and use the ln(KHSO4_T) 
        !of that temperature as that might be the best estimate for that case.
        T = T_K
        IF (T < 180.0D0) THEN !use T of 180 K for parameterization
            T = 180.0D0
        ELSE IF (T > 473.0D0) THEN !use T of 473 K
            T = 473.0D0
        ENDIF
        lnKHSO4atT = -4.54916799587D0 -fc1*(fc2*(1.0D0/T-fc3)-fc4*LOG(T/T0) -0.5D0*dcpdT*(T-T0))
        !..
        IF (mHSO4max-mscale > 0.0D0) THEN
            mHSO4max = mHSO4max-mscale
        ELSE
            mscale = mHSO4max*5.0D-1
            mHSO4max = mHSO4max-mscale
        ENDIF
    
        !Set stoichiometric limits such that none of ion molalities is exactly zero:
        IF (mHSO4 < mscale) THEN
            mHSO4 = mscale
            mH = MAX(mH -mscale, 0.0D0)
            mSulf = MAX(mSulf -mscale, 0.0D0)
        ENDIF
        IF (mH < mscale) THEN
            mscale = mscale*0.5D0
            mHSO4 = mHSO4 -mscale
            mH = mH +mscale
            mSulf = mSulf +mscale
        ENDIF
        IF (mSulf < mscale) THEN
            mscale = mscale*0.5D0
            mHSO4 = mHSO4 -mscale
            mH = mH +mscale
            mSulf = mSulf +mscale
        ENDIF
        !##! section to compute a good initial guess for the msulf1, msulf2 range based on a parameterization of alphaHSO4 as a function of mass frac. of HSO4(max) in the water portion of the solvent mixture.
        !(A) calculate the special normalized mass fraction of HSO4max in water plus HSO4max: wfHSO4max; only the water fractional amount associated with HSO4max is considered (excluding water associated with other ions in the system).
        mfHSO4maxinions = 2.0D0*mHSO4max/SumIonMolalities !mole fraction of (HSO4max plus counter cations) relative to all ion molar amounts.
        wfHSO4max = mHSO4max*MolarMassHSO4/(mHSO4max*MolarMassHSO4 + mfHSO4maxinions*wtf(1)/SUM(wtf(1:nneutral)))
        !(B) calculate the predicted degree of dissociation at this wfHSO4max value based on the parameterization:
        alphaHSO4pred = 1.0D0 - 1.0D0/(1.0D0 + (1.0D0/wfHSO4max**c1 - 1.0D0/wfHSO4max**c2))**c4 &
            & + c3*wfHSO4max**c5*(1.0D0-wfHSO4max)**1.75D0
        alphaHSO4pred = MIN(alphaHSO4pred, 1.0D0-mscale) !prevent values larger than 1.0D0
        alphaHSO4pred = MAX(alphaHSO4pred, mscale)
        !(C) calculate the lower and upper guesses for msulf based on the predicted dissociation degree:
        mHSO4Guess1 = (1.0D0 - MIN(alphaHSO4pred+0.15D0, 1.0D0))*mHSO4max
        mHSO4Guess2 = (1.0D0 - MAX(alphaHSO4pred-0.15D0, 0.0D0))*mHSO4max
        !##!
        success = .false.
        !(1a) search for root bracket for mSulf variable in guessed interval:
        ntry = 1
        geomscal = .false.
        ALLOCATE(mSb1(1:ntry+1), mSb2(1:ntry+1))
        CALL zerobracket_inwards(fHSO4dissoc, mHSO4Guess1, mHSO4Guess2, ntry, geomscal, nb, mSb1, mSb2, success)
        IF (.NOT. success) THEN !check lower bracket
            IF (mHSO4min < mHSO4Guess1-mscale) THEN !(1b) search in full interval from mSulfmin to mSulfmax:
                mHSO4Guess1 = mHSO4Guess1+mscale
                CALL zerobracket_inwards(fHSO4dissoc, mHSO4min, mHSO4Guess1, ntry, geomscal, nb, mSb1, mSb2, success)
            ENDIF
        ENDIF
        IF ((.NOT. success) ) THEN
            IF (mHSO4Guess2 < mHSO4max-mscale) THEN !check upper bracket !(1c) search in full interval from mSulfmin to mSulfmax:
                mHSO4Guess2 = mHSO4Guess2-mscale
                CALL zerobracket_inwards(fHSO4dissoc, mHSO4Guess2, mHSO4max, ntry, geomscal, nb, mSb1, mSb2, success)
            ENDIF
        ENDIF
        IF (success) THEN
            !(2) use Brent's method to find the root.
            mHSO4Guess1 = mSb1(1)
            mHSO4Guess2 = mSb2(1)
            tol = 2.2D0*DEPS
            mHSO4 = brentzero(mHSO4Guess1, mHSO4Guess2, DEPS, tol, fHSO4dissoc) !using Brent's method with fHSO4dissoc
            !--
            !check found zero and return best found mHSO4:
            CALL DiffKsulfuricDissoc(mHSO4, diffK)
        ELSE
            IF (alphaHSO4pred > 0.4D0) THEN
                mHSO4 = mHSO4min +alphaHSO4pred*mscale
                CALL DiffKsulfuricDissoc(mHSO4, diffK)
                dm = mSulf+mHSO4+mH
                mHSO4Guess1 = mHSO4
            ELSE
                mHSO4Guess1 = (1.0D0 - MIN(alphaHSO4pred, 1.0D0))*mHSO4max
                CALL DiffKsulfuricDissoc(mHSO4, diffK)
                dm = mSulf+mHSO4+mH
                mHSO4Guess1 = mHSO4
            ENDIF
            !now solve a quadratic equation, keeping the ion activity coeff. constant to determine the exact bisulfate-system ion molalities.
            IF (lngKbisulf > 300.0D0 .OR. lngKbisulf < -300.0D0) THEN !apply floating point overflow protection measures
                lngKbisulf = SIGN(300.0D0, lngKbisulf)
            ENDIF
            mr = EXP(lnKHSO4atT - lngKbisulf)
            b = -(SMC(idH) +SMA(idSO4) +mr)
            c = SMC(idH)*SMA(idSO4) -mr*SMA(idHSO4)
            tol = b*b -4.0D0*c
            IF (tol >= 0.0D0) THEN
                tol = SQRT(tol)
                !use a numerically stable approach to find the two roots:
                IF (b >= 0.0D0) THEN
                    dm = 0.5D0*(-b -tol)
                    dm2 = 2.0D0*c/(-b-tol)
                ELSE
                    dm = 2.0D0*c/(-b+tol)
                    dm2 = 0.5D0*(-b +tol)
                ENDIF
                IF (ABS(dm2) < ABS(dm)) THEN
                    dm = dm2
                ENDIF
                mSulf = SMA(idSO4)-dm
                mH = SMC(idH)-dm
                mHSO4 = SMA(idHSO4)+dm
                CALL DiffKsulfuricDissoc(mHSO4, diffK)
                alphaHSO4 = 1.0D0-(mHSO4/mHSO4max)
                IF ((ABS(diffK) > 1.0D-1 .AND. ABS(alphaHSO4 -alphaHSO4pred) > 0.4D0) &
                    & .OR. ABS(alphaHSO4 -alphaHSO4pred) > 0.99D0) THEN
                    CALL DiffKsulfuricDissoc(mHSO4Guess1, diffK)
                    mHSO4 = SMA(idHSO4)
                ENDIF
            ENDIF
        ENDIF
        alphaHSO4 = 1.0D0-(mHSO4/mHSO4max) !this is the generally valid expression for alphaHSO4
        diffKHSO4 = ABS(diffK)
        !--
    ELSE !case where HSO4 is part of the system components, but at present input composition, mHSO4max is at zero molality.
        alphaHSO4 = -9.999999D0
        CALL Gammas() !call Gammas once
    ENDIF !@@**##
    !================================================================================================= 

    CONTAINS
    
        !****************************************************************************************
        !*   :: Purpose ::                                                                      *
        !*   Subroutine to calculate the difference between the reference value (lnKHSO4atT)    *
        !*   and the estimated value of the dissociation constant of HSO4- <-> H+ + SO4--.      *
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
        SUBROUTINE DiffKsulfuricDissoc(mHSO4inp, diffK)
               
        !Public Variables:
        USE ModSystemProp, ONLY : idH, idHSO4, idSO4, topsubno

        IMPLICIT NONE

        REAL(8),INTENT(INOUT) :: mHSO4inp
        REAL(8),INTENT(OUT) :: diffK
        !local vars:
        REAL(8),PARAMETER :: dtiny = 9.1D0*EPSILON(1.0D0)
        REAL(8) :: lngammaH, lngammaSulf, lngammaHydSulf, mscale, mratio, madd, mSulf
        !................................................................

        !initial values:
        mscale = 1.0D-3*DEPS  !scaling factor

        !limit mHSO4 to the possible interval bounds:
        mHSO4 = mHSO4inp
        mHSO4 = MAX(mHSO4, mHSO4min)
        mHSO4 = MIN(mHSO4, mHSO4max)

        !Calculate the molalities of the variable ions:
        !HSO4-
        IF (mHSO4 < mscale) THEN !make sure that at least a tiny amount of HSO4- is present
            madd = MIN(mscale -mHSO4, mHSO4max)
            mHSO4 = mHSO4 +madd
        ENDIF
        SMA(idHSO4) = mHSO4
        !H+
        mH = MAX(mHmax-mHSO4, 0.0D0)
        IF (mH < mscale) THEN !make sure that at least a tiny amount of H+ is present
            madd = mscale - mH !the difference from the tiny mH value to mscale
            IF (mHSO4-madd > mHSO4min) THEN
                mHSO4 = mHSO4 -madd
                SMA(idHSO4) = mHSO4
                mH = mH +madd
            ENDIF
        ENDIF
        SMC(idH) = mH
        !SO4--
        mSulf = MAX(mSulfmax-mHSO4, 0.0D0) ![SO4--]
        IF (mSulf < mscale) THEN !make sure that at least a tiny amount of H+ is present
            madd = mscale - mSulf !the difference from the tiny mH value to mscale
            IF (mHSO4-madd > mHSO4min) THEN
                mHSO4 = mHSO4 -madd
                SMA(idHSO4) = mHSO4
                mSulf = mSulf +madd
            ENDIF
        ENDIF
        SMA(idSO4) = mSulf

        IF (mH*mSulf*mHSO4 < -dtiny) THEN  !check the values for mistakes in program
            !$OMP CRITICAL (DKD2)
            WRITE(*,*) "WARNING: negative molality values detected! wtf(water): ", wtf(1)
            WRITE(*,*) "mH, mSulf, mHSO4: ", mH, mSulf, mHSO4
            WRITE(*,*) "mHSO4inp ", mHSO4inp
            WRITE(*,*) ""
            !$OMP END CRITICAL (DKD2)
        ENDIF
        !*-*
        CALL Gammas() !compute the natural log of activity coefficients from SR, MR, & LR parts.
        !*-*
        !Use the mole fraction-based ion activity coeff. (or their natural log).
        IF (idH > 0) THEN
            lngammaH = -Tmolal +gcsrln(idH) +gcmrln(idH) +gclrln(idH)
        ELSE
            lngammaH = 0.0D0
        ENDIF
        IF (idHSO4 > 0) THEN
            lngammaHydSulf = -Tmolal +gasrln(idHSO4) +gamrln(idHSO4) +galrln(idHSO4)
        ELSE
            lngammaHydSulf = 0.0D0
        ENDIF
        IF (idSO4 > 0) THEN
            lngammaSulf = -Tmolal +gasrln(idSO4) +gamrln(idSO4) +galrln(idSO4)
        ELSE
            lngammaSulf = 0.0D0
        ENDIF
        !estimate the molalities of SO4--, HSO4- and H+ to get agreement 
        !with the reference value of the dissociation constant.
        lngKbisulf = lngammaH +lngammaSulf -lngammaHydSulf

        !comparison, logarithmic "rel." deviation:
        IF (mHSO4 > 0.0D0) THEN
            mratio = mH*mSulf/mHSO4
            IF (mratio > 0.0D0) THEN 
                diffK = (lngKbisulf +LOG(mratio)) -lnKHSO4atT
            ELSE
                mratio = DEPS/mHSO4
                diffK = (lngKbisulf +LOG(mratio)) -lnKHSO4atT
            ENDIF
        ELSE
            mratio = (mH*mSulf)/DEPS
            IF (mratio > 0.0D0) THEN 
                diffK = (lngKbisulf +LOG(mratio)) -lnKHSO4atT
            ELSE
                mratio = DEPS
                diffK = (lngKbisulf +LOG(mratio)) -lnKHSO4atT
            ENDIF
        ENDIF

        END SUBROUTINE DiffKsulfuricDissoc 
        !--------------------------------------------------------------------------------------------------
    
        FUNCTION fHSO4dissoc(mHSO4in) !internal wrapper function for the computation of Diffk as a function of mHSO4 input.
    
        IMPLICIT NONE
        !interface:
        REAL(8), INTENT(INOUT) :: mHSO4in
        REAL(8) :: fHSO4dissoc
        !......................................
        !call DiffKsulfuricDissoc to calculate the actual value of the dissociation (non)-equilibrium:
        CALL DiffKsulfuricDissoc(mHSO4in, fHSO4dissoc)
    
        END FUNCTION fHSO4dissoc
        !--------------------------------------------------------------------------------------------------

    END SUBROUTINE HSO4DissociationCheck
    !==========================================================================================================================
    
    
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
    !*   -> latest changes: 2019/09/24                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    SUBROUTINE DeltaActivities(xin, TKelvin, onlyDeltaVisc, dact, dactcoeff)

    USE ModSystemProp, ONLY : nindcomp
    USE ModAIOMFACvar, ONLY : deltaetamix

    IMPLICIT NONE

    !interface variables:
    REAL(8),DIMENSION(nindcomp),INTENT(IN) :: xin
    REAL(8),INTENT(IN) :: TKelvin
    LOGICAL(4),INTENT(IN) :: onlyDeltaVisc !default = .false.; this is set true when only the first component is varied; use only for delta(etamix) calculation.
    REAL(8),DIMENSION(nindcomp),INTENT(OUT) :: dact, dactcoeff
    !local variables:
    INTEGER(4) :: ni, j
    LOGICAL(4) :: xinputtype
    REAL(8),DIMENSION(nindcomp) :: xinit
    REAL(8),DIMENSION(:,:),ALLOCATABLE :: partdact, partdactcoeff
    REAL(8) :: Dtiny, dn, partdact_ji, partdactcoeff_ji
    PARAMETER (Dtiny = 1.79D1*EPSILON(1.0D0))
    PARAMETER (dn = 0.1D0*SQRT(EPSILON(1.0D0))) ![mol] a small, but not too small molar (change) for the numerical differentiation
    PARAMETER (xinputtype = .true.) !mole fraction input concentration
    !...........................................................
    !parameters:
    xinit = xin  !xinit is locally stored and unaffected by changes when reentrant code will be called that might feed back on xin!

    !calculate the finite differences (numerical "partial derivatives") of the component ni's activity at composition xinit while holding the moles of the other componentes fixed.
    !compositions for the forward / backward differences at the given point xinit (component nnvar is the one to be enhanced/diminished):
    IF (onlyDeltaVisc) THEN !e.g. for web-model version, only perform calculation with variation in component 1 (often water) and effect on component 1 only for deltaetamix etc (much faster).
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
            dact(ni) = SUM(ABS(partdact(ni,1:nindcomp))) !we skip here the *dn (which would allow calculating the total differential, as it leads to small numbers and sensitivity will be normalized later anyways)
            dactcoeff(ni) = SUM(ABS(partdactcoeff(ni,1:nindcomp)))
        ENDDO
    ENDIF
    !for log viscosity sensitivity estimate as a function of a change in mole fraction of component 1 (water):
    deltaetamix = ABS(deltaetamix)
    
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
    !*   -> latest changes: 2019/09/24                                                      *
    !*                                                                                      *
    !****************************************************************************************
    SUBROUTINE partialdactcoeff(xin, TKelvin, j, i, partdact_ji, partdactcoeff_ji)

    USE ModSystemProp, ONLY : nindcomp, nneutral, nd, Mmass, errorflagcalc
    USE ModCompScaleConversion
    USE ModAIOMFACvar, ONLY : activity, meanmolalactcoeff, actcoeff_n, deltaetamix, etamix
    
    IMPLICIT NONE

    !interface variables:
    INTEGER(4),INTENT(IN) :: j, i
    REAL(8),DIMENSION(nindcomp),INTENT(IN) :: xin
    REAL(8),INTENT(IN) :: TKelvin
    REAL(8),INTENT(OUT) :: partdact_ji, partdactcoeff_ji
    !local variables:
    LOGICAL(4) :: xinputtype
    REAL(8),DIMENSION(nindcomp) :: xinit, xplus, wtf
    REAL(8) :: ntplus, Dtiny, dn, actplus, actinit, actcoeffplus, actcoeffinit, etamixplus, etamixinit
    PARAMETER (Dtiny = 1.11D0*EPSILON(1.0D0))
    PARAMETER (dn = 0.1D0*SQRT(EPSILON(1.0D0))) ![mol] a small molar (change) for the numerical differentiation
    PARAMETER (xinputtype = .true.) !mole fraction input concentration
    !...........................................................
    !parameters:
    xinit = xin !xinit is locally stored and unaffected by changes when reentrant code will be called that might feed back on xin!
    
    !calculate the finite differences (numerical "partial derivatives") of component j's activity and activity coeff. with respect to component i,
    !at composition xinit, while holding the moles of the other components fixed.
    !compositions for the forward / backward differences at the given point xinit (component i is the one to be enhanced/diminished):
    ntplus = 1.0D0+dn !as the molar sum of all components can be set arbitrarily to 1.0
    !calculate the mole fractions at the slightly changed compositions:
    xplus(1:i-1) = xinit(1:i-1)/ntplus
    xplus(i) = (xinit(i)+dn)/ntplus
    xplus(i+1:nindcomp) = xinit(i+1:nindcomp)/ntplus
    !..
    !Call AIOMFAC_calc to calculate the activity at composition xplus:
    CALL MoleFrac2MassFrac(xplus, Mmass, wtf) 
    CALL AIOMFAC_calc(wtf, TKelvin)
    actplus = activity(j)
    IF (j <= nneutral) THEN !neutral component
        actcoeffplus = actcoeff_n(j)
    ELSE !electrolyte component
        actcoeffplus = meanmolalactcoeff(j-nneutral)
    ENDIF
    IF (j == 1 .AND. i == 1) THEN !save the etamix value for the sensitivity of viscosity
        etamixplus = etamix
    ENDIF
    !..
    !calculate the activity & activity coeff. at the initial composition xinit:
    CALL MoleFrac2MassFrac(xinit, Mmass, wtf) 
    CALL AIOMFAC_calc(wtf, TKelvin)
    actinit = activity(j)
    IF (j <= nneutral) THEN
        actcoeffinit = actcoeff_n(j)
    ELSE
        actcoeffinit = meanmolalactcoeff(j-nneutral)
    ENDIF
    IF (j == 1 .AND. i == 1) THEN !save the etamix value for the sensitivity of viscosity
        etamixinit = etamix
    ENDIF
    !..
    !calculate the numerical forward differences with respect to activity at point xinit 
    !(partial derivative of component j with resp. to a small molar change, dn, in component i)
    IF (actplus > 0.0D0 .AND. actinit > 0.0D0) THEN
        IF (ABS(actplus-actinit) < 1.0D-292) THEN !floating underflow risk
            IF (actplus-actinit < 0.0D0) THEN !"negative zero"
                partdact_ji = -1.0D-292 !almost zero
            ELSE !"positive zero"
                partdact_ji = 1.0D-292 !almost zero
            ENDIF
        ELSE !no underflow
            IF (actplus -actinit > 1.0D292) THEN !floating overflow risk
                partdact_ji = 1.0D292+(LOG(actplus)-LOG(actinit)) !set to huge positive number to allow following computations but prevent overflow
            ELSE IF (actplus -actinit < -1.0D292) THEN !floating overflow risk
                partdact_ji = -1.0D292-(LOG(actplus)-LOG(actinit)) !set to huge negative number to allow following computations but prevent overflow
            ELSE
                partdact_ji = (actplus -actinit)/dn !the numerical derivatives with respect to component i (while all other component's moles are kept const.) 
            ENDIF
        ENDIF
        IF (ABS(actcoeffplus -actcoeffinit) < 1.0D-292) THEN !floating underflow risk
            IF (actcoeffplus -actcoeffinit < 0.0D0) THEN !"negative zero"
                partdactcoeff_ji = -1.0D-292 !almost zero
            ELSE !"positive zero"
                partdactcoeff_ji = 1.0D-292 !almost zero
            ENDIF
        ELSE !ok, no underflow
            IF (actcoeffplus -actcoeffinit > 1.0D292) THEN !floating overflow risk
                partdactcoeff_ji = 1.0D292+(LOG(actcoeffplus)-LOG(actcoeffinit)) !set to huge positive number to allow following computations but prevent overflow
            ELSE IF (actcoeffplus -actcoeffinit < -1.0D292) THEN !floating overflow risk
                partdactcoeff_ji = -1.0D292-(LOG(actcoeffplus)-LOG(actcoeffinit)) !set to huge negative number to allow following computations but prevent overflow
            ELSE
                partdactcoeff_ji = (actcoeffplus -actcoeffinit)/dn
            ENDIF
        ENDIF
    ELSE IF (xinit(j) > Dtiny .AND. (actplus < Dtiny .OR. actinit < Dtiny)) THEN !there is some amount of a component in the mixture, but the model activity coefficient is too large (thus, very steep derivative)
        partdact_ji = 1.11111111D5 !a specific value indicating a floating point overflow in AIOMFAC_calc, while not suppressing the feedback of such a value.
        partdactcoeff_ji = 1.11111111D5
        errorflagcalc = 2
    ELSE !exceptions where problem occurred
        partdact_ji = 0.0D0
        partdactcoeff_ji = 0.0D0
        errorflagcalc = 3
    ENDIF
    !**--
    IF (j == 1 .AND. i == 1) THEN
        !partial forward difference for mixture viscosity:
        IF (etamixinit > 0.0D0) THEN
            deltaetamix = ((etamixplus -etamixinit)/dn)/etamixinit  != (LOG(etamixplus) -LOG(etamixinit))/dn;  relative partial differential error with respect to component i (while all other component's moles are kept const.)
        ELSE
            deltaetamix = 0.0D0
        ENDIF
    ENDIF
    !**--
    END SUBROUTINE partialdactcoeff
    !==========================================================================================================================
    
END MODULE ModCalcActCoeff