!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Submodule of module ModCalcActCoeff containing procedures used for the calculation *
!*   of aqueous electrolyte/ion dissociation equilibria, such as those involving        *
!*   sulfate--bisulfate ions and carbonate--bicarbonate--CO2(aq).                       *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend, Hang Yin,                                                              *
!*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2018    (originally as non-submodule part of ModCalcActCoeff)   *
!*   -> latest changes: 2021-12-04                                                      *
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
!*   :: List of subroutines and functions contained in this submodule:                  *
!*   --------------------------------------------------------------                     *
!*   -  SUBROUTINE HSO4_dissociation                                                    *
!*   -  SUBROUTINE DiffKsulfuricDissoc                                                  *
!*   -  FUNCTION   fHSO4dissoc                                                          *
!*   -  SUBROUTINE HSO4_and_HCO3_dissociation                                           *
!*   -  SUBROUTINE DiffK_carb_sulf                                                      *
!*   -  PURE ELEMENTAL SUBROUTINE rboundsCheck                                          *
!*   -  PURE FUNCTION sum_sorted                                                        *
!*                                                                                      *
!****************************************************************************************
SUBMODULE (ModCalcActCoeff) SubModDissociationEquil

USE ModSystemProp, ONLY : bisulfsyst, calcviscosity, errorflagcalc, frominpfile, idH, idHCO3, &
    & idCO3, idOH, idCO2, idHSO4, idSO4, idCa, nneutral, NGI 
!USE ModAIOMFACvar, ONLY : alphaHSO4, diffKHSO4, SMA, SMC, SumIonMolalities, T_K, wtf   !accessible via parent module use statement
USE ModCompScaleConversion, ONLY : MassFrac2SolvMolalities, MassFrac2IonMolalities, Moles2solvmass, &
    & Molality2SolvMoleFrac

IMPLICIT NONE
!submodule parameters and variables:
INTEGER(4),DIMENSION(:),ALLOCATABLE :: map_solver2molarinp      !array storing the indices for mapping between the active set of solver variables and the "molarinp" variables;
!..
REAL(8),PARAMETER :: deps = EPSILON(1.0D0), sqrtdeps = SQRT(deps)
REAL(8),PARAMETER :: ztiny = TINY(1.0D0)**0.3D0, lztiny = 4.0D0*deps
REAL(8),PARAMETER :: minlim = LOG(1.0D-6*deps)    
REAL(8),PARAMETER :: Rgas_atm = 8.2057366081D-5                 !Rgas in units of [m^3.atm.K^-1.mol^-1]
!..
REAL(8),PARAMETER :: T0 = 298.15D0                              ![K] reference temperature for dissociation constant temperature dependency calculation
REAL(8),PARAMETER :: MolarMassHSO4 = 0.097071D0                                                                      !the molar mass of the bisulfate ion in [kg/mol]
REAL(8),PARAMETER :: c1 = 2.63468D-01, c2 = 2.63204D-01, c3 = 1.72067D+00, c4 = 7.30039D+01, c5 = 7.45718D-01        !fitted parameters for the estimation of alphaHSO4pred based on scaled mass fraction of HSO4max
REAL(8),PARAMETER :: DH0 = -1.8554748598D4, cp0 = -2.05594443486D2, dcpdT = -9.94992864240D-1                        !parameters for the temperature dependent bisulfate dissociation constant calculation.
REAL(8),PARAMETER :: fe1 = 1.0D0/8.314472D0, fe2 = DH0-cp0*T0+0.5D0*dcpdT*T0*T0, fe3 = 1.0D0/T0, fe4 = cp0-dcpdT*T0  !terms/factors for the temperature dependent bisulfate dissociation constant calculation.
!..
REAL(8) :: lngbisulf, lngK1bicarb, lngK2bicarb, lngKbisulf, lngKOH, lnK1HCO3atT, lnK2HCO3atT, lnKCO2atT, lnKH2OatT, &
    & lnKHSO4atT, mHmax, mHSO4max, mHSO4min, mSulfmax, nCarbmax, nCO2_init, nCO2gas, nCO2max, nH2O_init, nHCO3max,  &
    & nHmax, nHSO4max, nOH_init, nOHmax, nSulfmax, ntiny, nVmax, target_CO2gas_ppm, V
REAL(8) :: mCarb, mCO2, mH, mHCO3, mHSO4, mOH, mSulf            !molalities for activity calculations
REAL(8) :: nCarb, nCO2, nH, nHCO3, nHSO4, nOH, nSulf            !molar amounts for composition determination
REAL(8),DIMENSION(8) :: ln_maxval_inp
REAL(8),DIMENSION(:),ALLOCATABLE :: mNeutral, molNeutral        !molalities and moles of neutral species
REAL(8),DIMENSION(:),ALLOCATABLE :: SNA, SNC                    !moles of ions
REAL(8),DIMENSION(:),ALLOCATABLE :: solve_var, solve_var_maxval, solve_var_saved
LOGICAL(4) :: use_CO2gas_equil

!$OMP THREADPRIVATE(use_CO2gas_equil, solve_var, solve_var_maxval, solve_var_saved, map_solver2molarinp, mCarb, mCO2,   &
    !$OMP & mH, mHCO3, mHSO4, mOH, mSulf, nCarb, nCO2, nH, nHCO3, nHSO4, nOH, nSulf, mNeutral, molNeutral, lngbisulf,   &
    !$OMP & lngK1bicarb, lngK2bicarb, lngKbisulf, lngKOH, lnK1HCO3atT, lnK2HCO3atT, lnKCO2atT, lnKH2OatT, lnKHSO4atT,   &
    !$OMP & ln_maxval_inp, mHmax, mHSO4max, mHSO4min, mSulfmax, nCarbmax, nCO2_init, nCO2gas, nCO2max, nH2O_init,       &
    !$OMP & nHCO3max, nHmax, nHSO4max, nOH_init, nOHmax, nSulfmax, ntiny, nVmax, SNA, SNC, target_CO2gas_ppm, V)

!======================
    CONTAINS
!======================  

    
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
    !*   -> latest changes: 2021-11-27                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    MODULE SUBROUTINE HSO4_dissociation()

    IMPLICIT NONE
    !Local variables and parameters:
    INTEGER(4) :: ntry, nb
    REAL(8) :: alphaHSO4pred, b, c, diffK, dm, dm2, mfHSO4maxinions, mHSO4Guess1, mHSO4Guess2, mr, &
        & mscale, T, tol, wfHSO4max 
    REAL(8),DIMENSION(:),ALLOCATABLE :: mSb1, mSb2
    REAL(8),EXTERNAL :: brentzero   !Brent's method for finding a zero of a univariate function (external function)
    LOGICAL(4) :: success, geomscal, calcviscosity_saved
    !................................................................................................

    !total [H+] molality: mH = 2*[H2SO4]+[HCl]+[HNO3]+[HBr] + ...; here this is the sum: [H+]+[HSO4-], 
    !because [H+] is already the sum of the molalities from all H+ containing electrolytes.
    mHmax = 0.0D0
    mHSO4max = 0.0D0
    mSulfmax = 0.0D0                        !highest possible molality of SO4-- 
    mSulf = 0.0D0
    mHSO4 = 0.0D0
    mH = 0.0D0
    
    !determine maximum possible molality amounts:
    !H+
    mHmax = mHmax + SMC(idH)
    mH = SMC(idH)
    !HSO4-
    mHmax = mHmax + SMA(idHSO4)
    mSulfmax = mSulfmax + SMA(idHSO4)       ![SO4--] from [HSO4-]
    mHSO4 = SMA(idHSO4)                     !initial molality of HSO4-
    !SO4--
    mSulfmax = mSulfmax + SMA(idSO4)        ![SO4--]
    mSulf = SMA(idSO4)
    mHSO4max = mHSO4 +MIN(mH, mSulf)
    mHSO4min = 0.0D0                        !this is always the case
    mscale = 1.0D-3*deps                    !minimal acceptable molality of any dissociation reagent (small scaling value).

    !If no proton or no sulfate or no sulfuric acid are in the solution at the same time, 
    !then step over, otherwise calculate dissociation:
    IF (mHSO4max > mscale) THEN !@@**## there is some dissociation, so proceed...
        !calculate the reference value of the 2nd dissociation constant of H2SO4 as a function of temperature 
        !using parameterization by Knopf et al. (2003) (corrected equation 16). 
        !The parameterization is valid in the T-range [179 k < T_K < 474 k], for lower temperatures, 
        !replace T with 180 k temporarily and use the ln(KHSO4_T) 
        !of that temperature as that might be the best estimate for that case.
        T = T_K
        IF (T < 180.0D0) THEN !use T of 180 k for parameterization
            T = 180.0D0
        ELSE IF (T > 473.0D0) THEN !use T of 473 k
            T = 473.0D0
        ENDIF
        lnKHSO4atT = -4.54916799587D0 -fe1*(fe2*(1.0D0/T-fe3)-fe4*LOG(T/T0) -0.5D0*dcpdT*(T-T0))
        calcviscosity_saved = calcviscosity
        calcviscosity = .false.
        !..
        IF (mHSO4max-mscale > 0.0D0) THEN
            mHSO4max = mHSO4max-mscale
        ELSE
            mscale = mHSO4max*5.0D-1
            mHSO4max = mHSO4max-mscale
        ENDIF
    
        !Set stoichiometric limits such that none of the ion molalities is exactly zero:
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
        !
        !##! section to compute a good initial guess for the msulf1, msulf2 range based on a parameterization of 
        !alphaHSO4 as a function of mass frac. of HSO4(max) in the water portion of the solvent mixture.
        !(A) calculate the special normalized mass fraction of HSO4max in water plus HSO4max: wfHSO4max; 
        !only the water fractional amount associated with HSO4max is considered (excluding water associated with other ions in the system).
        mfHSO4maxinions = 2.0D0*mHSO4max/SumIonMolalities !mole fraction of (HSO4max plus counter cations) relative to all ion molar amounts.
        wfHSO4max = mHSO4max*MolarMassHSO4/(mHSO4max*MolarMassHSO4 + mfHSO4maxinions*wtf(1)/SUM(wtf(1:nneutral)))
        !(B) calculate the predicted degree of dissociation at this wfHSO4max value based on the parameterization:
        alphaHSO4pred = 1.0D0 - 1.0D0/(1.0D0 + (1.0D0/wfHSO4max**c1 - 1.0D0/wfHSO4max**c2))**c4 &
            & + c3*wfHSO4max**c5*(1.0D0 -wfHSO4max)**1.75D0
        alphaHSO4pred = MIN(alphaHSO4pred, 1.0D0 -mscale)        !prevent values larger than 1.0D0
        alphaHSO4pred = MAX(alphaHSO4pred, mscale)
        !(C) calculate the lower and upper guesses for msulf based on the predicted dissociation degree:
        mHSO4Guess1 = (1.0D0 - MIN(alphaHSO4pred +0.15D0, 1.0D0))*mHSO4max
        mHSO4Guess2 = (1.0D0 - MAX(alphaHSO4pred -0.15D0, 0.0D0))*mHSO4max
        !##!
        success = .false.
        !(1a) search for bracketed root:
        ntry = 1
        geomscal = .false.
        ALLOCATE(mSb1(1:ntry+1), mSb2(1:ntry+1))
        CALL zerobracket_inwards(fHSO4dissoc, mHSO4Guess1, mHSO4Guess2, ntry, geomscal, nb, mSb1, mSb2, success)
        IF (.NOT. success) THEN                         !check lower bracket
            IF (mHSO4min < mHSO4Guess1-mscale) THEN     !(1b) search in full interval from mSulfmin to mSulfmax:
                mHSO4Guess1 = mHSO4Guess1+mscale
                CALL zerobracket_inwards(fHSO4dissoc, mHSO4min, mHSO4Guess1, ntry, geomscal, nb, mSb1, mSb2, success)
            ENDIF
        ENDIF
        IF ((.NOT. success) ) THEN
            IF (mHSO4Guess2 < mHSO4max-mscale) THEN     !check upper bracket !(1c) search in full interval from mSulfmin to mSulfmax:
                mHSO4Guess2 = mHSO4Guess2-mscale
                CALL zerobracket_inwards(fHSO4dissoc, mHSO4Guess2, mHSO4max, ntry, geomscal, nb, mSb1, mSb2, success)
            ENDIF
        ENDIF
        IF (success) THEN
            !(2) use Brent's method to find the root.
            mHSO4Guess1 = mSb1(1)
            mHSO4Guess2 = mSb2(1)
            tol = sqrtdeps
            mHSO4 = brentzero(mHSO4Guess1, mHSO4Guess2, deps, tol, fHSO4dissoc)     !using Brent's method with fHSO4dissoc
            !check found zero and return best found mHSO4:
            calcviscosity = calcviscosity_saved
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
            !now solve a quadratic equation, keeping the ion activity coeff. constant 
            !to determine the exact bisulfate-system ion molalities.
            IF (lngKbisulf > 300.0D0 .OR. lngKbisulf < -300.0D0) THEN       !apply floating point overflow protection measures
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
                    dm2 = 2.0D0*c/(-b -tol)
                ELSE
                    dm = 2.0D0*c/(-b +tol)
                    dm2 = 0.5D0*(-b +tol)
                ENDIF
                IF (ABS(dm2) < ABS(dm)) THEN
                    dm = dm2
                ENDIF
                mSulf = SMA(idSO4)-dm
                mH = SMC(idH)-dm
                mHSO4 = SMA(idHSO4)+dm
                calcviscosity = calcviscosity_saved
                CALL DiffKsulfuricDissoc(mHSO4, diffK)
                alphaHSO4 = 1.0D0-(mHSO4/mHSO4max)
                IF ((ABS(diffK) > 1.0D-1 .AND. ABS(alphaHSO4 -alphaHSO4pred) > 0.4D0) &
                    & .OR. ABS(alphaHSO4 -alphaHSO4pred) > 0.99D0) THEN
                    CALL DiffKsulfuricDissoc(mHSO4Guess1, diffK)
                    mHSO4 = SMA(idHSO4)
                ENDIF
            ENDIF
        ENDIF
        alphaHSO4 = 1.0D0 - (mHSO4/mHSO4max)    !this is a generally valid expression for alphaHSO4
        diffKHSO4 = ABS(diffK)
        !--
    ELSE    !case where HSO4 is part of the system components, but at present input composition, mHSO4max is at zero molality.
        alphaHSO4 = -9.999999D0
        CALL Gammas()                           !call Gammas once
    ENDIF !@@**##

    END SUBROUTINE HSO4_dissociation
    !================================================================================================= 

    
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
    !*   -> latest changes: 2021-12-02                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    MODULE SUBROUTINE DiffKsulfuricDissoc(mHSO4inp, diffK)
               
    USE ModSystemProp, ONLY : idH, idHSO4, idSO4, topsubno
    USE ModAIOMFACvar, ONLY : gasrln, gamrln, galrln, gcsrln, gcmrln, gclrln, SMA, SMC, wtf, Tmolal

    IMPLICIT NONE
    !interface arguments:
    REAL(8),INTENT(INOUT) :: mHSO4inp
    REAL(8),INTENT(OUT) :: diffK
    !local vars:
    REAL(8),PARAMETER :: dtiny = 9.1D0*EPSILON(1.0D0)    
    REAL(8) :: lngammaH, lngammaSulf, lngammaHydSulf, mscale, mratio, madd
    !................................................................

    mscale = 1.0D-3*deps  !scaling factor
    !limit mHSO4 to the possible interval bounds:
    mHSO4 = mHSO4inp
    mHSO4 = MAX(mHSO4, mHSO4min)
    mHSO4 = MIN(mHSO4, mHSO4max)

    !Calculate the molalities of the variable ions:
    !HSO4-
    IF (mHSO4 < mscale) THEN                !make sure that at least a tiny amount of HSO4- is present
        madd = MIN(mscale -mHSO4, mHSO4max)
        mHSO4 = mHSO4 +madd
    ENDIF
    SMA(idHSO4) = mHSO4
    !H+
    mH = MAX(mHmax-mHSO4, 0.0D0)
    IF (mH < mscale) THEN                   !make sure that at least a tiny amount of H+ is present
        madd = mscale - mH                  !the difference from the tiny mH value to mscale
        IF (mHSO4-madd > mHSO4min) THEN
            mHSO4 = mHSO4 -madd
            SMA(idHSO4) = mHSO4
            mH = mH +madd
        ENDIF
    ENDIF
    SMC(idH) = mH
    !SO4--
    mSulf = MAX(mSulfmax-mHSO4, 0.0D0)
    IF (mSulf < mscale) THEN                !make sure that at least a tiny amount of SO4-- is present
        madd = mscale - mSulf               !the difference from the tiny mSulf value to mscale
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
            mratio = deps/mHSO4
            diffK = (lngKbisulf +LOG(mratio)) -lnKHSO4atT
        ENDIF
    ELSE
        mratio = (mH*mSulf)/deps
        IF (mratio > 0.0D0) THEN 
            diffK = (lngKbisulf +LOG(mratio)) -lnKHSO4atT
        ELSE
            mratio = deps
            diffK = (lngKbisulf +LOG(mratio)) -lnKHSO4atT
        ENDIF
    ENDIF

    END SUBROUTINE DiffKsulfuricDissoc 
    !================================================================================================= 
    
    
    !=================================================================================================
    !internal wrapper function for the computation of Diffk as a function of mHSO4 input.
    MODULE FUNCTION fHSO4dissoc(mHSO4in)       
    
    IMPLICIT NONE
    !interface:
    REAL(8), INTENT(INOUT) :: mHSO4in
    REAL(8) :: fHSO4dissoc              !the function return value
    !......................................
    !calculate the actual value of the dissociation (non)-equilibrium:
    CALL DiffKsulfuricDissoc(mHSO4in, fHSO4dissoc)
    
    END FUNCTION fHSO4dissoc
    !=================================================================================================
    
        
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   The degree of dissociation of bicarbonate(HCO3-) and bisulfate(HSO4-) ions will    *
    !*   be determined iteratively for the given mixture and temperature.                   *
    !*                                                                                      *
    !*   This subroutine makes use of hybrd method for finding n roots of n equations.      *
    !*   This third-party code is distributed under the GNU LGPL license (see minpack.f90). *
    !*                                                                                      *
    !*   :: Authors & Copyright ::                                                          *
    !*   Hang Yin, Andi Zuend                                                               *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 2021 May                                                               *
    !*   -> latest changes: 2021-12-03                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    MODULE SUBROUTINE HSO4_and_HCO3_dissociation()
    
    USE Mod_MINPACK, ONLY : hybrd1
    
    IMPLICIT NONE
    !Local variables and parameters:
    INTEGER(4) :: k, n, nrd, info, iflag, iloop, maxiloop
    !real parameters:    
    REAL(8),PARAMETER :: T0 = 298.15D0              ![k] reference temperature for dissociation constant temperature dependency calculation
    REAL(8),PARAMETER :: fa1 = -8.204327D+02, fa2 = -1.4027266D-01, fa3 = 5.027549D+04, fa4 = 1.268339D+02, fa5 = -3.879660D+06 !terms/factors for the temperature dependent bicarbonate1 dissociation constant calculation.
    REAL(8),PARAMETER :: fb1 = -2.484192D+02, fb2 = -7.489962D-02, fb3 = 1.186243D+04, fb4 = 3.892561D+01, fb5 = -1.297999D+06  !terms/factors for the temperature dependent bicarbonate2 dissociation constant calculation.
    REAL(8),PARAMETER :: fc1 = 2.495691D+02, fc2 = 4.570806D-02, fc3 = -1.593281D+04, fc4 = -4.045154D+01, fc5 = 1.541270D+6    !terms/factors for the temperature dependent CO2 Henry's law constant
    REAL(8),PARAMETER :: fd1 = -1.4816780D+02, fd2 = 8.933802D-01, fd3 = -2.332199D-03, fd4 = 2.146860D-06
    !--
    REAL(8) :: alphaHSO4pred, dperturb, mfHSO4maxinions, nCa, nCaSO4, r, sumCat, T, tol, wfHSO4max
    REAL(8) :: ln_nCarbmax, ln_nCO2max, ln_nHCO3max, ln_nHmax, ln_nHSO4max, ln_nOHmax, wtf_water_init
    REAL(8),DIMENSION(8) :: molarinp
    REAL(8),DIMENSION(:),ALLOCATABLE :: diffK, randomval
    !--
    LOGICAL(4),PARAMETER :: memory_initial_guess = .true.   !if true, keeps determined solver values to be used for subsequent initial guess in calc.
    LOGICAL(4),PARAMETER :: HighPrecSolving = .true.        !switch to set mode of dissociation calculation between .false. (=faster) and high-precision solving (slower);
    LOGICAL(4) :: calcviscosity_saved
    !........................................................................
    
    !settings for calculations with CO2 equilibration to gas phase of constant mixing ratio.
    use_CO2gas_equil = .false.      !switch for use of pseudo CO2(g) <--> CO2(aq) equilibration to the "target_CO2gas_ppm" alongside the aqueous phase equilibria;
    target_CO2gas_ppm = 400.0D0     !400 ppm by volume CO2(g) in 1 m^3 of air at ptot = 1 atm;
    
    !Initial molar amounts of ions + neutrals in the mixture for 1 kg of solvent at input => equivalent to molality (initially)
    wtf_water_init = wtf(1)         !save input mass frac. of water (since that may change when water self-dissociation is accounted for)
    ALLOCATE( mNeutral(nneutral), molNeutral(nneutral), SNC(NGI), SNA(NGI) )
    molNeutral = 0.0D0
    SNA = SMA
    SNC = SMC
    CALL MassFrac2SolvMolalities(wtf, mNeutral)
    molNeutral = mNeutral
    
    sumCat = SUM(SNC)
    nH2O_init = mNeutral(1)         !store the initial amount of H2O for nH2O change
    nCO2_init = mNeutral(idCO2)     !store the initial amount of CO2(aq) for nH2O change
    ntiny = deps                    !minimal acceptable molar amount of any dissociation reagent (set to a small scaling value).
    
    !As water is involved in the reaction, the calculations should be molar amounts instead of molalities  
    nHCO3 = SNA(idHCO3)             !initial molar amount of [HCO3-]
    nCarb = SNA(idCO3)              !initial molar amount of [CO3--]
    nOH = SNA(idOH)                 !initial molar amount of [OH-]
    nCO2 = nCO2_init                !initial molar amount of [CO2]
    nH = SNC(idH)                   !initial molar amount of [H+]
    nOH_init = nOH
    IF (bisulfsyst) THEN
        nHSO4 = SNA(idHSO4)         !initial molar amount of [HSO4-]
        nSulf = SNA(idSO4)          !initial molar amount of [SO4--]
    ELSE
        nHSO4 = 0.0D0
        nSulf = 0.0D0
    ENDIF
    !for tests:  make sure that nOHmax is at least 3.5*sumCat value, if possible:
    r = 3.5D0*sumCat/nH2O_init
    r = MIN(r, 0.98D0)
    nOHmax = r*nH2O_init + nOH
    
    !calculate maximum possible molar amounts: 
    IF (use_CO2gas_equil) THEN
        n = 7
        nCO2gas = target_CO2gas_ppm*1.0D-6/(T_K*Rgas_atm)                           !initial molar amount of [CO2(g)]; 
        nHmax = nH + nHCO3 + nCO2 +nCO2 + nCO2gas +nCO2gas + (nOHmax -nOH) + nHSO4  !H+ from 2*min(nCO2,nH2O), usually nH2O is more abundant
        nCarbmax = nCarb + nHCO3 + nCO2 + nCO2gas
    ELSE
        n = 5
        nCO2gas = 0.0D0
        nHmax = nH + nHCO3 + nCO2 + nCO2 + (nOHmax -nOH) + nHSO4
        nCarbmax = nCarb + nHCO3 + nCO2
    ENDIF
    IF (bisulfsyst) THEN
        n = n + 1
    ENDIF
    nSulfmax = nHSO4 + nSulf
    nHSO4max = MIN(nHmax, nSulfmax)
    nHCO3max = MIN(nHmax, nCarbmax)
    nCO2max = MIN(0.5D0*nHmax, nCarbmax)
    nVmax = 1.0D9
    
    !n: the number of unknowns / equations to solve:
    ALLOCATE(solve_var(n), solve_var_maxval(n), map_solver2molarinp(n), randomval(n), diffK(n))
    map_solver2molarinp = 0                          
    IF (.NOT. ALLOCATED(solve_var_saved)) THEN
        ALLOCATE(solve_var_saved(n))                   !allocate submodule variable on first use, but not subsequently
        solve_var_saved = -9.999D280                   !initialized with an unfeasible value
    ENDIF

    !if Ca++ is part of the input, remove all the CaSO4(s) that could form;
    IF (idCa > 0) THEN
        nCa = SNC(idCa)
        nCaSO4 = MIN(nCa, nSulfmax)
        IF (nCa < nSulfmax) THEN
            SNC(idCa) = 0.0D0
            nSulf = nSulf + nHSO4 -nCa
            nH = nH + nHSO4
            nHSO4 = 0.0D0
            !WRITE(*,*) "WARNING: All the CaSO4 was precipitated out by default."
        ELSE
            !leave tiny amount (1000.0D0*deps) of SO4-- in the system  0.001% of Sulfmax or 1000.0D0*deps
            SNC(idCa) = nCa -nSulfmax + MIN(1.0D-5*nSulfmax, 1.0D3*deps)
            nSulf = MIN(1.0D-5*nSulfmax, 1.0D3*deps)
            nH = nH + nHSO4
            nHSO4 = 0.0D0            
            !WRITE(*,*) "WARNING: All the CaSO4 was precipitated out by default"
        ENDIF
        nSulfmax = nHSO4 + nSulf
        nHmax = nH + nHCO3 + nCO2+nCO2 + (nOHmax -nOH) + nHSO4
        nHSO4max = MIN(nHmax, nSulfmax)
        nHCO3max = MIN(nHmax, nCarbmax)
        nCO2max = MIN(0.5D0*nHmax, nCarbmax)
    ENDIF
        
    !if no H+ or no carbonate or no carbonate acid are in the solution at the same time, then step over, otherwise calculate dissociation:
    IF (nHCO3max > ntiny) THEN !@@**## there is some dissociation, so proceed...
        calcviscosity_saved = calcviscosity
        calcviscosity = .false.
        ln_nHCO3max = LOG(nHCO3max*(1.0D0 -lztiny))
        ln_nCarbmax = LOG(nCarbmax*(1.0D0 -lztiny))
        ln_nCO2max = LOG(nCO2max*(1.0D0 -lztiny))
        ln_nOHmax = LOG(nOHmax*(1.0D0 -lztiny))
        ln_nHmax = LOG(nHmax*(1.0D0 -lztiny))
        IF (bisulfsyst) THEN
            ln_nHSO4max = LOG(nHSO4max*(1.0D0 -lztiny))
        ELSE
            ln_nHSO4max = 0.0D0
        ENDIF
        ln_maxval_inp = [ ln_nHCO3max, ln_nCarbmax, ln_nCO2max, ln_nOHmax, ln_nHSO4max, &
                        & ln_nHmax, ln_nCO2max, LOG(nVmax*(1.0D0 -lztiny)) ]  
        
        !compute the known equilibrium constants for the dissociation equilibria at T_K:
        T = T_K
        lnK1HCO3atT = fa1 + fa2*T + fa3/T + fa4*LOG(T) + fa5/(T**2)
        lnK2HCO3atT = fb1 + fb2*T + fb3/T + fb4*LOG(T) + fb5/(T**2)
        lnKCO2atT = fc1 + fc2*T + fc3/T + fc4*LOG(T) + fc5/(T**2)
        lnKH2OatT = fd1 + fd2*T + fd3*(T**2) + fd4*(T**3)   !Valid from 273K to 323 k.
        IF (bisulfsyst) THEN
            !calculate the reference value of the 2nd dissociation constant of H2SO4 as a function of temperature 
            !using parameterization by Knopf et al. (2003) (corrected equation 16). 
            !The parameterization is valid in the T-range [179 k < T_K < 474 k], for lower temperatures, 
            !replace T with 180 k temporarily and use the ln(KHSO4_T) 
            !of that temperature as that might be the best estimate for that case.
            T = T_K
            IF (T < 180.0D0) THEN !use T of 180 k for parameterization
                T = 180.0D0
            ELSE IF (T > 473.0D0) THEN !use T of 473 k
                T = 473.0D0
            ENDIF
            lnKHSO4atT = -4.54916799587D0 -fe1*(fe2*(1.0D0/T-fe3)-fe4*LOG(T/T0) -0.5D0*dcpdT*(T-T0))
        ENDIF
        
        !##! section to compute a good initial guess for the msulf1, msulf2 range based on a parameterization of alphaHSO4 as a function of mass frac. of HSO4(max) in the water portion of the solvent mixture.
        IF (bisulfsyst) THEN
            !(A) calculate the special normalized mass fraction of HSO4max in water plus HSO4max: wfHSO4max; only the water fractional amount associated with HSO4max is considered (excluding water associated with other ions in the system).
            mfHSO4maxinions = 2.0D0*nSulfmax/SumIonMolalities !mole fraction of (HSO4max plus counter cations) relative to all ion molar amounts.
            wfHSO4max = nSulfmax*MolarMassHSO4/(nSulfmax*MolarMassHSO4 + mfHSO4maxinions*wtf(1)/SUM(wtf(1:nneutral)))
            !(B) calculate the predicted degree of dissociation at this wfHSO4max value based on the parameterization:
            alphaHSO4pred = 1.0D0 - 1.0D0/(1.0D0 + (1.0D0/wfHSO4max**c1 - 1.0D0/wfHSO4max**c2))**c4 &
                & + c3*wfHSO4max**c5*(1.0D0-wfHSO4max)**1.75D0
            alphaHSO4pred = MIN(alphaHSO4pred, 1.0D0-deps) !prevent values larger than 1.0D0
            alphaHSO4pred = MAX(alphaHSO4pred, deps)
        ENDIF

        CALL RANDOM_SEED (SIZE = nrd)
        CALL RANDOM_SEED(PUT=[(k*211839831, k=1,nrd)])     !initialize the pseudo-random number generator each time with the same seed (for debugging reproducability)
        IF (HighPrecSolving) THEN
            maxiloop = 40
            tol = 1.0D-2*sqrtdeps
        ELSE
            maxiloop = 9
            tol = sqrtdeps
        ENDIF
        !---------------------
        DO iloop = 1,maxiloop
            info = -1
            IF (iloop == 1) THEN
                !set the index mapping between molarinp entries and the array of solver variables via map_solver2molarinp:
                IF (bisulfsyst .AND. use_CO2gas_equil) THEN
                    map_solver2molarinp(1:n) = [1, 2, 3, 4, 5, 6, 7, 8]
                ELSE IF (bisulfsyst) THEN
                    map_solver2molarinp(1:n) = [1, 2, 3, 4, 5, 6]
                ELSE IF (use_CO2gas_equil) THEN
                    map_solver2molarinp(1:n) = [1, 2, 3, 4, 6, 7, 8]    !skip 5 in this case
                ELSE
                    map_solver2molarinp(1:n) = [1, 2, 3, 4, 6]          !skip 5
                ENDIF
                solve_var_maxval(1:n) = ln_maxval_inp(map_solver2molarinp(1:n))
                !initial guess
                IF (memory_initial_guess .AND. SUM(solve_var_saved(1:n)) > -9.9D270) THEN 
                    !use saved solution from prior succesful calculation (previous data point input) as initial guess:
                    solve_var(1:n) = solve_var_saved(1:n)
                    CALL rboundsCheck( solve_var(1:n), minlim, ntiny, solve_var_maxval(1:n) )
                ELSE
                    !use initial guess values; estimated from trial calculations distinguishing with/without bisulfate:
                    IF (bisulfsyst) THEN
                        molarinp(1) = 1.60D-6*nHCO3max      !nHCO3
                        molarinp(2) = 3.0D0*deps*nCarbmax   !nCarb
                        IF (use_CO2gas_equil) THEN
                            molarinp(3) = 2.0D-5*nCO2max    !nCO2(aq)
                        ELSE
                            molarinp(3) = 0.9999D0*nCO2max  !nCO2(aq)
                        ENDIF
                        molarinp(4) = 30.0D0*deps*nOHmax    !nOH
                        nSulf = alphaHSO4pred*nSulfmax      !nSulf
                        molarinp(5) = nSulfmax -nSulf       !nHSO4
                        molarinp(6) = molarinp(1) +2.0D0*molarinp(2) +molarinp(4) +molarinp(5) +2.0D0*nSulf     !nH; start with neutral condition
                        molarinp(7) = 0.9999D0*nCO2max      !nCO2gas
                        molarinp(8) = 40.0D0                !V
                    ELSE
                        molarinp(1) = 0.13D0*nHCO3max       !nHCO3
                        molarinp(2) = 0.6D0*nCarbmax        !nCarb
                        molarinp(3) = 1.0D-5*nCO2max        !nCO2
                        molarinp(4) = 2.0D-5*nOHmax         !nOH
                        nSulf = 0.0D0
                        molarinp(5) = 0.0D0                 !nHSO4
                        molarinp(6) = molarinp(1) +2.0D0*molarinp(2) +molarinp(4) +molarinp(5) +2.0D0*nSulf     !nH; start with neutral condition
                        molarinp(7) = 0.28D0*nCO2max        !nCO2gas
                        molarinp(8) = 20.0D0                !V
                    ENDIF
                    solve_var(1:n) = LOG(molarinp(map_solver2molarinp(1:n)))
                ENDIF
            ELSE IF (MOD(iloop,2) /= 0) THEN
                !perturb all molar amounts:
                CALL RANDOM_NUMBER(randomval)           !get random numbers in the interval [0.0, 1.0]
                randomval = -1.0D0 +2.0D0*randomval     !values in interval [-1.0, 1.0]
                dperturb = MIN( 0.8D0, 0.005D0 + 0.05D0**(1.0D0/REAL(iloop -1, KIND=8)) )
                solve_var = solve_var*(1.0D0-dperturb +2.0D0*randomval*dperturb)
                CALL rboundsCheck( solve_var(1:n), minlim, ntiny, solve_var_maxval(1:n) )
            ELSE
                !perturb all molar amounts slightly:
                CALL RANDOM_NUMBER(randomval)
                randomval = -1.0D0 +2.0D0*randomval
                dperturb = 2.0D-2*iloop
                solve_var = solve_var*(1.0D0-dperturb +2.0D0*randomval*dperturb)
                CALL rboundsCheck( solve_var(1:n), minlim, ntiny, solve_var_maxval(1:n) )
            ENDIF
            !--
            CALL hybrd1(DiffK_carb_sulf, n, solve_var, diffK, tol, info)
            !--
            IF (info == 1) THEN
                IF (memory_initial_guess) THEN
                    IF (SUM(ABS(diffK)) < sqrtdeps) THEN
                        solve_var_saved(1:n) = solve_var(1:n)   !save found solver values for potential use in subsequent calc.
                    ENDIF
                ENDIF
                !k = 1       !for debugging breakpoint
                EXIT
            ENDIF
        ENDDO !iloop
        !---------------------
            
        !BLOCK
        !    USE ModSystemProp, ONLY : nd, NGI, anionZ, cationZ
        !    REAL(8) :: t1, t2, t3
        !    t1 = SUM(SMA(1:NGI)*ABS(anionZ(1:NGI)))             !test the electrical charge neutrality condition in the solution
        !    t2 = SUM(SMC(1:NGI)*cationZ(1:NGI))
        !    t3 = t1 - t2
        !    IF (t3 < 1.0D6 .AND. ABS(t3) > 1.0D-9) THEN  !test
        !        WRITE(*,'(A, ES15.8)') "WARNING: Electrical charge neutrality condition violated! ", t3
        !        WRITE(*,'(A, I0,1X,ES15.8)') "nd, wtf(1) ", nd, wtf(1)
        !    ENDIF
        !END BLOCK
        
        !If success, check found zero and return best found nHCO3:
        calcviscosity = calcviscosity_saved
        CALL DiffK_carb_sulf(n, solve_var, diffK, iflag)
        molarinp(map_solver2molarinp(1:n)) = EXP(solve_var(1:n))
        nHCO3 = molarinp(1)
        nCarb = molarinp(2)
        nCO2 = molarinp(3)
        alphaHCO3 = 1.0D0 -nHCO3/nHCO3max       !this is a generally valid expression for alphaHCO3
        !predict the equil. partial pressure of CO2 [Pa:
        pCO2 = 1.013250D+05 * mNeutral(idCO2) * EXP( MIN(gnmrln(idCO2), logval_threshold) - lnKCO2atT )
        IF (bisulfsyst) THEN
            nHSO4 = molarinp(5)
            alphaHSO4 = 1.0D0 -nHSO4/nHSO4max       !this is a generally valid expression for alphaHSO4
        ELSE
            alphaHSO4 = -9.999999D0
        ENDIF
        !checks for debugging:
        IF (iloop > 9) THEN
            k = 88 !for breakpoint
            r = SUM(ABS(diffK))
            IF (r > sqrtdeps) THEN
                IF (errorflagcalc == 0) THEN    !report back issues in solving dissociation equilibrium
                    errorflagcalc = 17
                    IF (.NOT. frominpfile) THEN
                        !$OMP CRITICAL
                        WRITE(*,'(A)') "AIOMFAC ERROR 17: Issue with ion dissociation equilibria calculations."
                        WRITE(*,'(A)') "The numerical solution of electrolyte/ion dissociation equilibria was &
                            &not accomplished to the desired tolerance level. Model output for this data &
                            &point is unreliable and likely incorrect."
                        WRITE(*,'(A,I0,ES13.6)') "iloop, wtf_water_init = ", iloop, wtf_water_init
                        WRITE(*,'(A,ES13.6)') "SUM(ABS(diffK)) = ", r
                        WRITE(*,*)
                        !$OMP END CRITICAL
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
    ELSE    !case where HCO3 is part of the system components, but at present input composition, nHCO3max is at zero molar amount.
        alphaHCO3 = -9.999999D0
        alphaHSO4 = -9.999999D0
        pCO2 = 0.0D0
        CALL Gammas() !call Gammas once
    ENDIF !@@**##
    
    DEALLOCATE(solve_var, solve_var_maxval, map_solver2molarinp, randomval, diffK)
    DEALLOCATE(mNeutral, molNeutral, SNC, SNA)
    
    END SUBROUTINE HSO4_and_HCO3_dissociation
    !=================================================================================================
    
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to calculate the difference between the reference value (lnK1HCO3atT,   *
    !*   lnK2HCO3atT, lnKHSO4atT, and the estimated value of the dissociation constant of   *
    !*   CO2(aq) + H2O <-> H+ + HCO3-, HCO3- <-> H+ + CO3--, and HSO4- <-> H+ + SO4--       *
    !*                                                                                      *
    !*   :: Authors & Copyright ::                                                          *
    !*   Hang Yin, Andi Zuend                                                               *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 2021 May                                                               *
    !*   -> latest changes: 2021-12-04                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    SUBROUTINE DiffK_carb_sulf(n, solve_var, diffK, iflag)
    !subroutine containing the system of n equations with n unknowns to be solved using Minpack methods.

    !Public Variables:
    USE ModSystemProp, ONLY : nneutral, idH, idHCO3, idCO3, idOH, idCO2, idSO4, idHSO4, NGI, topsubno, &
        & waterpresent, Mmass, Ianion, Ication
    USE ModSubgroupProp, ONLY : SMWC, SMWA

    IMPLICIT NONE

    INTEGER(4),INTENT(IN) :: n
    REAL(8),DIMENSION(n),INTENT(INOUT) :: solve_var
    REAL(8),DIMENSION(n),INTENT(OUT) :: diffK
    INTEGER(4),INTENT(OUT) :: iflag
    !local vars:
    INTEGER(4) :: I, J, k
    REAL(8) :: deltaCO2, diffK1, diffK2, diffK3, diffK4, lngammaCarb, lngammaCO2, lngammaH, lngammaHydCarb, &
        & lngammaHydSulf, lngammaOH, lngammaSulf, lngammawv, mCarbcheck, mgeomean, mHCO3check, mOHcheck,    &
        & mratio1, mratio2, mratio3, mratio4, mSulfcheck, solvmass, totmass, xwater
    REAL(8),DIMENSION(n) :: exp_solve_var
    REAL(8),DIMENSION(nneutral) :: xneutral
    !................................................
        
    !make sure solver_var are consistent with upper limits given by ln_maxval_inp:
    iflag = 0
    CALL rboundsCheck( solve_var(1:n), minlim, ntiny, solve_var_maxval(1:n) )
    exp_solve_var = EXP(solve_var)
    nHCO3 = exp_solve_var(1)
    nCarb = exp_solve_var(2)
    nCO2 = exp_solve_var(3)
    nOH = exp_solve_var(4)
    IF (bisulfsyst .AND. use_CO2gas_equil) THEN
        nHSO4 = exp_solve_var(5)
        nH = exp_solve_var(6)
        nCO2gas = exp_solve_var(7)
        V = exp_solve_var(8)
    ELSE IF (bisulfsyst) THEN
        nHSO4 = exp_solve_var(5)
        nH = exp_solve_var(6)
    ELSE IF (use_CO2gas_equil) THEN
        nHSO4 = 0.0D0
        nH = exp_solve_var(5)
        nCO2gas = exp_solve_var(6)
        V = exp_solve_var(7)
    ELSE
        nHSO4 = 0.0D0
        nH = exp_solve_var(5)
    ENDIF
    nSulf = nSulfmax -nHSO4
    
    !update the moles of all species:
    deltaCO2 = nCO2 - nCO2_init
    molNeutral(1) = nH2O_init + deltaCO2 - (nOH - nOH_init)
    IF (molNeutral(1) < ntiny) THEN     !correction necessary
        molNeutral(1) = ntiny
        deltaCO2 = molNeutral(1) -nH2O_init + (nOH - nOH_init)
        nCO2 = deltaCO2 + nCO2_init
    ENDIF
    molNeutral(idCO2) = nCO2   !CO2(aq)
    SNA(idHCO3) = nHCO3
    SNA(idCO3) = nCarb
    SNA(idOH) = nOH
    SNC(idH) = nH
    IF (bisulfsyst) THEN
        SNA(idHSO4) = nHSO4
        SNA(idSO4) = nSulf
    ENDIF
        
    !update the molalities of all ions & solvents:
    CALL Moles2solvmass(molNeutral, solvmass)   
        
    totmass = solvmass
    DO I = 1,NGI
        k = ICation(I) - 200
        J = IAnion(I) - 240
        IF (k > 0) THEN
            totmass = totmass + SNC(I)*SMWC(k)*1.0D-3
        ENDIF 
        IF (J > 0) THEN 
            totmass = totmass + SNA(I)*SMWA(J)*1.0D-3 
        ENDIF
    ENDDO
    !update the mass fraction and molalities of neutrals:
    wtf(1:nneutral) = molNeutral*Mmass(1:nneutral)/totmass
    mNeutral = molNeutral/solvmass
    SMA = SNA/solvmass
    SMC = SNC/solvmass
    mH = SMC(idH)
    mHCO3 = SMA(idHCO3)
    mCarb = SMA(idCO3)
    mOH = SMA(idOH)
    mCO2 = mNeutral(idCO2)
    IF (bisulfsyst) THEN
        mHSO4 = SMA(idHSO4)
        mSulf = SMA(idSO4)
    ENDIF

    !Mole fraction of water
    CALL Molality2SolvMoleFrac(SMA, SMC, mNeutral, xneutral)
    xwater = xneutral(1)
    !--
    CALL Gammas()   !compute the natural log of activity coefficients from SR, MR, & LR interaction parts.
    !--
    !Use the mole fraction-based ion activity coeff. (or their natural log).
    IF (idH > 0) THEN
        lngammaH = -Tmolal +gcsrln(idH) +gcmrln(idH) +gclrln(idH)
    ELSE
        lngammaH = 0.0D0
    ENDIF
    IF (idHCO3 > 0) THEN
        lngammaHydCarb = -Tmolal +gasrln(idHCO3) +gamrln(idHCO3) +galrln(idHCO3)
    ELSE
        lngammaHydCarb = 0.0D0
    ENDIF
    IF (idCO3 > 0) THEN
        lngammaCarb = -Tmolal +gasrln(idCO3) +gamrln(idCO3) +galrln(idCO3)
    ELSE
        lngammaCarb = 0.0D0
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
    IF (idOH > 0) THEN
        lngammaOH = -Tmolal +gasrln(idOH) +gamrln(idOH) +galrln(idOH)
    ELSE
        lngammaOH = 0.0D0
    ENDIF
    IF (idCO2 > 0) THEN
        lngammaCO2 = gnmrln(idCO2)
    ELSE
        lngammaCO2 = 0.0D0
    ENDIF
    IF (waterpresent) THEN
        lngammawv = gnsrln(1) +gnmrln(1) +gnlrln(1)
    ELSE
        lngammawv = 0.0D0
    ENDIF
        
    !estimate the molalities of CO3--, HCO3- and H+ to get agreement
    !with the reference value of the dissociation constant.
    lngK1bicarb = lngammaH + lngammaHydCarb - lngammaCO2 - lngammawv
    lngK2bicarb = lngammaH + lngammaCarb - lngammaHydCarb
    lngbisulf = lngammaH + lngammaSulf - lngammaHydSulf
    lngKOH = lngammaH + lngammaOH - lngammawv

    !report resulting deviation vector diffK:
    diffK = 0.0D0
    IF (mCO2*xwater > 0.0D0) THEN
        mratio1 = mH*mHCO3/(mCO2*xwater)
    ELSE
        mratio1 = 0.0D0
    ENDIF
    IF (mHCO3 > 0.0D0) THEN
        mratio2 = mH*mCarb/mHCO3
    ELSE
        mratio2 = 0.0D0
    ENDIF
    IF (xwater > 0.0D0) THEN
        mratio3 = mH*mOH/xwater
    ELSE
        mratio3 = 0.0D0
    ENDIF
    IF (bisulfsyst) THEN
        IF (mHSO4 > 0.0D0) THEN
            mratio4 = mH*mSulf/mHSO4
        ELSE
            mratio4 = 0.0D0
        ENDIF
    ENDIF
    
    IF (mHCO3 > lztiny) THEN !normal case
        IF (mratio1 > 0.0D0) THEN
            diffK(1) = (lngK1bicarb +LOG(mratio1)) -lnK1HCO3atT
        ELSE
            mratio1 = ztiny*mHCO3/(mCO2*xwater)
            diffK(1) = (lngK1bicarb +LOG(mratio1)) -lnK1HCO3atT
        ENDIF
    ELSE !very low mHCO3, requires special treatment
        IF (mratio1 > 0.0D0) THEN
            !compute the correct value of mHCO3 needed to fulfill this equation 1:
            diffK1 = lnK1HCO3atT - lngK1bicarb
            diffK1 = MAX(diffK1, -logval_threshold)
            diffK1 = MIN(diffK1, logval_threshold)
            mHCO3check = (mCO2*xwater)/mH*EXP(diffK1)
            IF (mHCO3check*solvmass < 100.0D0*ntiny) THEN   !assume mHCO3check is essentially zero, so no equation to solve;
                diffK(1) = 0.0D0  !set for this exception
                nHCO3 = solvmass*mHCO3check
            ELSE
                mgeomean = SQRT(mHCO3*mHCO3check)
                mratio1 = mH*mgeomean/(mCO2*xwater)
                diffK(1) = (lngK1bicarb +LOG(mratio1)) -lnK1HCO3atT
                nHCO3 = solvmass*mgeomean
            ENDIF
            nHCO3 = MIN(nHCO3, nHCO3max)
        ENDIF
    ENDIF !mHCO3
    IF (mCarb > lztiny) THEN  !normal case
        IF (mratio2 > 0.0D0) THEN
            diffK(2) = (lngK2bicarb +LOG(mratio2)) -lnK2HCO3atT
        ELSE
            mratio2 = ztiny/mHCO3
            diffK(2) = (lngK2bicarb +LOG(mratio2)) -lnK2HCO3atT
        ENDIF
    ELSE !very low mCarb, requires special treatment
        IF (mratio2 > 0.0D0) THEN
            !compute the correct value of mCarb needed to fulfill this equation 2:
            diffK2 = lnK2HCO3atT - lngK2bicarb
            diffK2 = MAX(diffK2, -logval_threshold)
            diffK2 = MIN(diffK2, logval_threshold)
            mCarbcheck = mHCO3/mH*EXP(diffk2)
            IF (mCarbcheck*solvmass < 100.0D0*ntiny) THEN !assume mCarb is essentially zero, so no equation to solve;
                diffK(2) = 0.0D0  !set for this exception
                nCarb = solvmass*mCarbcheck
            ELSE
                mgeomean = SQRT(mCarb*mCarbcheck)
                mratio2 = mH*mgeomean/mHCO3
                diffK(2) = (lngK2bicarb +LOG(mratio2)) -lnK2HCO3atT
                nCarb = solvmass*mgeomean
            ENDIF
            nCarb = MIN(nCarb, nCarbmax)
        ELSE
            mratio2 = ztiny/mHCO3
            diffK(2) = (lngK2bicarb +LOG(mratio2)) -lnK2HCO3atT
        ENDIF
    ENDIF !nCarb
    IF (mOH > lztiny) THEN  !normal case
        IF (mratio3 > 0.0D0) THEN
            diffK(3) = (lngKOH +LOG(mratio3)) -lnKH2OatT
        ELSE
            mratio3 = mH*ztiny/xwater
            diffK(3) = (lngKOH +LOG(mratio3)) -lnKH2OatT
        ENDIF
    ELSE !very low mOH; requires special treatment
        !compute the correct value of mOH needed to fulfill equilibrium equation 3:
        IF (mratio3 > 0.0D0) THEN
            diffK3 = lnKH2OatT - lngKOH
            diffK3 = MAX(diffK3, -logval_threshold)
            diffK3 = MIN(diffK3, logval_threshold)
            mOHcheck = xwater/mH*EXP(diffK3)
            IF (mOHcheck*solvmass < 100.0D0*ntiny) THEN !assume mCarb is essentially zero, so no equation to solve;
                diffK(3) = 0.0D0
                nOH = solvmass*mOHcheck
            ELSE
                mgeomean = SQRT(mOH*mOHcheck)
                mratio3 = mH*mgeomean/xwater
                diffK(3) = (lngKOH +LOG(mratio3)) -lnKH2OatT
                nOH = solvmass*mgeomean
            ENDIF
            nOH = MIN(nOH, nOHmax)
        ENDIF
    ENDIF !mOH
    IF (bisulfsyst) THEN
        IF (mSulf > lztiny .AND. mHSO4 > lztiny) THEN
            IF (mratio4 > 0.0D0) THEN
                diffK(4) = (lngbisulf +LOG(mratio4)) -lnKHSO4atT
            ELSE
                mratio4 = mH*ztiny/xwater
                diffK(4) = (lngbisulf +LOG(mratio4)) -lnKHSO4atT
            ENDIF
        ELSE    !very low mSulf, requires special treatment
            !compute the "correct" value of mSulf needed to fulfill equilibrium equation 4:
            IF (mratio4 > 0.0D0) THEN
                diffK4 = lnKHSO4atT - lngbisulf
                diffK4 = MAX(diffK4, -logval_threshold)
                diffK4 = MIN(diffK4, logval_threshold)
                mSulfcheck = mHSO4/mH*EXP(diffk4)
                IF (mSulfcheck*solvmass < 100.0D0*ntiny) THEN
                    diffK(4) = 0.0D0
                    nSulf = solvmass*mSulfcheck
                ELSE
                    mgeomean = SQRT(mSulf*mSulfcheck)
                    mratio4 = mH*mgeomean/mHSO4
                    diffK(4) = (lngbisulf +LOG(mratio4)) -lnKHSO4atT
                    nSulf = solvmass*mgeomean
                ENDIF
                nSulf = MIN(nSulf, nSulfmax -2.0D0*deps*nSulfmax)
                nHSO4 = nSulfmax -nSulf
            ENDIF
        ENDIF !mSulf
    ENDIF

    IF (bisulfsyst .AND. use_CO2gas_equil) THEN
        !write(*,'(8(ES13.6,1X))') [nHCO3, nCarb, nCO2, nOH, nHSO4, nH, nCO2gas, V]
        solve_var = LOG([nHCO3, nCarb, nCO2, nOH, nHSO4, nH, nCO2gas, V])
        diffK(5) = mCO2*EXP(MIN(lngammaCO2, logval_threshold))*2.5D3 - EXP(lnKCO2atT)
        diffK(6) = 1.0D2*(sum_sorted([nHmax, -nHCO3, -nH, -nCO2, -nCO2, -nCO2gas, -nCO2gas, -nHSO4, -nOHmax, nOH]))
        diffK(7) = 1.0D2*(sum_sorted([nCarbmax, -nHCO3, -nCarb, -nCO2, -nCO2gas]))                   !molar balance
        diffK(8) = target_CO2gas_ppm*1.0D-6*V -nCO2gas*Rgas_atm*T_K
    ELSE IF (bisulfsyst) THEN
        solve_var = LOG([nHCO3, nCarb, nCO2, nOH, nHSO4, nH])
        diffK(5) = 1.0D2*(sum_sorted([nCarbmax, -nHCO3, -nCarb, -nCO2]))                            !nCarbmax -nHCO3 -nCarb -nCO2 !mass balance
        diffK(6) = 1.0D2*(sum_sorted([nHmax, -nHCO3, -nH, -nCO2, -nCO2, -nHSO4, -nOHmax, nOH]))     !nHmax -nHCO3 -nH -2.0D0*nCO2 -nHSO4 -nOHmax + nOH !mass balance
    ELSE IF (use_CO2gas_equil) THEN
        solve_var = LOG([nHCO3, nCarb, nCO2, nOH, nH, nCO2gas, V])
        diffK(4) = mCO2*EXP(MIN(lngammaCO2, logval_threshold))*2.5D3 - EXP(lnKCO2atT)
        diffK(5) = 1.0D2*(sum_sorted([nHmax, -nHCO3, -nH, -nCO2, -nCO2, -nCO2gas, -nCO2gas, -nOHmax, nOH]))
        diffK(6) = 1.0D2*(sum_sorted([nCarbmax, -nHCO3, -nCarb, -nCO2, -nCO2gas]))                   !molar balance
        diffK(7) = target_CO2gas_ppm*1.0D-6*V -nCO2gas*Rgas_atm*T_K
    ELSE
        solve_var = LOG([nHCO3, nCarb, nCO2, nOH, nH])
        diffK(4) = 1.0D2*(sum_sorted([nHmax, -nHCO3, -nH, -nCO2, -nCO2, -nOHmax, nOH]))
        diffK(5) = 1.0D2*(sum_sorted([nCarbmax, -nHCO3, -nCarb, -nCO2]))
    ENDIF

    END SUBROUTINE DiffK_carb_sulf
    !=================================================================================================
    
    
    !Utility subroutine to check, and if necessary adjust, values of variables limited to upper/lower bounds.
    !The upper bound is given by a slightly scaled maxlim value.
    !--------------------------------------------------------------------------------------------------
    PURE ELEMENTAL MODULE SUBROUTINE rboundsCheck(rset, lowerbound, ntiny, maxlim)
    !
    IMPLICIT NONE
    !interface arguments:
    REAL(8),INTENT(INOUT) :: rset
    REAL(8),INTENT(IN) :: lowerbound, ntiny, maxlim
    !local parameters:
    REAL(8),PARAMETER :: scaleval = 1.0D-3, deps = EPSILON(1.0D0), threedeps = 3.0D0*deps
    REAL(8) :: qa, qacorr, qaband, higherbound
    !.......................

    !use smart adjustments on values with upper/lower bounds if necessary:
    higherbound = (1.0D0 -ntiny)*maxlim
    qaband = 0.05D0*ABS(higherbound - lowerbound)
    IF (qaband < threedeps) THEN
        qaband = MAX(0.2D0*ABS(higherbound - lowerbound), deps)
    ENDIF

    IF (rset < lowerbound) THEN
        qacorr = (lowerbound-rset)*scaleval
        qa = lowerbound + qacorr
        DO WHILE (qa > (lowerbound +qaband)) !(qa > 1.0D0)
            qacorr = ABS(qacorr)*0.1D0
            qa = lowerbound + qacorr
        ENDDO
        rset = qa
    ELSE IF (rset > higherbound) THEN
        qacorr = (rset -higherbound)*scaleval
        qa = higherbound -qacorr
        DO WHILE (qa < (higherbound -qaband)) !(qa < 0.0D0)
            qacorr = ABS(qacorr)*0.1D0
            qa = higherbound -qacorr
        ENDDO
        rset = qa
    ENDIF

    END SUBROUTINE rboundsCheck
    !--------------------------------------------------------------------------------------------------
    
    
    !--------------------------------------------------------------------------------------------------
    !Utility function:
    !** A safer way for summing a "small" list of values that potentially differ substantially in magnitude ** 
    !** based on sorting them and summing the elements starting with the smallest absolute value.           **
    PURE MODULE FUNCTION sum_sorted(list) RESULT(summed)
        
    IMPLICIT NONE
    !interface arguments:
    REAL(8),DIMENSION(:),INTENT(IN) :: list
    REAL(8) :: summed
    !local arguments
    INTEGER(4) :: i, k
    LOGICAL(4),DIMENSION(SIZE(list)) :: available
    !...........................
        
    available = .true.
    !sort the list from smallest to largest in magnitude (ignoring neg/pos signs):
    summed = 0.0D0
    DO i = 1,SIZE(list)
        k = MINLOC(ABS(list(:)), MASK = available(:), DIM=1)   
        available(k) = .false.
        summed = summed + list(k) 
    ENDDO
        
    END FUNCTION sum_sorted
    !--------------------------------------------------------------------------------------------------
    

END SUBMODULE SubModDissociationEquil