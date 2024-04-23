!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   submodule of module ModCalcActCoeff containing procedures used for the calculation *
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
!*   -> latest changes: 2023-09-11                                                      *
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
!*   :: List of subroutines and functions contained in this submodule:                  *
!*   --------------------------------------------------------------                     *
!*   -  subroutine HSO4_dissociation                                                    *
!*   -  subroutine DiffKsulfuricDissoc                                                  *
!*   -  function   fHSO4dissoc                                                          *
!*   -  subroutine HSO4_and_HCO3_dissociation                                           *
!*   -  subroutine DiffK_carb_sulf                                                      *
!*   -  pure elemental subroutine rboundsCheck                                          *
!*   -  pure function sum_sorted                                                        *
!*                                                                                      *
!****************************************************************************************
submodule (ModCalcActCoeff) SubModDissociationEquil

use ModSystemProp, only : bisulfsyst, calcviscosity, errorflag_clist, frominpfile, idH, idHCO3, &
    & idCO3, idOH, idCO2, idHSO4, idSO4, idCa, nneutral, NGI 
use ModAIOMFACvar, only : alphaHSO4, diffKHSO4, SMA, SMC, SumIonMolalities, T_K, wtf
use ModCompScaleConversion, only : MassFrac2SolvMolalities, MassFrac2IonMolalities, Moles2solvmass, &
    & Molality2SolvMoleFrac

implicit none
!submodule parameters and variables:
integer,dimension(:),allocatable :: map_solver2molarinp             !array storing the indices for mapping between the active set of solver variables and the "molarinp" variables;
!..
real(wp),parameter :: deps = epsilon(1.0_wp), sqrtdeps = sqrt(deps)
real(wp),parameter :: ztiny = tiny(1.0_wp)**0.3_wp, lztiny = 4.0_wp*deps
real(wp),parameter :: minlim = log(1.0E-6_wp*deps)    
real(wp),parameter :: Rgas_atm = 8.2057366081E-5_wp                 !Rgas in units of [m^3.atm.K^-1.mol^-1]
!..
real(wp),parameter :: T0 = 298.15_wp                                ![K] reference temperature for dissociation constant temperature dependency calculation
real(wp),parameter :: MolarMassHSO4 = 0.097071_wp                                                                           !the molar mass of the bisulfate ion in [kg/mol]
real(wp),parameter :: c1 = 2.63468E-01_wp, c2 = 2.63204E-01_wp, c3 = 1.72067_wp, c4 = 7.30039E+01_wp, c5 = 7.45718E-01_wp   !fitted parameters for the estimation of alphaHSO4pred based on scaled mass fraction of HSO4max
real(wp),parameter :: DH0 = -1.8554748598E4_wp, cp0 = -2.05594443486E2_wp, dcpdT = -9.94992864240E-1_wp                     !parameters for the temperature dependent bisulfate dissociation constant calculation.
real(wp),parameter :: fe1 = 1.0_wp/8.314472_wp, fe2 = DH0-cp0*T0+0.5_wp*dcpdT*T0*T0, fe3 = 1.0_wp/T0, fe4 = cp0-dcpdT*T0    !terms/factors for the temperature dependent bisulfate dissociation constant calculation.
!..
real(wp) :: lngbisulf, lngK1bicarb, lngK2bicarb, lngKbisulf, lngKOH, lnK1HCO3atT, lnK2HCO3atT, lnKCO2atT, lnKH2OatT, &
    & lnKHSO4atT, mHmax, mHSO4max, mHSO4min, mSulfmax, nCarbmax, nCO2_init, nCO2gas, nCO2max, nH2O_init, nHCO3max,  &
    & nHmax, nHSO4max, nOH_init, nOHmax, nSulfmax, ntiny, nVmax, target_CO2gas_ppm, V
real(wp) :: mCarb, mCO2, mH, mHCO3, mHSO4, mOH, mSulf            !molalities for activity calculations
real(wp) :: nCarb, nCO2, nH, nHCO3, nHSO4, nOH, nSulf            !molar amounts for composition determination
real(wp),dimension(8) :: ln_maxval_inp
real(wp),dimension(:),allocatable :: mNeutral, molNeutral        !molalities and moles of neutral species
real(wp),dimension(:),allocatable :: SNA, SNC                    !moles of ions
real(wp),dimension(:),allocatable :: solve_var, solve_var_maxval, solve_var_saved
logical :: use_CO2gas_equil

!$OMP THREADPRIVATE(use_CO2gas_equil, solve_var, solve_var_maxval, solve_var_saved, map_solver2molarinp, mCarb, mCO2,   &
    !$OMP & mH, mHCO3, mHSO4, mOH, mSulf, nCarb, nCO2, nH, nHCO3, nHSO4, nOH, nSulf, mNeutral, molNeutral, lngbisulf,   &
    !$OMP & lngK1bicarb, lngK2bicarb, lngKbisulf, lngKOH, lnK1HCO3atT, lnK2HCO3atT, lnKCO2atT, lnKH2OatT, lnKHSO4atT,   &
    !$OMP & ln_maxval_inp, mHmax, mHSO4max, mHSO4min, mSulfmax, nCarbmax, nCO2_init, nCO2gas, nCO2max, nH2O_init,       &
    !$OMP & nHCO3max, nHmax, nHSO4max, nOH_init, nOHmax, nSulfmax, ntiny, nVmax, SNA, SNC, target_CO2gas_ppm, V)

!======================
    contains
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
    module subroutine HSO4_dissociation()

    use ModSystemProp, only : calcviscosity, idH, idHSO4, idSO4, nneutral

    implicit none
    !Local variables and parameters:
    integer :: ntry, nb
    real(wp) :: alphaHSO4pred, b, c, diffK, dm, dm2, mfHSO4maxinions, mHSO4Guess1, mHSO4Guess2, mr, &
        & mscale, T, tol, wfHSO4max 
    real(wp),dimension(:),allocatable :: mSb1, mSb2
    real(wp),external :: brentzero   !Brent's method for finding a zero of a univariate function (external function)
    logical :: success, geomscal, calcviscosity_saved
    !................................................................................................

    !total [H+] molality: mH = 2*[H2SO4]+[HCl]+[HNO3]+[HBr] + ...; here this is the sum: [H+]+[HSO4-], 
    !because [H+] is already the sum of the molalities from all H+ containing electrolytes.
    mHmax = 0.0_wp
    mHSO4max = 0.0_wp
    mSulfmax = 0.0_wp                        !highest possible molality of SO4-- 
    mSulf = 0.0_wp
    mHSO4 = 0.0_wp
    mH = 0.0_wp
    
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
    mHSO4max = mHSO4 +min(mH, mSulf)
    mHSO4min = 0.0_wp                       !this is always the case
    mscale = 1.0E-3_wp*deps                 !minimal acceptable molality of any dissociation reagent (small scaling value).

    !If no proton or no sulfate or no sulfuric acid are in the solution at the same time, 
    !then step over, otherwise calculate dissociation:
    if (mHSO4max > mscale) then !@@**## there is some dissociation, so proceed...
        !calculate the reference value of the 2nd dissociation constant of H2SO4 as a function of temperature 
        !using parameterization by Knopf et al. (2003) (corrected equation 16). 
        !The parameterization is valid in the T-range [179 k < T_K < 474 k], for lower temperatures, 
        !replace T with 180 k temporarily and use the ln(KHSO4_T) 
        !of that temperature as that might be the best estimate for that case.
        T = T_K
        if (T < 180.0_wp) then !use T of 180 k for parameterization
            T = 180.0_wp
        else if (T > 473.0_wp) then !use T of 473 k
            T = 473.0_wp
        endif
        lnKHSO4atT = -4.54916799587_wp -fe1*(fe2*(1.0_wp/T-fe3)-fe4*log(T/T0) -0.5_wp*dcpdT*(T-T0))
        calcviscosity_saved = calcviscosity
        calcviscosity = .false.
        !..
        if (mHSO4max-mscale > 0.0_wp) then
            mHSO4max = mHSO4max-mscale
        else
            mscale = mHSO4max*5.0E-1_wp
            mHSO4max = mHSO4max-mscale
        endif
    
        !Set stoichiometric limits such that none of the ion molalities is exactly zero:
        if (mHSO4 < mscale) then
            mHSO4 = mscale
            mH = max(mH -mscale, 0.0_wp)
            mSulf = max(mSulf -mscale, 0.0_wp)
        endif
        if (mH < mscale) then
            mscale = mscale*0.5_wp
            mHSO4 = mHSO4 -mscale
            mH = mH +mscale
            mSulf = mSulf +mscale
        endif
        if (mSulf < mscale) then
            mscale = mscale*0.5_wp
            mHSO4 = mHSO4 -mscale
            mH = mH +mscale
            mSulf = mSulf +mscale
        endif
        !
        !##! section to compute a good initial guess for the msulf1, msulf2 range based on a parameterization of 
        !alphaHSO4 as a function of mass frac. of HSO4(max) in the water portion of the solvent mixture.
        !(A) calculate the special normalized mass fraction of HSO4max in water plus HSO4max: wfHSO4max; 
        !only the water fractional amount associated with HSO4max is considered (excluding water associated with other ions in the system).
        mfHSO4maxinions = 2.0_wp*mHSO4max/SumIonMolalities !mole fraction of (HSO4max plus counter cations) relative to all ion molar amounts.
        wfHSO4max = mHSO4max*MolarMassHSO4/(mHSO4max*MolarMassHSO4 + mfHSO4maxinions*wtf(1)/sum(wtf(1:nneutral)))
        !(B) calculate the predicted degree of dissociation at this wfHSO4max value based on the parameterization:
        alphaHSO4pred = 1.0_wp - 1.0_wp/(1.0_wp + (1.0_wp/wfHSO4max**c1 - 1.0_wp/wfHSO4max**c2))**c4 &
            & + c3*wfHSO4max**c5*(1.0_wp -wfHSO4max)**1.75_wp
        alphaHSO4pred = min(alphaHSO4pred, 1.0_wp -mscale)        !prevent values larger than 1.0_wp
        alphaHSO4pred = max(alphaHSO4pred, mscale)
        !(C) calculate the lower and upper guesses for msulf based on the predicted dissociation degree:
        mHSO4Guess1 = (1.0_wp - min(alphaHSO4pred +0.15_wp, 1.0_wp))*mHSO4max
        mHSO4Guess2 = (1.0_wp - max(alphaHSO4pred -0.15_wp, 0.0_wp))*mHSO4max
        !##!
        success = .false.
        !(1a) search for bracketed root:
        ntry = 1
        geomscal = .false.
        allocate(mSb1(1:ntry+1), mSb2(1:ntry+1))
        call zerobracket_inwards(fHSO4dissoc, mHSO4Guess1, mHSO4Guess2, ntry, geomscal, nb, mSb1, mSb2, success)
        if (.NOT. success) then                         !check lower bracket
            if (mHSO4min < mHSO4Guess1-mscale) then     !(1b) search in full interval from mSulfmin to mSulfmax:
                mHSO4Guess1 = mHSO4Guess1+mscale
                call zerobracket_inwards(fHSO4dissoc, mHSO4min, mHSO4Guess1, ntry, geomscal, nb, mSb1, mSb2, success)
            endif
        endif
        if ((.NOT. success) ) then
            if (mHSO4Guess2 < mHSO4max-mscale) then     !check upper bracket !(1c) search in full interval from mSulfmin to mSulfmax:
                mHSO4Guess2 = mHSO4Guess2-mscale
                call zerobracket_inwards(fHSO4dissoc, mHSO4Guess2, mHSO4max, ntry, geomscal, nb, mSb1, mSb2, success)
            endif
        endif
        if (success) then
            !(2) use Brent's method to find the root.
            mHSO4Guess1 = mSb1(1)
            mHSO4Guess2 = mSb2(1)
            tol = min( max(mHSO4Guess1, mHSO4guess2)*sqrtdeps, 1.0E-2_wp*sqrtdeps )
            mHSO4 = brentzero(mHSO4Guess1, mHSO4Guess2, deps, tol, fHSO4dissoc)     !using Brent's method with fHSO4dissoc
            !check found zero and return best found mHSO4:
            calcviscosity = calcviscosity_saved
            call DiffKsulfuricDissoc(mHSO4, diffK)
            !mr = mHSO4      !just for debugging; breakpoint
        else
            if (alphaHSO4pred > 0.4_wp) then
                mHSO4 = mHSO4min +alphaHSO4pred*mscale
                call DiffKsulfuricDissoc(mHSO4, diffK)
                dm = mSulf+mHSO4+mH
                mHSO4Guess1 = mHSO4
            else
                mHSO4Guess1 = (1.0_wp - min(alphaHSO4pred, 1.0_wp))*mHSO4max
                call DiffKsulfuricDissoc(mHSO4, diffK)
                dm = mSulf+mHSO4+mH
                mHSO4Guess1 = mHSO4
            endif
            !now solve a quadratic equation, keeping the ion activity coeff. constant 
            !to determine the exact bisulfate-system ion molalities.
            if (lngKbisulf > 300.0_wp .OR. lngKbisulf < -300.0_wp) then       !apply floating point overflow protection measures
                lngKbisulf = sign(300.0_wp, lngKbisulf)
            endif
            mr = exp(lnKHSO4atT - lngKbisulf)
            b = -(SMC(idH) +SMA(idSO4) +mr)
            c = SMC(idH)*SMA(idSO4) -mr*SMA(idHSO4)
            tol = b*b -4.0_wp*c
            if (tol >= 0.0_wp) then
                tol = sqrt(tol)
                !use a numerically stable approach to find the two roots:
                if (b >= 0.0_wp) then
                    dm = 0.5_wp*(-b -tol)
                    dm2 = 2.0_wp*c/(-b -tol)
                else
                    dm = 2.0_wp*c/(-b +tol)
                    dm2 = 0.5_wp*(-b +tol)
                endif
                if (abs(dm2) < abs(dm)) then
                    dm = dm2
                endif
                mSulf = SMA(idSO4)-dm
                mH = SMC(idH)-dm
                mHSO4 = SMA(idHSO4)+dm
                calcviscosity = calcviscosity_saved
                call DiffKsulfuricDissoc(mHSO4, diffK)
                alphaHSO4 = 1.0_wp-(mHSO4/mHSO4max)
                if ((abs(diffK) > 0.1_wp .AND. abs(alphaHSO4 -alphaHSO4pred) > 0.4_wp) &
                    & .OR. abs(alphaHSO4 -alphaHSO4pred) > 0.99_wp) then
                    call DiffKsulfuricDissoc(mHSO4Guess1, diffK)
                    mHSO4 = SMA(idHSO4)
                endif
            endif
        endif
        alphaHSO4 = 1.0_wp - (mHSO4/mHSO4max)    !this is a generally valid expression for alphaHSO4
        diffKHSO4 = abs(diffK)
        !--
    else    !case where HSO4 is part of the system components, but at present input composition, mHSO4max is at zero molality.
        alphaHSO4 = -9.999999_wp
        call Gammas()                           !call Gammas once
    endif !@@**##

    end subroutine HSO4_dissociation
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
    module subroutine DiffKsulfuricDissoc(mHSO4inp, diffK)
               
    use ModSystemProp, only : idH, idHSO4, idSO4, topsubno
    use ModAIOMFACvar, only : gasrln, gamrln, galrln, gcsrln, gcmrln, gclrln, SMA, SMC, wtf, Tmolal

    implicit none
    !interface arguments:
    real(wp),intent(inout) :: mHSO4inp
    real(wp),intent(out) :: diffK
    !local vars:
    real(wp),parameter :: dtiny = 9.1_wp*epsilon(1.0_wp)    
    real(wp) :: lngammaH, lngammaSulf, lngammaHydSulf, mscale, mratio, madd
    !................................................................

    mscale = 1.0E-3_wp*deps  !scaling factor
    !limit mHSO4 to the possible interval bounds:
    mHSO4 = mHSO4inp
    mHSO4 = max(mHSO4, mHSO4min)
    mHSO4 = min(mHSO4, mHSO4max)

    !Calculate the molalities of the variable ions:
    !HSO4-
    if (mHSO4 < mscale) then                !make sure that at least a tiny amount of HSO4- is present
        madd = min(mscale -mHSO4, mHSO4max)
        mHSO4 = mHSO4 +madd
    endif
    SMA(idHSO4) = mHSO4
    !H+
    mH = max(mHmax-mHSO4, 0.0_wp)
    if (mH < mscale) then                   !make sure that at least a tiny amount of H+ is present
        madd = mscale - mH                  !the difference from the tiny mH value to mscale
        if (mHSO4-madd > mHSO4min) then
            mHSO4 = mHSO4 -madd
            SMA(idHSO4) = mHSO4
            mH = mH +madd
        endif
    endif
    SMC(idH) = mH
    !SO4--
    mSulf = max(mSulfmax-mHSO4, 0.0_wp)
    if (mSulf < mscale) then                !make sure that at least a tiny amount of SO4-- is present
        madd = mscale - mSulf               !the difference from the tiny mSulf value to mscale
        if (mHSO4-madd > mHSO4min) then
            mHSO4 = mHSO4 -madd
            SMA(idHSO4) = mHSO4
            mSulf = mSulf +madd
        endif
    endif
    SMA(idSO4) = mSulf

    if (mH*mSulf*mHSO4 < -dtiny) then  !check the values for mistakes in program
        !$OMP CRITICAL (DKD2)
        write(*,*) "WARNING: negative molality values detected! wtf(water): ", wtf(1)
        write(*,*) "mH, mSulf, mHSO4: ", mH, mSulf, mHSO4
        write(*,*) "mHSO4inp ", mHSO4inp
        write(*,*) ""
        !$OMP end CRITICAL (DKD2)
    endif
    !*-*
    call Gammas() !compute the natural log of activity coefficients from SR, MR, & LR parts.
    !*-*
    !Use the mole fraction-based ion activity coeff. (or their natural log).
    if (idH > 0) then
        lngammaH = -Tmolal +gcsrln(idH) +gcmrln(idH) +gclrln(idH)
    else
        lngammaH = 0.0_wp
    endif
    if (idHSO4 > 0) then
        lngammaHydSulf = -Tmolal +gasrln(idHSO4) +gamrln(idHSO4) +galrln(idHSO4)
    else
        lngammaHydSulf = 0.0_wp
    endif
    if (idSO4 > 0) then
        lngammaSulf = -Tmolal +gasrln(idSO4) +gamrln(idSO4) +galrln(idSO4)
    else
        lngammaSulf = 0.0_wp
    endif
    !estimate the molalities of SO4--, HSO4- and H+ to get agreement 
    !with the reference value of the dissociation constant.
    lngKbisulf = lngammaH +lngammaSulf -lngammaHydSulf

    !comparison, logarithmic "rel." deviation:
    if (mHSO4 > 0.0_wp) then
        mratio = mH*mSulf/mHSO4
        if (mratio > 0.0_wp) then 
            diffK = (lngKbisulf +log(mratio)) -lnKHSO4atT
        else
            mratio = deps/mHSO4
            diffK = (lngKbisulf +log(mratio)) -lnKHSO4atT
        endif
    else
        mratio = (mH*mSulf)/deps
        if (mratio > 0.0_wp) then 
            diffK = (lngKbisulf +log(mratio)) -lnKHSO4atT
        else
            mratio = deps
            diffK = (lngKbisulf +log(mratio)) -lnKHSO4atT
        endif
    endif

    end subroutine DiffKsulfuricDissoc 
    !================================================================================================= 
    
    
    !=================================================================================================
    !internal wrapper function for the computation of Diffk as a function of mHSO4 input.
    module function fHSO4dissoc(mHSO4in)       
    
    implicit none
    !interface:
    real(wp), intent(inout) :: mHSO4in
    real(wp) :: fHSO4dissoc              !the function return value
    !......................................
    !calculate the actual value of the dissociation (non)-equilibrium:
    call DiffKsulfuricDissoc(mHSO4in, fHSO4dissoc)
    
    end function fHSO4dissoc
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
    !*   -> latest changes: 2022-01-17                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    module subroutine HSO4_and_HCO3_dissociation()
    
    use Mod_MINPACK, only : hybrd1
    
    implicit none
    !Local variables and parameters:
    integer :: k, n, nrd, info, iflag, iloop, maxiloop
    !real parameters:    
    real(wp),parameter :: T0 = 298.15_wp              ![k] reference temperature for dissociation constant temperature dependency calculation
    real(wp),parameter :: fa1 = -8.204327E+02_wp, fa2 = -1.4027266E-01_wp, fa3 = 5.027549E+04_wp, fa4 = 1.268339E+02_wp, fa5 = -3.879660E+06_wp !terms/factors for the temperature dependent bicarbonate1 dissociation constant calculation.
    real(wp),parameter :: fb1 = -2.484192E+02_wp, fb2 = -7.489962E-02_wp, fb3 = 1.186243E+04_wp, fb4 = 3.892561E+01_wp, fb5 = -1.297999E+06_wp  !terms/factors for the temperature dependent bicarbonate2 dissociation constant calculation.
    real(wp),parameter :: fc1 = 2.495691E+02_wp, fc2 = 4.570806E-02_wp, fc3 = -1.593281E+04_wp, fc4 = -4.045154E+01_wp, fc5 = 1.541270E+06_wp   !terms/factors for the temperature dependent CO2 Henry's law constant
    real(wp),parameter :: fd1 = -1.4816780E+02_wp, fd2 = 8.933802E-01_wp, fd3 = -2.332199E-03_wp, fd4 = 2.146860E-06_wp
    !--
    real(wp) :: alphaHSO4pred, dperturb, mfHSO4maxinions, nCa, nCaSO4, r, sumCat, T, tol, wfHSO4max
    real(wp) :: ln_nCarbmax, ln_nCO2max, ln_nHCO3max, ln_nHmax, ln_nHSO4max, ln_nOHmax, wtf_water_init
    real(wp),dimension(8) :: molarinp
    real(wp),dimension(:),allocatable :: diffK, randomval
    !--
    logical,parameter :: memory_initial_guess = .true.   !if true, keeps determined solver values to be used for subsequent initial guess in calc.
    logical,parameter :: HighPrecSolving = .true.        !switch to set mode of dissociation calculation between .false. (=faster) and high-precision solving (slower);
    logical :: calcviscosity_saved
    !........................................................................
    
    !settings for calculations with CO2 equilibration to gas phase of constant mixing ratio.
    use_CO2gas_equil = .false.      !switch for use of pseudo CO2(g) <--> CO2(aq) equilibration to the "target_CO2gas_ppm" alongside the aqueous phase equilibria;
    target_CO2gas_ppm = 400.0_wp     !400 ppm by volume CO2(g) in 1 m^3 of air at ptot = 1 atm;
    
    !Initial molar amounts of ions + neutrals in the mixture for 1 kg of solvent at input => equivalent to molality (initially)
    wtf_water_init = wtf(1)         !save input mass frac. of water (since that may change when water self-dissociation is accounted for)
    allocate( mNeutral(nneutral), molNeutral(nneutral), SNC(NGI), SNA(NGI) )
    molNeutral = 0.0_wp
    SNA = SMA
    SNC = SMC
    call MassFrac2SolvMolalities(wtf, mNeutral)
    molNeutral = mNeutral
    
    sumCat = sum(SNC)
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
    if (bisulfsyst) then
        nHSO4 = SNA(idHSO4)         !initial molar amount of [HSO4-]
        nSulf = SNA(idSO4)          !initial molar amount of [SO4--]
    else
        nHSO4 = 0.0_wp
        nSulf = 0.0_wp
    endif
    !for tests:  make sure that nOHmax is at least 3.5*sumCat value, if possible:
    r = 3.5_wp*sumCat/nH2O_init
    r = min(r, 0.98_wp)
    nOHmax = r*nH2O_init + nOH
    
    !calculate maximum possible molar amounts: 
    if (use_CO2gas_equil) then
        n = 7
        nCO2gas = target_CO2gas_ppm*1.0E-6_wp/(T_K*Rgas_atm)                           !initial molar amount of [CO2(g)]; 
        nHmax = nH + nHCO3 + nCO2 +nCO2 + nCO2gas +nCO2gas + (nOHmax -nOH) + nHSO4  !H+ from 2*min(nCO2,nH2O), usually nH2O is more abundant
        nCarbmax = nCarb + nHCO3 + nCO2 + nCO2gas
    else
        n = 5
        nCO2gas = 0.0_wp
        nHmax = nH + nHCO3 + nCO2 + nCO2 + (nOHmax -nOH) + nHSO4
        nCarbmax = nCarb + nHCO3 + nCO2
    endif
    if (bisulfsyst) then
        n = n + 1
    endif
    nSulfmax = nHSO4 + nSulf
    nHSO4max = min(nHmax, nSulfmax)
    nHCO3max = min(nHmax, nCarbmax)
    nCO2max = min(0.5_wp*nHmax, nCarbmax)
    nVmax = 1.0E9_wp
    
    !n: the number of unknowns / equations to solve:
    allocate(solve_var(n), solve_var_maxval(n), map_solver2molarinp(n), randomval(n), diffK(n))
    map_solver2molarinp = 0                            
    if (.NOT. allocated(solve_var_saved)) then
        allocate(solve_var_saved(n))                    !allocate submodule variable on first use, but not subsequently
        solve_var_saved = -0.99_wp*huge(1.0_wp)         !initialized with an unfeasible value
    endif

    !if Ca++ is part of the input, remove all the CaSO4(s) that could form;
    if (idCa > 0) then
        nCa = SNC(idCa)
        nCaSO4 = min(nCa, nSulfmax)
        if (nCa < nSulfmax) then
            SNC(idCa) = 0.0_wp
            nSulf = nSulf + nHSO4 -nCa
            nH = nH + nHSO4
            nHSO4 = 0.0_wp
            !write(*,*) "WARNING: All the CaSO4 was precipitated out by default."
        else
            !leave tiny amount (1000.0_wp*deps) of SO4-- in the system  0.001% of Sulfmax or 1000.0_wp*deps
            SNC(idCa) = nCa -nSulfmax + min(1.0E-5_wp*nSulfmax, 1.0E3_wp*deps)
            nSulf = min(1.0E-5_wp*nSulfmax, 1.0E3_wp*deps)
            nH = nH + nHSO4
            nHSO4 = 0.0_wp            
            !write(*,*) "WARNING: All the CaSO4 was precipitated out by default"
        endif
        nSulfmax = nHSO4 + nSulf
        nHmax = nH + nHCO3 + nCO2+nCO2 + (nOHmax -nOH) + nHSO4
        nHSO4max = min(nHmax, nSulfmax)
        nHCO3max = min(nHmax, nCarbmax)
        nCO2max = min(0.5_wp*nHmax, nCarbmax)
    endif
        
    !if no H+ or no carbonate or no carbonate acid are in the solution at the same time, then step over, otherwise calculate dissociation:
    if (nHCO3max > ntiny) then !@@**## there is some dissociation, so proceed...
        calcviscosity_saved = calcviscosity
        calcviscosity = .false.
        ln_nHCO3max = log(nHCO3max*(1.0_wp -lztiny))
        ln_nCarbmax = log(nCarbmax*(1.0_wp -lztiny))
        ln_nCO2max = log(nCO2max*(1.0_wp -lztiny))
        ln_nOHmax = log(nOHmax*(1.0_wp -lztiny))
        ln_nHmax = log(nHmax*(1.0_wp -lztiny))
        if (bisulfsyst) then
            ln_nHSO4max = log(nHSO4max*(1.0_wp -lztiny))
        else
            ln_nHSO4max = 0.0_wp
        endif
        ln_maxval_inp = [ ln_nHCO3max, ln_nCarbmax, ln_nCO2max, ln_nOHmax, ln_nHSO4max, &
                        & ln_nHmax, ln_nCO2max, log(nVmax*(1.0_wp -lztiny)) ]  
        
        !compute the known equilibrium constants for the dissociation equilibria at T_K:
        T = T_K
        lnK1HCO3atT = fa1 + fa2*T + fa3/T + fa4*log(T) + fa5/(T**2)
        lnK2HCO3atT = fb1 + fb2*T + fb3/T + fb4*log(T) + fb5/(T**2)
        lnKCO2atT = fc1 + fc2*T + fc3/T + fc4*log(T) + fc5/(T**2)
        lnKH2OatT = fd1 + fd2*T + fd3*(T**2) + fd4*(T**3)   !Valid from 273K to 323 K.
        if (bisulfsyst) then
            !calculate the reference value of the 2nd dissociation constant of H2SO4 as a function of temperature 
            !using parameterization by Knopf et al. (2003) (corrected equation 16). 
            !The parameterization is valid in the T-range [179 K < T_K < 474 K], for lower temperatures, 
            !replace T with 180 K temporarily and use the ln(KHSO4_T) 
            !of that temperature as that might be the best estimate for that case.
            T = T_K
            if (T < 180.0_wp) then          !use T of 180 K for parameterization
                T = 180.0_wp
            else if (T > 473.0_wp) then     !use T of 473 K
                T = 473.0_wp
            endif
            lnKHSO4atT = -4.54916799587_wp -fe1*(fe2*(1.0_wp/T-fe3)-fe4*log(T/T0) -0.5_wp*dcpdT*(T-T0))
        endif
        
        !##! section to compute a good initial guess for the msulf1, msulf2 range based on a parameterization of alphaHSO4 as a function of mass frac. of HSO4(max) in the water portion of the solvent mixture.
        if (bisulfsyst) then
            !(A) calculate the special normalized mass fraction of HSO4max in water plus HSO4max: wfHSO4max; only the water fractional amount associated with HSO4max is considered (excluding water associated with other ions in the system).
            mfHSO4maxinions = 2.0_wp*nSulfmax/SumIonMolalities      !mole fraction of (HSO4max plus counter cations) relative to all ion molar amounts.
            wfHSO4max = nSulfmax*MolarMassHSO4/(nSulfmax*MolarMassHSO4 + mfHSO4maxinions*wtf(1)/sum(wtf(1:nneutral)))
            !(B) calculate the predicted degree of dissociation at this wfHSO4max value based on the parameterization:
            alphaHSO4pred = 1.0_wp - 1.0_wp/(1.0_wp + (1.0_wp/wfHSO4max**c1 - 1.0_wp/wfHSO4max**c2))**c4 &
                & + c3*wfHSO4max**c5*(1.0_wp-wfHSO4max)**1.75_wp
            alphaHSO4pred = min(alphaHSO4pred, 1.0_wp-deps)         !prevent values larger than 1.0_wp
            alphaHSO4pred = max(alphaHSO4pred, deps)
        endif

        call random_seed (size = nrd)
        call random_seed(PUT=[(k*211839831, k=1,nrd)])              !initialize the pseudo-random number generator each time with the same seed (for debugging reproducability)
        if (HighPrecSolving) then
            maxiloop = 40
            tol = 1.0E-2_wp*sqrtdeps
        else
            maxiloop = 9
            tol = sqrtdeps
        endif
        !---------------------
        do iloop = 1,maxiloop
            info = -1
            if (iloop == 1) then
                !set the index mapping between molarinp entries and the array of solver variables via map_solver2molarinp:
                if (bisulfsyst .AND. use_CO2gas_equil) then
                    map_solver2molarinp(1:n) = [1, 2, 3, 4, 5, 6, 7, 8]
                else if (bisulfsyst) then
                    map_solver2molarinp(1:n) = [1, 2, 3, 4, 5, 6]
                else if (use_CO2gas_equil) then
                    map_solver2molarinp(1:n) = [1, 2, 3, 4, 6, 7, 8]    !skip 5 in this case
                else
                    map_solver2molarinp(1:n) = [1, 2, 3, 4, 6]          !skip 5
                endif
                solve_var_maxval(1:n) = ln_maxval_inp(map_solver2molarinp(1:n))
                !initial guess
                if (memory_initial_guess .AND. solve_var_saved(1) > -0.98*huge(1.0_wp)) then 
                    !use saved solution from prior succesful calculation (previous data point input) as initial guess:
                    solve_var(1:n) = solve_var_saved(1:n)
                    call rboundsCheck( solve_var(1:n), minlim, ntiny, solve_var_maxval(1:n) )
                else
                    !use initial guess values; estimated from trial calculations distinguishing with/without bisulfate:
                    if (bisulfsyst) then
                        molarinp(1) = 1.60E-6_wp*nHCO3max       !nHCO3
                        molarinp(2) = 3.0_wp*deps*nCarbmax      !nCarb
                        if (use_CO2gas_equil) then
                            molarinp(3) = 2.0E-5_wp*nCO2max     !nCO2(aq)
                        else
                            molarinp(3) = 0.9999_wp*nCO2max     !nCO2(aq)
                        endif
                        molarinp(4) = 30.0_wp*deps*nOHmax       !nOH
                        nSulf = alphaHSO4pred*nSulfmax          !nSulf
                        molarinp(5) = nSulfmax -nSulf           !nHSO4
                        molarinp(6) = molarinp(1) +2.0_wp*molarinp(2) +molarinp(4) +molarinp(5) +2.0_wp*nSulf     !nH; start with neutral condition
                        molarinp(7) = 0.9999_wp*nCO2max         !nCO2gas
                        molarinp(8) = 40.0_wp                   !V
                    else
                        molarinp(1) = 0.13_wp*nHCO3max          !nHCO3
                        molarinp(2) = 0.6_wp*nCarbmax           !nCarb
                        molarinp(3) = 1.0E-5_wp*nCO2max         !nCO2
                        molarinp(4) = 2.0E-5_wp*nOHmax          !nOH
                        nSulf = 0.0_wp
                        molarinp(5) = 0.0_wp                    !nHSO4
                        molarinp(6) = molarinp(1) +2.0_wp*molarinp(2) +molarinp(4) +molarinp(5) +2.0_wp*nSulf     !nH; start with neutral condition
                        molarinp(7) = 0.28_wp*nCO2max           !nCO2gas
                        molarinp(8) = 20.0_wp                   !V
                    endif
                    solve_var(1:n) = log(molarinp(map_solver2molarinp(1:n)))
                endif
            else if (mod(iloop,2) /= 0) then
                !perturb all molar amounts:
                call random_number(randomval)                   !get random numbers in the interval [0.0_wp, 1.0_wp]
                randomval = -1.0_wp +2.0_wp*randomval           !values in interval [-1.0_wp, 1.0_wp]
                dperturb = min( 0.8_wp, 0.005_wp + 0.05_wp**(1.0_wp/real(iloop -1, kind=wp)) )
                solve_var = solve_var*(1.0_wp-dperturb +2.0_wp*randomval*dperturb)
                call rboundsCheck( solve_var(1:n), minlim, ntiny, solve_var_maxval(1:n) )
            else
                !perturb all molar amounts slightly:
                call random_number(randomval)
                randomval = -1.0_wp +2.0_wp*randomval
                dperturb = 2.0E-2_wp*iloop
                solve_var = solve_var*(1.0_wp-dperturb +2.0_wp*randomval*dperturb)
                call rboundsCheck( solve_var(1:n), minlim, ntiny, solve_var_maxval(1:n) )
            endif
            !--
            call hybrd1(DiffK_carb_sulf, n, solve_var, diffK, tol, info)
            !--
            if (info == 1) then
                if (memory_initial_guess) then
                    if (sum(abs(diffK)) < sqrtdeps) then
                        solve_var_saved(1:n) = solve_var(1:n)   !save found solver values for potential use in subsequent calc.
                    endif
                endif
                !k = 1                                          !for debugging breakpoint
                exit
            endif
        enddo !iloop
        !---------------------
            
        !BLOCK
        !    use ModSystemProp, only : nd, NGI, anionZ, cationZ
        !    real(wp) :: t1, t2, t3
        !    t1 = sum(SMA(1:NGI)*abs(anionZ(1:NGI)))             !test the electrical charge neutrality condition in the solution
        !    t2 = sum(SMC(1:NGI)*cationZ(1:NGI))
        !    t3 = t1 - t2
        !    if (t3 < 1.0E6_wp .AND. abs(t3) > 1.0E-9_wp) then  !test
        !        write(*,'(A, ES15.8)') "WARNING: Electrical charge neutrality condition violated! ", t3
        !        write(*,'(A, I0,1X,ES15.8)') "nd, wtf(1) ", nd, wtf(1)
        !    endif
        !end BLOCK
        
        !If success, check found zero and return best found nHCO3:
        calcviscosity = calcviscosity_saved
        call DiffK_carb_sulf(n, solve_var, diffK, iflag)
        molarinp(map_solver2molarinp(1:n)) = exp(solve_var(1:n))
        nHCO3 = molarinp(1)
        nCarb = molarinp(2)
        nCO2 = molarinp(3)
        alphaHCO3 = 1.0_wp -nHCO3/nHCO3max              !this is a generally valid expression for alphaHCO3
        !predict the equil. partial pressure of CO2 [Pa]:
        pCO2 = 1.013250E+05_wp * mNeutral(idCO2) * exp( min(gnmrln(idCO2), logval_threshold) - lnKCO2atT )
        if (bisulfsyst) then
            nHSO4 = molarinp(5)
            alphaHSO4 = 1.0_wp -nHSO4/nHSO4max          !this is a generally valid expression for alphaHSO4
        else
            alphaHSO4 = -9.999999_wp
        endif
        !checks for debugging:
        if (iloop > 9) then
            k = 88 !for breakpoint
            r = sum(abs(diffK))
            if (r > sqrtdeps) then
                errorflag_clist(17) = .true.
                !if (.NOT. frominpfile) then
                !    !$OMP CRITICAL
                !    write(*,'(A)') "AIOMFAC ERROR 17: Issue with ion dissociation equilibria calculations."
                !    write(*,'(A)') "The numerical solution of electrolyte/ion dissociation equilibria was &
                !        &not accomplished to the desired tolerance level. Model output for this data &
                !        &point is unreliable and likely incorrect."
                !    write(*,'(A,I0,ES13.6)') "iloop, wtf_water_init = ", iloop, wtf_water_init
                !    write(*,'(A,ES13.6)') "sum(abs(diffK)) = ", r
                !    write(*,*)
                !    !$OMP end CRITICAL
                !endif
            endif
        endif
    else    !case where HCO3 is part of the system components, but at present input composition, nHCO3max is at zero molar amount.
        alphaHCO3 = -9.999999_wp
        alphaHSO4 = -9.999999_wp
        pCO2 = 0.0_wp
        call Gammas() !call Gammas once
    endif !@@**##
    
    deallocate(solve_var, solve_var_maxval, map_solver2molarinp, randomval, diffK)
    deallocate(mNeutral, molNeutral, SNC, SNA)
    
    end subroutine HSO4_and_HCO3_dissociation
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
    subroutine DiffK_carb_sulf(n, solve_var, diffK, iflag)
    !subroutine containing the system of n equations with n unknowns to be solved using Minpack methods.

    !Public Variables:
    use ModSystemProp, only : nneutral, idH, idHCO3, idCO3, idOH, idCO2, idSO4, idHSO4, NGI, topsubno, &
        & waterpresent, Mmass, Ianion, Ication
    use ModSubgroupProp, only : SMWC, SMWA

    implicit none

    integer,intent(in) :: n
    real(wp),dimension(n),intent(inout) :: solve_var
    real(wp),dimension(n),intent(out) :: diffK
    integer,intent(out) :: iflag
    !local vars:
    integer :: I, J, k
    real(wp),parameter :: dtiny = 1.0E-2_wp*epsilon(1.0_wp)
    real(wp) :: deltaCO2, diffK1, diffK2, diffK3, diffK4, lngammaCarb, lngammaCO2, lngammaH, lngammaHydCarb, &
        & lngammaHydSulf, lngammaOH, lngammaSulf, lngammawv, mCarbcheck, mgeomean, mHCO3check, mOHcheck,    &
        & mratio1, mratio2, mratio3, mratio4, mSulfcheck, solvmass, totmass, xwater
    real(wp),dimension(n) :: exp_solve_var
    real(wp),dimension(nneutral) :: xneutral
    !................................................
        
    !make sure solver_var are consistent with upper limits given by ln_maxval_inp:
    iflag = 0
    call rboundsCheck( solve_var(1:n), minlim, ntiny, solve_var_maxval(1:n) )
    exp_solve_var = exp(solve_var)
    nHCO3 = exp_solve_var(1)
    nCarb = exp_solve_var(2)
    nCO2 = exp_solve_var(3)
    nOH = exp_solve_var(4)
    if (bisulfsyst .AND. use_CO2gas_equil) then
        nHSO4 = exp_solve_var(5)
        nH = exp_solve_var(6)
        nCO2gas = exp_solve_var(7)
        V = exp_solve_var(8)
    else if (bisulfsyst) then
        nHSO4 = exp_solve_var(5)
        nH = exp_solve_var(6)
    else if (use_CO2gas_equil) then
        nHSO4 = 0.0_wp
        nH = exp_solve_var(5)
        nCO2gas = exp_solve_var(6)
        V = exp_solve_var(7)
    else
        nHSO4 = 0.0_wp
        nH = exp_solve_var(5)
    endif
    nSulf = nSulfmax -nHSO4
    
    !update the moles of all species:
    deltaCO2 = nCO2 - nCO2_init
    molNeutral(1) = nH2O_init + deltaCO2 - (nOH - nOH_init)
    if (molNeutral(1) < ntiny) then     !correction necessary
        molNeutral(1) = ntiny
        deltaCO2 = molNeutral(1) -nH2O_init + (nOH - nOH_init)
        nCO2 = deltaCO2 + nCO2_init
    endif
    molNeutral(idCO2) = nCO2   !CO2(aq)
    SNA(idHCO3) = nHCO3
    SNA(idCO3) = nCarb
    SNA(idOH) = nOH
    SNC(idH) = nH
    if (bisulfsyst) then
        SNA(idHSO4) = nHSO4
        SNA(idSO4) = nSulf
    endif
        
    !update the molalities of all ions & solvents:
    call Moles2solvmass(molNeutral, solvmass)   
        
    totmass = solvmass
    do I = 1,NGI
        k = ICation(I) - 200
        J = IAnion(I) - 240
        if (k > 0) then
            totmass = totmass + SNC(I)*SMWC(k)*1.0E-3_wp
        endif 
        if (J > 0) then 
            totmass = totmass + SNA(I)*SMWA(J)*1.0E-3_wp 
        endif
    enddo
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
    if (bisulfsyst) then
        mHSO4 = SMA(idHSO4)
        mSulf = SMA(idSO4)
    endif

    !Mole fraction of water
    call Molality2SolvMoleFrac(SMA, SMC, mNeutral, xneutral)
    xwater = xneutral(1)
    !--
    call Gammas()   !compute the natural log of activity coefficients from SR, MR, & LR interaction parts.
    !--
    !Use the mole fraction-based ion activity coeff. (or their natural log).
    if (idH > 0) then
        lngammaH = -Tmolal +gcsrln(idH) +gcmrln(idH) +gclrln(idH)
    else
        lngammaH = 0.0_wp
    endif
    if (idHCO3 > 0) then
        lngammaHydCarb = -Tmolal +gasrln(idHCO3) +gamrln(idHCO3) +galrln(idHCO3)
    else
        lngammaHydCarb = 0.0_wp
    endif
    if (idCO3 > 0) then
        lngammaCarb = -Tmolal +gasrln(idCO3) +gamrln(idCO3) +galrln(idCO3)
    else
        lngammaCarb = 0.0_wp
    endif
    if (idHSO4 > 0) then
        lngammaHydSulf = -Tmolal +gasrln(idHSO4) +gamrln(idHSO4) +galrln(idHSO4)
    else
        lngammaHydSulf = 0.0_wp
    endif
    if (idSO4 > 0) then
        lngammaSulf = -Tmolal +gasrln(idSO4) +gamrln(idSO4) +galrln(idSO4)
    else
        lngammaSulf = 0.0_wp
    endif
    if (idOH > 0) then
        lngammaOH = -Tmolal +gasrln(idOH) +gamrln(idOH) +galrln(idOH)
    else
        lngammaOH = 0.0_wp
    endif
    if (idCO2 > 0) then
        lngammaCO2 = gnmrln(idCO2)
    else
        lngammaCO2 = 0.0_wp
    endif
    if (waterpresent) then
        lngammawv = gnsrln(1) +gnmrln(1) +gnlrln(1)
    else
        lngammawv = 0.0_wp
    endif
        
    !estimate the molalities of CO3--, HCO3- and H+ to get agreement
    !with the reference value of the dissociation constant.
    lngK1bicarb = lngammaH + lngammaHydCarb - lngammaCO2 - lngammawv
    lngK2bicarb = lngammaH + lngammaCarb - lngammaHydCarb
    lngbisulf = lngammaH + lngammaSulf - lngammaHydSulf
    lngKOH = lngammaH + lngammaOH - lngammawv

    !report resulting deviation vector diffK:
    diffK = 0.0_wp
    if (mCO2*xwater > 0.0_wp) then
        mratio1 = mH*mHCO3/(mCO2*xwater)
    else
        mratio1 = 0.0_wp
    endif
    if (mHCO3 > 0.0_wp) then
        mratio2 = mH*mCarb/mHCO3
    else
        mratio2 = 0.0_wp
    endif
    if (xwater > 0.0_wp) then
        mratio3 = mH*mOH/xwater
    else
        mratio3 = 0.0_wp
    endif
    if (bisulfsyst) then
        if (mHSO4 > 0.0_wp) then
            mratio4 = mH*mSulf/mHSO4
        else
            mratio4 = 0.0_wp
        endif
    endif
    
    if (mHCO3 > lztiny) then !normal case
        if (mratio1 > 0.0_wp) then
            diffK(1) = (lngK1bicarb +log(mratio1)) -lnK1HCO3atT
        else
            mratio1 = ztiny*mHCO3/(mCO2*xwater)
            diffK(1) = (lngK1bicarb +log(mratio1)) -lnK1HCO3atT
        endif
    else !very low mHCO3, requires special treatment
        if (mratio1 > 0.0_wp) then
            !compute the correct value of mHCO3 needed to fulfill this equation 1:
            diffK1 = lnK1HCO3atT - lngK1bicarb
            diffK1 = max(diffK1, -logval_threshold)
            diffK1 = min(diffK1, logval_threshold)
            mHCO3check = (mCO2*xwater)/mH*exp(diffK1)
            if (mHCO3check*solvmass < dtiny) then   !assume mHCO3check is essentially zero, so no equation to solve;
                diffK(1) = 0.0_wp                   !set for this exception
                nHCO3 = solvmass*mHCO3check
            else
                mgeomean = sqrt(mHCO3*mHCO3check)
                mratio1 = mH*mgeomean/(mCO2*xwater)
                diffK(1) = (lngK1bicarb +log(mratio1)) -lnK1HCO3atT
                nHCO3 = solvmass*mgeomean
            endif
            nHCO3 = min(nHCO3, nHCO3max)
        endif
    endif !mHCO3
    if (mCarb > lztiny) then  !normal case
        if (mratio2 > 0.0_wp) then
            diffK(2) = (lngK2bicarb +log(mratio2)) -lnK2HCO3atT
        else
            mratio2 = ztiny/mHCO3
            diffK(2) = (lngK2bicarb +log(mratio2)) -lnK2HCO3atT
        endif
    else !very low mCarb, requires special treatment
        if (mratio2 > 0.0_wp) then
            !compute the correct value of mCarb needed to fulfill this equation 2:
            diffK2 = lnK2HCO3atT - lngK2bicarb
            diffK2 = max(diffK2, -logval_threshold)
            diffK2 = min(diffK2, logval_threshold)
            mCarbcheck = mHCO3/mH*exp(diffk2)
            if (mCarbcheck*solvmass < dtiny) then   !assume mCarb is essentially zero, so no equation to solve;
                diffK(2) = 0.0_wp                   !set for this exception
                nCarb = solvmass*mCarbcheck
            else
                mgeomean = sqrt(mCarb*mCarbcheck)
                mratio2 = mH*mgeomean/mHCO3
                diffK(2) = (lngK2bicarb +log(mratio2)) -lnK2HCO3atT
                nCarb = solvmass*mgeomean
            endif
            nCarb = min(nCarb, nCarbmax)
        else
            mratio2 = ztiny/mHCO3
            diffK(2) = (lngK2bicarb +log(mratio2)) -lnK2HCO3atT
        endif
    endif !nCarb
    if (mOH > lztiny) then  !normal case
        if (mratio3 > 0.0_wp) then
            diffK(3) = (lngKOH +log(mratio3)) -lnKH2OatT
        else
            mratio3 = mH*ztiny/xwater
            diffK(3) = (lngKOH +log(mratio3)) -lnKH2OatT
        endif
    else !very low mOH; requires special treatment
        !compute the correct value of mOH needed to fulfill equilibrium equation 3:
        if (mratio3 > 0.0_wp) then
            diffK3 = lnKH2OatT - lngKOH
            diffK3 = max(diffK3, -logval_threshold)
            diffK3 = min(diffK3, logval_threshold)
            mOHcheck = xwater/mH*exp(diffK3)
            if (mOHcheck*solvmass < dtiny) then     !assume mCarb is essentially zero, so no equation to solve;
                diffK(3) = 0.0_wp
                nOH = solvmass*mOHcheck
            else
                mgeomean = sqrt(mOH*mOHcheck)
                mratio3 = mH*mgeomean/xwater
                diffK(3) = (lngKOH +log(mratio3)) -lnKH2OatT
                nOH = solvmass*mgeomean
            endif
            nOH = min(nOH, nOHmax)
        endif
    endif !mOH
    if (bisulfsyst) then
        if (mSulf > lztiny .AND. mHSO4 > lztiny) then
            if (mratio4 > 0.0_wp) then
                diffK(4) = (lngbisulf +log(mratio4)) -lnKHSO4atT
            else
                mratio4 = mH*ztiny/xwater
                diffK(4) = (lngbisulf +log(mratio4)) -lnKHSO4atT
            endif
        else    !very low mSulf, requires special treatment
            !compute the "correct" value of mSulf needed to fulfill equilibrium equation 4:
            if (mratio4 > 0.0_wp) then
                diffK4 = lnKHSO4atT - lngbisulf
                diffK4 = max(diffK4, -logval_threshold)
                diffK4 = min(diffK4, logval_threshold)
                mSulfcheck = mHSO4/mH*exp(diffk4)
                if (mSulfcheck*solvmass < dtiny) then
                    diffK(4) = 0.0_wp
                    nSulf = solvmass*mSulfcheck
                else
                    mgeomean = sqrt(mSulf*mSulfcheck)
                    mratio4 = mH*mgeomean/mHSO4
                    diffK(4) = (lngbisulf +log(mratio4)) -lnKHSO4atT
                    nSulf = solvmass*mgeomean
                endif
                nSulf = min(nSulf, nSulfmax -2.0_wp*deps*nSulfmax)
                nHSO4 = nSulfmax -nSulf
            endif
        endif !mSulf
    endif

    if (bisulfsyst .AND. use_CO2gas_equil) then
        !write(*,'(*(ES13.6,1X))') [nHCO3, nCarb, nCO2, nOH, nHSO4, nH, nCO2gas, V]
        solve_var = log([nHCO3, nCarb, nCO2, nOH, nHSO4, nH, nCO2gas, V])
        diffK(5) = mCO2*exp(min(lngammaCO2, logval_threshold))*2.5E3_wp - exp(lnKCO2atT)
        diffK(6) = 1.0E2_wp*(sum_sorted([nHmax, -nHCO3, -nH, -nCO2, -nCO2, -nCO2gas, -nCO2gas, -nHSO4, -nOHmax, nOH]))
        diffK(7) = 1.0E2_wp*(sum_sorted([nCarbmax, -nHCO3, -nCarb, -nCO2, -nCO2gas]))                   !molar balance
        diffK(8) = target_CO2gas_ppm*1.0E-6_wp*V -nCO2gas*Rgas_atm*T_K
    else if (bisulfsyst) then
        solve_var = log([nHCO3, nCarb, nCO2, nOH, nHSO4, nH])
        diffK(5) = 1.0E2_wp*(sum_sorted([nCarbmax, -nHCO3, -nCarb, -nCO2]))                            !nCarbmax -nHCO3 -nCarb -nCO2 !mass balance
        diffK(6) = 1.0E2_wp*(sum_sorted([nHmax, -nHCO3, -nH, -nCO2, -nCO2, -nHSO4, -nOHmax, nOH]))     !nHmax -nHCO3 -nH -2.0_wp*nCO2 -nHSO4 -nOHmax + nOH !mass balance
    else if (use_CO2gas_equil) then
        solve_var = log([nHCO3, nCarb, nCO2, nOH, nH, nCO2gas, V])
        diffK(4) = mCO2*exp(min(lngammaCO2, logval_threshold))*2.5E3_wp - exp(lnKCO2atT)
        diffK(5) = 1.0E2_wp*(sum_sorted([nHmax, -nHCO3, -nH, -nCO2, -nCO2, -nCO2gas, -nCO2gas, -nOHmax, nOH]))
        diffK(6) = 1.0E2_wp*(sum_sorted([nCarbmax, -nHCO3, -nCarb, -nCO2, -nCO2gas]))                   !molar balance
        diffK(7) = target_CO2gas_ppm*1.0E-6_wp*V -nCO2gas*Rgas_atm*T_K
    else
        solve_var = log([nHCO3, nCarb, nCO2, nOH, nH])
        diffK(4) = 1.0E2_wp*(sum_sorted([nHmax, -nHCO3, -nH, -nCO2, -nCO2, -nOHmax, nOH]))
        diffK(5) = 1.0E2_wp*(sum_sorted([nCarbmax, -nHCO3, -nCarb, -nCO2]))
    endif

    end subroutine DiffK_carb_sulf
    !=================================================================================================
    
    
    !Utility subroutine to check, and if necessary adjust, values of variables limited to upper/lower bounds.
    !The upper bound is given by a slightly scaled maxlim value.
    !--------------------------------------------------------------------------------------------------
    pure elemental module subroutine rboundsCheck(rset, lowerbound, ntiny, maxlim)
    !
    implicit none
    !interface arguments:
    real(wp),intent(inout) :: rset
    real(wp),intent(in) :: lowerbound, ntiny, maxlim
    !local parameters:
    real(wp),parameter :: scaleval = 1.0E-3_wp, deps = epsilon(1.0_wp), threedeps = 3.0_wp*deps
    real(wp) :: qa, qacorr, qaband, higherbound
    !.......................

    !use smart adjustments on values with upper/lower bounds if necessary:
    higherbound = (1.0_wp -ntiny)*maxlim
    qaband = 0.05_wp*abs(higherbound - lowerbound)
    if (qaband < threedeps) then
        qaband = max(0.2_wp*abs(higherbound - lowerbound), deps)
    endif

    if (rset < lowerbound) then
        qacorr = (lowerbound-rset)*scaleval
        qa = lowerbound + qacorr
        do while (qa > (lowerbound +qaband)) !(qa > 1.0_wp)
            qacorr = abs(qacorr)*0.1_wp
            qa = lowerbound + qacorr
        enddo
        rset = qa
    else if (rset > higherbound) then
        qacorr = (rset -higherbound)*scaleval
        qa = higherbound -qacorr
        do while (qa < (higherbound -qaband)) !(qa < 0.0_wp)
            qacorr = abs(qacorr)*0.1_wp
            qa = higherbound -qacorr
        enddo
        rset = qa
    endif

    end subroutine rboundsCheck
    !--------------------------------------------------------------------------------------------------
    
    
    !--------------------------------------------------------------------------------------------------
    !Utility function:
    !** A safer way for summing a "small" list of values that potentially differ substantially in magnitude ** 
    !** based on sorting them and summing the elements starting with the smallest absolute value.           **
    pure module function sum_sorted(list) result(summed)
        
    implicit none
    !interface arguments:
    real(wp),dimension(:),intent(in) :: list
    real(wp) :: summed
    !local arguments
    integer :: i, k
    logical,dimension(size(list)) :: available
    !...........................
        
    available = .true.
    !sort the list from smallest to largest in magnitude (ignoring neg/pos signs):
    summed = 0.0_wp
    do i = 1,size(list)
        k = minloc(abs(list(:)), mask = available(:), DIM=1)   
        available(k) = .false.
        summed = summed + list(k) 
    enddo
        
    end function sum_sorted
    !--------------------------------------------------------------------------------------------------
    

end submodule SubModDissociationEquil