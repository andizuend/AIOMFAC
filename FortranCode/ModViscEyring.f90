!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   This module contains subroutines to compute the viscosity of aqueous electrolyte   *
!*   solutions with treatment of single-ion and cation--anion pair contributions and    *
!*   associated parameters; see paper by Lilek and Zuend (2021, ACP).                   *
!*                                                                                      *
!*   Options for viscosity calculation modes for organic--inorganic mixtures are set in *
!*   the module declaration part below.                                                 *
!*                                                                                      *
!*   :: Authors & Copyright ::                                                          *
!*   Joseph Lilek, Andi Zuend,                                                          *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2020                                                            *
!*   -> latest changes: 2021-11-29                                                      *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  SUBROUTINE AqueousElecViscosity                                                 *
!*   -  SUBROUTINE WaterMolefracCorrection                                              *
!*   -  SUBROUTINE GoldsackViscEqn                                                      *
!*   -  SUBROUTINE MapIons2Fitpars3                                                     *
!*   -  SUBROUTINE ViscGelContribution                                                  *
!*                                                                                      *
!****************************************************************************************
MODULE ModViscEyring

USE ModAIOMFACvar, ONLY : galrln, gamrln, gasrln, gclrln, gcmrln, gcsrln, gnlrln, &
    & gnmrln, gnsrln, lnactcoeff_a, lnactcoeff_c, SMA, SMC, solvmixcorrMRc, Tmolal, &
    & TmolalSolvmix
USE ModSystemProp, ONLY : Ianion, Ication, Mmass, Nanion, Ncation, nd, solvmixrefnd, topsubno
USE ModSRparam, ONLY : SR_RR

IMPLICIT NONE

INTEGER(4),PRIVATE :: i
INTEGER(4),DIMENSION(1:69),PARAMETER,PRIVATE :: fitmapfwd = [(i, i = 1,69)]
INTEGER(4),DIMENSION(2,201:topsubno),PRIVATE :: ionfitmapfwd          
INTEGER(4),DIMENSION(2,201:240,241:topsubno),PRIVATE :: catanfitmapfwd
INTEGER(4),DIMENSION(:),ALLOCATABLE,PRIVATE :: ions  
REAL(8),DIMENSION(201:topsubno),PRIVATE :: Vion, Zion                 
REAL(8),DIMENSION(201:240,241:topsubno),PRIVATE :: nuecat, nuean      
REAL(8),PUBLIC :: delGstar_by_RT_w, ionicstrengthsave
!LOGICAL(4),PRIVATE :: gel_eqn1, gel_eqn2 

!private variables with initial values that may change:
LOGICAL(4),PRIVATE :: fitmap_complete = .true.                  !fitmap for viscosity model parameter init; needs to be .true. here; changed after first use
LOGICAL(4),PRIVATE :: calc_gel = .false.                        !false by default; currently gel effects are in development
!****
!**** set a few parameters viscosity model calculation options:
LOGICAL(4),PARAMETER,PUBLIC :: aquelec = .true.                 !true is default, false sets 'aquorg' mixing as mixing rule (except when ZSRvisc_on is true);
LOGICAL(4),PARAMETER,PUBLIC :: ZSRvisc_on = .false.             !false is default, true activates ZSR viscosity subroutine for org--inorg mixing.
LOGICAL(4),PARAMETER,PRIVATE :: calc_catan = .true.
LOGICAL(4),PARAMETER,PRIVATE :: newcatanPrime = .true.          !if true, use the non-normalized tau'. If false (default) use normalized tau.
!****
!****
REAL(8),DIMENSION(1:69),PARAMETER,PUBLIC :: fitpar = [ &
! 1 august 2021 (newcatan, newcatanprime; 46 input files)
&  1.55863D+00,  1.01039D-11,  4.98498D+00,  2.85945D-02,  5.81319D+00,  1.33891D-01,  7.27692D-01,  1.99747D-01, & 
&  5.20662D+00,  1.73959D+00,  1.52373D+01,  1.14135D-01,  1.98746D+00,  1.01039D-11,  1.45522D+01,  1.01039D-11, &
&  6.27643D+00,  2.06176D-01,  3.39886D+00,  7.30610D-01,  4.47385D+00,  1.01039D-11,  1.89144D+01,  3.58274D-01, &
&  5.20066D+00,  3.33591D-02,  2.47602D+00,  5.25846D-01,  1.90547D+00,  6.68723D-01,  8.80106D-01,  1.17183D+00, &
&  1.01039D-11,  5.23501D+00,  1.61773D+00,  1.64799D+00,  1.58764D+00,  1.67362D+00,  1.17183D+00,  1.18587D-01, &
&  5.23501D+00,  4.88160D-01,  1.01044D+00,  1.65761D+00,  2.01138D+00,  5.16958D+00,  6.11219D-01,  4.40650D+00, &
&  3.40868D+00,  6.02049D+00,  3.67356D+00,  4.06981D+00,  8.28882D+00,  9.75043D-01,  8.28882D+00,  2.04407D-01, &
&  2.40928D+00,  2.40928D+00,  2.40928D+00,  2.40928D+00,  2.40928D+00,  2.40928D+00,  1.54187D+00,  2.60910D+00, &
&  2.07340D+00,  2.41020D+00,  1.17183D+00,  1.01039D-11,  5.23501D+00 ]

!public procedures:
PUBLIC :: AqueousElecViscosity, WaterMolefracCorrection
PRIVATE     !default setting for procedures and variables

!make certain variables of this module threadprivate for use in parallel execution with openMP:
!$OMP THREADPRIVATE( ions, delGstar_by_RT_w, Zion, ionicstrengthsave)

!=====================================
    CONTAINS
!=====================================   
    
    !-------------------------------------------------------------------------------------
    SUBROUTINE AqueousElecViscosity(X_, lnGaSR_, RS_, ln_etaw, ln_eta_aquelec)
    
    USE ModAIOMFACvar, ONLY : galrln, gamrln, gasrln, gclrln, gcmrln, &
        & gcsrln, SMC, SMA, lnactcoeff_c, lnactcoeff_a, solvmixcorrMRc, solvmixcorrMRa, Tmolal, &
        & TmolalSolvmix, actcoeff_c, actcoeff_a, wtf, X
    USE ModSystemProp, ONLY : Ianion, Ication, Nanion, Ncation, topsubno, solvmixrefnd, nneutral, &
        & anionZ, cationZ    
    USE ModCompScaleConversion, ONLY : MassFrac2MoleFracMolality, zSolution2SpeciesMolality

    IMPLICIT NONE

    !interface variables:
    REAL(8),DIMENSION(:),INTENT(IN) :: X_, lnGaSR_, RS_
    REAL(8),INTENT(IN) :: ln_etaw
    REAL(8),INTENT(OUT) :: ln_eta_aquelec
    !local variables:
    INTEGER(4) :: I, NKNpNcat, nnp1, allocstat
    INTEGER(4) :: iion, icat, ian
    REAL(8),PARAMETER :: dtiny = EPSILON(1.0D0) 
    REAL(8),PARAMETER :: Vw = 1.39564D-5  ![m^3/mol]  Vw for water molec. volume;
    REAL(8) :: trunc, ionicstrength1, sumXion, sumX, ionicstrengthfactor, Vcomp, xw
    REAL(8),DIMENSION(201:topsubno) :: act_ion, X_ion, X_ion_organicfree
    !................................................
    
    !tasks that only need to be done once during program run:
    IF (fitmap_complete) THEN
        fitmap_complete = .false.
        Vion = 15.17D-06 * SR_RR(201:topsubno)
        CALL MapIons2Fitpars3()
    ENDIF
    
    !initialize:    
    Zion = 0.0D0
    X_ion = 0.0D0
    X_ion_organicfree = 0.0D0
    act_ion = 0.0D0
    
    NKNpNcat = nneutral + Ncation
    nnp1 = nneutral + 1

    gcsrln(1:Ncation) = lnGaSR_(nnp1:NKNpNcat)
    gasrln(1:Nanion)  = lnGaSR_(NKNpNcat+1:NKNpNcat+Nanion)
    
    IF (aquelec .AND. wtf(1) > 0.0D0) THEN    
        ionicstrengthfactor = SUM(wtf(1:nneutral))/wtf(1)
    ELSE
        ionicstrengthfactor = 1.0D0
    ENDIF
       
    !the activity coefficient contributions for LR and MR parts come via module access (USE ModAIOMFACvar, ONLY : ...)
    !initialize the output arrays:
    lnactcoeff_c = -9999.9D0
    lnactcoeff_a = -9999.9D0
    actcoeff_c = 0.0D0
    actcoeff_a = 0.0D0
    
    !ln of the activity coefficient for the cations:
    DO I = 1,Ncation
        IF (SMC(I) > dtiny) THEN
            lnactcoeff_c(I) = gcmrln(I) +gcsrln(I) +gclrln(I) -Tmolal  !this term converts to the molality scale (basis/scale conversion)
            IF (solvmixrefnd) THEN !correction terms for MR and SR part, because reference solution is the solvent mixture
                lnactcoeff_c(I) = lnactcoeff_c(I) +solvmixcorrMRc(I) +Tmolal -TmolalSolvmix    
            ENDIF
            IF (lnactcoeff_c(I) < -340.0D0 .OR. lnactcoeff_c(I) > 340.0D0) THEN 
                trunc = 1.0D-2*LOG( ABS(lnactcoeff_c(I)) )
                IF (lnactcoeff_c(I) < -340.0D0) THEN
                    actcoeff_c(I) = EXP(-340.0D0 - trunc)    
                ELSE
                    actcoeff_c(I) = EXP(340.0D0 + trunc)
                ENDIF
            ELSE
                actcoeff_c(I) = EXP(lnactcoeff_c(I))
            ENDIF
        ENDIF
    ENDDO
    
    !ln of the activity coefficient for the anions:
    DO I = 1,Nanion 
        IF (SMA(I) > dtiny) THEN
            lnactcoeff_a(I) = gamrln(I) +gasrln(I) +galrln(I) -Tmolal    !this term converts to the molality scale (basis/scale conversion)
            IF (solvmixrefnd) THEN  !correction terms for MR and SR because reference solution is the solvent mixture
                lnactcoeff_a(I) = lnactcoeff_a(I) +solvmixcorrMRa(I) +Tmolal -TmolalSolvmix 
            ENDIF
            IF (lnactcoeff_a(I) < -340.0D0 .OR. lnactcoeff_a(I) > 340.0D0) THEN 
                trunc = 1.0D-2*LOG( ABS(lnactcoeff_a(I)) )
                IF (lnactcoeff_a(I) < -340.0D0) THEN
                    actcoeff_a(I) = EXP(-340.0D0 - trunc)  
                ELSE
                    actcoeff_a(I) = EXP(340.0D0 + trunc) 
                ENDIF
            ELSE
                actcoeff_a(I) = EXP(lnactcoeff_a(I))
            ENDIF
        ENDIF
    ENDDO    
    
    ALLOCATE( ions(Ncation+Nanion), STAT=allocstat)
    icat = 0
    ian = 0
    
    !define mole fraction, activity of cations
    DO I = 1,Ncation
       !save activity and mole fraction of this ion in an array by actual AIOMFAC ion index:
        iion = Ication(I)   !iion is ion index number (201:topsubno)
        ions(I) = iion
        Zion(iion) = cationZ(I)
        act_ion(iion) = actcoeff_c(I)*SMC(I) * ionicstrengthfactor !ionicstrengthfactor to change solvent in molality to H2O instead of all neutral components
        X_ion(iion) = X_(nneutral+I)                            !these ion mole fractions are based on the full system (water, orgs, cation, anion)        
    ENDDO

    !define mole fraction, activity of anions
    DO I = 1,Nanion 
        iion = Ianion(I)
        ions(Ncation+I) = Ianion(I)
        Zion(iion) = ABS(anionZ(I))
        act_ion(iion) = actcoeff_a(I)*SMA(I)*ionicstrengthfactor
        X_ion(iion) = X_(nneutral+Ncation+I)
    ENDDO
    
    !ionicstrengthfactor  is used to change solvent to H2O instead of all neutral components:
    ionicstrength1 = 0.5D0*ionicstrengthfactor*( SUM(SMC(1:Ncation)*cationZ(1:Ncation)**2) + SUM(SMA(1:Nanion)*anionZ(1:Nanion)**2) )
    ionicstrengthsave = ionicstrength1
    
    IF (aquelec) THEN
        sumXion = SUM(X_ion(:))
        X_ion_organicfree = X_ion(:) / (X_(1) + sumXion)
        xw = 1.0D0 - SUM(X_ion_organicfree(:))
        CALL GoldsackViscEqn(X_ion_organicfree, act_ion, ionicstrength1, xw, Vw, ln_etaw, ln_eta_aquelec)
    ELSE !aquorg. in GoldsackViscEqn, xw = 1 - sum(X_ion), so mole fractions of organics are effectively counted as mole fraction of water
        xw = 1.0D0 - SUM(X_ion(:))
        sumX = SUM(X(1:nneutral))
        Vcomp = 15.17D-06 * SUM((X(1:nneutral)/sumX) * RS_(1:nneutral))   !use X (full system mole fractions to calc solvent volume)
        CALL GoldsackViscEqn(X_ion, act_ion, ionicstrength1, xw, Vcomp, ln_etaw, ln_eta_aquelec) !
    ENDIF
    
    DEALLOCATE( ions )

    END SUBROUTINE AqueousElecViscosity
    !-------------------------------------------------------------------------------------
    
    
    !-------------------------------------------------------------------------------------
    PURE SUBROUTINE GoldsackViscEqn(xin_, actin_, ionicstrength_, xw, Vw, ln_etaw, ln_etacalc)
    
    USE ModSystemProp, ONLY : Ianion, Ication, Nanion, Ncation, topsubno

    IMPLICIT NONE
    
    REAL(8),DIMENSION(201:topsubno),INTENT(IN) :: xin_, actin_
    REAL(8),INTENT(IN) :: ionicstrength_, ln_etaw, Vw, xw
    REAL(8),INTENT(OUT) :: ln_etacalc
    !local variables:
    INTEGER(4) :: iion, icat, ian
    ! Goldsack & Franchetto vars and parameters:
    INTEGER(4) :: cation, anion, numion                     !which value they refer to in the AIOMFAC ion numbering system (201:261)
    REAL(8) :: delGstar_over_RT, delGstar_by_RT_w, V, Garg_both, Garg_ion
    REAL(8) :: delGstar_ions, delGstar_catans, sumtauPrime, sumXz 
    REAL(8),DIMENSION(201:topsubno) :: chargefrac, delGstar_ion
    REAL(8),DIMENSION(201:240,241:topsubno) :: delGstar_catan
    REAL(8),DIMENSION(201:240,241:topsubno) :: tauPrime, tau
    REAL(8),PARAMETER :: deps = EPSILON(1.0D0)
    REAL(8),PARAMETER :: R = 8.3144598D0
    REAL(8),PARAMETER :: hPlanck = 6.62607015D-34           ![J s]        Planck's constant
    REAL(8),PARAMETER :: NA = 6.02214076D23                 ![#/mol]      Avogadro's constant
    REAL(8),PARAMETER :: MmassH2O = 1.801528D-02            !kg mol-1
    REAL(8),PARAMETER :: hNA = hPlanck*NA                   ![J s #/mol]
    REAL(8),PARAMETER :: ln_hNA = LOG(hNA)                  ![-] (normalized)
    !............................................
    
    !Note: xw = 1.0D0 - SUM(xin_(:))
    V = xw*Vw + fitpar(1)*SUM( xin_(:) * Vion(:) )          ![m^3/mol]  effective average hole volume; use fit parameter to affect ion vol. contributions (but not water);
        
    !(1) delGstar_ions (Garg_ion for all individual ions, summed)
    delGstar_ion = 0.0D0
    DO iion = 1,Ncation + Nanion
        numion = ions(iion)
        IF (ionfitmapfwd(1,numion) /= 0) THEN
            Garg_ion = actin_(numion)
            IF (Garg_ion < deps) THEN
                Garg_ion = 1.0D0
            ENDIF
            delGstar_ion(numion) = xin_(numion) * (fitpar(ionfitmapfwd(1,numion)) * LOG(Garg_ion) + fitpar(ionfitmapfwd(2,numion)))
        ELSE
            delGstar_ion(numion) = 0.0D0
            !add an error flag, indicating that this ion has not been considered
        ENDIF
    ENDDO
    delGstar_ions = SUM( delGstar_ion(:) )

    !(2) deltaG_catan (Garg_both of all cation--anion pairs, summed)
    IF (calc_catan) THEN
        !initialize
        delGstar_catan = 0.0D0
        tauPrime = 0.0D0
        tau = 0.0D0
        sumtauPrime = 0.0D0
        chargefrac = 0.0D0
        !use ion molalities
        sumXz = SUM(xin_(201:240)*ABS(Zion(201:240)))
        DO icat = 1,Ncation!ncat
            cation = Ication(icat)
            DO ian = 1,Nanion
                anion = Ianion(ian)
                !only one fitpar per catan pair
                IF (catanfitmapfwd(1,cation,anion) /= 0 ) THEN
                    Garg_both = fitpar(catanfitmapfwd(1,cation,anion)) * SQRT(ionicstrength_)
                    IF (Garg_both < deps) THEN
                        Garg_both = 1.0D0
                    ENDIF
                    !Implementation of new catan treatment (17 June 2021)
                    chargefrac(anion) = (xin_(anion)*ABS(Zion(anion))) / sumXz
                    tauPrime(cation,anion) = xin_(cation)/nuecat(cation,anion) * chargefrac(anion)
                    delGstar_catan(cation,anion) = Garg_both
                ELSE
                    delGstar_catan(cation,anion) = 0.0D0
                ENDIF
            ENDDO
        ENDDO
        sumtauPrime = SUM(tauPrime(:,:))
        IF (newcatanPrime .OR. sumtauPrime < deps) THEN
            WHERE (delGstar_catan(:,:) > 0.0D0)
                delGstar_catan(:,:) = delGstar_catan(:,:)*tauPrime(:,:)
            ENDWHERE
        ELSE
            tau = tauPrime/sumtauPrime
            WHERE (delGstar_catan(:,:) > 0.0D0)
                delGstar_catan(:,:) = delGstar_catan(:,:)*tau(:,:)
            ENDWHERE
        ENDIF
        delGstar_catans = SUM(delGstar_catan(:,:))
    ELSE
        delGstar_catans = 0.0D0
    ENDIF !calc_catan
    
    !*@#! delGstar_by_RT_w = xw * LOG( etaw * Vw / hNA )
    delGstar_by_RT_w = xw *(ln_etaw +LOG(Vw) -ln_hNA)
    !(4) deltaG_water + delGstar_ions + delGstar_catan
    delGstar_over_RT = &
        &   delGstar_ions &         !all cation & anion contributions multiplied by xin
        & + delGstar_catans &       !ion strength term multiplied by tau weighting
        & + delGstar_by_RT_w        !water contribution multiplied by xw

    !*@#! etacalc_ = (hNA / V) * EXP(delGstar_over_RT)  ![Pa s] the calculated mixture viscosity at this composition and T;
    ln_etacalc = ln_hNA -LOG(V) + delGstar_over_RT      !normalized ln(eta/[Pa s])

    END SUBROUTINE GoldsackViscEqn
    !-------------------------------------------------------------------------------------

    
    !-------------------------------------------------------------------------------------
    PURE SUBROUTINE WaterMolefracCorrection(X_, Xnew)
    
    USE ModAIOMFACvar, ONLY : wtf
    USE ModSystemProp, ONLY : nneutral, Mmass
        
    IMPLICIT NONE
        
    !interface vars
    REAL(8),DIMENSION(:),INTENT(IN) :: X_
    REAL(8),DIMENSION(:),INTENT(OUT) :: Xnew 
    !local vars
    REAL(8) :: Xaquorg, AquorgMmass, TotalAvgMmass, wtfaquorg
    !................................
    
    !(1) Calculate mass of water + organics:
    wtfaquorg = SUM(wtf(1:nneutral))
    !(2) Calculate avg molar mass of the aquorg mixture:
    AquorgMmass = SUM( (wtf(1:nneutral)/SUM(wtf(1:nneutral))) * Mmass(1:nneutral) )
    TotalAvgMmass = SUM(wtf * Mmass )       
        
    !(3) Calculate mole fraction amount for aquorg mixture:
    Xaquorg = (wtfaquorg / AquorgMmass) * TotalAvgMmass
        
    Xnew(1) = Xaquorg
    Xnew(2:nneutral) = 0.0D0
    Xnew(nneutral+1:) = X_(nneutral+1:)
    Xnew = Xnew / SUM(Xnew)
    
    END SUBROUTINE WaterMolefracCorrection
    !-------------------------------------------------------------------------------------
    
    
    !===============================================================================================    
    ! This subroutine maps ions to fit parameters.
    !    
    !*** AIOMFAC: ions: **************************************************************************************
    !    201 = Li+, 202 = Na+, 203 = K+, 204 = NH4+, 205 = H+, 221 = Ca2+, 222 = Ba2+
    !    223 = Mg2+, 224 = Sr2+, 225 = Co2+, 226 = Ni2+, 227 = Cu2+,
    !    228 = Zn2+, 229 = Hg2+, 241 = F-, 242 = Cl-, 243 = Br-, 244 = I-
    !    245 = NO3-, 246 = CH3COO-, 247 = SCN-, 248 = HSO4-, 249 = CH3SO3-, 261 = SO4--
    !*******************************************************************************************************    
    SUBROUTINE MapIons2Fitpars3()

    IMPLICIT NONE
    !......................
    
    !==========Ion fitpars=====================
    ionfitmapfwd(:,205) = fitmapfwd(2:3)   !H+
    ionfitmapfwd(:,201) = fitmapfwd(4:5)   !Li+
    ionfitmapfwd(:,203) = fitmapfwd(6:7)   !K+
    ionfitmapfwd(:,202) = fitmapfwd(8:9)   !Na+
    ionfitmapfwd(:,221) = fitmapfwd(10:11) !Ca++
    ionfitmapfwd(:,204) = fitmapfwd(12:13) !NH4+
    ionfitmapfwd(:,223) = fitmapfwd(14:15) !Mg++
    ionfitmapfwd(:,242) = fitmapfwd(16:17) !Cl-
    ionfitmapfwd(:,243) = fitmapfwd(18:19) !Br-
    ionfitmapfwd(:,245) = fitmapfwd(20:21) !NO3-
    ionfitmapfwd(:,261) = fitmapfwd(22:23) !SO4--
    ionfitmapfwd(:,248) = fitmapfwd(24:25) !HSO4-
    ionfitmapfwd(:,244) = fitmapfwd(26:27) !I-

    !==========Catan fitpars=====================
    nuecat = 1.0D0
    nuean = 1.0D0

    catanfitmapfwd(1,205,242) = fitmapfwd(28) !HCl
    catanfitmapfwd(1,201,242) = fitmapfwd(29) !LiCl
    catanfitmapfwd(1,203,242) = fitmapfwd(30) !KCl
    catanfitmapfwd(1,202,242) = fitmapfwd(31) !NaCl
    catanfitmapfwd(1,221,242) = fitmapfwd(32) !CaCl2
    nuean(221,242) = 2.0D0
    catanfitmapfwd(1,204,242) = fitmapfwd(33) !NH4Cl
    catanfitmapfwd(1,223,242) = fitmapfwd(34) !MgCl2
    nuean(223,242) = 2.0D0
    catanfitmapfwd(1,205,243) = fitmapfwd(35) !HBr
    catanfitmapfwd(1,201,243) = fitmapfwd(36) !LiBr
    catanfitmapfwd(1,203,243) = fitmapfwd(37) !KBr
    catanfitmapfwd(1,202,243) = fitmapfwd(38) !NaBr
    catanfitmapfwd(1,221,243) = fitmapfwd(39) !CaBr2
    nuean(221,243) = 2.0D0
    catanfitmapfwd(1,204,243) = fitmapfwd(40) !NH4Br
    catanfitmapfwd(1,223,243) = fitmapfwd(41) !MgBr2
    nuean(223,243) = 2.0D0
    catanfitmapfwd(1,205,245) = fitmapfwd(42) !HNO3
    catanfitmapfwd(1,201,245) = fitmapfwd(43) !LiNO3
    catanfitmapfwd(1,203,245) = fitmapfwd(44) !KNO3
    catanfitmapfwd(1,202,245) = fitmapfwd(45) !NaNO3
    catanfitmapfwd(1,221,245) = fitmapfwd(46) !CaNO32
    nuean(221,245) = 2.0D0
    catanfitmapfwd(1,204,245) = fitmapfwd(47) !NH4NO3
    catanfitmapfwd(1,223,245) = fitmapfwd(48) !MgNO32
    nuean(223,245) = 2.0D0
    catanfitmapfwd(1,205,261) = fitmapfwd(49) !H2SO4 (H+ and SO4--)
    nuecat(205,261) = 2.0D0
    catanfitmapfwd(1,201,261) = fitmapfwd(50) !Li2SO4
    nuecat(201,261) = 2.0D0
    catanfitmapfwd(1,203,261) = fitmapfwd(51) !K2SO4
    nuecat(203,261) = 2.0D0
    catanfitmapfwd(1,202,261) = fitmapfwd(52) !Na2SO4
    nuecat(202,261) = 2.0D0
    catanfitmapfwd(1,221,261) = fitmapfwd(53) !CaSO4
    catanfitmapfwd(1,204,261) = fitmapfwd(54) !NH42SO4
    nuecat(204,261) = 2.0D0
    catanfitmapfwd(1,223,261) = fitmapfwd(55) !MgSO4
    catanfitmapfwd(1,205,248) = fitmapfwd(56) !H2SO4 (H+ and HSO4-)
    catanfitmapfwd(1,201,248) = fitmapfwd(57) !LiHSO4
    catanfitmapfwd(1,203,248) = fitmapfwd(58) !KHSO4
    catanfitmapfwd(1,202,248) = fitmapfwd(59) !NaHSO4
    catanfitmapfwd(1,221,248) = fitmapfwd(60) !Ca(HSO4)2
    nuean(221,248) = 2.0D0
    catanfitmapfwd(1,204,248) = fitmapfwd(61) !NH4HSO4
    catanfitmapfwd(1,223,248) = fitmapfwd(62) !Mg(HSO4)2
    nuean(223,248) = 2.0D0
    catanfitmapfwd(1,205,244) = fitmapfwd(63) !HI
    catanfitmapfwd(1,201,244) = fitmapfwd(64) !LiI
    catanfitmapfwd(1,203,244) = fitmapfwd(65) !KI
    catanfitmapfwd(1,202,244) = fitmapfwd(66) !NaI
    catanfitmapfwd(1,221,244) = fitmapfwd(67) !CaI2
    nuean(221,244) = 2.0D0
    catanfitmapfwd(1,204,244) = fitmapfwd(68) !NH4I
    catanfitmapfwd(1,223,244) = fitmapfwd(69) !MgI2
    nuean(223,244) = 2.0D0

    END SUBROUTINE MapIons2Fitpars3
    !-------------------------------------------------------------------------------------

END MODULE ModViscEyring