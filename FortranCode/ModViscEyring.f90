!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   This module contains subroutines to compute the viscosity of aqueous electrolyte   *
!*   solutions with treatment of single-ion and cation--anion pair contributions and    *
!*   associated parameters; see paper by Lilek and Zuend (2022, ACP).                   *
!*                                                                                      *
!*   Options for viscosity calculation modes for organic--inorganic mixtures are set in *
!*   the module declaration part below.                                                 *
!*                                                                                      *
!*   :: Authors & Copyright ::                                                          *
!*   Joseph Lilek, Andi Zuend,                                                          *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2020                                                            *
!*   -> latest changes: 2022-02-11                                                      *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  SUBROUTINE AqueousElecViscosity                                                 *
!*   -  SUBROUTINE WaterMolefracCorrection                                              *
!*   -  SUBROUTINE GoldsackViscEqn                                                      *
!*   -  SUBROUTINE MapIons2Fitpars3                                                     *
!*                                                                                      *
!****************************************************************************************
MODULE ModViscEyring

USE ModAIOMFACvar, ONLY : galrln, gamrln, gasrln, gclrln, gcmrln, gcsrln, gnlrln, &
    & gnmrln, gnsrln, lnactcoeff_a, lnactcoeff_c, SMA, SMC, solvmixcorrMRc, Tmolal, &
    & TmolalSolvmix
USE ModSystemProp, ONLY : errorflag_clist, Ianion, Ication, Mmass, Nanion, Ncation, &
    & nd, solvmixrefnd, topsubno
USE ModSRparam, ONLY : SR_RR

IMPLICIT NONE

INTEGER(4),PRIVATE :: i
INTEGER(4),DIMENSION(1:105),PARAMETER,PRIVATE :: fitmapfwd = [(i, i = 1,105)]
INTEGER(4),DIMENSION(2,201:topsubno),PRIVATE :: ionfitmapfwd          
INTEGER(4),DIMENSION(2,201:240,241:topsubno),PRIVATE :: catanfitmapfwd
INTEGER(4),DIMENSION(:),ALLOCATABLE,PRIVATE :: ions  
REAL(8),DIMENSION(201:topsubno),PRIVATE :: Vion, Zion                 
REAL(8),DIMENSION(201:240,241:topsubno),PRIVATE :: nuecat, nuean      
REAL(8),PUBLIC :: delGstar_by_RT_w, ionicstrengthsave
!LOGICAL(4),PRIVATE :: gel_eqn1, gel_eqn2 

!private variables with initial values that may change:
LOGICAL(4),PRIVATE :: fitmap_complete = .true.                  !fitmap for viscosity model parameter init; needs to be .true. here; changes after first use
LOGICAL(4),PRIVATE :: calc_gel = .false.                        !false by default; currently gel effects are in development
!****
!**** set a few parameters viscosity model calculation options:
LOGICAL(4),PARAMETER,PUBLIC :: aquelec = .true.                 !true is default, false sets 'aquorg' mixing as mixing rule (except when ZSRvisc_on is true);
LOGICAL(4),PARAMETER,PUBLIC :: ZSRvisc_on = .false.             !false is default, true activates ZSR viscosity subroutine for org--inorg mixing.
LOGICAL(4),PARAMETER,PRIVATE :: calc_catan = .true.
LOGICAL(4),PARAMETER,PRIVATE :: newcatanPrime = .true.          !if true, use the non-normalized tau'. If false (default) use normalized tau.
!****
!****
REAL(8),DIMENSION(1:105),PARAMETER,PUBLIC :: fitpar = [ &
! 28 Dec 2021 (Lilek and Zuend revised, updated, DEFAULT)
& 1.67983D+00,  1.73760D-03,  8.20090D+00,  1.57456D-01,  9.57170D+00,  4.35501D-01,  3.52182D+00,  2.32051D-01,  8.70806D+00,  4.90251D-02, &
& 1.22486D+01,  2.20685D-01,  4.39636D+00,  3.25016D-02,  2.97547D+01,  2.58001D-02,  3.83450D+00,  9.70421D-03,  9.00033D-01,  1.42443D+00, &
& 1.63357D+00,  1.01039D-11,  1.84540D+01,  2.14728D-02,  7.64053D+00,  2.13062D-01,  1.92583D-02,  1.01039D-11,  2.42196D+01,  1.17449D+00, &
& 1.51967D+01,  1.01039D-11,  1.52129D+01,  1.60020D+01,  3.39381D+01,  3.26414D-01,  1.44456D+00,  3.39242D-01,  4.59481D-01,  4.69448D+00, &
& 1.01039D-11,  2.47812D+00,  1.43264D+00,  1.70274D+00,  1.42802D+00,  1.55456D+00,  4.69448D+00,  1.01039D-11,  2.47812D+00,  1.01039D-11, &
& 1.82254D-02,  1.74929D+00,  1.80425D+00,  6.77700D+00,  5.67534D-01,  8.19343D-01,  2.53147D-01,  3.70756D+00,  6.87103D-02,  3.17849D-01, &
& 3.39729D+00,  1.01039D-11,  3.39729D+00,  2.50423D-01,  1.44456D+00,  3.39242D-01,  1.31901D+00,  4.69448D+00,  1.01039D-11,  2.47812D+00, &
& 1.10588D+00,  1.34544D+00,  1.71202D+00,  1.85803D+00,  4.69448D+00,  1.01039D-11,  2.47812D+00,  2.30706D+01,  3.80327D+00,  1.52573D+00, &
& 3.80327D+00,  3.39729D+00,  1.52573D+00,  3.39729D+00,  2.37260D+01,  1.39357D-02,  1.01039D-11,  1.39357D-02,  6.77700D+00,  1.01039D-11, &
& 8.19343D-01,  2.45867D+01,  1.70520D+00,  3.45063D-02,  1.01039D-11,  4.69448D+00,  3.45063D-02,  2.47812D+00,  3.40645D+00,  1.82254D-02, &
& 1.74929D+00,  1.80425D+00,  6.77700D+00,  5.67534D-01,  8.19343D-01 ]

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
        CALL GoldsackViscEqn(X_ion_organicfree, act_ion, ionicstrength1, xw, Vw, ln_etaw, ln_eta_aquelec, errorflag_clist)
    ELSE !aquorg. in GoldsackViscEqn, xw = 1 - sum(X_ion), so mole fractions of organics are effectively counted as mole fraction of water
        xw = 1.0D0 - SUM(X_ion(:))
        sumX = SUM(X(1:nneutral))
        Vcomp = 15.17D-06 * SUM((X(1:nneutral)/sumX) * RS_(1:nneutral))   !use X (full system mole fractions to calc solvent volume)
        CALL GoldsackViscEqn(X_ion, act_ion, ionicstrength1, xw, Vcomp, ln_etaw, ln_eta_aquelec, errorflag_clist) !
    ENDIF
    
    DEALLOCATE( ions )

    END SUBROUTINE AqueousElecViscosity
    !-------------------------------------------------------------------------------------
    
    
    !-------------------------------------------------------------------------------------
    PURE SUBROUTINE GoldsackViscEqn(xin_, actin_, ionicstrength_, xw, Vw, ln_etaw, ln_etacalc, errorflag_clist)
    
    USE ModSystemProp, ONLY : Ianion, Ication, Nanion, Ncation, topsubno

    IMPLICIT NONE
    
    REAL(8),DIMENSION(201:topsubno),INTENT(IN) :: xin_, actin_
    REAL(8),INTENT(IN) :: ionicstrength_, ln_etaw, Vw, xw
    REAL(8),INTENT(OUT) :: ln_etacalc
    LOGICAL(4),DIMENSION(:),INTENT(INOUT) :: errorflag_clist
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
    REAL(8),PARAMETER :: R = 8.3144598D0                    !J K^-1 mol^-1 ideal gas constant
    REAL(8),PARAMETER :: hPlanck = 6.62607015D-34           ![J s]        Planck's constant
    REAL(8),PARAMETER :: NA = 6.02214076D23                 ![#/mol]      Avogadro's constant
    REAL(8),PARAMETER :: MmassH2O = 1.801528D-02            !kg mol^-1
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
                    errorflag_clist(18) = .true.
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
    
    delGstar_by_RT_w = xw *(ln_etaw +LOG(Vw) -ln_hNA)
    !(4) deltaG_water + delGstar_ions + delGstar_catan
    delGstar_over_RT = &
        &   delGstar_ions &         !all cation & anion contributions multiplied by xin
        & + delGstar_catans &       !ion strength term multiplied by tau weighting
        & + delGstar_by_RT_w        !water contribution multiplied by xw

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
    AquorgMmass = SUM( (X_(1:nneutral)/SUM(X_(1:nneutral))) * Mmass(1:nneutral) )
    TotalAvgMmass = SUM(X_ * Mmass )         
        
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

    SUBROUTINE MapIons2Fitpars3()

    IMPLICIT NONE
    !......................
    !initialize arrays (fitmapfwd is simply integers 1:105)
    ionfitmapfwd = 0
    catanfitmapfwd = 0
    
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
    
    !Add ions CO3--, HCO3-, OH-, IO3-
    ionfitmapfwd(:,262) = fitmapfwd(28:29) !CO3--
    ionfitmapfwd(:,250) = fitmapfwd(30:31) !HCO3-
    ionfitmapfwd(:,247) = fitmapfwd(32:33) !OH-
    ionfitmapfwd(:,246) = fitmapfwd(34:35) !IO3-

    !==========Catan fitpars=====================
    nuecat = 1.0D0
    nuean = 1.0D0
    
    catanfitmapfwd(1,205,242) = fitmapfwd(36) !HCl
    catanfitmapfwd(1,201,242) = fitmapfwd(37) !LiCl
    catanfitmapfwd(1,203,242) = fitmapfwd(38) !KCl
    catanfitmapfwd(1,202,242) = fitmapfwd(39) !NaCl
    catanfitmapfwd(1,221,242) = fitmapfwd(40) !CaCl2
    nuean(221,242) = 2.0D0
    catanfitmapfwd(1,204,242) = fitmapfwd(41) !NH4Cl
    catanfitmapfwd(1,223,242) = fitmapfwd(42) !MgCl2
    nuean(223,242) = 2.0D0
    !
    catanfitmapfwd(1,205,243) = fitmapfwd(43) !HBr
    catanfitmapfwd(1,201,243) = fitmapfwd(44) !LiBr
    catanfitmapfwd(1,203,243) = fitmapfwd(45) !KBr
    catanfitmapfwd(1,202,243) = fitmapfwd(46) !NaBr
    catanfitmapfwd(1,221,243) = fitmapfwd(47) !CaBr2
    nuean(221,243) = 2.0D0
    catanfitmapfwd(1,204,243) = fitmapfwd(48) !NH4Br
    catanfitmapfwd(1,223,243) = fitmapfwd(49) !MgBr2
    nuean(223,243) = 2.0D0
    !
    catanfitmapfwd(1,205,245) = fitmapfwd(50) !HNO3
    catanfitmapfwd(1,201,245) = fitmapfwd(51) !LiNO3
    catanfitmapfwd(1,203,245) = fitmapfwd(52) !KNO3
    catanfitmapfwd(1,202,245) = fitmapfwd(53) !NaNO3
    catanfitmapfwd(1,221,245) = fitmapfwd(54) !CaNO32
    nuean(221,245) = 2.0D0
    catanfitmapfwd(1,204,245) = fitmapfwd(55) !NH4NO3
    catanfitmapfwd(1,223,245) = fitmapfwd(56) !MgNO32
    nuean(223,245) = 2.0D0
    !
    catanfitmapfwd(1,205,261) = fitmapfwd(57) !H2SO4 (H+ and SO4--)
    nuecat(205,261) = 2.0D0
    catanfitmapfwd(1,201,261) = fitmapfwd(58) !Li2SO4
    nuecat(201,261) = 2.0D0
    catanfitmapfwd(1,203,261) = fitmapfwd(59) !K2SO4
    nuecat(203,261) = 2.0D0
    catanfitmapfwd(1,202,261) = fitmapfwd(60) !Na2SO4
    nuecat(202,261) = 2.0D0
    catanfitmapfwd(1,221,261) = fitmapfwd(61) !CaSO4
    catanfitmapfwd(1,204,261) = fitmapfwd(62) !NH42SO4
    nuecat(204,261) = 2.0D0
    !
    catanfitmapfwd(1,223,261) = fitmapfwd(63) !MgSO4
    catanfitmapfwd(1,205,248) = fitmapfwd(64) !H2SO4 (H+ and HSO4-)
    catanfitmapfwd(1,201,248) = fitmapfwd(65) !LiHSO4
    catanfitmapfwd(1,203,248) = fitmapfwd(66) !KHSO4
    catanfitmapfwd(1,202,248) = fitmapfwd(67) !NaHSO4
    catanfitmapfwd(1,221,248) = fitmapfwd(68) !Ca(HSO4)2
    nuean(221,248) = 2.0D0
    catanfitmapfwd(1,204,248) = fitmapfwd(69) !NH4HSO4
    catanfitmapfwd(1,223,248) = fitmapfwd(70) !Mg(HSO4)2
    nuean(223,248) = 2.0D0
    !
    catanfitmapfwd(1,205,244) = fitmapfwd(71) !HI
    catanfitmapfwd(1,201,244) = fitmapfwd(72) !LiI
    catanfitmapfwd(1,203,244) = fitmapfwd(73) !KI
    catanfitmapfwd(1,202,244) = fitmapfwd(74) !NaI
    catanfitmapfwd(1,221,244) = fitmapfwd(75) !CaI2
    nuean(221,244) = 2.0D0
    catanfitmapfwd(1,204,244) = fitmapfwd(76) !NH4I
    catanfitmapfwd(1,223,244) = fitmapfwd(77) !MgI2
    nuean(223,244) = 2.0D0
    !
    catanfitmapfwd(1,205,262) = fitmapfwd(78) !H2CO3 (H+ and CO3--)
    nuecat(205,262) = 2.0D0
    catanfitmapfwd(1,201,262) = fitmapfwd(79) !Li2CO3
    nuecat(201,262) = 2.0D0
    catanfitmapfwd(1,203,262) = fitmapfwd(80) !K2CO3
    nuecat(203,262) = 2.0D0
    catanfitmapfwd(1,202,262) = fitmapfwd(81) !Na2CO3
    nuecat(202,262) = 2.0D0
    catanfitmapfwd(1,221,262) = fitmapfwd(82) !CaCO3
    catanfitmapfwd(1,204,262) = fitmapfwd(83) !NH42CO3
    nuecat(204,262) = 2.0D0
    catanfitmapfwd(1,223,262) = fitmapfwd(84) !MgCO3
    !
    catanfitmapfwd(1,205,250) = fitmapfwd(85) !H2CO3 (H+ and HCO3-)
    catanfitmapfwd(1,201,250) = fitmapfwd(86) !LiHCO3
    catanfitmapfwd(1,203,250) = fitmapfwd(87) !KHCO3
    catanfitmapfwd(1,202,250) = fitmapfwd(88) !NaHCO3
    catanfitmapfwd(1,221,250) = fitmapfwd(89) !Ca(HCO3)2
    nuean(221,250) = 2.0D0
    catanfitmapfwd(1,204,250) = fitmapfwd(90) !NH4HCO3
    catanfitmapfwd(1,223,250) = fitmapfwd(91) !Mg(HCO3)2
    nuean(223,250) = 2.0D0
    !
    catanfitmapfwd(1,205,247) = fitmapfwd(92) !H2O (H+ and OH-)
    catanfitmapfwd(1,201,247) = fitmapfwd(93) !LiOH
    catanfitmapfwd(1,203,247) = fitmapfwd(94) !KOH
    catanfitmapfwd(1,202,247) = fitmapfwd(95) !NaOH
    catanfitmapfwd(1,221,247) = fitmapfwd(96) !Ca(OH)2
    nuean(221,247) = 2.0D0
    catanfitmapfwd(1,204,247) = fitmapfwd(97) !NH4OH
    catanfitmapfwd(1,223,247) = fitmapfwd(98) !Mg(OH)2
    nuean(223,247) = 2.0D0
    !
    catanfitmapfwd(1,205,246) = fitmapfwd(99) !HIO3 (H+ and IO3-)
    catanfitmapfwd(1,201,246) = fitmapfwd(100) !LiIO3
    catanfitmapfwd(1,203,246) = fitmapfwd(101) !KIO3
    catanfitmapfwd(1,202,246) = fitmapfwd(102) !NaIO3
    catanfitmapfwd(1,221,246) = fitmapfwd(103) !Ca(IO3)2
    nuean(221,246) = 2.0D0
    catanfitmapfwd(1,204,246) = fitmapfwd(104) !NH4IO3
    catanfitmapfwd(1,223,246) = fitmapfwd(105) !Mg(IO3)2
    nuean(223,246) = 2.0D0

    END SUBROUTINE MapIons2Fitpars3
    !-------------------------------------------------------------------------------------
 

END MODULE ModViscEyring