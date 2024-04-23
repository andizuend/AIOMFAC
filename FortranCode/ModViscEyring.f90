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
!*   -> latest changes: 2022-02-11                                                      *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  subroutine AqueousElecViscosity                                                 *
!*   -  subroutine WaterMolefracCorrection                                              *
!*   -  subroutine GoldsackViscEqn                                                      *
!*   -  subroutine MapIons2Fitpars3                                                     *
!*   -  subroutine ViscGelContribution                                                  *
!*                                                                                      *
!****************************************************************************************
module ModViscEyring

use Mod_kind_param, only : wp
use ModAIOMFACvar, only : galrln, gamrln, gasrln, gclrln, gcmrln, gcsrln, gnlrln, &
    & gnmrln, gnsrln, lnactcoeff_a, lnactcoeff_c, SMA, SMC, solvmixcorrMRc, Tmolal, &
    & TmolalSolvmix
use ModSystemProp, only : errorflag_clist, Ianion, Ication, Mmass, Nanion, Ncation, &
    & nd, solvmixrefnd, topsubno
use ModSRparam, only : SR_RR

implicit none

integer,private :: i
integer,dimension(1:105),parameter,private :: fitmapfwd = [(i, i = 1,105)]
integer,dimension(2,201:topsubno),private :: ionfitmapfwd          
integer,dimension(2,201:240,241:topsubno),private :: catanfitmapfwd
integer,dimension(:),allocatable,private :: ions  
real(wp),dimension(201:topsubno),private :: Vion, Zion                 
real(wp),dimension(201:240,241:topsubno),private :: nuecat, nuean      
real(wp),public :: delGstar_by_RT_w, ionicstrengthsave
!logical,private :: gel_eqn1, gel_eqn2 

!private variables with initial values that may change:
logical,private :: fitmap_complete = .true.                  !fitmap for viscosity model parameter init; needs to be .true. here; changed after first use
logical,private :: calc_gel = .false.                        !false by default; currently gel effects are in development
!****
!**** set a few parameters for viscosity model calculation options:
logical,parameter,public :: aquelec = .true.                 !true is default, false sets 'aquorg' mixing as mixing rule (except when ZSRvisc_on is true);
logical,parameter,public :: ZSRvisc_on = .false.             !false is default, true activates ZSR viscosity subroutine for org--inorg viscosity mixing
logical,parameter,private :: calc_catan = .true.
logical,parameter,private :: newcatanPrime = .true.          !if true, use the non-normalized tau'. If false use normalized tau.
!****
!****
real(wp),dimension(1:105),parameter,public :: fitpar = [ &

! 28 Dec 2021 (Lilek and Zuend revised, updated, default)
& 1.67983E+00_wp,  1.73760E-03_wp,  8.20090E+00_wp,  1.57456E-01_wp,  9.57170E+00_wp,  4.35501E-01_wp,  3.52182E+00_wp,  2.32051E-01_wp,  8.70806E+00_wp,  4.90251E-02_wp, &
& 1.22486E+01_wp,  2.20685E-01_wp,  4.39636E+00_wp,  3.25016E-02_wp,  2.97547E+01_wp,  2.58001E-02_wp,  3.83450E+00_wp,  9.70421E-03_wp,  9.00033E-01_wp,  1.42443E+00_wp, &
& 1.63357E+00_wp,  1.01039E-11_wp,  1.84540E+01_wp,  2.14728E-02_wp,  7.64053E+00_wp,  2.13062E-01_wp,  1.92583E-02_wp,  1.01039E-11_wp,  2.42196E+01_wp,  1.17449E+00_wp, &
& 1.51967E+01_wp,  1.01039E-11_wp,  1.52129E+01_wp,  1.60020E+01_wp,  3.39381E+01_wp,  3.26414E-01_wp,  1.44456E+00_wp,  3.39242E-01_wp,  4.59481E-01_wp,  4.69448E+00_wp, &
& 1.01039E-11_wp,  2.47812E+00_wp,  1.43264E+00_wp,  1.70274E+00_wp,  1.42802E+00_wp,  1.55456E+00_wp,  4.69448E+00_wp,  1.01039E-11_wp,  2.47812E+00_wp,  1.01039E-11_wp, &
& 1.82254E-02_wp,  1.74929E+00_wp,  1.80425E+00_wp,  6.77700E+00_wp,  5.67534E-01_wp,  8.19343E-01_wp,  2.53147E-01_wp,  3.70756E+00_wp,  6.87103E-02_wp,  3.17849E-01_wp, &
& 3.39729E+00_wp,  1.01039E-11_wp,  3.39729E+00_wp,  2.50423E-01_wp,  1.44456E+00_wp,  3.39242E-01_wp,  1.31901E+00_wp,  4.69448E+00_wp,  1.01039E-11_wp,  2.47812E+00_wp, &
& 1.10588E+00_wp,  1.34544E+00_wp,  1.71202E+00_wp,  1.85803E+00_wp,  4.69448E+00_wp,  1.01039E-11_wp,  2.47812E+00_wp,  2.30706E+01_wp,  3.80327E+00_wp,  1.52573E+00_wp, &
& 3.80327E+00_wp,  3.39729E+00_wp,  1.52573E+00_wp,  3.39729E+00_wp,  2.37260E+01_wp,  1.39357E-02_wp,  1.01039E-11_wp,  1.39357E-02_wp,  6.77700E+00_wp,  1.01039E-11_wp, &
& 8.19343E-01_wp,  2.45867E+01_wp,  1.70520E+00_wp,  3.45063E-02_wp,  1.01039E-11_wp,  4.69448E+00_wp,  3.45063E-02_wp,  2.47812E+00_wp,  3.40645E+00_wp,  1.82254E-02_wp, &
& 1.74929E+00_wp,  1.80425E+00_wp,  6.77700E+00_wp,  5.67534E-01_wp,  8.19343E-01_wp ]

!public procedures:
public :: AqueousElecViscosity, WaterMolefracCorrection
private     !default setting for procedures and variables

!make certain variables of this module threadprivate for use in parallel execution with openMP:
!$OMP THREADPRIVATE( ions, delGstar_by_RT_w, Zion, ionicstrengthsave)

!=====================================
    contains
!=====================================   
    
    !-------------------------------------------------------------------------------------
    subroutine AqueousElecViscosity(X_, lnGaSR_, RS_, ln_etaw, ln_eta_aquelec)
    
    use ModAIOMFACvar, only : galrln, gamrln, gasrln, gclrln, gcmrln, &
        & gcsrln, SMC, SMA, lnactcoeff_c, lnactcoeff_a, solvmixcorrMRc, solvmixcorrMRa, Tmolal, &
        & TmolalSolvmix, actcoeff_c, actcoeff_a, wtf, X
    use ModSystemProp, only : Ianion, Ication, Nanion, Ncation, topsubno, solvmixrefnd, nneutral, &
        & anionZ, cationZ    
    use ModCompScaleConversion, only : MassFrac2MoleFracMolality, zSolution2SpeciesMolality
    use ModNumericalTransformations, only : safe_exp

    implicit none

    !interface variables:
    real(wp),dimension(:),intent(in) :: X_, lnGaSR_, RS_
    real(wp),intent(in) :: ln_etaw
    real(wp),intent(out) :: ln_eta_aquelec
    !local variables:
    integer :: I, NKNpNcat, nnp1, allocstat
    integer :: iion, icat, ian
    real(wp),parameter :: lnhuge = 0.49_wp*log(huge(1.0_wp))
    real(wp),parameter :: logval_threshold = 0.4_wp*lnhuge
    real(wp),parameter :: dtiny = epsilon(1.0_wp) 
    real(wp),parameter :: Vw = 1.39564E-5_wp        ![m^3/mol]  Vw for water molec. volume;
    real(wp) :: ionicstrength1, sumXion, sumX, ionicstrengthfactor, Vcomp, xw
    real(wp),dimension(201:topsubno) :: act_ion, X_ion, X_ion_organicfree
    !................................................
    
    !tasks that only need to be done once during program run:
    if (fitmap_complete) then
        fitmap_complete = .false.
        Vion = 15.17E-06_wp * SR_RR(201:topsubno)
        call MapIons2Fitpars3()
    endif
    
    !initialize:    
    Zion = 0.0_wp
    X_ion = 0.0_wp
    X_ion_organicfree = 0.0_wp
    act_ion = 0.0_wp
    
    NKNpNcat = nneutral + Ncation
    nnp1 = nneutral + 1

    gcsrln(1:Ncation) = lnGaSR_(nnp1:NKNpNcat)
    gasrln(1:Nanion)  = lnGaSR_(NKNpNcat+1:NKNpNcat+Nanion)
    
    if (aquelec .AND. wtf(1) > 0.0_wp) then    
        ionicstrengthfactor = sum(wtf(1:nneutral))/wtf(1)
    else
        ionicstrengthfactor = 1.0_wp
    endif
       
    !the activity coefficient contributions for LR and MR parts come via module access (use ModAIOMFACvar, only : ...)
    !initialize the output arrays:
    lnactcoeff_c = -9999.9_wp
    lnactcoeff_a = -9999.9_wp
    actcoeff_c = 0.0_wp
    actcoeff_a = 0.0_wp
    
    !ln of the activity coefficient for the cations:
    do I = 1,Ncation
        if (SMC(I) > dtiny) then
            lnactcoeff_c(I) = gcmrln(I) +gcsrln(I) +gclrln(I) -Tmolal  !this term converts to the molality scale (basis/scale conversion)
            if (solvmixrefnd) then      !correction terms for MR and SR part, if  reference solution is the solvent mixture
                lnactcoeff_c(I) = lnactcoeff_c(I) +solvmixcorrMRc(I) +Tmolal -TmolalSolvmix    
            endif
            actcoeff_c(I) = safe_exp(lnactcoeff_c(I), logval_threshold)
        endif
    enddo
    
    !ln of the activity coefficient for the anions:
    do I = 1,Nanion 
        if (SMA(I) > dtiny) then
            lnactcoeff_a(I) = gamrln(I) +gasrln(I) +galrln(I) -Tmolal    !this term converts to the molality scale (basis/scale conversion)
            if (solvmixrefnd) then
                lnactcoeff_a(I) = lnactcoeff_a(I) +solvmixcorrMRa(I) +Tmolal -TmolalSolvmix 
            endif
            actcoeff_a(I) = safe_exp(lnactcoeff_a(I), logval_threshold)
        endif
    enddo    
    
    allocate( ions(Ncation+Nanion), STAT=allocstat)
    icat = 0
    ian = 0
    
    !define mole fraction, activity of cations
    do I = 1,Ncation
       !save activity and mole fraction of this ion in an array by actual AIOMFAC ion index:
        iion = Ication(I)   !iion is ion index number (201:topsubno)
        ions(I) = iion
        Zion(iion) = cationZ(I)
        act_ion(iion) = actcoeff_c(I)*SMC(I) * ionicstrengthfactor !ionicstrengthfactor to change solvent in molality to H2O instead of all neutral components
        X_ion(iion) = X_(nneutral+I)                            !these ion mole fractions are based on the full system (water, orgs, cation, anion)        
    enddo

    !define mole fraction, activity of anions
    do I = 1,Nanion 
        iion = Ianion(I)
        ions(Ncation+I) = Ianion(I)
        Zion(iion) = abs(anionZ(I))
        act_ion(iion) = actcoeff_a(I)*SMA(I)*ionicstrengthfactor
        X_ion(iion) = X_(nneutral+Ncation+I)
    enddo
    
    !ionicstrengthfactor  is used to change solvent to H2O instead of all neutral components:
    ionicstrength1 = 0.5_wp*ionicstrengthfactor*( sum(SMC(1:Ncation)*cationZ(1:Ncation)**2) + sum(SMA(1:Nanion)*anionZ(1:Nanion)**2) )
    ionicstrengthsave = ionicstrength1
    
    if (aquelec) then
        sumXion = sum(X_ion(:))
        X_ion_organicfree = X_ion(:) / (X_(1) + sumXion)
        xw = 1.0_wp - sum(X_ion_organicfree(:))
        call GoldsackViscEqn(X_ion_organicfree, act_ion, ionicstrength1, xw, Vw, ln_etaw, ln_eta_aquelec, errorflag_clist)
    else !aquorg. in GoldsackViscEqn, xw = 1 - sum(X_ion), so mole fractions of organics are effectively counted as mole fraction of water
        xw = 1.0_wp - sum(X_ion(:))
        sumX = sum(X(1:nneutral))
        Vcomp = 15.17E-06_wp * sum((X(1:nneutral)/sumX) * RS_(1:nneutral))   !use X (full system mole fractions to calc solvent volume)
        call GoldsackViscEqn(X_ion, act_ion, ionicstrength1, xw, Vcomp, ln_etaw, ln_eta_aquelec, errorflag_clist) !
    endif
    
    deallocate( ions )

    end subroutine AqueousElecViscosity
    !-------------------------------------------------------------------------------------
    
    
    !-------------------------------------------------------------------------------------
    pure subroutine GoldsackViscEqn(xin_, actin_, ionicstrength_, xw, Vw, ln_etaw, ln_etacalc, errorflag_clist)
    
    use ModSystemProp, only : Ianion, Ication, Nanion, Ncation, topsubno

    implicit none
    
    real(wp),dimension(201:topsubno),intent(in) :: xin_, actin_
    real(wp),intent(in) :: ionicstrength_, ln_etaw, Vw, xw
    real(wp),intent(out) :: ln_etacalc
    logical,dimension(:),intent(inout) :: errorflag_clist
    !local variables:
    integer :: iion, icat, ian
    ! Goldsack & Franchetto vars and parameters:
    integer :: cation, anion, numion                            !which value they refer to in the AIOMFAC ion numbering system (201:topsubno)
    real(wp) :: delGstar_over_RT, delGstar_by_RT_w, V, Garg_both, Garg_ion
    real(wp) :: delGstar_ions, delGstar_catans, sumtauPrime, sumXz 
    real(wp),dimension(201:topsubno) :: chargefrac, delGstar_ion
    real(wp),dimension(201:240,241:topsubno) :: delGstar_catan
    real(wp),dimension(201:240,241:topsubno) :: tauPrime, tau
    real(wp),parameter :: deps = epsilon(1.0_wp)
    real(wp),parameter :: R = 8.3144598_wp                      ![J K^-1 mol^-1]  ideal gas constant
    real(wp),parameter :: hPlanck = 6.62607015E-34_wp           ![J s]            Planck's constant
    real(wp),parameter :: NA = 6.02214076E+23_wp                ![#/mol]          Avogadro's constant
    real(wp),parameter :: MmassH2O = 1.801528E-02_wp            ![kg mol^-1]
    real(wp),parameter :: hNA = hPlanck*NA                      ![J s #/mol]
    real(wp),parameter :: ln_hNA = log(hNA)                     ![-] (normalized)
    !............................................
    
    !Note: xw = 1.0_wp - sum(xin_(:))
    V = xw*Vw + fitpar(1)*sum( xin_(:) * Vion(:) )          ![m^3/mol]  effective average hole volume; use fit parameter to affect ion vol. contributions (but not water);
        
    !(1) delGstar_ions (Garg_ion for all individual ions, summed)
    delGstar_ion = 0.0_wp
    do iion = 1,Ncation + Nanion
        numion = ions(iion)
        if (ionfitmapfwd(1,numion) /= 0) then
            Garg_ion = actin_(numion)
            if (Garg_ion < deps) then
                Garg_ion = 1.0_wp
            endif
            delGstar_ion(numion) = xin_(numion) * (fitpar(ionfitmapfwd(1,numion)) * log(Garg_ion) + fitpar(ionfitmapfwd(2,numion)))
        else
            delGstar_ion(numion) = 0.0_wp
        endif
    enddo
    delGstar_ions = sum( delGstar_ion(:) )

    !(2) deltaG_catan (Garg_both of all cation--anion pairs, summed)
    if (calc_catan) then
        !initialize
        delGstar_catan = 0.0_wp
        tauPrime = 0.0_wp
        tau = 0.0_wp
        sumtauPrime = 0.0_wp
        chargefrac = 0.0_wp
        !use ion molalities
        sumXz = sum(xin_(201:240)*abs(Zion(201:240)))
        do icat = 1,Ncation!ncat
            cation = Ication(icat)
            do ian = 1,Nanion
                anion = Ianion(ian)
                !only one fitpar per catan pair:
                if (catanfitmapfwd(1,cation,anion) /= 0 ) then
                    Garg_both = fitpar(catanfitmapfwd(1,cation,anion)) * sqrt(ionicstrength_)
                    if (Garg_both < deps) then
                        Garg_both = 1.0_wp
                    endif
                    !implementation of new catan treatment (17 June 2021):
                    chargefrac(anion) = (xin_(anion)*abs(Zion(anion))) / sumXz
                    tauPrime(cation,anion) = xin_(cation)/nuecat(cation,anion) * chargefrac(anion)
                    delGstar_catan(cation,anion) = Garg_both
                else
                    delGstar_catan(cation,anion) = 0.0_wp
                    errorflag_clist(18) = .true.
                endif
            enddo
        enddo
        sumtauPrime = sum(tauPrime(:,:))
        if (newcatanPrime .OR. sumtauPrime < deps) then
            where (delGstar_catan(:,:) > 0.0_wp)
                delGstar_catan(:,:) = delGstar_catan(:,:)*tauPrime(:,:)
            end where
        else
            tau = tauPrime/sumtauPrime
            where (delGstar_catan(:,:) > 0.0_wp)
                delGstar_catan(:,:) = delGstar_catan(:,:)*tau(:,:)
            end where
        endif
        delGstar_catans = sum(delGstar_catan(:,:))
    else
        delGstar_catans = 0.0_wp
    endif !calc_catan
    
    delGstar_by_RT_w = xw *(ln_etaw +log(Vw) -ln_hNA)
    !(4) deltaG_water + delGstar_ions + delGstar_catan
    delGstar_over_RT = &
        &   delGstar_ions &         !all cation & anion contributions multiplied by xin
        & + delGstar_catans &       !ionic strength term multiplied by tau weighting
        & + delGstar_by_RT_w        !water contribution multiplied by xw

    ln_etacalc = ln_hNA -log(V) + delGstar_over_RT      !normalized ln(eta/[Pa s])

    end subroutine GoldsackViscEqn
    !-------------------------------------------------------------------------------------

    
    !-------------------------------------------------------------------------------------
    pure subroutine WaterMolefracCorrection(X_, Xnew)
    
    use ModAIOMFACvar, only : wtf
    use ModSystemProp, only : nneutral, Mmass
        
    implicit none
        
    !interface vars
    real(wp),dimension(:),intent(in) :: X_
    real(wp),dimension(:),intent(out) :: Xnew 
    !local vars
    real(wp) :: Xaquorg, AquorgMmass, TotalAvgMmass, wtfaquorg
    !................................
    
    !(1) Calculate mass of water + organics:
    wtfaquorg = sum(wtf(1:nneutral))
    !(2) Calculate avg molar mass of the aquorg mixture:
    AquorgMmass = sum( (X_(1:nneutral)/sum(X_(1:nneutral))) * Mmass(1:nneutral) )
    TotalAvgMmass = sum(X_ * Mmass )         
        
    !(3) Calculate mole fraction amount for aquorg mixture:
    Xaquorg = (wtfaquorg / AquorgMmass) * TotalAvgMmass
        
    Xnew(1) = Xaquorg
    Xnew(2:nneutral) = 0.0_wp
    Xnew(nneutral+1:) = X_(nneutral+1:)
    Xnew = Xnew / sum(Xnew)
    
    end subroutine WaterMolefracCorrection
    !-------------------------------------------------------------------------------------
    
    
    !===============================================================================================    
    ! This subroutine maps ions to fit parameters.

    subroutine MapIons2Fitpars3()

    implicit none
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
    nuecat = 1.0_wp
    nuean = 1.0_wp
    
    catanfitmapfwd(1,205,242) = fitmapfwd(36) !HCl
    catanfitmapfwd(1,201,242) = fitmapfwd(37) !LiCl
    catanfitmapfwd(1,203,242) = fitmapfwd(38) !KCl
    catanfitmapfwd(1,202,242) = fitmapfwd(39) !NaCl
    catanfitmapfwd(1,221,242) = fitmapfwd(40) !CaCl2
    nuean(221,242) = 2.0_wp
    catanfitmapfwd(1,204,242) = fitmapfwd(41) !NH4Cl
    catanfitmapfwd(1,223,242) = fitmapfwd(42) !MgCl2
    nuean(223,242) = 2.0_wp
    !
    catanfitmapfwd(1,205,243) = fitmapfwd(43) !HBr
    catanfitmapfwd(1,201,243) = fitmapfwd(44) !LiBr
    catanfitmapfwd(1,203,243) = fitmapfwd(45) !KBr
    catanfitmapfwd(1,202,243) = fitmapfwd(46) !NaBr
    catanfitmapfwd(1,221,243) = fitmapfwd(47) !CaBr2
    nuean(221,243) = 2.0_wp
    catanfitmapfwd(1,204,243) = fitmapfwd(48) !NH4Br
    catanfitmapfwd(1,223,243) = fitmapfwd(49) !MgBr2
    nuean(223,243) = 2.0_wp
    !
    catanfitmapfwd(1,205,245) = fitmapfwd(50) !HNO3
    catanfitmapfwd(1,201,245) = fitmapfwd(51) !LiNO3
    catanfitmapfwd(1,203,245) = fitmapfwd(52) !KNO3
    catanfitmapfwd(1,202,245) = fitmapfwd(53) !NaNO3
    catanfitmapfwd(1,221,245) = fitmapfwd(54) !CaNO32
    nuean(221,245) = 2.0_wp
    catanfitmapfwd(1,204,245) = fitmapfwd(55) !NH4NO3
    catanfitmapfwd(1,223,245) = fitmapfwd(56) !MgNO32
    nuean(223,245) = 2.0_wp
    !
    catanfitmapfwd(1,205,261) = fitmapfwd(57) !H2SO4 (H+ and SO4--)
    nuecat(205,261) = 2.0_wp
    catanfitmapfwd(1,201,261) = fitmapfwd(58) !Li2SO4
    nuecat(201,261) = 2.0_wp
    catanfitmapfwd(1,203,261) = fitmapfwd(59) !K2SO4
    nuecat(203,261) = 2.0_wp
    catanfitmapfwd(1,202,261) = fitmapfwd(60) !Na2SO4
    nuecat(202,261) = 2.0_wp
    catanfitmapfwd(1,221,261) = fitmapfwd(61) !CaSO4
    catanfitmapfwd(1,204,261) = fitmapfwd(62) !NH42SO4
    nuecat(204,261) = 2.0_wp
    !
    catanfitmapfwd(1,223,261) = fitmapfwd(63) !MgSO4
    catanfitmapfwd(1,205,248) = fitmapfwd(64) !H2SO4 (H+ and HSO4-)
    catanfitmapfwd(1,201,248) = fitmapfwd(65) !LiHSO4
    catanfitmapfwd(1,203,248) = fitmapfwd(66) !KHSO4
    catanfitmapfwd(1,202,248) = fitmapfwd(67) !NaHSO4
    catanfitmapfwd(1,221,248) = fitmapfwd(68) !Ca(HSO4)2
    nuean(221,248) = 2.0_wp
    catanfitmapfwd(1,204,248) = fitmapfwd(69) !NH4HSO4
    catanfitmapfwd(1,223,248) = fitmapfwd(70) !Mg(HSO4)2
    nuean(223,248) = 2.0_wp
    !
    catanfitmapfwd(1,205,244) = fitmapfwd(71) !HI
    catanfitmapfwd(1,201,244) = fitmapfwd(72) !LiI
    catanfitmapfwd(1,203,244) = fitmapfwd(73) !KI
    catanfitmapfwd(1,202,244) = fitmapfwd(74) !NaI
    catanfitmapfwd(1,221,244) = fitmapfwd(75) !CaI2
    nuean(221,244) = 2.0_wp
    catanfitmapfwd(1,204,244) = fitmapfwd(76) !NH4I
    catanfitmapfwd(1,223,244) = fitmapfwd(77) !MgI2
    nuean(223,244) = 2.0_wp
    !
    catanfitmapfwd(1,205,262) = fitmapfwd(78) !H2CO3 (H+ and CO3--)
    nuecat(205,262) = 2.0_wp
    catanfitmapfwd(1,201,262) = fitmapfwd(79) !Li2CO3
    nuecat(201,262) = 2.0_wp
    catanfitmapfwd(1,203,262) = fitmapfwd(80) !K2CO3
    nuecat(203,262) = 2.0_wp
    catanfitmapfwd(1,202,262) = fitmapfwd(81) !Na2CO3
    nuecat(202,262) = 2.0_wp
    catanfitmapfwd(1,221,262) = fitmapfwd(82) !CaCO3
    catanfitmapfwd(1,204,262) = fitmapfwd(83) !NH42CO3
    nuecat(204,262) = 2.0_wp
    catanfitmapfwd(1,223,262) = fitmapfwd(84) !MgCO3
    !
    catanfitmapfwd(1,205,250) = fitmapfwd(85) !H2CO3 (H+ and HCO3-)
    catanfitmapfwd(1,201,250) = fitmapfwd(86) !LiHCO3
    catanfitmapfwd(1,203,250) = fitmapfwd(87) !KHCO3
    catanfitmapfwd(1,202,250) = fitmapfwd(88) !NaHCO3
    catanfitmapfwd(1,221,250) = fitmapfwd(89) !Ca(HCO3)2
    nuean(221,250) = 2.0_wp
    catanfitmapfwd(1,204,250) = fitmapfwd(90) !NH4HCO3
    catanfitmapfwd(1,223,250) = fitmapfwd(91) !Mg(HCO3)2
    nuean(223,250) = 2.0_wp
    !
    catanfitmapfwd(1,205,247) = fitmapfwd(92) !H2O (H+ and OH-)
    catanfitmapfwd(1,201,247) = fitmapfwd(93) !LiOH
    catanfitmapfwd(1,203,247) = fitmapfwd(94) !KOH
    catanfitmapfwd(1,202,247) = fitmapfwd(95) !NaOH
    catanfitmapfwd(1,221,247) = fitmapfwd(96) !Ca(OH)2
    nuean(221,247) = 2.0_wp
    catanfitmapfwd(1,204,247) = fitmapfwd(97) !NH4OH
    catanfitmapfwd(1,223,247) = fitmapfwd(98) !Mg(OH)2
    nuean(223,247) = 2.0_wp
    !
    catanfitmapfwd(1,205,246) = fitmapfwd(99) !HIO3 (H+ and IO3-)
    catanfitmapfwd(1,201,246) = fitmapfwd(100) !LiIO3
    catanfitmapfwd(1,203,246) = fitmapfwd(101) !KIO3
    catanfitmapfwd(1,202,246) = fitmapfwd(102) !NaIO3
    catanfitmapfwd(1,221,246) = fitmapfwd(103) !Ca(IO3)2
    nuean(221,246) = 2.0_wp
    catanfitmapfwd(1,204,246) = fitmapfwd(104) !NH4IO3
    catanfitmapfwd(1,223,246) = fitmapfwd(105) !Mg(IO3)2
    nuean(223,246) = 2.0_wp

    end subroutine MapIons2Fitpars3
    !-------------------------------------------------------------------------------------
     

    
!Interaction param numbers (= 27 + interaction index)
!       H, Li, K, Na, Ca, NH4, Mg
!Cl     1   2   3   4   5   6   7       
!Br     8   9   10  11  12  13  14
!NO3    15  16  17  18  19  20  21
!SO4    22  23  24  25  26  27  28
!HSO4   29  30  31  32  33  34  35
!I      36  37  38  39  40  41  42
    
!       H, Li, K, Na, Ca, NH4, Mg
!Cl     28  29  30  31  32  33  34       
!Br     35  36  37  38  39  40  41
!NO3    42  43  44  45  46  47  48
!SO4    49  50  51  52  53  54  55
!HSO4   56  57  58  59  60  61  62
!I      63  64  65  66  67  68  69
    
    
!       H,      Li,     K,          Na,     Ca,     NH4,    Mg
!Cl     28:29  30:31    32:33    34:35    36:37    38:39    40:41       
!Br     42:43  44:45    46:47    48:49    50:51    52:53    54:55
!NO3    56:57  58:59    60:61    62:63    64:65    66:67    68:69
!SO4    70:71  72:73    74:75    76:77    78:79    80:81    82:83
!HSO4   84:85  86:87    88:89    90:91    92:93    94:95    96:97
!I      98:99  100:101  102:103  104:105  106:107  108:109  110:111
!=========================================================    

end module ModViscEyring