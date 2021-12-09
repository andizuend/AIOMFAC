!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Submodule of module ModSystemProp containing subroutines to set system property    *
!*   variables, such as the list of components and subgroups present in the current     *
!*   mixture system.                                                                    *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2005                                                            *
!*   -> latest changes: 2021-10-03                                                      *
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
!*   -  SUBROUTINE definemixtures                                                       *
!*   -  SUBROUTINE defElectrolytes                                                      *
!*   -  SUBROUTINE SetMolarMass                                                         *
!*                                                                                      *
!****************************************************************************************
SUBMODULE (ModSystemProp) SubModDefSystem

IMPLICIT NONE
    
!==========================================================================================================================
    CONTAINS
!==========================================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Main subroutine to set and initialize system properties via the ITAB list of       *
    !*   'dataset_components' and an initial call of 'definemixtures'.                      *
    !*   When input data from a file (external source) is used, optional arguments must be. *
    !*   specified (such as the number of components and their subgroups).                  *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2018-05-24                                                      *
    !*   -> latest changes: 2021-09-30                                                      *
    !*                                                                                      *
    !****************************************************************************************
    MODULE SUBROUTINE SetSystem(ndi, datafromfile, ninp, cpnameinp, cpsubginp)

    USE ModSystemProp, ONLY : ninput, topsubno, bisulfsyst, definemixtures, frominpfile, cpname, &
        & bicarbsyst, noCO2input

    IMPLICIT NONE
    !interface variables:
    INTEGER(4),INTENT(IN) :: ndi            !dataset no. identifier (when not using dataset input from a file).
    LOGICAL(4),INTENT(IN) :: datafromfile   !set .true. if the system components are provided from an input file
    INTEGER(4),INTENT(IN) :: ninp           !value known at input only in the case of input from a file
    !optional input arguments:
    CHARACTER(LEN=200),DIMENSION(ninp),INTENT(IN),  OPTIONAL :: cpnameinp
    INTEGER(4),DIMENSION(ninp,topsubno),INTENT(IN), OPTIONAL :: cpsubginp
    !...
    !local variables
    INTEGER(4) :: AnFree, CatFree, i, ninputcomp, nnp1
    INTEGER(4),DIMENSION(:),ALLOCATABLE :: compID, compIDdat
    INTEGER(4),DIMENSION(:,:),ALLOCATABLE :: cpsubg, cpsubgdat
    LOGICAL(4) :: Hexists, HSO4exists, SO4exists, updbisulf
    LOGICAL(4) :: HCO3exists, CO3exists, updbicarb
    !........................................................

    frominpfile = datafromfile
    !(1b) allocate temporary arrays to contain the component information of the system:
    ninputcomp = ninp
    ALLOCATE( compID(ninputcomp+6), cpsubg(ninputcomp+6,topsubno), cpname(ninputcomp+4) ) 
    compID = 1500 !in this input file case
    cpsubg = 0
    cpname = "!?!"
    !(2b) transfer data from input via interface:
    DO i = 1,topsubno
        cpsubg(1:ninputcomp,i) = cpsubginp(1:ninputcomp, i)
    ENDDO
    cpname(1:ninputcomp) = cpnameinp(1:ninputcomp)
    ninput = ninputcomp !save original number of cpsubg-defined input components; i.e. prior to potential automatic changes to ITAB.

    !(3) check whether H+, HSO4- and SO4-- are part of the system or could be forming, 
    !    which may then require an adjustment to cpsubg.
    !    Check also whether H+, HCO3- and CO3-- are part of the system.
    Hexists = .false.
    HSO4exists = .false. 
    SO4exists = .false. 
    bisulfsyst = .false.
    updbisulf = .false.
    HCO3exists = .false. 
    CO3exists = .false. 
    bicarbsyst = .false.
    noCO2input = .false. !check if CO2 is an input or not
    updbicarb = .false.
    IF ( ANY(cpsubg(1:ninputcomp,205) > 0) ) THEN
        Hexists = .true.
    ENDIF
    IF ( ANY(cpsubg(1:ninputcomp,248) > 0) ) THEN
        HSO4exists = .true.
        bisulfsyst = .true.
    ENDIF
    IF ( ANY(cpsubg(1:ninputcomp,261) > 0) ) THEN
        SO4exists = .true.
    ENDIF
    IF ( ANY(cpsubg(1:ninputcomp,250) > 0) ) THEN
        HCO3exists = .true.
        bicarbsyst = .true.
    ENDIF
    IF ( ANY(cpsubg(1:ninputcomp,262) > 0) ) THEN
        CO3exists = .true.
    ENDIF
    IF (HSO4exists) THEN
        IF ( (.NOT. Hexists) .OR. (.NOT. SO4exists) ) THEN
            updbisulf = .true.  !HSO4- present but H+ and/or SO4-- not yet present at input
        ENDIF
    ELSE IF ( Hexists .AND. SO4exists ) THEN
        updbisulf = .true.      !H+ and SO4-- present but HSO4- not yet present at input
    ELSE IF ( HCO3exists .AND. SO4exists ) THEN
        updbisulf = .true.      !HCO3- and SO4-- present but HSO4- not yet present at input
    ENDIF
    IF (updbisulf) THEN
        bisulfsyst = .true.
        CatFree = 0
        AnFree = 0
        !determine the component number of the first electrolyte component:
        DO i = 1,ninputcomp 
            IF ( ANY(cpsubg(i,201:240) > 0) ) THEN !loop over all cations (as one must be part of first electrolyte component);
                nnp1 = i
                EXIT !leave loop
            ENDIF
        ENDDO
        IF (.NOT. Hexists) THEN !there was no H+ before the dissociation:
            !Find first "new" component in cpsubg that could contain the new cation:
            IF (ALL(cpsubg(ninputcomp+1,201:240) == 0)) THEN !found first cation-free cpsubg component
                CatFree = ninputcomp+1  !possibly anion free component...
                cpsubg(CatFree,205) = 1  !1x H+ for now, potentially corrected below if paired with SO4--
                ninputcomp = ninputcomp+1
            ENDIF
        ENDIF
        IF (.NOT. SO4exists) THEN !there was no SO4-- before the dissociation = > create a new anion number group:
            DO i = nnp1,ninputcomp+1
                !Find first "new" component in cpsubg that could contain the new anion:
                IF (ALL(cpsubg(i,241:topsubno) == 0)) THEN !found first anion-free cpsubg component
                    AnFree = i  ! anion-free component...
                    cpsubg(AnFree,261) = 1  !1x SO4--
                    ninputcomp = MAX(ninputcomp, i)
                    EXIT !leave the do-loop
                ENDIF
            ENDDO
        ENDIF
        IF (.NOT. HSO4exists) THEN !there was no HSO4- before the dissociation:
            DO i = nnp1,ninputcomp+1
                !Find first "new" component in cpsubg that could contain the new anion:
                IF (ALL(cpsubg(i,241:topsubno) == 0)) THEN !found first anion-free cpsubg component
                    AnFree = i  !possibly anion free component...
                    cpsubg(AnFree,248) = 1  !1x HSO4-
                    ninputcomp = MAX(ninputcomp, i)
                    EXIT !leave the do-loop
                ENDIF
            ENDDO
        ENDIF
        !check and guarantee a neutral ion pairing of a new elecrolyte component in cpsubg:
        IF (AnFree == CatFree) THEN
            !check stoichimetric number of H+ cations in electrolyte
            IF (cpsubg(CatFree,205) == 1) THEN
                IF (cpsubg(CatFree,261) > 0) THEN
                    cpsubg(CatFree,205) = 2 !it needs two H+ for one SO4--
                ENDIF
            ENDIF
        ELSE IF (AnFree > CatFree) THEN
            IF (cpsubg(AnFree,205) == 0) THEN
                IF (cpsubg(AnFree,248) > 0) THEN !248 is corresp. anion
                    cpsubg(AnFree,205) = 1
                ELSE !261 is corresponding anion
                    cpsubg(AnFree,205) = 2 !since it needs two H+ for one SO4--
                ENDIF
            ENDIF
        ELSE IF (CatFree > AnFree) THEN !need to put an anion to match one H+ (choose HSO4-)
            cpsubg(CatFree,248) = 1
        ENDIF
        !----------
    ENDIF !bisulfsyst
    
    !(4) H+, HCO3- and CO3-- are part of the system or could be forming, 
    !    which may then require an adjustment to cpsubg.
    IF (HCO3exists) THEN
        IF ( (.NOT. Hexists) .OR. (.NOT. CO3exists) ) THEN
            updbicarb = .true. !HCO3- present but H+ and/or CO3-- not yet present at input
        ENDIF
    ELSE IF ( Hexists .AND. CO3exists ) THEN
        updbicarb = .true. !H+ and CO3-- present but HCO3- not yet present at input
    ELSE IF ( HSO4exists .AND. CO3exists ) THEN
        updbicarb = .true. !HSO4- and CO3-- present but HCO3- not yet present at input
    ENDIF
    
    !determine the component number of the first electrolyte component:
    IF (bicarbsyst .OR. updbicarb) THEN
        DO i = 1,ninputcomp 
            IF ( ANY(cpsubg(i,201:240) > 0) ) THEN !loop over all cations (as one must be part of first electrolyte component);
                nnp1 = i
                EXIT !leave loop
            ENDIF
        ENDDO
    ENDIF
    IF (updbicarb) THEN
        bicarbsyst = .true.
        CatFree = 0
        AnFree = 0
        IF (.NOT. Hexists) THEN !there was no H+ before the dissociation:
            !Find first "new" component in cpsubg that could contain the new cation:
            IF (ALL(cpsubg(ninputcomp+1,201:240) == 0)) THEN !found first cation-free cpsubg component
                CatFree = ninputcomp+1  !possibly anion free component...
                cpsubg(CatFree,205) = 1  !1x H+ for now, potentially corrected below if paired with CO3--
                ninputcomp = ninputcomp +1
            ENDIF
        ENDIF
        IF (.NOT. CO3exists) THEN !there was no CO3-- before the dissociation = > create a new anion number group:
            DO i = nnp1,ninputcomp+1
                !Find first "new" component in cpsubg that could contain the new anion:
                IF (ALL(cpsubg(i,241:topsubno) == 0)) THEN !found first anion-free cpsubg component
                    AnFree = i  ! anion-free component...
                    cpsubg(AnFree,262) = 1  !1x CO3--
                    ninputcomp = MAX(ninputcomp, i)
                    EXIT !leave the do-loop
                ENDIF
            ENDDO
        ENDIF
        IF (.NOT. HCO3exists) THEN !there was no HCO3- before the dissociation:
            DO i = nnp1,ninputcomp+1
                !Find first "new" component in cpsubg that could contain the new anion:
                IF (ALL(cpsubg(i,241:topsubno) == 0)) THEN !found first anion-free cpsubg component
                    AnFree = i  !possibly anion free component...
                    cpsubg(AnFree,250) = 1  !1x HCO3-
                    ninputcomp = MAX(ninputcomp, i)
                    EXIT !leave the do-loop
                ENDIF
            ENDDO
        ENDIF
        !check and guarantee a neutral ion pairing of a new elecrolyte component in cpsubg:
        IF (AnFree == CatFree) THEN
            !check stoichimetric number of H+ cations in electrolyte
            IF (cpsubg(CatFree,205) == 1) THEN
                IF (cpsubg(CatFree,262) > 0) THEN
                    cpsubg(CatFree,205) = 2 !it needs two H+ for one CO3--
                ENDIF
            ENDIF
        ELSE IF (AnFree > CatFree) THEN
            IF (cpsubg(AnFree,205) == 0) THEN
                IF (cpsubg(AnFree,250) > 0) THEN !250 is corresp. anion
                    cpsubg(AnFree,205) = 1
                ELSE !262 is corresponding anion
                    cpsubg(AnFree,205) = 2 !since it needs two H+ for one CO3--
                ENDIF
            ENDIF
        ELSE IF (CatFree > AnFree) THEN !need to put an anion to match one H+ (choose HCO3-)
            cpsubg(CatFree,250) = 1
        ENDIF
        !----------
    ENDIF !updbicarb
    
    !Add CO2(aq) as an additional component if bicarbsys = .true.    
    IF (bicarbsyst) THEN
        IF (.NOT. ANY(cpsubg(1:ninputcomp,247) > 0)) THEN !Check for OH- existence
            ninputcomp = ninputcomp + 1
            cpsubg(ninputcomp,247) = 1  !1x OH-
            cpsubg(ninputcomp,205) = 1  !1x H+
        ENDIF
        IF (.NOT. ANY(cpsubg(1:nnp1,173) > 0)) THEN   
            noCO2input = .true.
            ALLOCATE(cpsubgdat(nnp1+1:ninputcomp+1,201:topsubno))
            cpsubgdat(nnp1+1:ninputcomp+1,201:topsubno) = cpsubg(nnp1:ninputcomp,201:topsubno)   !make a temporary array for the electrolytes
            cpsubg(nnp1:ninputcomp,201:topsubno) = 0
            cpsubg(nnp1,173) = 1        !add CO2(aq.) as the last neutral component
            compID(nnp1) = 402
            cpname(nnp1) = 'CO2(aq)'
            idCO2 = nnp1
            !transfer back the electrolyte component input information:
            cpsubg(nnp1+1:ninputcomp+1,201:topsubno) = cpsubgdat(nnp1+1:ninputcomp+1,201:topsubno)
            ninputcomp = ninputcomp +1
            DEALLOCATE(cpsubgdat)
        ELSE
            idCO2 = FINDLOC(cpsubg(1:nnp1,173), 1, DIM=1)
        ENDIF
    ELSE
        idCO2 = 0
    ENDIF
        
    !re-package data for transfer, since the number of input component and array sizes may have changed:
    ALLOCATE( compIDdat(ninputcomp), cpsubgdat(ninputcomp, topsubno) )
    compIDdat = 0
    cpsubgdat = 0
    compIDdat(1:ninputcomp) = compID(1:ninputcomp)
    DO i = 1,topsubno
        cpsubgdat(1:ninputcomp, i) = cpsubg(1:ninputcomp, i)
    ENDDO

    !(final step) define mixture properties and component properties of this system, including mapping of cpsubg to ITAB:
    CALL definemixtures(ndi, ninputcomp, compIDdat, cpsubgdat) 

    DEALLOCATE(compID, cpsubg, cpname, cpsubgdat, compIDdat)
    
    END SUBROUTINE SetSystem 
!==========================================================================================================================
    

    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   This subroutine defines the species (composition) of a system (nd) and calls       *
    !*   several subroutines to define SR, LR, and MR data that do *not* change with a      *
    !*   change in composition of the same mixture (e.g. molar masses, # ions, etc.).       *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich,                                                                *
    !*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA                *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2005                                                            *
    !*   -> latest changes: 2018/05/30                                                      *
    !*                                                                                      *
    !****************************************************************************************
    MODULE SUBROUTINE definemixtures(ndi, ninputcomp, compID, cpsubg)

    USE ModSubgroupProp, ONLY : cpsubgrstring, GroupMW, ioncharge, NKTAB
    USE ModAIOMFACvar, ONLY : AllocModAIOMFACvar
    USE ModMRpart, ONLY : MRinteractcoeff
    USE ModSRunifac, ONLY : SRsystm
    USE ModComponentNames, ONLY : names_mix

    IMPLICIT NONE
    !Interface variables:
    INTEGER(4),INTENT(IN) :: ndi, ninputcomp    
    INTEGER(4),DIMENSION(:),INTENT(IN) :: compID           
    INTEGER(4),DIMENSION(:,:),INTENT(IN) :: cpsubg                  !DIMENSION(ninputcomp,topsubno), list of component subgroups and corresponding quantity (e.g. from input file or from dataset_components)
    !Local variables:
    INTEGER(4) :: i, II, J, JJ, K, Icheck, countall, q, qq, ID
    INTEGER(4),DIMENSION(200) :: SolvSubs2                          !temporary array for solvent subgroups (we allow a maximum of 200)
    INTEGER(4),DIMENSION(ninputcomp*2) :: ElectSubs2
    INTEGER(4),DIMENSION(201:topsubno) :: ElectPos
    LOGICAL(4) :: already1, already2
    !--------------------------------------

    !initialize arrays and parameters:
    nd = 1          !for web-version (single data file / mixture only)
    solvmixrefnd = .false.
    errorflagmix = 0
    IF (ALLOCATED(ITAB)) THEN
        DEALLOCATE(ITAB, ITABsr, ITABMG, ITAB_dimflip, compN, Imaingroup, maingrindexofsubgr)
    ENDIF
    ALLOCATE( ITAB(ninputcomp, topsubno), ITAB_dimflip(topsubno, ninputcomp), compN(ninputcomp) )
    !assign the different mixture components from interface input data to ITAB:
    ITAB = cpsubg
    ITAB_dimflip = TRANSPOSE(ITAB)
    CompN = compID

    !count the number of neutral components, nneutral:
    nneutral = 0
    Icheck = 0
    DO i = 1,ninputcomp
        DO J = 1,200
            IF (ITAB_dimflip(J,i) /= 0) THEN
                IF (Icheck /= i) THEN
                    nneutral = nneutral+1
                    IF (i > nneutral) THEN
                        WRITE(*,*) "ERROR: there are skipped neutral component numbers! Please check your input."
                        READ(*,*)
                        STOP
                    ENDIF
                ENDIF
                Icheck = i
            ENDIF
        ENDDO
    ENDDO

    !count the number of electrolyte (input) components, nelectrol.
    !nelectrol will be updated below to account for additional potential electrolytes that may form from different cation-anion combinations!
    nelectrol = 0
    Icheck = 0
    DO i = nneutral+1,ninputcomp
        DO J = 201,topsubno
            IF (ITAB_dimflip(J,i) /= 0) THEN
                IF (Icheck /= i) THEN
                    nelectrol = nelectrol+1 !count the same ion only once when existing in different components
                    IF (i-nneutral > nelectrol) THEN
                        WRITE(*,*) "ERROR: there are skipped salt component numbers! Please check your input."
                        READ(*,*)
                        STOP
                        WRITE(*,*) ""
                    ENDIF
                ENDIF
                Icheck = i
            ENDIF
        ENDDO
    ENDDO

    !count the number of subgroups, NG, and other group numbers and properties:  
    Icheck = -1
    NG = 0
    NGN = 0
    NGS = 0
    SolvSubs2 = 0
    ElectSubs2 = 0
    countall = 0
    K = 1
    DO J = 1,200
        DO i = 1,ninputcomp
            IF (ITAB(i,J) /= 0) THEN
                countall = countall+1
                IF (Icheck /= J) THEN
                    SolvSubs2(K) = J !SolvSubs are the different solvent subgroups present in the mixture
                    K = K+1
                    NGN = NGN+1 !count the same subgroups only once, even if they are defined in different components
                ENDIF
                Icheck = J
            ENDIF
        ENDDO
    ENDDO
    IF (ALLOCATED(SolvSubs)) THEN
        DEALLOCATE(SolvSubs)
    ENDIF
    ALLOCATE( SolvSubs(NGN) )
    SolvSubs(1:NGN) = SolvSubs2(1:NGN)
    K = 1
    DO J = 201,topsubno
        DO i = 1,ninputcomp
            IF (ITAB(i,J) /= 0) THEN
                countall = countall+1
                IF (Icheck /= J) THEN
                    ElectSubs2(K) = J !the electrolyte subgroups (ions) present in the mixture
                    K = K+1
                    NGS = NGS+1 !count the same ions of salts / electrolytes only once, even if they exist in different components
                ENDIF
                Icheck = J
            ENDIF
        ENDDO
    ENDDO

    !check and set switch for presence/absence of water as a component in the system:
    IF (ANY(SolvSubs(1:NGN) == 16)) THEN !water is present
        waterpresent = .true.
        DO i = 1,ninputcomp
            IF (ITAB(i,16) == 1) THEN
                CompN(i) = 401 !use for cases where input files is used to assign components
                EXIT
            ENDIF
        ENDDO
    ELSE
        waterpresent = .false.
    ENDIF

    !now bring the subgroups into the right positions according to the components in ITAB (so that anions and cations of neutral electrolyte formula units are coupled):
    !the first anion in ElectSubs should correspond to the first cation in array ElectSubs.
    q = 1
    ElectPos = 0
    DO i = nneutral+1,ninputcomp
        DO J = 1,NGS
            K = ElectSubs2(J)
            IF (ITAB_dimflip(K,i) > 0 .AND. ElectPos(K) == 0) THEN  !to make sure that every elect. subgroup is only checked once.
                ElectPos(K) = q !the position in the later array of ElectSubs, so that the anions and cations in ElectSubs correspond.
                q = q+1
            ENDIF
        ENDDO
    ENDDO
    IF (ALLOCATED(ElectSubs)) THEN
        DEALLOCATE(ElectSubs, AllSubs)
    ENDIF
    ALLOCATE( ElectSubs(NGS), AllSubs(NGN+NGS) )
    ElectSubs = 0
    AllSubs = 0
    DO J = 1,NGS
        K = ElectSubs2(J)
        ElectSubs(ElectPos(K)) = K
    ENDDO

    !Define array AllSubs which contains all (different) subgroups present in the mixture:
    isPEGsystem = .false.
    DO K = 1,NGN
        AllSubs(K) = SolvSubs(K)
        IF (SolvSubs(K) == 154) THEN
            isPEGsystem = .true.
        ENDIF
    ENDDO
    DO K = 1,NGS
        AllSubs(NGN+K) = ElectSubs(K)
    ENDDO

    !------------------------------------------------------------------------------------------------------
    !Calculating Ication and Ianion: saves the different cations/anions in the arrays Ication/Ianion at the 
    !cation/anion Nr. array index position (in ascending order of the type nr.); Also assign cation and anion charges.
    NGI = MAX( COUNT(ElectSubs(1:NGS) < 241), COUNT(ElectSubs(1:NGS) > 240) )
    IF (ALLOCATED(Ication)) THEN
        DEALLOCATE(Ication, Ianion, cationZ, anionZ)
    ENDIF
    ALLOCATE( Ication(NGI), Ianion(NGI), cationZ(NGI), anionZ(NGI) )
    ICation = 0
    IAnion = 0
    CatNr = 0
    AnNr = 0
    cationZ = 0.0D0
    anionZ = 0.0D0
    DO J = 1,NGS
        already1 = .false.
        already2 = .false.
        K = ElectSubs(J)
        IF (K < 241 .AND. K > 200) THEN
            DO qq = 1,NGS
                IF (.NOT. already1) THEN
                    IF (Ication(qq) == 0 .OR. Ication(qq) == K) THEN
                        Ication(qq) = K
                        CatNr(K) = qq
                        cationZ(qq) = REAL(Ioncharge(K), kind=8)
                        already1 = .true.
                    ENDIF
                ENDIF
            ENDDO
        ELSE IF (K <= topsubno .AND. K > 240) THEN !K > 240
            DO qq = 1,NGS
                IF (.NOT. already2) THEN
                    IF (Ianion(qq) == 0 .OR. Ianion(qq) == K) THEN
                        Ianion(qq) = K
                        AnNr(K) = qq
                        anionZ(qq) = REAL(Ioncharge(K), kind=8)
                        already2 = .true.
                    ENDIF
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    
    !save index locations of sulfuric acid / bisulfate ions (of use later if present in the mixture):
    idH = CatNr(205)
    idHSO4 = AnNr(248)
    idSO4 = AnNr(261)
    idCa = CatNr(221)
    
    idHCO3 = AnNr(250)
    idCO3 = AnNr(262)
    idOH = AnNr(247)

    !------------------------------------------------------------
    !Calculation of the number of different cations and the number of different anions:
    Ncation = COUNT(Ication(:) > 0)
    Nanion = COUNT(Ianion(:) > 0)

    !define different binary cation-anion combinations as electrolyte components of the mixture:
    IF (nelectrol > 0) THEN
        elpresent = .true.
        !initialize KVLE(cation ID, anion ID) parameters for gas-liquid equilibria of different electrolyte components composed of binary cation+anion pairs
        !FROM LANGE's Handbook of Chemistry 15th Edition, Gibbs Energy of Formation Data; with R* = 8.314459 [J K^-1 mol^-1]
        KVLE_298K = 0.0D0
        KVLE_298K(205,243) = 7.236814145D8  !HBr  7.235752425D8
        KVLE_298K(205,242) = 1.986870884D6  !HCl
        KVLE_298K(205,241) = 5.749865093D8  !HF
        KVLE_298K(205,244) = 2.168075448D9  !HI
        KVLE_298K(205,245) = 4.191055023D6  !HNO3
        !KVLE_298K(205,204) = 1.07614415D11 !NH3  (i.e. NH3(g) + H+(aq) <--> NH4+(aq)
        !initialize the IAP vs. aw parameterization coefficients for semi-volatile electrolyte components considered:
        IAPcoeffs = 0.0D0
        IAPcoeffs(205,243,1:3) = [4.341292D-2, 1.119213D0, 5.644364D-1]  !HBr
        IAPcoeffs(205,242,1:3) = [4.391087D-2, 1.196585D0, 4.603556D-1]  !HCl
        IAPcoeffs(205,245,1:3) = [4.8857D-2, 1.3369D0, 3.423756D-1]      !HNO3
        !--
        CALL defElectrolytes(nneutral, NGS, nelectrol)
        !--
    ELSE
        elpresent = .false.
        IF (ALLOCATED(nuestoich)) DEALLOCATE(nuestoich)
        ALLOCATE( nuestoich(nneutral) )
        nuestoich = 1.0D0
    ENDIF

    NGI = MAX(Nanion, Ncation)  !max number of ion groups of either cation or anion kind for present input composition.
    NG = NGN + NGS !maximum number of molecular / ionic groups
    NKNpNGS = nneutral + NGS

    !------------------------------------------------------------
    !Calculation of ITABMG. ITABMG is the analogue to ITAB for the main groups of the neutral components.
    ALLOCATE( ITABMG(nneutral,Nmaingroups), Imaingroup(NGN), maingrindexofsubgr(NGN) )
    ITABMG = 0  !initialize
    Imaingroup = 0
    maingrindexofsubgr = 0
    DO i = 1,nneutral
        DO J = 1,NGN
            II = SolvSubs(J)
            JJ = NKTAB(II)
            ITABMG(i,JJ) = ITABMG(i,JJ) + ITAB_dimflip(II,i)
            Imaingroup(J) = JJ  !list of main groups (here not yet filtered and sorted)
        ENDDO
        !check input correctness with respect to OH groups and special CHn groups bonded to OH groups:
        IF (ITABMG(i,68) > 0) THEN
            IF (ITABMG(i,68) > ITABMG(i,69)) THEN !this would indicate an input error
                errorflagmix = 13
            ENDIF
        ENDIF
    ENDDO

    !------------------------------------------------------------
    !Calculation of the neutral subgroups molar masses; GroupMW is a list of the molecular masses of all neutral subgroups in the mixture.
    IF (ALLOCATED(SubGroupMW)) THEN
        DEALLOCATE(SubGroupMW)
    ENDIF
    ALLOCATE( SubGroupMW(NGN) )
    SubGroupMW(1:NGN) = GroupMW(SolvSubs(1:NGN))*1.0D-3 !assign subgroup molar mass in [kg/mol]

    !Calculation of Imaingroup, the sorted list of neutral main groups:
    DO i = 1,NGN-1
        DO J = i+1,NGN
            IF (Imaingroup(i) == Imaingroup(J)) THEN !combine the same main groups
                Imaingroup(J) = 0
            ENDIF
        ENDDO
    ENDDO
    !filter now the empty (Imaingroup = 0) positions out:
    DO i = 1,NGN-1
        IF (Imaingroup(i) == 0) THEN
            DO J = i+1,NGN
                IF (Imaingroup(J) /= 0) THEN
                    Imaingroup(i) = Imaingroup(J)
                    Imaingroup(J) = 0
                    EXIT !exit the inner J-loop
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    !assemble a list stating for each subgroup of the present mixture to which main group of the present mixture it belongs (e.g. used in subroutine LR_MR_activity):
    DO j = 1,NGN !j-loop over subgroups
        JJ = NKTAB(SolvSubs(j)) !static main group number associated with subgroup
        !find the index number in current mixture of the main group with static number JJ:
        DO i = 1,NGN
            IF (Imaingroup(i) == JJ) THEN
                maingrindexofsubgr(j) = i
                EXIT !leave inner i-loop
            ENDIF
        ENDDO
    ENDDO
    !------------------------------------------------------------

    !Calculation of ITABsr: ITABsr is similar to ITAB, but the salts are divided into cations and anions as components themselves for use within UNIFAC...
    !Transfer information of cations (first) then anions into ITABsr list
    ALLOCATE( ITABsr(NKNpNGS,topsubno) )
    ITABsr = 0 !initialize
    K = nneutral
    JJ = nneutral+Ncation
    DO J = 1,NGS
        ID = ElectSubs(J)
        !!i = SUM(MAXLOC(ITAB(nneutral+1:nneutral+nelectrol, ID))) +nneutral
        IF (ID > 200 .AND. ID < 240) THEN !cation
            K = K+1
            ITABsr(K,ID) = 1   !transfer cation info; just set to one ion, actual amount is coming from molality of the ion
        ELSE IF (ID > 240) THEN !anion
            JJ = JJ+1
            ITABsr(JJ,ID) = 1  !transfer anion info; just set to one ion, actual ion amount is coming from molality of the ion
        ENDIF
    ENDDO
    !Transfer of information for neutrals from ITAB to ITABsr:
    DO i = 1,nneutral
        DO J = 1,NGN
            JJ = SolvSubs(J)
            ITABsr(i,JJ) = ITAB_dimflip(JJ,i)  !update ITABsr list.
        ENDDO
    ENDDO

    !save the number of independent components (ions combined to salts as given by ITAB):
    nindcomp = nneutral +nelectrol !number of components of the system (input components + additional electrolyte components from potential additional electroneutral ion combinations)

    IF (ALLOCATED(OtoCratio)) THEN
        DEALLOCATE( OtoCratio, HtoCratio, compname, compnameTeX, ionname, ionnameTeX )
    ENDIF
    ALLOCATE( OtoCratio(nindcomp), HtoCratio(nindcomp), compname(nindcomp), compnameTeX(nindcomp), ionname(NGS), ionnameTeX(NGS) )
    !set the name strings for the different independent components and the O:C ratio of the neutrals (when set already):
    CALL names_mix(CompN, compname, compnameTeX, ionname, ionnameTeX, OtoCratio, HtoCratio)     
    !now compname contains the name strings of all independent components of the present mixture nd;

    CALL cpsubgrstring() !generate for each component a subgroup-string stored in compsubgroups, compsubgroupsTeX, and compsubgroupsHTML

    !define the middle range (MR) coefficients of the actual mixture:
    CALL MRinteractcoeff()

    !calculate some short range (SR) parameters which don't change for the same mixture:
    CALL SRsystm(NKNpNGS)

    !allocate and populate array Mmass with the molar masses in component order.
    IF (ALLOCATED(Mmass)) THEN
        DEALLOCATE(Mmass)
    ENDIF
    ALLOCATE(Mmass(nindcomp))
    !calculate the molar masses of the different mixture species:
    CALL SetMolarMass(Mmass)    !Mmass lists the molar mass in [kg/mol] of all mixture components

    !-- allocate several composition-dependent variables from module ModAIOMFACvar:
    CALL AllocModAIOMFACvar()

    END SUBROUTINE definemixtures
    !==========================================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to define all possible electrolytes of a system based on input of       *
    !*   cations and anions (from input electrolyte components).                            *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2014/07/21                                                      *
    !*   -> latest changes: 2018/05/24                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    MODULE SUBROUTINE defElectrolytes(nneutral, NGS, nelectrol)

    USE ModSubgroupProp, ONLY : Ioncharge, IonO2Cequiv

    IMPLICIT NONE
    !Interface variables:
    INTEGER(4),INTENT(IN) :: nneutral, NGS
    INTEGER(4),INTENT(INOUT) :: nelectrol
    !Local variables:
    INTEGER(4) :: i, J, JJ, K, KK, q, zc, za, ion, iel, nelin, nnp1
    LOGICAL(4) :: already1
    REAL(8),PARAMETER :: onethird = 1.0D0/3.0D0
    !------------------------------------------------------------------

    nnp1 = nneutral + 1
    !calculate all different electrolyte components that can be formed from the input electro-neutral
    !cation-anion combinations (e.g. for gas-particle partitioning of different electrolytes):
    iel = Ncation*Nanion !maximum number of possible electrolyte components just formed by one cation and anion each.
    IF (ALLOCATED(ElectComps)) THEN
        DEALLOCATE( ElectComps, ElectNues, ElectVolatile, ElectO2Cequiv, K_el )
    ENDIF
    ALLOCATE( ElectComps(iel,2), ElectNues(iel,2), ElectVolatile(iel), ElectO2Cequiv(iel), K_el(nnp1:nneutral+iel) )
    IF (ALLOCATED(nuestoich)) DEALLOCATE(nuestoich)
    ALLOCATE( nuestoich(nneutral+iel) )

    !initialize arrays for electrolytes:
    ElectComps = 0
    ElectNues = 0
    ElectVolatile = .false.
    ElectO2Cequiv = 0.0D0
    K_el = 0.0D0
    nuestoich = 1.0D0

    !First, list the input electrolyte components in the order of definition in ITAB:
    nelin = nelectrol !ITAB-defined electrolyte components;  ninput - nneutral
    iel = 0
    DO i = 1,nelin
        KK = 0
        JJ = 0
        DO J = 1,NGS
            !loop over ElectSubs(K)
            ion = ElectSubs(J)
            IF (KK == 0 .AND. ion > 200 .AND. ion < 240) THEN
                IF (ITAB_dimflip(ion,nneutral+i) > 0) THEN
                    KK = ion
                ENDIF
            ELSE IF (JJ == 0 .AND. ion > 240 .AND. ion <= topsubno) THEN
                IF (ITAB_dimflip(ion,nneutral+i) > 0) THEN
                    JJ = ion
                ENDIF
            ENDIF
            IF (KK > 0 .AND. JJ > 0) THEN
                iel = iel+1
                ElectComps(iel,1) = KK
                ElectComps(iel,2) = JJ
                K_el(nneutral+iel) = KVLE_298K(KK,JJ) !assign gas-particle equilibrium constants of the electrolyte components based on the involved cations and anions and tabulated values in KVLE_298K
                EXIT !leave inner J-loop
            ENDIF
        ENDDO
        IF (KK == JJ) THEN !then they are both == 0
            CYCLE !continue with next i loop
        ENDIF
        zc = Ioncharge(KK)
        za = ABS(Ioncharge(JJ))
        !now attribute the stoichiometric cation and anion numbers to creat an electroneutral component:
        IF (zc == za) THEN !1:1 electrolyte (both the same charge of +-1 or +-2)
            ElectNues(iel,1) = 1
            ElectNues(iel,2) = 1
            ElectO2Cequiv(iel) = (IonO2Cequiv(KK)*IonO2Cequiv(JJ))**0.5D0
            nuestoich(nneutral+iel) = 2.0D0
        ELSE
            SELECT CASE(zc)
            CASE(1) !2:1 electrolyte
                ElectNues(iel,1) = 2
                ElectNues(iel,2) = 1 !anion charge: -2
                ElectO2Cequiv(iel) = (IonO2Cequiv(KK)**2 *IonO2Cequiv(JJ))**onethird
                nuestoich(nneutral+iel) = 3.0D0
            CASE(2) !1:2 electrolyte
                ElectNues(iel,1) = 1 !cation charge: +2
                ElectNues(iel,2) = 2
                ElectO2Cequiv(iel) = (IonO2Cequiv(KK)*IonO2Cequiv(JJ)**2)**onethird
                nuestoich(nneutral+iel) = 3.0D0
            END SELECT
        ENDIF
        IF (KK == 205 .AND. JJ /= 261 .AND. JJ /= 262 .AND. JJ /= 250 .AND. JJ /= 248 &
        & .OR. (KK == 204 .AND. JJ == 245)) THEN
            ElectVolatile(iel) = .true.
        ENDIF
    ENDDO !i
    !Second, make all other possible electroneutral binary cation-anion combinations for electrolytes and add to the lists:
    i = iel
    DO K = 1,Ncation
        KK = Ication(K)
        DO J = 1,Nanion
            JJ = Ianion(J)
            !check if cation-anion combination is already assigned to a component:
            already1 = .false.
            DO q = 1,i
                IF (ElectComps(q,1) == KK .AND. ElectComps(q,2) == JJ) THEN
                    already1 = .true.
                ENDIF
            ENDDO
            IF (.NOT. already1) THEN !found new electrolyte component
                i = i+1
                ElectComps(i,1) = KK
                ElectComps(i,2) = JJ
                K_el(nneutral+i) = KVLE_298K(KK,JJ)
                zc = Ioncharge(KK)
                za = ABS(Ioncharge(JJ))
                !now attribute the stoichiometric cation and anion numbers to creat an electroneutral component:
                IF (zc == za) THEN !1:1 electrolyte (both the same charge of +-1 or +-2)
                    ElectNues(i,1) = 1
                    ElectNues(i,2) = 1
                    ElectO2Cequiv(i) = (IonO2Cequiv(KK)*IonO2Cequiv(JJ))**0.5D0
                    nuestoich(nneutral+i) = 2.0D0
                ELSE
                    SELECT CASE(zc)
                    CASE(1) !2:1 electrolyte
                        ElectNues(i,1) = 2
                        ElectNues(i,2) = 1 !anion charge: -2
                        ElectO2Cequiv(i) = (IonO2Cequiv(KK)**2 *IonO2Cequiv(JJ))**onethird
                        nuestoich(nneutral+i) = 3.0D0
                    CASE(2) !1:2 electrolyte
                        ElectNues(i,1) = 1 !cation charge: +2
                        ElectNues(i,2) = 2
                        ElectO2Cequiv(i) = (IonO2Cequiv(KK)*IonO2Cequiv(JJ)**2)**onethird
                        nuestoich(nneutral+i) = 3.0D0
                    END SELECT
                ENDIF
                IF (KK == 205 .AND. JJ /= 261 .AND. JJ /= 262 .AND. JJ /= 250 .AND. JJ /= 248 &
                & .OR. (KK == 204 .AND. JJ == 245)) THEN
                    ElectVolatile(i) = .true.
                ENDIF
            ENDIF
        ENDDO !J
    ENDDO !K
    nelectrol = i !the total number of different possible electrolyte components of the system (including both input electrolytes and additional, electroneutral ion combinations)
    NG = NGN + NGS !maximum number of molecular / ionic groups

    END SUBROUTINE defElectrolytes
    !==========================================================================================================================


    !************************************************************************
    !*                                                                      *
    !*  Set the molar masses of the given mixture components;               *
    !*  saved in array MolarM --> Mmass (neutrals first, then electrolytes) *
    !*                                                                      *
    !************************************************************************
    PURE MODULE SUBROUTINE SetMolarMass(MolarM)

    USE ModSubgroupProp, ONLY : GroupMW, SMWC, SMWA

    IMPLICIT NONE

    REAL(8),DIMENSION(:),INTENT(OUT) :: MolarM  !Molar mass in [kg/mol]
    INTEGER(4) :: i, J, K
    !...............................................
    !SMWA     : list of the molar masses of the anions
    !SMWC     : list of the molar masses of the cations
    !GroupMW  : list of molar masses of the organic/water subgroups
    MolarM = 0.0D0  !initialize array
    !calculate the neutral compounds' molar mass from the values of their subgroups listed in GroupMW array:
    DO i = 1,NGN    !loop over all solvent subgroups
        J = SolvSubs(i)
        MolarM(1:nneutral) = MolarM(1:nneutral) + ITAB(1:nneutral,J)*GroupMW(J)*1.0D-3  !compound molar mass in [kg/mol]
    ENDDO
    !electrolytes:
    DO i = 1,nelectrol
        J = ElectComps(i,1)-200 !the cation of this electrolyte in the list of cations as numbered in SMWC (1:40 instead of 201:240)
        K = ElectComps(i,2)-240 !the anion "
        MolarM(nneutral+i) = ( ElectNues(i,1)*SMWC(J) + ElectNues(i,2)*SMWA(K) )*1.0D-3 ! Molar mass of electrolyte i in [kg/mol]
    ENDDO

    END SUBROUTINE SetMolarMass
!==========================================================================================================================
    
END SUBMODULE SubModDefSystem