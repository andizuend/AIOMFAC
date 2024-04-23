!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   submodule of module ModSystemProp containing subroutines to set system property    *
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
!*   -> latest changes: 2023-03-17                                                      *
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
!*   -  subroutine definemixtures                                                       *
!*   -  subroutine defElectrolytes                                                      *
!*   -  subroutine SetMolarMass                                                         *
!*                                                                                      *
!****************************************************************************************
submodule (ModSystemProp) SubModDefSystem

implicit none

integer,dimension(:),allocatable :: compID, compIDdat
integer,dimension(:,:),allocatable :: cpsubg, cpsubgdat

!$OMP THREADPRIVATE (compID, compIDdat, cpsubg, cpsubgdat)
    
!==========================================================================================================================
    contains
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
    !*   -> latest changes: 2023-03-17                                                      *
    !*                                                                                      *
    !****************************************************************************************
    module subroutine SetSystem(ndi, datafromfile, ninp, cpnameinp, cpsubginp)

    use ModSystemProp, only : ninput, topsubno, bisulfsyst, definemixtures, frominpfile, cpname, &
        & bicarbsyst, noCO2input

    implicit none
    !interface variables:
    integer,intent(in) :: ndi            !dataset no. identifier (when not using dataset input from a file).
    logical,intent(in) :: datafromfile   !set .true. if the system components are provided from an input file
    integer,intent(in) :: ninp           !value known at input only in the case of input from a file
    !optional input arguments:
    character(len=200),dimension(ninp),intent(in),  OPTIONAL :: cpnameinp
    integer,dimension(ninp,topsubno),intent(in), OPTIONAL :: cpsubginp
    !...
    !local variables
    integer :: AnFree, CatFree, i, ninputcomp, nnp1
    logical :: Hexists, HSO4exists, SO4exists, updbisulf
    logical :: HCO3exists, CO3exists, updbicarb
    !........................................................

    frominpfile = datafromfile
    if (.NOT. datafromfile) then 
        !define dataset via call to dataset_components
    else
        !(1b) allocate temporary arrays to contain the component information of the system:
        ninputcomp = ninp
        allocate( compID(ninputcomp+6), cpsubg(ninputcomp+6,topsubno), cpname(ninputcomp+4) ) 
        compID = 1500 !in this input file case
        cpsubg = 0
        cpname = "!?!"
        !(2b) transfer data from input via interface:
        do i = 1,topsubno
            cpsubg(1:ninputcomp,i) = cpsubginp(1:ninputcomp, i)
        enddo
        cpname(1:ninputcomp) = cpnameinp(1:ninputcomp)
    endif
    ninput = ninputcomp !save original number of cpsubg-defined input components; i.e. prior to potential automatic changes to ITAB.
    
    !(4) check whether H+, HSO4- and SO4-- are part of the system or could be forming, 
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
    if ( any(cpsubg(1:ninputcomp,205) > 0) ) then
        Hexists = .true.
    endif
    if ( any(cpsubg(1:ninputcomp,248) > 0) ) then
        HSO4exists = .true.
        bisulfsyst = .true.
    endif
    if ( any(cpsubg(1:ninputcomp,261) > 0) ) then
        SO4exists = .true.
    endif
    if ( any(cpsubg(1:ninputcomp,250) > 0) ) then
        HCO3exists = .true.
        bicarbsyst = .true.
    endif
    if ( any(cpsubg(1:ninputcomp,262) > 0) ) then
        CO3exists = .true.
    endif
    if (HSO4exists) then
        if ( (.NOT. Hexists) .OR. (.NOT. SO4exists) ) then
            updbisulf = .true.  !HSO4- present but H+ and/or SO4-- not yet present at input
        endif
    else if ( Hexists .AND. SO4exists ) then
        updbisulf = .true.      !H+ and SO4-- present but HSO4- not yet present at input
    else if ( HCO3exists .AND. SO4exists ) then
        updbisulf = .true.      !HCO3- and SO4-- present but HSO4- not yet present at input
    endif
    if (updbisulf) then
        bisulfsyst = .true.
        CatFree = 0
        AnFree = 0
        !determine the component number of the first electrolyte component:
        do i = 1,ninputcomp 
            if ( any(cpsubg(i,201:240) > 0) ) then !loop over all cations (as one must be part of first electrolyte component);
                nnp1 = i
                exit !leave loop
            endif
        enddo
        if (.NOT. Hexists) then !there was no H+ before the dissociation:
            !Find first "new" component in cpsubg that could contain the new cation:
            if (all(cpsubg(ninputcomp+1,201:240) == 0)) then !found first cation-free cpsubg component
                CatFree = ninputcomp+1  !possibly anion free component...
                cpsubg(CatFree,205) = 1  !1x H+ for now, potentially corrected below if paired with SO4--
                ninputcomp = ninputcomp+1
            endif
        endif
        if (.NOT. SO4exists) then !there was no SO4-- before the dissociation = > create a new anion number group:
            do i = nnp1,ninputcomp+1
                !Find first "new" component in cpsubg that could contain the new anion:
                if (all(cpsubg(i,241:topsubno) == 0)) then !found first anion-free cpsubg component
                    AnFree = i  ! anion-free component...
                    cpsubg(AnFree,261) = 1  !1x SO4--
                    ninputcomp = max(ninputcomp, i)
                    exit !leave the do-loop
                endif
            enddo
        endif
        if (.NOT. HSO4exists) then !there was no HSO4- before the dissociation:
            do i = nnp1,ninputcomp+1
                !Find first "new" component in cpsubg that could contain the new anion:
                if (all(cpsubg(i,241:topsubno) == 0)) then !found first anion-free cpsubg component
                    AnFree = i  !possibly anion free component...
                    cpsubg(AnFree,248) = 1  !1x HSO4-
                    ninputcomp = max(ninputcomp, i)
                    exit !leave the do-loop
                endif
            enddo
        endif
        !check and guarantee a neutral ion pairing of a new elecrolyte component in cpsubg:
        if (AnFree == CatFree) then
            !check stoichimetric number of H+ cations in electrolyte
            if (cpsubg(CatFree,205) == 1) then
                if (cpsubg(CatFree,261) > 0) then
                    cpsubg(CatFree,205) = 2 !it needs two H+ for one SO4--
                endif
            endif
        else if (AnFree > CatFree) then
            if (cpsubg(AnFree,205) == 0) then
                if (cpsubg(AnFree,248) > 0) then !248 is corresp. anion
                    cpsubg(AnFree,205) = 1
                else !261 is corresponding anion
                    cpsubg(AnFree,205) = 2 !since it needs two H+ for one SO4--
                endif
            endif
        else if (CatFree > AnFree) then !need to put an anion to match one H+ (choose HSO4-)
            cpsubg(CatFree,248) = 1
        endif
        !----------
    endif !bisulfsyst
    
    !(5) H+, HCO3- and CO3-- are part of the system or could be forming, 
    !    which may then require an adjustment to cpsubg.
    if (HCO3exists) then
        if ( (.NOT. Hexists) .OR. (.NOT. CO3exists) ) then
            updbicarb = .true. !HCO3- present but H+ and/or CO3-- not yet present at input
        endif
    else if ( Hexists .AND. CO3exists ) then
        updbicarb = .true. !H+ and CO3-- present but HCO3- not yet present at input
    else if ( HSO4exists .AND. CO3exists ) then
        updbicarb = .true. !HSO4- and CO3-- present but HCO3- not yet present at input
    endif
    
    !determine the component number of the first electrolyte component:
    if (bicarbsyst .OR. updbicarb) then
        do i = 1,ninputcomp 
            if ( any(cpsubg(i,201:240) > 0) ) then !loop over all cations (as one must be part of first electrolyte component);
                nnp1 = i
                exit !leave loop
            endif
        enddo
    endif
    if (updbicarb) then
        bicarbsyst = .true.
        CatFree = 0
        AnFree = 0
        if (.NOT. Hexists) then !there was no H+ before the dissociation:
            !Find first "new" component in cpsubg that could contain the new cation:
            if (all(cpsubg(ninputcomp+1,201:240) == 0)) then !found first cation-free cpsubg component
                CatFree = ninputcomp+1  !possibly anion free component...
                cpsubg(CatFree,205) = 1  !1x H+ for now, potentially corrected below if paired with CO3--
                ninputcomp = ninputcomp +1
            endif
        endif
        if (.NOT. CO3exists) then !there was no CO3-- before the dissociation = > create a new anion number group:
            do i = nnp1,ninputcomp+1
                !Find first "new" component in cpsubg that could contain the new anion:
                if (all(cpsubg(i,241:topsubno) == 0)) then !found first anion-free cpsubg component
                    AnFree = i  ! anion-free component...
                    cpsubg(AnFree,262) = 1  !1x CO3--
                    ninputcomp = max(ninputcomp, i)
                    exit !leave the do-loop
                endif
            enddo
        endif
        if (.NOT. HCO3exists) then !there was no HCO3- before the dissociation:
            do i = nnp1,ninputcomp+1
                !Find first "new" component in cpsubg that could contain the new anion:
                if (all(cpsubg(i,241:topsubno) == 0)) then !found first anion-free cpsubg component
                    AnFree = i  !possibly anion free component...
                    cpsubg(AnFree,250) = 1  !1x HCO3-
                    ninputcomp = max(ninputcomp, i)
                    exit !leave the do-loop
                endif
            enddo
        endif
        !check and guarantee a neutral ion pairing of a new electrolyte component in cpsubg:
        if (AnFree == CatFree) then
            !check stoichimetric number of H+ cations in electrolyte
            if (cpsubg(CatFree,205) == 1) then
                if (cpsubg(CatFree,262) > 0) then
                    cpsubg(CatFree,205) = 2 !it needs two H+ for one CO3--
                endif
            endif
        else if (AnFree > CatFree) then
            if (cpsubg(AnFree,205) == 0) then
                if (cpsubg(AnFree,250) > 0) then !250 is corresp. anion
                    cpsubg(AnFree,205) = 1
                else !262 is corresponding anion
                    cpsubg(AnFree,205) = 2 !since it needs two H+ for one CO3--
                endif
            endif
        else if (CatFree > AnFree) then !need to put an anion to match one H+ (choose HCO3-)
            cpsubg(CatFree,250) = 1
        endif
        !----------
    endif !updbicarb
    
    !Add CO2(aq) as an additional component if bicarbsys = .true.    
    if (bicarbsyst) then
        if (.NOT. any(cpsubg(1:ninputcomp,247) > 0)) then !Check for OH- existence
            ninputcomp = ninputcomp + 1
            cpsubg(ninputcomp,247) = 1  !1x OH-
            cpsubg(ninputcomp,205) = 1  !1x H+
        endif
        if (.NOT. any(cpsubg(1:nnp1,173) > 0)) then   
            noCO2input = .true.
            allocate(cpsubgdat(nnp1+1:ninputcomp+1,201:topsubno))
            cpsubgdat(nnp1+1:ninputcomp+1,201:topsubno) = cpsubg(nnp1:ninputcomp,201:topsubno)   !make a temporary array for the electrolytes
            cpsubg(nnp1:ninputcomp,201:topsubno) = 0
            cpsubg(nnp1,173) = 1        !add CO2(aq.) as the last neutral component
            compID(nnp1) = 402
            cpname(nnp1) = 'CO2(aq)'
            idCO2 = nnp1
            !transfer back the electrolyte component input information:
            cpsubg(nnp1+1:ninputcomp+1,201:topsubno) = cpsubgdat(nnp1+1:ninputcomp+1,201:topsubno)
            ninputcomp = ninputcomp +1
            deallocate(cpsubgdat)
        else
            idCO2 = findloc(cpsubg(1:nnp1,173), 1, DIM=1)
        endif
    else
        idCO2 = 0
    endif

    !re-package data for transfer, since the number of input component and array sizes may have changed:
    allocate( compIDdat(ninputcomp), cpsubgdat(ninputcomp, topsubno) )
    compIDdat = 0
    cpsubgdat = 0
    compIDdat(1:ninputcomp) = compID(1:ninputcomp)
    do i = 1,topsubno
        cpsubgdat(1:ninputcomp, i) = cpsubg(1:ninputcomp, i)
    enddo

    !(final step) define mixture properties and component properties of this system, including mapping of cpsubg to ITAB:
    call definemixtures(ndi, ninputcomp, compIDdat, cpsubgdat) 

    deallocate(compID, cpsubg, cpname, cpsubgdat, compIDdat)
    
    end subroutine SetSystem 
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
    !*   -> latest changes: 2018-05-30                                                      *
    !*                                                                                      *
    !****************************************************************************************
    module subroutine definemixtures(ndi, ninputcomp, compID, cpsubg)

    use ModSubgroupProp, only : cpsubgrstring, GroupMW, ioncharge, NKTAB
    use ModAIOMFACvar, only : AllocModAIOMFACvar
    use ModMRpart, only : MRinteractcoeff
    use ModSRunifac, only : SRsystm
    use ModComponentNames, only : names_mix

    implicit none
    !Interface variables:
    integer,intent(in) :: ndi, ninputcomp    
    integer,dimension(:),intent(in) :: compID           
    integer,dimension(:,:),intent(in) :: cpsubg                  !dimension(ninputcomp,topsubno), list of component subgroups and corresponding quantity (e.g. from input file or from dataset_components)
    !Local variables:
    integer :: i, II, J, JJ, K, Icheck, countall, q, qq, ID
    integer,dimension(200) :: SolvSubs2                          !temporary array for solvent subgroups (we allow a maximum of 200)
    integer,dimension(ninputcomp*2) :: ElectSubs2
    integer,dimension(201:topsubno) :: ElectPos
    logical :: already1, already2
    !--------------------------------------

    !initialize arrays and parameters:
    nd = ndi  !nd = 1 for web-version (single data file / mixture only)
    solvmixrefnd = .false.
    errorflagmix = 0
    if (allocated(ITAB)) then
        deallocate(ITAB, ITABsr, ITABMG, ITAB_dimflip, compN, Imaingroup, maingrindexofsubgr)
    endif
    allocate( ITAB(ninputcomp, topsubno), ITAB_dimflip(topsubno, ninputcomp), compN(ninputcomp) )
    !assign the different mixture components from interface input data to ITAB:
    ITAB = cpsubg
    !!do i = 1,topsubno
    !!    ITAB(:,i) = cpsubg(:,i) !transfer data to module array ITAB
    !!enddo
    ITAB_dimflip = TRANSPOSE(ITAB)
    CompN = compID

    !count the number of neutral components, nneutral:
    nneutral = 0
    Icheck = 0
    do i = 1,ninputcomp
        do J = 1,200
            if (ITAB_dimflip(J,i) /= 0) then
                if (Icheck /= i) then
                    nneutral = nneutral+1
                    if (i > nneutral) then
                        write(*,*) "ERROR: there are skipped neutral component numbers! Please check your input."
                        read(*,*)
                        STOP
                    endif
                endif
                Icheck = i
            endif
        enddo
    enddo

    !count the number of electrolyte (input) components, nelectrol.
    !nelectrol will be updated below to account for additional potential electrolytes that may form from different cation-anion combinations!
    nelectrol = 0
    Icheck = 0
    do i = nneutral+1,ninputcomp
        do J = 201,topsubno
            if (ITAB_dimflip(J,i) /= 0) then
                if (Icheck /= i) then
                    nelectrol = nelectrol+1 !count the same ion only once when existing in different components
                    if (i-nneutral > nelectrol) then
                        write(*,*) "ERROR: there are skipped salt component numbers! Please check your input."
                        read(*,*)
                        STOP
                        write(*,*) ""
                    endif
                endif
                Icheck = i
            endif
        enddo
    enddo

    !count the number of subgroups, NG, and other group numbers and properties:  
    Icheck = -1
    NG = 0
    NGN = 0
    NGS = 0
    SolvSubs2 = 0
    ElectSubs2 = 0
    countall = 0
    K = 1
    do J = 1,200
        do i = 1,ninputcomp
            if (ITAB(i,J) /= 0) then
                countall = countall+1
                if (Icheck /= J) then
                    SolvSubs2(K) = J !SolvSubs are the different solvent subgroups present in the mixture
                    K = K+1
                    NGN = NGN+1 !count the same subgroups only once, even if they are defined in different components
                endif
                Icheck = J
            endif
        enddo
    enddo
    if (allocated(SolvSubs)) then
        deallocate(SolvSubs)
    endif
    allocate( SolvSubs(NGN) )
    SolvSubs(1:NGN) = SolvSubs2(1:NGN)
    K = 1
    do J = 201,topsubno
        do i = 1,ninputcomp
            if (ITAB(i,J) /= 0) then
                countall = countall+1
                if (Icheck /= J) then
                    ElectSubs2(K) = J !the electrolyte subgroups (ions) present in the mixture
                    K = K+1
                    NGS = NGS+1 !count the same ions of salts / electrolytes only once, even if they exist in different components
                endif
                Icheck = J
            endif
        enddo
    enddo

    !check and set switch for presence/absence of water as a component in the system:
    if (any(SolvSubs(1:NGN) == 16)) then !water is present
        waterpresent = .true.
        do i = 1,ninputcomp
            if (ITAB(i,16) == 1) then
                CompN(i) = 401 !use for cases where input files is used to assign components
                exit
            endif
        enddo
    else
        waterpresent = .false.
    endif

    !now bring the subgroups into the right positions according to the components in ITAB (so that anions and cations of neutral electrolyte formula units are coupled):
    !the first anion in ElectSubs should correspond to the first cation in array ElectSubs.
    q = 1
    ElectPos = 0
    do i = nneutral+1,ninputcomp
        do J = 1,NGS
            K = ElectSubs2(J)
            if (ITAB_dimflip(K,i) > 0 .AND. ElectPos(K) == 0) then  !to make sure that every elect. subgroup is only checked once.
                ElectPos(K) = q !the position in the later array of ElectSubs, so that the anions and cations in ElectSubs correspond.
                q = q+1
            endif
        enddo
    enddo
    if (allocated(ElectSubs)) then
        deallocate(ElectSubs, AllSubs)
    endif
    allocate( ElectSubs(NGS), AllSubs(NGN+NGS) )
    ElectSubs = 0
    AllSubs = 0
    do J = 1,NGS
        K = ElectSubs2(J)
        ElectSubs(ElectPos(K)) = K
    enddo

    !Define array AllSubs which contains all (different) subgroups present in the mixture:
    isPEGsystem = .false.
    do K = 1,NGN
        AllSubs(K) = SolvSubs(K)
        if (SolvSubs(K) == 154) then
            isPEGsystem = .true.
        endif
    enddo
    do K = 1,NGS
        AllSubs(NGN+K) = ElectSubs(K)
    enddo

    !------------------------------------------------------------------------------------------------------
    !Calculating Ication and Ianion: saves the different cations/anions in the arrays Ication/Ianion at the 
    !cation/anion Nr. array index position (in ascending order of the type nr.); Also assign cation and anion charges.
    NGI = max( count(ElectSubs(1:NGS) < 241), count(ElectSubs(1:NGS) > 240) )
    if (allocated(Ication)) then
        deallocate(Ication, Ianion, cationZ, anionZ)
    endif
    allocate( Ication(NGI), Ianion(NGI), cationZ(NGI), anionZ(NGI) )
    ICation = 0
    IAnion = 0
    CatNr = 0
    AnNr = 0
    cationZ = 0.0_wp
    anionZ = 0.0_wp
    do J = 1,NGS
        already1 = .false.
        already2 = .false.
        K = ElectSubs(J)
        if (K < 241 .AND. K > 200) then
            do qq = 1,NGS
                if (.NOT. already1) then
                    if (Ication(qq) == 0 .OR. Ication(qq) == K) then
                        Ication(qq) = K
                        CatNr(K) = qq
                        cationZ(qq) = real(Ioncharge(K), kind=wp)
                        already1 = .true.
                    endif
                endif
            enddo
        else if (K <= topsubno .AND. K > 240) then !K > 240
            do qq = 1,NGS
                if (.NOT. already2) then
                    if (Ianion(qq) == 0 .OR. Ianion(qq) == K) then
                        Ianion(qq) = K
                        AnNr(K) = qq
                        anionZ(qq) = real(Ioncharge(K), kind=wp)
                        already2 = .true.
                    endif
                endif
            enddo
        endif
    enddo
    
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
    Ncation = count(Ication(:) > 0)
    Nanion = count(Ianion(:) > 0)

    !define different binary cation--anion combinations as electrolyte components of the mixture:
    if (nelectrol > 0) then
        elpresent = .true.
        !initialize KVLE(cation ID, anion ID) parameters for gas-liquid equilibria of different electrolyte components composed of binary cation+anion pairs
        !FROM LANGE's Handbook of Chemistry 15th Edition, Gibbs Energy of Formation Data; with R* = 8.314459 [J K^-1 mol^-1]
        KVLE_298K = 0.0_wp
        KVLE_298K(205,243) = 7.236814145E8_wp  !HBr  7.235752425E8_wp
        KVLE_298K(205,242) = 1.986870884E6_wp  !HCl
        KVLE_298K(205,241) = 5.749865093E8_wp  !HF
        KVLE_298K(205,244) = 2.168075448E9_wp  !HI
        KVLE_298K(205,245) = 4.191055023E6_wp  !HNO3
        !KVLE_298K(205,204) = 1.07614415E11_wp !NH3  (i.e. NH3(g) + H+(aq) <--> NH4+(aq)
        !initialize the IAP vs. aw parameterization coefficients for semi-volatile electrolyte components considered:
        IAPcoeffs = 0.0_wp
        IAPcoeffs(205,243,1:3) = [4.341292E-2_wp, 1.119213E0_wp, 5.644364E-1_wp]    !HBr
        IAPcoeffs(205,242,1:3) = [4.391087E-2_wp, 1.196585E0_wp, 4.603556E-1_wp]    !HCl
        IAPcoeffs(205,245,1:3) = [4.8857E-2_wp,   1.3369E0_wp,   3.423756E-1_wp]    !HNO3
        !--
        call defElectrolytes(nneutral, NGS, nelectrol)
        !--
    else
        elpresent = .false.
        if (allocated(nuestoich)) deallocate(nuestoich)
        allocate( nuestoich(nneutral) )
        nuestoich = 1.0_wp
    endif

    NGI = max(Nanion, Ncation)  !max number of ion groups of either cation or anion kind for present input composition.
    NG = NGN + NGS !maximum number of molecular / ionic groups
    NKNpNGS = nneutral + NGS

    !------------------------------------------------------------
    !Calculation of ITABMG. ITABMG is the analogue to ITAB for the main groups of the neutral components.
    allocate( ITABMG(nneutral,Nmaingroups), Imaingroup(NGN), maingrindexofsubgr(NGN) )
    ITABMG = 0  !initialize
    Imaingroup = 0
    maingrindexofsubgr = 0
    do i = 1,nneutral
        do J = 1,NGN
            II = SolvSubs(J)
            JJ = NKTAB(II)
            ITABMG(i,JJ) = ITABMG(i,JJ) + ITAB_dimflip(II,i)
            Imaingroup(J) = JJ  !list of main groups (here not yet filtered and sorted)
        enddo
        !check input correctness with respect to OH groups and special CHn groups bonded to OH groups:
        if (ITABMG(i,68) > 0 .OR. ITABMG(i,52) > 0) then
            if (ITABMG(i,68) + ITABMG(i,52) > ITABMG(i,69)) then !this would indicate an input error
                errorflagmix = 13
            endif
        endif
    enddo

    !------------------------------------------------------------
    !Calculation of the neutral subgroups molar masses; GroupMW is a list of the molecular masses of all neutral subgroups in the mixture.
    if (allocated(SubGroupMW)) then
        deallocate(SubGroupMW)
    endif
    allocate( SubGroupMW(NGN) )
    SubGroupMW(1:NGN) = GroupMW(SolvSubs(1:NGN))*1.0E-3_wp !assign subgroup molar mass in [kg/mol]

    !Calculation of Imaingroup, the sorted list of neutral main groups:
    do i = 1,NGN-1
        do J = i+1,NGN
            if (Imaingroup(i) == Imaingroup(J)) then !combine the same main groups
                Imaingroup(J) = 0
            endif
        enddo
    enddo
    !filter now the empty (Imaingroup = 0) positions out:
    do i = 1,NGN-1
        if (Imaingroup(i) == 0) then
            do J = i+1,NGN
                if (Imaingroup(J) /= 0) then
                    Imaingroup(i) = Imaingroup(J)
                    Imaingroup(J) = 0
                    exit !exit the inner J-loop
                endif
            enddo
        endif
    enddo
    !assemble a list stating for each subgroup of the present mixture to which main group of the present mixture it belongs (e.g. used in subroutine LR_MR_activity):
    do j = 1,NGN !j-loop over subgroups
        JJ = NKTAB(SolvSubs(j)) !static main group number associated with subgroup
        !find the index number in current mixture of the main group with static number JJ:
        do i = 1,NGN
            if (Imaingroup(i) == JJ) then
                maingrindexofsubgr(j) = i
                exit !leave inner i-loop
            endif
        enddo
    enddo
    !------------------------------------------------------------

    !Calculation of ITABsr: ITABsr is similar to ITAB, but the salts are divided into cations and anions as components themselves for use within UNIFAC...
    !Transfer information of cations (first) then anions into ITABsr list
    allocate( ITABsr(NKNpNGS,topsubno) )
    ITABsr = 0 !initialize
    K = nneutral
    JJ = nneutral+Ncation
    do J = 1,NGS
        ID = ElectSubs(J)
        !!i = sum(maxloc(ITAB(nneutral+1:nneutral+nelectrol, ID))) +nneutral
        if (ID > 200 .AND. ID < 240) then !cation
            K = K+1
            ITABsr(K,ID) = 1   !transfer cation info; just set to one ion, actual amount is coming from molality of the ion
        else if (ID > 240) then !anion
            JJ = JJ+1
            ITABsr(JJ,ID) = 1  !transfer anion info; just set to one ion, actual ion amount is coming from molality of the ion
        endif
    enddo
    !Transfer of information for neutrals from ITAB to ITABsr:
    do i = 1,nneutral
        do J = 1,NGN
            JJ = SolvSubs(J)
            ITABsr(i,JJ) = ITAB_dimflip(JJ,i)  !update ITABsr list.
        enddo
    enddo

    !save the number of independent components (ions combined to salts as given by ITAB):
    nindcomp = nneutral +nelectrol !number of components of the system (input components + additional electrolyte components from potential additional electroneutral ion combinations)

    if (allocated(OtoCratio)) then
        deallocate( OtoCratio, HtoCratio, compname, compnameTeX, ionname, ionnameTeX )
    endif
    allocate( OtoCratio(nindcomp), HtoCratio(nindcomp), compname(nindcomp), compnameTeX(nindcomp), ionname(NGS), ionnameTeX(NGS) )
    !set the name strings for the different independent components and the O:C ratio of the neutrals (when set already):
    call names_mix(CompN, compname, compnameTeX, ionname, ionnameTeX, OtoCratio, HtoCratio)     
    !now compname contains the name strings of all independent components of the present mixture nd;

    call cpsubgrstring() !generate for each component a subgroup-string stored in compsubgroups, compsubgroupsTeX, and compsubgroupsHTML

    !define the middle range (MR) coefficients of the actual mixture:
    call MRinteractcoeff()

    !calculate some short range (SR) parameters which don't change for the same mixture:
    call SRsystm(NKNpNGS)

    !allocate and populate array Mmass with the molar masses in component order.
    if (allocated(Mmass)) then
        deallocate(Mmass)
    endif
    allocate(Mmass(nindcomp))
    !calculate the molar masses of the different mixture species:
    call SetMolarMass(Mmass)    !Mmass lists the molar mass in [kg/mol] of all mixture components
    if (any(Mmass(:) < 0.0_wp)) then
        !$OMP CRITICAL
        write(*,*) 'ERROR: A molar mass value is negative indicating a missing value set in ModSubgroupProp.'
        read(*,*)
        STOP
        !$OMP end CRITICAL
    endif

    !-- allocate several composition-dependent variables from module ModAIOMFACvar:
    call AllocModAIOMFACvar()

    end subroutine definemixtures
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
    !*   -> created:        2014-07-21                                                      *
    !*   -> latest changes: 2018-05-24                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    module subroutine defElectrolytes(nneutral, NGS, nelectrol)

    use ModSubgroupProp, only : Ioncharge, IonO2Cequiv

    implicit none
    !Interface variables:
    integer,intent(in) :: nneutral, NGS
    integer,intent(inout) :: nelectrol
    !Local variables:
    integer :: i, J, JJ, K, KK, q, zc, za, ion, iel, nelin, nnp1
    logical :: already1
    real(wp),parameter :: onethird = 1.0_wp/3.0_wp
    !------------------------------------------------------------------

    nnp1 = nneutral + 1
    !calculate all different electrolyte components that can be formed from the input electro-neutral
    !cation--anion combinations (e.g. for gas-particle partitioning of different electrolytes):
    iel = Ncation*Nanion !maximum number of possible electrolyte components just formed by one cation and anion each.
    if (allocated(ElectComps)) then
        deallocate( ElectComps, ElectNues, ElectVolatile, ElectO2Cequiv, K_el )
    endif
    allocate( ElectComps(iel,2), ElectNues(iel,2), ElectVolatile(iel), ElectO2Cequiv(iel), K_el(nnp1:nneutral+iel) )
    if (allocated(nuestoich)) deallocate(nuestoich)
    allocate( nuestoich(nneutral+iel) )

    !initialize arrays for electrolytes:
    ElectComps = 0
    ElectNues = 0
    ElectVolatile = .false.
    ElectO2Cequiv = 0.0_wp
    K_el = 0.0_wp
    nuestoich = 1.0_wp

    !First, list the input electrolyte components in the order of definition in ITAB:
    nelin = nelectrol !ITAB-defined electrolyte components;  ninput - nneutral
    iel = 0
    do i = 1,nelin
        KK = 0
        JJ = 0
        do J = 1,NGS
            !loop over ElectSubs(K)
            ion = ElectSubs(J)
            if (KK == 0 .AND. ion > 200 .AND. ion < 240) then
                if (ITAB_dimflip(ion,nneutral+i) > 0) then
                    KK = ion
                endif
            else if (JJ == 0 .AND. ion > 240 .AND. ion <= topsubno) then
                if (ITAB_dimflip(ion,nneutral+i) > 0) then
                    JJ = ion
                endif
            endif
            if (KK > 0 .AND. JJ > 0) then
                iel = iel+1
                ElectComps(iel,1) = KK
                ElectComps(iel,2) = JJ
                K_el(nneutral+iel) = KVLE_298K(KK,JJ) !assign gas-particle equilibrium constants of the electrolyte components based on the involved cations and anions and tabulated values in KVLE_298K
                exit !leave inner J-loop
            endif
        enddo
        if (KK == JJ) then !then they are both == 0
            cycle !continue with next i loop
        endif
        zc = Ioncharge(KK)
        za = abs(Ioncharge(JJ))
        !now attribute the stoichiometric cation and anion numbers to creat an electroneutral component:
        if (zc == za) then !1:1 electrolyte (both the same charge of +-1 or +-2)
            ElectNues(iel,1) = 1
            ElectNues(iel,2) = 1
            ElectO2Cequiv(iel) = (IonO2Cequiv(KK)*IonO2Cequiv(JJ))**0.5_wp
            nuestoich(nneutral+iel) = 2.0_wp
        else
            select case(zc)
            case(1) !2:1 electrolyte
                ElectNues(iel,1) = 2
                ElectNues(iel,2) = 1 !anion charge: -2
                ElectO2Cequiv(iel) = (IonO2Cequiv(KK)**2 *IonO2Cequiv(JJ))**onethird
                nuestoich(nneutral+iel) = 3.0_wp
            case(2) !1:2 electrolyte
                ElectNues(iel,1) = 1 !cation charge: +2
                ElectNues(iel,2) = 2
                ElectO2Cequiv(iel) = (IonO2Cequiv(KK)*IonO2Cequiv(JJ)**2)**onethird
                nuestoich(nneutral+iel) = 3.0_wp
            end select
        endif
        if (KK == 205 .AND. JJ /= 261 .AND. JJ /= 262 .AND. JJ /= 250 .AND. JJ /= 248 &
        & .OR. (KK == 204 .AND. JJ == 245)) then
            ElectVolatile(iel) = .true.
        endif
    enddo !i
    !Second, make all other possible electroneutral binary cation-anion combinations for electrolytes and add to the lists:
    i = iel
    do K = 1,Ncation
        KK = Ication(K)
        do J = 1,Nanion
            JJ = Ianion(J)
            !check if cation-anion combination is already assigned to a component:
            already1 = .false.
            do q = 1,i
                if (ElectComps(q,1) == KK .AND. ElectComps(q,2) == JJ) then
                    already1 = .true.
                endif
            enddo
            if (.NOT. already1) then !found new electrolyte component
                i = i+1
                ElectComps(i,1) = KK
                ElectComps(i,2) = JJ
                K_el(nneutral+i) = KVLE_298K(KK,JJ)
                zc = Ioncharge(KK)
                za = abs(Ioncharge(JJ))
                !now attribute the stoichiometric cation and anion numbers to creat an electroneutral component:
                if (zc == za) then !1:1 electrolyte (both the same charge of +-1 or +-2)
                    ElectNues(i,1) = 1
                    ElectNues(i,2) = 1
                    ElectO2Cequiv(i) = (IonO2Cequiv(KK)*IonO2Cequiv(JJ))**0.5_wp
                    nuestoich(nneutral+i) = 2.0_wp
                else
                    select case(zc)
                    case(1) !2:1 electrolyte
                        ElectNues(i,1) = 2
                        ElectNues(i,2) = 1 !anion charge: -2
                        ElectO2Cequiv(i) = (IonO2Cequiv(KK)**2 *IonO2Cequiv(JJ))**onethird
                        nuestoich(nneutral+i) = 3.0_wp
                    case(2) !1:2 electrolyte
                        ElectNues(i,1) = 1 !cation charge: +2
                        ElectNues(i,2) = 2
                        ElectO2Cequiv(i) = (IonO2Cequiv(KK)*IonO2Cequiv(JJ)**2)**onethird
                        nuestoich(nneutral+i) = 3.0_wp
                    end select
                endif
                if (KK == 205 .AND. JJ /= 261 .AND. JJ /= 262 .AND. JJ /= 250 .AND. JJ /= 248 &
                & .OR. (KK == 204 .AND. JJ == 245)) then
                    ElectVolatile(i) = .true.
                endif
            endif
        enddo !J
    enddo !K
    nelectrol = i !the total number of different possible electrolyte components of the system (including both input electrolytes and additional, electroneutral ion combinations)
    NG = NGN + NGS !maximum number of molecular / ionic groups

    end subroutine defElectrolytes
    !==========================================================================================================================


    !************************************************************************
    !*                                                                      *
    !*  Set the molar masses of the given mixture components;               *
    !*  saved in array MolarM --> Mmass (neutrals first, then electrolytes) *
    !*                                                                      *
    !************************************************************************
    pure module subroutine SetMolarMass(MolarM)

    use ModSubgroupProp, only : GroupMW, SMWC, SMWA

    implicit none

    real(wp),dimension(:),intent(out) :: MolarM  !Molar mass in [kg/mol]
    integer :: i, J, K
    !...............................................
    !SMWA     : list of the molar masses of the anions
    !SMWC     : list of the molar masses of the cations
    !GroupMW  : list of molar masses of the organic/water subgroups
    MolarM = 0.0_wp  !initialize array
    !calculate the neutral compounds' molar mass from the values of their subgroups listed in GroupMW array:
    do i = 1,NGN    !loop over all solvent subgroups
        J = SolvSubs(i)
        MolarM(1:nneutral) = MolarM(1:nneutral) + ITAB(1:nneutral,J)*GroupMW(J)*1.0E-3_wp  !compound molar mass in [kg/mol]
    enddo
    !electrolytes:
    do i = 1,nelectrol
        J = ElectComps(i,1)-200 !the cation of this electrolyte in the list of cations as numbered in SMWC (1:40 instead of 201:240)
        K = ElectComps(i,2)-240 !the anion "
        MolarM(nneutral+i) = ( ElectNues(i,1)*SMWC(J) + ElectNues(i,2)*SMWA(K) )*1.0E-3_wp ! Molar mass of electrolyte i in [kg/mol]
    enddo

    end subroutine SetMolarMass
!==========================================================================================================================
    
end submodule SubModDefSystem