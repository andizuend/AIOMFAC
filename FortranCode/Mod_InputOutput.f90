!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module containing input (file) data reading subroutines and output processing      *
!*   subroutines, including HTML output for the web model.                              *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andreas Zuend,                                                                     *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2021-07-26                                                      *
!*   -> latest changes: 2023-03-18                                                      *
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
!*   -  subroutine Output_TXT                                                           *
!*   -  subroutine Output_HTML                                                          *
!*   -  subroutine RepErrorWarning                                                      *
!*   -  subroutine ReadInputFile                                                        *
!*                                                                                      *
!****************************************************************************************
module Mod_InputOutput

use Mod_NumPrec, only : wp

implicit none
private

public :: Output_TXT, Output_HTML, ReadInputFile, RepErrorWarning

!============================================================================================
    contains
!============================================================================================

    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Read an input text file with information about the subgroups of each system        *
    !*   component and the compositions and temperatures for specific mixture calculations. *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
    !*                                                                                      *
    !*   -> created:        2011                                                            *
    !*   -> latest changes: 2019-10-17                                                      *
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
    !****************************************************************************************
    subroutine ReadInputFile(filepath, folderpathout, filename, ninpmax, maxpoints, unito, verbose, &
        & ncp, npoints, warningind, errorind, filevalid, cpnameinp, cpsubg, T_K, composition, xinputtype)

    use ModSystemProp, only : topsubno

    implicit none

    !interface variables:
    character(len=3000),intent(inout) :: filepath, folderpathout
    character(len=200),intent(inout) :: filename
    integer,intent(in) :: ninpmax, maxpoints
    integer,intent(inout) :: unito
    logical,intent(in) :: verbose
    integer,intent(out) :: ncp, npoints
    integer,intent(inout) :: warningind, errorind
    logical,intent(out) :: filevalid
    character(len=200),dimension(:),intent(out) :: cpnameinp    !list of assigned component names (from input file)
    integer,dimension(:,:),intent(out) :: cpsubg                !list of input component subgroups and corresponding subgroup quantities
    real(wp),dimension(:),intent(out) :: T_K                    !temperature of data points in Kelvin
    real(wp),dimension(:,:),intent(out) :: composition          !array of mixture composition points for which calculations should be run
    logical,intent(out) :: xinputtype
    !--
    !local variables:
    character(len=:),allocatable :: cn                          !this assumes a maximum four-digit component number in the system (max. 9999); to be adjusted otherwise.
    character(len=4) :: dashes, equalsigns, pluses
    character(len=20) :: dummy, cnformat
    character(len=50) :: txtcheck, inpfolder, outpfolder
    character(len=3000) :: errlogfile, fname
    character(len=20),dimension(ninpmax) :: txtarray
    integer :: cpno, i, inpfilesize, istat, k, kinpf, qty, subg, unitx
    logical :: fileopened, fileexists
    !...................................................................................

    !initialize variables
    ncp = 0
    cpsubg = 0
    cpnameinp = "none"
    composition = 0.0_wp
    T_K = 298.15_wp      !default temperature in [K]
    dashes = "----"
    equalsigns = "===="
    pluses = "++++"
    fileexists = .false.
    filevalid = .false.

    !==== read INPUT data ===========================================================

    !extract filepath from the input filepath (which could include a path to the directory):
    k = len_trim(filepath)
    if (k < 1) then !no file was provided
        write(*,*) ""
        write(*,*) "ERROR from AIOMFAC-web: no input file was provided!"
        write(*,*) "An input file path needs to be provided via command line argument 2"
        write(*,*) ""
        read(*,*)
        stop
    endif

    i = index(filepath, "/input_", BACK = .true.)  !find index position of input file, assuming that all input files start with "input_"
    if (i < k .and. i > 0) then !found a direcory path in UNIX file path style
        filename = trim(filepath(i+1:k))
        filepath = trim(filepath(1:k))
        kinpf = index(filepath(1:i-1), "/", BACK = .true.)  !find folder length of folder containing the input file
    else
        kinpf = 0
        i = index(filepath, "\input_", BACK = .true.) 
        if (i > 0 .and. i < k) then !found a direcory path in WINDOWS file path style
            filename = trim(filepath(i+1:k))
            filepath = trim(filepath(1:k))
            kinpf = index(filepath(1:i-1), "\", BACK = .true.)  !find folder length of folder containing the input file
        endif
    endif
    if (kinpf > i-1) then
        kinpf = 0
    endif
    inpfolder = trim(filepath(kinpf+1:i-1))
    kinpf = len_trim(inpfolder)

    !check / create for associated output directory:
    outpfolder = "Outputfiles/"  !initialize
    istat = index(inpfolder, "inp")
    if (istat > 0 .and. istat < kinpf) then
        outpfolder = "outp"//trim(inpfolder(4:))//"/"  
    else
        istat = index(inpfolder, "Inp")
        if (istat > 0) then
            outpfolder = "Outp"//trim(inpfolder(4:))//"/"   
        endif
    endif
    folderpathout = trim(filepath(1:i-(kinpf+1)))//trim(outpfolder)

    !use the name of the input file to create a corresponding output file name:
    i = index(filename, ".txt") !returns starting position of string input within string filename
    !create an error-logfile associated with the input file name:
    errlogfile = "Errorlog_"//filename(i-4:)
    fname = trim(folderpathout)//trim(errlogfile)
    OPEN (NEWUNIT = unito, FILE = fname, STATUS ='UNKNOWN') !unito is the error / logfile unit
    !-----
    !check if file exists and read its content if true:
    fname = trim(filepath)
    INQUIRE(FILE = fname, EXIST = fileexists, size = inpfilesize) !inpfilesize is the file size in [bytes]
    !Delete very large files that can only mean uploaded spam content and not actual input:
    if (fileexists) then
        if (real(inpfilesize, kind=wp) > 50.0_wp*(ninpmax +ninpmax*maxpoints) ) then !likely not a valid input file
            fname = trim(filepath)
            OPEN (NEWUNIT = unitx, FILE = fname, IOSTAT=istat, ACTION='read', STATUS='OLD')
            read(unitx,*) dummy, dummy, dummy, txtcheck
            close(unitx)
            if (.not. (txtcheck(1:11) == "AIOMFAC-web")) then !invalid file (likely spam)
                OPEN (NEWUNIT = unitx, FILE = fname, STATUS='OLD')
                close(unitx, STATUS='DELETE')  !close and delete the file
                fileexists = .false.
            endif
        endif
    endif

    !attempt to read a supposedly valid input file:
    if (fileexists) then
        if (verbose) then
            write(unito,*) ""
            write(unito,*) "MESSAGE from AIOMFAC: input file exists."
        endif
        fname = trim(filepath)
        OPEN (NEWUNIT = unitx, FILE = fname, IOSTAT=istat, ACTION='read', STATUS='OLD')
        if (istat /= 0) then ! an error occurred
            write(unito,*) ""
            write(unito,*) "MESSAGE from AIOMFAC: an error occurred while trying to open the file! IOSTAT: ", istat
        endif
        !validate file as a correct input text-file (no spam):
        read(unitx,*) dummy, dummy, dummy, txtcheck
        if (txtcheck(1:11) == "AIOMFAC-web") then
            filevalid = .true.
            if (verbose) then
                write(unito,*) ""
                write(unito,*) "MESSAGE from AIOMFAC: input file has passed the first line text validation and will be read."
            endif
            !read input data from file:
            BACKSPACE unitx !jump back to beginning of record (to the beginning of the line)
            read(unitx,*) dummy, dummy, dummy, dummy, dummy !read first line
            read(unitx,*)               !read empty line 2
            read(unitx,*) dummy, dummy  !read line 3
            read(unitx,*) dummy         !read line 4
        else !file not valid. It will be closed and deleted below
            filevalid = .false.
            if (verbose) then
                write(unito,*) ""
                write(unito,*) "MESSAGE from AIOMFAC: input file does not pass first line text validation and will be deleted."
            endif
        endif
        !loop over mixture components with variable numbers of subgroups to read (using inner loop):
        if (filevalid) then
            if (verbose) then
                write(unito,*) ""
                write(unito,*) "MESSAGE from AIOMFAC: reading component data from input file."
            endif
            k = max(2, CEILING(LOG10(real(ninpmax))) ) !determine order of magnitude digits
            write(dummy,'(I0)') k
            cnformat = "(I"//trim(dummy)//"."//trim(dummy)//")" !variable format specifier
            allocate(character(len=k) :: cn)
            do ncp = 1,ninpmax
                if (verbose .and. istat /= 0) then
                    write(unito,*) ""
                    write(unito,*) "MESSAGE from AIOMFAC: file end found in input file at ncp = ", ncp
                endif
                read(unitx,*,IOSTAT=istat) txtcheck !read only first argument on this line for subsequent check
                if (txtcheck(1:4) == pluses .or. istat /= 0) then !"++++" indicates end of this components definition part
                    EXIT !exit ncp do-loop
                else !in this case, argument 3 of txtcheck is the component no.:
                    BACKSPACE unitx !jump back to beginning of record (to the beginning of the line)
                    read(unitx,*) dummy, dummy, cpno !read line with component's number
                endif
                read(unitx,*) dummy, dummy, cpnameinp(cpno) !read line with component's name
                if (len_trim(cpnameinp(cpno)) < 1) then !no component name was assigned, generate a default name
                    write(cn,cnformat) cpno
                    cpnameinp(cpno) = "cp_"//cn
                endif
                if (verbose) then
                    write(unito,*) ""
                    write(unito,*) "MESSAGE from AIOMFAC: found a component cpno =", cpno
                endif
                do !until exit
                    read(unitx,*,IOSTAT=istat) txtcheck !read only first argument on this line
                    !check whether another subgroup is present or not
                    if (txtcheck(1:4) == dashes .or. istat /= 0) then !"----" indicates no more subgroups of this component
                        EXIT !leave the inner do-loop
                    else !else continue reading the next subgroup
                        BACKSPACE unitx
                        read(unitx,*) dummy, dummy, dummy, subg, qty  !continue reading line with subgroup no. and corresp. quantity
                        cpsubg(cpno,subg) = cpsubg(cpno,subg)+qty
                    endif
                enddo
            enddo
            ncp = ncp-1 !ncp-1 is the number of different components in this mixture
            if (verbose) then
                write(unito,*) ""
                write(unito,*) "MESSAGE from AIOMFAC: component data read."
                write(unito,*) "MESSAGE from AIOMFAC: number of components: ", ncp
            endif
            if (ncp == ninpmax) then
                read(unitx,*,IOSTAT=istat) txtcheck !read only first argument on this line for subsequent check
                if (txtcheck(1:4) /= pluses) then
                    errorind = 34
                    filevalid = .false.
                    write(*,*) "AIOMFAC ERROR 34: maximum number of input components reached while reading input file."
                    write(*,*) "AIOMFAC ERROR 34: check whether ninpmax value is too small."
                    write(unito,*) ""
                    write(unito,*) "AIOMFAC ERROR 34: maximum number of input components reached while reading input file."
                endif
            endif
            if (filevalid) then
                read(unitx,*) dummy, dummy     !read mixture composition line
                read(unitx,*) dummy, dummy, i  !read mass fraction? line
                if (i == 1) then !composition in mass fractions
                    xinputtype = .false.
                else !composition in mole fractions
                    xinputtype = .true.
                endif
                read(unitx,*) dummy, dummy, i !read mole fraction? line
                read(unitx,*) dummy !read dashes line
                !now read the lines with the mixture composition points:
                read(unitx,*) txtarray(1:ncp) !read header line of composition table
                do npoints = 1,maxpoints !or until exit
                    read(unitx,*,IOSTAT=istat) txtcheck !read only the first argument on this line
                    if (txtcheck(1:4) == equalsigns .or. istat /= 0) then !"====" indicates no more composition points (and last line of input file)
                        EXIT !leave the do-loop (normal exit point)
                    else if (IACHAR(txtcheck(1:1)) > 47 .and. IACHAR(txtcheck(1:1)) < 58) then !validate that the data is actual intended input and not some sort of text field spam).
                        BACKSPACE unitx
                        read(unitx,*) txtcheck, dummy !read only the first two arguments on this line
                        if (IACHAR(dummy(1:1)) > 47 .and. IACHAR(dummy(1:1)) < 58) then !it is a number
                            BACKSPACE unitx
                            read(unitx,*) i, T_K(npoints), composition(npoints,2:ncp) !read the temperature in [K] and composition values of the components [2:ncp] into the array
                        else
                            if (npoints == 1) then
                                filevalid = .false. !file is not valid due to incorrect input in text field
                                EXIT
                            else
                                warningind = 31
                                EXIT !abnormal exit point of loop due to non-number entries at a certain line in the text field
                            endif
                        endif
                    else
                        if (npoints == 1) then
                            filevalid = .false. !file is not valid due to incorrect input in text field
                            EXIT
                        else
                            warningind = 31
                            EXIT !abnormal exit point of loop due to non-number entries at a certain line in the text field
                        endif
                    endif
                enddo
                npoints = npoints-1 !the number of composition points
                if (verbose) then
                    write(unito,*) ""
                    write(unito,*) "MESSAGE from AIOMFAC: composition points read."
                    write(unito,*) "MESSAGE from AIOMFAC: number of points: ", npoints
                endif
                if (.not. filevalid) then
                    !close and delete file from server:
                    close(unitx, STATUS = 'DELETE')
                else
                    close(unitx)
                endif
                !assign component 01 its composition values based on the rest of the component's data:
                do i = 1,npoints
                    composition(i,1) = 1.0_wp - sum(composition(i,2:ncp))
                enddo
            endif !filevalid2
        endif !filevalid1
    else
        write(unito,*) ""
        write(unito,*) "ERROR in AIOMFAC: Input file does not exist at expected location: ", trim(filename)
        write(unito,*) ""
    endif !fileexists

    if (fileexists) then
        if (.not. filevalid) then !input file contains errors or is completely invalid (submitted spam etc.)
            if (verbose) then
                write(unito,*) ""
                write(unito,*) "MESSAGE from AIOMFAC: input file did not pass full validation and may be a spam file."
                write(unito,*) "MESSAGE from AIOMFAC: the input file will be deleted to prevent spam files and malicious code from occupying the server."
                write(unito,*) ""
            endif
            INQUIRE(FILE = fname, OPENED = fileopened)
            if (errorind == 0) then
                errorind = 32
                !close and delete file from server:
                if (fileopened) then
                    close(unitx, STATUS = 'DELETE')
                endif
            else
                if (fileopened) then
                    close(unitx)
                endif
            endif
        else
            !check for any valid points for calculations before initializing AIOMFAC:
            if (npoints < 1) then !no points input for calculations
                errorind = 33 !no input in text field..
                if (verbose) then
                    write(unito,*) ""
                    write(unito,*) "MESSAGE from AIOMFAC: no composition points have been entered in the text field. There is nothing to calculate."
                    write(unito,*) ""
                endif
                filevalid = .false.
            endif
        endif
    endif

    end subroutine ReadInputFile
    !========================================================================================
    
        
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Create an output ASCII text file with an overall mixture data header and           *
    !*   individual tables for all components / species (in case of ions).                  *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
    !*                                                                                      *
    !*   -> created:        2011                                                            *
    !*   -> latest changes: 2021-12-08                                                      *
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
    !****************************************************************************************
    subroutine Output_TXT(fname, VersionNo, nspecmax, npoints, watercompno, cpnameinp, T_K, &
        & px, out_data, out_viscdata)

    !module variables:
    use ModSystemProp, only : compname, compsubgroups, compsubgroupsTeX, idCO2, NGS, NKNpNGS, &
        & ninput, nneutral
    use ModSubgroupProp, only : subgrname, subgrnameTeX

    implicit none
    !interface variables:
    character(len=*),intent(in) :: fname 
    character(len=*),intent(in) :: VersionNo
    integer,intent(in) ::  nspecmax, npoints, watercompno
    character(len=200),dimension(nspecmax),intent(in) :: cpnameinp   !list of assigned component names (from input file)
    real(wp),dimension(npoints),intent(in) :: T_K
    integer,dimension(nspecmax),intent(inout) :: px
    real(wp),dimension(7,npoints,NKNpNGS),intent(in) :: out_data
    real(wp),dimension(3,npoints),intent(in) :: out_viscdata
    !--
    !local variables:
    character(len=:),allocatable :: cn
    character(len=5) :: tlen
    character(len=50) :: subntxt, Iformat
    character(len=150) :: cnformat, horizline, txtn, tablehead
    character(len=3000) :: txtsubs, mixturestring
    integer ::  i, k, kms, pointi, qty, unitx
    real(wp) :: RH
    !...................................................................................

    k = max(2, CEILING(LOG10(real(nspecmax))) )  
    write(tlen,'(I0)') k
    Iformat = "I"//trim(tlen)//"."//trim(tlen)  !variable integer specifier
    cnformat = "("//Iformat//")"                !constructed format specifier
    allocate(character(len=k) :: cn)

    !create a character string of the mixture as a series of its components:
    kms = len(mixturestring)
    mixturestring = trim(cpnameinp(1)) !first component
    !loop over all further components / species:
    do i = 2,nspecmax
        px(i) = max(1, px(i))
        if (INT(out_data(6,px(i),i)) == 0) then !neutral component
            txtn = trim(adjustl(cpnameinp(i)))
        else !ion (with its own link)
            if (INT(out_data(6,px(i),i)) > 1) then
                txtn = trim( adjustl(subgrname(INT(out_data(6,px(i),i) ))) )
            else
                txtn = "unknown_sub"
            endif
        endif
        k = len_trim(mixturestring) +len_trim(' + '//trim(txtn) )
        if (k < kms - 50) then
            mixturestring = trim(mixturestring)//' + '//trim(adjustl(txtn))
        else
            qty = nspecmax -i
            write(subntxt,'(I0)') qty
            mixturestring = trim(mixturestring)//' + '//trim(subntxt)//' additional components ...'
            EXIT
        endif
    enddo
    mixturestring = adjustl(mixturestring)

    OPEN (NEWUNIT = unitx, FILE = fname, STATUS = "UNKNOWN")
    write(unitx,'(A)') "==========================================================================================================="
    write(unitx,'(A)') "AIOMFAC-web, version "//VersionNo
    write(unitx,'(A)') "==========================================================================================================="
    write(unitx,*) ""
    mixturestring = "'"//trim(mixturestring)//"'"
    ![Note that Intel Fortran allows for variable format statements like A<len_trim(txtn)>, but that is a non-standard extension not supported by gfortran and other compilers - therefore avoided here.]
    write(unitx, '(A,A)') "Mixture name:  ", trim(mixturestring)
    write(unitx, '(A,I0.2)') "Number of independent input components: ", ninput
    write(unitx, '(A,I0.2)') "Number of different neutral components: ", nneutral
    write(unitx, '(A,I0.2)') "Number of different inorganic ions    : ", NGS
    write(unitx,*) ""
    write(unitx,'(A)') "The AIOMFAC output is tabulated for each component/species individually below."
    write(unitx,*) ""
    write(unitx,'(A)') '---  Table key  -------------------------------------------------------------------------------------------'
    write(unitx,'(A)') 'no.              : composition point number;                                                        '
    write(unitx,'(A)') 'T [K]            : absolute temperature;                                                            '
    write(unitx,'(A)') 'RH [%]           : relative humidity in equilibrium with the liquid mixture (bulk conditions);      '
    write(unitx,'(A)') 'w(j) [-]         : weight fraction (mass fraction) of species "j";                                  '
    write(unitx,'(A)') 'x_i(j) [-]       : mole fraction of species "j", calculated on the basis of completely              '
    write(unitx,'(A)') '                   dissociated inorganic ions; exception: the partial dissociation of bisulfate     '
    write(unitx,'(A)') '                   (HSO4- <--> H+ + SO4--) is explicitly considered when present in the mixture;    '
    write(unitx,'(A)') 'm_i(j) [mol/kg]  : molality of species "j" [mol(j)/(kg solvent mixture)], where "solvent mixture"   '
    write(unitx,'(A)') '                   refers to the electrolyte-free mixture (water + organics);                       '
    write(unitx,'(A)') 'a_coeff_x(j) [-] : activity coefficient of "j", defined on mole fraction basis (used for non-ionic  '
    write(unitx,'(A)') '                   components) with pure (liquid) component "j" reference state;                    '
    write(unitx,'(A)') 'a_coeff_m(j) [-] : activity coefficient of "j", defined on molality basis (used for inorg. ions)    '
    write(unitx,'(A)') '                   with reference state of infinite dilution of "j" in pure water;                  '
    write(unitx,'(A)') 'a_x(j) [-]       : activity of "j", defined on mole fraction basis (pure component reference state);'
    write(unitx,'(A)') 'a_m(j) [-]       : activity of "j", defined on molality basis (used for inorg. ions) with reference '
    write(unitx,'(A)') '                   state of infinite dilution of "j" in pure water;                                 '
    write(unitx,'(A)') 'log10(eta/[Pa.s]): base-10 log of the predicted dynamic viscosity of the mixture;                   '
    write(unitx,'(A)') '+/-log10(eta sens./[Pa.s]): base-10 log of the estimated sensitivity/uncertainty of the dynamic     '
    write(unitx,'(A)') '                   viscosity; see details on the "Hints & Examples" webpage;                        '
    write(unitx,'(A)') 'flag             : error/warning flag, a non-zero value (error/warning number) indicates that a     '
    write(unitx,'(A)') '                   numerical issue or a warning occurred at this data point, suggesting evaluation  '
    write(unitx,'(A)') '                   with caution (warnings) or exclusion (errors) of this data point.                '
    write(unitx,'(A)') '-----------------------------------------------------------------------------------------------------------'
    write(unitx,*) ''
    write(unitx,*) ''
    !--
    horizline = "-----------------------------------------------------------------------------------------------------------"
    !--
    !data table for mixture viscosity output
    write(unitx,'(A)') "Properties of this phase: mixture viscosity"
    !write table column headers:
    write(unitx,'(A)') adjustl(horizline)
    write(unitx,'(2X, A)') "no.   T_[K]     RH_[%]    log10(eta/[Pa.s])   +/-log10(eta sens./[Pa.s])     flag "
    write(unitx,'(A)') adjustl(horizline)
    !write data to viscosity table
    do pointi = 1,npoints  !loop over composition points
        if (watercompno > 0) then
            RH = out_data(5,pointi,watercompno)*100.0_wp !RH in %
            if (RH > 1000.0_wp .or. RH < 0.0_wp) then
                RH = -99.99_wp
            endif
        else
            RH = 0.0_wp
        endif
        write(unitx,'(I5.3,2X,F7.2,2X,F7.2,3X,2(ES12.5,10X),3X,I2)') pointi, T_K(pointi), RH, out_viscdata(1,pointi), out_viscdata(2,pointi), INT(out_viscdata(3,pointi))
    enddo !pointi
    write(unitx,'(A)') adjustl(horizline)
    write(unitx,*) ""

    !write individual data tables for each component / ionic species:
    do i = 1,nspecmax
        write(cn,cnformat) i !component / species number as character string
        write(unitx,*) ""
        !distinguish between neutral components and ionic species:
        if (INT(out_data(6,px(i),i)) == 0) then !neutral component
            write(unitx,'(A,I0.2)') "Mixture's component # : ", i
            txtn = "'"//trim(adjustl(compname(i)))//"'"
            write(unitx, '(A,A)')  "Component's name      : ", trim(txtn)
            txtsubs = "'"//trim(compsubgroups(i))//"'"
            write(unitx, '(A,A)') "Component's subgroups : ", trim(txtsubs)
            txtsubs = "'"//trim(compsubgroupsTeX(i))//"'"
            write(unitx, '(A,A)') "Subgroups, TeX format : ", trim(txtsubs)
            !write table column headers:
            write(unitx,'(A)') trim(horizline)
            if (i == idCO2) then  !check exception: CO2 has the activity defined on molality basis like ions:
                tablehead = "no.   T_[K]     RH_[%]   w("//cn//")          x_i("//cn//")        m_i("//cn//")        a_coeff_m("//cn//")   a_m("//cn//")        flag "
            else        !regular case
                tablehead = "no.   T_[K]     RH_[%]   w("//cn//")          x_i("//cn//")        m_i("//cn//")        a_coeff_x("//cn//")   a_x("//cn//")        flag "
            endif
            write(unitx, '(2X,A)') trim(tablehead)
            !--
        else if (INT(out_data(6,px(i),i)) < 240) then !cation
            write(unitx,'(A,I0.2)') "Mixture's species, #  : ", i
            subntxt = trim(adjustl(subgrname(INT(out_data(6,px(i),i)))))
            qty = len_trim(subntxt)
            txtn = adjustl(subntxt(2:qty-1)) !to print the ion name without the enclosing parenthesis ()
            txtn = "'"//trim(txtn)//"'"
            write(unitx, '(A,A)')  "Cation's name         : ", trim(txtn)
            txtsubs = "'"//trim(subntxt)//"'"
            write(unitx,'(A,A)') "Cation's subgroups    : ", trim(txtsubs)
            subntxt = trim(adjustl(subgrnameTeX(INT(out_data(6,px(i),i)))))
            txtsubs = "'"//trim(subntxt)//"'"
            write(unitx, '(A,A)') "Subgroups, TeX format : ", trim(txtsubs)
            !write table column headers:
            write(unitx, '(A)') trim(horizline)
            tablehead = "no.   T_[K]     RH_[%]   w("//cn//")          x_i("//cn//")        m_i("//cn//")        a_coeff_m("//cn//")   a_m("//cn//")        flag "
            write(unitx, '(2X,A)') trim(tablehead)
            !--
        else if (INT(out_data(6,px(i),i)) > 240) then !anion
            write(unitx,'(A,I0.2)') "Mixture's species, #  : ", i
            subntxt = trim( adjustl( subgrname(INT(out_data(6,px(i),i))) ) )
            qty = len_trim(subntxt)
            txtn = adjustl(subntxt(2:qty-1))
            txtn = "'"//trim(txtn)//"'"
            write(unitx, '(A,A)')  "Anion's name          : ", trim(txtn)
            txtsubs = "'"//trim(subntxt)//"'"
            write(unitx, '(A,A)') "Anion's subgroups     : ", trim(txtsubs)
            subntxt = trim( adjustl( subgrnameTeX(INT(out_data(6,px(i),i))) ) )
            txtsubs = "'"//trim(subntxt)//"'"
            write(unitx, '(A,A)') "Subgroups, TeX format : ", trim(txtsubs)
            !write table column headers:
            write(unitx,'(A)') trim(horizline)
            tablehead = "no.   T_[K]     RH_[%]   w("//cn//")          x_i("//cn//")        m_i("//cn//")        a_coeff_m("//cn//")   a_m("//cn//")        flag "
            write(unitx, '(2X, A)') trim(tablehead)
            !--
        else
            !error
        endif
        !--
        write(unitx,'(A)') adjustl(horizline)
        !write data to table:
        do pointi = 1,npoints  !loop over composition points
            if (watercompno > 0) then
                RH = out_data(5,pointi,watercompno)*100.0_wp !RH in %
                if (RH > 1000.0_wp .or. RH < 0.0_wp) then
                    RH = -99.99_wp
                endif
            else
                RH = 0.0_wp
            endif
            write(unitx,'(I5.3,2X,F7.2,2X,F7.2,3X,5(ES12.5,3X),3X,I2)') pointi, T_K(pointi), RH, out_data(1:5,pointi,i), INT(out_data(7,pointi,i))
        enddo !pointi
        write(unitx,'(A)') trim(horizline)
        write(unitx,*) ""
    enddo
    write(unitx,*) ""
    !--
    write(unitx,*) ""
    write(unitx,'(A)') "==========================================================================================================="
    close(unitx)
        
    end subroutine Output_TXT
    !========================================================================================
        
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Create an output file in HTML markup format with an overall mixture data header    *
    !*   and individual tables for all components / species (in case of ions).              *
    !*   This output format is used for showning the results on the AIOMFAC-web webpage.    *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
    !*                                                                                      *
    !*   -> created:        2011                                                            *
    !*   -> latest changes: 2021-12-05                                                      *
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
    !****************************************************************************************
    subroutine Output_HTML(fname, VersionNo, nspecmax, npoints, watercompno, cpnameinp, T_K, &
        & px, out_data, out_viscdata)

    !module variables:
    use ModSystemProp, only : compname, compsubgroupsHTML, idCO2, NGS, NKNpNGS, ninput, nneutral
    use ModSubgroupProp, only : subgrnameHTML

    implicit none
    !interface variables:
    character(len=*),intent(in) :: fname 
    character(len=*),intent(in) :: VersionNo
    integer,intent(in) ::  nspecmax, npoints, watercompno
    character(len=200),dimension(nspecmax),intent(in) :: cpnameinp      !list of assigned component names (from input file)
    real(wp),dimension(npoints),intent(in) :: T_K
    integer,dimension(nspecmax),intent(inout) :: px
    real(wp),dimension(7,npoints,NKNpNGS),intent(in) :: out_data
    real(wp),dimension(3,npoints),intent(in) :: out_viscdata
    !--
    !local variables:
    character(len=:),allocatable :: cn
    character(len=5) :: tlen
    character(len=50) :: subntxt, Iformat
    character(len=150) :: cnformat, tformat, txtn
    character(len=400) :: outtxtleft
    character(len=3000) :: txtsubs, mixturestring
    integer ::  i, k, kms, pointi, qty, unitx
    real(wp) :: RH
    !...................................................................................

    k = max(2, CEILING(LOG10(real(nspecmax))) )
    write(tlen,'(I0)') k
    Iformat = "I"//trim(tlen)//"."//trim(tlen)      !variable integer specifier
    cnformat = "("//Iformat//")"                    !constructed format specifier
    allocate(character(len=k) :: cn)
    !--
    OPEN (NEWUNIT = unitx, FILE = fname, STATUS = "UNKNOWN")
    outtxtleft = adjustl("<h3>AIOMFAC-web, version "//VersionNo//"</h3>")
    write(unitx,'(A)') outtxtleft
    !create a character string of the mixture as a series of its components and add links to the components:
    kms = len(mixturestring)
    mixturestring = '<a href="#'//trim(adjustl(compname(1)))//'">'//trim(adjustl(compname(1)))//'</a>' !first component
    !loop over all further components / species:
    do i = 2,nspecmax
        if (INT(out_data(6,px(i),i)) == 0) then !neutral component
            txtn = '<a href="#'//trim(adjustl(compname(i)))//'">'//trim(adjustl(compname(i)))//'</a>'
        else !ion (with its own link)
            txtn = '<a href="#'//trim(adjustl(subgrnameHTML(INT(out_data(6,px(i),i)))))//'">'//trim(adjustl(subgrnameHTML(INT(out_data(6,px(i),i) ))))//'</a>'
        endif
        k = len_trim(mixturestring) +len_trim(' + '//trim(txtn) )
        if (k < kms - 50) then
            mixturestring = trim(mixturestring)//' + '//trim(adjustl(txtn))
        else
            qty = nspecmax -i
            write(subntxt,'(I0)') qty
            mixturestring = trim(mixturestring)//' + '//trim(subntxt)//' additional components ...'
            EXIT
        endif
    enddo
    mixturestring = trim(mixturestring)
    write(unitx, '(A,A)') "<p>Mixture name: &nbsp;", adjustl(mixturestring)
    write(unitx, '(A,I0.2)') "<br> Number of independent input components: ", ninput
    write(unitx, '(A,I0.2)') "<br> Number of different neutral components: ", nneutral
    write(unitx, '(A,I0.2)') "<br> Number of different inorganic ions    : ", NGS
    write(unitx,*) "</p>"
    write(unitx,'(A)') adjustl('<p> The AIOMFAC output is tabulated for each component/species individually below. </p>')
    write(unitx,*) '<br>'
    !--
    write(unitx,'(A)') adjustl('<table class="tablekey">')
    write(unitx,'(A)') adjustl('<caption id="tablekeycapt"> Table key </caption>')
    write(unitx,'(A)') adjustl('<tr><td id="tablekeyc1">no.</td><td id="tablekeyc2">: </td><td id="tablekeyc3"> composition point number;</td></tr>')
    write(unitx,'(A)') adjustl('<tr><td>T [K]</td><td>: </td><td> absolute temperature;</td></tr>')
    write(unitx,'(A)') adjustl('<tr><td>RH [%]</td><td>: </td><td> relative humidity in equilibrium with the liquid mixture (bulk conditions);</td></tr>')
    write(unitx,'(A)') adjustl('<tr><td>w(j) [-]</td><td>: </td><td> weight fraction (mass fraction) of species "j";</td></tr>')
    write(unitx,'(A)') adjustl('<tr><td>x_i(j) [-]</td><td>: </td><td> mole fraction of species "j", calculated on the basis of completely dissociated &
        & inorganic ions; exception: the partial dissociation of bisulfate (<span class="chemf">HSO<sub>4</sub><sup>-</sup> &#8596; &
        & H<sup>+</sup> + SO<sub>4</sub><sup>2-</sup></span>) is explicitly considered when present in the mixture; </td></tr>')
    write(unitx,'(A)') adjustl('<tr><td>m_i(j) [mol/kg]</td><td>: </td><td> molality of species "j" [mol(j)/(kg solvent mixture)], where "solvent mixture" &
        & refers to the electrolyte-free mixture (water + organics); </td></tr>')
    write(unitx,'(A)') adjustl('<tr><td>a_coeff_x(j) [-]</td><td>: </td><td> activity coefficient of "j", defined on mole fraction basis (used for &
        & non-ionic components) with pure (liquid) component "j" reference state;</td></tr>')
    write(unitx,'(A)') adjustl('<tr><td>a_coeff_m(j) [-]</td><td>: </td><td> activity coefficient of "j", defined on molality basis (used for &
        & inorg. ions) with reference state of infinite dilution of "j" in pure water;</td></tr>')
    write(unitx,'(A)') adjustl('<tr><td>a_x(j) [-]</td><td>: </td><td> activity of "j", defined on mole fraction basis (pure component reference &
        & state);</td></tr>')
    write(unitx,'(A)') adjustl('<tr><td>a_m(j) [-]</td><td>: </td><td> activity of "j", defined on molality basis (used for inorg. ions) with reference &
        & state of infinite dilution of "j" in pure water;</td></tr>')
    write(unitx,'(A)') adjustl('<tr><td>log<sub>10</sub>(&eta;/[Pa.s])</td><td>: </td><td> base-10 log of the dynamic viscosity of the mixture;</td></tr>')
    write(unitx,'(A)') adjustl('<tr><td>&plusmn;log<sub>10</sub>(&eta;&nbsp;sens./[Pa.s])</td><td>: </td><td> base-10 log of the estimated sensitivity/uncertainty of &
        & the dynamic viscosity of the mixture; see details <a href="../help.html#viscosity">here</a>;</td></tr>')
    write(unitx,'(A)') adjustl('<tr><td>flag</td><td>: </td><td> error/warning flag, a non-zero value (error/warning number) indicates that a numerical issue &
        & or a warning occurred at this data point, suggesting evaluation with caution (warnings) or exclusion (errors) of this data point; see also &
        & the <a href="../AIOMFAC-model/List_of_Error_and_Warning_Codes.pdf">List_of_Error_and_Warning_Codes.pdf</a>.</td></tr>')
    write(unitx,'(A)') adjustl('</table>')
    write(unitx,'(A)') '<br>'
    write(unitx,'(A)') '<br>'

    !--
    !data table for mixture viscosity output
    write(unitx,'(A)') '<br>'
    !write component properties table above data table:
    write(unitx,'(A)') adjustl('<table class="datatabname">')
    write(unitx,'(A)') '<tbody>'
    write(unitx,'(A)') "<tr><td>Properties of this phase: mixture viscosity</td></tr>"
    write(unitx,'(A)') '</tbody>'
    write(unitx,'(A)') '</table>'
    !write table column headers:
    write(unitx,'(A)') adjustl('<table class="datatable">')
    write(unitx,'(A)') '<thead>'
    txtsubs = "<tr><th> no. </th><th>&nbsp;</th><th> T [K] </th><th>&nbsp;</th><th> RH [%] </th><th>&nbsp;</th><th> &
        & log<sub>10</sub>(&eta;/[Pa.s])</th><th>&nbsp;</th><th> &plusmn;log<sub>10</sub>(&eta;&nbsp;sens./[Pa.s]) </th><th>&nbsp;</th><th> flag </th></tr> "
    write(unitx,'(A)') adjustl(txtsubs) 
    write(unitx,'(A)') '</thead>'
    write(unitx,'(A)') '<tbody>'
    !--
    !write data to table:
    do pointi = 1,npoints  !loop over composition points
        if (watercompno > 0) then !water is present in mixture
            RH = out_data(5,pointi,watercompno)*100.0_wp !RH in %
            if (RH > 1000.0_wp .or. RH < 0.0_wp) then
                RH = -99.99_wp
            endif
        else
            RH = 0.0_wp
        endif
        tformat = '(A8,I5.3,"</td><td>&nbsp;</td><td>",F7.2,"</td><td>&nbsp;</td><td>",F7.2,"</td><td>&nbsp;</td><td>",2(ES12.5,"</td><td>&nbsp;</td><td>"),I2,A10)'
        write(unitx, tformat) "<tr><td>", pointi, T_K(pointi), RH, out_viscdata(1,pointi), out_viscdata(2,pointi), INT(out_viscdata(3,pointi)), "</td></tr>"
    enddo !pointi
    write(unitx,'(A)') '</tbody>'
    write(unitx,'(A)') '</table>'
    write(unitx,'(A)') '<a href="#top">&uarr; Top</a>'
    write(unitx,'(A)') "<br>"
    write(unitx,'(A)') "<br>"

    !write individual data tables for each component / ionic species:
    do i = 1,nspecmax
        write(cn,cnformat) i !component / species number as character string
        if (INT(out_data(6,px(i),i)) == 0) then !neutral component
            write(unitx,'(A)') '<br>'
            write(unitx,'(A)') adjustl('<a id="'//trim(adjustl(compname(i)))//'"></a>')
            !write component properties table above data table:
            write(unitx,'(A)') adjustl('<table class="datatabname">')
            write(unitx,'(A)') '<tbody>'
            write(unitx,'(A,I0.2,A)') "<tr><td> Mixture's component, #  </td> <td> &nbsp; : &nbsp; ", i ,"</td></tr>"
            txtn = trim(adjustl(compname(i)))
            txtn = trim(txtn)//"</td></tr>"
            write(unitx, '(A,A)') "<tr><td> Component's name </td> <td> &nbsp; : &nbsp; ", adjustl(txtn)
            txtsubs = '<span class="chemf">'//trim(compsubgroupsHTML(i))//"</span></td></tr>"
            write(unitx, '(A,A)') "<tr><td> Component's subgroups </td> <td> &nbsp; : &nbsp; ", adjustl(txtsubs)
            write(unitx,'(A)') '</tbody>'
            write(unitx,'(A)') '</table>'
            !write table column headers:
            write(unitx,'(A)') adjustl('<table class="datatable">')
            write(unitx,'(A)') '<thead>'
            if (i == idCO2) then  !check exception: CO2 has the activity defined on molality basis like ions:
                outtxtleft = adjustl("<tr><th> no. </th><th>&nbsp;</th><th> T [K] </th><th>&nbsp;</th><th> RH [%] </th><th>&nbsp;</th><th> &
                & w("//cn//")</th><th>&nbsp;</th><th> x_i("//cn//")</th><th>&nbsp;</th><th> m_i("//cn//")</th><th>&nbsp;</th><th> &
                & a_coeff_m("//cn//")</th><th>&nbsp;</th><th> a_m("//cn//")</th><th>&nbsp;</th><th> flag </th></tr> ")
            else        !regular case
                outtxtleft = adjustl("<tr><th> no. </th><th>&nbsp;</th><th> T [K] </th><th>&nbsp;</th><th> RH [%] </th><th>&nbsp;</th><th> &
                & w("//cn//")</th><th>&nbsp;</th><th> x_i("//cn//")</th><th>&nbsp;</th><th> m_i("//cn//")</th><th>&nbsp;</th><th> &
                & a_coeff_x("//cn//")</th><th>&nbsp;</th><th> a_x("//cn//")</th><th>&nbsp;</th><th> flag </th></tr> ")
            endif
            write(unitx,'(A)') outtxtleft
            !--
        else if (INT(out_data(6,px(i),i)) < 240) then !cation
            write(unitx,'(A)') '<br>'
            write(unitx,'(A)') '<a id="'//trim(adjustl(subgrnameHTML(INT(out_data(6,px(i),i)))))//'"></a>' !link target
            !---
            !write component properties table above data table:
            write(unitx,'(A)') adjustl('<table class="datatabname">')
            write(unitx,'(A)') '<tbody>'
            write(unitx,'(A,I0.2,A)') "<tr><td> Mixture's species, # </td> <td> &nbsp; : &nbsp; ", i ,"</td></tr>"
            subntxt = trim(adjustl(subgrnameHTML(INT(out_data(6,px(i),i)))))
            qty = len_trim(subntxt)
            txtn = adjustl(subntxt(2:qty-1)) !to print the ion name without the enclosing parenthesis ()
            txtn = '<span class="chemf">'//trim(txtn)//"</span></td></tr>"
            write(unitx, '(A,A)') "<tr><td> Cation's name </td> <td> &nbsp; : &nbsp; ", adjustl(txtn)
            txtsubs = '<span class="chemf">'//trim(subntxt)//"</span></td></tr>"
            write(unitx, '(A,A)') "<tr><td> Cation's subgroups </td> <td> &nbsp; : &nbsp; ", adjustl(txtsubs)
            write(unitx,'(A)') '</tbody>'
            write(unitx,'(A)') '</table>'
            !write data table column headers:
            write(unitx,'(A)') adjustl('<table class="datatable">')
            write(unitx,'(A)') '<thead>'
            outtxtleft = adjustl('<tr><th> no. </th><th>&nbsp;</th><th> T [K] </th><th>&nbsp;</th><th> RH [%] </th><th>&nbsp;</th><th> &
                & w('//cn//')</th><th>&nbsp;</th><th> x_i('//cn//')</th><th>&nbsp;</th><th> m_i('//cn//')</th><th>&nbsp;</th><th> &
                & a_coeff_m('//cn//')</th><th>&nbsp;</th><th> a_m('//cn//')</th><th>&nbsp;</th><th> flag </th></tr> ')
            write(unitx,'(A)') outtxtleft
            !--
        else if (INT(out_data(6,px(i),i)) > 240) then !anion
            write(unitx,'(A)') '<br>'
            write(unitx,'(A)') '<a id="'//trim(adjustl(subgrnameHTML(INT(out_data(6,px(i),i)))))//'"></a>'
            !---
            !write component properties table above data table:
            write(unitx,'(A)') adjustl('<table class="datatabname">')
            write(unitx,'(A)') '<tbody>'
            write(unitx,'(A,I0.2,A)') "<tr><td> Mixture's species, # </td> <td> &nbsp; : &nbsp; ", i ,"</td></tr>"
            subntxt = trim(adjustl(subgrnameHTML(INT(out_data(6,px(i),i)))))
            qty = len_trim(subntxt)
            txtn = adjustl(subntxt(2:qty-1)) !to print the ion name without the enclosing parenthesis ()
            txtn = '<span class="chemf">'//trim(txtn)//"</span></td></tr>"
            write(unitx, '(A,A)')  "<tr><td> Anion's name </td> <td> &nbsp; : &nbsp; ", adjustl(txtn)
            txtsubs = '<span class="chemf">'//trim(subntxt)//"</span></td></tr>"
            write(unitx, '(A,A)') "<tr><td> Anion's subgroups </td> <td> &nbsp; : &nbsp; ", adjustl(txtsubs)
            write(unitx,'(A)') '</tbody>'
            write(unitx,'(A)') '</table>'
            !write data table column headers:
            write(unitx,'(A)') adjustl('<table class="datatable">')
            write(unitx,'(A)') '<thead>'
            outtxtleft = adjustl("<tr><th> no. </th><th>&nbsp;</th><th> T [K] </th><th>&nbsp;</th><th> RH [%] </th><th>&nbsp;</th><th> &
                & w("//cn//")</th><th>&nbsp;</th><th> x_i("//cn//")</th><th>&nbsp;</th><th> m_i("//cn//")</th><th>&nbsp;</th><th> &
                & a_coeff_m("//cn//")</th><th>&nbsp;</th><th> a_m("//cn//")</th><th>&nbsp;</th><th> flag </th></tr> ")
            write(unitx,'(A)') outtxtleft
            !--
        else
            !error
        endif
        !--
        write(unitx,'(A)') '</thead>'
        write(unitx,'(A)') '<tbody>'
        !write data to table:
        do pointi = 1,npoints  !loop over composition points
            if (watercompno > 0) then !water is present in mixture
                RH = out_data(5,pointi,watercompno)*100.0_wp !RH in %
                if (RH > 1000.0_wp .or. RH < 0.0_wp) then
                    RH = -99.99_wp
                endif
            else
                RH = 0.0_wp
            endif
            tformat = '(A8,I5.3,"</td><td>&nbsp;</td><td>",F7.2,"</td><td>&nbsp;</td><td>",F7.2,"</td><td>&nbsp;</td><td>",5(ES12.5, &
                & "</td><td>&nbsp;</td><td>"),I2,A10)'
            write(unitx, tformat) "<tr><td>", pointi, T_K(pointi), RH, (out_data(k,pointi,i), k = 1,5), INT(out_data(7,pointi,i)), "</td></tr>"
        enddo !pointi
        write(unitx,'(A)') '</tbody>'
        write(unitx,'(A)') '</table>'
        write(unitx,'(A)') '<a href="#top">&uarr; Top</a>'
        write(unitx,'(A)') "<br>"
        write(unitx,'(A)') "<br>"
    enddo

    close(unitx)
        
    end subroutine Output_HTML
    !========================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Report error and warning messages to the error logfile from a list of defined      *
    !*   cases.                                                                             *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
    !*                                                                                      *
    !*   -> created:        2011                                                            *
    !*   -> latest changes: 2022-01-17                                                      *
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
    !****************************************************************************************
    subroutine RepErrorWarning(unito, errorflagmix, warningflag, errorflag_list, pointi, errorind, warningind)

    use ModSystemProp, only : errorflag_clist

    implicit none
    !interface variables:
    integer,intent(in) :: unito, errorflagmix, warningflag, pointi
    logical,dimension(size(errorflag_clist)),intent(in) :: errorflag_list
    integer,intent(out) :: errorind, warningind
    !local variables:
    integer :: en, n
    !...................................................................................

    if (errorflagmix /= 0) then !some mixture related error occurred:
        SELECT CASE(errorflagmix)
        CASE(1)
            write(unito,*) ""
            write(unito,*) "======================================================="
            write(unito,*) "AIOMFAC ERROR 1: mixture related."
            write(unito,*) "An organic main group <-> cation interaction parameter"
            write(unito,*) "is not defined for the requested mixture. "
            write(unito,*) "Please check your mixture components for available parameters"
            write(unito,*) "stated in the AIOMFAC interaction matrix."
            write(unito,*) "======================================================="
            write(unito,*) ""
        CASE(2)
            write(unito,*) ""
            write(unito,*) "======================================================="
            write(unito,*) "AIOMFAC ERROR 2: mixture related."
            write(unito,*) "An organic main group <-> anion interaction parameter"
            write(unito,*) "is not defined for the requested mixture. "
            write(unito,*) "Please check your mixture components for available parameters"
            write(unito,*) "stated in the AIOMFAC interaction matrix."
            write(unito,*) "======================================================="
            write(unito,*) ""
        CASE(9)
            write(unito,*) ""
            write(unito,*) "======================================================="
            write(unito,*) "AIOMFAC ERROR 9: mixture related."
            write(unito,*) "At least one cation <-> anion interaction parameter is "
            write(unito,*) "not defined for the requested mixture. "
            write(unito,*) "Please check all ion combinations for available parameters"
            write(unito,*) "stated in the AIOMFAC interaction matrix."
            write(unito,*) "======================================================="
            write(unito,*) ""
        CASE(13)
            write(unito,*) ""
            write(unito,*) "======================================================="
            write(unito,*) "AIOMFAC ERROR 13: Incorrect hydroxyl group assignment."
            write(unito,*) "At least one component containing (CH_n[(OH)]) groups"
            write(unito,*) "has been assigned an incorrect number of (OH) groups."
            write(unito,*) "Note that the notation of a CHn group bonded to an OH"
            write(unito,*) "group does not include the OH group; rather the hydroxyl"
            write(unito,*) "groups have to be defined separately.                   "
            write(unito,*) "======================================================="
            write(unito,*) ""
        CASE(14)
            write(unito,*) ""
            write(unito,*) "======================================================="
            write(unito,*) "AIOMFAC ERROR 14: Missing short-range ARR parameter."
            write(unito,*) "A neutral main group <-> main group interaction coeff."
            write(unito,*) "of this particular mixture is not available in the SR."
            write(unito,*) "part of the model."
            write(unito,*) "Check your organic components and their subgroups in"
            write(unito,*) "comparison to available subgroups in the AIOMFAC matrix."
            write(unito,*) "======================================================="
            write(unito,*) ""
        CASE(15)
            write(unito,*) ""
            write(unito,*) "======================================================="
            write(unito,*) "AIOMFAC ERROR 15: Missing short-range BRR parameter."
            write(unito,*) "A neutral main group <-> main group interaction coeff."
            write(unito,*) "of this particular mixture is not available in the SR."
            write(unito,*) "part of the model for 3-parameter temperature dependence."
            write(unito,*) "Check your organic components and their subgroups in"
            write(unito,*) "comparison to available subgroups in the AIOMFAC matrix."
            write(unito,*) "======================================================="
            write(unito,*) ""
        CASE DEFAULT
            write(unito,*) ""
            write(unito,*) "======================================================="
            write(unito,*) "AIOMFAC ERROR XX: an undefined mixture error occurred!"
            write(unito,*) "errorflagmix = ", errorflagmix
            write(unito,*) "======================================================="
            write(unito,*) ""
        end SELECT
        errorind = errorflagmix
    else 
        !check warnings and errors related to occurences during specific data point calculations:
        if (warningflag > 0) then
            SELECT CASE(warningflag)
            CASE(10)
                write(unito,*) ""
                write(unito,*) "======================================================="
                write(unito,*) "AIOMFAC WARNING 10: Temperature range related."
                write(unito,*) "At least one data point has a set temperature outside of"
                write(unito,*) "the recommended range for model calculations of "
                write(unito,*) "electrolyte-containing mixtures. This may be intended, "
                write(unito,*) "but caution is advised as AIOMFAC is not designed to  "
                write(unito,*) "perform well at this temperature."
                write(unito,*) "Data point no.: ", pointi
                write(unito,*) "======================================================="
                write(unito,*) ""
            CASE(11)
                write(unito,*) ""
                write(unito,*) "======================================================="
                write(unito,*) "AIOMFAC WARNING 11: Temperature range related."
                write(unito,*) "At least one data point has a set temperature outside of"
                write(unito,*) "the recommended range for model calculations of "
                write(unito,*) "electrolyte-free organic mixtures. This may be intended,"
                write(unito,*) "but caution is advised as AIOMFAC is not designed to "
                write(unito,*) "perform well at this temperature."
                write(unito,*) "Data point no.: ", pointi
                write(unito,*) "======================================================="
                write(unito,*) ""
            CASE(16)
                write(unito,*) ""
                write(unito,*) "======================================================="
                write(unito,*) "AIOMFAC-VISC WARNING 16: Mixture viscosity issue.      "
                write(unito,'(A)') "A problem occurred during the viscosity prediction, &
                    &likely related to a missing pure-component viscosity value. &
                    &Therefore, an unrealistic mixture viscosity of log_10(eta/[Pa.s]) = &
                    &-9999.9999 is output."
                write(unito,*) "Data point no.: ", pointi
                write(unito,*) "======================================================="
                write(unito,*) ""
            CASE DEFAULT
                write(unito,*) ""
                write(unito,*) "======================================================="
                write(unito,*) "AIOMFAC WARNING XX: an undefined WARNING occurred!"
                write(unito,*) "warningflag = ", warningflag
                write(unito,*) "======================================================="
                write(unito,*) ""
            end SELECT
            warningind = warningflag
        endif !warningflag
        if ( any(errorflag_list) ) then
            en = 0
            do n = 1,COUNT(errorflag_list)
                en = en + findloc(errorflag_list(en+1:), VALUE = .true., dim=1)
                SELECT CASE(en)
                CASE(3)
                    write(unito,*) ""
                    write(unito,*) "======================================================="
                    write(unito,*) "AIOMFAC ERROR 3: Mixture composition related."
                    write(unito,*) "Composition data for this point is missing or incorrect!"
                    write(unito,*) "The sum of the mole or mass fractions of all components"
                    write(unito,*) "has to be equal to 1.0 and individual mole or mass "
                    write(unito,*) "fractions have to be positive values <= 1.0!"
                    write(unito,*) "Composition point no.: ", pointi
                    write(unito,*) "======================================================="
                    write(unito,*) ""
                CASE(4,5)
                    write(unito,*) ""
                    write(unito,*) "======================================================="
                    write(unito,*) "AIOMFAC ERROR 4: Mixture composition related."
                    write(unito,*) "Composition data for this point is incorrect!"
                    write(unito,*) "The sum of the mole or mass fractions of all components"
                    write(unito,*) "has to be equal to 1.0 and individual mole or mass "
                    write(unito,*) "fractions have to be positive values <= 1.0!"
                    write(unito,*) "Composition point no.: ", pointi
                    write(unito,*) "======================================================="
                    write(unito,*) ""
                CASE(6,7)
                    write(unito,*) ""
                    write(unito,*) "======================================================="
                    write(unito,*) "AIOMFAC ERROR 6: Numerical issue."
                    write(unito,*) "A numerical issue occurred during computation of the"
                    write(unito,*) "data points flagged in the output tables."
                    write(unito,*) "This error was possibly caused due to input of very"
                    write(unito,*) "high electrolyte concentrations."
                    write(unito,*) "Composition point no.: ", pointi
                    write(unito,*) "======================================================="
                    write(unito,*) ""
                CASE(8)
                    write(unito,*) ""
                    write(unito,*) "======================================================="
                    write(unito,*) "AIOMFAC ERROR 8: Mixture composition related."
                    write(unito,*) "At least one neutral component must be present in the  "
                    write(unito,*) "system! "
                    write(unito,*) "Composition point no.: ", pointi
                    write(unito,*) "======================================================="
                    write(unito,*) ""
                CASE(12)
                    write(unito,*) ""
                    write(unito,*) "======================================================="
                    write(unito,*) "AIOMFAC ERROR 12: Charge neutrality violated."
                    write(unito,*) "The mixture violates the electrical charge neutrality  "
                    write(unito,*) "condition (moles cation*[cation charge] =              "
                    write(unito,*) "                          moles anion*[anion charge]). "
                    write(unito,*) "Make sure that selected integer amounts of cation and  "
                    write(unito,*) "anion 'subgroups' fulfill the charge balance (in the   "
                    write(unito,*) "inorganic component definition of the input file).     "
                    write(unito,*) "Composition point no.: ", pointi
                    write(unito,*) "======================================================="
                    write(unito,*) ""
                CASE(17)
                    write(unito,*) ""
                    write(unito,*) "======================================================="
                    write(unito,'(A)') "AIOMFAC ERROR 17: Issue with ion dissociation &
                                    &equilibria calculations."
                    write(unito,'(A)') "The numerical solution of electrolyte/ion dissociation &
                                    &equilibria was not accomplished to the desired tolerance &
                                    &level. Model output for this point is unreliable and likely &
                                    &incorrect."
                    write(unito,*) "Composition point no.: ", pointi
                    write(unito,*) "======================================================="
                    write(unito,*) ""
                CASE(18)
                    write(unito,*) ""
                    write(unito,*) "======================================================="
                    write(unito,'(A)') "AIOMFAC-VISC WARNING 18: Mixture viscosity issue."
                    write(unito,'(A)') "At least one of the cation–anion combinations included &
                        &in the mixture have not yet been supported for mixture viscosity &
                        &calculations. Therefore, the predicted viscosity is not reliable."
                    write(unito,*) "Composition point no.: ", pointi
                    write(unito,*) "======================================================="
                    write(unito,*) ""
                
                CASE DEFAULT
                    write(unito,*) ""
                    write(unito,*) "======================================================="
                    write(unito,*) "AIOMFAC ERROR XX: an undefined calculation error occurred!"
                    write(unito,*) "error flag = ", errorflag_list(en)
                    write(unito,*) "======================================================="
                    write(unito,*) ""
                end SELECT
            enddo
            errorind = findloc(errorflag_list(:), VALUE = .true., dim=1)
        endif
    endif !errorflagmix
        
    end subroutine RepErrorWarning
    !========================================================================================

end module Mod_InputOutput