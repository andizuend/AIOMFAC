!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Main program section to interface the AIOMFAC model with input parameters received * 
!*   from a webpage or command line.  The program gets the name of an input data file   *
!*   via a command line argument and calls AIOMFAC using these input parameters. Output *
!*   is produced in the form of an output txt-file and/or a related .html file for      *
!*   webpage display.                                                                   *
!*                                                                                      *
!*   Relative folder structure expected for input/output; directories 'Inputfiles' and  *
!*   'Outputfiles' need to be present at run time with read & write permissions:        *
!*   [location of executable AIOMFAC-web.out] -> [./Inputfiles/input_0???.txt]          *                       
!*   [location of executable AIOMFAC-web.out] -> [./Outputfiles/output_0???.txt]        *
!*                                                                                      *
!*   The AIOMFAC model expressions and parameters are described in Zuend et al. (2008,  * 
!*   Atmos. Chem. Phys.) and Zuend et al. (2011, Atmos. Chem. Phys.). Interaction       *
!*   parameters of Zuend et al. (2011) are used where they differ from the previous     *
!*   version. Additional parameters from Zuend and Seinfeld (2012), e.g. for peroxides, *
!*   are included as well.                                                              *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
!*                                                                                      *
!*   -> created:        2011  (this file)                                               *
!*   -> latest changes: 2018/08/10                                                      *
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
!****************************************************************************************
    
PROGRAM Main_IO_driver

!module variables:
!USE IFPORT  !for directory creation
USE ModSystemProp, ONLY : compname, compsubgroups, compsubgroupsHTML, compsubgroupsTeX, cpname, &
    & errorflagmix, NGS, nindcomp, ninput, NKNpNGS, nneutral, SetSystem, topsubno, waterpresent
USE ModSubgroupProp, ONLY : SubgroupNames, SubgroupAtoms, subgrname, subgrnameTeX, subgrnameHTML
USE ModMRpart, ONLY : MRdata
USE ModSRunifac, ONLY : SRdata

IMPLICIT NONE
!set preliminary input-related parameters:
INTEGER(4),PARAMETER :: maxpoints = 1000  !limit maximum number of composition points for web-version
INTEGER(4),PARAMETER :: ninpmax = 50      !set the maximum number of mixture components allowed (preliminary parameter)
!local variables:
CHARACTER(LEN=2) :: cn
CHARACTER(LEN=4) :: dashes, equalsigns, pluses, VersionNo
CHARACTER(LEN=20) :: dummy
CHARACTER(LEN=50) :: subntxt, txtcheck
CHARACTER(LEN=150) :: horizline, tablehead, tformat, txtn
CHARACTER(LEN=3000) :: errlogfile, filename, filepath, filepathout, fname, &
    &mixturestring, txtfilein, txtsubs  
CHARACTER(LEN=20),DIMENSION(ninpmax) :: txtarray
CHARACTER(LEN=60),DIMENSION(ninpmax) :: cpnameinp   !list of assigned component names (from input file)
CHARACTER(LEN=60),DIMENSION(:),ALLOCATABLE :: outnames
CHARACTER(LEN=400) :: outtxtleft
INTEGER(4) :: allocstat, cpno, errlogout, errorflagcalc, errorind, i, inpfilesize, &
    & istat, k, nc, ncp, npoints, nspecies, nspecmax, pointi, qty, standardout, &
    & subg, unito, unitx, warningflag, warningind, watercompno
INTEGER(4),DIMENSION(ninpmax) :: px
INTEGER(4),DIMENSION(ninpmax,topsubno) :: cpsubg   !list of input component subgroups and corresponding subgroup quantities
REAL(8) :: TKelvin, RH
REAL(8),DIMENSION(maxpoints) :: T_K
REAL(8),DIMENSION(:),ALLOCATABLE :: inputconc
REAL(8),DIMENSION(:,:),ALLOCATABLE :: composition, outputvars
REAL(8),DIMENSION(:,:,:),ALLOCATABLE :: out_data
LOGICAL(4) :: ext, fileexists, fileopened, filevalid, verbose, xinputtype
!...................................................................................
!
!==== INITIALIZATION section =======================================================
!
ncp = 0
nspecmax = 0
cpsubg = 0
standardout = 6     !standard output to terminal console window (Fortran standard unit *, 6, or 101)
errlogout = 71      !output unit for error log file
cpnameinp = "none"  !default component's name
T_K = 298.15D0      !default temperature in [K]
dashes = "----"
equalsigns = "===="
pluses = "++++"
VersionNo = "2.20"  !AIOMFAC-web version number (change here if minor or major changes require a version number change)
unito = errlogout
errorind = 0        !0 means no error found
warningind = 0      !0 means no warnings found
filevalid = .false.
fileexists = .false.
verbose = .true.    !if true, some debugging information will be printed to the unit "unito" (errorlog file)
!
!==== INPUT data section ===========================================================
!
!read command line for text-file name (which contains the input parameters to run the AIOMFAC progam):
CALL GET_COMMAND_ARGUMENT(1, txtfilein)
!---
!txtfilein = './Inputfiles/input_0985.txt' !just use this for debugging, otherwise comment out
!---
filepath = ADJUSTL(TRIM(txtfilein))
WRITE(*,*) ""
WRITE(*,*) "MESSAGE from AIOMFAC-web: program started, command line argument 1 = ", filepath
WRITE(*,*) ""
!extract filepath from the input filepath (which could include a path to the directory):
k = LEN_TRIM(filepath)
IF (k < 1) THEN !no file was provided
    WRITE(*,*) ""
    WRITE(*,*) "ERROR from AIOMFAC-web: no input file was provided!"
    WRITE(*,*) "An input file path needs to be provided via command line argument 2"
    WRITE(*,*) ""
    READ(*,*)
    STOP
ENDIF
i = INDEX(filepath, "Inputfiles/input")
IF (i > 0 .AND. i < k) THEN !found a direcory path in UNIX file path style
    filename = filepath(i+11:k)
    filepath = filepath(1:i+10)
ELSE
    i = INDEX(filepath, "Inputfiles\input")
    IF (i > 0 .AND. i < k) THEN !found a direcory path in WINDOWS file path style
        filename = filepath(i+11:k)
        filepath = filepath(1:i+10)
    ELSE
        filename = filepath(i+11:k)
        filepath = ""
    ENDIF
ENDIF
!-----
!check / create for associated output directory "Outputfiles":
i = LEN_TRIM(filepath)
filepathout = TRIM(filepath(1:i-11))//"Outputfiles/"
!use the name of the input file to create a corresponding output file name:
i = INDEX(filename, ".txt") !returns starting position of string input within string filename
!create an error-logfile associated with the input file name:
errlogfile = "Errorlog_"//filename(i-4:)
fname = TRIM(filepathout)//TRIM(errlogfile)
OPEN (NEWUNIT = unito, FILE = fname, STATUS ='UNKNOWN') !unito is the error / logfile unit
!-----
!check if file exists and read its content if true:
fname = TRIM(filepath)//TRIM(filename)
INQUIRE(FILE = fname, EXIST = fileexists, SIZE = inpfilesize) !inpfilesize is the file size in [bytes]
!Delete very large files that can only mean uploaded spam content and not actual input:
IF (fileexists) THEN
    IF (FLOAT(inpfilesize) > 2E6) THEN !cannot be a valid input file
        fname = TRIM(filepath)//TRIM(filename)
        OPEN (NEWUNIT = unitx, FILE = fname, STATUS='OLD')
        CLOSE(unitx, STATUS='DELETE')  !close and delete the file
        fileexists = .false.
    ENDIF
ENDIF

!attempt to read a supposedly valid input file:
IF (fileexists) THEN
    IF (verbose) THEN  
        WRITE(unito,*) ""
        WRITE(unito,*) "MESSAGE from AIOMFAC: input file exists."
    ENDIF
    fname = TRIM(filepath)//TRIM(filename)
    OPEN (NEWUNIT = unitx, FILE = fname, IOSTAT=istat, ACTION='READ', STATUS='OLD')
    IF (istat /= 0) THEN ! an error occurred
        WRITE(unito,*) ""
        WRITE(unito,*) "MESSAGE from AIOMFAC: an error occurred while trying to open the file! IOSTAT: ", istat
    ENDIF
    !validate file as a correct input text-file (no spam):
    READ(unitx,*) dummy, dummy, dummy, txtcheck
    IF (txtcheck(1:11) == "AIOMFAC-web") THEN
        filevalid = .true.
        IF (verbose) THEN
            WRITE(unito,*) ""
            WRITE(unito,*) "MESSAGE from AIOMFAC: input file has passed the first line text validation and will be read."
        ENDIF
        !read input data from file:
        BACKSPACE unitx !jump back to beginning of record (to the beginning of the line)
        READ(unitx,*) dummy, dummy, dummy, dummy, dummy !read first line
        READ(unitx,*)               !read empty line 2
        READ(unitx,*) dummy, dummy  !read line 3
        READ(unitx,*) dummy         !read line 4
    ELSE !file not valid. It will be closed and deleted below
        filevalid = .false.
        IF (verbose) THEN
            WRITE(unito,*) ""
            WRITE(unito,*) "MESSAGE from AIOMFAC: input file does not pass first line text validation and will be deleted."
        ENDIF
    ENDIF
    !loop over mixture components with variable numbers of subgroups to read (using inner loop):
    IF (filevalid) THEN
        IF (verbose) THEN
            WRITE(unito,*) ""
            WRITE(unito,*) "MESSAGE from AIOMFAC: reading component data from input file."
        ENDIF
        DO ncp = 1,ninpmax
            IF (verbose .AND. istat /= 0) THEN
                WRITE(unito,*) ""
                WRITE(unito,*) "MESSAGE from AIOMFAC: file end found in input file at ncp = ", ncp
            ENDIF
            READ(unitx,*,IOSTAT=istat) txtcheck !read only first argument on this line for subsequent check
            IF (txtcheck(1:4) == pluses .OR. istat /= 0) THEN !"++++" indicates end of this components definition part
                EXIT !exit ncp DO-loop
            ELSE !in this case, argument 3 of txtcheck is the component no.:
                BACKSPACE unitx !jump back to beginning of record (to the beginning of the line)
                READ(unitx,*) dummy, dummy, cpno !read line with component's number
            ENDIF
            READ(unitx,*) dummy, dummy, cpnameinp(cpno) !read line with component's name
            IF (LEN_TRIM(cpnameinp(cpno)) < 1) THEN !no component name was assigned, generate a default name
                WRITE(cn,'(I2.2)') cpno
                cpnameinp(cpno) = "cp_"//cn
            ENDIF
            IF (verbose) THEN
                WRITE(unito,*) ""
                WRITE(unito,*) "MESSAGE from AIOMFAC: found a component cpno =", cpno
            ENDIF
            DO !until exit
                READ(unitx,*,IOSTAT=istat) txtcheck !read only first argument on this line
                !check whether another subgroup is present or not
                IF (txtcheck(1:4) == dashes .OR. istat /= 0) THEN !"----" indicates no more subgroups of this component
                    EXIT !leave the inner DO-loop
                ELSE !else continue reading the next subgroup
                    BACKSPACE unitx
                    READ(unitx,*) dummy, dummy, dummy, subg, qty  !continue reading line with subgroup no. and corresp. quantity
                    cpsubg(cpno,subg) = cpsubg(cpno,subg)+qty
                ENDIF
            ENDDO
        ENDDO
        ncp = ncp-1 !ncp-1 is the number of different components in this mixture
        IF (verbose) THEN
            WRITE(unito,*) ""
            WRITE(unito,*) "MESSAGE from AIOMFAC: component data read."
            WRITE(unito,*) "MESSAGE from AIOMFAC: number of components: ", ncp
        ENDIF
        READ(unitx,*) dummy, dummy     !read mixture composition line
        READ(unitx,*) dummy, dummy, i  !read mass fraction? line
        IF (i == 1) THEN !composition in mass fractions
            xinputtype = .false.
        ELSE !composition in mole fractions
            xinputtype = .true.
        ENDIF
        READ(unitx,*) dummy, dummy, i !read mole fraction? line
        READ(unitx,*) dummy !read dashes line
        !now read the lines with the mixture composition points:
        !first allocate the composition array with a max of maxpoints points and the number of components present:
        ALLOCATE(composition(maxpoints,ncp), STAT=allocstat)
        composition = 0.0D0 !initialize array
        READ(unitx,*) txtarray(1:ncp) !read header line of composition table
        DO npoints = 1,maxpoints !or until exit
            READ(unitx,*,IOSTAT=istat) txtcheck !read only the first argument on this line
            IF (txtcheck(1:4) == equalsigns .OR. istat /= 0) THEN !"====" indicates no more composition points (and last line of input file)
                EXIT !leave the DO-loop (normal exit point)
            ELSE IF (IACHAR(txtcheck(1:1)) > 47 .AND. IACHAR(txtcheck(1:1)) < 58) THEN !validate that the data is actual intended input and not some sort of text field spam).
                BACKSPACE unitx
                READ(unitx,*) txtcheck, dummy !read only the first two arguments on this line
                IF (IACHAR(dummy(1:1)) > 47 .AND. IACHAR(dummy(1:1)) < 58) THEN !it is a number
                    BACKSPACE unitx
                    READ(unitx,*) i, T_K(npoints), composition(npoints,2:ncp) !read the temperature in [K] and composition values of the components [2:ncp] into the array
                ELSE
                    IF (npoints == 1) THEN
                        filevalid = .false. !file is not valid due to incorrect input in text field
                        EXIT
                    ELSE
                        warningind = 31
                        EXIT !abnormal exit point of loop due to non-number entries at a certain line in the text field
                    ENDIF
                ENDIF
            ELSE
                IF (npoints == 1) THEN
                    filevalid = .false. !file is not valid due to incorrect input in text field
                    EXIT
                ELSE
                    warningind = 31
                    EXIT !abnormal exit point of loop due to non-number entries at a certain line in the text field
                ENDIF
            ENDIF
        ENDDO
        npoints = npoints-1 !the number of composition points
        IF (verbose) THEN
            WRITE(unito,*) ""
            WRITE(unito,*) "MESSAGE from AIOMFAC: composition points read."
            WRITE(unito,*) "MESSAGE from AIOMFAC: number of points: ",npoints
        ENDIF
        IF (.NOT. filevalid) THEN
            !close and delete file from server:
            CLOSE(unitx, STATUS = 'DELETE')
        ELSE
            CLOSE(unitx)
        ENDIF
        !assign component 01 its composition values based on the rest of the component's data:
        DO i = 1,npoints
            composition(i,1) = 1.0D0-SUM(composition(i,2:ncp))
        ENDDO
    ENDIF !filevalid
ELSE
    WRITE(unito,*) ""
    WRITE(unito,*) "ERROR in AIOMFAC: Input file does not exist at expected location: ", TRIM(filename)
    WRITE(unito,*) ""
ENDIF !fileexists

IF (fileexists) THEN
    IF (.NOT. filevalid) THEN !input file contains errors or is completely invalid (submitted spam etc.)
        IF (verbose) THEN
            WRITE(unito,*) ""
            WRITE(unito,*) "MESSAGE from AIOMFAC: input file did not pass full validation and may be a spam file."
            WRITE(unito,*) "MESSAGE from AIOMFAC: the input file will be deleted to prevent spam files and malicious code from occupying the server."
            WRITE(unito,*) ""
        ENDIF
        errorind = 32
        !close and delete file from server:
        INQUIRE(FILE = fname, OPENED = fileopened)
        IF (fileopened) THEN
            CLOSE(unitx, STATUS = 'DELETE')
        ENDIF
    ELSE
        !check for any valid points for calculations before initializing AIOMFAC:
        IF (npoints < 1) THEN !no points input for calculations
            errorind = 33 !no input in text field..
            IF (verbose) THEN
                WRITE(unito,*) ""
                WRITE(unito,*) "MESSAGE from AIOMFAC: no composition points have been entered in the text field. There is nothing to calculate."
                WRITE(unito,*) ""
            ENDIF
            filevalid = .false.
        ENDIF
    ENDIF
ENDIF

IF (fileexists .AND. filevalid) THEN
    !
    !==== AIOMFAC initialization and calculation section ===============================
    !
    IF (verbose) THEN
        WRITE(unito,*) ""
        WRITE(unito,*) "MESSAGE from AIOMFAC: input file read, starting AIOMFAC mixture definitions and initialization... "
        WRITE(unito,*) ""
    ENDIF
    !load the MR and SR interaction parameter data:
    CALL MRdata()         !initialize the MRdata for the interaction coeff. calculations
    CALL SRdata()         !initialize data for the SR part coefficient calculations
    CALL SubgroupNames()  !initialize the subgroup names for the construction of component subgroup strings
    CALL SubgroupAtoms()
    !--
    !set mixture properties based on the data from the input file:
    CALL SetSystem(1, .true., ncp, cpnameinp(1:ncp), cpsubg(1:ncp,1:topsubno) )

    !check whether water is present in the mixture and as which component number:
    watercompno = 0
    IF (waterpresent) THEN
        watercompno = MAXLOC(cpsubg(1:ncp,16), DIM=1) !usually = 1
    ENDIF
    
    IF (errorflagmix /= 0) THEN !some mixture related error occured:
        CALL RepErrorWarning(unito, errorflagmix, warningflag, errorflagcalc, i, errorind, warningind)
    ENDIF

    IF (errorind == 0) THEN !go ahead and perform calculations, else jump to termination section
        !--
        ALLOCATE(inputconc(nindcomp), outputvars(6,NKNpNGS), outnames(NKNpNGS), out_data(7,npoints,NKNpNGS), STAT=allocstat)
        inputconc = 0.0D0
        out_data = 0.0D0
        !--
        IF (verbose) THEN
            WRITE(unito,*) ""
            WRITE(unito,*) "MESSAGE from AIOMFAC: mixture defined, calculating composition points... "
            WRITE(unito,*) ""
        ENDIF
        px = 0
        !set AIOMFAC input and call the main AIOMFAC subroutines for all composition points:
        DO pointi = 1,npoints !loop over points, changing composition / temperature
            inputconc(1:ncp) = composition(pointi,1:ncp)
            TKelvin = T_K(pointi)
            !--
            CALL AIOMFAC_inout(inputconc, xinputtype, TKelvin, nspecies, outputvars, outnames, errorflagcalc, warningflag)
            !--
            IF (warningflag > 0 .OR. errorflagcalc > 0) THEN
                !$OMP CRITICAL errwriting
                CALL RepErrorWarning(unito, errorflagmix, warningflag, errorflagcalc, pointi, errorind, warningind)
                !$OMP END CRITICAL errwriting
            ENDIF
            nspecmax = MAX(nspecmax, nspecies) !figure out the maximum number of different species in mixture (accounting for the 
                !possibility of HSO4- dissoc. and different species at different data points due to zero mole fractions).
            DO nc = 1,nspecmax !loop over species (ions dissociated and treated as individual species):
                out_data(1:6,pointi,nc) = outputvars(1:6,nc) !out_data general structure: | data columns 1:7 | data point | component no.|
                out_data(7,pointi,nc) = REAL(errorflagcalc, KIND=8)
                IF (errorflagcalc == 0 .AND. warningflag > 0) THEN !do not overwrite an errorflag if present!
                    out_data(7,pointi,nc) = REAL(warningflag, KIND=8) 
                ENDIF
                IF (px(nc) == 0 .AND. out_data(6,pointi,nc) >= 0.0D0) THEN
                    px(nc) = pointi !use a point for which this component's abundance is given, i.e. mole fraction(nc) > 0.0!
                ENDIF
            ENDDO !nc
        ENDDO !pointi
        !
        !==== OUTPUT data-to-file section ==================================================
        !
        !use the name of the input file to create a corresponding output file name from a string like "inputfile_0004.txt"
        i = INDEX(filename, ".txt")
        !replace "input" by "output":
        filename = "AIOMFAC_output_"//filename(i-4:)
        !-- for debugging
        WRITE(unito,'(A80)') "................................................................................"
        WRITE(unito,'(A60)') "MESSAGE from AIOMFAC: computations successfully performed."
        WRITE(dummy,'(I0)') LEN_TRIM(filepathout)
        tformat = '(A13, A24, A18, A'//TRIM(dummy)//')'  !dynamic format specifier
        WRITE(unito, tformat) "Output file, ", TRIM(filename), " created at path: ", TRIM(filepathout)
        WRITE(unito,'(A80)') "................................................................................"
        !create an output ASCII text file with an overall mixture header and individual tables for all components / species (in case of ions)
        fname = TRIM(filepathout)//TRIM(filename)
        CALL OutputTXT(fname, VersionNo, cpnameinp(1:nspecmax), nspecmax, npoints, watercompno, T_K(1:npoints), px(1:nspecmax), out_data)
        !--
        !>> write output HTML-file
        i = LEN_TRIM(filename)
        filename = filename(1:i-3)//"html"
        fname = TRIM(filepathout)//TRIM(filename)
        CALL OutputHTML(fname, VersionNo, cpnameinp(1:nspecmax), nspecmax, npoints, watercompno, T_K(1:npoints), px(1:nspecmax), out_data)
        !
        !==== TERMINATION section ==========================================================
        !
        DEALLOCATE(inputconc, outputvars, outnames, out_data, STAT=allocstat)
        ext = ALLOCATED(composition)
        IF (ext) THEN
            DEALLOCATE(composition, STAT=allocstat)
        ENDIF
    ENDIF !errorind
ENDIF !fileexists and valid

WRITE(unito,*) "+-+-+-+-+"
WRITE(unito,*) "Final warning indicator (an entry '00' means no warnings found):"
WRITE(unito,'(I2.2)') warningind
WRITE(unito,*) "+-+-+-+-+"
WRITE(unito,*) ""
WRITE(unito,*) "########"
WRITE(unito,*) "Final error indicator (an entry '00' means no errors found):"
WRITE(unito,'(I2.2)') errorind
WRITE(unito,*) "########"
CLOSE(unito) !close the error log-file

WRITE(*,*) ""
WRITE(*,*) "MESSAGE from AIOMFAC: End of program; final error indicator: ", errorind
WRITE(*,*) ""
!READ(*,*)  !Pause; just for debugging and testing.
!
!==== THE END ======================================================================
!
END PROGRAM Main_IO_driver