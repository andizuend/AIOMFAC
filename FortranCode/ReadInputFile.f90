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
!*   -> latest changes: 2019/10/17                                                      *
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
SUBROUTINE ReadInputFile(filepath, folderpathout, filename, ninpmax, maxpoints, unito, verbose, ncp, npoints, &
    & warningind, errorind, filevalid, cpnameinp, cpsubg, T_K, composition, xinputtype)

USE ModSystemProp, ONLY : topsubno

IMPLICIT NONE

!interface variables:
CHARACTER(LEN=3000),INTENT(INOUT) :: filepath, folderpathout
CHARACTER(LEN=200),INTENT(INOUT) :: filename
INTEGER(4),INTENT(IN) :: ninpmax, maxpoints
INTEGER(4),INTENT(INOUT) :: unito
LOGICAL(4),INTENT(IN) :: verbose
INTEGER(4),INTENT(OUT) :: ncp, npoints
INTEGER(4),INTENT(INOUT) :: warningind, errorind
LOGICAL(4),INTENT(OUT) :: filevalid
CHARACTER(LEN=200),DIMENSION(:),INTENT(OUT) :: cpnameinp    !list of assigned component names (from input file)
INTEGER(4),DIMENSION(:,:),INTENT(OUT) :: cpsubg             !list of input component subgroups and corresponding subgroup quantities
REAL(8),DIMENSION(:),INTENT(OUT) :: T_K                     !temperature of data points in Kelvin
REAL(8),DIMENSION(:,:),INTENT(OUT) :: composition           !array of mixture composition points for which calculations should be run
LOGICAL(4),INTENT(OUT) :: xinputtype
!--
!local variables:
CHARACTER(LEN=:),ALLOCATABLE :: cn   !this assumes a maximum four-digit component number in the system (max. 9999); to be adjusted otherwise.
CHARACTER(LEN=4) :: dashes, equalsigns, pluses
CHARACTER(LEN=20) :: dummy, cnformat
CHARACTER(LEN=50) :: txtcheck, inpfolder, outpfolder
CHARACTER(LEN=3000) :: errlogfile, fname
CHARACTER(LEN=20),DIMENSION(ninpmax) :: txtarray
INTEGER(4) :: cpno, i, inpfilesize, istat, k, kinpf, qty, subg, unitx
LOGICAL(4) :: fileopened, fileexists
!...................................................................................

!initialize variables
ncp = 0
cpsubg = 0
cpnameinp = "none"
composition = 0.0D0
T_K = 298.15D0      !default temperature in [K]
dashes = "----"
equalsigns = "===="
pluses = "++++"
fileexists = .false.
filevalid = .false.

!==== READ INPUT data ===========================================================

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

i = INDEX(filepath, "/input_", BACK = .true.)  !find index position of input file, assuming that all input files start with "input_"
IF (i < k .AND. i > 0) THEN !found a direcory path in UNIX file path style
    filename = TRIM(filepath(i+1:k))
    filepath = TRIM(filepath(1:k))
    kinpf = INDEX(filepath(1:i-1), "/", BACK = .true.)  !find folder length of folder containing the input file
ELSE
    kinpf = 0
    i = INDEX(filepath, "\input_", BACK = .true.) 
    IF (i > 0 .AND. i < k) THEN !found a direcory path in WINDOWS file path style
        filename = TRIM(filepath(i+1:k))
        filepath = TRIM(filepath(1:k))
        kinpf = INDEX(filepath(1:i-1), "\", BACK = .true.)  !find folder length of folder containing the input file
    ENDIF
ENDIF
IF (kinpf > i-1) THEN
    kinpf = 0
ENDIF
inpfolder = TRIM(filepath(kinpf+1:i-1))
kinpf = LEN_TRIM(inpfolder)

!check / create for associated output directory:
outpfolder = "Outputfiles/"  !initialize
istat = INDEX(inpfolder, "inp")
IF (istat > 0 .AND. istat < kinpf) THEN
    outpfolder = "outp"//TRIM(inpfolder(4:))//"/"  
ELSE
    istat = INDEX(inpfolder, "Inp")
    IF (istat > 0) THEN
        outpfolder = "Outp"//TRIM(inpfolder(4:))//"/"   
    ENDIF
ENDIF
folderpathout = TRIM(filepath(1:i-(kinpf+1)))//TRIM(outpfolder)

!use the name of the input file to create a corresponding output file name:
i = INDEX(filename, ".txt") !returns starting position of string input within string filename
!create an error-logfile associated with the input file name:
errlogfile = "Errorlog_"//filename(i-4:)
fname = TRIM(folderpathout)//TRIM(errlogfile)
OPEN (NEWUNIT = unito, FILE = fname, STATUS ='UNKNOWN') !unito is the error / logfile unit
!-----
!check if file exists and read its content if true:
fname = TRIM(filepath)
INQUIRE(FILE = fname, EXIST = fileexists, SIZE = inpfilesize) !inpfilesize is the file size in [bytes]
!Delete very large files that can only mean uploaded spam content and not actual input:
IF (fileexists) THEN
    IF (FLOAT(inpfilesize) > 50*(ninpmax +ninpmax*maxpoints) ) THEN !likely not a valid input file
        fname = TRIM(filepath)
        OPEN (NEWUNIT = unitx, FILE = fname, IOSTAT=istat, ACTION='READ', STATUS='OLD')
        READ(unitx,*) dummy, dummy, dummy, txtcheck
        CLOSE(unitx)
        IF (.NOT. (txtcheck(1:11) == "AIOMFAC-web")) THEN !invalid file (likely spam)
            OPEN (NEWUNIT = unitx, FILE = fname, STATUS='OLD')
            CLOSE(unitx, STATUS='DELETE')  !close and delete the file
            fileexists = .false.
        ENDIF
    ENDIF
ENDIF

!attempt to read a supposedly valid input file:
IF (fileexists) THEN
    IF (verbose) THEN
        WRITE(unito,*) ""
        WRITE(unito,*) "MESSAGE from AIOMFAC: input file exists."
    ENDIF
    fname = TRIM(filepath)
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
        k = MAX(2, CEILING(LOG10(REAL(ninpmax))) ) !determine order of magnitude digits
        WRITE(dummy,'(I0)') k
        cnformat = "(I"//TRIM(dummy)//"."//TRIM(dummy)//")" !variable format specifier
        ALLOCATE(CHARACTER(LEN=k) :: cn)
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
                WRITE(cn,cnformat) cpno
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
        IF (ncp == ninpmax) THEN
            READ(unitx,*,IOSTAT=istat) txtcheck !read only first argument on this line for subsequent check
            IF (txtcheck(1:4) /= pluses) THEN
                errorind = 34
                filevalid = .false.
                WRITE(*,*) "AIOMFAC ERROR 34: maximum number of input components reached while reading input file."
                WRITE(*,*) "AIOMFAC ERROR 34: check whether ninpmax value is too small."
                WRITE(unito,*) ""
                WRITE(unito,*) "AIOMFAC ERROR 34: maximum number of input components reached while reading input file."
            ENDIF
        ENDIF
        IF (filevalid) THEN
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
                WRITE(unito,*) "MESSAGE from AIOMFAC: number of points: ", npoints
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
        ENDIF !filevalid2
    ENDIF !filevalid1
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
        INQUIRE(FILE = fname, OPENED = fileopened)
        IF (errorind == 0) THEN
            errorind = 32
            !close and delete file from server:
            IF (fileopened) THEN
                CLOSE(unitx, STATUS = 'DELETE')
            ENDIF
        ELSE
            IF (fileopened) THEN
                CLOSE(unitx)
            ENDIF
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

END SUBROUTINE ReadInputFile