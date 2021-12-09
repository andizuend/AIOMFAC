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
!*   :: Methods and model references ::                                                 *
!*   The AIOMFAC model expressions and parameters are described in Zuend et al. (2008,  * 
!*   Atmos. Chem. Phys.) and Zuend et al. (2011, Atmos. Chem. Phys.). Interaction       *
!*   parameters of Zuend et al. (2011) are used where they differ from the previous     *
!*   version. Additional parameters from Zuend and Seinfeld (2012), e.g. for peroxides, *
!*   are included as well. Several ions and their interactions with other species,      *
!*   including I-, IO3-, CO3--, HCO3- and CO2(aq), are included based on Yin et al.     *
!*   (2021, Atmos. Chem. Phys.).                                                        *
!*   Viscosity predictions via AIOMFAC-VISC are included based on the articles by       *
!*   Gervasi et al. (2020) and Lilek and Zuend (2021).                                  *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
!*                                                                                      *
!*   -> created:        2011  (this file)                                               *
!*   -> latest changes: 2021-12-09                                                      *
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
USE ModSystemProp, ONLY : errorflagmix, idCO2, nindcomp, NKNpNGS, SetSystem, topsubno, waterpresent
USE ModSubgroupProp, ONLY : SubgroupAtoms, SubgroupNames
USE ModMRpart, ONLY : MRdata
USE ModSRparam, ONLY : SRdata

IMPLICIT NONE
!set preliminary input-related parameters:
INTEGER(4),PARAMETER :: maxpoints = 1001    !limit maximum number of composition points for web-version
INTEGER(4),PARAMETER :: ninpmax = 51        !set the maximum number of mixture components allowed (preliminary parameter)
!local variables:
CHARACTER(LEN=4) :: VersionNo
CHARACTER(LEN=200) :: filename
CHARACTER(LEN=3000) :: filepath, folderpathout, fname, txtfilein  
CHARACTER(LEN=200),DIMENSION(:),ALLOCATABLE :: cpnameinp   !list of assigned component names (from input file)
CHARACTER(LEN=200),DIMENSION(:),ALLOCATABLE :: outnames
INTEGER(4) :: allocstat, errorflagcalc, errorind, i, nc, ncp, npoints, nspecies, nspecmax, pointi, &
    & unito, warningflag, warningind, watercompno
INTEGER(4),DIMENSION(ninpmax) :: px
INTEGER(4),DIMENSION(:,:),ALLOCATABLE :: cpsubg   !list of input component subgroups and corresponding subgroup quantities
REAL(8) :: TKelvin
REAL(8),DIMENSION(:),ALLOCATABLE :: T_K
REAL(8),DIMENSION(:),ALLOCATABLE :: inputconc, outputviscvars
REAL(8),DIMENSION(:,:),ALLOCATABLE :: composition, compos2, outputvars, out_viscdata
REAL(8),DIMENSION(:,:,:),ALLOCATABLE :: out_data
LOGICAL(4) :: filevalid, verbose, xinputtype
!--
!explicit interfaces:
INTERFACE
    SUBROUTINE ReadInputFile(filepath, folderpathout, filename, ninpmax, maxpoints, unito, verbose, ncp, npoints, &
        & warningind, errorind, filevalid, cpnameinp, cpsubg, T_K, composition, xinputtype)
        USE ModSystemProp, ONLY : topsubno
        CHARACTER(LEN=3000),INTENT(INOUT) :: filepath, folderpathout
        CHARACTER(LEN=200),INTENT(INOUT) :: filename
        INTEGER(4),INTENT(IN) :: ninpmax, maxpoints
        INTEGER(4),INTENT(INOUT) :: unito
        LOGICAL(4),INTENT(IN) :: verbose
        INTEGER(4),INTENT(OUT) :: ncp, npoints
        INTEGER(4),INTENT(INOUT) :: warningind, errorind
        LOGICAL(4),INTENT(OUT) :: filevalid
        CHARACTER(LEN=200),DIMENSION(:),INTENT(OUT) :: cpnameinp !list of assigned component names (from input file)
        INTEGER(4),DIMENSION(:,:),INTENT(OUT) :: cpsubg          !list of input component subgroups and corresponding subgroup quantities
        REAL(8),DIMENSION(:),INTENT(OUT) :: T_K                  !temperature of data points in Kelvin
        REAL(8),DIMENSION(:,:),INTENT(OUT) :: composition        !array of mixture composition points for which calculations should be run
        LOGICAL(4),INTENT(OUT) :: xinputtype
    END SUBROUTINE ReadInputFile
END INTERFACE
!...................................................................................

!
!==== INITIALIZATION section =======================================================
!
VersionNo = "3.01"      !AIOMFAC-web version number (change here if minor or major changes require a version number change)
verbose = .true.        !if true, some debugging information will be printed to the unit "unito" (errorlog file)
nspecmax = 0
errorind = 0            !0 means no error found
warningind = 0          !0 means no warnings found
!
!==== INPUT data section ===========================================================
!
!read command line for text-file name (which contains the input parameters to run the AIOMFAC progam):
CALL GET_COMMAND_ARGUMENT(1, txtfilein)
IF (LEN_TRIM(txtfilein) < 4) THEN               !no command line argument; use specific input file for tests;
    txtfilein = './Inputfiles/input_0927.txt'   !just use this for debugging with a specific input file, otherwise comment out;
ENDIF
filepath = ADJUSTL(TRIM(txtfilein))
WRITE(*,*) ""
WRITE(*,'(A,A)') "MESSAGE from AIOMFAC-web: program started, command line argument 1 = ", TRIM(filepath)
WRITE(*,*) ""
ALLOCATE(cpsubg(ninpmax,topsubno), cpnameinp(ninpmax), composition(maxpoints,ninpmax), T_K(maxpoints), STAT=allocstat)
!--
CALL ReadInputFile(filepath, folderpathout, filename, ninpmax, maxpoints, unito, verbose, ncp, npoints, &
    & warningind, errorind, filevalid, cpnameinp, cpsubg, T_K, composition, xinputtype)
!--
IF (filevalid) THEN
    !
    !==== AIOMFAC initialization and calculation section ===============================
    !
    IF (verbose) THEN
        WRITE(unito,*) ""
        WRITE(unito,'(A)') "MESSAGE from AIOMFAC: input file read, starting AIOMFAC mixture definitions and initialization... "
        WRITE(unito,*) ""
    ENDIF
    !load the MR and SR interaction parameter data:
    CALL MRdata()         !initialize the MR data for the interaction coeff. calculations
    CALL SRdata()         !initialize data for the SR part coefficient calculations
    CALL SubgroupNames()  !initialize the subgroup names for the construction of component subgroup strings
    CALL SubgroupAtoms()
    !--
    !set mixture system properties (composition-independent ones) based on the data from the input file:
    CALL SetSystem(1, .true., ncp, cpnameinp(1:ncp), cpsubg(1:ncp,1:topsubno) )

    !check whether water is present in the mixture and as which component number:
    watercompno = 0
    IF (waterpresent) THEN
        watercompno = MAXLOC(cpsubg(1:ncp,16), DIM=1)   !usually = 1
    ENDIF
    IF (idCO2 > 0) THEN             !add CO2 as a (non-input) neutral component name
        cpnameinp(idCO2+1:ncp+1) = cpnameinp(idCO2:ncp)
        cpnameinp(idCO2) = 'CO2(aq)'   
    ENDIF
    !transfer composition data to adequate array size:
    ALLOCATE(compos2(npoints,ncp), STAT=allocstat)
    DO nc = 1,ncp
        compos2(1:npoints,nc) = composition(1:npoints,nc)
    ENDDO
    DEALLOCATE(cpsubg, composition, STAT=allocstat)
    
    IF (errorflagmix /= 0) THEN     !a mixture-related error occured:
        CALL RepErrorWarning(unito, errorflagmix, warningflag, errorflagcalc, i, errorind, warningind)
    ENDIF

    IF (errorind == 0) THEN         !perform AIOMFAC calculations; else jump to termination section
        !--
        ALLOCATE(inputconc(nindcomp), outputvars(6,NKNpNGS), outputviscvars(2), outnames(NKNpNGS), &
            & out_data(7,npoints,NKNpNGS), out_viscdata(3,npoints), STAT=allocstat)
        inputconc = 0.0D0
        out_data = 0.0D0
        out_viscdata = 0.0D0
        !--
        IF (verbose) THEN
            WRITE(unito,*) ""
            WRITE(unito,'(A)') "MESSAGE from AIOMFAC: mixture defined, calculating composition points... "
            WRITE(unito,*) ""
        ENDIF
        px = 0
        !set AIOMFAC input and call the main AIOMFAC subroutines for all composition points:
        DO pointi = 1,npoints       !loop over points, changing composition / temperature
            inputconc(1:ncp) = compos2(pointi,1:ncp)
            TKelvin = T_K(pointi)
            !--
            CALL AIOMFAC_inout(inputconc, xinputtype, TKelvin, nspecies, outputvars, outputviscvars, &
                & outnames, errorflagcalc, warningflag)
            !--
            IF (warningflag > 0 .OR. errorflagcalc > 0) THEN
                !$OMP CRITICAL errwriting
                CALL RepErrorWarning(unito, errorflagmix, warningflag, errorflagcalc, pointi, errorind, warningind)
                !$OMP END CRITICAL errwriting
            ENDIF
            nspecmax = MAX(nspecmax, nspecies)  !figure out the maximum number of different species in mixture (accounting for the 
                                                !possibility of HSO4- dissoc. and different species at different data points due to zero mole fractions).
            DO nc = 1,nspecmax                  !loop over species (ions dissociated and treated as individual species):
                out_data(1:6,pointi,nc) = outputvars(1:6,nc)            !out_data general structure: | data columns 1:7 | data point | component no.|
                out_data(7,pointi,nc) = REAL(errorflagcalc, KIND=8)
                out_viscdata(3,pointi) = REAL(errorflagcalc, KIND=8)
                IF (errorflagcalc == 0 .AND. warningflag > 0) THEN      !do not overwrite an errorflag if present!
                    IF (warningflag == 16) THEN                         !a warning that only affects viscosity calc.
                        out_viscdata(3,pointi) = REAL(warningflag, KIND=8)
                    ELSE
                        out_data(7,pointi,nc) = REAL(warningflag, KIND=8)
                        out_viscdata(3,pointi) = REAL(warningflag, KIND=8)
                    ENDIF
                ENDIF
                IF (px(nc) == 0 .AND. out_data(6,pointi,nc) >= 0.0D0) THEN
                    px(nc) = pointi     !use a point for which this component's abundance is given, i.e. mole fraction(nc) > 0.0!
                ENDIF
            ENDDO !nc
            out_viscdata(1:2,pointi) = outputviscvars(1:2)              !out_viscdata general structure: | data columns 1:2 | data point |
        ENDDO !pointi
        
        !
        !==== OUTPUT data-to-file section ==================================================
        !
        !use the name of the input file to create a corresponding output file name from a string like "inputfile_0004.txt"
        i = INDEX(filename, ".txt")
        filename = "AIOMFAC_output_"//filename(i-4:)
        !-- for debugging
        WRITE(unito,'(A)') "................................................................................"
        WRITE(unito,'(A)') "MESSAGE from AIOMFAC: computations successfully performed."
        WRITE(unito,'(4(A))') "Output file, ", TRIM(filename), " created at path: ", TRIM(folderpathout)
        WRITE(unito,'(A)') "................................................................................"
        !create an output ASCII text file with an overall mixture header and individual tables for all components / species (in case of ions)
        fname = TRIM(folderpathout)//TRIM(filename)
        CALL OutputTXT(fname, VersionNo,  nspecmax, npoints, watercompno, cpnameinp(1:nspecmax), T_K(1:npoints), &
            & px(1:nspecmax), out_data, out_viscdata)
        !--
        !>> write output HTML-file
        i = LEN_TRIM(filename)
        filename = filename(1:i-3)//"html"
        fname = TRIM(folderpathout)//TRIM(filename)
        CALL OutputHTML(fname, VersionNo, nspecmax, npoints, watercompno, cpnameinp(1:nspecmax), T_K(1:npoints), &
            & px(1:nspecmax), out_data, out_viscdata)
        !
        !==== TERMINATION section ==========================================================
        !
        DEALLOCATE(inputconc, outputvars, outputviscvars, outnames, out_data, out_viscdata, T_K, cpnameinp, STAT=allocstat)
        IF (ALLOCATED(compos2)) THEN
            DEALLOCATE(compos2, STAT=allocstat)
        ENDIF
    ENDIF !errorind
ENDIF !file valid

WRITE(unito,*) "+-+-+-+-+"
WRITE(unito,'(A)') "Final warning indicator (an entry '00' means no warnings found):"
WRITE(unito,'(I2.2)') warningind
WRITE(unito,*) "+-+-+-+-+"
WRITE(unito,*) ""
WRITE(unito,*) "########"
WRITE(unito,'(A)') "Final error indicator (an entry '00' means no errors found):"
WRITE(unito,'(I2.2)') errorind
WRITE(unito,*) "########"
CLOSE(unito)    !close the error log-file

WRITE(*,*) ""
WRITE(*,'(A,I0)') "MESSAGE from AIOMFAC: end of program; final error indicator: ", errorind
WRITE(*,*) ""
!READ(*,*)      !pause and wait for user action; just for debugging and testing.
!
!==== THE END ======================================================================
!
END PROGRAM Main_IO_driver