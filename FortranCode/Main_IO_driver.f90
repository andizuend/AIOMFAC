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
!*   Gervasi et al. (2020) and Lilek and Zuend (2022, Atmos. Chem. Phys.).              *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
!*                                                                                      *
!*   -> created:        2011  (this file)                                               *
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
!****************************************************************************************

program Main_IO_driver

!module variables:
use Mod_NumPrec, only : wp
use ModSystemProp, only : errorflag_clist, errorflagmix, idCO2, nindcomp, NKNpNGS, SetSystem, topsubno, waterpresent
use ModSubgroupProp, only : SubgroupAtoms, SubgroupNames
use ModMRpart, only : MRdata
use ModSRparam, only : SRdata
use Mod_InputOutput, only : Output_TXT, Output_HTML, RepErrorWarning, ReadInputFile

implicit none
!set preliminary input-related parameters:
integer,parameter :: maxpoints = 1001                       !limit maximum number of composition points for web-version
integer,parameter :: ninpmax = 51                           !set the maximum number of mixture components allowed (preliminary parameter)
!local variables:
character(len=4) :: VersionNo
character(len=200) :: filename
character(len=3000) :: filepath, folderpathout, fname, txtfilein  
character(len=200),dimension(:),allocatable :: cpnameinp    !list of assigned component names (from input file)
character(len=200),dimension(:),allocatable :: outnames
integer :: allocstat, errorind, i, nc, ncp, npoints, nspecies, nspecmax, pointi, &
    & unito, warningflag, warningind, watercompno
integer,dimension(ninpmax) :: px
integer,dimension(:,:),allocatable :: cpsubg                !list of input component subgroups and corresponding subgroup quantities
real(wp) :: TKelvin
real(wp),dimension(:),allocatable :: T_K
real(wp),dimension(:),allocatable :: inputconc, outputviscvars
real(wp),dimension(:,:),allocatable :: composition, compos2, outputvars, out_viscdata
real(wp),dimension(:,:,:),allocatable :: out_data
logical :: filevalid, verbose, xinputtype
logical,dimension(size(errorflag_clist)) :: errflag_list
!...................................................................................

!
!==== INITIALIZATION section =======================================================
!
VersionNo = "3.04"      !AIOMFAC-web version number (change here if minor or major changes require a version number change)
verbose = .true.        !if true, some debugging information will be printed to the unit "unito" (errorlog file)
nspecmax = 0
errorind = 0            !0 means no error found
warningind = 0          !0 means no warnings found
!
!==== INPUT data section ===========================================================
!
!read command line for text-file name (which contains the input parameters to run the AIOMFAC progam):
call get_command_argument(1, txtfilein)
if (len_trim(txtfilein) < 4) then               !no command line argument; use specific input file for tests;
    txtfilein = './Inputfiles/input_0986.txt'   !just use this for debugging with a specific input file, otherwise comment out;
endif
filepath = adjustl(trim(txtfilein))
write(*,*) ""
write(*,'(A,A)') "MESSAGE from AIOMFAC-web: program started, command line argument 1 = ", trim(filepath)
write(*,*) ""
allocate(cpsubg(ninpmax,topsubno), cpnameinp(ninpmax), composition(maxpoints,ninpmax), T_K(maxpoints), STAT=allocstat)
!--
call ReadInputFile(filepath, folderpathout, filename, ninpmax, maxpoints, unito, verbose, ncp, npoints, &
    & warningind, errorind, filevalid, cpnameinp, cpsubg, T_K, composition, xinputtype)
!--
if (filevalid) then
    !
    !==== AIOMFAC initialization and calculation section ===============================
    !
    if (verbose) then
        write(unito,*) ""
        write(unito,'(A)') "MESSAGE from AIOMFAC: input file read, starting AIOMFAC mixture definitions and initialization... "
        write(unito,*) ""
    endif
    !load the MR and SR interaction parameter data:
    call MRdata()         !initialize the MR data for the interaction coeff. calculations
    call SRdata()         !initialize data for the SR part coefficient calculations
    call SubgroupNames()  !initialize the subgroup names for the construction of component subgroup strings
    call SubgroupAtoms()
    !--
    !set mixture system properties (composition-independent properties) based on the data from the input file:
    call SetSystem(1, .true., ncp, cpnameinp(1:ncp), cpsubg(1:ncp,1:topsubno) )

    !check whether water is present in the mixture and as which component number:
    watercompno = 0
    if (waterpresent) then
        watercompno = maxloc(cpsubg(1:ncp,16), dim=1)   !usually = 1
    endif
    if (idCO2 > 0) then             !add CO2 as a (non-input) neutral component name
        cpnameinp(idCO2+1:ncp+1) = cpnameinp(idCO2:ncp)
        cpnameinp(idCO2) = 'CO2(aq)'   
    endif
    !transfer composition data to adequate array size:
    allocate(compos2(npoints,ncp), STAT=allocstat)
    do nc = 1,ncp
        compos2(1:npoints,nc) = composition(1:npoints,nc)
    enddo
    deallocate(cpsubg, composition, STAT=allocstat)
    
    if (errorflagmix /= 0) then     !a mixture-related error occurred:
        call RepErrorWarning(unito, errorflagmix, warningflag, errflag_list, i, errorind, warningind)
    endif

    if (errorind == 0) then         !perform AIOMFAC calculations; else jump to termination section
        !--
        allocate(inputconc(nindcomp), outputvars(6,NKNpNGS), outputviscvars(2), outnames(NKNpNGS), &
            & out_data(7,npoints,NKNpNGS), out_viscdata(3,npoints), STAT=allocstat)
        inputconc = 0.0_wp
        out_data = 0.0_wp
        out_viscdata = 0.0_wp
        !--
        if (verbose) then
            write(unito,*) ""
            write(unito,'(A)') "MESSAGE from AIOMFAC: mixture defined, calculating composition points... "
            write(unito,*) ""
        endif
        px = 0
        !set AIOMFAC input and call the main AIOMFAC subroutines for all composition points:
        do pointi = 1,npoints       
            !loop over points, changing composition / temperature
            inputconc(1:ncp) = compos2(pointi,1:ncp)
            TKelvin = T_K(pointi)
            !--
            call AIOMFAC_inout(inputconc, xinputtype, TKelvin, nspecies, outputvars, outputviscvars, &
                & outnames, errflag_list, warningflag)
            !--
            if (warningflag > 0 .or. any(errflag_list)) then
                !$OMP CRITICAL errwriting
                call RepErrorWarning(unito, errorflagmix, warningflag, errflag_list, pointi, errorind, warningind)
                !$OMP end CRITICAL errwriting
            endif
            nspecmax = max(nspecmax, nspecies)  !figure out the maximum number of different species in mixture (accounting for the 
                                                !possibility of HSO4- dissoc. and different species at different data points due to zero mole fractions).
            do nc = 1,nspecmax                  !loop over species (ions dissociated and treated as individual species):
                out_data(1:6,pointi,nc) = outputvars(1:6,nc)                !out_data general structure: | data columns 1:7 | data point | component no.|
                out_data(7,pointi,nc) = real(findloc(errflag_list(:), value = .true., dim=1), kind=wp)
                out_viscdata(3,pointi) = real(findloc(errflag_list(:), value = .true., dim=1), kind=wp)
                if (abs(out_viscdata(3,pointi) -18.0_wp) < 0.1_wp) then       !error 18 is only an error/warning for viscosity prediction, not activity coeff.
                    out_data(7,pointi,nc) = 0.0_wp
                endif
                if ((.not. any(errflag_list)) .and. warningflag > 0) then   !do not overwrite an errorflag if present!
                    if (warningflag == 16) then                             !a warning that only affects viscosity calc.
                        out_viscdata(3,pointi) = real(warningflag, kind=wp)
                    else
                        out_data(7,pointi,nc) = real(warningflag, kind=wp)
                        out_viscdata(3,pointi) = real(warningflag, kind=wp)
                    endif
                endif
                if (px(nc) == 0 .and. out_data(6,pointi,nc) >= 0.0_wp) then
                    px(nc) = pointi     !use a point for which this component's abundance is given, i.e. mole fraction(nc) > 0.0!
                endif
            enddo !nc
            out_viscdata(1:2,pointi) = outputviscvars(1:2)              !out_viscdata general structure: | data columns 1:2 | data point |
        enddo !pointi
        
        !
        !==== OUTPUT data-to-file section ==================================================
        !
        !use the name of the input file to create a corresponding output file name from a string like "inputfile_0004.txt"
        i = index(filename, ".txt")
        filename = "AIOMFAC_output_"//filename(i-4:)
        !-- for debugging
        write(unito,'(A)') "................................................................................"
        write(unito,'(A)') "MESSAGE from AIOMFAC: computations successfully performed."
        write(unito,'(4(A))') "Output file, ", trim(filename), " created at path: ", trim(folderpathout)
        write(unito,'(A)') "................................................................................"
        !create an output ASCII text file with an overall mixture header and individual tables for all components / species (in case of ions)
        fname = trim(folderpathout)//trim(filename)
        call Output_TXT(fname, VersionNo,  nspecmax, npoints, watercompno, cpnameinp(1:nspecmax), T_K(1:npoints), &
            & px(1:nspecmax), out_data, out_viscdata)
        !--
        !>> write output HTML-file
        i = len_trim(filename)
        filename = filename(1:i-3)//"html"
        fname = trim(folderpathout)//trim(filename)
        call Output_HTML(fname, VersionNo, nspecmax, npoints, watercompno, cpnameinp(1:nspecmax), T_K(1:npoints), &
            & px(1:nspecmax), out_data, out_viscdata)
        !
        !==== TERMINATION section ==========================================================
        !
        deallocate(inputconc, outputvars, outputviscvars, outnames, out_data, out_viscdata, T_K, cpnameinp, STAT=allocstat)
        if (allocated(compos2)) then
            deallocate(compos2, STAT=allocstat)
        endif
    endif !errorind
endif !file valid

write(unito,*) "+-+-+-+-+"
write(unito,'(A)') "Final warning indicator (an entry '00' means no warnings found):"
write(unito,'(I2.2)') warningind
write(unito,*) "+-+-+-+-+"
write(unito,*) ""
write(unito,*) "########"
write(unito,'(A)') "Final error indicator (an entry '00' means no errors found):"
write(unito,'(I2.2)') errorind
write(unito,*) "########"
close(unito)    !close the error log-file

write(*,*) ""
write(*,'(A,I0)') "MESSAGE from AIOMFAC: end of program; final error indicator: ", errorind
write(*,*) ""
!read(*,*)      !pause and wait for user action; just for debugging and testing.
!
!==== the end ======================================================================
!
end program Main_IO_driver