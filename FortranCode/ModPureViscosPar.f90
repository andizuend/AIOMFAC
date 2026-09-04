!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module providing access to the pure component viscosity parameter arrays.          *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend, Natalie Gervasi, Zixuan Shen, Joshua Penikis                           *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2012-09-07                                                      *
!*   -> latest changes: 2026-08-22                                                      *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  subroutine DeRieux_Tno_Est                                                      *
!*   -  subroutine Tg_ML_Armeli                                                         *
!*   -  subroutine VogelTemp                                                            *
!*   -  subroutine PureCompViscosity                                                    *
!*                                                                                      *
!****************************************************************************************
    
module ModPureCompViscos

use Mod_kind_param, only : wp

implicit none
!..................................

private     !as default

!public parameter arrays:
integer,dimension(1500),public :: CorrelEqNo
real(wp),dimension(1500,2),public :: CorrelTrange       !structure: (lower T limit, upper T limit) in [K]
real(wp),dimension(1500,5),public :: ViscosCorrelPar    !structure: (component no., param. A, B, C, D, E)
real(wp),public :: TgScale
real(wp),public :: TgScalePlot = -1.0_wp                !initialize with an unphysical value
real(wp),dimension(:),allocatable,private :: Tg_list    ![K] array storing the Tg values of the components for this system after the first look up / calculation
real(wp),dimension(:),allocatable,private :: TgML       ![K] array storing Tg values determined based on Machine Learning method by Armeli et al. (2023)

!module procedures:
public :: PureCompViscosity

private :: VogelTemp, DeRieux_Tno_Est, Tg_ML_Armeli

!$OMP THREADPRIVATE(Tg_list, TgML, TgScale, TgScalePlot)

!========================================================================================================== 
    contains
!========================================================================================================== 
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Calculation of Tg in [K] using the model from DeRieux et al. 2018 (ACP).           *
    !*   [cP] = 0.01 [Poise] = 0.01 [g/(cm.s)] = 0.001 [Pa.s] = 0.001 [N.s/m^2]             *
    !*                                                                                      *
    !*   :: Authors & Copyright ::                                                          *
    !*   Natalie Gervasi, Andi Zuend,                                                       *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *                       
    !*                                                                                      *
    !*   -> created:        2019-02-03                                                      *
    !*   -> latest changes: 2025-07-06                                                      *      
    !*                                                                                      *
    !****************************************************************************************
    pure subroutine DeRieux_Tno_Est(ind, expTg, Tg)

    use ModSubgroupProp, only : O2C_H2C_component

    implicit none
    !..................................
    !interface variables:
    integer,intent(in) :: ind           !index no. of component in ITAB
    real(wp),intent(in) :: expTg        !the experimental Tg (if available)
    real(wp),intent(out) :: Tg          !glass transition temp [K]
    !local variables:
    real(wp),parameter :: nCO = 12.13_wp, bC = 10.95_wp, bH = -41.82_wp, bCH = 21.61_wp, bO = 118.96_wp, bCO = -24.38_wp
    real(wp),parameter :: nCO_k = 1.96_wp, bC_k = 61.99_wp, bH_k = -113.33_wp, bCH_k = 28.74_wp
    real(wp) :: cpCarbon, cpHydrogen, cpOxygen, cpNitrogen, cpSulfur, OtoC, HtoC, NtoC, StoC
    !..................................
    
    call O2C_H2C_component(ind, cpCarbon, cpHydrogen, cpOxygen, cpNitrogen, cpSulfur, OtoC, HtoC, NtoC, StoC)

    ! Tg is either calculated or an experimentally determined value is used
    if (expTg > -999.0_wp) then ! use experimental value
        Tg = expTg
    else
        !choose the constants for the Tg model based on whether the compound contains oxygen or not
        if (cpHydrogen < 1.0_wp) then     !(AZ: added 2025-07-03); fix to allow for predictions when the compound contains no H (in those cases it often contains lots of N, which is ignored below, so the prediction may be poor anyways)
            cpHydrogen = 1.0_wp
        endif
        if (cpOxygen < 1.0_wp) then
            Tg = bC_k*(nCO_k + log(cpCarbon)) + bH_k*log(cpHydrogen) + bCH_k*log(cpCarbon)*log(cpHydrogen)
        else
            ! Shiraiwa et al. group model used to calculate Tg (DeRieux et al. 2018), eqn (2)
            Tg = bC*(nCO + log(cpCarbon)) + bH*log(cpHydrogen) + bCH*log(cpCarbon)*log(cpHydrogen) + &
                & bO*log(cpOxygen) + bCO*log(cpCarbon)*log(cpOxygen)
        endif
    endif

    end subroutine DeRieux_Tno_Est
    !========================================================================================================== 
    
    

    !***********************************************************************************************
    !*   :: Purpose ::                                                                             *
    !*   Calculation of Tg in [K] using machine learning prediction based on SMILES string of      *
    !*   components by Armeli et al. 2023 (ACS).                                                   *           
    !*                                                                                             *
    !*   :: Authors & Copyright ::                                                                 *
    !*   Zixuan Shen, Andi Zuend, Joshua Penikis                                                   *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                                 *         
    !*                                                                                             *
    !*   -> latest changes: 2026-08-24                                                             * 
    !*                                                                                             *
    !***********************************************************************************************
    subroutine Tg_ML_Armeli(SMILESfname, outfilename)
    
    use ModOScommands, only : isWindowsOS, f_replace_text
    use ModPureCompProp, only : RelPathModel
    use ModSystemProp, only : errorflagmix
    
    implicit none
    
    character(len=*),intent(in) :: SMILESfname, outfilename
    !local variables:
    character(len=200) :: pathTools
    character(len=300) :: command, command2
    integer :: Estat, Cstat
    !.................................
    
    write(*,'(A,/)') 'Note: updating pure-component Tg values via Python call of TgML_Armeli method. This could take a few seconds...'
    
    pathTools = RelPathModel//'TgML_Armeli/'      !path from AIOMFAC_Proj directory to the Tools for TgML folder
    
    command = trim(pathTools)//'.venv/Scripts/python.exe '//trim(pathTools)//'TgML_SMILES.py '//trim(pathTools)//'InputFiles/'//trim(SMILESfname)//' ' &
        & //trim(pathTools)//'OutputFiles/'//trim(outfilename)
    
    if (isWindowsOS) then
        command2 = f_replace_text(command, '/', '\')                                          !replace forward by backslashes for Windows commands
    else !Linux, MacOS, ...
        command2 = f_replace_text(command, '.venv/Scripts/python.exe', '.venv/bin/python')    !no need for executable file's extension on Linux
    endif
    
    call execute_command_line(trim(command2), exitstat=Estat, cmdstat=Cstat)
    if (Estat == 0) then
        !$OMP critical
        write(*,'(A,/)') 'Note from Tg_ML_Armeli: completed the pure-component Tg values update.'
        !$OMP end critical 
    else
        errorflagmix = 26
        !$OMP critical
        write(*,'(A,/)') ''
        write(*,'(A,/)') 'ERROR in Tg_ML_Armeli: Script access issue. Unsuccessful in running the "TgML_SMILES.py" &
            &script inside subfolder "TgML_Armeli". &
            &Ensure that the folder is present and that a virtual Python environment named ".venv" is present as a &
            &(hidden) subfolder. See also information in "requirements_console_script_AZ.txt".'
            call sleep(1)
        !$OMP end critical 
    endif
    
    end subroutine Tg_ML_Armeli
    !==========================================================================================================
    
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Calculation of Tno in [K] using the Vogel-Tamman-Fulcher (VTF) equation in the     *
    !*   form used by Angell (2002) and DeRieux et al. 2018.                                *
    !*   Fragility (D) is determined by the temperature of the run.                         *
    !*                                                                                      *
    !*   :: Authors & Copyright ::                                                          *
    !*   Natalie Gervasi, Andi Zuend,                                                       *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *                                  
    !*                                                                                      *
    !****************************************************************************************
    pure subroutine VogelTemp(Tg, TempK, D, Tno)

    implicit none

    !interface
    real(wp), intent(in) :: Tg, TempK       ! glass transition temperature, run temperature [K]
    real(wp), intent(out) :: Tno, D         ! Vogel temperature [K], fragility parameter [-]
    !........................  

    if (TempK < Tg) then                    ! pure-component substance below Tg that was once fragile will behave more strongly
        D = 30.0_wp
    else
        D = 10.0_wp                         ! reasonable guess for organics (Shiraiwa et al., 2017)
    endif
    Tno = (39.17_wp*Tg) / (D + 39.17_wp)    ! (DeRieux et al. 2018) eqn (7)

    end subroutine VogelTemp
    !==========================================================================================================
    

    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Calculation of pure-component dynamic viscosity (eta0) in Pascal seconds [Pa.s]    *
    !*   at a given temperature T [K]. Output is in form of ln(eta0/[Pa s]).                *
    !*   [cP] = 0.01 [Poise] = 0.01 [g/(cm.s)] = 0.001 [Pa.s] = 0.001 [N.s/m^2]             *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend, Natalie Gervasi, Joshua Penikis                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2012-09-07                                                      *
    !*   -> latest changes: 2026-08-24                                                      *
    !*                                                                                      *
    !****************************************************************************************
    subroutine PureCompViscosity(ind, TempK, ln_eta0, iflag, Tglass, fragility)
    
    use ModSystemProp, only : CompN, ITAB, nindcomp, waterpresent, errorflag_clist, maxsmileslength
    use Mod_InputOutput, only : cpsmiles, armeliON
    use ModPureCompProp, only : lookup_Tg, append_purecomp_entry, RelPathModel
    use ModOScommands, only : isWindowsOS, f_replace_text

    implicit none
    !..................................
    !interface variables:
    integer,intent(in) :: ind
    real(wp),intent(in) :: TempK
    real(wp),intent(out) :: Tglass, ln_eta0, fragility
    integer,intent(out) :: iflag
    !local variables:
    integer :: cpn, equationNo, i, ind_water, unsmiles, un1, nlines, istat, num, fsize, newfsize, &
        clockstart, clockcount, clockrate, exstat, cmdstat, unresult
    real(wp),parameter :: ln10 = log(10.0_wp), ln_bwater = log(1.3788E-4_wp)
    real(wp) :: a, b, c, d, e, ln_b, Tg, Tvog, Tgest, Tg_value, Tg_read, r, elapsed_s
    logical :: Tg_found, fileexists
    character(len=100) :: SMILES_input_file, TgML_output_file
    character(len=4) :: casenumber
    character(len=7) :: rdwr_status
    character(len=5) :: cmd_res
    character(len=500) :: watchdog_cmd, cmd_line, tmp_file, inp_file, outp_file, tmp2, inp2
    character(len=maxsmileslength) :: smiles_input, smiles_read, smiles_batch
    !..................................

    !initialize:
    if (ind > 0) then
        if (ITAB(ind,16) > 0) then          !is water
            cpn = 401
        else if (ITAB(ind,173) > 0) then    !is CO2(aq)
            cpn = 402
        else
            cpn = compN(ind)                !organic or water component ID (when defined)
        endif
    else if (ind == -1) then                !this should not happen
        cpn = 401                           !water
    endif

    select case(cpn)
    case(401) !water
        if (TempK >= 230.0_wp) then
            equationNo = 12
            a = 225.66_wp
            b = 1.3788E-4_wp
            ln_b = ln_bwater
            c = -1.6433_wp
            d = 0.0_wp
            e = 0.0_wp
            CorrelTrange(cpn,1) = 230.0_wp
            CorrelTrange(cpn,2) = 495.0_wp
        else
            equationNo = 10
            a = 136.0_wp                    !Tg of water from experiments (see, e.g. Koop et al., 2011, PCCP);
            b = -999.0_wp
            c = -999.0_wp
            d = -999.0_wp
            e = -999.0_wp
            CorrelTrange(cpn,1) = 136.0_wp
            CorrelTrange(cpn,2) = 230.0_wp
        endif
        
    case(402) !CO2(aq)
        equationNo = 10
        a = 130.0_wp
        b = -999.0_wp
        c = -999.0_wp
        d = -999.0_wp
        e = -999.0_wp
        CorrelTrange(cpn,1) = 130.0_wp
        CorrelTrange(cpn,2) = 1000.0_wp
        
    case(1500, 9999) !for system input of organics from file
        cpn = 1500
        equationNo = 10
        a = -999.0_wp
        b = -999.0_wp
        c = -999.0_wp
        d = -999.0_wp
        e = -999.0_wp
        CorrelTrange(cpn,1) = 1.0_wp
        CorrelTrange(cpn,2) = 1000.0_wp
        
    case default
        !assign parameter values to variables used below:
        equationNo = CorrelEqNo(cpn)
        a = ViscosCorrelPar(cpn,1)
        b = ViscosCorrelPar(cpn,2)
        c = ViscosCorrelPar(cpn,3)
        d = ViscosCorrelPar(cpn,4)
        e = ViscosCorrelPar(cpn,5)
        
    end select
    
    !initialize Tg_list:
    if (.not. allocated(Tg_list)) then
        allocate(Tg_list(nindcomp))
        Tg_list = -111.0_wp                                                    !initialize to unrealistic value
    endif
    
    if (armeliON .and. cpsmiles(ind) /= "") then
        !------------
        !initialize TgML array
        if (.not. allocated(TgML)) then
            allocate(TgML(size(cpsmiles)))
            TgML = -111.0_wp
        endif
        !--
        !compare requested temperature with temperature range of parameterization:
        if (TempK >= CorrelTrange(cpn,1) .and. TempK <= CorrelTrange(cpn,2)) then 
            !valid T range for the correlation equation selected
        else !use the default method with glass transition temperature estimation
            equationNo = 10
            a = -999.0_wp
            b = -999.0_wp
            c = -999.0_wp
            d = -999.0_wp
            e = -999.0_wp
        endif
        iflag = 0 !0 means no errors
    
        !determine apppropriate value in Tg range
        if (ind > 0) then
            Tg = Tg_list(ind)
            if (Tg < 0.0_wp) then                                                   !need to look up or compute Tg
                !---
                !check whether water is present 
                if (waterpresent) then
                    ind_water = 1
                else
                    ind_water = 0
                endif
            
                if ( (ind > ind_water) .and. ( .not. (errorflag_clist(19) .or. errorflag_clist(22)) ) ) then    !only enter when not water and no error with pure-component .txt file
                
                    smiles_input = trim(cpsmiles(ind))                              !pull individual smiles
                    
                    !check if Tg has already been predicted for SMILES
                    call lookup_Tg(smiles_input, Tg_value, Tg_found)
                
                    if (Tg_found) then                                              !Tg_value has valid value
                        !Tg already in table; continue below...
                    else
                        call random_seed()
                        call random_number(r)                                       !generate random number
                        num = floor(r * 10000.0_wp)                                 !make 4 digits
                        write(casenumber,'(I4.4)') num                              !assign unique filenumber
                        SMILES_input_file = "input_"//casenumber//"_SMILES.txt"
                        TgML_output_file = "output_"//casenumber//"_SMILES.txt"
                        inp_file = RelPathModel//"TgML_Armeli/InputFiles/"//trim(SMILES_input_file)
                        tmp_file = f_replace_text(inp_file, ".txt", ".tmp")           !also make a temporary file so that watchdog cannot process it before it is completed and sized.
                        tmp_file = f_replace_text(tmp_file, "InputFiles", "OutputFiles") 
                        open (newunit = unsmiles, file = trim(tmp_file), status = "new", action = "readwrite", iostat = istat)   !write a temporary SMILES input file
                        if (istat /= 0) write(*,*) "@PureCompViscosity: issue when opening temporary SMILES_input_file"
                        do i = 1, size(cpsmiles)                                    !batch all SMILES from this system
                            smiles_batch = trim(cpsmiles(i))
                            if (len_trim(smiles_batch) > 0) then                    !contains SMILES
                                call lookup_Tg(smiles_batch, Tg_value, Tg_found)    !check whether Tg already predicted for SMILES 
                                if (.not. Tg_found) then                            !write SMILES to batch file 
                                    write(unsmiles,'(A)') trim(smiles_batch) 
                                endif 
                            endif
                        enddo
                        close(unsmiles)
                    
                        inquire(file = trim(tmp_file), size = fsize, readwrite = rdwr_status)  !take initial size of temporary file (= initial size input file)

                        !set file access permissions to permissive, then 
                        !move the .tmp file to .txt so that watchdog may process it:
                        if (.not. isWindowsOS) then !Linux
                            cmd_line = 'chmod 777 '//trim(tmp_file)//' && mv '//trim(tmp_file)//' '//trim(inp_file)
                            call execute_command_line(trim(cmd_line), exitstat = exstat, cmdstat = cmdstat)
                            if (cmdstat /= 0) write(*,*) "ERROR in PureCompViscosity, moving file: cmdstat = ", cmdstat
                        else
                            inp2 = f_replace_text(inp_file, "/", "\")
                            tmp2 = f_replace_text(tmp_file, "/", "\")
                            cmd_line = 'icacls '//trim(tmp_file)//' /grant Users:F > NUL & move /y '//trim(tmp2)//' '//trim(inp2)//' > NUL'
                            call execute_command_line(trim(cmd_line), exitstat = exstat, cmdstat = cmdstat) 
                            if (cmdstat /= 0) write(*,*) "ERROR in PureCompViscosity, moving file: cmdstat = ", cmdstat
                        endif
                    
                        inquire(file = RelPathModel//"TgML_Armeli/TgML_SMILES_watchdog.py", exist = fileexists)          !inquire about watchdog existing
                        if (fileexists) then                                        !method *should* be running in background
                        
                            !inquire whether TgML watchdog process is running using command line call:
                            if (isWindowsOS) then
                                watchdog_cmd = 'powershell -NoProfile -Command "[bool](Get-CimInstance Win32_Process -Filter \"Name like'//" &
                                    & 'python%%'"//'\" | Where-Object { $_.ProcessId -ne $PID -and $_.CommandLine -like '// &
                                    & "'*TgML_SMILES_watchdog.py*'"//" }) | ForEach-Object { $_.ToString().ToLower() } | Out-File -FilePath &
                                    & "//RelPathModel//'Auxiliary/cmd_result.txt'//' -Encoding ascii"'
                            else    !Linux
                                watchdog_cmd = '(pgrep -f "TgML_SMILES_watchdog.py" > /dev/null && echo "true" || echo "false") > '//RelPathModel//'Auxiliary/cmd_result.txt'
                            endif
                            call execute_command_line(trim(watchdog_cmd), wait = .true., cmdstat = cmdstat)
                            !read the file written by the command line execution:
                            open(newunit = unresult, file = RelPathModel//"Auxiliary/cmd_result.txt", status = 'unknown', action="readwrite", iostat = istat)
                            read(unresult,'(A)') cmd_res
                            close(unresult, status = "delete")

                            if (cmdstat /= 0) then
                                errorflag_clist(20) = .true.                        !command couldn't run, contact admin
                            elseif (trim(cmd_res) == "false") then                  !watchdog isn't running, default to calling Tg_ML_Armeli
                                call Tg_ML_Armeli(trim(SMILES_input_file), trim(TgML_output_file))
                                if (.not. isWindowsOS) errorflag_clist(21) = .true. !on Linux, send the user a warning to contact admin (since watchdog is down)
                            else                                                    !cmd_res == "true", i.e. watchdog is running
                                call system_clock(count = clockstart, count_rate = clockrate)   !start timer
                                do  !until exit
                                    inquire(file = trim(inp_file), size = newfsize, readwrite = rdwr_status)  !check new file size
                                    if (newfsize > fsize + 1) then
                                        exit                                        !input file changed by > 1 storage unit (byte?) -> Tg method ran, output file should be generated
                                    endif
                                    call system_clock(count = clockcount)
                                    elapsed_s = real(clockcount - clockstart, wp) / real(clockrate, wp)
                                    if (elapsed_s > 5.0_wp) then
                                        exit                                        !time out after 5 seconds since likely there is another issue (access permissions)
                                    endif
                                enddo
                            endif
                        
                        else    !watchdog process is not running in background
                            !as slower alternative: run Tg prediction based on a direct call to ML method by Armeli et al.
                            call Tg_ML_Armeli(SMILES_input_file, TgML_output_file)  
                        endif

                        !check whether output file exists and read its content if true:
                        outp_file = RelPathModel//"TgML_Armeli/OutputFiles/"//trim(TgML_output_file)
                        inquire(file = trim(outp_file), exist = fileexists, readwrite = rdwr_status)
                    
                        if (fileexists .and. trim(rdwr_status) /= 'NO' .and. trim(rdwr_status) /= 'no') then
                        
                            open(newunit = un1, file = trim(outp_file), action="readwrite", status = "unknown")
                            !calculate the number of lines of Tg output
                            nlines = 0
                            do  !until exit
                                read(un1,*,iostat = istat)
                                if (istat /= 0) exit
                                nlines = nlines + 1
                            enddo
                            rewind(un1)                                             !read the Tg output again (from the beginning)
                            do i = 1, nlines                                        !read the value for Tg and keep the corresponding SMILES
                                read(un1,*,iostat = istat) Tg_read, smiles_read     !read the data line
                                if (Tg_read > -90.0_wp) then                        !a valid Tg was determined
                                    call append_purecomp_entry(smiles_read, Tg_read)!append to table and .txt file
                                endif
                                if (istat < 0) then                                 !end of file reached
                                    exit                                            !leave do-loop
                                elseif (istat > 0) then                             !error occurred
                                    write(un1,*) "an error occurred while reading the Tg output; istat = ", istat
                                    !read(un1,*)                                    !wait for user action
                                endif
                            enddo
                            close(un1, status = "delete")                           !delete the output file after reading
                        
                        else    !no output file was generated
                            Tg_value = -9999.999999_wp                              !indicate Tg could not be determined
                            errorflag_clist(23) = .true.                            !error attempting to read Armeli output file
                        endif
                    
                        !clean up: delete the input file after temporary use
                        inquire(file = trim(inp_file), exist = fileexists, readwrite = rdwr_status)
                        if (fileexists .and. trim(rdwr_status) /= 'NO' .and. trim(rdwr_status) /= 'no') then
                            open(newunit = unsmiles, file = trim(inp_file), action="readwrite", status = "old")
                            close(unsmiles, status = "delete")    
                        endif
                    
                        call lookup_Tg(trim(smiles_input), Tg_value, Tg_found)      !Tg now updated for the SMILES
                    endif
                
                    if (Tg_found .and. Tg_value > -90.0_wp) then
                        TgML(ind) = Tg_value                                        !assign Tg to TgML at ind
                    else    !Tg = -99.00 returned indicates that SMILES invalid
                        errorflag_clist(24) = .true.
                        TgML(ind) = 10.0_wp                                         !assign an incorrect but positive value nevertheless so the remaining program can run to completion and issue an error.
                    endif
                    Tg = TgML(ind)
                    !---
                else
                    !---
                    select case(equationNo)
                    case(10, 11)
                        Tgest = a
                    case default
                        if (cpn == 401) then
                            Tgest = 136.0_wp
                        else
                            Tgest = -9999.0_wp
                        endif
                    end select
                    call DeRieux_Tno_Est(ind, Tgest, Tg)
                endif
                Tg_list(ind) = Tg                                                   !save for next calcuation with this system
                !---
            endif !Tg check
            
            Tglass = Tg*TgScale
            TgScalePlot = Tg
            
        else if (ind == -1) then    !use water properties for ions
            Tg = 136.0_wp
            Tglass = Tg*TgScale
            TgScalePlot = Tg
        else
            TgScalePlot = -1.0_wp
        endif

        !compute ln_eta0 with given equation for the component
        select case(equationNo)
        case(1)                                                                 !Correlation Equation #1: Daubert and Danner; ln_eta0 in Pa.s, TempK in Kelvin:
            ln_eta0 = a + b/TempK + c*log(TempK) + d*TempK**e
        case(3)                                                                 !for estimated values at a single temperature value only.
            ln_eta0 = a*ln10
        case(4)                                                                 !for pure-component values found via experiment and quoted in the literature. eta in Pa.s
            ln_eta0 = log(a)
        case(5)                                                                 !Correlation equation for acetone. eta in Pa.s T = 293.15K-318.15K (data from Hafez)
            ln_eta0 = log(a*TempK + 0.0011_wp)
        case(6)                                                                 !Correlation equation for propanoic acid. eta in Pa.s T = 293.15K - 325.15K (data from Rattan)
            ln_eta0 = log(a*TempK + 0.0037_wp)
        case(7)                                                                 !Correlation equation of Vogel-Fulcher-Tammann (VFT); Ollett and Parker (1990)
            ln_eta0 = log(a) + b/(TempK - c)
        case(8)                                                                 !Correlation equation of Vogel-Fulcher-Tammann (VFT); Angell et al. (1982);
            ln_eta0 = log(a*1.0E-3_wp) + b/(TempK - c)                          !incl. conversion from centi-poise to Pa.s units
        case(9)                                                                 !Correlation equation of van Velzen et al. ; see Viswanath et al. book (2007) "Viscosity of liquids"
            ln_eta0 = log(1.0E-3_wp*10.0_wp**(a*((1.0_wp/TempK) - (1.0_wp/b)))) !(a is parameter B and b is T0 in van Velzen et al. eq.) !conversion from centi-poise to Pa.s units
        case(10, 11)                                                            !Vogel-Tammann-Fulcher (VFT), Angell (1991) using DeRieux et al. (2018) constants and DeRieux Tg (case 10) or experimental/predicted Tg values (case 11)                                                                                              
            call VogelTemp(Tglass, TempK, fragility, Tvog)
            if (Tvog >= TempK) then
                iflag = 1                                                       !1 = outside valid temperature range!
                ln_eta0 = 600.0_wp
            else
                ln_eta0 = ln10*( -5.0_wp + 0.434_wp*(fragility*Tvog/(TempK - Tvog)) )   !DeRieux et al. 2018, eqn (6)
            endif
        case(12)
            ln_eta0 = ln_b + c*log((TempK/a) - 1.0_wp)
        end select
        iflag = 0       !0 indicates no errors
        !------------
    else
        !------------
        !for use in AIOMFAC-web when only using DeRieux et al. estimation method:
        if (TempK >= CorrelTrange(cpn,1) .and. TempK <= CorrelTrange(cpn,2)) then   !valid
            iflag = 0           !0 means no errors
            if (ind > 0) then
                Tg = Tg_list(ind)
            endif
            if (Tg < 0.0_wp) then
                select case(equationNo)
                case(10, 11)
                    Tgest = a
                case default
                    if (cpn == 401) then
                        Tgest = 136.0_wp
                    else
                        Tgest = -9999.0_wp
                    endif
                end select
                call DeRieux_Tno_Est(ind, Tgest, Tg)
                Tg_list(ind) = Tg                                               !save for next calcuation with this system
            endif
            Tglass = Tg*TgScale
            TgScalePlot = Tg
            
            !compute eta with given equation for the component:
            select case(equationNo)
            case(12)
                ln_eta0 = ln_b + c*log((TempK/a) - 1.0_wp)
            case(10)     
                call VogelTemp(Tglass, TempK, fragility, Tvog)                  !Vogel-Tammann-Fulcher (VFT), Angell (1991) using DeRieux et al. (2018) constants and DeRieux Tg
                if (Tvog >= TempK) then
                    iflag = 1                                                   !1 = outside valid temperature range!
                    ln_eta0 = 600.0_wp
                else
                    ln_eta0 = ln10*( -5.0_wp + 0.434_wp*(fragility*Tvog/(TempK - Tvog)) )    !DeRieux et al. 2018, eqn (6)
                endif
            end select
            
        else
            iflag = 1                                                           !1 = outside valid temperature range!
        endif
    endif
        
    end subroutine PureCompViscosity
    !==========================================================================================================

end module ModPureCompViscos
