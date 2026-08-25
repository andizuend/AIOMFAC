    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Module with subroutines and functions to carry out different operating system      *
    !*   commands working both on Windows and Linux platforms. Main uses include copying of *
    !*   files via commands executed in a terminal command line.                            *
    !*   Some functions require Python to be installed and available from command line on   *
    !*   the local system.                                                                  *
    !*                                                                                      *
    !*   :: Authors & Copyright ::                                                          *
    !*   Andi Zuend                                                                         *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2020                                                            *
    !*   -> latest changes: 2026-08-24                                                      *
    !*                                                                                      *
    !*   :: List of subroutines and functions contained in this module:                     *
    !*   --------------------------------------------------------------                     *
    !*   -  SUBROUTINE copy_file                                                            *
    !*   -  FUNCTION   isWindowsPlatform                                                    *
    !*   -  FUNCTION   Replace_Text                                                         *
    !*                                                                                      *
    !****************************************************************************************
    module ModOScommands

    implicit none
    
    logical,public :: isWindowsOS

    public :: copy_file, isWindowsPlatform, Replace_Text

    contains

    !--------------------------------------------------------------------------------------
    !   This subroutine copies a file 'file_name' to 'file_name_new' via a system command.
    !   The file_name string can include a relative path on the sytem.
    !--------------------------------------------------------------------------------------
    subroutine copy_file(file_name, file_name_new)

    use ModStringFunctions

    implicit none

    character(len=*),intent(in) :: file_name, file_name_new
    !local:
    character(len = 100 + 2*max(len(file_name_new), len(file_name))) :: command, command2
    integer :: Estat, Cstat
    !.................................

    if (isWindowsOS) then
        command = 'copy "'//trim(file_name) //'" "'//trim(file_name_new)//'" > NUL'
        command2 = trim(Replace_Text(command, "/", "\"))    !replace forward- by backslashes for Windows commands
        call execute_command_line(trim(command2), exitstat=Estat, cmdstat=Cstat)            !first copy file
        if (Cstat /= 0) write(*,*) "ERROR in copy_file 1: Cstat = ", Cstat
        command = 'icacls '//trim(file_name_new)//' /grant Users:F'                         !then set access permissions on Windows to full control (permissive)
        call execute_command_line(trim(command), exitstat=Estat, cmdstat=Cstat)
        if (Cstat /= 0) write(*,*) "ERROR in copy_file 2: Cstat = ", Cstat
    else !on a LINUX OS?
        command2 = 'cp -p"'//trim(file_name) //'" "'//trim(file_name_new)//'" > NUL'        !copy while keeping access permissions as of the original file
        call execute_command_line(trim(command2), exitstat=Estat, cmdstat=Cstat)
        if (Cstat /= 0) write(*,*) "ERROR in copy_file 1: Cstat = ", Cstat
    endif

    end subroutine copy_file
    !------------------------------------------------------------------------------------


    !------------------------------------------------------------------------------------ 
    !A logical querry function that returns 'true' if the operating system / platform is a version of Windows. 
    function isWindowsPlatform() 
     
    implicit none 
     
    logical             :: isWindowsPlatform 
    character(len=32)   :: os_val 
    integer             :: val_len, status 
    !............................... 
     
    !Get the OS-specific 'OS' variable: 
    call get_environment_variable("OS", os_val, val_len, status) 
     
    if (status == 0) then 
        select case(os_val(1:3)) 
        case('Win', 'win', 'WIN') 
            isWindowsPlatform = .true. 
        case default    !otherwise likely Linux or Mac OS 
            isWindowsPlatform = .false. 
        end select 
    else  !otherwise likely Linux
        !write(*,'(A)') "ERROR: Could not determine operating system in function isWindowsPlatform!"
        !write(*,*) "status value was: ", status
        !read(*,*) !wait for user action
        isWindowsPlatform = .false. 
    endif 
     
    end function isWindowsPlatform 
    !------------------------------------------------------------------------------------ 
    
    
    !------------------------------------------------------------------------------------ 
    function Replace_Text (s, text, rep)  result(outs)
    
    implicit none
    
    character(*)            :: s, text, rep
    character(len(s)+100)   :: outs     ! provide outs with extra 100 char len
    integer                 :: i, nt, nr
    !...............................
    
    outs = s ; nt = len_trim(text) ; nr = len_trim(rep)
    do
        i = index(outs,text(:nt))
        if (i == 0) exit
        outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
    enddo
    
    end function Replace_Text
    !------------------------------------------------------------------------------------


    end module ModOScommands
