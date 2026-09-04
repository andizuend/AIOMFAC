!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module with subroutines and functions to carry out different operating system      *
!*   commands working both on Windows and Linux platforms. Main uses include copying of *
!*   files via commands executed in a terminal command line.                            *
!*                                                                                      *
!*   :: Authors & Copyright ::                                                          *
!*   Andi Zuend                                                                         *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2020                                                            *
!*   -> latest changes: 2026-09-03                                                      *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  subroutine copy_file                                                            *
!*   -  function   f_query_OS                                                           *
!*   -  function   f_replace_text                                                       *
!*                                                                                      *
!****************************************************************************************
module ModOScommands

implicit none
    
logical,public :: isWindowsOS
    
public :: copy_file, f_query_OS, f_replace_text

    contains

    !--------------------------------------------------------------------------------------
    !   This subroutine copies a file 'from_file' to 'to_file' via a system command.
    !   The file strings can include the (relative) path on the system.
    !--------------------------------------------------------------------------------------
    subroutine copy_file(from_file, to_file)

    !use ModStringFunctions

    implicit none

    character(len=*),intent(in) :: from_file, to_file
    !local:
    character(len = 100 + 2*max(len(from_file), len(to_file))) :: command, command2
    integer :: Estat, Cstat
    logical :: isWindowsOS
    !.................................

    isWindowsOS = f_query_OS()                                                       !function call
    
    if (isWindowsOS) then
        command = 'copy "'//trim(from_file) //'" "'//trim(to_file)//'" > NUL'
        command2 = trim(f_replace_text(command, "/", "\"))                                  !replace forward- by backslashes for Windows commands
        call execute_command_line(trim(command2), exitstat=Estat, cmdstat=Cstat)            !first copy file
        if (Cstat /= 0) write(*,*) "ERROR in copy_file 1: Cstat = ", Cstat
        
        command = 'icacls '//trim(to_file)//' /grant Users:M > NUL'                         !then set access permissions on Windows to permissive
        call execute_command_line(trim(command), exitstat=Estat, cmdstat=Cstat)
        if (Cstat /= 0) write(*,*) "ERROR in copy_file 2: Cstat = ", Cstat
        
    else !on a LINUX OS?
        command2 = 'cp -p"'//trim(from_file) //'" "'//trim(to_file)//'" > NUL'              !copy while keeping access permissions as of the original file
        call execute_command_line(trim(command2), exitstat=Estat, cmdstat=Cstat)
        if (Cstat /= 0) write(*,*) "ERROR in copy_file 1: Cstat = ", Cstat
    endif

    end subroutine copy_file
    !------------------------------------------------------------------------------------


    !------------------------------------------------------------------------------------
    !A logical querry function that returns 'true' if the operating system / platform is a version of Windows.
    function f_query_OS() result(isWindowsPlatform)

    implicit none

    logical             :: isWindowsPlatform        !output value
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
    else !it is likely Linux or we pretend it is
        !write(*,'(A)') "ERROR: Could not determine operating system in function isWindowsPlatform!"
        !read(*,*) !wait for user action
        isWindowsPlatform = .false.
    endif

    end function f_query_OS
    !------------------------------------------------------------------------------------
    
    
    !------------------------------------------------------------------------------------
    pure function f_replace_text(s, text, rep)  result(outs)
    
    implicit none
    !interface arguments:
    character(*),intent(in)       :: s          !text string to process
    character(*),intent(in)       :: text       !text characters to be searched for and replaced (repeatedly)
    character(*),intent(in)       :: rep        !replacement string of text
    character(len=:),allocatable  :: outs       !(output) the processed string
    !local variables:
    integer                 :: i, nt, nr
    !...............................

    allocate(character(len(s)) :: outs)
    outs = s
    nt = len_trim(text)
    nr = len_trim(rep)
    
    do
        i = index(outs, text(:nt))
        if (i == 0) exit
        outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
    enddo
    
    end function f_replace_text
    !------------------------------------------------------------------------------------


end module ModOScommands