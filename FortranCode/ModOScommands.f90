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
!*   -> latest changes: 2026-09-15                                                      *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  subroutine copy_file                                                            *
!*   -  function   f_query_OS                                                           *
!*   -  function   f_replace_text                                                       *
!*   -  function   f_epoch_time                                                         *
!*                                                                                      *
!****************************************************************************************
module ModOScommands

use Mod_kind_param, only : wp

implicit none
    
logical,public :: isWindowsOS
    
public :: copy_file, f_query_OS, f_replace_text, f_epoch_time

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
    
    
    !------------------------------------------------------------------------------------
    !-- Function to return the UNIX-style current epoch time in terms of the integer  
    !   number of seconds since January 1, 1970. The code interfaces via iso_c_binding   
    !   with the system's C library to provide a portable OS-independent implementation.
    function f_epoch_time()  result(time_sec)
    
    use, intrinsic :: iso_c_binding, only : c_int64_t
    
    implicit none
    !interface arguments:
    real(wp) :: time_sec                            ![s] (output) the time in seconds since January 1, 1970
    !local variables:
    integer(c_int64_t) :: epoch_seconds
    !...............................

    !interface block to bind to the standard C 'time' function
    interface
        function c_time(t) bind(C, name="time")
            import :: c_int64_t
            integer(c_int64_t), intent(in), value :: t
            integer(c_int64_t) :: c_time
        end function c_time
    end interface

    epoch_seconds = c_time(0_c_int64_t)             !passing 0_c_int64_t acts like passing NULL in C
    time_sec = real(epoch_seconds, kind=wp)

    end function f_epoch_time
    !------------------------------------------------------------------------------------


end module ModOScommands