!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module for records of SMILES-based pure-component properties.                      *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andreas Zuend, Joshua Penikis,                                                     *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2026-07-16                                                      *
!*   -> latest changes: 2026-08-24                                                      *
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
!*   -  subroutine load_purecomp_table                                                  *
!*   -  subroutine append_purecomp_entry                                                *
!*   -  function find_smiles_index                                                      *
!*   -  subroutine lookup_Tg                                                            *
!*                                                                                      *
!****************************************************************************************
module ModPureCompProp

use Mod_kind_param, only : wp
use ModSystemProp, only : maxsmileslength

implicit none
!module variables
character(len=:),allocatable,public :: PureCompFile, RelPathModel

!SMILES-based pure-component reference table type:
type :: purecomp_entry
    character(len=maxsmileslength) :: smiles
    real(wp) :: Tg
    !future: add variables holding other pure component properties
end type purecomp_entry

type(purecomp_entry),dimension(:),allocatable :: purecomp_table
!.....................................................

!public procedures
public :: load_purecomp_table, append_purecomp_entry, lookup_Tg

private     !as default

!============================================================================================
    contains
!============================================================================================

    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Load table containing SMILES and their associated pure-component properties        *
    !*   from CSV file located in Auxillary folder                                          *
    !*                                                                                      *
    !*   :: Authors & Copyright ::                                                          *
    !*   Joshua Penikis, Andi Zuend,                                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *                       
    !*                                                                                      *
    !*   -> created:        2026-07-16                                                      *
    !*   -> latest changes: 2026-08-24                                                      *      
    !*                                                                                      *
    !****************************************************************************************
    subroutine load_purecomp_table()
    
    use ModSystemProp, only : errorflag_clist
    
    implicit none
    !local variables
    character(len=20) :: dummy
    integer :: unitx, istat, nlines, i
    logical :: fileexists
    !.....................................
    
    !initialize variables
    if (allocated(purecomp_table)) deallocate(purecomp_table)
    
    !check if file exists and read its content if true:
    do i = 1,4
        select case(i)
        case(1)  
            RelPathModel = '../AIOMFAC-model/'
        case(2)
            RelPathModel = './'
        case(3)
            RelPathModel = '.././'
        case default
            write(*,'(A)') 'ERROR in load_purecomp_table: could not find file Pure_component_smiles_table.csv!'
        end select
        PureCompFile = RelPathModel//'Auxiliary/Pure_component_smiles_table.csv'
        inquire(file = PureCompFile, exist = fileexists)
        if (fileexists) then
            exit
        endif
    enddo
    
    nlines = 0
    if (fileexists) then                    !open file
        open (newunit = unitx, file = PureCompFile, iostat = istat, action = 'read', status = 'unknown')
        if (istat /= 0) then                !an error occurred
            errorflag_clist(19) = .true.    !error opening pure-component property file
        endif    
        !read file
        read(unitx,*) dummy                 !read SMILES,Tg header
        do
            read(unitx,*,iostat=istat) dummy, dummy
            if (istat /= 0 ) exit           !file end reached
            nlines = nlines + 1             !count smiles rows
        enddo
        rewind(unitx)                       !read file from beginning
        read(unitx,*)                       !skip past header
        allocate(purecomp_table(nlines))    !set amount of entries
        do i = 1, nlines
            read(unitx,*) purecomp_table(i)%smiles, purecomp_table(i)%Tg
        enddo
        close(unitx)
    else                                    !pure-component property file does not exist
        errorflag_clist(22) = .true.
    endif
    
    end subroutine load_purecomp_table
    !------------------------------------------------------------------------------------
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Append new pure-component property value(s) for respective SMILES to table and     *
    !*   CSV file for record keeping.                                                       *
    !*                                                                                      *
    !*   :: Authors & Copyright ::                                                          *
    !*   Joshua Penikis, Andi Zuend,                                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *                       
    !*                                                                                      *
    !*   -> created:        2026-07-16                                                      *
    !*   -> latest changes: 2026-07-23                                                      *      
    !*                                                                                      *
    !****************************************************************************************
    subroutine append_purecomp_entry(smiles_new, Tg_value)
    
    implicit none
    !local variables
    character(len=*),intent(in) :: smiles_new
    real(wp),intent(in) :: Tg_value
    type(purecomp_entry),dimension(:),allocatable :: temp     !temporary entry array
    integer :: unitx, n
    !.....................................
    
    !initialize variables
    n =  size(purecomp_table)               !current row count
    
    !append to custom type
    allocate(temp(n + 1))
    temp(1:n) = purecomp_table(1:n)
    temp(n + 1)%smiles = smiles_new         !append smiles
    temp(n + 1)%Tg = Tg_value               !append Tg value
    !future: append other PC properties
    call move_alloc(temp, purecomp_table)
    !---
    !append to CSV file
    open (newunit=unitx, file=PureCompFile, position='append', action='write', status='old')
    write(unitx, '("""", A, """", ",", F0.2)') trim(smiles_new), Tg_value    !smiles saved in quotes, Tg saved to 2 decimal places
    close(unitx)
    
    end subroutine append_purecomp_entry
    !------------------------------------------------------------------------------------
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Identify index of SMILES in pure-component properties table.                       *
    !*                                                                                      *
    !*   :: Authors & Copyright ::                                                          *
    !*   Joshua Penikis, Andi Zuend,                                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *                       
    !*                                                                                      *
    !*   -> created:        2026-07-16                                                      *
    !*   -> latest changes: 2026-08-23                                                      *      
    !*                                                                                      *
    !****************************************************************************************
    pure elemental function find_smiles_index(smiles_string)  result(ind)
    
    implicit none
    !local variables
    character(len=*),intent(in) :: smiles_string
    integer :: ind
    !.....................................
    
    !search for smiles match
    ind = findloc(purecomp_table(:)%smiles, value=smiles_string, dim=1)
    if (ind == 0) then
        !no match found
    endif
        
    end function find_smiles_index
    !------------------------------------------------------------------------------------
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Using SMILES, identify if there exists an entry containing its respective Tg       *
    !*   value. If found, return Tg value.                                                  *
    !*                                                                                      *
    !*   :: Authors & Copyright ::                                                          *
    !*   Joshua Penikis, Andi Zuend,                                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *                       
    !*                                                                                      *
    !*   -> created:        2026-07-16                                                      *
    !*   -> latest changes: 2026-07-23                                                      *      
    !*                                                                                      *
    !****************************************************************************************
    pure subroutine lookup_Tg(smiles_string, Tg_found, found)
    
    implicit none
    !local variables
    character(len=*),intent(in) :: smiles_string
    real(wp),intent(out) :: Tg_found
    logical,intent(out) :: found
    integer :: ind
    !..............................
    
    ind = find_smiles_index(trim(smiles_string))    !search for smiles in custom type
    
    found = (ind > 0)                               !set to .true. if ind not 0
    if (found) then
        Tg_found = purecomp_table(ind)%Tg
    else
        Tg_found = -95.0_wp                         !return an unrealistic value
    endif
    
    end subroutine lookup_Tg
    !------------------------------------------------------------------------------------
    
    
end module ModPureCompProp
