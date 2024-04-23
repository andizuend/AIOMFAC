!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module to define the character strings for the names of the different components   *
!*   or species in AIOMFAC and to compute the component names present in a mixture.     * 
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2006                                                            *
!*   -> latest changes: 2021-10-01                                                      *
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
!*   -  subroutine nametab                                                              *
!*   -  subroutine names_mix                                                            *
!*                                                                                      *
!****************************************************************************************

module ModComponentNames

implicit none

!module public vars:
character(len=60),dimension(:),allocatable,public :: NKname, NKnameTeX

!========================================================================================================== 
    contains
!========================================================================================================== 
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to set the character strings for the names of the different independent *
    !*   components in a mixture.                                                           *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, 2009                                                           *
    !*   Dept. Chem. Engineering, California Institute of Technology (2009 - 2012),         *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2009                                                            *
    !*   -> latest changes: 2023-08-02                                                      *
    !*                                                                                      *
    !****************************************************************************************

    subroutine nametab()

    implicit none
    !.....................................................
    
    allocate( NKname(1500), NKnameTeX(1500) )

    !list of neutral component names:
    NKname = "not_defined"
    NKnameTeX = "not_defined"
    
    !...
    NKname(1500) = "!fromInpFile!";  NKnameTeX(1500) = "!fromInpFile!";   !this ID number is used to indicate that the name is set from input via a file.

    end subroutine nametab
!==========================================================================================================
    

    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to set the character strings for the names of the different components  *
    !*   and ions in a mixture.                                                             * 
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
    !*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2006                                                            *
    !*   -> latest changes: 2021/10/01                                                      *
    !*                                                                                      *
    !****************************************************************************************
    pure subroutine names_mix(CompN, compname, compnameTeX, ionname, ionnameTeX, OtoCratio, HtoCratio)

    use Mod_kind_param, only : wp
    use ModSystemProp, only : ElectComps, ElectNues, ElectSubs, nneutral, Ncation, Nanion, NGS, frominpfile, cpname
    use ModSubgroupProp, only : O2C_H2C_component, subgrname, subgrnameTeX
    use ModStringFunctions, only : Replace_Text, Replace_Text_Advance

    implicit none

    !interface vars:
    integer,dimension(:),intent(in) :: CompN
    real(wp),dimension(:),intent(out) :: OtoCratio, HtoCratio
    character(len=*),dimension(:),intent(out) :: compname, compnameTeX
    character(len=*),dimension(:),intent(out) :: ionname, ionnameTeX
    !local vars:
    character(len=1) :: nchar
    character(len=16) :: catxt, antxt
    character(len=40) :: txt
    integer :: i, k, nn, cnt, cn, an, nue_c, nue_a
    real(wp) :: compC, compH, compO
    !...........................................................

    OtoCratio = -7.777777_wp  !set to an impossible value at initialization
    HtoCratio = -7.777777_wp
    !Get the names of the actual components in mixture nd
    !neutral components:
    !water is always component 1 (if there is water in the mixture):
    if (CompN(1) == 401) then !there is water
        compname(1) = "Water"
        compnameTeX(1) = "Water"
    endif

    do i = 1,nneutral
        k = CompN(i)
        if (k > 0 .AND. k /= 401) then !scan all components, except water
            if (frominpfile) then
                compname(i) = cpname(i) 
                compnameTeX(i) = cpname(i) 
            else
                compname(i) = trim(NKname(k)) !//"'"    !the "'" is added to make sure that for comma separated reading the commas in e.g. 1,2-Butanediol is not separating the species name.
                compnameTeX(i) = trim(NKnameTeX(k)) !//"'"
            endif
            !call the OtoC calculation subroutine which uses atom information from subgroups;  here just for the single components "i".
            call O2C_H2C_component(i, compC, compH, compO, OtoCratio(i), HtoCratio(i))
        endif
    enddo

    !construct electrolyte names from ion subgroup names:
    cnt = Ncation*Nanion
    do i = 1,cnt
        cn = ElectComps(i,1)
        an = ElectComps(i,2)  
        !remove the parathensis ( ) around the ion subgroup as well as the charge characters:
        catxt = trim( Replace_Text(trim(subgrname(cn)), '(', '') ) 
        catxt = trim( Replace_Text(catxt, ')', '') )
        catxt = trim( Replace_Text(catxt, '+', '') )
        antxt = trim( Replace_Text(trim(subgrname(an)), '(', '') ) 
        antxt = trim( Replace_Text(antxt, ')', '') )
        antxt = trim( Replace_Text(antxt, '-', '') )
        !construct 'neutral' electrolyte component name from ion info using
        !common naming convention from chemistry:
        nue_c = ElectNues(i,1)
        nue_a = ElectNues(i,2)
        if (nue_c /= nue_a) then
            if (nue_c > nue_a) then
                write(nchar,'(I0)') nue_c/nue_a
                k = len_trim(catxt)
                nn = iachar(catxt(k:k))
                if (nn > 49 .AND. nn < 57) then     !the last character is the number
                    compname(nneutral+i) = '('//trim(catxt)//')'//nchar//trim(antxt)
                else
                    compname(nneutral+i) = trim(catxt)//nchar//trim(antxt)
                endif   
            else
                write(nchar,'(I0)') nue_a/nue_c
                k = len_trim(antxt)
                nn = iachar(antxt(k:k))
                if (nn > 49 .AND. nn < 57) then     !the last character is the number
                    compname(nneutral+i) = trim(catxt)//'('//trim(antxt)//')'//nchar
                else
                    compname(nneutral+i) = trim(catxt)//trim(antxt)//nchar
                endif 
            endif
        else   
            compname(nneutral+i) = trim(catxt)//trim(antxt)        
        endif
        !generate the name also using TeX formatting:
        txt = trim(compname(nneutral+i))
        txt = trim( Replace_Text_Advance(trim(txt), '2', '$_2$') ) 
        txt = trim( Replace_Text_Advance(trim(txt), '3', '$_3$') ) 
        txt = trim( Replace_Text_Advance(trim(txt), '4', '$_4$') ) 
        compnameTeX(nneutral+i) = trim(txt)
    enddo
    !set the rest of the names to an empty string
    compname(nneutral+cnt+1:) = ""
    compnameTeX(nneutral+cnt+1:) = ""

    !set also ion names, especially for output use with systems containing several electrolytes:
    !cations:
    cnt = 0
    do i = 1,NGS !loop over ElectSubs
        k = ElectSubs(i)
        if (k < 240) then !cation
            cnt = cnt +1
            !remove the parathensis (...) around the ion subgroup:
            txt = trim( Replace_Text(trim(subgrname(k)), '(', '') ) 
            txt = trim( Replace_Text(txt, ')', '') )
            ionname(cnt) = trim(txt)
            txt = trim( Replace_Text(trim(subgrnameTeX(k)), '(', '') ) 
            txt = trim( Replace_Text(txt, ')', '') )
            ionnameTeX(cnt) = trim(txt)
        endif
    enddo
    !anions:
    do i = 2,NGS !loop over ElectSubs
        k = ElectSubs(i)
        if (k > 240) then !anion
            cnt = cnt +1
            !remove the parathensis (...) around the ion subgroup:
            txt = trim( Replace_Text(trim(subgrname(k)), '(', '') ) 
            txt = trim( Replace_Text(txt, ')', '') )
            ionname(cnt) = trim(txt)
            txt = trim( Replace_Text(trim(subgrnameTeX(k)), '(', '') ) 
            txt = trim( Replace_Text(txt, ')', '') )
            ionnameTeX(cnt) = trim(txt)
        endif
    enddo

    end subroutine names_mix
!========================================================================================================== 
    
end module ModComponentNames