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
!*   -> latest changes: 2021-12-08                                                      *
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
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  SUBROUTINE nametab                                                              *
!*   -  SUBROUTINE names_mix                                                            *
!*                                                                                      *
!****************************************************************************************

MODULE ModComponentNames

IMPLICIT NONE

!module public vars:
CHARACTER(LEN=60),DIMENSION(:),ALLOCATABLE,PUBLIC :: NKname, NKnameTeX

!========================================================================================================== 
    CONTAINS
!========================================================================================================== 
    
    SUBROUTINE nametab()  !not used specifically in web version of AIOMFAC

    IMPLICIT NONE
    !..............................................
    
    ALLOCATE( NKname(1500), NKnameTeX(1500) )
    
    !list of neutral component names:
    NKname = "not_defined"
    NKnameTeX = "not_defined"

    END SUBROUTINE nametab
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
    PURE SUBROUTINE names_mix(CompN, compname, compnameTeX, ionname, ionnameTeX, OtoCratio, HtoCratio)

    USE ModSystemProp, ONLY : ElectComps, ElectNues, ElectSubs, nneutral, Ncation, Nanion, NGS, frominpfile, cpname
    USE ModSubgroupProp, ONLY : O2C_H2C_component, subgrname, subgrnameTeX
    USE ModStringFunctions, ONLY : Replace_Text, Replace_Text_Advance

    IMPLICIT NONE

    !interface vars:
    INTEGER(4),DIMENSION(:),INTENT(IN) :: CompN
    REAL(8),DIMENSION(:),INTENT(OUT) :: OtoCratio, HtoCratio
    CHARACTER(LEN=*),DIMENSION(:),INTENT(OUT) :: compname, compnameTeX
    CHARACTER(LEN=*),DIMENSION(:),INTENT(OUT) :: ionname, ionnameTeX
    !local vars:
    CHARACTER(LEN=1) :: nchar
    CHARACTER(LEN=16) :: catxt, antxt
    CHARACTER(LEN=40) :: txt
    INTEGER(4) :: i, k, nn, cnt, cn, an, nue_c, nue_a
    REAL(8) :: compC, compH, compO
    !...........................................................

    OtoCratio = -7.777777D0  !set to an impossible value at initialization
    HtoCratio = -7.777777D0
    !Get the names of the actual components in mixture nd
    !neutral components:
    !water is always component 1 (if there is water in the mixture):
    IF (CompN(1) == 401) THEN !there is water
        compname(1) = "Water"
        compnameTeX(1) = "Water"
    ENDIF

    DO i = 1,nneutral
        k = CompN(i)
        IF (k > 0 .AND. k /= 401) THEN !scan all components, except water
            IF (frominpfile) THEN
                compname(i) = cpname(i) 
                compnameTeX(i) = cpname(i) 
            ELSE
                compname(i) = TRIM(NKname(k)) !//"'"    !the "'" is added to make sure that for comma separated reading the commas in e.g. 1,2-Butanediol is not separating the species name.
                compnameTeX(i) = TRIM(NKnameTeX(k)) !//"'"
            ENDIF
            !call the OtoC calculation subroutine which uses atom information from subgroups;  here just for the single components "i".
            CALL O2C_H2C_component(i, compC, compH, compO, OtoCratio(i), HtoCratio(i))
        ENDIF
    ENDDO

    !construct electrolyte names from ion subgroup names:
    cnt = Ncation*Nanion
    DO i = 1,cnt
        cn = ElectComps(i,1)
        an = ElectComps(i,2)  
        !remove the parathensis ( ) around the ion subgroup as well as the charge characters:
        catxt = TRIM( Replace_Text(TRIM(subgrname(cn)), '(', '') ) 
        catxt = TRIM( Replace_Text(catxt, ')', '') )
        catxt = TRIM( Replace_Text(catxt, '+', '') )
        antxt = TRIM( Replace_Text(TRIM(subgrname(an)), '(', '') ) 
        antxt = TRIM( Replace_Text(antxt, ')', '') )
        antxt = TRIM( Replace_Text(antxt, '-', '') )
        !construct 'neutral' electrolyte component name from ion info using
        !common naming convention from chemistry:
        nue_c = ElectNues(i,1)
        nue_a = ElectNues(i,2)
        IF (nue_c /= nue_a) THEN
            IF (nue_c > nue_a) THEN
                WRITE(nchar,'(I0)') nue_c/nue_a
                k = LEN_TRIM(catxt)
                nn = IACHAR(catxt(k:k))
                IF (nn > 49 .AND. nn < 57) THEN     !the last character is the number
                    compname(nneutral+i) = '('//TRIM(catxt)//')'//nchar//TRIM(antxt)
                ELSE
                    compname(nneutral+i) = TRIM(catxt)//nchar//TRIM(antxt)
                ENDIF   
            ELSE
                WRITE(nchar,'(I0)') nue_a/nue_c
                k = LEN_TRIM(antxt)
                nn = IACHAR(antxt(k:k))
                IF (nn > 49 .AND. nn < 57) THEN     !the last character is the number
                    compname(nneutral+i) = TRIM(catxt)//'('//TRIM(antxt)//')'//nchar
                ELSE
                    compname(nneutral+i) = TRIM(catxt)//TRIM(antxt)//nchar
                ENDIF 
            ENDIF
        ELSE   
            compname(nneutral+i) = TRIM(catxt)//TRIM(antxt)        
        ENDIF
        !generate the name also using TeX formatting:
        txt = TRIM(compname(nneutral+i))
        txt = TRIM( Replace_Text_Advance(TRIM(txt), '2', '$_2$') ) 
        txt = TRIM( Replace_Text_Advance(TRIM(txt), '3', '$_3$') ) 
        txt = TRIM( Replace_Text_Advance(TRIM(txt), '4', '$_4$') ) 
        compnameTeX(nneutral+i) = TRIM(txt)
    ENDDO
    !set the rest of the names to an empty string
    compname(nneutral+cnt+1:) = ""
    compnameTeX(nneutral+cnt+1:) = ""

    !set also ion names, especially for output use with systems containing several electrolytes:
    !cations:
    cnt = 0
    DO i = 1,NGS !loop over ElectSubs
        k = ElectSubs(i)
        IF (k < 240) THEN !cation
            cnt = cnt +1
            !remove the parathensis (...) around the ion subgroup:
            txt = TRIM( Replace_Text(TRIM(subgrname(k)), '(', '') ) 
            txt = TRIM( Replace_Text(txt, ')', '') )
            ionname(cnt) = TRIM(txt)
            txt = TRIM( Replace_Text(TRIM(subgrnameTeX(k)), '(', '') ) 
            txt = TRIM( Replace_Text(txt, ')', '') )
            ionnameTeX(cnt) = TRIM(txt)
        ENDIF
    ENDDO
    !anions:
    DO i = 2,NGS !loop over ElectSubs
        k = ElectSubs(i)
        IF (k > 240) THEN !anion
            cnt = cnt +1
            !remove the parathensis (...) around the ion subgroup:
            txt = TRIM( Replace_Text(TRIM(subgrname(k)), '(', '') ) 
            txt = TRIM( Replace_Text(txt, ')', '') )
            ionname(cnt) = TRIM(txt)
            txt = TRIM( Replace_Text(TRIM(subgrnameTeX(k)), '(', '') ) 
            txt = TRIM( Replace_Text(txt, ')', '') )
            ionnameTeX(cnt) = TRIM(txt)
        ENDIF
    ENDDO

    END SUBROUTINE names_mix
!========================================================================================================== 
    
END MODULE ModComponentNames