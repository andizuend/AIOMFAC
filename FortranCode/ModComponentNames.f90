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
!*   -> latest changes: 2018/05/22                                                      *
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
CHARACTER(LEN=60),DIMENSION(1500),PUBLIC :: NKname, NKnameTeX !neutral component names (alphabetical names and DISLIN TeX-code version)
CHARACTER(LEN=30),DIMENSION(40,40),PUBLIC :: electname, electnameTeX !ion combinations component names

!================================================================================================================================= 
    CONTAINS
!================================================================================================================================= 
    
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
    !*   -> latest changes: 2016/10/05                                                      *
    !*                                                                                      *
    !****************************************************************************************

    SUBROUTINE nametab()

    IMPLICIT NONE
    !...........................................................

    !list of neutral component names:
    NKname = "not_defined"
    NKnameTeX = "not_defined"

    END SUBROUTINE nametab

!=================================================================================================================================
    
    
    !**********************************************************************************************************
    !*                                                                                                        *
    !*  Subroutine to set the character strings for the names of the different components in a mixture.       *
    !*                                                                                                        *
    !*                                                                                                        *   
    !*                                      (c) Andi Zuend, IACETH, ETH Zurich, 2009                          *
    !*                         Div. Chem. Engineering, California Institute of Technology, 2009 - 2012        *
    !**********************************************************************************************************
    PURE SUBROUTINE names_mix(ncomp, nneutral, CompN, compname, compnameTeX, ionname, ionnameTeX, OtoCratio, HtoCratio)

    USE ModSystemProp, ONLY : ElectComps, ElectSubs, Ncation, Nanion, NGS, frominpfile, cpname
    USE ModSubgroupProp, ONLY : O2C_H2C_component, subgrname, subgrnameTeX

    IMPLICIT NONE

    !interface vars:
    INTEGER(4),INTENT(IN) :: ncomp, nneutral  !the component number (internal numbering as the compN in definemixtures, LRdata, etc.)
    INTEGER(4),DIMENSION(ncomp),INTENT(IN) :: CompN
    REAL(8),DIMENSION(ncomp),INTENT(OUT) :: OtoCratio, HtoCratio
    CHARACTER(LEN=*),DIMENSION(ncomp),INTENT(OUT) :: compname, compnameTeX
    CHARACTER(LEN=*),DIMENSION(NGS),INTENT(OUT) :: ionname, ionnameTeX
    !local vars:
    CHARACTER(LEN=16) :: txt
    INTEGER(4) :: i, k, cnt, cn, an
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

    !electrolyte names are set in 'nametab':
    cnt = Ncation*Nanion
    DO i = 1,cnt
        cn = ElectComps(i,1)
        an = ElectComps(i,2)  
        cn = cn -200 !index shift for cations
        an = an -240 !index shift for anions
        compname(nneutral+i) = TRIM(electname(cn,an))
        compnameTeX(nneutral+i) = TRIM(electnameTeX(cn,an)) 
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
            txt = TRIM(subgrname(k))
            cn = VERIFY(txt, "( )", BACK = .false.) !first character index that is not ( or ) or a space
            an = VERIFY(txt, "( )", BACK = .true.) !last character index that is not ( or ) or a space
            !remove the parathensis (...) around the ion subgroup:
            ionname(cnt) = TRIM(txt(cn:an))
            txt = TRIM(subgrnameTeX(k))
            cn = VERIFY(txt, "( )", BACK = .false.) 
            an = VERIFY(txt, "( )", BACK = .true.) 
            ionnameTeX(cnt) = TRIM(txt(cn:an))
        ENDIF
    ENDDO
    !anions:
    DO i = 2,NGS !loop over ElectSubs
        k = ElectSubs(i)
        IF (k > 240) THEN !anion
            cnt = cnt +1
            txt = TRIM(subgrname(k))
            cn = VERIFY(txt, "( )", BACK = .false.) 
            an = VERIFY(txt, "( )", BACK = .true.) 
            !remove the parathensis (...) around the ion subgroup:
            ionname(cnt) = TRIM(txt(cn:an))
            txt = TRIM(subgrnameTeX(k))
            cn = VERIFY(txt, "( )", BACK = .false.) 
            an = VERIFY(txt, "( )", BACK = .true.)
            ionnameTeX(cnt) = TRIM(txt(cn:an))
        ENDIF
    ENDDO

    END SUBROUTINE names_mix
!================================================================================================================================= 
    
END MODULE ModComponentNames