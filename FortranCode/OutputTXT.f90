!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Create an output ASCII text file with an overall mixture data header and           *
!*   individual tables for all components / species (in case of ions).                  *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
!*                                                                                      *
!*   -> created:        2011                                                            *
!*   -> latest changes: 2018/08/10                                                      *
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
SUBROUTINE OutputTXT(fname, VersionNo, cpnameinp, nspecmax, npoints, watercompno, T_K, px, out_data)

!module variables:
USE ModSystemProp, ONLY : compname, compsubgroups, compsubgroupsTeX, NGS, NKNpNGS, ninput, nneutral
USE ModSubgroupProp, ONLY : subgrname, subgrnameTeX

IMPLICIT NONE
!interface variables:
CHARACTER(LEN=*),INTENT(IN) :: fname 
CHARACTER(LEN=*),INTENT(IN) :: VersionNo
CHARACTER(LEN=60),DIMENSION(nspecmax),INTENT(IN) :: cpnameinp   !list of assigned component names (from input file)
INTEGER(4),INTENT(IN) ::  nspecmax, npoints, watercompno
REAL(8),DIMENSION(npoints),INTENT(IN) :: T_K
INTEGER(4),DIMENSION(nspecmax),INTENT(INOUT) :: px
REAL(8),DIMENSION(7,npoints,NKNpNGS),INTENT(IN) :: out_data
!--
!local variables:
CHARACTER(LEN=2) :: cn
CHARACTER(LEN=5) :: tlen
CHARACTER(LEN=50) :: subntxt
CHARACTER(LEN=150) :: horizline, tablehead, tformat, txtn
CHARACTER(LEN=3000) :: txtsubs, mixturestring
INTEGER(4) ::  i, k, pointi, qty, unitx
REAL(8) :: RH
!...................................................................................

!create a character string of the mixture as a series of its components:
mixturestring = TRIM(cpnameinp(1)) !first component
!loop over all further components / species:
DO i = 2,nspecmax
    px(i) = MAX(1, px(i))
    IF (INT(out_data(6,px(i),i)) == 0) THEN !neutral component
        txtn = TRIM(ADJUSTL(cpnameinp(i)))
    ELSE !ion (with its own link)
        IF (INT(out_data(6,px(i),i)) > 1) THEN
            txtn = TRIM( ADJUSTL(subgrname(INT(out_data(6,px(i),i) ))) )
        ELSE
            txtn = "unknown_sub"
        ENDIF
    ENDIF
    mixturestring = TRIM(mixturestring)//' + '//TRIM(ADJUSTL(txtn))
ENDDO
mixturestring = ADJUSTL(mixturestring)

OPEN (NEWUNIT = unitx, FILE = fname, STATUS = "UNKNOWN")
WRITE(unitx,'(A107)') "==========================================================================================================="
WRITE(unitx,'(A25)') "AIOMFAC-web, version "//VersionNo
WRITE(unitx,'(A107)') "==========================================================================================================="
WRITE(unitx,*) ""
mixturestring = "'"//TRIM(mixturestring)//"'"
![Note that Intel Fortran allows for variable format statements like A<LEN_TRIM(txtn)>, but that is a non-standard extension not supported by gfortran and other compilers - therefore avoided here.]
WRITE(tlen,'(I5.5)') LEN_TRIM(mixturestring)
tformat = '(A15, A'//tlen//')'  !dynamic format specifier
WRITE(unitx, tformat) "Mixture name:  ", ADJUSTL(mixturestring)
WRITE(unitx,'(A40, I2.2)') "Number of independent input components: ", ninput
WRITE(unitx,'(A40, I2.2)') "Number of different neutral components: ", nneutral
WRITE(unitx,'(A40, I2.2)') "Number of different inorganic ions    : ", NGS
WRITE(unitx,*) ""
WRITE(unitx,'(A78)') "The AIOMFAC output is tabulated for each component/species individually below."
WRITE(unitx,*) ""
WRITE(unitx,'(A107)') '---  Table key  -------------------------------------------------------------------------------------------'
WRITE(unitx,'(A100)') 'no.              : composition point number;                                                        '
WRITE(unitx,'(A100)') 'T [K]            : absolute temperature;                                                            '
WRITE(unitx,'(A100)') 'RH [%]           : relative humidity in equilibrium with the liquid mixture (bulk conditions);      '
WRITE(unitx,'(A100)') 'w(j) [-]         : weight fraction (mass fraction) of species "j";                                  '
WRITE(unitx,'(A100)') 'x_i(j) [-]       : mole fraction of species "j", calculated on the basis of completely              '
WRITE(unitx,'(A100)') '                   dissociated inorganic ions; exception: the partial dissociation of bisulfate     '
WRITE(unitx,'(A100)') '                   (HSO4- <--> H+ + SO4--) is explicitly considered when present in the mixture;    '
WRITE(unitx,'(A100)') 'm_i(j) [mol/kg]  : molality of species "j" [mol(j)/(kg solvent mixture)], where "solvent mixture"   '
WRITE(unitx,'(A100)') '                   refers to the electrolyte-free mixture (water + organics);                       '
WRITE(unitx,'(A100)') 'a_coeff_x(j) [-] : activity coefficient of "j", defined on mole fraction basis (used for non-ionic  '
WRITE(unitx,'(A100)') '                   components) with pure (liquid) component "j" reference state;                    '
WRITE(unitx,'(A100)') 'a_coeff_m(j) [-] : activity coefficient of "j", defined on molality basis (used for inorg. ions)    '
WRITE(unitx,'(A100)') '                   with reference state of infinite dilution of "j" in pure water;                  '
WRITE(unitx,'(A100)') 'a_x(j) [-]       : activity of "j", defined on mole fraction basis (pure component reference state);'
WRITE(unitx,'(A100)') 'a_m(j) [-]       : activity of "j", defined on molality basis (used for inorg. ions) with reference '
WRITE(unitx,'(A100)') '                   state of infinite dilution of "j" in pure water;                                 '
WRITE(unitx,'(A100)') 'flag             : error/warning flag, a non-zero value (error/warning number) indicates that a     '
WRITE(unitx,'(A100)') '                   numerical issue or a warning occurred at this data point, suggesting evaluation  '
WRITE(unitx,'(A100)') '                   with caution (warnings) or exclusion (errors) of this data point.                '
WRITE(unitx,'(A107)') '-----------------------------------------------------------------------------------------------------------'
WRITE(unitx,*) ''
WRITE(unitx,*) ''
!--
horizline = "-----------------------------------------------------------------------------------------------------------"
!write individual data tables for each component / ionic species:
DO i = 1,nspecmax
    WRITE(cn,'(I2.2)') i !component / species number as character string
    WRITE(unitx,*) ""
    !distinguish between neutral components and ionic species:
    IF (INT(out_data(6,px(i),i)) == 0) THEN !neutral component
        WRITE(unitx,'(A24,I2.2)') "Mixture's component # : ", i
        txtn = "'"//TRIM(ADJUSTL(compname(i)))//"'"
        WRITE(tlen,'(I5.5)') LEN_TRIM(txtn)
        tformat = '(A24, A'//tlen//')'
        WRITE(unitx, tformat)  "Component's name      : ", ADJUSTL(txtn)
        txtsubs = "'"//TRIM(compsubgroups(i))//"'"
        WRITE(tlen,'(I5.5)') LEN_TRIM(txtsubs)
        tformat = '(A24, A'//tlen//')'
        WRITE(unitx, tformat) "Component's subgroups : ", ADJUSTL(txtsubs)
        txtsubs = "'"//TRIM(compsubgroupsTeX(i))//"'"
        WRITE(tlen,'(I5.5)') LEN_TRIM(txtsubs)
        tformat = '(A24, A'//tlen//')'
        WRITE(unitx, tformat) "Subgroups, TeX format : ", ADJUSTL(txtsubs)
        !write table column headers:
        WRITE(unitx,'(A107)') ADJUSTL(horizline)
        tablehead = "no.   T_[K]     RH_[%]   w("//cn//")          x_i("//cn//")        m_i("//cn//")        a_coeff_x("//cn//")   a_x("//cn//")        flag "
        WRITE(unitx,'(2X, A105)') ADJUSTL(tablehead)
        !--
    ELSE IF (INT(out_data(6,px(i),i)) < 240) THEN !cation
        WRITE(unitx,'(A24,I2.2)') "Mixture's species, #  : ", i
        subntxt = TRIM(ADJUSTL(subgrname(INT(out_data(6,px(i),i)))))
        qty = LEN_TRIM(subntxt)
        txtn = ADJUSTL(subntxt(2:qty-1)) !to print the ion name without the enclosing parathesis ()
        txtn = "'"//TRIM(txtn)//"'"
        WRITE(tlen,'(I5.5)') LEN_TRIM(txtn)
        tformat = '(A24, A'//tlen//')'
        WRITE(unitx, tformat)  "Cation's name         : ", ADJUSTL(txtn)
        txtsubs = "'"//TRIM(subntxt)//"'"
        WRITE(tlen,'(I5.5)') LEN_TRIM(txtsubs)
        tformat = '(A24, A'//tlen//')'
        WRITE(unitx,tformat) "Cation's subgroups    : ", ADJUSTL(txtsubs)
        subntxt = TRIM(ADJUSTL(subgrnameTeX(INT(out_data(6,px(i),i)))))
        txtsubs = "'"//TRIM(subntxt)//"'"
        WRITE(tlen,'(I5.5)') LEN_TRIM(txtsubs)
        tformat = '(A24, A'//tlen//')'
        WRITE(unitx, tformat) "Subgroups, TeX format : ", ADJUSTL(txtsubs)
        !write table column headers:
        WRITE(unitx, '(A107)') ADJUSTL(horizline)
        tablehead = "no.   T_[K]     RH_[%]   w("//cn//")          x_i("//cn//")        m_i("//cn//")        a_coeff_m("//cn//")   a_m("//cn//")        flag "
        WRITE(unitx, '(2X,A105)') ADJUSTL(tablehead)
        !--
    ELSE IF (INT(out_data(6,px(i),i)) > 240) THEN !anion
        WRITE(unitx,'(A24,I2.2)') "Mixture's species, #  : ", i
        subntxt = TRIM( ADJUSTL( subgrname(INT(out_data(6,px(i),i))) ) )
        qty = LEN_TRIM(subntxt)
        txtn = ADJUSTL(subntxt(2:qty-1))
        txtn = "'"//TRIM(txtn)//"'"
        WRITE(tlen,'(I5.5)') LEN_TRIM(txtn)
        tformat = '(A24, A'//tlen//')'
        WRITE(unitx, tformat)  "Anion's name          : ", ADJUSTL(txtn)
        txtsubs = "'"//TRIM(subntxt)//"'"
        WRITE(tlen,'(I5.5)') LEN_TRIM(txtsubs)
        tformat = '(A24, A'//tlen//')'
        WRITE(unitx, tformat) "Anion's subgroups     : ", ADJUSTL(txtsubs)
        subntxt = TRIM( ADJUSTL( subgrnameTeX(INT(out_data(6,px(i),i))) ) )
        txtsubs = "'"//TRIM(subntxt)//"'"
        WRITE(tlen,'(I5.5)') LEN_TRIM(txtsubs)
        tformat = '(A24, A'//tlen//')'
        WRITE(unitx, tformat) "Subgroups, TeX format : ", ADJUSTL(txtsubs)
        !write table column headers:
        WRITE(unitx,'(A107)') ADJUSTL(horizline)
        tablehead = "no.   T_[K]     RH_[%]   w("//cn//")          x_i("//cn//")        m_i("//cn//")        a_coeff_m("//cn//")   a_m("//cn//")        flag "
        WRITE(unitx,'(2X,A105)') ADJUSTL(tablehead)
        !--
    ELSE
        !error
    ENDIF
    !--
    WRITE(unitx,'(A107)') ADJUSTL(horizline)
    !write data to table:
    DO pointi = 1,npoints  !loop over composition points
        IF (watercompno > 0) THEN
            RH = out_data(5,pointi,watercompno)*100.0D0 !RH in %
            IF (RH > 1000.0D0 .OR. RH < 0.0D0) THEN
                RH = -99.99D0
            ENDIF
        ELSE
            RH = 0.0D0
        ENDIF
        WRITE(unitx,'(I5.3,2X,F7.2,2X,F7.2,3X,5(ES12.5,3X),3X,I2)') pointi, T_K(pointi), RH, out_data(1:5,pointi,i), INT(out_data(7,pointi,i))
    ENDDO !pointi
    WRITE(unitx,'(A107)') ADJUSTL(horizline)
    WRITE(unitx,*) ""
ENDDO
WRITE(unitx,*) ""
WRITE(unitx,'(A107)') "==========================================================================================================="
CLOSE(unitx)
        
END SUBROUTINE OutputTXT