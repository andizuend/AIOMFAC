!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Subroutine providing an interface to AIOMFAC for computing specific liquid phase   *
!*   activity coefficients and phase compositions with consideration of a single liquid * 
!*   phase as solution.                                                                 *  
!*   Input is a specific composition point and temperature for a previously initialized *
!*   system. Output is listed in outputvars.                                            *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andreas Zuend,                                                                     *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2011                                                            *
!*   -> latest changes: 2018/06/28                                                      *
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
SUBROUTINE AIOMFAC_inout(inputconc, xinputtype, TKelvin, nspecies, outputvars, outnames, errorflag, warningflag)

USE ModSystemProp
USE ModAIOMFACvar
USE ModCompScaleConversion
USE ModCalcActCoeff, ONLY : AIOMFAC_calc
USE ModSubgroupProp, ONLY : SMWA, SMWC

IMPLICIT NONE
!interface variables declarations 
REAL(8),DIMENSION(nindcomp),INTENT(IN) :: inputconc      !inputconc = the concentration of a given input point (e.g., at an experimental data point)
REAL(8),DIMENSION(6,NKNpNGS),INTENT(OUT) :: outputvars   !2-D output array with computed compositions and activities for each species; structure is:  | mass-frac., mole-frac., molality, act.coeff., activity, ion-indicator | species-no |
REAL(8),INTENT(IN) :: TKelvin                            !the input temperature [K]
LOGICAL(4),INTENT(IN) :: xinputtype
INTEGER(4),INTENT(OUT) :: nspecies, errorflag, warningflag
CHARACTER(LEN=60),DIMENSION(NKNpNGS) :: outnames
!--
!local variable declarations:
CHARACTER(LEN=2) :: cn  !this assumes a maximum two-digit component number in the system (max. 99); to be adjusted otherwise.
CHARACTER(LEN=3) :: cino
INTEGER(4) :: i, ion_no, ion_indic, nc, NKSinput, NKSinputp1
REAL(8),PARAMETER :: DEPS = 1.1D1*(EPSILON(DEPS))
REAL(8) :: wtf_cp, xi_cp, mi_cp, actcoeff_cp, a_cp, sum_ms, sum_mions, sum_miMi
REAL(8),DIMENSION(nelectrol) :: mixingratio, wtfdry
!------------------------------------------------------------------------------------------- 

!check for debugging:
IF (ANY(ISNAN(inputconc(1:nindcomp)))) THEN
    WRITE(*,*) "inputconc is NaN: ", inputconc(1:nindcomp)
    WRITE(*,*) ""
    RETURN
ENDIF
      
! Set initial values of some array variables:
errorflag = 0  
warningflag = 0
outputvars = 0.0D0
NKSinput = ninput -nneutral
NKSinputp1 = NKSinput+1
wtfdry(1:NKSinput) = 1.0D0
wtfdry(NKSinputp1:) = 0.0D0
mixingratio(1:NKSinput) = 1.0D0
mixingratio(NKSinputp1:) = 0.0D0
nspecies = NKNpNGS

IF (nneutral < 1) THEN  !leave the subroutine and indicate a problem to the calling routine
    errorflag = 8 !there must be at least one neutral component in the mixture!
    RETURN
ENDIF

!weight (mass) fractions of the data point:
CALL Inputconc_to_wtf(inputconc, mixingratio, wtfdry, xinputtype, wtf)  !convert input concentrations to wtf
IF (nneutral > 0 .AND. SUM(wtf(1:nneutral)) < DEPS) THEN  !leave the subroutine and indicate a problem to the calling routine
    errorflag = 3
    RETURN
ENDIF
IF (nindcomp > 0 .AND. ANY(wtf(1:nindcomp) < -DEPS)) THEN  !leave the subroutine and indicate a problem to the calling routine
    errorflag = 4
    RETURN
ENDIF
IF (ANY(ISNaN(wtf(1:nindcomp)))) THEN  !leave the subroutine and indicate a problem to the calling routine
    errorflag = 5
    RETURN
ENDIF
!.....
CALL MassFrac2MoleFracMolality(wtf, XrespSalt, mrespSalt)
CALL AIOMFAC_calc(wtf, TKelvin) !calculate at given mass fraction and temperature
!.....

!Output of the AIOMFAC calculated values species-wise to array outputvars (ions separately):
sum_ms = SUM(mrespSalt(1:nneutral))
sum_miMi = 0.0D0
DO i = 1,NGI !calculate sum of ion molalities, mi, times ion molar mass, Mi; [kg]
    IF (Ication(i) > 0) THEN
        sum_miMi = sum_miMi + smc(i)*SMWC(Ication(i)-200)*1.0D-3 
    ENDIF
    IF (Ianion(i) > 0) THEN
        sum_miMi = sum_miMi + sma(i)*SMWA(Ianion(i)-240)*1.0D-3 
    ENDIF
ENDDO
DO nc = 1,nspecies !loop over components
    IF (nc <= nneutral) THEN !distinguish between neutral components and inorg. ions
        WRITE(cn,'(I2.2)') nc
        outnames(nc) = "comp_no_"//cn
        wtf_cp = wtf(nc)
        xi_cp = X(nc)
        mi_cp = mrespSalt(nc) !molality in solvent mixture [mol/(kg solvent mix)]      
        actcoeff_cp = actcoeff_n(nc)
        IF (wtf(nc) > 0.0D0) THEN
            actcoeff_cp = actcoeff_n(nc) 
        ELSE
            actcoeff_cp = 0.0D0
        ENDIF
        a_cp = activity(nc)
        ion_indic = 0
    ELSE
        ion_no = ElectSubs(nc-nneutral) !the current ion subgroup number to identify the ion as output component nc
        WRITE(cino,'(I3.3)') ion_no
        outnames(nc) = "ion_no_"//cino
        ion_indic = ion_no
        !detect whether it is a cation or an anion:
        IF (ion_no > 239) THEN !anion
            i = AnNr(ion_no) !the number i anion (storage location in sma(i) etc.)
            xi_cp = sma(i)/(sum_ms + SumIonMolalities)
            mi_cp = sma(i)
            wtf_cp = sma(i)*SMWA(Ianion(i)-240)*1.0D-3/(1.0D0+sum_miMi)
            IF (sma(i) > 0.0D0) THEN
                actcoeff_cp = actcoeff_a(i) !molal activity coeff.
            ELSE
                actcoeff_cp = 0.0D0
            ENDIF
            IF (actcoeff_a(i) >= 0.0D0) THEN
               a_cp = actcoeff_a(i)*sma(i) !molal activity
            ELSE
               a_cp = -9999.999999D0 !indicate a numerical problem
               errorflag = 7 
            ENDIF
        ELSE !cation
            i = CatNr(ion_no) 
            xi_cp = smc(i)/(sum_ms + SumIonMolalities)
            mi_cp = smc(i)
            wtf_cp = smc(i)*SMWC(Ication(i)-200)*1.0D-3/(1.0D0+sum_miMi)    
            IF (smc(i) > 0.0D0) THEN
                actcoeff_cp = actcoeff_c(i) !molal activity coeff.
            ELSE
                actcoeff_cp = 0.0D0
            ENDIF
            IF (actcoeff_c(i) >= 0.0D0) THEN
               a_cp = actcoeff_c(i)*smc(i) !molal activity
            ELSE
               a_cp = -9999.999999D0 !indicate a numerical problem 
               errorflag = 7
            ENDIF
        ENDIF 
    ENDIF
    outputvars(1,nc) = wtf_cp            !mass fraction of component nc
    outputvars(2,nc) = xi_cp             !mole fraction with respect to electrolytes (salts) dissociated into individual ions
    outputvars(3,nc) = mi_cp             !molality of nc (with respect to dissociated ions)
    outputvars(4,nc) = actcoeff_cp       !activity coefficient of the component / species
    outputvars(5,nc) = a_cp              !activity of the component
    outputvars(6,nc) = REAL(ion_indic, KIND=8) !the indicator if this species is an inorg. ion or not: 0 = not ion, a number > 200 indicates the ion ID from the list
ENDDO

!check applicable temperature range and state a warning "errorflag" if violated:
!applicable range for electrolyte-containing mixtures (approx.): 288.0 to 309.0 K (298.15 +- 10 K); (strictly valid range would be 298.15 K only)
!applicable range for electrlyte-free mixtures (approx.): 280.0 to 400.0 K
IF (NGS > 0) THEN !electrolyte-containing
    IF (TKelvin > 309.0D0 .OR. TKelvin < 288.0D0) THEN !set warning flag
        IF (warningflag == 0) THEN !don't overwrite an existing warning if non-zero
            warningflag = 10
        ENDIF
    ENDIF
ELSE
    IF (TKelvin > 400.0D0 .OR. TKelvin < 280.0D0) THEN !set warning flag
        IF (warningflag == 0) THEN !don't overwrite an existing warning if non-zero
            warningflag = 11
        ENDIF
    ENDIF
ENDIF

END SUBROUTINE AIOMFAC_inout
! ======================= END =======================================================
