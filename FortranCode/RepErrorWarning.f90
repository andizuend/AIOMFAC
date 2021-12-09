!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Report error and warning messages to the error logfile from a list of defined      *
!*   cases.                                                                             *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
!*                                                                                      *
!*   -> created:        2011                                                            *
!*   -> latest changes: 2021-12-06                                                      *
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
SUBROUTINE RepErrorWarning(unito, errorflagmix, warningflag, errorflagcalc, pointi, errorind, warningind)

IMPLICIT NONE
!interface variables:
INTEGER(4),INTENT(IN) :: unito, errorflagmix, warningflag, errorflagcalc, pointi
INTEGER(4),INTENT(OUT) :: errorind, warningind
!...................................................................................

IF (errorflagmix /= 0) THEN !some mixture related error occured:
    SELECT CASE(errorflagmix)
    CASE(1)
        WRITE(unito,*) ""
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) "AIOMFAC ERROR 1: mixture related."
        WRITE(unito,*) "An organic main group <-> cation interaction parameter"
        WRITE(unito,*) "is not defined for the requested mixture. "
        WRITE(unito,*) "Please check your mixture components for available parameters"
        WRITE(unito,*) "stated in the AIOMFAC interaction matrix."
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) ""
    CASE(2)
        WRITE(unito,*) ""
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) "AIOMFAC ERROR 2: mixture related."
        WRITE(unito,*) "An organic main group <-> anion interaction parameter"
        WRITE(unito,*) "is not defined for the requested mixture. "
        WRITE(unito,*) "Please check your mixture components for available parameters"
        WRITE(unito,*) "stated in the AIOMFAC interaction matrix."
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) ""
    CASE(9)
        WRITE(unito,*) ""
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) "AIOMFAC ERROR 9: mixture related."
        WRITE(unito,*) "At least one cation <-> anion interaction parameter is "
        WRITE(unito,*) "not defined for the requested mixture. "
        WRITE(unito,*) "Please check all ion combinations for available parameters"
        WRITE(unito,*) "stated in the AIOMFAC interaction matrix."
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) ""
    CASE(13)
        WRITE(unito,*) ""
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) "AIOMFAC ERROR 13: Incorrect hydroxyl group assignment."
        WRITE(unito,*) "At least one component containing (CH_n[(OH)]) groups"
        WRITE(unito,*) "has been assigned an incorrect number of (OH) groups."
        WRITE(unito,*) "Note that the notation of a CHn group bonded to an OH"
        WRITE(unito,*) "group does not include the OH group; rather the hydroxyl"
        WRITE(unito,*) "groups have to be defined separately.                   "
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) ""
    CASE(14)
        WRITE(unito,*) ""
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) "AIOMFAC ERROR 14: Missing short-range ARR parameter."
        WRITE(unito,*) "A neutral main group <-> main group interaction coeff."
        WRITE(unito,*) "of this particular mixture is not available in the SR."
        WRITE(unito,*) "part of the model."
        WRITE(unito,*) "Check your organic components and their subgroups in"
        WRITE(unito,*) "comparison to available subgroups in the AIOMFAC matrix."
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) ""
    CASE(15)
        WRITE(unito,*) ""
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) "AIOMFAC ERROR 15: Missing short-range BRR parameter."
        WRITE(unito,*) "A neutral main group <-> main group interaction coeff."
        WRITE(unito,*) "of this particular mixture is not available in the SR."
        WRITE(unito,*) "part of the model for 3-parameter temperature dependence."
        WRITE(unito,*) "Check your organic components and their subgroups in"
        WRITE(unito,*) "comparison to available subgroups in the AIOMFAC matrix."
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) ""
    CASE DEFAULT
        WRITE(unito,*) ""
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) "AIOMFAC ERROR XX: an undefined mixture error occurred!"
        WRITE(unito,*) "errorflagmix = ", errorflagmix
        WRITE(unito,*) "======================================================="
        WRITE(unito,*) ""
    END SELECT
    errorind = errorflagmix
ELSE 
    !check warnings and errors related to occurences during specific data point calculations:
    IF (warningflag > 0) THEN
        SELECT CASE(warningflag)
        CASE(10)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC WARNING 10: Temperature range related."
            WRITE(unito,*) "At least one data point has a set temperature outside of"
            WRITE(unito,*) "the recommended range for model calculations of "
            WRITE(unito,*) "electrolyte-containing mixtures. This may be intended, "
            WRITE(unito,*) "but caution is advised as AIOMFAC is not designed to  "
            WRITE(unito,*) "perform well at this temperature."
            WRITE(unito,*) "Data point no.: ", pointi
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE(11)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC WARNING 11: Temperature range related."
            WRITE(unito,*) "At least one data point has a set temperature outside of"
            WRITE(unito,*) "the recommended range for model calculations of "
            WRITE(unito,*) "electrolyte-free organic mixtures. This may be intended,"
            WRITE(unito,*) "but caution is advised as AIOMFAC is not designed to "
            WRITE(unito,*) "perform well at this temperature."
            WRITE(unito,*) "Data point no.: ", pointi
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE(16)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC-VISC WARNING 16: Mixture viscosity issue.      "
            WRITE(unito,'(A)') "A problem occurred during the viscosity prediction, &
                &likely related to a missing pure-component viscosity value. &
                &Therefore, an unrealistic mixture viscosity of log_10(eta/[Pa.s]) = &
                &-9999.9999 is output."
            WRITE(unito,*) "Data point no.: ", pointi
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE DEFAULT
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC WARNING XX: an undefined WARNING occurred!"
            WRITE(unito,*) "warningflag = ", warningflag
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        END SELECT
        warningind = warningflag
    ENDIF !warningflag
    IF (errorflagcalc > 0) THEN
        SELECT CASE(errorflagcalc)
        CASE(3)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC ERROR 3: Mixture composition related."
            WRITE(unito,*) "Composition data for this point is missing or incorrect!"
            WRITE(unito,*) "The sum of the mole or mass fractions of all components"
            WRITE(unito,*) "has to be equal to 1.0 and individual mole or mass "
            WRITE(unito,*) "fractions have to be positive values <= 1.0!"
            WRITE(unito,*) "Composition point no.: ", pointi
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE(4,5)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC ERROR 4: Mixture composition related."
            WRITE(unito,*) "Composition data for this point is incorrect!"
            WRITE(unito,*) "The sum of the mole or mass fractions of all components"
            WRITE(unito,*) "has to be equal to 1.0 and individual mole or mass "
            WRITE(unito,*) "fractions have to be positive values <= 1.0!"
            WRITE(unito,*) "Composition point no.: ", pointi
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE(6,7)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC ERROR 6: Numerical issue."
            WRITE(unito,*) "A numerical issue occurred during computation of the"
            WRITE(unito,*) "data points flagged in the output tables."
            WRITE(unito,*) "This error was possibly caused due to input of very"
            WRITE(unito,*) "high electrolyte concentrations."
            WRITE(unito,*) "Composition point no.: ", pointi
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE(8)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC ERROR 8: Mixture composition related."
            WRITE(unito,*) "At least one neutral component must be present in the  "
            WRITE(unito,*) "system! "
            WRITE(unito,*) "Composition point no.: ", pointi
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE(12)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC ERROR 12: Charge neutrality violated."
            WRITE(unito,*) "The mixture violates the electrical charge neutrality  "
            WRITE(unito,*) "condition (moles cation*[cation charge] =              "
            WRITE(unito,*) "                          moles anion*[anion charge]). "
            WRITE(unito,*) "Make sure that selected integer amounts of cation and  "
            WRITE(unito,*) "anion 'subgroups' fulfill the charge balance (in the   "
            WRITE(unito,*) "inorganic component definition of the input file).     "
            WRITE(unito,*) "Composition point no.: ", pointi
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE(17)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,'(A)') "AIOMFAC ERROR 17: Issue with ion dissociation &
                            &equilibria calculations."
            WRITE(unito,'(A)') "The numerical solution of electrolyte/ion dissociation &
                            &equilibria was not accomplished to the desired tolerance &
                            &level. Model output for this point is unreliable and likely &
                            &incorrect."
            WRITE(unito,*) "Composition point no.: ", pointi
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE DEFAULT
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC ERROR XX: an undefined calculation error occurred!"
            WRITE(unito,*) "errorflagcalc = ", errorflagcalc
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        END SELECT
        errorind = errorflagcalc
    ENDIF !errorflagcalc
ENDIF !errorflagmix
        
END SUBROUTINE RepErrorWarning