!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module containing several composition-dependent AIOMFAC variables of the present   * 
!*   mixture.                                                                           *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2009 (based on non-module version from 2005)                    *
!*   -> latest changes: 2018/05/29                                                      *
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
!*   -  SUBROUTINE AllocModAIOMFACcore                                                  *
!*                                                                                      *
!****************************************************************************************
MODULE ModAIOMFACvar

IMPLICIT NONE
!Public module variables:
REAL(8),PUBLIC :: alphaHSO4, deltaetamix, diffKHSO4, etamix, ionicstrength, lastTK, meanSolventMW, &
    & SumIonMolalities, T_K, Tmolal, TmolalSolvMix, Xwdissoc
REAL(8),DIMENSION(:),ALLOCATABLE,PUBLIC :: actcoeff_a, actcoeff_c, actcoeff_n, activity, eta0, lneta_cpn, &
    & fragil, galrln, gamrln, gasrln, gclrln, gcmrln, gcsrln, gnlrln, gnmrln, gnsrln, ionactivityprod, lnactcoeff_a, &
    & lnactcoeff_c, lnactcoeff_n, lnmeanmactcoeff, meanmolalactcoeff, mrespSalt, SMA, SMC, solvmixcorrMRa, &
    & solvmixcorrMRc, Tglass0, wtf, X, XN, XrespSalt
REAL(8),DIMENSION(201:261),PUBLIC :: actcoeff_ion, molality_ion
LOGICAL(4),PUBLIC :: DebyeHrefresh 
!..................................................

!make all variables of this module threadprivate for use in parallel execution with openMP:
!$OMP THREADPRIVATE( alphaHSO4, deltaetamix, diffKHSO4, etamix, lastTK, meanSolventMW, SumIonMolalities, &
    !$OMP & T_K, Xwdissoc, actcoeff_a, actcoeff_c, actcoeff_n, activity, eta0, lneta_cpn, fragil, galrln, gamrln, gasrln, &
    !$OMP & gclrln, gcmrln, gcsrln, gnlrln, gnmrln, gnsrln, ionactivityprod, ionicstrength, lnactcoeff_a, &
    !$OMP & lnactcoeff_c, lnactcoeff_n, lnmeanmactcoeff, meanmolalactcoeff, mrespSalt, SMA, SMC, &
    !$OMP & solvmixcorrMRa, solvmixcorrMRc, Tglass0, Tmolal, TmolalSolvMix, wtf, X, XN, XrespSalt, DebyeHrefresh )

!==========================================================================================================================
    CONTAINS
!==========================================================================================================================
    
    !utility subroutine to allocate/deallocate module variables after mixture parameters are known (from definemixtures).
    SUBROUTINE AllocModAIOMFACvar()
    
    USE ModSystemProp, ONLY : NGI, nindcomp, NKNpNGS, nneutral, nelectrol
    
    IMPLICIT NONE
    
    !-- allocate several composition-dependent variables:
    IF (ALLOCATED(wtf)) THEN
        DEALLOCATE ( sma, smc, wtf, X, XN, solvmixcorrMRc, solvmixcorrMRa, XrespSalt, mrespSalt, &
        & activity, meanmolalactcoeff, actcoeff_n, actcoeff_c, actcoeff_a, ionactivityprod, eta0, lneta_cpn, & 
        & fragil, gnlrln, gclrln, galrln, gnmrln, gcmrln, gamrln, gnsrln, gcsrln, gasrln, lnactcoeff_n, &
        & lnactcoeff_c, lnactcoeff_a, lnmeanmactcoeff, Tglass0 )
    ENDIF
    ALLOCATE( sma(NGI), smc(NGI), wtf(nindcomp), X(NKNpNGS), XN(NKNpNGS), XrespSalt(nindcomp), &
        & mrespSalt(nindcomp), activity(nindcomp), meanmolalactcoeff(nelectrol), actcoeff_n(nneutral), actcoeff_c(NGI), &
        & actcoeff_a(NGI), ionactivityprod(nelectrol), solvmixcorrMRc(NGI), solvmixcorrMRa(NGI), eta0(NKNpNGS), fragil(NKNpNGS), & 
        & lneta_cpn(NKNpNGS), Tglass0(NKNpNGS), & 
        & lnactcoeff_n(nneutral), gnlrln(nneutral), gnmrln(nneutral), gnsrln(nneutral), &
        & lnactcoeff_c(NGI), gclrln(NGI), gcmrln(NGI), gcsrln(NGI), lnmeanmactcoeff(nelectrol), &
        & lnactcoeff_a(NGI), galrln(NGI), gamrln(NGI), gasrln(NGI) )
    
    END SUBROUTINE AllocModAIOMFACvar
!==========================================================================================================================

END MODULE ModAIOMFACvar