!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module containing the main calculation subroutines of the SR (modified UNIFAC)     * 
!*   part of the AIOMFAC model.                                                         *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2018 (based on non-module version from 2004)                    *
!*   -> latest changes: 2021-09-17                                                      *
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
!*   -  SUBROUTINE SRunifac                                                             *
!*   -  SUBROUTINE SRcalcvisc                                                           *
!*   -  SUBROUTINE SRgref                                                               *
!*   -  SUBROUTINE SRgres                                                               *
!*   -  SUBROUTINE SRgcomb                                                              * 
!*   -  SUBROUTINE SRsystm                                                              *
!*                                                                                      *
!****************************************************************************************
MODULE ModSRunifac

USE ModAIOMFACvar, ONLY : lastTK, ln_eta0
USE ModSystemProp, ONLY : nd, nneutral, NG, NGN, Allsubs, ITABsr, Nmaingroups, topsubno, &
    & maingrindexofsubgr, Imaingroup, isPEGsystem, calcviscosity, solvmixrefnd, COMPN
USE ModSRparam, ONLY : ARR, BRR, CRR, SR_RR, SR_QQ

IMPLICIT NONE
!private module variables
INTEGER(4),PRIVATE :: grefcallID
INTEGER(4),DIMENSION(:,:),ALLOCATABLE,PRIVATE :: SRNY, SRNY_dimflip
REAL(8),DIMENSION(:),ALLOCATABLE,PRIVATE :: RS, QS, R, Q, XL, lnGaCinf, lnGaRref, XieRref
REAL(8),DIMENSION(:,:),ALLOCATABLE,PRIVATE :: parA, parB, parC, PsiT, PsiT_dimflip
REAL(8),DIMENSION(:,:),ALLOCATABLE,PRIVATE :: Nvis
LOGICAL(4),PRIVATE :: gcombrefresh

!$OMP THREADPRIVATE(parA, parB, parC, PsiT, PsiT_dimflip, grefcallID, SRNY, SRNY_dimflip, R, Q, RS, QS, XL, Nvis, &
   !$OMP & lnGaCinf, lnGaRref, XieRref, gcombrefresh)

!==========================================================================================================================
    CONTAINS
!==========================================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Short range SR (modified UNIFAC) part of the AIOMFAC model.                        *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
    !*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA                *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2004                                                            *
    !*   -> latest changes: 2019/09/17                                                      *
    !*                                                                                      *
    !****************************************************************************************
    
    !   lnGaC is the natural log of the combinatorial activity coefficient.
    !   lnGaRref is the natural log of the reference residual activity coefficient (pure component reference).
    !   lnGaR is the natural log of the residual activity coefficient of the components in the mixture.
    !   lnGaSR is the natural log of the total short-range activity coefficient.
    !   ITABsr lists the components and corresponding subgroups, sorted such that cations and anions are individual components (for SR part)
    !   XieC and XieR are the non-electrolyte component viscosity contributions of the combinatorial and residual parts, respectively.
    
    SUBROUTINE SRunifac(NK, T_K, X, XN, refreshgref, lnGaSR)

    USE ModPureCompViscos, ONLY : PureCompViscosity

    IMPLICIT NONE

    INTEGER(4),INTENT(IN) :: NK
    REAL(8),INTENT(IN) :: T_K
    REAL(8),DIMENSION(:),INTENT(IN) :: X, XN
    LOGICAL,INTENT(IN) :: refreshgref 
    REAL(8),DIMENSION(:),INTENT(OUT) :: lnGaSR
    !local variables:
    REAL(8),DIMENSION(NK) :: lnGaC, lnGaR, XieC, XieR
    !.......................................................
    
    !Determine the reference values (lnGaRref, XieRref) for the residual part:
    grefcallID = 0
    IF (refreshgref .OR. solvmixrefnd .OR. calcviscosity) THEN
        CALL SRgref(NK, T_K, XN, refreshgref, lnGaRref, XieRref)
    ENDIF

    CALL SRgres(X, lnGaR, NK, XieR)     !Residual part
    CALL SRgcomb(X, lnGaC, NK, XieC)    !Combinatorial part
    
    IF (calcviscosity) THEN
        lnGaSR(1:NK) = lnGaC(1:NK) + lnGaR(1:NK) - lnGaRref(1:NK)
        IF (NK > nneutral) THEN
            lnGaSR(nneutral+1:) = lnGaSR(nneutral+1:) - lnGaCinf(nneutral+1:)
        ENDIF
        CALL SRcalcvisc(NK, X, XN, lnGaSR, XieC, XieR)
    ELSE
        lnGaSR(1:NK) = lnGaC(1:NK) + lnGaR(1:NK) - lnGaRref(1:NK)
        IF (NK > nneutral) THEN
            lnGaSR(nneutral+1:) = lnGaSR(nneutral+1:) - lnGaCinf(nneutral+1:)
        ENDIF
    ENDIF
    
    END SUBROUTINE SRunifac
    !==========================================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Viscosity calculations based on contributions from SRunifac part for organics and  *
    !*   water as well as ModViscEyring in presence of ions.                                *
    !*                                                                                      *
    !*   :: Authors & Copyright ::                                                          *
    !*   Joseph Lilek, Natalie Gervasi, Andi Zuend,                                         *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2018                                                            *
    !*   -> latest changes: 2021-11-30                                                      *
    !*                                                                                      *
    !****************************************************************************************
    SUBROUTINE SRcalcvisc(NK, X, XN, lnGaSR, XieC, XieR)
    
    USE ModAIOMFACvar, ONLY : ln_etamix, ln_eta0, lneta_cpn, ln_eta_aquelec
    USE ModViscEyring, ONLY : AqueousElecViscosity, WaterMolefracCorrection, aquelec
    USE ModCompScaleConversion, ONLY : zSolution2SpeciesMolality
    
    IMPLICIT NONE
    !interface variables:
    INTEGER(4),INTENT(IN) :: NK
    REAL(8),DIMENSION(:),INTENT(IN) :: X, XN
    REAL(8),DIMENSION(:),INTENT(IN) :: lnGaSR, XieC, XieR
    !local variables:
    INTEGER(4) :: viscosity_mixmode
    REAL(8) :: RSS, ln_etaw
    REAL(8),DIMENSION(NK) :: lnGaC_loc, lnGaR_loc, XieC_loc, XieR_loc, phi, Xnew, lnGaSR_loc
    LOGICAL(4) :: ionspresent, orgspresent
    !...................................................
    
    orgspresent = .false. 
    ionspresent = .false.
    IF (nneutral > 1) THEN
        IF (SUM(X(2:nneutral)) > 0.0D0) THEN
            orgspresent = .true.
        ENDIF
    ENDIF
    IF (nneutral < NK) THEN
        IF (SUM(X(nneutral+1:NK)) > 0.0D0) THEN
            ionspresent = .true.
        ENDIF
    ENDIF
    
    !Select the mixing model affecting the viscosity calculation in mixed organic--inorganic (ion-containing) solutions.
    !Note that aquelec, aquorg or ZSR-based organic--inorganic mixing is set as option within 'ModViscEyring' (default is aquelec)
    viscosity_mixmode = 1
    IF (aquelec .AND. ionspresent) THEN     
        viscosity_mixmode = 1
    ELSE
        viscosity_mixmode = 2
    ENDIF
    IF (ionspresent) THEN
        CALL SRgres(XN, lnGaR_loc, NK, XieR_loc)        !Residual part using XN
        CALL SRgcomb(XN, lnGaC_loc, NK, XieC_loc)       !Combinatorial part using XN
    ENDIF

    SELECT CASE(viscosity_mixmode)
    CASE(1)     !"aquelec" org-inorg mixing used (and/or ions present in solution)
        IF (orgspresent) THEN
            Xnew = 0.0D0
            Xnew(1) = X(1) + SUM(X(nneutral+1:NK))  !use Xnew because XieC for aqueous elec solution becomes very small at low water content, counteracting the scaling of XieC by eta_aquelec
            Xnew(2:nneutral) = X(2:nneutral)
            Xnew = Xnew/SUM(Xnew)
            RSS = SUM(RS*Xnew)  
            phi = Xnew*RS/RSS
        ELSE
            RSS = SUM(RS*XN)
            phi = XN*RS/RSS        
            Xnew = XN
        ENDIF
        CALL AqueousElecViscosity(X, lnGaSR, RS, ln_eta0(1), ln_eta_aquelec)
        XieC_loc(1) = XieC_loc(1)*(ln_eta_aquelec/ln_eta0(1))*Xnew(1)/XN(1)     !the XieC of water corrected for by ion effects on the "pure" component viscosity of water.

        lneta_cpn(1:nneutral) = XieC_loc(1:nneutral) + phi(1:nneutral)*(XieR_loc(1:nneutral) - XieRref(1:nneutral))
        lneta_cpn(nneutral+1:) = 0.0D0
        ln_etamix = SUM(lneta_cpn(1:nneutral))
        
    CASE(2)     !"aquorg" org-inorg mixing is used (or ion-free mixture)          
        !Calculate the viscosity of the (electrolyte-free) aqueous organic mixture
        RSS = SUM(RS(:)*XN(:))              !calculating phi over neutrals only because the ion contributions are only impacting ln_eta0(1)
        phi = XN*RS/RSS      
        IF (ionspresent) THEN
            lnGaSR_loc(1:NK) = lnGaC_loc(1:NK) + lnGaR_loc(1:NK) - lnGaRref(1:NK)
            IF (NK > nneutral) THEN
                lnGaSR_loc(nneutral+1:) = lnGaSR_loc(nneutral+1:) - lnGaCinf(nneutral+1:)
            ENDIF
            lneta_cpn(1:nneutral) = XieC_loc(1:nneutral) + phi(1:nneutral)*(XieR_loc(1:nneutral) - XieRref(1:nneutral))
        ELSE
            lnGaSR_loc = lnGaSR                          !in this case, use the calculated values from SRunifac directly
            lneta_cpn(1:nneutral) = XieC(1:nneutral) + phi(1:nneutral)*(XieR(1:nneutral) - XieRref(1:nneutral))
        ENDIF
        lneta_cpn(nneutral+1:) = 0.0D0
        ln_etamix = SUM(lneta_cpn(1:nneutral))
        
        IF (ionspresent) THEN                           !the call to AqueousElecViscosity is only necessary in presence of ions
            ln_etaw = ln_etamix                         !assign the ln_etamix value from the electrolyte-free aqueous organic mixture as the "pure-water" value in this mode.
            IF (orgspresent) THEN
                CALL WaterMolefracCorrection(X, Xnew)   !adds mass of neutrals to water mass, then converts to mole fractions Xnew
            ELSE
                Xnew = X
            ENDIF
            CALL AqueousElecViscosity(Xnew, lnGaSR_loc, RS, ln_etaw, ln_etamix)
            !the returned ln_etamix will account for ion effects...
        ENDIF
    END SELECT
    
    END SUBROUTINE SRcalcvisc
    !==========================================================================================================================																																														 
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Calculation of the residual pure component reference activity coefficients and     *
    !*   reference mixture viscosity contributions.                                         *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
    !*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA                *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2004                                                            *
    !*   -> latest changes: 2021/01/31                                                     *
    !*                                                                                      *
    !****************************************************************************************
    SUBROUTINE SRgref(NK, T_K, XN, grefresh, lnGaRrefer, XieRrefer)
    
    USE ModAIOMFACvar, ONLY : Tglass0, fragil 
    USE ModPureCompViscos, ONLY : PureCompViscosity

    IMPLICIT NONE
    !interface variables:
    INTEGER(4),INTENT(IN) :: NK
    REAL(8),INTENT(IN) :: T_K
    REAL(8),DIMENSION(:),INTENT(IN) :: XN
    LOGICAL(4),INTENT(IN) :: grefresh
    REAL(8),DIMENSION(:),INTENT(OUT) :: lnGaRrefer, XieRrefer
    !local variables:
    INTEGER(4) :: I, iflag, K
    REAL(8),PARAMETER :: T0 = 298.15D0  !room temperature as the reference temperature for the calculation of the PsiT array
    REAL(8),DIMENSION(NK) :: Xrefsp, lnGaRX, XieX
    REAL(8) :: lneta, Tglass, D         !pure component viscosity, glass transition temperature, fragility
    !...................................................

    !Check / set temperature-dependent parameters:
    IF (grefresh) THEN !detected a change in temperature => PsiT needs to be updated
        SELECT CASE(nd)
        CASE(500:799, 2000:2434)    !(500:799,2000:2500)  2000:2434
            !new equation for UNIFAC/AIOMFAC 3-parameter version for larger temperature range; see Ganbavale et al. (2015, ACP, doi:10.5194/acp-15-447-2015)
            PsiT = EXP( -(parA/T_K) + parB*(1.0D0/T0 -1.0D0/T_K) + parC*( ((T0 - T_K)/T_K) +LOG(T_K/T0) ) ) 
        CASE DEFAULT
            PsiT = EXP(-parA/T_K)  !{original UNIFAC ("good" near room temperature: PsiT = EXP(-parA/T_K)} parameterization
        END SELECT
        DO I = 1,NG
            PsiT_dimflip(I,:) = PsiT(:,I) !store the PsiT data (in PsiT_dimflip) with the first dimension being the second and vice-versa.
        ENDDO
        !------
        !loop over neutral components and calculate the residual reference contributions from the different subgroups
        !reference residual contributions of pure components (and also hypothetical pure ions) are always 1.0 (i.e. for water) and
        !the calculations can therefore be omitted for components consisting of one subgroup only.
        Xrefsp = 0.0D0
        lnGaRrefer = 0.0D0
        XieRrefer = 0.0D0
        DO I = 1,nneutral
            !check if number of subgroups in component is greater than 1:
            K = SUM(SRNY_dimflip(1:NGN,I)) !SUM(SRNY(I,1:NGN))
            IF (K > 1) THEN !more than one subgroup
                grefcallID = I
                Xrefsp(1:nneutral) = 0.0D0
                Xrefsp(I) = 1.0D0 !set all other mole fractions but this (I) to 0.0. Thus reference state conditions of x(i) = 1
                CALL SRgres(Xrefsp, lnGaRX, NK, XieX)
                lnGaRrefer(I) = lnGaRX(I) !this line is necessary since the other values (not only I) get overwritten in call to SRgres!
                XieRrefer(I) = XieX(I)
                grefcallID = 0
            ENDIF
        ENDDO
    ENDIF

    !Check for solvent mixture reference state and potentially calculate the reference of the residual part for ions in solvent mixture reference state.
    IF (solvmixrefnd) THEN  !solvent mixture as reference
        grefcallID = 0
        DO I = nneutral+1,NK
            Xrefsp(1:nneutral) = XN(1:nneutral) !solvent
            Xrefsp(nneutral+1:) = 0.0D0 !ions
            CALL SRgres(Xrefsp, lnGaRX, NK, XieX)
            lnGaRrefer(I) = lnGaRX(I) !this line is necessary since the other values (not only I) get overwritten in call to SRgres!
            XieRrefer(I) = XieX(I)
        ENDDO
    ENDIF

    !---- calculate pure component viscosity of non-electrolytes at given temperature
    IF (calcviscosity) THEN
        DO I = 1,nneutral
            CALL PureCompViscosity(I, T_K, lneta, iflag, Tglass, D)   !calculate pure component dynamic viscosity eta
            IF (iflag == 0) THEN
                Tglass0(I) = Tglass
                fragil(I) = D
                ln_eta0(I) = lneta !eta in SI units of [Pa s]
            ELSE
                ln_eta0(I) = -7777.7D0  !to indicate a problem
                Tglass0(I) = -1.0D0
                fragil(I) = -1.0D0
                SELECT CASE(nd)
                CASE(500:799,1364:1366,2000:2500)
                    !do nothing for now....
                CASE DEFAULT
                    !$OMP CRITICAL
                    WRITE(*,*) ""
                    WRITE(*,*) "WARNING: a pure compound viscosity for component",I,"could not be set due to missing parameters or a temperature outside the available bounds!"
                    WRITE(*,*) "Dataset is nd = ", nd
                    !READ(*,*)
                    !$OMP END CRITICAL
                END SELECT
            ENDIF
        ENDDO
        IF (NK > nneutral) THEN !calculate viscosity of (hypothetical) individual ions
            !ions / electrolyte components:
            ln_eta0(nneutral+1:NK) = -3.0D0     !initialization; this value does not matter.
        ENDIF
    ELSE
        ln_eta0(1:NK) = -9999.9D0 !a tiny value to indicate uninitialized ln_eta0
    ENDIF

    END SUBROUTINE SRgref
!==========================================================================================================================

    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Calculation of residual activity coefficients [less the reference part stemming    *
    !*   from group k in pure component i (the latter is calculated in gref)].              *
    !*   Component residual contributions to mixture viscosity are also calculated.         *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
    !*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA                *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2004                                                            *
    !*   -> latest changes: 2018/05/18                                                      *
    !*                                                                                      *
    !****************************************************************************************
    PURE SUBROUTINE SRgres(Xr, lnGaR, NK, XieR)

    IMPLICIT NONE
    !interface variables:
    INTEGER(4),INTENT(IN) :: NK
    REAL(8),DIMENSION(:),INTENT(OUT) :: lnGaR, XieR
    REAL(8),DIMENSION(:),INTENT(IN) :: Xr
    !local variables:
    INTEGER(4) :: I, k
    REAL(8) :: S2, S3
    REAL(8),DIMENSION(NG) :: TH, GAML, S1, XG, S4, THbyS4, Xiek
    REAL(8),DIMENSION(NG,NG) :: THmn
    !...................................................

    IF (grefcallID > 0) THEN !this is only the case when called from SRgref
        DO I = 1,NG
            S1(I) = SRNY(grefcallID,I) !only this because all Xr other than that of grefcallID are zero in this reference case
        ENDDO
    ELSE !typical case
        DO I = 1,NG
            S1(I) = SUM(SRNY(1:NK,I)*Xr(1:NK))
        ENDDO
    ENDIF
    S2 = SUM(S1)
    XG = S1/S2
    S3 = SUM(Q(1:NG)*XG)
    TH = Q(1:NG)*XG/S3
    DO I = 1,NG
        S4(I) = SUM(TH*PsiT(1:NG,I))
    ENDDO
    THbyS4 = TH/S4
    DO I = 1,NG
        GAML(I) = Q(I)*( 1.0D0-LOG(S4(I)) -SUM(THbyS4*PsiT_dimflip(1:NG,I)) )
    ENDDO
    IF (grefcallID > 0) THEN !this is only the case when called from SRgref
        lnGaR(grefcallID) = SUM(SRNY_dimflip(1:NG,grefcallID)*GAML(1:NG))
    ELSE
        DO I = 1,NK
            lnGaR(I) = SUM(SRNY_dimflip(1:NG,I)*GAML(1:NG))
        ENDDO
    ENDIF

    !--------------------------------------------------------
    !Viscosity of mixture calculation; residual part.
    IF (calcviscosity) THEN
        DO k = 1,NG !k
            DO I = 1,NG !m
                S4(I) = SUM(TH(1:NG)*PsiT(1:NG,I))
                THmn(k,I) = TH(k)*PsiT_dimflip(I,k)/S4(I)
            ENDDO
        ENDDO
        !sum up the number of subgroups and their contributions to XieR in compound I:
        IF (grefcallID > 0) THEN !this is only the case when called from SRgref
            DO k = 1,NG
                Xiek(k) = (Q(k)/R(k))*Nvis(k,grefcallID)*SUM(THmn(1:NG,k)*LOG(PsiT(1:NG,k))) !Natalie mod 4
            ENDDO
            XieR(grefcallID) = SUM(SRNY_dimflip(1:NG,grefcallID)*Xiek(1:NG))
        ELSE
            DO I = 1,NK
                !calculate viscosity contributions by individual subgroups:
                DO k = 1,NG
                    Xiek(k) = (Q(k)/R(k))*Nvis(k,I)*SUM(THmn(1:NG,k)*LOG(PsiT(1:NG,k))) !Natalie mod 4
                ENDDO
                XieR(I) = SUM(SRNY_dimflip(1:NG,I)*Xiek(1:NG))
            ENDDO
        ENDIF
    ELSE !a default value indicating "NO viscosity calc."
        XieR = 0.0D0
    ENDIF
    !--------------------------------------------------------

    END SUBROUTINE SRgres
!==========================================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Calculation of combinatorial UNIFAC / AIOMFAC part for activity coefficients and   *
    !*   mixture viscosity contributions of components.                                     *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
    !*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA                *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2004                                                            *
    !*   -> latest changes: 2018/05/18                                                      *
    !*                                                                                      *
    !****************************************************************************************
    SUBROUTINE SRgcomb(X, lnGaC, NK, XieC)

    IMPLICIT NONE
    !interface variables:
    INTEGER(4),INTENT(IN) :: NK
    REAL(8),DIMENSION(:),INTENT(IN) :: X
    REAL(8),DIMENSION(:),INTENT(OUT) :: lnGaC, XieC
    !local variables:
    INTEGER(4) :: nnp1  
    REAL(8) :: QSSref, RSSref, QSS, RSS, Xsolvtot, XLSref
    REAL(8),DIMENSION(NK) :: B, Xsolv, ViV, phi, phiQ, QIQ, VQR !!, Xnew
    !.................................................................................
    nnp1 = nneutral+1
    QSS = SUM(QS(:)*X(:))
    RSS = SUM(RS(:)*X(:))
    ViV(:) = RS(:)/RSS ! the van der Waals volume of molecule to average volume ratio of (1:NK) molecules

    !orig. UNIFAC-equivalent formulation as described by Voutsas and Tassios (1997), Ind. Eng. Chem. Res. 1997, 36, 4965-4972
    VQR(:) = ViV(:)*QSS/QS(:)
    lnGaC(:) = LOG(ViV(:)) +1.0D0 -ViV(:) -5.0D0*QS(:)*(LOG(VQR(:)) +1.0D0 -VQR(:))
    
    !Calculation of "lnGaC infinite" for the ions. Water is taken as the reference solvent (usually):
    !decide whether water or the solvent mixture is the reference solvent:
    IF (solvmixrefnd) THEN  !solvent mixture as reference for solutes (ions)
        lnGaCinf = 0.0D0
        Xsolvtot = SUM(X(1:nneutral))
        Xsolv(1:nneutral) = X(1:nneutral)/Xsolvtot
        QSSref = SUM(QS(1:nneutral)*Xsolv(1:nneutral))
        RSSref = SUM(RS(1:nneutral)*Xsolv(1:nneutral))
        XLSref = SUM(XL(1:nneutral)*Xsolv(1:nneutral))
        B(nnp1:NK) = 5.0D0*QS(nnp1:NK)*LOG(QS(nnp1:NK)/QSSref*RSSref/RS(nnp1:NK)) +XL(nnp1:NK) -RS(nnp1:NK)/RSSref*XLSref
        lnGaCinf(nnp1:NK) = LOG(RS(nnp1:NK)/RSSref) +B(nnp1:NK)
    ELSE IF (gcombrefresh) THEN     !water is the reference solvent for solutes (ions); 
                                    !GAMCinf is a constant for the ions of the system and needs to be calculated only at the first call.
        lnGaCinf = 0.0D0
        B(nnp1:NK) = 5.0D0*QS(nnp1:NK)*LOG(QS(nnp1:NK)/1.40D0*0.92D0/RS(nnp1:NK)) +XL(nnp1:NK) -RS(nnp1:NK)/0.92D0*(-2.32D0)
        lnGaCinf(nnp1:NK) = LOG(RS(nnp1:NK)/0.92D0) +B(nnp1:NK)
        gcombrefresh = .false.
    ENDIF

    !--------------------------------------------------------
    !Viscosity of mixture calculation; combinatorial part.
    IF (calcviscosity) THEN
        phi = X*ViV !volume based component fraction
        QiQ = QS(1:NK)/QSS
        phiQ = X*QiQ
        !undefined ln_eta0 would contain negative values of large magnitude
        IF (ALL(ln_eta0(1:nneutral) > -7777.0D0)) THEN !pure-component values available (replaced NK with nneutral (JL november 6)
            XieC(1:nneutral) = (EXP(lnGaC(1:nneutral))*X(1:nneutral))*ln_eta0(1:nneutral)        
        ELSE
            XieC = -55555.5D0
        ENDIF
    ELSE !a default value indicating "NO viscosity calc."
        XieC = -99999.9D0
    ENDIF
    !--------------------------------------------------------

    END SUBROUTINE SRgcomb
!==========================================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine SRsystm puts together the compacted parameter matrix from the large       *
    !*   UNIFAC parameter tables on the basis of information stored in ITABsr for a given   *
    !*   system.                                                                            *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
    !*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA                *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2004                                                            *
    !*   -> latest changes: 2018/05/18                                                      *
    !*                                                                                      *
    !****************************************************************************************
    SUBROUTINE SRsystm(NK)

    USE ModSystemProp, ONLY : errorflagmix

    IMPLICIT NONE

    INTEGER(4),INTENT(IN) :: NK
    INTEGER(4) :: J, JJ, I, J1, L, L1 !, unitm
    REAL(8) :: Roxyethylene
    !.....................end of variable declarations.............................

    !--------------------------------------------------------------------------------------------------------
    IF (ALLOCATED(parA)) THEN
        DEALLOCATE(parA, parB, parC, PsiT, PsiT_dimflip, SRNY, SRNY_dimflip, Nvis, R, Q, XL, RS, QS, lnGaCinf, lnGaRref, XieRref)
    ENDIF
    ALLOCATE( parA(NG,NG), parB(NG,NG), parC(NG,NG), PsiT(NG,NG), PsiT_dimflip(NG,NG), SRNY(NK,NG), SRNY_dimflip(NG,NK), &
        & Nvis(NG,NK), R(NG), Q(NG), XL(NK), RS(NK), QS(NK), lnGaCinf(NK), lnGaRref(NK), XieRref(NK) ) !, XieCref(NK)

    parA = 0.0D0
    parB = 0.0D0
    parC = 0.0D0
    R = 0.0D0
    Q = 0.0D0
    XL = 0.0D0
    RS = 0.0D0  ! sum of the Van der Waals group volumes --> component molar volume
    QS = 0.0D0  ! sum of the Van der Waals group surface area --> component molar surface area
    SRNY = 0
    SRNY_dimflip = 0
    Nvis = 0.0D0
    lastTK = 0.0D0
    gcombrefresh = .true. !if .true., the infinite dilution reference value for the SR combinatorial part of ions will be calculated in SRgcomb and saved for subsequent use.

    !populate array SRNY for SR part
    DO J = 1,NG
        JJ = AllSubs(J)
        R(J) = SR_RR(JJ)
        Q(J) = SR_QQ(JJ)
        SRNY(:,J) = ITABsr(:,JJ)
        IF (JJ == 154) THEN !special PEG-group present in the mixture, so it is a "PEG system"
            Roxyethylene = SR_RR(JJ)
        ENDIF
        SRNY_dimflip(J,:) = SRNY(:,J) !same data but in flipped-dimension storage order for use with better memory alignment in certain equations
    ENDDO

    !assign interaction parameters between different subgroups from  main group interaction parameter matrix.
    !loop only over non-electrolyte (neutral solvent) groups as the interactions of neutral groups with ions is defined to be always 0.0D0 in SR-part.
    DO J = 1,NGN
        J1 = Imaingroup(maingrindexofsubgr(J)) !this maingrindexofsubgr array assigns the index of the main group (of the current mixture) associated with neutral subgroup index j
        DO L = 1,NGN
            L1 = Imaingroup(maingrindexofsubgr(L))
            parA(L,J) = ARR(L1,J1)
            parB(L,J) = BRR(L1,J1)
            parC(L,J) = CRR(L1,J1)
            !Apply check to avoid using SR inteaction parameters that are not assigned correctly or were not estimated yet:
            IF (ARR(L1,J1) < -888887.0D0) THEN  !this parameter has not yet been estimated nor is a fitparameter set, so don't use it!!
                errorflagmix = 14
                WRITE(*,*) ""
                WRITE(*,*) "======================================================="
                WRITE(*,*) "        WARNING from SRsystem:"
                WRITE(*,*) "A neutral main group <-> main group interaction coeff. is "
                WRITE(*,*) "not yet defined or no fit parameter was assigned. "
                WRITE(*,*) "System nd: ", nd
                WRITE(*,*) "L1, J1, ARR(L1,J1): ", L1, J1, ARR(L1,J1)
                WRITE(*,*) "======================================================="
                WRITE(*,*) ""
                READ(*,*)               !to READ(*,*) and wait for user interaction
            ENDIF
            SELECT CASE(nd)
            CASE(500:799,2000:2500)     !potentially used for 3-parameter AIOMFAC (UNIFAC) part for non-electrolyte systems
                IF (BRR(L1,J1) < -888887.0D0 .OR. CRR(L1,J1) < -888887.0D0) THEN    !this parameter has not yet been estimated nor is a fitparameter set, so don't use it!!
                    errorflagmix = 15
                    WRITE(*,*) ""
                    WRITE(*,*) "======================================================="
                    WRITE(*,*) "        WARNING from SRsystem:"
                    WRITE(*,*) "A neutral main group <-> main group interaction coeff. is "
                    WRITE(*,*) "not yet defined or no fit parameter was assigned. "
                    WRITE(*,*) "System nd: ", nd
                    WRITE(*,*) "L1, J1: ", L1, J1
                    WRITE(*,*) "BRR(L1,J1), CRR(L1,J1): ", BRR(L1,J1), CRR(L1,J1)
                    WRITE(*,*) "======================================================="
                    WRITE(*,*) ""
                    READ(*,*) !to READ(*,*) and wait for user interaction
                ENDIF
            END SELECT
        ENDDO
    ENDDO

    DO I = 1,NK
        DO J = 1,NG
            RS(I) = RS(I) +SRNY_dimflip(J,I)*R(J) !RS(I) +SRNY(I,J)*R(J)
            QS(I) = QS(I) +SRNY_dimflip(J,I)*Q(J) !QS(I) +SRNY(I,J)*Q(J)
            IF (isPEGsystem) THEN !PEG-polymer solution; use special subgroups and some scaled parameters
                IF (AllSubs(J) == 154 .AND. SRNY_dimflip(J,I) > 0) THEN  !oxyethylene group in PEG polymer, so use a special parametrisation for RS(I) in those cases (idea adapted from Ninni et al., (1999), but implemented in different form, basically to fit RS(I)/QS(I) of oxyethylene group rather than taking fixed Bondi (1964) values.):
                    RS(I) = RS(I) -SRNY_dimflip(J,I)*R(J) !reset above set value for this group.
                    R(J) = 1.381433094376D0     !fitSRparam(211)
                    RS(I) = RS(I) + SRNY_dimflip(J,I)*R(J)
                    !adjust Q(CH2OCH2) value according to the relation between z, RS(I) and QS(I) described in Vera et al. (1977) for straight chain molecules. The expression below for Q(CH2OCH2) is dependent on the overall chain length of the PEG polymer.
                    QS(I) = QS(I) -SRNY_dimflip(J,I)*Q(J) !reset above set value for this group.
                    Q(J) = 3.0D0    !fitSRparam(212) !as a test with limiting ratio of Q(J)/R(J) for a large number of monomer units
                    QS(I) = QS(I) + SRNY_dimflip(J,I)*Q(J)  !CH2OCH2[PEG]
                ENDIF
            ENDIF
        ENDDO
        XL(I) = 5.0D0*(RS(I)-QS(I)) -RS(I) +1.0D0  !corresponds to l(i) eq. (3) in UNIFAC Fredenslund et al. (1975) paper
        !...
        !calculate Nvis for viscosity calculations according to Cao et al. (1993), Ind. Eng. Chem. Res., 32, 2088--2092.
        DO J = 1,NG
            Nvis(J,I) = Q(J)*(0.5D0*(QS(I)-RS(I))-0.1D0*(1.0D0-RS(I)))
        ENDDO
    ENDDO

    END SUBROUTINE SRsystm
!==========================================================================================================================

END MODULE ModSRunifac