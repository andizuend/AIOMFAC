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
!*   -> latest changes: 2024-04-22                                                      *
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
!*   -  subroutine SRunifac                                                             *
!*   -  subroutine SRcalcvisc                                                           *
!*   -  subroutine SRgref                                                               *
!*   -  subroutine SRgres                                                               *
!*   -  subroutine SRgcomb                                                              * 
!*   -  subroutine SRsystm                                                              *
!*                                                                                      *
!****************************************************************************************
module ModSRunifac

use Mod_kind_param, only : wp
use ModAIOMFACvar, only : lastTK, ln_eta0
use ModSystemProp, only : nd, nneutral, NG, NGN, Allsubs, ITABsr, Nmaingroups, topsubno, &
    & maingrindexofsubgr, Imaingroup, isPEGsystem, calcviscosity, solvmixrefnd, COMPN
use ModSRparam, only : ARR, BRR, CRR, SR_RR, SR_QQ

implicit none
!private module variables
integer,private :: grefcallID
integer,dimension(:,:),allocatable,private :: SRNY, SRNY_dimflip
real(wp),dimension(:),allocatable,private :: RS, QS, R, Q, XL, lnGaCinf, lnGaRref, XieRref
real(wp),dimension(:,:),allocatable,private :: parA, parB, parC, PsiT, PsiT_dimflip
real(wp),dimension(:,:),allocatable,private :: Nvis
logical,private :: gcombrefresh

!$OMP THREADPRIVATE(parA, parB, parC, PsiT, PsiT_dimflip, grefcallID, SRNY, SRNY_dimflip, R, Q, RS, QS, XL, Nvis, &
   !$OMP & lnGaCinf, lnGaRref, XieRref, gcombrefresh)

!==========================================================================================================================
    contains
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
    !*   -> latest changes: 2024-04-22                                                      *
    !*                                                                                      *
    !****************************************************************************************
    
    !   lnGaC is the natural log of the combinatorial activity coefficient.
    !   lnGaRref is the natural log of the reference residual activity coefficient (pure component reference).
    !   lnGaR is the natural log of the residual activity coefficient of the components in the mixture.
    !   lnGaSR is the natural log of the total short-range activity coefficient.
    !   ITABsr lists the components and corresponding subgroups, sorted such that cations and anions are individual components (for SR part)
    !   XieC and XieR are the non-electrolyte component viscosity contributions of the combinatorial and residual parts, respectively.
    
    subroutine SRunifac(NK, T_K, X, XN, refreshgref, lnGaSR)

    use ModPureCompViscos, only : PureCompViscosity

    implicit none
    !interface:
    integer,intent(in) :: NK                                ![-] number of species, counting each ion separately
    real(wp),intent(in) :: T_K                              ![K] temperature
    real(wp),dimension(:),intent(in) :: X, XN               ![-] mole fraction array X; structure is: 
                                                            !   1) neutral components in component order,
                                                            !   2) ions: first the cations, then the anions
    logical,intent(in) :: refreshgref 
    real(wp),dimension(:),intent(out) :: lnGaSR             ![-] ln(total short-range activity coeff.)
    !local variables:
    real(wp),dimension(NK) :: lnGaC, lnGaR, XieC, XieR
    !.......................................................
    
    !Determine the reference values (lnGaRref, XieRref) for the residual part:
    grefcallID = 0
    if (refreshgref .OR. solvmixrefnd .OR. calcviscosity) then
        call SRgref(NK, T_K, XN, refreshgref, lnGaRref, XieRref)
    endif

    call SRgres(X, lnGaR, NK, XieR)                         !Residual part
    call SRgcomb(X, lnGaC, NK, XieC)                        !Combinatorial part
    
    if (calcviscosity) then
        lnGaSR(1:NK) = lnGaC(1:NK) + lnGaR(1:NK) - lnGaRref(1:NK)
        if (NK > nneutral) then
            lnGaSR(nneutral+1:) = lnGaSR(nneutral+1:) - lnGaCinf(nneutral+1:)
        endif
        call SRcalcvisc(NK, X, XN, lnGaSR, XieC, XieR)
    else
        lnGaSR(1:NK) = lnGaC(1:NK) + lnGaR(1:NK) - lnGaRref(1:NK)
        if (NK > nneutral) then
            lnGaSR(nneutral+1:) = lnGaSR(nneutral+1:) - lnGaCinf(nneutral+1:)
        endif
    endif
    
    end subroutine SRunifac
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
    subroutine SRcalcvisc(NK, X, XN, lnGaSR, XieC, XieR)
    
    use ModAIOMFACvar, only : ln_etamix, ln_eta0, lneta_cpn, ln_eta_aquelec, aquelecVar_save
    use ModViscEyring, only : AqueousElecViscosity, WaterMolefracCorrection, aquelec
    use ModCompScaleConversion, only : zSolution2SpeciesMolality
    
    implicit none
    !interface variables:
    integer,intent(in) :: NK
    real(wp),dimension(:),intent(in) :: X, XN
    real(wp),dimension(:),intent(in) :: lnGaSR, XieC, XieR
    !local variables:
    integer :: viscosity_mixmode
    real(wp) :: RSS, ln_etaw
    real(wp),dimension(NK) :: lnGaC_loc, lnGaR_loc, XieC_loc, XieR_loc, phi, Xnew, lnGaSR_loc
    logical :: ionspresent, orgspresent
    !...................................................
    
    orgspresent = .false. 
    ionspresent = .false.
    if (nneutral > 1) then
        if (sum(X(2:nneutral)) > 0.0_wp) then
            orgspresent = .true.
        endif
    endif
    if (nneutral < NK) then
        if (sum(X(nneutral+1:NK)) > 0.0_wp) then
            ionspresent = .true.
        endif
    endif
    
    !Select the mixing model affecting the viscosity calculation in mixed organic--inorganic (ion-containing) solutions.
    !Note that aquelec, aquorg or ZSR-based organic--inorganic mixing is set as option within 'ModViscEyring' (default is aquelec)
    if (aquelec .AND. ionspresent) then     
        viscosity_mixmode = 1
    else
        viscosity_mixmode = 2
    endif
    if (ionspresent) then
        call SRgres(XN, lnGaR_loc, NK, XieR_loc)        !Residual part using XN
        call SRgcomb(XN, lnGaC_loc, NK, XieC_loc)       !Combinatorial part using XN
    endif

    select case(viscosity_mixmode)
    case(1)     !"aquelec" org-inorg mixing used (and/or ions present in solution)
        if (orgspresent) then
            Xnew = 0.0_wp
            Xnew(1) = X(1) + sum(X(nneutral+1:NK))  !use Xnew because XieC for aqueous elec solution becomes very small at low water content, counteracting the scaling of XieC by eta_aquelec
            Xnew(2:nneutral) = X(2:nneutral)
            Xnew = Xnew/sum(Xnew)
            RSS = sum(RS*Xnew)  
            phi = Xnew*RS/RSS
        else
            RSS = sum(RS*XN)
            phi = XN*RS/RSS        
            Xnew = XN
        endif
        call AqueousElecViscosity(X, lnGaSR, RS, ln_eta0(1), ln_eta_aquelec)
        if (XN(1) > 0.0_wp) then
            XieC_loc(1) = XieC_loc(1)*(ln_eta_aquelec/ln_eta0(1))*Xnew(1)/XN(1)     !the XieC of water corrected for by ion effects on the "pure" component viscosity of water.
        else
            XieC_loc(1) = 0.0_wp 
        endif

        lneta_cpn(1:nneutral) = XieC_loc(1:nneutral) + phi(1:nneutral)*(XieR_loc(1:nneutral) - XieRref(1:nneutral))
        lneta_cpn(nneutral+1:) = 0.0_wp
        ln_etamix = sum(lneta_cpn(1:nneutral))
        
        !save these values (XieC_loc, phi, XieR_loc for aqueous electrolyte solution)
        aquelecVar_save = 0.0_wp
        aquelecVar_save(1) = XieC_loc(1)
        aquelecVar_save(3) = XN(1)
        if (nneutral > 1) then !data for debugging and diagnostics
            aquelecVar_save(2) = XieC_loc(2)
            aquelecVar_save(4) = XN(2)
        endif
        
    case(2)     !"aquorg" org-inorg mixing is used (or ion-free mixture)          
        !Calculate the viscosity of the (electrolyte-free) aqueous organic mixture
        RSS = sum(RS(:)*XN(:))              !calculating phi over neutrals only because the ion contributions are only impacting ln_eta0(1)
        phi = XN*RS/RSS      
        if (ionspresent) then
            lnGaSR_loc(1:NK) = lnGaC_loc(1:NK) + lnGaR_loc(1:NK) - lnGaRref(1:NK)
            if (NK > nneutral) then
                lnGaSR_loc(nneutral+1:) = lnGaSR_loc(nneutral+1:) - lnGaCinf(nneutral+1:)
            endif
            lneta_cpn(1:nneutral) = XieC_loc(1:nneutral) + phi(1:nneutral)*(XieR_loc(1:nneutral) - XieRref(1:nneutral))
        else
            lnGaSR_loc = lnGaSR                          !in this case, use the calculated values from SRunifac directly
            lneta_cpn(1:nneutral) = XieC(1:nneutral) + phi(1:nneutral)*(XieR(1:nneutral) - XieRref(1:nneutral))
        endif
        lneta_cpn(nneutral+1:) = 0.0_wp
        ln_etamix = sum(lneta_cpn(1:nneutral))
        
        if (ionspresent) then                           !the call to AqueousElecViscosity is only necessary in presence of ions
            ln_etaw = ln_etamix                         !assign the ln_etamix value from the electrolyte-free aqueous organic mixture as the "pure-water" value in this mode.
            if (orgspresent) then
                call WaterMolefracCorrection(X, Xnew)   !adds mass of neutrals to water mass, then converts to mole fractions Xnew
            else
                Xnew = X
            endif
            call AqueousElecViscosity(Xnew, lnGaSR_loc, RS, ln_etaw, ln_etamix)
            !the returned ln_etamix will account for ion effects...
        endif
    end select
    
    end subroutine SRcalcvisc
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
    !*   -> latest changes: 2021-12-15                                                     *
    !*                                                                                      *
    !****************************************************************************************
    subroutine SRgref(NK, T_K, XN, grefresh, lnGaRrefer, XieRrefer)
    
    use ModAIOMFACvar, only : Tglass0, fragil 
    use ModPureCompViscos, only : PureCompViscosity

    implicit none
    !interface variables:
    integer,intent(in) :: NK
    real(wp),intent(in) :: T_K
    real(wp),dimension(:),intent(in) :: XN
    logical,intent(in) :: grefresh
    real(wp),dimension(:),intent(out) :: lnGaRrefer, XieRrefer
    !local variables:
    integer :: I, iflag, K
    real(wp),parameter :: T0 = 298.15_wp  !room temperature as the reference temperature for the calculation of the PsiT array
    real(wp),dimension(NK) :: Xrefsp, lnGaRX, XieX
    real(wp) :: lneta, Tglass, D         !pure component viscosity, glass transition temperature, fragility
    !...................................................

    !Check / set temperature-dependent parameters:
    if (grefresh) then !detected a change in temperature => PsiT needs to be updated
        select case(nd)
        case(500:800, 2000:2434)    !(500:800,2000:2500)  2000:2434
            !new equation for UNIFAC/AIOMFAC 3-parameter version for larger temperature range; see Ganbavale et al. (2015, ACP, doi:10.5194/acp-15-447-2015)
            PsiT = exp( -(parA/T_K) + parB*(1.0_wp/T0 -1.0_wp/T_K) + parC*( ((T0 - T_K)/T_K) +log(T_K/T0) ) ) 
        case default
            PsiT = exp(-parA/T_K)  !{original UNIFAC ("good" near room temperature: PsiT = exp(-parA/T_K)} parameterization
        end select
        do I = 1,NG
            PsiT_dimflip(I,:) = PsiT(:,I) !store the PsiT data (in PsiT_dimflip) with the first dimension being the second and vice-versa.
        enddo
        !------
        !loop over neutral components and calculate the residual reference contributions from the different subgroups
        !reference residual contributions of pure components (and also hypothetical pure ions) are always 1.0 (i.e. for water) and
        !the calculations can therefore be omitted for components consisting of one subgroup only.
        Xrefsp = 0.0_wp
        lnGaRrefer = 0.0_wp
        XieRrefer = 0.0_wp
        do I = 1,nneutral
            !check if number of subgroups in component is greater than 1:
            K = sum(SRNY_dimflip(1:NGN,I)) !sum(SRNY(I,1:NGN))
            if (K > 1) then !more than one subgroup
                grefcallID = I
                Xrefsp(1:nneutral) = 0.0_wp
                Xrefsp(I) = 1.0_wp !set all other mole fractions but this (I) to 0.0. Thus reference state conditions of x(i) = 1
                call SRgres(Xrefsp, lnGaRX, NK, XieX)
                lnGaRrefer(I) = lnGaRX(I) !this line is necessary since the other values (not only I) get overwritten in call to SRgres!
                XieRrefer(I) = XieX(I)
                grefcallID = 0
            endif
        enddo
    endif

    !Check for solvent mixture reference state and potentially calculate the reference of the residual part for ions in solvent mixture reference state.
    if (solvmixrefnd) then  !solvent mixture as reference
        grefcallID = 0
        do I = nneutral+1,NK
            Xrefsp(1:nneutral) = XN(1:nneutral) !solvent
            Xrefsp(nneutral+1:) = 0.0_wp !ions
            call SRgres(Xrefsp, lnGaRX, NK, XieX)
            lnGaRrefer(I) = lnGaRX(I) !this line is necessary since the other values (not only I) get overwritten in call to SRgres!
            XieRrefer(I) = XieX(I)
        enddo
    endif

    !---- calculate pure component viscosity of non-electrolytes at given temperature
    if (calcviscosity) then
        do I = 1,nneutral
            call PureCompViscosity(I, T_K, lneta, iflag, Tglass, D)   !calculate pure component dynamic viscosity eta
            if (iflag == 0) then
                Tglass0(I) = Tglass
                fragil(I) = D
                ln_eta0(I) = lneta !eta in SI units of [Pa s]
            else
                ln_eta0(I) = -7777.7_wp  !to indicate a problem
                Tglass0(I) = -1.0_wp
                fragil(I) = -1.0_wp
                select case(nd)
                case(500:800, 1364:1366, 2000:2500)
                    !do nothing for now....
                case default
                    !$OMP CRITICAL
                    write(*,*) ""
                    write(*,*) "WARNING: a pure compound viscosity for component",I, "could not be set due to missing parameters or a temperature outside the available bounds!"
                    write(*,*) "Dataset is nd = ", nd
                    !read(*,*)
                    !$OMP end CRITICAL
                end select
            endif
        enddo
        if (NK > nneutral) then !calculate viscosity of (hypothetical) individual ions
            !ions / electrolyte components:
            ln_eta0(nneutral+1:NK) = -3.0_wp     !initialization
        endif
    else
        ln_eta0(1:NK) = -9999.9_wp !a tiny value to indicate uninitialized ln_eta0
    endif

    end subroutine SRgref
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
    pure subroutine SRgres(Xr, lnGaR, NK, XieR)

    implicit none
    !interface variables:
    integer,intent(in) :: NK
    real(wp),dimension(:),intent(out) :: lnGaR, XieR
    real(wp),dimension(:),intent(in) :: Xr
    !local variables:
    integer :: I, k
    real(wp) :: S2, S3
    real(wp),dimension(NG) :: TH, GAML, S1, XG, S4, THbyS4, Xiek
    real(wp),dimension(NG,NG) :: THmn
    !...................................................

    if (grefcallID > 0) then !this is only the case when called from SRgref
        do I = 1,NG
            S1(I) = SRNY(grefcallID,I) !only this because all Xr other than that of grefcallID are zero in this reference case
        enddo
    else !typical case
        do I = 1,NG
            S1(I) = sum(SRNY(1:NK,I)*Xr(1:NK))
        enddo
    endif
    S2 = sum(S1)
    XG = S1/S2
    S3 = sum(Q(1:NG)*XG)
    TH = Q(1:NG)*XG/S3
    do I = 1,NG
        S4(I) = sum(TH*PsiT(1:NG,I))
    enddo
    THbyS4 = TH/S4
    do I = 1,NG
        GAML(I) = Q(I)*( 1.0_wp-log(S4(I)) -sum(THbyS4*PsiT_dimflip(1:NG,I)) )
    enddo
    if (grefcallID > 0) then !this is only the case when called from SRgref
        lnGaR(grefcallID) = sum(SRNY_dimflip(1:NG,grefcallID)*GAML(1:NG))
    else
        do I = 1,NK
            lnGaR(I) = sum(SRNY_dimflip(1:NG,I)*GAML(1:NG))
        enddo
    endif

    !--------------------------------------------------------
    !Viscosity of mixture calculation; residual part.
    if (calcviscosity .OR. grefcallID > 0) then
        do k = 1,NG !k
            do I = 1,NG !m
                S4(I) = sum(TH(1:NG)*PsiT(1:NG,I))
                THmn(k,I) = TH(k)*PsiT_dimflip(I,k)/S4(I)
            enddo
        enddo
        !sum up the number of subgroups and their contributions to XieR in compound I:
        if (grefcallID > 0) then    !this is only the case when called from SRgref
            do k = 1,NG
                Xiek(k) = (Q(k)/R(k))*Nvis(k,grefcallID)*sum(THmn(1:NG,k)*log(PsiT(1:NG,k))) !Natalie mod 4
            enddo
            XieR(grefcallID) = sum(SRNY_dimflip(1:NG,grefcallID)*Xiek(1:NG))
            if (isPEGsystem) then  !presently special treatment for systems containing PEG oligomers (ignore residual contribution)
                XieR(grefcallID) = 0.0_wp
            endif
        else
            do I = 1,NK
                !calculate viscosity contributions by individual subgroups:
                do k = 1,NG
                    Xiek(k) = (Q(k)/R(k))*Nvis(k,I)*sum(THmn(1:NG,k)*log(PsiT(1:NG,k))) !Natalie mod 4
                enddo
                XieR(I) = sum(SRNY_dimflip(1:NG,I)*Xiek(1:NG))
                if (isPEGsystem) then  !presently special treatment for systems containing PEG oligomers (ignore residual contribution)
                    XieR(I) = 0.0_wp
                endif
            enddo
        endif
    else !a default value indicating "NO viscosity calc."
        XieR = 0.0_wp
    endif
    !--------------------------------------------------------

    end subroutine SRgres
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
    subroutine SRgcomb(X, lnGaC, NK, XieC)

    implicit none
    !interface variables:
    integer,intent(in) :: NK
    real(wp),dimension(:),intent(in) :: X
    real(wp),dimension(:),intent(out) :: lnGaC, XieC
    !local variables:
    integer :: nnp1  
    real(wp) :: QSSref, RSSref, QSS, RSS, Xsolvtot, XLSref
    real(wp),dimension(NK) :: B, Xsolv, ViV, phi, phiQ, QIQ, VQR !!, Xnew
    !.................................................................................
    nnp1 = nneutral+1
    QSS = sum(QS(:)*X(:))
    RSS = sum(RS(:)*X(:))
    ViV(:) = RS(:)/RSS ! the van der Waals volume of molecule to average volume ratio of (1:NK) molecules

    !orig. UNIFAC-equivalent formulation as described by Voutsas and Tassios (1997), Ind. Eng. Chem. Res. 1997, 36, 4965-4972
    VQR(:) = ViV(:)*QSS/QS(:)
    lnGaC(:) = log(ViV(:)) +1.0_wp -ViV(:) -5.0_wp*QS(:)*(log(VQR(:)) +1.0_wp -VQR(:))
    
    !Calculation of "lnGaC infinite" for the ions. Water is taken as the reference solvent (usually):
    !decide whether water or the solvent mixture is the reference solvent:
    if (solvmixrefnd) then  !solvent mixture as reference for solutes (ions)
        lnGaCinf = 0.0_wp
        Xsolvtot = sum(X(1:nneutral))
        Xsolv(1:nneutral) = X(1:nneutral)/Xsolvtot
        QSSref = sum(QS(1:nneutral)*Xsolv(1:nneutral))
        RSSref = sum(RS(1:nneutral)*Xsolv(1:nneutral))
        XLSref = sum(XL(1:nneutral)*Xsolv(1:nneutral))
        B(nnp1:NK) = 5.0_wp*QS(nnp1:NK)*log(QS(nnp1:NK)/QSSref*RSSref/RS(nnp1:NK)) +XL(nnp1:NK) -RS(nnp1:NK)/RSSref*XLSref
        lnGaCinf(nnp1:NK) = log(RS(nnp1:NK)/RSSref) +B(nnp1:NK)
    else if (gcombrefresh) then     !water is the reference solvent for solutes (ions); 
                                    !GAMCinf is a constant for the ions of the system and needs to be calculated only at the first call.
        lnGaCinf = 0.0_wp
        B(nnp1:NK) = 5.0_wp*QS(nnp1:NK)*log(QS(nnp1:NK)/1.40_wp*0.92_wp/RS(nnp1:NK)) +XL(nnp1:NK) -RS(nnp1:NK)/0.92_wp*(-2.32_wp)
        lnGaCinf(nnp1:NK) = log(RS(nnp1:NK)/0.92_wp) +B(nnp1:NK)
        gcombrefresh = .false.
    endif

    !--------------------------------------------------------
    !Viscosity of mixture calculation; combinatorial part.
    if (calcviscosity) then
        phi = X*ViV !volume based component fraction
        QiQ = QS(1:NK)/QSS
        phiQ = X*QiQ
        !undefined ln_eta0 would contain negative values of large magnitude
        if (all(ln_eta0(1:nneutral) > -7777.0_wp)) then  !pure-component values available (replaced NK with nneutral (JL november 6)
            XieC(1:nneutral) = (exp(lnGaC(1:nneutral))*X(1:nneutral))*ln_eta0(1:nneutral)     
            if (isPEGsystem) then 
                !with PEG, limit combinatorial activity coeff. to between 0 and 100 for viscosity contribution:
                XieC(1:nneutral) = min(exp(lnGaC(1:nneutral)), 1.0E2_wp)*X(1:nneutral)*ln_eta0(1:nneutral)
            endif
        else
            XieC = -55555.5_wp
        endif
    else !a default value indicating "no viscosity calc."
        XieC = -99999.9_wp
    endif
    !--------------------------------------------------------

    end subroutine SRgcomb
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
    subroutine SRsystm(NK)

    use ModSystemProp, only : errorflagmix

    implicit none

    integer,intent(in) :: NK
    integer :: J, JJ, I, J1, L, L1 !, unitm
    real(wp),parameter :: deps = epsilon(1.0_wp)
    real(wp) :: Roxyethylene
    !.....................end of variable declarations.............................

    !--------------------------------------------------------------------------------------------------------
    if (allocated(parA)) then
        deallocate(parA, parB, parC, PsiT, PsiT_dimflip, SRNY, SRNY_dimflip, Nvis, R, Q, XL, RS, QS, lnGaCinf, lnGaRref, XieRref)
    endif
    allocate( parA(NG,NG), parB(NG,NG), parC(NG,NG), PsiT(NG,NG), PsiT_dimflip(NG,NG), SRNY(NK,NG), SRNY_dimflip(NG,NK), &
        & Nvis(NG,NK), R(NG), Q(NG), XL(NK), RS(NK), QS(NK), lnGaCinf(NK), lnGaRref(NK), XieRref(NK) ) !, XieCref(NK)

    parA = 0.0_wp
    parB = 0.0_wp
    parC = 0.0_wp
    R = 0.0_wp
    Q = 0.0_wp
    XL = 0.0_wp
    RS = 0.0_wp  ! sum of the Van der Waals group volumes --> component molar volume
    QS = 0.0_wp  ! sum of the Van der Waals group surface area --> component molar surface area
    SRNY = 0
    SRNY_dimflip = 0
    Nvis = 0.0_wp
    lastTK = 0.0_wp
    gcombrefresh = .true. !if .true., the infinite dilution reference value for the SR combinatorial part of ions will be calculated in SRgcomb and saved for subsequent use.

    !populate array SRNY for SR part
    do J = 1,NG
        JJ = AllSubs(J)
        R(J) = SR_RR(JJ)
        Q(J) = SR_QQ(JJ)
        SRNY(:,J) = ITABsr(:,JJ)
        if (JJ == 154) then !special PEG-group present in the mixture, so it is a "PEG system"
            Roxyethylene = SR_RR(JJ)
        endif
        SRNY_dimflip(J,:) = SRNY(:,J) !same data but in flipped-dimension storage order for use with better memory alignment in certain equations
    enddo

    !assign interaction parameters between different subgroups from  main group interaction parameter matrix.
    !loop only over non-electrolyte (neutral solvent) groups as the interactions of neutral groups with ions is defined to be always 0.0_wp in SR-part.
    do J = 1,NGN
        J1 = Imaingroup(maingrindexofsubgr(J)) !this maingrindexofsubgr array assigns the index of the main group (of the current mixture) associated with neutral subgroup index j
        do L = 1,NGN
            L1 = Imaingroup(maingrindexofsubgr(L))
            parA(L,J) = ARR(L1,J1)
            parB(L,J) = BRR(L1,J1)
            parC(L,J) = CRR(L1,J1)
            !Apply check to avoid using SR inteaction parameters that are not assigned correctly or were not estimated yet:
            if (ARR(L1,J1) < -888887.0_wp) then  !this parameter has not yet been estimated nor is a fitparameter set, so don't use it!!
                errorflagmix = 14
                write(*,*) ""
                write(*,*) "======================================================="
                write(*,*) "        WARNING from SRsystem:"
                write(*,*) "A neutral main group <-> main group interaction coeff. is "
                write(*,*) "not yet defined or no fit parameter was assigned. "
                write(*,*) "System nd: ", nd
                write(*,*) "L1, J1, ARR(L1,J1): ", L1, J1, ARR(L1,J1)
                write(*,*) "======================================================="
                write(*,*) ""
                read(*,*)               !to read(*,*) and wait for user interaction
            endif
            select case(nd)
            case(500:800, 2000:2500)     !potentially used for 3-parameter AIOMFAC (UNIFAC) part for non-electrolyte systems
                if (BRR(L1,J1) < -888887.0_wp .OR. CRR(L1,J1) < -888887.0_wp) then    !this parameter has not yet been estimated nor is a fitparameter set, so don't use it!!
                    errorflagmix = 15
                    write(*,*) ""
                    write(*,*) "======================================================="
                    write(*,*) "        WARNING from SRsystem:"
                    write(*,*) "A neutral main group <-> main group interaction coeff. is "
                    write(*,*) "not yet defined or no fit parameter was assigned. "
                    write(*,*) "System nd: ", nd
                    write(*,*) "L1, J1: ", L1, J1
                    write(*,*) "BRR(L1,J1), CRR(L1,J1): ", BRR(L1,J1), CRR(L1,J1)
                    write(*,*) "======================================================="
                    write(*,*) ""
                    read(*,*) !to read(*,*) and wait for user interaction
                endif
            end select
        enddo
    enddo

    do I = 1,NK
        do J = 1,NG
            RS(I) = RS(I) +SRNY_dimflip(J,I)*R(J) !RS(I) +SRNY(I,J)*R(J)
            QS(I) = QS(I) +SRNY_dimflip(J,I)*Q(J) !QS(I) +SRNY(I,J)*Q(J)
            if (isPEGsystem) then !PEG-polymer solution; use special subgroups and some scaled parameters
                if (AllSubs(J) == 154 .AND. SRNY_dimflip(J,I) > 0) then  !oxyethylene group in PEG polymer, so use a special parametrisation for RS(I) in those cases (idea adapted from Ninni et al., (1999), but implemented in different form, basically to fit RS(I)/QS(I) of oxyethylene group rather than taking fixed Bondi (1964) values.):
                    RS(I) = RS(I) -SRNY_dimflip(J,I)*R(J) !reset above set value for this group.
                    R(J) = 1.381433_wp   !fitSRparam(211)    
                    RS(I) = RS(I) + SRNY_dimflip(J,I)*R(J)
                    !adjust Q(CH2OCH2) value according to the relation between z, RS(I) and QS(I) described in Vera et al. (1977) for straight chain molecules. The expression below for Q(CH2OCH2) is dependent on the overall chain length of the PEG polymer.
                    QS(I) = QS(I) -SRNY_dimflip(J,I)*Q(J) !reset above set value for this group.
                    Q(J) = 3.0_wp        !fitSRparam(212) !as a test with limiting ratio of Q(J)/R(J) for a large number of monomer units
                    QS(I) = QS(I) + SRNY_dimflip(J,I)*Q(J)  !CH2OCH2[PEG]
                endif
            endif
        enddo
        XL(I) = 5.0_wp*(RS(I)-QS(I)) -RS(I) +1.0_wp  !corresponds to l(i) eq. (3) in UNIFAC Fredenslund et al. (1975) paper
        !...
        !calculate Nvis for viscosity calculations according to Cao et al. (1993), Ind. Eng. Chem. Res., 32, 2088--2092.
        do J = 1,NG
            Nvis(J,I) = Q(J)*(0.5_wp*(QS(I)-RS(I))-0.1_wp*(1.0_wp-RS(I)))
        enddo
    enddo
    
    end subroutine SRsystm
!==========================================================================================================================

end module ModSRunifac