!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module providing access to the pure component viscosity parameter arrays.          *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend, Natalie Gervasi                                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2012-09-07                                                      *
!*   -> latest changes: 2021-11-29                                                      *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  subroutine DeRieux_Tno_Est                                                      *
!*   -  subroutine VogelTemp                                                            *
!*   -  subroutine PureCompViscosity                                                    *
!*                                                                                      *
!****************************************************************************************

module ModPureCompViscos

use Mod_kind_param, only : wp

implicit none
!..................................
!public parameter arrays:
real(wp),dimension(1500,2),public :: CorrelTrange    !structure: (lower T limit, upper T limit) in [K]
real(wp),dimension(1500,5),public :: ViscosCorrelPar !structure: (component no., param. A, B, C, D, E)
!public procedures:
public :: PureCompViscosity
private     !as default

!========================================================================================================== 
    contains
!========================================================================================================== 
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Calculation of Tg in [K] using the model from DeRieux et al. 2018 (ACP).           *
    !*   [cP] = 0.01 [Poise] = 0.01 [g/(cm.s)] = 0.001 [Pa.s] = 0.001 [N.s/m^2]             *
    !*                                                                                      *
    !*   :: Authors & Copyright ::                                                          *
    !*   Natalie Gervasi, Andi Zuend,                                                       *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *                         
    !*                                                                                      *
    !****************************************************************************************
    pure subroutine DeRieux_Tno_Est(ind, expTg, Tg)

    use ModSubgroupProp, only : O2C_H2C_component

    implicit none
    !..................................
    !interface variables:
    integer,intent(in) :: ind           !index no. of component in ITAB
    real(wp),intent(in) :: expTg        !the experimental Tg (if available)
    real(wp),intent(out) :: Tg          !glass transition temp [K]
    !local variables:
    real(wp),parameter :: nCO = 12.13_wp, bC = 10.95_wp, bH = -41.82_wp, bCH = 21.61_wp, bO = 118.96_wp, bCO = -24.38_wp
    real(wp),parameter :: nCO_k = 1.96_wp, bC_k = 61.99_wp, bH_k = -113.33_wp, bCH_k = 28.74_wp
    real(wp) :: sumC, sumH, sumO, OtoC, HtoC
    !..................................
    
    call O2C_H2C_component(ind, sumC, sumH, sumO, OtoC, HtoC)

    ! Tg is either calculated or an experimentally determined value is used
    if (expTg > -999.0_wp) then     !use experimental value
        Tg = expTg
    else
        !choose the constants for the Tg model based on whether the compound contains oxygen or not
        if (sumO < 1.0_wp) then
            Tg = bC_k*(nCO_k + log(sumC)) + bH_k*log(sumH) + bCH_k*log(sumC)*log(sumH)
        else
            ! Shiraiwa et al. group model used to calculate Tg (DeRieux et al. 2018), eqn (2)
            Tg = bC*(nCO + log(sumC)) + bH*log(sumH) + bCH*log(sumC)*log(sumH) + bO*log(sumO) + bCO*log(sumC)*log(sumO)
        endif
    endif

    end subroutine DeRieux_Tno_Est
    !========================================================================================================== 
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Calculation of Tno in [K] using the Vogel-Tamman-Fulcher (VTF) equation in the     *
    !*   form used by Angell (2002) and DeRieux et al. 2018.                                *
    !*   Fragility (D) is determined by the temperature of the run.                         *
    !*                                                                                      *
    !*   :: Authors & Copyright ::                                                          *
    !*   Natalie Gervasi, Andi Zuend,                                                       *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *                                  
    !*                                                                                      *
    !****************************************************************************************
    pure subroutine VogelTemp(Tg, TempK, D, Tno)

    implicit none

    !interface
    real(wp), intent(in) :: Tg, TempK   ! glass transition temperature, run temperature [K]
    real(wp), intent(out) :: Tno, D     ! Vogel temperature [K], fragility parameter [-]
    !........................  

    if (TempK < Tg) then                ! pure component substance below Tg that was once fragile will behave more strongly
        D = 30.0_wp
    else
        D = 10.0_wp                     ! reasonable guess for organics (Shiraiwa et al., 2017)
    endif
    Tno = (39.17_wp*Tg)/(D + 39.17_wp)  ! (DeRieux et al. 2018) eqn (7)

    end subroutine VogelTemp
    !==========================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Calculation of pure-component dynamic viscosity (eta0) in Pascal seconds [Pa.s]    *
    !*   at a given temperature T [K]. Output is in form of ln(eta0/[Pa s]).                *
    !*   [cP] = 0.01 [Poise] = 0.01 [g/(cm.s)] = 0.001 [Pa.s] = 0.001 [N.s/m^2]             *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend, Natalie Gervasi                                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2012-09-07                                                      *
    !*   -> latest changes: 2021-11-29                                                      *
    !*                                                                                      *
    !****************************************************************************************
    subroutine PureCompViscosity(ind, TempK, ln_eta0, iflag, Tglass, fragility)
    
    use ModSystemProp, only : CompN, ITAB
    
    implicit none
    !..................................
    !interface variables:
    integer,intent(in) :: ind
    real(wp),intent(in) :: TempK
    real(wp),intent(out) :: Tglass, ln_eta0, fragility
    integer,intent(out) :: iflag
    !local variables:
    integer :: cpn, equationNo
    real(wp),parameter :: ln10 = log(10.0_wp), ln_bwater = log(1.3788E-4_wp)
    real(wp) :: a, b, c, d, e, ln_b, Tg, Tvog
    !..................................

    !initialize:
    if (ind > 0) then
        if (ITAB(ind,16) > 0) then      !is water
            cpn = 401
        else
            cpn = compN(ind)            !organic or water component ID (when defined)
        endif
    else if (ind == -1) then
        cpn = 401                       !water
    endif

    select case(cpn)
    case(401) !water
        if (TempK >= 230.0_wp) then
            equationNo = 12
            a = 225.66_wp
            b = 1.3788E-4_wp
            ln_b = ln_bwater
            c = -1.6433_wp
            d = 0.0_wp
            e = 0.0_wp
            CorrelTrange(cpn,1) = 230.0_wp
            CorrelTrange(cpn,2) = 495.0_wp
        else
            equationNo = 10
            a = 136.0_wp                 !Tg of water from experiments (see, e.g. Koop et al., 2011, PCCP);
            b = -999.0_wp
            c = -999.0_wp
            d = -999.0_wp
            e = -999.0_wp
            CorrelTrange(cpn,1) = 136.0_wp
            CorrelTrange(cpn,2) = 230.0_wp
        endif
    case(402) !CO2(aq)
        equationNo = 10
        a = 130.0_wp
        b = -999.0_wp
        c = -999.0_wp
        d = -999.0_wp
        e = -999.0_wp
        CorrelTrange(cpn,1) = 130.0_wp
        CorrelTrange(cpn,2) = 1000.0_wp
    case default !(1500, 9999) !for system input of organics from file
        cpn = 1500
        equationNo = 10
        a = -999.0_wp
        b = -999.0_wp
        c = -999.0_wp
        d = -999.0_wp
        e = -999.0_wp
        CorrelTrange(cpn,1) = 1.0_wp
        CorrelTrange(cpn,2) = 1000.0_wp
    end select
    
    !for use in AIOMFAC-web instead of the if ... endif block above
    !compare requested temperature with temperature range of parameterization:
    if (TempK >= CorrelTrange(cpn,1) .AND. TempK <= CorrelTrange(cpn,2)) then   !valid
        iflag = 0           !0 means no errors
        !determine apppropriate value in Tg range
        call DeRieux_Tno_Est(ind, a, Tg)
        Tglass = Tg
        call VogelTemp(Tglass, TempK, fragility, Tvog)
        !compute eta with given equation for the component
        select case(equationNo)
        case(12)
            ln_eta0 = ln_b +c*log((TempK/a) - 1.0_wp)
        case(10)            !Vogel-Tammann-Fulcher (VFT), Angell (1991) using DeRieux et al. (2018) constants and DeRieux Tg
            if (Tvog >= TempK) then
                iflag = 1           !1 = outside valid temperature range!
                ln_eta0 = 1.0E300_wp
            else
                ln_eta0 = ln10*( -5.0_wp + 0.434_wp*(fragility*Tvog/(TempK - Tvog)) )    !DeRieux et al. 2018, eqn (6)
            endif
        end select
    else
        iflag = 1           !1 = outside valid temperature range!
    endif

    end subroutine PureCompViscosity
    !==========================================================================================================

end module ModPureCompViscos