!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module providing access to the pure component viscosity parameter arrays.          *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend, Natalie Gervasi                                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2012/09/07                                                      *
!*   -> latest changes: 2019/08/23                                                      *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  SUBROUTINE DeRieux_Tno_Est                                                      *
!*   -  SUBROUTINE VogelTemp                                                            *
!*                                                                                      *
!****************************************************************************************

MODULE Mod_PureViscosPar

IMPLICIT NONE
!..................................
!public parameter arrays:
REAL(8),DIMENSION(1500,2),PUBLIC :: CorrelTrange    !structure: (lower T limit, upper T limit) in [K]

!================================================================================================================================= 
    CONTAINS
!================================================================================================================================= 
    
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
    PURE SUBROUTINE DeRieux_Tno_Est(ind, expTg, Tg)

    USE ModSubgroupProp, ONLY : O2C_H2C_component

    IMPLICIT NONE
    !..................................
    !interface variables:
    INTEGER(4),INTENT(IN) :: ind  !index no. of component in ITAB
    REAL(8),INTENT(IN) :: expTg   !the experimental Tg (if available)
    REAL(8),INTENT(OUT) :: Tg     !glass transition temp [K]
    !local variables:
    REAL(8) :: nCO, bC, bH, bCH, bO, bCO, sumC, sumH, sumO, OtoC, HtoC
    !..................................
    
    CALL O2C_H2C_component(ind, sumC, sumH, sumO, OtoC, HtoC)

    ! Tg is either calculated or an experimentally determined value is used
    IF (expTg > -999.0D0) THEN ! use experimental value
        Tg = expTg
    ELSE
        !choose the constants for the Tg model based on if the compound has oxygens or not
        IF (sumO < 1.0D0) THEN
            nCO = 1.96D0
            bC = 61.99D0
            bH = -113.33D0
            bCH = 28.74D0
            Tg = bC*(nCO + LOG(sumC)) + bH*LOG(sumH) + bCH*LOG(sumC)*LOG(sumH)
        ELSE
            nCO = 12.13D0
            bC = 10.95D0
            bH = -41.82D0
            bCH = 21.61D0
            bO = 118.96D0
            bCO = -24.38D0
            ! Shiraiwa et al. group model used to calculate Tg (DeRieux et al. 2018), eqn (2)
            Tg = bC*(nCO + LOG(sumC)) + bH*LOG(sumH) + bCH*LOG(sumC)*LOG(sumH) + bO*LOG(sumO) + bCO*LOG(sumC)*LOG(sumO)
        ENDIF
    ENDIF

    END SUBROUTINE DeRieux_Tno_Est
    !================================================================================================================================= 
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Calculation of Tno in [K] using the Vogel-Tamman-Fulcher equation in the form      *
    !*   used by DeRieux et al. 2018. Fragility (D) is determined by the temperature of the *
    !*   run.                                                                               *
    !*                                                                                      *
    !*   :: Authors & Copyright ::                                                          *
    !*   Natalie Gervasi, Andi Zuend,                                                       *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *                                  
    !*                                                                                      *
    !****************************************************************************************
    PURE SUBROUTINE VogelTemp(Tg, TempK, D, Tno)

    IMPLICIT NONE

    !interface
    REAL(8), INTENT(IN) :: Tg, TempK ! glass transition temperature, run temperature [K]
    REAL(8), INTENT(OUT) :: Tno, D   ! vogel temperature, fragility
    !........................  

    IF (TempK < Tg) THEN ! pure component substance below Tg that was once fragile will behave more strongly
        D = 30.0D0
    ELSE
        D = 10.0D0 ! reasonable guess for organics (Shiraiwa et al., 2017)
    ENDIF
    Tno = (39.17D0*Tg)/(D + 39.17D0) ! (DeRieux et al. 2018) eqn (7)

    END SUBROUTINE VogelTemp
    !=================================================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Calculation of pure component dynamic viscosity (eta) in Pascal seconds [Pa.s] at  *
    !*   a given temperature T [K].                                                         *
    !*   [cP] = 0.01 [Poise] = 0.01 [g/(cm.s)] = 0.001 [Pa.s] = 0.001 [N.s/m^2]             *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend, Natalie Gervasi                                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2012/09/07                                                      *
    !*   -> latest changes: 2019/09/24                                                      *
    !*                                                                                      *
    !****************************************************************************************
    SUBROUTINE PureCompViscosity(ind, TempK, eta, iflag, Tglass, fragility)
    
    USE ModSystemProp, ONLY : CompN, ITAB
    
    IMPLICIT NONE
    !..................................
    !interface variables:
    INTEGER(4),INTENT(IN) :: ind
    REAL(8),INTENT(IN) :: TempK
    REAL(8),INTENT(OUT) :: Tglass, eta, fragility
    INTEGER(4),INTENT(OUT) :: iflag
    !local variables:
    INTEGER(4) :: cpn, equationNo
    REAL(8) :: a, b, c, d, e, Tg, Tvog
    !..................................

    !initialize:
    IF (ind > 0) THEN
        IF (ITAB(ind,16) > 0) THEN !is water
            cpn = 401
        ELSE
            cpn = compN(ind) !organic or water component ID (when defined)
        ENDIF
    ELSE IF (ind == -1) THEN
        cpn = 401 !water
    ENDIF
    
    SELECT CASE(cpn)
    CASE(401) !water
        equationNo = 12
        a = 225.66D0  
        b = 1.3788D-4
        c = -1.6433D0
        d = 0.0D0
        e = 0.0D0
        CorrelTrange(cpn,1) = 230.0D0
        CorrelTrange(cpn,2) = 495.0D0
    CASE DEFAULT !(1500, 9999) !for system input of any organics from file
        equationNo = 10
        a = -999.0D0
        b = -999.0D0
        c = -999.0D0
        d = -999.0D0
        e = -999.0D0
        CorrelTrange(cpn,1) = 1.0D0
        CorrelTrange(cpn,2) = 1000.0D0
    END SELECT
    
    !compare requested temperature with temperature range of parameterization:
    IF (TempK >= CorrelTrange(cpn,1) .AND. TempK <= CorrelTrange(cpn,2)) THEN !valid
        iflag = 0 !0 means no errors
        ! decide apppropriate value in Tg range
        IF (ind > 0) THEN
            CALL DeRieux_Tno_Est(ind, a, Tg)
            Tglass = Tg
            CALL VogelTemp(Tglass, TempK, fragility, Tvog)
        ELSE IF (ind == -1) THEN !use water properties for ions
            Tg = 225.66D0
            Tglass = Tg
        ENDIF
        !compute eta with given equation for the component
        SELECT CASE(equationNo)
        CASE(12)
            eta = b*((TempK/a) - 1)**c
        CASE(10) !Vogel-Fulcher-Tammann (VFT) using DeRieux et al. (2018) constants and DeRieux Tg
            eta = 10.0D0**(-5.0D0 + 0.434D0*(fragility*Tvog/(TempK - Tvog))) !DeRieux et al. 2018, eqn (6)
        END SELECT
        iflag = 0 !0 means no errors
    ELSE
        iflag = 1 !1 = outside valid temperature range!
    ENDIF

    END SUBROUTINE PureCompViscosity
    !=================================================================================================================================
    
END MODULE Mod_PureViscosPar