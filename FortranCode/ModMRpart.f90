!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module containing data and subroutines to define middle range interaction coeff.   * 
!*   overall and to set the specific interaction coeff. of the current mixture.         *
!*   Subroutine LR_MR_activity is provided to compute LR and MR activity coefficient    *
!*   contributions.                                                                     *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2004                                                            *
!*   -> latest changes: 2021-12-07                                                      *
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
!*   -  SUBROUTINE LR_MR_activity                                                       *
!*   -  SUBROUTINE MRinteractcoeff                                                      *
!*   -  SUBROUTINE MRdata                                                               *
!*   -  SUBROUTINE GammaCO2                                                             *
!*                                                                                      *
!****************************************************************************************
    
MODULE ModMRpart

USE ModSystemProp, ONLY : NGI, NG, topsubno, frominpfile, Nmaingroups, bicarbsyst

IMPLICIT NONE
!....................................................................................
!public variables:
REAL(8),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: Ksp
!** private module variables, which need to be set to public during parameter fitting...
REAL(8),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: bTABna, bTABnc, cTABna, cTABnc
REAL(8),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: Cn1TABAC, Cn2TABAC, omega2TAB, bTABAC, cTABAC
REAL(8),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: TABhighestWTF, TABKsp, omegaTAB
!private module variables:
INTEGER(4),PRIVATE :: jj
REAL(8),PRIVATE :: A_DebyeHw, B_DebyeHw  !the temperature-dependent Debye-Hueckel parameters for water as solvent
REAL(8),DIMENSION(:,:),ALLOCATABLE,PRIVATE :: BAC, CAC, omega, omega2, Cnac1, Cnac2, Raa, Rcc
REAL(8),DIMENSION(:,:),ALLOCATABLE,PRIVATE :: bnc, cnc, bna, cna, omegaNC, omegaNA
REAL(8),DIMENSION(:,:),ALLOCATABLE,PRIVATE :: RccTAB
REAL(8),DIMENSION(:,:,:),ALLOCATABLE,PRIVATE :: Qcca
REAL(8),DIMENSION(:,:,:),ALLOCATABLE,PRIVATE :: qcca1TAB  !REAL(8),DIMENSION(40,40,40),PRIVATE :: qcca1TAB
LOGICAL(4),PRIVATE :: QccaInteract, RccInteract
!set CO2(aq) <--> ion parameters:
REAL(8),DIMENSION(201:topsubno),PARAMETER,PRIVATE :: lambdaIN = [ &
    !Li+         Na+        K+        NH4+        H+, ...
& 0.0000D0, 0.1000D0, 0.0510D0, 0.01000D0, 0.0000D0, 0.000D0, 0.0000D0, 0.0000D0, 0.0000D0, 0.0000D0, & !(201:210)
& 0.0000D0, 0.0000D0, 0.0000D0, 0.00000D0, 0.0000D0, 0.000D0, 0.0000D0, 0.0000D0, 0.0000D0, 0.0000D0, & !(211:220)
& 0.1830D0, 0.0000D0, 0.1830D0, 0.00000D0, 0.0000D0, 0.000D0, 0.0000D0, 0.0000D0, 0.0000D0, 0.0000D0, & !(221:230)
& 0.0000D0, 0.0000D0, 0.0000D0, 0.00000D0, 0.0000D0, 0.000D0, 0.0000D0, 0.0000D0, 0.0000D0, 0.0000D0, & !(231:240)
& 0.0000D0, -0.005D0, 0.0000D0, 0.0000D0, -0.0457D0, 0.000D0, 0.0000D0, -0.003D0, 0.0000D0, 0.0000D0, & !(241:250)
& 0.0000D0, 0.0000D0, 0.0000D0, 0.00000D0, 0.0000D0, 0.000D0, 0.0000D0, 0.0000D0, 0.0000D0, 0.0000D0, & !(251:260)
& 0.0970D0, (0.0D0, jj = 262,topsubno) ] 
!....................................................................................

!$OMP THREADPRIVATE(A_DebyeHw, B_DebyeHw, BAC, CAC, omega, omega2, CnAC1, Cnac2, Raa, Rcc, &
    !$OMP & Ksp, bnc, cnc, bna, cna, omegaNC, omegaNA, Qcca, QccaInteract, RccInteract)
    
!==========================================================================================================================
    CONTAINS
!==========================================================================================================================
  
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to calculate the long range (LR) and middle range (MR) contributions to * 
    !*   the activity coefficients in ion-containing solutions.                             *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
    !*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA                *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2004                                                            *
    !*   -> latest changes: 2018/05/26                                                      *
    !*                                                                                      *
    !****************************************************************************************
    SUBROUTINE LR_MR_activity()

    USE ModSystemProp, ONLY : anionZ, cationZ, Imaingroup, ITAB, ITABMG, maingrindexofsubgr, Mmass, &
        & nelectrol, NGI, NGN, nneutral, solvmixrefnd, SolvSubs, SubGroupMW, errorflagcalc
    USE ModAIOMFACvar, ONLY : DebyeHrefresh, galrln, gamrln, gclrln, gcmrln, gnlrln, gnmrln,        &
        & Ionicstrength, meanSolventMW, SumIonMolalities, SMA, SMC, solvmixcorrMRa, solvmixcorrMRc, &
        & T_K, Tmolal, TmolalSolvmix, XN

    IMPLICIT NONE
    !Local Variables
    INTEGER(4) :: I, II, J, K
    REAL(8) :: A, b, bb, SI, SI2, Z, Test11, Test12, TermI, Sumkion1,  &
    & sumca1, sumq1, sumkion2, sumbsca, sumcnca, sumcsca, sumqc1, SumRc, sumxg, SumQa2, & 
    & lrw1, twoZ, ZSI, avgmaingrMW, Mavg
    REAL(8),PARAMETER :: dtiny = EPSILON(1.0D0)
    REAL(8),PARAMETER :: Mwater = 0.01801528D0  !the molar mass of water [kg/mol]
    REAL(8),PARAMETER :: cDens = 997.0D0        !density of water (kg/m^3) at 298.15 K
    REAL(8),PARAMETER :: densmix = cDens        !here by design the density of reference solvent water at 298.15 K is used
    REAL(8),PARAMETER :: sqrtDensmix = SQRT(densmix)
    REAL(8),PARAMETER :: cDiel = 78.54D0        !static dielectric const. (= relative static permittivity) of water at 298.15 K
    REAL(8),PARAMETER :: dielmix = cDiel
    REAL(8),DIMENSION(NGI) :: SumABca, SumBcx, SumBax, SumCBca, ZA, ZA2, ZC2, smcsma, oexp1, oexp2
    REAL(8),DIMENSION(NGN) :: subgroupxsf, mgroupxsf, mgroupMK, SumBM, gkmrln
    REAL(8),DIMENSION(NGI,NGI) :: BSca, CSca, Bca, Cnca
    REAL(8),DIMENSION(NGN,NGI) :: BKionc, BKiona, BSKionc, BSKiona
    !....................................................................

    !Calculation of the salt-free mole fractions of the neutral subgroups (combined contributions from different molecules if the same subgroup is present):
    DO I = 1,NGN
        II = SolvSubs(I)
        subgroupxsf(I) = SUM(ITAB(1:nneutral,II)*XN(1:nneutral))
    ENDDO
    !normalize salt-free basis subgroup mole fractions:
    sumxg = SUM(subgroupxsf)
    IF (sumxg > 0.0D0) THEN
        subgroupxsf = subgroupxsf/sumxg
    ENDIF

    !calculate the (salt-free basis) mole fractions, mgroupxsf, of the organic main groups present:
    mgroupxsf = 0.0D0
    DO j = 1,NGN !j-loop to cover subgroups
        I = maingrindexofsubgr(j) !this maingrindexofsubgr array assigns the index of the main group (of the current mixture) associated with subgroup index j
        mgroupxsf(I) = mgroupxsf(I) +subgroupxsf(j)
    ENDDO

    !Assignment of the molecular masses to subgroups:
    !mgroupMK: lists the molecular masses of the main groups (of Imaingroup) present in the mixture (as calculated based on the current contributions from different subgroups). 
    !calculate mole-fractional contributions of different subgroups to a certain main group for calculation of main group molar mass:
    mgroupMK = 0.0D0
    DO j = 1,NGN !j-loop to cover subgroups
        I = maingrindexofsubgr(j)
        IF (mgroupxsf(I) > 0.0D0) THEN !check for potential divison by zero problem
            mgroupMK(I) = mgroupMK(I) + SubGroupMW(j)*subgroupxsf(j)/mgroupxsf(I) !normalized contributions of subgroups to main group molar mass [kg/mol]
        ENDIF
    ENDDO
    Mavg = SUM(SubGroupMW(1:NGN)*subgroupxsf(1:NGN))

    !Calculation of the density of the salt free mixture
    !**--++1 water LR 
    !densmix = cDens !SUM(compVF(1:nneutral)*cDens(1:nneutral))
    !!Calculation of the dielectric constant of the salt free mixture:
    !dielmix = cDiel != SUM(compVF(1:nneutral)*cDiel(1:nneutral))
    IF (DebyeHrefresh) THEN !only calculate the first time and when T_K changes significantly
        A_DebyeHw = 1.327757D5*sqrtDensmix/((dielmix*T_K)**1.5D0)
        B_DebyeHw = 6.359696D0*sqrtDensmix/((dielmix*T_K)**0.5D0)
    ENDIF
    A = A_DebyeHw
    b = B_DebyeHw

    !Calculation of the ionic strength (SI): (the Ionic strength in the molality base)
    ZA = ABS(anionZ(1:NGI))
    ZA2 = ZA**2
    ZC2 = cationZ(1:NGI)**2
    SI = 0.5D0*SUM( SMA(1:NGI)*ZA2 + SMC(1:NGI)*ZC2 )
    SI2 = SQRT(SI)
    Ionicstrength = SI
    !Calculation of Z (number of charges per kg solvent):
    Z = SUM( SMA(1:NGI)*ZA(1:NGI) + SMC(1:NGI)*cationZ(1:NGI) )
    !for debugging tests:
    Test11 = SUM(SMA(1:NGI)*ZA(1:NGI))                      !Test of the electrical charge neutrality condition in the solution:  !test test
    Test12 = SUM(SMC(1:NGI)*cationZ(1:NGI))
    Test11 = Test11 - Test12
    IF (.NOT. ANY([bicarbsyst])) THEN
        IF (SI < 1.0D6 .AND. ABS(Test11) > 1.0D-9) THEN     !test
            errorflagcalc = 12                              !overall charge neutrality violation error
            IF (.NOT. frominpfile) THEN
                WRITE(*,*) "WARNING: Electrical charge neutrality condition in the solution is violated! ", Test11
            ENDIF
        ENDIF
    ENDIF

    !:: LR :: calculation of the long-range component of the activity coefficients of the neutrals:
    bb = 1.0D0+b*SI2
    TermI = bb-1.0D0/(bb)-2.0D0*LOG(bb) 
    !all LR properties set to those of pure water (except for the molar mass):
    lrw1 = 2.0D0*A/(b*b*b)*TermI
    gnlrln(1:nneutral) = Mmass(1:nneutral)*lrw1

    !:: LR :: calculation of the long-range component of the activity coefficients of the ions:
    lrw1 = A*SI2/bb
    gclrln = -ZC2*lrw1
    galrln = -ZA2*lrw1

    !################################################################################################################
    ! :: MR :: calculation of the middle-range contribution of the activity coefficients of the ions and neutrals:
    !################################################################################################################
    IF (nelectrol > 0 .AND. SI > dtiny) THEN !**--++
        SumCBca = 0.0D0
        SumABca = 0.0D0
        !Calculation of Bca, Cnca, Bkionc and Bkiona:
        lrw1 = -0.5D0/SI2
        DO J = 1,NGI
            oexp1(1:NGI) = omega(1:NGI,J)*SI2 
            WHERE (oexp1(1:NGI) > 300.0D0) !prevent floating point underflow as the second term becomes practically zero:
                Bca(1:NGI,J) = bac(1:NGI,J) !+cac(1:NGI,J)*EXP(-300.0D0)
            ELSEWHERE
                Bca(1:NGI,J) = bac(1:NGI,J) +cac(1:NGI,J)*EXP(-oexp1(1:NGI))
            ENDWHERE
            oexp2(1:NGI) = omega2(1:NGI,J)*SI2
            WHERE (oexp2(1:NGI) > 300.0D0) !prevent floating point underflow as the second term becomes practically zero:
                Cnca(1:NGI,J) = cnac1(1:NGI,J) !+cnac2(1:NGI,J)*EXP(-300.0D0)
            ELSEWHERE
                Cnca(1:NGI,J) = cnac1(1:NGI,J) +cnac2(1:NGI,J)*EXP(-oexp2(1:NGI))
            ENDWHERE
            ! Calculation of derivatives of Bca, Cnca with respect to the ionic strength SI. BSca: derivative of Bca.
            BSca(1:NGI,J) = lrw1*omega(1:NGI,J)*(Bca(1:NGI,J)-bac(1:NGI,J)) 
            CSca(1:NGI,J) = lrw1*omega2(1:NGI,J)*(Cnca(1:NGI,J)-cnac1(1:NGI,J))
        ENDDO

        DO J = 1,NGI
            IF (SI2 > 250.0D0) THEN !prevent floating point underflow as the second term becomes practically zero:
                Bkionc(1:NGN,J) = bnc(1:NGN,J)
                Bkiona(1:NGN,J) = bna(1:NGN,J)
            ELSE
                Bkionc(1:NGN,J) = bnc(1:NGN,J) +cnc(1:NGN,J)*EXP(-omegaNC(1:NGN,J)*SI2)
                Bkiona(1:NGN,J) = bna(1:NGN,J) +cna(1:NGN,J)*EXP(-omegaNA(1:NGN,J)*SI2)
            ENDIF
            !Calculation of derivatives of Bkionc and Bkiona with respect to the ionic strength SI: BSkionc: derivative of Bkionc; BSkiona: derivative of Bkiona
            BSkionc(1:NGN,J) = lrw1*omegaNC(1:NGN,J)*(Bkionc(1:NGN,J)-bnc(1:NGN,J))  != (-0.5D0*omegaNC(I,J)/SI2)*cnc(I,J)*EXP(-omegaNC(I,J)*SI2)
            BSkiona(1:NGN,J) = lrw1*omegaNA(1:NGN,J)*(Bkiona(1:NGN,J)-bna(1:NGN,J))  != (-0.5D0*omegaNA(I,J)/SI2)*cna(I,J)*EXP(-omegaNA(I,J)*SI2)
        ENDDO

        !the mean molecular mass of a main group in the solvent mixture.
        avgmaingrMW = SUM(mgroupMK*mgroupxsf) !this is the mean main group molar mass -- and not the mean solvent molar mass!

        !Calculation of the midrange part for the neutrals:
        DO I = 1,NGN
            SumBm(I) = SUM(Bkiona(I,1:NGI)*SMA(1:NGI) +Bkionc(I,1:NGI)*SMC(1:NGI))
        ENDDO
        Sumkion1 = 0.0D0
        DO J = 1,NGI
            Sumkion1 = Sumkion1 +SUM( (Bkiona(1:NGN,J)+SI*BSkiona(1:NGN,J))*mgroupxsf(1:NGN)*SMA(J) &
                                   & +(Bkionc(1:NGN,J)+SI*BSkionc(1:NGN,J))*mgroupxsf(1:NGN)*SMC(J) )
        ENDDO
        SumCA1 = 0.0D0
        twoZ = 2.0D0*Z
        ZSI = Z*SI
        DO J = 1,NGI
            SumCA1 = SumCA1 +SUM((Bca(1:NGI,J)+SI*BSca(1:NGI,J)+twoZ*Cnca(1:NGI,J)+ZSI*Csca(1:NGI,J))*SMC(1:NGI)*SMA(J))
        ENDDO

        !Calculation of the three ion interactions in the neutral part:      
        SumQ1 = 0.0D0
        IF (QccaInteract) THEN !so far only used for [NH4+|H+|HSO4-] containing mixtures
            DO I = 1,NGI
                DO J = I,NGI
                    SumQ1 = SumQ1 +SUM(2.0D0*Qcca(I,J,1:NGI)*SMC(I)*SMC(J)*SMA(1:NGI)) !term when ionic strength dependency is parameterized: SumQ1+(2.0D0*Qcca(I,J,K)+SI*QScca(I,J,K))*SMC(I)*SMC(J)*SMA(K)
                ENDDO
            ENDDO
        ENDIF
        ! Calculation of the two cation interactions in the neutral part:
        SumRc = 0.0D0
        IF (RccInteract) THEN !so far only used for [NH4+|H+|HSO4-] containing mixtures
            DO I = 1,NGI
                SumRc = SumRc +SUM(Rcc(I,I:NGI)*SMC(I)*SMC(I:NGI)) 
            ENDDO
        ENDIF

        !gkmrln(I): ln of the middle-range part of the activity coefficient of main group I
        gkmrln(1:NGN) = SumBm(1:NGN) -mgroupMK(1:NGN)/avgmaingrMW*Sumkion1 -mgroupMK(1:NGN)*SumCA1 &
            & -mgroupMK(1:NGN)*SumQ1 -mgroupMK(1:NGN)*SumRc

        !Calculation of the middle-range part of the activity coefficients of the neutral molecules from the contributions of main groups to individual molecules
        !gnmrln(I): ln of the middle-range part of the activity coefficient of neutral molecule I
        gnmrln = 0.0D0
        DO J = 1,NGN
            II = Imaingroup(J)
            IF (II > 0) THEN
                DO I = 1,nneutral
                    IF (ITABMG(I,II) > 0) THEN
                        DO K = 1,NGN
                            IF (Imaingroup(K) == II) THEN
                                gnmrln(I) = gnmrln(I) +ITABMG(I,II)*gkmrln(K)
                            ENDIF
                        ENDDO
                    ENDIF
                ENDDO
            ENDIF
        ENDDO

        !Calculation of the middle-range contribution part for the cations and anions: @@##
        SumBcx = 0.0D0
        SumBax = 0.0D0
        DO I = 1,NGI
            SumBCx(I) = SUM(Bkionc(1:NGN,I)*mgroupxsf(1:NGN)) 
            SumBax(I) = SUM(Bkiona(1:NGN,I)*mgroupxsf(1:NGN)) 
        ENDDO
        SumBsCA = 0.0D0
        SumCnca = 0.0D0
        SumCsca = 0.0D0
        DO J = 1,NGI
            SumCBca(1:NGI) = SumCBca(1:NGI) +(BCA(1:NGI,J) +Z*Cnca(1:NGI,J))*SMA(J)
            SumABca(1:NGI) = SumABca(1:NGI) +(BCA(J,1:NGI) +Z*Cnca(J,1:NGI))*SMC(J)
            smcsma(1:NGI) = SMC(1:NGI)*SMA(J)
            SumBsCA = sumBsCA +SUM(BSca(1:NGI,J)*smcsma(1:NGI))
            SumCnca = sumCnca +SUM(Cnca(1:NGI,J)*smcsma(1:NGI))
            SumCsca = sumCsca +SUM(Csca(1:NGI,J)*smcsma(1:NGI))
        ENDDO
        Sumkion2 = 0.0D0
        DO J = 1,NGI
            Sumkion2 = Sumkion2+SUM(BSkiona(1:NGN,J)*mgroupxsf(1:NGN)*SMA(J)+BSkionc(1:NGN,J)*mgroupxsf(1:NGN)*SMC(J))
        ENDDO

        !gcmrln(I): ln of the middle-range part of the activity coefficient of cation I
        !gamrln(I): ln of the middle-range part of the activity coefficient of anion I
        DO I = 1,NGI
            !cations:
            gcmrln(I) = sumbcx(I)/avgmaingrMW +ZC2(I)/(2.0D0*avgmaingrMW)*sumkion2 +sumCBca(I) &
                & +0.5D0*ZC2(I)*SumBsCA +cationZ(I)*SumCnca +0.5D0*ZC2(I)*Z*SumCsca
            IF (QccaInteract) THEN !so far only used for [NH4+|H+|HSO4-] containing mixtures
                SumQc1 = 0.0D0
                DO K = 1,NGI
                    SumQc1 = SumQc1 +SUM(SMC(1:NGI)*SMA(K)*Qcca(I,1:NGI,K))  !entspr. Qc1,c,a(I)
                ENDDO
                gcmrln(I) = gcmrln(I) +SumQc1 !the 3 ion interaction terms... 
            ENDIF
            !there might be an Interaction between two cations which is important for the system-fit (NH4+ <-> H+):
            IF (RccInteract) THEN !so far only used for [NH4+|H+] containing mixtures
                gcmrln(I) = gcmrln(I) +SUM(Rcc(I,1:NGI)*SMC(1:NGI))
            ENDIF 
            !-- anions:
            gamrln(I) = SumBax(I)/avgmaingrMW +ZA2(I)/(2.0D0*avgmaingrMW)*sumkion2 +sumABca(I) &
                & +0.5D0*ZA2(I)*SumBsCA+ZA(I)*SumCnca +0.5D0*ZA2(I)*Z*SumCsca
            IF (QccaInteract) THEN !so far only used for [NH4+|H+|HSO4-] containing mixtures
                SumQa2 = 0.0D0
                DO J = 1,NGI
                    SumQa2 = SumQa2 +SUM(SMC(J)*SMC(J:NGI)*Qcca(J,J:NGI,I)) 
                ENDDO
                gamrln(I) = gamrln(I)+SumQa2
            ENDIF
        ENDDO

        !Calculation of the term that is needed for the ions to convert from mole fraction to molality units.
        !The reference solvent is water for MR and SR and solvent mixture for LR (but LR is anyways set to water properties):
        Tmolal = LOG(Mwater/meanSolventMW +Mwater*SumIonMolalities)  !this is the molality basis conversion term with reference pure water.
        !----
        IF (solvmixrefnd) THEN  !calculate additional terms required to correct for different reference state in AIOMFAC_calc
            TmolalSolvmix = LOG(1.0D0 +SumIonMolalities*meanSolventMW)  !here reference state is the salt-free solvent mixture, because the LR-term is always with reference solvent mixture (unfortunately)!!
            solvmixcorrMRc = 0.0D0
            solvmixcorrMRa = 0.0D0
            DO I = 1,NGI
                solvmixcorrMRc(I) = solvmixcorrMRc(I)+SUM( (bnc(1:NGN,I) +cnc(1:NGN,I))*mgroupxsf(1:NGN) )
                solvmixcorrMRa(I) = solvmixcorrMRa(I)+SUM( (bna(1:NGN,I) +cna(1:NGN,I))*mgroupxsf(1:NGN) )
            ENDDO
            solvmixcorrMRc(1:NGI) = -solvmixcorrMRc(1:NGI)/meanSolventMW
            solvmixcorrMRa(1:NGI) = -solvmixcorrMRa(1:NGI)/meanSolventMW
        ENDIF
        !----
    ELSE  !**--++ (no ions present in mixture)
        gamrln = 0.0D0
        gcmrln = 0.0D0
        gnmrln = 0.0D0
        Tmolal = 0.0D0
        IF (solvmixrefnd) THEN
            TmolalSolvmix = 0.0D0
            solvmixcorrMRc = 0.0D0
            solvmixcorrMRa = 0.0D0
        ENDIF
    ENDIF !**--++

    END SUBROUTINE LR_MR_activity
!========================================================================================================================== 
      
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Middle-range (MR) part of the AIOMFAC model: set neutral group <--> ion and        * 
    !*   ion <--> ion interaction coeff. for current mixture.                               *
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
    SUBROUTINE MRinteractcoeff()
   
    !Public variables:  
    USE ModSystemProp, ONLY : NGN, NGI, NG, Ication, Ianion, Imaingroup, Ncation, isPEGsystem, &
        & CatNr, AnNr, errorflagmix, frominpfile 
    
    IMPLICIT NONE
    !Local variables:
    INTEGER(4) :: I, J, NN, NC, NA, NC2, K
    !.................................................................................
    
    !allocate public arrays for use in MR calculation part:
    IF (ALLOCATED(Ksp)) THEN
        DEALLOCATE(Ksp, BAC, CAC, omega, omega2, Cnac1, Cnac2, bnc, cnc, bna, cna, omegaNC, omegaNA)
    ENDIF
    ALLOCATE( Ksp(NGI,NGI), BAC(NGI,NGI), CAC(NGI,NGI), omega(NGI,NGI), omega2(NGI,NGI), Cnac1(NGI,NGI),  &
        & Cnac2(NGI,NGI), bnc(NG,NGI), cnc(NG,NGI), bna(NG,NGI), cna(NG,NGI), omegaNC(NG,NGI), omegaNA(NG,NGI) )
    
    !----------------------------------------------------------------   
    ! ASSIGNMENT OF THE MR-INTERACTION COEFFICIENTS
    !---------------------------------------------------------------- 
    !initialize neutral <-> ions interaction parameters:
    BNC = 0.0D0
    BNA = 0.0D0
    CNC = 0.0D0
    CNA = 0.0D0
    omegaNC = 1.2D0
    omegaNA = 1.2D0
      
    !calculate the neutral main group <-> ion interaction coefficients:        
    DO I = 1,NGN
        NN = Imaingroup(I)
        IF (NN > 0) THEN
            DO J = 1,NGI
                IF (Ication(J) > 0) THEN
                    NC = Ication(J)-200
                    BNC(I,J) = bTABnc(NN,NC)
                    CNC(I,J) = cTABnc(NN,NC)                !CNC(I,J) is CNC(main group I, cation J), thus interaction coeff. C: main group I <-> cation J
                    IF (isPEGsystem) THEN                   !use special parameters for PEG systems
                        IF (NN == 68 .AND. NC == 4) THEN    !CHn[OH,PEG] <--> NH4+
                            BNC(I,J) = 3.057678302653D-01   !fitparam(273) !bTABnc(NN,NC)
                            CNC(I,J) = 4.041113321682D-01   !fitparam(274) !cTABnc(NN,NC) 
                        ENDIF
                    ENDIF
                    !check to avoid using inteaction parameters that are not assigned correctly or were not estimated yet:
                    IF (BNC(I,J) < -1000.0D0 .OR. CNC(I,J) < -1000.0D0) THEN !this parameter is not yet estimated nor is a fitparameter set
                        errorflagmix = 1
                        IF (.NOT. frominpfile) THEN
                            !$OMP CRITICAL (printingtoscreen2)
                            WRITE(*,*) ""
                            WRITE(*,*) "======================================================="
                            WRITE(*,*) "        WARNING from MRinteractcoeff:"
                            WRITE(*,*) "This neutral main group <-> cation interaction coeff. is "
                            WRITE(*,*) "not yet defined or no fit parameter was assigned. "
                            WRITE(*,*) "NN, NC, bTABnc(NN,NC), cTABnc(NN,NC): ", NN, NC, bTABnc(NN,NC), cTABnc(NN,NC)
                            WRITE(*,*) "======================================================="
                            WRITE(*,*) ""
                            READ(*,*)
                            !$OMP END CRITICAL (printingtoscreen2)
                        ENDIF
                    ENDIF
                END IF
                IF (Ianion(J) > 0) THEN
                    NA = Ianion(J)-240
                    BNA(I,J) = bTABna(NN,NA)
                    CNA(I,J) = cTABna(NN,NA)
                    IF (isPEGsystem) THEN !use special parameters for PEG systems
                        IF (NN == 68 .AND. NA == 21) THEN !CHn[OH,PEG] <--> SO4--
                            BNA(I,J) = -2.1822703311941D-01 !fitparam(275) !bTABna(NN,NA)
                            CNA(I,J) = -1.2595743942958D-01 !fitparam(276) !cTABna(NN,NA) 
                        ENDIF
                    ENDIF
                    !check to avoid using inteaction parameters that are not assigned correctly or were not estimated yet:
                    IF (BNA(I,J) < -1000.0D0 .OR. CNA(I,J) < -1000.0D0) THEN !this parameter is not yet estimated nor is a fitparameter set
                        errorflagmix = 2
                        IF (.NOT. frominpfile) THEN
                            !$OMP CRITICAL (printingtoscreen3)
                            WRITE(*,*) ""
                            WRITE(*,*) "======================================================="
                            WRITE(*,*) "        WARNING from MRinteractcoeff:"
                            WRITE(*,*) "This neutral main group <-> anion interaction coeff. is "
                            WRITE(*,*) "not yet defined or no fit parameter was assigned. "
                            WRITE(*,*) "NN, NA, bTABna(NN,NA), cTABna(NN,NA): ", NN, NA, bTABna(NN,NA), cTABna(NN,NA)
                            WRITE(*,*) "======================================================="
                            WRITE(*,*) ""
                            READ(*,*)
                            !$OMP END CRITICAL (printingtoscreen3)
                        ENDIF
                    ENDIF
                END IF
            ENDDO
        END IF
    ENDDO
    !---------------------------------------------------------------- 

    !----------------------------------------------------------------  
    !initialize cation <-> anion interaction coefficients and input/output variables:
    BAC = 0.0D0
    CAC = 0.0D0
    CnAC1 = 0.0D0
    CnAC2 = 0.0D0
    Ksp = 0.0D0     !solubility coeff. @ 298 K
    omega = 0.8D0   !the same for all, will be overwritten in some specific cases
    omega2 = 0.6D0

    !----------------------------------------------------------------       
    !calculate the cation <-> anion interaction coefficients: 
    DO I = 1,NGI
        NC = Ication(I)-200
        DO J = 1,NGI
            NA = Ianion(J)-240
            IF (NC >= 1 .AND. NA >= 1) THEN
                BAC(I,J) = bTABAC(NC,NA)
                CAC(I,J) = cTABAC(NC,NA)
                CnAC1(I,J) = cn1TABAC(NC,NA)
                CnAC2(I,J) = cn2TABAC(NC,NA)
                omega2(I,J) = omega2TAB(NC,NA)
                omega(I,J) = omegaTAB(NC,NA)
                Ksp(I,J) = TABKsp(NC,NA)
                IF (ABS(bTABAC(NC,NA)) < 1.0D-10) THEN
                    IF (.NOT. (NC == 5 .AND. NA == 7)) THEN  !exclude H+ <--> OH- case
                        errorflagmix = 9
                        IF (.NOT. frominpfile) THEN
                            !$OMP CRITICAL (MRI1)
                            WRITE(*,*) "WARNING from MRinteractcoeff: interaction parameters of a cation--anion pair"
                            WRITE(*,*) "present in the mixture is not defined: cation-ID, anion-ID: ", NC+200, NA+240
                            READ(*,*)
                            !$OMP END CRITICAL (MRI1)
                        ENDIF
                    ENDIF
                ENDIF
            ENDIF
        ENDDO
    ENDDO

    !---------------------------------------------------------------- 
    !The interactions of three different ions are described with Qcca and Qcaa:
    !There is also an interaction parameter between NH4+ and H+ ions declared here:
    IF (ALLOCATED(Qcca)) THEN
        DEALLOCATE(Qcca, Rcc, Raa)
    ENDIF
    ALLOCATE( Qcca(NGI,NGI,NGI), Rcc(NGI,NGI), Raa(NGI,NGI) )
    Qcca = 0.0D0
    Rcc = 0.0D0
    Raa = 0.0D0

    !---------------------------------------------------------------- 
    !Build up the Qcca and Rcc1 for the actual composition:
    DO I = 1,NGI
        NC = Ication(I)-200
        DO J = 1,NGI
            NC2 = Ication(J)-200
            IF (NC >= 1 .AND. NC2 >= 1) THEN
                Rcc(I,J) = RccTAB(NC,NC2)
                DO K = 1,NGI
                    NA = Ianion(K)-240
                    IF (NA >= 1) THEN
                        Qcca(I,J,K) = qcca1TAB(NC,NC2,NA)
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    !---------------------------------------------------------------- 

    !set logical switches for use in LR_MR_activity:
    QccaInteract = .false.
    RccInteract = .false.
    IF (Ncation > 1) THEN  !there might be an Interaction between two cations, which is important for the system-fit (NH4+ <-> H+): 
        IF (CatNr(204) > 0 .AND. CatNr(205) > 0) THEN !so far only used for [NH4+|H+|HSO4-] containing mixtures
            RccInteract = .true.
            IF (AnNr(248) > 0) THEN
                QccaInteract = .true.
            ENDIF
        ENDIF
    ENDIF

    END SUBROUTINE MRinteractcoeff
!==========================================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Middle-range (MR) part of the AIOMFAC model: load data tables for organic groups   *
    !*   and ions and initialize/populate interaction coefficient data tables.              *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
    !*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA                *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2004                                                            *
    !*   -> latest changes: 2021-09-17                                                      *
    !*                                                                                      *
    !****************************************************************************************
    SUBROUTINE MRdata()

    IMPLICIT NONE
    !................................
    
    ALLOCATE( bTABnc(Nmaingroups, 23), cTABnc(Nmaingroups, 23), &
        & bTABna(Nmaingroups, topsubno-240), cTABna(Nmaingroups, topsubno-240) )
    
    !.... initialize main group <--> ion interaction parameter data ...................
    !main group <--> cations
    bTABnc = -7.777778D+04      !default value indicating that the parameter has not been determined
    cTABnc = -7.777778D+04
    !load data for main group <--> cation interactions (only for the existing cation IDs, e.g. cation 01 is Li+, 02 is Na+):
    !               main group_no         1           2                3             4               5              6             7              8            9              10              11            12              13           14              15          16               17              18              19              20           21              22             23           24              25            26             27              28             29             30             31             32             33            34             35             36             37             38            39              40             41             42             43             44             45             46             47             48             49             50             51             52             53             54             55             56             57             58             59             60             61             62             63             64             65             66             67             68             69             70             71             72             73             74             75                76
    !               main_groups  (CHn),            (C=C),         (ACHn),    (ACCHn),   (OH[stand_UNIFAC]), (CH3OH[methanol]), (H2O),        (ACOH),       (CHnCO),    (CHO[aldehyde]),   (CCOO),  (HCOO[formate]),  (CHnO[ether]),   (CNH2),        (CNH),       ((C)3N),         (ACNH2),     (PYRIDINE), (CCN[nitrile]), (COOH[stand_UNIFAC]), (CCl),     (CCl2),          (CCl3),         (CCl4),       (ACCl),         (CNO2),        (ACNO2),       (CS2),      (CH3SH),    (FURFURAL),         (DOH),          (I),           (BR),   (C-=C[triple_bond]), (DMSO),       (ACRY),        (ClCC),        (ACF),        (DMF),          (CF2),          (COO),          (SIH2),      (SIO),          (NMP),         (CClF),       (CON),       (OCCOH),           (CH2S),   (MORPHOLINE),   (THIOPHENE), (Inorg_Ions), (CHn[MingR_lo]), (OH[MingR_l]),(COOH[MingR_l]),(CHnCO[MingR_l]),(CHn[MingR_mo]),(CHnO[MingR_mo]),(OH[MingR_mo]),(CHn[MingR_hy]),(COOH[MingR_hy]),(OH[MingR_hy]),(CHn[MingR_di]),(COOH[MingR_di]),(OH[Peng]),(COOH),  (CHn[alc]),  (CHn[alc-tail]),  (CHn[OH]),       (OH),    (CH2OCH2[PEG]),  (CHnONO2),  (CHnOOH[perox]),(C(=O)OOH[perox]),(CHnOOCHm[perox]),(C(=O)OONO2[perox]), (CO2)  
    bTABnc(1:Nmaingroups,01) = [  7.961623D-02,  1.287780D-01,  6.208020D-02, -7.777778D+04,  5.877733D-03, -7.777778D+04,  0.000000D+00,  4.459015D-02,  1.117971D-01,  5.189833D-02,  1.117970D-01, -7.777778D+04,  6.128712D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  4.458747D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  4.458747D-02,  5.877733D-03, -7.777778D+04,  4.458747D-02,  5.877733D-03,  4.458747D-02,  7.961623D-02,  7.961623D-02,  6.188986D-02,  5.877733D-03, -7.777778D+04,  1.242037D-01,  6.716485D-02,  4.458747D-02,  1.225742D-01, -7.777778D+04,  0.000000D+00]
    cTABnc(1:Nmaingroups,01) = [  4.669577D-02,  1.362762D-01,  3.087262D-02, -7.777778D+04, -5.629984D-02, -7.777778D+04,  0.000000D+00, -6.965679D-03,  1.422573D-01, -1.682426D-03,  1.422573D-01, -7.777778D+04,  2.696984D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -6.968065D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -6.968065D-03, -5.629984D-02, -7.777778D+04, -6.968065D-03, -5.629984D-02, -6.968065D-03,  4.669577D-02,  4.669577D-02,  2.583107D-02, -5.629984D-02, -7.777778D+04,  3.972771D-02, -2.933000D-02, -6.968065D-03,  5.393967D-02, -7.777778D+04,  0.000000D+00]
    bTABnc(1:Nmaingroups,02) = [  1.058813D-01,  2.187202D-01,  1.058813D-01, -7.777778D+04,  7.649999D-03, -7.777778D+04,  0.000000D+00,  7.501697D-02,  1.628598D-01,  1.057398D-01,  1.628598D-01, -7.777778D+04,  1.018669D-01, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.556653D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.556653D-02,  7.649999D-03, -7.777778D+04,  5.556653D-02,  7.649999D-03,  5.556653D-02,  1.058813D-01,  1.058813D-01,  1.058805D-01,  7.649999D-03, -7.777778D+04,  1.614479D-01,  1.095169D-01,  5.556653D-02,  2.037338D-01, -7.777778D+04,  0.000000D+00]
    cTABnc(1:Nmaingroups,02) = [  2.266823D-02,  5.099094D-02,  2.266823D-02, -7.777778D+04, -7.485817D-03, -7.777778D+04,  0.000000D+00,  6.854487D-03,  5.045448D-02,  2.272655D-02,  5.045446D-02, -7.777778D+04,  1.485420D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -6.696041D-04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -6.696041D-04, -7.485817D-03, -7.777778D+04, -6.696041D-04, -7.485817D-03, -6.696041D-04,  2.266823D-02,  2.266823D-02,  9.730931D-03, -7.485817D-03, -7.777778D+04,  2.199863D-02,  7.368382D-03, -6.696041D-04,  2.970840D-02, -7.777778D+04,  0.000000D+00]
    bTABnc(1:Nmaingroups,03) = [  9.836417D-02,  1.593944D-01,  7.030872D-02, -7.777778D+04,  2.704778D-03, -7.777778D+04,  0.000000D+00,  4.748628D-02,  1.335596D-01,  4.748629D-02,  1.018410D-01, -7.777778D+04,  6.474752D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.240469D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.240469D-02,  2.704778D-03, -7.777778D+04,  3.240469D-02,  2.704778D-03,  3.240469D-02,  9.836417D-02,  9.836417D-02,  6.847925D-02,  2.704778D-03, -7.777778D+04,  1.307689D-01,  6.745230D-02,  3.240469D-02,  1.294950D-01, -7.777778D+04,  0.000000D+00]
    cTABnc(1:Nmaingroups,03) = [  4.203279D-02,  6.730633D-02,  3.245205D-02, -7.777778D+04, -3.237396D-02, -7.777778D+04,  0.000000D+00,  1.421022D-02,  5.039454D-02,  1.421023D-02, -4.268114D-03, -7.777778D+04,  6.472387D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  1.995250D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  1.995250D-02, -3.237396D-02, -7.777778D+04,  1.995250D-02, -3.237396D-02,  1.995250D-02,  4.203279D-02,  4.203279D-02,  3.207157D-02, -3.237396D-02, -7.777778D+04,  6.198529D-02,  3.234991D-02,  1.995250D-02,  1.294477D-01, -7.777778D+04,  0.000000D+00]
    bTABnc(1:Nmaingroups,04) = [  5.667445D-02,  9.381847D-02,  5.584332D-02, -7.777778D+04,  5.238242D-03, -7.777778D+04,  0.000000D+00,  4.060326D-02,  9.371399D-02,  5.622542D-02,  9.371399D-02, -7.777778D+04,  4.060311D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.534892D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.534892D-02,  5.238242D-03, -7.777778D+04, -2.534892D-02,  5.238242D-03, -2.534892D-02,  5.667445D-02,  5.667445D-02,  3.904649D-02,  5.238242D-03,  2.249583D-01,  3.132553D-02,  4.584135D-02, -2.534892D-02,  8.120622D-02, -7.777778D+04,  0.000000D+00]
    cTABnc(1:Nmaingroups,04) = [  4.316537D-02,  1.165552D-01,  4.299857D-02, -7.777778D+04,  5.684108D-03, -7.777778D+04,  0.000000D+00,  2.120821D-02,  5.169221D-02,  3.702557D-02,  5.169221D-02, -7.777778D+04,  4.505894D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -3.149018D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -3.149018D-02,  5.684108D-03, -7.777778D+04, -3.149018D-02,  5.684108D-03, -3.149018D-02,  4.316537D-02,  4.316537D-02,  4.444826D-02,  5.684108D-03,  3.997699D-01,  1.167519D-02,  5.074304D-02, -3.149018D-02,  9.011787D-02, -7.777778D+04,  0.000000D+00]
    bTABnc(1:Nmaingroups,05) = [  9.447874D-02,  1.981737D-01,  9.068728D-02, -7.777778D+04, -2.496397D-02, -7.777778D+04,  0.000000D+00,  5.500988D-02,  1.968493D-01,  1.023700D-01,  1.968490D-01, -7.777778D+04, -7.280049D-01, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -1.541675D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -1.541675D-02, -2.496397D-02, -7.777778D+04, -1.541675D-02, -2.496397D-02, -1.541675D-02,  9.447874D-02,  9.447874D-02,  7.969963D-02, -2.496397D-02, -7.777778D+04,  7.906199D-02, -7.529688D-01, -1.541675D-02, -1.456010D+00, -7.777778D+04,  0.000000D+00]
    cTABnc(1:Nmaingroups,05) = [  6.568965D-02,  1.381514D-01,  2.230661D-02, -7.777778D+04,  7.786872D-05, -7.777778D+04,  0.000000D+00,  1.413451D-03, -4.275766D-03, -6.996550D-02, -4.275770D-03, -7.777778D+04, -2.625517D-01, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.472152D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.472152D-02,  7.786872D-05, -7.777778D+04, -7.472152D-02,  7.786872D-05, -7.472152D-02,  6.568965D-02,  6.568965D-02,  3.134212D-02,  7.786872D-05, -7.777778D+04, -9.031866D-03, -2.624738D-01, -7.472152D-02, -5.251033D-01, -7.777778D+04,  0.000000D+00]
    bTABnc(1:Nmaingroups,21) = [  1.021409D-01,  2.062126D-01,  8.109440D-02, -7.777778D+04,  1.038726D-02, -7.777778D+04,  0.000000D+00,  7.100129D-02,  2.044581D-01,  1.021191D-01,  2.017796D-01, -7.777778D+04,  8.254596D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.605447D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.605447D-02,  1.038726D-02, -7.777778D+04,  2.605447D-02,  1.038726D-02,  2.605447D-02,  1.021409D-01,  1.021409D-01,  8.092567D-02,  1.038726D-02, -7.777778D+04,  1.281954D-01,  9.293322D-02,  2.605447D-02,  1.650919D-01, -7.777778D+04,  0.000000D+00]
    cTABnc(1:Nmaingroups,21) = [  7.474215D-02,  1.525628D-01,  5.931454D-02, -7.777778D+04, -5.465988D-04, -7.777778D+04,  0.000000D+00,  2.885781D-02,  1.524886D-01,  7.471416D-02, -1.077497D-01, -7.777778D+04,  7.077540D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.331719D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.331719D-02, -5.465988D-04, -7.777778D+04,  2.331719D-02, -5.465988D-04,  2.331719D-02,  7.474215D-02,  7.474215D-02,  5.926138D-02, -5.465988D-04, -7.777778D+04,  9.805934D-02,  7.022880D-02,  2.331719D-02,  1.415508D-01, -7.777778D+04,  0.000000D+00]
    bTABnc(1:Nmaingroups,23) = [  8.761896D-02,  1.465395D-01,  7.522760D-02, -7.777778D+04,  3.103042D-03, -7.777778D+04,  0.000000D+00,  6.831109D-02,  1.442190D-01,  6.831109D-02,  1.442190D-01, -7.777778D+04,  7.252455D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.112358D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.112358D-02,  3.103042D-03, -7.777778D+04,  5.112358D-02,  3.103042D-03,  5.112358D-02,  8.761896D-02,  8.761896D-02,  6.749006D-02,  3.103042D-03, -7.777778D+04,  1.387425D-01,  7.562759D-02,  5.112358D-02,  1.450491D-01, -7.777778D+04,  0.000000D+00]
    cTABnc(1:Nmaingroups,23) = [ -2.317547D-02, -1.395783D-02, -2.883615D-02, -7.777778D+04, -8.659318D-03, -7.777778D+04,  0.000000D+00, -4.267910D-02, -2.095079D-02, -4.267910D-02, -2.095079D-02, -7.777778D+04, -2.344034D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.495400D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.495400D-02, -8.659318D-03, -7.777778D+04, -2.495400D-02, -8.659318D-03, -2.495400D-02, -2.317547D-02, -2.317547D-02, -2.232628D-02, -8.659318D-03, -7.777778D+04, -4.812947D-02, -3.209966D-02, -2.495400D-02, -4.688069D-02, -7.777778D+04,  0.000000D+00]
  
    
    !main group <--> anions:
    bTABna = -7.777778D+04      !default value indicating that the parameter has not been determined
    cTABna = -7.777778D+04
    !               main group_no         1           2                3             4               5              6             7              8            9              10              11            12              13           14              15          16               17              18              19              20           21              22             23           24              25            26             27              28             29             30             31             32             33            34             35             36             37             38            39              40             41             42             43             44             45             46             47             48             49             50             51             52             53             54             55             56             57             58             59             60             61             62             63             64             65             66             67             68             69             70             71             72             73             74             75                76
    !               main_groups     (CHn),        (C=C),         (ACHn),    (ACCHn),   (OH[stand_UNIFAC]), (CH3OH[methanol]), (H2O),        (ACOH),       (CHnCO),    (CHO[aldehyde]),   (CCOO),  (HCOO[formate]),  (CHnO[ether]),   (CNH2),        (CNH),       ((C)3N),         (ACNH2),     (PYRIDINE), (CCN[nitrile]), (COOH[stand_UNIFAC]), (CCl),     (CCl2),          (CCl3),         (CCl4),       (ACCl),         (CNO2),        (ACNO2),       (CS2),      (CH3SH),    (FURFURAL),         (DOH),          (I),           (BR),   (C-=C[triple_bond]), (DMSO),       (ACRY),        (ClCC),        (ACF),        (DMF),          (CF2),          (COO),          (SIH2),      (SIO),          (NMP),         (CClF),       (CON),       (OCCOH),           (CH2S),   (MORPHOLINE),   (THIOPHENE), (Inorg_Ions), (CHn[MingR_lo]), (OH[MingR_l]),(COOH[MingR_l]),(CHnCO[MingR_l]),(CHn[MingR_mo]),(CHnO[MingR_mo]),(OH[MingR_mo]),(CHn[MingR_hy]),(COOH[MingR_hy]),(OH[MingR_hy]),(CHn[MingR_di]),(COOH[MingR_di]),(OH[Peng]),(COOH),  (CHn[alc]),  (CHn[alc-tail]),  (CHn[OH]),       (OH),    (CH2OCH2[PEG]),  (CHnONO2),  (CHnOOH[perox]),(C(=O)OOH[perox]),(CHnOOCHm[perox]),(C(=O)OONO2[perox]), (CO2)
    bTABna(1:Nmaingroups,02) = [  6.914307D-02,  1.329562D-01,  5.290417D-02, -7.777778D+04,  9.156391D-03, -7.777778D+04,  0.000000D+00,  4.028553D-02,  1.329562D-01,  6.904875D-02,  1.327976D-01, -7.777778D+04,  5.425750D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.885736D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.885736D-02,  9.156391D-03, -7.777778D+04,  3.885736D-02,  9.156391D-03,  3.885736D-02,  6.914307D-02,  6.914307D-02,  5.227910D-02,  9.156391D-03, -7.777778D+04,  1.080004D-01,  6.341389D-02,  3.885736D-02,  1.085150D-01, -7.777778D+04,  0.000000D+00]
    cTABna(1:Nmaingroups,02) = [  5.888996D-02,  1.363867D-01,  3.860389D-02, -7.777778D+04, -2.249250D-02, -7.777778D+04,  0.000000D+00,  3.767145D-02,  1.363866D-01,  5.886737D-02,  1.362981D-01, -7.777778D+04,  4.161487D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.009678D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.009678D-02, -2.249250D-02, -7.777778D+04,  3.009678D-02, -2.249250D-02,  3.009678D-02,  5.888996D-02,  5.888996D-02,  3.656014D-02, -2.249250D-02, -7.777778D+04,  8.898674D-02,  1.912237D-02,  3.009678D-02,  8.322975D-02, -7.777778D+04,  0.000000D+00]
    bTABna(1:Nmaingroups,03) = [  4.136791D-02,  8.229121D-02,  4.128208D-02, -7.777778D+04,  4.686682D-03, -7.777778D+04,  0.000000D+00,  2.998976D-02,  6.480397D-02,  3.717467D-02,  6.475417D-02, -7.777778D+04,  4.133287D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.998974D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.998974D-02,  4.686682D-03, -7.777778D+04,  2.998974D-02,  4.686682D-03,  2.998974D-02,  4.136791D-02,  4.136791D-02,  4.094891D-02,  4.686682D-03, -7.777778D+04,  7.135765D-02,  4.601955D-02,  2.998974D-02,  8.266574D-02, -7.777778D+04,  0.000000D+00]
    cTABna(1:Nmaingroups,03) = [  3.439384D-02,  5.607449D-02,  3.401668D-02, -7.777778D+04,  1.476176D-03, -7.777778D+04,  0.000000D+00,  1.237431D-02,  4.151743D-02,  1.807473D-02,  4.150937D-02, -7.777778D+04,  2.391176D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  1.237430D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  1.237430D-02,  1.476176D-03, -7.777778D+04,  1.237430D-02,  1.476176D-03,  1.237430D-02,  3.439384D-02,  3.439384D-02,  1.688331D-02,  1.476176D-03, -7.777778D+04,  4.676815D-02,  2.538794D-02,  1.237430D-02,  4.782352D-02, -7.777778D+04,  0.000000D+00]
    bTABna(1:Nmaingroups,04) = [  4.878929D-02,  7.134251D-02,  4.259999D-02, -7.777778D+04,  1.939024D-02, -7.777778D+04,  0.000000D+00,  4.133419D-02,  7.376797D-02,  3.983930D-02,  3.427790D-02, -7.777778D+04,  4.117653D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.427825D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.235025D-02,  1.939024D-02, -7.777778D+04,  3.235025D-02,  1.939024D-02,  3.235025D-02,  4.878929D-02,  4.878929D-02,  3.427807D-02,  1.939024D-02, -7.777778D+04,  8.306755D-02,  6.056677D-02,  3.427825D-02,  8.235307D-02, -7.777778D+04,  0.000000D+00]
    cTABna(1:Nmaingroups,04) = [ -8.516072D-03, -1.301376D-02, -1.896587D-02, -7.777778D+04, -1.636921D-02, -7.777778D+04,  0.000000D+00, -2.434093D-02, -1.520187D-02, -2.517013D-02,  2.169341D-02, -7.777778D+04, -1.417372D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -3.647636D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -1.342741D-02, -1.636921D-02, -7.777778D+04, -1.342741D-02, -1.636921D-02, -1.342741D-02, -8.516072D-03, -8.516072D-03, -6.082930D-03, -1.636921D-02, -7.777778D+04, -4.499244D-02, -3.054294D-02, -3.647636D-02, -2.834745D-02, -7.777778D+04,  0.000000D+00]
    bTABna(1:Nmaingroups,05) = [  4.233234D-02,  6.701021D-02,  4.209408D-02, -7.777778D+04, -2.707269D-02, -7.777778D+04,  0.000000D+00,  4.535367D-03,  4.542737D-02,  4.139841D-02,  4.542735D-02, -7.777778D+04,  4.773836D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.755684D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.755684D-03, -2.707269D-02, -7.777778D+04, -2.755684D-03, -2.707269D-02, -2.755684D-03,  4.233234D-02,  4.233234D-02,  3.000031D-02, -2.707269D-02, -7.777778D+04,  3.957665D-02, -2.229886D-02, -2.755684D-03,  9.547671D-03, -7.777778D+04,  0.000000D+00]
    cTABna(1:Nmaingroups,05) = [  4.853253D-02,  1.230290D-01,  4.848432D-02, -7.777778D+04,  8.326951D-03, -7.777778D+04,  0.000000D+00,  2.880167D-02,  6.624150D-02,  4.875210D-02,  6.624148D-02, -7.777778D+04,  8.524851D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  6.663517D-04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  6.663517D-04,  8.326951D-03, -7.777778D+04,  6.663517D-04,  8.326951D-03,  6.663517D-04,  4.853253D-02,  4.853253D-02,  3.967023D-02,  8.326951D-03, -7.777778D+04,  4.919888D-02,  9.357546D-02,  6.663517D-04,  1.704970D-01, -7.777778D+04,  0.000000D+00]
    bTABna(1:Nmaingroups,06) = [  4.233234D-02,  6.701021D-02,  4.209408D-02, -7.777778D+04, -2.707269D-02, -7.777778D+04,  0.000000D+00,  4.535367D-03,  4.542737D-02,  4.139841D-02,  4.542735D-02, -7.777778D+04,  4.773836D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.755684D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.755684D-03, -2.707269D-02, -7.777778D+04, -2.755684D-03, -2.707269D-02, -2.755684D-03,  4.233234D-02,  4.233234D-02,  3.000031D-02, -2.707269D-02, -7.777778D+04,  3.957665D-02, -2.229886D-02, -2.755684D-03,  9.547671D-03, -7.777778D+04,  0.000000D+00]
    cTABna(1:Nmaingroups,06) = [  4.853253D-02,  1.230290D-01,  4.848432D-02, -7.777778D+04,  8.326951D-03, -7.777778D+04,  0.000000D+00,  2.880167D-02,  6.624150D-02,  4.875210D-02,  6.624148D-02, -7.777778D+04,  8.524851D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  6.663517D-04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  6.663517D-04,  8.326951D-03, -7.777778D+04,  6.663517D-04,  8.326951D-03,  6.663517D-04,  4.853253D-02,  4.853253D-02,  3.967023D-02,  8.326951D-03, -7.777778D+04,  4.919888D-02,  9.357546D-02,  6.663517D-04,  1.704970D-01, -7.777778D+04,  0.000000D+00]
    bTABna(1:Nmaingroups,07) = [  4.878929D-02,  7.134251D-02,  4.259999D-02, -7.777778D+04,  1.939024D-02, -7.777778D+04,  0.000000D+00,  4.133419D-02,  7.376797D-02,  3.983930D-02,  3.427790D-02, -7.777778D+04,  4.117653D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.427825D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.235025D-02,  1.939024D-02, -7.777778D+04,  3.235025D-02,  1.939024D-02,  3.235025D-02,  4.878929D-02,  4.878929D-02,  3.427807D-02,  1.939024D-02, -7.777778D+04,  8.306755D-02,  6.056677D-02,  3.427825D-02,  8.235307D-02, -7.777778D+04,  0.000000D+00]
    cTABna(1:Nmaingroups,07) = [ -8.516072D-03, -1.301376D-02, -1.896587D-02, -7.777778D+04, -1.636921D-02, -7.777778D+04,  0.000000D+00, -2.434093D-02, -1.520187D-02, -2.517013D-02,  2.169341D-02, -7.777778D+04, -1.417372D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -3.647636D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -1.342741D-02, -1.636921D-02, -7.777778D+04, -1.342741D-02, -1.636921D-02, -1.342741D-02, -8.516072D-03, -8.516072D-03, -6.082930D-03, -1.636921D-02, -7.777778D+04, -4.499244D-02, -3.054294D-02, -3.647636D-02, -2.834745D-02, -7.777778D+04,  0.000000D+00]
    bTABna(1:Nmaingroups,08) = [  1.244526D-01,  2.299447D-01,  9.637726D-02, -7.777778D+04, -6.100228D-02, -7.777778D+04,  0.000000D+00,  4.375717D-02,  1.728215D-01,  7.159263D-02,  1.808274D-01, -7.777778D+04,  6.071562D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02, -6.100228D-02, -7.777778D+04,  7.999459D-02, -6.100228D-02,  7.999459D-02,  1.244526D-01,  1.244526D-01,  9.617478D-02, -6.100228D-02, -7.777778D+04,  2.044472D-01, -2.866620D-04,  7.999459D-02,  1.214312D-01, -7.777778D+04,  0.000000D+00]
    cTABna(1:Nmaingroups,08) = [  1.426790D-01,  1.962500D-01,  9.774756D-02, -7.777778D+04, -3.952179D-02, -7.777778D+04,  0.000000D+00,  3.042277D-03,  1.053198D-01,  2.249248D-02,  8.572575D-02, -7.777778D+04,  2.439568D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06, -3.952179D-02, -7.777778D+04,  5.620295D-06, -3.952179D-02,  5.620295D-06,  1.426790D-01,  1.426790D-01,  9.757109D-02, -3.952179D-02, -7.777778D+04,  1.426847D-01, -1.512611D-02,  5.620295D-06,  4.879137D-02, -7.777778D+04,  0.000000D+00]
    bTABna(1:Nmaingroups,09) = [  1.244526D-01,  2.299447D-01,  9.637726D-02, -7.777778D+04, -6.100228D-02, -7.777778D+04,  0.000000D+00, -7.777778D+04,  2.333990D-01,  1.089460D-01,  2.333990D-01, -7.777778D+04,  6.071562D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02, -6.100228D-02, -7.777778D+04,  7.999459D-02, -6.100228D-02,  7.999459D-02,  1.244526D-01,  1.244526D-01,  9.617478D-02, -6.100228D-02, -7.777778D+04,  2.044472D-01, -2.866600D-04,  7.999459D-02,  1.214312D-01, -7.777778D+04,  0.000000D+00]
    cTABna(1:Nmaingroups,09) = [  1.426790D-01,  1.962500D-01,  9.774756D-02, -7.777778D+04, -3.952179D-02, -7.777778D+04,  0.000000D+00, -7.777778D+04,  2.234230D-01,  8.074430D-02,  2.234230D-01, -7.777778D+04,  2.439568D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06, -3.952179D-02, -7.777778D+04,  5.620295D-06, -3.952179D-02,  5.620295D-06,  1.426790D-01,  1.426790D-01,  9.757109D-02, -3.952179D-02, -7.777778D+04,  1.426846D-01, -1.512611D-02,  5.620295D-06,  4.879136D-02, -7.777778D+04,  0.000000D+00]
    bTABna(1:Nmaingroups,10) = [  1.244526D-01,  2.299447D-01,  9.637726D-02, -7.777778D+04, -6.100228D-02, -7.777778D+04,  0.000000D+00,  4.375717D-02,  1.728215D-01,  7.159263D-02,  1.808274D-01, -7.777778D+04,  6.071562D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02, -6.100228D-02, -7.777778D+04,  7.999459D-02, -6.100228D-02,  7.999459D-02,  1.244526D-01,  1.244526D-01,  9.617478D-02, -6.100228D-02, -7.777778D+04,  2.044472D-01, -2.866620D-04,  7.999459D-02,  1.214312D-01, -7.777778D+04,  0.000000D+00]
    cTABna(1:Nmaingroups,10) = [  1.426790D-01,  1.962500D-01,  9.774756D-02, -7.777778D+04, -3.952179D-02, -7.777778D+04,  0.000000D+00,  3.042277D-03,  1.053198D-01,  2.249248D-02,  8.572575D-02, -7.777778D+04,  2.439568D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06, -3.952179D-02, -7.777778D+04,  5.620295D-06, -3.952179D-02,  5.620295D-06,  1.426790D-01,  1.426790D-01,  9.757109D-02, -3.952179D-02, -7.777778D+04,  1.426847D-01, -1.512611D-02,  5.620295D-06,  4.879137D-02, -7.777778D+04,  0.000000D+00]
    bTABna(1:Nmaingroups,21) = [  7.598433D-02,  1.241119D-01,  7.056449D-02, -7.777778D+04,  2.271222D-03, -7.777778D+04,  0.000000D+00,  3.663588D-02,  1.238882D-01,  3.904856D-02,  1.726755D-02, -7.777778D+04,  4.690451D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -4.482039D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -4.482039D-02,  2.271222D-03, -7.777778D+04, -4.482039D-02,  2.271222D-03, -4.482039D-02,  7.598433D-02,  7.598433D-02,  5.257356D-02,  2.271222D-03, -3.803980D-01,  3.116394D-02,  4.917573D-02, -4.482039D-02,  9.380902D-02, -7.777778D+04,  0.000000D+00]
    cTABna(1:Nmaingroups,21) = [  1.130960D-01,  2.717248D-01,  6.332340D-02, -7.777778D+04,  1.224391D-02, -7.777778D+04,  0.000000D+00,  6.143489D-02,  1.346749D-01,  7.196915D-02,  1.077373D-03, -7.777778D+04,  8.606047D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -5.001442D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -5.001442D-02,  1.224391D-02, -7.777778D+04, -5.001442D-02,  1.224391D-02, -5.001442D-02,  1.130960D-01,  1.130960D-01,  7.931526D-02,  1.224391D-02, -3.296446D-01,  6.308154D-02,  9.830437D-02, -5.001442D-02,  1.721209D-01, -7.777778D+04,  0.000000D+00]
    bTABna(1:Nmaingroups,22) = [  7.598433D-02,  1.241119D-01,  7.056449D-02, -7.777778D+04,  2.271222D-03, -7.777778D+04,  0.000000D+00,  3.663588D-02,  1.238882D-01,  3.904856D-02,  1.726755D-02, -7.777778D+04,  4.690451D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -4.482039D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -4.482039D-02,  2.271222D-03, -7.777778D+04, -4.482039D-02,  2.271222D-03, -4.482039D-02,  7.598433D-02,  7.598433D-02,  5.257356D-02,  2.271222D-03, -3.803980D-01,  3.116394D-02,  4.917573D-02, -4.482039D-02,  9.380902D-02, -7.777778D+04,  0.000000D+00]
    cTABna(1:Nmaingroups,22) = [  1.130960D-01,  2.717248D-01,  6.332340D-02, -7.777778D+04,  1.224391D-02, -7.777778D+04,  0.000000D+00,  6.143489D-02,  1.346749D-01,  7.196915D-02,  1.077373D-03, -7.777778D+04,  8.606047D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -5.001442D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -5.001442D-02,  1.224391D-02, -7.777778D+04, -5.001442D-02,  1.224391D-02, -5.001442D-02,  1.130960D-01,  1.130960D-01,  7.931526D-02,  1.224391D-02, -3.296446D-01,  6.308154D-02,  9.830437D-02, -5.001442D-02,  1.721209D-01, -7.777778D+04,  0.000000D+00]
    
    !------ cation <-->  anion interactions (in water) ----------------------
    !bkion and ckion interaction coefficients estimated by AIOMFAC model fit to experimental data:
    !further interaction parameters are set in MRfitpar, but only during fit of parameters.  

    ! Data for the ionic part of the midrange interaction term
    ! Determined by AZ in 2006, 2007 (4 - 5 parameter fit with omega = 0.8, omega2 = variable).
    ! bTABAC(I,J): contains the b values of the midrange terms. I: cations, J: anions.
    ! cTABAC(I,J): contains the c values of the midrange terms
    ! bAC(I,J): interaction parameter b between cation I and anion J.
    ! cAC(I,J): interaction parameter c between cation I and anion J.
    ALLOCATE( Cn1TABAC(40,40), Cn2TABAC(40,40), omega2TAB(40,40), bTABAC(40,40), cTABAC(40,40), &
        & TABhighestWTF(40,40), TABKsp(40,40), omegaTAB(40,40) )
    
    bTABAC = 0.0D0
    cTABAC = 0.0D0
    cn1TABAC = 0.0D0
    cn2TABAC = 0.0D0
    omega2TAB = 0.6D0 !arbitrary value
    omegaTAB = 0.8D0  !the usual value
    TABhighestWTF = 0.0D0
    TABKsp = 1.0D30 !default value on [molal basis] (indicating an extremely high solubility product as default)
  
    !Li+ <-> Cl-          4P AUG 2007      
    bTABAC(1,2) = 0.106554911884921D0        
    cTABAC(1,2) = 0.206369566797954D0
    cn1TABAC(1,2) = 0.0D0
    cn2TABAC(1,2) = 5.323941573688994D-002  
    omega2TAB(1,2) = 0.535548024591189D0     
    TABhighestwtf(1,2) = 0.45D0
    TABKsp(1,2) = 1.4815D6 ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !Li+ <-> Br-          4P AUG 2007      
    bTABAC(1,3) = 0.106384115653176D0 
    cTABAC(1,3) = 0.316479521480175D0    
    cn1TABAC(1,3) = 0.0D0
    cn2TABAC(1,3) = 5.760198325796455D-002    
    omega2TAB(1,3) = 0.464658335659553D0      
    TABhighestwtf(1,3) = 0.35D0
    
    !Li+ <-> I- 
    bTABAC(1,4) =  1.8588934753695888D-01
    cTABAC(1,4) =  6.0808982141071000D-01   
    cn1TABAC(1,4) = 0.0D0
    cn2TABAC(1,4) = 4.1814140083157633D-02
    omega2TAB(1,4) = 4.1341524514949507D-01
    TABKsp(1,4) = 1.107D06 ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !Li+ <-> NO3-         4P JUN 2007      
    bTABAC(1,5) = 7.631270368216216D-002
    cTABAC(1,5) = 0.300549960057599D0
    cn1TABAC(1,5) = 0.0D0
    cn2TABAC(1,5) = 4.670051954857116D-002
    omega2TAB(1,5) = 0.664928462062345
    TABhighestwtf(1,5) = 0.58D0               
    
    !Li+ <-> IO3-     May 2021
    bTABAC(1,6) = 7.6049825043886410D-02
    cTABAC(1,6) = -2.4930506321164203D-01
    omegaTAB(1,6) = 2.9689011671088517D-01
    cn1TABAC(1,6) = 0.0D0
    cn2TABAC(1,6) = -1.2852865367456584D-01
    omega2TAB(1,6) = 1.7423828299201640D+00
    TABhighestwtf(1,6) = 0.16D0 !NOt true value
    
    !Li+ <-> OH-      June 2021
    bTABAC(1,7) = 7.5665164617218825D-02
    cTABAC(1,7) = -9.1330575713229178D-01
    cn1TABAC(1,7) = 0.0D0
    cn2TABAC(1,7) = 6.9105686988206494D-01
    omega2TAB(1,7) =  2.3477855868588415D+00
    
    !Li+ <-> HSO4-    May 2021
    bTABAC(1,8) = 2.8878992634118333D-01
    cTABAC(1,8) = -2.5337206136122727D+00
    cn1TABAC(1,8) = 0.0D0 
    cn2TABAC(1,8) = 6.7541421521779388D-01
    omega2TAB(1,8) = 1.5684706739420728D+00
    
    !Li+ <-> HCO3-   estimated by K+ <-> HCO3-
    bTABAC(1,10) =  9.2903436457389355D-02
    cTABAC(1,10) =  -8.8237465618429356D-01
    omegaTAB(1,10) = 5.6552279057559085D-01
    cn1TABAC(1,10) = 0.0D0
    cn2TABAC(1,10) = 5.9250512478964340D-03
    omega2TAB(1,10) = 6.2492599123964832D-01

    !Li+ <-> SO4--        4P JUN 2007  
    bTABAC(1,21) = 0.114470419384503
    cTABAC(1,21) = 3.540074251345881D-002
    cn1TABAC(1,21) = 0.0D0
    cn2TABAC(1,21) = -0.263257805913055D0 
    omega2TAB(1,21) = 1.31696686093901D0
    TABhighestwtf(1,21) = 0.25D0
    TABKsp(1,21) = 3.25D0 ![molal basis] molal ion activity product at the solubility limit @ 298.15 K  

    !Li+ <-> CO3--        Oct 2020    
    bTABAC(1,22) = 7.5165264810580406D-04
    cTABAC(1,22) = -1.6045235322580549D0
    cn1TABAC(1,22) = 0.0D0
    cn2TABAC(1,22) = -9.5120991084252995D-03
    omega2TAB(1,22) =  3.6525829891663431D-01
    TABKsp(1,22) = 3.493D-03 ![molal basis] molal ion activity product at the solubility limit @ 298.15 K

    !Na+ <-> Cl-          4P JUN 2007         
    bTABAC(2,2) = 5.374079546760321D-002 
    cTABAC(2,2) = 7.977128474886985D-002   
    cn1TABAC(2,2) = 0.0D0
    cn2TABAC(2,2) = 2.455349637351519D-002  
    omega2TAB(2,2) = 0.562981194934520D0       
    TABhighestwtf(2,2) = 0.44D0
    TABKsp(2,2) = 37.83D0 ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !Na+ <-> Br-          4P JUN 2007        
    bTABAC(2,3) = 0.180807320471190D0
    cTABAC(2,3) = 0.273143813534931D0   
    cn1TABAC(2,3) = 0.0D0
    cn2TABAC(2,3) = -0.506597565983690D0 
    omega2TAB(2,3) = 2.20904992991060D0  
    TABhighestwtf(2,3) = 0.49D0 
    TABKsp(2,3) = 322.0D0 ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !Na+ <-> I-            
    bTABAC(2,4) = 6.834973D-02 
    cTABAC(2,4) = 3.796614D-01   
    cn1TABAC(2,4) = 0.0D0
    cn2TABAC(2,4) = 4.145436D-02
    omega2TAB(2,4) = 4.994560D-01  
    TABKsp(2,4) = 5.741D03 ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !Na+ <-> NO3-        5P JUN 2007
    bTABAC(2,5) = 1.164350105903698D-003
    cTABAC(2,5) = -0.102545526130049D0 
    omegaTAB(2,5) = 0.410452675064302D0 
    cn1TABAC(2,5) = 0.0D0
    cn2TABAC(2,5) = 2.534857925520951D-003 
    omega2TAB(2,5) = 0.512657343426316D0   
    TABhighestwtf(2,5) = 0.80D0
    TABKsp(2,5) = 11.82D0 ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !Na+ <-> IO3-     May 2021
    bTABAC(2,6) = 1.3437657607868917D-02
    cTABAC(2,6) = -9.6486635628713002D-01
    cn1TABAC(2,6) = 0.0D0
    cn2TABAC(2,6) = -1.3706943474700559D-01
    omega2TAB(2,6) =  1.3626195116594571D+00
    omegaTAB(2,6) =   1.3047771888155979D+00

    !Na+ <-> OH-      Mar 2020
    bTABAC(2,7) = 2.3840333396438568D-02
    cTABAC(2,7) = 6.6956349269406479D-01
    omegaTAB(2,7) = 2.2686902829375624D-01
    cn1TABAC(2,7) = 0.0D0
    cn2TABAC(2,7) = -5.2304938122407685D-01
    omega2TAB(2,7) = 1.3037546443208794D+00

    !Na+ <-> HSO4-       5P FEB 2011           
    bTABAC(2,8) = 1.5321434237810060D-02
    cTABAC(2,8) = 4.0D-01
    omegaTAB(2,8) = 4.2363461984697831D-01 
    cn1TABAC(2,8) = 0.0D0    
    cn2TABAC(2,8) = 3.5007219661690844D-03   
    omega2TAB(2,8) = 4.0D-01 
    TABhighestwtf(2,8) = 0.15D0   
    
    !Na+ <-> HCO3-  Jan 2021
    bTABAC(2,10) =  6.4998964390762845D-02
    cTABAC(2,10) =  -7.5189119274972557D-02
    cn1TABAC(2,10) = 0.0D0
    cn2TABAC(2,10) = -5.7300345498924488D-02
    omega2TAB(2,10) = 9.3276066630078436D-01
    
    !Na+ <-> SO4--       4P JUN 2007           
    bTABAC(2,21) = 1.890680686105647D-003
    cTABAC(2,21) = -0.424184309855867D0
    cn1TABAC(2,21) = 0.0D0
    cn2TABAC(2,21) = -0.223851238896391D0     
    omega2TAB(2,21) = 1.05361978864913D0  
    TABhighestwtf(2,21) = 0.60D0
    !TABKsp(2,21) = 0.1236D0 !@ 298 K; [molal basis] molal ion activity product at the solubility limit @ 298 K (these salt or hydrates have a strong temperature dependence!).   
        !value based on aqueous Na2SO4 solubility at 293 and 303 K from Okorafor (1999) (actually it forms a hydrate Na2SO4*10 H2O at T < 32 deg C!
    TABKsp(2,21) =  0.4877D0 !@ 308 K; [molal basis] molal ion activity product at the solubility limit @ 308 K, based on data of anhydrous Na2SO4 solubility at 308 K by Vener and Thompson (1949).                
    
    !Na+ <-> CO3--    Oct 2020          
    bTABAC(2,22) =  1.0695833656549962D-01
    cTABAC(2,22) =  -3.8465749714565167D-01
    cn1TABAC(2,22) = 0.0D0
    cn2TABAC(2,22) = -1.9818447835032355D-01
    omega2TAB(2,22) =  8.4635302771576537D-01
    
    !K+ <-> Cl-        4P JUL 2007
    bTABAC(3,2) = 1.656079143002206D-002       
    cTABAC(3,2) = -2.752209560273460D-003
    cn1TABAC(3,2) = 0.0D0
    cn2TABAC(3,2) = 2.083280655191561D-002
    omega2TAB(3,2) = 0.670529802215773D0
    TABhighestwtf(3,2) = 0.49D0
    TABKsp(3,2) = 8.0D0 ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !K+ <-> Br-          4P JUN 2007        
    bTABAC(3,3) = 3.368790641502646D-002 
    cTABAC(3,3) = 6.088180362627697D-002 
    cn1TABAC(3,3) = 0.0D0 
    cn2TABAC(3,3) = 1.529324768456540D-002 
    omega2TAB(3,3) = 0.565063111593592D0
    TABhighestwtf(3,3) = 0.40D0
    TABKsp(3,3) = 11.34D0 ![molal basis] molal ion activity product at the solubility limit @ 298 K
    
    !K+ <-> I-            
    bTABAC(3,4) = 5.442088D-02 
    cTABAC(3,4) = 1.882261D-01
    cn1TABAC(3,4) = 0.0D0
    cn2TABAC(3,4) = 2.534565D-02 
    omega2TAB(3,4) = 6.607733D-01 
    TABKsp(3,4) = 54.11D0 ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !K+ <-> NO3-          5P JUN 2007
    bTABAC(3,5) = 2.498205354734865D-005  
    cTABAC(3,5) =  -0.413171796795632D0   
    cn1TABAC(3,5) = 0.0D0
    cn2TABAC(3,5) = -4.551039248054616D-004 
    omega2TAB(3,5) = 0.242243640032790D0   
    omegaTAB(3,5) = 0.357227451106132D0  
    TABhighestwtf(3,5) = 0.26D0 
    TABKsp(3,5) = 0.9225D0 ![molal basis] molal ion activity product at the solubility limit @ 298 K 
    
    !K+ <-> IO3-    ! April 2021
    bTABAC(3,6) = 1.2531824958095275D-02
    cTABAC(3,6) = -1.8395581736986735D+00
    cn1TABAC(3,6) = 0.0D0
    cn2TABAC(3,6) = -1.4618851585020370D-01
    omega2TAB(3,6) =  1.3612241273448702D+00
    omegaTAB(3,6) = 1.2744128860114006D+00

    !K+ <-> OH-      ! Mar 2020
    bTABAC(3,7) = 2.5131496199052189D-01
    cTABAC(3,7) =  5.3647067717101127D-01
    cn1TABAC(3,7) = 0.0D0
    cn2TABAC(3,7) = -7.5910080641938760D-01
    omega2TAB(3,7) = 1.9363546657176669D+00
    TABKsp(3,7) = 999.0D0 ![not true value]! 
    
    !K+ <-> HSO4-    Revised May 2019
    bTABAC(3,8) = 8.4433898901815352D-02
    cTABAC(3,8) = -3.9301228443995739D-01
    cn1TABAC(3,8) = 0.0D0 
    cn2TABAC(3,8) = -9.0232965273929300D-01
    omega2TAB(3,8) = 1.8692081240631140D+00

    !K+ <-> SO4--     Revised May 2019 
    bTABAC(3,21) =  6.6714599718481538D-02
    cTABAC(3,21) = -1.0844118398860187D+00
    cn1TABAC(3,21) = 0.0D0
    cn2TABAC(3,21) = -6.6511262694399914D-02
    omega2TAB(3,21) = 5.8327209615736464D-01   
    TABhighestwtf(3,21) = 0.11D0  
    TABKsp(3,21) = 1.571D-02 ![molal basis]
    
    !K+ <-> HCO3-   Jan 2021
    bTABAC(3,10) =  9.2903436457389355D-02
    cTABAC(3,10) =  -8.8237465618429356D-01
    omegaTAB(3,10) = 5.6552279057559085D-01
    cn1TABAC(3,10) = 0.0D0
    cn2TABAC(3,10) = 5.9250512478964340D-03
    omega2TAB(3,10) = 6.2492599123964832D-01
    
    !K+ <-> CO3--   Oct 2020
    bTABAC(3,22) =  1.9268078927583168D-01
    cTABAC(3,22) =  1.3547830332027755D-01
    cn1TABAC(3,22) = 0.0D0
    cn2TABAC(3,22) = -1.8190589684462621D-01
    omega2TAB(3,22) = 9.4333743510918933D-01

    !NH4+ <-> Cl-      5P JUN 2007  
    bTABAC(4,2) = 1.520261933644218D-003
    cTABAC(4,2) = 4.907431846455088D-002  
    cn1TABAC(4,2) = 0.0D0    
    cn2TABAC(4,2) = 1.111220056439640D-002
    omega2TAB(4,2) = 0.653256017313056D0  
    TABhighestwtf(4,2) = 0.54D0 
    omegaTAB(4,2) = 0.116800561298130  
    TABKsp(4,2) = 17.2D0 ![molal basis] 

    !NH4+ <-> Br-     5P JUL 2007        
    bTABAC(4,3) = 2.498124477602936D-003
    cTABAC(4,3) = 8.151172908753296D-002 
    cn1TABAC(4,3) = 0.0D0
    cn2TABAC(4,3) = 1.379527431149357D-002  
    omega2TAB(4,3) = 0.728984320385571D0
    TABhighestwtf(4,3) = 0.43D0
    omegaTAB(4,3) = 0.143620623302905D0 
    TABKsp(4,3) = 24.0D0 ![molal basis]

    !NH4+ <-> I-            
    bTABAC(4,4) = 9.2191676199143202D-02
    cTABAC(4,4) = -4.6783666750062035D-02
    cn1TABAC(4,4) = 0.0D0
    cn2TABAC(4,4) =  8.1996177814860266D-02
    omega2TAB(4,4) =  1.3240421566062248D+00
    TABKsp(4,4) = 74.07D0  ![molal basis]

    !NH4+ <-> NO3-    5P  JUN 2007  
    bTABAC(4,5) = -5.735407413065387D-005
    cTABAC(4,5) = -0.171746228783622D0  
    cn1TABAC(4,5) = 0.0D0       
    cn2TABAC(4,5) = 5.510450416893165D-003   
    omega2TAB(4,5) = 0.529761886212099D0     
    TABhighestwtf(4,5) = 0.90D0 
    omegaTAB(4,5) = 0.26D0  
    TABKsp(4,5) = 12.0D0 ![molal basis]
    
    !NH4+ <-> IO3-   estimated by NH4+ <-> NO3- 
    bTABAC(4,6) = bTABAC(4,5) 
    cTABAC(4,6) = cTABAC(4,5)            
    cn1TABAC(4,6) = cn1TABAC(4,5)      
    cn2TABAC(4,6) = cn2TABAC(4,5)
    omega2TAB(4,6) = omega2TAB(4,5)
    omegaTAB(4,6) = omegaTAB(4,5)

    !NH4+ <-> OH-   estimated by K+ <-> OH-   
    bTABAC(4,7) = bTABAC(3,7)
    cTABAC(4,7) =  cTABAC(3,7)
    cn1TABAC(4,7) = cn1TABAC(3,7)
    cn2TABAC(4,7) = cn2TABAC(3,7)
    omega2TAB(4,7) = omega2TAB(3,7)

    !NH4+ <-> HSO4-    5P FEB 2011
    bTABAC(4,8) = 7.59734841D-03
    cTABAC(4,8) = 1.43011653D-01
    omegaTAB(4,8) = 2.03954033D-01
    cn1TABAC(4,8) = 0.0D0 
    cn2TABAC(4,8) = 6.31184288D-03
    omega2TAB(4,8) = 8.25385786D-01
    TABhighestwtf(4,8) = 0.90D0 

    !NH4+ <-> HCO3-    Jan 2021
    bTABAC(4,10) = 1.7140258422482357D-01
    cTABAC(4,10) = -2.1046493522524873D-01
    omegaTAB(4,10) = 4.3530915388705860D-01
    cn1TABAC(4,10) = 0.0D0 
    cn2TABAC(4,10) = -4.5619906309823675D-02
    omega2TAB(4,10) =  5.8023505996780234D-01

    !(NH4)2SO4  [resp. NH4+ <-> SO4--]  5P JUN 2007 
    bTABAC(4,21) = 3.729564498700757D-004
    cTABAC(4,21) = -0.906075204598434D0
    omegaTAB(4,21) = 0.54510882832D0
    cn1TABAC(4,21) = 0.0D0
    cn2TABAC(4,21) = -3.790519418639869D-004  
    omega2TAB(4,21) = 0.354206472573787D0
    TABhighestwtf(4,21) = 0.78D0
    TABKsp(4,21) = 1.03D0 ![molal basis]  

    !NH4+ <-> CO3--     Sept 2020
    bTABAC(4,22) = 7.8317937289243134E-03
    cTABAC(4,22) = -2.7115951835499486D0
    cn1TABAC(4,22) = 0.0D0 
    cn2TABAC(4,22) = -7.2643506216563458D-02
    omega2TAB(4,22) = 5.3625142404031889D-01

    !H+ <-> Cl-      4P JUN 2007
    bTABAC(5,2) = 0.182003163561796D0
    cTABAC(5,2) = 0.243340062772300D0  
    cn1TABAC(5,2) = 0.0D0
    cn2TABAC(5,2) = 3.331893957013864D-002 
    omega2TAB(5,2) = 0.504672383721131D0  
    TABhighestwtf(5,2) = 0.18D0

    !H+ <-> Br-      4P JUN 2007
    bTABAC(5,3) = 0.120325353682262D0 
    cTABAC(5,3) = 0.444858839503779D0 
    cn1TABAC(5,3) = 0.0D0
    cn2TABAC(5,3) = 8.076729905421279D-002
    omega2TAB(5,3) = 0.596775891883278D0
    TABhighestwtf(5,3) = 0.20D0  

    !H+ <-> I-  
    bTABAC(5,4) = 5.091034D-01 
    cTABAC(5,4) = 1.793665D-01
    cn1TABAC(5,4) = 0.0D0
    cn2TABAC(5,4) = -9.943891D-02
    omega2TAB(5,4) = 1.371292D0  

    !H+ <-> NO3-     4P JUN 2007
    bTABAC(5,5) = 0.210637923183262D0 
    cTABAC(5,5) = 0.122694135055374D0
    cn1TABAC(5,5) = 0.0D0
    cn2TABAC(5,5) = -0.101735827610270D0
    omega2TAB(5,5) = 1.67641985798924D0
    TABhighestwtf(5,5) = 0.16D0

    !H+ <-> IO3-     July 2021
    bTABAC(5,6) = 1.4981114516279127D-02
    cTABAC(5,6) = -3.2063552123769585D+00
    cn1TABAC(5,6) = 0.0D0
    cn2TABAC(5,6) = -4.2730539646547307D-02
    omega2TAB(5,6) = 5.3596305879882622D-01

    !H+ <-> OH- 
    bTABAC(5,7) = 0.0D0
    cTABAC(5,7) =  0.0D0
    cn1TABAC(5,7) = 0.0D0
    cn2TABAC(5,7) = 0.0D0
    omega2TAB(5,7) = 0.0D0
    omegaTAB(5,7) = 0.0D0

    !H2SO4 (-> HSO4- & H+) 5P FEB 2011
    bTABAC(5,8) = 2.15532299D-02  
    cTABAC(5,8) = 5.62965674D-01      
    omegaTAB(5,8) = 1.42442019D-01 !!   
    cn1TABAC(5,8) = 0.0D0
    cn2TABAC(5,8) = 7.03842440D-02
    omega2TAB(5,8) = 7.14194282D-01      
    TABhighestwtf(5,8) = 0.82D0

    !H2CO3 (-> H+ & HCO3- )  Jan 2021
    bTABAC(5,10) = 4.8650377643706250D-02
    cTABAC(5,10) = -4.7365194558795232D-01   
    cn1TABAC(5,10) = 0.0D0
    cn2TABAC(5,10) = 9.9088551465819994D-03
    omega2TAB(5,10) = 3.2804860454964357D-01  

    !HSO4-  (-> H+ & SO4--) 5P FEB 2011
    bTABAC(5,21) = 2.863428226D-01
    cTABAC(5,21) = -5.996148056D0   
    omegaTAB(5,21) = 1.368612639D0 
    cn1TABAC(5,21) = 0.0D0
    cn2TABAC(5,21) = -5.359772603D-01                 
    omega2TAB(5,21) = 9.071999765D-01 
    TABhighestwtf(5,21) = 0.82D0

    !HCO3-  (-> H+ & CO3--)  Jan 2021
    bTABAC(5,22) = -1.5220470327665464D-01 
    cTABAC(5,22) = 1.7020149536615228D-01  
    cn1TABAC(5,22) = 0.0D0
    cn2TABAC(5,22) = 8.6856221337292705D-02
    omega2TAB(5,22) = 5.9622219134548793D-01

    !Ca++ <-> Cl-         4P JUN 2007
    bTABAC(21,2) = 0.104920450619807D0     
    cTABAC(21,2) = 0.866923035893255D0    
    cn1TABAC(21,2) = 0.0D0
    cn2TABAC(21,2) = 7.206272036095257D-002    
    omega2TAB(21,2) = 0.365747420374652D0   
    TABhighestwtf(21,2) = 0.53D0

    !Ca++ <-> Br-         4P JUL 2010
    bTABAC(21,3) = 0.890929D0
    cTABAC(21,3) = 6.101342D-002
    cn1TABAC(21,3) = 0.0D0
    cn2TABAC(21,3) = -0.238788D0   
    omega2TAB(21,3) = 0.762961D0
    TABhighestwtf(21,3) = 0.54D0

    !Ca++ <-> I- 
    bTABAC(21,4) = 7.5976411109783837D-01
    cTABAC(21,4) = 1.1031941621629198D+00  
    cn1TABAC(21,4) = 0.0D0
    cn2TABAC(21,4) = -6.7243764463453726D-01
    omega2TAB(21,4) = 1.5180445738268933D+00

    !Ca++ <-> NO3-        4P MAR 2006
    bTABAC(21,5) = 0.163281976992298D0 
    cTABAC(21,5) = 0.203681108454362D0  
    cn1TABAC(21,5) = 0.0D0
    cn2TABAC(21,5) = -7.545167774717340D-002
    omega2TAB(21,5) = 1.21090585828713D0
    TABhighestwtf(21,5) = 0.50D0 
    TABKsp(21,5) = 2100.0D0 ![molal basis]  

    !Ca++ <-> IO3-        estimated by Ca++ <-> NO3-
    bTABAC(21,6) = bTABAC(21,5)
    cTABAC(21,6) = cTABAC(21,5) 
    cn1TABAC(21,6) = cn1TABAC(21,5)
    cn2TABAC(21,6) = cn2TABAC(21,5)
    omega2TAB(21,6) = omega2TAB(21,5)

    !Ca++ <-> OH-      estimated by Ca++ <-> I-
    bTABAC(21,7) = bTABAC(21,4) 
    cTABAC(21,7) = cTABAC(21,4)
    cn1TABAC(21,7) = cn1TABAC(21,4)
    cn2TABAC(21,7) = cn2TABAC(21,4)
    omega2TAB(21,7) =  omega2TAB(21,4)

    !Ca++ <-> HSO4-     AUG 2019
    bTABAC(21,8) = 8.2098127298847412D-02
    cTABAC(21,8) = -1.7166655056598561D+00
    cn1TABAC(21,8) = 0.0D0
    cn2TABAC(21,8) = 8.6726419889762707D-02
    omega2TAB(21,8) = 4.6756611045939778D-01

    !Ca++ <-> HCO3- !Jan 2021 
    bTABAC(21,10) = 4.0120219199493851D-01
    cTABAC(21,10) = -7.9562337248535087D-01
    cn1TABAC(21,10) = 0.0D0
    cn2TABAC(21,10) = -6.4427030445717304D-01
    omega2TAB(21,10) = 1.0886391158515725D+00

    !Ca++ <-> SO4--      4P JAN 2011
    bTABAC(21,21) = 1.2956723534082162D0
    cTABAC(21,21) = -6.9680634914651929D-01
    omegaTAB(21,21) = 0.8D0 !usually set to 0.8
    cn1TABAC(21,21) = 0.0D0
    cn2TABAC(21,21) = 1.5915889122259039D0
    omega2TAB(21,21) = 2.5621741189504110D-01
    TABhighestwtf(21,21) = 2.0D-3 !CaSO4   

    !Ca++ <-> CO3--  !Jan 2021 
    bTABAC(21,22) = 6.0377736237857693D-01
    cTABAC(21,22) = -8.0427942772938049D-03
    omegaTAB(21,22) = 1.0257268968493378D+00
    cn1TABAC(21,22) = 0.0D0
    cn2TABAC(21,22) = -4.9989429639370098D-01
    omega2TAB(21,22) = 1.5894433704988793D+00

    !Mg++ <-> Cl-        4P JUN 2007
    bTABAC(23,2) = 0.195908827051115D0
    cTABAC(23,2) = 0.332387167778681D0
    cn1TABAC(23,2) = 0.0D0
    cn2TABAC(23,2) = 7.865437200679259D-002    
    omega2TAB(23,2) = 0.397920452995359D0
    TABhighestwtf(23,2) = 0.32D0    !MgCl2
    TABKsp(23,2) = 1.8986D4  !hydrate-Ksp [molal basis] for saturated solution with respect to MgCl2.6H2O (hydrate); solubility value from Rard: 5.8101 mol MgCl2/kg water

    !Mg++ <-> Br-        4P JUL 2010
    bTABAC(23,3) = 0.260487D0
    cTABAC(23,3) = 1.017037D0
    cn1TABAC(23,3) = 0.0D0
    cn2TABAC(23,3) = 6.162636D-002    
    omega2TAB(23,3) = 0.299475D0
    TABhighestwtf(23,3) = 0.48D0    !MgBr2

    !Mg++ <-> I- 
    bTABAC(23,4) = 9.6129689010970776D-01
    cTABAC(23,4) = 7.7356425470356993D-01   
    cn1TABAC(23,4) = 0.0D0
    cn2TABAC(23,4) = -4.7285705361984376D-01
    omega2TAB(23,4) =  1.0827821509817128D+00

    !Mg++ <-> NO3-      4P JUN 2007
    bTABAC(23,5) = 0.430670745401918D0
    cTABAC(23,5) = 0.767241847531231D0  
    cn1TABAC(23,5) = 0.0D0
    cn2TABAC(23,5) = -0.511836366410133D0 
    omega2TAB(23,5) = 1.44093984791712D0 
    TABhighestwtf(23,5) = 0.54D0
    TABKsp(23,5) = 30383.0D0 ![molal basis]     

    !Mg++ <-> IO3-     JUN 2021
    bTABAC(23,6) = 1.2126052701130180D-01
    cTABAC(23,6) = 1.7147232434049486D-01
    cn1TABAC(23,6) = 0.0D0
    cn2TABAC(23,6) = -5.3215011337300289D-01
    omega2TAB(23,6) = 8.7729205968138912D-01

    !Mg++ <-> OH-      estimated by Mg++ <-> I-
    bTABAC(23,7) = bTABAC(23,4) 
    cTABAC(23,7) = cTABAC(23,4)
    cn1TABAC(23,7) = cn1TABAC(23,4)
    cn2TABAC(23,7) = cn2TABAC(23,4)
    omega2TAB(23,7) =  omega2TAB(23,4)

    !Mg++ <-> HSO4-     JULY 2019
    bTABAC(23,8) = 1.7141381771763914D-01
    cTABAC(23,8) =  2.7138143502966559D+00
    cn1TABAC(23,8) = 0.0D0
    cn2TABAC(23,8) = 2.0404160088725892D-01
    omega2TAB(23,8) = 5.0134675159022479D-01

    !Mg++ <-> HCO3-  !Jan 2021 
    bTABAC(23,10) = 1.8773950264361727D-01
    cTABAC(23,10) = -1.0243549517195567D+00
    omegaTAB(23,10) = 1.8303653537783598D-01
    cn1TABAC(23,10) = 0.0D0
    cn2TABAC(23,10) = 2.1693552210594391D-01
    omega2TAB(23,10) = 1.3960864277476452D+00

    !Mg++ <-> SO4--      4P JUN 2007
    bTABAC(23,21) = 0.122364180080768D0
    cTABAC(23,21) =  -3.42587571353539D0
    cn1TABAC(23,21) = 0.0D0
    cn2TABAC(23,21) = -0.738560814131576D0   
    omega2TAB(23,21) = 0.864380270721573D0
    TABhighestwtf(23,21) = 0.60D0 !MgSO4 
    TABKsp(23,21) = 999999.0D0 ![molal basis]   !this value is not correct, just as initial value here. MgSO4 forms a heptahydrate at room temp.: MgSO4.7H2O

    !Mg++ <-> CO3--  !Jan 2021 
    bTABAC(23,22) = 6.8615185959383962D-01
    cTABAC(23,22) = -2.2381706037879368D+00
    omegaTAB(23,22) = 1.0578308708725501D-01
    cn1TABAC(23,22) = 0.0D0
    cn2TABAC(23,22) = -5.6539053148654073D-01
    omega2TAB(23,22) = 1.2004054440397884D+00

    !H+ <-> CH3SO3-  5P JUN 2017;  [MSA = methanesulfonic acid]
    bTABAC(5,9) = 4.679002D-2
    cTABAC(5,9) = 4.562208D-1
    omegaTAB(5,9) = 0.24D0
    cn1TABAC(5,9) = 0.0D0
    cn2TABAC(5,9) = -1.702953D-1 
    omega2TAB(5,9) = 1.428179D0
    TABhighestwtf(5,9) = 0.791D0 !MSA

    !Na+ <-> CH3SO3-  5P JUN 2017;  [CH3SO3- is the methanesulfonate anion]
    bTABAC(2,9) = 6.9073215D-03
    cTABAC(2,9) =  3.166402726D-01
    omegaTAB(2,9) = 3.2567310229D-01
    cn1TABAC(2,9) = 0.0D0
    cn2TABAC(2,9) = -4.3807515687D-01
    omega2TAB(2,9) = 2.109695752D0
    TABhighestwtf(2,9) = 0.64D0

    !NH4+ <-> CH3SO3-  5P JUN 2017;  [CH3SO3- is the methanesulfonate anion]
    bTABAC(4,9) = 2.337520821D-02
    cTABAC(4,9) = 2.23469882D-02
    omegaTAB(4,9) = 2.20000D-01
    cn1TABAC(4,9) = 0.0D0
    cn2TABAC(4,9) = 2.472697276D-02
    omega2TAB(4,9) = 8.265735355D-01
    TABhighestwtf(4,9) = 0.64D0

    !---------------------------------------------------------------- 
    !The interactions of three different ions are described with Qcca and Qcaa:
    !There is also an interaction parameter between NH4+ and H+ ions (Rcc) declared here:
    !Allocate RccTAB with limited array size because these values are only used for a few selected ion interactions;
    ALLOCATE( qcca1TAB(23,23,40), RccTAB(23,23) )  !qcca1TAB(40,40,40)
    qcca1TAB = 0.0D0
    RccTAB = 0.0D0

    !Rcc interaction parameter between two cations: 
    !NH4+|H+ (= H+|NH4+)
    RccTAB(4,5) = -1.544864142D-01
    RccTAB(5,4) = -1.544864142D-01

    !Qcca parameters for some three body compositions: 
    !NH4+|H+|HSO4- (= H+|NH4+|HSO4-) 
    qcca1TAB(4,5,8) = 4.483540853D-04
    qcca1TAB(5,4,8) = 4.483540853D-04

    !""""""""""""""""""""" END OF  cation <-->  anion (in water) """"""""""""""""""""""""

    END SUBROUTINE MRdata
    !==========================================================================================================================
   
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to calculate the CO2 activity coefficients in the mixture.              *
    !*   This aqueous CO2 activity coefficient is defined on molality basis with reference  *
    !*   state of infinite dilution in pure water.                                          *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Hang Yin & Andi Zuend                                                              *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2020                                                            *
    !*   -> latest changes: 2020-05-26                                                      *
    !*                                                                                      *
    !****************************************************************************************
    SUBROUTINE GammaCO2()

    USE ModAIOMFACvar, ONLY : SMA, SMC, gnmrln, gnsrln, gnlrln
    USE ModSystemProp, ONLY : NGI, Ication, Ianion, idCO2

    IMPLICIT NONE

    INTEGER(4) :: I, NC, NA
    REAL(8) :: lngammaCO2
    !...............................
    
    lngammaCO2 = 0.0D0
    DO I = 1,NGI
        NC = Ication(I)
        IF (SMC(I) > 0.0D0) THEN
            lngammaCO2 = lngammaCO2 + 2.0D0*SMC(I)*lambdaIN(NC)
        ENDIF
        NA = Ianion(I)
        IF (SMA(I) > 0.0D0) THEN
            lngammaCO2 = lngammaCO2 + 2.0D0*SMA(I)*lambdaIN(NA)
        ENDIF
    ENDDO

    gnmrln(idCO2) = lngammaCO2
    gnsrln(idCO2) = 0.0D0
    gnlrln(idCO2) = 0.0D0

    END SUBROUTINE GammaCO2
    !=========================================================================================================================

END MODULE ModMRpart