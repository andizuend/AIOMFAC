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
!*   -> latest changes: 2024-04-10                                                      *
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
!*   -  subroutine LR_MR_activity                                                       *
!*   -  subroutine MRinteractcoeff                                                      *
!*   -  subroutine MRdata                                                               *
!*   -  subroutine GammaCO2                                                             *
!*                                                                                      *
!****************************************************************************************
    
module ModMRpart

use Mod_kind_param, only : wp
use ModSystemProp, only : NGI, NG, topsubno, frominpfile, Nmaingroups, bicarbsyst

implicit none
!....................................................................................
!public variables:
real(wp),dimension(:,:),allocatable,public :: Ksp
!** private module variables, which need to be set to public during parameter fitting...
real(wp),dimension(:,:),allocatable,public :: bTABna, bTABnc, cTABna, cTABnc
real(wp),dimension(:,:),allocatable,public :: Cn1TABAC, Cn2TABAC, omega2TAB, bTABAC, cTABAC
real(wp),dimension(:,:),allocatable,public :: TABhighestWTF, TABKsp, omegaTAB
!private module variables:
integer,private :: jj
real(wp),private :: A_DebyeHw, B_DebyeHw  !the temperature-dependent Debye-Hueckel parameters for water as solvent
real(wp),dimension(:,:),allocatable,private :: BAC, CAC, omega, omega2, Cnac1, Cnac2, Raa, Rcc
real(wp),dimension(:,:),allocatable,private :: bnc, cnc, bna, cna, omegaNC, omegaNA
real(wp),dimension(:,:),allocatable,private :: RccTAB
real(wp),dimension(:,:,:),allocatable,private :: Qcca
real(wp),dimension(:,:,:),allocatable,private :: qcca1TAB  !real(wp),dimension(40,40,40),private :: qcca1TAB
logical,private :: QccaInteract, RccInteract
!set CO2(aq) <--> ion parameters:
real(wp),dimension(201:topsubno),parameter,private :: lambdaIN = real([ &
  !Li+,       Na+,        K+,     NH4+,       H+, ...
& 0.0000D0, 0.1000D0, 0.0510D0, 0.01000D0, 0.0000D0, 0.000D0, 0.0000D0, 0.0000D0, 0.0000D0, 0.0000D0, & !(201:210)
& 0.0000D0, 0.0000D0, 0.0000D0, 0.00000D0, 0.0000D0, 0.000D0, 0.0000D0, 0.0000D0, 0.0000D0, 0.0000D0, & !(211:220)
& 0.1830D0, 0.0000D0, 0.1830D0, 0.00000D0, 0.0000D0, 0.000D0, 0.0000D0, 0.0000D0, 0.0000D0, 0.0000D0, & !(221:230)
& 0.0000D0, 0.0000D0, 0.0000D0, 0.00000D0, 0.0000D0, 0.000D0, 0.0000D0, 0.0000D0, 0.0000D0, 0.0000D0, & !(231:240)
& 0.0000D0, -0.005D0, 0.0000D0, 0.0000D0, -0.0457D0, 0.000D0, 0.0000D0, -0.003D0, 0.0000D0, 0.0000D0, & !(241:250)
& 0.0000D0, 0.0000D0, 0.0000D0, 0.00000D0, 0.0000D0, 0.000D0, 0.0000D0, 0.0000D0, 0.0000D0, 0.0000D0, & !(251:260)
& 0.0970D0, (0.0D0, jj = 262,topsubno) ], kind=wp)  
!....................................................................................

!$OMP THREADPRIVATE(A_DebyeHw, B_DebyeHw, BAC, CAC, omega, omega2, CnAC1, Cnac2, Raa, Rcc, &
    !$OMP & Ksp, bnc, cnc, bna, cna, omegaNC, omegaNA, Qcca, QccaInteract, RccInteract)
    
!==========================================================================================================================
    contains
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
    !*   -> latest changes: 2022-01-17                                                      *
    !*                                                                                      *
    !****************************************************************************************
    subroutine LR_MR_activity()

    use ModSystemProp, only : anionZ, cationZ, Imaingroup, ITAB, ITABMG, maingrindexofsubgr, Mmass, &
        & nd, nelectrol, NGI, NGN, nneutral, solvmixrefnd, SolvSubs, SubGroupMW, errorflag_clist
    use ModAIOMFACvar, only : DebyeHrefresh, galrln, gamrln, gclrln, gcmrln, gnlrln, gnmrln, &
        & Ionicstrength, meanSolventMW, SumIonMolalities, SMA, SMC, solvmixcorrMRa, solvmixcorrMRc, &
        & T_K, Tmolal, TmolalSolvmix, wtf, XN

    implicit none
    !Local Variables
    integer :: I, II, J, K
    real(wp) :: A, b, bb, SI, SI2, Z, Test11, Test12, TermI, Sumkion1,  &
    & sumca1, sumq1, sumkion2, sumbsca, sumcnca, sumcsca, sumqc1, SumRc, sumxg, SumQa2, & 
    & lrw1, twoZ, ZSI, avgmaingrMW, Mavg
    real(wp),parameter :: dtiny = epsilon(1.0_wp)
    real(wp),parameter :: Mwater = 0.01801528_wp    !the molar mass of water [kg/mol]
    real(wp),parameter :: cDens = 997.0_wp          !density of water (kg/m^3) at 298.15 K
    real(wp),parameter :: densmix = cDens           !here by design the density of reference solvent water at 298.15 K is used
    real(wp),parameter :: sqrtDensmix = sqrt(densmix)
    real(wp),parameter :: cDiel = 78.54E0_wp        !static dielectric const. (= relative static permittivity) of water at 298.15 K
    real(wp),parameter :: dielmix = cDiel
    real(wp),dimension(NGI) :: SumABca, SumBcx, SumBax, SumCBca, ZA, ZA2, ZC2, smcsma, oexp1, oexp2
    real(wp),dimension(NGN) :: subgroupxsf, mgroupxsf, mgroupMK, SumBM, gkmrln
    real(wp),dimension(NGI,NGI) :: BSca, CSca, Bca, Cnca
    real(wp),dimension(NGN,NGI) :: BKionc, BKiona, BSKionc, BSKiona
    !....................................................................

    !Calculation of the salt-free mole fractions of the neutral subgroups (combined contributions from different molecules if the same subgroup is present):
    do I = 1,NGN
        II = SolvSubs(I)
        subgroupxsf(I) = sum(ITAB(1:nneutral,II)*XN(1:nneutral))
    enddo
    !normalize salt-free basis subgroup mole fractions:
    sumxg = sum(subgroupxsf)
    if (sumxg > 0.0_wp) then
        subgroupxsf = subgroupxsf/sumxg
    endif

    !calculate the (salt-free basis) mole fractions, mgroupxsf, of the organic main groups present:
    mgroupxsf = 0.0_wp
    do j = 1,NGN !j-loop to cover subgroups
        I = maingrindexofsubgr(j) !this maingrindexofsubgr array assigns the index of the main group (of the current mixture) associated with subgroup index j
        mgroupxsf(I) = mgroupxsf(I) +subgroupxsf(j)
    enddo

    !Assignment of the molecular masses to subgroups:
    !mgroupMK: lists the molecular masses of the main groups (of Imaingroup) present in the mixture (as calculated based on the current contributions from different subgroups). 
    !calculate mole-fractional contributions of different subgroups to a certain main group for calculation of main group molar mass:
    mgroupMK = 0.0_wp
    do j = 1,NGN !j-loop to cover subgroups
        I = maingrindexofsubgr(j)
        if (mgroupxsf(I) > 0.0_wp) then !check for potential divison by zero problem
            mgroupMK(I) = mgroupMK(I) + SubGroupMW(j)*subgroupxsf(j)/mgroupxsf(I) !normalized contributions of subgroups to main group molar mass [kg/mol]
        endif
    enddo
    Mavg = sum(SubGroupMW(1:NGN)*subgroupxsf(1:NGN))

    !Calculation of the density of the salt free mixture
    !**--++1 water LR 
    !densmix = cDens !sum(compVF(1:nneutral)*cDens(1:nneutral))
    !!Calculation of the dielectric constant of the salt free mixture:
    !dielmix = cDiel != sum(compVF(1:nneutral)*cDiel(1:nneutral))
    if (DebyeHrefresh) then !only calculate the first time and when T_K changes significantly
        A_DebyeHw = 1.327757E5_wp*sqrtDensmix/((dielmix*T_K)**1.5_wp)
        B_DebyeHw = 6.359696_wp*sqrtDensmix/((dielmix*T_K)**0.5_wp)
    endif
    A = A_DebyeHw
    b = B_DebyeHw

    !Calculation of the ionic strength (SI): (the Ionic strength in the molality base)
    ZA = abs(anionZ(1:NGI))
    ZA2 = ZA**2
    ZC2 = cationZ(1:NGI)**2
    SI = 0.5_wp*sum( SMA(1:NGI)*ZA2 + SMC(1:NGI)*ZC2 )
    SI2 = sqrt(SI)
    Ionicstrength = SI
    !Calculation of Z (number of charges per kg solvent):
    Z = sum( SMA(1:NGI)*ZA(1:NGI) + SMC(1:NGI)*cationZ(1:NGI) )
    !for debugging tests:
    Test11 = sum(SMA(1:NGI)*ZA(1:NGI))                      !Test of the electrical charge neutrality condition in the solution:  !test test
    Test12 = sum(SMC(1:NGI)*cationZ(1:NGI))
    Test11 = Test11 - Test12
    if (.NOT. any([bicarbsyst])) then
        if (SI < 1.0E6_wp .AND. abs(Test11) > 1.0E-9_wp) then     !test
            errorflag_clist(12) = .true.                    !overall charge neutrality violation error
            if (.NOT. frominpfile) then
                write(*,*) "WARNING: Electrical charge neutrality condition in the solution is violated! ", Test11
            endif
        endif
    endif

    !:: LR :: calculation of the long-range component of the activity coefficients of the neutrals:
    bb = 1.0_wp+b*SI2
    TermI = bb-1.0_wp/(bb)-2.0_wp*log(bb) 
    !all LR properties set to those of pure water (except for the molar mass):
    lrw1 = 2.0_wp*A/(b*b*b)*TermI
    gnlrln(1:nneutral) = Mmass(1:nneutral)*lrw1

    !:: LR :: calculation of the long-range component of the activity coefficients of the ions:
    lrw1 = A*SI2/bb
    gclrln = -ZC2*lrw1
    galrln = -ZA2*lrw1

    !################################################################################################################
    ! :: MR :: calculation of the middle-range contribution of the activity coefficients of the ions and neutrals:
    !################################################################################################################
    if (nelectrol > 0 .AND. SI > dtiny) then !**--++
        SumCBca = 0.0_wp
        SumABca = 0.0_wp
        !Calculation of Bca, Cnca, Bkionc and Bkiona:
        lrw1 = -0.5_wp/SI2
        do J = 1,NGI
            oexp1(1:NGI) = omega(1:NGI,J)*SI2 
            where (oexp1(1:NGI) > 300.0_wp) !prevent floating point underflow as the second term becomes practically zero:
                Bca(1:NGI,J) = bac(1:NGI,J) !+cac(1:NGI,J)*exp(-300.0_wp)
            elsewhere
                Bca(1:NGI,J) = bac(1:NGI,J) +cac(1:NGI,J)*exp(-oexp1(1:NGI))
            end where
            oexp2(1:NGI) = omega2(1:NGI,J)*SI2
            where (oexp2(1:NGI) > 300.0_wp) !prevent floating point underflow as the second term becomes practically zero:
                Cnca(1:NGI,J) = cnac1(1:NGI,J) !+cnac2(1:NGI,J)*exp(-300.0_wp)
            elsewhere
                Cnca(1:NGI,J) = cnac1(1:NGI,J) +cnac2(1:NGI,J)*exp(-oexp2(1:NGI))
            end where
            ! Calculation of derivatives of Bca, Cnca with respect to the ionic strength SI. BSca: derivative of Bca.
            BSca(1:NGI,J) = lrw1*omega(1:NGI,J)*(Bca(1:NGI,J)-bac(1:NGI,J)) 
            CSca(1:NGI,J) = lrw1*omega2(1:NGI,J)*(Cnca(1:NGI,J)-cnac1(1:NGI,J))
        enddo

        do J = 1,NGI
            if (SI2 > 250.0_wp) then !prevent floating point underflow as the second term becomes practically zero:
                Bkionc(1:NGN,J) = bnc(1:NGN,J)
                Bkiona(1:NGN,J) = bna(1:NGN,J)
            else
                Bkionc(1:NGN,J) = bnc(1:NGN,J) +cnc(1:NGN,J)*exp(-omegaNC(1:NGN,J)*SI2)
                Bkiona(1:NGN,J) = bna(1:NGN,J) +cna(1:NGN,J)*exp(-omegaNA(1:NGN,J)*SI2)
            endif
            !Calculation of derivatives of Bkionc and Bkiona with respect to the ionic strength SI: BSkionc: derivative of Bkionc; BSkiona: derivative of Bkiona
            BSkionc(1:NGN,J) = lrw1*omegaNC(1:NGN,J)*(Bkionc(1:NGN,J)-bnc(1:NGN,J))  != (-0.5_wp*omegaNC(I,J)/SI2)*cnc(I,J)*exp(-omegaNC(I,J)*SI2)
            BSkiona(1:NGN,J) = lrw1*omegaNA(1:NGN,J)*(Bkiona(1:NGN,J)-bna(1:NGN,J))  != (-0.5_wp*omegaNA(I,J)/SI2)*cna(I,J)*exp(-omegaNA(I,J)*SI2)
        enddo

        !the mean molecular mass of a main group in the solvent mixture.
        avgmaingrMW = sum(mgroupMK*mgroupxsf) !this is the mean main group molar mass -- and not the mean solvent molar mass!

        !Calculation of the midrange part for the neutrals:
        do I = 1,NGN
            SumBm(I) = sum(Bkiona(I,1:NGI)*SMA(1:NGI) +Bkionc(I,1:NGI)*SMC(1:NGI))
        enddo
        Sumkion1 = 0.0_wp
        do J = 1,NGI
            Sumkion1 = Sumkion1 +sum( (Bkiona(1:NGN,J)+SI*BSkiona(1:NGN,J))*mgroupxsf(1:NGN)*SMA(J) &
                                   & +(Bkionc(1:NGN,J)+SI*BSkionc(1:NGN,J))*mgroupxsf(1:NGN)*SMC(J) )
        enddo
        SumCA1 = 0.0_wp
        twoZ = 2.0_wp*Z
        ZSI = Z*SI
        do J = 1,NGI
            SumCA1 = SumCA1 +sum((Bca(1:NGI,J)+SI*BSca(1:NGI,J)+twoZ*Cnca(1:NGI,J)+ZSI*Csca(1:NGI,J))*SMC(1:NGI)*SMA(J))
        enddo

        !Calculation of the three ion interactions in the neutral part:      
        SumQ1 = 0.0_wp
        if (QccaInteract) then !so far only used for [NH4+|H+|HSO4-] containing mixtures
            do I = 1,NGI
                do J = I,NGI
                    SumQ1 = SumQ1 +sum(2.0_wp*Qcca(I,J,1:NGI)*SMC(I)*SMC(J)*SMA(1:NGI)) !term when ionic strength dependency is parameterized: SumQ1+(2.0_wp*Qcca(I,J,K)+SI*QScca(I,J,K))*SMC(I)*SMC(J)*SMA(K)
                enddo
            enddo
        endif
        ! Calculation of the two cation interactions in the neutral part:
        SumRc = 0.0_wp
        if (RccInteract) then !so far only used for [NH4+|H+|HSO4-] containing mixtures
            do I = 1,NGI
                SumRc = SumRc +sum(Rcc(I,I:NGI)*SMC(I)*SMC(I:NGI)) 
            enddo
        endif

        !gkmrln(I): ln of the middle-range part of the activity coefficient of main group I
        gkmrln(1:NGN) = SumBm(1:NGN) -mgroupMK(1:NGN)/avgmaingrMW*Sumkion1 -mgroupMK(1:NGN)*SumCA1 &
            & -mgroupMK(1:NGN)*SumQ1 -mgroupMK(1:NGN)*SumRc

        !Calculation of the middle-range part of the activity coefficients of the neutral molecules from the contributions of main groups to individual molecules
        !gnmrln(I): ln of the middle-range part of the activity coefficient of neutral molecule I
        gnmrln = 0.0_wp
        do J = 1,NGN
            II = Imaingroup(J)
            if (II > 0) then
                do I = 1,nneutral
                    if (ITABMG(I,II) > 0) then
                        do K = 1,NGN
                            if (Imaingroup(K) == II) then
                                gnmrln(I) = gnmrln(I) +ITABMG(I,II)*gkmrln(K)
                            endif
                        enddo
                    endif
                enddo
            endif
        enddo

        !Calculation of the middle-range contribution part for the cations and anions: @@##
        SumBcx = 0.0_wp
        SumBax = 0.0_wp
        do I = 1,NGI
            SumBCx(I) = sum(Bkionc(1:NGN,I)*mgroupxsf(1:NGN)) 
            SumBax(I) = sum(Bkiona(1:NGN,I)*mgroupxsf(1:NGN)) 
        enddo
        SumBsCA = 0.0_wp
        SumCnca = 0.0_wp
        SumCsca = 0.0_wp
        do J = 1,NGI
            SumCBca(1:NGI) = SumCBca(1:NGI) +(BCA(1:NGI,J) +Z*Cnca(1:NGI,J))*SMA(J)
            SumABca(1:NGI) = SumABca(1:NGI) +(BCA(J,1:NGI) +Z*Cnca(J,1:NGI))*SMC(J)
            smcsma(1:NGI) = SMC(1:NGI)*SMA(J)
            SumBsCA = sumBsCA +sum(BSca(1:NGI,J)*smcsma(1:NGI))
            SumCnca = sumCnca +sum(Cnca(1:NGI,J)*smcsma(1:NGI))
            SumCsca = sumCsca +sum(Csca(1:NGI,J)*smcsma(1:NGI))
        enddo
        Sumkion2 = 0.0_wp
        do J = 1,NGI
            Sumkion2 = Sumkion2+sum(BSkiona(1:NGN,J)*mgroupxsf(1:NGN)*SMA(J)+BSkionc(1:NGN,J)*mgroupxsf(1:NGN)*SMC(J))
        enddo

        !gcmrln(I): ln of the middle-range part of the activity coefficient of cation I
        !gamrln(I): ln of the middle-range part of the activity coefficient of anion I
        do I = 1,NGI
            !cations:
            gcmrln(I) = sumbcx(I)/avgmaingrMW +ZC2(I)/(2.0_wp*avgmaingrMW)*sumkion2 +sumCBca(I) &
                & +0.5_wp*ZC2(I)*SumBsCA +cationZ(I)*SumCnca +0.5_wp*ZC2(I)*Z*SumCsca
            if (QccaInteract) then !so far only used for [NH4+|H+|HSO4-] containing mixtures
                SumQc1 = 0.0_wp
                do K = 1,NGI
                    SumQc1 = SumQc1 +sum(SMC(1:NGI)*SMA(K)*Qcca(I,1:NGI,K))  !entspr. Qc1,c,a(I)
                enddo
                gcmrln(I) = gcmrln(I) +SumQc1 !the 3 ion interaction terms... 
            endif
            !there might be an Interaction between two cations which is important for the system-fit (NH4+ <-> H+):
            if (RccInteract) then !so far only used for [NH4+|H+] containing mixtures
                gcmrln(I) = gcmrln(I) +sum(Rcc(I,1:NGI)*SMC(1:NGI))
            endif 
            !-- anions:
            gamrln(I) = SumBax(I)/avgmaingrMW +ZA2(I)/(2.0_wp*avgmaingrMW)*sumkion2 +sumABca(I) &
                & +0.5_wp*ZA2(I)*SumBsCA+ZA(I)*SumCnca +0.5_wp*ZA2(I)*Z*SumCsca
            if (QccaInteract) then !so far only used for [NH4+|H+|HSO4-] containing mixtures
                SumQa2 = 0.0_wp
                do J = 1,NGI
                    SumQa2 = SumQa2 +sum(SMC(J)*SMC(J:NGI)*Qcca(J,J:NGI,I)) 
                enddo
                gamrln(I) = gamrln(I)+SumQa2
            endif
        enddo

        !Calculation of the term that is needed for the ions to convert from mole fraction to molality units.
        !The reference solvent is water for MR and SR and solvent mixture for LR (but LR is anyways set to water properties):
        Tmolal = log(Mwater/meanSolventMW +Mwater*SumIonMolalities)  !this is the molality basis conversion term with reference solvent of pure water.
        !----
        if (solvmixrefnd) then  !calculate additional terms required to correct for different reference state in AIOMFAC_calc
            TmolalSolvmix = log(1.0_wp +SumIonMolalities*meanSolventMW)  !here reference state is the salt-free solvent mixture, because the LR-term is always with reference solvent mixture (unfortunately)!!
            solvmixcorrMRc = 0.0_wp
            solvmixcorrMRa = 0.0_wp
            do I = 1,NGI
                solvmixcorrMRc(I) = solvmixcorrMRc(I)+sum( (bnc(1:NGN,I) +cnc(1:NGN,I))*mgroupxsf(1:NGN) )
                solvmixcorrMRa(I) = solvmixcorrMRa(I)+sum( (bna(1:NGN,I) +cna(1:NGN,I))*mgroupxsf(1:NGN) )
            enddo
            solvmixcorrMRc(1:NGI) = -solvmixcorrMRc(1:NGI)/meanSolventMW
            solvmixcorrMRa(1:NGI) = -solvmixcorrMRa(1:NGI)/meanSolventMW
        endif
        !----
    else  !**--++ (no ions present in mixture)
        gamrln = 0.0_wp
        gcmrln = 0.0_wp
        gnmrln = 0.0_wp
        Tmolal = 0.0_wp
        if (solvmixrefnd) then
            TmolalSolvmix = 0.0_wp
            solvmixcorrMRc = 0.0_wp
            solvmixcorrMRa = 0.0_wp
        endif
    endif !**--++

    end subroutine LR_MR_activity
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
    !*   -> latest changes: 2018-05-18                                                      *
    !*                                                                                      *
    !****************************************************************************************
    subroutine MRinteractcoeff()
   
    use ModSystemProp, only : NGN, NGI, NG, Ication, Ianion, Imaingroup, Ncation, isPEGsystem, &
        & CatNr, AnNr, errorflagmix, frominpfile 
    
    implicit none
    !local variables:
    integer :: I, J, NN, NC, NA, NC2, K
    !.................................................................................
    
    !allocate public arrays for use in MR calculation part:
    if (allocated(Ksp)) then
        deallocate(Ksp, BAC, CAC, omega, omega2, Cnac1, Cnac2, bnc, cnc, bna, cna, omegaNC, omegaNA)
    endif
    allocate( Ksp(NGI,NGI), BAC(NGI,NGI), CAC(NGI,NGI), omega(NGI,NGI), omega2(NGI,NGI), Cnac1(NGI,NGI),  &
        & Cnac2(NGI,NGI), bnc(NG,NGI), cnc(NG,NGI), bna(NG,NGI), cna(NG,NGI), omegaNC(NG,NGI), omegaNA(NG,NGI) )
    
    !----------------------------------------------------------------   
    ! ASsignMENT OF THE MR-INTERACTION COEFFICIENTS
    !---------------------------------------------------------------- 
    !initialize neutral <-> ions interaction parameters:
    BNC = 0.0_wp
    BNA = 0.0_wp
    CNC = 0.0_wp
    CNA = 0.0_wp
    omegaNC = 1.2_wp
    omegaNA = 1.2_wp
      
    !calculate the neutral main group <-> ion interaction coefficients:        
    do I = 1,NGN
        NN = Imaingroup(I)
        if (NN > 0) then
            do J = 1,NGI
                if (Ication(J) > 0) then
                    NC = Ication(J)-200
                    BNC(I,J) = bTABnc(NN,NC)
                    CNC(I,J) = cTABnc(NN,NC)                    !CNC(I,J) is CNC(main group I, cation J), thus interaction coeff. C: main group I <-> cation J
                    if (isPEGsystem) then                       !use special parameters for PEG systems
                        if ((NN == 52 .OR. NN == 68) .AND. NC == 4) then    !CHn[OH,PEG] <--> NH4+
                            BNC(I,J) = 3.0576783E-01_wp         !fitparam(273) !bTABnc(NN,NC)
                            CNC(I,J) = 4.0411133E-01_wp         !fitparam(274) !cTABnc(NN,NC) 
                        endif
                    endif
                    !check to avoid using inteaction parameters that are not assigned correctly or were not estimated yet:
                    if (BNC(I,J) < -1000.0_wp .OR. CNC(I,J) < -1000.0_wp) then    
                        errorflagmix = 1                        !this parameter is not yet estimated nor is a fitparameter set
                        if (.NOT. frominpfile) then
                            !$OMP CRITICAL (printingtoscreen2)
                            write(*,*) ""
                            write(*,*) "======================================================="
                            write(*,*) "        WARNING from MRinteractcoeff:"
                            write(*,*) "This neutral main group <-> cation interaction coeff. is "
                            write(*,*) "not yet defined or no fit parameter was assigned. "
                            write(*,*) "NN, NC, bTABnc(NN,NC), cTABnc(NN,NC): ", NN, NC, bTABnc(NN,NC), cTABnc(NN,NC)
                            write(*,*) "======================================================="
                            write(*,*) ""
                            read(*,*)
                            !$OMP end CRITICAL (printingtoscreen2)
                        endif
                    endif
                end if
                if (Ianion(J) > 0) then
                    NA = Ianion(J)-240
                    BNA(I,J) = bTABna(NN,NA)
                    CNA(I,J) = cTABna(NN,NA)
                    if (isPEGsystem) then                                   !use special parameters for PEG systems
                        if ((NN == 52 .OR. NN == 68) .AND. NA == 21) then   !CHn[OH,PEG] <--> SO4--
                            BNA(I,J) = -2.1822703E-01_wp              !fitparam(275) !bTABna(NN,NA)
                            CNA(I,J) = -1.2595744E-01_wp              !fitparam(276) !cTABna(NN,NA) 
                        endif
                    endif
                    !check to avoid using inteaction parameters that are not assigned correctly or were not estimated yet:
                    if (BNA(I,J) < -1000.0_wp .OR. CNA(I,J) < -1000.0_wp) then !this parameter is not yet estimated nor is a fitparameter set
                        errorflagmix = 2
                        if (.NOT. frominpfile) then
                            !$OMP CRITICAL (printingtoscreen3)
                            write(*,*) ""
                            write(*,*) "======================================================="
                            write(*,*) "        WARNING from MRinteractcoeff:"
                            write(*,*) "This neutral main group <-> anion interaction coeff. is "
                            write(*,*) "not yet defined or no fit parameter was assigned. "
                            write(*,*) "NN, NA, bTABna(NN,NA), cTABna(NN,NA): ", NN, NA, bTABna(NN,NA), cTABna(NN,NA)
                            write(*,*) "======================================================="
                            write(*,*) ""
                            read(*,*)
                            !$OMP end CRITICAL (printingtoscreen3)
                        endif
                    endif
                end if
            enddo
        end if
    enddo
    !---------------------------------------------------------------- 

    !----------------------------------------------------------------  
    !initialize cation <-> anion interaction coefficients and input/output variables:
    BAC = 0.0_wp
    CAC = 0.0_wp
    CnAC1 = 0.0_wp
    CnAC2 = 0.0_wp
    Ksp = 0.0_wp     !solubility coeff. @ 298 K
    omega = 0.8_wp   !the same for all, will be overwritten in some specific cases
    omega2 = 0.6_wp

    !----------------------------------------------------------------       
    !calculate the cation <-> anion interaction coefficients: 
    do I = 1,NGI
        NC = Ication(I)-200
        do J = 1,NGI
            NA = Ianion(J)-240
            if (NC >= 1 .AND. NA >= 1) then
                BAC(I,J) = bTABAC(NC,NA)
                CAC(I,J) = cTABAC(NC,NA)
                CnAC1(I,J) = cn1TABAC(NC,NA)
                CnAC2(I,J) = cn2TABAC(NC,NA)
                omega2(I,J) = omega2TAB(NC,NA)
                omega(I,J) = omegaTAB(NC,NA)
                Ksp(I,J) = TABKsp(NC,NA)
                if (abs(bTABAC(NC,NA)) < 1.0E-10_wp) then
                    if (.NOT. (NC == 5 .AND. NA == 7)) then  !exclude H+ <--> OH- case
                        errorflagmix = 9
                        if (.NOT. frominpfile) then
                            !$OMP CRITICAL (MRI1)
                            write(*,*) "WARNING from MRinteractcoeff: interaction parameters of a cation--anion pair"
                            write(*,*) "present in the mixture is not defined: cation-ID, anion-ID: ", NC+200, NA+240
                            read(*,*)
                            !$OMP end CRITICAL (MRI1)
                        endif
                    endif
                endif
            endif
        enddo
    enddo

    !---------------------------------------------------------------- 
    !The interactions of three different ions are described with Qcca and Qcaa:
    !There is also an interaction parameter between NH4+ and H+ ions declared here:
    if (allocated(Qcca)) then
        deallocate(Qcca, Rcc, Raa)
    endif
    allocate( Qcca(NGI,NGI,NGI), Rcc(NGI,NGI), Raa(NGI,NGI) )
    Qcca = 0.0_wp
    Rcc = 0.0_wp
    Raa = 0.0_wp

    !---------------------------------------------------------------- 
    !Build up the Qcca and Rcc1 for the actual composition:
    do I = 1,NGI
        NC = Ication(I)-200
        do J = 1,NGI
            NC2 = Ication(J)-200
            if (NC >= 1 .AND. NC2 >= 1) then
                Rcc(I,J) = RccTAB(NC,NC2)
                do K = 1,NGI
                    NA = Ianion(K)-240
                    if (NA >= 1) then
                        Qcca(I,J,K) = qcca1TAB(NC,NC2,NA)
                    endif
                enddo
            endif
        enddo
    enddo
    !---------------------------------------------------------------- 

    !set logical switches for use in LR_MR_activity:
    QccaInteract = .false.
    RccInteract = .false.
    if (Ncation > 1) then  !there might be an Interaction between two cations, which is important for the system-fit (NH4+ <-> H+): 
        if (CatNr(204) > 0 .AND. CatNr(205) > 0) then !so far only used for [NH4+|H+|HSO4-] containing mixtures
            RccInteract = .true.
            if (AnNr(248) > 0) then
                QccaInteract = .true.
            endif
        endif
    endif

    end subroutine MRinteractcoeff
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
    subroutine MRdata()

    implicit none
    !................................
    
    allocate( bTABnc(Nmaingroups, 23), cTABnc(Nmaingroups, 23), &
        & bTABna(Nmaingroups, topsubno-240), cTABna(Nmaingroups, topsubno-240) )
    
    !.... initialize main group <--> ion interaction parameter data ...................
    !main group <--> cations
    bTABnc = -7.777778E+04_wp      !default value indicating that the parameter has not been determined
    cTABnc = -7.777778E+04_wp
    !load data for main group <--> cation interactions (only for the existing cation IEs, e.g. cation 01 is Li+, 02 is Na+):
    !               main group_no          1            2                3              4               5              6              7               8             9              10              11            12              13           14              15          16               17              18              19              20           21              22             23           24              25            26            27             28             29             30             31             32               33            34             35             36             37           38            39              40             41             42             43             44             45             46             47             48             49             50             51             52             53             54             55             56             57             58             59             60             61             62             63             64             65             66             67             68              69             70             71             72             73             74             75                76
    !               main_groups        (CHn),            (C=C),         (ACHn),    (ACCHn),   (OH[stand_UNIFAC]), (CH3OH[methanol]), (H2O),        (ACOH),       (CHnCO),  (CHO[aldehyde]),   (CCOO),  (HCOO[formate]),  (CHnO[ether]),   (CNH2),        (CNH),        ((C)3N),         (ACNH2),     (PYRIEINE),   (CCN[nitrile]), (COOH[stand_UNIFAC]), (CCl),     (CCl2),          (CCl3),    (CCl4),         (ACCl),         (CNO2),        (ACNO2),       (CS2),      (CH3SH),      (FURFURAL),         (EOH),          (I),           (BR),   (C#C[triple_bond]), (EMSO),       (ACRY),        (ClCC),        (ACF),        (EMF),          (CF2),          (COO),      (SIH2),         (SIO),          (NMP),         (CClF),       (CON),        (OCCOH),        (CH2S),     (MORPHOLINE),   (THIOPHENE),   (Inorg_Ions),   (CH[OH]), (OH[MingRL]), (COOH[MingRL]), (CHnCO[MingRL]), (CHn[Ming_mo]), (CHnO[Ming_mo]), (OH[Ming_mo]), (CHn[Ming_hy]), (COOH[Ming_hy]),(OH[Ming_hy]),(CHn[Ming_di]),(COOH[Ming_di]), (OH[Peng]),  (COOH),      (CHn[alc]),   (CHn[alc-tail]),  (CH2[OH]),       (OH),    (CH2OCH2[PEG]),  (CHnONO2),  (CHnOOH[perox]),(C(=O)OOH[perox]),(CHnOOCHm[perox]),(C(=O)OONO2[perox]), (CO2)  
    bTABnc(1:Nmaingroups,01) = real([ 7.961623D-02,  1.287780D-01,  6.208020D-02, -7.777778D+04,  5.877733D-03, -7.777778D+04,  0.000000D+00,  4.459015D-02,  1.117971D-01,  5.189833D-02,  1.117970D-01, -7.777778D+04,  6.128712D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  4.458747D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  4.458747D-02,  4.458747D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  6.188986D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  4.458747D-02,  5.877733D-03, -7.777778D+04,  4.458747D-02,  5.877733D-03,  4.458747D-02,  7.961623D-02,  7.961623D-02,  6.188986D-02,  5.877733D-03, -7.777778D+04,  1.242037D-01,  6.716485D-02,  4.458747D-02,  1.225742D-01,  1.242037D-01,  0.000000D+00], kind=wp)
    cTABnc(1:Nmaingroups,01) = real([ 4.669577D-02,  1.362762D-01,  3.087262D-02, -7.777778D+04, -5.629984D-02, -7.777778D+04,  0.000000D+00, -6.965679D-03,  1.422573D-01, -1.682426D-03,  1.422573D-01, -7.777778D+04,  2.696984D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -6.968065D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -6.968065D-03, -6.968065D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.583107D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -6.968065D-03, -5.629984D-02, -7.777778D+04, -6.968065D-03, -5.629984D-02, -6.968065D-03,  4.669577D-02,  4.669577D-02,  2.583107D-02, -5.629984D-02, -7.777778D+04,  3.972770D-02, -2.933000D-02, -6.968065D-03,  5.393968D-02,  3.972770D-02,  0.000000D+00], kind=wp)
    bTABnc(1:Nmaingroups,02) = real([ 1.058813D-01,  2.187202D-01,  1.058813D-01, -7.777778D+04,  7.649999D-03, -7.777778D+04,  0.000000D+00,  7.501697D-02,  1.628598D-01,  1.057398D-01,  1.628598D-01, -7.777778D+04,  1.018669D-01, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.556653D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.556653D-02,  5.556653D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  1.058805D-01, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.556653D-02,  7.649999D-03, -7.777778D+04,  5.556653D-02,  7.649999D-03,  5.556653D-02,  1.058813D-01,  1.058813D-01,  1.058805D-01,  7.649999D-03, -7.777778D+04,  1.614478D-01,  1.095169D-01,  5.556653D-02,  2.037338D-01,  1.614478D-01,  0.000000D+00], kind=wp)
    cTABnc(1:Nmaingroups,02) = real([ 2.266823D-02,  5.099094D-02,  2.266823D-02, -7.777778D+04, -7.485817D-03, -7.777778D+04,  0.000000D+00,  6.854487D-03,  5.045448D-02,  2.272655D-02,  5.045446D-02, -7.777778D+04,  1.485420D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -6.696041D-04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -6.696041D-04, -6.696041D-04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  9.730931D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -6.696041D-04, -7.485817D-03, -7.777778D+04, -6.696041D-04, -7.485817D-03, -6.696041D-04,  2.266823D-02,  2.266823D-02,  9.730931D-03, -7.485817D-03, -7.777778D+04,  2.199863D-02,  7.368383D-03, -6.696041D-04,  2.970840D-02,  2.199863D-02,  0.000000D+00], kind=wp)
    bTABnc(1:Nmaingroups,03) = real([ 9.836417D-02,  1.593944D-01,  7.030872D-02, -7.777778D+04,  2.704778D-03, -7.777778D+04,  0.000000D+00,  4.748628D-02,  1.335596D-01,  4.748629D-02,  1.018410D-01, -7.777778D+04,  6.474752D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.240469D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.240469D-02,  3.240469D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  6.847925D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.240469D-02,  2.704778D-03, -7.777778D+04,  3.240469D-02,  2.704778D-03,  3.240469D-02,  9.836417D-02,  9.836417D-02,  6.847925D-02,  2.704778D-03, -7.777778D+04,  1.307689D-01,  6.745230D-02,  3.240469D-02,  1.294950D-01,  1.307689D-01,  0.000000D+00], kind=wp)
    cTABnc(1:Nmaingroups,03) = real([ 4.203279D-02,  6.730633D-02,  3.245205D-02, -7.777778D+04, -3.237396D-02, -7.777778D+04,  0.000000D+00,  1.421022D-02,  5.039454D-02,  1.421023D-02, -4.268114D-03, -7.777778D+04,  6.472387D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  1.995250D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  1.995250D-02,  1.995250D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.207157D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  1.995250D-02, -3.237396D-02, -7.777778D+04,  1.995250D-02, -3.237396D-02,  1.995250D-02,  4.203279D-02,  4.203279D-02,  3.207157D-02, -3.237396D-02, -7.777778D+04,  6.198529D-02,  3.234991D-02,  1.995250D-02,  1.294477D-01,  6.198529D-02,  0.000000D+00], kind=wp)
    bTABnc(1:Nmaingroups,04) = real([ 5.667445D-02,  9.381847D-02,  5.584332D-02, -7.777778D+04,  5.238242D-03, -7.777778D+04,  0.000000D+00,  4.060326D-02,  9.371399D-02,  5.622542D-02,  9.371399D-02, -7.777778D+04,  4.060311D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.534892D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.534892D-02, -2.534892D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.904649D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.534892D-02,  5.238242D-03, -7.777778D+04, -2.534892D-02,  5.238242D-03, -2.534892D-02,  5.667445D-02,  5.667445D-02,  3.904649D-02,  5.238242D-03,  2.249583D-01,  3.132553D-02,  4.584135D-02, -2.534892D-02,  8.120622D-02,  3.132553D-02,  0.000000D+00], kind=wp)
    cTABnc(1:Nmaingroups,04) = real([ 4.316537D-02,  1.165552D-01,  4.299857D-02, -7.777778D+04,  5.684108D-03, -7.777778D+04,  0.000000D+00,  2.120821D-02,  5.169221D-02,  3.702557D-02,  5.169221D-02, -7.777778D+04,  4.505894D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -3.149018D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -3.149018D-02, -3.149018D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  4.444826D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -3.149018D-02,  5.684108D-03, -7.777778D+04, -3.149018D-02,  5.684108D-03, -3.149018D-02,  4.316537D-02,  4.316537D-02,  4.444826D-02,  5.684108D-03,  3.997699D-01,  1.167519D-02,  5.074305D-02, -3.149018D-02,  9.011788D-02,  1.167519D-02,  0.000000D+00], kind=wp)
    bTABnc(1:Nmaingroups,05) = real([ 9.447874D-02,  1.981737D-01,  9.068728D-02, -7.777778D+04, -2.496397D-02, -7.777778D+04,  0.000000D+00,  5.500988D-02,  1.968493D-01,  1.023700D-01,  1.968490D-01, -7.777778D+04,  9.357357D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -1.541675D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -1.541675D-02, -1.541675D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.969963D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -1.541675D-02, -2.496397D-02, -7.777778D+04, -1.541675D-02, -2.496397D-02, -1.541675D-02,  9.447874D-02,  9.447874D-02,  7.969963D-02, -2.496397D-02,  2.249583D-01,  7.906199D-02,  6.860960D-02, -1.541675D-02,  1.871471D-01,  7.906199D-02,  0.000000D+00], kind=wp)
    cTABnc(1:Nmaingroups,05) = real([ 6.568965D-02,  1.381514D-01,  2.230661D-02, -7.777778D+04,  7.786872D-05, -7.777778D+04,  0.000000D+00,  1.413451D-03, -4.275766D-03, -6.996550D-02, -4.275770D-03, -7.777778D+04,  6.145360D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.472152D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.472152D-02, -7.472152D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.134212D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.472152D-02,  7.786872D-05, -7.777778D+04, -7.472152D-02,  7.786872D-05, -7.472152D-02,  6.568965D-02,  6.568965D-02,  3.134212D-02,  7.786872D-05,  3.997699D-01, -9.031870D-03,  6.153147D-02, -7.472152D-02,  1.229072D-01, -9.031870D-03,  0.000000D+00], kind=wp)
    bTABnc(1:Nmaingroups,21) = real([ 1.021409D-01,  2.062126D-01,  8.109440D-02, -7.777778D+04,  1.038726D-02, -7.777778D+04,  0.000000D+00,  7.100129D-02,  2.044581D-01,  1.021191D-01,  2.017796D-01, -7.777778D+04,  8.254596D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.605447D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.605447D-02,  2.605447D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  8.092567D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.605447D-02,  1.038726D-02, -7.777778D+04,  2.605447D-02,  1.038726D-02,  2.605447D-02,  1.021409D-01,  1.021409D-01,  8.092567D-02,  1.038726D-02, -7.777778D+04,  1.281954D-01,  9.293322D-02,  2.605447D-02,  1.650919D-01,  1.281954D-01,  0.000000D+00], kind=wp)
    cTABnc(1:Nmaingroups,21) = real([ 7.474215D-02,  1.525628D-01,  5.931454D-02, -7.777778D+04, -5.465988D-04, -7.777778D+04,  0.000000D+00,  2.885781D-02,  1.524886D-01,  7.471416D-02, -1.077497D-01, -7.777778D+04,  7.077540D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.331719D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.331719D-02,  2.331719D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.926138D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.331719D-02, -5.465988D-04, -7.777778D+04,  2.331719D-02, -5.465988D-04,  2.331719D-02,  7.474215D-02,  7.474215D-02,  5.926138D-02, -5.465988D-04, -7.777778D+04,  9.805934D-02,  7.022880D-02,  2.331719D-02,  1.415508D-01,  9.805934D-02,  0.000000D+00], kind=wp)
    bTABnc(1:Nmaingroups,23) = real([ 8.761896D-02,  1.465395D-01,  7.522760D-02, -7.777778D+04,  3.103042D-03, -7.777778D+04,  0.000000D+00,  6.831109D-02,  1.442190D-01,  6.831109D-02,  1.442190D-01, -7.777778D+04,  7.252455D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.112358D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.112358D-02,  5.112358D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  6.749006D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.112358D-02,  3.103042D-03, -7.777778D+04,  5.112358D-02,  3.103042D-03,  5.112358D-02,  8.761896D-02,  8.761896D-02,  6.749006D-02,  3.103042D-03, -7.777778D+04,  1.387425D-01,  7.562759D-02,  5.112358D-02,  1.450491D-01,  1.387425D-01,  0.000000D+00], kind=wp)
    cTABnc(1:Nmaingroups,23) = real([-2.317547D-02, -1.395783D-02, -2.883615D-02, -7.777778D+04, -8.659318D-03, -7.777778D+04,  0.000000D+00, -4.267910D-02, -2.095079D-02, -4.267910D-02, -2.095079D-02, -7.777778D+04, -2.344034D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.495400D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.495400D-02, -2.495400D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.232628D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.495400D-02, -8.659318D-03, -7.777778D+04, -2.495400D-02, -8.659318D-03, -2.495400D-02, -2.317547D-02, -2.317547D-02, -2.232628D-02, -8.659318D-03, -7.777778D+04, -4.812947D-02, -3.209966D-02, -2.495400D-02, -4.688068D-02, -4.812947D-02,  0.000000D+00], kind=wp)

    
    !main group <--> anions:
    bTABna = -7.777778E+04_wp      !default value indicating that the parameter has not been determined
    cTABna = -7.777778E+04_wp
    !               main group_no         1             2                3             4               5              6             7              8            9              10              11            12              13           14              15          16               17              18              19              20           21              22             23           24              25            26             27              28             29             30             31             32               33            34             35             36             37             38            39              40             41             42             43             44             45             46             47             48             49             50             51             52             53             54             55             56             57             58             59             60             61             62             63             64             65             66             67             68             69             70             71             72             73             74             75                76
    !               main_groups        (CHn),            (C=C),         (ACHn),    (ACCHn),   (OH[stand_UNIFAC]), (CH3OH[methanol]), (H2O),        (ACOH),       (CHnCO),  (CHO[aldehyde]),   (CCOO),  (HCOO[formate]),  (CHnO[ether]),   (CNH2),        (CNH),        ((C)3N),         (ACNH2),     (PYRIEINE),   (CCN[nitrile]), (COOH[stand_UNIFAC]), (CCl),     (CCl2),          (CCl3),    (CCl4),         (ACCl),         (CNO2),        (ACNO2),       (CS2),      (CH3SH),      (FURFURAL),         (EOH),          (I),           (BR),   (C#C[triple_bond]), (EMSO),       (ACRY),        (ClCC),        (ACF),        (EMF),          (CF2),          (COO),      (SIH2),         (SIO),          (NMP),         (CClF),       (CON),        (OCCOH),        (CH2S),     (MORPHOLINE),   (THIOPHENE),   (Inorg_Ions),   (CH[OH]), (OH[MingRL]), (COOH[MingRL]), (CHnCO[MingRL]), (CHn[Ming_mo]), (CHnO[Ming_mo]), (OH[Ming_mo]), (CHn[Ming_hy]), (COOH[Ming_hy]),(OH[Ming_hy]),(CHn[Ming_di]),(COOH[Ming_di]), (OH[Peng]),  (COOH),      (CHn[alc]),   (CHn[alc-tail]),  (CH2[OH]),       (OH),    (CH2OCH2[PEG]),  (CHnONO2),  (CHnOOH[perox]),(C(=O)OOH[perox]),(CHnOOCHm[perox]),(C(=O)OONO2[perox]), (CO2)         
    bTABna(1:Nmaingroups,02) = real([ 6.914307D-02,  1.329562D-01,  5.290417D-02, -7.777778D+04,  9.156391D-03, -7.777778D+04,  0.000000D+00,  4.028553D-02,  1.329562D-01,  6.904875D-02,  1.327976D-01, -7.777778D+04,  5.425750D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.885736D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.885736D-02,  3.885736D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.227910D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.885736D-02,  9.156391D-03, -7.777778D+04,  3.885736D-02,  9.156391D-03,  3.885736D-02,  6.914307D-02,  6.914307D-02,  5.227910D-02,  9.156391D-03, -7.777778D+04,  1.080004D-01,  6.341389D-02,  3.885736D-02,  1.085150D-01,  1.080004D-01,  0.000000D+00], kind=wp)
    cTABna(1:Nmaingroups,02) = real([ 5.888996D-02,  1.363867D-01,  3.860389D-02, -7.777778D+04, -2.249250D-02, -7.777778D+04,  0.000000D+00,  3.767145D-02,  1.363866D-01,  5.886737D-02,  1.362981D-01, -7.777778D+04,  4.161487D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.009678D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.009678D-02,  3.009678D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.656014D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.009678D-02, -2.249250D-02, -7.777778D+04,  3.009678D-02, -2.249250D-02,  3.009678D-02,  5.888996D-02,  5.888996D-02,  3.656014D-02, -2.249250D-02, -7.777778D+04,  8.898674D-02,  1.912237D-02,  3.009678D-02,  8.322974D-02,  8.898674D-02,  0.000000D+00], kind=wp)
    bTABna(1:Nmaingroups,03) = real([ 4.136791D-02,  8.229121D-02,  4.128208D-02, -7.777778D+04,  4.686682D-03, -7.777778D+04,  0.000000D+00,  2.998976D-02,  6.480397D-02,  3.717467D-02,  6.475417D-02, -7.777778D+04,  4.133287D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.998974D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.998974D-02,  2.998974D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  4.094891D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  2.998974D-02,  4.686682D-03, -7.777778D+04,  2.998974D-02,  4.686682D-03,  2.998974D-02,  4.136791D-02,  4.136791D-02,  4.094891D-02,  4.686682D-03, -7.777778D+04,  7.135765D-02,  4.601955D-02,  2.998974D-02,  8.266574D-02,  7.135765D-02,  0.000000D+00], kind=wp)
    cTABna(1:Nmaingroups,03) = real([ 3.439384D-02,  5.607449D-02,  3.401668D-02, -7.777778D+04,  1.476176D-03, -7.777778D+04,  0.000000D+00,  1.237431D-02,  4.151743D-02,  1.807473D-02,  4.150937D-02, -7.777778D+04,  2.391176D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  1.237430D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  1.237430D-02,  1.237430D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  1.688331D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  1.237430D-02,  1.476176D-03, -7.777778D+04,  1.237430D-02,  1.476176D-03,  1.237430D-02,  3.439384D-02,  3.439384D-02,  1.688331D-02,  1.476176D-03, -7.777778D+04,  4.676814D-02,  2.538794D-02,  1.237430D-02,  4.782352D-02,  4.676814D-02,  0.000000D+00], kind=wp)
    bTABna(1:Nmaingroups,04) = real([ 4.878929D-02,  7.134253D-02,  4.259999D-02, -7.777778D+04,  1.939024D-02, -7.777778D+04,  0.000000D+00,  4.133419D-02,  7.376798D-02,  3.983930D-02,  3.427790D-02, -7.777778D+04,  4.117654D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.427825D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.427825D-02,  3.427825D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.427807D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.235025D-02,  1.939024D-02, -7.777778D+04,  3.235025D-02,  1.939024D-02,  3.235025D-02,  4.878929D-02,  4.878929D-02,  3.427807D-02,  1.939024D-02, -7.777778D+04,  8.306755D-02,  6.056677D-02,  3.427825D-02,  8.235307D-02,  8.306755D-02,  0.000000D+00], kind=wp)
    cTABna(1:Nmaingroups,04) = real([ 8.516072D-03, -1.301376D-02, -1.896587D-02, -7.777778D+04, -1.636921D-02, -7.777778D+04,  0.000000D+00, -2.434093D-02, -1.520187D-02, -2.517013D-02,  2.169341D-02, -7.777778D+04, -1.417372D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -3.647636D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -3.647636D-02, -3.647636D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -6.082930D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -1.342741D-02, -1.636921D-02, -7.777778D+04, -1.342741D-02, -1.636921D-02, -1.342741D-02, -8.516072D-03, -8.516072D-03, -6.082930D-03, -1.636921D-02, -7.777778D+04, -4.499244D-02, -3.054293D-02, -3.647636D-02, -2.834744D-02, -4.499244D-02,  0.000000D+00], kind=wp)
    bTABna(1:Nmaingroups,05) = real([ 4.233234D-02,  6.701021D-02,  4.209408D-02, -7.777778D+04, -2.707269D-02, -7.777778D+04,  0.000000D+00,  4.535367D-03,  4.542737D-02,  4.139841D-02,  4.542735D-02, -7.777778D+04,  4.773836D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.755684D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.755684D-03, -2.755684D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.000031D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.755684D-03, -2.707269D-02, -7.777778D+04, -2.755684D-03, -2.707269D-02, -2.755684D-03,  4.233234D-02,  4.233234D-02,  3.000031D-02, -2.707269D-02, -7.777778D+04,  3.957666D-02, -2.229885D-02, -2.755684D-03,  9.547672D-03,  3.957666D-02,  0.000000D+00], kind=wp)
    cTABna(1:Nmaingroups,05) = real([ 4.853253D-02,  1.230290D-01,  4.848432D-02, -7.777778D+04,  8.326951D-03, -7.777778D+04,  0.000000D+00,  2.880167D-02,  6.624150D-02,  4.875210D-02,  6.624148D-02, -7.777778D+04,  8.524851D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  6.663517D-04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  6.663517D-04,  6.663517D-04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.967023D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  6.663517D-04,  8.326951D-03, -7.777778D+04,  6.663517D-04,  8.326951D-03,  6.663517D-04,  4.853253D-02,  4.853253D-02,  3.967023D-02,  8.326951D-03, -7.777778D+04,  4.919888D-02,  9.357546D-02,  6.663517D-04,  1.704970D-01,  4.919888D-02,  0.000000D+00], kind=wp)
    bTABna(1:Nmaingroups,06) = real([ 4.233234D-02,  6.701021D-02,  4.209408D-02, -7.777778D+04, -2.707269D-02, -7.777778D+04,  0.000000D+00,  4.535367D-03,  4.542737D-02,  4.139841D-02,  4.542735D-02, -7.777778D+04,  4.773836D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.755684D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.755684D-03, -2.755684D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.000031D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -2.755684D-03, -2.707269D-02, -7.777778D+04, -2.755684D-03, -2.707269D-02, -2.755684D-03,  4.233234D-02,  4.233234D-02,  3.000031D-02, -2.707269D-02, -7.777778D+04,  3.957666D-02, -2.229885D-02, -2.755684D-03,  9.547672D-03,  3.957666D-02,  0.000000D+00], kind=wp)
    cTABna(1:Nmaingroups,06) = real([ 4.853253D-02,  1.230290D-01,  4.848432D-02, -7.777778D+04,  8.326951D-03, -7.777778D+04,  0.000000D+00,  2.880167D-02,  6.624150D-02,  4.875210D-02,  6.624148D-02, -7.777778D+04,  8.524851D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  6.663517D-04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  6.663517D-04,  6.663517D-04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.967023D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  6.663517D-04,  8.326951D-03, -7.777778D+04,  6.663517D-04,  8.326951D-03,  6.663517D-04,  4.853253D-02,  4.853253D-02,  3.967023D-02,  8.326951D-03, -7.777778D+04,  4.919888D-02,  9.357546D-02,  6.663517D-04,  1.704970D-01,  4.919888D-02,  0.000000D+00], kind=wp)
    bTABna(1:Nmaingroups,07) = real([ 4.878929D-02,  7.134253D-02,  4.259999D-02, -7.777778D+04,  1.939024D-02, -7.777778D+04,  0.000000D+00,  4.133419D-02,  7.376798D-02,  3.983930D-02,  3.427790D-02, -7.777778D+04,  4.117654D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.427825D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.427825D-02,  3.427825D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.427807D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  3.235025D-02,  1.939024D-02, -7.777778D+04,  3.235025D-02,  1.939024D-02,  3.235025D-02,  4.878929D-02,  4.878929D-02,  3.427807D-02,  1.939024D-02, -7.777778D+04,  8.306755D-02,  6.056677D-02,  3.427825D-02,  8.235307D-02,  8.306755D-02,  0.000000D+00], kind=wp)
    cTABna(1:Nmaingroups,07) = real([ 8.516072D-03, -1.301376D-02, -1.896587D-02, -7.777778D+04, -1.636921D-02, -7.777778D+04,  0.000000D+00, -2.434093D-02, -1.520187D-02, -2.517013D-02,  2.169341D-02, -7.777778D+04, -1.417372D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -3.647636D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -3.647636D-02, -3.647636D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -6.082930D-03, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -1.342741D-02, -1.636921D-02, -7.777778D+04, -1.342741D-02, -1.636921D-02, -1.342741D-02, -8.516072D-03, -8.516072D-03, -6.082930D-03, -1.636921D-02, -7.777778D+04, -4.499244D-02, -3.054293D-02, -3.647636D-02, -2.834744D-02, -4.499244D-02,  0.000000D+00], kind=wp)
    bTABna(1:Nmaingroups,08) = real([ 1.244526D-01,  2.299447D-01,  9.637726D-02, -7.777778D+04, -6.100228D-02, -7.777778D+04,  0.000000D+00,  4.375717D-02,  1.728215D-01,  7.159264D-02,  1.808274D-01, -7.777778D+04,  6.071562D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02,  7.999459D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  9.617478D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02, -6.100228D-02, -7.777778D+04,  7.999459D-02, -6.100228D-02,  7.999459D-02,  1.244526D-01,  1.244526D-01,  9.617478D-02, -6.100228D-02, -7.777778D+04,  2.044472D-01, -2.866600D-04,  7.999459D-02,  1.214312D-01,  2.044472D-01,  0.000000D+00], kind=wp)
    cTABna(1:Nmaingroups,08) = real([ 1.426790D-01,  1.962500D-01,  9.774756D-02, -7.777778D+04, -3.952179D-02, -7.777778D+04,  0.000000D+00,  3.042278D-03,  1.053198D-01,  2.249248D-02,  8.572575D-02, -7.777778D+04,  2.439568D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06,  5.620295D-06, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  9.757109D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06, -3.952179D-02, -7.777778D+04,  5.620295D-06, -3.952179D-02,  5.620295D-06,  1.426790D-01,  1.426790D-01,  9.757109D-02, -3.952179D-02, -7.777778D+04,  1.426846D-01, -1.512611D-02,  5.620295D-06,  4.879136D-02,  1.426846D-01,  0.000000D+00], kind=wp)     
    bTABna(1:Nmaingroups,09) = real([ 1.244526D-01,  2.299447D-01,  9.637726D-02, -7.777778D+04, -6.100228D-02, -7.777778D+04,  0.000000D+00, -7.777778D+04,  2.333990D-01,  1.089460D-01,  2.333990D-01, -7.777778D+04,  6.071562D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02,  7.999459D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  9.617478D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02, -6.100228D-02, -7.777778D+04,  7.999459D-02, -6.100228D-02,  7.999459D-02,  1.244526D-01,  1.244526D-01,  9.617478D-02, -6.100228D-02, -7.777778D+04,  2.044472D-01, -2.866600D-04,  7.999459D-02,  1.214312D-01,  2.044472D-01,  0.000000D+00], kind=wp)
    cTABna(1:Nmaingroups,09) = real([ 1.426790D-01,  1.962500D-01,  9.774756D-02, -7.777778D+04, -3.952179D-02, -7.777778D+04,  0.000000D+00, -7.777778D+04,  2.234230D-01,  8.074430D-02,  2.234230D-01, -7.777778D+04,  2.439568D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06,  5.620295D-06, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  9.757109D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06, -3.952179D-02, -7.777778D+04,  5.620295D-06, -3.952179D-02,  5.620295D-06,  1.426790D-01,  1.426790D-01,  9.757109D-02, -3.952179D-02, -7.777778D+04,  1.426846D-01, -1.512611D-02,  5.620295D-06,  4.879136D-02,  1.426846D-01,  0.000000D+00], kind=wp)
    bTABna(1:Nmaingroups,10) = real([ 1.244526D-01,  2.299447D-01,  9.637726D-02, -7.777778D+04, -6.100228D-02, -7.777778D+04,  0.000000D+00,  4.375717D-02,  1.728215D-01,  7.159264D-02,  1.808274D-01, -7.777778D+04,  6.071562D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02,  7.999459D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  9.617478D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.999459D-02, -6.100228D-02, -7.777778D+04,  7.999459D-02, -6.100228D-02,  7.999459D-02,  1.244526D-01,  1.244526D-01,  9.617478D-02, -6.100228D-02, -7.777778D+04,  2.044472D-01, -2.866600D-04,  7.999459D-02,  1.214312D-01,  2.044472D-01,  0.000000D+00], kind=wp)
    cTABna(1:Nmaingroups,10) = real([ 1.426790D-01,  1.962500D-01,  9.774756D-02, -7.777778D+04, -3.952179D-02, -7.777778D+04,  0.000000D+00,  3.042278D-03,  1.053198D-01,  2.249248D-02,  8.572575D-02, -7.777778D+04,  2.439568D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06,  5.620295D-06, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  9.757109D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.620295D-06, -3.952179D-02, -7.777778D+04,  5.620295D-06, -3.952179D-02,  5.620295D-06,  1.426790D-01,  1.426790D-01,  9.757109D-02, -3.952179D-02, -7.777778D+04,  1.426846D-01, -1.512611D-02,  5.620295D-06,  4.879136D-02,  1.426846D-01,  0.000000D+00], kind=wp)
    bTABna(1:Nmaingroups,21) = real([ 7.598433D-02,  1.241119D-01,  7.056449D-02, -7.777778D+04,  2.271222D-03, -7.777778D+04,  0.000000D+00,  3.663588D-02,  1.238882D-01,  3.904856D-02,  1.726755D-02, -7.777778D+04,  4.690451D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -4.482039D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -4.482039D-02, -4.482039D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.257356D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -4.482039D-02,  2.271222D-03, -7.777778D+04, -4.482039D-02,  2.271222D-03, -4.482039D-02,  7.598433D-02,  7.598433D-02,  5.257356D-02,  2.271222D-03, -3.803980D-01,  3.116394D-02,  4.917573D-02, -4.482039D-02,  9.380902D-02,  3.116394D-02,  0.000000D+00], kind=wp)
    cTABna(1:Nmaingroups,21) = real([ 1.130960D-01,  2.717248D-01,  6.332340D-02, -7.777778D+04,  1.224391D-02, -7.777778D+04,  0.000000D+00,  6.143489D-02,  1.346749D-01,  7.196915D-02,  1.077373D-03, -7.777778D+04,  8.606047D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -5.001442D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -5.001442D-02, -5.001442D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.931526D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -5.001442D-02,  1.224391D-02, -7.777778D+04, -5.001442D-02,  1.224391D-02, -5.001442D-02,  1.130960D-01,  1.130960D-01,  7.931526D-02,  1.224391D-02, -3.296446D-01,  6.308158D-02,  9.830438D-02, -5.001442D-02,  1.721209D-01,  6.308158D-02,  0.000000D+00], kind=wp)
    !CO3--
    bTABna(1:Nmaingroups,22) = real([ 7.598433D-02,  1.241119D-01,  7.056449D-02, -7.777778D+04,  2.271222D-03, -7.777778D+04,  0.000000D+00,  3.663588D-02,  1.238882D-01,  3.904856D-02,  1.726755D-02, -7.777778D+04,  4.690451D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -4.482039D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -4.482039D-02, -4.482039D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  5.257356D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -4.482039D-02,  2.271222D-03, -7.777778D+04, -4.482039D-02,  2.271222D-03, -4.482039D-02,  7.598433D-02,  7.598433D-02,  5.257356D-02,  2.271222D-03, -3.803980D-01,  3.116394D-02,  4.917573D-02, -4.482039D-02,  9.380902D-02,  3.116394D-02,  0.000000D+00], kind=wp)
    cTABna(1:Nmaingroups,22) = real([ 1.130960D-01,  2.717248D-01,  6.332340D-02, -7.777778D+04,  1.224391D-02, -7.777778D+04,  0.000000D+00,  6.143489D-02,  1.346749D-01,  7.196915D-02,  1.077373D-03, -7.777778D+04,  8.606047D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -5.001442D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -5.001442D-02, -5.001442D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04,  7.931526D-02, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -7.777778D+04, -5.001442D-02,  1.224391D-02, -7.777778D+04, -5.001442D-02,  1.224391D-02, -5.001442D-02,  1.130960D-01,  1.130960D-01,  7.931526D-02,  1.224391D-02, -3.296446D-01,  6.308158D-02,  9.830438D-02, -5.001442D-02,  1.721209D-01,  6.308158D-02,  0.000000D+00], kind=wp)  

    !------ cation <-->  anion interactions (in water) ----------------------
    !bkion and ckion interaction coefficients estimated by AIOMFAC model fit to experimental data:
    !further interaction parameters are set in MRfitpar, but only during fit of parameters.  

    ! Data for the ionic part of the midrange interaction term
    ! Determined by AZ in 2006, 2007 (4 - 5 parameter fit with omega = 0.8, omega2 = variable).
    ! bTABAC(I,J): contains the b values of the midrange terms. I: cations, J: anions.
    ! cTABAC(I,J): contains the c values of the midrange terms
    ! bAC(I,J): interaction parameter b between cation I and anion J.
    ! cAC(I,J): interaction parameter c between cation I and anion J.
    allocate( Cn1TABAC(40,40), Cn2TABAC(40,40), omega2TAB(40,40), bTABAC(40,40), cTABAC(40,40), &
        & TABhighestWTF(40,40), TABKsp(40,40), omegaTAB(40,40) )
    
    !initialize:
    bTABAC = 0.0E0_wp
    cTABAC = 0.0E0_wp
    cn1TABAC = 0.0E0_wp
    cn2TABAC = 0.0E0_wp
    omega2TAB = 0.6E0_wp !arbitrary value
    omegaTAB = 0.8E0_wp  !the usual value
    TABhighestWTF = 0.0E0_wp
    TABKsp = 1.0E30_wp !default value on [molal basis] (indicating an extremely high solubility product as default)
  
    !Li+ <-> Cl-          4P AUG 2007      
    bTABAC(1,2) = 0.106554911884921E0_wp        
    cTABAC(1,2) = 0.206369566797954E0_wp
    cn1TABAC(1,2) = 0.0E0_wp
    cn2TABAC(1,2) = 5.323941573688994E-02_wp  
    omega2TAB(1,2) = 0.535548024591189E0_wp     
    TABhighestwtf(1,2) = 0.45E0_wp
    TABKsp(1,2) = 1.4815E6_wp ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !Li+ <-> Br-          4P AUG 2007      
    bTABAC(1,3) = 0.106384115653176E0_wp 
    cTABAC(1,3) = 0.316479521480175E0_wp    
    cn1TABAC(1,3) = 0.0E0_wp
    cn2TABAC(1,3) = 5.760198325796455E-02_wp    
    omega2TAB(1,3) = 0.464658335659553E0_wp      
    TABhighestwtf(1,3) = 0.35E0_wp
    
    !Li+ <-> I- 
    bTABAC(1,4) =  1.8588934753695888E-01_wp
    cTABAC(1,4) =  6.0808982141071000E-01_wp   
    cn1TABAC(1,4) = 0.0E0_wp
    cn2TABAC(1,4) = 4.1814140083157633E-02_wp
    omega2TAB(1,4) = 4.1341524514949507E-01_wp
    TABKsp(1,4) = 1.107E06_wp ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !Li+ <-> NO3-         4P JUN 2007      
    bTABAC(1,5) = 7.631270368216216E-02_wp
    cTABAC(1,5) = 0.300549960057599E0_wp
    cn1TABAC(1,5) = 0.0E0_wp
    cn2TABAC(1,5) = 4.670051954857116E-02_wp
    omega2TAB(1,5) = 0.664928462062345
    TABhighestwtf(1,5) = 0.58E0_wp               
    
    !Li+ <-> IO3-     May 2021
    bTABAC(1,6) = 7.6049825043886410E-02_wp
    cTABAC(1,6) = -2.4930506321164203E-01_wp
    omegaTAB(1,6) = 2.9689011671088517E-01_wp
    cn1TABAC(1,6) = 0.0E0_wp
    cn2TABAC(1,6) = -1.2852865367456584E-01_wp
    omega2TAB(1,6) = 1.7423828299201640E0_wp
    TABhighestwtf(1,6) = 0.16E0_wp !NOt true value
    
    !Li+ <-> OH-      June 2021
    bTABAC(1,7) = 7.5665164617218825E-02_wp
    cTABAC(1,7) = -9.1330575713229178E-01_wp
    cn1TABAC(1,7) = 0.0E0_wp
    cn2TABAC(1,7) = 6.9105686988206494E-01_wp
    omega2TAB(1,7) =  2.3477855868588415E0_wp
    
    !Li+ <-> HSO4-    May 2021
    bTABAC(1,8) = 2.8878992634118333E-01_wp
    cTABAC(1,8) = -2.5337206136122727E0_wp
    cn1TABAC(1,8) = 0.0E0_wp 
    cn2TABAC(1,8) = 6.7541421521779388E-01_wp
    omega2TAB(1,8) = 1.5684706739420728E0_wp
    
    !Li+ <-> HCO3-   estimated by K+ <-> HCO3-
    bTABAC(1,10) =  9.2903436457389355E-02_wp
    cTABAC(1,10) =  -8.8237465618429356E-01_wp
    omegaTAB(1,10) = 5.6552279057559085E-01_wp
    cn1TABAC(1,10) = 0.0E0_wp
    cn2TABAC(1,10) = 5.9250512478964340E-03_wp
    omega2TAB(1,10) = 6.2492599123964832E-01_wp

    !Li+ <-> SO4--        4P JUN 2007  
    bTABAC(1,21) = 0.114470419384503
    cTABAC(1,21) = 3.540074251345881E-02_wp
    cn1TABAC(1,21) = 0.0E0_wp
    cn2TABAC(1,21) = -0.263257805913055E0_wp 
    omega2TAB(1,21) = 1.31696686093901E0_wp
    TABhighestwtf(1,21) = 0.25E0_wp
    TABKsp(1,21) = 3.25E0_wp ![molal basis] molal ion activity product at the solubility limit @ 298.15 K  

    !Li+ <-> CO3--        Oct 2020    
    bTABAC(1,22) = 7.5165264810580406E-04_wp
    cTABAC(1,22) = -1.6045235322580549E0_wp
    cn1TABAC(1,22) = 0.0E0_wp
    cn2TABAC(1,22) = -9.5120991084252995E-03_wp
    omega2TAB(1,22) =  3.6525829891663431E-01_wp
    TABKsp(1,22) = 3.493E-03_wp ![molal basis] molal ion activity product at the solubility limit @ 298.15 K

    !Na+ <-> Cl-          4P JUN 2007         
    bTABAC(2,2) = 5.374079546760321E-02_wp 
    cTABAC(2,2) = 7.977128474886985E-02_wp   
    cn1TABAC(2,2) = 0.0E0_wp
    cn2TABAC(2,2) = 2.455349637351519E-02_wp  
    omega2TAB(2,2) = 0.562981194934520E0_wp       
    TABhighestwtf(2,2) = 0.44E0_wp
    TABKsp(2,2) = 37.83E0_wp ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !Na+ <-> Br-          4P JUN 2007        
    bTABAC(2,3) = 0.180807320471190E0_wp
    cTABAC(2,3) = 0.273143813534931E0_wp   
    cn1TABAC(2,3) = 0.0E0_wp
    cn2TABAC(2,3) = -0.506597565983690E0_wp 
    omega2TAB(2,3) = 2.20904992991060E0_wp  
    TABhighestwtf(2,3) = 0.49E0_wp 
    TABKsp(2,3) = 322.0E0_wp ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !Na+ <-> I-            
    bTABAC(2,4) = 6.834973E-02_wp 
    cTABAC(2,4) = 3.796614E-01_wp   
    cn1TABAC(2,4) = 0.0E0_wp
    cn2TABAC(2,4) = 4.145436E-02_wp
    omega2TAB(2,4) = 4.994560E-01_wp  
    TABKsp(2,4) = 5.741E03_wp ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !Na+ <-> NO3-        5P JUN 2007
    bTABAC(2,5) = 1.164350105903698E-03_wp
    cTABAC(2,5) = -0.102545526130049E0_wp 
    omegaTAB(2,5) = 0.410452675064302E0_wp 
    cn1TABAC(2,5) = 0.0E0_wp
    cn2TABAC(2,5) = 2.534857925520951E-03_wp 
    omega2TAB(2,5) = 0.512657343426316E0_wp   
    TABhighestwtf(2,5) = 0.80E0_wp
    TABKsp(2,5) = 11.82E0_wp ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !Na+ <-> IO3-     May 2021
    bTABAC(2,6) = 1.3437657607868917E-02_wp
    cTABAC(2,6) = -9.6486635628713002E-01_wp
    cn1TABAC(2,6) = 0.0E0_wp
    cn2TABAC(2,6) = -1.3706943474700559E-01_wp
    omega2TAB(2,6) =  1.3626195116594571E0_wp
    omegaTAB(2,6) =   1.3047771888155979E0_wp

    !Na+ <-> OH-      Mar 2020
    bTABAC(2,7) = 2.3840333396438568E-02_wp
    cTABAC(2,7) = 6.6956349269406479E-01_wp
    omegaTAB(2,7) = 2.2686902829375624E-01_wp
    cn1TABAC(2,7) = 0.0E0_wp
    cn2TABAC(2,7) = -5.2304938122407685E-01_wp
    omega2TAB(2,7) = 1.3037546443208794E0_wp

    !Na+ <-> HSO4-       5P FEB 2011           
    bTABAC(2,8) = 1.5321434237810060E-02_wp
    cTABAC(2,8) = 4.0E-01_wp
    omegaTAB(2,8) = 4.2363461984697831E-01_wp 
    cn1TABAC(2,8) = 0.0E0_wp    
    cn2TABAC(2,8) = 3.5007219661690844E-03_wp   
    omega2TAB(2,8) = 4.0E-01_wp 
    TABhighestwtf(2,8) = 0.15E0_wp   
    
    !Na+ <-> HCO3-  Jan 2021
    bTABAC(2,10) =  6.4998964390762845E-02_wp
    cTABAC(2,10) =  -7.5189119274972557E-02_wp
    cn1TABAC(2,10) = 0.0E0_wp
    cn2TABAC(2,10) = -5.7300345498924488E-02_wp
    omega2TAB(2,10) = 9.3276066630078436E-01_wp
    
    !Na+ <-> SO4--       4P JUN 2007           
    bTABAC(2,21) = 1.890680686105647E-03_wp
    cTABAC(2,21) = -0.424184309855867E0_wp
    cn1TABAC(2,21) = 0.0E0_wp
    cn2TABAC(2,21) = -0.223851238896391E0_wp     
    omega2TAB(2,21) = 1.05361978864913E0_wp  
    TABhighestwtf(2,21) = 0.60E0_wp
    !TABKsp(2,21) = 0.1236E0_wp !@ 298 K; [molal basis] molal ion activity product at the solubility limit @ 298 K (these salt or hydrates have a strong temperature dependence!).   
        !value based on aqueous Na2SO4 solubility at 293 and 303 K from Okorafor (1999) (actually it forms a hydrate Na2SO4*10 H2O at T < 32 deg C!
    TABKsp(2,21) =  0.4877E0_wp !@ 308 K; [molal basis] molal ion activity product at the solubility limit @ 308 K, based on data of anhydrous Na2SO4 solubility at 308 K by Vener and Thompson (1949).                
    
    !Na+ <-> CO3--    Oct 2020          
    bTABAC(2,22) =  1.0695833656549962E-01_wp
    cTABAC(2,22) =  -3.8465749714565167E-01_wp
    cn1TABAC(2,22) = 0.0E0_wp
    cn2TABAC(2,22) = -1.9818447835032355E-01_wp
    omega2TAB(2,22) =  8.4635302771576537E-01_wp
    
    !K+ <-> Cl-        4P JUL 2007
    bTABAC(3,2) = 1.656079143002206E-02_wp       
    cTABAC(3,2) = -2.752209560273460E-03_wp
    cn1TABAC(3,2) = 0.0E0_wp
    cn2TABAC(3,2) = 2.083280655191561E-02_wp
    omega2TAB(3,2) = 0.670529802215773E0_wp
    TABhighestwtf(3,2) = 0.49E0_wp
    TABKsp(3,2) = 8.038E0_wp ![molal basis] molal ion activity product at the solubility limit @ 298 K (Pinho and Macedo, 2005)

    !K+ <-> Br-          4P JUN 2007        
    bTABAC(3,3) = 3.368790641502646E-02_wp 
    cTABAC(3,3) = 6.088180362627697E-02_wp 
    cn1TABAC(3,3) = 0.0E0_wp 
    cn2TABAC(3,3) = 1.529324768456540E-02_wp 
    omega2TAB(3,3) = 0.565063111593592E0_wp
    TABhighestwtf(3,3) = 0.40E0_wp
    TABKsp(3,3) = 11.34E0_wp ![molal basis] molal ion activity product at the solubility limit @ 298 K
    
    !K+ <-> I-            
    bTABAC(3,4) = 5.442088E-02_wp 
    cTABAC(3,4) = 1.882261E-01_wp
    cn1TABAC(3,4) = 0.0E0_wp
    cn2TABAC(3,4) = 2.534565E-02_wp 
    omega2TAB(3,4) = 6.607733E-01_wp 
    TABKsp(3,4) = 54.11E0_wp ![molal basis] molal ion activity product at the solubility limit @ 298 K

    !K+ <-> NO3-          5P JUN 2007
    bTABAC(3,5) = 2.498205354734865E-05_wp  
    cTABAC(3,5) =  -0.413171796795632E0_wp   
    cn1TABAC(3,5) = 0.0E0_wp
    cn2TABAC(3,5) = -4.551039248054616E-04_wp 
    omega2TAB(3,5) = 0.242243640032790E0_wp   
    omegaTAB(3,5) = 0.357227451106132E0_wp  
    TABhighestwtf(3,5) = 0.26E0_wp 
    TABKsp(3,5) = 0.9225E0_wp ![molal basis] molal ion activity product at the solubility limit @ 298 K 
    
    !K+ <-> IO3-    ! April 2021
    bTABAC(3,6) = 1.2531824958095275E-02_wp
    cTABAC(3,6) = -1.8395581736986735E0_wp
    cn1TABAC(3,6) = 0.0E0_wp
    cn2TABAC(3,6) = -1.4618851585020370E-01_wp
    omega2TAB(3,6) =  1.3612241273448702E0_wp
    omegaTAB(3,6) = 1.2744128860114006E0_wp

    !K+ <-> OH-      ! Mar 2020
    bTABAC(3,7) = 2.5131496199052189E-01_wp
    cTABAC(3,7) =  5.3647067717101127E-01_wp
    cn1TABAC(3,7) = 0.0E0_wp
    cn2TABAC(3,7) = -7.5910080641938760E-01_wp
    omega2TAB(3,7) = 1.9363546657176669E0_wp
    TABKsp(3,7) = 999.0E0_wp ![not true value]! 
    
    !K+ <-> HSO4-    Revised May 2019
    bTABAC(3,8) = 8.4433898901815352E-02_wp
    cTABAC(3,8) = -3.9301228443995739E-01_wp
    cn1TABAC(3,8) = 0.0E0_wp 
    cn2TABAC(3,8) = -9.0232965273929300E-01_wp
    omega2TAB(3,8) = 1.8692081240631140E0_wp

    !K+ <-> SO4--     Revised May 2019 
    bTABAC(3,21) =  6.6714599718481538E-02_wp
    cTABAC(3,21) = -1.0844118398860187E0_wp
    cn1TABAC(3,21) = 0.0E0_wp
    cn2TABAC(3,21) = -6.6511262694399914E-02_wp
    omega2TAB(3,21) = 5.8327209615736464E-01_wp   
    TABhighestwtf(3,21) = 0.11E0_wp  
    TABKsp(3,21) = 1.571E-02_wp ![molal basis]
    
    !K+ <-> HCO3-   Jan 2021
    bTABAC(3,10) =  9.2903436457389355E-02_wp
    cTABAC(3,10) =  -8.8237465618429356E-01_wp
    omegaTAB(3,10) = 5.6552279057559085E-01_wp
    cn1TABAC(3,10) = 0.0E0_wp
    cn2TABAC(3,10) = 5.9250512478964340E-03_wp
    omega2TAB(3,10) = 6.2492599123964832E-01_wp
    
    !K+ <-> CO3--   Oct 2020
    bTABAC(3,22) =  1.9268078927583168E-01_wp
    cTABAC(3,22) =  1.3547830332027755E-01_wp
    cn1TABAC(3,22) = 0.0E0_wp
    cn2TABAC(3,22) = -1.8190589684462621E-01_wp
    omega2TAB(3,22) = 9.4333743510918933E-01_wp

    !NH4+ <-> Cl-      5P JUN 2007  
    bTABAC(4,2) = 1.520261933644218E-03_wp
    cTABAC(4,2) = 4.907431846455088E-02_wp  
    cn1TABAC(4,2) = 0.0E0_wp    
    cn2TABAC(4,2) = 1.111220056439640E-02_wp
    omega2TAB(4,2) = 0.653256017313056E0_wp  
    TABhighestwtf(4,2) = 0.54E0_wp 
    omegaTAB(4,2) = 0.116800561298130  
    TABKsp(4,2) = 17.2E0_wp ![molal basis] 

    !NH4+ <-> Br-     5P JUL 2007        
    bTABAC(4,3) = 2.498124477602936E-03_wp
    cTABAC(4,3) = 8.151172908753296E-02_wp 
    cn1TABAC(4,3) = 0.0E0_wp
    cn2TABAC(4,3) = 1.379527431149357E-02_wp  
    omega2TAB(4,3) = 0.728984320385571E0_wp
    TABhighestwtf(4,3) = 0.43E0_wp
    omegaTAB(4,3) = 0.143620623302905E0_wp 
    TABKsp(4,3) = 24.0E0_wp ![molal basis]

    !NH4+ <-> I-            
    bTABAC(4,4) = 9.2191676199143202E-02_wp
    cTABAC(4,4) = -4.6783666750062035E-02_wp
    cn1TABAC(4,4) = 0.0E0_wp
    cn2TABAC(4,4) =  8.1996177814860266E-02_wp
    omega2TAB(4,4) =  1.3240421566062248E0_wp
    TABKsp(4,4) = 74.07E0_wp  ![molal basis]

    !NH4+ <-> NO3-    5P  JUN 2007  
    bTABAC(4,5) = -5.735407413065387E-05_wp
    cTABAC(4,5) = -0.171746228783622E0_wp  
    cn1TABAC(4,5) = 0.0E0_wp       
    cn2TABAC(4,5) = 5.510450416893165E-03_wp   
    omega2TAB(4,5) = 0.529761886212099E0_wp     
    TABhighestwtf(4,5) = 0.90E0_wp 
    omegaTAB(4,5) = 0.26E0_wp  
    TABKsp(4,5) = 12.0E0_wp ![molal basis]
    
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
    bTABAC(4,8) = 7.59734841E-03_wp
    cTABAC(4,8) = 1.43011653E-01_wp
    omegaTAB(4,8) = 2.03954033E-01_wp
    cn1TABAC(4,8) = 0.0E0_wp 
    cn2TABAC(4,8) = 6.31184288E-03_wp
    omega2TAB(4,8) = 8.25385786E-01_wp
    TABhighestwtf(4,8) = 0.90E0_wp 

    !NH4+ <-> HCO3-    Jan 2021
    bTABAC(4,10) = 1.7140258422482357E-01_wp
    cTABAC(4,10) = -2.1046493522524873E-01_wp
    omegaTAB(4,10) = 4.3530915388705860E-01_wp
    cn1TABAC(4,10) = 0.0E0_wp 
    cn2TABAC(4,10) = -4.5619906309823675E-02_wp
    omega2TAB(4,10) =  5.8023505996780234E-01_wp
    
    !(NH4)2SO4  [resp. NH4+ <-> SO4--]  5P JUN 2007 
    bTABAC(4,21) = 3.729564498700757E-04_wp
    cTABAC(4,21) = -0.906075204598434E0_wp
    omegaTAB(4,21) = 0.54510882832E0_wp
    cn1TABAC(4,21) = 0.0E0_wp
    cn2TABAC(4,21) = -3.790519418639869E-04_wp  
    omega2TAB(4,21) = 0.354206472573787E0_wp
    TABhighestwtf(4,21) = 0.78E0_wp
    TABKsp(4,21) = 1.03E0_wp ![molal basis]  

    !NH4+ <-> CO3--     Sept 2020
    bTABAC(4,22) = 7.8317937289243134E-03
    cTABAC(4,22) = -2.7115951835499486E0_wp
    cn1TABAC(4,22) = 0.0E0_wp 
    cn2TABAC(4,22) = -7.2643506216563458E-02_wp
    omega2TAB(4,22) = 5.3625142404031889E-01_wp

    !H+ <-> Cl-      4P JUN 2007
    bTABAC(5,2) = 0.182003163561796E0_wp
    cTABAC(5,2) = 0.243340062772300E0_wp  
    cn1TABAC(5,2) = 0.0E0_wp
    cn2TABAC(5,2) = 3.331893957013864E-02_wp 
    omega2TAB(5,2) = 0.504672383721131E0_wp  
    TABhighestwtf(5,2) = 0.18E0_wp

    !H+ <-> Br-      4P JUN 2007
    bTABAC(5,3) = 0.120325353682262E0_wp 
    cTABAC(5,3) = 0.444858839503779E0_wp 
    cn1TABAC(5,3) = 0.0E0_wp
    cn2TABAC(5,3) = 8.076729905421279E-02_wp
    omega2TAB(5,3) = 0.596775891883278E0_wp
    TABhighestwtf(5,3) = 0.20E0_wp  

    !H+ <-> I-  
    bTABAC(5,4) = 5.091034E-01_wp 
    cTABAC(5,4) = 1.793665E-01_wp
    cn1TABAC(5,4) = 0.0E0_wp
    cn2TABAC(5,4) = -9.943891E-02_wp
    omega2TAB(5,4) = 1.371292E0_wp  

    !H+ <-> NO3-     4P JUN 2007
    bTABAC(5,5) = 0.210637923183262E0_wp 
    cTABAC(5,5) = 0.122694135055374E0_wp
    cn1TABAC(5,5) = 0.0E0_wp
    cn2TABAC(5,5) = -0.101735827610270E0_wp
    omega2TAB(5,5) = 1.67641985798924E0_wp
    TABhighestwtf(5,5) = 0.16E0_wp

    !H+ <-> IO3-     July 2021
    bTABAC(5,6) = 1.4981114516279127E-02_wp
    cTABAC(5,6) = -3.2063552123769585E0_wp
    cn1TABAC(5,6) = 0.0E0_wp
    cn2TABAC(5,6) = -4.2730539646547307E-02_wp
    omega2TAB(5,6) = 5.3596305879882622E-01_wp

    !H+ <-> OH- 
    bTABAC(5,7) = 0.0E0_wp
    cTABAC(5,7) =  0.0E0_wp
    cn1TABAC(5,7) = 0.0E0_wp
    cn2TABAC(5,7) = 0.0E0_wp
    omega2TAB(5,7) = 0.0E0_wp
    omegaTAB(5,7) = 0.0E0_wp

    !H2SO4 (-> HSO4- & H+) 5P FEB 2011
    bTABAC(5,8) = 2.15532299E-02_wp  
    cTABAC(5,8) = 5.62965674E-01_wp      
    omegaTAB(5,8) = 1.42442019E-01_wp !!   
    cn1TABAC(5,8) = 0.0E0_wp
    cn2TABAC(5,8) = 7.03842440E-02_wp
    omega2TAB(5,8) = 7.14194282E-01_wp      
    TABhighestwtf(5,8) = 0.82E0_wp

    !H2CO3 (-> H+ & HCO3- )  Jan 2021
    bTABAC(5,10) = 4.8650377643706250E-02_wp
    cTABAC(5,10) = -4.7365194558795232E-01_wp   
    cn1TABAC(5,10) = 0.0E0_wp
    cn2TABAC(5,10) = 9.9088551465819994E-03_wp
    omega2TAB(5,10) = 3.2804860454964357E-01_wp  
    
    !HSO4-  (-> H+ & SO4--) 5P FEB 2011
    bTABAC(5,21) = 2.863428226E-01_wp
    cTABAC(5,21) = -5.996148056E0_wp   
    omegaTAB(5,21) = 1.368612639E0_wp 
    cn1TABAC(5,21) = 0.0E0_wp
    cn2TABAC(5,21) = -5.359772603E-01_wp                 
    omega2TAB(5,21) = 9.071999765E-01_wp 
    TABhighestwtf(5,21) = 0.82E0_wp

    !HCO3-  (-> H+ & CO3--)  Jan 2021
    bTABAC(5,22) = -1.5220470327665464E-01_wp 
    cTABAC(5,22) = 1.7020149536615228E-01_wp  
    cn1TABAC(5,22) = 0.0E0_wp
    cn2TABAC(5,22) = 8.6856221337292705E-02_wp
    omega2TAB(5,22) = 5.9622219134548793E-01_wp
    
    !Ca++ <-> Cl-         4P JUN 2007
    bTABAC(21,2) = 0.104920450619807E0_wp     
    cTABAC(21,2) = 0.866923035893255E0_wp    
    cn1TABAC(21,2) = 0.0E0_wp
    cn2TABAC(21,2) = 7.206272036095257E-02_wp    
    omega2TAB(21,2) = 0.365747420374652E0_wp   
    TABhighestwtf(21,2) = 0.53E0_wp

    !Ca++ <-> Br-         4P JUL 2010
    bTABAC(21,3) = 0.890929E0_wp
    cTABAC(21,3) = 6.101342E-02_wp
    cn1TABAC(21,3) = 0.0E0_wp
    cn2TABAC(21,3) = -0.238788E0_wp   
    omega2TAB(21,3) = 0.762961E0_wp
    TABhighestwtf(21,3) = 0.54E0_wp

    !Ca++ <-> I- 
    bTABAC(21,4) = 7.5976411109783837E-01_wp
    cTABAC(21,4) = 1.1031941621629198E0_wp  
    cn1TABAC(21,4) = 0.0E0_wp
    cn2TABAC(21,4) = -6.7243764463453726E-01_wp
    omega2TAB(21,4) = 1.5180445738268933E0_wp

    !Ca++ <-> NO3-        4P MAR 2006
    bTABAC(21,5) = 0.163281976992298E0_wp 
    cTABAC(21,5) = 0.203681108454362E0_wp  
    cn1TABAC(21,5) = 0.0E0_wp
    cn2TABAC(21,5) = -7.545167774717340E-02_wp
    omega2TAB(21,5) = 1.21090585828713E0_wp
    TABhighestwtf(21,5) = 0.50E0_wp 
    TABKsp(21,5) = 2100.0E0_wp ![molal basis]  

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
    bTABAC(21,8) = 8.2098127298847412E-02_wp
    cTABAC(21,8) = -1.7166655056598561E0_wp
    cn1TABAC(21,8) = 0.0E0_wp
    cn2TABAC(21,8) = 8.6726419889762707E-02_wp
    omega2TAB(21,8) = 4.6756611045939778E-01_wp

    !Ca++ <-> HCO3- !Jan 2021 
    bTABAC(21,10) = 4.0120219199493851E-01_wp
    cTABAC(21,10) = -7.9562337248535087E-01_wp
    cn1TABAC(21,10) = 0.0E0_wp
    cn2TABAC(21,10) = -6.4427030445717304E-01_wp
    omega2TAB(21,10) = 1.0886391158515725E0_wp

    !Ca++ <-> SO4--      4P JAN 2011
    bTABAC(21,21) = 1.2956723534082162E0_wp
    cTABAC(21,21) = -6.9680634914651929E-01_wp
    omegaTAB(21,21) = 0.8E0_wp !usually set to 0.8
    cn1TABAC(21,21) = 0.0E0_wp
    cn2TABAC(21,21) = 1.5915889122259039E0_wp
    omega2TAB(21,21) = 2.5621741189504110E-01_wp
    TABhighestwtf(21,21) = 2.0E-03_wp !CaSO4   

    !Ca++ <-> CO3--  !Jan 2021 
    bTABAC(21,22) = 6.0377736237857693E-01_wp
    cTABAC(21,22) = -8.0427942772938049E-03_wp
    omegaTAB(21,22) = 1.0257268968493378E0_wp
    cn1TABAC(21,22) = 0.0E0_wp
    cn2TABAC(21,22) = -4.9989429639370098E-01_wp
    omega2TAB(21,22) = 1.5894433704988793E0_wp

    !Mg++ <-> Cl-        4P JUN 2007
    bTABAC(23,2) = 0.195908827051115E0_wp
    cTABAC(23,2) = 0.332387167778681E0_wp
    cn1TABAC(23,2) = 0.0E0_wp
    cn2TABAC(23,2) = 7.865437200679259E-02_wp    
    omega2TAB(23,2) = 0.397920452995359E0_wp
    TABhighestwtf(23,2) = 0.32E0_wp    !MgCl2
    TABKsp(23,2) = 1.8986E04_wp  !hydrate-Ksp [molal basis] for saturated solution with respect to MgCl2.6H2O (hydrate); solubility value from Rard: 5.8101 mol MgCl2/kg water

    !Mg++ <-> Br-        4P JUL 2010
    bTABAC(23,3) = 0.260487E0_wp
    cTABAC(23,3) = 1.017037E0_wp
    cn1TABAC(23,3) = 0.0E0_wp
    cn2TABAC(23,3) = 6.162636E-02_wp    
    omega2TAB(23,3) = 0.299475E0_wp
    TABhighestwtf(23,3) = 0.48E0_wp    !MgBr2

    !Mg++ <-> I- 
    bTABAC(23,4) = 9.6129689010970776E-01_wp
    cTABAC(23,4) = 7.7356425470356993E-01_wp   
    cn1TABAC(23,4) = 0.0E0_wp
    cn2TABAC(23,4) = -4.7285705361984376E-01_wp
    omega2TAB(23,4) =  1.0827821509817128E0_wp

    !Mg++ <-> NO3-      4P JUN 2007
    bTABAC(23,5) = 0.430670745401918E0_wp
    cTABAC(23,5) = 0.767241847531231E0_wp  
    cn1TABAC(23,5) = 0.0E0_wp
    cn2TABAC(23,5) = -0.511836366410133E0_wp 
    omega2TAB(23,5) = 1.44093984791712E0_wp 
    TABhighestwtf(23,5) = 0.54E0_wp
    TABKsp(23,5) = 30383.0E0_wp ![molal basis]     

    !Mg++ <-> IO3-     JUN 2021
    bTABAC(23,6) = 1.2126052701130180E-01_wp
    cTABAC(23,6) = 1.7147232434049486E-01_wp
    cn1TABAC(23,6) = 0.0E0_wp
    cn2TABAC(23,6) = -5.3215011337300289E-01_wp
    omega2TAB(23,6) = 8.7729205968138912E-01_wp

    !Mg++ <-> OH-      estimated by Mg++ <-> I-
    bTABAC(23,7) = bTABAC(23,4) 
    cTABAC(23,7) = cTABAC(23,4)
    cn1TABAC(23,7) = cn1TABAC(23,4)
    cn2TABAC(23,7) = cn2TABAC(23,4)
    omega2TAB(23,7) =  omega2TAB(23,4)

    !Mg++ <-> HSO4-     JULY 2019
    bTABAC(23,8) = 1.7141381771763914E-01_wp
    cTABAC(23,8) =  2.7138143502966559E0_wp
    cn1TABAC(23,8) = 0.0E0_wp
    cn2TABAC(23,8) = 2.0404160088725892E-01_wp
    omega2TAB(23,8) = 5.0134675159022479E-01_wp

    !Mg++ <-> HCO3-  !Jan 2021 
    bTABAC(23,10) = 1.8773950264361727E-01_wp
    cTABAC(23,10) = -1.0243549517195567E0_wp
    omegaTAB(23,10) = 1.8303653537783598E-01_wp
    cn1TABAC(23,10) = 0.0E0_wp
    cn2TABAC(23,10) = 2.1693552210594391E-01_wp
    omega2TAB(23,10) = 1.3960864277476452E0_wp

    !Mg++ <-> SO4--      4P JUN 2007
    bTABAC(23,21) = 0.122364180080768E0_wp
    cTABAC(23,21) =  -3.42587571353539E0_wp
    cn1TABAC(23,21) = 0.0E0_wp
    cn2TABAC(23,21) = -0.738560814131576E0_wp   
    omega2TAB(23,21) = 0.864380270721573E0_wp
    TABhighestwtf(23,21) = 0.60E0_wp !MgSO4 
    TABKsp(23,21) = 999999.0E0_wp ![molal basis]   !this value is not correct, just as initial value here. MgSO4 forms a heptahydrate at room temp.: MgSO4.7H2O

    !Mg++ <-> CO3--  !Jan 2021 
    bTABAC(23,22) = 6.8615185959383962E-01_wp
    cTABAC(23,22) = -2.2381706037879368E0_wp
    omegaTAB(23,22) = 1.0578308708725501E-01_wp
    cn1TABAC(23,22) = 0.0E0_wp
    cn2TABAC(23,22) = -5.6539053148654073E-01_wp
    omega2TAB(23,22) = 1.2004054440397884E0_wp

    !H+ <-> CH3SO3-  5P JUN 2017;  [MSA = methanesulfonic acid]
    bTABAC(5,9) = 4.679002E-02_wp
    cTABAC(5,9) = 4.562208E-01_wp
    omegaTAB(5,9) = 0.24E0_wp
    cn1TABAC(5,9) = 0.0E0_wp
    cn2TABAC(5,9) = -1.702953E-01_wp 
    omega2TAB(5,9) = 1.428179E0_wp
    TABhighestwtf(5,9) = 0.791E0_wp !MSA

    !Na+ <-> CH3SO3-  5P JUN 2017;  [CH3SO3- is the methanesulfonate anion]
    bTABAC(2,9) = 6.9073215E-03_wp
    cTABAC(2,9) =  3.166402726E-01_wp
    omegaTAB(2,9) = 3.2567310229E-01_wp
    cn1TABAC(2,9) = 0.0E0_wp
    cn2TABAC(2,9) = -4.3807515687E-01_wp
    omega2TAB(2,9) = 2.109695752E0_wp
    TABhighestwtf(2,9) = 0.64E0_wp

    !NH4+ <-> CH3SO3-  5P JUN 2017;  [CH3SO3- is the methanesulfonate anion]
    bTABAC(4,9) = 2.337520821E-02_wp
    cTABAC(4,9) = 2.23469882E-02_wp
    omegaTAB(4,9) = 2.20000E-01_wp
    cn1TABAC(4,9) = 0.0E0_wp
    cn2TABAC(4,9) = 2.472697276E-02_wp
    omega2TAB(4,9) = 8.265735355E-01_wp
    TABhighestwtf(4,9) = 0.64E0_wp
    
    !---------------------------------------------------------------- 
    !The interactions of three different ions are described with Qcca and Qcaa:
    !There is also an interaction parameter between NH4+ and H+ ions (Rcc) declared here:
    ALLOCATE( qcca1TAB(23,23,40), RccTAB(23,23) )  !qcca1TAB(40,40,40)
    qcca1TAB = 0.0E0_wp
    RccTAB = 0.0E0_wp
    !RaaTAB = 0.0E0_wp

    !Rcc interaction parameter between two cations: 
    !NH4+|H+ (= H+|NH4+)
    RccTAB(4,5) = -1.544864142E-01_wp
    RccTAB(5,4) = -1.544864142E-01_wp

    !Qcca parameters for some three body compositions: 
    !NH4+|H+|HSO4- (= H+|NH4+|HSO4-) 
    qcca1TAB(4,5,8) = 4.483540853E-04_wp
    qcca1TAB(5,4,8) = 4.483540853E-04_wp

    !""""""""""""""""""""" end OF  cation <-->  anion (in water) """"""""""""""""""""""""

    end subroutine MRdata
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
    subroutine GammaCO2()

    use ModAIOMFACvar, only : SMA, SMC, gnmrln, gnsrln, gnlrln
    use ModSystemProp, only : NGI, Ication, Ianion, idCO2

    implicit none

    integer :: I, NC, NA
    real(wp) :: lngammaCO2
    !...............................
    
    lngammaCO2 = 0.0_wp
    do I = 1,NGI
        NC = Ication(I)
        if (SMC(I) > 0.0_wp) then
            lngammaCO2 = lngammaCO2 + 2.0_wp*SMC(I)*lambdaIN(NC)
        endif
        NA = Ianion(I)
        if (SMA(I) > 0.0_wp) then
            lngammaCO2 = lngammaCO2 + 2.0_wp*SMA(I)*lambdaIN(NA)
        endif
    enddo

    gnmrln(idCO2) = lngammaCO2
    gnsrln(idCO2) = 0.0_wp
    gnlrln(idCO2) = 0.0_wp

    end subroutine GammaCO2
    !=========================================================================================================================

end module ModMRpart