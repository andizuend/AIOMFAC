!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module offering different composition scale conversion subroutines.                *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2004 (non-module versions)                                      *
!*   -> latest changes: 2018/05/28                                                      *
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
!*   -  SUBROUTINE MassFrac2IonMolalities                                               *
!*   -  SUBROUTINE MoleFrac2MassFrac                                                    *
!*   -  SUBROUTINE Inputconc_to_wtf                                                     *
!*   -  SUBROUTINE SpecialInputConcConversion                                           *
!*   -  SUBROUTINE MassFrac2MoleFracMolality                                            *
!*   -  SUBROUTINE zSolution2SpeciesMolality                                            *
!*                                                                                      *
!****************************************************************************************   
MODULE ModCompScaleConversion

!Public Variables:
USE ModSystemProp, ONLY : AnNr, CatNr, ElectComps, ElectNues, Mmass, Ncation, nd, nelectrol, NGI, nindcomp, NKNpNGS, nneutral

IMPLICIT NONE
PUBLIC

!================================================================================================================================= 
    CONTAINS
!================================================================================================================================= 
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to calculate the molalities of the cations and anions from a given      *
    !*   input (phase) composition in mass fractions (wtf) of all components.               *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich (2004 - 2009),                                                  *
    !*   Dept. Chem. Engineering, California Institute of Technology (2009 - 2012),         *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2004                                                            *
    !*   -> latest changes: 2018/05/27                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    PURE SUBROUTINE MassFrac2IonMolalities(wtf, SMC, SMA)

    IMPLICIT NONE
    !interface variables:
    REAL(8),DIMENSION(nindcomp),INTENT(IN) :: wtf   !mass fractions of all components
    REAL(8),DIMENSION(NGI),INTENT(OUT) :: SMC, SMA  !cation and anion molalities
    !Local Variables and Parameters:
    INTEGER(4) :: I, J, K, II, ia, ic
    REAL(8) :: an, cn, sumWN
    !................................................................................................

    !calculate molalities of the cations and anions:
    sumWN = SUM(wtf(1:nneutral))
    SMA = 0.0D0
    SMC = 0.0D0
    ia = 0
    ic = 0
    DO K = 1,nelectrol
        ic = ElectComps(K,1) !cation identifier
        ia = ElectComps(K,2) !anion identifier
        cn = REAL(ElectNues(K,1), KIND=8)  !number of cations ic per electrolyte unit 
        an = REAL(ElectNues(K,2), KIND=8)  !number of anions ia per electrolyte unit 
        I = CatNr(ic)
        J = AnNr(ia)
        II = nneutral+K
        SMC(I) = SMC(I) + (wtf(II)/Mmass(II))*(cn/sumWN) !add molality contribution to cation I from electrolyte component K
        SMA(J) = SMA(J) + (wtf(II)/Mmass(II))*(an/sumWN) !add molality contribution to anion J from electrolyte component K
    ENDDO !K

    END SUBROUTINE MassFrac2IonMolalities
!================================================================================================================================= 
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutines to convert from mole fractions (xin) to mass fractions (wout) for      *
    !*   given input components and their molar masses. The applied procedure is designed   *
    !*   to avoid tiny rounding issues by identifying the most abundant component for       *
    !*   mass fraction summation to exactly 1.0D0 within machine precision.                 *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   Dept. Chem. Engineering, California Institute of Technology, 2009 - 2012           *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2012                                                            *
    !*   -> latest changes: 2018/06/11                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    PURE SUBROUTINE MoleFrac2MassFrac(xin, Mmass, wout) 

    IMPLICIT NONE
    !interface variables:
    REAL(8),DIMENSION(:),INTENT(IN) :: xin, Mmass
    REAL(8),DIMENSION(:),INTENT(OUT) :: wout
    !local variables:
    INTEGER(4) :: nc, nmaxindex
    REAL(8) :: totmass, sum1, sum2
    !...................................................
    !the conversion assumes that the mole fractions in xin sum to exactly 1.0D0 and makes sure 
    !that also the sum of the mass fractions equals 1.0D0, i.e., by using this constraint to 
    !avoid potential round-off inaccuracies within machine precision:
    nc = SIZE(xin)
    nmaxindex = MAXLOC(xin, DIM=1)
    totmass = SUM(xin*Mmass)
    IF (nmaxindex > 1) THEN
        wout(1:nmaxindex-1) = xin(1:nmaxindex-1)*Mmass(1:nmaxindex-1)/totmass
        sum1 = SUM(wout(1:nmaxindex-1))
    ELSE
        sum1 = 0.0D0
    ENDIF
    IF (nmaxindex < nc) THEN
        wout(nmaxindex+1:) = xin(nmaxindex+1:)*Mmass(nmaxindex+1:)/totmass
        sum2 = SUM(wout(nmaxindex+1:))
    ELSE
        sum2 = 0.0D0
    ENDIF
    !set the mass fraction value of the most abundant component, which is least sensitive to a tiny rounding error:
    wout(nmaxindex) = MAX(1.0D0-sum1-sum2, 0.0D0) !MAX() to ensure that wout is never a negative value.

    END SUBROUTINE MoleFrac2MassFrac
!================================================================================================================================= 
    

    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   This subroutine converts the input concentration given in mole or mass fraction to *
    !*   mass fractions (wtf).                                                              *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   Div. Chem. Engineering, California Institute of Technology, 2009 - 2012            *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2005                                                            *
    !*   -> latest changes: 2018/05/28                                                      *
    !*                                                                                      *
    !****************************************************************************************  
    SUBROUTINE Inputconc_to_wtf(inputconc, mixingratio, wtfdry, xinput, wtf)

    IMPLICIT NONE
    !interface variables:
    REAL(8),DIMENSION(nindcomp),INTENT(IN) :: inputconc   !the concentration of a given input point (e.g., at an experimental data point)
    REAL(8),DIMENSION(nelectrol),INTENT(IN) :: mixingratio, wtfdry
    LOGICAL(4),INTENT(IN) :: xinput !"true" indicates input is in mole fraction ("false" indicates mass fraction input)
    REAL(8),DIMENSION(nindcomp),INTENT(OUT) :: wtf
    !local variables:
    INTEGER(4), DIMENSION(1) :: minlocwtf, maxlocwtf
    REAL(8),PARAMETER :: lowval = 1.0D2*EPSILON(1.0D0)
    REAL(8),DIMENSION(nindcomp) :: x 
    LOGICAL(4) :: defaultcase
    !...................................
    wtf = 0.0D0
    defaultcase = .true.
    !===
    !consider special cases that need a scaling of input amounts (e.g. to distribute among multiple salts or PEG-oligomer components);
    !this is not needed for general customized input (e.g. remove for AIOMFAC-web version)
    CALL SpecialInputConcConversion(inputconc, mixingratio, wtfdry, xinput, wtf, defaultcase)
    !===
    IF (defaultcase) THEN  !(defaultcase should be set .true. if SpecialInputConcConversion is not used)
        IF (xinput) THEN
            x(2:nindcomp) = inputconc(2:nindcomp)  !mole fraction (with respect to salts not dissociated into ions) of other components including salts!
            x(1) = 1.0D0-SUM(x(2:nindcomp)) !for component water usually
            CALL MoleFrac2MassFrac(x, Mmass, wtf)
        ELSE
            wtf(2:nindcomp) = inputconc(2:nindcomp)
            wtf(1) = 1.0D0-SUM(wtf(2:nindcomp))
        ENDIF
    ENDIF
                    
    !check and correct mixture composition if necessary (avoiding floating point exceptions):
    IF (ANY(wtf(1:nindcomp) < 0.0D0)) THEN
        IF (ABS(MINVAL(wtf(1:nindcomp))) < 1.0D-8) THEN !correct floating point rounding problem
            minlocwtf(1) = MINLOC(wtf(1:nindcomp), DIM=1)
            maxlocwtf(1) = MAXLOC(wtf(1:nindcomp), DIM=1)
            wtf(maxlocwtf(1)) = wtf(maxlocwtf(1))+wtf(minlocwtf(1))
            wtf(minlocwtf(1)) = 0.0D0
        ELSE !there is something wrong...
            WRITE(*,*) ""
            WRITE(*,*) "WARNING from Inputconc_to_wtf: weight fraction of a component is less then 0.0 !!"
            WRITE(*,*) "nd, wtf(1:nindcomp): ", nd, wtf(1:nindcomp)
            WRITE(*,*) ""
            !  READ(*,*)
            RETURN
        ENDIF
    ENDIF
    IF (SUM(wtf(1:nneutral)) < lowval .AND. SUM(wtf(nneutral+1:nindcomp)) > lowval) THEN  !there has to be some water in the mixture or some organic solvent!!
        wtf(2:nindcomp) = wtf(2:nindcomp)*(1.0D0 -lowval)
        wtf(1) = 1.0D0-SUM(wtf(2:nindcomp))
    ENDIF

    END SUBROUTINE Inputconc_to_wtf
!================================================================================================================================= 
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to set molar mixing ratios and dry mass fractions of certain salt or    * 
    !*   PEG mixtures. This is used for model fits to experimental data and for             *
    !*   model output calculation of a list of special mixture systems.                     *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2005                                                            *
    !*   -> latest changes: 2018/05/28                                                      *
    !*                                                                                      *
    !****************************************************************************************   
    PURE SUBROUTINE SpecialInputConcConversion(inputconc, mixingratio, wtfdry, xinput, wtf, defaultcase)
    
    IMPLICIT NONE
    !interface:
    REAL(8),DIMENSION(nindcomp),INTENT(IN) :: inputconc   !inputconc = the concentration of a given input point (e.g., at an experimental data point)
    REAL(8),DIMENSION(nelectrol),INTENT(IN) :: mixingratio, wtfdry
    LOGICAL(4),INTENT(IN) :: xinput
    REAL(8),DIMENSION(nindcomp),INTENT(INOUT) :: wtf
    LOGICAL(4),INTENT(OUT) :: defaultcase
    !local:
    REAL(8) :: totalweight, SUMSalts
    REAL(8),DIMENSION(nindcomp) :: x 
    !..........................................

    defaultcase = .false. !initialize

    !weightfractions of the read in data:
    IF (.NOT. xinput) THEN !data is read in in mass fraction scale
        SELECT CASE(nd)
        CASE(26)
            wtf(2) = inputconc(2) !inputconc(2)
            wtf(1) = inputconc(3) !inputconc(3)
            wtf(3) = 0.0D0
        CASE(111)
            wtf(2) = inputconc(2) !inputconc(2)
            wtf(3) = inputconc(3) !inputconc(3)
            wtf(4) = inputconc(4) !inputconc(4) !second salt
            wtf(1) = 1.0D0-SUM(inputconc(2:4))  !water wtf
        CASE(117,122,123,134:135,199) !saltmixes
            wtf(nneutral+1:nindcomp) = inputconc(2)*wtfdry(1:nelectrol)
            wtf(1) = 1.0D0-inputconc(2)
        CASE(38,133,139:145,180:181,191,200:202) !saltmixes
            SUMSalts = SUM(Mmass(nneutral+1:nindcomp)*mixingratio(1:nelectrol))
            wtf(nneutral+1:nindcomp) = inputconc(2)*Mmass(nneutral+1:nindcomp)*mixingratio(1:nelectrol)/SUMSalts
            wtf(1) = 1.0D0-inputconc(2)
        CASE(203:204) !PEG-400 weight ratios for a 1:2 molar mixing ratio:
            wtf(2) = inputconc(2)*0.30885D0 !PEG-400-n7
            wtf(3) = inputconc(2)*0.69115D0 !PEG-400-n8
            wtf(1) = 1.0D0-inputconc(2) !water wtf
        CASE(205:209,219:221) !PEG-400 weight ratios for a 1:2 molar mixing ratio:
            wtf(2) = inputconc(2)*0.30885D0 !PEG-400-n7
            wtf(3) = inputconc(2)*0.69115D0 !PEG-400-n8
            wtf(4) = inputconc(3)
            wtf(1) = 1.0D0-inputconc(2)-inputconc(3) !water wtf
        !CASE(210,217,218) !PEG-1000 weight ratios for a 0.708913:0.291087 molar mixing ratio:
        !    wtf(2) = inputconc(2)*0.69982246682D0 !PEG-1000-n21
        !    wtf(3) = inputconc(2)*0.30017753058D0 !PEG-1000-n22
        !    wtf(4) = inputconc(3)
        !    wtf(1) = 1.0D0-inputconc(2)-inputconc(3) !water wtf
        CASE(211) !PEG-1000 weight ratios for a 0.708913:0.291087 molar mixing ratio:
            wtf(2) = inputconc(2)*0.69982246682D0 !PEG-1000-n21
            wtf(3) = inputconc(2)*0.30017753058D0 !PEG-1000-n22
            wtf(1) = 1.0D0-inputconc(2) !water wtf
        CASE(212) !PEG-600 weight ratios for a 0.708913:0.291087 molar mixing ratio:
            wtf(2) = inputconc(2)*0.7767D0 !PEG-1000-n12
            wtf(3) = inputconc(2)*0.2233D0 !PEG-1000-n13
            wtf(1) = 1.0D0-inputconc(2) !water wtf
        CASE(213) !PEG-1450 weight ratios for a 0.708913:0.291087 molar mixing ratio:
            wtf(2) = inputconc(2)*0.48630405D0 !PEG-1000-n31
            wtf(3) = inputconc(2)*0.51369595D0 !PEG-1000-n32
            wtf(1) = 1.0D0-inputconc(2) !water wtf
        !CASE(214) !PEG-1540 weight ratios for a 0.45089529:0.54910471 molar mixing ratio:
        !    wtf(2) = inputconc(2)*0.44381284D0 !PEG-1540-n33
        !    wtf(3) = inputconc(2)*0.55618716D0 !PEG-1540-n34
        !    wtf(4) = inputconc(3)
        !    wtf(1) = 1.0D0-inputconc(2)-inputconc(3) !water wtf
        CASE(1157) !PEG-200 as a mix of two PEG chain length; mass fractions are: 0.843854566*(PEG-200-n3) + 0.156145434*(PEG-200-n4)
            wtf(2) = inputconc(2)*0.843854566D0 !PEG-200-n3
            wtf(3) = inputconc(2)*0.156145434D0 !PEG-200-n4
            wtf(1) = 1.0D0-inputconc(2) !water wtf
        CASE DEFAULT
            defaultcase = .true.
        END SELECT
    ELSE !data is read in in mole fractions x => conversion to wtf is necessary:
        SELECT CASE(nd)
        CASE(38,133,139:145,180:181,191,200:202) !saltmixes
            x(nneutral+1:nindcomp) = inputconc(2)*mixingratio(1:nelectrol)
            x(1) = 1.0D0-SUM(x(nneutral+1:nindcomp))
            totalweight = SUM(Mmass(1:nneutral)*x(1:nneutral))
            totalweight = totalweight+SUM(Mmass(nneutral+1:nindcomp)*x(nneutral+1:nindcomp))
            wtf(nneutral+1:nindcomp) = x(nneutral+1:nindcomp)*Mmass(nneutral+1:nindcomp)/totalweight
            wtf(1) = 1.0D0-SUM(wtf(2:nindcomp))
        CASE(196:198)
            x(1) = 1.0D0-inputconc(2)-inputconc(3)
            x(nneutral+1:nindcomp) = inputconc(2)*mixingratio(1:nelectrol)/SUM(mixingratio(1:nelectrol))
            totalweight = SUM(Mmass(1:nneutral)*x(1:nneutral))
            totalweight = totalweight+SUM(Mmass(nneutral+1:nindcomp)*x(nneutral+1:nindcomp))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
            wtf(nneutral+1:nindcomp) = x(nneutral+1:nindcomp)*Mmass(nneutral+1:nindcomp)/totalweight
        CASE(203:204) !PEG-400 molar ratios
            x(1) = 1.0D0-inputconc(2)
            x(2) = inputconc(2)*(1.0D0/3.0D0) !PEG-400 n = 7
            x(3) = inputconc(2)*(2.0D0/3.0D0) !PEG-400 n = 8
            totalweight = SUM(Mmass(1:nneutral)*x(1:nneutral))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        CASE(205:209,219:221) !PEG-400 molar ratios
            x(1) = 1.0D0-inputconc(2)-inputconc(3)
            x(2) = inputconc(2)*(1.0D0/3.0D0) !PEG-400 n = 7
            x(3) = inputconc(2)*(2.0D0/3.0D0) !PEG-400 n = 8
            x(4) = inputconc(3) !AS
            totalweight = SUM(Mmass(1:nneutral)*x(1:nneutral))
            totalweight = totalweight+SUM(Mmass(nneutral+1:nindcomp)*x(nneutral+1:nindcomp))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
            wtf(nneutral+1:nindcomp) = x(nneutral+1:nindcomp)*Mmass(nneutral+1:nindcomp)/totalweight
        !CASE(210,217,218) !PEG-1000 molar ratios
        !    x(1) = 1.0D0-inputconc(2)-inputconc(3)
        !    x(2) = inputconc(2)*0.708913D0 !PEG-1000 n = 21
        !    x(3) = inputconc(2)*0.291087D0 !PEG-1000 n = 22
        !    x(4) = inputconc(3) !AS
        !    totalweight = SUM(Mmass(1:nneutral)*x(1:nneutral))
        !    totalweight = totalweight+SUMMmass(nneutral+1:nindcomp)*x(nneutral+1:nindcomp))
        !    wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        !    wtf(nneutral+1:nindcomp) = x(nneutral+1:nindcomp)*Mmass(nneutral+1:nindcomp)/totalweight
        CASE(211) !PEG-1000  0.708913:0.291087 molar mixing ratio:
            x(1) = 1.0D0-inputconc(2)-inputconc(3)
            x(2) = inputconc(2)*0.708913D0
            x(3) = inputconc(2)*0.291087D0
            totalweight = SUM(Mmass(1:nneutral)*x(1:nneutral))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        CASE(212) !PEG-600 as a mixture of 0.788926 PEG-600-n12 and 0.211074 PEG-600-n13
            x(1) = 1.0D0-inputconc(2)-inputconc(3)
            x(2) = inputconc(2)*0.788926D0
            x(3) = inputconc(2)*0.211074D0
            totalweight = SUM(Mmass(1:nneutral)*x(1:nneutral))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        CASE(213) !PEG-1450 (PEG-1450 as a mixture of 0.493898 PEG-1450-n31 and 0.506102 PEG-1450-n32)
            x(1) = 1.0D0-inputconc(2)-inputconc(3)
            x(2) = inputconc(2)*0.493898D0
            x(3) = inputconc(2)*0.506102D0
            totalweight = SUM(Mmass(1:nneutral)*x(1:nneutral))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        !CASE(214) !PEG-1540 molar ratios
        !    x(1) = 1.0D0-inputconc(2)-inputconc(3)
        !    x(2) = inputconc(2)*0.45089529 !PEG-1540 n = 33
        !    x(3) = inputconc(2)*0.54910471 !PEG-1540 n = 34
        !    x(4) = inputconc(3) !AS
        !    totalweight = SUM(Mmass(1:nneutral)*x(1:nneutral))
        !    totalweight = totalweight+SUM(Mmass(nneutral+1:nindcomp)*x(nneutral+1:nindcomp))
        !    wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        !    wtf(nneutral+1:nindcomp) = x(nneutral+1:nindcomp)*Mmass(nneutral+1:nindcomp)/totalweight
        CASE(510)  !ternary neutral mixture (no salts)
            x(2:nindcomp) = inputconc(2:nindcomp)  !mole fraction (with respect to salts not dissociated into ions) of other components including salts!
            x(1) = 1.0D0-SUM(x(2:nindcomp)) !for component water usually
            totalweight = SUM(Mmass(1:nneutral)*x(1:nneutral))
            wtf(2:nneutral) = x(2:nneutral)*Mmass(2:nneutral)/totalweight
            wtf(1) = 1.0D0-SUM(wtf(2:nindcomp))
        CASE(1157) !PEG-200; mole fractions of PEG-200-n3 and PEG-200-n4 to yield 200 g/mol:  0.8689391*(PEG-200-n3) + 0.1310609*(PEG-200-n4);
            x(1) = 1.0D0-inputconc(2)
            x(2) = inputconc(2)*0.8689391D0 !PEG-200 n = 3
            x(3) = inputconc(2)*0.1310609D0 !PEG-200 n = 4
            totalweight = SUM(Mmass(1:nneutral)*x(1:nneutral))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        CASE DEFAULT
            defaultcase = .true.
        END SELECT
    ENDIF

    END SUBROUTINE SpecialInputConcConversion
!=================================================================================================================================
    
    
    !********************************************************************************
    !*                        'MassFrac2MoleFracMolality'                           *
    !*                                                                              *
    !*  This routine calculates the mole fraction and the molality of the           *
    !*  components, including of the salts, from the given mass fractions (wtf).    *
    !*  In this case, mole fractions and molalities are calculated with respect to  *
    !*  undissociated (!) electrolytes of given ion pairings at input.              *
    !*                                                                              *
    !*            (c) Andi Zuend, IACETH, ETH Zurich, 2007 - 2009;                  *
    !*   Dept. Chem. Engineering, California Institute of Technology, 2009 - 2012   *
    !*                                                                              *
    !********************************************************************************
    PURE SUBROUTINE MassFrac2MoleFracMolality(wtf, XrespSalt, mrespSalt)

    IMPLICIT NONE
    !interface variables:
    REAL(8),DIMENSION(:),INTENT(IN) :: wtf          !mass fractions (input)
    REAL(8),DIMENSION(:),INTENT(OUT) :: XrespSalt   !mole fraction of the components with respect to undissociated salts/electrolytes (not individual ions)
    REAL(8),DIMENSION(:),INTENT(OUT) :: mrespSalt   !molality of the components with respect to undissociated electrolytes
    !local variables:
    INTEGER(4) :: i, nnp1
    REAL(8) :: totalmolenumber, sumsaltfreeWTF
    REAL(8),DIMENSION(SIZE(wtf)) :: saltfreeWTF, wtfbyMmass
    !................................................................

    nnp1 = nneutral+1
    XrespSalt = 0.0D0
    mrespSalt = 0.0D0
    wtfbyMmass = wtf/Mmass
    totalmolenumber = SUM(wtfbyMmass)

    saltfreeWTF = 0.0D0
    sumsaltfreeWTF = SUM(wtf(1:nneutral))
    IF (sumsaltfreeWTF > 0.0D0) THEN
        saltfreeWTF(1:nneutral) = wtf(1:nneutral)/sumsaltfreeWTF
    ENDIF
    !total number of moles of substances is known, now one can calculate the mole fraction:
    !for the neutrals:
    XrespSalt(1:nneutral) = wtfbyMmass(1:nneutral)/totalmolenumber
    mrespSalt(1:nneutral) = saltfreeWTF(1:nneutral)/Mmass(1:nneutral)
    !for the electrolytes
    IF (nelectrol > 0) THEN
        XrespSalt(nnp1:) = wtfbyMmass(nnp1:)/totalmolenumber
        IF (wtf(1) > 0.0D0) THEN
            mrespSalt(nnp1:) = wtfbyMmass(nnp1:)*(saltfreeWTF(1)/wtf(1))
        ELSE
            IF (nneutral > 1 .AND. ANY(wtf(2:nneutral) > 0.0D0)) THEN
                DO i = 2,nneutral
                    IF (wtf(i) > 0.0D0) THEN
                        EXIT !save the ith neutral
                    ENDIF
                ENDDO
                mrespSalt(nnp1:nindcomp) = wtfbyMmass(nnp1:nindcomp)*(saltfreewtf(i)/wtf(i))
            ELSE
                i = 777
                !!$OMP CRITICAL (MMA1)
                !WRITE(*,*) "WARNING from MassFrac2MoleFracMolality: wtf(1) = 0.0 "
                !WRITE(*,*) "There has to be some water in the mixture when inorganic salts are present!"
                !WRITE(*,*) "wtf(1:nindcomp): ", wtf(1:nindcomp)
                !WRITE(*,*) ""
                !READ(*,*)
                !!$OMP END CRITICAL (MMA1)
            ENDIF
        ENDIF
    ENDIF

    END SUBROUTINE MassFrac2MoleFracMolality
!================================================================================================================================= 
    
      
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to compute solvent and dissociated ion molalities from input of mixture * 
    !*   composition in mole fractions with respect to undissociated electrolytes (zl)      *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2005                                                            *
    !*   -> latest changes: 2018/05/28                                                      *
    !*                                                                                      *
    !****************************************************************************************  
    PURE SUBROUTINE zSolution2SpeciesMolality(zl, ml) 
    
    IMPLICIT NONE
    !interface variables:
    REAL(8),DIMENSION(nindcomp),INTENT(IN) :: zl !input mole fraction (undissociated electrolytes)
    REAL(8),DIMENSION(NKNpNGS),INTENT(OUT) :: ml !output molalities of neutral solvent components and dissociated ions (first cations then anions according to order in Ication, Ianion)
    !local variables:
    INTEGER(4) :: i, cn, an, cid, aid
    REAL(8) :: Msolv
    !................................
    ml = 0.0D0
    !liquid solvent mass:
    Msolv = SUM(zl(1:nneutral)*Mmass(1:nneutral))
    !molality of solvent components:
    ml(1:nneutral) = zl(1:nneutral)/Msolv
    !molality of individual ions:
    DO i = 1,nelectrol  
        cn = ElectComps(i,1) !the cation of this electrolyte
        an = ElectComps(i,2) !the anion
        cid = CatNr(cn) !the cation index ID within the cations of this mixture
        aid = AnNr(an)  !the anion index ID
        ml(nneutral+cid) = ml(nneutral+cid) + zl(nneutral+i)*REAL(ElectNues(i,1), KIND=8)/Msolv
        ml(nneutral+Ncation+aid) = ml(nneutral+Ncation+aid) + zl(nneutral+i)*REAL(ElectNues(i,2), KIND=8)/Msolv
    ENDDO
    
    END SUBROUTINE zSolution2SpeciesMolality
!==========================================================================================================================

END MODULE ModCompScaleConversion