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
!*   -> latest changes: 2021-10-01                                                      *
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
!*   -  subroutine MassFrac2IonMolalities                                               *
!*   -  subroutine MassFrac2SolvMolalities                                              *
!*   -  subroutine MoleFrac2MassFrac                                                    *
!*   -  subroutine Moles2solvmass                                                       *
!*   -  subroutine Molality2SolvMoleFrac                                                *
!*   -  subroutine Inputconc_to_wtf                                                     *
!*   -  subroutine SpecialInputConcConversion                                           *
!*   -  subroutine MassFrac2MoleFracMolality                                            *
!*   -  subroutine zSolution2SpeciesMolality                                            *
!*                                                                                      *
!****************************************************************************************   
module ModCompScaleConversion

!Public Variables:
use Mod_kind_param, only : wp
use ModSystemProp, only : AnNr, CatNr, ElectComps, ElectNues, Mmass, Ncation, nd, nelectrol, &
    & NGI, nindcomp, NKNpNGS, nneutral

implicit none
public

!========================================================================================================== 
    contains
!========================================================================================================== 
    
    
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
    !*   -> latest changes: 2018-05-27                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    pure subroutine MassFrac2IonMolalities(wtf, SMC, SMA)

    implicit none
    !interface variables:
    real(wp),dimension(nindcomp),intent(in) :: wtf       !mass fractions of all components
    real(wp),dimension(NGI),intent(out) :: SMC, SMA      !cation and anion molalities
    !Local Variables and Parameters:
    integer :: I, J, K, II, ia, ic
    real(wp) :: an, cn, sumWN
    !................................................................................................

    !calculate molalities of the cations and anions:
    sumWN = sum(wtf(1:nneutral))
    SMA = 0.0_wp
    SMC = 0.0_wp
    ia = 0
    ic = 0
    do K = 1,nelectrol
        ic = ElectComps(K,1)                !cation identifier
        ia = ElectComps(K,2)                !anion identifier
        cn = real(ElectNues(K,1), kind=wp)   !number of cations ic per electrolyte unit 
        an = real(ElectNues(K,2), kind=wp)   !number of anions ia per electrolyte unit 
        I = CatNr(ic)
        J = AnNr(ia)
        II = nneutral+K
        SMC(I) = SMC(I) + (wtf(II)/Mmass(II))*(cn/sumWN)    !add molality contribution to cation I from electrolyte component K
        SMA(J) = SMA(J) + (wtf(II)/Mmass(II))*(an/sumWN)    !add molality contribution to anion J from electrolyte component K
    enddo !K

    end subroutine MassFrac2IonMolalities
!========================================================================================================== 
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to calculate the molalities of the neutral (solvent) species from a     *
    !*   given input (phase) composition in mass fractions (wtf) of all components.         *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Hang Yin and Andi Zuend,                                                           *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2020                                                            *
    !*   -> latest changes: 2023-03-17                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    pure subroutine MassFrac2SolvMolalities(wtf, m_neutral)

    implicit none
    !interface variables:
    real(wp),dimension(nindcomp),intent(in) :: wtf           ![-]      mass fractions of all components
    real(wp),dimension(nneutral),intent(out) :: m_neutral    ![mol/kg] neutral component molalities
    !..........................................
    
    m_neutral = wtf(1:nneutral) / ( sum(wtf(1:nneutral)) * Mmass(1:nneutral) )

    end subroutine MassFrac2SolvMolalities
!==========================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *																							 
    !*   Subroutines to convert from mole fractions (xin) to mass fractions (wout) for      *
    !*   given input components and their molar masses. The applied procedure is designed   *
    !*   to avoid tiny rounding issues by identifying the most abundant component for       *
    !*   mass fraction summation to exactly 1.0_wp within machine precision.                 *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   Dept. Chem. Engineering, California Institute of Technology, 2009 - 2012           *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2012                                                            *
    !*   -> latest changes: 2019-10-29                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    pure subroutine MoleFrac2MassFrac(xin, Mmass, wout) 

    implicit none
    !interface variables:
    real(wp),dimension(:),intent(in) :: xin, Mmass
    real(wp),dimension(:),intent(out) :: wout
    !local variables:
    integer :: nc, nmaxindex
    real(wp) :: totmass, sum1, sum2
    !...................................................
    !the conversion assumes that the mole fractions in xin sum to exactly 1.0_wp and makes sure 
    !that also the sum of the mass fractions equals 1.0_wp, i.e., by using this constraint to 
    !avoid potential round-off inaccuracies within machine precision:
    nc = size(xin)
    nmaxindex = maxloc(xin, dim=1)
    totmass = sum(xin*Mmass)
    if (totmass > 0.0_wp) then
        if (nmaxindex > 1) then
            wout(1:nmaxindex-1) = xin(1:nmaxindex-1)*Mmass(1:nmaxindex-1)/totmass
            sum1 = sum(wout(1:nmaxindex-1))
        else
            sum1 = 0.0_wp
        endif
        if (nmaxindex < nc) then
            wout(nmaxindex+1:) = xin(nmaxindex+1:)*Mmass(nmaxindex+1:)/totmass
            sum2 = sum(wout(nmaxindex+1:))
        else
            sum2 = 0.0_wp
        endif
    else !zero total mass; assign artificial mass fraction distribution
        wout = 1.0_wp/real(nc, kind=wp)
        sum1 = sum(wout(1:nc-1))
        sum2 = 0.0_wp
    endif
    !set the mass fraction value of the most abundant component, which is least sensitive to a tiny rounding error:
    wout(nmaxindex) = max(1.0_wp-sum1-sum2, 0.0_wp) !max() to ensure that wout is never a negative value.

    end subroutine MoleFrac2MassFrac
!==========================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Utility subroutines to update solvent mass.                                        *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Hang Yin and Andi Zuend,                                                           *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2020                                                            *
    !*   -> latest changes: 2020-06-06                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    pure subroutine Moles2solvmass(moleNeutral, solvmass)

    implicit none
    !interface variables:
    real(wp),dimension(nneutral),intent(in) :: moleNeutral
    real(wp),intent(out) :: solvmass                         ![kg]
    !...................................................
    
    solvmass = sum(moleNeutral(1:nneutral)*Mmass(1:nneutral))

    end subroutine Moles2solvmass
!==========================================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to calculate the mole fractions of the neutral species from a given     *
    !*   input (phase) composition in molality of all components.                           *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Hang Yin and Andi Zuend,                                                           *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2020                                                            *
    !*   -> latest changes: 2023-03-17                                                      *
    !*                                                                                      *
    !**************************************************************************************** 
    pure subroutine Molality2SolvMoleFrac(SMA, SMC, m_neutral, xout)

    implicit none
    !interface variables:
    real(wp),dimension(:),intent(in) :: m_neutral        ![mol/kg] neutral molalities
    real(wp),dimension(:),intent(in) :: SMA, SMC         ![mol/kg] ion molalities
    real(wp),dimension(nneutral),intent(out) :: xout     ![-] neutral component mole fractions (on the basis of dissociated electrolytes) 
    !Local Variables and Parameters:
    real(wp) :: SumIonMolalities, sum_molal
    !.........................................
    
    SumIonMolalities = sum(SMA(1:NGI)) +sum(SMC(1:NGI))
    sum_molal = sum(m_neutral) +SumIonMolalities        !sum of all molalities
    xout = m_neutral / sum_molal                        !mole fraction of the neutral components 
                                                        !(on the basis of partially/fully dissociated electrolytes)    
    end subroutine Molality2SolvMoleFrac
!========================================================================================================== 
    

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
    !*   -> latest changes: 2018-05-28                                                      *
    !*                                                                                      *
    !****************************************************************************************  
    subroutine Inputconc_to_wtf(inputconc, mixingratio, wtfdry, xinput, wtf)

    implicit none
    !interface variables:
    real(wp),dimension(nindcomp),intent(in) :: inputconc             !the concentration of a given input point (e.g., at an experimental data point)
    real(wp),dimension(nelectrol),intent(in) :: mixingratio, wtfdry
    logical,intent(in) :: xinput                                 !"true" indicates input is in mole fraction ("false" indicates mass fraction input)
    real(wp),dimension(nindcomp),intent(out) :: wtf
    !local variables:
    integer :: minlocwtf, maxlocwtf
    real(wp),parameter :: lowval = 1.0E2_wp*epsilon(1.0_wp)
    real(wp),dimension(nindcomp) :: x 
    logical :: defaultcase
    !...................................
    wtf = 0.0_wp
    defaultcase = .true.
    !===
    !consider special cases that need a scaling of input amounts (e.g. to distribute among multiple salts or PEG-oligomer components);
    !this is not needed for general customized input (e.g. remove for AIOMFAC-web version)
    call SpecialInputConcConversion(inputconc, mixingratio, wtfdry, xinput, wtf, defaultcase)
    !===
    if (defaultcase) then                           !(defaultcase should be set .true. if SpecialInputConcConversion is not used)
        if (xinput) then
            x(2:nindcomp) = inputconc(2:nindcomp)   !mole fraction (with respect to salts not dissociated into ions) of other components including salts!
            x(1) = 1.0_wp-sum(x(2:nindcomp))         !for component water usually
            call MoleFrac2MassFrac(x, Mmass, wtf)
        else
            wtf(2:nindcomp) = inputconc(2:nindcomp)
            wtf(1) = 1.0_wp-sum(wtf(2:nindcomp))
        endif
    endif
                    
    !check and correct mixture composition if necessary (avoiding floating point exceptions):
    if (any(wtf(1:nindcomp) < 0.0_wp)) then
        if (abs(minval(wtf(1:nindcomp))) < 1.0E-8_wp) then     !correct floating point rounding problem
            minlocwtf = minloc(wtf(1:nindcomp), dim=1)
            maxlocwtf = maxloc(wtf(1:nindcomp), dim=1)
            wtf(maxlocwtf) = wtf(maxlocwtf)+wtf(minlocwtf)
            wtf(minlocwtf) = 0.0_wp
        else !there is something wrong...
            minlocwtf = minloc(wtf(1:nindcomp), dim=1)
            write(*,*) ""
            write(*,*) "WARNING from Inputconc_to_wtf: mass fraction of a component is less then 0.0 !!"
            write(*,*) "nd, wtf(minlocwtf): ", nd, wtf(minlocwtf)
            write(*,*) ""
            !  read(*,*)
            return
        endif
    endif
    if (sum(wtf(1:nneutral)) < lowval .and. sum(wtf(nneutral+1:nindcomp)) > lowval) then  !there has to be some water in the mixture or some organic solvent!!
        wtf(2:nindcomp) = wtf(2:nindcomp)*(1.0_wp -lowval)
        wtf(1) = 1.0_wp - sum(wtf(2:nindcomp))
    endif

    end subroutine Inputconc_to_wtf
!========================================================================================================== 
    
    
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
    !*   -> latest changes: 2018-05-28                                                      *
    !*                                                                                      *
    !****************************************************************************************   
    pure subroutine SpecialInputConcConversion(inputconc, mixingratio, wtfdry, xinput, wtf, defaultcase)
    
    implicit none
    !interface:
    real(wp),dimension(nindcomp),intent(in) :: inputconc   !inputconc = the concentration of a given input point (e.g., at an experimental data point)
    real(wp),dimension(nelectrol),intent(in) :: mixingratio, wtfdry
    logical,intent(in) :: xinput
    real(wp),dimension(nindcomp),intent(inout) :: wtf
    logical,intent(out) :: defaultcase
    !local:
    real(wp) :: totalweight, SUMSalts
    real(wp),dimension(nindcomp) :: x 
    !..........................................

    defaultcase = .false. !initialize

    !weightfractions of the read in data:
    if (.not. xinput) then !data is read in in mass fraction scale
        select case(nd)
        case(26)
            wtf(2) = inputconc(2) !inputconc(2)
            wtf(1) = inputconc(3) !inputconc(3)
            wtf(3) = 0.0_wp
        case(111)
            wtf(2) = inputconc(2) !inputconc(2)
            wtf(3) = inputconc(3) !inputconc(3)
            wtf(4) = inputconc(4) !inputconc(4) !second salt
            wtf(1) = 1.0_wp-sum(inputconc(2:4))  !water wtf
        case(117,122,123,134:135,199) !saltmixes
            wtf(nneutral+1:nindcomp) = inputconc(2)*wtfdry(1:nelectrol)
            wtf(1) = 1.0_wp-inputconc(2)
        case(38,133,139:145,180:181,191,200:202,1935,1966:1969,1978:1980,1988,3034:3036) !saltmixes
            SUMSalts = sum(Mmass(nneutral+1:nindcomp)*mixingratio(1:nelectrol))
            wtf(nneutral+1:nindcomp) = inputconc(2)*Mmass(nneutral+1:nindcomp)*mixingratio(1:nelectrol)/SUMSalts
            wtf(1) = 1.0_wp-inputconc(2)
        case(203:204) !PEG-400 weight ratios for a 1:2 molar mixing ratio:
            wtf(2) = inputconc(2)*0.30885_wp !PEG-400-n7
            wtf(3) = inputconc(2)*0.69115_wp !PEG-400-n8
            wtf(1) = 1.0_wp-inputconc(2) !water wtf
        case(205:209,219:221) !PEG-400 weight ratios for a 1:2 molar mixing ratio:
            wtf(2) = inputconc(2)*0.30885_wp !PEG-400-n7
            wtf(3) = inputconc(2)*0.69115_wp !PEG-400-n8
            wtf(4) = inputconc(3)
            wtf(1) = 1.0_wp-inputconc(2)-inputconc(3) !water wtf
        !case(210,217,218) !PEG-1000 weight ratios for a 0.708913:0.291087 molar mixing ratio:
        !    wtf(2) = inputconc(2)*0.69982246682_wp !PEG-1000-n21
        !    wtf(3) = inputconc(2)*0.30017753058_wp !PEG-1000-n22
        !    wtf(4) = inputconc(3)
        !    wtf(1) = 1.0_wp-inputconc(2)-inputconc(3) !water wtf
        case(211) !PEG-1000 weight ratios for a 0.708913:0.291087 molar mixing ratio:
            wtf(2) = inputconc(2)*0.69982246682_wp !PEG-1000-n21
            wtf(3) = inputconc(2)*0.30017753058_wp !PEG-1000-n22
            wtf(1) = 1.0_wp-inputconc(2) !water wtf
        case(212) !PEG-600 weight ratios for a 0.708913:0.291087 molar mixing ratio:
            wtf(2) = inputconc(2)*0.7767_wp !PEG-1000-n12
            wtf(3) = inputconc(2)*0.2233_wp !PEG-1000-n13
            wtf(1) = 1.0_wp-inputconc(2) !water wtf
        case(213) !PEG-1450 weight ratios for a 0.708913:0.291087 molar mixing ratio:
            wtf(2) = inputconc(2)*0.48630405_wp !PEG-1000-n31
            wtf(3) = inputconc(2)*0.51369595_wp !PEG-1000-n32
            wtf(1) = 1.0_wp-inputconc(2) !water wtf
        !case(214) !PEG-1540 weight ratios for a 0.45089529:0.54910471 molar mixing ratio:
        !    wtf(2) = inputconc(2)*0.44381284_wp !PEG-1540-n33
        !    wtf(3) = inputconc(2)*0.55618716_wp !PEG-1540-n34
        !    wtf(4) = inputconc(3)
        !    wtf(1) = 1.0_wp-inputconc(2)-inputconc(3) !water wtf
        case(1157) !PEG-200 as a mix of two PEG chain length; mass fractions are: 0.843854566*(PEG-200-n3) + 0.156145434*(PEG-200-n4)
            wtf(2) = inputconc(2)*0.843854566_wp !PEG-200-n3
            wtf(3) = inputconc(2)*0.156145434_wp !PEG-200-n4
            wtf(1) = 1.0_wp-inputconc(2) !water wtf
        case(1958:1961,1965,1970:1975,1985,1989,1996,3004,3022,3023:3030,3051,3053,3054,3056:3060)  !shift 2 to 3 for carbsyst with noCO2(aq)input
            SUMSalts = sum(Mmass(nneutral+1:nindcomp)*mixingratio(1:nelectrol))
            wtf(nneutral+1:nindcomp) = inputconc(2)*Mmass(nneutral+1:nindcomp)*mixingratio(1:nelectrol)/SUMSalts
            wtf(2) = 0.0_wp
            wtf(1) = 1.0_wp-inputconc(2)
        case(3050)  !shift 2 to 3 for malosyst with H2Malo input
            wtf(nneutral+1:nindcomp) = inputconc(2:nindcomp-1)
            wtf(2) = 0.0_wp
            wtf(1) = 1.0_wp-sum( wtf(nneutral+1:nindcomp))
        case default
            defaultcase = .true.
        end select
    else !data is read in in mole fractions x => conversion to wtf is necessary:
        select case(nd)
        case(38,133,139:145,180:181,191,200:202,1935,1966:1969,1978:1980,1988,3034:3036) !saltmixes
            x(nneutral+1:nindcomp) = inputconc(2)*mixingratio(1:nelectrol)
            x(1) = 1.0_wp-sum(x(nneutral+1:nindcomp))
            totalweight = sum(Mmass(1:nneutral)*x(1:nneutral))
            totalweight = totalweight+sum(Mmass(nneutral+1:nindcomp)*x(nneutral+1:nindcomp))
            wtf(nneutral+1:nindcomp) = x(nneutral+1:nindcomp)*Mmass(nneutral+1:nindcomp)/totalweight
            wtf(1) = 1.0_wp-sum(wtf(2:nindcomp))
        case(196:198)
            x(1) = 1.0_wp-inputconc(2)-inputconc(3)
            x(nneutral+1:nindcomp) = inputconc(2)*mixingratio(1:nelectrol)/sum(mixingratio(1:nelectrol))
            totalweight = sum(Mmass(1:nneutral)*x(1:nneutral))
            totalweight = totalweight+sum(Mmass(nneutral+1:nindcomp)*x(nneutral+1:nindcomp))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
            wtf(nneutral+1:nindcomp) = x(nneutral+1:nindcomp)*Mmass(nneutral+1:nindcomp)/totalweight
        case(203:204) !PEG-400 molar ratios
            x(1) = 1.0_wp-inputconc(2)
            x(2) = inputconc(2)*(1.0_wp/3.0_wp) !PEG-400 n = 7
            x(3) = inputconc(2)*(2.0_wp/3.0_wp) !PEG-400 n = 8
            totalweight = sum(Mmass(1:nneutral)*x(1:nneutral))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        case(205:209,219:221) !PEG-400 molar ratios
            x(1) = 1.0_wp-inputconc(2)-inputconc(3)
            x(2) = inputconc(2)*(1.0_wp/3.0_wp) !PEG-400 n = 7
            x(3) = inputconc(2)*(2.0_wp/3.0_wp) !PEG-400 n = 8
            x(4) = inputconc(3) !AS
            totalweight = sum(Mmass(1:nneutral)*x(1:nneutral))
            totalweight = totalweight+sum(Mmass(nneutral+1:nindcomp)*x(nneutral+1:nindcomp))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
            wtf(nneutral+1:nindcomp) = x(nneutral+1:nindcomp)*Mmass(nneutral+1:nindcomp)/totalweight
        !case(210,217,218) !PEG-1000 molar ratios
        !    x(1) = 1.0_wp-inputconc(2)-inputconc(3)
        !    x(2) = inputconc(2)*0.708913_wp !PEG-1000 n = 21
        !    x(3) = inputconc(2)*0.291087_wp !PEG-1000 n = 22
        !    x(4) = inputconc(3) !AS
        !    totalweight = sum(Mmass(1:nneutral)*x(1:nneutral))
        !    totalweight = totalweight+SUMMmass(nneutral+1:nindcomp)*x(nneutral+1:nindcomp))
        !    wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        !    wtf(nneutral+1:nindcomp) = x(nneutral+1:nindcomp)*Mmass(nneutral+1:nindcomp)/totalweight
        case(211) !PEG-1000  0.708913:0.291087 molar mixing ratio:
            x(1) = 1.0_wp-inputconc(2)-inputconc(3)
            x(2) = inputconc(2)*0.708913_wp
            x(3) = inputconc(2)*0.291087_wp
            totalweight = sum(Mmass(1:nneutral)*x(1:nneutral))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        case(212) !PEG-600 as a mixture of 0.788926 PEG-600-n12 and 0.211074 PEG-600-n13
            x(1) = 1.0_wp-inputconc(2)-inputconc(3)
            x(2) = inputconc(2)*0.788926_wp
            x(3) = inputconc(2)*0.211074_wp
            totalweight = sum(Mmass(1:nneutral)*x(1:nneutral))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        case(213) !PEG-1450 (PEG-1450 as a mixture of 0.493898 PEG-1450-n31 and 0.506102 PEG-1450-n32)
            x(1) = 1.0_wp-inputconc(2)-inputconc(3)
            x(2) = inputconc(2)*0.493898_wp
            x(3) = inputconc(2)*0.506102_wp
            totalweight = sum(Mmass(1:nneutral)*x(1:nneutral))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        !case(214) !PEG-1540 molar ratios
        !    x(1) = 1.0_wp-inputconc(2)-inputconc(3)
        !    x(2) = inputconc(2)*0.45089529 !PEG-1540 n = 33
        !    x(3) = inputconc(2)*0.54910471 !PEG-1540 n = 34
        !    x(4) = inputconc(3) !AS
        !    totalweight = sum(Mmass(1:nneutral)*x(1:nneutral))
        !    totalweight = totalweight+sum(Mmass(nneutral+1:nindcomp)*x(nneutral+1:nindcomp))
        !    wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        !    wtf(nneutral+1:nindcomp) = x(nneutral+1:nindcomp)*Mmass(nneutral+1:nindcomp)/totalweight
        case(510)  !ternary neutral mixture (no salts)
            x(2:nindcomp) = inputconc(2:nindcomp)  !mole fraction (with respect to salts not dissociated into ions) of other components including salts!
            x(1) = 1.0_wp-sum(x(2:nindcomp)) !for component water usually
            totalweight = sum(Mmass(1:nneutral)*x(1:nneutral))
            wtf(2:nneutral) = x(2:nneutral)*Mmass(2:nneutral)/totalweight
            wtf(1) = 1.0_wp-sum(wtf(2:nindcomp))
        case(1157) !PEG-200; mole fractions of PEG-200-n3 and PEG-200-n4 to yield 200 g/mol:  0.8689391*(PEG-200-n3) + 0.1310609*(PEG-200-n4);
            x(1) = 1.0_wp-inputconc(2)
            x(2) = inputconc(2)*0.8689391_wp !PEG-200 n = 3
            x(3) = inputconc(2)*0.1310609_wp !PEG-200 n = 4
            totalweight = sum(Mmass(1:nneutral)*x(1:nneutral))
            wtf(1:nneutral) = x(1:nneutral)*Mmass(1:nneutral)/totalweight
        !case () ! shift 2 to 3 for CO2(aq)
        !    x(3) = x(2)
        !    x(2) = 0.0_wp
        !    x(4:) = 0.0_wp
        case default
            defaultcase = .true.
        end select
    endif

    end subroutine SpecialInputConcConversion
!==========================================================================================================
    
    
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
    pure subroutine MassFrac2MoleFracMolality(wtf, XrespSalt, mrespSalt)

    implicit none
    !interface variables:
    real(wp),dimension(:),intent(in)  :: wtf         !mass fractions (input); typically these are mass fractions with electrolytes
                                                    !expressed as undissociated ions (matters for the output XrespSalt values)
    real(wp),dimension(:),intent(out) :: XrespSalt   !mole fraction of the components with respect to undissociated electrolytes
    real(wp),dimension(:),intent(out) :: mrespSalt   !molality of the (undissociated) components
    !local variables:
    integer :: i, nnp1
    real(wp) :: sum_wtfbyMmass, sum_saltfreeWTF
    real(wp),dimension(size(wtf)) :: saltfreeWTF, wtfbyMmass
    !................................................................

    nnp1 = nneutral+1
    XrespSalt = 0.0_wp
    mrespSalt = 0.0_wp
    wtfbyMmass = wtf/Mmass              !this is equivalent to sum(n_j) / sum(n_j M_j)
    sum_wtfbyMmass = sum(wtfbyMmass)

    saltfreeWTF = 0.0_wp
    sum_saltfreeWTF = sum(wtf(1:nneutral))
    if (sum_saltfreeWTF > 0.0_wp) then
        saltfreeWTF(1:nneutral) = wtf(1:nneutral) / sum_saltfreeWTF
    endif
    !total number of moles of substances is known, now one can calculate the mole fraction:
    !for the neutrals:
    XrespSalt(1:nneutral) = wtfbyMmass(1:nneutral) / sum_wtfbyMmass
    mrespSalt(1:nneutral) = wtfbyMmass(1:nneutral) / sum_saltfreeWTF
    !for the electrolytes
    if (nelectrol > 0) then
        XrespSalt(nnp1:) = wtfbyMmass(nnp1:) / sum_wtfbyMmass
        if (wtf(1) > 0.0_wp) then
            mrespSalt(nnp1:) = wtfbyMmass(nnp1:)*(saltfreeWTF(1) / wtf(1))
        else
            if (nneutral > 1 .and. any(wtf(2:nneutral) > 0.0_wp)) then
                do i = 2,nneutral
                    if (wtf(i) > 0.0_wp) then
                        exit !save the ith neutral
                    endif
                enddo
                mrespSalt(nnp1:nindcomp) = wtfbyMmass(nnp1:nindcomp)*(saltfreewtf(i) / wtf(i))
            else
                i = 777
                !!$OMP CRITICAL (MMA1)
                !write(*,*) "WARNING from MassFrac2MoleFracMolality: wtf(1) = 0.0 "
                !write(*,*) "There has to be some water in the mixture when inorganic salts are present!"
                !write(*,*) "wtf(1:nindcomp): ", wtf(1:nindcomp)
                !write(*,*) ""
                !read(*,*)
                !!$OMP end CRITICAL (MMA1)
            endif
        endif
    endif

    end subroutine MassFrac2MoleFracMolality
!========================================================================================================== 
    
      
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
    !*   -> latest changes: 2018-05-28                                                      *
    !*                                                                                      *
    !****************************************************************************************  
    pure subroutine zSolution2SpeciesMolality(zl, ml) 
    
    implicit none
    !interface variables:
    real(wp),dimension(nindcomp),intent(in) :: zl    !input mole fraction (undissociated electrolytes)
    real(wp),dimension(NKNpNGS),intent(out) :: ml    !output molalities of neutral solvent components and dissociated ions (first cations then anions according to order in Ication, Ianion)
    !local variables:
    integer :: i, cn, an, cid, aid
    real(wp) :: Msolv
    !................................
    ml = 0.0_wp
    !liquid solvent mass:
    Msolv = sum(zl(1:nneutral)*Mmass(1:nneutral))
    !molality of solvent components:
    ml(1:nneutral) = zl(1:nneutral)/Msolv
    !molality of individual ions:
    do i = 1,nelectrol  
        cn = ElectComps(i,1)    !the cation of this electrolyte
        an = ElectComps(i,2)    !the anion
        cid = CatNr(cn)         !the cation index ID within the cations of this mixture
        aid = AnNr(an)          !the anion index ID
        ml(nneutral+cid) = ml(nneutral+cid) + zl(nneutral+i)*ElectNues(i,1)/Msolv
        ml(nneutral+Ncation+aid) = ml(nneutral+Ncation+aid) + zl(nneutral+i)*ElectNues(i,2)/Msolv
    enddo
    
    end subroutine zSolution2SpeciesMolality
!==========================================================================================================================

end module ModCompScaleConversion