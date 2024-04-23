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
    nmaxindex = maxloc(xin, DIM=1)
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
            minlocwtf = minloc(wtf(1:nindcomp), DIM=1)
            maxlocwtf = maxloc(wtf(1:nindcomp), DIM=1)
            wtf(maxlocwtf) = wtf(maxlocwtf)+wtf(minlocwtf)
            wtf(minlocwtf) = 0.0_wp
        else !there is something wrong...
            minlocwtf = minloc(wtf(1:nindcomp), DIM=1)
            write(*,*) ""
            write(*,*) "WARNING from Inputconc_to_wtf: mass fraction of a component is less then 0.0 !!"
            write(*,*) "nd, wtf(minlocwtf): ", nd, wtf(minlocwtf)
            write(*,*) ""
            !  read(*,*)
            return
        endif
    endif
    if (sum(wtf(1:nneutral)) < lowval .AND. sum(wtf(nneutral+1:nindcomp)) > lowval) then  !there has to be some water in the mixture or some organic solvent!!
        wtf(2:nindcomp) = wtf(2:nindcomp)*(1.0_wp -lowval)
        wtf(1) = 1.0_wp - sum(wtf(2:nindcomp))
    endif

    end subroutine Inputconc_to_wtf
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
            if (nneutral > 1 .AND. any(wtf(2:nneutral) > 0.0_wp)) then
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