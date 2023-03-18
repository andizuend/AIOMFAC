!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Subroutine providing an interface to AIOMFAC for computing specific liquid phase   *
!*   activity coefficients and phase compositions with consideration of a single liquid *
!*   phase as solution.                                                                 *
!*   Input is a specific composition point and temperature for a previously initialized *
!*   system. Output is listed in outputvars.                                            *
!*                                                                                      *
!*   :: Methods and model references ::                                                 *
!*   The AIOMFAC model expressions and parameters are described in Zuend et al. (2008,  *
!*   Atmos. Chem. Phys.) and Zuend et al. (2011, Atmos. Chem. Phys.). Interaction       *
!*   parameters of Zuend et al. (2011) are used where they differ from the previous     *
!*   version. Additional parameters from Zuend and Seinfeld (2012), e.g. for peroxides, *
!*   are included as well. Several ions and their interactions with other species,      *
!*   including I-, IO3-, CO3--, HCO3- and CO2(aq), are included based on Yin et al.     *
!*   (2021, Atmos. Chem. Phys.).                                                        *
!*   Viscosity predictions via AIOMFAC-VISC are included based on the articles by       *
!*   Gervasi et al. (2020) and Lilek and Zuend (2021).                                  *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andreas Zuend,                                                                     *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2011                                                            *
!*   -> latest changes: 2023-03-18                                                      *
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
!****************************************************************************************
subroutine AIOMFAC_inout(inputconc, xinputtype, TKelvin, nspecies, outputvars, outputviscvars, &
    & outnames, errorflag_list, warningflag)

use Mod_NumPrec, only : wp
use ModSystemProp
use ModAIOMFACvar
use ModCompScaleConversion
use ModCalcActCoeff, only : AIOMFAC_calc
use ModFiniteDiffSens, only : DeltaActivities
use ModSubgroupProp, only : SMWA, SMWC
use, intrinsic :: IEEE_ARITHMETIC

implicit none
!interface variables: 
real(wp),dimension(nindcomp),intent(in)     :: inputconc        !inputconc = the concentration of a given input point (e.g., at an experimental data point)
real(wp),dimension(6,NKNpNGS),intent(out)   :: outputvars       !2-D output array with computed compositions and activities for each species; 
                                                                !structure is:  | mass-frac., mole-frac., molality, act.coeff., activity, ion-indicator | species-no |
real(wp),dimension(2),intent(out)           :: outputviscvars   !output array for viscosity related values: | viscosity | model sensitivity |
real(wp),intent(in)                         :: TKelvin          !the input temperature [K]
logical,intent(in)                          :: xinputtype
integer,intent(out)                         :: nspecies
character(len=*),dimension(NKNpNGS),intent(out) :: outnames
logical,dimension(size(errorflag_clist)),intent(out) :: errorflag_list
integer,intent(out)                         :: warningflag
!--
!local variables:
character(len=3) :: cn                      !this assumes a maximum three-digit component number in the system (max. 999); to be adjusted otherwise.
character(len=3) :: cino
integer :: i, ion_no, ion_indic, nc, NKSinput, NKSinputp1
logical :: onlyDeltaVisc
real(wp),parameter :: deps = 11.0_wp*(epsilon(deps)), ln10 = log(10.0_wp)
real(wp) :: wtf_cp, xi_cp, mi_cp, actcoeff_cp, a_cp, sum_ms, sum_miMi, xtolviscosity, w1perturb
real(wp),dimension(nelectrol) :: mixingratio, wtfdry
real(wp),dimension(nneutral) :: m_neutral
real(wp),dimension(nindcomp) :: inputconcZ, xinp, dact, dactcoeff, wfrac
!------------------------------------------------------------------------------------------- 
      
! Set initial values of some array variables:
errorflag_list = .false.  
warningflag = 0
outputvars = 0.0_wp
outputviscvars = 0.0_wp
NKSinput = ninput -nneutral
NKSinputp1 = NKSinput+1
inputconcZ = inputconc
if (noCO2input) then
    if (idCO2 > 0) then     !CO2 present, but not as input component;
        !shift input components since CO2 had been automatically added as the 
        !last neutral (non-electrolyte) component in SetSystem
        inputconcZ(idCO2+1:nindcomp) = inputconc(idCO2:nindcomp-1)
        inputconcZ(idCO2) = 0.0_wp
        NKSinput = NKSinput +1
        NKSinputp1 = NKSinputp1 +1
    endif
endif
wtfdry(1:NKSinput) = 1.0_wp
wtfdry(NKSinputp1:) = 0.0_wp
mixingratio(1:NKSinput) = 1.0_wp
mixingratio(NKSinputp1:) = 0.0_wp
nspecies = NKNpNGS
xtolviscosity = 0.0_wp
calcviscosity = .true.

if (nneutral < 1) then                      !leave the subroutine and indicate a problem to the calling routine
    errorflag_list(8) = .true.              !there must be at least one neutral component in the mixture!
else
    !weight (mass) fractions of the data point:
    call Inputconc_to_wtf(inputconcZ, mixingratio, wtfdry, xinputtype, wtf)
    if (nneutral > 0 .and. sum(wtf(1:nneutral)) < deps) then
        errorflag_list(3) = .true.
    endif
    if (nindcomp > 0 .and. any(wtf(1:nindcomp) < -deps)) then
        errorflag_list(4) = .true.
    endif
    if ( any( IEEE_IS_NAN(wtf(1:nindcomp)) ) ) then
        errorflag_list(5) = .true.
    endif
endif
if ( any(errorflag_list) ) then
    outnames = ""
    if (nspecies > nneutral) then
        outputvars(6,nneutral+1:) = real(ElectSubs(1:), kind=wp)
    endif
    return                                  !leave the subroutine and indicate a problem to the calling routine
endif
!.....
call MassFrac2MoleFracMolality(wtf, XrespSalt, mrespSalt)

if (calcviscosity) then
    onlyDeltaVisc = .true.
    xinp(1:nindcomp) = XrespSalt(1:nindcomp)
    call DeltaActivities(xinp, TKelvin, onlyDeltaVisc, dact, dactcoeff)     !will also call AIOMFAC_calc and compute activity coeff.
    w1perturb = 0.02_wp
    wfrac = wtf
    wfrac(1) = wfrac(1) + w1perturb
    wfrac = wfrac/(1.0_wp + w1perturb)
    call MassFrac2MoleFracMolality(wfrac, XrespSalt, mrespSalt)
    xtolviscosity = XrespSalt(1) - xinp(1)
else
    call AIOMFAC_calc(wtf, TKelvin)                                         !calculate at given mass fraction and temperature
endif
!.....

!Output of the AIOMFAC calculated values species-wise to array outputvars (ions separately):
call MassFrac2SolvMolalities(wtf, m_neutral)
sum_ms = sum(m_neutral(1:nneutral))
sum_miMi = 0.0_wp
do i = 1,NGI                                !calculate sum of ion molalities, mi, times ion molar mass, Mi; [kg]
    if (Ication(i) > 0) then
        sum_miMi = sum_miMi + smc(i)*SMWC(Ication(i)-200)*1.0E-3_wp 
    endif
    if (Ianion(i) > 0) then
        sum_miMi = sum_miMi + sma(i)*SMWA(Ianion(i)-240)*1.0E-3_wp 
    endif
enddo
do nc = 1,nspecies                          !loop over components
    if (nc <= nneutral) then                !distinguish between neutral components and inorg. ions
        write(cn,'(I0.2)') nc
        outnames(nc) = "comp_no_"//trim(cn)
        wtf_cp = wtf(nc)
        xi_cp = X(nc)
        mi_cp = mrespSalt(nc)               !molality in solvent mixture [mol/(kg solvent mix)]      
        actcoeff_cp = actcoeff_n(nc)
        if (wtf(nc) > 0.0_wp) then
            actcoeff_cp = actcoeff_n(nc) 
        else
            actcoeff_cp = 0.0_wp
        endif
        a_cp = activity(nc)
        ion_indic = 0
    else
        ion_no = ElectSubs(nc-nneutral)     !the current ion subgroup number to identify the ion as output component nc
        write(cino,'(I3.3)') ion_no
        outnames(nc) = "ion_no_"//cino
        ion_indic = ion_no
        !detect whether it is a cation or an anion:
        if (ion_no > 239) then              !anion
            i = AnNr(ion_no)                !the number i anion (storage location in sma(i) etc.)
            xi_cp = sma(i)/(sum_ms + SumIonMolalities)
            mi_cp = sma(i)
            wtf_cp = sma(i)*SMWA(Ianion(i)-240)*1.0E-3_wp/(1.0_wp + sum_miMi)
            if (sma(i) > 0.0_wp) then
                actcoeff_cp = actcoeff_a(i) !molal activity coeff.
            else
                actcoeff_cp = 0.0_wp
            endif
            if (actcoeff_a(i) >= 0.0_wp) then
               a_cp = actcoeff_a(i)*sma(i)  !molal activity
            else
               a_cp = -9999.999999_wp        !indicate a numerical problem
               errorflag_list(7) = .true. 
            endif
        else !cation
            i = CatNr(ion_no) 
            xi_cp = smc(i)/(sum_ms + SumIonMolalities)
            mi_cp = smc(i)
            wtf_cp = smc(i)*SMWC(Ication(i)-200)*1.0E-3_wp/(1.0_wp + sum_miMi)    
            if (smc(i) > 0.0_wp) then
                actcoeff_cp = actcoeff_c(i) !molal activity coeff.
            else
                actcoeff_cp = 0.0_wp
            endif
            if (actcoeff_c(i) >= 0.0_wp) then
               a_cp = actcoeff_c(i)*smc(i)  !molal activity
            else
               a_cp = -9999.999999_wp        !indicate a numerical problem 
               errorflag_list(7) = .true.
            endif
        endif 
    endif
    outputvars(1,nc) = wtf_cp                       ![-]   mass fraction of component nc
    outputvars(2,nc) = xi_cp                        ![-]   mole fraction with respect to electrolytes (salts) partially (e.g. HSO4-) or completely dissociated 
                                                    !      into individual ions (this is potentially different from the input mole fractions, 
                                                    !      which are usually with respect to undissociated electrolyte components!)
    outputvars(3,nc) = mi_cp                        ![mol/kg] molality of nc (with respect to dissociated ions); molality is defined as 
                                                    !         moles of species per kg of solvent mixture, where the mass of the solvent mixture means 
                                                    !         the cumulative mass of neutral compounds (water + organics) present in solution phase;
    outputvars(4,nc) = actcoeff_cp                  ![-]   activity coefficient of the component / species (see remarks below about mole fraction or molality basis)
    outputvars(5,nc) = a_cp                         ![-]   activity of the component (mole fraction basis for neutral compounds, 
                                                    !      molality basis for ions (with infinite dilution in water as reference state)
    outputvars(6,nc) = real(ion_indic, kind=wp)     !the indicator if this species is an inorg. ion or not: 0 = not ion, a number > 200 indicates the ion ID from the AIOMFAC list
enddo

! viscosity output
if (calcviscosity) then
    outputviscvars(1) = ln_etamix/ln10              !conversion by 1.0/ln(10) for output as log10(eta/[Pa s])
    outputviscvars(2) = xtolviscosity*partial_log10_etamix  !this is the log10-scale +/- error value to be added to the log10 viscosity value.
else
    outputviscvars(1) = -9999.9999999999_wp         !unrealistic values to signal "property not calculated"
    outputviscvars(2) = -9999.9999999999_wp
endif

!transfer certain error flags raised during calc.:
where(.not. errorflag_list)
    errorflag_list = errorflag_clist
endwhere

!check applicable temperature range and state a warning "errorflag" if violated:
!applicable range for electrolyte-containing mixtures (approx.): 288.0 to 309.0 K (298.15 +- 10 K); (strictly valid range would be 298.15 K only)
!applicable range for electrlyte-free mixtures (approx.): 280.0 to 400.0 K
if (NGS > 0) then !electrolyte-containing
    if (TKelvin > 309.0_wp .or. TKelvin < 288.0_wp) then    !set warning flag
        if (warningflag == 0) then                          !do not overwrite an existing warning when non-zero
            warningflag = 10
        endif
    endif
else
    if (TKelvin > 400.0_wp .or. TKelvin < 280.0_wp) then    !set warning flag
        if (warningflag == 0) then
            warningflag = 11
        endif
    endif
endif
if (.not. calcviscosity) then
    if (warningflag == 0) then
        warningflag = 16
    endif
endif

end subroutine AIOMFAC_inout
! ======================= end =======================================================