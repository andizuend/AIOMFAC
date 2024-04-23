!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module to provide public variables and parameters characterizing the composition-  *
!*   -independent properties of the mixture system;                                     *
!*   Note: composition-dependent variables are declared elsewhere (e.g. ModAIOMFACvar). *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2005                                                            *
!*   -> latest changes: 2022-01-17                                                      *
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
!*   -  submodule SubModDefSystem : SetSystem, definemixtures, defElectrolytes,         *
!*                                  SetMolarMass;                                       *
!*                                                                                      *
!****************************************************************************************

!nneutral = number of neutral components
!nelectrol = number of (inorganic) electrolyte components
!NGN: number of different neutral groups
!NGS: number of different ions
!NG:  number of different groups (NGN + NGS)
!Ncation: number of cations present in the mixture
!Nanion:  number of anions present in the mixture
!NGI: the max(Ncation, Nanion)
!ITAB: indicates the subgroups and their amounts that are defining each component; 
!In ITAB, electrolytes are taken as components (but with differen ions as subgroups).
!ITABsr: similar to ITAB, yet electrolytes (salts) are divided into ions: each ion is listed as a separate species.
    
module ModSystemProp

use Mod_kind_param, only : wp

implicit none

integer,public,parameter :: Nmaingroups = 76     !total number of SR (UNIFAC) main groups
integer,public,parameter :: topsubno = 265       !index no. of the last subgroup considered (top boundary of ITAB subgroups)
!--
integer,public :: nd, NGS, NGN, NGI, NG, NKNpNGS, nindcomp, nneutral, nelectrol, ninput, errorflagmix
integer,public :: Nanion, Ncation                !Nanion = number of anions, Ncation = number of cations
integer,public :: idH, idHSO4, idSO4, idHCO3, idCO3, idOH, idCO2, idCa !the index locations of these ions in the cation and anion arrays (e.g., in SMC, and SMA); idH = CatNr(205), idHSO4 = AnNr(248), idSO4 = AnNr(261)
integer,dimension(:),allocatable,public :: CompN
integer,dimension(:,:),allocatable,public :: ElectComps, ElectNues
integer,dimension(:),allocatable,public :: Ication, Ianion, ElectSubs, SolvSubs, AllSubs
integer,dimension(:),allocatable,public :: Imaingroup, maingrindexofsubgr
integer,dimension(201:topsubno),public :: CatNr, AnNr    !CatNr, AnNr: saves present mixture index entry of a certain ion related to Ication or Ianion list (whether it is cation 1, 2, 3,...).
integer,dimension(:,:),allocatable,public :: ITAB, ITABsr, ITABMG, ITAB_dimflip
!--
real(wp),parameter,public :: Rgas = 8.3144598_wp              !the universal gas constant in J/(K*mol); 8.314 4598  according to NIST (2015)
real(wp),dimension(:),allocatable,public :: cationZ, anionZ  !cationZ and anionZ are the integer charges of the cations and anions in current mixture ion order (as in Ication, Ianion);
real(wp),dimension(:),allocatable,public :: OtoCratio, HtoCratio, ElectO2Cequiv 
real(wp),dimension(201:240, 241:topsubno, 1:3),public :: IAPcoeffs
real(wp),dimension(201:240, 241:topsubno),public :: KVLE_298K
real(wp),dimension(:),allocatable :: K_el, nuestoich, SubGroupMW
real(wp),dimension(:),allocatable :: Mmass                    !component molar mass in order of neutrals then electrolytes
real(wp) :: wtf_saved
!--
character(len=200),dimension(:),allocatable,public :: cpname, compname, compnameTeX  !component names in order of mixture components
character(len=16),dimension(:),allocatable,public :: ionname, ionnameTeX
character(len=3000),dimension(:),allocatable,public :: compsubgroups, compsubgroupsTeX, compsubgroupsHTML
!--
logical,public :: bisulfsyst, calcviscosity, elpresent, frominpfile, waterpresent, solvmixrefnd
logical,public :: bicarbsyst, noCO2input
logical,public :: incl_bisulfate = .false.
logical,public :: isPEGsystem                            !to mark a special case: systems containing a PEG polymer
logical,dimension(:),allocatable,public :: ElectVolatile
logical,dimension(50),public :: errorflag_clist
!----
!interfaces to subroutines in submodule SubModDefSystem:
interface
    module subroutine SetSystem(ndi, datafromfile, ninp, cpnameinp, cpsubginp)
        integer,intent(in) :: ndi 
        logical,intent(in) :: datafromfile
        integer,intent(in) :: ninp
        !optional arguments at call:
        character(len=200),dimension(ninp),intent(in),  OPTIONAL :: cpnameinp
        integer,dimension(ninp,topsubno),intent(in), OPTIONAL :: cpsubginp
    end subroutine SetSystem
    !--
    module subroutine definemixtures(ndi, ninputcomp, compID, cpsubg)
        integer,intent(in) :: ndi, ninputcomp
        integer,dimension(:),intent(in) :: compID
        integer,dimension(:,:),intent(in) :: cpsubg 
    end subroutine definemixtures
    !--
    module subroutine defElectrolytes(nneutral, NGS, nelectrol)
        integer,intent(in) :: nneutral, NGS
        integer,intent(inout) :: nelectrol
    end subroutine defElectrolytes
    !--
    pure module subroutine SetMolarMass(MolarM)
        real(wp),dimension(:),intent(out) :: MolarM
    end subroutine SetMolarMass
    !--
    end interface
!....................................................................................

!$OMP THREADPRIVATE(nd, nindcomp, nelectrol, nneutral, ninput, NGS, NGN, NGI, NG, NKNpNGS, Nanion,  &
    !$OMP & Ncation, Mmass, idH, idHSO4, idSO4, CompN, Ication, Ianion, ElectSubs, SolvSubs, AllSubs, &
    !$OMP & Imaingroup, CatNr, AnNr, ITAB, ITAB_dimflip, ITABsr, ITABMG, OtoCratio, HtoCratio, cpname, compname,  &
    !$OMP & compnameTeX, compsubgroups, compsubgroupsTeX, compsubgroupsHTML, ionname, ionnameTeX, solvmixrefnd,  &
    !$OMP & frominpfile, bisulfsyst, waterpresent, calcviscosity, elpresent, isPEGsystem, maingrindexofsubgr,   &
    !$OMP & ElectComps, ElectNues, ElectVolatile, IAPcoeffs, KVLE_298K, K_el, SubGroupMW, ElectO2Cequiv, cationZ,  &
    !$OMP & anionZ, errorflagmix, errorflag_clist, nuestoich, idHCO3, idCO3, idOH, idCO2, idCa, bicarbsyst, &
    !$OMP & noCO2input, incl_bisulfate)
    
end module ModSystemProp