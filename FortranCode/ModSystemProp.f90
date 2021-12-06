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
!*   -> latest changes: 2021-09-22                                                      *
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
!*   -  SUBMODULE SubModDefSystem : SetSystem, definemixtures, defElectrolytes,         *
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
    
MODULE ModSystemProp

IMPLICIT NONE

INTEGER(4),PUBLIC,PARAMETER :: Nmaingroups = 76     !total number of SR (UNIFAC) main groups
INTEGER(4),PUBLIC,PARAMETER :: topsubno = 264       !index no. of the last subgroup considered (top boundary of ITAB subgroups)
!--
INTEGER(4),PUBLIC :: nd, NGS, NGN, NGI, NG, NKNpNGS, nindcomp, nneutral, nelectrol, ninput, errorflagmix, errorflagcalc
INTEGER(4),PUBLIC :: Nanion, Ncation                !Nanion = number of anions, Ncation = number of cations
INTEGER(4),PUBLIC :: idH, idHSO4, idSO4, idHCO3, idCO3, idOH, idCO2, idCa   !the index locations of these ions in the cation and anion arrays (e.g., in SMC, and SMA); idH = CatNr(205), idHSO4 = AnNr(248), idSO4 = AnNr(261)
INTEGER(4),DIMENSION(:),ALLOCATABLE,PUBLIC :: CompN
INTEGER(4),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: ElectComps, ElectNues
INTEGER(4),DIMENSION(:),ALLOCATABLE,PUBLIC :: Ication, Ianion, ElectSubs, SolvSubs, AllSubs
INTEGER(4),DIMENSION(:),ALLOCATABLE,PUBLIC :: Imaingroup, maingrindexofsubgr
INTEGER(4),DIMENSION(201:topsubno),PUBLIC :: CatNr, AnNr    !CatNr, AnNr: saves present mixture index entry of a certain ion related to Ication or Ianion list (whether it is cation 1, 2, 3,...).
INTEGER(4),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: ITAB, ITABsr, ITABMG, ITAB_dimflip
!--
REAL(8),PARAMETER,PUBLIC :: Rgas = 8.3144598D0              !the universal gas constant in J/(K*mol); 8.314 4598  according to NIST (2015)
REAL(8),DIMENSION(:),ALLOCATABLE,PUBLIC :: cationZ, anionZ  !cationZ and anionZ are the integer charges of the cations and anions in current mixture ion order (as in Ication, Ianion);
REAL(8),DIMENSION(:),ALLOCATABLE,PUBLIC :: OtoCratio, HtoCratio, ElectO2Cequiv 
REAL(8),DIMENSION(201:240, 241:topsubno, 1:3),PUBLIC :: IAPcoeffs
REAL(8),DIMENSION(201:240, 241:topsubno),PUBLIC :: KVLE_298K
REAL(8),DIMENSION(:),ALLOCATABLE :: K_el, nuestoich, SubGroupMW
REAL(8),DIMENSION(:),ALLOCATABLE :: Mmass                    !component molar mass in order of neutrals then electrolytes
!--
CHARACTER(LEN=200),DIMENSION(:),ALLOCATABLE,PUBLIC :: cpname, compname, compnameTeX  !component names in order of mixture components
CHARACTER(LEN=16),DIMENSION(:),ALLOCATABLE,PUBLIC :: ionname, ionnameTeX
CHARACTER(LEN=3000),DIMENSION(:),ALLOCATABLE,PUBLIC :: compsubgroups, compsubgroupsTeX, compsubgroupsHTML
!--
LOGICAL(4),PUBLIC :: bisulfsyst, calcviscosity, elpresent, frominpfile, waterpresent, solvmixrefnd
LOGICAL(4),PUBLIC :: bicarbsyst, noCO2input
LOGICAL(4),PUBLIC :: isPEGsystem                            !to mark a special case: systems containing a PEG polymer
LOGICAL(4),DIMENSION(:),ALLOCATABLE,PUBLIC :: ElectVolatile
!----
!interfaces to subroutines in submodule SubModDefSystem:
INTERFACE
    MODULE SUBROUTINE SetSystem(ndi, datafromfile, ninp, cpnameinp, cpsubginp)
        INTEGER(4),INTENT(IN) :: ndi 
        LOGICAL(4),INTENT(IN) :: datafromfile
        INTEGER(4),INTENT(IN) :: ninp
        !optional arguments at call:
        CHARACTER(LEN=200),DIMENSION(ninp),INTENT(IN),  OPTIONAL :: cpnameinp
        INTEGER(4),DIMENSION(ninp,topsubno),INTENT(IN), OPTIONAL :: cpsubginp
    END SUBROUTINE SetSystem
    !--
    MODULE SUBROUTINE definemixtures(ndi, ninputcomp, compID, cpsubg)
        INTEGER(4),INTENT(IN) :: ndi, ninputcomp
        INTEGER(4),DIMENSION(:),INTENT(IN) :: compID
        INTEGER(4),DIMENSION(:,:),INTENT(IN) :: cpsubg 
    END SUBROUTINE definemixtures
    !--
    MODULE SUBROUTINE defElectrolytes(nneutral, NGS, nelectrol)
        INTEGER(4),INTENT(IN) :: nneutral, NGS
        INTEGER(4),INTENT(INOUT) :: nelectrol
    END SUBROUTINE defElectrolytes
    !--
    PURE MODULE SUBROUTINE SetMolarMass(MolarM)
        REAL(8),DIMENSION(:),INTENT(OUT) :: MolarM
    END SUBROUTINE SetMolarMass
    !--
END INTERFACE
!....................................................................................

!$OMP THREADPRIVATE(nd, nindcomp, nelectrol, nneutral, ninput, NGS, NGN, NGI, NG, NKNpNGS, Nanion,  &
    !$OMP & Ncation, Mmass, idH, idHSO4, idSO4, CompN, Ication, Ianion, ElectSubs, SolvSubs, AllSubs, &
    !$OMP & Imaingroup, CatNr, AnNr, ITAB, ITAB_dimflip, ITABsr, ITABMG, OtoCratio, HtoCratio, cpname, compname,  &
    !$OMP & compnameTeX, compsubgroups, compsubgroupsTeX, compsubgroupsHTML, ionname, ionnameTeX, solvmixrefnd,  &
    !$OMP & frominpfile, bisulfsyst, waterpresent, calcviscosity, elpresent, isPEGsystem, maingrindexofsubgr,   &
    !$OMP & ElectComps, ElectNues, ElectVolatile, IAPcoeffs, KVLE_298K, K_el, SubGroupMW, ElectO2Cequiv, cationZ,  &
    !$OMP & anionZ, errorflagmix, errorflagcalc, nuestoich, idHCO3, idCO3, idOH, idCO2, idCa, bicarbsyst, noCO2input)

END MODULE ModSystemProp