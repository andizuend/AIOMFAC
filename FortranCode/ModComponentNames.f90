!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module to define the character strings for the names of the different components   *
!*   or species in AIOMFAC and to compute the component names present in a mixture.     * 
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend,                                                                        *
!*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2006                                                            *
!*   -> latest changes: 2024-08-13                                                      *
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
!*   -  subroutine nametab                                                              *
!*   -  subroutine names_mix                                                            *
!*                                                                                      *
!****************************************************************************************

module ModComponentNames

implicit none

!module public vars:
character(len=60),dimension(:),allocatable,public :: NKname, NKnameTeX
character(len=500),dimension(:),allocatable,public :: NKsmiles

!========================================================================================================== 
    contains
!========================================================================================================== 
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to set the character strings for the names of the different independent *
    !*   components in a mixture.                                                           *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, 2009                                                           *
    !*   Dept. Chem. Engineering, California Institute of Technology (2009 - 2012),         *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2009                                                            *
    !*   -> latest changes: 2026-09-17                                                      *
    !*                                                                                      *
    !****************************************************************************************

    subroutine nametab()

    implicit none
    !.....................................................
    
    allocate( NKname(1500), NKnameTeX(1500), NKsmiles(1500) )

    !list of neutral component names:
    NKname = "not_defined"
    NKnameTeX = "not_defined"
    NKsmiles = "not_defined"    !not predefined list

    !Alkanes
    NKname(1) = "Methane";        NKnameTeX(1) = "Methane";        NKsmiles(1) = "C";
    NKname(2) = "Ethane";         NKnameTeX(2) = "Ethane";         NKsmiles(2) = "CC";
    NKname(3) = "Propane";        NKnameTeX(3) = "Propane";        NKsmiles(3) = "CCC";
    NKname(4) = "Butane";         NKnameTeX(4) = "Butane";         NKsmiles(4) = "CCCC";
    NKname(5) = "Pentane";        NKnameTeX(5) = "Pentane";        NKsmiles(5) = "CCCCC";
    NKname(6) = "Hexane";         NKnameTeX(6) = "Hexane";         NKsmiles(6) = "CCCCCC";
    NKname(7) = "Heptane";        NKnameTeX(7) = "Heptane";        NKsmiles(7) = "CCCCCCC";
    NKname(8) = "Octane";         NKnameTeX(8) = "Octane";         NKsmiles(8) = "CCCCCCCC";
    NKname(9) = "Nonane";         NKnameTeX(9) = "Nonane";         NKsmiles(9) = "CCCCCCCCC";
    NKname(10) = "Decane";        NKnameTeX(10) = "Decane";        NKsmiles(10) = "CCCCCCCCCC";
    NKname(11) = "Undecane";      NKnameTeX(11) = "Undecane";      NKsmiles(11) = "CCCCCCCCCCC";
    NKname(12) = "Dodecane";      NKnameTeX(12) = "Dodecane";      NKsmiles(12) = "CCCCCCCCCCCC";
    NKname(13) = "Cyclohexane";   NKnameTeX(13) = "Cyclohexane";   NKsmiles(13) = "C1CCCCC1";
    NKname(14) = "Tetradecane";   NKnameTeX(14) = "Tetradecane";   NKsmiles(14) = "CCCCCCCCCCCCCC";
    NKname(15) = "Pentadecane";   NKnameTeX(15) = "Pentadecane";   NKsmiles(15) = "CCCCCCCCCCCCCCC";
    NKname(16) = "Hexadecane";    NKnameTeX(16) = "Hexadecane";    NKsmiles(16) = "CCCCCCCCCCCCCCCC";
    NKname(17) = "Heptadecane";   NKnameTeX(17) = "Heptadecane";   NKsmiles(17) = "CCCCCCCCCCCCCCCCC";
    NKname(18) = "Octadecane";    NKnameTeX(18) = "Octadecane";    NKsmiles(18) = "CCCCCCCCCCCCCCCCCC";
    NKname(19) = "Nonadecane";    NKnameTeX(19) = "Nonadecane";    NKsmiles(19) = "CCCCCCCCCCCCCCCCCCC";
    NKname(20) = "Icosane";       NKnameTeX(20) = "Icosane";       NKsmiles(20) = "CCCCCCCCCCCCCCCCCCCC";
    NKname(21) = "Henicosane";    NKnameTeX(21) = "Henicosane";    NKsmiles(21) = "CCCCCCCCCCCCCCCCCCCCC";
    NKname(22) = "Docosane";      NKnameTeX(22) = "Docosane";      NKsmiles(22) = "CCCCCCCCCCCCCCCCCCCCCC";
    NKname(23) = "Tricosane";     NKnameTeX(23) = "Tricosane";     NKsmiles(23) = "CCCCCCCCCCCCCCCCCCCCCCC";
    NKname(24) = "Tetracosane";   NKnameTeX(24) = "Tetracosane";   NKsmiles(24) = "CCCCCCCCCCCCCCCCCCCCCCCC";
    NKname(25) = "Pentacosane";   NKnameTeX(25) = "Pentacosane";   NKsmiles(25) = "CCCCCCCCCCCCCCCCCCCCCCCCC";
    NKname(26) = "Hexacosane";    NKnameTeX(26) = "Hexacosane";    NKsmiles(26) = "CCCCCCCCCCCCCCCCCCCCCCCCCC";
    
    NKname(30) = "Squalane";      NKnameTeX(30) = "Squalane";      NKsmiles(30) = "CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C";
    NKname(33) = "Tridecane";     NKnameTeX(33) = "Tridecane";     NKsmiles(33) = "CCCCCCCCCCCCC";

    !...
    !Alkenes
    NKname(31) = "Squalene";      NKnameTeX(31) = NKname(31);      NKsmiles(31) = "CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCCC=C(C)CCC=C(C)C"
    
    !Squalene ozonolysis products
    NKname(51) = "SqualeneO3p01";  NKnameTeX(51) = NKname(51);  !C/C(CCC1OC(C)(C)OO1)=C\CC/C=C(CC/C=C(C)/CCC=O)\C
    NKname(52) = "SqualeneO3p02";  NKnameTeX(52) = NKname(52);
    NKname(53) = "SqualeneO3p03";  NKnameTeX(53) = NKname(53);
    NKname(54) = "SqualeneO3p04";  NKnameTeX(54) = NKname(54);
    NKname(55) = "SqualeneO3p05";  NKnameTeX(55) = NKname(55);
    NKname(56) = "SqualeneO3p06";  NKnameTeX(56) = NKname(56);
    NKname(57) = "SqualeneO3p07";  NKnameTeX(57) = NKname(57);
    NKname(58) = "SqualeneO3p08";  NKnameTeX(58) = NKname(58);
    NKname(59) = "SqualeneO3p09";  NKnameTeX(59) = NKname(59);
    NKname(60) = "SqualeneO3p10";  NKnameTeX(60) = NKname(60);
    NKname(61) = "SqualeneO3p11";  NKnameTeX(61) = NKname(61);
    NKname(62) = "SqualeneO3p12";  NKnameTeX(62) = NKname(62);
    NKname(63) = "SqualeneO3p13";  NKnameTeX(63) = NKname(63);
    NKname(64) = "SqualeneO3p14";  NKnameTeX(64) = NKname(64);
    NKname(65) = "SqualeneO3p15";  NKnameTeX(65) = NKname(65);
    NKname(66) = "SqualeneO3p16";  NKnameTeX(66) = NKname(66);
    NKname(67) = "SqualeneO3p17";  NKnameTeX(67) = NKname(67);
    NKname(68) = "SqualeneO3p18";  NKnameTeX(68) = NKname(68);
    NKname(69) = "SqualeneO3p19";  NKnameTeX(69) = NKname(69);
    
    !...
    !Alcohols
    NKname(101) = "Methanol";            NKnameTeX(101) = "Methanol";            NKsmiles(101) = "CO";
    NKname(102) = "Ethanol";             NKnameTeX(102) = "Ethanol";             NKsmiles(102) = "CCO";
    NKname(103) = "1-Propanol";          NKnameTeX(103) = "1-Propanol";          NKsmiles(103) = "CCCO";
    NKname(104) = "1-Butanol";           NKnameTeX(104) = "1-Butanol";           NKsmiles(104) = "CCCCO";
    NKname(105) = "1-Pentanol";          NKnameTeX(105) = "1-Pentanol";          NKsmiles(105) = "CCCCCO";
    NKname(106) = "1-Hexanol";           NKnameTeX(106) = "1-Hexanol";           NKsmiles(106) = "CCCCCCO";
    NKname(131) = "2-Propanol";          NKnameTeX(131) = "2-Propanol";          NKsmiles(131) = "CC(C)O";
    NKname(132) = "2-Butanol";           NKnameTeX(132) = "2-Butanol";           NKsmiles(132) = "CCC(C)O";
    NKname(133) = "Isobutanol";          NKnameTeX(133) = "Isobutanol";          NKsmiles(133) = "CC(C)CO";
    NKname(134) = "tert-Butanol";        NKnameTeX(134) = "$tert$-Butanol";      NKsmiles(134) = "CC(C)(C)O";
    NKname(135) = "2-Pentanol";          NKnameTeX(135) = "2-Pentanol";          NKsmiles(135) = "CCCC(C)O";
    NKname(136) = "3-Pentanol";          NKnameTeX(136) = "3-Pentanol";          NKsmiles(136) = "CCC(CC)O";
    NKname(137) = "2-Methyl-2-butanol";  NKnameTeX(137) = "2-Methyl-2-butanol";  NKsmiles(137) = "CCC(C)(C)O";
    NKname(138) = "3-Methyl-1-butanol";  NKnameTeX(138) = "3-Methyl-1-butanol";  NKsmiles(138) = "CC(C)CCO";
    NKname(139) = "Cyclopentanol";       NKnameTeX(139) = "Cyclopentanol";       NKsmiles(139) = "OC1CCCC1";
    NKname(140) = "2-Hexanol";           NKnameTeX(140) = "2-Hexanol";           NKsmiles(140) = "CC(O)CCCC";
    NKname(141) = "3-Hexanol";           NKnameTeX(141) = "3-Hexanol";           NKsmiles(141) = "CCC(O)CCC";
    NKname(143) = "Cyclohexanol";        NKnameTeX(143) = "Cyclohexanol";        NKsmiles(143) = "OC1CCCCC1";
    NKname(147) = "2-Octanol";           NKnameTeX(147) = "2-Octanol";           NKsmiles(147) = "CC(O)CCCCCC";

    !142 = 4-methyl-2-pentanol, 144 = 2-heptanol, 
    !    145 = 3-heptanol, 146 = 4-heptanol, 147 = 2-octanol, 148 = 3-octanol,
    !    149 = 4-octanol, 
    NKname(150) = "2-Methyl-1-butanol";  NKnameTeX(150) = "2-Methyl-1-butanol";  NKsmiles(150) = "CCC(C)CO";

    !...
    !(mono) carboxylic acids
    NKname(201) = "Formic_acid";             NKnameTeX(201) = "Formic_acid";             NKsmiles(201) = "C(=O)O";
    NKname(202) = "Acetic_acid";             NKnameTeX(202) = "Acetic_acid";             NKsmiles(202) = "CC(=O)O";
    NKname(203) = "Propanoic_acid";          NKnameTeX(203) = "Propanoic_acid";          NKsmiles(203) = "CCC(=O)O";
    NKname(204) = "Butyric_acid";            NKnameTeX(204) = "Butyric_acid";            NKsmiles(204) = "CCCC(=O)O";
    NKname(210) = "2-Methylpropanoic_acid";  NKnameTeX(210) = "2-Methylpropanoic_acid";  NKsmiles(210) = "CC(C)C(=O)O";
    NKname(220) = "Pyruvic_acid";            NKnameTeX(220) = "Pyruvic_acid";            NKsmiles(220) = "CC(=O)C(=O)O";
    NKname(221) = "Methacrylic_acid";        NKnameTeX(221) = "Methacrylic_acid";        NKsmiles(221) = "CC(=C)C(=O)O";
    NKname(222) = "Palmitic_acid";           NKnameTeX(222) = "Palmitic_acid";           NKsmiles(222) = "CCCCCCCCCCCCCCCC(=O)O";
    NKname(223) = "Stearic_acid";            NKnameTeX(223) = "Stearic_acid";            NKsmiles(223) = "CCCCCCCCCCCCCCCCCC(=O)O";
    NKname(224) = "Oleic_acid";              NKnameTeX(224) = "Oleic_acid";              NKsmiles(224) = "CCCCCCCCC=CCCCCCCCC(=O)O";
    NKname(225) = "Docosanoic_acid";         NKnameTeX(225) = "Docosanoic_acid";         NKsmiles(225) = "CCCCCCCCCCCCCCCCCCCCCC(=O)O";
    NKname(226) = "Tetracosanoic_acid";      NKnameTeX(226) = "Tetracosanoic_acid";      NKsmiles(226) = "CCCCCCCCCCCCCCCCCCCCCCCC(=O)O";
    NKname(227) = "Sapienic_acid";           NKnameTeX(227) = "Sapienic_acid";           NKsmiles(227) = "CCCCCCCC=CCCCCCCCC(=O)O";
    NKname(228) = "Gluconic_acid";           NKnameTeX(228) = "Gluconic_acid";           NKsmiles(228) = "OCC(O)C(O)C(O)C(O)C(=O)O";          
    NKname(229) = "Levulinic_acid";          NKnameTeX(229) = "Levulinic_acid";          NKsmiles(229) = "CC(=O)CCC(=O)O";

    !...
    !Amino acids
    NKname(241) = "Serine";   NKnameTeX(241) = "Serine";   NKsmiles(241) = "NCC(O)C(=O)O";      != 2-Amino-3-hydroxypropanoic acid
    NKname(242) = "Glycine";  NKnameTeX(242) = "Glycine";  NKsmiles(242) = "NCC(=O)O";          != Aminoacetic acid
    NKname(243) = "Alanine";  NKnameTeX(243) = "Alanine";  NKsmiles(243) = "CC(N)C(=O)O";       != 2-Aminopropanoic acid
    
    !...
    !Dicarboxylic acids
    NKname(301) = "Oxalic_acid";              NKnameTeX(301) = "Oxalic_acid";              NKsmiles(301) = "O=C(O)C(=O)O";
    NKname(302) = "Malonic_acid";             NKnameTeX(302) = "Malonic_acid";             NKsmiles(302) = "O=C(O)CC(=O)O";
    NKname(303) = "Succinic_acid";            NKnameTeX(303) = "Succinic_acid";            NKsmiles(303) = "O=C(O)CCC(=O)O";
    NKname(304) = "Glutaric_acid";            NKnameTeX(304) = "Glutaric_acid";            NKsmiles(304) = "O=C(O)CCCC(=O)O";
    NKname(305) = "Adipic_acid";              NKnameTeX(305) = "Adipic_acid";              NKsmiles(305) = "O=C(O)CCCCC(=O)O";
    NKname(306) = "Pimelic_acid";             NKnameTeX(306) = "Pimelic_acid";             NKsmiles(306) = "O=C(O)CCCCCC(=O)O";
    NKname(307) = "Suberic_acid";             NKnameTeX(307) = "Suberic_acid";             NKsmiles(307) = "O=C(O)CCCCCCC(=O)O";
    NKname(308) = "Methylsuccinic_acid";      NKnameTeX(308) = "Methylsuccinic_acid";      NKsmiles(308) = "CC(C(=O)O)CC(=O)O";
    NKname(309) = "Dimethylmalonic_acid";     NKnameTeX(309) = "Dimethylmalonic_acid";     NKsmiles(309) = "CC(C)(C(=O)O)C(=O)O";
    NKname(310) = "2-Methylglutaric_acid";    NKnameTeX(310) = "2-Methylglutaric_acid";    NKsmiles(310) = "CC(C(=O)O)CCC(=O)O";
    NKname(311) = "3-Methylglutaric_acid";    NKnameTeX(311) = "3-Methylglutaric_acid";    NKsmiles(311) = "O=C(O)CC(C)CC(=O)O";
    NKname(312) = "2,2-Dimethylsuccinic_acid";NKnameTeX(312) = "2,2-Dimethylsuccinic_acid";NKsmiles(312) = "CC(C)(C(=O)O)CC(=O)O";
    NKname(313) = "3-Methyladipic_acid";      NKnameTeX(313) = "3-Methyladipic_acid";      NKsmiles(313) = "O=C(O)CC(C)CCC(=O)O";
    NKname(314) = "3,3-Dimethylglutaric_acid";NKnameTeX(314) = "3,3-Dimethylglutaric_acid";NKsmiles(314) = "O=C(O)CC(C)(C)CC(=O)O";
    NKname(315) = "Diethylmalonic_acid";      NKnameTeX(315) = "Diethylmalonic_acid";      NKsmiles(315) = "CCC(C(CC)C(=O)O)C(=O)O";
    NKname(316) = "Tartaric_acid";            NKnameTeX(316) = "Tartaric_acid";            NKsmiles(316) = "O=C(O)C(O)C(O)C(=O)O";
    NKname(317) = "Methylmalonic_acid";       NKnameTeX(317) = "Methylmalonic_acid";       NKsmiles(317) = "CC(C(=O)O)C(=O)O";
    NKname(318) = "Tartronic_acid";           NKnameTeX(318) = "Tartronic_acid";           NKsmiles(318) = "O=C(O)C(O)C(=O)O";

    NKname(320) = "Azelaic_acid";             NKnameTeX(320) = "Azelaic_acid";             NKsmiles(320) = "O=C(O)CCCCCCCC(=O)O";
    NKname(321) = "Sebaic_acid";              NKnameTeX(321) = "Sebaic_acid";              NKsmiles(321) = "O=C(O)CCCCCCCCC(=O)O";
    NKname(322) = "Dodecanedioic_acid";       NKnameTeX(322) = "Dodecanedioic_acid";       NKsmiles(322) = "O=C(O)CCCCCCCCCCC(=O)O";

    !...
    !Water
    NKname(401) = "H2O";    NKnameTeX(401) = "H$_2$O";    NKsmiles(401) = "O"
    NKname(402) = "CO2";    NKnameTeX(402) = "CO$_2$";    NKsmiles(402) = "O=C=O" !402 for carbon dioxide CO2(aq)
    !...
    !hydroperoxides
    NKname(411) = "hydroperoxyundecane";       NKnameTeX(411) = "hydroperoxyundecane";       NKsmiles(411) = "CCCCCCCCCCCOO";
    NKname(412) = "hydroperoxydodecane";       NKnameTeX(412) = "hydroperoxydodecane";       NKsmiles(412) = "CCCCCCCCCCCCOO";
    NKname(413) = "hydroperoxytridecane";      NKnameTeX(413) = "hydroperoxytridecane";      NKsmiles(413) = "CCCCCCCCCCCCCOO";
    NKname(414) = "hydroperoxytetradecane";    NKnameTeX(414) = "hydroperoxytetradecane";    NKsmiles(414) = "CCCCCCCCCCCCCCOO";
    NKname(415) = "hydroperoxypentadecane";    NKnameTeX(415) = "hydroperoxypentadecane";    NKsmiles(415) = "CCCCCCCCCCCCCCCOO";
    NKname(416) = "hydroperoxyhexadecane";     NKnameTeX(416) = "hydroperoxyhexadecane";     NKsmiles(416) = "CCCCCCCCCCCCCCCCOO";
    NKname(417) = "hydroperoxyheptadecane";    NKnameTeX(417) = "hydroperoxyheptadecane";    NKsmiles(417) = "CCCCCCCCCCCCCCCCCOO";
    NKname(418) = "hydroperoxyoctadecane";     NKnameTeX(418) = "hydroperoxyoctadecane";     NKsmiles(418) = "CCCCCCCCCCCCCCCCCCOO";
    NKname(419) = "hydroperoxynonadecane";     NKnameTeX(419) = "hydroperoxynonadecane";     NKsmiles(419) = "CCCCCCCCCCCCCCCCCCCOO";
    NKname(420) = "hydroperoxyicosane";        NKnameTeX(420) = "hydroperoxyicosane";        NKsmiles(420) = "CCCCCCCCCCCCCCCCCCCCOO";
    NKname(421) = "hydroperoxyhenicosane";     NKnameTeX(421) = "hydroperoxyhenicosane";     NKsmiles(421) = "CCCCCCCCCCCCCCCCCCCCCOO";
    NKname(422) = "hydroperoxydocosane";       NKnameTeX(422) = "hydroperoxydocosane";       NKsmiles(422) = "CCCCCCCCCCCCCCCCCCCCCCOO";
    NKname(423) = "hydroperoxytricosane";      NKnameTeX(423) = "hydroperoxytricosane";      NKsmiles(423) = "CCCCCCCCCCCCCCCCCCCCCCCOO";
    NKname(424) = "hydroperoxytetracosane";    NKnameTeX(424) = "hydroperoxytetracosane";    NKsmiles(424) = "CCCCCCCCCCCCCCCCCCCCCCCCOO";
    NKname(425) = "hydroperoxypentacosane";    NKnameTeX(425) = "hydroperoxypentacosane";    NKsmiles(425) = "CCCCCCCCCCCCCCCCCCCCCCCCCOO";

    !...
    !dihydroperoxides
    NKname(426) = "dihydroperoxyundecane";  NKnameTeX(426) = "dihydroperoxyundecane";  
    NKname(427) = "dihydroperoxydodecane";  NKnameTeX(427) = "dihydroperoxyundecane";  
    NKname(428) = "dihydroperoxytridecane";  NKnameTeX(428) = "dihydroperoxytridecane";  
    NKname(429) = "dihydroperoxytetradecane";  NKnameTeX(429) = "dihydroperoxytetradecane";   
    NKname(430) = "dihydroperoxypentadecane";  NKnameTeX(430) = "dihydroperoxypentadecane";  
    NKname(431) = "dihydroperoxyhexadecane";  NKnameTeX(431) = "dihydroperoxyhexadecane";   
    NKname(432) = "dihydroperoxyheptadecane";  NKnameTeX(432) = "dihydroperoxyheptadecane";  
    NKname(433) = "dihydroperoxyoctadecane";  NKnameTeX(433) = "dihydroperoxyoctadecane";  
    NKname(434) = "dihydroperoxynonadecane";  NKnameTeX(434) = "dihydroperoxynonadecane";  
    NKname(435) = "dihydroperoxyicosane";  NKnameTeX(435) = "dihydroperoxyicosane";  
    NKname(436) = "dihydroperoxyhenicosane";  NKnameTeX(436) = "dihydroperoxyhenicosane";  
    NKname(437) = "dihydroperoxydocosane";  NKnameTeX(437) = "dihydroperoxydocosane";  
    NKname(438) = "dihydroperoxytricosane";  NKnameTeX(438) = "dihydroperoxytricosane";  

    !hydroperoxyketones
    NKname(441) = "hydroperoxyundecanone";  NKnameTeX(441) = "hydroperoxyundecanone";  
    NKname(442) = "hydroperoxydodecanone";  NKnameTeX(442) = "hydroperoxydodecanone";  
    NKname(443) = "hydroperoxytridecanone";  NKnameTeX(443) = "hydroperoxytridecanone";  
    NKname(444) = "hydroperoxytetradecanone";  NKnameTeX(444) = "hydroperoxytetradecanone";  
    NKname(445) = "hydroperoxypentadecanone";  NKnameTeX(445) = "hydroperoxypentadecanone";  
    NKname(446) = "hydroperoxyhexadecanone";  NKnameTeX(446) = "hydroperoxyhexadecanone";  
    NKname(447) = "hydroperoxyheptadecanone";  NKnameTeX(447) = "hydroperoxyheptadecanone";  
    NKname(448) = "hydroperoxyoctadecanone";  NKnameTeX(448) = "hydroperoxyoctadecanone";  
    NKname(449) = "hydroperoxynonadecanone";  NKnameTeX(449) = "hydroperoxynonadecanone";  
    NKname(450) = "hydroperoxyicosanone";  NKnameTeX(450) = "hydroperoxyicosanone";  
    NKname(451) = "hydroperoxyhenicosanone";  NKnameTeX(451) = "hydroperoxyhenicosanone";  
    NKname(452) = "hydroperoxydocosanone";  NKnameTeX(452) = "hydroperoxydocosanone";  
    NKname(453) = "hydroperoxytricosanone";  NKnameTeX(453) = "hydroperoxytricosanone";  
    !...
    !hydroperoxydiketones
    NKname(454) = "hydroperoxyundecanedione";  NKnameTeX(454) = "hydroperoxyundecanedione";  
    NKname(455) = "hydroperoxydodecanedione";  NKnameTeX(455) = "hydroperoxydodecanedione";  
    NKname(456) = "hydroperoxytridecanedione";  NKnameTeX(456) = "hydroperoxytridecanedione";  
    NKname(457) = "hydroperoxytetradecanedione";  NKnameTeX(457) = "hydroperoxytetradecanedione";  
    NKname(458) = "hydroperoxypentadecanedione";  NKnameTeX(458) = "hydroperoxypentadecanedione";  
    NKname(459) = "hydroperoxyhexadecanedione";  NKnameTeX(459) = "hydroperoxyhexadecanedione";  
    NKname(460) = "hydroperoxyheptadecanedione";  NKnameTeX(460) = "hydroperoxyheptadecanedione";  
    NKname(461) = "hydroperoxyoctadecanedione";  NKnameTeX(461) = "hydroperoxyoctadecanedione";  
    NKname(462) = "hydroperoxynonadecanedione";  NKnameTeX(462) = "hydroperoxynonadecanedione";  
    NKname(463) = "hydroperoxyicosanedione";  NKnameTeX(463) = "hydroperoxyicosanedione";  
    NKname(464) = "hydroperoxyhenicosanedione";  NKnameTeX(464) = "hydroperoxyhenicosanedione";  
    NKname(465) = "hydroperoxydocosanedione";  NKnameTeX(465) = "hydroperoxydocosanedione";  
    NKname(466) = "hydroperoxytricosanedione";  NKnameTeX(466) = "hydroperoxytricosanedione";  

    !dihydroperoxytetraketones
    NKname(467) = "dihydroperoxytridecanetetraone";  NKnameTeX(467) = "dihydroperoxytridecanetetraone";  
    NKname(468) = "dihydroperoxytetradecanetetraone";  NKnameTeX(468) = "dihydroperoxytetradecanetetraone";  
    NKname(469) = "dihydroperoxypentadecanetetraone";  NKnameTeX(469) = "dihydroperoxypentadecanetetraone";  
    NKname(470) = "dihydroperoxyhexadecanetetraone";  NKnameTeX(470) = "dihydroperoxyhexadecanetetraone";  
    NKname(471) = "dihydroperoxyheptadecanetetraone";  NKnameTeX(471) = "dihydroperoxyheptadecanetetraone";  
    NKname(472) = "dihydroperoxyoctadecanetetraone";  NKnameTeX(472) = "dihydroperoxyoctadecanetetraone";  
    NKname(473) = "dihydroperoxynonadecanetetraone";  NKnameTeX(473) = "dihydroperoxynonadecanetetraone";  
    NKname(474) = "dihydroperoxyicosanetetraone";  NKnameTeX(474) = "dihydroperoxyicosanetetraone";  
    NKname(475) = "dihydroperoxyhenicosanetetraone";  NKnameTeX(475) = "dihydroperoxyhenicosanetetraone";  
    NKname(476) = "dihydroperoxydocosanetetraone";  NKnameTeX(476) = "dihydroperoxydocosanetetraone";  
    NKname(477) = "dihydroperoxytricosanetetraone";  NKnameTeX(477) = "dihydroperoxytricosanetetraone";  
    !...
    !set of bitumen with extra oxidations to account for longer carbon chains than C12
    !hydroperoxytriketones
    NKname(478) = "hydroperoxyundecanetrione";  NKnameTeX(478) = "hydroperoxyundecanetrione";  
    NKname(479) = "hydroperoxydodecanetrione";  NKnameTeX(479) = "hydroperoxydodecanetrione";  
    NKname(480) = "hydroperoxytridecanetrione";  NKnameTeX(480) = "hydroperoxytridecanetrione";  
    NKname(481) = "hydroperoxytetradecanetrione";  NKnameTeX(481) = "hydroperoxytetradecanetrione";  
    NKname(482) = "hydroperoxypentadecanetrione";  NKnameTeX(482) = "hydroperoxypentadecanetrione";  
    NKname(483) = "hydroperoxyhexadecanetrione";  NKnameTeX(483) = "hydroperoxyhexadecanetrione";  
    NKname(484) = "hydroperoxyheptadecanetrione";  NKnameTeX(484) = "hydroperoxyheptadecanetrione";  
    NKname(485) = "hydroperoxyoctadecanetrione";  NKnameTeX(485) = "hydroperoxyoctadecanetrione";  
    NKname(486) = "hydroperoxynonadecanetrione";  NKnameTeX(486) = "hydroperoxynonadecanetrione";  
    NKname(487) = "hydroperoxyicosanetrione";  NKnameTeX(487) = "hydroperoxyicosanetrione";  
    NKname(488) = "hydroperoxyhenicosanetrione";  NKnameTeX(488) = "hydroperoxyhenicosanetrione";  
    NKname(489) = "hydroperoxydocosanetrione";  NKnameTeX(489) = "hydroperoxydocosanetrione";  
    NKname(490) = "hydroperoxytricosanetrione";  NKnameTeX(490) = "hydroperoxytricosanetrione";  

    !...
    !Polycarboxylic acids, functionalized aliphatic acids 
    NKname(501) = "Maleic_acid";              NKnameTeX(501) = "Maleic_acid";              NKsmiles(501) = "O=C(O)C=CC(=O)O";
    NKname(502) = "Fumaric_acid";             NKnameTeX(502) = "Fumaric_acid";             NKsmiles(502) = "O=C(O)C=CC(=O)O";           !NOTE: ignoring stereoisomerism, fumaric and maleic are equivalent
    NKname(503) = "Malic_acid";               NKnameTeX(503) = "Malic_acid";               NKsmiles(503) = "O=C(O)CC(O)C(=O)O";
    NKname(504) = "Tartaric_acid";            NKnameTeX(504) = "Tartaric_acid";            NKsmiles(504) = "O=C(O)C(O)C(O)C(=O)O";
    NKname(505) = "Citric_acid";              NKnameTeX(505) = "Citric_acid";              NKsmiles(505) = "O=C(O)CC(O)(CC(=O)O)C(=O)O";
    NKname(506) = "alpha-Ketoglutaric_acid";  NKnameTeX(506) = "$\alpha$-Ketoglutaric_acid"; NKsmiles(506) = "O=C(O)CCC(=O)C(=O)O";
    NKname(507) = "Glycolic_acid";            NKnameTeX(507) = "Glycolic_acid";            NKsmiles(507) = "OCC(=O)O";
    NKname(508) = "Methoxyacetic_acid";       NKnameTeX(508) = "Methoxyacetic_acid";       NKsmiles(508) = "COCC(=O)O";
    NKname(509) = "Shikimic_acid";            NKnameTeX(509) = "Shikimic_acid";            NKsmiles(509) = "O=C(O)C1=CC(O)C(O)C1O";
    NKname(510) = "Isopimaric_acid";          NKnameTeX(510) = "Isopimaric_acid";          NKsmiles(510) = "CC(C)C1CCC2(C)C3CCC(C)(C)C3=CCC2C1C(=O)O";
    NKname(511) = "Abietic_acid";             NKnameTeX(511) = "Abietic_acid";             NKsmiles(511) = "CC(C)C1=CC2CC(C)(C)C3CCC(C)(C)C3C2CC1C(=O)O";

    !...
    !dihydroperoxydiketones
    NKname(512) = "dihydroperoxydodecanedione";  NKnameTeX(512) = "dihydroperoxydodecanedione";  
    NKname(513) = "dihydroperoxytridecanedione";  NKnameTeX(513) = "dihydroperoxytridecanedione";  
    NKname(514) = "dihydroperoxytetradecanedione";  NKnameTeX(514) = "dihydroperoxytetradecanedione";  
    NKname(515) = "dihydroperoxypentadecanedione";  NKnameTeX(515) = "dihydroperoxypentadecanedione";  
    NKname(516) = "dihydroperoxyhexadecanedione";  NKnameTeX(516) = "dihydroperoxyhexadecanedione";  
    NKname(517) = "dihydroperoxyheptadecanedione";  NKnameTeX(517) = "dihydroperoxyheptadecanedione";  
    NKname(518) = "dihydroperoxyoctadecanedione";  NKnameTeX(518) = "dihydroperoxyoctadecanedione";  
    NKname(519) = "dihydroperoxynonadecanedione";  NKnameTeX(519) = "dihydroperoxynonadecanedione";  
    NKname(520) = "dihydroperoxyicosanedione";  NKnameTeX(520) = "dihydroperoxyicosanedione";  
    NKname(521) = "dihydroperoxyhenicosanedione";  NKnameTeX(521) = "dihydroperoxyhenicosanedione";  
    NKname(522) = "dihydroperoxydocosanedione";  NKnameTeX(523) = "dihydroperoxydocosanedione";  
    NKname(523) = "dihydroperoxytricosanedione";  NKnameTeX(524) = "dihydroperoxytricosanedione";  
    !...
    !dihydroperoxytriketones
    NKname(524) = "dihydroperoxytetradecanetrione";  NKnameTeX(524) = "dihydroperoxytetradecanetrione";  
    NKname(525) = "dihydroperoxypentadecanetrione";  NKnameTeX(525) = "dihydroperoxypentadecanetrione";  
    NKname(526) = "dihydroperoxyhexadecanetrione";  NKnameTeX(526) = "dihydroperoxyhexadecanetrione";  
    NKname(527) = "dihydroperoxyheptadecanetrione";  NKnameTeX(527) = "dihydroperoxyheptadecanetrione";  
    NKname(528) = "dihydroperoxyoctadecanetrione";  NKnameTeX(528) = "dihydroperoxyoctadecanetrione";  
    NKname(529) = "dihydroperoxynonadecanetrione";  NKnameTeX(529) = "dihydroperoxynonadecanetrione";  
    NKname(530) = "dihydroperoxyicosanetrione";  NKnameTeX(530) = "dihydroperoxyicosanetrione";  
    NKname(531) = "dihydroperoxyhenicosanetrione";  NKnameTeX(531) = "dihydroperoxyhenicosanetrione";  
    NKname(532) = "dihydroperoxydocosanetrione";  NKnameTeX(532) = "dihydroperoxydocosanetrione";  
    NKname(533) = "dihydroperoxytricosanetrione";  NKnameTeX(533) = "dihydroperoxytricosanetrione";  
    
    !...
    !Diols, Triols, Polyols, Sugars
    NKname(601) = "1,2-Ethanediol";         NKnameTeX(601) = "1,2-Ethanediol";         NKsmiles(601) = "OCCO";
    NKname(602) = "1,2-Propanediol";        NKnameTeX(602) = "1,2-Propanediol";        NKsmiles(602) = "CC(O)CO";
    NKname(603) = "1,3-Propanediol";        NKnameTeX(603) = "1,3-Propanediol";        NKsmiles(603) = "OCCCO";
    NKname(604) = "Glycerol";               NKnameTeX(604) = "Glycerol";               NKsmiles(604) = "OCC(O)CO";
    NKname(605) = "1,2-Butanediol";         NKnameTeX(605) = "1,2-Butanediol";         NKsmiles(605) = "CCC(O)CO";
    NKname(606) = "1,3-Butanediol";         NKnameTeX(606) = "1,3-Butanediol";         NKsmiles(606) = "CC(O)CCO";
    NKname(607) = "1,4-Butanediol";         NKnameTeX(607) = "1,4-Butanediol";         NKsmiles(607) = "OCCCCO";
    NKname(608) = "2,3-Butanediol";         NKnameTeX(608) = "2,3-Butanediol";         NKsmiles(608) = "CC(O)C(O)C";
    NKname(609) = "1,2,3-Butanetriol";      NKnameTeX(609) = "1,2,3-Butanetriol";      NKsmiles(609) = "CC(O)C(O)CO";
    NKname(610) = "1,2,4-Butanetriol";      NKnameTeX(610) = "1,2,4-Butanetriol";      NKsmiles(610) = "OCC(O)CCO";
    NKname(611) = "Erythritol";             NKnameTeX(611) = "Erythritol";             NKsmiles(611) = "OCC(O)C(O)CO";
    NKname(612) = "1,2-Pentanediol";        NKnameTeX(612) = "1,2-Pentanediol";        NKsmiles(612) = "CCCC(O)CO";
    NKname(613) = "1,3-Pentanediol";        NKnameTeX(613) = "1,3-Pentanediol";        NKsmiles(613) = "CCC(O)CCO";
    NKname(614) = "1,4-Pentanediol";        NKnameTeX(614) = "1,4-Pentanediol";        NKsmiles(614) = "CC(O)CCCO";
    NKname(615) = "1,5-Pentanediol";        NKnameTeX(615) = "1,5-Pentanediol";        NKsmiles(615) = "OCCCCCO";
    NKname(616) = "2,3-Pentanediol";        NKnameTeX(616) = "2,3-Pentanediol";        NKsmiles(616) = "CC(O)C(O)CC";
    NKname(617) = "2,4-Pentanediol";        NKnameTeX(617) = "2,4-Pentanediol";        NKsmiles(617) = "CC(O)CC(O)C";
    NKname(618) = "1,2-Hexanediol";         NKnameTeX(618) = "1,2-Hexanediol";         NKsmiles(618) = "CCCCC(O)CO";
    NKname(619) = "1,3-Hexanediol";         NKnameTeX(619) = "1,3-Hexanediol";         NKsmiles(619) = "CCCC(O)CCO";
    NKname(620) = "1,4-Hexanediol";         NKnameTeX(620) = "1,4-Hexanediol";         NKsmiles(620) = "CCC(O)CCCO";
    NKname(621) = "1,5-Hexanediol";         NKnameTeX(621) = "1,5-Hexanediol";         NKsmiles(621) = "CC(O)CCCCO";
    NKname(622) = "1,6-Hexanediol";         NKnameTeX(622) = "1,6-Hexanediol";         NKsmiles(622) = "OCCCCCCO";
    NKname(623) = "2,3-Hexanediol";         NKnameTeX(623) = "2,3-Hexanediol";         NKsmiles(623) = "CC(O)C(O)CCC";
    NKname(624) = "2,4-Hexanediol";         NKnameTeX(624) = "2,4-Hexanediol";         NKsmiles(624) = "CC(O)CC(O)CC";
    NKname(625) = "2,5-Hexanediol";         NKnameTeX(625) = "2,5-Hexanediol";         NKsmiles(625) = "CC(O)CCC(O)C";
    NKname(626) = "1,2,5-Hexanetriol";      NKnameTeX(626) = "1,2,5-Hexanetriol";      NKsmiles(626) = "CC(O)CCC(O)CO";
    NKname(627) = "1,2,6-Hexanetriol";      NKnameTeX(627) = "1,2,6-Hexanetriol";      NKsmiles(627) = "OCCCCC(O)CO";
    NKname(628) = "2,3,4-Hexanetriol";      NKnameTeX(628) = "2,3,4-Hexanetriol";      NKsmiles(628) = "CC(O)C(O)C(O)CC";
    NKname(629) = "Sorbitol";               NKnameTeX(629) = "Sorbitol";               NKsmiles(629) = "OCC(O)C(O)C(O)C(O)CO";
    NKname(630) = "Mannitol";               NKnameTeX(630) = "Mannitol";               NKsmiles(630) = "OCC(O)C(O)C(O)C(O)CO";
    NKname(631) = "1,4-Heptanediol";        NKnameTeX(631) = "1,4-Heptanediol";        NKsmiles(631) = "CCCC(O)CCCO";
    NKname(632) = "1,5-Heptanediol";        NKnameTeX(632) = "1,5-Heptanediol";        NKsmiles(632) = "CCC(O)CCCCO";
    NKname(633) = "1,7-Heptanediol";        NKnameTeX(633) = "1,7-Heptanediol";        NKsmiles(633) = "OCCCCCCCO";
    NKname(634) = "2,4-Heptanediol";        NKnameTeX(634) = "2,4-Heptanediol";        NKsmiles(634) = "CC(O)CC(O)CCC";
    NKname(635) = "Xylitol";                NKnameTeX(635) = "Xylitol";                NKsmiles(635) = "OCC(O)C(O)C(O)CO";
    NKname(638) = "1,3-Nonanediol";         NKnameTeX(638) = "1,3-Nonanediol";         NKsmiles(638) = "CCCCCCC(O)CCO";
    NKname(639) = "1,4-Dihydroxy-2-butene"; NKnameTeX(639) = "1,4-Dihydroxy-2-butene"; NKsmiles(639) = "OCC=CCO";
    NKname(640) = "Levoglucosan";           NKnameTeX(640) = "Levoglucosan";           NKsmiles(640) = "OC1C2COC(O2)C(O)C1O";
    NKname(641) = "D-Fructopyranose";       NKnameTeX(641) = "D-Fructopyranose";       NKsmiles(641) = "C1C(C(C(C(O1)(CO)O)O)O)O";
    NKname(642) = "D-Mannopyranose";        NKnameTeX(642) = "D-Mannopyranose";        NKsmiles(642) = "C1C(C(C(C(C(O1)O)O)O)O)O";
    !..
    NKname(643) = "1,2,10-Decanetriol";      NKnameTeX(643) = "1,2,10-Decanetriol";      NKsmiles(643) = "OCCCCCCCCC(O)CO";
    NKname(644) = "1,2,5,8-Octanetetrol";    NKnameTeX(644) = "1,2,5,8-Octanetetrol";    NKsmiles(644) = "OCCCC(O)CCC(O)CO";
    NKname(645) = "1,2,7,8-Octanetetrol";    NKnameTeX(645) = "1,2,7,8-Octanetetrol";    NKsmiles(645) = "OCC(O)CCCCC(O)CO";
    
    NKname(646) = "D-Ribofuranose";          NKnameTeX(646) = "D-Ribofuranose";          NKsmiles(646) = "C([CH]1[CH)O";
    NKname(647) = "1,2,6-Hexanetriol";       NKnameTeX(647) = "1,2,6-Hexanetriol";       NKsmiles(647) = "OCCCCC(O)CO";
    NKname(648) = "2-Methylerythritol";      NKnameTeX(648) = "2-Methylerythritol";      NKsmiles(648) = "CC(CO)(O)C(O)CO";
    
    NKname(651) = "Xylose";                  NKnameTeX(651) = "Xylose";                  NKsmiles(651) = "OCC(O)C(O)C(O)C=O";
    NKname(652) = "Glucose";                 NKnameTeX(652) = "Glucose";                 NKsmiles(652) = "OCC(O)C(O)C(O)C(O)C=O";
    NKname(653) = "Fructose";                NKnameTeX(653) = "Fructose";                NKsmiles(653) = "OCC(O)C(O)C(O)C(O)C(=O)CO";
    NKname(654) = "Sucrose";                 NKnameTeX(654) = "Sucrose";                 NKsmiles(654) = "OCC1OC(OCC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O";
    NKname(655) = "D-Glucopyranose";         NKnameTeX(655) = "D-Glucopyranose";         NKsmiles(655) = "C1C(C(C(C(C(O1)O)O)O)O)O";
    NKname(656) = "Raffinose";               NKnameTeX(656) = "Raffinose";               NKsmiles(656) = "OCC1OC(OC2C(O)C(O)C(OC3OC(CO)C(O)C(O)C3O)C(O)O2)C(O)C(O)C1O";
    NKname(657) = "2-Chloro-D-Glucose";      NKnameTeX(657) = "2-Chloro-D-Glucose";      NKsmiles(657) = "C(C1C(C(C(C(O1)O)Cl)O)O)O";
    
    !...
    !Poly (ethylene glycols),  PEG
    NKname(658) = "Triethylene_glycol";     NKnameTeX(658) = "Triethylene_glycol";      NKsmiles(658) = "OCCOCCOCCO"; 
    NKname(659) = "PEG-300-n5";  NKnameTeX(659) = "PEG-300-n5";    NKsmiles(659) = "OCCOCCOCCOCCOCCOCCO"; 
    NKname(660) = "PEG-300-n6";  NKnameTeX(660) = "PEG-300-n6";    NKsmiles(660) = "OCCOCCOCCOCCOCCOCCOCCO";  
    NKname(661) = "PEG-3400-n76";  NKnameTeX(661) = "PEG-3400-n76";  
    NKname(662) = "PEG-2000-n44";  NKnameTeX(662) = "PEG-2000-n44";  
    NKname(663) = "PEG-200-n3";  NKnameTeX(663) = "PEG-200-n3";    NKsmiles(663) = "OCCOCCOCCOCCO";  
    NKname(664) = "PEG-200-n4";  NKnameTeX(664) = "PEG-200-n4";    NKsmiles(664) = "OCCOCCOCCOCCOCCO";   
    NKname(665) = "PEG-8000-n180";  NKnameTeX(665) = "PEG-8000-n180";  
    NKname(666) = "PEG-20000-n453";  NKnameTeX(666) = "PEG-20000-n453";  
    NKname(667) = "PEG-6000-n135";  NKnameTeX(667) = "PEG-6000-n135";  
    NKname(668) = "PEG-3350-n75";  NKnameTeX(668) = "PEG-3350-n75";  
    NKname(669) = "PEG-4000-n89";  NKnameTeX(669) = "PEG-4000-n89";  
    NKname(670) = "PEG-10000-n226";  NKnameTeX(670) = "PEG-10000-n226";  
    NKname(671) = "PEG-400-n7";  NKnameTeX(671) = "PEG-400-n7";    NKsmiles(671) = "OCCOCCOCCOCCOCCOCCOCCOCCO";  !(PEG-400 as a mixture of polyethylene glycols with n = number of oxyethylene groups. 1/3 n=7 and 2/3 n=8 for the weight of approx. 400 g/mol) 
    NKname(672) = "PEG-400-n8";  NKnameTeX(672) = "PEG-400-n8";    NKsmiles(672) = "OCCOCCOCCOCCOCCOCCOCCOCCOCCO";
    NKname(673) = "PEG-1000-n21";  NKnameTeX(673) = "PEG-1000-n21";    NKsmiles(673) = "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO";  
    NKname(674) = "PEG-1000-n22";  NKnameTeX(674) = "PEG-1000-n22";    NKsmiles(674) = "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO";  
    NKname(675) = "PEG-600-n12";  NKnameTeX(675) = "PEG-600-n12";    NKsmiles(675) = "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO";  
    NKname(676) = "PEG-600-n13";  NKnameTeX(676) = "PEG-600-n13";    NKsmiles(676) = "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO";  
    NKname(677) = "PEG-1450-n31";  NKnameTeX(677) = "PEG-1450-n31";  NKsmiles(677) = "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO";
    NKname(678) = "PEG-1450-n32";  NKnameTeX(678) = "PEG-1450-n32";  NKsmiles(678) = "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO";
    NKname(679) = "PEG-1540-n33";  NKnameTeX(679) = "PEG-1540-n33";  NKsmiles(679) = "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO";
    NKname(680) = "PEG-1540-n34";  NKnameTeX(680) = "PEG-1540-n34";  NKsmiles(680) = "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO";

    !...
    NKname(681) = "C18H29(OH)9";  NKnameTeX(681) = "C$_{18}$H$_{29}$(OH)$_9$";    ![C18H29(OH)9] 1,3,5,7,9,11,13,15,17-octadecanol "super-polyol"
    NKname(682) = "2,2,6,6-Tetrakis(hydroxymethyl)cyclohexanol";  NKnameTeX(682) = "2,2,6,6-Tetrakis(hydroxymethyl)cyclohexanol";    !(=2,2,6,6-tetrakis(hydroxylmethyl)cyclohexanol, C10)
    NKname(683) = "D-Glucuronic_acid";  NKnameTeX(683) = "D-Glucuronic_acid";  NKsmiles(683) = "C1C(C(C(C(C(O1)O)O)O)O)C(=O)O";
    NKname(684) = "Triton_X-100";  NKnameTeX(684) = "Triton_X-100";   NKsmiles(684) = "CC(C)(C)CC(C)(C)C1=CC=C(C=C1)OCCOCCOCCOCCOCCOCCOCCOCCOCCO";  
    NKname(685) = "Tween20";  NKnameTeX(685) = "Tween20";   NKsmiles(685) = "CCCCCCCCCCCC(=O)OC1CC(O)C(OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO)OC1COCCOCCOCCOCCOCCOCCOCCOCCOCCO";  
    !...
    !Polcyclic aromatic hydrocarbons (PAH)
    NKname(701) = "Benzene";                        NKnameTeX(701) = "Benzene";                        NKsmiles(701) = "c1ccccc1";
    NKname(702) = "Naphthalene";                    NKnameTeX(702) = "Naphthalene";                    NKsmiles(702) = "c1ccc2ccccc2c1";
    NKname(703) = "Anthracene";                     NKnameTeX(703) = "Anthracene";                     NKsmiles(703) = "c1ccc2cc3ccccc3cc2c1";
    NKname(704) = "Phenanthrene";                   NKnameTeX(704) = "Phenanthrene";                   NKsmiles(704) = "c1ccc2c(c1)ccc1ccccc12";
    NKname(705) = "Fluoranthene";                   NKnameTeX(705) = "Fluoranthene";                   NKsmiles(705) = "c1ccc-2c(c1)-c3cccc4c3c2ccc4";
    NKname(706) = "Pyrene";                         NKnameTeX(706) = "Pyrene";                         NKsmiles(706) = "c1cc2cccc3c2c4c1cccc4cc3";
    NKname(707) = "Chrysene";                       NKnameTeX(707) = "Chrysene";                       NKsmiles(707) = "c1ccc2c(c1)ccc1c3ccccc3ccc21";
    NKname(708) = "Perylene";                       NKnameTeX(708) = "Perylene";                       NKsmiles(708) = "c1ccc5cccc4c5c1c2cccc3cccc4c23";
    NKname(709) = "Coronene";                       NKnameTeX(709) = "Coronene";                       NKsmiles(709) = "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67";
    NKname(710) = "Toluene";                        NKnameTeX(710) = "Toluene";                        NKsmiles(710) = "Cc1ccccc1";
    NKname(711) = "2-Methylnaphthalene";            NKnameTeX(711) = "2-Methylnaphthalene";            NKsmiles(711) = "Cc1ccc2ccccc2c1";
    NKname(712) = "7-Isopropyl-1-methylphenanthrene"; NKnameTeX(712) = "7-Isopropyl-1-methylphenanthrene"; NKsmiles(712) = "Cc1cccc3c1ccc2cc(C(C)C)ccc23";   != Retene
    NKname(713) = "Benz[a]anthracene";              NKnameTeX(713) = "Benz[a]anthracene";              NKsmiles(713) = "c1ccc2c(c1)ccc3c2cc4ccccc4c3";
    NKname(714) = "Benzo[b]fluoranthene";           NKnameTeX(714) = "Benzo[b]fluoranthene";           NKsmiles(714) = "c5c3c1cccc2c1c(ccc2)c3c4ccccc4c5";
    NKname(715) = "Indeno[1,2,3-cd]pyrene";         NKnameTeX(715) = "Indeno[1,2,3-cd]pyrene";         NKsmiles(715) = "c1ccc2c(c1)c1c3ccc4ccccc4c3c3ccccc23";
    NKname(716) = "Benzo[a]pyrene";                 NKnameTeX(716) = "Benzo[a]pyrene";                 NKsmiles(716) = "c1ccc2c(c1)cc3ccc4cccc5c4c3c2cc5";
    
    NKname(717) = "6,12-dihydroxy-BaP";             NKnameTeX(717) = "6,12-dihydroxy-BaP";             NKsmiles(717) = "";
    NKname(718) = "BaP-6,12-dione";                 NKnameTeX(718) = "BaP-6,12-dione";                 NKsmiles(718) = "";
    NKname(719) = "BaP-derived_diCOOH";             NKnameTeX(719) = "BaP-derived_diCOOH";             NKsmiles(719) = "";   !7-oxo-7H-benzo[de]anthracene-3,4-dicarboxylic acid
    
    
    !...
    !Multifunctional aromatic hydrocarbons, phenolic compounds 
    NKname(720) = "Phenol";                    NKnameTeX(720) = "Phenol";                    NKsmiles(720) = "Oc1ccccc1";
    NKname(721) = "Protocatechuic_acid";       NKnameTeX(721) = "Protocatechuic_acid";       NKsmiles(721) = "O=C(O)c1ccc(O)c(O)c1";
    NKname(722) = "Vanillin";                  NKnameTeX(722) = "Vanillin";                  NKsmiles(722) = "O=Cc1ccc(O)c(OC)c1";      !4-Hydroxy-3-methoxybenzaldehyde
    NKname(723) = "Vanillic_acid";             NKnameTeX(723) = "Vanillic_acid";             NKsmiles(723) = "O=C(O)c1ccc(O)c(OC)c1";
    NKname(724) = "Gallic_acid";               NKnameTeX(724) = "Gallic_acid";               NKsmiles(724) = "O=C(O)c1cc(O)c(O)c(O)c1";
    NKname(725) = "Ferulic_acid";              NKnameTeX(725) = "Ferulic_acid";              NKsmiles(725) = "O=C(O)C=Cc1ccc(O)c(OC)c1";
    NKname(726) = "Syringic_acid";             NKnameTeX(726) = "Syringic_acid";             NKsmiles(726) = "O=C(O)c1cc(OC)c(O)c(OC)c1";
    NKname(727) = "2-Hydroxybenzoic_acid";     NKnameTeX(727) = "2-Hydroxybenzoic_acid";     NKsmiles(727) = "O=C(O)c1ccccc1O";         !(= Salicylic acid)
    NKname(728) = "3-Hydroxybenzoic_acid";     NKnameTeX(728) = "3-Hydroxybenzoic_acid";     NKsmiles(728) = "O=C(O)c1cccc(O)c1";
    NKname(729) = "4-Hydroxybenzoic_acid";     NKnameTeX(729) = "4-Hydroxybenzoic_acid";     NKsmiles(729) = "O=C(O)c1ccc(O)cc1";
    NKname(730) = "Phthalic_acid";             NKnameTeX(730) = "Phthalic_acid";             NKsmiles(730) = "O=C(O)c1ccccc1C(=O)O";    ! (= benzene-1,2-dicarboxylic acid)
    NKname(731) = "2,4-Dihydroxybenzaldehyde"; NKnameTeX(731) = "2,4-Dihydroxybenzaldehyde"; NKsmiles(731) = "O=Cc1ccc(O)cc1O";
    NKname(732) = "Vanillylmandelic_acid";     NKnameTeX(732) = "Vanillylmandelic_acid";     NKsmiles(732) = "COc1cc(CC(O)C(=O)O)ccc1O";
    NKname(733) = "3,5-Dihydroxybenzoic_acid"; NKnameTeX(733) = "3,5-Dihydroxybenzoic_acid"; NKsmiles(733) = "O=C(O)c1cc(O)cc(O)c1";
    NKname(734) = "Mandelic_acid";             NKnameTeX(734) = "Mandelic_acid";             NKsmiles(734) = "O=C(O)C(O)c1ccccc1";
    NKname(735) = "Dimethyl_phthalate";        NKnameTeX(735) = "Dimethyl_phthalate";        NKsmiles(735) = "COC(=O)c1ccccc1C(=O)OC";
    NKname(736) = "2,5-Dihydroxybenzoic_acid"; NKnameTeX(736) = "2,5-Dihydroxybenzoic_acid"; NKsmiles(736) = "O=C(O)c1cc(O)ccc1O";
    NKname(737) = "Resorcinol";                NKnameTeX(737) = "Resorcinol";                NKsmiles(737) = "Oc1cccc(O)c1";            !(= 1,3-benzenediol, = m-benzenediol)
    NKname(738) = "p-Cresol";                  NKnameTeX(738) = "p-Cresol";                  NKsmiles(738) = "Cc1ccc(O)cc1";
    !m-Cresol and o-Cresol are listed under 895, 896
    
    NKname(739) = "p-Benzenediol";              NKnameTeX(739) = "p-Benzenediol";              NKsmiles(739) = "Oc1ccc(O)cc1";                    !(= 1,4-benzenediol)
    NKname(740) = "4-Methylguaiacol";           NKnameTeX(740) = "4-Methylguaiacol";           NKsmiles(740) = "COc1cc(C)ccc1O";
    NKname(741) = "4-Propylguaiacol";           NKnameTeX(741) = "4-Propylguaiacol";           NKsmiles(741) = "CCCc1ccc(O)c(OC)c1";
    NKname(742) = "Coniferaldehyde";            NKnameTeX(742) = "Coniferaldehyde";            NKsmiles(742) = "COc1cc(C=CC=O)ccc1O";
    NKname(743) = "4-Methylsyringol";           NKnameTeX(743) = "4-Methylsyringol";           NKsmiles(743) = "COc1cc(C)c(O)c(OC)c1";
    NKname(744) = "Syringyl_acetone";           NKnameTeX(744) = "Syringyl_acetone";           NKsmiles(744) = "CC(=O)CCCc1cc(OC)c(O)c(OC)c1";      !(= 4-(4-hydroxy-3,5-dimethoxyphenyl)butan-2-one)
    NKname(745) = "p-Tolualdehyde";             NKnameTeX(745) = "p-Tolualdehyde";             NKsmiles(745) = "Cc1ccc(C=O)cc1";   !4-Methylbenzaldehyde
    NKname(746) = "2,5-Dimethylbenzaldehyde";   NKnameTeX(746) = "2,5-Dimethylbenzaldehyde";   NKsmiles(746) = "Cc1cc(C)ccc1C=O";
    NKname(747) = "Indan-1-one";                NKnameTeX(747) = "Indan-1-one";                NKsmiles(747) = "O=C1Cc2ccccc2C1";
    NKname(748) = "1H-phenalen-1-one";          NKnameTeX(748) = "1H-phenalen-1-one";          NKsmiles(748) = "O=c1c2cccc3cccc(c23)cc1";
    NKname(749) = "3,4-Dimethoxytoluene";       NKnameTeX(749) = "3,4-Dimethoxytoluene";       NKsmiles(749) = "Cc1cc(OC)c(OC)cc1";
    NKname(750) = "Veratric_acid";              NKnameTeX(750) = "Veratric_acid";              NKsmiles(750) = "COc1cc(C(=O)O)ccc1OC";
    NKname(751) = "Benzoic_acid";               NKnameTeX(751) = "Benzoic_acid";               NKsmiles(751) = "O=C(O)c1ccccc1";
    NKname(752) = "2-Naphthol";                 NKnameTeX(752) = "2-Naphthol";                 NKsmiles(752) = "Oc1ccc2ccccc2c1";
    NKname(753) = "1,4-Naphthalenedione";       NKnameTeX(753) = "1,4-Naphthalenedione";       NKsmiles(753) = "O=C1C=CC(=O)c2ccccc12";
    NKname(754) = "5-Hydroxy-1,4-naphthalenedione"; NKnameTeX(754) = "5-Hydroxy-1,4-naphthalenedione"; NKsmiles(754) = "O=C1C=CC(=O)c2cc(O)cc12";
    NKname(755) = "2-Carboxycinnamic_acid";     NKnameTeX(755) = "2-Carboxycinnamic_acid";     NKsmiles(755) = "O=C(O)c1ccccc1C=CC(=O)O";
    NKname(756) = "4-Hydroxyphthalic_acid";     NKnameTeX(756) = "4-Hydroxyphthalic_acid";     NKsmiles(756) = "O=C(O)c1cc(O)cc(C(=O)O)c1";
    NKname(757) = "Carminic_acid";              NKnameTeX(757) = "Carminic_acid";              NKsmiles(757) = "O=C(O)c2c(c3C(=O)c1c(O)c(c(O)c(O)c1C(=O)c3cc2O)[C@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)CO)C";   !complex anthraquinone glycoside; recommend database-derived structure
    NKname(758) = "Triolein";                   NKnameTeX(758) = "Triolein";                   NKsmiles(758) = "CCCCCCCCC=CCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCC=CCCCCCCCC)OC(=O)CCCCCCCCC=CCCCCCCCC";
    NKname(759) = "Linoleic_acid";              NKnameTeX(759) = "Linoleic_acid";              NKsmiles(759) = "CCCCCC=CCC=CCCCCCCCC(=O)O";
        
    !other complex functionalyzed compounds, MCM names for additions to the toluene system by NG
    NKname(760) = "C3DIALOOH";  NKnameTeX(760) = "C3DIALOOH";  
    NKname(761) = "C33CO";      NKnameTeX(761) = "C33CO";
    NKname(762) = "C23O3CCHO";  NKnameTeX(762) = "C23O3CCHO";
    NKname(763) = "C535OOH";    NKnameTeX(763) = "C535OOH";
    NKname(764) = "C534OOH";    NKnameTeX(764) = "C534OOH";
    !...
    !Carbonyls (Ketones and Aldehydes)
    NKname(802) = "Acetone";        NKnameTeX(802) = "Acetone";         NKsmiles(802) = "CC(=O)C";
    NKname(803) = "2-Butanone";     NKnameTeX(803) = "2-Butanone";      NKsmiles(803) = "CCC(=O)C";
    NKname(804) = "3-Methyl-2-butanone";    NKnameTeX(804) = "3-Methyl-2-butanone";     NKsmiles(804) = "CC(C)C(=O)C"; !Methyl Isopropyl Ketone
    NKname(805) = "4-Methyl-2-pentanone";   NKnameTeX(805) = "4-Methyl-2-pentanone";    NKsmiles(805) = "CC(C)CC(=O)C";
    NKname(806) = "2-Pentanone";    NKnameTeX(806) = "2-Pentanone";     NKsmiles(806) = "CCCC(=O)C";
    NKname(807) = "2-Hexanone";     NKnameTeX(807) = "2-Hexanone";      NKsmiles(807) = "CCCCC(=O)C";
    NKname(808) = "2-Heptanone";    NKnameTeX(808) = "2-Heptanone";     NKsmiles(808) = "CCCCCC(=O)C";
    NKname(809) = "3-Heptanone";    NKnameTeX(809) = "3-Heptanone";     NKsmiles(809) = "CCCCC(=O)CC";
    NKname(810) = "3-Pentanone";    NKnameTeX(810) = "3-Pentanone";     NKsmiles(810) = "CCC(=O)CC";
    NKname(811) = "2-Octanone";     NKnameTeX(811) = "2-Octanone";      NKsmiles(811) = "CC(=O)CCCCCC";
    NKname(812) = "4-Heptanone";    NKnameTeX(812) = "4-Heptanone";     NKsmiles(812) = "CCCC(=O)CCC";
    NKname(813) = "Acetaldehyde";   NKnameTeX(813) = "Acetaldehyde";    NKsmiles(813) = "CC=O";
    NKname(814) = "Propionaldehyde";NKnameTeX(814) = "Propionaldehyde"; NKsmiles(814) = "CCC=O";
    NKname(815) = "Butyraldehyde";  NKnameTeX(815) = "Butyraldehyde";   NKsmiles(815) = "CCCC=O";
    NKname(816) = "2-Nonanone";     NKnameTeX(816) = "2-Nonanone";      NKsmiles(816) = "CC(=O)CCCCCCC";
    NKname(817) = "Dodecanone";     NKnameTeX(817) = "Dodecanone";      NKsmiles(817) = "CC(=O)CCCCCCCCCC";
    NKname(818) = "Tridecanone";    NKnameTeX(818) = "Tridecanone";     NKsmiles(818) = "CC(=O)CCCCCCCCCCC";
    NKname(819) = "Tetradecanone";  NKnameTeX(819) = "Tetradecanone";   NKsmiles(819) = "CC(=O)CCCCCCCCCCCC";
    NKname(820) = "Pentadecanone";  NKnameTeX(820) = "Pentadecanone";   NKsmiles(820) = "CC(=O)CCCCCCCCCCCCC";
    NKname(821) = "Hexadecanone";   NKnameTeX(821) = "Hexadecanone";    NKsmiles(821) = "CC(=O)CCCCCCCCCCCCCC";
    NKname(822) = "Heptadecanone";  NKnameTeX(822) = "Heptadecanone";   NKsmiles(822) = "CC(=O)CCCCCCCCCCCCCCC";
    NKname(823) = "Octadecanone";   NKnameTeX(823) = "Octadecanone";    NKsmiles(823) = "CC(=O)CCCCCCCCCCCCCCCC";
    NKname(824) = "Nonadecanone";   NKnameTeX(824) = "Nonadecanone";    NKsmiles(824) = "CC(=O)CCCCCCCCCCCCCCCCC";
    NKname(825) = "Icosanone";      NKnameTeX(825) = "Icosanone";       NKsmiles(825) = "CC(=O)CCCCCCCCCCCCCCCCCC";
    NKname(826) = "Henicosanone";   NKnameTeX(826) = "Henicosanone";    NKsmiles(826) = "CC(=O)CCCCCCCCCCCCCCCCCCC";
    NKname(827) = "Docosanone";     NKnameTeX(827) = "Docosanone";      NKsmiles(827) = "CC(=O)CCCCCCCCCCCCCCCCCCCC";
    NKname(828) = "Tricosanone";    NKnameTeX(828) = "Tricosanone";     NKsmiles(828) = "CC(=O)CCCCCCCCCCCCCCCCCCCCC";
    NKname(829) = "Tetracosanone";  NKnameTeX(829) = "Tetracosanone";   NKsmiles(829) = "CC(=O)CCCCCCCCCCCCCCCCCCCCCC";
    NKname(830) = "Pentacosanone";  NKnameTeX(830) = "Pentacosanone";   NKsmiles(830) = "CC(=O)CCCCCCCCCCCCCCCCCCCCCCC";
    NKname(831) = "Undecanone";     NKnameTeX(831) = "Undecanone";      NKsmiles(831) = "CC(=O)CCCCCCCCC";

    !...
    !Esters
    NKname(837) = "Cholesteryl_oleate";        NKnameTeX(837) = NKname(837);    NKsmiles(837) = "CCCCCCCCC=CCCCCCCCC(=O)OC1CCC2C3CCC4=CC(CC[C@]4(C)C3CC[C@]12C)C";
    NKname(838) = "Triacontyl_palmitate";      NKnameTeX(838) = NKname(838);    NKsmiles(838) = "CCCCCCCCCCCCCCCC(=O)OCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    NKname(839) = "1-Palmitoyl-2-oleoyl-3-linoleoyl_glycerol"; NKnameTeX(839) = NKname(839);    NKsmiles(839) = "CCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCC=CCCCCCCCC)OC(=O)CCCCCC=CCC=CCCCCCCCC";
    NKname(840) = "bis(2-ethylhexyl)_sebacate"; NKnameTeX(840) = "bis(2-ethylhexyl)_sebacate";  NKsmiles(840) = "CCCCC(CC)COC(=O)CCCCCCCCCC(=O)OCC(CC)CCCCC";
    NKname(841) = "Methyl_acetate";             NKnameTeX(841) = "Methyl_acetate";          NKsmiles(841) = "CC(=O)OC";
    NKname(842) = "Ethyl_acetate";              NKnameTeX(842) = "Ethyl_acetate";           NKsmiles(842) = "CC(=O)OCC";
    NKname(843) = "1-Propyl_acetate";           NKnameTeX(843) = "1-Propyl_acetate";        NKsmiles(843) = "CC(=O)OCCC";
    NKname(844) = "1-Butyl_acetate";            NKnameTeX(844) = "1-Butyl_acetate";         NKsmiles(844) = "CC(=O)OCCCC";
    NKname(845) = "Isobutyl_acetate";           NKnameTeX(845) = "Isobutyl_acetate";        NKsmiles(845) = "CC(=O)OCC(C)C";
    NKname(846) = "2-Butyl_acetate";            NKnameTeX(846) = "2-Butyl_acetate";         NKsmiles(846) = "CC(=O)OC(C)CC";
    NKname(847) = "tert-Butyl_acetate";         NKnameTeX(847) = "tert-Butyl_acetate";      NKsmiles(847) = "CC(=O)OC(C)(C)C";
    NKname(848) = "1-Pentyl_acetate";           NKnameTeX(848) = "1-Pentyl_acetate";        NKsmiles(848) = "CC(=O)OCCCCC";
    NKname(849) = "1-Hexyl_acetate";            NKnameTeX(849) = "1-Hexyl_acetate";         NKsmiles(849) = "CC(=O)OCCCCCC";
    NKname(850) = "2-Ethoxyethyl_acetate";      NKnameTeX(850) = "2-Ethoxyethyl_acetate";   NKsmiles(850) = "CC(=O)OCCOCC";
    NKname(851) = "Octadecanoic_acid_methyl_ester";   NKnameTeX(851) = "Octadecanoic_acid_methyl_ester";    NKsmiles(851) = "CCCCCCCCCCCCCCCCCC(=O)OC";     !(= methyl stearate)
    NKname(852) = "Octadecanoic_acid_ethyl_ester";    NKnameTeX(852) = "Octadecanoic_acid_ethyl_ester";     NKsmiles(852) = "CCCCCCCCCCCCCCCCCC(=O)OCC";    !(= ethyl stearate)

    !Diketones
    NKname(853) = "Undecanedione";      NKnameTeX(853) = "Undecanedione";  
    NKname(854) = "Dodecanedione";      NKnameTeX(854) = "Dodecanedione";  
    NKname(855) = "Tridecanedione";     NKnameTeX(855) = "Tridecanedione";  
    NKname(856) = "Tetradecanedione";   NKnameTeX(856) = "Tetradecanedione";  
    NKname(857) = "Pentadecanedione";   NKnameTeX(857) = "Pentadecanedione";  
    NKname(858) = "Hexadecanedione";    NKnameTeX(858) = "Hexadecanedione";  
    NKname(859) = "Heptadecanedione";   NKnameTeX(859) = "Heptadecanedione";  
    NKname(860) = "Octadecanedione";    NKnameTeX(860) = "Octadecanedione";  
    NKname(861) = "Nonadecanedione";    NKnameTeX(861) = "Nonadecanedione";  
    NKname(862) = "Icosanedione";       NKnameTeX(862) = "Icosanedione";  
    NKname(863) = "Henicosanedione";    NKnameTeX(863) = "Henicosanedione";  
    NKname(864) = "Docosanedione";      NKnameTeX(864) = "Docosanedione";  
    NKname(865) = "Tricosanedione";     NKnameTeX(865) = "Tricosanedione";  

    !...
    !Ethers
    NKname(881) = "2-Methoxy-2-methylpropane";   NKnameTeX(881) = "2-Methoxy-2-methylpropane";   NKsmiles(881) = "CC(C)(C)OC";
    NKname(882) = "2-Methoxyethanol";            NKnameTeX(882) = "2-Methoxyethanol";            NKsmiles(882) = "COCCO";
    NKname(883) = "2-Ethoxyethanol";             NKnameTeX(883) = "2-Ethoxyethanol";             NKsmiles(883) = "CCOCCO";
    NKname(884) = "1-Methoxy-2-propanol";        NKnameTeX(884) = "1-Methoxy-2-propanol";        NKsmiles(884) = "CC(O)COC";
    NKname(885) = "2-Isopropoxyethanol";         NKnameTeX(885) = "2-Isopropoxyethanol";         NKsmiles(885) = "CC(C)OCCO";
    NKname(886) = "2-Butoxyethanol";             NKnameTeX(886) = "2-Butoxyethanol";             NKsmiles(886) = "CCCCOCCO";
    NKname(887) = "2-Methoxypropanol";           NKnameTeX(887) = "2-Methoxypropanol";           NKsmiles(887) = "CC(CO)OC";
    NKname(888) = "1-(2-methoxypropoxy)-2-propanol"; NKnameTeX(888) = "1-(2-methoxypropoxy)-2-propanol"; NKsmiles(888) = "CC(O)COC(C)COC";
    NKname(889) = "2-(2-methoxyethoxy)ethanol";  NKnameTeX(889) = "2-(2-methoxyethoxy)ethanol";  NKsmiles(889) = "COCCOCCO";
    NKname(890) = "2-(2-ethoxyethoxy)ethanol";   NKnameTeX(890) = "2-(2-ethoxyethoxy)ethanol";   NKsmiles(890) = "CCOCCOCCO";
    NKname(891) = "1,4-Dioxane";                 NKnameTeX(891) = "1,4-Dioxane";                 NKsmiles(891) = "C1COCCO1";
    NKname(892) = "Tetrahydrofuran";             NKnameTeX(892) = "Tetrahydrofuran";             NKsmiles(892) = "C1CCOC1";
    NKname(893) = "2-Ethoxy-2-methylpropane";    NKnameTeX(893) = "2-Ethoxy-2-methylpropane";    NKsmiles(893) = "CCOC(C)(C)C";
    NKname(894) = "Diethyl_ether";               NKnameTeX(894) = "Diethyl_ether";               NKsmiles(894) = "CCOCC"; 
    !...
    !other complex functionalyzed compounds, MCM names:
    NKname(895) = "m-Cresol";       NKnameTeX(895) = "m-Cresol";       NKsmiles(895) = "Cc1cccc(O)c1";
    NKname(896) = "o-Cresol";       NKnameTeX(896) = "o-Cresol";       NKsmiles(896) = "Cc1ccccc1O";
    NKname(897) = "Anisaldehyde";   NKnameTeX(897) = "Anisaldehyde";   NKsmiles(897) = "COc1ccc(C=O)cc1";   ! usually p-Anisaldehyde (4-methoxybenzaldehyde)
    NKname(898) = "Trehalose";      NKnameTeX(898) = "Trehalose";      NKsmiles(898) = "OC1C(O)C(O)C(OC2OC(CO)C(O)C(O)C2O)OC1CO";
    NKname(899) = "Maltose";        NKnameTeX(899) = "Maltose";        NKsmiles(899) = "OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O";

    NKname(901) = "alpha-Pinene";  NKnameTeX(901) = "$\alpha$-Pinene";     !alpha-Pinene = 2,6,6-trimethylbicyclo[3.1.1]hept-2-ene
    NKname(902) = "PINONIC";  NKnameTeX(902) = "PINONIC";  
    NKname(903) = "PINAL";  NKnameTeX(903) = "PINAL";  
    NKname(904) = "PINIC";  NKnameTeX(904) = "PINIC";  
    NKname(905) = "Norpinonic_acid";  NKnameTeX(905) = "Norpinonic_acid";  
    NKname(906) = "NORPINIC";  NKnameTeX(906) = "NORPINIC";  
    NKname(907) = "C89CO2H";  NKnameTeX(907) = "C89CO2H";    !Pinalic_acid
    NKname(908) = "HOPINONIC";  NKnameTeX(908) = "HOPINONIC";  
    NKname(909) = "C721CHO";  NKnameTeX(909) = "C721CHO";   !Norpinalic_acid
    NKname(910) = "8-Hydroxypinalic_acid";  NKnameTeX(910) = "8-Hydroxypinalic_acid";  
    NKname(911) = "C108OOH";  NKnameTeX(911) = "C108OOH";   !C10H16O5 (a hydroperoxyketoaldehyde)
    NKname(912) = "C97OOH";  NKnameTeX(912) = "C97OOH";   !C97OOH (MCM name), C9H16O4
    NKname(913) = "ESTER_dimer";  NKnameTeX(913) = "ESTER_dimer";   !C19H28O7 ketoesterdiacid, dimer from ester formation
    NKname(914) = "ALDOL_dimer";  NKnameTeX(914) = "ALDOL_dimer";   !C19H28O7 hydroperoxyketoaldehydeacid, dimer from Aldol condensation
    NKname(915) = "C813OOH";  NKnameTeX(915) = "C813OOH";   !C813OOH (MCM name);  C8H14O6 
    NKname(916) = "C107OOH";  NKnameTeX(916) = "C107OOH";   
    NKname(917) = "APINBOH";  NKnameTeX(917) = "APINBOH";  
    NKname(918) = "C107OH";  NKnameTeX(918) = "C107OH";  
    NKname(919) = "APINAOOH";  NKnameTeX(919) = "APINAOOH";   
    NKname(920) = "C108OH";  NKnameTeX(920) = "C108OH";   
    NKname(921) = "C98OOH";  NKnameTeX(921) = "C98OOH";   
    NKname(922) = "APINCOOH";  NKnameTeX(922) = "APINCOOH";   
    NKname(923) = "C921OOH";  NKnameTeX(923) = "C921OOH";   
    NKname(924) = "C97OH";  NKnameTeX(924) = "C97OH";   
    NKname(925) = "C812OH";  NKnameTeX(925) = "C812OH";   
    NKname(926) = "APINBCO";  NKnameTeX(926) = "APINBCO";   
    NKname(927) = "C811OH";  NKnameTeX(927) = "C811OH";  
    NKname(928) = "PINALOOH";  NKnameTeX(928) = "PINALOOH";  
    NKname(929) = "C109OOH";  NKnameTeX(929) = "C109OOH";  
    NKname(930) = "C812OOH";  NKnameTeX(930) = "C812OOH";   
    NKname(931) = "C109OH";  NKnameTeX(931) = "C109OH";   
    NKname(932) = "Pinolic_acid";  NKnameTeX(932) = "Pinolic_acid";   
    NKname(933) = "2-Hydroxypinane-3-nitrate";  NKnameTeX(933) = "2-Hydroxypinane-3-nitrate";   
    NKname(934) = "KETALDOOH";  NKnameTeX(934) = "KETALDOOH";   
    NKname(935) = "KETALDOH";  NKnameTeX(935) = "KETALDOH";  
    NKname(936) = "DECOMP6";  NKnameTeX(936) = "DECOMP6";  
    NKname(937) = "DECOMP1";  NKnameTeX(937) = "DECOMP1";  
    NKname(938) = "DECOMP2";  NKnameTeX(938) = "DECOMP2";  
    NKname(939) = "DECOMP3";  NKnameTeX(939) = "DECOMP3";   != 5-methyl-2(5H)-furanone
    NKname(940) = "DECOMP4";  NKnameTeX(940) = "DECOMP4";   != 2(5H)-furanone
    NKname(941) = "DECOMP5";  NKnameTeX(941) = "DECOMP5";  
    NKname(942) = "RBCOHOHOOH";  NKnameTeX(942) = "RBCOHOHOOH";  
    NKname(943) = "ROHOOH";  NKnameTeX(943) = "ROHOOH";  
    NKname(944) = "REPOX";  NKnameTeX(944) = "REPOX";  
    NKname(945) = "RBCOHOOH";  NKnameTeX(945) = "RBCOHOOH";  
    NKname(946) = "ROHOH";  NKnameTeX(946) = "ROHOH";  
    NKname(947) = "ROOH";  NKnameTeX(947) = "ROOH";  
    NKname(948) = "CARB";  NKnameTeX(948) = "CARB";  
    NKname(949) = "OHROOH";  NKnameTeX(949) = "OHROOH";  
    NKname(950) = "CARBROOH";  NKnameTeX(950) = "CARBROOH";  
    NKname(951) = "pCARBROOHplusOH";  NKnameTeX(951) = "pCARBROOHplusOH";  
    NKname(952) = "C11ALDEHYDE";  NKnameTeX(952) = "C11ALDEHYDE";  
    NKname(953) = "pCARBROOHplushv";  NKnameTeX(953) = "pCARBROOHplushv";  
    NKname(954) = "OHCARB";  NKnameTeX(954) = "OHCARB";  
    NKname(955) = "PHA_dimer";  NKnameTeX(955) = "PHA_dimer";  
    NKname(956) = "Hydroxymethylfurfural";  NKnameTeX(956) = "Hydroxymethylfurfural";  
    NKname(957) = "3-methyl-1,2,3-butane-tricarboxylic_acid";  NKnameTeX(957) = "3-methyl-1,2,3-butane-tricarboxylic_acid";  
    NKname(958) = "Hydroperoxy_dimethoxy_diol";  NKnameTeX(958) = "Hydroperoxy_dimethoxy_diol";   !complete name: 4-hydroperoxy-2,6-dimethoxycyclohexa-1,5-diene-1,3-diol
    
    !.... Ampritta's isoprene surrogate system ................
    NKName(959) = "Glyoxal"; NKnameTeX(959) = "Glyoxal";
    NKName(960) = "Methylglyoxal"; NKnameTeX(960) = "Methylglyoxal";
    NKName(961) = "Methyl_vinyl_ketone"; NKnameTeX(961) = "Methyl_vinyl_ketone";
    NKName(962) = "MACO3H"; NKnameTeX(962) = "MACO3H";
    NKName(963) = "MACROOH"; NKnameTeX(963) = "MACROOH";
    NKName(964) = "HMACROOH"; NKnameTeX(964) = "HMACROOH";
    NKName(965) = "2-Methyltetrol"; NKnameTeX(965) = "2-Methyltetrol";
    NKName(966) = "2-Methylglyceric_acid"; NKnameTeX(966) = "2-Methylglyceric_acid";
    NKName(967) = "2-Hydroxy-dihydroperoxide"; NKnameTeX(967) = "2-Hydroxy-dihydroperoxide";
    NKName(968) = "C5-alkene_triol"; NKnameTeX(968) = "C5-alkene_triol"; 
    NKName(969) = "C10-hemiacetal_dimer"; NKnameTeX(969) = "C10-hemiacetal_dimer";

    !!!add new organics for MT system used in project 2016-05-13 by JM
    !!NKname(960) = "C813OH";  NKnameTeX(960) = "C813OH";   !complete name: ?
    !!NKname(961) = "C922OOH";  NKnameTeX(961) = "C922OOH";   !complete name: ?
    !!NKname(962) = "C98OH";  NKnameTeX(962) = "C98OH";   !complete name: ?
    !!NKname(963) = "C813O2";  NKnameTeX(963) = "C813O2";   !complete name: ?
    !!NKname(964) = "C922O2";  NKnameTeX(964) = "C922O2";   !complete name: ?
    !!NKname(965) = "C811PAN";  NKnameTeX(965) = "C811PAN";   !complete name: ?
    !!NKname(966) = "C920PAN";  NKnameTeX(966) = "C920PAN";   !complete name: ?

    !add new organics for IS system used in project 2016-05-14 by JM
    NKname(970) = "IEB1OOH";  NKnameTeX(970) = "IEB1OOH";   !complete name: ?
    NKname(971) = "IEB2OOH";  NKnameTeX(971) = "IEB2OOH";   !complete name: ?
    NKname(972) = "C59OOH";  NKnameTeX(972) = "C59OOH";   !complete name: ?
    NKname(973) = "IEC1OOH";  NKnameTeX(973) = "IEC1OOH";   !complete name: ?
    NKname(974) = "C58OOH";  NKnameTeX(974) = "C58OOH";   !complete name: ?
    NKname(975) = "IEPOXA";  NKnameTeX(975) = "IEPOXA";   !complete name: ?
    NKname(976) = "C57OOH";  NKnameTeX(976) = "C57OOH";   !complete name: ?
    NKname(977) = "IEPOXC";  NKnameTeX(977) = "IEPOXC";   !complete name: ?
    NKname(978) = "HIEB1OOH";  NKnameTeX(978) = "HIEB1OOH";   !complete name: ?
    NKname(979) = "INDOOH";  NKnameTeX(979) = "INDOOH";   !complete name: ?
    NKname(980) = "IEACO3H";  NKnameTeX(980) = "IEACO3H";   !complete name: ?
    NKname(981) = "C525OOH";  NKnameTeX(981) = "C525OOH";   !complete name: ?
    NKname(982) = "HIEB2OOH";  NKnameTeX(982) = "HIEB2OOH";   !complete name: ?
    NKname(983) = "IEC2OOH";  NKnameTeX(983) = "IEC2OOH";   !complete name: ?
    NKname(984) = "INAOOH";  NKnameTeX(984) = "INAOOH";   !complete name: ?
    NKname(985) = "C510OOH";  NKnameTeX(985) = "C510OOH";   !complete name: ?
    NKname(986) = "INB1OOH";  NKnameTeX(986) = "INB1OOH";   !complete name: ?
    NKname(987) = "IECCO3H";  NKnameTeX(987) = "IECCO3H";   !complete name: ?
    NKname(988) = "INCOOH";  NKnameTeX(988) = "INCOOH";   !complete name: ?
    NKname(989) = "INB2OOH";  NKnameTeX(989) = "INB2OOH";   !complete name: ?
    NKname(990) = "2-Methyltetrol_dimer";  NKnameTeX(990) = "2-Methyltetrol_dimer";   !complete name: 2-methyl-4-((1,3,4-trihydroxy-2-methylbutan-2-yl)oxy)butane-1,2,3-triol
    NKname(991) = "MBTCA";  NKnameTeX(991) = "MBTCA";   !complete name: 3-methylbutane-1,2,3-tricarboxylic acid
    !compounds used in Havala Pye's SOAS mixture
    NKname(992) = "MO-OOA";  NKnameTeX(992) = "MO-OOA";   !complete name: ?
    NKname(993) = "BBOA";  NKnameTeX(993) = "BBOA";   !complete name: ?
    NKname(994) = "IEPOXOA";  NKnameTeX(994) = "IEPOXOA";   !complete name: ?
    NKname(995) = "LO-OOA";  NKnameTeX(995) = "LO-OOA";   !complete name: ?
    NKname(996) = "2-Methyltetrol";  NKnameTeX(996) = "2-Methyltetrol";    !complete name: ?
    NKname(997) = "C5-alkenetriol";  NKnameTeX(997) = "C5-alkenetriol";    !complete name: a C5-alkene_triol from isoprene oxidation
    NKname(998) = "MGA";  NKnameTeX(998) = "MGA";   !2,3-dihydroxy-2-methylpropanoic acid;   CC(C(O)=O)(CO)O
    NKname(999) = "Hydroxyglutaric_acid";  NKnameTeX(999) = "Hydroxyglutaric_acid";  
    
    NKname(1000) = "C15H18O7_syringol_aqSOA";  NKnameTeX(1000) = "C$_{15}$H$_{18}$O$_{7}$_syringol_aqSOA";   !complete name: 6-(2,4-dihydroxy-3,5-dimethoxyphenyl)-5-methoxyhexa-3,5-dienoic acid
    
    !NKname(1002) -- 1010) temporarily reserved for Antoine
    
    !organosulfates
    NKname(1011) = "Methyl_sulfate";  NKnameTeX(1011) = "Methyl_sulfate";  
    NKname(1012) = "Ethyl_sulfate";   NKnameTeX(1012) = "Ethyl_sulfate";  
    NKname(1013) = "Isoprene_OS_1";  NKnameTeX(1013) = "Isoprene_OS_1";  
    NKname(1014) = "Isoprene_OS_2";  NKnameTeX(1014) = "Isoprene_OS_2";  
    NKname(1015) = "Isoprene_OS_3";  NKnameTeX(1015) = "Isoprene_OS_3";  
    NKname(1016) = "Isoprene_OS_4";  NKnameTeX(1016) = "Isoprene_OS_4";

    NKname(1020) = "IEPOX";  NKnameTeX(1020) = "IEPOX"; NKsmiles(1020) = "OCC1OC1(C)CO"; !Isoprene-derived epoxydiol;
    
    !...
    NKname(1500) = "!fromInpFile!";  NKnameTeX(1500) = "!fromInpFile!";   !this ID number is used to indicate that the name is set from input via a file.

    end subroutine nametab
!==========================================================================================================
    
    

    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to set the character strings for the names of the different components  *
    !*   and ions in a mixture.                                                             * 
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend,                                                                        *
    !*   IACETH, ETH Zurich, (2004 - 2009)                                                  *
    !*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2006                                                            *
    !*   -> latest changes: 2026-09-16                                                      *
    !*                                                                                      *
    !****************************************************************************************
    subroutine names_mix(CompN, compname, compnameTeX, cpsmiles, ionname, ionnameTeX, OtoCratio, HtoCratio)

    use Mod_kind_param, only : wp
    use ModSystemProp, only : ElectComps, ElectNues, ElectSubs, nneutral, Ncation, Nanion, NGS, frominpfile, cpname
    use ModSubgroupProp, only : O2C_H2C_component, subgrname, subgrnameTeX
    use ModStringFunctions, only : replace_text, replace_text_advance

    implicit none

    !interface vars:
    integer,dimension(:),intent(in) :: CompN
    real(wp),dimension(:),intent(out) :: OtoCratio, HtoCratio
    character(len=*),dimension(:),intent(out) :: compname, compnameTeX
    character(len=*),dimension(:),intent(inout) :: cpsmiles
    character(len=*),dimension(:),intent(out) :: ionname, ionnameTeX
    !local vars:
    character(len=1) :: nchar
    character(len=16) :: catxt, antxt
    character(len=40) :: txt
    integer :: i, k, nn, cnt, cn, an, nue_c, nue_a
    real(wp) :: cpCarbon, cpHydrogen, cpOxygen, cpNitrogen, cpSulfur, NtoCratio, StoCratio
    !...........................................................

    OtoCratio = -7.777777_wp  !set to an impossible value at initialization
    HtoCratio = -7.777777_wp
    !Get the names of the actual components in mixture nd
    !neutral components:
    !water is always component 1 (if there is water in the mixture):
    if (CompN(1) == 401) then !there is water
        compname(1) = "Water"
        compnameTeX(1) = "Water"
        cpsmiles(1) = "O"
    endif

    do i = 1,nneutral
        k = CompN(i)
        if (k > 0 .and. k /= 401) then                  !scan all components, except water
            if (frominpfile) then
                compname(i) = cpname(i) 
                compnameTeX(i) = cpname(i) 
            else
                compname(i) = trim(NKname(k)) !//"'"    !the "'" is added to make sure that for comma separated reading the commas in e.g. 1,2-Butanediol is not separating the species name.
                compnameTeX(i) = trim(NKnameTeX(k)) !//"'"
            endif
            !call the OtoC calculation subroutine which uses atom information from subgroups;  here just for the single components "i".
            call O2C_H2C_component(i, cpCarbon, cpHydrogen, cpOxygen, cpNitrogen, cpSulfur, OtoCratio(i), HtoCratio(i), NtoCratio, StoCratio)
            
            !also assign SMILES to those components with a defined one:
            if (len_trim(cpsmiles(i)) < 1) then
                if (NKsmiles(k) /= "not_defined") then
                    cpsmiles(i) = trim(NKsmiles(k))
                endif
            endif
        endif
        
    enddo

    !construct electrolyte names from ion subgroup names:
    cnt = Ncation*Nanion
    do i = 1,cnt
        cn = ElectComps(i,1)
        an = ElectComps(i,2)  
        !remove the parathensis ( ) around the ion subgroup as well as the charge characters:
        catxt = trim( replace_text(trim(subgrname(cn)), '(', '') ) 
        catxt = trim( replace_text(catxt, ')', '') )
        catxt = trim( replace_text(catxt, '+', '') )
        antxt = trim( replace_text(trim(subgrname(an)), '(', '') ) 
        antxt = trim( replace_text(antxt, ')', '') )
        antxt = trim( replace_text(antxt, '-', '') )
        !construct 'neutral' electrolyte component name from ion info using
        !common naming convention from chemistry:
        nue_c = ElectNues(i,1)
        nue_a = ElectNues(i,2)
        if (nue_c /= nue_a) then
            if (nue_c > nue_a) then
                write(nchar,'(I0)') nue_c/nue_a
                k = len_trim(catxt)
                nn = iachar(catxt(k:k))
                if (nn > 49 .and. nn < 57) then     !the last character is the number
                    compname(nneutral+i) = '('//trim(catxt)//')'//nchar//trim(antxt)
                else
                    compname(nneutral+i) = trim(catxt)//nchar//trim(antxt)
                endif   
            else
                write(nchar,'(I0)') nue_a/nue_c
                k = len_trim(antxt)
                nn = iachar(antxt(k:k))
                if (nn > 49 .and. nn < 57) then     !the last character is the number
                    compname(nneutral+i) = trim(catxt)//'('//trim(antxt)//')'//nchar
                else
                    compname(nneutral+i) = trim(catxt)//trim(antxt)//nchar
                endif 
            endif
        else   
            compname(nneutral+i) = trim(catxt)//trim(antxt)        
        endif
        !generate the name also using TeX formatting:
        txt = trim(compname(nneutral+i))
        txt = trim( replace_text_advance(trim(txt), '2', '$_2$') ) 
        txt = trim( replace_text_advance(trim(txt), '3', '$_3$') ) 
        txt = trim( replace_text_advance(trim(txt), '4', '$_4$') ) 
        compnameTeX(nneutral+i) = trim(txt)
    enddo
    !set the rest of the names to an empty string
    compname(nneutral+cnt+1:) = ""
    compnameTeX(nneutral+cnt+1:) = ""

    !set also ion names, especially for output use with systems containing several electrolytes:
    !cations:
    cnt = 0
    do i = 1,NGS !loop over ElectSubs
        k = ElectSubs(i)
        if (k < 240) then !cation
            cnt = cnt +1
            !remove the parathensis (...) around the ion subgroup:
            txt = trim( replace_text(trim(subgrname(k)), '(', '') ) 
            txt = trim( replace_text(txt, ')', '') )
            ionname(cnt) = trim(txt)
            txt = trim( replace_text(trim(subgrnameTeX(k)), '(', '') ) 
            txt = trim( replace_text(txt, ')', '') )
            ionnameTeX(cnt) = trim(txt)
        endif
    enddo
    !anions:
    do i = 2,NGS !loop over ElectSubs
        k = ElectSubs(i)
        if (k > 240) then !anion
            cnt = cnt +1
            !remove the parathensis (...) around the ion subgroup:
            txt = trim( replace_text(trim(subgrname(k)), '(', '') ) 
            txt = trim( replace_text(txt, ')', '') )
            ionname(cnt) = trim(txt)
            txt = trim( replace_text(trim(subgrnameTeX(k)), '(', '') ) 
            txt = trim( replace_text(txt, ')', '') )
            ionnameTeX(cnt) = trim(txt)
        endif
    enddo

    end subroutine names_mix
!========================================================================================================== 
    
end module ModComponentNames
