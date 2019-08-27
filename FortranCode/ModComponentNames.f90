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
!*   -> latest changes: 2018/05/22                                                      *
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
!*   -  SUBROUTINE nametab                                                              *
!*   -  SUBROUTINE names_mix                                                            *
!*                                                                                      *
!****************************************************************************************

MODULE ModComponentNames

IMPLICIT NONE

!module public vars:
CHARACTER(LEN=60),DIMENSION(1500),PUBLIC :: NKname, NKnameTeX !neutral component names (alphabetical names and DISLIN TeX-code version)
CHARACTER(LEN=30),DIMENSION(40,40),PUBLIC :: electname, electnameTeX !ion combinations component names

!================================================================================================================================= 
    CONTAINS
!================================================================================================================================= 
    
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
    !*   -> latest changes: 2016/10/05                                                      *
    !*                                                                                      *
    !****************************************************************************************

    SUBROUTINE nametab()

    IMPLICIT NONE
    !...........................................................

    !list of neutral component names:
    NKname = "not_defined"
    NKnameTeX = "not_defined"

    !Alkanes
    NKname(1) = "Methane";  NKnameTeX(1) = "Methane";  
    NKname(2) = "Ethane";  NKnameTeX(2) = "Ethane";  
    NKname(3) = "Propane";  NKnameTeX(3) = "Propane";  
    NKname(4) = "Butane";  NKnameTeX(4) = "Butane";  
    NKname(5) = "Pentane";  NKnameTeX(5) = "Pentane";  
    NKname(6) = "Hexane";  NKnameTeX(6) = "Hexane";  
    NKname(7) = "Heptane";  NKnameTeX(7) = "Heptane";  
    NKname(8) = "Octane";  NKnameTeX(8) = "Octane";  
    NKname(9) = "Nonane";  NKnameTeX(9) = "Nonane";  
    NKname(10) = "Decane";  NKnameTeX(10) = "Decane";  
    NKname(11) = "Undecane";  NKnameTeX(11) = "Undecane";   
    NKname(12) = "Dodecane";  NKnameTeX(12) = "Dodecane";   
    NKname(13) = "Cyclohexane";  NKnameTeX(13) = "Cyclohexane";   
    NKname(14) = "Tetradecane";  NKnameTeX(14) = "Tetradecane";   
    NKname(15) = "Pentadecane";  NKnameTeX(15) = "Pentadecane";  
    NKname(16) = "Hexadecane";  NKnameTeX(16) = "Hexadecane";   
    NKname(17) = "Heptadecane";  NKnameTeX(17) = "Heptadecane";  
    NKname(18) = "Octadecane";  NKnameTeX(18) = "Octadecane";  
    NKname(19) = "Nonadecane";  NKnameTeX(19) = "Nonadecane";  
    NKname(20) = "Icosane";  NKnameTeX(20) = "Icosane";  
    NKname(21) = "Henicosane";  NKnameTeX(21) = "Henicosane";  
    NKname(22) = "Docosane";  NKnameTeX(22) = "Docosane";  
    NKname(23) = "Tricosane";  NKnameTeX(23) = "Tricosane";  
    NKname(24) = "Tetracosane";  NKnameTeX(24) = "Tetracosane";  
    NKname(25) = "Pentacosane";  NKnameTeX(25) = "Pentacosane";  
    NKname(26) = "Hexacosane";  NKnameTeX(26) = "Hexacosane";  
    NKname(33) = "Tridecane";  NKnameTeX(33) = "Tridecane";   !could put into 13			
    NKname(30) = "Squalane";  NKnameTeX(30) = "Squalane";  

    !...
    !Alkenes
    
    !...
    !Alcohols
    NKname(101) = "Methanol";   NKnameTeX(101) = "Methanol";  
    NKname(102) = "Ethanol";  NKnameTeX(102) = "Ethanol";  
    NKname(103) = "1-Propanol";  NKnameTeX(103) = "1-Propanol";  
    NKname(104) = "1-Butanol";  NKnameTeX(104) = "1-Butanol";  
    NKname(105) = "1-Pentanol";  NKnameTeX(105) = "1-Pentanol";  
    NKname(106) = "1-Hexanol";  NKnameTeX(106) = "1-Hexanol";  
    NKname(131) = "2-Propanol";  NKnameTeX(131) = "2-Propanol";  
    NKname(132) = "2-Butanol";  NKnameTeX(132) = "2-Butanol";  
    NKname(133) = "Isobutanol";  NKnameTeX(133) = "Isobutanol";  
    NKname(134) = "tert-Butanol";  NKnameTeX(134) = "$tert$-Butanol";  
    NKname(135) = "2-Pentanol";  NKnameTeX(135) = "2-Pentanol";  
    NKname(136) = "3-Pentanol";  NKnameTeX(136) = "3-Pentanol";  
    NKname(137) = "2-Methyl-2-butanol";  NKnameTeX(137) = "2-Methyl-2-butanol";  
    NKname(138) = "3-Methyl-1-butanol";  NKnameTeX(138) = "3-Methyl-1-butanol";  
    NKname(139) = "Cyclopentanol";  NKnameTeX(139) = "Cyclopentanol";  
    NKname(140) = "2-Hexanol";  NKnameTeX(140) = "2-Hexanol";  
    NKname(141) = "3-Hexanol";  NKnameTeX(141) = "3-Hexanol";  
    NKname(143) = "Cyclohexanol";  NKnameTeX(143) = "Cyclohexanol";  
    NKname(147) = "2-Octanol";  NKnameTeX(147) = "2-Octanol";  
    !142 = 4-methyl-2-pentanol, 144 = 2-heptanol, 
    !    145 = 3-heptanol, 146 = 4-heptanol, 147 = 2-octanol, 148 = 3-octanol,
    !    149 = 4-octanol, 
    NKname(150) = "2-Methyl-1-butanol";  NKnameTeX(150) = "2-Methyl-1-butanol";  

    !...
    !(mono) Carboxylic acids
    NKname(201) = "Formic_acid";  NKnameTeX(201) = "Formic_acid";  
    NKname(202) = "Acetic_acid";  NKnameTeX(202) = "Acetic_acid";  
    NKname(203) = "Propanoic_acid";  NKnameTeX(203) = "Propanoic_acid";    != propionic acid
    NKname(204) = "Butyric_acid";  NKnameTeX(204) = "Butyric_acid";  
    NKname(210) = "2-Methylpropanoic_acid";  NKnameTeX(210) = "2-Methylpropanoic_acid";    != isobutyric acid
    NKname(220) = "Pyruvic_acid";  NKnameTeX(220) = "Pyruvic_acid";    !(= 2-oxopropanoic acid)
    NKname(221) = "Methacrylic_acid";  NKnameTeX(221) = "Methacrylic_acid";    !(= 2-methyl-2-propenoic acid)
    NKname(222) = "Palmitic_acid";  NKnameTeX(222) = "Palmitic_acid";    !(= hexadecanoic acid)
    NKname(223) = "Stearic_acid";  NKnameTeX(223) = "Stearic_acid";    !(= Octadecanoic acid)
    NKname(224) = "Oleic_acid";  NKnameTeX(224) = "Oleic_acid";   !(= (9Z)-Octadec-9-enoic acid )
    NKname(225) = "Docosanoic_acid";  NKnameTeX(225) = "Docosanoic_acid";   !Behenic acid
    NKname(226) = "Tetracosanoic_acid";  NKnameTeX(226) = "Tetracosanoic_acid";  
    !...
    !Dicarboxylic acids
    NKname(301) = "Oxalic_acid";  NKnameTeX(301) = "Oxalic_acid";  
    NKname(302) = "Malonic_acid";  NKnameTeX(302) = "Malonic_acid";  
    NKname(303) = "Succinic_acid";  NKnameTeX(303) = "Succinic_acid";  
    NKname(304) = "Glutaric_acid";  NKnameTeX(304) = "Glutaric_acid";  
    NKname(305) = "Adipic_acid";  NKnameTeX(305) = "Adipic_acid";  
    NKname(306) = "Pimelic_acid";  NKnameTeX(306) = "Pimelic_acid";  
    NKname(307) = "Suberic_acid";  NKnameTeX(307) = "Suberic_acid";  
    NKname(308) = "Methylsuccinic_acid";  NKnameTeX(308) = "Methylsuccinic_acid";  
    NKname(309) = "Dimethylmalonic_acid";  NKnameTeX(309) = "Dimethylmalonic_acid";  
    NKname(310) = "2-Methylglutaric_acid";  NKnameTeX(310) = "2-Methylglutaric_acid";  
    NKname(311) = "3-Methylglutaric_acid";  NKnameTeX(311) = "3-Methylglutaric_acid";  
    NKname(312) = "2,2-Dimethylsuccinic_acid";  NKnameTeX(312) = "2,2-Dimethylsuccinic_acid";  
    NKname(313) = "3-Methyladipic_acid";  NKnameTeX(313) = "3-Methyladipic_acid";  
    NKname(314) = "3,3-Dimethylglutaric_acid";  NKnameTeX(314) = "3,3-Dimethylglutaric_acid";  
    NKname(315) = "Diethylmalonic_acid";  NKnameTeX(315) = "Diethylmalonic_acid";  
    NKname(316) = "Tartaric_acid";  NKnameTeX(316) = "Tartaric_acid";  
    NKname(317) = "Methylmalonic_acid";  NKnameTeX(317) = "Methylmalonic_acid";  
    NKname(318) = "Tartronic_acid";  NKnameTeX(318) = "Tartronic_acid";  

    NKname(320) = "Azelaic_acid";  NKnameTeX(320) = "Azelaic_acid";  
    NKname(321) = "Sebaic_acid";  NKnameTeX(321) = "Sebaic_acid";  
    NKname(322) = "Dodecanedioic_acid";  NKnameTeX(322) = "Dodecanedioic_acid";  
    !...
    !Water
    NKname(401) = "H2O";  NKnameTeX(401) = "H$_2$O";  
    !...
    !hydroperoxides
    NKname(411) = "hydroperoxyundecane";  NKnameTeX(411) = "hydroperoxyundecane";  
    NKname(412) = "hydroperoxydodecane";  NKnameTeX(412) = "hydroperoxydodecane";  
    NKname(413) = "hydroperoxytridecane";  NKnameTeX(413) = "hydroperoxytridecane";  
    NKname(414) = "hydroperoxytetradecane";  NKnameTeX(414) = "hydroperoxytetradecane";   
    NKname(415) = "hydroperoxypentadecane";  NKnameTeX(415) = "hydroperoxypentadecane";  
    NKname(416) = "hydroperoxyhexadecane";  NKnameTeX(416) = "hydroperoxyhexadecane";   
    NKname(417) = "hydroperoxyheptadecane";  NKnameTeX(417) = "hydroperoxyheptadecane";  
    NKname(418) = "hydroperoxyoctadecane";  NKnameTeX(418) = "hydroperoxyoctadecane";  
    NKname(419) = "hydroperoxynonadecane";  NKnameTeX(419) = "hydroperoxynonadecane";  
    NKname(420) = "hydroperoxyicosane";  NKnameTeX(420) = "hydroperoxyicosane";  
    NKname(421) = "hydroperoxyhenicosane";  NKnameTeX(421) = "hydroperoxyhenicosane";  
    NKname(422) = "hydroperoxydocosane";  NKnameTeX(422) = "hydroperoxydocosane";  
    NKname(423) = "hydroperoxytricosane";  NKnameTeX(423) = "hydroperoxytricosane";  
    NKname(424) = "hydroperoxytetracosane";  NKnameTeX(424) = "hydroperoxytetracosane";  
    NKname(425) = "hydroperoxypentacosane";  NKnameTeX(425) = "hydroperoxypentacosane";  
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
    NKname(501) = "Maleic_acid";  NKnameTeX(501) = "Maleic_acid";  
    NKname(502) = "Fumaric_acid";  NKnameTeX(502) = "Fumaric_acid";  
    NKname(503) = "Malic_acid";  NKnameTeX(503) = "Malic_acid";  
    NKname(504) = "Tartaric_acid";  NKnameTeX(504) = "Tartaric_acid";  
    NKname(505) = "Citric_acid";  NKnameTeX(505) = "Citric_acid";  
    NKname(506) = "alpha-Ketoglutaric_acid";  NKnameTeX(506) = "$\alpha$-Ketoglutaric_acid";  
    NKname(507) = "Glycolic_acid";  NKnameTeX(507) = "Glycolic_acid";  
    NKname(508) = "Methoxyacetic_acid";  NKnameTeX(508) = "Methoxyacetic_acid";  
    NKname(509) = "Shikimic_acid";  NKnameTeX(509) = "Shikimic_acid";  
    NKname(510) = "Isopimaric_acid";  NKnameTeX(510) = "Isopimaric_acid";  
    NKname(511) = "Abietic_acid";  NKnameTeX(511) = "Abietic_acid";  
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
    NKname(601) = "1,2-Ethanediol";  NKnameTeX(601) = "1,2-Ethanediol";  
    NKname(602) = "1,2-Propanediol";  NKnameTeX(602) = "1,2-Propanediol";  
    NKname(603) = "1,3-Propanediol";  NKnameTeX(603) = "1,3-Propanediol";  
    NKname(604) = "Glycerol";  NKnameTeX(604) = "Glycerol";  
    NKname(605) = "1,2-Butanediol";  NKnameTeX(605) = "1,2-Butanediol";  
    NKname(606) = "1,3-Butanediol";  NKnameTeX(606) = "1,3-Butanediol";  
    NKname(607) = "1,4-Butanediol";  NKnameTeX(607) = "1,4-Butanediol";  
    NKname(608) = "2,3-Butanediol";  NKnameTeX(608) = "2,3-Butanediol";  
    NKname(609) = "1,2,3-Butanetriol";  NKnameTeX(609) = "1,2,3-Butanetriol";  
    NKname(610) = "1,2,4-Butanetriol";  NKnameTeX(610) = "1,2,4-Butanetriol";  
    NKname(611) = "Erythritol";  NKnameTeX(611) = "Erythritol";  
    NKname(612) = "1,2-Pentanediol";  NKnameTeX(612) = "1,2-Pentanediol";  
    NKname(613) = "1,3-Pentanediol";  NKnameTeX(613) = "1,3-Pentanediol";  
    NKname(614) = "1,4-Pentanediol";  NKnameTeX(614) = "1,4-Pentanediol";  
    NKname(615) = "1,5-Pentanediol";  NKnameTeX(615) = "1,5-Pentanediol";  
    NKname(616) = "2,3-Pentanediol";  NKnameTeX(616) = "2,3-Pentanediol";  
    NKname(617) = "2,4-Pentanediol";  NKnameTeX(617) = "2,4-Pentanediol";  
    NKname(618) = "1,2-Hexanediol";  NKnameTeX(618) = "1,2-Hexanediol";  
    NKname(619) = "1,3-Hexanediol";  NKnameTeX(619) = "1,3-Hexanediol";  
    NKname(620) = "1,4-Hexanediol";  NKnameTeX(620) = "1,4-Hexanediol";  
    NKname(621) = "1,5-Hexanediol";  NKnameTeX(621) = "1,5-Hexanediol";  
    NKname(622) = "1,6-Hexanediol";  NKnameTeX(622) = "1,6-Hexanediol";  
    NKname(623) = "2,3-Hexanediol";  NKnameTeX(623) = "2,3-Hexanediol";  
    NKname(624) = "2,4-Hexanediol";  NKnameTeX(624) = "2,4-Hexanediol";  
    NKname(625) = "2,5-Hexanediol";  NKnameTeX(625) = "2,5-Hexanediol";  
    NKname(626) = "1,2,5-Hexanetriol";  NKnameTeX(626) = "1,2,5-Hexanetriol";  
    NKname(627) = "1,2,6-Hexanetriol";  NKnameTeX(627) = "1,2,6-Hexanetriol";  
    NKname(628) = "2,3,4-Hexanetriol";  NKnameTeX(628) = "2,3,4-Hexanetriol";  
    NKname(629) = "Sorbitol";  NKnameTeX(629) = "Sorbitol";  
    NKname(630) = "Mannitol";  NKnameTeX(630) = "Mannitol";  
    NKname(631) = "1,4-Heptanediol";  NKnameTeX(631) = "1,4-Heptanediol";  
    NKname(632) = "1,5-Heptanediol";  NKnameTeX(632) = "1,5-Heptanediol";  
    NKname(633) = "1,7-Heptanediol";  NKnameTeX(633) = "1,7-Heptanediol";  
    NKname(634) = "2,4-Heptanediol";  NKnameTeX(634) = "2,4-Heptanediol";  
    NKname(638) = "1,3-Nonanediol";  NKnameTeX(638) = "1,3-Nonanediol";  
    NKname(639) = "1,4-Dihydroxy-2-butene";  NKnameTeX(639) = "1,4-Dihydroxy-2-butene";   
    NKname(640) = "Levoglucosan";  NKnameTeX(640) = "Levoglucosan";  
    NKname(641) = "D-Fructopyranose";  NKnameTeX(641) = "D-Fructopyranose";  
    NKname(642) = "D-Mannopyranose";  NKnameTeX(642) = "D-Mannopyranose";  
    !.. 
    NKname(643) = "1,2,10-Decanetriol";  NKnameTeX(643) = "1,2,10-Decanetriol";  
    NKname(644) = "1,2,5,8-Octanetetrol";  NKnameTeX(644) = "1,2,5,8-Octanetetrol";  
    NKname(645) = "1,2,7,8-Octanetetrol";  NKnameTeX(645) = "1,2,7,8-Octanetetrol";  
    !...
    NKname(646) = "D-Ribofuranose";  NKnameTeX(646) = "D-Ribofuranose";  
    NKname(647) = "1,2,6-Hexanetriol";  NKnameTeX(647) = "1,2,6-Hexanetriol";  
    NKname(648) = "2-Methylerythritol";  NKnameTeX(648) = "2-Methylerythritol";  
    !...
    NKname(651) = "Xylose";  NKnameTeX(651) = "Xylose";  
    NKname(652) = "Glucose";  NKnameTeX(652) = "Glucose";  
    NKname(653) = "Fructose";  NKnameTeX(653) = "Fructose";  
    NKname(654) = "Sucrose";  NKnameTeX(654) = "Sucrose";  
    NKname(655) = "D-Glucopyranose";  NKnameTeX(655) = "D-Glucopyranose";  
    NKname(656) = "Raffinose";  NKnameTeX(656) = "Raffinose";  
    !...
    !Poly (ethylene glycols),  PEG
    NKname(659) = "PEG-300-n5";  NKnameTeX(659) = "PEG-300-n5";  
    NKname(660) = "PEG-300-n6";  NKnameTeX(660) = "PEG-300-n6";  
    NKname(661) = "PEG-3400-n76";  NKnameTeX(661) = "PEG-3400-n76";  
    NKname(662) = "PEG-2000-n44";  NKnameTeX(662) = "PEG-2000-n44";  
    NKname(663) = "PEG-200-n3";  NKnameTeX(663) = "PEG-200-n3";  
    NKname(664) = "PEG-200-n4";  NKnameTeX(664) = "PEG-200-n4";  
    NKname(665) = "PEG-8000-n180";  NKnameTeX(665) = "PEG-8000-n180";  
    NKname(666) = "PEG-20000-n453";  NKnameTeX(666) = "PEG-20000-n453";  
    NKname(667) = "PEG-6000-n135";  NKnameTeX(667) = "PEG-6000-n135";  
    NKname(668) = "PEG-3350-n75";  NKnameTeX(668) = "PEG-3350-n75";  
    NKname(669) = "PEG-4000-n89";  NKnameTeX(669) = "PEG-4000-n89";  
    NKname(670) = "PEG-10000-n226";  NKnameTeX(670) = "PEG-10000-n226";  
    NKname(671) = "PEG-400-n7";  NKnameTeX(671) = "PEG-400-n7";    !(PEG-400 as a mixture of polyethylene glycols with n = number of oxyethylene groups. 1/3 n=7 and 2/3 n=8 for the weight of approx. 400 g/mol) 
    NKname(672) = "PEG-400-n8";  NKnameTeX(672) = "PEG-400-n8";  
    NKname(673) = "PEG-1000-n21";  NKnameTeX(673) = "PEG-1000-n21";  
    NKname(674) = "PEG-1000-n22";  NKnameTeX(674) = "PEG-1000-n22";  
    NKname(675) = "PEG-600-n12";  NKnameTeX(675) = "PEG-600-n12";  
    NKname(676) = "PEG-600-n13";  NKnameTeX(676) = "PEG-600-n13";  
    NKname(677) = "PEG-1450-n31";  NKnameTeX(677) = "PEG-1450-n31";  
    NKname(678) = "PEG-1450-n32";  NKnameTeX(678) = "PEG-1450-n32";  
    NKname(679) = "PEG-1540-n33";  NKnameTeX(679) = "PEG-1540-n33";  
    NKname(680) = "PEG-1540-n34";  NKnameTeX(680) = "PEG-1540-n34";  
    !...
    NKname(681) = "C18H29(OH)9";  NKnameTeX(681) = "C$_{18}$H$_{29}$(OH)$_9$";    ![C18H29(OH)9] 1,3,5,7,9,11,13,15,17-octadecanol "super-polyol"
    NKname(682) = "2,2,6,6-Tetrakis(hydroxymethyl)cyclohexanol";  NKnameTeX(682) = "2,2,6,6-Tetrakis(hydroxymethyl)cyclohexanol";    !(=2,2,6,6-tetrakis(hydroxylmethyl)cyclohexanol, C10)
    !...
    !Polcyclic aromatic hydrocarbons (PAH)
    NKname(701) = "Benzene";  NKnameTeX(701) = "Benzene";  
    NKname(702) = "Naphthalene";  NKnameTeX(702) = "Naphthalene";  
    NKname(703) = "Anthracene";  NKnameTeX(703) = "Anthracene";  
    NKname(704) = "Phenanthrene";  NKnameTeX(704) = "Phenanthrene";  
    NKname(705) = "Fluoranthene";  NKnameTeX(705) = "Fluoranthene";  
    NKname(706) = "Pyrene";  NKnameTeX(706) = "Pyrene";  
    NKname(707) = "Chrysene";  NKnameTeX(707) = "Chrysene";  
    NKname(708) = "Perylene";  NKnameTeX(708) = "Perylene";  
    NKname(709) = "Coronene";  NKnameTeX(709) = "Coronene";  
    NKname(710) = "Toluene";  NKnameTeX(710) = "Toluene";  
    NKname(711) = "2-Methylnaphthalene";  NKnameTeX(711) = "2-Methylnaphthalene";  
    NKname(712) = "7-Isopropyl-1-methylphenanthrene";  NKnameTeX(712) = "7-Isopropyl-1-methylphenanthrene";   != Retene
    NKname(713) = "Benz[a]anthracene";  NKnameTeX(713) = "Benz[a]anthracene";  
    NKname(714) = "Benzo[b]fluoranthene";  NKnameTeX(714) = "Benzo[b]fluoranthene";  
    NKname(715) = "Indeno[1,2,3-cd]pyrene";  NKnameTeX(715) = "Indeno[1,2,3-cd]pyrene";  
    NKname(716) = "Benzo[a]pyrene";  NKnameTeX(716) = "Benzo[a]pyrene";  
    NKname(717) = "6,12-dihydroxy-BaP";  NKnameTeX(717) = "6,12-dihydroxy-BaP";  
    NKname(718) = "BaP-6,12-dione";  NKnameTeX(718) = "BaP-6,12-dione";  
    NKname(719) = "BaP-derived_diCOOH";  NKnameTeX(719) = "BaP-derived_diCOOH";   !7-oxo-7H-benzo[de]anthracene-3,4-dicarboxylic acid
    !...
    !Multifunctional aromatic hydrocarbons, phenolic compounds
    NKname(720) = "Phenol";  NKnameTeX(720) = "Phenol";  
    NKname(721) = "Protocatechuic_acid";  NKnameTeX(721) = "Protocatechuic_acid";  
    NKname(722) = "Vanillin";  NKnameTeX(722) = "Vanillin";   !4-Hydroxy-3-methoxybenzaldehyde
    NKname(723) = "Vanillic_acid";  NKnameTeX(723) = "Vanillic_acid";  
    NKname(724) = "Gallic_acid";  NKnameTeX(724) = "Gallic_acid";  
    NKname(725) = "Ferulic_acid";  NKnameTeX(725) = "Ferulic_acid";  
    NKname(726) = "Syringic_acid";  NKnameTeX(726) = "Syringic_acid";  
    NKname(727) = "2-Hydroxybenzoic_acid";  NKnameTeX(727) = "2-Hydroxybenzoic_acid";   !(= Salicylic acid)
    NKname(728) = "3-Hydroxybenzoic_acid";  NKnameTeX(728) = "3-Hydroxybenzoic_acid";  
    NKname(729) = "4-Hydroxybenzoic_acid";  NKnameTeX(729) = "4-Hydroxybenzoic_acid";  
    NKname(730) = "Phthalic_acid";  NKnameTeX(730) = "Phthalic_acid";    ! (= benzene-1,2-dicarboxylic acid)
    NKname(731) = "2,4-Dihydroxybenzaldehyde";  NKnameTeX(731) = "2,4-Dihydroxybenzaldehyde";  
    NKname(732) = "Vanillylmandelic_acid";  NKnameTeX(732) = "Vanillylmandelic_acid";  
    NKname(733) = "3,5-Dihydroxybenzoic_acid";  NKnameTeX(733) = "3,5-Dihydroxybenzoic_acid";  
    NKname(734) = "Mandelic_acid";  NKnameTeX(734) = "Mandelic_acid";  
    NKname(735) = "Dimethyl_phthalate";  NKnameTeX(735) = "Dimethyl_phthalate";  
    NKname(736) = "2,5-Dihydroxybenzoic_acid";  NKnameTeX(736) = "2,5-Dihydroxybenzoic_acid";  
    NKname(737) = "Resorcinol";  NKnameTeX(737) = "Resorcinol";   !(= 1,3-benzenediol, = m-benzenediol)
    NKname(738) = "p-Cresol";  NKnameTeX(738) = "p-Cresol";  
    !m-Cresol and o-Cresol are listed under 895, 896
    NKname(739) = "p-Benzenediol";  NKnameTeX(739) = "p-Benzenediol";   !(= 1,4-benzenediol)
    NKname(740) = "4-Methylguaiacol";  NKnameTeX(740) = "4-Methylguaiacol";  
    NKname(741) = "4-Propylguaiacol";  NKnameTeX(741) = "4-Propylguaiacol";  
    NKname(742) = "Coniferaldehyde";  NKnameTeX(742) = "Coniferaldehyde";  
    NKname(743) = "4-Methylsyringol";  NKnameTeX(743) = "4-Methylsyringol";   
    NKname(744) = "Syringyl_acetone";  NKnameTeX(744) = "Syringyl_acetone";   !(= 4-(4-hydroxy-3,5-dimethoxyphenyl)butan-2-one )
    NKname(745) = "p-Tolualdehyde";  NKnameTeX(745) = "p-Tolualdehyde";   !4-Methylbenzaldehyde
    NKname(746) = "2,5-Dimethylbenzaldehyde";  NKnameTeX(746) = "2,5-Dimethylbenzaldehyde";  
    NKname(747) = "Indan-1-one";  NKnameTeX(747) = "Indan-1-one";  
    NKname(748) = "1H-phenalen-1-one";  NKnameTeX(748) = "1H-phenalen-1-one";  
    NKname(749) = "3,4-Dimethoxytoluene";  NKnameTeX(749) = "3,4-Dimethoxytoluene";  
    NKname(750) = "Veratric_acid";  NKnameTeX(750) = "Veratric_acid";  
    NKname(751) = "Benzoic_acid";  NKnameTeX(751) = "Benzoic_acid";  
    NKname(752) = "2-Naphthol";  NKnameTeX(752) = "2-Naphthol";  
    NKname(753) = "1,4-Naphthalenedione";  NKnameTeX(753) = "1,4-Naphthalenedione";  
    NKname(754) = "5-Hydroxy-1,4-naphthalenedione";  NKnameTeX(754) = "5-Hydroxy-1,4-naphthalenedione";  
    NKname(755) = "2-Carboxycinnamic_acid";  NKnameTeX(755) = "2-Carboxycinnamic_acid";  
    NKname(756) = "4-Hydroxyphthalic_acid";  NKnameTeX(756) = "4-Hydroxyphthalic_acid";  
    NKname(757) = "Carminic_acid";  NKnameTeX(757) = "Carminic_acid";  
    NKname(758) = "Triolein";  NKnameTeX(758) = "Triolein";  
    NKname(759) = "Linoleic_acid";  NKnameTeX(759) = "Linoleic_acid";  
    !...
    !Carbonyls (Ketones and Aldehydes)
    NKname(802) = "Acetone";  NKnameTeX(802) = "Acetone";  
    NKname(803) = "2-Butanone";  NKnameTeX(803) = "2-Butanone";  
    NKname(804) = "3-Methyl-2-butanone";  NKnameTeX(804) = "3-Methyl-2-butanone";   !Methyl Isopropyl Ketone
    NKname(805) = "4-Methyl-2-pentanone";  NKnameTeX(805) = "4-Methyl-2-pentanone";  
    NKname(806) = "2-Pentanone";  NKnameTeX(806) = "2-Pentanone";  
    NKname(807) = "2-Hexanone";  NKnameTeX(807) = "2-Hexanone";  
    NKname(808) = "2-Heptanone";  NKnameTeX(808) = "2-Heptanone";  
    NKname(809) = "3-Heptanone";  NKnameTeX(809) = "3-Heptanone";  
    NKname(810) = "3-Pentanone";  NKnameTeX(810) = "3-Pentanone";  
    NKname(811) = "2-Octanone";  NKnameTeX(811) = "2-Octanone";  
    NKname(812) = "4-Heptanone";  NKnameTeX(812) = "4-Heptanone";  
    NKname(813) = "Acetaldehyde";  NKnameTeX(813) = "Acetaldehyde";  
    NKname(814) = "Propionaldehyde";  NKnameTeX(814) = "Propionaldehyde";  
    NKname(815) = "Butyraldehyde";  NKnameTeX(815) = "Butyraldehyde";  
    NKname(816) = "2-Nonanone";  NKnameTeX(816) = "2-Nonanone";  
    NKname(817) = "Dodecanone";  NKnameTeX(817) = "Dodecanone";  
    NKname(818) = "Tridecanone";  NKnameTeX(818) = "Tridecanone";  
    NKname(819) = "Tetradecanone";  NKnameTeX(819) = "Tetradecanone";  
    NKname(820) = "Pentadecanone";  NKnameTeX(820) = "Pentadecanone";  
    NKname(821) = "Hexadecanone";  NKnameTeX(821) = "Hexadecanone";  
    NKname(822) = "Heptadecanone";  NKnameTeX(822) = "Heptadecanone";  
    NKname(823) = "Octadecanone";  NKnameTeX(823) = "Octadecanone";  
    NKname(824) = "Nonadecanone";  NKnameTeX(824) = "Nonadecanone";  
    NKname(825) = "Icosanone";  NKnameTeX(825) = "Icosanone";  
    NKname(826) = "Henicosanone";  NKnameTeX(826) = "Henicosanone";  
    NKname(827) = "Docosanone";  NKnameTeX(827) = "Docosanone";  
    NKname(828) = "Tricosanone";  NKnameTeX(828) = "Tricosanone";  
    NKname(829) = "Tetracosanone";  NKnameTeX(829) = "Tetracosanone";  
    NKname(830) = "Pentacosanone";  NKnameTeX(830) = "Pentacosanone";  
    NKname(831) = "Undecanone";  NKnameTeX(831) = "Undecanone";  

    !...
    !Esters
    NKname(840) = "bis(2-ethylhexyl)_sebacate";  NKnameTeX(840) = "bis(2-ethylhexyl)_sebacate";   !BES
    NKname(841) = "Methyl_acetate";  NKnameTeX(841) = "Methyl_acetate";  
    NKname(842) = "Ethyl_acetate";  NKnameTeX(842) = "Ethyl_acetate";  
    NKname(843) = "1-Propyl_acetate";  NKnameTeX(843) = "1-Propyl_acetate";  
    NKname(844) = "1-Butyl_acetate";  NKnameTeX(844) = "1-Butyl_acetate";  
    NKname(845) = "Isobutyl_acetate";  NKnameTeX(845) = "Isobutyl_acetate";  
    NKname(846) = "2-Butyl_acetate";  NKnameTeX(846) = "2-Butyl_acetate";  
    NKname(847) = "tert-Butyl_acetate";  NKnameTeX(847) = "tert-Butyl_acetate";   
    NKname(848) = "1-Pentyl_acetate";  NKnameTeX(848) = "1-Pentyl_acetate";  
    NKname(849) = "1-Hexyl_acetate";  NKnameTeX(849) = "1-Hexyl_acetate";  
    NKname(850) = "2-Ethoxyethyl_acetate";  NKnameTeX(850) = "2-Ethoxyethyl_acetate";  
    NKname(851) = "Octadecanoic_acid_methyl_ester";  NKnameTeX(851) = "Octadecanoic_acid_methyl_ester";   !(= methyl Stearate)
    NKname(852) = "Octadecanoic_acid_ethyl_ester";  NKnameTeX(852) = "Octadecanoic_acid_ethyl_ester";   !(= ethyl Stearate)

    !Diketones
    NKname(853) = "Undecanedione";  NKnameTeX(853) = "Undecanedione";  
    NKname(854) = "Dodecanedione";  NKnameTeX(854) = "Dodecanedione";  
    NKname(855) = "Tridecanedione";  NKnameTeX(855) = "Tridecanedione";  
    NKname(856) = "Tetradecanedione";  NKnameTeX(856) = "Tetradecanedione";  
    NKname(857) = "Pentadecanedione";  NKnameTeX(857) = "Pentadecanedione";  
    NKname(858) = "Hexadecanedione";  NKnameTeX(858) = "Hexadecanedione";  
    NKname(859) = "Heptadecanedione";  NKnameTeX(859) = "Heptadecanedione";  
    NKname(860) = "Octadecanedione";  NKnameTeX(860) = "Octadecanedione";  
    NKname(861) = "Nonadecanedione";  NKnameTeX(861) = "Nonadecanedione";  
    NKname(862) = "Icosanedione";  NKnameTeX(862) = "Icosanedione";  
    NKname(863) = "Henicosanedione";  NKnameTeX(863) = "Henicosanedione";  
    NKname(864) = "Docosanedione";  NKnameTeX(864) = "Docosanedione";  
    NKname(865) = "Tricosanedione";  NKnameTeX(865) = "Tricosanedione";  

    !...
    !Ethers
    NKname(881) = "2-Methoxy-2-methylpropane";  NKnameTeX(881) = "2-Methoxy-2-methylpropane";  
    NKname(882) = "2-Methoxyethanol";  NKnameTeX(882) = "2-Methoxyethanol";  
    NKname(883) = "2-Ethoxyethanol";  NKnameTeX(883) = "2-Ethoxyethanol";  
    NKname(884) = "1-Methoxy-2-propanol";  NKnameTeX(884) = "1-Methoxy-2-propanol";  
    NKname(885) = "2-Isopropoxyethanol";  NKnameTeX(885) = "2-Isopropoxyethanol";   
    NKname(886) = "2-Butoxyethanol";  NKnameTeX(886) = "2-Butoxyethanol";   
    NKname(887) = "2-Methoxypropanol";  NKnameTeX(887) = "2-Methoxypropanol";   
    NKname(888) = "1-(2-methoxypropoxy)-2-propanol";  NKnameTeX(888) = "1-(2-methoxypropoxy)-2-propanol";   
    NKname(889) = "2-(2-methoxyethoxy)ethanol";  NKnameTeX(889) = "2-(2-methoxyethoxy)ethanol";   
    NKname(890) = "2-(2-ethoxyethoxy)ethanol";  NKnameTeX(890) = "2-(2-ethoxyethoxy)ethanol";   
    NKname(891) = "1,4-Dioxane";  NKnameTeX(891) = "1,4-Dioxane";   
    NKname(892) = "Tetrahydrofuran";  NKnameTeX(892) = "Tetrahydrofuran";  
    NKname(893) = "2-Ethoxy-2-methylpropane";  NKnameTeX(893) = "2-Ethoxy-2-methylpropane";  
    NKname(894) = "Diethyl_ether";  NKnameTeX(894) = "Diethyl_ether";  
    !...
    !other complex functionalyzed compounds, MCM names:
    NKname(895) = "m-Cresol";  NKnameTeX(895) = "m-Cresol";  
    NKname(896) = "o-Cresol";  NKnameTeX(896) = "o-Cresol";  
    NKname(897) = "Anisaldehyde";  NKnameTeX(897) = "Anisaldehyde";  
    NKname(898) = "Trehalose";  NKnameTeX(898) = "Trehalose";  
    NKname(899) = "Maltose";  NKnameTeX(899) = "Maltose";  

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
    NKname(959) = "C15H18O7_syringol_aqSOA";  NKnameTeX(959) = "C$_{15}$H$_{18}$O$_{7}$_syringol_aqSOA";   !complete name: 6-(2,4-dihydroxy-3,5-dimethoxyphenyl)-5-methoxyhexa-3,5-dienoic acid
    !....................................................

    !add new organics for MT system used in project 2016-05-13 by JM
    NKname(960) = "C813OH";  NKnameTeX(960) = "C813OH";   !complete name: ?
    NKname(961) = "C922OOH";  NKnameTeX(961) = "C922OOH";   !complete name: ?
    NKname(962) = "C98OH";  NKnameTeX(962) = "C98OH";   !complete name: ?
    NKname(963) = "C813O2";  NKnameTeX(963) = "C813O2";   !complete name: ?
    NKname(964) = "C922O2";  NKnameTeX(964) = "C922O2";   !complete name: ?
    NKname(965) = "C811PAN";  NKnameTeX(965) = "C811PAN";   !complete name: ?
    NKname(966) = "C920PAN";  NKnameTeX(966) = "C920PAN";   !complete name: ?

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

    
    !...
    NKname(1500) = "!fromInpFile!";  NKnameTeX(1500) = "!fromInpFile!";   !this ID number is used to indicate that the name is set from input via a file.

    
    !list of electrolyte component names:
    electname = "not_defined"
    electnameTeX = "not_defined"

    !Li+ <-> Cl-
    electname(1,2) = "LiCl"
    electnameTeX(1,2) = "LiCl"
    !Li+ <-> Br-  
    electname(1,3) = "LiBr"
    electnameTeX(1,3) = "LiBr"
    !Li+ <-> NO3-   
    electname(1,5) = "LiNO3"
    electnameTeX(1,5) = "LiNO$_3$"                   
    !Li+ <-> SO4--   
    electname(1,21) = "Li2SO4"
    electnameTeX(1,21) = "Li$_2$SO$_4$"    
    !Na+ <-> Cl-   
    electname(2,2) = "NaCl"
    electnameTeX(2,2) = "NaCl"         	
    !Na+ <-> Br-   
    electname(2,3) = "NaBr"
    electnameTeX(2,3) = "NaBr"
    !Na+ <-> NO3-   
    electname(2,5) = "NaNO3"
    electnameTeX(2,5) = "NaNO$_3$"  
    !Na+ <-> HSO4-  
    electname(2,8) = "NaHSO4"
    electnameTeX(2,8) = "NaHSO$_4$"  
    !NaCH3SO3 (-> CH3SO3- & Na+) != Na-MSA = sodium methanesulfonate
    electname(2,9) = "NaCH3SO3"
    electnameTeX(2,9) = "NaCH$_3$SO$_3$"   
    !Na+ <-> SO4-- 
    electname(2,21) = "Na2SO4"
    electnameTeX(2,21) = "Na$_2$SO$_4$"    
    !K+ <-> Cl-  
    electname(3,2) = "KCl"
    electnameTeX(3,2) = "KCl"
    !K+ <-> Br-  
    electname(3,3) = "KBr"
    electnameTeX(3,3) = "KBr"
    !K+ <-> NO3-   
    electname(3,5) = "KNO3"
    electnameTeX(3,5) = "KNO$_3$"   
    !K+ <-> SO4--  
    electname(3,21) = "K2SO4"
    electnameTeX(3,21) = "K$_2$SO$_4$" 
    !NH4+ <-> Cl-   
    electname(4,2) = "NH4Cl"
    electnameTeX(4,2) = "NH$_4$Cl"            
    !NH4+ <-> Br- 
    electname(4,3) = "NH4Br"
    electnameTeX(4,3) = "NH$_4$Br"  
    !NH4+ <-> NO3- 
    electname(4,5) = "NH4NO3"
    electnameTeX(4,5) = "NH$_4$NO$_3$"  
    !NH4+ <-> HSO4- 
    electname(4,8) = "NH4HSO4"
    electnameTeX(4,8) = "NH$_4$HSO$_4$" 
    !NH4CH3SO3 (-> CH3SO3- & NH4+) != NH4-MSA = sodium methanesulfonate
    electname(4,9) = "NH4CH3SO3"
    electnameTeX(4,9) = "NH$_4$CH$_3$SO$_3$"   
    !(NH4)2SO4  
    electname(4,21) = "(NH4)2SO4"
    electnameTeX(4,21) = "(NH$_4$)$_2$SO$_4$" 
    !H+ <-> Cl-  
    electname(5,2) = "HCl"
    electnameTeX(5,2) = "HCl"  
    !H+ <-> Br-  
    electname(5,3) = "HBr"
    electnameTeX(5,3) = "HBr"   
    !H+ <-> NO3- 
    electname(5,5) = "HNO3"
    electnameTeX(5,5) = "HNO$_3$" 
    !H2SO4 (- > HSO4- & H+) 
    electname(5,8) = "H2SO4"
    electnameTeX(5,8) = "H$_2$SO$_4$"    
    !HCH3SO3 (- > CH3SO3- & H+) != MSA = methanesulfonic acid
    electname(5,9) = "MSA"
    electnameTeX(5,9) = "MSA"   
    !HSO4-  (- >2x H+ & SO4--) 
    electname(5,21) = "H2SO4"
    electnameTeX(5,21) = "H$_2$SO$_4$"     
    !Ca++ <-> Cl- 
    electname(21,2) = "CaCl2"
    electnameTeX(21,2) = "CaCl$_2$"     
    !Ca++ <-> Br- 
    electname(21,3) = "CaBr2"
    electnameTeX(21,3) = "CaBr$_2$" 
    !Ca++ <-> NO3-  
    electname(21,5) = "Ca(NO3)2"
    electnameTeX(21,5) = "Ca(NO$_3$)$_2$"                        
    !Ca++ <-> SO4-- 
    electname(21,21) = "CaSO4"
    electnameTeX(21,21) = "CaSO$_4$"         
    ! Mg++ <-> Cl-  
    electname(23,2) = "MgCl2"
    electnameTeX(23,2) = "MgCl$_2$"      
    ! Mg++ <-> Br-  
    electname(23,3) = "MgBr2"
    electnameTeX(23,3) = "MgBr$_2$" 
    ! Mg++ <-> NO3-  
    electname(23,5) = "Mg(NO3)2"
    electnameTeX(23,5) = "Mg(NO$_3$)$_2$"                      
    ! Mg++ <-> SO4-- 
    electname(23,21) = "MgSO4"
    electnameTeX(23,21) = "MgSO$_4$" 
    !...................................................

    END SUBROUTINE nametab

!=================================================================================================================================
    
    
    !**********************************************************************************************************
    !*                                                                                                        *
    !*  Subroutine to set the character strings for the names of the different components in a mixture.       *
    !*                                                                                                        *
    !*                                                                                                        *   
    !*                                      (c) Andi Zuend, IACETH, ETH Zurich, 2009                          *
    !*                         Div. Chem. Engineering, California Institute of Technology, 2009 - 2012        *
    !**********************************************************************************************************
    PURE SUBROUTINE names_mix(ncomp, nneutral, CompN, compname, compnameTeX, ionname, ionnameTeX, OtoCratio, HtoCratio)

    USE ModSystemProp, ONLY : ElectComps, ElectSubs, Ncation, Nanion, NGS, frominpfile, cpname
    USE ModSubgroupProp, ONLY : subgrname, subgrnameTeX
    USE ModSubgroupProp, ONLY : OtoCandHtoCmix

    IMPLICIT NONE

    !interface vars:
    INTEGER(4),INTENT(IN) :: ncomp, nneutral  !the component number (internal numbering as the compN in definemixtures, LRdata, etc.)
    INTEGER(4),DIMENSION(ncomp),INTENT(IN) :: CompN
    REAL(8),DIMENSION(ncomp),INTENT(OUT) :: OtoCratio, HtoCratio
    CHARACTER(LEN=*),DIMENSION(ncomp),INTENT(OUT) :: compname, compnameTeX
    CHARACTER(LEN=*),DIMENSION(NGS),INTENT(OUT) :: ionname, ionnameTeX
    !local vars:
    CHARACTER(LEN=16) :: txt
    INTEGER(4) :: i, k, cnt, cn, an
    REAL(8) :: OtoC, HtoC
    REAL(8),DIMENSION(nneutral) :: xfrac
    !...........................................................

    OtoCratio = -7.777777D0  !set to an impossible value at initialization
    HtoCratio = -7.777777D0
    !Get the names of the actual components in mixture nd
    !neutral components:
    !water is always component 1 (if there is water in the mixture):
    IF (CompN(1) == 401) THEN !there is water
        compname(1) = TRIM(NKname(401)) !//"'" 
        compnameTeX(1) = TRIM(NKnameTeX(401)) !//"'"
    ENDIF

    DO i = 1,nneutral
        k = CompN(i)
        IF (k > 0 .AND. k /= 401) THEN !scan all components, except water
            IF (frominpfile) THEN
                compname(i) = cpname(i) 
                compnameTeX(i) = cpname(i) 
            ELSE
                compname(i) = TRIM(NKname(k)) !//"'"    !the "'" is added to make sure that for comma separated reading the commas in e.g. 1,2-Butanediol is not separating the species name.
                compnameTeX(i) = TRIM(NKnameTeX(k)) !//"'"
            ENDIF
            !call the OtoC calculation subroutine which uses atom information from subgroups;  here just for the single components.
            xfrac = 0.0D0
            xfrac(i) = 1.0D0 !set only the current component to 1, in order to compute only its O:C ratio.
            CALL OtoCandHtoCmix(nneutral, xfrac, OtoC, HtoC)
            OtoCratio(i) = OtoC
            HtoCratio(i) = HtoC
        ENDIF
    ENDDO

    !electrolyte names are set in 'nametab':
    cnt = Ncation*Nanion
    DO i = 1,cnt
        cn = ElectComps(i,1)
        an = ElectComps(i,2)  
        cn = cn -200 !index shift for cations
        an = an -240 !index shift for anions
        compname(nneutral+i) = TRIM(electname(cn,an))
        compnameTeX(nneutral+i) = TRIM(electnameTeX(cn,an)) 
    ENDDO
    !set the rest of the names to an empty string
    compname(nneutral+cnt+1:) = ""
    compnameTeX(nneutral+cnt+1:) = ""

    !set also ion names, especially for output use with systems containing several electrolytes:
    !cations:
    cnt = 0
    DO i = 1,NGS !loop over ElectSubs
        k = ElectSubs(i)
        IF (k < 240) THEN !cation
            cnt = cnt +1
            txt = TRIM(subgrname(k))
            cn = VERIFY(txt, "( )", BACK = .false.) !first character index that is not ( or ) or a space
            an = VERIFY(txt, "( )", BACK = .true.) !last character index that is not ( or ) or a space
            !remove the parathensis (...) around the ion subgroup:
            ionname(cnt) = TRIM(txt(cn:an))
            txt = TRIM(subgrnameTeX(k))
            cn = VERIFY(txt, "( )", BACK = .false.) !first character index that is not ( or ) or a space
            an = VERIFY(txt, "( )", BACK = .true.) !last character index that is not ( or ) or a space
            ionnameTeX(cnt) = TRIM(txt(cn:an))
        ENDIF
    ENDDO
    !anions:
    DO i = 2,NGS !loop over ElectSubs
        k = ElectSubs(i)
        IF (k > 240) THEN !anion
            cnt = cnt +1
            txt = TRIM(subgrname(k))
            cn = VERIFY(txt, "( )", BACK = .false.) !first character index that is not ( or ) or a space
            an = VERIFY(txt, "( )", BACK = .true.) !last character index that is not ( or ) or a space
            !remove the parathensis (...) around the ion subgroup:
            ionname(cnt) = TRIM(txt(cn:an))
            txt = TRIM(subgrnameTeX(k))
            cn = VERIFY(txt, "( )", BACK = .false.) !first character index that is not ( or ) or a space
            an = VERIFY(txt, "( )", BACK = .true.) !last character index that is not ( or ) or a space
            ionnameTeX(cnt) = TRIM(txt(cn:an))
        ENDIF
    ENDDO

    END SUBROUTINE names_mix
!================================================================================================================================= 
    
END MODULE ModComponentNames