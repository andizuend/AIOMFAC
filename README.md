# AIOMFAC-web, Public Model Code Repository
This public repository provides the AIOMFAC model Fortran code (AIOMFAC-web, version 2.20 and newer) and additional information about building and running the model on your own system. May it be of use to you.

## About AIOMFAC
AIOMFAC stands for Aerosol Inorganic&ndash;Organic Mixtures Functional groups Activity Coefficients; it is a thermodynamic group-contribution model to describe non-ideal mixing in liquid solutions (phases). If you are unfamiliar with the purpose and applications of AIOMFAC, please visit the [AIOMFAC website](https://aiomfac.lab.mcgill.ca) for more information.

----
> [!TIP]
> Click on the <a href="#"><img src="./Auxiliary/outline_icon.jpg" alt="outline" style="width:3ex"/></a> icon at the top right of this readme file to show the *table of contents* of this file with links to specific sections.

## AIOMFAC-web model versions and related code and feature updates
Information associated with specific model versions, including comments on major AIOMFAC-web changes and new features, are provided on the page accompanying that release; see under [releases](https://github.com/andizuend/AIOMFAC/releases). 

## Applications, Modifications and Citation
If you use any of our AIOMFAC code in your own projects / code, following the GNU license restrictions, we would appreciate hearing about it. In scientific or other publications, we also request that you reference the main peer-reviewed publications which describe the theoretical underpinning of the AIOMFAC model and its parameterizations, as described in more detail on the AIOMFAC website under: https://aiomfac.lab.mcgill.ca/citation.html.

### License
All files presented here are covered under the GNU GPL license v3.0. For more information, please read the license file. A brief overview of the viable permissions can be found here: https://choosealicense.com/licenses/.

## Dependencies
- Starting with AIOMFAC-web v3.14, which introduced the support of pure-component viscosity predictions via a machine learning method by [Armeli et al. (2023)](https://dx.doi.org/10.1021/acsomega.2c08146) implemented in Python, there are several specific Python packages that will need to be installed alongside the Fortran program using a dedicated virtual environment (follow the installation instructions provided below).
- The Fortran code itself is dependency-free. It requires a compiler supporting the Fortran 2008 standard (or newer). For example, the gfortran v9 and newer (v12, v15) and the Intel oneAPI ifx compilers have been tested and work as expected for our code, but any other recent Fortran compiler should be able to compile the Fortran sources.

## Installation instructions
> [!NOTE] 
> The following steps are first outlined for a Windows 64-bit installation (denoted by steps tagged as [Windows]). Equivalent steps are also shown for installation on a Linux machine (denoted by tag [Linux]). The Linux steps were tested with RHEL v8.1; the details for other Linux distributions may differ slightly.

### (1) Relative folder structure
Copy/clone the AIOMFAC folders and contained files from this repository to your local project (e.g. from command terminal when in your desired parent directory enter `git clone https://github.com/andizuend/AIOMFAC.git`).
On Linux, the main folder structure should look as illustrated below (not showing all subfolders of the `.venv` directory). On Windows the structure is the same but the subfolders inside `.venv` differ. The `.venv` content will get generated automatically; see step (2) below.

```
AIOMFAC
├───Auxiliary
├───FortranCode
├───Inputfiles
├───Outputfiles
└───TgML_Armeli
    ├───.venv
    │   ├───bin
    │   ├───lib
    │   └───include
    ├───InputFiles
    ├───OutputFiles
    └───pickle
```

> [!NOTE] 
> Within folder `TgML_Armeli`, the subfolders `InputFiles` and `OutputFiles` need to exist (with read and write permissions set for the current user). During normal operation of the AIOMFAC program with use of the glass transition temperature prediction based on the machine learning method by [Armeli et al. (2023)](https://dx.doi.org/10.1021/acsomega.2c08146), temporary files may be created in those folders and deleted a moment later. That's why they will look unused, but are needed for the proper functioning of the setup.

### (2) Generate a (virtual) Python environment
For reasons of compatibility with the machine learning methods run in the background (called from the AIOMFAC Fortran program), it is necessary to install Python v3.9, e.g. specific version 3.9.13, in a virtual environment together with the specific Python packages outlined in the following steps:
- In a command prompt run on [Windows]  `py --list` or on [Linux]  `compgen -c python | grep -E '^python[0-9.]+$' ` to see the Python versions already installed on the system. 
- If Python 3.9 if not among them, install it on the system (consult a guide for your operating system if it is unclear to you how to do this correctly).
- Create a virtual environment inside the `TgML_Armeli` folder. In a command prompt (or terminal), navigate to the `TgML_Armeli` folder and execute the command:
    - [Windows]    `py -3.9 -m venv .venv`
    - [Linux]    `python3.9 -m venv .venv`
- Activate the virtual environment using the command:
    - [Windows]    `.\.venv\Scripts\activate.bat`
    - [Linux]    `source .venv/bin/activate`
- Given the activated Python environment in the command prompt, use pip to install the specific package versions listed in the following:
    -  `pip install numpy==1.22.4`
	-  `pip install deepchem==2.5.0`
	-  `pip install rdkit-pypi==2022.3.2.1`
	-  `pip install rdkit==2022.9.1`
	-  `pip install scikit-learn==1.1.1`
	-  `pip install tensorflow-cpu==2.9.0`

### (3) Test the TgML_Armeli Python code execution
- From a command prompt when navigated to the `TgML_Armeli` folder, execute the command:
	- [Windows] 	`.venv\Scripts\python.exe TgML_SMILES.py` 
	- [Linux] 	`.venv/bin/python TgML_SMILES.py`
- Running the above may take a few seconds since large Python packages are first imported. If the test was successful you should see a message in the terminal stating "done with processing 1 SMILES..." and "Note: all SMILES were confirmed to be valid.". Further, in folder `\OutputFiles`, you will find a new file `output_1000_Tg.txt`. If this test was unsuccessful, check the error message issued and investigate whether all the above listed Python packages were installed successfully into the `.venv`.

### (4) Compile and link the AIOMFAC Fortran program
Building the AIOMFAC program from the Fortran source code can be done in a few distinct ways outlined in the following. If all you wish to do is to generate the executable AIOMFAC program on your system to subsequently  run your customized input files / cases, it is recommended to build the program using the instructions provided in the file `build_command_line.txt` included in folder `FortranCode`. We recommend using either [GNU's gfortran](https://gcc.gnu.org/fortran/) or [Intel's oneAPI ifx](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) compiler, both of which have been confirmed to successfully compile the Fortran source files. Other modern Fortran compilers should work as well (untested). 
Briefly, for command line compilation, the following steps need to be completed (examples described below apply to AIOMFAC-web v3.14 and later):
- On [Windows]:
	- open a dedicated Intel oneAPI terminal (which provides the necessary environment variable settings) or, alternatively, use the Windows subsystem for Linux (WSL) from a regular terminal -- in that case follow the instructions for building the program as described for Linux;
	- navigate to your local AIOMFAC Fortran source code directory;
	- copy & paste, then execute the following command line (for compilation with /O3 optimization set using ifx, example works for AIOMFAC-web v3.14 and later):
	```
	ifx /o AIOMFAC-web.exe /O3 Mod_kind_param.f90 ModStringFunctions.f90 ModSystemProp.f90 Mod_MINPACK.f90 ModSubgroupProp.f90 ModCompScaleConversion.f90 ModSRparam.f90 ModAIOMFACvar.f90 ModMRpart.f90 ModOScommands.f90 ModPureCompProp.f90 ModComponentNames.f90 ModNumericalTransformations.f90 Mod_InputOutput.f90 ModViscEyring.f90 ModPureViscosPar.f90 ModSRunifac.f90 SubModDefSystem.f90 ModCalcActCoeff.f90 ModZSRvisc.f90 SubModDissociationEquil.f90 ModFiniteDiffSens.f90 zerobracket_inwards.f90 brent.f90 AIOMFAC_inout.f90 Main_IO_driver.f90
	```
 	- the generated executable file named `AIOMFAC-web.exe` will be placed into the Fortran code folder.

- On [Linux]:
	- open a terminal and make sure that a recent version of gfortran is available (check with `gfortran --version`);
 	- navigate to your local AIOMFAC Fortran source code directory;
	- copy & paste, then execute the following command line (for compilation with -O3 optimization set when using gfortran, example works for AIOMFAC-web v3.14 and later):
   	```
	gfortran -o AIOMFAC-web.out -O3 -ffree-line-length-none -fstack-protector-strong -fbounds-check Mod_kind_param.f90  ModStringFunctions.f90 ModSystemProp.f90 Mod_MINPACK.f90 ModSubgroupProp.f90 ModCompScaleConversion.f90 ModSRparam.f90 ModAIOMFACvar.f90 ModMRpart.f90 ModOScommands.f90 ModPureCompProp.f90 ModComponentNames.f90 ModNumericalTransformations.f90 Mod_InputOutput.f90 ModViscEyring.f90 ModPureViscosPar.f90 ModSRunifac.f90 SubModDefSystem.f90 ModCalcActCoeff.f90 ModZSRvisc.f90 SubModDissociationEquil.f90 ModFiniteDiffSens.f90 zerobracket_inwards.f90 brent.f90 AIOMFAC_inout.f90 Main_IO_driver.f90
	```
    - the generated executable file named AIOMFAC-web.out will be placed into the Fortran code folder.
    - The `build_command_line.txt` file includes alternative command lines for debug-mode compilation as well as information on how to activate a recent gfortran version on RedHat and CentOS Linux.

- Alternatively, on [Linux] one can use the included makefile to build the code (on command line, navigate to the FortranCode folder and enter `make`). You could also re-generate a makefile by running the attached Perl script `mkmf.pl` (developed by V. Balaji, v.balaji@noaa.gov); that requires Perl to be installed and available from command line. I slightly modified an older version of the make make `mkmf` application to enable Fortran submodules to help establishing the correct dependencies of modules, submodules and subroutines/functions; see also [information here](https://github.com/NOAA-GFDL/mkmf/tree/main) and read the instructions provided in the `maketarget_commands_info_mkmf_Perl.txt` file included under `FortranCode`.
   
- Moreover, for in-depth code editing, debugging and development purposes, on [Windows] I recommend using MS Visual Studio Community with Intel's oneAPI Fortran compiler integration.

### (5) Test the Fortran program
-To be added...


----
## Quick guide to running AIOMFAC from a command prompt
To run the AIOMFAC program for your own system of components, this is a relatively straightforward task. The following inputs need to be provided.

-To be added...
