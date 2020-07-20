## About this code folder
This folder contains all Fortran files (.f90) needed to compile and link the modular AIOMFAC Fortran program. 
For a list of releases, inlcuding the latest release of the Fortran code, see https://github.com/andizuend/AIOMFAC/releases.
### Compiling and linking
You can compile and link using your compiler of choice (a recent compiler is necessary; one with support of submodules, e.g. `gfortran` version 7 or newer). On Linux, e.g. with gfortran, using the attached Makefile - or by building a new makefile using the attached Perl script mkmf.pl (developed by V. Balaji (v.balaji@noaa.gov), which I slightly modified to deal with the relatively recent feature of Fortran submodules, see also http://www.gfdl.noaa.gov/~vb/mkmf.html) and the [maketarget_commands_info_mkmf_Perl.txt](./maketarget_commands_info_mkmf_Perl.txt) file to help building the correct dependencies of modules, submodules and subroutines/functions from command line). 
On MS Windows, I recommend using MS Visual Studio with the Intel Fortran Composer (Parallel studio; `ifort`) integration (if you have it available - it is not free, except in special cases, such as for qualified students).

## Main program
File <code>Main_IO_driver.f90</code> is the entry point to the program and likely the only file you may wish to adjust if you would like to modify the AIOMFAC interface for your purposes.
