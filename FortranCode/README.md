## About this code folder
This folder contains all Fortran files (.f90) needed to compile and link the modular AIOMFAC Fortran program. 
### Compiling and linking
You can compile and link using your compiler of choice (a recent compiler is necessary; one with support of submodules, e.g. gfortran version 7 or newer). On Linux, e.g. with gfortran, using the attached Makefile - or by building a new makefile using the attached Perl script mkmf.pl (developed by V. Balaji (v.balaji@noaa.gov), which I slightly modified to deal with the relatively recent feature of Fortran submodules, see also http://www.gfdl.noaa.gov/~vb/mkmf.html) and command files to help building the correct dependencies of modules, submodules and subroutines/functions). On MS Windows, I recommend using MS Visual Studio with the Intel Fortran Composer (Parallel studio) integration (if you have it available - it is not free, except in special cases, such as for certain students).

## Main program
File <code>Main_IO_driver.f90</code> is the entry point to the program and likely the only file you may wish to adjust if you would like to modify the AIOMFAC interface for your purposes.
