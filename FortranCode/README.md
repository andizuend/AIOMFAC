This folder contains all Fortran files (.f90) needed to compile and link the modular AIOMFAC Fortran program. 
You can compile and link using your compiler of choice. On Linux, e.g. with gfortran, using the attaced Makefile - or building a new makefile using the attached Perl script mkmf.pl and command files to help building the correct dependencies of modules, submodules and subroutines/functions). On MS Windows, I recommend using MS Visual Studio with the Intel Fortran Composer (Parallel studio) integration (if you have it available - it is not free, except in special cases).

File <code>Main_IO_driver.f90</code> is the entry point to the program and likely the only file you may wish to adjust if you would like to modify the AIOMFAC interface for your purposes.
