## About this code folder
This folder contains all Fortran files (.f90) needed to compile and link the modular AIOMFAC Fortran program. 
For a list of releases and associated descriptions, including the latest release of the Fortran code, see https://github.com/andizuend/AIOMFAC/releases.
### Compiling and linking
The file [build_command_line.txt](./build_command_line.txt) provides information about how to compile and link the source code with `gfortran` or `ifx` (specifically for v3.00 and more recent).

Alternatively, you can compile and link using your compiler of choice (a recent Fortran compiler is necessary; one with support of submodules, e.g. `gfortran` version 9 or newer). On Linux, e.g. with gfortran, using the attached Makefile - or by building a new makefile using the attached Perl script mkmf.pl (developed by V. Balaji (v.balaji@noaa.gov). I slightly modified an older version of the mkmf application to deal to enable Fortran submodules, see also https://www.gfdl.noaa.gov/~vb/mkmf.html) and the [maketarget_commands_info_mkmf_Perl.txt](./maketarget_commands_info_mkmf_Perl.txt) file to help building the correct dependencies of modules, submodules and subroutines/functions from command line).
On MS Windows, I recommend using either (i) MS Visual Studio Community and Intel's oneAPI Fortran compiler integration (which is free software) or (ii) an Intel oneAPI command prompt for direct compilation from the terminal or (iii) the use of Windows subsystem for Linux (WSL) -- in the latter case follow then the instructions for building the program described for Linux.

## Main program
File <code>Main_IO_driver.f90</code> is the entry point to the program and likely the only file you may wish to adjust if you would like to modify the AIOMFAC interface for your purposes.
