@ECHO OFF
SETLOCAL

ECHO  ************************************************************
ECHO *                                                            *
ECHO * This is a Windows batch file to run an executable (here    *
ECHO * AIOMFAC-web.exe) with an input file located in the stated  *
ECHO * folder given by the relative path.                         *
ECHO *                                      (script by A. Zuend)  *
ECHO *                                                            *
ECHO  ************************************************************

SET _input=.\Inputfiles\input_0008.txt

:: run the executable with the chosen input file:

AIOMFAC-web.exe %_input%

::------------------------------ The End ------------------------------------------------


ECHO -------------------------------------
ECHO The executable batch job is finished. 
ECHO -------------------------------------
ECHO.

pause