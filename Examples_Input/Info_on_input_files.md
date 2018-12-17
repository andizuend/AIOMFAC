Running the AIOMFAC-web program from your own system with custom input requires an input file. 
You can run your compiled and built executable (from the Fortran code) with the name of the input file as a command line argument. Having  the executable, e.g. <code>AIOMFAC-web.out</code>, in a folder of your choice, you will need to create two more directories within that parent directory, namely folder <code>Inputfiles</code> and folder <code>Outputfiles</code> - these two folders will not be created automatically. Then, depending on your operating system, you can run a script to run the AIOMFAC program (such as the attached Run_Model_batchfile.bat on Windows; edit that file in notepad (or your text-editor of choice) to add your input file name) or by running it from a command terminal. For example, from a Linux terminal, type: <code>./AIOMFAC-web.out ./Inputfiles/input_0008.txt </code>. The model output will be in a corresponding text-file in the folder "Outputfiles". 

Examples for input files are included in this folder. For more information, please also read the hints on https://aiomfac.lab.mcgill.ca/help.html.

You can also use the input form on the AIOMFAC website (https://aiomfac.lab.mcgill.ca/model.html) to create valid input files, which can then be downloaded via the text-file link on the "Results" page of the website.

Information about the identification codes of valid AIOMFAC subgroups is provided via the spreadsheet in the folder AIOMFAC_Documents on this site.
