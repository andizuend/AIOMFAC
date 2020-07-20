# Makefile created by mkmf.pl $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include maketemplate.txt


.DEFAULT:
	-echo $@ does not exist.
all: AIOMFAC-web.out
AIOMFAC_inout.o: ./AIOMFAC_inout.f90 ModSystemProp.o ModAIOMFACvar.o ModCompScaleConversion.o ModCalcActCoeff.o ModSubgroupProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./AIOMFAC_inout.f90
Main_IO_driver.o: ./Main_IO_driver.f90 ModSystemProp.o ModSubgroupProp.o ModMRpart.o ModSRunifac.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./Main_IO_driver.f90
ModAIOMFACvar.o: ./ModAIOMFACvar.f90 ModSystemProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModAIOMFACvar.f90
ModCalcActCoeff.o: ./ModCalcActCoeff.f90 ModAIOMFACvar.o ModSystemProp.o ModCompScaleConversion.o ModSRunifac.o ModMRpart.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModCalcActCoeff.f90
ModCompScaleConversion.o: ./ModCompScaleConversion.f90 ModSystemProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModCompScaleConversion.f90
ModComponentNames.o: ./ModComponentNames.f90 ModSystemProp.o ModSubgroupProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModComponentNames.f90
ModMRpart.o: ./ModMRpart.f90 ModSystemProp.o ModAIOMFACvar.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModMRpart.f90
ModSRunifac.o: ./ModSRunifac.f90 ModAIOMFACvar.o ModSystemProp.o Mod_PureViscosPar.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModSRunifac.f90
ModSubgroupProp.o: ./ModSubgroupProp.f90 ModSystemProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModSubgroupProp.f90
ModSystemProp.o: ./ModSystemProp.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModSystemProp.f90
Mod_PureViscosPar.o: ./Mod_PureViscosPar.f90 ModSubgroupProp.o ModSystemProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./Mod_PureViscosPar.f90
OutputHTML.o: ./OutputHTML.f90 ModSystemProp.o ModSubgroupProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./OutputHTML.f90
OutputTXT.o: ./OutputTXT.f90 ModSystemProp.o ModSubgroupProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./OutputTXT.f90
ReadInputFile.o: ./ReadInputFile.f90 ModSystemProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ReadInputFile.f90
RepErrorWarning.o: ./RepErrorWarning.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./RepErrorWarning.f90
SubModDefSystem.o: ./SubModDefSystem.f90 ModSystemProp.o ModSubgroupProp.o ModAIOMFACvar.o ModMRpart.o ModSRunifac.o ModComponentNames.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./SubModDefSystem.f90
brent.o: ./brent.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./brent.f90
zerobracket_inwards.o: ./zerobracket_inwards.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./zerobracket_inwards.f90
SRC = ./zerobracket_inwards.f90 ./ModCalcActCoeff.f90 ./AIOMFAC_inout.f90 ./RepErrorWarning.f90 ./ModAIOMFACvar.f90 ./ModMRpart.f90 ./ModSubgroupProp.f90 ./ReadInputFile.f90 ./OutputHTML.f90 ./ModSystemProp.f90 ./ModCompScaleConversion.f90 ./brent.f90 ./SubModDefSystem.f90 ./ModSRunifac.f90 ./OutputTXT.f90 ./Main_IO_driver.f90 ./Mod_PureViscosPar.f90 ./ModComponentNames.f90
OBJ = zerobracket_inwards.o ModCalcActCoeff.o AIOMFAC_inout.o RepErrorWarning.o ModAIOMFACvar.o ModMRpart.o ModSubgroupProp.o ReadInputFile.o OutputHTML.o ModSystemProp.o ModCompScaleConversion.o brent.o SubModDefSystem.o ModSRunifac.o OutputTXT.o Main_IO_driver.o Mod_PureViscosPar.o ModComponentNames.o
clean: neat
	-rm -f .AIOMFAC-web.out.cppdefs $(OBJ) AIOMFAC-web.out
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
AIOMFAC-web.out: $(OBJ) 
	$(LD) $(OBJ) -o AIOMFAC-web.out  $(LDFLAGS)
