# Makefile created by mkmf.pl $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include maketemplate.txt


.DEFAULT:
	-echo $@ does not exist.
all: AIOMFAC-web.out
AIOMFAC_inout.o: ./AIOMFAC_inout.f90 Mod_NumPrec.o ModSystemProp.o ModAIOMFACvar.o ModCompScaleConversion.o ModCalcActCoeff.o ModFiniteDiffSens.o ModSubgroupProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./AIOMFAC_inout.f90
Main_IO_driver.o: ./Main_IO_driver.f90 Mod_NumPrec.o ModSystemProp.o ModSubgroupProp.o ModMRpart.o ModSRparam.o Mod_InputOutput.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./Main_IO_driver.f90
ModAIOMFACvar.o: ./ModAIOMFACvar.f90 ModSystemProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModAIOMFACvar.f90
ModCalcActCoeff.o: ./ModCalcActCoeff.f90 ModAIOMFACvar.o ModSystemProp.o ModCompScaleConversion.o ModNumericalTransformations.o ModSRunifac.o ModMRpart.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModCalcActCoeff.f90
ModCompScaleConversion.o: ./ModCompScaleConversion.f90 ModSystemProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModCompScaleConversion.f90
ModComponentNames.o: ./ModComponentNames.f90 ModSystemProp.o ModSubgroupProp.o ModStringFunctions.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModComponentNames.f90
ModFiniteDiffSens.o: ./ModFiniteDiffSens.f90 ModSystemProp.o ModAIOMFACvar.o ModCompScaleConversion.o ModCalcActCoeff.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModFiniteDiffSens.f90
ModMRpart.o: ./ModMRpart.f90 ModSystemProp.o ModAIOMFACvar.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModMRpart.f90
ModNumericalTransformations.o: ./ModNumericalTransformations.f90 ModSystemProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModNumericalTransformations.f90
ModPureViscosPar.o: ./ModPureViscosPar.f90 ModSubgroupProp.o ModSystemProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModPureViscosPar.f90
ModSRparam.o: ./ModSRparam.f90 ModSystemProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModSRparam.f90
ModSRunifac.o: ./ModSRunifac.f90 ModAIOMFACvar.o ModSystemProp.o ModSRparam.o ModPureViscosPar.o ModViscEyring.o ModCompScaleConversion.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModSRunifac.f90
ModStringFunctions.o: ./ModStringFunctions.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModStringFunctions.f90
ModSubgroupProp.o: ./ModSubgroupProp.f90 ModSystemProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModSubgroupProp.f90
ModSystemProp.o: ./ModSystemProp.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModSystemProp.f90
ModViscEyring.o: ./ModViscEyring.f90 ModAIOMFACvar.o ModSystemProp.o ModSRparam.o ModCompScaleConversion.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ModViscEyring.f90
Mod_InputOutput.o: ./Mod_InputOutput.f90 Mod_NumPrec.o ModSystemProp.o ModSubgroupProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./Mod_InputOutput.f90
Mod_MINPACK.o: ./Mod_MINPACK.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./Mod_MINPACK.f90
Mod_NumPrec.o: ./Mod_NumPrec.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./Mod_NumPrec.f90
SubModDefSystem.o: ./SubModDefSystem.f90 ModSystemProp.o ModSubgroupProp.o ModAIOMFACvar.o ModMRpart.o ModSRunifac.o ModComponentNames.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./SubModDefSystem.f90
SubModDissociationEquil.o: ./SubModDissociationEquil.f90 ModSystemProp.o ModCompScaleConversion.o ModAIOMFACvar.o Mod_MINPACK.o ModSubgroupProp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./SubModDissociationEquil.f90
brent.o: ./brent.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./brent.f90
zerobracket_inwards.o: ./zerobracket_inwards.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./zerobracket_inwards.f90
SRC = ./ModComponentNames.f90 ./Main_IO_driver.f90 ./SubModDefSystem.f90 ./ModNumericalTransformations.f90 ./ModSubgroupProp.f90 ./ModFiniteDiffSens.f90 ./Mod_MINPACK.f90 ./ModSRunifac.f90 ./brent.f90 ./Mod_NumPrec.f90 ./ModCalcActCoeff.f90 ./ModMRpart.f90 ./AIOMFAC_inout.f90 ./ModStringFunctions.f90 ./ModSystemProp.f90 ./ModViscEyring.f90 ./ModSRparam.f90 ./SubModDissociationEquil.f90 ./ModCompScaleConversion.f90 ./ModPureViscosPar.f90 ./Mod_InputOutput.f90 ./ModAIOMFACvar.f90 ./zerobracket_inwards.f90
OBJ = ModComponentNames.o Main_IO_driver.o SubModDefSystem.o ModNumericalTransformations.o ModSubgroupProp.o ModFiniteDiffSens.o Mod_MINPACK.o ModSRunifac.o brent.o Mod_NumPrec.o ModCalcActCoeff.o ModMRpart.o AIOMFAC_inout.o ModStringFunctions.o ModSystemProp.o ModViscEyring.o ModSRparam.o SubModDissociationEquil.o ModCompScaleConversion.o ModPureViscosPar.o Mod_InputOutput.o ModAIOMFACvar.o zerobracket_inwards.o
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
