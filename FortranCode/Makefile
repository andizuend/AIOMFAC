# Makefile created by mkmf.pl 2022.01.00


MK_TEMPLATE = maketemplate.txt
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: AIOMFAC-web.out
AIOMFAC_inout.o: ./AIOMFAC_inout.f90 Mod_kind_param.o ModSystemProp.o ModAIOMFACvar.o ModCompScaleConversion.o ModCalcActCoeff.o ModFiniteDiffSens.o ModSubgroupProp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./AIOMFAC_inout.f90
Main_IO_driver.o: ./Main_IO_driver.f90 Mod_kind_param.o ModSystemProp.o ModSubgroupProp.o ModMRpart.o ModSRparam.o Mod_InputOutput.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./Main_IO_driver.f90
ModAIOMFACvar.o: ./ModAIOMFACvar.f90 Mod_kind_param.o ModSystemProp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModAIOMFACvar.f90
ModCalcActCoeff.o: ./ModCalcActCoeff.f90 Mod_kind_param.o ModAIOMFACvar.o ModSystemProp.o ModCompScaleConversion.o ModNumericalTransformations.o ModSRunifac.o ModMRpart.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModCalcActCoeff.f90
ModCompScaleConversion.o: ./ModCompScaleConversion.f90 Mod_kind_param.o ModSystemProp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModCompScaleConversion.f90
ModComponentNames.o: ./ModComponentNames.f90 Mod_kind_param.o ModSystemProp.o ModSubgroupProp.o ModStringFunctions.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModComponentNames.f90
ModFiniteDiffSens.o: ./ModFiniteDiffSens.f90 Mod_kind_param.o ModSystemProp.o ModAIOMFACvar.o ModCompScaleConversion.o ModCalcActCoeff.o ModViscEyring.o ModZSRvisc.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModFiniteDiffSens.f90
ModMRpart.o: ./ModMRpart.f90 Mod_kind_param.o ModSystemProp.o ModAIOMFACvar.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModMRpart.f90
ModNumericalTransformations.o: ./ModNumericalTransformations.f90 Mod_kind_param.o ModSystemProp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModNumericalTransformations.f90
ModPureViscosPar.o: ./ModPureViscosPar.f90 Mod_kind_param.o ModSubgroupProp.o ModSystemProp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModPureViscosPar.f90
ModSRparam.o: ./ModSRparam.f90 Mod_kind_param.o ModSystemProp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModSRparam.f90
ModSRunifac.o: ./ModSRunifac.f90 Mod_kind_param.o ModAIOMFACvar.o ModSystemProp.o ModSRparam.o ModPureViscosPar.o ModViscEyring.o ModCompScaleConversion.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModSRunifac.f90
ModStringFunctions.o: ./ModStringFunctions.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModStringFunctions.f90
ModSubgroupProp.o: ./ModSubgroupProp.f90 Mod_kind_param.o ModSystemProp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModSubgroupProp.f90
ModSystemProp.o: ./ModSystemProp.f90 Mod_kind_param.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModSystemProp.f90
ModViscEyring.o: ./ModViscEyring.f90 Mod_kind_param.o ModAIOMFACvar.o ModSystemProp.o ModSRparam.o ModCompScaleConversion.o ModNumericalTransformations.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModViscEyring.f90
ModZSRvisc.o: ./ModZSRvisc.f90 Mod_kind_param.o ModSystemProp.o ModAIOMFACvar.o ModCalcActCoeff.o Mod_MINPACK.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModZSRvisc.f90
Mod_InputOutput.o: ./Mod_InputOutput.f90 Mod_kind_param.o ModSystemProp.o ModSubgroupProp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./Mod_InputOutput.f90
Mod_MINPACK.o: ./Mod_MINPACK.f90 Mod_kind_param.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./Mod_MINPACK.f90
Mod_kind_param.o: ./Mod_kind_param.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./Mod_kind_param.f90
SubModDefSystem.o: ./SubModDefSystem.f90 ModSystemProp.o ModSubgroupProp.o ModAIOMFACvar.o ModMRpart.o ModSRunifac.o ModComponentNames.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./SubModDefSystem.f90
SubModDissociationEquil.o: ./SubModDissociationEquil.f90 ModCalcActCoeff.o ModSystemProp.o ModAIOMFACvar.o ModCompScaleConversion.o Mod_MINPACK.o ModSubgroupProp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./SubModDissociationEquil.f90
brent.o: ./brent.f90 Mod_kind_param.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./brent.f90
zerobracket_inwards.o: ./zerobracket_inwards.f90 Mod_kind_param.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./zerobracket_inwards.f90
SRC = ./ModSRparam.f90 ./SubModDefSystem.f90 ./Mod_kind_param.f90 ./ModPureViscosPar.f90 ./ModSRunifac.f90 ./Main_IO_driver.f90 ./ModSubgroupProp.f90 ./ModCompScaleConversion.f90 ./ModSystemProp.f90 ./ModViscEyring.f90 ./ModAIOMFACvar.f90 ./SubModDissociationEquil.f90 ./Mod_InputOutput.f90 ./Mod_MINPACK.f90 ./ModNumericalTransformations.f90 ./ModComponentNames.f90 ./AIOMFAC_inout.f90 ./brent.f90 ./ModFiniteDiffSens.f90 ./ModStringFunctions.f90 ./ModZSRvisc.f90 ./ModMRpart.f90 ./ModCalcActCoeff.f90 ./zerobracket_inwards.f90
OBJ = ModSRparam.o SubModDefSystem.o Mod_kind_param.o ModPureViscosPar.o ModSRunifac.o Main_IO_driver.o ModSubgroupProp.o ModCompScaleConversion.o ModSystemProp.o ModViscEyring.o ModAIOMFACvar.o SubModDissociationEquil.o Mod_InputOutput.o Mod_MINPACK.o ModNumericalTransformations.o ModComponentNames.o AIOMFAC_inout.o brent.o ModFiniteDiffSens.o ModStringFunctions.o ModZSRvisc.o ModMRpart.o ModCalcActCoeff.o zerobracket_inwards.o
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

