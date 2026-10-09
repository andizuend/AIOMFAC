# Makefile created by mkmf.pl 2022.01.00


MK_TEMPLATE = maketemplate.txt
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: AIOMFAC-web.out
AIOMFAC_inout.o: ./AIOMFAC_inout.f90 Mod_kind_param.o ModSystemProp.o ModAIOMFACvar.o ModCompScaleConversion.o ModCalcActCoeff.o ModFiniteDiffSens.o ModSubgroupProp.o ModPureViscosPar.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./AIOMFAC_inout.f90
Main_IO_driver.o: ./Main_IO_driver.f90 Mod_kind_param.o ModSystemProp.o ModSubgroupProp.o ModMRpart.o ModSRparam.o Mod_InputOutput.o ModPureCompProp.o ModComponentNames.o ModOScommands.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./Main_IO_driver.f90
ModAIOMFACvar.o: ./ModAIOMFACvar.f90 Mod_kind_param.o ModSystemProp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModAIOMFACvar.f90
ModCalcActCoeff.o: ./ModCalcActCoeff.f90 Mod_kind_param.o ModAIOMFACvar.o ModSystemProp.o ModCompScaleConversion.o ModNumericalTransformations.o ModSRunifac.o ModMRpart.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModCalcActCoeff.f90
ModCompScaleConversion.o: ./ModCompScaleConversion.f90 Mod_kind_param.o ModSystemProp.o ModNumericalTransformations.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModCompScaleConversion.f90
ModComponentNames.o: ./ModComponentNames.f90 ModSystemProp.o Mod_kind_param.o ModSubgroupProp.o ModStringFunctions.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModComponentNames.f90
ModFiniteDiffSens.o: ./ModFiniteDiffSens.f90 Mod_kind_param.o ModSystemProp.o ModAIOMFACvar.o ModCompScaleConversion.o ModCalcActCoeff.o ModViscEyring.o ModZSRvisc.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModFiniteDiffSens.f90
ModMRpart.o: ./ModMRpart.f90 Mod_kind_param.o ModSystemProp.o ModAIOMFACvar.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModMRpart.f90
ModNumericalTransformations.o: ./ModNumericalTransformations.f90 Mod_kind_param.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModNumericalTransformations.f90
ModOScommands.o: ./ModOScommands.f90 Mod_kind_param.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModOScommands.f90
ModPureCompProp.o: ./ModPureCompProp.f90 Mod_kind_param.o ModSystemProp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ModPureCompProp.f90
ModPureViscosPar.o: ./ModPureViscosPar.f90 Mod_kind_param.o ModSubgroupProp.o ModOScommands.o ModPureCompProp.o ModSystemProp.o Mod_InputOutput.o
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
Mod_InputOutput.o: ./Mod_InputOutput.f90 Mod_kind_param.o ModSystemProp.o ModComponentNames.o ModSubgroupProp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./Mod_InputOutput.f90
Mod_MINPACK.o: ./Mod_MINPACK.f90 Mod_kind_param.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./Mod_MINPACK.f90
Mod_NumPrec.o: ./Mod_NumPrec.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./Mod_NumPrec.f90
Mod_kind_param.o: ./Mod_kind_param.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./Mod_kind_param.f90
SubModDefSystem.o: ./SubModDefSystem.f90 ModSystemProp.o Mod_InputOutput.o ModSubgroupProp.o ModAIOMFACvar.o ModMRpart.o ModSRunifac.o ModComponentNames.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./SubModDefSystem.f90
SubModDissociationEquil.o: ./SubModDissociationEquil.f90 ModCalcActCoeff.o ModSystemProp.o ModAIOMFACvar.o ModCompScaleConversion.o Mod_MINPACK.o ModSubgroupProp.o ModNumericalTransformations.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./SubModDissociationEquil.f90
aiomfac_inout__genmod.o: ./aiomfac_inout__genmod.f90 ModSystemProp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./aiomfac_inout__genmod.f90
brent.o: ./brent.f90 Mod_kind_param.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./brent.f90
brentzero__genmod.o: ./brentzero__genmod.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./brentzero__genmod.f90
glomin__genmod.o: ./glomin__genmod.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./glomin__genmod.f90
local_min__genmod.o: ./local_min__genmod.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./local_min__genmod.f90
local_min_rc__genmod.o: ./local_min_rc__genmod.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./local_min_rc__genmod.f90
timestampbrent__genmod.o: ./timestampbrent__genmod.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./timestampbrent__genmod.f90
zero_rc__genmod.o: ./zero_rc__genmod.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./zero_rc__genmod.f90
zerobracket_inwards.o: ./zerobracket_inwards.f90 Mod_kind_param.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./zerobracket_inwards.f90
zerobracket_inwards__genmod.o: ./zerobracket_inwards__genmod.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./zerobracket_inwards__genmod.f90
SRC = ./ModSubgroupProp.f90 ./ModPureCompProp.f90 ./zerobracket_inwards__genmod.f90 ./ModCalcActCoeff.f90 ./ModStringFunctions.f90 ./ModSRunifac.f90 ./ModViscEyring.f90 ./ModPureViscosPar.f90 ./ModSystemProp.f90 ./ModZSRvisc.f90 ./ModCompScaleConversion.f90 ./brentzero__genmod.f90 ./ModMRpart.f90 ./glomin__genmod.f90 ./ModComponentNames.f90 ./brent.f90 ./SubModDefSystem.f90 ./Mod_InputOutput.f90 ./SubModDissociationEquil.f90 ./Mod_NumPrec.f90 ./Mod_MINPACK.f90 ./ModNumericalTransformations.f90 ./ModAIOMFACvar.f90 ./zero_rc__genmod.f90 ./local_min_rc__genmod.f90 ./Main_IO_driver.f90 ./AIOMFAC_inout.f90 ./timestampbrent__genmod.f90 ./ModFiniteDiffSens.f90 ./ModOScommands.f90 ./local_min__genmod.f90 ./Mod_kind_param.f90 ./aiomfac_inout__genmod.f90 ./zerobracket_inwards.f90 ./ModSRparam.f90
OBJ = ModSubgroupProp.o ModPureCompProp.o zerobracket_inwards__genmod.o ModCalcActCoeff.o ModStringFunctions.o ModSRunifac.o ModViscEyring.o ModPureViscosPar.o ModSystemProp.o ModZSRvisc.o ModCompScaleConversion.o brentzero__genmod.o ModMRpart.o glomin__genmod.o ModComponentNames.o brent.o SubModDefSystem.o Mod_InputOutput.o SubModDissociationEquil.o Mod_NumPrec.o Mod_MINPACK.o ModNumericalTransformations.o ModAIOMFACvar.o zero_rc__genmod.o local_min_rc__genmod.o Main_IO_driver.o AIOMFAC_inout.o timestampbrent__genmod.o ModFiniteDiffSens.o ModOScommands.o local_min__genmod.o Mod_kind_param.o aiomfac_inout__genmod.o zerobracket_inwards.o ModSRparam.o
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

