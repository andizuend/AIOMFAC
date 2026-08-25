        !COMPILER-GENERATED INTERFACE MODULE: Tue Aug 25 15:47:28 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AIOMFAC_INOUT__genmod
          INTERFACE 
            SUBROUTINE AIOMFAC_INOUT(INPUTCONC,XINPUTTYPE,TKELVIN,      &
     &NSPECIES,OUTPUTVARS,OUTPUTVISCVARS,OUTNAMES,ERRORFLAG_LIST,       &
     &WARNINGFLAG)
              USE MODSYSTEMPROP
              REAL(KIND=8), INTENT(IN) :: INPUTCONC(NINDCOMP)
              LOGICAL(KIND=4), INTENT(IN) :: XINPUTTYPE
              REAL(KIND=8), INTENT(IN) :: TKELVIN
              INTEGER(KIND=4), INTENT(OUT) :: NSPECIES
              REAL(KIND=8), INTENT(OUT) :: OUTPUTVARS(6,NKNPNGS)
              REAL(KIND=8), INTENT(OUT) :: OUTPUTVISCVARS(2)
              CHARACTER(*), INTENT(OUT) :: OUTNAMES(NKNPNGS)
              LOGICAL(KIND=4), INTENT(OUT) :: ERRORFLAG_LIST(50)
              INTEGER(KIND=4), INTENT(OUT) :: WARNINGFLAG
            END SUBROUTINE AIOMFAC_INOUT
          END INTERFACE 
        END MODULE AIOMFAC_INOUT__genmod
