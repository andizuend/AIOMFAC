        !COMPILER-GENERATED INTERFACE MODULE: Tue Aug 25 15:47:28 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZEROBRACKET_INWARDS__genmod
          INTERFACE 
            SUBROUTINE ZEROBRACKET_INWARDS(FX,XLOWER,XUPPER,NTRY,       &
     &GEOMSCAL,NB,XBLOW,XBUP,SUCCESS)
              INTEGER(KIND=4), INTENT(IN) :: NTRY
              REAL(KIND=8) :: FX
              EXTERNAL FX
              REAL(KIND=8), INTENT(INOUT) :: XLOWER
              REAL(KIND=8), INTENT(INOUT) :: XUPPER
              LOGICAL(KIND=4), INTENT(IN) :: GEOMSCAL
              INTEGER(KIND=4), INTENT(OUT) :: NB
              REAL(KIND=8), INTENT(OUT) :: XBLOW(1:NTRY+1)
              REAL(KIND=8), INTENT(OUT) :: XBUP(1:NTRY+1)
              LOGICAL(KIND=4), INTENT(OUT) :: SUCCESS
            END SUBROUTINE ZEROBRACKET_INWARDS
          END INTERFACE 
        END MODULE ZEROBRACKET_INWARDS__genmod
