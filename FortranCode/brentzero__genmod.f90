        !COMPILER-GENERATED INTERFACE MODULE: Tue Aug 25 15:47:28 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BRENTZERO__genmod
          INTERFACE 
            FUNCTION BRENTZERO(A,B,MACHEP,T,F)
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: MACHEP
              REAL(KIND=8) :: T
              REAL(KIND=8) :: F
              EXTERNAL F
              REAL(KIND=8) :: BRENTZERO
            END FUNCTION BRENTZERO
          END INTERFACE 
        END MODULE BRENTZERO__genmod
