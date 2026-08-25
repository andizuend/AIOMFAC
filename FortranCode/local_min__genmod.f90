        !COMPILER-GENERATED INTERFACE MODULE: Tue Aug 25 15:47:28 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LOCAL_MIN__genmod
          INTERFACE 
            FUNCTION LOCAL_MIN(A,B,EPS,T,F,X)
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: EPS
              REAL(KIND=8) :: T
              REAL(KIND=8) :: F
              EXTERNAL F
              REAL(KIND=8) :: X
              REAL(KIND=8) :: LOCAL_MIN
            END FUNCTION LOCAL_MIN
          END INTERFACE 
        END MODULE LOCAL_MIN__genmod
