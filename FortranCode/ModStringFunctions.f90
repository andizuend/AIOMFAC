! -----------------------------------------------
MODULE ModStringFunctions   ! by David Frank  dave_frank@hotmail.com
                            ! http://home.earthlink.net/~dave_gemini/strings.f90
                            ! with a few minor modifications (PURE functions, INTENT(IN),...) by Andi Zuend (McGill U.)
IMPLICIT NONE   

PUBLIC
! Copy (generic) char array to string or string to char array
! Clen           returns same as LEN      unless last non-blank char = null
! Clen_trim      returns same as LEN_TRIM    "              "
! Ctrim          returns same as TRIM        "              "
! Count_Items    in string that are blank or comma separated
! Reduce_Blanks  in string to 1 blank between items, last char not blank
! Replace_Text   in all occurances in string with replacement string
! Spack          pack string's chars == extract string's chars
! Tally          occurances in string of text arg
! Translate      text arg via indexed code table
! Upper/Lower    case the text arg
! Replace_Text_Advance ::   replace a string with replacement text, then advance, such that the 
!                           replacement text itself is not recursively processed. (added by A. Zuend)

INTERFACE Copy    ! generic
    MODULE PROCEDURE copy_a2s, copy_s2a
END INTERFACE Copy

    CONTAINS
    
    ! ------------------------
    PURE FUNCTION Copy_a2s(a)  RESULT (s)    ! copy char array to string
    CHARACTER,INTENT(IN) :: a(:)
    CHARACTER(SIZE(a)) :: s
    INTEGER(4) :: i
    DO i = 1,SIZE(a)
        s(i:i) = a(i)
    END DO
    END FUNCTION Copy_a2s

    ! ------------------------
    PURE FUNCTION Copy_s2a(s)  RESULT (a)   ! copy s(1:Clen(s)) to char array
    CHARACTER(*),INTENT(IN) :: s
    CHARACTER :: a(LEN(s))
    INTEGER(4) :: i
    DO i = 1,LEN(s)
        a(i) = s(i:i)
    END DO
    END FUNCTION Copy_s2a

    ! ------------------------
    PURE INTEGER(4) FUNCTION Clen(s)    ! returns same result as LEN unless:
    CHARACTER(*),INTENT(IN) :: s        ! last non-blank char is null
    INTEGER(4) :: i
    Clen = LEN(s)
    i = LEN_TRIM(s)
    IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
    END FUNCTION Clen

    ! ------------------------
    PURE INTEGER(4) FUNCTION Clen_trim(s)   ! returns same result as LEN_TRIM unless:
    CHARACTER(*),INTENT(IN) :: s            ! last char non-blank is null, if true:
    INTEGER(4) :: i                         ! then len of C string is returned, note:
    ! Ctrim is only user of this function
    i = LEN_TRIM(s) ; Clen_trim = i
    IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
    END FUNCTION Clen_trim

    ! ----------------
    PURE FUNCTION Ctrim(s1)  RESULT(s2)     ! returns same result as TRIM unless:
    CHARACTER(*),INTENT(IN)  :: s1     ! last non-blank char is null in which
    CHARACTER(Clen_trim(s1)) :: s2     ! case trailing blanks prior to null
    s2 = s1                            ! are output
    END FUNCTION Ctrim

    ! --------------------
    PURE INTEGER(4) FUNCTION Count_Items(s1)  ! in string or C string that are blank or comma separated
    CHARACTER(*),INTENT(IN) :: s1
    CHARACTER(Clen(s1)) :: s
    INTEGER(4) :: i, k

    s = s1                            ! remove possible last char null
    k = 0  ; IF (s /= ' ') k = 1      ! string has at least 1 item
    DO i = 1,LEN_TRIM(s)-1
        IF (s(i:i) /= ' '.AND.s(i:i) /= ',' &
            .AND. s(i+1:i+1) == ' '.OR. s(i+1:i+1) == ',') k = k+1
    END DO
    Count_Items = k
    END FUNCTION Count_Items

    ! --------------------
    PURE FUNCTION Reduce_Blanks(s)  RESULT (outs)
    CHARACTER(*),INTENT(IN)      :: s
    CHARACTER(LEN_TRIM(s)) :: outs
    INTEGER(4)           :: i, k, n

    n = 0  ; k = LEN_TRIM(s)          ! k=index last non-blank (may be null)
    DO i = 1,k-1                      ! dont process last char yet
        n = n+1 ; outs(n:n) = s(i:i)
        IF (s(i:i+1) == '  ') n = n-1  ! backup/discard consecutive output blank
    END DO
    n = n+1  ; outs(n:n)  = s(k:k)    ! last non-blank char output (may be null)
    IF (n < k) outs(n+1:) = ' '       ! pad trailing blanks
    END FUNCTION Reduce_Blanks

    ! ------------------
    PURE FUNCTION Replace_Text (s, text, rep)  RESULT(outs)
    CHARACTER(*),INTENT(IN) :: s, text, rep
    CHARACTER(LEN(s)+100)   :: outs     ! provide outs with extra 100 char len
    INTEGER(4)              :: i, nt, nr
    !..........................
    outs = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
    DO
        i = INDEX(outs, text(:nt)) 
        IF (i == 0) EXIT
        outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
    END DO
    END FUNCTION Replace_Text
    
    
    ! ------------------
    PURE FUNCTION Replace_Text_Advance(s, text, rep)  RESULT(outs)
    CHARACTER(*),INTENT(IN) :: s, text, rep
    CHARACTER(LEN(s)+100)   :: outs     ! provide outs with extra 100 char len
    INTEGER(4)              :: i, nt, nr, k
    !..........................
    outs = s
    nt = LEN_TRIM(text) 
    nr = LEN_TRIM(rep)
    k = 1
    DO
        i = k -1 + INDEX( outs(k:), text(:nt) ) 
        IF (i == k -1) THEN
            EXIT
        ELSE
            outs = outs(:i-1)//rep(:nr)//outs(i+nt:)
            k = i + nr
        ENDIF
    END DO
    
    END FUNCTION Replace_Text_Advance

    
    ! ---------------------------------
    PURE FUNCTION Spack (s,ex)  RESULT (outs)
    CHARACTER(*),INTENT(IN) :: s,ex
    CHARACTER(LEN(s)) :: outs
    CHARACTER :: aex(LEN(ex))   ! array of ex chars to extract
    INTEGER(4)   :: i, n

    n = 0  ;  aex = Copy(ex)
    DO i = 1,LEN(s)
        IF (.NOT.ANY(s(i:i) == aex)) CYCLE   ! dont pack char
        n = n+1 ; outs(n:n) = s(i:i)
    END DO
    outs(n+1:) = ' '     ! pad with trailing blanks
    END FUNCTION Spack

    ! --------------------
    PURE INTEGER(4) FUNCTION Tally (s,text)
    CHARACTER(*),INTENT(IN) :: s, text
    INTEGER(4) :: i, nt

    Tally = 0 ; nt = LEN_TRIM(text)
    DO i = 1,LEN(s)-nt+1
        IF (s(i:i+nt-1) == text(:nt)) Tally = Tally+1
    END DO
    END FUNCTION Tally

    ! ---------------------------------
    PURE FUNCTION Translate(s1,codes)  RESULT (s2)
    CHARACTER(*),INTENT(IN) :: s1, codes(2)
    CHARACTER(LEN(s1)) :: s2
    CHARACTER          :: ch
    INTEGER(4)            :: i, j

    DO i = 1,LEN(s1)
        ch = s1(i:i)
        j = INDEX(codes(1),ch) ; IF (j > 0) ch = codes(2)(j:j)
        s2(i:i) = ch
    END DO
    END FUNCTION Translate

    ! ---------------------------------
    PURE FUNCTION Upper(s1)  RESULT (s2)
    CHARACTER(*),INTENT(IN) :: s1
    CHARACTER(LEN(s1)) :: s2
    CHARACTER          :: ch
    INTEGER(4),PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
    INTEGER(4)            :: i

    DO i = 1,LEN(s1)
        ch = s1(i:i)
        IF (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
        s2(i:i) = ch
    END DO
    END FUNCTION Upper

    ! ---------------------------------
    PURE FUNCTION Lower(s1)  RESULT (s2)
    CHARACTER(*),INTENT(IN) :: s1
    CHARACTER(LEN(s1)) :: s2
    CHARACTER          :: ch
    INTEGER(4),PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
    INTEGER(4)            :: i

    DO i = 1,LEN(s1)
        ch = s1(i:i)
        IF (ch >= 'A'.AND.ch <= 'Z') ch = CHAR(ICHAR(ch)-DUC)
        s2(i:i) = ch
    END DO
    END FUNCTION Lower

END MODULE ModStringFunctions