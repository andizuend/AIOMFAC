! -----------------------------------------------
module ModStringFunctions   ! by David Frank  dave_frank@hotmail.com
                            ! http://home.earthlink.net/~dave_gemini/strings.f90
                            ! with a few minor modifications (pure functions, intent(in),...) by Andi Zuend (McGill U.)
implicit none   

public
! Copy (generic) char array to string or string to char array
! Clen           returns same as len      unless last non-blank char = null
! Clen_trim      returns same as len_trim    "              "
! Ctrim          returns same as trim        "              "
! Count_Items    in string that are blank or comma separated
! Reduce_Blanks  in string to 1 blank between items, last char not blank
! Replace_Text   in all occurances in string with replacement string
! Spack          pack string's chars == extract string's chars
! Tally          occurances in string of text arg
! Translate      text arg via indexed code table
! Upper/Lower    case the text arg
! Replace_Text_Advance ::   replace a string with replacement text, then advance, such that the 
!                           replacement text itself is not recursively processed. (added by A. Zuend)

interface Copy    ! generic
    module procedure copy_a2s, copy_s2a
end interface Copy

    contains
    
    ! ------------------------
    pure function Copy_a2s(a)  result (s)    ! copy char array to string
    character,intent(in) :: a(:)
    character(size(a)) :: s
    integer :: i
    do i = 1,size(a)
        s(i:i) = a(i)
    end do
    end function Copy_a2s

    ! ------------------------
    pure function Copy_s2a(s)  result (a)   ! copy s(1:Clen(s)) to char array
    character(*),intent(in) :: s
    character :: a(len(s))
    integer :: i
    do i = 1,len(s)
        a(i) = s(i:i)
    end do
    end function Copy_s2a

    ! ------------------------
    pure integer function Clen(s)    ! returns same result as len unless:
    character(*),intent(in) :: s        ! last non-blank char is null
    integer :: i
    Clen = len(s)
    i = len_trim(s)
    if (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
    end function Clen

    ! ------------------------
    pure integer function Clen_trim(s)   ! returns same result as len_trim unless:
    character(*),intent(in) :: s            ! last char non-blank is null, if true:
    integer :: i                         ! then len of C string is returned, note:
    ! Ctrim is only user of this function
    i = len_trim(s) ; Clen_trim = i
    if (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
    end function Clen_trim

    ! ----------------
    pure function Ctrim(s1)  result(s2)     ! returns same result as trim unless:
    character(*),intent(in)  :: s1     ! last non-blank char is null in which
    character(Clen_trim(s1)) :: s2     ! case trailing blanks prior to null
    s2 = s1                            ! are output
    end function Ctrim

    ! --------------------
    pure integer function Count_Items(s1)  ! in string or C string that are blank or comma separated
    character(*),intent(in) :: s1
    character(Clen(s1)) :: s
    integer :: i, k

    s = s1                            ! remove possible last char null
    k = 0  ; if (s /= ' ') k = 1      ! string has at least 1 item
    do i = 1,len_trim(s)-1
        if (s(i:i) /= ' '.AND.s(i:i) /= ',' &
            .AND. s(i+1:i+1) == ' '.OR. s(i+1:i+1) == ',') k = k+1
    end do
    Count_Items = k
    end function Count_Items

    ! --------------------
    pure function Reduce_Blanks(s)  result (outs)
    character(*),intent(in)      :: s
    character(len_trim(s)) :: outs
    integer           :: i, k, n

    n = 0  ; k = len_trim(s)          ! k=index last non-blank (may be null)
    do i = 1,k-1                      ! dont process last char yet
        n = n+1 ; outs(n:n) = s(i:i)
        if (s(i:i+1) == '  ') n = n-1  ! backup/discard consecutive output blank
    end do
    n = n+1  ; outs(n:n)  = s(k:k)    ! last non-blank char output (may be null)
    if (n < k) outs(n+1:) = ' '       ! pad trailing blanks
    end function Reduce_Blanks

    ! ------------------
    pure function Replace_Text (s, text, rep)  result(outs)
    character(*),intent(in) :: s, text, rep
    character(len(s)+100)   :: outs     ! provide outs with extra 100 char len
    integer              :: i, nt, nr
    !..........................
    outs = s ; nt = len_trim(text) ; nr = len_trim(rep)
    do
        i = index(outs, text(:nt)) 
        if (i == 0) exit
        outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
    end do
    end function Replace_Text
    
    
    ! ------------------
    pure function Replace_Text_Advance(s, text, rep)  result(outs)
    character(*),intent(in) :: s, text, rep
    character(len(s)+100)   :: outs     ! provide outs with extra 100 char len
    integer              :: i, nt, nr, k
    !..........................
    outs = s
    nt = len_trim(text) 
    nr = len_trim(rep)
    k = 1
    do
        i = k -1 + index( outs(k:), text(:nt) ) 
        if (i == k -1) then
            exit
        else
            outs = outs(:i-1)//rep(:nr)//outs(i+nt:)
            k = i + nr
        endif
    end do
    
    end function Replace_Text_Advance

    
    ! ---------------------------------
    pure function Spack (s,ex)  result (outs)
    character(*),intent(in) :: s,ex
    character(len(s)) :: outs
    character :: aex(len(ex))   ! array of ex chars to extract
    integer   :: i, n

    n = 0  ;  aex = Copy(ex)
    do i = 1,len(s)
        if (.NOT.any(s(i:i) == aex)) cycle   ! dont pack char
        n = n+1 ; outs(n:n) = s(i:i)
    end do
    outs(n+1:) = ' '     ! pad with trailing blanks
    end function Spack

    ! --------------------
    pure integer function Tally (s,text)
    character(*),intent(in) :: s, text
    integer :: i, nt

    Tally = 0 ; nt = len_trim(text)
    do i = 1,len(s)-nt+1
        if (s(i:i+nt-1) == text(:nt)) Tally = Tally+1
    end do
    end function Tally

    ! ---------------------------------
    pure function Translate(s1,codes)  result (s2)
    character(*),intent(in) :: s1, codes(2)
    character(len(s1)) :: s2
    character          :: ch
    integer            :: i, j

    do i = 1,len(s1)
        ch = s1(i:i)
        j = index(codes(1),ch) ; if (j > 0) ch = codes(2)(j:j)
        s2(i:i) = ch
    end do
    end function Translate

    ! ---------------------------------
    pure function Upper(s1)  result (s2)
    character(*),intent(in) :: s1
    character(len(s1)) :: s2
    character          :: ch
    integer,parameter  :: DUC = ICHAR('A') - ICHAR('a')
    integer            :: i

    do i = 1,len(s1)
        ch = s1(i:i)
        if (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
        s2(i:i) = ch
    end do
    end function Upper

    ! ---------------------------------
    pure function Lower(s1)  result (s2)
    character(*),intent(in) :: s1
    character(len(s1)) :: s2
    character          :: ch
    integer,parameter  :: DUC = ICHAR('A') - ICHAR('a')
    integer            :: i

    do i = 1,len(s1)
        ch = s1(i:i)
        if (ch >= 'A'.AND.ch <= 'Z') ch = CHAR(ICHAR(ch)-DUC)
        s2(i:i) = ch
    end do
    end function Lower

end module ModStringFunctions