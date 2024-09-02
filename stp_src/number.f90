INTEGER FUNCTION number(string)
!
! convert a string into an integer
! first character in string can not be a blank space
!
!
IMPLICIT NONE
! dummy argument
CHARACTER(len=*), INTENT(IN) :: string
! local variables
INTEGER (kind=4) nchar,i,n
CHARACTER ch

nchar = LEN_TRIM(string)   !string length
number=0    !initializes answer
DO i=1,nchar                  !for each character
  ch = string(i:i)            !isolate character
  IF(ch /= ' ') THEN          !if not blank space
    n = ICHAR(ch)-48          !numerical value of character
    IF(n >= 0 .AND. n <= 9) THEN   !check it is a digit
      number = number*10 + n       !add to answer
    ELSE                      !error if a non-digit present
      WRITE(6,*)'ERROR:String cannot be converted into an integer'
      STOP
    END IF
  ELSE IF(number == 0) THEN   !ignore trailing blanks
    CYCLE
  ELSE                        !end of integer found
    EXIT
  END IF
END DO

RETURN
END FUNCTION number
