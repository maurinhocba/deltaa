 SUBROUTINE readln (c,lin,word,value,flags,names,cname,units)
 !***********************************************************************
 !
 !***reads a string and interprets it as key-words a parameter, flags
 !                variable name and component names
 !
 !   - maximum number of flags = 2
 !   - maximum number of component names  = 6
 !   - only the first five characters of each word are decoded.
 !   - the underline characters '_' are discarted.
 !   - lower case letters are converted to upper case.
 !   - variable name must be sorrounded by { }
 !   - component names must be sorrounded by [ ]
 !   - each component name must be separated by ','
 !
 !***********************************************************************

 IMPLICIT NONE

 !***  arguments
 CHARACTER(len=1 ) ,INTENT(IN OUT) :: c
 CHARACTER(len=99) ,INTENT(IN OUT) :: lin
 CHARACTER(len=5 ) ,INTENT(OUT) :: word
 CHARACTER(len=1 ) ,INTENT(OUT) :: flags(2)
 CHARACTER(len=15) ,INTENT(OUT) :: names(2)
 CHARACTER(len=15) ,INTENT(OUT) :: cname(6,2)
 CHARACTER(len=15) ,INTENT(OUT) :: units
 REAL(kind=8) ,INTENT(OUT) :: value


 REAL(kind=8) ::  digit
 INTEGER (kind=4) :: k,first,last,leng,ioflag,flag,ib,ie,ic
 CHARACTER (len=100) :: card
 CHARACTER :: letra
 LOGICAL :: nread,cread,cont

 INTERFACE
   SUBROUTINE upcase(string,n)
   IMPLICIT NONE
   INTEGER (kind=4), INTENT(IN) :: n
   CHARACTER(len=*), INTENT(IN OUT) :: string
   END SUBROUTINE upcase
 END INTERFACE

 !***  begin

 card  = c//lin    ! concatenate First char and rest of the line
 word  = '     '   ! initializes WORDS
 names = '               '   ! initializes NAMES
 cname = '               '   ! initializes Component NAMEs
 flags = ' '
 value = 0.0       ! initializes value
 units = ''        ! initializes units

 flag  = 0          ! inializes number of flags
 nread = .FALSE.    ! inializes to name read
 cread = .FALSE.    ! inializes to component read

 leng = LEN_TRIM(card)                ! compute the length.
 cont = ( card(leng:leng) == '/' ) .OR. ( card(leng:leng) == '\' )

 ! see if last character is a continuation one

 IF(leng == 0 )RETURN               ! blank line

 ! first find the keyword

 !  block to find the FIRST character of the Key-word
 first = 1                   ! first non-processed character
 DO                          ! Skip over auxiliar characters
   IF(first > leng) RETURN            ! No valid char RETURN
   letra = card(first:first)          ! present character
   IF( letra /= ' ' ) EXIT
   first = first+1                    ! new first non-blank char
 END DO
 !         read key word
 !  loop to find the LAST character of the Key-word
 last = first                         !initializes
 ! find first auxiliar characters
 DO
   letra = card(last:last)            ! present character
   IF( letra == ' ' .OR.  &           ! Look for separator ( =:,)
       letra == '=' .OR.  &
       letra == ':' .OR.  &
       letra == ',' ) EXIT
   last = last+1                      ! point to next character
   IF( last > leng )EXIT              ! look for END of line
 END DO

 last = last - 1                      ! correct last position
 k = 0
 DO ib=first,last
   letra = card(ib:ib)
   IF( letra == '_' )CYCLE
   k = k + 1
   word(k:k) = letra        ! keyword
   IF( k == 5 )EXIT         ! no more than 5 chars
 END DO
 CALL upcase(word,k)       ! convert to upper case.

 outer : DO                            ! decode the rest of the card.
   first = last + 1
   IF( first > leng )THEN
     IF( cont ) THEN
       READ(7,'(a1,a99)',IOSTAT=ioflag)c,lin   !read a line
       card  = c//lin    ! concatenate First char and rest of the line
       leng = LEN_TRIM(card)                ! compute the length.
       cont = ( card(leng:leng) == '/' )
       first = 1
     ELSE
       EXIT
     END IF
   END IF
   DO                          ! Skip over auxiliar characters
     IF(first > leng) RETURN            ! No valid char RETURN
     letra = card(first:first)          ! present character
     IF( letra /= ' ' ) EXIT
     first = first+1                    ! new first non-blank char
   END DO
   SELECT CASE (letra)
   CASE ('A':'Z','a':'z')               ! Read a Flag
     flag = flag+1
     IF( flag > 2 )RETURN !too many flags
     IF(ICHAR(letra) >= 97) letra = CHAR(ICHAR(letra)-32)      ! convert to upper case.
     !CALL upcase(letra,1)                     ! convert to upper case.
     flags(flag) = letra
     nread = .TRUE.
     last = first+1
     DO
       letra = card(last:last)            ! present character
       IF( letra == ' ' )EXIT             ! Look for separator ( )
       last = last+1                      ! point to next character
       IF( last > leng )EXIT              ! look for END of line
     END DO

   CASE ('{')                           ! Read Variable Name
     last = first + 1
     DO
       letra = card(last:last)            ! present character
       IF( letra == '}' )EXIT             ! Look for separator }
       last = last+1                      ! point to next character
       IF( last > leng )THEN              ! look for END of line
         !error reading a name
       END IF
     END DO
     IF(nread)THEN
       names(flag) = card(first+1:last-1)
       nread = .FALSE.
       cread = .TRUE.
     END IF
   CASE ('[')                           ! Read component Names
     last = first + 1
     DO
       letra = card(last:last)            ! present character
       IF( letra == ']' )EXIT             ! Look for separator }
       last = last+1                      ! point to next character
       IF( last > leng )THEN              ! look for END of line
         !error reading components
         PRINT *, card
         STOP 'lines in STP.CFG files must not exceed 90 chars'
       END IF
     END DO
     IF(cread)THEN
       ib = first+1
       ic = 0
       DO
         ie = ib
         DO
           letra = card(ie:ie)
           IF( letra == ']' .OR. letra == ',' )EXIT
           ie = ie + 1
         END DO
         ic = ic + 1
         IF( ic > 6 )EXIT
         cname(ic,flag) = card(ib:ie-1)
         ib = ie+1
         IF( ib >= last )EXIT
       END DO
       cread = .FALSE.
     ELSE
       units = card(first+1:last-1)
     END IF
   CASE ('0':'9','.','-')
     last = first+1
     DO
       letra = card(last:last)            ! present character
       IF( letra == ' ' )EXIT             ! Look for separator ( )
       last = last+1                      ! point to next character
       IF( last > leng )EXIT              ! look for END of line
     END DO

     READ(card(first:last),'(G20.0)',IOSTAT = ioflag) digit  !read a number
     IF(ioflag == 0) THEN              ! O.K. it is a parameter.
       value = digit                   ! store param
     ELSE                              ! it is a word.
       !error reading a number
     END IF
   CASE DEFAULT
     last = first
   END SELECT
 END DO outer

 RETURN

 END SUBROUTINE readln

 SUBROUTINE upcase(string,n)
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: n
 CHARACTER, INTENT(IN OUT) :: string(n)

 INTEGER (kind=4) i,k

 DO i=1,n
   k = ICHAR(string(i))
   IF(k >= 97 .AND. k <= 122)string(i)=char(k-32)
 END DO

 RETURN
 END SUBROUTINE upcase
