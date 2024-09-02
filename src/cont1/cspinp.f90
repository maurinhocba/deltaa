SUBROUTINE cspinp(p)

!.... read contact pair data

USE param_db, ONLY : mich
USE c_input
USE cont1_db
IMPLICIT NONE
!     arguments
TYPE (pair1_db), POINTER :: p  !INTENT(OUT)

!     local variables
INTEGER (kind=4) j
CHARACTER (len=6) :: surpos
CHARACTER(len=mich):: inttoch


!.... read master and slave contact surfaces for each pair
IF (exists('IPAIR ',j)) THEN     !check Key-word
  IF( param(j) == 0)THEN         !if no associated number
    p%pname = words(j+1)(1:6)    !label is the next word
  ELSE                           !label is a two digit number
    p%pname = '    '//TRIM(inttoch(INT(param(j)),2)) !fill with blanks
  END IF
  WRITE(lures,"(/'  PAIR NAME:',a8,/)")p%pname
ELSE                             !If kew-word not found
  CALL runend ('CSPINP: Pair Name not given !      ')
END IF

IF (exists('SLAVE ',j)) THEN            !check Key-word
  IF( param(j) == 0)THEN                !if no associated number
    p%slave = words(j+1)(1:6)           !label is the next word
    surpos  = words(j+2)(1:6)           !position (top center bottom)
  ELSE                                  !label is a two digit number
    p%slave = '    '//TRIM(inttoch(INT(param(j)),2)) !fill with blanks
    surpos  = words(j+1)(1:6)           !position (top center bottom)
  END IF
  SELECT CASE (surpos)
  CASE ('BOTTOM')
    p%slsur = -2
  CASE ('REVERS')
    p%mtsur = -1
  CASE ('TOP   ')
    p%slsur =  2
  CASE DEFAULT
    p%slsur =  1
    surpos = 'CENTER'
  END SELECT
  WRITE(lures,"('         Slave surface :',a8, &
               &'         Surface posit.:',a8)")p%slave,surpos
ELSE
  CALL runend ('CSPINP: Slave surface not defined !')
END IF
IF (exists('MASTER',j)) THEN           !check Key-word
  IF( param(j) == 0)THEN               !if no associated number
    p%master = words(j+1)(1:6)         !label is the next word
    surpos  = words(j+2)(1:6)          !position (top center bottom)
  ELSE                                 !label is a two digit number
    p%master = '    '//TRIM(inttoch(INT(param(j)),2)) !fill with blanks
    surpos  = words(j+1)(1:6)          !position (top center bottom)
  END IF
  SELECT CASE (surpos)
  CASE ('BOTTOM')
    p%mtsur = -2
  CASE ('REVERS')
    p%mtsur = -1
  CASE ('TOP   ')
    p%mtsur =  2
  CASE DEFAULT
    p%mtsur =  1
    surpos = 'CENTER'
  END SELECT
  WRITE(lures,"('         Master surface:',a8, &
               &'         Surface posit.:',a8)")p%master,surpos
  IF( p%master == p%slave ) &
   CALL runend ('CSPINP: Slave == Master not allowed')
ELSE
  CALL runend ('CSPINP: Master surface not defined!')
END IF

!.... read contact pair properties data

p%indcon =getint('INDCON',0,'!Contact forces flag ..............')
SELECT CASE (p%indcon)
CASE (0)
  WRITE(lures,"(15x,'Contact Forces applied on BOTH surfaces')")
CASE (1)
  WRITE(lures,"(15x,'Contact Forces applied on SLAVE surface ONLY')")
CASE (2)
  WRITE(lures,"(15x,'Contact Forces applied on MASTER surface ONLY')")
CASE DEFAULT
  p%indcon = 0
  WRITE(lures,"(15x,'Contact Forces applied on BOTH surfaces')")
END SELECT

p%npenal =getrea('NPENAL',0d0,'!Normal penalty coeff .............')
p%tpenal =getrea('TPENAL',0d0,'!Tangential penalty coeff..........')
p%static =getrea('STATIC',0d0,'!Friction coefficient .............')
p%cutoff=-getrea('CUTOFF',1d9,' CUT OFF value for penetration.....')
IF( p%cutoff == 0d0 ) p%cutoff = -1d9
p%bhforc =getrea('BHFORC',0d0,' Blank Holder Force ...............')
p%start  =getrea('START ',0d0,' Activation time ..................')
p%end    =getrea('END   ',1d9,' Deactivation time  ...............')
p%freq  = getint('FREQ  ',100,' Frequency of global search .......')

p%prev = .FALSE.  !not existence in previous step
WRITE(lures,"(//)")
RETURN
END SUBROUTINE cspinp
