SUBROUTINE cinpu4(maxcc,iwrit,label,npoin,coord)

!.... contact input routine (NEW problem)

USE c_input
USE cont4_db  !INTENT(OUT)

IMPLICIT NONE
!     arguments
INTEGER (kind=4),INTENT (IN) :: iwrit,npoin,label(:),maxcc
REAL (kind=8), INTENT (IN) :: coord(:,:)

!     local variables
TYPE (pair4_db), POINTER :: pair
INTEGER (kind=4) tsurf,oldsr,oldpr,surdi
INTEGER (kind=4), ALLOCATABLE :: codes(:,:)
CHARACTER (len=6), ALLOCATABLE :: surname(:)

INTERFACE
  INCLUDE 'cspin4.h'
  INCLUDE 'csinp4.h'
END INTERFACE

!...  Initializes Data Base
npair = 0            !number of pairs
nsurf = 0            !number of surface
oldis = 0d0          !latest maximum displacement increment
CALL ini_cont4(headp,tailp)   !initializes pair list
CALL ini_srf4 (shead,stail)   !initializes surface list

!.... Title print and first line

IF(iwrit /= 0)WRITE(lures, &
  "( //10x,'*** D A T A   F O R   C O N T A C T ***',//, &
      & 10x,'*** A L G O R I T H M   4  ***',//)")

CALL listen('CINPU4')  !Read first line

ctime = getrea('DTIME ',0d0,' Dtime to compute penalty param. ..')
IF( ctime == 0d0) ctime = -1d0
wear  = exists('WEARIN')
IF( wear ) THEN
  WRITE(lures,"(15x,'Work due to Friction will be computed')")
  ALLOCATE( wwear(npoin) )    !reserve space
  wwear = 0d0                 !initialize
END IF

!.... input contact surfaces pairs DATA

CALL listen('CINPU4')  !read second line
IF (.NOT.exists('PAIRIN'))   &  ! pair information expected
    CALL runend('CSPINP:PAIR_INFO CARD EXPECTED     ')

WRITE(lures, &
  "(/,'  C O N T A C T   S U R F A C E S   P A I R S   D A T A '/)")
npair = 0   !initializes
DO
  CALL listen ('CINPU4')        !read a line
  IF (exists('ENDPAI')) EXIT    !end of pair data detected, EXIT
  CALL new_pair4 (pair)         !Initializes new pair
  CALL cspin4(pair)             !Read pair data
  npair = npair + 1             !increase pair counter
  CALL add_pair4 (pair, headp, tailp) !add to end of the list
END DO
                         !set maximum number of surfaces
surdi = 2*npair          !twice npair seems to be enough
ALLOCATE ( codes(2,surdi), surname(surdi) )  !allocate auxiliar arrays
CALL chksu4(surdi,codes,surname,tsurf)       !generate auxiliar arrays
!tsurf = number of expected surfaces
!surname = surface names
!codes = (1)times as master (2) times as slave

!.... input and generate contact surfaces DATA

oldsr = 0  !number of old surfaces
oldpr = 0  !number of old pairs
CALL csinp4(maxcc,iwrit,label,npoin,coord,oldsr,oldpr, &
            codes,surname,tsurf,surdi)

DEALLOCATE( codes , surname )  !release auxiliar arrays
ALLOCATE (surtf(3,npair) )     !reserve space for total contact forces

WRITE(*,"('+ CONTAC: End of initial data        ')")

RETURN
END SUBROUTINE cinpu4
