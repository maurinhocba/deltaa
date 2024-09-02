SUBROUTINE celmnt(ttime,iwric,coord,emass,fcont,coorb,coort)

!.... perform contact searching & computes contact forces

USE cont1_db    !INTENT(IN OUT)
IMPLICIT NONE
!     Dummy arguments
INTEGER (kind=4), INTENT(IN) :: iwric                      !INTENT(IN)
REAL (kind=8), INTENT(IN) :: ttime,emass(:,:)              !INTENT(IN)
REAL (kind=8), POINTER :: coord(:,:),coorb(:,:),coort(:,:) !INTENT(IN)
REAL (kind=8), INTENT(IN OUT) :: fcont(:,:)                !INTENT(OUT)

!     local variables
LOGICAL :: cmptf
INTEGER (kind=4) ipair,ims,iss,i
INTEGER (kind=4), SAVE :: istep = 0
REAL (kind=8), POINTER :: coors(:,:),coorm(:,:)
TYPE (pair1_db), POINTER :: pair
TYPE (surf1_db), POINTER :: master,slave

INTERFACE
  INCLUDE 'csrfc1.h'
END INTERFACE

istep = istep + 1                   !increments counter

pair => headp           !pointer to first pair
DO ipair = 1, npair     !for each pair
  ! if pair is active for this time
  IF( pair%start <= ttime .AND. pair%end >= ttime)THEN
    !.... identify master and slave surfaces numbers
    ims = pair%imast    !master surface order in list
    master => shead     !point to first surface
    DO i=1,ims-1        !loop until pointer is correct
      master => master%next
    END DO
    IF( .NOT.pair%prev )THEN
      pair%prev  = .TRUE.    !set Previous Step as TRUE
      master%cxc = .TRUE.    !set to recompute coordinates
    END IF
    !If necessary set to recompute coordinates of segments centers
    IF( MOD(istep,pair%freq) == 0 )master%cxc = .TRUE.
    iss = pair%islav    !slave surface order in list
    slave => shead      !point to first surface
    DO i=1,iss-1        !loop until pointer is correct
      slave => slave%next
    END DO
    !select coordinates to use for each surface
    SELECT CASE (pair%mtsur)       !for master surface
    CASE(-2)
      coorm => coorb               !bottom
    CASE(-1,1)
      coorm => coord               !central
    CASE( 2)
      coorm => coort               !top
    END SELECT
    SELECT CASE (pair%slsur)       !for slave surface
    CASE(-2)
      coors => coorb               !bottom
    CASE(-1,1)
      coors => coord               !central
    CASE( 2)
      coors => coort               !top
    END SELECT

    cmptf = pair%bhforc > 0d0 .OR. iwric > 0    ! compute total force?
    !.... go into next routine
    CALL csrfc1(ipair,coors,coorm,pair,master,slave,cmptf,emass,fcont)
  END IF
  pair => pair%next     !point to next pair
END DO
RETURN
END SUBROUTINE celmnt
