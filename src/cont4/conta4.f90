SUBROUTINE conta4(itask,ttime,iwrit,veloc,dtime,nsymm,force,stiff,ustif,sname,&
                  maxve)

!     main contac routine (ALGORITHM 4)

USE lispa0
USE npoi_db, ONLY : ndofn, npoin, label, oldlb, emass, coora, coord, coor1
USE kin1_db
USE kin2_db
USE cont_db ! iwric, bottom, top, coorb(:,:), coort(:,:) , fcont(:,:), maxcc
USE cont4_db

IMPLICIT NONE
   !        Dummy arguments
CHARACTER(len=6),INTENT(IN) :: itask  !task to perform
CHARACTER(len=6),INTENT(IN), OPTIONAL :: sname  !surface name
INTEGER (kind=4),INTENT(IN), OPTIONAL :: iwrit,nsymm,maxve
REAL (kind=8),INTENT(IN), OPTIONAL :: dtime,ttime,veloc(:)
REAL (kind=8), INTENT(IN OUT), OPTIONAL :: force(:),stiff(:),ustif(:)

   !        variable for statistical record
INTEGER (kind=4) :: nnn(41,3)
COMMON /proj3/ nnn

   !        Local variables
INTEGER (kind=4) :: ipair,i
TYPE (pair4_db), POINTER :: p

INTERFACE
  INCLUDE 'cinpu4.h'
  !INCLUDE 'ninpu4.h'
  INCLUDE 'celmn4.h'
  INCLUDE 'prsur4.h'
  INCLUDE 'surms4.h'
END INTERFACE

!....  PERFORM SELECT TASK

SELECT CASE (itask)

CASE('ACTUAL')
  !....    updates internal (friction) variables into last converged
  disma = ABS( MAXVAL(veloc(1:neq)) - MINVAL(veloc(1:neq)) )
  IF(nvfix > 0)disma = MAX(disma,(MAX(ABS(MAXVAL(velor(:,nvelr+1))),  &
                                  ABS(MINVAL(velor(:,nvelr+1)))))*dtime)
  disma = 19.d0*disma
  IF(disma == 0d0) disma = 1d0
  CALL conac4(ttime)

CASE ('NEW   ')   !Read Input Data for a NEW problem
  !from npo_db INTENT(IN) :: npoin,coord,label
  !            INTENT(IN OUT) :: bottom,top,emass
  nnn = 0   ! extended check statistics
  ! .. Initializes data base and input element DATA
  CALL cinpu4(maxve,iwrit,label(:),npoin,coord)

  ! compute mass and detect surface to use
  CALL surms4(emass,coord)

!CASE ('NSTRAT') !Modifies Data
  !!from npo_db INTENT(IN) :: npoin,coora,label
  !!            INTENT(IN OUT) :: bottom,top,emass
  !CALL ninpu4(maxve,iwrit,label(:),npoin,coora,   &
  !            oldlb(:),itask,esets,eset)
  !! compute mass and detect surface to use
  !CALL surms4(emass,coord)

CASE ('FORCES')  !Performs contact search and computes contact forces
  !from npo_db INTENT(IN) :: npoin,coora,coorb,coort,label,emass
  !            INTENT(IN OUT) :: fcont

  ! ....  Perform Contact Searching & Compute contact forces
  CALL celmn4(ttime,iwric,coor1,emass,fcont,coorb,coort)

  !DO i=1,npoin
  !  IF(fcont(3,i) /=0 )WRITE(55,"(i10,e15.5)")label(i),fcont(3,i)
  !END DO

CASE ('STIFFM')
  !....    compute contact stiffness contribution
  CALL cstif4(ttime,force(1),stiff(1))

CASE ('DUMPIN')   !write data to a re-start file
  CALL cdump4 (npoin)
  WRITE(55,"(i5,3i20)")(i,nnn(i,1:3),i=1,41)  !statistic to debug file

CASE ('RESTAR')   !read data from a re-start file
  !from npo_db INTENT(IN) :: bottom,top
  CALL crest4 (npoin)
  ALLOCATE ( surtf(3,npair) )              !get memory for total forces

CASE ('OUTPUT')   !Writes contact forces between pairs
  IF( iwric > 0 )THEN
    WRITE(41)ttime                        !control variable
    p => headp
    DO ipair=1,npair                      !for each pair
      WRITE(41)surtf(1:3,ipair)           !average contact forces
      IF(p%bhforc > 0d0 .AND. p%start < ttime .AND. p%end  > ttime) &
        WRITE(lures,"(' Normal penalty for pair ',a6,' is',E12.5)") &
                                                  p%pname, p%npenal
      p => p%next
    END DO
  END IF

CASE ('WRTPOS')
  ! OPEN files for postprocess
  CALL prsur4(iwric,bottom,top)

END SELECT

RETURN
END SUBROUTINE conta4
