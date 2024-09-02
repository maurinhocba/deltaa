SUBROUTINE conta1(itask,ttime,iwrit,veloc,dtime,nsymm,force,stiff,ustif, &
               sname,maxve)
!     main contac routine (ALGORITHM 1) 2D problems

USE lispa0
USE npoi_db, ONLY : ndofn, npoin, label, oldlb, emass, coora, coord, &
                    coor1
USE kin1_db
USE kin2_db
USE cont_db   ! bottom, top, fcont(:,:), coorb(:,:), coort(:,:)
USE cont1_db

IMPLICIT NONE
   !        Dummy arguments
CHARACTER(len=6),INTENT(IN) :: itask  !task to perform
CHARACTER(len=6),INTENT(IN), OPTIONAL :: sname ! surface name
INTEGER (kind=4),INTENT(IN), OPTIONAL :: iwrit,nsymm,maxve
REAL (kind=8),INTENT(IN), OPTIONAL :: dtime,ttime,veloc(:)
REAL (kind=8), INTENT(IN OUT), OPTIONAL :: force(:),stiff(:),ustif(:)

   !        Local variables
INTEGER (kind=4) :: ipair
TYPE (pair1_db), POINTER :: p

INTERFACE
  INCLUDE 'cinput.h'
  !INCLUDE 'ninput.h'
  INCLUDE 'celmnt.h'
  INCLUDE 'prsurf.h'
  INCLUDE 'surms1.h'
END INTERFACE

!....  PERFORM SELECT TASK

SELECT CASE (itask)

CASE('ACTUAL')
  !....    updates internal (friction) variables into last converged
  disma = ABS( MAXVAL(veloc(1:neq)) - MINVAL(veloc(1:neq)) )
  IF(nvfix > 0)disma = MAX(disma,(MAX(ABS(MAXVAL(velor(:,nvelr+1))),  &
                                  ABS(MINVAL(velor(:,nvelr+1)))))*dtime)
  disma = 16.d0*disma
  IF(disma == 0d0) disma = 1d0
  CALL conact(ttime)

CASE ('NEW   ')   !Read Input Data for a NEW problem
  !from npo_db INTENT(IN) :: npoin,coord,label
  !            INTENT(IN OUT) :: bottom,top,emass
  ! .. Initializes data base and input element DATA
  CALL cinput(maxve,iwrit,label(:),npoin,coord)

  ! compute mass and detect surface to use
  CALL surms1(emass,coord)

!CASE ('NSTRAT') !Modifies Data
  !!from npo_db INTENT(IN) :: npoin,coora,label
  !!            INTENT(IN OUT) :: bottom,top,emass
  !CALL ninput(maxve,iwrit,label(:,2),npoin,coora,   &
  !            label(:,1),itask,esets,eset)
  !! compute mass and detect surface to use
  !CALL surms1(emass,coord)

CASE ('FORCES')  !Performs contact search and computes contact forces
  !from npo_db INTENT(IN) :: npoin,coora,coorb,coort,label,emass
  !            INTENT(IN OUT) :: fcont

  ! ....  Perform Contact Searching & Compute contact forces
  CALL celmnt(ttime,iwric,coor1,emass,fcont,coorb,coort)

CASE ('STIFFM')
  !....    compute contact stiffness contribution
  CALL cstiff(nsymm,ttime,force(1),stiff(1),ustif(1))

CASE ('DUMPIN')   !write data to a re-start file
  CALL cdump1 (npoin)

CASE ('RESTAR')   !read data from a re-start file
  !from npo_db INTENT(IN) :: bottom,top
  CALL crest1 (npoin)
  ALLOCATE ( surtf(2,npair) )              !get memory for total forces

!CASE ('AREMES')   !Modifies data for MODIFIED surface SNAME
!CALL chsur1(sname)

CASE ('OUTPUT')   !Writes contact forces between pairs
  IF( iwric > 0 )THEN
    WRITE(41)ttime                        !control variable
    p => headp
    DO ipair=1,npair                      !for each pair
      WRITE(41)surtf(1:2,ipair)           !average contact forces
      IF(p%bhforc > 0d0 .AND. p%start < ttime .AND. p%end  > ttime) &
        WRITE(lures,"(' Normal penalty for pair ',a6,' is',E12.5)") &
                                                  p%pname, p%npenal
      p => p%next
    END DO
  END IF
CASE ('WRTPOS')    ! OPEN files for postprocess
  CALL prsurf(iwric,bottom,top)

END SELECT

RETURN
END SUBROUTINE conta1
