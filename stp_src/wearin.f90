SUBROUTINE wearin(npoin,ttime,label,ip)
! reads contact friction work
USE data_db, ONLY : loadstep
USE cont_db
IMPLICIT NONE
INTEGER(kind=4), INTENT(IN) :: npoin,label(npoin),ip
REAL(kind=8), INTENT(IN) :: ttime

INTEGER(kind=4) :: i
REAL(kind=8) :: vmax(1)

 INTERFACE
   INCLUDE 'cab_gid.h'
   INCLUDE 'cab_gid_bin.h'
 END INTERFACE

  READ(45) (wwork(i),i=1,npoin) !friction work for each point in the mesh
  ! uses associated area to consider
  DO i=1,npoin
    IF(areas(i) /= 0d0) THEN
      wwork(i) = wwork(i)/ABS(areas(i))*wear_f         !includes wear factor to change units
    ELSE
      wwork(i) = 0d0
    END IF
  END DO
  IF(LEN(wear_l(2)) == 0 )THEN
    vmax = MAXVAL(wwork)   !compute maximum value to normalize
    IF( vmax(1) > 0d0 )wwork = wwork/vmax(1)
  END IF
  IF (ip == 2) THEN
    CALL cab_gid(1,1,wear_l(2),wear_l(2),1,loadstep,ttime,units=wear_u)
    DO i=1,npoin
      IF(areas(i) < 0d0) WRITE(13,"(i8,e13.5)")label(i),wwork(i)
    END DO
    WRITE(13,"('End Values')")
    CALL cab_gid(1,1,wear_l(3),wear_l(3),1,loadstep,ttime,units=wear_u)
    DO i=1,npoin
      IF(areas(i) > 0d0) WRITE(13,"(i8,e13.5)")label(i),wwork(i)
    END DO
    WRITE(13,"('End Values')")

  ELSE IF (ip == 4 .OR. ip == 5) THEN
    CALL cab_gid_bin(1,1,wear_l(2),wear_l(2),1,loadstep,ttime,units=wear_u)
    DO i=1,npoin          !for each point
      IF(areas(i) < 0d0) CALL GID_WRITESCALAR(label(i),wwork(i))
    END DO
    CALL GID_ENDRESULT()
    CALL cab_gid_bin(1,1,wear_l(3),wear_l(3),1,loadstep,ttime,units=wear_u)
    DO i=1,npoin          !for each point
      IF(areas(i) > 0d0) CALL GID_WRITESCALAR(label(i),wwork(i))
    END DO
    CALL GID_ENDRESULT()

  END IF

RETURN
END SUBROUTINE wearin
