SUBROUTINE prtcob(nnode,dim,sname,nodes,iset)
! print coordinates for GiD
USE data_db,ONLY: npoin, ndime, wsdisp, wtdisp, label, coord, coors, meshn, tdisp_f, sdisp_f
IMPLICIT NONE

  INTEGER(kind=4),INTENT(IN):: nnode,dim,nodes(npoin),iset
  CHARACTER(len=32),INTENT(IN):: sname
  REAL(kind=8) :: x(3)

  INTEGER (kind=4) i, elmt

  SELECT CASE (dim)           !element dimension
  CASE (0)
    IF (nnode == 1) elmt=1
  CASE (1) !linear elements
    IF (nnode == 2) elmt=2
    IF (nnode == 3) elmt=2
  CASE (2)
    IF (nnode == 2) elmt=2
    IF (nnode == 3) elmt=3
    IF (nnode == 4 .OR. nnode == 8) elmt=4
  CASE (3)
    IF (nnode == 4) elmt=5
    IF (nnode == 8 .OR. nnode == 20) elmt=6
    IF (nnode == 6 .OR. nnode == 15) elmt=7
  END SELECT


  CALL GID_BEGINMESH(TRIM(sname),ndime,elmt,nnode)
  CALL GID_BEGINCOORDINATES()
  IF (ndime==2) x(3) = 0d0
  IF( .NOT.wsdisp .AND. wtdisp )THEN  !write original coordinates
    DO i=1,npoin
      IF (nodes(i) == iset) THEN
        x(1:ndime) = coord(:,i)*tdisp_f
        CALL GID_WRITECOORDINATES(label(i),x(1),x(2),x(3))
      END IF
      IF (nodes(i) > 0) meshn(i)=.TRUE.
    END DO
  ELSE !Write stage coordinates
    DO i=1,npoin
      IF (nodes(i) == iset) THEN
        x(1:ndime) = coors(:,i)*sdisp_f
        CALL GID_WRITECOORDINATES(label(i),x(1),x(2),x(3))
      END IF
      IF (nodes(i) > 0) meshn(i)=.TRUE.
    END DO
  END IF
  CALL GID_ENDCOORDINATES()

RETURN
END SUBROUTINE prtcob
