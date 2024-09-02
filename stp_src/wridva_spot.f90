 SUBROUTINE wridva_spot(d,f)
 !
 !  updates and prints auxiliar nodes data
 !
 USE data_db
 IMPLICIT NONE

 REAL(kind=8), POINTER :: d(:,:)
 REAL(kind=8) :: f

 TYPE (spot), POINTER :: eset
 INTEGER :: ielem,j,n,m,nnode,last
 REAL (kind=8) :: val(3)

 eset => spot_head  !for first set, point the head
 last = label(npoin)

 DO
   nnode = eset%nnode

   IF( nnode == 6 )THEN
     DO ielem=1,eset%nelem               !process all elements
       j = 0
       DO m=1,2
         val = 0d0
         DO n=1,3
           j = j + 1
           val = val + d(:,eset%lnods(j,ielem))*eset%lc(n,m,ielem)
         END DO
         last = last+1
         val = val*f
         IF(ip == 2 )THEN
           WRITE(13,2003)last,val
         ELSE IF(ip == 4)THEN
           CALL GID_WRITEVECTOR(last,val(1),val(2),val(3))
         END IF
       END DO
     END DO
   ELSE
     val = 0d0
     DO ielem=1,eset%nelem
       IF(eset%lnods(2,ielem) > 0 )CYCLE
       last = last+1
       val = val*f
       IF(ip == 2 )THEN
         WRITE(13,2003)last,val
       ELSE IF(ip == 4)THEN
         CALL GID_WRITEVECTOR(last,val(1),val(2),val(3))
       END IF
     END DO
   END IF
   eset => eset%next
   IF( .NOT.ASSOCIATED(eset)) EXIT
 END DO
 RETURN
 2003 FORMAT(i8,3e13.5)

 END SUBROUTINE wridva_spot
