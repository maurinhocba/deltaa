 SUBROUTINE wridva_bpr(d,f,type)
 !
 !  updates and prints auxiliar nodes data
 !
 USE data_db
 IMPLICIT NONE

 REAL(kind=8), POINTER :: d(:,:)
 REAL(kind=8) :: f
 INTEGER(kind=4) :: type ! 0:scalar, 1:vector, 2:tensor

 TYPE (sol3D), POINTER :: eset
 INTEGER :: ielem,j,n,m,nnode,last,dim
 REAL (kind=8) :: val(6)

 eset => sol3d_head  !for first set, point the head

 DO
   IF( eset%etype /= 27 .OR. eset%nnode /= 15 ) THEN
     eset => eset%next
     IF( .NOT.ASSOCIATED(eset)) EXIT
     CYCLE
   END IF

   last = eset%first_l
   nnode = eset%nnode   !15
   dim = SIZE(d,1)

   DO ielem=1,eset%nelem               !process all elements
     DO m=1,3
       val(1:dim) = (d(:,eset%lnods(m,ielem)) + d(:,eset%lnods(m+3,ielem)))/2d0
       last = last+1
       val = val*f
       IF(ip == 2 )THEN
         WRITE(13,2003)last,val(1:dim)
       ELSE IF(ip == 4)THEN
         SELECT CASE (type)
         CASE (0)
           CALL GID_WRITESCALAR(last,val(1))
         CASE (1)
           CALL GID_WRITEVECTOR(last,val(1),val(2),val(3))
         CASE (2)
         END SELECT
       END IF
     END DO
   END DO

   eset => eset%next
   IF( .NOT.ASSOCIATED(eset)) EXIT
 END DO
 RETURN
 2003 FORMAT(i8,3e15.7)

 END SUBROUTINE wridva_bpr
