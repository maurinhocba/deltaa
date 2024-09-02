 SUBROUTINE ensvec(nvarl,lm,locvc,glovc)
 !*************************************************************************
 !
 !     assemble a local vector into a global vector
 !
 !*************************************************************************
 USE kinc_db, ONLY : nn,npsdf,nesdf,ftsdf,ftsd0
 USE Ctrl_db, ONLY : inverse
 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN) :: nvarl,lm(nvarl)
 REAL (kind=8),INTENT(IN) :: locvc(nvarl)
 REAL (kind=8),INTENT(IN OUT) :: glovc(*)

 INTEGER (kind=4) :: i,j,k,nec,ib,ie

 DO i = 1,nvarl                                           !for each value
   nec = lm(i)                                            !assoc. equation
   SELECT CASE (nec)
   CASE (1:)                                              !if active dof
     glovc(nec) = glovc(nec) + locvc(i)                   !sums on global
   CASE (-nn:-1)                                          !if slave dof
     ib = npsdf(-nec)                                     !first post.
     ie = npsdf(-nec+1)-1                                 !last post.
     DO j=ib,ie                                           !for each master
       k = nesdf(j)                                       !assoc. equation
       IF(k > 0)THEN
         IF(inverse) THEN
           glovc(k) = glovc(k) + locvc(i)*ftsd0(j)  !sums on global
         ELSE
           glovc(k) = glovc(k) + locvc(i)*ftsdf(j)  !sums on global
         END IF
       END IF
     END DO
   END SELECT
 END DO
 RETURN

 END SUBROUTINE ensvec
