 SUBROUTINE ensmal(nvarl,lm,locst,glost)
 !*************************************************************************
 !
 !     assembles a local diagonal matrix into a global diagonal matrix
 !     both matrices are stored as arrays
 !
 !*************************************************************************
 USE kinc_db, ONLY: nn,npsdf,nesdf,ftsdf   !(IN) npsdf(:),nesdf(:),ftsdf(:)
 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN) :: nvarl,lm(nvarl)
 REAL (kind=8),INTENT(IN) :: locst(nvarl)
 REAL (kind=8),INTENT(IN OUT) :: glost(*)

 INTEGER (kind=4) :: i,m,ib,ie,neci,necm
 REAL    (kind=8) :: fm,stk

 DO i = 1,nvarl                                   !for each column
   neci = lm(i)                                   !assoc. equation
   stk  = locst(i)                                !value to assemble
   IF(neci > 0) THEN                              !if active dof
     glost(neci) = glost(neci) + stk              !sums on global matrix
   ELSE IF(neci < 0 .AND. neci > -nn) THEN        !if slave dof
     ib = npsdf(-neci)                            !first positon in array
     ie = npsdf(-neci+1)-1                        !last position in array
     DO m=ib,ie
       fm   = ftsdf(m)                            !assoc. factor
       necm = nesdf(m)                            !assoc. equation
       IF(necm > 0) &                             !if active dof
         glost(necm) = glost(necm) + stk*fm*fm    !sums on global matrix
     END DO                                       !m=ib,ie
   END IF
 END DO                                           !i=1,nvarl
 RETURN

 END SUBROUTINE ensmal
