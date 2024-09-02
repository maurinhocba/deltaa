 SUBROUTINE smom05(nnode,ngaus,emass,dvolu,shape)
 !**********************************************************************
 !
 !****this routine evaluates the element smoothing matrix and inverts it
 !    for element 05
 !**********************************************************************
 IMPLICIT NONE

 !     routine variables

 INTEGER (kind=4) nnode,ngaus
 REAL   (kind=8) emass(nnode,nnode),dvolu(ngaus),shape(nnode,ngaus)

 !     local variables

 INTEGER (kind=4) i,j,g,err
 REAL    (kind=8) a(nnode,nnode),auxva,denom

 INTERFACE
   INCLUDE 'findinv.h'
 END INTERFACE
 !***consistent smoothing matrix


 a = 0d0

 DO g=1,ngaus
   DO i=1,nnode
     auxva = shape(i,g)*dvolu(g)
     DO j=i,nnode
       a(i,j) = a(i,j) + auxva*shape(j,g)
     END DO
   END DO
 END DO

 DO i=1,nnode
   DO j=1,i-1
     a(i,j) = a(j,i)
   END DO
 END DO

 ! inverse of an n x n matrix
 IF( nnode == 6 )THEN
   denom = a(1,1)*a(4,4) - a(1,4)**2
   auxva = a(4,4)/denom
   emass(1:3,1:3) = auxva
   auxva = a(1,1)/denom
   emass(4:6,4:6) = auxva
   auxva =-a(1,4)/denom
   emass(1:3,4:6) = auxva
   emass(4:6,1:3) = auxva
 ELSE
   IF( ngaus == 4 )THEN
   ! 4 gauss gives as poor consistent mass matrix
     emass = 0d0
     DO i=1,nnode
       denom = SUM( a(1:nnode,i) )
       emass(i,i) = 1d0/denom
     END DO
   ELSE
     CALL findinv(a,emass,nnode,err)
     ! de momento inversa de la matriz diagonalizada
     IF( err == -1 )THEN
       emass = 0d0
       DO i=1,nnode
         denom = SUM( a(1:nnode,i) )
         emass(i,i) = 1d0/denom
       END DO
     END IF
   END IF
 END IF
 RETURN
 END SUBROUTINE smom05
