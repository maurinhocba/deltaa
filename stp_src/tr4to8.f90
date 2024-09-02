 SUBROUTINE tr4to8(gausv,newgv,iv,nv,ng,mv)
 !          For 8-node elements
 ! transfer gauss point information from 4-Gauss points to 8-Gauss points(GiD)
 !
 ! Original gauss points (a = 1/sqrt(3)
 !           In progs                                     Inverse
 !       4-G           8-G      GiD                   GiD   4G    8G
 ! G1 = (/+a,-a,-a /)    2       2                  !  1   N1     1
 ! G2 = (/-a,+a,-a /)    3       4 x                !  2   G1     2
 ! G3 = (/-a,-a,+a /)    5       5                  !  3   N2     4
 ! G4 = (/+a,+a,+a /)    8       7 x                !  4   G2     3
 ! N1 = (/-a,-a,-a /)    1       1                  !  5   G3     5
 ! N2 = (/+a,+a,-a /)    4       3 x                !  6   N3     6
 ! N3 = (/+a,-a,+a /)    6       6                  !  7   G4     8
 ! N4 = (/-a,+a,+a /)    7       8 x                !  8   N4     7
 !
 ! ALSO     For 8 and 6-node elements
 ! transfer gauss point information from external faces to 6/8-Gauss points (GID)
 !
 IMPLICIT NONE
 INTEGER(kind=4), INTENT(IN) :: iv, & !last position to transfer in array GAUSV
                                nv, & !number of variables to transfer
                                ng       !number of original gauss points
 REAL(kind=8) :: gausv(:,:)           !original variables
 REAL(kind=8) :: newgv(6,8)           !transfered variables
 REAL(kind=8), OPTIONAL :: mv         !minimum value

 INTEGER(kind=4) :: i,g,m,n
 REAL(kind=8) :: avera(nv)
 INTEGER(kind=4), PARAMETER :: j(4) = (/ 2,4,5,7 /), & !direct positions
                               k(4) = (/ 1,3,6,8 /), & !extrapolated points
                               l(4) = (/ 4,3,2,1 /)    !auxiliar
 REAL (kind=8), PARAMETER :: f(2,2) = (/0.788675134594813, 0.211324865405187, 0.211324865405187, 0.788675134594813/)

 i = iv-nv  !pointer to first value in array GAUSV
 IF( ng > 0 )THEN
   !loop over original values to compute average values (twice)
   DO g=1,nv
     avera(g) = SUM(gausv(i+g,1:4))/2d0
   END DO
   ! loop to transfer
   DO g=1,4
     newgv(1:nv,k(g)) = avera - gausv(i+1:iv,l(g))  !new interpolated values
     newgv(1:nv,j(g)) = gausv(i+1:iv,g)             !transfered values
   END DO
   ! if a minimum value exists, check values (used for positive values)
   IF( PRESENT (mv ) ) THEN
     DO g=1,4        !for each interpolated point
       DO i=1,nv       !for each variable
         IF( newgv(i,k(g)) < mv ) newgv(i,k(g)) = mv
       END DO
     END DO
   END IF
                          !interpolate from external surfaces
 ELSE !IF (ng <= 0) THEN  !interpolate to Gauss points
   n = 4                  !Hexahedra
   IF( ng ==-1 ) n = 3    !Prism
   DO g=1,2
     DO m=1,n
       newgv(1:nv,m+n*(g-1)) = f(1,g)*gausv(i+1:iv,1)+f(2,g)*gausv(i+1:iv,2)
     END DO
   END DO
 END IF

 RETURN
 END SUBROUTINE tr4to8
