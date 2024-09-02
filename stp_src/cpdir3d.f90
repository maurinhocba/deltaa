 SUBROUTINE  cpdir3d(ndime,nnode,nelem,lnods,dirt,coord)

 ! computes thicknes direction for each element

 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: ndime, & !=2 problem dimension
                                 nnode, & !3 or 4 number of nodes per element
                                 nelem, & !number of elements in the mesh
                                 lnods(:,:)     !connectivities
 REAL(kind=8), INTENT(IN) :: coord(:,:)        !mesh coordinates
 REAL(kind=8), INTENT(OUT) :: dirt(:,:)        !thickness direction
                                                    !one for each element
 ! local variables
 INTEGER (kind=4) :: ne,i,j,k,is,ms(1),n1,n2,k1,k2,m1,m2
 INTEGER (kind=4), POINTER :: bound(:,:)
 INTEGER (kind=4), ALLOCATABLE :: side(:),ns(:)
 REAL (kind=8) :: t(ndime),n(ndime),l,l1,l2
 REAL (kind=8), ALLOCATABLE :: norm(:,:),top(:,:,:),bot(:,:,:),angle(:)
 REAL (kind=8), PARAMETER :: pi4 = 0.7D0
 LOGICAL :: first,chang

! INTERFACE
!   SUBROUTINE surf05 ( lnods, nnode, nelem, bound, ne)
!   !******************************************************************
!   !
!   !*** get the boundary definition from the element set
!   !    segments are oriented (outward normal)
!   !
!   !******************************************************************
!   IMPLICIT NONE
!   !dummy arguments
!   INTEGER (kind=4), INTENT(IN) :: nelem      !number of elements
!   INTEGER (kind=4), INTENT(IN) :: nnode      !number of nodes/element
!   INTEGER (kind=4), INTENT(IN) :: lnods(nnode,nelem) !connectivities
!   INTEGER (kind=4), POINTER :: bound(:,:)
!   INTEGER (kind=4), ALLOCATABLE :: side(:)
!   INTEGER (kind=4), INTENT(OUT) :: ne
!   END SUBROUTINE surf05
! END INTERFACE

 !compute solid boundary ==> bound(2,ne)  ne=number of segments defining boundary
! CALL surf05 ( lnods, nnode, nelem, bound, ne)
PRINT *,'***    ERROR DE PROGRAMACION    ***'
PRINT *,'*** FALTA LA SUBRUTINA ''SURF05'' ***'
STOP
 ALLOCATE( norm(ndime,ne), angle(ne), side(ne) )
 !compute side normals ==> norm
 DO i=1,ne           !for each segment
   t = coord(:,bound(2,i)) - coord(:,bound(1,i))  !side vector
   n = (/ t(2), -t(1) /)                          !outward normal
   CALL vecuni(ndime,n,l)                         !unit vector
   norm(:,i) = n                                  !assign
 END DO
 !compute angles with previous segment
 is = 1                       !initializes number of physical sides
 t = norm(:,ne)               !normal of last segment
 DO i=1,ne                    !for each segment
   angle(i) =  DOT_PRODUCT( t,norm(:,i) )    !Angle Cosine
   IF( angle(i) < pi4 .AND. i > 1)is = is+1  !increase number of sides
   side(i) = is               !side
   t = norm(:,i)              !keep normal of previos segment
 END DO
 IF( angle(1) > pi4 )THEN     !correct if first was not a corner
   DO i=ne,1,-1               !downward loop on segments
     IF( side(i) /= is )EXIT  !end of side found, exit correction loop
     side(i) = 1              !side
   END DO
   is = is-1                  !decrement number of sides
 END IF
 ! count nunber of segments in each side
 ALLOCATE( ns(is) )
 ns = 0                       !initializes
 DO i=1,ne                    !for each segment
   j = side(i)                !side
   ns(j) = ns(j) + 1          !sum
 END DO
 ! discard the two sides with the lowest number of segment
 DO k=1,2                      !for each side to discard
   ms = MINLOC( ns(1:is+1-k) ) !side to discard
   DO i=1,ne                   !for each segment
     j = side(i)               !side number
     IF( j == ms(1) )THEN      !if segment belongs to side
       side(i) = 0             !set to zero
     ELSE IF( j > ms(1) )THEN  !if segment belongs to an upper side number
       side(i) = side(i) - 1   !diminish side number
     END IF
   END DO
   DO j=ms(1)+1,is             !compact array ns
     ns(j-1) = ns(j)
   END DO
 END DO
 is = is - 2                   !number of remaining sides
 IF( is > 2 )THEN              !compact sides if neccessary
   first = .TRUE.              !set to side 1
   chang = .TRUE.              !change when transversal cut found
   DO i=1,ne                   !for each segment
     IF( side(i) /= 0 )THEN    !a side segment
       IF( first )THEN         !if side 1
         side(i) = 1           !set to side 1
       ELSE                    !else side 2
         side(i) = 2           !set to side 2
       END IF
       chang = .TRUE.          !always in a side set to change
     ELSE IF( chang )THEN      !for a transversal cut, change
       first = .NOT.first      !change to the other side
       chang = .FALSE.         !set to not change until cut ended
     END IF
   END DO
 END IF
 n1 = COUNT( side == 1 )  !number of segments on 'TOP' surface
 n2 = COUNT( side == 2 )  !number of segments on 'BOTTON' surface
 ! compute coordinates of the center of each side
 ALLOCATE( top(ndime,2,n1), bot(ndime,2,n2) )
 m1 = 0           !initializes counter of segments on side 1
 m2 = 0           !initializes counter of segments on side 2
 DO i=1,ne        !for each segment
   IF( side(i) == 0 )CYCLE  !a transversal cut
   IF( side(i) == 1)THEN    !for top surface
     m1 = m1+1              !increase number of segments
     top(:,1,m1) = (coord(:,bound(2,i)) + coord(:,bound(1,i)))/2d0 !center
     top(:,2,m1) = norm(:,i)                                       !normal
   ELSE                     !for bottom surface
     m2 = m2+1              !increase number of segments
     bot(:,1,m2) = (coord(:,bound(2,i)) + coord(:,bound(1,i)))/2d0 !center
     bot(:,2,m2) = -norm(:,i)                                      !normal
   END IF
 END DO
 DEALLOCATE( norm,side,ns,angle,bound )  !release auxiliar memory
 ! for each element compute normal direction
 DO i=1,nelem               !for each element in the mesh
   ! compute center of element
   DO j=1,ndime
     t(j) = SUM(coord(j,lnods(:,i)))/nnode
   END DO
   ! search for the nearest side on both TOP and BOTTOM surfaces
   l1 = 1e9                 !initializes to a large number
   k1 = 0                   !initializes position of nearest top segment
   DO j=1,n1                !for each segment on top surface
     n = t - top(:,1,j)     !distance from element center to segment center
     CALL vecuni(ndime,n,l) !compute distance ==> l
     IF( l < l1 )THEN       !compare distances
       l1 = l               !new minimum distance
       k1 = j               !new closest segment
     END IF
   END DO
   l2 = 1e9                 !initializes to a large number
   k2 = 0                   !initializes position of nearest bottom segment
   DO j=1,n2                !for each segment on bottom surface
     n = t - bot(:,1,j)     !distance from element center to segment center
     CALL vecuni(ndime,n,l) !compute distance ==> l
     IF( l < l2 )THEN       !compare distances
       l2 = l               !new minimum distance
       k2 = j               !new closest segment
     END IF
   END DO
   n = top(:,2,k1)*l2 + bot(:,2,k2)*l1  !average thickness direction
   CALL vecuni(ndime,n,l)           !unit vector
   dirt(:,i) = REAL(n,kind=8)       !assign
 END DO
 DEALLOCATE (top,bot)       !release auxiliar arrays
 RETURN
 END SUBROUTINE cpdir3d
