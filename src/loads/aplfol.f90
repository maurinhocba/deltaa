 SUBROUTINE aplfol (ndime,coora,force,ntype,headf)
 !     applies follower load
 USE loa_db, ONLY: foll_seg
 USE npo_db, ONLY : ifpre
 IMPLICIT NONE

   !--- Dummy variables
   INTEGER(kind=4):: ndime, & !problem dimension
                     ntype    !problem type for 2D
   REAL(kind=8), INTENT (IN) :: coora(:,:)   !present coordinates
   REAL(kind=8), INTENT (IN OUT) :: force(:) !assembled forces
   TYPE(foll_seg),POINTER:: headf            !pointer to first segment
   !--- Local variables
   REAL (kind=8), PARAMETER :: pi=3.141592653589793
   INTEGER(kind=4):: n,j,nnode ,k(4)
   !REAL(kind=8):: pforc,press,x(ndime,4),t(ndime,2),vn(ndime),px(ndime) !px(ndime,4)
   REAL(kind=8):: pforc,press,x(ndime,4),t(ndime,2),vn(ndime),px(ndime,4)
   TYPE (foll_seg),POINTER :: seg
   !---------------------------------------------------------------
   IF(ndime == 2 .AND. ntype == 3 )THEN
     k = (/ 1,2,0,0 /)
   ELSE
     k = 1
   END IF

   seg => headf
   DO
     IF (.NOT.ASSOCIATED (seg) ) EXIT
     press = -seg%fload !change sign


     nnode = seg%nnode
     x(:,1:nnode) = coora(:,seg%lnofl(1:nnode))
     IF ( nnode == 2) THEN  !2D problems
       vn(1) =  x(2,1)-x(2,2)  !normal vector x length
       vn(2) =  x(1,2)-x(1,1)
       pforc = press/nnode
       IF(ntype == 3) THEN
         px(:,1) = pforc*pi*(2d0*x(1,1)+x(1,2))/1.5d0*vn
         px(:,2) = pforc*pi*(x(1,1)+2d0*x(1,2))/1.5d0*vn
       ELSE
         px(:,1) = pforc*vn
       END IF
     ELSE IF ( nnode == 3) THEN !3D with triangles
       t(:,1) = x(:,2) - x(:,1)    ! side 3
       t(:,2) = x(:,3) - x(:,1)    !-side 2
       !*** cross product  T3 = T1 x T2 = 2 area * t3
       vn(1) = t(2,1)*t(3,2) - t(3,1)*t(2,2)
       vn(2) = t(3,1)*t(1,2) - t(1,1)*t(3,2)
       vn(3) = t(1,1)*t(2,2) - t(2,1)*t(1,2)
       pforc = press/2d0/nnode
       px(:,1)  = pforc*vn
     ELSE IF (nnode == 4) THEN  !3D with quads
       t(:,1) = (-x(:,1)+x(:,3))  !d1
       t(:,2) = (-x(:,2)+x(:,4))  !d2
       !*** cross product  2A T = d1 x d2
       vn(1) = t(2,1)*t(3,2) - t(3,1)*t(2,2)
       vn(2) = t(3,1)*t(1,2) - t(1,1)*t(3,2)
       vn(3) = t(1,1)*t(2,2) - t(2,1)*t(1,2)
       pforc = press/2d0/nnode
       px(:,1)  = pforc*vn
     END IF

     DO j=1,nnode
       n = seg%lnofl(j) !internal nodes
       CALL ensvec(ndime,ifpre(1,n),px(1,k(j)),force(1))
     END DO
     seg => seg%next
   END DO


 RETURN
 END SUBROUTINE aplfol
