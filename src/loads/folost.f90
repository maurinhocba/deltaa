 SUBROUTINE folost (ttime,nsymm,force,stiff,ustif,coora,headf,factor,acurv)
 !     applies follower load
 USE loa_db, ONLY: foll_seg
 USE ctrl_db, ONLY : ndime,ntype
 IMPLICIT NONE

   !--- Dummy variables
   REAL (kind=8), INTENT(IN) :: ttime,     & !to compute load factor
                                force(:),  & !assembled forces
                                stiff(:),  & !stiffness matrix (symmetric)
                                ustif(:),  & !stiffness matrix (antisymmetric)
                                coora(:,:),& !present coordinates
                                factor       !load factor
   INTEGER(kind=4):: nsymm, & !0: symmetric matrix only, 1:include unsymmetric part
                     acurv    !assigned curve
   TYPE(foll_seg),POINTER:: headf            !pointer to first segment

   !--- Local variables
   REAL (kind=8), PARAMETER :: pi=3.141592653589793
   INTEGER(kind=4):: n,inode,nnode,onode,nv,ns,na
   REAL(kind=8):: press,x(ndime,4),t1(ndime),t2(ndime),t3(ndime),fs(2),ro
   REAL(kind=8), ALLOCATABLE :: st(:),as(:)
   TYPE (foll_seg),POINTER :: seg

   INTERFACE
     INCLUDE 'functs.h'
   END INTERFACE
   !---------------------------------------------------------------

   fs(1:2) = functs(acurv,ttime)*factor  !present value of load function
   IF( fs(1) == 0d0 )RETURN
   IF( ndime == 2 .AND. ntype == 3 )fs(1) = fs(1)*pi  !for axilsymmetric problems
   onode = 0              !initializes
   seg => headf           !point to first segment in the list
   DO             !loop over all the segments
     IF (.NOT.ASSOCIATED (seg) ) EXIT   !if last segment processed, EXIT
     press = -seg%fload*fs(1)    !change sign for stiffness matrix
     nnode = seg%nnode           !number of nodes in the segment
     IF( nnode /= onode )THEN  !if nnode have changed
       nv = nnode*ndime          !number of DOF
       ns = (nv*(nv+1))/2        !size of symetric stiffness
       na = ns-nv                !size of antisymetric stiffness
       IF( ALLOCATED(st) ) DEALLOCATE(st,as) !check
       ALLOCATE(st(ns),as(na))   !get memory for matrices
       st = 0d0                  !initializes
       as = 0d0
       onode = nnode             !keep number of nodes
     END IF
     x(1:ndime,1:nnode) = coora(1:ndime,seg%lnofl(1:nnode))  !nodal coordinates
     IF ( nnode == 2) THEN  !2D problems (2-node segments)
       press = press/2d0
       !  symmetric part for axilsymmetric problems only
       IF(ntype == 3)THEN            !for axilsymmetric problems
         t2(1) = (x(2,1)-x(2,2))*press        ! t1 * L  normal vector x length
         t2(2) = (x(1,2)-x(1,1))*press/2d0    ! t2 * L/2
         ro = (x(1,1)+x(1,2))*press           ! r1+r2
         st = (/ t2(1),t2(2),t2(1)   ,t2(2)-ro, &
                         0d0,t2(2)+ro,     0d0, &
                                t2(1),   t2(2), &
                                           0d0 /)

         !  antisymmetric part
         IF( nsymm == 1 ) &
           as = (/ t2(2)-ro,   0d0,    t2(2), &
                            -t2(2),      0d0, &
                                    t2(2)+ro /)
       ELSE
         !st = (/   0d0,   0d0,   0d0,-press, &
         !                 0d0, press,   0d0, &
         !                        0d0,   0d0, &
         !                               0d0 /)
         st(4) = -press; st(6) = press
         !  antisymmetric part
         IF( nsymm == 1 ) THEN
         !  as = (/ -press,  0d0,   0d0, &
         !                   0d0,   0d0, &
         !                        press /)
           as(1) = -press; as(6) = press

         END IF
       END IF

     ELSE IF ( nnode == 3) THEN !3D with triangles
       press = press/12d0
       t1 = (x(1:3,3) - x(1:3,2))*press    ! side 1
       t2 = (x(1:3,1) - x(1:3,3))*press    ! side 2
       t3 = (x(1:3,2) - x(1:3,1))*press    ! side 3
       x(:,1) = t1-t2
       x(:,2) = t1-t3
       x(:,3) = t2-t3
       st = (/ 0d0,0d0,0d0,    0d0,+x(3,1),-x(2,1),    0d0,+x(3,2),-x(2,2), &
                   0d0,0d0,-x(3,1),    0d0,+x(1,1),-x(3,2),    0d0,+x(1,2), &
                       0d0,+x(2,1),-x(1,1),    0d0,+x(2,2),-x(1,2),    0d0, &
                               0d0,    0d0,    0d0,    0d0,+x(3,3),-x(2,3), &
                                       0d0,    0d0,-x(3,3),    0d0,+x(1,3), &
                                               0d0,+x(2,3),-x(1,3),    0d0, &
                                                       0d0,    0d0,    0d0, &
                                                               0d0,    0d0, &
                                                                       0d0 /)
       IF( nsymm == 1 ) THEN
         t1 = 2d0*t1
         t2 = 2d0*t2
         t3 = 2d0*t3
         as = (/ +t1(3),-t1(2),   0d0,+t3(3),-t3(2),   0d0,+t2(3),-t2(2), &
                        +t1(1),-t3(3),   0d0,+t3(1),-t2(3),   0d0,+t2(1), &
                               +t3(2),-t3(1),   0d0,+t2(2),-t2(1),   0d0, &
                                      +t2(3),-t2(2),   0d0,+t1(3),-t1(2), &
                                             +t2(1),-t1(3),   0d0,+t1(1), &
                                                    +t1(2),-t1(1),   0d0, &
                                                           +t3(3),-t3(2), &
                                                                  +t3(1) /)
       END IF

     ELSE IF (nnode == 4) THEN  !3D with quads
       press = press/4d0
       t1 = (-x(:,1)+x(:,3))*press  !d1     *p/4
       t2 = (-x(:,2)+x(:,4))*press  !d2     *p/4
       x(:,1) = t1+t2               !d1+d2  *p/4
       x(:,2) =-t1+t2               !-d1+d2 *p/4
       t1 = 2d0*t1                  !2*d1   *p/4
       t2 = 2d0*t2                  !2*d2   *p/4
       st = (/ 0d0,0d0,0d0,    0d0,+x(3,1),-x(2,1),    0d0,+ t2(3),- t2(2),    0d0,+x(3,2),-x(2,2), &
                   0d0,0d0,-x(3,1),    0d0,+x(1,1),- t2(3),    0d0,+ t2(1),-x(3,2),    0d0,+x(1,2), &
                       0d0,+x(2,1),-x(1,1),    0d0,+ t2(2),- t2(1),    0d0,+x(2,2),-x(1,2),    0d0, &
                               0d0,    0d0,    0d0,    0d0,+x(3,2),-x(2,2),    0d0,- t1(3),+ t1(2), &
                                       0d0,    0d0,-x(3,2),    0d0,+x(1,2),+ t1(3),    0d0,- t1(1), &
                                               0d0,+x(2,2),-x(1,2),    0d0,- t1(2),+ t1(1),    0d0, &
                                                       0d0,    0d0,    0d0,    0d0,-x(3,1),+x(2,1), &
                                                               0d0,    0d0,+x(3,1),    0d0,-x(1,1), &
                                                                       0d0,-x(2,1),+x(1,1),    0d0, &
                                                                               0d0,    0d0,    0d0, &
                                                                                       0d0,    0d0, &
                                                                                               0d0/)
       IF( nsymm == 1 ) as =(/ &
           + t2(3),- t2(2),    0d0,+x(3,2),-x(2,2),    0d0,    0d0,    0d0,    0d0,+x(3,1),-x(2,1), &
                   + t2(1),-x(3,2),    0d0,+x(1,2),    0d0,    0d0,    0d0,-x(3,1),    0d0,+x(1,1), &
                           +x(2,2),-x(1,2),    0d0,    0d0,    0d0,    0d0,+x(2,1),-x(1,1),    0d0, &
                                   - t1(3),+ t1(2),    0d0,-x(3,1), x(2,1),    0d0,    0d0,    0d0, &
                                           - t1(1), x(3,1),    0d0,-x(1,1),    0d0,    0d0,    0d0, &
                                                   -x(2,1), x(1,1),    0d0,    0d0,    0d0,    0d0, &
                                                           - t2(3),+ t2(2),    0d0,-x(3,2),+x(2,2), &
                                                                   - t2(1),+x(3,2),    0d0,-x(1,2), &
                                                                           -x(2,2),+x(1,2),    0d0, &
                                                                                   + t1(3),- t1(2), &
                                                                                           + t1(1) /)
     END IF
     ! add into stiffness matrix
     CALL addlhs(seg%lnofl(1),ndime,nnode,force(1),stiff(1),ustif(1),nsymm,st(1),as(1),nv)
     seg => seg%next
   END DO
   DEALLOCATE(st,as)

 RETURN
 END SUBROUTINE folost
