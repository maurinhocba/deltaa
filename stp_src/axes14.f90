 SUBROUTINE axes14(x,a,b,sides,t,angle,a2,locax)
 !***********************************************************************
 !
 !    this routine compute the element local axes system
 !    for the 3 node element, and for the adjacent elements
 !    (local x-axis is directed along fiber at an Angle with
 !    standard direction (intersection with X-Y plane)
 !***********************************************************************
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8),INTENT(IN) :: x(3,6), &  !nodal coordinates
                             angle      !angle between standard X1 and local X1
 REAL (kind=8),INTENT(OUT) :: t(3,3),a(3,4),b(3,4),a2
 LOGICAL, INTENT(IN) :: sides(*)
 INTEGER (kind=4), INTENT(IN) :: locax  !local system option

 ! local variables
 INTEGER (kind=4) ii,jj,i,j,k
 INTEGER (kind=4), SAVE :: kk(3,3) = RESHAPE((/ 4,3,2, 5,1,3, 6,2,1 /), (/3,3/) )
 REAL (kind=8) l1(3),l2(3),t1(3),t2(3),t3(3),ll1,ll2,lt, &
               cosa,sina,cosb,sinb,cosg,sing,area2(4)

 !*** evaluate the first two side vectors

 l1 = x(1:3,3) - x(1:3,2)                             !side 1
 l2 = x(1:3,1) - x(1:3,3)                             !side 2

 !*** evaluate the cross product => plane normal

 t3(1) = l1(2)*l2(3) - l1(3)*l2(2)                    !normal * area2
 t3(2) = l1(3)*l2(1) - l1(1)*l2(3)
 t3(3) = l1(1)*l2(2) - l1(2)*l2(1)

 area2(1) = SQRT(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3)) !computes twice area
 a2 = area2(1)

 t3 = t3/area2(1)                                     !t3 (unit length)

 SELECT CASE (locax)
 CASE (1)
   lt = (t3(2)*t3(2)+t3(3)*t3(3)) !component in th Y-Z plane
   IF( lt  < 1.0d-5) THEN         !If t3 is almost orthogonal to  Y-Z plane
     t2 = (/ -t3(3), 0d0, t3(1) /) !choose t2 orthogonal to global Y direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t3,t1)
   ELSE       !         SELECT local y=t1 in the global YZ plane
     t1 = (/ 0d0, -t3(3), t3(2)  /)
     t2 = (/ lt, -t3(2)*t3(1), -t3(3)*t3(1) /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 CASE (2)
   lt = (t3(3)*t3(3)+t3(1)*t3(1)) !component in th Z-X plane
   IF( lt  < 1.0d-5) THEN         !If t3 is almost orthogonal to  Z-Y plane
     t2 = (/ t3(2), -t3(1), 0d0 /) !choose t2 orthogonal to global Z direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t3,t1)
   ELSE       !         SELECT local z=t1 in the global ZX plane
     t1 = (/ t3(3), 0d0, -t3(1)  /)
     t2 = (/  -t3(1)*t3(2), lt, -t3(3)*t3(2) /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 CASE (3)
   lt = (t3(1)*t3(1)+t3(2)*t3(2)) !component in th X-Y plane
   IF( lt  < 1.0d-5) THEN         !If t3 is almost orthogonal to  X-Y plane
     t2 = (/ 0d0, t3(3), -t3(2) /) !choose t2 orthogonal to global X direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t3,t1)
   ELSE       !         SELECT local x=t1 in the global xy plane
     t1 = (/ -t3(2), t3(1) , 0d0 /)
     t2 = (/ -t3(1)*t3(3), -t3(2)*t3(3), lt /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 END SELECT

 cosa = COS(angle)                            !angle to compute
 sina = SIN(angle)                            !local X1 direction

 t(1:3,1) = t1*cosa + t2*sina                 !local X1 direction
 t(1:3,2) =-t1*sina + t2*cosa                 !local X2 direction
 t(1:3,3) = t3                                !Normal direction

 !*** find the local coordinates

 a(1,1) = l1(1)*t(1,1)+l1(2)*t(2,1)+l1(3)*t(3,1)      ! l1 . t1
 a(2,1) = l2(1)*t(1,1)+l2(2)*t(2,1)+l2(3)*t(3,1)      ! l2 . t1
 a(3,1) = -a(1,1)-a(2,1)
 b(1,1) = l1(1)*t(1,2)+l1(2)*t(2,2)+l1(3)*t(3,2)      ! l1 . t2
 b(2,1) = l2(1)*t(1,2)+l2(2)*t(2,2)+l2(3)*t(3,2)      ! l2 . t2
 b(3,1) = -b(1,1)-b(2,1)

 DO ii=1,3
   IF(.NOT.sides(ii))CYCLE                             !no adjacent elemt.
   jj = ii+1                                           !next element
   i = kk(1,ii)                                        !local first node
   j = kk(2,ii)                                        !local second node
   k = kk(3,ii)                                        !local third node
   l1 = x(1:3,k) - x(1:3,j)                            !side J-K  (I)
   l2 = x(1:3,i) - x(1:3,k)                            !side K-I  (J)

   ll1 = SQRT(l1(1)*l1(1)+l1(2)*l1(2)+l1(3)*l1(3))     !length of side (I)
   ll2 = SQRT(l2(1)*l2(1)+l2(2)*l2(2)+l2(3)*l2(3))     !length of side (J)

   t3(1) = l1(2)*l2(3) - l1(3)*l2(2)                   !normal of elem (I)
   t3(2) = l1(3)*l2(1) - l1(1)*l2(3)
   t3(3) = l1(1)*l2(2) - l1(2)*l2(1)
   area2(jj) = SQRT(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3))!Area of Elem (I)
   !t(1:3,i) = t3/area2(jj)                             !normal of elem (I)

   a(1,jj) = -a(ii,1)                                  !projec. of side (I)
   b(1,jj) = -b(ii,1)
   cosa = a(1,jj)/ll1                                  !angle of side (I)
   sina = b(1,jj)/ll1

   cosb = (l1(1)*l2(1)+l1(2)*l2(2)+l1(3)*l2(3))/ll1/ll2!angle between (I-J)
   sinb = area2(jj)/ll1/ll2

   cosg = cosa*cosb - sina*sinb                !angle of side (J) g = a+b
   sing = cosa*sinb + sina*cosb

   a(2,jj) = ll2*cosg                                  !projec. of side (J)
   b(2,jj) = ll2*sing
   a(3,jj) = -a(1,jj)-a(2,jj)                          !projec. of side (K)
   b(3,jj) = -b(1,jj)-b(2,jj)

 END DO
 DO i=1,4
   IF(i > 1 )THEN
     IF( sides(i-1) )THEN
       a(1:3,i) = 0d0
       b(1:3,i) = 0d0
     END IF
   ELSE
     a(1:3,i) = a(1:3,i)/area2(i)
     b(1:3,i) = b(1:3,i)/area2(i)
    END IF
 END DO
 RETURN
 END SUBROUTINE axes14
