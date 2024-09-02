SUBROUTINE axes12(x,a,b,sides,tt,a2)
!***********************************************************************
!
!    this routine compute the element local axes system
!    for the 3 node element, and for the adjacent elements
!    (local x-axis is directed from node 1 to opposite mid-side
!***********************************************************************
IMPLICIT NONE

REAL (kind=8),INTENT(IN) :: x(3,6)    !nodal coordinates
REAL (kind=8),INTENT(OUT) :: a(3,4),b(3,4),tt(3,3),a2
LOGICAL, INTENT(IN) :: sides(3)


INTEGER (kind=4) ii,jj,i,j,k
INTEGER (kind=4), SAVE :: &
         kk(3,3) = RESHAPE((/ 4,3,2, 5,1,3, 6,2,1 /), (/3,3/) )
REAL (kind=8) l1(3),l2(3),t3(3),t1(3),t2(3),ll1,ll2,lt, &
              cosa,sina,cosb,sinb,cosg,sing,area2(4)

!*** evaluate the local axes

t1 = x(1:3,2) - x(1:3,1)                             !side 3 ==> t1
lt = SQRT(DOT_PRODUCT(t1,t1))                        !length of l3
t1 = t1/lt                                           !t1 (unit length)
l1 = x(1:3,3) - x(1:3,2)                             !side 1
l2 = x(1:3,1) - x(1:3,3)                             !side 2
t2 = l1                                              !candidate for t2
t2 = t2 - DOT_PRODUCT(t1,t2)*t1                      !t2 direction
lt = SQRT(DOT_PRODUCT(t2,t2))                        !lenght of t2
t2 = t2/lt                                           !t2 (unit length)

!*** evaluate the cross product

t3(1) = l1(2)*l2(3) - l1(3)*l2(2)                    !normal * area2
t3(2) = l1(3)*l2(1) - l1(1)*l2(3)
t3(3) = l1(1)*l2(2) - l1(2)*l2(1)

area2(1) = SQRT(DOT_PRODUCT(t3,t3))                  !computes twice area
a2 = area2(1)
t3 = t3/area2(1)                                     !t3 (unit length)

tt(:,1) = t1
tt(:,2) = t2
tt(:,3) = t3

!*** find the local coordinates

a(1,1) = DOT_PRODUCT(l1,t1)                          ! l1 . t1
a(2,1) = DOT_PRODUCT(l2,t1)                          ! l2 . t1
a(3,1) = -a(1,1)-a(2,1)
b(1,1) = DOT_PRODUCT(l1,t2)                          ! l1 . t2
b(2,1) = DOT_PRODUCT(l2,t2)                          ! l2 . t2
b(3,1) = -b(1,1)-b(2,1)

DO ii=1,3
  IF(.NOT.sides(ii))CYCLE                             !no adjacent elemt.
  jj = ii+1                                           !next element
  i = kk(1,ii)                                        !local first node
  j = kk(2,ii)                                        !local second node
  k = kk(3,ii)                                        !local third node
  l1 = x(1:3,k) - x(1:3,j)                            !side J-K  (I)
  l2 = x(1:3,i) - x(1:3,k)                            !side K-I  (J)

  ll1 = SQRT(DOT_PRODUCT(l1,l1))                      !length of side (I)
  ll2 = SQRT(DOT_PRODUCT(l2,l2))                      !length of side (J)

  t3(1) = l1(2)*l2(3) - l1(3)*l2(2)                   !normal of elem (I)
  t3(2) = l1(3)*l2(1) - l1(1)*l2(3)
  t3(3) = l1(1)*l2(2) - l1(2)*l2(1)
  area2(jj) = SQRT(DOT_PRODUCT(t3,t3))                !Area of Elem (I)

  a(1,jj) = -a(ii,1)                                  !projec. of side (I)
  b(1,jj) = -b(ii,1)
  cosa = a(1,jj)/ll1                                  !angle of side (I)
  sina = b(1,jj)/ll1

  cosb = DOT_PRODUCT(l1,l2)/ll1/ll2                   !angle between (I-J)
  sinb = area2(jj)/ll1/ll2

  cosg = cosa*cosb - sina*sinb                !angle of side (J) g = a+b
  sing = cosa*sinb + sina*cosb

  a(2,jj) = ll2*cosg                                  !projec. of side (J)
  b(2,jj) = ll2*sing
  a(3,jj) = -a(1,jj)-a(2,jj)                          !projec. of side (K)
  b(3,jj) = -b(1,jj)-b(2,jj)

END DO
a(1:3,1) = a(1:3,1)/area2(1)
b(1:3,1) = b(1:3,1)/area2(1)
DO jj=2,4
  IF( sides(jj-1) )THEN
    a(1:3,jj) = a(1:3,jj)/area2(jj)
    b(1:3,jj) = b(1:3,jj)/area2(jj)
  END IF
END DO

RETURN
END SUBROUTINE axes12
