 SUBROUTINE axes25(x0,x,a,b,t,angle,a2)
 !***********************************************************************
 !
 !    this routine compute the element local axes system
 !    for the 3 node element, and for the adjacent elements
 !    (local x-axis is directed along fiber at an Angle with
 !    standard direction (intersection with X-Y plane)
 !***********************************************************************
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8),INTENT(IN) :: x0(3,4), &  !original nodal coordinates
                             x(3,4),  &  !actual nodal coordinates
                             angle       !angle between standard X1 and local X1
 REAL (kind=8),INTENT(OUT) :: t(3,2),a(4),b(4),a2

 ! local variables
 REAL (kind=8) l1(3),l2(3),t1(3),t2(3),t3(3),lt,cosa,sina,jac(2,2),ji(2,2)
 REAL (kind=8), PARAMETER :: deriv(4,2) = RESHAPE( &  !shape function derivatives
   (/-0.25d0, 0.25d0, 0.25d0,-0.25d0,-0.25d0,-0.25d0, 0.25d0, 0.25d0 /), &
   (/4,2/))

  !*** evaluate the local system at the center

  l1 = MATMUL(x0,deriv(:,1))  !x_xita_0
  l2 = MATMUL(x0,deriv(:,2))  !x_eta_0

  !*** evaluate the cross product => plane normal

  CALL vecpro(l1,l2,t3)               !normal direction
  a2 = SQRT(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3)) !length = area/4
  t3 = t3/a2                                     !t3 (unit length)

  lt = (t3(1)*t3(1)+t3(2)*t3(2)) !component in th X-Y plane
  IF( lt  < 1.0d-5) THEN         !If t3 is almost orthogonal to  X-Y plane
    t1 = (/  1d0, 0d0, 0d0 /)
    IF(t3(3) > 0d0 )THEN         !if t3 = Z
      t2 = (/  0d0, 1d0, 0d0 /)
    ELSE                         !if t3 = -Z
      t2 = (/  0d0,-1d0, 0d0 /)
    END IF
  ELSE
    !         SELECT local x=t1 in the global xy plane
    t1 = (/ -t3(2), t3(1) , 0d0 /)
    !         of course local y = t2 = t3 x t1
    t2 = (/ -t3(1)*t3(3), -t3(2)*t3(3), lt /)
    !           normalizes t1 & t2
    CALL vecuni(3,t1,lt)
    CALL vecuni(3,t2,lt)
  END IF

  cosa = COS(angle)                            !angle to compute
  sina = SIN(angle)                            !local X1 direction

  t(1:3,1) = t1*cosa + t2*sina                 !local X1 direction
  t(1:3,2) =-t1*sina + t2*cosa                 !local X2 direction

  !*** find the cartesian derivative

  jac(1,1) = DOT_PRODUCT(t(:,1),l1)         !x1,xi
  jac(1,2) = DOT_PRODUCT(t(:,1),l2)         !x1,eta
  jac(2,1) = DOT_PRODUCT(t(:,2),l1)         !x2,xi
  jac(2,2) = DOT_PRODUCT(t(:,2),l2)         !x2,eta
  lt = jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)     !determinant = area(0)
  ji(1,1)  =  jac(2,2)/lt                      !xi,x1
  ji(1,2)  = -jac(1,2)/lt                      !xi,x2
  ji(2,1)  = -jac(2,1)/lt                      !eta,x1
  ji(2,2)  =  jac(1,1)/lt                      !eta,x2
  a(:) = MATMUL(deriv(:,:),ji(:,1))            !N(j)_x1
  b(:) = MATMUL(deriv(:,:),ji(:,2))            !N(j)_x2
  a2 = 4d0*a2      !element area
  t(:,1) = MATMUL(x,a)
  t(:,2) = MATMUL(x,b)
 RETURN
 END SUBROUTINE axes25
