SUBROUTINE stra14(a,b,x,sides,stran)
!
!     Compute first and second fundamental forms for element BST
!
IMPLICIT NONE
REAL (kind=8), INTENT(IN) :: b(3,4),x(3,6),a(3,4)
REAL (kind=8), INTENT(OUT) :: stran(6)
LOGICAL, INTENT(IN) :: sides(3)

INTEGER (kind=4) :: i,j,k,n,ns,nc
REAL (kind=8) :: lb,ls(3),t(3,3),tr(3,2),fs(3),no(2,2),det,nd(2,2)

INTEGER (kind=4), SAVE ::  &
         kk(3,3) = RESHAPE((/ 4,3,2, 5,1,3, 6,2,1 /), (/3,3/) ), &
         nn(2,3) = RESHAPE((/ 2,3, 3,1, 1,2 /), (/2,3/) )

ns = 0
DO i=1,4
  IF( i == 1)THEN
    t(1:3,1) = -MATMUL(x(1:3,1:3),b(1:3,1))     ! x(1)
    t(1:3,2) = +MATMUL(x(1:3,1:3),a(1:3,1))     ! x(2)
    t(1,3) = t(2,1)*t(3,2) - t(3,1)*t(2,2)      !normal * rA
    t(2,3) = t(3,1)*t(1,2) - t(1,1)*t(3,2)
    t(3,3) = t(1,1)*t(2,2) - t(2,1)*t(1,2)

    !computes twice area
    lb = 1d0/SQRT(t(1,3)*t(1,3)+t(2,3)*t(2,3)+t(3,3)*t(3,3))
    t(1:3,3) = t(1:3,3)*lb                   ! normalizes normal vector

  ELSE
    n = i-1
    ls(n) = SQRT( b(n,1)**2 + a(n,1)**2 )  !side length
    IF(sides(n))THEN                           !if side element exists
      ns = ns+1                                !number of sides
      tr = 0d0                                 !initializes
      DO j=1,3                                 !for each node in the el.
        k = kk(j,n)                            !associated node
        tr(1:3,1) = tr(1:3,1) - b(j,i)*x(1:3,k)   !x(1)(i)
        tr(1:3,2) = tr(1:3,2) + a(j,i)*x(1:3,k)   !x(2)(i)
      END DO
      ! computes gradient along the normal and proyection over t3
      tr(:,1) = ( b(n,1)*tr(:,1) - a(n,1)*tr(:,2) )/ls(n)
      fs(n)   = DOT_PRODUCT(tr(:,1),t(:,3))
    ELSE
      nc = n
    END IF
  END IF
END DO

! Membrane part
stran(1) = DOT_PRODUCT(t(1:3,1),t(1:3,1))
stran(2) = DOT_PRODUCT(t(1:3,2),t(1:3,2))
stran(3) = DOT_PRODUCT(t(1:3,1),t(1:3,2))

! Bending part
SELECT CASE (ns)
CASE (1)
  DO i=1,3
    IF( sides(i) )THEN
      nc = i
      EXIT
    END IF
  END DO
  DO i=1,3
    IF( sides(i) )CYCLE
    fs(i) = -fs(nc)*(b(nc,1)*b(i,1)+a(nc,1)*a(i,i))/ls(nc)/ls(i)
  END DO
CASE (2)
  no(1:2,1) = (/ b(nn(1,nc),1), -a(nn(1,nc),1) /)/ls(nn(1,nc))
  no(1:2,2) = (/ b(nn(2,nc),1), -a(nn(2,nc),1) /)/ls(nn(2,nc))
  det = no(1,1)*no(2,2) - no(1,2)*no(2,1)
  nd(1,1) =  no(2,2)/det
  nd(2,1) = -no(1,2)/det
  nd(1,2) = -no(2,1)/det
  nd(2,2) =  no(1,1)/det
  fs(nc) = (-b(nc,1)*nd(1,1)+a(nc,1)*nd(2,1))/ls(nc)*fs(nn(1,nc)) &
         + (-b(nc,1)*nd(1,2)+a(nc,1)*nd(2,2))/ls(nc)*fs(nn(2,nc))
END SELECT
stran(4:6) = 0d0
DO i=1,3
 ! IF(sides(i))THEN      !side element exist
    stran(4) = stran(4) - b(i,1)**2*fs(i)/ls(i)         !k11
    stran(5) = stran(5) - a(i,1)**2*fs(i)/ls(i)         !k22
    stran(6) = stran(6) + a(i,1)*b(i,1)*fs(i)/ls(i)     !k12
 ! ELSE
 !   !nothing yet
 ! END IF
END DO

RETURN
END SUBROUTINE stra14
