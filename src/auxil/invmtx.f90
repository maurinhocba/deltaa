SUBROUTINE invmtx(a,b,deter,n)
!***********************************************************************
!
!***this routine inverts a square matrix "a"
!
!     a  -  given n*n matrix
!     b  -  inverse of "a"
!     n  -  matrix size
!
!***********************************************************************
IMPLICIT NONE
INTEGER (kind=4),INTENT(IN) :: n
REAL (kind=8),INTENT(IN OUT) :: a(n,n)
REAL (kind=8),INTENT(OUT):: b(n,n),deter

REAL (kind=8) :: denom,t1,t2,t3,t4
integer :: i,j,k,l,m,irow
real:: big,dum

SELECT CASE (n)
CASE (1)
  !***inverse of a 1*1 matrix
  deter=a(1,1)
  IF(deter == 0d0) RETURN
  b(1,1) = 1d0/a(1,1)

CASE (2)
  !***inverse of a 2*2 matrix
  deter = a(1,1)*a(2,2)-a(2,1)*a(1,2)
  IF(deter == 0d0) RETURN
  b(1,1) = a(2,2)/deter
  b(2,2) = a(1,1)/deter
  b(2,1) =-a(2,1)/deter
  b(1,2) =-a(1,2)/deter

CASE (3)
  !*** inverse of a 3*3 matrix
  t1  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
  t2  =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
  t3  = a(2,1)*a(3,2) - a(3,1)*a(2,2)
  deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
  IF(deter == 0d0) RETURN
  b(1,1) = t1/deter
  b(2,1) = t2/deter
  b(3,1) = t3/deter
  b(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))/deter
  b(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))/deter
  b(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))/deter
  b(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))/deter
  b(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))/deter
  b(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))/deter

CASE (4)
  !*** inverse of a 4*4 matrix
  t1= a(2,2)*a(3,3)*a(4,4) + a(2,3)*a(3,4)*a(4,2) +               &
      a(2,4)*a(3,2)*a(4,3) - a(2,3)*a(3,2)*a(4,4) -               &
      a(2,2)*a(3,4)*a(4,3) - a(2,4)*a(3,3)*a(4,2)
  t2=-a(2,1)*a(3,3)*a(4,4) - a(2,3)*a(3,4)*a(4,1) -               &
      a(2,4)*a(3,1)*a(4,3) + a(2,4)*a(3,3)*a(4,1) +               &
      a(2,3)*a(3,1)*a(4,4) + a(2,1)*a(3,4)*a(4,3)
  t3=+a(2,1)*a(3,2)*a(4,4) + a(2,2)*a(3,4)*a(4,1) +               &
      a(2,4)*a(3,1)*a(4,2) - a(2,4)*a(3,2)*a(4,1) -               &
      a(2,2)*a(3,1)*a(4,4) - a(2,1)*a(3,4)*a(4,2)
  t4=-a(2,1)*a(3,2)*a(4,3) - a(2,2)*a(3,3)*a(4,1) -               &
      a(2,3)*a(3,1)*a(4,2) + a(2,3)*a(3,2)*a(4,1) +               &
      a(2,2)*a(3,1)*a(4,3) + a(2,1)*a(3,3)*a(4,2)
  deter= a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3 + a(1,4)*t4
  IF(deter == 0d0) RETURN
  denom = 1.d0/deter
  b(1,1) = t1*denom
  b(2,1) = t2*denom
  b(3,1) = t3*denom
  b(4,1) = t4*denom
  b(1,2) =(- a(1,2)*a(3,3)*a(4,4) - a(1,3)*a(3,4)*a(4,2)          &
         - a(1,4)*a(3,2)*a(4,3) + a(1,3)*a(3,2)*a(4,4)            &
         + a(1,2)*a(3,4)*a(4,3) + a(1,4)*a(3,3)*a(4,2))*denom
  b(2,2) =(  a(1,1)*a(3,3)*a(4,4) + a(1,3)*a(3,4)*a(4,1)          &
           + a(1,4)*a(3,1)*a(4,3) - a(1,4)*a(3,3)*a(4,1)          &
           - a(1,3)*a(3,1)*a(4,4) - a(1,1)*a(3,4)*a(4,3))*denom
  b(3,2) =(- a(1,1)*a(3,2)*a(4,4) - a(1,2)*a(3,4)*a(4,1)          &
           - a(1,4)*a(3,1)*a(4,2) + a(1,4)*a(3,2)*a(4,1)          &
           + a(1,2)*a(3,1)*a(4,4) + a(1,1)*a(3,4)*a(4,2))*denom
  b(4,2) =(  a(1,1)*a(3,2)*a(4,3) + a(1,2)*a(3,3)*a(4,1)          &
           + a(1,3)*a(3,1)*a(4,2) - a(1,3)*a(3,2)*a(4,1)          &
           - a(1,2)*a(3,1)*a(4,3) - a(1,1)*a(3,3)*a(4,2))*denom
  b(1,3) =(  a(1,2)*a(2,3)*a(4,4) + a(1,3)*a(2,4)*a(4,2)          &
           + a(1,4)*a(2,2)*a(4,3) - a(1,3)*a(2,2)*a(4,4)          &
           - a(1,2)*a(2,4)*a(4,3) - a(1,4)*a(2,3)*a(4,2))*denom
  b(2,3) =(- a(1,1)*a(2,3)*a(4,4) - a(1,3)*a(2,4)*a(4,1)          &
           - a(1,4)*a(2,1)*a(4,3) + a(1,4)*a(2,3)*a(4,1)          &
           + a(1,3)*a(2,1)*a(4,4) + a(1,1)*a(2,4)*a(4,3))*denom
  b(3,3) =(  a(1,1)*a(2,2)*a(4,4) + a(1,2)*a(2,4)*a(4,1)          &
           + a(1,4)*a(2,1)*a(4,2) - a(1,4)*a(2,2)*a(4,1)          &
           - a(1,2)*a(2,1)*a(4,4) - a(1,1)*a(2,4)*a(4,2))*denom
  b(4,3) =(- a(1,1)*a(2,2)*a(4,3) - a(1,2)*a(2,3)*a(4,1)          &
           - a(1,3)*a(2,1)*a(4,2) + a(1,3)*a(2,2)*a(4,1)          &
           + a(1,2)*a(2,1)*a(4,3) + a(1,1)*a(2,3)*a(4,2))*denom
  b(1,4) =(- a(1,2)*a(2,3)*a(3,4) - a(1,3)*a(2,4)*a(3,2)          &
           - a(1,4)*a(2,2)*a(3,3) + a(1,4)*a(2,3)*a(3,2)          &
           + a(1,3)*a(2,2)*a(3,4) + a(1,2)*a(2,4)*a(3,3))*denom
  b(2,4) =(  a(1,1)*a(2,3)*a(3,4) + a(1,3)*a(2,4)*a(3,1)          &
           + a(1,4)*a(2,1)*a(3,3) - a(1,4)*a(2,3)*a(3,1)          &
           - a(1,3)*a(2,1)*a(3,4) - a(1,1)*a(2,4)*a(3,3))*denom
  b(3,4) =(- a(1,1)*a(2,2)*a(3,4) - a(1,2)*a(2,4)*a(3,1)          &
           - a(1,4)*a(2,1)*a(3,2) + a(1,4)*a(2,2)*a(3,1)          &
           + a(1,2)*a(2,1)*a(3,4) + a(1,1)*a(2,4)*a(3,2))*denom
  b(4,4) =(  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)          &
           + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)          &
           - a(1,2)*a(2,1)*a(3,3) - a(1,1)*a(2,3)*a(3,2))*denom
  CASE DEFAULT
    ! subroutine to calculate the inverse of a matrix using Gauss-Jordan elimination
    ! the inverse of matrix a(n,n) is calculated and stored in the matrix b(n,n)

    !build the identity matrix
    b = 0d0
    DO i = 1,n
      b(i,i) = 1.0d0
    END DO

    DO i = 1,n ! this is the big loop over all the columns of a(n,n)
               ! in case the entry a(i,i) is zero, we need to find a good pivot; this pivot
               ! is chosen as the largest value on the column i from a(j,i) with j = 1,n
      big = a(i,i)
      DO j = i,n
        IF (ABS(a(j,i)) >= ABS(big)) THEN
          big = a(j,i)
          irow = j
        END IF
      END DO
      ! interchange lines i with irow for both a() and b() matrices
      IF (ABS(big) >= ABS(a(i,i))) THEN
        DO k = 1,n
          dum = a(i,k) ! matrix a()
          a(i,k) = a(irow,k)
          a(irow,k) = dum
          dum = b(i,k) ! matrix b()
          b(i,k) = b(irow,k)
          b(irow,k) = dum
        END DO
      END IF
      ! divide all entries in line i from a(i,j) by the value a(i,i);
      ! same operation for the identity matrix
      dum = a(i,i)
      DO j = 1,n
        a(i,j) = a(i,j)/dum
        b(i,j) = b(i,j)/dum
      END DO
      ! make zero all entries in the column a(j,i); same operation for indent()
      DO j = i+1,n
        dum = a(j,i)
        DO k = 1,n
          a(j,k) = a(j,k) - dum*a(i,k)
          b(j,k) = b(j,k) - dum*b(i,k)
        END DO
      END DO
    END DO

    ! substract appropiate multiple of row j from row j-1
    DO i = 1,n-1
      DO j = i+1,n
        dum = a(i,j)
        DO l = 1,n
          a(i,l) = a(i,l)-dum*a(j,l)
          b(i,l) = b(i,l)-dum*b(j,l)
        END DO
      END DO
    END DO
  END SELECT

RETURN

END SUBROUTINE invmtx
