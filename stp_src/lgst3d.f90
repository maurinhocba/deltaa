SUBROUTINE lgst3d(u2,lb)
!
!  compute log-strains for a 3d problem
!
IMPLICIT NONE
REAL(kind=8), INTENT(IN)  :: u2(6)  !U^2 in vector form (11-22-33-12-23-31)
REAL(kind=8), INTENT(OUT) :: lb(3)  !princial log strains

 ! local variables
 INTEGER (kind=4) :: i
 REAL (kind=8) i1,i2,i3,sb,sbb,phi,a12,a13,a23,a11,a22,a33
 REAL (kind=8),PARAMETER :: s23  = 1.15470053837925, pi6  = 0.523598775598299,  &
                            pi23 = 2.09439510239320, pi43 = 4.18879020478639

 !
 !     compute principal values
 !

 a11 = u2(1) - 1d0  !to improve numerical accuracy for small strains
 a22 = u2(2) - 1d0
 a33 = u2(3) - 1d0
 a12 = u2(4)
 a13 = u2(5)
 a23 = u2(6)

 IF( ABS(a12) + ABS(a13) + ABS(a23) > 1d-16)THEN
   i1  = a11 + a22 + a33                                        !first inv.
   i2  = a11*a22 + a11*a33 + a22*a33 -a12**2 - a13**2 - a23**2  !second
   i3  = a11*a22*a33 + 2d0*a12*a13*a23              &           !third inv
        -a11*a23**2 - a22*a13**2 - a33*a12**2
   sb  = SQRT((i1*i1 - 3d0*i2)/3d0)
   sbb = i1*i2 - 2d0/9d0*i1**3 - 3d0*i3
   phi = sbb/(s23*sb**3)
   IF( phi > 1d0 )THEN         !check range for round errors
     phi = pi6
   ELSE IF (phi < -1d0) THEN
     phi = -pi6
   ELSE                        !normaly
     phi = ASIN(phi)/3d0
   END IF
   i1 = i1/3d0  !mean value
   sb = s23*sb
   lb(1) = sb*SIN(phi+pi23) + i1         !largest eigenvalue
   lb(2) = sb*SIN(phi)      + i1         !mediun eigenvalue
   lb(3) = sb*SIN(phi+pi43) + i1         !smallest eigenvalue
 ELSE
   !main directions concides with cartesian axes
   !lb = (/ a11, a22, a33 /)   !assume
   lb(1) = MAX( a11,a22,a33 )
   lb(3) = MIN( a11,a22,a33 )
   lb(2) = a11+a22+a33-lb(1)-lb(3)
 END IF
 !Compute eigenvalues and principal strains
 lb = lb + 1d0           !restore eigenvalues of U^2
 DO i=1,3
   IF( lb(i) > 0d0 )THEN
     lb(i) = SQRT(lb(i)) !previous values where eigenvalues of U^2
   ELSE
     WRITE(*,"('Negative eigenvalue found',e15.4)")lb(i)
     lb(i) = 6.737947d-3
   END IF
   lb(i) = LOG(lb(i))    !log strains
 END DO
 RETURN
END SUBROUTINE lgst3d
