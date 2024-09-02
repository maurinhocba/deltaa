SUBROUTINE lgst2d(stran,r1,r2,lb)

! compute log strains from C = U^2 = F^T F

IMPLICIT NONE
! dummy arguments
REAL(kind=8), INTENT (IN OUT) :: stran(3) !(IN) C, (OUT) log strains
REAL(kind=8), INTENT (OUT) :: r1,r2, & !components of first eigevector
                              lb(2)    !eigenvalues of U

! local variables
REAL (kind=8) :: c,r,l1,l2

!compute eigenvalues
c = (stran(1)+stran(2))/2d0                        !center of circle
r = SQRT((stran(1)-stran(2))**2/4d0+stran(3)**2)   !circle radius
l1 = c+r                                 !first (maximum) eigenvalue
l2 = c-r                                 !second (minimum) eigenvalue
!compute eigenvectors
c = stran(1) - l1                        !first diagonal element
IF( c < -1d-15)THEN                      !check
  r = -stran(3)/c                        ! v = (r,1)
  c = SQRT(r**2+1d0)                     ! eigenvector length
  r1 = r/c                               !first component of eigevector
  r2 = 1d0/c                             !first component of eigevector
ELSE
  r1 = 1d0                               !local direction is the vector
  r2 = 0d0
END IF
!Compute eigenvalues and principal strains
!IF( l2 < 0d0 )THEN
!  WRITE(55,*)" Too distorted mesh, negative eigen-value detected"
!  WRITE(55,*)' U^2 = ',stran
!  WRITE(55,*)' l1 = ',l1,' l2 = ',l2
!  WRITE(55,*)' CALLED FROM ',sname
!  ierr = 1
!  RETURN
!END IF
lb(1) = SQRT(l1)     !previous values where eigenvalues of U^2
lb(2) = SQRT(l2)
l1 = LOG(lb(1))      !ln of principal stretches
l2 = LOG(lb(2))
!Compute 2-D strains with the spectral decomposition
stran(1) = l1*r1*r1 + l2*r2*r2
stran(2) = l1*r2*r2 + l2*r1*r1
stran(3) = 2d0*(l1-l2)*r1*r2   !twice the strain
lb(1) = l1
lb(2) = l2
RETURN
END SUBROUTINE lgst2d
