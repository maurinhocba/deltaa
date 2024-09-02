SUBROUTINE jacob5(cartd,deriv,jac,elcod,nnode)
!*********************************************************************
!
!**** this SUBROUTINE evaluates the jacobian matrix and the cartesian
!     shape FUNCTION derivatives
!
!*********************************************************************
IMPLICIT NONE
INTEGER (kind=4), PARAMETER :: ndime = 3
INTEGER (kind=4), INTENT(IN) :: nnode   !number of nodes per element
REAL (kind=8), INTENT(IN)  :: deriv(nnode,ndime),elcod(ndime,nnode)
REAL (kind=8), INTENT(OUT) :: cartd(nnode,ndime),jac

REAL    (kind=8) jaci(ndime,ndime), jacm(ndime,ndime)

!     create jacobian matrix jacm

jacm = MATMUL(elcod,deriv)

!     calculate determinant and inverse of jacobian matrix

jac =  jacm(1,1)*jacm(2,2)*jacm(3,3)                              &
      +jacm(1,3)*jacm(2,1)*jacm(3,2)                              &
      +jacm(3,1)*jacm(1,2)*jacm(2,3)                              &
      -jacm(3,1)*jacm(2,2)*jacm(1,3)                              &
      -jacm(3,3)*jacm(1,2)*jacm(2,1)                              &
      -jacm(1,1)*jacm(2,3)*jacm(3,2)

jaci(1,1) =  (jacm(2,2)*jacm(3,3)-jacm(2,3)*jacm(3,2))/jac
jaci(2,1) = -(jacm(2,1)*jacm(3,3)-jacm(3,1)*jacm(2,3))/jac
jaci(3,1) =  (jacm(2,1)*jacm(3,2)-jacm(2,2)*jacm(3,1))/jac
jaci(1,2) = -(jacm(1,2)*jacm(3,3)-jacm(1,3)*jacm(3,2))/jac
jaci(2,2) =  (jacm(1,1)*jacm(3,3)-jacm(3,1)*jacm(1,3))/jac
jaci(3,2) = -(jacm(1,1)*jacm(3,2)-jacm(1,2)*jacm(3,1))/jac
jaci(1,3) =  (jacm(1,2)*jacm(2,3)-jacm(1,3)*jacm(2,2))/jac
jaci(2,3) = -(jacm(1,1)*jacm(2,3)-jacm(1,3)*jacm(2,1))/jac
jaci(3,3) =  (jacm(1,1)*jacm(2,2)-jacm(2,1)*jacm(1,2))/jac

!     calculate cartesian derivatives

cartd = MATMUL(deriv,jaci)

RETURN
END SUBROUTINE jacob5
