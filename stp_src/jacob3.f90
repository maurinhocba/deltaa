SUBROUTINE jacob3(cartd,deriv,detj,elcod,nnode)
!*********************************************************************
!
!**** this SUBROUTINE evaluates the jacobian matrix and the cartesian
!     shape function derivatives
!
!*********************************************************************
IMPLICIT NONE

INTEGER (kind=4), PARAMETER :: ndime=2
INTEGER (kind=4) nnode
REAL    (kind=8) cartd(nnode,ndime),deriv(nnode,ndime),detj,      &
                 elcod(ndime,nnode)

REAL    (kind=8) jaci(ndime,ndime) ,jacm(ndime,ndime)

!     ***  create jacobian matrix jacm

jacm = MATMUL(elcod,deriv) !jacm(i,j) = d[x(i)]/d[xita(j)]

!     ***  calculate determinant and inverse of jacobian matrix

detj = jacm(1,1)*jacm(2,2) - jacm(1,2)*jacm(2,1)
!IF(detj <= 0) THEN
!  WRITE(*,"(//, '  PROGRAM halted in SUBROUTINE jacob3')")
!  WRITE(*,"(11x,' zero or negative area')")
!  WRITE(*,"(11x,' element number ',i5)") elem
!  WRITE(*,"(2e15.4)") elcod,detj
!  RETURN
!END IF
jaci(1,1) =  jacm(2,2)/detj       !jaci(i,j) = d[xita(i)]/d[x(j)]
jaci(1,2) = -jacm(1,2)/detj
jaci(2,1) = -jacm(2,1)/detj
jaci(2,2) =  jacm(1,1)/detj

!     ***  calculate cartesian derivatives

cartd = MATMUL(deriv,jaci)       ! cartd(i,j)= d[N(i)]/d[x(j)]
RETURN
END SUBROUTINE jacob3
