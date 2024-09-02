 SUBROUTINE shape4 (nnode, shape, deriv ,xita, eta, zeta, bezier, order )
 !********************************************************************
 !
 !***  calculates shape functions and their derivatives for 3d
 !     prismatic elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER(Kind=4) :: nnode
 REAL (kind=8)    deriv(nnode,3),shape(nnode),xita,eta,zeta
 LOGICAL, INTENT(IN), OPTIONAL :: bezier,order

 END SUBROUTINE shape4
