FUNCTION functs (iload,ttime)
!********************************************************************
!
!***  heaviside(1), harmonic(2), multi-linear, etc.  time FUNCTIONs
!
!OUTPUT: functs(1) = function value,  functs(2) = function derivative
!********************************************************************
IMPLICIT NONE
REAL (kind=8) :: functs(2)
INTEGER (kind=4),INTENT(IN) :: iload
REAL (kind=8),INTENT(IN) :: ttime
END FUNCTION functs
