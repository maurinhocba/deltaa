 SUBROUTINE check_orientation8(lnods,x,t)
 ! check element connectivities
 IMPLICIT NONE
 !dummy arguments
 INTEGER(kind=4), INTENT(IN OUT) :: lnods(8) !connectivities
 REAL(kind=8), INTENT(IN OUT) :: x(3,8)      !coordinates
 REAL(kind=8), INTENT(IN) :: t(3)            !shell normal
 END SUBROUTINE check_orientation8
