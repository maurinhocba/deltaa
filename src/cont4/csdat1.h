SUBROUTINE csdat1(ncnod,nsegm,ncnxx,nsexx,nnseg,lcnod,lcseg,iwrit, &
                  label,npoin,coord,sname)

!.... READ node numbers and segment connectivities

USE surf_db
IMPLICIT NONE
!     arguments
INTEGER (kind=4), INTENT(IN) :: nnseg,iwrit,npoin,                 &
                  label(:),ncnxx,nsexx
INTEGER (kind=4), INTENT(IN OUT) :: ncnod,nsegm
INTEGER (kind=4), INTENT(OUT) :: lcnod(:),lcseg(:,:)
REAL(kind=8), INTENT(IN) :: coord(:,:)
CHARACTER (len=6), INTENT(IN) :: sname
END SUBROUTINE csdat1
