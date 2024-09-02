      SUBROUTINE csdata(maxve,iprint,label,npoin,coord,surf) !,esets,eset)

!.... read contact surface data and generate database

      USE cont1_db

      IMPLICIT NONE
!     arguments
      INTEGER (kind=4), INTENT(IN) :: iprint,npoin,label(:) !,esets,eset(esets)
      INTEGER (kind=4), INTENT(IN) :: maxve
      REAL (kind=8), INTENT(IN) :: coord(:,:)
      TYPE (surf1_db), POINTER :: surf  !INTENT(OUT)
      END SUBROUTINE csdata
