SUBROUTINE csrfc1(ipair,coors,coorm,pair,master,slave,cmptf,emass,fcont)

!.... searching procedure for a given master-slave 2-d surfaces pair

USE cont1_db
IMPLICIT NONE
!     Dummy arguments
LOGICAL, INTENT(IN) :: cmptf           !compute total contact force?
INTEGER (kind=4), INTENT(IN) :: ipair  !pair position
REAL (kind=8), INTENT(IN) :: emass(:,:), &! nodal masses
                             coors(:,:), &! coordinates of slave surface
                             coorm(:,:)   ! coordinates of master surface
REAL (kind=8), INTENT(IN OUT) :: fcont(:,:)  !contact forces
TYPE (pair1_db), POINTER :: pair             !contact pair
TYPE (surf1_db), POINTER :: master,slave     !contact surfaces

END SUBROUTINE csrfc1
