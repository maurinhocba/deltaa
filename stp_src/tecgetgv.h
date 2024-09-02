 SUBROUTINE tecgetgv(flag,gv,nodes,nps,ngv)
 !USE data_db, ONLY : npoin,wsdisp,coord,coors
 IMPLICIT NONE
 ! dummy arguments
 LOGICAL, INTENT(IN) :: flag  !.FALSE. coordinates  .TRUE. all nodal variables
 INTEGER(kind=4), INTENT(IN) :: nps  !number of active nodes
 INTEGER(kind=4), POINTER :: nodes(:,:)  !number of active nodes
 INTEGER(kind=4), INTENT(OUT) :: ngv  !number of global variables
 REAL(kind=4), POINTER :: gv(:,:)     !global variables
 END SUBROUTINE tecgetgv
