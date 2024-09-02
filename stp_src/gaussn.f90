SUBROUTINE gaussn ( )
USE data_db
IMPLICIT NONE

INTEGER (kind=4) ::  etype, iset, &
                     spot_iset, truss_iset, sol2d_iset, sol3d_iset, shl3d_iset, &
                     beame_iset, shrev_iset, ele11_iset, bst_iset,   &
                     ele15_iset, spher_iset

spot_iset  = 0       !initializes set counters
truss_iset = 0       !initializes set counters
sol2d_iset = 0
sol3d_iset = 0
shl3d_iset = 0
beame_iset = 0
shrev_iset = 0
ele11_iset = 0
bst_iset   = 0
ele15_iset = 0
spher_iset = 0

DO iset=1,nsets

  etype = setyp(iset)     !element type

  SELECT CASE(etype)

  CASE (1)      ! skip Gauss data for TRUSS elements
    CALL gausn1 (spot_iset)

  CASE (2)      ! skip Gauss data for TRUSS elements
    CALL gausn2 (truss_iset)

  CASE (3,17,20)  ! skip Gauss data for 2-D SOLID elements
    CALL gausn3 (sol2d_iset)

  CASE (4,5,12,16,18)      ! skip Gauss data for 3-D SOLID elements
    CALL gausn5 (sol3d_iset)

  CASE (6,7)    ! skip Gauss data for 3-D SHELL (Simo Theory) elements
    CALL gausn6 (shl3d_iset)

  CASE (8)      ! skip Gauss data for 3-D BEAM (Simo Theory) element
    CALL gausn8 (shrev_iset)

  CASE (9,11)      ! skip Gauss data for 2-D SHELL/BEAM (Simo Theory) element
    CALL gausn9 (shrev_iset)

  CASE (13:15,25)  ! skip Gauss data for SHELL (BST) element
    CALL gaun12 (bst_iset)

  CASE (19)     ! skip Gauss data for SPHERE elements
    !CALL gaun16 (spher_iset)

  END SELECT
END DO
IF( drawb_sets > 0 )THEN
  READ(43)  !ttime
  DO iset=1,drawb_sets
    READ(43)  !daz
  END DO
END IF
RETURN
END SUBROUTINE gaussn
