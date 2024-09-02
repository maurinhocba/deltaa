SUBROUTINE gidresbin(etype )
!  Print result file for GiD binary
USE data_db,ONLY: spot_sets, truss_nps, sol2d_nps, sol3d_nps, shl3d_nps, beame_nps,   &
                  shrev_nps, bst_nps, drawb_sets
IMPLICIT NONE
  INTEGER (kind=4), INTENT(IN) :: etype

  IF (spot_sets > 0) CALL prgib01( )
  IF (truss_nps > 0) CALL prgib02( )
  IF (sol2d_nps > 0) CALL prgib03( )
  IF (sol3d_nps > 0) CALL prgib05( )
  IF (shl3d_nps > 0) CALL prgib06( )
  IF (beame_nps > 0) CALL prgib08( )
  IF (shrev_nps > 0) CALL prgib09(etype)
  IF (bst_nps   > 0) CALL prgib12( )
  IF (drawb_sets > 0)CALL prgidb ( )

RETURN
END SUBROUTINE gidresbin
