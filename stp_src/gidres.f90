 SUBROUTINE gidres (etype )
 !
 !  Print result file for GiD
 !
 USE data_db
 IMPLICIT NONE
  INTEGER (kind=4), INTENT(IN) :: etype

  IF(spot_sets > 0)CALL prgi01 ( )
  IF(truss_nps > 0)CALL prgi02 ( )
  IF(sol2d_nps > 0)CALL prgi03 ( )
  IF(sol3d_nps > 0)CALL prgi05 ( )
  IF(shl3d_nps > 0)CALL prgi06 ( )
  IF(beame_nps > 0)CALL prgi08 ( )
  IF(shrev_nps > 0)THEN
    IF( ntype /= 4 )THEN
      CALL prgi09 (etype)
    ELSE
      CALL prgi09c( )
    END IF
  END IF
  IF(bst_nps    > 0)CALL prgi12 ( )
  IF(drawb_sets > 0)CALL prgidb ( )

 RETURN
 END SUBROUTINE gidres
