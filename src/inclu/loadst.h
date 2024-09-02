 SUBROUTINE loadst (ttime,nsymm,force,stiff,ustif,coora)

 ! calculates and applies stiffness due to follower loads

 !USE load_db, ONLY : nload,loa_set,headls
 !USE hydf_db
 IMPLICIT NONE
 ! dummy arguments
   REAL (kind=8), INTENT(IN) :: ttime,     & !to compute load factor
                                force(:),  & !assembled forces
                                stiff(:),  & !stiffness matrix (symmetric)
                                ustif(:),  & !stiffness matrix (antisymmetric)
                                coora(:,:)   !present coordinates
   INTEGER(kind=4):: nsymm    !0: symmetric matrix only, 1:include unsymmetric part


 END SUBROUTINE loadst
