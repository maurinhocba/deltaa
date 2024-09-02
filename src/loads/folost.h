 SUBROUTINE folost (ttime,nsymm,force,stiff,ustif,coora,headf,factor,acurv)
 !     applies follower load
 USE loa_db, ONLY: foll_seg
 !USE ctrl_db, ONLY : ndime,ntype
 IMPLICIT NONE

   !--- Dummy variables
   REAL (kind=8), INTENT(IN) :: ttime,     & !to compute load factor
                                force(:),  & !assembled forces
                                stiff(:),  & !stiffness matrix (symmetric)
                                ustif(:),  & !stiffness matrix (antisymmetric)
                                coora(:,:),& !present coordinates
                                factor       !load factor
   INTEGER(kind=4):: nsymm, & !0: symmetric matrix only, 1:include unsymmetric part
                     acurv    !assigned curve
   TYPE(foll_seg),POINTER:: headf            !pointer to first segment
 END SUBROUTINE folost
