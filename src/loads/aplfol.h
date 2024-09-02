 SUBROUTINE aplfol (ndime,coora,force,ntype,headf)
 !     applies follower load
 USE loa_db, ONLY: foll_seg
 !USE kin2_db, ONLY : ifpre
 IMPLICIT NONE

   !--- Dummy variables
   INTEGER(kind=4):: ndime, & !problem dimension
                     ntype    !problem type for 2D
   REAL(kind=8), INTENT (IN) :: coora(:,:)   !present coordinates
   REAL(kind=8), INTENT (IN OUT) :: force(:) !assembled forces
   TYPE(foll_seg),POINTER:: headf            !pointer to first segment

 END SUBROUTINE aplfol
