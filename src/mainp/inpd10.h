 SUBROUTINE inpd10 (task, ndime, nel, iwrit,elsnam,nelms)
 !*****************************************************************
 !
 !*** READ control DATA for RIGID elements
 !
 !*****************************************************************

 !USE ctrl_db, ONLY: ndime,ndofn ! problem dimension  & number of DOFs per node
 !USE ele10_db
 IMPLICIT NONE
 !     routine parameters

 CHARACTER(len=*),INTENT(IN) :: elsnam
 CHARACTER(len=*),INTENT(IN) :: task
 INTEGER (kind=4) ndime,        & ! problem dimensio
                  nel,          & ! number of elements in the set
                  nelms,        & ! number of element sets
                  iwrit           ! flag to echo data input

 END SUBROUTINE inpd10
