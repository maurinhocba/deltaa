 SUBROUTINE inpda2(task,ndime,nelem,iwrit,elsnam,nelms)
 !******************************************************************
 !
 !*** READ control DATA for 2-node 2/3-d truss element
 !
 !******************************************************************

 IMPLICIT NONE

 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: ndime,   & ! problem dimension
&                    iwrit,   & ! flag to echo data input
&                    nelem,   & ! number of elements
&                    nelms      ! number of element sets of this type
 END SUBROUTINE inpda2
