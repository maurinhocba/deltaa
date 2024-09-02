 SUBROUTINE inpd18 (task, nelem, iwrit, elsnam, nelms)

 !   READ control DATA for element number 18 (TL QUAD 4)

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets
                     nelem,   & ! number of elements in the set
                     iwrit      ! flag to echo data input
 END SUBROUTINE inpd18
