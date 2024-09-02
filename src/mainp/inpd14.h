 SUBROUTINE inpd14 (task, nelem, iwrit, elsnam, nelms)

 !   READ control DATA for element number 14 (TL BST++)

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets of this type
                     nelem,   & ! number of elements in the set
                     iwrit      ! flag to echo data input

 END SUBROUTINE inpd14
