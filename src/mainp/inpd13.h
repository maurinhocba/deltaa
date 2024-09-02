      SUBROUTINE inpd13 (task, nelem, iwrit, elsnam, nelms)

      !   READ control DATA for element number 13 (TL BST++)

      USE lispa0
      USE ele13_db

      IMPLICIT NONE
      ! dummy arguments
      CHARACTER(len=*) :: elsnam   ! element set name
      CHARACTER(len=* ) :: task     ! requested task
      INTEGER (kind=4) :: nelms,  & ! number of element sets of this type
                          nelem,   & ! number of elements in the set
                          iwrit     ! flag to echo data input

      END SUBROUTINE inpd13
