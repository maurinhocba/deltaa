 SUBROUTINE inpd25 (task, nelem, iwrit, elsnam, nelms)
 !   READ control DATA for element number 25 (TL BSQ)

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*) :: elsnam    ! element set name
 CHARACTER(len=*) :: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets of this type
                     nelem,   & ! number of elements in the set
                     iwrit      ! flag to echo data input

 END SUBROUTINE inpd25
