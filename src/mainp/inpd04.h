 SUBROUTINE inpd04 (task, nelem, iwrit, elsnam, nelms)

 !   READ control DATA for element number 04 (TL CST++)

 !USE ele04_db

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*) :: elsnam,  &  ! element set name
                     task        ! requested task
 INTEGER (kind=4) :: iwrit,   & ! flag to echo data input
                     nelem,   & ! number of elements
                     nelms      ! number of element sets

 END SUBROUTINE inpd04
