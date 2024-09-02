SUBROUTINE csinp4(maxve,iwrit,label,npoin,coord,oldsr,oldpr, &
                  codes,surname,tsurf,surdi)

!.... input and generate contact surfaces data

!USE cont4_db
IMPLICIT NONE
!     arguments
INTEGER (kind=4), INTENT(IN) :: iwrit,npoin,label(:),oldpr,maxve,surdi
INTEGER (kind=4), INTENT(IN OUT) :: oldsr,tsurf,codes(2,surdi)
REAL (kind=8), INTENT(IN) :: coord(:,:)
CHARACTER (len=6), INTENT (IN OUT) :: surname(surdi)
END SUBROUTINE csinp4
