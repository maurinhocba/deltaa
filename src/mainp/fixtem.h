SUBROUTINE fixtem(iwrit,npoin,iftmp,label)

  !***  APPLY Temperature prescribed values

  IMPLICIT NONE
  INTEGER (kind=4), INTENT(IN) :: iwrit,npoin,label(:)
  INTEGER (kind=4), INTENT(IN OUT) :: iftmp(:,:)

END SUBROUTINE fixtem
