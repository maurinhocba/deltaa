MODULE dyna_db
  !USE ctrl_db, ONLY :   ndyna    ! analysis type
  !                      = 0      !static
  !                      = 1      !damping + consistent mass
  !                      = 2      !damping + lumped mass
  !                      = 3      !undamped + consistent mass
  !                      = 4      !undamped + lumped mass
  ! Information association to dynamic analysis
  IMPLICIT NONE

  REAL (kind=8) :: &
    alpha, &  ! alpha coefficient for Newmark integration scheme
    beta,  &  ! beta coefficient for Newmark integration scheme
    gamma, &  ! gamma coefficient for Newmark integration scheme
    ccm,   &  ! mass coefficient for proportional damping
    cck       ! stiffness coefficient for proportional damping
  LOGICAL :: eigen                 !.TRUE. compute first neigen eigenvalues and eigenvectors


CONTAINS

  SUBROUTINE dump_dyna(neq,maxa)
    IMPLICIT NONE
    INTEGER (kind=4), INTENT(IN) :: neq,maxa
    INTEGER (kind=4) :: n1,m1,i

    WRITE(50,ERR=9999) alpha, beta, gamma, ccm, cck
    !IF( ndyna > 0)THEN
    !  WRITE(50,ERR=9999) (accel(i),i=1,neq),(veloc(i),i=1,neq)
    !  SELECT CASE (ndyna)
    !  CASE (1)
    !    n1 = maxa
    !    m1 = maxa
    !  CASE (2)
    !    n1 = neq
    !    m1 = maxa
    !  CASE (3)
    !    n1 = maxa
    !    m1 = 1
    !  CASE (4)
    !    n1 = neq
    !    m1 = 1
    !  END SELECT
    !
    !  WRITE(50,ERR=9999) (mass(i),i=1,n1),(damp(i),i=1,m1)
    !END IF
    RETURN
    9999 CALL runen2(' ')
  END SUBROUTINE dump_dyna

  SUBROUTINE rest_dyna(neq,maxa)
    IMPLICIT NONE
    INTEGER (kind=4), INTENT(IN) :: neq,maxa

    INTEGER (kind=4) :: i

    READ(51) alpha, beta, gamma, ccm, cck
    !IF( ndyna > 0)THEN
    !  ALLOCATE( accel(neq), veloc(neq) )
    !  READ (51) (accel(i),i=1,neq),(veloc(i),i=1,neq)
    !  SELECT CASE (ndyna)
    !  CASE (1)
    !    ALLOCATE( mass(maxa), damp(maxa) )
    !    READ (51) (mass(i),i=1,maxa),(damp(i),i=1,maxa)
    !  CASE (2)
    !    ALLOCATE( mass(neq), damp(maxa) )
    !    READ (51) (mass(i),i=1,neq),(damp(i),i=1,maxa)
    !  CASE (3)
    !    ALLOCATE( mass(maxa), damp(1) )
    !    READ (51) (mass(i),i=1,maxa),(damp(i),i=1,1)
    !  CASE (4)
    !    ALLOCATE( mass(neq), damp(1) )
    !    READ (51) (mass(i),i=1,neq),(damp(i),i=1,1)
    !  END SELECT
    !
    !ELSE
    !  ALLOCATE( accel(1), veloc(1), damp(1), mass(1) )
    !
    !END IF

  END SUBROUTINE rest_dyna

END MODULE dyna_db
