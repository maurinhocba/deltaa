SUBROUTINE conini(actio, maxve)

  !*** contact data READ

  USE c_input
  USE ctrl_db, ONLY: ctype, mconts, ndime, npoin, nconp, numct
  USE npo_db,  ONLY: fcont
  USE outp_db, ONLY: nreqc, nprqc, iwrit
  IMPLICIT NONE

  CHARACTER(len=*),INTENT(IN):: actio
  INTEGER(kind=4),INTENT(INOUT):: maxve

  INTEGER (kind=4) i, ncont

  !INTERFACE
  !  INCLUDE 'coninp.h'
  !END INTERFACE

  CALL listen('CONINI')
  IF (.NOT.exists('CONTAC')) THEN
    IF ( numct == 0 ) THEN
      ALLOCATE( fcont(1,1) )
      IF (nreqc > 0) THEN
        DEALLOCATE (nprqc)
        ALLOCATE (nprqc(1))
      END IF
    END IF
    backs = .TRUE.

  ELSE

    WRITE (lures,"(/,'  C O N T A C T   P A R A M E T E R S ',/)",ERR=9999)
    DO
      CALL listen('CONINI')
      IF (exists('ENDCON')) EXIT
      ncont = getint('NCONT ',0,'!Contact Type  ....................')
      IF (ncont < 0 .OR. ncont > mconts) CALL runend('INPDAT:NCONT MUST BE 1-7')
      !CALL coninp (actio, ncont, iwrit=iwrit, maxve=maxve)
    END DO
    IF (TRIM(actio) /= 'NEW') DEALLOCATE ( fcont )
    ALLOCATE( fcont(ndime,npoin) )
    fcont = 0d0
  END IF

  numct = 0     !initializes
  ctype = 0

  DO i=1,mconts                 !for each contact algorithm
    IF (nconp(i) > 0) THEN      !if there are contact pairs using this alg.
      numct = numct + 1         !increase the number of cont algors used
      ctype (numct) = i         !keep for this set the element type
    END IF
  END DO

  IF (TRIM(actio) == 'NSTRA2' .AND. numct > 0) THEN
    DEALLOCATE ( fcont )
    ALLOCATE( fcont(ndime,npoin) )
  END IF

  RETURN
 9999 CALL runen2('')
END SUBROUTINE conini
