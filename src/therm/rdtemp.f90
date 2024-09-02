SUBROUTINE rdtemp (iwrit,actio)

  !***  READ fixed and prescribed temperatures

  USE param_db,ONLY: mnam
  USE ift_db       !temperature prescribed data_base
  USE c_input
  USE nsets_db
  USE curv_db, ONLY : del_cur

  IMPLICIT NONE
  CHARACTER(len=*) :: actio
  INTEGER (kind=4) :: iwrit    !for echo into .RSN

  ! Local
  LOGICAL set, found
  CHARACTER(len=mnam) :: sname,lbl
  INTEGER :: i,j,k,n,ki,kf,icod,nnods,nod,nrve,nrpt
  INTEGER (kind=4), ALLOCATABLE :: nods(:)
  REAL (kind=8) :: v
  TYPE (rpt_set), POINTER :: rves,rpts,sant,spos
  TYPE (rpt_nod), POINTER :: rven
  TYPE (nset), POINTER :: ns

  INTERFACE
     INCLUDE 'rdpret.h'
  END INTERFACE

 !------------------------------------------------------------

  IF (iwrit == 1) WRITE(lures, "(' PRESCRIBED TEMPERATURES IN TIME ',/)",ERR=9999)

  IF (npret > 0) THEN              !if previous prescribed values
    CALL listen('RDTEMP')          !read a card
    IF (exists('DELETE')) THEN     !Key word DELETE is present
      IF (iwrit == 1) WRITE(lures,"(' Deleting prescribed temperature sets ',/)",ERR=9999)
      IF (TRIM(actio) == 'NSTRA0') actio = 'NSTRA1'

      DO                           ! loop over the sets to delete
        CALL listen('RDTEMP')      ! read a card
        IF (exists('ENDDEL') ) EXIT     !kew-word END_DELETE found => exit loop
        lbl = get_name('TEMSET',stype='CURV',        &   !assoc curpt
              texts='!PRESCRIBED TEMPERATURE SET .......')
        CALL srch_rpts (headv, sant, rpts, lbl, found)  !search for the curpt
        IF (.NOT.found) THEN                    !If not found error in data input
          WRITE(lures, "(' Warning! Temperature set using curpt', &
     &              a,' does not exist')",ERR=9999)TRIM(lbl)
        ELSE
          CALL del_rpts (headv, tailv, sant, rpts)   !delete set assoc to LC
          CALL del_cur (lbl)                         !delete curve data associated LC
          npret = npret - 1                          !correct counter
        END IF
      END DO

    ELSE              !nothing to delete
      backs = .TRUE.                       !one line back
    ENDIF

  END IF

  IF (iwrit == 1) WRITE (lures, "(//)",ERR=9999)

  DO
    ! loop over temperature sets - reading
    CALL listen('RDTEMP')              !read a card
    IF (.NOT.exists('TEMSET')) THEN    ! if key-word TEM_SET not found
      backs = .TRUE.                       !one line back
      EXIT                             ! exit loop
    END IF
    ALLOCATE (rpts)                    !get memory for a list of nodes (SET)
    rpts%lc = get_name('TEMSET',stype='CURV',        &   !assoc curpt
              texts='!PRESCRIBED TEMPERATURE SET ...')
    rpts%factor=getrea('FACTOR',1d0, ' Participation factor of this set .')

    !check if associated function LC already used (if TRUE stop)

    CALL srch_rpts (headv, sant, spos, rpts%lc, found)
    IF (found) CALL runend ('RDTEMP: Set using this name already')
    CALL rdcurv('VELOC',rpts%lc)

    IF(iwrit == 1) THEN    !echo title
      WRITE(lures, "(//,5X,' PRESCRIBED TEMPERATURES',//)",ERR=9999)
      IF ( ndoft == 1) WRITE(lures,  &
         "(6X,'NODE',5X,'TEMP.'/)",ERR=9999)
      IF ( ndoft == 2) WRITE(lures,  &
         "(6X,'NODE',5X,'TEMP-I',9X,'TEMP-S.',/)",ERR=9999)
      IF ( ndoft == 3) WRITE(lures,  &
         "(6X,'NODE',5X,'TEMP-N',8X,'TEMP-I',8X,'TEMP-S',/)",ERR=9999)
    END IF
    IF (TRIM(actio) == 'NSTRA0') actio = 'NSTRA1'

    CALL rdpret (iwrit, rpts%head, rpts%tail,nrpt) !read prescribed temp sets

    rpts%nrv=nrpt                       !keep number of nodes in the set
    CALL add_rpts (rpts, headv, tailv)  !add set of nodes to the end of the list
    npret = npret + 1                   !increment number of sets

  END DO  ! loop for sets of prescribed temperatures
  CALL listen('RDTEMP')  ! read a line (END_PRESCRIBED line expected)
  IF (.NOT. exists('ENDPRE'))CALL runend ('RDTEMP:end_prescribed_temper expc.')
  RETURN
  ! the convention: the last read overrides the previous data
  ! prescribed temperature will be applied to the nodes previously
  ! fixed o released


 RETURN
 9999 CALL runen2('')
END SUBROUTINE rdtemp
