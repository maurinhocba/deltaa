      SUBROUTINE histor( )
!***********************************************************************
!
! ** initial values and output requested DATA
!
!***********************************************************************
      USE param_db,ONLY: mnam
      USE ctrl_db, ONLY: memo,ndyna
      USE c_input
      USE npo_db, ONLY : label
      USE outp_db

      IMPLICIT NONE


      INTEGER (kind=4) :: i,j,n
      INTEGER (kind=4), ALLOCATABLE ::  m(:)

      INTEGER (kind=4) chnode


      ALLOCATE ( m(memo) )     ! auxiliar array for temporary storage

!-----------------     Read requested output for time history

      CALL listen('INTIME')                       !read a card

      IF (.NOT.exists('HISTOR')) THEN             !search for key-word HISTORY

        backs = .TRUE.                          !One line back

      ELSE             ! New data for nodal output will be read

        DO                            !loop for each output variable requested
          IF (exists('ENDHIS')) EXIT  !key-word END_HISTORY found, exit loop

          ! ***  READ selected nodes for output (displacements)

          IF (exists('DISPLA')) THEN              !key-word DISPLACEMENTS
            CALL rdfrin('INTIME',m,n,memo)   !read the list of labels
            nreqd = 0
            DO i=1,n               !change from labels to internal numeration
              j = chnode(m(i))
              IF( j > 0 )THEN
                nreqd = nreqd + 1
                m(nreqd) = j
              END IF
            END DO
            IF(ASSOCIATED(nprqd)) DEALLOCATE( nprqd)  !release previous memory
            IF( nreqd > 0 )THEN
              ALLOCATE( nprqd(nreqd) )              !get memory for the list
              nprqd = m(1:nreqd)                    !store the list
            END IF
            CALL listen('INTIME')                 !read key-word END_DISPLAC
            IF (.NOT.exists('ENDDIS'))          & !Error if not present
     &          CALL runend('INTIME:END_DISPLACEMENT EXPECTED   ')

          END IF

          ! ***  READ selected nodes for output (velocities)

          IF (exists('VELOCI')) THEN             !key-word VELOCITIES
            CALL rdfrin('INTIME',m,n,memo)   !read the list of labels
            nreqv = 0
            IF( ndyna > 0 )THEN
              DO i=1,n               !change from labels to internal numeration
                j = chnode(m(i))
                IF( j > 0 )THEN
                  nreqv = nreqv + 1
                  m(nreqv) = j
                END IF
              END DO
              IF(ASSOCIATED(nprqv)) DEALLOCATE( nprqv)
              IF( nreqv > 0 )THEN
                ALLOCATE( nprqv(nreqv) )              !get memory for the list
                nprqv = m(1:nreqv)                    !store the list
              END IF
            END IF
            CALL listen('INTIME')                 !read key-word END_VELOCITIES
            IF (.NOT.exists('ENDVEL'))          & !Error if not present
     &        CALL runend('INTIME: END_VELOCITIES EXPECTED    ')
          END IF

          ! ***  READ selected nodes for output (accelerations)

          IF (exists('ACCELE')) THEN             !key-word ACCELERATIONS
            CALL rdfrin('INTIME',m,n,memo)   !read the list of labels
            nreqa = 0
            IF( ndyna > 0 )THEN
              DO i=1,n               !change from labels to internal numeration
                j = chnode(m(i))
                IF( j > 0 )THEN
                  nreqa = nreqa + 1
                  m(nreqa) = j
                END IF
              END DO
              IF(ASSOCIATED(nprqa)) DEALLOCATE( nprqa)
              IF( nreqa > 0 )THEN
                ALLOCATE( nprqa(nreqa) )              !get memory for the list
                nprqa = m(1:nreqa)                    !store the list
              END IF
            END IF
            CALL listen('INTIME')                !read key-word END_ACCELERATIONS
            IF (.NOT.exists('ENDACC'))         & !Error if not present
     &        CALL runend('INTIME: END_ACCELERATIONS EXPECTED ')

          END IF

          !     ***  READ selected nodes for output (internal loads)

          IF (exists('INTERN')) THEN             !key-word INTERNAL
            CALL rdfrin('INTIME',m,n,memo)   !read the list of labels
            nreql = 0
            DO i=1,n               !change from labels to internal numeration
              j = chnode(m(i))
              IF( j > 0 )THEN
                nreql = nreql + 1
                m(nreql) = j
              END IF
            END DO
            IF(ASSOCIATED(nprql)) DEALLOCATE( nprql)
            IF( nreql > 0 )THEN
              ALLOCATE( nprql(nreql) )              !get memory for the list
              nprql = m(1:nreql)                    !store the list
            END IF

            CALL listen('INTIME')               !read key-word END_INTERNAL
            IF (.NOT.exists('ENDINT'))        & !Error if not present
     &        CALL runend('INTIME:END_INTERNAL_FORCES EXPECTED')

          END IF

          ! ***  READ selected nodes for output (contact forces)

          IF (exists('CONTAC')) THEN             !key-word CONTACT
            CALL rdfrin('INTIME',m,n,memo)   !read the list of labels
            nreqc = 0
            DO i=1,n               !change from labels to internal numeration
              j = chnode(m(i))
              IF( j > 0 )THEN
                nreqc = nreqc + 1
                m(nreqc) = j
              END IF
            END DO
            IF(ASSOCIATED(nprqc)) DEALLOCATE( nprqc)
            IF( nreqc > 0 )THEN
              ALLOCATE( nprqc(nreqc) )              !get memory for the list
              nprqc = m(1:nreqc)                    !store the list
            END IF

            CALL listen('INTIME')                !read key-word END_CONTACT
            IF (.NOT.exists('ENDCON'))         & !Error if not present
              CALL runend('INTIME: END_CONTACT_FORCES EXPECTED')

          END IF

          CALL listen('INTIME')

        END DO

      END IF

      IF(.NOT.ASSOCIATED(nprqd)) ALLOCATE( nprqd(1) )      !allocate 1 space for
      IF(.NOT.ASSOCIATED(nprqv)) ALLOCATE( nprqv(1) )      !non requested output
      IF(.NOT.ASSOCIATED(nprqa)) ALLOCATE( nprqa(1) )
      IF(.NOT.ASSOCIATED(nprql)) ALLOCATE( nprql(1) )
      IF(.NOT.ASSOCIATED(nprqc)) ALLOCATE( nprqc(1) )

      IF(iwrit == 1) THEN      !ECHO the requested OUTPUT

        IF (nreqd > 0) THEN
          WRITE(lures,"(//,5x,' Selective Output Requested for Displacements ')",ERR=9999)
          WRITE(lures,"(5x,10i5)",ERR=9999) (label(nprqd(i)),i=1,nreqd)
        END IF

        IF (nreqv > 0) THEN
          WRITE(lures,"(//,5x,' Selective Output Requested for Velocities ')",ERR=9999)
          WRITE(lures,"(5x,10i5)",ERR=9999) (label(nprqv(i)),i=1,nreqv)
        END IF

        IF (nreqa > 0) THEN
          WRITE(lures,"(//,5x,' Selective Output Requested for Accelerations ')",ERR=9999)
          WRITE(lures,"(5x,10i5)",ERR=9999) (label(nprqa(i)),i=1,nreqa)
        END IF

        IF (nreql > 0) THEN
          WRITE(lures,"(//,5x,' Selective Output Requested for Internal Loads')",ERR=9999)
          WRITE(lures,"(5x,10i5)",ERR=9999) (label(nprql(i)),i=1,nreql)
        END IF

        IF (nreqc > 0) THEN
          WRITE(lures,"(//,5x,' Selective output requested for Contac Forces ')",ERR=9999)
          WRITE(lures,"(5x,10i5)",ERR=9999) (label(nprqc(i)),i=1,nreqc)
        END IF

      END IF

      DEALLOCATE ( m )     !release auxiliar space

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE histor
