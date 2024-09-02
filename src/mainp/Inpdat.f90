 SUBROUTINE inpdat (actio)
 !-----------------------------------------------------------------------
 !***  input data routine
 !-----------------------------------------------------------------------
  USE param_db,ONLY: midn, mnam, mich
  USE ctrl_db, ONLY: ndime, ndofn, nrotd, neulr, npoin, ndoft, itemp
  USE esets_db,ONLY : nelms, nel, elsets, elty, msets, add_name, rot_free
  USE ifx_db, ONLY : nifx
  USE kinc_db
  USE c_input
  USE nes_db,  ONLY : nnes,reasd0        !slave nodes
  USE npo_db, ONLY: tempe,dtemp
  USE ndp_db, ONLY : nndp
  USE mat_dba, ONLY : mat_inp
  USE outp_db, ONLY : iwrit
  IMPLICIT NONE

   !--- Dummy variable
   CHARACTER(len=*),INTENT(IN OUT):: actio !NEW, NSTRA0
   !--- Functions
   CHARACTER(len=mich):: inttoch
   !--- Local variable
   CHARACTER(len=midn) :: elmty
   CHARACTER(len=mnam) :: elsnam
   INTEGER (kind=4) :: i,j,n
   LOGICAL          :: found

  INTERFACE
    INCLUDE 'arisi0.h'
    INCLUDE 'elmdat.h'
    INCLUDE 'fixva0.h'
    INCLUDE 'rdvel0.h'
    INCLUDE 'readp0.h'
  END INTERFACE


  CALL import( actio)  ! import sets from finary files

  CALL rdcoor (actio)  ! read coordinates (& nodal systems)

  ! input of material properties
  CALL mat_inp ( )

   !******* Input of Element data SETs
   !        if added element sets
   !         actio = NSTRA0 => actio = NSTRA1
   CALL listen('INPDAT')

   IF(exists('SETDEF'))THEN

     WRITE (lures,"(/,80('-'),//,' ELEMENT SETS DEFINITION',/)",ERR=9999)

     CALL listen('INPDAT')

     IF (TRIM(actio) == 'NEW') THEN
       nelms = 0
       nel = 0
     ELSE
       IF (TRIM(actio) == 'NSTRA0') actio='NSTRA2'
     END IF

     DO
       IF (exists('ENDSET') )EXIT
       IF (exists('ELMTY ',i) ) THEN
         elmty = words(i+1)(1:midn)
         i = 1
         DO
           IF( TRIM(elmty) == TRIM(elty(i)) )EXIT    !found element type
           i = i+1
           IF( i <= msets)CYCLE
           WRITE(lures,*,ERR=9999)' Unrecognizable element type ',TRIM(elmty)
           CALL runend ('INPDAT:A NON-VALID ELEMNT TYPE READ')
         END DO
         elsnam = get_name('ELSNAM',found,stype='ESET')
         IF (found) THEN
           WRITE(lures, "(' Element set named: ',A,'  will be read')",ERR=9999) TRIM(elsnam)
         ELSE
           elsnam = 'ES'//TRIM(inttoch(i,2))//TRIM(inttoch(nelms(i),2))
           WRITE(lures, "(' Warning:',/  &
             & ' An element set without name has been defined.', &
             & ' Name ',A,' will be used for future reference')",ERR=9999) TRIM(elsnam)
         END IF
         CALL add_name(elsnam,2)
         WRITE(lures,"(/'   Elements type (',A,')'/)",ERR=9999) TRIM(elty(i))
         IF((i > 5 .AND. i <10) .AND. (neulr == 0))THEN
           WRITE(lures, "(' *ERROR, set NEULR=1 if beam/shell elements')",ERR=9999)
           !WRITE(*    , "(' *ERROR, set NEULR=1 if beam/shell elements')")
           STOP
         END IF

         CALL elmdat('INPUT ',elsnam,i)

       ELSE
         CALL runend('INPDAT: ELMT_TYPE CARD EXPECTED    ')
       END IF

       CALL listen('INPDAT')
     END DO
     CALL elsets ( )  !computes ESETS & ESET(:)

   ELSE IF (TRIM(actio) == 'NEW') THEN
     CALL runend('INPDAT: SET_DEFINITION EXPECTED    ')

   ELSE
     backs = .TRUE.                       !one line back
   END IF

!***  READ the kinematic restrictions

  CALL listen('INPDAT')              ! read a lines
  IF (.NOT.exists('KINEMA') ) THEN   ! check key-word
    backs = .TRUE.                       !one line back
  ELSE
    WRITE(lures, "(//' KINEMATIC CONDITIONS read in this strategy',/)",ERR=9999)
    DO                               !loop for constraints
      CALL listen('INPDAT')          !read a line

      SELECT CASE (words(1)(1:midn))         !according to first key-word

      CASE ('NESCVN')                ! read multi-point constraint data
        CALL reasd0 (iwrit,actio)

      CASE ('NDEPDN')                ! read dependent nodes
        CALL readp0 (iwrit,actio)

       CASE ('SLIDIN')                ! read sliding nodes
         !IF( order > 3 )CALL runend('INPDAT: Kinematic constrains not in order, check please     ')
         CALL rdsldn (iwrit,actio)

      ! CASE ('RFNDEP')                ! read Rotation Free Nodal DEPendencies
      !   CALL rdrfd0 (iwrit,actio)

      CASE ('NARISN')                ! read dependent nodes on a side
        CALL arisi0 (naris, iwrit,actio)

      CASE ('BOUNDA')                ! read fixed and released node data
        CALL fixva0 (iwrit,ndime,nrotd,actio,rot_free,nvelr)

      CASE ('PRESCR')                ! read prescribed velocities
        CALL rdvel0 (iwrit,ndime,ndofn,actio,nvelr)

      CASE ('ENDKIN')                ! key-word for end of data
        EXIT
      CASE DEFAULT
        CALL runend('INPDAT: UNEXPECTED DATA READ.      ')
      END SELECT
    END DO
    ndepd = nndp + ndumn
    IF( neulr == 0 .AND. (ndepd+naris) > 0) THEN  !Neulr > 0 for dependant nodes
      WRITE(lures, "(' WARNING, set NEULR=1 if ndepd >0 or naris >0 ')",ERR=9999)
      WRITE(*    , "(' WARNING, set NEULR=1 if ndepd >0 or naris >0 ')")
    END IF
    IF(actio == 'NSTRA0')actio = 'NSTRA1'       !if only kinematic conditions modified
  END IF ! IF (.NOT.exists('KINEMA') )

 !    READ Heat transfer essential boundary conditions
  CALL listen('INPDAT')              ! read a lines
  IF (.NOT.exists('TEMPER') ) THEN   ! check key-word
    backs = .TRUE.                       !one line back
  ELSE
    WRITE(lures, "(//' TEMPERATURE CONDITIONS read in this strategy',/)",ERR=9999)

    CALL listen('INPDAT')          !read a line

    IF( exists('PRESCR') )THEN          ! read prescribed temperature in time
      CALL rdtemp (iwrit,actio)
    ELSE
      CALL runend('INPDAT: UNEXPECTED DATA READ.      ')
    END IF

    itemp = .TRUE.
    IF( .NOT.ASSOCIATED (tempe) )THEN
      ALLOCATE (tempe(ndoft,npoin), dtemp(ndoft,npoin))
      tempe = 0d0
      dtemp = 0d0
    END IF
    CALL listen('INPDAT')              ! read a lines
    IF (.NOT.exists('ENDTEM') ) CALL runend('INPDAT: UNEXPECTED DATA READ.      ')  ! check key-word
  END IF ! IF (.NOT.exists('TEMPER') )

  CALL rdconm (actio,nrotd,iwrit)      !read concentrated masses

  CALL flushf(lures)                     !flush output file

  RETURN
  9999 CALL runen2('')
END SUBROUTINE inpdat
