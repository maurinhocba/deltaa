 SUBROUTINE inpd04 (task, nel, iwrit, elsnam, nelms)

 !   READ control DATA for element number 04 (TL CST++)

 USE ele04_db

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*) :: elsnam,  &  ! element set name
                     task        ! requested task
 INTEGER (kind=4) :: iwrit,   & ! flag to echo data input
                     nel,     & ! number of elements
                     nelms      ! number of element sets

 ! local variables
 LOGICAL :: oldset
 INTEGER (kind=4) :: nreqs, narch, nelem

 TYPE (ele04_set), POINTER, SAVE  :: elset, anter

 IF (TRIM(task) == 'INPUT ') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele04 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm04 (1, nelem,  nreqs, narch, elsnam, elset)
     elset%lside = .NOT.elset%linear  !initializes flag to compute LSIDE
     nel = nel - nelem
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     CALL listen('INPD04')  !read a line
     nreqs=getint('NREQS ',0,' Gauss pt for stress time history..')
     elset%angdf(1) =getrea('ALPHA ',0d0,' First Euler Angle from X and Ortho')
     elset%angdf(2) =getrea('BETA  ',0d0,' Second Euler Angle from X and Orth')
     elset%angdf(3) =getrea('GAMMA ',0d0,' Third Euler Angle from X and Ortho')
     IF( elset%linear )THEN
       WRITE(lures,"(/,5X,'Volume Approximation based on neighbours will be used')")
       !elset%vfact =getrea('VOLFAC',0.2d0,' Volume factor for neighbours      ')
     ELSE
       WRITE(lures,"(/,5X,'STANDARD Approximation Strain) will be used')")
     END IF
     elset%lside = .NOT.elset%linear  !initializes flag to compute LSIDE

     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements

     !Initialize empty list Point both pointer to nothing
     CALL ini_ele04e (elset%head, elset%tail)
   END IF
   !  read new data or add to previous data
   CALL elmd04( nelem, elset%head, elset%tail, iwrit)
   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs (ngaus,nreqs, elset%ngrqs, iwrit )

   CALL comm04(0, nelem,  nreqs, narch, elsnam, elset)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele04 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF
   nel = nel + nelem

 !ELSE IF (TRIM(task) == 'RESTAR') THEN
 !
 !  ALLOCATE (elset)          !initializes a list
 !  NULLIFY(elset%head)       !nullify head pointer
 !  !  read control parameters
 !  elset%sname = elsnam
 !  READ (51) elset%nelem, elset%nreqs, elset%narch, elset%gauss, &
 !      elset%plstr, elset%angdf, elset%btscal, elset%small, elset%linear, elset%vfact
 !  CALL rest04 (elset%nelem,  elset%nreqs, elset%head, elset%tail, &
 !               elset%ngrqs, elset%linear)
 !  ! add to list of elements
 !  CALL add_ele04 (elset, head, tail)
 !
 ELSE
   CALL runend('INPD04: NON-EXISTENT TASK .        ')
 END IF
 RETURN

 END SUBROUTINE inpd04
