 SUBROUTINE inpd29 (task, nel, iwrit, elsnam, nelms)

 !   READ control DATA for element number 29: Bezier PRISM solid-shell element

 USE ele29_db
 USE gvar_db,ONLY : fimpo,overw
 USE npo_db, ONLY : cpx,coord

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets
                     nel,     & ! number of elements in the set
                     iwrit      ! flag to echo data input

 ! local variables
 LOGICAL :: oldset
 INTEGER (kind=4) :: nreqs, narch, ngaus, nelem, ansmm, anssh, nassp, i

 TYPE (ele29_set), POINTER, SAVE  :: elset, anter

 IF (TRIM(task) == 'INPUT') THEN

   IF( nelms == 0 )THEN  ! for first set only
     ! use Mid-Side in-plane Gauss Points (only for first set)
     midsidegausspoints = exists('MSGP   ') .OR. .NOT.exists('GPINT ')
     IF ( midsidegausspoints )THEN
        WRITE(lures,"(' Gauss points at mid-side points used')")
        psg  = psg_ms
        pa2  = pa2_ms
        pa2b = pa2b_ms
        pag  = pag_ms
     ELSE
        WRITE(lures,"(' Gauss points at inner points used')")
        psg  = psg_in
        pa2  = pa2_in
        pa2b = pa2b_in
        pag  = pag_in
     END IF
   END IF
   ! check if list of sets and set exists and initializes
   CALL srch_ele29 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm29 (1, nelem,  nreqs, narch, elsnam, elset, ngaus, ansmm, anssh, nassp)
     nel = nel - nelem
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     CALL listen('INPD29')  !read a line
     ngaus = getint('NGAUS ',2,' NUMBER OF GAUSS POINTS TTT (>= 2) .')
     IF( ngaus < 2 )THEN
       WRITE(lures,"(' Minimumn number of gauss points is 2 ........')",ERR=9999)
       elset%ngaus = 2
     END IF

     elset%locax = 3       !default
     IF( exists('LOCAL  ',i)) THEN   !local axis definition
       IF( INT(param(i)) > 0 .AND. INT(param(i)) < 4 )THEN  !if a valid code
         elset%locax =  INT(param(i))
       ELSE  !keep default value and p
         WRITE(lures,"(/,5X,'Error in the definition of local system', / &
                         5X,'Invalid axis: ',i3,' Default value used (3)')") INT(param(i))
       END IF
     END IF
     elset%angdf = getrea('ALPHA ',0d0,' First Euler Angle from X and Ortho')
     ansmm =  2 !Use mid-side points of each subtriangle
     IF( exists('ANSMM ',i)) THEN   !Assumed Natural Strain model for Mem
       IF( INT(param(i)) == 0 .OR. INT(param(i)) == 2 )THEN  !if a valid code
         ansmm =  INT(param(i))
         WRITE(lures,"(/,5X,'ANS Membrane Model',i3)") ansmm
       ELSE  !keep default value
         WRITE(lures,"(/,5X,'Error in the option for ANS membrane model', / &
                         5X,'Invalid model: ',i3,' Default value used (2)')") INT(param(i))
       END IF
     END IF
     anssh =  2 !Use incomplete quadratic
     nassp = nnas2
     IF( exists('ANSSH ',i)) THEN   !Assumed Natural Strain model for shear
       IF( INT(param(i)) >= 0 .AND. INT(param(i)) < 3 )THEN  !if a valid code
         anssh =  INT(param(i))
         WRITE(lures,"(/,5X,'ANS Shear Model',i3)") anssh
         IF( anssh /= 2 )THEN
           nassp = nnas1
           IF ( midsidegausspoints )THEN
             pag(:,1:nnas1,:)  = pag0_ms    !(2,nnas1,ngaup)
           ELSE
             pag(:,1:nnas1,:)  = pag0_in    !(2,nnas1,ngaup)
           END IF
         END IF
       ELSE  !keep default value
         WRITE(lures,"(/,5X,'Error in the option for ANS Shear model', / &
                         5X,'Invalid model: ',i3,' Default value used (2)')") INT(param(i))
       END IF
     END IF
     nreqs=getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     elset%small = exists('SMALL')
     IF(elset%small) WRITE(lures,"(' Green strains will be used if possible')")
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements

     !Initialize empty list Point both pointer to nothing
     CALL ini_ele29e (elset%head, elset%tail)
   END IF
   !  read new data or add data to previous data
   CALL elmd29(nelem, elset%head, elset%tail, iwrit, ngaus, ansmm, anssh, nassp )
   ! compute coordinates of control points
   CALL nodx29(nelem,elset%head,coord,cpx)
   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs (ngaus,nreqs, elset%ngrqs, iwrit )

   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele29 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF
   nel = nel + nelem

! ELSE IF (TRIM(task) == 'RESTAR') THEN
!
!   ALLOCATE (elset)          !initializes a list
!   NULLIFY(elset%head)       !nullify head pointer
!   ! read control parameters
!   READ (51) nelem, nreqs, narch, ngaus, &
!             elset%gauss, elset%small, elset%plstr, elset%angdf, &
!             elset%locax
!   ! restore list of elements
!   CALL rest29 (nelem, nreqs, elset%head, elset%tail, &
!                elset%ngrqs, ngaus )
!   ! add to list of elements
!   CALL add_ele29 (elset, head, tail)
!
! ELSE IF (TRIM(task) == 'IMPORT') THEN
!   ! check if list of sets and set exists and initializes
!   IF( overw )THEN
!     CALL srch_ele29 (head, anter, elset, elsnam, oldset)
!     IF (oldset) THEN    !if set exists
!       CALL comm29(1, nelem,  nreqs, narch, elsnam, elset, ngaus)
!     ELSE
!       CALL runen2(' Old set to overwrite does not exist')
!     END IF
!     READ (fimpo)
!   ELSE
!     ALLOCATE (elset)       !reserve memory for set
!     CALL ini_ele29e (elset%head, elset%tail)
!     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
!     READ (fimpo) nelem,ngaus,elset%angdf,elset%small, &
!                  elset%locax
!     nreqs = 0
!     narch = 0
!   END IF
!   ! restore list of elements
!   CALL impo29 ( nelem, ngaus, elset%head, elset%tail)
!   ! add to list of elements
!   IF( .NOT.overw )THEN
!     CALL add_ele29 (elset, head, tail)
!     nelms = nelms + 1 ! increased set counter for this element type
!   END IF
!
 ELSE
   CALL runend('INPD29: NON-EXISTENT TASK .        ')
 END IF
 CALL comm29(0, nelem, nreqs, narch, elsnam, elset, ngaus, ansmm, anssh, nassp)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpd29
