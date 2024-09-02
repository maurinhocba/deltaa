 SUBROUTINE inpd30 (task, nel , iwrit, elsnam, nelms)

 !   READ control DATA for element number 30 (TL BST++)

 USE ele30_db
 USE gvar_db

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets of this type
                     nel,     & ! number of elements in the set
                     iwrit      ! flag to echo data input

 ! local variables
 LOGICAL :: oldset, logst
 INTEGER (kind=4) :: nreqs, narch, nelem, nnode, i
 CHARACTER(len=mnam) :: sname     ! element set name

 TYPE (ele30_set), POINTER, SAVE  :: elset, anter

 sname = elsnam
 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele30 (head, anter, elset, sname, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm30 (1, nelem,  nreqs, narch, sname, elset, logst)
     elset%lside = .FALSE.  !initializes flag to compute LSIDE
     nel = nel - nelem
   ELSE                !set ELSNAM does not exist
     CALL new_ele30(elset)
     CALL listen('INPD30')  !read a line
     nreqs = getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     nnode = getint('NNODE ',3,' Nodes included in connectivities..')
     elset%angdf = getrea('ANGLE ',0d0,' Default angle between X1 and Ort_1')
     IF( exists('LOCAL  ',i)) THEN   !local axis definition
       elset%locax =  INT(param(i))
       IF( elset%locax < 1 .OR. elset%locax > 3 )THEN
         WRITE(lures,"(/,5X,'Error in the definition of local system', / &
                         5X,'Invalid axis: ',i3,' Default value used (3)')") elset%locax
         elset%locax = 3
       END IF
     END IF
     logst = .NOT.exists('SMALL ')
     narch = 0           !to check
     nelem = 0           !new set, initializes number of elements

     IF(iwrit == 1)WRITE(lures,"(/,5X,'CONTROL PARAMETERS FOR SHELL ELEMENT' &
                   &       //,5X,'REQUIRED STRESS (NREQS) =',I10,  &
                   &       //,5X,'NODES in CONNS  (NNODE) =',I10,  &
                   &       //,5X,'USE G-L strains (SMALL) =',L10,/)",ERR=9999)&
                   nreqs,nnode,.NOT.logst

   END IF
   !  read new data or add to previous data
   CALL elmd30(nelem, elset%head, elset%tail, iwrit, nnode)
   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs ( 1 ,nreqs, elset%ngrqs, iwrit )
   CALL comm30(0, nelem,  nreqs, narch, sname, elset, logst)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele30 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF
   IF( ASSOCIATED(elset%stint) )DEALLOCATE(elset%stint)
   ALLOCATE(elset%stint(13,nelem))
   elset%stint = 0d0
   nel = nel + nelem

! ELSE IF (TRIM(task) == 'RESTAR') THEN
!
!   CALL new_ele30(elset)
!   ! read control parameters
!   elset%sname = sname
!   READ (51) elset%nelem, elset%nreqs, elset%narch, elset%logst, elset%lside, &
!             elset%gauss, elset%plstr, elset%angdf, elset%nonrg, elset%locax
!   ! restore list of elements
!   ALLOCATE (elset%stint(13,elset%nelem))         !initializes a list
!   CALL rest30 (elset%nelem, elset%nreqs, elset%head, elset%tail, &
!                elset%ngrqs, elset%stint)
!   ! add to list of elements
!   CALL add_ele30 (elset, head, tail)
!
! ELSE IF (TRIM(task) == 'IMPORT') THEN
!   ! check if list of sets and set exists and initializes
!   IF( overw )THEN
!     CALL srch_ele30 (head, anter, elset, elsnam, oldset)
!     IF (oldset) THEN    !if set exists
!       CALL comm30 (1, nelem,  nreqs, narch, sname, elset, logst)
!       !elset%lside = .FALSE.   !initializes flag to compute LSIDE
!     ELSE
!       CALL runen2(' Old set to overwrite does not exist')
!     END IF
!     READ (fimpo)
!   ELSE
!     CALL new_ele30(elset)
!     elset%sname = sname
!     READ (fimpo) nelem,elset%angdf,logst,elset%nonrg,elset%locax
!     nreqs = 0
!     narch = 0
!     ALLOCATE(elset%stint(13,nelem))
!     elset%stint = 0d0
!   END IF
!   ! restore list of elements
!   CALL impo30 ( nelem,elset%head, elset%tail)
!   CALL comm30 (0, nelem,  nreqs, narch, sname, elset, logst)
!   ! add to list of elements
!   IF( .NOT.overw )THEN
!     CALL add_ele30 (elset, head, tail)
!     nelms = nelms + 1 ! increased set counter for this element type
!     nel = nel + nelem
!   END IF
!
 ELSE
   CALL runend('INPD30: NON-EXISTENT TASK .        ')
 END IF

 RETURN
 9999 CALL runen2('')

 END SUBROUTINE inpd30
