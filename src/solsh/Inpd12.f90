 SUBROUTINE inpd12 (task, nel, iwrit, elsnam, nelms)

 !   READ control DATA for element number 12: SOLSH laminated solid-shell element

 USE ele12_db
 USE gvar_db,ONLY : fimpo,overw
 USE ctrl_db,ONLY : ndofn
 USE npo_db,ONLY : eule0,euler,coord

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets
                     nel,     & ! number of elements in the set
                     iwrit      ! flag to echo data input
 INTEGER (kind=4), ALLOCATABLE  :: lnods(:,:)
 TYPE (ele12), POINTER :: el
 INTERFACE
   INCLUDE 'nodnor.h'
 END INTERFACE

 ! local variables
 LOGICAL :: oldset,quad,zigzag,cmpea
 INTEGER (kind=4) :: nreqs, narch, nnode, nelem, nnb, i

 TYPE (ele12_set), POINTER, SAVE  :: elset, anter

 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele12 (head, anter, elset, elsnam, oldset)
   cmpea = .FALSE.
   IF (oldset) THEN    !if set exists
     CALL comm12 (1, nnode, nelem,  nreqs, narch, elsnam, elset, zigzag, quad, nnb)
     nel = nel - nelem
   ELSE                !set ELSNAM does not exist
     CALL new_ele12(elset)       !reserve memory for set
     CALL listen('INPD12')  !read a line
     nnb=getint('NNODE ',6,' NUMBER OF BASIC NODES PER ELEMENT.')
     IF( nnb /=6 .AND. nnb/=8 ) CALL runend('INPD12: NNODE must be 6 or 8       ')
     IF(exists('QUAD ').AND.nnb==6) THEN
       WRITE(lures,"(' quaratic approach for in-plane components ...')",ERR=9999)
       quad = .TRUE.
       nnode = 12
     ELSE
       quad = .FALSE.
       nnode = nnb
       elset%lface = .TRUE.
     END IF

     IF( exists('LOCAL  ',i)) THEN   !local axis definition
       IF( elset%locax > 0 .AND. elset%locax < 4 )THEN  !if a valid code
         elset%locax =  INT(param(i))
       ELSE  !keep default value and p
         WRITE(lures,"(/,5X,'Error in the definition of local system', / &
                         5X,'Invalid axis: ',i3,' Default value used (3)')") INT(param(i))
       END IF
     END IF
     elset%check= exists('CHECK ')
     nreqs=getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     elset%angdf = getrea('ALPHA ',0d0,' First Euler Angle from X and Ortho')
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements

     zigzag = .FALSE.
     IF( exists('ZIGZAG') ) THEN !local axis definition
       IF( ndofn == 8 )THEN
         WRITE(lures,"(/,5X,'Additional in-plane DOFs considered for this Set',/)")
         zigzag = .TRUE.
         elset%nstra  = 16        !modifies default value of 8
         IF( exists('EULER') ) cmpea = .TRUE.
       END IF
     END IF
     IF( nnb == 8 )THEN
       IF( exists('BETA  ',i)) THEN   !stabilization factors
         elset%beta = param(i:i+1)
         WRITE(lures,"(' Stabilization factors :',2f8.5)")elset%beta
       ELSE
         elset%beta = (/ 0.15d0, 0.15d0 /)
       END IF
     END IF
   END IF
   !  read new data or data to previous data
   CALL elmd12(nelem, elset%head, elset%tail, iwrit, nnode, elset%nstra, nnb)
   IF (.NOT.oldset) CALL rdreqs (ngaus,nreqs, elset%ngrqs, iwrit )
   IF(cmpea)THEN
     ALLOCATE (lnods(3,nelem))
     el => elset%head
     DO i=1,nelem
       lnods(:,i) = el%lnods(1:3)
       el => el%next
     END DO
     CALL nodnor(nelem,3,lnods,coord,eule0,euler)
     DEALLOCATE (lnods)
   END IF

   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele12 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF
   nel = nel + nelem

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)          !initializes a list
   NULLIFY(elset%head)       !nullify head pointer
   ! read control parameters
   elset%sname = elsnam
   READ (51) elset%nelem, elset%nnode, elset%nreqs, elset%narch, &
             elset%gauss, elset%lface, elset%angdf,              &
             elset%quad,  elset%locax, elset%zigzag
   ! restore list of elements
   CALL rest12 (elset%nelem, elset%nnb, elset%nreqs, elset%head, elset%tail, &
                elset%ngrqs, elset%quad  )
   ! add to list of elements
   CALL add_ele12 (elset, head, tail)

! ELSE IF (TRIM(task) == 'IMPORT') THEN
!   ! check if list of sets and set exists and initializes
!   IF( overw )THEN
!     CALL srch_ele12 (head, anter, elset, elsnam, oldset)
!     IF (oldset) THEN    !if set exists
!       CALL comm12(1, nnode, nelem,  nreqs, narch, elsnam, elset, zigzag, quad,nnb)
!     ELSE
!       CALL runen2(' Old set to overwrite does not exist')
!     END IF
!     READ (fimpo)
!   ELSE
!     ALLOCATE (elset)       !reserve memory for set
!     CALL ini_ele12e (elset%head, elset%tail)
!     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
!     elset%lface = .FALSE.  !initializes flag to compute extended connectivities
!     READ (fimpo) nelem,ngaus,elset%angdf,elset%quad,elset%locax
!     nreqs = 0
!     narch = 0
!   END IF
!   ! restore list of elements
!   CALL impo12 ( nelem, nnb,nnode, elset%nvar, elset%head, elset%tail)
!   ! add to list of elements
!   IF( .NOT.overw )THEN
!     CALL add_ele12 (elset, head, tail)
!     nelms = nelms + 1 ! increased set counter for this element type
!   END IF
!
 ELSE
   CALL runend('INPD12: NON-EXISTENT TASK .        ')
 END IF
 CALL comm12(0, nnode, nelem,  nreqs, narch, elsnam, elset, zigzag, quad, nnb)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpd12
