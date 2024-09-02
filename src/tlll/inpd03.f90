 SUBROUTINE inpd03 (task,nel, eule0,euler,coord, iwrit,elsnam,nelms)
 !******************************************************************
 !
 !*** READ control DATA for 6-node triangular shell element
 !
 !******************************************************************

 USE ele03_db
 USE ctrl_db, ONLY : ndofn
 USE gvar_db

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam
 CHARACTER(len=*),INTENT(IN):: task
 INTEGER (kind=4) :: nel,nelms,iwrit
 REAL    (kind=8) :: coord(:,:),eule0(:,:),euler(:,:)


 ! local variables
 LOGICAL :: cmpea, &   !CoMPute Euler Angles
            oldset     !
 INTEGER (kind=4) :: nreqs, narch,i,ng,stype,nelem
 INTEGER (kind=4), ALLOCATABLE  :: lnods(:,:)
 TYPE (ele03_set), POINTER :: elset,anter
 TYPE (ele03), POINTER :: el
 CHARACTER(len=mnam) :: sname     ! element set name

 INTERFACE
   INCLUDE 'nodn03.h'
 END INTERFACE

  sname = elsnam
 ! ALLOCATE (elset)

 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele03 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm03 (1, nelem,  nreqs, narch, elsnam, elset)
     cmpea = .FALSE.
     nel = nel - nelem
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     elset%lside = .FALSE.  !     "     flag to compute LSIDE
     CALL listen('INPDA7')
     nreqs = getint('NREQS ',0,'!GAUSS PT FOR STRESS TIME HISTORY .')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     cmpea = exists('EULER ')
     IF(cmpea)WRITE(lures,"(/,5X,' Automatic Evaluation of Local ', &
                          &      'Nodal Systems Selected'//)",ERR=9999)
     elset%angdf=getrea('ANGLE ',0d0,' Default angle between X1 and Ort_1')
     IF( exists('LOCAL  ',i)) THEN   !local axis definition
       elset%locax =  INT(param(i))
       IF( elset%locax < 1 .OR. elset%locax > 3 )THEN
         WRITE(lures,"(/,5X,'Error in the definition of local system', / &
                         5X,'Invalid axis: ',i3,' Default value used (3)')") elset%locax
         elset%locax = 3
       END IF
     ELSE
       elset%locax = 3
     END IF
     elset%stabq=getrea('STABQ ',1d0,' Stabilization factor for shear    ')
     elset%quad = exists('QUAD  ')
     IF(elset%quad)THEN
       WRITE(lures,"(/,5X,' Quadratic approach for membrane ', &
                          &      'behavior Selected'//)",ERR=9999)
       elset%nnode = 9
     ELSE
       elset%nnode = 6
     END IF
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements
     !Initialize empty list Point both pointer to nothing
     CALL ini_ele03e (elset%head, elset%tail)
     IF( exists('ZIGZAG') ) THEN !local axis definition
       IF( ndofn == 8 )THEN
         elset%zigzag = .TRUE.
         elset%nstre  = 14
       END IF
     END IF
   END IF

   CALL elmd03(nelem,elset%head,elset%tail,iwrit,elset%quad,elset%nnode,elset%nstre)

   IF(cmpea)THEN
     ALLOCATE (lnods(6,nelem))
     el => elset%head
     DO i=1,nelem
       lnods(:,i) = el%lnods(1:6)
       el => el%next
     END DO
     CALL nodn03(nelem,lnods,coord,eule0,euler)  !locax
     DEALLOCATE (lnods)
   END IF

   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs ( 1 ,nreqs, elset%ngrqs, iwrit )

   CALL comm03(0, nelem,  nreqs, narch, elsnam, elset)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele03 (elset, head, tail) ! add to the list of sets
     nelms = nelms + 1 ! increased set counter for this element type
   END IF
   nel = nel + nelem

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)          !initializes a list
   NULLIFY(elset%head)       !nullify head pointer
   ! read control parameters
   elset%sname = elsnam
   READ (51) elset%nelem, elset%nreqs, elset%narch, elset%gauss, elset%lside, &
             elset%quad, elset%plstr, elset%angdf , elset%stabq, elset%nnode, elset%zigzag, elset%nstre
   ! restore list of elements
   CALL rest03 (elset%nelem, elset%nreqs, elset%head, elset%tail, &
                elset%ngrqs, elset%quad, elset%nnode, elset%nstre )
   ! add to list of elements
   CALL add_ele03 (elset, head, tail)

 ELSE IF (TRIM(task) == 'IMPORT') THEN
   ! check if list of sets and set exists and initializes
   IF( overw )THEN
     CALL srch_ele03 (head, anter, elset, elsnam, oldset)
     IF (oldset) THEN    !if set exists
       CALL comm03(1, nelem,  nreqs, narch, elsnam, elset)
       !elset%lside = .FALSE.   !initializes flag to compute LSIDE
     ELSE
       CALL runen2(' Old set to overwrite does not exist')
     END IF
     READ (fimpo)
   ELSE
     ALLOCATE (elset)       !reserve memory for set
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     elset%lside = .FALSE.  !     "     flag to compute LSIDE
     elset%plstr = 0     ! do not compute plastic strains
     elset%locax = 3        ! must be modified lateer
     CALL ini_ele03e(elset%head,elset%tail)
     elset%sname = sname
     READ (fimpo) nelem,elset%nnode,ng,stype,elset%angdf,elset%stabq
     elset%quad = elset%nnode == 9

     nreqs = 0
     narch = 0
   END IF
   ! restore list of elements
   CALL impo03 ( nelem,elset%head, elset%tail, elset%nnode)
   CALL comm03 (0, nelem,  nreqs, narch, elsnam, elset)
   ! add to list of elements
   IF( .NOT.overw )THEN
     CALL add_ele03 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
     nel = nel + nelem
   END IF
 ELSE
   CALL runend('INPD03: NON-EXISTENT TASK .        ')
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpd03
