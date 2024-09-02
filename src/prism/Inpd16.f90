 SUBROUTINE inpd16 (task, nel, iwrit, elsnam, nelms)

 !   READ control DATA for element number 16 (TL PRISM 6-15)

 USE ele16_db
 USE gvar_db,ONLY : fimpo,overw
 USE npo_db, ONLY : cpx,coord
 USE ctrl_db, ONLY : npoin

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets
                     nel,     & ! number of elements in the set
                     iwrit      ! flag to echo data input

 ! local variables
 LOGICAL :: oldset,quad,shell,bbar,bezier
 INTEGER (kind=4) :: nreqs, narch, nn, ngaus, nassp, nelem, i,ngaup
 TYPE (ele16_set), POINTER, SAVE  :: elset, anter

 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele16 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm16 (1, nelem,  nreqs, narch, elsnam, elset, nn, ngaus,nassp,quad,shell,bbar,bezier)
     nel = nel - nelem
   ELSE                !set ELSNAM does not exist
     nassp = 0
     ALLOCATE (elset)       !reserve memory for set
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     CALL listen('INPD16')  !read a line
     nn    = getint('NNODE ',6,' NUMBER OF ELEMENT NODES (6/12/15/18)')
     IF( nn /= 6 .AND. nn /= 12 .AND. nn /= 15 .AND. nn /= 18 )CALL runend('PRISM: NNODE must be 6, 15 or 18    ')
     shell = exists('SHELL')
     IF( nn == 6 )THEN
       ngaus = 2
       ngaup = 1
       quad = exists('QUAD')
       IF(quad) WRITE(lures,"(9X,'QUADratic approximation will be used.')")
       elset%lface = .FALSE.
       bbar  = exists('BBAR ')
       IF(bbar )WRITE(lures,"(' Average volumetric strain will be used')",ERR=9999)
       IF(shell)WRITE(lures,"(' Asumed transverse shear strains will be used')",ERR=9999)
       nassp = 3
     ELSE
       quad = .FALSE.
       bbar  = exists('BBAR ') .AND. nn == 12
       IF(exists('SHQUAD')) THEN
         nassp = 8
         WRITE(lures,"(' Asumed transverse shear QUAD-strains will be used')",ERR=9999)
         shell = .TRUE.
       ELSE IF( shell) THEN
         IF(shell)WRITE(lures,"(' Asumed transverse Linear-shear strains will be used')",ERR=9999)
         nassp = 6
       END IF
       bezier = exists('BEZIER')
       IF( bezier ) THEN
         IF( .NOT.ASSOCIATED (cpx) )THEN
           ALLOCATE (cpx(3,npoin))
           cpx = 0
         END IF
         WRITE(lures,"(' BEZIER quadratic polynomials will be used')",ERR=9999)
       END IF
       IF( nn == 12 )ngaus = 6
       IF( nn == 15 )ngaus = 6
       IF( nn == 18 )ngaus = 8
       ngaus = getint('NGAUS ',ngaus,' NUMBER OF GAUSS POINTS  (6/7/8/9 ONLY)')
       IF( ngaus < 6 .OR. ngaus > 9)CALL runend('PRISM: NGAUS must be 6..9 for quadratic prism')
       IF( shell .AND. ngaus == 7 ) ngaus = 6
       ngaup = 3
       IF( MOD(ngaus,4) == 0 )ngaup = 4
     END IF
     IF( shell )THEN
       ALLOCATE( elset%gpa(2,nassp),elset%nfnda(nn,nassp,2),elset%pag(2,nassp,ngaup),elset%amat(nassp,nassp))
       IF( nassp == 3 ) elset%gpa = gpa3
       IF( nassp == 6 .OR. nassp == 8 ) elset%gpa = gpa8(:,1:nassp)
       IF( nassp == 6 )THEN
         elset%amat = amal
       ELSE IF (nassp == 8 )THEN
         elset%amat = amaq
       END IF
     END IF
     IF(quad .OR. shell ) THEN
       elset%locax = 3
     ELSE
       elset%locax = 0
     END IF
     IF( exists('LOCAL  ',i)) THEN   !local axis definition
       IF( param(i)  > 0 .AND. param(i) < 4 )THEN  !if a valid code
         elset%locax =  INT(param(i))
       ELSE  !keep default value and p
         WRITE(lures,"(/,5X,'Error in the definition of local system', / &
                         5X,'Invalid axis: ',i3,' Default value used (3)')") INT(param(i))
       END IF
     END IF
     nreqs=getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     elset%angdf(1) =getrea('ALPHA ',0d0,' First Euler Angle from X and Ortho')
     elset%angdf(2) =getrea('BETA  ',0d0,' Second Euler Angle from X and Orth.')
     elset%angdf(3) =getrea('GAMMA ',0d0,' Third Euler Angle from X and Ortho')
     elset%small = exists('SMALL')
     IF(elset%small) WRITE(lures,"(' Green strains will be used if possible')")
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements

     !Initialize empty list Point both pointer no nothing
     CALL ini_ele16e (elset%head, elset%tail)
   END IF
   !  read new data or add to previous data
   CALL elmd16(nelem, nn, elset%head, elset%tail, iwrit, ngaus, quad, shell, nassp)
   IF(nn > 6) CALL nodx16(nelem,nn,elset%head,coord,cpx,bezier)
   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs (ngaus,nreqs, elset%ngrqs, iwrit )

   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele16 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF
   nel = nel + nelem

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)          !initializes a list
   NULLIFY(elset%head)       !nullify head pointer
   ! read control parameters
   elset%sname = elsnam
   READ (51) elset%nelem, elset%nnode, elset%nreqs, elset%narch, elset%gauss, &
             elset%plstr, elset%angdf, elset%ngaus, elset%small, elset%shell,&
             elset%quad, elset%locax , elset%bbar
   elset%lface = .TRUE.
   ! restore list of elements
   CALL rest16 (elset%nelem, elset%nnode, elset%nreqs, elset%head, elset%tail, &
                elset%ngrqs, elset%ngaus, elset%quad, elset%shell  )
   ! add to list of elements
   CALL add_ele16 (elset, head, tail)

 ELSE IF (TRIM(task) == 'IMPORT') THEN
   ! check if list of sets and set exists and initializes
   IF( overw )THEN
     CALL srch_ele16 (head, anter, elset, elsnam, oldset)
     IF (oldset) THEN    !if set exists
       CALL comm16(1, nelem,  nreqs, narch, elsnam, elset, nn, ngaus,nassp,quad,shell,bbar,bezier)
     ELSE
       CALL runen2(' Old set to overwrite does not exist')
     END IF
     READ (fimpo)
   ELSE
     ALLOCATE (elset)       !reserve memory for set
     CALL ini_ele16e (elset%head, elset%tail)
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     elset%lface = .FALSE.  !initializes flag to compute extended connectivities
     READ (fimpo) nelem,nn,ngaus,elset%angdf,elset%small,elset%quad, &
                  elset%bbar,elset%shell,elset%locax
     nreqs = 0
     narch = 0
   END IF
   ! restore list of elements
   CALL impo16 ( nelem,nn,ngaus, elset%head, elset%tail, elset%shell, elset%quad)
   ! add to list of elements
   IF( .NOT.overw )THEN
     CALL add_ele16 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE
   CALL runend('INPD16: NON-EXISTENT TASK .        ')
 END IF
 CALL comm16(0, nelem,  nreqs, narch, elsnam, elset, nn, ngaus,nassp,quad,shell,bbar,bezier)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpd16
