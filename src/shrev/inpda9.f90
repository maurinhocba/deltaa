 SUBROUTINE inpda9 (task,nel, eule0,euler,coord, iwrit,elsnam,nelms)
!******************************************************************
!
!*** READ control DATA for 2-3-node 2-d beam/shell element
!
!******************************************************************

 USE ctrl_db, ONLY : ntype,ndofn
 USE ele09_db

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam
 CHARACTER(len=*),INTENT(IN):: task
 INTEGER (kind=4) :: nel,nelms,iwrit
 REAL    (kind=8) :: coord(:,:),eule0(:,:),euler(:,:)

 INTEGER (kind=4) :: nnode, ngaus, nstre, axesc, nreqs, narch, nelem, ip

 CHARACTER (len=12) :: ptype(4) =(/ 'Plane stress', &
                                    'Plane strain', &
                                    'Axisymmetric', &
                                    '2-D C-beam  '  /)

 LOGICAL ::  oldset
 TYPE (ele09_set), POINTER :: elset,anter
 TYPE (ele09), POINTER :: elem

 !ALLOCATE (elset)

 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele09 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL commv9 (1, nelem, nnode, ngaus, nstre, axesc, &
                     nreqs, narch, elsnam, elset)
     nel = nel - nelem
   ELSE                !set ELSNAM does not exist
     CALL new_ele09 (elset)       !reserve memory for set
   !     READ the set DATA card
     CALL listen('INPDA9')
     WRITE(lures,"(/,5x,'Control parameters for beam element'//)",ERR=9999)
     nnode=getint('NNODE ',2,' NUMBER OF NODES PER ELEMENT ......')
     ngaus=getint('NGAUS ',nnode-1,' Integration points in the set ....')
     axesc=getint('AXESCO',0,' Local axes code ..................')
     elset%ndofe = 3
     IF( exists('ZIGZAG') ) THEN !local axis definition
       IF( ndofn == 4 )THEN
         elset%zigzag = .TRUE.
         elset%ndofe = 4
       END IF
     ELSE IF (exists('ZIGZ++'))THEN
       IF( ndofn == 7 )THEN
         elset%zigzpp = .TRUE.
         elset%ndofe = 7
       END IF
     END IF
     SELECT CASE (ntype)
     CASE (1)
       nstre = 3
     CASE (2:3)
       nstre = 5
     CASE (4)
       nstre = 5
       IF( nnode /= 2 )THEN
         WRITE(lures,"(/,5x,'Linear beam must have 2 nodes',/,5x,'Changed')",ERR=9999)
         nnode = 2
       END IF
       IF( elset%zigzag .OR. elset%zigzpp )THEN
         WRITE(lures,"(/,5x,'Linear beam Cannot include Zigzag',/,5x,'Changed')",ERR=9999)
         elset%zigzag = .FALSE.
         elset%ndofe = 3
       END IF
       IF( ngaus < 3 ) ngaus = 5
       axesc = 0
     END SELECT
     nreqs=getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )

     IF(ABS(axesc) >= 2) axesc= nnode*axesc/ABS(axesc)
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements
     !Initialize empty list Point both pointer to nothing
     CALL ini_ele09e (elset%head, elset%tail)
     ALLOCATE( elset%weigh(ngaus),elset%posgp(ngaus), &
               elset%shape(nnode,ngaus),elset%deriv(nnode,ngaus))
   END IF
   CALL elmda9(nelem,nnode,nstre,ngaus,axesc,elset%head,elset%tail,iwrit)
   WRITE(lures,"(10x,A,/,10x,'No of stresses  (NSTRE) =',i3,/)",ERR=9999)    &
         TRIM(ptype(ntype)),nstre
   IF( elset%zigzag .OR. elset%zigzpp )THEN
     ALLOCATE( elset%estr(nstre*ngaus,nelem))  !keeps element strains to compute and print element profiles
     ALLOCATE( elset%esecs(nelem))
     elem => elset%head
     DO ip=1,nelem
       elset%esecs(ip) = elem%matno
       elem => elem%next
     END DO
   ELSE
     NULLIFY(elset%esecs, elset%estr )
   END IF

   IF(nnode == 3)CALL nodxy9(nelem,elset%head,coord,eule0,euler)
   CALL locla9(nelem,nnode,axesc,elset%head,iwrit,eule0,coord)

   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) THEN
     CALL rdreqs ( 1 ,nreqs, elset%ngrqs, iwrit )
     CALL add_ele09 (elset, head, tail) ! add to the list of sets
     nelms = nelms + 1 ! increased the counter of element sets of this type
   END IF
   nel = nel + nelem

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   READ (51) nnode, ngaus, nstre, axesc, nreqs, narch

   ALLOCATE( elset%weigh(ngaus),elset%posgp(ngaus), &
             elset%shape(nnode,ngaus),elset%deriv(nnode,ngaus))
   CALL rest09(nelem,nreqs,nnode,ngaus,nstre,axesc,elset%head,elset%ngrqs, &
               elset%posgp,elset%shape,elset%deriv,elset%weigh)
   CALL add_ele09 (elset, head, tail)

 ELSE
   CALL runend('INPDA9: NON-EXISTENT TASK .        ')
 ENDIF

 CALL commv9 (0, nelem, nnode, ngaus, nstre, axesc, &
                 nreqs, narch, elsnam, elset)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpda9
