 SUBROUTINE inpda7 (task,nel, eule0,euler,coord, iwrit,elsnam,nelms)
 !******************************************************************
 !
 !*** READ control DATA for 6-node triangular shell element
 !
 !******************************************************************

 USE npo_db, ONLY : cpx
 USE ctrl_db, ONLY : ndofn
 USE ele07_db

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam
 CHARACTER(len=*),INTENT(IN):: task
 INTEGER (kind=4) :: nel,nelms,iwrit
 REAL    (kind=8) :: coord(:,:),eule0(:,:),euler(:,:)


 ! local variables
 LOGICAL :: cmpea, &   !CoMPute Euler Angles
            oldset     !
 INTEGER (kind=4) :: nreqs, narch,i,ansmm,nnass,g,nelem
 INTEGER (kind=4), ALLOCATABLE  :: lnods(:,:)
 TYPE (ele07_set), POINTER :: elset,anter
 TYPE (ele07), POINTER :: el
 REAL (kind=8) :: a,b

 INTERFACE
   INCLUDE 'nodnor.h'
 END INTERFACE

 ! ALLOCATE (elset)

 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele07 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL commv7 (1, nelem,  nreqs, narch, ansmm, nnass, elsnam, elset)
     cmpea = .FALSE.
     nel = nel - nelem
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
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
     IF( exists('GPINT  ',i)) THEN   !Internal Gauss point
       WRITE(lures,"(/,5X,' Interior Gauss points will be used')",ERR=9999)
       elset%gpint = .TRUE.
       a = 2d0/3d0
       b = 1d0/6d0
       elset%posgp = RESHAPE ((/  b , b , a , b , b , a /),(/2,ngaus/)) ! Gauss point position
     ELSE
       WRITE(lures,"(/,5X,' Mid-side Gauss points will be used')",ERR=9999)
       elset%gpint = .FALSE.
       elset%posgp = RESHAPE ((/ 0.5d0, 0.0d0, 0.5d0, 0.5d0, 0.0d0, 0.5d0  /),(/2,ngaus/)) ! Gauss point position
     END IF

     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements
     !Initialize empty list Point both pointer to nothing
     CALL ini_ele07e (elset%head, elset%tail)

     ! assumed membrane model

     ansmm = getint('ANSMM ',0,' ANS Membrane Model for this set ..')
     IF( ansmm < 0  .OR. ansmm > 2 )ansmm = 0 !check
     IF( ansmm /= 0 ) CALL omat07(ansmm,elset%posgp,elset%omat)
     IF( exists('SHQUAD') ) THEN !MITC6 quadratic approach for shear
       WRITE(lures,"(/,5X,' Quadratic approach for shear strains will be used')",ERR=9999)
       nnass = nnas2
     ELSE IF( exists('STDSH') ) THEN !use STandarD displacement formulation for SHear
       WRITE(lures,"(/,5X,' Standard approach for shear strains (locks) will be used')",ERR=9999)
       nnass = ngaus
     ELSE
       WRITE(lures,"(/,5X,' Linear approach for shear strains will be used')",ERR=9999)
       nnass = nnas1
     END IF
     ALLOCATE( elset%ap1(2,nnass,ngaus) )
     IF( ansmm < 0  .OR. ansmm > 2 )ansmm = 0 !check
     IF( ansmm /= 0 ) CALL omat07(ansmm,elset%posgp,elset%omat)
     IF( exists('ZIGZAG') ) THEN !local axis definition
       IF( ndofn == 8 )THEN
         elset%zigzag = .TRUE.
         elset%nstre  = 14
       END IF
     END IF
   END IF

   !  ******************************

   CALL elmda7(nelem,elset%head,elset%tail,ansmm,nnass,iwrit,elset%nstre)
   IF(cmpea)THEN !This requires mid-side coordinates are given as data
     ALLOCATE (lnods(nnode,nelem))
     el => elset%head
     DO i=1,nelem
       lnods(:,i) = el%lnods
       el => el%next
     END DO
     CALL nodnor(nelem,nnode,lnods,coord,eule0,euler)
     DEALLOCATE (lnods)
   END IF
   CALL nodxy7(nelem,elset%head,coord,eule0,euler,cpx)

   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs ( 1 ,nreqs, elset%ngrqs, iwrit )

   CALL commv7(0, nelem,  nreqs, narch, ansmm, nnass, elsnam, elset)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele07 (elset, head, tail) ! add to the list of sets
     nelms = nelms + 1 ! increased set counter for this element type
   END IF
   nel = nel + nelem

! ELSE IF (TRIM(task) == 'RESTAR') THEN
!
!   ALLOCATE (elset)          !initializes a list
!   NULLIFY(elset%head)       !nullify head pointer
!   ! read control parameters
!   elset%sname = elsnam
!   READ (51) elset%nelem, elset%nreqs, elset%narch, elset%ansmm, &
!             elset%gauss, elset%plstr, elset%angdf, elset%zigzag, elset%nstre
!   ! restore list of elements
!   CALL rest07 (elset%nelem, elset%nreqs, elset%head, elset%tail, &
!                elset%ngrqs, elset%ansmm, elset%nstre)
!   ! add to list of elements
!   CALL add_ele07 (elset, head, tail)
!
 ELSE
   CALL runend('INPD07: NON-EXISTENT TASK .        ')
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpda7
