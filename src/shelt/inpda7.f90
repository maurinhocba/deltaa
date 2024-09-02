 SUBROUTINE inpda7 (task,nel, eule0,euler,coord, iwrit,elsnam,nelms)
 !******************************************************************
 !
 !*** READ control DATA for 4-node cuadrilateral shell element
 !
 !******************************************************************

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
 INTEGER (kind=4) :: nreqs, narch,i,stype,g,nelem
 INTEGER (kind=4), ALLOCATABLE  :: lnods(:,:)
 TYPE (ele07_set), POINTER :: elset,anter
 TYPE (ele07), POINTER :: el
 REAL (kind=8) :: a,b,deriv(nnode,2)

 INTERFACE
   INCLUDE 'nodnor.h'
 END INTERFACE

 ! ALLOCATE (elset)

 !       gauss points in local coordinates and weigths

 ! assumed membrane strain points

 a = 1d0/3d0
 CALL shape7(0d0,0d0,deriv,dn(1,1,1,1))    !node 1 (vertex)
 CALL shape7(1d0,0d0,deriv,dn(1,1,2,1))    !node 2 (vertex)
 CALL shape7(0d0,1d0,deriv,dn(1,1,3,1))    !node 3 (vertex)
 CALL shape7(0.5d0,0.0d0,deriv,dn(1,1,1,2))    !
 CALL shape7(0.5d0,0.5d0,deriv,dn(1,1,2,2))    !
 CALL shape7(0.0d0,0.5d0,deriv,dn(1,1,3,2))    !

 ! ******************************
 ! gauss points in local coordinates and weigths
 a = 1d0/6d0                   !
 weigp = a                     !weight
 ! version for Integration points at mid-side nodes
 !b = 1d0/2d0                   !coordinate at mid point
 !posgp(:,1) = (/ b  ,0d0 /)    !
 !posgp(:,2) = (/ b  ,b   /)    !
 !posgp(:,3) = (/ 0d0,b   /)    !
 !! version for Integration points at interior points
 b = 2d0/3d0                   !coordinate at interior point
 posgp(:,1)= (/ a,a /)
 posgp(:,2)= (/ b,a /)
 posgp(:,3)= (/ a,b /)

 ! gauss points shape functions and derivatives of nodal functions

 DO g=1,ngaus
   CALL shape7(posgp(1,g),posgp(2,g),shape(1,g),deriv) !at integration points
   CALL ap1tm7(ap1(:,:,g),posgp(1,g),posgp(2,g))       !(AP)^T
 END DO

 !  ******************************
 !  assumed shear strain points

 a = (1d0-1d0/SQRT(3d0))/2d0                         !first point along side
 b = (1d0+1d0/SQRT(3d0))/2d0                         !second point along side
 CALL shape7(a  ,0d0,deriv(1,1),nfdas(1,1,1))        !A
 CALL shape7(b  ,0d0,deriv(1,2),nfdas(1,1,2))        !B
 CALL shape7(b  ,a  ,deriv(1,1),nfdas(1,1,3))        !C
 CALL shape7(a  ,b  ,deriv(1,2),nfdas(1,1,4))        !D
 CALL shape7(0d0,b  ,deriv(1,1),nfdas(1,1,5))        !E
 CALL shape7(0d0,a  ,deriv(1,2),nfdas(1,1,6))        !F

 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele07 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL commv7 (1, nelem,  nreqs, narch, stype, elsnam, elset)
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
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements
     !Initialize empty list Point both pointer to nothing
     CALL ini_ele07e (elset%head, elset%tail)
     stype = getint('STYPE ',0,' Formulation type for this set ....')
     IF( stype /= 2 .AND. stype /= 3 )stype = 0 !check
     IF( exists('ZIGZAG') ) THEN !local axis definition
       IF( ndofn == 8 )THEN
         elset%zigzag = .TRUE.
         elset%nstre  = 14
       END IF
     END IF
   END IF

   CALL elmda7(nelem,elset%head,elset%tail,stype,iwrit,elset%nstre)
   CALL nodxy7(nelem,elset%head,coord,eule0,euler)

   IF(cmpea)THEN
     ALLOCATE (lnods(nnode,nelem))
     el => elset%head
     DO i=1,nelem
       lnods(:,i) = el%lnods
       el => el%next
     END DO
     CALL nodnor(nelem,nnode,lnods,coord,eule0,euler)
     DEALLOCATE (lnods)
   END IF

   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs ( 1 ,nreqs, elset%ngrqs, iwrit )

   CALL commv7(0, nelem,  nreqs, narch, stype, elsnam, elset)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele07 (elset, head, tail) ! add to the list of sets
     nelms = nelms + 1 ! increased set counter for this element type
   END IF
   nel = nel + nelem

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)          !initializes a list
   NULLIFY(elset%head)       !nullify head pointer
   ! read control parameters
   elset%sname = elsnam
   READ (51) elset%nelem, elset%nreqs, elset%narch, elset%stype, &
             elset%gauss, elset%plstr, elset%angdf, elset%zigzag, elset%nstre
   ! restore list of elements
   CALL rest07 (elset%nelem, elset%nreqs, elset%head, elset%tail, &
                elset%ngrqs, elset%stype, elset%nstre)
   ! add to list of elements
   CALL add_ele07 (elset, head, tail)

 ELSE
   CALL runend('INPD07: NON-EXISTENT TASK .        ')
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpda7
