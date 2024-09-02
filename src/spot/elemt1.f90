 SUBROUTINE elemt1(task ,ttime,iload,iforce,                       &
                   coord,euler,sumat,mass,emass,resid,gstif,       &
                   istop,flag1,flag2,iwrit)

 USE ctrl_db, ONLY : ndime,neulr,ndofn
 USE loa_db, ONLY : igrav,gravy,gv
 USE npo_db, ONLY : ifpre,loadv
 USE ele01_db
 IMPLICIT NONE
 ! Dummy arguments
 CHARACTER(len=*) :: task
 INTEGER (kind=4), INTENT(IN), OPTIONAL  :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN), OPTIONAL  :: coord(:,:),euler(:,:)
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),emass(:,:),resid(:,:),gstif(:)
 LOGICAL, OPTIONAL  :: flag1,flag2

 ! local variables
 INTEGER (kind=4) nelem,nreqs,narch,nnode
 CHARACTER(len=mnam) elsnam
 TYPE (ele01_set), POINTER :: elset, anter

 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO

   CALL commv1 (1,   nnode,nelem,nreqs,narch,elsnam,elset)

   SELECT CASE (TRIM(task))

   CASE ('NEW   ','NSTRA1','NSTRA2')

     CALL acvdf1(ndime,nnode,ndofn,nelem,ifpre,elset%head)

!   CASE ('ACTUAL')
!     CALL actua1(nelem,elset%head)

   !CASE ('CLOSEF')
   !CALL close1(nreqs,narch)

   CASE ('DUMPIN','dumpin')
     CALL dumpi1 ( elset )

   CASE ('GAUSSC')
     CALL gauss1(ndime,nnode,neulr,nelem,elset%head,coord,euler)

   CASE ('LOADPL')
     IF( nnode == 2 )CALL loadp1(ndime,nelem,igrav,loadv(:,:,iload),gv,gravy,elset%head)

   CASE ('MASMTX')
     IF( nnode == 2 )CALL masmt1(ndime,ndofn,nelem,elset%head,emass,flag1,mass,sumat,ifpre)

   CASE ('OUTPUT')
     IF(flag1.OR.flag2)                                              &
       CALL outpu1(flag1,flag2,nnode,nelem,nreqs,narch,iwrit,elset%head,   &
                   ndime,elset%ngrqs,ttime)

   CASE ('RESVPL')
     IF( nnode == 2 )CALL resvp1(ndime,ndofn,nelem,coord,euler,resid,elset%head,ttime)

   CASE ('WRTPOS')
     CALL masel1(nreqs,nnode,nelem,narch,elset%head,elset%ngrqs,elsnam)
     elset%narch = narch

   CASE ('STIFFM')
     IF( nnode == 2 )CALL stiff1(ndime,ndofn,nelem,elset%head,coord,euler,gstif,   &
                 resid(:,iforce),ttime)

   END SELECT
   IF ( ASSOCIATED (elset%next) ) THEN
     elset => elset%next
   ELSE
     EXIT
   END IF

 END DO
 RETURN
 END SUBROUTINE elemt1
