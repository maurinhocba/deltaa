 SUBROUTINE elemt2(task ,ttime,iload,iforce,                       &
                   coora,sumat,mass,emass,resid,gstif,             &
                   istop,flag1,flag2,iwrit)

 USE ctrl_db, ONLY : ndime
 USE npo_db, ONLY : coord,loadv, ifpre
 USE loa_db, ONLY : igrav,gravy,gv
 USE ele02_db
 IMPLICIT NONE
 ! Dummy arguments
 CHARACTER(len=*) :: task
 INTEGER (kind=4), INTENT(IN), OPTIONAL  :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN), OPTIONAL  :: coora(:,:)
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),emass(:,:),resid(:,:),gstif(:)
 LOGICAL, OPTIONAL  :: flag1,flag2

 ! local variables
 INTEGER (kind=4) nelem,nreqs,narch
 CHARACTER(len=mnam) elsnam
 TYPE (ele02_set), POINTER :: elset, anter

 
 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO

   CALL commv2 (1,   nelem,nreqs,narch,elsnam,elset)

   SELECT CASE (TRIM(task))

   CASE ('NEW   ','NSTRA1','NSTRA2')

     CALL acvdf2(ndime,nelem,elset%head,ifpre)

   CASE ('ACTUAL')
     CALL actua2(nelem,elset%head)

   !CASE ('CLOSEF')
   !CALL close1(nreqs,narch)

   CASE ('DUMPIN','dumpin')
     CALL dumpi2 ( elset )

   CASE ('GAUSSC')
     CALL gauss2(ndime,nelem,elset%head, coora,flag1,istop)

   CASE ('LOADPL')
     CALL loadp2(ndime,nelem,igrav,loadv(:,:,iload),gv,gravy,elset%head)

   CASE ('MASMTX')
     CALL masmt2(ndime,nelem,elset%head,emass,flag1,mass,sumat,ifpre)

   CASE ('OUTPUT')
     IF(flag1.OR.flag2)                                              &
  &    CALL outpu2(flag1,flag2,nelem,nreqs,narch,iwrit,elset%head,  &
  &                ndime,elset%ngrqs,ttime)

   CASE ('RESVPL')
     CALL resvp2(ndime,nelem,coora,resid,elset%head,ttime,coord)

   CASE ('WRTPOS')
     CALL masel2(nreqs,nelem,narch,elset%head,elset%ngrqs,elsnam)
     elset%narch = narch

   CASE ('STIFFM')
     CALL stiff2(ndime,nelem,elset%head,coora,gstif,resid(:,iforce),   &
  &              emass,ttime)

   END SELECT
   IF ( ASSOCIATED (elset%next) ) THEN
     elset => elset%next
   ELSE
     EXIT
   END IF

 END DO
 RETURN
 END SUBROUTINE elemt2
