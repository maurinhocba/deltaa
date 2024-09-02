 SUBROUTINE elemt8(task ,ttime,iload,iforce,                       &
                   coora,locsy,sumat,mass,emass,resid,gstif,       &
                   istop,flag1,flag2,iwrit)

 USE ele08_db
 USE ctrl_db, ONLY: ndime, ndofn, npoin
 USE npo_db, ONLY : naeul, ifpre, coord, euler, coort, coorb, ifact, loadv, eule0, velnp
 USE loa_db, ONLY : gravy, gv, igrav

 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN) :: task
 INTEGER (kind=4), INTENT(IN), OPTIONAL  :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime,coora(:,:),locsy(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),resid(:,:),gstif(:),     &
                                  emass(:,:)
 LOGICAL, OPTIONAL  ::  flag1,flag2


 !REAL (kind=8),SAVE :: energ(8)
 TYPE (ele08_set), POINTER :: elset, anter
 INTEGER (kind=4) nelem,nnode,ngaus,axesc,narch,nreqs
 CHARACTER(len=mnam) :: sname

 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO

   CALL commv8 (1,     nelem, nnode, ngaus, axesc,               &
                nreqs, narch, sname, elset)

   SELECT CASE (TRIM(task))

   CASE ('NEW','NSTRA1','NSTRA2')
   CALL acvdf8(nelem,nnode,ifpre,elset%head,naeul)

   CASE('ACTUAL')
     CALL actua8(elset%head)

   CASE ('GAUSSC')
   CALL gauss8(ndime,nelem,nnode,ngaus,elset%axesc,elset%head,coord,   &
               eule0,elset%posgp,elset%shape,elset%deriv,        &
               elset%weigh,flag1,istop)

   CASE ('LOADPL')
   CALL loadp8(ndime,ndofn,nelem,loadv(:,:,iload),               &
               gv,gravy,nnode,ngaus,axesc,elset%head,           &
               elset%shape,elset%weigh)

   CASE ('MASMTX')
   CALL masmt8(ndime,nelem,iwrit,nnode,ngaus,axesc,elset%head,            &
               elset%weigh,emass,elset%shape,sumat,mass,flag1,ifpre)

   CASE ('OUTPUT')
   IF(flag1.OR.flag2)                                            &
     CALL outpu8(flag1,flag2,nelem,nreqs,narch,iwrit,ngaus,      &
                 elset%ngrqs,ttime,elset%head)

   CASE ('RESVPL')
   CALL resvp8(ndime,nelem,nnode,ngaus,axesc,coora,              &
               locsy,velnp,resid,elset%weigh,elset%shape,        &
               elset%deriv,elset%head,istop)

   CASE ('DUMPIN')
     CALL dumpi8 (sname, nelem, nnode, ngaus, axesc, nreqs,         &
                  narch, elset%head, elset%ngrqs,                   &
                  elset%posgp, elset%shape, elset%deriv, elset%weigh)

   CASE ('WRTPOS')
   CALL masel8(nnode,ngaus,nreqs,nelem,narch,elset%head,              &
               elset%ngrqs,elset%sname)
    elset%narch = narch

   CASE('STIFFM')
     CALL stiff8(ndime,nelem,nnode,ngaus,axesc,elset%shape,  &
                 elset%weigh,elset%deriv,elset%head,               &
                 coora, locsy, gstif, resid(:,iforce))

   !CASE ('DELETE')
   !  flag2 = TRIM(elsnam) == TRIM(sname)
   !  IF (flag2) THEN
   !    CALL del_ele08 (head, tail, anter, elset)
   !    nelms(8) = nelms(8) - 1
   !    EXIT
   !  END IF
   !  IF ( ASSOCIATED (elset%next) ) anter => elset
   !
   !CASE ('SURFAC')      !compute contact surface
   !  IF( flag2 )EXIT
   !  flag2 = TRIM(elsnam) == TRIM(sname)
   !  IF (flag2) THEN
   !    ! get surface definition from the element set
   !    CALL surf08 (elset%head,elset%nelem)
   !    EXIT
   !  END IF

   END SELECT

   IF ( ASSOCIATED (elset%next) ) THEN
     elset => elset%next
   ELSE
     EXIT
   END IF

 END DO

 RETURN

 END SUBROUTINE elemt8
