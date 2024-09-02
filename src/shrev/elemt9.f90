 SUBROUTINE elemt9(task ,ttime,iload,iforce,                       &
                   coora,locsy,sumat,mass,emass,resid,gstif,       &
                   istop,flag1,flag2,iwrit)

 USE ele09_db
 USE ctrl_db, ONLY: ndime, npoin, top, bottom, ntype
 USE npo_db, ONLY : naeul, ifpre, coort, coorb, ifact, coord, euler, loadv, eule0
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

 REAL (kind=8),SAVE :: energ(9)
 TYPE (ele09_set), POINTER :: elset, anter
 INTEGER (kind=4) :: nelem, nnode, ngaus, nstre, axesc, nreqs, narch
 CHARACTER (len=mnam) :: sname


 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL commv9 (1, nelem, nnode, ngaus, nstre, axesc, &
                   nreqs, narch, sname, elset)


   SELECT CASE (TRIM(task))

   CASE ('NEW','NSTRA1','NSTRA2')
     CALL acvdf9(nnode,ifpre,elset%head,naeul,elset%zigzag,elset%zigzpp)

   CASE('ACTUAL')
     CALL actua9(elset%head,ntype)

   !CASE ('CLOSEF')
   !
   !  CALL close1(nreqs,narch)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump09 ( elset )

   CASE ('GAUSSC')
     IF( ntype <= 3 )THEN
     CALL gauss9(ndime,ntype,nstre,nnode,ngaus,axesc,coord,eule0,      &
                 istop,elset%head,elset%gauss, &
                 elset%posgp,elset%weigh,elset%shape,elset%deriv,elset%zigzag,elset%zigzpp)
     ELSE
      CALL gauss9c(ndime,nnode,ngaus,coord,            &
                   istop,elset%head,elset%gauss,             &
                   elset%posgp,elset%weigh,elset%shape,elset%deriv)
     END IF

   CASE ('LOADPL')
     CALL loadp9 (ntype,nelem,loadv(:,:,iload),gv,gravy,nnode,ngaus, &
                  elset%head,elset%shape,elset%weigh,elset%posgp)


   CASE ('MASMTX')

     CALL masmt9(flag1,nelem,nnode,ngaus,elset%shape,elset%weigh, &
                 emass,mass,sumat,elset%head,ntype,ifpre,elset%ndofe,elset%zigzag)

   CASE ('OUTPUT')

     IF(flag1 .OR. flag2) THEN     !flag1 = b1    flag2 = b2   argm3 = ttime
       CALL outpu9(flag1,flag2,nelem,ngaus,nreqs,narch,iwrit,ntype, &
                   nstre,elset%ngrqs,ttime,elset%head,elset%zigzag,elset%zigzpp,elset%estr)
       !IF(flag1) WRITE(55,'(8e12.4)')energ/2d0   !print components
     END IF

   CASE ('WRTPOS')

     CALL masel9 (ntype,nnode,nstre,ngaus,nreqs,nelem,narch, &
                  elset%ngrqs,elset%sname,elset%head,elset%zigzag)
     elset%narch = narch

   CASE ('RESVPL')
     energ = 0d0
     IF( ntype /= 4 )THEN
     CALL resvp9(ntype,nelem,nnode,elset%ndofe,ngaus,axesc,nstre,elset%head,        &
                 coora,locsy,resid,elset%weigh,elset%shape,elset%deriv,elset%estr,  &
                 energ,istop,bottom,top,coorb,coort,ifact,elset%zigzag,elset%zigzpp,ttime,elset%esecs)
     ELSE
     CALL resvp9c(nelem,ngaus,elset%head,coord,euler,       &
                 coora,locsy,resid,elset%posgp)
     END IF

   CASE('STIFFM')
     CALL stiff9(nelem,nnode,elset%ndofe,ngaus,nstre,axesc,ntype,elset%shape,  &
                 elset%weigh,elset%deriv,elset%head,               &
                 coora, locsy, gstif, resid(:,iforce),elset%zigzag,elset%zigzpp,elset%esecs)

   END SELECT
   elset => elset%next
 END DO
 RETURN
 END SUBROUTINE elemt9
