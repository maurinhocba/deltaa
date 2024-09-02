 SUBROUTINE elem11(task ,ttime,iload,iforce,                &
                   coora,sumat,mass,emass,resid,gstif,       &
                   istop,flag1,flag2,iwrit)

 !master routine for element 11 rotation free 2-D beam/shell element
 USE ele11_db
 USE loa_db, ONLY : igrav,gravy,gv
 USE ctrl_db, ONLY : ndime,npoin,ntype,top,bottom
 USE npo_db,  ONLY : label,oldlb,coord,coort,coorb,ifact, &
                     iffix,ifpre,loadv

 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN) :: task

 INTEGER (kind=4), INTENT(IN), OPTIONAL :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime,coora(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),resid(:,:),gstif(:),     &
                                  emass(:,:)
 LOGICAL, OPTIONAL ::  flag1,flag2

 TYPE (ele11_set), POINTER :: elset, anter
 INTEGER (kind=4) :: nelem, nstre, nbn, nreqs, narch, ngaus
 CHARACTER (len=mnam) :: sname

 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL comm11 (1, nelem, nstre, nbn, ngaus, &
                   nreqs, narch, sname, elset)


   SELECT CASE (TRIM(task))

   CASE ('NEW','NSTRA1','NSTRA2')
     CALL acvd11(ifpre,elset%head,elset%lside,nelem,elset%nbn,elset%nhead)

   CASE('ACTUAL')
     CALL actu11(elset%head,nstre,ngaus)

   !CASE ('CLOSEF')
   !
   !  CALL close1(nreqs,narch)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump11 ( elset )

   CASE ('GAUSSC')
     CALL gaus11(nstre,nbn,ntype,ngaus,iffix,coord,      &
                 istop,elset%shap,elset%head,elset%gauss, &
                 elset%nhead,elset%strai)

   CASE ('LOADPL')
     CALL load11 (nelem,loadv(:,:,iload),gv,gravy,ntype,elset%head,coord)


   CASE ('MASMTX')

     CALL masm11(ifpre,flag1,nelem,emass,mass,sumat,elset%head,ntype,coord)

   CASE ('OUTPUT')

     IF(flag1 .OR. flag2) THEN     !flag1 = b1    flag2 = b2   argm3 = ttime
       CALL outp11(flag1,flag2,nelem,nstre,nreqs,narch,ntype,ngaus,iwrit, &
                   elset%ngrqs,ttime,elset%head,elset%stint,elset%shap)

     END IF

   CASE ('WRTPOS')

     CALL mase11 (nstre,nreqs,nelem,narch,ntype,ngaus, &
                  elset%ngrqs,elset%sname,elset%head)
     elset%narch = narch

   CASE ('RESVPL')

     CALL resv11(nelem,nbn,nstre,ntype,ngaus,elset%head,elset%nhead,elset%shap,   &
                 iffix,coora,resid,istop,ttime,bottom,top,coorb,coort,ifact,      &
                 elset%stabs,elset%stint)

   CASE('STIFFM')
     CALL stif11(elset%head, coora, gstif, resid(:,iforce), elset%nhead, &
                 ntype, nstre,ngaus,elset%shap,elset%stabs,elset%stint)

   END SELECT
   elset => elset%next
 END DO
 RETURN
 END SUBROUTINE elem11
