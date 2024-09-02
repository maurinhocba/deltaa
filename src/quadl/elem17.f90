 SUBROUTINE elem17(task ,ttime,iload,iforce,                       &
                   coora,sumat,mass,emass,resid,gstif,ustif,       &
                   istop,flag1,flag2,iwrit)

 USE ctrl_db, ONLY : ntype,inverse
 USE loa_db, ONLY : igrav,gravy,gv
 USE npo_db, ONLY : label,coord,loadv,ifpre
 USE ele17_db
 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN) :: task
 INTEGER (kind=4), INTENT(IN), OPTIONAL  :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime,coora(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),resid(:,:),gstif(:),     &
                                ustif(:),emass(:,:)
 LOGICAL, OPTIONAL  ::  flag1,flag2

 INTEGER (kind=4), PARAMETER :: ndime=2

 TYPE (ele17_set), POINTER :: elset, anter
 INTEGER (kind=4) nelem, nreqs, narch, ngaus, nnode
 CHARACTER (len=mnam) :: sname


 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL comm17 (1,nelem,nreqs,narch,sname,elset,ngaus,nnode)

   SELECT CASE (TRIM(task))

   CASE ('NEW   ','NSTRA1','NSTRA2')
     CALL acvd17(ifpre,elset%head,nnode)

   CASE ('ACTUAL')
     CALL actu17(elset%head,ngaus)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump17 ( elset )

   CASE ('GAUSSC')
     CALL gaus17(elset%head,coord,istop,coora,elset%gauss, &
                 elset%plstr,elset%angdf,ntype,ngaus,nnode)

   CASE ('LOADPL')
     CALL load17 (igrav,loadv(:,:,iload),gv,gravy,elset%head,ngaus,nnode)

   !CASE ('CLOSEF')
   !  CALL close1(nreqs,narch)

   CASE ('MASMTX')
     CALL masm17(ndime,elset%head,emass,mass,sumat,flag1,ngaus,ifpre,nnode,coord,ntype)

!   CASE ('OUTDYN')
!     CALL outd17(nelem,m(nmatno),m(nindex),gprop(1),m(ngausv),nquad,
!  &              sname)

   CASE ('OUTPUT')
     IF(flag1.OR.flag2) CALL outp17 (flag1,flag2, iwrit, &
            elset%head, nreqs, narch, elset%ngrqs, ttime, ngaus)

   CASE ('RESVPL')
     IF( elset%nocom )THEN
       CALL resv17a(elset%head, coora, resid, istop, ngaus, nnode, elset%lcar , elset%epsi )
     ELSE
       IF( inverse) THEN
         CALL resv17i (elset%head, ntype, coora, resid,  &
                      istop, ttime, coord, ngaus, nnode)
       ELSE
       CALL resv17 (elset%head, ntype, coora, resid,  &
                    istop, ttime, coord, ngaus, nnode)
       END IF
     END IF

   CASE ('WRTPOS')
     CALL mase17 (ntype, nreqs, nelem, elset%head, elset%ngrqs, narch, &
                  elset%angdf, elset%sname, ngaus ,nnode)
     elset%narch = narch

   CASE ('STIFFM')
     IF( elset%nocom )THEN
       CALL stif17a(elset%head,gstif,resid(:,iforce),coora,ngaus,nnode,elset%epsi)
     ELSE
       IF( inverse) THEN
         CALL stif17i(elset%head,coord,gstif,ustif,resid(:,iforce),ntype,coora,ngaus,nnode)
       ELSE
         CALL stif17(elset%head,coord,gstif,resid(:,iforce),ntype,coora,ngaus,nnode)
       END IF
     END IF
!
!   CASE ('SURFAC')
!       ! get surface definition from the element set
!     CALL surf12 (m(nlnods),nelem,flag2,sname)

   END SELECT
   elset => elset%next
 END DO

 RETURN
 END SUBROUTINE elem17
