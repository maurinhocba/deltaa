 SUBROUTINE elem04(task ,ttime,iload,iforce,                       &
                   coora,sumat,mass,emass,resid,gstif,             &
                   istop,flag1,flag2,iwrit)

 USE ctrl_db, ONLY : npoin
 USE loa_db, ONLY : igrav,gravy,gv
 USE npo_db, ONLY : label,coord,ifpre,loadv !,oldlb
 USE ele04_db
 IMPLICIT NONE
 CHARACTER(len=6), INTENT(IN) :: task

 INTEGER (kind=4), INTENT(IN), OPTIONAL  :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime,coora(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),resid(:,:),gstif(:),     &
                                  emass(:,:)
 LOGICAL, OPTIONAL  ::  flag1,flag2


 INTEGER (kind=4) nelem, nreqs, narch
 CHARACTER (len=mnam) :: sname
 TYPE (ele04_set), POINTER :: elset, anter


 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL comm04 (1,nelem,nreqs,narch,sname,elset)

   SELECT CASE (task)

   CASE ('NEW   ','NSTRA1','NSTRA2')
     CALL acvd04(ifpre,elset%head,elset%lside,elset%linear)

   CASE ('GAUSSC')
     CALL gaus04(elset%head,coord,istop,coora, &
                 elset%gauss,elset%plstr,elset%angdf)

   CASE ('LOADPL')
     CALL load04 (igrav,loadv(:,:,iload),gv,gravy,elset%head,coord)


   !CASE ('CLOSEF')
   !  CALL close1(nreqs,narch)

   CASE ('MASMTX')
     CALL masm04(elset%head,emass,mass,sumat,flag1,ifpre)

!   CASE ('OUTDYN')
!     CALL outd04(nelem,m(nmatno),m(nindex),gprop(1),m(ngausv),nquad,
!  &              sname)

   CASE ('OUTPUT')
     IF(flag1.OR.flag2) CALL outp04 (flag1,flag2, iwrit, &
            elset%head, nreqs, narch, elset%ngrqs, ttime)

   CASE ('RESVPL')
     CALL resv04 (elset%head, coora, resid, istop, ttime)

   CASE ('WRTPOS')
     CALL mase04 (nreqs, nelem, elset%head, elset%ngrqs, narch, &
                  elset%angdf, elset%sname )
     elset%narch = narch

   CASE ('STIFFM')
     CALL stif04(elset%head,gstif,resid(:,iforce),coora)
!
!   CASE ('SURFAC')
!       ! get surface definition from the element set
!     CALL surf12 (m(nlnods),nelem,flag2,sname)

   END SELECT
   elset => elset%next
 END DO

 RETURN
 END SUBROUTINE elem04
