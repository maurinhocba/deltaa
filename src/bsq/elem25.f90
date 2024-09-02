 SUBROUTINE elem25(task ,ttime,iload,iforce,                       &
                   coora,sumat,mass,emass,resid,gstif,             &
                   istop,flag1,flag2,iwrit)

 USE loa_db, ONLY : igrav,gravy,gv
 USE npo_db, ONLY : label,coord,iffix,ifpre,coort,coorb,ifact,loadv
 USE ctrl_db, ONLY : npoin,top,bottom

 USE ele25_db
 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN) :: task
 INTEGER (kind=4), INTENT(IN), OPTIONAL :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL :: ttime,coora(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL :: mass(:),resid(:,:),gstif(:),     &
                                  emass(:,:)
 LOGICAL, OPTIONAL ::  flag1,flag2

 INTEGER (kind=4),PARAMETER :: ndiem=3, nnode=4, nstre=3

 TYPE (ele25_set), POINTER :: elset, anter
 INTEGER (kind=4) nelem, nreqs, narch
 LOGICAL :: logst
 CHARACTER (len=mnam) :: sname

 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL comm25 (1, nelem, nreqs, narch, sname, elset, logst)

   SELECT CASE (TRIM(task))

   CASE ('NEW   ','NSTRA1','NSTRA2')
     CALL acvd25(elset%head,elset%lside, &
                 nelem,ifpre,elset%nbs,elset%bhead)

   CASE ('ACTUAL')
     CALL actu25(elset%head)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump25 ( elset )

   CASE ('GAUSSC')
   IF( ASSOCIATED(elset%stint) )DEALLOCATE(elset%stint)
     ALLOCATE(elset%stint(10,nelem))
     elset%stint = 0d0
     CALL gaus25(elset%head,coord,iffix,istop,elset%gauss, &
          elset%angdf,elset%nbs,elset%bhead,nelem,elset%locax)

   CASE ('LOADPL')
     CALL load25 (igrav, loadv(:,:,iload), gv, gravy, elset%head)


   !CASE ('CLOSEF')
   !  CALL close1(nreqs,narch)

   CASE ('MASMTX')
     CALL masm25(ifpre,elset%head,emass,mass,sumat,flag1)

!   CASE ('OUTDYN')
!     CALL outd25(nelem,m(nmatno),m(nindex),gprop(1),m(ngausv),nquad,
!  &              sname)

   CASE ('OUTPUT')
      IF(flag1.OR.flag2) CALL outp25 (flag1,flag2,elset%logst, iwrit, &
   &           elset%head, nreqs, narch, elset%ngrqs, ttime, elset%stint)

   CASE ('RESVPL')
     CALL resv25 (elset%head, iffix, coora, resid, logst, istop, ttime, &
                  bottom, top, coorb, coort, ifact, elset%nbs, elset%bhead, &
                  elset%stabs, elset%stabb, elset%stint)

   CASE ('WRTPOS')
     CALL mase25 (nreqs, nelem, elset%head, elset%ngrqs, narch, &
   &              elset%angdf, elset%sname )
     elset%narch = narch

   CASE ('STIFFM')
     CALL stif25(elset%head,coora,gstif,resid(:,iforce),logst, elset%stabs, &
                 elset%stabb, elset%stint, iffix)
!
!   CASE ('SURFAC')
!       ! get surface definition from the element set
!     CALL surf12 (m(nlnods),nelem,flag2,sname)

   END SELECT
   elset => elset%next
 END DO

 RETURN
 END SUBROUTINE elem25
