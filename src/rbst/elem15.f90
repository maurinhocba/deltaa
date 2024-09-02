 SUBROUTINE elem15(task ,ttime,iload,iforce,                       &
                   coora,sumat,mass,emass,resid,gstif,             &
                   istop,flag1,flag2,iwrit)

 USE loa_db, ONLY : igrav,gravy,gv
 USE ctrl_db, ONLY : npoin,top,bottom
 USE npo_db, ONLY : label,oldlb,coord,loadv,iffix,coort,coorb,ifact,ifpre
 USE ele15_db
 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN) :: task
 INTEGER (kind=4), INTENT(IN), OPTIONAL :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime,coora(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),resid(:,:),gstif(:),     &
                                  emass(:,:)
 LOGICAL, OPTIONAL  ::  flag1,flag2

 INTEGER (kind=4), PARAMETER :: ndime=3, nnode=3, nstre=3

 TYPE (ele15_set), POINTER :: elset, anter
 INTEGER (kind=4) nelem, nreqs, narch
 LOGICAL :: logst
 CHARACTER (len=mnam) :: sname


 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL comm15 (1, nelem, nreqs, narch, sname, elset, logst)

   SELECT CASE (TRIM(task))

   CASE ('NEW   ','NSTRA1','NSTRA2')
     CALL acvd15(ifpre,elset%head,elset%lside, &
                 nelem,elset%nrf,elset%rhead)

   CASE ('ACTUAL')
     CALL actu15(elset%head)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump15 ( elset )

   CASE ('GAUSSC')
     IF( ASSOCIATED(elset%stint) )DEALLOCATE(elset%stint)
     ALLOCATE(elset%stint(10,nelem))
     elset%stint = 0d0
     CALL gaus15(elset%head,coord,iffix,istop,elset%gauss,  &
             elset%angdf,elset%nrf,elset%rhead,nelem, &
             elset%shear,elset%shears,elset%factors,elset%ninv)

   CASE ('LOADPL')
     CALL load15 (igrav, loadv(:,:,iload), gv, gravy, elset%head, elset%rhead)


   !CASE ('CLOSEF')
   !  CALL close1(nreqs,narch)

   CASE ('MASMTX')
     CALL masm15(ndime,nnode,elset%head,emass,mass,sumat,flag1,elset%nrf,elset%rhead,ifpre)

!   CASE ('OUTDYN')
!     CALL outd15(nelem,m(nmatno),m(nindex),gprop(1),m(ngausv),nquad,
!  &              sname)

   CASE ('OUTPUT')
      IF(flag1.OR.flag2) CALL outp15 (flag1,flag2,elset%logst, iwrit, &
   &           elset%head, nreqs, narch, elset%ngrqs, ttime, elset%stint, &
               elset%rhead)

   CASE ('RESVPL')
     CALL resv15 (elset%head, iffix, coora, resid, logst, istop, ttime, &
                  bottom, top, coorb, coort, ifact, elset%nrf, elset%rhead, &
                  elset%stint,elset%shear,elset%shears,elset%factors,elset%ninv)

   CASE ('WRTPOS')
     CALL mase15 (nreqs, nelem, elset%head, elset%ngrqs, narch, &
   &              elset%angdf, elset%sname, elset%nrf, elset%rhead )
     elset%narch = narch

   CASE ('STIFFM')
     CALL stif15(elset%head,coora,gstif,resid(:,iforce),logst, elset%stint, iffix)
!
!   CASE ('SURFAC')
!       ! get surface definition from the element set
!     CALL surf12 (m(nlnods),nelem,flag2,sname)

   END SELECT
   elset => elset%next
 END DO

 RETURN
 END SUBROUTINE elem15
