 SUBROUTINE elem19(task ,elsnam, ttime,iload,iforce,               &
&                  coora,sumat,mass,emass,resid,gstif,             &
&                  istop,flag1,flag2,iwrit)

 USE loa_db, ONLY : igrav,gravy,gv
 USE ctrl_db, ONLY : npoin,ntype
 USE npo_db, ONLY : label,coord,coorc,loadv,ifpre,auxi
 USE ele19_db
 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN) :: task
 CHARACTER (len=*), OPTIONAL :: elsnam
 INTEGER (kind=4), INTENT(IN), OPTIONAL  :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime,coora(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),resid(:,:),gstif(:),     &
                                 emass(:,:)
 LOGICAL, OPTIONAL  ::  flag1,flag2

 INTEGER (kind=4), PARAMETER :: ndime=2, nstre=3

 TYPE (ele19_set), POINTER :: elset, anter
 INTEGER (kind=4) nelem, nreqs, narch
 CHARACTER (len=mnam) :: sname


 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL comm19 (1,nelem,nreqs,narch,sname,elset)

   SELECT CASE (TRIM(task))

   CASE ('NEW   ','NSTRA1','NSTRA2')
       CALL acvd19(ifpre, elset%head)

   CASE ('ACTUAL')
     CALL actu19(elset%head)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump19 ( elset )

   CASE ('GAUSSC')
     CALL gaus19(elset%head,coord,elset%gauss,elset%angdf,ntype,elset%ver)

   CASE ('LOADPL')
     CALL load19 (igrav,loadv(:,:,iload),gv,gravy,elset%head,elset%ver)

   !CASE ('CLOSEF')
   !  CALL close1(nreqs,narch)

   CASE ('MASMTX')
     CALL masm19(ndime,elset%head,emass,mass,sumat,flag1,elset%ver)

   CASE ('OUTPUT')
     IF(flag1.OR.flag2) CALL outp19 (flag1,flag2, iwrit, &
            elset%head, nreqs, narch, elset%ngrqs, ttime)

   CASE ('RESVPL','INTERN')
     CALL resv19(nelem,elset%head,ntype,coora,resid,istop,ttime,elset%angdf,coorc,elset%eulrf)

   CASE ('WRTPOS')
     CALL mase19 (ntype, nreqs, nelem, elset%head, elset%ngrqs, narch, &
                  elset%angdf, elset%sname, elset%ver )
     elset%narch = narch

   CASE ('STIFFM')
     CALL stif19(elset%head,gstif,resid(:,iforce),ntype,coora,elset%eulrf)

   CASE ('NODSET')      !compute nodes set from element set
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! extracts nodes used in the discretization to nodset
       CALL nods19( elset%head, label)
       EXIT
     END IF

   CASE ('SECDAT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       !  extract SECTIONs used in the set
       CALL secd19 (elset%head,elset%nelem,auxi)
       EXIT
     END IF

   CASE ('EXPORT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! export set data
       CALL expo19 (elset,flag1,istop)
       EXIT
     END IF

!
!   CASE ('SURFAC')
!       ! get surface definition from the element set
!     CALL surf12 (m(nlnods),nelem,flag2,sname)

   END SELECT
   elset => elset%next
 END DO

 RETURN
 END SUBROUTINE elem19
