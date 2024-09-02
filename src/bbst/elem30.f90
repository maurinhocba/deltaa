 SUBROUTINE elem30(task , elsnam, ttime,iload,iforce,              &
                   coora,sumat,mass,emass,resid,gstif,             &
                   istop,flag1,flag2,iwrit)

 USE loa_db, ONLY : igrav,gravy,gv
 USE ctrl_db, ONLY : npoin,top,bottom
 USE npo_db, ONLY : label,coord,coort,coorb,ifact,iffix,ifpre,loadv,auxi
 USE ele30_db
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

 INTEGER (kind=4), PARAMETER :: ndime=3, nnode=3, nstre=3

 TYPE (ele30_set), POINTER :: elset, anter
 INTEGER (kind=4) nelem, nreqs, narch
 LOGICAL :: logst
 CHARACTER (len=mnam) :: sname

 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL comm30 (1, nelem, nreqs, narch, sname, elset, logst)

   SELECT CASE (TRIM(task))

   CASE ('NEW','NSTRA1','NSTRA2')
      CALL acvd30 ( ifpre, elset%head, nelem, elset%lside, npoin)

!   CASE ('ACTUAL')
!     CALL actu30(elset%head)
!
!   CASE ('DUMPIN')      !dumps variables for restart
!     CALL dump30 ( elset )
!
   CASE ('GAUSSC')
     CALL gaus30(elset%head,coord,iffix,istop, &
                 elset%gauss,elset%angdf,elset%locax)

    CASE ('LOADPL')
      CALL load30 (igrav, loadv(:,:,iload), gv, gravy, elset%head)

!   !CASE ('CLOSEF')
!   !  CALL close1(nreqs,narch)
!
    CASE ('MASMTX')
      CALL masm30(ndime,nnode,ifpre,elset%head,emass,mass,sumat,flag1)
!
!!   CASE ('OUTDYN')
!!     CALL outd30(nelem,m(nmatno),m(nindex),gprop(1),m(ngausv),nquad,
!!  &              sname)
!
    CASE ('OUTPUT')
       IF(flag1.OR.flag2) CALL outp30 (flag1,flag2,elset%logst, iwrit, &
    &           elset%head, nreqs, narch, elset%ngrqs, ttime, elset%stint)

   CASE ('RESVPL','INTERN')
     CALL resv30 (nelem, elset%head, iffix, coora, resid, logst, istop, ttime, bottom, &
                  top, coorb, coort, ifact, elset%stint)

    CASE ('WRTPOS')
      CALL mase30 (nreqs, nelem, elset%head, elset%ngrqs, narch, &
    &              elset%angdf, elset%sname, elset%locax )
      elset%narch = narch

    CASE ('STIFFM')
      CALL stif30(elset%head,coora,gstif,resid(:,iforce),logst, elset%stint, iffix)
!
!   CASE ('NODSET')      !compute nodes set from element set
!     IF( flag2 )EXIT
!     flag2 = TRIM(elsnam) == TRIM(sname)
!     IF (flag2) THEN
!       ! extracts nodes used in the discretization to nodset
!       CALL nods30(nelem, elset%head, label)
!       EXIT
!     END IF
!
!!   CASE ('INIGAU')      !read initial stresses or internal variables
!!     IF( flag2 )EXIT
!!     flag2 = TRIM(elsnam) == TRIM(sname)
!!     IF (flag2) THEN
!!       ! read initial internal variables for the element set
!!       CALL inig30 (elset%head,nelem,iwrit,elset%plstr,coora,iffix,elset%nonrg)
!!       EXIT
!!     END IF
!!
!   CASE ('SECDAT')
!     IF( flag2 )EXIT
!     flag2 = TRIM(elsnam) == TRIM(sname)
!     IF (flag2) THEN
!       !  extract SECTIONs used in the set
!       CALL secd30 (elset%head,elset%nelem,auxi)
!       EXIT
!     END IF
!
!   CASE ('EXPORT')
!     IF( flag2 )EXIT
!     flag2 = TRIM(elsnam) == TRIM(sname)
!     IF (flag2) THEN
!       ! export set data
!       CALL expo30 (elset,flag1,istop)
!       EXIT
!     END IF
!
!!
!!   CASE ('SURFAC')
!!       ! get surface definition from the element set
!!     CALL surf12 (m(nlnods),nelem,flag2,sname)

   END SELECT
   elset => elset%next
 END DO

 RETURN
 END SUBROUTINE elem30
