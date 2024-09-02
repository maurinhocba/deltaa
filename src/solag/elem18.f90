 SUBROUTINE elem18(task ,elsnam,ttime,iload,iforce,                &
                   coora,sumat,mass,emass,resid,gstif,ustif,       &
                   istop,flag1,flag2,iwrit)

 USE ctrl_db, ONLY : inverse
 USE loa_db, ONLY : igrav,gravy,gv
 USE npo_db, ONLY : label,coord,loadv,ifpre,auxi
 USE ele18_db
 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN) :: task
 CHARACTER (len=*),OPTIONAL :: elsnam
 INTEGER (kind=4), INTENT(IN), OPTIONAL  :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime,coora(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),resid(:,:),gstif(:),     &
                                  ustif(:),emass(:,:)
 LOGICAL, OPTIONAL  ::  flag1,flag2

 ! Local variables
 INTEGER (kind=4), PARAMETER :: ndime = 3
 CHARACTER (len=mnam) :: sname
 INTEGER (kind=4) :: nreqs, narch, nelem, newel, nnode, ngaus
 TYPE (ele18_set), POINTER :: elset, anter

 IF ( .NOT.ASSOCIATED (head) ) RETURN  !no element sets

 NULLIFY (anter)      !nullify pointer to previous set
 elset => head        !point to first set of elements

 DO                   !loop over the element sets
   ! recover control variables for present element set
   CALL comm18 (1,nelem,nreqs,narch,sname,elset,nnode,ngaus)

   SELECT CASE (TRIM(task)) !according to the requested task

   CASE ('NEW   ','NSTRA1','NSTRA2')  !release DOFs and other tasks
     CALL acvd18(nnode,ifpre,elset%head)

   CASE ('ACTUAL')
     CALL actu18(elset%head,ngaus)

   CASE ('GAUSSC')      !compute initial constants
     CALL gaus18(elset%head,coord,istop,nnode,elset%gauss,elset%angdf, &
                 ngaus,elset%shell,elset%gpc,elset%check,elset%locax)

   CASE ('LOADPL')      !compute gravity load vector
    CALL load18 (igrav, loadv(:,:,iload), gv, gravy, elset%head, nnode, ngaus)

   CASE ('MASMTX')
     IF( nnode == 8 )THEN
       CALL masm18(ndime,elset%head,emass,mass,sumat,flag1,nnode, ngaus,ifpre)
     ELSE
       CALL masm18q(ndime,elset%head,emass,mass,sumat,flag1,nnode,ifpre)
     END IF

   CASE ('OUTPUT')      !output variables for post-processing

     IF(flag1.OR.flag2) CALL outp18 (flag1,flag2, iwrit, &
            elset%head, nreqs, narch, elset%ngrqs, ttime, ngaus)

   CASE ('RESVPL','INTERN')      !compute internal nodal forces
     IF( inverse) THEN
       CALL resv18i(nelem, elset%head, coora, resid, istop, ttime, elset%small, nnode, coord, ngaus, &
                    elset%shell, elset%gpc,elset%bbar)
     ELSE
       CALL resv18 (nelem, elset%head, coora, resid, istop, ttime, elset%small, nnode, ngaus, &
                    elset%shell, elset%gpc,elset%bbar)
     END IF

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump18 ( elset )

   CASE ('WRTPOS')      !writes geometric variables for post-process
     CALL mase18 ( nreqs, nelem, elset%head, elset%ngrqs, narch, &
                   elset%angdf, elset%sname, nnode, ngaus, elset%locax )
     elset%narch = narch              !keep output file

   CASE ('STIFFM')
     IF( inverse) THEN
       CALL stif18i(elset%head,gstif,ustif,resid(:,iforce),coora,coord,nnode,ngaus,elset%shell,elset%gpc,elset%bbar)
     ELSE
       CALL stif18 (elset%head,gstif,resid(:,iforce),coora,nnode,ngaus,elset%shell,elset%gpc,elset%bbar)
     END IF

   CASE ('NODSET')      !compute nodes set from element set
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! extracts nodes used in the discretization to nodset
       CALL nods18(nelem,nnode, elset%head, label)
       EXIT
     END IF

   !CASE('SLNODS')
   !  IF( flag2 )EXIT
   !  flag2 = TRIM(elsnam) == TRIM(sname)
   !  IF (flag2) THEN
   !    esta = PRESENT(flag1)
   !    IF( esta ) esta = flag1
   !    ! extracts element connectivities into array IMAT
   !    CALL slno18(nelem,nnode,elset%head,lnod,esta)
   !    EXIT
   !  END IF

   CASE ('SECDAT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       !  extract SECTIONs used in the set
       CALL secd18 (elset%head,elset%nelem,auxi)
       EXIT
     END IF

   CASE ('EXPORT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! export set data
       CALL expo18 (elset,flag1,istop)
       EXIT
     END IF

   END SELECT

   IF ( ASSOCIATED (elset%next) ) THEN   !more sets to process
     elset => elset%next                 !point to next set
   ELSE
     EXIT
   END IF

 END DO

 RETURN

 END SUBROUTINE elem18
