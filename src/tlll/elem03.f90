 SUBROUTINE elem03(task ,ttime,iload,iforce,                       &
                   coora,locsy,sumat,mass,emass,resid,gstif,       &
                   istop,flag1,flag2,iwrit)

 USE ctrl_db, ONLY: ndime, npoin, top, bottom
 USE ele03_db
 USE npo_db, ONLY : ifpre, coord, euler, coort, coorb, ifact, loadv, eule0
 USE loa_db, ONLY : gravy, gv, igrav
 IMPLICIT NONE
 CHARACTER(len=6), INTENT(IN) :: task
 INTEGER (kind=4), INTENT(IN), OPTIONAL  :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime,coora(:,:),locsy(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),resid(:,:),gstif(:),     &
                                  emass(:,:)
 LOGICAL, OPTIONAL  ::  flag1,flag2

 REAL (kind=8),SAVE :: energ(8)
 TYPE (ele03_set), POINTER :: elset, anter
 INTEGER (kind=4) nelem, nreqs, narch
 CHARACTER (len=30) :: sname


 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head
 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL comm03 (1, nelem, nreqs, narch, sname, elset )


   SELECT CASE (task)

   CASE ('NEW   ','NSTRA1','NSTRA2')
     CALL acvd03(ifpre,elset%head,elset%quad,elset%lside,nelem,elset%nnode,elset%zigzag)

   CASE('ACTUAL')
     CALL actu03(elset%head)

   !CASE ('CLOSEF')
   !
   !  CALL close1(nreqs,narch)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump03 ( elset )

   CASE ('GAUSSC')
     CALL gaus03(coord,eule0,istop,  &
                 elset%head,elset%gauss,elset%angdf,elset%locax,elset%quad,elset%zigzag)

   CASE ('LOADPL')
     CALL load03 (igrav, loadv(:,:,iload), gv, gravy, elset%head)


   CASE ('MASMTX')

     CALL masm03(elset%head,coord,emass,mass,sumat,flag1,ifpre,elset%zigzag)

   CASE ('OUTPUT')

     IF(flag1 .OR. flag2) THEN     !flag1 = b1    flag2 = b2   argm3 = ttime
       CALL outp03(flag1,flag2,nreqs,elset%head, &
                   narch,iwrit,elset%ngrqs,ttime,elset%nstre)
       !IF(flag1) WRITE(55,'(8e12.4)')energ/2d0   !print components
     END IF

   CASE ('WRTPOS')

     CALL mase03 (nreqs, nelem, elset%head, elset%ngrqs, narch, &
                  elset%angdf, elset%sname, elset%nstre )
     elset%narch = narch

   CASE ('RESVPL','INTERN')
     energ = 0d0
     CALL resv03(elset%head,coora,locsy,resid,energ,istop, &
                 bottom,top,coorb,coort,ifact,ttime,elset%stabq,elset%quad,elset%nstre,elset%zigzag)

   CASE('STIFFM')
     CALL stif03(elset%head, coora, locsy, gstif, resid(:,iforce),elset%stabq, &
                 elset%quad,elset%nnode,elset%nstre,elset%zigzag)

   !CASE ('NODSET')      !compute nodes set from element set
   !  IF( flag2 )EXIT
   !  flag2 = TRIM(elsnam) == TRIM(sname)
   !  IF (flag2) THEN
   !    ! extracts nodes used in the discretization to nodset
   !    CALL nods03(nelem, elset%head, label)
   !    EXIT
   !  END IF

   !          ivect not defined
   !CASE ('SECDAT')
   !  IF( flag2 )EXIT
   !  flag2 = TRIM(elsnam) == TRIM(sname)
   !  IF (flag2) THEN
   !    !  extract SECTIONs used in the set
   !    CALL secd03 (elset%head,elset%nelem,ivect)
   !    EXIT
   !  END IF

   !CASE ('EXPORT')
   !  IF( flag2 )EXIT
   !  flag2 = TRIM(elsnam) == TRIM(sname)
   !  IF (flag2) THEN
   !    ! export set data
   !    CALL expo03 (elset,flag1,istop)
   !    EXIT
   !  END IF

   END SELECT
   elset => elset%next
 END DO
 RETURN
 END SUBROUTINE elem03
