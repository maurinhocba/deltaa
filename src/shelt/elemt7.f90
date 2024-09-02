 SUBROUTINE elemt7(task ,ttime,iload,iforce,                       &
                   coora,locsy,sumat,mass,emass,resid,gstif,       &
                   istop,flag1,flag2,iwrit)

 USE ctrl_db, ONLY: ndime, npoin, top, bottom
 USE ele07_db
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
 TYPE (ele07_set), POINTER :: elset, anter
 INTEGER (kind=4) nelem, nreqs, narch, stype
 CHARACTER (len=30) :: sname


 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL commv7 (1, nelem, nreqs, narch, stype, sname, elset )


   SELECT CASE (task)

   CASE ('NEW   ','NSTRA1','NSTRA2')
     CALL acvdf7(ifpre,elset%head,elset%zigzag)

   CASE('ACTUAL')
     CALL actua7(elset%head)

   !CASE ('CLOSEF')
   !
   !  CALL close1(nreqs,narch)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump07 ( elset )

   CASE ('GAUSSC')
     CALL gauss7(stype,coord,eule0,istop,  &
                 elset%head,elset%gauss,elset%angdf,elset%locax,elset%zigzag)

   CASE ('LOADPL')
     CALL loadp7 (igrav, loadv(:,:,iload), gv, gravy, elset%head)


   CASE ('MASMTX')

     CALL masmt7(elset%head,coord,emass,mass,sumat,flag1,ifpre,5)   !ndofe = 5

   CASE ('OUTPUT')

     IF(flag1 .OR. flag2) THEN     !flag1 = b1    flag2 = b2   argm3 = ttime
       CALL outpu7(flag1,flag2,nreqs,elset%head, &
                   narch,iwrit,elset%ngrqs,ttime,elset%nstre)
       !IF(flag1) WRITE(55,'(8e12.4)')energ/2d0   !print components
     END IF

   CASE ('WRTPOS')

     CALL masel7 (nreqs, nelem, elset%head, elset%ngrqs, narch, &
                  elset%angdf, elset%sname, elset%nstre)
     elset%narch = narch

   CASE ('RESVPL')
     energ = 0d0
     CALL resvp7(elset%head,coora,locsy,resid,energ,istop, &
                 bottom,top,coorb,coort,ifact,ttime,stype,elset%nstre,elset%zigzag)

   CASE('STIFFM')
     CALL stiff7(elset%head, coora, locsy, gstif, resid(:,iforce), stype,elset%nstre,elset%zigzag)

   END SELECT
   elset => elset%next
 END DO
 RETURN
 END SUBROUTINE elemt7
