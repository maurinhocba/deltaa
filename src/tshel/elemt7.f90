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
 INTEGER (kind=4) nelem, nreqs, narch, ansmm, nnass
 CHARACTER (len=30) :: sname


 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL commv7 (1, nelem, nreqs, narch, ansmm, nnass, sname, elset )


   SELECT CASE (task)

   CASE ('NEW   ','NSTRA1','NSTRA2')
     CALL acvdf7(ifpre,elset%head,elset%zigzag)

   CASE('ACTUAL')
     CALL actua7(elset%head)

   !CASE ('CLOSEF')
   !
   !  CALL close1(nreqs,narch)

   CASE ('DUMPIN')      !dumps variables for restart
   !  CALL dump07 ( elset )

   CASE ('GAUSSC')
     CALL gauss7(coord,eule0,istop,ansmm,nnass,elset%locax,elset%posgp,elset%shape,elset%ap1,  &
                 elset%omat,elset%head,elset%gauss,elset%angdf,elset%zigzag)

   CASE ('LOADPL')
     CALL loadp7 (igrav, loadv(:,:,iload), gv, gravy, elset%head, elset%shape)


   CASE ('MASMTX')

     CALL masmt7(elset%head,emass,mass,sumat,flag1,ifpre,5,elset%shape)   !ndofe = 5 flag1=lumpd

   CASE ('OUTPUT')

     IF(flag1 .OR. flag2) THEN     !flag1 = b1    flag2 = b2   argm3 = ttime
       CALL outpu7(flag1,flag2,nreqs,elset%head, &
                   narch,iwrit,elset%ngrqs,ttime,elset%nstre)
       !IF(flag1) WRITE(55,'(8e12.4)')energ/2d0   !print components
     END IF

   CASE ('WRTPOS')

     CALL masel7 (nreqs, nelem, elset%head, elset%ngrqs, narch, &
                  elset%angdf, elset%sname, elset%nstre, elset%gpint )
     elset%narch = narch

   CASE ('RESVPL')
     energ = 0d0
     CALL resv07(elset%head,coora,locsy,resid,energ,istop,bottom,top,coorb,coort,ifact, &
                 ttime,ansmm,nnass,elset%nstre,elset%zigzag,elset%shape,elset%ap1,elset%omat)

   CASE('STIFFM')
     CALL stiff7(elset%head, coora, locsy, gstif, resid(:,iforce), ansmm, nnass, elset%nstre, &
                 elset%zigzag,elset%posgp, elset%shape, elset%omat, elset%ap1 )

   END SELECT
   elset => elset%next
 END DO
 RETURN
 END SUBROUTINE elemt7
