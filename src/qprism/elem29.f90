 SUBROUTINE elem29(task , elsnam, ttime,iload,iforce,              &
                   coora,sumat,mass,emass,resid,gstif,             &
                   istop,flag1,flag2,iwrit)

 USE loa_db, ONLY : igrav,gravy,gv
 !USE ctrl_db, ONLY : npoin
 USE npo_db, ONLY : label,coord,loadv,ifpre
 USE ele29_db
 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN) :: task
 CHARACTER (len=*),OPTIONAL :: elsnam
 INTEGER (kind=4), INTENT(IN), OPTIONAL  :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime,coora(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),resid(:,:),gstif(:),     &
                                  emass(:,:)
 LOGICAL, OPTIONAL  ::  flag1,flag2

 ! Local variables
 INTEGER (kind=4), PARAMETER :: ndime = 3
 CHARACTER (len=mnam) :: sname
 INTEGER (kind=4) :: nreqs, narch, nelem, newel, ngaus, ansmm, anssh, nassp
 TYPE (ele29_set), POINTER :: elset, anter

 IF ( .NOT.ASSOCIATED (head) ) RETURN  !no element sets

 NULLIFY (anter)      !nullify pointer to previous set
 elset => head        !point to first set of elements

 DO                   !loop over the element sets
   ! recover control variables for present element set
   CALL comm29 (1,nelem,nreqs,narch,sname,elset,ngaus,ansmm,anssh,nassp)

   SELECT CASE (TRIM(task)) !according to the requested task

   CASE ('NEW   ','NSTRA1','NSTRA2')  !release DOFs and other tasks
     CALL acvd29(ifpre,elset%head,nelem)

   CASE ('ACTUAL')
   !  CALL actu29(elset%head,ngaus)

   CASE ('GAUSSC')      !compute initial constants
     IF( .NOT.elset%gauss)THEN
       CALL gaus29(elset%head,coord,istop,elset%angdf,elset%locax, &
                   elset%isgo,elset%psgo,ngaus,ansmm,anssh,nassp,elset%small)
       elset%gauss = .TRUE.
     END IF

   CASE ('LOADPL')      !compute gravity load vector
    CALL load29 (igrav, loadv(:,:,iload), gv, gravy, elset%head)

!   CASE ('CLOSEF')      !close output file
!     CALL close1(nreqs,narch)

   CASE ('MASMTX')
     CALL masm29(ndime,elset%head,emass,mass,sumat,flag1,ifpre)

   CASE ('OUTPUT')      !output variables for post-processing

     IF(flag1.OR.flag2) CALL outp29 (flag1,flag2, iwrit, &
            elset%head, nreqs, narch, elset%ngrqs, ttime, ngaus, &
            elset%isgo,elset%psgo)

   CASE ('RESVPL')      !compute internal nodal forces
     CALL resv29 (elset%head, coora, resid, istop, ttime, elset%small, ngaus,  &
                  ansmm, anssh, nassp)

   CASE ('DUMPIN')      !dumps variables for restart
     !CALL dump29 ( elset )

   CASE ('WRTPOS')      !writes geometric variables for post-process
     CALL mase29 ( nreqs, nelem, elset%head, elset%ngrqs, narch, &
                   elset%angdf, sname, elset%locax )
     elset%narch = narch              !keep output file

   CASE ('STIFFM')
      CALL stif29(elset%head, gstif, resid(:,iforce), coora, elset%small, ngaus, &
                  ansmm,anssh,nassp)

   !CASE ('SPBACK')
   !    CALL outs29 (elset, label)


!   CASE ('SEARCH')      !search if set named SNAME exists
!     IF( flag2 )EXIT
!     flag2 = elsnam == sname
!     IF (flag2) EXIT

!   CASE ('SURFAC','BOUNDA')
!     IF( flag2 )EXIT
!     flag2 = elsnam == sname
!     IF (flag2) THEN
!       ! get surface definition from the element set
!       CALL surf29 (elset%head,elset%nelem)
!       EXIT
!     END IF
!

   !CASE ('SMOOTH')      !performs smoothing for remeshing and transfer
   !  IF( flag2 )EXIT
   !  flag2 = elsnam == sname
   !  IF (flag2) THEN
   !    ! smooth gaussian variables and writes
   !    CALL smth29(coord,coora,strnd,elset%head,iffix,npoin, &
   !               m_last,r_last,r_elm_ratio,r_elm_size,r_min_size)
   !    EXIT
   !  END IF

   !CASE ('NODSET')      !compute nodes set from element set
   !  IF( flag2 )EXIT
   !  flag2 = elsnam == sname
   !  IF (flag2) THEN
   !    ! extracts nodes used in the discretization to nodset
   !    CALL nods29 (npoin, nelem, numpo, nodset, elset%head, label, l_old)
   !    EXIT
   !  END IF

   !CASE ('INIGAU')      !read initial stresses or internal variables
   !  IF( flag2 )EXIT
   !  flag2 = elsnam == sname
   !  IF (flag2) THEN
   !    ! read initial internal variables for the element set
   !    CALL inig29 (elset%head,nelem,iwrit,elset%plstr)
   !    EXIT
   !  END IF

   !CASE ('REALLO')
   !
   !  IF( elsnam == sname ) THEN
   !    ! reallocates the element set after remeshing
   !    ! exchange mesh - deallocate old mesh and put new mesh,
   !    ! interpolate all the Gaussian variables from nodes to Gauss points
   !    CALL real29 (npoin, ne_new, maxnn, l_new, elset,            &
   !                 s_new, iwrit, label, coord, r_last)
   !
   !  ELSE IF( r_last )THEN
   !    ! updates local node numbering in not remeshed sets
   !    CALL acvd29 (task, npoin, ifpre, &
   !                 elset%head, label, oldlb)
   !  END IF

   CASE ('EXPORT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       !export set data
       !CALL expo29 (elset,flag1,istop)
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

 END SUBROUTINE elem29
