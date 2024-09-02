 SUBROUTINE output(inicia, foutp, ecdis )
 !***********************************************************************
 !
 !** output routine
 !
 !***********************************************************************
 USE ctrl_db,ONLY: istep, numct, ndime, ndofn, nrotd, ndyna, &
                   neulr, nload, npoin, dtime, ttime, &
                   begtm, neq, itemp, ndoft
 USE kinc_db, ONLY : nvelr,nn,nm,velor,npsdf,nesdf,ftsdf,nsdof
 USE lispa0
 USE npo_db, ONLY : label,coord,coora,ifpre,resid,fcont,force,euler,scale, &
                    veloc,accel,loadv,naeul,psia,tempe,cpx
 USE solv_db, ONLY : lambd
 USE outp_db
 USE c_input,ONLY: openfi
 IMPLICIT NONE

   !--- Dummy variables
   LOGICAL,INTENT(IN):: foutp, & ! .TRUE. write compulsory
                        inicia   ! .TRUE. first call, write initial values
   INTEGER(kind=4),INTENT(IN):: ecdis  !internal equation number controlling output

   !--- Local varibles
   INTEGER(kind=4):: ipoin,ireq,idofn,ieq,iq,ib,ie,i,j,k,nv1
   REAL(kind=8):: cputm,auxil,value,angls(3),rm(3,3),fac,xx(ndime,2)
   LOGICAL:: b1,b2,br
   ! auxiliar arrays to write down information
   REAL(kind=8),POINTER,SAVE:: disp1(:,:),velo1(:,:),acce1(:,:), disp3(:,:)
   TYPE (slave_list), POINTER :: sl_d
   CHARACTER (len=130) :: react

   INTERFACE
     INCLUDE 'angeul.h'
     INCLUDE 'elemnt.h'
     INCLUDE 'contac.h'
   END INTERFACE


   ! first determine if OUTPUT must be done at this step

 b1 = .FALSE.
 b2 = .FALSE.
 CALL timuse(cputm)
 cputm = cputm-cpui
 IF(cputm < 0) cputm = cputm + 86400d0

 IF(ndyna == 0) THEN
   value = lambd
 ELSE
   value = ttime
 END IF

 IF(MOD(istep,noutd) == 0 .OR. foutp .OR. inicia) THEN
   b1 = .TRUE.

   !***  OUTPUT FOR SELECTED VALUES

   !         D I S P L A C E M E N T S
   IF(nreqd > 0) THEN
     ALLOCATE ( disp1(ndofn,nreqd) )
     DO ireq = 1,nreqd
       ipoin = nprqd(ireq)
       DO idofn = 1,ndime
         disp1(idofn,ireq) = (coora(idofn,ipoin)-coord(idofn,ipoin))
       END DO
       ! for local systems
       IF(neulr > 0) THEN
         IF(ndime == 2) THEN                   !2-D problems (1 DOF)
           disp1(nrotd,ireq) = euler(1,ipoin)   !present angle
         ELSE                                  !3-D problems (3 DOFs)
           rm = RESHAPE( euler(1:9,ipoin),(/3,3/))          !rotation matrix
           angls = 0d0                    !angles
           CALL angeul(rm,angls,.TRUE.)   !returns Euler angles (in rads)
           disp1(ndime+1:nrotd,ireq) =   angls(1:ndofn-ndime)
         END IF
         IF( ndofn == 4 )disp1(4,ireq) = psia(1,ipoin)
         IF( ndofn == 8 )disp1(7:8,ireq) = psia(:,ipoin)
       END IF
     END DO
     WRITE(11,ERR=9999) value,((disp1(i,ireq),i=1,ndofn),ireq=1,nreqd)
     DEALLOCATE ( disp1 )
   END IF
   !    V E L O C I T I E S
   IF(nreqv > 0) THEN
     nv1 = nvelr + 1
     ALLOCATE ( velo1(nrotd,nreqv) )
     DO ireq = 1,nreqv
       ipoin = nprqv(ireq)
       DO idofn = 1,nrotd
         ieq = ifpre(idofn,ipoin)
         IF(ieq > 0 )THEN                        !active DOFs
           velo1(idofn,ireq) = veloc(ieq)
         ELSE IF(ieq < -nn)THEN                  !fixed DOFs
           velo1(idofn,ireq) = velor(-ieq-nn,nv1)
         ELSE IF(ieq < 0 .AND. ieq > -nn)THEN    !slave DOFs
           ib = npsdf(-ieq)
           ie = npsdf(-ieq+1)-1
           auxil = 0d0
           DO i = ib,ie
             iq = nesdf(i)
             IF(iq > 0)THEN
               auxil = auxil + veloc(iq)*ftsdf(i)
             ELSE IF(iq < -nn)THEN
               auxil = auxil + velor(-iq-nn,nv1)*ftsdf(i)
             END IF
           END DO
           velo1(idofn,ireq) = auxil
         ELSE
           velo1(idofn,ireq) = 0d0               !non-existent DOFs
         END IF
       END DO
     END DO
     WRITE(19,ERR=9999) value,((velo1(i,ireq),i=1,nrotd),ireq=1,nreqv)
     DEALLOCATE ( velo1 )
   END IF
   !   A C C E L E R A T I O N S
   IF(nreqa > 0) THEN
     ALLOCATE ( acce1(nrotd,nreqa) )
     DO ireq = 1,nreqa
       ipoin = nprqa(ireq)
       DO idofn = 1,nrotd
         ieq = ifpre(idofn,ipoin)
         IF(ieq > 0 )THEN                        !active DOFs
           acce1(idofn,ireq) = accel(ieq)
         ELSE
           acce1(idofn,ireq) = 0d0               !otherwise 0
         END IF
       END DO
     END DO
     WRITE(20,ERR=9999) value,((acce1(i,ireq),i=1,nrotd),ireq=1,nreqa)
     DEALLOCATE ( acce1 )
   END IF
   !  N O D A L   E Q U I V A L E N T   F O R C E S
   IF(nreql > 0) THEN
     res = resid(:,nprql(1:nreql)) !direct contributions
     !search for slave dependencies
     IF( nsdof > 0 )THEN
       sl_d => sl_head
       DO ipoin=1,nreql
         DO idofn=1,ndofn
           DO i=1,sl_d%nvalues
             j  = sl_d%deps(1,i)
             k  = sl_d%deps(2,i)
             ib = sl_d%deps(3,i)
             res(idofn,ipoin) = res(idofn,ipoin) + ftsdf(ib)*resid(j,k)
           END DO
           sl_d => sl_d%next
         END DO
       END DO

     END IF
     WRITE(12,ERR=9999) value,res(1:ndofn,1:nreql)
   END IF
   !  N O D A L   C O N T A C T   F O R C E S
   IF(numct > 0 .AND. nreqc > 0) WRITE(14,ERR=9999) value,fcont(1:ndime,nprqc(1:nreqc))
   !  P A I R S   C O N T A C T   F O R C E S
   IF(numct > 0) CALL contac('OUTDY1', ttime=value, dtime = dtime)
   !!  E N E R G Y   V A L U E S
   !IF(iener /=  0) THEN
   !  cinet = 0d0
   !  IF(nvelr > 0)THEN
   !    DO ieq = 1,neq
   !      cinet = cinet + ymass(ieq)*veloc(ieq)**2
   !    END DO
   !  END IF
   !  poter = 0d0
   !  enkint=0.
   !  DO ipoin=1,npoin
   !    DO idofn=1,ndime
   !      ieq = ifpre(idofn,ipoin)
   !      IF(ieq > 0)THEN
   !        tdisp = coora(idofn,ipoin) - coord(idofn,ipoin)
   !        IF(nload > 0)poter = poter +tdisp*force(ieq,nload+1)
   !        enkint= enkint+ymass(ieq)*veloc(ieq)**2
   !      END IF
   !    END DO
   !  END DO
   !  IF(neulr > 0)THEN
   !    enkinr = cinet - enkint
   !  ELSE
   !    enkinr = 0d0
   !    cinet  = enkint
   !  END IF
   !  velcr  = SQRT(enkint/sumat)
   !  enkint = 0.5*enkint
   !  enkinr = 0.5*enkinr
   !  cinet  = 0.5*cinet
   !  fract = 0.
   !  IF (energ(2)+cinet-energ(3).ne.0.) fract=1d0-energ(1)/(energ(2)+cinet-energ(3))
   !  WRITE(15,ERR=9999) value,cinet,velcr,dtime,enkint,enkinr,   &
   !    energ(1),energ(2),fract,cinet-energ(3),energ(5),poter
   !  WRITE(55,"(' Time=',e12.4,' Kin. E=',e12.5,' Pot. E=',e12.5)",ERR=9999) ttime,cinet,poter
   !END IF
   !! volume & pressure of volume dependent follower loads
   !CALL wrtfl1 (ttime)
 END IF

 !*** GLOBAL OUTPUT

 IF((MOD(istep,noutp) == 0 .AND. istep > 0 ) .OR. foutp ) THEN
   b2 = .TRUE.
   !message in .RSN file, if IWRIT == 1, results are printed to ASCII file
   WRITE(lures,"(5x,//'Results are reported for postprocess at:',/, &
     10x,'Istep= ',i8,5x,'Ttime= ',e15.7,/)",ERR=9999) istep,value
   WRITE(17,ERR=9999) istep,value,postype              !heading
   !** Total Displacements
   IF( .NOT.ASSOCIATED (cpx) )THEN
     DO ipoin = 1,npoin
       WRITE(17,ERR=9999) (coora(1:ndime,ipoin) - coord(1:ndime,ipoin))
     END DO
   ELSE
     DO ipoin=1,npoin
       IF( cpx(1,ipoin) == 0 )THEN
         WRITE(17,ERR=9999) (coora(1:ndime,ipoin) - coord(1:ndime,ipoin))
       ELSE ! (+/-1/2/3)
         IF( ABS(cpx(1,ipoin)) == 1  )THEN  ! .OR. ABS(cpx(1,ipoin)) == 2
           xx(:,1) = (2d0*coord(:,ipoin) + coord(:,cpx(2,ipoin))+coord(:,cpx(3,ipoin)))/4d0
           xx(:,2) = (2d0*coora(:,ipoin) + coora(:,cpx(2,ipoin))+coora(:,cpx(3,ipoin)))/4d0
         ELSE IF( ABS(cpx(1,ipoin)) == 2 )THEN
           xx(:,1) = (2d0*coord(:,ipoin) + 3d0*(coord(:,cpx(2,ipoin))+coord(:,cpx(3,ipoin))))/8d0
           xx(:,2) = (2d0*coora(:,ipoin) + 3d0*(coora(:,cpx(2,ipoin))+coora(:,cpx(3,ipoin))))/8d0
         ELSE IF( ABS(cpx(1,ipoin)) == 3 )THEN !use a couple of diagonal nodes
           xx(:,1) = (2d0*coord(:,ipoin) + coord(:,cpx(2,ipoin))+coord(:,cpx(3,ipoin)))/4d0
           xx(:,2) = (2d0*coora(:,ipoin) + coora(:,cpx(2,ipoin))+coora(:,cpx(3,ipoin)))/4d0
         END IF
         WRITE(17,ERR=9999) xx(1:ndime,2)-xx(1:ndime,1)
       END IF
     END DO
   END IF
   ! for local systems
   IF(neulr > 0) THEN
     fac = 45d0/ATAN(1d0)
     IF(ndime == 2) THEN                   !2-D problems (1 DOF)
       DO ipoin = 1,npoin
         WRITE(17,ERR=9999) euler(1,ipoin)*fac !present angle
       END DO
     ELSE                                  !3-D problems (3 DOFs)
       DO ipoin = 1,npoin
         rm = RESHAPE( euler(1:9,ipoin),(/3,3/))          !rotation matrix
         angls = 0d0                  !angles
         CALL angeul(rm,angls,.TRUE.) !returns Euler angles (in rads)
         WRITE (17,ERR=9999) angls
       END DO
     END IF
     IF( ndofn > nrotd )THEN
       DO ipoin = 1,npoin
         WRITE (17,ERR=9999) psia(:,ipoin)
       END DO
     END IF
   END IF
   IF(iwrit == 1) THEN   !write to ASCII file (.RSN)
     WRITE(lures,"(//5x,'Displacements at time step ',i10,5x,'Time ',e20.11/,    &
    &                5x,'Nnode',3x,'X-disp',6x,'Y-disp',6x,'Z-disp'/)",ERR=9999) &
    &                istep,value
     DO ipoin = 1,npoin  !displacements
       WRITE(lures,902,ERR=9999) label(ipoin),(coora(1:ndime,ipoin)-coord(1:ndime,ipoin))
     END DO
     IF(neulr > 0) THEN  !for local systems
       IF(ndime == 2) THEN             !2-D problems
         WRITE(lures,930,ERR=9999)
         DO i=1,npoin
           WRITE(lures,933,ERR=9999) label(i),euler(1,i)*fac
         END DO
       ELSE                            !3-D problemas
         WRITE(lures,931,ERR=9999)
         DO ipoin = 1,npoin
           IF(.NOT.naeul(ipoin) )CYCLE
           rm = RESHAPE( euler(1:9,ipoin),(/3,3/))          !rotation matrix
           angls = 0d0                  !angles
           CALL angeul(rm,angls,.TRUE.) !returns Euler angles (in rads)
           WRITE(lures,932,ERR=9999) label(ipoin),angls*fac
         END DO
       END IF
       IF( ndofn > nrotd  )THEN
         WRITE(lures,941,ERR=9999)
         DO ipoin = 1,npoin
           IF(ifpre(ndofn,ipoin) /= 0)WRITE(lures,902,ERR=9999) label(ipoin),psia(:,ipoin)
         END DO
       END IF
     END IF
   END IF
   IF(ndyna > 0)THEN
     !**  velocity vector
     nv1 = nvelr + 1
     IF(.NOT.ASSOCIATED(disp3)) ALLOCATE ( disp3(nrotd,npoin) )
     DO ipoin = 1,npoin
       DO idofn = 1,nrotd
         ieq = ifpre(idofn,ipoin)
         SELECT CASE (ieq)
         CASE (1:)                                 !active dof
           disp3(idofn,ipoin) = veloc(ieq)
         CASE (0)                                  !null dof
           disp3(idofn,ipoin) = 0.
         CASE ( -nm:-1)                            !slave DOFs
           ib = npsdf(-ieq)
           ie = npsdf(-ieq+1)-1
           auxil = 0d0
           DO i = ib,ie
             iq = nesdf(i)
             IF(iq > 0)THEN
               auxil = auxil + veloc(iq)*ftsdf(i)
             ELSE IF(iq < -nn)THEN
               auxil = auxil + velor(-iq-nn,nv1)*ftsdf(i)
             END IF
           END DO
           disp3(idofn,ipoin) = auxil
         CASE ( :-nn)                           !fixed DOFs
           disp3(idofn,ipoin) = velor(-ieq-nn,nv1)
         END SELECT
       END DO
     END DO

     DO ipoin=1,npoin
        WRITE(17,ERR=9999) disp3(1:nrotd,ipoin)
     END DO

     IF(iwrit == 1) THEN  !write to ASCII file
       WRITE(lures,950,ERR=9999) istep,value
       IF(ndime == 2)  WRITE(lures,952,ERR=9999)
       IF(ndime == 3)  WRITE(lures,953,ERR=9999)
       DO ipoin = 1,npoin
         WRITE(lures,902,ERR=9999) label(ipoin),disp3(1:nrotd,ipoin)
       END DO
     END IF
     !**  aceleration vector
     DO ipoin = 1,npoin   ! rearrange first
       DO idofn = 1,nrotd
         ieq = ifpre(idofn,ipoin)
         IF(ieq > 0) THEN
           disp3(idofn,ipoin) = accel(ieq)
         ELSE
           disp3(idofn,ipoin) = 0d0
         END IF
       END DO
       WRITE(17,ERR=9999) disp3(1:nrotd,ipoin)
     END DO
     IF(iwrit == 1) THEN  ! write to ASCII file
       WRITE(lures,960,ERR=9999) istep,value
       IF(ndime == 2)  WRITE(lures,962,ERR=9999)
       IF(ndime == 3)  WRITE(lures,963,ERR=9999)
       DO ipoin=1,npoin
         WRITE(lures,902,ERR=9999) label(ipoin),disp3(1:nrotd,ipoin)
       END DO
     END IF
     DEALLOCATE( disp3 ) ! release memory
   END IF
   !** CONTACT values (gaps, pressures, friction-work)
   IF (numct > 0) CALL contac('OUTDY2',ttime=value,dtime=dtime)
   !** temperatures at  nodal points
   IF (itemp) THEN
     DO ipoin=1,npoin
       WRITE(17,ERR=9999) tempe(:,ipoin)
     END DO
     !tstra = (ttime - begtm)*tscal  ! time from the strategy start
     !WRITE(17,ERR=9999) (bqgen(ipoin)/tstra,ipoin=1,npoin) ! average heat generation rate
     ! the above is for wear algorithm
     IF(iwrit==1) THEN
       WRITE(3,961,ERR=9999) istep,ttime
       DO ipoin=1,npoin
         WRITE(3,902,ERR=9999) label(ipoin),(tempe(:,ipoin))
       END DO
     END IF
   END IF
 END IF

 !*** elemental variables (Gauss-points)
 IF(iwrit == 1 .AND. b2) WRITE(lures,"(/,10x,'Stresses ',/)",ERR=9999)
 IF(b1 .OR. b2) THEN
   CALL elemnt ('OUTPUT', dtime=dtime, ttime=value, flag1=b1, flag2=b2)
   IF( ireac > 0 .AND. b2)THEN
     WRITE(lures,904,ERR=9999)
     DO ipoin=1,npoin
       WRITE(react,"(130(' '))")
       WRITE(react(1:8),"(i8)")label(ipoin)
       i = 9
       br = .FALSE.
       DO idofn=1,ndofn
         IF(ifpre(idofn,ipoin) < -nn)THEN ! next line is not rigth
           WRITE(react(i:i+14),"(E15.7)")resid(idofn,ipoin) !- loadv(idofn,ipoin,1)
           br = .TRUE.
         ELSE
           WRITE(react(i:i+14),"(15(' '))")
         END IF
         i = i+15
       END DO
       IF(br) WRITE(lures,"(a130)",ERR=9999)react
     END DO
   END IF
   IF(iwrit == 1) CALL flushf(lures)
   IF(b2) THEN
     CALL flushf(16)
     CALL flushf(17)
     WRITE(*,"(5x,'Results have been reported for postprocess')")  !screen
   END IF
   !IF (b2) CALL openfi(0,flag=0)   !erase GiD files
 END IF

 RETURN
   902 FORMAT(5x,i7,6e23.13)
   904 FORMAT(/,10x,'Nodal Reactions ',/)
   930 FORMAT(//'  Nodal angles '/)
   931 FORMAT(//'  Nodal cartesyan systems (Euler angles in rads) '/)
   932 FORMAT(1x,i10,1x,3f18.12)
   933 FORMAT(i7,3x,e15.7)
   941 FORMAT(//'  Nodal ZigZag amplitudes '/)
   950 FORMAT(//5x,'Velocity at time step ',i10,5x,'Time ',e20.11/)
   952 FORMAT(  5x,'Nnode',3x,'X-vel',7x,'Y-vel',7x,'Omega'/)
   953 FORMAT(  5x,'Nnode',3x,'X-vel',7x,'Y-vel',7x,'Z-vel',6x, &
      &            'Omega-1',5x,'Omega-2',5x,'Omega-3'/)
   960 FORMAT(//5x,'Acceleration at time step ',i10,5x,'Time ',e20.11,/)
   961 FORMAT(//5X,'Temperatures at time step ',I10,5X,'Time ',e20.11,/ &
      &         5X,'Nnode',3X,'Temperatures')
   962 FORMAT(  5x,'Nnode',3x,'X-accel',6x,'Y-accel',6x,'Alpha',/)
   963 FORMAT(  5x,'nnode',3x,'X-accel',6x,'Y-accel',6x,'Z-accel',6x, &
      &            'Alpha-1',6x,'Alpha-2',6x,'Alpha-3'/)
   965 FORMAT(//5x,'Hydrostatic pressure at time step ',i10, &
      &         5x,'Time ',e20.11/)
   966 FORMAT(i10,e13.5)
  9999 CALL runen2('')

 END SUBROUTINE output
