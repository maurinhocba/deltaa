SUBROUTINE eqnums (actio )

  !***  calculates equation numbers,
  !***  applies boundary conditions
  !***  computes matrices profile
  !***  allocates memory for arrays
  !***  generates CURve-PArameters array

  USE lispa0, ONLY : lures
  USE ctrl_db, ONLY: ndime, ndofn, nrotd, neulr, npoin,  neq, memo, numct, nload,  &
                     maxa, top, bottom, ndyna, ndoft, itemp, inverse
  USE curv_db
  USE esets_db,ONLY : maxnn, melem, gnods, bestn, gelem, nelms, rot_free
  USE kinc_db, ONLY : nescv,ndepd,naris,ndumn,distd,nndpd, &
                      nsdof,nesdf,ftsdf,npsdf,ftsd0,maxav, &
                      nvfix,nvelr,lcvel,velor,nn,sldnd,sldvl
  USE ndp_db, ONLY : nndp
  USE nsld_db, ONLY : nnsld
  USE nes_db, ONLY : nnes
  USE npo_db, ONLY : ifpre,iffix,label,coord,coora,coorc,naeul,   &
                     coor1,euler,eule1,locs1,emass,ifact,coort,   &
                     coorb,fcont,id,resid,locsy,force,loadv,loass,&
                     veloc,velnp,accel,mass,damp,tempe,dtemp,iftmp
  USE outp_db, ONLY : iwrit,nreql
  USE solv_db, ONLY : nbuck,nsymm,displ,ddisp,disax,gvect,stiff

  IMPLICIT NONE

  CHARACTER(len=*),INTENT(IN):: actio

  ! local
  LOGICAL          :: error
  INTEGER (kind=4) :: i,j,k,n,ndf,isize,newnn
  INTEGER (kind=4), ALLOCATABLE :: dbase(:)
  REAL (kind=8) :: ang(3)


  INTERFACE
    INCLUDE 'elemnt.h'
    INCLUDE 'inrotm.h'
    INCLUDE 'fixval.h'
    INCLUDE 'ubicmx.h'
    INCLUDE 'fixtem.h'
    INCLUDE 'rdvelr.h'
  END INTERFACE


  CALL gen_curpa ( )  !generates Curve parameters array
  IF (TRIM(actio) /= 'NSTRA0') THEN   !if a new problem or important changes

    ALLOCATE ( gnods(maxnn,melem), bestn(npoin) )     !auxiliar temporary array
    gelem = 0                           !initializes
    gnods = 0

    ! deallocate arrays that need redimensioning
    IF( ASSOCIATED(resid) )DEALLOCATE(resid)
    IF( ASSOCIATED(ifpre) )DEALLOCATE(ifpre)
    IF( ASSOCIATED(iffix) )DEALLOCATE(iffix)
    IF( ASSOCIATED(emass) )DEALLOCATE(emass)
    ALLOCATE( ifpre(ndofn,npoin), emass(ndofn,npoin),resid(ndofn,npoin))
    resid = 0d0
    IF (rot_free) THEN
      ALLOCATE( iffix(npoin) )
    ELSE
      ALLOCATE( iffix(1) )
    END IF
    !  arrays associated to contac
    IF( ASSOCIATED(fcont))DEALLOCATE ( fcont )
    IF (numct > 0) THEN
      ALLOCATE( fcont(ndime,npoin) )
      IF ( top .OR. bottom ) THEN       !if bottom or top surfaces will be used
        ! get memory for auxiliar arrays (Coordinates and factors)
        IF (ASSOCIATED(ifact)) DEALLOCATE(ifact)
        ALLOCATE( ifact(npoin) )
        IF (bottom )THEN
          IF(ASSOCIATED(coorb)) DEALLOCATE(coorb)
          ALLOCATE(coorb(ndime,npoin))
        END IF
        IF (top) THEN
          IF(ASSOCIATED(coort)) DEALLOCATE(coort)
          ALLOCATE(coort(ndime,npoin))
        END IF
      END IF
    ELSE
      ALLOCATE( fcont(1,1) ) !necessary?
    END IF

    ifpre = 1        !initializes the ifpre array
    IF(itemp)THEN
      IF( ASSOCIATED(iftmp) )DEALLOCATE(iftmp)
      ALLOCATE( iftmp(ndoft,npoin) )
      iftmp = 0        !initializes the iftmp array
    END IF

    ndumn = 0        !initializes number of dummy nodes (dependant nodes in rigid bodies)
    IF( neulr )THEN
      IF( ASSOCIATED(naeul) )DEALLOCATE(naeul)
      ALLOCATE( naeul(npoin) )
      naeul = .FALSE. !initializes
    END IF

    CALL elemnt (actio)    !release DOFs associated to elements & compute NDUMN

    ! after this loop, NAEUL kwows the nodes associated to Shear deformable shells and beams
    ! generates rotation matrix for the nodal cartesian systems

    IF (TRIM(actio) /= 'NSTRA1') THEN   !if a new problem or important changes
      IF(neulr == 1) THEN      !ndime == 3
        DO n=1,npoin                !for each node
          locsy(:,n) = euler(:,n)
        END DO
      END IF
      IF(neulr == 9) THEN      !ndime == 3
        DO n=1,npoin                !for each node
          IF( .NOT. naeul(n) )CYCLE
          ang = euler(1:3,n)              !Recover the Euler angles (radians)
          CALL inrotm(ang,euler(1:9,n))   !compute rotation matrix
          locsy(:,n) = euler(:,n)
        END DO
      END IF
    END IF
    ! and master nodes in rigid bodies
    ndepd = nndp+ndumn     !total number of dependant nodes
    IF (ASSOCIATED(ftsdf)) DEALLOCATE ( ftsdf,nesdf,npsdf,nndpd,distd )
    ALLOCATE( distd(3,MAX(ndepd,1)), nndpd(3,MAX(ndepd+naris,1)) )
    distd = 0d0
    ! sliding nodes
    IF( ASSOCIATED (sldnd) )DEALLOCATE( sldnd, sldvl )
    IF( nnsld > 0 )THEN
        ALLOCATE( sldnd(3,nnsld), sldvl(2,nnsld) )
    END IF
    WRITE(lures, "(//' KINEMATIC CONDITIONS applied in this strategy',/)",ERR=9999)

    !   *** APPLY linear constraints on translational degrees of freedom
    !auxiliar pointers (divide auxiliar array DBASE in three arrays)
    ALLOCATE ( dbase(memo) )           !auxiliar temporary array
    isize = memo/13     !first  pointer 1 >  npsdf    pointers
    j = 1 + isize       !second pointer j >  nesdf    master equations
    i = j + 4*isize     !third  pointer i >  ftsdf    factors
    dbase(1) = 1
    nescv = 0   !initializes number of slave constrains
    IF(nnes > 0) CALL reasdf(nescv,dbase(1),dbase(j),dbase(i))
    !consider nodes on an 'arista'
    nsdof = nescv+1          !first position after already computed part
    IF( naris > 0) CALL arisin(naris,dbase(1),dbase(j),dbase(i),nndpd(1,ndepd+1),nsdof,isize)
    !consider dependent nodes defined as a pair or rigid bodies
    IF(ndepd > 0) CALL readpd(dbase(1),dbase(j),dbase(i),nsdof,isize)
    !!consider rotation-free dependent nodes
    !IF(nrfdp > 0) CALL rdrfdp(dbase(1),dbase(j),dbase(i),nsdof,isize)
    ! the following lines transfer the slave DOF data
    !consider sliding nodes
    IF(nnsld > 0) CALL inisld(dbase(1),dbase(j),dbase(i),nsdof,isize)
    ALLOCATE( npsdf(nsdof) )
    npsdf(:) = dbase(1:nsdof) !assign pointers already computed
    !compute Maximum total length of arrays NESDF & FTSDF
    n   = npsdf(nsdof)
    ALLOCATE ( nesdf(n), ftsdf(n) )  !reserve Maximum space
    CALL vecasi(n,dbase(i),ftsdf(1)) !assign factors already computed
    nesdf(1:n) = dbase(j:j+n-1)      !assign equations already computed
    nsdof = nsdof - 1
    IF( inverse .AND. nsdof > 0 )THEN
      ALLOCATE ( ftsd0(n) )  !reserve Maximum space
      ftsd0 = ftsdf
    END IF

    DEALLOCATE (dbase)
    !   ***    APPLY the fixed values.

    IF(rot_free) iffix = 0      !initializes boundary conditions for BST

    CALL fixval(iwrit,ndime,ndofn,ifpre,nvfix,rot_free,iffix,label)

    !   APPLY constrained velocities

    CALL rdvelr(nvelr,ndime,ndofn,npoin,iwrit,label, &
                ifpre,nvfix,lcvel,velor)

    !     renumber nodes for optimum matrix storage

    CALL renumn(npoin,newnn)

    !       *** assign number to active equations

    neq = 0                         !initializes
    DO n=1,npoin                    !for each node in the mesh
      j = bestn(n)
      DO i=1,ndofn                  !for each DOFs
        IF(ifpre(i,j) == 0) THEN    !IF an active DOF
          neq = neq+1               !increase number of active DOFs
          ifpre(i,j) = neq          !assign equation number
        ELSE IF(ifpre(i,j) == 1) THEN  !if a non existent DOF
          ifpre(i,j) = 0            !put zero
        END IF
      END DO
    END DO
    ! check if any local system must be activated (may be unnecessary)
    IF( neulr )THEN
      ndf = ndime+1
      DO i=1,npoin  !released DOF or nonzero-prescribed DOF
        IF( naeul(i) )CYCLE
        naeul(i) =  ANY(ifpre(ndf:nrotd,i) > 0) .OR.  ANY(ifpre(ndf:nrotd,i) < -nn )
      END DO
    END IF
    ! modify NESDF array to keep the real DOF
    error = .FALSE.
    DO j=1,npsdf(nsdof+1)-1      !for any kind of slave DOFs
      k = nesdf(j)               !compound value
      IF( k > 0)THEN             !excluding nodes in rigid bodies
        i = MOD(k,10)            !DOF
        n = k/10                 !node (internal)
        k = ifpre(i,n)           !associated equation
        IF( k >= 0 .OR. k < -nn)THEN   !check that it is not SLAVE-SLAVE
          nesdf(j) = k      !transfer equation (active, inactive or prescribed)
        ELSE
          WRITE(lures,"(' ERROR, slave DOF depending on slave DOF',/ &
              & ,'  NODE',i6,' DOF',i3,' Must NOT be slave')",ERR=9999)label(n),i
          WRITE(*,"(' ERROR, slave DOF depending on slave DOF')")label(n),i
          error = .TRUE.
        END IF
      END IF
    END DO
    IF( error ) CALL runen3(' INVALID SLAVE-SLAVE association')
    ! modify IFFIX for nodes prescribed on a side
    IF( rot_free )THEN        !when rotation free elements exist
      DO i=ndepd+1,ndepd+naris            !for each node of this type
        j = ABS(nndpd(1,i)) !slave node
        iffix(j) = -i  !clamped condition is assumed
        !IF( iffix(j) > 0 ) iffix(j) = -i  !keep position in array NNDPD
      END DO
    END IF

    IF (ASSOCIATED(force)) DEALLOCATE (force, loadv, loass)
    IF( nload /= 0 )ALLOCATE( force(neq+1,nload+1), loadv(ndofn,npoin,nload), loass(MAX(nload,1)) )

    CALL flushf(lures)
    IF(iwrit == 1) THEN
      IF(ndime == 2)WRITE(lures,"(//'  Information Relative to the ', &
         & 'Equations Numbers ',/,6x,'Node',7x,'X',7x,'Y',4x,'Alfa'/)",ERR=9999)
      IF(ndime == 3)WRITE(lures,"(//'  Information Relative to the ', &
         & 'Equations Numbers ',/,6x,'NODE',7x,'X',7x,'Y',7x,'Z',4x, &
         &  'Alfa    Beta   Gamma  Orig')",ERR=9999)
      DO n=1,npoin
        WRITE(lures,"(11i8)",ERR=9999) label(n),ifpre(1:ndofn,n),n,bestn(n)
      END DO
    END IF
    CALL flushf(lures)

    !     ***  writes the final coordinates

    IF(iwrit == 1) THEN
      IF(ndime == 2) THEN
        IF(neulr == 0) WRITE(lures,"(/'  Coordinates After Generation',/,      &
                 &              ' Node     X',9X,'Y'/,(I5,2E15.5))",ERR=9999)  &
                 &             (label(n),coord(1:ndime,n),n=1,npoin)
        IF(neulr == 1) WRITE(lures,"(/' Coordinates After Generation',/,       &
                 &              ' Node     X',9X,'Y',8X,'Alfa'/,(I5,3E15.5))",ERR=9999)  &
                 &             (label(n),coord(1:ndime,n),euler(1,n),n=1,npoin)
      ELSE
        WRITE(lures,"(/'  Coordinates After Generation',/, &
                 &     '  Node     X',9X,'Y',9X,'Z')",ERR=9999)
        DO n=1,npoin
          WRITE(lures,"(I6,3E15.5)",ERR=9999)label(n),(coord(k,n),k=1,ndime)
        END DO
        IF (neulr == 9) THEN
          WRITE(lures,"(/'  Nodal Systems After Generation',/, &
             &  '   Node  X1     X2     X3',7x,'Y1     Y2     Y3',7x, &
             &  'Z1     Z2     Z3')",ERR=9999)
          DO n=1,npoin
            IF(.NOT.naeul(n) )CYCLE
            WRITE(lures,"(i7,3(2x,3f7.4))",ERR=9999) label(n),(euler(k,n),k=1,neulr)
          END DO
        END IF
      END IF
    END IF
    IF( itemp ) CALL fixtem(iwrit,npoin,iftmp,label)

    CALL flushf(lures)


    !***  Allocate space for arrays associated to NEQ

    ALLOCATE ( maxav(neq+1) )
    CALL ubicmx ( ndofn,melem,newnn,neq,maxa,ifpre,gnods,  &
                      maxav,npsdf,nesdf)
    DEALLOCATE ( gnods, bestn )

    ALLOCATE( displ(neq),ddisp(neq),gvect(neq) )
    displ = 0d0
    ddisp = 0d0
    resid = 0d0
    ALLOCATE(  disax(neq+1,7) )
    disax = 0d0              !initializes

    i = 1
    IF(nsymm > 0) i = i+1
    IF(nbuck > 0) i = i+1
    ALLOCATE( stiff(maxa,i) , resid(ndofn,npoin)  )
    resid = 0d0
    SELECT CASE (ndyna)
    CASE (0)
      ALLOCATE( mass(1), damp(1), coor1(ndime,npoin) )
      IF( neulr > 0 )THEN
        ALLOCATE( locs1(neulr,npoin),eule1(neulr,npoin) )
      ELSE
        ALLOCATE( locs1(1,1), eule1(1,1) )
      END IF
    CASE (1)
      ALLOCATE( mass(maxa), damp(maxa) )
    CASE (2)
      ALLOCATE( mass(neq), damp(maxa) )
    CASE (3)
      ALLOCATE( mass(maxa), damp(1) )
    CASE (4)
      ALLOCATE( mass(neq), damp(1) )
    END SELECT
    IF(ndyna == 0)THEN
      ALLOCATE( veloc(1), accel(1), velnp(1,1)  )
    ELSE
      ALLOCATE( veloc(neq), accel(neq)  )
      accel = 0d0
      veloc = 0d0
      !restore values for velocity
      DO i = 1,npoin
        DO j = 1,ndofn
          k = ifpre(j,i)
          IF (k > 0) veloc(k) = velnp(j,i)
        END DO
      END DO
    END IF
    CALL flushf(lures)

  END IF ! IF (actio /= 'NSTRA0')
  IF( nreql > 0 ) CALL cmp_slist( )

  RETURN
  9999 CALL runen2(' ')
END SUBROUTINE eqnums
   SUBROUTINE cmp_slist ( )
      !     generate list of dependencies to compute
      !     internal forces over master nodes

     USE ctrl_db, ONLY:  ndofn, npoin
     USE npo_db, ONLY: ifpre
     USE kinc_db, ONLY : nsdof,npsdf,nesdf,nn
     USE outp_db

     IMPLICIT NONE

     INTEGER (kind=4) :: i,j,k,l,m1,m2,n,ii,jj
     TYPE (slist), POINTER :: slave_d,slave_h,slave_t
     TYPE (slave_list), POINTER :: sl_d


     IF(ASSOCIATED(res))DEALLOCATE(res) !deallocate array if exist
     ALLOCATE(res(ndofn,nreql))         !get memory
     IF( nsdof == 0) RETURN             !If no slave DOF nothing to do
     !release previous memory if allocated
     DO
       IF( .NOT.ASSOCIATED(sl_head))EXIT
       IF( sl_head%nvalues > 0 )DEALLOCATE(sl_head%deps)
       sl_d => sl_head%next
       DEALLOCATE(sl_head)
       sl_head => sl_d
     END DO
     !
     DO i=1,nreql         !for each required node
       DO j=1,ndofn         !for each DOF of the node
         ALLOCATE(sl_d)       !a pointer is required
         sl_d%nvalues = 0     !initializes number of slave DOFs
         n = ifpre(j,nprql(i))       !master type of DOF
         IF( n > 0 .OR. n < -nn ) THEN     !if free or prescribed
           ii = 1                            !initializes counter
           DO k=1,nsdof                      !search in list of slave DOFs
             m1 = npsdf(k)                     !first position of this slave DOF
             m2 = npsdf(k+1)-1                 !last position of this slave DOF
             DO l=m1,m2                        !search  in list
               IF( nesdf(l) == n )THEN           !see if DOF is present
                 sl_d%nvalues = sl_d%nvalues + 1   !increase number of dependent DOFs
                 ALLOCATE(slave_d)                 !get a pointer
                 ext : DO               !loop over all the nodes
                   DO jj=1,ndofn                      !and all the DOFs
                     IF( ifpre(jj,ii) == -k )EXIT ext    !until found
                   END DO
                   ii= ii+1
                   IF( ii > npoin )ii=1
                 END DO ext
                 slave_d%dof  = jj                 !keep DOF
                 slave_d%node = ii                 !and NODE
                 slave_d%pos  = l                  !and POSition in NESDF
                 IF( sl_d%nvalues == 1 )THEN       !if the first time
                   slave_h => slave_d                !set as head
                 ELSE                              !already occured
                   slave_t%next => slave_d           !add at the end
                 END IF
                 slave_t => slave_d                !point to the last value
               END IF
             END DO
           END DO
         END IF
         IF( sl_d%nvalues > 0 )THEN                !if the DOF is master of any DOF
           ALLOCATE(sl_d%deps(3,sl_d%nvalues))        !get memory for dependencies
           slave_d => slave_h                         !point to first slave
           DO k=1,sl_d%nvalues                        !loop to generate list
             sl_d%deps(:,k) = (/ slave_d%dof,slave_d%node,slave_d%pos /)   !add values
             slave_t => slave_d                       !point to this value to delete
             slave_d => slave_t%next                  !point to next value before deleting
             DEALLOCATE(slave_t)                      !release memory
           END DO
         END IF
         IF( i == 1 .AND. j == 1)THEN              !for the first node and dof
           sl_head => sl_d                            !remember first position in the list
         ELSE
           sl_tail%next => sl_d                       !add to the list
         END IF
         sl_tail => sl_d                              !update list tail
       END DO
     END DO

     RETURN
   END SUBROUTINE cmp_slist
