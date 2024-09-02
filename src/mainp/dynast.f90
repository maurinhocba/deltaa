 SUBROUTINE dynast(istop,actio)
 !********************************************************************
 !
 ! *** iterative dynamic process routine
 !
 !********************************************************************
 USE lispa0, ONLY : lures
 USE c_input, ONLY : openfi
 USE ctrl_db, ONLY : ndime,neulr,ndofn,nrotd,npoin,neq,maxa,nload,istep,ndyna,dtime,ttime
 USE npo_db, ONLY  : coora,coorc,coori,euler,locsy,ifpre,psia,psic, &
                     force,loadv,loass,resid,veloc,velnp,accel,mass,damp
 USE solv_db
 USE dyna_db
 USE kinc_db, ONLY : ndepd,naris, maxav,nvelr,velor,nvfix,lcvel,nesdf,npsdf,ftsdf

 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN OUT) :: istop
 CHARACTER(len=*),INTENT(IN):: actio

 LOGICAL ::  convg,dampg
 INTEGER(kind=4) :: itera,i,first,every,ittol,niter,n1,nld,nsu,j,k,l
 LOGICAL :: kstif

 REAL (kind=8) :: initw,deltw,normf,normr,normd,nordi,             &
                  cm,cv,c1,c2,c3,c4,ca0,cb0,poter,ntvel,           &
                  facts(nload,2),factv(nvelr,2),ddt,ff(nload)
 REAL (kind=8),SAVE :: ntdis

 INTERFACE
   INCLUDE 'actual.h'
   INCLUDE 'colsol.h'
   INCLUDE 'coract.h'
   INCLUDE 'functs.h'
   INCLUDE 'mulmat.h'
   INCLUDE 'outbuc.h'
   INCLUDE 'resvpl.h'
   INCLUDE 'stiffm.h'
   INCLUDE 'timing.h'
   INCLUDE 'velnpo.h'
 END INTERFACE

 !        COMPUTE constants and factors

 first = MOD(nralg,10)          !first iteration to update stiffness matrix
 every = MOD(nralg,100)/10      !frequency of stiffness matrix evaluation
 ittol = MOD(nralg,1000)/100    !type of convergence check
 niter = nralg/10000000                                  !two digits
 IF(every == 0) every = niter + 1  !do not update
 ! compute factors used in computation
 cm = gamma/beta/dtime**2/alpha    !
 cv = gamma/beta/dtime
 dampg = (ndyna <= 2 ) .AND. (ccm > 0 .OR. cck > 0)  !damping present
 IF(dampg) THEN
   ca0 = cm +cv*ccm
   cb0 = 1d0+cv*cck
 ELSE
   ca0 = cm
   cb0 = 1d0
 END IF

 nsu  = 1 + ABS(nsymm)             !position of non-symmetric part
 n1 = nload+1                      !position of conformed force vector
 nld = n1
 IF(nvelr > 0) nld = nld+1         !position of prescrived velocities force vector

!  EXTRACTED FROM EXPLIT, BUT SOMETHING LIKE THIS MUST BE PRESENT
! ! update lumped mass matrix and force vector when dependant nodes exists
! IF (((ndepd+naris+nrfdp+nnsld ) > 0) .AND. (MOD(istep,20) == 0)) THEN !
!   IF (nload >0) force(1:neq,1:nload)=0d0
!   DO i=1,nload
!     CALL ensve1(ndofn*npoin,ifpre(1,1),loadv(1,1,i),force(1,i),npsdf(1),nesdf(1),ftsdf(1))
!   END DO
!   IF( lumped )THEN
!     ymass = 0d0
!     CALL ensmal(ndofn*npoin,neq,ifpre(1,1),emass(1,1),ymass(1),npsdf(1),nesdf(1),ftsdf(1))
!     IF( selective_mass_scaling )CALL fibre_mass(emass)
!     DO i=1,neq
!       ymass(i) = 1d0/ymass(i)
!     END DO
!   END IF
! END IF

 !        NATURAL MODES AND FREQUENCIES
 !        if required and only for first step
 !        this could be changed so that eigen-pairs are calculated
 !        for deformed configurations
 IF(eigen .AND. istep == 1) THEN
   !Compute stiffness matrix
   CALL resvpl (istop,ttime,resid,gvect)  ! ,.TRUE.)
   CALL stiffm (istop,ttime,.TRUE.,nld,stiff(:,1),force,nsymm,stiff(:,nsu))
   !number of eigenvectors to use in computation  => i
   i = MAX(2*neigen, neigen+8)
   i = MIN(i,neq)
   l =i*(i+1)/2       !size of square condensed matrix
   j = neq                        !size of lumped mass matrix
   IF(MOD(ndyna,2) /= 0)j = maxa  !size of consistent mass matrix
   ! reserve memory for auxiliar arryas
   ALLOCATE (r(0:neq,i), eigv(i), ar(l), br(l),vec(l,l), d(i), rtolv(i))
   ! initializes (perhaps unnecessary)
   r = 0d0 ; eigv = 0d0 ; ar = 0d0 ; br = 0d0 ; vec = 0d0 ; d= 0d0 ; rtolv = 0d0
   k = 1  !verifys sturm sequence
   ! compute eigen-pairs usin Subspace iteration algorithm
   CALL openfi(47)
   CALL  sspace (stiff(1,nsu),mass(1),maxav(1),r(0,1),eigv,disax(1,1),disax(1,2),      &
  &               ar(1),br(1),vec(1,1),d(1),rtolv(1),disax(1,3),disax(1,4),disax(1,5), &
  &               neq,neq+1,maxa,j,neigen,1e-8,i,l,32,k,1,47,55,0)
   CLOSE(47)
   ! output eigen-pairs
   DO k=1,neigen
     CALL outbuc(npoin,ndofn,k,ifpre,r(1:neq,k),0d0,eigv(k), 3)
   END DO
   !release auxiliar memory
   DEALLOCATE (r, eigv, ar, br,vec, d, rtolv)
 END IF

 !     ********************  INITIAL ACCELERATIONS  *********************
 IF(istep == 1) THEN    !for first step only
   kstif = .TRUE. !compute stiffness matrix the first time
   CALL resvpl (istop,ttime,resid,gvect)   !compute residual forces at t=0
   ! IF load present compute initial load vector
   IF(nload > 0)THEN
     DO i=1,nload
       facts(i,1:2) = functs(loass(i),ttime)*force(neq+1,i)
     END DO
     DO i=1,neq
       force(i,n1)= DOT_PRODUCT(force(i,1:nload),facts(1:nload,1))
       gvect(i) = force(i,n1) + gvect(i)  !residual forces
     END DO
   END IF
   IF(dampg) THEN  ! if damping compute C matrix
     ! stiffness matrix to compute contribution to C
     CALL stiffm (istop,ttime,.TRUE.,nld,stiff(:,1),force,nsymm,stiff(:,nsu))
     IF(ndyna == 1) THEN      !consistent mass matrix
       damp(1:maxa) = ccm*mass(1:maxa)+cck*stiff(1:maxa,1)
     ELSE  ! ndyna = 2       !lumped mass matrix
       DO i=1,maxa
         damp(i) = cck*stiff(i,1)
       END DO
       !damp(1:maxa) = cck*stiff(1:maxa,1)
       DO i=1,neq
         damp(maxav(i)) = damp(maxav(i))+ccm*mass(i)
       END DO
     END IF
     CALL mulmat(damp,maxav,ddisp,veloc,neq)   !contribution to equilibrium
     gvect = gvect - ddisp
   END IF
   CALL timing(5,.TRUE.)
   IF(MOD(ndyna,2) /= 0) THEN  !ndyna =1 or 3 consistent mass matrix
     accel = gvect                    !umbalanced forces
     DO i=1,maxa
        stiff(i,1) = mass(i)   !use stiff array as auxiliar
     END DO
     CALL colsol(stiff(:,1),maxav,neq,1,55,1,0)           !mass matrix must be symmetric
     CALL colsol(stiff(:,1),maxav,neq,2,55,1,0,v=accel)
   ELSE             !ndyna = 2 or 4 lumped mass matrix
     DO i = 1,neq
       accel(i) = gvect(i)/mass(i)
     END DO
   END IF
   CALL timing(5,.FALSE.)
   ntdis = 0d0
 END IF
 !     ****  END of initial accelerations  *****

 !     **************************     PREDICTION     ***********************

 c1 = alpha*dtime              !fraction time where equilibrium is checked
 c2 = c1*dtime*(0.5d0-beta)    !h^2(1/2 - beta) * alpha
 c4 = 1d0-gamma                !1-gamma
 c3 = c1*c4                    !alpha*h*(1-gamma)
 ! compute values at alpha*h
 DO i = 1,neq
   disax(i,1) =    veloc(i) + c3*accel(i)   !velocities         v_n + alpha*h*(1-gamma) a_n
   disax(i,2) = c4*accel(i)                 !mean acceleration  (1-gamma) a_n
   disax(i,3) = c1*veloc(i) + c2*accel(i)   !incremental displ. alpha*h v_n + h^2(1/2 - beta) * alpha * a_n
 END DO

 ! compute prescribed veloctities at equilibrium time
 ddt   = dtime*alpha           !initializes to  h*alpha
 IF(nvelr > 0)THEN
   DO i=1,nvelr
     factv(i,1:2) = functs(lcvel(i),ttime+ddt)*velor(nvfix+1 ,i)
   END DO
   velor(1:nvfix,nvelr+1) = MATMUL(velor(1:nvfix,1:nvelr),factv(1:nvelr,1))
 END IF
 ! update coordinates & nodal system from prediction
 DO i=1,npoin             !loop over each node
   coori(:,i) = coora(:,i)
 END DO
 CALL coract(ndime,npoin,neulr,ndofn,ddt,coora,disax(:,3),coora,locsy)

 !  UNBALANCED NODAL FORCES
 !  1- contribution from external forces
 IF(nload > 0)THEN
   IF((ndepd+naris) > 0 )THEN !regenerate force vector
     force = 0d0
     DO i=1,nload
       CALL ensvec (nrotd*npoin,ifpre,loadv(1,1,i),force(1,i))
     END DO
   END IF
   ff = 0
   DO i=1,nload
     facts(i,1:2) = functs(loass(i),ttime)*force(neq+1,i)
     ff(i) = facts(i,1)
     facts(i,1:2) = functs(loass(i),ttime+dtime)*force(neq+1,i)
     ff(i) = ff(i)*(1-alpha) + facts(i,1)*alpha
   END DO
   force(1:neq,n1)= MATMUL(force(1:neq,1:nload),ff )
 END IF
 !  2- contribution from internal forces
 CALL resvpl (istop,ttime,resid,gvect)
 IF(istop == 1)RETURN
 ttime = ttime + dtime       !update time (end of internal)
 !  3- inertia forces
 IF(MOD(ndyna,2) == 1) THEN  !consistent mass matrix
    CALL mulmat(mass,maxav,ddisp,disax(:,2),neq)
 ELSE                       !lumped mass matrix
   DO i=1,neq
     ddisp(i) = mass(i)*disax(i,2)
   END DO
 END IF
 !     internal + external + inertia forces
 gvect = gvect - ddisp
 IF(nload > 0) gvect = gvect + force(1:neq,n1)
 !  damping forces are added once damping matrix is computed
 !     *****  END prediction  ***
 !     set initial values
 convg = .FALSE.
 itera = 0
 !     *********************   iteration loop   ***********************
 DO
   itera = itera + 1     !update iteration counter
   !       ************   new stiffness matrix if necessary   ************
   IF(itera >= first .AND. MOD(itera-first,every) == 0)kstif = .TRUE.
   IF(kstif) THEN
     !   evaluates stiffness matrix and force vector for presc. displac.
     CALL stiffm (istop,ttime,.TRUE.,nld,stiff(:,1),force,nsymm,stiff(:,nsu))
     IF(dampg) THEN   !damping matrix evaluation
       IF(ndyna == 1) THEN   !consistent mass matrix
         damp(1:maxa) = ccm*mass(1:maxa)+cck*stiff(1:maxa,1)
       ELSE  ! ndyna = 2       !lumped mass matrix
         DO i=1,maxa
           damp(i) = cck*stiff(i,1)
         END DO
         !damp(1:maxa) = cck*stiff(1:maxa,1)
         DO i=1,neq
           damp(maxav(i))=damp(maxav(i))+ccm*mass(i)
         END DO
       END IF               !consistent or lumped
     END IF          !end of damping matrix evaluation
     IF(MOD(ndyna,2) == 0) THEN !lumped mass matrix
       stiff(1:maxa,1) = cb0*stiff(1:maxa,1)
       DO i=1,neq
         stiff(maxav(i),1) = stiff(maxav(i),1) + ca0*mass(i)
       END DO
     ELSE !consistent mass matrix
       DO i=1,maxa
         stiff(i,1) = ca0*mass(i)+cb0*stiff(i,1) !modified stiffness
       END DO
     END IF                    !consistent or lumped
     !         L D U factorization of effective stiffness matrix
     CALL timing(5,.TRUE.)
     CALL colsol(stiff(:,1),maxav,neq,1,55,1,nsymm,u=stiff(:,nsu))
     CALL timing(5,.FALSE.)
     kstif = .FALSE.
   END IF
   !       *** END of new stiffness matrix   ***
   !       ***** damping and prescribed displ. forces in first iteration   *****
   IF(itera == 1) THEN        !for first iteration
     IF(dampg) THEN               !sums damping forces
       CALL mulmat(damp,maxav,ddisp,disax(:,1),neq)
       gvect = gvect - ddisp
     END IF
     IF(nvelr > 0 ) THEN  !sums prescribed displ. forces
       gvect =  gvect + ddt*force(1:neq,nld)
     END IF
   END IF
   !       **** end of correction of residual forces in first iteration ***
   ddisp = gvect          !
   !       evaluates incremental displacements due to residual forces
   CALL timing(5,.TRUE.)
   !WRITE(58,"(20e14.4,/)")ddisp
   CALL colsol(stiff(:,1),maxav,neq,2,55,1,nsymm,v=ddisp,u=stiff(:,nsu))
   !WRITE(58,"(20e14.5,/)")ddisp
   CALL timing(5,.FALSE.)
   CALL screed(istep,itera,ttime,disax(ecdis,3))
   !       values to check convergence
   deltw = DOT_PRODUCT(ddisp,gvect)
   nordi = SQRT(DOT_PRODUCT(ddisp,ddisp))
   IF(itera == 1) THEN
     normf = MAX(SQRT(DOT_PRODUCT(gvect,gvect)),toler)
     initw = MAX(ABS(deltw),toler)
     normd = MAX(nordi,toler)
     normr = SQRT(DOT_PRODUCT(gvect,gvect))          !residual norm
   END IF
   IF(ittol <= 2.AND.ABS(deltw/initw) > 1./toler)initw=deltw*toler
   IF(ittol <= 2.AND.ABS(normr/normf) > 1./toler)normf=normr*toler
   IF(ittol == 3.AND.ABS(nordi/normd) > 1./toler)normd=nordi*toler
   !       updates coordinates and nodal local systems
   DO i=1,npoin             !loop over each node
     coori(:,i) = coora(:,i)
   END DO
   CALL coract(ndime,npoin,neulr,ndofn,0d0,coora,ddisp,coora,locsy)
   !       updates (n+alpha) velocities, mean acceler. and incr. displacements,==> disax
   DO i = 1,neq
     disax(i,1)=disax(i,1)+cv*ddisp(i)        !velocity at n+alpha
     disax(i,2)=disax(i,2)+cm*ddisp(i)        !mean acceleration
     disax(i,3)=disax(i,3)+   ddisp(i)        !incremental displacements
   END DO
   !Rearrange velocity vector (necessary for some element routines) at n+alpha
   CALL velnpo(nrotd, npoin, nvelr, ifpre, nesdf, npsdf, ftsdf, velor, disax(:,1), velnp)
   !       ******  evaluates residual of the motion equations  ****
   CALL resvpl(istop,ttime,resid,gvect)            !internal forces
   IF(istop == 1)EXIT
   IF(MOD(ndyna,2) /= 0) THEN   !consistent mass
     CALL mulmat(mass,maxav,ddisp,disax(:,2),neq)  !inertia forces
   ELSE                        !lumped mass
     DO i=1,neq
       ddisp(i) = mass(i)*disax(i,2)               !inertia forces
     END DO
   END IF                      !end of inertia forces
   gvect = gvect - ddisp
   IF( nload > 0)gvect = gvect + force(1:neq,n1)   !external forces
   IF(dampg) THEN               !damping forces
     CALL mulmat(damp,maxav,ddisp,disax(:,1),neq)
     gvect = gvect - ddisp
   END IF
   !       *** END of evaluation of residual forces
   !       convergence checks
   normr = SQRT(DOT_PRODUCT(gvect,gvect))          !residual norm
   SELECT CASE (ittol)
   CASE (1)        !check incremental work
     IF( ABS(deltw/initw) < toler) convg = .TRUE.
   CASE (2)        !check umbalanced forces
     IF( ABS(normr/normf) < toler) convg = .TRUE.
   CASE (3)        !check incremental displ.
     IF( ABS(nordi/normd) < toler) convg = .TRUE.
   END SELECT
   !       PRINT convergence checks
   IF(itera == 1) THEN
     WRITE(lures,100)initw,normf,normd
     WRITE(55,100,ERR=9999)initw,normf,normd
   ELSE
     WRITE(lures,102)deltw/initw,normr/normf,nordi/normd,ittol,convg
     WRITE(55,102)deltw/initw,normr/normf,nordi/normd,ittol,convg
   END IF
  IF(itera == niter .AND. (.NOT.convg)) THEN
     IF(ittol == 2 .AND. ABS(deltw/initw) < toler) THEN
       convg = .TRUE.
     ELSE               !stops when no convergence
       WRITE(lures,200)istep,itera,ntvel,ttime
       istop = 1
       EXIT
     END IF
  END IF
   IF( convg ) EXIT
 END DO   !end of iterative procedure
 IF(convg)THEN
 !       computes kinetic and potential energy
 !        IF(MOD(ndyna,2) == 1) THEN  !ndyna =1 or 3 consistent mass matrix
 !         CALL mulmat(mass,maxav,disax(1,6),disax(1,1),neq)
 !          kinet = DOT_PRODUCT(disax(1:neq,6),disax(1:neq,1))
 !        ELSE  ! ndyna = 2 or 4  lumped mass matrix
 !          kinet = 0d0
 !          DO i = 1,neq
 !            kinet = kinet + mass(i)*disax(i,1)**2
 !          END DO
 !        END IF
   disax(1:neq,4) = disax(1:neq,4) + disax(1:neq,3)  !total displac.
   ! poter = - lambd*DOT_PRODUCT(force(1:neq,1),disax(1:neq,4))
   ntdis = disax(ecdis,4)
   ntvel = disax(ecdis,1)
   ! WRITE(UNIT=4,FMT="(' K=',e13.5,' V=',e13.5,' D=',  &
   !      &        e13.5,' v=',e13.5)") kinet/2d0, poter, ntdis, ntvel
   !       actualizes converged values
   WRITE(lures,101)itera
   WRITE(55,101)itera
   c1 = (1d0-alpha)/alpha
   !  c4 = (1d0-gamma)/gamma
   DO i=1,neq
     ddisp(i) = disax(i,3) * c1
     veloc(i) = veloc(i) + dtime*disax(i,2)
     accel(i) = (disax(i,2)-c4*accel(i))/gamma
     disax(i,4) = disax(i,4) + ddisp(i)
   END DO
   !   actualizes last configuration used to evaluate internal forces
   coorc = coora
   euler = locsy
   IF( ASSOCIATED (psia) )  psic = psia  !psi functions
   ddt = ddt*c1
   CALL coract(ndime,npoin,neulr,ndofn,ddt,coora,ddisp,coora,locsy)
   !!Rearrange velocity vector (necessary for some element routines)
   !!CALL velnpo(nrotd, npoin, nvelr, ifpre, nesdf, npsdf, ftsdf, velor, veloc, velnp)
   CALL actual(disax(:,4),ttime,dtime)

 ELSE
   istep = istep-1
   ttime = ttime - dtime
   coora = coorc
   locsy = euler
 END IF

RETURN
 9999 CALL runen2(' ')

 100 FORMAT('  iw=',e9.3,' nf=',e9.3,' id=',e9.3)
 101 FORMAT(' convergence in ',i3,' iterations ')
 102 FORMAT(' dw/iw=',e9.3,' nr/nf=',e9.3,' dd/id=',e9.3,              &
    & ' ittol=',i2,' convg=',l3)
 200 FORMAT( //,' non-convergence in step =',i5,' in ',i3,' iterations'&
    &/,' Control Displacement ',e12.4,' at time',e12.4,/,              &
    &' PROGRAM  s t o p p e d ')

 END SUBROUTINE dynast
