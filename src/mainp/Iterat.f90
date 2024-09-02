 SUBROUTINE iterat(istop,actio)
 !********************************************************************
 !
 ! *** iterative static process routine
 !
 !********************************************************************
 USE lispa0, ONLY : lures
 USE ctrl_db, ONLY : ndime,neulr,ndofn,nrotd,npoin,neq,maxa,nload,istep,itemp,ndoft
 USE npo_db, ONLY  : coord,coora,coorc,coori,coor1,euler,locsy,locs1,ifpre, &
                     force,loadv,loass,resid,psia,psic,psii,tempe,dtemp,iftmp
 USE solv_db
 USE gvar_db, ONLY : ksnxt           !computes Stiffness matrix in next iteration
 USE kinc_db, ONLY : ndepd,naris, maxav,nvelr,velor,nvfix,lcvel,nn,nn1
 USE ift_db, ONLY : npret,nprev,lctmp,prtmp
 !USE loa_db, ONLY : loass
 USE c_input, ONLY : openfi
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN OUT) :: istop

 CHARACTER(len=*),INTENT(IN):: actio

 LOGICAL, SAVE :: recompute = .TRUE.
 LOGICAL :: convg,crita
 INTEGER (kind=4),SAVE :: ephi,newtv
 INTEGER (kind=4) :: i,ittol,itera,first,every,karcl,jj,kstif,&
                     lauto,nsb,nsu,nld,n1,niter,l,j,k,ii,lmi(1),lms(1)
 REAL (kind=8),SAVE :: piter = 4d0
 REAL (kind=8) :: initw,deltw,normf,normr,ddlmb,normd,nordi,       &
                  diter,oldlb,aux1,facts(nload,2),factv(nvelr,2),factt(npret,2)
 REAL (kind=8), ALLOCATABLE :: aux(:)    ,sl(:),su(:)  ! for COLSOL verification
 !REAL (kind=8),SAVE :: poter = 0
 REAL (kind=8) :: tdisp
 REAL (kind=8) :: r0,r1,h,a,h1,dl

 INTERFACE
   INCLUDE 'actual.h'
   INCLUDE 'colsol.h'
   INCLUDE 'coract.h'
   INCLUDE 'functs.h'
   INCLUDE 'iteinv.h'
   INCLUDE 'mulmat.h'
   INCLUDE 'mulmatu.h'
   INCLUDE 'outbuc.h'
   INCLUDE 'resvpl.h'
   INCLUDE 'stiffm.h'
   INCLUDE 'timing.h'
   INCLUDE 'update.h'
 END INTERFACE

 IF((ndepd+naris) > 0 )THEN !regenerate force vector
   IF(nload > 0) force(1:neq,1:nload) = 0d0
   DO i=1,nload
     CALL ensvec (ndofn*npoin,ifpre(1,1),loadv(1,1,i),force(1,i))
   END DO
 END IF

 n1  = nload+1
 nld = n1
 IF( nvelr > 0 .OR. npret > 0 ) nld = nld + 1  !position in array
 ! compute prescribed velocities
 IF(nvelr > 0)THEN
   DO i=1,nvelr
     factv(i,1:2) = functs(lcvel(i),lambd)*velor(nvfix+1,i)
   END DO
   velor(1:nvfix,nvelr+1) = MATMUL(velor(1:nvfix,1:nvelr),factv(1:nvelr,1))
 END IF
 !   compute prescribed temperature velocities
 IF(npret > 0)THEN
   DO i=1,npret
     factt(i,1:2) = functs (lctmp(i),lambd)*prtmp(nprev+1,i)
   END DO

   !  update prescribed temperatures

   DO i=1,nprev
     prtmp(i,npret+1)= DOT_PRODUCT(prtmp(i,1:npret),factt(1:npret,1))
   END DO

   DO i=1,npoin             !loop over each node
     DO j=1,ndoft
       l = iftmp(j,i)       !DOF
       IF( l < -nn) tempe(j,i) = prtmp(-l-nn,npret+1)
     END DO
   END DO
   dtemp = tempe                      !keep original values

   !  update temperature velocity
   DO i=1,nprev
     prtmp(i,npret+1)= DOT_PRODUCT(prtmp(i,1:npret),factt(1:npret,2))
   END DO

 END IF

 oldlb = lambd

 first = MOD(nralg,10)
 every = MOD(nralg,100)/10
 ittol = MOD(nralg,1000)/100
 lauto = MOD(nralg,10000)/1000
 karcl = MOD(nralg,100000)/10000
 diter = REAL(INT(MOD(nralg,10000000)/100000),kind=8)    !two digits
 niter = nralg/10000000                                  !two digits
 IF( every == 0 .AND.first == 0 .AND.niter == 98 ) niter =1200
 IF(every == 0)every = niter + 1
 !
 IF( ecdis > 0 )  tdisp = coora(MOD(ncdis,10),ncdis/10)-coord(MOD(ncdis,10),ncdis/10)

 nsu  = 1 + ABS(nsymm)
 jj   = 1 + nsu                   !position of DSTIF
 nsb  = 0                         !stiffness symmetry for critical search
 IF(TRIM(actio) /= 'STEP' ) THEN
   piter = diter           !initializes
   ephi  = 0               !initializes
   kstif = 1               !compute stiffness matrix
   ksnxt = .TRUE.          !compute stiffness matrix in next iteration
   DO i=1,npoin             !loop over each node
     coori(:,i) = coora(:,i)
   END DO
   CALL resvpl (istop,lambd,resid,gvect,.TRUE.)
   normr = MAX(SQRT(DOT_PRODUCT(gvect,gvect)),toler)
   IF(nload > 0)THEN
     DO i=1,nload
       facts(i,1:2) = functs(loass(i),lambd)*force(neq+1,i)
     END DO
     DO i=1,neq
       force(i,n1)= DOT_PRODUCT(force(i,1:nload),facts(1:nload,1))
       gvect(i) = force(i,n1) + gvect(i)  !residual forces
     END DO
   END IF
 ELSE
   kstif = 0
 END IF
 !                      Initializes values for iteration
 convg = .FALSE.
 itera = 0
 !     **************** BEGIN ITERATION LOOP *******************
 DO
   itera = itera + 1
   IF(itera >= first .AND. MOD(itera-first,every) == 0 .OR. recompute)kstif = 1
   IF(kstif == 1) THEN
     kstif = 0 ; recompute = .FALSE.
     crita = (itera == first .AND. nbuck > 0 )
     IF(crita) crita = (MOD(istep,MOD(nbuck,1000)) == 0)
     ! Evaluates stiffness matrix and force vector for presc. displac
     IF( nvelr > 0 .OR. itemp) force(:,nld) = 0d0
     CALL stiffm (istop,lambd,.TRUE.,nld,stiff(:,1),force,nsymm,    &
                  stiff(:,nsu),.TRUE.)
!    prints STIFFNESS MATRIX
!    WRITE(58,"(2x,42(6x,i2,13x))",ERR=9999)(l,l=1,neq)
!    ALLOCATE(aux(neq))
!    DO i=1,neq
!      ! diagonal and below
!      j = maxav(i)
!      k = maxav(i+1)
!      DO l=i,1,-1
!        IF( j < k ) THEN
!          aux(l) = stiff(j,1)
!          j = j + 1
!        ELSE
!          aux(l) = 0d0
!        END IF
!      END DO
!      ! above diagonal
!      DO l=i+1,neq
!        j = maxav(l)            !diagonal position
!        k = maxav(l+1)          !next diagonal position
!        IF( k - j >  l-i) THEN
!          aux(l) = stiff(j+l-i,2)
!        ELSE
!          aux(l) = 0d0
!        END IF
!      END DO
!      WRITE(58,"(i2,12e20.8)",ERR=9999)i,(aux(l),l=1,neq)
!    END DO
!    DEALLOCATE(aux)

 !    IF(crita) stiff(:,jj) = - stiff(:,1) !stores actual Kt
     IF(crita)THEN !eigenvalue analysis to be performed
       DO i=1,maxa    !present stiffness matrix (sign reversed)
         stiff(i,jj) = -stiff(i,1)
       END DO
       WRITE(58,"(10E20.6))",ERR=9999)stiff(:,jj)
     END IF
     ! L D U factorization of stiffness matrix
     CALL timing(8,.TRUE.)
!     WRITE(58,"(12e15.5)",ERR=9999)(stiff(i,1:2),i=1,144)
     CALL colsol(stiff(:,1),maxav,neq,1,55,1,nsymm,u=stiff(:,nsu))
!!    prints STIFFNESS MATRIX
!     DO i=1,maxav(neq+1)
!       WRITE(58,"(i2,2e15.5)",ERR=9999)i,stiff(i,1:2)
!     END DO
     CALL timing(8,.FALSE.)
     ! evaluates tangent solution vector due to external forces and b.c.
     !external forces
     IF(nload > 0)THEN
       DO i=1,nload
         facts(i,1:2) = functs(loass(i),lambd)*force(neq+1,i)
       END DO
       !force(1:neq,n1)= MATMUL(force(1:neq,1:nload),facts(1:nload,2))
       !disax(1:neq,1) = force(1:neq,n1)
       DO i=1,neq
         force(i,n1)= DOT_PRODUCT(force(i,1:nload),facts(1:nload,2))
         disax(i,1) = force(i,n1)
       END DO
       !WRITE(55,"(f12.0)",ERR=9999)(disax(l,1),l=1,neq)
     ELSE
       disax(1:neq,1) = 0d0
     END IF
     !force vector due to prescribed displacements
     IF(nvelr > 0 .OR. itemp) THEN
       DO i=1,neq
         disax(i,1) = disax(i,1) - force(i,nld)
       END DO
     END IF
     CALL timing(8,.TRUE.)
     IF( crita .AND. neigen > 0 )THEN !buckling analysis will be performed
       CALL openfi(47)  ! only for neigen > in fact ==> SSPACE
       REWIND 47
       WRITE(47) stiff(:,1)           !keeps tangent matrix
     END IF
     !WRITE(58,"('f-vector')",ERR=9999)
     !WRITE(58,"(24e18.8)",ERR=9999)(disax(l,1),l=1,neq)
     CALL colsol(stiff(:,1),maxav,neq,2,55,1,nsymm,v=disax(1:neq,1), &  !computes Ut
                 u=stiff(:,nsu))
     newtv = 1
     CALL timing(8,.FALSE.)
     ! critical point search if desired
     IF(crita) THEN
       ! generates spatial configuration at (v + Ut*dl/10) ==> coor1,locs1
       aux1  = dlamb /1000d0          ! auxiliar load increment
       DO i=1,neq
         ddisp(i) = disax(i,1)*aux1
       END DO
       ! computes stresses and plastic consistency parameters and them
       ! stiffness matrix for a point along the tangent path
       ! generates spatial configuration with incremental displac. ddisp
       IF( neulr > 0 )THEN
         DO i=1,npoin
           locs1(:,i) = locsy(:,i)
         END DO
       END IF
       IF( ASSOCIATED(psia) )psii=psia !does it make sense? no load-geometry contributions
       CALL coract(ndime,npoin,neulr,ndofn,0d0,coora,ddisp,coor1,locs1)
       IF(npret > 0)THEN
         DO i=1,npret
           factt(i,1:2) = functs (lctmp(i),lambd+aux1)*prtmp(nprev+1,i)
         END DO
         !  update prescribed temperatures
         DO i=1,nprev
           prtmp(i,npret+1)= DOT_PRODUCT(prtmp(i,1:npret),factt(1:npret,1))
         END DO
         DO i=1,npoin             !loop over each node
           DO j=1,ndoft
             l = iftmp(j,i)       !DOF
             IF( l < -nn) tempe(j,i) = prtmp(-l-nn,npret+1)
           END DO
         END DO
         !  back to previous temperature velocities
         DO i=1,nprev
           prtmp(i,npret+1)= DOT_PRODUCT(prtmp(i,1:npret),factt(1:npret,2))
         END DO
       END IF
       ! computes stresses, plastic consistency pars. and residual forces
       ksnxt = .TRUE.
       CALL resvpl(istop,lambd,resid,ddisp) ! ==> ddisp
       CALL stiffm(istop,lambd,.FALSE.,4,stiff(:,jj),disax,nsymm, stiff(:,nsu))   !forces due to displacements stored in disax4
       ! IF( ASSOCIATED(psia) )psia=psii ! may be here it should restore psia
       IF(npret > 0) tempe = dtemp !restore temperatures
       ! Numerical derivate of stiffness matrix, symmetric part only
       aux1  = -1d0/aux1              ! inverse of auxiliar load increment
       DO i=1,maxa
         stiff(i,jj) = stiff(i,jj)*aux1 ! derivative of stiffness matrix
       END DO
       buckl = lambd                  ! base load factor
       IF( neigen == 1 )THEN
         CALL iteinv(stiff(:,1),stiff(:,jj),maxav(:),disax(:,4),     &
                     buckl,1d-9,neq,mbite,3,ephi)
         ephi = 1
         CALL outbuc(npoin,nrotd,istep,ifpre,disax(:,4),lambd,buckl, -1)
       ELSE
         !number of eigenvectors to use in computation  => i
         i = MAX(2*neigen, neigen+8)
         i = MIN(i,neq)
         l =i*(i+1)/2       !size of square condensed matrix
         j = maxa           !size derivative matrix
         ! reserve memory for auxiliar arryas
         IF( .NOT.ALLOCATED(r) )THEN
            ALLOCATE (r(0:neq,i), eigv(i), ar(l), br(l),vec(l,l), d(i), rtolv(i))
           ! initializes (perhaps unnecessary)
            r = 0d0 ; eigv = 0d0 ; ar = 0d0 ; br = 0d0 ; vec = 0d0 ; d= 0d0 ; rtolv = 0d0
         END IF
         k = 1  !verifys sturm sequence
         ! compute eigen-pairs usin Subspace iteration algorithm
         !WRITE(47) stiff(:,1)
         WRITE(47) disax
         CALL  sspace (stiff(1,1),stiff(1,jj),maxav(1),r(0,1),eigv,disax(1,1),disax(1,2),    &
        &               ar(1),br(1),vec(1,1),d(1),rtolv(1),disax(1,3),disax(1,4),disax(1,5), &
        &               neq,neq+1,maxa,maxa,neigen,1e-8,i,l,32,k,1,47,58,1)
         !READ(47) stiff(:,1)
         READ(47) disax
         CLOSE(47)
         ! output eigen-pairs
         DO k=1,neigen
           CALL outbuc(npoin,nrotd,istep,ifpre,r(1:neq,k),lambd,eigv(k),-k)
         END DO
       END IF
     END IF    !critical point search
   END IF    !kstif == 1
   ! evaluates incremental displacements due to residual forces
   DO i=1,neq
     ddisp(i) = gvect(i)
   END DO
   CALL timing(8,.TRUE.)
   lms = MAXLOC(ddisp)
   !WRITE(55,"(i5,e15.5)",ERR=9999)lms,ddisp(lms(1))
   !WRITE(55,"(10e15.5)",ERR=9999)(ddisp(l),l=1,neq,4)
   CALL colsol(stiff(:,1),maxav,neq,2,55,1,nsymm,v=ddisp,u=stiff(:,nsu))
   CALL timing(8,.FALSE.)
   ! updates load factor and incre. displ. considering selected algorithm
   CALL update(istep,itera,lauto,ecdis,neq,arcln,lambd,dlamb,      &
               disax,displ,ddisp,ddlmb,karcl,piter,diter,newtv,ncdis)
   IF( npret > 0 )THEN !  update temperatures
     DO i=1,npoin             !loop over each node
       DO j=1,ndoft             !loop over each node
         l = iftmp(j,i)       !DOF
         IF( l <= -nn1) tempe(j,i) = dtemp(j,i) + prtmp(-l-nn,npret+1)*dlamb
       END DO
     END DO
   END IF

   IF( ecdis > 0  )THEN                            !control displacement present
     tdisp = tdisp + ddisp(ecdis)                  !total displacement
     CALL scree2(istep,itera,ddlmb,lambd,tdisp,ddisp(ecdis))
   ELSE
     CALL screen(istep,itera,ddlmb,lambd)
   END IF
   ! evaluates total residual forces for this iteration
   IF(nload > 0 )THEN
     DO i=1,neq
       gvect(i) = gvect(i) + force(i,n1)*ddlmb
     END DO
   END IF
   deltw = DOT_PRODUCT(ddisp,gvect)
   nordi = SQRT(DOT_PRODUCT(ddisp,ddisp))
   IF(itera == 1) THEN
     normf = MAX(SQRT(DOT_PRODUCT(gvect,gvect)),toler)
     initw = MAX(ABS(deltw),toler)
     normd = MAX(nordi,toler)
   END IF
   IF(ittol <= 2.AND.ABS(deltw/initw) > 1./toler)initw=deltw*toler
   IF(ittol <= 2.AND.ABS(normr/normf) > 1./SQRT(toler))normf=normr*SQRT(toler)
   IF(ittol == 3.AND.ABS(nordi/normd) > 1./toler)normd=nordi*toler
   ! updates coordinates and nodal local systems
   DO i=1,npoin             !loop over each node
     coori(:,i) = coora(:,i)
   END DO
   IF( ASSOCIATED(psia) )psii=psia
   CALL coract(ndime,npoin,neulr,ndofn,ddlmb,coora,ddisp,coora,locsy)
   ! updates step displacements

   DO i=1,neq
     displ(i) = displ(i) + ddisp(i)
   END DO
   !WRITE(55,"(6e15.5)",ERR=9999)(ddisp(l),l=1,neq)

   ksnxt = (itera+1 >= first .AND. MOD(itera+1-first,every) == 0)
   CALL resvpl (istop,lambd,resid,gvect,.TRUE.)

   IF( istop == 1 )EXIT
   IF(nload > 0)THEN
     DO i=1,nload
       facts(i,1:2) = functs(loass(i),lambd)*force(neq+1,i)
     END DO
     DO i=1,neq
       gvect(i) = DOT_PRODUCT(force(i,1:nload),facts(1:nload,1))+gvect(i)
     END DO
     !WRITE(58,"(8e13.5)")-gvect
   END IF
   !r0 = MINVAL(gvect)*2d0/3d0
   !r1 = MAXVAL(gvect)*2d0/3d0
   !k = 1
   !DO i=1,neq
   !  IF( gvect(i) > r1 .OR. gvect(i) < r0)THEN
   !    loop : DO
   !      DO j=1,ndofn
   !        IF( ifpre(j,k) == i )THEN
   !          WRITE(58,"(i7,i3,e15.5)")k,j,gvect(i)
   !          EXIT loop
   !        END IF
   !      END DO
   !      k = k + 1
   !    END DO loop
   !  END IF
   !END DO
   ! convergence checks
   normr = SQRT(DOT_PRODUCT(gvect,gvect))
   SELECT CASE (ittol)
   CASE (1)
     IF(ABS(deltw/initw) < toler .AND. &
     (ABS(normr/normf) < SQRT(toler) .OR. ABS(nordi/normd) < toler )) &
        convg = .TRUE.
   CASE (2)
     IF(ABS(normr/normf) < toler) convg = .TRUE.
   CASE (3)
     IF(ABS(nordi/normd) < toler) convg = .TRUE.
   END SELECT
   IF(.NOT.convg) newtv = 0

   IF(itera == 1) WRITE(lures,100,ERR=9999)deltw,initw,normr,normf,nordi,normd
   WRITE(55,100,ERR=9999)deltw,initw,normr,normf,nordi,normd
   IF(itera == niter .AND. (.NOT.convg)) THEN
     IF(ittol == 2 .AND. ABS(deltw/initw) < toler) THEN
       convg = .TRUE.
     ELSE
       IF(niter >= diter) THEN
         WRITE(lures,200,ERR=9999)istep,itera,lambd,ddlmb
         EXIT
       ELSE                  !go to next step but recompute matrix
         convg = .TRUE.
         recompute = .TRUE.
       END IF
     END IF
   END IF
   IF( convg ) EXIT
   ! line search ***********************
   IF( itera == 1 )CYCLE
   IF(.NOT.lsearch) CYCLE  !line search disabled
   !IF( karcl /= 0 )CYCLE
   r0 = deltw                          !R(0)
   l = 1
   DO
     r1 = DOT_PRODUCT(gvect,ddisp)       !R(1)
     a = r0/r1                           !R(0)/R(1)
     IF( ABS(a) > 2.0 )EXIT
     IF( a < 0d0 )THEN        !a root exists
       h = a/2d0 + SQRT(a**2/4d0 - a )     !root approximation
     ELSE                                  !no root
       h = a/2d0                           !null tangent point
     END IF
     WRITE(55,"('L',i2,' R(0):',e12.4,' R(1):',e12.4,' h:',f6.2)")l,r0,r1,h
     WRITE(*,"('L',i2,' R(0):',e12.4,' R(1):',e12.4,' h:',f6.2)")l,r0,r1,h
     IF( h <= 0.5d0 ) recompute = .TRUE.     !improve matrix
     h1 = h-1d0
     disax(1:neq,7) = h1*ddisp               !change in displacements
     ddisp = h*ddisp                         !modify incremental displacement
     dl = h1*ddlmb                           !change in lambda
     ddlmb = h*ddlmb                         !modify incremental lambda
     displ = displ + disax(1:neq,7)       !modify step displacements
     lambd = lambd + dl                   !modify load factor
     DO i=1,npoin             !loop over each node
       coori(:,i) = coora(:,i)
     END DO
     CALL coract(ndime,npoin,neulr,ndofn,dl,coora,disax(1:neq,7),coora,locsy) !update coordinates
     CALL resvpl (istop,lambd,resid,gvect,.TRUE.) !recompute internal forces
     IF( istop == 1 )EXIT
     IF(nload > 0)THEN
       DO i=1,nload
         facts(i,1:2) = functs(loass(i),lambd)*force(neq+1,i)
       END DO
       gvect = MATMUL(force(1:neq,1:nload),facts(1:nload,1))+gvect
     END IF
     l = l+1
     IF( l > 3 )EXIT
   END DO
   ! end line search ******************
 END DO
 !     actualizes converged values
 IF ( convg )THEN
   WRITE(lures,"(' convergence in ',i3,' iterations ')",ERR=9999)itera
   WRITE(lures,100,ERR=9999)deltw,initw,normr,normf,nordi,normd
   !WRITE(55,100,ERR=9999)deltw,initw,normr,normf,nordi,normd
   WRITE(55,"(' convergence for step ',i4,' in ',i3,' iterations ')",ERR=9999)istep,itera
   piter = REAL(itera,kind=8)
   DO i=1,npoin           !updates converged configuration
     coorc(:,i) = coora(:,i)
   END DO
   IF( neulr > 0 )THEN
     DO i=1,npoin           !updates converged configuration
       euler(:,i) = locsy(:,i)
     END DO
   END IF
   IF( ASSOCIATED (psia) )  psic = psia  !psi functions
   dlamb = lambd - oldlb
   !CALL actual(displ,lambd,dlamb)
   !   WRITE(55,"(' convergence, itera',i5,' l =',e12.4)",ERR=9999)itera,lambd
   !   poter = poter - lambd*DOT_PRODUCT(force(1:neq,1),displ)
   !   tdisp = tdisp + displ(ecdis)
   !   WRITE(UNIT=55,ADVANCE='no',FMT="(' D=',e13.5,' l =',e12.4,' V=', &
   !     &      e13.5)",ERR=9999)tdisp,lambd,poter
   IF( npret > 0 )THEN !  update temperatures
     DO i=1,npoin             !loop over each node
       DO j=1,ndoft             !loop over each node
         l = iftmp(j,i)       !DOF
         IF( l <= -nn1) tempe(j,i) = dtemp(j,i) + prtmp(-l-nn,npret+1)*dlamb
       END DO
     END DO
   END IF
 ELSE
   istop = 1
   WRITE(lures,200,ERR=9999)istep,itera,lambd,ddlmb
   istep = istep - 1
   lambd = oldlb
   coora = coorc
   locsy = euler
 END IF
 RETURN
 9999 CALL runen2(' ')
 100 FORMAT(' dw=',e9.3,' iw=',e9.3,' nr=',e9.3,' nf=',e9.3,           &
    &       ' dd=',e9.3,' id=',e9.3)
 200 FORMAT( //,' non-convergence in step =',i5,' in ',i3,' iterations'&
    &/,' actual load PARAMETER ',e12.4,' last increment in load',      &
    &e12.4,/,' PROGRAM  s t o p p e d ')

 END SUBROUTINE iterat
