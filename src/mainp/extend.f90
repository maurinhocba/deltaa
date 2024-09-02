 SUBROUTINE extend(istop,actio)
 !********************************************************************
 !
 ! *** Extended system iterative process routine
 !
 !********************************************************************
 USE lispa0
 USE curv_db
 USE ctrl_db, ONLY : ndime,neulr,ndofn,npoin,nload,istep
 USE npo_db, ONLY : coord,coora,coorc,coor1,euler,locsy,locs1,ifpre, &
                    force,loadv,loass,resid
 USE solv_db
 USE kinc_db, ONLY : ndepd,naris,          &
                     neq,maxa,maxav,       &
                     nvelr,velor,nvfix,lcvel
 USE ift_db, ONLY : npret,nprev,lctmp,prtmp
 !USE loa_db, ONLY :
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN OUT) :: istop
 CHARACTER(len=*),INTENT(IN):: actio

 ! local variables
 LOGICAL :: convg

 INTEGER (kind=4) :: itera,i,j,jj,ecphi(1),diter,ittol,nsb,        &
                     nsu,nld,n1,niter
 REAL (kind=8) :: initw,deltw,normf,normr,dl,normd,g1,g2,dmu,      &
                  nordi,oldlb,aux1,a11,a12,a21,a22,mu,             &
                  factv(nvelr,2),facts(nload,2)

 INTERFACE
   INCLUDE 'actual.h'
   INCLUDE 'colsol.h'
   INCLUDE 'coract.h'
   INCLUDE 'functs.h'
   INCLUDE 'mulmat.h'
   INCLUDE 'outbuc.h'
   INCLUDE 'resvpl.h'
   INCLUDE 'stiffm.h'
 END INTERFACE

 jj    = 2+ABS(nsymm)              !first position in DSTIF
 oldlb = lambd                     !keeps original load factor
 ittol = MOD(nralg,1000)/100       !convergence type
 diter = REAL(INT(MOD(nralg,10000000)/100000),kind=8)    !two digits
 niter = nralg/10000000                                  !two digits
 nsb   = 0
 nsu   = 1 + ABS(nsymm)

 IF((ndepd+naris) > 0 )THEN !regenerate force vector
   IF(nload > 0) force(1:neq,1:nload) = 0d0
   DO j=1,nload
     CALL ensvec (ndofn*npoin,ifpre(1,1),loadv(1,1,j),force(1,j))
   END DO
 END IF

 n1  = nload+1
 nld = n1
 IF( nvelr > 0 .OR. npret > 0 ) nld = nld + 1  !position in array
 ! compute prescribed velocities
 IF(nvelr > 0)THEN
   DO j=1,nvelr
     factv(j,1:2) = functs(lcvel(j),lambd)*velor(nvfix+1,j)
   END DO
   velor(1:nvfix,nvelr+1) = MATMUL(velor(1:nvfix,1:nvelr),factv(1:nvelr,1))
 END IF

 !     most relevant DOF in eigenvector  ==> i
 ecphi = MAXLOC(disax(1:neq,4))    !maximum value = +1
 i = ecphi(1)
 WRITE(lures,"(' Spring on DOF=',i5,' Ut(i)=',e12.4,' Phi(i)=',    &
              &  e12.4)") i, disax(i,1), disax(i,4)
 !*******************
 !     prediction of critical state
 dl    = buckl - lambd             !lambda increment
 displ = dl * disax(1:neq,1)       !displacement increment
 lambd = buckl                     !prediction for critical load
 mu    = displ(i)                  !increment in most relevant DOF
 !     updates coordinates, nodal local systems and computes RESID
 CALL coract(ndime,npoin,neulr,ndofn,dl,coora,displ,coora,locsy) !=> coora & locsy
 CALL resvpl(istop,lambd,resid,gvect,flag=.TRUE.)                !=> resid & gvect
 IF(nload > 0)THEN
   DO j=1,nload        !present load factor and derivatives
     facts(j,1:2) = functs(loass(j),lambd)*force(neq+1,j)
   END DO
   force(1:neq,n1)=MATMUL(force(1:neq,1:nload),facts(1:nload,2)) !incremental (derivative) force vector
   gvect = MATMUL(force(1:neq,1:nload),facts(1:nload,1))+gvect   !umbalanced equivalent nodal forces
 END IF
 ! gvect(i) = gvect(i) + gamm1*(mu-displ(i))     !useless
 !     computes comparison values to check convergence
 initw = MAX(ABS(DOT_PRODUCT(displ,gvect)),toler)  !energy norm
 normf = MAX(SQRT(DOT_PRODUCT(gvect,gvect)),toler) !residual norm
 normd = MAX(SQRT(DOT_PRODUCT(displ,displ)),toler) !displacement norm
 !     Initializes values for iteration
 convg = linear
 itera = 0
 CALL screen(istep,itera,dl,lambd)
 !     begin iteration loop with Full Newton-Raphson
 DO
   IF( convg ) EXIT
   itera = itera + 1
   WRITE(lures,"(/'iter',i2,' L',e15.5,' dl',e15.5)")itera,lambd,dl
   WRITE(58,"(/,'iter',i2,' L',e15.5,' dl',e15.5)")itera,lambd,dl
   !    evaluates stiffness matrix and force vector for presc. displac
   CALL stiffm(istop,lambd,.TRUE.,nld,stiff(:,1),force,nsymm,stiff(:,nsu)) !=> Stiff(1) & force(nld)
   !    generates spatial configuration at (v + alpha * phi) ==> coor1, locs1
   ddisp =alph1*disax(1:neq,4)       !ALPH1 in SOLV_DB  read in CONTOL
   !    computes stresses and plastic consistency parameters and them
   !    stiffness matrix for a point along the bifurcation path
   locs1 = locsy
   CALL coract(ndime,npoin,neulr,ndofn,0d0,coora,ddisp,coor1,locs1)  ! => coor1 & locs1 (default arrays for RESVPL)
   ! computes stresses, plastic consistency pars. and residual forces
   CALL resvpl(istop,lambd,resid,ddisp)  ! => resid & ddisp
   ! stiffness matrix at altered configuration
   CALL stiffm(istop,lambd,.TRUE.,1,stiff(:,jj),disax,nsymm,stiff(:,nsu)) ! => stiff(jj) & disax(1)
   !    numerical derivative of (symmetric) stiffness along bifurcation path
   DO j = 1,maxa
     stiff(j,jj) = (stiff(j,1)-stiff(j,jj))/alph1  ! - Du(Kt.Phi)
   END DO
   stiff(maxav(i),1) = stiff(maxav(i),1) + gamm1 ! modify stiffness matrix (diagonal value associated with largest dof in eigenvalue)
   !    L D L factorization of stiffness matrix
   CALL colsol(stiff(:,1),maxav,neq,1,58,1,0,u=stiff(:,nsu))
   WRITE(58,"('DOF',i6,'D(i)',e12.4)")i,stiff(maxav(i),1)
   !    evaluates tangent solution vector due to external forces and b.c
   !external forces
   IF(nload > 0)THEN
     DO j=1,nload      !present load factor and derivatives
       facts(j,1:2) = functs(loass(j),lambd)*force(neq+1,j)
     END DO
     force(1:neq,n1)=MATMUL(force(1:neq,1:nload),facts(1:nload,2))  !incremental (derivative) force vector
     disax(1:neq,1) = force(1:neq,n1) + disax(1:neq,1)              !equivalent nodal forces due to loads and prescribed displacements (derivative)
   END IF
   !incremental displacements due to increments in loads and prescribed displacements
   CALL colsol(stiff(:,1),maxav,neq,2,58,1,0,v=disax(:,1),u=stiff(:,nsu)) !=> disax(1)
   !    evaluates incremental displacements due to residual forces
   disax(1:neq,2) = gvect
   CALL colsol(stiff(:,1),maxav,neq,2,58,1,0,v=disax(:,2),u=stiff(:,nsu)) !=> disax(2)
   !    evaluates incremental displacements due to vector gamma*e
   disax(1:neq,3) = 0d0
   disax(i,3) = gamm1
   CALL colsol(stiff(:,1),maxav,neq,2,58,1,0,v=disax(:,3),u=stiff(:,nsu))
   WRITE(58,"('disax',3e12.4)")disax(i,1:3)
   !    computes vectors h
   CALL mulmat(stiff(:,jj),maxav,disax(:,5),disax(:,1),neq)
   CALL mulmat(stiff(:,jj),maxav,disax(:,6),disax(:,2),neq)
   CALL mulmat(stiff(:,jj),maxav,disax(:,7),disax(:,3),neq)
   !    computes vectors p
   CALL colsol(stiff(:,1),maxav,neq,2,58,1,0,v=disax(:,5),u=stiff(:,nsu))
   CALL colsol(stiff(:,1),maxav,neq,2,58,1,0,v=disax(:,6),u=stiff(:,nsu))
   CALL colsol(stiff(:,1),maxav,neq,2,58,1,0,v=disax(:,7),u=stiff(:,nsu))
   !    solves for dLambda y dMu
   a11 = disax(i,5)
   a12 = disax(i,7)
   a21 = disax(i,1)
   a22 = disax(i,3) - 1d0
   g1 = 1d0 - disax(i,3)  - disax(i,6)
   g2 = mu  - displ(i)    - disax(i,2)
   aux1 = a11*a22 - a12*a21
   IF(ABS(aux1) > 1d-6 .AND. ABS(a22) > 1e-8) THEN
     dl  = ( a22*g1 - a12*g2)/aux1  ! incremental load factor
     dmu = (-a21*g1 + a11*g2)/aux1  ! incremental most relevant DOF
   ELSE
     ! dl  = (1d0 - disax(i,6))/disax(i,5)
     ! dmu = - gamm1*disax(i,3)/disax(i,7)
     IF( a11 /= 0d0)THEN  !what if very low
       dl  = g1/a11
     ELSE
       dl = 0d0
     END IF
     dmu = disax(i,2)+dl*disax(i,1)
   END IF
   WRITE(58,"(2e12.4,4x,e12.4)") a11,a12,g1
   WRITE(58,"(2e12.4,4x,e12.4)") a21,a22,g2
   WRITE(58,"(3e12.4)")aux1,dl,dmu
   ! incremental displacements
   ddisp = dl*disax(1:neq,1) + disax(1:neq,2) + dmu*disax(1:neq,3)
   ! updates eigenvector and critical load
   disax(1:neq,4) =    dl*disax(1:neq,5) + disax(1:neq,6)          &
                     + dmu*disax(1:neq,7) + disax(1:neq,3)
   lambd = lambd + dl
   mu    = mu + dmu
   CALL screen(istep,itera,dl,lambd)
   !    evaluates total residual forces for this iteration (for checks only)
   IF(nload > 0 ) gvect = gvect + force(1:neq,n1)*dl
   deltw =      DOT_PRODUCT(ddisp,gvect)  ! energy in this iteration
   nordi = SQRT(DOT_PRODUCT(ddisp,ddisp)) ! norm of iterative displacements
   IF(ittol <= 2.AND.ABS(deltw/initw) > 1./toler)initw=deltw*toler
   IF(ittol == 3.AND.ABS(nordi/normd) > 1./toler)normd=nordi*toler
   !    updates coordinates, nodal local systems & step displacements
   CALL coract(ndime,npoin,neulr,ndofn,dl,coora,ddisp,coora,locsy)
   displ = displ + ddisp                  ! total incremental displacement
   !    computes residual forces
   CALL resvpl(istop,lambd,resid,gvect,.TRUE.)
   IF( istop == 1 )EXIT
   IF(nload > 0)THEN
     DO j=1,nload
       facts(j,1:2) = functs(loass(j),lambd)*force(neq+1,j)
     END DO
     gvect = MATMUL(force(1:neq,1:nload),facts(1:nload,1))+gvect
   END IF
   !    convergence checks
   normr = SQRT(DOT_PRODUCT(gvect,gvect))
   IF(ittol <= 2.AND.ABS(normr/normf) > 1./toler)normf=normr*toler
   SELECT CASE (ittol)
   CASE (1)
     IF(ABS(deltw/initw) < tol1 .AND. ABS(normr/normf) < SQRT(tol1) ) &
          convg = .TRUE.
   CASE (2)
     IF(ABS(normr/normf) < tol1 ) convg = .TRUE.
   CASE (3)
     IF(ABS(nordi/normd) < tol1 ) convg = .TRUE.
   END SELECT

   WRITE(lures,100)deltw,initw,normr,normf,nordi,normd
   WRITE(58,100)deltw,initw,normr,normf,nordi,normd
   IF(itera == niter.AND.(.NOT.convg)) THEN
     IF(diter > niter) EXIT
     IF(ittol == 2 .AND. ABS(deltw/initw) < tol1 ) THEN
       convg = .TRUE.
     ELSE
       WRITE(lures,200)istep,itera,lambd,dl
       STOP
     END IF
   END IF
 END DO
 WRITE(lures,"(' convergence',l5,' in ',i3,' iterations ')")convg,itera
 !     actualizes converged values
 coorc = coora           !updates converged configuration
 euler = locsy
 dlamb = lambd - oldlb
 CALL actual(displ,lambd,dlamb)
 CALL outbuc(npoin,ndime,istep,ifpre,disax(:,4),oldlb,lambd, 0 )
 nbuck = -nbuck
 ! poter = -lambd*DOT_PRODUCT(force(1:neq,1),tdisp)
 ! aux1 = tdisp(ecdis)
 ! WRITE(UNIT=58,ADVANCE='no',FMT="(' D=',e13.5,' l =',e12.4,' V=',  &
 !           &           e13.5)")aux1,lambd,poter

 RETURN
 100 FORMAT(' dw=',e9.3,' iw=',e9.3,' nr=',e9.3,' nf=',e9.3,           &
           &' dd=',e9.3,' id=',e9.3)
 200 FORMAT( //,' non-convergence in step =',i5,' in ',i3,' iterations'&
    &/,' actual load PARAMETER ',e12.4,' last increment in load',      &
    &e12.4,/,' PROGRAM  s t o p p e d ')

 END SUBROUTINE extend
