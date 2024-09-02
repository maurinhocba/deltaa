SUBROUTINE switch(istop,actio)
!********************************************************************
!
! *** Branch switching iterative routine
!
!********************************************************************
USE c_input
USE curv_db
USE npoi_db, ONLY : ndime,neulr,ndofn,npoin, &
                    coord,coora,coorc,coor1,euler,locsy,locs1
USE solv_db
USE kin0_db, ONLY : ndepd,naris
USE kin1_db, ONLY : neq,maxa,maxav
USE kin2_db, ONLY : nvelr,ifpre,velor,nvfix
USE load_db, ONLY : nload,force,loadv,loass
IMPLICIT NONE
INTEGER (kind=4), INTENT(IN OUT) :: istop
CHARACTER(len=*),INTENT(IN):: actio

LOGICAL :: convg, error
INTEGER(kind=4) :: itera,i,jj,neq1,ecphi(1),karcl,ittol,j,n, &
                   n1,nld,l,nsu,niter
REAL (kind=8) :: initw,deltw,normf,normr,dl,normd,                &
                 nordi,oldlb,fphi,a,b,c,d,e,aux1,eps,l1,l2,       &
                 facts(nload,2)

INTERFACE
  INCLUDE 'timing.h'
  INCLUDE 'coract.h'
  INCLUDE 'resvpl.h'
  !INCLUDE 'screen.h'
  INCLUDE 'stiffm.h'
  INCLUDE 'colsol.h'
  INCLUDE 'update.h'
  INCLUDE 'mulmat.h'
  INCLUDE 'outbuc.h'
  INCLUDE 'actual.h'
  INCLUDE 'functs.h'
END INTERFACE

error = .FALSE.
neq1  = neq+1
maxa  = maxav(neq1)               !size of stiffness matrix
jj    = 2 +ABS(nsymm)             !first position in DSTIF
nsu   = 1 +ABS(nsymm)             !position of the asymmetric part
n1    = nload + 1                 !position of the load vector
nld   = n1
IF(nvelr > 0 ) nld = n1+1         !position of the load vector

!     C O M P U T E   A S Y M P O T I C   D E R I V A T I V E S
!     most relevant DOF in eigenvector  ==> i
ecphi = MAXLOC(disax(1:neq,4))    !maximum value = +1
i = ecphi(1)
!     correct tangent vector so that q(1)i = 1d0  ( y(i) = 0d0 )
disax(1:neq,1) = disax(1:neq,1)-disax(i,1)*disax(1:neq,4)  !
!     computes epsilon (increment for numerical derivative of residual)
eps   = SQRT(alph1)
!     proyection of critical mode over total displacements and force vector
IF(nload > 0)THEN
  disax(1:neq,5) = force(1:neq,n1)
ELSE
  disax(1:neq,5) = 0d0
END IF
IF(nvelr > 0) disax(1:neq,5) = disax(1:neq,5) + force(1:neq,nld)
normd = SQRT( DOT_PRODUCT( disax(1:neq,4),disax(1:neq,4)) )
coor1 = coora - coord                            !total displacements
aux1  = 0d0                                      !initializes
nordi = 0d0                                      !initializes
DO n=1,npoin                                     !for every node
  DO j=1,ndime                                   !for every direction in space
    l = ifpre(j,n)                               !active DOF
    IF( l > 0 )THEN                              !if an active DOF
      aux1 = aux1 + coor1(j,n) * disax(l,5)      ! TDISP * FORCE
      nordi= nordi+ coor1(j,n) * coor1(j,n)      ! TDISP * TDISP
    END IF
  END DO
END DO
fphi  = DOT_PRODUCT( disax(1:neq,4),disax(1:neq,5) )   ! PHI * FORCE
IF(nordi > 0d0 .AND. aux1 /= 0d0)THEN
  aux1  = fphi*nordi/aux1/normd
ELSE
  aux1  = fphi/normd
END IF
WRITE(lures,"(' normalized coefficient F.Phi=',e12.4)")aux1
!     notice that in STIFF(1,jj) we have  - Du(Kt.Phi) == - Vijk Xk
CALL mulmat(stiff(:,jj),maxav,disax(:,2),disax(:,4),neq)
c  = - DOT_PRODUCT(disax(1:neq,2),disax(1:neq,4)) !Vijk Xi Xj Xk
!     results will be stored as  u(1) = disax(4)  , u(2) = disax(5)
IF(ABS(aux1) > 1d-2 .AND. .NOT.linear )THEN          ! Limit point
  WRITE(lures,"(//' LIMIT POINT ')")
  l1 = 0d0             !First load derivative
  l2 = c/fphi          !Second Load derivative
  CALL colsol(stiff(:,1),maxav,neq,2,55,1,0,v=disax(:,2),u=stiff(:,1))
  disax(1:neq,5) = disax(1:neq,2) + l2*disax(1:neq,1)
ELSE                              ! Bifurcation
  b = - DOT_PRODUCT(disax(1:neq,2),disax(1:neq,1))
  CALL mulmat(stiff(:,jj),maxav,disax(:,6),disax(:,1),neq)
  a = - DOT_PRODUCT(disax(1:neq,6),disax(1:neq,1))
  IF( b /= 0)THEN
    aux1 =c/b*SQRT(DOT_PRODUCT(disax(1:neq,1),disax(1:neq,1)))/normd
  ELSE
    b = -0.5d0
    aux1 = 1d0
  END IF
  WRITE(lures,"(//' BIFURCATION'/' a,b,c c/b coefficients =',     &
               &   4e12.4)")       a,b,c,aux1
  IF( a /= 0d0 )THEN        !quadratic term not NULL
    d  = b*b-a*c            !discriminant of the quadratic equation
    IF( d < 0 )THEN         !if roots are imaginary
      WRITE(lures,"(' error in quadratic equation, d=',e15.4)")d
      l1   = 0d0            !set as symmetric
      l2   = 1d0            !any value is O.K. here
      aux1 = 0d0
      error = .TRUE.        !flag not to compute the second derivative
    ELSE
      d = SQRT(d)           !root of the discriminant (radius)
      IF( b > 0 ) THEN      !if B coefficient is positive
        l1 = (-b+d)/a       !smallest root
        l2 = (-b-d)/a       !biggest root
      ELSE
        l2 = (-b+d)/a       !biggest root
        l1 = (-b-d)/a       !smallest root
      END IF
    END IF
  ELSE
    l1 = 0d0                !first root
    l2 = -c/(2d0*b)         !root of the linear remaining problem
  END IF
  WRITE(lures,"(' Roots =',2e12.4)")l1,l2
  IF(( ABS(aux1) > 1d-3 .OR.  ABS(l1/l2) > 1d-5) .AND. .NOT.linear ) THEN    ! Asymmetric
    WRITE(lures,"(//' ASYMMETRIC BIFURCATION')")
    ! u(1) = disax(4) + l1 * disax(1)
    IF( l1 == 0d0) l1=l2
    disax(1:neq,4) = disax(1:neq,4) + l1* disax(1:neq,1)
    ! no values for second order derivatives yet
    l2 = 0d0
    disax(1:neq,5) = 0d0
  ELSE                                        ! Symmetric Bifurcation
    WRITE(lures,"(//' SYMMETRIC BIFURCATION')")
    l1 = 0d0
    IF(nload > 0)THEN
      DO j=1,nload
       facts(j,1:2) = functs(loass(j),lambd)*force(neq+1,j)
      END DO
      disax(1:neq,7) = MATMUL(force(1:neq,1:nload),facts(1:nload,1))+ddisp
    END IF
    ! First computes Vijkl Xi Xj Xk Xl as a numerical derivative
    ddisp = -2d0*eps*disax(1:neq,4)  ! (v - 2eps * phi)
    locs1 = locsy
    CALL coract(ndime,npoin,neulr,ndofn,0d0,coora,ddisp,coor1,locs1)
    ! computes stresses, plastic consistency pars. and residual forces
    CALL resvpl(istop,lambd,gvect,ddisp)
    IF(nload > 0) ddisp = ddisp + disax(1:neq,7)

    d = + DOT_PRODUCT(ddisp,disax(1:neq,4))/2d0
    ddisp = -eps*disax(1:neq,4)      ! (v - eps * phi)
    locs1 = locsy
    CALL coract(ndime,npoin,neulr,ndofn,0d0,coora,ddisp,coor1,locs1)
    ! computes stresses, plastic consistency pars. and residual forces
    CALL resvpl(istop,lambd,gvect,ddisp)
    IF(nload > 0) ddisp = ddisp + disax(1:neq,7)
    d = d - DOT_PRODUCT(ddisp,disax(1:neq,4))
    ddisp = eps*disax(1:neq,4)       ! (v + eps * phi)
    locs1 = locsy
    CALL coract(ndime,npoin,neulr,ndofn,0d0,coora,ddisp,coor1,locs1)
    ! computes stresses, plastic consistency pars. and residual forces
    CALL resvpl(istop,lambd,gvect,ddisp)
    IF(nload > 0) ddisp = ddisp + disax(1:neq,7)
    d = d + DOT_PRODUCT(ddisp,disax(1:neq,4))
    ddisp = 2d0*eps*disax(1:neq,4)   ! (v + 2*eps * phi)
    locs1 = locsy
    CALL coract(ndime,npoin,neulr,ndofn,0d0,coora,ddisp,coor1,locs1)
    ! computes stresses, plastic consistency pars. and residual forces
    CALL resvpl(istop,lambd,gvect,ddisp)
    IF(nload > 0) ddisp = ddisp + disax(1:neq,7)
    d = (d - DOT_PRODUCT(ddisp,disax(1:neq,4))/2d0)/eps**3
    WRITE(lures,"(' Vijkl Xi Xj Xk Xl=Phi.Du[Du(Kt.Phi).Phi].Phi',&
              e12.4)")d
    ! Compute Vijk Xi Xj Zk = Phi.Du(Kt.Phi).Z  as a numerical derivative
    disax(1:neq,7) = disax(1:neq,2)
    CALL colsol(stiff(:,1),maxav,neq,2,55,1,0,v=disax(:,7))
    e = - DOT_PRODUCT(disax(1:neq,2),disax(1:neq,7))
    WRITE(lures,"(' Vijk Xi Xj Zk = Phi.Du(Kt.Phi).Z ',e12.4)")e
    IF( error .OR. linear )THEN                   ! if D < 0
      l2 = 0d0                     !set second derivative to 0
    ELSE
      l2 = -( d + 3d0 * e)/3d0/b
    END IF
    !         correction to account for U(2)i = 0
    ! aux1 = l2*disax(i,1) + disax(i,7)
    disax(1:neq,5) =   l2*disax(1:neq,1) + disax(1:neq,7)       ! &
    !                  - aux1*disax(1:neq,4)
  END IF
END IF
WRITE(lures,"(/,' Load derivatives, l1, l2 ',2e15.5)")l1,l2
WRITE(lures,"(' u(2)i=',e12.4,' u(2)M=',e12.4,' u(2)m=',e12.4)")  &
      &    disax(i,5),MAXVAL(disax(1:neq,5)),MINVAL(disax(1:neq,5))
CALL outbuc(npoin,ndofn,istep,ifpre,disax(:,4),lambd,lambd, 1)

CALL outbuc(npoin,ndofn,istep,ifpre,disax(:,5),lambd,lambd,2)

!                  S W I T C H I N G   P R O C E D U R E
oldlb = lambd                     !keep critical load
karcl = 4                         !adopt prescribed displacement
ittol = MOD(nralg,1000)/100       !convergence type
WRITE(lures,"(' Most relevant DOF ',i5)")i

!      read control data for Switching

IF (arcl1 /= 0d0 ) arcln = arcl1
IF( scdis > 0  )THEN
  j = scdis / 10
  jj= MOD( scdis, 10 )
  jj = ifpre(jj,j)      !Control DOF
  ncdis = scdis
ELSE             !use largest DOF in eigenvalue
  jj = i              !keep DOF
  !      search asociated NODE and DOF
  ncdis = 0           !initializes to not found
  outer : DO n=1,npoin        !search for the node and DOF
    DO j=1,ndofn              !search for control DOF
      IF(ifpre(j,n) == i)THEN             !DOF found
        IF( j <= ndime ) ncdis = j+10*n   !if a translational DOF OK.
      EXIT outer                        !exit search
      END IF
    END DO
  END DO outer
END IF
arcl1 = disax(i,4)/disax(jj,4)*arcln    !scale arc-length
i = jj                           !control DOF

!     prediction values for incremental displacements and load factor
aux1  = arcl1**2/2d0
displ = disax(1:neq,4)*arcl1 + disax(1:neq,5)*aux1
dl    =  l1*arcl1 + l2*aux1
lambd = lambd + dl
disax(1:neq,2) = arcl1*disax(1:neq,4)  !comparison vector
!    updates coordinates and nodal local systems and computes residual forces
CALL coract(ndime,npoin,neulr,ndofn,dl,coora,displ,coora,locsy)
CALL resvpl(istop,lambd,gvect,resid,.TRUE.)

IF((ndepd+naris) > 0 )THEN !regenerate force vector
  IF(nload > 0) force(1:neq,1:nload) = 0d0
  DO j=1,nload
    CALL ensvec (ndofn*npoin,ifpre(1,1),loadv(1,1,j),force(1,j))
  END DO
END IF
IF(nload > 0)THEN
  DO j=1,nload
    facts(j,1:2) = functs(loass(j),lambd)*force(neq+1,j)
  END DO
  force(1:neq,n1)=MATMUL(force(1:neq,1:nload),facts(1:nload,2))
  resid = MATMUL(force(1:neq,1:nload),facts(1:nload,1))+resid
END IF
!                    Initializes values for iteration
niter = nralg/10000000                                  !two digits
convg = .FALSE.
itera = 1
initw = MAX( ABS(DOT_PRODUCT(displ,resid)),toler)
normf = MAX(SQRT(DOT_PRODUCT(resid,resid)),toler)
normd = MAX(SQRT(DOT_PRODUCT(displ,displ)),toler)
WRITE(lures,100)initw,initw,normf,normf,normd,normd
CALL screen(istep,itera,dl,lambd)
!                  begin iteration loop
DO
  itera = itera + 1
  !       evaluates stiffness matrix and force vector for presc. displac
  CALL stiffm(istop,lambd,.TRUE.,nld,stiff(:,1),force,nsymm, &
              stiff(:,nsu))

  !       L D U factorization of stiffness matrix
  CALL timing(5,.TRUE.)
  CALL colsol(stiff(:,1),maxav,neq,1,55,1,nsymm)
  CALL timing(5,.FALSE.)
  !       evaluates tangent solution vector due to external forces and b.c.
  IF(nload > 0)THEN
    DO j=1,nload
      facts(j,1:2) = functs(loass(j),lambd)*force(neq+1,j)
    END DO
   force(1:neq,n1)=MATMUL(force(1:neq,1:nload),facts(1:nload,2))
   disax(1:neq,1) = force(1:neq,n1)
  ELSE
    disax(1:neq,1) = 0d0
  END IF
  ! force vector due to prescribed displacements
  IF(nvelr > 0) disax(1:neq,1) = disax(1:neq,1) + force(1:neq,nld)
  CALL timing(5,.TRUE.)
  CALL colsol(stiff(:,1),maxav,neq,2,55,1,nsymm,v=disax(:,1),u=stiff(:,nsu))
  CALL timing(5,.FALSE.)
  !       evaluates incremental displacements due to residual forces
  ddisp = resid
  CALL timing(5,.TRUE.)
  CALL colsol(stiff(:,1),maxav,neq,2,4,1,nsymm,v=ddisp,u=stiff(:,1))
  CALL timing(5,.FALSE.)
  !       updates load factor and incre. displ. considering selected algorithm
  CALL update(istep,itera,  0  ,  i  ,neq,arcln,lambd,dlamb,      &
              disax,displ,ddisp,dl,karcl, 1d0 , 1d0 ,  0  ,ncdis)
  CALL screen(istep,itera,dl,lambd)
  !       evaluates total residual forces for this iteration
  IF(nload > 0 ) resid = resid + force(1:neq,n1)*dl
  deltw = DOT_PRODUCT(ddisp,resid)
  nordi = SQRT(DOT_PRODUCT(ddisp,ddisp))
  IF(ittol <= 2.AND.ABS(deltw/initw) > 1./toler) initw = deltw*toler
  IF(ittol <= 2.AND.ABS(normr/normf) > 1./toler) normf = normr*toler
  IF(ittol == 3.AND.ABS(nordi/normd) > 1./toler) normd = nordi*toler
  !       updates coordinates and nodal local systems
  CALL coract(ndime,npoin,neulr,ndofn,dl,coora,ddisp,coora,locsy)
  ! updates step displacements
  displ = displ + ddisp
  ! WRITE(55,"(i7,3e15.6)")(i,coora(1:3,i)-coorc(1:3,i),i=1,npoin)
  CALL resvpl(istop,lambd,gvect,resid,.TRUE.)
  IF( istop == 1 )EXIT
  IF(nload > 0)THEN
    DO j=1,nload
      facts(j,1:2) = functs(loass(j),lambd)*force(neq+1,j)
    END DO
    resid = MATMUL(force(1:neq,1:nload),facts(1:nload,1))+resid
  END IF
  !       convergence checks
  normr = SQRT(DOT_PRODUCT(resid,resid))
  SELECT CASE (ittol)
  CASE (1)
    IF( ABS(deltw/initw) < toler .AND. ABS(normr/normf) < SQRT(toler)) &
      convg = .TRUE.
  CASE (2)
    IF( ABS(normr/normf) < toler) convg = .TRUE.
  CASE (3)
    IF( ABS(nordi/normd) < toler) convg = .TRUE.
  END SELECT
  WRITE(lures,100)deltw,initw,normr,normf,nordi,normd
  !       final decision about convergence
  IF(itera == niter.AND.(.NOT.convg)) THEN
    IF(ittol == 2 .AND. ABS(deltw/initw) < toler) THEN
      convg = .TRUE.
    ELSE
      !IF(niter > diter) THEN
        WRITE(lures,200)istep,itera,lambd,dl
        STOP
      !ELSE
      !  convg = .TRUE.
      !END IF
    END IF
  END IF
  IF( convg )EXIT
END DO
!     actualize converged values
WRITE(lures,"(' convergence in ',i3,' iterations ')")itera
!     tdisp = tdisp + displ
coorc = coora           !updates converged configuration
euler = locsy
dlamb = lambd - oldlb
CALL actual(displ,lambd,dlamb)
karcl = MOD(nralg,100000)/10000
!     modify advance parameters for next standard step
SELECT CASE (karcl)
CASE (0)
   nralg = nralg+20000   !updated tangent plane
   arcln = normd/2d0
CASE (1:3)
  arcln = normd/2d0
CASE (4)
  arcln = arcl1/2d0
CASE (5)
  !!!!!!!!! no sense bifurcation with tools
END SELECT
nbuck = 0
disax(1:neq+1,5:7) = 0d0
!     poter = -lambd*DOT_PRODUCT(force(1:neq,1),tdisp)
ecdis = i           !position with largest value in eigenvector or selected DOF
IF( ncdis == 0 )THEN  !IF not a Translational DOF search the largest
  aux1 = 0d0          !initializes the largest displacement increment
  DO n=1,npoin        !search for the node and DOF
    DO j=1,ndime
      i = ifpre(j,n)            !associated DOF
      IF(i > 0)THEN             !If an active DOF
        IF( ABS(displ(i)) >= aux1 ) THEN     !if displacement is greater
          aux1 = ABS(displ(i))               !updates maximum Translational DOF
          ncdis = 10*n+j                     !updates node and DOF
          ecdis = i                          !updates ACTIVE DOF
        END IF
      END IF
    END DO
  END DO
  IF( karcl == 4 ) arcln = displ(ecdis)/2d0    !new ARC-LENGTH
END IF

! aux1 = tdisp(ecdis)
! WRITE(UNIT=4,ADVANCE='no',FMT="(' D=',e13.5,' l =',e12.4,' V=',   &
! &             e13.5)")aux1,lambd,poter

RETURN
100 FORMAT(' dw=',e9.3,' iw=',e9.3,' nr=',e9.3,' nf=',e9.3,           &
    &       ' dd=',e9.3,' id=',e9.3)
200 FORMAT( //,' non-convergence in step =',i5,' in ',i3,' iterations'&
    &/,' actual load PARAMETER ',e12.4,' last increment in load',      &
    &e12.4,/,' PROGRAM  s t o p p e d ')

END SUBROUTINE switch
