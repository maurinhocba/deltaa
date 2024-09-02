 SUBROUTINE j2smt9(nstre,nn,yistr,ishar,khard,f1,coefn,coefm,coenm,&
&                  ddmat,gkh,dmatx,strsg,h,pymat)
 !------------------------------------------------------------------------------
 !
 !     this routine compute the constitutive matrix for the j2 shell plasticity
 !
 !  input
 !     strsg: actual stress (n+1)
 !     h: (1,2=consistency parameters)
 !     dmatx: elastic constitutive matrix
 !  output
 !     dmatx: constitutive matrix
 !
 !------------------------------------------------------------------------------
 IMPLICIT NONE
 !***  parameters
 INTEGER (kind=4) nstre,nn
 REAL    (kind=8) dmatx(nstre,nstre),strsg(nstre),h(2),yistr,ishar,&
&                 khard,f1,coefn,coefm,coenm,ddmat(3),gkh,pymat(nn)
 !***  local variables
 LOGICAL ih
 INTEGER (kind=4) i,j,k,iyiel,jyiel,nq
 REAL    (kind=8) valu1,valu2,valu3,spcon,dpcon,thick,signo(2),    &
&                 sgefe(nstre),jmati(2,2),jmatd(2,2),a(4,4),b(4,4),&
&                 flows(nstre,2),nflow(nstre,2),ematp(nstre,nstre),&
&                 emasp(nstre,nstre)
 SAVE signo
 DATA signo /1d0, -1d0/
 !
 spcon = h(1) + h(2)
 dpcon = h(1) - h(2)
 IF (spcon == 0) RETURN
 thick = SQRT(1d0/coefn)

 !     compute the efective stresses  in diagonalized base

 nq = 2*nn+1
 sgefe = strsg
 CALL j2sdia(nn,sgefe)
 IF(khard > 0) THEN
   ih = .TRUE.
   ematp = 0d0
   emasp = 0d0
 ELSE
   ih = .FALSE.
 END IF

 !     compute algorithmic moduli

 dmatx = 0d0

 DO i=1,nn
    valu1 = spcon*coefn*pymat(i)
    valu2 = spcon*coefm*pymat(i)
    valu3 = dpcon*coenm*pymat(i)
 !           e(s,s)
    a(1,1) = 1.0d0/ddmat(i)    + valu1
    a(2,1) =                   + valu3
    a(1,2) =                   + valu3
    a(2,2) = 1.0d0/ddmat(i)/f1 + valu2
    IF (ih) THEN
 !           e(s,p)
      a(1,3) =                   - valu1
      a(1,4) =                   - valu3
      a(2,3) =                   - valu3
      a(2,4) =                   - valu2
 !           e(p,s)
      a(3,1) =                   - valu1
      a(3,2) =                   - valu3
      a(4,1) =                   - valu3
      a(4,2) =                   - valu2
 !           e(p,p)
      a(3,3) =  pymat(i)/khard   + valu1
      a(3,4) =                   + valu3
      a(4,3) =                   + valu3
      a(4,4) = pymat(i)/khard/f1 + valu2
      CALL invmtx(a,b,valu1,4)
    ELSE
      valu1  =  a(1,1)*a(2,2)-a(1,2)*a(2,1)
      b(1,1) =  a(2,2)/valu1
      b(1,2) = -a(1,2)/valu1
      b(2,1) = -a(2,1)/valu1
      b(2,2) =  a(1,1)/valu1
    END IF
    DO j=0,1
       DO k=0,1
          dmatx(i+nn*j,i+nn*k) = b(j+1,k+1)
          IF(ih) THEN
            ematp(i+nn*j,i+nn*k) = b(j+3,k+3)
            emasp(i+nn*j,i+nn*k) = b(j+1,k+3)
          END IF
       END DO
    END DO
 END DO
 valu1 = 2d0*spcon*coefn
 a(1,1)= 1.0d0/gkh   + valu1
 IF(ih) THEN
   a(2,1) =             - valu1
   a(3,1) =             - valu1
   a(4,1) = 1.0d0/khard + valu1
   CALL invmtx(a,b,valu1,2)
 ELSE
   b(1,1) = 1d0/a(1,1)
 END IF
 DO i=nq,nstre
    dmatx(i,i) = b(1,1)
    IF(ih) THEN
      emasp(i,i) = b(2,1)
      ematp(i,i) = b(4,1)
    END IF
 END DO

 !     compute actual jacobian matrix  "jmati"

 !         flow rule
 CALL j2sflw(nstre,nn,nq,flows,sgefe,pymat,coefn,coefm,coenm,      &
&            signo)
 nflow = 0d0
 DO i = 1,nstre
   DO j = 1,nstre
     IF(ih) THEN
       valu1 = dmatx(i,j) - 2d0*emasp(i,j) + ematp(i,j)
     ELSE
       valu1 = dmatx(i,j)
     END IF
     DO iyiel=1,2
       nflow(i,iyiel) = nflow(i,iyiel) + valu1*flows(j,iyiel)
     END DO
   END DO
 END DO
 jmati(1,1) = DOT_PRODUCT(flows(1:nstre,1),nflow(1:nstre,1))
 jmati(1,2) = DOT_PRODUCT(flows(1:nstre,1),nflow(1:nstre,2))
 jmati(2,1) = jmati(1,2)
 jmati(2,2) = DOT_PRODUCT(flows(1:nstre,2),nflow(1:nstre,2))
 !     INCLUDE the effect of isotropic hardening
 IF (ishar > 0) THEN
   valu1 = 4d0/9d0*yistr**2*ishar/thick
   jmati(1,1) = jmati(1,1)+valu1
   jmati(1,2) = jmati(1,2)+valu1
   jmati(2,1) = jmati(2,1)+valu1
   jmati(2,2) = jmati(2,2)+valu1
 END IF
 !
 valu1 = jmati(1,1)*jmati(2,2)-jmati(1,2)*jmati(2,1)
 IF (valu1 /= 0 .AND. h(1) > 0 .AND. h(2) > 0) THEN
   jmatd(1,1) =   jmati(2,2)/valu1
   jmatd(2,2) =   jmati(1,1)/valu1
   jmatd(1,2) = - jmati(1,2)/valu1
   jmatd(2,1) = - jmati(2,1)/valu1
 ELSE  !  select the case diagonal
   jmatd(1,1) = h(1)/(jmati(1,1)*h(1)+jmati(1,2)*h(2))
   jmatd(1,2) = 0.0d0
   jmatd(2,1) = 0.0d0
   jmatd(2,2) = h(2)/(jmati(2,1)*h(1)+jmati(2,2)*h(2))
 END IF

 !     compute elasto-platic matrix in diagonalized base

 IF(ih) THEN
   nflow = 0d0
   DO iyiel=1,2
     IF (h(iyiel) > 0) THEN
       DO i = 1,nstre
         DO j = 1,nstre
           nflow(i,iyiel) = nflow(i,iyiel)+                        &
&                          ( dmatx(j,i)-emasp(j,i) )*flows(j,iyiel)
         END DO
       END DO
     END IF
   END DO
 END IF

 DO iyiel=1,2
   IF(h(iyiel) > 0) THEN
     DO jyiel=1,2
       IF(h(jyiel) > 0) THEN
         DO i=1,nstre
           DO j=1,nstre
             dmatx(i,j) = dmatx(i,j) -                             &
&                nflow(i,iyiel)*jmatd(iyiel,jyiel)*nflow(j,jyiel)
           END DO
         END DO
       END IF
     END DO
   END IF
 END DO

 !     transform to non-diagonalized base

 DO i = 1,nstre
   CALL j2scar(nn,dmatx(1,i))
 END DO
 DO i = 1,nstre
   DO j = i+1,nstre
     valu1 = dmatx(i,j)
     dmatx(i,j) = dmatx(j,i)
     dmatx(j,i) = valu1
   END DO
 END DO
 DO i = 1,nstre
   CALL j2scar(nn,dmatx(1,i))
 END DO

 END SUBROUTINE j2smt9
