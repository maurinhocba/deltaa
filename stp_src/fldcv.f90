SUBROUTINE CalDisFLD(npoif,fldcur,xp,yp,fld)
!*****************************************************************************
!     Calculate the distant to FLD curve if point is out of the
!     admisible state. In other case, the distant is cero.
!*****************************************************************************
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: npoif
  REAL(kind=8),INTENT(IN):: fldcur(2,npoif), xp, yp
  REAL(kind=8),INTENT(OUT):: fld
  !--- Local variables
  INTEGER(kind=4):: nsegm, isegm, ipoif, nsigno
  INTEGER(kind=4),ALLOCATABLE:: nyorc(:)
  REAL(kind=8):: dref, xdif, ydif, xnorm, xtang, dmin, d
  REAL(kind=8),ALLOCATABLE:: c1(:), c2(:), xlist(:)

  nsegm = npoif - 1
  ALLOCATE(nyorc(nsegm),xlist(nsegm),c1(nsegm),c2(nsegm))

  !**** CALCULO DE C1(ISEGM);C2(ISEGM);NYORC(ISEGM) PARA CADA SEGMENTO DE LA FLD
  !
  DO isegm=1,nsegm
    ipoif = isegm
    xdif = fldcur(1,ipoif+1) - fldcur(1,ipoif)
    ydif = fldcur(2,ipoif+1) - fldcur(2,ipoif)
    !
    !**** SI EL PUNTO CONSECUTIVO ESTA MUY PROXIMO AL ANTERIOR SE PARA
    !
    IF ((DABS(ydif).LE.0.0001) .AND. (DABS(xdif).LE.0.0001)) THEN
      WRITE(*,"(' ERROR: FLD CURVE WRONG DESCRIBED: 2 POINTS CLOSE ',2I5)") ipoif,ipoif+1
      STOP
    !
    !**** CALCULO DE NYORC
    !
    ELSE
      nyorc(isegm) = 3
      IF (DABS(ydif).LE.0.0001) nyorc(isegm)=2
      IF (DABS(xdif).LE.0.0001 .AND. ydif.GT.0.0001)  nyorc(isegm)=1
      IF (DABS(xdif).LE.0.0001 .AND. ydif.LT.-0.0001) nyorc(isegm)=4
    END IF

    SELECT CASE (nyorc(isegm))
    CASE (1,4)
      !**** PARA SEGMENTOS CUASIVERTICALES C1,C2 SE DETERMINAN EN LAS ORDENADAS
      xnorm = -xdif/ydif
      xlist(isegm) = xnorm
      c1(isegm) = fldcur(2,ipoif)-xlist(isegm)*fldcur(1,ipoif)
      c2(isegm) = fldcur(2,ipoif+1)-xlist(isegm)*fldcur(1,ipoif+1)
    CASE (2)
      !**** PARA SEGMENTOS CUASI HORIZONTALES C1,C2 SE DETERMINAN EN "X" CON XTANG
      xtang = ydif/xdif
      xlist(isegm) = xtang
      c1(isegm) = fldcur(1,ipoif) + xlist(isegm)*fldcur(2,ipoif)
      c2(isegm) = fldcur(1,ipoif+1) + xlist(isegm)*fldcur(2,ipoif+1)
    CASE (3)
      !**** PARA SEGMENTOS EN LA ZONA 3 SE UTILIZA XNORM#0
      xnorm = -xdif/ydif
      xlist(isegm) = xnorm
      c1(isegm) = fldcur(1,ipoif) - fldcur(2,ipoif)/xlist(isegm)
      c2(isegm) = fldcur(1,ipoif+1) - fldcur(2,ipoif+1)/xlist(isegm)
    END SELECT
  END DO
  !
  !**** LOOKING FOR DMAX(MAXIMUM AMONG MINIMUM) FOR A LINEAR TRANSFORMATION
  !
  CALL CalDisCurv(0d0,0d0,npoif,nsegm,nyorc,fldcur,c1,c2,xlist,nsigno,dref)
  !
  !**** CALCULO DE LA DISTANCIA A LA FLC
  !
  CALL CalDisCurv(xp,yp,npoif,nsegm,nyorc,fldcur,c1,c2,xlist,nsigno,dmin)
  !
  !**** LINEAR TRANSFORMATION
  !
  IF (nsigno.EQ.0) THEN
    d  = - dmin/dref
  ELSE
    d  = dmin/dref
  END IF
  fld = 1d2*(1d0+d)

  DEALLOCATE (nyorc,xlist,c1,c2)

RETURN
END SUBROUTINE CalDisFLD

SUBROUTINE CalDisCurv(xp,yp,npoif,nsegm,nyorc,fldcur,c1,c2,xlist,nsigno,dmin)
!*****************************************************************************
!     Calculate the distant to a curve defined by segments if point is
!     out of the admisible state. In other case, the distant is cero.
!*****************************************************************************
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: npoif, nsegm, nyorc(nsegm)
  INTEGER(kind=4),INTENT(OUT):: nsigno
  REAL(kind=8),INTENT(IN):: xp, yp, fldcur(2,npoif), c1(nsegm), c2(nsegm), xlist(nsegm)
  REAL(kind=8),INTENT(OUT):: dmin
  !--- Local variables
  INTEGER(kind=4):: isegm, ipoif, msign, msegm
  REAL(kind=8):: dt, fpx, fpy
  REAL(kind=8),ALLOCATABLE:: cp(:)

  ALLOCATE (cp(nsegm))

  dmin = HUGE(1d0)
!  dmin=10**4
  msign = -1
  !
  !**** CP(ISEGM) DETERMINATION
  DO isegm=1,nsegm
    ipoif=isegm
    SELECT CASE (nyorc(isegm))
      CASE (1)
        cp(isegm)=yp-xlist(isegm)*xp
      CASE (2)
        cp(isegm)=xp+xlist(isegm)*yp
      CASE (3)
        cp(isegm)=xp-yp/xlist(isegm)
      CASE (4)
        cp(isegm)=yp-xlist(isegm)*xp
    END SELECT
    !
    !**** RECTANGULAR ZONES ALGORITHM
    !
    IF(nyorc(isegm).LE.3)THEN
      IF(c1(isegm).LE.cp(isegm).AND.cp(isegm).LE.c2(isegm))THEN
        CALL distseg(dt,npoif,isegm,xp,yp,fldcur,nyorc(isegm))
        CALL signo(nsigno,npoif,isegm,xp,yp,fldcur)
        IF (dt.LE.dmin) THEN
          dmin = dt
          msign = nsigno
        END IF
      END IF
    ELSE
      IF (c2(isegm).LE.cp(isegm) .AND. cp(isegm).LE.c1(isegm)) THEN
        CALL distseg(dt,npoif,isegm,xp,yp,fldcur,nyorc(isegm))
        CALL signo(nsigno,npoif,isegm,xp,yp,fldcur)
        IF (dt.LE.dmin) THEN
          dmin = dt
          msign = nsigno
        END IF
      END IF
    END IF
  END DO
  !
  !**** TRIANGULAR ZONES ALGORITHM FOR POINTS WHOSE DISTANCE
  !     TO INTERMEDIATE POINTS IS SHORTER
  !
  msegm=nsegm-1
  DO isegm=1,msegm
    ipoif=isegm
    IF (nyorc(isegm).LE.3 .AND. nyorc(isegm+1).LE.3) THEN
      IF (c2(isegm).LE.cp(isegm) .AND. cp(isegm+1).LE.c1(isegm+1)) THEN
        fpx=xp-fldcur(1,ipoif+1)
        fpy=yp-fldcur(2,ipoif+1)
        dt=SQRT(fpx**2+fpy**2)
        nsigno=0
        IF (fpy.GE.0.0) nsigno=1
        IF (dt.LE.dmin) THEN
          dmin = dt
          msign = nsigno
        END IF
      END IF
    END IF
    IF (nyorc(isegm).LE.3 .AND. nyorc(isegm+1).EQ.4) THEN
      IF (c2(isegm).LE.cp(isegm) .AND. cp(isegm+1).GE.c1(isegm+1)) THEN
        fpx=xp-fldcur(1,ipoif+1)
        fpy=yp-fldcur(2,ipoif+1)
        dt=SQRT(fpx**2+fpy**2)
        nsigno=1
        IF (dt.LE.dmin) THEN
          dmin = dt
          msign = nsigno
        END IF
      END IF
    END IF
    IF (nyorc(isegm).EQ.4 .AND. nyorc(isegm+1).LE.3) THEN
      IF (c2(isegm).GE.cp(isegm) .AND. cp(isegm+1).LE.c1(isegm+1)) THEN
        fpx=xp-fldcur(1,ipoif+1)
        fpy=yp-fldcur(2,ipoif+1)
        dt=SQRT(fpx**2+fpy**2)
        nsigno=0
        IF (dt.LE.dmin) THEN
          dmin = dt
          msign = nsigno
        END IF
      END IF
    END IF
    IF (nyorc(isegm).EQ.4 .AND. nyorc(isegm+1).EQ.4) THEN
      IF (c2(isegm).GE.cp(isegm) .AND. cp(isegm+1).GE.c1(isegm+1)) THEN
        fpx=xp-fldcur(1,ipoif+1)
        fpy=yp-fldcur(2,ipoif+1)
        dt=SQRT(fpx**2+fpy**2)
        nsigno=0
        IF (dt.LE.dmin) THEN
          dmin = dt
          msign = nsigno
        END IF
      END IF
    END IF
  END DO

  DO isegm=1,nsegm
    IF (isegm.EQ.1 .OR. isegm.EQ.nsegm) THEN
      IF (isegm.EQ.1) ipoif=1
      IF(isegm.EQ.nsegm) ipoif=npoif
      fpx=xp-fldcur(1,ipoif)
      fpy=yp-fldcur(2,ipoif)
      dt=SQRT(fpx**2+fpy**2)
      IF (dt.LT.dmin) THEN
        dmin=dt
        CALL distseg(dt,npoif,isegm,xp,yp,fldcur,nyorc(isegm))
        CALL signo(nsigno,npoif,isegm,xp,yp,fldcur)
        IF (dt.LE.dmin) THEN
          dmin = dt
          msign = nsigno
        END IF
      END IF
    END IF
  END DO

  DEALLOCATE (cp)

  nsigno = msign

RETURN
END SUBROUTINE CalDisCurv

SUBROUTINE distseg(dt,npoif,isegm,xp,yp,fldcur,nyork)
!*****************************************************************************
!**** DISTANCE FROM POINT P TO IANY SEGMENT
!*****************************************************************************
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: npoif, isegm, nyork
  REAL(kind=8),INTENT(IN):: xp, yp, fldcur(2,npoif)
  REAL(kind=8),INTENT(OUT):: dt
  !--- Local variables
  INTEGER(kind=4):: ipoif
  REAL(kind=8):: xdif, ydif, xtang, xnorm, eval, fval, xt, yt

  !
  !***  DISTANCE SIMPLIFIED FOR CASES 1,2,4
  !
  ipoif=isegm
  SELECT CASE (nyork)
    CASE (1)
      dt=DABS(xp-fldcur(1,ipoif))
    CASE (2)
      dt=DABS(yp-fldcur(2,ipoif))
    CASE (3)
      xdif=fldcur(1,ipoif+1)-fldcur(1,ipoif)
      ydif=fldcur(2,ipoif+1)-fldcur(2,ipoif)
      xtang=ydif/xdif
      xnorm=-1/xtang
      eval=yp-xnorm*xp-fldcur(2,ipoif)+xtang*fldcur(1,ipoif)
      fval=xtang-xnorm
      xt=eval/fval
      yt=yp+xnorm*(xt-xp)
      dt=SQRT((xt-xp)**2+(yt-yp)**2)
    CASE (4)
      dt=DABS(xp-fldcur(1,ipoif))
    END SELECT

RETURN
END SUBROUTINE distseg

SUBROUTINE signo(msign,npoif,isegm,xp,yp,fldcur)
!*****************************************************************************
!     VECTORIAL PRODUCT IF POINT REMAINS TO THE SEGMENT'S RIGHT.
!     THEN, MSIGN=0 (BREAKOUT)
!*****************************************************************************
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: npoif
  INTEGER(kind=4),INTENT(OUT):: msign
  REAL(kind=8),INTENT(IN):: xp, yp, fldcur(2,npoif)
  !--- Local variables
  INTEGER(kind=4):: isegm, ipoif
  REAL(kind=8):: xdif, ydif, xpx, ypy, cval

  !
  !**** VECTORIAL PRODUCT OF (XDIF,YDIF,0)^(XPX,YPY,0)
  !
  ipoif=isegm
  xdif=fldcur(1,ipoif+1)-fldcur(1,ipoif)
  ydif=fldcur(2,ipoif+1)-fldcur(2,ipoif)
  xpx=xp-fldcur(1,ipoif)
  ypy=yp-fldcur(2,ipoif)
  cval=xdif*ypy-xpx*ydif

  msign=1
  IF (cval.LT.0.0) msign=0

RETURN
END SUBROUTINE signo

SUBROUTINE FLDdg(LmM,LmLS,LmPS,LmT,xp,yp,fld,SafZn,ForZn)
!*****************************************************************************
!  This subroutine calculates 'Safety Zone' diagram and 'Forming Zone' diagram.
!*****************************************************************************
IMPLICIT NONE

  !Dummy variables
  REAL(kind=8),INTENT(IN):: LmM, LmLS, LmPS, LmT, xp, yp, fld
  REAL(kind=8),INTENT(OUT):: SafZn,ForZn
  !Local Variables
  REAL(kind=8):: pi, AngBS, AngPSi, AngPSs, AngUT, AngPS, AngUC, LimThk, radio, theta

  pi = 4d0*DATAN(1d0)

  !Limits values
  AngBS = 0.25d0*pi          !Biaxial Stretch Angle
  AngPSi = 0.5d0*pi - LmPS   !Plane Strain Upper Angle
  AngPSs = 0.5d0*pi + LmPS   !Plane Strain Lower Angle
  AngUT = pi - DATAN(2d0)    !Uniaxial Tension Angle
  AngPS = 0.75d0*pi          !Pure Shear Angle
  AngUC = pi - DATAN(0.5d0)  !Uniaxial Compresion Angle
  LimThk = DLOG(1d0-LmT)     !Excessive thinning in logaritmic strains

  CALL CarToPol(xp,yp,radio,theta)

  !*** SAFETY ZONE DIAGRAM ***
  IF (radio <= LmLS) THEN          !LOW STRAIN
    SafZn = 3d0
  ELSE IF (fld > 1d2) THEN         !FAIL
    SafZn = 7d0
  ELSE IF (xp+yp > -LimThk) THEN   !EXCESS. THINNING
    SafZn = 6d0
  ELSE IF (fld > 1d2-LmM) THEN     !MARGINAL
    SafZn = 5d0
  ELSE IF (theta >= AngUC) THEN    !STRONG WRINKLING
    SafZn = 1d0
  ELSE IF (theta >= AngPS) THEN    !WRINKLING
    SafZn = 2d0
  ELSE                             !SAFE
    SafZn = 4d0
  END IF

  !*** FORMING ZONE DIAGRAM ***
  IF (radio <= LmLS) THEN         !LOW STRAIN
    ForZn = 5d0
  ELSE IF (theta >= AngUC) THEN   !STRONG WRINKLING
    ForZn = 7d0
  ELSE IF (theta >= AngPS) THEN   !WRINKLING
    ForZn = 6d0
  ELSE IF (theta >= AngUT) THEN   !LOOSE
    ForZn = 4d0
  ELSE IF (theta >= AngPSs) THEN  !SEMI TIGHT
    ForZn = 3d0
  ELSE IF (theta >= AngPSi) THEN  !PLANE STRAIN
    ForZn = 2d0
  ELSE IF (theta >= AngBS) THEN   !TIGHT
    ForZn = 1d0
  END IF

RETURN
END SUBROUTINE FLDdg

SUBROUTINE CarToPol(x,y,radio,theta)
!*****************************************************************************
!  This subroutine transform cartesian coordinates to polar coordiantes
!*****************************************************************************
IMPLICIT NONE

  !Dummy variables
  REAL (kind=8) :: x, y, radio, theta
  !Local Variables
  REAL (kind=8) :: pi

  pi = 4d0 * ATAN(1.0d0)

  radio = SQRT(x*x + y*y)

  IF (radio > 0d0) THEN
    IF (x /= 0d0) THEN
      theta = ATAN(ABS(y/x))
    ELSE IF (y > 0d0) THEN
      theta = 0.5d0*pi
    ELSE
      theta = 1.5d0*pi
    END IF

    IF ( (x < 0d0) .AND. (y >= 0d0) ) THEN
      theta = pi - theta
    ELSE IF ( (x < 0d0) .AND. (y < 0d0) ) THEN
      theta = pi + theta
    ELSE IF ( (x > 0d0) .AND. (y < 0d0) ) THEN
      theta = 2d0*pi - theta
    END IF
  ELSE
    theta = 0d0
  END IF

RETURN
END SUBROUTINE CarToPol
