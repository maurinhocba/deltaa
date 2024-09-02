 FUNCTION functs (iload,ttime)
 !********************************************************************
 !
 !***  heaviside(1), harmonic(2), multi-linear, etc.  time FUNCTIONs
 !
 !OUTPUT: functs(1) = function value,  functs(2) = function derivative
 !********************************************************************
 USE curv_db,ONLY: curpar => curpa
 USE ctrl_db,ONLY: ndyna
 IMPLICIT NONE
 REAL (kind=8) :: functs(2)
 INTEGER (kind=4),INTENT(IN) :: iload
 REAL (kind=8),INTENT(IN) :: ttime

 INTEGER (kind=4) ifunc,n1,n2
 REAL(kind=8):: argum,argu1,a,b,phi,w,t0,te,tc,time
 REAL (kind=8),PARAMETER :: pi = 3.1415926535898_8

 functs = 0d0    !initializes
 time = ttime    !initializes
 IF( iload == 0 ) RETURN     !no associated curve

 n1    = INT(curpar(iload))  !pointer to first parameter
 ifunc = INT(curpar(n1))     !curve type (ltype)
 IF( ifunc == 0) RETURN      !accelerations fixed
 t0 = curpar(n1+2)           !START
 IF( time < t0 )THEN         !function not started yet
   IF( ndyna > 0 .OR. ifunc /= 12 )RETURN  !for dynamical problems
 END IF
 te = curpar(n1+3)           !END time
 IF( time > te) RETURN       !function ended

 n1 = n1+4                   !pointer to parameters

 SELECT CASE (ifunc)         !according to function type
 CASE (1)                                     !constant
   functs = (/ curpar(n1) , 0d0 /)

 CASE (2)                                     !sino
   a  = curpar(n1  )
   b  = curpar(n1+1)
   w  = curpar(n1+2)
   functs(1) = a + b*  SIN(w*(time-t0))
   functs(2) =     b*w*COS(w*(time-t0))

 CASE (3)                                     !multi-linear
   n2 = MAX(n1+1,INT(curpar(n1)))                 !position in list
    DO         !loop to find interval
      IF(time >= curpar(n2) )THEN
        IF(time < curpar(n2+2))THEN
          argum = curpar(n2+1)
          argu1 = curpar(n2+3) - argum
          te = curpar(n2+2) - curpar(n2)
          functs(2) = argu1/te
          functs(1) = argum + functs(2)*(time-curpar(n2))   !multi-linear
          curpar(n1) = n2      !keep position in curve
          EXIT
        ELSE
          n2 = n2+2
        END IF
      ELSE
        n2 = n2-2
      END IF
    END DO

 CASE (4)                                     !cosine
   b  = curpar(n1)
   w  = curpar(n1+1)
   functs(1) = b*(1d0-COS(w*(time-t0)))/2d0
   functs(2) = b*w*SIN(w*(time-t0))/2d0

 CASE (5)                                     !cosine until maximum
   b  = curpar(n1)
   w  = curpar(n1+1)
   argum = w*(time-t0)
   IF(argum < pi) THEN
     functs(1) = b*(1d0-COS(argum))/2d0
     functs(2) = b*w*SIN(argum)/2d0
   ELSE
     functs(1) = b                            !then constant
     functs(2) = 0d0
   END IF

 !CASE (6)                                      !cosine until maximum
 !  b  = curpar(n1  )
 !  w  = curpar(n1+1)
 !  tc = curpar(n1+2)
 !  argum = w*(time-t0)
 !  argu1 = w*(tc-time)
 !  IF(argum <= pi) THEN
 !    functs(1) = b*(1d0-COS(argum))/2d0
 !    functs(2) = b*w*SIN(argum)/2d0
 !  ELSE IF(argu1 <= pi) THEN
 !    functs(1) = b*(1d0-COS(2d0*pi-argu1))/2d0 ! ???
 !    functs(2) = b*w*SIN(2d0*pi-argu1)/2d0
 !  ELSE
 !    functs(1) = b                             !then constant
 !    functs(2) = 0d0                           !then constant
 !  END IF

 CASE (7)                                      !slope-step
   a = curpar(n1  )
   w = curpar(n1+1)
   functs = (/ a , 0d0 /)
   argum  = time - t0
   IF(argum < w) functs = (/ a*argum/w, a/w /)

 CASE (8)                                      !cosine
   a = curpar(n1)
   b = curpar(n1+1)
   w = curpar(n1+2)
   functs(1) = a + b*(1d0-COS(w*(time-t0)))
   functs(2) = b*w*SIN(w*(time-t0))

 CASE (9) !  Cosine increment over a short period over a base
   a  = curpar(n1  )
   b  = curpar(n1+1)
   w  = curpar(n1+2)
   tc = curpar(n1+3)
   functs = (/ a , 0d0 /)
   argum = time - tc
   IF( argum > 0d0 .AND. argum < w )THEN
     functs(1) = a*(1d0+b*(1d0-COS( 6.283185307d0*argum/w )))
     functs(2) = a*b*6.283185307*argum/w*SIN(6.283185307*argum/w)
   END IF

 CASE (10)        !   Linear Decay after an initial delay
   a  = curpar(n1  )
   w  = curpar(n1+1)
   argum = time - t0
   IF( argum >= 0d0 .AND. argum <= w ) functs = (/ a*(1d0-argum/w), -a/w /)

 CASE (11)                                !hat
   a  = curpar(n1  )
   b  = curpar(n1+1)
   w  = curpar(n1+2)
   argum = MODULO((time-t0)/w,1d0)
   IF(argum <= 0.5) functs = (/ a + 2d0*b*argum, 2d0*b/w /)
   IF(argum >  0.5) functs = (/ a + 2d0*b*(1d0-argum), -2d0*b/w /)

 CASE (12)               !linear function
   a = curpar(n1)
   argum = time - t0
   functs = (/ a*argum, a /)

 CASE (13)
   !  add here any special FUNCTION depending on "n" parameters

 END SELECT
 RETURN

 END FUNCTION functs
