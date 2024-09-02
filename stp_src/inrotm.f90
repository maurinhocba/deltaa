 SUBROUTINE inrotm(a,b,g,euler)
 !***************************************************************
 !     calculates the initial nodal coordinate systems
 !***************************************************************
 IMPLICIT NONE

 REAL (kind=8),INTENT(IN) :: a,b,g
 REAL (kind=8),INTENT(OUT) :: euler(9)

 REAL (kind=8) :: ca,cb,cg,sa,sb,sg

 IF(a == 0 .AND. b == 0 .AND. g == 0) THEN

   euler = (/ 1d0 , 0d0 , 0d0 ,  &
              0d0 , 1d0 , 0d0 ,  &
              0d0 , 0d0 , 1d0 /)

 ELSE

   ca = COS(a)
   sa = SIN(a)
   cb = COS(b)
   sb = SIN(b)
   cg = COS(g)
   sg = SIN(g)

   euler = (/   ca*cg - cb*sa*sg  , ca*cb*sg + cg*sa , sb*sg , &
             - (ca*sg + cb*cg*sa) , ca*cb*cg - sa*sg , cg*sb , &
                    sa*sb         ,      - ca*sb     ,   cb   /)
 END IF

 RETURN
 END SUBROUTINE inrotm

 SUBROUTINE rotate( sig, r)

 IMPLICIT NONE

 REAL (kind=8),INTENT(IN) :: r(3,3)
 REAL (kind=8),INTENT(IN OUT) :: sig(6)

 !local variable
 REAL (kind=8) :: str(6)

   str = sig

   sig(1) = r(1,1)*str(1)*r(1,1)+r(1,2)*str(2)*r(1,2)+r(1,3)*str(3)*r(1,3)+ &
       2d0*(r(1,1)*str(4)*r(1,2)+r(1,1)*str(5)*r(1,3)+r(1,2)*str(6)*r(1,3))
   sig(2) = r(2,1)*str(1)*r(2,1)+r(2,2)*str(2)*r(2,2)+r(2,3)*str(3)*r(2,3)+ &
       2d0*(r(2,1)*str(4)*r(2,2)+r(2,1)*str(5)*r(2,3)+r(2,2)*str(6)*r(2,3))
   sig(3) = r(3,1)*str(1)*r(3,1)+r(3,2)*str(2)*r(3,2)+r(3,3)*str(3)*r(3,3)+ &
       2d0*(r(3,1)*str(4)*r(3,2)+r(3,1)*str(5)*r(3,3)+r(3,2)*str(6)*r(3,3))
   sig(4) = r(1,1)*str(1)*r(2,1)+r(1,2)*str(2)*r(2,2)+r(1,3)*str(3)*r(2,3)+ &
            r(1,1)*str(4)*r(2,2)+r(1,1)*str(5)*r(2,3)+r(1,2)*str(6)*r(2,3)+ &
            r(1,2)*str(4)*r(2,1)+r(1,3)*str(5)*r(2,1)+r(1,3)*str(6)*r(2,2)
   sig(5) = r(1,1)*str(1)*r(3,1)+r(1,2)*str(2)*r(3,2)+r(1,3)*str(3)*r(3,3)+ &
            r(1,1)*str(4)*r(3,2)+r(1,1)*str(5)*r(3,3)+r(1,2)*str(6)*r(3,3)+ &
            r(1,2)*str(4)*r(3,1)+r(1,3)*str(5)*r(3,1)+r(1,3)*str(6)*r(3,2)
   sig(6) = r(2,1)*str(1)*r(3,1)+r(2,2)*str(2)*r(3,2)+r(2,3)*str(3)*r(3,3)+ &
            r(2,1)*str(4)*r(3,2)+r(2,1)*str(5)*r(3,3)+r(2,2)*str(6)*r(3,3)+ &
            r(2,2)*str(4)*r(3,1)+r(2,3)*str(5)*r(3,1)+r(2,3)*str(6)*r(3,2)

!   sig(1) = r(1,1)*str(1)*r(1,1)+r(2,1)*str(2)*r(2,1)+r(3,1)*str(3)*r(3,1)+ &
!       2d0*(r(1,1)*str(4)*r(2,1)+r(1,1)*str(5)*r(3,1)+r(2,1)*str(6)*r(3,1))
!   sig(2) = r(1,2)*str(1)*r(1,2)+r(2,2)*str(2)*r(2,2)+r(3,2)*str(3)*r(3,2)+ &
!       2d0*(r(1,2)*str(4)*r(2,2)+r(1,2)*str(5)*r(3,2)+r(2,2)*str(6)*r(3,2))
!   sig(3) = r(1,3)*str(1)*r(1,3)+r(2,3)*str(2)*r(2,3)+r(3,3)*str(3)*r(3,3)+ &
!       2d0*(r(1,3)*str(4)*r(2,3)+r(1,3)*str(5)*r(3,3)+r(2,3)*str(6)*r(3,3))
!   sig(4) = r(1,1)*str(1)*r(1,2)+r(2,1)*str(2)*r(2,2)+r(3,1)*str(3)*r(3,2)+ &
!            r(1,1)*str(4)*r(2,2)+r(1,1)*str(5)*r(3,2)+r(2,1)*str(6)*r(3,2)+ &
!            r(2,1)*str(4)*r(1,2)+r(3,1)*str(5)*r(1,2)+r(3,1)*str(6)*r(2,2)
!   sig(5) = r(1,1)*str(1)*r(1,3)+r(2,1)*str(2)*r(2,3)+r(3,1)*str(3)*r(3,3)+ &
!            r(1,1)*str(4)*r(2,3)+r(1,1)*str(5)*r(3,3)+r(2,1)*str(6)*r(3,3)+ &
!            r(2,1)*str(4)*r(1,3)+r(3,1)*str(5)*r(1,3)+r(3,1)*str(6)*r(2,3)
!   sig(6) = r(1,2)*str(1)*r(1,3)+r(2,2)*str(2)*r(2,3)+r(3,2)*str(3)*r(3,3)+ &
!            r(1,2)*str(4)*r(2,3)+r(1,2)*str(5)*r(3,3)+r(2,2)*str(6)*r(3,3)+ &
!            r(2,2)*str(4)*r(1,3)+r(3,2)*str(5)*r(1,3)+r(3,2)*str(6)*r(2,3)
 RETURN
 END SUBROUTINE rotate
