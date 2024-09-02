 SUBROUTINE shapri(shape,ngaus,nnode)
 !
 !  compute Shape functions
 !
 IMPLICIT NONE
 !     arguments
 INTEGER(kind=4), INTENT(IN) :: ngaus,nnode
 REAL (kind=8), INTENT(OUT) :: shape(nnode,ngaus)
 !     local variables
 !  data value for 15-node prism with 6 integration poins
 REAL (kind=8), PARAMETER :: a = 0.166666666666667, b = 0.333333333333333, c = 0.666666666666667, &
                             e = 0.500000000000000, f = 1.666666666666667,                        &
                             l1 =-1.267949192431123,l2 = 1.267949192431123
 REAL (kind=8), PARAMETER :: pt(3,15) = RESHAPE((/                    &
                -b,-b,l1,    f,-b,l1,   -b, f, l1,                    & !lower vertex nodes
                -b,-b,l2,    f,-b,l2,   -b, f, l2,                    & !upper vertex nodes
                 e,-b,l1,    c, c,l1,   -b, c, l1,                    & !lower mid-side nodes
                -b,-b,0.,    f,-b,0.,   -b, f, 0.,                    & !arista mid-side nodes
                 e,-b,l2,    c, c,l2,   -b, c, l2 /),(/3,15/) )         !upper mid-side nodes

 !  data value for 20-node hexahedra with 8 integration poins
 REAL (kind=8), PARAMETER :: g = 1.7320508075688772935274463415059
 REAL (kind=8), PARAMETER :: pq(3,20) = RESHAPE((/                    &
                -g,-g,-g,    g,-g,-g,    g, g,-g,   -g, g,-g,         & !lower vertex nodes
                -g,-g, g,    g,-g, g,    g, g, g,   -g, g, g,         & !upper vertex nodes
                 0,-g,-g,    g, 0,-g,    0, g,-g,   -g, 0,-g,         & !lower mid-side nodes
                -g,-g, 0,    g,-g, 0,    g, g, 0,   -g, g, 0,         & !arista mid-side nodes
                 0,-g, g,    g, 0, g,    0, g, g,   -g, 0, g /),(/3,20/) )         !upper mid-side nodes

 INTEGER (kind=4) :: n
 REAL(kind=8) :: s,t,q,zp,sp,tp,sm,tm,zm

   IF( nnode == 15 )THEN
     DO n=1,nnode ! for each node

       s = pt(1,n)
       t = pt(2,n)
       q = pt(3,n)
       sp = (1d0-q)/2d0   !N1
       tp = (1d0+q)/2d0   !N2
       zp = 1d0-s-t

       shape(n,1) =   zp * sp
       shape(n,2) =   s  * sp
       shape(n,3) =   t  * sp
       shape(n,4) =   zp * tp
       shape(n,5) =   s  * tp
       shape(n,6) =   t  * tp

     END DO
   ELSE IF( nnode == 20 )THEN

     DO n=1,nnode
       s = pq(1,n)
       t = pq(2,n)
       q = pq(3,n)
       sm = 0.5d0*(1d0-s)
       tm = 0.5d0*(1d0-t)
       zm = 0.5d0*(1d0-q)
       sp = 0.5d0*(1d0+s)
       tp = 0.5d0*(1d0+t)
       zp = 0.5d0*(1d0+q)

       !   shape functions

       shape(n,1) = sm*tm*zm
       shape(n,2) = sp*tm*zm
       shape(n,3) = sp*tp*zm
       shape(n,4) = sm*tp*zm
       shape(n,5) = sm*tm*zp
       shape(n,6) = sp*tm*zp
       shape(n,7) = sp*tp*zp
       shape(n,8) = sm*tp*zp
     END DO
   END IF
 RETURN
 END SUBROUTINE shapri
