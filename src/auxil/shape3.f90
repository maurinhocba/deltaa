      SUBROUTINE shape3(deriv,shape,s,t,nnode )
      !********************************************************************
      !
      !***  calculates shape functions and their derivatives for 2d elements
      !
      !********************************************************************
      IMPLICIT NONE
      INTEGER (kind=4), INTENT(IN) :: nnode
      REAL (kind=8), INTENT(IN) :: s,t
      REAL (kind=8), INTENT(OUT) ::   deriv(nnode,2),shape(nnode)
      ! local auxiliar variables
      REAL (kind=8)    st,ss,tt,s2,t2,sst,stt,st2

      IF (nnode == 4) THEN

         st = s*t
         !  shape functions for 4 noded element

         shape(1) = (1D0 - s - t + st)/4D0
         shape(2) = (1D0 + s - t - st)/4D0
         shape(3) = (1D0 + s + t + st)/4D0
         shape(4) = (1D0 - s + t - st)/4D0

         ! and derivatives

         deriv(1,1) = (-1d0 + t)/4d0
         deriv(2,1) = -deriv(1,1)
         deriv(3,1) = ( 1d0 + t)/4d0
         deriv(4,1) = -deriv(3,1)
         deriv(1,2) = (-1d0 + s)/4d0
         deriv(2,2) = (-1d0 - s)/4d0
         deriv(3,2) = -deriv(2,2)
         deriv(4,2) = -deriv(1,2)

      ELSE IF (nnode == 3) THEN

         !  shape functions for 3-node element

         shape(1) = 1.d0 - s - t
         shape(2) = s
         shape(3) = t

         !  and derivatives

         deriv(1,1) = -1.d0
         deriv(1,2) = -1.d0
         deriv(2,1) =  1.d0
         deriv(2,2) =  0.d0
         deriv(3,1) =  0.d0
         deriv(3,2) =  1.d0

      ELSE IF (nnode == 8) THEN
        !  my version
        !ss = s*s
        !st = s*t
        !tt = t*t
        !!  shape functions for 8 noded serendipity element
        !! corner nodes
        !shape(1) = (1D0 - s - t + st)*(-s-t-1d0)/4D0
        !shape(2) = (1D0 + s - t - st)*(+s-t-1d0)/4D0
        !shape(3) = (1D0 + s + t + st)*(+s+t-1d0)/4D0
        !shape(4) = (1D0 - s + t - st)*(-s+t-1d0)/4D0
        !! mid-side nodes
        !shape(5) = (1D0 - ss)*(1d0-t)/2D0
        !shape(6) = (1D0 - tt)*(1d0+s)/2D0
        !shape(7) = (1D0 - ss)*(1d0+t)/2D0
        !shape(8) = (1D0 - tt)*(1d0-s)/2D0
        !
        !! and derivatives
        !deriv(1,1) = ( 2D0*s +t -2d0*st -tt)/4D0
        !deriv(1,2) = (+s +2d0*t -ss -2d0*st)/4D0
        !deriv(2,1) = ( 2D0*s -t -2d0*st +tt)/4D0
        !deriv(2,2) = (-s +2d0*t -ss +2d0*st)/4D0
        !deriv(3,1) = ( 2D0*s +t +2d0*st +tt)/4D0
        !deriv(3,2) = (+s +2d0*t +ss +2d0*st)/4D0
        !deriv(4,1) = ( 2D0*s -t +2d0*st -tt)/4D0
        !deriv(4,2) = (-s +2d0*t +ss -2d0*st)/4D0
        !deriv(5,1) = -s*(1d0 -t)
        !deriv(5,2) = -(1d0 - ss)/2d0
        !deriv(6,1) = +(1d0 - tt)/2d0
        !deriv(6,2) = -t*(1d0 +s)
        !deriv(7,1) = -s*(1d0 +t)
        !deriv(7,2) = +(1d0 - ss)/2d0
        !deriv(8,1) = -(1d0 - tt)/2d0
        !deriv(8,2) = -t*(1d0 -s)

        !  STP version
         st = s*t
         s2 = s*2d0
         t2 = t*2d0
         ss = s*s
         tt = t*t
         sst = ss*t
         stt = st*t
         st2 = st*2d0

         !     ***  shape functions

         shape(1) = (-1d0+ st+ ss+ tt- sst- stt)/4d0
         shape(2) = (-1d0- st+ ss+ tt- sst+ stt)/4d0
         shape(3) = (-1d0+ st+ ss+ tt+ sst+ stt)/4d0
         shape(4) = (-1d0- st+ ss+ tt+ sst- stt)/4d0
         shape(5) = ( 1d0- t- ss+ sst)/2d0
         shape(6) = ( 1d0+ s- tt- stt)/2d0
         shape(7) = ( 1d0+ t- ss- sst)/2d0
         shape(8) = ( 1d0- s- tt+ stt)/2d0

         !     *** ant terivatives

         deriv(1,1) = ( t+ s2- st2- tt)/4d0
         deriv(2,1) = (-t+ s2- st2+ tt)/4d0
         deriv(3,1) = ( t+ s2+ st2+ tt)/4d0
         deriv(4,1) = (-t+ s2+ st2- tt)/4d0
         deriv(5,1) = -s+ st
         deriv(6,1) = (1d0- tt)/2d0
         deriv(7,1) = -s- st
         deriv(8,1) = (-1d0+ tt)/2d0
         deriv(1,2) = ( s+ t2- ss- st2)/4d0
         deriv(2,2) = (-s+ t2- ss+ st2)/4d0
         deriv(3,2) = ( s+ t2+ ss+ st2)/4d0
         deriv(4,2) = (-s+ t2+ ss- st2)/4d0
         deriv(5,2) = (-1d0+ ss)/2d0
         deriv(6,2) =  -t- st
         deriv(7,2) = (1d0- ss)/2d0
         deriv(8,2) =  -t+ st
      END IF

      END SUBROUTINE shape3
