 SUBROUTINE smotpr6(ngaus,nstre,dvolu,strsg,strea,shape,facta,nnode)
 !***********************************************************************
 !
 !****this routine computes the elemental contribution to nodal stresses
 !       -----> strea
 !***********************************************************************
 IMPLICIT NONE

 !     routine parameters

 INTEGER (kind=4), INTENT(IN) :: ngaus, & !number of Gauss points
                                 nstre, & !number of variables to smooth
                                 nnode    !number of nodes per element

 REAL(kind=8), INTENT(IN) :: strsg(nstre,ngaus), & !Gauss point values
                             dvolu(ngaus),       & !Gauss point volumes
                             shape(nnode,ngaus)    !Nodal shape functions
 REAL(kind=8), INTENT(OUT) :: facta,  &          !nodal factor
                              strea(nstre,nnode) !nodal contributions

 !     local variables

 INTEGER (kind=4) i,ng1
 REAL(kind=8) :: f(2,2)


 ng1 = ngaus - 1
 !     compute extrapolation factors
 f(2,1) =    -(0.333333333333d0-shape(1,1))/(shape(1,1)-shape(1,2))
 f(1,1) = 1d0-f(2,1)
 f(1,2) =    -(0.333333333333d0-shape(4,ngaus))/(shape(4,ngaus)-shape(4,ng1))
 f(2,2) = 1d0-f(1,2)

 facta = 1d0/SUM(dvolu)

 !     extrapolate gauss-point values to nodes

 DO i=1,nstre
   strea(i,1:3) = f(1,1)*strsg(i,1)   + f(2,1)*strsg(i,2)
   strea(i,4:6) = f(1,2)*strsg(i,ng1) + f(2,2)*strsg(i,ngaus)
 END DO

 RETURN
 END SUBROUTINE smotpr6
