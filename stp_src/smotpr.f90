 SUBROUTINE smotpr(ngaus,nstre,dvolu,strsg,strea,shape,facta,nnode)
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

 INTEGER (kind=4) i,g,n


 !     evaluates the area

 facta = 1d0/SUM(dvolu)

 !     extrapolate gauss-point values to nodes

 strea = 0d0

 DO n=1,nnode
   DO g=1,ngaus
     DO i=1,nstre
       strea(i,n) = strea(i,n) + shape(n,g) *strsg(i,g)
     END DO
   END DO
 END DO

 RETURN
 END SUBROUTINE smotpr
