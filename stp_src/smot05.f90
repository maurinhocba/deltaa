 SUBROUTINE smot05(ngaus,nstre,dvolu,strsg,strea,shape,facta,nnode)
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

 REAL (kind=8)  auxva,emass(nnode*nnode),streb(nstre,nnode)

 !     evaluates the volume

 facta = 1d0/SUM(dvolu)

 CALL smom05(nnode,ngaus,emass,dvolu,shape)
 !     integrate nodal weighted stresses and compute nodal contributions

 streb = 0d0

 DO n=1,nnode
   DO g=1,ngaus
     auxva = shape(n,g)*dvolu(g)
     DO i=1,nstre
       streb(i,n) = streb(i,n) + auxva*strsg(i,g)
     END DO
   END DO
 END DO

 CALL proma2(strea,streb,emass,nstre,nnode,nnode)

 RETURN
 END SUBROUTINE smot05
