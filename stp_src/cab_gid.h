 SUBROUTINE cab_gid(vtype,where,vname,vcomp,nc,loadstep,ttime,gpname,rrt,units)
 !
 !
 !
 IMPLICIT NONE
 INTEGER(kind=4),INTENT(IN):: vtype  !1=Scalar,2=Vector,3=matrix,4=2d-matrix
 INTEGER(kind=4),INTENT(IN):: where  !1=OnNodes, 2=OnGaussPoints, 3=CotourRanges
 INTEGER(kind=4),INTENT(IN):: nc     !number of components
 CHARACTER(len=*),INTENT(IN):: vname !varname
 CHARACTER(len=*),INTENT(IN):: vcomp(nc) !varname
 CHARACTER(len=*),INTENT(IN):: loadstep  !load step type
 REAL(kind=8),INTENT(IN):: ttime          !step time
 CHARACTER(len=*),INTENT(IN),OPTIONAL:: units     !units
 CHARACTER(len=*),INTENT(IN),OPTIONAL:: gpname    !Gauss point name
 CHARACTER(len=*),INTENT(IN),OPTIONAL:: rrt       !result range table

 END SUBROUTINE cab_gid
