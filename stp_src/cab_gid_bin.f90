SUBROUTINE cab_gid_bin(vtype,where,vname,vcomp,nc,loadstep,ttime,gpname,rrt,units)
IMPLICIT NONE

  INTEGER(kind=4),INTENT(IN):: vtype  !1=Scalar,2=Vector,3=matrix,4=2d-matrix
  INTEGER(kind=4),INTENT(IN):: where  !1=OnNodes, 2=OnGaussPoints, 3=CotourRanges
  INTEGER(kind=4),INTENT(IN):: nc     !number of components
  CHARACTER(len=*),INTENT(IN):: vname !varname
  CHARACTER(len=*),INTENT(IN):: vcomp(nc) !components of varname
  CHARACTER(len=*),INTENT(IN):: loadstep  !load step type
  REAL(kind=8),INTENT(IN):: ttime          !step time
  CHARACTER(len=*),INTENT(IN),OPTIONAL:: gpname    !Gauss point name
  CHARACTER(len=*),INTENT(IN),OPTIONAL:: rrt       !Result Range Table
  CHARACTER(len=*),INTENT(IN),OPTIONAL:: units     !variable units

  CHARACTER(len=4):: NULL
  CHARACTER(len=25):: lrrt,var_u,var_c(6),vu
  INTEGER :: l,i

  NULL = CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
  IF( PRESENT(rrt)) THEN
    l = LEN_TRIM(rrt)
    lrrt(1:l) = TRIM(rrt)
  ELSE
    l = 4
    lrrt(1:4) = NULL
  END IF

  IF( PRESENT(units)) THEN
    var_u = vu(vname,units)
    DO i=1,nc
      var_c(i) = vu(vcomp(i),units)
    END DO
  ELSE
    var_u = vname
    DO i=1,nc
      var_c(i) = vcomp(i)
    END DO
  END IF

  IF( where == 2 )THEN  !for Gauss points
    SELECT CASE(vtype)
    CASE (1)   !Scalar
      IF( PRESENT(rrt) )THEN
        CALL GID_BEGINSCALARRESULT(TRIM(var_u),TRIM(loadstep),ttime,1,TRIM(gpname),      &
             TRIM(lrrt),NULL )
      ELSE IF( TRIM(vname) /= TRIM(vcomp(1)) .AND. LEN_TRIM(vcomp(1)) > 0 )THEN
        CALL GID_BEGINSCALARRESULT(TRIM(var_u),TRIM(loadstep),ttime,1,TRIM(gpname),      &
             NULL,TRIM(var_c(1)))
      ELSE
        CALL GID_BEGINSCALARRESULT(TRIM(var_u),TRIM(loadstep),ttime,1,TRIM(gpname),      &
             NULL,NULL)
      END IF

    CASE (2)   !Vector
      IF (nc == 3) THEN
        CALL GID_BEGINVECTORRESULT(TRIM(var_u),TRIM(loadstep),ttime,1,TRIM(gpname),NULL, &
          TRIM(var_c(1)),TRIM(var_c(2)),TRIM(var_c(3)),NULL)
      ELSE ! IF (nc == 2) THEN
        CALL GID_BEGINVECTORRESULT(TRIM(var_u),TRIM(loadstep),ttime,1,TRIM(gpname),NULL, &
          TRIM(var_c(1)),TRIM(var_c(2)),NULL,NULL)
      END IF

    CASE (3)   !Matrix
      IF (nc == 6) THEN
        CALL GID_BEGIN3DMATRESULT(TRIM(var_u),TRIM(loadstep),ttime,1,TRIM(gpname),NULL,  &
          TRIM(var_c(1)),TRIM(var_c(2)),TRIM(var_c(3)),TRIM(var_c(4)),TRIM(var_c(5)),    &
          TRIM(var_c(6)))
      ELSE ! IF (nc == 3) THEN
        CALL GID_BEGIN3DMATRESULT(TRIM(var_u),TRIM(loadstep),ttime,1,TRIM(gpname),NULL,  &
          TRIM(var_c(1)),TRIM(var_c(2)),TRIM(var_c(3)),NULL,NULL,NULL)
      END IF

    CASE (4)   !Matrix
      IF (nc == 4) THEN
        CALL GID_BEGIN3DMATRESULT(TRIM(var_u),TRIM(loadstep),ttime,1,TRIM(gpname),NULL,  &
          TRIM(var_c(1)),TRIM(var_c(2)),TRIM(var_c(3)),TRIM(var_c(4)),NULL,NULL)
        !CALL GiD_BeginPDMMatResult(TRIM(var_u),TRIM(loadstep),ttime,1,TRIM(gpname),NULL,  &
        !  TRIM(var_c(1)),TRIM(var_c(2)),TRIM(var_c(3)),TRIM(var_c(4)),NULL,NULL)

      ELSE ! IF (nc == 3) THEN
        CALL GID_BEGIN3DMATRESULT(TRIM(var_u),TRIM(loadstep),ttime,1,TRIM(gpname),NULL,  &
          TRIM(var_c(1)),TRIM(var_c(2)),'"V_zz=0"',TRIM(var_c(3)),NULL,NULL)
        ! CALL GiD_Begin2DMatResult(TRIM(var_u),TRIM(loadstep),ttime,1,TRIM(gpname),NULL,  &
        !  TRIM(var_c(1)),TRIM(var_c(2)),TRIM(var_c(3)))

      END IF

    END SELECT

  ELSE                  !for nodes
    SELECT CASE(vtype)
    CASE (1)   !Scalar
      IF(vname /= vcomp(1) .AND. vcomp(1) /= ' ') THEN
        CALL GID_BEGINSCALARRESULT(TRIM(var_u),TRIM(loadstep),ttime,0,   &
             NULL,NULL,NULL)
      ELSE
        CALL GID_BEGINSCALARRESULT(TRIM(var_u),TRIM(loadstep),ttime,0,   &
             NULL,NULL,TRIM(var_c(1)))
      END IF

    CASE (2)   !Vector
      IF (nc == 3) THEN
        CALL GID_BEGINVECTORRESULT(TRIM(var_u),TRIM(loadstep),ttime,0,NULL,NULL,         &
          TRIM(var_c(1)),TRIM(var_c(2)),TRIM(var_c(3)),NULL)
      ELSE ! IF (nc == 2) THEN
        CALL GID_BEGINVECTORRESULT(TRIM(var_u),TRIM(loadstep),ttime,0,NULL,NULL,         &
          TRIM(var_c(1)),TRIM(var_c(2)),NULL,NULL)
      END IF

    CASE (3)   !Matrix
      IF (nc == 6) THEN
        CALL GID_BEGIN3DMATRESULT(TRIM(var_u),TRIM(loadstep),ttime,0,NULL,NULL,          &
          TRIM(var_c(1)),TRIM(var_c(2)),TRIM(var_c(3)),TRIM(var_c(4)),TRIM(var_c(5)),    &
          TRIM(var_c(6)))
      ELSE ! IF (nc == 3) THEN
        CALL GID_BEGIN3DMATRESULT(TRIM(var_u),TRIM(loadstep),ttime,0,NULL,NULL,          &
          TRIM(var_c(1)),TRIM(var_c(2)),TRIM(var_c(3)),NULL,NULL,NULL)
      END IF

    CASE (4)   !Matrix
      IF (nc == 4) THEN
        CALL GID_BEGIN3DMATRESULT(TRIM(var_u),TRIM(loadstep),ttime,0,NULL,NULL,          &
          TRIM(var_c(1)),TRIM(var_c(2)),TRIM(var_c(3)),TRIM(var_c(4)),NULL,NULL)
      ELSE ! IF (nc == 3) THEN
        CALL GID_BEGIN3DMATRESULT(TRIM(var_u),TRIM(loadstep),ttime,0,NULL,NULL,          &
          TRIM(var_c(1)),TRIM(var_c(2)),'"V_zz=0"',TRIM(var_c(3)),NULL,NULL)
      END IF

    END SELECT
  END IF

RETURN
END SUBROUTINE cab_gid_bin
