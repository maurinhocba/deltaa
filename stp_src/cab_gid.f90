 SUBROUTINE cab_gid(vtype,where,vname,vcomp,nc,loadstep,ttime,gpname,rrt,units)
 !
 !   prints results heading for GiD in ASCII form
 !

 ! HAY QUE DECIDIR EL ORDEN DE LAS VARIABLES OPCIONALES, VER PORQUE HE DEJADO '  ' PARA GPNAME EN LAS LLAMADAS !!
 IMPLICIT NONE
 INTEGER(kind=4),INTENT(IN):: vtype  !1=Scalar,2=Vector,3=matrix,4=2d-matrix
 INTEGER(kind=4),INTENT(IN):: where  !1=OnNodes, 2=OnGaussPoints
 INTEGER(kind=4),INTENT(IN):: nc     !number of components
 CHARACTER(len=*),INTENT(IN):: vname !varname
 CHARACTER(len=*),INTENT(IN):: vcomp(nc) !varname
 CHARACTER(len=*),INTENT(IN):: loadstep  !load step type
 REAL(kind=8),INTENT(IN):: ttime          !step time
 CHARACTER(len=*),INTENT(IN),OPTIONAL:: units     !units
 CHARACTER(len=*),INTENT(IN),OPTIONAL:: gpname    !Gauss point name
 CHARACTER(len=*),INTENT(IN),OPTIONAL:: rrt       !result range table

 CHARACTER(len=22),PARAMETER:: v_type(4) = (/'Scalar                ','Vector                ', &
                                             'Matrix                ','PlainDeformationMatrix' /)
 CHARACTER(len=7 ),PARAMETER:: where1 = 'OnNodes'
 CHARACTER(len=15),PARAMETER:: where2 = 'OnGaussPoints "'
 CHARACTER(len=80):: line
 CHARACTER(len=25):: var_u,var_c(6),vu

 INTEGER :: lv,lc(nc),i,ib,ie

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
 lv = LEN_TRIM(var_u)


   IF( where == 2 )THEN  !for Gauss points
     WRITE(13,"('Result ""',a,'"" ""',a10,'"" ',e12.4,' ',a,' ',a16,a,'""')") &
        &  var_u(1:lv),loadstep,ttime,TRIM(v_type(vtype)),where2,TRIM(gpname)
     IF( PRESENT(rrt) ) WRITE(13,"('ResultRangesTable ""',A,'""')")TRIM(rrt)
   ELSE                  !for nodes
     WRITE(13,"('Result ""',a,'"" ""',a10,'"" ',e12.4,' ',a,' ',a8)") &
        &  var_u(1:lv),loadstep,ttime,TRIM(v_type(vtype)),where1
   END IF

   IF(vtype /= 1)THEN    !for vector and matrix variables
     ib = 1
     DO i=1,nc
       lc(i) = LEN_TRIM(var_c(i))
       ie = ib+lc(i)+1
       line(ib:ie) = '"'//var_c(i)(1:lc(i))//'"'
       IF( i < nc) THEN
         ie = ie+1
         line(ie:ie) = ','
         ib = ie+1
       END IF
       IF( vtype == 4 .AND. nc == 3 .AND. i == 3 )THEN
         ie = ie+1
         line(ie:ie) = ','
         ib = ie+1
         ie = ib+7
         line(ib:ie) = '"V_zz=0"'
       END IF
     END DO
     WRITE(13,"('ComponentNames ',a)")line(1:ie)

   ELSE IF( (vname /= vcomp(1)) .AND. (vcomp(1) /= ' ')) THEN  !for scalar variables
     ie = LEN_TRIM(var_c(1)) + 2
     line(1:ie) = '"'//var_c(1)(1:ie-2)//'"'
     WRITE(13,"('ComponentNames ',a)")line(1:ie)
   END IF
   WRITE(13,"('Values')")
 RETURN
 END SUBROUTINE cab_gid
 !
 FUNCTION vu(var_name,units)
 ! joints var_name with units
 IMPLICIT NONE
 CHARACTER(len=25) vu
 CHARACTER(len=*), INTENT(IN) :: var_name,units
 INTEGER(kind=4) :: lv,lu
 lv = LEN_TRIM(var_name)
 lu = LEN_TRIM(units)
 IF( lu == 0 )THEN !no units
   vu = TRIM(var_name)
 ELSE
   vu = TRIM(var_name)//'_('//TRIM(units)//')'
 END IF
 RETURN
 END FUNCTION vu
