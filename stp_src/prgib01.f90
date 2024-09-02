SUBROUTINE prgib01( )
!  print results for GiD for spot elements
USE data_db,ONLY: npoin, loadstep, ttime, label, spot, spot_head,   &
    spot_nodes, spot_force, spot_force_l, spot_shear, spot_shear_l,   &
    spot_eqpst, spot_eqpst_l, &
    spot_force_u, spot_shear_u, spot_eqpst_u, spot_force_f, spot_shear_f, spot_eqpst_f
IMPLICIT NONE

  TYPE(spot),POINTER:: eset
  INTEGER(kind=4):: ielem,g,iv,is,i
  INTEGER(kind=4),POINTER:: nodes(:,:)
  REAL(kind=8),POINTER:: vargs(:,:)
  CHARACTER(len=34):: gpname

  INTERFACE
    INCLUDE 'cab_gid_bin.h'
  END INTERFACE



  eset => spot_head          !point to first set

  !        print Gaussian variables
  DO                          !loop for all the sets
    gpname = 'GP'//TRIM(eset%sname)   !gauss point name string "GPsname"
    iv = 0  !initializes pointer to Gauss values

    IF (spot_force > 1) THEN
      CALL cab_gid_bin(1,2,spot_force_l(1),spot_force_l(2),1,loadstep,ttime,gpname,units=spot_force_u)
      iv = iv+1
      DO ielem=1,eset%nelem
        DO g=1,eset%ngaus
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,g,ielem)*spot_force_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (spot_shear > 1) THEN
      CALL cab_gid_bin(1,2,spot_shear_l(1),spot_shear_l(2),1,loadstep,ttime,gpname,units=spot_shear_u)
      iv = iv+1
      DO ielem=1,eset%nelem
        DO g=1,eset%ngaus
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,g,ielem)*spot_shear_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (spot_eqpst > 1) THEN
      CALL cab_gid_bin(1,2,spot_eqpst_l(2),spot_eqpst_l(2),1,loadstep,ttime,gpname,units=spot_eqpst_u)
      iv = iv+1
      DO ielem=1,eset%nelem
        DO g=1,eset%ngaus
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,g,ielem)*spot_eqpst_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    eset => eset%next
    IF (.NOT.ASSOCIATED(eset)) EXIT
  END DO

RETURN
END SUBROUTINE prgib01
