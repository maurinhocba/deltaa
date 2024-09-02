SUBROUTINE prgib06( )
!  print results for GiD for shl3d elements
USE data_db,ONLY: npoin, loadstep, ttime, label, shl3d, shl3d_head, shl3d_nvarn,         &
    shl3d_nodes, shl3d_vargs, shl3d_force, shl3d_force_l, shl3d_momen, shl3d_momen_l,    &
    shl3d_shear, shl3d_shear_l, shl3d_logst, shl3d_logst_l, shl3d_eqpst, shl3d_eqpst_l,  &
    shl3d_vmise, shl3d_vmise_l, shl3d_thrat, shl3d_thrat_l,                              &
    shl3d_force_u, shl3d_force_f, shl3d_momen_u, shl3d_momen_f, shl3d_shear_u,           &
    shl3d_shear_f, shl3d_logst_u, shl3d_logst_f, shl3d_eqpst_u, shl3d_eqpst_f,           &
    shl3d_vmise_u, shl3d_vmise_f, shl3d_thrat_u, shl3d_thrat_f
IMPLICIT NONE

  TYPE(shl3d),POINTER:: eset
  INTEGER(kind=4):: ielem,g,iv,is,i,ng,ngaus,nelem
  INTEGER(kind=4),POINTER:: nodes(:,:)
  REAL(kind=8),POINTER:: vargs(:,:)
  CHARACTER(len=34):: gpname
  REAL(kind=8) :: aux(3)
  !         non-standard Gauss positions for GiD
  INTEGER, PARAMETER :: o(4,3) = RESHAPE( (/ 1, 0, 0, 0, 1, 2, 3, 0,  1, 3, 4, 2 /), &
                                          (/ 4, 3 /) )


  INTERFACE
    INCLUDE 'cab_gid_bin.h'
  END INTERFACE

  eset => shl3d_head          !point to first set

  !        print Gaussian variables
  DO                          !loop for all the sets
    gpname = 'GP'//TRIM(eset%sname)   !gauss point name string "GPsname"
    iv = 0  !initializes pointer to Gauss values
    ngaus = eset%ngaus    !number of Gauss points per element
    nelem = eset%nelem    !number of elements
    ng = MAX(ngaus-1,1)   !1 - 2 for triangles 3 for quads

    IF (shl3d_force > 1) THEN
      CALL cab_gid_bin(4,2,shl3d_force_l(5),shl3d_force_l(6),3,loadstep,ttime,gpname,units=shl3d_force_u)
      iv = iv+3
      DO ielem=1,nelem
        DO g=1,ngaus
          aux = eset%elvar(iv-2:iv,o(g,ng),ielem)*shl3d_force_f
          CALL GID_WRITE3DMATRIX(ielem,aux(1),aux(2),aux(3),0d0,0d0,0d0)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (shl3d_momen > 1) THEN
      CALL cab_gid_bin(4,2,shl3d_momen_l(5),shl3d_momen_l(6),3,loadstep,ttime,gpname,units=shl3d_momen_u)
      iv = iv+3
      DO ielem=1,nelem
        DO g=1,ngaus
          aux = eset%elvar(iv-2:iv,o(g,ng),ielem)*shl3d_momen_f
          CALL GID_WRITE3DMATRIX(ielem,aux(1),aux(2),aux(3),0d0,0d0,0d0)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (shl3d_shear > 1) THEN
      CALL cab_gid_bin(1,2,shl3d_shear_l(4),shl3d_shear_l(5),1,loadstep,ttime,gpname,units=shl3d_shear_u)
      iv = iv+1
      DO ielem=1,nelem
        DO g=1,ngaus
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*shl3d_shear_f)
        END DO
      END DO
      CALL GID_ENDRESULT()

      CALL cab_gid_bin(1,2,shl3d_shear_l(4),shl3d_shear_l(6),1,loadstep,ttime,gpname,units=shl3d_shear_u)
      iv = iv+1
      DO ielem=1,nelem
        DO g=1,ngaus
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*shl3d_shear_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (shl3d_logst > 1) THEN
      CALL cab_gid_bin(2,2,shl3d_logst_l(5),shl3d_logst_l(6),3,loadstep,ttime,gpname,units=shl3d_logst_u)
      iv = iv+3
      DO ielem=1,nelem
        DO g=1,ngaus
          aux = eset%elvar(iv-2:iv,o(g,ng),ielem)*shl3d_logst_f
          CALL GID_WRITEVECTOR(ielem,aux(1),aux(2),aux(3))
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (shl3d_eqpst > 1) THEN
      CALL cab_gid_bin(1,2,shl3d_eqpst_l(3),shl3d_eqpst_l(4),1,loadstep,ttime,gpname,units=shl3d_eqpst_u)
      iv = iv+1
      DO ielem=1,nelem
        DO g=1,ngaus
          CALL GID_WRITESCALAR(ielem,MAX(eset%elvar(iv,o(g,ng),ielem),0d0)*shl3d_eqpst_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (shl3d_vmise > 1) THEN
      CALL cab_gid_bin(1,2,shl3d_vmise_l(3),shl3d_vmise_l(4),1,loadstep,ttime,gpname,units=shl3d_vmise_u)
      iv = iv+1
      DO ielem=1,nelem
        DO g=1,ngaus
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*shl3d_vmise_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (shl3d_thrat > 1) THEN
      CALL cab_gid_bin(1,2,shl3d_thrat_l(3),shl3d_thrat_l(4),1,loadstep,ttime,gpname,units=shl3d_thrat_u)
      iv = iv+1
      DO ielem=1,nelem
        DO g=1,ngaus
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*shl3d_thrat_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    eset => eset%next
    IF (.NOT.ASSOCIATED(eset)) EXIT
  END DO

  !        print Nodal variables

  IF (shl3d_nvarn == 0) RETURN    !no variables to print, => exit

  nodes => shl3d_nodes
  vargs => shl3d_vargs
  iv = 0

  IF (MOD(shl3d_force,2) == 1) THEN
    CALL cab_gid_bin(4,1,shl3d_force_l(1),shl3d_force_l(2),3,loadstep,ttime,gpname,units=shl3d_force_u)
    iv = iv+3
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      aux = vargs(iv-2:iv,is)*shl3d_force_f
      CALL GID_WRITE3DMATRIX(label(i),aux(1),aux(2),aux(3),0d0,0d0,0d0)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(shl3d_momen,2) == 1) THEN
    CALL cab_gid_bin(4,1,shl3d_momen_l(1),shl3d_momen_l(2),3,loadstep,ttime,gpname,units=shl3d_momen_u)
    iv = iv+3
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      aux = vargs(iv-2:iv,is)*shl3d_momen_f
      CALL GID_WRITE3DMATRIX(label(i),aux(1),aux(2),aux(3),0d0,0d0,0d0)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(shl3d_shear,2) == 1) THEN
    CALL cab_gid_bin(1,1,shl3d_shear_l(2),shl3d_shear_l(2),1,loadstep,ttime,gpname,units=shl3d_shear_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*shl3d_shear_f)
    END DO
    CALL GID_ENDRESULT()

    CALL cab_gid_bin(1,1,shl3d_shear_l(3),shl3d_shear_l(3),1,loadstep,ttime,gpname,units=shl3d_shear_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*shl3d_shear_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(shl3d_logst,2) == 1) THEN
    CALL cab_gid_bin(2,1,shl3d_logst_l(1),shl3d_logst_l(2),3,loadstep,ttime,gpname,units=shl3d_logst_u)
    iv = iv+3
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      aux = vargs(iv-2:iv,is)*shl3d_logst_f
      CALL GID_WRITEVECTOR(label(i),aux(1),aux(2),aux(3))
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(shl3d_eqpst,2) ==  1) THEN
    CALL cab_gid_bin(1,1,shl3d_eqpst_l(1),shl3d_eqpst_l(2),1,loadstep,ttime,gpname,units=shl3d_eqpst_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*shl3d_eqpst_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(shl3d_vmise,2) ==  1) THEN
    CALL cab_gid_bin(1,1,shl3d_vmise_l(1),shl3d_vmise_l(2),1,loadstep,ttime,gpname,units=shl3d_vmise_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*shl3d_vmise_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(shl3d_thrat,2) ==  1) THEN
    CALL cab_gid_bin(1,1,shl3d_thrat_l(1),shl3d_thrat_l(2),1,loadstep,ttime,gpname,units=shl3d_thrat_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*shl3d_thrat_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

RETURN
END SUBROUTINE prgib06
