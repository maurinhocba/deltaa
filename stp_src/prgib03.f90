SUBROUTINE prgib03( )
!  print results for GiD for
!  2-D solid elements
USE data_db,ONLY: npoin, loadstep, ttime, label, sol2d, sol2d_head, sol2d_nvarn,         &
    sol2d_nodes, sol2d_vargs, sol2d_stres, sol2d_stres_l, sol2d_logst, sol2d_logst_l,    &
    sol2d_shtst, sol2d_shtst_l,sol2d_thrat, sol2d_thrat_l, sol2d_eqpst, sol2d_eqpst_l,   &
    sol2d_vmise, sol2d_vmise_l, sol2d_fldma, sol2d_fldma_l, sol2d_wfldFZ, sol2d_fldFZ_l, &
    sol2d_wfldSZ, sol2d_fldSZ_l, rrt_fr, rrt_sf, sol2d_stres_u, sol2d_stres_f,           &
    sol2d_logst_u, sol2d_logst_f, sol2d_shtst_u, sol2d_shtst_f, sol2d_thrat_u,           &
    sol2d_thrat_f, sol2d_eqpst_u, sol2d_eqpst_f, sol2d_vmise_u, sol2d_vmise_f
IMPLICIT NONE

  ! local variables
  TYPE(sol2d),POINTER:: eset
  INTEGER(kind=4):: ielem,g,i,iv,nelem,nnode,ngaus,is,ng
  INTEGER(kind=4),POINTER:: nodes(:,:)
  REAL(kind=8),POINTER:: vargs(:,:)
  CHARACTER(len=34):: gpname
  REAL(kind=8) :: aux(4)
  !     non standart GiD Gauss point positions (INTERNAL)
  INTEGER, PARAMETER :: o(9,3) = RESHAPE( (/ 1, 0, 0, 0, 0, 0, 0, 0, 0,   &
                                             1, 3, 4, 2, 0, 0, 0, 0, 0,   &
                                             1, 7, 9, 3, 4, 8, 6, 2, 5 /),&
                                             (/ 9,3 /) )

  INTERFACE
    INCLUDE 'cab_gid_bin.h'
  END INTERFACE

  eset => sol2d_head

  !        print Gaussian variables
  DO
    gpname = 'GP'//TRIM(eset%sname)
    nelem = eset%nelem
    nnode = eset%nnode
    ngaus = eset%ngaus
    ng = NINT(SQRT(REAL(ngaus)))
    iv = 0       !initializes pointer

    IF (sol2d_stres > 1) THEN
      CALL cab_gid_bin(4,2,sol2d_stres_l(6),sol2d_stres_l(7),4,loadstep,ttime,gpname,units=sol2d_stres_u)
      iv = iv+4
      aux = 0d0
      DO ielem=1,nelem
        DO g=1,ngaus
          aux(1:4) = eset%elvar(iv-3:iv,o(g,ng),ielem)*sol2d_stres_f
          CALL GID_WRITE3DMATRIX(ielem,aux(1),aux(2),aux(3),aux(4),0d0,0d0)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (sol2d_logst > 1) THEN
      CALL cab_gid_bin(2,2,sol2d_logst_l(5),sol2d_logst_l(6),3,loadstep,ttime,gpname,units=sol2d_logst_u)
      iv = iv+3
      aux = 0d0
      DO ielem=1,nelem
        DO g=1,ngaus
          aux(1:3) = eset%elvar(iv-2:iv,o(g,ng),ielem)*sol2d_logst_f
          CALL GID_WRITEVECTOR(ielem,aux(1),aux(2),aux(3))
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (sol2d_shtst > 1) THEN
      CALL cab_gid_bin(2,2,sol2d_shtst_l(5),sol2d_shtst_l(6),3,loadstep,ttime,gpname,units=sol2d_shtst_u)
      iv = iv+3
      aux = 0d0
      DO ielem=1,nelem
        DO g=1,ngaus
          aux(1:3) = eset%elvar(iv-2:iv,o(g,ng),ielem)*sol2d_shtst_f
          CALL GID_WRITEVECTOR(ielem,aux(1),aux(2),aux(3))
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (sol2d_thrat > 1) THEN
      CALL cab_gid_bin(1,2,sol2d_thrat_l(3),sol2d_thrat_l(4),1,loadstep,ttime,gpname,units=sol2d_thrat_u)
      iv = iv+1
      DO ielem=1,nelem
        DO g=1,ngaus
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*sol2d_thrat_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (sol2d_eqpst > 1) THEN
      CALL cab_gid_bin(1,2,sol2d_eqpst_l(3),sol2d_eqpst_l(4),1,loadstep,ttime,gpname,units=sol2d_eqpst_u)
      iv = iv+1
      DO ielem=1,nelem
        DO g=1,ngaus
          CALL GID_WRITESCALAR(ielem,MAX(eset%elvar(iv,o(g,ng),ielem)*sol2d_eqpst_f,0d0))
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (sol2d_vmise > 1) THEN
      CALL cab_gid_bin(1,2,sol2d_vmise_l(3),sol2d_vmise_l(4),1,loadstep,ttime,gpname,units=sol2d_vmise_u)
      iv = iv+1
      DO ielem=1,nelem
        DO g=1,ngaus
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*sol2d_vmise_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (sol2d_fldma > 1) THEN
      CALL cab_gid_bin(1,2,sol2d_fldma_l(3),sol2d_fldma_l(4),1,loadstep,ttime,gpname)
      iv = iv+1
      DO ielem=1,nelem
        DO g=1,ngaus
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem))
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (sol2d_wfldFZ) THEN
      CALL cab_gid_bin(1,2,sol2d_fldFZ_l(1),sol2d_fldFZ_l(1),1,loadstep,ttime,gpname,TRIM(rrt_fr%name))
      iv = iv+1
      DO ielem=1,nelem
        DO g=1,ngaus
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem))
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (sol2d_wfldSZ) THEN
      CALL cab_gid_bin(1,2,sol2d_fldSZ_l(1),sol2d_fldSZ_l(1),1,loadstep,ttime,gpname,TRIM(rrt_sf%name))
      iv = iv+1
      DO ielem=1,nelem
        DO g=1,ngaus
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem))
        END DO
      END DO
      CALL GID_ENDRESULT()

    END IF

    eset => eset%next
    IF (.NOT.ASSOCIATED(eset)) EXIT
  END DO

  !        print Nodal variables

  IF (sol2d_nvarn == 0) RETURN    !no variables to print, => exit

  nodes => sol2d_nodes
  vargs => sol2d_vargs
  iv = 0       !initializes pointer to vargs

  IF (MOD(sol2d_stres,2) == 1) THEN
    CALL cab_gid_bin(4,1,sol2d_stres_l(1),sol2d_stres_l(2),4,loadstep,ttime,gpname,units=sol2d_stres_u)
    iv = iv+4
    is = 0
    aux = 0d0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      aux(1:4) = vargs(iv-3:iv,is)*sol2d_stres_f
      CALL GID_WRITE3DMATRIX(label(i),aux(1),aux(2),aux(3),aux(4),0d0,0d0)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(sol2d_logst,2) == 1) THEN
    CALL cab_gid_bin(2,1,sol2d_logst_l(1),sol2d_logst_l(2),3,loadstep,ttime,gpname,units=sol2d_logst_u)
    iv = iv+3
    is = 0
    aux = 0d0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      aux(1:3) = vargs(iv-2:iv,is)*sol2d_logst_f
      CALL GID_WRITEVECTOR(label(i),aux(1),aux(2),aux(3))
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(sol2d_shtst,2) == 1) THEN
    CALL cab_gid_bin(2,1,sol2d_shtst_l(1),sol2d_shtst_l(2),3,loadstep,ttime,gpname,units=sol2d_shtst_u)
    iv = iv+3
    is = 0
    aux = 0d0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      aux(1:3) = vargs(iv-2:iv,is)*sol2d_shtst_f
      CALL GID_WRITEVECTOR(label(i),aux(1),aux(2),aux(3))
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(sol2d_thrat,2) == 1) THEN
    CALL cab_gid_bin(1,1,sol2d_thrat_l(1),sol2d_thrat_l(2),1,loadstep,ttime,gpname,units=sol2d_thrat_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*sol2d_thrat_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(sol2d_eqpst,2) == 1) THEN
    CALL cab_gid_bin(1,1,sol2d_eqpst_l(1),sol2d_eqpst_l(2),1,loadstep,ttime,gpname,units=sol2d_eqpst_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*sol2d_eqpst_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(sol2d_vmise,2) == 1) THEN
    CALL cab_gid_bin(1,1,sol2d_vmise_l(1),sol2d_vmise_l(2),1,loadstep,ttime,gpname,units=sol2d_vmise_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*sol2d_vmise_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(sol2d_fldma,2) == 1) THEN
    CALL cab_gid_bin(1,1,sol2d_fldma_l(1),sol2d_fldma_l(2),1,loadstep,ttime,gpname)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is))
    END DO
    CALL GID_ENDRESULT()
  END IF

RETURN
END SUBROUTINE prgib03
