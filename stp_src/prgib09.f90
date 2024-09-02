SUBROUTINE prgib09(etype )
!  print results for GiD for shrev elements
USE data_db,ONLY: npoin, ntype, loadstep, ttime, label, shrev, shrev_head, shrev_nvarn,  &
    shrev_nodes, shrev_vargs, shrev_force, shrev_force_l, shrev_shear, shrev_shear_l,    &
    shrev_momen, shrev_momen_l, shrev_eqpst, shrev_eqpst_l, shrev_vmise, shrev_vmise_l,  &
    shrev_thrat, shrev_thrat_l, shrev_force_u, shrev_force_f, shrev_shear_u, shrev_shear_f,    &
    shrev_momen_u, shrev_momen_f, shrev_eqpst_u, shrev_eqpst_f, shrev_vmise_u, shrev_vmise_f,  &
    shrev_thrat_u, shrev_thrat_f
IMPLICIT NONE
  INTEGER, INTENT(IN) :: etype   !element type

  TYPE(shrev),POINTER:: eset
  INTEGER(kind=4):: ielem,g,iv,is,i,ng,nst
  INTEGER(kind=4),POINTER:: nodes(:,:)
  REAL(kind=8),POINTER:: vargs(:,:)
  CHARACTER(len=34):: gpname,auxvl
  CHARACTER(len=4):: NULL
  INTEGER(kind=4),PARAMETER:: o(3,3) = RESHAPE( (/ 1, 0, 0,  1, 2, 0, 1, 3, 2 /), &
                                                (/ 3, 3 /) )


  INTERFACE
    INCLUDE 'cab_gid_bin.h'
  END INTERFACE

  NULL = CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)

  nst = 2
  IF( ntype == 1 )nst = 1

  eset => shrev_head          !point to first set

  !        print Gaussian variables
  DO                          !loop for all the sets
    gpname = 'GP'//TRIM(eset%sname)   !gauss point name string "GPsname"
    iv = 0  !initializes pointer to Gauss values
    ng = eset%ngaus

    ! There are problems with values written at Gauss points
    IF(shrev_force > 1)THEN
     !CALL cab_gid_bin(1,2,shrev_force_l(4),shrev_force_l(5),2,loadstep,ttime,gpname,units=shrev_force_u)
      CALL GID_BEGINSCALARRESULT(TRIM(shrev_force_l(5)),TRIM(loadstep),ttime,1,TRIM(gpname),NULL,NULL)
      iv = iv+1
      DO ielem=1,eset%nelem
        DO g=1,ng
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*shrev_force_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
      IF(ntype > 1)THEN
        CALL GID_BEGINSCALARRESULT(TRIM(shrev_force_l(6)),TRIM(loadstep),ttime,1,TRIM(gpname),NULL,NULL)
        iv = iv+1
        DO ielem=1,eset%nelem
          DO g=1,ng
            CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*shrev_force_f)
          END DO
        END DO
        CALL GID_ENDRESULT()
      END IF
    END IF

    IF(shrev_momen > 1)THEN
     !CALL cab_gid_bin(1,2,shrev_momen_l(4),shrev_momen_l(5),2,loadstep,ttime,gpname,units=shrev_momen_u)
      CALL GID_BEGINSCALARRESULT(TRIM(shrev_momen_l(5)),TRIM(loadstep),ttime,1,TRIM(gpname),NULL,NULL)
      iv = iv+1
      DO ielem=1,eset%nelem
        DO g=1,ng
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*shrev_momen_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
      IF(ntype > 1)THEN
        CALL GID_BEGINSCALARRESULT(shrev_momen_l(6),TRIM(loadstep),ttime,1,TRIM(gpname),NULL,NULL)
        iv = iv+1
        DO ielem=1,eset%nelem
          DO g=1,ng
            CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*shrev_momen_f)
          END DO
        END DO
        CALL GID_ENDRESULT()
      END IF
    END IF

    IF(shrev_shear > 1)THEN
     !CALL cab_gid_bin(1,2,shrev_shear_l(3),shrev_shear_l(4),1,loadstep,ttime,gpname,units=shrev_shear_u)
      CALL GID_BEGINSCALARRESULT(TRIM(shrev_shear_l(3)),TRIM(loadstep),ttime,1,TRIM(gpname),NULL,NULL)
      iv = iv+1
      DO ielem=1,eset%nelem
        DO g=1,ng
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*shrev_shear_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF(shrev_eqpst > 1)THEN
      CALL cab_gid_bin(1,2,shrev_eqpst_l(3),shrev_eqpst_l(4),1,loadstep,ttime,gpname,units=shrev_eqpst_u)
      iv = iv+1
      DO ielem=1,eset%nelem
        DO g=1,ng
          CALL GID_WRITESCALAR(ielem,MAX(eset%elvar(iv,o(g,ng),ielem),0d0)*shrev_eqpst_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF(shrev_vmise > 1)THEN
      CALL cab_gid_bin(1,2,shrev_vmise_l(3),shrev_vmise_l(4),1,loadstep,ttime,gpname,units=shrev_vmise_u)
      iv = iv+1
      DO ielem=1,eset%nelem
        DO g=1,ng
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*shrev_vmise_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF(shrev_thrat > 1)THEN
      CALL cab_gid_bin(1,2,shrev_thrat_l(3),shrev_thrat_l(4),1,loadstep,ttime,gpname,units=shrev_thrat_u)
      iv = iv+1
      DO ielem=1,eset%nelem
        DO g=1,ng
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*shrev_thrat_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    eset => eset%next
    IF( .NOT.ASSOCIATED (eset) )EXIT
  END DO

  !        print Nodal variables

  IF( shrev_nvarn == 0 )RETURN    !no variables to print, => exit

  nodes => shrev_nodes
  vargs => shrev_vargs
  iv = 0

  IF(MOD(shrev_force,2) == 1)THEN
    CALL cab_gid_bin(nst,1,shrev_force_l(1),shrev_force_l(2),nst,loadstep,ttime,gpname,units=shrev_force_u)
    iv = iv+nst
    is = 0
    DO i=1,npoin
      IF(nodes(i,1) == 0)CYCLE
      is = is + 1
      IF( nst == 1 )THEN
        CALL GID_WRITESCALAR(label(i),vargs(iv,is)*shrev_force_f)
      ELSE
        CALL GID_WRITEVECTOR(label(i),vargs(iv-1,is)*shrev_force_f,vargs(iv,is)*shrev_force_f,0d0)
      END IF
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF(MOD(shrev_momen,2) == 1)THEN
    CALL cab_gid_bin(nst,1,shrev_momen_l(1),shrev_momen_l(2),nst,loadstep,ttime,gpname,units=shrev_momen_u)
    iv = iv+nst
    is = 0
    DO i=1,npoin
      IF(nodes(i,1) == 0)CYCLE
      is = is + 1
      IF( nst == 1 )THEN
        CALL GID_WRITESCALAR(label(i),vargs(iv,is)*shrev_momen_f)
      ELSE
        CALL GID_WRITEVECTOR(label(i),vargs(iv-1,is)*shrev_momen_f,vargs(iv,is)*shrev_momen_f,0d0)
      END IF
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF(MOD(shrev_shear,2) == 1)THEN
   !CALL cab_gid_bin(1,1,shrev_shear_l(1),shrev_shear_l(2),1,loadstep,ttime,gpname,units=shrev_shear_u)
   CALL GID_BEGINSCALARRESULT(TRIM(shrev_shear_l(2)),TRIM(loadstep),ttime,0,NULL,NULL,NULL)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF(nodes(i,1) == 0)CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*shrev_shear_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF(MOD(shrev_eqpst,2) ==  1)THEN
    auxvl = shrev_eqpst_l(1)
    IF( etype == 11 ) auxvl = TRIM(auxvl)//'_1'
    CALL cab_gid_bin(1,1,auxvl,shrev_eqpst_l(2),1,loadstep,ttime,gpname,units=shrev_eqpst_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF(nodes(i,1) == 0)CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*shrev_eqpst_f)
    END DO
    CALL GID_ENDRESULT()
    IF( etype == 11 )THEN
      auxvl = shrev_eqpst_l(1)
      IF( etype == 11 ) auxvl = TRIM(auxvl)//'_n'
      CALL cab_gid_bin(1,1,auxvl,shrev_eqpst_l(2),1,loadstep,ttime,gpname,units=shrev_eqpst_u)
      iv = iv+1
      is = 0
      DO i=1,npoin
        IF(nodes(i,1) == 0)CYCLE
        is = is + 1
        CALL GID_WRITESCALAR(label(i),vargs(iv,is)*shrev_eqpst_f)
      END DO
      CALL GID_ENDRESULT()
    END IF
  END IF

  IF(MOD(shrev_vmise,2) ==  1 .AND. etype == 9 )THEN
    CALL cab_gid_bin(1,1,shrev_vmise_l(1),shrev_vmise_l(2),1,loadstep,ttime,gpname,units=shrev_vmise_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF(nodes(i,1) == 0)CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*shrev_vmise_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF(MOD(shrev_thrat,2) ==  1)THEN
    CALL cab_gid_bin(1,1,shrev_thrat_l(1),shrev_thrat_l(2),1,loadstep,ttime,gpname,units=shrev_thrat_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF(nodes(i,1) == 0)CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*shrev_thrat_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

RETURN
END SUBROUTINE prgib09
