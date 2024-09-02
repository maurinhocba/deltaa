SUBROUTINE prgib08( )
!  print results for GiD for beame elements
USE data_db,ONLY: npoin, loadstep, ttime, label, beame, beame_head, beame_nvarn,          &
    beame_nodes, beame_vargs, beame_force, beame_force_l, beame_momen, beame_momen_l,     &
    beame_eqpst, beame_eqpst_l, beame_arrat, beame_arrat_l, beame_force_u, beame_force_f, &
    beame_momen_u, beame_momen_f, beame_eqpst_u, beame_eqpst_f, beame_arrat_u, beame_arrat_f
IMPLICIT NONE

  TYPE(beame),POINTER:: eset
  INTEGER(kind=4):: ielem,g,iv,is,i,ng
  INTEGER(kind=4),POINTER:: nodes(:,:)
  REAL(kind=8),POINTER:: vargs(:,:)
  CHARACTER(len=34):: gpname
  INTEGER(kind=4),PARAMETER:: o(3,3) = RESHAPE( (/ 1, 0, 0,  1, 2, 0, 1, 3, 2 /), &
                                                (/ 3, 3 /) )
  REAL(kind=8):: aux(3)

  INTERFACE
    INCLUDE 'cab_gid_bin.h'
  END INTERFACE

  eset => beame_head          !point to first set

  !        print Gaussian variables
  DO                          !loop for all the sets
    gpname = 'GP'//TRIM(eset%sname)   !gauss point name string "GPsname"
    iv = 0  !initializes pointer to Gauss values
    ng = eset%ngaus

    IF(beame_force > 1)THEN
      CALL cab_gid_bin(2,2,beame_force_l(5),beame_force_l(6),3,loadstep,ttime,gpname,units=beame_force_u)
      iv = iv+3
      DO ielem=1,eset%nelem
        DO g=1,ng
          aux = eset%elvar(iv-2:iv,o(g,ng),ielem)*beame_force_f
          CALL GID_WRITEVECTOR(ielem,aux(1),aux(2),aux(3))
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF(beame_momen > 1)THEN
      CALL cab_gid_bin(2,2,beame_momen_l(5),beame_momen_l(6),3,loadstep,ttime,gpname,units=beame_momen_u)
      iv = iv+3
      DO ielem=1,eset%nelem
        DO g=1,ng
          aux = eset%elvar(iv-2:iv,o(g,ng),ielem)*beame_momen_f
          CALL GID_WRITEVECTOR(ielem,aux(1),aux(2),aux(3))
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF(beame_eqpst > 1)THEN
      CALL cab_gid_bin(1,2,beame_eqpst_l(3),beame_eqpst_l(4),1,loadstep,ttime,gpname,units=beame_eqpst_u)
      iv = iv+1
      DO ielem=1,eset%nelem
        DO g=1,ng
          CALL GID_WRITESCALAR(ielem,MAX(eset%elvar(iv,o(g,ng),ielem),0d0)*beame_eqpst_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    !IF(beame_vmise > 1)THEN
    ! CALL cab_gid_bin(1,2,beame_vmise_l(3),beame_vmise_l(4),1,loadstep,ttime,gpname,units=beame_vmise_u)
    !  iv = iv+1
    !  DO ielem=1,eset%nelem
    !    DO g=1,ng
    !      CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*beame_vmise_f)
    !    END DO
    !  END DO
    ! CALL GID_ENDRESULT()
    !END IF

    IF(beame_arrat > 1)THEN
      CALL cab_gid_bin(1,2,beame_arrat_l(3),beame_arrat_l(4),1,loadstep,ttime,gpname,units=beame_arrat_u)
      iv = iv+1
      DO ielem=1,eset%nelem
        DO g=1,ng
          CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*beame_arrat_f)
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    eset => eset%next
    IF( .NOT.ASSOCIATED (eset) )EXIT
  END DO

  !        print Nodal variables

  IF( beame_nvarn == 0 )RETURN    !no variables to print, => exit

  nodes => beame_nodes
  vargs => beame_vargs
  iv = 0

  IF(MOD(beame_force,2) == 1)THEN
    CALL cab_gid_bin(2,1,beame_force_l(1),beame_force_l(2),3,loadstep,ttime,gpname,units=beame_force_u)
    iv = iv+3
    is = 0
    DO i=1,npoin
      IF(nodes(i,1) == 0)CYCLE
      is = is + 1
      aux = vargs(iv-2:iv,is)*beame_force_f
      CALL GID_WRITEVECTOR(label(i),aux(1),aux(2),aux(3))
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF(MOD(beame_momen,2) == 1)THEN
    CALL cab_gid_bin(2,1,beame_momen_l(1),beame_momen_l(2),3,loadstep,ttime,gpname,units=beame_momen_u)
    iv = iv+3
    is = 0
    DO i=1,npoin
      IF(nodes(i,1) == 0)CYCLE
      is = is + 1
      aux = vargs(iv-2:iv,is)*beame_momen_f
      CALL GID_WRITEVECTOR(label(i),aux(1),aux(2),aux(3))
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF(MOD(beame_eqpst,2) ==  1)THEN
    CALL cab_gid_bin(1,1,beame_eqpst_l(1),beame_eqpst_l(2),1,loadstep,ttime,gpname,units=beame_eqpst_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF(nodes(i,1) == 0)CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*beame_eqpst_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  !IF(MOD(beame_vmise,2) ==  1)THEN
  ! CALL cab_gid_bin(1,1,beame_vmise_l(1),beame_vmise_l(2),1,loadstep,ttime,gpname,units=beame_vmise_u))
  !  iv = iv+1
  !  is = 0
  !  DO i=1,npoin
  !   IF(nodes(i,1) == 0)CYCLE
  !    is = is + 1
  !    CALL GID_WRITESCALAR(label(i),vargs(iv,is)*beame_vmise_f)
  !  END DO
  !  CALL GID_ENDRESULT()
  !END IF

  IF(MOD(beame_arrat,2) ==  1)THEN
    CALL cab_gid_bin(1,1,beame_arrat_l(1),beame_arrat_l(2),1,loadstep,ttime,gpname,units=beame_arrat_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF(nodes(i,1) == 0)CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*beame_arrat_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

RETURN
END SUBROUTINE prgib08
