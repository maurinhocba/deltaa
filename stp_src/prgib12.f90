SUBROUTINE prgib12( )
!  print results for GiD for
!  3-D BST Shell elements
USE data_db,ONLY: npoin, loadstep, ttime, label, bst, bst_head, bst_nvarn,     &
    bst_nodes, bst_vargs, bst_force, bst_force_l, bst_momen, bst_momen_l,      &
    bst_shear, bst_shear_l, &
    bst_logst, bst_logst_l, bst_curva, bst_curva_l, bst_thrat, bst_thrat_l,    &
    bst_eqpst, bst_eqpst_l, bst_vmise, bst_vmise_l, bst_fldma, bst_fldma_l,    &
    bst_wfldFZ, bst_fldFZ_l, bst_wfldSZ, bst_fldSZ_l,                          &
    udmat,udmat_head, rrt_fr, rrt_sf, bst_force_u, bst_force_f, bst_momen_u, bst_momen_f,   &
    bst_shear_u, bst_shear_f, &
    bst_logst_u, bst_logst_f, bst_curva_u, bst_curva_f, bst_thrat_u, bst_thrat_f,    &
    bst_eqpst_u, bst_eqpst_f, bst_vmise_u, bst_vmise_f

IMPLICIT NONE

  ! local variables
  TYPE(bst),POINTER:: eset
  TYPE (udmat), POINTER :: postd
  INTEGER(kind=4):: ielem,i,iv,nelem,nnode,nstre,ngaus,is,dim,naux,matno,nvar
  INTEGER(kind=4),POINTER:: nodes(:,:)
  REAL(kind=8),POINTER:: vargs(:,:)
  REAL(kind=8) :: aux(3)
  CHARACTER(len=34):: gpname
  CHARACTER(len=17):: saux,caux

  INTERFACE
    INCLUDE 'cab_gid_bin.h'
  END INTERFACE

  eset => bst_head

  !        print Gaussian variables
  DO
    gpname = 'GP'//TRIM(eset%sname)
    nelem = eset%nelem
    nnode = eset%nnode
    nstre = eset%nstre
    ngaus = eset%ngaus
    iv = 0       !initializes pointer

    IF (bst_force > 1) THEN
      CALL cab_gid_bin(2,2,bst_force_l(5),bst_force_l(6),3,loadstep,ttime,gpname,units=bst_force_u)
      iv = iv+3
      DO ielem=1,nelem
        aux = eset%elvar(iv-2:iv,1,ielem)*bst_force_f
        CALL GID_WRITEVECTOR(ielem,aux(1),aux(2),aux(3))
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (bst_momen > 1) THEN
      CALL cab_gid_bin(2,2,bst_momen_l(5),bst_momen_l(6),3,loadstep,ttime,gpname,units=bst_momen_u)
      !CALL cab_gid_bin(4,2,bst_momen_l(5),bst_momen_l(6),3,loadstep,ttime,gpname,units=bst_momen_u)
      !CALL cab_gid_bin(4,2,bst_momen_l(5),bst_momen_l(6),4,loadstep,ttime,gpname,units=bst_momen_u)
      iv = iv+3
      DO ielem=1,nelem
        aux = eset%elvar(iv-2:iv,1,ielem)*bst_momen_f
        CALL GID_WRITEVECTOR(ielem,aux(1),aux(2),aux(3))
        !CALL GiD_Write2DMatrix(ielem,aux(1),aux(2),aux(3))
        !CALL GiD_WritePlainDefMatrix(ielem,aux(1),aux(2),aux(3),0d0)

      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (bst_shear > 1) THEN
      CALL cab_gid_bin(2,2,bst_shear_l(4),bst_shear_l(5),2,loadstep,ttime,gpname,units=bst_shear_u)
      iv = iv+2
      DO ielem=1,nelem
        aux(1:2) = eset%elvar(iv-1:iv,1,ielem)*bst_shear_f
        CALL GID_WRITEVECTOR(ielem,aux(1),aux(2),0d0)
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (bst_logst > 1) THEN
      CALL cab_gid_bin(2,2,bst_logst_l(5),bst_logst_l(6),3,loadstep,ttime,gpname,units=bst_logst_u)
      iv = iv+3
      DO ielem=1,nelem
        aux = eset%elvar(iv-2:iv,1,ielem)*bst_logst_f
        CALL GID_WRITEVECTOR(ielem,aux(1),aux(2),aux(3))
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (bst_curva > 1) THEN
      CALL cab_gid_bin(2,2,bst_curva_l(5),bst_curva_l(6),3,loadstep,ttime,gpname,units=bst_curva_u)
      iv = iv+3
      DO ielem=1,nelem
        aux = eset%elvar(iv-2:iv,1,ielem)*bst_curva_f
        CALL GID_WRITEVECTOR(ielem,aux(1),aux(2),aux(3))
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (bst_thrat > 1) THEN
      CALL cab_gid_bin(1,2,bst_thrat_l(3),bst_thrat_l(4),1,loadstep,ttime,gpname,units=bst_thrat_u)
      iv = iv+1
      DO ielem=1,nelem
        CALL GID_WRITESCALAR(ielem,eset%elvar(iv,1,ielem)*bst_thrat_f)
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (bst_eqpst > 1) THEN
      saux = TRIM(bst_eqpst_l(3))//'_1'
      caux = TRIM(bst_eqpst_l(4))//'_1'
      CALL cab_gid_bin(1,2,saux,caux,1,loadstep,ttime,gpname,units=bst_eqpst_u)
      iv = iv+1
      DO ielem=1,nelem
        CALL GID_WRITESCALAR(ielem,MAX(eset%elvar(iv,1,ielem),0d0)*bst_eqpst_f)
      END DO
      CALL GID_ENDRESULT()

      saux = TRIM(bst_eqpst_l(3))//'_N'
      caux = TRIM(bst_eqpst_l(4))//'_N'
      CALL cab_gid_bin(1,2,saux,caux,1,loadstep,ttime,gpname,units=bst_eqpst_u)
      iv = iv+1
      DO ielem=1,nelem
        CALL GID_WRITESCALAR(ielem,MAX(eset%elvar(iv,1,ielem),0d0)*bst_eqpst_f)
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (bst_vmise > 1) THEN
      saux = TRIM(bst_vmise_l(3))//'_1'
      caux = TRIM(bst_vmise_l(4))//'_1'
      CALL cab_gid_bin(1,2,saux,caux,1,loadstep,ttime,gpname,units=bst_vmise_u)
      iv = iv+1
      DO ielem=1,nelem
        CALL GID_WRITESCALAR(ielem,eset%elvar(iv,1,ielem)*bst_vmise_f)
      END DO
      CALL GID_ENDRESULT()

      saux = TRIM(bst_vmise_l(3))//'_N'
      caux = TRIM(bst_vmise_l(4))//'_N'
      CALL cab_gid_bin(1,2,saux,caux,1,loadstep,ttime,gpname,units=bst_vmise_u)
      iv = iv+1
      DO ielem=1,nelem
        CALL GID_WRITESCALAR(ielem,eset%elvar(iv,1,ielem)*bst_vmise_f)
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (bst_fldma > 1) THEN
      CALL cab_gid_bin(1,2,bst_fldma_l(3),bst_fldma_l(4),1,loadstep,ttime,gpname)
      iv = iv+1
      DO ielem=1,nelem
        CALL GID_WRITESCALAR(ielem,eset%elvar(iv,1,ielem))
      END DO
      CALL GID_ENDRESULT()
    END IF
    IF (bst_wfldFZ) THEN
      CALL cab_gid_bin(1,2,bst_fldFZ_l(1),bst_fldFZ_l(1),1,loadstep,ttime,gpname,TRIM(rrt_fr%name))
      iv = iv+1
      DO ielem=1,nelem
        CALL GID_WRITESCALAR(ielem,eset%elvar(iv,1,ielem))
      END DO
      CALL GID_ENDRESULT()
    END IF
    IF (bst_wfldSZ) THEN
      CALL cab_gid_bin(1,2,bst_fldSZ_l(1),bst_fldSZ_l(1),1,loadstep,ttime,gpname,TRIM(rrt_sf%name))
      iv = iv+1
      DO ielem=1,nelem
        CALL GID_WRITESCALAR(ielem,eset%elvar(iv,1,ielem))
      END DO
      CALL GID_ENDRESULT()
    END IF

   !  user defined internal variables
   naux = eset%nstre - 13
   IF( naux > 0 )THEN
     matno = eset%matno(1)
     postd => udmat_head
     DO
       IF( postd%matno == matno )EXIT
       postd => postd%next
     END DO
     nvar = postd%nvar
     DO i=1,nvar
       SELECT CASE (postd%type(i))
       CASE (0) !scalar
         CALL cab_gid_bin(1,2,postd%name(1,i),postd%name(2,i),1,            &
                      loadstep,ttime,gpname)
         iv = iv+1
         DO ielem=1,nelem
           CALL GID_WRITESCALAR(ielem,eset%elvar(iv,1,ielem))
         END DO
       CASE (1) !vector
         dim = postd%dim(i) !must be three
         CALL cab_gid_bin(2,2,postd%name(1,i),postd%name(2,i),dim, &
                      loadstep,ttime,gpname)
         iv = iv+dim        !must be three
         DO ielem=1,nelem
           CALL GID_WRITEVECTOR(ielem,eset%elvar(iv-2,1,ielem),                   &
               eset%elvar(iv-1,1,ielem),eset%elvar(iv,1,ielem))
         END DO
       CASE (2) !tensor
         dim = postd%dim(i)  !must be three (number of components)
         iv = iv+dim           !three components included in a vector
         !CALL cab_gid_bin(4,2,postd%name(1,i),postd%name(2,i),dim,loadstep,ttime,gpname)
         CALL cab_gid_bin(2,2,postd%name(1,i),postd%name(2,i),3,loadstep,ttime,gpname)
         DO ielem=1,nelem
           CALL GID_WRITEVECTOR(ielem,eset%elvar(iv-2,1,ielem),                   &
               eset%elvar(iv-1,1,ielem),eset%elvar(iv,1,ielem))
           !CALL GiD_Write2DMatrix(ielem,eset%elvar(iv-2,1,ielem),                   &
           !     eset%elvar(iv-1,1,ielem),eset%elvar(iv,1,ielem))

         END DO
       END SELECT
       CALL GID_ENDRESULT()
     END DO
   END IF

    eset => eset%next
    IF (.NOT.ASSOCIATED(eset)) EXIT
  END DO

  !        print Nodal variables

  IF (bst_nvarn == 0) RETURN    !no variables to print, => exit

  nodes => bst_nodes
  vargs => bst_vargs
  iv = 0       !initializes pointer to vargs

  IF (MOD(bst_force,2) == 1) THEN
    CALL cab_gid_bin(2,1,bst_force_l(1),bst_force_l(2),3,loadstep,ttime,gpname,units=bst_force_u)
    iv = iv+3
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      aux = vargs(iv-2:iv,is)*bst_force_f
      CALL GID_WRITEVECTOR(label(i),aux(1),aux(2),aux(3))
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(bst_momen,2) == 1) THEN
    CALL cab_gid_bin(2,1,bst_momen_l(1),bst_momen_l(2),3,loadstep,ttime,gpname,units=bst_momen_u)
    iv = iv+3
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      aux = vargs(iv-2:iv,is)*bst_momen_f
      CALL GID_WRITEVECTOR(label(i),aux(1),aux(2),aux(3))
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(bst_shear,2) == 1) THEN
    CALL cab_gid_bin(2,1,bst_shear_l(1),bst_shear_l(2),2,loadstep,ttime,gpname,units=bst_shear_u)
    iv = iv+2
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      aux(1:2) = vargs(iv-1:iv,is)*bst_momen_f
      CALL GID_WRITEVECTOR(label(i),aux(1),aux(2),0d0)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(bst_logst,2) == 1) THEN
    CALL cab_gid_bin(2,1,bst_logst_l(1),bst_logst_l(2),3,loadstep,ttime,gpname,units=bst_logst_u)
    iv = iv+3
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      aux = vargs(iv-2:iv,is)*bst_logst_f
      CALL GID_WRITEVECTOR(label(i),aux(1),aux(2),aux(3))
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(bst_curva,2) == 1) THEN
    CALL cab_gid_bin(2,1,bst_curva_l(1),bst_curva_l(2),3,loadstep,ttime,gpname,units=bst_curva_u)
    iv = iv+3
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      aux = vargs(iv-2:iv,is)*bst_curva_f
      CALL GID_WRITEVECTOR(label(i),aux(1),aux(2),aux(3))
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(bst_thrat,2) == 1) THEN
    CALL cab_gid_bin(1,1,bst_thrat_l(1),bst_thrat_l(2),1,loadstep,ttime,gpname,units=bst_thrat_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*bst_thrat_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(bst_eqpst,2) == 1) THEN
    saux = TRIM(bst_eqpst_l(1))//'_1'
    caux = TRIM(bst_eqpst_l(2))//'_1'
    CALL cab_gid_bin(1,1,saux,caux,1,loadstep,ttime,gpname,units=bst_eqpst_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*bst_eqpst_f)
    END DO
    CALL GID_ENDRESULT()

    saux = TRIM(bst_eqpst_l(1))//'_N'
    caux = TRIM(bst_eqpst_l(2))//'_N'
    CALL cab_gid_bin(1,1,saux,caux,1,loadstep,ttime,gpname,units=bst_eqpst_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*bst_eqpst_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(bst_vmise,2) == 1) THEN
    saux = TRIM(bst_vmise_l(1))//'_1'
    caux = TRIM(bst_vmise_l(2))//'_1'
    CALL cab_gid_bin(1,1,saux,caux,1,loadstep,ttime,gpname,units=bst_vmise_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*bst_vmise_f)
    END DO
    CALL GID_ENDRESULT()

    saux = TRIM(bst_vmise_l(1))//'_N'
    caux = TRIM(bst_vmise_l(2))//'_N'
    CALL cab_gid_bin(1,1,saux,caux,1,loadstep,ttime,gpname,units=bst_vmise_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*bst_vmise_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(bst_fldma,2) == 1) THEN
    CALL cab_gid_bin(1,1,bst_fldma_l(1),bst_thrat_l(2),1,loadstep,ttime,gpname)
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
END SUBROUTINE prgib12
