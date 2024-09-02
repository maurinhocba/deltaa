SUBROUTINE prgib05( )
!  print results for GiD for
!  3-D solid elements
USE data_db,ONLY: npoin, loadstep, ttime, label, sol3d, sol3d_head, sol3d_nvarn,         &
    sol3d_nodes, sol3d_vargs, sol3d_stres, sol3d_stres_l, sol3d_logst, sol3d_logst_l,    &
    sol3d_eqpst, sol3d_eqpst_l, sol3d_vmise, sol3d_vmise_l, sol3d_thrat, sol3d_thrat_l,  &
    sol3d_fldma, sol3d_fldma_l, given, udmat , udmat_head,                               &
    sol3d_stres_u, sol3d_stres_f, sol3d_logst_u, sol3d_logst_f,                          &
    sol3d_eqpst_u, sol3d_eqpst_f, sol3d_vmise_u, sol3d_vmise_f, sol3d_thrat_u, sol3d_thrat_f

IMPLICIT NONE

  ! local variables
  TYPE(sol3d),POINTER:: eset
  TYPE (udmat), POINTER :: postd
  INTEGER(kind=4):: i,j,iv,is,nelem,nnode,ielem,ngaus,g,ng,nvar,naux,matno,dim
  INTEGER(kind=4),POINTER:: nodes(:,:)
  REAL(kind=8),POINTER:: vargs(:,:)
  CHARACTER(len=34):: gpname = 'no name yet '
  LOGICAL :: direct
 !     non standart GiD Gauss point positions (INTERNAL)
 INTEGER, PARAMETER :: o(8,6) = RESHAPE( (/ 1, 0, 0, 0, 0, 0, 0, 0,   &  !hexa-tetra-prism
                                            1, 2, 0, 0, 0, 0, 0, 0,   &  !prism
                                            1, 2, 3, 4, 0, 0, 0, 0,   &  !hexa
                                            1, 1, 1, 2, 2, 2, 0, 0,   &  !prism
                                            1, 2, 3, 4, 5, 6, 0, 0,   &  !prism
                                            1, 2, 4, 3, 5, 6, 8, 7/), &  !hexa
                               (/ 8,6 /) )
 INTEGER (kind=4) :: etype
 REAL(kind=8) :: auxil(6,8),aux(6)
  INTERFACE
    INCLUDE 'tr4to8.h'
    INCLUDE 'cab_gid_bin.h'
  END INTERFACE

  eset => sol3d_head


  !        print Gaussian variables
  DO
    gpname = 'GP'//TRIM(eset%sname)
    nelem = eset%nelem
    nnode = eset%nnode
    ngaus = eset%ngaus
    etype = eset%etype
    direct = .TRUE.
    SELECT CASE (ngaus)  !set position in arrya "o" => ng
    CASE(1)              !TETRAHEDRA
      ng = 1             !position in array "O"
    CASE (2)             !PRISMA
      IF( given )THEN    !print in actual positions
        ng = 2
      ELSE               !print in GiD positions
        IF( etype == 4 )THEN ! Hexahedra
          ng = 0           !flag to interpolate
          ngaus = 8
          direct = .FALSE.
        ELSE IF( etype == 5 )THEN ! SPRISM
          ng = -1          !flag to interpolate
          ngaus = 6
          direct = .FALSE.
        ELSE                  ! PRISMA
          ng = 4           !position in array "O"
          ngaus = 6
        END IF
      END IF
    CASE (4)             !HEXAHEDRA
      IF( given )THEN    !print in actual positions
        ng = 3           !position in array "O"
      ELSE               !print in GiD positions
        ng = 6           !position in array "O"
        ngaus = 8
        direct = .FALSE.
      END IF
    CASE (6)             !PRISMA
      ng = 5             !position in array "O"
    CASE (8)             !HEXAHEDRA
      ng = 6             !position in array "O"
    END SELECT

    iv = 0       !initializes pointer

    IF (sol3d_stres > 1) THEN
      CALL cab_gid_bin(3,2,sol3d_stres_l(8),sol3d_stres_l(9),6,loadstep,ttime,gpname,units=sol3d_stres_u)
      iv = iv+6
      DO ielem=1,nelem
        IF( .NOT. direct ) CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,6,ng)
        DO g=1,ngaus
          IF( direct )THEN
            aux = eset%elvar(iv-5:iv,o(g,ng),ielem)*sol3d_stres_f
          ELSE
            aux = auxil(:,g)*sol3d_stres_f
          END IF
          CALL GID_WRITE3DMATRIX(ielem,aux(1),aux(2),aux(3),aux(4),aux(5),aux(6))
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (sol3d_logst > 1) THEN
      CALL cab_gid_bin(2,2,sol3d_logst_l(5),sol3d_logst_l(6),3,loadstep,ttime,gpname,units=sol3d_logst_u)
      iv = iv+3
      DO ielem=1,nelem
        IF( .NOT. direct ) CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,3,ng)
        DO g=1,ngaus
          IF( direct )THEN
            aux(1:3) = eset%elvar(iv-2:iv,o(g,ng),ielem)*sol3d_logst_f
          ELSE
            aux(1:3) = auxil(1:3,g)*sol3d_logst_f
          END IF
          CALL GID_WRITEVECTOR(ielem,aux(1),aux(2),aux(3))
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (sol3d_thrat > 1) THEN
      CALL cab_gid_bin(1,2,sol3d_thrat_l(3),sol3d_thrat_l(4),1,loadstep,ttime,gpname,units=sol3d_thrat_u)
      iv = iv+1
      DO ielem=1,nelem
        IF( .NOT. direct ) CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,1,ng,0d0)
        DO g=1,ngaus
          IF( direct )THEN
            CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*sol3d_thrat_f)
          ELSE
            CALL GID_WRITESCALAR(ielem,auxil(1,g)*sol3d_thrat_f)
          END IF
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (sol3d_eqpst > 1) THEN
      CALL cab_gid_bin(1,2,sol3d_eqpst_l(3),sol3d_eqpst_l(4),1,loadstep,ttime,gpname,units=sol3d_eqpst_u)
      iv = iv+1
      DO ielem=1,nelem
        IF( .NOT. direct ) CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,1,ng,0d0)
        DO g=1,ngaus
          IF( direct )THEN
            CALL GID_WRITESCALAR(ielem,MAX(eset%elvar(iv,o(g,ng),ielem),0d0)*sol3d_eqpst_f)
          ELSE
            CALL GID_WRITESCALAR(ielem,auxil(1,g)*sol3d_eqpst_f)
          END IF
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (sol3d_vmise > 1) THEN
      CALL cab_gid_bin(1,2,sol3d_vmise_l(3),sol3d_vmise_l(4),1,loadstep,ttime,gpname,units=sol3d_vmise_u)
      iv = iv+1
      DO ielem=1,nelem
        IF( .NOT. direct ) CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,1,ng,0d0)
        DO g=1,ngaus
          IF( direct )THEN
            CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem)*sol3d_vmise_f)
          ELSE
            CALL GID_WRITESCALAR(ielem,auxil(1,g)*sol3d_vmise_f)
          END IF
        END DO
      END DO
      CALL GID_ENDRESULT()
    END IF

    IF (sol3d_fldma > 1) THEN
      CALL cab_gid_bin(1,2,sol3d_fldma_l(3),sol3d_fldma_l(4),1,loadstep,ttime,gpname)
      iv = iv+1
      DO ielem=1,nelem
        IF( .NOT. direct ) CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,1,ng,0d0)
        DO g=1,ngaus
          IF( direct )THEN
            CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem))
          ELSE
            CALL GID_WRITESCALAR(ielem,auxil(1,g))
          END IF
        END DO
      END DO
      CALL GID_ENDRESULT()

    END IF
    ! ************** special variables **********************
    !  user defined internal variables
    naux = eset%nstre - 8
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
          CALL cab_gid_bin(1,2,postd%name(1,i),postd%name(2,i),1,loadstep,ttime,gpname)
          iv = iv+1
          DO ielem=1,nelem
            IF( .NOT. direct ) CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,1,ng,0d0)
            DO g=1,ngaus
              IF( direct )THEN
                CALL GID_WRITESCALAR(ielem,eset%elvar(iv,o(g,ng),ielem))
              ELSE
                CALL GID_WRITESCALAR(ielem,auxil(1,g))
              END IF
            END DO
          END DO
          CALL GID_ENDRESULT()

        CASE (1) !vector
          CALL cab_gid_bin(2,2,postd%name(1,i),postd%name(2,i),3, &
                       loadstep,ttime,gpname)
          dim = postd%dim(i)
          j  = iv+1
          iv = iv+dim
          DO ielem=1,nelem
            IF( .NOT. direct ) CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,dim,ng)
            DO g=1,ngaus
               IF( direct )THEN
                 CALL GID_WRITEVECTOR(ielem,eset%elvar(j,o(g,ng),ielem),          &
                 eset%elvar(j+1,o(g,ng),ielem),eset%elvar(iv,o(g,ng),ielem))
               ELSE
                  CALL GID_WRITEVECTOR(ielem,auxil(1,g),auxil(2,g),auxil(3,g))
               END IF
            END DO
           END DO
           CALL GID_ENDRESULT()

        CASE (2) !tensor
          dim = postd%dim(i)
          CALL cab_gid_bin(3,2,postd%name(1,i),postd%name(2,i),6, &
                         loadstep,ttime,gpname)
          j  = iv+1
          iv = iv+dim
          DO ielem=1,nelem
            IF( .NOT. direct ) CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,dim,ng)
            DO g=1,ngaus
              IF( direct )THEN
                CALL GID_WRITE3DMATRIX(ielem,eset%elvar(j,o(g,ng),ielem),         &
                  eset%elvar(j+1,o(g,ng),ielem),eset%elvar(j+2,o(g,ng),ielem),     &
                  eset%elvar(j+3,o(g,ng),ielem),eset%elvar(iv-1,o(g,ng),ielem),     &
                  eset%elvar(iv,o(g,ng),ielem))
              ELSE
                CALL GID_WRITE3DMATRIX(ielem,auxil(1,g),auxil(2,g),auxil(3,g),     &
                                             auxil(4,g),auxil(5,g),auxil(6,g))
              END IF
            END DO
          END DO
          CALL GID_ENDRESULT()

        END SELECT
      END DO
    END IF


    eset => eset%next
    IF (.NOT.ASSOCIATED(eset)) EXIT
  END DO

  !        print Nodal variables
  IF (sol3d_nvarn == 0) RETURN    !no variables to print, => exit

  nodes => sol3d_nodes
  vargs => sol3d_vargs
  iv = 0       !initializes pointer to vargs

  IF (MOD(sol3d_stres,2) == 1) THEN
    CALL cab_gid_bin(3,1,sol3d_stres_l(1),sol3d_stres_l(2),6,loadstep,ttime,gpname,units=sol3d_stres_u)
    iv = iv+6
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      aux = vargs(iv-5:iv,is)*sol3d_stres_f
      CALL GID_WRITE3DMATRIX(label(i),aux(1),aux(2),aux(3),aux(4),aux(5),aux(6))
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(sol3d_logst,2) == 1) THEN
    CALL cab_gid_bin(2,1,sol3d_logst_l(1),sol3d_logst_l(2),3,loadstep,ttime,gpname,units=sol3d_logst_u)
    iv = iv+3
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      aux(1:3) = vargs(iv-2:iv,is)*sol3d_logst_f
      CALL GID_WRITEVECTOR(label(i),aux(1),aux(2),aux(3))
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(sol3d_thrat,2) == 1) THEN
    CALL cab_gid_bin(1,1,sol3d_thrat_l(1),sol3d_thrat_l(2),1,loadstep,ttime,gpname,units=sol3d_thrat_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*sol3d_thrat_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(sol3d_eqpst,2) == 1) THEN
    CALL cab_gid_bin(1,1,sol3d_eqpst_l(1),sol3d_eqpst_l(2),1,loadstep,ttime,gpname,units=sol3d_eqpst_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*sol3d_eqpst_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(sol3d_vmise,2) == 1) THEN
    CALL cab_gid_bin(1,1,sol3d_vmise_l(1),sol3d_vmise_l(2),1,loadstep,ttime,gpname,units=sol3d_vmise_u)
    iv = iv+1
    is = 0
    DO i=1,npoin
      IF (nodes(i,1) == 0) CYCLE
      is = is + 1
      CALL GID_WRITESCALAR(label(i),vargs(iv,is)*sol3d_vmise_f)
    END DO
    CALL GID_ENDRESULT()
  END IF

  IF (MOD(sol3d_fldma,2) == 1) THEN
    CALL cab_gid_bin(1,1,sol3d_fldma_l(1),sol3d_fldma_l(2),1,loadstep,ttime,gpname)
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
END SUBROUTINE prgib05
