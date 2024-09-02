SUBROUTINE sgaus12()
! Process special Gauss information and compute values for 3-D BST Shell elements
USE data_db,ONLY: npoin, bst, bst_head, bst_nvarn, bst_nodes, bst_vargs, bst_force,      &
                  bst_momen, bst_shear, bst_logst, bst_curva, bst_thrat, bst_eqpst,      &
                  bst_vmise, bst_fldma, bst_wfldFZ, bst_wfldSZ
USE flc_db,ONLY: flc_tp, hpflc, srch_flc, lblflc
IMPLICIT NONE

  !--- Local variables
  INTEGER(kind=4):: iw, iv, is, ipoin, ielem, nelem, inode, lbl, olbl
  INTEGER,POINTER:: nodes(:,:)
  REAL(kind=8):: aux, smtfac, fld, SafZn, ForZn, elstr(2)
  REAL(kind=8),POINTER:: vargs(:,:), ndstr(:,:)
  LOGICAL:: found
  TYPE(bst),POINTER:: eset
  TYPE(flc_tp),POINTER:: flc

  IF (bst_nvarn == 0) RETURN    !no variables to print, => exit

  nodes => bst_nodes   !Local number for bst nodes
  vargs => bst_vargs   !Smoothed variables
  ALLOCATE(ndstr(2,npoin))

  !Smoothed Nodal variables: "Principal strain"
  iw = 0       !initializes pointer to vargs
  IF (MOD(bst_force,2) == 1) iw=iw+3
  IF (MOD(bst_momen,2) == 1) iw=iw+3
  IF (MOD(bst_shear,2) == 1) iw=iw+2
  IF (MOD(bst_logst,2) == 1) iw=iw+3
  IF (MOD(bst_curva,2) == 1) iw=iw+3
  IF (MOD(bst_thrat,2) == 1) iw=iw+1
  IF (MOD(bst_eqpst,2) == 1) iw=iw+2
  IF (MOD(bst_vmise,2) == 1) iw=iw+2
  IF (MOD(bst_fldma,2) == 1) iw=iw+1
  IF (bst_wfldFZ .OR. bst_wfldSZ) THEN
    iw = iw + 2
    is = 0
    DO ipoin=1,npoin
      IF (nodes(ipoin,1) == 0) CYCLE
      is = is + 1
      ndstr(2,ipoin) = vargs(iw-1,is)  !E1
      ndstr(1,ipoin) = vargs(iw,is)    !E2
    END DO
  END IF

  olbl = -1
  eset => bst_head
  DO
    IF( .NOT.ASSOCIATED(eset) )EXIT
    nelem = eset%nelem      !number of elements in the set

    !Gaussian variables
    iv = 0       !initializes pointer
    IF (bst_force > 1) iv=iv+3
    IF (bst_momen > 1) iv=iv+3
    IF (bst_shear > 1) iv=iv+2
    IF (bst_logst > 1) iv=iv+3
    IF (bst_curva > 1) iv=iv+3
    IF (bst_thrat > 1) iv=iv+1
    IF (bst_eqpst > 1) iv=iv+2
    IF (bst_vmise > 1) iv=iv+2
    IF (bst_fldma > 1) iv=iv+1
    IF (bst_wfldFZ .OR. bst_wfldSZ) THEN
      iv = iv + 2
      smtfac  = 1d0/3d0    !Smooth factor
      DO ielem=1,nelem
        IF( lblflc(eset%matno(ielem)) /= olbl)THEN
          olbl = lblflc(eset%matno(ielem))
          CALL srch_flc(hpflc,olbl,found,flc)
        END IF
        elstr(1:2) = 0d0   !Init. smoothing elements output ("second smoothing")
        DO inode=1,3
          ipoin = eset%lnods(inode,ielem)
          elstr(1:2) = elstr(1:2) + ndstr(1:2,ipoin)  !Calculate smoothing values
        END DO
        elstr(1:2) = smtfac*elstr(1:2)   !Second smoothing over elements outputs
        IF (elstr(1) > elstr(2)) THEN    !Can be TRUE in some cases motivated by
          aux = elstr(2)                 !numerical errors during smoothing
          elstr(2) = elstr(1)
          elstr(1) = aux
        END IF
        !Calculate 'Safety Zone' diagram and 'Forming Zone' diagram
        CALL CalDisFLD(flc%npt,flc%cv,elstr(1),elstr(2),fld)
        CALL FLDdg(flc%LmM,flc%LmLS,flc%LmPS,flc%LmT,elstr(1),elstr(2),fld,SafZn,ForZn)
        !Update Safety Zone and Forming Zone
        eset%elvar(iv-1,1,ielem) = ForZn   !assign
        eset%elvar(iv,1,ielem) = SafZn   !assign
      END DO
    END IF

    eset => eset%next
  END DO

  DEALLOCATE(ndstr)

RETURN
END SUBROUTINE sgaus12
