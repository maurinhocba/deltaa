 SUBROUTINE modp20(dmatx,ntype,young,poiss,prop,efpst,defps, &
                   bbar,elast,newmt)
 !*********************************************************************
 !
 !**** this SUBROUTINE evaluates the elastic d-matrix (upper part only)
 !
 !
 !**********************************************************************
 IMPLICIT NONE
 LOGICAL, INTENT(IN) :: bbar,elast
 LOGICAL, INTENT(IN OUT) :: newmt
 INTEGER (kind=4), INTENT(IN) :: ntype
 REAL (kind=8), INTENT(IN) :: prop(5),young,poiss,defps,efpst
 REAL (kind=8), INTENT(OUT) :: dmatx(4,4)

 END SUBROUTINE modp20
