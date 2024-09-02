 SUBROUTINE modp18(dmatx,young,poiss,prop,efpst,defps, &
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
 REAL (kind=8), INTENT(IN) :: prop(5),young,poiss,defps,efpst
 REAL (kind=8), INTENT(OUT) :: dmatx(6,6)

 END SUBROUTINE modp18
