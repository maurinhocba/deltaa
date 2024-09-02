  SUBROUTINE streb2(stran,sigma,cm,prop,chi,varin,ierr,newmt,plast,elast)

   ! computes stresses for plane stress model

     IMPLICIT NONE
   ! dummy arguments
   INTEGER (kind=4), INTENT(OUT) :: ierr     !flag to indicate error
   REAL (kind=8), INTENT(IN) :: stran(:),  & !(2) log strains
                                cm(:),     & !(3) elasticity coeffs.
                                chi(:),    & !(12) Hill coeffs.
                                prop(:)      !(5) material properties
   REAL (kind=8), INTENT(IN OUT) :: varin(:) !(3) internal variables
   REAL (kind=8), INTENT(OUT) :: sigma(:)    !(2) stress
   LOGICAL, INTENT(IN) :: plast,elast
   LOGICAL, INTENT(IN OUT) :: newmt

  END SUBROUTINE streb2
