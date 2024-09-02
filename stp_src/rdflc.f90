 SUBROUTINE rdflc()
 !=======================================================================
 ! Read FLC curves
 !=======================================================================
 USE flc_db  !,ONLY: nflc, flc_tp, hpflc,tpflc
 USE param_db,ONLY: mlng
 IMPLICIT NONE

   !Local variables
   INTEGER(kind=4):: i,j
   REAL(kind=8):: pi
   CHARACTER(len=mlng):: chdum
   TYPE(flc_tp),POINTER:: flc

   pi = 4d0*DATAN(1d0)
   !Read used FLC curves
   DO i=1,nflc
     CALL new_flc(flc)
     READ(17) flc%lbl,flc%sttyp,flc%npt                !Number of points of the FLC curve
     READ(17) flc%LmM, flc%LmPS, flc%LmLS, flc%LmT     !read  strains limits
     ALLOCATE(flc%cv(2,flc%npt))
     READ(17) (flc%cv(1:2,j),j=1,flc%npt)              !read  strains coordinates of the FLC
     CALL add_flc(flc,hpflc,tpflc)
   END DO

 RETURN
 END SUBROUTINE rdflc
