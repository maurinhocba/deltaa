MODULE solv_db
  !stiffness matrices, iterative vectors & control parameters
  IMPLICIT NONE
  INTEGER (kind=4) :: &
    nsymm, & ! flag to use symmetric or non-symmetric matrices
    nbuck, & ! frequency for buckling analysis
    nralg, & ! encrypted advance algorithm parameters
    ncdis, & ! node & DOF controlling advance
    ecdis, & ! equation controlling advance
    scdis    ! node & DOF controlling switching

  REAL (kind=8) :: &
    arcln,  & ! generalized arc length
    dlamb,  & ! increment in load parameter
    lambd=0d0,  & ! actual load parameter
    buckl,  & ! last computed buckling load parameter
    toler,  & ! convergence tolerance
    tol1 ,  & ! convergence tolerance for extended analysis
    alph1,  & ! alpha parameter for extended systems
    gamm1,  & ! gamma parameter for extended systems
    arcl1     ! arc-length parameter for switching paths

  LOGICAL :: linear  !.TRUE. for Linear switching along eigenvalue
  LOGICAL :: lsearch !.TRUE. for Linear search


  REAL (kind=8), POINTER :: &
    stiff(:,:), & !(maxa,?) stiffness matrix and derivative
    displ(:),   & !(neq) incremental displacement
    ddisp(:),   & !(neq) iterative incremental displacement
    gvect(:),   & !(neq) assembled internal forces
    disax(:,:)    !(neq+1,7) auxiliar vectors
                  !(1) incremental load vector to compute ==> Ut
                  !(2) step incremental displacements
                  !(4) forces due to incremental prescribed displacements ==> buckling mode
                  !(5:7) last 3 steps increments to predict a quadratic approach

                  !(1) du_1 = K^-1 f  = U_t
                  !(2) du_2
                  !(3) du_3
                  !(4) buckling mode
                  !(5)  - Du(Kt.Phi)* du_1
                  !(6)  - Du(Kt.Phi)* du_2
                  !(7)  - Du(Kt.Phi)* du_3
  INTEGER (kind=4) :: neigen,mbite
  REAL(kind=8), ALLOCATABLE :: r(:,:), eigv(:), ar(:), br(:),vec(:,:), d(:), rtolv(:)

CONTAINS

  SUBROUTINE dump_solv (neq)
    IMPLICIT NONE
    INTEGER (kind=4), INTENT(IN) :: neq

    INTEGER (kind=4) :: i,j !,nds

    WRITE(50) nsymm, nbuck, nralg, ncdis, ecdis, scdis, linear
    WRITE(50) arcln, dlamb, lambd, buckl, toler, tol1,  alph1, gamm1, arcl1, neigen, mbite
    !nds = 1
    !IF(nsymm > 0) nds = nds+1
    !IF(nbuck > 0) nds = nds+1
    !WRITE(50) ((stiff(i,j),i=1,maxa),j=1,nds)
    !WRITE(50) (gvect(i),i=1,neq)
    !WRITE(50) (ddisp(i),i=1,neq)
    WRITE(50) (displ(i),i=1,neq)
    WRITE(50) ((disax(i,j),i=1,neq+1),j=1,7)

    RETURN
  END SUBROUTINE dump_solv

  SUBROUTINE rest_solv (maxa,neq)
    IMPLICIT NONE
    INTEGER (kind=4), INTENT(IN) :: maxa,neq

    INTEGER (kind=4) :: i,j ,nds

    READ(51) nsymm, nbuck, nralg, ncdis, ecdis, scdis, linear
    READ(51) arcln, dlamb, lambd, buckl, toler, tol1,  alph1, gamm1, arcl1, neigen, mbite
    nds = 1
    IF(nsymm > 0) nds = nds+1
    IF(nbuck > 0) nds = nds+1
    ALLOCATE( stiff(maxa,nds), gvect(neq), displ(neq), ddisp(neq))
    !READ(51) ((stiff(i,j),i=1,maxa),j=1,nds)
    !READ(51) (gvect(i),i=1,neq)
    !READ(51) (ddisp(i),i=1,neq)
    READ(51) (displ(i),i=1,neq)
    ALLOCATE( disax(neq+1,7) )
    READ(51) ((disax(i,j),i=1,neq+1),j=1,7)

    RETURN
  END SUBROUTINE rest_solv

END MODULE solv_db
