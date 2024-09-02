 MODULE ctrl_db
 ! global control paramters
 USE param_db, ONLY: midn, mttl, mnam
 IMPLICIT NONE

   CHARACTER (len=midn) :: ptype='IMPLCT' ! problem type
   CHARACTER (len=mttl) :: text           ! problem title

   ! size (dimension) parameters
   INTEGER (kind=4) ::  &
             ndime,        & !problem dimension
             ndofn,        & !number of DOFs per node
             nrotd,        & !number of DOFs per node, exclusing extra DOFs
             neulr,        & !if local systems used
             ntype,        & !problem type for 2-D 1:plane stress 2:plane strain 3:axilsymmetric
             nload=0,      & !number of load sets
             npoin=0,      & !number of points in the mesh
             npoio=0         !number of points in previous mesh

   INTEGER (kind=4) :: neq,    & !number of active equations
                       maxa      !number of elements in stiffness matrix
   INTEGER (kind=4) :: neqt=0    !number of active temperature equations

   LOGICAL       ::  inverse= .FALSE.      ! compute original configuration from deformed one
   ! stepping control parameters
   INTEGER (kind=4)  :: istep=0,   & ! step number
                        nstep,     & ! number of steps in the analysis
                        nstra=1      ! counter of strategies read from data file

   ! time integration parameters
   REAL (kind=8) :: &
                    begtm,     & ! Time of Begining of the new strategy
                    endtm,     & ! Final Time
                    dtime,     & ! time step
                    ttime=0d0, & ! accumulated time of the process
                    mscal        ! Mass SCALing factor (affecting dtime)

   ! contact/drawbead codes
   INTEGER (kind=4), PARAMETER :: mconts = 8 ! number of impl. contact algors.
                                             ! 1-6: contacts, 7-8: drawbeads
   INTEGER (kind=4) :: nconp(mconts) = 0,  & ! numbers of contact pairs
                       numct = 0,          & ! number of contact algors. used
                       ctype(mconts) = 0     ! contact types (algors.) used
   LOGICAL :: bottom = .FALSE.,            & ! bottom surface is used
              top    = .FALSE.               ! top surface is used

 !   temperature analysis parameters
   LOGICAL       :: &
                    itemp= .FALSE.      ! compute deformations due to temperature mod

   INTEGER (kind=4) :: &
                    ndoft      !number of temperature DOFs per node

   ! other
   LOGICAL :: echo_chnode = .TRUE.
   LOGICAL :: initial_displacements = .FALSE.
   LOGICAL :: initial_rotations = .FALSE.

  ! Information association to dynamic analysis
  INTEGER (kind=4) :: &
    ndyna    ! analysis type

  REAL (kind=8) :: &
    alpha, &  ! alpha coefficient for Newmark integration scheme
    beta,  &  ! beta coefficient for Newmark integration scheme
    gamma, &  ! gamma coefficient for Newmark integration scheme
    ccm,   &  ! mass coefficient for proportional damping
    cck       ! stiffness coefficient for proportional damping

  INTEGER (kind=4)   ::  memo        !size of auxiliar array for input


 CONTAINS

    SUBROUTINE dump_ctrl
    !   dumps database of control parameters
    IMPLICIT NONE

      WRITE (50,ERR=9999) ndime, ndofn, nrotd, neulr, ntype, nload, npoin, npoio
      WRITE (50,ERR=9999) neq, maxa, memo, neqt
      WRITE (50,ERR=9999) ptype, text
      WRITE (50,ERR=9999) istep, nstep, nstra
      WRITE (50,ERR=9999) begtm, endtm, dtime, ttime, mscal
      WRITE (50,ERR=9999) nconp, numct, ctype, bottom, top, itemp, ndoft
      WRITE (50,ERR=9999) ndyna
      WRITE (50,ERR=9999) alpha, beta, gamma, ccm, cck

    RETURN
    9999 CALL runen2('')
    END SUBROUTINE dump_ctrl

    SUBROUTINE rest_ctrl

      !   restores database of control parameters

      READ (51) ndime, ndofn, nrotd, neulr, ntype, nload, npoin, npoio
      READ (51) neq, maxa, memo, neqt
      READ (51) ptype, text
      READ (51) istep, nstep, nstra
      READ (51) begtm, endtm, dtime, ttime, mscal
      READ (51) nconp, numct, ctype, bottom, top, itemp, ndoft
      READ (51) ndyna
      READ (51) alpha, beta, gamma, ccm, cck

    END SUBROUTINE rest_ctrl

 END MODULE ctrl_db
