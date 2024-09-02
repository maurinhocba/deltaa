MODULE npo_db
  USE esets_db, ONLY: nel
  ! nodal point information
  IMPLICIT NONE

  REAL(kind=8):: scale   !Scale factor for nodal coordinates

  INTEGER (kind=4),POINTER :: &
            ifpre(:,:),   & !(ndofn,npoin) Associated equation to each DOF's
            iffix(:),     & !(npoin)       rotation constraint (shell elements)
            iftmp(:,:),   & !(ndoft,npoin) Associated equation to each Temp DOF
            label(:),     & !(npoin)       node labels
            loass(:),     & !(nload)       curve assigment for each load set
            oldlb(:),     & !(npoio)       node labels in previous strategy
            id(:,:),      & !(ndime,npoin) array to assemble contact forces
            ifact(:),     & !(npoin)       factors to compute coorb & coort
            cpx(:,:)        !(3,npoin) mid-side control points

  REAL (kind=8),POINTER :: &
          coora(:,:),     & !(ndime,npoin) actual coordinates
          coorc(:,:),     & !(ndime,npoin) coordinates at previous step
          coord(:,:),     & !(ndime,npoin) original coordinates
          coori(:,:),     & !(ndime,npoin) coordinates a previous iteration
          coors(:,:),     & !(ndime,npoin) coordinates at beginning of stage
          coor1(:,:),     & !(ndime,npoin) auxiliar coordinates
          coorb(:,:),     & !(ndime,npoin) coordinates of bottom surface
          coort(:,:),     & !(ndime,npoin) coordinates of top surface
          emass(:,:),     & !(nrotd,npoin) nodal mass
          euler(:,:),     & !(neulr,npoin) local coordinate system (actual)
          eule0(:,:),     & !(1-3,npoin) local coordinate system (initial)
          eulei(:,:),     & !(1-3,npoin) local coordinate system (initial)
          locsy(:,:),     & !(neulr,npoin) local coordinate system (actual)
          eule1(:,:),     & !(neulr,npoin) local coordinate system auxiliar
          locs1(:,:),     & !(neulr,npoin) local coordinate system auxiliar
          fcont(:,:),     & !(ndime,npoin) contact forces
          force(:,:),     & !(neq,nload+1) external forces
          loadv(:,:,:),   & !(ndofn,npoin,nload) external forces
          resid(:,:),     & !(ndofn,npoin) internal forces
          psic(:,:),      & !(2,npoin) psi functions (converged)
          psii(:,:),      & !(2,npoin) psi functions (iterative)
          psia(:,:),      & !(2,npoin) psi functions (present)
          dtemp(:,:),     & !(ndoft,npoin) prescribed temperature velocities
          tempe(:,:)        ! (2,npoin) temperature

  LOGICAL ,POINTER :: naeul(:)    !(npoin) if local system is time dependent

  REAL (kind=8), POINTER :: &
    damp(:),   & !(0, neq, or maxa) damping matrix
    accel(:),  & !(neq) accelerations
    veloc(:),  & !(neq) velocities
    velnp(:,:),& !(nrotd,npoin) velocities
    mass(:)      !(neq) equivalent DOF mass

  ! scratch
  INTEGER (kind=4), POINTER :: nodset(:),auxi(:)   ! (numnp), nodes of the eset, auxiliar array
  INTEGER (kind=4) :: numpo                        ! number of points in a set

         ! TO BE COMPLETED

CONTAINS

  SUBROUTINE dump_npo ( )
    USE ctrl_db, ONLY : ndime, ndofn, nrotd, neulr, nload, npoin, neq,maxa,ndyna, &
                        numct, bottom, top
    USE esets_db,ONLY:  rot_free,nel
    IMPLICIT NONE
    INTEGER (kind=4) :: i,j,k,n1,m1
    LOGICAL :: psi

    WRITE (50,ERR=9999) scale
    WRITE (50,ERR=9999) ((ifpre(i,k),i=1,ndofn),k=1,npoin)
    IF( rot_free )THEN
      WRITE(50,ERR=9999)( iffix(k),k=1,npoin)
    ELSE
      WRITE(50,ERR=9999) iffix(1)
    END IF

    WRITE(50,ERR=9999) (label(i),i=1,npoin)
    !WRITE(50,ERR=9999) (oldlb(i),i=1,npoio)
    WRITE(50,ERR=9999) ((coora(i,j),i=1,ndime),j=1,npoin), &
              ((coorc(i,j),i=1,ndime),j=1,npoin), &
              ((coord(i,j),i=1,ndime),j=1,npoin), &
              ((emass(i,j),i=1,ndofn),j=1,npoin)
    IF(neulr > 0)THEN
      WRITE(50,ERR=9999)(eule0(:,k),k=1,npoin)
      WRITE(50,ERR=9999)((euler(i,k),i=1,neulr),k=1,npoin),(naeul(k),k=1,npoin)
      WRITE(50,ERR=9999)((locsy(i,j),i=1,neulr),j=1,npoin)
      psi = ASSOCIATED(psia)
      WRITE(50,ERR=9999)psi
      IF( psi )THEN
        WRITE(50,ERR=9999)(psic(:,j),j=1,npoin),(psia(:,j),j=1,npoin)
      END IF
    ELSE
      WRITE(50,ERR=9999) euler(1,1)
    END IF
    IF(ASSOCIATED(cpx))WRITE(50,ERR=9999)(cpx(:,j),j=1,npoin)

    IF( ndyna > 0)THEN
      WRITE(50,ERR=9999) (accel(i),i=1,neq),(veloc(i),i=1,neq), (velnp(:,i),i=1,npoin)
      SELECT CASE (ndyna)
      CASE (1)
        n1 = maxa
        m1 = maxa
      CASE (2)
        n1 = neq
        m1 = maxa
      CASE (3)
        n1 = maxa
        m1 = 1
      CASE (4)
        n1 = neq
        m1 = 1
      END SELECT

      WRITE(50,ERR=9999) (mass(i),i=1,n1),(damp(i),i=1,m1)
    END IF
  RETURN
  9999 CALL runen2(' ')
  END SUBROUTINE dump_npo

  SUBROUTINE rest_npo ( )
    USE ctrl_db, ONLY : ndime, ndofn, nrotd, neulr, nload, npoin, neq,maxa,ndyna, &
                        numct, bottom, top
    IMPLICIT NONE
    INTEGER (kind=4) :: i,j,k,neu
    LOGICAL :: psi

    READ (51) scale
    ALLOCATE( ifpre(ndofn,npoin))
    READ (51) ((ifpre(i,k),i=1,ndofn),k=1,npoin)
    READ (51) k
    IF( k > 0 )THEN
       ALLOCATE( iffix(k))
       READ (51)( iffix(i),i=1,k)
    END IF

    ALLOCATE( label(npoin)) !, oldlb(npoio) )
    READ(51) (label(i),i=1,npoin)
    !READ(51) (oldlb(i),i=1,npoio)

    ALLOCATE( coora(ndime,npoin), coorc(ndime,npoin), coord(ndime,npoin), &
              emass(ndofn,npoin))
    READ(51)  ((coora(i,j),i=1,ndime),j=1,npoin), &
              ((coorc(i,j),i=1,ndime),j=1,npoin), &
              ((coord(i,j),i=1,ndime),j=1,npoin), &
              ((emass(i,j),i=1,ndofn),j=1,npoin)

    IF( neulr > 0)THEN
      neu = 2*ndime-3
      ALLOCATE( eule0(neu,npoin),euler(neulr,npoin), locsy(neulr,npoin), naeul(npoin) )
      READ(51) ((eule0(i,k),i=1,neu),k=1,npoin)
      READ(51) ((euler(i,k),i=1,neulr),k=1,npoin),(naeul(k),k=1,npoin)
      READ(51) ((locsy(i,j),i=1,neulr),j=1,npoin)
      READ(51) psi
      IF( psi )THEN
        IF( ndofn == 8 )THEN
          ALLOCATE(psic(2,npoin),psia(2,npoin))
        ELSE
          ALLOCATE(psic(1,npoin),psia(1,npoin))
        END IF
        READ(51)(psic(:,j),j=1,npoin),(psia(:,j),j=1,npoin)
      END IF
    ELSE
       ALLOCATE( euler (1,1))
       READ (51) euler(1,1)
    END IF
    IF(nel(19) > 0)READ(51)(cpx(:,j),j=1,npoin)

    IF( ndyna > 0)THEN
      ALLOCATE( accel(neq), veloc(neq), velnp(nrotd,npoin) )
      READ (51) (accel(i),i=1,neq),(veloc(i),i=1,neq), (velnp(:,i),i=1,npoin)
      SELECT CASE (ndyna)
      CASE (1)
        ALLOCATE( mass(maxa), damp(maxa) )
        READ (51) (mass(i),i=1,maxa),(damp(i),i=1,maxa)
      CASE (2)
        ALLOCATE( mass(neq), damp(maxa) )
        READ (51) (mass(i),i=1,neq),(damp(i),i=1,maxa)
      CASE (3)
        ALLOCATE( mass(maxa), damp(1) )
        READ (51) (mass(i),i=1,maxa),(damp(i),i=1,1)
      CASE (4)
        ALLOCATE( mass(neq), damp(1) )
        READ (51) (mass(i),i=1,neq),(damp(i),i=1,1)
      END SELECT

    END IF

  RETURN
  END SUBROUTINE rest_npo

END MODULE npo_db
