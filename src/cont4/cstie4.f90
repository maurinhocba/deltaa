SUBROUTINE cstie4(rssdb,stiff,cprop,issdb,facto)

!.... compute stiffness matrix for a four-node 3-d contact element

IMPLICIT NONE
!     arguments
INTEGER (kind=4), INTENT(IN) :: issdb(:)
REAL (kind=8), INTENT(IN) :: rssdb(:),cprop(:),facto
REAL (kind=8), INTENT(OUT) :: stiff(:)
!     local variables
INTEGER (kind=4) i,j,k,ii,jj
REAL (kind=8) gap,sh(4),pnlty,vns(3),cofri,pnltc,iter
REAL (kind=8) nn(3,3),mn(3,3)

!     form shape functions
sh(1)   =  1d0
sh(3:4) = -rssdb(5:6)
sh(2)   = -1d0 - sh(3) - sh(4)
!     normal vector
vns = rssdb(2:4)

!     restore surface properties
pnlty = cprop(1) * facto
iter = DBLE(issdb(4)+1)/4d0
IF( iter < 1d0 ) pnlty = pnlty*iter

cofri = cprop(3)
IF(cofri > 0d0 .AND. issdb(3) == 1 )THEN     !If contact in previous step
  pnltc = cprop(2) * facto
  gap = ABS(rssdb(5)-rssdb(7)) + ABS(rssdb(6)-rssdb(8))
  IF(gap > .1e-11)  pnltc = pnltc*cofri
ELSE
  pnltc = 0d0
END IF

!     form nn(T)  &  1 - nn(T) operators

mn = 0d0
DO j = 1,3  ! 3 = ndime
  mn(j,j) = 1d0
  DO i = 1,3    ! 3 = ndime
    nn(i,j) = vns(i)*vns(j)    !standard tensorial product
  END DO
END DO
mn = (mn - nn)*pnltc              !orthogonal projection matrix (scaled)
nn =  nn*pnlty + mn               !tensorial product (scaled)
!     compute stiffness matrix
k = 0                          !index on stiffness matrix
DO i = 1,4                     !for each NODE
  DO ii=1,3                    !for each direction in SPACE
    DO j=i,4                   !for each node
      DO jj=1,3                !for each direction in SPACE
        IF(i == j .AND. jj < ii )CYCLE  !if under the diagonal skip
        k = k+1                !increase counter
        ! compute tangent operator for linearized kinematics
        ! (symmetric only) consistent non symmetric matrix not available
        stiff(k) = nn(ii,jj)*sh(i)*sh(j)
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE cstie4
