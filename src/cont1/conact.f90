SUBROUTINE conact(ttime)

!.... updates internal friction variables into last converged

USE cont1_db
IMPLICIT NONE
!     arguments
REAL (kind=8), INTENT(IN) :: ttime

!     local variables
INTEGER (kind=4) :: ia,icnod,ipair,nearn,np
REAL (kind=4) :: bhf,gap
INTEGER (kind=4), POINTER :: issdb(:,:)
REAL (kind=8), POINTER :: rssdb(:,:)
TYPE (pair1_db), POINTER :: p

IF(oldis < disma)THEN
  ffdis = oldis/disma
  IF(ffdis > 0.975) ffdis = 1d0
ELSE
  ffdis = 1d0
END IF
oldis = disma

p => headp
DO ipair = 1, npair
  IF( p%start <= ttime .AND. p%end >= ttime)THEN
    issdb => p%issdb
    rssdb => p%rssdb
    DO icnod=1,p%ncnod
      !       issdb(1,icnod) = nearest segment
      nearn = issdb(1,icnod)
      gap = rssdb(1,icnod)
      IF( nearn > 0 .AND. gap < 0d0 )THEN   !if proyected and penetrated
        issdb(2,icnod) = issdb(1,icnod)     !target segment
        issdb(3,icnod) =  1   !penetration is true
        IF(p%static > 0d0)THEN
          IF( rssdb(5,icnod) /= rssdb(4,icnod) )THEN
            rssdb(5,icnod) = rssdb(4,icnod)    !local coordinates of onset
            issdb(3,icnod) =  2                !sliding
          END IF
        END IF
        rssdb(6,icnod) = rssdb(1,icnod)     !keep penetration of step
      ELSE
        IF( nearn > -1000000 .AND. gap > disma )THEN !searched and too far
          issdb(1,icnod) = -ABS(nearn) - 10000000  !do not search for 10 steps
        ELSE IF( nearn < -1000000) THEN !not searched
          nearn = -nearn
          np = (nearn / 1000000)*ffdis !remaining steps for next search
          ia = MAX(INT(np-1),0)        !update remaining steps
          nearn = MOD(nearn, 1000000) + ia * 1000000  !new value for nearn
          issdb(1,icnod) = -nearn    !store new value
        END IF
        issdb(3,icnod) = MIN(issdb(3,icnod),1) - 1
        rssdb(6,icnod) = 0d0              ! no effective penetration
      END IF
      issdb(4,icnod) = 0                  ! initializes iterative penetration

      IF( p%wrink )THEN
        IF( rssdb(1,icnod) < p%mingp(icnod) )THEN
          IF( rssdb(1,icnod) > 0d0 ) p%mingp(icnod) = rssdb(1,icnod)
        END IF
      END IF

    END DO
    ! Modifies penalty to account for constant B_H_forces
    IF( p%bhforc > 0d0 )THEN
      bhf = ABS(surtf(3,ipair))                    !total force
      p%npenal = p%npenal * p%bhforc /bhf !modifies penalty
    END IF
  END IF

  p => p%next
END DO

RETURN
END SUBROUTINE conact
