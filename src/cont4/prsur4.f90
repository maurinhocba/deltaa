      SUBROUTINE prsur4(iwric,bottom,top)

!.... write rigid surface data for visualization
!.... open file for contact forces between surfaces
!.... update flags BOTTOM and TOP

      USE c_input, ONLY : openfi
      USE cont4_db
      IMPLICIT NONE
!     dummy arguments
      INTEGER (kind=4), INTENT(IN) :: iwric
      LOGICAL, INTENT(OUT) :: bottom,top
!     local variables
      INTEGER (kind=4) :: i,j,isurf,nsegm

      TYPE (surf4_db), POINTER :: surf
      TYPE (pair4_db), POINTER :: pair
      LOGICAL :: openp

!.... Open file por post-process
      IF(iwric > 0)THEN           !if contact forces are desired
        CALL openfi(41)           !open file for contact forces
        WRITE(41)npair,3          !number of forces and vector size
        pair => headp             !point to first pair
        DO i=1,npair              !for each pair
          WRITE(41) i,pair%pname  !contact pair order and name
          pair => pair%next       !point to next pair
        END DO
      END IF

!.... WRITE contact surface for post-proccessing
      CALL openfi(18)               !open file
      surf => shead                 !pointer to first surface
      DO isurf=1,nsurf              !for each surface
        nsegm  = surf%nsegm         !number of segments defining the surface
        IF( surf%iwrit )THEN        !if visualization desired
          WRITE(18) nsegm,surf%sname  ! No of segments
          WRITE(18) 3             ! No of nodes per segment
          WRITE(18) surf%lcseg(1:3,1:nsegm) !connectivities
        END IF
        surf => surf%next       !pointer to next surface
      END DO

      CLOSE (18)                    !close file

!.... Update global flags BOTTOM and TOP

      openp = .FALSE.          !for surface pressure
      pair => headp             !point to first pair
      DO i=1,npair              !for each pair
        IF( pair%mtsur == -2 .OR. pair%slsur == -2 )bottom = .TRUE.
        IF( pair%mtsur ==  2 .OR. pair%slsur ==  2 )top    = .TRUE.
        IF( pair%press .OR. pair%wrink )THEN    !if surface pressure or wrinkles
          IF( .NOT.openp )THEN
            CALL openfi(44)    !file for surface data
            openp = .TRUE.
          END IF
          surf => shead                 !pointer to first surface
          isurf = 1
          DO
            IF(isurf == pair%islav)EXIT
            isurf = isurf+1
            surf => surf%next
          END DO

          WRITE(44)surf%sname,pair%press,pair%wrink !name, press and wrinkles
          WRITE(44)surf%ncnod,surf%nsegm            !No of nodes & segments in the surface
          WRITE(44)(surf%lcnod(j),j=1,surf%ncnod)   !internal numeration of nodes
          WRITE(44)(surf%lcseg(:,j),j=1,surf%nsegm) !surface connectivities

        END IF
        pair => pair%next       !point to next pair
      END DO
      IF( openp ) WRITE(44)'NONAME',.FALSE.,.FALSE. !mark to end

      IF( wear )THEN    !if friction work
        CALL openfi(45)
        surf => shead                 !pointer to first surface

        DO j=1,nsurf
          IF(surf%nsegm > 0)THEN
            WRITE(45)surf%nsegm            !No of segments in the surface
            WRITE(45)(surf%lcseg(:,i),i=1,surf%nsegm) !surface connectivities
          END IF
          surf => surf%next
        END DO
        WRITE(45) 0  !mark to end, associated to nsegm
      END IF
      RETURN
      END SUBROUTINE prsur4
