SUBROUTINE mastd4(nsegm,nnseg,lcseg,nhseg,confor,x)
!
!   Generates data base for surfaces that will act as master in some pairs.
!   This data includes the neighbours segments
!
!.... input
!....   nsegm = number of segments of the master surface considered
!....   nnseg = number of nodes per segment of the surface considered
!....   lcseg(inseg,icseg) = global node number for the local element node
!....                        [inseg] of the segment [icseg]
!....   confor : .TRUE.   conforming surface is considered
!....            .FALSE.  non-conforming surface is considered
!....   x     = nodal coordinates
!.... output
!....   nhseg(nnseg,nsegm) = neighbours segments to each side

IMPLICIT NONE
!     Dummy arguments
INTEGER (kind=4), INTENT(IN) :: nsegm,nnseg,lcseg(:,:)
INTEGER (kind=4), INTENT(OUT) :: nhseg(:,:)
LOGICAL, INTENT(IN), OPTIONAL :: confor
REAL (kind=8), INTENT(IN), OPTIONAL :: x(:,:)
!     local variables
INTEGER (kind=4) jm,km,i,j,jn,kn,icseg,inseg,ns
INTEGER (kind=4), PARAMETER :: nextn(3) = (/2,3,1/)
INTEGER (kind=4), ALLOCATABLE :: ni(:,:)
!      LOGICAL, ALLOCATABLE :: nj(:,:)

INTERFACE
  INCLUDE 'mastda.h'
END INTERFACE

!     Generate standard Data Base for conforming elements

nhseg = 0   !initializes connected segments
DO icseg = 1,nsegm                  !for each segment
  outer : DO inseg = 1,nnseg                !for each node in the segment
    IF( nhseg(inseg,icseg) /= 0)CYCLE       !if already considered cycle
    jm = lcseg(nextn(inseg),icseg)          !first node of the side
    km = lcseg(nextn(nextn(inseg)),icseg)   !second node of the side
    ! search the TWIN side
    DO i=icseg+1,nsegm                      !for the rest of the segments
      DO j=1,nnseg                          !for each side
        IF(nhseg(j,i) /= 0)CYCLE            !if already considered cycle
        jn = lcseg(nextn(j),i)              !first node of the side
        IF( jn == km )THEN                  !if second = first
          kn = lcseg(nextn(nextn(j)),i)     !second node of the side
          IF( kn == jm )THEN                !if first = second
            nhseg(inseg,icseg) = i          !this is the side
            nhseg(j,i) = icseg              !
            CYCLE outer                     !Go to next side
          END IF
        END IF
      END DO
    END DO
  END DO outer
END DO

IF( .NOT.PRESENT(confor) .OR. .NOT.PRESENT(x) )RETURN
IF( confor )RETURN

!....  computed 'best' neighbour segment for boundary sides
!....  in Non-conforming meshes

ns = COUNT(nhseg(1:3,1:nsegm) == 0) !number of non-paired sides
ALLOCATE ( ni(4,ns) )  !Auxiliar array
ni = 0                 !initializes array

! search for boundary nodes (sides), make a list
ns = 0                   !initializes position in array
DO icseg=1,nsegm         !for each segment
  DO inseg=1,nnseg          !for each side
    IF(nhseg(inseg,icseg) /= 0 )CYCLE   !if paired side cycle
    ns = ns+1                           !increase counter
    ni(1,ns) = icseg                    !assign segment number
    ni(2,ns) = inseg                    !assign side number
    i = nextn(inseg)                    !first local node of segment
    j = nextn(i)                        !second local node of segmnet
    ni(3,ns) = lcseg(i,icseg)           !first global node
    ni(4,ns) = lcseg(j,icseg)           !second global node
  END DO
END DO

! check if all sides are included in closed loops
! this is not strictly necessary but advisable to check meshes

!      ALLOCATE ( nj(2,ns) )  !Auxiliar array
!      nj = .FALSE.
!      DO i=1,ns             !for each side
!        jn = ni(3,i)        !first node of the side
!        kn = ni(4,i)        !second node of the side
!        DO j=i+1,ns         !for each remaining side
!          ! if both sides found, go to next side
!          IF( nj(1,i) .AND. nj(2,i) )EXIT
!          jm = ni(3,j)      !first node of the side
!          km = ni(4,j)      !second node of the side
!          IF( jn == km )THEN
!            ! previous side found
!            IF( nj(1,i) .OR. nj(2,j))CYCLE
!            nj(1,i) = .TRUE.
!            nj(2,j) = .TRUE.
!          ELSE IF( kn == jm )THEN
!            ! next side found
!            IF( nj(2,i) .OR. nj(1,j) )CYCLE
!            nj(2,i) = .TRUE.
!            nj(1,j) = .TRUE.
!          END IF
!        END DO
!      END DO

!       IF open lines found STOP
!      IF( ns > 0 .AND. .NOT.ALL(nj) )THEN
!        write(55,"(4i6)")ni
!        CALL runend('MASTD4: Error in SEG connectivities')
!      END IF
!      DEALLOCATE ( nj )
!     compare boundaries and assign alternate nodes

CALL mastda (ns,ni,x,nhseg)
DEALLOCATE ( ni )

RETURN
END SUBROUTINE mastd4
