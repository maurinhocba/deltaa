SUBROUTINE cdump4(npoin)
!
!     writes contact data for restar
!
USE cont4_db   !INTENT(IN)
IMPLICIT NONE
INTEGER (kind=4), INTENT(IN) :: npoin
INTEGER (kind=4) :: nn(41,3)     !variable for statistic records
COMMON /proj3/ nn

!local variables
TYPE (pair4_db), POINTER :: pair
TYPE (surf4_db), POINTER :: surf
INTEGER (kind=4) :: ipair,isurf,i,j


WRITE(50) nsurf,npair,oldis,disma,ffdis,ctime,wear   !Control variables

IF( wear )WRITE(50) (wwear(i),i=1,npoin)

!     Dumps pair data base

pair => headp
DO  ipair=1,npair
  WRITE(50) pair%pname , pair%master, pair%slave ,              &
            pair%imast , pair%islav , pair%indcon, pair%ncnod , &
            pair%mtsur , pair%slsur , pair%freq ,  pair%bhforc, &
            pair%npenal, pair%tpenal, pair%static, pair%kinet,  &
            pair%cutoff, pair%gapinc, pair%start , pair%end,    &
            pair%prev,   pair%press , pair%wrink

  WRITE(50) ((pair%issdb(i,j),i=1,nisdb),j=1,pair%ncnod)
  WRITE(50) ((pair%rssdb(i,j),i=1,nrsdb),j=1,pair%ncnod)
  IF(pair%press) WRITE(50) (pair%presn(j),j=1,pair%ncnod)
  IF(pair%wrink) WRITE(50) (pair%mingp(j),j=1,pair%ncnod)

  pair => pair%next
END DO

!     Dumps surface data base

surf => shead
DO isurf=1,nsurf
  WRITE(50) surf%sname, surf%cxc,    surf%bottom, surf%confor,  &
            surf%ncnod, surf%curved, surf%nsegm,  surf%iwrit,   &
            surf%imcod, surf%iscod,  surf%press,  surf%density
  IF( surf%iscod )WRITE(50) (surf%lcnod(i),i=1,surf%ncnod)
  IF( surf%nsegm > 0 )WRITE(50) ((surf%lcseg(i,j),i=1,3),j=1,surf%nsegm)
  IF( surf%imcod )THEN
    WRITE(50) ((surf%nhseg(i,j),i=1,3),j=1,surf%nsegm)
    WRITE(50) ((surf%xc(i,j),i=1,3),j=1,surf%nsegm)
  END IF
  IF(surf%bottom )THEN
    WRITE(50) ((surf%lcseb(i,j),i=1,3),j=1,surf%nsegm)
    WRITE(50) ((surf%nhseb(i,j),i=1,3),j=1,surf%nsegm)
  END IF
  IF(surf%curved .AND. surf%iscod ) &
                     WRITE(50) ((surf%tn(i,j),i=1,3),j=1,surf%ncnod)
  IF(surf%curved .AND. surf%imcod ) &
                     WRITE(50) ((surf%cu(i,j),i=1,3),j=1,surf%nsegm)

  surf => surf%next
END DO

WRITE(50)((nn(i,j),i=1,41),j=1,3)       !statistical records

RETURN

END SUBROUTINE cdump4
