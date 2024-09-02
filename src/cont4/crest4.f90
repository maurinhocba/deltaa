SUBROUTINE crest4( npoin )
!
!     Reads contact data from a restar file
!
USE cont4_db  !INTENT(OUT)
IMPLICIT NONE
INTEGER (kind=4), INTENT(IN) :: npoin
INTEGER (kind=4) :: nn(41,3)     !variable for statistic records
COMMON /proj3/ nn

!Local variables
TYPE (pair4_db), POINTER :: pair
TYPE (surf4_db), POINTER :: surf
INTEGER (kind=4) :: ipair,isurf,i,j


READ(51) nsurf,npair,oldis,disma,ffdis,ctime,wear    !control variables

IF( wear )THEN
  ALLOCATE( wwear(npoin) )
  READ(51) (wwear(i),i=1,npoin)
END IF

!...  Read Pair Data Base

CALL ini_cont4(headp,tailp)              !initializes pair list

DO  ipair=1,npair
  CALL new_pair4 (pair)
  READ(51) pair%pname , pair%master, pair%slave ,              &
           pair%imast , pair%islav , pair%indcon, pair%ncnod , &
           pair%mtsur , pair%slsur , pair%freq ,  pair%bhforc, &
           pair%npenal, pair%tpenal, pair%static, pair%kinet,  &
           pair%cutoff, pair%gapinc, pair%start , pair%end,    &
           pair%prev,   pair%press , pair%wrink
  ALLOCATE ( pair%issdb(nisdb,pair%ncnod) )
  ALLOCATE ( pair%rssdb(nrsdb,pair%ncnod) )
  READ(51) ((pair%issdb(i,j),i=1,nisdb),j=1,pair%ncnod)
  READ(51) ((pair%rssdb(i,j),i=1,nrsdb),j=1,pair%ncnod)
  IF(pair%press) THEN
    ALLOCATE( pair%presn(pair%ncnod) )
    READ(51) (pair%presn(j),j=1,pair%ncnod)
  END IF
  IF(pair%wrink) THEN
    ALLOCATE( pair%mingp(pair%ncnod) )
    READ(51) (pair%mingp(j),j=1,pair%ncnod)
  END IF
  CALL add_pair4 (pair, headp, tailp)
END DO

!...  Read Surface Database

CALL ini_srf4 (shead,stail)              !initializes surface list

DO isurf=1,nsurf
  CALL new_surf4 (surf)
  READ(51)  surf%sname, surf%cxc,    surf%bottom, surf%confor,  &
            surf%ncnod, surf%curved, surf%nsegm,  surf%iwrit,   &
            surf%imcod, surf%iscod,  surf%press,  surf%density
  IF(surf%iscod)THEN
    ALLOCATE ( surf%lcnod(surf%ncnod) )
    READ(51) (surf%lcnod(i),i=1,surf%ncnod)
  END IF
  IF(surf%nsegm > 0)THEN
    ALLOCATE ( surf%lcseg(3,surf%nsegm) )
    READ(51) ((surf%lcseg(i,j),i=1,3),j=1,surf%nsegm)
  END IF
  IF(surf%imcod)THEN
    ALLOCATE ( surf%nhseg(3,surf%nsegm) , surf%xc(3,surf%nsegm) )
    READ(51) ((surf%nhseg(i,j),i=1,3),j=1,surf%nsegm)
    READ(51) ((surf%xc(i,j),i=1,3),j=1,surf%nsegm)
  END IF
  IF(surf%bottom)THEN
    ALLOCATE ( surf%lcseb(3,surf%nsegm) , surf%nhseb(3,surf%nsegm) )
    READ(51) ((surf%lcseb(i,j),i=1,3),j=1,surf%nsegm)
    READ(51) ((surf%nhseb(i,j),i=1,3),j=1,surf%nsegm)
  END IF
  IF(surf%curved)THEN
    IF(surf%iscod == 1)THEN
      ALLOCATE( surf%tn(3,surf%ncnod) )
      READ(51) ((surf%tn(i,j),i=1,3),j=1,surf%ncnod)
    END IF
    IF(surf%imcod)THEN
      ALLOCATE( surf%cu(3,surf%nsegm) )
      READ(51) ((surf%cu(i,j),i=1,3),j=1,surf%nsegm)
    END IF
  END IF

  CALL add4_srf (surf, shead, stail)

END DO

READ (51)((nn(i,j),i=1,41),j=1,3)  !read statistics

RETURN

END SUBROUTINE crest4
