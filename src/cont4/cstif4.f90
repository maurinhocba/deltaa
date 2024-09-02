SUBROUTINE cstif4(ttime,force,stiff)
!.... contact stiffness routine
USE cont4_db
IMPLICIT NONE
!     arguments
REAL (kind=8), INTENT(IN OUT) :: force(*),stiff(*)
REAL (kind=8), INTENT(IN) :: ttime
!     local variables
LOGICAL :: prev
INTEGER (kind=4) ::ipair,i,iss,ims,ien(4),icnod,nearn,iter
INTEGER (kind=4), POINTER :: lcnod(:),lcseg(:,:)
REAL (kind=8) :: xfact,dummy,cprop(3),lstif(78)
TYPE (pair4_db), POINTER :: p
TYPE (surf4_db), POINTER :: master,slave

INTERFACE
  INCLUDE 'cstie4.h'
  INCLUDE 'cfrmlm.h'
END INTERFACE

p => headp

DO ipair=1,npair

  IF( p%start <= ttime .AND. p%end >= ttime) THEN
    !.... identify master and slave surfaces numbers
    ims = p%imast       !master surface order in list
    master => shead     !point to first surface
    DO i=1,ims-1        !loop until pointer is correct
      master => master%next
    END DO
    IF( p%mtsur < 0)THEN      ! if master surface uses the bottom surface
      lcseg => master%lcseb
    ELSE                      ! else uses top surface
      lcseg => master%lcseg
    END IF

    iss = p%islav         !slave surface order in list
    slave => shead        !point to first surface
    DO i=1,iss-1          !loop until pointer is correct
      slave => slave%next
    END DO
    lcnod => slave%lcnod

    cprop = (/ p%npenal, p%tpenal, p%static /)

    DO icnod = 1, p%ncnod
      iter  = p%issdb(4,icnod)             !actual active contact node
      prev  = p%rssdb(9,icnod) < 0d0       !previous active contact node

      IF( iter > 0 .OR. prev ) THEN
        !....     active contact element    nen : projected segment
        nearn = p%issdb(1,icnod)
        ien = (/ lcnod(icnod), lcseg(:,nearn) /)
        CALL cfrmlm(ien, 4 ,ctime,xfact,p%rssdb(5,icnod),p%indcon)
        ! compute element left-hand-side
        ! 4-node (s,m1,m2,m3) 3D contact element (node-triangle penetration)
        CALL cstie4(p%rssdb(:,icnod),lstif,cprop,p%issdb(:,icnod),xfact)
        !  form lm array  issdb(2-4): isn, imn1, imn2 & add element left-hand-side
        CALL addlhs(ien,3,4,force(1),stiff(1),dummy,0,lstif(1),dummy,12)
        !         3=ndime, 4=nen, dummy=ustif 0=nsymm, dummy=astif, 12=ndime*nen
      END IF
    END DO

  END IF
  p => p%next
END DO

RETURN
END SUBROUTINE cstif4
