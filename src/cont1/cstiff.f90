 SUBROUTINE cstiff(nsymm,ttime,force,stiff,ustif)
 !.... contact stiffness routine
 USE cont1_db
 IMPLICIT NONE
 !     arguments
 INTEGER (kind=4), INTENT(IN) :: nsymm
 REAL (kind=8), INTENT(IN OUT) :: force(*),stiff(*),ustif(*)
 REAL (kind=8), INTENT(IN) :: ttime
 !     local variables
 LOGICAL :: prev
 INTEGER (kind=4) :: ipair,i,iss,ims,ien(3),icnod,nearn,iter
 INTEGER (kind=4), POINTER :: lcnod(:),lcseg(:,:)
 REAL (kind=8) :: xfact,cprop(4),lstif(21),astif(15)
 TYPE (pair1_db), POINTER :: p
 TYPE (surf1_db), POINTER :: master,slave

 INTERFACE
   INCLUDE 'cstie1.h'
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

     cprop = (/ p%npenal, p%tpenal, p%static, p%kinet /)

     DO icnod = 1, p%ncnod
       iter  = p%issdb(4,icnod)             !actual active contact node
       prev  = p%rssdb(6,icnod) < 0d0       !previous active contact node

       IF( iter > 0 .OR. prev ) THEN
         !....     active contact element    nen : projected segment
         nearn = p%issdb(1,icnod)
         ien = (/ lcnod(icnod), lcseg(:,nearn) /)
         CALL cfrmlm(ien, 3 ,ctime,xfact,p%rssdb(4,icnod),p%indcon)
         ! compute element left-hand-side
         ! 3-node (s,m1,m2) 2D contact element (node-segment penetration)
         CALL cstie1(p%rssdb(:,icnod),lstif,cprop,p%issdb(:,icnod),xfact,  &
                     nsymm,astif)
         !  form lm array  issdb(2-3): isn, imn1,imn2 & add element left-hand-side
         CALL addlhs(ien,2,3,force(1),stiff(1),ustif(1),nsymm,lstif(1),  &
                     astif(1),6)         !   2=ndime, 3=nen, 6=ndime*nen
       END IF
     END DO

   END IF
   p => p%next
 END DO

 RETURN
 END SUBROUTINE cstiff
