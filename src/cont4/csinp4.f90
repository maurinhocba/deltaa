SUBROUTINE csinp4(maxve,iwrit,label,npoin,coord,oldsr,oldpr, &
                  codes,surname,tsurf,surdi)

!.... input and generate contact surfaces data

USE param_db, ONLY : mich
USE c_input
USE cont4_db
IMPLICIT NONE
!     arguments
INTEGER (kind=4), INTENT(IN) :: iwrit,npoin,label(:),oldpr,maxve,surdi
INTEGER (kind=4), INTENT(IN OUT) :: oldsr,tsurf,codes(2,surdi)
REAL (kind=8), INTENT(IN) :: coord(:,:)
CHARACTER (len=6), INTENT (IN OUT) :: surname(surdi)

!     local variables
LOGICAL :: found
INTEGER (kind=4) :: j,k,isurf,deflt,ipair
TYPE (pair4_db), POINTER :: pair
TYPE (surf4_db), POINTER :: surf, anter, posic
CHARACTER (len=6) :: sname
CHARACTER(len=mich):: inttoch

INTERFACE
  INCLUDE 'csdat4.h'
END INTERFACE

!.... Check if Old Surfaces will be used by New Pairs
 pair => headp                  !point to first pair
 DO ipair=1,oldpr               !skip OLD pair
   pair => pair%next
 END DO
 DO ipair=oldpr+1,npair         !for each new-pair
   DO j=1,2                     !for each surface
     IF( j == 1 )THEN           ! master surface (j=1)
       sname = pair%master                          !master surface label
     ELSE                       ! slave surface  (j=2)
       sname = pair%slave                           !slave surface label
     END IF
     k = 0                      !initializes surface position in list
     surf => shead              !point to first (old) surface
     DO
       k=k+1                    !increase counter
       IF( k > oldsr)EXIT       !all existing surfaces considered EXIT
       IF( surf%sname == sname )THEN  !if surface found
         IF( j == 1)THEN      !master surface
           IF( surf%imcod )THEN     !surface was supposed master
             pair%imast = k            !O.K., assign internal position
             ! if BOTTOM surface and this was not previously assigned
             IF(pair%mtsur <  0 .AND. .NOT.ASSOCIATED(surf%lcseb))THEN
               ALLOCATE(surf%lcseb(3,surf%nsegm),surf%nhseb(3,surf%nsegm))
               surf%lcseb(1,:) = surf%lcseg(1,:) !reversed nodes order
               surf%lcseb(2,:) = surf%lcseg(3,:)
               surf%lcseb(3,:) = surf%lcseg(2,:)
               surf%nhseb(1,:) = surf%nhseg(1,:) !reversed segment neighbour
               surf%nhseb(2,:) = surf%nhseg(3,:)
               surf%nhseb(3,:) = surf%nhseg(2,:)
               surf%bottom = .TRUE.              !set to .TRUE.
             END IF
             EXIT  !exit surface search
           ELSE
             WRITE(lures,"('Surface',6a,' Cannot be MASTER')")sname
             CALL runend('CSINPT: Surface Cannot be MASTER   ')
           END IF
         ELSE                 !slave surface
           IF( surf%iscod )THEN      !surface was supposed slave
             pair%islav = k             !O.K., assign internal position
             pair%ncnod = surf%ncnod    !O.K.
             ALLOCATE ( pair%issdb(nisdb,surf%ncnod) ) !slave Int. DB
             ALLOCATE ( pair%rssdb(nrsdb,surf%ncnod) ) !slave Real DB
             pair%issdb = 0                            !initializes
             pair%rssdb = 0d0                          !initializes
             IF( pair%press )THEN                      !compute normal force
               ALLOCATE( pair%presn(surf%ncnod) )      !reserve space
               pair%presn = 0d0                        !initializes
               surf%press = .TRUE.                     !set to true
             END IF
             IF( pair%wrink )THEN                      !wrinkle control
               ALLOCATE( pair%mingp(surf%ncnod) )      !reserve space
               pair%mingp = 1d9                        !initializes
             END IF
             EXIT  !exit surface search
           ELSE
             WRITE(lures,"('Surface',a6,' Cannot be SLAVE')")sname
             CALL runend('CSINPT: Surface Cannot be SLAVE    ')
           END IF
         END IF
       END IF
       surf => surf%next    !point to next surface
     END DO  !search in surface Data Base
   END DO ! j=1,2
   pair => pair%next        !consider next pair
 END DO !ipair=oldpr+1,npair

 !.... loop over new contact surfaces (READ)

 DO
   ! loop over surfaces

   CALL listen('CSINPT')
   IF (.NOT.exists('ISURF ',k)) THEN
     backs = .TRUE.
     EXIT
   END IF

   CALL new_surf4 (surf)              !allocate new surface
   IF( param(k) == 0d0)THEN           !if no associated number
     surf%sname = words(k+1)(1:6)     !next label is the name
   ELSE                               !name is a two digit number
     surf%sname = '    '//TRIM(inttoch(INT(param(k)),2))  !fill with blanks
   END IF
   WRITE(lures,"(/' CONTACT SURFACE:',a8,/)")surf%sname  !print

   !check if surface label already used
   CALL srch4_srf (shead, anter, posic, surf%sname, found)
   IF (found) CALL runend ('RDSURF: Surface already defined !!!')

   isurf= 0        !search for surface label in CODES array
   DO
     isurf = isurf+1
     IF( isurf > tsurf .AND. isurf < surdi)THEN   !not defined but possible
       surname(isurf) = surf%sname
       WRITE(lures,"(' WARNING: surface ',6a,' not defined in any ', &
              & 'pair')")surf%sname
     ELSE IF(surname(isurf) == surf%sname)THEN    !found
       EXIT
     ELSE IF(isurf == surdi)THEN        !space exhausted
       CALL runend('CSINPT: More surfaces than expected')
     END IF
   END DO

   CALL listen('CSINPT')
   !                          read flag WRITE
   deflt = 0
   IF(codes(1,isurf) > 0 .AND. codes(2,isurf) == 0)deflt = 1
   surf%iwrit = getint('WRITE ',deflt,' Save surface for Post-Processing..') == 1
   IF( iwrit > 0)THEN
     IF( surf%iwrit )THEN
       WRITE(lures,"(12x,'WILL be saved for post-Processing ')")
     ELSE
       WRITE(lures,"(12x,'WILL NOT be saved for post-Processing')")
     END IF
   END IF

   !                          read flag NCNOD
   deflt = 0
   IF(codes(1,isurf) == 0 .AND. codes(2,isurf) /= 0)deflt = 1
   surf%ncnod =getint('NCNOD ',deflt,' Nodes in the surface will be read ')
   IF(iwrit > 0)THEN
     IF( surf%ncnod > 0)THEN
       WRITE(lures,"(12x,'List of node will be read ')")
     ELSE
       WRITE(lures,"(12x,'No list of nodes expected')")
     END IF
   END IF

   !                          read flag NSEGM
   deflt = 0
   IF( surf%iwrit ) deflt = 1
   IF( codes(1,isurf) > 0 .OR. surf%ncnod == 0) deflt = 1
   surf%nsegm =getint('NSEGM ',deflt,' Segments of surface will be read..')
   IF(iwrit > 0)THEN
     IF( surf%nsegm > 0)THEN
       WRITE(lures,"(12x,'Connectivities will be read ')")
     ELSE
       WRITE(lures,"(12x,'No connectivities expected')")
     END IF
   END IF

   !                          read flag IMCOD
   deflt = 0
   IF( codes(1,isurf) > 0 ) deflt = 1
   surf%imcod = getint('IMCOD ',deflt,' Code for master surface ..........') == 1
   IF( surf%imcod .AND. surf%nsegm == 0 )CALL runend &
      ('CSINPT:IMCOD=1 requires nsegm > 0  ')
   IF(iwrit > 0)THEN
     IF( surf%imcod )THEN
       WRITE(lures,"(12x,'Will act as MASTER for some pair')")
     ELSE
       WRITE(lures,"(12x,'Will not act as MASTER for any pair')")
     END IF
   END IF
   IF( .NOT.surf%imcod .AND. deflt == 1 )THEN
     WRITE(lures,"(' Bad Master Code for surface',a7)") surf%sname
     CALL runend('CSINP4: Master Code Must be 1 !!!   ')
   ELSE IF(surf%imcod .AND. deflt == 0 )THEN
     WRITE(lures,"(' Warning Master Code for surface',a7,' set= 1', &
     & /,t15,'but is not as master in any pair')")surf%sname
   END IF

   !                          read flag ISCOD
   deflt = 0
   IF( codes(2,isurf) > 0 ) deflt = 1
   surf%iscod =getint('ISCOD ',deflt,' Code for slave surface ...........') == 1
   IF(iwrit > 0)THEN
     IF( surf%iscod )THEN
       WRITE(lures,"(12x,'Will act as SLAVE for some pair')")
     ELSE
       WRITE(lures,"(12x,'Will not act as SLAVE for any pair')")
     END IF
   END IF

   IF( .NOT.surf%iscod .AND. deflt == 1 )THEN
     WRITE(lures,"(' Bad SLAVE Code for surface',a7)")surf%sname
     CALL runend('CSINPT: SLAVE Code Must be 1 !!!    ')
   ELSE IF( surf%iscod .AND. deflt == 0 )THEN
     WRITE(lures,"(' Warning SLAVE Code for surface',a7,' set= 1', &
     & /,t15,'but is not as SLAVE in any pair')")surf%sname
   END IF

   ! warning if WRITE = 1 and NSEGM = 0
   IF( surf%iwrit .AND. surf%nsegm == 0)THEN
     WRITE(lures,"(' Sorry, NSEGM must be > 0 to save surface for ' &
                 &  'visualization')")
     surf%iwrit = .FALSE.
   END IF

   surf%density =getrea('DENSIT', 0d0 ,' Surface density for mass comput...')

   surf%curved = exists('CURVED')
   IF( surf%curved ) WRITE(lures,"(12x,'Surface will be treated as CURVED')")

   surf%confor = exists('CONFOR')
   IF( surf%confor ) WRITE(lures,"(12x,'Surface will be supposed CONFORMING')")

   ! read surface connectivities and generate data base
   CALL csdat4(maxve,iwrit,label,npoin,coord,surf)

   nsurf = 1 + nsurf    !increase number of surfaces

   !Check if the surface will be used as SLAVE for some pairs

   pair => headp               !point to first pair
   DO ipair=1,codes(2,isurf)   !codes(2) = times used as slave
     DO
       IF( pair%slave == surname(isurf) )EXIT    !found
       pair => pair%next                         !point to next pair
     END DO
     pair%islav = nsurf               !assign internal position
     pair%ncnod = surf%ncnod          !assign number of nodes
     ALLOCATE( pair%issdb(nisdb,surf%ncnod),  &   !reserve space for
               pair%rssdb(nrsdb,surf%ncnod) )     !slave data base
     pair%issdb = 0                   !initializes Int. array
     pair%rssdb = 0d0                 !initializes Real array
     IF( pair%press )THEN
       ALLOCATE( pair%presn(pair%ncnod) )   !reserve space
       surf%press = .TRUE.                  !set to true
       pair%presn = 0d0                !initialize
     END IF
     IF( pair%wrink )THEN                      !wrinkle control
       ALLOCATE( pair%mingp(pair%ncnod) )      !reserve space
       pair%mingp = 1d9                        !initializes
     END IF
     pair => pair%next                !next pair
   END DO

   !Check if the surface will be used as MASTER for some pairs

   pair => headp               !point to first pair
   DO ipair=1,codes(1,isurf)   !codes(1) = times used as master
     DO
       IF( pair%master == surname(isurf) )EXIT   !found
       pair => pair%next                         !point to next pair
     END DO
     pair%imast = nsurf               !assign internal position

     ! if BOTTOM surface and this was not previously assigned
     IF(pair%mtsur < 0 .AND. .NOT.ASSOCIATED(surf%lcseb))THEN
       ALLOCATE(surf%lcseb(3,surf%nsegm),surf%nhseb(3,surf%nsegm))
       surf%lcseb(1,:) = surf%lcseg(1,:) !reversed nodes order
       surf%lcseb(2,:) = surf%lcseg(3,:)
       surf%lcseb(3,:) = surf%lcseg(2,:)
       surf%nhseb(1,:) = surf%nhseg(1,:) !reversed segment neighbour
       surf%nhseb(2,:) = surf%nhseg(3,:)
       surf%nhseb(3,:) = surf%nhseg(2,:)
       surf%bottom = .TRUE.  !set to .TRUE.
     END IF
     pair => pair%next                   !point to next pair
   END DO
   CALL add4_srf (surf, shead, stail)    !add surface to data_base

 END DO         !surface reading

 ! check if all the pairs have surface correctly defined
 pair => headp               !point to first pair
 found = .FALSE.
 DO ipair=1,npair            !for each pair
   IF( pair%imast < 1 .OR. pair%imast > nsurf )THEN
     WRITE(lures,"(' master surface', a7,' asssociated to pair',a7,  &
&                  ' does NOT exists')")pair%master,pair%pname
     WRITE(55,"(' master surface', a7,' asssociated to pair',a7,     &
&                  ' does NOT exists')")pair%master,pair%pname
     found = .TRUE.
   END IF
   IF( pair%islav < 1 .OR. pair%islav > nsurf )THEN
     WRITE(lures,"(' slave surface ', a7,' asssociated to pair',a7,  &
&                  ' does NOT exists')")pair%slave,pair%pname
     WRITE(55,"(' slave surface ', a7,' asssociated to pair',a7,     &
&                  ' does NOT exists')")pair%slave,pair%pname
     found = .TRUE.
   END IF

   pair => pair%next                   !point to next pair
 END DO
 IF( found )CALL runend('CSINPT: contact surface missing   ')
 RETURN
 END SUBROUTINE csinp4
