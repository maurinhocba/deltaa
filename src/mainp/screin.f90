 SUBROUTINE screin(task,cpui)

 USE name_db,ONLY : prognm
 REAL(kind=8) :: cpui
 CHARACTER(len=*),INTENT(IN) :: task !not used

 INTEGER (kind=4) i,lng
 CHARACTER(len=8) :: version = '7.0.0   ', verdate = '01/06/20'

 INTERFACE
   INCLUDE 'timuse.h'
 END INTERFACE

 CALL timuse(cpui)

 !--------------------write out headers
!IF(TRIM(prognm) == 'ALPHA')THEN
!  WRITE(6,"(/                                             &
!           & t22,'  оо   оо    ооооо  оо  оо   оо  '/     &
!           & t22,' оооо  оо    оо  оо оо  оо  оооо '/     &
!           & t22,'оо  оо оо    оо  оо оо  оо оо  оо'/     &
!           & t22,'оо  оо оо    оо  оо оо  оо оо  оо'/     &
!           & t22,'оооооо оо    ооооо  оооооо оооооо'/     &
!           & t22,'оо  оо оо    оо     оо  оо оо  оо'/     &
!           & t22,'оо  оо ооооо оо     оо  оо оо  оо'/)")
 IF(TRIM(prognm) == 'DELTA')THEN
   WRITE(6,"(/                                             &
            & t22,'оооо   оооооо оо    оооооо   оо  '/     &
            & t22,'оо оо  оо     оо      оо    оооо '/     &
            & t22,'оо  оо оо     оо      оо   оо  оо'/     &
            & t22,'оо  оо оооо   оо      оо   оо  оо'/     &
            & t22,'оо  оо оо     оо      оо   оооооо'/     &
            & t22,'оо оо  оо     оо      оо   оо  оо'/     &
            & t22,'оооо   оооооо ооооо   оо   оо  оо'/)")
 ELSE
 !                    PROGRAM name
    lng = LEN_TRIM(prognm)
    WRITE(6,"(///////33X,<lng>(A,1X),/)") (prognm(i:i),i=1,lng)
 END IF
 !                    PROGRAM version number and date
 WRITE (6,"(33X,'version:',A9/33X,'date   :',A9,/)")version,verdate
 !      WRITE (17)prognm,version,verdate
 WRITE (6,"(t20,'STATIC/DYNAMIC ANALYSIS OF STRUCTURES',//                &
           &t22,'     BY THE FINITE ELEMENT METHOD',//)")
 WRITE(6,"(                                                               &
      & t16,'GROUP OF NUMERICAL METHODS IN STRUCTURAL MECHANICS',/,       &
      & t16,'             DEPARTMENT OF STRUCTURES',//,                   &
      & t16,'          NATIONAL UNIVERSITY OF CORDOBA',//,                &
      & t16,'             (C) G M N M E   CopyRight',//,                  &
      & t16,'                CORDOBA - ARGENTINA')")
 WRITE (6,"('+ WORKING ON PRELIMINARY CALCULATIONS')")

 RETURN
 END SUBROUTINE screin

 SUBROUTINE screen(istep,itera,dlamb,lambd)
 USE outp_db, ONLY : cpui
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: istep,itera
 REAL (kind=8),INTENT(IN) ::   dlamb,lambd

 REAL (kind=8)    ::  cpuc
 INTEGER (kind=4) :: seg,iha,imi,isg
 INTERFACE
   INCLUDE 'timuse.h'
 END INTERFACE

 CALL timuse(cpuc)
 cpuc = cpuc - cpui
 seg = INT(cpuc)
 isg = MOD(seg,60)
 imi = MOD(seg/60,60)
 iha = seg/3600
 WRITE(6,1) istep,itera,lambd,dlamb,iha,imi,isg

 RETURN
 1 FORMAT('+',' STEP=',I5,' ITER=',I3,'  LAMBDA=',E11.5,             &
         &    '   DLAMB=',E11.4,'     CPU ',I2,':',I2,':',I2)

 END SUBROUTINE screen

 SUBROUTINE screed(istep,itera,ttime,tdisp)
 USE outp_db, ONLY : cpui
 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN) :: istep,itera
 REAL (kind=8),INTENT(IN)  :: ttime,tdisp

 REAL (kind=8)    :: cpuc
 INTEGER (kind=4) :: seg,iha,imi,isg
 INTERFACE
   INCLUDE 'timuse.h'
 END INTERFACE

 CALL timuse(cpuc)
 cpuc = cpuc - cpui
 seg = INT(cpuc)
 isg = MOD(seg,60)
 imi = MOD(seg/60,60)
 iha = seg/3600
 WRITE(6,1) istep,itera,ttime,tdisp,iha,imi,isg

 RETURN
 1 FORMAT('+',' STEP=',I5,' ITER=',I3,'    TTIME=',E11.4,            &
        & '   TDISP =',E11.5,'    CPU ',I2,':',I2,':',I2)

 END SUBROUTINE screed

 SUBROUTINE scree2(istep,itera,dlamb,lambd,du,ddu)
 USE outp_db, ONLY : cpui
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: istep,itera
 REAL (kind=8),INTENT(IN) ::   dlamb,lambd,du,ddu

 REAL (kind=8)    ::  cpuc
 INTEGER (kind=4) :: seg,iha,imi,isg
 INTERFACE
   INCLUDE 'timuse.h'
 END INTERFACE

 CALL timuse(cpuc)
 cpuc = cpuc - cpui
 seg = INT(cpuc)
 isg = MOD(seg,60)
 imi = MOD(seg/60,60)
 iha = seg/3600
 WRITE(6,1) istep,itera,lambd,dlamb,du,ddu,iha,imi,isg

 RETURN
 1 FORMAT(I5,I3,' L=',E11.4,' DL=',E11.4,             &
         & ' U=',E11.3,' DDU=',E11.3,' ',I2,':',I2,':',I2)

 END SUBROUTINE scree2
SUBROUTINE flushf(n)
INTEGER (kind = 4) n
!     dummy routine (in main-frame fortran it clears buffers for file n)
!      CALL flush(n)           !Not ALL FORTRANs
END SUBROUTINE flushf
