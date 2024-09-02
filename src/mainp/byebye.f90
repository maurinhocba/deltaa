      SUBROUTINE byebye (time)
!*********************************************************************
!     writes cpu times
!*********************************************************************
      USE c_input, ONLY : openfi
      USE gvar_db, ONLY : actchk
      IMPLICIT NONE
      REAL (kind=8),INTENT(IN OUT) :: time(:)

      INTEGER (kind=4) :: i,iha,imi,isg

      IF(actchk)RETURN

      time(10) = time(10)-SUM(time(6:9))
      DO i=2,11 !% times
        time(20+i) = time(i)*100d0/time(1)
      END DO

      CALL openfi(2)
      WRITE (2,40,ERR=9999) (time(i),time(i+20),i=2,11)
      iha = time(1)/3600.d0                      !hours
      imi =(time(1)-iha*3600.0)/60.0             !minutes
      isg = time(1)-iha*3600.0-imi*60.0          !seconds
      WRITE (2,42,ERR=9999) iha,imi,isg
      CLOSE(2,STATUS='keep')
      RETURN
   40 FORMAT(///,'T I M I N G   I N F O R M A T I O N',//,&
     &15X,'ACCION                           CPU(SEC)         %',//,&
     &15X,'DATA INPUT ...................',E14.4,2X,F6.2,' %',//,&     !1
     &15X,'RESTART INPUT-OUTPUT .........',E14.4,2X,F6.2,' %',//,&     !2
     &15X,'INITIAL COMPUTATIONS..........',E14.4,2X,F6.2,' %',//,&     !3
     &15X,'MASS MATRIX COMPUTATION ......',E14.4,2X,F6.2,' %',//,&     !4
     &15X,'STIFFNESS MATRIX COMPUTATION .',E14.4,2X,F6.2,' %',//,&     !5
     &15X,'RESIDUAL FORCES COMPUTATION ..',E14.4,2X,F6.2,' %',//,&     !6
     &15X,'CONTACT COMPUTATIONS .........',E14.4,2X,F6.2,' %',//,&     !7
     &15X,'LDL DECOMPOSITION + SOLVING ..',E14.4,2X,F6.2,' %',//,&     !8
     &15X,'ITERATIVE PROCESS ............',E14.4,2X,F6.2,' %',//,&     !9
     &15X,'RESULTS OUTPUT ...............',E14.4,2X,F6.2,' %',//)      !10
   41 FORMAT(5X,'T O T A L                    ',E14.4,//)
   42 FORMAT(////,15X,'HOURS: ',I6,8X,'MINUTES:',I2,8X,'SECONDS:',I2,//)
 9999 CALL runen2('')
      END SUBROUTINE byebye
