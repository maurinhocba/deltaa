 SUBROUTINE wrtpos(flag)
 !
 ! write initial data for post-process
 ! updates stage coordinates for a new strategy
 !
 USE ctrl_db, ONLY: ndime,npoin,ndofn,nrotd,numct,text,ndyna,ndoft,itemp
 USE c_input, ONLY: openfi
 USE npo_db, ONLY: label,coord,coors,cpx
 USE outp_db
 USE mat_dba, ONLY : nusect
 IMPLICIT NONE

   ! Dummy variables
   INTEGER(kind=4):: flag
   ! Local variables
   INTEGER (kind=4) :: i , iflag, nada1 = 1, nada0 = 0
   INTEGER (kind=4), PARAMETER :: cero = 0
   REAL(kind=8) :: xx(ndime,2)

   INTERFACE
     INCLUDE 'contac.h'
     INCLUDE 'elemnt.h'
   END INTERFACE

   IF (flag /= 1) THEN  !Open files and write initial data for post-processing
     CALL openfi(16)    ! open post-process files
     CALL openfi(17)
     ! open file for especific output (for HISTORY or CURVES)
     IF (nreqd > 0) THEN              !displacement
       CALL openfi(11) ! ,oflag=iflag)
       !IF (iflag == 0) THEN
         WRITE(11,ERR=9999) postype,nreqd,ndofn,(label(nprqd(i)),i=1,nreqd)
       !END IF
     END IF
     IF (nreql > 0) THEN              !equivalent nodal forces
       CALL openfi(12) !,oflag=iflag)
       !IF (iflag == 0) THEN
         WRITE(12,ERR=9999) postype,nreql,ndofn,(label(nprql(i)),i=1,nreql)
       !END IF
     END IF
     IF (iener == 1) THEN             !energies
       CALL openfi(15) !,oflag=iflag)
       !IF (iflag == 0) THEN
         WRITE(15,ERR=9999) postype,11   !number of values to print per step
       !END IF
     END IF
     IF (nreqv > 0) THEN              !velocities
       CALL openfi(19) !,oflag=iflag)
       !IF (iflag == 0) THEN
         WRITE(19,ERR=9999) postype,nreqv,nrotd,(label(nprqv(i)),i=1,nreqv)
       !END IF
     END IF
     IF (nreqa > 0) THEN              !accelerations
       CALL openfi(20) !,oflag=iflag)
       !IF (iflag == 0) THEN
         WRITE(20,ERR=9999) postype,nreqa,nrotd,(label(nprqa(i)),i=1,nreqa)
       !END IF
     END IF
     !IF (nreqt > 0) THEN              !temperatures
     !  CALL openfi(70,oflag=iflag)
     !  IF (iflag == 0) THEN
     !    WRITE(70,ERR=9999)nreqt,ndoft,(label(nprqt(i)),i=1,nreqt)
     !  END IF
     !END IF
      !writes information for post-processing associated to follower loads
     CALL wrtfl0()
     CALL actcd()  !updates stage coordinates
     ! write some control parameter
     i = 0
     IF( itemp ) i = ndoft
     WRITE(17,ERR=9999) ndime,ndofn,i !i = number of termal DOFs
     WRITE(17,ERR=9999) npoin,MIN(ndyna,1),text
     IF( .NOT.ASSOCIATED (cpx) )THEN
       DO i=1,npoin
         WRITE(17,ERR=9999) label(i),coord(1:ndime,i),coors(1:ndime,i)
       END DO
     ELSE
       DO i=1,npoin
         IF( cpx(1,i) == 0 )THEN
           WRITE(17,ERR=9999) label(i),coord(1:ndime,i),coors(1:ndime,i)
         ELSE ! (+/-1)
           xx(:,1) = (2d0*coord(:,i) + coord(:,cpx(2,i))+coord(:,cpx(3,i)))/4d0
           xx(:,2) = (2d0*coors(:,i) + coors(:,cpx(2,i))+coors(:,cpx(3,i)))/4d0
           WRITE(17,ERR=9999) label(i),xx(1:ndime,1),xx(1:ndime,2)
         END IF
       END DO
     END IF
     ! additional data for sections (flc no implemented)
     WRITE(17,ERR=9999) nusect
     WRITE(17,ERR=9999) (nada0,i=1,nusect),nada0     !lblfld, nfc
     !
     CALL elemnt ('WRTPOS') !write element data (mesh connectivities)
     WRITE(17,ERR=9999) nada0,nada0,'                              '       !to indicate end of elements sets

     IF (numct > 0) THEN                 !contact
       IF (nreqc > 0) THEN               !contact forces at selected points
         CALL openfi(14) !,oflag=iflag)
         !IF (iflag == 0) THEN
           WRITE(14,ERR=9999)postype, nreqc,ndime,(label(nprqc(i)),i=1,nreqc)
         !END IF
       END IF
       ! writing rigid (contact) surfaces for postprocess
       CALL contac ('WRTSUR')
     END IF

   ELSE !IF (flag == 1) THEN   !Open files to add data for post-processing
 ! ***
     CALL openfi(1016)      ! open post-process files
     CALL openfi(1017)

     IF (nreqd > 0) CALL openfi(11)
     IF (nreql > 0) CALL openfi(12)
     IF (iener == 1) CALL openfi(15)
     IF (nreqv > 0) CALL openfi(19)
     IF (nreqa > 0) CALL openfi(20)
     IF (numct > 0) THEN
       IF (nreqc > 0) CALL openfi(14)
       CALL contac ('WRTSUR')
     END IF
     !CALL wrtfl0 !(52)
 !        IF( thickc )CALL openfi(46,oflag=iflag)
 !        IF(nreqp > 0) CALL openfi(150,oflag=iflag)
 !        IF(nreqt > 0) CALL openfi(151,oflag=iflag)

   END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE wrtpos
