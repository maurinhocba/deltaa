 SUBROUTINE gidmesh(input,leng)
 !
 !        writes mesh data for GiD in ASCII form
 !
 USE data_db
 IMPLICIT NONE
 CHARACTER (len=*) :: input
 INTEGER :: leng,i

 OPEN (11,FILE=input(1:leng)//'.post.msh',STATUS='UNKNOWN',  &
       ACCESS='SEQUENTIAL',FORM='FORMATTED')
 OPEN (13,FILE=input(1:leng)//'.post.res',STATUS='UNKNOWN',  &
       ACCESS='SEQUENTIAL',FORM='FORMATTED')


 ! heading for new format
 WRITE(13,"('GiD Post Results File 1.0',/)")

 ! writes mesh at file .msh and Gauss points information at file .res
 IF( spot_sets  > 0 )  CALL pcgi01( ) ! conect DATA for element type 1 (SPOT)
 IF( truss_sets > 0 )  CALL pcgi02( ) ! conect DATA for element type 2 (TRUSS)
 IF( sol2d_sets > 0 )  CALL pcgi03( ) ! conect DATA for element type 3, 17, 20
 IF( sol3d_sets > 0 )  CALL pcgi05( ) ! conect DATA for element type 4,5,16,18,27 (SOLID)
 IF( shl3d_sets > 0 )  CALL pcgi06( ) ! conect DATA for element type 6-7(SHELL)
 IF( beame_sets > 0 )  CALL pcgi08( ) ! conect DATA for element type 8 (BEAME)
 IF( shrev_sets > 0 )  CALL pcgi09( ) ! conect DATA for element type 9 (SHREV)
 IF(   bst_sets > 0 )  CALL pcgi12( ) ! conect DATA for element type 12:15,25(BST,BSQ)
 IF( rigid_sets > 0 )  CALL pcgi10( ) ! conect DATA for element type 10(RIGID)
 !IF( drawb_sets > 0 .AND. drawb_shn > 0) CALL pcgid_db(0) ! connect DATA for drawbeads non-associated surfaces
 WRITE(11,"('End Elements')")
 CLOSE (11)

 ! writes tables for FLD maps at file .res
 IF (bst_sets > 0 .OR. sol2d_sets > 0) THEN
   WRITE(13,*)
   !Definition of table associated to 'FORMING ZONE' chart
   IF (bst_wfldFZ .OR. sol2d_wfldFZ) THEN
     WRITE(13,"('ResultRangesTable ""',a,'""')")TRIM(rrt_fr%name)
     DO i=1,rrt_fr%nval
       WRITE(13,"(f5.1,' -',f4.1,': ""',a,'""')")rrt_fr%vals(i),rrt_fr%vals(i+1),TRIM(rrt_fr%lab(i))
     END DO
     WRITE(13,"('End ResultRangesTable',/)")
     WRITE(13,*)
   END IF
   !Definition of table associated to 'SAFETY ZONE' chart
   IF (bst_wfldSZ .OR. sol2d_wfldSZ) THEN
     WRITE(13,"('ResultRangesTable ""',a,'""')")TRIM(rrt_sf%name)
     DO i=1,rrt_sf%nval
       WRITE(13,"(f5.1,' -',f4.1,': ""',a,'""')")rrt_sf%vals(i),rrt_sf%vals(i+1),TRIM(rrt_sf%lab(i))
     END DO
     WRITE(13,"('End ResultRangesTable',/)")
   END IF
 END IF
 IF (drawb_sets > 0 .AND. drawb_wafz) THEN
   WRITE(13,*)
   !Definition of table associated to 'DRAWBEAD AFFECTED ZONE' chart
   WRITE(13,"('ResultRangesTable ""',a,'""')")TRIM(rrt_db%name)
   DO i=1,rrt_db%nval
     WRITE(13,"(f5.1,' -',f4.1,': ""',a,'""')")rrt_db%vals(i),rrt_db%vals(i+1),TRIM(rrt_db%lab(i))
   END DO
   WRITE(13,"('End ResultRangesTable',/)")
 END IF

 RETURN
 END SUBROUTINE gidmesh
