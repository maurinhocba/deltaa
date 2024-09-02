 SUBROUTINE gidbinmesh(input,leng)
 !
 !  writes mesh data for GiD in Binary form
 !
 USE data_db
 IMPLICIT NONE
 CHARACTER(len=*):: input
 INTEGER:: leng,i

 !  open GID POST RESULT FILE
 IF( ip == 5 )THEN !ASCII type using GiD routines
   CALL GID_OPENPOSTMESHFILE(input(1:leng)//'.post.msh',0)
   CALL GID_OPENPOSTRESULTFILE(input(1:leng)//'.post.res',0)
 ELSE
   CALL GID_OPENPOSTRESULTFILE(input(1:leng)//'.post.bin',2)
 END IF

 ! writes mesh at GiD binary file for post
 IF (spot_sets  > 0) CALL pcmgib01() !connect DATA for element type 1 (SPOT)
 IF (truss_sets > 0) CALL pcmgib02() !connect DATA for element type 2 (TRUSS)
 IF (sol2d_sets > 0) CALL pcmgib03() !connect DATA for element type 1, 3 or 4
 IF (sol3d_sets > 0) CALL pcmgib05() !connect DATA for element type 5 (SOLID)
 IF (  bst_sets > 0) CALL pcmgib12() !connect DATA for element type 12-14(BST)
 IF (shl3d_sets > 0) CALL pcmgib06() !connect DATA for element type 6-7(SHELL)
 IF (beame_sets > 0) CALL pcmgib08() !connect DATA for element type 8 (BEAME)
 IF (shrev_sets > 0) CALL pcmgib09() !connect DATA for element type 9 (SHREV)
 IF (rigid_sets > 0) CALL pcmgib10() !connect DATA for element type 10(RIGID)
 !IF( drawb_sets > 0 .AND. drawb_shn > 0) CALL pcgid_db(1) ! connect DATA for drawbeads non-associated surfaces

 ! writes Gauss points information at GiD binary file for post
 IF (spot_sets  > 0) CALL pcggib01() !GP DATA for element type 1 (SPOT)
 IF (truss_sets > 0) CALL pcggib02() !GP DATA for element type 2 (TRUSS)
 IF (sol2d_sets > 0) CALL pcggib03() !GP DATA for element type 1, 3 or 4
 IF (sol3d_sets > 0) CALL pcggib05() !GP DATA for element type 5 (SOLID)
 IF (  bst_sets > 0) CALL pcggib12() !GP DATA for element type 12-14(BST)
 IF (shl3d_sets > 0) CALL pcggib06() !GP DATA for element type 6-7(SHELL)
 IF (beame_sets > 0) CALL pcggib08() !GP DATA for element type 8 (BEAME)
 IF (shrev_sets > 0) CALL pcggib09() !GP DATA for element type 9 (SHREV)
 !IF( drawb_sets > 0 .AND. drawb_shn > 0) CALL pcgid_db(2) ! GP DATA for drawbeads non-associated surfaces

 ! writes tables for FLD maps
 IF (bst_sets > 0 .OR. sol2d_sets > 0) THEN
   !Definition of table associated to 'FORMING ZONE' chart
   IF (bst_wfldFZ .OR. sol2d_wfldFZ) THEN
     CALL GID_BEGINRANGETABLE(TRIM(rrt_fr%name))
     DO i=1,rrt_fr%nval
       CALL GID_WRITERANGE(rrt_fr%vals(i),rrt_fr%vals(i+1),TRIM(rrt_fr%lab(i)))
     END DO
     CALL GID_ENDRANGETABLE()
   END IF
   !Definition of table associated to 'SAFETY ZONE' chart
   IF (bst_wfldSZ .OR. sol2d_wfldSZ) THEN
     CALL GID_BEGINRANGETABLE(TRIM(rrt_sf%name))
     DO i=1,rrt_sf%nval
       CALL GID_WRITERANGE(rrt_sf%vals(i),rrt_sf%vals(i+1),TRIM(rrt_sf%lab(i)))
     END DO
     CALL GID_ENDRANGETABLE()
   END IF
 END IF
 !Definition of table associated to 'DRAWBEAD AFFECTED ZONE' chart
 IF (drawb_sets > 0 .AND. drawb_wafz) THEN
     CALL GID_BEGINRANGETABLE(TRIM(rrt_db%name))
     DO i=1,rrt_db%nval
       CALL GID_WRITERANGE(rrt_db%vals(i),rrt_db%vals(i+1),TRIM(rrt_db%lab(i)))
     END DO
     CALL GID_ENDRANGETABLE()
   END IF

 RETURN
 END SUBROUTINE gidbinmesh
