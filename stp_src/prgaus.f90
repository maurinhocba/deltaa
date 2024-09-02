 SUBROUTINE prgaus (input,leng)
 !
 !  Process Gaussian and Nodal variables (smoothed from Gaussian)
 !
 USE data_db
 IMPLICIT NONE
 CHARACTER (len=25) :: input
 INTEGER :: leng

 INTEGER (kind=4) ::  etype, iset, &
                      spot_iset, truss_iset, sol2d_iset, sol3d_iset, shl3d_iset, &
                      beame_iset, shrev_iset, bst_iset
 LOGICAL :: quads

 spot_iset  = 0       !initializes set counters
 truss_iset = 0
 sol2d_iset = 0
 sol3d_iset = 0
 shl3d_iset = 0
 beame_iset = 0
 shrev_iset = 0
 bst_iset   = 0
 !spher_iset = 0

 ! first loop to read Gaussian Variables

 DO iset=1,nsets

   etype = setyp(iset)     !element type

   SELECT CASE(etype)

   CASE (1)      ! Process Gauss data for TRUSS elements
     CALL gauss1 (spot_iset)

   CASE (2)      ! Process Gauss data for TRUSS elements
     CALL gauss2 (truss_iset)

   CASE (17,19,20)  ! Process Gauss data for 2-D SOLID elements
     CALL gauss3(sol2d_iset)

   CASE (4,5,12,16,18,27,29)      ! Process Gauss data for 3-D SOLID elements
     CALL gauss5 (sol3d_iset)

   CASE (6,7)    ! Process Gauss data for 3-D SHELL (Simo Theory) elements
     CALL gauss6 (shl3d_iset)

   CASE (8)      ! Process Gauss data for 3-D BEAM (Simo Theory) element
     CALL gauss8 (shrev_iset)

   CASE (9,11)   ! Process Gauss data for 2-D SHELL/BEAM element
     IF( ntype /= 4 )THEN
       CALL gauss9 (shrev_iset,etype)
     ELSE
       CALL gauss9c(shrev_iset)
     END IF

   CASE (13:15,25)  ! Process Gauss data for SHELL (BST) element
     CALL gaus12 (bst_iset,etype)

   END SELECT
 END DO

 ! smoothing (averaging)

 IF(truss_nps > 0 .AND. truss_nvarn > 0 )  &
   CALL averag(truss_nps,truss_nvarn,truss_vargs,truss_accpn,quads)
 IF(sol2d_nps > 0 .AND. sol2d_nvarn > 0 )THEN
   CALL averag(sol2d_nps,sol2d_nvarn,sol2d_vargs,sol2d_accpn,quads)
   IF( quads )CALL aves2d()
 END IF
 IF(sol3d_nps > 0 .AND. sol3d_nvarn > 0 )  THEN
   CALL averag(sol3d_nps,sol3d_nvarn,sol3d_vargs,sol3d_accpn,quads)
   IF( quads )CALL aves3d()
 END IF
 IF(shl3d_nps > 0 .AND. shl3d_nvarn > 0 )THEN
   CALL averag(shl3d_nps,shl3d_nvarn,shl3d_vargs,shl3d_accpn,quads)
   IF( quads )CALL aveshl()
 END IF
 IF(beame_nps > 0 .AND. beame_nvarn > 0 )  &
   CALL averag(beame_nps,beame_nvarn,beame_vargs,beame_accpn,quads)
 IF(shrev_nps > 0 .AND. shrev_nvarn > 0 .AND. ntype /= 4 ) &
   CALL averag(shrev_nps,shrev_nvarn,shrev_vargs,shrev_accpn,quads)
 IF(bst_nps   > 0 .AND. bst_nvarn   > 0 )  &
   CALL averag(bst_nps,  bst_nvarn,  bst_vargs,  bst_accpn  ,quads)
 ! Process special Gauss data after smoothing
 CALL spgaus()
 CALL readbr()
 ! print results

 SELECT CASE (ip)
 CASE (1)
   !CALL tecbinres(input,leng)
 CASE (2)
   CALL gidres ( etype )
 CASE (3)
   CALL tecres(input,leng)
 CASE (4,5)
   CALL gidresbin( etype )
 END SELECT

 RETURN
 END SUBROUTINE prgaus
