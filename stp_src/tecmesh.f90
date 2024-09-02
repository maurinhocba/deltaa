 SUBROUTINE tecmesh(input,leng)
 !
 !    Print original mesh in ASCII format
 !
 USE data_db
 IMPLICIT NONE
 ! dummy arguments
 CHARACTER (len=*) :: input  !file root name
 INTEGER :: leng             !length of file root name
 ! local variables
 CHARACTER (len=22) :: line    !string
 CHARACTER (len=14) :: linev   !string
 LOGICAL, PARAMETER :: flag = .FALSE.      !mesh
 INTEGER(kind=4), PARAMETER :: iload = 0   !mesh


 ! open file
 OPEN(11,FILE=input(1:leng)//'l000'//'.dat',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='FORMATTED')
 ! writes header
 line ='"'//TRIM(text)//'"'  !  header
 linev = '"X" "Y"' !  coordinates names
 IF( ndime == 3 )  linev(8:11) = ' "Z"'
 WRITE(11,"('title = ',a)")TRIM(line)
 WRITE(11,"('FILETYPE=GRID')")
 WRITE(11,"('VARIABLES = ',a)")TRIM(linev)

 sh_var ='VARSHARELIST=([1,2])'
 IF( ndime == 3 )sh_var(19:22)= ',3])'

  IF(truss_nps > 0 )CALL pcte02 (flag,iload)
  IF(sol2d_nps > 0 )CALL pcte03 (flag,iload)
  IF(sol3d_nps > 0 )CALL pcte05 (flag,iload)
  IF(shl3d_nps > 0 )CALL pcte06 (flag,iload)
  IF(beame_nps > 0 )CALL pcte08 (flag,iload)
  IF(shrev_nps > 0 )CALL pcte09 (flag,iload)
  IF(rigid_nps > 0 )CALL pcte10 (flag,iload)
  IF(bst_nps   > 0 )CALL pcte12 (flag,iload)

 RETURN
 END SUBROUTINE tecmesh
