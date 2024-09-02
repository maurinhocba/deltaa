 SUBROUTINE tecres (input,leng)
 !
 !  Print result file for TecPlot 360
 !  presently for ASCII option in BLOCK FORMAT
 !
 !  UNDER DEVELOPMENT: Hints for a new format
 !  1-Generate a GRID file with original coordinates and connectivities
 !    separate different element types: Each element type must have its own numeration
 !    but different sets of same element type share numeration (for smoothing)
 !    3 node lineseg are divided in two
 !    use BLOCK format instead of POINT format
 !  2-Generate SOLUTION files with:
 !    a)global nodal variables: present coordinates
 !                              displacements, velocities, acceleration
 !    b)specific element type variables, use PASSIVEVARLIST not to include unused variables
 !    Gauss point variables are be averaged and passed as cell-centered
 !
 !    Tasks missing: binary files
 !                   user defined variables

 USE data_db
 IMPLICIT NONE
 ! dummy arguments
 CHARACTER (len=*) :: input !file root name
 INTEGER :: leng            !length of root name
 ! local variables
 CHARACTER(len=6 ) :: name2    !composed file name for TECPLOT
 CHARACTER(len=30) :: sname    !set name
 CHARACTER(len=3 ) :: ch,comp  !function to convert a 3 digit integer into a string
 CHARACTER (len=55) :: line    !string
 INTEGER :: l                  !length of root name
 INTEGER(kind=4), SAVE ::  iload=0  !counts number of printed steps
 LOGICAL, PARAMETER :: flag = .TRUE.      !solution

 CALL tecvar()            !generate string of variables

 iload = iload + 1        !results for step ILOAD

 name2 = 'l'//ch(iload)     !to be add to root name
 ! open file
 OPEN(11,FILE=input(1:leng)//name2(1:4)//'.dat',STATUS='UNKNOWN', &
         ACCESS='SEQUENTIAL',FORM='FORMATTED')
 ! writes header
 line = '"'//TRIM(text)//' At time: '
 l = LEN_TRIM(line) + 2
 WRITE(line(l:l+13),"(e13.5,a)")ttime,'"'
 WRITE(11,"('title = ',a)")TRIM(line)
 WRITE(11,"('FILETYPE = SOLUTION')")
 WRITE(11,"('VARIABLES = ',a)")TRIM(t_vars(1:t_p))

 !IF( spot_nps > 0 )CALL pcte01 (flag,iload)    !  1
  IF(truss_nps > 0 )CALL pcte02 (flag,iload)    !  2
  IF(sol2d_nps > 0 )CALL pcte03 (flag,iload)    !  3
  IF(sol3d_nps > 0 )CALL pcte05 (flag,iload)    !  4
  IF(shl3d_nps > 0 )CALL pcte06 (flag,iload)    !  5
  IF(beame_nps > 0 )CALL pcte08 (flag,iload)    !  6
  IF(shrev_nps > 0 )CALL pcte09 (flag,iload)    !  7
  IF(rigid_nps > 0 )CALL pcte10 (flag,iload)    !  8
  IF(bst_nps   > 0 )CALL pcte12 (flag,iload)    !  9
 !IF(?????_nps > 0 )CALL pcte?? (flag,iload)    ! 10

 CLOSE(11)
 RETURN

 END SUBROUTINE tecres
