SUBROUTINE bgraph( bgout, nfile )

  !  read list of connected (contacting) pairs of spheres
  !  write contact bond graph file

  USE data_db
  USE bonds_db
  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: bgout    ! bond graph output file name
  INTEGER,          INTENT(IN) :: nfile    ! number of bond file

  CHARACTER (len=1), SAVE :: ni(0:9)=(/'0','1','2','3','4','5','6','7','8','9'/)
  CHARACTER (len=13) :: elmt
  INTEGER :: i,j,j0,j1,j2,leng,nelem,nn(2)
  REAL (kind=8) :: f(2)
  TYPE (list), POINTER :: ls        ! list of 2-point elements

  CALL ls_ini (ls)                  !initializes list
  DO
    ! read list of connected (contacting) pairs of spheres
    READ(16) nn, f
    IF (nn(1) == 0) EXIT
    CALL add_at_end (ls, nn, f)     !add node to SET
  END DO

  leng = LEN_TRIM(bgout)            !length of the file name
  j0 = MOD(nfile,1000)/100          !digits for hundreds, tens & units
  j1 = MOD(nfile,100)/10
  j2 = MOD(nfile,10)
  ! open file for mesh definition
  OPEN (99,FILE=bgout(1:leng)//'_bg'//ni(j0)//ni(j1)//ni(j2)// &
        '_f.msh',STATUS='UNKNOWN',FORM='FORMATTED')
  elmt = 'Linear       '
  WRITE(99,"('MESH    dimension =',i2,' ElemType ',a13, &
           & ' Nnode = ',i2,/,'Coordinates')")ndime,elmt,2
  DO i=1,npoin
    WRITE(99,"(i8,3e18.10)")label(i),coorf(1:ndime,i)
  END DO
  WRITE(99,"('End coordinates'/'Elements')")
  nelem = get_length (ls)           !get number of elements in the list
  CALL ls_head (ls)                 !go to top of the list of elems
  DO i=1,nelem                      !for each node in the list
    nn = get_nodes (ls)             !get nodal data
    WRITE(99,'(15i8)') i,(label(nn(j)),j=1,2),1
    CALL ls_next (ls)               !next node in the list
  END DO
  WRITE(99,"('End Elements')")
  CLOSE (99)

  ! open file for results (contact forces)
  OPEN (99,FILE=bgout(1:leng)//'_bg'//ni(j0)//ni(j1)//ni(j2)// &
        '_f.res',STATUS='UNKNOWN',FORM='FORMATTED')
  WRITE(99,2005)
  WRITE(99,2004)ttime
  WRITE(99,"('FORCE        ')")
  CALL ls_head (ls)                 !go to top of the list of elems
  DO i=1,nelem                      !for each element in the list
    f = get_f (ls)                  !get list element data (force)
    WRITE(99,'(i8,e18.10)')i,f(1)
    CALL ls_next (ls)               !next element in the list
  END DO
  CALL dalloc_list (ls)             !deallocate the list
  CLOSE (99)

  RETURN

 2004 FORMAT('TRUSS_FORCES    1 ',e12.4,'   1   2   1 "Chain"')
 2005 FORMAT('GaussPoints "Chain" ElemType Linear'/                    &
     &       '  Number of Gauss Points: 1',/ '  Nodes not Included ',/ &
     &       '  Natural Coordinates: Internal',/ 'End Gausspoints')

END SUBROUTINE bgraph
