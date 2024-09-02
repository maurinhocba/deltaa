 SUBROUTINE readbr ( )
  ! read drawbead affected zones for all drawbeads
  USE data_db

  REAL(kind=8) :: ftime          !present time (not used)
  TYPE (drawb), POINTER :: eset  !pointer to a drawbead
  INTEGER(kind=4) :: iset        !counter

  IF( drawb_sets > 0 )THEN      !check if drawbeads exist
    READ(43)  ftime              !read present time (could be compared to check)
    eset => drawb_head             !point to first drawbead
    DO iset=1,drawb_sets             !loop all the drawbeads
      READ(43)  eset%daz               !read affected zone variable
      eset => eset%next                !point to next drawbead
    END DO
  END IF

  RETURN
  END SUBROUTINE readbr
