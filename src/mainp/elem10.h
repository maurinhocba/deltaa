 SUBROUTINE elem10(TASK, elsnam, ttime, iload, iforce, &
                   coora, locsy, sumat, mass, emass,   &
                   istop,flag1,flag2,iwrit)

 ! Tasks Routine for element RIGID

 IMPLICIT NONE
 CHARACTER(len=*),INTENT(IN):: TASK

 ! optional parameters

 CHARACTER (len=*),OPTIONAL :: elsnam
 LOGICAL, OPTIONAL :: flag1,flag2
 INTEGER (kind=4), OPTIONAL  :: istop, iload, iforce, iwrit
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime,coora(:,:),locsy(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),emass(:,:)

 END SUBROUTINE elem10
