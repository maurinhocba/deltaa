 SUBROUTINE elemnt (TASK, name, dtime, ttime, istop, flag1, flag2, &
                    sumat, iload, iforce, mass, stiff, ustif, resid)
 !********************************************************************
 !
 !*** standard ELEMENT routine
 !
 !********************************************************************

 IMPLICIT NONE
 !                        Dummy (compulsory) arguments
 CHARACTER (len=*) :: TASK
 !                        Dummy (optional) arguments
 CHARACTER(len=*),OPTIONAL:: name
 INTEGER (kind=4), OPTIONAL :: istop,iload,iforce
 REAL (kind=8), OPTIONAL :: dtime,ttime,sumat,resid(:,:),stiff(:),ustif(:),mass(:)
 LOGICAL, OPTIONAL :: flag1,flag2

 END SUBROUTINE elemnt
