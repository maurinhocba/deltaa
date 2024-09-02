 SUBROUTINE actual (displ,lambd,dlamb)
 !********************************************************************
 !
 !***   actualization of current values into converged values
 !
 !********************************************************************
 !USE cont_db, ONLY : ncont
 IMPLICIT NONE

 REAL (kind=8),INTENT(IN) :: displ(:),lambd,dlamb

 INTERFACE
   INCLUDE 'elemnt.h'
   !INCLUDE 'contac.h'
 END INTERFACE

 !***   update internal variables for elements

 CALL elemnt('ACTUAL')

 !IF(ncont > 0) CALL contac('ACTUAL',ncont,ttime=lambd,dtime=dlamb,veloc=displ)

 RETURN
 END SUBROUTINE actual
