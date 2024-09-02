 SUBROUTINE conend ( )
 !********************************************************************
 !
 !*** contact definition card READ
 !
 !********************************************************************
 USE c_input

 CALL listen('CONEND')
 IF (.NOT.exists('ENDCON')) &
   CALL runend('CONEND: END_CONTACT_DEFINITION CARD')

 RETURN
 END SUBROUTINE conend
