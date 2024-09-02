 SUBROUTINE pcgi01 ( )
 !
 !  Write mesh information for GiD for SPOT elements
 !
 USE data_db
 IMPLICIT NONE
 TYPE( spot ), POINTER :: e
 INTEGER :: iset,iel,n,k
 LOGICAL :: gauss
 CHARACTER (len=13) :: etype
 CHARACTER (len=32) :: sname
 CHARACTER (len=34) :: gpname

 INTERFACE
   SUBROUTINE prtcoo_spot(nelem,nnode,iset,sname,x,fl,lnods,matno)
     INTEGER (kind=4), INTENT(IN) :: nelem,nnode,fl,iset,matno
     CHARACTER (len=32), INTENT(IN) :: sname
     INTEGER(kind=4), POINTER :: lnods(:,:)
     REAL (kind=8), POINTER :: x(:,:,:)
   END SUBROUTINE prtcoo_spot
 END INTERFACE

 gauss = spot_nvarg > 0         !Gauss variables to be print

 e => spot_head                 !point to first element set

 DO iset=1,spot_sets            !for each element set
   k = LEN_TRIM(e%sname)
   sname = '"'//e%sname(1:k)//'"'
   gpname = '"GP'//sname(2:32)
   k = k + 4
   CALL prtcoo_spot(e%nelem,e%nnode,iset,sname,e%a_x,e%first_l,e%lnods,e%matno)        !write coordinates and headers
   IF( gauss )THEN
     etype = 'Linear       '
     WRITE(13,"('GaussPoints  ',a,' Elemtype  ',a13,3x,a)") &
                gpname(1:k),etype,sname
     WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
             &  '  Natural Coordinates: internal',/,     &
             &  'End GaussPoints')")e%ngaus
   END IF
   e => e%next
 END DO


 RETURN

 END SUBROUTINE pcgi01
