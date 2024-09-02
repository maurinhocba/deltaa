 SUBROUTINE check_orientation8(lnods,x,t)
 USE npo_db, ONLY : label
 ! check element connectivities
 IMPLICIT NONE
 !dummy arguments
 INTEGER(kind=4), INTENT(IN OUT) :: lnods(8) !connectivities
 REAL(kind=8), INTENT(IN OUT) :: x(3,8)      !coordinates
 REAL(kind=8), INTENT(IN) :: t(3)            !shell normal
 !local variables
 REAL(kind=8) :: t1(3),t2(3),t3(3),ca,z(3),y(3,8)
 INTEGER(kind=4) :: ln(8)
 INTERFACE
   INCLUDE 'vecuni.h'
   INCLUDE 'vecpro.h'
 END INTERFACE

 z = t                ! normal shell direction
 CALL vecuni(3,z,ca)  !unit shell normal
 t1 =  -x(:,1)+x(:,2)+x(:,3)-x(:,4)   !x_xita
 t2 =  -x(:,1)-x(:,2)+x(:,3)+x(:,4)   !x_eta
 CALL vecpro(t1,t2,t3)                !x_xita x x_eta
 CALL vecuni(3,t3,ca)                 !element normal
 ca = ABS( DOT_PRODUCT(z,t3) )        !COS between normals
 IF( ca > 0.866 )RETURN !O.K   angle is less than 30 degrees
 ! explore other options
 ln = lnods                           !assign connectivities to auxiliar array
 y = x                                !assign coordinates to auxiliar array
 ! set face 1 4 8 5 as bottom face
 t1 =  -x(:,1)+x(:,4)+x(:,8)-x(:,5)   !x_xita
 t2 =  -x(:,1)-x(:,4)+x(:,5)+x(:,8)   !x_eta
 CALL vecpro(t1,t2,t3)                !x_xita x x_eta
 CALL vecuni(3,t3,ca)                 !element normal
 ca = ABS( DOT_PRODUCT(z,t3) )        !COS between shell normal and first arista
 IF( ca > 0.866 )THEN                   !if angle is less than 30 degrees
   lnods = ln((/ 1,4,8,5, 2,3,7,6 /))   !swap connectivities
   x     = y(:,(/ 1,4,8,5, 2,3,7,6 /))  !swap coordinates
   WRITE(58,"(8i6)")label(lnods)
   RETURN !O.K
 END IF
 ! set face 1 5 6 2 as bottom face
 t1 =  -x(:,1)+x(:,2)+x(:,6)-x(:,5)   !x_xita
 t2 =  -x(:,1)-x(:,2)+x(:,5)+x(:,6)   !x_eta
 CALL vecpro(t1,t2,t3)                !x_xita x x_eta
 CALL vecuni(3,t3,ca)                 !element normal
 ca = ABS( DOT_PRODUCT(z,t3) )        !COS between shell normal and first arista
 IF( ca > 0.866 )THEN                   !if angle is less than 30 degrees
   lnods = ln((/ 1,5,6,2, 4,8,7,3 /))   !swap connectivities
   x     = y(:,(/ 1,5,6,2, 4,8,7,3 /))  !swap coordinates
   WRITE(58,"(8i6)")label(lnods)
   RETURN !O.K
 END IF
 WRITE(55,"(8i6,'WARNING: ORIENTATION COULD NOT BE CORRECTED')")lnods

 RETURN
 END SUBROUTINE check_orientation8
