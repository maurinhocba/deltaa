 SUBROUTINE prgi05( )
 !
 !  print results for GiD for
 !  3-D solid elements
 !
 USE data_db
 IMPLICIT NONE

 ! local variables
 TYPE (sol3d), POINTER :: eset
 TYPE (udmat), POINTER :: postd
 INTEGER :: i,j,k,iv,is,nelem,nnode,ielem,ngaus,g,ng,nvar,naux,matno,dim,etype,npg
 INTEGER, POINTER :: nodes(:,:)
 REAL(kind=8), POINTER :: vargs(:,:)
 CHARACTER (len=34) :: gpname = 'no name yet '
 LOGICAL :: direct
 !     non standart GiD Gauss point positions (INTERNAL)
 INTEGER, PARAMETER :: o(9,6) = RESHAPE( (/ 1, 0, 0, 0, 0, 0, 0, 0, 0,   & !1 hexa-tetra-prism
                                            1, 2, 0, 0, 0, 0, 0, 0, 0,   &  !2 prism
                                            1, 2, 3, 4, 0, 0, 0, 0, 0,   &  !3 hexa
                                            1, 1, 1, 2, 2, 2, 0, 0, 0,   &  !4 prism 2 ==> 6
                                            1, 2, 3, 4, 5, 6, 7, 8, 9,   &  !5 prism normal
                                            1, 2, 4, 3, 5, 6, 8, 7, 0 /),&  !6 hexa
                                            (/ 9,6 /) )
! or this                                   1, 5, 7, 3, 2, 6, 8, 4/), &  !hexa
 REAL(kind=8) :: press
 REAL(kind=8) :: auxil(6,8)
 CHARACTER (len=27) :: f_3dten = "(i8,6e12.4,  (/,8x,6e12.4))" , &
                       f_3dvec = "(i8,3e15.5,  (/,8x,3e15.5))" , &
                       f_scala = "(i8, e15.5,  (/,8x, e15.5))"
 INTERFACE
   INCLUDE 'tr4to8.h'
   INCLUDE 'cab_gid.h'
 END INTERFACE

 eset => sol3d_head

 !        print Gaussian variables

 DO
   gpname = 'GP'//eset%sname
   nelem = eset%nelem
   nnode = eset%nnode
   ngaus = eset%ngaus
   etype = eset%etype
   direct = .TRUE.
   SELECT CASE (nnode)
   CASE (4)      !Tetrahedra
     ng = 1             !position in array "O"
   CASE (6)      !Linear Prism
     SELECT CASE (ngaus)
     CASE (1)
       ! etype == 12 SOLSH
         ng =  1    !position in array "O"
                    ! GiD does not work with "Given" gauss points
         ngaus = 1  ! so the only option is to use "Internal",
         direct = .TRUE.
     CASE (2)
       IF( etype == 16 )THEN  !PRISM
         ng = 5             !position in array "O"
         ngaus = 6
         direct = .TRUE.
       ELSE IF( etype == 5 )THEN  !SPRISM
         ng = -1            !position in array "O"
         ngaus = 6
         direct = .FALSE.
       END IF
     CASE (6)
        ng = 5             !position in array "O"
        direct = .TRUE.
     END SELECT
   CASE (8,20)   !Hexahedra
     SELECT CASE (ngaus)
     CASE (1)
       ng = 1
     CASE (2)             !bottom and top surfaces
       ng = 0             !flag to interpolate
       ngaus = 8          !print in GiD positinos
       direct = .FALSE.
     CASE (4)             !selected 4 points
       IF( given )THEN    !print in actual positions
         ng = 3           !position in array "O"
       ELSE               !print in GiD positions
         ng = 6           !position in array "O"
         ngaus = 8
         direct = .FALSE.
       END IF
     CASE (8)
        ng = 6
     CASE (27)
        RETURN
     END SELECT
   CASE (15,18)    !Quadratic Prism
     ng = 5
     direct = .TRUE.
   END SELECT
   npg   = ngaus              ! the standard
   ! formats for Gauss-point WRITES
   WRITE(f_3dten(12:13),"(i2)")npg-1
   WRITE(f_3dvec(12:13),"(i2)")npg-1
   WRITE(f_scala(12:13),"(i2)")npg-1
   !
   iv = 0       !initializes pointer

   k = LEN_TRIM(gpname)
   IF(sol3d_stres > 1)THEN
     iv = iv+6
     CALL cab_gid(3,2,sol3d_stres_l(8),sol3d_stres_l(9),6,loadstep,ttime,gpname,units=sol3d_stres_u)
     DO ielem=1,nelem
       IF( direct )THEN
         WRITE(13,f_3dten)ielem,(eset%elvar(iv-5:iv,o(g,ng),ielem)*sol3d_stres_f, g=1,npg  )
       ELSE
         CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,6,ng)
         WRITE(13,f_3dten)ielem,(auxil(1:6,g)*sol3d_stres_f, g=1,npg  )
       END IF
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(sol3d_logst > 1)THEN
     iv = iv+3
     CALL cab_gid(2,2,sol3d_logst_l(5),sol3d_logst_l(6),3,loadstep,ttime,gpname,units=sol3d_logst_u)
     DO ielem=1,nelem
       IF( direct )THEN
         WRITE(13,f_3dvec)ielem,(eset%elvar(iv-2:iv,o(g,ng),ielem)*sol3d_logst_f, g=1,npg  )
       ELSE
         CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,3,ng)
         WRITE(13,f_3dvec)ielem,(auxil(1:3,g)*sol3d_logst_f, g=1,npg  )
       END IF
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(sol3d_thrat > 1)THEN
     iv = iv+1
     CALL cab_gid(1,2,sol3d_thrat_l(3),sol3d_thrat_l(4),1,loadstep,ttime,gpname,units=sol3d_thrat_u)
     DO ielem=1,nelem
       IF( direct )THEN
         WRITE(13,f_scala)ielem,(eset%elvar(iv,o(g,ng),ielem)*sol3d_thrat_f, g=1,npg  )
       ELSE
         CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,1,ng,0d0)
         WRITE(13,f_scala)ielem,(auxil(1,g)*sol3d_thrat_f, g=1,npg  )
       END IF
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(sol3d_eqpst > 1)THEN
     iv = iv+1
     CALL cab_gid(1,2,sol3d_eqpst_l(3),sol3d_eqpst_l(4),1,loadstep,ttime,gpname,units=sol3d_eqpst_u)
     DO ielem=1,nelem
       IF( direct )THEN
         WRITE(13,f_scala)ielem,(eset%elvar(iv,o(g,ng),ielem)*sol3d_eqpst_f, g=1,npg  )
       ELSE
         CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,1,ng,0d0)
         WRITE(13,f_scala)ielem,(auxil(1,g)*sol3d_eqpst_f, g=1,npg  )
       END IF
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(sol3d_vmise > 1)THEN
     iv = iv+1
     CALL cab_gid(1,2,sol3d_vmise_l(3),sol3d_vmise_l(4),1,loadstep,ttime,gpname,units=sol3d_vmise_u)
     DO ielem=1,nelem
       IF( direct )THEN
         WRITE(13,f_scala)ielem,(eset%elvar(iv,o(g,ng),ielem)*sol3d_vmise_f, g=1,npg  )
       ELSE
         CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,1,ng,0d0)
         WRITE(13,f_scala)ielem,(auxil(1,g)*sol3d_vmise_f, g=1,npg  )
       END IF
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(sol3d_fldma > 1)THEN
     iv = iv+1
     CALL cab_gid(1,2,sol3d_fldma_l(3),sol3d_fldma_l(4),1,loadstep,ttime,gpname)
     DO ielem=1,nelem
       IF( direct )THEN
         WRITE(13,f_scala)ielem,(eset%elvar(iv,o(g,ng),ielem), g=1,npg  )
       ELSE
         CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,1,ng,0d0)
         WRITE(13,f_scala)ielem,(auxil(1,g), g=1,npg)
       END IF
     END DO
     WRITE(13,"('End Values')")
   END IF

   !  user defined internal variables
   naux = eset%nstre - 8
   IF( naux > 0 )THEN
     matno = eset%matno(1)
     postd => udmat_head
     DO
       IF( postd%matno == matno )EXIT
       postd => postd%next
     END DO
     nvar = postd%nvar
     DO i=1,nvar
       SELECT CASE (postd%type(i))
       CASE (0) !scalar
         CALL cab_gid(1,2,postd%name(1,i),postd%name(2,i),1,            &
                      loadstep,ttime,gpname)
         iv = iv+1
         DO ielem=1,nelem
           IF( direct )THEN
             WRITE(13,f_scala)ielem,(eset%elvar(iv,o(g,ng),ielem),g=1,npg)
           ELSE
             CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,1,ng,0d0)
             WRITE(13,f_scala)ielem,(auxil(1,g),g=1,npg)
           END IF
         END DO
       CASE (1) !vector
         CALL cab_gid(2,2,postd%name(1,i),postd%name(2,i),postd%dim(i), &
                      loadstep,ttime,gpname)
         dim = postd%dim(i)
         j  = iv+1
         iv = iv+dim
         DO ielem=1,nelem
           IF( direct )THEN
             WRITE(13,"(i8,(3e15.5))")ielem, eset%elvar(j:iv,o(1,ng),ielem)
             DO g=2,npg
               WRITE(13,"(8x,3e15.5)")eset%elvar(j:iv,o(g,ng),ielem)
             END DO
           ELSE
             CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,dim,ng)
             WRITE(13,"(i8,(3e15.5))")ielem,auxil(1:dim,1)
             DO g=2,npg
               WRITE(13,"(8x,(3e15.5))")auxil(1:dim,g)
             END DO
           END IF
         END DO
       CASE (2) !tensor
         dim = postd%dim(i)
         j  = iv+1
         iv = iv+dim
         IF( dim == 3 )THEN ! a 2-d matrix
           CALL cab_gid(3,2,postd%name(1,i),postd%name(2,i),dim, &
                        loadstep,ttime,gpname)
           DO ielem=1,nelem
             WRITE(13,"(i8,6e15.5)")ielem, eset%elvar(j:j+1,1,ielem),0d0, eset%elvar(iv,1,ielem)
             !DO g=2,npg
             !  WRITE(13,"(8x,6e15.5)") eset%elvar(j:j+1,g,ielem),0d0, eset%elvar(iv,g,ielem)
             !END DO
           END DO
         ELSE
           CALL cab_gid(4,2,postd%name(1,i),postd%name(2,i),dim, &
                        loadstep,ttime,gpname)
           DO ielem=1,nelem
             IF( direct )THEN
               WRITE(13,"(i8,6e15.5)")ielem, eset%elvar(j:iv,o(1,ng),ielem)
               DO g=2,npg
                 WRITE(13,"(8x,6e15.5)")eset%elvar(j:iv,o(g,ng),ielem)
               END DO
             ELSE
               CALL tr4to8(eset%elvar(:,:,ielem),auxil,iv,dim,ng)
               WRITE(13,"(i8,6e15.5)")ielem,auxil(1:dim,1)
               DO g=2,npg
                 WRITE(13,"(8x,6e15.5)")auxil(1:dim,g)
               END DO
             END IF
           END DO
         END IF
       END SELECT
       WRITE(13,"('End Values')")
     END DO
   END IF

   eset => eset%next
   IF( .NOT.ASSOCIATED(eset) )EXIT
 END DO

 !        print Nodal variables

 IF( sol3d_nvarn == 0 )RETURN    !no variables to print, => exit

 nodes => sol3d_nodes
 vargs => sol3d_vargs
 iv = 0       !initializes pointer to vargs

 IF(MOD(sol3d_stres,2) == 1)THEN
   CALL cab_gid(3,1,sol3d_stres_l(1),sol3d_stres_l(2),6,loadstep,ttime,gpname,units=sol3d_stres_u)
   iv = iv+6
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(6e12.4))")label(i),vargs(iv-5:iv,is)*sol3d_stres_f
   END DO
   WRITE(13,"('End Values')")
   !  first invariant
   IF(MOD(sol3d_press,2) == 1)THEN
     CALL cab_gid(1,1,sol3d_press_l(1),sol3d_press_l(2),1,loadstep,ttime,gpname,units=sol3d_press_u)
     is = 0
     DO i=1,npoin
       IF(nodes(i,1) == 0)CYCLE
       is = is + 1
       press = -SUM(vargs(iv-5:iv-3,is))/3d0
       WRITE(13,"(i8,(1e12.4))")label(i),press*sol3d_stres_f
     END DO
     WRITE(13,"('End Values')")
   END IF
 END IF

 IF(MOD(sol3d_logst,2) == 1)THEN
   CALL cab_gid(2,1,sol3d_logst_l(1),sol3d_logst_l(2),3,loadstep,ttime,gpname,units=sol3d_logst_u)
   iv = iv+3
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv-2:iv,is)*sol3d_logst_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(sol3d_thrat,2) == 1)THEN
   CALL cab_gid(1,1,sol3d_thrat_l(1),sol3d_thrat_l(2),1,loadstep,ttime,gpname,units=sol3d_thrat_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv,is)*sol3d_thrat_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(sol3d_eqpst,2) == 1)THEN
   CALL cab_gid(1,1,sol3d_eqpst_l(1),sol3d_eqpst_l(2),1,loadstep,ttime,gpname,units=sol3d_eqpst_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv,is)*sol3d_eqpst_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(sol3d_vmise,2) == 1)THEN
   CALL cab_gid(1,1,sol3d_vmise_l(1),sol3d_vmise_l(2),1,loadstep,ttime,gpname,units=sol3d_vmise_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv,is)*sol3d_vmise_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(sol3d_fldma,2) == 1)THEN
   CALL cab_gid(1,1,sol3d_fldma_l(1),sol3d_fldma_l(2),1,loadstep,ttime,gpname)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv,is)
   END DO
   WRITE(13,"('End Values')")
 END IF

 RETURN
 END SUBROUTINE prgi05
