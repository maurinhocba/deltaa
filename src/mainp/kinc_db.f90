 MODULE kinc_db
   USE ctrl_db, ONLY : ndime,neq
   ! kinematic constraints data base
   IMPLICIT NONE
   SAVE
   INTEGER (kind=4) :: &
     naris=0, & !number of nodes constrained to remain over a side (2-D)
     ndepd, & !number of total dependent nodes (defined as pairs + rigid bodies)
     ndumn, & !number of dependent nodes in rigid bodies
     nescv=0, & !number of slave constrains
     nsdof, & !total number of slave DOFs
     nvelr=0, & !number of constrained velocity sets
     nvfix=0, & !number of fixed values
     maxa     !size or stiffness matrix

   INTEGER (kind=4), POINTER :: &
     nndpd(:,:), & !(3,ndepd) slave-master full dependent nodes
     sldnd(:,:)    !(3,nsldn) slave-master & flags sliding nodes
   REAL (kind=8), POINTER :: &
     distd(:,:), & !(3,ndepd) distance vector between dependant nodes
     sldvl(:,:)    !(2,nsldn) distance and rotation values for sliding nodes

   INTEGER (kind=4),POINTER :: &
     lcvel(:), & !
     nesdf(:), & !( ) active equations associated to slave DOFs
     npsdf(:)    !(nsdof) pointer to array NESDF & FTSDF
   REAL (kind=8),POINTER :: &
     velor(:,:), & !(nvfix+1,nvelr+1)  prescribed velocities
     ftsdf(:),   & !( ) factors associated to slave DOFs
     ftsd0(:)      !( ) factors associated to slave DOFs
   INTEGER (kind=4), POINTER :: &
     maxav(:)      !(neq+1) position of diagonal terms on stiffness matrix

   INTEGER (kind=4), PARAMETER :: nn = 1000000, nn1=1000001, nm=999999
   INTEGER (kind=4) :: &
     nv1      !position of actual velocity set


 CONTAINS
   SUBROUTINE dump_kinc (nnsld)
   IMPLICIT NONE
   INTEGER(Kind=4) :: nnsld
   INTEGER(Kind=4) :: i,j,k

   j = SIZE(nesdf)
   WRITE(50,ERR=9999) naris, ndepd, ndumn, nescv, nsdof, &
                      nvelr, nvfix, nnsld, j
   WRITE(50,ERR=9999) lcvel, nndpd, npsdf    !kinc_db
   WRITE(50,ERR=9999) (nesdf(i),i=1,j)
   WRITE(50,ERR=9999) ((distd(i,k),i=1,3),k=1,MAX(ndepd,1)),ftsdf
   IF( nnsld > 0 )WRITE(50,ERR=9999) sldnd,sldvl
   WRITE(50,ERR=9999) ((velor(i,k),i=1,nvfix+1),k=1,nvelr+1)
   WRITE(50,ERR=9999) (maxav(i),i=1,neq+1)
   RETURN
   9999 CALL runen2('')
   END SUBROUTINE dump_kinc

   SUBROUTINE rest_kinc (nnsld)
   IMPLICIT NONE
   INTEGER(Kind=4) :: nnsld
   INTEGER(Kind=4) :: i,j,k

   READ (51) naris, ndepd, ndumn, nescv, nsdof, nvelr, nvfix, nnsld, j
   ALLOCATE( lcvel(nvelr+1),nndpd(3,MAX(ndepd+naris,1)),npsdf(nsdof+1),nesdf(j), &
             ftsdf(j),     &
             distd(3,MAX(ndepd,1)),velor(nvfix+1,nvelr+1) )

   READ (51) lcvel, nndpd, npsdf
   READ (51) (nesdf(i),i=1,j)
   READ (51) ((distd(i,k),i=1,3),k=1,MAX(ndepd,1)),(ftsdf(i),i=1,j)
   IF( nnsld > 0 )THEN
     ALLOCATE( sldnd(3,nnsld),sldvl(2,nnsld))
     READ(51) sldnd,sldvl
   END IF
   READ (51) ((velor(i,k),i=1,nvfix+1),k=1,nvelr+1)
   ALLOCATE( maxav(neq+1))
   READ(51) (maxav(i),i=1,neq+1)

   END SUBROUTINE rest_kinc
 END MODULE kinc_db

