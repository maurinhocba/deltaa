SUBROUTINE fixtem(iwrit,npoin,iftmp,label)

  !***  APPLY Temperature prescribed values

  USE c_input
  USE ift_db
  USE curv_db, ONLY: getcun
  IMPLICIT NONE
  INTEGER (kind=4), INTENT(IN) :: iwrit,npoin,label(:)
  INTEGER (kind=4), INTENT(IN OUT) :: iftmp(:,:)

  INTEGER (kind=4) :: ifix(3),g,i,j,l,n,ipoin,chnode,np,nposn,nrve
  INTEGER (kind=4), PARAMETER :: nn = 1000000, nn1 = 1000001
  REAL (kind=8) xg(ndoft),xf(ndoft)
  REAL (kind=8), ALLOCATABLE :: auxil(:,:)
  TYPE (rpt_set), POINTER :: rves
  TYPE (rpt_nod), POINTER :: rven


  IF(iwrit == 1) THEN
    IF(ndoft == 1) WRITE(lures,"(/,' Boundary Conditions', &
                               & /,'   Node    T')",ERR=9999)
    IF(ndoft == 2) WRITE(lures,"(/,' Boundary Conditions', &
                               & /,'   Node    TI TS')",ERR=9999)
    IF(ndoft == 3) WRITE(lures,"(/,' Boundary Conditions', &
                               & /,'   Node    TN TI TS')",ERR=9999)
  END IF


  !***  Applies prescribed temperatures

  nprev = 0    !initializes

  IF (npret > 0) THEN  !if prescribed temperatures exist

    ALLOCATE( lctmp(npret+1) )                !space for associated curves
    ALLOCATE( auxil(0:ndoft*npoin,npret) )    !auxiliar space
    auxil = 0d0                               !initializes

    rves => headv                             !point to first set

    DO l=1,npret                   !loop on each set

      lctmp(l)  = getcun (rves%lc)  ! get associated curve number for this set
      auxil(0,l)= rves%factor       !associated factor for this set

      IF(iwrit == 1) THEN
        WRITE(lures, "(//, &
        & 5X,'Curve scaling this temperature set .',i10/ &
        & 5X,'Scaling factor .....................',e14.7,/)",ERR=9999) &
          lctmp(l),auxil(0,l)

        WRITE(lures, "(//,5X,'Prescribed temperatures',//)",ERR=9999)
        IF (ndoft == 1 ) WRITE(lures, &
           "(6X,'NODE',6X,'Temp.')",ERR=9999)
        IF (ndoft == 2 ) WRITE(lures, &
           "(6X,'NODE',6X,'Temp-S',8X,'Temp-I')",ERR=9999)
        IF (ndoft == 3 ) WRITE(lures, &
           "(6X,'NODE',6X,'Temp-M',8X,'Temp-I',8X,'Temp-S')",ERR=9999)
      END IF

      nrve = rves%nrv         !number of nodes in this set
      rven => rves%head       !point to first node in the set
      DO j=1,nrve             !loop for each node in the set
        g = rven%node         !node label
        n = chnode(g)         !internal number

        xg(1:ndoft) = rven%v(1:ndoft)  !prescribed temperatures components
        xf = xg                        !copy onto XF
        DO i=1,ndoft                   !for each DOF
          nposn = iftmp(i,n)           !restriction code
          SELECT CASE (nposn)
          CASE (:-nn1)                        !prescribed values
            auxil(-nposn+nn,l) = xg(i)        !assign value
          CASE (0)                            !active DOF or non existent
            IF(xg(i) /= 0d0)THEN              !if a non-zero value
              nprev = nprev+1                 !increase number of fixed values
              iftmp(i,n) = -nprev-nn          !modify restriction code
              auxil(nprev,l) = xg(i)          !assign value
            END IF
          END SELECT
        END DO
        ! echo effectively assigned values
        IF(iwrit == 1)  WRITE(lures,"(i10,3e14.5)",ERR=9999) label(n), xf(1:ndoft)
        rven => rven%next        !point to next node
      END DO
      rves => rves%next       !point to next set
    END DO

    ALLOCATE( prtmp(nprev+1,npret+1) )     !reserve space for definitive array
    prtmp(1:nprev,1:npret) = auxil(1:nprev,1:npret)  !transfer temperature data
    prtmp(nprev+1,1:npret) = auxil(0,1:npret)        !transfer factors
    prtmp(1:nprev,npret+1) = 0.d0  !initializes compounded prescribed temperatures
    DEALLOCATE ( auxil )           !release auxiliar array

  ELSE  !no prescribed temperatures

    ALLOCATE( prtmp(1,1) ) !to avoid null pointers only
    ALLOCATE( lctmp(npret+1) )
    prtmp = 0d0

  END IF

RETURN
 9999 CALL runen2('')
END SUBROUTINE fixtem
