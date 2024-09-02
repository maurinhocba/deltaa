 SUBROUTINE elmdat (task,elsnam,itype)

 !     read element sets

 USE ctrl_db, ONLY: ndime, neulr, npoin
 USE outp_db, ONLY:  iwrit
 USE esets_db, ONLY: nelms,nel
 USE npo_db, ONLY:   label, coord, eule0, euler, cpx
 IMPLICIT NONE

 CHARACTER(len=*),INTENT(IN) :: elsnam
 CHARACTER(len=*),INTENT(IN) :: task
 INTEGER (kind=4), INTENT(IN) :: itype

 INTERFACE
   INCLUDE 'inpda1.h'
   INCLUDE 'inpda2.h'
   INCLUDE 'inpd03.h'
   INCLUDE 'inpd04.h'
   INCLUDE 'inpd05.h'
   INCLUDE 'inpda6.h'
   INCLUDE 'inpda7.h'
   INCLUDE 'inpda8.h'
   INCLUDE 'inpda9.h'
   INCLUDE 'inpd10.h'
   INCLUDE 'inpd11.h'
   INCLUDE 'inpd12.h'
   INCLUDE 'inpd13.h'
   INCLUDE 'inpd14.h'
   INCLUDE 'inpd15.h'
   INCLUDE 'inpd16.h'
   INCLUDE 'inpd17.h'
   INCLUDE 'inpd18.h'
   INCLUDE 'inpd19.h'
   INCLUDE 'inpd20.h'
   INCLUDE 'inpd24.h'
   INCLUDE 'inpd25.h'
   INCLUDE 'inpd26.h'
   INCLUDE 'inpd27.h'
   INCLUDE 'inpd29.h'
   INCLUDE 'inpd30.h'
 END INTERFACE

 SELECT CASE (itype)
 CASE (1)
   CALL inpda1(task,ndime,neulr,nel(1),iwrit,elsnam,nelms(1))
 CASE (2)
   CALL inpda2(task,ndime,nel(2),iwrit,elsnam,nelms(2))
 CASE (3)
   CALL inpd03(task,nel(3),eule0,euler,coord,iwrit,elsnam,nelms(3))
 CASE (4)
   CALL inpd04(task,nel(4),iwrit,elsnam,nelms(4))
 CASE (5)
   CALL inpd05(task,nel(5),iwrit,elsnam,nelms(5))
 CASE (6)
   CALL inpda6(task,nel(6),eule0,euler,coord,iwrit,elsnam,nelms(6))
 CASE (7)
   IF( .NOT.ASSOCIATED (cpx) )THEN
     ALLOCATE (cpx(3,npoin))
     cpx = 0
   END IF
   CALL inpda7(task,nel(7),eule0,euler,coord,iwrit,elsnam,nelms(7))
 CASE (8)
   CALL inpda8(task,nel(8),eule0,euler,coord,iwrit,elsnam,nelms(8))
 CASE (9)
   CALL inpda9(task,nel(9), eule0,euler,coord, iwrit,elsnam,nelms(9))
 CASE (10)
   CALL inpd10(task,ndime, nel(10), iwrit,elsnam,nelms(10))
 CASE (11)
   CALL inpd11(task,nel(11),iwrit,elsnam,nelms(11))
 CASE (12)
   CALL inpd12(task,nel(12),iwrit,elsnam,nelms(12))
 CASE (13)
   CALL inpd13(task,nel(13),iwrit,elsnam,nelms(13))
 CASE (14)
   CALL inpd14(task,nel(14),iwrit,elsnam,nelms(14))
 CASE (15)
   CALL inpd15(task,nel(15),iwrit,elsnam,nelms(15))
 CASE (16)
   CALL inpd16(task,nel(16),iwrit,elsnam,nelms(16))
 CASE (17)
   CALL inpd17(task,nel(17),iwrit,elsnam,nelms(17))
 CASE (18)
   CALL inpd18(task,nel(18),iwrit,elsnam,nelms(18))
 CASE (19)
   IF( .NOT.ASSOCIATED (cpx) )THEN
     ALLOCATE (cpx(3,npoin))
     cpx = 0
   END IF
   CALL inpd19(task,nel(19),iwrit,elsnam,nelms(19))
 CASE (20)
   CALL inpd20(task,nel(20),iwrit,elsnam,nelms(20))
 CASE (24)
   CALL inpd24(task,nel(24),iwrit,elsnam,nelms(24))
 CASE (25)
   CALL inpd25(task,nel(25),iwrit,elsnam,nelms(25))
 CASE (26)
   CALL inpd26(task,nel(26),iwrit,elsnam,nelms(26))
 CASE (27)
   IF( .NOT.ASSOCIATED (cpx) )THEN
     ALLOCATE (cpx(3,npoin))
     cpx = 0
   END IF
   CALL inpd27(task,nel(27),iwrit,elsnam,nelms(27))
 CASE (29)
   IF( .NOT.ASSOCIATED (cpx) )THEN
     ALLOCATE (cpx(3,npoin))
     cpx = 0
   END IF
   CALL inpd29(task,nel(29),iwrit,elsnam,nelms(29))
 CASE (30)
   CALL inpd30(task,nel(30),iwrit,elsnam,nelms(30))
 CASE DEFAULT
   CALL runend('ELMDAT: ELEMENT_TYPE NOT EXISTENT  ')
 END SELECT

 RETURN
 END SUBROUTINE elmdat
