  SUBROUTINE dmat06(dm,d)
  IMPLICIT NONE
  REAL(kind=8) :: dm(:),d(:)
  dm( 1) = d( 1); dm( 2) = d( 2); dm( 3) = d( 3); dm( 4) = d(16); dm( 5) = d(17); dm( 6) = d(18)
                  dm(15) = d( 4); dm(16) = d( 5); dm(17) = d(19); dm(18) = d(20); dm(19) = d(21)
                                  dm(28) = d( 6); dm(29) = d(22); dm(30) = d(23); dm(31) = d(24)
                                                  dm(40) = d( 7); dm(41) = d( 8); dm(42) = d( 9)
                                                                  dm(51) = d(10); dm(52) = d(11)
                                                                                  dm(61) = d(12)
                   dm( 7: 8) = 0d0; dm( 9) = d(35); dm(10) = d(36); dm(11) = d(37); dm(12) = d(38); dm(13:14) = 0d0
                   dm(20:21) = 0d0; dm(22) = d(39); dm(23) = d(40); dm(24) = d(41); dm(25) = d(42); dm(26:27) = 0d0
                   dm(32:33) = 0d0; dm(34) = d(43); dm(35) = d(44); dm(36) = d(45); dm(37) = d(46); dm(38:39) = 0d0
                   dm(43:44) = 0d0; dm(45) = d(47); dm(46) = d(48); dm(47) = d(49); dm(48) = d(50); dm(49:50) = 0d0
                   dm(53:54) = 0d0; dm(55) = d(51); dm(56) = d(52); dm(57) = d(53); dm(58) = d(54); dm(59:60) = 0d0
                   dm(62:63) = 0d0; dm(64) = d(55); dm(65) = d(56); dm(66) = d(57); dm(67) = d(58); dm(68:69) = 0d0
  dm(70) = d(13); dm(71) = d(14); dm(72:75) = 0d0;                                                dm(76) = d(62); dm(77) = d(63)
                  dm(78) = d(15); dm(79:82) = 0d0;                                                dm(83) = d(64); dm(84) = d(65)
                                  dm(85) = d(25); dm(86) = d(26); dm(87) = d(27); dm(88) = d(28); dm(89:89) = 0d0
                                                  dm(91) = d(29); dm(92) = d(30); dm(93) = d(31); dm(94:95) = 0d0
                                                                  dm(96) = d(32); dm(97) = d(33); dm(98:99) = 0d0
                                                                                  dm(100)= d(34); dm(101:102)= 0d0
                                                                                                  dm(103)= d(59); dm(104)= d(60)
                                                                                                                  dm(105)= d(61)
  RETURN
  END SUBROUTINE dmat06

! dm( 1) = d( 1)
! dm( 2) = d( 2)
! dm( 3) = d( 3)
! dm( 4) = d(16)
! dm( 5) = d(17)
! dm( 6) = d(18)
! dm(7:8)= 0d0
! dm( 9) = d(35)
! dm(10) = d(36)
! dm(11) = d(37)
! dm(12) = d(38)
! dm(13:14)= 0d0
! dm(15) = d( 4)
! dm(16) = d( 5)
! dm(17) = d(19)
! dm(18) = d(20)
! dm(19) = d(21)
! dm(20:21)= 0d0
! dm(22) = d(39)
! dm(23) = d(40)
! dm(24) = d(41)
! dm(25) = d(42)
! dm(26:27)= 0d0
! dm(28) = d( 6)
! dm(29) = d(22)
! dm(30) = d(23)
! dm(31) = d(24)
! dm(32:33)= 0d0
! dm(34) = d(43)
! dm(35) = d(44)
! dm(36) = d(45)
! dm(37) = d(46)
! dm(38:39) = 0d0
! dm(40) = d( 7)
! dm(41) = d( 8)
! dm(42) = d( 9)
! dm(43:44) = 0d0
! dm(45) =  d(47)
! dm(46) =  d(48)
! dm(47) =  d(49)
! dm(48) =  d(50)
! dm(49:50) = 0d0
! dm(51) = d(10)
! dm(52) = d(11)
! dm(53:54) = 0d0
! dm(55) =  d(51)
! dm(56) =  d(52)
! dm(57) =  d(53)
! dm(58) =  d(54)
! dm(59:60) = 0d0
! dm(61) = d(12)
! dm(62:63) = 0d0
! dm(64) =  d(55)
! dm(65) =  d(56)
! dm(66) =  d(57)
! dm(67) =  d(58)
! dm(68:69) = 0d0
! dm(70) =  d(13)
! dm(71) =  d(14)
! dm(72:75) = 0d0
! dm(76) =  d(62)
! dm(77) =  d(63)
! dm(78) =  d(15)
! dm(79:82) = 0d0
! dm(83) =  d(64)
! dm(84) =  d(65)
! dm(85) = d(25)
! dm(86) = d(26)
! dm(87) = d(27)
! dm(88) = d(28)
! dm(89:89)= 0d0
! dm(91) = d(29)
! dm(92) = d(30)
! dm(93) = d(31)
! dm(94:95)= 0d0
! dm(96) = d(32)
! dm(97) = d(33)
! dm(98:99)= 0d0
! dm(100)= d(34)
! dm(101:102)= 0d0
! dm(103)= d(59)
! dm(104)= d(60)
! dm(105)= d(61)
