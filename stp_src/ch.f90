FUNCTION ch(istep)
!
! function to convert a 3 digit integer into a string
!
IMPLICIT NONE
CHARACTER*3 ch
INTEGER (kind=4) istep
CHARACTER*1 first,second,third

first  = char(MOD(istep/100,10)+48)
second = char(MOD(istep/10,10)+48)
third  = char(MOD(istep,10)+48)

ch = first//second//third

END FUNCTION ch
