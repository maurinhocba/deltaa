MODULE param_db
! Module where are defined global parameters of Stampack
IMPLICIT NONE

  !--- Constants for character lenght
  INTEGER(kind=4),PARAMETER:: mlen = 256,    & !Maximum length for input parameters and file names
                              mlng = 1024      !Maximum length for long strings

END MODULE
