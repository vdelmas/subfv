module precision_module
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: QUAD = real128
  integer, parameter :: DOUBLE = real64
  integer, parameter :: ENTIER = int32
  integer, parameter :: ENTIER_D = int64
  integer, parameter :: ENTIER_8 = int8
  integer, parameter :: ENTIER_16 = int16
end module precision_module
