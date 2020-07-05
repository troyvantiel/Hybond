module version
  use chm_kinds

  !
  !     VERNUM  - version number (simple integer); 0 < VERNUM < 10000 is assumed 
  !               for endian check. /LNI February 2013
  !     VERNMC  - version number (character string)
  !
  INTEGER, PARAMETER :: VERNUM=44
  CHARACTER(len=24), PARAMETER :: VERNMC='44b2   February 15, 2020'
  !                                       123456789+123456789+1234
  logical, parameter :: free_version = .true.
end module version

