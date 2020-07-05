! Placeholder for Fortran wrapper interface provided by OpenMM.
! https://simtk.org/home/openmm

#if KEY_OPENMM==1

! pgf95 wants 'end subroutine', not just 'end' as in OpenMM 4.1 wrapper

INCLUDE 'OpenMMFortranModule.f90'
INCLUDE 'CharmmOpenMMFortranModule.f90'
INCLUDE 'OpenMMGBSWFortranModule.f90'

#else /* KEY_OPENMM */

! humor setmk.com

MODULE OpenMM_Types
END MODULE OpenMM_Types

MODULE OpenMM
END MODULE OpenMM

MODULE OpenMM_Charmm
END MODULE OpenMM_Charmm

MODULE OpenMMGBSW_Types
END MODULE OpenMMGBSW_Types

MODULE OpenMMGBSW
END MODULE OpenMMGBSW
#endif
