module cstran_mod
  use chm_kinds
  implicit none

  !
  !  Center of mass harmonic restraints
  !     QHMCM      Flag indicating that potential is in use
  !     HMCMX      Reference centers of mass for each atom subset
  !     HMCMY
  !     HMCMZ
  !     KHMCM      Force constants for each restraint subset
  !     RHMCM      Reference distance for each restraint subset
  !     NHMCM      Number of restraint subsets
  !     IHMCM      Pointer to index list for restraint subset
  !     INHMCM     Number of atoms in a restraint subset
  !     PHMCM      Index list for restraint subset
  !     NPHMCM     Number of indices in PHMCM
  !     LHMCMM     Logical flag for mass weighting of the center of mass
  !
  !     NHMCMR     Number of relative distance restraints
  !     KHMCMR     Force constants for relative restraints
  !     HMCMR      Reference distance between centers of mass
  !     IHMCM1     First center of mass
  !     IHMCM2     Second center of mass
  !
  logical :: qhmcm
#if KEY_HMCOM==1
!  integer,parameter :: MAXHCM=500
!  integer,parameter :: MAXPCM=5000
!  integer,parameter :: MAXRCM=2000

  integer,parameter,public :: HMCM_X=1,HMCM_Y=2,HMCM_Z=4,HMCM_XY=3,HMCM_XZ=5,HMCM_YZ=6,HMCM_XYZ=7 
  integer,save,public ::  nhmcm, nphmcm, nhmcmr, maxhcm, maxpcm, maxrcm
  integer,allocatable, dimension(:),save,public :: ihmcm,inhmcm,phmcm,ihmcm1,ihmcm2,hmcmtype
       
  real(chm_real),allocatable, dimension(:),save,public :: hmcmx,hmcmy,hmcmz,khmcm,khmcmr,hmcmr,rhmcm

  logical,allocatable, dimension(:),save,public :: lhmcmm

!  integer :: nhmcm,ihmcm(maxhcm),inhmcm(maxhcm),phmcm(maxpcm), &
!       nphmcm,nhmcmr,ihmcm1(maxrcm),ihmcm2(maxrcm), &
!       hmcmtype(maxhcm)
       

!  real(chm_real) :: hmcmx(maxhcm),hmcmy(maxhcm),hmcmz(maxhcm), &
!       khmcm(maxhcm),khmcmr(maxrcm),hmcmr(maxrcm), &
!       rhmcm(maxhcm)

!  logical :: lhmcmm(maxhcm)

#endif 

character(len=*), private,parameter :: file_name="cstran.src"

contains

#if KEY_HMCOM==1 /*hmcm_init*/
  subroutine hmcm_init
    !     center of mass constraints
    qhmcm=.false.
    nhmcm=0
    nphmcm=0
    nhmcmr=0
    maxhcm=0
    maxpcm=0
    maxrcm=0
    return
  end subroutine hmcm_init

  subroutine allocate_hmcm
    use memory
    use number

    integer, allocatable, dimension(:) :: iwork
    logical, allocatable, dimension(:) :: lwork
    integer  :: olddim
    character(len=*), parameter :: routine_name="allocate_hmcm"

    real(chm_real), allocatable, dimension(:) :: work

    if(.not. allocated(ihmcm).and. nhmcm>maxhcm) then  ! allocate space for hmcm restraints
       call chmalloc(file_name,routine_name,'ihmcm',nhmcm,intg=ihmcm)
       call chmalloc(file_name,routine_name,'inhmcm',nhmcm,intg=inhmcm)
       call chmalloc(file_name,routine_name,'hmcmtype',nhmcm,intg=hmcmtype)
       call chmalloc(file_name,routine_name,'lhmcmm',nhmcm,log=lhmcmm)
       call chmalloc(file_name,routine_name,'khmcm',nhmcm,crl=khmcm)
       call chmalloc(file_name,routine_name,'rhmcm',nhmcm,crl=rhmcm)
       call chmalloc(file_name,routine_name,'hmcmx',nhmcm,crl=hmcmx)
       call chmalloc(file_name,routine_name,'hmcmy',nhmcm,crl=hmcmy)
       call chmalloc(file_name,routine_name,'hmcmz',nhmcm,crl=hmcmz)
       maxhcm = nhmcm
       ihmcm(1:maxhcm)=0
       inhmcm(1:maxhcm)=0
       hmcmtype(1:maxhcm)=0
       lhmcmm(1:maxhcm)=.false.
       khmcm(1:maxhcm)=zero
       rhmcm(1:maxhcm)=zero
       hmcmx(1:maxhcm)=zero
       hmcmy(1:maxhcm)=zero
       hmcmz(1:maxhcm)=zero
    elseif(nhmcm>maxhcm) then  !resize space for hmcm restraints
       olddim=size(ihmcm)
       call chmalloc(file_name,routine_name,'iwork',olddim,intg=iwork)
       call chmalloc(file_name,routine_name,'lwork',olddim,log=lwork)
       call chmalloc(file_name,routine_name,'work',olddim,crl=work)

       iwork=ihmcm
       call chmdealloc(file_name,routine_name,'ihmcm',olddim,intg=ihmcm)
       call chmalloc(file_name,routine_name,'ihmcm',nhmcm,intg=ihmcm)
       ihmcm(1:olddim)=iwork(1:olddim)

       iwork=inhmcm
       call chmdealloc(file_name,routine_name,'inhmcm',olddim,intg=inhmcm)
       call chmalloc(file_name,routine_name,'inhmcm',nhmcm,intg=inhmcm)
       inhmcm(1:olddim)=iwork(1:olddim)

       iwork=hmcmtype
       call chmdealloc(file_name,routine_name,'hmcmtype',olddim,intg=hmcmtype)
       call chmalloc(file_name,routine_name,'hmcmtype',nhmcm,intg=hmcmtype)
       hmcmtype(1:olddim)=iwork(1:olddim)

       lwork=lhmcmm
       call chmdealloc(file_name,routine_name,'lhmcmm',olddim,log=lhmcmm)
       call chmalloc(file_name,routine_name,'lhmcmm',nhmcm,log=lhmcmm)
       lhmcmm(1:olddim)=lwork(1:olddim)

       work=khmcm
       call chmdealloc(file_name,routine_name,'khmcm',olddim,crl=khmcm)
       call chmalloc(file_name,routine_name,'khmcm',nhmcm,crl=khmcm)
       khmcm(1:olddim)=work(1:olddim)

       work=rhmcm
       call chmdealloc(file_name,routine_name,'rhmcm',olddim,crl=rhmcm)
       call chmalloc(file_name,routine_name,'rhmcm',nhmcm,crl=rhmcm)
       rhmcm(1:olddim)=work(1:olddim)

       work=hmcmx
       call chmdealloc(file_name,routine_name,'hmcmx',olddim,crl=hmcmx)
       call chmalloc(file_name,routine_name,'hmcmx',nhmcm,crl=hmcmx)
       hmcmx(1:olddim)=work(1:olddim)

       work=hmcmy
       call chmdealloc(file_name,routine_name,'hmcmy',olddim,crl=hmcmy)
       call chmalloc(file_name,routine_name,'hmcmy',nhmcm,crl=hmcmy)
       hmcmy(1:olddim)=work(1:olddim)

       work=hmcmz
       call chmdealloc(file_name,routine_name,'hmcmz',olddim,crl=hmcmz)
       call chmalloc(file_name,routine_name,'hmcmz',nhmcm,crl=hmcmz)
       hmcmz(1:olddim)=work(1:olddim)

       call chmdealloc(file_name,routine_name,'iwork',olddim,intg=iwork)
       call chmdealloc(file_name,routine_name,'lwork',olddim,log=lwork)
       call chmdealloc(file_name,routine_name,'work',olddim,crl=work)
       maxhcm=nhmcm

    elseif(.not.allocated(phmcm)) then !allocate space for CM atom pointers
       maxpcm=nphmcm
       call chmalloc(file_name,routine_name,'phmcm',maxpcm,intg=phmcm)
       phmcm=0
    elseif(nphmcm>maxpcm) then  !rsize space for CM atom pointers
       olddim=size(phmcm)
       call chmalloc(file_name,routine_name,'iwork',olddim,intg=iwork)
       iwork=phmcm
       call chmdealloc(file_name,routine_name,'phmcm',olddim,intg=phmcm)
       call chmalloc(file_name,routine_name,'phmcm',nphmcm,intg=phmcm)
       phmcm(1:olddim)=iwork(1:olddim)
       call chmdealloc(file_name,routine_name,'iwork',olddim,intg=iwork)
       maxpcm=nphmcm

    elseif(.not.allocated(khmcmr).and.nhmcmr>maxrcm) then  !allocate space for hmcr restraints
       call chmalloc(file_name,routine_name,'khmcmr',nhmcmr,crl=khmcmr)
       call chmalloc(file_name,routine_name,'hmcmr',nhmcmr,crl=hmcmr)
       call chmalloc(file_name,routine_name,'ihmcm1',nhmcmr,intg=ihmcm1)
       call chmalloc(file_name,routine_name,'ihmcm2',nhmcmr,intg=ihmcm2)
       khmcmr=zero
       hmcmr=zero
       ihmcm1=0
       ihmcm2=0
       maxrcm=nhmcmr
    elseif(nhmcmr>maxrcm) then  !resize hmcr restraints
       olddim=size(khmcmr)
       call chmalloc(file_name,routine_name,'iwork',olddim,intg=iwork)
       call chmalloc(file_name,routine_name,'work',olddim,crl=work)

       work=khmcmr
       call chmdealloc(file_name,routine_name,'khmcmr',olddim,crl=khmcmr)
       call chmalloc(file_name,routine_name,'khmcmr',nhmcmr,crl=khmcmr)
       khmcmr(1:olddim)=work(1:olddim)

       work=hmcmr
       call chmdealloc(file_name,routine_name,'hmcmr',olddim,crl=hmcmr)
       call chmalloc(file_name,routine_name,'hmcmr',nhmcmr,crl=hmcmr)
       hmcmr(1:olddim)=work(1:olddim)

       iwork=ihmcm1
       call chmdealloc(file_name,routine_name,'ihmcm1',olddim,intg=ihmcm1)
       call chmalloc(file_name,routine_name,'ihmcm1',nhmcmr,intg=ihmcm1)
       ihmcm1(1:olddim)=iwork(1:olddim)

       iwork=ihmcm2
       call chmdealloc(file_name,routine_name,'ihmcm2',olddim,intg=ihmcm2)
       call chmalloc(file_name,routine_name,'ihmcm2',nhmcmr,intg=ihmcm2)
       ihmcm2(1:olddim)=iwork(1:olddim)

       maxrcm=nhmcmr
       call chmdealloc(file_name,routine_name,'iwork',olddim,intg=iwork)
       call chmdealloc(file_name,routine_name,'work',olddim,crl=work)

    endif
    return
  end  subroutine allocate_hmcm

  subroutine deallocate_hmcm
    use memory
    integer :: olddim
    character(len=*), parameter :: routine_name="deallocate_hmcm"

    if(allocated(ihmcm).and..not.qhmcm) then  !process hcm clear command
       olddim=size(ihmcm)
       call chmdealloc(file_name,routine_name,'ihmcm',olddim,intg=ihmcm)
       call chmdealloc(file_name,routine_name,'inhmcm',olddim,intg=inhmcm)
       call chmdealloc(file_name,routine_name,'hmcmtype',olddim,intg=hmcmtype)
       call chmdealloc(file_name,routine_name,'lhmcmm',olddim,log=lhmcmm)
       call chmdealloc(file_name,routine_name,'khmcm',olddim,crl=khmcm)
       call chmdealloc(file_name,routine_name,'rhmcm',olddim,crl=rhmcm)
       call chmdealloc(file_name,routine_name,'hmcmx',olddim,crl=hmcmx)
       call chmdealloc(file_name,routine_name,'hmcmy',olddim,crl=hmcmy)
       call chmdealloc(file_name,routine_name,'hmcmz',olddim,crl=hmcmz)
       maxhcm=0

       if(allocated(phmcm)) then
          olddim=size(phmcm)
          call chmdealloc(file_name,routine_name,'phmcm',olddim,intg=phmcm)
          maxpcm=0
       endif

       if(allocated(khmcmr)) then
          olddim=size(khmcmr)
          call chmdealloc(file_name,routine_name,'khmcmr',olddim,crl=khmcmr)
          call chmdealloc(file_name,routine_name,'hmcmr',olddim,crl=hmcmr)
          call chmdealloc(file_name,routine_name,'ihmcm1',olddim,intg=ihmcm1)
          call chmdealloc(file_name,routine_name,'ihmcm2',olddim,intg=ihmcm2)
          maxrcm=0
       endif
    elseif(allocated(ihmcm).and.allocated(khmcmr)) then  !process hcmr clear command
       olddim=size(khmcmr)
       call chmdealloc(file_name,routine_name,'khmcmr',olddim,crl=khmcmr)
       call chmdealloc(file_name,routine_name,'hmcmr',olddim,crl=hmcmr)
       call chmdealloc(file_name,routine_name,'ihmcm1',olddim,intg=ihmcm1)
       call chmdealloc(file_name,routine_name,'ihmcm2',olddim,intg=ihmcm2)
       maxrcm=0
    endif
    return
  end   subroutine deallocate_hmcm

#endif /* (hmcm_init)*/


  SUBROUTINE CSTRAN(COMLYN,COMLEN)
    !
    !     Process the CONSTRAINT command.
    !

  use cons_rmsd,only: setrmsd0
  use cnst_fcm,only: allocate_cnst

  use chm_kinds
  use dimens_fcm
  use code
  use bases_fcm
  use coord
  use coordc
#if KEY_EMAP==1
  use emapmod,only:emaprest
#endif
  use psf
  use memory
  use stream
  use string

#if KEY_CFF==1 || KEY_MMFF==1
  use ffieldm

#endif 
#if KEY_CPATH==1
  use deriv    
#endif

    implicit none

    CHARACTER COMLYN*(*)
    INTEGER   COMLEN
    ! . Local variables.
    LOGICAL   LCOMP, MSTUP
    integer,allocatable,dimension(:),target :: ispace
    integer,pointer,dimension(:) :: islct,jslct


    call chmalloc('cstran.src','CSTRAN','ispace',2*NATOM,intg=ispace)
    islct => ispace(1:natom)
    jslct => ispace(natom+1:2*natom)
    MSTUP = .FALSE.

    LCOMP = INDXA(COMLYN, COMLEN, 'COMP')  >  0

    if(indxa(comlyn, comlen, 'RMSD')  >  0)then
       call allocate_cnst(natom)
       call setrmsd0(islct,jslct,lcomp,x,y,z,wmain,xcomp,ycomp,zcomp,wcomp)
#if KEY_EMAP==1
    else if(indxa(comlyn, comlen, 'EMAP')  >  0)then
       call emaprest(COMLYN,COMLEN,islct,jslct,lcomp,x,y,z,wmain,xcomp,ycomp,zcomp,wcomp)
#endif
    else

       IF(LCOMP) THEN
          CALL CSTRN2(COMLYN,COMLEN,XCOMP,YCOMP,ZCOMP,WCOMP,LCOMP, &
               ISLCT,JSLCT, &
               ICB,ICT,ICP,ICI, &
#if KEY_CMAP==1
               ICCT, &     
#endif
#if KEY_CPATH==1
               DX,DY,DZ, &    
#endif
               MSTUP,MUSTUP)
       ELSE
          CALL CSTRN2(COMLYN,COMLEN,X,Y,Z,WMAIN,LCOMP,ISLCT, &
               JSLCT,ICB,ICT,ICP,ICI, &
#if KEY_CMAP==1
               ICCT, &    
#endif
#if KEY_CPATH==1
               DX,DY,DZ, &  
#endif
               MSTUP,MUSTUP)
       ENDIF

    ENDIF !IF(INDXA(COMLYN, COMLEN, 'RMSD')  >  0)THEN

    call chmdealloc('cstran.src','CSTRAN','ispace',2*NATOM,intg=ispace)
    IF (MSTUP) CALL PSFSUM(OUTU)

    RETURN
  END SUBROUTINE CSTRAN

  SUBROUTINE CSTRN2(COMLYN,COMLEN,X,Y,Z,WT,QCOMP,ISLCT,JSLCT, &
       ICB,ICT,ICP,ICI, &
#if KEY_CMAP==1
       ICCT, & 
#endif
#if KEY_CPATH==1
       DX,DY,DZ, &       
#endif
       MSTUP,MUSTUP)
    !
    !             THIS ROUTINE SETS UP THE ARRAY IMOVE TO FIX ATOMS, OR
    !     THE NECESSARY ARRAYS TO APPLY HARMONIC OR DIHEDRAL CONSTRAINTS. IF
    !     IMOVE IS SET, THE MINIMIZATION AND DYNAMICS ROUTINES WILL FIX THE
    !     ATOM. THE ENERGY CALCULATION WILL EXCLUDE ALL CONTRIBUTIONS
    !     BETWEEN FIXED ATOMS. THE HARMONIC CONSTRAINT APPLYS A CONSTRAINT
    !     TO AN ATOM'S LOCATION IN SPACE, AND THE DIHEDRAL CONSTRAINT APPLYS
    !     AN IMPROPER DIHEDRAL POTENTIAL TO THE SPECIFIED TORSION.
    !
    !             THE SYNTAX FOR CONSTRAINTS IS AS FOLLOWS:
    !
    ! CONStraint FIX  purge-spec selection-spec
    !
    !       purge-spec ::= [PURGE] [BOND] [THETA] [PHI] [IMPHI] [CMAP]
    !
    ! CONS HARMonic { [ABSOlute] absolute-specs force-const-spec   coordinate-spec }
    !               {  BESTfit   bestfit-specs  force-const-spec   coordinate-spec }
    !               {  RELAtive  bestfit-specs  force-const-spec 2nd-atom-selection}
    !               {  CLEAr                                                       }
    !
    !      force-const-spec ::= { FORCE real } atom-selection [MASS]
    !                           { WEIGhting  }
    !
    !      absolute-specs ::= [EXPOnent int] [XSCAle real] [YSCAle real] [ZSCAle real]
    !
    !      bestfit-specs  ::= [ NOROtation ] [ NOTRanslation ]
    !
    !      coordinate-spec::= { [MAIN] }
    !                         {  COMP  }
    !                         {  KEEP  }
    !
    ! CONStraint DIHEdral [FORC real] [MIN real] dihedral-spec
    ! CONStraint CLDH
    !
    ! CONStraint IC       [BOND real] [EXPOnent int] [UPPEr]
    !                                           [ANGL real] [DIHE real]
    !
    ! CONStraint DROPlet  [FORCe real] [EXPOnent integer] [NOMAss]
    !
    ! CONStraint HMCM     [FORCe real] [WEIGhting] [RDISt real]
    !                     reference-spec coordinate-spec
    !            HMCR     [FORCe real] REFR real CM1 integer CM2 integer
    !            HMCM     CLEAr
    !
    !      reference-spec ::= REFX real REFY real REFZ real
    !
    ! CONStraint PATH     spline-pts [COMP]
    !                     FORCe real TZERo real selection-spec
    !                     FORCe real TZERo real selection-spec
    !                     ...
    !
    !            PATH     CLEAr
    !
    !            PATH     POST
    !
    !      spline-pts ::= { SPLIne selection-spec }
    !                     { SPLUnit unit }
    !
    !  where:
    !      selection-spec is an atom selection specification
    !
    !  The order of terms in the selection and dihedral specs may be
    !  important, but the  force, angle, and refrence units are picked up
    !  independent of order in the command line.
    !
    !  FIX is handled through setting IMOVE.  Each CONS FIX ... command line
    !  clears out the old set of fixed atoms and starts over.  The purge
    !  spec will allow some (by terms) or all (PURGE) of the terms in the
    !  energy expression that involve only fixed atoms to be deleted.
    !
    !  HARMonic constraints are applied through ECNSTR as isotropic harmonic
    !  constraints.  Each CONS HARM ... command line sets the constraints
    !  for the specified atoms without changing the constraints applied
    !  previously to other atoms.  To clear out the HARMonic constraints
    !  all the constraints can be set to zero with the command;
    !              CONS HARM FORCE 0.0 SELE ALL END
    !  If COMP is specified, the reference coordinates are copied from
    !  the comparison set.
    !
    !  For dihedrals
    !
    !       dihe-spec ::== [atom-spec  atom-spec  atom-spec atom-spec]
    !                      [BYNUM iatom jatom katom latom]
    !       atom-spec ::= SEGID RESID IUPAC
    !  The two real numbers are the force (in kcal/mole/radian**2)
    !  and minimum angle (in degrees) for dihedral constraint.  Each
    !  CONS DIHE ...  command line adds a new constraint dihedral without
    !  touching the old ones.
    !  CLDH clears the dihedral angle constraints.
    !
    !     First version written by Robert Bruccoleri
    !     Modifications by David States
    !      Overhauled by Bernard R. Brooks   1983
    !
    use chm_kinds
    use dimens_fcm
    use number
    use memory
    use psf
    use cnst_fcm
    use conshelix_m
    use consta
    use select
    use stream
    use string
    use fourdm
#if KEY_OPENMM==1
    use omm_ctrl, only : omm_system_changed
#endif
#if KEY_ACE==1
    use ace_module,only:qfiace           
#endif
#if KEY_CFF==1 || KEY_MMFF==1
    use cff_fcm
    use mmffm
    use ffieldm
#endif 
    use parallel
    use repdstr
    use intcor2,only:geticv
    use chutil
    use pucker_mod,only: setup_pucker_constraint,ncspuck
#if KEY_DOMDEC
    use domdec_common,only:domdec_system_changed
#endif
    use param_store, only: set_param

    use shake, only: qshake
    
    implicit none
    !
    INTEGER ICB(*),ICT(*),ICP(*),ICI(*)
#if KEY_CMAP==1
    INTEGER ICCT(*)                   
#endif
    INTEGER ISLCT(*),JSLCT(*)
    INTEGER COMLEN
    character(len=*) COMLYN
    real(chm_real) X(*),Y(*),Z(*),WT(*)
    LOGICAL QCOMP,MSTUP,MUSTUP
    !
    INTEGER I,J,IGRP,IS,IQ,QAT(6),NQAT,HRTYP,ISET,JSET
    real(chm_real) AMIN
    integer :: NINSETI(NUMHSETS), NINSETF(NUMHSETS+3)
    character(len=4) WRD
    character(len=8), dimension(6) :: SID, RID, REN, AC
    LOGICAL FIXED,LCMASS,QWEIG,QMAIN,QHMOD
    LOGICAL QBOND,QANGLE,QPHI,QIMP,ERR
    logical dihe_ok
#if KEY_CMAP==1
    LOGICAL QCMAP                   
#endif
    real(chm_real) KC,KV,RJ(4)
#if KEY_HMCOM==1
    real(chm_real) RT               
#endif
    character(len=4) CCOMP
#if KEY_CFF==1 || KEY_MMFF==1
    INTEGER K
#endif 
#if KEY_CPATH==1
    INTEGER SPLUNIT,TESTUNIT
    real(chm_real)  CPATHENER,DX(*),DY(*),DZ(*)
#endif 
#if KEY_CONSHELIX==1
    LOGICAL LCOM
#endif 
    integer,allocatable,dimension(:) :: mbnd,mtheta

    call allocate_cnst(natom)
    !
    !     THE FIRST WORD IN THE COMMAND LINE DETERMINES WHAT IS BEING DONE.
    !
#if KEY_MMFF==1 || KEY_CFF==1 /*ffield_fcm*/
    IF (FFIELD == CFF .OR. FFIELD == MMFF) THEN
       call chmalloc('cstran.src','CSTRAN','MBND',NBOND,intg=MBND)
       call chmalloc('cstran.src','CSTRAN','MTHETA',NTHETA,intg=MTHETA)
    else
       call chmalloc('cstran.src','CSTRAN','MBND',1,intg=MBND)
       call chmalloc('cstran.src','CSTRAN','MTHETA',1,intg=MTHETA)
    ENDIF
#else /* (ffield_fcm)*/
    call chmalloc('cstran.src','CSTRAN','MBND',1,intg=MBND)
    call chmalloc('cstran.src','CSTRAN','MTHETA',1,intg=MTHETA)
#endif /* (ffield_fcm)*/

    MSTUP=.FALSE.
    QHMOD=.FALSE.
    CCOMP='MAIN'
    IF(QCOMP) CCOMP='COMP'
    WRD=NEXTA4(COMLYN,COMLEN)
    !
    !     BRANCH TO THE APPROPRIATE SECTION.
    !
    !-----------------------------------------------------------------------
    IF(WRD == 'FIX ') THEN
       ! fix-atom-positions
#if KEY_ACE==1
       !        self energy contributions of fixed atoms will be recalculated
       QFIACE=0
#endif

       if (qshake) then
          call wrndie(-1, '<CSTRAN>', &
               'shake already active and may affect fixed atoms')
       end if

       ! First check the command line to see if internal coordinates should
       ! be purged for fixed atoms.
       !
       QBOND=INDXA(COMLYN,COMLEN,'BOND') /= 0
       QANGLE=INDXA(COMLYN,COMLEN,'THET') /= 0
       QPHI=INDXA(COMLYN,COMLEN,'PHI') /= 0
       QIMP=INDXA(COMLYN,COMLEN,'IMPH') /= 0
#if KEY_CMAP==1
       QCMAP=INDXA(COMLYN,COMLEN,'CMAP') /= 0
#endif 
       IF (INDXA(COMLYN,COMLEN,'PURG') /= 0) THEN
          QBOND=.TRUE.
          QANGLE=.TRUE.
          QPHI=.TRUE.
          QIMP=.TRUE.
#if KEY_CMAP==1
          QCMAP=.TRUE.
#endif 
       ENDIF
       IF (QBOND.OR.QANGLE.OR.QPHI.OR.QIMP &
#if KEY_CMAP==1
            .OR. QCMAP &        
#endif
            ) THEN
          IF(WRNLEV >= 2) WRITE(OUTU,1456)
1456      FORMAT(' CSTRAN: Irreversible compression of', &
               ' PSF is performed.')
          MSTUP=.TRUE.
       ENDIF
       !
       ! Now find the atoms that are to be fixed.
       !
       CALL SELCTA(COMLYN,COMLEN,ISLCT,X,Y,Z,WT,.TRUE.)
       !
       DO I=1,NATOM
          IF(IMOVE(I) >= 0) IMOVE(I)=ISLCT(I)
       ENDDO
#if KEY_OPENMM==1
       call omm_system_changed()  
#endif
       !
       ! Modify fixed internal coordinate lists (bonds,thetas,phis and imphis)
       !
#if KEY_CFF==1 || KEY_MMFF==1
       !      Initialize MBND and MTHETA for all bonds and angles fixed
       !
       IF (FFIELD == CFF .OR. FFIELD == MMFF) THEN
          DO I=1,NBOND
             MBND(I)=0
          ENDDO
          DO I=1,NTHETA
             MTHETA(I)=0
          ENDDO
       ENDIF
#endif 
       IF(QPHI) THEN
          J=0
          DO I=1,NPHI
             FIXED=IMOVE(IP(I)) > 0.AND.IMOVE(JP(I)) > 0.AND. &
                  IMOVE(KP(I)) > 0.AND.IMOVE(LP(I)) > 0
             IF(.NOT.FIXED) THEN
                J=J+1
                IP(J)=IP(I)
                JP(J)=JP(I)
                KP(J)=KP(I)
                LP(J)=LP(I)
                ICP(J)=ICP(I)
#if KEY_CFF==1
                !      Mark the two angles in this torsion as moving
                !
                IF (FFIELD == CFF) THEN
                   MTHETA(IPTW(1,I))=1
                   MTHETA(IPTW(2,I))=1
                ENDIF
#endif 
             ENDIF
          ENDDO
          NPHI=J
#if KEY_CFF==1
          IF (FFIELD == CFF) THEN
             J=0
             DO I=1,NBB
                K = IBBW(I)
                FIXED=IMOVE(IP(K)) > 0.AND.IMOVE(JP(K)) > 0.AND. &
                     IMOVE(KP(K)) > 0.AND.IMOVE(LP(K)) > 0
                IF (.NOT.FIXED) THEN
                   J=J+1
                   IBBW(J)=IBBW(I)
                ENDIF
             ENDDO
             NBB=J
          ENDIF
#endif 
       ENDIF
       IF(QIMP) THEN
          J=0
          DO I=1,NIMPHI
             FIXED=IMOVE(IM(I)) > 0.AND.IMOVE(JM(I)) > 0.AND. &
                  IMOVE(KM(I)) > 0.AND.IMOVE(LM(I)) > 0
             IF (.NOT.FIXED) THEN
                J=J+1
                IM(J)=IM(I)
                JM(J)=JM(I)
                KM(J)=KM(I)
                LM(J)=LM(I)
                ICI(J)=ICI(I)
             ENDIF
          ENDDO
          NIMPHI=J
       ENDIF

#if KEY_CMAP==1
       IF(QCMAP) THEN
          J=0
          DO I=1,NCRTERM
             FIXED=IMOVE(I1CT(I)) > 0.AND. &
                  IMOVE(J1CT(I)) > 0.AND. &
                  IMOVE(K1CT(I)) > 0.AND. &
                  IMOVE(L1CT(I)) > 0.AND. &
                  IMOVE(I2CT(I)) > 0.AND. &
                  IMOVE(J2CT(I)) > 0.AND. &
                  IMOVE(K2CT(I)) > 0.AND. &
                  IMOVE(L2CT(I)) > 0
             IF (.NOT.FIXED) THEN
                J=J+1
                I1CT(J)=I1CT(I)
                J1CT(J)=J1CT(I)
                K1CT(J)=K1CT(I)
                L1CT(J)=L1CT(I)
                I2CT(J)=I2CT(I)
                J2CT(J)=J2CT(I)
                K2CT(J)=K2CT(I)
                L2CT(J)=L2CT(I)
                ICCT(J)=ICCT(I)
             ENDIF
          ENDDO
          NCRTERM=J
       ENDIF
#endif 

       IF(QANGLE) THEN
#if KEY_CFF==1
          !      Loop through the angle-angle cross-terms.  If any atom in a term
          !      is moving then mark both angles as moving.  Otherwise, remove it.
          !
          IF (FFIELD == CFF) THEN
             J=0
             DO I=1,NTHTH
                FIXED=IMOVE(ITTW(1,I)) > 0.AND.IMOVE(ITTW(2,I)) > 0.AND. &
                     IMOVE(ITTW(3,I)) > 0.AND.IMOVE(ITTW(4,I)) > 0
                IF (.NOT.FIXED) THEN
                   J=J+1
                   ITTW(1,J)=ITTW(1,I)
                   ITTW(2,J)=ITTW(2,I)
                   ITTW(3,J)=ITTW(3,I)
                   ITTW(4,J)=ITTW(4,I)
                   ICTT(J)=ICTT(I)
                   MTHETA(ITTTW(1,I))=1
                   MTHETA(ITTTW(2,I))=1
                ENDIF
             ENDDO
             NTHTH=J
          ENDIF
#endif 
#if KEY_CFF==1 || KEY_MMFF==1
          !      Loop through all the angles.  If they are still marked as fixed
          !      then check to see if any atoms are moving.  This will catch all
          !      moving angles for MMFF.  For CFF it will catch 3 atom molecules
          !      which aren't part of a cross-term.
          !
          IF (FFIELD == CFF .OR. FFIELD == MMFF) THEN
             DO I=1,NTHETA
                IF (MTHETA(I) == 0) THEN
                   FIXED=IMOVE(IT(I)) > 0.AND.IMOVE(JT(I)) > 0.AND. &
                        IMOVE(KT(I)) > 0
                   IF (.NOT.FIXED) MTHETA(I)=1
                ENDIF
             ENDDO
             !      Now remove all fixed angles and mark the bonds in the moving
             !      angles as moving.
             J=0
             DO I=1,NTHETA
                IF (MTHETA(I) == 1) THEN
                   J=J+1
                   IT(J)=IT(I)
                   JT(J)=JT(I)
                   KT(J)=KT(I)
                   ICT(J)=ICT(I)
#if KEY_CFF==1
                   IF (FFIELD == CFF) THEN
                      MBND(ITBW(1,I))=1
                      MBND(ITBW(2,I))=1
                   ENDIF
#endif 
#if KEY_MMFF==1
                   IF (FFIELD == MMFF) THEN
                      MBND(STRBLIST(I,1))=1
                      MBND(STRBLIST(I,2))=1
                   ENDIF
#endif 
                   MTHETA(I)=J
                ENDIF
             ENDDO
          ELSE
#endif 
             J=0
             DO I=1,NTHETA
                FIXED=IMOVE(IT(I)) > 0.AND.IMOVE(JT(I)) > 0.AND. &
                     IMOVE(KT(I)) > 0
                IF(.NOT.FIXED) THEN
                   J=J+1
                   IT(J)=IT(I)
                   JT(J)=JT(I)
                   KT(J)=KT(I)
                   ICT(J)=ICT(I)
                ENDIF
             ENDDO
#if KEY_CFF==1 || KEY_MMFF==1
          ENDIF
#endif 
          NTHETA=J
       ENDIF
       IF(QBOND) THEN
#if KEY_CFF==1 || KEY_MMFF==1
          !      Loop through all the bonds and if they aren't already marked as
          !      moving then check whether either atom is moving.  This will catch
          !      any 2 atom molecules which aren't part of an angle.
          !
          IF (FFIELD == CFF .OR. FFIELD == MMFF) THEN
             DO I=1,NBOND
                IF (MBND(I) == 0) THEN
                   FIXED=IMOVE(IB(I)) > 0.AND.IMOVE(JB(I)) > 0
                   IF (.NOT.FIXED) MBND(I)=1
                ENDIF
             ENDDO
             !      Now remove the fixed bonds.
             J=0
             DO I=1,NBOND
                IF (MBND(I) == 1) THEN
                   J=J+1
                   IB(J)=IB(I)
                   JB(J)=JB(I)
                   ICB(J)=ICB(I)
                   MBND(I)=J
                ENDIF
             ENDDO
             NBOND=J
#if KEY_MMFF==1
             !      Since bond numbers have changed it is also necessary
             !      to change the pointers to the bonds in each angle.
             IF (FFIELD == MMFF) THEN
                DO I=1,NTHETA
                   STRBLIST(I,1)=MBND(STRBLIST(I,1))
                   STRBLIST(I,2)=MBND(STRBLIST(I,2))
                ENDDO
             ENDIF
#endif 
#if KEY_CFF==1
             !      Since bond numbers have changed it is also necessary
             !      to change the pointers to the bonds in each angle.
             IF (FFIELD == CFF) THEN
                DO I=1,NTHETA
                   ITBW(1,I)=MBND(ITBW(1,I))
                   ITBW(2,I)=MBND(ITBW(2,I))
                ENDDO
                !       Change the pointers to the bonds and angles in each torsion.
                DO I=1,NPHI
                   IPBW(1,I)=MBND(IPBW(1,I))
                   IPBW(2,I)=MBND(IPBW(2,I))
                   IPBW(3,I)=MBND(IPBW(3,I))
                   IPTW(1,I)=MTHETA(IPTW(1,I))
                   IPTW(2,I)=MTHETA(IPTW(2,I))
                ENDDO
                !      Change the pointers to the angles in each angle-angle term.
                DO I=1,NTHTH
                   ITTTW(1,I)=MTHETA(ITTTW(1,I))
                   ITTTW(2,I)=MTHETA(ITTTW(2,I))
                ENDDO
             ENDIF
#endif 
          ELSE
#endif 
             DO I=1,NBOND
                FIXED=IMOVE(IB(I)) > 0.AND.IMOVE(JB(I)) > 0
                IF (FIXED) ICB(I)=0
             ENDDO
#if KEY_CFF==1 || KEY_MMFF==1
          ENDIF
#endif 
       ENDIF
       !
       ! Set up group fixed array
       DO IGRP=1,NGRP
          IS=IGPBS(IGRP)+1
          IQ=IGPBS(IGRP+1)
          J=IMOVE(IS)
          DO I=IS,IQ
             IF(IMOVE(I) <= 0) J=0
          ENDDO
          IMOVEG(IGRP)=J
       ENDDO
       !
       MUSTUP=.TRUE.
       !
       !     REMOVE/RESTORE CONSTRAINT DIHEDRALS FOR FIXED/FREED ATOMS
       !
       !     CALL EPHI(ETERM(CDIHE),ICS,JCS,KCS,LCS,ICCS,NCSPHI,
       DO I=1,NCSPHI
          FIXED=IMOVE(ICS(I)) > 0.AND.IMOVE(JCS(I)) > 0
          FIXED=FIXED.AND.IMOVE(KCS(I)) > 0.AND.IMOVE(LCS(I)) > 0
          IF (FIXED) THEN
             ICCS(I)=0
          ELSE
             ICCS(I)=I
          ENDIF
       ENDDO
       !
       !-----------------------------------------------------------------------
    ELSE IF(WRD == 'FIX4') THEN
       !
#if KEY_FOURD==1
       ! Now find the atoms that are to be fixed in 4D.
       !
       CALL SELCTA(COMLYN,COMLEN,ISLCT,X,Y,Z,WT,.TRUE.)
       !
       DO I=1,NATOM
          IF(IMOVE4(I) >= 0) IMOVE4(I)=ISLCT(I)
       ENDDO
#else /**/
       CALL WRNDIE(-1,'<CSTRAN>','No 4D code compiled')
       call cstrn2_RETURN
       return
#endif 
       !-----------------------------------------------------------------------
       ! CONStraint HARMonic { [ABSOlute]  absolute-specs    } force-const-spec
       !                     {  BESTfit    coordinate-spec   }
       !                     {  RELAtive  2nd-atom-selection }
       !                     {  CLEAr                        }
       !
       !       force-const-spec ::=  { FORCE real } atom-selection [MASS]
       !                             { WEIGhting  }
       !
       !       absolute-specs ::= coordinate-spec [EXPOnent int]
       !                              [XSCAle real] [YSCAle real] [ZSCAle real]
       !
       !       coordinate-spec::= { [MAIN] }
       !                          {  COMP  }
       !                          {  KEEP  }
       !
    ELSE IF(WRD == 'HARM') THEN
       ! set-up-harmonic-constraints
       !
       ! The harmonic position restraints are parsed here.
       !
       ! Determine what type of restraint set is requested.
       !
       WRD=NEXTA4(COMLYN,COMLEN)
       IF(WRD == 'CLEA') THEN
          ! Clear the harmonic restraints
          QCNSTR=.FALSE.
          NUMHSETS=0
          call hset_clear()
          DO I=1,NATOM
             IHSET(I)=0
             KCNSTR(I)=0.0
             REFX(I)=0.0
             REFY(I)=0.0
             REFZ(I)=0.0
          ENDDO
          call cstrn2_RETURN
          return
       ELSE IF(WRD == 'BEST') THEN
          ! process the best-fit option (add a new set)
          HRTYP=1
          QHMOD=.TRUE.
       ELSE IF(WRD == 'RELA') THEN
          ! process the relative option (add two new sets)
          HRTYP=2
       ELSE IF(WRD == 'ABSO') THEN
          ! The old absolute approach
          QHMOD=.TRUE.
          HRTYP=0
       ELSE
          ! Assume the old approach (put the command back)
          CALL JOINWD(COMLYN,MXCMSZ,COMLEN,WRD,4)
          HRTYP=0
          QHMOD=.TRUE.
       ENDIF
       !
       ! make sure the data is cleared if disabled.
       IF(.NOT.QCNSTR) THEN
          NUMHSETS=0
          call hset_clear()
          DO I=1,NATOM
             KCNSTR(I)=0.0
             IHSET(I)=0
          ENDDO
       ENDIF
       !
       ! Determine active sets
       NINSETI = 0
       NINSETF = 0
       DO I=1,NATOM
          J=IHSET(I)
          IF(J > NUMHSETS .OR. J < 0) THEN
             CALL WRNDIE(-4,'<CSTRAN>', &
                  'Bad atom set value: coding error')
             call cstrn2_RETURN
             return
          ENDIF
          IF(J > 0) NINSETI(J)=NINSETI(J)+1
       ENDDO
       !
       ! parse the atom selection(s)
       ERR=.FALSE.
       IF(HRTYP == 2) THEN
          CALL SELCTD(COMLYN,COMLEN,ISLCT,JSLCT,X,Y,Z,WT,.TRUE.,ERR)
          IF(ERR) THEN
             CALL WRNDIE(-1,'<CSTRAN>', &
                  'Double atom selection error - IGNORED')
             call cstrn2_RETURN
             return
          ENDIF
       ELSE
          CALL SELCTA(COMLYN,COMLEN,ISLCT,X,Y,Z,WT,.TRUE.)
       ENDIF
       !
       ! Check to see if any active sets will be disrupted
       !
       DO I=1,NATOM
          J=IHSET(I)
          IF(ISLCT(I) == 1) THEN
             J=-1
             IHSET(I)=-1
          ENDIF
          IF(HRTYP == 2) THEN
             IF(JSLCT(I) == 1) THEN
                IF(J < 0) ERR=.TRUE.
                IHSET(I)=-2
             ENDIF
          ENDIF
          IF(J > 0) NINSETF(J)=NINSETF(J)+1
       ENDDO
       IF(ERR) THEN
          CALL WRNDIE(-1,'<CSTRAN>', &
               'Double atom selection overlap - IGNORED')
          call cstrn2_RETURN
          return
       ENDIF
       DO I=1,NUMHSETS
          IF(NINSETI(I) /= NINSETF(I) .AND. NINSETF(I) > 0) THEN
             !           this set is disrupted - issue error message
             IF(WRNLEV > 2) WRITE(OUTU,235) I,NINSETI(I),NINSETF(I)
235          FORMAT(' *****  WARNING  ***** ',/, &
                  ' CSTRAN: Harmonic restraint set',I3,' is disrupted.',/ &
                  '     Previous atom count:',I7,'  New atom count:',I7,/ &
                  '     Some atoms have been removed from this set.')
             CALL WRNDIE(-1,'<CSTRAN>', &
                  'Harmonic restraint atom selection overlap')
          ENDIF
       ENDDO
       !
       ! Get a new set number (look for an empty set to use)
       !
       ISET=-1
       J=1
       IF(HRTYP == 2) THEN
          J=2
          IF(NUMHSETS == 0) THEN
             NUMHSETS=1
             call hset_ensure_capacity(NUMHSETS)
             KCEXPN(1)=2
             XHSCALE(1)=ONE
             YHSCALE(1)=ONE
             ZHSCALE(1)=ONE
             TYPHSET(1)=0
             QHNORT(1)=.FALSE.
             QHNOTR(1)=.FALSE.
          ENDIF
       ENDIF
       DO I=NUMHSETS,J,-1
          IF(NINSETF(I) == 0) ISET=I
       ENDDO
       IF(ISET < 0) THEN
          NUMHSETS=NUMHSETS+1
          call hset_ensure_capacity(NUMHSETS)
          ISET=NUMHSETS
       ENDIF
       !
       IF(HRTYP == 0) THEN
          KCEXPN(ISET)=GTRMI(COMLYN,COMLEN,'EXPO',2)
          XHSCALE(ISET)=GTRMF(COMLYN,COMLEN,'XSCA',ONE)
          YHSCALE(ISET)=GTRMF(COMLYN,COMLEN,'YSCA',ONE)
          ZHSCALE(ISET)=GTRMF(COMLYN,COMLEN,'ZSCA',ONE)
          QHNORT(ISET)=.TRUE.
          QHNOTR(ISET)=.TRUE.
#if KEY_OPENMM==1
          call omm_system_changed()  
#endif
       ELSE
          KCEXPN(ISET)=2
          XHSCALE(ISET)=ONE
          YHSCALE(ISET)=ONE
          ZHSCALE(ISET)=ONE
          QHNORT(ISET)=INDXA(COMLYN,COMLEN,'NORO') > 0
          QHNOTR(ISET)=INDXA(COMLYN,COMLEN,'NOTR') > 0
       ENDIF
       TYPHSET(ISET)=HRTYP
       !
       ! Get a second set number (look for an empty set to use)
       !
       IF(HRTYP == 2) THEN
          JSET=-1
          DO I=NUMHSETS,2,-1
             IF(NINSETF(I) == 0 .AND. I /= ISET) JSET=I
          ENDDO
          IF(JSET < 0) THEN
             NUMHSETS=NUMHSETS+1
             call hset_ensure_capacity(NUMHSETS)
             JSET=NUMHSETS
          ENDIF
          !
          KCEXPN(JSET) = KCEXPN(ISET)
          XHSCALE(JSET)= XHSCALE(ISET)
          YHSCALE(JSET)= YHSCALE(ISET)
          ZHSCALE(JSET)= ZHSCALE(ISET)
          QHNORT(JSET) = QHNORT(ISET)
          QHNOTR(JSET) = QHNOTR(ISET)
          TYPHSET(JSET)= ISET
          TYPHSET(ISET)= JSET
          NINSETF(JSET)=0
       ENDIF
       !
       NINSETF(ISET)=0
       DO I=1,NATOM
          IF(IHSET(I) == -1) THEN
             IHSET(I)=ISET
             NINSETF(ISET)=NINSETF(ISET)+1
          ELSE IF(IHSET(I) == -2) THEN
             IHSET(I)=JSET
             NINSETF(JSET)=NINSETF(JSET)+1
          ENDIF
       ENDDO
       IF(HRTYP == 2) THEN
          IF(NINSETF(ISET) /= NINSETF(JSET)) CALL WRNDIE(-1,'<CSTRAN>', &
               'Number of selected atoms does not match - Extras ignored')
          NINSETF(JSET)=0
       ENDIF
       !
       ! Parse the force constant specs
       !
       LCMASS=INDXA(COMLYN,COMLEN,'MASS') > 0
       QWEIG=(INDXA(COMLYN,COMLEN,'WEIG') > 0)
       IF(.NOT.QWEIG) THEN
          KC=GTRMF(COMLYN,COMLEN,'FORC',MINSIX)
          IF(KC == MINSIX) THEN
             CALL WRNDIE(2,'<CSTRN2>', &
                  'No force constant specified. Set to zero')
             KC=ZERO
          ENDIF
          IF(KC < 0.0) THEN
             CALL WRNDIE(2,'<CSTRN2>', &
                  'Negative force constant specified. Set to zero')
             KC=ZERO
          ENDIF
       ENDIF
       !
       DO I=1,NATOM
          IF(ISLCT(I) == 1) THEN
             IF(QWEIG) THEN
                KV=WT(I)
             ELSE
                KV=KC
             ENDIF
             IF(LCMASS) KV=KV*AMASS(I)
             KCNSTR(I)=KV
          ENDIF
       ENDDO
       !
       IF(QHMOD) THEN
          IF(.NOT.QCOMP) THEN
             QMAIN=(INDXA(COMLYN,COMLEN,'MAIN') > 0)
             IF(.NOT.QMAIN) QHMOD=(INDXA(COMLYN,COMLEN,'KEEP') == 0)
          ENDIF
       ENDIF
       !
       IF(QHMOD) THEN
          DO I=1,NATOM
             IF(ISLCT(I) == 1) THEN
                REFX(I)=X(I)
                REFY(I)=Y(I)
                REFZ(I)=Z(I)
             ENDIF
          ENDDO
          IF(INDXA(COMLYN,COMLEN,'KEEP') > 0) CALL WRNDIE(0, &
               '<CSTRAN>','KEEP option not valid')
       ELSE
          IF(QCOMP) CALL WRNDIE(0,'<CSTRAN>', &
               'COMP option not valid - KEEP used')
          IF(INDXA(COMLYN,COMLEN,'MAIN') > 0) CALL WRNDIE(0, &
               '<CSTRAN>','MAIN option not valid - KEEP used')
       ENDIF
       !
       ! Print out the status of the new harmonic restraints
       !
       IF(PRNLEV > 2) THEN
          WRITE(OUTU,1204) 'CSTRAN: Harmonic Restraints'
          !
          IF(HRTYP == 0) THEN
             WRITE(OUTU,1206) &
                  'ABSOlute type as set number',ISET, &
                  '.  Number of selected atoms:',NINSETF(ISET)
          ELSE IF(HRTYP == 1) THEN
             WRITE(OUTU,1206) &
                  'BESTfit type as set number',ISET, &
                  '.  Number of selected atoms:',NINSETF(ISET)
          ELSE IF(HRTYP == 2) THEN
             WRITE(OUTU,1206) &
                  'RELAtive type as set number',ISET, &
                  '.  Number of selected atoms:',NINSETF(ISET)
             WRITE(OUTU,1206) &
                  'bestfit to atoms set number',JSET, &
                  '.  Number of selected atoms:',NINSETF(JSET)
          ENDIF
          !
          IF (QHMOD) THEN
             IF(QCOMP) THEN
                WRITE(OUTU, 1205) &
                     'Reference coordinates set to comparison coordinates.'
             ELSE
                WRITE(OUTU, 1205) &
                     'Reference coordinates set to main coordinates.'
             ENDIF
          ELSE
             WRITE(OUTU, 1205) &
                  'No reference coordinates used.'
          ENDIF
          !
          IF(LCMASS) THEN
             WRITE(OUTU,1205) &
                  'Mass weighting will be used for new restraints.'
          ELSE
             WRITE(OUTU,1205) &
                  'Mass weighting will NOT be used for new restraints.'
          ENDIF
          IF(QWEIG) THEN
             WRITE(OUTU,1205) &
                  'The weighting array will be used for force constants.'
          ELSE
             WRITE(OUTU,1207) 'The force constant of',KC, &
                  ' will be used.'
          ENDIF
          !
          IF(HRTYP == 0) THEN
             WRITE(OUTU,1206) 'An exponent of',KCEXPN(ISET), &
                  ' will be used.'
             WRITE(OUTU,1208) XHSCALE(ISET),YHSCALE(ISET),ZHSCALE(ISET)
          ELSE
             IF(QHNORT(ISET)) WRITE(OUTU,1205) &
                  'The bestfit will be done with no rotation.'
             IF(QHNOTR(ISET)) WRITE(OUTU,1205) &
                  'The bestfit will be done with no translation.'
          ENDIF
          !
          J=0
          DO I=1,NATOM
             IF(KCNSTR(I) /= ZERO .AND. IHSET(I) > 0) J=J+1
          ENDDO
          WRITE(OUTU,1209) 'A total of',J,' atoms are restrained.'
1204      FORMAT(1X,A)
1205      FORMAT(10X,A)
1206      FORMAT(10X,A,I3,A,I7)
1207      FORMAT(10X,A,F14.5,A)
1208      FORMAT(10X,'The XYZ scale factors are:',3F14.5)
1209      FORMAT(10X,A,I7,A,I7)


          ! Set the QCNSTR flag to indicate the presence of constraints if
          ! any have been set.
          !
       ENDIF
       !
       QCNSTR=.FALSE.
       DO I=1,NATOM
          IF(KCNSTR(I) /= 0.0 .AND. IHSET(I) > 0) THEN
             QCNSTR=.TRUE.
             call cstrn2_RETURN
             return
          ENDIF
       ENDDO
       call cstrn2_RETURN
       return
       !-----------------------------------------------------------------------
    ELSE IF(WRD == 'DIHE') THEN
       ! Add a dihedral restraint and check it
       dihe_ok = .false.
       NCSPHI=NCSPHI+1
       CALL set_param('NCSP',NCSPHI)
       call csphi_ensure_capacity(NCSPHI)
       ICCS(NCSPHI)=NCSPHI
       CCSC(NCSPHI)=GTRMF(COMLYN,COMLEN,'FORC',ZERO)
#if KEY_REPDSTR==1
       IF (.not. QREPDSTR) THEN  
#endif
          IF(CCSC(NCSPHI) == 0.0) THEN
             IF(WRNLEV >= 2) WRITE(OUTU,401)
401          FORMAT(' *****  WARNING  ***** ',/, &
                  ' Constraint dihedral has no force', &
                  ' constant specified.')
             GOTO 800
          ENDIF
#if KEY_REPDSTR==1
       ENDIF  
#endif
       AMIN=GTRMF(COMLYN,COMLEN,'MIN',FMARK)
       IF(AMIN /= FMARK) THEN
          CCSB(NCSPHI)=AMIN*PI/180.0
       ELSE
          QMAIN=(INDXA(COMLYN,COMLEN,'MAIN') > 0)
          IF(QMAIN.OR.QCOMP) THEN
             IF(PRNLEV > 2) WRITE(OUTU,416) CCOMP
416          FORMAT(' CSTRAN: Restrained dihedral minimum angle value', &
                  ' is obtained from the ',A4,' coordinates.')
          ELSE
             IF(WRNLEV >= 2) WRITE(OUTU,402)
402          FORMAT(' *****  ERROR  ***** ',/,' Angle not', &
                  ' specified for a constraint dihedral.')
             GOTO 800
          ENDIF
       ENDIF
       ! Choice of the functional form. AB.
       CCSD(NCSPHI)=GTRMI(COMLYN,COMLEN,'PERI',-9999)
       CCSW(NCSPHI)=GTRMF(COMLYN,COMLEN,'WIDT',ZERO)
       !
       IF (CCSD(NCSPHI) == -9999) THEN
          IF(WRNLEV > 2) WRITE(OUTU,403)
403       FORMAT(' No periodicity specified: Econs=K*(phi-phi0)**2', &
               ' will be used to constraint.')
          CCSD(NCSPHI)=0
       ENDIF
       IF (CCSD(NCSPHI) < 0) THEN
          CCSD(NCSPHI)=-CCSD(NCSPHI)
          IF(WRNLEV > 2) WRITE(OUTU,404) CCSD(NCSPHI)
404       FORMAT(' Multiple constraint dihedral not allowed.', &
               '  Periodicity set to :',I6)
       ENDIF
       !
       IF(CCSW(NCSPHI) < ZERO) THEN
          IF(WRNLEV > 2) WRITE(OUTU,405)
405       FORMAT(' Negative widths not allowed.')
          CCSW(NCSPHI)=ZERO
       ENDIF
       !
       IF(CCSW(NCSPHI) > ZERO .AND. CCSD(NCSPHI) > 0) THEN
          IF(WRNLEV > 2) WRITE(OUTU,406)
406       FORMAT(' Width specification allowed only for improper type.')
          CCSW(NCSPHI)=ZERO
       ENDIF
       !
       CALL NXTATM(QAT,NQAT,4,COMLYN,COMLEN,ISLCT, &
            SEGID,RESID,ATYPE,IBASE,NICTOT,NSEG,RES,NATOM)
       IF(NQAT /= 4) GOTO 800
       ICS(NCSPHI)=QAT(1)
       I=ICS(NCSPHI)
       IF(I <= 0.OR.I > NATOMT) GOTO 800
       CALL ATOMID(I,SID(1),RID(1),REN(1),AC(1))
       JCS(NCSPHI)=QAT(2)
       I=JCS(NCSPHI)
       IF(I <= 0.OR.I > NATOMT) GOTO 800
       CALL ATOMID(I,SID(2),RID(2),REN(2),AC(2))
       KCS(NCSPHI)=QAT(3)
       I=KCS(NCSPHI)
       IF(I <= 0.OR.I > NATOMT) GOTO 800
       CALL ATOMID(I,SID(3),RID(3),REN(3),AC(3))
       LCS(NCSPHI)=QAT(4)
       I=LCS(NCSPHI)
       IF(I <= 0.OR.I > NATOMT) GOTO 800
       CALL ATOMID(I,SID(4),RID(4),REN(4),AC(4))
       dihe_ok = .true.

800    if (dihe_ok) then
          ! Dihedral restraint OK, keep it
          IF(AMIN == FMARK) THEN
             CALL GETICV(QAT(1),QAT(2),QAT(3),QAT(4),.FALSE., &
                  RJ(1),RJ(2),AMIN,RJ(3),RJ(4),X,Y,Z)
             CCSB(NCSPHI)=AMIN*PI/180.0
          ENDIF

          IF (CCSD(NCSPHI) == 0) THEN
             CCSCOS(NCSPHI)=COS(CCSB(NCSPHI))
             CCSSIN(NCSPHI)=SIN(CCSB(NCSPHI))
          ELSE
             CCSCOS(NCSPHI)=COS(CCSB(NCSPHI)*CCSD(NCSPHI)+PI)
             CCSSIN(NCSPHI)=SIN(CCSB(NCSPHI)*CCSD(NCSPHI)+PI)
          ENDIF

          IF(PRNLEV >= 2) WRITE(OUTU,425) NCSPHI, &
               RID(1)(1:idleng),SID(1)(1:idleng),AC(1)(1:idleng), &
               RID(2)(1:idleng),SID(2)(1:idleng),AC(2)(1:idleng), &
               RID(3)(1:idleng),SID(3)(1:idleng),AC(3)(1:idleng), &
               RID(4)(1:idleng),SID(4)(1:idleng),AC(4)(1:idleng), &
               CCSC(NCSPHI),AMIN,CCSD(NCSPHI),CCSW(NCSPHI)
425       FORMAT(' NEW DIHEDRAL RESTRAINT ADDED'/I5,':',4X, &
               3(A,1X,A,1X,A,'/  '),(A,1X,A,1X,A),/, &
               6X,'FORCe=',F9.2,'  MIN=',F9.2,'  PERIod=',I4, &
               '  WIDTh=',F9.2)
#if KEY_DOMDEC==1
          call domdec_system_changed()
#endif
#if KEY_OPENMM==1
          call omm_system_changed()  
#endif
       else
          ! Dihedral restraint bad, remove it
          NCSPHI=NCSPHI-1
          CALL set_param('NCSP',NCSPHI)
          IF(WRNLEV >= 2) CALL WRNDIE(0,'<CSTRAN>', &
               'Some component of the dihedral restraint was invalid')
       endif
       !-----------------------------------------------------------------------
    else if (wrd == 'PUCK') then
       call setup_pucker_constraint(comlyn,comlen)
       !-----------------------------------------------------------------------
    ELSE IF(WRD == 'CLPC') THEN
       ! Clear all puckering restraints
       ncspuck = 0
       call set_param('ncsu',ncspuck)
       call csphi_clear()

       !-----------------------------------------------------------------------
    ELSE IF(WRD == 'CLDH') THEN
       ! Clear all dihedral restraints
       NCSPHI=0
       CALL set_param('NCSP',NCSPHI)
       call csphi_clear()
#if KEY_DOMDEC
       call domdec_system_changed()
#endif
#if KEY_OPENMM==1
       call omm_system_changed()  
#endif
       !-----------------------------------------------------------------------
    ELSE IF(WRD == 'IC  ') THEN
       ! set-up-ic-constraints
       LCIC=.TRUE.
       CCBIC=GTRMF(COMLYN,COMLEN,'BOND',ZERO)
       CCTIC=GTRMF(COMLYN,COMLEN,'ANGL',ZERO)
       CCPIC=GTRMF(COMLYN,COMLEN,'DIHE',ZERO)
       CCIIC=GTRMF(COMLYN,COMLEN,'IMPR',ZERO)
       KBEXPN=GTRMI(COMLYN,COMLEN,'EXPO',2)
       LUPPER=(INDXA(COMLYN,COMLEN,'UPPE') > 0)
       IF(MOD(KBEXPN,2) /= 0) KBEXPN=KBEXPN+1
       IF(KBEXPN < 2) KBEXPN=2
       IF(CCBIC == ZERO .AND. CCTIC == ZERO .AND. &
            CCPIC == ZERO .AND. CCIIC == ZERO) LCIC=.FALSE.
       IF(PRNLEV >= 2) THEN
          IF(LCIC) THEN
             WRITE(OUTU,132) CCBIC,CCTIC,CCPIC,CCIIC,KBEXPN,LUPPER
132          FORMAT(' INTERNAL COORDINATE CONSTRAINTS ACTIVE.', &
                  ' FORCE CONSTANTS:'/, &
                  '   BONDS =',F14.6,'  ANGLES =',F14.6, &
                  '  DIHEDRALS =',F14.6,/'  IMPROPERS =',F14.6, &
                  '   EXPON =',I3,13X, &
                  '(BONDS AS UPPER LIMIT: ',L6,')')
          ELSE
             WRITE(OUTU,134)
134          FORMAT(' INTERNAL COORDINATE CONSTRAINTS TURNED OFF')
          ENDIF
       ENDIF
       !-----------------------------------------------------------------------
    ELSE IF(WRD == 'DROP') THEN
       ! set-up-droplet-constraints
       KQCNST=GTRMF(COMLYN,COMLEN,'FORC',ZERO)
       KQEXPN=GTRMI(COMLYN,COMLEN,'EXPO',4)
       !        LQMASS=INDXA(COMLYN,COMLEN,'MASS') > 0
       LQMASS=.NOT.(INDXA(COMLYN,COMLEN,'NOMA') > 0)
       QQCNST=(KQCNST > 0.0)
       IF(PRNLEV >= 2) THEN
          IF(QQCNST) THEN
             IF(LQMASS) THEN
                WRITE(OUTU,142) KQCNST,KQEXPN
             ELSE
                WRITE(OUTU,143) KQCNST,KQEXPN
             ENDIF
          ELSE
             WRITE(OUTU,144)
          ENDIF
142       FORMAT(' DROPLET CONTRAINT ACTIVE WITH MASS WEIGHTING.'/, &
               '   FORCE CONSTANT =',E14.4,'  EXPONENT =',I4)
143       FORMAT(' DROPLET CONTRAINT ACTIVE. NO MASS WEIGHTING.'/, &
               '   FORCE CONSTANT =',E14.4,'  EXPONENT =',I4)
144       FORMAT(' DROPLET CONSTRAINTS REMOVED')
       ENDIF
#if KEY_HMCOM==1
       ! ---------------------------------------------------------------------
    else if(wrd == 'HMCM') then
       ! set up center of mass harmonic constraints
       if (indxa(comlyn,comlen,'CLEA') > 0) then
          nphmcm=0
          nhmcm=0
          qhmcm=.false.
          call deallocate_hmcm
#if KEY_DOMDEC==1
          call domdec_system_changed()
#endif
       else
          nhmcm=nhmcm+1
          if (nhmcm >= maxhcm) call allocate_hmcm
          khmcm(nhmcm)=gtrmf(comlyn,comlen,'FORC',zero)
          rhmcm(nhmcm)=gtrmf(comlyn,comlen,'RDIS',zero)
          lhmcmm(nhmcm)= ( indxa(comlyn,comlen,'WEIG') > 0 ) &
               .or. ( indxa(comlyn,comlen,'MASS') > 0 ) 
          hmcmx(nhmcm)=gtrmf(comlyn,comlen,'REFX',fmark)
          hmcmy(nhmcm)=gtrmf(comlyn,comlen,'REFY',fmark)
          hmcmz(nhmcm)=gtrmf(comlyn,comlen,'REFZ',fmark)
          hmcmtype(nhmcm)=0
          if(hmcmx(nhmcm) /= fmark) hmcmtype(nhmcm)=hmcmtype(nhmcm)+1
          if(hmcmy(nhmcm) /= fmark) hmcmtype(nhmcm)=hmcmtype(nhmcm)+2
          if(hmcmz(nhmcm) /= fmark) hmcmtype(nhmcm)=hmcmtype(nhmcm)+4  ! why is this 4?
          
          
          if (hmcmtype(nhmcm) == 0 ) then
             if (wrnlev >= 2) call wrndie(0,'<cstran>', &
                  'At least one X,Y, or Z value of reference center of mass must ' &
                  //'be specified.')
             write(outu,'(" HMCM>>>> last restraint will be ignored")')
             nhmcm=nhmcm-1
          else
             ihmcm(nhmcm)=nphmcm+1
             call selcta(comlyn,comlen,islct,x,y,z,wt,.true.)
             inhmcm(nhmcm)=0
             do i=1,natom
                if (islct(i) > 0) then
                   nphmcm=nphmcm+1
                   inhmcm(nhmcm)=inhmcm(nhmcm)+1
                endif
             enddo
             if(.not.allocated(phmcm) .or. nphmcm > maxpcm) call allocate_hmcm
             nphmcm=nphmcm-inhmcm(nhmcm)
             do i=1,natom
                if (islct(i) > 0) then
                   nphmcm=nphmcm+1
                   phmcm(nphmcm)=i
                endif
             enddo
             if(inhmcm(nhmcm) > 0) then
                qhmcm=.true.
                if (prnlev >= 2) then
                   if(lhmcmm(nhmcm)) then
                      write(outu,191) 'MASS', &
                           nhmcm,khmcm(nhmcm),rhmcm(nhmcm), &
                           hmcmx(nhmcm),hmcmy(nhmcm),hmcmz(nhmcm), &
                           inhmcm(nhmcm),nphmcm
                   else
                      write(outu,191) 'Geometry', &
                           nhmcm,khmcm(nhmcm),rhmcm(nhmcm), &
                           hmcmx(nhmcm),hmcmy(nhmcm),hmcmz(nhmcm), &
                           inhmcm(nhmcm),nphmcm
                   endif
191                format(' Center of ',a,' Harmonic Constraint #',i3,/ &
                        ' force constant     = ',e14.4,/ &
                        ' reference distance = ',e14.4,/ &
                        ' reference point: (',3f14.4,')',/ &
                        ' for ',i4,' out of ',i5,' total selected atoms.')
                endif
#if KEY_DOMDEC==1
                call domdec_system_changed()
#endif
             else
                if (wrnlev >= 2) call wrndie(0,'<cstran>', &
                     ' No atoms selected.')
                write(outu,'(" HMCM>>>> last restraint will be ignored")')
                nhmcm=nhmcm-1
             endif
          endif
       endif
       ! ----------------------------------------------------------------
    else if(wrd == 'HMCR' .and. qhmcm) then
       ! set up relative center of mass harmonic constraints
       if (indxa(comlyn,comlen,'CLEA') > 0) then
          nhmcmr=0
          call deallocate_hmcm
       else
          nhmcmr=nhmcmr+1
          if (nhmcmr >= maxrcm) call allocate_hmcm
          khmcmr(nhmcmr)=gtrmf(comlyn,comlen,'FORC',zero)
          hmcmr(nhmcmr)=gtrmf(comlyn,comlen,'REFR',fmark)

          if (hmcmr(nhmcmr) <= fmark) then
             if (wrnlev >= 2) call wrndie(0,'<cstran>', &
                  ' Reference distance between centers of mass must ' &
                  //'be specified.')
             nhmcmr=nhmcmr-1
          else
             ihmcm1(nhmcmr)=gtrmi(comlyn,comlen,'CM1',1)
             ihmcm2(nhmcmr)=gtrmi(comlyn,comlen,'CM2',1)
             write(outu,192) &
                  nhmcmr,ihmcm1(nhmcmr),ihmcm2(nhmcmr), &
                  khmcmr(nhmcmr),hmcmr(nhmcmr)
192          format(' Relative center of mass/geometry ', &
                  'Harmonic restraint #',i3, &
                  ' between ',i5,' and ',i5, &
                  ' Force constant = ',e14.4, &
                  ' Reference distance: ',f14.4)
          endif
       endif
#endif 
!-----------------------------------------------------------------------
       !-----------------------------------------------------------------------
       ! Taehoon Kim/Wonpil Im, KU, 2011
#if KEY_CONSHELIX==1
      ELSE IF(WRD.EQ.'HELI') THEN
         lcom=.false.
         call conshelix(lcom)
      ELSE IF(WRD.EQ.'COM ') THEN
         lcom=.true.
         call conshelix(lcom)
#endif 
!-----------------------------------------------------------------------
       !-----------------------------------------------------------------------
       ! Sean Law/Michael Feig, MSU, 2007
#if KEY_CPATH==1
    ELSE IF (WRD == 'PATH') THEN
       !     Set up PATH Constraints
       if(.not. allocated(PATHNATMS)) then
          call cpath_init
       end if
       IF (IOLEV > 0) THEN
          IF (INDXA(COMLYN, COMLEN, 'CLEA')  >  0) THEN
             if(allocated(PATHNATMS)) then
                call cpath_dealloc
             end if
             PATHN=0
             PATHNCOUNT=0
          ELSEIF ((INDXA(COMLYN,COMLEN, 'POST') > 0).AND. &
               (PATHN > 0)) THEN
             CALL PATHEPROJ (CPATHENER,X,Y,Z,DX,DY,DZ)
             IF (PRNLEV > 2)  &
                  WRITE(OUTU,*) "CSTRAN>     PATH COM projection"
             CALL PATHTEST (OUTU)
             IF (PRNLEV > 2)  &
                  WRITE(OUTU,*) "CSTRAN>     PATH Energy"
             WRITE(OUTU,*) 'CPATHENERGY = ',CPATHENER
          ELSE
             IF (PATHN >= MAXNPATHS) THEN
                IF (PRNLEV > 2) WRITE(OUTU,1787) PATHN,MAXNPATHS
1787            FORMAT("PATHN=",I3,", MAXNPATHS=",I3)
                IF(WRNLEV >= 2) CALL WRNDIE (0,'<CSTRAN>', &
                     'Too many path constraints. Increase MAXNPATHS.')
             ELSE
                PATHN=PATHN+1
                PATHNCOUNT=PATHN
                PATHNATMS(PATHN)=0
                PATHNCOM(PATHN)=0
                TESTUNIT=GTRMI(COMLYN,COMLEN,'TEST',-1)
                PROJUNIT=GTRMI(COMLYN,COMLEN,'PROJ',-1)
                IF (INDXA(COMLYN,COMLEN,'SPLI') > 0) THEN
                   CALL SELCTA(COMLYN,COMLEN,ISLCT,X,Y,Z,WT,.TRUE.)
                   DO I=1,NATOM
                      IF (ISLCT(I) > 0) THEN
                         PATHNATMS(PATHN)=PATHNATMS(PATHN)+1
                         IF (PATHNATMS(PATHN) >= MAXNATMS) THEN
                            WRITE(OUTU,3113) PATHN,PATHNATMS,MAXNATMS
3113                        FORMAT ("SPLINE=",I3,", PATHNATMS=",I3 &
                                 ,", MAXNATMS=",I3)
                            CALL WRNDIE (0,'<CSTRAN>', &
                                 'Too many atoms defining a path. Increase MAXNATMS.')
                         ENDIF
                         PATHD(1,PATHNATMS(PATHN),PATHN)=X(I)
                         PATHD(2,PATHNATMS(PATHN),PATHN)=Y(I)
                         PATHD(3,PATHNATMS(PATHN),PATHN)=Z(I)
                      ENDIF
                   ENDDO
                ELSE

                   !     IF UNIT IS SPECIFIED, READ COORDINATES FROM FILE
                   !     THE FILE FORMAT IS:
                   !     3F8.3 WITH NO LEADING SPACES AT THE START OF EACH LINE
                   !
                   !     EXAMPLE:
                   !     100.123  -23.099  16.442
                   !     99.201  -83.841  10.506
                   !
                   SPLUNIT=GTRMI(COMLYN,COMLEN,'SPLU',-1)
                   IF (SPLUNIT < 0) THEN
                      CALL WRNDIE(-3,'<CSTRAN>', &
                           'Need to specify either SPLIne or SPLUnit.')
                   ENDIF

                   IF (PRNLEV > 2) WRITE(OUTU,*)  &
                        "CSTRAN>    Reading spline points from unit ", &
                        SPLUNIT
                   PATHIOSTAT=0
                   DO WHILE(PATHIOSTAT  == 0)
                      READ(SPLUNIT,*,IOSTAT=PATHIOSTAT) PATHFILEX, &
                           PATHFILEY,PATHFILEZ
                      IF (PATHIOSTAT == 0) THEN
                         PATHNATMS(PATHN)=PATHNATMS(PATHN)+1
                         !     WRITE(*,*) PATHFILEX,PATHFILEY,PATHFILEZ
                         IF (PATHNATMS(PATHN) >= MAXNATMS) THEN
                            IF (WRNLEV >= 2) THEN
                               WRITE(OUTU,3114) PATHN,PATHNATMS, &
                                    MAXNATMS
3114                           FORMAT ("SPLINE=",I3,", PATHNATMS=",I3 &
                                    ,", MAXNATMS=",I3)
                               CALL WRNDIE (0,'<CSTRAN>', &
                                    'Too many atoms defining a path.  Increase MAXNATMS.')
                            ENDIF
                         ENDIF
                         PATHD(1,PATHNATMS(PATHN),PATHN)=PATHFILEX
                         PATHD(2,PATHNATMS(PATHN),PATHN)=PATHFILEY
                         PATHD(3,PATHNATMS(PATHN),PATHN)=PATHFILEZ
                      ENDIF
                   ENDDO
                ENDIF
                IF (PATHNATMS(PATHN) < 4.0) THEN
                   CALL WRNDIE (0,'<CSTRAN>', &
                        'A minimum of four atoms is required to define a spline.')
                ENDIF
                !     CALL SPLINTOLD ()
                CALL SPLINT (PATHN,PATHNATMS(PATHN) &
                     ,PATHA,PATHB,PATHC,PATHD,MAXNATMS,MAXNPATHS)

                !            DO I=1,PATHNATMS(PATHN)-1
                !            WRITE(*,*) I,PATHN
                !               WRITE (*,77) PATHA(1,I,PATHN),PATHB(1,I,PATHN)
                !     $              ,PATHC(1,I,PATHN),PATHD(1,I,PATHN)
                ! 77            FORMAT (4F8.3)
                !               WRITE (*,77) PATHA(2,I,PATHN),PATHB(2,I,PATHN)
                !     $              ,PATHC(2,I,PATHN),PATHD(2,I,PATHN)
                !               WRITE (*,77) PATHA(3,I,PATHN),PATHB(3,I,PATHN)
                !     $              ,PATHC(3,I,PATHN),PATHD(3,I,PATHN)
                !            ENDDO
                DO WHILE (CURRA4(COMLYN,COMLEN)  ==  'FORC')
                   PATHNCOM(PATHN)=PATHNCOM(PATHN)+1
                   !     Get Force constant from command
                   PATHK(PATHN,PATHNCOM(PATHN)) &
                        =GTRMF(COMLYN,COMLEN,'FORC',ZERO)
                   !     Get TVAL from command
                   PATHTZERO(PATHN,PATHNCOM(PATHN))= &
                        GTRMF(COMLYN,COMLEN,'TZER',ZERO)
                   !     Get COM selection from command
                   CALL SELCTA(COMLYN,COMLEN,ISLCT,X,Y,Z,WT,.TRUE.)
                   PATHNCOMATMS=0
                   PATHCOMI(PATHN,PATHNCOM(PATHN))=PATHCOMLYN+1
                   PATHCOMILEN(PATHN,PATHNCOM(PATHN))=0
                   PATHCOMX=0
                   PATHCOMY=0
                   PATHCOMZ=0
                   DO I=1,NATOM
                      IF(ISLCT(I) > 0) THEN
                         IF(PATHNCOMATMS >= MAXCOMATMS) THEN
                            IF(WRNLEV >= 2) CALL WRNDIE(0,'<CSTRAN>', &
                                 'Too many atoms defining a path. Increase MAXCOMATMS.')
                         ELSE
                            PATHNCOMATMS=PATHNCOMATMS+1
                            PATHCOMLYN=PATHCOMLYN+1
                            PATHINDX(PATHCOMLYN)=I
                         ENDIF
                      ENDIF
                      !     PATHCOMILEN should be equal to the number
                      !     of atoms selected for a given COM
                      PATHCOMILEN(PATHN,PATHNCOM(PATHN))=PATHNCOMATMS
                   ENDDO
                ENDDO
                !                  CALL PATHINIPROJ(CCOMP)
                IF (PRNLEV > 2)  &
                     WRITE(OUTU,*) "CSTRAN>     Begin TINI projection"
                CALL PATHINIPROJ(CCOMP)
                IF (PRNLEV > 2) &
                     WRITE(OUTU,*) "CSTRAN>     TINI projection" &
                     //" completed"
                CCOMP='PRYM'
                IF (PRNLEV > 2) &
                     WRITE(OUTU,*) "CSTRAN>     Begin TPRIME " &
                     //"projection"
                !     PATHTPRIME is important for keeping track of the initial projection
                !     of the COM from the MAIN set when dealing with more than one
                !     umbrella sampling window.  PATHTZERO should be close to PATHTPRIME
                !     in order to give a good approximation of the tangent at PATHTZERO.
                CALL PATHINIPROJ(CCOMP)
                !                  CALL PATHEPROJ (CPATHENER,X,Y,Z,DX,DY,DZ)
                IF (PRNLEV > 2) &
                     WRITE(OUTU,*) "CSTRAN>     TPRIME projection"  &
                     //" completed"
                !     Generate output of spline interpolation and initial projection
                !     in order to visualize the user defined path(s) and projection(s).
                IF (TESTUNIT > 0) THEN
                   IF (PRNLEV > 2) &
                        WRITE(OUTU,*) "CSTRAN>     Writing spline "  &
                        //"interpolation and initial projection" &
                        //" to unit ",TESTUNIT
                   CALL PATHTEST (TESTUNIT)
                ENDIF
             ENDIF
          ENDIF
       ENDIF
#endif 
       !-----------------------------------------------------------------------
    ELSE
       ! WRD is none of the above
       IF(WRNLEV >= 2) CALL WRNDIE(0,'<CSTRAN>', &
            ' Unknown constraint/restraint specified')
       !-----------------------------------------------------------------------
    ENDIF
    !
    call cstrn2_RETURN
    return

  contains
    subroutine cstrn2_return
#if KEY_MMFF==1 || KEY_CFF==1 /*ffield_fcm*/
      IF (FFIELD == CFF .OR. FFIELD == MMFF) THEN
         call chmdealloc('cstran.src','CSTRAN','MBND',NBOND,intg=MBND)
         call chmdealloc('cstran.src','CSTRAN','MTHETA',NTHETA,intg=MTHETA)
      else
         call chmdealloc('cstran.src','CSTRAN','MBND',1,intg=MBND)
         call chmdealloc('cstran.src','CSTRAN','MTHETA',1,intg=MTHETA)
      ENDIF
#else /* (ffield_fcm)*/
      call chmdealloc('cstran.src','CSTRAN','MBND',1,intg=MBND)
      call chmdealloc('cstran.src','CSTRAN','MTHETA',1,intg=MTHETA)
#endif /* (ffield_fcm)*/
      return
    end subroutine cstrn2_return
  END SUBROUTINE CSTRN2

  SUBROUTINE MKFRAT(IMOVE,NATOM,FREEAT,NFREAT)
    !
    !     GIVEN A MARKER ARRAY FOR ALL THE ATOMS, IMOVE, WHICH SPECIFIES THE
    !     ALLOWED MOBILITY OF AN ATOM, MKFRAT (MAKE FREE ATOMS) CONSTRUCTS
    !     AN ARRAY OF ALL ATOM NUMBERS, FREEAT, WHICH ARE FREE TO MOVE.
    !     lone pairs are not included in the FREEAT listing
    !
    !     Authors: Robert Bruccoleri
    !              Bernie Brooks
    !
  use chm_kinds
    implicit none
    INTEGER NATOM,NFREAT
    INTEGER IMOVE(*)
    INTEGER FREEAT(*)
    !
    INTEGER I
    !
    NFREAT=0
    DO I=1,NATOM
       IF(IMOVE(I) <= 0) THEN
          NFREAT=NFREAT+1
          FREEAT(NFREAT)=I
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE MKFRAT

  SUBROUTINE RDCNST(UNIT,TITLE,NTITL,NCSPHI,ICS,JCS,KCS,LCS,ICCS, &
       CCSC,CCSB,CCSD,CCSW,CCSCOS,CCSSIN)
    !
    ! Read dihedral restraints
    !
  use chm_kinds
  use dimens_fcm
  use stream
  use consta
    implicit none
    !
    INTEGER UNIT,NTITL
    character(len=*) TITLE(*)
    INTEGER NCSPHI
    INTEGER ICS(*),JCS(*),KCS(*),LCS(*),CCSD(*),ICCS(*)
    real(chm_real)  CCSC(*),CCSB(*),CCSW(*),CCSCOS(*),CCSSIN(*)
    !
    INTEGER I,ITEMP,IRNAT
    !
    IF(IOLEV > 0) THEN
       CALL RDTITL(TITLE,NTITL,UNIT,0)
       CALL WRTITL(TITLE,NTITL,OUTU,1)
       READ(UNIT,110) IRNAT,NCSPHI
110    FORMAT(3I5)
       !
       IF(IRNAT > 0) THEN
          CALL WRNDIE(-3,'<RDCNST>', &
               'Harmonic restraint I/O is no longer supported')
          RETURN
       ENDIF
       ! Dihedral constraint io changed.
       !
       call csphi_ensure_capacity(NCSPHI)
       DO I=1,NCSPHI
          ! Try to help machdep. AB.
          CCSD(I)=0
          READ(UNIT,140) ITEMP,ICS(I),JCS(I),KCS(I),LCS(I),ICCS(I), &
               CCSC(I),CCSB(I),CCSD(I),CCSW(I)
140       FORMAT(6I5,2F10.5,I5,F10.5)
          ! Update cos and sin. AB.
          IF (CCSD(I) == 0) THEN
             CCSCOS(I)=COS(CCSB(I))
             CCSSIN(I)=SIN(CCSB(I))
          ELSE
             CCSCOS(I)=COS(CCSB(I)*CCSD(I)+PI)
             CCSSIN(I)=SIN(CCSB(I)*CCSD(I)+PI)
          ENDIF
       ENDDO
    ENDIF
    !
#if KEY_PARALLEL==1
    !
    CALL PSND4(NCSPHI,1)
    CALL PSND4(ICS,NCSPHI)
    CALL PSND4(JCS,NCSPHI)
    CALL PSND4(KCS,NCSPHI)
    CALL PSND4(LCS,NCSPHI)
    CALL PSND4(ICCS,NCSPHI)
    CALL PSND8(CCSC,NCSPHI)
    CALL PSND8(CCSB,NCSPHI)
    ! AB.
    CALL PSND8(CCSCOS,NCSPHI)
    CALL PSND8(CCSSIN,NCSPHI)
    ! AB.
#endif 
    !
    RETURN
  END SUBROUTINE RDCNST

  SUBROUTINE WRCNST(UNIT,TITLE,NTITL,QCNSTR,LCIC,QQCNST, &
       NCSPHI,ICS,JCS,KCS,LCS,ICCS,CCSC,CCSB,CCSD,CCSW)
    !
    ! Write out diheral restraints
    !
    use chm_kinds
    use dimens_fcm
    use stream
    implicit none
    !
    INTEGER UNIT,NTITL
    character(len=*) TITLE(*)
    LOGICAL QCNSTR,LCIC,QQCNST
    INTEGER NCSPHI,ICS(*),JCS(*),KCS(*),LCS(*),CCSD(*),ICCS(*)
    real(chm_real)  CCSC(*),CCSB(*),CCSW(*)
    !
    INTEGER I
    !
    IF(QCNSTR) CALL WRNDIE(0,'<WRCNST>', &
         'Harmonic restrain I/O is no longer supported')
    IF (LCIC) CALL WRNDIE(0,'<WRCNST>','Cannot write IC restraints')
    IF (QQCNST) CALL WRNDIE(0,'<WRCNST>', &
         'Cannot write droplet restraint')
    !
    IF(IOLEV < 0) RETURN
    CALL WRTITL(TITLE,NTITL,UNIT,0)
    WRITE(UNIT,110) 0,NCSPHI,0
110 FORMAT(3I5)
    DO I=1,NCSPHI
       WRITE(UNIT,140) I,ICS(I),JCS(I),KCS(I),LCS(I),ICCS(I), &
            CCSC(I),CCSB(I),CCSD(I),CCSW(I)
140    FORMAT(6I5,2F10.5,I5,F10.5)
    ENDDO
    RETURN
  END SUBROUTINE WRCNST

  SUBROUTINE PRCNST(UNIT, &
       NCSPHI,ICS,JCS,KCS,LCS,CCSC,CCSB,CCSD,CCSW, &
       NATOM,QCNSTR,NUMHSETS,IHSET,TYPHSET,KCNSTR, &
       KCEXPN,LCIC,CCBIC,CCTIC,CCPIC,CCIIC, &
       QQCNST,LQMASS,KQCNST,KQEXPN)
    !
    ! Print restrain energy terms
    !
    use chm_kinds
    use consta
    use stream
    use chutil,only:atomid
    use pucker_mod,only: ncspuck,print_pucker_constraints
    implicit none
    !
    INTEGER UNIT,NCSPHI,ICS(*),JCS(*),KCS(*),LCS(*)
    real(chm_real)  CCSC(*),CCSB(*),CCSW(*)
    INTEGER CCSD(*),NATOM
    LOGICAL QCNSTR
    INTEGER NUMHSETS,IHSET(*),TYPHSET(*)
    real(chm_real)  KCNSTR(*)
    INTEGER KCEXPN(*)
    LOGICAL LCIC
    real(chm_real)  CCBIC,CCTIC,CCPIC,CCIIC
    LOGICAL QQCNST,LQMASS
    real(chm_real)  KQCNST
    INTEGER KQEXPN
    !
    real(chm_real) AMIN
    INTEGER ICSPHI
    LOGICAL FOUND,FOUNDX
    character(len=8), dimension(6) :: SID,RID,REN,AC
    character(len=8) SETTYP
    INTEGER NINSET(NUMHSETS), NPAIR
    INTEGER I,J,K,ISET,JSET
    !
    IF(OUTU /= UNIT .AND. IOLEV < 0) RETURN
    IF(OUTU == UNIT .AND. PRNLEV < 3) RETURN
    !
    FOUND=.FALSE.
    FOUNDX=.FALSE.
    !
    !-----------------------------------------------------------------------
    ! print diherdal restraints
    IF(NCSPHI > 0) THEN
       FOUNDX=.TRUE.
       WRITE(UNIT,209)
209    FORMAT(/6X,'FORCE CONSTANTS, ANGLES AND PERIODICITY', &
            ' FOR DIHEDRAL CONSTRAINTS'/)
       DO ICSPHI=1,NCSPHI
          AMIN=CCSB(ICSPHI)*180.0/PI
          ! Seemed to be a bug in the following. Corrected AB.
          I=ICS(ICSPHI)
          CALL ATOMID(I,SID(1),RID(1),REN(1),AC(1))
          I=JCS(ICSPHI)
          CALL ATOMID(I,SID(2),RID(2),REN(2),AC(2))
          I=KCS(ICSPHI)
          CALL ATOMID(I,SID(3),RID(3),REN(3),AC(3))
          I=LCS(ICSPHI)
          CALL ATOMID(I,SID(4),RID(4),REN(4),AC(4))
          WRITE(UNIT,425) ICSPHI, &
               SID(1)(1:idleng),RID(1)(1:idleng),AC(1)(1:idleng), &
               SID(2)(1:idleng),RID(2)(1:idleng),AC(2)(1:idleng), &
               SID(3)(1:idleng),RID(3)(1:idleng),AC(3)(1:idleng), &
               SID(4)(1:idleng),RID(4)(1:idleng),AC(4)(1:idleng), &
               CCSC(ICSPHI),AMIN,CCSD(ICSPHI),CCSW(ICSPHI)
425       FORMAT(I4,':',3X,3(A,' ',A,' ',A,'/ '),(A,' ',A,' ',A), &
               /,12X,'FORCE=',F9.2,',  MIN=',F9.2, &
               ',  PERIod=',I4,',  WIDTh=',F9.2)
       ENDDO
    ENDIF
    !
    !-----------------------------------------------------------------------
    ! print puckering restraints
    if(ncspuck > 0) call print_pucker_constraints(unit)
    !
    !-----------------------------------------------------------------------
    IF(QCNSTR) THEN
       WRITE(UNIT,310) NUMHSETS
310    FORMAT(' Harmonic restraints active using',I4,' atom sets.'/)
       DO ISET=1,NUMHSETS
          NINSET(ISET)=0
       ENDDO
       ! Determine which sets are active and check for errors.
       DO I=1,NATOM
          J=IHSET(I)
          IF(J > NUMHSETS .OR. J < 0) THEN
             CALL WRNDIE(-4,'<PRCNST>','Bad atom set value:coding error')
             RETURN
          ENDIF
          IF(J > 0) NINSET(J)=NINSET(J)+1
       ENDDO

       DO ISET=1,NUMHSETS
          IF(NINSET(ISET) > 0) THEN
             ! print selection from active set.
             FOUNDX=.TRUE.
             IF(TYPHSET(ISET) == 0) THEN
                SETTYP='ABSOLUTE'
             ELSE IF(TYPHSET(ISET) == 1) THEN
                SETTYP='BESTFIT'
             ELSE IF(TYPHSET(ISET) > 1) THEN
                SETTYP='RELATIVE'
             ELSE
                CALL WRNDIE(-4,'<PRCNST>','Bad set type:coding error')
                RETURN
             ENDIF
             WRITE(OUTU,320) ISET,SETTYP,KCEXPN(ISET)
320          FORMAT('     Atoms of harmonic restraint set',I4,' type: ', &
                  A,', exponent:',I3)
             !
             ! print absolute and best fit in a simple column
             IF(TYPHSET(ISET) <= 1) THEN
                DO I=1,NATOM
                   IF(IHSET(I) == ISET) THEN
                      IF(KCNSTR(I) /= 0.0) THEN
                         CALL ATOMID(I,SID(1),RID(1),REN(1),AC(1))
                         WRITE(OUTU,345) SID(1)(1:idleng),RID(1)(1:idleng), &
                              REN(1)(1:idleng),AC(1)(1:idleng),KCNSTR(I)
345                      FORMAT(10X,4(1X,A),F12.5)
                      ENDIF
                   ENDIF
                ENDDO
                !
                ! print relative as matched pairs.
             ELSE
                JSET=TYPHSET(ISET)
                IF(TYPHSET(JSET) /= ISET) CALL WRNDIE(-4,'<PRCNST>', &
                     'Error in relative set types - coding error')
                IF(ISET < JSET) THEN  ! only process once (ISET<JSET)
                   NPAIR=NINSET(ISET)
                   IF(NPAIR /= NINSET(JSET)) CALL WRNDIE(-2,'<PRCNST>', &
                        'Number of selected atoms in two sets do not match')
                   IF(NPAIR > NINSET(JSET)) NPAIR=NINSET(JSET)
                   !
                   I=1
                   J=1
                   DO K=1,NPAIR
                      DO WHILE(IHSET(I) /= ISET)
                         I=I+1
                      ENDDO
                      DO WHILE(IHSET(J) /= JSET)
                         J=J+1
                      ENDDO
                      IF(KCNSTR(I) /= 0.0) THEN
                         CALL ATOMID(I,SID(1),RID(1),REN(1),AC(1))
                         CALL ATOMID(J,SID(2),RID(2),REN(2),AC(2)) ! brb fix 24-Jul-2003
                         WRITE(OUTU,355) SID(1)(1:idleng),RID(1)(1:idleng), &
                              REN(1)(1:idleng),AC(1)(1:idleng), &
                              SID(2)(1:idleng),RID(2)(1:idleng), &
                              REN(2)(1:idleng),AC(2)(1:idleng),KCNSTR(I)
355                      FORMAT(10X,4(1X,A),' with ',4(1X,A),F12.5)
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF
          ENDIF
       ENDDO
    ENDIF
    !-----------------------------------------------------------------------
    ! Print internal coordinate restraint parameters
    IF(LCIC) THEN
       FOUNDX=.TRUE.
       WRITE(UNIT,132) CCBIC,CCTIC,CCPIC,CCIIC
132    FORMAT(' Internal Coordinate restraints active.', &
            ' FORCE CONSTANTS:'/,'   BONDS =',F14.6, &
            '  ANGLES =',F14.6,'  DIHEDRALS =',F14.6, &
            '  IMPROPERS =',F14.6)
    ENDIF
    !
    !-----------------------------------------------------------------------
    ! Print quartic droplet restraint parameters
    IF(QQCNST) THEN
       FOUNDX=.TRUE.
       IF(LQMASS) THEN
          WRITE(UNIT,142) KQCNST,KQEXPN
142       FORMAT(' Droplet restraint active with mass weighting.'/, &
               '   FORCE CONSTANT =',E14.4,'  EXPONENT =',I4)
       ELSE
          WRITE(UNIT,143) KQCNST,KQEXPN
143       FORMAT(' Droplet restraint active. No mass weighting.'/, &
               '   FORCE CONSTANT =',E14.4,'  EXPONENT =',I4)
       ENDIF
    ENDIF
    !
    !-----------------------------------------------------------------------
    ! Print shake constraints
    CALL PRNSHK(FOUND)
    !
    !-----------------------------------------------------------------------
    IF (.NOT.(FOUND.OR.FOUNDX)) WRITE(UNIT,240)
240 FORMAT(' There are no constraints/restraints.')
    IF(.NOT.FOUND) RETURN
    !-----------------------------------------------------------------------
    IF(QCNSTR) THEN
       WRITE(UNIT,250)
250    FORMAT(' The harmonic restraint flag is turned ON.')
    ELSE
       WRITE(UNIT,260)
260    FORMAT(' The harmonic restraint flag is turned OFF.')
    ENDIF
    !-----------------------------------------------------------------------
    RETURN
  END SUBROUTINE PRCNST

  subroutine csphi_ensure_capacity(wantcap)
    use cnst_fcm
    use memory
    integer, intent(in) :: wantcap
    integer, parameter :: mincap = 20
    real, parameter :: pad = 1.25
    integer :: newcap

    if (allocated(ICCS)) then
       if (wantcap > size(ICCS)) then
          newcap = ceiling(pad * wantcap)
          call chmrealloc('cstran.src','csphi_ensure_capacity','ICCS',newcap,intg=ICCS)
          call chmrealloc('cstran.src','csphi_ensure_capacity','ICS',newcap,intg=ICS)
          call chmrealloc('cstran.src','csphi_ensure_capacity','JCS',newcap,intg=JCS)
          call chmrealloc('cstran.src','csphi_ensure_capacity','KCS',newcap,intg=KCS)
          call chmrealloc('cstran.src','csphi_ensure_capacity','LCS',newcap,intg=LCS)
          call chmrealloc('cstran.src','csphi_ensure_capacity','CCSD',newcap,intg=CCSD)
          call chmrealloc('cstran.src','csphi_ensure_capacity','CCSC',newcap,crl=CCSC)
          call chmrealloc('cstran.src','csphi_ensure_capacity','CCSB',newcap,crl=CCSB)
          call chmrealloc('cstran.src','csphi_ensure_capacity','CCSW',newcap,crl=CCSW)
          call chmrealloc('cstran.src','csphi_ensure_capacity','CCSCOS',newcap,crl=CCSCOS)
          call chmrealloc('cstran.src','csphi_ensure_capacity','CCSSIN',newcap,crl=CCSSIN)
       endif
    else
       newcap = max(mincap, ceiling(pad * wantcap))
       call chmalloc('cstran.src','csphi_ensure_capacity','ICCS',newcap,intg=ICCS)
       call chmalloc('cstran.src','csphi_ensure_capacity','ICS',newcap,intg=ICS)
       call chmalloc('cstran.src','csphi_ensure_capacity','JCS',newcap,intg=JCS)
       call chmalloc('cstran.src','csphi_ensure_capacity','KCS',newcap,intg=KCS)
       call chmalloc('cstran.src','csphi_ensure_capacity','LCS',newcap,intg=LCS)
       call chmalloc('cstran.src','csphi_ensure_capacity','CCSD',newcap,intg=CCSD)
       call chmalloc('cstran.src','csphi_ensure_capacity','CCSC',newcap,crl=CCSC)
       call chmalloc('cstran.src','csphi_ensure_capacity','CCSB',newcap,crl=CCSB)
       call chmalloc('cstran.src','csphi_ensure_capacity','CCSW',newcap,crl=CCSW)
       call chmalloc('cstran.src','csphi_ensure_capacity','CCSCOS',newcap,crl=CCSCOS)
       call chmalloc('cstran.src','csphi_ensure_capacity','CCSSIN',newcap,crl=CCSSIN)
    endif
  end subroutine csphi_ensure_capacity

  subroutine csphi_clear()
    use cnst_fcm
    use memory
    integer :: oldcap
    if (allocated(ICCS)) then
       oldcap = size(ICCS)
       call chmdealloc('cstran.src','csphi_clear','ICCS',oldcap,intg=ICCS)
       call chmdealloc('cstran.src','csphi_clear','ICS',oldcap,intg=ICS)
       call chmdealloc('cstran.src','csphi_clear','JCS',oldcap,intg=JCS)
       call chmdealloc('cstran.src','csphi_clear','KCS',oldcap,intg=KCS)
       call chmdealloc('cstran.src','csphi_clear','LCS',oldcap,intg=LCS)
       call chmdealloc('cstran.src','csphi_clear','CCSD',oldcap,intg=CCSD)
       call chmdealloc('cstran.src','csphi_clear','CCSC',oldcap,crl=CCSC)
       call chmdealloc('cstran.src','csphi_clear','CCSB',oldcap,crl=CCSB)
       call chmdealloc('cstran.src','csphi_clear','CCSW',oldcap,crl=CCSW)
       call chmdealloc('cstran.src','csphi_clear','CCSCOS',oldcap,crl=CCSCOS)
       call chmdealloc('cstran.src','csphi_clear','CCSSIN',oldcap,crl=CCSSIN)
    endif
  end subroutine csphi_clear

  subroutine hset_ensure_capacity(wantcap)
    use cnst_fcm
    use memory
    integer, intent(in) :: wantcap
    integer, parameter :: mincap = 20
    real, parameter :: pad = 1.25
    integer :: newcap

    if (allocated(KCEXPN)) then
       if (wantcap > size(KCEXPN)) then
          newcap = ceiling(pad * wantcap)
          call chmrealloc('cstran.src','hset_ensure_capacity','KCEXPN',newcap,intg=KCEXPN)
          call chmrealloc('cstran.src','hset_ensure_capacity','XHSCALE',newcap,crl=XHSCALE)
          call chmrealloc('cstran.src','hset_ensure_capacity','YHSCALE',newcap,crl=YHSCALE)
          call chmrealloc('cstran.src','hset_ensure_capacity','ZHSCALE',newcap,crl=ZHSCALE)
          call chmrealloc('cstran.src','hset_ensure_capacity','TYPHSET',newcap,intg=TYPHSET)
          call chmrealloc('cstran.src','hset_ensure_capacity','QHNORT',newcap,log=QHNORT)
          call chmrealloc('cstran.src','hset_ensure_capacity','QHNOTR',newcap,log=QHNOTR)
       endif
    else
       newcap = max(mincap, ceiling(pad * wantcap))
       call chmalloc('cstran.src','hset_ensure_capacity','KCEXPN',newcap,intg=KCEXPN)
       call chmalloc('cstran.src','hset_ensure_capacity','XHSCALE',newcap,crl=XHSCALE)
       call chmalloc('cstran.src','hset_ensure_capacity','YHSCALE',newcap,crl=YHSCALE)
       call chmalloc('cstran.src','hset_ensure_capacity','ZHSCALE',newcap,crl=ZHSCALE)
       call chmalloc('cstran.src','hset_ensure_capacity','TYPHSET',newcap,intg=TYPHSET)
       call chmalloc('cstran.src','hset_ensure_capacity','QHNORT',newcap,log=QHNORT)
       call chmalloc('cstran.src','hset_ensure_capacity','QHNOTR',newcap,log=QHNOTR)
    endif
  end subroutine hset_ensure_capacity

  subroutine hset_clear()
    use cnst_fcm
    use memory
#if KEY_OPENMM==1
    use omm_ctrl, only: omm_system_changed
#endif
    integer :: oldcap
    if (allocated(KCEXPN)) then
       oldcap = size(KCEXPN)
       call chmdealloc('cstran.src','hset_clear','KCEXPN',oldcap,intg=KCEXPN)
       call chmdealloc('cstran.src','hset_clear','XHSCALE',oldcap,crl=XHSCALE)
       call chmdealloc('cstran.src','hset_clear','YHSCALE',oldcap,crl=YHSCALE)
       call chmdealloc('cstran.src','hset_clear','ZHSCALE',oldcap,crl=ZHSCALE)
       call chmdealloc('cstran.src','hset_clear','TYPHSET',oldcap,intg=TYPHSET)
       call chmdealloc('cstran.src','hset_clear','QHNORT',oldcap,log=QHNORT)
       call chmdealloc('cstran.src','hset_clear','QHNOTR',oldcap,log=QHNOTR)
#if KEY_OPENMM==1
       call omm_system_changed()  
#endif
    endif
  end subroutine hset_clear

end module cstran_mod

