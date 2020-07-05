!-----------------------------------------------------------------------
!     Generate replicas: one per process or group of processes,
!     and maybe many replicas per process.
!     By Paul Maragakis and Milan Hodoscek, 2005
!     Support for parallel/parallel and I/O setup: Milan Hodoscek, May 2007
!
!     Milan Hodoscek, April 2009:
!      - replica exchange printout and fixes
!      - SGLD replica exchange
!      - support for equal temperatures of selected/all replicas
!      - initial implementation of TIGER method.
!      TIGER is based on the code from Satoru Itoh
!
MODULE REPDSTRMOD
  use chm_kinds
  use parallel, only: maxnode
#if KEY_REPDSTR==1
  use repdstr, only: maxrepdstr
#endif
#if KEY_OPENMM==1
  use omm_glblopts, only : qtor_repex, torsion_lambda
#endif

  use param_store, only: set_param

  implicit none

  integer,save                                         :: reptag

  !
  ! NB curtemps maps indices -> temperatures,
  real(chm_real),save                                  :: tempcurrent               ! for fast repexch
  real(chm_real),allocatable,dimension(:)              :: curtemps                  ! also for fast repexch
  logical                                              :: ffrun                     ! fast repexch logical for first run
  integer,save                                         :: loglevel                  ! fastrepexch how many exchanges be written
  
  real(chm_real),allocatable,dimension(:),save,public  :: sgtemprx,sgftrx
  real(chm_real),allocatable,dimension(:),save         :: rhener, rlener            ! energies of what is in the reservoir
  real(chm_real4),allocatable,dimension(:),save        :: rescrdx,rescrdy,rescrdz   ! coordinates from reservoir
                                                                                    ! traj stuff in single prec

  logical,save                                         :: qecor

  ! Track and export exchange probabilities
  real(chm_real):: EXRATUP,EXRATUP2,EXRATDN,EXRATDN2,RERAT,RERAT2,REPROB,REPROB2
  integer,save  :: noppup,noppdn,nsucup,nsucdn,noppup2,noppdn2,nsucup2,nsucdn2 

  !!!!!!!!!!!!!!!!!!!!!!
  !! Begin 2D-REX     !!
  !!!!!!!!!!!!!!!!!!!!!!
  logical,save                                         :: q2drex,q2ditemp
  integer,save                                         :: d1freq,d2freq,ndim1,ndim2
  integer,save                                         :: myrepd1,myrepd2,irexd1,irexd2
  integer,save                                         :: nbrup1,nbrdn1,nbrup2,nbrdn2,cnbrup1,cnbrdn1,cnbrup2,cnbrdn2
  logical,save,dimension(2)                            :: q2dtemp,q2dham,q2dph
  !! End 2D-REX

  !
  integer :: repdsynccount = 0
  real(chm_real),save :: rhtemp, rltemp ! temperatures of the reservoirs
  logical,save,public :: qsump,qrxsgld,qrxtham,qreservoir,qreshigh,qreslow,qfastrepdstr
#if KEY_PHMD==1
  logical,save,public :: qphrx 
#endif
#if KEY_BLOCK==1
  logical,save,public :: QMSPHRX    ! GG MSLD-compatibility
#endif

  ! discrete-state constant pH related variables
  logical,save,public :: qphrex
  real(chm_real),save,public  :: phrx(maxnode)
  real(chm_real)              :: rhph,rlph      ! pH values of the reservoir
  integer,allocatable,target,dimension(:)   :: nprothigh,nprotlow ! store number of protonated groups in each reservoir entry
  integer,allocatable,target,dimension(:,:) :: reshtstate,resltstate
  real(chm_real4),allocatable,dimension(:)  :: rescg

  logical,save,public :: qresboltz,qresnobo ! what exchange criterion to use for reservoir
  integer,save,public :: irex = 0     ! number of exchanges
  integer,save,public :: isuc = 0     ! number of successful exchanges
  integer,save,public :: iresexch = 0 ! number of reservoir exchanges
  integer,save,public :: rhunit,rlunit ! units for high and low resorvoirs
  integer,save,public :: highressz,lowressz,maxressz ! current and maximum size(s) of the 
                                                     ! reservoir(s) in adaptive calculations
  integer,save,public :: repdid     ! replica id at current replica
  integer,save,public :: nrepeat    ! number of times to repeat exchange procedure
  !
  !     TIGER stuff:
  !     (Temperature Intervals with Global Energy Reassignment)
  ! 
  !                     __,,,,_
  !          _ __..-;''`--/'/ /.',-`-.
  !      (`/' ` |  \ \ \\ / / / / .-'/`,_
  !     /'`\ \   |  \ | \| // // / -.,/_,'-,
  !    /<7' ;  \ \  | ; ||/ /| | \/    |`-/,/-.,_,/')
  !   /  _.-, `,-\,__|  _-| / \ \/|_/  |    '-/.;.\'
  !   `-`  f/ ;      / __/ \__ `/ |__/ |
  !        `-'      |  -| =|\_  \  |-' |
  !              __/   /_..-' `  ),'  //
  !             ((__.-'((___..-'' \__.'
  !
  !
  !     qrxtiger  - logical: are we using TIGER
  !     qpxtiger  - logical: do we need more dynamics before exchange
  !     qrxtmin   - logical: do we need to perform minimization at this step
  !     tigergr   - real: gradient tolerance for the minimization step
  !     tigerit   - integer: how many mini&equil cycles are needed
  !     tigeriti  - integer: how many cycles already performed
  !     tigernm   - integer: number of minimization steps
  !     tigerneq  - integer: number of equilibration steps
  !
  logical,save,public :: qrxtiger,qpxtiger=.false.,qrxtmin=.false.
  real(chm_real),save,public :: TIGERGR
  integer,save,public :: TIGERIT,TIGERNM,TIGERNEQ,TIGERITI=0
  !
CONTAINS
  !
#if KEY_REPDSTR==1 /* repdstr_main */
  SUBROUTINE REPDSTRMAIN
    !-----------------------------------------------------------------------
    !
#if KEY_PARALLEL==1 /* pll */
    !-----------------------------------------------------------------------
  use number
  use dimens_fcm
  use comand
  use cstuff, only: getpid
  use psf
  use coord
  use parallel
  use repdstr
  use stream
  use string
  use memory
  use deriv
  use block_ltm
  use lambdam    !GG MSLD-compatibility
  use consta     !GG MSLD-compatibility

  implicit none

    logical :: want_openmm
    LOGICAL QFIN,QFINSYNC,QFINONE,QSYNC
    !
    real(chm_real) STEMP,DTEMP,MTEMP,SGTEMP,DSGTEMP,MSGTEMP,SGFT,DSGFT, ttemp
    !
    character(len=100) FNAME
    INTEGER J !GG New MSLD pH-REX commands added
    INTEGER I,IOERR,MYPID,ENEUN,PROBUNIT,NAT3,RDIM
    !
    QFIN=(INDXA(COMLYN,COMLEN,'RESE').GT.0)
    !
    IF(QFIN)THEN
       !        Must be called before qrepdstr=.false.
       !        Restore the I/O; can be called individually! No global comm.!
       CALL DREPRESIO(IOLEV,PRNLEV,WRNLEV)
       QREPDSTR = .FALSE.
       QRDQTT = .FALSE.
       QWRQTT = .FALSE.
       !!         IF(MYNODG.EQ.0)QRDQTT=.TRUE.
       CALL PSETGLOB
       !
       ! is iparpt OK?
       ! is ippmap OK?
       ! is INODE() OK?
       !
       QFINSYNC=(INDXA(COMLYN,COMLEN,'SYNC').GT.0)
       IF(QFINSYNC)CALL PSYNC()
       QFINONE=(INDXA(COMLYN,COMLEN,'PONE').GT.0)
       IF(QFINONE)THEN
          !           IF(MYNODG.EQ.0)THEN
          MYNOD=0
          NUMNOD=1
          CALL CUBE(mynod,numnod,ippmap)
          !           ENDIF
       ENDIF
       !!         write(outu,'(A,L)')'REPDSTR>after sync:qfinsync=',qfinsync
       if(qrxsgld)deallocate(sgtemprx)
       if(qrxsgld)deallocate(sgftrx)

       IF(QFASTREPDSTR) &
          call chmdealloc('repdstr.src','REPDSTRMAIN','curtemps',NREPDSTR,crl=curtemps)
       IF(QRESERVOIR) THEN
          call chmdealloc('repdstr.src','REPDSTRMAIN','rescrdx',NATOM,cr4=rescrdx)
          call chmdealloc('repdstr.src','REPDSTRMAIN','rescrdy',NATOM,cr4=rescrdy)
          call chmdealloc('repdstr.src','REPDSTRMAIN','rescrdz',NATOM,cr4=rescrdz)
          call chmdealloc('repdstr.src','REPDSTRMAIN','rescg',NATOM,cr4=rescg)
       ENDIF
       RETURN
    ENDIF

    !
    !
    !     This is global operation so every input script
    !     should have MATCHING repd sync!!!
    !
    QSYNC=(INDXA(COMLYN,COMLEN,'SYNC').GT.0)

    IF(QSYNC)THEN
       repdsynccount=repdsynccount+1
       if(prnlev.gt.2) then
          write(outu,'(a,i5)')'REPDSTR>sync count=',repdsynccount
       endif
       CALL PSETGLOB()
       CALL PSYNC()
       CALL PSETLOC()
       RETURN
    ENDIF

    !   sometimes we want this in our scripts :-)
    if (indxa(comlyn,comlen,'IORES').gt.0) then
       qrepioset = .false.
       return
    endif

    if (indxa(comlyn,comlen,'IOSET').gt.0) then
       qrepioset = .true.
       return
    endif
    !
    NREPDSTR = GTRMI(COMLYN,COMLEN,'NREP',1)
    !
    NATREPCMD = GTRMI(COMLYN,COMLEN,'NATR',-1)
    !
    QREPDSTR = .TRUE.
    !
    QFASTREPDSTR = (INDXA(COMLYN,COMLEN,'FAST').GT.0)
    IF(QFASTREPDSTR) THEN
       CALL CHMALLOC('repdstr.src','REPDSTRMAIN','CURTEMPS',NREPDSTR,crl=curtemps)
       FFRUN=.TRUE.

       LOGLEVEL=GTRMI(COMLYN,COMLEN,'LOGL',1)
    ENDIF

    qrepioset = .true.
    !
    !     This is the default situation. This changes when stream is called!
    QRDQTT = .FALSE.
    !     This is the default situation. This changes when outu is called!
    QWRQTT = .FALSE.

    NREPEAT=GTRMI(COMLYN,COMLEN,'REPE',1)
    IF(PRNLEV.GT.2) WRITE(OUTU,'(A,I3)') 'REPDSTR> EXCHANGE REPEAT = ', NREPEAT


    Q2DREX = (INDXA(COMLYN,COMLEN,'TWOD').GT.0)
    IF(Q2DREX) THEN
       CALL SETUP_2D(COMLYN,COMLEN)
       RETURN
    ENDIF

    ! Tim Miller: June, 2011
    ! Made default feature: October, 2014
    ! Decide if we need to do discrete state PH based replica exchange
    QPHREX=(INDXA(COMLYN,COMLEN,'PHRE').GT.0)

    IF(QPHREX) THEN
       !    repd nrep <int> phrex freq <int> phval <real> phval <real> ...

       IF(QFASTREPDSTR) &
          CALL WRNDIE(-4,'<REPDSTR>','FAST PH REX IS NOT SUPPORTED')

       REPSEED=123+IREPDSTR
       IUNREX=GTRMI(COMLYN,COMLEN,'UNIT',OUTU)
       IREXFQ=GTRMI(COMLYN,COMLEN,'FREQ',1)
      
       DO I=1,NREPDSTR
          PHRX(I)=GTRMF(COMLYN,COMLEN,'PHVA',-ONE)
          !! Fix for Ana -- allow negative pH values
          !!IF(PHRX(I).LT.0) &
          !!   CALL WRNDIE(-5,'<REPDSTR>','NEGATIVE PH VALUES ARE NOT ALLOWED!')
       ENDDO
    ENDIF

    QREXCHG = (INDXA(COMLYN,COMLEN,'EXCH').GT.0)
    QREXCHGL = (INDXA(COMLYN,COMLEN,'EXLM').GT.0)

    ! Tim Miller: June, 2011
    ! This code takes care of reservoirs for Asim Okur's reservoir relica
    ! exchange method
    QRESERVOIR=(INDXA(COMLYN,COMLEN,'RSRV').GT.0)
    QRESHIGH=(INDXA(COMLYN,COMLEN,'RESH').GT.0)
    QRESLOW=(INDXA(COMLYN,COMLEN,'RESL').GT.0)
    IF(QRESERVOIR) THEN

       IF(QFASTREPDSTR) &
          CALL WRNDIE(-4,'<REPDSTR>','FAST RESERVOIR REX IS NOT SUPPORTED')

       QRESBOLTZ=(INDXA(COMLYN,COMLEN,'BOLT').GT.0)
       QRESNOBO=(INDXA(COMLYN,COMLEN,'NOBO').GT.0)
       IF(.not. (QRESHIGH .or. QRESLOW)) &
          CALL WRNDIE(-3,'<REPDSTRMAIN>', 'RESERVOIR NEEDS RESHIGH OR RESLOW')

       IF(QREXCHGL) THEN
          ! Only allow "Boltzmann" for Hamiltonian REX
          QRESNOBO=.FALSE.
          QRESBOLTZ=.TRUE.
          IF(PRNLEV.GE.3) WRITE(OUTU,'(A)') &
             'REPDSTRMAIN> HAMILTONIAN REX IN USE. BOLTZMANN EXCHANGE CRITERION ACTIVATED.'
       ELSE
          ! check to make sure there aren't multiple specs
          IF(QRESBOLTZ .and. QRESNOBO) THEN
             CALL WRNDIE(-4,'<REPDSTRMAIN>', 'CONFLICTING EXCHANGE CRITERIA SET FOR RESERVOIR REX')
          ENDIF
          IF(.not. ( QRESBOLTZ .or. QRESNOBO )) THEN
             CALL WRNDIE(0,'<REPDSTRMAIN>', &
                         'NO EXCHANGE CRITERIA SET: RESERVOIR REX WILL USE BOLTZMANN CRITERION')
             QRESBOLTZ=.true.
          ENDIF
       ENDIF
       IF(QRESBOLTZ.OR.QRESNOBO) THEN
          QECOR=(INDXA(COMLYN,COMLEN,'ECOR').GT.0)
       ELSE
          QECOR=.FALSE.
       ENDIF

       call chmalloc('repdstr.src','REPDSTRMAIN','rescrdx',NATOM,cr4=rescrdx)
       call chmalloc('repdstr.src','REPDSTRMAIN','rescrdy',NATOM,cr4=rescrdy)
       call chmalloc('repdstr.src','REPDSTRMAIN','rescrdz',NATOM,cr4=rescrdz)
       call chmalloc('repdstr.src','REPDSTRMAIN','rescg',NATOM,cr4=rescg)

       ! Get info on the trajectory files containing the reservoirs      
       IF(QRESHIGH) THEN
          HIGHRESSZ=GTRMI(COMLYN,COMLEN,'RHSZ',-1)
          RHUNIT=GTRMI(COMLYN,COMLEN,'RHUN',-1)
          IF(RHUNIT.LT.0) &
             CALL WRNDIE(-3,'<REPDSTRMAIN>','BAD UNIT FOR TOP RESERVOIR')

          IF(HIGHRESSZ.LT.0) CALL WRNDIE(-3,'<REPDSTRMAIN>','MUST SPECIFY A VALID HIGH RESERVOIR SIZE')
          

          IF(QREXCHGL) THEN
             RHTEMP=GTRMF(COMLYN,COMLEN,'RHTE',-1.0)
             IF(RHTEMP.LE.0) &
                CALL WRNDIE(-3,'<REPDSTRMAIN>','BAD BOLTZMANN TEMP FOR HIGH RESERVOIR')
             ENEUN=GTRMI(COMLYN,COMLEN,'FHEN',-1)
             IF(ENEUN.LT.1) &
                CALL WRNDIE(-3,'<REPDSTRMAIN>','RESERVOIR H-REX CANNOT PRECALC ENERGIES.')
             CALL GET_ENE_VAL(ENEUN,QECOR,RLTEMP,.TRUE.)
          ELSE
             RHTEMP=GTRMF(COMLYN,COMLEN,'RHTE',-1.0)
             IF(QRESBOLTZ) THEN
                IF(QPHREX) THEN
                   ! handle boltzmann pH rex, hoorah
                   RHPH=GTRMF(COMLYN,COMLEN,'RHPH',-1.0)
                ELSE
                   ! Get temp of high reservoir
                   IF(RHTEMP.LE.0) &
                      CALL WRNDIE(-3,'<REPDSTRMAIN>','BAD BOLTZMANN TEMP FOR HIGH RESERVOIR')
                ENDIF 
             ENDIF

             IF(QRESBOLTZ.OR.QRESNOBO) THEN
                ENEUN=GTRMI(COMLYN,COMLEN,'FHEN',-1)
                IF(ENEUN.GT.0) THEN
                   CALL GET_ENE_VAL(ENEUN,QECOR,RHTEMP,.TRUE.)
                ELSE
                   CALL PRECALCENE(.TRUE.,X,Y,Z,QECOR,RHTEMP)
                ENDIF
             ENDIF
          ENDIF
       ENDIF
       IF(QRESLOW) THEN
          LOWRESSZ=GTRMI(COMLYN,COMLEN,'RLSZ',-1)
          RLUNIT=GTRMI(COMLYN,COMLEN,'RLUN',-1)          
          IF(RLUNIT.LT.0) &
             CALL WRNDIE(-3,'<REPDSTRMAIN>','BAD UNIT FOR BOTTOM RESERVOIR')

          IF(LOWRESSZ.LT.0) CALL WRNDIE(-3,'<REPDSTRMAIN>','MUST SPECIFY A VALID LOW RESERVOIR SIZE')


          IF(QREXCHGL) THEN
             RLTEMP=GTRMF(COMLYN,COMLEN,'RLTE',-1.0)
             IF(RLTEMP.LE.0) &
                CALL WRNDIE(-3,'<REPDSTRMAIN>','BAD BOLTZMANN TEMP FOR LOW RESERVOIR')
             ENEUN=GTRMI(COMLYN,COMLEN,'FLEN',-1)
             IF(ENEUN.LT.1) &
                CALL WRNDIE(-3,'<REPDSTRMAIN>','RESERVOIR H-REX CANNOT PRECALC ENERGIES.')
             CALL GET_ENE_VAL(ENEUN,QECOR,RLTEMP,.FALSE.)
          ELSE
             RLTEMP=GTRMF(COMLYN,COMLEN,'RLTE',-1.0)
             IF(QRESBOLTZ) THEN
                ! Get temp of high reservoir
                IF(QPHREX) THEN
                   ! handle boltzmann pH rex, hoorah
                   RLPH=GTRMF(COMLYN,COMLEN,'RLPH',-1.0)
                ELSE
                   IF(RLTEMP.LE.0) &
                      CALL WRNDIE(-3,'<REPDSTRMAIN>','BAD BOLTZMANN TEMP FOR LOW RESERVOIR')
                ENDIF
             ENDIF


             IF(QPHREX) THEN
                IF(QFASTREPDSTR) &
                   CALL WRNDIE(-4,'<REPDSTR>','FAST PH REX IS NOT SUPPORTED')
             ELSE
                IF(QRESBOLTZ.OR.QRESNOBO) THEN
                   ENEUN=GTRMI(COMLYN,COMLEN,'FLEN',-1)
                   IF(ENEUN.GT.0) THEN
                      CALL GET_ENE_VAL(ENEUN,QECOR,RLTEMP,.FALSE.)
                   ELSE
                      CALL PRECALCENE(.FALSE.,X,Y,Z,QECOR,RLTEMP)
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
       ENDIF
    ELSE
       QRESBOLTZ=.false.
       QRESNOBO=.false.
    ENDIF

    !OLD:C
    !OLD:C     Redistribute the parallel setup
    !OLD:C     Current limitation (number of replica)/(number of processes) = 1
    !OLD:C     this way we know which replica we are dealing with, without changing
    !OLD:C     charmm PSF, COOR etc structures
    !OLD:C     If you change this rule do it also elsewhere like for example here:
    !OLD:C        - VOPEN() in machio.src
    !
    !     This is initialization. Use PSETLOC(), PSETGLOB() elsewhere
    !     This needs to be generalized!!! Also first two lines are now in paral1.src!!!
    !     No need to be here... IREPDSTR must be calculated from them
    !
    !      MYNODG=MYNOD
    !      NUMNODG=NUMNOD
    NUMNOD=NUMNODG/NREPDSTR
    IF(PRNLEV.GT.2)WRITE(OUTU,'(A,2I5)') &
         ' REPD> Number of processes and NREP=',NUMNODG,NREPDSTR
    IF(NUMNOD*NREPDSTR.NE.NUMNODG)THEN
       CALL WRNDIE(-5,'<REPDSTR>', &
            'Wrong combination of NREP and number of processes.')
    ENDIF
    IREPDSTR=MYNODG/NUMNOD
    REPDID=IREPDSTR
    !
    !      write(70+mynodg,'(a,4i7)')'REPDSTR-0>me,meg,np,npg=',
    !     $     mynod,mynodg,numnod,numnodg
    CALL PSETLOC
    !
    !      write(70+mynodg,'(a,4i7)')'REPDSTR-1>me,meg,np,npg=',
    !     $     mynod,mynodg,numnod,numnodg
    !      write(70+mynodg,'(a,4i7)')'REPDSTR-2>irepdstr,iolev,prnlev=',
    !     $     irepdstr,iolev,prnlev
    !
    !     Maybe here we could reassign mynode,numnode, too!?
    CALL set_param('MYREP',IREPDSTR)
    CALL set_param('NREP',NREPDSTR)
    !
    !C      write(*,*)'me,temp(1-5)=',mynodg,(temprx(i),i=1,5)
    !
    !C      write(*,*)'REPDSTR-1>me,meg,np,npg=',mynod,mynodg,numnod,numnodg
    !
    !     This is old, now we hope we deal with everything!!!
    !     At this point we could redefine IOLEV but there are case that it
    !     wouldn't work, so we deal with them individualy and currently it
    !     works for:
    !     1. OPEN      (machio.src: maybe we want also: VINQRE, VCLOSE)
    !     2. READ COOR (coorio.src: binary??)
    !     3. WRITE COOR
    !     4. TRAJIO
    !
    !     define the I/O flags:
    CALL DREPSETIO(IOLEV,PRNLEV,WRNLEV)
    !
    !     Replica exchange data: get the temperatures...
    QEXPT = (INDXA(COMLYN,COMLEN,'EXPT').GT.0)
    QEX2D = (INDXA(COMLYN,COMLEN,'EX2D').GT.0)
    QEXBK = (INDXA(COMLYN,COMLEN,'EXBK').GT.0)

    IF(QREXCHG) THEN
       !
       !     repd nrep <int> exch freq <int> temp <real> temp <real> ...
       !
       !     could be also:
       !     repd nrep <int> exch stemp <real> dtemp <real>
       !     where stemp is starting temperature and dtemp is the
       !     interval between the temperatures.
       !
       !         REPSEED=123+MYNODG
       REPSEED=123+IREPDSTR
       QSUMP=(INDXA(COMLYN,COMLEN,'SUMP').GT.0)
       QRXSGLD=(INDXA(COMLYN,COMLEN,'SGLD').GT.0)
       QRXTIGER=(INDXA(COMLYN,COMLEN,'TIGE').GT.0)
       QRXTHAM=(INDXA(COMLYN,COMLEN,'THAM').GT.0)
#if KEY_PHMD==1
       QPHRX=(INDXA(COMLYN,COMLEN,'PHMD').GT.0) !JAW. Flag for using PHMD and PH exchange
#endif
       IUNREX=GTRMI(COMLYN,COMLEN,'UNIT',OUTU)
       IREXFQ=GTRMI(COMLYN,COMLEN,'FREQ',1)
#if KEY_OPENMM==1
       want_openmm = (indxa(comlyn,comlen,'OMM')>0)
       qtor_repex = (indxa(comlyn,comlen,'TORS')>0)
       TTEMP=GTRMF(COMLYN,COMLEN,'TTEM',-ONE)
#endif
       STEMP=GTRMF(COMLYN,COMLEN,'STEM',-ONE)
       IF(QRXTHAM.AND.QRXSGLD) &
          CALL WRNDIE(-4,'<REPDSTR>', &
          'TEMP-HAMILTONIAN REPLICA EXCHANGE IS INCOMPATIBLE WITH SGLD REM')
       IF(STEMP>ZERO) THEN
          MTEMP=GTRMF(COMLYN,COMLEN,'MTEM',-ONE)
          IF(MTEMP>ZERO)THEN
             DTEMP=EXP(DLOG(MTEMP/STEMP)/(NREPDSTR-1))
             DO I=1,NREPDSTR
                TEMPRX(I)=STEMP*DTEMP**(I-1)
                IF(QFASTREPDSTR) CURTEMPS(I)=TEMPRX(I)
          ENDDO
          ELSE
             DTEMP=GTRMF(COMLYN,COMLEN,'DTEM',-ONE)
             IF(DTEMP.LT.ZERO) CALL WRNDIE(-5,'<REPDSTR>', &
                  'replica EXCHange needs interval between temperatures.')
             DO I=1,NREPDSTR
                TEMPRX(I)=STEMP+(I-1)*DTEMP
                IF(QFASTREPDSTR) CURTEMPS(I)=TEMPRX(I)
             ENDDO
          ENDIF
#if KEY_OPENMM==1
       else if(ttemp<=0) then
#else
       ELSE
#endif
          DO I=1,NREPDSTR
             TEMPRX(I)=GTRMF(COMLYN,COMLEN,'TEMP',-ONE)
             IF(TEMPRX(I).LT.ZERO) THEN
                write(outu,'(a,i4,a,i4,a)')'Number of specifed temperatures ',i, &
                     ', but ',nrepdstr,' are needed'
                CALL WRNDIE(-5,'<REPDSTR>', &
                     'replica EXCHange needs all temperatures.')
             ENDIF
             IF(QFASTREPDSTR) CURTEMPS(I)=TEMPRX(I)
          ENDDO
       ENDIF
#if KEY_OPENMM==1
       if(qtor_repex) then
          if(ttemp <=0) call wrndie(-5,'<REPDSTR>', &
               'Torsion space replica exchange needs varible TTEM set.')
          DO I=1,NREPDSTR
             TEMPRX(I)=GTRMF(COMLYN,COMLEN,'PHIS',-ONE)
             IF(TEMPRX(I)<=0) THEN
                write(outu,'(a,i4,a,i4,a)')'Number of specifed torsion scaling ',i, &
                     ', but ',nrepdstr,' are needed'
                CALL WRNDIE(-5,'<REPDSTR>', &
                     'replica EXCHange needs all torsion scalings (PHIS).')
             ENDIF
             IF(QFASTREPDSTR) CURTEMPS(I)=TEMPRX(I)
          enddo
            temprx(1:nrepdstr) = ttemp
        if(prnlev>=2) then
             write(outu,'(a,i4,a)') &
                  ' Using OpenMM to run Replica Exchance in torsion space using ', nrepdstr,' replicas'
             write(outu,'(a,f6.2)') ' Dynamics will be run at single temperature ', ttemp
             write(outu,'(a,20(f6.2))') ' Each replica dihedral potentials scaled by', (curtemps(i), i=1, nrepdstr)
          endif
       endif
#endif
       IF(QRXSGLD)THEN
          IF(QFASTREPDSTR) &
             CALL WRNDIE(-4,'<REPDSTR>','FAST SGLD REX IS NOT SUPPORTED')

          allocate(sgtemprx(nrepdstr))
          allocate(sgftrx(nrepdstr))
          SGTEMP=GTRMF(COMLYN,COMLEN,'SGTE',ZERO)
          SGFT=GTRMF(COMLYN,COMLEN,'SGFT',ZERO)
          DSGFT=GTRMF(COMLYN,COMLEN,'DSGF',ZERO)
          DO I=1,NREPDSTR
                SGFTRX(I)=SGFT+(I-1)*DSGFT
          ENDDO
          MSGTEMP=GTRMF(COMLYN,COMLEN,'MSGT',-ONE)
          IF(SGTEMP>ZERO)THEN
            IF(MSGTEMP>ZERO)THEN
              DSGTEMP=EXP(DLOG(MSGTEMP/SGTEMP)/(NREPDSTR-1))
              DO I=1,NREPDSTR
                SGTEMPRX(I)=SGTEMP*DSGTEMP**(I-1)
              ENDDO
            ELSE 
              DSGTEMP=GTRMF(COMLYN,COMLEN,'DSGT',ZERO)
              DO I=1,NREPDSTR
                SGTEMPRX(I)=SGTEMP+(I-1)*DSGTEMP
              ENDDO
            ENDIF
          ELSE
            DO I=1,NREPDSTR
              SGTEMPRX(I)=GTRMF(COMLYN,COMLEN,'SGTT',-ONE)
              IF(SGTEMPRX(I).LT.ZERO) THEN
                write(outu,'(a,i4,a,f6.2)')'No SGTT input on stage ',i, &
                     ', set to simulation temperature: ',TEMPRX(I)
                SGTEMPRX(I)=TEMPRX(I)
              ENDIF
            ENDDO
          ENDIF
       ENDIF
       !
       !     TIGER parameters parsed here:
       !
       IF(QRXTIGER)THEN
          TIGERGR=GTRMF(COMLYN,COMLEN,'TOLG',ZERO)
          TIGERIT=GTRMI(COMLYN,COMLEN,'ITER',1)
          TIGERNM=GTRMI(COMLYN,COMLEN,'NMIN',100)
          TIGERNEQ=GTRMI(COMLYN,COMLEN,'NEQU',1000)
       ENDIF
    ENDIF

    IF(QREXCHGL)THEN
       IF(QFASTREPDSTR) &
          CALL WRNDIE(-4,'<REPDSTR>','FAST LAMBDA REX IS NOT SUPPORTED')

       QSUMP=(INDXA(COMLYN,COMLEN,'SUMP').GT.0)
       EWRITU=GTRMI(COMLYN,COMLEN,'UEWR',IUNREX)

       REPSEED=123+MYNODG
       IUNREX=GTRMI(COMLYN,COMLEN,'UNIT',OUTU)
       IREXFQ=GTRMI(COMLYN,COMLEN,'FREQ',1)
       !!IF(PRNLEV.GT.2) WRITE(OUTU,'(A,I5)') 'TIM DEBUG> IREXFQ = ',IREXFQ

       IF((.NOT.QEXPT).AND.(.NOT.QEX2D)) THEN
          !     repd nrep <int> exlm freq <int>
#if KEY_PHMD==1
          QPHRX=(INDXA(COMLYN,COMLEN,'PHMD').GT.0) !JAW. Flag for using PHMD and PH exchange
#endif
#if KEY_BLOCK==1
          QMSPHRX=(INDXA(COMLYN,COMLEN,'MSPH').GT.0) !GG: Flag for using pH-REX in CPHMD^MSLD
          IF(QMSPHRX)THEN
             !GG: Usage example "repd nrep <int> exlm freq <int> msph sph <int> mph <int>"
             IF (NREPDSTR .EQ. NREPLICA) THEN        !GG: Check number of replicas are the same in BLOCK and REPDSTR
                call msld_phrex(comlyn,comlen,nrepdstr,1,nblock)
             ELSE
                CALL WRNDIE(-5,'<REPDSTR>', &
                'Number of replicas declared in BLOCK and REPDSTR do not match!')
             ENDIF        !GG: For "NREPDSTR .EQ. NREPLICA" loop
          ENDIF           !GG: For "QMSPHRX" loop
#endif

          IF(QEXBK)IRBK=GTRMI(COMLYN,COMLEN,'REBK',1)
       ELSE IF (QEXPT) THEN
          !     repd nrep <int> exlm expt nrpt <int> freq <int>
          NREPT=GTRMI(COMLYN,COMLEN,'NRPT',1)
          IF(QEXBK)IRBK=GTRMI(COMLYN,COMLEN,'REBK',1)
       ELSE
          ! repd nrep <int> exlm ex2d nrpx <int> freq <int>
          NREPX=GTRMI(COMLYN,COMLEN,'NRPX',1)
          IF(QEXBK)IRBK=GTRMI(COMLYN,COMLEN,'REBK',1)
       ENDIF
    ENDIF
    
    ! initialize the starting temperature of this replica
    IF(QFASTREPDSTR) THEN
       TEMPCURRENT=CURTEMPS(IREPDSTR+1)
#if KEY_OPENMM==1
       if(qtor_repex) torsion_lambda = tempcurrent
#endif
       IF(PRNLEV.GE.6) &
          WRITE(OUTU,'(A,I3,A,F10.3)') 'FREX DEBUG> REPL ', IREPDSTR, ' INIT TEMPCURRENT = ', TEMPCURRENT
    ENDIF
    REPTAG=IREPDSTR
    NOPPUP=0
    NOPPUP2=0
    NOPPDN=0
    NOPPDN2=0
    NSUCUP=0
    NSUCUP2=0
    NSUCDN=0
    NSUCDN2=0
    REPROB=0.0
    REPROB2=0.0

    ! Initialize in case someone decides to try to use these"
    call set_param('REPROB',REPROB)
    call set_param('REPROB2',REPROB2)
    call set_param('EXRUP',ZERO)
    call set_param('EXRUP2',ZERO)
    call set_param('EXRDN',ZERO)
    call set_param('EXRDN2',ZERO)

    CALL FLUSH(OUTU)
    mypid = getpid()
    CALL PSYNC()

#else /* pll */
    CALL WRNDIE(-5,'<REPDSTR>', &
         'REPlica DiSTRibuted runs only in parallel.')
#endif /* pll */
    RETURN
  END SUBROUTINE REPDSTRMAIN
  !

  !
  SUBROUTINE SETUP_2D(COMLYN,COMLEN)
     use repdstr
     use parallel
     use string
     use stream
     use number
     use param_store, only: set_param

     ! Passed variables
     character(len=*)      :: comlyn
     integer               :: comlen

     ! Local variables
     integer               :: i,j
     character(len=4)      :: scratch

     ndim1=gtrmi(comlyn,comlen,'DIM1',0)
     ndim2=gtrmi(comlyn,comlen,'DIM2',0)
     d1freq=gtrmi(comlyn,comlen,'D1FR',-1)
     d2freq=gtrmi(comlyn,comlen,'D2FR',-1)
     iunrex=gtrmi(comlyn,comlen,'UNIT',-1)

     irexd1=0
     irexd2=0
     irexfq=min(d1freq,d2freq)
     qrepdstr=.true.
     qrexchg=.true.

     if(ndim1 <= 0 .or. ndim2 <= 0) call wrndie(-5,'<SETUP_2D>','DIMENSIONS CANNOT BE NEGATIVE!')
     if(d1freq < 1 .or. d2freq < 1) call wrndie(-5,'<SETUP_2D>','2D REX REQUIRES POSITIVE EXCHANGE FREQUENCIES')

     scratch=gtrma(comlyn,comlen,'D1CR')
     if(scratch == 'TEMP') then
        q2dtemp(1)=.true.
        q2dham(1)=.false.
        q2dph(1)=.false.
        q2ditemp=(indxa(comlyn,comlen,'ITEM').gt.0)
     else if(scratch == 'HAM') then
        q2dtemp(1)=.false.
        q2dham(1)=.true.
        q2dph(1)=.false.
     else if(scratch == 'PH') then
        q2dtemp(1)=.false.
        q2dham(1)=.false.
        q2dph(1)=.true.
     else
        call wrndie(-5,'<SETUP_2D>','UNRECOGNIZED EXCHANGE CRITERIA FOR DIMENSION 1.')
     endif

     scratch=gtrma(comlyn,comlen,'D2CR')
     if(scratch == 'TEMP') then
        q2dtemp(2)=.true.
        q2dham(2)=.false.
        q2dph(2)=.false.
        q2ditemp=(indxa(comlyn,comlen,'ITEM').gt.0)
     else if(scratch == 'HAM') then
        q2dtemp(2)=.false.
        q2dham(2)=.true.
        q2dph(2)=.false.
     else if(scratch == 'PH') then
        q2dtemp(2)=.false.
        q2dham(2)=.false.
        q2dph(2)=.true.
     else
        call wrndie(-5,'<SETUP_2D>','UNRECOGNIZED EXCHANGE CRITERIA FOR DIMENSION 2.')
     endif

     if(q2dtemp(1).and.q2dtemp(2)) call wrndie(-5,'<SETUP_2D>','ONLY ONE DIMENSION CAN USE TEMPERATURE EXCHANGE!')
     if(q2dph(1).and.q2dph(2)) call wrndie(-5,'<SETUP_2D>','ONLY ONE DIMENSION CAN USE PH EXCHANGE!')

     do i=1,ndim1
        if(q2dtemp(1)) then
           temprx(i)=gtrmf(comlyn,comlen,'TEMP',-1.0)
        else if(q2dph(1)) then
           phrx(i)=gtrmf(comlyn,comlen,'PHVA',-1.0)
        endif
     enddo
     do i=1,ndim2
        if(q2dtemp(2)) then
           temprx(i)=gtrmf(comlyn,comlen,'TEMP',-1.0)
        else if(q2dph(2)) then
           phrx(i)=gtrmf(comlyn,comlen,'PH',-1.0)
        endif
     enddo

     ! redistribute parallel set-up
     nrepdstr=ndim1*ndim2
     numnod=numnodg/nrepdstr

     if(prnlev.gt.2) &
        write(outu,'(a,3i4)') 'SETUP_2D> Number of processors, NREPDSTR, and NUMNOD = ',numnodg,nrepdstr,numnod
     if(numnod*nrepdstr /= numnodg) call wrndie(-5,'<SETUP_2D>','WRONG COMBINATION OF REPLICAS AND PROCESSORS')

     call psetloc
     call drepsetio(iolev,prnlev,wrnlev)

     irepdstr=mynodg/numnod
     repdid=irepdstr
     repseed=irepdstr+123

     myrepd1=mod(irepdstr,ndim1)
     myrepd2=irepdstr/ndim1 

     if(myrepd1 < ndim1-1) then
        nbrup1=irepdstr+1
        cnbrup1=mynodg+numnod
     else
        nbrup1=-1
        cnbrup1=-1
     endif
     if(myrepd1 > 0) then
        nbrdn1=irepdstr-1
        cnbrdn1=mynodg-numnod
     else
        nbrdn1=-1
        cnbrdn1=-1
     endif

     if(myrepd2 < ndim2-1) then
        nbrup2=irepdstr+ndim1
        cnbrup2=mynodg+(ndim1*numnod)
     else
        nbrup2=-1
        cnbrup2=-1
     endif
     if(myrepd2 > 0) then
        nbrdn2=irepdstr-ndim1
        cnbrdn2=mynodg-(ndim1*numnod)
     else 
        nbrdn2=-1
        cnbrdn2=-1
     endif

     call set_param('MYREP',irepdstr)
     call set_param('MYREPD1',myrepd1)
     call set_param('MYREPD2',myrepd2)
     call set_param('NREP',nrepdstr)
     call set_param('NREPD1',ndim1)
     call set_param('NREPD2',ndim2)

     call set_param('REPROB',zero)
     call set_param('REPROB2',zero)
     call set_param('EXRUP',zero)
     call set_param('EXRUP2',zero)
     call set_param('EXRDN',zero)
     call set_param('EXRDN2',zero)

     if(mynod.eq.0) then
        write(outu,'(a,8i5)') 'SETUP_2D> MYNODG,MYNOD,IREPDSTR,MYREPD1,NBRUP1,NBRDN1,CNBRUP1,CNBRDN1 = ', &
                              mynodg,mynod,irepdstr,myrepd1,nbrup1,nbrdn1,cnbrup1,cnbrdn1
        write(outu,'(a,8i5)') 'SETUP_2D> MYNODG,MYNOD,IREPDSTR,MYREPD2,NBRUP2,NBRDN2,CNBRUP2,CNBRDN2 = ', &
                              mynodg,mynod,irepdstr,myrepd2,nbrup2,nbrdn2,cnbrup2,cnbrdn2
     endif

  END SUBROUTINE SETUP_2D

  !
  SUBROUTINE DO2DEXCH(X,Y,Z,WMAIN,VX,VY,VZ,XOLD,YOLD,ZOLD,MYEPOT,TTEMP,ISTART,JHSTRT, &
                     ISEED,IASVEL,IGVOPT,CALLSEQ &
#if KEY_TSM==1
                     ,BACKLS & 
#endif
                     ,IDIDPHREX & 
                    )
    
    use psf
    use number
    use stream
    use parallel
    use repdstr
    use consta
    use memory
    use image
    use bases_fcm
    use energym, only: energy,eprop,epot
    use deriv,only: dx,dy,dz
    use consph,only: tstate
    use clcg_mod,only: random
    use imgup,only: upimag
    use param_store, only: set_param, get_param

    !
    ! passed-in variables
    real(chm_real) :: X(:), Y(:), Z(:),MYEPOT,VX(*),VY(*),VZ(*),WMAIN(*),TTEMP
    real(chm_real) XOLD(*),YOLD(*),ZOLD(*)
    INTEGER ISTART,JHSTRT,ISEED,IASVEL,IGVOPT,CALLSEQ
#if KEY_TSM==1
    INTEGER BACKLS(*) 
#endif

    !
    ! Local variables
    real(chm_real),allocatable,dimension(:)     :: w
    real(chm_real),allocatable,dimension(:,:,:) :: transf
    real(chm_real),dimension(6)                 :: oldxtlabc
    real(chm_real)                              :: eneigh,p,rn,nbrtemp,scalef,ourtemp
    real(chm_real)                              :: epot1p,epot1q,epot2p,epot2q,myexr
    integer                                     :: exchdim,me,neighbor,cneighbor,step,i
    logical                                     :: qexc,qcrys,qdim1,qdim2,qup,qdn,lexattempt

    integer                           :: iproto, jproto, ididphrex, phresstruct
    integer, allocatable,dimension(:) :: nstate
    real(chm_real)                    :: ph_l, ph_m, ph_delta

    ididphrex = 0

    qdim1 = .false.
    qdim2 = .false.
    qup = .false.
    qdn = .false.
    qcrys = (xtltyp /= '    ')
    lexattempt = .false.

    if(mod(istart-1,d2freq).eq.0) then
       qdim2 = .true.
       irexd2=irexd2+1

       exchdim=2
       step=mod(irexd2,2)
       me=myrepd2
       if(step==1) then
          if(mod(me,2) == 0) then
             qup=.true.
             neighbor=nbrup2
             cneighbor=cnbrup2
          else
             qdn=.true.
             neighbor=nbrdn2
             cneighbor=cnbrdn2
          endif
       else
          if(mod(me,2) == 0) then
             qdn=.true.
             neighbor=nbrdn2
             cneighbor=cnbrdn2
          else
             qup=.true.
             neighbor=nbrup2
             cneighbor=cnbrup2
          endif
       endif
    else if(mod(istart-1,d1freq).eq.0) then
       qdim1 = .true.
       irexd1=irexd1+1

       exchdim=1
       step=mod(irexd1,2)
       me=myrepd1
       if(step==1) then
          if(mod(me,2) == 0) then
             qup=.true.
             neighbor=nbrup1
             cneighbor=cnbrup1
          else
             qdn=.true.
             neighbor=nbrdn1
             cneighbor=cnbrdn1
          endif
       else
          if(mod(me,2) == 0) then
             qdn=.true.
             neighbor=nbrdn1
             cneighbor=cnbrdn1
          else
             qup=.true.
             neighbor=nbrup1
             cneighbor=cnbrup1
          endif
       endif
    else
       return
    endif

    if(cneighbor > -1) then
       lexattempt = .true.
       if(mynod == 0) then
          qexc = .false.
          write(outu,'(a,i3,a,i6,a,i5,a,i5)') 'DO2DEXCH> REPLICA ',IREPDSTR,' STEP ',istart-1, &
                                              ' EXCHANGE WITH NEIGHBOR ',neighbor, &
                                              ' PROCESSOR ',cneighbor

          call chmalloc('repdstr.src','DO2DEXCH','w',natom,crl=w)

          if(q2dtemp(exchdim)) then
             write(outu,'(a)') 'DO2DEXCH> DO TEMP EXCH'

             ! exchange temperature with neighbor and calculate if exchange succeeded
             call grecsen(cneighbor,1,eneigh,1,myepot,1)
             if(q2ditemp) then
                ourtemp=ttemp
             else
                ourtemp=temprx(me+1)
             endif
             call grecsen(cneighbor,4,nbrtemp,1,ourtemp,1)
             p=exp(-(one/(kboltz*ourtemp))-(one/(kboltz*nbrtemp))*(eneigh-myepot))
             p=min(p,one)
             rn=random(repseed)
             qexc=rn.lt.p

             if(irepdstr < neighbor) then
                ! we control the exchange
                call gsen(cneighbor,3,qexc,4)
                call gsen(cneighbor,4,p,8)
                call gsen(cneighbor,5,rn,8)
             else
                ! get our info from our natural superior
                call grec(cneighbor,3,qexc,4)
                call grec(cneighbor,4,p,8)
                call grec(cneighbor,5,rn,8)
             endif

             if(qexc) then
                call grecsen(cneighbor,2,w,natom,x,natom)
                x(1:natom)=w(1:natom)
                call grecsen(cneighbor,2,w,natom,y,natom)
                y(1:natom)=w(1:natom)
                call grecsen(cneighbor,2,w,natom,z,natom)
                z(1:natom)=w(1:natom)
                call grecsen(cneighbor,2,w,natom,vx,natom)
                vx(1:natom)=w(1:natom)
                call grecsen(cneighbor,2,w,natom,vy,natom)
                vy(1:natom)=w(1:natom)
                call grecsen(cneighbor,2,w,natom,vz,natom)
                vz(1:natom)=w(1:natom)

                scalef=sqrt(ttemp/nbrtemp)
                do i=1,natom
                   vx(i)=vx(i)*scalef
                   vy(i)=vy(i)*scalef
                   vz(i)=vz(i)*scalef
                enddo

                if(qcrys.and.xdim.gt.0) then
                   call grecsen(cneighbor,4,W,6,XTLABC,6)
                   xtlabc(1:6) = w(1:6)
                endif
                call grecsen(cneighbor,11,w,natom,wmain,natom)
                wmain(1:natom) = w(1:natom)

             endif
             write(iunrex,'(a,i10,a,f7.2,a,i3,a,f7.2,a,f8.2,a,f8.2,a,f5.3,a,f5.3,a,l1)') &
                   'T-REX>',istart-1,' TEMP ',ourtemp,' NBR ',neighbor,' NBRTMP ', &
                   nbrtemp, ' OURENE ',myepot,' NBRENE ',eneigh,' P ',p,' RN ', &
                   rn,' SUC? ',qexc

          else if(q2dham(exchdim)) then

             if(prnlev > 3) write(outu,'(a)') 'DO2DEXCH> DO HAMIL EXCH'
             qexc = .false.

             epot1p = myepot
             call grecsen(cneighbor,1,epot2q,1,epot1p,1)

             ! Do test coordinate exchange --actual work is
             ! done down below.
             call grecsen(cneighbor,2,w,natom,x,natom)
             x(1:natom)=w(1:natom)
             call grecsen(cneighbor,2,w,natom,y,natom)
             y(1:natom)=w(1:natom)
             call grecsen(cneighbor,2,w,natom,z,natom)
             z(1:natom)=w(1:natom)

             if(qcrys.and.xdim.gt.0) then
                call grecsen(cneighbor,4,W,6,XTLABC,6)
                xtlabc(1:6) = w(1:6)
             endif
             call grecsen(cneighbor,11,w,natom,wmain,natom)
             wmain(1:natom) = w(1:natom)
             if(prnlev > 3) write(outu,'(a)') 'DO2DEXCH> DONE PART A'

          else if(q2dph(exchdim)) then

             if(prnlev > 3) write(outu,'(a)') 'DO2DEXCH> DO PH EXCH'
             ididphrex = 1
             qexc = .false.
             ph_l = phrx(me+1)

             do i=1,nres
                if(tstate(i).eq.1) iproto = iproto + 1
             enddo
             call grecsen(cneighbor,6,jproto,1,iproto,1)
             call grecsen(cneighbor,7,ph_m,1,ph_l,1)

             ph_delta = log(10.0)*(ph_m - ph_l)*(iproto - jproto)
             if(ph_delta.le.zero) then
                p=one
             else
                p=min(one,exp(-ph_delta))
             endif
             rn=random(repseed)
             qexc=rn.lt.p

             if(irepdstr < neighbor) then
                ! we control the exchange
                call gsen(cneighbor,3,qexc,4)
                call gsen(cneighbor,4,p,8)
                call gsen(cneighbor,5,rn,8)
             else
                ! get our info from our natural superior
                call grec(cneighbor,3,qexc,4)
                call grec(cneighbor,4,p,8)
                call grec(cneighbor,5,rn,8)
             endif

             if(qexc) then
                call chmalloc('repdstr.src','DO2DEXCH','nstate',nres,intg=nstate)

                call grecsen(cneighbor,2,w,natom,x,natom)
                x(1:natom)=w(1:natom)
                call grecsen(cneighbor,2,w,natom,y,natom)
                y(1:natom)=w(1:natom)
                call grecsen(cneighbor,2,w,natom,z,natom)
                z(1:natom)=w(1:natom)
                call grecsen(cneighbor,2,w,natom,vx,natom)
                vx(1:natom)=w(1:natom)
                call grecsen(cneighbor,2,w,natom,vy,natom)
                vy(1:natom)=w(1:natom)
                call grecsen(cneighbor,2,w,natom,vz,natom)
                vz(1:natom)=w(1:natom)

                if(qcrys.and.xdim.gt.0) then
                   call grecsen(cneighbor,4,W,6,XTLABC,6)
                   xtlabc(1:6) = w(1:6)
                endif
                call grecsen(cneighbor,11,w,natom,wmain,natom)
                wmain(1:natom) = w(1:natom)

                call grecsen(cneighbor,14,w,natom,cg,natom)
                cg(1:natom)=w(1:natom)
                call grecsen(cneighbor,15,nstate,nres,tstate,nres)
                tstate(1:natom)=nstate(1:natom)

                call chmdealloc('repdstr.src','DO2DEXCH','nstate',nres,intg=nstate)
             endif
             write(iunrex,'(a,i10,a,f5.2,a,i3,a,f5.2,a,f5.3,a,f5.3,a,l1)') &
                   'PHREX>',istart-1,' PH ',ph_m,' NBR ',neighbor,' NBRPH ', &
                   ph_l,' P ',p,' RN ',rn,' SUC? ',qexc

          endif
          call chmdealloc('repdstr.src','DO2DEXCH','w',natom,crl=w)

       endif

       ! this part is executed on all nodes in da group

       call psnd4(qexc,1)
       call psnd4(exchdim,1)

#if KEY_CONSPH==1
       if(q2dph(exchdim)) then
          call psnd4(tstate,natom)
          call psnd8(cg,natom)
       endif
#endif

       if(qexc .or. q2dham(exchdim)) then
          call psnd8(x,natom)
          call psnd8(y,natom)
          call psnd8(z,natom)
          call psnd8(wmain,natom)
          if(qcrys) then
             oldxtlabc(1:6)=xtlabc(1:6)
             call psnd8(xtlabc,6)
             call xtllat(xucell,xtlabc)
             call chmalloc('repdstr.src','REPXCHG','TRANSF',3,4,XNSYMM,crl=transf)
             call imfill(transf,.false.)
             call chmdealloc('repdstr.src','REPXCHGL','TRANSF',3,4,XNSYMM,crl=transf)

             ! if we've done pH replica exchange, the charges might have changed.
             if(q2dph(exchdim)) call psnd8(cg,natom)

             call upimag(x,y,z,wmain,0,x,y,z,vx,vy,vz)
          endif
          call nbonds(x,y,z,bnbnd,bimag)
          call energy(x,y,z,dx,dy,dz,bnbnd,bimag,0)
       endif

       if(q2dham(exchdim)) then

          if(prnlev > 3) write(outu,'(a)') 'DO2DEXCH> DO PART B'

          ! This is unpleasant, because we have to do actual work here
          call psnd8(epot1p,1)
          call psnd8(epot2q,1)
          epot1q = eprop(epot) ! we calculated energy above
          
          ! Decide if any of this crap worked.
          if(mynod == 0) then
             call grecsen(cneighbor,5,epot2p,1,epot1q,1)

             if(irepdstr < neighbor) then
                p=exp(-(one/(kboltz*ttemp))*(epot2p+epot1q-epot2q-epot1p))
                p=min(p,one)
                rn=random(repseed)
                qexc=rn.lt.p

                call gsen(cneighbor,1,rn,8)
                call gsen(cneighbor,2,p,8)
                call gsen(cneighbor,3,qexc,4)
             else
                call grec(cneighbor,1,rn,8)
                call grec(cneighbor,2,p,8)
                call grec(cneighbor,3,qexc,4)
             endif

             if(prnlev > 3) write(outu,'(a)') 'DO2DEXCH> DONE PART B'
          endif

          call psnd4(qexc,1)
          if(qexc) then
             if(prnlev > 3) then
                write(outu,'(a)') 'DO2DEXCH> START PART C1'
                call flush(outu)
             endif
             call chmalloc('repdstr.src','DO2DEXCH','w',natom,crl=w)
             call grecsen(cneighbor,2,w,natom,vx,natom)
             vx(1:natom)=w(1:natom)
             call grecsen(cneighbor,2,w,natom,vy,natom)
             vy(1:natom)=w(1:natom)
             call grecsen(cneighbor,2,w,natom,vz,natom)
             vz(1:natom)=w(1:natom)

             call psnd8(vx,natom)
             call psnd8(vy,natom)
             call psnd8(vz,natom)
             call chmdealloc('repdstr.src','DO2DEXCH','w',natom,crl=w)
             if(prnlev > 3) write(outu,'(a)') 'DO2DEXCH> DONE PART C1'
          else
          
             ! FML: we have to undo all of the hard work we just did
             if(prnlev > 3) then 
                write(outu,'(a)') 'DO2DEXCH> START PART C2'
                call flush(outu)
             endif
             if(mynod == 0) then
                call chmalloc('repdstr.src','DO2DEXCH','w',natom,crl=w)
                call grecsen(cneighbor,2,w,natom,x,natom)
                x(1:natom)=w(1:natom)
                call grecsen(cneighbor,2,w,natom,y,natom)
                y(1:natom)=w(1:natom)
                call grecsen(cneighbor,2,w,natom,z,natom)
                z(1:natom)=w(1:natom)

                if(qcrys.and.xdim.gt.0) then
                   call grecsen(cneighbor,4,W,6,XTLABC,6)
                   xtlabc(1:6) = w(1:6)
                endif
                call grecsen(cneighbor,11,w,natom,wmain,natom)
                wmain(1:natom) = w(1:natom)
                call chmdealloc('repdstr.src','DO2DEXCH','w',natom,crl=w)
             endif

             call psnd8(x,natom)
             call psnd8(y,natom)
             call psnd8(z,natom)
             call psnd8(wmain,natom)
             if(qcrys) then
                oldxtlabc(1:6)=xtlabc(1:6)
                call psnd8(xtlabc,6)
                call xtllat(xucell,xtlabc)
                call chmalloc('repdstr.src','REPXCHG','TRANSF',3,4,XNSYMM,crl=transf)
                call imfill(transf,.false.)
                call chmdealloc('repdstr.src','REPXCHGL','TRANSF',3,4,XNSYMM,crl=transf)
                call upimag(x,y,z,wmain,0,x,y,z,vx,vy,vz)
             endif 
             call nbonds(x,y,z,bnbnd,bimag)
             call energy(x,y,z,dx,dy,dz,bnbnd,bimag,0)
             if(prnlev > 3) then
                write(outu,'(a)') 'DO2DEXCH> DONE PART C2'
                call flush(outu)
             endif

          endif

          write(iunrex,'(a,i10,a,f8.2,a,f8.2,a,f8.2,a,f8.2,a,f5.3,a,f5.3,a,l1)') &
                'H-REX>',istart-1,' EPOT1(P)',epot1p,' EPOT1(Q) ',epot1q,' EPOT2(P) ', &
                epot2p,' EPOT2(Q) ',epot2q,' P ',p,' RN ',rn,' SUC? ',qexc

       endif

       if(qexc) then
          ! Yay; our work is done! Just reset the dynamics with new velocities and go...

          jhstrt=0
          igvopt=2

          call psnd4(jhstrt,1)
          call psnd4(igvopt,1)
          call psnd8(vx,natom)
          call psnd8(vy,natom)
          call psnd8(vz,natom)
       endif

    else
       write(iunrex,'(a,i9,a)') '2DREX> ',istart-1,' SKIP EXCH'
    endif

    if(lexattempt) then

       ! statistics update time.

       if(qdim1.and.qdim2) call wrndie(-4,'<DO2DEXCH>','LOGIC BUG QDIM1 QDIM2')
       if(qup.and.qdn) call wrndie(-4,'<DO2DEXCH>','LOGIC BUG QUP QDN')

       if(qdim1.and.qup) then
          noppup=noppup+1
          if(qexc) nsucup=nsucup+1

          myexr=float(nsucup)/float(noppup)
          call set_param('EXRUP',myexr)
       else if(qdim2.and.qup) then
          noppup2=noppup2+1
          if(qexc) nsucup2=nsucup2+1

          myexr=float(nsucup2)/float(noppup2)
          call set_param('EXRUP2',myexr)
       else if(qdim1.and.qdn) then
          noppdn=noppdn+1
          if(qexc) nsucdn=nsucdn+1

          myexr=float(nsucdn)/float(noppdn)
          call set_param('EXRDN',myexr)
       else if(qdim2.and.qdn) then
          noppdn2=noppdn2+1
          if(qexc) nsucdn2=nsucdn2+1

          myexr=float(nsucdn2)/float(noppdn2)
          call set_param('EXRDN2',myexr)
       endif

       if((noppup+noppdn.gt.0).and.qdim1.and.lexattempt) then
          reprob=(min(p,one)+((noppup+noppdn-1)*reprob))/float(noppup+noppdn)
          call set_param('REPROB',reprob)
          if(prnlev > 2) then
             write(iunrex,'(a,4i3)') 'REPEXCH> noppup, noppdn, nsucup, nsucdn = ',noppup,noppdn,nsucup,nsucdn
             write(iunrex,'(a,f10.5)') 'REPEXCH> P = ',min(p,one)
             call get_param('EXRUP', exratup)
             call get_param('EXRDN', exratdn)
             write(iunrex,'(a,3f8.5)') 'REPEXCH> EXRUP EXRDN REPROB = ',exratup,exratdn,reprob
          endif
       endif
       if((noppup2+noppdn2.gt.0).and.qdim2.and.lexattempt) then
          reprob2=(min(p,one)+((noppup2+noppdn2-1)*reprob2))/float(noppup2+noppdn2)
          call set_param('REPROB2',reprob2)
          if(prnlev > 2) then
             write(iunrex,'(a,4i3)') 'REPEXCH> noppup2, noppdn2, nsucup2, nsucdn2 = ',noppup2,noppdn2,nsucup2,nsucdn2
             write(iunrex,'(a,f10.5)') 'REPEXCH> P = ',min(p,one)
             call get_param('EXRUP2', exratup2)
             call get_param('EXRDN2', exratdn2)
             write(iunrex,'(a,3f8.5)') 'REPEXCH> EXRUP2 EXRDN2 REPROB2 = ',exratup2,exratdn2,reprob2
          endif
       endif
    endif

  END SUBROUTINE DO2DEXCH
  
  !
  SUBROUTINE FASTREPEXCHG(WMAIN,EPOT,ISTEP,JHSTRT,IGVOPT,VX,VY,VZ,XOLD,YOLD,ZOLD &
#if KEY_TSM==1
                          ,BACKLS &
#endif
                         )
     use mpi
     use memory
     use parallel
     use repdstr 
     use stream
     use number
     use consta
     use clcg_mod,only: random

     ! Arguments
     REAL(CHM_REAL),INTENT(IN)    :: WMAIN(*),EPOT
     REAL(CHM_REAL),INTENT(INOUT) :: VX(*),VY(*),VZ(*),XOLD(*),YOLD(*),ZOLD(*)
     INTEGER,INTENT(IN)           :: ISTEP
     INTEGER,INTENT(OUT)          :: JHSTRT
     INTEGER,INTENT(INOUT)        :: IGVOPT
#if KEY_TSM==1
     INTEGER BACKLS(*)
#endif

     ! Local variables
     INTEGER                                 :: I,J,X,IERR,REP1,REP2,TMPUN,NPERREP
     INTEGER                                 :: OURSELVES,NEIGHBOR,REPORDER(MAXREPDSTR)
     LOGICAL                                 :: QEXC,QMASTER,QDOLOG
     REAL(CHM_REAL)                          :: EPOTARR(MAXREPDSTR),OURTEMP
     REAL(CHM_REAL)                          :: PROBARR(MAXREPDSTR)
     REAL(CHM_REAL)                          :: NBRTEMP,P,PROB,LOWTEMP,LASTLOW,REXP
     REAL(CHM_REAL),ALLOCATABLE,DIMENSION(:) :: SCARR

     IF(MOD(ISTEP-1,IREXFQ).NE.0) RETURN
     CALL PSETGLOB()
    
     IF(PRNLEV.GE.6) WRITE(OUTU,'(A,I6)') 'FREX DEBUG> IN FASTREPEXCHG AT STEP ', ISTEP
     IF(NREPDSTR.GT.MAXREPDSTR) &
        CALL WRNDIE(-5,'<FASTREPEXCHG>','TOO MANY REPLICAS FOR FAST EXCHANGING')
     DO I=1,MAXREPDSTR
        QEXC=.FALSE.
     ENDDO
     PROBARR(1:NREPDSTR)=ZERO

     NPERREP=NUMNOD/NREPDSTR
     CALL CHMALLOC('repdstr.src','FASTREPEXCHG','scarr',NUMNOD,crl=SCARR)

     CALL MPI_GATHER(EPOT,1,MPI_REAL8,SCARR,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
     IF(IERR.NE.MPI_SUCCESS) &
        CALL WRNDIE(-5,'<FASTREPEXCHG>','BUNGLED MPI COMMUNICATION')
     IF(MYNOD.EQ.0) THEN
        J=1
        DO I=1,NREPDSTR
           EPOTARR(I)=SCARR(J)
           J=J+NPERREP
        ENDDO

        IF(FFRUN) THEN
           WRITE(IUNREX,'(A)') '# replica temp. ener. neighbor ntemp nene prob p success? newrep'
           FFRUN=.FALSE.
        ENDIF

        ! We need to actually figure out how the exchanges are going to
        ! happen, to do so, apply the formula and temps to the energies.
        DO X=1,NREPEAT
           IF(X.EQ.1.OR.X.EQ.NREPEAT) THEN
              QDOLOG=.TRUE.
           ELSE IF(LOGLEVEL.EQ.0) THEN
              QDOLOG=.FALSE.
           ELSE IF(MOD(X,LOGLEVEL).EQ.0) THEN
              QDOLOG=.TRUE.
           ELSE
              QDOLOG=.FALSE.
           ENDIF

           IREX=IREX+1
           IF(QDOLOG) WRITE(IUNREX,'(A,I15,A,I12,A,I5)') '# Exchange ', IREX, ': STEP ', ISTEP-1, ': REPEAT ', X

           ! put the reservoirs in order of temperature ... a bit of a hacky
           ! bubble sort for now.
           LASTLOW=9999.0
           DO I=1,NREPDSTR
              LOWTEMP=9999.0
              DO J=1,NREPDSTR
                 IF(I.GT.1) THEN
                    IF(CURTEMPS(J).GT.LASTLOW) THEN
                       IF(CURTEMPS(J).LT.LOWTEMP) THEN
                          REPORDER(I)=J
                          LOWTEMP=CURTEMPS(J)
                       ENDIF
                    ENDIF
                 ELSE
                    IF(CURTEMPS(J).LT.LOWTEMP) THEN
                       REPORDER(I)=J
                       LOWTEMP=CURTEMPS(J)
                    ENDIF
                 ENDIF
              ENDDO
              LASTLOW=LOWTEMP
           ENDDO

           QMASTER=(MOD(IREX,2).EQ.0)
           DO I=1,NREPDSTR
              QEXC=.FALSE.
              OURSELVES=REPORDER(I)
              OURTEMP=CURTEMPS(OURSELVES)
              
              ! find our neighbor (next highest temperature replica)
              IF(QMASTER) THEN
                 IF(I.LE.NREPDSTR-1) THEN
                    NEIGHBOR=REPORDER(I+1)
                    NBRTEMP=CURTEMPS(NEIGHBOR)
                 ELSE 
                    NEIGHBOR=-ONE
                    NBRTEMP=9999.0
                 ENDIF
              ELSE
                 ! Special case to make sure that there's a print out for the first
                 ! replica.
                 IF(I.EQ.1) THEN
                    IF(QDOLOG) WRITE(IUNREX,'(I2,x,F12.6,x,F15.6,x,I2,x,F12.6,x,F15.6,x,F5.3,x,F5.3,x,L,x,I2)') &
                               OURSELVES,OURTEMP,EPOTARR(OURSELVES),-1,9999.0,0.00,0.00,0.00,.FALSE.,OURSELVES
                 ENDIF
                 QMASTER=.NOT.QMASTER
                 CYCLE
              ENDIF

              !!WRITE(IUNREX,'(A,I3,A,F10.3,A,I3,A,F10.3)') &
              !!      '# FREP DEBUG> OURSELVES = ', OURSELVES, ' OURTEMP = ', OURTEMP, &
              !!      ' NEIGHBOR = ', NEIGHBOR, ' NBRTEMP = ', NBRTEMP  

              IF(NEIGHBOR.GT.0) THEN
                 ! We control the exchange; our neighbor has the next
                 ! highest temperature.
                 IF(OURTEMP.EQ.NBRTEMP) THEN
                    PROB=ONE
                 ELSE
#if KEY_OPENMM==1
                    if(qtor_repex) then
                       prob = min(one, exp(-(one/(kboltz*temprx(1)))       &
                            * ( (nbrtemp/ourtemp-one)*epotarr(ourselves)   & 
                            +   (ourtemp/nbrtemp-one)*epotarr(neighbor) )  &
                            ))
                       p = -(one/(kboltz*temprx(1)))       &
                            * ( (nbrtemp/ourtemp-one)*epotarr(ourselves)   & 
                            +   (ourtemp/nbrtemp-one)*epotarr(neighbor) )  

                    else
#endif
                       PROB=MIN(ONE,EXP(-(ONE/(KBOLTZ*OURTEMP) &
                            -ONE/(KBOLTZ*NBRTEMP))*(EPOTARR(NEIGHBOR)-EPOTARR(OURSELVES))))
#if KEY_OPENMM==1
                    endif
#endif
                 ENDIF
                 P=RANDOM(REPSEED)
                
                 IF(P.LE.PROB) THEN
                    QEXC=.TRUE.
                    CURTEMPS(NEIGHBOR)=OURTEMP
                    CURTEMPS(OURSELVES)=NBRTEMP
                    REP1=NEIGHBOR
                    REP2=OURSELVES
                 ELSE
                    QEXC=.FALSE.
                    REP1=OURSELVES
                    REP2=NEIGHBOR
                 ENDIF
                 PROBARR(OURSELVES)=(((X-1)*PROBARR(OURSELVES))+MIN(PROB,ONE))/FLOAT(X)
                 PROBARR(NEIGHBOR)=(((X-1)*PROBARR(NEIGHBOR))+MIN(PROB,ONE))/FLOAT(X)
                 !!!WRITE(IUNREX,'(A,2I3)') 'UPDATING PROBARR OF REPLICAS ',OURSELVES,NEIGHBOR
                 IF(QDOLOG) WRITE(IUNREX,'(I2,x,F12.6,x,F15.6,x,I2,x,F12.6,x,F15.6,x,F5.3,x,F5.3,x,L,x,I2)') &
                            OURSELVES,OURTEMP,EPOTARR(OURSELVES),NEIGHBOR,NBRTEMP,EPOTARR(NEIGHBOR),PROB,P,QEXC,REP1

                 ! write out a line for the neighboring replica, as well
                 IF(QDOLOG) WRITE(IUNREX,'(I2,x,F12.6,x,F15.6,x,I2,x,F12.6,x,F15.6,x,F5.3,x,F5.3,x,L,x,I2)') &
                            NEIGHBOR,NBRTEMP,EPOTARR(NEIGHBOR),OURSELVES,OURTEMP,EPOTARR(OURSELVES),PROB,P,QEXC,REP2

              ELSE
                 ! Give this a 0 exchange probability for records-keeping purposes
                 PROBARR(OURSELVES)=((X-1)*PROBARR(OURSELVES))/FLOAT(X)
                 !!!WRITE(IUNREX,'(A,I2,A,F6.4)') 'PROBDEBUG> NO EXCHANGE FOR ',OURSELVES,' PROBARR = ',PROBARR(OURSELVES)

                 IF(QDOLOG) WRITE(IUNREX,'(I2,x,F12.6,x,F15.6,x,I2,x,F12.6,x,F15.6,x,F5.3,x,F5.3,x,L,x,I2)') &
                            OURSELVES,OURTEMP,EPOTARR(OURSELVES),-1,9999.0,0.00,0.00,0.00,.FALSE.,OURSELVES

              ENDIF ! END IF(OURSELVES.GT.0)
              QMASTER=.NOT.QMASTER
           ENDDO
        ENDDO

        ! ----
        ! Copy exchange probability data to a point where it can be blasted out
        ! one per node. Mike C's COMM_MASTER will remove the need for hijinks
        ! like this.
        ! ----
        WRITE(IUNREX,'(A)',ADVANCE='NO') 'PROBARR> '
        J=1
        DO I=1,NUMNOD
           WRITE(IUNREX,'(A,I1,A,F6.4,X)',ADVANCE='NO') 'PROBARR(', J, ') =', PROBARR(J)
           SCARR(I)=PROBARR(J)
           IF(MOD(I,NPERREP).EQ.0) J=J+1
        ENDDO
        WRITE(IUNREX,'(A)') ' '

     ENDIF ! END part only executed on processor 0

     ! Now that we have figured out the final temperatures at each state,
     ! broadcast curtemps and call dofastexchg to actually make the exchange.
     CALL MPI_BCAST(CURTEMPS,NREPDSTR,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
     CALL MPI_BCAST(IREX,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
     CALL MPI_BCAST(NREPEAT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
     CALL MPI_SCATTER(SCARR,1,MPI_REAL8,REXP,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

     CALL CHMDEALLOC('repdstr.src','FASTREPEXCHG','scarr',NUMNOD,crl=SCARR)
     CALL PSETLOC()
     CALL DOFASTEXCHG(CURTEMPS(IREPDSTR+1),JHSTRT,IGVOPT,VX,VY,VZ,XOLD,YOLD,ZOLD,REXP &
#if KEY_TSM==1
                      ,BACKLS &
#endif
                     )
     CALL FLUSH(OUTU)
  END SUBROUTINE FASTREPEXCHG

  SUBROUTINE DOFASTEXCHG(TEMPNEW,JHSTRT,IGVOPT,VX,VY,VZ,XOLD,YOLD,ZOLD,REXP &
#if KEY_TSM==1
                         ,BACKLS &
#endif
                        )
     use mpi
     use psf
     use parallel
     use stream
     use repdstr
     use coord,only: x,y,z
!     use omm_main, only : omm_change_lambda
  
     ! Arguments
     REAL(CHM_REAL),INTENT(IN)    :: TEMPNEW,REXP
     REAL(CHM_REAL),INTENT(INOUT) :: VX(*),VY(*),VZ(*),XOLD(*),YOLD(*),ZOLD(*)
     INTEGER,INTENT(OUT)          :: JHSTRT
     INTEGER,INTENT(INOUT)        :: IGVOPT
#if KEY_TSM==1
     INTEGER BACKLS(*) 
#endif    

     ! Local variables
     LOGICAL        :: QUPVELOC
     REAL(CHM_REAL) :: SCALEF
     INTEGER        :: I,IERR

     IF(MYNOD.EQ.0) THEN
        REPROB=(((IREX-NREPEAT)*REPROB)+(NREPEAT*REXP))/FLOAT(IREX)
        call set_param('REPROB',REPROB)

        ! In FAST replica exchange, we don't track things explicitly
        ! because all decisions are made on the main processor and it
        ! would be a bit costly to send them to each replica; so up and
        ! down are rather arbitrary in this context...
        IF(MOD(IREX,2)==0) THEN
           NOPPDN=NOPPDN+1
        ELSE
           NOPPUP=NOPPUP+1
        ENDIF
        IF(PRNLEV.GT.5) THEN
           WRITE(IUNREX,'(A,2I5)')  'FREX PROB> NOPPUP,NOPPDN = ',NOPPUP,NOPPDN
           WRITE(IUNREX,'(A,F8.6)') 'FREX PROB> REPROB = ',REPROB
        ENDIF

        IF(PRNLEV.GT.3) &
           WRITE(OUTU,'(A,2F10.3)') 'FREX DEBUG> FORMER AND NEW TEMPS: ', TEMPCURRENT, TEMPNEW
        IF(TEMPCURRENT.NE.TEMPNEW) THEN

           IF(TEMPCURRENT.GT.TEMPNEW) THEN
              NSUCUP=NSUCUP+1
              call set_param('EXRUP', REAL(NSUCUP, chm_real) / REAL(NOPPUP, chm_real))
           ELSE
              NSUCDN=NSUCDN+1
              call set_param('EXRDN', REAL(NSUCDN, chm_real) / REAL(NOPPDN, chm_real))
           ENDIF
           IF(PRNLEV.GT.5) THEN
              WRITE(IUNREX,'(A,2I5)')   'FREX PROB> NSUCUP,NSUCDN = ',NSUCUP,NSUCDN
              WRITE(IUNREX,'(A,2F9.6)') 'FREX PROB> EXRUP,EXRDN = ', &
                  REAL(NSUCUP, chm_real) / REAL(NOPPUP, chm_real), &
                  REAL(NSUCDN, chm_real) / REAL(NOPPDN, chm_real)
           ENDIF

#if KEY_OPENMM==1
           if(.not.qtor_repex) then
#endif
              QUPVELOC=.TRUE.
              SCALEF=SQRT(TEMPNEW/TEMPCURRENT)
              IF(PRNLEV.GE.6) WRITE(OUTU,'(A,F7.4)') 'FREX DEBUG> SCALEF: ', SCALEF
              DO I=1,NATOM
                 VX(I)=VX(I)*SCALEF
                 VY(I)=VY(I)*SCALEF
                 VZ(I)=VZ(I)*SCALEF
                 XOLD(I)=XOLD(I)*SCALEF
                 YOLD(I)=YOLD(I)*SCALEF
                 ZOLD(I)=ZOLD(I)*SCALEF
              ENDDO
              TEMPCURRENT=TEMPNEW
#if KEY_OPENMM==1
           else
              tempcurrent = tempnew
              qupveloc = .false.
           endif
#endif
        ELSE
           QUPVELOC=.FALSE.
        ENDIF
     ENDIF
 
     CALL PSND4(QUPVELOC,1)
     CALL PSND8(TEMPCURRENT,1)

#if KEY_OPENMM==1
     if(qtor_repex) then
        torsion_lambda = tempcurrent
!        call omm_change_lambda(tempcurrent)
     endif
#endif
     IF(QUPVELOC) THEN
        JHSTRT=0
        IGVOPT=2
        
        ! BTM: I am not sure if we still need to send xold, yold, and zold,
        ! but I am going to do it anyway to be safe.
        CALL PSND4(IGVOPT,1)
        CALL PSND4(JHSTRT,1)
        CALL PSND8(VX,NATOM)
        CALL PSND8(VY,NATOM)
        CALL PSND8(VZ,NATOM)
        CALL PSND8(XOLD,NATOM)
        CALL PSND8(YOLD,NATOM)
        CALL PSND8(ZOLD,NATOM)
     ENDIF               
  END SUBROUTINE DOFASTEXCHG

  SUBROUTINE REPEXCHG(X,Y,Z,WMAIN,VX,VY,VZ,XOLD,YOLD,ZOLD,MYEPOT,TTEMP,ISTART,JHSTRT, &
                      ISEED,IASVEL,IGVOPT,CALLSEQ &
#if KEY_TSM==1
                      ,BACKLS &
#endif
                      ,IDIDPHREX &
                     )

    !-----------------------------------------------------------------------
    !     Perform necessary communication and calculate new data...
    !
  use number
  use consta
  use comand
  use dimens_fcm
  use psf
  use stream
  use parallel
  use repdstr
  use memory
  use phmd !JAW
  use bases_fcm
  use deriv
  use image
  use energym, only: energy,eprop,eterm,lenent,epot
  use clcg_mod,only: random
  use cnst_fcm,only: fbeta
  use reawri,only: delta
  use imgup,only: upimag
  use sgld,only: EPOTLF,EPOTHF,AVGEFLF,AVGEFHF,AVGCFLF,AVGCFHF,AVGTLF,TREFLF,TRXLF, &
        SGVX,SGVY,SGVZ,SGFX,SGFY,SGFZ,SGGX,SGGY,SGGZ,SGHX,SGHY,SGHZ,SGKX,SGKY,SGKZ,QSGLD
  use consph,only: tstate
  use block_ltm
  use lambdam    !GG: MSLD-compatibility

    !
    real(chm_real) :: X(:), Y(:), Z(:),MYEPOT,VX(*),VY(*),VZ(*),WMAIN(*),TTEMP
    real(chm_real) XOLD(*),YOLD(*),ZOLD(*)

    INTEGER ISTART,JHSTRT,ISEED,IASVEL
#if KEY_TSM==1
    INTEGER BACKLS(*) 
#endif

    !
    LOGICAL QEXC,QCRYS,LUSED,LEXATTEMPT
    INTEGER STEP,IDECIDE,ME,NEIGHBOR,IGVOPT,I,J,CNEIGHBOR,FNEIGH,CALLSEQ
    real(chm_real) ENEIGH,SCALED,P,srate,rn,ttx,ttsgx,tfsgx,sgldarg,hrexarg
    REAL(CHM_REAL),allocatable,dimension(:) :: w,comar,dd,timx,timy,timz
    REAL(CHM_REAL),allocatable,dimension(:) :: oldwmain
    logical qhes, qdidrsvr
    real(chm_real) SGARRAY(10),SGARRAYN(10),oldxtlabc(6)
    real(chm_real) SCALSG,DTEFLF,DTCFLF,FACT
    real(chm_real) ecsum,ourpoti,ourpotj,nbrpoti,nbrpotj
    real(chm_real),allocatable,dimension(:,:,:) :: TRANSF
#if KEY_PHMD==1
    real(chm_real) l(ntitr) !JAW
#endif
#if KEY_BLOCK==1
    integer K                                      !GG: MSLD-compatibility
    real(chm_real) n(nsitemld,nblock)              !GG: MSLD-compatibility
    real(chm_real) m(nsitemld*nblock)              !GG: MSLD-compatibility
    real(chm_real) THETAVMLDS(nsitemld*nblock)     !GG: MSLD-compatibility
    real(chm_real) THETAMLDS(nsitemld*nblock)      !GG: MSLD-compatibility
    real(chm_real) THETAMLDOLDS(nsitemld*nblock)   !GG: MSLD-compatibility
    real(chm_real) THETAFMLDS(nsitemld*nblock)     !GG: MSLD-compatibility
#endif
    integer n6,ical,oldrep

    integer                           :: iproto, jproto, ididphrex, phresstruct
    integer, allocatable,dimension(:) :: nstate
    real(chm_real)                    :: ph_l, ph_m, ph_delta

    ididphrex = 0

    lexattempt=.false.
    oldrep=reptag

    IF(Q2DREX) THEN
       CALL DO2DEXCH(X,Y,Z,WMAIN,VX,VY,VZ,XOLD,YOLD,ZOLD,MYEPOT,TTEMP,ISTART,JHSTRT, &
                      ISEED,IASVEL,IGVOPT,CALLSEQ &
#if KEY_TSM==1
                      ,BACKLS &
#endif
                      ,IDIDPHREX & 
                     )
       RETURN
    ENDIF


    !
    !     When qpxtiger is on we are in a tiger pre-exchange state
    IF(QPXTIGER)THEN
       IF(MOD(ISTART-1,TIGERNEQ).NE.0) RETURN
       qrxtmin=.true.   ! we can do mini first now
    ELSE
       IF(MOD(ISTART-1,IREXFQ).NE.0) RETURN
    ENDIF
    QCRYS = (XTLTYP.NE.'    ')

    !
    !     If we are about to do the exchange in the case of TIGER method
    !     we need some preparation:
    !       1. run number of iteration cycles of
    !          the pair of minimizer and equlibration
    !       2. then try for the exchange
    !
    !
    IF(QRXTIGER)THEN
       IF(QRXTMIN)THEN
          write(comlyn, &
               '(''mini abnr nstep '',i6,'' tolg '',f12.8,'' nprint 100'')') &
               tigernm,tigergr
          comlen=51
          call maincomx(comlyn,comlen,lused)
          qrxtmin=.false.
       ENDIF
       IF(QPXTIGER)THEN
          TIGERITI=TIGERITI+1
          IF(TIGERITI.LT.TIGERIT)THEN
             return  ! maybe
          ELSE
             qpxtiger=.false.
             tigeriti=0
          ENDIF
       ENDIF
    ENDIF
    IF(QRXSGLD)THEN
       IF(TREFLF<RSMALL)THEN
          CALL PSETGLOB()
          CALL PSYNC()
          TRXLF=AVGTLF
          CALL PSND8(TRXLF,1)
          TRXLF=TRXLF*temprx(irepdstr+1)/temprx(1)
          CALL PSETLOC()
       ENDIF
    ENDIF
    !
    !     Prepare the variable NEIGHBOR.
    !     NEIGHBOR=-1 if no communication is needed on this process

    STEP=MOD(IREX,2)
    ME=MOD(IREPDSTR,2)

    IF(STEP.EQ.1)THEN
       NEIGHBOR=IREPDSTR+1
       IF(ME.NE.0)NEIGHBOR=IREPDSTR-1
    ELSE
       NEIGHBOR=IREPDSTR-1
       IF(ME.NE.0)NEIGHBOR=IREPDSTR+1
    ENDIF
    IF(NEIGHBOR.GE.NREPDSTR)NEIGHBOR=-1

    if(neighbor >= 0) then
       if(neighbor > irepdstr) then
          noppup=noppup+1
       else
          noppdn=noppdn+1
       endif
    endif

    eneigh=zero                ! printout looks better this way
    ourpoti=zero
    ourpotj=zero
    nbrpoti=zero
    nbrpotj=zero
    qdidrsvr=.false.

    QEXC=.FALSE.
    RN=RANDOM(REPSEED)
    scaled=one
    scalsg=one

    IF(NEIGHBOR.GE.0) THEN

      LEXATTEMPT=.TRUE.
      CNEIGHBOR=NEIGHBOR*NUMNOD
      IF(QRXSGLD.AND.MYNOD.EQ.0)THEN
          ! Let's take care of SGLD before we start mucking
          ! about with the energies in the H-REX code below.
          SGARRAY(1)=EPOTLF
          SGARRAY(2)=EPOTHF+EPOTLF
          SGARRAY(3)=(AVGEFLF*AVGCFLF-AVGEFHF*AVGCFHF)/(kboltz*temprx(irepdstr+1))
          SGARRAY(4)=AVGEFHF*AVGCFHF/(kboltz*temprx(irepdstr+1))
          SGARRAY(5)=REPDID
          SGARRAY(6)=AVGTLF
          SGARRAY(7)=AVGEFLF
          SGARRAY(8)=AVGCFLF
          CALL GRECSEN(CNEIGHBOR,1,SGARRAYN,8,SGARRAY,8)
          ENEIGH=SGARRAYN(2)
          !p=min(one,exp(-(SGARRAY(3)-SGARRAYN(3))*(SGARRAYN(1)-SGARRAY(1))  &
          ! -(SGARRAY(4)-SGARRAYN(4))*(SGARRAYN(2)-SGARRAY(2))  &
          ! -(SGARRAY(3)*SGARRAYN(7)-SGARRAYN(3)*SGARRAY(7))*(SGARRAY(6)-SGARRAYN(6))))
          !p=min(one,exp(-(SGARRAY(3)-SGARRAYN(3))*(SGARRAYN(1)-SGARRAY(1))  &
          ! -(SGARRAY(4)-SGARRAYN(4))*(SGARRAYN(2)-SGARRAY(2))))

          if(.not.QRXTHAM) then
             sgldarg = -((SGARRAY(3)-SGARRAYN(3))*(SGARRAYN(1)-SGARRAY(1)) &
                       -(SGARRAY(4)-SGARRAYN(4))*(SGARRAYN(2)-SGARRAY(2)))
          endif
      ELSE IF(QRXTHAM) THEN

          ! ouch -- perform a hamiltonian coordinate swap INSIDE temperature REX
          !
          ! We're doing this up here because we need to perform the test swap and
          ! energy eval BEFORE telling other processors to get lost.
          if(mynod.eq.0) then
             CALL GRECSEN(CNEIGHBOR,1,ENEIGH,1,MYEPOT,1)
             ourpoti=myepot
             nbrpotj=eneigh

             !     Perform a test coordinate exchange and calculate the new energies
             call chmalloc('repdstr.src','REPXCHG','W',NATOM,crl=w)
             call chmalloc('repdstr.src','REPXCHG','TIMX',NATOM,crl=timx)
             call chmalloc('repdstr.src','REPXCHG','TIMY',NATOM,crl=timy)
             call chmalloc('repdstr.src','REPXCHG','TIMZ',NATOM,crl=timz)
             call chmalloc('repdstr.src','REPXCHG','OLDWMAIN',NATOM,crl=oldwmain)

             timx(1:natom)=x(1:natom)
             timy(1:natom)=y(1:natom)
             timz(1:natom)=z(1:natom)
             oldwmain(1:natom)=wmain(1:natom)

             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,X,NATOM)
             x(1:natom) = w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,Y,NATOM)
             y(1:natom) = w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,Z,NATOM)
             z(1:natom) = w(1:natom)
             if(qcrys.and.xdim.gt.0) then
                call grecsen(cneighbor,4,W,6,XTLABC,6)
                XTLABC(1:6) = W(1:6)
             endif
             CALL GRECSEN(CNEIGHBOR,11,W,NATOM,WMAIN,NATOM)
             wmain(1:natom) = w(1:natom)
          endif

          call psnd8(x,natom)
          call psnd8(y,natom)
          call psnd8(z,natom)
          call psnd8(wmain,natom)

          if(qcrys) then
             oldxtlabc(1:6)=xtlabc(1:6)
             call psnd8(xtlabc,6)
             call xtllat(xucell,xtlabc)
             call chmalloc('repdstr.src','REPXCHG','TRANSF',3,4,XNSYMM,crl=transf)
             call imfill(transf,.false.)
             call chmdealloc('repdstr.src','REPXCHGL','TRANSF',3,4,XNSYMM,crl=transf)
             call upimag(x,y,z,wmain,0,x,y,z,vx,vy,vz)
          endif
          call nbonds(x,y,z,bnbnd,bimag)
          call energy(x,y,z,dx,dy,dz,bnbnd,bimag,0)
          ecsum=zero
          do i=1,lenent
             ecsum=eterm(i)+ecsum
          enddo
          eprop(epot)=ecsum

          if(mynod.eq.0) then

             CALL GRECSEN(CNEIGHBOR,1,ENEIGH,1,ECSUM,1)
             ourpotj=ecsum
             nbrpoti=eneigh

             ! calculate ham REX argument taking different temperatures
             ! into account, from Sugita, Kitao, and Okamoto, J. Chem.
             ! Phys. 113, 6042 (2000).

             write(outu,'(a,4f12.4)') 'DD> NBRPOTJ,NBRPOTI,OURPOTJ,OURPOTI = ',nbrpotj,nbrpoti,ourpotj,ourpoti
             write(outu,'(a,f11.4)')  'DD> k  = ',kboltz
             write(outu,'(a,f11.4)')  'DD> T (us)    = ',temprx(irepdstr+1)
             write(outu,'(a,f11.4)')  'DD> T (them)  = ',temprx(neighbor+1)
             write(outu,'(a,f11.4)')  'DD> 1/kT (us) = ',(one/(kboltz*temprx(irepdstr+1)))
             write(outu,'(a,f11.4)')  'DD> 1/kT (them) = ',(one/(kboltz*temprx(neighbor+1)))
             !hrexarg=-((one/(kboltz*temprx(irepdstr+1)))*(nbrpotj-nbrpoti) &
             !         -(one/(kboltz*temprx(neighbor+1)))*(ourpotj-ourpoti))
             !write(outu,'(a,f11.4)') 'hrexarg = ',hrexarg

             hrexarg=-((one/(kboltz*temprx(irepdstr+1)))*(ourpotj-ourpoti) &
                      -(one/(kboltz*temprx(neighbor+1)))*(nbrpotj-nbrpoti))

             ! We have no clue whether this will succeed or fail, so for now, just put
             ! things back the way that they were and let the exchange code down
             ! below take care of things.

             x(1:natom)=timx(1:natom)
             y(1:natom)=timy(1:natom)
             z(1:natom)=timz(1:natom)
             wmain(1:natom)=oldwmain(1:natom)

             call chmdealloc('repdstr.src','REPXCHG','W',NATOM,crl=w)
             call chmdealloc('repdstr.src','REPXCHG','TIMX',NATOM,crl=timx)
             call chmdealloc('repdstr.src','REPXCHG','TIMY',NATOM,crl=timy)
             call chmdealloc('repdstr.src','REPXCHG','TIMZ',NATOM,crl=timz)
             call chmdealloc('repdstr.src','REPXCHG','OLDWMAIN',NATOM,crl=oldwmain)

          endif
          call psnd8(x,natom)
          call psnd8(y,natom)
          call psnd8(z,natom)
          call psnd8(wmain,natom)

          if(qcrys) then
             xtlabc(1:6)=oldxtlabc(1:6)
             call xtllat(xucell,xtlabc)
             call chmalloc('repdstr.src','REPXCHG','TRANSF',3,4,XNSYMM,crl=transf)
             call imfill(transf,.false.)
             call chmdealloc('repdstr.src','REPXCHG','TRANSF',3,4,XNSYMM,crl=transf)
             call upimag(x,y,z,wmain,0,x,y,z,vx,vy,vz)
          endif
          call nbonds(x,y,z,bnbnd,bimag)
          call energy(x,y,z,dx,dy,dz,bnbnd,bimag,0)
          ecsum=zero
          do i=1,lenent
             ecsum=eterm(i)+ecsum
          enddo
          eprop(epot)=ecsum

      ENDIF
    ENDIF

    IF(MYNOD.NE.0) GOTO 99

#if KEY_BLOCK==1
    call msld_checkvariables(1) !GG Variables Before MC Exchange
#endif 
    !
    !     FIXME: Maybe for protection we can call psetloc() here
    !            This code relies on the fact that psetloc()
    !            values are in NUMNOD and MYNOD. It is important
    !            that every psetglob() is immediately followed
    !            by psetloc() in any part of the code!
    !
    !

    IF(NEIGHBOR.EQ.-1) THEN
       IRESEXCH=IRESEXCH+1
       IF(QRESERVOIR) THEN

          IF((IREPDSTR.EQ.0).AND.QRESLOW) THEN
             QDIDRSVR=.TRUE.
             LEXATTEMPT=.TRUE.
             NOPPDN=NOPPDN+1
             IF(PRNLEV.GE.6) &
                WRITE(IUNREX,*) 'REPEXCH> BOT REPL WOULD EXCH W/ RESERVOIR'
             IF(QPHREX) THEN
                IDIDPHREX=1 
                IPROTO=0
                DO I=1,NRES
                   IF(TSTATE(I).EQ.1) IPROTO=IPROTO+1
                ENDDO
                CALL RESEXCH_PH(TTEMP,.FALSE.,X,Y,Z,VX,VY,VZ,PHRX(IREPDSTR+1),IPROTO, &
                                ISEED,IASVEL,IGVOPT,P,RN,JPROTO,PH_M,QEXC,FNEIGH &
#if KEY_TSM==1
                                ,BACKLS &
#endif
                               )
             ELSE
                CALL RESEXCH(.FALSE.,X,Y,Z,VX,VY,VZ,QEXC,TEMPRX(IREPDSTR+1),MYEPOT, &
                             ISEED,IASVEL,IGVOPT,P,RN,ENEIGH,FNEIGH &
#if KEY_TSM==1
                             ,BACKLS &
#endif
                            )
             ENDIF
             IF(PRNLEV.GE.6) WRITE(IUNREX,'(A,I6,3F12.6,A,l1)') &
                'REPEXCH> PHREX: FNEIGH,ENEIGH,P,RN = ', FNEIGH, ENEIGH, P, RN, ' SUCCESS = ', QEXC
          ENDIF
          IF((IREPDSTR.EQ.NREPDSTR-1).AND.QRESHIGH) THEN
             NOPPUP=NOPPUP+1
             LEXATTEMPT=.TRUE.
             QDIDRSVR=.TRUE.
             IF(PRNLEV.GE.6) &
                WRITE(IUNREX,*) 'REPEXCH> TOP REPL WOULD EXCH W/ RESERVOIR'
             IF(QPHREX) THEN
                IDIDPHREX=1
                IPROTO=0
                DO I=1,NRES
                   IF(TSTATE(I).EQ.1) IPROTO=IPROTO+1
                ENDDO
                CALL RESEXCH_PH(TTEMP,.TRUE.,X,Y,Z,VX,VY,VZ,PHRX(IREPDSTR+1),IPROTO, & 
                                ISEED,IASVEL,IGVOPT,P,RN,JPROTO,PH_M,QEXC,FNEIGH &
#if KEY_TSM==1
                                ,BACKLS & 
#endif
                               )
             ELSE
                CALL RESEXCH(.TRUE.,X,Y,Z,VX,VY,VZ,QEXC,TEMPRX(IREPDSTR+1),MYEPOT, &
                             ISEED,IASVEL,IGVOPT,P,RN,ENEIGH,FNEIGH & 
#if KEY_TSM==1
                             ,BACKLS &
#endif
                            )
             ENDIF

             IF(PRNLEV.GE.6) WRITE(IUNREX,'(A,I6,3F12.6,A,l1)') &
                'TIM DBG> FNEIGH,ENEIGH,P,RN = ', FNEIGH, ENEIGH, P, RN, ' SUCCESS = ', QEXC
          ENDIF
       ENDIF
       GOTO 99
    ENDIF

    IF(QPHREX) THEN
      ! we are state i and have a pH value of pH_l, our neighbor
      ! is state j and has a pH value of pH_m

      ! count number of prototnated residues (in state 1) and swap with
      ! the neighbor
      ididphrex = 1
      iproto = 0
      do i=1,nres
         if(tstate(i).eq.1) iproto = iproto + 1 
      enddo
      if(irepdstr.gt.neighbor) then
         call grec(cneighbor,5,jproto,4)
         call gsen(cneighbor,6,iproto,4)
      else
         call gsen(cneighbor,5,iproto,4)
         call grec(cneighbor,6,jproto,4)
      endif

      ph_l = phrx(irepdstr+1)
      ph_m = phrx(neighbor+1)

      ph_delta = log(10.0)*(ph_m - ph_l)*(iproto - jproto)
      if(ph_delta.le.zero) then
         p=one
      else
         p=min(one,exp(-ph_delta))
      endif

      IF(QSGLD.AND.MYNOD.EQ.0)THEN
        ! we should be swapping SGLD stuffs too...
        SGARRAY(1)=EPOTLF
        SGARRAY(2)=EPOTHF+EPOTLF
        SGARRAY(3)=(AVGEFLF*AVGCFLF-AVGEFHF*AVGCFHF)/(kboltz*temprx(irepdstr+1))
        SGARRAY(4)=AVGEFHF*AVGCFHF/(kboltz*temprx(irepdstr+1))
        SGARRAY(5)=REPDID
        SGARRAY(6)=AVGTLF
        SGARRAY(7)=AVGEFLF
        SGARRAY(8)=AVGCFLF
        CALL GRECSEN(CNEIGHBOR,1,SGARRAYN,8,SGARRAY,8)
        ENEIGH=SGARRAYN(2)
      ENDIF
    ELSE
       IF(QRXSGLD) THEN
          p=min(one,exp(sgldarg))
       ELSE
          IF(QRXTHAM) THEN
             p=min(one,exp(hrexarg))
          ELSE
             CALL GRECSEN(CNEIGHBOR,1,ENEIGH,1,MYEPOT,1)
             if(temprx(irepdstr+1).eq.temprx(neighbor+1)) then
                p=ONE
             else
                P=MIN(ONE,EXP(-(ONE/(KBOLTZ*TEMPRX(IREPDSTR+1)) &
                     -ONE/(KBOLTZ*TEMPRX(NEIGHBOR+1)))*(ENEIGH-MYEPOT)))
             endif
          ENDIF
       ENDIF

    ENDIF

    !
    !     QEXC would be the result of the probability test...
    !
    QEXC=P.GT.RN
    !      if(prnlev.ge.2)write(IUNREX,'(a,i5,l5)')
    !     $     'REPEXCHG>me,qexc=',mynodg,qexc
    !
    IF(MOD(IREPDSTR,2).EQ.0)THEN
       CALL GREC(CNEIGHBOR,3,QEXC,4)  !GG: Updates QEXEC to all other replicas?
       CALL GREC(CNEIGHBOR,4,RN,8)
    ELSE
       CALL GSEN(CNEIGHBOR,3,QEXC,4)
       CALL GSEN(CNEIGHBOR,4,RN,8)
    ENDIF
    !      if(prnlev.ge.2)write(IUNREX,'(a,i5,l5)')
    !     $     'REPEXCHG>me,qexc=',mynodg,qexc
    !
    !     Perform the exchange of coordinates and velocities:
    IF(QEXC) THEN
       if(irepdstr.gt.neighbor) then
          nsucdn=nsucdn+1
          call grec(cneighbor,5,reptag,4)
          call gsen(cneighbor,6,oldrep,4)
       else
          nsucup=nsucup+1
          call gsen(cneighbor,5,oldrep,4)
          call grec(cneighbor,6,reptag,4)
       endif
       !
       call chmalloc('repdstr.src','REPXCHG','W',NATOM,crl=w)
       if(qphrex) then
          call chmalloc('repdstr.src','REPEXCHG','NSTATE',nres,intg=nstate)
          scaled=one
       else
          scaled=sqrt(temprx(irepdstr+1)/temprx(neighbor+1))
       endif

       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,X,NATOM)
       x(1:natom)=w(1:natom)
       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,Y,NATOM)
       y(1:natom)=w(1:natom)
       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,Z,NATOM)
       z(1:natom)=w(1:natom)
       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,VX,NATOM)
       vx(1:natom)=w(1:natom)
       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,VY,NATOM)
       vy(1:natom)=w(1:natom)
       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,VZ,NATOM)
       vz(1:natom)=w(1:natom)

       if(qrxtham) then
          CALL GRECSEN(CNEIGHBOR,11,W,NATOM,WMAIN,NATOM)
          wmain(1:natom) = w(1:natom)
          if(qcrys.and.xdim.gt.0) then
             call grecsen(cneighbor,4,W,6,XTLABC,6)
             XTLABC(1:6) = W(1:6)
          endif
       endif

#if KEY_PHMD==1
       IF(QPHRX)THEN !JAW - if exchange then swap theta and velocities
         CALL GRECSEN(CNEIGHBOR,2,L,NTITR,PH_THETA,NTITR)
         PH_THETA(1:ntitr) = l(1:ntitr)
         CALL GRECSEN(CNEIGHBOR,2,L,NTITR,VPH_THETA,NTITR)
         VPH_THETA(1:ntitr) = l(1:ntitr)
         CALL GRECSEN(CNEIGHBOR,2,L,NTITR,THETAOLD,NTITR)
         THETAOLD(1:ntitr) = l(1:ntitr)
         IF(PHBETA .gt. ZERO) THEN ! Langevin PHMD
            CALL GRECSEN(CNEIGHBOR,2,L,NTITR,VPHOLD,NTITR)
            VPHOLD(1:ntitr) = l(1:ntitr)
            CALL GRECSEN(CNEIGHBOR,2,L,NTITR,DPHOLD,NTITR)
            DPHOLD(1:ntitr) = l(1:ntitr)
         ENDIF
       ENDIF !JAW
#endif
#if KEY_BLOCK==1
       IF (QMSPHRX) THEN
          write(outu,'(a)') 'Transferring Theta variables'
          K = nsitemld*nblock                          !GG: Calc total no. of elements in MSLD array
          THETAMLDS = RESHAPE(THETAMLD,SHAPE(M))       !GG: Reshape array for coordinates to 1D array
          CALL GRECSEN(CNEIGHBOR,2,M,K,THETAMLDS,K)    !GG: Transfer array
          THETAMLDS(1:K) = M(1:K)
          THETAMLD = RESHAPE(THETAMLDS,SHAPE(N))       !GG: Reshape array back to nD array
          THETAMLDOLDS = RESHAPE(THETAMLDOLD,SHAPE(M)) !GG: Processing previous theta coordinates
          CALL GRECSEN(CNEIGHBOR,2,M,K,THETAMLDOLDS,K)
          THETAMLDOLDS(1:K) = M(1:K)
          THETAMLDOLD = RESHAPE(THETAMLDOLDS,SHAPE(N))
          THETAVMLDS = RESHAPE(THETAVMLD,SHAPE(M))     !GG: Processing current velocity coordinates
          CALL GRECSEN(CNEIGHBOR,2,M,K,THETAVMLDS,K)
          THETAVMLDS(1:K) = M(1:K)
          THETAVMLD = RESHAPE(THETAVMLDS,SHAPE(N))
          THETAFMLDS = RESHAPE(THETAFMLD,SHAPE(M))     !GG: Processing current force coordinates
          CALL GRECSEN(CNEIGHBOR,2,M,K,THETAFMLDS,K)
          THETAFMLDS(1:K) = M(1:K)
          THETAFMLD = RESHAPE(THETAFMLDS,SHAPE(N))
       ENDIF
#endif

       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,XOLD,NATOM)
       xOLD(1:natom)=w(1:natom)
       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,YOLD,NATOM)
       yOLD(1:natom)=w(1:natom)
       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,ZOLD,NATOM)
       zOLD(1:natom)=w(1:natom)

       IF(QPHREX) THEN
          call grecsen(cneighbor,2,w,natom,cg,natom)
          cg(1:natom)=w(1:natom)
          if(irepdstr.gt.neighbor) then
             call grec(cneighbor,5,nstate,4*nres)
             call gsen(cneighbor,6,tstate,4*nres)
          else
             call gsen(cneighbor,5,tstate,4*nres)
             call grec(cneighbor,6,nstate,4*nres)
          endif
          tstate(1:nres)=nstate(1:nres)

          IF(QSGLD)THEN
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGVX,NATOM)
             sgvx(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGVY,NATOM)
             sgvy(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGVZ,NATOM)
             sgvz(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGFX,NATOM)
             sgfx(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGFY,NATOM)
             sgfy(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGFZ,NATOM)
             sgfz(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGGX,NATOM)
             sggx(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGGY,NATOM)
             sggy(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGGZ,NATOM)
             sggz(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGHX,NATOM)
             sghx(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGHY,NATOM)
             sghy(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGHZ,NATOM)
             sghz(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGKX,NATOM)
             sgkx(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGKY,NATOM)
             sgky(1:natom)=w(1:natom)
             CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGKZ,NATOM)
             sgkz(1:natom)=w(1:natom)
             scalsg=sqrt(SGARRAY(6)/SGARRAYN(6))
             REPDID=NINT(SGARRAYN(5))
             EPOTLF=SGARRAYN(1)
             DTEFLF=SGARRAY(7)-ONE-(SGARRAYN(7)-ONE)*SCALSG
             DTCFLF=SGARRAY(8)-ONE-(SGARRAYN(8)-ONE)*SCALSG
             DO I = 1, NATOM
                SGVX(I)=SGVX(I)*SCALSG
                SGVY(I)=SGVY(I)*SCALSG
                SGVZ(I)=SGVZ(I)*SCALSG
                SGGX(I)=SGGX(I)*SCALSG+DTEFLF*SGFX(I)
                SGGY(I)=SGGY(I)*SCALSG+DTEFLF*SGFY(I)
                SGGZ(I)=SGGZ(I)*SCALSG+DTEFLF*SGFZ(I)
                FACT=DTCFLF*TIMFAC*FBETA(I)*AMASS(I)/DELTA
                SGHX(I)=SGHX(I)*SCALSG+FACT*SGVX(I)
                SGHY(I)=SGHY(I)*SCALSG+FACT*SGVY(I)
                SGHZ(I)=SGHZ(I)*SCALSG+FACT*SGVZ(I)
             ENDDO
          ENDIF
       ENDIF

       IF(QRXSGLD)THEN
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGVX,NATOM)
          sgvx(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGVY,NATOM)
          sgvy(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGVZ,NATOM)
          sgvz(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGFX,NATOM)
          sgfx(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGFY,NATOM)
          sgfy(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGFZ,NATOM)
          sgfz(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGGX,NATOM)
          sggx(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGGY,NATOM)
          sggy(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGGZ,NATOM)
          sggz(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGHX,NATOM)
          sghx(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGHY,NATOM)
          sghy(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGHZ,NATOM)
          sghz(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGKX,NATOM)
          sgkx(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGKY,NATOM)
          sgky(1:natom)=w(1:natom)
          CALL GRECSEN(CNEIGHBOR,2,W,NATOM,SGKZ,NATOM)
          sgkz(1:natom)=w(1:natom)
       ENDIF
       !
       call chmdealloc('repdstr.src','REPXCHG','W',NATOM,crl=w)
       if(qphrex) call chmdealloc('repdstr.src','REPEXCHG','NSTATE',nres,intg=nstate)
       !
       !         if(prnlev.ge.2)write(IUNREX,'(a,i5,2f20.8)')
       !     $ 'REPEXCHG>me,temps=',mynodg,temprx(irepdstr+1),temprx(neighbor+1)

       DO I = 1, NATOM
          VX(I)=VX(I)*SCALED
          VY(I)=VY(I)*SCALED
          VZ(I)=VZ(I)*SCALED
          XOLD(I)=XOLD(I)*SCALED
          YOLD(I)=YOLD(I)*SCALED
          ZOLD(I)=ZOLD(I)*SCALED
       ENDDO

       IF(QRXSGLD)THEN
          scalsg=sqrt(SGARRAY(6)/SGARRAYN(6))
          REPDID=NINT(SGARRAYN(5))
          EPOTLF=SGARRAYN(1)
          DTEFLF=SGARRAY(7)-ONE-(SGARRAYN(7)-ONE)*SCALSG
          DTCFLF=SGARRAY(8)-ONE-(SGARRAYN(8)-ONE)*SCALSG
          DO I = 1, NATOM
             SGVX(I)=SGVX(I)*SCALSG
             SGVY(I)=SGVY(I)*SCALSG
             SGVZ(I)=SGVZ(I)*SCALSG
             SGGX(I)=SGGX(I)*SCALSG+DTEFLF*SGFX(I)
             SGGY(I)=SGGY(I)*SCALSG+DTEFLF*SGFY(I)
             SGGZ(I)=SGGZ(I)*SCALSG+DTEFLF*SGFZ(I)
             FACT=DTCFLF*TIMFAC*FBETA(I)*AMASS(I)/DELTA
             SGHX(I)=SGHX(I)*SCALSG+FACT*SGVX(I)
             SGHY(I)=SGHY(I)*SCALSG+FACT*SGVY(I)
             SGHZ(I)=SGHZ(I)*SCALSG+FACT*SGVZ(I)
          ENDDO
       ENDIF
#if KEY_PHMD==1
       IF(QPHRX)THEN
          DO I=1,NTITR           
             VPH_THETA(I)=VPH_THETA(I)*SCALED           
          ENDDO
          IF(PHBETA .gt. ZERO)THEN
             VPHOLD(I)=VPHOLD(I)*SCALED
          ENDIF
       ENDIF
#endif
#if KEY_BLOCK==1
       IF (QMSPHRX) THEN
          K = nsitemld*nblock
          DO I=1,K
             !GG: Scale the 1D array directly, since that is the array that will be broadcasted
             THETAVMLDS(I)=THETAVMLDS(I)*SCALED
             THETAVMLD = RESHAPE(THETAVMLDS,SHAPE(N))
          ENDDO
       ENDIF
#endif
       !ENDIF
       !         write(50+mynodg,'(3f20.10)')(x(i),y(i),z(i),i=1,natom)
       !         close(50+mynodg)
       !         call psync
       !         call mpi_finalize(i)
       !         stop
 

    ENDIF !GG: End loop for QEXC = .TRUE.
    !
    !     We need to broadcast the data within the same replica group
    !
99  CONTINUE

    CALL PSND8(ENEIGH,1)
    CALL PSND8(P,1)
    CALL PSND8(RN,1)
    CALL PSND4(QEXC,1)
    IF(QPHREX) THEN
       CALL PSND4(ididphrex,1)
       CALL PSND4(iproto,1)
       CALL PSND4(jproto,1)
    ENDIF

    !     but we need these:
    !

    IF(QEXC)THEN
#if KEY_BLOCK==1
       call msld_checkvariables(2) !GG: Before transmission
#endif

       ! reset the dynamics algorithm if there's an exchange
       ! Do we need to set igvopt to 2 since this is what ASSVEL does?
       ! assvel does
       JHSTRT=0
       IGVOPT=2
       CALL PSND4(IGVOPT,1)
       CALL PSND4(JHSTRT,1)
       CALL PSND8(X,NATOM)
       CALL PSND8(Y,NATOM)
       CALL PSND8(Z,NATOM)
       CALL PSND8(VX,NATOM)
       CALL PSND8(VY,NATOM)
       CALL PSND8(VZ,NATOM)
       CALL PSND8(XOLD,NATOM)
       CALL PSND8(YOLD,NATOM)
       CALL PSND8(ZOLD,NATOM)

       IF(QRXTHAM) THEN
          CALL PSND8(WMAIN,NATOM)
          if(qcrys) then
             xtlabc(1:6)=oldxtlabc(1:6)
             call xtllat(xucell,xtlabc)
             call chmalloc('repdstr.src','REPXCHG','TRANSF',3,4,XNSYMM,crl=transf)
             call imfill(transf,.false.)
             call chmdealloc('repdstr.src','REPXCHGL','TRANSF',3,4,XNSYMM,crl=transf)
             call upimag(x,y,z,wmain,0,x,y,z,vx,vy,vz)
          endif
          ! not sure if we need these calls, but I am leaving them in now
          ! to be safe
          call nbonds(x,y,z,bnbnd,bimag)
          call energy(x,y,z,dx,dy,dz,bnbnd,bimag,0)
       ENDIF

       IF(QRXSGLD) THEN
          ! Xiongwu doesn't do this, but it probably needs
          ! to be done -- ask about it!
          CALL PSND8(SGVX,NATOM)
          CALL PSND8(SGVY,NATOM)
          CALL PSND8(SGVZ,NATOM)
          CALL PSND8(SGFX,NATOM)
          CALL PSND8(SGFY,NATOM)
          CALL PSND8(SGFZ,NATOM)
          CALL PSND8(SGGX,NATOM)
          CALL PSND8(SGGY,NATOM)
          CALL PSND8(SGGZ,NATOM)
          CALL PSND8(SGKX,NATOM)
          CALL PSND8(SGKY,NATOM)
          CALL PSND8(SGKZ,NATOM)
          CALL PSND8(SGHX,NATOM)
          CALL PSND8(SGHY,NATOM)
          CALL PSND8(SGHZ,NATOM)
       ENDIF

       IF(QPHREX) THEN
          CALL PSND8(CG,NATOM)
          IF(QRXSGLD) THEN
             ! Xiongwu doesn't do this, but it probably needs
             ! to be done -- ask about it!
             CALL PSND8(SGVX,NATOM)
             CALL PSND8(SGVY,NATOM)
             CALL PSND8(SGVZ,NATOM)
             CALL PSND8(SGFX,NATOM)
             CALL PSND8(SGFY,NATOM)
             CALL PSND8(SGFZ,NATOM)
             CALL PSND8(SGGX,NATOM)
             CALL PSND8(SGGY,NATOM)
             CALL PSND8(SGGZ,NATOM)
             CALL PSND8(SGKX,NATOM)
             CALL PSND8(SGKY,NATOM)
             CALL PSND8(SGKZ,NATOM)
             CALL PSND8(SGHX,NATOM)
             CALL PSND8(SGHY,NATOM)
             CALL PSND8(SGHZ,NATOM)
          ENDIF
       ENDIF

#if KEY_PHMD==1
       IF(QPHRX)THEN
          CALL PSND8(PH_THETA,NTITR)
          CALL PSND8(VPH_THETA,NTITR) 
          IF(PHBETA .gt. ZERO)THEN
             CALL PSND8(VPHOLD,NTITR)
             CALL PSND8(DPHOLD,NTITR)
          ENDIF
       ENDIF
#endif

#if KEY_BLOCK==1
       IF (QMSPHRX) THEN   !GG: Broadcast swapped data to other replicas
          K = nsitemld*nblock
          CALL PSND8(THETAMLDS,K)
          CALL PSND8(THETAMLDOLDS,K)
          CALL PSND8(THETAVMLDS,K)
          CALL PSND8(THETAFMLDS,K)
          IF (MYNOD.EQ.0) GOTO 991
          THETAMLD = RESHAPE(THETAMLDS,SHAPE(N))
          THETAMLDOLD = RESHAPE(THETAMLDOLDS,SHAPE(N))
          THETAVMLD = RESHAPE(THETAVMLDS,SHAPE(N))
          THETAFMLD = RESHAPE(THETAFMLDS,SHAPE(N))
991       CONTINUE
       ENDIF
#endif
    ENDIF  !GG: End loop for QEXC = .TRUE.

    !
    !     Reinitialize the stuff for TIGER
    !
    IF(QRXTIGER)THEN
       qpxtiger=.false.
       qrxtmin=.false.
       tigeriti=0
    ENDIF
    !
    ! Printout the results:
    !
    ! on each exchange during dynamics
    irex=irex+1
 
    if(noppup > 0) then
       exratup=float(nsucup)/float(noppup)
       call set_param('EXRUP',exratup)
    endif
    if(noppdn > 0) then
       exratdn=float(nsucdn)/float(noppdn)
       call set_param('EXRDN',exratdn)
    endif
    if((noppup+noppdn.gt.0).and.lexattempt) then
       reprob=(min(p,one)+((noppup+noppdn-1)*reprob))/float(noppup+noppdn)
       call set_param('REPROB',reprob)
       if(prnlev > 5) then
          write(iunrex,'(a,4i3)') 'REPEXCH> noppup, noppdn, nsucup, nsucdn = ',noppup,noppdn,nsucup,nsucdn
          write(iunrex,'(a,f10.5)') 'REPEXCH> P = ',min(p,one)
          write(iunrex,'(a,f10.5)') 'REPEXCH> NUMERATOR = ',(min(p,one)+((noppup+noppdn-1)*reprob))
          write(iunrex,'(a,f10.5)') 'REPEXCH> DENOMINATOR = ',float(noppup+noppdn)
          write(iunrex,'(a,3f6.3)') 'REPEXCH> EXRUP EXRDN REPROB = ',exratup,exratdn,reprob
       endif
    endif

    if(prnlev.ge.2) then
       !
       if(qexc)isuc=isuc+1
       srate=real(isuc)/real(irex)

       if(qphrex) then
          write(iunrex,'(a)') &
               '------------- pH Replica Exchange ------------'
          write(iunrex,'(a,i10,a,i10)') &
               'REX>EXCHANGE = ', irex, '  Step =', istart-1
          write(iunrex,'(a,i5,a,f7.3,a,i5)') &
               'REX>REPL     = ',irepdstr, &
               '  pH = ', phrx(irepdstr+1), &
               ' nproto = ', iproto
          write(iunrex,'(a,i5,a,f7.3,a,i5)') &           
               'REX>NEIGHBOR = ',neighbor, &
               '  pH = ', ph_m,     &
               ' nproto = ', jproto
       else
          ttx=zero
          if(neighbor.ge.0)ttx=temprx(neighbor+1)
          write(iunrex,'(a)') &
               '------------- Replica Exchange ------------'
          write(iunrex,'(a,i10,a,i10)') &
               'REX>EXCHANGE = ', irex, '  Step =', istart-1
          write(iunrex,'(a,i5,a,f7.3,a,f20.8)') &
               'REX>REPL     = ',irepdstr,      &
               '  Temp = ', temprx(irepdstr+1), &
               '  Epot = ', myepot
          write(iunrex,'(a,i5,a,f7.3,a,f20.8)') &
               'REX>NEIGHBOR = ',neighbor,      &
               '  Temp = ', ttx,                &
               '  Epot = ', eneigh
          if(qrxsgld) then
             !write(iunrex,'(a6,i6,1x,f12.2,f12.2,f6.3,f6.3,f6.3,1x,f7.3)') &
             !    'RXSG> ',irex,myepot,eneigh,scaled,scalsg,srate,sgldarg

             write(iunrex,'(a6,i6,1x,i10,i4,i4,i4,f12.2,f12.2,f6.3,f6.3,f6.3,1x,l1)') &
                  'RXSG> ',irex,istart-1,irepdstr,neighbor,repdid,  &
                   myepot,eneigh,scaled,scalsg,srate,qexc

          endif
          if(qrxtham) then
             write(iunrex,'(a6,4f12.2,x,f11.3)') &
                 'THAM> ',nbrpotj,nbrpoti,ourpoti,ourpotj,hrexarg
          endif
       endif
       write(iunrex,'(a,i5,a,i5)') 'REX>ORIGINAL TAG ',oldrep,' NEW TAG ',reptag

       if(qdidrsvr) &
          write(iunrex,'(a,i7)') 'REX>RESERVOIR STRUCT = ', fneigh
       write(iunrex,'(a,f8.5,a,f8.5,a,f7.4,a,l1)') &
            'REX>PROB     = ', p, ' Rand = ', rn,  &
            ' Tscale = ', scaled,                  &
            ' Success = ', qexc

       !
       ! Printout the summary results if sump flag specified:
       !
       if(qsump) then
          allocate(comar(6))
          comar(1)=scaled
          comar(2)=srate
          comar(3)=temprx(irepdstr+1)
          comar(4)=temprx(neighbor+1)
          comar(5)=ttemp
          comar(6)=myepot
          write(iunrex,'(a,a,a)')'#REXSUM Rep# Tscale', &
               '  Sratio   Temp      NewTemp',          &
               '  CurrTemp       Epot'
          write(iunrex,'(a,i3,2f8.3,3f10.3,f15.4)')'REXSUM> ', &
               irepdstr,(comar(j),j=1,6)
          !
          ! Lets do some extra communication for a nice rexsum> printout
          !
          if(irepdstr.eq.0)then
             do i = 1, nrepdstr-1
                call grec(i*numnod,5,comar,6*8)
                write(iunrex,'(a,i3,2f8.3,3f10.3,f15.4)')'REXSUM> ', &
                     i,(comar(j),j=1,6)
             enddo
          else
             call gsen(0,5,comar,6*8)
          endif
          deallocate(comar)
       endif

       if(.not.qrxsgld.and..not.qphrex) then
          write(iunrex,'(a)') &
            '------------- Replica Exchange End --------'
       else if(qphrex) then
          write(iunrex,'(a)') &
            '------------- pH Replica Exchange End --------'
       endif

    endif
    !
    !     Since we exchanged the coordinates we should perform
    !     image and non-bond update. NOTE! Order is critical...
    !
    CALL UPIMAG(X,Y,Z,WMAIN,0,X,Y,Z,VX,VY,VZ)
    call nbonds(x,y,z,bnbnd,bimag)

#if KEY_PHMD==1
    IF(QPHRX)THEN !JAW
      call UpdatePHMD(1, 1, 0, 1)
    ENDIF
#endif

#if KEY_BLOCK==1
    IF (QMSPHRX) THEN
       call msld_swapcoeff_mod     !GG: Toggle msld_setblcoef_fnexp to use thetamold values
       call msld_checkvariables(3) !GG: Variables after Swap, before Call Energy
    ENDIF
#endif
    CALL ENERGY(X,Y,Z,DX,DY,DZ,BNBND,BIMAG,0) !JAW
#if KEY_BLOCK==1
    IF (QMSPHRX) THEN
       call msld_swapcoeff_norm    !GG: Toggle msld_setblcoef_fnexp to use thetamold values
       call msld_checkvariables(4) !GG: Variables after Swap, after Call Energy
    ENDIF
#endif

    RETURN
  END SUBROUTINE REPEXCHG
  !

  SUBROUTINE GARRETTSPACER
    RETURN   !GG:  Subroutine Spacer to get bookmarks working correctly
  END SUBROUTINE GARRETTSPACER


  SUBROUTINE REPEXCHGL(X,Y,Z,WMAIN,VX,VY,VZ,EPTT,TEMNEW, &
       ISTART,ISEED,IASVEL,IGVOPT,JHSTRT &
#if KEY_TSM==1
       ,BACKLS                    &            
#endif
       )

    !-----------------------------------------------------------------------
    !     Perform necessary communication and calculate new data...
    !
  use chm_kinds
  use chm_types
  use number
  use consta
  use dimens_fcm
  use energym
  use deriv
  use image
  use psf
  use stream
  use parallel
  use repdstr
  use memory
  use phmd !JW
  use bases_fcm
  use imgup
#if KEY_GCMC==1
  use gcmc          
#endif
  use clcg_mod,only: random
  use block_ltm
  use lambdam    !GG: MSLD-compatibility
  use sgld,only: SGVX,SGVY,SGVZ,SGFX,SGFY,SGFZ,SGGX,SGGY,SGGZ,SGHX,SGHY,SGHZ,SGKX,SGKY,SGKZ  !GG: SGLD-compatibility

    !
    implicit none

    REAL(chm_real) :: X(:), Y(:), Z(:), WMAIN(*),EPTT,VX(*),VY(*),VZ(*)
    REAL(chm_real) LAMBDA,LAMBDA2,ECSUM, EMTS(LENENT)
    INTEGER ISTART,OLDREP,J
    REAL(chm_real) TEMNEW,RN
    INTEGER ISEED,IASVEL,IGVOPT,JHSTRT

#if KEY_TSM==1
    INTEGER BACKLS(*)
#endif
    !
    LOGICAL QEXC,QCRYS,LEXATTEMPT
    INTEGER STEP,IDECIDE,ME,MEG,NEIGHBOR,I,CNEIGHBOR,IREPR,IREPL
    REAL(chm_real) ENEIGH,ENEIGHI,ENEIGHJ,EPOTI,EPOTJ,P
    real(chm_real),allocatable,dimension(:) :: W
#if KEY_PHMD==1
    real(chm_real) :: L(NTitr) !JW
#endif
#if KEY_BLOCK==1
    integer K                                      !GG: MSLD-compatibility
    real(chm_real) n(nsitemld,nblock)              !GG: MSLD-compatibility
    real(chm_real) m(nsitemld*nblock)              !GG: MSLD-compatibility
    real(chm_real) THETAVMLDS(nsitemld*nblock)     !GG: MSLD-compatibility
    real(chm_real) THETAMLDS(nsitemld*nblock)      !GG: MSLD-compatibility
    real(chm_real) THETAMLDOLDS(nsitemld*nblock)   !GG: MSLD-compatibility
    real(chm_real) THETAFMLDS(nsitemld*nblock)     !GG: MSLD-compatibility
    real(chm_real) EDIFF                           !GG: MSLD-compatibility
#endif
    logical,allocatable,dimension(:) :: WL
    real(chm_real),allocatable,dimension(:,:,:) :: TRANSF
    !

    IF(MOD(ISTART-1,IREXFQ).NE.0) RETURN
    OLDREP=REPTAG
    LEXATTEMPT=.FALSE.

    !
    !     Prepare the variable NEIGHBOR.
    !     NEIGHBOR=-1 if no communication is needed on this process
    IF((.NOT.QEXPT).AND.(.NOT.QEX2D))THEN !GG: CPHMD^MSLD uses this protocol (EXLM keyword only)
       STEP=MOD(ISTART/IREXFQ,2)          !GG: What is the current MC step
       ME=MOD(IREPDSTR,2)                 !GG: Is the current replica no. even (STEP=0) or odd (STEP=1)?
       IF(STEP.EQ.1)THEN
          NEIGHBOR=IREPDSTR+1
          IF(ME.NE.0)NEIGHBOR=IREPDSTR-1
       ELSE
          NEIGHBOR=IREPDSTR-1
          IF(ME.NE.0)NEIGHBOR=IREPDSTR+1
       ENDIF
    ELSE IF (QEXPT)THEN
       STEP=MOD(ISTART/IREXFQ,4)
       ME=MOD(IREPDSTR,2)
       MEG=MOD(IREPDSTR/NREPT,2)
       IF(STEP.EQ.1)THEN
          NEIGHBOR=IREPDSTR+1
          IF(ME.NE.0)NEIGHBOR=IREPDSTR-1
          IF((ME.EQ.0).AND.(MOD(NEIGHBOR,NREPT).EQ.0))NEIGHBOR=-1
          IF((ME.NE.0).AND.(MOD((NEIGHBOR+1),NREPT).EQ.0))NEIGHBOR=-1
       ELSE IF(STEP.EQ.2)THEN
          NEIGHBOR=IREPDSTR-1
          IF(ME.NE.0)NEIGHBOR=IREPDSTR+1
          IF((ME.NE.0).AND.(MOD(NEIGHBOR,NREPT).EQ.0))NEIGHBOR=-1
          IF((ME.EQ.0).AND.(MOD((NEIGHBOR+1),NREPT).EQ.0))NEIGHBOR=-1
       ELSE IF(STEP.EQ.3)THEN
          IF(MOD(IREPDSTR,NREPT).EQ.0)THEN
             NEIGHBOR=IREPDSTR+NREPT
             IF(MEG.NE.0)NEIGHBOR=IREPDSTR-NREPT
          ELSE
             NEIGHBOR=-1
          ENDIF
       ELSE
          IF(MOD(IREPDSTR,NREPT).EQ.0)THEN
             NEIGHBOR=IREPDSTR-NREPT
             IF(MEG.NE.0)NEIGHBOR=IREPDSTR+NREPT
          ELSE
             NEIGHBOR=-1
          ENDIF
       ENDIF
    ELSE                              !GG: EX2D uses this???
      STEP=MOD(ISTART/IREXFQ,4)       !GG: What is the current MC step
      ME=MOD(IREPDSTR,2)              !GG: Is the current replica even (STEP=0) or odd (STEP=1)?
      MEG=MOD(IREPDSTR/NREPX,2)
      IF(STEP.EQ.1)THEN
         NEIGHBOR=IREPDSTR+1
         IF(ME.NE.0)NEIGHBOR=IREPDSTR-1
         IF((ME.EQ.0).AND.(MOD(NEIGHBOR,NREPX).EQ.0))NEIGHBOR=-1
         IF((ME.NE.0).AND.(MOD((NEIGHBOR+1),NREPX).EQ.0))NEIGHBOR=-1
      ELSE IF(STEP.EQ.2)THEN
         NEIGHBOR=IREPDSTR-1
         IF(ME.NE.0)NEIGHBOR=IREPDSTR+1
         IF((ME.NE.0).AND.(MOD(NEIGHBOR,NREPX).EQ.0))NEIGHBOR=-1
         IF((ME.EQ.0).AND.(MOD((NEIGHBOR+1),NREPX).EQ.0))NEIGHBOR=-1
      ELSE IF(STEP.EQ.3)THEN
         NEIGHBOR=IREPDSTR+NREPX
         IF(MEG.NE.0)NEIGHBOR=IREPDSTR-NREPX
      ELSE
         NEIGHBOR=IREPDSTR-NREPX
         IF(MEG.NE.0)NEIGHBOR=IREPDSTR+NREPX
      ENDIF
    ENDIF
    IF(NEIGHBOR.LT.0)NEIGHBOR=-1
    IF(NEIGHBOR.GE.NREPDSTR)NEIGHBOR=-1
    IF(QEXBK)THEN
       IREPR=IREPDSTR+1
       IREPL=IREPDSTR-1
       IF((IREPR.EQ.IRBK).AND.(NEIGHBOR.EQ.IRBK))NEIGHBOR=-1
       IF((IREPDSTR.EQ.IRBK).AND.(NEIGHBOR.EQ.IREPL))NEIGHBOR=-1
    ENDIF

    IF(NEIGHBOR.EQ.-1) THEN
       IF(QRESERVOIR) THEN
          IF((IREPDSTR.EQ.0).AND.QRESLOW) THEN
             NOPPDN=NOPPDN+1
             LEXATTEMPT=.TRUE.
             CALL RESEXCHL(.FALSE.,X,Y,Z,VX,VY,VZ,TEMNEW,RHTEMP,EPTT, &
                           ISEED,IASVEL,IGVOPT,P,ISTART-1,JHSTRT &
#if KEY_TSM==1
                           ,BACKLS & 
#endif
                          )
          ENDIF
          IF((IREPDSTR.EQ.NREPDSTR-1).AND.QRESHIGH) THEN
             NOPPUP=NOPPUP+1
             LEXATTEMPT=.TRUE.
             CALL RESEXCHL(.TRUE.,X,Y,Z,VX,VY,VZ,TEMNEW,RLTEMP,EPTT, &
                           ISEED,IASVEL,IGVOPT,P,ISTART-1,JHSTRT &
#if KEY_TSM==1
                           ,BACKLS &
#endif
                          )
          ENDIF
       ENDIF
       GOTO 222
    ENDIF

    IF(IREXFQ .GE. NWREX)THEN
#if KEY_PHMD==1
       IF(.NOT. QPHRX)THEN
#endif
       !!write(IUNREX,*)'REPEXCHGL>me,irep,istart,step,neighbor=', &
       !!     mynodg,irepdstr,istart-1,step,neighbor
#if KEY_PHMD==1
       ENDIF
#endif
    ENDIF

    QCRYS = (XTLTYP.NE.'    ')
    !
    !     FIXME: Maybe for protection we can call psetloc() here
    CNEIGHBOR=NEIGHBOR*NUMNOD
    IF(prnlev.ge.7) write(outu,'(a5, i2)') 'MYNOD', MYNOD !GG: Prints out current node

    ! -----------------------------GG: Start of MYNOD .EQ. 0 processes -----------------------------------------

    IF (MYNOD.NE.0) GOTO 99

#if KEY_BLOCK==1
    call msld_checkvariables(1) !GG Variables Before MC Exchange
#endif

    LEXATTEMPT=.TRUE.
    CALL GRECSEN(CNEIGHBOR,1,ENEIGH,1,EPTT,1)
    ENEIGHI = ENEIGH
    EPOTI = EPTT

    !     Perform a test coordinate exchange and calculate the new energies
    call chmalloc('repdstr.src','REPXCHGL','W',NATOM,crl=w)
    CALL GRECSEN(CNEIGHBOR,2,W,NATOM,X,NATOM)
    x(1:natom) = w(1:natom)
    CALL GRECSEN(CNEIGHBOR,2,W,NATOM,Y,NATOM)
    y(1:natom) = w(1:natom)
    CALL GRECSEN(CNEIGHBOR,2,W,NATOM,Z,NATOM)
    z(1:natom) = w(1:natom)
#if KEY_PHMD==1
    IF (QPHRX) THEN ! JAW. EXCHANGE THETA VALUES BEFORE TESTING EXCHANGE  
       CALL GRECSEN(CNEIGHBOR,2,L,NTITR,PH_THETA,NTITR)
       PH_THETA(1:ntitr) = l(1:ntitr)
    ENDIF
#endif
#if KEY_BLOCK==1
    IF (QMSPHRX) THEN
       write(outu,'(a)') 'Transferring Theta variables'
       K = nsitemld*nblock                          !GG: Calc total no. of elements in MSLD array
       THETAMLDS = RESHAPE(THETAMLD,SHAPE(M))       !GG: Reshape array for coordinates to 1D array
       CALL GRECSEN(CNEIGHBOR,2,M,K,THETAMLDS,K)    !GG: Transfer array
       THETAMLDS(1:K) = M(1:K)
       THETAMLD = RESHAPE(THETAMLDS,SHAPE(N))       !GG: Reshape array back to nD array
       THETAMLDOLDS = RESHAPE(THETAMLDOLD,SHAPE(M)) !GG: Processing previous theta coordinates
       CALL GRECSEN(CNEIGHBOR,2,M,K,THETAMLDOLDS,K)
       THETAMLDOLDS(1:K) = M(1:K)
       THETAMLDOLD = RESHAPE(THETAMLDOLDS,SHAPE(N))
       THETAVMLDS = RESHAPE(THETAVMLD,SHAPE(M))     !GG: Processing current velocity coordinates
       CALL GRECSEN(CNEIGHBOR,2,M,K,THETAVMLDS,K)
       THETAVMLDS(1:K) = M(1:K)
       THETAVMLD = RESHAPE(THETAVMLDS,SHAPE(N))
       THETAFMLDS = RESHAPE(THETAFMLD,SHAPE(M))     !GG: Processing current force coordinates
       CALL GRECSEN(CNEIGHBOR,2,M,K,THETAFMLDS,K)
       THETAFMLDS(1:K) = M(1:K)
       THETAFMLD = RESHAPE(THETAFMLDS,SHAPE(N))
    ENDIF
#endif

    call chmdealloc('repdstr.src','REPXCHGL','W',NATOM,crl=w)
    !Exchange cell dimension information
      IF (QCRYS .AND. XDIM.GT.0) THEN
         call chmalloc('repdstr.src','REPXCHGL','W',6,crl=w)
         CALL GRECSEN(CNEIGHBOR,4,W,6,XTLABC,6)
         XTLABC(1:6) = W(1:6)
         call chmdealloc('repdstr.src','REPXCHGL','W',6,crl=w)
      ENDIF
      call chmalloc('repdstr.src','REPXCHGL','W',NATOM,crl=w)
      CALL GRECSEN(CNEIGHBOR,11,W,NATOM,WMAIN,NATOM)
      WMAIN(1:natom) = w(1:natom)
      call chmdealloc('repdstr.src','REPXCHGL','W',NATOM,crl=w)
#if KEY_GCMC==1
    IF (QGCMC) THEN
       call chmalloc('repdstr.src','REPXCHGL','WL',MAXA,log=wl)
       CALL GRECSEN(CNEIGHBOR,9,WL,MAXA,GCMCON,MAXA)
       gcmcon(1:maxa) = WL(1:MAXA)
       call chmdealloc('repdstr.src','REPXCHGL','WL',MAXA,log=wl)
    ENDIF
#endif
    ! -----------------------------GG: Start of MYNOD .EQ. ALL processes -----------------------------------------
99  CONTINUE
#if KEY_BLOCK==1
    call msld_checkvariables(2) !GG: Variables Before Transmission
#endif
    CALL PSND8(X,NATOM)
    CALL PSND8(Y,NATOM)
    CALL PSND8(Z,NATOM)
    CALL PSND8(WMAIN, NATOM)
      IF (QCRYS) THEN
         CALL PSND8(XTLABC,6)
         CALL XTLLAT(XUCELL,XTLABC)
         call chmalloc('repdstr.src','REPXCHGL','TRANSF',3,4,XNSYMM,crl=transf) 
         CALL IMFILL(TRANSF,.FALSE.)
         call chmdealloc('repdstr.src','REPXCHGL','TRANSF',3,4,XNSYMM,crl=transf) 
         CALL UPIMAG(X,Y,Z,WMAIN,0,X,Y,Z,VX,VY,VZ)
      ENDIF
#if KEY_BLOCK==1
      IF (QMSPHRX) THEN   !GG: Broadcast swapped data to other replicas
         K = nsitemld*nblock
         CALL PSND8(THETAMLDS,K)
         CALL PSND8(THETAMLDOLDS,K)
         CALL PSND8(THETAVMLDS,K)
         CALL PSND8(THETAFMLDS,K)
         IF (MYNOD.EQ.0) GOTO 991
         THETAMLD = RESHAPE(THETAMLDS,SHAPE(N))
         THETAMLDOLD = RESHAPE(THETAMLDOLDS,SHAPE(N))
         THETAVMLD = RESHAPE(THETAVMLDS,SHAPE(N))
         THETAFMLD = RESHAPE(THETAFMLDS,SHAPE(N))
991      CONTINUE
      ENDIF
#endif
#if KEY_GCMC==1
    IF (QGCMC) THEN
       CALL PSND4(GCMCON,MAXA)
    ENDIF
#endif
#if KEY_PHMD==1
    IF(QPHRX)THEN
       CALL UpdatePHMD(1, 1, 0, 1) ! JAW. After Exchanging Coordinates and Lambda it is necesarry to update charges and
                                   ! total system charge so energy calculations are correct. 
    ENDIF
#endif
#if KEY_BLOCK==1
    !GG No need for BLOCK version of charge update since charges are stored individually
#endif
    CALL NBONDS(X,Y,Z,BNBND,BIMAG)  ! LNBND,BIMAG,LIMAG)
#if KEY_BLOCK==1
    IF (QMSPHRX) THEN
       call msld_swapcoeff_mod      !GG: Toggle msld_setblcoef_fnexp to use thetamold values
       call msld_checkvariables(3) !GG: Variables After Swap, Before Call Energy
    ENDIF
#endif
    CALL ENERGY(X,Y,Z,DX,DY,DZ,BNBND,BIMAG,0)
#if KEY_BLOCK==1
    IF (QMSPHRX) THEN
       call msld_swapcoeff_norm       !GG: Toggle msld_setblcoef_fnexp to use thetam values
       call msld_checkvariables(4) !GG: Variables  After Swap, After Call Energy
    ENDIF
#endif
    IF (QMSPHRX) THEN
       !GG: Not summing ETERM(I), since PE from MSLD biases have no individual ETERM(I) value
       IF(prnlev.ge.8) THEN
          write(outu,'(a)') 'Using MSLD-PHREX potential energy calculations'
       ENDIF
    ELSE
       !GG: Jana's implementation of traditional pH-REX CPHMD
       ECSUM=ZERO
       DO I=1,LENENT
          ECSUM=ETERM(I)+ECSUM
       ENDDO
       EPROP(EPOT)=ECSUM
    ENDIF

    IF (MYNOD.NE.0) GOTO 999

       CALL GRECSEN(CNEIGHBOR,1,ENEIGH,1,EPROP(EPOT),1)
       ENEIGHJ = ENEIGH
       EPOTJ = EPROP(EPOT)

    !
    IF (QMSPHRX) THEN
       !GG: Somehow TEMNEW is not parsed in correctly in MSLD-CPHMD, using MSLD temp instead
       !GG: I/J here refer to before/after swap
       P=MIN(ONE,EXP(-(ONE/(KBOLTZ*TBLD))*(ENEIGHJ+EPOTJ-ENEIGHI-EPOTI)))
       !P=ZERO
    ELSE
       P=MIN(ONE,EXP(-(ONE/(KBOLTZ*TEMNEW))*(ENEIGHJ+EPOTJ-ENEIGHI-EPOTI)))
    ENDIF
    RN=RANDOM(REPSEED)
    QEXC=P.GT.RN

#if KEY_PHMD==1
    IF (QPHRX) THEN
111 format(a16,a4,i5,a9,i5,a5,i10) 
112 format(a15,f16.4,a15,f16.4)
113 format(a15,f16.4,a15,f16.4)
         IF (QEXC) THEN
            write(iunrex,'(a)') &
            '-------------------------------------------------------------'
            write(iunrex,111) &
            ' PH-REX> ACCEPT ','REP',ME,' NEIGHBOR',CNEIGHBOR,'STEP',ISTART-1
         ELSE 
            write(iunrex,'(a)') &
            '-------------------------------------------------------------'
            write(iunrex,111) &
            ' PH-REX> REJECT ','REP',ME,'NEIGHBOR',CNEIGHBOR,'STEP',ISTART-1
         ENDIF
         write(iunrex,112) ' Ei(REP)=',EPOTI,' Ej(Neighbor)=',ENEIGHI
         write(iunrex,113) ' Ei(Neighbor)=',EPOTJ,'Ej(REP)=',ENEIGHJ       
     ELSEIF ( .NOT. QMSPHRX) THEN !GG: CPHMD^MSLD PH-REX printout is written out later
     ELSE
#endif
1111     format('H-REX> REPL ',I3,' NEIGHBOR ',I3,' STEP ',I8)
1112     format('H-REX> EPOT1(P) ',f12.2,' EPOT1(Q) ',f12.2,' EPOT2(Q) ',f12.2,' EPOT2(P) ',f12.2)
         if(prnlev.gt.0) then
            write(IUNREX,'(A)') '------------------'
            write(IUNREX,1111) irepdstr,neighbor,istart-1
            write(IUNREX,1112) epoti,epotj,eneighi,eneighj
         endif

#if KEY_PHMD==1
    ENDIF
#endif
    !
    !     QEXC would be the result of the probability test...
    IF(QEXPT.OR.QEX2D)THEN
       IF((STEP.EQ.1).OR.(STEP.EQ.2))THEN
          IF(MOD(IREPDSTR,2).EQ.0)THEN
             CALL GREC(CNEIGHBOR,3,QEXC,4)
          ELSE
             CALL GSEN(CNEIGHBOR,3,QEXC,4)
          ENDIF
       ELSE
          !         ALLOCATE(W(1))
          !         CALL GRECSEN(CNEIGHBOR,9,W,1,QEXC,1)
          IF(MEG.EQ.0)THEN

             CALL GREC(CNEIGHBOR,3,QEXC,4)
          ELSE
             CALL GSEN(CNEIGHBOR,3,QEXC,4)
          ENDIF
       ENDIF
    ELSE
       IF(MOD(IREPDSTR,2).EQ.0)THEN
          CALL GREC(CNEIGHBOR,3,QEXC,4)
          CALL GREC(CNEIGHBOR,9,RN,8)
       ELSE
          CALL GSEN(CNEIGHBOR,3,QEXC,4)
          CALL GSEN(CNEIGHBOR,9,RN,8)
       ENDIF
    ENDIF

    IF(NEIGHBOR.GT.IREPDSTR) THEN
       NOPPUP=NOPPUP+1
       IF(QEXC) NSUCUP=NSUCUP+1
    ELSE
       NOPPDN=NOPPDN+1
       IF(QEXC) NSUCDN=NSUCDN+1
    ENDIF

#if KEY_BLOCK==1
    !GG: CPHMD^MSLD PH-REX printout done after the Call GREC/GSEN(CNEIGHBOR,3,QEXC,4) commands
    !GG  If done before, QEXC will not be consistent since not updated?
    IF (QMSPHRX) THEN
711 format(a16,a4,i5,a9,i5,a5,i10)     !!GG: writes out exchange results
712 format(a15,f16.4,a15,f16.4)
713 format(a15,f16.4,a15,f16.4)
714 format(a8,f18.14)
         IF (QEXC) THEN
            write(outu,'(a)') &
            '------------Replica Exchange (PH-REX) for CPHMD^MSLD------------'
            write(outu,711) &
            ' PH-REXv2> ACCEPT ','REP', IREPDSTR,' NEIGHBOR',NEIGHBOR,'STEP',ISTART-1
         ELSE
            write(outu,'(a)') &
            '------------Replica Exchange (PH-REX) for CPHMD^MSLD------------'
            write(outu,711) &
            ' PH-REXv2> REJECT ','REP', IREPDSTR,' NEIGHBOR',NEIGHBOR,'STEP',ISTART-1
         ENDIF
         write(outu,712) ' Ei(REP)=',EPOTI,' Ej(Neighbor)=',ENEIGHI  !GG: I/J here refers to pH
         write(outu,713) ' Ei(Neighbor)=',EPOTJ,'Ej(REP)=',ENEIGHJ
         write(outu,714) ' Ediff', EDIFF
         write(outu,714) ' Prob', P
         write(outu,*) ' Success', QEXC
    ENDIF
#endif

    !
    !     If acceptance, keep the test exchange of coordinate;if rejection,
    !     perform the exchange again(go back to the configuration before exchange):
999 CONTINUE
    CALL PSND4(QEXC,1)
    CALL PSND4(NOPPUP,1)
    CALL PSND4(NOPPDN,1)
    CALL PSND4(NSUCUP,1)
    CALL PSND4(NSUCDN,1)
    CALL PSND8(P,1)

1113    format('H-REX> PROB ',f6.4,' EXCH ',L1)

    if(prnlev.gt.0) then
       write(IUNREX,1113) p,qexc
       call flush(iunrex)
    endif

    IF(QEXC)THEN   !GG: ------------------- IF EXCHANGE WAS ACCEPTED -------------------------
       !     exchange velocities - more stable than assigning new velocities

       ! BTM -- this code is a gory mess, but since we've exchanged velocities, we need
       ! to restart the dynamics algorithm
       JHSTRT=0
       IGVOPT=2
       IF(MYNOD.NE.0) GOTO 888
       call chmalloc('repdstr.src','REPXCHGL','W',NATOM,crl=w)
       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,VX,NATOM)
       vx(1:natom)=w(1:natom)
       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,VY,NATOM)
       vy(1:natom)=w(1:natom)
       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,VZ,NATOM)
       vz(1:natom)=w(1:natom)
       call chmdealloc('repdstr.src','REPXCHGL','W',NATOM,crl=w) 
       if(irepdstr.gt.neighbor) then
          call grec(cneighbor,5,reptag,4)
          call gsen(cneighbor,6,oldrep,4)
       else
          call gsen(cneighbor,5,oldrep,4)
          call grec(cneighbor,6,reptag,4)
       endif
#if KEY_PHMD==1
       IF(QPHRX)THEN ! JAW. EXCHANGE THETA VELOCITY TERMS, EXCHANGE WAS ACCEPTED
          CALL GRECSEN(CNEIGHBOR,2,L,NTITR,VPH_THETA,NTITR)
          VPH_THETA(1:ntitr) = l(1:ntitr)
          CALL GRECSEN(CNEIGHBOR,2,L,NTITR,THETAOLD,NTITR)
          THETAOLD(1:ntitr) = l(1:ntitr)
          IF(PHBETA .gt. ZERO) THEN ! Langevin PHMD
             CALL GRECSEN(CNEIGHBOR,2,L,NTITR,VPHOLD,NTITR)
             VPHOLD(1:ntitr) = l(1:ntitr)
             CALL GRECSEN(CNEIGHBOR,2,L,NTITR,DPHOLD,NTITR)
             DPHOLD(1:ntitr) = l(1:ntitr)
          ENDIF
       ENDIF
#endif
#if KEY_BLOCK==1
       IF (QMSPHRX) THEN
          IF(prnlev.ge.6) THEN
             write(outu,'(a)') "EXCHANGE ACCEPTED, not transferring theta variables back" !GG: All theta variables prev transferred
          ENDIF
       ENDIF
       !
706 format(a2,i4,a5,f20.8)
#endif
888    CONTINUE
       CALL PSND8(VX,NATOM)
       CALL PSND8(VY,NATOM)
       CALL PSND8(VZ,NATOM)
       CALL PSND4(REPTAG,1)
    ELSE   !GG: ------------------- IF EXCHANGE WAS REJECTED -------------------------

#if KEY_PHMD==1
       !JAW - Exchange was rejected--> we need to reset the non-bond lists. also an energy call resets GB. THIS IS NECESSARY!!! 
       IF(QPHRX)THEN
          CALL UpdatePHMD(1, 1, 0, 1) ! After Exchanging Coordinates and Lambda it is necesarry to update charges and
                                      ! total system charge so energy calculations are correct. JW

          !???cb3 are the next two calls supposed to be inside te if(qphrx?)
          !???btm I believe that they are...
          CALL NBONDS(X,Y,Z,BNBND,BIMAG)  ! LNBND,BIMAG,LIMAG)
          CALL ENERGY(X,Y,Z,DX,DY,DZ,BNBND,BIMAG,0)
       ENDIF
#endif

       !JAW
       IF (MYNOD.NE.0) GOTO 9999
       call chmalloc('repdstr.src','REPXCHGL','W',NATOM,crl=w)
       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,X,NATOM)
       x(1:natom) = w(1:natom)
       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,Y,NATOM)
       y(1:natom) = w(1:natom)
       CALL GRECSEN(CNEIGHBOR,2,W,NATOM,Z,NATOM)
       z(1:natom) = w(1:natom)
       call chmdealloc('repdstr.src','REPXCHGL','W',NATOM,crl=w)

#if KEY_PHMD==1
       IF(QPHRX)THEN ! JW: SWITCH THETA VALUES BACK, EXCHANGE WAS REJECTED 
          CALL GRECSEN(CNEIGHBOR,2,L,NTITR,PH_THETA,NTITR)
          PH_THETA(1:ntitr) = l(1:ntitr)
       ENDIF 
#endif
#if KEY_BLOCK==1
       IF (QMSPHRX) THEN
          IF(prnlev.ge.6) THEN
             write(outu,'(a)') "EXCHANGE REJECTED, transferring theta variables back"
          ENDIF
          K = nsitemld*nblock                          !GG: Calc total no. of elements in MSLD array
          THETAMLDS = RESHAPE(THETAMLD,SHAPE(M))       !GG: Reshape array for coordinates to 1D array
          CALL GRECSEN(CNEIGHBOR,2,M,K,THETAMLDS,K)    !GG: Transfer array
          THETAMLDS(1:K) = M(1:K)
          THETAMLD = RESHAPE(THETAMLDS,SHAPE(N))       !GG: Reshape array back to nD array
          THETAMLDOLDS = RESHAPE(THETAMLDOLD,SHAPE(M)) !GG: Processing previous theta coordinates
          CALL GRECSEN(CNEIGHBOR,2,M,K,THETAMLDOLDS,K)
          THETAMLDOLDS(1:K) = M(1:K)
          THETAMLDOLD = RESHAPE(THETAMLDOLDS,SHAPE(N))
          THETAVMLDS = RESHAPE(THETAVMLD,SHAPE(M))     !GG: Processing current velocity coordinates
          CALL GRECSEN(CNEIGHBOR,2,M,K,THETAVMLDS,K)
          THETAVMLDS(1:K) = M(1:K)
          THETAVMLD = RESHAPE(THETAVMLDS,SHAPE(N))
          THETAFMLDS = RESHAPE(THETAFMLD,SHAPE(M))     !GG: Processing current force coordinates
          CALL GRECSEN(CNEIGHBOR,2,M,K,THETAFMLDS,K)
          THETAFMLDS(1:K) = M(1:K)
          THETAFMLD = RESHAPE(THETAFMLDS,SHAPE(N))
       ENDIF
#endif

    !Exchange cell dimension information
      IF (QCRYS .AND. XDIM.GT.0) THEN
         call chmalloc('repdstr.src','REPXCHGL','W',6,crl=w)
         CALL GRECSEN(CNEIGHBOR,4,W,6,XTLABC,6)
         XTLABC(1:6) = W(1:6)
         call chmdealloc('repdstr.src','REPXCHGL','W',6,crl=w)
      ENDIF
      call chmalloc('repdstr.src','REPXCHGL','W',NATOM,crl=w)
      CALL GRECSEN(CNEIGHBOR,11,W,NATOM,WMAIN,NATOM)
      WMAIN(1:natom) = w(1:natom)
      call chmdealloc('repdstr.src','REPXCHGL','W',NATOM,crl=w)

#if KEY_GCMC==1
       IF (QGCMC) THEN
          call chmalloc('repdstr.src','REPXCHGL','WL',MAXA,log=wl)
          CALL GRECSEN(CNEIGHBOR,9,W,MAXA,GCMCON,MAXA)
          gcmcon(1:maxa) = WL(1:MAXA)
          call chmdealloc('repdstr.src','REPXCHGL','WL',MAXA,log=wl)
       ENDIF
#endif
9999   CONTINUE
#if KEY_BLOCK==1
       call msld_checkvariables(2) !GG: Before transmission
#endif
       !     We need to broadcast the data within the same replica group
       CALL PSND8(X,NATOM)
       CALL PSND8(Y,NATOM)
       CALL PSND8(Z,NATOM)
       CALL PSND8(VX,NATOM)
       CALL PSND8(VY,NATOM)
       CALL PSND8(VZ,NATOM)
       CALL PSND8(WMAIN,NATOM)
#if KEY_GCMC==1
       IF (QGCMC) THEN
          CALL PSND4(GCMCON,MAXA)
       ENDIF
#endif
      IF (QCRYS) THEN
         CALL PSND8(XTLABC,6)
         CALL XTLLAT(XUCELL,XTLABC)
         call chmalloc('repdstr.src','REPXCHGL','TRANSF',3,4,XNSYMM, &
                        crl=TRANSF)
         CALL IMFILL(TRANSF,.FALSE.)
         call chmdealloc('repdstr.src','REPXCHGL','TRANSF',3,4,XNSYMM, &
                          crl=TRANSF)
         !!CALL UPIMAG(X,Y,Z,WMAIN,0,X,Y,Z,VX,VY,VZ)
      ENDIF

#if KEY_BLOCK==1
      IF (QMSPHRX) THEN   !GG: Broadcast swapped data to other replicas
         K = nsitemld*nblock
         CALL PSND8(THETAMLDS,K)
         CALL PSND8(THETAMLDOLDS,K)
         CALL PSND8(THETAVMLDS,K)
         CALL PSND8(THETAFMLDS,K)
         IF (MYNOD.EQ.0) GOTO 992
         THETAMLD = RESHAPE(THETAMLDS,SHAPE(N))
         THETAMLDOLD = RESHAPE(THETAMLDOLDS,SHAPE(N))
         THETAVMLD = RESHAPE(THETAVMLDS,SHAPE(N))
         THETAFMLD = RESHAPE(THETAFMLDS,SHAPE(N))
992      CONTINUE
     ENDIF
     IF (QMSPHRX) THEN
         call msld_swapcoeff_mod      !GG: Toggle msld_setblcoef_fnexp to use thetamold values
         call msld_checkvariables(3) !GG: Variables  After Swap, Before Call Energy
     ENDIF
     CALL ENERGY(X,Y,Z,DX,DY,DZ,BNBND,BIMAG,0)
     IF (QMSPHRX) THEN
        call msld_swapcoeff_norm      !GG: Toggle msld_setblcoef_fnexp to use thetamold values
        call msld_checkvariables(4) !GG: Variables  After Swap, After Call Energy
     ENDIF
#endif

      !!CALL NBONDS(X,Y,Z,BNBND,bimag)   !LNBND,BIMAG,LIMAG)

    ENDIF

    IF(QCRYS) CALL UPIMAG(X,Y,Z,WMAIN,0,X,Y,Z,VX,VY,VZ)
    CALL NBONDS(X,Y,Z,BNBND,bimag)
    CALL ENERGY(X,Y,Z,DX,DY,DZ,BNBND,BIMAG,0)
#if KEY_BLOCK==1
    IF (QMSPHRX) THEN
       call msld_checkvariables(5) !GG: Variables  After Swap, After MC Exchange
    ENDIF
#endif
222 CONTINUE

1114    format('SUMHREX> ',I8,' NBR ',I2,' ORTAG ',I2,' NWTAG ',I2,' epot1(p) ',F12.2,' epot1(q) ',F12.2, &
               ' epot2(q) ',F12.2,' epot2(p) ',f12.2,' prob ',f5.3,' rn ',f5.3,' exc ',L1)
1115    format('HREX UP> OPPORTUNITIES ',I8,' SUCCESSES ',I8,' RATIO ',F6.4)
1116    format('HREX DN> OPPORTUNITIES ',I8,' SUCCESSES ',I8,' RATIO ',F6.4)
1117    format(I8,x,ES25.16,x,ES25.16,x,ES25.16,x,ES25.16)

    if(noppup > 0) exratup=real(nsucup)/real(noppup)
    if(noppdn > 0) exratdn=real(nsucdn)/real(noppdn)
    if((noppup+noppdn.gt.0).and.lexattempt) then
       reprob=(min(p,one)+((noppup+noppdn-1)*reprob))/float(noppup+noppdn)
       if(prnlev > 5) then
          write(iunrex,'(a,4i3)') 'REPEXCHL> noppup, noppdn, nsucup, nsucdn = ',noppup,noppdn,nsucup,nsucdn
          write(iunrex,'(a,f10.5)') 'REPEXCHL> P = ',min(p,one)
          write(iunrex,'(a,3f6.3)') 'REPEXCHL> EXRUP EXRDN REPROB = ',exratup,exratdn,reprob
       endif
    endif

    call set_param('EXRUP',exratup)
    call set_param('EXRDN',exratdn)
    call set_param('REPROB',reprob)

    if(prnlev.gt.0) then
       if(cneighbor > mynodg .and. ewritu > 0) &
          write(EWRITU,1117) istart-1,epoti,eneighj,epotj,eneighi

       write(IUNREX,1114) istart-1,neighbor,oldrep,reptag,epoti,epotj,eneighi,eneighj,p,rn,qexc
       if(qsump) then
          if(noppup.gt.0) then
             write(IUNREX,1115) noppup,nsucup,exratup
          endif
          if(noppdn.gt.0) then 
             write(IUNREX,1116) noppdn,nsucdn,exratdn
          endif
          if(nsucup.eq.0.and.nsucdn.eq.0) then
             write(IUNREX,'(a)') 'WARNING! NO SUCCESSFUL EXCHANGES'
          else if(nsucup.gt.0) then
             rerat=exratdn/exratup
             if(noppdn.gt.0.and.(rerat.gt.2.0.or.rerat.lt.0.5)) then
                write(IUNREX,'(a)') 'WARNING! UNBALANCED EXCHANGE RATIO; CONSIDER MOVING REPLICA OVERLAPS!'
             endif
          else if(nsucdn.gt.0) then ! defensive programming, FTW!
             rerat=exratup/exratdn
             if(noppup.gt.0.and.(rerat.gt.2.0.or.rerat.lt.0.5)) then
                write(IUNREX,'(a)') 'WARNING! UNBALANCED EXCHANGE RATIO; CONSIDER MOVING REPLICA OVERLAPS!'
             endif
          endif
       endif
    endif

    RETURN
  END SUBROUTINE REPEXCHGL

  !
  SUBROUTINE DREPSETIO(IOL,PRNL,WRNL)
    !-----------------------------------------------------------------------
    !     Set the IOLEV information from global values, ie all system
    !
  use parallel
  use repdstr
    !
    INTEGER IOL,PRNL,WRNL
    !
    !
    !??? not used here anymore      IF (IOON.EQ.1) RETURN
    !
    !     this must be here since it is easier to protect it here then
    !     in the calling routines
    !
    IF (.NOT.QREPDSTR) RETURN
    !
    IOON=1
    !
    IOLORIG=IOL
    PRNLORIG=PRNL
    WRNLORIG=WRNL
    !      write(70+mynodg,'(a,5i6)')'DREPSETIO-0>me,meg,iolev=',
    !     $     mynod,mynodg,iol,prnl
    !
    !     There is no problem here with IOLEV, ie always like this
    !     even when QRDQTT=.TRUE., since we deal with this separately
    IF(MYNOD.EQ.0)IOL=1
    !
    !     For QWRQTT=.TRUE. we can deal with it here:
    IF(QWRQTT.AND.(MYNOD.EQ.0))PRNL=5
    IF(QWRQTT.AND.(MYNOD.EQ.0))WRNL=5
    !      
    !      write(70+mynodg,'(a,5i6)')'DREPSETIO-1>me,meg,iolev=',
    !     $     mynod,mynodg,iol,prnl
    !
    RETURN
  END SUBROUTINE DREPSETIO
  !
  SUBROUTINE DREPRESIO(IOL,PRNL,WRNL)
    !-----------------------------------------------------------------------
    !     Restores the IOLEV information from global values, ie all system
    !
    !     Generalize this too!!!
    !
  use parallel
  use repdstr
    !
    INTEGER IOL,PRNL,WRNL
    !
    !! not used anymore:       IF (IOON.EQ.0) RETURN
    IF (.NOT.QREPDSTR) RETURN
    !
    IOON=0
    !      write(50+mynodg,*)'DREPRESIO-0>me,meg,iolev=',
    !     $     mynod,mynodg,iol,prnl
    IOL=IOLORIG
    PRNL=PRNLORIG
    WRNL=WRNLORIG
    !      write(50+mynodg,*)'DREPRESIO-1>me,meg,iolev=',
    !     $     mynod,mynodg,iol,prnl
    !
    RETURN
  END SUBROUTINE DREPRESIO
  !
  SUBROUTINE PSETGLOB
    !-----------------------------------------------------------------------
    !     Set the parallel information from global values, ie all system
    !
  use parallel
    !
    integer:: i
    MYNOD=MYNODG
    NUMNOD=NUMNODG
    MYNODP=MYNOD+1
    NODDIM=NPTWO()
    CALL CUBE(MYNOD,NUMNOD,IPPMAP)
    !
    !     Lets do the INODE array, too:
    !
    DO I=1,MAXNODE
       INODE(I)=MOD(I,NUMNOD)
    ENDDO
    !
    RETURN
  END SUBROUTINE PSETGLOB
  !
  SUBROUTINE PSETLOC
    !-----------------------------------------------------------------------
    !     Set the local parallel information.
    !     Put PSETLOC,PSETGLOB into paral1.src ??
    !
  use parallel
  use repdstr
    !
    INTEGER I
    !
    NUMNOD=NUMNODG/NREPDSTR
    MYNOD=MOD(MYNODG,NUMNOD)
    MYNODP=MYNOD+1
    NODDIM=NPTWO()
    DO I = 0, NUMNOD
       IPPMAP(I)=I+NUMNOD*IREPDSTR
    ENDDO
    !
    !     Lets do the INODE array, too:
    !
    DO I=1,MAXNODE
       INODE(I)=MOD(I,NUMNOD)
    ENDDO
    !
    !     Also iparpt(), nparpt() have to be localized??? !!!! Pressure problems???
    !     Maybe not so trivial??? Better save
    !     both cases when generated and then copy them as needed ???
    !  IPARPT and NPARPT are more tricky, since they depend on
    !  groups. It can't be dealt generally since when repd is started we
    !  might don't have the PSF yet. It is maybe better to use separate routine
    !  and call it before vdgsum & vdgbr, or whatever routines use them.
    !  or use nbonds command in the repd script after PSFs are set up.
    !

    RETURN
  END SUBROUTINE PSETLOC

  SUBROUTINE RESEXCH(QHIGH,X,Y,Z,VX,VY,VZ,QEXC,CURTMP,CURENE,ISEED,IASVEL,IGVOPT, &
                     P,RN,ENEIGH,FNEIGH &
#if KEY_TSM==1
                     ,BACKLS & 
#endif
                    )
     use psf
     use clcg_mod,only: random
     use stream
     use consta,only: kboltz
     use number,only: one
     use repdstr
     use energym
     use bases_fcm
     use memory
     use dynutil, only: assvel

     LOGICAL, INTENT(IN)             :: QHIGH
     LOGICAL, INTENT(OUT)            :: QEXC
     INTEGER, INTENT(OUT)            :: FNEIGH
     REAL(CHM_REAL), INTENT(OUT)     :: P,RN,ENEIGH
     REAL(CHM_REAL), INTENT(INOUT)   :: X(*),Y(*),Z(*),VX(*),VY(*),VZ(*)
     REAL(CHM_REAL), INTENT(INOUT)   :: CURTMP,CURENE
     INTEGER,INTENT(INOUT)           :: ISEED,IASVEL,IGVOPT
#if KEY_TSM==1
     INTEGER,INTENT(INOUT)           :: BACKLS(*) 
#endif

     INTEGER                                 :: UNUM,RESSZ,TRGT,I
     REAL(CHM_REAL)                          :: C,S,NEWE,NEWT,TMP,PCUR,PGEN,DCUR,DGEN
     REAL(CHM_REAL),DIMENSION(NATOM)         :: RX,RY,RZ,ORX,ORY,ORZ,RDX,RDY,RDZ
     REAL(CHM_REAL),ALLOCATABLE,DIMENSION(:) :: DXYZ
     REAL(CHM_REAL)                          :: EBACKUP

     IF(QHIGH) THEN
        UNUM=RHUNIT
        RESSZ=HIGHRESSZ
     ELSE
        UNUM=RLUNIT
        RESSZ=LOWRESSZ
     ENDIF
     TRGT=CEILING(RANDOM(REPSEED)*RESSZ)
     IF(TRGT.EQ.0) TRGT=1 ! Just in case the RNG returns 0
     FNEIGH=TRGT

     ! Decide if we want to swap with the reservoir
     QEXC = .FALSE.
     IF(QRESBOLTZ) THEN
        IF(QHIGH) THEN
           NEWT=RHTEMP        
           NEWE=RHENER(TRGT)
        ELSE
           NEWT=RLTEMP
           NEWE=RLENER(TRGT)
        ENDIF

        IF(CURTMP.EQ.NEWT) THEN
           P=MIN(ONE,EXP(-(NEWE-CURENE)/(KBOLTZ*CURTMP)))
        ELSE
           P=MIN(ONE,EXP(-(ONE/(KBOLTZ*CURTMP) &
            -ONE/(KBOLTZ*NEWT))*(NEWE-CURENE)))
        ENDIF
        RN=RANDOM(REPSEED)
        QEXC=RN.LT.P
        ENEIGH=NEWE
     ELSE IF(QRESNOBO) THEN
        IF(QHIGH) THEN
           NEWE=RHENER(TRGT)
        ELSE
           NEWE=RLENER(TRGT)
        ENDIF
        P=MIN(ONE,EXP(-(ONE/(KBOLTZ*CURTMP))*(NEWE-CURENE)))
        RN=RANDOM(REPSEED)
        QEXC=RN.LT.P
        ENEIGH=NEWE
     ENDIF

     ! go ahead and swap in the new atomic coordinates, saving
     ! the old coordinates in rescrdx,y,z for possible adding to 
     ! the reservoir
     IF(QEXC) THEN

        IF(QHIGH) THEN
           NSUCUP=NSUCUP+1
        ELSE
           NSUCDN=NSUCDN+1
        ENDIF

        ! actually read the new coordinates from the file
        I=((TRGT-1)*3)+1
        READ(UNUM,REC=I)   RESCRDX
        READ(UNUM,REC=I+1) RESCRDY
        READ(UNUM,REC=I+2) RESCRDZ

        IF(PRNLEV.GE.6) &
          WRITE(IUNREX,'(A,I9,A)') '-- RESEXCH: SWAPPED IN ELT ', TRGT, ' OF THE RESERVOIR --'
        DO I=1,NATOM
           IF(PRNLEV.GE.7) &
              WRITE(IUNREX,'(A,I5,A,3F10.6)') 'RESEXCH> CRD ', I, ' = ', RESCRDX(I), RESCRDY(I), RESCRDZ(I)
           ! X
           TMP=X(I)
           X(I)=RESCRDX(I)
           RESCRDX(I)=TMP
           ! Y
           TMP=Y(I)
           Y(I)=RESCRDY(I)
           RESCRDY(I)=TMP
           ! Z
           TMP=Z(I)
           Z(I)=RESCRDZ(I)
           RESCRDZ(I)=TMP
        ENDDO

        ! Since we've reassigned the coordinates, we need to call assvel to
        ! reassign the velocities
        CALL ASSVEL(CURTMP,X,Y,Z,VX,VY,VZ,AMASS,ISEED,IASVEL,IGVOPT,NATOM,IMOVE &
#if KEY_TSM==1
                    ,BACKLS & 
#endif
                   )
       
    ELSE IF(PRNLEV.GE.6) THEN
        WRITE(IUNREX,'(A,I9,A)') '--- RESEXCH: SWAP WITH RESERVOIR NOT ACCEPTED ---'
    ENDIF

  END SUBROUTINE RESEXCH

  SUBROUTINE RESEXCHL(QHIGH,X,Y,Z,VX,VY,VZ,TEMP,RESTEMP,OUREPOT,ISEED,IASVEL,IGVOPT,P,STEP,JHSTRT &
#if KEY_TSM==1
                             ,BACKLS & 
#endif
                     )

     use psf
     use consta
     use number
     use repdstr
     use parallel
     use stream
     use clcg_mod,only: random
     use memory
     use dynutil, only: assvel

     ! passed in variables
     logical, intent(in)          :: qhigh
     integer, intent(in)          :: iseed,iasvel,step
     integer, intent(inout)       :: igvopt,jhstrt,backls(*)
     real(chm_real),intent(in)    :: temp,restemp,ourepot
     real(chm_real),intent(inout) :: x(*),y(*),z(*),vx(*),vy(*),vz(*)
     real(chm_real),intent(out)   :: p

     ! local variables
     integer                                 :: unum,ressz,trgt,i
     logical                                 :: qexc
     real(chm_real)                          :: resepot,rn

     if(mynod.eq.0) then
        ! OK, we have a problem here. If we're replica I w/ cords Q and the reservoir is replica
        ! J w/ coords Q', it's easy to calculate E_I(Q) and E_I(Q'), but it's not so easy to
        ! calculate E_J(Q) and E_J(Q'), since we don't have an explicit replica to ship the coordinates
        ! off to. One idea is to treat the reservoir as Boltzmann, which is what I've implemented.

        if(qhigh) then
           unum=rhunit
           ressz=highressz
        else
           unum=rlunit   
           ressz=lowressz
        endif
        trgt=ceiling(random(iseed)*ressz)
        if(trgt.eq.0) trgt=1 ! Just in case the RNG returns 0
        if(prnlev.gt.3) write(iunrex,'(a,i6)') 'RESEXCHL> TRGT = ', trgt
        call flush(iunrex)
        if(qhigh) then
           resepot=rhener(trgt)
        else
           resepot=rlener(trgt)
        endif

        ! Decide if we want to swap with the reservoir
        rn=random(repseed)
        if(temp.eq.restemp) then
           p=one
        else
           p=min(one,exp(-(one/(kboltz*temp) &
                -one/(kboltz*restemp))*(ourepot-resepot)))
        endif
        qexc=p.gt.rn

        if(qexc) then
           ! do the swap
           i=(3*trgt)-2

           read(unum,rec=i)   rescrdx
           read(unum,rec=i+1) rescrdy
           read(unum,rec=i+2) rescrdz
           x(1:natom)=rescrdx(1:natom)
           y(1:natom)=rescrdy(1:natom)
           z(1:natom)=rescrdz(1:natom)
        endif
     endif

     call psnd4(qexc,1)
     if(qexc) then     
        if(qhigh) then
           nsucup=nsucup+1
        else
           nsucdn=nsucdn+1
        endif

       jhstrt=0
       igvopt=2
       
       call psnd8(x,natom)
       call psnd8(y,natom)
       call psnd8(z,natom)
       call assvel(temp,x,y,z,vx,vy,vz,amass,iseed,iasvel,igvopt,natom,imove &
#if KEY_TSM==1
                   ,backls & 
#endif
                   )
     endif

     ! almost done -- report output to the user
123  format('H-R-REX> ',i8,' REP ENE ',f10.4,' RES ENE ',f10.4,' PROB ',f6.4,' RN ',f6.4,' SUCCESS? ',L1)
     if(prnlev.gt.3) then
        write(iunrex,123) step,ourepot,resepot,p,rn,qexc
        call flush(iunrex)
     endif

  END SUBROUTINE RESEXCHL


  SUBROUTINE RESEXCH_PH(CURTMP,QHIGH,X,Y,Z,VX,VY,VZ,OURPH,IPROTO, &
                        ISEED,IASVEL,IGVOPT,P,RN,JPROTO,PH_M,QEXC,RESNUM &
#if KEY_TSM==1
                        ,BACKLS & 
#endif
                       )

     use psf
     use clcg_mod,only: random
     use stream
     use consta,only: kboltz
     use number,only: one
     use consph,only: tstate 
     use repdstr
     use image
     use bases_fcm
     use dynutil, only: assvel

     LOGICAL, INTENT(IN)            :: QHIGH
     LOGICAL, INTENT(OUT)           :: QEXC
     INTEGER, INTENT(IN)            :: IPROTO
     INTEGER, INTENT(OUT)           :: RESNUM,JPROTO
     REAL(CHM_REAL), INTENT(OUT)    :: P,RN,PH_M
     REAL(CHM_REAL), INTENT(INOUT)  :: X(*),Y(*),Z(*),VX(*),VY(*),VZ(*)
     INTEGER,INTENT(INOUT)          :: ISEED,IASVEL,IGVOPT
#if KEY_TSM==1
     INTEGER,INTENT(INOUT)          :: BACKLS(*) 
#endif
     REAL(CHM_REAL), INTENT(IN)     :: CURTMP,OURPH

     INTEGER                        :: UNUM,RESSZ,TRGT,I,J,K,PIDX
     REAL(CHM_REAL)                 :: PH_DELTA,TMP
     INTEGER,DIMENSION(NATOM)       :: TMPTSTATE
     INTEGER,DIMENSION(:,:),POINTER :: RESTSTATE

     IF(QHIGH) THEN
        UNUM=RHUNIT
        RESSZ=HIGHRESSZ
        PH_M=RHPH
        RESTSTATE=>RESHTSTATE
     ELSE
        UNUM=RLUNIT
        RESSZ=LOWRESSZ
        PH_M=RLPH
        RESTSTATE=>RESLTSTATE
     ENDIF
     TRGT=CEILING(RANDOM(REPSEED)*RESSZ)
     IF(TRGT.EQ.0) TRGT=1 ! Just in case the RNG returns 0
     RESNUM=TRGT

     ! Decide if we want to swap with the reservoir
     QEXC = .FALSE.

     ! For the Boltzmann case at least, we need to read in the reservoir file
     ! so we can calculate JPROTO.
     I=((TRGT-1)*5)+1
     READ(UNUM,REC=I)   RESCRDX
     READ(UNUM,REC=I+1) RESCRDY
     READ(UNUM,REC=I+2) RESCRDZ
     READ(UNUM,REC=I+3) RESCG
     READ(UNUM,REC=I+4) TMPTSTATE
     JPROTO=0
     DO J=1,NRES
        IF(TMPTSTATE(J).EQ.1) JPROTO=JPROTO+1
     ENDDO
     !!!WRITE(OUTU,'(A,I5)') 'TIM DEBUG> COORDINATES OF RESERVOIR STRUCTURE ',TRGT
     !!!DO J=1,NATOM
     !!!   WRITE(OUTU,'(3F8.3)') RESCRDX(J),RESCRDY(J),RESCRDZ(J)
     !!!ENDDO
     WRITE(OUTU,'(A,I5)') 'TIM DEBUG> PROTONATION COUNT OF RESERVOIR STRUCT = ', JPROTO
     CALL FLUSH(OUTU)

     IF(QRESBOLTZ) THEN
        PH_DELTA=LOG(10.0)*(PH_M-OURPH)*(IPROTO-JPROTO)
        IF(PH_DELTA.LE.0) THEN
           P=ONE
        ELSE
           P=MIN(ONE,EXP(-PH_DELTA))
        ENDIF
     ELSE IF(QRESNOBO) THEN
        CALL WRNDIE(-3,'<RESEXCH_PH>','NONBOLTZMANN REX NOT IMPLEMENTED FOR PH')
     ENDIF
     RN=RANDOM(REPSEED)
     QEXC=RN.LT.P

     IF(QEXC) THEN
        IF(PRNLEV.GE.6) &
          WRITE(IUNREX,'(A,I9,A)') '-- RESEXCH: SWAP IN ELT ', TRGT, ' OF THE RESERVOIR --'


        IF(QHIGH) THEN
           NSUCUP=NSUCUP+1
        ELSE
           NSUCDN=NSUCDN+1
        ENDIF

        DO I=1,NATOM
           !!WRITE(IUNREX,'(A,I5,A,4XF10.3)') 'RESEXCH> CRD ', I, ' = ', RESCRDX(I), RESCRDY(I), RESCRDZ(I), RESCG(I)
           ! X
           TMP=X(I)
           X(I)=RESCRDX(I)
           RESCRDX(I)=TMP
           ! Y
           TMP=Y(I)
           Y(I)=RESCRDY(I)
           RESCRDY(I)=TMP
           ! Z
           TMP=Z(I)
           Z(I)=RESCRDZ(I)
           RESCRDZ(I)=TMP
           ! CHARGE
           TMP=CG(I)
           CG(I)=RESCG(I)
           RESCG(I)=TMP

           IF(CG(I).NE.RESCG(I)) THEN
              if(ntrans > 0) then
                 do k=natom+1,natim
                    pidx=bimag%IMATTR(k)
                    if(pidx > 0) then
                       cg(k)=cg(i)
                    endif
                 enddo
              endif
           ENDIF
        ENDDO
        ! now we need to set the protonation states of each of the residues
        ! (via the tstate array)
        DO I=1,NRES
           TSTATE(I)=TMPTSTATE(I)
           !!WRITE(IUNREX,'(A,I5,A,I3)') 'RESEXCH> SET TSTATE(', I, ') = ', TSTATE(I)
        ENDDO

        ! Since we've reassigned the coordinates, we need to call assvel to
        ! reassign the velocities
        CALL ASSVEL(CURTMP,X,Y,Z,VX,VY,VZ,AMASS,ISEED,IASVEL,IGVOPT,NATOM,IMOVE &
#if KEY_TSM==1
                    ,BACKLS & 
#endif
                   )
     ENDIF

  END SUBROUTINE RESEXCH_PH


  SUBROUTINE PRECALCENE(QHIGH,X,Y,Z,QECOR,ECORTEMP)
     use psf
     use number
     use memory
     use stream
     use bases_fcm
     use energym
     use deriv
     use shake
     use consta
     use parallel

     LOGICAL,INTENT(IN)           :: QHIGH,QECOR
     REAL(CHM_REAL),INTENT(IN)    :: ECORTEMP
     REAL(CHM_REAL),INTENT(INOUT) :: X(NATOM),Y(NATOM),Z(NATOM)
     INTEGER                      :: UNUM,RESSZ,I,J,R
     
     REAL(CHM_REAL), ALLOCATABLE, DIMENSION(:) :: TMPX, TMPY, TMPZ
     REAL(CHM_REAL)                            :: CORRECTION
     INTEGER                                   :: NDEGF

     IF(QECOR) THEN
        IF(PRNLEV.GT.1) WRITE(OUTU,'(A,F10.4)') 'PRECALCENE> ADDING ENERGY CORRECTION TERM AT TEMP ', ECORTEMP
        CALL GET_ECOR_ADJUST(CORRECTION,ECORTEMP)
        IF(PRNLEV.GT.1) WRITE(OUTU,'(A,F10.4)') 'PRECALCENE> ENERGY CORRECTION = ', CORRECTION
     ELSE
        NDEGF=-1
        CORRECTION=ZERO
     ENDIF

     IF(QHIGH) THEN
        UNUM=RHUNIT
        RESSZ=HIGHRESSZ
        call chmalloc('repdstr.src','PRECALCENE','rhener',RESSZ,crl=RHENER)
     ELSE
        UNUM=RLUNIT
        RESSZ=LOWRESSZ
        call chmalloc('repdstr.src','PRECALCENE','rlener',RESSZ,crl=RLENER)
     ENDIF

     call chmalloc('repdstr.src','PRECALCENE','tmpx',NATOM,crl=TMPX)
     call chmalloc('repdstr.src','PRECALCENE','tmpy',NATOM,crl=TMPY)
     call chmalloc('repdstr.src','PRECALCENE','tmpz',NATOM,crl=TMPZ)

     DO I=1,NATOM
        TMPX(I)=X(I)
        TMPY(I)=Y(I)
        TMPZ(I)=Z(I)
     ENDDO
     DO I=1,RESSZ
        ! NB, we have to read these into rescrd{x,y,z} because the sdcd
        ! stores them as 4 byte floats, not eight byte...

        RESCRDX(:)=ZERO
        RESCRDY(:)=ZERO
        RESCRDZ(:)=ZERO

        IF(MYNOD.EQ.0) THEN
           R=(I*3)-2
           READ(UNUM,REC=R)   RESCRDX
           READ(UNUM,REC=R+1) RESCRDY
           READ(UNUM,REC=R+2) RESCRDZ

           DO J=1,NATOM
              X(J)=RESCRDX(J)
              Y(J)=RESCRDY(J)
              Z(J)=RESCRDZ(J)
           ENDDO
        ENDIF

        ! Update coordinates and non-bond lists
        ! FIXME: We should call UPDATE here instead
        CALL PSND8(X,NATOM)
        CALL PSND8(Y,NATOM)
        CALL PSND8(Z,NATOM)
        CALL NBONDS(X,Y,Z,BNBND,BIMAG)

        ! Get the energy
        ! The potential energy will be stored in the third element
        ! of the ETERM array
        !CALL GETE(X,Y,Z,X,Y,Z,0)
        CALL ENERGY(X,Y,Z,DX,DY,DZ,BNBND,BIMAG,1)
        IF(PRNLEV.GT.3) WRITE(OUTU,'(A,I6,A,F12.6)') 'PRECALCENE> POTENTIAL OF STRUCT ', &
                                     I, ' = ', EPROP(EPOT)
        !CALL PRINTE(OUTU, EPROP, ETERM, 'ENER', 'ENR', .TRUE., 0, 0, 0, .TRUE.)       

        IF(QECOR) THEN
           IF(PRNLEV.GT.3) WRITE(OUTU,'(A,F12.6)') 'PRECALCENE> CORRECTED ENE = ', EPROP(EPOT)+CORRECTION
           IF(QHIGH) THEN
              RHENER(I)=EPROP(EPOT)+CORRECTION
           ELSE
              RLENER(I)=EPROP(EPOT)+CORRECTION
           ENDIF
        ELSE
           IF(QHIGH) THEN
              RHENER(I)=EPROP(EPOT)
           ELSE
              RLENER(I)=EPROP(EPOT)
           ENDIF
        ENDIF
     ENDDO
     ! Restore original coordinates
     DO I=1,NATOM
        X(I)=TMPX(I)
        Y(I)=TMPY(I)
        Z(I)=TMPZ(I)
     ENDDO
     call chmdealloc('repdstr.src','PRECALCENE','tmpx',NATOM,crl=TMPX)
     call chmdealloc('repdstr.src','PRECALCENE','tmpy',NATOM,crl=TMPY)
     call chmdealloc('repdstr.src','PRECALCENE','tmpz',NATOM,crl=TMPZ)
 
     ! reset the energy and force arrays to whatever they were before
     ! I'm not sure if this matters, but I'm doing it to be safe...
     CALL PSND8(X,NATOM)
     CALL PSND8(Y,NATOM)
     CALL PSND8(Z,NATOM)
     CALL NBONDS(X,Y,Z,BNBND,BIMAG)
     !CALL GETE(X,Y,Z,X,Y,Z,0)
     CALL ENERGY(X,Y,Z,DX,DY,DZ,BNBND,BIMAG,1)

  END SUBROUTINE PRECALCENE

  SUBROUTINE GET_ENE_VAL(ENEUN,QECOR,TEMP,QHIGH)
     use psf
     use shake
     use consta
     use number
     use memory
     use stream

     integer, intent(in)       :: ENEUN
     logical, intent(in)       :: QHIGH,QECOR
     real(chm_real),intent(in) :: TEMP
     integer                   :: i,sz,ndegf
     real(chm_real)            :: correction

     IF(QECOR) THEN
        IF(PRNLEV.GT.1) WRITE(OUTU,'(A,F10.4)') 'GET_ENE_VAL> ADDING ENERGY CORRECTION TERM AT TEMP ', TEMP
        CALL GET_ECOR_ADJUST(CORRECTION,TEMP)
        IF(PRNLEV.GT.4) WRITE(OUTU,'(A,F10.4)') 'GET_ENE_VAL> ENERGY CORRECTION = ', CORRECTION
     ELSE
        NDEGF=-1
        CORRECTION=ZERO
     ENDIF

     if(qhigh) then
        sz=highressz
        call chmalloc('repdstr.src','GET_ENE_VAL','rhener',sz,crl=RHENER)
     else
        sz=lowressz
        call chmalloc('repdstr.src','GET_ENE_VAL','rlener',sz,crl=RLENER)
     endif
     do i=1,sz
        if(qhigh) then
           read(eneun,'(f9.4)') rhener(i)
           rhener(i)=rhener(i)+correction
           if(prnlev.gt.3) write(outu,'(a,i4,a,f10.4)') 'GET_ENE_VAL> ENERGY OF HIGH RES STRUCT ', i, ' = ', rhener(i)
        else
           read(eneun,'(f9.4)') rlener(i)
           rlener(i)=rlener(i)+correction
           if(prnlev.gt.3) write(outu,'(a,i4,a,f10.4)') 'GET_ENE_VAL> ENERGY OF LOW RES STRUCT ', i, ' = ', rlener(i)
        endif
     enddo

  END SUBROUTINE GET_ENE_VAL

  SUBROUTINE GET_ECOR_ADJUST(C,TEMP)
     USE PSF
     USE SHAKE
     USE CONSTA,ONLY:KBOLTZ
     USE STREAM

     REAL(CHM_REAL),INTENT(OUT) :: C
     REAL(CHM_REAL),INTENT(IN)  :: TEMP
     INTEGER                    :: NDEGF,I

     NDEGF=0 
     DO I=1,NATOM
        IF(IMOVE(I).EQ.0) NDEGF=NDEGF+3
     ENDDO
     IF(QSHAKE) THEN
       DO I=1,NCONST
          IF((IMOVE(SHKAPR(1,I)) == 0).OR.(IMOVE(SHKAPR(2,I)) == 0)) THEN
             NDEGF=NDEGF-1  
          ENDIF
       ENDDO
     ENDIF  
     IF(PRNLEV.GT.1) WRITE(OUTU,'(A,I4)') 'GET_ECOR_ADJUST> NUMBER OF DEGREES OF FREEDOM = ', NDEGF

     IF(TEMP.LE.0) &
        CALL WRNDIE(-3,'<GET_ECOR_ADJUST>','ENERGY CORRECTION REQUESTED WITH INVALID TEMPERATURE')
     C=(KBOLTZ*TEMP*NDEGF)/2.0
     IF(PRNLEV.GT.4) WRITE(OUTU,'(A,F10.4)') 'GET_ENE_VAL> ENERGY CORRECTION = ', C

  END SUBROUTINE GET_ECOR_ADJUST

#else /* repdstr_main */

  SUBROUTINE REPDSTRMAIN
    CALL WRNDIE(-1, &
         '<REPDSTR>','REPlica DiSTRibuted code not compiled.')
    RETURN
  END SUBROUTINE REPDSTRmain
#endif /* repdstr_main */


END MODULE REPDSTRMOD
