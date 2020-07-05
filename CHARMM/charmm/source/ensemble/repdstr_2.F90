module repdstrmod2
  use chm_kinds
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

  use parallel, only: maxnode
#if KEY_REPDSTR2==1
  use repdstr, only: maxrepdstr,qrepioset

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

  logical,save                                         :: qecor, qrepdverbose

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
  logical,save,public :: qphrex
#if KEY_CONSPH==1
  real(chm_real),save,public  :: phrx(maxnode)
  real(chm_real)              :: rhph,rlph      ! pH values of the reservoir
  integer,allocatable,target,dimension(:)   :: nprothigh,nprotlow ! store number of protonated groups in each reservoir entry
  integer,allocatable,target,dimension(:,:) :: reshtstate,resltstate
  real(chm_real4),allocatable,dimension(:)  :: rescg
#endif
  logical,save,public :: qresboltz,qresnobo ! what exchange criterion to use for reservoir
  integer,save,public :: irex = 0     ! number of exchanges
  integer,save,public :: isuc = 0     ! number of successful exchanges
  integer,save,public :: iresexch = 0 ! number of reservoir exchanges
  integer,save,public :: rhunit,rlunit ! units for high and low resorvoirs
  integer,save,public :: highressz,lowressz,maxressz ! current and maximum size(s) of the 
                                                     ! reservoir(s) in adaptive calculations
  integer,save,public :: repdid     ! replica id at current replica
  integer,save,public :: nrepeat    ! number of times to repeat exchange procedure
  integer,save        :: noppup,noppdn,nsucup,nsucdn ! keep track of successes and failures
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
#endif
  !
! RLH ->
  integer,allocatable,dimension(:) :: map_node2rep, map_rep2node, map_buffer
! <- RLH
  !------------------END of DECLARATIONS------------------------------


contains
  !*******************************************************************
  !                SUBROUTINES
  !*******************************************************************
  !

#if KEY_REPDSTR2==1 /* repdstr_main */
  !********************************************************
  !        REPDSTRMAIN
  !********************************************************
  subroutine repdstrmain
    !-----------------------------------------------------------------------
    !
#ifndef KEY_PARALLEL /* pll */
    CALL WRNDIE(-5,'<REPDSTR2>', &
         'REPlica DiSTRibuted runs only in parallel.')
#endif /* pll */
    !-----------------------------------------------------------------------
    use number
    use dimens_fcm
    use comand
    use cstuff, only: getpid
    use psf
    use coord
    use parallel
!    use repdstr
    use stream
    use string
    use memory
    use deriv
    use block_ltm
    use lambdam    !GG MSLD-compatibility
    use consta     !GG MSLD-compatibility
    use repd_ensemble
    use repdstr       !,only:qrepdnames,qrepioset,nrepdstr,natrepcmd
    use param_store,only:set_param

    implicit none
    
    logical qfin,qfinsync,qfinone,qsync,qsynce

    real(chm_real) stemp,dtemp,mtemp,sgtemp,dsgtemp,msgtemp,sgft,dsgft

    character(len=100) fname
    character(len=6) repnamestyle
    integer :: repnamestylelen,repnamestylemax=6
    integer j !gg new msld ph-rex commands added
    integer i,ioerr,mypid,eneun,probunit,nat3,rdim
    integer ierror

    !=====================================================================
    !   *******  RESET **********
    if(indxa(comlyn,comlen,'RESE') > 0)then
       call repd2_fin
       return
    end if

    !=====================================================================
    !   *******  SYNC **********
    !     This is global operation so every input script
    !     should have MATCHING repd sync!!!
    !
    qsync=(indxa(comlyn,comlen,'SYNC') > 0)
    if(qsync)then
       repdsynccount=repdsynccount+1
       call repd_reps_barrier(ierror)
       if(prnlev > 2) then
          write(outu,'(a,i5)')'REPDSTR>sync count=',repdsynccount
       endif

       return
    endif

    !=====================================================================
    !   *******  REPIO = NORMAL or REPD **********
    !   NORMAL: setting REPIO to normal will not add the _# at the end of file names
    !           Users must use their own method of keeping track of files and filenames
    !   REPD:   <default> Original REPD method of adding _# to each input or output file.
    !           where # is the replica number of this rep.
    qrepdnames=.true.
!    if (indxa(comlyn,comlen,'REPIO') > 0) then
       call gtrmwd(comlyn,comlen,'REPIO',5,repnamestyle,repnamestylemax,repnamestylelen)
       qrepdnames = .not. (repnamestyle == "NORMAL")
       if(qrepdnames)then
          if(mynod==0)write(outu,'(a)')"REPD using old repd file names with trailing _#"
       else
          if(mynod==0)write(outu,'(a)')"REPD using normal file names without trailing _#"
       endif
!    endif

    !=====================================================================
    !   *******  IORES, IOSET **********
    !   sometimes we want this in our scripts :-)
    if (indxa(comlyn,comlen,'IORES') > 0) then
       qrepioset = .false.
       return
    endif

    if (indxa(comlyn,comlen,'IOSET') > 0) then
       qrepioset = .true.
       return
    endif

    !===================================================================
    !   *******  NREP,NATR **********
    nrepdstr = gtrmi(comlyn,comlen,'NREP',1)
    natrepcmd = gtrmi(comlyn,comlen,'NATR',-1)
    qrepdverbose= (indxa(comlyn,comlen,'VERBOSE') > 0)

    qrepdstr = .true.
    !===================================================================
    !   *******  OUTU **********
    outu_repdstr = gtrmi(comlyn,comlen,'OUTU',2)

    !=========================================================================
    !   ***************  COMMUNICATOR SETUP  **************************
    !=========================================================================
    !---- Set up communicators and file I/O
    call repd_ensnensem(nrepdstr,qrepdverbose)

    !===================================================================
    !   *******  FAST **********
    qfastrepdstr = (indxa(comlyn,comlen,'FAST') > 0)
    if(qfastrepdstr) then
       call chmalloc('repdstr.src','repdstrmain','curtemps',nrepdstr,crl=curtemps)
       ffrun=.true.
       loglevel=gtrmi(comlyn,comlen,'LOGL',1)
    endif

    ! qrepioset = .true.   !MFC why is this here, it unsets what was done with iores/ioset
    !     This is the default situation. This changes when stream is called!
    qrdqtt = .false.
    !     This is the default situation. This changes when outu is called!
    qwrqtt = .false.

    !=========================================================================
    !   *******  REPE **********
    nrepeat=gtrmi(comlyn,comlen,'REPE',1)
    if(prnlev > 2) write(outu,'(a,i3)') 'REPDSTR> EXCHANGE REPEAT = ', NREPEAT


    !future    !=========================================================================
    !future    !   *******  TWOD **********
    !future    q2drex = (indxa(comlyn,comlen,'TWOD') > 0)
    !future    if(q2drex) then
    !future       call setup_2d(comlyn,comlen)
    !future       return
    !future    endif

#if KEY_CONSPH==1
    !=============================================================
    !   *******  pH Replica Echange (PHRE) **********
    ! Tim Miller: June, 2011
    ! Decide if we need to do PH based replica exchange
    qphrex = (indxa(comlyn,comlen,'PHRE') > 0)

    if(qphrex) then
       !    repd nrep <int> phrex freq <int> phval <real> phval <real> ...

       if(qfastrepdstr) call wrndie(-4,'<REPDSTR>','FAST PH REX IS NOT SUPPORTED')

       repseed=123+irepdstr
       iunrex=gtrmi(comlyn,comlen,'UNIT',outu)
       irexfq=gtrmi(comlyn,comlen,'FREQ',1)
      
       do i=1,nrepdstr
          phrx(i)=gtrmf(comlyn,comlen,'PHVA',-one)
       enddo
    endif
#endif

    !=========================================================================
    !   *******  EXCH, EXLM **********
    qrexchg  = (indxa(comlyn,comlen,'EXCH') > 0)
    if(qrexchg .and. .not. qfastrepdstr) &
         call wrndie(-4,"<repdstr_2.src>repdstrmain", &
         "repd exchange not supported, try fast version")
    qrexchgl = (indxa(comlyn,comlen,'EXLM') > 0)

    qreservoir=(indxa(comlyn,comlen,'RSRV') > 0)
    qreshigh  =(indxa(comlyn,comlen,'RESH') > 0)
    qreslow   =(indxa(comlyn,comlen,'RESL') > 0)
    call repd_reservoir_setup   !-- contained subroutine at end of this routine

    !     This is initialization. Use PSETLOC(), PSETGLOB() elsewhere
    !     This needs to be generalized!!! Also first two lines are now in paral1.src!!!


    
    !=========================================================================
    !   *******  Startup for REPD/ENSEMBLE  **********
    !=========================================================================
    ! 
    numnod=numnodg/nrepdstr
    if(prnlev > 2)write(outu,'(a,2i5)') &
         ' REPD Startup> Number of processes and NREP=',numnodg,nrepdstr
    call repd_check_numnod
    irepdstr=mynodg/numnod
    repdid=irepdstr

    !=========================================================================
    !        Set scripting parameters myrep and nrep
    call set_param('MYREP',irepdstr)
    call set_param('NREP',nrepdstr)
    call drepsetio(iolev,prnlev,wrnlev)

    !=========================================================================
    !     Replica exchange data: get the temperatures...
    qexpt = (indxa(comlyn,comlen,'expt') > 0)
    qex2d = (indxa(comlyn,comlen,'ex2d') > 0)
    qexbk = (indxa(comlyn,comlen,'exbk') > 0)

    call repd_exchange_setup  !-- contained subroutine at end of this routine

    if(qrexchgl)then
       ! call wrndie(-4,'<REPDSTR>','FAST LAMBDA REX IS NOT SUPPORTED')

       qsump=(indxa(comlyn,comlen,'SUMP') > 0)

       repseed=123+mynodg
       iunrex=gtrmi(comlyn,comlen,'UNIT',outu)
       irexfq=gtrmi(comlyn,comlen,'FREQ',1)
       !!IF(PRNLEV > 2) WRITE(OUTU,'(A,I5)') 'TIM DEBUG> IREXFQ = ',IREXFQ

       if((.not.qexpt).and.(.not.qex2d)) then
          !     repd nrep <int> exlm freq <int>
#if KEY_PHMD==1
          qphrx=(indxa(comlyn,comlen,'PHMD') > 0) !JAW. Flag for using PHMD and PH exchange
#endif
#if KEY_BLOCK==1
          qmsphrx=(indxa(comlyn,comlen,'MSPH') > 0) !GG: Flag for using pH-REX in CPHMD^MSLD
          if(qmsphrx)then
             !gg: usage example "repd nrep <int> exlm freq <int> msph sph <int> mph <int>"
             if (nrepdstr  ==  nreplica) then        !gg: check number of replicas are the same in block and repdstr
                call msld_phrex(comlyn,comlen,nrepdstr,1,nblock)
             else
                call wrndie(-5,'<REPDSTR>', &
                'Number of replicas declared in BLOCK and REPDSTR do not match!')
             endif        !gg: for "nrepdstr  ==  nreplica" loop
          endif           !gg: for "qmsphrx" loop
#endif

          if(qexbk)irbk=gtrmi(comlyn,comlen,'REBK',1)
       else if (qexpt) then
          !     repd nrep <int> exlm expt nrpt <int> freq <int>
          nrept=gtrmi(comlyn,comlen,'NRPT',1)
          if(qexbk)irbk=gtrmi(comlyn,comlen,'REBK',1)
       ELSE
          ! repd nrep <int> exlm ex2d nrpx <int> freq <int>
          nrepx=gtrmi(comlyn,comlen,'NRPX',1)
          if(qexbk)irbk=gtrmi(comlyn,comlen,'REBK',1)
       endif
    endif

    ! initialize the starting temperature of this replica
    if(qfastrepdstr) then
       tempcurrent=curtemps(irepdstr+1)
       if(prnlev >= 6) &
            write(outu,'(a,i3,a,f10.3)') 'FREX DEBUG> REPL ', IREPDSTR, ' INIT TEMPCURRENT = ', TEMPCURRENT
    endif
    reptag=irepdstr
    noppup=0
    noppdn=0
    nsucup=0
    nsucdn=0

    call flush(outu)
    mypid = getpid()
    call repd_reps_barrier(ierror)
    return

  contains
    subroutine repd_check_numnod
      if(numnod*nrepdstr == numnodg) return
      write(outu,*)" *** ERROR *** nrep must be a factor of number of mpi processes"
      call wrndie(-5,'<REPDSTR>', &
           'Wrong combination of NREP and number of processes.')
    end subroutine repd_check_numnod

    subroutine repd_exchange_setup  
      if(qrexchg) then
         !     repd nrep <int> exch freq <int> temp <real> temp <real> ...
         !
         !     could be also:
         !     repd nrep <int> exch stemp <real> dtemp <real>
         !     where stemp is starting temperature and dtemp is the
         !     interval between the temperatures.
         repseed=123+irepdstr
         qsump   =(indxa(comlyn,comlen,'SUMP') > 0)
         qrxsgld =(indxa(comlyn,comlen,'SGLD') > 0)
         qrxtiger=(indxa(comlyn,comlen,'TIGE') > 0)
         qrxtham =(indxa(comlyn,comlen,'THAM') > 0)
#if KEY_PHMD==1
         qphrx=(indxa(comlyn,comlen,'PHMD') > 0) !jaw. flag for using phmd and ph exchange
#endif
         iunrex=gtrmi(comlyn,comlen,'UNIT',outu)
         irexfq=gtrmi(comlyn,comlen,'FREQ',1)
         stemp=gtrmf(comlyn,comlen,'STEM',-one)
         if(qrxtham.and.qrxsgld) &
              call wrndie(-4,'<REPDSTR>', &
              'TEMP-HAMILTONIAN REPLICA EXCHANGE IS INCOMPATIBLE WITH SGLD REM')
         if(stemp > zero) then
            mtemp=gtrmf(comlyn,comlen,'MTEM',-one)
            if(mtemp>zero)then
               dtemp=exp(dlog(mtemp/stemp)/(nrepdstr-1))
               do i=1,nrepdstr
                  temprx(i)=stemp*dtemp**(i-1)
               enddo
            else
               dtemp=gtrmf(comlyn,comlen,'DTEM',-one)
               if(dtemp < zero) call wrndie(-5,'<REPDSTR>', &
                    'replica EXCHange needs interval between temperatures.')
               do i=1,nrepdstr
                  temprx(i)=stemp+(i-1)*dtemp
               enddo
            endif
         else
            do i=1,nrepdstr
               temprx(i)=gtrmf(comlyn,comlen,'TEMP',-ONE)
               if(temprx(i) < zero) then
                  write(outu,*)'Number of specifed temperatures ',i, &
                       ', but ',nrepdstr,' are needed'
                  CALL WRNDIE(-5,'<REPDSTR>', &
                       'replica EXCHange needs all temperatures.')
               endif
               if(qfastrepdstr) curtemps(i)=temprx(i)
            enddo
         endif
         if(qrxsgld)then
            if(qfastrepdstr) &
                 CALL WRNDIE(-4,'<REPDSTR>','FAST SGLD REX IS NOT SUPPORTED')
            
            allocate(sgtemprx(nrepdstr))
            allocate(sgftrx(nrepdstr))
            sgtemp=gtrmf(comlyn,comlen,'SGTE',zero)
            sgft=gtrmf(comlyn,comlen,'SGFT',zero)
            dsgft=gtrmf(comlyn,comlen,'DSGF',zero)
            do i=1,nrepdstr
               sgftrx(i)=sgft+(i-1)*dsgft
            enddo
            msgtemp=gtrmf(comlyn,comlen,'MSGT',-one)
            if(sgtemp>zero)then
               if(msgtemp>zero)then
                  dsgtemp=exp(dlog(msgtemp/sgtemp)/(nrepdstr-1))
                  do i=1,nrepdstr
                     sgtemprx(i)=sgtemp*dsgtemp**(i-1)
                  enddo
               else 
                  dsgtemp=gtrmf(comlyn,comlen,'DSGT',zero)
                  do i=1,nrepdstr
                     sgtemprx(i)=sgtemp+(i-1)*dsgtemp
                  enddo
               endif
            else
               do i=1,nrepdstr
                  sgtemprx(i)=gtrmf(comlyn,comlen,'SGTT',-one)
                  if(sgtemprx(i) < zero) then
                     write(outu,*)'No SGTT input on stage ',i, &
                          ', set to simulation temperature: ',temprx(i)
                     sgtemprx(i)=temprx(i)
                  endif
               enddo
            endif
         endif
         !
         !     TIGER parameters parsed here:
         !
         if(qrxtiger)then
            tigergr=gtrmf(comlyn,comlen,'TOLG',zero)
            tigerit=gtrmi(comlyn,comlen,'ITER',1)
            tigernm=gtrmi(comlyn,comlen,'NMIN',100)
            tigerneq=gtrmi(comlyn,comlen,'NEQU',1000)
         endif
      endif
! RLH ->
      call chmalloc('repdstr.src','repd_exchange_setup','map_node2rep',nrepdstr,intg=map_node2rep)
      call chmalloc('repdstr.src','repd_exchange_setup','map_rep2node',nrepdstr,intg=map_rep2node)
      call chmalloc('repdstr.src','repd_exchange_setup','map_buffer',nrepdstr,intg=map_buffer)
      do i=1,nrepdstr
        map_rep2node(i)=i-1
        map_node2rep(i)=i-1
      enddo
! RLHFIX - never free mapping arrays
! <- RLH
      return
    end subroutine repd_exchange_setup


    subroutine repd_reservoir_setup
      reservoir_setup: if(qreservoir) then

         if(qfastrepdstr) &
              call wrndie(-4,'<REPDSTR>','FAST RESERVOIR REX IS NOT SUPPORTED')

         qresboltz=(indxa(comlyn,comlen,'BOLT') > 0)
         qresnobo =(indxa(comlyn,comlen,'NOBO') > 0)
         if(.not. (qreshigh .or. qreslow)) &
              call wrndie(-3,'<REPDSTRMAIN>', 'RESERVOIR NEEDS RESHIGH OR RESLOW')

         if(qrexchgl) then
            ! Only allow "Boltzmann" for Hamiltonian REX
            qresnobo=.false.
            qresboltz=.true.
            if(prnlev >= 3) write(outu,'(a)') &
                 'REPDSTRMAIN> HAMILTONIAN REX IN USE. BOLTZMANN EXCHANGE CRITERION ACTIVATED.'
         else
            ! check to make sure there aren't multiple specs
            if(qresboltz .and. qresnobo) then
               call wrndie(-4,'<REPDSTRMAIN>', 'CONFLICTING EXCHANGE CRITERIA SET FOR RESERVOIR REX')
            endif
            if(.not. ( qresboltz .or. qresnobo )) then
               call wrndie(0,'<REPDSTRMAIN>', &
                    'NO EXCHANGE CRITERIA SET: RESERVOIR REX WILL USE BOLTZMANN CRITERION')
               qresboltz=.true.
            endif
         endif
         if(qresboltz.or.qresnobo) then
            qecor=(indxa(comlyn,comlen,'ECOR') > 0)
         else
            qecor=.false.
         endif

         call chmalloc('repdstr.src','repdstrmain','rescrdx',natom,cr4=rescrdx)
         call chmalloc('repdstr.src','repdstrmain','rescrdy',natom,cr4=rescrdy)
         call chmalloc('repdstr.src','repdstrmain','rescrdz',natom,cr4=rescrdz)
#if KEY_CONSPH==1
         call chmalloc('repdstr.src','repdstrmain','rescg',natom,cr4=rescg)
#endif

         ! Get info on the trajectory files containing the reservoirs      
         if(qreshigh) then
            highressz=gtrmi(comlyn,comlen,'RHSZ',-1)
            rhunit=gtrmi(comlyn,comlen,'RHUN',-1)
            if(rhunit < 0) &
                 call wrndie(-3,'<REPDSTRMAIN>','BAD UNIT FOR TOP RESERVOIR')

            if(highressz < 0) call wrndie(-3,'<REPDSTRMAIN>','MUST SPECIFY A VALID HIGH RESERVOIR SIZE')


            if(qrexchgl) then
               rhtemp=gtrmf(comlyn,comlen,'RHTE',-1.0)
               if(rhtemp <= 0) &
                    call wrndie(-3,'<REPDSTRMAIN>','BAD BOLTZMANN TEMP FOR HIGH RESERVOIR')
               eneun=gtrmi(comlyn,comlen,'FHEN',-1)
               if(eneun < 1) &
                    call wrndie(-3,'<REPDSTRMAIN>','RESERVOIR H-REX CANNOT PRECALC ENERGIES.')
               call get_ene_val(eneun,qecor,rltemp,.true.)
            else
               rhtemp=gtrmf(comlyn,comlen,'RHTE',-1.0)
               if(qresboltz) then
#if KEY_CONSPH==1
                  if(qphrex) then
                     ! handle boltzmann ph rex, hoorah
                     rhph=gtrmf(comlyn,comlen,'RHPH',-1.0)
                  else
#endif
                     ! Get temp of high reservoir
                     if(rhtemp <= 0) &
                          call wrndie(-3,'<REPDSTRMAIN>','BAD BOLTZMANN TEMP FOR HIGH RESERVOIR')
#if KEY_CONSPH==1
                  endif
#endif
               endif

               if(qresboltz.or.qresnobo) then
                  eneun=gtrmi(comlyn,comlen,'FHEN',-1)
                  if(eneun > 0) then
                     call get_ene_val(eneun,qecor,rhtemp,.true.)
                  else
                     call precalcene(.true.,x,y,z,qecor,rhtemp)
                  endif
               endif
            endif
         endif
         if(qreslow) then
            lowressz=gtrmi(comlyn,comlen,'RLSZ',-1)
            rlunit=gtrmi(comlyn,comlen,'RLUN',-1)          
            if(rlunit < 0) &
                 call wrndie(-3,'<REPDSTRMAIN>','BAD UNIT FOR BOTTOM RESERVOIR')

            if(lowressz < 0) call wrndie(-3,'<REPDSTRMAIN>','MUST SPECIFY A VALID LOW RESERVOIR SIZE')


            if(qrexchgl) then
               rltemp=gtrmf(comlyn,comlen,'RLTE',-1.0)
               if(rltemp <= 0) &
                    call wrndie(-3,'<REPDSTRMAIN>','BAD BOLTZMANN TEMP FOR LOW RESERVOIR')
               eneun=gtrmi(comlyn,comlen,'FLEN',-1)
               if(eneun < 1) &
                    call wrndie(-3,'<REPDSTRMAIN>','RESERVOIR H-REX CANNOT PRECALC ENERGIES.')
               call get_ene_val(eneun,qecor,rltemp,.false.)
            else
               rltemp=gtrmf(comlyn,comlen,'RLTE',-1.0)
               if(qresboltz) then
                  ! Get temp of high reservoir
#if KEY_CONSPH==1
                  if(qphrex) then
                     ! handle boltzmann ph rex, hoorah
                     rlph=gtrmf(comlyn,comlen,'RLPH',-1.0)
                  else
#endif
                     if(rltemp <= 0) &
                          call wrndie(-3,'<REPDSTRMAIN>','BAD BOLTZMANN TEMP FOR LOW RESERVOIR')
#if KEY_CONSPH==1
                  endif
#endif
               endif


#if KEY_CONSPH==1
               if(qphrex) then
                  if(qfastrepdstr) &
                       call wrndie(-4,'<REPDSTR>','FAST PH REX IS NOT SUPPORTED')
               else
#endif
                  if(qresboltz.or.qresnobo) then
                     eneun=gtrmi(comlyn,comlen,'FLEN',-1)
                     if(eneun > 0) then
                        call get_ene_val(eneun,qecor,rltemp,.false.)
                     else
                        call precalcene(.false.,x,y,z,qecor,rltemp)
                     endif
                  endif
#if KEY_CONSPH==1
               endif
#endif
            endif
         endif
      else reservoir_setup
         qresboltz=.false.
         qresnobo=.false.
      endif reservoir_setup
      return
    end subroutine repd_reservoir_setup


    subroutine repd2_fin
      qfin=.true.
      !        Must be called before qrepdstr=.false.
      !        Restore the I/O; can be called individually! No global comm.!
      call drepresio(iolev,prnlev,wrnlev)
      qrepdstr = .false.
      qrdqtt = .false.
      qwrqtt = .false.
      ! is iparpt OK?
      ! is ippmap OK?
      ! is INODE() OK?
      qfinsync=(indxa(comlyn,comlen,'SYNC') > 0)
      if(qfinsync)call psync()
      qfinone=(indxa(comlyn,comlen,'PONE') > 0)
      if(qfinone)then
         mynod=0
         numnod=1
         call cube(mynod,numnod,ippmap)
      endif
      if(qrxsgld)deallocate(sgtemprx)
      if(qrxsgld)deallocate(sgftrx)
      
      if(qfastrepdstr) &
           call chmdealloc('repdstr.src','repdstrmain','curtemps',nrepdstr,crl=curtemps)
      if(qreservoir) then
         call chmdealloc('repdstr.src','repdstrmain','rescrdx',natom,cr4=rescrdx)
         call chmdealloc('repdstr.src','repdstrmain','rescrdy',natom,cr4=rescrdy)
         call chmdealloc('repdstr.src','repdstrmain','rescrdz',natom,cr4=rescrdz)
#if KEY_CONSPH==1
         call chmdealloc('repdstr.src','repdstrmain','rescg',natom,cr4=rescg)
#endif
         
      endif
      call repd2_ensfin
!      call repd2_clear_communicators
      return
    end subroutine repd2_fin
  end subroutine repdstrmain




  !*******************************************************************
  !                FASTREPEXCHG
  !*******************************************************************
  subroutine fastrepexchg(wmain,epot,istep,jhstrt,igvopt,vx,vy,vz, &
       xold,yold,zold &
#if KEY_TSM==1
       ,backls &
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
     use repd_ensemble

     ! Arguments
     real(chm_real),intent(in)    :: wmain(*),epot
     real(chm_real),intent(inout) :: vx(*),vy(*),vz(*),xold(*),yold(*),zold(*)
     integer,intent(in)           :: istep
     integer,intent(out)          :: jhstrt
     integer,intent(inout)        :: igvopt
#if KEY_TSM==1
     integer backls(*)
#endif

     ! Local variables
     integer                      :: i,j,x,ierr,rep1,rep2,tmpun
     integer                      :: ourselves,neighbor,reporder(maxrepdstr)
     logical                      :: qexc,qmaster,qdolog
     real(chm_real)               :: epotarr(maxrepdstr),nperrep,ourtemp
     real(chm_real)               :: nbrtemp,p,prob,lowtemp,lastlow
     real(chm_real),allocatable,dimension(:) :: scarr

     if(mod(istep-1,irexfq) /= 0) return
    
     if(prnlev >= 6) write(outu,'(a,i6)') 'FREX DEBUG> IN FASTREPEXCHG AT STEP ', istep
     if(nrepdstr > maxrepdstr) &
        call wrndie(-5,'<FASTREPEXCHG>','TOO MANY REPLICAS FOR FAST EXCHANGING')
     do i=1,maxrepdstr
        qexc=.false.
     enddo

     nperrep=numnodg/nrepdstr
     if(mynodg == 0) &
        call chmalloc('repdstr.src','fastrepexchg','scarr',nrepdstr,crl=scarr)
     if(lmasternode)then
        call mpi_gather(epot,1,mpi_real8,scarr,1,mpi_real8,0,comm_master,ierr)
        if(ierr /= mpi_success) &
             call wrndie(-5,'<fastrepexchg>','bungled mpi communication')
     endif

     !=================================================================
     !   Master of all masters does this work
     !=================================================================
     masterwork: if(mynodg == 0) then
        epotarr(1:nrepdstr)=scarr(1:nrepdstr)
        call chmdealloc('repdstr.src','fastrepexchg','scarr',nrepdstr,crl=scarr)
        if(ffrun) then
           write(iunrex,'(a)') '# replica temp. ener. neighbor ntemp nene prob p success? newrep'
           ffrun=.false.
        endif

        ! We need to actually figure out how the exchanges are going to
        ! happen, to do so, apply the formula and temps to the energies.
        repeatloop: do x=1,nrepeat
           if(x == 1.or.x == nrepeat) then
              qdolog=.true.
           else if(loglevel == 0) then
              qdolog=.false.
           else if(mod(x,loglevel) == 0) then
              qdolog=.true.
           else
              qdolog=.false.
           endif

           irex=irex+1
           if(qdolog) write(iunrex,'(a,i15,a,i12,a,i5)') '# Exchange ', IREX, ': STEP ', ISTEP-1, ': REPEAT ', X

           ! put the reservoirs in order of temperature ... a bit of a hacky
           ! bubble sort for now.
           lastlow=9999.0
           do i=1,nrepdstr
              lowtemp=9999.0
              do j=1,nrepdstr
                 if(i > 1) then
                    if(curtemps(j) > lastlow) then
                       if(curtemps(j) < lowtemp) then
                          reporder(i)=j
                          lowtemp=curtemps(j)
                       endif
                    endif
                 else
                    if(curtemps(j) < lowtemp) then
                       reporder(i)=j
                       lowtemp=curtemps(j)
                    endif
                 endif
              enddo
              lastlow=lowtemp
           enddo

           qmaster=(mod(irex,2) == 0)
           mloop2: do i=1,nrepdstr
              qexc=.false.
              ourselves=reporder(i)
              ourtemp=curtemps(ourselves)
              
              ! find our neighbor (next highest temperature replica)
              if(qmaster) then
                 if(i <= nrepdstr-1) then
                    neighbor=reporder(i+1)
                    nbrtemp=curtemps(neighbor)
                 else 
                    neighbor=-one
                    nbrtemp=9999.0
                 endif
              else
                 ! special case to make sure that there's a print out for the first
                 ! replica.
                 if(i == 1) then
                    if(qdolog) write(iunrex,'(i2,x,f12.6,x,f15.6,x,i2,x,f12.6,x,f15.6,x,f6.3,x,f6.3,x,l,x,i2)') &
                               ourselves,ourtemp,epotarr(ourselves),-1,9999.0,0.00,0.00,0.00,.false.,ourselves
                 endif
                 qmaster=.not.qmaster
                 cycle
              endif

              if(neighbor > 0) then
                 ! we control the exchange; our neighbor has the next
                 ! highest temperature.
                 if(ourtemp == nbrtemp) then
                    prob=one
                 else
                    prob=min(one,exp(-(one/(kboltz*ourtemp) &
                             -one/(kboltz*nbrtemp))*(epotarr(neighbor)-epotarr(ourselves))))
                 endif
                 p=random(repseed)
                
                 if(p <= prob) then
                    qexc=.true.
                    curtemps(neighbor)=ourtemp
                    curtemps(ourselves)=nbrtemp
                    rep1=neighbor
                    rep2=ourselves
                 else
                    qexc=.false.
                    rep1=ourselves
                    rep2=neighbor
                 endif
                 if(qdolog) write(iunrex,'(i2,x,f12.6,x,f15.6,x,i2,x,f12.6,x,f15.6,x,f6.3,x,f6.3,x,l,x,i2)') &
                            ourselves,ourtemp,epotarr(ourselves),neighbor,nbrtemp,epotarr(neighbor),prob,p,qexc,rep1

                 ! write out a line for the neighboring replica, as well
                 if(qdolog) write(iunrex,'(i2,x,f12.6,x,f15.6,x,i2,x,f12.6,x,f15.6,x,f6.3,x,f6.3,x,l,x,i2)') &
                            neighbor,nbrtemp,epotarr(neighbor),ourselves,ourtemp,epotarr(ourselves),prob,p,qexc,rep2

              else
                 if(qdolog) write(iunrex,'(i2,x,f12.6,x,f15.6,x,i2,x,f12.6,x,f15.6,x,f6.3,x,f6.3,x,l,x,i2)') &
                      ourselves,ourtemp,epotarr(ourselves), &
                      -1,9999.0,0.00,0.00,0.00,.false.,ourselves

              endif ! end if(ourselves > 0)
              qmaster=.not.qmaster
           enddo mloop2
        enddo repeatloop

     endif masterwork ! end part only executed on processor 0

     ! Now that we have figured out the final temperatures at each state,
     ! broadcast curtemps and call dofastexchg to actually make the exchange.
     ! only masternodes will do this communication.
     if(lmasternode) &
          call mpi_bcast(curtemps,nrepdstr,mpi_real8,0,comm_master,ierr)

     call dofastexchg(curtemps(irepdstr+1),jhstrt,igvopt,vx,vy,vz, &
          xold,yold,zold)
     call flush(outu)
  end subroutine fastrepexchg

  !********************************************************
  !        DOFASTEXCHG
  !********************************************************
  subroutine dofastexchg(tempnew,jhstrt,igvopt,vx,vy,vz,xold,yold,zold)
     use mpi
     use psf
     use parallel
     use stream
     use repdstr
     use coord,only: x,y,z
     use repd_ensemble
     ! Arguments
     real(chm_real),intent(in)    :: tempnew
     real(chm_real),intent(inout),dimension(natom) :: vx,vy,vz,xold,yold,zold
     integer,intent(out)          :: jhstrt
     integer,intent(inout)        :: igvopt

     ! Local variables
     logical        :: qupveloc
     real(chm_real) :: scalef
     integer        :: i,ierr,comm_save,mynod_save

     if(lmasternode) then
        write(outu,'(a,2f10.3)') 'FREX DEBUG> FORMER AND NEW TEMPS: ', &
             tempcurrent, tempnew
        if(tempcurrent /= tempnew) then
           qupveloc=.true.
           scalef=sqrt(tempnew/tempcurrent)
           if(prnlev >= 6) write(outu,'(a,f7.4)')'FREX DEBUG> SCALEF: ', SCALEF
           vx  (1:natom)=vx  (1:natom)*scalef
           vy  (1:natom)=vy  (1:natom)*scalef
           vz  (1:natom)=vz  (1:natom)*scalef
           xold(1:natom)=xold(1:natom)*scalef
           yold(1:natom)=yold(1:natom)*scalef
           zold(1:natom)=zold(1:natom)*scalef
           tempcurrent=tempnew
        else
           qupveloc=.false.
        endif
     endif
 
     call mpi_bcast(qupveloc,1,mpi_logical,0,comm_charmm,ierr)
     call mpi_bcast(tempcurrent,1,mpi_double_precision,0,comm_charmm,ierr)
     if(qupveloc) then
        jhstrt=0
        igvopt=2
        
        ! BTM: I am not sure if we still need to send xold, yold, and zold,
        ! but I am going to do it anyway to be safe.
     call mpi_bcast(igvopt,1,mpi_int,0,comm_charmm,ierr)
     call mpi_bcast(jhstrt,1,mpi_int,0,comm_charmm,ierr)
     call mpi_bcast(vx,natom,mpi_double_precision,0,comm_charmm,ierr)
     call mpi_bcast(vy,natom,mpi_double_precision,0,comm_charmm,ierr)
     call mpi_bcast(vz,natom,mpi_double_precision,0,comm_charmm,ierr)
     call mpi_bcast(xold,natom,mpi_double_precision,0,comm_charmm,ierr)
     call mpi_bcast(yold,natom,mpi_double_precision,0,comm_charmm,ierr)
     call mpi_bcast(zold,natom,mpi_double_precision,0,comm_charmm,ierr)
     endif               
  end subroutine dofastexchg

!future  SUBROUTINE SETUP_2D(COMLYN,COMLEN)
!future     use repdstr
!future     use parallel
!future     use string
!future     use stream
!future
!future     ! Passed variables
!future     character(len=*)      :: comlyn
!future     integer               :: comlen
!future
!future     ! Local variables
!future     integer               :: i,j
!future     character(len=4)      :: scratch
!future
!future     ndim1=gtrmi(comlyn,comlen,'DIM1',0)
!future     ndim2=gtrmi(comlyn,comlen,'DIM2',0)
!future     d1freq=gtrmi(comlyn,comlen,'D1FR',-1)
!future     d2freq=gtrmi(comlyn,comlen,'D2FR',-1)
!future     iunrex=gtrmi(comlyn,comlen,'UNIT',-1)
!future
!future     irexd1=0
!future     irexd2=0
!future     irexfq=min(d1freq,d2freq)
!future     qrepdstr=.true.
!future     qrexchg=.true.
!future
!future     if(ndim1 <= 0 .or. ndim2 <= 0) call wrndie(-5,'<SETUP_2D>','BAD DIMENSIONS, DUMBASS!')
!future     if(d1freq < 1 .or. d2freq < 1) call wrndie(-5,'<SETUP_2D>','2D REX REQUIRES POSITIVE EXCHANGE FREQUENCIES')
!future
!future     scratch=gtrma(comlyn,comlen,'D1CR')
!future     if(scratch == 'TEMP') then
!future        q2dtemp(1)=.true.
!future        q2dham(1)=.false.
!future        q2dph(1)=.false.
!future        q2ditemp=(indxa(comlyn,comlen,'ITEM') > 0)
!future     else if(scratch == 'HAM') then
!future        q2dtemp(1)=.false.
!future        q2dham(1)=.true.
!future        q2dph(1)=.false.
!future#if KEY_CONSPH==1
!future     else if(scratch == 'PH') then
!future        q2dtemp(1)=.false.
!future        q2dham(1)=.false.
!future        q2dph(1)=.true.
!future#endif
!future     else
!future        call wrndie(-5,'<SETUP_2D>','UNRECOGNIZED EXCHANGE CRITERIA FOR DIMENSION 1.')
!future     endif
!future
!future     scratch=gtrma(comlyn,comlen,'D2CR')
!future     if(scratch == 'TEMP') then
!future        q2dtemp(2)=.true.
!future        q2dham(2)=.false.
!future        q2dph(2)=.false.
!future        q2ditemp=(indxa(comlyn,comlen,'ITEM') > 0)
!future     else if(scratch == 'HAM') then
!future        q2dtemp(2)=.false.
!future        q2dham(2)=.true.
!future        q2dph(2)=.false.
!future#if KEY_CONSPH==1
!future     else if(scratch == 'PH') then
!future        q2dtemp(2)=.false.
!future        q2dham(2)=.false.
!future        q2dph(2)=.true.
!future#endif
!future     else
!future        call wrndie(-5,'<SETUP_2D>','UNRECOGNIZED EXCHANGE CRITERIA FOR DIMENSION 2.')
!future     endif
!future
!future     if(q2dtemp(1).and.q2dtemp(2)) call wrndie(-5,'<SETUP_2D>','ONLY ONE TEMP. DIMENSION ALLOWED!')
!future     if(q2dph(1).and.q2dph(2)) call wrndie(-5,'<SETUP_2D>','ONLY ONE PH DIMENSION ALLOWED!')
!future
!future     do i=1,ndim1
!future        if(q2dtemp(1)) then
!future           temprx(i)=gtrmf(comlyn,comlen,'TEMP',-1.0)
!future
!future#if KEY_CONSPH==1
!future        else if(q2dph(1)) then
!future           phrx(i)=gtrmf(comlyn,comlen,'PHVA',-1.0)
!future#endif
!future        endif
!future     enddo
!future     do i=1,ndim2
!future        if(q2dtemp(2)) then
!future           temprx(i)=gtrmf(comlyn,comlen,'TEMP',-1.0)
!future
!future#if KEY_CONSPH==1
!future        else if(q2dph(2)) then
!future           phrx(i)=gtrmf(comlyn,comlen,'PH',-1.0)
!future#endif
!future        endif
!future     enddo
!future
!future     ! redistribute parallel set-up
!future     nrepdstr=ndim1*ndim2
!future     numnod=numnodg/nrepdstr
!future
!future     if(prnlev > 2) &
!future        write(outu,'(a,3i4)') 'SETUP_2D> Number of processors, NREPDSTR, and NUMNOD = ',numnodg,nrepdstr,numnod
!future     if(numnod*nrepdstr /= numnodg) call wrndie(-5,'<SETUP_2D>','WRONG COMBINATION OF REPLIUCAS AND PROCESSORS')
!future
!future     call psetloc
!future     call drepsetio(iolev,prnlev,wrnlev)
!future
!future     irepdstr=mynodg/numnod
!future     repdid=irepdstr
!future     repseed=irepdstr+123
!future
!future     myrepd1=mod(irepdstr,ndim1)
!future     myrepd2=irepdstr/ndim1 
!future
!future     if(myrepd1 < ndim1-1) then
!future        nbrup1=irepdstr+1
!future        cnbrup1=mynodg+numnod
!future     else
!future        nbrup1=-1
!future        cnbrup1=-1
!future     endif
!future     if(myrepd1 > 0) then
!future        nbrdn1=irepdstr-1
!future        cnbrdn1=mynodg-numnod
!future     else
!future        nbrdn1=-1
!future        cnbrdn1=-1
!future     endif
!future
!future     if(myrepd2 < ndim2-1) then
!future        nbrup2=irepdstr+ndim1
!future        cnbrup2=mynodg+(ndim1*numnod)
!future     else
!future        nbrup2=-1
!future        cnbrup2=-1
!future     endif
!future     if(myrepd2 > 0) then
!future        nbrdn2=irepdstr-ndim1
!future        cnbrdn2=mynodg-(ndim1*numnod)
!future     else 
!future        nbrdn2=-1
!future        cnbrdn2=-1
!future     endif
!future
!future     call set_param('MYREP',irepdstr)
!future     call set_param('MYREPD1',myrepd1)
!future     call set_param('MYREPD2',myrepd2)
!future     call set_param('NREP',nrepdstr)
!future     call set_param('NREPD1',ndim1)
!future     call set_param('NREPD2',ndim2)
!future
!future     if(mynod == 0) then
!future        write(outu,'(a,8i5)') 'SETUP_2D> MYNODG,MYNOD,IREPDSTR,MYREPD1,NBRUP1,NBRDN1,CNBRUP1,CNBRDN1 = ', &
!future                              mynodg,mynod,irepdstr,myrepd1,nbrup1,nbrdn1,cnbrup1,cnbrdn1
!future        write(outu,'(a,8i5)') 'SETUP_2D> MYNODG,MYNOD,IREPDSTR,MYREPD2,NBRUP2,NBRDN2,CNBRUP2,CNBRDN2 = ', &
!future                              mynodg,mynod,irepdstr,myrepd2,nbrup2,nbrdn2,cnbrup2,cnbrdn2
!future     endif
!future 
!future  END SUBROUTINE SETUP_2D
!future
!future  !
!future  SUBROUTINE DO2DEXCH(X,Y,Z,WMAIN,VX,VY,VZ,XOLD,YOLD,ZOLD,MYEPOT,TTEMP,ISTART,JHSTRT, &
!future                     ISEED,IASVEL,IGVOPT,CALLSEQ &
!future#if KEY_TSM==1
!future                     ,BACKLS & 
!future#endif
!future#if KEY_CONSPH==1
!future                     ,IDIDPHREX & 
!future#endif
!future                    )
!future    
!future    use psf
!future    use number
!future    use stream
!future    use parallel
!future    use repdstr
!future    use consta
!future    use memory
!future    use image
!future    use bases_fcm
!future    use energym, only: energy,eprop,epot
!future    use deriv,only: dx,dy,dz
!future    use consph,only: tstate !##CONSPH
!future    use clcg_mod,only: random
!future    use imgup,only: upimag
!future
!future    !
!future    ! passed-in variables
!future    real(chm_real) :: X(:), Y(:), Z(:),MYEPOT,VX(*),VY(*),VZ(*),WMAIN(*),TTEMP
!future    real(chm_real) XOLD(*),YOLD(*),ZOLD(*)
!future    INTEGER ISTART,JHSTRT,ISEED,IASVEL,IGVOPT,CALLSEQ
!future    INTEGER BACKLS(*) !##TSM
!future
!future    !
!future    ! Local variables
!future    real(chm_real),allocatable,dimension(:)     :: w
!future    real(chm_real),allocatable,dimension(:,:,:) :: transf
!future    real(chm_real),dimension(6)                 :: oldxtlabc
!future    real(chm_real)                              :: eneigh,p,rn,nbrtemp,scalef,ourtemp
!future    real(chm_real)                              :: epot1p,epot1q,epot2p,epot2q
!future    integer                                     :: exchdim,me,neighbor,cneighbor,step,i
!future    logical                                     :: qexc,qcrys
!future
!future#if KEY_CONSPH==1
!future    integer                           :: iproto, jproto, ididphrex, phresstruct
!future    integer, allocatable,dimension(:) :: nstate
!future    real(chm_real)                    :: ph_l, ph_m, ph_delta
!future
!future    ididphrex = 0
!future#endif
!future
!future    qcrys = (xtltyp /= '    ')
!future
!future    if(mod(istart-1,d2freq) == 0) then
!future       irexd2=irexd2+1
!future
!future       exchdim=2
!future       step=mod(irexd2,2)
!future       me=myrepd2
!future       if(step==1) then
!future          if(mod(me,2) == 0) then
!future             neighbor=nbrup2
!future             cneighbor=cnbrup2
!future          else
!future             neighbor=nbrdn2
!future             cneighbor=cnbrdn2
!future          endif
!future       else
!future          if(mod(me,2) == 0) then
!future             neighbor=nbrdn2
!future             cneighbor=cnbrdn2
!future          else
!future             neighbor=nbrup2
!future             cneighbor=cnbrup2
!future          endif
!future       endif
!future    else if(mod(istart-1,d1freq) == 0) then
!future       irexd1=irexd1+1
!future
!future       exchdim=1
!future       step=mod(irexd1,2)
!future       me=myrepd1
!future       if(step==1) then
!future          if(mod(me,2) == 0) then
!future             neighbor=nbrup1
!future             cneighbor=cnbrup1
!future          else
!future             neighbor=nbrdn1
!future             cneighbor=cnbrdn1
!future          endif
!future       else
!future          if(mod(me,2) == 0) then
!future             neighbor=nbrdn1
!future             cneighbor=cnbrdn1
!future          else
!future             neighbor=nbrup1
!future             cneighbor=cnbrup1
!future          endif
!future       endif
!future    else
!future       return
!future    endif
!future
!future    if(cneighbor > -1) then
!future       if(mynod == 0) then
!future          qexc = .false.
!future          write(outu,'(a,i3,a,i6,a,i5,a,i5)') 'DO2DEXCH> REPLICA ',IREPDSTR,' STEP ',istart-1, &
!future                                              ' EXCHANGE WITH NEIGHBOR ',neighbor, &
!future                                              ' PROCESSOR ',cneighbor
!future
!future          call chmalloc('repdstr.src','DO2DEXCH','w',natom,crl=w)
!future
!future          if(q2dtemp(exchdim)) then
!future             write(outu,'(a)') 'DO2DEXCH> DO TEMP EXCH'
!future
!future             ! exchange temperature with neighbor and calculate if exchange succeeded
!future             call masterrecsenrl(neighbor,1,eneigh,1,myepot,1)
!future             if(q2ditemp) then
!future                ourtemp=ttemp
!future             else
!future                ourtemp=temprx(me+1)
!future             endif
!future             call masterrecsenrl(neighbor,4,nbrtemp,1,ourtemp,1)
!future             p=exp(-(one/(kboltz*ourtemp))-(one/(kboltz*nbrtemp))*(eneigh-myepot))
!future             p=min(p,one)
!future             rn=random(repseed)
!future             qexc=rn < p
!future
!future             if(irepdstr < neighbor) then
!future                ! we control the exchange
!future                call mastersenint(neighbor,3,qexc,4)
!future                call mastersenint(neighbor,4,p,8)
!future                call mastersenint(neighbor,5,rn,8)
!future             else
!future                ! get our info from our natural superior
!future                call grec(cneighbor,3,qexc,4)
!future                call grec(cneighbor,4,p,8)
!future                call grec(cneighbor,5,rn,8)
!future             endif
!future
!future             if(qexc) then
!future                call masterrecsenrl(neighbor,2,w,natom,x,natom)
!future                x(1:natom)=w(1:natom)
!future                call masterrecsenrl(neighbor,2,w,natom,y,natom)
!future                y(1:natom)=w(1:natom)
!future                call masterrecsenrl(neighbor,2,w,natom,z,natom)
!future                z(1:natom)=w(1:natom)
!future                call masterrecsenrl(neighbor,2,w,natom,vx,natom)
!future                vx(1:natom)=w(1:natom)
!future                call masterrecsenrl(neighbor,2,w,natom,vy,natom)
!future                vy(1:natom)=w(1:natom)
!future                call masterrecsenrl(neighbor,2,w,natom,vz,natom)
!future                vz(1:natom)=w(1:natom)
!future
!future                scalef=sqrt(ttemp/nbrtemp)
!future                do i=1,natom
!future                   vx(i)=vx(i)*scalef
!future                   vy(i)=vy(i)*scalef
!future                   vz(i)=vz(i)*scalef
!future                enddo
!future
!future                if(qcrys.and.xdim > 0) then
!future                   call masterrecsenrl(neighbor,4,W,6,XTLABC,6)
!future                   xtlabc(1:6) = w(1:6)
!future                endif
!future                call masterrecsenrl(neighbor,11,w,natom,wmain,natom)
!future                wmain(1:natom) = w(1:natom)
!future
!future             endif
!future             write(iunrex,'(a,i10,a,f7.2,a,i3,a,f7.2,a,f8.2,a,f8.2,a,f6.3,a,f6.3,a,l1)') &
!future                   'T-REX>',istart-1,' TEMP ',ourtemp,' NBR ',neighbor,' NBRTMP ', &
!future                   nbrtemp, ' OURENE ',myepot,' NBRENE ',eneigh,' P ',p,' RN ', &
!future                   rn,' SUC? ',qexc
!future
!future          else if(q2dham(exchdim)) then
!future
!future             if(prnlev > 3) write(outu,'(a)') 'DO2DEXCH> DO HAMIL EXCH'
!future             qexc = .false.
!future
!future             epot1p = myepot
!future             call masterrecsenrl(neighbor,1,epot2q,1,epot1p,1)
!future
!future             ! Do test coordinate exchange --actual work is
!future             ! done down below.
!future             call masterrecsenrl(neighbor,2,w,natom,x,natom)
!future             x(1:natom)=w(1:natom)
!future             call masterrecsenrl(neighbor,2,w,natom,y,natom)
!future             y(1:natom)=w(1:natom)
!future             call masterrecsenrl(neighbor,2,w,natom,z,natom)
!future             z(1:natom)=w(1:natom)
!future
!future             if(qcrys.and.xdim > 0) then
!future                call masterrecsenrl(neighbor,4,W,6,XTLABC,6)
!future                xtlabc(1:6) = w(1:6)
!future             endif
!future             call masterrecsenrl(neighbor,11,w,natom,wmain,natom)
!future             wmain(1:natom) = w(1:natom)
!future             if(prnlev > 3) write(outu,'(a)') 'DO2DEXCH> DONE PART A'
!future
!future#if KEY_CONSPH==1
!future          else if(q2dph(exchdim)) then
!future
!future             if(prnlev > 3) write(outu,'(a)') 'DO2DEXCH> DO PH EXCH'
!future             ididphrex = 1
!future             qexc = .false.
!future             ph_l = phrx(me+1)
!future
!future             do i=1,nres
!future                if(tstate(i) == 1) iproto = iproto + 1
!future             enddo
!future             call masterrecsen(neighbor,6,jproto,1,iproto,1)
!future             call masterrecsenrl(neighbor,7,ph_m,1,ph_l,1)
!future
!future             ph_delta = log(10.0)*(ph_m - ph_l)*(iproto - jproto)
!future             if(ph_delta <= zero) then
!future                p=one
!future             else
!future                p=min(one,exp(-ph_delta))
!future             endif
!future             rn=random(repseed)
!future             qexc=rn < p
!future
!future             if(irepdstr < neighbor) then
!future                ! we control the exchange
!future                call mastersenint(neighbor,3,qexc,4)
!future                call mastersenint(neighbor,4,p,8)
!future                call mastersenint(neighbor,5,rn,8)
!future             else
!future                ! get our info from our natural superior
!future                call grec(cneighbor,3,qexc,4)
!future                call grec(cneighbor,4,p,8)
!future                call grec(cneighbor,5,rn,8)
!future             endif
!future
!future             if(qexc) then
!future                call chmalloc('repdstr.src','DO2DEXCH','nstate',nres,intg=nstate)
!future
!future                call masterrecsenrl(neighbor,2,w,natom,x,natom)
!future                x(1:natom)=w(1:natom)
!future                call masterrecsenrl(neighbor,2,w,natom,y,natom)
!future                y(1:natom)=w(1:natom)
!future                call masterrecsenrl(neighbor,2,w,natom,z,natom)
!future                z(1:natom)=w(1:natom)
!future                call masterrecsenrl(neighbor,2,w,natom,vx,natom)
!future                vx(1:natom)=w(1:natom)
!future                call masterrecsenrl(neighbor,2,w,natom,vy,natom)
!future                vy(1:natom)=w(1:natom)
!future                call masterrecsenrl(neighbor,2,w,natom,vz,natom)
!future                vz(1:natom)=w(1:natom)
!future
!future                if(qcrys.and.xdim > 0) then
!future                   call masterrecsenrl(neighbor,4,W,6,XTLABC,6)
!future                   xtlabc(1:6) = w(1:6)
!future                endif
!future                call masterrecsenrl(neighbor,11,w,natom,wmain,natom)
!future                wmain(1:natom) = w(1:natom)
!future
!future                call masterrecsenrl(neighbor,14,w,natom,cg,natom)
!future                cg(1:natom)=w(1:natom)
!future                call masterrecsen(neighbor,15,nstate,nres,tstate,nres)
!future                tstate(1:natom)=nstate(1:natom)
!future
!future                call chmdealloc('repdstr.src','DO2DEXCH','nstate',nres,intg=nstate)
!future             endif
!future             write(iunrex,'(a,i10,a,f5.2,a,i3,a,f5.2,a,f6.3,a,f6.3,a,l1)') &
!future                   'PHREX>',istart-1,' PH ',ph_m,' NBR ',neighbor,' NBRPH ', &
!future                   ph_l,' P ',p,' RN ',rn,' SUC? ',qexc
!future
!future#endif
!future          endif
!future          call chmdealloc('repdstr.src','DO2DEXCH','w',natom,crl=w)
!future       endif
!future
!future       ! this part is executed on all nodes in da group
!future
!future       call psnd4(qexc,1)
!future       call psnd4(exchdim,1)
!future
!future       if(q2dph(exchdim)) then
!future          ! Broadcast the charges, T-State, and other such shit
!future       endif
!future
!future       if(qexc .or. q2dham(exchdim)) then
!future          call psnd8(x,natom)
!future          call psnd8(y,natom)
!future          call psnd8(z,natom)
!future          call psnd8(wmain,natom)
!future          if(qcrys) then
!future             oldxtlabc(1:6)=xtlabc(1:6)
!future             call psnd8(xtlabc,6)
!future             call xtllat(xucell,xtlabc)
!future             call chmalloc('repdstr.src','REPXCHG','TRANSF',3,4,XNSYMM,crl=transf)
!future             call imfill(transf,.false.)
!future             call chmdealloc('repdstr.src','REPXCHGL','TRANSF',3,4,XNSYMM,crl=transf)
!future
!future             ! if we've done pH replica exchange, the charges might have changed.
!future             if(q2dph(exchdim)) call psnd8(cg,natom)
!future
!future             call upimag(x,y,z,wmain,0,x,y,z,vx,vy,vz)
!future          endif
!future          call nbonds(x,y,z,bnbnd,bimag)
!future          call energy(x,y,z,dx,dy,dz,bnbnd,bimag,0)
!future       endif
!future
!future       if(q2dham(exchdim)) then
!future
!future          if(prnlev > 3) write(outu,'(a)') 'DO2DEXCH> DO PART B'
!future
!future          ! This is unpleasant, because we have to do actual work here
!future          call psnd8(epot1p,1)
!future          call psnd8(epot2q,1)
!future          epot1q = eprop(epot) ! we calculated energy above
!future          
!future          ! Decide if any of this crap worked.
!future          if(mynod == 0) then
!future             call masterrecsenrl(neighbor,5,epot2p,1,epot1q,1)
!future
!future             if(irepdstr < neighbor) then
!future                p=exp(-(one/(kboltz*ttemp))*(epot2p+epot1q-epot2q-epot1p))
!future                p=min(p,one)
!future                rn=random(repseed)
!future                qexc=rn < p
!future
!future                call mastersenint(neighbor,1,rn,8)
!future                call mastersenint(neighbor,2,p,8)
!future                call mastersenint(neighbor,3,qexc,4)
!future             else
!future                call grec(cneighbor,1,rn,8)
!future                call grec(cneighbor,2,p,8)
!future                call grec(cneighbor,3,qexc,4)
!future             endif
!future
!future             if(prnlev > 3) write(outu,'(a)') 'DO2DEXCH> DONE PART B'
!future          endif
!future
!future          call psnd4(qexc,1)
!future          if(qexc) then
!future             if(prnlev > 3) then
!future                write(outu,'(a)') 'DO2DEXCH> START PART C1'
!future                call flush(outu)
!future             endif
!future             call chmalloc('repdstr.src','DO2DEXCH','w',natom,crl=w)
!future             call masterrecsenrl(neighbor,2,w,natom,vx,natom)
!future             vx(1:natom)=w(1:natom)
!future             call masterrecsenrl(neighbor,2,w,natom,vy,natom)
!future             vy(1:natom)=w(1:natom)
!future             call masterrecsenrl(neighbor,2,w,natom,vz,natom)
!future             vz(1:natom)=w(1:natom)
!future
!future             call psnd8(vx,natom)
!future             call psnd8(vy,natom)
!future             call psnd8(vz,natom)
!future             call chmdealloc('repdstr.src','DO2DEXCH','w',natom,crl=w)
!future             if(prnlev > 3) write(outu,'(a)') 'DO2DEXCH> DONE PART C1'
!future          else
!future          
!future             ! FML: we have to undo all of the hard work we just did
!future             if(prnlev > 3) then 
!future                write(outu,'(a)') 'DO2DEXCH> START PART C2'
!future                call flush(outu)
!future             endif
!future             if(mynod == 0) then
!future                call chmalloc('repdstr.src','DO2DEXCH','w',natom,crl=w)
!future                call masterrecsenrl(neighbor,2,w,natom,x,natom)
!future                x(1:natom)=w(1:natom)
!future                call masterrecsenrl(neighbor,2,w,natom,y,natom)
!future                y(1:natom)=w(1:natom)
!future                call masterrecsenrl(neighbor,2,w,natom,z,natom)
!future                z(1:natom)=w(1:natom)
!future
!future                if(qcrys.and.xdim > 0) then
!future                   call masterrecsenrl(neighbor,4,W,6,XTLABC,6)
!future                   xtlabc(1:6) = w(1:6)
!future                endif
!future                call masterrecsenrl(neighbor,11,w,natom,wmain,natom)
!future                wmain(1:natom) = w(1:natom)
!future                call chmdealloc('repdstr.src','DO2DEXCH','w',natom,crl=w)
!future             endif
!future
!future             call psnd8(x,natom)
!future             call psnd8(y,natom)
!future             call psnd8(z,natom)
!future             call psnd8(wmain,natom)
!future             if(qcrys) then
!future                oldxtlabc(1:6)=xtlabc(1:6)
!future                call psnd8(xtlabc,6)
!future                call xtllat(xucell,xtlabc)
!future                call chmalloc('repdstr.src','REPXCHG','TRANSF',3,4,XNSYMM,crl=transf)
!future                call imfill(transf,.false.)
!future                call chmdealloc('repdstr.src','REPXCHGL','TRANSF',3,4,XNSYMM,crl=transf)
!future                call upimag(x,y,z,wmain,0,x,y,z,vx,vy,vz)
!future             endif 
!future             call nbonds(x,y,z,bnbnd,bimag)
!future             call energy(x,y,z,dx,dy,dz,bnbnd,bimag,0)
!future             if(prnlev > 3) then
!future                write(outu,'(a)') 'DO2DEXCH> DONE PART C2'
!future                call flush(outu)
!future             endif
!future
!future          endif
!future          write(iunrex,'(a,i10,a,f8.2,a,f8.2,a,f8.2,a,f8.2,a,f6.3,a,f6.3,a,l1)') &
!future                'H-REX>',istart-1,' EPOT1(P)',epot1p,' EPOT1(Q) ',epot1q,' EPOT2(P) ', &
!future                epot2p,' EPOT2(Q) ',epot2q,' P ',p,' RN ',rn,' SUC? ',qexc
!future
!future       endif
!future
!future       if(qexc) then
!future          ! Yay; our work is done! Just reset the dynamics with new velocities and go...
!future
!future          jhstrt=0
!future          igvopt=2
!future
!future          call psnd4(jhstrt,1)
!future          call psnd4(igvopt,1)
!future          call psnd8(vx,natom)
!future          call psnd8(vy,natom)
!future          call psnd8(vz,natom)
!future       endif
!future
!future    else
!future       write(iunrex,'(a,i9,a)') '2DREX> ',istart-1,' SKIP EXCH'
!future    endif
!future
!future    ! broadcast everything within replica group
!future
!future  END SUBROUTINE DO2DEXCH
  
  subroutine repexchg(x,y,z,wmain,vx,vy,vz,xold,yold,zold,myepot,ttemp,istart,jhstrt, &
       iseed,iasvel,igvopt,callseq &
#if KEY_TSM==1
    ,backls &
#endif
    ,ididphrex &
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
    use parallel
    use mpi
    use repd_ensemble
    !
    real(chm_real) :: x(:), y(:), z(:),myepot,vx(*),vy(*),vz(*),wmain(*),ttemp
    real(chm_real) xold(*),yold(*),zold(*)

    integer istart,jhstrt,iseed,iasvel
#if KEY_TSM==1
    integer backls(*) 
#endif

    !
    logical qexc,qcrys,lused,lexattempt
    integer step,idecide,me,neighbor,igvopt,i,j,cneighbor,fneigh,callseq
    real(chm_real) eneigh,scaled,p,srate,rn,ttx,ttsgx,tfsgx,sgldarg,hrexarg
    real(chm_real),allocatable,dimension(:) :: w,comar,dd,timx,timy,timz
    real(chm_real),allocatable,dimension(:) :: oldwmain
    logical qhes, qdidrsvr
    real(chm_real) sgarray(10),sgarrayn(10),oldxtlabc(6)
    real(chm_real) scalsg,dteflf,dtcflf,fact
    real(chm_real) ecsum,ourpoti,ourpotj,nbrpoti,nbrpotj
    real(chm_real),allocatable,dimension(:,:,:) :: transf
#if KEY_PHMD==1
    real(chm_real) l(ntitr) !JAW
#endif
#if KEY_BLOCK==1
    integer k                                      !gg: msld-compatibility
    real(chm_real) n(nsitemld,nblock)              !gg: msld-compatibility
    real(chm_real) m(nsitemld*nblock)              !gg: msld-compatibility
    real(chm_real) thetavmlds(nsitemld*nblock)     !gg: msld-compatibility
    real(chm_real) thetamlds(nsitemld*nblock)      !gg: msld-compatibility
    real(chm_real) thetamldolds(nsitemld*nblock)   !gg: msld-compatibility
    real(chm_real) thetafmlds(nsitemld*nblock)     !gg: msld-compatibility
#endif
    integer n6,ical,oldrep

    integer                           :: iproto, jproto, ididphrex, phresstruct
    integer, allocatable,dimension(:) :: nstate
    real(chm_real)                    :: ph_l, ph_m, ph_delta
    integer :: ierror

    lexattempt=.false.
    ididphrex = 0
!#endif

    oldrep=reptag

    !future    if(q2drex) then
    !future       call do2dexch(x,y,z,wmain,vx,vy,vz,xold,yold,zold,myepot,ttemp,istart,jhstrt, &
    !future            iseed,iasvel,igvopt,callseq &
    !future#if KEY_TSM==1
    !future       ,backls &
    !future#endif
    !future#if KEY_CONSPH==1
    !future       ,ididphrex & 
    !future#endif
    !future       )
    !future       return
    !future    endif


    !
    !     When qpxtiger is on we are in a tiger pre-exchange state
    if(qpxtiger)then
       if(mod(istart-1,tigerneq) /= 0) return
       qrxtmin=.true.   ! we can do mini first now
    else
       if(mod(istart-1,irexfq) /= 0) return
    endif
    qcrys = (xtltyp /= '    ')

    !
    !     If we are about to do the exchange in the case of TIGER method
    !     we need some preparation:
    !       1. run number of iteration cycles of
    !          the pair of minimizer and equlibration
    !       2. then try for the exchange
    !
    !
    if(qrxtiger)then
       if(qrxtmin)then
          write(comlyn, &
               '(''mini abnr nstep '',i6,'' tolg '',f12.8,'' nprint 100'')'), &
               tigernm,tigergr
          comlen=51
          call maincomx(comlyn,comlen,lused)
          qrxtmin=.false.
       endif
       if(qpxtiger)then
          tigeriti=tigeriti+1
          if(tigeriti < tigerit)then
             return  ! maybe
          else
             qpxtiger=.false.
             tigeriti=0
          endif
       endif
    endif
    if(qrxsgld)then
       if(treflf<rsmall)then
!          call psetglob()
!          call psync()
!          call mpi_barrier(comm_charmm,ierror) 
!          call mpi_barrier(comm_master,ierror)
          trxlf=avgtlf
          call mpi_bcast(trxlf,1,mpi_double_precision,0,comm_master,ierror) !communicate between reps
          call mpi_bcast(trxlf,1,mpi_double_precision,0,comm_charmm,ierror) !communicate within rep
!          call psnd8(trxlf,1)
          trxlf=trxlf*temprx(irepdstr+1)/temprx(1)
!          call psetloc()
       endif
    endif
    !
    !     Prepare the variable NEIGHBOR.
    !     NEIGHBOR=-1 if no communication is needed on this process

    step=mod(irex,2)
    me=mod(irepdstr,2)

    if(step == 1)then
       neighbor=irepdstr+1
       if(me /= 0)neighbor=irepdstr-1
    else
       neighbor=irepdstr-1
       if(me /= 0)neighbor=irepdstr+1
    endif
    if(neighbor >= nrepdstr)neighbor=-1
    eneigh=zero                ! printout looks better this way
    ourpoti=zero
    ourpotj=zero
    nbrpoti=zero
    nbrpotj=zero
    qdidrsvr=.false.

    qexc=.false.
    rn=random(repseed)
    scaled=one
    scalsg=one

    exchanging: if(neighbor >= 0) then

       cneighbor=neighbor     !*numnod
       if(qrxsgld.and.mynod == 0)then
          ! let's take care of sgld before we start mucking
          ! about with the energies in the h-rex code below.
          sgarray(1)=epotlf
          sgarray(2)=epothf+epotlf
          sgarray(3)=(avgeflf*avgcflf-avgefhf*avgcfhf)/(kboltz*temprx(irepdstr+1))
          sgarray(4)=avgefhf*avgcfhf/(kboltz*temprx(irepdstr+1))
          sgarray(5)=repdid
          sgarray(6)=avgtlf
          sgarray(7)=avgeflf
          sgarray(8)=avgcflf
          call masterrecsenrl(neighbor,1,sgarrayn,8,sgarray,8)
          eneigh=sgarrayn(2)
          !p=min(one,exp(-(SGARRAY(3)-SGARRAYN(3))*(SGARRAYN(1)-SGARRAY(1))  &
          ! -(SGARRAY(4)-SGARRAYN(4))*(SGARRAYN(2)-SGARRAY(2))  &
          ! -(SGARRAY(3)*SGARRAYN(7)-SGARRAYN(3)*SGARRAY(7))*(SGARRAY(6)-SGARRAYN(6))))
          !p=min(one,exp(-(SGARRAY(3)-SGARRAYN(3))*(SGARRAYN(1)-SGARRAY(1))  &
          ! -(SGARRAY(4)-SGARRAYN(4))*(SGARRAYN(2)-SGARRAY(2))))

          if(.not.qrxtham) then
             sgldarg = -((sgarray(3)-sgarrayn(3))*(sgarrayn(1)-sgarray(1)) &
                  -(sgarray(4)-sgarrayn(4))*(sgarrayn(2)-sgarray(2)))
          endif
       else if(qrxtham) then

          ! ouch -- perform a hamiltonian coordinate swap inside temperature rex
          !
          ! We're doing this up here because we need to perform the test swap and
          ! energy eval BEFORE telling other processors to get lost.
          if(mynod == 0) then
             call masterrecsenrl(neighbor,1,eneigh,1,myepot,1)
             ourpoti=myepot
             nbrpotj=eneigh

             !     perform a test coordinate exchange and calculate the new energies
             call chmalloc('repdstr.src','repxchg','w',natom,crl=w)
             call chmalloc('repdstr.src','repxchg','timx',natom,crl=timx)
             call chmalloc('repdstr.src','repxchg','timy',natom,crl=timy)
             call chmalloc('repdstr.src','repxchg','timz',natom,crl=timz)
             call chmalloc('repdstr.src','repxchg','oldwmain',natom,crl=oldwmain)

             timx(1:natom)=x(1:natom)
             timy(1:natom)=y(1:natom)
             timz(1:natom)=z(1:natom)
             oldwmain(1:natom)=wmain(1:natom)

             call masterrecsenrl(neighbor,2,w,natom,x,natom)
             x(1:natom) = w(1:natom)
             call masterrecsenrl(neighbor,2,w,natom,y,natom)
             y(1:natom) = w(1:natom)
             call masterrecsenrl(neighbor,2,w,natom,z,natom)
             z(1:natom) = w(1:natom)
             if(qcrys.and.xdim > 0) then
                call masterrecsenrl(neighbor,4,w,6,xtlabc,6)
                xtlabc(1:6) = w(1:6)
             endif
             call masterrecsenrl(neighbor,11,w,natom,wmain,natom)
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

          if(mynod == 0) then

             CALL MASTERRECSENRL(NEIGHBOR,1,ENEIGH,1,ECSUM,1)
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

             call chmdealloc('repdstr.src','repxchg','w',natom,crl=w)
             call chmdealloc('repdstr.src','repxchg','timx',natom,crl=timx)
             call chmdealloc('repdstr.src','repxchg','timy',natom,crl=timy)
             call chmdealloc('repdstr.src','repxchg','timz',natom,crl=timz)
             call chmdealloc('repdstr.src','repxchg','oldwmain',natom,crl=oldwmain)

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

       endif
    endif exchanging

    masterwork1: if(mynod /= 0) then

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


       neighborminus: IF(NEIGHBOR == -1) THEN
          IRESEXCH=IRESEXCH+1
          IF(QRESERVOIR) THEN

             IF((IREPDSTR == 0).AND.QRESLOW) THEN
                QDIDRSVR=.TRUE.
                IF(PRNLEV >= 6) &
                     WRITE(IUNREX,*) 'REPEXCH> BOT REPL WOULD EXCH W/ RESERVOIR'
#if KEY_CONSPH==1
                IF(QPHREX) THEN
                   IDIDPHREX=1 
                   IPROTO=0
                   DO I=1,NRES
                      IF(TSTATE(I) == 1) IPROTO=IPROTO+1
                   ENDDO
                   CALL RESEXCH_PH(TTEMP,.FALSE.,X,Y,Z,VX,VY,VZ,PHRX(IREPDSTR+1),IPROTO, &
                        ISEED,IASVEL,IGVOPT,P,RN,JPROTO,PH_M,QEXC,FNEIGH &
#if KEY_TSM==1
                   ,BACKLS &
#endif
                   )
                ELSE
#endif
                   CALL RESEXCH(.FALSE.,X,Y,Z,VX,VY,VZ,QEXC,TEMPRX(IREPDSTR+1),MYEPOT, &
                        ISEED,IASVEL,IGVOPT,P,RN,ENEIGH,FNEIGH &
#if KEY_TSM==1
                   ,BACKLS &
#endif
                   )
#if KEY_CONSPH==1
                ENDIF
#endif
                IF(PRNLEV >= 6) WRITE(IUNREX,'(A,I6,3F12.6,A,l1)') &
                     'TIM DBG> FNEIGH,ENEIGH,P,RN = ', FNEIGH, ENEIGH, P, RN, ' SUCCESS = ', QEXC
             ENDIF
             IF((IREPDSTR == NREPDSTR-1).AND.QRESHIGH) THEN
                QDIDRSVR=.TRUE.
                IF(PRNLEV >= 6) &
                     WRITE(IUNREX,*) 'REPEXCH> TOP REPL WOULD EXCH W/ RESERVOIR'
#if KEY_CONSPH==1
                IF(QPHREX) THEN
                   IDIDPHREX=1
                   IPROTO=0
                   DO I=1,NRES
                      IF(TSTATE(I) == 1) IPROTO=IPROTO+1
                   ENDDO
                   CALL RESEXCH_PH(TTEMP,.TRUE.,X,Y,Z,VX,VY,VZ,PHRX(IREPDSTR+1),IPROTO, & 
                        ISEED,IASVEL,IGVOPT,P,RN,JPROTO,PH_M,QEXC,FNEIGH &
#if KEY_TSM==1
                        ,BACKLS & 
#endif
                   )
                ELSE
#endif
                   CALL RESEXCH(.TRUE.,X,Y,Z,VX,VY,VZ,QEXC,TEMPRX(IREPDSTR+1),MYEPOT, &
                        ISEED,IASVEL,IGVOPT,P,RN,ENEIGH,FNEIGH & 
#if KEY_TSM==1
                   ,BACKLS &
#endif
                   )
#if KEY_CONSPH==1
                ENDIF
#endif

                IF(PRNLEV >= 6) WRITE(IUNREX,'(A,I6,3F12.6,A,l1)') &
                     'TIM DBG> FNEIGH,ENEIGH,P,RN = ', FNEIGH, ENEIGH, P, RN, ' SUCCESS = ', QEXC
             ENDIF
          ENDIF
          !  GOTO 99
       else neighborminus

#if KEY_CONSPH==1
          IF(QPHREX) THEN
             ! we are state i and have a pH value of pH_l, our neighbor
             ! is state j and has a pH value of pH_m

             ! count number of prototnated residues (in state 1) and swap with
             ! the neighbor
             ididphrex = 1
             iproto = 0
             do i=1,nres
                if(tstate(i) == 1) iproto = iproto + 1 
             enddo
             if(irepdstr > neighbor) then
                call masterrecint(neighbor,5,jproto,4)
                call mastersenint(neighbor,6,iproto,4)
             else
                call mastersenint(neighbor,5,iproto,4)
                call masterrecint(neighbor,6,jproto,4)
             endif

             ph_l = phrx(irepdstr+1)
             ph_m = phrx(neighbor+1)

             ph_delta = log(10.0)*(ph_m - ph_l)*(iproto - jproto)
             if(ph_delta <= zero) then
                p=one
             else
                p=min(one,exp(-ph_delta))
             endif

             IF(QSGLD.AND.MYNOD == 0)THEN
                ! we should be swapping SGLD stuffs too...
                SGARRAY(1)=EPOTLF
                SGARRAY(2)=EPOTHF+EPOTLF
                SGARRAY(3)=(AVGEFLF*AVGCFLF-AVGEFHF*AVGCFHF)/(kboltz*temprx(irepdstr+1))
                SGARRAY(4)=AVGEFHF*AVGCFHF/(kboltz*temprx(irepdstr+1))
                SGARRAY(5)=REPDID
                SGARRAY(6)=AVGTLF
                SGARRAY(7)=AVGEFLF
                SGARRAY(8)=AVGCFLF
                CALL MASTERRECSENRL(NEIGHBOR,1,SGARRAYN,8,SGARRAY,8)
                ENEIGH=SGARRAYN(2)
                !p=min(one,exp(-(SGARRAY(3)-SGARRAYN(3))*(SGARRAYN(1)-SGARRAY(1))  &
                ! -(SGARRAY(4)-SGARRAYN(4))*(SGARRAYN(2)-SGARRAY(2))  &
                ! -(SGARRAY(3)*SGARRAYN(7)-SGARRAYN(3)*SGARRAY(7))*(SGARRAY(6)-SGARRAYN(6))))
                !p=min(one,exp(-(SGARRAY(3)-SGARRAYN(3))*(SGARRAYN(1)-SGARRAY(1))  &
                ! -(SGARRAY(4)-SGARRAYN(4))*(SGARRAYN(2)-SGARRAY(2))))

                !! WE don't actually need this next line here

                !sgldarg = -((SGARRAY(3)-SGARRAYN(3))*(SGARRAYN(1)-SGARRAY(1)) &
                !          -(SGARRAY(4)-SGARRAYN(4))*(SGARRAYN(2)-SGARRAY(2)))
             ENDIF
          ELSE
#endif

             IF(QRXSGLD) THEN
                p=min(one,exp(sgldarg))
             ELSE
                IF(QRXTHAM) THEN
                   p=min(one,exp(hrexarg))
                ELSE
                   CALL MASTERRECSENRL(NEIGHBOR,1,ENEIGH,1,MYEPOT,1)
                   if(temprx(irepdstr+1) == temprx(neighbor+1)) then
                      p=ONE
                   else
                      P=MIN(ONE,EXP(-(ONE/(KBOLTZ*TEMPRX(IREPDSTR+1)) &
                           -ONE/(KBOLTZ*TEMPRX(NEIGHBOR+1)))*(ENEIGH-MYEPOT)))
                   endif
                ENDIF
             ENDIF

#if KEY_CONSPH==1
          ENDIF
#endif
          !
          !     QEXC would be the result of the probability test...
          !
          QEXC=P > RN
          !      if(prnlev >= 2)write(IUNREX,'(a,i5,l5)')
          !     $     'REPEXCHG>me,qexc=',mynodg,qexc
          !
          IF(MOD(IREPDSTR,2) == 0)THEN
             CALL MASTERRECINT(NEIGHBOR,3,QEXC,4)  !GG: Updates QEXEC to all other replicas?
             CALL MASTERRECINT(NEIGHBOR,4,RN,8)
          ELSE
             CALL MASTERSENINT(NEIGHBOR,3,QEXC,4)
             CALL MASTERSENINT(NEIGHBOR,4,RN,8)
          ENDIF
          !      if(prnlev >= 2)write(IUNREX,'(a,i5,l5)')
          !     $     'REPEXCHG>me,qexc=',mynodg,qexc
          !
          !     Perform the exchange of coordinates and velocities:
          IF(QEXC) THEN

             if(irepdstr > neighbor) then
                call masterrecint(neighbor,5,reptag,4)
                call mastersenint(neighbor,6,oldrep,4)
             else
                call mastersenint(neighbor,5,oldrep,4)
                call masterrecint(neighbor,6,reptag,4)
             endif
             !
             call chmalloc('repdstr.src','REPXCHG','W',NATOM,crl=w)
#if KEY_CONSPH==1
             if(qphrex) then
                call chmalloc('repdstr.src','REPEXCHG','NSTATE',nres,intg=nstate)
                scaled=one
             else
                scaled=sqrt(temprx(irepdstr+1)/temprx(neighbor+1))
             endif
#else
             scaled=sqrt(temprx(irepdstr+1)/temprx(neighbor+1))
#endif

             CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,X,NATOM)
             x(1:natom)=w(1:natom)
             CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,Y,NATOM)
             y(1:natom)=w(1:natom)
             CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,Z,NATOM)
             z(1:natom)=w(1:natom)
             CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,VX,NATOM)
             vx(1:natom)=w(1:natom)
             CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,VY,NATOM)
             vy(1:natom)=w(1:natom)
             CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,VZ,NATOM)
             vz(1:natom)=w(1:natom)

             if(qrxtham) then
                CALL MASTERRECSENRL(NEIGHBOR,11,W,NATOM,WMAIN,NATOM)
                wmain(1:natom) = w(1:natom)
                if(qcrys.and.xdim > 0) then
                   call masterrecsenrl(neighbor,4,W,6,XTLABC,6)
                   XTLABC(1:6) = W(1:6)
                endif
             endif

#if KEY_PHMD==1
             IF(QPHRX)THEN !JAW - if exchange then swap theta and velocities
                CALL masterrecsenrl(NEIGHBOR,2,L,NTITR,PH_THETA,NTITR)
                PH_THETA(1:ntitr) = l(1:ntitr)
                CALL masterrecsenrl(NEIGHBOR,2,L,NTITR,VPH_THETA,NTITR)
                VPH_THETA(1:ntitr) = l(1:ntitr)
                CALL masterrecsenrl(NEIGHBOR,2,L,NTITR,THETAOLD,NTITR)
                THETAOLD(1:ntitr) = l(1:ntitr)
                IF(PHBETA  >  ZERO) THEN ! Langevin PHMD
                   CALL masterrecsenrl(NEIGHBOR,2,L,NTITR,VPHOLD,NTITR)
                   VPHOLD(1:ntitr) = l(1:ntitr)
                   CALL masterrecsenrl(NEIGHBOR,2,L,NTITR,DPHOLD,NTITR)
                   DPHOLD(1:ntitr) = l(1:ntitr)
                ENDIF
             ENDIF !JAW
#endif
#if KEY_BLOCK==1
             IF (QMSPHRX) THEN
                write(outu,'(a)') 'Transferring Theta variables'
                K = nsitemld*nblock                          !GG: Calc total no. of elements in MSLD array
                THETAMLDS = RESHAPE(THETAMLD,SHAPE(M))       !GG: Reshape array for coordinates to 1D array
                CALL MASTERRECSENRL(NEIGHBOR,2,M,K,THETAMLDS,K)    !GG: Transfer array
                THETAMLDS(1:K) = M(1:K)
                THETAMLD = RESHAPE(THETAMLDS,SHAPE(N))       !GG: Reshape array back to nD array
                THETAMLDOLDS = RESHAPE(THETAMLDOLD,SHAPE(M)) !GG: Processing previous theta coordinates
                CALL MASTERRECSENRL(NEIGHBOR,2,M,K,THETAMLDOLDS,K)
                THETAMLDOLDS(1:K) = M(1:K)
                THETAMLDOLD = RESHAPE(THETAMLDOLDS,SHAPE(N))
                THETAVMLDS = RESHAPE(THETAVMLD,SHAPE(M))     !GG: Processing current velocity coordinates
                CALL MASTERRECSENRL(NEIGHBOR,2,M,K,THETAVMLDS,K)
                THETAVMLDS(1:K) = M(1:K)
                THETAVMLD = RESHAPE(THETAVMLDS,SHAPE(N))
                THETAFMLDS = RESHAPE(THETAFMLD,SHAPE(M))     !GG: Processing current force coordinates
                CALL MASTERRECSENRL(NEIGHBOR,2,M,K,THETAFMLDS,K)
                THETAFMLDS(1:K) = M(1:K)
                THETAFMLD = RESHAPE(THETAFMLDS,SHAPE(N))
             ENDIF
#endif

             CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,XOLD,NATOM)
             xOLD(1:natom)=w(1:natom)
             CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,YOLD,NATOM)
             yOLD(1:natom)=w(1:natom)
             CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,ZOLD,NATOM)
             zOLD(1:natom)=w(1:natom)
#if KEY_CONSPH==1
             IF(QPHREX) THEN
                call masterrecsenrl(neighbor,2,w,natom,cg,natom)
                cg(1:natom)=w(1:natom)
                if(irepdstr > neighbor) then
                   call masterrecint(neighbor,5,nstate,4*nres)
                   call mastersenint(neighbor,6,tstate,4*nres)
                else
                   call mastersenint(neighbor,5,tstate,4*nres)
                   call masterrecint(neighbor,6,nstate,4*nres)
                endif
                tstate(1:nres)=nstate(1:nres)

                IF(QSGLD)THEN
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGVX,NATOM)
                   sgvx(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGVY,NATOM)
                   sgvy(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGVZ,NATOM)
                   sgvz(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGFX,NATOM)
                   sgfx(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGFY,NATOM)
                   sgfy(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGFZ,NATOM)
                   sgfz(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGGX,NATOM)
                   sggx(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGGY,NATOM)
                   sggy(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGGZ,NATOM)
                   sggz(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGHX,NATOM)
                   sghx(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGHY,NATOM)
                   sghy(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGHZ,NATOM)
                   sghz(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGKX,NATOM)
                   sgkx(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGKY,NATOM)
                   sgky(1:natom)=w(1:natom)
                   CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGKZ,NATOM)
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
#endif
             IF(QRXSGLD)THEN
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGVX,NATOM)
                sgvx(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGVY,NATOM)
                sgvy(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGVZ,NATOM)
                sgvz(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGFX,NATOM)
                sgfx(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGFY,NATOM)
                sgfy(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGFZ,NATOM)
                sgfz(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGGX,NATOM)
                sggx(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGGY,NATOM)
                sggy(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGGZ,NATOM)
                sggz(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGHX,NATOM)
                sghx(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGHY,NATOM)
                sghy(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGHZ,NATOM)
                sghz(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGKX,NATOM)
                sgkx(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGKY,NATOM)
                sgky(1:natom)=w(1:natom)
                CALL MASTERRECSENRL(NEIGHBOR,2,W,NATOM,SGKZ,NATOM)
                sgkz(1:natom)=w(1:natom)
             ENDIF
             !
             call chmdealloc('repdstr.src','REPXCHG','W',NATOM,crl=w)
#if KEY_CONSPH==1
             if(qphrex) call chmdealloc('repdstr.src','REPEXCHG','NSTATE',nres,intg=nstate)
#endif
             !
             !         if(prnlev >= 2)write(IUNREX,'(a,i5,2f20.8)')
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
                IF(PHBETA  >  ZERO)THEN
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
       ENDIF neighborminus
    endif masterwork1

    CALL PSND8(ENEIGH,1)
    CALL PSND8(P,1)
    CALL PSND8(RN,1)
    CALL PSND4(QEXC,1)
#if KEY_CONSPH==1
    IF(QPHREX) THEN
       CALL PSND4(ididphrex,1)
       CALL PSND4(iproto,1)
       CALL PSND4(jproto,1)
    ENDIF
#endif
    !     but we need these:
    !

    IF(QEXC)THEN
#if KEY_BLOCK==1
       call msld_checkvariables(2) !GG: Before transmission
#endif

       ! reset the dynamics algorithm if there's an exchange
       ! Do we need to set igvopt to 2 since this is what ASSVEL does?
       ! assvel does
       jhstrt=0
       igvopt=2
!       call mpi_barrier(comm_master,ierror)
!       call mpi_barrier(comm_master,ierror)
       call psnd4(igvopt,1)
       call psnd4(jhstrt,1)
       call psnd8(x,natom)
       call psnd8(y,natom)
       call psnd8(z,natom)
       call psnd8(vx,natom)
       call psnd8(vy,natom)
       call psnd8(vz,natom)
       call psnd8(xold,natom)
       call psnd8(yold,natom)
       call psnd8(zold,natom)

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

#if KEY_CONSPH==1
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
#endif

#if KEY_PHMD==1
       IF(QPHRX)THEN
          CALL PSND8(PH_THETA,NTITR)
          CALL PSND8(VPH_THETA,NTITR) 
          IF(PHBETA  >  ZERO)THEN
             CALL PSND8(VPHOLD,NTITR)
             CALL PSND8(DPHOLD,NTITR)
          ENDIF
       ENDIF
#endif

#if KEY_BLOCK==1
       if (qmsphrx) then   !gg: broadcast swapped data to other replicas
          k = nsitemld*nblock
          call psnd8(thetamlds,k)
          call psnd8(thetamldolds,k)
          call psnd8(thetavmlds,k)
          call psnd8(thetafmlds,k)
          if (.not. lmasternode) then
             thetamld = reshape(thetamlds,shape(n))
             thetamldold = reshape(thetamldolds,shape(n))
             thetavmld = reshape(thetavmlds,shape(n))
             thetafmld = reshape(thetafmlds,shape(n))
          endif
       endif
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
    if(prnlev >= 2) then
       !
       if(qexc)isuc=isuc+1
       srate=real(isuc)/real(irex)
#if KEY_CONSPH==1
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
#endif
          ttx=zero
          if(neighbor >= 0)ttx=temprx(neighbor+1)
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
#if KEY_CONSPH==1
       endif
#endif
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
          if(irepdstr == 0)then
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

#if KEY_CONSPH==1
       if(.not.qrxsgld.and..not.qphrex) then
#else
          if(.not.qrxsgld) then
#endif
             write(iunrex,'(a)') &
                  '------------- Replica Exchange End --------'
#if KEY_CONSPH==1
          else if(qphrex) then
             write(iunrex,'(a)') &
                  '------------- pH Replica Exchange End --------'
#endif
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


  subroutine garrettspacer
    return   !gg:  subroutine spacer to get bookmarks working correctly
  end subroutine garrettspacer


  subroutine fastrepexchgl(istart,iseed)
    use chm_kinds
    use number
    use consta
    use stream
    use repdstr
    use memory
    use clcg_mod,only: random
    use block_ltm, only: nblock,ninter,blcoep
    use lambdam
    use mpi
    use parallel
    use repd_ensemble
    implicit none
    integer :: istart
    integer :: iseed
    integer :: i
    integer,dimension(nrepdstr) :: prev_rep2node
    ! Energies before exchange, during exchange, and after rejection
    real(chm_real) :: energy1, energy2, energyn1, energyn2
    ! Variables in setneighbors
    integer :: me, step, neighbor, rep, neighborrep
    ! Variables in repexchgld_swap
    integer, allocatable, dimension(:) :: intbufferlong
    real(chm_real), allocatable, dimension(:) :: realbuffer, realbufferlong
    ! Variables in repexchgld_test
    real(chm_real) :: eneigh,eneighi,eneighj,epoti,epotj,p,rn
    logical :: qexc
    ! Variables in repexchgld_print
    character(len=120) :: fmt711,fmt712,fmt713,fmt714
    ! Variables in repexchgld_bcast
    integer :: ierr

    ! Before we do anything, fix the lambda values
    ! Lambda values not updated to match theta values until next energy call
    call msld_setblcoef(nblock,ninter,bixlam,blcoep)

    if(mod(istart-1,irexfq) /= 0) return

    call setneighbors   !contained routine

    if (lmasternode) then
      if(neighbor /= -1) then
        energy1 = 0
        call msld_add_potentialenergy(nblock,bixlam,energy1)

        call repexchgld_swap

        energy2 = 0
        call msld_add_potentialenergy(nblock,bixlam,energy2)

        call repexchgld_test

        if (.not. qexc) then
          call repexchgld_swap
        endif
      endif

      call repexchgld_remap

      if (irepdstr==0) then
        call repexchgld_print
      endif
    endif

    ! RLHFIX Still need to broadcast new biasing potential to other nodes
    call repexchgld_bcast

    return

  contains

    !------------------------------------------------------
    !       Set Neighbors
    !------------------------------------------------------
    subroutine setneighbors
      irex=irex+1
      !
      !     Prepare the variable NEIGHBOR.
      !     NEIGHBOR=-1 if no communication is needed on this process
      step=mod(irex,2)         !gg: what is the current mc step
      me=mod(map_node2rep(irepdstr+1),2)  !gg: is the current replica no. 
                                        !    even (step=0) or odd (step=1)?
      rep=map_node2rep(irepdstr+1)
      if(step == me)then
        neighborrep=rep+1
      else
        neighborrep=rep-1
      endif
      if (neighborrep == nrepdstr) then
        neighborrep=-1
      endif
      if (neighborrep >= 0) then
        neighbor=map_rep2node(neighborrep+1)
      else
        neighbor=-1
      endif

      ! print '(a,i10,a,i3,a,i3,a,i3,a,i3,a,i3)',"istart",istart,"step",step,"me",me,"irepdstr",irepdstr,"neighborrep",neighborrep,"neighbor",neighbor

      prev_rep2node=map_rep2node

      return
    end subroutine setneighbors

    !------------------------------------------------------
    !       Swap biasing potential with neighbor
    !------------------------------------------------------
    subroutine repexchgld_swap
      call chmalloc('repdstr.src','repxchgld','intbufferlong',nbiasv,intg=intbufferlong)
      call chmalloc('repdstr.src','repxchgld','realbuffer',nblock,crl=realbuffer)
      call chmalloc('repdstr.src','repxchgld','realbufferlong',nbiasv,crl=realbufferlong)

      call masterrecsenrl(neighbor,1,realbuffer,nblock,bielam,nblock)
      bielam=realbuffer
      
      call masterrecsenin(neighbor,2,intbufferlong,nbiasv,ipbias,nbiasv)
      ipbias=intbufferlong
      call masterrecsenin(neighbor,3,intbufferlong,nbiasv,ibclas,nbiasv)
      ibclas=intbufferlong

      call masterrecsenrl(neighbor,4,realbufferlong,nbiasv,irreup,nbiasv)
      irreup=realbufferlong
      call masterrecsenrl(neighbor,5,realbufferlong,nbiasv,irrlow,nbiasv)
      irrlow=realbufferlong
      call masterrecsenrl(neighbor,6,realbufferlong,nbiasv,ikbias,nbiasv)
      ikbias=realbufferlong

      call chmdealloc('repdstr.src','repxchgld','intbufferlong',nbiasv,intg=intbufferlong)
      call chmdealloc('repdstr.src','repxchgld','realbuffer',nblock,crl=realbuffer)
      call chmdealloc('repdstr.src','repxchgld','realbufferlong',nbiasv,crl=realbufferlong)
    end subroutine repexchgld_swap

    !-----------------------------------------------------------
    !        Masters decide on exchange
    !-----------------------------------------------------------
    subroutine repexchgld_test
      call masterrecsenrl(neighbor,1,energyn1,1,energy1,1)
      call masterrecsenrl(neighbor,1,energyn2,1,energy2,1)
      
      !GG: Somehow TEMNEW is not parsed in correctly in MSLD-CPHMD, using MSLD temp instead
      p=min(one,exp(-(one/(kboltz*tbld))*(energyn2+energy2-energyn1-energy1)))
      rn=random(repseed)
      qexc=p > rn

      if(me == 0)then
        call masterreclog(neighbor,3,qexc,1)
      else
        call mastersenlog(neighbor,3,qexc,1)
      endif

      if(me == step) then
         noppup=noppup+1
         if(qexc) nsucup=nsucup+1
      else
         noppdn=noppdn+1
         if(qexc) nsucdn=nsucdn+1
      endif
      
      fmt711    = '(a16,a4,i5,a9,i5,a5,i10)'     !!GG: writes out exchange results
      fmt712    = '(a15,f16.4,a15,f16.4)'
      fmt713    = '(a15,f16.4,a15,f16.4)'
      fmt714    = '(a8,f18.14)'
      IF (QEXC) THEN
        write(outu,'(a)') &
          '------------Replica Exchange (PH-REX) for CPHMD^MSLD------------'
        write(outu,fmt711) &
          ' PH-REXv2> ACCEPT ','REP', IREPDSTR,' NEIGHBOR',NEIGHBOR,'STEP',ISTART-1
        write(outu,fmt711) &
          ' PH-REXv2> ACCEPT ','REP', rep,' NEIGHBOR',neighborrep,'STEP',ISTART-1
      ELSE
        write(outu,'(a)') &
          '------------Replica Exchange (PH-REX) for CPHMD^MSLD------------'
        write(outu,fmt711) &
          ' PH-REXv2> REJECT ','REP', IREPDSTR,' NEIGHBOR',NEIGHBOR,'STEP',ISTART-1
        write(outu,fmt711) &
          ' PH-REXv2> REJECT ','REP', rep,' NEIGHBOR',neighborrep,'STEP',ISTART-1
      ENDIF
      write(outu,fmt712) ' Ei(REP)=',energy1,' Ej(Neighbor)=',energyn1  !GG: I/J here refers to pH
      write(outu,fmt713) ' Ei(neighbor)=',energy2,'Ej(rep)=',energyn2
      write(outu,fmt714) ' Ediff', (energyn2+energy2-energyn1-energy1)
      write(outu,fmt714) ' Prob', p
      write(outu,*) ' Success', qexc

      if(qexc) then
        rep=neighborrep
      endif

      return
    end subroutine repexchgld_test

    !-----------------------------------------------------------
    !        Send mapping to master of masters
    !-----------------------------------------------------------
    subroutine repexchgld_remap
      if(irepdstr == 0) then
        map_node2rep(1)=rep
        map_rep2node(rep+1)=0
      endif

      do i=1,nrepdstr-1
        if(irepdstr == 0) then
          call masterrecint(i,i,map_node2rep(i+1),1)
          map_rep2node(map_node2rep(i+1)+1)=i
        else if (irepdstr == i) then
          call mastersenint(0,i,rep,1)
        endif
      enddo

      do i=1,nrepdstr-1
        if(irepdstr == 0) then
          call mastersenint(i,1,map_node2rep,nrepdstr)
          call mastersenint(i,2,map_rep2node,nrepdstr)
        else if (irepdstr == i) then
          call masterrecint(0,1,map_node2rep,nrepdstr)
          call masterrecint(0,2,map_rep2node,nrepdstr)
        endif
      enddo

      return
    end subroutine repexchgld_remap

    !-----------------------------------------------------------
    !        Master of masters prints exchange file
    !-----------------------------------------------------------
    subroutine repexchgld_print
      if(ffrun) then
        write(iunrex,'(a)') '# replica temp. ener. neighbor ntemp nene prob p success? newrep'
        ffrun=.false.
      endif

      ! irex=irex+1 ! done earlier
      write(iunrex,'(a,i15,a,i12,a,i5)') '# Exchange ', IREX, ': STEP ', istart-1, ': REPEAT ', 1

      do i=0,nrepdstr-1
        step=mod(irex,2)
        me=mod(i,2)
        rep=i
        if(step == me)then
          neighborrep=rep+1
        else
          neighborrep=rep-1
        endif
        if (neighborrep == nrepdstr) then
          neighborrep=-1
        endif
        if (neighborrep >= 0) then
          neighbor=prev_rep2node(neighborrep+1)
        else
          neighbor=-1
        endif

        write(iunrex,'(i2,x,f12.6,x,f15.6,x,i2,x,f12.6,x,f15.6,x,f5.3,x,f5.3,x,l,x,i2)') &
          prev_rep2node(rep+1)+1,i+1.0,0.0,neighbor+1,neighborrep+1.0,0.0,0.0,0.0,(prev_rep2node(rep+1) /= map_rep2node(rep+1)),map_rep2node(rep+1)+1
          ! ourselves,ourtemp,epotarr(ourselves),neighbor,nbrtemp,epotarr(neighbor),prob,p,qexc,rep1
      enddo

      return
    end subroutine repexchgld_print

    !-----------------------------------------------------------
    !        Master of each replica broadcasts the new potential
    !-----------------------------------------------------------
    subroutine repexchgld_bcast
      call mpi_bcast(map_node2rep,nrepdstr,mpi_integer,0,comm_charmm,ierr)
      call mpi_bcast(map_rep2node,nrepdstr,mpi_integer,0,comm_charmm,ierr)

      call mpi_bcast(bielam,nblock,mpi_real8,0,comm_charmm,ierr)
      
      call mpi_bcast(ipbias,nbiasv,mpi_integer,0,comm_charmm,ierr)
      call mpi_bcast(ibclas,nbiasv,mpi_integer,0,comm_charmm,ierr)

      call mpi_bcast(irreup,nbiasv,mpi_real8,0,comm_charmm,ierr)
      call mpi_bcast(irrlow,nbiasv,mpi_real8,0,comm_charmm,ierr)
      call mpi_bcast(ikbias,nbiasv,mpi_real8,0,comm_charmm,ierr)

      return
    end subroutine repexchgld_bcast

  end subroutine fastrepexchgl


  subroutine fastrepexchgl_read_restart(u)
    use stream
    use repdstr
    use repd_ensemble
    implicit none
    integer, intent(in) :: u
    integer :: i,rep,status
    character(len=128) :: line

    ! Read the checkpoint file
    read(u,'(/A)',iostat=status) line
    if ( status/= 0 ) call wrndie(-3,'<fastrepexchgl>','EOF during read')
    read(u,'(1I12)') rep

    if(lmasternode)then
      ! Send target state to master of masters
      if(irepdstr == 0) then
        map_buffer(1)=rep  !! target map_node2rep
      endif

      do i=1,nrepdstr-1
        if(irepdstr == 0) then
          call masterrecint(i,i,map_buffer(i+1),1)
        else if (irepdstr == i) then
          call mastersenint(0,i,rep,1)
        endif
      enddo

      ! Send back from master of masters to all masters
      do i=1,nrepdstr-1
        if(irepdstr == 0) then
          call mastersenint(i,1,map_buffer,nrepdstr)
        else if (irepdstr == i) then
          call masterrecint(0,1,map_buffer,nrepdstr)
        endif
      enddo
    endif

  end subroutine fastrepexchgl_read_restart


  subroutine fastrepexchgl_restart_broadcast()
    use stream
    use repdstr
    use memory
    use block_ltm, only: nblock
    use lambdam
    use mpi
    use parallel
    use repd_ensemble
    implicit none
    integer :: i,j
    integer, allocatable, dimension(:) :: intbufferlong
    real(chm_real), allocatable, dimension(:) :: realbuffer, realbufferlong
    integer :: ierr

    if(lmasternode)then
      do i=nrepdstr-1,1,-1
        ! swap i
        if(map_node2rep(i+1) .ne. map_buffer(i+1))then
          if(irepdstr == i) then
            call repexchgld_swap(map_rep2node(map_buffer(i+1)+1))
          else if (irepdstr == map_rep2node(map_buffer(i+1)+1)) then
            call repexchgld_swap(i)
          endif
          map_node2rep(map_rep2node(map_buffer(i+1)+1)+1)=map_node2rep(i+1)
          map_node2rep(i+1)=map_buffer(i+1)
          do j=0,nrepdstr-1
            map_rep2node(map_node2rep(j+1)+1)=j
          enddo
        endif
      enddo

      do i=0,nrepdstr-1
        write(outu,*) 'REPEXCHGLD> Replica', i, 'is starting with biasing potential', map_node2rep(i+1)
        write(outu,*) 'REPEXCHGLD> as specified in checkpoint file'
      enddo
    endif

    call mpi_bcast(map_node2rep,nrepdstr,mpi_integer,0,comm_charmm,ierr)
    call mpi_bcast(map_rep2node,nrepdstr,mpi_integer,0,comm_charmm,ierr)

    call mpi_bcast(bielam,nblock,mpi_real8,0,comm_charmm,ierr)

    call mpi_bcast(ipbias,nbiasv,mpi_integer,0,comm_charmm,ierr)
    call mpi_bcast(ibclas,nbiasv,mpi_integer,0,comm_charmm,ierr)

    call mpi_bcast(irreup,nbiasv,mpi_real8,0,comm_charmm,ierr)
    call mpi_bcast(irrlow,nbiasv,mpi_real8,0,comm_charmm,ierr)
    call mpi_bcast(ikbias,nbiasv,mpi_real8,0,comm_charmm,ierr)

    return

  contains

    !------------------------------------------------------
    !       Swap biasing potential n
    !------------------------------------------------------
    subroutine repexchgld_swap(n)
      integer :: n
      
      call chmalloc('repdstr.src','repxchgld','intbufferlong',nbiasv,intg=intbufferlong)
      call chmalloc('repdstr.src','repxchgld','realbuffer',nblock,crl=realbuffer)
      call chmalloc('repdstr.src','repxchgld','realbufferlong',nbiasv,crl=realbufferlong)

      call masterrecsenrl(n,1,realbuffer,nblock,bielam,nblock)
      bielam=realbuffer
      
      call masterrecsenin(n,2,intbufferlong,nbiasv,ipbias,nbiasv)
      ipbias=intbufferlong
      call masterrecsenin(n,3,intbufferlong,nbiasv,ibclas,nbiasv)
      ibclas=intbufferlong

      call masterrecsenrl(n,4,realbufferlong,nbiasv,irreup,nbiasv)
      irreup=realbufferlong
      call masterrecsenrl(n,5,realbufferlong,nbiasv,irrlow,nbiasv)
      irrlow=realbufferlong
      call masterrecsenrl(n,6,realbufferlong,nbiasv,ikbias,nbiasv)
      ikbias=realbufferlong

      call chmdealloc('repdstr.src','repxchgld','intbufferlong',nbiasv,intg=intbufferlong)
      call chmdealloc('repdstr.src','repxchgld','realbuffer',nblock,crl=realbuffer)
      call chmdealloc('repdstr.src','repxchgld','realbufferlong',nbiasv,crl=realbufferlong)
    end subroutine repexchgld_swap

  end subroutine fastrepexchgl_restart_broadcast


!*******************************************************************************
!     REPEXCHGL   lambda exchange
!*******************************************************************************
  subroutine repexchgl(x,y,z,wmain,vx,vy,vz,eptt,temnew, &
       istart,iseed,iasvel,igvopt,jhstrt &
#if KEY_TSM==1
       ,backls                    &            
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
    use mpi
    use parallel
    use repd_ensemble

    !
    implicit none

    real(chm_real),dimension(natom) :: x,y,z,wmain,vx,vy,vz
    real(chm_real) :: lambda,lambda2,ecsum, eptt
    real(chm_real),dimension(lenent) :: emts(lenent)
    integer :: istart,oldrep,j
    real(chm_real) :: temnew,rn,rntmp,exratup,exratdn,rerat
    integer :: iseed,iasvel,igvopt,jhstrt

#if KEY_TSM==1
    integer :: backls(*)
#endif

    logical :: qexc,qcrys
    integer :: step,idecide,me,meg,neighbor,i,cneighbor,irepr,irepl
    real(chm_real) :: eneigh,eneighi,eneighj,epoti,epotj,p
    real(chm_real),allocatable,dimension(:) :: w
#if KEY_PHMD==1
    real(chm_real) :: l(ntitr) 
#endif
#if KEY_BLOCK==1
    !GG: MSLD-compatibility
    integer :: K                                  
    real(chm_real),dimension(nsitemld,nblock) :: n
    real(chm_real),dimension(nsitemld*nblock) :: &
         m,thetavmlds,thetamlds,thetamldolds,thetafmlds
    real(chm_real) ::ediff
#endif
    logical,allocatable,dimension(:) :: wl
    real(chm_real),allocatable,dimension(:,:,:) :: transf
    integer :: nbr,ierr,tag,nmsg
    logical :: qsend
    character(len=120) :: fmt711,fmt712,fmt713,fmt714

    if(mod(istart-1,irexfq) /= 0) return
    oldrep=reptag

    !-------------------------------------------------------------
    !        NEIGHBOR SETUP
    call setneighbors   !contained routine
    call set_right_left !contained routine

    if(neighbor == -1) then
       call end_rep_work ! when  rep=0 or nrep  and has no partner
       call finishup
       return
    endif

    qcrys = (xtltyp /= '    ')

    ! -----------------------------GG: Start of MYNOD  ==  0 processes -------

    masterwork_00: if (lmasternode) then
       
       call msld_checkvariables(1) !GG Variables Before MC Exchange only for block
       
       call masterrecsenrl(neighbor,1,eneigh,1,eptt,1)
       eneighi = eneigh
       epoti = eptt
       call chmalloc('repdstr.src','repxchgl','w',natom,crl=w)
       call repexchgl_transfer_coords    ! contained subroutine
       call repexchgl_transfer_theta_stuff    ! contained subroutine
       if (qcrys .and. xdim > 0) call repexchgl_transfer_xtlabc
       call chmdealloc('repdstr.src','repxchgl','w',natom,crl=w)
#if KEY_GCMC==1
       if (qgcmc) call repexchgl_transfer_gcmcon
#endif
    endif masterwork_00

#if KEY_BLOCK==1
    call msld_checkvariables(2) !GG: Variables Before Transmission
#endif
    call repexchgl_spread_coords_in_rep
#if KEY_PHMD==1
    if(qphrx)then
       ! jaw. after exchanging coordinates and lambda it is necesarry 
       !      to update charges and
       !      total system charge so energy calculations are correct. 
       call updatephmd(1, 1, 0, 1) 
    endif
#endif

    !============================================================
    !     Compute Energy
    !
    call repexchl_compute_energy_for_exchange   !contained subroutine

    !============================================================
    !     Test for Exchange criteria
    !
    if (lmasternode)   call repexchgl_master_exchange_test
    !
    !     If acceptance, keep the test exchange of coordinate;if rejection,
    !     perform the exchange again(go back to the configuration before exchange):

    call spread_result_to_slaves

    if(qexc) then   
       !------- EXCHANGE ACCEPTED ---------
       !     exchange velocities - more stable than assigning new velocities
       ! BTM -- this code is a gory mess, but since we've exchanged velocities, we need
       ! to restart the dynamics algorithm
       ! 
       jhstrt=0
       igvopt=2
       if(lmasternode) call masterwork_1
       call Distribute_after_exchange_to_slaves
    else            
       !------- EXCHANGE REJECTED -----------------
#if KEY_PHMD==1
       !JAW - Exchange was rejected--> we need to reset the non-bond lists. 
       !               also an energy call resets GB. THIS IS NECESSARY!!! 
       if(qphrx)then
          ! after exchanging coordinates and lambda it is necesarry 
          !   to update charges and
          !   total system charge so energy calculations are correct. JW
          call updatephmd(1, 1, 0, 1) 
          
          !???cb3 are the next two calls supposed to be inside te if(qphrx?)
          !???btm I believe that they are...
          call nbonds(x,y,z,bnbnd,bimag)  ! lnbnd,bimag,limag)
          call energy(x,y,z,dx,dy,dz,bnbnd,bimag,0)
       endif
#endif

       !------- EXCHANGE REJECTED go back to original state-----------------
       masterwork_3: IF (lmasternode) then
          call chmalloc('repdstr.src','repxchgl','w',natom,crl=w)
          call repexchgl_transfer_coords
          call repexchgl_transfer_theta_stuff_2
          if (qcrys .and. xdim > 0)call repexchgl_transfer_xtlabc
          call chmdealloc('repdstr.src','repxchgl','w',natom,crl=w)
#if KEY_GCMC==1
          if(qgcmc)call repexchgl_transfer_gcmcon
#endif
       endif masterwork_3

#if KEY_BLOCK==1
       call msld_checkvariables(2) !GG: Before transmission
#endif
       !     We need to broadcast the data within the same replica group

       call repexhgl_restore_within_replica
#if KEY_BLOCK==1
       IF (QMSPHRX) THEN
          !-- Toggle msld_setblcoef_fnexp to use thetamold values
          call msld_swapcoeff_mod     
          !-- Variables  After Swap, Before Call Energy
          call msld_checkvariables(3) 
       endif
#endif
       call energy(x,y,z,dx,dy,dz,bnbnd,bimag,0)
#if KEY_BLOCK==1
       if (qmsphrx) then
          call msld_swapcoeff_norm      !GG: Toggle msld_setblcoef_fnexp to use thetamold values
          call msld_checkvariables(4) !GG: Variables  After Swap, After Call Energy
       ENDIF
#endif
    endif

    if(qcrys) call upimag(x,y,z,wmain,0,x,y,z,vx,vy,vz)
    call nbonds(x,y,z,bnbnd,bimag)
    call energy(x,y,z,dx,dy,dz,bnbnd,bimag,0)
#if KEY_BLOCK==1
    if (qmsphrx) then
       call msld_checkvariables(5) !gg: variables  after swap, after mc exchange
    endif
#endif
    call finishup    !contained routine
    return

    !==========================================================
    !     Contained Subroutines finishup and setneighbors
    !==========================================================
  contains

    !-----------------------------------------------------------
    !      Compute energy for exchange criterion check
    !-----------------------------------------------------------
    subroutine repexchl_compute_energy_for_exchange
      call nbonds(x,y,z,bnbnd,bimag)  ! lnbnd,bimag,limag)
#if KEY_BLOCK==1
      !GG: Toggle msld_setblcoef_fnexp to use thetamold values
      if (qmsphrx) then
         call msld_swapcoeff_mod
         call msld_checkvariables(3) !gg: variables after swap, before call energy
      endif
#endif
      call energy(x,y,z,dx,dy,dz,bnbnd,bimag,0)
#if KEY_BLOCK==1
      IF (QMSPHRX) THEN
         call msld_swapcoeff_norm       !GG: Toggle msld_setblcoef_fnexp to use thetam values
         call msld_checkvariables(4) !GG: Variables  After Swap, After Call Energy
      endif
#endif
      if (qmsphrx) then
         !GG: Not summing ETERM(I), since PE from MSLD biases have no 
         !     individual ETERM(I) value
         if(prnlev >= 8) then
            write(outu,'(a)') 'Using MSLD-PHREX potential energy calculations'
         endif
      else
         !GG: Jana's implementation of traditional pH-REX CPHMD
         ecsum=zero
         do i=1,lenent
            ecsum=eterm(i)+ecsum
         enddo
         eprop(epot)=ecsum
      endif
    end subroutine repexchl_compute_energy_for_exchange

    !-----------------------------------------------------------
    !        Masters decide on exchange
    !-----------------------------------------------------------
    subroutine repexchgl_master_exchange_test
      nbr=neighbor
      call masterrecsenrl(nbr,1,eneigh,1,eprop(epot),1)
      eneighj = eneigh
      epotj = eprop(epot)
      
      if (qmsphrx) then
         !GG: Somehow TEMNEW is not parsed in correctly in MSLD-CPHMD, using MSLD temp instead
         !GG: I/J here refer to before/after swap
         p=min(one,exp(-(one/(kboltz*tbld))*(eneighj+epotj-eneighi-epoti)))
         !p=zero
      else
         p=min(one,exp(-(one/(kboltz*temnew))*(eneighj+epotj-eneighi-epoti)))
      endif
      rn=random(repseed)
      qexc=p > rn
      call write_result_of_exchange_test ! contained subroutine
      
      !-----------------------------------------------------------
      !     QEXC would be the result of the probability test...
      if(qexpt.or.qex2d)then
         if((step == 1).or.(step == 2))then
            if(mod(irepdstr,2) == 0)then
               nbr=neighbor
               call masterreclog(nbr,3,qexc,1)
            else
               nbr=neighbor
               call mastersenlog(nbr,3,qexc,1)
            endif
         else
            !         allocate(w(1))
            !         call masterrecsen(neighbor,9,w,1,qexc,1)
            tag=3
            nbr=neighbor
            if(meg == 0)then
               !  call mpi_recv(qexc,1,mpi_logical,nbr,tag,comm_master, ierr)
               call masterreclog(nbr,3,qexc,1)
            else
               ! call mpi_send(qexc,1,mpi_logical,nbr,tag,comm_master,ierr)
               call mastersenlog(nbr,3,qexc,1)
            endif
         endif
      else
         tag=3
         if(mod(irepdstr,2) == 0)then
            !call mpi_recv(qexc,1,mpi_byte,nbr,tag,comm_master, ierr)
            call masterreclog(neighbor,3,qexc,1)
         else
            call mastersenlog(neighbor,3,qexc,1)
         endif
         
         nbr=neighbor
         tag=9
         nmsg=1
         rntmp=rn
         qsend=(mod(irepdstr,2)==1)
         call masterrecsenrl_1(nbr,tag,rntmp,nmsg,rntmp,nmsg,qsend)
      endif
      
      if(neighbor > irepdstr) then
         noppup=noppup+1
         if(qexc) nsucup=nsucup+1
      else
         noppdn=noppdn+1
         if(qexc) nsucdn=nsucdn+1
      endif
      
#if KEY_BLOCK==1
      !GG: CPHMD^MSLD PH-REX printout done after the Call GREC/MASTERSENINT(NEIGHBOR,3,QEXC,4) commands
      !GG  If done before, QEXC will not be consistent since not updated?
      outresult_0: if (qmsphrx) then
         fmt711    = '(a16,a4,i5,a9,i5,a5,i10)'     !!GG: writes out exchange results
         fmt712    = '(a15,f16.4,a15,f16.4)'
         fmt713    = '(a15,f16.4,a15,f16.4)'
         fmt714    = '(a8,f18.14)'
         IF (QEXC) THEN
            write(outu,'(a)') &
                 '------------Replica Exchange (PH-REX) for CPHMD^MSLD------------'
            write(outu,fmt711) &
                 ' PH-REXv2> ACCEPT ','REP', IREPDSTR,' NEIGHBOR',NEIGHBOR,'STEP',ISTART-1
         ELSE
            write(outu,'(a)') &
                 '------------Replica Exchange (PH-REX) for CPHMD^MSLD------------'
            write(outu,fmt711) &
                 ' PH-REXv2> REJECT ','REP', IREPDSTR,' NEIGHBOR',NEIGHBOR,'STEP',ISTART-1
         ENDIF
         write(outu,fmt712) ' Ei(REP)=',EPOTI,' Ej(Neighbor)=',ENEIGHI  !GG: I/J here refers to pH
         write(outu,fmt713) ' Ei(neighbor)=',epotj,'Ej(rep)=',eneighj
         write(outu,fmt714) ' Ediff', ediff
         write(outu,fmt714) ' Prob', p
         write(outu,*) ' Success', qexc
      endif outresult_0
#endif
      return
    end subroutine repexchgl_master_exchange_test
    
    !-----------------------------------------------------------
    !        Communicate among masters after exchange decision
    !-----------------------------------------------------------
    subroutine masterwork_1
      call chmalloc('repdstr.src','repxchgl','w',natom,crl=w)
      call masterrecsenrl(neighbor,22,w,natom,vx,natom)
      vx(1:natom)=w(1:natom)
      call masterrecsenrl(neighbor,23,w,natom,vy,natom)
      vy(1:natom)=w(1:natom)
      call masterrecsenrl(neighbor,24,w,natom,vz,natom)
      vz(1:natom)=w(1:natom)
      call chmdealloc('repdstr.src','repxchgl','w',natom,crl=w) 
      if(irepdstr > neighbor) then
         call masterrecint(neighbor,25,reptag,4)
         call mastersenint(neighbor,26,oldrep,4)
      else
         call mastersenint(neighbor,25,oldrep,4)
         call masterrecint(neighbor,26,reptag,4)
      endif
#if KEY_PHMD==1
      if(qphrx)then ! jaw. exchange theta velocity terms, exchange was accepted
         call masterrecsenrl(neighbor,2,l,ntitr,vph_theta,ntitr)
         vph_theta(1:ntitr) = l(1:ntitr)
         call masterrecsenrl(neighbor,2,l,ntitr,thetaold,ntitr)
         thetaold(1:ntitr) = l(1:ntitr)
         if(phbeta  >  zero) then ! langevin phmd
            call masterrecsenrl(neighbor,2,l,ntitr,vphold,ntitr)
            vphold(1:ntitr) = l(1:ntitr)
            call masterrecsenrl(neighbor,2,l,ntitr,dphold,ntitr)
            dphold(1:ntitr) = l(1:ntitr)
         endif
      endif
#endif
#if KEY_BLOCK==1
      if (qmsphrx) then
         if(prnlev >= 6) then
            write(outu,'(a)') "EXCHANGE ACCEPTED, not transferring theta variables back" !GG: All theta variables prev transferred
         endif
      endif
#endif
      return
    end subroutine masterwork_1

    !-----------------------------------------------------------
    !       Distribute to slaves after exchange decision
    !-----------------------------------------------------------
    subroutine Distribute_after_exchange_to_slaves
      call mpi_bcast(vx,natom,mpi_double_precision,0,comm_charmm,ierr)
      call mpi_bcast(vy,natom,mpi_double_precision,0,comm_charmm,ierr)
      call mpi_bcast(vz,natom,mpi_double_precision,0,comm_charmm,ierr)
      call mpi_bcast(reptag,natom,mpi_integer,0,comm_charmm,ierr)
      !       CALL PSND8(VX,NATOM)
      !       CALL PSND8(VY,NATOM)
      !       CALL PSND8(VZ,NATOM)
      !       CALL PSND4(REPTAG,1)
      return
    end subroutine Distribute_after_exchange_to_slaves
    
    !-----------------------------------------------------------
    !            Transfer coords
    !-----------------------------------------------------------
    subroutine repexchgl_transfer_coords
      !     Perform a test coordinate exchange and calculate the new energies
      call masterrecsenrl(neighbor,2,w,natom,x,natom)
      x(1:natom) = w(1:natom)
      call masterrecsenrl(neighbor,2,w,natom,y,natom)
      y(1:natom) = w(1:natom)
      call masterrecsenrl(neighbor,2,w,natom,z,natom)
      z(1:natom) = w(1:natom)
#if KEY_PHMD==1
      if (qphrx) then ! jaw. exchange theta values before testing exchange  
         call masterrecsenrl(neighbor,2,l,ntitr,ph_theta,ntitr)
         ph_theta(1:ntitr) = l(1:ntitr)
      endif
#endif
      call masterrecsenrl(neighbor,11,w,natom,wmain,natom)
      wmain(1:natom) = w(1:natom)
      return
    end subroutine repexchgl_transfer_coords

    !-----------------------------------------------------------
    !            Transfer coords II
    !-----------------------------------------------------------
    subroutine repexchgl_transfer_coords_2
      call masterrecsenrl(neighbor,2,w,natom,x,natom)
      x(1:natom) = w(1:natom)
      call masterrecsenrl(neighbor,2,w,natom,y,natom)
      y(1:natom) = w(1:natom)
      call masterrecsenrl(neighbor,2,w,natom,z,natom)
      z(1:natom) = w(1:natom)
      call chmdealloc('repdstr.src','repxchgl','w',natom,crl=w)
      
#if KEY_PHMD==1
      if(qphrx)then ! jw: switch theta values back, exchange was rejected 
         call masterrecsenrl(neighbor,2,l,ntitr,ph_theta,ntitr)
         ph_theta(1:ntitr) = l(1:ntitr)
      endif
#endif
      return
    end subroutine repexchgl_transfer_coords_2
       
    !-----------------------------------------------------------
    !             Transfer theta stuff II
    !-----------------------------------------------------------
    subroutine repexchgl_transfer_theta_stuff_2
#if KEY_BLOCK==1
      if (qmsphrx) then
         if(prnlev >= 6) then
            write(outu,'(a)') &
                 "EXCHANGE REJECTED, transferring theta variables back"
         endif
         k = nsitemld*nblock                          
         !calc total no. of elements in msld array
         thetamlds = reshape(thetamld,shape(m))
         !reshape array for coordinates to 1d array
         call masterrecsenrl(neighbor,2,m,k,thetamlds,k)    
         !transfer array
         thetamlds(1:k) = m(1:k)
         thetamld = reshape(thetamlds,shape(n))       
         !reshape array back to nd array
         thetamldolds = reshape(thetamldold,shape(m)) 
         !processing previous theta coordinates
         call masterrecsenrl(neighbor,2,m,k,thetamldolds,k)
         thetamldolds(1:k) = m(1:k)
         thetamldold = reshape(thetamldolds,shape(n))
         thetavmlds = reshape(thetavmld,shape(m))     
         !processing current velocity coordinates
         call masterrecsenrl(neighbor,2,m,k,thetavmlds,k)
         thetavmlds(1:k) = m(1:k)
         thetavmld = reshape(thetavmlds,shape(n))
         thetafmlds = reshape(thetafmld,shape(m))     
         !processing current force coordinates
         call masterrecsenrl(neighbor,2,m,k,thetafmlds,k)
         thetafmlds(1:k) = m(1:k)
         thetafmld = reshape(thetafmlds,shape(n))
      endif
#endif
      return
    end subroutine repexchgl_transfer_theta_stuff_2


    !-----------------------------------------------------------
    !             Transfer theta stuff
    !-----------------------------------------------------------
    subroutine repexchgl_transfer_theta_stuff
#if KEY_BLOCK==1
      if (qmsphrx) then
         k = nsitemld*nblock         !gg: calc total no. of elements in msld array
         thetamlds = reshape(thetamld,shape(m))       !gg: reshape array for coordinates to 1d array
         nbr=neighbor
         call masterrecsenrl(nbr,2,m,k,thetamlds,k)    !gg: transfer array
         thetamlds(1:k) = m(1:k)
         thetamld = reshape(thetamlds,shape(n))       !gg: reshape array back to nd array
         thetamldolds = reshape(thetamldold,shape(m)) !gg: processing previous theta coordinates
         nbr=neighbor
         call masterrecsenrl(nbr,2,m,k,thetamldolds,k)
         thetamldolds(1:k) = m(1:k)
         thetamldold = reshape(thetamldolds,shape(n))
         thetavmlds = reshape(thetavmld,shape(m))     !gg: processing current velocity coordinates
         nbr=neighbor
         call masterrecsenrl(nbr,2,m,k,thetavmlds,k)
         thetavmlds(1:k) = m(1:k)
         thetavmld = reshape(thetavmlds,shape(n))
         thetafmlds = reshape(thetafmld,shape(m))     !gg: processing current force coordinates
         nbr=neighbor
         call masterrecsenrl(nbr,2,m,k,thetafmlds,k)
         thetafmlds(1:k) = m(1:k)
         thetafmld = reshape(thetafmlds,shape(n))
      endif
#endif
      return
    end subroutine repexchgl_transfer_theta_stuff

    !-----------------------------------------------------------
    !          Spread coords locally in rep
    !-----------------------------------------------------------
    subroutine repexchgl_spread_coords_in_rep
      call psnd8(x,natom)
      call psnd8(y,natom)
      call psnd8(z,natom)
      call psnd8(wmain, natom)
      if (qcrys) then
         call psnd8(xtlabc,6)
         call xtllat(xucell,xtlabc)
         call chmalloc('repdstr.src','repxchgl','transf',3,4,xnsymm,crl=transf) 
         call imfill(transf,.false.)
         call chmdealloc('repdstr.src','repxchgl','transf',3,4,xnsymm,crl=transf) 
         call upimag(x,y,z,wmain,0,x,y,z,vx,vy,vz)
      endif
#if KEY_BLOCK==1
      if (qmsphrx) then   !gg: broadcast swapped data to other replicas
         k = nsitemld*nblock
         call psnd8(thetamlds,k)
         call psnd8(thetamldolds,k)
         call psnd8(thetavmlds,k)
         call psnd8(thetafmlds,k)
         if (mynod == 0) then
            thetamld = reshape(thetamlds,shape(n))
            thetamldold = reshape(thetamldolds,shape(n))
            thetavmld = reshape(thetavmlds,shape(n))
            thetafmld = reshape(thetafmlds,shape(n))
         endif
      endif
#endif
#if KEY_GCMC==1
      if (qgcmc) then
         call psnd4(gcmcon,maxa)
      endif
#endif
      return
    end subroutine repexchgl_spread_coords_in_rep

    !-----------------------------------------------------------
    !          Spread result of exchange decision locally in rep
    !-----------------------------------------------------------
    subroutine spread_result_to_slaves    !contained
      call mpi_bcast(qexc,1,mpi_logical,0,comm_charmm,ierr)
      call mpi_bcast(noppup,1,mpi_integer,0,comm_charmm,ierr)
      call mpi_bcast(noppdn,1,mpi_integer,0,comm_charmm,ierr)
      call mpi_bcast(nsucup,1,mpi_integer,0,comm_charmm,ierr)
      call mpi_bcast(nsucdn,1,mpi_integer,0,comm_charmm,ierr)
      call mpi_bcast(p,1,mpi_double_precision,0,comm_charmm,ierr)
      
      if(prnlev > 0) then
         write(iunrex,'("H-REX> PROB ",f7.4," EXCH ",L1)') p,qexc
         call flush(iunrex)
      endif
      !    CALL PSND4(QEXC,1)
      !    CALL PSND4(NOPPUP,1)
      !    CALL PSND4(NOPPDN,1)
      !    CALL PSND4(NSUCUP,1)
      !    CALL PSND4(NSUCDN,1)
      !    CALL PSND8(P,1)
      return
    end subroutine spread_result_to_slaves
    
    !-----------------------------------------------------------
    !     Restoring state within replica after exchange reject
    !-----------------------------------------------------------
    subroutine repexhgl_restore_within_replica
      call psnd8(x,natom)
      call psnd8(y,natom)
      call psnd8(z,natom)
      call psnd8(vx,natom)
      call psnd8(vy,natom)
      call psnd8(vz,natom)
      call psnd8(wmain,natom)
#if KEY_GCMC==1
      if (qgcmc) then
         call psnd4(gcmcon,maxa)
      endif
#endif
      if (qcrys) then
         call psnd8(xtlabc,6)
         call xtllat(xucell,xtlabc)
         call chmalloc('repdstr.src','repxchgl','transf',3,4,xnsymm, &
              crl=transf)
         call imfill(transf,.false.)
         call chmdealloc('repdstr.src','repxchgl','transf',3,4,xnsymm, &
              crl=transf)
         !!call upimag(x,y,z,wmain,0,x,y,z,vx,vy,vz)
      endif
      
#if KEY_BLOCK==1
      if (qmsphrx) then   !gg: broadcast swapped data to other replicas
         k = nsitemld*nblock
         call psnd8(thetamlds,k)
         call psnd8(thetamldolds,k)
         call psnd8(thetavmlds,k)
         call psnd8(thetafmlds,k)
         if (mynod /= 0) then
            thetamld = reshape(thetamlds,shape(n))
            thetamldold = reshape(thetamldolds,shape(n))
            thetavmld = reshape(thetavmlds,shape(n))
            thetafmld = reshape(thetafmlds,shape(n))
         endif
      endif
#endif
      return
    end subroutine repexhgl_restore_within_replica

    !-----------------------------------------------------------
    !             Transfer xtlabc
    !-----------------------------------------------------------
    subroutine repexchgl_transfer_xtlabc
      nbr=neighbor
      call masterrecsenrl(nbr,4,w,6,xtlabc,6)
      xtlabc(1:6) = w(1:6)
      return
    end subroutine repexchgl_transfer_xtlabc
    !-----------------------------------------------------------
    !             Transfer gcmcon
    !-----------------------------------------------------------
    subroutine repexchgl_transfer_gcmcon
      call chmalloc('repdstr.src','repxchgl','wl',maxa,log=wl)
      nbr=neighbor
      call masterrecsenrl(nbr,9,wl,maxa,gcmcon,maxa)
      gcmcon(1:maxa) = wl(1:maxa)
      call chmdealloc('repdstr.src','repxchgl','wl',maxa,log=wl)
      return
    end subroutine repexchgl_transfer_gcmcon
   
   
    !-----------------------------------------------------------
    !             Write result of exchange test
    !-----------------------------------------------------------
    subroutine write_result_of_exchange_test
#if KEY_PHMD==1
      if (qphrx) then
         if (qexc) then
            write(iunrex,'(a)') &
                 '-------------------------------------------------------------'
            write(iunrex,'(a16,a4,i5,a9,i5,a5,i10) ') &
                 ' PH-REX> ACCEPT ','REP',ME,' NEIGHBOR',CNEIGHBOR,'STEP',ISTART-1
         ELSE 
            write(iunrex,'(a)') &
                 '-------------------------------------------------------------'
            write(iunrex,'(a16,a4,i5,a9,i5,a5,i10) ') &
                 ' PH-REX> REJECT ','REP',ME,'NEIGHBOR',CNEIGHBOR,'STEP',ISTART-1
         ENDIF
         write(iunrex,'(a15,f16.4,a15,f16.4)') ' Ei(REP)=',EPOTI,' Ej(Neighbor)=', &
              ENEIGHI
         write(iunrex,'(a15,f16.4,a15,f16.4)') ' Ei(Neighbor)=',EPOTJ,'Ej(REP)=', &
              ENEIGHJ       
      elseif ( .not. qmsphrx) then !gg: cphmd^msld ph-rex printout written later
      else
#endif
         if(prnlev > 0) then
            write(iunrex,'(a)') '------------------'
            write(IUNREX,"('H-REX> REPL ',I3,' NEIGHBOR ',I3,' STEP ',I8)") &
                 irepdstr,neighbor,istart-1
            write(IUNREX,'(a,f12.2,a,f12.2,a,f12.2,a,f12.2)') &
                 'H-REX> EPOT1(P)' ,epoti,' EPOT1(Q) ',epotj, &
                 ' EPOT2(Q) ',eneighi,' EPOT2(P) ',eneighj
         endif
#if KEY_PHMD==1
      endif
#endif
      return
    end subroutine write_result_of_exchange_test
    
    !-----------------------------------------------------------
    !         Set right left
    !-----------------------------------------------------------
    subroutine set_right_left
      if(qexbk)then
         irepr=irepdstr+1
         irepl=irepdstr-1
         if((irepr == irbk).and.(neighbor == irbk))neighbor=-1
         if((irepdstr == irbk).and.(neighbor == irepl))neighbor=-1
      endif
      return
    end subroutine set_right_left

    !-----------------------------------------------------------
    !         END_REP_WORK
    !-----------------------------------------------------------
    subroutine end_rep_work
      if(qreservoir) then
         if((irepdstr == 0).and.qreslow) then
            call resexchl(.false.,x,y,z,vx,vy,vz,temnew,rhtemp,eptt, &
                 iseed,iasvel,igvopt,p,istart-1,jhstrt &
#if KEY_TSM==1
                 ,backls & 
#endif
                 )
         endif
         if((irepdstr == nrepdstr-1).and.qreshigh) then
            call resexchl(.true.,x,y,z,vx,vy,vz,temnew,rltemp,eptt, &
                 iseed,iasvel,igvopt,p,istart-1,jhstrt &
#if KEY_TSM==1
                 ,backls &
#endif
                 )
         endif
      endif
      return
    end subroutine end_rep_work
      
    subroutine finishup
      if(prnlev > 0) then
         write(IUNREX,&
              '(a,I8,a,I2,a,I2,a,I2,a,F12.2,a,F12.2,a,F12.2,a,f12.2,a,f6.3,a,f6.3,a,L1)') &
              'SUMHREX> ',istart-1,' NBR ',neighbor, &
              ' ORTAG ',oldrep,' NWTAG ',reptag, &
              ' epot1(p) ',epoti,' epot1(q) ',epotj, &
              ' epot2(q) ',eneighi,' epot2(p) ',eneighj, &
              ' prob ',p,' rn ',rn,' exc ',qexc
         if(qsump) then
            if(noppup > 0) then
               exratup=real(nsucup)/real(noppup)
               write(IUNREX,"(a,I8,a,I8,a,F7.4)") &
                    'HREX UP> OPPORTUNITIES ',noppup, &
                    ' SUCCESSES ',nsucup,' RATIO ',exratup
            endif
            if(noppdn > 0) then 
               exratdn=real(nsucdn)/real(noppdn)
               write(IUNREX,'(a,I8,a,I8,a,F7.4)') &
                    'HREX DN> OPPORTUNITIES ',noppdn,' SUCCESSES ',nsucdn,' RATIO ',exratdn
            endif
            if(nsucup == 0.and.nsucdn == 0) then
               write(IUNREX,'(a)') 'WARNING! NO SUCCESSFUL EXCHANGES'
            else if(nsucup > 0) then
               rerat=exratdn/exratup
               if(noppdn > 0.and.(rerat > 2.0.or.rerat < 0.5)) then
                  write(IUNREX,'(a)') 'WARNING! UNBALANCED EXCHANGE RATIO; CONSIDER MOVING REPLICA OVERLAPS!'
               endif
            else if(nsucdn > 0) then ! defensive programming, FTW!
               rerat=exratup/exratdn
               if(noppup > 0.and.(rerat > 2.0.or.rerat < 0.5)) then
                  write(IUNREX,'(a)') 'WARNING! UNBALANCED EXCHANGE RATIO; CONSIDER MOVING REPLICA OVERLAPS!'
               endif
            endif
         endif
      endif
      return
    end subroutine finishup

    !------------------------------------------------------
    !       Set Neighbors
    !------------------------------------------------------
    subroutine setneighbors
      !
      !     Prepare the variable NEIGHBOR.
      !     NEIGHBOR=-1 if no communication is needed on this process
      !gg: cphmd^msld uses this protocol (exlm keyword only)
      type_a: if((.not.qexpt).and.(.not.qex2d))then 
         step=mod(istart/irexfq,2)     !gg: what is the current mc step
         me=mod(irepdstr,2)            !gg: is the current replica no. 
                                       !    even (step=0) or odd (step=1)?
         if(step == 1)then
            neighbor=irepdstr+1
            if(me /= 0)neighbor=irepdstr-1
         else
            neighbor=irepdstr-1
            if(me /= 0)neighbor=irepdstr+1
         endif
      else if (qexpt)then    type_a
         step=mod(istart/irexfq,4)
         me=mod(irepdstr,2)
         meg=mod(irepdstr/nrept,2)
         if(step == 1)then
            neighbor=irepdstr+1
            if(me /= 0)neighbor=irepdstr-1
            if((me == 0).and.(mod(neighbor,nrept) == 0))neighbor=-1
            if((me /= 0).and.(mod((neighbor+1),nrept) == 0))neighbor=-1
         else if(step == 2)then
            neighbor=irepdstr-1
            if(me /= 0)neighbor=irepdstr+1
            if((me /= 0).and.(mod(neighbor,nrept) == 0))neighbor=-1
            if((me == 0).and.(mod((neighbor+1),nrept) == 0))neighbor=-1
         else if(step == 3)then
            if(mod(irepdstr,nrept) == 0)then
               neighbor=irepdstr+nrept
               if(meg /= 0)neighbor=irepdstr-nrept
            else
               neighbor=-1
            endif
         else
            if(mod(irepdstr,nrept) == 0)then
               neighbor=irepdstr-nrept
               if(meg /= 0)neighbor=irepdstr+nrept
            else
               neighbor=-1
            endif
         endif
      else   type_a                     !gg: ex2d uses this???
         step=mod(istart/irexfq,4)       !gg: what is the current mc step
         me=mod(irepdstr,2)              !gg: is the current replica even (step=0) or odd (step=1)?
         meg=mod(irepdstr/nrepx,2)
         if(step == 1)then
            neighbor=irepdstr+1
            if(me /= 0)neighbor=irepdstr-1
            if((me == 0).and.(mod(neighbor,nrepx) == 0))neighbor=-1
            if((me /= 0).and.(mod((neighbor+1),nrepx) == 0))neighbor=-1
         else if(step == 2)then
            neighbor=irepdstr-1
            if(me /= 0)neighbor=irepdstr+1
            if((me /= 0).and.(mod(neighbor,nrepx) == 0))neighbor=-1
            if((me == 0).and.(mod((neighbor+1),nrepx) == 0))neighbor=-1
         else if(step == 3)then
            neighbor=irepdstr+nrepx
            if(meg /= 0)neighbor=irepdstr-nrepx
         else
            neighbor=irepdstr-nrepx
            if(meg /= 0)neighbor=irepdstr+nrepx
         endif
      endif type_a

      if(neighbor < 0)neighbor=-1
      if(neighbor >= nrepdstr)neighbor=-1
      return
    end subroutine setneighbors
    !-------------- End of contained subroutines -------------------
  end subroutine repexchgl


  !=========================================================
  !        DREPSETIO
  !=========================================================
  subroutine drepsetio(iol,prnl,wrnl)
    !-----------------------------------------------------------------------
    !     Set the IOLEV information from global values, ie all system
  use parallel
  use repdstr
    integer iol,prnl,wrnl

    !     this must be here since it is easier to protect it here then
    !     in the calling routines
    if (.not.qrepdstr) return

    ioon=1
    iolorig=iol
    prnlorig=prnl
    wrnlorig=wrnl
    !      write(70+mynodg,'(a,5i6)')'DREPSETIO-0>me,meg,iolev=',
    !     $     mynod,mynodg,iol,prnl
    !
    !     There is no problem here with IOLEV, ie always like this
    !     even when QRDQTT=.TRUE., since we deal with this separately
    if(mynod == 0)iol=1

    !     For QWRQTT=.TRUE. we can deal with it here:
    if(qwrqtt.and.(mynod == 0))prnl=5
    if(qwrqtt.and.(mynod == 0))wrnl=5
    !      write(70+mynodg,'(a,5i6)')'DREPSETIO-1>me,meg,iolev=',
    !     $     mynod,mynodg,iol,prnl
    return
  end subroutine drepsetio

  SUBROUTINE DREPRESIO(IOL,PRNL,WRNL)
    !-----------------------------------------------------------------------
    !     Restores the IOLEV information from global values, ie all system
    !
    !     Generalize this too!!!
    !
  use parallel
  use repdstr
    integer iol,prnl,wrnl

    if (.not.qrepdstr) return

    ioon=0
    !      write(50+mynodg,*)'DREPRESIO-0>me,meg,iolev=',
    !     $     mynod,mynodg,iol,prnl
    iol=iolorig
    prnl=prnlorig
    wrnl=wrnlorig
    !      write(50+mynodg,*)'DREPRESIO-1>me,meg,iolev=',
    !     $     mynod,mynodg,iol,prnl
    !
    return
  end subroutine drepresio

!MFC to be deleted  subroutine psetglob
!MFC to be deleted    !-----------------------------------------------------------------------
!MFC to be deleted    !     Set the parallel information from global values, ie all system
!MFC to be deleted  use parallel
!MFC to be deleted    integer:: i
!MFC to be deleted    mynod=mynodg
!MFC to be deleted    numnod=numnodg
!MFC to be deleted    mynodp=mynod+1
!MFC to be deleted    noddim=nptwo()
!MFC to be deleted    call cube(mynod,numnod,ippmap)
!MFC to be deleted    !     Lets do the INODE array, too:
!MFC to be deleted    do i=1,maxnode
!MFC to be deleted       inode(i)=mod(i,numnod)
!MFC to be deleted    enddo
!MFC to be deleted    return
!MFC to be deleted  end subroutine psetglob
!MFC to be deleted
!MFC to be deleted  subroutine psetloc
!MFC to be deleted    !-----------------------------------------------------------------------
!MFC to be deleted    !     Set the local parallel information.
!MFC to be deleted    !     Put PSETLOC,PSETGLOB into paral1.src ??
!MFC to be deleted  use parallel
!MFC to be deleted  use repdstr
!MFC to be deleted    integer i
!MFC to be deleted
!MFC to be deleted    numnod=numnodg/nrepdstr
!MFC to be deleted    mynod=mod(mynodg,numnod)
!MFC to be deleted    mynodp=mynod+1
!MFC to be deleted    noddim=nptwo()
!MFC to be deleted    do i = 0, numnod
!MFC to be deleted       ippmap(i)=i+numnod*irepdstr
!MFC to be deleted    enddo
!MFC to be deleted    !     Lets do the INODE array, too:
!MFC to be deleted    do i=1,maxnode
!MFC to be deleted       inode(i)=mod(i,numnod)
!MFC to be deleted    enddo
!MFC to be deleted    !  Also iparpt(), nparpt() have to be localized??? !!!! Pressure problems???
!MFC to be deleted    !     Maybe not so trivial??? Better save
!MFC to be deleted    !     both cases when generated and then copy them as needed ???
!MFC to be deleted    !  IPARPT and NPARPT are more tricky, since they depend on
!MFC to be deleted    !  groups. It can't be dealt generally since when repd is started we
!MFC to be deleted    !  might don't have the PSF yet. It is maybe better to use separate routine
!MFC to be deleted    !  and call it before vdgsum & vdgbr, or whatever routines use them.
!MFC to be deleted    !  or use nbonds command in the repd script after PSFs are set up.
!MFC to be deleted    !
!MFC to be deleted
!MFC to be deleted    return
!MFC to be deleted  end subroutine psetloc

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
     use dynutil,only: assvel

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
     IF(TRGT == 0) TRGT=1 ! Just in case the RNG returns 0
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

        IF(CURTMP == NEWT) THEN
           P=MIN(ONE,EXP(-(NEWE-CURENE)/(KBOLTZ*CURTMP)))
        ELSE
           P=MIN(ONE,EXP(-(ONE/(KBOLTZ*CURTMP) &
            -ONE/(KBOLTZ*NEWT))*(NEWE-CURENE)))
        ENDIF
        RN=RANDOM(REPSEED)
        QEXC=RN < P
        ENEIGH=NEWE
     ELSE IF(QRESNOBO) THEN
        IF(QHIGH) THEN
           NEWE=RHENER(TRGT)
        ELSE
           NEWE=RLENER(TRGT)
        ENDIF
        P=MIN(ONE,EXP(-(ONE/(KBOLTZ*CURTMP))*(NEWE-CURENE)))
        RN=RANDOM(REPSEED)
        QEXC=RN < P
        ENEIGH=NEWE
     ENDIF

     ! go ahead and swap in the new atomic coordinates, saving
     ! the old coordinates in rescrdx,y,z for possible adding to 
     ! the reservoir
     IF(QEXC) THEN

        ! actually read the new coordinates from the file
        I=((TRGT-1)*3)+1
        READ(UNUM,REC=I)   RESCRDX
        READ(UNUM,REC=I+1) RESCRDY
        READ(UNUM,REC=I+2) RESCRDZ

        IF(PRNLEV >= 6) &
          WRITE(IUNREX,'(A,I9,A)') '-- RESEXCH: SWAPPED IN ELT ', TRGT, ' OF THE RESERVOIR --'
        DO I=1,NATOM
           IF(PRNLEV >= 7) &
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
       
    IF(PRNLEV >= 6) &
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
     use dynutil,only: assvel

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

     if(mynod == 0) then
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
        if(trgt == 0) trgt=1 ! Just in case the RNG returns 0
        if(prnlev > 3) write(iunrex,'(a,i6)') 'RESEXCHL> TRGT = ', trgt
        call flush(iunrex)
        if(qhigh) then
           resepot=rhener(trgt)
        else
           resepot=rlener(trgt)
        endif

        ! Decide if we want to swap with the reservoir
        rn=random(repseed)
        if(temp == restemp) then
           p=one
        else
           p=min(one,exp(-(one/(kboltz*temp) &
                -one/(kboltz*restemp))*(ourepot-resepot)))
        endif
        qexc=p > rn

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
123  format('H-R-REX> ',i8,' REP ENE ',f10.4,' RES ENE ',f10.4,' PROB ',f7.4,' RN ',f7.4,' SUCCESS? ',L1)
     if(prnlev > 3) then
        write(iunrex,'(a)') '------------------'
        write(iunrex,123) step,ourepot,resepot,p,rn,qexc
        call flush(iunrex)
     endif

  END SUBROUTINE RESEXCHL


#if KEY_CONSPH==1
   subroutine resexch_ph(curtmp,qhigh,x,y,z,vx,vy,vz,ourph,iproto, &
                        iseed,iasvel,igvopt,p,rn,jproto,ph_m,qexc,resnum &
#if KEY_TSM==1
                        ,backls & 
#endif
                       )

     use psf
     use clcg_mod,only: random
     use stream
     use consta,only: kboltz
     use number,only: one
#if KEY_CONSPH==1
     use consph,only: tstate 
#endif
     use repdstr
     use image
     use bases_fcm
     use dynutil,only:assvel

     logical, intent(in)            :: qhigh
     logical, intent(out)           :: qexc
     integer, intent(in)            :: iproto
     integer, intent(out)           :: resnum,jproto
     real(chm_real), intent(out)    :: p,rn,ph_m
     real(chm_real), intent(inout)  :: x(*),y(*),z(*),vx(*),vy(*),vz(*)
     integer,intent(inout)          :: iseed,iasvel,igvopt
#if KEY_TSM==1
     integer,intent(inout)          :: backls(*) 
#endif
     real(chm_real), intent(in)     :: curtmp,ourph

     integer                        :: unum,ressz,trgt,i,j,k,pidx
     real(chm_real)                 :: ph_delta,tmp
     integer,dimension(natom)       :: tmptstate
     integer,dimension(:,:),pointer :: reststate

     if(qhigh) then
        unum=rhunit
        ressz=highressz
        ph_m=rhph
        reststate=>reshtstate
     else
        unum=rlunit
        ressz=lowressz
        ph_m=rlph
        reststate=>resltstate
     endif
     trgt=ceiling(random(repseed)*ressz)
     if(trgt == 0) trgt=1 ! just in case the rng returns 0
     resnum=trgt

     ! decide if we want to swap with the reservoir
     qexc = .false.

     ! for the boltzmann case at least, we need to read in the reservoir file
     ! so we can calculate jproto.
     i=((trgt-1)*5)+1
     read(unum,rec=i)   rescrdx
     read(unum,rec=i+1) rescrdy
     read(unum,rec=i+2) rescrdz
     read(unum,rec=i+3) rescg
     read(unum,rec=i+4) tmptstate
     jproto=0
     do j=1,nres
        if(tmptstate(j) == 1) jproto=jproto+1
     enddo
     !!!WRITE(OUTU,'(A,I5)') 'TIM DEBUG> COORDINATES OF RESERVOIR STRUCTURE ',TRGT
     !!!DO J=1,NATOM
     !!!   WRITE(OUTU,'(3F8.3)') RESCRDX(J),RESCRDY(J),RESCRDZ(J)
     !!!ENDDO
     write(outu,'(a,i5)') 'TIM DEBUG> PROTONATION COUNT OF RESERVOIR STRUCT = ', jproto
     call flush(outu)

     if(qresboltz) then
        ph_delta=log(10.0)*(ph_m-ourph)*(iproto-jproto)
        if(ph_delta <= 0) then
           p=one
        else
           p=min(one,exp(-ph_delta))
        endif
     else if(qresnobo) then
        call wrndie(-3,'<RESEXCH_PH>','NONBOLTZMANN REX NOT IMPLEMENTED FOR PH')
     endif
     rn=random(repseed)
     qexc=rn < p

     if(qexc) then
        if(prnlev >= 6) &
          WRITE(IUNREX,'(A,I9,A)') '-- RESEXCH: SWAP IN ELT ', TRGT, ' OF THE RESERVOIR --'

        do i=1,natom
           !!WRITE(IUNREX,'(A,I5,A,4XF10.3)') 'RESEXCH> CRD ', I, ' = ', RESCRDX(I), RESCRDY(I), RESCRDZ(I), RESCG(I)
           ! X
           tmp=x(i)
           x(i)=rescrdx(i)
           rescrdx(i)=tmp
           ! y
           tmp=y(i)
           y(i)=rescrdy(i)
           rescrdy(i)=tmp
           ! z
           tmp=z(i)
           z(i)=rescrdz(i)
           rescrdz(i)=tmp
           ! charge
           tmp=cg(i)
           cg(i)=rescg(i)
           rescg(i)=tmp

           if(cg(i) /= rescg(i)) then
              if(ntrans > 0) then
                 do k=natom+1,natim
                    pidx=bimag%imattr(k)
                    if(pidx > 0) then
                       cg(k)=cg(i)
                    endif
                 enddo
              endif
           endif
        enddo
        ! now we need to set the protonation states of each of the residues
        ! (via the tstate array)
        tstate(1:nres)=tmptstate(1:nres)
        !!write(iunrex,'(a,i5,a,i3)') 'resexch> set tstate(', i, ') = ', tstate(i)

        ! since we've reassigned the coordinates, we need to call assvel to
        ! reassign the velocities
        call assvel(curtmp,x,y,z,vx,vy,vz,amass,iseed,iasvel,igvopt,natom,imove &
#if KEY_TSM==1
                    ,backls & 
#endif
                   )
     endif
  end subroutine resexch_ph
#endif


  subroutine precalcene(qhigh,x,y,z,qecor,ecortemp)
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
     use mpi

     logical,intent(in)           :: qhigh,qecor
     real(chm_real),intent(in)    :: ecortemp
     real(chm_real),intent(inout) :: x(natom),y(natom),z(natom)
     integer                      :: unum,ressz,i,j,r
     
     real(chm_real), allocatable, dimension(:) :: tmpx, tmpy, tmpz
     real(chm_real)                            :: correction
     integer                                   :: ndegf, status

     if(qecor) then
        if(prnlev > 1) write(outu,'(a,f10.4)') 'PRECALCENE> ADDING ENERGY CORRECTION TERM AT TEMP ', ecortemp
        call get_ecor_adjust(correction,ecortemp)
        if(prnlev > 1) write(outu,'(a,f10.4)') 'PRECALCENE> ENERGY CORRECTION = ', correction
     else
        ndegf=-1
        correction=zero
     endif

     if(qhigh) then
        unum=rhunit
        ressz=highressz
        call chmalloc('repdstr.src','precalcene','rhener',ressz,crl=rhener)
     else
        unum=rlunit
        ressz=lowressz
        call chmalloc('repdstr.src','precalcene','rlener',ressz,crl=rlener)
     endif

     call chmalloc('repdstr.src','precalcene','tmpx',natom,crl=tmpx)
     call chmalloc('repdstr.src','precalcene','tmpy',natom,crl=tmpy)
     call chmalloc('repdstr.src','precalcene','tmpz',natom,crl=tmpz)

     do i=1,natom
        tmpx(i)=x(i)
        tmpy(i)=y(i)
        tmpz(i)=z(i)
     enddo
     do i=1,ressz
        ! NB, we have to read these into rescrd{x,y,z} because the sdcd
        ! stores them as 4 byte floats, not eight byte...
        !MFC why are rescrd being stored in a file?s
        r=(i*3)-2
        read(unum,rec=r)   rescrdx
        read(unum,rec=r+1) rescrdy
        read(unum,rec=r+2) rescrdz

        x(1:natom)=rescrdx(1:natom)
        y(1:natom)=rescrdy(1:natom)
        z(1:natom)=rescrdz(1:natom)

        ! Update coordinates and non-bond lists
        ! FIXME: We should call UPDATE here instead
        call mpi_bcast(x,natom,mpi_double_precision,0,comm_charmm,status)
        call mpi_bcast(y,natom,mpi_double_precision,0,comm_charmm,status)
        call mpi_bcast(z,natom,mpi_double_precision,0,comm_charmm,status)
        !        call psnd8(x,natom)
        !        call psnd8(y,natom)
        !        call psnd8(z,natom)
        !MFC does nbonds need to be called if energy will be called right after it?
        call nbonds(x,y,z,bnbnd,bimag)

        ! Get the energy
        ! The potential energy will be stored in the third element of the ETERM array
        call energy(x,y,z,dx,dy,dz,bnbnd,bimag,1)
        if(prnlev > 3) write(outu,'(a,i6,a,f12.6)') 'PRECALCENE> POTENTIAL OF STRUCT ', &
                                     i, ' = ', eprop(epot)

        if(qecor) then
           if(prnlev > 3) write(outu,'(a,f12.6)') 'PRECALCENE> CORRECTED ENE = ', eprop(epot)+correction
           if(qhigh) then
              rhener(i)=eprop(epot)+correction
           else
              rlener(i)=eprop(epot)+correction
           endif
        else
           if(qhigh) then
              rhener(i)=eprop(epot)
           else
              rlener(i)=eprop(epot)
           endif
        endif
     enddo
     ! Restore original coordinates
     x(1:natom)=tmpx(1:natom)
     y(1:natom)=tmpy(1:natom)
     z(1:natom)=tmpz(1:natom)
     call chmdealloc('repdstr.src','precalcene','tmpx',natom,crl=tmpx)
     call chmdealloc('repdstr.src','precalcene','tmpy',natom,crl=tmpy)
     call chmdealloc('repdstr.src','precalcene','tmpz',natom,crl=tmpz)
 
     ! reset the energy and force arrays to whatever they were before
     ! I'm not sure if this matters, but I'm doing it to be safe...
     !MFC check if the above statement has any validity. Do the coords need to be broadcast?
     !MFC     Does nbonds need to be called if energy is called (doesnt energy call nbonds)?
     !MFC     Does energy need to be called at all here?
     call mpi_bcast(x,natom,mpi_double_precision,0,comm_charmm,status)
     call mpi_bcast(y,natom,mpi_double_precision,0,comm_charmm,status)
     call mpi_bcast(z,natom,mpi_double_precision,0,comm_charmm,status)
!     CALL PSND8(X,NATOM)
!     CALL PSND8(Y,NATOM)
!     CALL PSND8(Z,NATOM)
     call nbonds(x,y,z,bnbnd,bimag)
     call energy(x,y,z,dx,dy,dz,bnbnd,bimag,1)

  end subroutine precalcene

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
        IF(PRNLEV > 1) WRITE(OUTU,'(A,F10.4)') 'GET_ENE_VAL> ADDING ENERGY CORRECTION TERM AT TEMP ', TEMP
        CALL GET_ECOR_ADJUST(CORRECTION,TEMP)
        IF(PRNLEV > 4) WRITE(OUTU,'(A,F10.4)') 'GET_ENE_VAL> ENERGY CORRECTION = ', CORRECTION
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
           if(prnlev > 3) write(outu,'(a,i4,a,f10.4)') 'GET_ENE_VAL> ENERGY OF HIGH RES STRUCT ', i, ' = ', rhener(i)
        else
           read(eneun,'(f9.4)') rlener(i)
           rlener(i)=rlener(i)+correction
           if(prnlev > 3) write(outu,'(a,i4,a,f10.4)') 'GET_ENE_VAL> ENERGY OF LOW RES STRUCT ', i, ' = ', rlener(i)
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
        IF(IMOVE(I) == 0) NDEGF=NDEGF+3
     ENDDO
     IF(QSHAKE) THEN
       DO I=1,NCONST
          IF((IMOVE(SHKAPR(1,I)) == 0).OR.(IMOVE(SHKAPR(2,I)) == 0)) THEN
             NDEGF=NDEGF-1  
          ENDIF
       ENDDO
     ENDIF  
     IF(PRNLEV > 1) WRITE(OUTU,'(A,I4)') 'GET_ECOR_ADJUST> NUMBER OF DEGREES OF FREEDOM = ', NDEGF

     IF(TEMP <= 0) &
        CALL WRNDIE(-3,'<GET_ECOR_ADJUST>','ENERGY CORRECTION REQUESTED WITH INVALID TEMPERATURE')
     C=(KBOLTZ*TEMP*NDEGF)/2.0
     IF(PRNLEV > 4) WRITE(OUTU,'(A,F10.4)') 'GET_ENE_VAL> ENERGY CORRECTION = ', C

   END SUBROUTINE GET_ECOR_ADJUST

#else /* repdstr_main */

  SUBROUTINE REPDSTRMAIN
    CALL WRNDIE(-1, &
         '<REPDSTR>','REPlica DiSTRibuted code not compiled.')
    RETURN
  END SUBROUTINE REPDSTRmain
#endif /* repdstr_main */


END MODULE REPDSTRMOD2
