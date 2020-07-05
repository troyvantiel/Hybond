PROGRAM CHARMM
  use new_timer,only:init_timers,timer_start,t_total              
  !
  !      Chemistry at HARvard Macromolecular Mechanics
  !      -            ---     -              -
  !
  !      Version 44 - Developmental Version (c44b2) - February 15, 2020
  !
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !                                                                      C
  !      COPYRIGHT(c) 1984-2020                                          C
  !      President and Fellows of Harvard College                        C
  !                                                                      C
  !      All rights reserved                                             C
  !                                                                      C
  !      This copyright includes all source files and documentation.     C
  !                                                                      C
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !
  !
  !      Present and past developers:
  !      Ioan Andricioaei
  !      Georgios Archontis
  !      Jay L. Banks
  !      Christian Bartels
  !      Robert Best
  !      Arnaud Blondel
  !      Stefan Boresch
  !      John Brady
  !      Bernard R. Brooks
  !      Charles L. Brooks III
  !      Robert E. Bruccoleri
  !      Axel Brunger
  !      Amedeo Caflisch
  !      Leo Caves
  !      Jhih-Wei Chu
  !      Michael Crowley
  !      Qiang Cui
  !      Ryszard Czerminski
  !      Aaron R. Dinner
  !      Ron Elber
  !      Marcus Elstner
  !      Jeff Evanseck
  !      Michael Feig
  !      Scott Feller
  !      Martin J. Field
  !      Stefan Fischer
  !      Stephen H. Fleischman
  !      Jiali Gao
  !      Michael Garrahan
  !      Bruce Gelin
  !      Garrett Goh
  !      Urs Haberthuer
  !      Thomas A. Halgren
  !      Sergio A. Hassan
  !      Milan Hodoscek
  !      Antti-Pekka Hynninen
  !      Toshiko Ichiye
  !      Wonpil Im
  !      Mary E. Karpen
  !      Jennifer Knight
  !      Jeyapandian Kottalam
  !      Krzysztof Kuczera
  !      Themis Lazaridis
  !      Michael S. Lee
  !      Paul Lyne
  !      Jianpeng Ma
  !      Alex Mackerell
  !      Paul Maragakis
  !      Markus Meuwly
  !      Robert Nagle
  !      Benjamin Tim Miller
  !      Tom Ngo
  !      Lennart Nilsson
  !      Barry D. Olafson
  !      Victor Ovchinnikov
  !      Emanuele Paci
  !      Richard W. Pastor
  !      David Perahia
  !      Robert J. Petrella
  !      B. Montgomery Pettitt
  !      Carol B. Post
  !      Zingzhi Pu
  !      Walter E. Reiher III
  !      Benoit Roux
  !      Michael Schaefer
  !      Jana Shen
  !      Paul Sherwood
  !      Tom Simonson
  !      Jeremy Smith
  !      David J. States
  !      Peter J. Steinbach
  !      Roland Stote
  !      John Straub
  !      Sundaramoothi Swaminathan
  !      Walter Thiel
  !      Bruce Tidor
  !      Douglas J. Tobias
  !      Don G. Truhlar
  !      Arjan van der Vaart
  !      Richard M. Venable
  !      Herman van Vlijmen
  !      Joanna Wiorkiewicz
  !      Masa Watanabe
  !      Youngdo Won
  !      Lee H. Woodcock
  !      Thomas B. Woolf
  !      Xiongwu Wu
  !      Wei Yang
  !      Darrin M. York
  !      William S. Young
  !
  !      Developed under the overall direction of Martin Karplus,
  !      Department of Chemistry & Chemical Biology, Harvard University,
  !      12 Oxford Street, Cambridge, MA 02138
  !
  !      Before including in a publication any data calculated using
  !      CHARMM, please contact Martin Karplus at the address above to get
  !      the appropriate publication(s) to be referenced.
  !
  !      Refer to the documentation for information and usage.
  !
  !      DIMENSIONING INFORMATION - SEE PARTICULAR .FCM FILES
  !-----------------------------------------------------------------------
  !
  use chm_kinds
  use chm_types
  use dimens_fcm
  use number
  use bases_fcm
  use comand
  use ctitla
  use image
  use inbnd
  use param_store, only: param_store_init, set_param
  use param
  use pert
  use psf
  use startup
  use stream
  use string
  use timerm
  use parallel
#ifdef KEY_MPI
  use mpi
#endif
  use repdstr
#if KEY_REPDSTR==1
  use REPDSTRMOD                           
#endif
#if KEY_REPDSTR2==1
  use REPDSTRMOD2
  use repd_ensemble
#endif
  use intcor_module,only:initialize_icr_structs
  use machutil,only:Initialize_timers
  use machio,only:vopen
  use vangle_mm, only: ptrini
#if KEY_CFF==1
  use rtf,only: ucase    
#endif
  use usermod,only: usrini
  use cmdpar,only:cmdpar_init
#if KEY_ENSEMBLE==1
  use ensemble     
#endif
#if KEY_GAMUS==1
  use gamusmodule,only:gamusinit 
#endif
#if KEY_DHDGB==1
!AP/MF
  use dhdgb
  use derivdhdgb
#endif

  ! work-around for intel 18 compiler bug affecting fsystem
#if (__INTEL_COMPILER == 1800) && defined(_OPENMP)
  use ifport, only: setenvqq
#endif
  
  implicit none

  logical quantaok    

  !     The following are local variables.
  integer   istart
  logical   eof,error
  logical   lused,ok,oktmp
  logical   wantqu
  logical qrdcmd,get_next_cmd

  !--mfc temporary variables
  integer :: j,ierr
  character(len=4) :: winit

  ! work-around for intel 18 compiler bug affecting fsystem
#if (__INTEL_COMPILER == 1800) && defined(_OPENMP)
  error = setenvqq("KMP_INIT_AT_FORK=FALSE")
#endif
  
  !
  !*********************************************************
  !     END DECLARATIONS -- EXECUTION BEGINS
  !*********************************************************
  !
  !     Set I/O units and zero mscpar totals; default to lower case file names
  outu=poutu
  prnlev=5
  lower=.true.
  bomlev=0
  iolev=1
  wrnlev=5

  call param_store_init()

  call set_param('BOMLEV',bomlev)
  call set_param('WRNLEV',wrnlev)
  call set_param('PRNLEV',prnlev)
  call set_param('IOLEV',iolev)
#if KEY_CFF==1
  ucase = .true. 
#endif

  !     Start times and do machine specific startup.
  call Startup_machine_dependent_code !used to be called from jobini
  call Initialize_timers   ! used to be jobini

  !     Get the CHARMM command line arguments and initialize for different
  !     platforms different variables. Check if CHARMM should communicate
  !     with QUANTA or not
  call cmdpar_init
  call argumt(wantqu)
  call set_dimens()

  call init_timers()
  call timer_start(T_total)

  !     open standard input on unit 5 and standard output on unit OUTU.
  call vopen(5,'$$INITIAL$$',' ',' ',error,0)

  !     print a header for the output.
  call header

  !     Initialize the index common blocks

  quantaok=.false.

  !     Open the input file and read title for run
  !     attention: 'call header' have already been done by now
  if (quantaok) then
     nstrm=0
     istrm=-1
     eof=.false.
  else
     nstrm=1
     istrm=5
     jstrm(nstrm)=istrm
     eof=.false.
     call tryoro(istrm,'FORMATTED')
     ntitla=0
     call rdtitl(titlea,ntitla,istrm,0)
  endif

  call initialize_icr_structs()
#if KEY_TSM==1
  call tsminit(.false.)       
#endif
  call iniall

  !     Initialize local variables
  call getpref()
  comlen=0
  altlen=0
  istart=1
  ! checking for dimension chsize command or related max.. sizes
  do while(comlen == 0)
     call rdcmnd(comlyn,mxcmsz,comlen,istrm,eof,.true.,.true., 'CHARMM> ')
#if KEY_ENSEMBLE==1
     call sav_comlyn(comlyn,comlen)     
#endif
  enddo
  call parse_size(comlyn,comlen,qrdcmd)



  call allocate_all
  !     Initialize the rest, simply because I cannot figure out right 
  !     now whether the stuff in gtnbct needs the allocation first.
  !     Eventually the following lines will probably move into  iniall
  !     above.
  j=4
  winit='INIT'
  call gtnbct(winit,j,bnbnd)


  !=======================================================================
  ! call user defined startup routine -- mfc pulled this from iniall
  call usrini
  !=======================================================================

  !     Initialize pointers
  call ptrini

  !=======================================================================
  !     START OF MAIN COMMAND LOOP
  !=======================================================================
  !     before reading the next command, make sure that the variables have not
  !     exceeded their bounds.

  main_cmd_loop: do while(.true.)
     ok = check_for_exceed_bounds()  !internal function
     !     Check for unparsed junk of the last command.
     if(qrdcmd) CALL XTRANE(COMLYN,COMLEN,'CHARMM')

     get_next_cmd=.true.
     miscom_cmd_loop:do while(get_next_cmd)

        !     If a timelimit has passed, execute the alternate commandline
        call chklim
        if (atlim) then
           !     Reset deadlines (otherwise we could have an infinite loop here)
              cpulim=0.0
              deadhr=-1
              if(prnlev >= 2) write(outu,*) &
                   'CHARMM: ',LIMTYP,' timelimit(s) exceeded'
              if (altlen > 0) then
                 if(prnlev >= 2) then
                    WRITE(OUTU,*) 'CHARMM: Executing alternate command'
                    WRITE(OUTU,*) 'ALTCOM: ',ALTCOM(1:ALTLEN)
                 endif
                 call copyst(comlyn,mxcmsz,comlen,altcom,altlen)
                 goto 120
              endif
              if(prnlev >= 2) write(outu,*) 'CHARMM: Terminating execution'
              call stopch('limit reached')
        endif

        !---------------------------------------------------------------
        !     main loop command reader
        if(qrdcmd) then
#if KEY_ENSEMBLE==1
           call ensprint(" CHM>> main loop"," ") 
#endif
           call rdcmnd(comlyn,mxcmsz,comlen,istrm,eof,.true.,.true., &
                'CHARMM> ')
#if KEY_ENSEMBLE==1
           call sav_comlyn(comlyn,comlen)     
#endif
#if KEY_ENSEMBLE==1
           call ensprint(" CHM>>>>> ",comlyn(1:len_trim(comlyn))) 
#endif
!           call psync_world
        endif

        qrdcmd = .true.
        
120     continue
        ! See if it is one of the miscellaneous commands...
        call miscom(comlyn,mxcmsz,comlen,lused)

        get_next_cmd = ( lused .and. (.not. eof) ) 
     end do miscom_cmd_loop

     eoftest: if (eof) then
        !     If we run out of stuff on a particular stream, pop to the
        !     previous stream. quit when there are no more streams.
140     call ppstrm(ok)

#if KEY_REPDSTR==1
        !     This is excuted by everyone, so we can broadcast here if process 0
        !     wants to finish. Always restore the local parallel setup: this is safe
        if(qrepdstr)then
           call psetglob
           call psnd4(ok,1)
           call psetloc
        endif
#endif 
#if KEY_REPDSTR2==1
        !     This is excuted by everyone, so we can broadcast here if process 0
        !     wants to finish. Always restore the local parallel setup: this is safe
        if(qrepdstr)then
           call MPI_ALLREDUCE (mpi_in_place,ok,1,MPI_LOGICAL,MPI_LOR,comm_master,ierr) 
        endif
#endif 

        if (.not.(ok)) then
           call stopch('END OF FILE')
        endif
        eof=.false.
        !     Don't loop back to read from a script if no script
        if(quantaok.and.istrm == 0) goto 140
        cycle main_cmd_loop
     endif eoftest

     ! parse the command
     call maincomx(comlyn,comlen,lused)

  end do main_cmd_loop
  !--------------------------------------------------------------------

contains
  logical function check_for_exceed_bounds() result(ok)
    ok=nseg <= maxseg .and. natom <= maxa .and. nbond <= maxb .and. &
         ntheta <= maxt .and. nphi <= maxp .and. nimphi <= maximp .and. &
         nnb <= maxnb .and. ndon <= maxpad .and. nacc <= maxpad .and. &
         nres <= maxres .and. natc <= maxatc .and. ncb <= maxcb .and. &
         nct <= maxct .and. ncp <= maxcp .and. nci <= maxci .and. &
         nch <= maxch .and. ncn <= maxcn &
#if KEY_CMAP==1
         .and. ncrterm <= maxcrt                              
#endif
#if KEY_CMAP==0
    ;                                                    
#endif
    IF (.NOT.(OK)) THEN
       IF(WRNLEV.GE.2) THEN
          ! Write out the list of conditions so that one can see what
          ! actually caused the error.
          write(outu,'(a17,l5,2i10)')  'NSEG <= MAXSEG   ',NSEG <= MAXSEG   ,NSEG,MAXSEG 
          write(outu,'(a17,l5,2i10)')  'NATOM <= MAXA    ',NATOM <= MAXA    ,NATOM,MAXA 
          write(outu,'(a17,l5,2i10)')  'NBOND <= MAXB    ',NBOND <= MAXB    ,NBOND,MAXB 
          write(outu,'(a17,l5,2i10)')  'NTHETA <= MAXT   ',NTHETA <= MAXT   ,NTHETA,MAXT 
          write(outu,'(a17,l5,2i10)')  'NPHI <= MAXP     ',NPHI <= MAXP     ,NPHI,MAXP 
          write(outu,'(a17,l5,2i10)')  'NIMPHI <= MAXIMP ',NIMPHI <= MAXIMP ,NIMPHI,MAXIMP 
          write(outu,'(a17,l5,2i10)')  'NNB <= MAXNB     ',NNB <= MAXNB     ,NNB,MAXNB 
          write(outu,'(a17,l5,2i10)')  'NDON <= MAXPAD   ',NDON <= MAXPAD   ,NDON,MAXPAD 
          write(outu,'(a17,l5,2i10)')  'NACC <= MAXPAD   ',NACC <= MAXPAD   ,NACC,MAXPAD 
          write(outu,'(a17,l5,2i10)')  'NRES <= MAXRES   ',NRES <= MAXRES   ,NRES,MAXRES 
          write(outu,'(a17,l5,2i10)')  'NATC <= MAXATC   ',NATC <= MAXATC   ,NATC,MAXATC 
          write(outu,'(a17,l5,2i10)')  'NCB <= MAXCB     ',NCB <= MAXCB     ,NCB,MAXCB 
          write(outu,'(a17,l5,2i10)')  'NCT <= MAXCT     ',NCT <= MAXCT     ,NCT,MAXCT 
          write(outu,'(a17,l5,2i10)')  'NCP <= MAXCP     ',NCP <= MAXCP     ,NCP,MAXCP 
          write(outu,'(a17,l5,2i10)')  'NCI <= MAXCI     ',NCI <= MAXCI     ,NCI,MAXCI 
          write(outu,'(a17,l5,2i10)')  'NCH <= MAXCH     ',NCH <= MAXCH     ,NCH,MAXCH 
          write(outu,'(a17,l5,2i10)')  'NCN <= MAXCN     ',NCN <= MAXCN     ,NCN,MAXCN
#if KEY_CMAP==1
          write(outu,'(a17,l5,2i10)')'NCRTERM <= MAXCRT  ', &      
#endif
#if KEY_CMAP==1
               NCRTERM <= MAXCRT,NCRTERM,MAXCRT                    
#endif
       ENDIF
       CALL wrndie(-4,'<charmm_main.src>charmm','Bounds exceeded, try DIMENSION command.')
    ENDIF
70  FORMAT(//10X,'****  ERROR  **** A COUNTER VARIABLE HAS', &
         ' EXCEEDED ITS MAXIMUM ALLOWABLE VALUE.', &
         /16X,'EXECUTION WILL BE TERMINATED.', &
         /16X,'THE CURRENT COUNTER VARIABLES AND THEIR ', &
         'MAXIMUM ALLOWED VALUES ARE:', &
         /14X,'  NSEG = ',I8,' MAXSEG = ',I8,'  NATOM = ',I8, &
         '   MAXA = ',I8,'   NBOND = ',I8,'  MAXB = ',I8, &
         /14X,' NTHETA = ',I8,'   MAXT = ',I8,'   NPHI = ',I8, &
         '   MAXP = ',I8,'  NIMPHI = ',I8,' MAXIMP = ',I8, &
#if KEY_CMAP==1
         /14X,' NCRTERM = ',I8,' MAXCRT = ',I8,                 & 
#endif
         /14X,'    NNB = ',I8,'  MAXNB = ',I8,'   NDON = ',I8, &
         ' MAXPAD = ',I8,'    NACC = ',I8,' MAXPAD = ',I8, &
         /14X,'   NRES = ',I8,'  MAXRES = ',I8,'  NATC = ',I8, &
         ' MAXATC = ',I8,'     NCB = ',I8,'  MAXCB = ',I8, &
         /14X,'    NCT = ',I8,'   MAXCT = ',I8,'   NCP = ',I8, &
         '  MAXCP = ',I8,'     NCI = ',I8,'  MAXCI = ',I8, &
         /14X,'     NCH = ',I8,'   MAXCH = ',I8,'   NCN = ',I8, &
         '  MAXCN = ',I8)
  end function check_for_exceed_bounds

END PROGRAM CHARMM

recursive subroutine MAINCOMX(COMLYN,COMLEN,LUSED)
  !
  !  This is CHARMM's main command parser
  !
  !  Note: Miscellaneous commands (READ, WRITe,...) are parsed in MISCOM
  !        It should be called prior to calling this routine (if desired).
  !        Miscellaneous commands are available in subcommand parsers.
  !        The commands parsed here are not.    - BRB
  !
  !  Note: All commands in CHARMM are abbreviated to first 4 characters
  !        (except for graphics subcommands which are 3 characters).
  !
  use chm_kinds
  use chm_types
  use dimens_fcm
  use number
  use bases_fcm
  use contrl
  use coord
  use coordc
  use ctitla
  use hbondm
  use image
  use psf
  use modpsf
  use genpsf_m
  use stream
  use string
  use timerm
  use dcntrl_mod, only: dynopt
  use estats_mod,only: anal
#if KEY_ESTATS==1
  use estats_mod,only: estats        
#endif
  use cveloci_mod,only: cveloci
  use mltcanon,only:mltcanon_prs
  use mmfp,only:mmfp0
  use cadpac_mod,only:cadini
  use cstran_mod,only:cstran
#if KEY_QMMMSEMI==1
  use qmmmsemi,only: qmmm_startup         
#endif
#if KEY_EPMF==1
  use epmf, only: epmf_set                
#endif
#if KEY_PRIMO==1
  use primomodule, only: primo_set        
#endif
  use genborn,only: genborn_set
  use gbim, only: SetGBIM
  use gbmv, only:gbmv_set                 
  use gbsw, only: gbsw_set                


#if KEY_DENBIAS==1
  use denbias, only: denbias_set
#endif

  use phmd
  use gnn, only: gnncall
  use tamdmodule, only: tamd
#if KEY_MSCALE==1
  use mscalemod, only: mscale             
#endif
#if KEY_CSA==1 || KEY_DISTENE==1
  use csacommmod                          
#endif
#if KEY_REPDSTR2==1
  use repdstrmod2, only: repdstrmain
#else
  use repdstrmod, only: repdstrmain
#endif
#if KEY_OVERLAP==1
  use olapmod,only:olapcmd                
#endif
#if KEY_RPATH==1
  use EPATHMOD, ONLY: PATHPS              
#endif
#if KEY_FLUCQ==1
  use flucqm, only: fqinit                
#endif
#if KEY_TSALLIS==1
  use tsallis_module,only: iascale,qttsall 
#endif
#if KEY_FACTS==1
  use facts_module,only:fctini            
#endif
#if KEY_AFM==1
  use afm_module,only: afmini              
#endif
#if KEY_AXD==1
  use axd_module,only: axdini              
#endif
  use gukini_mod,only: gukini
  use gamess_fcm,only: qmused_qchem,qmused_g09,qmused_turbo,qmused_gamess
  use gamess_fcm,only: qmused_squantm,qmused_sccdftb,qmused_mndo97
  use gamess_fcm,only: qmused_nwchem
  use traj_mod,only:reatrj,wrttrj,trajio
  use block_fcm, only : block
#if KEY_SCPISM==1
  use scpismm,only:scpparse                
#endif
  use corman_mod,only:corcom
#if KEY_DMCONS==1
  use dmcons,only:dmcset                   
#endif
#if KEY_HQBM==1
  use hqbmm, only: hqbmini                 
#endif
  use pull_mod, only: pull_setup
  use rmsdyn_mod,only:rmsdyn
#if KEY_GENETIC==1
  use galgor,only:genetic_alg             
#endif
#if KEY_PATHINT==1
  use mpathint, only: pint_init            
#endif
#if KEY_POLAR==1
  use polarm, only: polar0                 
#endif
#if KEY_RGYCONS==1
  use rgym, only: rgyset                   
#endif
  use rush_mod, only: rush                 
#if KEY_RXNCOR==1
  use lup, only: lupopt                    
#endif
  use grid_dock, only: gridset
#if KEY_ASPENER==1
  use eef1_mod, only: eef1                 
#endif
#if KEY_PBEQ==1
  use pbeq, only: pbeq0                    
#endif
#if KEY_TORQUE==1
  use torque, only: torque_parse           
#endif
  use consph, only: getphresidues          
#if KEY_EDS==1
  use edsmod, only: process_eds            
#endif
#if KEY_CHEQ==1
  use cheq,only: cheqprep 
#endif
  use scalar_module, only: scalar
#if KEY_PROTO==1
  use proto_mod,only:proto     
#endif
  use intcor_module,only: intcor
  use rdfsol_mod,only:rdfsol
  use rxcons,only:rxconsps
#if KEY_CORSOL==1
  use corsol_mod,only:corsol     
#endif
  use shell,only:shlini
  use correl_mod,only: correl
  use nmrm,only:nmr
  use tmd,only: tmdinit
  use dims,only: dimsinit
  use travelmain,only: trek
  use nbndcc_util
#if KEY_RXNCOR==1
  use rxpath,only:pathc            
#endif
#if KEY_RXNCOR==1
  use rxdefs,only:rxpars           
#endif
  use replica_mod,only:replica
#if KEY_DYNVV2==1
  use tpvv2,only: tpcontrol         
#endif
  use resdist,only:redset
  use machutil,only:jobdat
#if KEY_CHARMMRATE==1
  use charmmrate_mod,only:charmmrate   
#endif
#if KEY_MC==1
  use mcc, only: mccall               
#endif
#if KEY_MC==1
  use mcmoveio, only: moverd, movewr  
#endif
#if KEY_MC==1
  use mcmoveln, only: moveln          
#endif
#if KEY_MC==1
  use mcmvad                          
#endif
# if KEY_ACTBOND==1
  use eintern_fast, only: bactiv
#endif
  use lonepr
  use mtp_fcm, only: mtp
  use mtpl_fcm
  use ediff, only: ediffop
  use enbond_mod
  use eutil
  use fitchg
#if KEY_MCMA==1
  use mcmamod                         
#endif
#if KEY_OPENMM==1
  use omm_ctrl, only : omm_command, omm_system_changed
#endif
  use molvco_m
  use pert_mod
  use qub_m
  use quantm,only:qmused_quantum
  use testch_m
  use select
  use shapes
  use tbmts
  use varcutm
  use vibran_m
#if KEY_DOMDEC==1
  use domdec,only:domdec_com          
#endif
  use linkatom
  use drude
#if KEY_PARALLEL==1
  use parallel,only:mynod,mynodg
  use paral4,only: setcpustruc, test_po
#endif
#if KEY_ENSEMBLE==1
  use ensemble,only:nensem,old_mynod,whoiam,ensprint     
#endif
#if KEY_ABPO==1
  use abpo_ltm  
#endif
#if KEY_ABPO==1
  use abpo,only:abpo_cntrl 
#endif
  use prssre
#if KEY_GAMUS==1
  use gamusmodule,only:gamusinit 
#endif
#if KEY_MULTICOM==1 /*  VO : stringm v */
  use multicom, only: multicom_main
  use ifstack                          ! VO:  parsing conditionals in parallel
#endif
  use zmodule
#if KEY_PNM==1 /* VO: Plasic Network Model */
  use pnm, only : pnm_main
#endif
#if KEY_LARMORD==1
use larmord, only : setup_larmord
#endif
#if KEY_SSNMR==1
  use ssnmr
#endif
#if KEY_RDC==1
  use rdc, only : rdcset
#endif
use gopair, only : GoPair_setup
use freeene_calc, only: frencalc, qmfix, mkdummy
use keywords, only: print_keys
use cstuff, only: unbuffer_stdout
  
  implicit none
  character(len=*) comlyn
  integer comlen
  logical lused

  !     The following are local variables.

#if KEY_TSALLIS==1
  integer   err  
#endif
  integer   freeat
  integer  :: i,ierror
  integer   idm1,idm2,idm3
  integer   istart,itemp
  integer   natiml
  integer   nsgiml

  real(chm_real)    xtmp

  logical    lcomp
  logical   qform,qopen,qwrite

  character(len=4)  wrd
  character(len=6)  access
  character(len=12) form

  integer  tlen
  integer,parameter :: tmax=50
  character(len=tmax) tname

  integer :: icycle=0
  real(chm_real),dimension(1) :: ddv_dummy

  lused=.true.

#if KEY_MULTICOM==1 /*  VO stringm : conditional evaluation in parallel */
  if (.not.peek_if()) then
   comlen=0
   return ! this if_then_else_endif block disabled on this node
  endif
#endif /* VO */

  !     Shorten word to four characters or less and pad with blanks
  !     main conditional for processing commands

  wrd=nexta4(comlyn,comlen)

  cmds: select case(wrd)
  case('    ') cmds
  case('ANAL') cmds
     call anal(comlyn,comlen)
  case('ATLI') cmds
     call trime(comlyn,comlen)
     call copyst(altcom,mxalsz,altlen,comlyn,comlen)
     comlen=0
  case('AUTO') cmds
     call autogen(comlyn,comlen)
  case('BLOC') cmds
     call block(comlyn,comlen)
#if KEY_TSALLIS==1
  case('TSAL') cmds
     wrd = nexta4(comlyn, comlen)
     if (wrd  ==  'TORS') THEN
        qttsall = .true.
        allocate(iascale(natom), stat=err)
        if (err /= 0) then
           call wrndie(-4,'<charmm_main>', &
                'Abort: unable to allocate IASCALE')
        endif
        wrd = nexta4(comlyn,comlen)
        if (wrd  ==  'SELE') THEN
           ! Assign selection
           call joinwd(comlyn,mxcmsz,comlen,wrd,4)
           call selcta(comlyn,comlen,iascale, &
                x,y,z,wmain,.true.)
        endif
     endif
#endif 
  case('NOSE') cmds
     call nosect(comlyn,comlen)
  case('MTS ') cmds
     call mts(comlyn,comlen)
  case('BOUN') cmds
     call bound(comlyn,comlen)
  case('BUIL') cmds
     call intcor(comlyn,comlen)
  case('CADP') cmds
     call cadini(comlyn,comlen)
  case('CONS') cmds
     call cstran(comlyn,comlen)
  case('COOR') cmds
     call corcom(comlyn,comlen)
  case('CORR') cmds
     call correl(0,ddv_dummy)
  case('CRYS') cmds
     call crystl
  case('DIME','RESI' ) cmds
     call wrndie(-1,'<CHARMM>', &
          'Use DIMEnsion change only as first CHARMM command')
#if KEY_DIMS==1
  case('DIMS') cmds
     call dimsinit(comlyn,comlen)
#endif 
  case('ENSE') cmds    
     call enscmd(comlyn,comlen)

#if KEY_CHARMMRATE==1
  case('POLY','RATE') cmds
     call charmmrate(comlyn,comlen)
#endif 
  case('NBAC') cmds
     call nbactv(comlyn,comlen)
#if KEY_ACTBOND==1
  case('BACT') cmds
     call bactiv(comlyn,comlen)
#endif 
  case('MKCL') cmds
     call mkclust(comlyn,comlen)
#if KEY_ESTATS==1
  case('ESTA') cmds
     call estats(comlyn,comlen)
#endif /* */
     !**********************************************************
  case('MKPR') cmds
     call mkpres0
  case('CVEL') cmds
     if(prnlev >= 2) write(outu, &
          '(/,1X,''CVELOCI>  SETING UP VECTOR DIRECTIVES'',/)')
     call cveloci(comlyn,comlen)
  case('DEAD') cmds
     !---- Procedure PROCESS-DEADLINE-COMMAND
     cpulim=gtrmf(comlyn,comlen,'CPU',zero)*60.d0
     if(cpulim > 0.0) call jobdat(xtmp,cpuini,idm1,idm2,idm3)
     xtmp=gtrmf(comlyn,comlen,'CLOC',minone)
     if (xtmp >= 0.0) then
        deadhr=xtmp
        deadmn=100.*(xtmp-deadhr)
        pasmid=0
        call chklim
        if(atlim) pasmid=pasmid+1
     else
        deadhr=-1
     endif
     !---- End Procedure PROCESS-DEADLINE-COMMAND
  case('DELE') cmds
     !     call from Subroutine mmff_io.
     call deltic(comlyn,comlen)
  case('DONO','ACCE') cmds
     call edhbatoms(comlyn,comlen,wrd)
#if KEY_TMD==1
  case('TMDI') cmds
     call tmdinit(comlyn,comlen)
#endif 
  case('DYNA') cmds
#if KEY_ABPO==1
     if (q_abpo) then
        call abpo_cntrl(comlyn,comlen)
     else
        call dynopt(comlyn,comlen)
     endif
#else /**/
     call dynopt(comlyn,comlen)
#endif /* ABPO*/
  case('EDIF') cmds
     call ediffop(comlyn,comlen,icycle)
  case('ENER') cmds
     call gete0('ENER', comlyn, comlen)
#if KEY_CHEQ==1
  case('FQBA') cmds
     call nosefq(comlyn,comlen)
  case('CHEQ') cmds
     call cheqprep(comlyn,comlen)
#endif 
  case('FLUC') cmds
#if KEY_FLUCQ==1
     call fqinit(comlyn,comlen)
#else /**/
     CALL WRNDIE(-1,'<CHARMM>','FLUCQ code is not compiled.')
#endif 
  case('FOUR') cmds
     call parse4d(comlyn, comlen)
#if KEY_FACTS==1
  case('FACT') cmds
     call fctini(comlyn,comlen)
#endif 
  case('GBOR') cmds
     call genborn_set(comlyn, comlen)
  case('GBIM') cmds
     call setgbim(comlyn, comlen)
  case('GBMV') cmds
     call gbmv_set(comlyn, comlen)
  case('GBSW') cmds
     call gbsw_set(comlyn, comlen)
  case('GOPA') cmds
     call GoPair_setup(comlyn, comlen)

#if KEY_DENBIAS==1
  case('DBIA') cmds
     call denbias_set(comlyn, comlen)
#endif

#if KEY_EPMF==1
  case('EPMF') cmds
     call epmf_set(comlyn,comlen)
#endif 
#if KEY_PRIMO==1
   case('PRIM') cmds
     call primo_set(comlyn,comlen)
#endif 
  case('GENE') cmds
     call genpsf(comlyn, comlen, istart)
  case('GETE') cmds
     call gete0('GETE', comlyn, comlen)
  case('GAME') cmds
     qmused_gamess=.true.
     call gukini(comlyn,comlen)
  case('NWCH') cmds
     qmused_nwchem=.true.
     call gukini(comlyn,comlen)
  case('QCHE') cmds
     qmused_qchem=.true.
     call gukini(comlyn,comlen)
! Guanhua_QC_UW1111 / JZ_UW12
  case('GAUS') cmds
     qmused_g09=.true.
     call gukini(comlyn,comlen)
  case('QTUR') cmds
     qmused_turbo=.true.
     call gukini(comlyn,comlen)

#if KEY_NOGRAPHICS==0
  case('GRAP') cmds
     call graphx
#endif 
  case('HBON') cmds
     if(ihbfrq == 0) ihbfrq=999
     call update(comlyn,comlen,x,y,z,wmain,.true., &
          .false.,.false.,.true.,.false.,0,0,0,0,0,0,0)
     if(ihbfrq == 999) ihbfrq=0
  case('HBTR') cmds
     call hbtrim
  case('HBUI') cmds
     call hbuild(comlyn,comlen)
  case('IC  ') cmds
     call intcor(comlyn,comlen)
  case('INQU') cmds
     !---- Procedure PROCESS-INQUIRE-FILE-COMMAND
     !     inquire all open files
     if(prnlev >= 2) write(outu,'(a)') &
          ' CHARMM: list of open files:'
     do i=1,99
        call vinqre('UNIT',tname,tmax,tlen,qopen,qform,qwrite,i)
        if (qopen) then
           if (qform) then
              form=' formatted'
           else
              form=' unformatted'
           endif
           if (qwrite) then
              access=' write'
           else
              access=' read'
           endif
           if(prnlev >= 2) write(outu,'(9x,i3,1x,3a)') &
                i,tname(1:tlen),access,form
        endif
     enddo
     !---- End Procedure PROCESS-INQUIRE-FILE-COMMAND
  case('INTE') cmds
     call inteop(comlyn,comlen,icycle)
  case('JOIN') cmds
     call joinsg(comlyn,comlen)
  case('LONE') cmds
     call loneprs(comlyn,comlen)
#if KEY_RXNCONS==1
  case('RCON') cmds
     call rxconsps(comlyn, comlen)
#endif 
#if KEY_RXNCOR==1
  case('LUPO') cmds
     call lupopt(comlyn,comlen)
#endif 
#if KEY_CSA==1 || KEY_DISTENE==1
  case('MAST') cmds
     call masterdstr(comlyn,comlen)
  case('ETRA') cmds
     call etraj(comlyn,comlen)
  case('RECE') cmds
     call calcrece(comlyn,comlen)
  case('TRAN') cmds
     call calctran(comlyn,comlen)
#if KEY_CSA==1
  case('CSA ') cmds
     call csacntrl(comlyn,comlen)
#endif 
#endif 
  case('MC  ') cmds
#if KEY_MC==1 /*mc1*/
     call mccall(comlyn,comlen)
#else /*  (mc1)*/
     CALL WRNDIE(-1,'<CHARMM>','MC code is not compiled.')
#endif /* (mc1)*/
  case('GNN ') cmds
     call gnncall(comlyn, comlen)
  case('MERG') cmds
     CALL TMERGE(COMLYN,COMLEN,TITLEB,NTITLB)
  case('MINI') cmds
     CALL MINMIZ(COMLYN,COMLEN)
  case('MLTE') cmds           !bt_080627
     CALL MULTE0(COMLYN,COMLEN)
  case('MMFP') cmds
#if KEY_OPENMM==1
     call omm_system_changed()
#endif
     CALL MMFP0
     !----namkh 02/03/03
  case('MNDO') cmds
#if KEY_MNDO97==1
     qmused_mndo97=.true.
     CALL MNDINI(COMLYN,COMLEN)
#endif 
  case('MOLV') cmds
     CALL MOLVCO(COMLYN,COMLEN)
  case('MONI') cmds
     !---- Procedure PROCESS-MONITOR-COMMANDS
     WRD=NEXTA4(COMLYN,COMLEN)
     IF (WRD == 'DIHE') THEN
        CALL MONDIH(COMLYN,COMLEN)
     ELSE
        CALL WRNDIE(0,'<CHARMM>','UNRECOGNIZED MONITOR OPTION')
     ENDIF
     !---- End Procedure PROCESS-MONITOR-COMMANDS
  case('MOVE') cmds
#if KEY_MC==1 /*mc2*/
     ! ----- Procedure PROCESS-MOVE-COMMANDS
     wrd=nexta4(comlyn,comlen)
     movecmds: select case(wrd)
     case('ADD ') movecmds
        call movead(comlyn,comlen)
     case('LINK') movecmds
        call moveln(comlyn,comlen)
     case('DELE') movecmds
        call movedl(comlyn,comlen)
     case('EDIT') movecmds
        call moveed(comlyn,comlen)
     case('READ') movecmds
        call moverd(comlyn,comlen)
     case('WRIT')  movecmds
        call movewr(comlyn,comlen)
     case default movecmds
        CALL WRNDIE(0,'<CHARMM>','UNRECOGNIZED MOVE OPTION')
     end select movecmds
     ! ----- End Procedure PROCESS-MOVE-COMMANDS
#else /*   (mc2)*/
     CALL WRNDIE(-1,'<CHARMM>','MC code is not compiled.')
#endif /*  (mc2)*/
  case('MLTC') cmds
     call mltcanon_prs(comlyn,comlen)
#if KEY_MSCALE==1
  case('MSCA') cmds
     call mscale(0,comlyn,comlen)
  case('SERV') cmds
     call mscale(1,comlyn,comlen)
#endif 
  case('NBON') cmds
     if(inbfrq == 0) inbfrq=999
     call update(comlyn,comlen,x,y,z,wmain,.true., &
          .false.,.true.,.false.,.true.,0,0,0,0,0,0,0)
     if(inbfrq == 999) inbfrq=0
  case('NOE ') cmds
     call noeset
  case('RESD') cmds
     call redset
  case('NMR ') cmds
     call nmr
  case('TAMD') cmds
     call tamd
  case('ZMAT') cmds
     call zmtrx
#if KEY_MCMA==1 /*mcma*/
  case('MCMA') cmds
     call mcma
#endif /* (mcma)*/
  case('MMQM') cmds
#if KEY_NOMISC==1 /*mmqm*/
     CALL WRNDIE(-1,'<CHARMM>','MMQM code is not compiled.')
#else /* (mmqm)*/
     call g94ini
#endif /* (mmqm)*/
  case('MTP ') cmds
     !  atomic multipole moments mtp module
     call mtp
  case('MTPL') cmds
     !  atomic multipole moments mtpl module--uses local axis systems for
     !  arbitrarily large molecules.
#if KEY_MTPL==1
     call mtpl
#else
     CALL WRNDIE(-1,'<CHARMM>','MTPL code is not compiled.')
#endif /* mtpl */
#if KEY_OVERLAP==1
  case('OLAP') cmds
     call olapcmd(comlyn,comlen)
#endif 
#if KEY_OPENMM==1
  case('OMM ') cmds
     call omm_command(comlyn, comlen)
#endif 
#if KEY_PARCMD==1
  case('PARA') cmds
     call parcmd(comlyn,comlen)
#endif 
  case('PATC') cmds
     call patch(comlyn,comlen)

! phmd dependent on gbsw or gbsv
#if KEY_PHMD==1
  case('PHMD') cmds
     call startphmd(comlyn,comlen)
  case('PHTE') cmds
     call dophmdtest(comlyn, comlen)
#endif 
#if KEY_PBEQ==1
  case('PBEQ') cmds
     call pbeq0
#endif 
#if KEY_GRID==1
  case('GRID') cmds
     call gridset(comlyn, comlen)
#endif
  case('PERT') cmds
     call perts(comlyn,comlen)
#if KEY_POLAR==1
  case('POLA') cmds
     call polar0
#endif 
#if KEY_PATHINT==1
  case('PINT') cmds
     call pint_init(comlyn,comlen)
#endif 
#if KEY_PNM==1
  case('PNM') cmds
     call pnm_main(comlyn,comlen)
#endif 
  case('PRES') cmds
     call getprs(comlyn,comlen)
  case('PRIN') cmds
     call mainio(wrd)
  case('PULL') cmds
     call pull_setup(comlyn,comlen)
  case('QUAN') cmds
! Because of this we cannot compile SQUANTM & QUANTUM at the same time
! fixing this later ...
#if KEY_SQUANTM==1
     qmused_squantm=.true.
     call sqmini(comlyn,comlen)
#else /**/
     qmused_quantum=.true.
     call qmdefn(comlyn,comlen)
#endif 
#if KEY_QMMMSEMI==1
  case('IFQN') cmds
     call qmmm_startup(comlyn,comlen,natom,x,y,z)
#endif 
  case('READ') cmds
     call mainio(wrd)
#if KEY_GENETIC==1
     !---- 12-Jan-98, CLBIII
  case('GALG') cmds
     IF (IndxA(comLyn,comLen,'SETU') > 0) then
        CALL Genetic_Alg(.true., .false.)
     ELSEIF (IndxA(comLyn,comLen,'EVOL') > 0) then
        CALL Genetic_Alg(.false.,.true.)
     endif
#endif 
  case('RENA') cmds
     !---- Procedure PROCESS-RENAME-COMMAND
     call crename(comlyn,comlen)
  case('REPL') cmds
     call replica
  case('REPD') cmds
     call repdstrmain
  case('RESE') cmds
     !---- Procedure PROCESS-RESET-COMMAND
     call iniall
     if(prnlev >= 2) write(outu, &
          '(/,1X,''**** CHARMM DATA STRUCTURES RESET ****'',/)')
     !---- End Procedure PROCESS-RESET-COMMAND
     !-----cb3 added prefx control for compilation
#if KEY_HQBM==1
  case('HQBM') cmds  ! 02-Jul-1997
     call hqbmini
#endif 
#if KEY_AFM==1
  case('AFM') cmds   ! 20-Jun-2003
     call afmini
#endif 
#if KEY_AXD==1
  case('AXD') cmds
     call axdini(natom)
#endif 
  case('RISM') cmds
     call rismcmd
#if KEY_RPATH==1
  case('RPAT') cmds
     call pathps(comlyn,comlen)
#endif 
#if KEY_RXNCOR==1
  case('PATH') cmds
     call pathc(x,y,z,xcomp,ycomp,zcomp,wmain,comlyn,comlen)
  case('RXNC') cmds
     call rxpars
#endif 
  case('SBOU') cmds
     !---- Procedure PROCESS-SOLVENT-BOUNDARY-COMMANDS
#if KEY_NOMISC==1
     CALL WRNDIE(-1,'<CHARMM>','SBOUND code is not compiled.')
#else /**/
     wrd=nexta4(comlyn,comlen)
     solvboun: select case(wrd)
     case('POTE') solvboun
        call sbintg
     case('SET ') solvboun
        call bndset
     case('READ') solvboun
        call sbread
     end select solvboun
#endif 
     !---- End Procedure PROCESS-SOLVENT-BOUNDARY-COMMANDS
  case('SCAL') cmds
     call scalar
#if KEY_SCCDFTB==1
     !---- SCCDFTB CODE (QCME)
  case('SCCD') cmds
     qmused_sccdftb=.true.
     call scctbini(comlyn,comlen)
#endif 
  case('SHAK') cmds
     call shkcom(comlyn, comlen)
#if KEY_SHAPES==1 /*shpcom*/
  case('SHAP') cmds
     call shpcom(comlyn, comlen)
#endif /* (shpcom)*/
  case('SHOW') cmds
     comlen=0
  case('TEST') cmds
     call testch(comlyn,comlen)
#if KEY_PARALLEL==1
  case('TSPO') cmds
     call test_po
  case('CPUS') cmds
     call setcpustruc
#endif

  case('TRAJ') cmds
     !---- Procedure PROCESS-TRAJECTORY-COMMAND
     if (indxa(comlyn,comlen,'READ') > 0) THEN
        if (indxa(comlyn,comlen,'COMP') <= 0) THEN
           call reatrj(x,y,z)
        else
           call reatrj(xcomp,ycomp,zcomp)
        endif
     else if (indxa(comlyn,comlen,'WRIT') > 0) THEN
        if (indxa(comlyn,comlen,'COMP') <= 0) THEN
           call wrttrj(x,y,z)
        else
           call wrttrj(xcomp,ycomp,zcomp)
        endif
     else
        call trajio
     endif
     !---- End Procedure PROCESS-TRAJECTORY-COMMAND
#if KEY_ADUMB==1
  case('UMBR') cmds
     call umban(comlyn,comlen)
#endif 
#if KEY_GAMUS==1
  case('GAMU','GUMB')
     CALL GAMUSINIT(COMLYN,COMLEN)
#endif 
     !mf switch stdout to unbuffered I/O
  case('UNBU') cmds
#if KEY_UNIX==1
     ierror = unbuffer_stdout()
     if (ierror .ne. 0) then
        call wrndie(0, '<CHARMM>', 'UNBU command failed')
     end if
#endif
  case('TREK','TRAV') cmds
     call trek(comlyn,comlen)
  case('TSM ') cmds
     call tsms
  case('WRIT') cmds
     call mainio(wrd)
  case('UPDA') cmds
     !---- Procedure PROCESS-UPDATE-COMMAND
     !
     !     UPDATE NBOND, HBOND, AND CODES LISTS
     !
     lcomp=(indxa(comlyn,comlen,'COMP') > 0)
     if (lcomp) then
        call update(comlyn,comlen,xcomp,ycomp,zcomp,wcomp,.true., &
             .true.,.true.,.true.,.true.,0,0,0,0,0,0,0)
     else
        call update(comlyn,comlen,x,y,z,wmain,.true., &
             .true.,.true.,.true.,.true.,0,0,0,0,0,0,0)
     endif
     !---- End Procedure PROCESS-UPDATE-COMMAND
  case('VARC') cmds
     call setup_varcut(comlyn,comlen,natom)
  case('VIBR') cmds
     call vibopt(comlyn,comlen)
  case('WHAM') cmds
     call wham0
  case('IMAG') cmds
     !---- Procedure PROCESS-IMAGE-SPECIFY-COMMAND
     !
     !     SPECIFY IMAGE CENTERING OPTIONS
     if(ntrans <= 0) then
        CALL WRNDIE(-2,'<CHARMM>','IMAGES NEED TO BE PRESENT')
        return
     endif
     call imspec(comlyn,comlen,bimag%imcenf, &
          limcen,imxcen,imycen,imzcen,natom,x,y,z,wmain)
     !---- End Procedure PROCESS-IMAGE-SPECIFY-COMMAND
  case('IMPA') cmds
     call impatc(comlyn,comlen)
#if KEY_NOGRAPHICS==0 && KEY_NODISPLAY==0
  case('DRAW') cmds
     !       For XDISPLAY
     call drawit(-1,x,y,z,xcomp,ycomp,zcomp)
#endif 
#if KEY_NOMISC==0
  case('BARI') cmds
     !---- Procedure PROCESS-BARRIER-COMMAND
     IF(WRNLEV >= 2) WRITE(OUTU,'(A)') &
          ' CHARMM> BARRier command is no longer supported.'
     IF(WRNLEV >= 2) WRITE(OUTU,'(A)') &
          '         Use the GRADient option of MINImize commands.'
     CALL WRNDIE(-5,'<CHARMM>','Unsupported command: '//WRD)
     !---- End Procedure PROCESS-BARRIER-COMMAND
  case('RMSD') cmds
     !---- Procedure PROCESS-RMSDYNAMICS-COMMAND
     natiml=natom
     nsgiml=nseg
     if(indxa(comlyn,comlen,'IMAG') > 0) THEN
        natiml=natomt
        nsgiml=nsegt
     endif

     call rmsdyn('main',x,y,z,wmain,xcomp,ycomp,zcomp,wcomp,natiml,amass)

     !---- End Procedure PROCESS-RMSDYNAMICS-COMMAND
#endif 
#if KEY_QUANTUM==1 || KEY_GAMESS==1 || KEY_GAMESSUK==1 || KEY_CADPAC==1 || \
    KEY_SCCDFTB==1 || KEY_QCHEM==1 || KEY_MNDO97==1 || KEY_SQUANTM==1 || \
    KEY_QTURBO==1 || KEY_G09==1
  case('ADDL') cmds
     call addlnat(outu)
  case('RELL') cmds
     call rellnat(outu)
#endif 
#if KEY_QUANTUM==1
  case('MULL') cmds
     call mullik(comlyn,comlen)
#endif 
#if KEY_QUANTUM==1 || KEY_SCCDFTB==1 || KEY_SQUANTM==1 || KEY_QCHEM==1 || KEY_GAMESSUK==1 || KEY_QTURBO==1
  case('QUB ') cmds
     call qub(comlyn,comlen)
#endif 
#if KEY_SQUANTM==1
  case('SQUA') cmds
     qmused_squantm=.true.
     call sqmini(comlyn,comlen)
#endif 
#if KEY_DMCONS==1
  case('DMCO') cmds
     call dmcset
#endif 
#if KEY_RGYCONS==1
  case('RGYR') cmds
     call rgyset
#endif 
  case('ETEN') cmds
     call etenset
  case('ETSR') cmds
     call etsrset
#if KEY_ASPENER==1
  case('EEF1') cmds
     call eef1
#endif 
#if KEY_SCPISM==1
  case('SCPI') cmds
     call scpparse
#endif 
#if KEY_SASAE==1
  case('SASA') cmds
     call sasini(comlyn,comlen)
#endif 
  case('PREF') cmds
     call print_keys()
#if KEY_EMAP==1
  case('EMAP') cmds
     call emapopt(comlyn, comlen)
#endif 
     ! BEGIN DRUDE (B. Roux and G. Lamoureux)
  case('DRUD') cmds
     IF (Indx(comlyn,comlen,'L_WALL',6) > 0) then
        ! Lei Huang, impose hard wall constraint on drude bond length
        L_WALL = GTRMF(COMLYN,COMLEN,'L_WALL',0.2)
     
        IF(L_WALL .gt. 0.0) THEN
          QHARDWALL = .true.
#if KEY_PARALLEL==1
          IF(mynod .eq. 0) then
#endif 
            WRITE(outu,*)   &
            'Hard wall constraint on drude bond length is turned ON. L_WALL = ',L_WALL
#if KEY_PARALLEL==1
          endif
#endif
        ELSE
          QHARDWALL = .false.
#if KEY_PARALLEL==1
          IF(mynod .eq. 0) then
#endif
            WRITE(outu,*)   &
            'Hard wall constraint on drude bond length is turned OFF.'
#if KEY_PARALLEL==1
          endif
#endif
        ENDIF
     ELSE
        call drude0(comlyn,comlen)
     ENDIF

  case('TPCO','NOS2') cmds
#if KEY_DYNVV2==1
     call tpcontrol(comlyn,comlen)                     
#endif
#if KEY_DYNVV2==0
     CALL WRNDIE(0,'<CHARMM>','DYNA VV2 not compiled') 
#endif

  case('FITC') cmds
     call fitcharge

  case('FITP') cmds
     call fitparam

#if KEY_RDFSOL==1
  case('RDFS') cmds
     !       new solvent radial distribution function module
     call rdfsol
#endif 
#if KEY_SHELL==1
  case('SHEL') cmds
     !       shell decomposition module
     call shlini
#endif 

  case('ZMOD') cmds
     call zerom

#if KEY_PROTO==1 /*proto*/
  case('PROT') cmds
     !       prototype functions
     call proto
#endif /*   (proto)*/

  case('PIPF','PFBA') cmds
#if KEY_PIPF==1
     if (wrd == 'PIPF') then
        call pfprep(comlyn,comlen)
     else if (wrd == 'PFBA') then
        call nosepf(comlyn,comlen)
     end if
#else /**/
     CALL WRNDIE(-1,'<CHARMM>','PIPF code is not compiled.')
#endif 
  case('RUSH') cmds
     call rush(comlyn,comlen)
#if KEY_CORSOL==1 /*corsol*/
  case('CORS') cmds
     !        solvent correlation functions
     call corsol
#endif /*     (corsol)*/
#if KEY_TORQUE==1 /*torque*/
  case('TORQ') cmds
     !       process TORQue command
     call torque_parse
#endif /* (torque)*/
  case('CNSP') cmds
     ! process CNSPh command
     call getphresidues(comlyn,comlen)
#if KEY_LARMORD==1
  case('LARM') cmds
     ! set-up calculation of chemical shifts with LarmorD
     call setup_larmord()
#endif
#if KEY_EDS==1
  case('EDS') cmds
     ! process EDS command
     call process_eds(comlyn,comlen)
#endif 
#if KEY_DOMDEC==1
  case('DOMD') cmds
     call domdec_com(comlyn, comlen)
#endif
#if KEY_STRINGM==1 /*  VO stringm v */
  case('STRI') cmds
    CALL SM_MAIN(COMLYN,COMLEN) ! string method parser
#endif
#if KEY_MULTICOM==1 /*  VO stringm communication v */
  case('MCOM') cmds
    CALL MULTICOM_MAIN(COMLYN,COMLEN) ! manipulate communicators for stringm
#endif /* VO ^ */
  case('FREN') cmds
    CALL FRENCALC(COMLYN,COMLEN)
  case('QMFI') cmds
    CALL QMFIX(COMLYN,COMLEN)
  case('MKDU') cmds
    CALL MKDUMMY(COMLYN,COMLEN) ! Write a new parameter file that turns selected atoms to dummy atoms

#if KEY_RDC==1
   case('RDC ') cmds
      call rdcset
#endif
#if KEY_SSNMR==1
  case('CCS ') cmds
    CALL CSSET
#endif
     !
  case default cmds
     lused=.false.
#if KEY_ENSEMBLE==1
     if(nensem > 1 ) write(outu,'(3(a,i3),2a)') &     
#endif
#if KEY_ENSEMBLE==1
          ">>> Ensemble ",whoiam," Node ",mynod," Worldnod ",old_mynod, &    
#endif
#if KEY_ENSEMBLE==1
          " ---- cmd problem = ",comlyn(1:comlen)       
#endif
     call wrndie(0,'<CHARMM>','Unrecognized command: '//WRD)
  end select cmds

  return

end subroutine maincomx

