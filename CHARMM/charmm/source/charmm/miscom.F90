SUBROUTINE MISCOM(COMLYN,MXCMS2,COMLEN,LUSED)
  !-----------------------------------------------------------------------
  !     This routine processes many of the miscellaneous commands
  !     it may be called from any routine that processes commands.
  !
  !      Note: length of all character strings used is 4.
  !
  !     Control flow commands by David States
  !     Overhauled by Bernard R. Brooks and Axel Brunger, 1-MAR-83
  !     Further modifications by Leo Caves, 17-JAN-94 (parameter table)
  !


  ! 12/12/07 SL: use module cross to allocate memory for RMD arrays
#if KEY_RMD==1
  use cross, only: CROSSINIT, QCROS     
#endif
#if KEY_REPDSTR==1
  use repdstrmod                        
#endif
#if KEY_REPDSTR2==1
  use repdstrmod2                        
#endif

  use chm_kinds
  use dimens_fcm
  use number
  use cmdpar
  use ctitla
  use eutil
#if KEY_BLOCK==1
  use lambdam       
#endif
  use fast
#if KEY_LOOKUP==1
  use lookup,only:wwsetup        
#endif
  use machdep
  use param_store, only: write_real_params, write_int_params, write_str_params
  use select
  use stream
  use string
  use timerm
  use parallel
  use usermod,only: usersb
  use quick,only:quicka
  use aidxmod,only:aidx
  use repdstr
  use calc,only:calculator
  use clcg_mod,only: randspec, irandom
  use machutil,only:timre,timrb,daytim,wcpu,csystem
#if KEY_VALBOND==1
  use valbond, only: vbcomm   
#endif
  use mtp_fcm
  use mtpl_fcm
#if KEY_MMPT==1
  use mmpt_fcm,only:mmptinit                              
#endif
#if KEY_MRMD==1
  use mrmd_fcm,only:mrmd_init
#endif
#if KEY_MULTICOM==1
  use ifstack    ! VO: if/else/endif conditionals in srtingm communication
  use multicom_aux, only: MPI_COMM_PARSER, ME_PARSER, SIZE_PARSER
#endif

#if KEY_UNIX==1
  use cstuff, only: setenv
  use, intrinsic :: iso_c_binding, only: C_NULL_CHAR
#endif /* UNIX */
  
  implicit none

  character(len=*) comlyn
  integer   mxcms2,comlen
  logical   lused

  ! following temp. string sizes taken from cmdpar.f90
  character(len=mxtlen) toktmp
  character(len=mxvlen) valtmp
  ! a variable name LENVAL is declared downstream!
  integer lentok
  logical   eof, done, ok
#if KEY_MULTICOM==1 /*   VO : conditionals  */
  logical :: ok2                
#endif
  character(len=4)   wrd,junk
  integer lenvar,lenval,ipt
  integer,parameter :: maxenv=120
  character(len=maxenv) envval, envvar
  integer wd2len,wd3len
  integer,parameter :: wdmax=50
  character(len=wdmax) wrd2,wrd3
  integer   i, j, ipar, iunit, idummy, irwait
  integer   k,l,m,ii,jj,kk
  real(chm_real)    value, temp, increm
  ! for numerical parameter modifying operations
  integer   parop
  integer,parameter :: incr=1,decr=2,mult=3,divi=4,expo=5

  integer   lablen
  integer,parameter :: labmax=80
  character(len=labmax) lablyn
  integer ilevel

  character(len=1) :: sdblq='"'
  integer idum(1)

  lused=.true.
  done=.false.
  eof=.false.
  ilevel=0

  loop101: do while(.not. done)
     done=.true.
     parOp = -1 
     wrd=nexta4(comlyn,comlen)

     !     main conditional for processing commands
     !
#if KEY_MULTICOM==1 /* (mcom)  VO : conditionals for stringm */
     !     if/else/endif followed by 'ELSE IF (peek_if)' MUST be processed first; 
     !     this is to permit the evaluation of multi-line nested conditionals in parallel
     !======================================================================================
      IF (WRD.EQ.'ELSE') THEN
        OK=peek_if() ! 'peek_if()' flags whether current if block (or main) is executed on this node
        OK2=pop_if() ! pop the 'if' stack
        if (.not.OK2) then ! errors popping stack?
          call wrndie(-2,' MISCOM>','UNEXPECTED ELSE STATEMENT.')
        else
        OK2=peek_if()
         IF(PRNLEV.GE.2.and.(OK.and.OK2)) & ! if previous statements executed, then the following cannot!
&             WRITE(OUTU,'(2A)') ' Skip to ELSE or ENDIF'
         call push_if(.not.OK.and.OK2) ! the opposite of the previous condition is true
        endif
     !
      ELSE IF (WRD.EQ.'ENDI') THEN
        OK=pop_if() ! pop the 'if' stack
        if (.not.OK) then 
          call wrndie(0,' MISCOM>','UNEXPECTED ENDIF STATEMENT.')
        endif
        COMLEN=0
     !
      ELSE IF (WRD.EQ.'IF  ') THEN
       OK2=peek_if()          ! execution flag at current level of nesting
       if (OK2) then          ! evaluate condition only if loop active, oherwise just push a 'false' (below)
!---- PROCESS-IF-COMMAND
        OK=.FALSE.
        CALL NEXTWD(COMLYN,COMLEN,TOKTMP,MXTLEN,LENTOK)
!     Do text replacement for a single character first operand.
!# <caves>-Aug-11-1993 (Leo Caves) if string (any size) is a symbol,
!                                  then substitute
! don't double substitute
        IF((FFOUR.NE. 'IF @').AND.(FFOUR.NE.'IF ?')) THEN
          IPAR=PARNUM(TOKTMP,LENTOK)
          IF(IPAR.GT.0) THEN
            IF(VALLEN(IPAR).GT.0) &
            CALL COPYST(TOKTMP,MXTLEN,LENTOK,VALNAM(IPAR),VALLEN(IPAR))
          ENDIF
        ENDIF
     !
        WRD=NEXTA4(COMLYN,COMLEN)
        CALL NEXTWD(COMLYN,COMLEN,VALTMP,MXVLEN,LENVAL)
     !
        IF(PRNLEV.GE.2) WRITE(OUTU,92) TOKTMP(1:LENTOK),&
                                       VALTMP(1:LENVAL)
  92    FORMAT(' Comparing "',A,'" and "',A,'".')
     !
        IF (WRD.EQ.'EQ  ' .OR. WRD.EQ.'.EQ.') THEN
          OK=EQST(TOKTMP,LENTOK,VALTMP,LENVAL)
        ELSE IF (WRD.EQ.'NE  ' .OR. WRD.EQ.'.NE.') THEN
          OK=.NOT.EQST(TOKTMP,LENTOK,VALTMP,LENVAL)
        ELSE
          VALUE=DECODF(TOKTMP,LENTOK)
          TEMP=DECODF(VALTMP,LENVAL)
          IF (WRD.EQ.'GT  ' .OR. WRD.EQ.'.GT.') THEN
            OK=VALUE.GT.TEMP
          ELSE IF (WRD.EQ.'GE  ' .OR. WRD.EQ.'.GE.') THEN
            OK=VALUE.GE.TEMP
          ELSE IF (WRD.EQ.'LT  ' .OR. WRD.EQ.'.LT.') THEN
            OK=VALUE.LT.TEMP
          ELSE IF (WRD.EQ.'LE  ' .OR. WRD.EQ.'.LE.') THEN
            OK=VALUE.LE.TEMP
          ELSE IF (WRD.EQ.'AE  ' .OR. WRD.EQ.'.AE.') THEN
            OK=ABS(VALUE-TEMP).LT.PT0001
          ELSE
! Better to trap this, an error could make affect a simulation
! without being noticed. L Nilsson
            call wrndie(-2,'<MISCOM>','Illegal IF test')
          ENDIF
        ENDIF
       else
     !       remove conditional without processing 
        CALL NEXTWD(COMLYN,COMLEN,TOKTMP,MXTLEN,LENTOK) ! remove token
        WRD=NEXTA4(COMLYN,COMLEN)                       ! remove operand
        CALL NEXTWD(COMLYN,COMLEN,VALTMP,MXVLEN,LENVAL) ! remove value
       endif
     ! for parallel : all nodes (not just root) scan the body of the conditional; the condition will in
     ! general, be false on some nodes (and true on others).
     !
       OK=OK.and.OK2           ! i.e. conditions at both levels need to be valid
       call push_if(OK) ! add level to the `if' stack 
     !
       IF (OK) THEN ! execute commands
          IF(PRNLEV.GE.2) WRITE(OUTU,94)
  94      FORMAT(' IF test evaluated as true.  Performing command')
          DONE=.false. ! tells miscom to re-start (see end of routine)
          CALL TRIMA(COMLYN,COMLEN)
          IF(COMLYN(1:4).EQ.'THEN') then 
            WRD=NEXTA4(COMLYN,COMLEN) ! i.e. remove 'THEN' from the line; this is a multi-line conditional terminated by ENDIF
          ELSE
           OK=pop_if() ! this is a one-line conditional without ELSE/ENDIF, so pop the stack (reuse 'OK' logical)
          ENDIF 
       ELSE ! I.E. (.NOT.OK) THEN ! will skip loop body
          WRD=NEXTA4(COMLYN,COMLEN)
          JUNK=NEXTA4(COMLYN,COMLEN)
     !
          IF((WRD.EQ.'THEN').AND.(JUNK.eq.'    '))THEN
              IF(PRNLEV.GE.2.and.OK2) WRITE(OUTU,'(A)') & ! write only if prev. level execution valid 
&               ' IF test evaluated as false.  Skip to ELSE or ENDIF'
          ELSE
             IF(PRNLEV.GE.2.and.OK2) WRITE(OUTU,'(A)') &
&               ' IF test evaluated as false.  Skipping command'
             COMLEN=0
             OK=pop_if() ! this is a one-line conditional without ELSE/ENDIF, so pop the stack
          ENDIF
       ENDIF ! OK/NOT OK
     !
      ELSE IF (.not.peek_if()) THEN
         COMLEN=0                ! ignore command line since execution is disabled on this node at this level of nesting
     !============================= DONE WITH CONDITIONALS ===================================================
      else & ! goes with "if" below
#endif /*(mcom) */
     if (wrd == 'BANN') then
        call header
        !
     else if (wrd == 'BOMB' .or. wrd.eq.'BOML') then
        kk=nexti(comlyn,comlen)
        if(kk <  -2)then
#if KEY_MULTICOM==1 /*  VO */
         if (mynod.eq.0) then                   
#endif
           if (prnlev >= 2) write(outu,'(a)') &
                ' MISCOM> Setting BOMLev < -2 is NOT a good idea.'
#if KEY_MULTICOM==1
         endif                                  
#endif
        endif
        bomlev=kk
        call set_param('BOMLEV',bomlev)
        !
     else if (wrd == 'CLOS') then
        call clolgu(comlyn,comlen,idummy)
        !
     ELSE IF (WRD == 'RXMD') THEN
#if KEY_RMD==1
        qcros = .true.
        call crossinit
#else /**/
        call wrndie(-1,'<MISCOM>', &
             'CROSS NOT AVAILABLE - NOT COMPILED WITH RMD FLAG')
#endif 
     ELSE IF (WRD == 'MRMD') THEN
#if KEY_MRMD==1
#if KEY_PARALLEL==1
        IF(MYNOD.EQ.0) &
#endif
        call mrmd_init
#else
        call wrndie(-1,'<MISCOM>', &
             'MRMD NOT AVAILABLE - NOT COMPILED WITH MRMD FLAG')
#endif
   ELSE IF (WRD.EQ.'VALB') THEN
#if KEY_VALBOND==1
         CALL VBCOMM()
#else /**/
         CALL WRNDIE(-1,'<MISCOM>', &
                'VALB NOT AVAILABLE - NOT COMPILED WITH VALBOND FLAG')
#endif 
     else if (wrd == 'MMPT') then
#if KEY_MMPT==1
        CALL MMPTINIT
#else /**/
        CALL wrndie(-1,'<MISCOM>', &
             'MMPT NOT AVAILABLE - NOT COMPILED WITH MMPT FLAG')
#endif 
     else if (wrd == 'DATE') then
        if(prnlev >= 2) then
           call daytim(k,l,m,ii,jj,kk)
           write(outu,54) k,l,m,ii,jj,kk
54         FORMAT(17x,'Current date: ',i2,'/',i2,'/',i2,' at ', &
                i2,':',i2,':',i2)
        endif

     else if (wrd == 'DECR') then
        !---- PROCESS-DECREMENT-COMMAND
        parop = decr

     else if (wrd == 'DEFI') then
        !---- process-define-command
        call filsky(comlyn,comlen,lused,.false.,idum)
        if (.not.(lused)) then
           call joinwd(comlyn,mxcms2,comlen,wrd,4)
           done=.true.
        endif
        !
     else if (wrd == 'DIVI') then
        !---- process-divide-command
        parop  = divi
        !
     else if (wrd == 'ECHU') then
        ! get output unit for echo command
        iecho=nexti(comlyn,comlen)
        if(iecho == 0) iecho=outu
     else if (wrd == 'ECHO') then
        ! just echo command line, without concern for its length, and skip first blank
        if(prnlev >= 2)then 
           if(comlen >= 2)then 
              write(iecho,'(a)') comlyn(2:comlen)
           else
              write(iecho,'(1x)')
           endif
        endif
        comlen=0

#if KEY_UNIX==1 /*unix*/
     else if (wrd == 'ENVI') then
        !---- process-environment-command
        call nextwd(comlyn, comlen, envvar, maxenv - 1, lenvar)
        call nextwd(comlyn, comlen, envval, maxenv - 1, lenval)

        ipt = 0
        do i = 1, lenvar
           if (envvar(i:i) /= '"') then  ! clean double quote
              ipt = ipt + 1
              envvar(ipt:ipt) = envvar(i:i)
           end if
        end do
        lenvar = ipt
        
        ipt = 0
        do i = 1, lenval
           if (envval(i:i) /= '"') then  ! clean double quote
              ipt = ipt + 1
              envval(ipt:ipt) = envval(i:i)
           end if
        end do
        lenval = ipt

        envvar = envvar(1:lenvar) // C_NULL_CHAR
        envval = envval(1:lenval) // C_NULL_CHAR
        i = setenv(envvar, envval, 1)

        if (i .ne. 0) then
           call wrndie(0, '<MISCOM>', &
                'failed to change environment variable')
        end if
#endif /* (unix)*/

     else if (wrd == 'EXPO') then
        !---- process-exponent-command
        parop  = expo

        !
     else if (wrd == 'FAST') then
        !---- process-fast-command
        faster=0
        wrd=nexta4(comlyn,comlen)
        fasteroption: select case(wrd)
        case('OFF ' ) fasteroption
           faster=-1
        case('GENE'  ) fasteroption
           faster=1
        case('ON  '  ) fasteroption
           faster=2
        case('DEFA') fasteroption
           faster=0
        case('EXPA') fasteroption
           faster=2
        case('SCAL') fasteroption
           CALL WRNDIE(0,'<MISCOM>','FAST option SCALar not available.')
        CASE('VECT') fasteroption
           CALL WRNDIE(0,'<MISCOM>','FAST option VECTor not available.')
        CASE('PARV') fasteroption
           CALL WRNDIE(0,'<MISCOM>','FAST option PARVec not available.')
        CASE('VPAR') fasteroption
           CALL WRNDIE(0,'<MISCOM>','FAST option VPAR not available.')
        CASE('CRAY') fasteroption
           CALL WRNDIE(0,'<MISCOM>','FAST option CRAY not available.')
#if KEY_PARALLEL==1
        case('DPCO') fasteroption
           dpcomm=1
        case('SPCO') fasteroption
           dpcomm=-1
#endif 
        case('BUCK') fasteroption
           IF(WRNLEV >= 2) WRITE (OUTU,'(A)') &
                ' MISCOM> Bucket fast routine is removed.'
        case default
           i=4
           faster=nexti(wrd,i)
           if(faster == 5) then
              if(wrnlev >= 2) write (outu,'(a)') &
                   ' MISCOM> Bucket fast routine is removed.'
           endif
           if(faster < -1 .or. faster > 3) then
              if(wrnlev >= 2) write (outu,'(a,i3)') &
                   ' MISCOM> Unknown FAST Level ',FASTER
              faster=0
              if(wrnlev >= 2) write(outu,'(a,i3)') &
                   '         FASTER will be set to ',FASTER
           endif
        end select fasteroption
        !
        !     set up LFAST value and inform the user the FASTER level set
        IF (FASTER < 0) THEN
           FASTER=-1
           LFAST =-1
#if KEY_PBOUND==1
           !     also turn of the hardwired periodic boundary in this case
           WRD='OFF'
           I=3
           CALL BOUND(WRD,I)
#endif 
        ENDIF
        LFAST=FASTER
        !
        IF (FASTER == -1) THEN
           IF(PRNLEV >= 2) WRITE (OUTU,'(A)') &
                ' MISCOM> FAST option: OFF (full feature routines)'
        ELSE IF (FASTER == 0) THEN
           IF(PRNLEV >= 2) WRITE (OUTU,'(A)')  &
                ' MISCOM> FAST option: DEFAULT (use best available)'
        ELSE IF (FASTER == 1) THEN
           IF(PRNLEV >= 2) WRITE (OUTU,'(A)') &
                ' MISCOM> FAST option: ON (generic fast routines)'
        ELSE IF (FASTER == 2) THEN
           IF(PRNLEV >= 2) WRITE (OUTU,'(A)') &
                ' MISCOM> FAST option: EXPANDED (limited fast routines)'
        ELSE
           IF(PRNLEV >= 2) WRITE (OUTU,'(A)')  &
                ' MISCOM> Bad FAST option: (coding error)'
        ENDIF
        CALL set_param('LFAST',LFAST)
        CALL set_param('FASTER',FASTER)
        !
     ELSE IF (WRD == 'FORM') THEN
        !---- PROCESS-FORM-COMMAND
        CALL NEXTWD(COMLYN,COMLEN,FMTWD,FMTMAX,FMTLEN)
        !
        !yw++
     ELSE IF (WRD == 'IOFO') THEN
        !---- process IOFOrmat command
        !     IOFOrmat [EXTEnded  ]
        !              [NOEXtended]
        qoldfmt=.false.
        qoldfmt=indxa(comlyn,comlen,'NOEX') > 0
        if (qoldfmt) then
           qextfmt=.false.
           idleng=4
           if (prnlev >= 2) write(outu,'(a)') &
                ' MISCOM> Expanded I/O format is NOT used.'
        endif
        qnewfmt=.false.
        qnewfmt=indxa(comlyn,comlen,'EXTE') > 0
        if (qnewfmt) then
           qextfmt=.true.
           idleng=8
           if (prnlev >= 2) write(outu,'(a)') &
                ' MISCOM> Expanded I/O format is used.'
        endif
        !yw--
     ELSE IF (WRD == 'GET ') THEN
        !---- GET-COMMAND-PARAMETER
        IUNIT=GTRMI(COMLYN,COMLEN,'UNIT',-1)
        CALL NEXTWD(COMLYN,COMLEN,TOKTMP,MXTLEN,LENTOK)
        IF (IUNIT  >  0) THEN
           CALL RDCMND(COMLYN,MXCMS2,COMLEN,IUNIT,EOF, &
#if KEY_MULTICOM==0 /*  VO broadcast value of get to all local nodes */
                .FALSE.,&
#else
                .TRUE., &
#endif
                .FALSE.,' ')
#if KEY_PARALLEL==1 || KEY_ENSEMBLE==1
           !     Do not use parallel code in RDCMND since we should not compress
           !     command line parameters.
#if KEY_MULTICOM==0 /*  VO skip broadcast; it must be done within RDCMND for some communicators */
           CALL PSND4(COMLEN,1)
           CALL PSND4(EOF,1)
#endif
           IF(COMLEN > 0) CALL PSNDC(COMLYN(1:COMLEN),1)
#endif 
           IF(EOF) THEN
              COMLYN='END-OF-FILE'
              COMLEN=11
              IF(WRNLEV >= 2) WRITE(OUTU,3320) IUNIT
3320          FORMAT(' MISCOM WARNING:', &
                   ' EOF when getting parameter from unit',I4)
           ENDIF
           CALL TRIMA(COMLYN,COMLEN)
           IPAR = PARINS(TOKTMP,LENTOK,COMLYN,COMLEN)
           IF(IPAR < 0) THEN
              CALL WRNDIE(0,'<MISCOM>','Failed to install parameter.')
           ENDIF
        ELSE
           CALL WRNDIE(0,'<MISCOM>','Illegal or missing UNIT for GET')
        ENDIF
        COMLEN=0
        !
        !
     ELSE IF (WRD == 'GOTO') THEN
        !---- PROCESS-GOTO-COMMAND
        CALL NEXTWD(COMLYN,COMLEN,WRD2,WDMAX,WD2LEN)
#if KEY_MULTICOM==0 /*  VO there could be several parser communicators (even with iolev < 0) */
        IF(IOLEV > 0) THEN
#else
        if (ME_PARSER.eq.0) then
#endif
           REWIND ISTRM
           OK=.FALSE.
           EOF=.FALSE.
111        CONTINUE
           READ(ISTRM,'(A)',ERR=945) COMLYN(1:80)
           COMLEN=80
           CALL CNVTUC(COMLYN,COMLEN)
           WRD=NEXTA4(COMLYN,COMLEN)
           IF (WRD == 'LABE') THEN
              CALL NEXTWD(COMLYN,COMLEN,WRD3,WDMAX,WD3LEN)
              OK=EQST(WRD2,WD2LEN,WRD3,WD3LEN)
           ENDIF
           IF (.NOT.(OK.OR.EOF)) GOTO 111
           !
           IF (.NOT.EOF) GOTO 845
945        CALL WRNDIE(-2,'<MISCOM>','Unable to find GOTO label')
845        CONTINUE
        ENDIF
        !
#if KEY_MULTICOM==0 /* (mcom)  VO : stringm conditionals already processed above */
     ELSE IF (WRD == 'IF  ') THEN
        !---- PROCESS-IF-COMMAND
        OK=.FALSE.
        CALL NEXTWD(COMLYN,COMLEN,TOKTMP,MXTLEN,LENTOK)
        !     Do text replacement for a single character first operand.
        !# <caves>-Aug-11-1993 (Leo Caves) if string (any size) is a symbol,
        !                                  then substitute
        ! don't double substitute
        IF((FFOUR /=  'IF @').AND.(FFOUR.NE.'IF ?')) THEN
           IPAR=PARNUM(TOKTMP,LENTOK)
           IF(IPAR > 0) THEN
              IF(VALLEN(IPAR) > 0) &
                   CALL COPYST(TOKTMP,MXTLEN,LENTOK,VALNAM(IPAR),VALLEN(IPAR))
           ENDIF
        ENDIF
        !
        WRD=NEXTA4(COMLYN,COMLEN)
        CALL NEXTWD(COMLYN,COMLEN,VALTMP,MXVLEN,LENVAL)
        !
        IF(PRNLEV >= 2) WRITE(OUTU,92) TOKTMP(1:LENTOK), &
             VALTMP(1:LENVAL)
92      FORMAT(' Comparing "',A,'" and "',A,'".')
        !
        IF (WRD == 'EQ  ' .OR. WRD.EQ.'.EQ.') THEN
           OK=EQST(TOKTMP,LENTOK,VALTMP,LENVAL)
        ELSE IF (WRD == 'NE  ' .OR. WRD.EQ.'.NE.') THEN
           OK=.NOT.EQST(TOKTMP,LENTOK,VALTMP,LENVAL)
        ELSE
           VALUE=DECODF(TOKTMP,LENTOK)
           TEMP=DECODF(VALTMP,LENVAL)
           IF (WRD == 'GT  ' .OR. WRD.EQ.'.GT.') THEN
              OK=VALUE > TEMP
           ELSE IF (WRD == 'GE  ' .OR. WRD.EQ.'.GE.') THEN
              OK=VALUE >= TEMP
           ELSE IF (WRD == 'LT  ' .OR. WRD.EQ.'.LT.') THEN
              OK=VALUE < TEMP
           ELSE IF (WRD == 'LE  ' .OR. WRD.EQ.'.LE.') THEN
              OK=VALUE <= TEMP
           ELSE IF (WRD == 'AE  ' .OR. WRD.EQ.'.AE.') THEN
              OK=ABS(VALUE-TEMP) < PT0001
           ELSE
! Better to trap this, an error could make affect a simulation
! without being noticed. L Nilsson
            call wrndie(-2,'<MISCOM>','Illegal IF test')
           ENDIF
        ENDIF
        !
        IF (OK) THEN
           IF(PRNLEV >= 2) WRITE(OUTU,94)
94         FORMAT(' IF test evaluated as true.  Performing command')
           DONE=.FALSE.
           CALL TRIMA(COMLYN,COMLEN)
           IF(COMLYN(1:4) == 'THEN') WRD=NEXTA4(COMLYN,COMLEN)
        ENDIF
        IF (.NOT.OK) THEN
           WRD=NEXTA4(COMLYN,COMLEN)
           JUNK=NEXTA4(COMLYN,COMLEN)
           IF((WRD == 'THEN').AND.(JUNK.eq.'    '))THEN
              IF(PRNLEV >= 2) WRITE(OUTU,'(A)')  &
                   ' IF test evaluated as false.  Skip to ELSE or ENDIF'
              IF(IOLEV > 0) THEN
                 ILEVEL=1
                 OK=.FALSE.
112              CONTINUE
                 READ(ISTRM,'(A)',END=946,ERR=946) COMLYN(1:80)
                 COMLEN=80
                 CALL CNVTUC(COMLYN,COMLEN)
                 WRD=NEXTA4(COMLYN,COMLEN)
                 IF (WRD == 'IF  '.and. INDEX(COMLYN,'THEN') /= 0)THEN
                    ILEVEL=ILEVEL+1
                 ENDIF
                 IF (WRD == 'ENDI')THEN
                    ILEVEL=ILEVEL-1
                 ENDIF
                 OK=(ILEVEL == 0)
                 IF ((WRD == 'ELSE').AND.(ILEVEL <= 1)) OK=.TRUE.
                 IF (.NOT.OK) GOTO 112
946              IF (.NOT.OK) THEN    
                    CALL WRNDIE(-2,'<MISCOM>', &
                         'Unable to find ELSE or ENDIF')
                 ENDIF
              ENDIF
           ELSE
              IF(PRNLEV >= 2) WRITE(OUTU,'(A)')  &
                   ' IF test evaluated as false.  Skipping command'
              COMLEN=0
           ENDIF
        ENDIF
        !
#endif /*(mcom) */
     ELSE IF (WRD == 'INCR') THEN
        !---- PROCESS-INCREMENT-COMMAND
        parOp = INCR
        !
#if KEY_MULTICOM==0 /* (mcom)  VO stringm conditionals already processed above */
     ELSE IF (WRD == 'ENDI') THEN
        COMLEN=0
        !
#endif /*(mcom) */
     ELSE IF (WRD == 'LABE') THEN
        COMLEN=0
!        FLUSH (OUTU)  ! for MMTSB
        !
#if KEY_MULTICOM==0 /* (mcom)  VO stringm conditionals already processed */
     ELSE IF (WRD == 'ELSE') THEN
        IF(PRNLEV >= 2)THEN
           WRITE(OUTU,'(2A)')  &
                ' Skip commands until next ENDIF'
        ENDIF
        IF(IOLEV > 0) THEN
           OK=.FALSE.
           ILEVEL=1
113        CONTINUE
           READ(ISTRM,'(A)',END=947,ERR=947) COMLYN(1:80)
           COMLEN=80
           CALL CNVTUC(COMLYN,COMLEN)
           WRD=NEXTA4(COMLYN,COMLEN)
           IF (WRD == 'IF  '.and. INDEX(COMLYN,'THEN') /= 0)THEN
              ILEVEL=ILEVEL+1
           ENDIF
           IF (WRD == 'ENDI')THEN
              ILEVEL=ILEVEL-1
           ENDIF
           OK=(ILEVEL == 0)
           IF (.NOT.OK) GOTO 113
947        IF (.NOT.OK) THEN
              CALL WRNDIE(-2,'<MISCOM>','Unable to find ENDIF')
           ENDIF
        ENDIF
        !
#endif /*(mcom) */
     ELSE IF (WRD == 'CALC') THEN
        ! Dec-01-1994 (Benoit Roux ) arithmetic parameter manipulation
        ! get the token
        CALL NEXTWD(COMLYN,COMLEN,TOKTMP,MXTLEN,LENTOK)
        CALL TRIMA(COMLYN,COMLEN)
        ! grab the parameter which is to be modified
        IPAR = PARNUM(TOKTMP,LENTOK)
        IF(IPAR <= 0) THEN
           ! install
           !         CALL WrnDie(0,'<MISCOM>',
           !    ,              'Non-existent parameter. '//swdtch(1:swdlen))
           IF(NUMPAR+1 > MAXPAR)THEN
              CALL WrnDie(0,'<MISCOM>','Parameter table full. Ignored.')
           ENDIF
           NUMPAR = NUMPAR + 1
           IPAR = NUMPAR
           TOKNAM(IPAR) = TOKTMP(1:LENTOK)
           TOKLEN(IPAR) = LENTOK
        ENDIF
        !
        ! dispense with an assignment operator here.
        IF (COMLYN(1:2) == '= ') THEN
           COMLYN(1:2) = '  '
           CALL TRIMB(COMLYN,COMLEN)
        ENDIF
        ! evaluate the arithmetic expression
        call calculator(VALUE,COMLYN)
        ! Set the new value to the parameter
        CALL ENCODF(value,VALNAM(IPAR),MXVLEN,VALLEN(IPAR))
        IF(PRNLEV >= 2) CALL PARWRI(OUTU,IPAR,1)
        COMLEN = 0
        !
     ELSE IF (WRD == 'LET ') THEN
        !# <caves>-Aug-13-1993 (Leo Caves) Limited parameter manipulation
        IPAR =  PARMAN(COMLYN,COMLEN)
        IF(IPAR < 0) THEN
           CALL WrnDie(0,'<MISCOM>','Error in parameter modification') 
           COMLEN = 0
        ENDIF
        !
     ELSE IF (WRD == 'LOWE') THEN
        LOWER=.TRUE.
        !
     ELSE IF (WRD == 'MULT') THEN
        !---- PROCESS-MULTIPLY-COMMAND
        parOp = MULT
        !
     ELSE IF (WRD.EQ.'MTP') THEN
        !  atomic multipole moments MTP module
        CALL MTP
        !
     ELSE IF (WRD.EQ.'MTPL') THEN
        !  atomic multipole moments MTPL module--uses local axis systems for
        !  arbitrarily large molecules
#if KEY_MTPL==1
        CALL MTPL
#else
        CALL WrnDie(0,'<MISCOM>','MTPL module not compiled.')
#endif /* mtpl */
        !
     ELSE IF (WRD == 'NOBO') THEN
        BOMLEV=-1
        CALL set_param('BOMLEV',BOMLEV)
     ELSE IF (WRD == 'OPEN') THEN
        CALL OPNLGU(COMLYN,COMLEN,IDUMMY)
     ELSE IF (WRD == 'OUTU') THEN
        OUTU  =NEXTI(COMLYN,COMLEN)
        CALL set_param('OUTU',OUTU)
#if KEY_REPDSTR==1 || KEY_REPDSTR2==1
        QWRQTT=.TRUE.                                  
        CALL DREPSETIO(IOLEV,PRNLEV,WRNLEV)                   
#endif
     ELSE IF (WRD == 'PRLE' .OR. WRD.EQ.'PRNL') THEN
        I=NEXTI(COMLYN,COMLEN)
#if KEY_PARALLEL==1
        J=GTRMI(COMLYN,COMLEN,'NODE',-100)
        IF(J == -100) THEN
           PRNLEV=I
        ELSE IF(J < 0 .OR. J > NUMNOD) THEN
           CALL WRNDIE(0,'<MISCOM>','Bad node number - ignored')
        ELSE
           IF(MYNOD == J) PRNLEV=I
        ENDIF
        plnod0=prnlev
        IF(MYNOD /= 0) plnod0=0
        call gbor(plnod0,1)

#else /**/
        PRNLEV=I
#endif 
        CALL set_param('PRNLEV',PRNLEV)
#if KEY_NOMISC==0
     ELSE IF (WRD == 'QUIC' .OR. WRD.EQ.'Q') THEN
        CALL QUICKA()
     ELSE IF (WRD == 'AIDX') THEN
        CALL AIDX()
#endif 
     ELSE IF (WRD == 'RAND') THEN
        CALL RANDSPEC()
        !
     ELSE IF (WRD == 'IRAN') THEN
        CALL IRANDOM()
        !
     ELSE IF (WRD == 'RETU') THEN
        CALL PPSTRM(OK)
        !.....mf050629
     ELSE IF (WRD == 'REWF') THEN
        reallow=.false.
     ELSE IF (WRD == 'REWT') THEN
        reallow=.true.
        !
     ELSE IF (WRD == 'REWI') THEN
        !
        ! Rewind file (useful for repetitive stream changing)
        !
        IUNIT=GTRMI(COMLYN,COMLEN,'UNIT',-1)
        IF (IUNIT > 0) THEN
           IF(IOLEV > 0) REWIND IUNIT
           IF(PRNLEV >= 2) WRITE(OUTU,1220) IUNIT
1220       FORMAT(20X,'REWINDING UNIT ',I5)
        ELSE
           CALL WRNDIE(0,'<MISCOM>','Logical unit number incorrect.')
        ENDIF
        !
     ELSE IF (WRD == 'SET ') THEN
        !---- SET-COMMAND-PARAMETER
        !# <caves>-Aug-10-1993 (Leo Caves) New parameter table
        ! get the token
        CALL NEXTWD(COMLYN,COMLEN,TOKTMP,MXTLEN,LENTOK)
        CALL TRIMA(COMLYN,COMLEN)
        !
        ! dispense with an assignment operator here.
        IF (COMLYN(1:2) == '= ') THEN
           COMLYN(1:2) = '  '
           CALL TRIMB(COMLYN,COMLEN)
        ENDIF
        ! install token
        IPAR = PARINS(TOKTMP,LENTOK,COMLYN,COMLEN)
        IF(IPAR < 0) THEN
           CALL WrnDie(0,'<MISCOM>','Failed to install parameter.')
        ENDIF
        COMLEN = 0
        !
        !
     ELSE IF (WRD == 'SHOW') THEN
        !---- process-show-command
        WRD = NEXTA4(COMLYN,COMLEN)
        ! this command not relevant if I/O level not appropriate
        IF((WRD == '   ').OR.(WRD.EQ.'BUIL'))THEN
           IF(PRNLEV >= 2) THEN
              WRITE(OUTU,143)
              call write_real_params(outu)
              call write_int_params(outu)
              call write_str_params(outu)
           ENDIF
        ELSE IF (WRD(1:3) == 'PAR') THEN
           ! dump parameter table
           IUNIT=GTRMI(COMLYN,COMLEN,'UNIT',OUTU)
           IF(PRNLEV >= 2) THEN
              WRITE(OUTU,'(A,2(I5,A))') &
                   ' MISCOM: No. of parameters= ',NUMPAR,'  (max.=',MAXPAR,')'
              ! verbose gets you string limits
              IF(INDXA(COMLYN,COMLEN,'VERB') > 0)  &
                   WRITE(OUTU,'(2(A,I5))')  &
                   ' MISCOM: Maximum string sizes: tokens= ',MXTLEN, &
                   ' values= ',MXVLEN
              IF(NUMPAR >= 1) THEN
                 ! here's the table
                 DO IPAR = 1, NUMPAR
                    CALL PARWRI(IUNIT,IPAR,-1)
                 ENDDO
              ENDIF
           ENDIF
        ELSE
           ! if the show command keyword was unrecognized,
           ! perhaps this is the CORREL SHOW command.  Reset the command line.
           CALL JOINWD(COMLYN,MXCMS2,COMLEN,WRD,4)
           WRD2='SHOW '
           CALL JOINWD(COMLYN,MXCMS2,COMLEN,WRD2,5)
           LUSED=.FALSE.
           DONE=.TRUE.
        ENDIF

143     FORMAT(' MISCOM: List of command substitution parameter values')

     ELSE IF (WRD == 'STRE') THEN
        !
#if KEY_REPDSTR==1 
!                   || KEY_REPDSTR2==1
        IF(QREPDSTR) THEN
!!!           CALL DREPSETIO(IOLEV,PRNLEV)
           !!           IF(MYNOD == 0) QRDQTT=.TRUE.
           QRDQTT=.TRUE.
        ENDIF
#endif 
        CALL PUSTRM(COMLYN,MXCMS2,COMLEN)
        !
        ! Process automatic set options: 
        !  STREam filename   parm1 parm2 parm3 ...
        !     this becomes:   IN1   IN2   IN3  ... 
        IF(COMLEN > 0) THEN
           I=0
           DO WHILE(COMLEN > 0 .AND. I < 9)
              I=I+1
              WRITE(TOKTMP,166) I
166           FORMAT('IN',I1)
              LENTOK=3
              CALL TRIMA(COMLYN,COMLEN)
              CALL NEXTWD(COMLYN,COMLEN,WRD2,WDMAX,WD2LEN)
              ! install token
              IPAR = PARINS(TOKTMP,LENTOK,WRD2,WD2LEN)
              IF(IPAR < 0) THEN
                 CALL WrnDie(0,'<MISCOM>','Failed to install parameter.')
              ENDIF
           ENDDO
        ENDIF
        !
     ELSE IF (WRD == 'SKIP') THEN
        CALL SKIPE(COMLYN,COMLEN)
        !
     ELSE IF (WRD == 'TIME') THEN
        !---- PROCESS-TIME-COMMAND
        IF (INDXA(COMLYN,COMLEN,'NOW')  >  0) THEN
           CALL TIMRE
        ELSE IF (INDXA(COMLYN,COMLEN,'DIFF')  >  0) THEN
           CALL TIMRE
           CALL TIMRB
        ELSE
           TIMER=NEXTI(COMLYN,COMLEN)
           CALL set_param('TIMER',TIMER)
        ENDIF
        !
#if KEY_LOOKUP==1
     ELSE IF (WRD == 'NBSO'.OR.WRD.EQ.'LOOK')THEN
        CALL WWSETUP(COMLYN,COMLEN) 
#endif 
        !
     ELSE IF (WRD == 'TITL') THEN
        !---- PROCESS-TITLE-COMMAND
        IF (INDXA(COMLYN,COMLEN,'COPY') > 0) THEN
           NTITLA=NTITLB
           DO I=1,NTITLB
              TITLEA(I)=TITLEB(I)
           ENDDO
        ENDIF
        CALL RDTITL(TITLEA,NTITLA,ISTRM,0)
        !
#if KEY_BLOCK==1
     ELSE IF (WRD == 'LDTI') THEN
        !---- PROCESS-TITLE-COMMAND FOR Lambda Dynamics
        CALL RDTITL(TITLEL,NTITLL,ISTRM,0)
#endif 
        !
     ELSE IF (WRD == 'TRIM') THEN
        !---- PROCESS-TRIM-COMMAND
        CALL NEXTWD(COMLYN,COMLEN,SWDTCH,SWDMAX,SWDLEN)
        IPAR=PARNUM(SWDTCH,SWDLEN)
        IF (IPAR > 0) THEN
           I=VALLEN(IPAR)
           DO WHILE(I > 1 .AND. VALNAM(IPAR)(I:I) == ' ')
              I=I-1
           ENDDO
           J=1
           DO WHILE(J < I .AND. VALNAM(IPAR)(J:J) == ' ')
              J=J+1
           ENDDO
           I=GTRMI(COMLYN,COMLEN,'TO',I)
           J=GTRMI(COMLYN,COMLEN,'FROM',J)
           IF(I > MXVLEN)I=MXVLEN
           IF(I > VALLEN(IPAR)) &
                CALL FILSPC(VALNAM(IPAR),I,VALLEN(IPAR))
           CALL COPSUB(VALNAM(IPAR),MXVLEN, &
                VALLEN(IPAR),VALNAM(IPAR),J,I)
           IF(PRNLEV >= 2) CALL PARWRI(OUTU,IPAR,1) 
        ELSE
           CALL WrnDie(0,'<MISCOM>', &
                ' Parameter not found. Nothing to trim.')
        ENDIF
        !
     ELSE IF (WRD == 'UPPE') THEN
        LOWER=.FALSE.
     ELSE IF (WRD == 'USE ') THEN
        CALL SETUPFF(COMLYN,COMLEN)
        !
     ELSE IF (WRD == 'USER') THEN
        CALL USERSB
     ELSE IF (WRD == 'WRNL') THEN
        I=NEXTI(COMLYN,COMLEN)
#if KEY_PARALLEL==1
        J=GTRMI(COMLYN,COMLEN,'NODE',-100)
        IF(J == -100) THEN
           WRNLEV=I
#if KEY_MULTICOM==0 /*  VO stringm v */
        ELSE IF(J < 0 .OR. J >= NUMNOD) THEN
#else
        ELSE IF(J < 0 .OR. J >= SIZE_PARSER) THEN
#endif
           CALL WRNDIE(0,'<MISCOM>','Bad node number - ignored')
        ELSE
#if KEY_MULTICOM==0
           IF(MYNOD == J) WRNLEV=I
#else
           IF(ME_PARSER == J) WRNLEV=I
#endif /* VO stringm v */
        ENDIF
#else /**/
        WRNLEV=I
#endif 
        CALL set_param('WRNLEV',WRNLEV)
     ELSE IF (WRD == 'IOLE') THEN
        I=NEXTI(COMLYN,COMLEN)
#if KEY_PARALLEL==1
        J=GTRMI(COMLYN,COMLEN,'NODE',-100)
        IF(J == -100) THEN
           IOLEV=I
#if KEY_MULTICOM==0 /*  VO stringm v */
        ELSE IF(J < 0 .OR. J >= NUMNOD) THEN
#else
        ELSE IF(J < 0 .OR. J >= SIZE_PARSER) THEN
#endif
           CALL WRNDIE(0,'<MISCOM>','Bad node number - ignored')
        ELSE
#if KEY_MULTICOM==0
           IF(MYNOD == J) IOLEV=I
#else
           IF(ME_PARSER == J) IOLEV=I
#endif /* VO stringm v */
        ENDIF
#else /**/
        IOLEV=I
#endif 
        CALL set_param('IOLEV',IOLEV)
     ELSE IF (WRD == 'LONG') THEN
        QLONGL=.TRUE.
     ELSE IF (WRD == 'SHOR') THEN
        QLONGL=.FALSE.
        !
     ELSE IF (WRD == 'SYST') THEN
        CALL CSYSTEM(COMLYN,COMLEN)

     ELSE IF (WRD == 'FRCD') THEN
        iunit = GTRMI(COMLYN,COMLEN,'UNIT',outu)
        call DumpEnFrc(iunit)

     ELSE IF (WRD == 'STOP') THEN
        CALL STOPCH('NORMAL STOP')
        !
     ELSE
        CALL JOINWD(COMLYN,MXCMS2,COMLEN,WRD,4)
        LUSED=.FALSE.
        DONE=.TRUE.
     ENDIF
     !
     !# <caves>-Aug-13-1993 (Leo Caves)
     ! support existing parameter operations. prefer that any further operations
     ! be placed in a separate evaluation function. 
     ! eg LET command: SUBROUTINE PARMAN
     IF (parOp > 0) THEN
        INCREM=GTRMF(COMLYN,COMLEN,'BY',ONE)
        CALL NEXTWD(COMLYN,COMLEN,SWDTCH,SWDMAX,SWDLEN)
        ! lookup
        IPAR = PARNUM(SWDTCH,SWDLEN)
        IF (IPAR >= 0) THEN
           VALUE=DECODF(VALNAM(IPAR),VALLEN(IPAR))
           IF      (parOp == INCR) THEN
              VALUE = VALUE + INCREM
           ELSE IF (parOp == DECR) THEN
              VALUE = VALUE - INCREM
           ELSE IF (parOp == MULT) THEN
              VALUE = VALUE * INCREM
           ELSE IF (parOp == DIVI) THEN
              IF (INCREM /= ZERO) THEN
                 VALUE = VALUE / INCREM
              ELSE
                 CALL WRNDIE(0,'<MISCOM>','Divide by zero. IGNORED.')
              ENDIF
           ELSE IF (parOp == EXPO) THEN
              VALUE = EXP(VALUE)
           ENDIF
           !
           CALL ENCODF(VALUE,VALNAM(IPAR),MXVLEN,VALLEN(IPAR))
           IF (PRNLEV >= 2) CALL PARWRI(OUTU,IPAR,1)
        ELSE
           CALL WrnDie(0,'<MISCOM>',  &
                'Parameter not found. Nothing to modify.')
        ENDIF ! IPAR
     ENDIF ! parOp
     !
  enddo loop101
  !
  if (lused) call xtrane(comlyn,comlen,'miscom')
  return
end subroutine miscom

!---------------------------------------------------------------------
!          SETUPFF
!---------------------------------------------------------------------
subroutine setupff(comlyn,comlen)
  !----------------------------------------------------------------------
  use chm_kinds
  use dimens_fcm
#if KEY_DOMDEC==1
  use domdec_common, only: q_domdec
#endif
  use consta
  use energym
  use ffieldm
  use inbnd
  use machdep
#if KEY_PERT==1
  use pert, only: qpert
#endif
  use psf
  use param_store, only: set_param
#if KEY_CFF==1
  use cff_fcm  
#endif
#if KEY_MMFF==1
  use mmffm  
#endif
  use stream
  use string

  implicit none

  character(len=*) comlyn
  integer   comlen
  integer   i

  if(indxa(comlyn,comlen,'AMBER') > 0) then
     ffield=amberffn
     ccelec=ccelec_amber
     call set_param('CCELEC',ccelec)
     if(prnlev >= 2) write(outu,'(" AMBER Force Field will be used")')
  endif
     
  if(indxa(comlyn,comlen,'CHARMM') > 0) then
     ffield=charmm
     ccelec=ccelec_charmm
     call set_param('CCELEC',ccelec)
     if(prnlev >= 2) write(outu,'(" CHARMM Force Field will be used")')
  endif
     
#if KEY_MMFF==1 || KEY_CFF==1 /*mmff_cff*/
#if KEY_CFF==1
  if(indxa(comlyn,comlen,'CFF') > 0) then
     ffield=cff
     ! Allocate data structures for CFF force fields. cb3
     if(.not. allocated(itflg)) call allocate_cff
  endif
#endif 
#if KEY_MMFF==1
  if(indxa(comlyn,comlen,'MMFF') > 0) then

#if KEY_DOMDEC==1
     if (q_domdec) then
        call wrndie(-2,'<SETUPFF>','MMFF is not compatible with DOMDEC')
     end if
#endif /* domdec */

#if KEY_PERT==1
     if (qpert) then
        call wrndie(-2,'<SETUPFF>','MMFF is not compatible with PERT')
     end if     
#endif /* pert */
     
     e14fac=0.75
     v14fac=1.
     ffield=mmff
     ! Allocate data structures for CFF force fields. cb3
     if(.not. allocated(bondtype)) then
        call allocate_mmff
        call mmff_init
     endif
     if(indxa(comlyn,comlen,'TYPE') > 0) then
        call mmff_setup
        i=indxa(comlyn,comlen,'ATOM')
     endif
  endif
#endif 

  if(indxa(comlyn,comlen,'CHARMM') > 0) then
     ffield=charmm
     if(prnlev >= 2) write(outu,'(a)') &
          ' CHARMM Force Field will be used'
  endif
#else /* (mmff_cff)*/
  if(indxa(comlyn,comlen,'CFF') > 0) then
     call wrndie(-2,'<SETUPFF>','CFF code not compiled')
  endif
  if(indxa(comlyn,comlen,'MMFF') > 0) then
     call wrndie(-2,'<SETUPFF>','MMFF code not compiled')
  endif
  if(indxa(comlyn,comlen,'CHARMM') > 0) then
     if(prnlev >= 2) write(outu,'(a)') &
          ' CHARMM Force Field will be used'
  endif
#endif /* (mmff_cff)*/
  i=indxa(comlyn,comlen,'FORCE')
  i=indxa(comlyn,comlen,'FIELD')
  return
end subroutine setupff

