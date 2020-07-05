!
!     THIS GROUP OF STRING COMMANDS DEALS WITH I/O
!
SUBROUTINE RDCMND(COMLYN,MXCMS2,COMLEN,UNIT,EOF, &
     QCMPRS,QPRINT,ECHOST)
  !-----------------------------------------------------------------------
  !     This routine reads a command in from UNIT. The command
  !     field on a line may extend to column 100 or may be terminated by
  !     an exclamation mark anywhere on the line. A command is one command
  !     field of information unless the last nonblank character in that
  !     field is a hyphen. In that case, the command field on the next
  !     record is appended to the characters preceding the hyphen.
  !     Trailing blanks are removed from all command fields, and lowercase
  !     letters are converted to uppercase. If an end of file is seen
  !     while reading a command, EOF is turned on. If EOF is turned on
  !     when called, the routine returns immediately. If QPRINT is on,
  !     the records are printed on unit OUTU as they are read. Each line
  !     printed on output device is prefixed by the ECHOST so that it will
  !     be easier to find it when scanning the output.
  !
  !     36 parameters can be specified by the syntax @n where n is a
  !     alphanumeric characher {0-9,a-z,A-Z}.  The specified parameter
  !     string n (stored in CMDPAR.FMC) will be substituted for @n.
  !
  !     The flag QCMPRS tells whether to compress the string to remove
  !     extra blanks.  This should only be done from standard parsing
  !     sections.  For the parallel versions, this flag is used to control
  !     communication broadcasting (e.g. the string is broadcast only if
  !     this flag is set). - BRB
  !
  !      Authors: Robert Bruccoleri
  !               David States
  !               Bernie Brooks
  !               Youngdo Won
  !      Victor Anisimov, 2004
  !         Disable HYPEN mechanism when RDCMND is invoked from G94INP module,
  !         which creates Gaussian input file.
  !
  !      Modified by Rick Lapp  5 Aug 99 to turn off upper case conversion
  !        while reading CFF forcefield.
  !
#if KEY_REPDSTR==1
  use repdstrmod                           
#endif
#if KEY_REPDSTR2==1
  use repdstrmod2                           
  use repd_ensemble
#endif
  !
  use chm_kinds
  use dimens_fcm
  use exfunc
  use cmdpar,only:parse1
  use string
  use stream
  use parallel
  use repdstr
#if KEY_PARALLEL==1
  use mpi       
#endif
  ! VO begin
#if KEY_MULTICOM==1
  use mpi                                 
  use multicom_aux                        
  use ifstack
#endif
  ! VO end
#if KEY_ENSEMBLE==1
  use ensemble,only:ensprint     
#endif
#if KEY_CFF==1
  use rtf,only:ucase
#endif 
  implicit none
  !
  integer mxcms2,comlen,unit
  character(len=*) comlyn
  !
  !
  integer enrlen
  integer,parameter :: enrmax=20, mxcard=200
  character(len=enrmax) enrst
  integer cardln,ierror
  character(len=mxcard) card
  integer wdlen,iend,ipar,i,j
  logical eof,qprint,qcmprs,usehyp
  character(len=*) echost
  character(len=80) chint
#if KEY_REPDSTR==1 || KEY_REPDSTR2==1
  integer mynodt,numnodt  ,ierr                   
#endif
  logical qmonoscript
  character(len=1) :: HYPHEN='-', EXCLMK='!', ATMARK='@', SDBLQ='"', &
       SQUES='?',stemp
  !
#if KEY_ENSEMBLE==1
  write(chint,*)qcmprs,iolev                       
#endif
#if KEY_ENSEMBLE==1
  call ensprint("RDCMND>>> qcmprs,iolev",chint)    
#endif
#if KEY_PARALLEL==1
#if KEY_REPDSTR==1 
!                   || KEY_REPDSTR2==1
  if(.not.qrepdstr) then                     
#endif
     if(.not.qcmprs .and. iolev < 0) then
        comlen=0
        return
     endif
#if KEY_REPDSTR==1 
!                    || KEY_REPDSTR2==1
  endif                                      
#endif
#endif 
  if(eof) return

  if(index(echost,'G94INP>') > 0) then
     ! G94 input file case: do not use HYPHEN mechanism
     usehyp=.false.
  else
     ! Standard case: use HYPEN mechanism
     usehyp=.true.
  endif

  !--- Clean up the previous command line in case of something left
  comlyn=' '
  comlen=0
  if(unit < 0) goto 9
  if(qprint .and. prnlev >= 3) write(outu,201)
201 format('  ')



1 continue
  !
  card(1:mxcard)=' '
  !     In case of distributed I/O we read from "every" node
#if KEY_REPDSTR==1 
!------------------|| KEY_REPDSTR2==1
  !
  !   In case of repdstr=.true. we break the rule of IOLEV
  !   since we have two different cases to parse charmm input scripts:
  !   1. main input script: qrdqtt=.f.
  !   2. stream, ie each group has its own charmm commands: qrdqtt=.t.
  !
  !
  qmonoscript=qrepdstr.and.(.not.qrdqtt)

  !! The next two lines completely break subcommand parsing. I have
  !! no idea why Milan did this -- Tim Miller

  qmonoscript=qmonoscript.and.(echost == 'CHARMM>' .or. echost == 'MSCALE> ' &
              .or. echost == 'BLOCK> ')

  if(qmonoscript)then
     ! one input CHARMM script only for all independent groups:
     ! each mynod=0 has iolev>0, so we need the test on mynodg here!
     if(mynodg == 0) read(unit,'(a)',end=9) card
  else
     ! we are in the streams or no QREPD: each mynod=0 has to do it:
     if(iolev > 0) read(unit,'(a)',end=9) card
  endif
#else /**/
#if KEY_MULTICOM==1 /*  VO string v */
  if (MPI_COMM_PARSER.eq.MPI_COMM_NULL) return ! should not be necessary, but added for clarity
  if ((ME_PARSER.eq.0).or.(iolev.gt.0)) then   ! can do without this by setting iolev>=0 where MPI_PARSER=0, but kept for clarity
#else /*  VO string end ^ */
  if(iolev > 0) then
#endif
   read(unit,'(a)',end=9) card
  endif
#endif /* repdstr */
  cardln=len_trim(card)

#if KEY_REPDSTR==1
  ! this seems to be OK, because the PSNDC after this one below does it
  ! for the stream case. We need special case only for non stream
  if(qmonoscript)then
     !     this has to be global!
     call psetglob
     call psndc(card,1)
     call psetloc
  endif
#endif 
!mfc  #if KEY_REPDSTR2==1
!mfc    ! this seems to be OK, because the PSNDC after this one below does it
!mfc    ! for the stream case. We need special case only for non stream
!mfc    if(qmonoscript)then
!mfc       !     this has to be global!
!mfc       call repd_reps_barrier(ierr)
!mfc       if(lmasternode)then
!mfc          call mpi_comm_rank(comm_master,i,ierr)
!mfc          call mpi_comm_size(comm_master,i,ierr)
!mfc          call mpi_bcast(cardln,1,mpi_integer,0,comm_master,ierr)
!mfc          call mpi_bcast(card,cardln,mpi_character,0,comm_master,ierr)
!mfc       endif
!mfc    endif
!mfc  #endif 

#if KEY_PARALLEL==1
  if(qcmprs) then
#if KEY_ENSEMBLE==1
     call ensprint("RDCMND broadcast"," ")          
#endif
     ! VO string v
#if KEY_MULTICOM==1
     CALL MPI_BCAST(CARD,MXCARD,MPI_BYTE,0,MPI_COMM_PARSER,J) 
#else
     CALL PSNDC(CARD,1)                                       
#endif
     ! VO string ^
#if KEY_ENSEMBLE==1
     call ensprint("RDCMND command",card)           
#endif
     if(card == 'END-OF-FILE') then
        eof=.true.
        return
     endif
  endif
#endif 

  cardln=mxcard
  call trime(card,cardln)
  if(cardln == 0) cardln=1
  if(qprint.and.prnlev >= 3 &
#if KEY_MULTICOM==1 /*  VO stringm : conditional evaluation in parallel */
 &          .and. peek_if() & 
#endif
 &                          ) write(outu,200) echost,card(1:cardln)
200 format(1x,a8,3x,a)

  iend = index(card(1:cardln), exclmk(1:1))
  if (iend == 1) goto 1
  if (iend /= 0) then
     cardln=iend-1
     call trime(card,cardln)
  endif

#if KEY_CFF==1
  if (ucase) &     
#endif
       call cnvtuc(card,cardln)
  if(qcmprs) call cmprst(card,cardln)
  if(cardln == 0) goto 2
  if(card(cardln:cardln) == hyphen.and.usehyp) then
     if(cardln == 1) goto 1
     if(comlen+cardln-1 > mxcms2) then
        CALL WRNDIE(-1,'<RDCMND>','Command line too long: truncated.')
     ENDIF
     CALL ADDST(COMLYN,MXCMS2,COMLEN,CARD,CARDLN-1)
     GOTO 1
  ENDIF
  !
2 CONTINUE
  IF(COMLEN+CARDLN > MXCMS2) THEN
     CALL WRNDIE(-1,'<RDCMND>','Command line too long: truncated.')
  ENDIF
  CALL ADDST(COMLYN,MXCMS2,COMLEN,CARD,CARDLN)
  !
  !     Before returning the string make any parameter
  !     substitutions that may be required.
  !
  ! 2000 FORMAT(' ****  WARNING  **** RDCMND Invalid',
  !     $       ' parameter subsitution. Command ignored.')
  ! 2010 FORMAT(' RDCMND substituted parameter ',I2,
  !     $       ' (',A1,') : "',80A1)
  !
  IF(QCMPRS) FFOUR=COMLYN(1:4)
#if KEY_MULTICOM==1 /* VO stringm conditional execution in parallel */
  if (peek_if()) & 
#endif
  CALL PARSE1(COMLYN,MXCMS2,COMLEN,QPRINT)
  IF(COMLEN == 0)RETURN
  !
  !     Before returning the string make any energy
  !     substitutions that may be required.
  !
2020 FORMAT(' RDCMND substituted energy or value "',80A1)
2030 FORMAT(' RDCMND: can not substitute energy "',80A1)

  ipar = index(comlyn(1:(comlen - 1)), sques(1:1))
  IF (IPAR > 0 &
#if KEY_MULTICOM==1 /*  VO stringm : conditional evaluation in parallel */
&     .and. peek_if() & 
#endif
&              ) THEN
400  CONTINUE
     CALL COPSUB(ENRST,ENRMAX,ENRLEN,COMLYN,IPAR+1, &
          MIN(COMLEN,ENRMAX+IPAR))
     CALL SUBENR(WDLEN,ENRLEN,ENRST,ENRMAX)
     IF (ENRLEN > 0) THEN
        IF (QPRINT .AND. PRNLEV >= 3) THEN
           WRITE(OUTU,2020) (COMLYN(J:J),J=IPAR,IPAR+WDLEN),SDBLQ, &
                ' ','t','o',' ',SDBLQ,(ENRST(J:J),J=1,ENRLEN),SDBLQ
        ENDIF
        CALL COPSUB(SCRTCH,SCRMAX,SCRLEN,COMLYN,IPAR+WDLEN+1,COMLEN)
        COMLEN=IPAR-1
        CALL ADDST(COMLYN,MXCMS2,COMLEN,ENRST,ENRLEN)
        CALL ADDST(COMLYN,MXCMS2,COMLEN,SCRTCH,SCRLEN)
        ipar = index(comlyn(1:(comlen - 1)), sques(1:1))
     ELSE
        !  WE WANT THE WHOLE PARAMETER NAME
        IF(WRNLEV >= 2) WRITE(OUTU,2030) &
             (COMLYN(J:J),J=IPAR,IPAR+WDLEN),SDBLQ
        IPAR=0
     ENDIF
     IF (IPAR > 0) GOTO 400
  ENDIF

  CALL TRIMA(COMLYN,COMLEN)
  RETURN
9 EOF=.TRUE.

#if KEY_PARALLEL==1
  IF(QCMPRS) THEN
     CARD='END-OF-FILE'
#if KEY_REPDSTR==1 
     IF(QREPDSTR.AND.(.NOT.QRDQTT))CALL PSETGLOB    
#endif
     cardln = len_trim(card)
#if KEY_ENSEMBLE==1
     call ensprint("RDCMND broadcast"," ")     
#endif
     ! VO string v
#if KEY_MULTICOM==0
     CALL PSNDC(CARD,1)                                       
#endif
#if KEY_MULTICOM==1
     CALL MPI_BCAST(CARD,MXCARD,MPI_BYTE,0,MPI_COMM_PARSER,J) 
#endif
     ! VO string ^
#if KEY_ENSEMBLE==1
     call ensprint("RDCMND command",card)     
#endif
#if KEY_REPDSTR==1 
     IF(QREPDSTR.AND.(.NOT.QRDQTT))CALL PSETLOC     
#endif
  ENDIF
#endif 
  !
  RETURN
END SUBROUTINE RDCMND

SUBROUTINE SUBENR(WDLEN,ENRLEN,ENRST,ENRMAX)
  !-----------------------------------------------------------------------
  !     THIS ROUTINE GENERATES AN ENERGY STRING SUBSTITUTION
  !
  !     This routine is to remain internal to STRING.FLX  - BRB
  !
  !      Bernard R. Brooks   9/9/83
  !
  use chm_kinds
  use exfunc
  use energym
  use param_store
  use rndnum
  use clcg_mod,only:ranumb
  use string
  !
  implicit none
  !
  INTEGER WDLEN,ENRLEN,ENRMAX
  character(len=*) ENRST
  character(len=8) WRD
  INTEGER I
  !
  real(chm_real) R
  INTEGER IVAL
  logical :: found
  found = .false.
  !
  WRD=NEXTA8(ENRST,ENRLEN)
  IF(WRD == ' ') RETURN
  WDLEN=8
  CALL TRIME(WRD,WDLEN)
  !
  !     Do energy value substitutions
  ENRLEN=0
  DO I=1,LENENP
     IF(WRD == CEPROP(I)) THEN
        R=EPROP(I)
        CALL ENCODF(R,ENRST,ENRMAX,ENRLEN)
        RETURN
     ENDIF
  ENDDO
  DO I=1,LENENT
     IF(WRD == CETERM(I)) THEN
        R=ETERM(I)
        CALL ENCODF(R,ENRST,ENRMAX,ENRLEN)
        RETURN
     ENDIF
  ENDDO
  DO I=1,LENENV
     IF(WRD == CEPRSS(I)) THEN
        R=EPRESS(I)
        CALL ENCODF(R,ENRST,ENRMAX,ENRLEN)
        RETURN
     ENDIF
  ENDDO

  !     Do random number substitution
  IF(WRD == 'RAND' .OR. WRD.EQ.'RANDOM') THEN
     R=RANUMB()
     CALL ENCODF(R,ENRST,ENRMAX,ENRLEN)
     RETURN
  ENDIF

  !     Do random number substitution
  IF(WRD == 'ISEE' .OR. WRD.EQ.'ISEED') THEN
     IVAL = IRNDSD
     CALL ENCODI(IVAL,ENRST,ENRMAX,ENRLEN)
     RETURN
  ENDIF

  !     Do miscellaneous real substitutions
  call find_param(wrd, R, found)
  IF (found) THEN
    CALL ENCODF(R,ENRST,ENRMAX,ENRLEN)
    RETURN
  END IF

  !     Do miscellaneous integer substitutions
  call find_param(wrd, ival, found)
  IF (found) THEN
    CALL ENCODI(ival, ENRST, ENRMAX, ENRLEN)
    RETURN
  END IF

  !     Do miscellaneous character substitutions
  call find_param(wrd, enrst, found)
  if (found) then
    ENRLEN = 8
    CALL TRIME(ENRST, ENRLEN)
    RETURN
  end if

  RETURN
END SUBROUTINE SUBENR

SUBROUTINE XTRANE(ST,STLEN,IDST)
  !-----------------------------------------------------------------------
  !     Checks to see if ST contains extraneous characters, ie. sees if it
  !     non-blank and prints a message including the character string,
  !     IDST, if it does.
  !
  !      Author: Robert Bruccoleri
  !
  use chm_kinds
  use stream
  use string
  implicit none
  character(len=*) ST
  INTEGER STLEN
  character(len=*) IDST
  !
  CALL TRIMA(ST,STLEN)
  IF(STLEN > 0) THEN
     IF(WRNLEV >= 2) THEN
        if (prnlev >= 2) WRITE(OUTU,90) IDST
        if (prnlev >= 2) CALL PRNTST(OUTU,ST,STLEN,1,80)
     ENDIF
     STLEN=0
  ENDIF
90 FORMAT(' **** Warning ****  The following extraneous characters', &
       /,' were found while command processing in ',A)
  !
  RETURN
END SUBROUTINE XTRANE
