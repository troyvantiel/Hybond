SUBROUTINE MINMIZ(COMLYN,COMLEN)
  !-----------------------------------------------------------------------
  !     MINMIZ controls the minimization options.
  !
  !     The input arguments are:
  !
  !     COMLYN           - The command line.
  !     COMLEN           - The length of the command line.
  !
  use abnerm,only:abner

#if KEY_CHEQ==1
  use cheq,only:qcg,cgmodel,qcginv,qcginvf,qpbeq1,   & 
#endif
#if KEY_CHEQ==1
     minnorm,qcgmine,qnoco,qpolar1,ipolar1,     & 
#endif
#if KEY_CHEQ==1
     checketa,checkqnorm,qpartbin,allocate_cheq   
#endif

  use chm_kinds
  use dimens_fcm
  use bases_fcm
  use contrl
  use coord
  use euler
  use eutil
  use hbondm
  use image
  use nraph_m
  use stream
  use string
  use parallel, only: numnod
  use pert
  use pert_mod
  use reawri
  use pathm
  use powell_mod, only: powell
  use psf, only: natom, ngrp
#if KEY_TMD==1
  use tmd,only:inrt   
#endif
#if KEY_DOMDEC==1
  use domdec_common, only: q_domdec  
#endif
#if KEY_DHDGB==1
!AP/MF
  use dhdgb
#endif

  implicit none
  !
  !     Passed variables.
  !
  character(len=*) COMLYN
  INTEGER       COMLEN
  !
  !     Local variables.
  !
  character(len=4) MINOPT
  LOGICAL QIMCEN
  !
#if KEY_CHEQ==1
  LOGICAL QCHEQRDPRM,QCHEQNORM
#endif /* */
  !
#if KEY_TNPACK==1
  !yw...TNPACK, 28-Jul-95
  LOGICAL QSTART
  !
  QLOC1  = .FALSE.
#endif 
  !
  !
  MINXYZ = .TRUE.
  LMINUC = .FALSE.
  QIMCEN = LIMCEN


#if KEY_CHEQ==1
  QCGMIN=(INDXA(COMLYN,COMLEN,'CHEQ') > 0)
  IF (.not.QCG .and. QCGMIN) CALL WRNDIE(-3,'<MINMIZ>', &
       'Fluctuating charges not set-up, use CHEQ ON command first!')
  IF (QCG.AND.QCGMIN)  THEN 
     if(.not.allocated(qpartbin)) then
        call wrndie(-1,'<minmiz> CHEQ not set-up')
     elseif(natim>natom) then
        call allocate_cheq(natim,ngrp)
     endif
     CALL CHECKETA(QCHEQRDPRM)
     IF (.not.QCHEQRDPRM) THEN
        if(prnlev > 1)write(outu,'(2a)') &
             'CHEQ DYNAMICS HAS BEEN REQUESTED BUT', &
             ' CORRECT PARAMETERS HAVE NOT BEEN READ'
        CALL WRNDIE(-1,'<MINMIZ>', &
             'CHEQ PARAMETERS HAVE NOT BEEN READ')
     ENDIF
     CALL CHECKQNORM(QCHEQNORM) 
     IF (.not.QCHEQNORM) THEN
        if(prnlev > 1)write(outu,'(2a)') &
             'CHEQ MINIMIZATION HAS BEEN REQUESTED BUT', &
             ' CORRECT NORMALIZATION SPECIFICATIONS LACKING'
        CALL WRNDIE(-1,'<MINMIZ>', &
             'CHEQ NORMALIZATION ASSIGNMENT INCORRECT')
     ENDIF
     !  ----   SINCE USER WANTS CHEQ, MAKE SURE NORMALIZATION IS SET UP AND CORRECT

     CGMODEL=GTRMI(COMLYN,COMLEN,'CGMD',CGMODEL)
     QCGINV = INDXA(COMLYN, COMLEN, 'CGIN')  >  0
     QCGINVF = INDXA(COMLYN, COMLEN, 'CGFC')  >  0
     QPBEQ1 = INDXA(COMLYN, COMLEN, 'PBEQ')  >  0

     !   FOR CHARGE NORMALIZATION
     MINNORM=.TRUE.   ! NORMALIZE w/out using different masses for minimization
     !  
     IF (PRNLEV > 0) THEN
        IF (QCGMIN) WRITE(OUTU,'(a)')"CHEQ HAS BEEN READ"
        IF (QPBEQ1) WRITE(OUTU,'(a)')"PBEQ is requested"
        IF (QCGINV) WRITE(OUTU,'(a)')"<MINI>: CG INVERSION REQUESTED"
        IF (QCGINVF) WRITE(OUTU,'(a)')"<MINI>: CG FORCE REQUESTED"
     ENDIF
     QCGMINE=QCGMIN
     MINXYZ=INDXA(COMLYN, COMLEN, 'NOCO') <= 0
     QNOCO=.NOT.MINXYZ
     QPOLAR1=INDXA(COMLYN, COMLEN, 'QPOL')  >  0
     IPOLAR1=GTRMI(COMLYN,COMLEN,'IPOL',IPOLAR1)
  ENDIF
#endif 
  !
  MINOPT = NEXTA4(COMLYN,COMLEN)
#if KEY_DOMDEC==1
  if (q_domdec .and. .not. ( minopt == 'ABNR' .or. minopt == 'SD  ' ) ) &
       call wrndie(-1, '<MINMIZ>', &
           'Cannot minimize after enabling DOMDEC.')

  if (q_domdec .and. (numnod .gt. 1)) &
       call wrndie(-1, '<MINMIZ>', &
           'Cannot minimize with more than one MPI process after enabling DOMDEC.')
#endif

  !
  !     Process update commands.
  !
  CALL UPDATE(COMLYN,COMLEN,X,Y,Z,WMAIN,.TRUE., &
       .TRUE.,.TRUE.,.TRUE.,.TRUE.,0,0,0,0,0,0,0)
  CALL FINCYC(NUPFRQ,0,0,0,0,INBFRQ,IHBFRQ,0,IMGFRQ,0,0,0 &
#if KEY_TMD==1
       ,inrt &  
#endif
       )
  !
  ! turn off image centering during minimization (but restore it later)
  LIMCEN = .FALSE.
  !
  !-----------------------------------------------------------------------
  ! Parse general minimization options.
  !
  !     Parse the saddle code option (i.e. minimize the gradient**2)
  QSADLE=(INDXA(COMLYN,COMLEN,'GRAD') > 0)
  IF(QSADLE .AND. PRNLEV >= 2) WRITE(OUTU,255)
255 FORMAT(' CHARMM> Energy will be the mean squared gradient', &
       ' during minimizations.')
  !
  !     Parse the numerical derivatives option.
  QNUMER = INDXA(COMLYN,COMLEN,'NUME')  >  0
  IF(QNUMER .AND. PRNLEV >= 2) WRITE(OUTU,256)
256 FORMAT(' CHARMM> Forces will be determined by finite differences', &
       ' during minimizations.')
  !
  ! Parse trajectory options
  IUNCRD = GTRMI(COMLYN,COMLEN,'IUNC',-1)
  NSAVC  = GTRMI(COMLYN,COMLEN,'NSAVC',1)
  IUNXYZ = GTRMI(COMLYN,COMLEN,'IUNX',-1)
  NSAVX  = GTRMI(COMLYN,COMLEN,'NSAVX',1)
  MXYZ   = GTRMI(COMLYN,COMLEN,'MXYZ',1)
  !
  ! SAPATEL
#if KEY_CHEQ==1
  !      IF ( .not. (MINOPT  ==  'CGSD' .OR. MINOPT .EQ. 'CONJ'
  !     $    .or. MINOPT  ==  'SD  ') .and. QCGMIN )
  !     $   CALL WRNDIE(-3,'<MINMIZ>',
  !     $   'CHEQ only supported for CONJ and SD Minimizers')
  IF ( QSADLE .AND. QCG ) CALL WRNDIE(-3,'<MINMIZ>', &
       'SADLE option not currently supported with CHEQ.')
#endif 
  ! SAPATEL
  !-----------------------------------------------------------------------
#if KEY_PERT==1
  IF(QPERT) CALL PERTDF
#endif 
  !
  !     Branch on the minimization option.
  !
  IF (MINOPT  ==  'ABNR') THEN
     CALL ABNER(COMLYN,COMLEN)
     !
  ELSE IF (MINOPT  ==  'TN  ') THEN
     !yw...TNPACK: updated 28-Jul-95 and 12-Aug-95
#if KEY_TNPACK==1
     QSTART=.TRUE.
     CALL TNDRIV(COMLYN,COMLEN,QSTART)
#else /**/
     CALL WRNDIE(-3,'<MINMIZ>','TN minmizer NOT compiled.')
#endif 
     !
  ELSE IF (MINOPT  ==  'CGSD' .OR. MINOPT .EQ. 'CONJ') THEN
     CALL CONJUG(COMLYN,COMLEN)
     !
  ELSE IF (MINOPT  ==  'NRAP') THEN
     !CC      IF (NTRANS  >  0) CALL XNBLST
     CALL NRAPH(COMLYN,COMLEN,X,Y,Z)
     !
  ELSE IF (MINOPT  ==  'POWE') THEN
     ! GAMESS has it too
!--mfc-- ##IF GAMESS
!--mfc--      CALL POWELLC(COMLYN,COMLEN,X,Y,Z)
!--mfc-- ##ELSE
     CALL POWELL(COMLYN,COMLEN,X,Y,Z)
!--mfc-- ##ENDIF
     !
  ELSE IF (MINOPT  ==  'SD  ') THEN
#if KEY_RPATH==1 /*rpath*/
     IF (QPNEB) THEN 
        CALL STEEPDNEB(COMLYN,COMLEN)
        !         write(*,*)'IM IN STEEPDNEB'
     ELSE
#endif /*    (rpath) */
        CALL STEEPD(COMLYN,COMLEN)
        !         write(*,*)'IM IN NORMAL STEEPD'
#if KEY_RPATH==1
     ENDIF     
#endif
     !
  ELSE
     CALL WRNDIE(-3,'<MINMIZ>', &
          'Unrecognised minimization procedure.')
  ENDIF
  !
  MINXYZ = .TRUE.
  LMINUC = .FALSE.
  LIMCEN = QIMCEN
  !
#if KEY_PERT==1
  IF(QPERT) CALL PERTAN(.FALSE.)
#endif 
  ! SAPATEL
#if KEY_CHEQ==1
  QCGINV=.FALSE.
  QCGINVF=.FALSE.
  QCGMINE=.FALSE.
  QPBEQ1 =.FALSE.
  QNOCO =.FALSE.
  CGMODEL=0
  QPOLAR1=.FALSE.
  IPOLAR1=-1
  MINNORM=.FALSE.
#endif 
  ! SAPATEL
  !
  RETURN
END SUBROUTINE MINMIZ

SUBROUTINE MINTRJ(NCALLS,NSTEPS,STEP)
  !-----------------------------------------------------------------------
  !     Writes trajectory frame of minimization steps
  !
#if KEY_CHEQ==1
  use cheq,only:qcg     
#endif

  use chm_kinds
  use dimens_fcm
  use reawri
  use psf
  use coord
  use cvio
  use deriv
  use dynio
  use number
  use ctitla
  !
  implicit none
  !
  INTEGER NCALLS,NSTEPS
  real(chm_real) STEP,T(1)
  !
  INTEGER NAT3
  !
  T(1)=ZERO
  IF((NSAVC > 0).AND.(IUNCRD >= 0) &
       .AND.(MOD(NCALLS,NSAVC) == 0)) THEN
     NAT3=3*NATOM
     CALL WRITCV(X,Y,Z, &
#if KEY_CHEQ==1
          CG,QCG,                                & 
#endif
          NATOM, (/ 0 /), NATOM,NCALLS,NCALLS,NAT3,STEP,NSAVC, &
          NSTEPS,TITLEA,NTITLA,IUNCRD,.FALSE.,.FALSE., (/ 0 /), .FALSE., (/ ZERO /))
  ENDIF
  IF((NSAVX > 0).AND.(IUNXYZ >= 0).AND.(MOD(NCALLS,NSAVX) == 0)) &
       CALL WRXYZ(IUNXYZ,X,Y,Z,T,T,T,DX,DY,DZ,NCALLS)
  !
  RETURN
END SUBROUTINE MINTRJ

