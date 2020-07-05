module genpsf_m
  use chm_kinds
  implicit none

contains

SUBROUTINE GENPSF(COMLYN,COMLEN,ISTART)
  !
  !     Generate the PSF.
  !
  !     Jay Banks 23 Oct 95: add MMFF codes
  !     Rick Lapp 25 May 98: add CFF code
  !

#if KEY_CHEQ==1
  use cheq,only:qcgrp,molbl,molsrt,allocate_cheq  
#endif
  use intcor_module
  use dimens_fcm
  use bases_fcm
  use memory
  use psf
  use rtf,only:autod,autot
  use stream
  use string
  use ffieldm
  use io
  use lonepr
  use aniso_fcm
#if KEY_REPLICA==1
  use replica_mod         
#endif
  use number

  integer,allocatable,dimension(:,:) :: IATBON
  integer,allocatable,dimension(:) :: I_sv
  integer,allocatable,dimension(:) :: NATBON
#if KEY_CHEQ==1
  integer,allocatable,dimension(:,:) :: IWORK1
#endif 
  ! . Passed variables.
  CHARACTER(len=*) COMLYN
  INTEGER   COMLEN, ISTART
  ! . Local variables.
  CHARACTER(len=8) PATF, PATL, WRD            ! yw A4->A8
  CHARACTER(len=8) WRDTMP
  INTEGER     WTLEN
  INTEGER   I, IDUMMY, ISTOP, NBONDL, NSEGM
  LOGICAL   CMPLTD, LANGLE, LPHI, LSETIC, LWARN
  !
  LOGICAL   LDRUDE, LSHOW
  real(chm_real)    DMASS
  INTEGER   NATOLD
  !
  WRD = GTRMA(COMLYN, COMLEN, 'DUPL')
#if KEY_CHEQ==1
  ! Default is to use real molecule as molecules
  ! has to specify to use group as molecules
  QCGRP = INDXA(COMLYN, COMLEN, 'GROU').GT.0
  QCGRP = (INDXA(COMLYN, COMLEN, 'GRP ').GT.0).OR.QCGRP
#endif 
  IF(WRD .EQ. ' ') THEN
     ! Generate a segment from a sequence (or Merck file).
     !
#if KEY_MMFF==1
     LDRUDE = .false.
     if(DATA_READ.eq.SEQUENCE) then
#endif 
        NSEGM  = NSEG
        NSEG   = NSEG + 1
        NBONDL = NBOND + 1
        LSETIC = INDXA(COMLYN, COMLEN, 'SETU') .GT. 0
        LWARN  = INDXA(COMLYN, COMLEN, 'WARN') .GT. 0
        !
        PATF = GTRMA(COMLYN, COMLEN, 'FIRS')
        IF (PATF .EQ. ' ') PATF = 'DEFA'
        PATL = GTRMA(COMLYN, COMLEN, 'LAST')
        IF (PATL .EQ. ' ') PATL = 'DEFA'
        !yw A4 -> A8
        SEGID(NSEG) = NEXTA8(COMLYN, COMLEN)
        IF (SEGID(NSEG) .EQ. ' ') &
             CALL ENCODI(NSEG, SEGID(NSEG), 8, IDUMMY)
        !yw--
        ! . Check for duplicate segment names.
        DO I = 1,NSEGM
           IF(SEGID(I) .EQ. SEGID(NSEG)) THEN
              NSEG = NSEGM
              CALL WRNDIE(-2,'<CHARMM>','DUPLICATE SEGMENT NAMES')
              RETURN
           ENDIF
        ENDDO
        !
        LDRUDE = (INDXA(COMLYN, COMLEN, 'DRUD') .GT. 0)
        LSHOW  = (INDXA(COMLYN, COMLEN, 'SHOW') .GT. 0)
        IF(LDRUDE)THEN
           !
           ! Drude uses lonepair facility. Make sure arrays are
           ! allocated, if not, allocate them. cb3
#if KEY_LONEPAIR==1
           if(.not.allocated(lpnhost)) call allocate_lonepair  
#endif
           !
           ! Drude also needs aniso data structures, check if present
           ! and allocate if not. cb3
           if(.not.allocated(lstani1)) call allocate_aniso
           !
           !
           DMASS  = GTRMF(COMLYN,COMLEN,'DMAS',ZERO)
           if (prnlev > 2) write(outu,'(1x,3A,f10.4)')  &
                'Drude polarizability will be setup for SEGID: ', &
                SEGID(NSEG), &
                'mass of Drudes particles =',DMASS

           ! Haibo Yu Hyper-Polarizability
           QHYPER=(INDXA(COMLYN,COMLEN,'HYPE').GT.0)
           IF(QHYPER)THEN
              HORDER = GTRMI(COMLYN,COMLEN,'HORD',0)
              KHYPER = GTRMF(COMLYN,COMLEN,'KHYP',ZERO)
              RHYPER = GTRMF(COMLYN,COMLEN,'RHYP',ZERO)
              IF(PRNLEV.GT.2)THEN  
                 WRITE(OUTU,'(A/)') ' DRUDE>HyperPolar included'
              ENDIF
           ENDIF
        ENDIF
        !
        ISTART = NICTOT(NSEG) + 1
        ISTOP  = NRES
        call chmalloc('genpsf.src','GENPSF','I_sv',MAXRES,intg=I_sv)
        NATOLD = NATOM
        CALL GENIC(ISTART, ISTOP, LWARN, LSETIC, PATF, PATL, I_sv, &
             LDRUDE)
        call chmdealloc('genpsf.src','GENPSF','I_sv',MAXRES,intg=I_sv)
#if KEY_MMFF==1
        ! . Set (needed for MMFF) and print out the structure file counters.
        if(FFIELD.eq.MMFF) CALL PSFSUM(OUTU)
#endif 
        ! . Check for auto generation requests.
        LANGLE = AUTOT
        LPHI   = AUTOD
        IF (INDXA(COMLYN, COMLEN, 'ANGL') .GT. 0) LANGLE = .TRUE.
        IF (INDXA(COMLYN, COMLEN, 'NOAN') .GT. 0) LANGLE = .FALSE.
        IF (INDXA(COMLYN, COMLEN, 'DIHE') .GT. 0) LPHI = .TRUE.
        IF (INDXA(COMLYN, COMLEN, 'NODI') .GT. 0) LPHI = .FALSE.
        IF (LANGLE .OR. LPHI) THEN
           call chmalloc('genpsf.src','GENPSF','NATBON',NATOM,intg=NATBON)
           call chmalloc('genpsf.src','GENPSF','IATBON',IATBMX,NATOM,intg=IATBON)
           CALL AUTGEN(NBONDL,NATBON,IATBON,LANGLE,LPHI, &
                LDRUDE)
           call chmdealloc('genpsf.src','GENPSF','NATBON',NATOM,intg=NATBON)
           call chmdealloc('genpsf.src','GENPSF','IATBON',IATBMX,NATOM,intg=IATBON)
        ENDIF
#if KEY_MMFF==1
     ENDIF ! if(DATA_READ.eq.SEQUENCE) then
     IF(FFIELD.EQ.MMFF) CALL MMFF_SETUP
#endif 
     !
     IF(LDRUDE)THEN
        QDRUDE = .TRUE.
        CALL MKDRUD(NATOLD,NATOM,ATYPE,IAC,CG,AMASS,DMASS, &
             THOLEI,ALPHADP,ISDRUDE, &
             nbond,ib,jb,NDRUDE,LSHOW)
        nbdrude=nbond
     ENDIF

  ELSE
     ! . Duplicate an existing segment.
     NSEGM  = NSEG
     NSEG   = NSEG + 1
     LSETIC = INDXA(COMLYN, COMLEN, 'SETU') .GT. 0
     !yw++ A4->A8
     SEGID(NSEG) = NEXTA8(COMLYN, COMLEN)
     IF (SEGID(NSEG) .EQ. ' ') &
          CALL ENCODI(NSEG, SEGID(NSEG), 8, IDUMMY)
     !yw--
     ! . Check for duplicate segment names.
     DO I = 1,NSEGM
        IF(SEGID(I) .EQ. SEGID(NSEG)) THEN
           NSEG = NSEGM
           CALL WRNDIE(-2,'<CHARMM>','DUPLICATE SEGMENT NAMES')
           RETURN
        ENDIF
     ENDDO
     CALL DUPSEG(LSETIC, WRD, CMPLTD)
     IF(.NOT.CMPLTD) THEN
        NSEG = NSEGM
        CALL WRNDIE(0,'<CHARMM>','Could not duplicate segment')
        RETURN
     ENDIF
  ENDIF
  !
#if KEY_CHEQ==1
#if KEY_REPLICA==1
  IF(.NOT.QREP) THEN      
#endif
     !  find molecule label
     if(allocated(molbl)) then
        ! psf changed, check whether arrays are still right size
        call allocate_cheq(natom, ngrp)
        CALL MOLSRT(NATOM,IB,JB,NBOND,QCGRP,NGRP,IGPBS)
        IF (QCGRP.AND.PRNLEV.GE.2) WRITE(OUTU,'(A,/,A)') &
             'GENPSF> Treat GRP as Molecules in CG dynamics.', &
             '  This is true for the whole PSF.'
     endif
#if KEY_REPLICA==1
  ENDIF                   
#endif
#endif 

  WRDTMP=SEGID(NSEG)
  WTLEN=LEN(WRDTMP)
  CALL TRIMA(WRDTMP,WTLEN)
  IF(PRNLEV.GE.2) WRITE (OUTU,'(A,I3,A,A,A)')   & ! yw A4->A
       ' GENPSF> Segment ', NSEG, &
       ' has been generated. Its identifier is ',WRDTMP(1:WTLEN), '.'
  ! . Print out the structure file counters.
  CALL PSFSUM(OUTU)
  !
  NBDRUDE=NBOND
  RETURN
END SUBROUTINE GENPSF

SUBROUTINE GENIC(ISTART,ISTOP,LWARN,LSETIC,PATFX,PATLX,ITPL, &
     LDRUDE)
  !
  !     GENERATES A SEGMENT OF A PROTEIN STRUCTURE FILE.
  !
  !     WRITTEN BY B. R. BROOKS  APRIL-1982
  !
#if KEY_CHEQ==1
  use cheq, only:ech,eha           
#endif
  use intcor_module
  use dimens_fcm
  use number
  use consta
  use bases_fcm
  use psf
  use rtf, only: nrtrs, nic, ftp, mac, armass, chg, alph, thol, &
    grpr, aa, rtrtyp, mib, mjb, mit, mjt, mkt, mip, mjp, mkp, &
    mlp, mim, mjm, mkm, mlm, &
    mlp1ct, mlp2ct, mlp3ct, mlp4ct, mlp0ct, mlpct1d,rlpct, tlpct, plpct, mianis, &
    maxcent_hosts, &
    mjanis, mkanis, mlanis, a22anis, a11anis, mh, md, ma, maa, bari, barj, &
    bark, barl, bart, aaf, aal, delat, mxn, mnb, icb1, icb2, &
    icth1, icth2, icphi

#if KEY_MMFF==1
  use rtf, only: atnumr, mbtype   
#endif
#if KEY_CHEQ==1
  use rtf, only: echh, ehah       
#endif
#if KEY_CMAP==1
  use rtf, only: mi1ct, mj1ct, mk1ct, ml1ct, mi2ct,mj2ct, mk2ct, ml2ct 
#endif

  use stream
  use string
  use lonepr
  use aniso_fcm
#if KEY_MMFF==1 || KEY_CFF==1
  use ffieldm
  use mmffm
#endif 
#if KEY_CFF==1
  use cff_fcm
#endif 
  use machutil,only:die

  INTEGER ISTART,ISTOP
  LOGICAL LWARN,LSETIC
  CHARACTER(len=*) PATFX,PATLX
  INTEGER ITPL(*)
  LOGICAL LDRUDE
  !
  !
  real(chm_real) QTOT
  real(chm_real) A33
  INTEGER NATOLD,IPATCH,I,ITEMP,IPATF,IPATL
  INTEGER NGRPFS,IGRLST,IRES,K,NATP,NATPT,NAT,IPAT,IGRP
  INTEGER II,JJ,IPT,N,NL,KK,LL,J,JK,NATPT2,NAT2,I2,IMATCH,JMATCH
  INTEGER, DIMENSION(maxcent_hosts):: IIJ  ! The dimension is equal to maxcent_hosts: Maximum number of hosts for defining lone pair center.
#if KEY_CMAP==1
  INTEGER II1,JJ1,KK1,LL1,II2,JJ2,KK2,LL2
#endif 
  INTEGER IRSZ
  INTEGER I400,I500,I600,I700
  CHARACTER(len=8) PATF,PATL,TEMP
  INTEGER WTLEN
  LOGICAL FOUND,OK,LW,LADDED
  INTEGER ID1
  CHARACTER(len=1) ID2
  CHARACTER(len=8) :: NONE='NONE',BLANK='    ',DEFA='DEFA'
  INTEGER, parameter :: MARK = -99999999
  !
  ! Define ITPL array  -- type(number) of residue
  LW=.FALSE.
  !
  NATOLD=NATOM
  !
  IPATCH=0
  I=ISTART
50 CONTINUE
  TEMP=RES(I)
  !        find-this-residue
  call sub500
505 CONTINUE
  ITPL(I)=ITEMP
  I=I+1
  IF(I.LE.ISTOP) GOTO 50
  !
  ! Find first and last patch residues
  !
243 FORMAT(' NO PATCHING WILL BE DONE ON THE FIRST RESIDUE')
245 FORMAT(' THE PATCH ''',A,''' WILL BE USED FOR THE FIRST RESIDUE')
253 FORMAT(' NO PATCHING WILL BE DONE ON THE LAST  RESIDUE')
255 FORMAT(' THE PATCH ''',A,''' WILL BE USED FOR THE LAST  RESIDUE')
  PATF=PATFX
  IPATCH=1
  IF(PATF.EQ.DEFA) PATF=AAF(ITPL(ISTART))
  IF(PATF.EQ.NONE) THEN
     IPATF=0
     IF(PRNLEV.GE.2) WRITE(OUTU,243)
  ELSE
     WTLEN=LEN(PATF)
     CALL TRIME(PATF,WTLEN)
     IF(PRNLEV.GE.2) WRITE(OUTU,245) PATF(1:WTLEN)
     TEMP=PATF
     !        find-this-residue
     call sub500

507  CONTINUE
     IPATF=ITEMP
  ENDIF
  !
  PATL=PATLX
  IF(PATL.EQ.DEFA) PATL=AAL(ITPL(ISTOP))
  IF(PATL.EQ.NONE) THEN
     IPATL=0
     IF(PRNLEV.GE.2) WRITE(OUTU,253)
  ELSE
     WTLEN=LEN(PATL)
     CALL TRIME(PATL,WTLEN)
     IF(PRNLEV.GE.2) WRITE(OUTU,255) PATL(1:WTLEN)
     TEMP=PATL
     !        find-this-residue
     call sub500

509  CONTINUE
     IPATL=ITEMP
  ENDIF
  !
  ! Generate atom sequence for new segment. added atoms of the
  ! first patch go at the begining. added atoms of the last patch
  ! are added at the end.
  !
  !
  NGRPFS=NGRP+1
  IGRLST=-1
  IRES=ISTART-1
  !     do ires=istart,istop
300 CONTINUE
  IRES=IRES+1
  K=ITPL(IRES)
  NATP=0
  IBASE(IRES)=NATOM
  !
  ! Process initial patch for first residue
  IF(IRES.EQ.ISTART.AND.IPATF.GT.0) THEN
     NATPT=0
     IF(IPATF.GT.1) NATPT=NIC(1,IPATF-1)
     NAT=NIC(1,IPATF)-NATPT
     LADDED=.FALSE.
     I=0
     !           do i=1,nat
320  CONTINUE
     I=I+1
     IF(.NOT.DELAT(I+NATPT)) THEN
        OK=.TRUE.
        IMATCH=-1
        IF(IRES.EQ.ISTOP.AND.IPATL.GT.0) THEN
           IPAT=IPATL
           !                    check-this-atom
           call sub700

705        CONTINUE
        ENDIF
        IF(.NOT.OK) OK=.NOT.DELAT(IMATCH)
        IF(OK) THEN
           !                    add-this-atom
           call sub600

605        CONTINUE
           JMATCH=IMATCH
           !
           ! Check for charge augmenting for first patch
           IF(ABS(CG(NL)).GT. NINETY) THEN
              IF(CG(NL).GT.0.0) THEN
                 CG(NL)=CG(NL)-HUNDRD
              ELSE
                 CG(NL)=CG(NL)+HUNDRD
              ENDIF
              IMATCH=-1
              IPAT=K
              !                       check-this-atom
              call sub700

707           CONTINUE
              IF(IMATCH.GT.0) THEN
                 CG(NL)=CG(NL)+CHG(IMATCH)
              ELSE
                 CALL WRNDIE(-2,'<GENPSF>', &
                      'CHARGE VALUE OUT OF RANGE')
              ENDIF
           ENDIF

           ! Check for more charge augmenting with one residue segments
           IF(JMATCH.GT.0) THEN
              IF(ABS(CHG(JMATCH)).GT. NINETY) THEN
                 IF(CHG(JMATCH).GT.0.0) THEN
                    CG(NL)=CG(NL)-HUNDRD+CHG(JMATCH)
                 ELSE
                    CG(NL)=CG(NL)+HUNDRD+CHG(JMATCH)
                 ENDIF
              ELSE
                 CG(NL)=CHG(JMATCH)
              ENDIF
              IBLO(NL)=JMATCH
              ITEMP=MAC(JMATCH)
              IAC(NL)=ITEMP
              AMASS(NL)=ARMASS(ITEMP)
              IF(LDRUDE)THEN
                 !                       write(*,*) 'JMATCH? ',NL,I+NATPT,JMATCH,
                 !    &                  CG(NL),ALPH(JMATCH),THOL(JMATCH)
                 ALPHADP(NL)=ALPH(JMATCH)
                 THOLEI(NL)=THOL(JMATCH)
              ENDIF
              RSCLF(NL)=1.0
#if KEY_WCA==1
              WCA(NL)=1.0            
#endif
#if KEY_MMFF==1
              !
              !RCZ940119 - this modification should allow to use RTF with MMFF ???
              !
              if(ffield == mmff) AtNum(NL)=AtNumR(JMATCH)! atomic number
#endif 
           ENDIF
           !
        ENDIF
     ENDIF
     IF(I.LT.NAT) GOTO 320
     IF (LADDED) IGRLST=0
  ENDIF
  !
  ! Process main body of residue
  NATPT=0
  IF(K.GT.1) NATPT=NIC(1,K-1)
  NAT=NIC(1,K)-NATPT
  I=0
  !        do i=1,nat
340 CONTINUE
  I=I+1
  OK=.TRUE.
  IMATCH=-1
  IF(IRES.EQ.ISTOP.AND.IPATL.GT.0) THEN
     IPAT=IPATL
     !              check-this-atom
     call sub700

709  CONTINUE
  ENDIF
  IF(.NOT.OK) OK=.NOT.DELAT(IMATCH)
  JMATCH=IMATCH
  IMATCH=-1
  IF(IRES.EQ.ISTART.AND.IPATF.GT.0) THEN
     IPAT=IPATF
     !              check-this-atom
     call sub700

711  CONTINUE
  ENDIF
  IF(IMATCH.GT.0) THEN
     IF(IGRLST.NE.GRPR(I+NATPT)) IGRLST=GRPR(I+NATPT)
  ENDIF
  IF (OK) THEN
     !              add-this-atom
     call sub600

607  CONTINUE
     IF(ABS(CG(NL)).GT.900.0) &
          CALL WRNDIE(-2,'<GENPSF>','CHARGE VALUE OUT OF RANGE')
     !
     ! Check for charge augmenting with last residue patch
     IF(JMATCH.GT.0) THEN
        IF(ABS(CHG(JMATCH)).GT. NINETY) THEN
           IF(CHG(JMATCH).GT.0.0) THEN
              CG(NL)=CG(NL)-HUNDRD+CHG(JMATCH)
           ELSE
              CG(NL)=CG(NL)+HUNDRD+CHG(JMATCH)

           ENDIF
        ELSE
           CG(NL)=CHG(JMATCH)
#if KEY_CHEQ==1
           ECH(NL)=ECHH(JMATCH)
           EHA(NL)=EHAH(JMATCH)
#endif 
        ENDIF
        IBLO(NL)=JMATCH
        ITEMP=MAC(JMATCH)
        IAC(NL)=ITEMP
        AMASS(NL)=ARMASS(ITEMP)
        IF(LDRUDE)THEN
           !                 write(*,*) 'augmenting last?',NL,I+NATPT,JMATCH,
           !    &            CG(NL),ALPH(JMATCH),THOL(JMATCH)
           ALPHADP(NL)=ALPH(JMATCH)
           THOLEI(NL)=THOL(JMATCH)
        ENDIF
        RSCLF(NL)=1.0
#if KEY_WCA==1
        WCA(NL)=1.0               
#endif
#if KEY_MMFF==1
        if(ffield == mmff) AtNum(NL)=AtNumR(JMATCH)  
#endif
     ENDIF
     !
  ENDIF
  IF(I.LT.NAT) GOTO 340
  IF (LADDED) IGRLST=0
  !
  ! Process final patch for last residue
  IF(IRES.EQ.ISTOP.AND.IPATL.GT.0) THEN
     NATPT=0
     IF(IPATL.GT.1) NATPT=NIC(1,IPATL-1)
     NAT=NIC(1,IPATL)-NATPT
     I=0
     !           do i=1,nat
360  CONTINUE
     I=I+1
     OK=.TRUE.
     IF(IRES.EQ.ISTART.AND.IPATF.GT.0) THEN
        IPAT=IPATF
        !                 check-this-atom
        call sub700

713     CONTINUE
     ENDIF
     IPAT=K
     !              check-this-atom
     call sub700

715  CONTINUE
     IF(.NOT.OK .AND. .NOT.DELAT(I+NATPT)) THEN
        IF(IGRLST.NE.GRPR(I+NATPT)) IGRLST=GRPR(I+NATPT)
     ENDIF
     IF(OK .AND. .NOT.DELAT(I+NATPT)) THEN
        !                 add-this-atom
        call sub600

609     CONTINUE
        !
        ! Check for charge augmenting for final patch
        IF(ABS(CG(NL)).GT.NINETY) CALL WRNDIE(-2,'<GENPSF>', &
             'CHARGE VALUE OUT OF RANGE')
        !
     ENDIF
     IF(I.LT.NAT) GOTO 360
  ENDIF
  !
  NATOM=NATOM+NATP
  IF(IRES.LT.ISTOP) GOTO 300
  !
  IBASE(ISTOP+1)=NATOM
  IGPBS(NGRP+1)=NATOM
  !
  ! Find group types
  DO IGRP=NGRPFS,NGRP
     II=IGPBS(IGRP)+1
     JJ=IGPBS(IGRP+1)
#if KEY_NOST2==0
     IF(IGPTYP(IGRP).EQ.0) THEN
#endif 
        QTOT=0.0
        DO I=II,JJ
           QTOT=QTOT+CG(I)
           IF(CG(I).NE.0.0) IGPTYP(IGRP)=1     ! MBA Something seems wrong here.
        ENDDO
!        IF(ABS(QTOT).LE.0.0001) IGPTYP(IGRP)=1  !MBA
        IF(ABS(QTOT).GT.0.0001) IGPTYP(IGRP)=2
#if KEY_NOST2==0
     ELSE
        IF((JJ-II).NE.4) &
             CALL WRNDIE(-4,'<GENIC>','ST2 group must have 5 atoms')
        NST2=NST2+1
        !           set imove for lone pairs
        IMOVE(JJ-1)=-1
        IMOVE(JJ)=-1
     ENDIF
#endif 
  ENDDO
  !
  ! Add in internal coordinates using final atom sequence
  !
  IRES=ISTART
  !     do ires=istart,istop

380 CONTINUE
  IF(IRES.EQ.ISTART.AND.IPATF.GT.0) THEN
     K=IPATF
     !           add-in-psf-elements-for-this-residue
     call sub400


405  CONTINUE
  ENDIF
  K=ITPL(IRES)
  !        add-in-psf-elements-for-this-residue
  call sub400

407 CONTINUE
  IF(IRES.EQ.ISTOP.AND.IPATL.GT.0) THEN
     K=IPATL
     !           add-in-psf-elements-for-this-residue
     call sub400

409  CONTINUE
  ENDIF
  !        put in nonbonded exclusions
  IPT=IBASE(IRES)
  N=IBASE(IRES+1)-IPT
  DO I=1,N
     IPT=IPT+1
     NL=IBLO(IPT)
     II=0
     IF(NL.GT.1) II=MXN(NL-1)
     JJ=MXN(NL)-II
     DO KK=1,JJ
        II=II+1
        CALL PATOM(LL,ID1,ID2,MNB(II),ATYPE,IBASE,[IRES],0,ISTART, &
             ISTOP,LW)
        IF(LL.GT.0) THEN
           NNB=NNB+1
           INB(NNB)=LL
        ELSE
           IF(LWARN .AND. WRNLEV.GE.2 .and. prnlev >= 2) WRITE(OUTU,180) IRES, &
                AA(K),MNB(II),ATYPE(IPT)(1:idleng)
        ENDIF
     ENDDO
     IBLO(IPT)=NNB
  ENDDO
  !
  IRES=IRES+1
  IF(IRES.LE.ISTOP) GOTO 380
  !

  ! Initialize atom coordinates, etc.
  !
  CALL ATMINI(NATOLD+1,NATOM)
  !
  RETURN
  !
180 FORMAT(' ** WARNING ** NONB-EXCLUSION NOT FOUND FOR RESIDUE ', &
       I3,1X,A6/,' EXCLUSION "',A6,'"', &
       ' WAS REQUESTED FOR ATOM "',A,'"')
  !
  !

contains
  !
  !=======================================================================
  ! to add-in-psf-elements-for-this-residue
  subroutine sub400
    use intcor_module
    use intcor2,only:puticel
    logical:: lok
    !
    ! put in bonds

    IPT=0
    IF(K.GT.1) IPT=NIC(2,K-1)
    N=NIC(2,K)-IPT
    DO I=1,N
       IPT=IPT+1
       CALL PATOM(II,ID1,ID2,MIB(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(JJ,ID1,ID2,MJB(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       IF(II.GT.0.AND.JJ.GT.0) THEN
            FOUND = .false.
            DO KK=1,NBOND
              IF(((II.EQ.IB(KK)).AND.(JJ.EQ.JB(KK))) &
                 .OR. &
              ((II.EQ.JB(KK)).AND.(JJ.EQ.IB(KK))))THEN
                 FOUND = .true.
              ENDIF
            ENDDO
            IF (.not. FOUND) THEN
              NBOND=NBOND+1
              IB(NBOND)=II
              JB(NBOND)=JJ
            ENDIF

#if KEY_MMFF==1
            if(ffield == mmff) BondType(NBOND)=MBTYPE(IPT) 
#endif
       ELSE
          IF(LWARN .AND. WRNLEV.GE.2 .and. prnlev >= 2) &
               WRITE(OUTU,140) IRES,AA(K),MIB(IPT),MJB(IPT)
       ENDIF
    ENDDO
    !
    ! put in angles
    IPT=0
    IF(K.GT.1) IPT=NIC(3,K-1)
    N=NIC(3,K)-IPT
    DO I=1,N
       IPT=IPT+1
       CALL PATOM(II,ID1,ID2,MIT(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(JJ,ID1,ID2,MJT(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(KK,ID1,ID2,MKT(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       IF(II.GT.0.AND.JJ.GT.0.AND.KK.GT.0) THEN
          NTHETA=NTHETA+1
          IT(NTHETA)=II
          JT(NTHETA)=JJ
          KT(NTHETA)=KK
       ELSE
          IF(LWARN .AND. WRNLEV.GE.2 .and. prnlev >= 2) &
               WRITE(OUTU,150) IRES,AA(K),MIT(IPT),MJT(IPT),MKT(IPT)
       ENDIF
    ENDDO
    !
    ! put in dihedrals
    IPT=0
    IF(K.GT.1) IPT=NIC(4,K-1)
    N=NIC(4,K)-IPT
    DO I=1,N
       IPT=IPT+1
       CALL PATOM(II,ID1,ID2,MIP(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(JJ,ID1,ID2,MJP(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(KK,ID1,ID2,MKP(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(LL,ID1,ID2,MLP(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       IF(II.GT.0.AND.JJ.GT.0.AND.KK.GT.0.AND.LL.GT.0) THEN
          NPHI=NPHI+1
          IP(NPHI)=II
          JP(NPHI)=JJ
          KP(NPHI)=KK
          LP(NPHI)=LL
       ELSE
          IF(LWARN .AND. WRNLEV.GE.2 .and. prnlev >= 2) WRITE(OUTU,160) IRES,AA(K), &
               MIP(IPT),MJP(IPT),MKP(IPT),MLP(IPT)
       ENDIF
    ENDDO
    !
    ! put in improper dihedrals
    IPT=0
    IF(K.GT.1) IPT=NIC(5,K-1)
    N=NIC(5,K)-IPT
    DO I=1,N
       IPT=IPT+1
       CALL PATOM(II,ID1,ID2,MIM(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(JJ,ID1,ID2,MJM(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(KK,ID1,ID2,MKM(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(LL,ID1,ID2,MLM(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       IF(II.GT.0.AND.JJ.GT.0.AND.KK.GT.0.AND.LL.GT.0) THEN
          NIMPHI=NIMPHI+1
#if KEY_CFF==1
          IF (FFIELD.EQ.CFF) THEN
             IM(NIMPHI)=JJ
             JM(NIMPHI)=II
          ELSE
#endif 
             IM(NIMPHI)=II
             JM(NIMPHI)=JJ
#if KEY_CFF==1
          ENDIF
#endif 
          KM(NIMPHI)=KK
          LM(NIMPHI)=LL
       ELSE
          IF(LWARN .AND. WRNLEV.GE.2 .and. prnlev >= 2) WRITE(OUTU,170) IRES,AA(K), &
               MIM(IPT),MJM(IPT),MKM(IPT),MLM(IPT)
170       FORMAT(' ** WARNING ** IMPROPER NOT FOUND FOR RESIDUE ',I3,1X,A6, &
               '.'/,' ATOMS',4(1X,'"',A6,'"'),' WERE REQUESTED.')
       ENDIF
    ENDDO

#if KEY_CMAP==1
    ! put in cross-term maps
    IPT=0
    IF(K.GT.1) IPT=NIC(11,K-1)
    N=NIC(11,K)-IPT

    DO I=1,N
       IPT=IPT+1
       CALL PATOM(II1,ID1,ID2,MI1CT(IPT),ATYPE,IBASE,[IRES],0, &
            ISTART,ISTOP,LW)
       CALL PATOM(JJ1,ID1,ID2,MJ1CT(IPT),ATYPE,IBASE,[IRES],0, &
            ISTART,ISTOP,LW)
       CALL PATOM(KK1,ID1,ID2,MK1CT(IPT),ATYPE,IBASE,[IRES],0, &
            ISTART,ISTOP,LW)
       CALL PATOM(LL1,ID1,ID2,ML1CT(IPT),ATYPE,IBASE,[IRES],0, &
            ISTART,ISTOP,LW)
       CALL PATOM(II2,ID1,ID2,MI2CT(IPT),ATYPE,IBASE,[IRES],0, &
            ISTART,ISTOP,LW)
       CALL PATOM(JJ2,ID1,ID2,MJ2CT(IPT),ATYPE,IBASE,[IRES],0, &
            ISTART,ISTOP,LW)
       CALL PATOM(KK2,ID1,ID2,MK2CT(IPT),ATYPE,IBASE,[IRES],0, &
            ISTART,ISTOP,LW)
       CALL PATOM(LL2,ID1,ID2,ML2CT(IPT),ATYPE,IBASE,[IRES],0, &
            ISTART,ISTOP,LW)

       IF(II1.GT.0.AND.JJ1.GT.0.AND.KK1.GT.0.AND.LL1.GT.0 &
            .AND.II2.GT.0.AND.JJ2.GT.0.AND.KK2.GT.0.AND.LL2.GT.0) THEN
          NCRTERM=NCRTERM+1
          I1CT(NCRTERM)=II1
          J1CT(NCRTERM)=JJ1
          K1CT(NCRTERM)=KK1
          L1CT(NCRTERM)=LL1
          I2CT(NCRTERM)=II2
          J2CT(NCRTERM)=JJ2
          K2CT(NCRTERM)=KK2
          L2CT(NCRTERM)=LL2
       ELSE
          IF(LWARN .AND. WRNLEV.GE.2 .and. prnlev >= 2) WRITE(OUTU,171) IRES,AA(K), &
               MI1CT(IPT),MJ1CT(IPT),MK1CT(IPT),ML1CT(IPT), &
               MI2CT(IPT),MJ2CT(IPT),MK2CT(IPT),ML2CT(IPT)
       ENDIF
    ENDDO
#endif 
    !------------------------
    ! Put in lone pairs (this works only for relative and bisector type, so far...)
    ! colinear type is now also supported, 03.2013
    !     LONE PAIR
    !     MLP0CT           type of LP (geometric information)
    !     MLP1CT           first atom
    !     MLP2CT           second atom
    !     MLP3CT           third atom
    !     MLP4CT           fourth atom
    !     DELLPCT          lone-pairs true if term is to be deleted
    !     RLPCT            distance  (positive if RELATIVE, 0 if CENTER, negative if BISE)
    !     TLPCT            angle
    !     PLPCT            dihedral
    !
#if KEY_LONEPAIR==1
    IPT=0
    IF(K.GT.1) THEN
       IPT=NIC(12,K-1)
    END IF
    N=NIC(12,K)-IPT
    DO I=1,N
       IPT=IPT+1
       CALL PATOM(II,ID1,ID2,MLP1CT(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(JJ,ID1,ID2,MLP2CT(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(KK,ID1,ID2,MLP3CT(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(LL,ID1,ID2,MLP4CT(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       LOK=.FALSE.
       j=1
       do while (mlpct1d((ipt-1)*maxcent_hosts+j)/=blank)
            CALL PATOM(IIJ(J),ID1,ID2,MLPCT1D((IPT-1)*MAXCENT_HOSTS+J),ATYPE,IBASE,[IRES],0,ISTART,ISTOP,LW)  
            LOK=(IIJ(J).gt.0)
            j=j+1
       END DO

       IF(II.GT.0.AND.JJ.GT.0.AND.KK.GT.0.AND.LL.GT.0.OR.LOK) THEN
          IF (MLP0CT(IPT).eq.'CENT') THEN         ! MBA  
             NUMLP = NUMLP + 1
             if(.not.allocated(lpnhost)) call allocate_lonepair
             LPHPTR(NUMLP) = NUMLPH + 1
             NUMLPH = NUMLPH + J - 1  
             LPNHOST(NUMLP)   = J - 2 
             DO JK = 1, J-1 ! MBA
                LPHOST(LPHPTR(NUMLP)+(JK-1))   = IIJ(JK)
             END DO
             LPWGHT(NUMLP) = .FALSE.
             LPVALUE(1,NUMLP) = RLPCT(IPT)
             LPVALUE(2,NUMLP) = TLPCT(IPT)
             LPVALUE(3,NUMLP) = PLPCT(IPT)
             JK=IIJ(1)
             IMOVE(JK) = -1
          ELSEIF((MLP0CT(IPT).eq.'RELA').or.(MLP0CT(IPT).eq.'BISE')) THEN           
             NUMLP = NUMLP + 1
             if(.not.allocated(lpnhost)) call allocate_lonepair
             LPHPTR(NUMLP) = NUMLPH + 1
             NUMLPH = NUMLPH + 4
             LPNHOST(NUMLP)   = 3
             LPHOST(LPHPTR(NUMLP))   = II
             LPHOST(LPHPTR(NUMLP)+1) = JJ
             LPHOST(LPHPTR(NUMLP)+2) = KK
             LPHOST(LPHPTR(NUMLP)+3) = LL
             LPWGHT(NUMLP) = .FALSE.
             LPVALUE(1,NUMLP) = RLPCT(IPT)
             IF(MLP0CT(IPT).eq.'BISE') &
                  LPVALUE(1,NUMLP) = -RLPCT(IPT)
             LPVALUE(2,NUMLP) = TLPCT(IPT)
             LPVALUE(3,NUMLP) = PLPCT(IPT)
             IMOVE(II) = -1

          ELSE IF((MLP0CT(IPT).eq.'COLI')) THEN     ! MBA
             NUMLP = NUMLP + 1
             if(.not.allocated(lpnhost)) call allocate_lonepair
             LPHPTR(NUMLP) = NUMLPH + 1
             NUMLPH = NUMLPH + 3
             LPNHOST(NUMLP)   = 2
             LPHOST(LPHPTR(NUMLP))   = II
             LPHOST(LPHPTR(NUMLP)+1) = JJ
             LPHOST(LPHPTR(NUMLP)+2) = KK
             LPWGHT(NUMLP) = .FALSE.
             LPVALUE(1,NUMLP) = RLPCT(IPT)
             LPVALUE(2,NUMLP) = TLPCT(IPT)
             IMOVE(II) = -1
          ELSE 
             CALL WRNDIE(-4,'<GENPSF>','Unknown lone pair type')
          ENDIF
       ENDIF
    ENDDO
#if KEY_PARALLEL==1
    if (numlp > 0) CALL PSND8(LPVALUE,3*NUMLP)         
#endif
#endif 
    !
    !------------------------
    !
    ! Put in anisotropic polarizability terms
    !     MIANIS           first atom
    !     MJANIS           second atom
    !     MKANIS           third atom
    !     MLANIS           fourth atom
    !     DELANIS          true if anistropic term is to be deleted
    !     A11ANIS          principal axis of the tensor
    !     A22ANIS
    !
    IPT=0
    IF(K.GT.1) IPT=NIC(13,K-1)
    N=NIC(13,K)-IPT
    DO I=1,N
       IPT=IPT+1
       CALL PATOM(II,ID1,ID2,MIANIS(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(JJ,ID1,ID2,MJANIS(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(KK,ID1,ID2,MKANIS(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(LL,ID1,ID2,MLANIS(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       IF(II.GT.0.AND.JJ.GT.0.AND.KK.GT.0.AND.LL.GT.0) THEN

          !         write(*,'(1x,5(A,2x),2x,3(A,f8.3))')
          !    &    'Anisotropy: ',
          !    &    MIANIS(IPT),MJANIS(IPT),MKANIS(IPT),MLANIS(IPT),
          !    &    ' K11 ', A11ANIS(IPT),
          !    &    ' K22 ',A22ANIS(IPT)
          !         write(*,*) 'atom number ',II,JJ,KK,LL

          !     NANISO  - Maximum number of anisotropic terms in the list
          !     LSTANI1 - 1st component in the list
          !     LSTANI2 - 2nd component in the list
          !     LSTANI3 - 3rd component in the list
          !     LSTANI4 - 4th component in the list
          !     Fraction of the total isotropic tensor along principal axis
          !     A11    - principal axis 11
          !     A22    - principal axis 22
          !     A33    - principal axis 33
          !     K11    - force constant for parallel component
          !     K22    - force constant for perpendicular component
          !     K33    - force constant orthogonal to 11 & 22

          IF( (A11ANIS(IPT).ne.ZERO).and. &
               (A22ANIS(IPT).ne.ZERO) ) THEN
             NANISO = NANISO + 1
             LSTANI1(NANISO) = II
             LSTANI2(NANISO) = JJ
             LSTANI3(NANISO) = KK
             LSTANI4(NANISO) = LL
             A33 = THREE-A11ANIS(IPT)-A22ANIS(IPT)
             IF(A33.le.zero) CALL WRNDIE(5,'<RTFRDR>', &
                  'ANISO: A33 negative ')

             K11(NANISO) = ONE/A11ANIS(IPT)
             K22(NANISO) = ONE/A22ANIS(IPT)
             K33(NANISO) = ONE/A33

             if(PRNLEV.GE.8)then
                write(*,'(1x,a,5i6,3f10.3)')  &
                     'Add anisotropy for one atom: ', &
                     naniso,ii,jj,kk,ll, &
                     K11(NANISO),K22(NANISO),K33(NANISO)
             endif
          ELSE
             CALL WRNDIE(-4,'<GENPSF>','Zero aniso force constant')
          ENDIF

       ENDIF
    ENDDO
    !
    !------------------------
    !
    ! put in donors
    IPT=0
    IF(K.GT.1) IPT=NIC(7,K-1)
    N=NIC(7,K)-IPT
    DO I=1,N
       IPT=IPT+1
       CALL PATOM(II,ID1,ID2,MH(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(JJ,ID1,ID2,MD(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       IF(II.EQ.MARK) II=0
       IF(JJ.GT.0 .AND. (II.GT.0 .OR. MH(IPT).EQ.BLANK)) THEN
          NDON=NDON+1
          IHD1(NDON)=II
          IDON(NDON)=JJ
       ELSE
          IF(LWARN .AND. WRNLEV.GE.2 .and. prnlev >= 2) &
               WRITE(OUTU,190) IRES,AA(K),MH(IPT),MD(IPT)
       ENDIF
    ENDDO
    !
    ! put in acceptors
    IPT=0
    IF(K.GT.1) IPT=NIC(8,K-1)
    N=NIC(8,K)-IPT
    DO I=1,N
       IPT=IPT+1
       CALL PATOM(II,ID1,ID2,MA(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       CALL PATOM(JJ,ID1,ID2,MAA(IPT),ATYPE,IBASE,[IRES],0,ISTART, &
            ISTOP,LW)
       IF(JJ.EQ.MARK) JJ=0
       IF(II.GT.0) THEN
          NACC=NACC+1
          IACC(NACC)=II
          IAC1(NACC)=JJ
       ELSE
          IF(LWARN .AND. WRNLEV.GE.2 .and. prnlev >= 2) &
               WRITE(OUTU,200) IRES,AA(K),MA(IPT),MAA(IPT)
       ENDIF
    ENDDO
    !



    ! put in build/ic elements
    IF(LSETIC) THEN
       IPT=0
       IF(K.GT.1) IPT=NIC(9,K-1)
       N=NIC(9,K)-IPT
       DO I=1,N
          IPT=IPT+1
          CALL PATOM(II,ID1,ID2,BARI(IPT),ATYPE,IBASE,[IRES],0, &
               ISTART,ISTOP,LW)
          CALL PATOM(JJ,ID1,ID2,BARJ(IPT),ATYPE,IBASE,[IRES],0, &
               ISTART,ISTOP,LW)
          CALL PATOM(KK,ID1,ID2,BARK(IPT),ATYPE,IBASE,[IRES],0, &
               ISTART,ISTOP,LW)
          CALL PATOM(LL,ID1,ID2,BARL(IPT),ATYPE,IBASE,[IRES],0, &
               ISTART,ISTOP,LW)
          IF(II.GT.0.OR.JJ.GT.0.OR.KK.GT.0.OR.LL.GT.0) THEN
             IF(icr_struct%intlen.LE.icr_struct%lenic) THEN  ! get more space...
                IRSZ=icr_struct%lenic+200
                CALL reintc_new(irsz,icr_struct)
             ENDIF
             icr_struct%lenic = icr_struct%lenic+1
             CALL PUTICEL(II,JJ,KK,LL,BART(IPT),ICB1(IPT),ICB2(IPT), &
                  ICTH1(IPT),ICTH2(IPT),ICPHI(IPT),&
                  icr_struct%lenic, &
                  icr_struct%B1ic,icr_struct%B2ic, &
                  icr_struct%T1ic,icr_struct%T2ic, &
                  icr_struct%PIC, icr_struct%IAR, &
                  icr_struct%JAR, icr_struct%KAR, &
                  icr_struct%LAR, icr_struct%TAR)
          ELSE
             IF(LWARN .AND. WRNLEV.GE.2 .and. prnlev >= 2) WRITE(OUTU,210) IRES,AA(K), &
                  BARI(IPT),BARJ(IPT),BARK(IPT),BARL(IPT)
          ENDIF
       ENDDO
    ENDIF
140 FORMAT(' ** WARNING ** BOND NOT FOUND FOR RESIDUE ',I3,1X,A6, &
         '.'/,' ATOMS',2(1X,'"',A6,'"'),' WERE REQUESTED.')
150 FORMAT(' ** WARNING ** ANGLE NOT FOUND FOR RESIDUE ',I3,1X,A6, &
         '.'/,' ATOMS',3(1X,'"',A6,'"'),' WERE REQUESTED.')
160 FORMAT(' ** WARNING ** DIHEDRAL NOT FOUND FOR RESIDUE ',I3,1X,A6, &
         '.'/,' ATOMS',4(1X,'"',A6,'"'),' WERE REQUESTED.')
#if KEY_CMAP==1
171 FORMAT(' ** WARNING ** CROSSTERM NOT FOUND FOR RESIDUE ',I3,1X,A6, &
         '.'/,' ATOMS',4(1X,'"',A6,'"'), &
         /,' ATOMS',4(1X,'"',A6,'"'),' WERE REQUESTED.')
#endif 
190 FORMAT(' ** WARNING ** DONOR NOT FOUND FOR RESIDUE ',I3,1X,A6, &
         '.'/,' ATOMS',2(1X,'"',A6,'"'),' WERE REQUESTED.')
200 FORMAT(' ** WARNING ** ACCEPTOR NOT FOUND FOR RESIDUE ',I3,1X,A6, &
         '.'/,' ATOMS',2(1X,'"',A6,'"'),' WERE REQUESTED.')
210 FORMAT(' ** WARNING ** IC NOT FOUND FOR RESIDUE ',I3,1X,A6, &
         '.'/,' ATOMS',4(1X,'"',A6,'"'),' WERE REQUESTED.')


    return
  end subroutine sub400
  !
  !=======================================================================
  ! to find-this-residue

  subroutine sub500
    FOUND=.FALSE.
    DO J=1,NRTRS
       IF(TEMP.EQ.AA(J)) THEN
          IF(RTRTYP(J).NE.IPATCH) THEN
             IF(WRNLEV.GE.2) WRITE(OUTU,231) TEMP
             CALL DIE
          ELSE
             ITEMP=J
             FOUND=.TRUE.
          ENDIF
       ENDIF
    ENDDO
    IF(.NOT.FOUND) THEN
       IF(WRNLEV.GE.2) WRITE(OUTU,233) TEMP
       CALL DIE
    ENDIF
231 FORMAT(/' ***** ERROR in GENIC ***** Residue ''',A6, &
         ''' has the wrong patch type')
233 FORMAT(/' ***** ERROR in GENIC ***** Residue ''',A6,''' was', &
         ' not found.')

    return
  end subroutine sub500

  !
  !=======================================================================
  ! to add-this-atom

  subroutine sub600
    LADDED=.TRUE.
    NATP=NATP+1
    NL=NATP+NATOM
    IBLO(NL)=I+NATPT
    ITEMP=MAC(I+NATPT)
    IAC(NL)=ITEMP
    AMASS(NL)=ARMASS(ITEMP)
    RSCLF(NL)=ONE
#if KEY_WCA==1
    WCA(NL)=ONE               
#endif
#if KEY_MMFF==1
    if(ffield == mmff) AtNum(NL)=AtNumR(I+NATPT) 
#endif
    ATYPE(NL)=FTP(I+NATPT)
    CG(NL)=CHG(I+NATPT)
    IF(LDRUDE)THEN
       !        write(*,*) '600 continue ',NL,I+NATPT,
       !    &   CG(NL),AMASS(NL),ALPH(I+NATPT),THOL(I+NATPT)
       ALPHADP(NL)=ALPH(I+NATPT)
       THOLEI(NL)=THOL(I+NATPT)
    ENDIF
#if KEY_CHEQ==1
    ECH(NL)=ECHH(I+NATPT)
    EHA(NL)=EHAH(I+NATPT)
#endif 
    IMOVE(NL)=0
    IF(IGRLST.NE.GRPR(I+NATPT)) THEN
       IGRLST=GRPR(I+NATPT)
       NGRP=NGRP+1
       IGPBS(NGRP)=NL-1
       IF(AA(K).EQ.'ST2 ') THEN
#if KEY_NOST2==0
          IGPTYP(NGRP)=3
#else /**/
          CALL WRNDIE(-1,'<GENPSF>','ST2 code is not compiled.')
          IGPTYP(NGRP)=0
#endif 
       ELSE
          IGPTYP(NGRP)=0
       ENDIF
       IMOVEG(NGRP)=0
    ENDIF

    return
  end subroutine sub600

  !
  !=======================================================================
  ! to check-this-atom

  subroutine sub700
    NATPT2=0
    IF(IPAT.GT.1) NATPT2=NIC(1,IPAT-1)
    NAT2=NIC(1,IPAT)-NATPT2
    DO I2=1,NAT2
       IF(FTP(I2+NATPT2).EQ.FTP(I+NATPT)) THEN
          OK=.FALSE.
          IMATCH=I2+NATPT2
       ENDIF
    ENDDO

    return
  end subroutine sub700

  !=======================================================================
  !
END SUBROUTINE GENIC

SUBROUTINE AUTGEN(NBONDL,NATBON,IATBON,LTHETA,LPHI,LDRUDE)
  !
  !     This routine automatically generates the angle and dihedral
  !     lists for a set of bonds. This may be done in lieu of specifying
  !     angles and dihedrals explicitly in the topology file. This
  !     procedure only adds terms, there is no check for duplicate terms
  !     from the existing angle or dihedral list. The multiple dihedral
  !     potential may be obtained by specifying a dihedral once in the
  !     topology file and the second will be generated here.
  !
  !     For angles, IT(i)<KT(i) For Dihedrals, IP(i)<LP(i)
  !
  !RCZ - 90/03/22 Ryszard Czerminski
  !     subroutine modified to avoid generation of dihedrals
  !     between linear bonds
  !RCZ
  !
  !     By Bernard R. Brooks  15-FEB-1984
  !
  use exfunc, only: order5
  use stream
  use dimens_fcm
  use psf
  use code
  use param

  INTEGER NBONDL
  INTEGER NATBON(*),IATBON(IATBMX,*)
  LOGICAL LTHETA,LPHI,LDRUDE, OK
  !
  !
  INTEGER I,IBT,JBT,J,IJ,IJA,K,IK,IKA
  LOGICAL SKIP
  EXTERNAL  EXCH5
  !
  ! Construct the cross reference bond list
  !
  DO I=1,NATOM
     NATBON(I)=0
  ENDDO
  DO I=NBONDL,NBOND
     IBT=IB(I)
     JBT=JB(I)
     IF(IBT.GT.0 .AND. JBT.GT.0) THEN
        NATBON(IBT)=NATBON(IBT)+1
        IF(NATBON(IBT).GT.IATBMX+2) THEN
           IF(WRNLEV.GE.2)then
              WRITE(OUTU,335) NATBON(IBT),IBT
335           FORMAT(' <AUTGEN>: ',I5, &
                   ' is too  many bonds for atom',I5, &
                   ' Check code')
              CALL DIEWRN(-4)
           ENDIF
        ENDIF
        IATBON(NATBON(IBT),IBT)=I
        NATBON(JBT)=NATBON(JBT)+1
        IF(NATBON(JBT).GT.IATBMX) THEN
           IF(WRNLEV.GE.2) WRITE(OUTU,335) JBT
           CALL DIEWRN(-4)
        ENDIF
        IATBON(NATBON(JBT),JBT)=-I
     ENDIF
  ENDDO
  !
  ! Now make the unsorted list of 1-3 interactions by taking bonds
  ! and extending in a direction one bond.
  !
  IF(LTHETA) THEN
     DO I=1,NBOND
        IBT=IB(I)
        JBT=JB(I)
        DO J=1,NATBON(IBT)
           IJ=IATBON(J,IBT)
           IF(IJ.GT.0) THEN
              IJA=JB(IJ)
           ELSE
              IJA=IB(ABS(IJ))
           ENDIF
           IF(IJA.GT.JBT) THEN
              OK = (ATYPE(JBT)(1:1) == 'D') .OR. &
                   (ATYPE(IBT)(1:1) == 'D') .OR. &
                   (ATYPE(IJA)(1:1) == 'D') .OR. &
                   (ATYPE(JBT)(1:2) == 'LP') .OR. &
                   (ATYPE(IBT)(1:2) == 'LP') .OR. &
                   (ATYPE(IJA)(1:2) == 'LP')
              IF(OK)THEN
                 IF(PRNLEV.GE.8) &
                      write(outu,*) 'skip drude angle',jbt,ibt,ija
              ELSE
                 NTHETA=NTHETA+1
                 IT(NTHETA)=JBT
                 JT(NTHETA)=IBT
                 KT(NTHETA)=IJA
              ENDIF
           ENDIF
        ENDDO
        DO J=1,NATBON(JBT)
           IJ=IATBON(J,JBT)
           IF(IJ.GT.0) THEN
              IJA=JB(IJ)
           ELSE
              IJA=IB(ABS(IJ))
           ENDIF
           IF(IJA.GT.IBT) THEN
              OK = (ATYPE(IBT)(1:1) == 'D') .OR. &
                   (ATYPE(JBT)(1:1) == 'D') .OR. &
                   (ATYPE(IJA)(1:1) == 'D') .OR. &
                   (ATYPE(IBT)(1:2) == 'LP') .OR. &
                   (ATYPE(JBT)(1:2) == 'LP') .OR. &
                   (ATYPE(IJA)(1:2) == 'LP')
              IF(OK)THEN
                 IF(PRNLEV.GE.8) &
                      write(outu,*) 'skip drude angle',ibt,jbt,ija
              ELSE
                 NTHETA=NTHETA+1
                 IT(NTHETA)=IBT
                 JT(NTHETA)=JBT
                 KT(NTHETA)=IJA
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  !
  !
  ! Now make the unsorted list of 1-4 interactions by taking bonds
  ! and extending in each direction one bond.
  !
  IF(.NOT.LPHI) RETURN
  !RCZ
  ! to avoid autogeneration of dihedral terms in which three
  ! atoms are colinear we have to call here CODES subroutine to
  ! determine for which triples equilibrium valence angle is 180
  ! degrees.
  CALL CODES(0,ICT,0,0, &
       NATOM,IMOVE,IAC,0,0,0, &
       NTHETA,IT,JT,KT,0,0,0,0,0,0,0,0,0,0, &
       .FALSE.,0,                                   & ! DRUDE
#if KEY_CMAP==1
       0,0,0,0,0,0,0,0,0,0,                         & 
#endif
       .FALSE.,.FALSE.)
  !RCZ
  DO I=1,NBOND
     IBT=IB(I)
     JBT=JB(I)
     DO J=1,NATBON(IBT)
        IJ=IATBON(J,IBT)
        IF(ABS(IJ).NE.I) THEN
           IF(IJ.GT.0) THEN
              IJA=JB(IJ)
           ELSE
              IJA=IB(ABS(IJ))
           ENDIF
           DO K=1,NATBON(JBT)
              IK=IATBON(K,JBT)
              IF(ABS(IK).NE.I) THEN
                 IF(IK.GT.0) THEN
                    IKA=JB(IK)
                 ELSE
                    IKA=IB(ABS(IK))
                 ENDIF
                 !

                 OK = (ATYPE(IJA)(1:1) == 'D') .OR. &
                      (ATYPE(IBT)(1:1) == 'D') .OR. &
                      (ATYPE(JBT)(1:1) == 'D') .OR. &
                      (ATYPE(IKA)(1:1) == 'D') .OR. &
                      (ATYPE(IJA)(1:2) == 'LP') .OR. &
                      (ATYPE(IBT)(1:2) == 'LP') .OR. &
                      (ATYPE(JBT)(1:2) == 'LP') .OR. &
                      (ATYPE(IKA)(1:2) == 'LP')
                 IF(OK)THEN
                    IF(PRNLEV.GE.8) &
                         write(outu,*) 'skip drude dihedral',ija,ibt,jbt,ika
                 ELSE

                    IF(IKA.EQ.IBT .OR. IJA.EQ.JBT) THEN
                       CALL WRNDIE(-2,'<AUTGEN>', &
                            'DOUBLE BOND SPECS FOR SOME ATOM PAIRS')
                    ELSE IF(IKA.EQ.IJA .AND.  &
                         AMASS(IKA).GT.0.002.AND.AMASS(IBT).GT.0.002.AND.  &
                         AMASS(JBT).GT.0.002) THEN
                       !va Skip this part if Lone-Pairs are involved in 3-member ring
                       CALL WRNDIE(0,'<AUTGEN>', &
                            'THREE MEMBER RING FOUND, NO DIHEDRAL ADDED')
                    ELSE
                       !RCZ/BRB
                       IF(IKA.LT.IJA) THEN
                          CALL CHECKDH(NTHETA,IT,JT,KT,IKA,JBT,IBT,IJA, &
                               SKIP,AMASS)
                          IF(.NOT.SKIP) THEN
                             NPHI=NPHI+1
                             IP(NPHI)=IKA
                             JP(NPHI)=JBT
                             KP(NPHI)=IBT
                             LP(NPHI)=IJA
                          ENDIF
                       ELSE
                          CALL CHECKDH(NTHETA,IT,JT,KT,IJA,IBT,JBT,IKA, &
                               SKIP,AMASS)
                          IF(.NOT.SKIP) THEN
                             NPHI=NPHI+1
                             IP(NPHI)=IJA
                             JP(NPHI)=IBT
                             KP(NPHI)=JBT
                             LP(NPHI)=IKA
                          ENDIF
                       ENDIF
                       !RCZ/BRB
                    ENDIF
                    !
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
  !  Now sort the new dihedral list to allow for double potential terms.
  !
  CALL SORT(NPHI,EXCH5,ORDER5,IP,JP,KP,LP,0,0,0,4)
  !
  RETURN
END SUBROUTINE AUTGEN

SUBROUTINE CHECKDH(NTHETA,IT,JT,KT,ID,JD,KD,LD,SKIP,AMASS)
  !
  !RCZ 90/03/22 checking if there are equilibrium valence angles
  !             equal 180 degrees between ID,JD,KD,LD.
  !BRB 91/10/16 Modified for efficiency.
  !
  use consta
  use dimens_fcm
  use code
  use param
  use stream

  INTEGER NTHETA,IT(*),JT(*),KT(*)
  INTEGER ID,JD,KD,LD
  LOGICAL SKIP
  real(chm_real) AMASS(*)
  INTEGER N,ITH,II,JJ,KK
  real(chm_real), parameter :: EPS=0.01
  !
  SKIP=.TRUE.
  N=0
  !
  ! Test first three atoms
  II=MIN(ID,KD)
  JJ=JD
  KK=MAX(ID,KD)
  DO ITH=1,NTHETA
     IF(II.EQ.IT(ITH)) THEN
        IF(JJ.EQ.JT(ITH)) THEN
           IF(KK.EQ.KT(ITH)) THEN
              N=N+1
              IF(ABS(CTB(ICT(ITH))-PI).LT.EPS .OR. &
                   ABS(CTB(ICT(ITH))).LT.EPS) THEN
                 !va Skip warning message in case of Lone-Pairs unless high PRNLEV is defined
                 IF(PRNLEV.GE.7 .OR. (PRNLEV.GE.2 .AND. &
                      AMASS(ID).GT.0.002.AND.AMASS(JD).GT.0.002.AND. &
                      AMASS(KD).GT.0.002.AND.AMASS(LD).GT.0.002)) THEN
                    WRITE(OUTU,80) ID,JD,KD,LD
80                  FORMAT(' CHECKDH> dihedral :',4I5, &
                         ' will NOT be generated')
                 ENDIF
                 RETURN
              ENDIF
              GOTO 105
           ENDIF
        ENDIF
     ENDIF
  ENDDO
105 CONTINUE
  !
  ! Test last three atoms
  II=MIN(JD,LD)
  JJ=KD
  KK=MAX(JD,LD)
  DO ITH=1,NTHETA
     IF(II.EQ.IT(ITH)) THEN
        IF(JJ.EQ.JT(ITH)) THEN
           IF(KK.EQ.KT(ITH)) THEN
              N=N+1
              IF(ABS(CTB(ICT(ITH))-PI).LT.EPS .OR. &
                   ABS(CTB(ICT(ITH))).LT.EPS) THEN
                 IF(PRNLEV.GE.7 .OR. (PRNLEV.GE.2 .AND. &
                      AMASS(ID).GT.0.002.AND.AMASS(JD).GT.0.002.AND. &
                      AMASS(KD).GT.0.002.AND.AMASS(LD).GT.0.002)) THEN
                    WRITE(OUTU,80) ID,JD,KD,LD
                 ENDIF
                 RETURN
              ENDIF
              GOTO 205
           ENDIF
        ENDIF
     ENDIF
  ENDDO
205 CONTINUE
  !
  !     IF(N.GT.1) WRITE(OUTU,200) N,ID,JD,KD,LD
  ! 200 FORMAT(' CHECKDH> WARNING; N=',I5,' ID,JD,KD,LD=',4I5)
  !.....CHECKDH
  SKIP=.FALSE.
  RETURN
END SUBROUTINE CHECKDH


SUBROUTINE DUPSEG(LSETIC,SEGD,CMPLTD)
  !
  !     This routine duplicates all entries of a given segment into
  !     a new segment. This is especially usefull in setting up a crystal
  !     for viewing or other analysis.
  !     Refernces from PSF, COORD, COORDC and INTCR are copied.
  !
  !     Bernard R. Brooks    10/22/83
  !

#if KEY_CHEQ==1
  use cheq, only: ech, eha  
#endif
  use intcor_module
  use intcor2,only:geticel,puticel
  use dimens_fcm
  use number
  use bases_fcm
  use psf
  use coord
  use coordc
  use mmffm
  use ffieldm

  LOGICAL LSETIC,CMPLTD
  CHARACTER(len=*) SEGD
  !
  INTEGER NSEGM,ISEG,I,IRESI,IRESF,IATMI,IATMF,IATOFF,IRES,NGRPO
  INTEGER IGRP,NNBI,IAT,NNBF,NBONDO,NTHETO
  INTEGER NPHIO,NIMPHO,NDONO,NACCO,LENICO
#if KEY_CMAP==1
  INTEGER NCRT0
#endif 
  CHARACTER(len=8) CHX              ! yw A4->A8
  INTEGER II,JJ,KK,LL,IRSZ
  LOGICAL TT
  real(chm_real)  BB1,BB2,TT1,TT2,PP
  !
  ! Do the psf arrays and numbers
  !
  NSEGM=NSEG-1
  ISEG=-1
  DO I=1,NSEGM
     IF(SEGID(I).EQ.SEGD) ISEG=I
  ENDDO
  IF(ISEG.LT.0) THEN
     CMPLTD=.FALSE.
     CALL WRNDIE(-1,'<DUPSEG>','COPY SEGMENT NOT PRESENT')
     RETURN
  ENDIF
  !
  IRESI=NICTOT(ISEG)+1
  IRESF=NICTOT(ISEG+1)
  IATMI=IBASE(IRESI)+1
  IATMF=IBASE(IRESF+1)
  IATOFF=NATOM+1-IATMI
  !
  DO IRES=IRESI,IRESF
     NRES=NRES+1
     CHX=RES(IRES)
     RES(NRES)=CHX
     CHX=RESID(IRES)
     RESID(NRES)=CHX
     IBASE(NRES)=IATOFF+IBASE(IRES)
  ENDDO
  NICTOT(NSEG+1)=NRES
  !
  NGRPO=NGRP
  DO IGRP=1,NGRPO
     IF((IGPBS(IGRP)+1.LE.IATMF).AND.(IGPBS(IGRP+1).GE.IATMI)) THEN
        NGRP=NGRP+1
        IGPTYP(NGRP)=IGPTYP(IGRP)
        IMOVEG(NGRP)=IMOVEG(IGRP)
        IGPBS(NGRP)=IGPBS(IGRP)+IATOFF
#if KEY_NOST2==0
        IF(IGPTYP(IGRP).EQ.3) NST2=NST2+1
#endif 
     ENDIF
  ENDDO
  IGPBS(NGRP+1)=IATOFF+IATMF
  !
  IF(IATMI.EQ.1) THEN
     NNBI=0
  ELSE
     NNBI=IBLO(IATMI-1)
  ENDIF
  DO IAT=IATMI,IATMF
     NATOM=NATOM+1
     CHX=ATYPE(IAT)
     ATYPE(NATOM)=CHX
     CG(NATOM)=CG(IAT)
     ISDRUDE(NATOM)=ISDRUDE(IAT)
     AMASS(NATOM)=AMASS(IAT)
#if KEY_CHEQ==1
     ECH(NATOM)=ECH(IAT)
     EHA(NATOM)=EHA(IAT)
#endif 
     RSCLF(NATOM)=ONE
#if KEY_WCA==1
     WCA(NATOM)=ONE           
#endif
#if KEY_MMFF==1
     if(ffield == mmff) AtNum(NATOM)=AtNum(IAT)  
#endif
     IAC(NATOM)=IAC(IAT)
     IMOVE(NATOM)=IMOVE(IAT)
     IBLO(NATOM)=IBLO(IAT)-NNBI+NNB
  ENDDO
  IBASE(NRES+1)=NATOM
  !
  NNBF=IBLO(IATMF)
  NNBI=NNBI+1
  DO I=NNBI,NNBF
     NNB=NNB+1
     INB(NNB)=INB(I)+IATOFF
  ENDDO
  !
  NBONDO=NBOND
  DO I=1,NBONDO
     IF(IB(I).LE.IATMF.AND.IB(I).GE.IATMI) THEN
        IF(JB(I).LE.IATMF.AND.JB(I).GE.IATMI) THEN
           NBOND=NBOND+1
           IB(NBOND)=IB(I)+IATOFF
           JB(NBOND)=JB(I)+IATOFF
#if KEY_MMFF==1
           if(ffield == mmff) BondType(NBOND)=BondType(I) 
#endif
        ENDIF
     ENDIF
  ENDDO
  !
  NTHETO=NTHETA
  DO I=1,NTHETO
     IF(IT(I).LE.IATMF.AND.IT(I).GE.IATMI) THEN
        IF(JT(I).LE.IATMF.AND.JT(I).GE.IATMI) THEN
           IF(KT(I).LE.IATMF.AND.KT(I).GE.IATMI) THEN
              NTHETA=NTHETA+1
              IT(NTHETA)=IT(I)+IATOFF
              JT(NTHETA)=JT(I)+IATOFF
              KT(NTHETA)=KT(I)+IATOFF
           ENDIF
        ENDIF
     ENDIF
  ENDDO
  !
  NPHIO=NPHI
  DO I=1,NPHIO
     IF(IP(I).LE.IATMF.AND.IP(I).GE.IATMI) THEN
        IF(JP(I).LE.IATMF.AND.JP(I).GE.IATMI) THEN
           IF(KP(I).LE.IATMF.AND.KP(I).GE.IATMI) THEN
              IF(LP(I).LE.IATMF.AND.LP(I).GE.IATMI) THEN
                 NPHI=NPHI+1
                 IP(NPHI)=IP(I)+IATOFF
                 JP(NPHI)=JP(I)+IATOFF
                 KP(NPHI)=KP(I)+IATOFF
                 LP(NPHI)=LP(I)+IATOFF
              ENDIF
           ENDIF
        ENDIF
     ENDIF
  ENDDO
  !
  NIMPHO=NIMPHI
  DO I=1,NIMPHO
     IF(IM(I).LE.IATMF.AND.IM(I).GE.IATMI) THEN
        IF(JM(I).LE.IATMF.AND.JM(I).GE.IATMI) THEN
           IF(KM(I).LE.IATMF.AND.KM(I).GE.IATMI) THEN
              IF(LM(I).LE.IATMF.AND.LM(I).GE.IATMI) THEN
                 NIMPHI=NIMPHI+1
                 IM(NIMPHI)=IM(I)+IATOFF
                 JM(NIMPHI)=JM(I)+IATOFF
                 KM(NIMPHI)=KM(I)+IATOFF
                 LM(NIMPHI)=LM(I)+IATOFF
              ENDIF
           ENDIF
        ENDIF
     ENDIF
  ENDDO

#if KEY_CMAP==1
  NCRT0=NCRTERM
  DO I=1,NCRT0
     IF(I1CT(I).LE.IATMF.AND.I1CT(I).GE.IATMI) THEN
        IF(J1CT(I).LE.IATMF.AND.J1CT(I).GE.IATMI) THEN
           IF(K1CT(I).LE.IATMF.AND.K1CT(I).GE.IATMI) THEN
              IF(L1CT(I).LE.IATMF.AND.L1CT(I).GE.IATMI) THEN
                 IF(I2CT(I).LE.IATMF.AND.I2CT(I).GE.IATMI) THEN
                    IF(J2CT(I).LE.IATMF.AND.J2CT(I).GE.IATMI) THEN
                       IF(K2CT(I).LE.IATMF.AND.K2CT(I).GE.IATMI) THEN
                          IF(L2CT(I).LE.IATMF.AND.L2CT(I).GE.IATMI) THEN
                             NCRTERM=NCRTERM+1
                             I1CT(NCRTERM)=I1CT(I)+IATOFF
                             J1CT(NCRTERM)=J1CT(I)+IATOFF
                             K1CT(NCRTERM)=K1CT(I)+IATOFF
                             L1CT(NCRTERM)=L1CT(I)+IATOFF
                             I2CT(NCRTERM)=I2CT(I)+IATOFF
                             J2CT(NCRTERM)=J2CT(I)+IATOFF
                             K2CT(NCRTERM)=K2CT(I)+IATOFF
                             L2CT(NCRTERM)=L2CT(I)+IATOFF
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
     ENDIF
  ENDDO
#endif 
  !
  NDONO=NDON
  DO I=1,NDONO
     IF(IDON(I).LE.IATMF.AND.IDON(I).GE.IATMI) THEN
        NDON=NDON+1
        IDON(NDON)=IDON(I)+IATOFF
        IF(IHD1(I).GT.0) IHD1(NDON)=IHD1(I)+IATOFF
     ENDIF
  ENDDO
  !
  NACCO=NACC
  DO I=1,NACCO
     IF(IACC(I).LE.IATMF.AND.IACC(I).GE.IATMI) THEN
        NACC=NACC+1
        IACC(NACC)=IACC(I)+IATOFF
        IF(IAC1(I).GT.0) IAC1(NACC)=IAC1(I)+IATOFF
     ENDIF
  ENDDO
  !
  ! Reset the coordinates
  !
  CALL ATMINI(IATMI+IATOFF,IATMF+IATOFF)
  !
  ! Do internal coordinates if requested
  !
  IF(LSETIC) THEN
     LENICO=icr_struct%lenic
     DO I=1,LENICO
        CALL GETICEL(II,JJ,KK,LL,TT,BB1,BB2,TT1,TT2,PP,I, &
             icr_struct%B1ic,icr_struct%B2ic, &
             icr_struct%T1ic,icr_struct%T2ic, &
             icr_struct%PIC, icr_struct%IAR, &
             icr_struct%JAR, icr_struct%KAR, &
             icr_struct%LAR, icr_struct%TAR)
        !
        IF(II.LE.IATMF.AND.II.GE.IATMI) THEN
           IF(JJ.LE.IATMF.AND.JJ.GE.IATMI) THEN
              IF(KK.LE.IATMF.AND.KK.GE.IATMI) THEN
                 IF(LL.LE.IATMF.AND.LL.GE.IATMI) THEN
                    !
                    IF(icr_struct%intlen.LE.icr_struct%lenic) THEN  ! get more space...
                       IRSZ=icr_struct%lenic+200
                       CALL REINTC_NEW(IRSZ,icr_struct)
                    ENDIF
                    icr_struct%lenic=icr_struct%lenic+1
                    !
                    II=II+IATOFF
                    JJ=JJ+IATOFF
                    KK=KK+IATOFF
                    LL=LL+IATOFF
                    CALL PUTICEL(II,JJ,KK,LL,TT,BB1,BB2,TT1,TT2,PP,&
                         icr_struct%lenic, &
                         icr_struct%B1ic,icr_struct%B2ic, &
                         icr_struct%T1ic,icr_struct%T2ic, &
                         icr_struct%PIC, icr_struct%IAR, &
                         icr_struct%JAR, icr_struct%KAR, &
                         icr_struct%LAR, icr_struct%TAR)
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDIF
  !
  CMPLTD=.TRUE.
  RETURN
END SUBROUTINE DUPSEG

SUBROUTINE PATOM(PATATM,PATRES,PTYP,ATOM,ATYPE,IBASE,IPR, &
     NPR,ISTART,ISTOP,LWARN)
  !
  !     This routine returns the atom number given a residue number
  !     preceding the atom name. the residue number is found by a lookup
  !     into IPR. If no number precedes the atom iupac name or if NPR
  !     equal zero the first residue in IPR is assumed.
  !
  !     23-JAN-83, Axel Brunger
  !
  use stream

  INTEGER   PATATM, PATRES
  CHARACTER(len=*)   PTYP, ATOM, ATYPE(*)
  INTEGER   IBASE(*)
  INTEGER   IPR(*), NPR, ISTART, ISTOP
  LOGICAL   LWARN
  !
  INTEGER   I, J
  !
  !     character definitions following
  !
  INTEGER, parameter :: MARK = -99999999
  CHARACTER(len=1) CHAR
  CHARACTER(len=10) WORD,WORD2                  ! yw A4->A8
  CHARACTER(len=1) IN(9)
  DATA IN/'1','2','3','4','5','6','7','8','9'/
  !
  PATRES=0
  CHAR=ATOM(1:1)
  WORD=ATOM
  !
  IF(NPR.EQ.0) THEN
     J=1
  ELSE
     J=1
     DO WHILE(CHAR.NE.IN(J).AND.J.LT.9)
        J=J+1
     ENDDO
     IF(CHAR.EQ.IN(J)) THEN
        IF(J.GT.NPR) GOTO 900
        WORD2=WORD(2:)
        WORD=WORD2
        CHAR=WORD(1:1)
     ELSE
        J=1
     ENDIF
  ENDIF
  !
  PATRES=IPR(J)
100 CONTINUE
  IF(CHAR.EQ.'-') THEN
     PATRES=PATRES-1
     WORD2=WORD(2:)
     WORD=WORD2
     CHAR=WORD(1:1)
     GOTO 100
  ELSE IF(CHAR.EQ.'+') THEN
     PATRES=PATRES+1
     WORD2=WORD(2:)
     WORD=WORD2
     CHAR=WORD(1:1)
     GOTO 100
  ELSE IF(CHAR.EQ.'#') THEN
     PATRES=PATRES+2
     WORD2=WORD(2:)
     WORD=WORD2
     CHAR=WORD(1:1)
     GOTO 100
  ENDIF
  !
  IF (PATRES.LT.ISTART.OR.PATRES.GT.ISTOP) THEN
     GOTO 900
  END IF
  I=IBASE(PATRES)+1
  DO WHILE (WORD /= ATYPE(I) .AND. I < IBASE(PATRES+1))
     I=I+1
  ENDDO
  IF (WORD == ATYPE(I)) THEN
     PATATM=I
     PTYP=WORD
  ELSE
     GOTO 900
  ENDIF
  RETURN
  !
900 CONTINUE
  PTYP=WORD
  PATATM=MARK
  IF(LWARN.AND. WORD.NE.' ') THEN
     !
     ! we dont want to have the following message if word=blank
     !
     IF(WRNLEV.GE.2) WRITE(OUTU,1000) ATOM(1:idleng),PATRES
1000 FORMAT(' ERROR in PATOM: ',A, &
          ' Atom type not found for residue',I5)
     CALL DIEWRN(2)
  ENDIF
  RETURN
  !
END SUBROUTINE PATOM

SUBROUTINE ATMINI(ISTART,ISTOP)
  !
  !     routine sets coordinate arrays and harmonic constraint array
  !     to default values for atom numbers between ISTART and ISTOP.
  !
  !     Axel Brunger, 17-MAY-83
  !
  use dimens_fcm
  use number
  use cnst_fcm
  use coord
  use coordc

  INTEGER ISTART,ISTOP
  INTEGER I
  !
  DO I=ISTART,ISTOP
     X(I)=ANUM
     Y(I)=ANUM
     Z(I)=ANUM
     WMAIN(I)=0.0
     XCOMP(I)=ANUM
     YCOMP(I)=ANUM
     ZCOMP(I)=ANUM
     WCOMP(I)=0.0
#if KEY_COMP2==1 /*comp2_1*/
     XCOMP2(I)=ANUM
     YCOMP2(I)=ANUM
     ZCOMP2(I)=ANUM
     WCOMP2(I)=0.0
#endif /* (comp2_1)*/
  ENDDO
  RETURN
END SUBROUTINE ATMINI

SUBROUTINE MKDRUD(NATOLD,NATOM,ATYPE,IAC,CG,AMASS,DMASS, &
     THOLEI,ALPHADP,ISDRUDE,nbond,ib,jb,NDRUDE,LSHOW)
  use dimens_fcm
  use param
  use number
  use consta
  use stream
  use aniso_fcm

  INTEGER I, NATOLD, NATOM, IAC(*)
  CHARACTER(len=*) ATYPE(*)
  real(chm_real) CG(*), AMASS(*), DMASS, THOLEI(*),ALPHADP(*)
  LOGICAL ISDRUDE(*), LSHOW, OK1, OK2, OK3, OK4
  integer j,nbond,ib(*),jb(*)
  integer NDRUDE
  integer numberl, i1,j1, jjj
  real(chm_real) KDRUDE, CGDRUDE
  integer icount
  integer IANISO, ITEMP

  icount = 0

  !     write(*,*) 'NCB ',ncb
  !     do i=1,ncb
  !     I1=SQRT(two*KCB(i))+half
  !     J1=KCB(i)-I1*(I1-1)/2
  !     write(*,*) i,kcb(j),cbb(j),cbc(j)
  !     WRITE(*,70) I,ATC(i1),ATC(j1),KCB(I),CBC(I),CBB(I),ICBCNT(I)
  ! 70  FORMAT(6X,I4,2X,A4,' - ',A4,3X,I7,F10.1,F10.3,I8)
  !     enddo

  if(prnlev.ge.5)then
     write(outu,*)
     write(outu,'(1x,a)') &
          'MKDRUDE generate list and setup for drude polarizability'
     write(outu,*)
  endif

  do i = natold+1,natom
     if(i.eq.1) cycle  ! because  of ALPHADP(i-1) argument

     ! By default atoms are not drudes...
     isdrude(i)=.false.

     !     OK1 = (atype(i)(1:1).eq.'D').and.(amass(i).eq.0.0)
     !    &     .and.(atc(iac(i))(1:3).eq.'DRU')

     OK1 = (atype(i)(1:1).eq.'D').and.(amass(i).eq.0.0) &
          .and.(abs(ALPHADP(i-1)).ne.0.0)

     if(OK1)then
        isdrude(i)=.true.

        i1=iac(i)
        j1=iac(i-1)

        if(prnlev.gt.5) write(outu,'(2(a,i5))') &
             'search for bond between atom ',i,' and ',i-1

        OK2 = .false.
        do j=1,nbond
           if(ib(j).eq.i .and. jb(j).eq.i-1) OK2=.true.
           if(ib(j).eq.i-1 .and. jb(j).eq.i) OK2=.true.
        enddo
        if(prnlev.gt.5) then
           if(OK2) write(outu,'(a)') 'bond found in psf '
        endif

        if(i1.gt.j1)then
           numberl = ((i1-1)*i1)/2 + j1
        else
           numberl = ((j1-1)*j1)/2 + i1
        endif
        jjj = 0
        do j=1,ncb
           if(numberl.eq.kcb(j)) jjj = j
        enddo

        OK3= jjj.ne.0

        if(prnlev.gt.8) &
             write(outu,'(3(a,2x),2f8.3)') &
             'bond type ',ATC(i1),ATC(j1), &
             cbc(jjj),cbb(jjj)

        KDRUDE = cbc(jjj)

        OK4 = (KDRUDE .ne. 0.0) .and. (cbb(jjj) .eq. 0.0)

        if (OK1 .and. OK2 .and. OK3 .and. OK4) then
           CGDRUDE = int(SQRT(2*ABS(ALPHADP(i-1))*KDRUDE/CCELEC)*10000.0 + 0.5)/10000.0
!           CGDRUDE  = SQRT(2*ABS(ALPHADP(i-1))*KDRUDE/CCELEC)
           CGDRUDE  = SIGN(CGDRUDE,ALPHADP(i-1)) ! copy sign of alpha
           CG(i)    = CGDRUDE                    ! set drude charge
           CG(i-1)  = CG(i-1)-CGDRUDE            ! set heavy atom charge
           AMASS(i) = DMASS
           AMASS(i-1) = AMASS(i-1)-DMASS

           !           CONVERTS K11,K22,K33 TO KPAR,KPERP AND KISO FOR USE IN
           !           EANISOTROPY ROUTINE (E. Harder 2007)

           do IANISO = 1, NANISO
              if (LSTANI1(IANISO)+1 .eq. i) then
                 K11(IANISO) = KDRUDE*K11(IANISO)
                 K22(IANISO) = KDRUDE*K22(IANISO)
                 K33(IANISO) = KDRUDE*K33(IANISO)
                 K33(IANISO) = K33(IANISO) - KDRUDE
                 K11(IANISO) = K11(IANISO) - KDRUDE - K33(IANISO)
                 K22(IANISO) = K22(IANISO) - KDRUDE - K33(IANISO)
                 if(prnlev.ge.8)then
                    write(outu,'(1x,a,4i6,3f10.3)') &
                         'anisotropy: ', &
                         lstani1(ianiso), &
                         lstani2(ianiso), &
                         lstani3(ianiso), &
                         lstani4(ianiso), &
                         k11(ianiso), &
                         k22(ianiso), &
                         k33(ianiso)
                 endif
              endif
           enddo
           !
           !

           if(prnlev.ge.8)then
              write(outu,'(1x,a,i6,2(a,f8.3))') &
                   'All okay for drude particle ',i, &
                   ' set its charge with KDRUDE =',kdrude, &
                   ' and its mass to ',DMASS
           endif
           icount = icount + 1

        else

           CALL WRNDIE(-5,'<MKDRUD>','COULD NOT SETUP DRUDES')

        endif

     endif

  enddo

#if KEY_PARALLEL==1
  CALL PSND8(AMASS,NATOM)
  CALL PSND8(CG,NATOM)
#endif 
  if(LSHOW.or.(prnlev.ge.4))then
     if (prnlev > 0) write(outu,'(1x,a,i6,a,/)') &
          'All okay for ',icount,' added Drude particles '
  endif

  !     NDRUDE used in fitcharge.src dcntrl.src and dynamvv2.src
  ndrude = ndrude + icount

  IF(LSHOW.or.(prnlev.gt.8))then
     do i = natold+1,natom
        if (prnlev > 0) write(outu,100) i, atype(i), iac(i), atc(iac(i)),cg(i),  &
             amass(i),ALPHADP(i), isdrude(i)
100     format(i4,3x,a5,i4,2x,a5,3f10.3,l3)
     enddo
     do i = natold+1,natom
        do IANISO = 1, NANISO
           if (LSTANI1(IANISO)+1 .eq. i) then
              if (prnlev > 0) write(outu,'(1x,a,5i6,3f10.3)') &
                   'anisotropy: ',ianiso, &
                   lstani1(ianiso), &
                   lstani2(ianiso), &
                   lstani3(ianiso), &
                   lstani4(ianiso), &
                   k11(ianiso), &
                   k22(ianiso), &
                   k33(ianiso)
           endif
        enddo
     enddo
     if (prnlev > 0) write(outu,*)
  ENDIF
  RETURN
END SUBROUTINE MKDRUD

end module genpsf_m

