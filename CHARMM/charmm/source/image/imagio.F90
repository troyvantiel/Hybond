SUBROUTINE IMREAD(COMLYN,COMLEN,IUNIT,LPRINT)
  !
  !     THIS ROUTINE READS THE IMAGE DATA FILE FROM CARDS
  !
  !     By Bernard R. Brooks    9/83
  !
  use chm_kinds
  use dimens_fcm
  use vector
  use image
  use stream
  use ctitla
  use corsubs,only:fndrot,fndu
  use string
  implicit none
  !
  character(len=*) COMLYN
  INTEGER COMLEN,IUNIT
  LOGICAL LPRINT
  !
  real(chm_real) SCALE(3),U(9)
  LOGICAL DONE,EOF,INVERT
  !
  real(chm_real) PHI,PHIX
  real(chm_real) UA(3,3),UB(3,3),UC(3,3)
  real(chm_real) TA(3),TB(3),TC(3),TBX(3)
  INTEGER I,ITRANS,IPT
  real(chm_real) UNCHK,YT,ZT,XT,UX,DIST
  INTEGER JTRANS,JPT,K,J,L,II,JJ,IPTX,NT
  LOGICAL LOK
  character(len=4) WRD
  !
  !
  !     READ TRANSFORMATION INFORMATION FROM CARDS
  !
  CALL TRYORO(IUNIT,'FORMATTED')
  !
  CALL RDTITL(TITLEB,NTITLB,IUNIT,0)
  IF(PRNLEV.GE.2) WRITE(OUTU,13) IUNIT
13 FORMAT(/12X,' READING IMAGE FILE FROM UNIT',I5/)
  !
  NTRANS=0
  SCALE(1)=1.0
  SCALE(2)=1.0
  SCALE(3)=1.0
  !
  EOF=.FALSE.
  DONE=.FALSE.
5000 CONTINUE
  !
  CALL RDCMND(COMLYN,MXCMSZ,COMLEN,IUNIT,EOF,.TRUE., &
       LPRINT,'IMREAD> ')
  IF(EOF) DONE=.TRUE.
  !
  WRD=NEXTA4(COMLYN,COMLEN)
  IF (WRD.EQ.'IMAG') THEN
     !       Begin Procedure PROCESS-IMAGE-COMMAND
     !
     IF(LPRINT .AND. PRNLEV.GE.2 .AND. NTRANS.GT.0) THEN
        !         print out the transformation just finished
        WRITE(OUTU,39) NTRANS,IMNAME(NTRANS)
39      FORMAT(/' TRANSFORMATION',I4,'  NAME:',A8)
        IPT=(NTRANS-1)*12+1
        DO J=1,4
           IPTX=IPT+2
           WRITE(OUTU,6) (IMTRNS(K),K=IPT,IPTX)
6          FORMAT(20X,3F10.6)
           IPT=IPT+3
        ENDDO
     ENDIF
     !
     NT=NTRANS+1
     IF(NT.GT.MAXTRN) CALL WRNDIE(-2,'<IMREAD>', &
          'OVERFLOW IN NUMBER OF TRANSFORMATIONS')
     WRD=NEXTA4(COMLYN,COMLEN)
     DO I=1,NTRANS
        IF(WRD.EQ.IMNAME(I)) &
             CALL WRNDIE(-2,'<IMREAD>','DUPLICATE TRANSFORMATION NAME')
     ENDDO
     IMNAME(NT)=WRD
     IPT=(NT-1)*12
     DO I=1,12
        IMTRNS(IPT+I)=0.0
     ENDDO
     IMTRNS(IPT+1)=1.0
     IMTRNS(IPT+5)=1.0
     IMTRNS(IPT+9)=1.0
     NTRANS=NT
     !       End Procedure PROCESS-IMAGE-COMMAND
  ELSE IF (WRD.EQ.'DEFI') THEN
     !       Begin Procedure PROCESS-DEFINE-COMMAND
     INVERT=.FALSE.
6000 CONTINUE
     WRD=NEXTA4(COMLYN,COMLEN)
     IF (WRD.EQ.'   ') THEN
     ELSE IF (WRD.EQ.'INVE') THEN
        INVERT=.TRUE.
     ELSE 
        !
        !         FIND FORMER TRANSFORMATION
        !
        DO ITRANS=1,NTRANS
           IF(WRD.EQ.IMNAME(ITRANS)) GOTO 200
        ENDDO
        IF(WRNLEV.GE.2) WRITE(OUTU,192) WRD
192     FORMAT(/' TRANSFORMATION NAME ''',A4,''' NOT FOUND')
        CALL WRNDIE(-2,'<IMREAD>', &
             'NAME NOT FOUND, RESET TO FIRST ONE')
        ITRANS=1
200     CONTINUE
        IPT=(ITRANS-1)*12
        DO I=1,3
           DO J=1,3
              IPT=IPT+1
              UB(J,I)=IMTRNS(IPT)
           ENDDO
        ENDDO
        !
        DO I=1,3
           IPT=IPT+1
           TB(I)=IMTRNS(IPT)
        ENDDO
        IF(INVERT) THEN
           !
           !           INVERT TRANSLATION VECTOR (AND ROTATE)
           DO I=1,3
              TC(I)=0.0
              DO J=1,3
                 TC(I)=TC(I)-TB(J)*UB(I,J)
              ENDDO
           ENDDO
           !           TC(1)=-(TB(1)*UB(1,1)+TB(2)*UB(1,2)+TB(3)*UB(1,3))
           !           TC(2)=-(TB(1)*UB(2,1)+TB(2)*UB(2,2)+TB(3)*UB(2,3))
           !           TC(3)=-(TB(1)*UB(3,1)+TB(2)*UB(3,2)+TB(3)*UB(3,3))
           TB(1)=TC(1)
           TB(2)=TC(2)
           TB(3)=TC(3)
           !
           !           INVERT ROTATION MATRIX
           DO I=1,3
              DO J=1,I
                 UX=UB(J,I)
                 UB(J,I)=UB(I,J)
                 UB(I,J)=UX
              ENDDO
           ENDDO
        ENDIF
        !
        !         MULTIPLY-TRANSFORMATIONS
        CALL MULTTR(TA,TB,TC,UA,UB,UC)
        !
        INVERT=.FALSE.
     ENDIF
     IF (WRD.NE.'    ') GOTO 6000
     !       End Procedure PROCESS-DEFINE-COMMAND
  ELSE IF (WRD.EQ.'SCAL') THEN
     !       Begin Procedure PROCESS-SCALE-COMMAND
     SCALE(1)=NEXTF(COMLYN,COMLEN)
     SCALE(2)=NEXTF(COMLYN,COMLEN)
     SCALE(3)=NEXTF(COMLYN,COMLEN)
     !       End Procedure PROCESS-SCALE-COMMAND
  ELSE IF (WRD.EQ.'ROTA') THEN
     !       Begin Procedure PROCESS-ROTA-COMMAND
     DO I=1,3
        TB(I)=NEXTF(COMLYN,COMLEN)
     ENDDO
     CALL NORMALL(TB,3)
     PHI=NEXTF(COMLYN,COMLEN)
     CALL FNDU(UB,TB,PHI,LOK)
     IF (.NOT.LOK) THEN
        !         Begin Procedure CRAP-OUT
        IF(WRNLEV.GE.2) WRITE(OUTU,925)
925     FORMAT(/' ***** ERROR ***** <IMREAD> PARSING ERROR.'/)
        CALL DIEWRN(0)
        RETURN
        !         End Procedure CRAP-OUT
     ENDIF
     !
     TB(1)=0.0
     TB(2)=0.0
     TB(3)=0.0
     !       MULTIPLY-TRANSFORMATIONS
     CALL MULTTR(TA,TB,TC,UA,UB,UC)
     IF(LPRINT .AND. PRNLEV.GE.2) THEN
        WRITE(OUTU,61) UB
61      FORMAT(' ROTATION MATRIX'/,3(3F12.6/))
        CALL FNDROT(UB,TBX,PHIX,LPRINT)
     ENDIF
     !       End Procedure PROCESS-ROTA-COMMAND
  ELSE IF (WRD.EQ.'TRAN') THEN
     !       Begin Procedure PROCESS-TRAN-COMMAND
     TB(1)=NEXTF(COMLYN,COMLEN)
     TB(2)=NEXTF(COMLYN,COMLEN)
     TB(3)=NEXTF(COMLYN,COMLEN)
     !
     CALL TRIME(COMLYN,COMLEN)
     IF(COMLEN.GT.0) THEN
        DIST=NEXTF(COMLYN,COMLEN)
        CALL NORMALL(TB,3)
        TB(1)=TB(1)*DIST
        TB(2)=TB(2)*DIST
        TB(3)=TB(3)*DIST
     ENDIF
     !
     TB(1)=TB(1)*SCALE(1)
     TB(2)=TB(2)*SCALE(2)
     TB(3)=TB(3)*SCALE(3)
     IPT=(NTRANS-1)*12+9
     IMTRNS(IPT+1)=IMTRNS(IPT+1)+TB(1)
     IMTRNS(IPT+2)=IMTRNS(IPT+2)+TB(2)
     IMTRNS(IPT+3)=IMTRNS(IPT+3)+TB(3)
     IF(LPRINT .AND. PRNLEV.GE.2) WRITE(OUTU,71) TB
71   FORMAT(' TRANSLATION VECTOR ',3F12.6)
     !       End Procedure PROCESS-TRAN-COMMAND
  ELSE IF (WRD.EQ.'NEGA') THEN
     !       Begin Procedure PROCESS-NEGATE-COMMAND
     DO I=1,3
        DO J=1,3
           UB(I,J)=0.0
        ENDDO
        UB(I,I)=-1.0
        TB(I)=0.0
     ENDDO
     !       MULTIPLY-TRANSFORMATIONS
     CALL MULTTR(TA,TB,TC,UA,UB,UC)
     !       End Procedure PROCESS-NEGATE-COMMAND
  ELSE IF (WRD.EQ.'END ') THEN
     !       Begin Procedure PROCESS-END-COMMAND
     DONE=.TRUE.
     !       End Procedure PROCESS-END-COMMAND
  ELSE IF (WRD.EQ.'    ') THEN
  ELSE
     CALL WRNDIE(0,'<IMREAD>','Unrecognized command')
  ENDIF
  !
  CALL XTRANE(COMLYN,COMLEN,'IMREAD')
  !
  IF (.NOT.DONE) GOTO 5000
  !
  !     Begin Procedure CHECK-ALL-TRANSFORMATIONS
  !
  DO ITRANS=1,NTRANS
     IPT=(ITRANS-1)*12
     UNCHK=0.0
     UNCHK=UNCHK+IMTRNS(IPT+1)*IMTRNS(IPT+5)*IMTRNS(IPT+9)
     UNCHK=UNCHK-IMTRNS(IPT+1)*IMTRNS(IPT+6)*IMTRNS(IPT+8)
     UNCHK=UNCHK+IMTRNS(IPT+2)*IMTRNS(IPT+6)*IMTRNS(IPT+7)
     UNCHK=UNCHK-IMTRNS(IPT+2)*IMTRNS(IPT+4)*IMTRNS(IPT+9)
     UNCHK=UNCHK+IMTRNS(IPT+3)*IMTRNS(IPT+4)*IMTRNS(IPT+8)
     UNCHK=UNCHK-IMTRNS(IPT+3)*IMTRNS(IPT+5)*IMTRNS(IPT+7)
     IF(ABS(ABS(UNCHK)-1.0).GT.1.D-4) THEN
        IF(WRNLEV.GE.2) WRITE(OUTU,85) ITRANS,UNCHK
85      FORMAT(/' WARNING: TRANSFORMATION MATRIX',I4, &
             ' IS NOT UNITARY', &
             /'      DETERMINANT IS  ',F14.8/)
        CALL DIEWRN(-2)
     ENDIF
     IF(UNCHK.LT.0.0) THEN
        IF(PRNLEV.GE.2) WRITE(OUTU,23) ITRANS,IMNAME(ITRANS)
23      FORMAT(' TRANSFORMATION',I5,2X,A4,' IS A MIRROR REFLECTION')
     ENDIF
  ENDDO
  !
  IF(PRNLEV.GE.2) WRITE(OUTU,125) NTRANS
125 FORMAT(/I6,'  TRANSFORMATIONS HAVE BEEN READ'/)
  NOROT=.TRUE.
  !
  loop9k:DO ITRANS=1,NTRANS
     IPT=(ITRANS-1)*12
     loop8k:DO JTRANS=1,NTRANS
        JPT=(JTRANS-1)*12
        YT=IMTRNS(IPT+10)*IMTRNS(JPT+4)+IMTRNS(IPT+11)*IMTRNS(JPT+5)+ &
             IMTRNS(IPT+12)*IMTRNS(JPT+6)+IMTRNS(JPT+11)
        ZT=IMTRNS(IPT+10)*IMTRNS(JPT+7)+IMTRNS(IPT+11)*IMTRNS(JPT+8)+ &
             IMTRNS(IPT+12)*IMTRNS(JPT+9)+IMTRNS(JPT+12)
        XT=IMTRNS(IPT+10)*IMTRNS(JPT+1)+IMTRNS(IPT+11)*IMTRNS(JPT+2)+ &
             IMTRNS(IPT+12)*IMTRNS(JPT+3)+IMTRNS(JPT+10)
        UNCHK=XT*XT+YT*YT+ZT*ZT
        IF(UNCHK.LE.1.D-4) THEN
           !
           K=0
           DO I=1,3
              DO J=1,3
                 K=K+1
                 U(K)=0.0
                 DO L=1,3
                    II=IPT+3*J+L-3
                    JJ=JPT+3*L+I-3
                    U(K)=U(K)+IMTRNS(II)*IMTRNS(JJ)
                 ENDDO
              ENDDO
           ENDDO
           !
           U(1)=U(1)-1.0
           U(5)=U(5)-1.0
           U(9)=U(9)-1.0
           UNCHK=0.0
           DO I=1,9
              UNCHK=UNCHK+U(I)*U(I)
           ENDDO
           IF(UNCHK.LE.0.001) THEN
              IMINV(ITRANS)=JTRANS
              GOTO 170
           ENDIF
        ENDIF
     enddo loop8k
     !
     IF(WRNLEV.GE.2) WRITE(OUTU,185) ITRANS
185  FORMAT(/' WARNING: TRANSFORMATION NUMBER',I4, &
          '  HAS NO INVERSE'/)
     IMINV(ITRANS)=0
     CALL DIEWRN(-3)
     !
170  CONTINUE
     DO I=1,9
        J=I+IPT
        U(I)=IMTRNS(J)
     ENDDO
     U(1)=U(1)-1.0
     U(5)=U(5)-1.0
     U(9)=U(9)-1.0
     UNCHK=0.0
     DO I=1,9
        UNCHK=UNCHK+U(I)*U(I)
     ENDDO
     IF(UNCHK.GT.0.001) NOROT=.FALSE.
     !
  enddo loop9k
  !
  IF(PRNLEV.GE.2) THEN
     IF (NOROT) THEN
        WRITE(OUTU,'(20X,A)') &
             ' THERE ARE NO ROTATIONS FOR THIS TRANSFORMATION SET'
     ELSE
        WRITE(OUTU,'(20X,A)') &
             ' THERE ARE ROTATIONS FOR THIS TRANSFORMATION SET'
     ENDIF
  ENDIF
  !
  !     End Procedure CHECK-ALL-TRANSFORMATIONS
  RETURN
  !
END SUBROUTINE IMREAD

SUBROUTINE MULTTR(TA,TB,TC,UA,UB,UC)
  !
  !     TO MULTIPLY-TRANSFORMATIONS
  !
  !     UA  =  UB  X  UC
  !
  use chm_kinds
  use dimens_fcm
  use number
  !
  use image
  implicit none
  !
  real(chm_real) UA(3,3),UB(3,3),UC(3,3)
  real(chm_real) TA(3),TB(3),TC(3)
  !
  !     local
  !
  INTEGER I,IPT,J,K
  !
  IPT=(NTRANS-1)*12
  DO I=1,3
     DO J=1,3
        IPT=IPT+1
        UC(J,I)=IMTRNS(IPT)
     ENDDO
  ENDDO
  !
  DO I=1,3
     IPT=IPT+1
     TC(I)=IMTRNS(IPT)
  ENDDO
  !
  DO I=1,3
     DO J=1,3
        UA(J,I)=ZERO
        DO K=1,3
           UA(J,I)=UA(J,I)+UB(K,I)*UC(J,K)
        ENDDO
     ENDDO
  ENDDO
  !
  DO I=1,3
     TA(I)=TB(I)
     DO K=1,3
        TA(I)=TA(I)+UB(K,I)*TC(K)
     ENDDO
  ENDDO
  !
  IPT=IPT-12
  DO I=1,3
     DO J=1,3
        IPT=IPT+1
        IMTRNS(IPT)=UA(J,I)
     ENDDO
  ENDDO
  DO I=1,3
     IPT=IPT+1
     IMTRNS(IPT)=TA(I)
  ENDDO
  !
  RETURN
END SUBROUTINE MULTTR

SUBROUTINE IMWRIT(IUNIT,COMLYN,COMLEN,BIMAG,IPRINT)
  !
  !     THIS ROUTINE WRITES IN CARD IMAGE OR PRINTS THE IMAGE DATA FILE
  !
  !     By Bernard R. Brooks    9/83
  !
  use chm_kinds
  use chm_types
  use dimens_fcm
  use exfunc
  use psf
  use image
  use stream
  use ctitla
  use string

  !
  implicit none
  !
  character(len=*) COMLYN
  INTEGER IUNIT,COMLEN,IPRINT
  type(imageDataStructure) BIMAG
  !
  INTEGER I,ITRANS,IPT,J,IPTX,K
  character(len=4) WRD
  real(chm_real) TOTF(3),TOTT(3)
  !
  WRD=NEXTA4(COMLYN,COMLEN)
  IF (WRD.EQ.'PSF ') THEN
     !       Begin Procedure PRINT-IMAGE-PSF
     CALL IMPSFW(IUNIT,IPRINT,NTRANS,IMNAME,IMINV, &
          BIMAG%IMATPT, &
          RES,RESID,IBASE,ATYPE,CG,AMASS,IAC, &
          IMOVE,IB(NBOND+1),JB(NBOND+1), &
          IT(NTHETA+1),JT(NTHETA+1),KT(NTHETA+1), &
          IP(NPHI+1),JP(NPHI+1),KP(NPHI+1),LP(NPHI+1), &
          IM(NIMPHI+1),JM(NIMPHI+1),KM(NIMPHI+1),LM(NIMPHI+1), &
#if KEY_CMAP==1
          I1CT(NCRTERM+1),J1CT(NCRTERM+1),K1CT(NCRTERM+1),L1CT(NCRTERM+1), & 
#endif
#if KEY_CMAP==1
          I2CT(NCRTERM+1),J2CT(NCRTERM+1),K2CT(NCRTERM+1),L2CT(NCRTERM+1), & 
#endif
#if KEY_CMAP==1
          NIMCRT, &   
#endif
          BIMAG%NIMINB,BIMAG%IMINB,BIMAG%IMIBLO, &
          NIMRES,NIMGRP,NATIM,NIMBON,NIMANG,NIMDIH,NIMIMP, &
          NRES,NATOM,NGRP,IGPBS,IGPTYP,IMOVEG)
     !       End Procedure PRINT-IMAGE-PSF
  ELSE IF (WRD.EQ.'FORC') THEN
     !       Begin Procedure PRINT-IMAGE-FORCES
     !
     DO I=1,3
        TOTT(I)=0.0
        TOTF(I)=0.0
     ENDDO
     !
     DO ITRANS=1,NTRANS
        IPT=(ITRANS-1)*3
        WRITE(IUNIT,734) ITRANS,IMNAME(ITRANS), &
             (IMFORC(IPT+J),J=1,3),(IMTORQ(IPT+J),J=1,3)
734     FORMAT(' TRANSFORMATION',I5,2X,A4/, &
             '   FORCES:',3F16.6/'   TORQUE:',3F16.6)
        DO J=1,3
           TOTT(J)=TOTT(J)+IMTORQ(IPT+J)
           TOTF(J)=TOTF(J)+IMFORC(IPT+J)
        ENDDO
     ENDDO
     !
     WRITE(IUNIT,735) (TOTF(J),J=1,3),(TOTT(J),J=1,3)
735  FORMAT(' TOTAL'/,'   FORCES:',3F16.6/'   TORQUE:',3F16.6)
     !       End Procedure PRINT-IMAGE-FORCES
  ELSE IF (WRD.EQ.'TRAN') THEN
     !       Begin Procedure PRINT-IMAGE-TRANSFORMATIONS
     CALL WRTITL(TITLEA,NTITLA,IUNIT,1)
     WRITE(IUNIT,37) NTRANS
37   FORMAT(/' NUMBER OF TRANSFORMATIONS IS',I4)
     IPT=1
     DO I=1,NTRANS
        WRITE(IUNIT,39) I,IMNAME(I),IMINV(I)
39      FORMAT(/' TRANSFORMATION',I4,'  NAME:',A8,'  INVERSE',I4)
        DO J=1,4
           IPTX=IPT+2
           WRITE(IUNIT,6) (IMTRNS(K),K=IPT,IPTX)
6          FORMAT(20X,3F10.6)
           IPT=IPT+3
        ENDDO
     ENDDO
     !       End Procedure PRINT-IMAGE-TRANSFORMATIONS
  ELSE IF (WRD.EQ.'    ') THEN
     CONTINUE
  ELSE
     IF(WRNLEV.GE.0) WRITE(OUTU,77) WRD
77   FORMAT('IMWRIT: Unable to write IMAGE file type: "',A4,'".')
     CALL WRNDIE(0,'<IMWRIT>','Unrecognized IMAGE file type')
  ENDIF
  !
  CALL XTRANE(COMLYN,COMLEN,'IMWRIT')
  RETURN
  !
END SUBROUTINE IMWRIT

SUBROUTINE IMPSFW(IUNIT,IPRINT,NTRANS,IMNAME,IMINV,IMATPT, &
     RES,RESID,IBASE,ATYPE,CG,AMASS,IAC, &
     IMOVE,IB,JB,IT,JT,KT,IP,JP,KP,LP,IM,JM,KM,LM, &
#if KEY_CMAP==1
     I1CT,J1CT,K1CT,L1CT,I2CT,J2CT,K2CT,L2CT,NIMCRT, &   
#endif
     NIMINB,IMINB,IMIBLO, &
     NIMRES,NIMGRP,NATIM,NIMBON,NIMANG,NIMDIH,NIMIMP, &
     NRES,NATOM,NGRP,IGPBS,IGPTYP,IMOVEG)
  !
  !     THIS ROUTINE WRITES THE PSF FOR IMAGE ATOMS.
  !
  use chm_kinds
  use exfunc
  use stream
  implicit none
  !
  INTEGER IUNIT,IPRINT,NTRANS,IMINV(*)
  character(len=*) IMNAME(*)
  INTEGER     IMATPT(*)
  real(chm_real)      CG(*),AMASS(*)
  character(len=*) RES(*),RESID(*),ATYPE(*)
  INTEGER     IBASE(*),IAC(*)
  INTEGER     IMOVE(*),IB(*),JB(*),IT(*),JT(*),KT(*)
  INTEGER     IP(*),JP(*),KP(*),LP(*),IM(*),JM(*),KM(*),LM(*)
  INTEGER     IMINB(*),IMIBLO(*)
#if KEY_CMAP==1
  INTEGER I1CT(*),J1CT(*),K1CT(*),L1CT(*)
  INTEGER I2CT(*),J2CT(*),K2CT(*),L2CT(*)
  INTEGER NIMCRT
#endif 
  INTEGER NIMRES(*),NIMGRP,NATIM,NIMBON,NIMANG,NIMDIH,NIMIMP
  INTEGER NRES,NATOM,NGRP,NIMINB
  INTEGER IGPBS(*),IGPTYP(*),IMOVEG(*)
  !
  INTEGER ITRANS,NEXTT,IRES,IGRP,I,INV,IS,IW,J

  ! IF(IPRINT.EQ.0) CALL WRNDIE(-5,'<IMPSFW>','Bad print option')

  WRITE(IUNIT,40) NATIM,NIMBON,NIMANG,NIMDIH,NIMIMP, &
#if KEY_CMAP==1
       NIMCRT, &   
#endif
       NIMINB,NIMRES(NTRANS),NTRANS,NIMGRP

#if KEY_CMAP==1
40 FORMAT(/9X,'NATIM  NIMBON NIMANG NIMDIH NIMIMP NIMCRT ', &
       'NIMINB NIMRES NTRANS NIMGRP'/7X,13I7)
#else /**/
40 FORMAT(/9X,'NATIM  NIMBON NIMANG NIMDIH NIMIMP NIMINB ', &
       'NIMRES NTRANS NIMGRP'/7X,12I7)
#endif 

  WRITE(IUNIT,105)
105 FORMAT(/10X,'ATOM CHARACTERISTICS :'/19X,'ATOM TYPE      ', &
       'CHARGE    ATOM CODE   COUNT OF    MOVEMENT FLAG', &
       '        MASS'/56X,'EXCLUSIONS'/)

  ITRANS=0
  NEXTT=0
  IRES=NRES+1
  IGRP=NGRP+1

50 FORMAT(' TRANSFORMATION',I4,2X,A4,'  INVERSE',I4,2X,A4)
108 FORMAT('     RESIDUE',I4,2X,A,2X,A,8X,'  TO',I6)
109 FORMAT(10X,'GROUP',I4,'  TYPE',I3,' MOVE',I2,'  TO',I6)
110 FORMAT(20X,I5,2X,A,F12.4,I5,I6,I5,1PG14.6)
  DO I=NATOM+1,NATIM
     IF(I.GT.NEXTT) THEN
        ITRANS=ITRANS+1
        NEXTT=IMATPT(ITRANS)
        INV=IMINV(ITRANS)
        WRITE(IUNIT,50) ITRANS,IMNAME(ITRANS),INV,IMNAME(INV)
     ENDIF
     IF(IBASE(IRES)+1.EQ.I) THEN
        WRITE(IUNIT,108) IRES,RESID(IRES)(1:idleng), &
             RES(IRES)(1:idleng),IBASE(IRES+1)
        IRES=IRES+1
     ENDIF
     IF(IGPBS(IGRP)+1.EQ.I) THEN
        WRITE(IUNIT,109) IGRP,IGPTYP(IGRP),IMOVEG(IGRP),IGPBS(IGRP+1)
        IGRP=IGRP+1
     ENDIF
     WRITE(IUNIT,110) I,ATYPE(I)(1:idleng),CG(I),IAC(I),IMIBLO(I), &
          IMOVE(I),AMASS(I)
  ENDDO

  WRITE(IUNIT,'(/10X,A)') 'BOND ARRAY (BY COLUMNS) :'
  DO I=1,NIMBON,20
     IS=I
     IW=I+19
     IF(IW.GT.NIMBON) IW=NIMBON
     WRITE(IUNIT,70) I,(IB(J),J=IS,IW)
     WRITE(IUNIT,160) (JB(J),J=IS,IW)
  ENDDO

  WRITE(IUNIT,'(/10X,A)') 'THETA ARRAY (BY COLUMNS) :'
  DO I=1,NIMANG,20
     IS=I
     IW=I+19
     IF(IW.GT.NIMANG) IW=NIMANG
     WRITE(IUNIT,70) I,(IT(J),J=IS,IW)
     WRITE(IUNIT,160) (JT(J),J=IS,IW)
     WRITE(IUNIT,160) (KT(J),J=IS,IW)
  ENDDO

  WRITE(IUNIT,'(/10X,A)') 'PHI ARRAY (BY COLUMNS) :'
  DO I=1,NIMDIH,20
     IS=I
     IW=I+19
     IF(IW.GT.NIMDIH) IW=NIMDIH
     WRITE(IUNIT,70) I,(IP(J),J=IS,IW)
     WRITE(IUNIT,160) (JP(J),J=IS,IW)
     WRITE(IUNIT,160) (KP(J),J=IS,IW)
     WRITE(IUNIT,160) (LP(J),J=IS,IW)
  ENDDO
  !
  IF(NIMIMP.EQ.0) GOTO 235
  WRITE(IUNIT,'(/10X,A)') 'IMPROPER TORSION ARRAY (BY COLUMNS) :'
  DO I=1,NIMIMP,20
     IS=I
     IW=I+19
     IF(IW.GT.NIMIMP) IW=NIMIMP
     WRITE(IUNIT,70) I,(IM(J),J=IS,IW)
     WRITE(IUNIT,160) (JM(J),J=IS,IW)
     WRITE(IUNIT,160) (KM(J),J=IS,IW)
     WRITE(IUNIT,160) (LM(J),J=IS,IW)
  ENDDO

#if KEY_CMAP==1
  IF(NIMCRT.EQ.0) GOTO 235
  WRITE(IUNIT,'(/10X,A)') 'CROSSTERM ARRAY (BY COLUMNS) :'
  DO I=1,NIMCRT,20
     IS=I
     IW=I+19
     IF(IW.GT.NIMCRT) IW=NIMCRT
     WRITE(IUNIT,70) I,(I1CT(J),J=IS,IW)
     WRITE(IUNIT,160) (J1CT(J),J=IS,IW)
     WRITE(IUNIT,160) (K1CT(J),J=IS,IW)
     WRITE(IUNIT,160) (L1CT(J),J=IS,IW)
     WRITE(IUNIT,160) (I2CT(J),J=IS,IW)
     WRITE(IUNIT,160) (J2CT(J),J=IS,IW)
     WRITE(IUNIT,160) (K2CT(J),J=IS,IW)
     WRITE(IUNIT,160) (L2CT(J),J=IS,IW)
  ENDDO
#endif 

70 FORMAT(10X,I5,5X,20I5)
160 FORMAT(20X,20I5)

235 WRITE(IUNIT,'(/10X,A)') 'NON-BONDED EXCLUSION ARRAY :'
  DO I=1,NIMINB,20
     IS=I
     IW=I+19
     IF(IW.GT.NIMINB) IW=NIMINB
     WRITE(IUNIT,70) I,(IMINB(J),J=IS,IW)
  ENDDO

  RETURN
END SUBROUTINE IMPSFW


