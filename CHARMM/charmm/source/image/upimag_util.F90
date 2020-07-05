module upimag_util
  implicit none
  private

  public upimnb, mkimnb

contains

  SUBROUTINE UPIMNB(BIMAG)
    !-----------------------------------------------------------------------
    !     GENERATE EXCLUSION LISTS FOR IMAGE-PRIMARY ATOMS IF PRESENT
    !
    use chm_kinds
    use chm_types
    !  use dimens_fcm
    use number
    !
    use machdep
    use psf
    use inbnd
    use image
    use memory
    use stream
    use pbound
    use datstr, only: imgrow
    use nbexcl_util,only:makgrp
    implicit none
    !
    !-----------------------------------------------------------------------
    !  Miscellaneous:
    !
    !  MAXING - The maximum number of atoms in any electrostatic group.
    !
    integer,parameter :: MAXING=1000
    !-----------------------------------------------------------------------
    integer,allocatable,dimension(:,:) :: IATBON
    integer,allocatable,dimension(:) :: NATBON
    integer,allocatable,dimension(:) :: IATTR
    integer,allocatable,dimension(:) :: IPRA
    integer,allocatable,dimension(:) :: IIMA
    integer,allocatable,dimension(:) :: IPRB
    integer,allocatable,dimension(:) :: IIMB
    integer,allocatable,dimension(:) :: IGR, IPT, ISTOP

    type(imageDataStructure) BIMAG
    !
    INTEGER I,MAXNNG,IATMXB
    LOGICAL CMPLTD
    !
    !
    IF(NTRANS.EQ.0) RETURN
    IF(NBXMOD.LT.-5 .OR. NBXMOD.GT.5) RETURN
    !
#if KEY_PBOUND==1
    If (qBoun) Return
#endif
    !
    !     GENERATE EXCLUSION LISTS FOR IMAGE-PRIMARY ATOMS IF PRESENT
    !
    NBONDT=NBOND+NIMBON
    I=3*NBONDT
    IF(IABS(NBXMOD).GT.3) I=I*2
    IF(NBXMOD.GT.0) I=I+NNB+2

    call IMGROW(bimag%IMINB, I)

    call chmalloc('upimag_util.src','UPIMNB','NATBON',NATOMT,intg=NATBON)
    IATMXB=8
    IF(NATOMT.LT.200) IATMXB=20
    call chmalloc('upimag_util.src','UPIMNB','IATBON',IATMXB,NATOMT,intg=IATBON)
    !
    call chmalloc('upimag_util.src','UPIMNB','IATTR',NATOM,intg=IATTR)
    call chmalloc('upimag_util.src','UPIMNB','IPRA',IATMXB,intg=IPRA)
    call chmalloc('upimag_util.src','UPIMNB','IIMA',IATMXB,intg=IIMA)
    call chmalloc('upimag_util.src','UPIMNB','IPRB',IATMXB*IATMXB,intg=IPRB)
    call chmalloc('upimag_util.src','UPIMNB','IIMB',IATMXB*IATMXB,intg=IIMB)
    !
    DO
       CALL MKIMNB(NATOM,IB,JB,NBOND,NBONDT,NATIM, &
            NTRANS,BIMAG%IMATPT,BIMAG%IMATTR,IATTR, &
            IPRA,IIMA,IPRB,IIMB, &
            BIMAG%IMINB, &
            BIMAG%IMIBLO,BIMAG%NIMINB,(I), &
            NBXMOD,NATBON,IATBON,IATMXB, &
            CMPLTD)
       IF (CMPLTD) EXIT
       call IMGROW(bimag%IMINB, 1.5, 10)
    ENDDO

    call chmdealloc('upimag_util.src','UPIMNB','IATBON',IATMXB,NATOMT,intg=IATBON)

    ! generate image group exclusion lists

    maxnng = BIMAG%NIMINB / 2 + 10
    call IMGROW(bimag%IMING, maxnng)

    do
       call chmalloc('upimag_util.src','UPIMNB','IGR', natom, intg=IGR)
       call chmalloc('upimag_util.src','UPIMNB','IPT', MAXING, intg=IPT)
       call chmalloc('upimag_util.src','UPIMNB','ISTOP', MAXING, intg=ISTOP)

       CALL MAKGRP(NGRP+1,NGRPT,NGRP,BIMAG%IMIGLO, &
            BIMAG%IMING,BIMAG%NIMING, size(bimag%iming), &
            BIMAG%IMIBLO,BIMAG%IMINB,IGPBS,IGR, &
            CMPLTD,IPT,ISTOP)

       call chmdealloc('upimag_util.src','UPIMNB','IGR', natom, intg=IGR)
       call chmdealloc('upimag_util.src','UPIMNB','IPT',MAXING,intg=IPT)
       call chmdealloc('upimag_util.src','UPIMNB','ISTOP',MAXING,intg=ISTOP)

       IF (CMPLTD) EXIT
       call IMGROW(bimag%IMING, 1.5, 10)
    end do

    call chmdealloc('upimag_util.src','UPIMNB','NATBON',NATOMT,intg=NATBON)
    call chmdealloc('upimag_util.src','UPIMNB','IATTR',NATOM,intg=IATTR)
    call chmdealloc('upimag_util.src','UPIMNB','IPRA',IATMXB,intg=IPRA)
    call chmdealloc('upimag_util.src','UPIMNB','IIMA',IATMXB,intg=IIMA)
    call chmdealloc('upimag_util.src','UPIMNB','IPRB',IATMXB*IATMXB,intg=IPRB)
    call chmdealloc('upimag_util.src','UPIMNB','IIMB',IATMXB*IATMXB,intg=IIMB)
    RETURN
  END SUBROUTINE UPIMNB

  SUBROUTINE MKIMNB(NATOM,IB,JB,NBOND,NBONDT,NATIM,NTRANS, &
       IMATPT,IMATTR,IATTR,IPRA,IIMA,IPRB,IIMB, &
       INB14,IBLO14,NNB14,MXNB14,MODE,NATBON, &
       IATBON,IATBMX,CMPLTD)
    !-----------------------------------------------------------------------
    !     THIS ROUTINE TAKES THE PSF INFORMATION AND GENERATES
    !     A NONBONDED EXCLUSION LIST (INB14 and IBLO14).
    !
    !     MODE =    0        LEAVE THE EXISTING INB14/IBLO14 LISTS
    !     MODE = +- 1        INCLUDE NOTHING
    !     MODE = +- 2        INCLUDE ONLY 1-2 (BOND) INTERACTIONS
    !     MODE = +- 3        INCLUDE 1-2 AND 1-3 (BOND AND ANGLE)
    !     MODE = +- 4        INCLUDE 1-2 1-3 AND 1-4's
    !     MODE = +- 5        INCLUDE 1-2 1-3 AND SPECIAL 1-4 interactions
    !     A POSITIVE MODE VALUE CAUSES THE EXISTING INB ARRAY TO BE ADDED.
    !
    !
    !
    !     By Bernard R. Brooks    August 1983
    !
    use chm_kinds
    use stream
    use exfunc
    implicit none
    !
    INTEGER NATOM,NBOND,NBONDT,NATIM,NTRANS
    INTEGER NNB14,MXNB14,MODE,IATBMX
    INTEGER IB(*),JB(*)
    INTEGER INB14(*),IBLO14(*)
    INTEGER NATBON(*),IATBON(IATBMX,*)
    INTEGER IMATTR(*),IMATPT(*),IATTR(*)
    INTEGER IPRA(*),IIMA(*),IPRB(*),IIMB(*)
    LOGICAL LEX14,CMPLTD,LFOUND
    EXTERNAL  EXCH5
    !
    !     LOCAL STORAGE
    INTEGER IPK(MXNB14),JPK(MXNB14),IPK14(MXNB14),JPK14(MXNB14)
    INTEGER MODEX,NPAIR,NPAIR4,MAXWRK,MXWRK4,NIMA,NPRA,NIMB,NPRB
    INTEGER NBONDP,IBOND,ITRANS,IATOM,IS,IQ,NX14,NEXT14
    INTEGER I,J,IJ,IPR,IBT,JBT,IIM,IA,NA,JA,I2,J2,I2L,J2L
    !
    !
    CMPLTD=.FALSE.
    MODEX=IABS(MODE)
    !
    IF(NBOND.EQ.NBONDT) THEN
       MODEX=0
       NNB14=0
    ENDIF
    !
    IF(MODEX.EQ.0) THEN
       IF(NNB14.EQ.0) THEN
          DO I=1,NATIM
             IBLO14(I)=0
          ENDDO
       ENDIF
       IF(NNB14.NE.IBLO14(NATIM)) THEN
          CALL WRNDIE(-1,'<MKIMNB>', &
               'NBXMod=0 not allowed without existing exclusion list')
       ENDIF
       CMPLTD=.TRUE.
       RETURN
    ENDIF
    IF(MODEX.GT.5) MODEX=5
    !
    NPAIR=0
    MAXWRK=MXNB14
    !
    IF(MODEX.EQ.5) THEN
       NPAIR4=0
       MXWRK4=MXNB14
    ENDIF
    !
    !     FOR MODE GREATER THAN ZERO INCLUDE THE EXISTING INB/IBLO LISTS,
    !     BUT NO IMAGE INTERACTIONS ARE ON EXCLUDED LIST.
    !
    !     COMPILE A LIST OF ALL THE SPECIFIED INTERACTIONS
    !
    !     This is based on the bond list.  First make a list of all bonds
    !     for every atom.
    !
    DO I=1,NATOM
       NATBON(I)=0
    ENDDO
    DO I=1,NBOND
       IBT=IB(I)
       JBT=JB(I)
       IF(IBT.GT.0 .AND. JBT.GT.0) THEN
          NATBON(IBT)=NATBON(IBT)+1
          IF(NATBON(IBT).GT.IATBMX) THEN
             IF(WRNLEV.GE.2) WRITE(OUTU,35) IBT
             CALL DIEWRN(-4)
          ENDIF
          IATBON(NATBON(IBT),IBT)=I
          NATBON(JBT)=NATBON(JBT)+1
          IF(NATBON(JBT).GT.IATBMX) THEN
             IF(WRNLEV.GE.2) WRITE(OUTU,35) JBT
             CALL DIEWRN(-4)
          ENDIF
          IATBON(NATBON(JBT),JBT)=-I
       ENDIF
    ENDDO
35  FORMAT(' <MKIMNB>: Too many bonds for atom',I5,' Check code')
    !
    !
    NBONDP=NBOND+1
    DO IBOND=NBONDP,NBONDT
       IF(IB(IBOND).GT.JB(IBOND)) THEN
          I=IB(IBOND)
          IB(IBOND)=JB(IBOND)
          JB(IBOND)=I
       ENDIF
    ENDDO
    !
    DO ITRANS=1,NTRANS
       IF (ITRANS.EQ.1) THEN
          IS=NATOM+1
       ELSE
          IS=IMATPT(ITRANS-1)+1
       ENDIF
       IQ=IMATPT(ITRANS)
       LFOUND=.FALSE.
       DO IBOND=NBONDP,NBONDT
          IF(JB(IBOND).LE.IQ .AND. JB(IBOND).GE.IS) THEN
             !
             !           THIS BOND CONNECT ITRANS. PROCESS IT.
             !
             IF (.NOT.(LFOUND)) THEN
                LFOUND=.TRUE.
                DO J=1,NATOM
                   IATTR(J)=0
                ENDDO
                DO I=IS,IQ
                   J=IMATTR(I)
                   IATTR(J)=I
                ENDDO
             ENDIF
             !
             !           SET UP LINKAGE ARRAYS
             !
             IPR=IB(IBOND)
             IIM=IMATTR(JB(IBOND))
             NPRA=NATBON(IPR)
             DO J=1,NPRA
                IJ=IATBON(J,IPR)
                IF (IJ.GT.0) THEN
                   IPRA(J)=JB(IJ)
                ELSE
                   IPRA(J)=IB(-IJ)
                ENDIF
             ENDDO
             !
             NIMA=NATBON(IIM)
             DO J=1,NIMA
                IJ=IATBON(J,IIM)
                IF (IJ.GT.0) THEN
                   IIMA(J)=JB(IJ)
                ELSE
                   IIMA(J)=IB(-IJ)
                ENDIF
             ENDDO
             !
             NPRB=0
             DO I=1,NPRA
                IA=IPRA(I)
                NA=NATBON(IA)
                DO J=1,NA
                   IJ=IATBON(J,IA)
                   IF (IJ.GT.0) THEN
                      JA=JB(IJ)
                   ELSE
                      JA=IB(-IJ)
                   ENDIF
                   IF(JA.NE.IPR) THEN
                      NPRB=NPRB+1
                      IPRB(NPRB)=JA
                   ENDIF
                ENDDO
             ENDDO
             !
             NIMB=0
             DO I=1,NIMA
                IA=IIMA(I)
                NA=NATBON(IA)
                DO J=1,NA
                   IJ=IATBON(J,IA)
                   IF (IJ.GT.0) THEN
                      JA=JB(IJ)
                   ELSE
                      JA=IB(-IJ)
                   ENDIF
                   IF(JA.NE.IIM) THEN
                      NIMB=NIMB+1
                      IIMB(NIMB)=JA
                   ENDIF
                ENDDO
             ENDDO
             !
             !
             IF (MODEX.GT.1) THEN
                I2=IPR
                J2=IATTR(IIM)
                IF(I2.GT.0 .AND. J2.GT.NATOM) THEN
                   NPAIR=NPAIR+1
                   IF(NPAIR.GT.MAXWRK) THEN
                      !                 RAN-OUT-OF-SPACE
                      IF(WRNLEV.GE.2) WRITE(OUTU,988)
                      RETURN
                   ENDIF
                   IPK(NPAIR)=J2
                   JPK(NPAIR)=I2
                ENDIF
             ENDIF
             !
             IF (MODEX.GT.2) THEN
                I2=IPR
                DO I=1,NIMA
                   J2=IATTR(IIMA(I))
                   IF(I2.GT.0 .AND. J2.GT.NATOM) THEN
                      NPAIR=NPAIR+1
                      IF(NPAIR.GT.MAXWRK) THEN
                         !                   RAN-OUT-OF-SPACE
                         IF(WRNLEV.GE.2) WRITE(OUTU,988)
                         RETURN
                      ENDIF
                      IPK(NPAIR)=J2
                      JPK(NPAIR)=I2
                   ENDIF
                ENDDO
                J2=IATTR(IIM)
                DO I=1,NPRA
                   I2=IPRA(I)
                   IF(I2.GT.0 .AND. J2.GT.NATOM) THEN
                      NPAIR=NPAIR+1
                      IF(NPAIR.GT.MAXWRK) THEN
                         !                   RAN-OUT-OF-SPACE
                         IF(WRNLEV.GE.2) WRITE(OUTU,988)
                         RETURN
                      ENDIF
                      IPK(NPAIR)=J2
                      JPK(NPAIR)=I2
                   ENDIF
                ENDDO
             ENDIF
             !
             !
             IF (MODEX.GT.3) THEN
                I2=IPR
                DO I=1,NIMB
                   J2=IATTR(IIMB(I))
                   IF(I2.GT.0 .AND. J2.GT.NATOM) THEN
                      NPAIR=NPAIR+1
                      IF(NPAIR.GT.MAXWRK) THEN
                         !                   RAN-OUT-OF-SPACE
                         IF(WRNLEV.GE.2) WRITE(OUTU,988)
                         RETURN
                      ENDIF
                      IPK(NPAIR)=J2
                      JPK(NPAIR)=I2
                      IF(MODEX.EQ.5) THEN
                         NPAIR4=NPAIR4+1
                         IF(NPAIR4.GT.MXWRK4) THEN
                            !                     RAN-OUT-OF-SPACE
                            IF(WRNLEV.GE.2) WRITE(OUTU,988)
                            RETURN
                         ENDIF
                         IPK14(NPAIR4)=J2
                         JPK14(NPAIR4)=I2
                      ENDIF
                   ENDIF
                ENDDO
                J2=IATTR(IIM)
                DO I=1,NPRB
                   I2=IPRB(I)
                   IF(I2.GT.0 .AND. J2.GT.NATOM) THEN
                      NPAIR=NPAIR+1
                      IF(NPAIR.GT.MAXWRK) THEN
                         !                   RAN-OUT-OF-SPACE
                         IF(WRNLEV.GE.2) WRITE(OUTU,988)
                         RETURN
                      ENDIF
                      IPK(NPAIR)=J2
                      JPK(NPAIR)=I2
                      IF(MODEX.EQ.5) THEN
                         NPAIR4=NPAIR4+1
                         IF(NPAIR4.GT.MXWRK4) THEN
                            !                     RAN-OUT-OF-SPACE
                            IF(WRNLEV.GE.2) WRITE(OUTU,988)
                            RETURN
                         ENDIF
                         IPK14(NPAIR4)=J2
                         JPK14(NPAIR4)=I2
                      ENDIF
                   ENDIF
                ENDDO
                DO I=1,NPRA
                   I2=IPRA(I)
                   DO J=1,NIMA
                      J2=IATTR(IIMA(J))
                      IF(I2.GT.0 .AND. J2.GT.NATOM) THEN
                         NPAIR=NPAIR+1
                         IF(NPAIR.GT.MAXWRK) THEN
                            !                     RAN-OUT-OF-SPACE
                            IF(WRNLEV.GE.2) WRITE(OUTU,988)
                            RETURN
                         ENDIF
                         IPK(NPAIR)=J2
                         JPK(NPAIR)=I2
                         IF(MODEX.EQ.5) THEN
                            NPAIR4=NPAIR4+1
                            IF(NPAIR4.GT.MXWRK4) THEN
                               !                       RAN-OUT-OF-SPACE
                               IF(WRNLEV.GE.2) WRITE(OUTU,988)
                               RETURN
                            ENDIF
                            IPK14(NPAIR4)=J2
                            JPK14(NPAIR4)=I2
                         ENDIF
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
             !
          ENDIF
       ENDDO
    ENDDO
    !
    !     Next sort the list of all the possible 1-4 interactions.
    !
    IF(MODEX.EQ.5) THEN
       NPAIR4=NPAIR4+1
       IF(NPAIR4.GT.MXWRK4) THEN
          !         RAN-OUT-OF-SPACE
          IF(WRNLEV.GE.2) WRITE(OUTU,988)
          RETURN
       ENDIF
       IPK14(NPAIR4)=NATIM
       JPK14(NPAIR4)=NATIM
       CALL SORT(NPAIR4,EXCH5,ORDER5,IPK14,JPK14,0,0,0,0,0,2)
    ENDIF
    !
    !     SORT THE PAIR LIST.
    !
    CALL SORT(NPAIR,EXCH5,ORDER5,IPK,JPK,0,0,0,0,0,2)
    !
    !     PROCESS THE SORTED PAIR LIST TO MAKE INB. CHECK THAT THERE ARE NOT
    !     MULTIPLE ENTRIES.
    !
    NNB14=0
    NX14=0
    I2L=0
    J2L=0
    IATOM=1
    NEXT14=1
    LEX14=.FALSE.
    DO I=1,NPAIR
       I2=IPK(I)
       J2=JPK(I)
       IF(MODEX.EQ.5) THEN
          LEX14=I2.EQ.IPK14(NEXT14) .AND. J2.EQ.JPK14(NEXT14)
          IF(LEX14) NEXT14=NEXT14+1
       ENDIF
       IF (I2.EQ.I2L .AND. J2.EQ.J2L) THEN
          !         THIS IS A REPITION. REMOVE 1-4 SIGN IF APPLIED.
          IF(.NOT.LEX14 .AND. INB14(NNB14).LT.0) THEN
             INB14(NNB14)=-INB14(NNB14)
             NX14=NX14-1
          ENDIF
       ELSE
          !         THIS IS A NEW ONE.
          I2L=I2
          J2L=J2
          IF (I2.GT.IATOM) THEN
             DO J=IATOM,I2-1
                IBLO14(J)=NNB14
             ENDDO
             IATOM=I2
          ENDIF
          NNB14=NNB14+1
          !
          IF (NNB14.GT.MXNB14) THEN
             !           RAN-OUT-OF-SPACE
             IF(WRNLEV.GE.2) WRITE(OUTU,988)
             RETURN
          ENDIF
          !
          IF (LEX14) THEN
             INB14(NNB14)=-J2
             NX14=NX14+1
          ELSE
             INB14(NNB14)=J2
          ENDIF
       ENDIF
    ENDDO
    !
    DO J=IATOM,NATIM
       IBLO14(J)=NNB14
    ENDDO
    CMPLTD=.TRUE.
    I=NNB14-NX14
    !
    IF(PRNLEV.GE.5) WRITE(OUTU,45) MODE,I,NX14
45  FORMAT(' <MKIMNB> with mode',I4,' found',I7,' exclusions and',I7, &
         ' interactions(1-4)')
988 FORMAT(' <MKIMNB>: Ran out of space. RESIZING')
    !
    RETURN
  END SUBROUTINE MKIMNB

end module upimag_util
