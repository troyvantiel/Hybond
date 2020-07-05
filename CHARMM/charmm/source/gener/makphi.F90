SUBROUTINE MAKPHI
  use chm_kinds
  use dimens_fcm
  use number
  use consta
  use psf
  use param
  use stream
  use vangle_mm
  use code
  use memory
  implicit none
  ! local
  INTEGER IPHI,IC
  LOGICAL GO,ERR
  !
  !  Generate arrays for fast vector ephi routines.
  !  Author: Stephen Fleischman
  !
  !  Trigonometric tables introduced to support
  !   the new dihedral energy routines.
  !
  !   by: Arnaud Blondel
  !
  !     determine memory requirement
  NPHIV=0
  DO IPHI = 1,NPHI
     IC=ICP(IPHI)
     IF (IC.NE.0) THEN
        NPHIV = NPHIV+1
        GO =.TRUE.
        DO WHILE(GO)
           IF (CPD(IC).LT.0) THEN
              IC=IC+1
              NPHIV = NPHIV+1
           ELSE
              GO = .FALSE.
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  !     allocate memory for the vector phi arrays
  IF(ALLOCATED(IPV)) THEN
     CALL CHMREALLOC('makphi.src','MAKPHI','IPV',newsz=NPHIV,intg=IPV)
  ELSE
     CALL CHMALLOC('makphi.src','MAKPHI','IPV',size=NPHIV,intg=IPV)
  ENDIF
  IF(ALLOCATED(JPV)) THEN
     CALL CHMREALLOC('makphi.src','MAKPHI','JPV',newsz=NPHIV,intg=JPV)
  ELSE
     CALL CHMALLOC('makphi.src','MAKPHI','JPV',size=NPHIV,intg=JPV)
  ENDIF
  IF(ALLOCATED(KPV)) THEN
     CALL CHMREALLOC('makphi.src','MAKPHI','KPV',newsz=NPHIV,intg=KPV)
  ELSE
     CALL CHMALLOC('makphi.src','MAKPHI','KPV',size=NPHIV,intg=KPV)
  ENDIF
  IF(ALLOCATED(LPV)) THEN
     CALL CHMREALLOC('makphi.src','MAKPHI','LPV',newsz=NPHIV,intg=LPV)
  ELSE
     CALL CHMALLOC('makphi.src','MAKPHI','LPV',size=NPHIV,intg=LPV)
  ENDIF
  IF(ALLOCATED(VIND)) THEN
     CALL CHMREALLOC('makphi.src','MAKPHI','VIND',newsz=NPHIV,intg=VIND)
  ELSE
     CALL CHMALLOC('makphi.src','MAKPHI','VIND',size=NPHIV,intg=VIND)
  ENDIF
  IF (associated(VCPC)) THEN
     call chmrealloc('makphi.src','MAKPHI','VCPC',NPHIV,crlp=VCPC)
  ELSE
     call chmalloc('makphi.src','MAKPHI','VCPC',NPHIV,crlp=VCPC)
  ENDIF
  IF (associated(VCPD)) THEN
     call chmrealloc('makphi.src','MAKPHI','VCPD',NPHIV,intgp=VCPD)
  ELSE
     call chmalloc('makphi.src','MAKPHI','VCPD',NPHIV,intgp=VCPD)
  ENDIF
  IF (associated(VCPB)) THEN
     CALL chmrealloc('makphi.src','MAKPHI','VCPB',NPHIV,crlp=VCPB)
  ELSE
     CALL chmalloc('makphi.src','MAKPHI','VCPB',NPHIV,crlp=VCPB)
  ENDIF
  ! Trigonometric tables. AB.
  IF (associated(VCPCOS)) THEN
     CALL chmrealloc('makphi.src','MAKPHI','VCPCOS',NPHIV,crlp=VCPCOS)
  ELSE
     CALL chmalloc('makphi.src','MAKPHI','VCPCOS',NPHIV,crlp=VCPCOS)
  ENDIF
  IF (associated(VCPSIN)) THEN
     CALL chmrealloc('makphi.src','MAKPHI','VCPSIN',NPHIV,crlp=VCPSIN)
  ELSE
     CALL chmalloc('makphi.src','MAKPHI','VCPSIN',NPHIV,crlp=VCPSIN)
  ENDIF
  ! AB.
  !
  CALL MAKPHI2(NPHIV,IPV,JPV,KPV,LPV,VIND, &
       VCPC,VCPD,VCPB, &
       VCPCOS,VCPSIN,ERR)
  RETURN
END SUBROUTINE MAKPHI
!
SUBROUTINE MAKPHI2(NPHIV,IPV,JPV,KPV,LPV,VIND,VCPC,VCPD,VCPB, &
     VCPCOS,VCPSIN,ERR)
  use chm_kinds
  use dimens_fcm
  use number
  use consta
  use psf
  use param
  use stream
  use code
  implicit none
  !
  INTEGER IPV(*),JPV(*),KPV(*),LPV(*),VIND(*),VCPD(*),NPHIV
  real(chm_real) VCPC(*),VCPB(*),VCPCOS(*),VCPSIN(*)
  LOGICAL ERR
  !
  !     indices
  INTEGER IC,IPHI
  !     difference vectors
  LOGICAL GO
  !
  ERR = .FALSE.
  IF (NPHI.EQ.0) RETURN
  !
  NPHIV=0
  DO IPHI = 1,NPHI
     IC=ICP(IPHI)
     IF (IC.NE.0) THEN
        NPHIV=NPHIV+1
        ! Limitation to 6 kept for the parallel routines. AB.
        IF (ABS(CPD(IC)).GT.6 .OR. CPD(IC).EQ.0) THEN
           ! AB.
           CALL WRNDIE(-5,'<MAKPHI>', &
                'Bad periodicity in list for dihedral angles.')
           ERR = .TRUE.
           RETURN
        ENDIF
        GO =.TRUE.
        DO WHILE(GO)
           IPV(NPHIV) = IP(IPHI)
           JPV(NPHIV) = JP(IPHI)
           KPV(NPHIV) = KP(IPHI)
           LPV(NPHIV) = LP(IPHI)
           VIND(NPHIV)=IPHI
           VCPC(NPHIV)=CPC(IC)
           VCPD(NPHIV)=CPD(IC)
           VCPB(NPHIV)=CPB(IC)
           VCPCOS(NPHIV)=COS(VCPB(NPHIV))
           VCPSIN(NPHIV)=SIN(VCPB(NPHIV))
           IF (VCPD(NPHIV).LT.0) THEN
              VCPD(NPHIV)=-VCPD(NPHIV)
              IC=IC+1
              NPHIV=NPHIV+1
           ELSE
              GO = .FALSE.
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  !
  IF(PRNLEV.GT.7) WRITE(OUTU,845) NPHI,NPHIV
845 FORMAT(' MAKPHI:: From',I6,' dihedrals, a total of',I6, &
       ' fast primitive dihedral elements were produced.')
  !
  RETURN
END SUBROUTINE MAKPHI2

SUBROUTINE AUTOGEN(COMLYN,COMLEN)
  !
  ! This routine embodied the AUTOgenerate command which
  ! conditionally deletes all angles and/or dihedrals
  ! and regenerates them for the entire PSF.
  ! This command is intended for the conclusion of a series
  ! of patch commands.
  !                       Bernard R. Brooks, NIH, 18-JUL-1995
  !
  use chm_kinds
  use stream
  use string
  use dimens_fcm
  use psf
  use genpsf_m, only: autgen
  use rtf,only:autop,autod,autot
  use memory
  implicit none

  integer,allocatable,dimension(:,:) :: IATBON
  integer,allocatable,dimension(:) :: NATBON
  CHARACTER(len=*) COMLYN
  INTEGER COMLEN
  LOGICAL LANGLE,LPHI,LPATCH,LOFF,LDRUDE
  !
  LANGLE=(INDXA(COMLYN, COMLEN, 'ANGL') .GT. 0)
  LPHI  =(INDXA(COMLYN, COMLEN, 'DIHE') .GT. 0)
  LPATCH=(INDXA(COMLYN, COMLEN, 'PATC') .GT. 0)
  LOFF  =(INDXA(COMLYN, COMLEN, 'OFF') .GT. 0)
  LDRUDE=(INDXA(COMLYN, COMLEN, 'DRUD') .GT. 0)
  !
  IF(LANGLE) THEN
     NTHETA=0
     NTHETT=0
     AUTOT=.TRUE.
     !mf 7/14/10 only print if prnlev is 2 or larger
     IF (PRNLEV.GE.2) WRITE(OUTU,45) 'angles'
  ENDIF
  IF(LPHI) THEN
     NPHI=0
     NPHIT=0
     AUTOD=.TRUE.
     !mf 7/14/10 only print if prnlev is 2 or larger
     IF (PRNLEV.GE.2) WRITE(OUTU,45) 'dihedrals'
  ENDIF
45 FORMAT(' AUTOGEN: All ',A,' are removed and regenerated.')
  IF(LPATCH) THEN
     AUTOP=.TRUE.
  ENDIF
  IF(LDRUDE) THEN
     QDRUDE=.TRUE.
  ENDIF
  IF(LANGLE .OR. LPHI .OR. LDRUDE) THEN
     call chmalloc('makphi.src','AUTOGEN','NATBON',NATOM,intg=NATBON)
     call chmalloc('makphi.src','AUTOGEN','IATBON',IATBMX,NATOM,intg=IATBON)

     !yw... Jul-07-2007 fix eh070703 add .FLASE. for LDRUDE
     CALL AUTGEN(1,NATBON,IATBON,LANGLE,LPHI,QDRUDE)
     !yw...
     call chmdealloc('makphi.src','AUTOGEN','NATBON',NATOM,intg=NATBON)
     call chmdealloc('makphi.src','AUTOGEN','IATBON',IATBMX,NATOM,intg=IATBON)
     CALL PSFSUM(OUTU)
  ELSE IF (LOFF) THEN
     AUTOT=.FALSE.
     AUTOD=.FALSE.
     AUTOP=.FALSE.
     IF (PRNLEV.GE.2) WRITE(OUTU,*) ' AUTOGEN: turned off.'
  ELSE
     CALL WRNDIE(0,'<AUTOGEN>', &
          'Neither ANGLe nor DIHEdral was specified')
  ENDIF
  !
  RETURN
END SUBROUTINE AUTOGEN

