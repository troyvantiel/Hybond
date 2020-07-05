module sgld
  use chm_kinds
  use dimens_fcm
  implicit none
#if KEY_SGLD==1 /*vars*/

!  Self Guided Langevin Dynamics (SGLD)
!
!        By Xiongwu Wu and Bernard Brooks, NIH, July 2003
!                                       Modified Dec 2005
!                                       Modified July 2010
!                                       Modified Jan 2013
!                                       Modified July 2015
!
!      variables and commom areas for SGLD simulation
!
!     QSGLD - logical flag for SGLD
!     QSGMD - logical flag for SGMDp simulation
!     QSGBZ - logical flag to enforce Boltzmann distribution
!     QSGSTOPT - logical flag to remove translation of the center of mass
!     QSGSTOPR - logical flag to remove rotation of the center of mass
!     QSGNONET - logical flag to remove guiding force and torque
!     QCTSG - logical flag for constant guiding temperature
!     QSGLDG - logical flag for SGLDG to maintain ensemble distribution
!
!    SGLD ARRAYS (HEAP INDEX):
!     SGVX,SGVY,SGVZ   ! velocity averages
!     SGFX,SGFY,SGFZ   ! force averages
!     SGDX,SGDY,SGDZ   ! Interaction forces
!     SGGX,SGGY,SGGZ   ! Storage arrays
!     SGHX,SGHY,SGHZ   ! Storage arrays
!     SGKX,SGKY,SGKZ   ! Storage arrays
!     SHKX,SHKY,SHKZ   ! SHAKE forces
!     SGWT             ! weighting factor for guiding forces. Using SCALAR to manipulate.
!
!    SGLD APPLYING ATOM RANGE
!      ISGSTA     ! Begining atom index applying sgld
!      ISGEND     ! Ending atom index applying sgld
!
!    SGLD VARIABLES
!     SGAVG0  ! Local average remains
!     SGAVG1  ! Local average factor, SGAVG1=1-SGAVG0
!     SGAVP0  ! Convergency control average remains
!     SGAVP1  ! Convergency control average factor, SGAVP1=1-SGAVP0
!     TSGAVG  ! Local average time, ps
!     TSGAVP  ! Convergency control average time, ps
!     SGFT    ! Momentum guiding factor
!     SGFTI   ! Current momentum guiding factor
!     SGFF    ! force guiding factor
!     SGFFI   ! Current force guiding factor
!     TEMPSG  ! Target guiding temperature, K
!     TEMPLF  ! Low frequency  temperature, K
!     TREFLF  ! Reference low frequency temperature, K
!     TEMPLFI ! Current low frequency  temperature, K
!     TREFLFI ! Current reference low frequency  temperature, K
!     TEMPHFI ! Current high frequency  temperature, K
!     TREFHFI ! Current reference high frequency  temperature, K
!     FRCHF   ! High frequency force parameter
!     FRCLF   ! Low frequency force parameter
!     VIRSG   ! Virial from the guiding forces
!
  real(chm_real),allocatable,dimension(:) :: SGWT,SGGAMMA
  real(chm_real),allocatable,dimension(:) :: SGVX,SGVY,SGVZ,SGFX,SGFY,SGFZ,SGDX,SGDY,SGDZ
  real(chm_real),allocatable,dimension(:) :: SGGX,SGGY,SGGZ,SGHX,SGHY,SGHZ,SGKX,SGKY,SGKZ
  real(chm_real),allocatable,dimension(:) :: SHKDX,SHKDY,SHKDZ
  integer ISGSTA,ISGEND,NDEGFSG,NSGSUM,NSGOUT

  real(chm_real) TSGAVG,TSGAVP,SGAVG0,SGAVG1,SGAVP0,SGAVP1
  real(chm_real) SGFT,SGFTI,SGFF,SGFFI,SGFD,SGFDI,SGFDJ,FRCLF,FRCHF,SGSCAL,FSGLDG,PSGLDG
  real(chm_real) TSGSET,TEMPSG,TREFLF,TREFHF,TEMPLF,TEMPHF,TOTMASS
  real(chm_real) TEMPSGI,TREFLFI,TREFHFI,TEMPLFI,TEMPHFI,TRXLF,TRXHF,TEMPLLI
  real(chm_real) EPOTLF,EPOTLF2,EPOTHF,VIRSG,VPRESG(9),EREFLF
  real(chm_real) SGEXPLF,SGEXPHF,SGEXP,DEPDT
  real(chm_real) ACWTSG,ACEREF,ACTEMPLF,ACTEMPHF,ACTREFLF,ACTREFHF,ACFRCLF,ACFRCHF
  real(chm_real) AVGEFLF,AVGEFHF,AVGCFLF,AVGCFHF,AVGTLF,AVGSGEXP
  real(chm_real) SGSUM(20),SGSUM2(20)
  real(chm_real) AVGV2,AVGG2,AVGFF,AVGDD,AVGGF,AVGGD,AVGFD
  real(chm_real) AVGPP,AVGGP,AVGFP,AVGQQ,AVGHQ,AVGDQ
  real(chm_real) SGCFLF,SGCFHF,SGCGLF,SGCGHF,SGEFLF,SGEFHF

  LOGICAL QSGLD,QSGMD,QSGBZ,QSGSTOPT,QSGSTOPR,QSGNONET,QCTSG,QSGLDG
#endif /* (vars)*/

contains
#if KEY_SGLD==1 /*sgld_main*/
  SUBROUTINE PSGLD(NATOM,TIMEST,VX,VY,VZ)
    !
    !    This routine allocate memory for Average arrays
    !
  use chm_kinds
  use cnst_fcm
  use memory
  use stream, only:prnlev,outu

    INTEGER NATOM
    real(chm_real) TIMEST
    real(chm_real) VX(*),VY(*),VZ(*)
    !
    if(prnlev>0)write(OUTU,'("     SGMD-SGLD VERSION: 07-15-2015     ")')
    !  Allocate memory for SGLD arrays
    call chmalloc('sglds.src','PSGLD','SGVX',NATOM,crl=SGVX)
    call chmalloc('sglds.src','PSGLD','SGVY',NATOM,crl=SGVY)
    call chmalloc('sglds.src','PSGLD','SGVZ',NATOM,crl=SGVZ)
    call chmalloc('sglds.src','PSGLD','SGDX',NATOM,crl=SGDX)
    call chmalloc('sglds.src','PSGLD','SGDY',NATOM,crl=SGDY)
    call chmalloc('sglds.src','PSGLD','SGDZ',NATOM,crl=SGDZ)
    call chmalloc('sglds.src','PSGLD','SGFX',NATOM,crl=SGFX)
    call chmalloc('sglds.src','PSGLD','SGFY',NATOM,crl=SGFY)
    call chmalloc('sglds.src','PSGLD','SGFZ',NATOM,crl=SGFZ)
    call chmalloc('sglds.src','PSGLD','SGGX',NATOM,crl=SGGX)
    call chmalloc('sglds.src','PSGLD','SGGY',NATOM,crl=SGGY)
    call chmalloc('sglds.src','PSGLD','SGGZ',NATOM,crl=SGGZ)
    call chmalloc('sglds.src','PSGLD','SGHX',NATOM,crl=SGHX)
    call chmalloc('sglds.src','PSGLD','SGHY',NATOM,crl=SGHY)
    call chmalloc('sglds.src','PSGLD','SGHZ',NATOM,crl=SGHZ)
    call chmalloc('sglds.src','PSGLD','SGKX',NATOM,crl=SGKX)
    call chmalloc('sglds.src','PSGLD','SGKY',NATOM,crl=SGKY)
    call chmalloc('sglds.src','PSGLD','SGKZ',NATOM,crl=SGKZ)
    call chmalloc('sglds.src','PSGLD','SHKDX',NATOM,crl=SHKDX)
    call chmalloc('sglds.src','PSGLD','SHKDY',NATOM,crl=SHKDY)
    call chmalloc('sglds.src','PSGLD','SHKDZ',NATOM,crl=SHKDZ)
    call chmalloc('sglds.src','PSGLD','SGGAMMA',NATOM,crl=SGGAMMA)

    CALL SGINIT(NATOM,ISGSTA,ISGEND,TIMEST,TSGAVG, &
         FBETA,VX,VY,VZ)

    !     Caculate initial constants
    RETURN
  END SUBROUTINE PSGLD

  SUBROUTINE SGFREE
    !
    !      This routinefREE memory space for Average arrays
    !
  use chm_kinds
  use psf
  use memory
    !WXW    free storage for SGLD substructure data
    call chmdealloc('sglds.src','SGFREE','SGVX',NATOM,crl=SGVX)
    call chmdealloc('sglds.src','SGFREE','SGVY',NATOM,crl=SGVY)
    call chmdealloc('sglds.src','SGFREE','SGVZ',NATOM,crl=SGVZ)
    call chmdealloc('sglds.src','SGFREE','SGDX',NATOM,crl=SGDX)
    call chmdealloc('sglds.src','SGFREE','SGDY',NATOM,crl=SGDY)
    call chmdealloc('sglds.src','SGFREE','SGDZ',NATOM,crl=SGDZ)
    call chmdealloc('sglds.src','SGFREE','SGFX',NATOM,crl=SGFX)
    call chmdealloc('sglds.src','SGFREE','SGFY',NATOM,crl=SGFY)
    call chmdealloc('sglds.src','SGFREE','SGFZ',NATOM,crl=SGFZ)
    call chmdealloc('sglds.src','SGFREE','SGGX',NATOM,crl=SGGX)
    call chmdealloc('sglds.src','SGFREE','SGGY',NATOM,crl=SGGY)
    call chmdealloc('sglds.src','SGFREE','SGGZ',NATOM,crl=SGGZ)
    call chmdealloc('sglds.src','SGFREE','SGHX',NATOM,crl=SGHX)
    call chmdealloc('sglds.src','SGFREE','SGHY',NATOM,crl=SGHY)
    call chmdealloc('sglds.src','SGFREE','SGHZ',NATOM,crl=SGHZ)
    call chmdealloc('sglds.src','SGFREE','SGKX',NATOM,crl=SGKX)
    call chmdealloc('sglds.src','SGFREE','SGKY',NATOM,crl=SGKY)
    call chmdealloc('sglds.src','SGFREE','SGKZ',NATOM,crl=SGKZ)
    call chmdealloc('sglds.src','SGFREE','SHKDX',NATOM,crl=SHKDX)
    call chmdealloc('sglds.src','SGFREE','SHKDY',NATOM,crl=SHKDY)
    call chmdealloc('sglds.src','SGFREE','SHKDZ',NATOM,crl=SHKDZ)
    call chmdealloc('sglds.src','SGFREE','SGGAMMA',NATOM,crl=SGGAMMA)
    RETURN
  END SUBROUTINE SGFREE

  SUBROUTINE ALLOCATE_SGWT(NATOM)
    !
    !    This routine allocate memory for Average arrays
    !
  use chm_kinds
  use cnst_fcm
  use memory
  use stream, only:prnlev,outu

    INTEGER NATOM
    !
    !  Allocate memory for SGLD array SGWT
    call chmalloc('sglds.src','ALLOCATE_SGWT','SGWT',NATOM,crl=SGWT)
    SGWT=1.0D0
    RETURN
  END SUBROUTINE ALLOCATE_SGWT

 SUBROUTINE SGINIT(NATOM,ISGSTA,ISGEND,TIMEST,TSGAVG, &
       FBETA,VX,VY,VZ)  !  ,AVGVX,AVGVY,AVGVZ,AVGV2,AVGG2)
    !
    !    This routine perform initiation for average arrays
    !
  use chm_kinds
  use dimens_fcm
  use consta
  use number
  use psf, only: amass
  use stream
  use shake
    INTEGER NATOM,ISGSTA,ISGEND
    INTEGER I,II,NDEGFT
    real(chm_real) VX(*),VY(*),VZ(*),FBETA(*)
!    real(chm_real) SGVX(*),SGVY(*),SGVZ(*),SGV2(*),SGG2(*)
    real(chm_real) VXI,VYI,VZI
    real(chm_real) TIMEST,TSGAVG,SGBETA
    real(chm_real) EGG,GAMM,DELTA,AMASSI,RFD
    ! allocate SGWT
    IF(.not.allocated(sgwt))call allocate_sgwt(natom)
    !WXW    Assign initial values
    NSGOUT=7
    DELTA=TIMEST/TIMFAC
    GAMM=SQRT(TIMEST/TSGAVG)
    IF(TREFLF>RSMALL)GAMM=SQRT(TREFLF/TSGSET)
    NDEGFSG=0
    NDEGFT=0
    TOTMASS=ZERO
    DO I=1,NATOM
       II=0
       IF(I>=ISGSTA .AND. I<=ISGEND)II=1
       IF(QSHAKE)THEN
          NDEGFT=NDEGFT+IDGF2(I)
          NDEGFSG=NDEGFSG+IDGF2(I)*II
       ELSE
          NDEGFT=NDEGFT+3
          NDEGFSG=NDEGFSG+3*II
       ENDIF
       AMASSI=AMASS(I)
       TOTMASS=TOTMASS+AMASSI
       SGGAMMA(I)=ZERO
       !WXW    initialize arrays
       VXI=VX(I)*DELTA
       VYI=VY(I)*DELTA
       VZI=VZ(I)*DELTA
       !SGVX(I)=VXI*GAMM
       !SGVY(I)=VYI*GAMM
       !SGVZ(I)=VZI*GAMM
       SGVX(I)=ZERO
       SGVY(I)=ZERO
       SGVZ(I)=ZERO
       SGFX(I)=ZERO
       SGFY(I)=ZERO
       SGFZ(I)=ZERO
       SGGX(I)=ZERO
       SGGY(I)=ZERO
       SGGZ(I)=ZERO
       SGHX(I)=ZERO
       SGHY(I)=ZERO
       SGHZ(I)=ZERO
       SGKX(I)=ZERO
       SGKY(I)=ZERO
       SGKZ(I)=ZERO
       SHKDX(I)=ZERO
       SHKDY(I)=ZERO
       SHKDZ(I)=ZERO
       AVGV2=AVGV2+AMASSI*(VXI*VXI+VYI*VYI+VZI*VZI)
    ENDDO
    IF(QSHAKE)NDEGFSG=NDEGFSG/2
    AVGG2=AVGV2*TREFLF/TSGSET
    SGFTI=SGFT
    SGFFI=ZERO
    SGFDI=ZERO
    IF(SGFF*SGFF>RSMALL)SGFFI=SGFF
    IF(SGFD*SGFD>RSMALL)SGFDI=SGFD
    IF(TSGSET<RSMALL)TSGSET=MAX(TEMPLF,300.0D0)
    TEMPSGI=TEMPSG
    TREFLFI=TREFLF
    IF(TREFLF<RSMALL)TREFLFI=GAMM*TSGSET*NDEGFSG/NDEGFT
    TEMPLFI=TREFLFI
    TREFHFI=TSGSET-TREFLFI
    IF(TREFHF<RSMALL)TREFHFI=TSGSET-TREFLFI
    TEMPHFI=TREFHFI
    SGEXPLF=ONE
    SGEXPHF=ONE
    IF(TEMPSG<RSMALL)THEN
       TEMPSGI=TSGSET
       TEMPLF=TREFLF
       TEMPHF=TSGSET-TEMPLF
    ELSE
      TEMPLF=TEMPSG*TREFLF*TSGSET/(TSGSET*TEMPSG-(TSGSET-TREFLF)*(TEMPSG-TSGSET))
      TEMPHF=TSGSET-TEMPLF
      TEMPSGI=TEMPSG
      IF(TREFLF>RSMALL)THEN
         SGEXPLF=TREFLFI/TEMPLF
         SGEXPHF=TREFHFI/TEMPHF
      ENDIF
    ENDIF
    IF(QSGLDG)THEN
      IF(TEMPSG<RSMALL)TEMPSG=TSGSET
      TEMPLLI=TEMPLFI
      IF(QSGMD)THEN
        IF(SGFTI*SGFTI>RSMALL)THEN
          SGFFI=SQRT(ONE+SGFTI)-ONE
        ELSE
          SGFTI=(ONE+SGFFI)*(ONE+SGFFI)-ONE
        ENDIF
        FSGLDG=ZERO
        PSGLDG=ZERO
      ELSE
        IF(SGFT>ONE)SGFT=ONE-RSMALL
        FSGLDG=-ONE+SQRT(ONE-SGFTI)
        PSGLDG=(ONE+SGFFI)*(ONE+SGFFI)-ONE
      ENDIF
    ENDIF
    TRXLF=ZERO
    EPOTLF=RBIG*TWO
    SGEXP=ZERO
    VIRSG=ZERO
    VPRESG=ZERO
    FRCLF=ONE
    FRCHF=ONE
    ACWTSG=ZERO
    ACEREF=ZERO
    ACTEMPLF=ZERO
    ACTEMPHF=ZERO
    ACTREFLF=ZERO
    ACTREFHF=ZERO
    ACFRCLF=ZERO
    ACFRCHF=ZERO
    AVGFF=RSMALL
    AVGDD=RSMALL
    AVGFD=RSMALL
    AVGGF=ZERO
    AVGGD=ZERO
    AVGPP=RSMALL
    AVGGP=ZERO
    AVGFP=ZERO
    AVGQQ=RSMALL
    AVGHQ=ZERO
    AVGDQ=ZERO
    SGEFLF=ONE
    SGEFHF=ONE
    SGCFLF=ONE
    SGCFHF=ONE
    SGCGLF=ZERO
    SGCGHF=ZERO
    NSGSUM=0
    DO I=1,2*NSGOUT
      SGSUM(I)=ZERO
      SGSUM2(I)=ZERO
    ENDDO
    IF(QSGSTOPT)THEN
      IF(NDEGFSG.EQ.NDEGFT)NDEGFSG=NDEGFSG-3
      NDEGFT=NDEGFT-3
    ENDIF
    IF(QSGSTOPR)THEN
      IF(NDEGFSG.EQ.NDEGFT)NDEGFSG=NDEGFSG-3
      NDEGFT=NDEGFT-3
    ENDIF
    IF(QSGNONET)THEN
      IF(.NOT.QSGSTOPT)NDEGFSG=NDEGFSG-3
      IF(.NOT.QSGSTOPR)NDEGFSG=NDEGFSG-3
    ENDIF
    IF(PRNLEV>2)WRITE(OUTU,'(" SGLD/SGMD NDEGFSG> ",I6)')NDEGFSG
    RETURN
  END SUBROUTINE SGINIT

  SUBROUTINE SGFAVG(ATFRST,ATLAST,EPOT,VIRINT,VOLUME,IMOVE,AMASS,X,Y,Z,DX,DY,DZ)
    !
    !     This routine calculates the local average forces
    !
  use chm_kinds
  use dimens_fcm
  use consta
  use number
  use reawri
  use econtmod, only:qecont,econt
  use machutil,only:eclock
#if KEY_DOMDEC==1
  use domdec_common,only:q_domdec,natoml,atoml
#endif
#if KEY_PARALLEL==1
  use parallel
  real(chm_real) TIMMER
#endif
    INTEGER ATFRST,ATLAST
    real(chm_real)  EPOT,VIRINT,VOLUME
    INTEGER IMOVE(*)
    real(chm_real)  AMASS(*),X(*),Y(*),Z(*),DX(*),DY(*),DZ(*)
    !
    !
    INTEGER I,IFRST,ILAST,IA
    real(chm_real) EPSUM,DXI,DYI,DZI,TMP,TMPS(20)
    !
    IFRST=ATFRST
    ILAST=ATLAST
    IF(ISGSTA > IFRST) IFRST=ISGSTA
    IF(ISGEND < ILAST) ILAST=ISGEND
    !
    do ia = atfrst,atlast
#if KEY_DOMDEC==1
       if (q_domdec) then
          i = atoml(ia)
       else
#endif
          i = ia
#if KEY_DOMDEC==1
       endif
#endif
       IF(I<ISGSTA .OR. I>ISGEND)CYCLE
       IF(IMOVE(I) == 0) THEN
          DXI=DX(I)
          DYI=DY(I)
          DZI=DZ(I)
          SGDX(I)=DXI
          SGDY(I)=DYI
          SGDZ(I)=DZI
       ENDIF
    ENDDO
    IF(QECONT)THEN
      ! when sgld is applied to only a part of a system
      EPSUM=ZERO
      DO I=IFRST,ILAST
        EPSUM=EPSUM+ECONT(I)
      ENDDO
#if KEY_PARALLEL==1
      CALL PSYNC()
      TIMMER=ECLOCK()
      CALL GCOMB(EPSUM,1)
      TMERI(TIMGCOMM)=TMERI(TIMGCOMM)+ECLOCK()-TIMMER
      TIMMER=ECLOCK()
#endif
    ELSE
      EPSUM=EPOT
    ENDIF
    IF(EPOTLF>RBIG)THEN
       EPOTLF=EPSUM
       EPOTHF=ZERO
       AVGSGEXP=ZERO
    ELSE
       EPOTLF=SGAVG0*EPOTLF+SGAVG1*(EPSUM)
       EPOTHF=EPSUM-EPOTLF
    ENDIF
    ACWTSG=ACWTSG+ONE
    ACEREF=ACEREF+EPOTLF
    ACTEMPLF=ACTEMPLF+TEMPLFI
    IF(TREFHF>RSMALL)THEN
      ACTEMPHF=ACTEMPHF+TEMPHFI
    ELSE
      ACTEMPHF=ACTEMPHF+TSGSET-TEMPLFI
    ENDIF
    ACTREFLF=ACTREFLF+TREFLFI
    ACTREFHF=ACTREFHF+TREFHFI
    ACFRCLF=ACFRCLF+FRCLF
    ACFRCHF=ACFRCHF+FRCHF
    AVGEFLF=SGAVP0*AVGEFLF+SGAVP1*FRCLF
    AVGEFHF=SGAVP0*AVGEFHF+SGAVP1*FRCHF
    AVGCFLF=SGAVP0*AVGCFLF+SGAVP1*TREFLFI/TEMPLFI
    AVGCFHF=SGAVP0*AVGCFHF+SGAVP1*TREFHFI/(TSGSET-TEMPLFI)
    AVGTLF=SGAVP0*AVGTLF+SGAVP1*TEMPLFI
    AVGTLF=SGAVP0*AVGTLF+SGAVP1*TEMPLFI
    IF(QSGLDG)THEN
      !SGEXP=EPOTLF*((ONE-SGFF)*TSGSET/TEMPSG-ONE)/(KBOLTZ*TSGSET)
      !SGEXP=((SGCFLF-ONE)*(EPOTLF-ACEREF/ACWTSG)+  &
      !   (SGCFHF-ONE)*EPOTHF+VIRSG)/(KBOLTZ*TSGSET)
      !SGEXP=-(TEMPSG-TSGSET)*(TSGSET*EPOTLF-AVGTLF*EPOT)/  &
      !   ((TSGSET-AVGTLF)*KBOLTZ*TSGSET*TEMPSG)
      SGEXP=-(TEMPSG-TSGSET)*(TSGSET*(EPOTLF-ACEREF/ACWTSG) &
      -AVGTLF*(EPOT-ACEREF/ACWTSG))/((TSGSET-AVGTLF)*KBOLTZ*TSGSET*TEMPSG)
    ELSE
      SGEXP=((AVGEFLF*AVGCFLF-ONE)*(EPOTLF-ACEREF/ACWTSG)+  &
         (AVGEFHF*AVGCFHF-ONE)*EPOTHF+VIRSG)/(KBOLTZ*TSGSET)
    ENDIF
    AVGSGEXP=SGAVP0*AVGSGEXP+SGAVP1*SGEXP
    !
    ! Guiding variable averages
    NSGSUM=NSGSUM+1
    TMPS(1)=SGFTI
    TMPS(2)=TEMPSGI
    TMPS(3)=TEMPLFI
    TMPS(4)=TREFLFI
    TMPS(5)=FRCLF
    TMPS(6)=EPOTLF
    TMPS(7)=SGEXP
    TMPS(8)=SGFFI
    TMPS(9)=SGSCAL
    TMPS(10)=TEMPHFI
    TMPS(11)=TREFHFI
    TMPS(12)=FRCHF
    TMPS(13)=EPOTHF
    TMPS(14)=VIRSG
    DO I=1,2*NSGOUT
       TMP=TMPS(I)
       SGSUM(I)=SGSUM(I)+TMP
       SGSUM2(I)=SGSUM2(I)+TMP*TMP
    ENDDO
    RETURN
  END SUBROUTINE SGFAVG


  SUBROUTINE SGFSHK(ATFRST,ATLAST,IMOVE,DX,DY,DZ)
    !
    !     This routine calculates the local average forces
    !
  use chm_kinds
  use dimens_fcm
  use consta
  use number
#if KEY_DOMDEC==1
  use domdec_common,only:q_domdec,natoml,atoml
#endif
    INTEGER ATFRST,ATLAST
    INTEGER IMOVE(*)
    real(chm_real)  DX(*),DY(*),DZ(*)
    !
    !
    INTEGER I,IFRST,ILAST,IA
    real(chm_real) DXI,DYI,DZI
    !
    !RETURN
    IFRST=ATFRST
    ILAST=ATLAST
    IF(ISGSTA > IFRST) IFRST=ISGSTA
    IF(ISGEND < ILAST) ILAST=ISGEND
    !
    do ia = atfrst,atlast
#if KEY_DOMDEC==1
       if (q_domdec) then
          i = atoml(ia)
       else
#endif
          i = ia
#if KEY_DOMDEC==1
       endif
#endif
       IF(I<ISGSTA .OR. I>ISGEND)CYCLE
       IF(IMOVE(I) == 0) THEN
          SHKDX(I)=DX(I)
          SHKDY(I)=DY(I)
          SHKDZ(I)=DZ(I)
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE SGFSHK

  subroutine sgld_ave_params()
    use chm_kinds, only: chm_real
    use param_store, only: set_param

    implicit none

    integer i
    real(chm_real) tmps(20)

    do i = 1, 2 * nsgout
       tmps(i) = sgsum(i) / nsgsum
    end do

    call set_param('SGFT', TMPS(1))
    call set_param('TEMPSG', TMPS(2))
    call set_param('TEMPLF', TMPS(3))
    call set_param('TREFLF', TMPS(4))
    call set_param('FRCLF', TMPS(5))
    call set_param('EPOTLF', TMPS(6))
    call set_param('SGEXP', TMPS(7))
    call set_param('SGFF', TMPS(8))
    call set_param('SGFD', TMPS(9))
    call set_param('TEMPHF', TMPS(10))
    call set_param('TREFHF', TMPS(11))
    call set_param('FRCHF', TMPS(12))
    call set_param('EPOTHF', TMPS(13))
    call set_param('VIRSG', TMPS(14))
  end subroutine sgld_ave_params

  SUBROUTINE PRNTSG(OUTU,QHEAD,SGTYPE,QINIT)
    !
    !     This routine calculates the guiding forces and adds them
    !     to the systematic forces.
    !         Assum SGFT*gamma=constant for all atoms
    !
  use chm_kinds
  use dimens_fcm
  use number
  use param_store, only: set_param

  implicit none

    INTEGER OUTU,I
    LOGICAL QHEAD,QINIT
    CHARACTER(len=4) SGTYPE
    CHARACTER(len=10) MARK
    CHARACTER(len=70) LABEL
    real(chm_real) TMP,TMPS(20)
  !WXW    determine data type to be printed
    IF(QSGLD)THEN
       MARK=SGTYPE//' SGLD'
    ELSE IF(QSGMD)THEN
       MARK=SGTYPE//' SGMD'
    ELSE
       RETURN
    ENDIF
    IF(QSGBZ.OR.QSGLDG)SGSCAL=SGFDI
    IF(SGTYPE.EQ.'AVER')THEN
       DO I=1,2*NSGOUT
          TMPS(I)=SGSUM(I)/NSGSUM
       ENDDO
    ELSE IF(SGTYPE.EQ.'FLUC')THEN
       DO I=1,2*NSGOUT
          TMP=SGSUM(I)/NSGSUM
          TMPS(I)=SGSUM2(I)/NSGSUM-TMP*TMP
          IF(TMPS(I)>ZERO)THEN
              TMPS(I)=SQRT(TMPS(I))
          ELSE
              TMPS(I)=ZERO
          ENDIF
       ENDDO
    ELSE IF(SGTYPE.EQ.'DYNA')THEN
       TMPS(1)=SGFTI
       TMPS(2)=TEMPSGI
       TMPS(3)=TEMPLFI
       TMPS(4)=TREFLFI
       TMPS(5)=FRCLF
       TMPS(6)=EPOTLF
       TMPS(7)=SGEXP
       TMPS(8)=SGFFI
       TMPS(9)=SGSCAL
       TMPS(10)=TEMPHFI
       TMPS(11)=TREFHFI
       TMPS(12)=FRCHF
       TMPS(13)=EPOTHF
       TMPS(14)=VIRSG
    ELSE
       TMPS(1)=SGFTI
       TMPS(2)=TEMPSGI
       TMPS(3)=TEMPLFI
       TMPS(4)=TREFLFI
       TMPS(5)=FRCLF
       TMPS(6)=EPOTLF
       TMPS(7)=SGEXP
       TMPS(8)=SGFFI
       TMPS(9)=SGSCAL
       TMPS(10)=TEMPHFI
       TMPS(11)=TREFHFI
       TMPS(12)=FRCHF
       TMPS(13)=EPOTHF
       TMPS(14)=VIRSG
    ENDIF
    IF(QHEAD)THEN
      WRITE(OUTU,'(A4," SGLF: ","    SGFT ","  TEMPSG ",&
&          "   TEMPLF","   TREFLF ", "FRCLF  ", "       EPOTLF   ",  "   SGEXP")')SGTYPE
      WRITE(OUTU,'(A4," SGHF: ","    SGFF ","    SGFD  ",&
&           "  TEMPHF ", "  TREFHF ", "FRCHF ", "        EPOTHF  ", "    VIRSG")')SGTYPE
    ENDIF
    WRITE(OUTU,'(A4," SGLF> ",F8.4,X,F8.2,X,F8.2,X,F8.2,X,F6.4,X,F14.4,X,F10.4)')&
&               SGTYPE,(TMPS(I),I=1,NSGOUT)
    WRITE(OUTU,'(A4," SGHF> ",F8.4,X,F8.4,X,F8.2,X,F8.2,X,F6.4,X,F14.4,X,F10.4)')&
&               SGTYPE,(TMPS(I),I=NSGOUT+1,2*NSGOUT)
    IF(QINIT)THEN
       NSGSUM=0
       DO I=1,2*NSGOUT
          SGSUM(I)=ZERO
          SGSUM2(I)=ZERO
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE PRNTSG

  SUBROUTINE SGLDW(ATFRST,ATLAST,DELTA,NDEGF,TEMPI, &
                   IMOVE,AMASS,X,Y,Z,VX,VY,VZ,DX,DY,DZ)
    !
    !     This routine calculates the guiding forces and adds them
    !     to the systematic forces.
    !         Assum SGFT*gamma=constant for all atoms
    !
    use chm_kinds
    use dimens_fcm
    use number
    use stream
    use consta
    use cnst_fcm
#if KEY_PARALLEL==1
    use parallel
#endif
    use machutil,only:eclock
    use prssre
#if KEY_DOMDEC==1
    use domdec_common,only:q_domdec,natoml,atoml
#endif
    !
    INTEGER ATFRST,ATLAST,NDEGF
    real(chm_real)  DELTA,TEMPI
    INTEGER IMOVE(*)
    real(chm_real)  AMASS(*),X(*),Y(*),Z(*),VX(*),VY(*),VZ(*), &
         DX(*),DY(*),DZ(*)
    !
#if KEY_PARALLEL==1
    real(chm_real) GCARR(31),TIMMER
#endif
    !
    INTEGER I,J,K,IFRST,ILAST,IA
    real(chm_real) XI,YI,ZI,VXI,VYI,VZI,VXT,VYT,VZT
    real(chm_real) FXI,FYI,FZI,GXI,GYI,GZI,DXI,DYI,DZI
    real(chm_real) SGXI,SGYI,SGZI,SGXJ,SGYJ,SGZJ,SGXK,SGYK,SGZK
    real(chm_real) GX2I,GY2I,GZ2I,GX2J,GY2J,GZ2J,GX2K,GY2K,GZ2K
    real(chm_real) GAM,GAMM,FACT0,FACT1,FACT2,DTEMP,SGFTJ
    real(chm_real) AMASSI,DELTA2
    real(chm_real) SUMGV,SUMVV,SUMV2,SUMG2,EKIN,EKINSG,EKINSGG
    real(chm_real) SUMFF,SUMDD,SUMGF,SUMGD
    real(chm_real) SUMPP,SUMGP,SUMFP,SUMQQ,SUMHQ,SUMDQ
    real(chm_real) XCM,YCM,ZCM,XX,XY,XZ,YY,YZ,ZZ
    real(chm_real) FXT,FYT,FZT,PXT,PYT,PZT
    real(chm_real) LXT,LYT,LZT,TXT,TYT,TZT
    real(chm_real) VXCM,VYCM,VZCM,AXCM,AYCM,AZCM
    real(chm_real) VXG,VYG,VZG,OXG,OYG,OZG
    real(chm_real) FXG,FYG,FZG,AXG,AYG,AZG
    real(chm_real) TXG,TYG,TZG,GXCM,GYCM,GZCM
    real(chm_real) OXCM,OYCM,OZCM,QXCM,QYCM,QZCM
    real(chm_real) TCM(3,3),U(3,3),SCR(24),AMOM(6)
    !
    DELTA2=DELTA*DELTA
    !
    IFRST=ATFRST
    ILAST=ATLAST
    IF(ISGSTA > IFRST) IFRST=ISGSTA
    IF(ISGEND < ILAST) ILAST=ISGEND
    !
    !WXW Calculate constraint factor to eliminate energy input from guiding force
    SUMGV=ZERO
    SUMVV=ZERO
    XCM=ZERO
    YCM=ZERO
    ZCM=ZERO
    XX=ZERO
    XY=ZERO
    XZ=ZERO
    YY=ZERO
    YZ=ZERO
    ZZ=ZERO
    PXT=ZERO
    PYT=ZERO
    PZT=ZERO
    LXT=ZERO
    LYT=ZERO
    LZT=ZERO
    FXT=ZERO
    FYT=ZERO
    FZT=ZERO
    TXT=ZERO
    TYT=ZERO
    TZT=ZERO
    FXG=ZERO
    FYG=ZERO
    FZG=ZERO
    TXG=ZERO
    TYG=ZERO
    TZG=ZERO
    do ia = atfrst,atlast
#if KEY_DOMDEC==1
       if (q_domdec) then
          i = atoml(ia)
       else
#endif
          i = ia
#if KEY_DOMDEC==1
       endif
#endif
       IF(IMOVE(I) == 0) THEN
          AMASSI=AMASS(I)
          XI=X(I)
          YI=Y(I)
          ZI=Z(I)
          VXI=VX(I)
          VYI=VY(I)
          VZI=VZ(I)
          DXI=DX(I)
          DYI=DY(I)
          DZI=DZ(I)
          IF(I>=ISGSTA .AND. I<=ISGEND)THEN
            SGFTJ=SGFTI*SGWT(I)
            GAM=TIMFAC*FBETA(I)
            GAMM=TWO/(TWO+GAM*DELTA)
            GXI=SGAVG0*SGVX(I)+SGAVG1*VXI
            GYI=SGAVG0*SGVY(I)+SGAVG1*VYI
            GZI=SGAVG0*SGVZ(I)+SGAVG1*VZI
            SGVX(I)=GXI
            SGVY(I)=GYI
            SGVZ(I)=GZI
            FACT1=SGFTJ*GAM*AMASSI/DELTA
            FXI=SGAVG0*SGFX(I)+SGAVG1*(SGDX(I)+SHKDX(I))
            FYI=SGAVG0*SGFY(I)+SGAVG1*(SGDY(I)+SHKDY(I))
            FZI=SGAVG0*SGFZ(I)+SGAVG1*(SGDZ(I)+SHKDZ(I))
            SGFX(I)=FXI
            SGFY(I)=FYI
            SGFZ(I)=FZI
            IF(QSGLDG)THEN
              SGXJ=SGAVG0*SGGX(I)+SGAVG1*FXI
              SGYJ=SGAVG0*SGGY(I)+SGAVG1*FYI
              SGZJ=SGAVG0*SGGZ(I)+SGAVG1*FZI
              SGGX(I)=SGXJ
              SGGY(I)=SGYJ
              SGGZ(I)=SGZJ
              SGXK=SGAVG0*SGKX(I)+SGAVG1*(DXI-SGDX(I))
              SGYK=SGAVG0*SGKY(I)+SGAVG1*(DYI-SGDY(I))
              SGZK=SGAVG0*SGKZ(I)+SGAVG1*(DZI-SGDZ(I))
              SGKX(I)=SGXK
              SGKY(I)=SGYK
              SGKZ(I)=SGZK
              SGXI=(SGDX(I)+SHKDX(I)-FXI)
              SGYI=(SGDY(I)+SHKDY(I)-FYI)
              SGZI=(SGDZ(I)+SHKDZ(I)-FZI)
              SGHX(I)=SGAVP0*SGHX(I)+SGAVP1*(SGXI*VXI+SGYI*VYI+SGZI*VZI)
              SGHY(I)=SGAVP0*SGHY(I)+SGAVP1*AMASSI*TIMFAC*(VXI*VXI+VYI*VYI+VZI*VZI)
              SGGAMMA(I)=SGAVP0*SGGAMMA(I)+SGAVP1*SGHX(I)/SGHY(I)
              FACT2=PSGLDG*SGGAMMA(I)*TIMFAC*AMASSI*SGWT(I)
              SGXJ=FXI-SGXJ
              SGYJ=FYI-SGYJ
              SGZJ=FZI-SGZJ
              SGXJ=FACT2*GXI-SGFFI*SGXJ
              SGYJ=FACT2*GYI-SGFFI*SGYJ
              SGZJ=FACT2*GZI-SGFFI*SGZJ
              SGXI=FACT1*GXI+SGXJ-FSGLDG*SGXK
              SGYI=FACT1*GYI+SGYJ-FSGLDG*SGYK
              SGZI=FACT1*GZI+SGZJ-FSGLDG*SGZK
            ELSE
              FXI=(SGFFI-SGFDI)*FXI+SGFDI*SGDX(I)
              FYI=(SGFFI-SGFDI)*FYI+SGFDI*SGDY(I)
              FZI=(SGFFI-SGFDI)*FZI+SGFDI*SGDZ(I)
              SGXI=FACT1*GXI-FXI
              SGYI=FACT1*GYI-FYI
              SGZI=FACT1*GZI-FZI
              SGXJ=SGXI
              SGYJ=SGYI
              SGZJ=SGZI
            ENDIF
            DXI=DXI-SGXI
            DYI=DYI-SGYI
            DZI=DZI-SGZI
            FXG=FXG-SGXI
            FYG=FYG-SGYI
            FZG=FZG-SGZI
            TXG=TXG-(YI*SGZI-ZI*SGYI)
            TYG=TYG-(ZI*SGXI-XI*SGZI)
            TZG=TZG-(XI*SGYI-YI*SGXI)
            DX(I)=DXI
            DY(I)=DYI
            DZ(I)=DZI
            FACT2=DELTA2/TWO/AMASSI
            VXT =VX(I)-DXI*FACT2
            VYT =VY(I)-DYI*FACT2
            VZT =VZ(I)-DZI*FACT2
            SUMGV=SUMGV+(SGXJ*VXT+SGYJ*VYT+SGZJ*VZT)*GAMM
            SUMVV=SUMVV+AMASSI*(VXT*VXT+VYT*VYT+VZT*VZT)*GAMM*GAMM
          ENDIF
          XCM=XCM+XI*AMASSI
          YCM=YCM+YI*AMASSI
          ZCM=ZCM+ZI*AMASSI
          XX=XX+XI*XI*AMASSI
          XY=XY+XI*YI*AMASSI
          XZ=XZ+XI*ZI*AMASSI
          YY=YY+YI*YI*AMASSI
          YZ=YZ+YI*ZI*AMASSI
          ZZ=ZZ+ZI*ZI*AMASSI
          PXT=PXT+AMASSI*VXI
          PYT=PYT+AMASSI*VYI
          PZT=PZT+AMASSI*VZI
          LXT=LXT+(YI*VZI-ZI*VYI)*AMASSI
          LYT=LYT+(ZI*VXI-XI*VZI)*AMASSI
          LZT=LZT+(XI*VYI-YI*VXI)*AMASSI
          FXT=FXT+DXI
          FYT=FYT+DYI
          FZT=FZT+DZI
          TXT=TXT+(YI*DZI-ZI*DYI)
          TYT=TYT+(ZI*DXI-XI*DZI)
          TZT=TZT+(XI*DYI-YI*DXI)
       ENDIF
    ENDDO
    !
#if KEY_PARALLEL==1
    CALL PSYNC()
    TIMMER=ECLOCK()
    GCARR(1)=SUMGV
    GCARR(2)=SUMVV
    GCARR(3) = XCM
    GCARR(4) = YCM
    GCARR(5) = ZCM
    GCARR(6) = XX
    GCARR(7) = XY
    GCARR(8) = XZ
    GCARR(9) = YY
    GCARR(10) = YZ
    GCARR(11) = ZZ
    GCARR(12)=PXT
    GCARR(13)=PYT
    GCARR(14)=PZT
    GCARR(15)=LXT
    GCARR(16)=LYT
    GCARR(17)=LZT
    GCARR(18)=FXT
    GCARR(19)=FYT
    GCARR(20)=FZT
    GCARR(21) = TXT
    GCARR(22) = TYT
    GCARR(23) = TZT
    GCARR(24)=FXG
    GCARR(25)=FYG
    GCARR(26)=FZG
    GCARR(27) = TXG
    GCARR(28) = TYG
    GCARR(29) = TZG
    CALL GCOMB(GCARR,29)
    SUMGV=GCARR(1)
    SUMVV=GCARR(2)
    XCM=GCARR(3)
    YCM=GCARR(4)
    ZCM=GCARR(5)
    XX=GCARR(6)
    XY=GCARR(7)
    XZ=GCARR(8)
    YY=GCARR(9)
    YZ=GCARR(10)
    ZZ=GCARR(11)
    PXT=GCARR(12)
    PYT=GCARR(13)
    PZT=GCARR(14)
    LXT=GCARR(15)
    LYT=GCARR(16)
    LZT=GCARR(17)
    FXT=GCARR(18)
    FYT=GCARR(19)
    FZT=GCARR(20)
    TXT=GCARR(21)
    TYT=GCARR(22)
    TZT=GCARR(23)
    FXG=GCARR(24)
    FYG=GCARR(25)
    FZG=GCARR(26)
    TXG=GCARR(27)
    TYG=GCARR(28)
    TZG=GCARR(29)
    TMERI(TIMGCOMM)=TMERI(TIMGCOMM)+ECLOCK()-TIMMER
    TIMMER=ECLOCK()
#endif
    ! rotation inertia
    IF(QSGSTOPR.OR.QSGNONET)THEN
      XCM=XCM/TOTMASS
      YCM=YCM/TOTMASS
      ZCM=ZCM/TOTMASS
      LXT=LXT-(YCM*PZT-ZCM*PYT)
      LYT=LYT-(ZCM*PXT-XCM*PZT)
      LZT=LZT-(XCM*PYT-YCM*PXT)
      TXT=TXT-(YCM*FZT-ZCM*FYT)
      TYT=TYT-(ZCM*FXT-XCM*FZT)
      TZT=TZT-(XCM*FYT-YCM*FXT)
      XX=XX-XCM*XCM*TOTMASS
      XY=XY-XCM*YCM*TOTMASS
      XZ=XZ-XCM*ZCM*TOTMASS
      YY=YY-YCM*YCM*TOTMASS
      YZ=YZ-YCM*ZCM*TOTMASS
      ZZ=ZZ-ZCM*ZCM*TOTMASS
    !
      AMOM(1)=YY+ZZ
      AMOM(2)=-XY
      AMOM(3)=-XZ
      AMOM(4)=XX+ZZ
      AMOM(5)=-YZ
      AMOM(6)=XX+YY
    !
      CALL DIAGQ(3,3,AMOM,U,SCR(4),SCR(7),SCR(10),SCR(13),SCR(1),  &
         SCR(16),SCR(19),SCR(22),0)
    !
      DO I=1,3
         IF(ABS(SCR(I)) < TENM5) SCR(I)=SIGN(TENM5,SCR(I))
         SCR(I)=ONE/SCR(I)
      ENDDO
    !
      DO I=1,3
         DO J=1,3
            TCM(I,J)=0.0
            DO K=1,3
               TCM(I,J)=TCM(I,J)+U(I,K)*U(J,K)*SCR(K)
            ENDDO
         ENDDO
      ENDDO
    ENDIF
    ! net translational and rotational guiding forces
    IF(QSGNONET)THEN
      !  net accelation
      AXG=FXG/TOTMASS
      AYG=FYG/TOTMASS
      AZG=FZG/TOTMASS
      TXG=TXG-(YCM*FZG-ZCM*FYG)
      TYG=TYG-(ZCM*FXG-XCM*FZG)
      TZG=TZG-(XCM*FYG-YCM*FXG)
    ! angular acceleration of COM
      GXCM=TXG*TCM(1,1)+TYG*TCM(1,2)+TZG*TCM(1,3)
      GYCM=TXG*TCM(2,1)+TYG*TCM(2,2)+TZG*TCM(2,3)
      GZCM=TXG*TCM(3,1)+TYG*TCM(3,2)+TZG*TCM(3,3)
      !  net velocity for local average
      VXG=PXT/TOTMASS
      VYG=PYT/TOTMASS
      VZG=PZT/TOTMASS
    ! angular velocity of COM for local average
      OXG=LXT*TCM(1,1)+LYT*TCM(1,2)+LZT*TCM(1,3)
      OYG=LXT*TCM(2,1)+LYT*TCM(2,2)+LZT*TCM(2,3)
      OZG=LXT*TCM(3,1)+LYT*TCM(3,2)+LZT*TCM(3,3)
    ELSE
      AXG=ZERO
      AYG=ZERO
      AZG=ZERO
      GXCM=ZERO
      GYCM=ZERO
      GZCM=ZERO
      !  net velocity for local average
      VXG=ZERO
      VYG=ZERO
      VZG=ZERO
    ! angular velocity of COM for local average
      OXG=ZERO
      OYG=ZERO
      OZG=ZERO
    ENDIF
    ! Translation
    IF(QSGSTOPT)THEN
      !  net velocity
      VXCM=PXT/TOTMASS
      VYCM=PYT/TOTMASS
      VZCM=PZT/TOTMASS
      ! net accelation
      AXCM=FXT/TOTMASS
      AYCM=FYT/TOTMASS
      AZCM=FZT/TOTMASS
      !  net velocity for local average
      VXG=ZERO
      VYG=ZERO
      VZG=ZERO
    ELSE
      VXCM=ZERO
      VYCM=ZERO
      VZCM=ZERO
      AXCM=AXG
      AYCM=AYG
      AZCM=AZG
    ENDIF
    ! Rotation
    IF(QSGSTOPR)THEN
    ! angular velocity of COM
      OXCM=LXT*TCM(1,1)+LYT*TCM(1,2)+LZT*TCM(1,3)
      OYCM=LXT*TCM(2,1)+LYT*TCM(2,2)+LZT*TCM(2,3)
      OZCM=LXT*TCM(3,1)+LYT*TCM(3,2)+LZT*TCM(3,3)
    ! angular acceleration of COM
      QXCM=TXT*TCM(1,1)+TYT*TCM(1,2)+TZT*TCM(1,3)
      QYCM=TXT*TCM(2,1)+TYT*TCM(2,2)+TZT*TCM(2,3)
      QZCM=TXT*TCM(3,1)+TYT*TCM(3,2)+TZT*TCM(3,3)
    ! angular velocity of COM for local average
      OXG=ZERO
      OYG=ZERO
      OZG=ZERO
    ELSE
      OXCM=ZERO
      OYCM=ZERO
      OZCM=ZERO
      QXCM=GXCM
      QYCM=GYCM
      QZCM=GZCM
    ENDIF
    !
    SGSCAL=DELTA2*SUMGV/(SUMVV-HALF*DELTA2*SUMGV)
    !IF(QSGLDG .and. QSGBZ)SGSCAL=ZERO
    !
    EKIN=ZERO
    EKINSG=ZERO
    EKINSGG=ZERO
    SUMFF=ZERO
    SUMDD=ZERO
    SUMGF=ZERO
    SUMGD=ZERO
    SUMPP=ZERO
    SUMGP=ZERO
    do ia = atfrst,atlast
#if KEY_DOMDEC==1
       if (q_domdec) then
          i = atoml(ia)
       else
#endif
          i = ia
#if KEY_DOMDEC==1
       endif
#endif
       IF(IMOVE(I) == 0) THEN
          AMASSI=AMASS(I)
          XI=X(I)-XCM
          YI=Y(I)-YCM
          ZI=Z(I)-ZCM
          VXI=VX(I)-VXCM-OYCM*ZI+OZCM*YI
          VYI=VY(I)-VYCM-OZCM*XI+OXCM*ZI
          VZI=VZ(I)-VZCM-OXCM*YI+OYCM*XI
          VX(I)=VXI
          VY(I)=VYI
          VZ(I)=VZI
          DXI=DX(I)-AMASSI*(AXCM+QYCM*ZI-QZCM*YI)
          DYI=DY(I)-AMASSI*(AYCM+QZCM*XI-QXCM*ZI)
          DZI=DZ(I)-AMASSI*(AZCM+QXCM*YI-QYCM*XI)
          IF(I<ISGSTA .OR. I>ISGEND)THEN
            DX(I)=DXI
            DY(I)=DYI
            DZ(I)=DZI
            CYCLE
          ENDIF
          GAM=TIMFAC*FBETA(I)*DELTA
          GAMM=ONE+HALF*(GAM+SGSCAL)
          SGFTJ=SGFTI*SGWT(I)
          FACT1=(ONE+HALF*GAM)/GAMM
          FACT2=AMASSI*SGSCAL/GAMM/DELTA2
          DXI=FACT1*DXI+FACT2*VXI
          DYI=FACT1*DYI+FACT2*VYI
          DZI=FACT1*DZI+FACT2*VZI
          DX(I)=DXI
          DY(I)=DYI
          DZ(I)=DZI
          FACT2=HALF*DELTA2/GAMM/AMASSI
          VXT =FACT1*VXI-DXI*FACT2
          VYT =FACT1*VYI-DYI*FACT2
          VZT =FACT1*VZI-DZI*FACT2
          ! local average of velocities
          VXI=VXI-VXG-OYG*ZI+OZG*YI
          VYI=VYI-VYG-OZG*XI+OXG*ZI
          VZI=VZI-VZG-OXG*YI+OYG*XI
          GXI=SGVX(I)
          GYI=SGVY(I)
          GZI=SGVZ(I)
          EKIN=EKIN+AMASSI*(VXI*VXI+VYI*VYI+VZI*VZI)
          EKINSG=EKINSG+AMASSI*(GXI*GXI+GYI*GYI+GZI*GZI)
         IF(.NOT.QSGLDG)THEN
          ! Accumulation for SGLD parameter calculation
          FXI=SGFX(I)
          FYI=SGFY(I)
          FZI=SGFZ(I)
          DXI=SGDX(I)
          DYI=SGDY(I)
          DZI=SGDZ(I)
          FACT1=AMASSI/DELTA2
          SGXI=FACT1*((SGSCAL)*VXT-SGFTJ*GAM*GXI)
          SGYI=FACT1*((SGSCAL)*VYT-SGFTJ*GAM*GYI)
          SGZI=FACT1*((SGSCAL)*VZT-SGFTJ*GAM*GZI)
          SGXK=SGXI+SGFDI*DXI+(SGFFI-SGFDI)*FXI
          SGYK=SGYI+SGFDI*DYI+(SGFFI-SGFDI)*FYI
          SGZK=SGZI+SGFDI*DZI+(SGFFI-SGFDI)*FZI
          SGXJ=SGXK-SGFDI*(DXI-FXI)
          SGYJ=SGYK-SGFDI*(DYI-FYI)
          SGZJ=SGZK-SGFDI*(DZI-FZI)
          GX2I=SGAVG0*SGGX(I)+SGAVG1*SGXI
          GY2I=SGAVG0*SGGY(I)+SGAVG1*SGYI
          GZ2I=SGAVG0*SGGZ(I)+SGAVG1*SGZI
          SGGX(I)=GX2I
          SGGY(I)=GY2I
          SGGZ(I)=GZ2I
          GX2K=SGAVG0*SGHX(I)+SGAVG1*SGXK
          GY2K=SGAVG0*SGHY(I)+SGAVG1*SGYK
          GZ2K=SGAVG0*SGHZ(I)+SGAVG1*SGZK
          SGHX(I)=GX2K
          SGHY(I)=GY2K
          SGHZ(I)=GZ2K
          GX2J=SGAVG0*SGKX(I)+SGAVG1*FXI
          GY2J=SGAVG0*SGKY(I)+SGAVG1*FYI
          GZ2J=SGAVG0*SGKZ(I)+SGAVG1*FZI
          SGKX(I)=GX2J
          SGKY(I)=GY2J
          SGKZ(I)=GZ2J
          VXI=VXT-GXI
          VYI=VYT-GYI
          VZI=VZT-GZI
          DXI=DXI-FXI
          DYI=DYI-FYI
          DZI=DZI-FZI
          SUMFF=SUMFF+(FXI*FXI+FYI*FYI+FZI*FZI)
          SUMDD=SUMDD+(DXI*DXI+DYI*DYI+DZI*DZI)
          SUMGF=SUMGF+(FXI*GX2I+FYI*GY2I+FZI*GZ2I)
          SUMGD=SUMGD+(DXI*(SGXI-GX2I)+DYI*(SGYI-GY2I)+DZI*(SGZI-GZ2I))
          FACT1=GAM*AMASSI/DELTA2
          FACT2=FACT1*FACT1
          SUMPP=SUMPP+FACT2*(GXI*GXI+GYI*GYI+GZI*GZI)
          SUMGP=SUMGP+FACT1*(GX2K*GXI+GY2K*GYI+GZ2K*GZI)
         ENDIF
       ENDIF
    ENDDO
#if KEY_PARALLEL==1
    CALL PSYNC()
    TIMMER=ECLOCK()
    GCARR(1)=SUMFF
    GCARR(2)=SUMDD
    GCARR(3)=SUMGF
    GCARR(4)=SUMGD
    GCARR(5)=SUMPP
    GCARR(6)=SUMGP
    GCARR(7)=VIRSG
    GCARR(8)=EKIN
    GCARR(9)=EKINSG
    CALL GCOMB(GCARR,9)
    SUMFF=GCARR(1)
    SUMDD=GCARR(2)
    SUMGF=GCARR(3)
    SUMGD=GCARR(4)
    SUMPP=GCARR(5)
    SUMGP=GCARR(6)
    VIRSG=GCARR(7)
    EKIN=GCARR(8)
    EKINSG=GCARR(9)
    TMERI(TIMGCOMM)=TMERI(TIMGCOMM)+ECLOCK()-TIMMER
    TIMMER=ECLOCK()
#endif
    SGSCAL=SGSCAL/TIMFAC/DELTA
    ! Estimate momentum guiding factor and guiding temperature
    TEMPLFI=SGAVP0*TEMPLFI+SGAVP1*EKINSG/(KBOLTZ*NDEGFSG)/DELTA/DELTA
    TEMPHFI=SGAVP0*TEMPHFI+SGAVP1*((EKIN-EKINSG)/(KBOLTZ*NDEGFSG)/DELTA/DELTA)
    IF(QSGLDG)THEN
      FRCLF=ONE+SGFFI
      FRCHF=ONE
      SGCFLF=TSGSET/TEMPSG
      SGCFHF=TSGSET/(TSGSET-TEMPLFI)
    ! Estimate reference temperatures
      TREFLFI=TREFLF
      IF(TREFLF<RSMALL)THEN
        TREFLFI=TEMPLFI*SGCFLF
#if KEY_REPDSTR==1 /*rexsgld*/
        if(TRXLF>ZERO)TREFLFI=TRXLF
#endif /* (rexsgld)*/
      ENDIF
      IF(TREFHF<RSMALL)TREFHFI=TSGSET-TREFLFI
      TEMPSGI=SGAVP0*TEMPSGI+SGAVP1*TEMPI*TREFHFI*TEMPLFI/TREFLFI/TEMPHFI
      RETURN
    ENDIF
    AVGFF=SGAVP0*AVGFF+SGAVP1*SUMFF
    AVGDD=SGAVP0*AVGDD+SGAVP1*SUMDD
    AVGGF=SGAVP0*AVGGF+SGAVP1*SUMGF
    AVGGD=SGAVP0*AVGGD+SGAVP1*SUMGD
    AVGPP=SGAVP0*AVGPP+SGAVP1*SUMPP
    AVGGP=SGAVP0*AVGGP+SGAVP1*SUMGP
    ! Estimate SGLD factors
    SGEFLF=(AVGGF/AVGFF+ONE)
    SGEFHF=(AVGGD/(AVGDD)+ONE)
    SGCFLF=(AVGGP/AVGPP)+ONE
    ! Estimate reference temperatures
    TREFLFI=TREFLF
    IF(TREFLF<RSMALL)THEN
      TREFLFI=MAX(TEMPLFI*SGCFLF,TEMPLFI/TEN)
#if KEY_REPDSTR==1 /*rexsgld*/
      if(TRXLF>ZERO)TREFLFI=TRXLF
#endif /* (rexsgld)*/
    ENDIF
    IF(TREFHF<RSMALL)TREFHFI=TSGSET-TREFLFI
    TEMPSGI=SGAVP0*TEMPSGI+SGAVP1*TEMPI*TREFHFI*TEMPLFI/TREFLFI/TEMPHFI
    FRCLF=SGEFLF+SGFFI
    FRCHF=SGEFHF+SGFDI
    SGCFLF=TREFLFI/TEMPLFI
    SGCFHF=TREFHFI/TEMPHFI
    ! Adjust guiding factors if wanted
    IF(QCTSG)THEN
      SGFTI=SGFTI+SGAVP1*(TEMPSG-TEMPSGI)*ABS(TEMPSG-TSGSET)/TSGSET/(TSGSET+TEMPSG)
      SGFTI=MIN(SGFTI,TEN)
      SGFTI=MAX(SGFTI,-TEN)
    ENDIF
    IF(QSGBZ)THEN
    ! Estimate force guiding factor
      IF(QCTSG)THEN
        IF(SGFD*SGFD<RSMALL)SGFDI=SGAVP0*SGFDI+SGAVP1*((TEMPI-TEMPLF)/(TEMPI-TREFLFI)-SGEFHF)
        IF(SGFF*SGFF<RSMALL)SGFFI=SGAVP0*SGFFI+SGAVP1*(TEMPLF/TREFLFI-SGEFLF)
      ELSE
        IF(SGFD*SGFD<RSMALL)THEN
           IF(TREFHF>RSMALL)THEN
              FACT1=TEMPHFI
           ELSE
              FACT1=TSGSET-TEMPLFI
           ENDIF
           SGFDI=SGAVP0*SGFDI+SGAVP1*(FACT1/(TSGSET-TREFLFI)-SGEFHF)
        ENDIF
        IF(SGFF*SGFF<RSMALL)SGFFI=SGAVP0*SGFFI+SGAVP1*(TEMPLFI/TREFLFI-SGEFLF)
      ENDIF
    ENDIF
    RETURN
  END SUBROUTINE SGLDW

  SUBROUTINE SGMDW(ATFRST,ATLAST,DELTA,NDEGF,TEMPI, &
                   IMOVE,AMASS,X,Y,Z,VX,VY,VZ,DX,DY,DZ)
    !
    !     This routine calculates the guiding forces and adds them
    !     to the systematic forces.
    !         Assum SGFT*gamma=constant for all atoms
    !
    use chm_kinds
    use clcg_mod,only:random
    use dimens_fcm
    use number
    use stream
    use consta
    use parallel
    use machutil,only:eclock
    use prssre
    use rndnum
#if KEY_DOMDEC==1
    use domdec_common,only:q_domdec,natoml,atoml
#endif
    !
    INTEGER ATFRST,ATLAST,NDEGF
    real(chm_real)  DELTA,TEMPI
    INTEGER IMOVE(*)
    real(chm_real)  AMASS(*)
    real(chm_real)  X(*),Y(*),Z(*),VX(*),VY(*),VZ(*), &
         DX(*),DY(*),DZ(*)
    !
    INTEGER IG,IA
#if KEY_PARALLEL==1
    real(chm_real) GCARR(30),TIMMER
#endif
    !
    INTEGER I,J,K,IFRST,ILAST
    real(chm_real) DELTA2,GAM
    real(chm_real) VXI,VYI,VZI,VXT,VYT,VZT,GXI,GYI,GZI,DXI,DYI,DZI
    real(chm_real) FXI,FYI,FZI,SGXI,SGYI,SGZI,SGXJ,SGYJ,SGZJ,SGXK,SGYK,SGZK
    real(chm_real) GX2I,GY2I,GZ2I,GX2J,GY2J,GZ2J,GX2K,GY2K,GZ2K
    real(chm_real) FACT0,FACT1,FACT2,AMASSI,SGFTJ
    real(chm_real) A,B,RDX,RDY,RDZ
    real(chm_real) SUMGV,SUMVV,SUMV2,SUMG2,EKIN,EKINSG,EKINSGG
    real(chm_real) SUMFF,SUMDD,SUMGF,SUMGD
    real(chm_real) SUMPP,SUMGP,SUMFP,SUMQQ,SUMHQ,SUMDQ
    real(chm_real) XI,YI,ZI,XCM,YCM,ZCM,XX,XY,XZ,YY,YZ,ZZ
    real(chm_real) FXT,FYT,FZT,PXT,PYT,PZT
    real(chm_real) LXT,LYT,LZT,TXT,TYT,TZT
    real(chm_real) VXCM,VYCM,VZCM,AXCM,AYCM,AZCM
    real(chm_real) OXCM,OYCM,OZCM,QXCM,QYCM,QZCM
    real(chm_real) TCM(3,3),U(3,3),SCR(24),AMOM(6)
    !
    DELTA2=DELTA*DELTA
    !
    IFRST=ATFRST
    ILAST=ATLAST
    IF(ISGSTA > IFRST) IFRST=ISGSTA
    IF(ISGEND < ILAST) ILAST=ISGEND
    !
    !WXW Calculate constraint factor to eliminate energy input from guiding force
    SUMGV=ZERO
    SUMVV=ZERO
    EKIN=ZERO
    EKINSG=ZERO
    EKINSGG=ZERO
    XCM=ZERO
    YCM=ZERO
    ZCM=ZERO
    XX=ZERO
    XY=ZERO
    XZ=ZERO
    YY=ZERO
    YZ=ZERO
    ZZ=ZERO
    PXT=ZERO
    PYT=ZERO
    PZT=ZERO
    LXT=ZERO
    LYT=ZERO
    LZT=ZERO
    FXT=ZERO
    FYT=ZERO
    FZT=ZERO
    TXT=ZERO
    TYT=ZERO
    TZT=ZERO
    K=0
    do ia = atfrst,atlast
#if KEY_DOMDEC==1
       if (q_domdec) then
          i = atoml(ia)
       else
#endif
          i = ia
#if KEY_DOMDEC==1
       endif
#endif
       IF(IMOVE(I) == 0) THEN
          AMASSI=AMASS(I)
          XI=X(I)
          YI=Y(I)
          ZI=Z(I)
          VXI=VX(I)
          VYI=VY(I)
          VZI=VZ(I)
          DXI=DX(I)
          DYI=DY(I)
          DZI=DZ(I)
          IF(I>=ISGSTA .AND. I<=ISGEND)THEN
            SGFTJ=SGFTI*SGWT(I)
            GXI=SGAVG0*SGVX(I)+SGAVG1*VXI
            GYI=SGAVG0*SGVY(I)+SGAVG1*VYI
            GZI=SGAVG0*SGVZ(I)+SGAVG1*VZI
            SGVX(I)=GXI
            SGVY(I)=GYI
            SGVZ(I)=GZI
            FACT1=TIMFAC*AMASSI/DELTA
            FXI=SGAVG0*SGFX(I)+SGAVG1*(SGDX(I)+SHKDX(I))
            FYI=SGAVG0*SGFY(I)+SGAVG1*(SGDY(I)+SHKDY(I))
            FZI=SGAVG0*SGFZ(I)+SGAVG1*(SGDZ(I)+SHKDZ(I))
            SGFX(I)=FXI
            SGFY(I)=FYI
            SGFZ(I)=FZI
            SGXI=SGAVG0*SGKX(I)+SGAVG1*GXI
            SGYI=SGAVG0*SGKY(I)+SGAVG1*GYI
            SGZI=SGAVG0*SGKZ(I)+SGAVG1*GZI
            SGKX(I)=SGXI
            SGKY(I)=SGYI
            SGKZ(I)=SGZI
            EKINSGG=EKINSGG+AMASSI*(SGXI*SGXI+SGYI*SGYI+SGZI*SGZI)
            !SGFTJ=(TFRIC-SGFT)*SGWT(I)
            SGFTJ=(SGFTI)*SGWT(I)
            SGXJ=SGAVG0*SGGX(I)+SGAVG1*FXI
            SGYJ=SGAVG0*SGGY(I)+SGAVG1*FYI
            SGZJ=SGAVG0*SGGZ(I)+SGAVG1*FZI
            SGGX(I)=SGXJ
            SGGY(I)=SGYJ
            SGGZ(I)=SGZJ
            SGXI=(SGDX(I)+SHKDX(I)-FXI)
            SGYI=(SGDY(I)+SHKDY(I)-FYI)
            SGZI=(SGDZ(I)+SHKDZ(I)-FZI)
            SGHX(I)=SGAVP0*SGHX(I)+SGAVP1*(SGXI*VXI+SGYI*VYI+SGZI*VZI)
            SGHY(I)=SGAVP0*SGHY(I)+SGAVP1*AMASSI*TIMFAC*(VXI*VXI+VYI*VYI+VZI*VZI)
            SGGAMMA(I)=SGAVP0*SGGAMMA(I)+  &
              SGAVP1*SGHX(I)/SGHY(I)
            FACT1=SGGAMMA(I)*TIMFAC*AMASSI
            SGXK=FXI-SGXJ
            SGYK=FYI-SGYJ
            SGZK=FZI-SGZJ
            SGXI=(SGFTJ*FACT1)*GXI-SGFFI*SGXK
            SGYI=(SGFTJ*FACT1)*GYI-SGFFI*SGYK
            SGZI=(SGFTJ*FACT1)*GZI-SGFFI*SGZK
            SGXJ=SGXI
            SGYJ=SGYI
            SGZJ=SGZI
            SGDX(I)=-SGXJ
            SGDY(I)=-SGYJ
            SGDZ(I)=-SGZJ
            !write(*,'("SGLDG:",I6,8F10.4)')I,SGWT(I),SGGAMMA(I),RDX,RDY,RDZ,SGXI,SGYI,SGZI
            DXI=DXI-SGXI
            DYI=DYI-SGYI
            DZI=DZI-SGZI
            DX(I)=DXI
            DY(I)=DYI
            DZ(I)=DZI
            FACT2=DELTA2/TWO/AMASSI
            VXT =VXI-DXI*FACT2
            VYT =VYI-DYI*FACT2
            VZT =VZI-DZI*FACT2
            SUMGV=SUMGV+(SGXJ*VXT+SGYJ*VYT+SGZJ*VZT)
            SUMVV=SUMVV+AMASSI*(VXT*VXT+VYT*VYT+VZT*VZT)
            EKIN=EKIN+AMASSI*(VXI*VXI+VYI*VYI+VZI*VZI)
            EKINSG=EKINSG+AMASSI*(GXI*GXI+GYI*GYI+GZI*GZI)
          ENDIF
          XCM=XCM+XI*AMASSI
          YCM=YCM+YI*AMASSI
          ZCM=ZCM+ZI*AMASSI
          XX=XX+XI*XI*AMASSI
          XY=XY+XI*YI*AMASSI
          XZ=XZ+XI*ZI*AMASSI
          YY=YY+YI*YI*AMASSI
          YZ=YZ+YI*ZI*AMASSI
          ZZ=ZZ+ZI*ZI*AMASSI
          PXT=PXT+AMASSI*VXI
          PYT=PYT+AMASSI*VYI
          PZT=PZT+AMASSI*VZI
          LXT=LXT+(YI*VZI-ZI*VYI)*AMASSI
          LYT=LYT+(ZI*VXI-XI*VZI)*AMASSI
          LZT=LZT+(XI*VYI-YI*VXI)*AMASSI
          FXT=FXT+DXI
          FYT=FYT+DYI
          FZT=FZT+DZI
          TXT=TXT+(YI*DZI-ZI*DYI)
          TYT=TYT+(ZI*DXI-XI*DZI)
          TZT=TZT+(XI*DYI-YI*DXI)
       ENDIF
    ENDDO
    !
#if KEY_PARALLEL==1
    CALL PSYNC()
    TIMMER=ECLOCK()
    GCARR(1)=SUMGV
    GCARR(2)=SUMVV
    GCARR(3)=EKIN
    GCARR(4)=EKINSG
    GCARR(5) = XCM
    GCARR(6) = YCM
    GCARR(7) = ZCM
    GCARR(8) = XX
    GCARR(9) = XY
    GCARR(10) = XZ
    GCARR(11) = YY
    GCARR(12) = YZ
    GCARR(13) = ZZ
    GCARR(14)=PXT
    GCARR(15)=PYT
    GCARR(16)=PZT
    GCARR(17)=LXT
    GCARR(18)=LYT
    GCARR(19)=LZT
    GCARR(20)=FXT
    GCARR(21)=FYT
    GCARR(22)=FZT
    GCARR(23) = TXT
    GCARR(24) = TYT
    GCARR(25) = TZT
    GCARR(26)=EKINSGG
    CALL GCOMB(GCARR,26)
    SUMGV=GCARR(1)
    SUMVV=GCARR(2)
    EKIN=GCARR(3)
    EKINSG=GCARR(4)
    XCM=GCARR(5)
    YCM=GCARR(6)
    ZCM=GCARR(7)
    XX=GCARR(8)
    XY=GCARR(9)
    XZ=GCARR(10)
    YY=GCARR(11)
    YZ=GCARR(12)
    ZZ=GCARR(13)
    PXT=GCARR(14)
    PYT=GCARR(15)
    PZT=GCARR(16)
    LXT=GCARR(17)
    LYT=GCARR(18)
    LZT=GCARR(19)
    FXT=GCARR(20)
    FYT=GCARR(21)
    FZT=GCARR(22)
    TXT=GCARR(23)
    TYT=GCARR(24)
    TZT=GCARR(25)
    EKINSGG=GCARR(26)
    TMERI(TIMGCOMM)=TMERI(TIMGCOMM)+ECLOCK()-TIMMER
    TIMMER=ECLOCK()
#endif
    !
    ! Translation
    IF(QSGSTOPT)THEN
      !  net velocity
      VXCM=PXT/TOTMASS
      VYCM=PYT/TOTMASS
      VZCM=PZT/TOTMASS
      !  net accelation
      AXCM=FXT/TOTMASS
      AYCM=FYT/TOTMASS
      AZCM=FZT/TOTMASS
    ELSE
      VXCM=ZERO
      VYCM=ZERO
      VZCM=ZERO
      AXCM=ZERO
      AYCM=ZERO
      AZCM=ZERO
    ENDIF
    ! Rotation
    IF(QSGSTOPR)THEN
      XCM=XCM/TOTMASS
      YCM=YCM/TOTMASS
      ZCM=ZCM/TOTMASS
      XX=XX-XCM*XCM*TOTMASS
      XY=XY-XCM*YCM*TOTMASS
      XZ=XZ-XCM*ZCM*TOTMASS
      YY=YY-YCM*YCM*TOTMASS
      YZ=YZ-YCM*ZCM*TOTMASS
      ZZ=ZZ-ZCM*ZCM*TOTMASS
      LXT=LXT-(YCM*PZT-ZCM*PYT)
      LYT=LYT-(ZCM*PXT-XCM*PZT)
      LZT=LZT-(XCM*PYT-YCM*PXT)
      TXT=TXT-(YCM*FZT-ZCM*FYT)
      TYT=TYT-(ZCM*FXT-XCM*FZT)
      TZT=TZT-(XCM*FYT-YCM*FXT)
    !
      AMOM(1)=YY+ZZ
      AMOM(2)=-XY
      AMOM(3)=-XZ
      AMOM(4)=XX+ZZ
      AMOM(5)=-YZ
      AMOM(6)=XX+YY
    !
      CALL DIAGQ(3,3,AMOM,U,SCR(4),SCR(7),SCR(10),SCR(13),SCR(1),  &
         SCR(16),SCR(19),SCR(22),0)
    !
      DO I=1,3
         IF(ABS(SCR(I)) < TENM5) SCR(I)=SIGN(TENM5,SCR(I))
         SCR(I)=ONE/SCR(I)
      ENDDO
    !
      DO I=1,3
         DO J=1,3
            TCM(I,J)=0.0
            DO K=1,3
               TCM(I,J)=TCM(I,J)+U(I,K)*U(J,K)*SCR(K)
            ENDDO
         ENDDO
      ENDDO
    ! angular velocity of COM
      OXCM=LXT*TCM(1,1)+LYT*TCM(1,2)+LZT*TCM(1,3)
      OYCM=LXT*TCM(2,1)+LYT*TCM(2,2)+LZT*TCM(2,3)
      OZCM=LXT*TCM(3,1)+LYT*TCM(3,2)+LZT*TCM(3,3)
    ! angular acceleration of COM
      QXCM=TXT*TCM(1,1)+TYT*TCM(1,2)+TZT*TCM(1,3)
      QYCM=TXT*TCM(2,1)+TYT*TCM(2,2)+TZT*TCM(2,3)
      QZCM=TXT*TCM(3,1)+TYT*TCM(3,2)+TZT*TCM(3,3)
    ELSE
      OXCM=ZERO
      OYCM=ZERO
      OZCM=ZERO
      QXCM=ZERO
      QYCM=ZERO
      QZCM=ZERO
    ENDIF
    !
    SGSCAL=DELTA2*SUMGV/(SUMVV-HALF*DELTA2*SUMGV)
    !
    FACT0=TWO/(TWO+SGSCAL)
    !
    SUMFF=ZERO
    SUMDD=ZERO
    SUMGF=ZERO
    SUMGD=ZERO
    SUMPP=ZERO
    SUMGP=ZERO
    SUMQQ=ZERO
    do ia = atfrst,atlast
#if KEY_DOMDEC==1
       if (q_domdec) then
          i = atoml(ia)
       else
#endif
          i = ia
#if KEY_DOMDEC==1
       endif
#endif
       IF(IMOVE(I) == 0) THEN
          AMASSI=AMASS(I)
          XI=X(I)-XCM
          YI=Y(I)-YCM
          ZI=Z(I)-ZCM
          VXI=VX(I)-VXCM-OYCM*ZI+OZCM*YI
          VYI=VY(I)-VYCM-OZCM*XI+OXCM*ZI
          VZI=VZ(I)-VZCM-OXCM*YI+OYCM*XI
          VX(I)=VXI
          VY(I)=VYI
          VZ(I)=VZI
          DXI=DX(I)-AMASSI*(AXCM+QYCM*ZI-QZCM*YI)
          DYI=DY(I)-AMASSI*(AYCM+QZCM*XI-QXCM*ZI)
          DZI=DZ(I)-AMASSI*(AZCM+QXCM*YI-QYCM*XI)
          IF(I<ISGSTA .OR. I>ISGEND)THEN
            DX(I)=DXI
            DY(I)=DYI
            DZ(I)=DZI
            CYCLE
          ENDIF
          FACT1=AMASSI/DELTA2
          FACT2=SGSCAL*FACT0*FACT1
          DXI=FACT0*DXI+FACT2*VXI
          DYI=FACT0*DYI+FACT2*VYI
          DZI=FACT0*DZI+FACT2*VZI
          DX(I)=DXI
          DY(I)=DYI
          DZ(I)=DZI
          FACT2=HALF*FACT0/FACT1
          VXT =FACT0*VXI-DXI*FACT2
          VYT =FACT0*VYI-DYI*FACT2
          VZT =FACT0*VZI-DZI*FACT2
          IF(QSGLDG)CYCLE
          ! Accumulation for SGMD parameter calculation
          FXI=SGFX(I)
          FYI=SGFY(I)
          FZI=SGFZ(I)
          DXI=SGDX(I)
          DYI=SGDY(I)
          DZI=SGDZ(I)
          GXI=SGVX(I)
          GYI=SGVY(I)
          GZI=SGVZ(I)
          GAM=SGFTI*SGWT(I)*TIMFAC*DELTA
          SGXI=FACT1*((SGSCAL)*VXT-GAM*GXI)
          SGYI=FACT1*((SGSCAL)*VYT-GAM*GYI)
          SGZI=FACT1*((SGSCAL)*VZT-GAM*GZI)
          SGXK=SGXI+SGFDI*DXI+(SGFFI-SGFDI)*FXI
          SGYK=SGYI+SGFDI*DYI+(SGFFI-SGFDI)*FYI
          SGZK=SGZI+SGFDI*DZI+(SGFFI-SGFDI)*FZI
          SGXJ=SGXK-SGFDI*(DXI-FXI)
          SGYJ=SGYK-SGFDI*(DYI-FYI)
          SGZJ=SGZK-SGFDI*(DZI-FZI)
          ! save guiding forces for virial calculation
          SGDX(I)=SGXK
          SGDY(I)=SGYK
          SGDZ(I)=SGZK
          VXI=VXT-GXI
          VYI=VYT-GYI
          VZI=VZT-GZI
          DXI=DXI-FXI
          DYI=DYI-FYI
          DZI=DZI-FZI
          SUMFF=SUMFF+(FXI*FXI+FYI*FYI+FZI*FZI)
          SUMDD=SUMDD+(DXI*DXI+DYI*DYI+DZI*DZI)
          SUMGF=SUMGF+(FXI*SGXI+FYI*SGYI+FZI*SGZI)
          SUMGD=SUMGD+(DXI*SGXJ+DYI*SGYJ+DZI*SGZJ)
         ! lamda derivatibes
          GX2K=SGKX(I)
          GY2K=SGKY(I)
          GZ2K=SGKZ(I)
          GX2J=SGHX(I)-SGAVP1*SGCGHF*GX2K
          GY2J=SGHY(I)-SGAVP1*SGCGHF*GY2K
          GZ2J=SGHZ(I)-SGAVP1*SGCGHF*GZ2K
          GX2I=SGGX(I)-SGCGHF*GX2K
          GY2I=SGGY(I)-SGCGHF*GY2K
          GZ2I=SGGZ(I)-SGCGHF*GZ2K
          GX2I=SGAVP0*((TWO-SGSCAL)*GX2I+  &
               GAM*(ONE+SGAVG0)*GX2J-TWO*DELTA*SGXK)/(TWO+SGSCAL-SGAVG1*GAM)
          GY2I=SGAVP0*((TWO-SGSCAL)*GY2I+  &
               GAM*(ONE+SGAVG0)*GY2J-TWO*DELTA*SGYK)/(TWO+SGSCAL-SGAVG1*GAM)
          GZ2I=SGAVP0*((TWO-SGSCAL)*GZ2I+  &
               GAM*(ONE+SGAVG0)*GZ2J-TWO*DELTA*SGZK)/(TWO+SGSCAL-SGAVG1*GAM)
          GX2K=SGAVP0*SGKX(I)+SGAVP1*GX2I
          GY2K=SGAVP0*SGKY(I)+SGAVP1*GY2I
          GZ2K=SGAVP0*SGKZ(I)+SGAVP1*GZ2I
          GX2J=SGAVP0*(SGAVG0*GX2J+SGAVG1*(GX2I))
          GY2J=SGAVP0*(SGAVG0*GY2J+SGAVG1*(GY2I))
          GZ2J=SGAVP0*(SGAVG0*GZ2J+SGAVG1*(GZ2I))
          SGKX(I)=GX2K
          SGKY(I)=GY2K
          SGKZ(I)=GZ2K
          SGHX(I)=GX2J
          SGHY(I)=GY2J
          SGHZ(I)=GZ2J
          SGGX(I)=GX2I
          SGGY(I)=GY2I
          SGGZ(I)=GZ2I
          ! lmda derivative of temperature
          SUMQQ=SUMQQ+(GXI*GX2J+GYI*GY2J+GZI*GZ2J)
          SUMPP=SUMPP+(VXT*GX2I+VYT*GY2I+VZT*GZ2I)
          SUMGP=SUMGP+(VXT*GX2K+VYT*GY2K+VZT*GZ2K)
       ENDIF
    ENDDO
#if KEY_PARALLEL==1
    CALL PSYNC()
    TIMMER=ECLOCK()
    GCARR(1)=SUMFF
    GCARR(2)=SUMDD
    GCARR(3)=SUMGF
    GCARR(4)=SUMGD
    GCARR(5)=SUMPP
    GCARR(6)=SUMGP
    GCARR(7)=SUMQQ
    GCARR(8)=VIRSG
    CALL GCOMB(GCARR,8)
    SUMFF=GCARR(1)
    SUMDD=GCARR(2)
    SUMGF=GCARR(3)
    SUMGD=GCARR(4)
    SUMPP=GCARR(5)
    SUMGP=GCARR(6)
    SUMQQ=GCARR(7)
    VIRSG=GCARR(8)
    TMERI(TIMGCOMM)=TMERI(TIMGCOMM)+ECLOCK()-TIMMER
    TIMMER=ECLOCK()
#endif
    ! Estimate momentum guiding factor and guiding temperature
    TEMPLFI=SGAVP0*TEMPLFI+SGAVP1*EKINSG/(KBOLTZ*NDEGFSG)/DELTA/DELTA
    TEMPHFI=SGAVP0*TEMPHFI+SGAVP1*((EKIN-EKINSG)/(KBOLTZ*NDEGFSG)/DELTA/DELTA)
    IF(QSGLDG)THEN
      FRCLF=ONE+SGFFI
      FRCHF=ONE
      SGCFLF=TSGSET/TEMPSG
      SGCFHF=TSGSET/(TSGSET-TEMPLFI)

    ! Estimate reference temperatures
      TREFLFI=TREFLF
      IF(TREFLF<RSMALL)THEN
        TREFLFI=TEMPLFI*SGCFLF
#if KEY_REPDSTR==1 /*rexsgld*/
        if(TRXLF>ZERO)TREFLFI=TRXLF
#endif /* (rexsgld)*/
      ENDIF
      IF(TREFHF<RSMALL)TREFHFI=TSGSET-TREFLFI
      TEMPSGI=SGAVP0*TEMPSGI+SGAVP1*TEMPI*TREFHFI*TEMPLFI/TREFLFI/TEMPHFI
      RETURN
    ENDIF
    AVGFF=SGAVP0*AVGFF+SGAVP1*SUMFF
    AVGDD=SGAVP0*AVGDD+SGAVP1*SUMDD
    AVGGF=SGAVP0*AVGGF+SGAVP1*SUMGF
    AVGGD=SGAVP0*AVGGD+SGAVP1*SUMGD
    AVGPP=SGAVP0*AVGPP+SGAVP1*SUMPP
    AVGGP=SGAVP0*AVGGP+SGAVP1*SUMGP
    AVGQQ=SGAVP0*AVGQQ+SGAVP1*SUMQQ*TWO/(KBOLTZ*NDEGFSG)/DELTA
    ! Estimate SGLD factors
    SGEFLF=(AVGGF/AVGFF+ONE)
    SGEFHF=(AVGGD/(AVGDD)+ONE)
    ! Estimate SGMD factors
    IF(SUMGP*SUMGP>RSMALL)THEN
      SGCGHF=SGAVG1*SUMPP/SUMGP
    ELSE
      SGCGHF=SGAVP0*SGCGHF
    ENDIF
    ! Estimate reference temperatures
    TREFLFI=TREFLF
    TREFHFI=TREFHF
    IF(TREFLF<RSMALL)THEN
      TREFLFI=TEMPLFI
#if KEY_REPDSTR==1 /*rexsgld*/
       if(TRXLF>ZERO)TREFLFI=TRXLF
#endif /* (rexsgld)*/
    ENDIF
    IF(TREFHF<RSMALL)TREFHFI=TSGSET-TREFLFI
    IF(TREFHF<RSMALL)TEMPHFI=TSGSET-TEMPLFI
    TEMPSGI=SGAVP0*TEMPSGI+SGAVP1*TEMPI*TREFHFI*TEMPLFI/TREFLFI/TEMPHFI
    FRCLF=SGEFLF+SGFFI
    FRCHF=SGEFHF+SGFDI
    ! Adjust guiding factors if wanted
    IF(QCTSG)THEN
      SGFTI=SGFTI+SGAVP1*(TEMPSG-TEMPSGI)*ABS(TEMPSG-TSGSET)/TSGSET/(TSGSET+TEMPSG)
      SGFTI=MIN(SGFTI,TEN)
      SGFTI=MAX(SGFTI,-TEN)
    ENDIF
    IF(QSGBZ)THEN
    ! Estimate force guiding factor
      IF(QCTSG)THEN
        IF(SGFD*SGFD<RSMALL)SGFDI=SGAVP0*SGFDI+SGAVP1*((TEMPI-TEMPLF)/(TEMPI-TREFLFI)-SGEFHF)
        IF(SGFF*SGFF<RSMALL)SGFFI=SGAVP0*SGFFI+SGAVP1*(TEMPLF/TREFLFI-SGEFLF)
      ELSE
        IF(SGFD*SGFD<RSMALL)THEN
           IF(TREFHF>RSMALL)THEN
              FACT1=TEMPHFI
           ELSE
              FACT1=TSGSET-TEMPLFI
           ENDIF
           SGFDI=SGAVP0*SGFDI+SGAVP1*(FACT1/(TSGSET-TREFLFI)-SGEFHF)
        ENDIF
        IF(SGFF*SGFF<RSMALL)SGFFI=SGAVP0*SGFFI+SGAVP1*(TEMPLFI/TREFLFI-SGEFLF)
      ENDIF
      SGSCAL=SGFDI
    ENDIF
    RETURN
  END SUBROUTINE SGMDW

  SUBROUTINE SGLDG(ATFRST,ATLAST,DELTA,NDEGF, &
               INLCKP,QNHLANG,isDrude,NHGAMMA,NHGAMMAD,  &
                   IMOVE,AMASS,X,Y,Z,VX,VY,VZ,DX,DY,DZ)
    !
    !     This routine calculates the guiding forces and adds them
    !     to the systematic forces.
    !         Assum SGFT*gamma=constant for all atoms
    !
    use chm_kinds
    use dimens_fcm
    use number
    use stream
    use consta
    use cnst_fcm
#if KEY_PARALLEL==1
    use parallel
#endif
    use machutil,only:eclock
    use prssre
#if KEY_DOMDEC==1
    use domdec_common,only:q_domdec,natoml,atoml
#endif
    !
    INTEGER ATFRST,ATLAST,NDEGF
    real(chm_real)  DELTA,NHGAMMA,NHGAMMAD
    INTEGER INLCKP(*)
    LOGICAL QNHLANG(*),isDrude(*)
    INTEGER IMOVE(*)
    real(chm_real)  AMASS(*),X(*),Y(*),Z(*),VX(*),VY(*),VZ(*), &
         DX(*),DY(*),DZ(*)
    !
#if KEY_PARALLEL==1
    real(chm_real) GCARR(31),TIMMER
#endif
    !
    INTEGER I,J,K,IFRST,ILAST,IA,TI
    real(chm_real) XI,YI,ZI,VXI,VYI,VZI,VXT,VYT,VZT
    real(chm_real) FXI,FYI,FZI,GXI,GYI,GZI,DXI,DYI,DZI
    real(chm_real) SGXI,SGYI,SGZI,SGXJ,SGYJ,SGZJ,SGXK,SGYK,SGZK
    real(chm_real) GX2I,GY2I,GZ2I,GX2J,GY2J,GZ2J,GX2K,GY2K,GZ2K
    real(chm_real) GAM,GAMM,FACT0,FACT1,FACT2,DTEMP,SGFTJ
    real(chm_real) AMASSI,DELTA2
    real(chm_real) SUMGV,SUMVV,SUMV2,SUMG2,EKIN,EKINSG,EKINSGG
    real(chm_real) SUMFF,SUMDD,SUMGF,SUMGD
    real(chm_real) SUMPP,SUMGP,SUMFP,SUMQQ,SUMHQ,SUMDQ
    real(chm_real) XCM,YCM,ZCM,XX,XY,XZ,YY,YZ,ZZ
    real(chm_real) FXT,FYT,FZT,PXT,PYT,PZT
    real(chm_real) LXT,LYT,LZT,TXT,TYT,TZT
    real(chm_real) VXCM,VYCM,VZCM,AXCM,AYCM,AZCM
    real(chm_real) VXG,VYG,VZG,OXG,OYG,OZG
    real(chm_real) FXG,FYG,FZG,AXG,AYG,AZG
    real(chm_real) TXG,TYG,TZG,GXCM,GYCM,GZCM
    real(chm_real) OXCM,OYCM,OZCM,QXCM,QYCM,QZCM
    real(chm_real) TCM(3,3),U(3,3),SCR(24),AMOM(6)
    !
    DELTA2=DELTA*DELTA
    !
    IFRST=ATFRST
    ILAST=ATLAST
    IF(ISGSTA > IFRST) IFRST=ISGSTA
    IF(ISGEND < ILAST) ILAST=ISGEND
    !
    !WXW Calculate constraint factor to eliminate energy input from guiding force
    SUMGV=ZERO
    SUMVV=ZERO
    XCM=ZERO
    YCM=ZERO
    ZCM=ZERO
    XX=ZERO
    XY=ZERO
    XZ=ZERO
    YY=ZERO
    YZ=ZERO
    ZZ=ZERO
    PXT=ZERO
    PYT=ZERO
    PZT=ZERO
    LXT=ZERO
    LYT=ZERO
    LZT=ZERO
    FXT=ZERO
    FYT=ZERO
    FZT=ZERO
    TXT=ZERO
    TYT=ZERO
    TZT=ZERO
    FXG=ZERO
    FYG=ZERO
    FZG=ZERO
    TXG=ZERO
    TYG=ZERO
    TZG=ZERO
    do ia = atfrst,atlast
#if KEY_DOMDEC==1
       if (q_domdec) then
          i = atoml(ia)
       else
#endif
          i = ia
#if KEY_DOMDEC==1
       endif
#endif
       IF(IMOVE(I) == 0) THEN
          AMASSI=AMASS(I)
          XI=X(I)
          YI=Y(I)
          ZI=Z(I)
          VXI=VX(I)
          VYI=VY(I)
          VZI=VZ(I)
          DXI=DX(I)
          DYI=DY(I)
          DZI=DZ(I)
          IF(I>=ISGSTA .AND. I<=ISGEND)THEN
            SGFTJ=SGFTI*SGWT(I)
            ti = INLCKP(i)
            j = i + 1
            if (QNHLANG(ti).and.(isDrude(i).or.isDrude(j)))then
              GAM=NHGAMMAD
            else
              GAM=NHGAMMA
            endif
            GAMM=TWO/(TWO+GAM*DELTA)
            GXI=SGAVG0*SGVX(I)+SGAVG1*VXI
            GYI=SGAVG0*SGVY(I)+SGAVG1*VYI
            GZI=SGAVG0*SGVZ(I)+SGAVG1*VZI
            SGVX(I)=GXI
            SGVY(I)=GYI
            SGVZ(I)=GZI
            FACT1=SGFTJ*GAM*AMASSI
            FXI=SGAVG0*SGFX(I)+SGAVG1*(SGDX(I)+SHKDX(I))
            FYI=SGAVG0*SGFY(I)+SGAVG1*(SGDY(I)+SHKDY(I))
            FZI=SGAVG0*SGFZ(I)+SGAVG1*(SGDZ(I)+SHKDZ(I))
            SGFX(I)=FXI
            SGFY(I)=FYI
            SGFZ(I)=FZI
            IF(QSGLDG)THEN
              SGXJ=SGAVG0*SGGX(I)+SGAVG1*FXI
              SGYJ=SGAVG0*SGGY(I)+SGAVG1*FYI
              SGZJ=SGAVG0*SGGZ(I)+SGAVG1*FZI
              SGGX(I)=SGXJ
              SGGY(I)=SGYJ
              SGGZ(I)=SGZJ
              SGXK=SGAVG0*SGKX(I)+SGAVG1*(DXI-SGDX(I))
              SGYK=SGAVG0*SGKY(I)+SGAVG1*(DYI-SGDY(I))
              SGZK=SGAVG0*SGKZ(I)+SGAVG1*(DZI-SGDZ(I))
              SGKX(I)=SGXK
              SGKY(I)=SGYK
              SGKZ(I)=SGZK
              SGXI=(SGDX(I)+SHKDX(I)-FXI)
              SGYI=(SGDY(I)+SHKDY(I)-FYI)
              SGZI=(SGDZ(I)+SHKDZ(I)-FZI)
              SGHX(I)=SGAVP0*SGHX(I)+SGAVP1*(SGXI*VXI+SGYI*VYI+SGZI*VZI)
              SGHY(I)=SGAVP0*SGHY(I)+SGAVP1*AMASSI*TIMFAC*(VXI*VXI+VYI*VYI+VZI*VZI)
              SGGAMMA(I)=SGAVP0*SGGAMMA(I)+SGAVP1*SGHX(I)/SGHY(I)
              FACT2=PSGLDG*SGGAMMA(I)*TIMFAC*AMASSI*SGWT(I)
              SGXJ=FXI-SGXJ
              SGYJ=FYI-SGYJ
              SGZJ=FZI-SGZJ
              SGXJ=FACT2*GXI-SGFFI*SGXJ
              SGYJ=FACT2*GYI-SGFFI*SGYJ
              SGZJ=FACT2*GZI-SGFFI*SGZJ
              SGXI=FACT1*GXI+SGXJ-FSGLDG*SGXK
              SGYI=FACT1*GYI+SGYJ-FSGLDG*SGYK
              SGZI=FACT1*GZI+SGZJ-FSGLDG*SGZK
            ELSE
              FXI=(SGFFI-SGFDI)*FXI+SGFDI*SGDX(I)
              FYI=(SGFFI-SGFDI)*FYI+SGFDI*SGDY(I)
              FZI=(SGFFI-SGFDI)*FZI+SGFDI*SGDZ(I)
              SGXI=FACT1*GXI-FXI
              SGYI=FACT1*GYI-FYI
              SGZI=FACT1*GZI-FZI
              SGXJ=SGXI
              SGYJ=SGYI
              SGZJ=SGZI
            ENDIF
            DXI=DXI-SGXI
            DYI=DYI-SGYI
            DZI=DZI-SGZI
            FXG=FXG-SGXI
            FYG=FYG-SGYI
            FZG=FZG-SGZI
            TXG=TXG-(YI*SGZI-ZI*SGYI)
            TYG=TYG-(ZI*SGXI-XI*SGZI)
            TZG=TZG-(XI*SGYI-YI*SGXI)
            DX(I)=DXI
            DY(I)=DYI
            DZ(I)=DZI
            !FACT2=DELTA/TWO/AMASSI
            !VXT =VX(I)-DXI*FACT2
            !VYT =VY(I)-DYI*FACT2
            !VZT =VZ(I)-DZI*FACT2
            VXT =VXI
            VYT =VYI
            VZT =VZI
            SUMGV=SUMGV+(SGXJ*VXT+SGYJ*VYT+SGZJ*VZT)*GAMM
            SUMVV=SUMVV+AMASSI*(VXT*VXT+VYT*VYT+VZT*VZT)*GAMM*GAMM
          ENDIF
          XCM=XCM+XI*AMASSI
          YCM=YCM+YI*AMASSI
          ZCM=ZCM+ZI*AMASSI
          XX=XX+XI*XI*AMASSI
          XY=XY+XI*YI*AMASSI
          XZ=XZ+XI*ZI*AMASSI
          YY=YY+YI*YI*AMASSI
          YZ=YZ+YI*ZI*AMASSI
          ZZ=ZZ+ZI*ZI*AMASSI
          PXT=PXT+AMASSI*VXI
          PYT=PYT+AMASSI*VYI
          PZT=PZT+AMASSI*VZI
          LXT=LXT+(YI*VZI-ZI*VYI)*AMASSI
          LYT=LYT+(ZI*VXI-XI*VZI)*AMASSI
          LZT=LZT+(XI*VYI-YI*VXI)*AMASSI
          FXT=FXT+DXI
          FYT=FYT+DYI
          FZT=FZT+DZI
          TXT=TXT+(YI*DZI-ZI*DYI)
          TYT=TYT+(ZI*DXI-XI*DZI)
          TZT=TZT+(XI*DYI-YI*DXI)
       ENDIF
    ENDDO
    !
#if KEY_PARALLEL==1
    CALL PSYNC()
    TIMMER=ECLOCK()
    GCARR(1)=SUMGV
    GCARR(2)=SUMVV
    GCARR(3) = XCM
    GCARR(4) = YCM
    GCARR(5) = ZCM
    GCARR(6) = XX
    GCARR(7) = XY
    GCARR(8) = XZ
    GCARR(9) = YY
    GCARR(10) = YZ
    GCARR(11) = ZZ
    GCARR(12)=PXT
    GCARR(13)=PYT
    GCARR(14)=PZT
    GCARR(15)=LXT
    GCARR(16)=LYT
    GCARR(17)=LZT
    GCARR(18)=FXT
    GCARR(19)=FYT
    GCARR(20)=FZT
    GCARR(21) = TXT
    GCARR(22) = TYT
    GCARR(23) = TZT
    GCARR(24)=FXG
    GCARR(25)=FYG
    GCARR(26)=FZG
    GCARR(27) = TXG
    GCARR(28) = TYG
    GCARR(29) = TZG
    CALL GCOMB(GCARR,29)
    SUMGV=GCARR(1)
    SUMVV=GCARR(2)
    XCM=GCARR(3)
    YCM=GCARR(4)
    ZCM=GCARR(5)
    XX=GCARR(6)
    XY=GCARR(7)
    XZ=GCARR(8)
    YY=GCARR(9)
    YZ=GCARR(10)
    ZZ=GCARR(11)
    PXT=GCARR(12)
    PYT=GCARR(13)
    PZT=GCARR(14)
    LXT=GCARR(15)
    LYT=GCARR(16)
    LZT=GCARR(17)
    FXT=GCARR(18)
    FYT=GCARR(19)
    FZT=GCARR(20)
    TXT=GCARR(21)
    TYT=GCARR(22)
    TZT=GCARR(23)
    FXG=GCARR(24)
    FYG=GCARR(25)
    FZG=GCARR(26)
    TXG=GCARR(27)
    TYG=GCARR(28)
    TZG=GCARR(29)
    TMERI(TIMGCOMM)=TMERI(TIMGCOMM)+ECLOCK()-TIMMER
    TIMMER=ECLOCK()
#endif
    ! rotation inertia
    IF(QSGSTOPR.OR.QSGNONET)THEN
      XCM=XCM/TOTMASS
      YCM=YCM/TOTMASS
      ZCM=ZCM/TOTMASS
      LXT=LXT-(YCM*PZT-ZCM*PYT)
      LYT=LYT-(ZCM*PXT-XCM*PZT)
      LZT=LZT-(XCM*PYT-YCM*PXT)
      TXT=TXT-(YCM*FZT-ZCM*FYT)
      TYT=TYT-(ZCM*FXT-XCM*FZT)
      TZT=TZT-(XCM*FYT-YCM*FXT)
      XX=XX-XCM*XCM*TOTMASS
      XY=XY-XCM*YCM*TOTMASS
      XZ=XZ-XCM*ZCM*TOTMASS
      YY=YY-YCM*YCM*TOTMASS
      YZ=YZ-YCM*ZCM*TOTMASS
      ZZ=ZZ-ZCM*ZCM*TOTMASS
    !
      AMOM(1)=YY+ZZ
      AMOM(2)=-XY
      AMOM(3)=-XZ
      AMOM(4)=XX+ZZ
      AMOM(5)=-YZ
      AMOM(6)=XX+YY
    !
      CALL DIAGQ(3,3,AMOM,U,SCR(4),SCR(7),SCR(10),SCR(13),SCR(1),  &
         SCR(16),SCR(19),SCR(22),0)
    !
      DO I=1,3
         IF(ABS(SCR(I)) < TENM5) SCR(I)=SIGN(TENM5,SCR(I))
         SCR(I)=ONE/SCR(I)
      ENDDO
    !
      DO I=1,3
         DO J=1,3
            TCM(I,J)=0.0
            DO K=1,3
               TCM(I,J)=TCM(I,J)+U(I,K)*U(J,K)*SCR(K)
            ENDDO
         ENDDO
      ENDDO
    ENDIF
    ! net translational and rotational guiding forces
    IF(QSGNONET)THEN
      !  net accelation
      AXG=FXG/TOTMASS
      AYG=FYG/TOTMASS
      AZG=FZG/TOTMASS
      TXG=TXG-(YCM*FZG-ZCM*FYG)
      TYG=TYG-(ZCM*FXG-XCM*FZG)
      TZG=TZG-(XCM*FYG-YCM*FXG)
    ! angular acceleration of COM
      GXCM=TXG*TCM(1,1)+TYG*TCM(1,2)+TZG*TCM(1,3)
      GYCM=TXG*TCM(2,1)+TYG*TCM(2,2)+TZG*TCM(2,3)
      GZCM=TXG*TCM(3,1)+TYG*TCM(3,2)+TZG*TCM(3,3)
      !  net velocity for local average
      VXG=PXT/TOTMASS
      VYG=PYT/TOTMASS
      VZG=PZT/TOTMASS
    ! angular velocity of COM for local average
      OXG=LXT*TCM(1,1)+LYT*TCM(1,2)+LZT*TCM(1,3)
      OYG=LXT*TCM(2,1)+LYT*TCM(2,2)+LZT*TCM(2,3)
      OZG=LXT*TCM(3,1)+LYT*TCM(3,2)+LZT*TCM(3,3)
    ELSE
      AXG=ZERO
      AYG=ZERO
      AZG=ZERO
      GXCM=ZERO
      GYCM=ZERO
      GZCM=ZERO
      !  net velocity for local average
      VXG=ZERO
      VYG=ZERO
      VZG=ZERO
    ! angular velocity of COM for local average
      OXG=ZERO
      OYG=ZERO
      OZG=ZERO
    ENDIF
    ! Translation
    IF(QSGSTOPT)THEN
      !  net velocity
      VXCM=PXT/TOTMASS
      VYCM=PYT/TOTMASS
      VZCM=PZT/TOTMASS
      ! net accelation
      AXCM=FXT/TOTMASS
      AYCM=FYT/TOTMASS
      AZCM=FZT/TOTMASS
      !  net velocity for local average
      VXG=ZERO
      VYG=ZERO
      VZG=ZERO
    ELSE
      VXCM=ZERO
      VYCM=ZERO
      VZCM=ZERO
      AXCM=AXG
      AYCM=AYG
      AZCM=AZG
    ENDIF
    ! Rotation
    IF(QSGSTOPR)THEN
    ! angular velocity of COM
      OXCM=LXT*TCM(1,1)+LYT*TCM(1,2)+LZT*TCM(1,3)
      OYCM=LXT*TCM(2,1)+LYT*TCM(2,2)+LZT*TCM(2,3)
      OZCM=LXT*TCM(3,1)+LYT*TCM(3,2)+LZT*TCM(3,3)
    ! angular acceleration of COM
      QXCM=TXT*TCM(1,1)+TYT*TCM(1,2)+TZT*TCM(1,3)
      QYCM=TXT*TCM(2,1)+TYT*TCM(2,2)+TZT*TCM(2,3)
      QZCM=TXT*TCM(3,1)+TYT*TCM(3,2)+TZT*TCM(3,3)
    ! angular velocity of COM for local average
      OXG=ZERO
      OYG=ZERO
      OZG=ZERO
    ELSE
      OXCM=ZERO
      OYCM=ZERO
      OZCM=ZERO
      QXCM=GXCM
      QYCM=GYCM
      QZCM=GZCM
    ENDIF
    !
    !SGSCAL=DELTA*SUMGV/(SUMVV-HALF*DELTA*SUMGV)
    SGSCAL=DELTA*SUMGV/(SUMVV)
    !IF(QSGLDG .OR. QSGBZ)SGSCAL=ZERO
    !
    EKIN=ZERO
    EKINSG=ZERO
    EKINSGG=ZERO
    SUMFF=ZERO
    SUMDD=ZERO
    SUMGF=ZERO
    SUMGD=ZERO
    SUMPP=ZERO
    SUMGP=ZERO
    do ia = atfrst,atlast
#if KEY_DOMDEC==1
       if (q_domdec) then
          i = atoml(ia)
       else
#endif
          i = ia
#if KEY_DOMDEC==1
       endif
#endif
       IF(IMOVE(I) == 0) THEN
          AMASSI=AMASS(I)
          XI=X(I)-XCM
          YI=Y(I)-YCM
          ZI=Z(I)-ZCM
          VXI=VX(I)-VXCM-OYCM*ZI+OZCM*YI
          VYI=VY(I)-VYCM-OZCM*XI+OXCM*ZI
          VZI=VZ(I)-VZCM-OXCM*YI+OYCM*XI
          VX(I)=VXI
          VY(I)=VYI
          VZ(I)=VZI
          DXI=DX(I)-AMASSI*(AXCM+QYCM*ZI-QZCM*YI)
          DYI=DY(I)-AMASSI*(AYCM+QZCM*XI-QXCM*ZI)
          DZI=DZ(I)-AMASSI*(AZCM+QXCM*YI-QYCM*XI)
          IF(I<ISGSTA .OR. I>ISGEND)THEN
            DX(I)=DXI
            DY(I)=DYI
            DZ(I)=DZI
            CYCLE
          ENDIF
          ti = INLCKP(i)
          j = i + 1
          if (QNHLANG(ti).and.(isDrude(i).or.isDrude(j)))then
             GAM=NHGAMMAD*DELTA
          else if(QNHLANG(ti))then
             GAM=NHGAMMA*DELTA
          endif
          GAMM=ONE+HALF*(GAM+SGSCAL)
          SGFTJ=SGFTI*SGWT(I)
          FACT1=(ONE+HALF*GAM)/GAMM
          FACT2=AMASSI*SGSCAL/GAMM/DELTA
          DXI=FACT1*DXI+FACT2*VXI
          DYI=FACT1*DYI+FACT2*VYI
          DZI=FACT1*DZI+FACT2*VZI
          DX(I)=DXI
          DY(I)=DYI
          DZ(I)=DZI
          FACT2=HALF*DELTA/GAMM/AMASSI
          VXT =FACT1*VXI-DXI*FACT2
          VYT =FACT1*VYI-DYI*FACT2
          VZT =FACT1*VZI-DZI*FACT2
          ! local average of velocities
          VXI=VXI-VXG-OYG*ZI+OZG*YI
          VYI=VYI-VYG-OZG*XI+OXG*ZI
          VZI=VZI-VZG-OXG*YI+OYG*XI
          GXI=SGVX(I)
          GYI=SGVY(I)
          GZI=SGVZ(I)
          EKIN=EKIN+AMASSI*(VXI*VXI+VYI*VYI+VZI*VZI)
          EKINSG=EKINSG+AMASSI*(GXI*GXI+GYI*GYI+GZI*GZI)
         IF(.NOT.QSGLDG)THEN
          ! Accumulation for SGLD parameter calculation
          FXI=SGFX(I)
          FYI=SGFY(I)
          FZI=SGFZ(I)
          DXI=SGDX(I)
          DYI=SGDY(I)
          DZI=SGDZ(I)
          FACT1=AMASSI/DELTA
          SGXI=FACT1*((SGSCAL)*VXT-SGFTJ*GAM*GXI)
          SGYI=FACT1*((SGSCAL)*VYT-SGFTJ*GAM*GYI)
          SGZI=FACT1*((SGSCAL)*VZT-SGFTJ*GAM*GZI)
          SGXK=SGXI+SGFDI*DXI+(SGFFI-SGFDI)*FXI
          SGYK=SGYI+SGFDI*DYI+(SGFFI-SGFDI)*FYI
          SGZK=SGZI+SGFDI*DZI+(SGFFI-SGFDI)*FZI
          SGXJ=SGXK-SGFDI*(DXI-FXI)
          SGYJ=SGYK-SGFDI*(DYI-FYI)
          SGZJ=SGZK-SGFDI*(DZI-FZI)
          ! save guiding forces for virial calculation
          SGDX(I)=SGXK
          SGDY(I)=SGYK
          SGDZ(I)=SGZK
          GX2I=SGAVG0*SGGX(I)+SGAVG1*SGXI
          GY2I=SGAVG0*SGGY(I)+SGAVG1*SGYI
          GZ2I=SGAVG0*SGGZ(I)+SGAVG1*SGZI
          SGGX(I)=GX2I
          SGGY(I)=GY2I
          SGGZ(I)=GZ2I
          GX2K=SGAVG0*SGHX(I)+SGAVG1*SGXK
          GY2K=SGAVG0*SGHY(I)+SGAVG1*SGYK
          GZ2K=SGAVG0*SGHZ(I)+SGAVG1*SGZK
          SGHX(I)=GX2K
          SGHY(I)=GY2K
          SGHZ(I)=GZ2K
          GX2J=SGAVG0*SGKX(I)+SGAVG1*FXI
          GY2J=SGAVG0*SGKY(I)+SGAVG1*FYI
          GZ2J=SGAVG0*SGKZ(I)+SGAVG1*FZI
          SGKX(I)=GX2J
          SGKY(I)=GY2J
          SGKZ(I)=GZ2J
          VXI=VXT-GXI
          VYI=VYT-GYI
          VZI=VZT-GZI
          DXI=DXI-FXI
          DYI=DYI-FYI
          DZI=DZI-FZI
          SUMFF=SUMFF+(FXI*FXI+FYI*FYI+FZI*FZI)
          SUMDD=SUMDD+(DXI*DXI+DYI*DYI+DZI*DZI)
          SUMGF=SUMGF+(FXI*GX2I+FYI*GY2I+FZI*GZ2I)
          SUMGD=SUMGD+(DXI*(SGXI-GX2I)+DYI*(SGYI-GY2I)+DZI*(SGZI-GZ2I))
          FACT1=GAM*AMASSI/DELTA
          FACT2=FACT1*FACT1
          SUMPP=SUMPP+FACT2*(GXI*GXI+GYI*GYI+GZI*GZI)
          SUMGP=SUMGP+FACT1*(GX2K*GXI+GY2K*GYI+GZ2K*GZI)
         ENDIF
       ENDIF
    ENDDO
#if KEY_PARALLEL==1
    CALL PSYNC()
    TIMMER=ECLOCK()
    GCARR(1)=SUMFF
    GCARR(2)=SUMDD
    GCARR(3)=SUMGF
    GCARR(4)=SUMGD
    GCARR(5)=SUMPP
    GCARR(6)=SUMGP
    GCARR(7)=VIRSG
    GCARR(8)=EKIN
    GCARR(9)=EKINSG
    GCARR(10)=EKINSGG
    CALL GCOMB(GCARR,10)
    SUMFF=GCARR(1)
    SUMDD=GCARR(2)
    SUMGF=GCARR(3)
    SUMGD=GCARR(4)
    SUMPP=GCARR(5)
    SUMGP=GCARR(6)
    VIRSG=GCARR(7)
    EKIN=GCARR(8)
    EKINSG=GCARR(9)
    EKINSGG=GCARR(10)
    TMERI(TIMGCOMM)=TMERI(TIMGCOMM)+ECLOCK()-TIMMER
    TIMMER=ECLOCK()
#endif
    SGSCAL=SGSCAL/TIMFAC/DELTA
    ! Estimate momentum guiding factor and guiding temperature
    TEMPLFI=SGAVP0*TEMPLFI+SGAVP1*EKINSG/(KBOLTZ*NDEGFSG)
    TEMPHFI=SGAVP0*TEMPHFI+SGAVP1*((EKIN-EKINSG)/(KBOLTZ*NDEGFSG))
    IF(QSGLDG)THEN
      FRCLF=ONE+SGFFI
      FRCHF=ONE
      SGCFLF=TSGSET/TEMPSG
      SGCFHF=TSGSET/(TSGSET-TEMPLFI)
    ! Estimate reference temperatures
      TREFLFI=TREFLF
      IF(TREFLF<RSMALL)THEN
        TREFLFI=TEMPLFI*SGCFLF
#if KEY_REPDSTR==1 /*rexsgld*/
        if(TRXLF>ZERO)TREFLFI=TRXLF
#endif /* (rexsgld)*/
      ENDIF
      IF(TREFHF<RSMALL)TREFHFI=TSGSET-TREFLFI
      TEMPSGI=SGAVP0*TEMPSGI+SGAVP1*TSGSET*TREFHFI*TEMPLFI/TREFLFI/TEMPHFI
      RETURN
    ENDIF
    AVGFF=SGAVP0*AVGFF+SGAVP1*SUMFF
    AVGDD=SGAVP0*AVGDD+SGAVP1*SUMDD
    AVGGF=SGAVP0*AVGGF+SGAVP1*SUMGF
    AVGGD=SGAVP0*AVGGD+SGAVP1*SUMGD
    AVGPP=SGAVP0*AVGPP+SGAVP1*SUMPP
    AVGGP=SGAVP0*AVGGP+SGAVP1*SUMGP
    ! Estimate SGLD factors
    SGEFLF=(AVGGF/AVGFF+ONE)
    SGEFHF=(AVGGD/(AVGDD)+ONE)
    SGCFLF=(AVGGP/AVGPP)+ONE
    ! Estimate reference temperatures
    TREFLFI=TREFLF
    IF(TREFLF<RSMALL)THEN
      TREFLFI=MAX(TEMPLFI*SGCFLF,TEMPLFI/TEN)
#if KEY_REPDSTR==1 /*rexsgld*/
      if(TRXLF>ZERO)TREFLFI=TRXLF
#endif /* (rexsgld)*/
    ENDIF
    IF(TREFHF<RSMALL)TREFHFI=TSGSET-TREFLFI
    TEMPSGI=SGAVP0*TEMPSGI+SGAVP1*TSGSET*TREFHFI*TEMPLFI/TREFLFI/TEMPHFI
    FRCLF=SGEFLF+SGFFI
    FRCHF=SGEFHF+SGFDI
    SGCFLF=TREFLFI/TEMPLFI
    SGCFHF=TREFHFI/TEMPHFI
    ! Adjust guiding factors if wanted
    IF(QCTSG)THEN
      SGFTI=SGFTI+SGAVP1*(TEMPSG-TEMPSGI)*ABS(TEMPSG-TSGSET)/TSGSET/(TSGSET+TEMPSG)
      SGFTI=MIN(SGFTI,TEN)
      SGFTI=MAX(SGFTI,-TEN)
    ENDIF
    IF(QSGBZ)THEN
    ! Estimate force guiding factor
      IF(QCTSG)THEN
        IF(SGFD*SGFD<RSMALL)SGFDI=SGAVP0*SGFDI+SGAVP1*((TSGSET-TEMPLF)/(TSGSET-TREFLFI)-SGEFHF)
        IF(SGFF*SGFF<RSMALL)SGFFI=SGAVP0*SGFFI+SGAVP1*(TEMPLF/TREFLFI-SGEFLF)
      ELSE
        IF(SGFD*SGFD<RSMALL)THEN
           IF(TREFHF>RSMALL)THEN
              FACT1=TEMPHFI
           ELSE
              FACT1=TSGSET-TEMPLFI
           ENDIF
           SGFDI=SGAVP0*SGFDI+SGAVP1*(FACT1/(TSGSET-TREFLFI)-SGEFHF)
        ENDIF
        IF(SGFF*SGFF<RSMALL)SGFFI=SGAVP0*SGFFI+SGAVP1*(TEMPLFI/TREFLFI-SGEFLF)
      ENDIF
    ENDIF
    RETURN
  END SUBROUTINE SGLDG

  SUBROUTINE SGMDG(ATFRST,ATLAST,DELTA,NDEGF,MAXNOS,INLCKP,RTMPR,SNHV, &
                   IMOVE,AMASS,X,Y,Z,VX,VY,VZ,DX,DY,DZ)
    !
    !     This routine calculates the guiding forces and adds them
    !     to the systematic forces.
    !         Assum SGFT*gamma=constant for all atoms
    !
    use chm_kinds
    use clcg_mod,only:random
    use dimens_fcm
    use number
    use stream
    use consta
    use parallel
    use machutil,only:eclock
    use rndnum
#if KEY_DOMDEC==1
    use domdec_common,only:q_domdec,natoml,atoml
#endif
    !
    integer MAXNOS            ! INPUT   Maximum value for NOBL
    INTEGER ATFRST,ATLAST,NDEGF
    real(chm_real)  DELTA
    INTEGER IMOVE(*)
    integer,dimension(:) :: INLCKP  ! INPUT
    real(chm_real)  SNHV(MAXNOS)      ! I/O     NH velocities (zeta)
    real(chm_real)  RTMPR(MAXNOS)     ! INPUT   NH target temperatures
    real(chm_real)  AMASS(*)
    real(chm_real)  X(*),Y(*),Z(*),VX(*),VY(*),VZ(*), &
         DX(*),DY(*),DZ(*)
    !
    INTEGER IG,TI
#if KEY_PARALLEL==1
    real(chm_real) GCARR(30),TIMMER
#endif
    !
    INTEGER I,J,K,IA,IFRST,ILAST
    real(chm_real) DELTA2,GAM
    real(chm_real) VXI,VYI,VZI,VXT,VYT,VZT,GXI,GYI,GZI,DXI,DYI,DZI
    real(chm_real) FXI,FYI,FZI,SGXI,SGYI,SGZI,SGXJ,SGYJ,SGZJ,SGXK,SGYK,SGZK
    real(chm_real) GX2I,GY2I,GZ2I,GX2J,GY2J,GZ2J,GX2K,GY2K,GZ2K
    real(chm_real) FACT0,FACT1,FACT2,AMASSI,SGFTJ
    real(chm_real) A,B,RDX,RDY,RDZ
    real(chm_real) SUMGV,SUMVV,SUMV2,SUMG2,EKIN,EKINSG,EKINSGG
    real(chm_real) SUMFF,SUMDD,SUMGF,SUMGD
    real(chm_real) SUMPP,SUMGP,SUMFP,SUMQQ,SUMHQ,SUMDQ
    real(chm_real) XI,YI,ZI,XCM,YCM,ZCM,XX,XY,XZ,YY,YZ,ZZ
    real(chm_real) FXT,FYT,FZT,PXT,PYT,PZT
    real(chm_real) LXT,LYT,LZT,TXT,TYT,TZT
    real(chm_real) VXCM,VYCM,VZCM,AXCM,AYCM,AZCM
    real(chm_real) OXCM,OYCM,OZCM,QXCM,QYCM,QZCM
    real(chm_real) TCM(3,3),U(3,3),SCR(24),AMOM(6)
    !
    DELTA2=DELTA*DELTA
    !
    IFRST=ATFRST
    ILAST=ATLAST
    IF(ISGSTA > IFRST) IFRST=ISGSTA
    IF(ISGEND < ILAST) ILAST=ISGEND
    !
    !WXW Calculate constraint factor to eliminate energy input from guiding force
    SUMGF=ZERO
    SUMGD=ZERO
    SUMGV=ZERO
    SUMVV=ZERO
    EKIN=ZERO
    EKINSG=ZERO
    EKINSGG=ZERO
    XCM=ZERO
    YCM=ZERO
    ZCM=ZERO
    XX=ZERO
    XY=ZERO
    XZ=ZERO
    YY=ZERO
    YZ=ZERO
    ZZ=ZERO
    PXT=ZERO
    PYT=ZERO
    PZT=ZERO
    LXT=ZERO
    LYT=ZERO
    LZT=ZERO
    FXT=ZERO
    FYT=ZERO
    FZT=ZERO
    TXT=ZERO
    TYT=ZERO
    TZT=ZERO
    K=0
    do ia = atfrst,atlast
#if KEY_DOMDEC==1
       if (q_domdec) then
          i = atoml(ia)
       else
#endif
          i = ia
#if KEY_DOMDEC==1
       endif
#endif
       IF(IMOVE(I) == 0) THEN
          AMASSI=AMASS(I)
          XI=X(I)
          YI=Y(I)
          ZI=Z(I)
          VXI=VX(I)
          VYI=VY(I)
          VZI=VZ(I)
          DXI=DX(I)
          DYI=DY(I)
          DZI=DZ(I)
          IF(I>=ISGSTA .AND. I<=ISGEND)THEN
            TI = INLCKP(I)
            FACT0=SNHV(TI)*AMASSI  ! instantaneous friction
            SGFTJ=SGFTI*SGWT(I)
            GXI=SGAVG0*SGVX(I)+SGAVG1*VXI
            GYI=SGAVG0*SGVY(I)+SGAVG1*VYI
            GZI=SGAVG0*SGVZ(I)+SGAVG1*VZI
            SGVX(I)=GXI
            SGVY(I)=GYI
            SGVZ(I)=GZI
            FXI=SGAVG0*SGFX(I)+SGAVG1*(SGDX(I)+SHKDX(I))
            FYI=SGAVG0*SGFY(I)+SGAVG1*(SGDY(I)+SHKDY(I))
            FZI=SGAVG0*SGFZ(I)+SGAVG1*(SGDZ(I)+SHKDZ(I))
            SGFX(I)=FXI
            SGFY(I)=FYI
            SGFZ(I)=FZI
            SGXI=SGAVG0*SGKX(I)+SGAVG1*GXI
            SGYI=SGAVG0*SGKY(I)+SGAVG1*GYI
            SGZI=SGAVG0*SGKZ(I)+SGAVG1*GZI
            SGKX(I)=SGXI
            SGKY(I)=SGYI
            SGKZ(I)=SGZI
            EKINSGG=EKINSGG+AMASSI*(SGXI*SGXI+SGYI*SGYI+SGZI*SGZI)
            !SGFTJ=(TFRIC-SGFT)*SGWT(I)
            SGFTJ=(SGFTI)*SGWT(I)
            SGXJ=SGAVG0*SGGX(I)+SGAVG1*FXI
            SGYJ=SGAVG0*SGGY(I)+SGAVG1*FYI
            SGZJ=SGAVG0*SGGZ(I)+SGAVG1*FZI
            SGGX(I)=SGXJ
            SGGY(I)=SGYJ
            SGGZ(I)=SGZJ
            SGXI=(SGDX(I)+SHKDX(I)-FXI)
            SGYI=(SGDY(I)+SHKDY(I)-FYI)
            SGZI=(SGDZ(I)+SHKDZ(I)-FZI)
            SGHX(I)=SGAVP0*SGHX(I)+SGAVP1*(SGXI*VXI+SGYI*VYI+SGZI*VZI)
            SGHY(I)=SGAVP0*SGHY(I)+SGAVP1*AMASSI*TIMFAC*(VXI*VXI+VYI*VYI+VZI*VZI)
            SGGAMMA(I)=SGAVP0*SGGAMMA(I)+  &
              SGAVP1*SGHX(I)/SGHY(I)
            FACT1=SGGAMMA(I)*TIMFAC*AMASSI
            SGXK=FXI-SGXJ
            SGYK=FYI-SGYJ
            SGZK=FZI-SGZJ
            SGXI=(SGFTJ*FACT1)*GXI-SGFFI*SGXK
            SGYI=(SGFTJ*FACT1)*GYI-SGFFI*SGYK
            SGZI=(SGFTJ*FACT1)*GZI-SGFFI*SGZK
            SGXJ=SGXI
            SGYJ=SGYI
            SGZJ=SGZI
            SGDX(I)=-SGXJ
            SGDY(I)=-SGYJ
            SGDZ(I)=-SGZJ
            !write(*,'("SGLDG:",I6,8F10.4)')I,SGWT(I),SGGAMMA(I),RDX,RDY,RDZ,SGXI,SGYI,SGZI
            DXI=DXI-SGXI
            DYI=DYI-SGYI
            DZI=DZI-SGZI
            DX(I)=DXI
            DY(I)=DYI
            DZ(I)=DZI
            !IF(AMASSI>RSMALL)THEN
            !  FACT2=DELTA/TWO/AMASSI
            !ELSE
            !  FACT2=ZERO
            !ENDIF
            !VXT =VXI-DXI*FACT2
            !VYT =VYI-DYI*FACT2
            !VZT =VZI-DZI*FACT2
            VXT =VXI
            VYT =VYI
            VZT =VZI
            SUMGV=SUMGV+(SGXJ*VXT+SGYJ*VYT+SGZJ*VZT)
            SUMVV=SUMVV+AMASSI*(VXT*VXT+VYT*VYT+VZT*VZT)
            EKIN=EKIN+AMASSI*(VXI*VXI+VYI*VYI+VZI*VZI)
            EKINSG=EKINSG+AMASSI*(GXI*GXI+GYI*GYI+GZI*GZI)
          ENDIF
          XCM=XCM+XI*AMASSI
          YCM=YCM+YI*AMASSI
          ZCM=ZCM+ZI*AMASSI
          XX=XX+XI*XI*AMASSI
          XY=XY+XI*YI*AMASSI
          XZ=XZ+XI*ZI*AMASSI
          YY=YY+YI*YI*AMASSI
          YZ=YZ+YI*ZI*AMASSI
          ZZ=ZZ+ZI*ZI*AMASSI
          PXT=PXT+AMASSI*VXI
          PYT=PYT+AMASSI*VYI
          PZT=PZT+AMASSI*VZI
          LXT=LXT+(YI*VZI-ZI*VYI)*AMASSI
          LYT=LYT+(ZI*VXI-XI*VZI)*AMASSI
          LZT=LZT+(XI*VYI-YI*VXI)*AMASSI
          FXT=FXT+DXI
          FYT=FYT+DYI
          FZT=FZT+DZI
          TXT=TXT+(YI*DZI-ZI*DYI)
          TYT=TYT+(ZI*DXI-XI*DZI)
          TZT=TZT+(XI*DYI-YI*DXI)
       ENDIF
    ENDDO
    !write(*,'("SGLDG GAMMA: ",8E15.6)')SGFTI,SGFFI,SUM(SGGAMMA)/SUM(SGWT),SUMGD,SUMGF
    !write(*,'("SGLDG GAMMA: ",4E15.6)')SUM(SGGAMMA),SUM(SGHX),SUM(SGHY),SUM(SGGAMMA)/SUM(SGWT)
    !
#if KEY_PARALLEL==1
    CALL PSYNC()
    TIMMER=ECLOCK()
    GCARR(1)=SUMGV
    GCARR(2)=SUMVV
    GCARR(3)=EKIN
    GCARR(4)=EKINSG
    GCARR(5) = XCM
    GCARR(6) = YCM
    GCARR(7) = ZCM
    GCARR(8) = XX
    GCARR(9) = XY
    GCARR(10) = XZ
    GCARR(11) = YY
    GCARR(12) = YZ
    GCARR(13) = ZZ
    GCARR(14)=PXT
    GCARR(15)=PYT
    GCARR(16)=PZT
    GCARR(17)=LXT
    GCARR(18)=LYT
    GCARR(19)=LZT
    GCARR(20)=FXT
    GCARR(21)=FYT
    GCARR(22)=FZT
    GCARR(23) = TXT
    GCARR(24) = TYT
    GCARR(25) = TZT
    GCARR(26)=EKINSGG
    CALL GCOMB(GCARR,26)
    SUMGV=GCARR(1)
    SUMVV=GCARR(2)
    EKIN=GCARR(3)
    EKINSG=GCARR(4)
    XCM=GCARR(5)
    YCM=GCARR(6)
    ZCM=GCARR(7)
    XX=GCARR(8)
    XY=GCARR(9)
    XZ=GCARR(10)
    YY=GCARR(11)
    YZ=GCARR(12)
    ZZ=GCARR(13)
    PXT=GCARR(14)
    PYT=GCARR(15)
    PZT=GCARR(16)
    LXT=GCARR(17)
    LYT=GCARR(18)
    LZT=GCARR(19)
    FXT=GCARR(20)
    FYT=GCARR(21)
    FZT=GCARR(22)
    TXT=GCARR(23)
    TYT=GCARR(24)
    TZT=GCARR(25)
    EKINSGG=GCARR(26)
    TMERI(TIMGCOMM)=TMERI(TIMGCOMM)+ECLOCK()-TIMMER
    TIMMER=ECLOCK()
#endif
    !
    ! Translation
    IF(QSGSTOPT)THEN
      !  net velocity
      VXCM=PXT/TOTMASS
      VYCM=PYT/TOTMASS
      VZCM=PZT/TOTMASS
      !  net accelation
      AXCM=FXT/TOTMASS
      AYCM=FYT/TOTMASS
      AZCM=FZT/TOTMASS
    ELSE
      VXCM=ZERO
      VYCM=ZERO
      VZCM=ZERO
      AXCM=ZERO
      AYCM=ZERO
      AZCM=ZERO
    ENDIF
    ! Rotation
    IF(QSGSTOPR)THEN
      XCM=XCM/TOTMASS
      YCM=YCM/TOTMASS
      ZCM=ZCM/TOTMASS
      XX=XX-XCM*XCM*TOTMASS
      XY=XY-XCM*YCM*TOTMASS
      XZ=XZ-XCM*ZCM*TOTMASS
      YY=YY-YCM*YCM*TOTMASS
      YZ=YZ-YCM*ZCM*TOTMASS
      ZZ=ZZ-ZCM*ZCM*TOTMASS
      LXT=LXT-(YCM*PZT-ZCM*PYT)
      LYT=LYT-(ZCM*PXT-XCM*PZT)
      LZT=LZT-(XCM*PYT-YCM*PXT)
      TXT=TXT-(YCM*FZT-ZCM*FYT)
      TYT=TYT-(ZCM*FXT-XCM*FZT)
      TZT=TZT-(XCM*FYT-YCM*FXT)
    !
      AMOM(1)=YY+ZZ
      AMOM(2)=-XY
      AMOM(3)=-XZ
      AMOM(4)=XX+ZZ
      AMOM(5)=-YZ
      AMOM(6)=XX+YY
    !
      CALL DIAGQ(3,3,AMOM,U,SCR(4),SCR(7),SCR(10),SCR(13),SCR(1),  &
         SCR(16),SCR(19),SCR(22),0)
    !
      DO I=1,3
         IF(ABS(SCR(I)) < TENM5) SCR(I)=SIGN(TENM5,SCR(I))
         SCR(I)=ONE/SCR(I)
      ENDDO
    !
      DO I=1,3
         DO J=1,3
            TCM(I,J)=0.0
            DO K=1,3
               TCM(I,J)=TCM(I,J)+U(I,K)*U(J,K)*SCR(K)
            ENDDO
         ENDDO
      ENDDO
    ! angular velocity of COM
      OXCM=LXT*TCM(1,1)+LYT*TCM(1,2)+LZT*TCM(1,3)
      OYCM=LXT*TCM(2,1)+LYT*TCM(2,2)+LZT*TCM(2,3)
      OZCM=LXT*TCM(3,1)+LYT*TCM(3,2)+LZT*TCM(3,3)
    ! angular acceleration of COM
      QXCM=TXT*TCM(1,1)+TYT*TCM(1,2)+TZT*TCM(1,3)
      QYCM=TXT*TCM(2,1)+TYT*TCM(2,2)+TZT*TCM(2,3)
      QZCM=TXT*TCM(3,1)+TYT*TCM(3,2)+TZT*TCM(3,3)
    ELSE
      OXCM=ZERO
      OYCM=ZERO
      OZCM=ZERO
      QXCM=ZERO
      QYCM=ZERO
      QZCM=ZERO
    ENDIF
    !
    !SGSCAL=DELTA*SUMGV/(SUMVV-HALF*DELTA*SUMGV)
    SGSCAL=DELTA*SUMGV/(SUMVV)
    !SGSCAL=ZERO
    !
    FACT0=TWO/(TWO+SGSCAL)
    !
    SUMFF=ZERO
    SUMDD=ZERO
    SUMGF=ZERO
    SUMGD=ZERO
    SUMPP=ZERO
    SUMGP=ZERO
    SUMQQ=ZERO
    do ia = atfrst,atlast
#if KEY_DOMDEC==1
       if (q_domdec) then
          i = atoml(ia)
       else
#endif
          i = ia
#if KEY_DOMDEC==1
       endif
#endif
       IF(IMOVE(I) == 0) THEN
          AMASSI=AMASS(I)
          XI=X(I)-XCM
          YI=Y(I)-YCM
          ZI=Z(I)-ZCM
          VXI=VX(I)-VXCM-OYCM*ZI+OZCM*YI
          VYI=VY(I)-VYCM-OZCM*XI+OXCM*ZI
          VZI=VZ(I)-VZCM-OXCM*YI+OYCM*XI
          VX(I)=VXI
          VY(I)=VYI
          VZ(I)=VZI
          DXI=DX(I)-AMASSI*(AXCM+QYCM*ZI-QZCM*YI)
          DYI=DY(I)-AMASSI*(AYCM+QZCM*XI-QXCM*ZI)
          DZI=DZ(I)-AMASSI*(AZCM+QXCM*YI-QYCM*XI)
          IF(I<ISGSTA .OR. I>ISGEND)THEN
            DX(I)=DXI
            DY(I)=DYI
            DZ(I)=DZI
            CYCLE
          ENDIF
          FACT1=SGSCAL*AMASSI/DELTA
          DXI=DXI+FACT1*VXI
          DYI=DYI+FACT1*VYI
          DZI=DZI+FACT1*VZI
          DX(I)=DXI
          DY(I)=DYI
          DZ(I)=DZI
       ENDIF
    ENDDO
    ! Estimate momentum guiding factor and guiding temperature
    TEMPLFI=SGAVP0*TEMPLFI+SGAVP1*EKINSG/(KBOLTZ*NDEGFSG)
    TEMPHFI=SGAVP0*TEMPHFI+SGAVP1*((EKIN-EKINSG)/(KBOLTZ*NDEGFSG))
      FRCLF=ONE+SGFFI
      FRCHF=ONE
      SGCFLF=TSGSET/TEMPSG
      SGCFHF=TSGSET/(TSGSET-TEMPLFI)
    ! Estimate reference temperatures
      TREFLFI=TREFLF
      IF(TREFLF<RSMALL)THEN
        TREFLFI=TEMPLFI*SGCFLF
#if KEY_REPDSTR==1 /*rexsgld*/
        if(TRXLF>ZERO)TREFLFI=TRXLF
#endif /* (rexsgld)*/
      ENDIF
      IF(TREFHF<RSMALL)TREFHFI=TSGSET-TREFLFI
      TEMPSGI=SGAVP0*TEMPSGI+SGAVP1*TSGSET*TREFHFI*TEMPLFI/TREFLFI/TEMPHFI
    RETURN
  END SUBROUTINE SGMDG

#endif /* (sgld_main)*/
  SUBROUTINE NULL_SGLD
    RETURN
  END SUBROUTINE NULL_SGLD
end module sgld
