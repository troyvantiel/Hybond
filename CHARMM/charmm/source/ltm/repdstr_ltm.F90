module repdstr
  use chm_kinds
  use dimens_fcm
  use parallel, only: maxnode
  implicit none

#if KEY_REPDSTR==0 && KEY_REPDSTR2 == 0 /*repdstr_fcm*/
      logical, parameter :: QREPDSTR = .false.
#else /* (repdstr_fcm)*/
#if KEY_PARALLEL==1 /*check_dependency*/
!-----------------------------------------------------------------------
!
!     QREPDSTR         are distributed replicas in use at present ?
!     QRDQTT           flag for use of RDCMND (true=no global broadcast)
!     QWRQTT           flag for use of OUTU (true=every mynod=0 writes)
!     NREPDSTR         number of distributed replicas in use at present ?
!     MAXREPDSTR       maximum number of distributed replicas
!     QREXCHG          flag if this is replica exchange run
!     QREXCHGL         flag if this is lambda replica exchange run for free
!                      energy purturbation or unbrella sampling
!     QEXPT            flag if this is parallel tempering replica exchange
!                      for internal barrier of protein side chain or 
!                      backbone(?)
!     QEXBK            flag if replica exchange blocker is applied.
!     IRBK             IREPDSTR of replica exchange blocker   
!
!-----------------------------------------------------------------------
      LOGICAL QREPDSTR,QRDQTT,QWRQTT,QREXCHG,QREXCHGL,QEXPT,QEX2D,QEXBK
      logical :: qrepioset, qrepdnames = .true.
      INTEGER NREPDSTR,IREPDSTR,IREXFQ
      INTEGER ITDRB,REPSEED,IOLORIG,PRNLORIG,WRNLORIG,IOON,IUNREX,NATREPCMD
      integer NREPT,NREPX,IRBK,EWRITU
      integer,parameter  :: NWREX=1000
      integer,parameter  :: MAXREPDSTR=64

      real(chm_real) TEMPRX(MAXNODE)

#endif /* (check_dependency)*/

contains
  subroutine repdstr_iniall
    ioon=0
    qrdqtt=.false.
    qwrqtt=.false.
    qrepioset=.false.
    return
  end subroutine repdstr_iniall

#endif /* (repdstr_fcm)*/

end module repdstr

