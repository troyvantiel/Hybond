module rtf
  use chm_kinds
  use dimens_fcm
  use parallel,only:mynod
  implicit none
  !-----------------------------------------------------------------------
  !     The Residue Topology File (RTF)
  !
  !     Purpose: To store information for generating macromolecular energy
  !     expressions on a monomer unit basis. Note: a very good description
  !     of the various fields is available in IO.DOC (get there via INFO).
  !
  !     I/O: Formatted and unformatted read: RTFRDR in io/rtfio.src
  !
  !
  !     Variable  Index  Purpose
  !
  !     NRTRS            Number of residues in RTF
  !     RTFTYP           USED TO CODE AUTO GENERATION FLAGS
  !     RTCNTL           Control array read from the file
  !     NATCT            Number of atom type codes in RTF
  !     DEFF             Name of default first patch residue
  !     DEFL             Name of default last patch residue
  !     AUTOT            Logical flag to request auto angle generation
  !     AUTOD            Logical flag to request auto dihedral generation
  !     AUTOP            Logical flag to request PATCh command auto generation
  !
  !     ARMASS    A-type Mass of atoms indexed by atom type code
  !     ATCT      A-type Names of parameter types. Corresponds to ARMASS
  !
  !     NIC       Res.   Number of atoms, bonds, angles, phi's, imphi's,
  !                      exclusions, hydrogen bond donors, acceptors,
  !                      builder internal coordinates, and groups;
  !                      respectively.
  !     RTRTYP    Res.   Type of the residue (normal or patch)
  !     AA        Res.   Name of the residue (6 characters)
  !     AAF       Res.   Name of first patch residue
  !     AAL       Res.   Name of last patch residue
  !
  !     CHG       Atom   Charge of each atom
  !     ALPH      Atom   Atomic polarizability of each atom (with Drude method)
  !     THOL      Atom   Thole screening parameter of each atom (with Drude method)
  !     FTP       Atom   IUPAC name of each atom
  !     MAC       Atom   Atom type code (index into ATCT)
  !     MXN       Atom   Number of non-bonded exclusions for each atom
  !     GRPR      Atom   Group number of atom
  !     DELAT     Atom   true if atom to be deleted (in patch residues)
  !
  !     MIB       Bond   First atom of bond
  !     MJB       Bond   Second atom of bond
  !     DELBD     Bond   true if bond to be deleted (in patch residues)
  !
  !     MIT       Angle  First atom of bond angle
  !     MJT       Angle  Second atom of bond angle
  !     MKT       Angle  Third atom of bond angle
  !     DELAN     Angle  true if angle to be deleted (in patch residues)
  !
  !
  !     MIP       Phi    First atom of torsion
  !     MJP       Phi    Second atom of torsion
  !     MKP       Phi    Third atom of torsion
  !     MLP       Phi    Fourth atom of torsion
  !     DELPT     Phi    true if torsion to be deleted (in patch residues)
  !
  !     MIM       Imphi  First atom of improper torsion
  !     MJM       Imphi  Second atom of improper torsi
  !     MKM       Imphi  Third atom of improper torsion
  !     MLM       Imphi  Fourth atom of improper torsion
  !     DELMT     Imphi  true if improper torsion to be deleted (in patch)
  !
#if KEY_CMAP==1
  !     MI1CT     Cmap   First atom of first torsion
  !     MJ1CT     Cmap   Second atom of first torsion
  !     MK1CT     Cmap   Third atom of first torsion
  !     ML1CT     Cmap   Fourth atom of first torsion
  !     MI2CT     Cmap   First atom of second torsion
  !     MJ2CT     Cmap   Second atom of second torsion
  !     MK2CT     Cmap   Third atom of second torsion
  !     ML2CT     Cmap   Fourth atom of second torsion
  !     DELPCT    Cmap   true if cross-term map is to be deleted
  !
  !     LONE PAIR
  !     MLP0CT           type of LP (geometric information)
  !     MLP1CT           first atom
  !     MLP2CT           second atom
  !     MLP3CT           third atom
  !     MLP4CT           fourth atom
  !     MLPCT1D          1-dimensional array of host of lp center
  !     DELLPCT          lone-pairs true if term is to be deleted
  !     RLPCT            distance
  !     TLPCT            angle
  !     PLPCT            dihedral
  !
  !     Anisotropic polarizability for drudes
  !     MIANIS           first atom
  !     MJANIS           second atom
  !     MKANIS           third atom
  !     MLANIS           fourth atom
  !     DELANIS          true if anistropic term is to be deleted
  !     A11ANIS          principal axis 11 of the tensor
  !     A22ANIS          principal axis 22 of the tensor

#endif
  !
  !
  !     MNB       NBex   Non-bonded exclusions for the atoms (index with MXN)
  !
  !     MA        Acc.   Hydrogen bond acceptor
  !     MAA       Acc.   First antecedent to the hydrogen bond acceptor
  !     MAB       Acc.   !NOT USED ANYMORE BUT STILL PRESENT
  !     DELAC     Acc.   true if acceptor to be deleted (in patch residues)
  !
  !     MD        Donor  Hydrogen bond heavy atom donor
  !     MDA       Donor  !NOT USED ANYMORE BUT STILL PRESENT
  !     MDB       Donor  !NOT USED ANYMORE BUT STILL PRESENT
  !     MH        Donor  Hydrogen of hydrogen bond donor
  !     DELDN     Donor  true if donor to be deleted (in patch residues)
  !
  !     BARI      ic     First atom of internal coordinate.
  !     BARJ      ic     Second atom of internal coordinate.
  !     BARK      ic     Third atom of internal coordinate.
  !     BARL      ic     Fourth atom of internal coordinate.
  !     BART      ic     Logical flag specifying an improper type
  !     ICPHI     ic     Value of torsion angle for internal coordinate
  !     ICB1      ic     First bond length for internal coordinate
  !     ICB2      ic     Second bond length for internal coordinate
  !     ICTH1     ic     First bond angle for internal coordinate
  !     ICTH2     ic     Second bond angle for internal coordinate
  !     DELIC     ic     true if internal coordinate to be deleted (in patch)
  !
  !     dim_NICM             Number of different types of information stored
  !                      in the RTF
  !
  !      MAXIMUMS IN USE AS OF 5/21/82 IN TOPH9
  !     314,312,462,160,152,42,47,55,244,149
  !

  integer,parameter :: dim_nicm=13,initial_dim=100 !,initial_dim2=1000
  integer,dimension(dim_nicm) :: dim_rt
  integer,parameter ::  &
       ind_at  = 1 , &
       ind_bo  = 2 , &
       ind_th  = 3 , &
       ind_ph  = 4 , &
       ind_im  = 5 , &
       ind_xx  = 6 , &
       ind_hd  = 7 , &
       ind_ha  = 8 , &
       ind_ic  = 9 , &
       ind_cmap= 11, &
       ind_lp  = 12, &
       ind_anis= 13

  integer :: dim_natct=100
  integer :: dim_rtrs=100
  real(chm_real) :: factor_rt=1.3


  logical ucase

  character(len=6),allocatable,dimension(:)   :: aa, aaf, aal
  integer,         allocatable,dimension(:,:) :: nic
  integer,         allocatable,dimension(:)   :: rtrtyp

  character(len=6),allocatable,dimension(:) :: ftp
  integer,allocatable,dimension(:)   :: mac,mxn,grpr
#if KEY_MMFF==1
  integer,allocatable,dimension(:)   :: atnumr
#endif
  real(chm_real),allocatable,dimension(:)   :: chg,alph,thol
#if KEY_CHEQ==1
  real(chm_real),allocatable,dimension(:)   :: echh,ehah
#endif
  logical,allocatable,dimension(:)          :: delat

  character(len=6),allocatable,dimension(:) :: mib,mjb
#if KEY_MMFF==1
  integer,allocatable,dimension(:)   :: mbtype
#endif
  logical,allocatable,dimension(:)          :: delbd

  character(len=6),allocatable,dimension(:) :: mit
  character(len=6),allocatable,dimension(:) :: mjt
  character(len=6),allocatable,dimension(:) :: mkt
  logical,allocatable,dimension(:)          :: delan

  character(len=6),allocatable,dimension(:) :: mip
  character(len=6),allocatable,dimension(:) :: mjp
  character(len=6),allocatable,dimension(:) :: mkp
  character(len=6),allocatable,dimension(:) :: mlp
  logical,         allocatable,dimension(:) :: delpt

  character(len=6),allocatable,dimension(:) :: mim
  character(len=6),allocatable,dimension(:) :: mjm
  character(len=6),allocatable,dimension(:) :: mkm
  character(len=6),allocatable,dimension(:) :: mlm
  logical,         allocatable,dimension(:) :: delmt

  character(len=6),allocatable,dimension(:) :: bari
  character(len=6),allocatable,dimension(:) :: barj
  character(len=6),allocatable,dimension(:) :: bark
  character(len=6),allocatable,dimension(:) :: barl
  real(chm_real),  allocatable,dimension(:) :: icb1,icb2,icth1,icth2,icphi
  logical,         allocatable,dimension(:) :: delic,bart

  character(len=6),allocatable,dimension(:) :: mnb

  character(len=6),allocatable,dimension(:) :: ma
  character(len=6),allocatable,dimension(:) :: maa
  character(len=6),allocatable,dimension(:) :: mab
  logical,         allocatable,dimension(:) :: delac

  character(len=6),allocatable,dimension(:) :: md
  character(len=6),allocatable,dimension(:) :: mh
  character(len=6),allocatable,dimension(:) :: mda
  character(len=6),allocatable,dimension(:) :: mdb
  logical,         allocatable,dimension(:) :: deldn
  character(len=6),allocatable,dimension(:) :: atct

  character(len=6) :: deff, defl

  character(len=4),allocatable,dimension(:) :: mlp0ct
  character(len=6),allocatable,dimension(:) :: mlp1ct,mlp2ct,mlp3ct,mlp4ct, &
       mianis, mjanis, mkanis, mlanis
  !character(len=6),allocatable,dimension(:,:) :: mlpct
  character(len=6),allocatable,dimension(:) :: mlpct1d 
  integer :: maxcent_hosts=100  

#if KEY_CMAP==1
  character(len=6),allocatable,dimension(:) :: mi1ct_test, &
       mi1ct,mj1ct,mk1ct,ml1ct, &
       mi2ct,mj2ct,mk2ct,ml2ct
#endif

#if KEY_MMFF==1
  !     Modified by Jay Banks 12 Oct 95: added bond types (SINGLE, DOUBLE,
  !     and atomic numbers for MMFF.  Also (13 Oct) included "bond.fcm"
  !     from Merck/MSI version
  !
  ! from source/fcm/bond.fcm 1.3  (MSI)
  integer dummy,single,double,triple,aromatic
  parameter(dummy=0,single=1,double=2,triple=3,aromatic=4)
  integer,allocatable,dimension(:) :: atnumt ! atomic numbers for atom types
#endif

  real(chm_real),allocatable,dimension(:) :: armass

  integer nrtrs,rtftyp, rtcntl(20), natct, icenthst
  integer, allocatable,dimension(:):: ncenthst
  logical autot, autod, autop

#if KEY_CMAP==1
  logical,allocatable,dimension(:) ::  delpct,dellpct
#endif

  ! Lone Pairs
  real(chm_real),allocatable,dimension(:) :: rlpct, tlpct, plpct

  ! Anisotropic polarizability for drudes
  real(chm_real),allocatable,dimension(:) :: a11anis, a22anis

  logical,allocatable,dimension(:) ::  delanis

#if KEY_MMFF==1
  character(len=4),allocatable,dimension(:) :: ctype   !??? symbolic atom types
  integer,allocatable,dimension(:) ::    intdef  !??? integer atom types
#endif

contains

  subroutine rtf_iniall()
    nrtrs=0
    natct=0
    defl='NONE'
    deff='NONE'
    autot=.false.
    autod=.false.
    autop=.false.

    return
  end subroutine rtf_iniall

  subroutine rtfrdr(unit,title,ntitl,icard,qprint,lappe)
    !-----------------------------------------------------------------------
    !     THIS ROUTINE READS INFORMATION ABOUT THE INTERNAL COORDINATES
    !     AND THE GEOMETRY OF RESIDUES
    !     First symbolic input version: Robert Bruccoleri
    !     Last overhaul -- 3/82 -- BRB
    !     MODIFICATIONS FOR APPEND AND DELETE (IN PATCHES) -- 1/83 -- AB
    !
    !     The indices into the NIC arrays are as follows:
    !     IAT = 1, IBO = 2, ITH = 3, IPH = 4, IIM = 5,
    !     IX = 6,  ID  = 7, IA  = 8, IBL=9, IGRP=10
    !
    !     Modified by Rick Lapp 4 Aug 99
    !       Allow lower case atom type names for CFF.
#if KEY_MMFF==1
    !     MODIFIED BY JAY BANKS 12 OCT 95:
    !       Add atomic numbers and bond-order specification for MMFF.
#endif
    !

    use exfunc,only:srchws
    use consta,only:ccelec
    use number,only:zero,one
#if KEY_CFF==1
    use cff_fcm
#endif
#if KEY_MMFF==1
    use mmffm
#endif
#if KEY_CFF==1 || KEY_MMFF==1
    use ffieldm
#endif
    use comand
    use stream
    use string
    use machutil,only:die
    !
#if KEY_STRINGM==1 /*  VO stringm */
    use machio, only: ifreeu
    use multicom_aux
#endif
    !
    implicit none

    integer, parameter :: mxcms2=mxcmsz
    character(len=mxcms2) comly2
    integer comle2

    integer unit,ntitl,icard
    character(len=*) title(*)
    character(len=6) patf,patl
    character(len=4) hdr
    logical qprint,prnflg,qdelet,lappe,qfirst,newres

    integer j2,idel
    logical ok,eof,addok, ldrude
    logical, allocatable, dimension(:) :: wok
    real(chm_real)  tholea
    real(chm_real)  totchg
    integer i,j,ipatch,inx,ind,iatc,iat,ires
    integer nrx,nat,nbo,nth,nph,nim,nd,na,ngr,nbl
    integer nrxpt,natpt,nbopt,nthpt,nphpt,nimpt,ndpt,napt,nblpt
    integer errcnt,wrncnt,wrncto
    character(len=4) wrd
    character(len=6) w1,w2,w3,w4,wtemp
#if KEY_CMAP==1
    integer npct,  nict
    character(len=4) w5
    character(len=6) w6,w7,w8
    character(len=6), dimension(maxcent_hosts):: w
!   Dimension of maxcent_hosts = 100 is provided to w because number of lonepair host
!   is not known apriori. Thus provided a large number
    integer nlpct, nilpt
    integer nanisct, nanis
#endif
    integer     iw1,iw4
    integer start,stop,rnicm,nrtrs0
#if KEY_CHEQ==1
    real(chm_real) tech,teca
    real(chm_real) :: factkcalr
#endif
    integer, parameter :: mxtabl=50
    character(len=6) tblkey(mxtabl),keytab
    integer ntabl

    character(len=6) typenm
    integer typecd
    character(len=6) name
    real(chm_real) DIST, ANGL, DIHE, SCALE
    real(chm_real) a11, a22

#if KEY_MMFF==1
    integer  AtomicNumber
    external AtomicNumber
#endif
    !
#if KEY_STRINGM==1 /*  VO stringm v */
    integer :: oldiol
    logical :: qstr
    common /replicaio/ qstr ! need global variable
#endif
    character(len=4) :: srtf='RTF '
    character(len=4) :: blank='    '




#if KEY_CHEQ==1
    factkcalr=one/ccelec
#endif
    ldrude = .false.

    errcnt=0
    wrncnt=0
    wrncto=0
    if(.not.lappe) then
       autot=.false.
       autod=.false.
       autop=.false.
    endif
    !
#if KEY_STRINGM==1 /*  VO stringm */
    qstr=unit.le.size(ifreeu) ! protect against OOB on ifreeu
    if (qstr) qstr=((MOD(IFREEU(UNIT),8).eq.0).and.(IFREEU(UNIT).ne.0).and.(ME_LOCAL.eq.0))
    if (qstr) then
     oldiol=iolev
     iolev=1
    endif
#endif /* VO stringm */
    !
    if(iolev <= 0) then
#if KEY_PARALLEL==1
       call  rtf_parallel_bcast
#endif
       return
    end if

    formattedfile: if (icard /= 0) then
       call tryoro(unit,'FORMATTED')
       call rdtitl(title,ntitl,unit,0)

       call get_rtf_version

       if(.not.lappe) then
          if(allocated(atct))atct(1:dim_natct)=blank
          if(allocated(armass))armass(1:dim_natct)=zero
       endif
       if (rtcntl(1) == 17) call wrndie(-1,'<RTFRDR>', &
            'This version of the rtf is no longer supported')
       !           try to read it anyway

       if (rtcntl(1) >= 17) then
          !           begin procedure read-free-format-version-18-topology-file
          !
          ntabl=1
          tblkey(1)=blank
          qfirst=.true.
          if (lappe) then
             natpt=nic(ind_at,nrtrs)
             nbopt=nic(ind_bo,nrtrs)
             nthpt=nic(ind_th,nrtrs)
             nphpt=nic(ind_ph,nrtrs)
             nimpt=nic(ind_im,nrtrs)
#if KEY_CMAP==1
             npct=nic(ind_cmap,nrtrs)
             nlpct=nic(ind_lp,nrtrs)
             nanisct=nic(ind_anis,nrtrs)
#endif
             nrxpt=nic(ind_xx,nrtrs)
             ndpt=nic (ind_hd,nrtrs)
             napt=nic (ind_ha,nrtrs)
             nblpt=nic(ind_ic,nrtrs)
          else
             defl='NONE'
             deff='NONE'
             autot=.false.
             autod=.false.
             autop=.false.
             natpt=0
             nbopt=0
             nthpt=0
             nphpt=0
             nimpt=0
#if KEY_CMAP==1
             npct=0
             nlpct=0
             nanisct=0
             nanis=0
#endif
             nrxpt=0
             napt=0
             ndpt=0
             nblpt=0
             !
             natct=0
             !# when not appending, set initial valued of residue counter to
             !             zero
             nrtrs=0
          endif
#if KEY_CFF==1
          !       For CFF forcefields turn off the automatic upper case conversion
          !       done by PARSE.
          if (ffield == cff) ucase=.false.
#endif
          call reallocate_rtrs(nrtrs)

          readingfile: do while(.not.eof)
             if(errcnt+wrncnt /= wrncto) then
                wrncto=errcnt+wrncnt
                if (.not.(prnflg)) then
                   call prntst(outu,comly2,comle2,10,79)
                endif
             endif
             call rdcmnd(comlyn,mxcmsz,comlen,unit,eof,.false., &
                  prnflg,'RTF>    ')

             if (.not.(eof)) then
                ipatch=0
                call copyst(comly2,mxcms2,comle2,comlyn,comlen)
                qdelet=(indxa(comlyn,comlen,'DELE') > 0)
                !-MH: This is not optimal, but nrtrs should not be zero this time, but later
                !-MH: it is tested as zero ???
                !MFC 12-2012 has this ever been addressed?

                nrtrs0=nrtrs
                if(nrtrs0 == 0)nrtrs0=1
                if(rtrtyp(nrtrs0) /= 1.and.qdelet) then
                   call wrndie(0,'<RTFIO>', &
                        'DELEte only for Patch RESidues allowed. Will be ' &
                        //'ignored')
                   qdelet=.false.
                endif
                wrd=nexta4(comlyn,comlen)
                notblank: if (.not.(wrd == blank)) then
#if KEY_CFF==1
                   !     For CFF forcefields take care of most upper case conversion here.
                   !     Keywords ATOM and MASS need to get lower case atom type names.
                   if (.not.ucase) then
                      call cnvtuc(wrd,4)
                      if (wrd /= 'ATOM'.and.wrd /= 'MASS') &
                           call cnvtuc(comlyn,comlen)
                   endif
#endif
                   select case(wrd)
                   case( 'MASS' )
                      call rtfrdr_process_mass
                   case ( 'DECL')
                      call rtfrdr_process_declare
                   case ( 'RESI')
                      call rtfrdr_process_residue
                   case ( 'ATOM')
                      call rtfrdr_process_atoms
                   case ( 'BOND', 'SING', 'DOUB', 'TRIP', 'AROM' )
                      call rtfrdr_process_bond

                   case ( 'ANGL', 'THET')
                      call rtfrdr_process_angle
                   case ( 'DIHE', 'PHI ')
                      call rtfrdr_process_dihedral
                   case ( 'IMPH', 'IMPR')
                      call rtfrdr_process_improper
                   case ( 'ANIS')
                      call rtfrdr_process_anisotropic
                   case ( 'LONE')
                      call rtfrdr_process_lonepairs
#if KEY_CMAP==1
                   case ( 'CMAP')
                      call rtfrdr_process_cmap
#endif
                   case ( 'DONO')
                      call rtfrdr_process_donor
                   case ( 'ACCE')
                      call rtfrdr_process_acceptor
                   case ( 'IC  ', 'BUIL', 'BILD')
                      call rtfrdr_process_ic
                   case ( 'GROU')
                      call rtfrdr_process_group
                   case ( 'END ')
                      eof=.true.
                   case ( 'PRIN')
                      call rtfrdr_process_print
                   case ( 'DEFA')
                      call rtfrdr_process_default_patch
                   case ( 'AUTO')
                      call rtfrdr_process_auto_command
                   case ( 'PATC')
                      call rtfrdr_process_patch
                   case ( 'PRES')
                      call rtfrdr_process_pres
                   case ( '21  ')
                      call wrndie(2,'<RTFRDR>','Version number 21 not ' &
                           //'needed.')
                   case ( '22  ')
                      call wrndie(2,'<RTFRDR>','Version number 22 not ' &
                           //'needed.')
                   case default
                      call wrndie(2,'<RTFRDR>','Unrecognized command')
                      errcnt=errcnt+1
                   end select
                endif notblank
             endif
          enddo readingfile
#if KEY_CFF==1
          !     Turn automatic upper case conversion back on.
          if (ffield == cff) ucase=.true.
#endif
          if (.not.qfirst) then
             !             finish-current-residue
             call fresid(wrncnt,patf,patl,totchg,ind, &
                  tblkey,ntabl,errcnt, &
                  natpt,nat,nbopt,nbo,nthpt,nth,nphpt,nph, &
                  nimpt,nim, &
#if KEY_CMAP==1
                  npct,nict,nlpct,nilpt, nanisct, nanis, &
#endif
                  nrxpt,nrx,napt,na,ndpt,nd,nblpt,nbl,ngr, &
                  ok,newres,idel)
          endif
          if(wrncnt > 0 .and. wrnlev >= 2) write(outu,630) wrncnt
          if(errcnt > 0) then
             if(wrnlev >= 2) write(outu,640) errcnt
             call die
          endif
          !           End Procedure READ-FREE-FORMAT-VERSION-18-TOPOLOGY-FILE
       else
          call wrndie(-1,'<RTFRDR>', &
               'This version of the rtf is no longer supported')
       endif
    else formattedfile
       !---------------------------------------------------------
       !   Binary File Reading
       !---------------------------------------------------------
       if (lappe) then
          call wrndie(0,'<RTFRDR>', &
               'APPEND not allowed for binary files.')
       else
          call tryoro(unit,'UNFORMATTED')
          if (reallow) then
             rewind(unit=unit)
          endif
          read(unit) hdr,rtcntl
          if(hdr /= srtf) then
             if(wrnlev >= 2) write(outu,201) srtf,hdr
201          FORMAT(' EXPECTED = "',A4,'" FOUND = "',A4,'"')
             call wrndie(-1,'<RTFIO>','HEADERS DONT MATCH')
          endif
          call rdtitl(title,ntitl,unit,-1)
          call wrtitl(title,ntitl,outu,1)
          if (rtcntl(1) >= 20) then
             !             begin procedure read-version-20-binary-topology-file
             read (unit) nrtrs,rtftyp,rnicm,natct
             call reallocate_rtrs(nrtrs)
             call reallocate_natct(natct)

             if(rnicm /= dim_nicm) then
                if(wrnlev >= 2) write(outu,305) rnicm,dim_nicm
305             FORMAT(' **** ERROR in RTFRDR **** NICM on file = ',I12, &
                     ' different from current value = ',I3)
                CALL WRNDIE(0,'<RTFRDR>', &
                     'BAD NICM VALUE IN TOPOLOGY FILE')
             endif
             autot=.false.
             autod=.false.
             autop=.false.
             if(rtcntl(1) == 17) rtftyp=0
             if (mod(rtftyp,2) == 1) autot=.true.
             if (mod(rtftyp,4) >= 2) autod=.true.
             if (mod(rtftyp,8) >= 4) autop=.true.
             !
             if (rnicm <= dim_nicm) then
                read (unit) ((nic(i,ires),i=1,rnicm),aa(ires), &
                     rtrtyp(ires),aaf(ires),aal(ires),ires=1,nrtrs)
             else
                rnicm=rnicm-dim_nicm
                read (unit) ((nic(i,ires),i=1,dim_nicm),(j2,j=1,rnicm), &
                     aa(ires),rtrtyp(ires),aaf(ires),aal(ires), &
                     ires=1,nrtrs)
             endif
             read (unit) (atct(i),armass(i),i=1,natct)
             call reallocate_rta(nic(ind_at,nrtrs))
             read (unit) (chg(i),ftp(i),mac(i),mxn(i),grpr(i),delat(i), &
                  i=1,nic(ind_at,nrtrs))
             call reallocate_rtb(nic(ind_bo,nrtrs))
             read (unit) (mib(i),mjb(i),delbd(i),i=1,nic(ind_bo,nrtrs))
             call reallocate_rtt(nic(ind_th,nrtrs))
             read (unit) (mit(i),mjt(i),mkt(i),delan(i), &
                  i=1,nic(ind_th,nrtrs))
             call reallocate_rtp(nic(ind_ph,nrtrs))
             read (unit) (mip(i),mjp(i),mkp(i),mlp(i),delpt(i), &
                  i=1,nic(ind_ph,nrtrs))
             call reallocate_rti(nic(ind_im,nrtrs))
             read (unit) (mim(i),mjm(i),mkm(i),mlm(i),delmt(i), &
                  i=1,nic(ind_im,nrtrs))
             call reallocate_rtx(nic(ind_xx,nrtrs))
             read (unit) (mnb(i),i=1,nic(ind_xx,nrtrs))
             call reallocate_rthd(nic(ind_hd,nrtrs))
             read (unit) (md(i),mh(i),mda(i),mdb(i),deldn(i), &
                  i=1,nic(ind_hd,nrtrs))
             call reallocate_rtha(nic(ind_ha,nrtrs))
             read (unit) (ma(i),maa(i),mab(i),delac(i), &
                  i=1,nic(ind_ha,nrtrs))
             call reallocate_rtic(nic(ind_ic,nrtrs))
             read (unit) (icb2(i),icb1(i),icth2(i),icth1(i),icphi(i), &
                  bari(i),barj(i),bark(i),barl(i),bart(i), &
                  delic(i),i=1,nic(ind_ic,nrtrs))
             !             end procedure read-version-20-binary-topology-file
          else
             call wrndie(-1,'<RTFRDR>', &
                  'This version of the rtf is no longer supported')
          endif
       endif
    endif formattedfile
    !
#if KEY_STRINGM==1 /*  VO : restore iolev */
    if (qstr) iolev=oldiol
#endif
    !
#if KEY_PARALLEL==1
   call   rtf_parallel_bcast
#endif

    return

630 FORMAT(' There were ',I3,' warning(s) from RTFRDR.')
640 FORMAT(' There were ',I4,' error(s) from RTFRDR.' &
         /' Execution will be terminated.')
    !

  contains
#if KEY_PARALLEL==1 /*par_ens*/
    subroutine rtf_parallel_bcast
      integer :: status
      integer :: i
      integer dim_rt_at, dim_rt_bo, dim_rt_th, dim_rt_ph, dim_rt_im
      integer dim_rt_xx, dim_rt_hd, dim_rt_ha, dim_rt_ic, dim_rt_cmap
      integer dim_rt_lp, dim_rt_natct, dim_rt_anis
      
      call psnd4(nrtrs,1)
      call reallocate_rtrs(nrtrs)
      call psnd4(natct,1)
      call psnd4(dim_natct,1)
      call psnd4(rtftyp,1)
      call psnd4(rtcntl,20)
      if(nrtrs>0) call psnd4(nic,dim_nicm*nrtrs)
      call psnd4(dim_rt,dim_nicm)
      call psync
      dim_rt_at   = dim_rt(ind_at)
      dim_rt_bo   = dim_rt(ind_bo)
      dim_rt_th   = dim_rt(ind_th)
      dim_rt_ph   = dim_rt(ind_ph)
      dim_rt_im   = dim_rt(ind_im)
      dim_rt_xx   = dim_rt(ind_xx)
      dim_rt_hd   = dim_rt(ind_hd)
      dim_rt_ha   = dim_rt(ind_ha)
      dim_rt_ic   = dim_rt(ind_ic)
#if KEY_CMAP==1
      dim_rt_cmap = dim_rt(ind_cmap)
#endif
#if KEY_CMAP==1
      dim_rt_lp   = dim_rt(ind_lp)
#endif
!      dim_rt_lp   = dim_rt(ind_anis)
      dim_rt_natct= dim_natct
!      dim_rt_anis = 0
      dim_rt_anis = dim_rt(ind_anis)

      call reallocate_natct(dim_rt_natct)
      call reallocate_rta( dim_rt_at)
      call reallocate_rtb( dim_rt_bo)
      call reallocate_rtt( dim_rt_th)
      call reallocate_rtp( dim_rt_ph)
      call reallocate_rti( dim_rt_im)
      call reallocate_rtx( dim_rt_xx)
      call reallocate_rthd(dim_rt_hd)
      call reallocate_rtha(dim_rt_ha)
      call reallocate_rtic(dim_rt_ic)
      call reallocate_natct(dim_rt_ic)
#if KEY_CMAP==1
      call reallocate_cmap(dim_rt_cmap)
#endif
#if KEY_CMAP==1
      call reallocate_lp(dim_rt_lp)
#endif
      call reallocate_anis(dim_rt_anis)
      call psync

      call psnd8(armass,natct)
      call psnd4(autot,1)
      call psnd4(autod,1)
      call psnd4(autop,1)
      if(nrtrs>0) then
         call psndc(aa,nrtrs)
         call psndc(aaf,nrtrs)
         call psndc(aal,nrtrs)
      endif
      call psndc(atct,natct)
      call psndc(deff,1)
      call psndc(defl,1)
      if(nrtrs>0) call psnd4(rtrtyp,nrtrs)
      call psync

     if(nrtrs>0) then

      if(nic(ind_at,nrtrs) > 0) then
         call psnd4(mac,   nic(ind_at,nrtrs))
         call psnd4(mxn,   nic(ind_at,nrtrs))
         call psnd4(grpr,  nic(ind_at,nrtrs))
#if KEY_MMFF==1
         call psnd4(atnumr,nic(ind_at,nrtrs))
#endif
         call psnd8(chg,   nic(ind_at,nrtrs))
         call psnd4(delat, nic(ind_at,nrtrs))
         call psndc(ftp,   nic(ind_at,nrtrs))
        ! begin drude (is this all)
         call psnd8(alph,  nic(ind_at,nrtrs))
         call psnd8(thol,  nic(ind_at,nrtrs))
#if KEY_CHEQ==1
         call psnd8(echh,  nic(ind_at,nrtrs))
#endif
#if KEY_CHEQ==1
         call psnd8(ehah,  nic(ind_at,nrtrs))
#endif
      endif
      call psync

      if(nic(ind_bo,nrtrs) > 0) then
         call psnd4(delbd,nic(ind_bo,nrtrs))
         call psndc(mib  ,nic(ind_bo,nrtrs))
         call psndc(mjb  ,nic(ind_bo,nrtrs))
      endif
      call psync

      if(nic(ind_th,nrtrs) > 0) then
         call psnd4(delan,nic(ind_th,nrtrs))
         call psndc(mit  ,nic(ind_th,nrtrs))
         call psndc(mjt  ,nic(ind_th,nrtrs))
         call psndc(mkt  ,nic(ind_th,nrtrs))
      endif
      call psync

      if(nic(ind_ph,nrtrs) > 0) then
           call psnd4(delpt,nic(ind_ph,nrtrs))
           call psndc(mip  ,nic(ind_ph,nrtrs))
           call psndc(mjp  ,nic(ind_ph,nrtrs))
           call psndc(mkp  ,nic(ind_ph,nrtrs))
           call psndc(mlp  ,nic(ind_ph,nrtrs))
      endif
      call psync

      if(nic(ind_im,nrtrs) > 0) then
         call psnd4(delmt,nic(ind_im,nrtrs))
         call psndc(mim  ,nic(ind_im,nrtrs))
         call psndc(mjm  ,nic(ind_im,nrtrs))
         call psndc(mkm  ,nic(ind_im,nrtrs))
         call psndc(mlm  ,nic(ind_im,nrtrs))
      endif
      call psync

      if(nic(ind_xx,nrtrs) > 0) then
         call psndc(mnb ,nic(ind_xx,nrtrs))
      endif
      call psync

      if(nic(ind_hd,nrtrs) > 0) then
         call psndc(md   ,nic(ind_hd,nrtrs))
         call psndc(mh   ,nic(ind_hd,nrtrs))
         call psndc(mda  ,nic(ind_hd,nrtrs))
         call psndc(mdb  ,nic(ind_hd,nrtrs))
         call psnd4(deldn,nic(ind_hd,nrtrs))
      endif
      call psync

      if(nic(ind_ha,nrtrs) > 0) then
         call psndc(ma   ,nic(ind_ha,nrtrs))
         call psndc(maa  ,nic(ind_ha,nrtrs))
         call psndc(mab  ,nic(ind_ha,nrtrs))
         call psnd4(delac,nic(ind_ha,nrtrs))
      endif
      call psync

      if(nic(ind_ic,nrtrs) > 0) then
         call psnd8(icb1, nic(ind_ic,nrtrs))
         call psnd8(icb2, nic(ind_ic,nrtrs))
         call psnd8(icth1,nic(ind_ic,nrtrs))
         call psnd8(icth2,nic(ind_ic,nrtrs))
         call psnd8(icphi,nic(ind_ic,nrtrs))
         call psnd4(bart, nic(ind_ic,nrtrs))
         call psndc(bari, nic(ind_ic,nrtrs))
         call psndc(barj, nic(ind_ic,nrtrs))
         call psndc(bark, nic(ind_ic,nrtrs))
         call psndc(barl, nic(ind_ic,nrtrs))
         call psnd4(delic,nic(ind_ic,nrtrs))
      endif
      call psync

#if KEY_CMAP==1
      if(nic(ind_cmap,nrtrs) > 0) then
         call psndc(mi1ct, nic(ind_cmap,nrtrs))
         call psndc(mj1ct, nic(ind_cmap,nrtrs))
         call psndc(mk1ct, nic(ind_cmap,nrtrs))
         call psndc(ml1ct, nic(ind_cmap,nrtrs))
         call psndc(mi2ct, nic(ind_cmap,nrtrs))
         call psndc(mj2ct, nic(ind_cmap,nrtrs))
         call psndc(mk2ct, nic(ind_cmap,nrtrs))
         call psndc(ml2ct, nic(ind_cmap,nrtrs))
         call psnd4(delpct,nic(ind_cmap,nrtrs))
      endif
      if(nic(ind_lp,nrtrs) > 0) then
         call psndc(mlp0ct, nic(ind_lp,nrtrs))
         call psndc(mlp1ct, nic(ind_lp,nrtrs))
         call psndc(mlp2ct, nic(ind_lp,nrtrs))
         call psndc(mlp3ct, nic(ind_lp,nrtrs))
         call psndc(mlp4ct, nic(ind_lp,nrtrs))
         call psndc(mlpct1d,nic(ind_lp,nrtrs)*maxcent_hosts)   ! Look for comments in machdep/paral1.src
         call psnd4(dellpct,nic(ind_lp,nrtrs))
      endif
      if(nic(ind_anis,nrtrs) > 0) then
         call psnd4(delanis,nic(ind_anis,nrtrs))
         call psnd4(mianis ,nic(ind_anis,nrtrs))
         call psnd4(mjanis ,nic(ind_anis,nrtrs))
         call psnd4(mkanis ,nic(ind_anis,nrtrs))
         call psnd4(mlanis ,nic(ind_anis,nrtrs))
         call psnd8(a11anis,nic(ind_anis,nrtrs))
         call psnd8(a22anis,nic(ind_anis,nrtrs))
      endif
#endif
   endif
      return
    end subroutine rtf_parallel_bcast
#endif /* (par_ens)*/


    !**********************************************************
    subroutine get_rtf_version
      eof=.false.
      prnflg=(qprint .and. prnlev >= 2)
      if(nrtrs == 0.and.lappe) &
           CALL WRNDIE(-1,'<RTFRDR>','Cant append to zero RTF')

      call rdcmnd(comlyn,mxcmsz,comlen,unit,eof,.false., &
           prnflg,'rtf>    ')
      if(eof) call wrndie(0,'<rtfrdr>', &
           'Unexpected end of file during read')
      rtcntl(1)=nexti(comlyn,comlen)

      !         No version number specified.
      if (rtcntl(1) == 0) then
         backspace (unit=unit)
         rtcntl(1)=18
         if(wrnlev >= 2) write (outu,'(a)') &
              ' RTFRDR> WARNING: Version number is NOT specified.'
      else
         i=1
         do while (comlen > 0)
            i=i+1
            rtcntl(i)=nexti(comlyn,comlen)
         enddo
      endif
      return
    end subroutine get_rtf_version



    subroutine rtfrdr_process_mass
      !                 Begin Procedure PROCESS-MASS
      !
      !                 MASS atom-type-code atom-type-name mass element-name
      !
      iatc=nexti(comlyn,comlen)
      if (iatc < 1 .and. iatc /= -1) then
         if(wrnlev >= 2) write(outu,'(a,i12,a)') &
              ' **** ERROR in RTFRDR **** Bad atom type code = ', &
              IATC,' was specified. It was ignored.'
         errcnt=errcnt+1
      else
         if(iatc == -1) then
            if(prnlev>7) write(outu,'(10x,a)') &
                 'CHARMM> RTFIO: expanding natct because iatc = -1'
            iatc = natct+1
         endif
         call reallocate_natct(iatc)
         if(atct(iatc) /= blank.and..not.lappe) then
            call wrndie(2,'<RTFRDR>', &
                 'Duplicate atom type index. replacing')
            wrncnt=wrncnt+1
         endif
         natct=max(iatc,natct)
         name=nexta6(comlyn,comlen)
         if(srchws(atct,natct,name) /= 0.and..not.lappe) then
            call wrndie(2,'<RTFRDR>', &
                 'Duplicate atom type name')
            errcnt=errcnt+1
         endif
         atct(iatc)=name
         armass(iatc)=nextf(comlyn,comlen)
         if(armass(iatc) < 0.0) then
            call wrndie(2,'<RTFRDR>', &
                 'Negative mass specified')
            errcnt=errcnt+1
         endif
         name=nexta6(comlyn,comlen)
#if KEY_CFF==1
         if (.not.ucase) call cnvtuc(name,6)
#endif
#if KEY_MMFF==1
         if(name(1:1) == '-') then
            AtNumT(IATC)=-AtomicNumber(name(2:))
            if(prnlev >= 2) write(OUTU,'(a,a,2i5)') &
                 'RTFRDR> dummy atom: name,iatc,AtNumT=', &
                 name,iatc,AtNumT(IATC)
         else
            AtNumT(IATC)=AtomicNumber(name)
         endif
         if(AtNumT(IATC) == 0 .and. FFIELD.EQ.MMFF) then
            SCRTCH='Unknown atomic symbol >'//name//'<'// &
                 ' for atom type '//ATCT(IATC)
            call wrndie(-1,'<RTFRDR>', &
                 SCRTCH(:LEN_TRIM(SCRTCH)))
         endif
#endif
         call xtrane(comlyn,comlen,'rtf reader')
      endif
      !                 end procedure process-mass
    end subroutine rtfrdr_process_mass


    subroutine rtfrdr_process_declare
      !                 Begin Procedure PROCESS-DECLARE
      !
      !                 DECLARE name value
      !
      keytab=nexta6(comlyn,comlen)
      if (srchws(tblkey,ntabl,keytab) /= 0) then
         call wrndie(2,'<RTFRDR>', &
              'Symbol redeclared. New declaration ignored')
         errcnt=errcnt+1
      else
         ntabl=ntabl+1
         if(ntabl > mxtabl) then
            call wrndie(0,'<RTFRDR>', &
                 'Too many declarations - MXTABL exceeded')
         endif
         tblkey(ntabl)=keytab
         call xtrane(comlyn,comlen,'RTF reader')
      endif
      !                 End Procedure PROCESS-DECLARE
    end subroutine rtfrdr_process_declare

    subroutine rtfrdr_process_residue
      if (.not.qfirst) then
         !                   finish-current-residue
         call fresid(wrncnt,patf,patl,totchg,ind, &
              tblkey,ntabl,errcnt, &
              natpt,nat,nbopt,nbo,nthpt,nth,nphpt,nph, &
              nimpt,nim, &
#if KEY_CMAP==1
              npct,nict, nlpct,nilpt, nanisct, nanis, &
#endif
              nrxpt,nrx,napt,na,ndpt, &
              nd,nblpt,nbl,ngr, &
              ok,newres,idel)
      else
         qfirst=.false.
      endif
      !                 PROCESS-RESIDUE
      call prores(errcnt,newres,idel,wrncnt,totchg, &
           ipatch,patf,patl, &
           nat,nbo,nth,nph,nim, &
#if KEY_CMAP==1
           nict,nilpt,nanis,  &
#endif
           nrx,na,nd,nbl,ngr)
    end    subroutine rtfrdr_process_residue

    subroutine rtfrdr_process_atoms

      !                 Begin Procedure PROCESS-ATOM
      !
      !               ATOM iupac atom-type-code charge repeat(exclusion-names)
      !               DELETE ATOM iupac  !only for patch residues.
      !
      if (nrtrs == 0) then
         call wrndie(2,'<RTFRDR>', &
              'A residue must be specified before anything ' &
              //'is put into it.')
         errcnt=errcnt+1
      else
         nat=nat+1
         call reallocate_rta(nat+natpt)
         name=nexta6(comlyn,comlen)
         if (name == blank) then
            call wrndie(2,'<RTFRDR>','Atoms must have names')
            errcnt=errcnt+1
         else
#if KEY_CFF==1
            if (.not.ucase) call cnvtuc(name,6)
#endif
            if (srchws(ftp(natpt+1),nat-1,name) > 0) then
               call wrndie(2,'<RTFRDR>','Duplicate atom name')
               errcnt=errcnt+1
            else
               ftp(nat+natpt)=name
               delat(nat+natpt)=qdelet
               if (qdelet) then
                  !
                  !   process the combine option for the delete atom command
                  !
                  alph(nat+natpt)= &
                       gtrmf(comlyn,comlen,'ALPH',zero)
                  name=nexta4(comlyn,comlen)
                  if (name == 'COMB') then
                     name=nexta6(comlyn,comlen)
                     nrx=nrx+1
                     call reallocate_rtx(nrx+nrxpt)
                     if(nrx+nrxpt > dim_rt(ind_xx)) then
                        call wrndie(0,'<RTFRDR>rtfio.src', &
                             'Limit of non-bonded exclusions ' &
                             //'dim_rt(ind_xx) been exceeded')
                     endif
                     mnb(nrx+nrxpt)=name
                     if(prnflg) write(outu,'(a,a6,a,a,a6,a)') &
                          ' Atom "',ftp(nat+natpt), &
                          '" is to be deleted and references moved to', &
                          ' atom "',name,'".'
                  elseif (name == blank) then
                  else
                     call wrndie(2,'<RTFRDR>','Error parsing ' &
                          //'the DELEte ATOM option.')
                     errcnt=errcnt+1
                  endif

                  chg(nat+natpt)=zero
                  mxn(nat+natpt)=nrx+nrxpt

                  if(ldrude)then
                     if(alph(nat+natpt) /= zero)then
                        name='d'//ftp(nat+natpt)
                        alph(nat+natpt)=zero
                        thol(nat+natpt)=zero
                        
                        nat=nat+1
                        call reallocate_rta(nat+natpt)
                        ftp(nat+natpt)=name
                        delat(nat+natpt)=qdelet
                        chg(nat+natpt)=zero
                        alph(nat+natpt)=zero
                        thol(nat+natpt)=zero
                        if(prnlev > 5) &
                             write(outu,'(4a,3(a,f8.3))')  &
                             ' delete drude in   ',aa(nrtrs), &
                             ' atom type=',FTP(NAT+NATPT)
                     endif
                  endif

               else

                  call reallocate_rta(nat+natpt)
                  typenm=nexta6(comlyn,comlen)
                  if (typenm == blank) then
                     typecd=0
                  else
                     typecd=srchws(atct,natct,typenm)
                  endif
                  if (typecd == 0) then
                     if(wrnlev >= 2) write(outu,'(a,a6,a)') &
                          ' **** ERROR in RTFRDR ****  The atom type code, ', &
                          typenm,' is unknown. The atom will be ignored.'
                     errcnt=errcnt+1
                  else
                     mac(nat+natpt)=typecd
                     chg(nat+natpt)=nextf(comlyn,comlen)
#if KEY_CHEQ==1
                     tech= 0.d0
                     teca=-1.d0
                     tech=gtrmf(comlyn,comlen,'ECHH',tech)
                     teca=gtrmf(comlyn,comlen,'ECHA',teca)
                     echh(nat+natpt)=tech
                     ehah(nat+natpt)=teca
                     ehah(nat+natpt)=ehah(nat+natpt)*factkcalr
#endif
                     grpr(nat+natpt)=ngr

                     ifdrude: if(ldrude)then
                        ! build a drude only on atoms that have a non-zero alpha
                        alph(nat+natpt)= &
                             gtrmf(comlyn,comlen,'ALPH',zero)
                        if(alph(nat+natpt) /= zero)then
                           ! plan for patches involving multiple residues...
                           if(ftp(nat+natpt)(1:1) == '1')then
                              name='1D'//ftp(nat+natpt)(2:)
                           elseif(ftp(nat+natpt)(1:1) == '2')then
                              name='2D'//ftp(nat+natpt)(2:)
                           elseif(ftp(nat+natpt)(1:1) == '3')then
                              name='3D'//ftp(nat+natpt)(2:)
                           else
                              name='D'//ftp(nat+natpt)
                           endif
                           thol(nat+natpt)= &
                                gtrmf(comlyn,comlen,'THOL',tholea)
                           nat=nat+1
                           call reallocate_rta(nat+natpt)
                           ftp(nat+natpt)=name
                           delat(nat+natpt)=qdelet
                           chg(nat+natpt)=0.0d0
                           alph(nat+natpt)=0.0d0
                           thol(nat+natpt)=0.0d0
                           grpr(nat+natpt)=ngr
                           typenm='DRUD'
                           call gtrmwd(comlyn,comlen,'TYPE',4, &
                                typenm,6,j2)
                           typecd=srchws(atct,natct,typenm)
                           mac(nat+natpt)=typecd
                           if(prnlev > 5) &
                                write(outu,'(6a,a,i4,4(a,f8.3))')  &
                                ' Add drude atom in ', AA(NRTRS), &
                                ' atom type=',FTP(NAT+NATPT), &
                                ' chem type =',ATCT(MAC(NAT+NATPT)), &
                                ' chem type number=',MAC(NAT+NATPT), &
                                ' charge=',CHG(NAT+NATPT), &
                                ' alpha=',ALPH(NAT+NATPT), &
                                ' thole=',ALPH(NAT+NATPT), &
                                ' mass=',ARMASS(MAC(NAT+NATPT))
                           nbo=nbo+1                 ! add bond immediately
                           call reallocate_rtb(nbo+nbopt)
                           w1=ftp(nat+natpt-1)
                           w2=ftp(nat+natpt)
                           mib(nbo+nbopt)=w1
                           mjb(nbo+nbopt)=w2
                           if(prnlev > 5) &
                                write(outu,'(a,i4,2a)')  &
                                ' Add bond ',NBO,W1,W2
                        endif
                     endif ifdrude

#if KEY_MMFF==1
                     AtNumR(NAT+NATPT)=AtNumT(TYPECD)
#endif
                     name=nexta6(comlyn,comlen)
                     do while (name /= blank)
                        nrx=nrx+1
                        call reallocate_rtx(nrx+nrxpt)
                        if(nrx+nrxpt > dim_rt(ind_xx)) &
                             CALL WRNDIE(0,'<RTFRDR>rtfio.src', &
                             'Limit of non-bonded exclusions ' &
                             //'dim_rt(ind_xx) been exceeded')
                        mnb(nrx+nrxpt)=name
                        name=nexta6(comlyn,comlen)
                     enddo
                     mxn(nat+natpt)=nrx+nrxpt
                  endif
               endif
            endif
         endif
      endif
      !                 End Procedure PROCESS-ATOM
    end subroutine rtfrdr_process_atoms

    subroutine rtfrdr_process_bond
      !                 Begin Procedure PROCESS-BOND
      !
      !                 [DELETE] BOND repeat(iupac iupac)
      !
      if (nrtrs == 0) then
         call wrndie(2,'<RTFRDR>', &
              'A residue must be specified before anything ' &
              //'is put into it.')
         errcnt=errcnt+1
      else
         w1=nexta6(comlyn,comlen)
         do while (w1 /= blank)
            w2=nexta6(comlyn,comlen)
            if (w2 == blank) then
               call wrndie(2,'<RTFRDR>', &
                    'A lone atom was specified as a bond. It ' &
                    //'will be ignored.')
               errcnt=errcnt+1
            else
               name=w1
               !                       lookup-name-into-ind
               call lookpn(ind,natpt,nat,name,ok,tblkey, &
                    ntabl,errcnt)
               if(ok) then
                  name=w2
                  !                         lookup-name-into-ind
                  call lookpn(ind,natpt,nat,name,ok,tblkey, &
                       ntabl,errcnt)
                  if(ok) then
                     if (w1 == w2) then
                        call wrndie(2,'<RTFRDR>', &
                             'Two atoms in bond are identical. ' &
                             //'They will be ignored.')
                        errcnt=errcnt+1
                     else
                        nbo=nbo+1
                        call reallocate_rtb(nbo+nbopt)
                        delbd(nbo+nbopt)=qdelet
                        mib(nbo+nbopt)=w1
                        mjb(nbo+nbopt)=w2
#if KEY_MMFF==1
                        if(WRD == 'BOND' .or. WRD.eq.'SING') then
                           MBTYPE(NBO+NBOPT)=SINGLE
                        elseif(WRD == 'DOUB') then
                           MBTYPE(NBO+NBOPT)=DOUBLE
                        elseif(WRD == 'TRIP') then
                           MBTYPE(NBO+NBOPT)=TRIPLE
                        elseif(WRD == 'AROM') then
                           MBTYPE(NBO+NBOPT)=AROMATIC
                        else
                           MBTYPE(NBO+NBOPT)=SINGLE
                           CALL WRNDIE(-2,'<RTFRDR>', &
                                'Unknown bond type')
                        endif
#endif
                     endif
                  endif
               endif
            endif
            w1=nexta6(comlyn,comlen)
         enddo
         call xtrane(comlyn,comlen,'RTF reader')
      endif
      !                 End Procedure PROCESS-BOND
    end subroutine rtfrdr_process_bond


    subroutine rtfrdr_process_angle
      !                 Begin Procedure PROCESS-ANGLE
      !
      !                 [DELETE] ANGLE repeat(iupac iupac iupac)
      !
      if (nrtrs == 0) then
         call wrndie(2,'<RTFRDR>', &
              'A residue must be specified before anything ' &
              //'is put into it.')
         errcnt=errcnt+1
      else
         w1=nexta6(comlyn,comlen)
         do while (w1 /= blank)
            w2=nexta6(comlyn,comlen)
            w3=nexta6(comlyn,comlen)
            if (w2 == blank.or.w3.eq.blank) then
               call wrndie(2,'<RTFRDR>', &
                    'Atoms were missing in the last angle. It ' &
                    //'will be ignored.')
               errcnt=errcnt+1
            else
               name=w1
               !                       LOOKUP-NAME-INTO-IND
               call lookpn(ind,natpt,nat,name,ok,tblkey, &
                    ntabl,errcnt)
               if(ok) then
                  name=w2
                  !                         lookup-name-into-ind
                  call lookpn(ind,natpt,nat,name,ok,tblkey, &
                       ntabl,errcnt)
                  if(ok) then
                     name=w3
                     !                           lookup-name-into-ind
                     call lookpn(ind,natpt,nat,name,ok,tblkey, &
                          ntabl,errcnt)
                     if(ok) then
                        ok=w1 /= w2.and.w1 /= w3.and.w2 /= w3
                        if (.not.ok) then
                           call wrndie(2,'<RTFRDR>', &
                                'The atoms in an angle are not all ' &
                                //'different. Will be ignored')
                           errcnt=errcnt+1
                        else
                           nth=nth+1
                           call reallocate_rtt(nth+nthpt)
                           delan(nth+nthpt)=qdelet
                           mit(nth+nthpt)=w1
                           mjt(nth+nthpt)=w2
                           mkt(nth+nthpt)=w3
                        endif
                     endif
                  endif
               endif
            endif
            w1=nexta6(comlyn,comlen)
         enddo
         call xtrane(comlyn,comlen,'RTF reader')
      endif
      !                 End Procedure PROCESS-ANGLE
    end subroutine rtfrdr_process_angle


    subroutine rtfrdr_process_dihedral
      !                 Begin Procedure PROCESS-TORSION
      !
      !                 [DELETE] { TORSION  } repeat(iupac iupac iupac iupac)
      !                 { DIHEDRAL }
      !
      if (nrtrs == 0) then
         call wrndie(2,'<RTFRDR>', &
              'A residue must be specified before anything ' &
              //'is put into it.')
         errcnt=errcnt+1
      else
         w1=nexta6(comlyn,comlen)
         do while (w1 /= blank)
            w2=nexta6(comlyn,comlen)
            w3=nexta6(comlyn,comlen)
            w4=nexta6(comlyn,comlen)
            if (w2 == blank.or.w3.eq.blank.or.w4.eq.blank) &
                 then
               call wrndie(2,'<rtfrdr>', &
                    'Atoms were missing in the last torsion. ' &
                    //'It will be ignored.')
               errcnt=errcnt+1
            else
               name=w1
               !                       lookup-name-into-ind
               call lookpn(ind,natpt,nat,name,ok,tblkey, &
                    ntabl,errcnt)
               iw1=ind
               if(ok) then
                  name=w2
                  !                         lookup-name-into-ind
                  call lookpn(ind,natpt,nat,name,ok,tblkey, &
                       ntabl,errcnt)
                  if(ok) then
                     name=w3
                     !                           lookup-name-into-ind
                     call lookpn(ind,natpt,nat,name,ok,tblkey, &
                          ntabl,errcnt)
                     if(ok) then
                        name=w4
                        !                             lookup-name-into-ind
                        call lookpn(ind,natpt,nat,name,ok, &
                             tblkey,ntabl,errcnt)
                        iw4=ind
                        if(ok) then
                           ok=w1 /= w2.and.w1 /= w3.and. &
                                w1 /= w4.and.w2 /= w3.and. &
                                w2 /= w4.and.w3 /= w4
                           if (.not.ok) then
                              call wrndie(2,'<rtfrdr>', &
                                   'the atoms in a torsion are not ' &
                                   //'all different. will be ignored')
                              errcnt=errcnt+1
                           else
                              if (iw1 > iw4) then
                                 wtemp=w1
                                 w1=w4
                                 w4=wtemp
                                 wtemp=w2
                                 w2=w3
                                 w3=wtemp
                              endif
                              addok=.true.
                              do i=nphpt+1,nphpt+nph
                                 if (mip(i) == w1.and.mjp(i).eq.w2.and. &
                                      mkp(i) == w3.and.mlp(i).eq.w4) &
                                      addok=.false.
                              enddo
                              if(addok) then
                                 nph=nph+1
                                 call reallocate_rtp(nph+nphpt)
                                 delpt(nph+nphpt)=qdelet
                                 mip(nph+nphpt)=w1
                                 mjp(nph+nphpt)=w2
                                 mkp(nph+nphpt)=w3
                                 mlp(nph+nphpt)=w4
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
            w1=nexta6(comlyn,comlen)
         enddo
         call xtrane(comlyn,comlen,'RTF reader')
      endif
      !                 end procedure process-torsion
    end subroutine rtfrdr_process_dihedral

    !**********************************************************
    subroutine rtfrdr_process_improper
      !                 Begin Procedure PROCESS-IMPROPER
      !
      !                 [DELETE] { IMPROPER } repeat(iupac iupac iupac iupac)
      !                 {  IMPHI   }
      !
      if (nrtrs == 0) then
         call wrndie(2,'<RTFRDR>', &
              'A residue must be specified before anything ' &
              //'is put into it.')
         errcnt=errcnt+1
      else
         w1=nexta6(comlyn,comlen)
         do while (w1 /= blank)
            w2=nexta6(comlyn,comlen)
            w3=nexta6(comlyn,comlen)
            w4=nexta6(comlyn,comlen)
            if (w2 == blank.or.w3.eq.blank.or.w4.eq.blank) &
                 then
               call wrndie(2,'<RTFRDR>', &
                    'Atoms were missing in the last improper ' &
                    //'torsion. Ignored.')
               errcnt=errcnt+1
            else
               name=w1
               !                       lookup-name-into-ind
               call lookpn(ind,natpt,nat,name,ok,tblkey, &
                    ntabl,errcnt)
               if(ok) then
                  name=w2
                  !                         lookup-name-into-ind
                  call lookpn(ind,natpt,nat,name,ok,tblkey, &
                       ntabl,errcnt)
                  if(ok) then
                     name=w3
                     !                           lookup-name-into-ind
                     call lookpn(ind,natpt,nat,name,ok,tblkey, &
                          ntabl,errcnt)
                     if(ok) then
                        name=w4
                        !                             lookup-name-into-ind
                        call lookpn(ind,natpt,nat,name,ok, &
                             tblkey,ntabl,errcnt)
                        if(ok) then
                           ok &
                                =w1 /= w2.and.w1 /= w3.and.w1 /= w4.and. &
                                w2 /= w3.and.w2 /= w4.and.w3 /= w4
                           if (.not.ok) then
                              call wrndie(2,'<RTFRDR>', &
                                   'The atoms in an improper ' &
                                   //'torsion are not all different.')
                              errcnt=errcnt+1
                           else
                              nim=nim+1
                              call reallocate_rti(nim+nimpt)
                              delmt(nim+nimpt)=qdelet
                              mim(nim+nimpt)=w1
                              mjm(nim+nimpt)=w2
                              mkm(nim+nimpt)=w3
                              mlm(nim+nimpt)=w4
                           endif
                        endif
                     endif
                  endif
               endif
            endif
            w1=nexta6(comlyn,comlen)
         enddo
         call xtrane(comlyn,comlen,'RTF reader')
      endif
      !                 End Procedure PROCESS-IMPROPER
      !
      return
    end subroutine rtfrdr_process_improper


    subroutine rtfrdr_process_anisotropic
      !                 Begin Procedure PROCESS anisotropic polarizability
      !
      !                 [DELETE] { ANISOTROPY } 4-atoms  A11 {real}   A22 {real}
      !
      if (nrtrs == 0) then
         call wrndie(2,'<RTFRDR>', &
              'A residue must be specified before anything ' &
              //'is put into it.')
         errcnt=errcnt+1
      else

         a11=gtrmf(comlyn,comlen,'A11',one)
         a22=gtrmf(comlyn,comlen,'A22',one)

         w1=nexta6(comlyn,comlen)
         do while (w1 /= blank)
            w2=nexta6(comlyn,comlen)
            w3=nexta6(comlyn,comlen)
            w4=nexta6(comlyn,comlen)
            if (w2 == blank.or.w3.eq.blank.or.w4.eq.blank) &
                 then
               call wrndie(2,'<rtfrdr>', &
                    'Atoms were missing in the anisotropic' &
                    //'alpha. Ignored.')
               errcnt=errcnt+1
            else
               name=w1
               !                       lookup-name-into-ind
               call lookpn(ind,natpt,nat,name,ok,tblkey, &
                    ntabl,errcnt)
               if(ok) then
                  name=w2
                  !                         lookup-name-into-ind
                  call lookpn(ind,natpt,nat,name,ok,tblkey, &
                       ntabl,errcnt)
                  if(ok) then
                     name=w3
                     !                           lookup-name-into-ind
                     call lookpn(ind,natpt,nat,name,ok,tblkey, &
                          ntabl,errcnt)
                     if(ok) then
                        name=w4
                        !                             lookup-name-into-ind
                        call lookpn(ind,natpt,nat,name,ok, &
                             tblkey,ntabl,errcnt)
                        if(ok) then
                           ok &
                                =w1 /= w2.and.w1 /= w3.and.w1 /= w4.and. &
                                w2 /= w3.and.w2 /= w4.and.w3 /= w4
                           if (.not.ok) then
                              call wrndie(2,'<rtfrdr>', &
                                   'The atoms in anisotropic ' &
                                   //'term are not all different.')
                              errcnt=errcnt+1
                           else
#if KEY_CMAP==1
                              nanis=nanis+1
                              call reallocate_anis(nanis+nanisct)
                              delanis(nanis+nanisct)=qdelet
                              mianis(nanis+nanisct)=w1
                              mjanis(nanis+nanisct)=w2
                              mkanis(nanis+nanisct)=w3
                              mlanis(nanis+nanisct)=w4
                              a11anis(nanis+nanisct)=a11
                              a22anis(nanis+nanisct)=a22
#endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
            w1=nexta6(comlyn,comlen)
         enddo
         call xtrane(comlyn,comlen,'RTF reader')
      endif
      !                 end procedure process-anisotropic-polarizability
      !
    end subroutine rtfrdr_process_anisotropic


    subroutine rtfrdr_process_lonepairs
      !                 Begin Procedure PROCESS-lone-pairs
      !
      !                 [DELETE] { LONEPAIR } dist {real} angle {real} dihe {real}  4-atoms
      !                 or
      !                 [DELETE] COLI dist {real} scal {real}  3-atoms
      !  Modified to introduce Lone pair center from rtf: Alexander D. MacKerell and Anmol Kumar
      integer:: i,j
      if (nrtrs == 0) then
         call wrndie(2,'<rtfrdr>', &
              'A residue must be specified before anything ' &
              //'is put into it.')
         errcnt=errcnt+1
      else
! Check if I have to remove KEY_CMAP from here and why.
#if KEY_CMAP==1
      w5=nexta4(comlyn,comlen)
#endif /* KEY_CMAP */
      dist=gtrmf(comlyn,comlen,'DIST',zero)
      angl=gtrmf(comlyn,comlen,'ANGL',zero)
      dihe=gtrmf(comlyn,comlen,'DIHE',zero)
      SCALE=GTRMF(COMLYN,COMLEN,'SCAL',ZERO)
#if KEY_CMAP==1
      if (w5 == 'CENT') then
         icenthst=0
         i = 0
123      i = i + 1
         w(i) = nexta6(comlyn,comlen)
         if (w(i) == blank.and.i==1 ) then
             call wrndie(2,'<rtfrdr>', &
                    'Atoms were missing in the last cross-term ' &
                    //'map. Ignored.')
             errcnt=errcnt+1
         else if (w(i) == blank.and.i.eq.2 ) then
             call wrndie(2,'<rtfrdr>', &
                    'Lone pair center requires more than one host')
             errcnt=errcnt+1
         else if (i.ge.maxcent_hosts) then
             call wrndie(2,'<rtfrdr>', &
                    'Lone pair center cannot have more than 100 host atoms')
             errcnt=errcnt+1
         else if (w(i) == blank.and.i.gt.2 ) then
              go to 124
         else if (w(i)/= blank) then
             name=w(i)
             ! lookup-name-into-ind
             call lookpn(ind,natpt,nat,name,ok,tblkey,ntabl,errcnt)
             if (ok) then
                if (i.gt.1) then
                  do j=1,i
                    if (j.lt.i) then
                     if (w(j).eq.w(j+1)) then
                      call wrndie(2,'<rtfrdr>', &
                           'The atoms map ' &
                           //'are not all different.')
                      errcnt=errcnt+1
                     end if
                    end if
                  end do
                end if
             else
                call wrndie(2,'<rtfrdr>', &
                           'The atoms ' &
                           //'are not present in residue.')
                errcnt=errcnt+1
             end if
             go to 123
          end if
124       continue
          icenthst = i   ! The loop comes here when it finds a blank host name, thus one greater than actual number of hosts. This helps in genpsf
          nilpt=nilpt+1
          call reallocate_lp(nilpt+nlpct)
          dellpct(nilpt+npct)=qdelet
          mlp0ct(nilpt+nlpct)=w5
          mlp1ct(nilpt+nlpct)=blank
          mlp2ct(nilpt+nlpct)=blank
          mlp3ct(nilpt+nlpct)=blank
          mlp4ct(nilpt+nlpct)=blank
     !     if (nilpt+nlpct.ge.1000) then
     !        call wrndie(2,'<rtfrdr>', &
     !               'Lone pair center command fails if it is read after 1000th entry &
     !                of lone pair among all the toppar files. Code needs modification')
     !        errcnt=errcnt+1
     !     end if
          !do j= 1,icenthst
          !     mlpct(nilpt+nlpct,j)=w(j)
          !end do
          do j = 1, icenthst
             mlpct1d((nilpt+nlpct-1)*maxcent_hosts + j) = w(j) 
          enddo
          ncenthst(nilpt+nlpct)=icenthst-1    ! Storing the actual number of host for ith lone pair
          if (ncenthst(nilpt+nlpct) < maxcent_hosts) then
             do j = ncenthst(nilpt+nlpct)+1, maxcent_hosts
                 mlpct1d((nilpt+nlpct-1)*maxcent_hosts + j) = blank
             enddo
          endif
          rlpct(nilpt+nlpct)= 0
          tlpct(nilpt+nlpct)= 0
          plpct(nilpt+nlpct)= 0
      else if(w5 == 'RELA'.or. w5 == 'BISE')then
            w1=nexta6(comlyn,comlen)
            w2=nexta6(comlyn,comlen)
            w3=nexta6(comlyn,comlen)
            w4=nexta6(comlyn,comlen)
            if (w1 == blank.or.w2.eq.blank.or.w3.eq.blank.or. &
                 w4 == blank) &
                 then
               call wrndie(2,'<rtfrdr>', &
                    'Atoms were missing in the last cross-term ' &
                    //'map. Ignored.')
               errcnt=errcnt+1
            else
               name=w1
               !                       lookup-name-into-ind
               call lookpn(ind,natpt,nat,name,ok,tblkey, &
                    ntabl,errcnt)
               if(ok) then
                  name=w2
                  !                         lookup-name-into-ind
                  call lookpn(ind,natpt,nat,name,ok,tblkey, &
                       ntabl,errcnt)
                  if(ok) then
                     name=w3
                     !                           lookup-name-into-ind
                     call lookpn(ind,natpt,nat,name,ok,tblkey, &
                          ntabl,errcnt)
                     if(ok) then
                        name=w4
                        !                             lookup-name-into-ind
                        call lookpn(ind,natpt,nat,name,ok, &
                             tblkey,ntabl,errcnt)

                        if(ok) then
                           ok &
                                =w1 /= w2.and.w1 /= w3.and.w1 /= w4.and. &
                                w2 /= w3.and.w2 /= w4.and.w3 /= w4

                           if (.not.ok) then
                              call wrndie(2,'<rtfrdr>', &
                                   'The atoms map ' &
                                   //'are not all different.')
                              errcnt=errcnt+1
                           else
                              nilpt=nilpt+1
                              icenthst=0
                              call reallocate_lp(nilpt+nlpct)
                              dellpct(nilpt+npct)=qdelet
                              mlp0ct(nilpt+nlpct)=w5
                              mlp1ct(nilpt+nlpct)=w1
                              mlp2ct(nilpt+nlpct)=w2
                              mlp3ct(nilpt+nlpct)=w3
                              mlp4ct(nilpt+nlpct)=w4
                              do j = 1, size(w)
                                 mlpct1d((nilpt+nlpct-1)*maxcent_hosts + j) = blank
                              enddo    
                             ! do j=1,size(w)
                             !    mlpct(nilpt+nlpct,j)=blank   ! For sanity
                             ! end do
                              ncenthst(nilpt+nlpct)= 0
                              rlpct(nilpt+nlpct)=dist
                              tlpct(nilpt+nlpct)=angl
                              plpct(nilpt+nlpct)=dihe
                              if(prnlev > 5) write(outu, &
                                   '(1x,a,2i6,1x,5(a,2x),3f10.3)')  &
                                   'Lone Pair ', &
                                   nilpt,nlpct, &
                                   mlp0ct(nilpt+nlpct),  &
                                   mlp1ct(nilpt+nlpct),  &
                                   mlp2ct(nilpt+nlpct),  &
                                   mlp3ct(nilpt+nlpct),  &
                                   mlp4ct(nilpt+nlpct),  &
                                   rlpct(nilpt+nlpct),  &
                                   tlpct(nilpt+nlpct),  &
                                   plpct(nilpt+nlpct)
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         !-------- for colinear case  ------------------------------
         elseif(w5 == 'COLI')then
            w1=nexta6(comlyn,comlen)
            w2=nexta6(comlyn,comlen)
            w3=nexta6(comlyn,comlen)
            if (w1 == blank.or.w2.eq.blank.or.w3.eq.blank) &
                 then
               call wrndie(2,'<RTFRDR>', &
                    'Atoms were missing in the last cross-term ' &
                    //'map. Ignored.')
               errcnt=errcnt+1
            else
               name=w1
               !                       lookup-name-into-ind
               call lookpn(ind,natpt,nat,name,ok,tblkey,ntabl,errcnt)
               if(ok) then
                  name=w2
                  !                         lookup-name-into-ind
                  call lookpn(ind,natpt,nat,name,ok,tblkey,ntabl,errcnt)
                  if(ok) then
                     name=w3
                     !                           lookup-name-into-ind
                     call lookpn(ind,natpt,nat,name,ok,tblkey,ntabl,errcnt)
                     if(ok) then
                        ok = w1 /= w2.and.w1 /= w3.and. w2 /= w3
                        if (.not.ok) then
                           CALL WRNDIE(2,'<RTFRDR>', &
                                'The atoms map ' &
                                //'are not all different.')
                           errcnt=errcnt+1
                        else  ! still use the data structure of
                           !    RELATI and BISECT case for consistency HJ
                           nilpt=nilpt+1
                           icenthst=0
                           call reallocate_lp(nilpt+nlpct)
                           dellpct(nilpt+npct)=qdelet
                           mlp0ct(nilpt+nlpct)=w5
                           mlp1ct(nilpt+nlpct)=w1
                           mlp2ct(nilpt+nlpct)=w2
                           mlp3ct(nilpt+nlpct)=w3
                           mlp4ct(nilpt+nlpct)=w1
                           do j = 1, size(w)
                              mlpct1d((nilpt+nlpct-1)*maxcent_hosts + j) = blank
                           enddo    
                           !do j=1,size(w)
                           !   mlpct(nilpt+nlpct,j)=blank   ! For sanity
                           !end do
                           ncenthst(nilpt+nlpct)= 0
                           rlpct(nilpt+nlpct)=dist
                           tlpct(nilpt+nlpct)=scale
                           plpct(nilpt+nlpct)=0
                           if(prnlev > 5) write(outu, &
                                '(1x,a,2i6,1x,4(a,2x),2f10.3)')  &
                                'Lone Pair ', &
                                nilpt,nlpct, &
                                mlp0ct(nilpt+nlpct),  &
                                mlp1ct(nilpt+nlpct),  &
                                mlp2ct(nilpt+nlpct),  &
                                mlp3ct(nilpt+nlpct),  &
                                rlpct(nilpt+nlpct),  &
                                tlpct(nilpt+nlpct)
                        endif
                     endif
                  endif
               endif
            endif
         endif
#endif

         call xtrane(comlyn,comlen,'RTF reader')
      endif
      !                 End Procedure PROCESS-LONE-PAIRS
    end subroutine rtfrdr_process_lonepairs


#if KEY_CMAP==1
    subroutine rtfrdr_process_cmap
      !                 Begin Procedure PROCESS-CROSSTERM-MAP
      !
      !                 [DELETE] { CMAP } repeat(8 x iupac)
      !
      if (nrtrs == 0) then
         call wrndie(2,'<rtfrdr>', &
              'A residue must be specified before anything ' &
              //'is put into it.')
         errcnt=errcnt+1
      else
         w1=nexta6(comlyn,comlen)
         do while (w1 /= blank)
            w2=nexta6(comlyn,comlen)
            w3=nexta6(comlyn,comlen)
            w4=nexta6(comlyn,comlen)
            w5=nexta6(comlyn,comlen)
            w6=nexta6(comlyn,comlen)
            w7=nexta6(comlyn,comlen)
            w8=nexta6(comlyn,comlen)
            if (w2 == blank.or.w3.eq.blank.or.w4.eq.blank.or. &
                 w5 == blank.or.w6.eq.blank.or.w7.eq.blank.or. &
                 w8 == blank) &
                 then
               call wrndie(2,'<rtfrdr>', &
                    'Atoms were missing in the last cross-term ' &
                    //'map. Ignored.')
               errcnt=errcnt+1
            else
               name=w1
               !                       lookup-name-into-ind
               call lookpn(ind,natpt,nat,name,ok,tblkey, &
                    ntabl,errcnt)
               if(ok) then
                  name=w2
                  !                         lookup-name-into-ind
                  call lookpn(ind,natpt,nat,name,ok,tblkey, &
                       ntabl,errcnt)
                  if(ok) then
                     name=w3
                     !                           lookup-name-into-ind
                     call lookpn(ind,natpt,nat,name,ok,tblkey, &
                          ntabl,errcnt)
                     if(ok) then
                        name=w4
                        !                             lookup-name-into-ind
                        call lookpn(ind,natpt,nat,name,ok, &
                             tblkey,ntabl,errcnt)

                        if(ok) then
                           name=w5
                           !                             lookup-name-into-ind
                           call lookpn(ind,natpt,nat,name,ok, &
                                tblkey,ntabl,errcnt)
                           if(ok) then
                              name=w6
                              !                             lookup-name-into-ind
                              call lookpn(ind,natpt,nat,name,ok, &
                                   tblkey,ntabl,errcnt)
                              if(ok) then
                                 name=w7
                                 !                             lookup-name-into-ind
                                 call lookpn(ind,natpt,nat,name,ok, &
                                      tblkey,ntabl,errcnt)
                                 if(ok) then
                                    name=w8
                                    !                             lookup-name-into-ind
                                    call lookpn(ind,natpt,nat,name,ok, &
                                         tblkey,ntabl,errcnt)

                                    if(ok) then
                                       ok = &
                                            w1 /= w2.and.w1 /= w3.and. &
                                            w1 /= w4.and. &
                                            w2 /= w3.and.w2 /= w4.and. &
                                            w3 /= w4.and. &
                                            w5 /= w6.and.w5 /= w7.and. &
                                            w5 /= w8.and. &
                                            w6 /= w7.and.w6 /= w8.and.w7 /= w8

                                       if (.not.ok) then
                                          call wrndie(2,'<rtfrdr>', &
                                               'The atoms in a cross-term map ' &
                                               //'are not all different.')
                                          errcnt=errcnt+1
                                       else
                                          nict=nict+1
                                          call reallocate_cmap(nict+npct)
                                          delpct(nict+npct)=qdelet
                                          mi1ct(nict+npct)=w1
                                          mj1ct(nict+npct)=w2
                                          mk1ct(nict+npct)=w3
                                          ml1ct(nict+npct)=w4
                                          mi2ct(nict+npct)=w5
                                          mj2ct(nict+npct)=w6
                                          mk2ct(nict+npct)=w7
                                          ml2ct(nict+npct)=w8
                                       endif
                                    endif
                                 endif
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
            w1=nexta6(comlyn,comlen)
         enddo
         call xtrane(comlyn,comlen,'RTF reader')
      endif
      !                 end procedure process-crossterm-map
      return
    end subroutine rtfrdr_process_cmap
#endif


    subroutine rtfrdr_process_donor
      !                 Begin Procedure PROCESS-DONOR
      !
      !                 [DELETE] DONOR [ hydrogen ] heavy-atom
      !                 ['BLNK'    ]
      !
      !
      !                 if 'BLNK' is specified instead of a hydrogen
      !                 an extended donor is assumed.
      !
      !
      if (nrtrs == 0) then
         call wrndie(2,'<rtfrdr>', &
              'A residue must be specified before anything ' &
              //'is put into it.')
         errcnt=errcnt+1
      else
         w1=nexta6(comlyn,comlen)
         if (w1 == 'BLNK') then
            w1=blank
            ok=.true.
            w2=nexta6(comlyn,comlen)
         else
            name=w1
            !                     lookup-name-into-ind
            call lookpn(ind,natpt,nat,name,ok,tblkey,ntabl,errcnt)
            w2=nexta6(comlyn,comlen)
         endif
         if(ok) then
            if (w2 == blank) then
               call wrndie(2,'<rtfrdr>', &
                    'The heavy atom in a hydrogen bond donor ' &
                    //'must be specified.')
               errcnt=errcnt+1
               ok=.false.
            else
               name=w2
               !                       lookup-name-into-ind
               call lookpn(ind,natpt,nat,name,ok,tblkey,ntabl,errcnt)
            endif
            if(ok) then
               nd=nd+1
               call reallocate_rthd(nd+ndpt)
               mh(nd+ndpt)=w1
               deldn(nd+ndpt)=qdelet
               md(nd+ndpt)=w2
               call xtrane(comlyn,comlen,'RTF reader')
            endif
         endif
      endif
      !                 end procedure process-donor
    end subroutine rtfrdr_process_donor


    subroutine rtfrdr_process_acceptor
      !                 begin procedure process-acceptor
      !
      !                 [delete] acceptor iupac [iupac]
      !
      if (nrtrs == 0) then
         call wrndie(2,'<rtfrdr>', &
              'A residue must be specified before anything ' &
              //'is put into it.')
         errcnt=errcnt+1
      else
         w1=nexta6(comlyn,comlen)
         w2=nexta6(comlyn,comlen)
         name=w1
         !                   lookup-name-into-ind
         call lookpn(ind,natpt,nat,name,ok,tblkey,ntabl,errcnt)
         if(ok) then
            name=w2
            !                     lookup-name-into-ind
            call lookpn(ind,natpt,nat,name,ok,tblkey,ntabl,errcnt)
            if(ok) then
               na=na+1
               call reallocate_rtha(na+napt)
               delac(na+napt)=qdelet
               ma(na+napt)=w1
               maa(na+napt)=w2
               call xtrane(comlyn,comlen,'RTF reader')
            endif
         endif
      endif
      !                 end procedure process-acceptor
    end subroutine rtfrdr_process_acceptor


    subroutine rtfrdr_process_ic
      !                 begin procedure process-ic
      !
      !    { build }
      !    [delete] { bild  } name name [*]name name bond angle phi angle bond
      !    { ic    }
      !
      if (nrtrs == 0) then
         call wrndie(2,'<rtfrdr>', &
              'A residue must be specified before anything ' &
              //'is put into it.')
         errcnt=errcnt+1
      else
         nbl=nbl+1
         call reallocate_rtic(nbl+nblpt)
         if(nbl+nblpt > dim_rt(ind_ic)) call wrndie(-5,'<rtfrdr>', &
              'Too many IC.')
         delic(nbl+nblpt)=qdelet
         bari(nbl+nblpt)=nexta6(comlyn,comlen)
         barj(nbl+nblpt)=nexta6(comlyn,comlen)
         w3=nexta6(comlyn,comlen)
         if (w3(1:1) == '*') then
            bart(nbl+nblpt)=.true.
            bark(nbl+nblpt)=w3(2:)
         else
            bart(nbl+nblpt)=.false.
            bark(nbl+nblpt)=w3
         endif
         barl(nbl+nblpt)=nexta6(comlyn,comlen)
         if(bari(nbl+nblpt) == 'BLNK') bari(nbl+nblpt)=blank
         if(barj(nbl+nblpt) == 'BLNK') barj(nbl+nblpt)=blank
         if(bark(nbl+nblpt) == 'BLNK') bark(nbl+nblpt)=blank
         if(barl(nbl+nblpt) == 'BLNK') barl(nbl+nblpt)=blank
         if(.not.qdelet) then
            icb2(nbl+nblpt)=nextf(comlyn,comlen)
            icth2(nbl+nblpt)=nextf(comlyn,comlen)
            icphi(nbl+nblpt)=nextf(comlyn,comlen)
            icth1(nbl+nblpt)=nextf(comlyn,comlen)
            icb1(nbl+nblpt)=nextf(comlyn,comlen)
         endif
         call xtrane(comlyn,comlen,'RTF reader')
      endif
      !                 end procedure process-ic
    end subroutine rtfrdr_process_ic

    subroutine rtfrdr_process_group
      !                 begin procedure process-group
      !
      if (nrtrs == 0) then
         call wrndie(2,'<rtfrdr>', &
              'A residue must be specified before anything ' &
              //'is put into it.')
         errcnt=errcnt+1
      else
         ngr=ngr+1
      endif
      !                 end procedure process-group
    end subroutine rtfrdr_process_group


    subroutine rtfrdr_process_print
      !                 begin procedure process-print
      wrd=nexta4(comlyn,comlen)
      if (wrd == 'ON  ') then
         prnflg=(prnlev >= 2)
      else if (wrd == 'OFF ') then
         prnflg=.false.
      else
         call wrndie(2,'<RTFRDR>','Unrecognized PRINT ' &
              //'option')
         wrncnt=wrncnt+1
      endif
      !                 end procedure process-print
    end subroutine rtfrdr_process_print


    subroutine rtfrdr_process_default_patch
      !                 begin procedure process-default-patch
      !
      !                 DEFAULT [ FIRSt patch-residue ] [ LAST patch-residue ]
      !                 patch-residue:== {   NONE  }
      !                 { name of one PRES }
      !
      wrd=nexta4(comlyn,comlen)
      w2=nexta6(comlyn,comlen)
      if (wrd == 'FIRS') then
         deff=w2
         if(w2 == blank) then
            call wrndie(2,'<RTFRDR>','No default patch ' &
                 //'residue specified')
            wrncnt=wrncnt+1
         endif
      else if (wrd == 'LAST') then
         defl=w2
         if(w2 == blank) then
            call wrndie(2,'<RTFRDR>','No default patch ' &
                 //'residue specified')
            wrncnt=wrncnt+1
         endif
      else
         call wrndie(2,'<RTFRDR>','Illegal default patch ' &
              //'option')
         wrncnt=wrncnt+1
      endif
      wrd=nexta4(comlyn,comlen)
      if(wrd /= blank) then
         w2=nexta6(comlyn,comlen)
         if (wrd == 'FIRS') then
            deff=w2
            if(w2 == blank) then
               call wrndie(2,'<RTFRDR>','No default patch ' &
                    //'residue specified')
               wrncnt=wrncnt+1
            endif
         else if (wrd == 'LAST') then
            defl=w2
            if(w2 == blank) then
               call wrndie(2,'<RTFRDR>','No default patch ' &
                    //'residue specified')
               wrncnt=wrncnt+1
            endif
         else
            call wrndie(2,'<RTFRDR>','Illegal default patch ' &
                 //'option')
            wrncnt=wrncnt+1
         endif
         call xtrane(comlyn,comlen,'RTF reader')
      endif
      !                 end procedure process-default-patch
    end subroutine rtfrdr_process_default_patch


    subroutine rtfrdr_process_auto_command
      !                 begin procedure process-auto-command
      !
      !                 auto [ angle ] [ dihedral ] [ patches ]
      !
      do while (comlen > 0)
         wrd=nexta4(comlyn,comlen)
         call trima(comlyn,comlen)
         if (wrd == 'ANGL') then
            autot=.true.
         else if (wrd == 'NOAN') then
            autot=.false.
         else if (wrd == 'DIHE') then
            autod=.true.
         else if (wrd == 'NODI') then
            autod=.false.
         else if (wrd == 'PATC') then
            autop=.true.
         else if (wrd == 'NOPAT') then
            autop=.false.
         else if (wrd == 'OFF') then
            autot=.false.
            autod=.false.
            autop=.false.
         else if (wrd == 'DRUD') then
            ldrude = .true.
            write(outu,'(/,A)') &
                 ' DRUDES PARTICLES WILL BE GENERATED'// &
                 ' AUTOMATICALLY FOR ALL ATOMS WITH'// &
                 ' NON-ZERO ALPHA'
            call parsethol(tholea)
         else
            call wrndie(2,'<rtfrdr>', &
                 'Illegal auto-gen option')
            wrncnt=wrncnt+1
         endif
      enddo
      call xtrane(comlyn,comlen,'RTF reader')
      !                 end procedure process-auto-command
    end subroutine rtfrdr_process_auto_command


    subroutine rtfrdr_process_patch
      !                 begin procedure process-patch
      !
      !                 PATCH [ FIRSt patch-residue ] [ LAST patch-residue ]
      !                 patch-residue:== {   NONE  }
      !                 { name of one PRES }
      !
      if (nrtrs == 0) then
         call wrndie(2,'<RTFRDR>', &
              'A residue must be specified before anything ' &
              //'is put into it.')
         errcnt=errcnt+1
      else
         wrd=nexta4(comlyn,comlen)
         w2=nexta6(comlyn,comlen)
         if (wrd == 'FIRS') then
            patf=w2
            if(w2 == blank) then
               call wrndie(2,'<RTFRDR>','No patch residue ' &
                    //'was specified')
               wrncnt=wrncnt+1
            endif
         else if (wrd == 'LAST') then
            patl=w2
            if(w2 == blank) then
               call wrndie(2,'<RTFRDR>','No patch residue ' &
                    //'was specified')
               wrncnt=wrncnt+1
            endif
         else
            call wrndie(2,'<RTFRDR>','Illegal patch option')
            wrncnt=wrncnt+1
         endif
         wrd=nexta4(comlyn,comlen)
         if(wrd /= blank) then
            w2=nexta6(comlyn,comlen)
            if (wrd == 'FIRS') then
               patf=w2
               if(w2 == blank) then
                  call wrndie(2,'<RTFRDR>','No patch residue ' &
                       //'was specified')
                  wrncnt=wrncnt+1
               endif
            else if (wrd == 'LAST') then
               patl=w2
               if(w2 == blank) then
                  call wrndie(2,'<RTFRDR>','No patch residue ' &
                       //'was specified')
                  wrncnt=wrncnt+1
               endif
            else
               call wrndie(2,'<RTFRDR>','Illegal patch ' &
                    //'option')
               wrncnt=wrncnt+1
            endif
            call xtrane(comlyn,comlen,'RTF reader')
         endif
      endif
      !                 end procedure process-patch
    end subroutine rtfrdr_process_patch


    subroutine rtfrdr_process_pres
      ipatch=1
      if (.not.qfirst) then
         !                   finish-current-residue
         call fresid(wrncnt,patf,patl,totchg,ind, &
              tblkey,ntabl,errcnt, &
              natpt,nat,nbopt,nbo,nthpt,nth,nphpt,nph, &
              nimpt,nim, &
#if KEY_CMAP==1
              npct,nict,nlpct,nilpt,nanisct, nanis, &
#endif
              nrxpt,nrx,napt,na,ndpt, &
              nd,nblpt,nbl,ngr, &
              ok,newres,idel)
      else
         qfirst=.false.
      endif
      !                 process-residue
      call prores(errcnt,newres,idel,wrncnt,totchg, &
           ipatch,patf,patl, &
           nat,nbo,nth,nph,nim, &
#if KEY_CMAP==1
           nict,nilpt, nanis, &
#endif
           NRX,NA,ND,NBL,NGR)
    end subroutine rtfrdr_process_pres

  end subroutine rtfrdr

  subroutine reallocate_rtrs(n)
    use memory
    integer n,dim_rtrs_new
    character(len=15) :: subname="reallocate_rtrs"
    character(len=9) :: fname="rtfio.src"

    if(.not.allocated(aa))then
       dim_rtrs_new = max(n,dim_rtrs)
       dim_rtrs = 0
       call chmalloc(fname,subname,'nic'   ,dim_nicm,dim_rtrs_new,intg=nic)
       call chmalloc(fname,subname,'aa'    ,dim_rtrs_new,ch6=aa)
       call chmalloc(fname,subname,'aaf'   ,dim_rtrs_new,ch6=aaf)
       call chmalloc(fname,subname,'aal'   ,dim_rtrs_new,ch6=aal)
       call chmalloc(fname,subname,'rtrtyp',dim_rtrs_new,intg=rtrtyp)
    else
       if(n <= size(aa)) return
       dim_rtrs_new = max(n, int(dim_rtrs*factor_rt) )
       call chmrealloc(fname,subname,'nic',   dim_nicm,dim_rtrs_new,intg=nic)
       call chmrealloc(fname,subname,'aa',    dim_rtrs_new,ch6=aa)
       call chmrealloc(fname,subname,'aaf',   dim_rtrs_new,ch6=aaf)
       call chmrealloc(fname,subname,'aal',   dim_rtrs_new,ch6=aal)
       call chmrealloc(fname,subname,'rtrtyp',dim_rtrs_new,intg=rtrtyp)
    endif
    aa    (           dim_rtrs+1:dim_rtrs_new) = '      '
    aaf   (           dim_rtrs+1:dim_rtrs_new) = '      '
    aal   (           dim_rtrs+1:dim_rtrs_new) = '      '
    nic   (1:dim_nicm,dim_rtrs+1:dim_rtrs_new) = 0
    rtrtyp(           dim_rtrs+1:dim_rtrs_new) = 0
    dim_rtrs = dim_rtrs_new
    return
  end subroutine reallocate_rtrs

  subroutine reallocate_rta(n)
    use memory
    use number,only:zero
    character(len=14) :: subname="reallocate_rta"
    character(len=9) :: fname="rtfio.src"

    integer n,dim_rta_new
    if(.not.allocated(ftp))then

       dim_rta_new=max(initial_dim,n)
       dim_rt(ind_at)=0
       call   chmalloc(fname,subname,'ftp'   ,dim_rta_new,ch6=ftp)
#if KEY_MMFF==1
       call   chmalloc(fname,subname,'atnumr',dim_rta_new,intg=atnumr)
#endif
       call   chmalloc(fname,subname,'mac'   ,dim_rta_new,intg=mac)
       call   chmalloc(fname,subname,'mxn'   ,dim_rta_new,intg=mxn)
       call   chmalloc(fname,subname,'grpr'  ,dim_rta_new,intg=grpr)
#if KEY_CHEQ==1
       call   chmalloc(fname,subname,'echh'  ,dim_rta_new,crl=echh)
#endif
#if KEY_CHEQ==1
       call   chmalloc(fname,subname,'ehah'  ,dim_rta_new,crl=ehah)
#endif
       call   chmalloc(fname,subname,'chg'   ,dim_rta_new,crl=chg)
       call   chmalloc(fname,subname,'alph'  ,dim_rta_new,crl=alph)
       call   chmalloc(fname,subname,'thol'  ,dim_rta_new,crl=thol)
       call   chmalloc(fname,subname,'delat' ,dim_rta_new,log=delat)
    else
       if ( n <= size(ftp))return
       dim_rta_new = max(n,int(dim_rt(ind_at)*factor_rt))
       call chmrealloc(fname,subname,'ftp'   ,dim_rta_new,ch6=ftp)
#if KEY_MMFF==1
       call chmrealloc(fname,subname,'atnumr',dim_rta_new,intg=atnumr)
#endif
       call chmrealloc(fname,subname,'mac'   ,dim_rta_new,intg=mac)
       call chmrealloc(fname,subname,'mxn'   ,dim_rta_new,intg=mxn)
       call chmrealloc(fname,subname,'grpr'  ,dim_rta_new,intg=grpr)
#if KEY_CHEQ==1
       call chmrealloc(fname,subname,'echh'  ,dim_rta_new,crl=echh)
#endif
#if KEY_CHEQ==1
       call chmrealloc(fname,subname,'ehah'  ,dim_rta_new,crl=ehah)
#endif
       call chmrealloc(fname,subname,'chg'   ,dim_rta_new,crl=chg)
       call chmrealloc(fname,subname,'alph'  ,dim_rta_new,crl=alph)
       call chmrealloc(fname,subname,'thol'  ,dim_rta_new,crl=thol)
       call chmrealloc(fname,subname,'delat' ,dim_rta_new,log=delat)
    endif
    ftp   (dim_rt(ind_at)+1:dim_rta_new) = '      '
#if KEY_MMFF==1
    atnumr(dim_rt(ind_at)+1:dim_rta_new) = 0
#endif
    mac   (dim_rt(ind_at)+1:dim_rta_new) = 0
    mxn   (dim_rt(ind_at)+1:dim_rta_new) = 0
    grpr  (dim_rt(ind_at)+1:dim_rta_new) = 0
#if KEY_CHEQ==1
    echh  (dim_rt(ind_at)+1:dim_rta_new) = zero
#endif
#if KEY_CHEQ==1
    ehah  (dim_rt(ind_at)+1:dim_rta_new) = zero
#endif
    chg   (dim_rt(ind_at)+1:dim_rta_new) = zero
    alph  (dim_rt(ind_at)+1:dim_rta_new) = zero
    thol  (dim_rt(ind_at)+1:dim_rta_new) = zero
    delat (dim_rt(ind_at)+1:dim_rta_new) = .false.

    dim_rt(ind_at)=dim_rta_new
    return
  end subroutine reallocate_rta

  subroutine reallocate_rtb(n)
    use memory
    character(len=14) :: subname="reallocate_rtb"
    character(len=9) :: fname="rtfio.src"
    integer n,dim_rtb_new

    if(.not.allocated(mib))then
       dim_rtb_new=max(initial_dim,n)
       dim_rt(ind_bo)=0
       call   chmalloc(fname,subname,'mib',   dim_rtb_new,ch6=mib)
       call   chmalloc(fname,subname,'mjb',   dim_rtb_new,ch6=mjb)
       call   chmalloc(fname,subname,'delbd', dim_rtb_new,log=delbd)
#if KEY_MMFF==1
       call   chmalloc(fname,subname,'mbtype',dim_rtb_new,intg=mbtype)
#endif
    else
       if ( n < size(mib))return
       dim_rtb_new = max(n,int(dim_rt(ind_bo) * factor_rt))
       call chmrealloc(fname,subname,'mib',   dim_rtb_new,ch6=mib)
       call chmrealloc(fname,subname,'mjb',   dim_rtb_new,ch6=mjb)
       call chmrealloc(fname,subname,'delbd', dim_rtb_new,log=delbd)
#if KEY_MMFF==1
       call chmrealloc(fname,subname,'mbtype',dim_rtb_new,intg=mbtype)
#endif
    endif
    mib   (dim_rt(ind_bo)+1:dim_rtb_new) = '      '
    mjb   (dim_rt(ind_bo)+1:dim_rtb_new) = '      '
    delbd (dim_rt(ind_bo)+1:dim_rtb_new) = .false.
#if KEY_MMFF==1
    mbtype(dim_rt(ind_bo)+1:dim_rtb_new) = 0
#endif
    dim_rt(ind_bo) = dim_rtb_new
    return
  end subroutine reallocate_rtb

  subroutine reallocate_rtt(n)
    use memory
    character(len=14) :: subname="reallocate_rtt"
    character(len=9) :: fname="rtfio.src"
    integer n,dim_rtt_new

    if(.not.allocated(mit))then
       dim_rtt_new = max(n,initial_dim)
       dim_rt(ind_th) = 0
       call   chmalloc(fname,subname,'mit'  ,dim_rtt_new,ch6=mit)
       call   chmalloc(fname,subname,'mjt'  ,dim_rtt_new,ch6=mjt)
       call   chmalloc(fname,subname,'mkt'  ,dim_rtt_new,ch6=mkt)
       call   chmalloc(fname,subname,'delan',dim_rtt_new,log=delan)
    else
       if ( n < size(mit))return
       dim_rtt_new = max(int(dim_rt(ind_th) * factor_rt), n)
       call chmrealloc(fname,subname,'mit'  ,dim_rtt_new,ch6=mit)
       call chmrealloc(fname,subname,'mjt'  ,dim_rtt_new,ch6=mjt)
       call chmrealloc(fname,subname,'mkt'  ,dim_rtt_new,ch6=mkt)
       call chmrealloc(fname,subname,'delan',dim_rtt_new,log=delan)
    endif
    mit  (dim_rt(ind_th)+1:dim_rtt_new) = '      '
    mjt  (dim_rt(ind_th)+1:dim_rtt_new) = '      '
    mkt  (dim_rt(ind_th)+1:dim_rtt_new) = '      '
    delan(dim_rt(ind_th)+1:dim_rtt_new) = .false.
    dim_rt(ind_th) = dim_rtt_new
    return
  end subroutine reallocate_rtt

  subroutine reallocate_rtp(n)
    use memory
    character(len=14) :: subname="reallocate_rtp"
    character(len=9) :: fname="rtfio.src"
    integer n,dim_rtp_new

    if(.not.allocated(mip))then
       dim_rtp_new = max(n,initial_dim)
       dim_rt(ind_ph)=0
       call   chmalloc(fname,subname,'mip'  ,dim_rtp_new,ch6=mip)
       call   chmalloc(fname,subname,'mjp'  ,dim_rtp_new,ch6=mjp)
       call   chmalloc(fname,subname,'mkp'  ,dim_rtp_new,ch6=mkp)
       call   chmalloc(fname,subname,'mlp'  ,dim_rtp_new,ch6=mlp)
       call   chmalloc(fname,subname,'delpt',dim_rtp_new,log=delpt)
    else
       if ( n < size(mip))return
       dim_rtp_new = max(n, int(dim_rt(ind_ph) * factor_rt) )
       call chmrealloc(fname,subname,'mip'  ,dim_rtp_new,ch6=mip)
       call chmrealloc(fname,subname,'mjp'  ,dim_rtp_new,ch6=mjp)
       call chmrealloc(fname,subname,'mkp'  ,dim_rtp_new,ch6=mkp)
       call chmrealloc(fname,subname,'mlp'  ,dim_rtp_new,ch6=mlp)
       call chmrealloc(fname,subname,'delpt',dim_rtp_new,log=delpt)
    endif
    mip  (dim_rt(ind_ph)+1:dim_rtp_new) = '      '
    mjp  (dim_rt(ind_ph)+1:dim_rtp_new) = '      '
    mkp  (dim_rt(ind_ph)+1:dim_rtp_new) = '      '
    mlp  (dim_rt(ind_ph)+1:dim_rtp_new) = '      '
    delpt(dim_rt(ind_ph)+1:dim_rtp_new) = .false.
    dim_rt(ind_ph) = dim_rtp_new
    return
  end subroutine reallocate_rtp

  subroutine reallocate_rti(n)
    use memory
    character(len=14) :: subname="reallocate_rti"
    character(len=9) :: fname="rtfio.src"
    integer n,dim_rti_new

    if(.not.allocated(mim))then
       dim_rti_new = max(initial_dim,n)
       dim_rt(ind_im) = 0
       call   chmalloc(fname,subname,'mim'  ,dim_rti_new,ch6=mim)
       call   chmalloc(fname,subname,'mjm'  ,dim_rti_new,ch6=mjm)
       call   chmalloc(fname,subname,'mkm'  ,dim_rti_new,ch6=mkm)
       call   chmalloc(fname,subname,'mlm'  ,dim_rti_new,ch6=mlm)
       call   chmalloc(fname,subname,'delmt',dim_rti_new,log=delmt)
    else
       if ( n < size(mim))return
       dim_rti_new = max(n, int(dim_rt(ind_im) * factor_rt) )
       call chmrealloc(fname,subname,'mim'  ,dim_rti_new,ch6=mim)
       call chmrealloc(fname,subname,'mjm'  ,dim_rti_new,ch6=mjm)
       call chmrealloc(fname,subname,'mkm'  ,dim_rti_new,ch6=mkm)
       call chmrealloc(fname,subname,'mlm'  ,dim_rti_new,ch6=mlm)
       call chmrealloc(fname,subname,'delmt',dim_rti_new,log=delmt)
    endif
    mim  (dim_rt(ind_im)+1:dim_rti_new) = '      '
    mjm  (dim_rt(ind_im)+1:dim_rti_new) = '      '
    mkm  (dim_rt(ind_im)+1:dim_rti_new) = '      '
    mlm  (dim_rt(ind_im)+1:dim_rti_new) = '      '
    delmt(dim_rt(ind_im)+1:dim_rti_new) = .false.
    dim_rt(ind_im) = dim_rti_new
    return
  end subroutine reallocate_rti

  subroutine reallocate_rtic(n)
    use memory
    use number,only:zero
    character(len=15) :: subname="reallocate_rtic"
    character(len=9) :: fname="rtfio.src"
    integer n,dim_rtic_new
    if(.not.allocated(bari))then
       dim_rtic_new = max(n,initial_dim)
       dim_rt(ind_ic) = 0
       call chmalloc(  fname,subname,'bari' ,dim_rtic_new,ch6=bari)
       call chmalloc(  fname,subname,'barj' ,dim_rtic_new,ch6=barj)
       call chmalloc(  fname,subname,'bark' ,dim_rtic_new,ch6=bark)
       call chmalloc(  fname,subname,'barl' ,dim_rtic_new,ch6=barl)
       call chmalloc(  fname,subname,'icb1' ,dim_rtic_new,crl=icb1)
       call chmalloc(  fname,subname,'icb2' ,dim_rtic_new,crl=icb2)
       call chmalloc(  fname,subname,'icth1',dim_rtic_new,crl=icth1)
       call chmalloc(  fname,subname,'icth2',dim_rtic_new,crl=icth2)
       call chmalloc(  fname,subname,'icphi',dim_rtic_new,crl=icphi)
       call chmalloc(  fname,subname,'bart' ,dim_rtic_new,log=bart)
       call chmalloc(  fname,subname,'delic',dim_rtic_new,log=delic)
    else
       if ( n <= size(bari))return
       dim_rtic_new = max(n, int(dim_rt(ind_ic) * factor_rt) )
       call chmrealloc(fname,subname,'bari' ,dim_rtic_new,ch6=bari)
       call chmrealloc(fname,subname,'barj' ,dim_rtic_new,ch6=barj)
       call chmrealloc(fname,subname,'bark' ,dim_rtic_new,ch6=bark)
       call chmrealloc(fname,subname,'barl' ,dim_rtic_new,ch6=barl)
       call chmrealloc(fname,subname,'icb1' ,dim_rtic_new,crl=icb1)
       call chmrealloc(fname,subname,'icb2' ,dim_rtic_new,crl=icb2)
       call chmrealloc(fname,subname,'icth1',dim_rtic_new,crl=icth1)
       call chmrealloc(fname,subname,'icth2',dim_rtic_new,crl=icth2)
       call chmrealloc(fname,subname,'icphi',dim_rtic_new,crl=icphi)
       call chmrealloc(fname,subname,'bart' ,dim_rtic_new,log=bart)
       call chmrealloc(fname,subname,'delic',dim_rtic_new,log=delic)
    endif
    bari (dim_rt(ind_ic)+1:dim_rtic_new) = '      '
    barj (dim_rt(ind_ic)+1:dim_rtic_new) = '      '
    bark (dim_rt(ind_ic)+1:dim_rtic_new) = '      '
    barl (dim_rt(ind_ic)+1:dim_rtic_new) = '      '
    icb1 (dim_rt(ind_ic)+1:dim_rtic_new) = zero
    icb2 (dim_rt(ind_ic)+1:dim_rtic_new) = zero
    icth1(dim_rt(ind_ic)+1:dim_rtic_new) = zero
    icth2(dim_rt(ind_ic)+1:dim_rtic_new) = zero
    icphi(dim_rt(ind_ic)+1:dim_rtic_new) = zero
    bart (dim_rt(ind_ic)+1:dim_rtic_new) = .false.
    delic(dim_rt(ind_ic)+1:dim_rtic_new) = .false.
    dim_rt(ind_ic) = dim_rtic_new
    return
  end subroutine reallocate_rtic

  subroutine reallocate_rtx(n)
    use memory
    integer n,dim_rtx_new
    character(len=14) :: subname="reallocate_rtx"
    character(len=9) :: fname="rtfio.src"

    if(.not.allocated(mnb))then
       dim_rtx_new = max(n,initial_dim)
       dim_rt(ind_xx) = 0
       call chmalloc(fname,subname,'mnb'    ,dim_rtx_new,ch6=mnb)
    else
       if(n <= size(mnb)) return
       if(n <= dim_rt(ind_xx))then
          dim_rtx_new = dim_rt(ind_xx)
       else
          dim_rtx_new = int(dim_rt(ind_xx)*factor_rt)
       endif
       call chmrealloc(fname,subname,'mnb'    ,dim_rtx_new,ch6=mnb)
    endif
    mnb(dim_rt(ind_xx)+1:dim_rtx_new) = '      '
    dim_rt(ind_xx) = dim_rtx_new
    if(dim_rt(ind_xx) /= size(mnb)) then
       call wrndie(-5,"<rtfio.src> reallocate_rtx", &
            "size of array mnb does not match saved dimension")
    endif
    return
  end subroutine reallocate_rtx

  subroutine reallocate_rtha(n)
    use memory
    integer n,dim_rtha_new
    character(len=15) :: subname="reallocate_rtha"
    character(len=9) :: fname="rtfio.src"

    if(.not.allocated(ma))then
       dim_rtha_new = max(n,initial_dim)
       dim_rt(ind_ha) = 0
       call chmalloc(fname,subname,'ma'    ,dim_rtha_new,ch6=ma)
       call chmalloc(fname,subname,'maa'   ,dim_rtha_new,ch6=maa)
       call chmalloc(fname,subname,'mab'   ,dim_rtha_new,ch6=mab)
       call chmalloc(fname,subname,'delac' ,dim_rtha_new,log=delac)
    else
       if(n <= size(ma)) return
       dim_rtha_new = max(n, int(dim_rt(ind_ha)*factor_rt))
       call chmrealloc(fname,subname,'ma'    ,dim_rtha_new,ch6=ma)
       call chmrealloc(fname,subname,'maa'   ,dim_rtha_new,ch6=maa)
       call chmrealloc(fname,subname,'mab'   ,dim_rtha_new,ch6=mab)
       call chmrealloc(fname,subname,'delac' ,dim_rtha_new,log=delac)
    endif
    ma (dim_rt(ind_ha)+1:dim_rtha_new) = '      '
    maa(dim_rt(ind_ha)+1:dim_rtha_new) = '      '
    mab(dim_rt(ind_ha)+1:dim_rtha_new) = '      '
    delac(dim_rt(ind_ha)+1:dim_rtha_new) = .false.
    dim_rt(ind_ha) = dim_rtha_new
    return
  end subroutine reallocate_rtha

  subroutine reallocate_rthd(n)
    use memory
    integer n,dim_rthd_new
    character(len=15) :: subname="reallocate_rthd"
    character(len=9) :: fname="rtfio.src"

    if(.not.allocated(md))then
       dim_rthd_new = max(n,initial_dim)
       dim_rt(ind_hd) = 0
       call chmalloc(fname,subname,'md'    ,dim_rthd_new,ch6=md)
       call chmalloc(fname,subname,'mh'    ,dim_rthd_new,ch6=mh)
       call chmalloc(fname,subname,'mda'   ,dim_rthd_new,ch6=mda)
       call chmalloc(fname,subname,'mdb'   ,dim_rthd_new,ch6=mdb)
       call chmalloc(fname,subname,'deldn' ,dim_rthd_new,log=deldn)
    else
       if(n <= size(md)) return
       dim_rthd_new = max(n, int(dim_rt(ind_hd)*factor_rt) )
       call chmrealloc(fname,subname,'md'    ,dim_rthd_new,ch6=md)
       call chmrealloc(fname,subname,'mh'    ,dim_rthd_new,ch6=mh)
       call chmrealloc(fname,subname,'mda'   ,dim_rthd_new,ch6=mda)
       call chmrealloc(fname,subname,'mdb'   ,dim_rthd_new,ch6=mdb)
       call chmrealloc(fname,subname,'deldn' ,dim_rthd_new,log=deldn)
    endif
    md (dim_rt(ind_hd)+1:dim_rthd_new) = '      '
    mh (dim_rt(ind_hd)+1:dim_rthd_new) = '      '
    mda(dim_rt(ind_hd)+1:dim_rthd_new) = '      '
    mdb(dim_rt(ind_hd)+1:dim_rthd_new) = '      '
    deldn(dim_rt(ind_hd)+1:dim_rthd_new) = .false.
    dim_rt(ind_hd) = dim_rthd_new
    return
  end subroutine reallocate_rthd



!------------------------------------------------------------------
!         REALLOCATE_NATCT
!------------------------------------------------------------------
  subroutine reallocate_natct(n)
    use memory
    use number,only:zero
    character(len=14) :: subname="reallocate_natct"
    character(len=9) :: fname="rtfio.src"

    integer n,dim_natct_new
    character(len=4) :: blank='    '

    if(.not.allocated(atct))then
       dim_natct_new=max(initial_dim,n)
       dim_natct=0

       call   chmalloc(fname,subname,'atct'   ,dim_natct_new,ch6=atct)
       call   chmalloc(fname,subname,'armass' ,dim_natct_new,crl=armass)

#if KEY_MMFF==1
       call   chmalloc(fname,subname,'atnumt' ,dim_natct_new,intg=atnumt)
#endif
#if KEY_MMFF==1
       call   chmalloc(fname,subname,'ctype'  ,dim_natct_new,ch4=ctype)
#endif
#if KEY_MMFF==1
       call   chmalloc(fname,subname,'intdef' ,dim_natct_new,intg=intdef)
#endif
    else
       if ( n <= size(atct))return
       dim_natct_new = max(n,int(dim_natct*factor_rt))
       call chmrealloc(fname,subname,'atct'   ,dim_natct_new,ch6=atct)
       call chmrealloc(fname,subname,'armass' ,dim_natct_new,crl=armass)

#if KEY_MMFF==1
       call chmrealloc(fname,subname,'atnumt' ,dim_natct_new,intg=atnumt)
#endif
#if KEY_MMFF==1
       call chmrealloc(fname,subname,'ctype'  ,dim_natct_new,ch4=ctype)
#endif
#if KEY_MMFF==1
       call chmrealloc(fname,subname,'intdef' ,dim_natct_new,intg=intdef)
#endif
    endif

    atct  (dim_natct+1:dim_natct_new) = blank
    armass(dim_natct+1:dim_natct_new) = zero

#if KEY_MMFF==1
    atnumt(dim_natct+1:dim_natct_new) = 0
#endif
#if KEY_MMFF==1
    ctype (dim_natct+1:dim_natct_new) = '    '
#endif
#if KEY_MMFF==1
    intdef(dim_natct+1:dim_natct_new) = 0
#endif


    dim_natct=dim_natct_new

   return
  end subroutine reallocate_natct

#if KEY_CMAP==1 /*cmap*/
  !------------------------------------------------------------------
  !         REALLOCATE_CMAP
  !------------------------------------------------------------------
  subroutine reallocate_cmap(n)
    use memory
    use number,only:zero
    character(len=14) :: subname="reallocate_cmap"
    character(len=9) :: fname="rtfio.src"

    integer n,dim_cmap_new
    character(len=4) :: blank='    '

    if(.not.allocated(mi1ct_test))then
       dim_cmap_new=max(initial_dim,n)
       dim_rt(ind_cmap)=0

       call   chmalloc(fname,subname,'mi1ct_test' ,dim_cmap_new,ch6=mi1ct_test)
       call   chmalloc(fname,subname,'mi1ct' ,dim_cmap_new,ch6=mi1ct)
       call   chmalloc(fname,subname,'mj1ct' ,dim_cmap_new,ch6=mj1ct)
       call   chmalloc(fname,subname,'mk1ct' ,dim_cmap_new,ch6=mk1ct)
       call   chmalloc(fname,subname,'ml1ct' ,dim_cmap_new,ch6=ml1ct)
       call   chmalloc(fname,subname,'mi2ct' ,dim_cmap_new,ch6=mi2ct)
       call   chmalloc(fname,subname,'mj2ct' ,dim_cmap_new,ch6=mj2ct)
       call   chmalloc(fname,subname,'mk2ct' ,dim_cmap_new,ch6=mk2ct)
       call   chmalloc(fname,subname,'ml2ct' ,dim_cmap_new,ch6=ml2ct)
       call   chmalloc(fname,subname,'delpct',dim_cmap_new,log=delpct)
    else
       if ( n <= size(mi1ct_test))return
       dim_cmap_new = max(n,int(dim_rt(ind_cmap)*factor_rt))
       call chmrealloc(fname,subname,'mi1ct_test' ,dim_cmap_new,ch6=mi1ct_test)
       call chmrealloc(fname,subname,'mi1ct' ,dim_cmap_new,ch6=mi1ct)
       call chmrealloc(fname,subname,'mj1ct' ,dim_cmap_new,ch6=mj1ct)
       call chmrealloc(fname,subname,'mk1ct' ,dim_cmap_new,ch6=mk1ct)
       call chmrealloc(fname,subname,'ml1ct' ,dim_cmap_new,ch6=ml1ct)
       call chmrealloc(fname,subname,'mi2ct' ,dim_cmap_new,ch6=mi2ct)
       call chmrealloc(fname,subname,'mj2ct' ,dim_cmap_new,ch6=mj2ct)
       call chmrealloc(fname,subname,'mk2ct' ,dim_cmap_new,ch6=mk2ct)
       call chmrealloc(fname,subname,'ml2ct' ,dim_cmap_new,ch6=ml2ct)
       call chmrealloc(fname,subname,'delpct',dim_cmap_new,log=delpct)
    endif

!--    mi1ct(dim_rt(ind_cmap)+1:dim_cmap_new) = blank
!--    mj1ct(dim_rt(ind_cmap)+1:dim_cmap_new) = blank
!--    mk1ct(dim_rt(ind_cmap)+1:dim_cmap_new) = blank
!--    ml1ct(dim_rt(ind_cmap)+1:dim_cmap_new) = blank
!--    mi2ct(dim_rt(ind_cmap)+1:dim_cmap_new) = blank
!--    mj2ct(dim_rt(ind_cmap)+1:dim_cmap_new) = blank
!--    mk2ct(dim_rt(ind_cmap)+1:dim_cmap_new) = blank
!--    ml2ct(dim_rt(ind_cmap)+1:dim_cmap_new) = blank

    dim_rt(ind_cmap)=dim_cmap_new

    return
  end subroutine reallocate_cmap

  !------------------------------------------------------------------
  !         REALLOCATE_LP
  !------------------------------------------------------------------
  subroutine reallocate_lp(n)
    use memory
    use number,only:zero
    character(len=14) :: subname="reallocate_lp"
    character(len=9) :: fname="rtfio.src"

    integer n,dim_lp_new
    character(len=4) :: blank='    '
    if(.not.allocated(dellpct))then
       dim_lp_new=max(initial_dim,n)
       dim_rt(ind_lp)=0
       call   chmalloc(fname,subname,'mlp0ct' ,dim_lp_new,ch4=mlp0ct)
       call   chmalloc(fname,subname,'mlp1ct' ,dim_lp_new,ch6=mlp1ct)
       call   chmalloc(fname,subname,'mlp2ct' ,dim_lp_new,ch6=mlp2ct)
       call   chmalloc(fname,subname,'mlp3ct' ,dim_lp_new,ch6=mlp3ct)
       call   chmalloc(fname,subname,'mlp4ct' ,dim_lp_new,ch6=mlp4ct)
    !   call   chmalloc(fname,subname,'mlpct'  ,dim_lp_new,maxcent_hosts,ch6=mlpct)
       call   chmalloc(fname,subname,'mlpct1d',dim_lp_new*maxcent_hosts,ch6=mlpct1d)
       call   chmalloc(fname,subname,'ncenthst',dim_lp_new,intg=ncenthst)
       call   chmalloc(fname,subname,'rlpct'  ,dim_lp_new,crl=rlpct)
       call   chmalloc(fname,subname,'rlpct'  ,dim_lp_new,crl=tlpct)
       call   chmalloc(fname,subname,'rlpct'  ,dim_lp_new,crl=plpct)
       call   chmalloc(fname,subname,'dellpct',dim_lp_new,log=dellpct)
    else
       if ( n <= size(dellpct))return
       dim_lp_new = max(n,int(dim_rt(ind_lp)*factor_rt))
       call chmrealloc(fname,subname,'mlp0ct' ,dim_lp_new,ch4=mlp0ct)
       call chmrealloc(fname,subname,'mlp1ct' ,dim_lp_new,ch6=mlp1ct)
       call chmrealloc(fname,subname,'mlp2ct' ,dim_lp_new,ch6=mlp2ct)
       call chmrealloc(fname,subname,'mlp3ct' ,dim_lp_new,ch6=mlp3ct)
       call chmrealloc(fname,subname,'mlp4ct' ,dim_lp_new,ch6=mlp4ct)
   !    call chmrealloc(fname,subname,'mlpct'  ,dim_lp_new,maxcent_hosts,ch6=mlpct)
       call chmrealloc(fname,subname,'mlpct1d',dim_lp_new*maxcent_hosts,ch6=mlpct1d)
       call chmrealloc(fname,subname,'ncenthst' ,dim_lp_new,intg=ncenthst)
       call chmrealloc(fname,subname,'rlpct'  ,dim_lp_new,crl=rlpct)
       call chmrealloc(fname,subname,'rlpct'  ,dim_lp_new,crl=tlpct)
       call chmrealloc(fname,subname,'rlpct'  ,dim_lp_new,crl=plpct)
       call chmrealloc(fname,subname,'dellpct',dim_lp_new,log=dellpct)
    endif

!--    mi1ct(dim_rt(ind_lp)+1:dim_lp_new) = blank
!--    mj1ct(dim_rt(ind_lp)+1:dim_lp_new) = blank
!--    mk1ct(dim_rt(ind_lp)+1:dim_lp_new) = blank
!--    ml1ct(dim_rt(ind_lp)+1:dim_lp_new) = blank
!--    mi2ct(dim_rt(ind_lp)+1:dim_lp_new) = blank
!--    mj2ct(dim_rt(ind_lp)+1:dim_lp_new) = blank
!--    mk2ct(dim_rt(ind_lp)+1:dim_lp_new) = blank
!--    ml2ct(dim_rt(ind_lp)+1:dim_lp_new) = blank

    dim_rt(ind_lp)=dim_lp_new

    return
  end subroutine reallocate_lp
#endif /* (cmap)*/

  !------------------------------------------------------------------
  !         REALLOCATE_ANIS
  !------------------------------------------------------------------
  subroutine reallocate_anis(n)
    use memory
    use number,only:zero
    character(len=14) :: subname="reallocate_anis"
    character(len=9) :: fname="rtfio.src"

    integer n,dim_anis_new
    character(len=4) :: blank='    '

    if(.not.allocated(delanis))then
       dim_anis_new=max(initial_dim,n)
       dim_rt(ind_anis)=0
       call   chmalloc(fname,subname,'mianis' ,dim_anis_new,ch6=mianis)
       call   chmalloc(fname,subname,'mjanis' ,dim_anis_new,ch6=mjanis)
       call   chmalloc(fname,subname,'mkanis' ,dim_anis_new,ch6=mkanis)
       call   chmalloc(fname,subname,'mlanis' ,dim_anis_new,ch6=mlanis)
       call   chmalloc(fname,subname,'a11anis'  ,dim_anis_new,crl=a11anis)
       call   chmalloc(fname,subname,'a22anis'  ,dim_anis_new,crl=a22anis)
       call   chmalloc(fname,subname,'delanis',dim_anis_new,log=delanis)
    else
       if ( n <= size(delanis))return
       dim_anis_new = max(n,int(dim_rt(ind_anis)*factor_rt))
       call chmrealloc(fname,subname,'mianis' ,dim_anis_new,ch6=mianis)
       call chmrealloc(fname,subname,'mjanis' ,dim_anis_new,ch6=mjanis)
       call chmrealloc(fname,subname,'mkanis' ,dim_anis_new,ch6=mkanis)
       call chmrealloc(fname,subname,'mlanis' ,dim_anis_new,ch6=mlanis)
       call chmrealloc(fname,subname,'a11anis'  ,dim_anis_new,crl=a11anis)
       call chmrealloc(fname,subname,'a22anis'  ,dim_anis_new,crl=a22anis)
       call chmrealloc(fname,subname,'delanis',dim_anis_new,log=delanis)
    endif

    dim_rt(ind_anis)=dim_anis_new

    return
  end subroutine reallocate_anis

!------------------------------------------------------------------
!         Parse Thol
!------------------------------------------------------------------
  SUBROUTINE PARSETHOL(THOLEA)
    use exfunc,only:srchws
    use number,only:zero,one
    use stream,only:outu
    use string
    use comand
    use psf,only:qthole,tholes

    implicit none
    real(chm_real) tholea,tdef

    tdef = 1.3*one
    qthole = .false.
    tholea = gtrmf(comlyn,comlen,'THOL',tdef)
    tholea = abs(tholea)
    tholes = gtrmi(comlyn,comlen,'SHAP',1)

    if (tholea > zero)then
       qthole = .true.

       if (tholes == 1)then
          write(outu,'(1x,3a,f10.6,1x,a,/)')  &
               'Thole-type dipole screening, ', &
               'Slater-Delta shape {S(u) = 1 - (1+u/2)*exp(-u)}, ', &
               'default radius =',THOLEA
       else if (tholes == 2)then
          write(outu,'(1x,3a,f10.6,1x,a,/)')  &
               'Thole-type dipole screening, ', &
               'Gaussian shape {S(u) = erf(u)}, ', &
               'default radius =',THOLEA
       else
          call wrndie(-3,'<DRUDE>','INVALID SHAPE')
       endif

    endif

    return
  end subroutine parsethol


  subroutine lookpn(ind,natpt,nat,name,ok,tblkey,ntabl,errcnt)
    !     to lookup-name-into-ind
    use stream,only:outu
    use exfunc,only:srchws
    implicit none
    integer ind,natpt,nat,ntabl,errcnt
    character(len=*) name,tblkey(*)
    logical ok

    !     local
    integer :: iok=999
    save iok
    ind=srchws(ftp(natpt+1),nat,name)
    if (ind > 0) then
       ok=.true.
    else
       ind=srchws(tblkey,ntabl,name)
       if (ind == 0) then
          if (rtrtyp(nrtrs) == 0) then
             if(wrnlev >= 2) write(outu,410) name
410          FORMAT(' **** ERROR in RTFRDR ****  The atom name, ',A6, &
                  ' was not found in this residue or in the DECLARE table.' &
                  /' It''s use will be ignored.')
             errcnt=errcnt+1
             ok=.false.
          else
             ind=iok
             iok=iok+1
             ok=.true.
          endif
       else
          ok=.true.
       endif
    endif
    !
    return
  end subroutine lookpn

  subroutine prores(errcnt,newres,idel,wrncnt,totchg, &
       ipatch,patf,patl, &
       nat,nbo,nth,nph,nim, &
#if KEY_CMAP==1
       nict,nilpt,nanis, &
#endif
       nrx,na,nd,nbl,ngr)

    !     TO PROCESS-RESIDUE

    !
    use comand
    use stream
    use string

    implicit none
    character(len=6) patf,patl
    integer     errcnt,wrncnt
    integer     idel,ipatch,nrx,nat,nbo,nth,nph,nim,nd,na,ngr,nbl
#if KEY_CMAP==1
    integer     nict
    integer     nilpt
    integer     nanis
#endif
    logical     newres
    real(chm_real)      totchg

    !     local
    character(len=6) :: blank='      ',name


    !     RESIDUE name [total-charge]

    name=nexta6(comlyn,comlen)
    if (name == blank) then
       call wrndie(2,'<RTFRDR>','Residue must have a name')
       errcnt=errcnt+1
    endif

    newres=.true.
    idel=0
    do while (newres .and. idel < nrtrs)
       idel=idel+1
       if(aa(idel) == name) then
          newres=.false.
          if(wrnlev >= 2) write(outu,341) aa(idel)
341       FORMAT('*** WARNING **** residue ',A4,' already exists ', &
               '(old one deleted)')
          wrncnt=wrncnt+1
          CALL WRNDIE(1,'<RTFRDR>', 'Residue already exists.')
       endif
    enddo
    nrtrs=nrtrs+1
    call reallocate_rtrs(nrtrs)

    aa(nrtrs)=name

    totchg=nextf(comlyn,comlen)
    call xtrane(comlyn,comlen,'RTF reader')
    rtrtyp(nrtrs)=ipatch
    patf=deff
    patl=defl
    nat=0
    nbo=0
    nth=0
    nph=0
    nim=0
#if KEY_CMAP==1
    nict=0
    nilpt=0
    nanis=0
#endif
    nrx=0
    na=0
    nd=0
    nbl=0
    ngr=0

    return
  end subroutine prores

  subroutine fresid(wrncnt,patf,patl,totchg,ind, &
       tblkey,ntabl,errcnt, &
       natpt,nat,nbopt,nbo,nthpt,nth,nphpt,nph, &
       nimpt,nim, &
#if KEY_CMAP==1
       npct,nict,nlpct,nilpt,nanisct, nanis, &
#endif
       nrxpt,nrx,napt,na,ndpt,nd,nblpt,nbl, &
       ngr,ok,newres,idel)

    !     to finish-current-residue
    use comand
    use stream
    implicit none
    character(len=*) patf,patl,tblkey(*)
    integer errcnt,wrncnt,ind,ntabl
    integer natpt,nat,nbopt,nbo,nthpt,nth,nphpt,nph,idel
    integer nimpt,nim,nrxpt,nrx,napt,na,ndpt,nd,nblpt,nbl,ngr
#if KEY_CMAP==1
    integer npct,nict
    integer nlpct,nilpt
    integer nanisct, nanis
#endif
    logical ok,newres
    real(chm_real)  totchg

    !     local
    character(len=6) name
    integer iat,inx,start,stop
    real(chm_real)  sum
    sum=0.0
    do iat=1,nat
       if(.not.delat(iat+natpt)) then
          sum=sum+chg(iat+natpt)
       endif
    enddo
    if(abs(sum-totchg) > 1.0e-4) then
       if(wrnlev >= 2) write(outu,430) aa(nrtrs),sum,totchg
       wrncnt=wrncnt+1
    endif
430 FORMAT(' **** WARNING from RTFRDR **** The total charge of the', &
         ' residue, ',A4,', ',F11.7,',' &
         /' does not equal the expected charge, ',F11.7,'.')
    aaf(nrtrs)=patf
    aal(nrtrs)=patl

    !     Here we convert the non-bonded exclusion lists stored temporarily
    iat=0
    do while (iat < nat)
       iat=iat+1
       if (iat+natpt == 1) then
          start=1
       else
          start=mxn(iat+natpt-1)+1
       endif
       stop=mxn(iat+natpt)
       inx=start
       do while (inx <= stop)
          name=mnb(inx)
          !         lookup-name-into-ind
          call lookpn(ind,natpt,nat,name,ok,tblkey,ntabl,errcnt)
          inx=inx+1
       enddo
    enddo
    !
    natpt=natpt+nat
    nbopt=nbopt+nbo
    nthpt=nthpt+nth
    nphpt=nphpt+nph
    nimpt=nimpt+nim
#if KEY_CMAP==1
    npct=npct+nict
    nlpct=nlpct+nilpt
    nanisct=nanisct+nanis
#endif
    nrxpt=nrxpt+nrx
    napt=napt+na
    ndpt=ndpt+nd
    nblpt=nblpt+nbl
    nic(ind_at,nrtrs)=natpt
    nic(ind_bo,nrtrs)=nbopt
    nic(ind_th,nrtrs)=nthpt
    nic(ind_ph,nrtrs)=nphpt
    nic(ind_im,nrtrs)=nimpt

#if KEY_CMAP==1
    nic(11,nrtrs)=npct
    nic(ind_lp,nrtrs)=nlpct
    nic(ind_anis,nrtrs)=nanisct
#endif
    nic(ind_xx,nrtrs)=nrxpt
    nic(ind_hd,nrtrs)=ndpt
    nic(ind_ha,nrtrs)=napt
    nic(ind_ic,nrtrs)=nblpt
    nic(10,nrtrs)=ngr
    if(.not.newres) then
       call delete(idel)
       natpt=nic(ind_at,nrtrs)
       nbopt=nic(ind_bo,nrtrs)
       nthpt=nic(ind_th,nrtrs)
       nphpt=nic(ind_ph,nrtrs)
       nimpt=nic(ind_im,nrtrs)
#if KEY_CMAP==1
       npct=nic(11,nrtrs)
       nlpct=nic(ind_lp,nrtrs)
       nanisct=nic(ind_anis,nrtrs)
#endif
       nrxpt=nic(ind_xx,nrtrs)
       ndpt =nic(ind_hd,nrtrs)
       napt =nic(ind_ha,nrtrs)
       nblpt=nic(ind_ic,nrtrs)
       ngr=nic(10,nrtrs)
    endif
    ok = natpt <= dim_rt(ind_at) .and. &
         nbopt <= dim_rt(ind_bo) .and. &
         nthpt <= dim_rt(ind_th) .and. &
         nphpt <= dim_rt(ind_ph) .and. &
         nimpt <= dim_rt(ind_im) .and. &
         nrxpt <= dim_rt(ind_xx) .and. &
#if KEY_CMAP==1
         npct  <= dim_rt(ind_cmap) .and. &
#endif
         napt  <= dim_rt(ind_ha) .and. &
         ndpt  <= dim_rt(ind_hd) .and. &
         nblpt <= dim_rt(ind_ic)
    if (.not.(ok)) then
       if(wrnlev >= 2) write(outu,590) &
            natpt, dim_rt(ind_at), &
            nbopt, dim_rt(ind_bo), &
            nthpt, dim_rt(ind_th), &
            nphpt, dim_rt(ind_ph), &
            nimpt, dim_rt(ind_im), &
#if KEY_CMAP==1
            npct,  dim_rt(ind_xx), &
#endif
#if KEY_CMAP==1
            nrxpt, dim_rt(ind_cmap), &
#endif
            napt,  dim_rt(ind_ha), &
            ndpt,  dim_rt(ind_hd), &
            nblpt, dim_rt(ind_ic)
590    FORMAT(13X,'IN USE   MAXIMUM',/ &
            ' ATOMS      ',2I8,/ &
            ' BONDS      ',2I8,/ &
            ' ANGLES     ',2I8,/ &
            ' DIHEDRALS  ',2I8,/ &
            ' IMPROPERS  ',2I8,/ &
#if KEY_CMAP==1
            ' CROSSTERMS ',2I8,/ &
#endif
            ' EXCLUSIONS ',2I8,/ &
            ' ACCEPTORS  ',2I8,/ &
            ' DONORS     ',2I8,/ &
            ' BUILD/IC   ',2I8)
       call wrndie(-4,'<RTFRDR>','LIMIT EXCEEDED')
    endif
    return
  end subroutine fresid

  subroutine rtwrit(title,ntitl,unit,iprint,len)
    !
    !     THIS ROUTINE WRITES THE RESIDUE TOPOLGY INFORMATION TO DISK
    !     OR TO OUTPUT.  WRITES OUT BOND, ANGLES, DIHEDRALS, AND IMPROPERS
    !     BY CODE NUMBER AND ATOM NAME.
    !     24-May-1981  DJS
    !     26-Oct-1981  19:42:42 BRUC
    !     25-Mar-1982  BRB
    !     21-JAN-1983  AB
    !
    !     The indices into the NIC arrays are as follows:
    !     IAT = 1, IBO = 2, ITH = 3, IPH = 4, IIM = 5,
    !     IX = 6,  ID  = 7, IA  = 8, IBILD=9, IGRP=10
    !

    use stream
    use number,only:zero
    implicit none
    !
    character(len=*) title(*)
    integer ntitl,iprint,unit,len,numbe
    integer iout(len)
    real(chm_real)  out(len)
    character(len=8) cout(len)
    integer i,ires,kxstrt,k,kx,j
    integer natpt,nbopt,nthpt,nphpt,nimpt,nrxpt,napt,ndpt,nblpt, &
         ngrppt
    integer nat,nbo,nth,nph,nim,nrx,na,nd,nbl,ngr
#if KEY_CMAP==1
    integer npct,nict
    integer nlpct,nilpt
    integer nanisct, nanis
#endif
    real(chm_real)  cg, alpht
    character(len=4) :: hdr='RTF ',blank='    ',junk='JUNK'

    if(iolev < 0) return

    if (iprint == 0) then
       !       begin procedure write-binary-topology-file
       rtcntl(1)=20
       write(unit) hdr,rtcntl
       call wrtitl(title,ntitl,unit,-1)

       !       set unused arrays
       do i=1,nic(ind_hd,nrtrs)
          mda(i)=junk
       enddo
       do i=1,nic(ind_hd,nrtrs)
          mdb(i)=junk
       enddo
       do i=1,nic(ind_ha,nrtrs)
          mab(i)=junk
       enddo
       rtftyp=0
       if(autot) rtftyp=1
       if(autod) rtftyp=rtftyp+2
       if(autop) rtftyp=rtftyp+4
       !
       write(unit) nrtrs,rtftyp,dim_nicm,natct
       write(unit) ((nic(i,ires),i=1,dim_nicm),aa(ires),rtrtyp(ires), &
            aaf(ires),aal(ires),ires=1,nrtrs)
       write(unit) (atct(i),armass(i),i=1,natct)
       write(unit) (chg(i),ftp(i),mac(i),mxn(i),grpr(i),delat(i), &
            i=1,nic(ind_at,nrtrs))
       write(unit) (mib(i),mjb(i),delbd(i),i=1,nic(ind_bo,nrtrs))
       write(unit) (mit(i),mjt(i),mkt(i),delan(i),i=1,nic(ind_th,nrtrs))
       write(unit) (mip(i),mjp(i),mkp(i),mlp(i),delpt(i), &
            i=1,nic(ind_ph,nrtrs))
       write(unit) (mim(i),mjm(i),mkm(i),mlm(i),delmt(i), &
            i=1,nic(ind_im,nrtrs))
       write(unit) (mnb(i),i=1,nic(ind_xx,nrtrs))
       write(unit) (md(i),mh(i),mda(i),mdb(i),deldn(i), &
            i=1,nic(ind_hd,nrtrs))
       write(unit) (ma(i),maa(i),mab(i),delac(i),i=1,nic(ind_ha,nrtrs))
       write(unit) (icb2(i),icb1(i),icth2(i),icth1(i),icphi(i),bari(i), &
            barj(i),bark(i),barl(i),bart(i),delic(i), &
            i=1,nic(ind_ic,nrtrs))
       !       end procedure write-binary-topology-file
    else
       !       begin procedure print-topology-file
       !
       natpt=0
       nbopt=0
       nthpt=0
       nphpt=0
       nimpt=0
#if KEY_CMAP==1
       npct=0
       nlpct=0
       nanisct=0
#endif
       nrxpt=0
       napt=0
       ndpt=0
       nblpt=0
       ngrppt=0
       kxstrt=1

       write(unit,35) rtcntl
35     FORMAT(/10X,'RESIDUE TOPOLOGY FILE MODULE' &
            /10X,'CONTROL ARRAY :',20I5)
       call wrtitl(title,ntitl,unit,1)
       write(unit,30) nrtrs
30     FORMAT(/10X,'  Number of residues is ',I3)


       write(unit,39)
39     FORMAT(/10X,'ATOM TYPE CODES : MASSES')
       do k=1,natct
          if(atct(k) /= blank) &
               write(unit,'(4X,I5,A10,F16.4)') k,atct(k),armass(k)
       enddo
       do i=1,nrtrs
          write(unit,'(//,80("-"))')
          write(unit,'("1",10X,"- ",I2," -",4X,A4/)') i, aa(i)
          if(rtrtyp(i) > 0) write(unit,'(" * NOTE *  THIS IS A PATCH"/)')

          nat=nic(ind_at,i)-natpt
          nbo=nic(ind_bo,i)-nbopt
          nth=nic(ind_th,i)-nthpt
          nph=nic(ind_ph,i)-nphpt
          nim=nic(ind_im,i)-nimpt
#if KEY_CMAP==1
          nict=nic(11,i)-npct
          nilpt=nic(ind_lp,i)-nlpct
          nanis=nic(ind_anis,i)-nanisct
#endif
          nrx=nic(ind_xx,i)-nrxpt
          nd=nic(ind_hd,i)-ndpt
          na=nic(ind_ha,i)-napt
          nbl=nic(ind_ic,i)-nblpt
          ngr=nic(10,i)-ngrppt

          write(unit,70) nat,(k,k=1,nat)
70        FORMAT (1X,I2,' ATOMS, WITH FOLLOWING ATTRIBUTES:'/ &
               (1X,'ATOM NUMBER', 17I9)/(12X,17I7))
          write(unit,90) (ftp(k+natpt),k=1,nat)
90        FORMAT ((1X,'ATOM NAME  ',17(3X,A6))/(12X,17(3X,A6)))
          do k=1,nat
             if (delat(k+natpt)) then
                cout(k)='dele'
             else
                cout(k)=atct(mac(k+natpt))
             endif
          enddo
          write(unit,'((1X,"ATOM TYPE  ",17(3X,A6))/(12X,17(3X,A6)))') &
               (cout(k),k=1,nat)
          do k=1,nat
             if (delat(k+natpt)) then
                iout(k)=0
             else
                iout(k)=mac(k+natpt)
             endif
          enddo
          write(unit,'((1X,"TYPE CODE  ",17(I7,2X))/(12X,17(I5,2X)))') &
               (iout(k),k=1,nat)

          do k=1,nat
             if (delat(k+natpt)) then
                iout(k)=0
             else
                iout(k)=grpr(k+natpt)
             endif
          enddo
          write(unit,'((1X,"GROUP NO   ",17(I7,2X))/(12X,17(I5,2X)))' ) &
               (iout(k),k=1,nat)

          do k=1,nat
             if (delat(k+natpt)) then
                out(k)=0.0
             else
                out(k)=chg(k+natpt)
             endif
          enddo
          write(unit,'((1X,"CHARGE     ",17F9.3)/(12X,17F7.3))') &
               (out(k),k=1,nat)

          do k=1,nat
             if (delat(k+natpt)) then
                out(k)=0.0
             else
                out(k)=alph(k+natpt)
             endif
          enddo
          write(unit,'((1X,"ALPHA      ",17F9.3)/(12X,17F7.3))') &
               (out(k),k=1,nat)

          do k=1,nat
             if (delat(k+natpt)) then
                out(k)=0.0
             else
                out(k)=thol(k+natpt)
             endif
          enddo
          write(unit,'((1X,"THOLE      ",17F9.3)/(12X,17F7.3))') &
               (out(k),k=1,nat)

          cg=zero
          do k=1,nat
             if(.not.delat(k+natpt)) then
                cg=cg + chg(k+natpt)
             endif
          enddo
          write(unit,'(/20X,"TOTAL CHARGE = ", F7.4)') cg

          alpht=0.0
          do k=1,nat
             if(.not.delat(k+natpt)) then
                alpht=alpht + alph(k+natpt)
             endif
          enddo
          write(unit,'(20X,"TOTAL ALPHA  = ", F7.4)') alpht


          write(unit, &
               '(/1X,I2," BONDS:",10(1X,A6," : ",A6)/(10X,10(1X,A6," : ",A6)))') &
               nbo, (mib(k+nbopt),mjb(k+nbopt),k=1,nbo)

          write(unit,56) nth,(mit(k+nthpt),mjt(k+nthpt),mkt(k+nthpt), &
               k=1,nth)
56        FORMAT(/1X,I2,' ANGLES:',5(4X,A6,' : ',A6,' : ',A6)/ &
               (11X,5(4X,A6,' : ',A6,' : ',A6)))
          write(unit,61) nph,(mip(k+nphpt),mjp(k+nphpt),mkp(k+nphpt), &
               mlp(k+nphpt),k=1,nph)
61        FORMAT(/1X,I2,' DIHEDRALS:',4(3X,A6,' : ',A6,' : ',A6,' : ', &
               A6)/(14X,4(3X,A6,' : ',A6,' : ',A6,' : ',A6)))
          write(unit,71) nim,(mim(k+nimpt),mjm(k+nimpt),mkm(k+nimpt), &
               mlm(k+nimpt),k=1,nim)
71        FORMAT(/1X,I2,' IMPROPERS:',4(3X,A6,' : ',A6,' : ',A6,' : ', &
               A6)/(14X,4(3X,A6,' : ',A6,' : ',A6,' : ',A6)))
#if KEY_CMAP==1
          write(unit,72) nict,(mi1ct(k+npct),mj1ct(k+npct), &
               mk1ct(k+npct),ml1ct(k+npct), &
               mi2ct(k+npct),mj2ct(k+npct), &
               mk2ct(k+npct),ml2ct(k+npct), &
               k=1,nict)
72        FORMAT(/1X,I2,' CROSS-TERM:',4(3X,A6,' : ',A6,' : ',A6,' : ', &
               A6)/(14X,4(3X,A6,' : ',A6,' : ',A6,' : ',A6)))

!          write(unit,73) nilpt,(mlp0ct(k+nlpct), &
!               mlp1ct(k+nlpct),mlp2ct(k+nlpct), &
!               mlp3ct(k+nlpct),mlp4ct(k+nlpct), &
!               K=1,NILPT)
!73        FORMAT(/1X,I2,' LONEPAIRS:',8(3X,(A6,' { ',4(A6,2x),'}')))
           DO K=1,NILPT
              IF (MLP0CT(K+NLPCT).NE."CENT") THEN
                write(unit,73) nilpt,mlp0ct(k+nlpct), &
                mlp1ct(k+nlpct),mlp2ct(k+nlpct), &
                mlp3ct(k+nlpct),mlp4ct(k+nlpct)
              ELSE
                NUMBE=ncenthst(k+nlpct)
                write(unit,*) nilpt,"LONEPAIRS: ",mlp0ct(k+nlpct),"{", (mlpct1d((k+nlpct-1)*maxcent_hosts+j),j=1,numbe),"}"
              ENDIF
           END DO

73        FORMAT(/1X,I2,' LONEPAIRS:',8(3X,(A4,' { ',4(A6,2x),'}')))

          write(unit,74) nanis, (mianis(k+nanisct), &
               mjanis(k+nanisct),mkanis(k+nanisct), &
               mlanis(k+nanisct),k=1,nanis)
74        FORMAT(/1X,I2,' ANISOTROPY:',8(4x,4(A,2x)))

#endif
          write(unit,101) nd,(mh(k+ndpt),md(k+ndpt),k=1,nd)
101       FORMAT(/1X,I2,' HB DONORS:',4(3X,A6,' : ',A6)/ &
               (14X,4(3X,A6,' : ',A6)))
          write(unit,106) na,(ma(k+napt),maa(k+napt),k=1,na)
106       FORMAT(/1X,I2,' HB ACCEPTORS:',4(3X,A6,' : ',A6)/ &
               (17X,4(3X,A6,' : ',A6)))
          write(unit,'(/a)') ' NONBOND EXCLUSIONS BY ATOM:'
          do k=1,nat
             kx=mxn(k+natpt)-kxstrt+1
             write(unit, &
                  '(1X,I2," EXCLUSIONS FROM ",A6,4X,15(2X,A6),10(/28X,15(2X,A6)))') &
                  kx,ftp(k+natpt),(mnb(kx), &
                  kx=kxstrt,mxn(k+natpt))
             kxstrt=mxn(k+natpt)+1
          enddo
          write(unit,150) (bari(k+nblpt),barj(k+nblpt),bark(k+nblpt), &
               barl(k+nblpt),icb2(k+nblpt),icth2(k+nblpt),icphi(k+nblpt), &
               icth1(k+nblpt),icb1(k+nblpt),bart(k+nblpt),k=1,nbl)
150       FORMAT(/1X,'RESIDUE IC INFORMATION'/ &
               '  I    J    K    L      BOND IJ    ANGLE IJK  DIHEDRAL', &
               ' ANGLE JKL    BOND KL  IMP'/(1X,4A5,F11.5,3F11.3, &
               F11.5,L3))

          natpt=natpt+nat
          nbopt=nbopt+nbo
          nthpt=nthpt+nth
          nphpt=nphpt+nph
          nimpt=nimpt+nim
#if KEY_CMAP==1
          npct=npct+nict
          nlpct=nlpct+nilpt
          nanisct=nanisct+nanis
#endif
          nrxpt=nrxpt+nrx
          napt=napt+na
          ndpt=ndpt+nd
          nblpt=nblpt+nbl
          ngrppt=ngrppt+ngr
       enddo
       write(unit,590) natpt,nbopt,nthpt,nphpt,nimpt, &
#if KEY_CMAP==1
            npct, nlpct, nanisct, &
#endif
            nrxpt,napt,ndpt,nblpt,ngrppt
590    FORMAT('1',10X,'TOTALS IN USE',/ &
            ' ATOMS      ',I5,/ &
            ' BONDS      ',I5,/ &
            ' ANGLES     ',I5,/ &
            ' DIHEDRALS  ',I5,/ &
            ' IMPROPERS  ',I5,/ &
#if KEY_CMAP==1
            ' CROSSTERMS ',I5,/ &
#endif
#if KEY_CMAP==1
            ' LONE PAIRS ',I5,/ &
#endif
#if KEY_CMAP==1
            ' ANISOTROP  ',I5,/ &
#endif
            ' EXCLUSIONS ',I5,/ &
            ' ACCEPTORS  ',I5,/ &
            ' DONORS     ',I5,/ &
            ' BUILD/IC   ',I5,/ &
            ' GROUPS     ',I5)
       !       End Procedure PRINT-TOPOLOGY-FILE
    endif
    return
  end subroutine rtwrit

  subroutine delete(ires)
    use stream
    implicit none

    integer ires

    !     local
    integer first,i,iptot,j,last,shift,numbe
    !
    !     This subroutine deletes residue name from the list
    !     and trims all associated tables
    !
    !.....delete residue name                : AA
    !.....delete residue type                : RTRTYP
    !.....delete residue first patch         : AAF
    !.....delete residue last  patch         : AAL
    !.....delete (1) atoms                   : CHG,FTP,MAC,MXN,GRPR,DELAT,ALPH,THOL  
    !.....delete (2) bonds                   : MIB,MJB,DELBD
    !.....delete (3) angles                  : MIT,MJT,MKT,DELAN
    !.....delete (4) dihedrals               : MIP,MJP,MKP,MLP,DELPT
    !.....delete (5) impropers               : MIM,MJM,MKM,MLM,DELMT
#if KEY_CMAP==1
    !.....delete (11) cross-term map         : MI1CT,MJ1CT,MK1CT,ML1CT,
    !                                          MI2CT,MJ2CT,MK2CT,ML2CT,
    !                                          DELPCT
#endif
    !.....delete (6) exclusions              : MNB
    !.....delete (7) hydrogen bond donors    : MD,MH,MDA,MDB,DELDN
    !.....delete (8) hydrogen bond acceptors : MA,MAA,MAB,DELAC
    !.....delete (9) internal coordinates    : ICB2,ICB1,ICTH2,ICTH1,ICPHI,B
    !     BARJ,BARK,BARL,BART,DELIC
    !
    !.....delete residue name, type,   first patch, last patch
    !     AA    RTFTYP  AAF          AAL
    !
    do i=ires,nrtrs-1
       aa(i)=aa(i+1)
       rtrtyp(i)=rtrtyp(i+1)
       aaf(i)=aaf(i+1)
       aal(i)=aal(i+1)
    enddo

    !.....delete (1) atoms                   : CHG,FTP,MAC,MXN,GRPR,DELAT
    call shiftx(1,dim_nicm,ires,nrtrs,nic,first,last,shift)

    iptot=nic(ind_xx,ires)
    if(ires > 1) iptot=nic(ind_xx,ires)-nic(ind_xx,ires-1)
    do i=first,last
       chg(i-shift)=chg(i)
       alph(i-shift)=alph(i)
       thol(i-shift)=thol(i)
       ftp(i-shift)=ftp(i)
       mac(i-shift)=mac(i)
       mxn(i-shift)=mxn(i)-iptot
       grpr(i-shift)=grpr(i)
       delat(i-shift)=delat(i)
    enddo

    !.....delete (2) bonds                   : MIB,MJB,DELBD
    call shiftx(2,dim_nicm,ires,nrtrs,nic,first,last,shift)
    do i=first,last
       mib(i-shift)=mib(i)
       mjb(i-shift)=mjb(i)
       delbd(i-shift)=delbd(i)
    enddo

    !.....delete (3) angles                  : MIT,MJT,MKT,DELAN
    call shiftx(3,dim_nicm,ires,nrtrs,nic,first,last,shift)
    do i=first,last
       mit(i-shift)=mit(i)
       mjt(i-shift)=mjt(i)
       mkt(i-shift)=mkt(i)
       delan(i-shift)=delan(i)
    enddo

    !.....delete (4) dihedrals               : MIP,MJP,MKP,MLP,DELPT
    call shiftx(4,dim_nicm,ires,nrtrs,nic,first,last,shift)
    do i=first,last
       mip(i-shift)=mip(i)
       mjp(i-shift)=mjp(i)
       mkp(i-shift)=mkp(i)
       mlp(i-shift)=mlp(i)
       delpt(i-shift)=delpt(i)
    enddo

    !.....delete (5) impropers               : MIM,MJM,MKM,MLM,DELMT
    call shiftx(5,dim_nicm,ires,nrtrs,nic,first,last,shift)
    do i=first,last
       mim(i-shift)=mim(i)
       mjm(i-shift)=mjm(i)
       mkm(i-shift)=mkm(i)
       mlm(i-shift)=mlm(i)
       delmt(i-shift)=delmt(i)
    enddo

#if KEY_CMAP==1
    !
    !.....delete (11) cross-term maps
    !
    call shiftx(11,dim_nicm,ires,nrtrs,nic,first,last,shift)
    do i=first,last
       mi1ct(i-shift)=mi1ct(i)
       mj1ct(i-shift)=mj1ct(i)
       mk1ct(i-shift)=mk1ct(i)
       ml1ct(i-shift)=ml1ct(i)
       mi2ct(i-shift)=mi2ct(i)
       mj2ct(i-shift)=mj2ct(i)
       mk2ct(i-shift)=mk2ct(i)
       ml2ct(i-shift)=ml2ct(i)
       delpct(i-shift)=delpct(i)
    enddo

    !.....delete (12) lone-pairs
    call shiftx(12,dim_nicm,ires,nrtrs,nic,first,last,shift)
    do i=first,last
       mlp0ct(i-shift)=mlp0ct(i)
       mlp1ct(i-shift)=mlp1ct(i)
       mlp2ct(i-shift)=mlp2ct(i)
       mlp3ct(i-shift)=mlp3ct(i)
       mlp4ct(i-shift)=mlp4ct(i)
       numbe=ncenthst(i)
       do j= 1,numbe
          mlpct1d((i-shift-1)*maxcent_hosts+j)=mlpct1d((i-1)*maxcent_hosts+j)
       end do
       dellpct(i-shift)=dellpct(i)
       rlpct(i-shift)= rlpct(i)
       tlpct(i-shift)=tlpct(i)
       plpct(i-shift)=plpct(i)
    enddo
#endif

    !.....delete (6) exclusions              : MNB
    call shiftx(6,dim_nicm,ires,nrtrs,nic,first,last,shift)
    do i=first,last
       mnb(i-shift)=mnb(i)
    enddo

    !.....delete (7) hydrogen bond donors    : MD,MH,MDA,MDB,DELDN
    call shiftx(7,dim_nicm,ires,nrtrs,nic,first,last,shift)
    do i=first,last
       md(i-shift)=md(i)
       mh(i-shift)=mh(i)
       mda(i-shift)=mda(i)
       mdb(i-shift)=mdb(i)
       deldn(i-shift)=deldn(i)
    enddo

    !.....delete (8) hydrogen bond acceptors : MA,MAA,MAB,DELAC
    call shiftx(8,dim_nicm,ires,nrtrs,nic,first,last,shift)
    do i=first,last
       ma(i-shift)=ma(i)
       maa(i-shift)=maa(i)
       mab(i-shift)=mab(i)
       delac(i-shift)=delac(i)
    enddo

    !.....delete (9) internal coordinates    : ICB2,ICB1,ICTH2,ICTH1,ICPHI,B
    !     BARJ,BARK,BARL,BART,DELIC
    call shiftx(9,dim_nicm,ires,nrtrs,nic,first,last,shift)
    do i=first,last
       icb1(i-shift)=icb1(i)
       icb2(i-shift)=icb2(i)
       icth1(i-shift)=icth1(i)
       icth2(i-shift)=icth2(i)
       icphi(i-shift)=icphi(i)
       bari(i-shift)=bari(i)
       barj(i-shift)=barj(i)
       bark(i-shift)=bark(i)
       barl(i-shift)=barl(i)
       bart(i-shift)=bart(i)
       delic(i-shift)=delic(i)
    enddo

    !.....finally shift NIC array
    do j=1,dim_nicm
       shift=nic(j,ires)
       if(ires > 1) shift=shift-nic(j,ires-1)
       !
       !       nic(dim_nicm,*) is not cumulative therefore shift is 0 for j == dim_nicm
       !
       if(j == dim_nicm) shift=0
       do i=ires,nrtrs-1
          nic(j,i)=nic(j,i+1)-shift
       enddo
       nic(j,nrtrs)=0
    enddo
    nrtrs=nrtrs-1
    return
  end subroutine delete


  subroutine shiftx(iprop,ndim,ires,nrtrs,nic,first,last,shift)
    integer iprop,ndim,ires,nrtrs,first,last,shift
    integer nic(ndim,*)
    first=nic(iprop,ires)+1
    last =nic(iprop,nrtrs)
    if(ires == 1) shift=first-1
    if(ires > 1) shift=first-nic(iprop,ires-1)-1
    return
  end subroutine shiftx

end module rtf
