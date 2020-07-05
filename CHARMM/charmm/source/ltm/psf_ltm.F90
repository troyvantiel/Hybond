module psf
  use chm_kinds
  use dimens_fcm
  character(len=*),private,parameter :: file_name   ="psf_ltm.src"
  !CHARMM Element source/fcm/psf.fcm 1.1
  !
  !     The Protein Structure File (PSF)
  !
  !     Purpose:
  !
  !     To store the lists of bonds, angles, torsions, and improper
  !     torsions for the entire structure, as well as information needed
  !     to generate the non-bonded and hydrogen bond lists.
  !
  !     I/O:    Unformatted I/O and Print: PSFIO in [MK.PROT]PSFRES.FOR
  !
  !     Notes:    The index entry gives the domain of the data.
  !
  !     Variable  Index    Purpose
  !
  !     NATOM              Number of atoms
  !     NRES               Number of residues
  !     NSEG               Number of segments
  !     NBOND              Number of bonds
  !     NTHETA             Number of angles
  !     NPHI               Number of dihedrals
  !     NIMPHI             Number of improper dihedrals
#if KEY_CMAP==1
  !     NCRTERM            Number of cross-term maps             
#endif
  !     NNB                Number of non-bonded exclusions
  !     NDON               Number of donors
  !     NACC               Number of acceptors
  !     NGRP               Number of groups
  !     NST2               Number of ST2 waters
  !....................................................................
  !     Polarizable force field with classical Drude oscillators
  !     QDRUDE             .TRUE. if there are some Drude oscillators
  !     ISDRUDE   Atom     .TRUE. if this atom is a Drude oscillator
  !     NDRUDE             Number of drude particles
  !     NBDRUDE            Pointer to the list of drude bonds (soon to be obsolete)
  !     QTHOLE             .TRUE. if the dipoles have Thole-type interactions
  !     THOLES             Shape of the dipole screening function
  !     ALPHADP   Atom     ALPHADP(I)=alpha 
  !     THOLEI    Atom     thole radius parameters
  !....................................................................
  !     AMASS     Atom     Mass of each atom
  !     CG        Atom     Charge of each atom
  !     IAC       Atom     Parameter type code
  !     IBLO      Atom     Points to INB giving last non-bonded exclusion
  !                        for this atom
  !     IMOVE     Atom     Flag indicating whether this atom can move.
  !                            1 = fixed atom
  !                            0 = mobile
  !                           -1 = a lonepair (may move, but no degfs)
  !     ATYPE     Atom     IUPAC Name for each atom
  !     RSCLF     Atom     Radius Scale Factor for nonbonded (vdw)
  !
  !     IB        Bond     First atom of a bond
  !     JB        Bond     Second atom of bond
  !
  !     IT        Angle    First atom of bond angle
  !     JT        Angle    Second atom of bond angle
  !     KT        Angle    Third atom of bond angle
  !
  !     IP        Phi      First atom of torsion
  !     JP        Phi      Second atom of torsion angle
  !     KP        Phi      Third atom of torsion angle
  !     JP        Phi      Fourth atom of torsion angle
  !
  !     IM        Imphi    First atom of improper torsion
  !     JM        Imphi    Second atom of improper torsion
  !     KM        Imphi    Third atom of improper torsion
  !     LM        Imphi    Fourth atom of improper torsion
#if KEY_CMAP==1
  !     I1CT      Cmap     First atom of first torsion
  !     J1CT      Cmap     Second atom of first torsion
  !     K1CT      Cmap     Third atom of first torsion
  !     L1CT      Cmap     Fourth atom of first torsion
  !     I2CT      Cmap     First atom of second torsion
  !     J2CT      Cmap     Second atom of second torsion
  !     K2CT      Cmap     Third atom of second torsion
  !     L2CT      Cmap     Fourth atom of second torsion
#endif 
  !
  !     IAC1      Acceptor First antecedent to a hydrogen bond acceptor
  !     IACC      Acceptor Hydrogen bond acceptor
  !
  !     IDON      Donor    Heavy atom of a hydrogen bond donor
  !     IHD1      Donor    Hydrogen of a hydrogen bond donor. Zero if
  !                        hydrogens are not used.
  !     INB       NBex     Non-bonded exclusions, indexed by IBLO
  !
  !     IMOVEG    Groups   Used if entire group has fixed atoms
  !     IGPBS     Groups   Base pointer to first atom in each group.
  !     IGPTYP    Groups   Coded type of each group
  !                         (0-no charges,1-neutral group,2-charged,3-ST2)
  !
  !     IBASE     Residue  IBASE(1)=0, IBASE(IRES+1) gives last atom of
  !                        IRESth residue
  !     RES       Residue  Residue Name (e.g. TYR, GLY, ADE,..)
  !     RESID     Residue  Residue Identifier (e.g. 1, 23, 45B,...)
  !     SEGID     Segment  Segment Identifier (e.g. A1, MAIN,...)
  !     NICTOT    Segment  Partition of structure in segments.
  !                        segments where NICTOT(1)=0 and
  !                        NICTOT(NSEG+1) gives NRES.
  !
  !     Both of the following data structures (Z14, Z14G) are built in
  !     NBEXCL.SRC and used in NBNDGC.SRC.
  !
  !     Other Parameters in the Common block
  !     Note added later:  These parameters are defined in DIMENS.FCM,
  !     actually. 
  !
  !     MAXA       Maximum number of atoms
  !     MAXB       Maximum number of bonds
  !     MAXIMP     Maximum number of improper torsions
#if KEY_CMAP==1
  !     MAXCRT     Maximum number of cross terms
#endif 
  !     MAXNB      Maximum number of non-bonded exclusions
  !     MAXP       Maximum number of torsions
  !     MAXPAD     Maximum number of hydrogen bond donors or acceptors
  !     MAXRES     Maximum number of residues
  !     MAXSEG     Maximum number of segments
  !     MAXGRP     Maximum number of groups
  !     MAXT       Maximum number of bond angles
  !

  !
  !====================================================================
  !             Free Format Version
  !====================================================================
  ! integers

  !--mfc--       INTEGER NATOM, NRES, NSEG, NBOND, NTHETA, NPHI, NIMPHI, &
  !--mfc--              NGRP, NATOMT, NREST, NSEGT, NBONDT, NTHETT, NPHIT, &
  !--mfc--              NIMPHT, NGRPT, NNB, NDON, NACC, NST2, NICTOT, IBASE, &
  !--mfc-- ##IF CMAP
  !--mfc--              NCRTERM,NCRTT, &
  !--mfc-- ##ENDIF
  !--mfc--              IGPBS, IGPTYP, IMOVEG, IAC, IMOVE, IB, JB, &
  !--mfc--              IT, JT, KT, IP, JP, KP, LP, IM, JM, KM, LM, &
  !--mfc-- ##IF CMAP
  !--mfc--              I1CT,J1CT,K1CT,L1CT,I2CT,J2CT,K2CT,L2CT, &
  !--mfc-- ##ENDIF
  !--mfc--              IDON, IHD1, IACC, IAC1, INB, IBLO, &
  !--mfc--              NDRUDE, NBDRUDE, THOLES

  integer :: NATOM, NRES, NSEG, NBOND, NTHETA, NPHI, &
       NIMPHI, NGRP, NATOMT, NREST, NSEGT, NBONDT, &
#if KEY_CMAP==1
       NCRTERM,NCRTT, &     
#endif
       NTHETT, NPHIT, NIMPHT, NGRPT, NNB, NDON, NACC, NST2
  integer,allocatable,dimension(:) :: NICTOT, IBASE,IGPBS, IGPTYP, &
       IMOVEG, IAC, IMOVE, &
       IB, JB, IT, JT, KT, &
       IP, JP, KP, LP, IM, JM, KM, LM, IDON, IHD1, &
#if KEY_CMAP==1
       I1CT,J1CT,K1CT,L1CT, &
       I2CT,J2CT,K2CT,L2CT, &
#endif 
       IACC, IAC1, INB, IBLO

  integer,save:: NDRUDE, NBDRUDE,THOLES

  !--mfc-- 
  !--mfc-- ! Weeks Chandler Anderson decomposition scalar
#if KEY_WCA==1
  !--mfc--       real(chm_real) WCA                                
#endif
  !--mfc-- 
  !--mfc-- 

  real(chm_real),allocatable,dimension(:) :: CG, AMASS, RSCLF, ALPHADP, THOLEI
#if KEY_WCA==1
  real(chm_real),allocatable,dimension(:) :: WCA  
#endif
  real(chm_real) :: CGTOT

  ! characters
  character(len=8),allocatable,dimension(:) :: SEGID, RES, RESID, ATYPE

  logical,allocatable,dimension(:) :: ISDRUDE
  logical :: QDRUDE, QTHOLE
  logical QHYPER ! hyperpolarizability
  logical QHARDWALL ! hard wall constraint on drude bond length
  INTEGER HORDER
  REAL(chm_real)  KHYPER, RHYPER, L_WALL
  !
contains

  subroutine psf_iniall()
    implicit none
    natom=0
    nres=0
    nseg=0
    nbond=0
    ntheta=0
    nphi=0
    nimphi=0
#if KEY_CMAP==1
    ncrterm=0
#endif 
    ngrp=0
    nnb=0
    ndon=0
    nacc=0
    nst2=0
    ! BEGIN DRUDE (G. Lamoureux, E. Harder)
    ndrude=0
    nbdrude=0
    qdrude=.false.
    QHARDWALL=.false.
    ! END DRUDE (G. Lamoureux)
    return
  end subroutine psf_iniall

  subroutine allocate_psf_ltm()
    use memory
    character(len=*),parameter :: routine_name="allocate_psf_ltm"
    call chmalloc(file_name,routine_name,'nictot',maxseg+1,intg=nictot)
    NICTOT(1)=0
    call chmalloc(file_name,routine_name,'ibase ',maxres+1,intg=ibase)
    call chmalloc(file_name,routine_name,'igpbs ',maxgrp+1,intg=igpbs)
    call chmalloc(file_name,routine_name,'igptyp',maxgrp  ,intg=igptyp)
    call chmalloc(file_name,routine_name,'imoveg',maxgrp  ,intg=imoveg)
    call chmalloc(file_name,routine_name,'iac   ',maxaim  ,intg=iac)
    iac=0
    call chmalloc(file_name,routine_name,'imove ',maxaim  ,intg=imove)
    imove=0
    call chmalloc(file_name,routine_name,'ib    ',maxb    ,intg=ib)
    call chmalloc(file_name,routine_name,'jb    ',maxb    ,intg=jb)
    call chmalloc(file_name,routine_name,'it    ',maxt    ,intg=it)
    call chmalloc(file_name,routine_name,'jt    ',maxt    ,intg=jt)
    call chmalloc(file_name,routine_name,'kt    ',maxt    ,intg=kt)
    call chmalloc(file_name,routine_name,'ip    ',maxp    ,intg=ip)
    call chmalloc(file_name,routine_name,'jp    ',maxp    ,intg=jp)
    call chmalloc(file_name,routine_name,'kp    ',maxp    ,intg=kp)
    call chmalloc(file_name,routine_name,'lp    ',maxp    ,intg=lp)
    call chmalloc(file_name,routine_name,'im    ',maximp  ,intg=im)
    call chmalloc(file_name,routine_name,'jm    ',maximp  ,intg=jm)
    call chmalloc(file_name,routine_name,'km    ',maximp  ,intg=km)
    call chmalloc(file_name,routine_name,'lm    ',maximp  ,intg=lm)
    call chmalloc(file_name,routine_name,'idon  ',maxpad  ,intg=idon)
    call chmalloc(file_name,routine_name,'ihd1  ',maxpad  ,intg=ihd1)
#if KEY_CMAP==1
    call chmalloc(file_name,routine_name,'i1ct  ',maxcrt  ,intg=i1ct)       
#endif
#if KEY_CMAP==1
    call chmalloc(file_name,routine_name,'j1ct  ',maxcrt  ,intg=j1ct)       
#endif
#if KEY_CMAP==1
    call chmalloc(file_name,routine_name,'k1ct  ',maxcrt  ,intg=k1ct)       
#endif
#if KEY_CMAP==1
    call chmalloc(file_name,routine_name,'l1ct  ',maxcrt  ,intg=l1ct)       
#endif
#if KEY_CMAP==1
    call chmalloc(file_name,routine_name,'i2ct  ',maxcrt  ,intg=i2ct)       
#endif
#if KEY_CMAP==1
    call chmalloc(file_name,routine_name,'j2ct  ',maxcrt  ,intg=j2ct)       
#endif
#if KEY_CMAP==1
    call chmalloc(file_name,routine_name,'k2ct  ',maxcrt  ,intg=k2ct)       
#endif
#if KEY_CMAP==1
    call chmalloc(file_name,routine_name,'l2ct  ',maxcrt  ,intg=l2ct)       
#endif
    call chmalloc(file_name,routine_name,'iacc  ',maxpad  ,intg=iacc)
    call chmalloc(file_name,routine_name,'iac1  ',maxpad  ,intg=iac1)
    call chmalloc(file_name,routine_name,'inb   ',maxnb   ,intg=inb)
    call chmalloc(file_name,routine_name,'iblo  ',maxaim  ,intg=iblo)
    call chmalloc(file_name,routine_name,'cg    ',maxaim  ,crl=cg)
    call chmalloc(file_name,routine_name,'amass ',maxaim  ,crl=amass)
    call chmalloc(file_name,routine_name,'rsclf ',maxaim  ,crl=rsclf)
    call chmalloc(file_name,routine_name,'alphadp',maxaim  ,crl=alphadp)
    call chmalloc(file_name,routine_name,'tholei',maxaim  ,crl=tholei)
#if KEY_WCA==1
    call chmalloc(file_name,routine_name,'wca   ',maxaim  ,crl=wca)          
#endif
    call chmalloc(file_name,routine_name,'segid ',maxseg  ,ch8=segid)
    call chmalloc(file_name,routine_name,'res   ',maxres  ,ch8=res)
    call chmalloc(file_name,routine_name,'resid ',maxres  ,ch8=resid)
    call chmalloc(file_name,routine_name,'atype ',maxaim  ,ch8=atype)
    call chmalloc(file_name,routine_name,'isdrude',maxaim  ,log=isdrude)

    ! initialize some module variables
    isdrude=.false.
    iblo = 0
    return
  end subroutine allocate_psf_ltm

  logical function fixed_atoms() result(LFIXED)
    
    LFIXED = any(IMOVE /= 0)
    
  end function fixed_atoms

  logical function pair_fixed_atoms(i,j) result(LFIXED)
    integer, intent(in) :: i, j
    LFIXED = (IMOVE(i) /= 0 .or. IMOVE(j) /=0)
    
  end function pair_fixed_atoms

end module psf

