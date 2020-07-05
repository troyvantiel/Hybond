SUBROUTINE psf_read_formatted(U,TITLE,NTITLE,QAPPE,qxplor)
  !
  ! new card format
  ! ===============
  ! This routine reads the protein structure file module (formatted file)
  ! Axel Brunger, 13-AUG-84
  ! =======================
  !

#if KEY_CHEQ==1
  use cheq,only:ech,eha,molt,molnt          
#endif

  use chm_kinds
  use dimens_fcm
  !     input/output
  use psf
  use param,only:natc,atc
  use lonepr
  use stream
  use string
#if KEY_MMFF==1 || KEY_CFF==1
  use ffieldm
  use mmffm
#endif 
#if KEY_CFF==1
  use cff_fcm
#endif 
  use exfunc
  use aniso_fcm
  use comand
  use parallel,only:mynod
  !
#if KEY_STRINGM==1 /*   VO string method */
    use machio, only: ifreeu
    use multicom_aux
#endif
  !
  implicit none
#if KEY_CHEQ==1
  INTEGER MOLNTX
#else /**/
  integer molntx_dummy,molnt_dummy,molt_dummy
  real(chm_real) ech_dummy,eha_dummy
#endif 
  logical,intent(inout) :: qxplor
  !
  character(len=*) title(*)
  integer ntitle
  integer u
  logical qappe
  ! local
  logical eof
  integer i, ii, j, iseg, ires, iatmtype
  integer natomx,nbondx,nthetx,nphix, &
       nimphx,nnbx,ndonx,naccx,ngrpx,nst2x, &
       nat,numlpx,numlphx,nanisox


  character(len=4) hdr
  character(len=8) lsegid, lresid, lres, ltype
  character(len=8) osegid, oresid, ores

#if KEY_CMAP==1
  integer ncrtx             
#endif

  character(len=MXCMSZ) line
  integer linelen
  logical qcmap_psf
  logical qcheq_psf
  logical qdrude_psf
  !
!MFC-- Needs to go away if separate io works.
!-- ##IF ENSEMBLE
!--   integer flen
!--   integer, parameter :: mxfile=128
!--   character(len=mxfile) junknm
!--   logical qopen,qform,qwrite,qens
!-- ##ENDIF
  !
#if KEY_STRINGM==1 /*  VO stringm v */
  integer :: oldiol
  logical :: qstr
  common /replicaio/ qstr ! need global variable
#endif
  !
  !yw  24-Jan-2003 expansion code
  character(len=10) :: fmt00,fmt03,fmt04,fmt05
  character(len=20) :: fmt06,fmt07
  character(len=60) :: fmt01,fmt01a,fmt01b,fmt02,fmt02a,fmt02b
  character(len=6) atc_x
  ! begin
  !
!MFC-- Needs to go away if separate io works.
!-- ##IF ENSEMBLE
!--   !     check whether it is serial or parallel read 
!--   !     for serial, root will read one file and broadcast to all
!--   !     for parallel, each node reads its own file
!--   call ensinq('unit',junknm,mxfile,flen,qopen, &
!--        qform,qwrite,qens,u)
!--   if(iolev > 0) then
!--      if (qens) then
!--         write(outu,*) ' psfrd2>  reading separate psfs'
!--      else
!--         write(outu,*) ' psfrd2>  reading single psf'
!--      endif
!--   endif
!-- ##ENDIF
#if KEY_STRINGM==1 /*  VO stringm */
  qstr=u.le.size(ifreeu) ! protect against OOB on ifreeu
  if (qstr) qstr=((MOD(IFREEU(U),8).eq.0).and.(IFREEU(U).ne.0).and.(ME_LOCAL.eq.0))
  if (qstr) then
   oldiol=iolev
   iolev=1
  endif
#endif
  !
  if(qxplor) then
     CALL WRNDIE(1,'<psf_read_formatted>', &
          'Attempting to read XPLOR psf format ')
  endif

  masterio: IF( (IOLEV > 0)  )THEN

     EOF=.FALSE.
     !
     !       read fixed-field formatted file
     !
     CALL RDCMND(line,MXCMSZ,lineLEN,u,EOF,.FALSE.,.FALSE.,' ')
     !
#if KEY_STRINGM==1 /*  VO : restore iolev */
     if (qstr) iolev=oldiol 
#endif
     !
     hdr=nexta4(line,linelen)
     if (hdr /= 'PSF') &
          call wrndie(-1,'<PSFRD2>','Headers dont match (PSF missing)')
     qextfmt=indxa(line,linelen,'EXT') > 0
     qcmap_psf = (INDXA(line, linelen, 'CMAP')  >  0)
     qcheq_psf = (INDXA(line, linelen, 'CHEQ')  >  0) 
     qdrude_psf = (INDXA(line, linelen, 'DRUDE')  >  0)
     qxplor = qxplor .or.  (INDXA(line, linelen, 'XPLOR')  >  0)
     if (qdrude_psf) then
        QDRUDE = .TRUE.
#if KEY_LONEPAIR==1
        if (.not. allocated(lpnhost)) call allocate_lonepair  
#endif
        if (.not. allocated(lstani1)) call allocate_aniso
     endif

#if KEY_CHEQ==0
     if(qcheq_psf)then
        write(outu,'(a/a)') &
             "<PSFRDR2> This psf has CHEQ information ", &
             "          Charmm not compiled with CHEQ "
        !           CALL WRNDIE(-4,'<PSFRD2)','Cannot read psf.')
     endif
#endif 

     fmt00(1:10) = ' '
     fmt01(1:60) = ' '
     fmt01a(1:60)= ' '
     fmt01b(1:60)= ' '
     fmt02(1:60) = ' '
     fmt02a(1:60)= ' '
     fmt02b(1:60)= ' '
     fmt03(1:10) = ' '
     fmt04(1:10) = ' '
     fmt05(1:10) = ' '
     fmt06(1:20) = ' '
     if (qextfmt) then
        idleng=8
        if (prnlev >= 2) write(outu,'(a)') &
             ' PSFRD2> Reading PSF in the expanded format.'
        fmt00='(/I10)'
        fmt01=  '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)'
        fmt01a= '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6)'
        fmt01b= '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6,L1)'
        fmt02=  '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6)'
        fmt02a= '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6)'
        fmt02b= '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6,l1)'
        fmt03='(8I10)'
        fmt04='(9I10)'
        fmt05='(/2I10)'
        fmt06='(2I10,3X,L1,3G14.6)'
        fmt07='(10X,3G14.6)'
     else
        idleng=4
        fmt00='(/I8)'
        fmt01='(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)'
        fmt01a= &
             '(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6)'
        fmt02='(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)'
        fmt02a= &
             '(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6)'
        fmt03='(8I8)'
        fmt04='(9I8)'
        fmt05='(/2I8)'
        fmt06='(2I8,3X,L1,3G14.6)'
     endif

     read(u,fmt00,end=6) ntitle
     read(u,'(a)',end=6) (title(j),j=1,ntitle)
     call wrtitl(title,ntitle,outu,+1)
     read(u,fmt00,end=6) natomx
     !     
     if(.not.qappe) then
        nseg=0
        nres=0
        natom=0
     endif
     nat=natom
     i=natom+1
     natom=natom+natomx
     !
     xplorfmt:IF(qxplor)THEN
        cheqon:IF(qcheq_psf)THEN
#if KEY_CHEQ==1
           if(qdrude_psf)then
              read(u,fmt02b,end=6) &
                   ii,lsegid,lresid,lres,atype(i),atc_x, &
                   cg(i),amass(i),imove(i),alphadp(i),tholei(i)
              if((amass(i) < 1.0).and.(amass(i) > 0.0001))then
                 isdrude(i)=.true.
                 ndrude=ndrude+1
              endif
           else
              read(u,fmt02a,end=6) &
                   ii,lsegid,lresid,lres,atype(i),atc_x, &
                   cg(i),amass(i),imove(i),ech(i),eha(i)
           endif
#else /**/
           read(u,fmt02a,end=6) &
                ii,lsegid,lresid,lres,atype(i),atc_x, &
                cg(i),amass(i),imove(i),ech_dummy,eha_dummy
#endif 
        else !cheq
           read(u,fmt02,end=6) &
                ii,lsegid,lresid,lres,atype(i),atc_x, &
                cg(i),amass(i),imove(i)
        endif cheqon
        iac(i) = -1
        iatmtype = srchws(atc,natc,atc_x)
        if(iatmtype>0) iac(i) = iatmtype
        if(iac(i) == -1) then
           write(outu,'("Atom ",i10,i10," ",a," problem in psf")') i, iatmtype, atc_x
           CALL WRNDIE(-3,'<psf_read_formatted>', &
                'atom type not found for atom')
        endif
     else ! not xplorfmt
        cheqon2: IF(qcheq_psf)THEN
#if KEY_CHEQ==1
           if(qdrude_psf)then
              read(u,fmt01b,end=6) &
                   ii,lsegid,lresid,lres,atype(i),iac(i), &
                   cg(i),amass(i),imove(i),alphadp(i),tholei(i),isdrude(i)
              if(isdrude(i)) ndrude=ndrude+1
           else
              read(u,fmt01a,end=6) &
                   ii,lsegid,lresid,lres,atype(i),iac(i), &
                   cg(i),amass(i),imove(i),ech(i),eha(i)
           endif
#else /**/
           read(u,fmt01a,end=6) &
                ii,lsegid,lresid,lres,atype(i),iac(i), &
                cg(i),amass(i),imove(i),ech_dummy,eha_dummy
#endif 
        else   !cheq
           read(u,fmt01,end=6) &
                ii,lsegid,lresid,lres,atype(i),iac(i), &
                cg(i),amass(i),imove(i)
        endif cheqon2
     endif xplorfmt
     rsclf(i)=1.0
#if KEY_WCA==1
     wca(i)=1.0                                
#endif
     !RCZ 92/06/17 - convert to upper case
     call cnvtuc(lsegid,idleng)
     call cnvtuc(lresid,idleng)
     call cnvtuc(lres,idleng)
     call cnvtuc(atype(i),idleng)
     !
     !       generate NSEG, NICTOT, LSEGID and NRES, IBASE, LRES, LRESID
     !       linked lists
     nseg=nseg+1
     nictot(nseg)=nres
     segid(nseg)=lsegid
     nres=nres+1
     ibase(nres)=nat
     resid(nres)=lresid
     res(nres)=lres
     osegid=lsegid
     oresid=lresid
     ores=lres
     ndrude=0
     isdrude(nat+1)= .false.

     !-------------------------------------------------
     !      Atom Loop
     !-------------------------------------------------
     atomsloop: do i=nat+2,natom
        isdrude(i)= .false.  

        xplorfmt2:IF(qxplor)THEN
           cheqon3:IF(qcheq_psf)THEN
#if KEY_CHEQ==1
              if(qdrude_psf)then
                 read(u,fmt02b,end=6) &
                      ii,lsegid,lresid,lres,atype(i),atc_x, &
                      cg(i),amass(i),imove(i), &
                      alphadp(i),tholei(i)
                 if((amass(i) < 1.0).and.(amass(i) > 0.0001))then
                    isdrude(i)=.true.
                    ndrude=ndrude+1
                 endif
              else
                 read(u,fmt02a,end=6) &
                      ii,lsegid,lresid,lres,atype(i),atc_x,cg(i),amass(i), &
                      imove(i),ech(i),eha(i)
              endif
#else /**/
              read(u,fmt02a,end=6) &
                   ii,lsegid,lresid,lres,atype(i),atc_x,cg(i),amass(i), &
                   imove(i),ech_dummy,eha_dummy
#endif 
           else !cheqon3
              read(u,fmt02,end=6) &
                   ii,lsegid,lresid,lres,atype(i),atc_x,cg(i),amass(i), &
                   imove(i)
           endif cheqon3
           iac(i) = -1
           iatmtype = srchws(atc,natc,atc_x)
           if(iatmtype>0) iac(i) = iatmtype
           if(iac(i) == -1)then
              write(outu,'("Atom first ",i10,i10," ",a," problem in psf")') i, iatmtype, atc_x
              CALL WRNDIE(-3,'<psf_read_formatted>', &
                   'atom type not found for atom')
           endif
        else   !xplorfmt2
           cheqon4:IF(qcheq_psf)THEN
#if KEY_CHEQ==1
              if(qdrude_psf)then
                 read(u,fmt01b,end=6) &
                      ii,lsegid,lresid,lres,atype(i),iac(i), &
                      cg(i),amass(i),imove(i), &
                      alphadp(i),tholei(i),isdrude(i)
                 if(isdrude(i)) ndrude=ndrude+1
              else
                 read(u,fmt01a,end=6) &
                      ii,lsegid,lresid,lres,atype(i),iac(i),cg(i),amass(i), &
                      imove(i),ech(i),eha(i)
              endif
#else /**/
              read(u,fmt01a,end=6) &
                   ii,lsegid,lresid,lres,atype(i),iac(i),cg(i),amass(i), &
                   imove(i),ech_dummy,eha_dummy
#endif 
           else
              read(u,fmt01,end=6) &
                   ii,lsegid,lresid,lres,atype(i),iac(i),cg(i),amass(i), &
                   imove(i)
           endif cheqon4
        endif xplorfmt2
        rsclf(i)=1.0
#if KEY_WCA==1
        wca(i)=1.0             
#endif
        !RCZ 92/06/17 - convert to upper case
        call cnvtuc(lsegid,idleng)
        call cnvtuc(lresid,idleng)
        call cnvtuc(lres,idleng)
        call cnvtuc(atype(i),idleng)
        !
        if (lresid /= oresid.or.lres.ne.ores.or.lsegid.ne.osegid) then
           oresid=lresid
           ores=lres
           nres=nres+1
           resid(nres)=lresid
           res(nres)=lres
           ibase(nres)=i-1
        endif
        if (lsegid /= osegid) then
           osegid=lsegid
           nseg=nseg+1
           segid(nseg)=lsegid
           nictot(nseg)=nres-1
        endif
     enddo atomsloop

     nictot(nseg+1)=nres
     ibase(nres+1)=natom
     !
     !       read the rest of the psf file
     if(qappe) then
        !
        read(u,fmt00,end=6) nbondx
        if(nbond+nbondx > maxb) then
           CALL WRNDIE(-3,'<PSFRD2)', &
                'Maximum number of bonds exceeded')
           return
        endif

        read(u,fmt03,end=6) (ib(i+nbond),jb(i+nbond),i=1,nbondx)
        do i=1,nbondx
           ib(i+nbond)=ib(i+nbond)+nat
           jb(i+nbond)=jb(i+nbond)+nat
        enddo
        nbond=nbond+nbondx
        nbdrude = nbond
        !
        read(u,fmt00,end=6) nthetx
        if(ntheta+nthetx > maxt) then
           CALL WRNDIE(-3,'<PSFRD2)', &
                'Maximum number of angles exceeded')
           return
        endif
        read(u,fmt04,end=6) (it(i+ntheta),jt(i+ntheta), &
             kt(i+ntheta),i=1,nthetx)
        do i=1,nthetx
           it(i+ntheta)=it(i+ntheta)+nat
           jt(i+ntheta)=jt(i+ntheta)+nat
           kt(i+ntheta)=kt(i+ntheta)+nat
        enddo
        ntheta=ntheta+nthetx
        !
        read(u,fmt00,end=6) nphix
        if(nphi+nphix > maxp) then
           CALL WRNDIE(-3,'<PSFRD2)', &
                'Maximum number of dihedrals exceeded')
           return
        endif
        read(u,fmt03,end=6) (ip(i+nphi),jp(i+nphi), &
             kp(i+nphi),lp(i+nphi),i=1,nphix)
        do i=1,nphix
           ip(i+nphi)=ip(i+nphi)+nat
           jp(i+nphi)=jp(i+nphi)+nat
           kp(i+nphi)=kp(i+nphi)+nat
           lp(i+nphi)=lp(i+nphi)+nat
        enddo
        nphi=nphi+nphix
        !
        read(u,fmt00,end=6) nimphx
        if(nimphi+nimphx > maximp) then
           CALL WRNDIE(-3,'<PSFRD2)', &
                'Maximum number of improper dihedrals exceeded')
           return
        endif
#if KEY_CFF==1
        if(ffield == cff)then
           read(u,fmt03,end=6)(jm(i+nimphi),im(i+nimphi), &
                km(i+nimphi),lm(i+nimphi),i=1,nimphx)
        else
#endif 
           read(u,fmt03,end=6)(im(i+nimphi),jm(i+nimphi), &
                km(i+nimphi),lm(i+nimphi),i=1,nimphx)
#if KEY_CFF==1
        endif    
#endif
        do i=1,nimphx
           im(i+nimphi)=im(i+nimphi)+nat
           jm(i+nimphi)=jm(i+nimphi)+nat
           km(i+nimphi)=km(i+nimphi)+nat
           lm(i+nimphi)=lm(i+nimphi)+nat
        enddo
        nimphi=nimphi+nimphx
        !
        read(u,fmt00,end=6) ndonx
        if(ndon+ndonx > maxpad) then
           CALL WRNDIE(-3,'<PSFRD2)', &
                'Maximum number of donors exceeded')
           return
        endif
        read(u,fmt03,end=6) (idon(i+ndon),ihd1(i+ndon),i=1,ndonx)
        do i=1,ndonx
           idon(i+ndon)=idon(i+ndon)+nat
           ihd1(i+ndon)=ihd1(i+ndon)+nat
        enddo
        ndon=ndon+ndonx
        !
        read(u,fmt00,end=6) naccx
        if(nacc+naccx > maxpad) then
           CALL WRNDIE(-3,'<PSFRD2)', &
                'Maximum number of acceptors exceeded')
           return
        endif
        read(u,fmt03,end=6) (iacc(i+nacc),iac1(i+nacc),i=1,naccx)
        do i=1,naccx
           iacc(i+nacc)=iacc(i+nacc)+nat
           iac1(i+nacc)=iac1(i+nacc)+nat
        enddo
        nacc=nacc+naccx
        !
        read(u,fmt00,end=6) nnbx
        if(nnb+nnbx > maxnb) then
           CALL WRNDIE(-3,'<PSFRD2)', &
                'Maximum number of nonbond exclusions exceeded')
           return
        endif
        read(u,fmt03,end=6) (inb(i+nnb),i=1,nnbx)
        do i=1,nnbx
           inb(i+nnb)=inb(i+nnb)+nat
        enddo
        read(u,fmt03,end=6) (iblo(i+nat),i=1,natomx)
        do i=1,natomx
           iblo(i+nat)=iblo(i+nat)+nnb
        enddo
        nnb=nnb+nnbx
        !
        read(u,fmt05,end=6) ngrpx,nst2x
        !
#if KEY_NOST2==1
        IF(NST2 > 0) &
             CALL WRNDIE(-2,'<PSFIO>','ST2 code is not compiled.')
#endif 
        !
        read(u,fmt04,end=6) (igpbs(i+ngrp),igptyp(i+ngrp), &
             imoveg(i+ngrp),i=1,ngrpx)
        do i=1,ngrpx
           igpbs(i+ngrp)=igpbs(i+ngrp)+nat
        enddo
        ngrp=ngrp+ngrpx
        nst2=nst2+nst2x
        !
#if KEY_CHEQ==1
        ! zz12-05-95
        if(qcheq_psf)then
           read(u,fmt00,end=6) molntx
           read(u,fmt03,end=6) (molt(i),i=1+nat,natom)
           do i=1+nat,natom
              molt(i)=molt(i)+molnt
           enddo
           molnt=molnt+molntx
        endif
#else /**/
        if(qcheq_psf)then
           read(u,fmt00,end=6) molntx_dummy
           read(u,fmt03,end=6) (molt_dummy,i=1+nat,natom)
        endif
#endif 

#if KEY_LONEPAIR==1
        ! Read lone pair stuff
        read(u,fmt05,end=35) numlpx,numlphx
        if(numlpx > 0) then
           if (.not. allocated(lpnhost)) call allocate_lonepair ! allocate memory before reading, HJ
           do i=numlp+1,numlp+numlpx
              read(u,fmt06,end=6) lpnhost(i),lphptr(i), &
                   lpwght(i),lpvalue(1,i),lpvalue(2,i),lpvalue(3,i)
              lphptr(i)=lphptr(i)+numlph
           enddo
           read(u,fmt03,end=6) (lphost(i+numlph),i=1,numlphx)
           do i=1,numlphx
              lphost(i+numlph)=lphost(i+numlph)+nat
           enddo
           numlp=numlp+numlpx
           numlph=numlph+numlphx
        endif
35      continue
#else /**/
        read(u,fmt05,end=35) numlpx,numlphx
        if(numlpx > 0) then
           call wrndie(-1,'<psfio>', &
                'lonepairs present, but lonepair code not compiled')
           do i=1,numlpx
              read(u,*,end=6)
           enddo
           read(u,*,end=6)
        endif
35      continue
#endif 

#if KEY_CMAP==1
        if(qcmap_psf)then
           read(u,fmt00,end=36,err=75) ncrtx
           if(ncrterm+ncrtx > maxcrt) then
              call wrndie(-3,'<psfrd2)', &
                   'maximum number of cross-term maps exceeded')
              return
           endif
           read(u,fmt03,end=6)(i1ct(i+ncrterm),j1ct(i+ncrterm), &
                k1ct(i+ncrterm),l1ct(i+ncrterm), &
                i2ct(i+ncrterm),j2ct(i+ncrterm), &
                k2ct(i+ncrterm),l2ct(i+ncrterm), &
                i=1,ncrtx)
           do i=1,ncrtx
              i1ct(i+ncrterm)=i1ct(i+ncrterm)+nat
              j1ct(i+ncrterm)=j1ct(i+ncrterm)+nat
              k1ct(i+ncrterm)=k1ct(i+ncrterm)+nat
              l1ct(i+ncrterm)=l1ct(i+ncrterm)+nat
              i2ct(i+ncrterm)=i2ct(i+ncrterm)+nat
              j2ct(i+ncrterm)=j2ct(i+ncrterm)+nat
              k2ct(i+ncrterm)=k2ct(i+ncrterm)+nat
              l2ct(i+ncrterm)=l2ct(i+ncrterm)+nat
           enddo
           ncrterm=ncrterm+ncrtx
        endif
36      continue
#endif 
        !
     else
        read(u,fmt00,end=6) nbond
        read(u,fmt03,end=6) (ib(i),jb(i),i=1,nbond)
        nbdrude = nbond
        read(u,fmt00,end=6) ntheta
        read(u,fmt04,end=6) (it(i),jt(i),kt(i),i=1,ntheta)
        read(u,fmt00,end=6) nphi
        read(u,fmt03,end=6) (ip(i),jp(i),kp(i),lp(i),i=1,nphi)
        read(u,fmt00,end=6) nimphi
#if KEY_CFF==1
        if(ffield == cff)then
           read(u,fmt03,end=6)(jm(i),im(i),km(i),lm(i),i=1,nimphi)
        else
#endif 
           read(u,fmt03,end=6)(im(i),jm(i),km(i),lm(i),i=1,nimphi)
#if KEY_CFF==1
        ENDIF    
#endif
        read(u,fmt00,end=6) ndon
        read(u,fmt03,end=6) (idon(i),ihd1(i),i=1,ndon)
        read(u,fmt00,end=6) nacc
        read(u,fmt03,end=6) (iacc(i),iac1(i),i=1,nacc)
        read(u,fmt00,end=6) nnb
        read(u,fmt03,end=6) (inb(i),i=1,nnb)
        read(u,fmt03,end=6) (iblo(i),i=1,natom)
        read(u,fmt05,end=6) ngrp,nst2
#if KEY_NOST2==1
        IF(NST2 > 0) &
             CALL WRNDIE(-2,'<PSFIO>','ST2 code is not compiled.')
#endif 
        read(u,fmt04,end=6) (igpbs(i),igptyp(i),imoveg(i),i=1,ngrp)
        !
        if(qcheq_psf)then
#if KEY_CHEQ==1
           read(u,fmt00,end=6) molnt
           read(u,fmt03,end=6) (molt(i),i=1,natom)
#else /**/
           read(u,fmt00,end=6) molnt_dummy
           read(u,fmt03,end=6) (molt_dummy,i=1,natom)
#endif 
        endif
#if KEY_LONEPAIR==1
        ! Read lone pair stuff
        numlp=0
        numlph=0
        read(u,fmt05,end=45) numlpx,numlphx
        if(numlpx > 0) then
           if (.not. allocated(lpnhost)) call allocate_lonepair ! allocate memory before reading, HJ
           do i=1,numlpx
              read(u,fmt06,end=6) lpnhost(i),lphptr(i), &
                   lpwght(i),lpvalue(1,i),lpvalue(2,i),lpvalue(3,i)
           enddo
           read(u,fmt03,end=6) (lphost(i),i=1,numlphx)
        endif
        numlp=numlpx
        numlph=numlphx
45      continue
        ! read anisotropic terms, if any
        if(qdrude) then
           naniso=0
           read(u,fmt00,end=37) nanisox
           if(nanisox > 0) then
              do i=1,nanisox
                 read(u,fmt07,end=6) k11(i),k22(i),k33(i)
              enddo
              read(u,fmt03,end=6) (lstani1(i),lstani2(i),lstani3(i), &
                   lstani4(i),i=1,nanisox)
           endif
37         continue
           naniso=nanisox
        endif
#else /**/
        read(u,fmt05,end=45) numlpx,numlphx
        if(numlpx > 0) then
           CALL WRNDIE(-1,'<PSFIO>', &
                'Lonepairs present, but LONEPAIR code not compiled')
           do i=1,numlpx
              read(u,*,end=6)
           enddo
           read(u,*,end=6)
        endif
45      continue
#endif 

#if KEY_CMAP==1
        if(qcmap_psf)then
           read(u,fmt00,end=46,err=75) ncrtx
           read(u,fmt03,end=6)(i1ct(i),j1ct(i),k1ct(i),l1ct(i), &
                i2ct(i),j2ct(i),k2ct(i),l2ct(i), &
                i=1,ncrtx)
           ncrterm=ncrtx
        endif
46      continue
#endif 
        !
     endif
     !
     igpbs(ngrp+1)=natom
     !
     !
     !-eof-label----
     goto 77
#if KEY_CMAP==1
75   write(OUTU,'(//2a)')' %PSFRD2-ERR: CMAP in psf header ', &
          'But no CMAP info in psf'
     CALL WRNDIE(-4,'<PSFRD2>','No CMAP entries in psf')
     goto 77
#endif 
6    eof=.true.
     call wrndie(-1,'<psfrd2>','eof during read')
77   continue
     !-----------------
     !
     !       make dimension check
     if(wrnlev >= 2) then
        if (natom > maxa) then
           WRITE(OUTU,'(A)') ' %PSFRD2-ERR: MAXA exceeded'
        ELSEIF (NBOND > MAXB) THEN
           WRITE(OUTU,'(A)') ' %PSFRD2-ERR: MAXB exceeded'
        ELSEIF (NTHETA > MAXT) THEN
           WRITE(OUTU,'(A)') ' %PSFRD2-ERR: MAXT exceeded'
        ELSEIF (NPHI > MAXP) THEN
           WRITE(OUTU,'(A)') ' %PSFRD2-ERR: NPHI exceeded'
        ELSEIF (NIMPHI > MAXIMP) THEN
           WRITE(OUTU,'(A)')' %PSFRD2-ERR: NIMPHI exceeded'
#if KEY_CMAP==1
        ELSEIF (NCRTERM > MAXCRT) THEN
           WRITE(OUTU,'(A)')' %PSFRD2-ERR: NCRTERM exceeded'
#endif 
        ELSEIF (NDON > MAXPAD) THEN
           WRITE(OUTU,'(A)') ' %PSFRD2-ERR: MAXPAD exceeded'
        ELSEIF (NACC > MAXPAD) THEN
           WRITE(OUTU,'(A)') ' %PSFRD2-ERR: MAXPAD exceeded'
        ELSEIF (NNB > MAXNB) THEN
           WRITE(OUTU,'(A)') ' %PSFRD2-ERR: MAXNB exceeded'
        ELSEIF (NGRP > MAXGRP) THEN
           WRITE(OUTU,'(A)') ' %PSFRD2-ERR: MAXGRP exceeded'
        ELSEIF (NRES > MAXRES) THEN
           WRITE(OUTU,'(A)') ' %PSFRD2-ERR: MAXRES exceeded'
        ELSEIF (NSEG > MAXSEG) THEN
           WRITE(OUTU,'(A)') ' %PSFRD2-ERR: MAXSEG exceeded'
        ENDIF
     ENDIF
  ENDIF masterio

#if KEY_PARALLEL==1 /*par_ens*/
  call psnd4(qdrude,1)
  if(qdrude) then 
     if( .not. allocated(lstani1))  call allocate_aniso
#if KEY_LONEPAIR==1
     if (.not. allocated(lpnhost)) call allocate_lonepair  
#endif
  endif
  call psnd4(natom,1)
  call psndc(atype,natom)
  call psnd4(iac,natom)
  call psnd8(cg,natom)
  call psnd8(amass,natom)
  call psnd8(rsclf,natom)
#if KEY_WCA==1
  call psnd8(wca,natom)               
#endif
  call psnd4(imove,natom)
  !     drude infromation
  call psnd4(isdrude,natom)
  call psnd4(ndrude,1)
  call psnd4(qdrude,1)
  if(qdrude)then
     call psnd8(alphadp,natom)
     call psnd8(tholei,natom)
  endif
  !
  call psnd4(nseg,1)
  call psnd4(nictot,nseg+1)
  call psndc(segid,nseg)
  !
  call psnd4(nres,1)
  call psnd4(ibase,nres+1)
  call psndc(resid,nres)
  call psndc(res,nres)
  !
  call psnd4(ngrp,1)
  call psnd4(igpbs,ngrp+1)
  call psnd4(igptyp,ngrp)
  call psnd4(imoveg,ngrp)
  !
  call psnd4(nbond,1)
  call psnd4(ib,nbond)
  call psnd4(jb,nbond)
  if(qdrude)then
     call psnd4(nbdrude,1)
  endif
  !
  call psnd4(ntheta,1)
  call psnd4(it,ntheta)
  call psnd4(jt,ntheta)
  call psnd4(kt,ntheta)
  !
  call psnd4(nphi,1)
  call psnd4(ip,nphi)
  call psnd4(jp,nphi)
  call psnd4(kp,nphi)
  call psnd4(lp,nphi)
  !
  call psnd4(nimphi,1)
  call psnd4(im,nimphi)
  call psnd4(jm,nimphi)
  call psnd4(km,nimphi)
  call psnd4(lm,nimphi)
#if KEY_CMAP==1 /*cmap*/
  call psnd4(ncrterm,1)
  call psnd4(i1ct,ncrterm)
  call psnd4(j1ct,ncrterm)
  call psnd4(k1ct,ncrterm)
  call psnd4(l1ct,ncrterm)
  call psnd4(i2ct,ncrterm)
  call psnd4(j2ct,ncrterm)
  call psnd4(k2ct,ncrterm)
  call psnd4(l2ct,ncrterm)
#endif /* (cmap)*/
  call psnd4(ndon,1)
  call psnd4(idon,ndon)
  call psnd4(ihd1,ndon)
  call psnd4(nacc,1)
  call psnd4(iacc,nacc)
  call psnd4(iac1,nacc)

  call psnd4(nnb,1)
  call psnd4(inb,nnb)
  call psnd4(iblo,natom)
  
  call psnd4(nst2,1)
  
#if KEY_LONEPAIR==1 /*send_lonepairs*/
     call psnd4(numlp,1)
     call psnd4(numlph,1)
     if(.not.allocated(lpnhost))call allocate_lonepair
     call psnd4(lpnhost,numlp)
     call psnd4(lphptr,numlp)
     call psnd4(lpwght,numlp)
     call psnd8(lpvalue,numlp*3)
     call psnd4(lphost,numlph)
     if(qdrude)then
        if(.not.allocated(lstani1))call allocate_aniso
        call psnd4(naniso,1)
        call psnd8(k11,naniso)
        call psnd8(k22,naniso)
        call psnd8(k33,naniso)
        call psnd4(lstani1,naniso)
        call psnd4(lstani2,naniso)
        call psnd4(lstani3,naniso)
        call psnd4(lstani4,naniso)
     endif
#endif /* (send_lonepairs)*/
#endif /* (par_ens)*/

#if KEY_MMFF==1
  if(ffield == mmff) then
     do i=1,nbond
        bondtype(i)=1
     enddo
     bondtype(3)=2
     bondtype(7)=2
  endif
#endif 
  return
end subroutine psf_read_formatted


        
SUBROUTINE PSFWR2(U,TITLE,NTITLE)
  !-----------------------------------------------------------------------
  !     This routine writes the protein structure file module
  !     Axel Brunger, 13-AUG-84
  !     =======================
  !
#if KEY_CHEQ==1
  use cheq,only:ech,eha,molt,molnt      
#endif
  use chm_kinds
  use dimens_fcm
  use exfunc
  !     input/output
  !RCZ 10/24/91 to write XPLOR compatible PSF file
  use comand
  use param
  !RCZ
  use psf
  use lonepr
  use aniso_fcm
  use stream
  use string
  !
#if KEY_STRINGM==1 /*   VO string method */
    use machio, only: ifreeu
    use multicom_aux
#endif
  !
  implicit none
  INTEGER U
  CHARACTER(len=*) TITLE(*)
  INTEGER NTITLE
  !     local
  INTEGER I, IRES, ISEG
  !RCZ 91/04/26
  LOGICAL  QXPLOR
  !yw  24-Jan-2003 expansion code
  integer wtlen
  character(len=10) :: fmt00,fmt03,fmt04,fmt05,fmt08
  character(len=30) :: fmt06,fmt07,heading
  character(len=60) :: fmt01,fmt01a,fmt01b,fmt02,fmt02a,fmt02b
  !
#if KEY_STRINGM==1 /*  VO stringm v */
  logical :: qstr
  common /replicaio/ qstr ! need global variable
  qstr=u.le.size(ifreeu) ! protect against OOB on ifreeu
  if (qstr) qstr=((MOD(IFREEU(U),8).eq.0).and.(IFREEU(U).ne.0).and.(ME_LOCAL.eq.0))
#endif

  IF(IOLEV<0 &
#if KEY_STRINGM==1 /*  VO */
    &         .and..not.qstr &  
#endif
    &                         ) RETURN

  fmt00 (1:10)=' '
  fmt01 (1:60)=' '
  fmt01a(1:60)=' '
  fmt02 (1:60)=' '
  fmt02a(1:60)=' '
  fmt02b(1:60)=' '
  fmt03 (1:10)=' '
  fmt04 (1:10)=' '
  fmt05 (1:10)=' '
  fmt06 (1:20)=' '
  fmt07 (1:20)=' '
  fmt08 (1:10)=' '
  heading(1:20)=' '
  qextfmt=qxform()
  if (qextfmt) then
     !23456789     +123456789+123456789+123456789+123456789+123456789+123456789+
     heading='PSF EXT'
     fmt00='(/I10,A)'
     fmt01= '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)'
     fmt01a='(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6)'
     fmt01b='(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6,L1)'
     fmt02= '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6)'
     fmt02a='(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6)'
     fmt02b='(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6)'
     fmt03='(8I10)'
     fmt04='(9I10)'
     fmt05='(/2I10,A)'
     fmt06='(2I10,3X,L1,3G14.6)'
     fmt07='(10X,3G14.6)'
     fmt08='(/I10,A)'
  else
     heading='PSF'
     fmt00='(/I8,A)'
     fmt01='(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)'
     fmt01a='(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6)'
     fmt02='(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)'
     fmt02a='(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6)'
     fmt03='(8I8)'
     fmt04='(9I8)'
     fmt05='(/2I8,A)'
     fmt06='(2I8,3X,L1,3G14.6)'
     fmt07='(8X,3G14.6)'
     fmt08='(/I8,A)'
  endif
  QXPLOR=INDXA(COMLYN,COMLEN,'XPLO') > 0
  ! Always write an xplor psf unless specifically asked not to
  qxplor = qxplor .or. .not. (indxa(comlyn,comlen,'OLDPSF')>0)
  !RCZ 91/04/26
  !
  wtlen=len(heading)
  call trime(heading,wtlen)
#if KEY_CMAP==1
  heading(wtlen+1:wtlen+5)=' CMAP'
  wtlen=len(heading)
  call trime(heading,wtlen)
#endif 
#if KEY_CHEQ==1
  heading(wtlen+1:wtlen+5)=' CHEQ'
  wtlen=len(heading)
  call trime(heading,wtlen)
#endif 
  IF(NDRUDE > 0) THEN
     heading(wtlen+1:wtlen+6)=' DRUDE'
     wtlen=len(heading)
     call trime(heading,wtlen)
  ENDIF
  if(qxplor) then
     heading(wtlen+1:wtlen+6)=' XPLOR'
     wtlen=len(heading)
     call trime(heading,wtlen)
  endif
  WRITE(U,'(A)') heading(1:wtlen)
  CALL WRTITL(TITLE,NTITLE,0,+2)
  WRITE(U,fmt00) NTITLE,' !NTITLE'
  WRITE(U,'(A)') (TITLE(I),I=1,NTITLE)
  WRITE(U,fmt00) NATOM, ' !NATOM'
  DO ISEG=1,NSEG
     DO IRES=NICTOT(ISEG)+1,NICTOT(ISEG+1)
        DO I=IBASE(IRES)+1,IBASE(IRES+1)
           !RCZ 91/04/26 use QXPLOR for Xplor format
#if KEY_CHEQ==1
           IF (QXPLOR) THEN
              IF(NDRUDE > 0) THEN
                 WRITE(U,fmt02b) &
                      I,SEGID(ISEG),RESID(IRES),RES(IRES),ATYPE(I),ATC(IAC(I)), &
                      CG(I),AMASS(I),IMOVE(I),ALPHADP(I),THOLEI(I)
              ELSE
                 WRITE(U,fmt02a) &
                      I,SEGID(ISEG),RESID(IRES),RES(IRES),ATYPE(I),ATC(IAC(I)), &
                      CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)
              ENDIF
           ELSE 
              IF(NDRUDE > 0) THEN
                 WRITE(U,fmt01b) &
                      I,SEGID(ISEG),RESID(IRES),RES(IRES),ATYPE(I),IAC(I), &
                      CG(I),AMASS(I),IMOVE(I),ALPHADP(I),THOLEI(I),ISDRUDE(I)
              ELSE              
                 WRITE(U,fmt01a) &
                      I,SEGID(ISEG),RESID(IRES),RES(IRES),ATYPE(I),IAC(I), &
                      CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)
              ENDIF
           ENDIF
#else /**/
           IF (QXPLOR) THEN
              WRITE(U,fmt02) &
                   I,SEGID(ISEG),RESID(IRES),RES(IRES),ATYPE(I),ATC(IAC(I)), &
                   CG(I),AMASS(I),IMOVE(I)
           ELSE 
              WRITE(U,fmt01) &
                   I,SEGID(ISEG),RESID(IRES),RES(IRES),ATYPE(I),IAC(I), &
                   CG(I),AMASS(I),IMOVE(I)
           ENDIF
#endif 
        ENDDO
     ENDDO
  ENDDO
  WRITE(U,fmt00) NBOND, ' !NBOND: bonds'
  WRITE(U,fmt03) (IB(I),JB(I),I=1,NBOND)
  WRITE(U,fmt00) NTHETA,' !NTHETA: angles'
  WRITE(U,fmt04) (IT(I),JT(I),KT(I),I=1,NTHETA)
  WRITE(U,fmt00) NPHI,' !NPHI: dihedrals'
  WRITE(U,fmt03) (IP(I),JP(I),KP(I),LP(I),I=1,NPHI)
  WRITE(U,fmt00) NIMPHI,' !NIMPHI: impropers'
  WRITE(U,fmt03)(IM(I),JM(I),KM(I),LM(I),I=1,NIMPHI)
  WRITE(U,fmt00) NDON,' !NDON: donors'
  WRITE(U,fmt03) (IDON(I),IHD1(I),I=1,NDON)
  WRITE(U,fmt00) NACC,' !NACC: acceptors'
  WRITE(U,fmt03) (IACC(I),IAC1(I),I=1,NACC)
  WRITE(U,fmt00) NNB,' !NNB'
  WRITE(U,fmt03) (INB(I),I=1,NNB)
  WRITE(U,fmt03) (IBLO(I),I=1,NATOM)
  WRITE(U,fmt05) NGRP,NST2,' !NGRP NST2'
  WRITE(U,fmt04) (IGPBS(I),IGPTYP(I),IMOVEG(I),I=1,NGRP)

#if KEY_CHEQ==1
  WRITE(U,fmt00) MOLNT,' !MOLNT'
  WRITE(U,fmt03) (MOLT(I),I=1,NATOM)
#endif 

#if KEY_LONEPAIR==1
  ! integers
  WRITE(U,fmt05) NUMLP,NUMLPH,' !NUMLP NUMLPH'
  IF(NUMLP > 0) THEN
     DO I=1,NUMLP
        WRITE(U,fmt06) LPNHOST(I),LPHPTR(I), &
             LPWGHT(I),LPVALUE(1,I),LPVALUE(2,I),LPVALUE(3,I)
     ENDDO
     WRITE(U,fmt03) (LPHOST(I),I=1,NUMLPH)
  ENDIF

  !c36a3      IF (QXPLOR) THEN
  IF(NDRUDE > 0) THEN
     WRITE(U,fmt00) NANISO,' !NUMANISO'
     IF(NANISO > 0) THEN
        DO I=1,NANISO
           WRITE(U,fmt07) K11(I),K22(I),K33(I)
        ENDDO
        WRITE(U,fmt03) (LSTANI1(I),LSTANI2(I),LSTANI3(I),LSTANI4(I), &
             I=1,NANISO)
     ENDIF
  ENDIF
  !c36a3      ENDIF
#else /**/
  WRITE(U,fmt05) 0,0,' !NUMLP NUMLPH'
#endif 

#if KEY_CMAP==1
  WRITE(U,fmt00) NCRTERM,' !NCRTERM: cross-terms'
  WRITE(U,fmt03)(I1CT(I),J1CT(I),K1CT(I),L1CT(I), &
       I2CT(I),J2CT(I),K2CT(I),L2CT(I), &
       I=1,NCRTERM)
#endif 
  !
  RETURN
END SUBROUTINE PSFWR2
        

!----------------------- Not integer*4 binaries ---------------------------
! ------ use native integers in binaries, binary files not portable -MFC===
!
SUBROUTINE PSFRDR(IUNIT,ICNTRL,TITLE,NTITL,QAPPE)
  !
  !  THIS ROUTINE READS A BINARY FILE PSF FROM THE DISK
  !
  !     binary format
  !     ===============
  !
  !     Authors: Robert Bruccoleri
  !     Bernie Brooks
  !
#if KEY_CHEQ==1
  use cheq,only:eha,ech,molnt,molt           
#endif
  use chm_kinds
  use dimens_fcm
  use psf
  use lonepr
  use stream
  use string
  use comand
  use param
  use memory
  !
#if KEY_STRINGM==1 /*   VO string method */
  use machio, only: ifreeu
  use multicom_aux
#endif
  !
  use exfunc, only : srchws
  implicit none
  INTEGER ICNTRL(20)
  CHARACTER(len=*) TITLE(*)
  INTEGER NTITL,IUNIT
  LOGICAL QAPPE
  !
  CHARACTER(len=4) :: PSFHDR='PSF ',HDR
  INTEGER NSEG1,IGRP,I,J,K
  INTEGER NSEGX,NRESX,NATOMX,NAT,NBONDX,NTHETX,NPHIX, &
       NIMPHX,NNBX,NDONX,NACCX,NGRPX,NST2X,NUMLPX,NUMLPHX
  INTEGER IJUNK,IS,IQ,oleng
  real(chm_real)  CGT
  logical :: qxplor
  character(len=8), allocatable, dimension(:) :: atmcode
#if KEY_CHEQ==1
  integer molntx
  LOGICAL qcheq_psf
#endif 
  !
  !
#if KEY_STRINGM==1 /*  VO stringm v */
  logical :: qstr
  common /replicaio/ qstr ! need global variable
  qstr=iunit.le.size(ifreeu) ! protect against OOB on ifreeu
  if (qstr) qstr=((MOD(IFREEU(IUNIT),8).eq.0).and.(IFREEU(IUNIT).ne.0).and.(ME_LOCAL.eq.0))
#endif
  !
  IF(IOLEV > 0 &
#if KEY_STRINGM==1 /*  VO */
    &         .or.qstr &  
#endif
    &                         ) THEN
     CALL TRYORO(IUNIT,'UNFORMATTED')
     IF (reallow) THEN      
        REWIND (UNIT=IUNIT)
     ENDIF                  
     READ(IUNIT) HDR,ICNTRL
     IF(HDR /= PSFHDR) THEN
        IF(WRNLEV >= 2) WRITE(OUTU,201) PSFHDR,HDR
201     FORMAT(' EXPECTED = "',A4,'" FOUND = "',A4,'"')
        CALL WRNDIE(-1,'<PSFIO>','HEADERS DONT MATCH')
     ENDIF
     !
     !       CHECK THE VERSION OF THE PSF MODULE. OLD VERSIONS (ICNTRL(1)=1) DID
     !       NOT PUT OUT A RECORD FOR THE IMPROPERS IF THERE WERE NO IMPROPERS,
     !       BUT DID PUT OUT A RECORD FOR ALL OTHER ZERO LENGTH ARRAYS. ALSO
     !       FORTRAN 77 DO LOOPS CHANGE THESE FROM ONE TO ZERO LENGTH RECORDS.
     !
     IF (ICNTRL(1) < 2 .AND. WRNLEV >= 2) WRITE(OUTU,10) ICNTRL(1)
10   FORMAT(' *****  WARNING  *****  ERRORS WILL OCCUR READING', &
          ' OLD PSFS WITH NO IMPROPER TORSIONS IN THE STRUCTURE'/ &
          ' ICNTRL(1)=',I2,' THIS VERSION OF PSFIO EXPECTS ICNTRL(1)>=2')
#if KEY_CHEQ==1
     qcheq_psf=(ICNTRL(13) == 1) 
#else /**/
     if (ICNTRL(13) == 1) then
        IF(WRNLEV >= 2) WRITE(OUTU,202) 
202     FORMAT(' PSF written using CHEQ (c30+ version) with CHEQ', &
             /'    CHEQ not compiled for this charmm')
        CALL WRNDIE(-1,'<PSFIO>', &
             'Cannot read binary psf containg CHEQ')
     endif
#endif 
     READ(IUNIT) NTITL,NSEGX
     CALL RDTITL(TITLE,NTITL,IUNIT,-2)
     CALL WRTITL(TITLE,NTITL,OUTU,1)
     !
     oleng=idleng
     if (icntrl(1) >= 30) then
        idleng=8
        if (prnlev >= 2) write(outu,'(a)') &
             ' PSFRDR> Reading extended format binary PSF file.'
     else
        idleng=4
     endif

     qxplor = (icntrl(20)==1)

     IF(QAPPE) THEN
        IF(NSEG+NSEGX > MAXSEG) THEN
           CALL WRNDIE(-3,'<PSFRDR)', &
                'Maximum number of segments exceeded')
           RETURN
        ENDIF
        NSEG1=NSEGX+1
        READ(IUNIT) (NICTOT(J+NSEG),NATOMX,NBONDX,NTHETX,NPHIX, &
             NIMPHX,NNBX,NDONX,NACCX,K,J=1,NSEG1)
        NRESX=NICTOT(NSEG+NSEG1)
        DO J=1,NSEG1
           NICTOT(J+NSEG)=NICTOT(J+NSEG)+NRES
        ENDDO
        READ(IUNIT) (SEGID(I+NSEG)(1:idleng),I=1,NSEGX)
        NSEG=NSEGX+NSEG
        !
        IF(NRES+NRESX > MAXRES) THEN
           CALL WRNDIE(-3,'<PSFRDR)', &
                'Maximum number of residues exceeded')
           RETURN
        ENDIF
        READ(IUNIT) (RES(I+NRES)(1:idleng),RESID(I+NRES)(1:idleng), &
             IJUNK,I=1,NRESX)
        READ(IUNIT) (IBASE(I+NRES),I=1,NRESX+1)
        DO I=1,NRESX+1 
           IBASE(I+NRES)=IBASE(I+NRES)+NATOM
        ENDDO
        NRES=NICTOT(NSEG+1)
        !
        NAT=NATOM
        IF(NATOM+NATOMX > MAXA) THEN
           CALL WRNDIE(-3,'<PSFRDR)', &
                'Maximum number of atoms exceeded')
           RETURN
        ENDIF
        !
        if(qxplor) &
             call chmalloc('psfres.src','psfrdr','atmcode',natom+natomx,ch8=atmcode)
        !
#if KEY_CHEQ==1
        if(qcheq_psf)then
           if(qxplor) then
              READ(IUNIT) (ATYPE(I+NAT)(1:idleng),CG(I+NAT),AMASS(I+NAT), &
                   atmcode(I+NAT),IMOVE(I+NAT),ECH(I+NAT),EHA(I+NAT), &
                   I=1,NATOMX)
           else
              READ(IUNIT) (ATYPE(I+NAT)(1:idleng),CG(I+NAT),AMASS(I+NAT), &
                   IAC(I+NAT),IMOVE(I+NAT),ECH(I+NAT),EHA(I+NAT), &
                   I=1,NATOMX)
           endif
        else
#endif 
           if(qxplor) then
              READ(IUNIT) (ATYPE(I+NAT)(1:idleng),CG(I+NAT),AMASS(I+NAT), &
                   atmcode(I+nat),IMOVE(I+nat),I=1,NATOMX)
           else
              READ(IUNIT) (ATYPE(I+NAT)(1:idleng),CG(I+NAT),AMASS(I+NAT), &
                   IAC(I+nat),IMOVE(I+nat),I=1,NATOMX)
           endif
#if KEY_CHEQ==1
        endif                                             
#endif
        DO I=1,NATOMX
           RSCLF(I+NAT)=1.0
#if KEY_WCA==1
           WCA(I+NAT)=1.0                                  
#endif
           if(qxplor) iac(i+nat) = srchws(atc,natc,atmcode(i+nat))
        ENDDO
        NATOM=NATOM+NATOMX
        !
        IF(NBOND+NBONDX > MAXB) THEN
           CALL WRNDIE(-3,'<PSFRDR)', &
                'Maximum number of bonds exceeded')
           RETURN
        ENDIF
        READ(IUNIT) (IB(I+NBOND),JB(I+NBOND),I=1,NBONDX)
        DO I=1,NBONDX
           IB(I+NBOND)=IB(I+NBOND)+NAT
           JB(I+NBOND)=JB(I+NBOND)+NAT
        ENDDO
        NBOND=NBOND+NBONDX
        !
        IF(NTHETA+NTHETX > MAXT) THEN
           CALL WRNDIE(-3,'<PSFRDR)', &
                'Maximum number of angles exceeded')
           RETURN
        ENDIF
        READ(IUNIT)(IT(I+NTHETA),JT(I+NTHETA),KT(I+NTHETA),I=1,NTHETX)
        DO I=1,NTHETX
           IT(I+NTHETA)=IT(I+NTHETA)+NAT
           JT(I+NTHETA)=JT(I+NTHETA)+NAT
           KT(I+NTHETA)=KT(I+NTHETA)+NAT
        ENDDO
        NTHETA=NTHETA+NTHETX
        !
        IF(NPHI+NPHIX > MAXP) THEN
           CALL WRNDIE(-3,'<PSFRDR)', &
                'Maximum number of dihedrals exceeded')
           RETURN
        ENDIF
        READ(IUNIT) (IP(I+NPHI),JP(I+NPHI),KP(I+NPHI), &
             LP(I+NPHI),I=1,NPHIX)
        DO I=1,NPHIX
           IP(I+NPHI)=IP(I+NPHI)+NAT
           JP(I+NPHI)=JP(I+NPHI)+NAT
           KP(I+NPHI)=KP(I+NPHI)+NAT
           LP(I+NPHI)=LP(I+NPHI)+NAT
        ENDDO
        NPHI=NPHI+NPHIX
        !
        IF(NIMPHI+NIMPHX > MAXIMP) THEN
           CALL WRNDIE(-3,'<PSFRDR)', &
                'Maximum number of improper dihedrals exceeded')
           RETURN
        ENDIF
        READ(IUNIT) (IM(I+NIMPHI),JM(I+NIMPHI),KM(I+NIMPHI), &
             LM(I+NIMPHI),I=1,NIMPHX)
        DO I=1,NIMPHX
           IM(I+NIMPHI)=IM(I+NIMPHI)+NAT
           JM(I+NIMPHI)=JM(I+NIMPHI)+NAT
           KM(I+NIMPHI)=KM(I+NIMPHI)+NAT
           LM(I+NIMPHI)=LM(I+NIMPHI)+NAT
        ENDDO
        NIMPHI=NIMPHI+NIMPHX
        !
        IF(NDON+NDONX > MAXPAD) THEN
           CALL WRNDIE(-3,'<PSFRDR)', &
                'Maximum number of donors exceeded')
           RETURN
        ENDIF
        READ(IUNIT) (IDON(I+NDON),IJUNK,IJUNK,IHD1(I+NDON),I=1,NDONX)
        DO I=1,NDONX
           IDON(I+NDON)=IDON(I+NDON)+NAT
           IHD1(I+NDON)=IHD1(I+NDON)+NAT
        ENDDO
        NDON=NDON+NDONX
        !
        IF(NACC+NACCX > MAXPAD) THEN
           CALL WRNDIE(-3,'<PSFRDR)', &
                'Maximum number of acceptors exceeded')
           RETURN
        ENDIF
        READ(IUNIT) (IACC(I+NACC),IAC1(I+NACC),IJUNK,I=1,NACCX)
        DO I=1,NACCX
           IACC(I+NACC)=IACC(I+NACC)+NAT
           IAC1(I+NACC)=IAC1(I+NACC)+NAT
        ENDDO
        NACC=NACC+NACCX
        !
        IF(NNB+NNBX > MAXNB) THEN
           CALL WRNDIE(-3,'<PSFRDR)', &
                'Maximum number of nonbond exclusions exceeded')
           RETURN
        ENDIF
        READ(IUNIT) (INB(I+NNB),I=1,NNBX)
        DO I=1,NNBX
           INB(I+NNB)=INB(I+NNB)+NAT
        ENDDO
        READ(IUNIT) (IBLO(I+NAT),I=1,NATOMX)
        DO I=1,NATOMX
           IBLO(I+NAT)=IBLO(I+NAT)+NNB
        ENDDO
        NNB=NNB+NNBX
        !
        if(qxplor) &
             call chmdealloc('psfres.src','psfrdr','atmcode',natom+natomx,ch8=atmcode)
        !

        !
     ELSE
#if KEY_CHEQ==1
        NAT=0                                                 
#endif
#if KEY_CHEQ==1
        MOLNT=0                                               
#endif
        NSEG=NSEGX
        NSEG1=NSEG+1
        READ(IUNIT) (NICTOT(J),NATOM,NBOND,NTHETA,NPHI, &
             NIMPHI,NNB,NDON,NACC,K,J=1,NSEG1)
        NRES=NICTOT(NSEG1)
        READ(IUNIT) (SEGID(I)(1:idleng),I=1,NSEG)
        READ(IUNIT) (RES(I)(1:idleng),RESID(I)(1:idleng), &
             IJUNK,I=1,NRES)
        READ(IUNIT) (IBASE(I),I=1,NRES+1)
        !
        if(qxplor) &
             call chmalloc('psfres.src','psfrdr','atmcode',natom,ch8=atmcode)
        !
        if(qxplor) then
           READ(IUNIT) (ATYPE(I)(1:idleng),CG(I),AMASS(I),atmcode(I), &
                IMOVE(I), &
#if KEY_CHEQ==1
                ECH(I),EHA(I),                            & 
#endif
                I=1,NATOM)
        else
           READ(IUNIT) (ATYPE(I)(1:idleng),CG(I),AMASS(I),IAC(I), &
                IMOVE(I), &
#if KEY_CHEQ==1
                ECH(I),EHA(I),                            & 
#endif
                I=1,NATOM)
        endif
        DO I=1,NATOM
           RSCLF(I)=1.0
#if KEY_WCA==1
           WCA(I)=1.0                                          
#endif
           if(qxplor) iac(i) = srchws(atc,natc,atmcode(i))
        ENDDO
        !
        if(qxplor) &
             call chmdealloc('psfres.src','psfrdr','atmcode',natom,ch8=atmcode)
        !
        READ(IUNIT) (IB(I),JB(I),I=1,NBOND)
        READ(IUNIT) (IT(I),JT(I),KT(I),I=1,NTHETA)
        READ(IUNIT) (IP(I),JP(I),KP(I),LP(I),I=1,NPHI)
        READ(IUNIT) (IM(I),JM(I),KM(I),LM(I),I=1,NIMPHI)
        READ(IUNIT) (IDON(I),IJUNK,IJUNK,IHD1(I),I=1,NDON)
        READ(IUNIT) (IACC(I),IAC1(I),IJUNK,I=1,NACC)
        READ(IUNIT) (INB(I),I=1,NNB)
        READ(IUNIT) (IBLO(I),I=1,NATOM)
     ENDIF
     !
     !       GET GROUP INFORMATION
     IF (ICNTRL(1) < 17) THEN
        !         old psf, just copy from residue arrays
        IF(QAPPE) CALL WRNDIE(-2,'<PSFRDR>', &
             'Group information is lost when reading an old PSF with APPEnd')
        NGRP=NRES
        NST2=0

        IGPBS(NRES+1)=IBASE(NRES+1)
        DO IGRP=1,NGRP
           IGPBS(IGRP)=IBASE(IGRP)
           IS=IBASE(IGRP)+1
           IQ=IBASE(IGRP+1)
           IF (RES(IGRP) == 'ST2 ') THEN
              !      ----- Regenerate the group information -----------
              !      ----- count st2 
              !      ----- get group types assigned as if for new type psf
#if KEY_NOST2==0
              NST2=NST2+1
              IGPTYP(IGRP)=3
              IMOVE(IQ-1)=-1
              IMOVE(IQ)=-1
#else /**/
              CALL WRNDIE(-1,'<PSFIO>','ST2 code is not compiled.')
#endif 
           ELSE
              IGPTYP(IGRP)=0
              CGT=0.0
              DO I=IS,IQ
                 CGT=CGT+CG(I)
                 IF(CG(I) /= 0.0) IGPTYP(IGRP)=1
              ENDDO
              IF(ABS(CGT) > 0.00001) IGPTYP(IGRP)=2
           ENDIF
           !
           !           FIX UP IMOVEG ARRAY
           J=IMOVE(IS)
           DO I=IS,IQ
              IF(IMOVE(I) == 0) J=0
           ENDDO
           IMOVEG(IGRP)=J
        ENDDO
        !
     ELSE
        !         this is a new psf, just read group stuff
        IF(QAPPE) THEN
           READ(IUNIT) NGRPX,NST2X
           READ(IUNIT) (IGPBS(I+NGRP),IGPTYP(I+NGRP), &
                IMOVEG(I+NGRP),I=1,NGRPX)
           DO I=1,NGRPX
              IGPBS(I+NGRP)=IGPBS(I+NGRP)+NAT
           ENDDO
           NGRP=NGRP+NGRPX
           NST2=NST2+NST2X
        ELSE
           READ(IUNIT) NGRP,NST2
           READ(IUNIT) (IGPBS(I),IGPTYP(I),IMOVEG(I),I=1,NGRP)
        ENDIF
        IGPBS(NGRP+1)=NATOM
     ENDIF
#if KEY_NOST2==1
     IF(NST2 > 0)CALL WRNDIE(-4,'<PSFIO>', &
          'ST2 code is not compiled.')
#endif 
     !
#if KEY_CHEQ==1
     READ(IUNIT) MOLNTX
     READ(IUNIT) (MOLT(NAT+I),I=1,NATOM)
     DO I=1+NAT,NATOM
        MOLT(I)=MOLT(I)+MOLNT
     ENDDO
     MOLNT=MOLNT+MOLNTX
#endif 

#if KEY_LONEPAIR==1
     ! Read lone pair stuff
     READ(IUNIT,END=35) NUMLPX,NUMLPHX
     IF(NUMLPX > 0) THEN
        READ(IUNIT) (LPNHOST(I),LPHPTR(I),LPWGHT(I),LPVALUE(1,I), &
             LPVALUE(2,I),LPVALUE(3,I),I=NUMLP+1,NUMLP+NUMLPX)
        DO I=NUMLP+1,NUMLP+NUMLPX
           LPHPTR(I)=LPHPTR(I)+NUMLPH
        ENDDO
        READ(IUNIT) (LPHOST(I+NUMLPH),I=1,NUMLPHX)
        DO I=1,NUMLPHX
           LPHOST(I+NUMLPH)=LPHOST(I+NUMLPH)+NAT
        ENDDO
        NUMLP=NUMLP+NUMLPX
        NUMLPH=NUMLPH+NUMLPHX
     ENDIF
35   CONTINUE
#else /**/
     READ(IUNIT,END=35) NUMLPX,NUMLPHX
     IF(NUMLPX > 0) THEN
        CALL WRNDIE(-1,'<PSFIO>', &
             'Lonepairs present, but LONEPAIR code not compiled')
        DO I=1,NUMLPX
           READ(IUNIT,END=35)
        ENDDO
        READ(IUNIT,END=35)
     ENDIF
35   CONTINUE
#endif 

  ENDIF                              ! misplaced ENDIF fixed (c30a2)

#if KEY_PARALLEL==1
  CALL PSND4(NATOM,1)
  CALL PSNDC(ATYPE,NATOM)
  CALL PSND4(IAC,NATOM)
  CALL PSND8(CG,NATOM)
  CALL PSND8(AMASS,NATOM)
  CALL PSND4(IMOVE,NATOM)
  !
  CALL PSND4(NSEG,1)
  CALL PSND4(NICTOT,NSEG+1)
  CALL PSNDC(SEGID,NSEG)
  !
  CALL PSND4(NRES,1)
  CALL PSND4(IBASE,NRES+1)
  CALL PSNDC(RESID,NRES)
  CALL PSNDC(RES,NRES)
  !
  CALL PSND4(NGRP,1)
  CALL PSND4(IGPBS,NGRP+1)
  CALL PSND4(IGPTYP,NGRP)
  CALL PSND4(IMOVEG,NGRP)
  !
  CALL PSND4(NBOND,1)
  CALL PSND4(IB,NBOND)
  CALL PSND4(JB,NBOND)
  !
  CALL PSND4(NTHETA,1)
  CALL PSND4(IT,NTHETA)
  CALL PSND4(JT,NTHETA)
  CALL PSND4(KT,NTHETA)
  !
  CALL PSND4(NPHI,1)
  CALL PSND4(IP,NPHI)
  CALL PSND4(JP,NPHI)
  CALL PSND4(KP,NPHI)
  CALL PSND4(LP,NPHI)
  !
  CALL PSND4(NIMPHI,1)
  CALL PSND4(IM,NIMPHI)
  CALL PSND4(JM,NIMPHI)
  CALL PSND4(KM,NIMPHI)
  CALL PSND4(LM,NIMPHI)
  !
  CALL PSND4(NDON,1)
  CALL PSND4(IDON,NDON)
  CALL PSND4(IHD1,NDON)
  CALL PSND4(NACC,1)
  CALL PSND4(IACC,NACC)
  CALL PSND4(IAC1,NACC)
  !
  CALL PSND4(NNB,1)
  CALL PSND4(INB,NNB)
  CALL PSND4(IBLO,NATOM)
  !
  CALL PSND4(NST2,1)
  !
#if KEY_LONEPAIR==1 /*send_lonepairs*/
  CALL PSND4(NUMLP,1)
  CALL PSND4(NUMLPH,1)
  CALL PSND4(LPNHOST,NUMLP)
  CALL PSND4(LPHPTR,NUMLP)
  CALL PSND4(LPWGHT,NUMLP)
  CALL PSND8(LPVALUE,NUMLP*3)
  CALL PSND4(LPHOST,NUMLPH)
#endif /* (send_lonepairs)*/
  !
#endif 

  !yw...28-Jan-2003, Check if expanded format is needed
  if (oleng /= idleng) qextfmt=qxform()

  RETURN
END SUBROUTINE PSFRDR

SUBROUTINE PSFWRT(IUNIT,IPRINT,ICNTRL,TITLE,NTITL)
  !
  !     THIS ROUTINE WRITES THE PROTEIN STRUCTURE FILE MODULE. IPRINT
  !     DETERMINES WHETHER THE DATA IS TO BE PRINTED OR WRITTEN OUT
  !     UNFORMATTED.
  !
  !     binary format
  !     ===============
  !
  !     Authors: Robert Bruccoleri
  !     Bernie Brooks
  !
#if KEY_CHEQ==1
  use cheq,only:ech,eha,molnt,molt       
#endif
  use chm_kinds
  use dimens_fcm
  use psf
  use lonepr
  use aniso_fcm
  use stream
  use string
#if KEY_MMFF==1
  use mmffm
  use ffieldm
#endif 
  !
#if KEY_STRINGM==1 /*   VO string method */
  use machio, only: ifreeu
  use multicom_aux
#endif
  !
  use comand
  use param
  implicit none
  INTEGER ICNTRL(20)
  CHARACTER(len=*) TITLE(*)
  INTEGER NTITL,IUNIT,IPRINT
  !
  CHARACTER(len=4) :: PSFHDR='PSF '
  character(len=8) :: atmcode(natom)
  INTEGER NSEG1,IGRP,IRES,IW,I,J,K,IS,IQ
  INTEGER IJUNK
  integer iptlen
  character(len=20) fmt1,fmt2
  logical :: qxplor
  !
#if KEY_STRINGM==1 /*  VO stringm v */
  logical :: qstr
  common /replicaio/ qstr ! need global variable
#endif
  !
  NSEG1=NSEG+1
  ICNTRL(1)=18
  DO I=2,20
     ICNTRL(I)=0
  ENDDO

  QXPLOR=INDXA(COMLYN,COMLEN,'XPLO') > 0
  ! Always write an xplor psf unless specifically asked not to
  qxplor = qxplor .or. .not. (indxa(comlyn,comlen,'OLDPSF')>0)
  qextfmt=qxform()
  if (qextfmt) then
     icntrl(1)=30
     iptlen=10
     fmt1='(10X,I5,5X,10I10)'
     fmt2='(20X,10I10)'
  else
     iptlen=20
     fmt1='(10X,I5,5X,20I5)'
     fmt2='(20X,20I5)'
  endif
  !
#if KEY_STRINGM==1 /*  VO stringm i/o */
    qstr=iunit.le.size(ifreeu) ! protect against OOB on ifreeu
    if (qstr) qstr=((MOD(IFREEU(IUNIT),8).eq.0).and.(IFREEU(IUNIT).ne.0).and.(ME_LOCAL.eq.0))
#endif
  !
#if KEY_CHEQ==1
  !    SPatel/MFC Adding flag to mark this psf as including
  !      CHEQ data so the reader can know whether to read it or not.
  !
  icntrl(13)=1
#endif 

  if(qxplor) icntrl(20) = 1

  IF(IPRINT /= 0) GOTO 25
  if(ndrude>0) &
       call wrndie(-1,'<PSFWRT>','Binary PSF not supported with drude model')
  !
  !     The call to WRTITL should update the date on the title
  !
  CALL WRTITL(TITLE,NTITL,OUTU,2)
  !
  IF(IOLEV<0 &
#if KEY_STRINGM==1 /*  VO */
    &         .and..not.qstr &  
#endif
    &                         ) RETURN
  !
  !     the following initializes unused arrays:
  IJUNK=-9999
  K=0
  !
  WRITE(IUNIT) PSFHDR,ICNTRL
  WRITE(IUNIT) NTITL,NSEG
  CALL WRTITL(TITLE,NTITL,IUNIT,-2)
  WRITE(IUNIT) (NICTOT(J),NATOM,NBOND,NTHETA,NPHI, &
       NIMPHI,NNB,NDON,NACC,K,J=1,NSEG1)
  WRITE(IUNIT) (SEGID(I)(1:idleng),I=1,NSEG)
  WRITE(IUNIT) (RES(I)(1:idleng),RESID(I)(1:idleng),IJUNK,I=1,NRES)
  WRITE(IUNIT) (IBASE(I),I=1,NRES+1)
  ! From this point into future, we use the atc specification for atom 
  ! types and reconstruct the iac array. This allows one to increment atom
  ! type codes automatically. clb3, 08/05/14
  if(qxplor) then
     do i = 1, natom
        atmcode(i) = atc(iac(i))
     enddo
     WRITE(IUNIT) (ATYPE(I)(1:idleng),CG(I),AMASS(I),atmcode(i),IMOVE(I), &
#if KEY_CHEQ==1
          ECH(I),EHA(I),                             & 
#endif
          I=1,NATOM)
  else
     WRITE(IUNIT) (ATYPE(I)(1:idleng),CG(I),AMASS(I),IAC(I),IMOVE(I), &
#if KEY_CHEQ==1
          ECH(I),EHA(I),                             & 
#endif
          I=1,NATOM)
  endif
  WRITE(IUNIT) (IB(I),JB(I),I=1,NBOND)
  WRITE(IUNIT) (IT(I),JT(I),KT(I),I=1,NTHETA)
  WRITE(IUNIT) (IP(I),JP(I),KP(I),LP(I),I=1,NPHI)
  WRITE(IUNIT) (IM(I),JM(I),KM(I),LM(I),I=1,NIMPHI)
  WRITE(IUNIT) (IDON(I),IJUNK,IJUNK,IHD1(I),I=1,NDON)
  WRITE(IUNIT) (IACC(I),IAC1(I),IJUNK,I=1,NACC)
  WRITE(IUNIT) (INB(I),I=1,NNB)
  WRITE(IUNIT) (IBLO(I),I=1,NATOM)
  !
  !     WRITE GROUP STUFF
  WRITE(IUNIT) NGRP,NST2
  WRITE(IUNIT) (IGPBS(I),IGPTYP(I),IMOVEG(I),I=1,NGRP)

#if KEY_CHEQ==1
  WRITE(IUNIT) MOLNT
  WRITE(IUNIT) (MOLT(I),I=1,NATOM)
#endif 

#if KEY_LONEPAIR==1
  ! Write lone pair stuff
  WRITE(IUNIT) NUMLP,NUMLPH
  IF(NUMLP > 0) THEN
     WRITE(IUNIT) (LPNHOST(I),LPHPTR(I),LPWGHT(I),LPVALUE(1,I), &
          LPVALUE(2,I),LPVALUE(3,I),I=1,NUMLP)
     WRITE(IUNIT) (LPHOST(I),I=1,NUMLPH)
  ENDIF
#else /**/
  WRITE(IUNIT) 0,0
#endif 

  RETURN
  !
25 CONTINUE
  WRITE(IUNIT,30) ICNTRL
30 FORMAT(/10X,'PSF FILE MODULE'/10X,'CONTROL ARRAY :',20I5)
  CALL WRTITL(TITLE,NTITL,IUNIT,0)
  WRITE(IUNIT,40) NATOM,NBOND,NTHETA,NPHI,NIMPHI,NNB,NDON,NACC, &
       NRES,NSEG,NGRP,NST2
40 FORMAT(/9X,'NATOM NBOND NTHETA NPHI NIMPHI  NNB  NDON', &
       '  NACC  NRES  NSEG  NGRP  NST2'/7X,12I6)
  !
  WRITE(IUNIT,45)
45 FORMAT(//10X,'PARTITION OF SEGMENTS - NICTOT ARRAY AND SEGMENT', &
       ' IDENTIFIERS :', &
       /' SEGMENT NUMBER  ID  ',6X,'NRES',5X,'NATOM',5X,'NBOND', &
       4X,'NTHETA',6X,'NPHI',4X,'NIMPHI',7X,'NNB',6X,'NDON', &
       6X,'NACC',6X,'TYPE')
  DO I=1,NSEG
     WRITE(IUNIT,50) I,SEGID(I)(1:idleng),NICTOT(I+1)
50   FORMAT(1X,I10,6X,A,10I10)
  ENDDO
  !
  WRITE(IUNIT,105)
105 FORMAT(/10X,'ATOM CHARACTERISTICS :'/19X,'ATOM TYPE      ', &
       'CHARGE    ATOM CODE   COUNT OF    MOVEMENT FLAG', &
       '        MASS'/56X,'EXCLUSIONS'/)
  !
  IRES=1
  IGRP=1
  DO I=1,NATOM
     IF(IBASE(IRES)+1 == I) THEN
        WRITE(IUNIT,108) IRES,RESID(IRES)(1:idleng), &
             RES(IRES)(1:idleng),IBASE(IRES+1)
108     FORMAT(' RESIDUE',I4,2X,A,2X,A,8X,'  TO',I6)
        IRES=IRES+1
     ENDIF
     IF(IGPBS(IGRP)+1 == I) THEN
        WRITE(IUNIT,109) IGRP,IGPTYP(IGRP),IMOVEG(IGRP),IGPBS(IGRP+1)
109     FORMAT(10X,'GROUP',I4,'  TYPE',I3,' MOVE',I2,'  TO',I6)
        IGRP=IGRP+1
     ENDIF
     WRITE(IUNIT,110) I,ATYPE(I)(1:idleng),CG(I),IAC(I),IBLO(I), &
          IMOVE(I),AMASS(I)
110  FORMAT(20X,I5,2X,A,F12.4,I5,I6,I5,1PG14.6)
  ENDDO
  !
  WRITE(IUNIT,120)
120 FORMAT(/10X,'BOND ARRAY (BY COLUMNS) :')
  DO I=1,NBOND,iptlen
     IS=I
     IW=I+iptlen-1
     IF(IW > NBOND) IW=NBOND
     WRITE(IUNIT,fmt1) I,(IB(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (JB(J),J=IS,IW)
#if KEY_MMFF==1
     IF(FFIELD == MMFF) WRITE(IUNIT,fmt2) (BondType(J),J=IS,IW)
#endif 
  ENDDO
  !
  WRITE(IUNIT,130)
130 FORMAT(/10X,'THETA ARRAY (BY COLUMNS) :')
  DO I=1,NTHETA,iptlen
     IS=I
     IW=I+iptlen-1
     IF(IW > NTHETA) IW=NTHETA
     WRITE(IUNIT,fmt1) I,(IT(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (JT(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (KT(J),J=IS,IW)
  ENDDO
  !
  WRITE(IUNIT,140)
140 FORMAT(/10X,'PHI ARRAY (BY COLUMNS) :')
  DO I=1,NPHI,iptlen
     IS=I
     IW=I+iptlen-1
     IF(IW > NPHI) IW=NPHI
     WRITE(IUNIT,fmt1) I,(IP(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (JP(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (KP(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (LP(J),J=IS,IW)
  ENDDO
  !
  WRITE(IUNIT,150)
150 FORMAT(/10X,'IMPROPER TORSION ARRAY (BY COLUMNS) :')
  DO I=1,NIMPHI,iptlen
     IS=I
     IW=I+iptlen-1
     IF(IW > NIMPHI) IW=NIMPHI
     WRITE(IUNIT,fmt1) I,(IM(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (JM(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (KM(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (LM(J),J=IS,IW)
  ENDDO
  !

#if KEY_CMAP==1
  WRITE(IUNIT,155)
155 FORMAT(/10X,'CROSSTERM TORSION ARRAY (BY COLUMNS) :')
  DO I=1,NCRTERM,iptlen
     IS=I
     IW=I+iptlen-1
     IF (IW > NCRTERM) IW=NCRTERM
     WRITE(IUNIT,fmt1) I,(I1CT(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (J1CT(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (K1CT(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (L1CT(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (I2CT(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (J2CT(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (K2CT(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (L2CT(J),J=IS,IW)
  ENDDO
#endif 

  WRITE(IUNIT,165)
165 FORMAT(/10X,'HYDROGEN DONOR ARRAYS :')
  DO I=1,NDON,iptlen
     IS=I
     IW=I+iptlen-1
     IF(IW > NDON) IW=NDON
     WRITE(IUNIT,fmt1) I,(IDON(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (IHD1(J),J=IS,IW)
  ENDDO
  !
  WRITE(IUNIT,205)
205 FORMAT(/10X,'HYDROGEN ACCEPTOR ARRAYS :')
  DO I=1,NACC,iptlen
     IS=I
     IW=I+iptlen-1
     IF(IW > NACC) IW=NACC
     WRITE(IUNIT,fmt1) I,(IACC(J),J=IS,IW)
     WRITE(IUNIT,fmt2) (IAC1(J),J=IS,IW)
  ENDDO
  !
  WRITE(IUNIT,240)
240 FORMAT(/10X,'NON-BONDED EXCLUSION ARRAY :')
  DO I=1,NNB,iptlen
     IS=I
     IW=I+iptlen-1
     IF(IW > NNB) IW=NNB
     WRITE(IUNIT,fmt1) I,(INB(J),J=IS,IW)
  ENDDO
  !
#if KEY_LONEPAIR==1
  ! integers
  IF(NUMLP > 0) THEN
     WRITE(IUNIT,250) NUMLP,NUMLPH
250  FORMAT(/10X,'LONEPAIR ARRAYS: NUMLP,NUMLPH=',2I8)
     DO I=1,NUMLP
        WRITE(IUNIT,'(3I8,3X,L1,3G14.6)') I,LPNHOST(I),LPHPTR(I), &
             LPWGHT(I),LPVALUE(1,I),LPVALUE(2,I),LPVALUE(3,I)
     ENDDO
     DO I=1,NUMLPH,iptlen
        IS=I
        IW=I+iptlen-1
        IF(IW > NUMLPH) IW=NUMLPH
        WRITE(IUNIT,fmt1) I,(LPHOST(J),J=IS,IW)
     ENDDO
  ENDIF
#endif 
  !
  IF(NANISO > 0) THEN
     WRITE(IUNIT,260) NANISO
260  FORMAT(/10X,'ANISOTROPIC ALPHA: NANISO=',I8)
     DO I=1,NANISO
        WRITE(IUNIT,'(5I8,3x,3(a,f8.3,3x))') &
             i,LSTANI1(i),LSTANI2(i),LSTANI3(i),LSTANI4(i), &
             ' KPAR   ',K11(i), &
             ' KPERP ',K22(i), &
             ' KISO  ',K33(i)
     ENDDO
  ENDIF

  RETURN
END SUBROUTINE PSFWRT

SUBROUTINE SEQRDR(COMLYN,COMLEN,MXCMSZ,IUNIT,TITLE,NTITL,MAXTIT, &
     RES,NRES,RESID,MODE,ISTART,CHAIN,SEGIDX,NCHAIN,NSKIP,SKIP,NALI,ALIAS, &
     LATOM,LHETATM,LSEQRES,IFIRST)

  !
  !     Reads the sequence for a piece of chain from IUNIT.
  !     MODE = 0  The normal reading is done free field.
  !     2  Read the sequence from a protein system coordinate file.
  !     3  Read sequence and resid's from coordinate file
  !     4  as 3 for Brookhaven coordinate files
  !     5  QUANTA interface
  !     Authors: Robert Bruccoleri
  !     David States
  !     and others
  !
  use chm_kinds
  use exfunc
  use stream
  use string, only:decodi,encodi,nexti,indxa,nextwd,nexta4,filspc,cnvtuc,trima
#if KEY_STRINGM==1 /*   VO string method */
  use machio, only: ifreeu
  use multicom_aux
#endif
  use param_store, only: set_param

  implicit none

  CHARACTER(len=*) COMLYN
  INTEGER COMLEN,MXCMSZ,IUNIT,NTITL,MAXTIT
  CHARACTER(len=*) TITLE(*)
  CHARACTER(len=*) RES(*),RESID(*),SEGIDX
  CHARACTER(len=1) CHAIN
  INTEGER NRES,MODE,ISTART,NCHAIN,NSKIP,NALI,IFIRST
  CHARACTER(len=8) ALIAS(2,NALI),SKIP(NSKIP)
  LOGICAL LATOM,LHETATM,LSEQRES

  INTEGER WDLEN
  INTEGER, PARAMETER :: WDMAX=20
  CHARACTER(len=WDMAX) WD
  LOGICAL DONE,EOF
  !
  INTEGER, PARAMETER :: NCHECK=8
  CHARACTER(len=4) :: CHECK(NCHECK)= &
    (/'HIS ','HSD ','HSE ','HSP ','ASP ','GLU ','LYS ','TYR '/),IBLANK='    '
  INTEGER CHKNUM(NCHECK)
  INTEGER ERRCNT,RDRES,CNTRES,NLINE,CHKTOT,LRESO,LRES
  INTEGER START,I,J,K,DUMMY,JUNK
  CHARACTER(len=8) RESI,ARES     !yw
  !
  !     PDB reader variables
  INTEGER IJ, ISEQ, pdblen, ilen,SLEN
  CHARACTER(len=80) PDBLIN
  character(len=40) fmt1100
  CHARACTER(len=8) :: SQSEGID = ' '
  CHARACTER(len=8) ORID, SID, RID, ATOMIN !yw
  CHARACTER(len=1) ICODE
  real(chm_real) XIN, YIN, ZIN, WIN
  LOGICAL QATOM, lextfmt, QSEARCH,QTER
  !
#if KEY_STRINGM==1 /*  VO stringm v */
  integer :: oldiol
  logical :: qstr
  common /replicaio/ qstr ! need global variable
  qstr=iunit.le.size(ifreeu) ! protect against OOB on ifreeu
  if (qstr) qstr=((MOD(IFREEU(IUNIT),8).eq.0).and.(IFREEU(IUNIT).ne.0).and.(ME_LOCAL.eq.0))
  if (qstr) then
   oldiol=iolev
   iolev=1
  endif
#endif /* VO stringm */
  !
  IF(IOLEV > 0) THEN
     ERRCNT=0
     START=NRES+1
     CALL TRYORO(IUNIT,'FORMATTED')
     IF (MODE == 0) THEN
        CALL RDTITL(TITLE,NTITL,IUNIT,0)
        EOF=.FALSE.
        CALL RDCMND(COMLYN,MXCMSZ,COMLEN,IUNIT,EOF,.FALSE., &
             .TRUE.,'SEQRDR> ')
        IF (EOF) THEN
           IF(WRNLEV >= 2) WRITE(OUTU,200)
200        FORMAT(' ***** Error in SEQRDR ***** End of file when', &
                ' reading a sequence.')
           CALL DIEWRN(0)
           GOTO 900
           !
        ENDIF
        CALL NEXTWD(COMLYN,COMLEN,WD,WDMAX,WDLEN)
        IF(WDLEN == 0) THEN
           IF(WRNLEV >= 2) WRITE(OUTU,210)
210        FORMAT(' ***** Error in SEQRDR ***** Number of residues', &
                ' not specified.')
           CALL DIEWRN(0)
           GOTO 900
           !
        ENDIF
        RDRES=DECODI(WD,WDLEN)
        IF(RDRES == 0 .AND. WRNLEV >= 2) WRITE(OUTU,215)
215     FORMAT(' ***** Warning from SEQRDR ***** A zero residue count', &
             ' was specified.'/' You may have neglected to specify a', &
             ' a residue count for your sequence, '/' so the first', &
             ' residue may have been eaten.')
        CNTRES=0
        DONE=.FALSE.
        DO WHILE (.NOT.DONE)
           WDLEN=1
           DO WHILE (WDLEN /= 0)
              CALL NEXTWD(COMLYN,COMLEN,WD,WDMAX,WDLEN)
              IF (WDLEN /= 0) THEN
                 IF (WDLEN > 8) THEN
                    IF(WRNLEV >= 2) WRITE(OUTU,220)
220                 FORMAT(' ***** Error in SEQRDR ***** Residue name longer than', &
                         ' eight characters.')
                    ERRCNT=ERRCNT+1
                 ENDIF
                 CALL FILSPC(WD,8,WDLEN)
                 CNTRES=CNTRES+1
                 RES(NRES+CNTRES)=WD
              ENDIF
           ENDDO
           IF (RDRES > 0.AND.CNTRES >= RDRES) THEN
              DONE=.TRUE.
              IF(CNTRES > RDRES .AND. WRNLEV >= 2) WRITE(OUTU,230)  &
                   CNTRES,RDRES
230           FORMAT(' ***** Warning from SEQRDR ***** Number of residues', &
                   ' read = ',I10/' is greater than number specified =',I10)
           ELSE
              CALL RDCMND(COMLYN,MXCMSZ,COMLEN,IUNIT,EOF,.FALSE., &
                   .TRUE.,'SEQRDR> ')
              IF(EOF.OR.COMLEN == 0) DONE=.TRUE.
           ENDIF
        ENDDO
        !         ENDDO: DO WHILE (.NOT.DONE)
        NRES=NRES+CNTRES
     ELSEIF (MODE == 2) THEN
        CALL RDTITL(TITLE,NTITL,IUNIT,0)
        !yw++ 28-Jan-2003 use PDBLIN to process a line
        !yw       READ(IUNIT,1100) NLINE
        !yw 1100  FORMAT(2I5,1X,A4,40X,1X,A4)
        READ(IUNIT,'(A)') PDBLIN
        pdblen=len(pdblin)
        nline=nexti(pdblin,pdblen)
        lextfmt=indxa(pdblin,pdblen,'EXT') > 0
        if (nline >= 100000) lextfmt=.true.
        if (lextfmt) then
           fmt1100='(2I10,2X,A8,70X,2X,A8,2X,A8)'
           ilen=8
        else
           fmt1100='(2I5,1X,A4,35X,1X,A4,1X,A4)'
           ilen=4
        endif
        !yw--
        LRESO=-999
        DO I=1,NLINE
           READ(IUNIT,fmt1100) JUNK,LRES,RESI,SID
           IF(SEGIDX ==' ' .OR. SEGIDX(1:ilen) == SID)THEN
             IF (LRES /= LRESO) THEN
                NRES=NRES+1
                RES(NRES)=RESI
                LRESO=LRES
             ENDIF
           ENDIF
        ENDDO
        IF (reallow) THEN  
           REWIND IUNIT
        ENDIF              
     ELSEIF (MODE == 3) THEN
        CALL RDTITL(TITLE,NTITL,IUNIT,0)
        !yw       READ(IUNIT,1100) NLINE
        READ(IUNIT,'(A)') PDBLIN
        pdblen=len(pdblin)
        nline=nexti(pdblin,pdblen)
        lextfmt=indxa(pdblin,pdblen,'EXT') > 0
        if (nline >= 100000) lextfmt=.true.
        if (lextfmt) then
           fmt1100='(2I10,2X,A8,70X,2X,A8,2X,A8)'
           ilen=8
        else
           fmt1100='(2I5,1X,A4,35X,1X,A4,1X,A4)'
           ilen=4
        endif
        LRESO=-999
        DO I=1,NLINE
           READ(IUNIT,fmt1100) JUNK,LRES,RESI,SID,ARES
           IF(SEGIDX ==' ' .OR. SEGIDX(1:ilen) == SID)THEN
             IF (LRES /= LRESO) THEN
                DO J=START,NRES
                   IF(ARES == RESID(J)) THEN
                      IF(WRNLEV >= 2) WRITE(OUTU,1083) J,NRES+1,ARES(1:ilen)
                      CALL DIEWRN(0)
                   ENDIF
                ENDDO
          
                !
                NRES=NRES+1
                RES(NRES)=RESI
                IF(ARES == IBLANK) THEN
                   IF(WRNLEV >= 2) WRITE(OUTU,1089) JUNK,LRES,RESI(1:ilen)
                   CALL DIEWRN(0)
                   CALL ENCODI(NRES-ISTART,ARES,8,DUMMY)
                ENDIF
                RESID(NRES)=ARES
                LRESO=LRES
             ENDIF
           ENDIF
        ENDDO
        IF (reallow) THEN  
           REWIND IUNIT
        ENDIF              
1083    FORMAT(' ** ERROR IN SEQRDR ** REPEATED RESID FOR ',2I5,1X,A)
1089    FORMAT(' **** ERROR **** NO RESID SPECIFIED.',2I5,1X,A)
     ELSEIF (MODE == 4) THEN
        !
        !         PDB format (read with RESID)
        !         see also module COORIO !!!
        !
        ORID=' '
        !         read PDB title
        NTITL=0
        DONE=.FALSE.
        DO WHILE (.NOT.DONE)
           READ(IUNIT,'(A)',END=61,ERR=61) PDBLIN
           SLEN=LEN(PDBLIN)
           CALL CNVTUC(PDBLIN,SLEN)
           DONE=PDBLIN(1:6) /= 'REMARK'
           IF (.NOT.DONE) THEN
              NTITL=NTITL+1
              TITLE(NTITL)=PDBLIN(8:)
              DONE= (NTITL >= MAXTIT)
           ENDIF 
        ENDDO
        CALL WRTITL(TITLE,NTITL,OUTU,+1)
        !
        !         now read other records
        IF(LSEQRES)THEN
        ! Use SEQRES-records instead of the ATOM/HETATM records
        ! How are residue insertion codes handled in SEQRES? Not at all?
        ! Find first SEQRES  
          DO WHILE(PDBLIN(1:6) /= 'SEQRES')
            READ(IUNIT,'(A)',END=61,ERR=61) PDBLIN
            SLEN=LEN(PDBLIN)
            CALL CNVTUC(PDBLIN,SLEN)
           ENDDO
           
          IF(CHAIN /= ' ') THEN
            DO WHILE(CHAIN /= PDBLIN(12:12) )
              READ(IUNIT,'(A)',END=61,ERR=61) PDBLIN
              SLEN=LEN(PDBLIN)
              CALL CNVTUC(PDBLIN,SLEN)
            ENDDO

            I=IFIRST-1
            WRITE(OUTU,'(A)') 'Reading SEQRES sequence for CHAIN '//CHAIN
            WRITE(OUTU,'(A,I6)') 'Assigning RESIDs starting from',IFIRST
            DO WHILE(CHAIN == PDBLIN(12:12))
              PDBLIN(1:19)=' '
              CALL TRIMA(PDBLIN,SLEN)
              DO WHILE(SLEN > 0) 
                ! insert new residue
                NRES=NRES+1
                I=I+1
                RES(NRES)=NEXTA4(PDBLIN,SLEN)
                CALL ENCODI(I,RESID(NRES),IDLENG, DUMMY)
              ENDDO
              READ(IUNIT,'(A)',END=61,ERR=61) PDBLIN
              SLEN=LEN(PDBLIN)
              CALL CNVTUC(PDBLIN,SLEN)
            ENDDO            
            GOTO 61
          ELSEIF( NCHAIN > 0 )THEN
            I=1
            CHAIN=PDBLIN(12:12)
            DO WHILE(I < NCHAIN)
              DO WHILE(PDBLIN(12:12) == CHAIN)
                READ(IUNIT,'(A)',END=61,ERR=61) PDBLIN
                SLEN=LEN(PDBLIN)
                CALL CNVTUC(PDBLIN,SLEN)
              ENDDO
              CHAIN=PDBLIN(12:12)
              I=I+1
            ENDDO 
            WRITE(OUTU,'(A,I5,A)') &
               'Reading SEQRES sequence for CHAIN #',NCHAIN,' (CHAIN ID '//CHAIN//')'
            SQSEGID=CHAIN(1:)
            I=IFIRST-1
            DO WHILE(CHAIN == PDBLIN(12:12))
              PDBLIN(1:19)=' '
              CALL TRIMA(PDBLIN,SLEN)
              DO WHILE(SLEN > 0) 
                ! insert new residue
                NRES=NRES+1
                I=I+1
                RES(NRES)=NEXTA4(PDBLIN,SLEN)
                CALL ENCODI(I,RESID(NRES),IDLENG, DUMMY)
              ENDDO
              READ(IUNIT,'(A)',END=61,ERR=61) PDBLIN
              SLEN=LEN(PDBLIN)
              CALL CNVTUC(PDBLIN,SLEN)
            ENDDO            
            CHAIN=' '
            GOTO 61 
          ELSE
        ! No CHAIN is specified, just gulp up all the residues
            WRITE(OUTU,'(A)') 'Reading SEQRES sequence'
            I=0
            DO WHILE(PDBLIN(1:6) == 'SEQRES')
              PDBLIN(1:19)=' '
              CALL TRIMA(PDBLIN,SLEN)
              DO WHILE(SLEN > 0) 
                ! insert new residue
                NRES=NRES+1
                I=I+1
                RES(NRES)=NEXTA4(PDBLIN,SLEN)
                CALL ENCODI(I,RESID(NRES),IDLENG, DUMMY)
              ENDDO
              READ(IUNIT,'(A)',END=61,ERR=61) PDBLIN
              SLEN=LEN(PDBLIN)
              CALL CNVTUC(PDBLIN,SLEN)
            ENDDO            
            CHAIN=' '
            GOTO 61
          ENDIF
! Done reading SEQRES records, and should not ever be here
        ENDIF
        QTER=.FALSE.
        IF(NCHAIN > 1) THEN  ! skip the first NCHAIN-1 TER-terminated chains
           DONE=.FALSE.
           I=0
           DO WHILE(.NOT. DONE)
             ! I don't think we need to worry about END records here? LNI
             IF(PDBLIN(1:3) == 'TER')THEN
               I=I+1
               IF(I == NCHAIN-1) DONE=.TRUE.
             ENDIF  
             READ(IUNIT,'(A)',END=61,ERR=61) PDBLIN
             SLEN=LEN(PDBLIN)
             CALL CNVTUC(PDBLIN,SLEN)
           ENDDO
        ENDIF
        IF(NCHAIN > 0) THEN 
           IF(PRNLEV > 2 ) &
             WRITE(OUTU,'(A,I5,A)') 'Reading sequence for CHAIN # ',I+1, &
             ' (CHAIN ID '//PDBLIN(12:12)//')'
           SQSEGID=PDBLIN(12:12)
           QTER=.TRUE. ! We now have to respect the next TER-line
        ENDIF        
        !
        ! Now we should be at the specified CHAIN number (or at the beginning)
        QSEARCH= .NOT.(CHAIN == ' ' .AND. SEGIDX == ' ')
7777    CONTINUE
        IF (PDBLIN(1:3) == 'END' ) THEN
           GOTO 61
        ELSEIF (PDBLIN(1:3) == 'TER' .AND. QTER) THEN
           GOTO 61
        ELSEIF (LATOM .AND. PDBLIN(1:4) == 'ATOM') THEN
!LNI, December 2013. There are unusual things in HETATM records, so we may
! want to ignore those when reading the sequence from a PDB file.
        ELSEIF (LHETATM .AND. PDBLIN(1:4) == 'HETA') THEN
        ELSE 
           ! keep reading until reaching ATOM/HETATM or END/TER
           QATOM=.FALSE.
           DO WHILE (.NOT.QATOM)
              READ(IUNIT,'(A)',END=61,ERR=61) PDBLIN
              SLEN=LEN(PDBLIN)
              CALL CNVTUC(PDBLIN,SLEN)
              IF (PDBLIN(1:3) == 'END' .OR. QTER.AND.PDBLIN(1:3) == 'TER') GOTO 61
              QATOM=( (LATOM .AND. PDBLIN(1:4) == 'ATOM') .OR. &
                       (LHETATM .AND. PDBLIN(1:4).EQ.'HETA') )
           ENDDO
        ENDIF
        !
        ! If CHAIN/SEGI has been specified skip lines that do not match CHAIN/SEGI   
        IF( .NOT. QSEARCH .OR.  &
         (CHAIN /=  ' ' .AND. CHAIN  == PDBLIN(22:22)) .OR. &
         (SEGIDX /= ' ' .AND. SEGIDX == PDBLIN(73:76)) ) THEN
           !         process ATOM line
           ! Format correction and residue insertion check. L.Nilsson, November 07
           ! Note: This still deviates sligthly from the official PDB, in that four
           ! characters are allowed for the residue name (RESI), where PDB only allows three.
           ! Previous format: (6X,I5,1X,A4,1X,A4,2X,A4,3X,3F8.3,6X,F6.2,6X,A4)
           !
           IF(QSEARCH .AND. .NOT. QTER) THEN
              IF(PRNLEV > 2 )THEN
                IF( CHAIN == PDBLIN(22:22) )THEN
                  WRITE(OUTU,'(A)') 'Reading sequence for CHAIN '//CHAIN
                  SQSEGID=CHAIN
                ELSE
                  WRITE(OUTU,'(A)') 'Reading sequence for SEGID '//SEGIDX
                  SQSEGID=SEGIDX
                ENDIF
              ENDIF
           ! chains are contiguous so if we found it we can turn on QTER, so
           ! that a TER record will terminate reading 
              QTER=.TRUE.
           ENDIF
           READ(PDBLIN, &
                '(6X,I5,1X,A4,1X,A4,1X,A4,A1,3X,3F8.3,6X,F6.2,6X,A4)') &
                ISEQ,ATOMIN,RESI,RID,ICODE,XIN,YIN,ZIN,WIN,SID
           !
           ! make RID left-justified
           IJ=4
           CALL TRIMA(RID,IJ)
           ! and append the insertion code 
           IF(ICODE  /=  ' ') RID=RID(1:IJ)//ICODE
           ! make RESI left-justified
           IJ=4
           CALL TRIMA(RESI,IJ)
           !
           ! test whether it is a new residue
           IF(ORID /= RID)THEN
              !
              IF(ICODE /= ' ' .AND. WRNLEV >= 5) WRITE(OUTU,'(A)') &
                   'SEQRDR> Residue insertion at RESID "'//RID(1:idleng)//'"'
              ! test for duplications
              DO J=START,NRES
                 IF(RID == RESID(J)) THEN
                    IF(WRNLEV >= 2) WRITE(OUTU,1083)J,NRES+1,RID(1:idleng)
                    CALL DIEWRN(0)
                 ENDIF
              ENDDO
              !
              ! insert new residue
              NRES=NRES+1
              RES(NRES)=RESI
              RESID(NRES)=RID
              ORID=RID
           ENDIF
        ENDIF
        !
        !         read next PDB line and go back
        READ(IUNIT,'(A)',ERR=61,END=61) PDBLIN
        SLEN=LEN(PDBLIN)
        CALL CNVTUC(PDBLIN,SLEN)
        GOTO 7777
        !
        !         PDB termination label:
61      CONTINUE
        IF (reallow) THEN  
           REWIND IUNIT
        ENDIF              
     ELSE 
        IF(WRNLEV >= 2) WRITE(OUTU,2131)
2131    FORMAT(' *****  FATAL  ***** SEQRDR illegal mode.')
        CALL DIEWRN(-3)
        GOTO 900
        !
     ENDIF
     IF(MODE == 4) THEN
       ! Process ALIAS
       IF(NALI >0) THEN
         DO I=START,NRES
           DO J=1,NALI
             IF(RES(I) == ALIAS(1,J)) THEN
               RES(I)=ALIAS(2,J)
               EXIT
             ENDIF
           ENDDO
         ENDDO 
       ENDIF
       ! Process SKIP
       IF(NSKIP >0) THEN
         I=START 
         DO WHILE (I<=NRES)
           DO J=1,NSKIP
             IF(RES(I) == SKIP(J)) THEN
               ! Remove this entry and contract from end of lists
               DO K=I,NRES-1
                 RES(K)=RES(K+1)
                 RESID(K)=RESID(K+1)
               ENDDO
               NRES=NRES-1
               ! We will have to reprocess position I
               I=I-1 
               EXIT
             ENDIF
           ENDDO
           I=I+1
         ENDDO 
       ENDIF
     ENDIF
     !
     CNTRES=NRES-START+1
     IF(PRNLEV >= 2) WRITE(OUTU,65) CNTRES, &
          (RES(I)(1:idleng),I=START,NRES)
65   FORMAT(/10X,'RESIDUE SEQUENCE -- ',I5,' RESIDUES',/,(10X,20A))
     CHKTOT=0
     DO J=1,NCHECK
        CHKNUM(J)=0
     ENDDO
     DO I=START,NRES
        DO J=1,NCHECK
           IF (RES(I) == CHECK(J)) THEN
              CHKNUM(J)=CHKNUM(J)+1
              CHKTOT=CHKTOT+1
           ENDIF
        ENDDO
     ENDDO
     IF (CHKTOT > 0 .AND. WRNLEV >= 2) &
          WRITE(OUTU,67) CHKTOT,(CHECK(J),CHKNUM(J),J=1,NCHECK)
67   FORMAT(' ***** Message from SEQRDR ***** THE SYSTEM CONTAINS',I3, &
          ' TITRATABLE GROUPS'/ &
          ' THE USER MUST PREDETERMINE THE PROTONATION STATE', &
          ' THROUGH THE SEQUENCE AND RTF'/(10(1X,A4,'-',I3,1X)))
     IF (ERRCNT > 0) THEN
        IF(WRNLEV >= 2) WRITE(OUTU,250)
250     FORMAT(' Execution terminated due to errors in reading ', &
             'the sequence.')
        CALL DIEWRN(0)
        GOTO 900
     ENDIF
     !
     !       DEFINE RESID ARRAY (NOT FOR RESID OR PDB)
     !
     IF(MODE < 3) THEN
        DO I=START,NRES
           CALL ENCODI(I-ISTART,RESID(I),8,DUMMY)
        ENDDO
     ENDIF
     !
900  CONTINUE
  ENDIF
  !
#if KEY_STRINGM==1 /*  VO : restore iolev */
  if (qstr) iolev=oldiol 
#endif
  !
#if KEY_PARALLEL==1
!-- ##IF .not.PARALLEL .not.ENSEMBLE
  CALL PSND4(CNTRES,1)
  CALL PSND4(NRES,1)
  CALL PSNDC(SQSEGID,1)
  CALL PSNDC(RES,NRES)
  CALL PSNDC(RESID,NRES)
!-- ##ELSE
!-- !MFC-- Needs to go away if separate io works.
!-- !MFC--- I think this is wrong, parallel ensemble should never send to world
!--   CALL PSND4_WORLD(NRES,1)
!--   CALL PSNDC_WORLD(RES,NRES)
!--   CALL PSNDC_WORLD(RESID,NRES)
!-- ##ENDIF
#endif 
     CALL set_param('SQNRES',CNTRES)
     CALL set_param('SQSEGID', SQSEGID)
  !
  RETURN
END SUBROUTINE SEQRDR

subroutine fillint4(ifill,isrc,nfill,index)
  use chm_kinds
  implicit none
  integer nfill,index
  integer isrc(*)
  integer*4 ifill(*)
  integer i

  do i=1,nfill
     ifill(index)=isrc(i)
     index=index+1
  enddo
  return
end subroutine fillint4

subroutine extractint4(dest,src,index,num)
  use chm_kinds
  implicit none
  integer index,num
  integer*4 src(index+num-1)
  integer dest(num)
  integer i

  do i=0,num-1
     dest(i+1)=src(index)
     index=index+1
  enddo
  return
end subroutine extractint4

subroutine extractint8(dest,src,index,num)
  use chm_kinds
  implicit none
  integer index,num
  integer*8 src(index+num-1)
  integer dest(num)
  integer i

  do i=0,num-1
     dest(i+1)=src(index)
     index=index+1
  enddo
  return
end subroutine extractint8

subroutine writeint4(iunit,iarray,n)
  use chm_kinds
  implicit none
  integer iunit,n
  integer*4 iarray(n)
  integer i

  write(iunit)(iarray(i),i=1,n)
  return
end subroutine writeint4

subroutine readint4(iunit,iarray,n)
  use chm_kinds
  implicit none
  integer n
  integer*4 iarray(n)
  integer iunit,i
  read(iunit)(iarray(i),i=1,n)
  return
end subroutine readint4

subroutine readint8(iunit,iarray,n)
  use chm_kinds
  implicit none
  integer n
  integer*8 iarray(n)
  integer iunit,i
  read(iunit)(iarray(i),i=1,n)
  return
end subroutine readint8

