#if KEY_REPDSTR2==1 
!----------------------------------------------------------------------
!         REPD ENSNENSEM
!----------------------------------------------------------------------
subroutine repd_ensnensem(nrep,qverb)
  !-----------------------------------------------------------------------
  ! Set the number of replicas NENSEM
  !-----------------------------------------------------------------------
  use chm_kinds
  use repd_ensemble
  use stream
  use dimens_fcm
  use string
  use memory
  use parallel
  use parallel_groups
  use mpi
  use param_store,only:set_param
  implicit none

  integer,intent(in) :: nrep
  logical,intent(in) :: qverb
  integer i,j,k,l,node, startnode, endnode, numnod_orig
  integer ens_group! , old_comm_charmm
  integer,save :: orig_comm_charmm,orig_charmm_group,orig_mynod
  integer new_charmm_group, new_comm_charmm, new_comm_charmm_index
  integer master_group,     master_comm,      comm_master_index
  integer ierror,status, everb
  integer charmm_group

  integer, allocatable, dimension(:) :: nodes, masternodes
  logical, save :: already_called = .false.
!  character(len=14),dimension(0:maxens) :: chmout
  character(len=80) :: chint
  logical :: eof

  orig_comm_charmm=comm_charmm
  orig_mynod=mynod
  repd_ensemble_verbose=qverb
  !-------- For purposes of being able to turn ensemble off, need a new communicator
  !-------- Make a new comm_charmm complete copy of current comm_charmm (identical)
  !-------- When turning ensemble off, need to kill all child communicators
  !-------- Killing this copy with the machinery in paralgroups, we can eliminate
  !--------   all child communicators and be left where we started before ensemble
  call set_repd_nensem_and_verbosity         !contained sub
  call setup_repd_comm_charmm     !contained sub
  call create_repd_communicators  !contained sub
  call open_chmout               !contained sub
  call repd2_setup_stream_input
  !---- Reset the printlevel/warnlevel/iolev for new group
  if(mynod /= 0) plnod0=0
  if(mynod /= 0) then
     prnlev = -1
     wrnlev = -5
     iolev  = -1
  else
     call mpi_bcast(prnlevr,1,mpi_integer,0,comm_master,ierror)
     call mpi_bcast(wrnlevr,1,mpi_integer,0,comm_master,ierror)
     call mpi_bcast( iolevr,1,mpi_integer,0,comm_master,ierror)
     call mpi_bcast(plnod0r,1,mpi_integer,0,comm_master,ierror)
  endif

  call set_param('WHOIAM',whoiam)
  call set_param('NENSEM',nensem)
  call set_param('MYNODE',mynod)
  call set_param('NUMNODE',numnod)
  do i=1,maxnode
     inode(i)=mod(i,numnod)
  enddo
  call rensprint("   ENSNENSEM>> Testing outu unit","6")
  call cube(mynod,numnod,ippmap)

!  call ens_setup_stream_input()
  return

contains
  !----------------------------------------------------------------------
  subroutine setup_repd_comm_charmm
    !--- Save original comm_charmm stuff to reset when finishing repd
    old_comm_charmm=comm_charmm  ! dont know what this is for
    old_mynod=mynod              ! dont know what this is for 
    old_numnod=numnod


    !--- Make a new group and communicator as a duplicate of the current
    !---   comm_charmm and charmm_group
    !--- get the group number of comm_charmm, call it orig_charmm_group
    call mpi_comm_group(COMM_CHARMM, orig_charmm_group, ierror)

    !--- fill a node table to make a new group and commumicator
    call chmalloc('ensemble.src','ensnensem','nodes',numnod,intg=nodes)
    do i=1,numnod
       nodes(i)=i-1
    enddo

    !--- form a new group (duplicate) as a subset of the nodes
    !---  in orig_charmm_group, contains all the nodes.
    !---  Call it new_charmm_group
    call mpi_group_incl(orig_charmm_group, numnod, nodes, new_charmm_group, ierror)

    !--- Make a communicator for the new group
    call mpi_comm_create(orig_comm_charmm, new_charmm_group, new_COMM_CHARMM, ierror)

    !------------------------------------------------------------------------
    ! mfc will use this for embedded replicas later
    !    call comm_save(orig_comm_charmm,new_comm_charmm,new_charmm_group, &
    !         mynod,numnod,comm_charmm_index,new_comm_charmm_index)
    !    comm_charmm_index = new_comm_charmm_index
    !------------------------------------------------------------------------
 
    !--- clean up and print out result if verbose
    call chmdealloc('ensemble.src','ensnensem','nodes',numnod,intg=nodes)
    write(chint,'(i12)')new_comm_charmm
    call rensprint("setup_repd_comm_charmmcomm_ensemble",chint)

    !--- global values from original communicator: mynodg, numnodg
    mynodg = mynod
    numnodg = numnod
    numnod_orig = numnod

    !--- numnod will now be number of nodes per ensemble replica
    numnod = numnod/nensem
    !--- whoiam is which replica in the current scheme this processor is in.
    whoiam = mynod/numnod

    !---  Sanity check and get ready to make replica communicators and groups
    call repd_check_and_allocate_node_tables   !internal subroutine

    !------------------------------------------------------------------------
    !--- mfc next few commented lines are for future embedded replicas
    !    current_ens_layer = ensemble_layers
    !    comm_ensemble(current_ens_layer) = new_comm_charmm
    !    comm_ens_index(current_ens_layer) = new_comm_charmm_index
    !------------------------------------------------------------------------

    !------------ Now set the comm_charmm to this child copy of the original
    comm_charmm = new_comm_charmm

    call mpi_comm_group(COMM_CHARMM, ens_group, ierror)

    return
  end subroutine setup_repd_comm_charmm

  !----------------------------------------------------------------------
  subroutine create_repd_communicators
    !-- Set up groups. Make group for each set of processors (master group and member rep).
    call setup_repd_groups          !contained sub

    !--- Create COMM_CHARMM communicator
    !---    Redefine comm_charmm to be the one replica set of processors that
    !---    this processor belongs to, as defined above in setup_repd_groups

    call mpi_comm_create(comm_charmm, new_charmm_group, new_COMM_CHARMM, ierror)
    !--- mfc for embedded replicas:
    !    call comm_save(old_comm_charmm,new_comm_charmm,new_charmm_group, &
    !         mynod,numnod,comm_charmm_index,new_comm_charmm_index)

    write(chint,'(i4)')new_comm_charmm
    call rensprint("setup_ens_comm_charmm comm_charmm rep",chint)

    !--- Create COMM_MASTER communicator is the communicator 
    !---        between the master node only of each replica
    call mpi_comm_create(comm_charmm, master_group, COMM_MASTER, ierror)
    write(chint,'(i4)')comm_master
    call rensprint("setup_ens_comm_charmm comm_master",chint)

    ensmasternod=-1
    if(lmasternode) then
       call mpi_comm_rank(comm_master,ensmasternod,ierror)
       print *,mynodg," ensmasternod ",ensmasternod
    endif

    call chmdealloc('ensemble.src','ensnensem','nodes',numnod,intg=nodes)
    call chmdealloc('ensemble.src','ensnensem','masternodes',nensem,intg=masternodes)

    comm_charmm = new_comm_charmm

    if (numnod > 1) then
       call mpi_comm_rank(comm_charmm, mynod, ierror)
    else
       mynod = 0
    endif

    write(chint,'(2i4)')mynod,comm_charmm
    call rensprint("ENSNENSEM create_ens_communicators mynod, comm_charmm rep",chint)

    if (mynod /= 0) lslavenode = .true.
    mynodp = mynod + 1
    plnod0=prnlev

    mastermaster = (whoiam == 0 .and. mynod == 0)
    return
  end subroutine create_repd_communicators

  !----------------------------------------------------------------------
  subroutine setup_repd_groups
    !-- Set up groups. Make group for each set of processors.
    !-- First fill arrays with list of processors for each group and the master group
    node = 0
    lmasternode=.false.
    do i=1,nensem
       startnode = node
       do j=1,numnod
          nodes(j) = node
          node = node + 1
       enddo
       endnode = node - 1
       !-- Make a new group for a replica if this node is part 
       !--  of the group responsible for that replica
       ! include this mynod to group charmm_group, 
       !     if startnode <= mynod <= endnode
       if (mynod >= startnode .and. mynod <= endnode ) then
          call mpi_group_incl(ens_group, numnod, nodes, new_charmm_group, ierror)
       endif
       masternodes(i) = startnode
       lmasternode = lmasternode .or. (startnode == mynod)
    enddo
    master_node = 0
    
    call mpi_group_incl(ens_group, nensem, masternodes, master_group, ierror)
    return
  end subroutine setup_repd_groups

  !----------------------------------------------------------------------
  subroutine open_chmout
    logical opn
return
    !--- fill in default stdout file names------
    outloop: do i=0,9
       do j = 0,9
          do k = 0,9
             l = 100*i + 10*j + k
!             if( l > maxens)exit outloop
             write(chmout(l),'("charmm.out.",3i1)')i,j,k
          enddo
       enddo
    enddo outloop

    if(mynod == 0 .and. whoiam > 0 ) then
       call rensprint("   ENSEMBLE>> opening file ",chmout(whoiam))
       inquire(unit=6,opened=opn)
       if(opn)then
          !          inquire(unit=6,name=chmout(whoiam))
          call rensprint("   ENSEMBLE>> unit 6 already open: ",chmout(whoiam))
          close(6)
       endif
       open (unit=6, file=chmout(whoiam),status="REPLACE")
       call rensprint("   ENSEMBLE>> opened file ",chmout(whoiam))
    endif
    !    write(chint,'("6",3i6)')whoiam,mynod,numnod
    !    call ensprint("   ENSNENSEM>> Testing outu unit",chint(1:len_trim(chint)))
    if(iolev > 0 )write(outu,*)"Testing outu ",whoiam,mynod
    return
  end subroutine open_chmout

  !----------------------------------------------------------------------
  subroutine repd_check_and_allocate_node_tables
    if (numnod*nensem /= numnod_orig) then
       write(outu,*)"ENSEMBLE>> ERROR numnod,nensem ",numnod,nensem
       CALL WRNDIE(-3,'<ENSNENSEM>', &
            'nensem must be factor of number of nodes (mod(numnod,nensem)=0).')
    endif
    
    IF (numnod < 1) CALL WRNDIE(-3,'<ENSNENSEM>', 'NUMNOD CANNOT BE SMALLER THAN ONE.')

    !--- mfc for opening default chmout files, can probably be deleted.
    !    IF (nensem > maxens) then
    !       write(outu,*)"ENSEMBLE>  Maxens = ",maxens
    !       CALL WRNDIE(-3,'<ENSNENSEM>', &
    !            'Need to redimension chmout for more than maxens ensemble members')
    !    endif
    !--- mfc embedded replicas future development
    !    ensemble_layers = ensemble_layers + 1
    !    if(ensemble_layers > maxensemble_layers) &
    !         CALL WRNDIE(-3,'<ENSNENSEM>ensemble.src', &
    !         'Too many ensemble layers (ensemble nested too many times).')
         
    call chmalloc('ensemble.src','ensnensem','nodes',numnod,intg=nodes)
    call chmalloc('ensemble.src','ensnensem','masternodes',nensem,intg=masternodes)
  end subroutine repd_check_and_allocate_node_tables

  !----------------------------------------------------------------------
  subroutine set_repd_nensem_and_verbosity
    if (iolev >= 0) then
       nensem = nrep
       if (nensem < 1) then
          write(outu,*)" ERROR number of reps < 1    nrep = ",nrep
          call wrndie(-3,'<ENSNENSEM>', 'Invalid number of ensemble copies.')
       endif
    endif
    everb=0
    if(repd_ensemble_verbose)everb=1
    call mpi_bcast(nensem,1,mpi_integer,0,comm_charmm,status)
    call mpi_bcast(everb ,1,mpi_integer,0,comm_charmm,status)
    repd_ensemble_verbose = (everb == 1)
  end subroutine set_repd_nensem_and_verbosity

end subroutine repd_ensnensem

!----------------------------------------------------------------------
!          REPD2_ENSFIN
!----------------------------------------------------------------------
subroutine repd2_ensfin
  !----------------------------------------------------------------------
  ! clean up on exit
  !----------------------------------------------------------------------
  use chm_kinds
  use dimens_fcm
  use stream
  use repd_ensemble
  use psf
  use mpi
  use parallel,only:comm_charmm,mynod,numnod
!  use parallel_groups,only:current_comm_index,allcomms
  implicit none

  integer ierror
  character(len=ensbfl) tbuf

!  do while(allcomms(comm_ens_index(current_ens_layer))%parent_index > 0 ) 
!     current_comm_index = allcomms(comm_ens_index(current_ens_layer))%parent_index
!  enddo
 ! do while(allcomms(current_comm_index)%parent_index > 0 ) 
 !    current_comm_index = allcomms(current_comm_index)%parent_index
 ! enddo
 ! comm_charmm=allcomms(current_comm_index)%comm
 ! mynod      =allcomms(current_comm_index)%mynod
 ! numnod     =allcomms(current_comm_index)%numnod
  comm_charmm=old_comm_charmm
  mynod = old_mynod
  numnod = old_numnod
  if(mynod > 0) close(5)
  call mpi_barrier(comm_charmm,ierror)

  !------------ Close open charmm.out files --------------------------
  

  !------------ Shut down communicators -------------------------------
  !--- comm_ensemble will be the one to shut down and all its children


!  !------------ All done, charmm_comm is reset, if verbose, print message ----
!  if(repd_ensemble_verbose)then
!     write(tbuf, '(a,i5,a)') &
!          ' ENSEMBLE>    REPLICA NODE ', WHOIAM, ' STOPPING!'
!     call ensprn(outu,tbuf,len_trim(tbuf))
!  endif
  return 
end subroutine repd2_ensfin

!----------------------------------------------------------------------
!          REPD_SYNC
!----------------------------------------------------------------------
subroutine repd_sync(comlyn,comlen)
  !-----------------------------------------------------------------------
  ! Set the number of replicas NENSEM
  !-----------------------------------------------------------------------
  use chm_kinds
  use repd_ensemble
  use stream
  use dimens_fcm
  use string
  use memory
  use parallel
  use parallel_groups
  use mpi
  use comand,only:comlyn_save,comlen_save
  implicit none

  character(len=*) comlyn
  character(len=4) wrd
  integer comlen,ierror
  character(len=80) :: chint

  wrd=nexta4(comlyn,comlen)
  
  call repd_reps_barrier(ierror)

  return

end subroutine repd_sync

!----------------------------------------------------------------------
!            REPD_REPS_BARRIER
!----------------------------------------------------------------------
subroutine repd_reps_barrier(ierror)
  use mpi
  use repd_ensemble
  use parallel,only:comm_charmm,mynod,mynodg
  implicit none
  integer,intent(out) :: ierror
  integer :: ier
  character(len=80) :: chint

  !----- First sync all the slaves of this master ------------
  !  write(chint,'(2i4)')comm_charmm,mynod
  !  call rensprint("   ENS_SYNC>> calling rep barrier with communicator",chint)
  call mpi_barrier(comm_charmm,ierror)
  !  call rensprint("   ENS_SYNC>>  "," after rep barrier")

  !----- Now sync all the masters of the reps ----------------
  if(lmasternode) then
     call mpi_barrier(comm_master,ierror)
     if(ierror/=0)print *,mynod,mynodg,"REPD_REPS_BARRIER ERROR"
  endif
  !  call rensprint("   ENS_SYNC>>  "," after master barrier")
  ierror = ier + ierror
  return
end subroutine repd_reps_barrier

!----------------------------------------------------------------------
!          REPD_BCAST_MASTERS
!----------------------------------------------------------------------
subroutine REPD_bcast_masters(array, length)
  !-----------------------------------------------------------------------
  ! This routine performs on node 0 broadcast to all other MASTER nodes
  ! and receive from node 0 on all other nodes.
  ! Usually called after read on node 0. For single precision arrays.

  use chm_kinds
  use dimens_fcm
  use repd_ensemble
  use mpi

  implicit none
  integer,intent(in) :: length
  integer,dimension(length),intent(inout) :: array(*)

  integer status

  if(nensem == 1) return
  if (.not.lmasternode) return
  call mpi_bcast(array,4*length,mpi_byte,0,comm_master,status)
  return
end subroutine REPD_bcast_masters

!----------------------------------------------------------------------
!          REPD_BCAST_ALL
!----------------------------------------------------------------------
subroutine repd_bcast_all(array, length)
  !-----------------------------------------------------------------------
  ! This routine performs on node 0 broadcast to all other MASTER nodes
  ! Then broadcast to all nodes in each copy.
  use chm_kinds
  use dimens_fcm
  use repd_ensemble
  use mpi    
  use parallel,only:comm_charmm
  implicit none
  integer,intent(in) :: length
  integer,intent(inout),dimension(length) :: array
  
  integer status
  
  if(nensem == 1) return
  if(lmasternode) then
     call mpi_bcast(array,4*length,mpi_byte,0,comm_master,status)
  endif
  call mpi_bcast(array,4*length,mpi_byte,0,comm_charmm,status)
  return
end subroutine repd_bcast_all

!----------------------------------------------------------------------
!          PSND8_repd
!----------------------------------------------------------------------
subroutine psnd8_repd(array, length)
  !-----------------------------------------------------------------------
  ! This routine performs on node 0 broadcast to all other nodes
  ! and receive from node 0 on all other nodes.
  ! Usually called after read on node 0. For real(chm_real) precision arrays.

  use chm_kinds
  use dimens_fcm
  use repd_ensemble
  use mpi
  implicit none
  real(chm_real) array(*)
  integer length

  integer status

  if(nensem == 1) return
  if (lslavenode) return
  call mpi_bcast(array,8*length,mpi_byte,0,comm_master,status)
  return
end subroutine psnd8_repd



!----------------------------------------------------------------------
!          ENS_SETUP_STREAM_INPUT
!----------------------------------------------------------------------
subroutine repd2_setup_stream_input()
  use chm_kinds
  use dimens_fcm,only:mxcmsz
  use stream
  use string,only:cnvtuc
  use comand,only:comlyn_save,comlen_save
  use repd_ensemble,only:lmasternode,Ensmasternod,comm_master
  use mpi
#ifdef __INTEL_COMPILER
  use ifport, only: ftell
#endif

  implicit none
  character(len=mxcmsz) :: line,linesav
  character(len=80) :: chint
  integer :: linelen,offset_cinfile
  logical keepreading,eof
  integer ierror

#ifdef __PGI
  interface
     integer function ftell(lu) bind(C)
       integer :: lu
     end function ftell
  end interface
#endif /* PGI */

  
  !--- Reset input file reading -- Open file for each master and read through
  !---     the ensemble command
  if(ensmasternod == 0)then
     if(.not. qcmdin) &
          call wrndie("ENS_setup_stream<ensemble.src>", &
          "Input must NOT be stdin, use -i or -input")
     offset_cinfile=ftell(5)
     print *,"MASTER has offset ",offset_cinfile
  endif
  close(5)
  if(lmasternode) then
     call mpi_bcast(offset_cinfile,1,mpi_integer,0,comm_master,ierror)
     call mpi_bcast(lcinfil,1,mpi_integer,0,comm_master,ierror)
     call mpi_bcast(cinfile,lcinfil,mpi_character,0,comm_master,ierror)
  
     write(chint,'(i4)')ensmasternod
     !call ensprint("   ENS_setup_stream>> ensmasternod ",chint )
     if(lmasternode)then     !(ensmasternod>0)then
        print *,ensmasternod," opening input file ",cinfile(1:lcinfil)
        open(unit=5,file=cinfile(1:lcinfil))
        print *,ensmasternod,"  has seek point of ",offset_cinfile
        !mfc     rewind(5)
        !mfc     call ensprint("   ENS_setup_stream>> opening input ",cinfile(1:lcinfil) )
        !mfc     write(chint,'(i6)')istrm
        !mfc     call ensprint("   ENS_setup_stream>> istrm ",chint )
        !mfc     line(1:mxcmsz) = " "
        !mfc     linesav(1:comlen_save) = comlyn_save(1:comlen_save)
        !mfc     call cnvtuc(linesav,comlen_save)
        !mfc     keepreading = .true.
        !mfc     do while(keepreading)
        !mfc        eof=.false.
        !mfc!MFC        read(5,'(a)')line
        !mfc        !-- setting qcmprs false so it will not broadcast ---
        !mfc        call rdcmnd(line,mxcmsz,linelen,istrm,eof,.false.,.false., &
        !mfc             'ENS_setup_stream read> ')
        !mfc        write(chint,*)eof
        !mfc        call ensprint("   ENS_setup_stream>> eof ",chint(1:len_trim(chint)) )
        !mfc!MFC        if(eof)call mpi_abort(mpi_comm_world,-1,status)
        !mfc        
        !mfc        call ensprint("   ENS_setup_stream>> line returned",line(1:len_trim(line)) )
        !mfc        call cnvtuc(line,len_trim(line))
        !mfc        !        keepreading = .not. ( indx(line,len_trim(line),"ENSE",4) > 0 )
        !mfc        keepreading = .not. ( line(1:comlen_save) == linesav(1:comlen_save) )
        !mfc     enddo
        !mfc          call ensprint("   ENS_setup_stream>> found ensemble line ",line(1:len_trim(line)))
     endif
#if KEY_IFORT==1
     ierror = fseek(5,offset_cinfile,0)
#else
     call fseek(5,offset_cinfile,0,ierror)
#endif
     !call ensprint("   ENS_setup_stream>> leaving "," ")
  endif


  return
end subroutine repd2_setup_stream_input



!MFC !----------------------------------------------------------------------
!MFC !            ENSINI
!MFC !----------------------------------------------------------------------
!MFC subroutine ensini
!MFC   !----------------------------------------------------------------------
!MFC   ! Initialize ensemble when charmm is started
!MFC   !----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use stream
!MFC   use ensemble
!MFC   use parallel,only:mynod,numnod,comm_charmm
!MFC   use mpi
!MFC   implicit none
!MFC   integer l
!MFC   integer ierror
!MFC   character(len=ensbfl) tbuf
!MFC 
!MFC   
!MFC   ensemble_verbose = .false.
!MFC   !ensemble_verbose = .true.
!MFC   slave_ens = .false.
!MFC   old_mynod=mynod
!MFC   ensemble_layers=0
!MFC   ! Set WHOIAM and NENSEM for now
!MFC   !   these will be modified when/if the user calls ENSEMBLE NENSEM command
!MFC   whoiam = mynod
!MFC   nensem = numnod
!MFC   lmasternode=.true.
!MFC   call mpi_barrier(comm_charmm, ierror)
!MFC 
!MFC   comm_master = comm_charmm
!MFC   if(ensemble_verbose)then
!MFC      write(tbuf, '(a,i5,a)') &
!MFC           ' ENSEMBLE>    REPLICA NODE ', WHOIAM, ' STARTED!'
!MFC      call ensprn(outu,tbuf,len_trim(tbuf))
!MFC 
!MFC      write(tbuf,'(3i5)')prnlev,wrnlev,iolev
!MFC      call ensprint("  ENSINI>>  plev,wlev,iolev: ",tbuf)
!MFC   endif
!MFC   if (nensem > maxens) then
!MFC      if (iolev > 0) write (outu,'(a)') &
!MFC           ' ENSINI> TOO MANY REPLICAS!', &
!MFC           ' ENSINI> INCREASE MAXENS IN ensemble.src'
!MFC      CALL WRNDIE(0, '<ENSINI>', 'TOO MANY REPLICAS REQUESTED')
!MFC   endif
!MFC 
!MFC   call set_param('WHOIAM',whoiam)
!MFC   call set_param('NENSEM',nensem)
!MFC   jrex = .false.
!MFC   jrswap = .false.
!MFC   qensexp = .false.
!MFC   ensswmethod = 0
!MFC   ensas = .false.
!MFC 
!MFC   !     this is ignored:
!MFC !MFC-- Huh?
!MFC   swapc = .true.
!MFC   t2repo = -1
!MFC   rep2to = -1
!MFC   ensexpu = -1
!MFC 
!MFC   if (iolev > 0 .and. ensemble_verbose) then
!MFC #if KEY_STRINGM==1
!MFC #if KEY_ABPO==1
!MFC      write(outu,'(a)')  &
!MFC           ' ENSINI> ALL REPLICAS PRESENTLY USE ONE TEMPERATURE', &
!MFC           ' ENSINI> USE "ENSEmble EXCHange"', &
!MFC           ' ENSINI> TO SET UP REPLICA EXCHANGE.', &
!MFC           ' ENSINI> USE ENSEmble STRIng [...]', &
!MFC           ' ENSINI> TO USE THE STRING METHOD.', &
!MFC           ' ENSINI> USE ENSEmble ABPO [...]', &
!MFC           ' ENSINI> TO USE THE ABPO METHOD.'
!MFC #else /**/
!MFC      write(outu,'(a)')  &
!MFC           ' ENSINI> ALL REPLICAS PRESENTLY USE ONE TEMPERATURE', &
!MFC           ' ENSINI> USE "ENSEmble EXCHange"', &
!MFC           ' ENSINI> TO SET UP REPLICA EXCHANGE.', &
!MFC           ' ENSINI> USE ENSEmble STRIng [...]', &
!MFC           ' ENSINI> TO USE THE STRING METHOD.'
!MFC #endif 
!MFC #else /**/
!MFC #if KEY_ABPO==1
!MFC      write(outu,'(a)')  &
!MFC           ' ENSINI> ALL REPLICAS PRESENTLY USE ONE TEMPERATURE', &
!MFC           ' ENSINI> USE "ENSEmble EXCHange"', &
!MFC           ' ENSINI> TO SET UP REPLICA EXCHANGE.', &
!MFC           ' ENSINI> USE ENSEmble ABPO [...]', &
!MFC           ' ENSINI> TO USE THE ABPO METHOD.'
!MFC #else /**/
!MFC      write(outu,'(a)')  &
!MFC           ' ENSINI> ALL REPLICAS PRESENTLY USE ONE TEMPERATURE', &
!MFC           ' ENSINI> USE "ENSEmble EXCHange"', &
!MFC           ' ENSINI> TO SET UP REPLICA EXCHANGE.'
!MFC #endif 
!MFC #endif 
!MFC   endif
!MFC   return
!MFC end subroutine ensini
!MFC 
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !            ENSPRN
!MFC !----------------------------------------------------------------------
!MFC subroutine ensprn(unum,message,l)
!MFC   !-----------------------------------------------------------------------
!MFC   ! Print a message for each process, _via_ root. 
!MFC   !     UNUM: unit to print to
!MFC   !     MESSAGE: message to write for this node
!MFC   !     L: length of message
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use stream
!MFC   use ensemble
!MFC   use parallel, only:mynod,comm_charmm
!MFC   use mpi
!MFC   implicit none
!MFC 
!MFC   character(len=*) message
!MFC   integer,intent(in) :: unum, l
!MFC   ! local
!MFC   integer i
!MFC   integer, parameter :: biglen=ensbfl*maxens
!MFC   character(len=biglen) msgarray
!MFC   integer ierror,lenst,stat(mpi_status_size)
!MFC 
!MFC   if(l == 0) return
!MFC   if(lmasternode)then
!MFC      if(nensem == 1) then
!MFC         write(unum, '(a)') message
!MFC      else
!MFC         lenst = l   !len_trim(message(1:1))
!MFC      endif
!MFC      if (whoiam == 0 ) write(unum, '(a)') message(1:lenst)
!MFC      do i=1,nensem-1
!MFC         if(whoiam == 0) then
!MFC            call mpi_recv(message,lenst,mpi_byte,i,i, &
!MFC                 comm_master,stat,ierror)
!MFC            write(unum, '(a)') message(1:lenst)
!MFC            
!MFC         else
!MFC            call mpi_send(message,lenst,mpi_byte,0,i, &
!MFC                 comm_master,ierror)
!MFC         endif
!MFC      enddo
!MFC   endif
!MFC   call gflush(unum)
!MFC   return
!MFC end subroutine ensprn
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !            ENSCMD
!MFC !----------------------------------------------------------------------
!MFC subroutine enscmd(comlyn,comlen)
!MFC   !----------------------------------------------------------------------
!MFC   !     parse command line and call subcommands, namely
!MFC   !           OPEN
!MFC   !           CLOSE
!MFC   !           SEED
!MFC   !     etc...
!MFC   !----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use string
!MFC   use stream
!MFC   use dimens_fcm
!MFC   use coord
!MFC   use parallel,only:mynod
!MFC #if KEY_ABPO==1
!MFC   use abpo, only:abpo_setup  
!MFC #endif
!MFC   implicit none
!MFC 
!MFC   character(len=*) comlyn
!MFC   integer comlen,ierror
!MFC 
!MFC   if (indxa(comlyn,comlen,'OPEN')  >  0) then
!MFC      call ensopn(comlyn,comlen)
!MFC   else if (indxa(comlyn,comlen,'CLOSE')  >  0) then
!MFC      call ensclo(comlyn,comlen)
!MFC   else if (indxa(comlyn,comlen,'SEED')  >  0) then
!MFC      call enseed(comlyn,comlen)
!MFC   else if (indxa(comlyn,comlen,'PRSEED')  >  0) then
!MFC      call ensdpr(comlyn,comlen)
!MFC   else if (indxa(comlyn,comlen,'EXPAVG')  >  0) then
!MFC      call eavgsetup(comlyn,comlen)
!MFC   else if (indxa(comlyn,comlen,'EXCH')  >  0) then
!MFC      call ensrex(comlyn,comlen)
!MFC   else if (indxa(comlyn,comlen,'WRIT')  >  0) then
!MFC      call enswri(comlyn,comlen)
!MFC #if KEY_STRINGM==1
!MFC   else if (indxa(comlyn,comlen,'STRI')  >  0) then
!MFC      CALL parse_string_commands(comlyn,comlen) ! parse string commands
!MFC #endif 
!MFC #if KEY_ABPO==1
!MFC   else if (indxa(comlyn,comlen,'ABPO')  >  0) then
!MFC      call abpo_setup(comlyn,comlen) ! parse abpo commands
!MFC #endif 
!MFC   else if (indxa(comlyn,comlen,'SWON')  >  0) then
!MFC      jrswap=.true.
!MFC      if (iolev > 0)  &
!MFC           write(outu,'(a)') &
!MFC           'SWAPPING OF TEMPERATURE REPLICAS ENABLED'
!MFC   else if (indxa(comlyn,comlen,'SWOFF')  >  0) then
!MFC      jrswap=.false.
!MFC      if (iolev > 0)  &
!MFC           write(outu,'(a)')  &
!MFC           'SWAPPING OF TEMPERATURE REPLICAS DISABLED'
!MFC   else if (indxa(comlyn,comlen,'INFO')  >  0) then
!MFC      call ensinf()
!MFC 
!MFC   else if (indxa(comlyn,comlen,'NENSEM')  >  0) then
!MFC !     call ensini
!MFC      call ensnensem(comlyn,comlen)
!MFC 
!MFC   elseif (indxa(comlyn,comlen,'SYNC')  >  0) then
!MFC      call ens_sync(comlyn,comlen)
!MFC 
!MFC   else
!MFC      call wrndie(0, '<ENSPAR>', 'UNRECOGNIZED SUBCOMMAND')
!MFC   endif
!MFC   call ensprint("  ENSCMD>> returning "," ")
!MFC   return
!MFC end subroutine enscmd
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSDR
!MFC !----------------------------------------------------------------------
!MFC subroutine ensdpr(comlyn,comlen)
!MFC   !----------------------------------------------------------------------
!MFC   !     PRINT RANDOM SEED FOR EACH NODE
!MFC   !----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use dimens_fcm
!MFC   use consta
!MFC   use stream
!MFC   use string
!MFC   use number
!MFC   use reawri
!MFC   use parallel, only: mynod
!MFC   implicit none
!MFC 
!MFC   character(len=*) comlyn
!MFC   integer comlen
!MFC   ! local ...
!MFC   integer rseed,i,j
!MFC   integer seeds(maxens),myseed
!MFC   character(len=ensbfl) tbuf
!MFC   real(chm_real) r
!MFC 
!MFC   if (mynod == 0) then
!MFC      rseed=gtrmi(comlyn,comlen,'ROOT',-1)
!MFC      if (iolev > 0) then
!MFC         write(outu,'(a)') ' ENSEED> RANDOM SEEDS FROM EACH NODE: '
!MFC      end if
!MFC      write(tbuf,'(a,i3,a,i20)') ' ENSEMBLE> REPLICA NODE ', WHOIAM, ' HAS ISEED = ', ISEED
!MFC      call ensprn(outu,tbuf,len(tbuf))
!MFC   endif
!MFC 
!MFC   return
!MFC end subroutine ensdpr
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSEED
!MFC !----------------------------------------------------------------------
!MFC subroutine enseed(comlyn,comlen)
!MFC   !----------------------------------------------------------------------
!MFC   !     The root node takes the random seed given on the command line
!MFC   !     and generates random seeds for all the other nodes, which are
!MFC   !     placed in the ?myseed variable on each node
!MFC   !     WARNING: THIS IS NOT A GOOD IDEA!
!MFC   !----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use dimens_fcm
!MFC   use consta
!MFC   use stream
!MFC   use string
!MFC   use number
!MFC   use clcg_mod
!MFC   use parallel, only: mynod,comm_charmm
!MFC   implicit none
!MFC 
!MFC   character(len=*) comlyn
!MFC   integer comlen
!MFC   ! local ...
!MFC   integer rseed,i,j
!MFC   integer seeds(maxens),myseed
!MFC   character(len=ensbfl) tbuf
!MFC   real(chm_real) r
!MFC 
!MFC   rseed=gtrmi(comlyn,comlen,'ROOT',-1)
!MFC   if (iolev > 0) then
!MFC      WRITE(OUTU,'(A)') ' ENSEED> ASSIGNING SEEDS AUTOMATICALLY'
!MFC      WRITE(OUTU,'(A)') ' ENSEED> THIS IS A BAD IDEA'
!MFC      WRITE(OUTU,'(A)') ' ENSEED> RATHER SUPPLY YOUR OWN SEEDS!'
!MFC   end if
!MFC   if (whoiam == 0 .and. mynod == 0) then
!MFC      seeds(1) = rseed
!MFC      do i=2,nensem
!MFC         r = random(rseed)
!MFC         seeds(i) = rseed
!MFC      end do
!MFC   endif
!MFC 
!MFC   call ens_bcast_all(seeds,nensem)
!MFC   myseed = seeds(whoiam+1)
!MFC   call set_param('MYSEED',myseed)
!MFC   write(tbuf,'(a,i3,a,i20)') &
!MFC        ' ENSEMBLE> REPLICA NODE ', WHOIAM, ' HAS SEED ', MYSEED
!MFC   call ensprn(outu,tbuf,len(tbuf))
!MFC   return
!MFC end subroutine enseed
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSINQ
!MFC !----------------------------------------------------------------------
!MFC subroutine ensinq(mode,name,maxlen,length,qopen,qform,qwrite, &
!MFC      qens,unit)
!MFC   !----------------------------------------------------------------------
!MFC   !     taken from VINQRE()
!MFC   !
!MFC   !     file inquiry by file name or FORTRAN unit
!MFC   !     Flag QOPEN indicates whether file or unit is "open".
!MFC   !     Flag QFORM indicates whether file was opened formatted.
!MFC   !     Flag QWRITe indicates whether file was opened write-access.
!MFC   !     For inquiry by unit MAXLEN has to be specified (max length of NAME
!MFC   !     and LENGTH returns with the length of NAME.
!MFC   !     For inquiry by file the two names INPUT and OUTPUT are reserved
!MFC   !     for the standard input and output channels 5 and OUTU.
!MFC   !
!MFC   use chm_kinds
!MFC   use stream
!MFC   use dimens_fcm
!MFC   use machio
!MFC   use string
!MFC   use parallel
!MFC   use ensemble,only:comm_master,lmasternode
!MFC   implicit none
!MFC   character(len=mxcmsz) ename
!MFC   integer elen
!MFC   character(len=*) mode, name
!MFC   integer maxlen, length
!MFC   logical qopen, qform, qwrite, qens
!MFC   integer unit
!MFC   integer qshare(4),i
!MFC 
!MFC   qopen=.true.
!MFC   length=0
!MFC 
!MFC   if (mode == 'FILE'.and.name == 'INPUT') then
!MFC      unit=5
!MFC   else if (mode == 'FILE' .and. name == 'OUTPUT') then
!MFC      unit=outu
!MFC   else if (mode == 'FILE') then
!MFC      ename = name
!MFC      call expnam(ename,elen,qopen)
!MFC      qopen = .not.qopen
!MFC      unit = -1
!MFC      if(qopen) then
!MFC         inquire(file=ename(1:elen),opened=qopen,number=unit)
!MFC         if (unit <= 0.or.unit > 99) qopen=.false.
!MFC      endif
!MFC   else if (mode == 'UNIT') then
!MFC      inquire(unit=unit,opened=qopen,name=name)
!MFC      length=maxlen
!MFC      call trime(name,length)
!MFC 
!MFC   endif
!MFC 
!MFC   !     if file is open then get QFORM and QWRITE flags
!MFC   !     ... and QENS
!MFC   if (qopen.and.iolev > 0) then
!MFC      if ((ifreeu(unit) == +1).or.(ifreeu(unit) == +7)) then
!MFC         qform=.true.
!MFC         qwrite=.true.
!MFC      else if ((ifreeu(unit) == +10).or.(ifreeu(unit) == +70)) then
!MFC         qform=.true.
!MFC         qwrite=.false.
!MFC      else if ((ifreeu(unit) == -1).or.(ifreeu(unit) == -7)) then
!MFC         qform=.false.
!MFC         qwrite=.true.
!MFC      else if ((ifreeu(unit) == -10).or.(ifreeu(unit) == -70)) then
!MFC         qform=.false.
!MFC         qwrite=.false.
!MFC      else if (ifreeu(unit) == 0) then
!MFC         qform=.false.
!MFC         qwrite=.false.
!MFC      endif
!MFC      if ((mod(ifreeu(unit),7) == 0).and.(ifreeu(unit) /= 0)) then
!MFC         qens = .true.
!MFC      else
!MFC         qens = .false.
!MFC      endif
!MFC   endif
!MFC   !     if single rather than parallel input only root knows file
!MFC   !     status
!MFC   qshare(1:4) = 0
!MFC 
!MFC   if (iolev > 0) then
!MFC      if (qopen) qshare(1)  = 1
!MFC      if (qform) qshare(2)  = 1
!MFC      if (qwrite) qshare(3) = 1
!MFC      if (qens) qshare(4)  = 1
!MFC   endif
!MFC !MFC-- Need to check this one out, probably should not be world.
!MFC !  if(lmasternode)call psnd4_comm(comm_master,qshare,4)
!MFC !  call psnd4_comm(comm_charmm,qshare,4)
!MFC !HH: the following line causes error when ensinq is not called on all
!MFC !  processors, commented out by HH
!MFC !  call ens_bcast_all(qshare,4)
!MFC   qopen=.false.
!MFC   qform=.false.
!MFC   qwrite=.false.
!MFC   qens=.false.
!MFC   if (qshare(1) == 1) qopen = .true.
!MFC   if (qshare(2) == 1) qform = .true.
!MFC   if (qshare(3) == 1) qwrite = .true.
!MFC   if (qshare(4) == 1) qens = .true.
!MFC   return
!MFC end subroutine ensinq
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSDMP
!MFC !----------------------------------------------------------------------
!MFC subroutine ensdmp()
!MFC   !-----------------------------------------------------------------------
!MFC   !     dumps coordinate files to "node0.pdb" "node1.pdb" etc on a 
!MFC   !     crash
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use dimens_fcm
!MFC   use stream
!MFC   use string
!MFC   use coord
!MFC   use psf
!MFC   use ctitla
!MFC   use machio
!MFC   use memory
!MFC   use coorio_mod,only:cwrite
!MFC   use parallel, only: mynod
!MFC 
!MFC   implicit none
!MFC   integer,allocatable,dimension(:) :: sel
!MFC   integer, parameter :: mxfile=128
!MFC   character(len=mxfile) filenm, junknm, tmp
!MFC   integer crashu,flen
!MFC   integer icntrl(20)
!MFC   integer i
!MFC   logical qopen,qform,qwrite,qens
!MFC 
!MFC   if (mynod == 0) then
!MFC      crashu = 50
!MFC      write(filenm, '(a,i2,a)') &
!MFC           'node',whoiam,'.pdb'
!MFC      if(filenm(5:5) == ' ') filenm(5:5)='0'
!MFC      flen=strlng(filenm)
!MFC      call ensinq('UNIT',junknm,mxfile,flen,qopen,qform,qwrite,qens, &
!MFC           crashu)
!MFC      if (qopen) then
!MFC         if(wrnlev >= 2.and.iolev > 0) write(outu,'(2a)') &
!MFC              ' ENSDMP> Unit already open.', &
!MFC              ' The old file will be closed first.'
!MFC         close(unit=crashu)   
!MFC         ifreeu(crashu) = 0
!MFC      endif
!MFC      OPEN(UNIT=CRASHU,FILE=FILENM,FORM='FORMATTED',STATUS='UNKNOWN', &
!MFC           ACCESS='SEQUENTIAL')
!MFC      call chmalloc('ensemble.src','ensdmp','sel',natom+1,intg=sel)
!MFC      sel(1:natom) = 1
!MFC      call cwrite(crashu,titlea,ntitla,icntrl,x,y,z,wmain, &
!MFC           res,atype,ibase,nres,natom,sel,4,0,0,.false.)
!MFC      
!MFC      close(unit=crashu)
!MFC      call chmdealloc('ensemble.src','ensdmp','sel',natom+1,intg=sel)
!MFC   endif
!MFC 
!MFC   return 
!MFC end subroutine ensdmp
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSOPN
!MFC !----------------------------------------------------------------------
!MFC subroutine ensopn(comlyn,comlen)
!MFC   !-----------------------------------------------------------------------
!MFC   ! Opens 1 file per replica process
!MFC   ! TODO: Include checking to make sure each replica's file name is
!MFC   ! different (when opening with write access).
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use dimens_fcm
!MFC   use consta
!MFC   use stream
!MFC   use string
!MFC   use number
!MFC   use machio
!MFC   use parallel,only:mynod
!MFC   implicit none
!MFC 
!MFC   character(len=*) comlyn
!MFC   integer comlen
!MFC 
!MFC   integer mxpref, mxsuff
!MFC   integer flen, l, unum, unum2, ierror, i
!MFC   logical err, qopen, qform, qwrite, qens
!MFC #if KEY_UNIX==1
!MFC   integer, parameter :: mxfile=128
!MFC #else /**/
!MFC   integer, parameter :: mxfile=256
!MFC #endif 
!MFC   integer ipt
!MFC   integer, parameter :: fmlen=30,stslen=20,faccl=20
!MFC   character(len=mxfile) filex, junknm
!MFC   character(len=fmlen) formt
!MFC   character(len=stslen) fstat
!MFC   character(len=faccl) facc
!MFC   character(len=ensbfl) tbuf
!MFC 
!MFC  if(iolev < 0 )then
!MFC      call ensprint("   ENSOPN>> "," not opening file")
!MFC      comlyn=""
!MFC      comlen=0
!MFC 
!MFC   else
!MFC      unum=gtrmi(comlyn,comlen,'UNIT',-1)
!MFC      if (unum < 0) then
!MFC         call wrndie(0,'<OPNLGU>','NO UNIT NUMBER SPECIFIED')
!MFC         return
!MFC      endif
!MFC      
!MFC      ! filename
!MFC      call gtrmwa(comlyn, comlen, 'NAME', 4, filex, mxfile, flen)
!MFC      if (flen <= 0) then
!MFC         call wrndie(0, '<ENSOPN>', 'NO FILE NAME GIVEN')
!MFC      endif
!MFC      call ensprint("  ENSOPN>> " ,filex)
!MFC 
!MFC      ! remove quotes
!MFC      ipt=0
!MFC      do i=1,flen
!MFC         if(filex(i:i) /= '"') then
!MFC            ipt=ipt+1
!MFC            filex(ipt:ipt)=filex(i:i)
!MFC         endif
!MFC      enddo
!MFC      filex = filex(1:ipt)
!MFC      flen=ipt
!MFC      
!MFC      !     convert file name to lowercase
!MFC      if(lower) call cnvtlc(filex,flen)
!MFC      call ensprint("   ENSOPN>> filename",filex)
!MFC      
!MFC      !     unit already open?
!MFC      call ensinq('UNIT',junknm,mxfile,flen,qopen,qform,qwrite,qens, &
!MFC           unum)
!MFC      if (qopen) then
!MFC         if(wrnlev >= 2.and.iolev > 0) write(outu,'(2a)') &
!MFC              ' ENSOPN> Unit already open.', &
!MFC              ' The old file will be closed first.'
!MFC         close(unit=unum)   
!MFC         ifreeu(unum) = 0
!MFC      endif
!MFC      
!MFC      !     file already in use
!MFC      call ensinq('FILE',filex,i,i,qopen,qform,qwrite,qens,unum2)
!MFC      if (qopen) then
!MFC         if(wrnlev >= 2.and.iolev > 0) write(outu,'(a,/,2a)') &
!MFC              ' ENSOPN> ***** WARNING ***** another unit is already ', &
!MFC              '         assigned to the file -', &
!MFC              ' it will be disconnected first.'
!MFC         close(unit=unum2)   
!MFC         ifreeu(unum2) = 0
!MFC      endif
!MFC      
!MFC      ! format
!MFC      if (indxa(comlyn,comlen,'UNFO') > 0) then
!MFC         formt='UNFORMATTED'
!MFC      else if (indxa(comlyn,comlen,'FILE') > 0) then
!MFC         FORMT='UNFORMATTED'
!MFC      ELSE IF (INDXA(COMLYN,COMLEN,'FORM') > 0) THEN
!MFC         FORMT='FORMATTED'
!MFC      ELSE IF (INDXA(COMLYN,COMLEN,'CARD') > 0) THEN
!MFC         FORMT='FORMATTED'
!MFC      else
!MFC         call wrndie(1,'<OPNLGU>', &
!MFC              'NO FORMATTING SPECIFICATION, IT WILL BE OPENED UNFORMATTED')
!MFC         FORMT='UNFORMATTED'
!MFC      endif
!MFC      
!MFC      ! status
!MFC      call gtrmwa(comlyn, comlen, 'STAT', 4, fstat, stslen, l)
!MFC      if (l <= 0) then
!MFC         fstat = 'UNKNOWN'
!MFC      endif
!MFC      
!MFC      ! access
!MFC      facc='READ'
!MFC      if (indxa(comlyn,comlen,'APPE') > 0) then
!MFC         facc='APPEND'
!MFC      ELSE IF (INDXA(COMLYN,COMLEN,'READ') > 0) THEN
!MFC         FACC='READ'
!MFC      ELSE IF (INDXA(COMLYN,COMLEN,'WRIT') > 0) THEN
!MFC         FACC='WRITE'
!MFC      ELSE
!MFC         FACC='READ'
!MFC      endif
!MFC 
!MFC      tbuf(1:len(tbuf))=" "
!MFC      write (tbuf, '(a,I3)') ' ENSEMBLE>   REPLICA NODE ', WHOIAM
!MFC      call ensprint(tbuf(1:len_trim(tbuf))," ") 
!MFC      tbuf(1:len(tbuf))=" "
!MFC      write (tbuf, '(a26,a54)') ' ENSEMBLE>   OPENING FILE ', FILEX
!MFC      call ensprint(tbuf(1:len_trim(tbuf))," ") 
!MFC      tbuf(1:len(tbuf))=" "
!MFC      WRITE (tbuf, '(A,I3)')  ' ENSEMBLE>   ON UNIT ', UNUM
!MFC      call ensprint(tbuf(1:len_trim(tbuf))," ") 
!MFC      tbuf(1:len(tbuf))=" "
!MFC      WRITE (tbuf, '(A,A15,A,A10)')  ' ENSEMBLE>   WITH FORMAT ',  &
!MFC           FORMT, ' AND ACCESS ', FACC
!MFC      call ensprint(tbuf(1:len_trim(tbuf))," ") 
!MFC 
!MFC      ! open it
!MFC      IF (FACC == 'APPEND') THEN
!MFC         OPEN(UNIT=UNUM,FILE=FILEX,FORM=FORMT,STATUS='OLD', &
!MFC              ACCESS='SEQUENTIAL')
!MFC      ELSE IF (FACC == 'READ') THEN
!MFC         OPEN(UNIT=UNUM,FILE=FILEX,FORM=FORMT,STATUS='OLD', &
!MFC              ACCESS='SEQUENTIAL')
!MFC      ELSE IF (FACC == 'WRITE') THEN
!MFC         OPEN(UNIT=UNUM,FILE=FILEX,FORM=FORMT,STATUS='UNKNOWN', &
!MFC              ACCESS='SEQUENTIAL')
!MFC      END IF
!MFC 
!MFC      ! update ifreeu array
!MFC      INQUIRE(FILE=FILEX,OPENED=QOPEN,NUMBER=UNUM)
!MFC      IF (.NOT. QOPEN) THEN
!MFC         CALL WRNDIE(0, '<ENSOPN>', 'Could not open file')
!MFC      ELSE
!MFC         !
!MFC         !     put appropriate code in IFREEU array to play nicely
!MFC         !     with rest of charmm:
!MFC         !        +70 ensemble read formatted
!MFC         !        +10 read formatted
!MFC         !        +7  ensemble write/append formatted
!MFC         !        +1  write/append formatted
!MFC         !        -1  write/append unformatted
!MFC         !        -7  ensemble write/append unformatted
!MFC         !        -10 read unformatted
!MFC         !        -70 ensemble read unformatted
!MFC         !       i.e. ifreeu(unum)%7 tells whether we have ensemble file
!MFC         IF (FORMT == 'FORMATTED') THEN
!MFC            IFREEU(UNUM)=7
!MFC         ELSE
!MFC            IFREEU(UNUM)=-7
!MFC         ENDIF
!MFC         IF (FACC == 'READ') IFREEU(UNUM)=IFREEU(UNUM)*10
!MFC      endif
!MFC   endif
!MFC   call ensprint("   ENSOPN>> returning "," ")
!MFC 
!MFC   return
!MFC end subroutine ensopn
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSCLO
!MFC !----------------------------------------------------------------------
!MFC subroutine ensclo(comlyn,comlen)
!MFC   !-----------------------------------------------------------------------
!MFC   ! Closes replica output
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use stream
!MFC   use dimens_fcm
!MFC   use string
!MFC   use machio
!MFC   use mpi
!MFC   use parallel,only:comm_charmm
!MFC   implicit none
!MFC 
!MFC   character(len=*) comlyn
!MFC   integer comlen
!MFC   integer unum, ierror
!MFC   character(len=ensbfl) tbuf
!MFC 
!MFC   ! words
!MFC   if (iolev > 0) then
!MFC      write(outu,'(a80)') comlyn
!MFC   endif
!MFC   unum=gtrmi(comlyn,comlen,'UNIT',-1)
!MFC   write (tbuf, '(a,i3,a,i3)') ' ENSEMBLE>   REPLICA NODE ',  &
!MFC        whoiam,  &
!MFC        ' CLOSING UNIT ', UNUM
!MFC   call ensprn(outu,tbuf,len(tbuf))
!MFC 
!MFC   ! action
!MFC   close(unit=unum)   
!MFC   ifreeu(unum) = 0
!MFC   call mpi_barrier(comm_charmm, ierror)
!MFC 
!MFC   return
!MFC end subroutine ensclo
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSREX
!MFC !----------------------------------------------------------------------
!MFC subroutine ensrex(comlyn,comlen)
!MFC   !-----------------------------------------------------------------------
!MFC   ! Set up replica exchange -parses CHARMM command options for
!MFC   ! "ENSEMBLE EXCHANGE ..."
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use stream
!MFC   use dimens_fcm
!MFC   use string
!MFC   use consta
!MFC   use number,only: zero
!MFC   use memory
!MFC   use psf
!MFC   use parallel
!MFC 
!MFC   implicit none
!MFC   !
!MFC   character(len=*) comlyn
!MFC   integer comlen
!MFC   !
!MFC   integer i,j,mapu,ii,swrules,tmpi
!MFC   logical jautot
!MFC   real(chm_real) lowt,tgrad,pswap,tmpt,a,b,c,deff
!MFC 
!MFC   ! words
!MFC   deff=zero
!MFC   jrex=.true.
!MFC   jrswap=.true.
!MFC   jensc=.true.
!MFC   ensas=.false.
!MFC   t2repo=gtrmi(comlyn,comlen,'T2RE',-1)
!MFC   rep2to=gtrmi(comlyn,comlen,'REP2',-1)
!MFC   ensfrq=gtrmi(comlyn,comlen,'FREQ',0)
!MFC   ensswmethod=gtrmi(comlyn,comlen,'SWME',0)
!MFC   call ens_bcast_all(ensswmethod,1)
!MFC 
!MFC   mapu=gtrmi(comlyn,comlen,'MAPU',-1)
!MFC   swrules=gtrmi(comlyn,comlen,'RULE',-1)
!MFC   if (iolev > 0) then
!MFC      write(outu,'(a,i5,a)') ' ENSREX> RUNNING', &
!MFC           nensem, ' REPLICA(S)'
!MFC      write(outu,'(a)') ' ENSREX> SWITCHING REPLICAS EVERY '
!MFC      write(outu,'(a,i8,a)') ' ENSREX> ', ENSFRQ, ' DYNAMICS STEPS'
!MFC   endif
!MFC   if (indxa(comlyn, comlen, 'NOSC')  >  0) jensc=.false.
!MFC   if (indxa(comlyn, comlen, 'ASVE')  >  0) ensas=.true.
!MFC   i = indxa(comlyn, comlen, 'SWAPC')
!MFC   if (iolev > 0) write(outu,'(a)') &
!MFC        ' ENSREX> WILL SWAP COORDINATES BETWEEN REPLICAS'
!MFC   call chmalloc('ensemble.src','ENSREX','ENSBUF',NATOM,crl=ENSBUF)
!MFC   !
!MFC   jautot = (INDXA(COMLYN, COMLEN, 'AUTO')  >  0)
!MFC   if (jautot) then
!MFC      LOWT=GTRMF(COMLYN,COMLEN,'LOWT',DEFF)
!MFC      TGRAD=GTRMF(COMLYN,COMLEN,'TGRAD',DEFF)
!MFC      PSWAP=GTRMF(COMLYN,COMLEN,'PSWAP',DEFF)
!MFC      !          ... set up temperature series automagically...  
!MFC      enstem(1)=lowt
!MFC      do i=2,nensem
!MFC         a=tgrad
!MFC         b=(-2.0*tgrad+log(pswap)*kboltz)*enstem(i-1)
!MFC         c=tgrad*enstem(i-1)**2
!MFC         enstem(i)=(-b+sqrt(b**2-4.0*a*c))/(2.0*a)
!MFC      enddo
!MFC   else
!MFC      do i=1,nensem
!MFC         tmpt = nextf(comlyn,comlen)
!MFC         enstem(i)=tmpt
!MFC      enddo
!MFC   endif
!MFC 
!MFC   if (mapu > 0) then
!MFC      if (iolev > 0) then
!MFC         do i=1,nensem
!MFC            read(mapu,*) ii,rep2t(i)
!MFC            t2rep(rep2t(i))=i
!MFC         enddo
!MFC      endif
!MFC      call ens_bcast_all(rep2t,nensem)
!MFC      call ens_bcast_all(t2rep,nensem)
!MFC   else 
!MFC      do i=1,nensem
!MFC         t2rep(i)=i
!MFC         rep2t(i)=i
!MFC      enddo
!MFC   endif
!MFC 
!MFC   if (swrules > 0) then
!MFC      !        read allowed swap rules from unit swrules
!MFC      if (iolev > 0) then
!MFC         read(swrules,*) ensnsw
!MFC         do i=1,ensnsw
!MFC            read(swrules,*) ensisw(i),ensjsw(i)
!MFC            if (ensisw(i) > ensjsw(i)) then
!MFC               tmpi=ensisw(i)
!MFC               ensisw(i)=ensjsw(i)
!MFC               ensjsw(i)=tmpi
!MFC            endif
!MFC         enddo
!MFC      endif
!MFC      call ens_bcast_all(ensnsw,1)
!MFC      call ens_bcast_all(ensisw,ensnsw)
!MFC      call ens_bcast_all(ensjsw,ensnsw)
!MFC   else 
!MFC      !        setup default swaps between 'adjacent' replicas
!MFC      ensnsw=nensem-1
!MFC      do i=1,ensnsw
!MFC         ensisw(i)=i
!MFC         ensjsw(i)=i+1
!MFC      enddo
!MFC   endif
!MFC 
!MFC   if (iolev > 0) then
!MFC      write(outu,'(a)')' ENSREX> THE FOLLOWING ARE ALLOWED EXCHANGES'
!MFC      do i=1,ensnsw
!MFC         write(outu,'(a,i4,a,i4)') ' ENSREX> ',ENSISW(I),' <---> ', &
!MFC              ENSJSW(I)
!MFC      enddo
!MFC   endif
!MFC 
!MFC   ensatt(1:ensnsw)=0
!MFC   enssuc(1:ensnsw)=0
!MFC 
!MFC   repswp(1:nensem)=0
!MFC   do i=1,nensem
!MFC      if (whoiam == (i-1)) then
!MFC         ensmyt=enstem(rep2t(i)) 
!MFC         call set_param('ensmyt',ensmyt)
!MFC      endif
!MFC   enddo
!MFC   ensmyt=enstem(rep2t(whoiam+1))  !MFC does this do the same thing as above loop? 
!MFC         call set_param('ENSMYT',ensmyt)
!MFC 
!MFC   return
!MFC end subroutine ensrex
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSWRI
!MFC !----------------------------------------------------------------------
!MFC subroutine enswri(comlyn,comlen)
!MFC   !     ... now redundant
!MFC   !-----------------------------------------------------------------------
!MFC   ! Write out replica target temperatures (for restart)
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use stream
!MFC   use dimens_fcm
!MFC   use string
!MFC   use parallel, only:mynod
!MFC   implicit none
!MFC 
!MFC   character(len=*) comlyn
!MFC   integer comlen
!MFC 
!MFC   integer unum,i
!MFC 
!MFC   if (mynod == 0) then
!MFC      UNUM=GTRMI(COMLYN,COMLEN,'UNIT',-1)
!MFC      if (iolev <= 0) return
!MFC      if (unum < 1)  call wrndie(-3,'<ENSWRI>', &
!MFC           'COULD NOT OPEN UNIT FOR OUTPUT.')
!MFC      do i=1,nensem
!MFC         write(unum,'(i5,i5)') i, rep2t(i)
!MFC      enddo
!MFC   endif
!MFC 
!MFC   return 
!MFC end subroutine enswri
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSCHK
!MFC !----------------------------------------------------------------------
!MFC subroutine enschk()
!MFC   !-----------------------------------------------------------------------
!MFC   ! Check that all replicas set up OK before starting dynamics
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use stream
!MFC   use dimens_fcm
!MFC   use parallel, only:mynod
!MFC   implicit none
!MFC   integer i
!MFC 
!MFC   if (mynod == 0) then
!MFC      if (.not.jrex) return
!MFC      do i=1,nensem
!MFC         if (enstem(i) < 1.0) then
!MFC            call wrndie(-3,'<enschk>', &
!MFC                 'TEMPERATURE NOT SET FOR SOME REPLICAS, QUITTING.')
!MFC         endif
!MFC      enddo
!MFC   endif
!MFC 
!MFC   return
!MFC end subroutine enschk
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSINF
!MFC !----------------------------------------------------------------------
!MFC subroutine ensinf()
!MFC   !-----------------------------------------------------------------------
!MFC   ! Print replica exchange information
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use stream
!MFC   use dimens_fcm
!MFC   use parallel, only:mynod
!MFC   implicit none
!MFC   !
!MFC   integer i
!MFC   integer unum, ierror
!MFC   character(len=ensbfl) tbuf
!MFC   !
!MFC   if (iolev > 0) then
!MFC      write(outu,'(a)') &
!MFC           ' ENSINF> **************************************************'
!MFC      write(outu,'(a)') &
!MFC           ' ENSINF> CURRENT REPLICA EXCHANGE INFO:'
!MFC      write(outu,'(a,i5)') &
!MFC           ' ENSINF>           NUMBER OF REPLICAS = ', NENSEM
!MFC      if (jrswap) then
!MFC         write(outu,'(a)') &
!MFC              ' ENSINF>   REPLICA SWAPPING ENABLED'
!MFC      else
!MFC         write(outu,'(a)') &
!MFC              ' ENSINF>   REPLICA SWAPPING DISABLED'
!MFC      endif
!MFC      if (jensc) then
!MFC         write(outu,'(a)') &
!MFC              ' ENSINF>   TEMPERATURE SCALING ENABLED'
!MFC      else
!MFC         write(outu,'(a)') &
!MFC              ' ENSINF>   TEMPERATURE SCALING DISABLED'
!MFC      endif
!MFC      if (ensas) then
!MFC         write(outu,'(a)') &
!MFC              ' ENSINF>   TEMPERATURE ASSIGNING ENABLED'
!MFC      else
!MFC         write(outu,'(a)') &
!MFC              ' ENSINF>   TEMPERATURE ASSIGNING DISABLED'
!MFC      endif
!MFC      write(outu,'(a,i8)') &
!MFC           ' ENSINF>   FREQUENCY FOR ATTEMPTING SWAPS = ', ENSFRQ
!MFC      write(outu,'(a,i3)') &
!MFC           ' ENSINF>     WRITING REPLICA MAP TO UNIT = ', T2REPO
!MFC      write(outu,'(a,i3)') &
!MFC           ' ENSINF> WRITING TEMPERATURE MAP TO UNIT = ', REP2TO
!MFC      do i=1,nensem
!MFC         write(outu,'(a,i5,a,f8.3)') ' ENSINF> REPLICA ',I-1, &
!MFC              ' HAS TEMPERATURE ', enstem(rep2t(i))
!MFC      enddo
!MFC      write(outu,'(a)') &
!MFC           ' ENSINF> **************************************************'
!MFC      write(outu,'(a)') &
!MFC           ' ENSINF> AND NOW, A WORD FROM EACH REPLICA ...'
!MFC   endif
!MFC   if (mynod == 0) &
!MFC        write (tbuf, '(a,i5,a,f8.3)') ' ENSINF>   REPLICA NODE ',  &
!MFC             whoiam,  &
!MFC             ': MY TEMPERATURE IS ', ENSMYT
!MFC   call ensprn(outu,tbuf,len(tbuf))
!MFC 
!MFC   return
!MFC end subroutine ensinf
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSTAT
!MFC !----------------------------------------------------------------------
!MFC subroutine enstat()
!MFC   !-----------------------------------------------------------------------
!MFC   ! Print replica exchange statistics at end of run
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use stream
!MFC   use dimens_fcm
!MFC   use mpi
!MFC   use number,only:zero
!MFC   implicit none
!MFC 
!MFC   integer i
!MFC   integer unum, ierror
!MFC   character(len=ensbfl) tbuf
!MFC   real(chm_real) ratio
!MFC 
!MFC   if (jrex.and.(iolev > 0)) then
!MFC      write(outu,'(a)') &
!MFC           ' ENSTAT> **************************************************'
!MFC      write(outu,'(a)') &
!MFC           ' ENSTAT> REPLICA EXCHANGE STATISTICS:'
!MFC      write(outu,'(a,i5)') &
!MFC           ' ENSINF>         NUMBER OF REPLICAS = ', NENSEM
!MFC      write(outu,'(a12,a10,a10,a10)') &
!MFC           'TRANSITION', 'ATTEMPTS', 'SUCCESSES', 'RATIO'
!MFC      do i=1,ensnsw
!MFC         if (ensatt(i) /= 0) then
!MFC            ratio = real(enssuc(i))/real(ensatt(i))
!MFC         else
!MFC            ratio = zero
!MFC         endif
!MFC         write(outu,'(2x,i3,a4,i3,2x,i10,i10,f10.3)') &
!MFC              ensisw(i),'<-->',ensjsw(i),ensatt(i),enssuc(i), &
!MFC              ratio
!MFC      enddo
!MFC      write(outu,'(a)') &
!MFC           ' ENSINF> **************************************************'
!MFC      call gflush(outu)
!MFC   endif
!MFC   return
!MFC end subroutine enstat
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSSWL
!MFC !----------------------------------------------------------------------
!MFC subroutine ensswl(xold,yold,zold,xnew,ynew,znew,xcomp,ycomp,zcomp,natom)
!MFC   !-----------------------------------------------------------------------
!MFC   ! Attempt replica swap -- swapping coordinates between nodes
!MFC   ! for leapfrog dynamics (dynamc.src integrator)
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use chm_types
!MFC   use ensemble
!MFC   use stream
!MFC   use dimens_fcm
!MFC   use number
!MFC   use reawri
!MFC   use nose_mod
!MFC   use energym
!MFC   use consta
!MFC   use coord
!MFC   use deriv
!MFC   use bases_fcm
!MFC   use image
!MFC   use comand
!MFC   use clcg_mod
!MFC   use parallel
!MFC   use mpi
!MFC   implicit none
!MFC   integer natom
!MFC   real(chm_real),dimension(natom) :: xold,yold,zold,xnew,ynew,znew,xcomp,ycomp,zcomp
!MFC   !
!MFC   real(chm_real) tmp6(6),pist(3),pistmp(3)
!MFC   integer ierror,i,j,k,ri,rj,partner,tmpi,ii,s
!MFC   integer, parameter :: myenel=4
!MFC   real(chm_real) myene(myenel),ensene(myenel*maxens)
!MFC   real(chm_real) eii,eij,ejj,eji,cdel,rnum,scafac,oldt,db
!MFC   real(chm_real) p_i,v_i,p_j,v_j
!MFC   real(chm_real) pivi,pivj,pjvi,pjvj
!MFC   logical qswp,sendfirst,debug
!MFC   integer atfrst,atlast
!MFC   integer swtmpi(maxswp), swtmpj(maxswp), ssave(maxswp), kk
!MFC   real(chm_real) randtbl(maxswp), tempd
!MFC   integer tempi,status
!MFC   logical swapped
!MFC 
!MFC   !  If I is this node and J is potential swap partner, then:
!MFC   !  MYENE(1) = EI(XI)
!MFC   !  MYENE(2) = EI(XJ)
!MFC 
!MFC   debug=.false.
!MFC   !  debug=.true.
!MFC 
!MFC   call vdgbr(xcomp,ycomp,zcomp,1)
!MFC 
!MFC   !     ROOT decides which nodes will attempt a swap
!MFC   setswap: if (mastermaster) then
!MFC      write(outu,'(a)') ' ******************************************************'
!MFC      write(outu,'(a)') ' **   Testing whether replicas need to be swapped ...**'
!MFC 
!MFC      if (ensswmethod == 0) then
!MFC         s=int(random(iseed)*(ensnsw))+1
!MFC         i=ensisw(s)
!MFC         j=ensjsw(s)
!MFC         write(outu,'(a,i3,a,i3)')  &
!MFC              ' ATTEMPTING TO SWAP REPLICAS ', I, ' AND ', J 
!MFC         call flush(outu)
!MFC         repswp(1:nensem)=0
!MFC         repswp(i)=j
!MFC         repswp(j)=i
!MFC         ensatt(s) = ensatt(s)+1
!MFC      else
!MFC         
!MFC         do k=1,ensnsw
!MFC            randtbl(k) = random(iseed)
!MFC            ssave(k) = k
!MFC         enddo
!MFC         
!MFC         swapped = .true.
!MFC         do while (swapped)
!MFC            do k=1,ensnsw-1
!MFC               if (randtbl(k) > randtbl(k+1)) then
!MFC                  tempd = randtbl(k)
!MFC                  randtbl(k) = randtbl(k+1)
!MFC                  randtbl(k+1) = tempd
!MFC                  tempi = ssave(k)
!MFC                  ssave(k) = ssave(k+1)
!MFC                  ssave(k+1) = tempi
!MFC                  swapped = .true.
!MFC               else
!MFC                  swapped = .false.
!MFC               endif
!MFC            enddo
!MFC         enddo
!MFC 
!MFC         do k=1,ensnsw
!MFC            swtmpi(k) = ensisw(ssave(k))
!MFC            swtmpj(k) = ensjsw(ssave(k))
!MFC            ensatt(ssave(k)) = ensatt(ssave(k)) + 1
!MFC         enddo
!MFC 
!MFC      endif
!MFC   endif setswap
!MFC 
!MFC   kk = ensnsw
!MFC   if (ensswmethod /= 0) kk = 1
!MFC 
!MFC   swaploop: do while (kk <= ensnsw)
!MFC      if (ensswmethod == 1) then
!MFC         s = ssave(kk)
!MFC      
!MFC         if (mastermaster) then
!MFC            repswp(1:nensem) = 0
!MFC            i = swtmpi(kk)
!MFC            j = swtmpj(kk)
!MFC            repswp(i)=j
!MFC            repswp(j)=i
!MFC            write(outu,'(a,i3,a,i3)') ' ATTEMPTING TO SWAP REPLICAS ', I, ' AND ', J
!MFC            call flush(outu)
!MFC         endif
!MFC      endif
!MFC 
!MFC      !     ROOT broadcasts swap array to other nodes
!MFC      call ens_bcast_all(repswp,nensem)
!MFC      !     TRIAL NODES evaluate old and new energies 
!MFC      if (repswp(whoiam+1) /= 0) then
!MFC         partner=repswp(whoiam+1)
!MFC         sendfirst=(partner > (whoiam+1))
!MFC         !        must evaluate energy in main coordinates to 
!MFC         !        accommmodate fast energy routines
!MFC         atfrst=1
!MFC         atlast=natom
!MFC         x(atfrst:atlast)=xcomp(atfrst:atlast)
!MFC         y(atfrst:atlast)=ycomp(atfrst:atlast)
!MFC         z(atfrst:atlast)=zcomp(atfrst:atlast)
!MFC         !        energy of local coords in local energy function
!MFC         !        (just to be on the safe side, re-evaluate energy)
!MFC         call update(comlyn,comlen,x,y,z,wmain,.false.,.false.,.true., &
!MFC              .false.,.true.,0,0,0,0,0,0,0)
!MFC         call energy(x,y,z,dx,dy,dz,bnbnd,bimag,0)
!MFC         myene(1)=eprop(epot)
!MFC         !        swap current coordinates
!MFC         if (debug) then
!MFC            write (outu,'(a14,2i3," x(1) = ",f12.6)') 'before: rep ',whoiam+1,mynod,x(1)
!MFC            call flush(outu)
!MFC         endif
!MFC         !        swap main coordinates
!MFC         call ensscv(x,ensbuf,natom,partner-1,sendfirst)
!MFC         call ensscv(y,ensbuf,natom,partner-1,sendfirst)
!MFC         call ensscv(z,ensbuf,natom,partner-1,sendfirst)
!MFC         if (debug) then
!MFC            write (outu,'(a14,2i3," x(1) = ",f12.6)') 'after: rep ',whoiam+1,mynod,x(1)
!MFC            call flush(outu)
!MFC         endif
!MFC         !        energy of remote coords in local energy function
!MFC         call update(comlyn,comlen,x,y,z,wmain,.false.,.false.,.true., &
!MFC              .false.,.true.,0,0,0,0,0,0,0)
!MFC         call energy(x,y,z,dx,dy,dz,bnbnd,bimag,0)
!MFC         myene(2)=eprop(epot)
!MFC         if (debug) then
!MFC            write (outu,'(a,i3,a,e12.5)') 'node ', whoiam+1,  &
!MFC                 ' olde = ', myene(1)
!MFC            write (outu,'(a,i3,a,e12.5)') 'node ', whoiam+1,  &
!MFC                 ' newe = ', myene(2)
!MFC            call flush(outu)
!MFC         endif
!MFC      endif
!MFC 
!MFC      !     ALL NODES share energy arrays (though only above data will be used)
!MFC 
!MFC      call ens_global_barrier(ierror)
!MFC      if (mynod == 0) then
!MFC         call mpi_allgather( &
!MFC              myene, myenel,mpi_real8, &
!MFC              ensene,myenel,mpi_real8, &
!MFC              comm_master,ierror)
!MFC      endif
!MFC      !     ROOT decides if TRIAL NODES will swap
!MFC      if (ensmasternod == 0) then
!MFC         qswp=.false.
!MFC         eii=ensene((i-1)*myenel+1)
!MFC         eij=ensene((i-1)*myenel+2)
!MFC         ejj=ensene((j-1)*myenel+1)
!MFC         eji=ensene((j-1)*myenel+2)
!MFC         write(outu,'(a,e12.5)') 'eii (kcal/mol) = ', eii
!MFC         write(outu,'(a,e12.5)') 'eij (kcal/mol) = ', eij
!MFC         write(outu,'(a,e12.5)') 'ejj (kcal/mol) = ', ejj
!MFC         write(outu,'(a,e12.5)') 'eji (kcal/mol) = ', eji
!MFC         call flush(outu)
!MFC 
!MFC         db = (eij-eii)/(enstem(i)*kboltz) &
!MFC              +(eji-ejj)/(enstem(j)*kboltz)
!MFC         write(outu,'(a,f12.6)') 'db = ', db
!MFC         write(outu,'(a,f12.6)') 'exp(-db) = ', exp(-db)
!MFC         call flush(outu)
!MFC         if (db <= zero) then
!MFC            qswp=.true.
!MFC         else
!MFC            rnum = random(iseed)
!MFC            write(outu,'(a,f12.6)') 'rnum = ', rnum
!MFC            if (rnum <= dexp(-db)) then
!MFC               qswp=.true.
!MFC            else
!MFC               qswp=.false.
!MFC            endif
!MFC         endif
!MFC         if (qswp) then
!MFC            write(outu,'(a,i3,a,i3)')  &
!MFC                 'swapping replicas for nodes ', i, &
!MFC                 ' and ', j
!MFC            call flush(outu)
!MFC            !              this just tracks which rep is at which temp for post-analysis
!MFC            ri=t2rep(i)
!MFC            rj=t2rep(j)
!MFC            t2rep(i)=rj
!MFC            t2rep(j)=ri
!MFC            rep2t(ri)=j
!MFC            rep2t(rj)=i
!MFC            enssuc(s)=enssuc(s)+1
!MFC            repswp(i)=-j
!MFC            repswp(j)=-i
!MFC         else
!MFC            repswp(i)=j
!MFC            repswp(j)=i
!MFC         endif
!MFC      endif
!MFC 
!MFC      !     ROOT broadcasts swap array to other nodes
!MFC      if(lmasternode)then 
!MFC         call psnd4_comm(comm_master,repswp,nensem)
!MFC         call psnd4_comm(comm_master,t2rep,nensem)
!MFC         call psnd4_comm(comm_master,rep2t,nensem)
!MFC         write(outu,'(a,10i4)') 'rep2t ', rep2t(1:nensem)
!MFC         write(outu,'(a,10i4)') 't2rep ', t2rep(1:nensem)
!MFC         call flush(outu)
!MFC      endif
!MFC      call psnd4_comm(comm_charmm,repswp,nensem)
!MFC 
!MFC      call psnd4_comm(comm_charmm,rep2t,nensem)
!MFC      call psnd4_comm(comm_charmm,t2rep,nensem)
!MFC 
!MFC      !     TRIAL NODES swap coordinates if success
!MFC      !     else also swap other dynamics arrays
!MFC      if (repswp(whoiam+1) /= 0) then
!MFC         partner=abs(repswp(whoiam+1))
!MFC         if (partner > (whoiam+1)) then
!MFC            sendfirst=.true.
!MFC         else
!MFC            sendfirst=.false.
!MFC         endif
!MFC         if (repswp(whoiam+1) < 0) then
!MFC            !           ... success
!MFC            !           swap coordinates used in leapfrog algorithm
!MFC            call ensscv(xold,ensbuf,natom,partner-1,sendfirst)
!MFC            call ensscv(yold,ensbuf,natom,partner-1,sendfirst)
!MFC            call ensscv(zold,ensbuf,natom,partner-1,sendfirst)
!MFC            call ensscv(xnew,ensbuf,natom,partner-1,sendfirst)
!MFC            call ensscv(ynew,ensbuf,natom,partner-1,sendfirst)
!MFC            call ensscv(znew,ensbuf,natom,partner-1,sendfirst)
!MFC            call ensscv(xcomp,ensbuf,natom,partner-1,sendfirst)
!MFC            call ensscv(ycomp,ensbuf,natom,partner-1,sendfirst)
!MFC            call ensscv(zcomp,ensbuf,natom,partner-1,sendfirst)
!MFC            !           update nb
!MFC            call update(comlyn,comlen,xcomp,ycomp,zcomp,wmain, &
!MFC                 .false.,.false.,.true., &
!MFC                 .false.,.true.,0,0,0,0,0,0,0)
!MFC            oldt=enstem(partner)
!MFC            scafac=sqrt(ensmyt/oldt)
!MFC            if (jensc) then
!MFC               do k=1,natom
!MFC                  xnew(k)=(two*scafac-one)*xold(k)+scafac*(xnew(k)-xold(k))
!MFC                  ynew(k)=(two*scafac-one)*yold(k)+scafac*(ynew(k)-yold(k))
!MFC                  znew(k)=(two*scafac-one)*zold(k)+scafac*(znew(k)-zold(k))
!MFC               enddo
!MFC            endif
!MFC 
!MFC         else
!MFC            call update(comlyn,comlen,xcomp,ycomp,zcomp,wmain, &
!MFC                 .false.,.false., &
!MFC                 .true.,.false.,.true.,0,0,0,0,0,0,0)
!MFC         endif
!MFC      endif
!MFC 
!MFC      kk = kk + 1
!MFC   enddo swaploop
!MFC   
!MFC 
!MFC   return
!MFC end subroutine ensswl
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSSCV
!MFC !----------------------------------------------------------------------
!MFC subroutine ensscv(array,tmpv,natom,partner,sendfirst)
!MFC   !-----------------------------------------------------------------------
!MFC   !     client-to-client array swapping for coordinate version of
!MFC   !     replica exchange
!MFC   !-----------------------------------------------------------------------
!MFC   ! In Parallel Ensemble, Send/Receive is between nodes with same mynod:
!MFC   ! 0   1   2  ... numnod-1
!MFC   ! |   |   |         |
!MFC   ! 0   1   2  ... numnod-1
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use stream
!MFC   use dimens_fcm
!MFC   use number
!MFC   use reawri
!MFC   use nose_mod
!MFC   use energym
!MFC   use parallel, only:mynod,numnod,comm_charmm
!MFC   use mpi
!MFC   implicit none
!MFC   integer,intent(in) :: natom,partner
!MFC   real(chm_real),dimension(natom) :: array,tmpv
!MFC   logical,intent(in) :: sendfirst
!MFC 
!MFC   integer ierror,stat(mpi_status_size)
!MFC 
!MFC   if (sendfirst) then
!MFC      if(lmasternode)then
!MFC         call mpi_send(array,natom,mpi_real8,partner,10, &
!MFC              comm_master,ierror)
!MFC         call mpi_recv(tmpv,natom,mpi_real8,partner,11, &
!MFC              comm_master,stat,ierror)
!MFC      endif
!MFC   ELSE
!MFC      if(lmasternode)then
!MFC         call mpi_recv(tmpv,natom,mpi_real8,partner,10, &
!MFC              comm_master,stat,ierror)
!MFC         call mpi_send(array,natom,mpi_real8,partner,11, &
!MFC              comm_master,ierror)
!MFC      endif
!MFC   endif
!MFC   call psnd8_comm(comm_charmm,tmpv,natom)
!MFC 
!MFC      array=tmpv
!MFC      return
!MFC end subroutine ensscv
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSOUT
!MFC !----------------------------------------------------------------------
!MFC SUBROUTINE ENSOUT()
!MFC   !-----------------------------------------------------------------------
!MFC   ! write rex output
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use stream
!MFC   use dimens_fcm
!MFC   implicit none
!MFC   integer i
!MFC   if (mastermaster .and. iolev >= 0) then
!MFC      if (rep2to > 0) then
!MFC         write(rep2to,'(10i5)') rep2t(1:nensem)
!MFC         call gflush(rep2to)
!MFC      endif
!MFC      if (t2repo > 0) then 
!MFC         write(t2repo,'(10i5)') t2rep(1:nensem)
!MFC         call gflush(t2repo)
!MFC      endif
!MFC   endif
!MFC   return
!MFC end subroutine ensout
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSAVE
!MFC !----------------------------------------------------------------------
!MFC SUBROUTINE ENSAVE(PARMS,BUFF,PARML)
!MFC   !-----------------------------------------------------------------------
!MFC   ! Simple linear ensemble averaging ...
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use mpi
!MFC   implicit none
!MFC 
!MFC   real(chm_real) PARMS(*), BUFF(*)
!MFC   INTEGER PARML,IERROR
!MFC   INTEGER I,M
!MFC   real(chm_real) AVE
!MFC   IF (NENSEM == 1) RETURN
!MFC   CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!MFC   CALL MPI_ALLGATHER(PARMS(1),PARML,MPI_REAL8,BUFF(1), &
!MFC        PARML,MPI_REAL8,MPI_COMM_WORLD,IERROR)
!MFC   DO I=1,PARML
!MFC      AVE = 0.0
!MFC      DO M=1,NENSEM
!MFC         AVE = AVE + BUFF((M-1)*PARML+I)
!MFC      ENDDO
!MFC      PARMS(I) = AVE/REAL(NENSEM)
!MFC   ENDDO
!MFC   RETURN
!MFC END SUBROUTINE ENSAVE
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSAV3
!MFC !----------------------------------------------------------------------
!MFC SUBROUTINE ENSAV3(PARMS,BUFF,PARML)
!MFC   !-----------------------------------------------------------------------
!MFC   ! Third power (noe) averaging
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use mpi
!MFC   implicit none
!MFC   real(chm_real) PARMS(*), BUFF(*)
!MFC   INTEGER PARML,IERROR
!MFC   INTEGER I,M
!MFC   real(chm_real) AVE
!MFC   IF (NENSEM == 1) RETURN
!MFC   CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!MFC   CALL MPI_ALLGATHER(PARMS(1),PARML,MPI_REAL8,BUFF(1), &
!MFC        PARML,MPI_REAL8,MPI_COMM_WORLD,IERROR)
!MFC   DO I=1,PARML
!MFC      AVE = 0.0
!MFC      DO M=1,NENSEM
!MFC         AVE = AVE + (BUFF((M-1)*PARML+I))**(-3.0)
!MFC      ENDDO
!MFC      PARMS(I) = (AVE/REAL(NENSEM))**(-1.0/3.0)
!MFC   ENDDO
!MFC   RETURN
!MFC END SUBROUTINE ENSAV3
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSAV6
!MFC !----------------------------------------------------------------------
!MFC SUBROUTINE ENSAV6(PARMS,BUFF,PARML)
!MFC   !-----------------------------------------------------------------------
!MFC   ! Sixth power (noe) averaging
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use mpi
!MFC   implicit none
!MFC   real(chm_real) PARMS(*), BUFF(*)
!MFC   INTEGER PARML,IERROR
!MFC   INTEGER I,M
!MFC   real(chm_real) AVE
!MFC   IF (NENSEM == 1) RETURN
!MFC   CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!MFC   CALL MPI_ALLGATHER(PARMS(1),PARML,MPI_REAL8,BUFF(1), &
!MFC        PARML,MPI_REAL8,MPI_COMM_WORLD,IERROR)
!MFC   DO I=1,PARML
!MFC      AVE = 0.0
!MFC      DO M=1,NENSEM
!MFC         AVE = AVE + (BUFF((M-1)*PARML+I))**(-6.0)
!MFC      ENDDO
!MFC      PARMS(I) = (AVE/REAL(NENSEM))**(-1.0/6.0)
!MFC   ENDDO
!MFC   RETURN
!MFC END SUBROUTINE ENSAV6
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSS2
!MFC !----------------------------------------------------------------------
!MFC SUBROUTINE ENSS2(BXIJ,BYIJ,BZIJ,XIJBUF,YIJBUF,ZIJBUF,NPARM,RIJ,S2)
!MFC   !-----------------------------------------------------------------------
!MFC   ! Calculate order parameter over ensemble
!MFC   ! Takes r_ij vector components (BXIJ etc for this replica) and calculates
!MFC   ! S2 over ensemble. Also calculates components of forces, and returns
!MFC   ! this in BXIJ etc.
!MFC   !-----------------------------------------------------------------------
!MFC   use chm_kinds
!MFC   use ensemble
!MFC   use mpi
!MFC   implicit none
!MFC   real(chm_real) BXIJ(*),BYIJ(*),BZIJ(*),XIJBUF(*),YIJBUF(*),ZIJBUF(*)
!MFC   real(chm_real) RIJ(*),S2(*)
!MFC   INTEGER NPARM
!MFC   !     local vbls
!MFC   INTEGER IERROR,I,M
!MFC   real(chm_real) SX2,SY2,SZ2,SXY,SXZ,SYZ,XIJ,YIJ,ZIJ
!MFC   !
!MFC   IF (NENSEM == 1) RETURN
!MFC   CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!MFC   CALL MPI_ALLGATHER(BXIJ(1),NPARM,MPI_REAL8,XIJBUF(1), &
!MFC        NPARM,MPI_REAL8,MPI_COMM_WORLD,IERROR)
!MFC   CALL MPI_ALLGATHER(BYIJ(1),NPARM,MPI_REAL8,YIJBUF(1), &
!MFC        NPARM,MPI_REAL8,MPI_COMM_WORLD,IERROR)
!MFC   CALL MPI_ALLGATHER(BZIJ(1),NPARM,MPI_REAL8,ZIJBUF(1), &
!MFC        NPARM,MPI_REAL8,MPI_COMM_WORLD,IERROR)
!MFC   do i=1,nparm
!MFC      sx2 = 0.0
!MFC      sy2 = 0.0
!MFC      sz2 = 0.0
!MFC      sxy = 0.0
!MFC      sxz = 0.0
!MFC      syz = 0.0
!MFC      do m=1,nensem
!MFC         xij = xijbuf((m-1)*nparm+i)
!MFC         yij = yijbuf((m-1)*nparm+i)
!MFC         zij = zijbuf((m-1)*nparm+i)
!MFC         sx2 = sx2 + xij**2
!MFC         sy2 = sy2 + yij**2
!MFC         sz2 = sz2 + zij**2
!MFC         sxy = sxy + xij*yij
!MFC         sxz = sxz + xij*zij
!MFC         syz = syz + yij*zij
!MFC      end do
!MFC      sx2 = sx2 / real(nensem)
!MFC      sy2 = sy2 / real(nensem)
!MFC      sz2 = sz2 / real(nensem)
!MFC      sxy = sxy / real(nensem)
!MFC      sxz = sxz / real(nensem)
!MFC      syz = syz / real(nensem)
!MFC      !        order parameter:
!MFC      s2(i) = (1.5/rij(i)**4)*(sx2**2+sy2**2+sz2**2+2.0*sxy**2 &
!MFC           + 2.0*sxz**2 + 2.0*syz**2) - 0.5
!MFC      !        components of force:
!MFC      xij = bxij(i)
!MFC      yij = byij(i)
!MFC      zij = bzij(i)
!MFC      bxij(i) = sx2*xij+sxy*yij+sxz*zij
!MFC      byij(i) = sy2*yij+sxy*xij+syz*zij
!MFC      bzij(i) = sz2*zij+sxz*xij+syz*yij
!MFC   enddo
!MFC   return
!MFC end subroutine enss2
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          EAVGSETUP
!MFC !----------------------------------------------------------------------
!MFC subroutine eavgsetup(comlyn,comlen)
!MFC   !-----------------------------------------------------------------------
!MFC   ! the nodes share their forces and energy and and do exp(-H)=exp(-HA)+exp(-HB)
!MFC   ! averaging
!MFC   !
!MFC   ! ... setup & allocate storage at beginning
!MFC   !
!MFC   use chm_kinds
!MFC   use dimens_fcm
!MFC   use string
!MFC   use stream
!MFC   use psf
!MFC   use ensemble
!MFC   use number
!MFC   use memory
!MFC   implicit none
!MFC   CHARACTER(len=*) COMLYN
!MFC   INTEGER COMLEN
!MFC   !
!MFC   ! local:
!MFC   INTEGER BUFLEN,I
!MFC   LOGICAL JOFFSET
!MFC   !
!MFC   QENSEXP = .true.
!MFC   BUFLEN = NATOM*NENSEM
!MFC   IF (IOLEV > 0) THEN
!MFC      WRITE(OUTU,'(A)') ' WILL DO EXPONENTIAL AVERAGING OF ENERGY'
!MFC   ENDIF
!MFC   !     write output?
!MFC   ENSEXPU=GTRMI(COMLYN,COMLEN,'UNIT',-1)
!MFC   IF ((IOLEV > 0).AND.(ENSEXPU > 0)) THEN
!MFC      WRITE(OUTU,'(A,I5)') ' WRITING EXP ENERGY OUTPUT TO UNIT', &
!MFC           ENSEXPU
!MFC   ENDIF
!MFC   !     set up beta for averaging
!MFC   ENSEBETA=GTRMF(COMLYN,COMLEN,'BETA',MINONE)
!MFC   IF (IOLEV > 0) THEN
!MFC      WRITE(OUTU,'(A,F12.5)') ' WEIGHTING ENERGIES USING BETA = ', &
!MFC           ENSEBETA
!MFC   ENDIF
!MFC   IF (INDXA(COMLYN, COMLEN, 'OFFSET')  >  0) THEN
!MFC      JOFFSET =.TRUE.
!MFC   ELSE
!MFC      JOFFSET =.FALSE.
!MFC   ENDIF
!MFC   IF (JOFFSET) THEN
!MFC      DO I=1,NENSEM
!MFC         EXPOFF(I) = NEXTF(COMLYN,COMLEN)
!MFC      ENDDO
!MFC   ELSE
!MFC      DO I=1,NENSEM
!MFC         EXPOFF(I) = 0.0
!MFC      ENDDO
!MFC   ENDIF
!MFC   IF (IOLEV > 0) THEN
!MFC      WRITE(OUTU,'(A)') ' ENERGY OFFSETS: '
!MFC      DO I=1,NENSEM
!MFC         WRITE(OUTU,'(A,I3,A,F12.5)') ' NODE = ', I-1,  &
!MFC              ' OFFSET = ', EXPOFF(I)
!MFC      ENDDO
!MFC   ENDIF
!MFC   IF (INDXA(COMLYN, COMLEN, 'BOLTZ')  >  0) THEN
!MFC      QEXPBW =.TRUE.
!MFC      IF (IOLEV > 0) THEN
!MFC         WRITE(OUTU,'(A)') 'USING ALTERNATIVE EXPONENTIAL AVE:'
!MFC         WRITE(OUTU,'(A)') '         exp(-2BUa)+exp(-2BUb)'
!MFC         WRITE(OUTU,'(A)') 'exp(-BU)=---------------------'
!MFC         WRITE(OUTU,'(A)') '         exp(-BUa)+exp(-BUb)'
!MFC      ENDIF
!MFC   ELSE
!MFC      QEXPBW =.FALSE.
!MFC      IF (IOLEV > 0) THEN
!MFC         WRITE(OUTU,'(A)') 'USING STANDARD EXPONENTIAL AVE:'
!MFC         WRITE(OUTU,'(A)') 'exp(-BU)=exp(-BUa)+exp(-BUb)'
!MFC      ENDIF
!MFC   ENDIF
!MFC   !     allocate space for MPI receive buffers
!MFC   call chmalloc('ensemble.src','EAVGSETUP','ENSDX',BUFLEN,crl=ENSDX)
!MFC   call chmalloc('ensemble.src','EAVGSETUP','ENSDY',BUFLEN,crl=ENSDy)
!MFC   call chmalloc('ensemble.src','EAVGSETUP','ENSDZ',BUFLEN,crl=ENSDZ)
!MFC   call chmalloc('ensemble.src','EAVGSETUP','ENSH',BUFLEN,crl=ENSH)
!MFC   !
!MFC   return
!MFC end subroutine eavgsetup
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSEXPAVG
!MFC !----------------------------------------------------------------------
!MFC subroutine ensexpavg(dxk,dyk,dzk,hk,natom)
!MFC   !-----------------------------------------------------------------------
!MFC   ! the nodes share their forces and energy and and do 
!MFC   ! exp(-H)=exp(-HA)+exp(-HB) averaging
!MFC   !
!MFC   ! ... calculate energy & forces at each time step
!MFC 
!MFC   use chm_kinds
!MFC   use dimens_fcm
!MFC   use stream
!MFC   use ensemble
!MFC   implicit none
!MFC 
!MFC   real(chm_real) dxk(*),dyk(*),dzk(*),hk
!MFC   integer natom
!MFC   ! local:
!MFC   integer buflen
!MFC   !
!MFC   buflen = natom*nensem
!MFC   call ensexpavg2(dxk,dyk,dzk,hk,natom,ensdx,ensdy,ensdz,ensh,buflen)
!MFC   return
!MFC end subroutine ensexpavg
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ENSEXPAVG2
!MFC !----------------------------------------------------------------------
!MFC subroutine ensexpavg2(dxk,dyk,dzk,hk,natom, &
!MFC      rcvdx,rcvdy,rcvdz,rcvh,buflen)
!MFC   !-----------------------------------------------------------------------
!MFC   ! the nodes share their forces and energy and and do exp(-H)=exp(-HA)+exp(-HB)
!MFC   ! averaging
!MFC   !
!MFC   ! ... calculate energy & forces at each time step
!MFC   !
!MFC   use chm_kinds
!MFC   use stream
!MFC   use ensemble
!MFC   use contrl
!MFC   use number
!MFC   use mpi
!MFC   implicit none
!MFC   !
!MFC   real(chm_real) DXK(*),DYK(*),DZK(*),RCVDX(*),RCVDY(*),RCVDZ(*)
!MFC   real(chm_real) RCVH(*),HK
!MFC   INTEGER NATOM,BUFLEN
!MFC   !     local
!MFC   INTEGER IERROR,I,M
!MFC   real(chm_real) EXPSUM,EXP2SUM,DXI,DYI,DZI,MINE
!MFC   real(chm_real) DXI2,DYI2,DZI2
!MFC   real(chm_real) TMPH(1)
!MFC   real(chm_real) EXPWT(MAXENS)
!MFC   !
!MFC   TMPH(1) = HK
!MFC   CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!MFC   CALL MPI_ALLGATHER(DXK(1),NATOM,MPI_REAL8,RCVDX(1), &
!MFC        NATOM,MPI_REAL8,MPI_COMM_WORLD,IERROR)
!MFC   CALL MPI_ALLGATHER(DYK(1),NATOM,MPI_REAL8,RCVDY(1), &
!MFC        NATOM,MPI_REAL8,MPI_COMM_WORLD,IERROR)
!MFC   CALL MPI_ALLGATHER(DZK(1),NATOM,MPI_REAL8,RCVDZ(1), &
!MFC        NATOM,MPI_REAL8,MPI_COMM_WORLD,IERROR)
!MFC   CALL MPI_ALLGATHER(TMPH(1),1,MPI_REAL8,RCVH(1), &
!MFC        1,MPI_REAL8,MPI_COMM_WORLD,IERROR)
!MFC   !     first calculate mean to use as an offset to avoid numerical
!MFC   !     errors in exp(-beta U) for very small (large negative) U.
!MFC   MINE=RCVH(1)+EXPOFF(1)
!MFC   DO I=2,NENSEM
!MFC      MINE = MIN(MINE,RCVH(I)+EXPOFF(I))
!MFC   ENDDO
!MFC   !     convert entries HK(I) to EXP(-beta*HK(I))
!MFC   EXPSUM=0.0
!MFC   EXP2SUM=0.0
!MFC   DO I=1,NENSEM
!MFC      EXPWT(I)=DEXP(-ENSEBETA*(RCVH(I)+EXPOFF(I)-MINE))
!MFC      EXPSUM=EXPSUM+EXPWT(I)
!MFC      IF (QEXPBW) EXP2SUM=EXP2SUM+(EXPWT(I))**2
!MFC   ENDDO
!MFC   !     calculate forces
!MFC   IF (QEXPBW) THEN
!MFC      !        write (outu,*) 'DOING EXPBW'
!MFC      !        call flush(outu)
!MFC      DO I=1,NATOM
!MFC         DXI = 0.0
!MFC         DYI = 0.0
!MFC         DZI = 0.0
!MFC         DXI2 = 0.0
!MFC         DYI2 = 0.0
!MFC         DZI2 = 0.0
!MFC         DO M=1,NENSEM
!MFC            DXI2=DXI2+RCVDX((M-1)*NATOM+I)*(EXPWT(M)**2)
!MFC            DYI2=DYI2+RCVDY((M-1)*NATOM+I)*(EXPWT(M)**2)
!MFC            DZI2=DZI2+RCVDZ((M-1)*NATOM+I)*(EXPWT(M)**2)
!MFC            DXI =DXI +RCVDX((M-1)*NATOM+I)*EXPWT(M)
!MFC            DYI =DYI +RCVDY((M-1)*NATOM+I)*EXPWT(M)
!MFC            DZI =DZI +RCVDZ((M-1)*NATOM+I)*EXPWT(M)
!MFC         ENDDO
!MFC         DXK(I) = TWO*DXI2/EXP2SUM-DXI/EXPSUM
!MFC         DYK(I) = TWO*DYI2/EXP2SUM-DYI/EXPSUM
!MFC         DZK(I) = TWO*DZI2/EXP2SUM-DZI/EXPSUM
!MFC      ENDDO
!MFC      HK = DLOG(EXP2SUM/EXPSUM)/(-ENSEBETA)+MINE
!MFC   ELSE
!MFC      DO I=1,NATOM
!MFC         DXI = 0.0
!MFC         DYI = 0.0
!MFC         DZI = 0.0
!MFC         DO M=1,NENSEM
!MFC            DXI=DXI+RCVDX((M-1)*NATOM+I)*EXPWT(M)
!MFC            DYI=DYI+RCVDY((M-1)*NATOM+I)*EXPWT(M)
!MFC            DZI=DZI+RCVDZ((M-1)*NATOM+I)*EXPWT(M)
!MFC         ENDDO
!MFC         DXK(I) = DXI/EXPSUM
!MFC         DYK(I) = DYI/EXPSUM
!MFC         DZK(I) = DZI/EXPSUM
!MFC      ENDDO
!MFC      HK = DLOG(EXPSUM)/(-ENSEBETA)+MINE
!MFC   ENDIF
!MFC   !     write output ...
!MFC   IF ((IOLEV >= 0).AND.(ENSEXPU > 0)) THEN
!MFC      IF(MOD(MDSTEP,NPRINT) == 0) THEN
!MFC         WRITE(ENSEXPU,'(I10$)') MDSTEP
!MFC         DO I=1,NENSEM
!MFC            WRITE(ENSEXPU,'(1X,F12.5$)') RCVH(I)
!MFC         ENDDO
!MFC         WRITE(ENSEXPU,'(1X,F12.5)') HK
!MFC      ENDIF
!MFC   ENDIF
!MFC 
!MFC   return
!MFC end subroutine ensexpavg2
!MFC 
!MFC !----------------------------------------------------------------------
!MFC !          ERDCMD
!MFC !----------------------------------------------------------------------
!MFC subroutine erdcmd(comlyn,mxcms2,comlen,unit,eof, &
!MFC      qcmprs,qprint,echost)
!MFC   !-----------------------------------------------------------------------
!MFC   !     parallel input version of rdcmnd
!MFC   use chm_kinds
!MFC   use dimens_fcm
!MFC   use cmdpar,only:parse1
!MFC   use string
!MFC   use stream
!MFC   use rtf,only:ucase
!MFC 
!MFC 
!MFC   implicit none
!MFC   integer mxcms2,comlen,unit
!MFC   character(len=*) comlyn
!MFC   !
!MFC   integer enrlen
!MFC   integer, parameter :: enrmax=20
!MFC   character(len=enrmax) enrst
!MFC   integer cardln
!MFC   integer, parameter :: mxcard=200
!MFC   character(len=mxcard) card
!MFC   character(len=1) stemp
!MFC   integer wdlen,iend,ipar,i,j
!MFC   logical eof,qprint,qcmprs
!MFC   character(len=*) echost
!MFC   character(len=1) hyphen/'-'/,exclmk/'!'/,atmark/'@'/,sdblq/'"'/
!MFC   character(len=1) sques/'?'/
!MFC   !
!MFC   if(eof) return
!MFC 
!MFC   !--- clean up the previous command line in case of something left
!MFC   comlyn=' '
!MFC   comlen=0
!MFC   if(unit < 0) then
!MFC      eof=.true.
!MFC      return
!MFC   endif
!MFC 
!MFC   if(qprint .and. prnlev >= 3) write(outu,'("  ")')
!MFC 
!MFC   loop0: do while(.true.)
!MFC      iend=1
!MFC      loop1: do while(iend == 1)
!MFC         read(unit,'(a)',end=9) card
!MFC         cardln=mxcard
!MFC         call trime(card,cardln)
!MFC         if(cardln == 0) cardln=1
!MFC         if(qprint.and.prnlev >= 3) write(outu,'(1x,a8,3x,a)') echost,card(1:cardln)
!MFC         iend=indx(card,cardln,exclmk,1)
!MFC      enddo loop1    ! if (iend == 1) goto 1
!MFC      
!MFC      if (iend /= 0) then
!MFC         cardln=iend-1
!MFC         call trime(card,cardln)
!MFC      endif
!MFC      
!MFC #if KEY_CFF==1
!MFC      if (ucase) then      
!MFC #endif
!MFC         call cnvtuc(card,cardln)
!MFC #if KEY_CFF==1
!MFC      endif      
!MFC #endif
!MFC 
!MFC      if(qcmprs) call cmprst(card,cardln)
!MFC      if(cardln == 0) exit loop0
!MFC      if(card(cardln:cardln) /= hyphen) exit loop0
!MFC      if(cardln == 1) cycle loop0
!MFC      if(comlen+cardln-1 > mxcms2) then
!MFC         call wrndie(-1,'<RDCMND>','Command line too long: truncated.')
!MFC      endif
!MFC      call addst(comlyn,mxcms2,comlen,card,cardln-1)
!MFC   enddo loop0
!MFC 
!MFC   if(comlen+cardln > mxcms2) then
!MFC      call wrndie(-1,'<RDCMND>','Command line too long: truncated.')
!MFC   endif
!MFC   call addst(comlyn,mxcms2,comlen,card,cardln)
!MFC 
!MFC   if(qcmprs) ffour=comlyn(1:4)
!MFC   call parse1(comlyn,mxcms2,comlen,qprint)
!MFC   if(comlen == 0)return
!MFC 
!MFC   !     Before returning the string make any energy
!MFC   !     substitutions that may be required.
!MFC 
!MFC 2020 FORMAT(' ERDCMD substituted energy or value "',80A1)
!MFC 2030 FORMAT(' ERDCMD: can not substitute energy "',80A1)
!MFC 
!MFC   ipar=indx(comlyn,comlen-1,sques,1)
!MFC   loop400: do while(ipar > 0)
!MFC      call copsub(enrst,enrmax,enrlen,comlyn,ipar+1, &
!MFC           min(comlen,enrmax+ipar))
!MFC      call subenr(wdlen,enrlen,enrst,enrmax)
!MFC      if (enrlen > 0) then
!MFC         if (qprint .and. prnlev >= 3) then
!MFC            write(outu,2020) (comlyn(j:j),j=ipar,ipar+wdlen),sdblq, &
!MFC                 ' ','t','o',' ',sdblq,(enrst(j:j),j=1,enrlen),sdblq
!MFC         endif
!MFC         call copsub(scrtch,scrmax,scrlen,comlyn,ipar+wdlen+1,comlen)
!MFC         comlen=ipar-1
!MFC         call addst(comlyn,mxcms2,comlen,enrst,enrlen)
!MFC         call addst(comlyn,mxcms2,comlen,scrtch,scrlen)
!MFC         ipar=indx(comlyn,comlen-1,sques,1)
!MFC      else
!MFC         !          fix may-90/ln, we rather want the whole parameter name
!MFC         if(wrnlev >= 2) write(outu,2030) &
!MFC              (comlyn(j:j),j=ipar,ipar+wdlen),sdblq
!MFC         ipar=0
!MFC      endif
!MFC   enddo loop400
!MFC 
!MFC   !--- Trim the command line before return to the caller
!MFC   call trima(comlyn,comlen)
!MFC   return
!MFC 9 eof=.true.
!MFC 
!MFC   return
!MFC end subroutine erdcmd

!MFC !=======================================================================
!MFC !     The routines below were lifted from the parallel code
!MFC !=======================================================================
!MFC !----------------------------------------------------------------------
!MFC !          PSYNC_ENS
!MFC !----------------------------------------------------------------------
!MFC subroutine psync_ens()
!MFC   !-----------------------------------------------------------------------
!MFC   !
!MFC   !     This is a wrapper routine for global sync.
!MFC   !
!MFC   use chm_kinds
!MFC   !
!MFC   use dimens_fcm
!MFC   use ensemble
!MFC #if KEY_CMPI==0
!MFC   use mpi    
!MFC #endif
!MFC   implicit none
!MFC   integer status
!MFC 
!MFC #if KEY_PARINFNTY==1
!MFC   if(qinfinity) return                 
!MFC #endif
!MFC   if(nensem == 1) return
!MFC   if (slave_ens) return
!MFC #if KEY_CMPI==1
!MFC   call cmpi_barrier(comm_master,status)
!MFC #endif 
!MFC 
!MFC   return
!MFC end subroutine psync_ens
!MFC 


!MFC #else /* (ensemble_main)*/
!MFC 
!MFC 
!MFC !====================== NO ENSEMBLES ... ==============================
!MFC subroutine enscmd(comlyn,comlen)
!MFC   !----------------------------------------------------------------------
!MFC   use stream
!MFC   character(len=*) comlyn
!MFC   integer comlen
!MFC   write (outu, '(a)') 'NO ENSEMBLE CODE COMPILED'
!MFC   return
!MFC end subroutine enscmd
!MFC 
!MFC subroutine ensfin
!MFC   return
!MFC end subroutine ensfin
!MFC 
!MFC subroutine enschk
!MFC   return
!MFC end subroutine enschk
!MFC 

subroutine masterrecsenin(tonode,tag,rbuf,rlen,sbuf,slen,mtype)
  use chm_kinds
  use mpi
  use parallel
  use repd_ensemble
  implicit none
  integer,intent(in) :: mtype,tonode,tag,rlen,slen
  integer,dimension(rlen),intent(inout) :: rbuf
  integer,dimension(slen),intent(inout) :: sbuf
  integer req(2),ierr,wtall,status_array(mpi_status_size,2)


  wtall=2
  ! print '("masterrecsen ",i3,a,i3,a,i3,e15.3,i8)',ensmasternod," sending to ",tonode," tag ",tag,sbuf(1),slen
  call mpi_isend(sbuf,slen,mpi_integer,tonode,tag,comm_master,req(1),ierr)
  if(ierr /= mpi_success)  call wrndie(-4,'<masterrecsen>','mpi isend error')
  call mpi_irecv(rbuf,rlen,mpi_integer,tonode,tag,comm_master,req(2),ierr)
  if(ierr /= mpi_success)  call wrndie(-4,'<masterrecsen>','mpi irecv error')
  call mpi_waitall(wtall,req,status_array,ierr)
  if(ierr /= mpi_success)  call wrndie(-4,'<masterrecsen>','mpi waitall error')
  ! print '("masterrecsen ",i3,a,i3,a,i3,e15.3,i8)',ensmasternod," receive from ",tonode," tag ",tag,rbuf(1),rlen
! print *,ensmasternod," received",tag,rbuf(1),rlen," from ", tonode
  return
end subroutine masterrecsenin

subroutine masterrecsenrl(tonode,tag,rbuf,rlen,sbuf,slen,mtype)
  use chm_kinds
  use mpi
  use parallel
  use repd_ensemble
  implicit none
  integer,intent(in) :: mtype,tonode,tag,rlen,slen
  real(chm_real),dimension(rlen),intent(inout) :: rbuf
  real(chm_real),dimension(slen),intent(inout) :: sbuf
  integer req(2),ierr,wtall,status_array(mpi_status_size,2)


  wtall=2
  ! print '("masterrecsen ",i3,a,i3,a,i3,e15.3,i8)',ensmasternod," sending to ",tonode," tag ",tag,sbuf(1),slen
  call mpi_isend(sbuf,slen,mpi_real8,tonode,tag,comm_master,req(1),ierr)
  if(ierr /= mpi_success)  call wrndie(-4,'<masterrecsen>','mpi isend error')
  call mpi_irecv(rbuf,rlen,mpi_real8,tonode,tag,comm_master,req(2),ierr)
  if(ierr /= mpi_success)  call wrndie(-4,'<masterrecsen>','mpi irecv error')
  call mpi_waitall(wtall,req,status_array,ierr)
  if(ierr /= mpi_success)  call wrndie(-4,'<masterrecsen>','mpi waitall error')
  ! print '("masterrecsen ",i3,a,i3,a,i3,e15.3,i8)',ensmasternod," receive from ",tonode," tag ",tag,rbuf(1),rlen
! print *,ensmasternod," received",tag,rbuf(1),rlen," from ", tonode
  return
end subroutine masterrecsenrl

subroutine masterrecsenrl_1(tonode,tag,rbuf,rlen,sbuf,slen,qsend)
  use chm_kinds
  use mpi
  use parallel
  use repd_ensemble
  implicit none
  integer,intent(in) :: tonode,tag,rlen,slen
  real(chm_real),dimension(rlen),intent(inout) :: rbuf
  real(chm_real),dimension(slen),intent(inout) :: sbuf
  logical,intent(in) :: qsend
  ! integer req,ierr,status_array(mpi_status_size,2) ! Josh/RLH
  integer req,ierr,status_array(mpi_status_size)

!  print '("masterrecsenrl_1 ",i3,a,i3,a,i3,e15.3,i8)',ensmasternod," sending to ",tonode," tag ",tag,sbuf(1),slen
  if(qsend) then
     call mpi_isend(sbuf,slen,mpi_real8,tonode,tag,comm_master,req,ierr)
     if(ierr /= mpi_success)  call wrndie(-4,'<masterrecsen>','mpi isend error')
  else
     call mpi_irecv(rbuf,rlen,mpi_real8,tonode,tag,comm_master,req,ierr)
     if(ierr /= mpi_success)  call wrndie(-4,'<masterrecsen>','mpi irecv error')
  endif
  call mpi_wait(req,status_array,ierr)
!  call mpi_waitall(wtall,req,status_array,ierr)
  if(ierr /= mpi_success)  call wrndie(-4,'<masterrecsen>','mpi waitall error')
  ! print '("masterrecsen ",i3,a,i3,a,i3,e15.3,i8)',ensmasternod," receive from ",tonode," tag ",tag,rbuf(1),rlen
  ! print *,ensmasternod," received",tag,rbuf(1),rlen," from ", tonode
  return
end subroutine masterrecsenrl_1

!===========================================================
!             Integer
!===========================================================
subroutine mastersenint(to,tag,buf,len)
  use chm_kinds
  use mpi                          
  use repd_ensemble
  implicit none
  
  integer,intent(in) :: to, tag, len
  integer,intent(in),dimension(len) ::  buf
  
  integer ierr,req,m_stat(mpi_status_size)
  if(len <= 0) return

  ! print '("senint ",i3,a,i3,a,i3,i10,i8)',ensmasternod," send to ",to," tag ",tag,buf(1),len

  call mpi_isend(buf,len,mpi_integer,to,tag,comm_master,req,ierr)
  if(ierr /= mpi_success) call wrndie(-4,'repd_ensemble.src<mastersenint>', &
       'mpi isend error in mastersenint')
  call mpi_wait(req,m_stat,ierr)
  ! print '("senint ",i3,a,i3,a,i3,"  success ")',ensmasternod," send to ",to," tag ",tag
  if(ierr /= mpi_success) call wrndie(-4,'repd_ensemble.src<mastersenint>', &
       'mpi wait error in mastersenint')
  return
end subroutine mastersenint

subroutine masterrecint(from,tag,buf,len)
  use chm_kinds
  use mpi                          
  use repd_ensemble
  implicit none

  integer,intent(in) :: from, tag, len
  integer,intent(out),dimension(len) :: buf

  integer m_stat(mpi_status_size),ierr,req
  if(len <= 0) return

  ! print '("recint ",i3,a,i3,a,i3,a)',ensmasternod," rec frm ",from," tag ",tag,"posted"

  call mpi_irecv(buf,len,mpi_integer,from,tag,comm_master,req, ierr)
  if(ierr /= mpi_success)  &
       call wrndie(-4,'repd-ensemble.src<masterrecint>','mpi irecv error in masterrecint')
  ! print '("recint ",i3,a,i3,a,i3,i10,i8)',ensmasternod," waiting   ",from," tag ",tag
  call mpi_wait(req,m_stat,ierr)
  ! print '("recint ",i3,a,i3,a,i3,i10,i8)',ensmasternod," done    ",from," tag ",tag,buf(1),len
  if(ierr /= mpi_success)  &
       call wrndie(-4,'repd-ensemble.src<masterrecint>','mpi wait error in masterrecint')
  return
end subroutine masterrecint

!===========================================================
!             Logical
!===========================================================
subroutine mastersenlog(to,tag,buf,len)
  use chm_kinds
  use mpi                          
  use repd_ensemble
  implicit none
  
  integer,intent(inout) :: to, tag, len
  logical,intent(inout),dimension(len) ::  buf
  
  integer ierr,req,m_stat(mpi_status_size)
  if(len <= 0) return

  ! print '("senint ",i3,a,i3,a,i3,i10,i8)',ensmasternod," send to ",to," tag ",tag,buf(1),len

  call mpi_isend(buf,len,mpi_logical,to,tag,comm_master,req,ierr)
  if(ierr /= mpi_success) call wrndie(-4,'repd_ensemble.src<mastersenint>', &
       'mpi isend error in mastersenint')
  call mpi_wait(req,m_stat,ierr)
  ! print '("senint ",i3,a,i3,a,i3,"  success ")',ensmasternod," send to ",to," tag ",tag
  if(ierr /= mpi_success) call wrndie(-4,'repd_ensemble.src<mastersenint>', &
       'mpi wait error in mastersenint')
  return
end subroutine mastersenlog

subroutine masterreclog(from,tag,buf,len)
  !-----------------------------------------------------------------------
  !
  !     general receive routine. its purpose is to ease replacement
  !     for various platforms (intel, cm-5, sp1, tcp/ip, ...)
  !
  use chm_kinds
  use parallel,only:comm_charmm    
  use mpi                          
  use repd_ensemble
  implicit none

  integer,intent(in) :: from, tag, len
  logical,intent(inout),dimension(len) :: buf

  integer m_stat(mpi_status_size),ierr,req
  if(len <= 0) return

  ! print '("recint ",i3,a,i3,a,i3,a)',ensmasternod," rec frm ",from," tag ",tag,"posted"

  call mpi_irecv(buf,len,mpi_logical,from,tag,comm_master,req, ierr)
  if(ierr /= mpi_success)  &
       call wrndie(-4,'repd-ensemble.src<masterrecint>','mpi irecv error in masterrecint')
  ! print '("recint ",i3,a,i3,a,i3,i10,i8)',ensmasternod," waiting   ",from," tag ",tag
  call mpi_wait(req,m_stat,ierr)
  ! print '("recint ",i3,a,i3,a,i3,i10,i8)',ensmasternod," done    ",from," tag ",tag,buf(1),len
  if(ierr /= mpi_success)  &
       call wrndie(-4,'repd-ensemble.src<masterrecint>','mpi wait error in masterrecint')
  return
end subroutine masterreclog

!===========================================================
!             Real
!===========================================================
subroutine mastersenreal(to,tag,buf,len)
  use chm_kinds
  use mpi                          
  use repd_ensemble
  implicit none
  integer,intent(inout) :: to, tag, len
  real(chm_real),intent(inout),dimension(len) ::  buf
  integer ierr,req,m_stat(mpi_status_size)
  if(len <= 0) return
  ! print '("senreal ",i3,a,i3,a,i3,e15.3,i8)',ensmasternod," send to ",to," tag ",tag,buf(1),len

  call mpi_send(buf,len,mpi_real,to,tag,comm_master,ierr)
  !--   call mpi_isend(buf,len,mpi_real8,to,tag,comm_master,req,ierr)
  !-- 
  !--   if(ierr /= mpi_success) call wrndie(-4,'repd_ensemble.src<mastersenint>', &
  !--        'mpi isend error in mastersenint')
  !--   call mpi_wait(req,m_stat,ierr)
  !-- print '("senreal ",i3,a,i3,a,i3,"  success ")',ensmasternod," send to ",to," tag ",tag
  if(ierr /= mpi_success) call wrndie(-4,'repd_ensemble.src<mastersenint>', &
       'mpi wait error in mastersenint')
  return
end subroutine mastersenreal

subroutine masterrecreal(from,tag,buf,len)
  use chm_kinds
  use mpi                          
  use repd_ensemble
  implicit none

  integer,intent(inout) :: from, tag, len
  real(chm_real),intent(inout),dimension(len) :: buf
  integer m_stat(mpi_status_size),ierr,req

  if(len <= 0) return
  ! print '("recreal ",i3,a,i3,a,i3,a)',ensmasternod," rec frm ",from," tag ",tag,"posted"

  ! RLH call mpi_recv(buf,len,mpi_real,from,tag,comm_master, ierr)
  call mpi_recv(buf,len,mpi_real,from,tag,comm_master,m_stat,ierr)
  !--   call mpi_irecv(buf,len,mpi_real8,from,tag,comm_master,req, ierr)
  !-- 
  !--   if(ierr /= mpi_success)  &
  !--        call wrndie(-4,'repd-ensemble.src<masterrecint>','mpi irecv error in masterrecint')
  !-- print '("recreal ",i3,a,i3,a,i3,e15.3,i8)',ensmasternod," waiting   ",from," tag ",tag
  !--   call mpi_wait(req,m_stat,ierr)
  !-- print '("recreal ",i3,a,i3,a,i3,e15.5,i8)',ensmasternod," done    ",from," tag ",tag,buf(1),len
  if(ierr /= mpi_success)  &
       call wrndie(-4,'repd-ensemble.src<masterrecint>','mpi wait error in masterrecint')
  return
end subroutine masterrecreal

#endif /* KEY_REPDSTR2 */
