module repd_ensemble
  use chm_kinds
  use dimens_fcm

! General ENSEMBLE variables
! --------------------------
! WHOIAM: number of each node according to MPI
! NENSEM: number of copies started by MPI

#if KEY_REPDSTR2==1

  integer,save :: outu_repdstr
  integer :: whoiam,nensem
  integer,parameter :: maxensemble_layers=1
  integer comm_master
  integer prnlevr,wrnlevr,iolevr,plnod0r
  logical lmasternode,mastermaster,lslavenode
  integer,save :: old_mynod,old_comm_charmm,old_numnod
  integer,save :: ensmasternod
  logical,save :: repd_ensemble_verbose
  integer,save :: ensbfl=120
  character(len=80),dimension(:),allocatable :: chmout

contains

  subroutine rensprint(a1,a2)
    use parallel,only:mynod
    character(len=*),intent(in) :: a1,a2
    character(len=80) :: fmt
    character(len=80) :: filnam
    if(.not. repd_ensemble_verbose) return
    write(filnam,'("prn",i3.3)')old_mynod
    open(action="write",file=filnam,unit=80,position="APPEND")

    fmt='(a," e",i2," n",i2," w",i2," >>> ",a)'
    write(80,fmt) &
!    print fmt, &
         a1,whoiam,mynod,old_mynod,a2
    close(unit=80)
    return
  end subroutine rensprint
#endif 
end module repd_ensemble

