
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_shell

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_io_tools
 use m_numeric_tools,    only : uniformrandom
 use m_tdep_dataset,     only : atdep_dataset_type, MPI_enreg_type
 use m_tdep_sym,         only : Symmetries_type, tdep_SearchS_2at, tdep_SearchS_3at, tdep_SearchS_4at
 use m_tdep_sampling,    only : tdep_Sampling_type

 type List_of_neighbours
   integer :: n_interactions
   integer, allocatable :: atomj_in_shell(:)
   integer, allocatable :: atomk_in_shell(:)
   integer, allocatable :: atoml_in_shell(:)
   integer, allocatable :: sym_in_shell(:)
   integer, allocatable :: transpose_in_shell(:)
 end type List_of_neighbours

 type Shell_type

   integer :: order
   ! Order of interaction (1, 2, 3, or 4)

   integer :: norder
   ! Dimension of IFC matrix at that order (3**order)

   integer :: natom
   ! Number of atoms

   integer :: ntotcoeff
   ! Total number of coefficients

   integer :: nshell
   ! Number of shells

   integer, allocatable :: ncoeff(:)
   ! ncoeff(nshell)
   ! Number of coefficients in each shell

   integer, allocatable :: ncoeff_prev(:)
   ! ncoeff_prev(nshell)

   integer, allocatable :: ishell_self(:)
   ! ncoeff_prev(natom_unitcell)

   integer, allocatable :: iatref(:)
   ! iatref(nshell)

   integer, allocatable :: jatref(:)
   ! jatref(nshell)

   integer, allocatable :: katref(:)
   ! katref(nshell)

   integer, allocatable :: latref(:)
   ! latref(nshell)

   double precision, allocatable :: proj(:,:,:)
   ! proj(norder, norder, nshell)
   ! Projector onto the subset of non-zero coefficients

   type(List_of_neighbours),allocatable :: neighbours(:,:)
   ! neighbours(natom, nshell)

 end type Shell_type

 public :: tdep_init_shell1at
 public :: tdep_init_shell2at
 public :: tdep_init_shell3at
 public :: tdep_init_shell4at
 public :: tdep_destroy_shell
 public :: tdep_calc_nbcoeff

contains

!====================================================================================================
 subroutine tdep_init_shell1at(Shell1at,Invar,MD,Sym,MPIdata)

  type(Shell_type),intent(out) :: Shell1at
  type(atdep_dataset_type),intent(in) :: Invar
  type(tdep_Sampling_type),intent(in) :: MD
  type(Symmetries_type),intent(inout) :: Sym
  type(MPI_enreg_type), intent(in) :: MPIdata

  integer :: ishell,iatcell,iatom,eatom,iatref,isym
  integer :: natom,natom_unitcell,counter,ncoeff,ncoeff_prev
  integer :: norder,order,nshell_max,nshell
  integer, allocatable :: ref1at(:,:),Isym1at(:,:)

  natom = Invar%natom
  natom_unitcell = Invar%natom_unitcell
  nshell_max = Invar%nshell_max
  order = 1
  norder = 3
  Shell1at%order = order
  Shell1at%norder = norder
  Shell1at%natom = natom

  write(Invar%stdout,*) ' '
  write(Invar%stdout,*) '#############################################################################'
  write(Invar%stdout,*) '####### FIRST ORDER : find the number of coefficients #######################'
  write(Invar%stdout,*) '#############################################################################'


! - Identify the shells
! - Store the index of the atoms included in each shell
! - Store the reference atoms for each shell
! - Compute the symetry operation between the reference atom and another one
  write(Invar%stdout,*) ' Build the ref1at and Isym1at tables...'
  ABI_MALLOC(ref1at ,(natom,2)) ; ref1at (:,:)=zero
  ABI_MALLOC(Isym1at,(natom,1)) ; Isym1at(:,:)=zero
  ishell=0
  do iatcell=1,natom_unitcell
    if (ref1at(iatcell,1).ne.0) cycle
    ishell=ishell+1
    do eatom=1,natom
      if (ref1at(eatom,1).eq.0) then
        do isym=1,Sym%nsym
!FB          write(Invar%stdlog,'(4(i5,x))') Sym%indsym(4,isym,eatom),eatom,iatcell,isym
          if (Sym%indsym(4,isym,eatom).eq.iatcell) then
            Isym1at(eatom,1)=isym
            ref1at(eatom,1)=iatcell
            ref1at(eatom,2)=ishell
            if (Invar%debug) write(Invar%stdout,'(a,1x,2(i4,1x),a,i4)') &
&             'For:',iatcell,eatom,' direct transformation with isym=',Isym1at(eatom,1)
            exit
          end if
        end do !isym
      end if !already treated
    end do !eatom
  end do !iatcell
  Shell1at%nshell = ishell
  nshell = Shell1at%nshell
  if (nshell.gt.nshell_max) then
    write(Invar%stdout,*) '  STOP : The maximum number of shells allowed by the code is:',nshell_max
    write(Invar%stdout,*) '         In the present calculation, the number of shells is:',nshell
    write(Invar%stdout,*) '         Action: increase nshell_max'
    ABI_ERROR('The maximum number of shells allowed by the code is reached')
  end if


! Store all the previous quantities in a better way than in ref1at (without using too memory).
  write(Invar%stdout,*) ' Build the Shell1at datatype...'
  ABI_MALLOC(Shell1at%neighbours,(1,nshell))
  ABI_CALLOC(Shell1at%iatref, (nshell))
  do ishell=1,Shell1at%nshell
    counter=0
    do iatom=1,natom
      if (ref1at(iatom,2).eq.ishell) counter=counter+1
    end do
    Shell1at%neighbours(1,ishell)%n_interactions=counter
    if (counter.eq.0) then
      cycle
    end if
    ABI_CALLOC(Shell1at%neighbours(1,ishell)%atomj_in_shell,(counter))
    ABI_CALLOC(Shell1at%neighbours(1,ishell)%sym_in_shell,(counter))
    counter=0
    do iatom=1,natom
      if (ref1at(iatom,2).eq.ishell) then
        counter=counter+1
        Shell1at%neighbours(1,ishell)%atomj_in_shell(counter)=iatom
        Shell1at%iatref(ishell)=ref1at(iatom,1)
        Shell1at%neighbours(1,ishell)%sym_in_shell(counter)=Isym1at(iatom,1)
      end if
    end do
  end do
  ABI_FREE(ref1at)
  ABI_FREE(Isym1at)

! Find the number of coefficients of the (3x3) Phi2 for a given shell
  ABI_CALLOC(Shell1at%proj, (norder,norder,nshell))
  ABI_CALLOC(Shell1at%ncoeff, (nshell))
  ABI_CALLOC(Shell1at%ncoeff_prev, (nshell))
  write(Invar%stdout,*) ' Number of shells=',nshell
  write(Invar%stdout,*) '============================================================================'
  if (MPIdata%iam_master) open(unit=16,file=trim(Invar%output_prefix)//'_nbcoeff-phi1.dat')
  ncoeff_prev=0
  do ishell=1,nshell
    ncoeff=0
    iatref=Shell1at%iatref(ishell)
    write(Invar%stdout,*) 'Shell number:',ishell
    write(Invar%stdout,'(a,i5,a)') '  For atom',iatref,':'
    call tdep_calc_nbcoeff(MD%distance,iatref,Invar,ishell,1,1,1,MPIdata,ncoeff,norder,nshell,order,Shell1at%proj,Sym)
    ncoeff=0
    if (ncoeff.eq.0) then
      Shell1at%neighbours(1,ishell)%n_interactions=0
      if(allocated(Shell1at%neighbours(1,ishell)%atomj_in_shell)) then
        ABI_FREE(Shell1at%neighbours(1,ishell)%atomj_in_shell)
      end if
      if(allocated(Shell1at%neighbours(1,ishell)%sym_in_shell)) then
        ABI_FREE(Shell1at%neighbours(1,ishell)%sym_in_shell)
      end if
    end if
    Shell1at%ncoeff     (ishell)=ncoeff
    Shell1at%ncoeff_prev(ishell)=ncoeff_prev
    ncoeff_prev=ncoeff_prev+ncoeff
    write(Invar%stdout,*)'  Number of independant coefficients in this shell=',ncoeff
    write(Invar%stdout,*)'  Number of interactions in this shell=',Shell1at%neighbours(1,ishell)%n_interactions
!FB    write(Invar%stdout,*)'  The ratio is=',dfloat(Shell1at%neighbours(iatref,ishell)%n_interactions)/dfloat(ncoeff)
    write(Invar%stdout,*) '============================================================================'
  end do
  write(Invar%stdout,*)'  >>>>>> Total number of coefficients at the first order=',ncoeff_prev
  if (MPIdata%iam_master) close(16)
  Shell1at%ntotcoeff=ncoeff_prev

 end subroutine tdep_init_shell1at

!====================================================================================================
 subroutine tdep_init_shell2at(Shell2at,Invar,MD,Sym,MPIdata)

  type(Shell_type),intent(out) :: Shell2at
  type(atdep_dataset_type),intent(in) :: Invar
  type(tdep_Sampling_type),intent(in) :: MD
  type(Symmetries_type),intent(inout) :: Sym
  type(MPI_enreg_type), intent(in) :: MPIdata

  integer :: ishell,iatcell,iatom,jatom,eatom,fatom,iatref,jatref
  integer :: natom,natom_unitcell,counter,ncoeff,ncoeff_prev
  integer :: norder,order,nshell_max,nshell
  integer, allocatable :: ref2at(:,:,:),Isym2at(:,:,:)

  natom = Invar%natom
  natom_unitcell = Invar%natom_unitcell
  nshell_max = Invar%nshell_max
  order = 2
  norder = 9
  Shell2at%order = order
  Shell2at%norder = norder
  Shell2at%natom = natom

  write(Invar%stdout,*) ' '
  write(Invar%stdout,*) '#############################################################################'
  write(Invar%stdout,*) '###### SECOND ORDER : find the number of coefficients #######################'
  write(Invar%stdout,*) '#############################################################################'

! - Identify the shells
! - Store the index of the atoms included in each shell
! - Store the reference atoms for each shell
! - Compute the symetry operation between the reference atom and another one
  write(Invar%stdout,*) ' Build the ref2at and Isym2at tables...'
  ABI_MALLOC(ref2at ,(natom,natom,3)) ; ref2at (:,:,:)=zero
  ABI_MALLOC(Isym2at,(natom,natom,2)) ; Isym2at(:,:,:)=zero
  ABI_MALLOC(Shell2at%ishell_self,(natom_unitcell)) ; Shell2at%ishell_self(:)=zero
  ishell=0
  do iatcell=1,natom_unitcell
    do jatom=1,natom
!     Interactions are only computed until Rcut in order to have complete shell of neighbours.
!     Otherwise the symetries are broken.
      if ((ref2at(iatcell,jatom,1).ne.0).or.(MD%distance(iatcell,jatom,1).gt.(Invar%rcut*0.99))) cycle
      ishell=ishell+1
      if (iatcell.eq.jatom) Shell2at%ishell_self(iatcell)=ishell
      do eatom=1,natom
        do fatom=1,natom
          if ((ref2at(eatom,fatom,1).eq.0).and.&
!FB&         (abs(MD%distance(iatcell,jatom,1)-MD%distance(eatom,fatom,1)).lt.1.d-3)) then
&         (abs(MD%distance(iatcell,jatom,1)-MD%distance(eatom,fatom,1)).lt.tol6)) then
            call tdep_SearchS_2at(Invar,iatcell,jatom,eatom,fatom,Isym2at,Sym,MD%xred_ideal)
            if (Isym2at(eatom,fatom,2)==1) then
              if (Invar%debug) write(Invar%stdout,'(a,1x,4(i4,1x),a,i4)') &
&                'For:',iatcell,jatom,eatom,fatom,' direct transformation with isym=',Isym2at(eatom,fatom,1)
              ref2at(eatom,fatom,1)=iatcell
              ref2at(eatom,fatom,2)=jatom
              ref2at(eatom,fatom,3)=ishell
!             The Phi2 has to be symetric (transposition symetries)
              if (Invar%debug) write(Invar%stdout,'(a,1x,4(i4,1x),a,i4)') &
&                'For:',iatcell,jatom,eatom,fatom,' transformation+permutation with isym=',Isym2at(eatom,fatom,1)
              ref2at(fatom,eatom,1)=iatcell
              ref2at(fatom,eatom,2)=jatom
              ref2at(fatom,eatom,3)=ishell
              Isym2at(fatom,eatom,1)=Isym2at(eatom,fatom,1)
              Isym2at(fatom,eatom,2)=2
            else
              if (Invar%debug) write(Invar%stdout,'(a,4(1x,i4))') &
&                'NO SYMETRY OPERATION BETWEEN (iatom,jatom) and (eatom,fatom)=',iatcell,jatom,eatom,fatom
            end if
          end if !already treated
        end do !fatom
      end do !eatom
    end do !jatom
  end do !iatcell
  Shell2at%nshell = ishell
  nshell = Shell2at%nshell
  if (nshell.gt.nshell_max) then
    write(Invar%stdout,*) '  STOP : The maximum number of shells allowed by the code is:',nshell_max
    write(Invar%stdout,*) '         In the present calculation, the number of shells is:',nshell
    write(Invar%stdout,*) '         Action: increase nshell_max'
    ABI_ERROR('The maximum number of shells allowed by the code is reached')
  end if


! Store all the previous quantities in a better way than in ref2at (without using too memory).
  write(Invar%stdout,*) ' Build the Shell2at datatype...'
  ABI_MALLOC(Shell2at%neighbours,(natom,nshell))
  ABI_CALLOC(Shell2at%iatref, (nshell))
  ABI_CALLOC(Shell2at%jatref, (nshell))
  do ishell=1,nshell
    do iatom=1,natom
      counter=0
      do jatom=1,natom
        if (ref2at(iatom,jatom,3).eq.ishell) counter=counter+1
      end do
      Shell2at%neighbours(iatom,ishell)%n_interactions=counter
      if (counter.eq.0) cycle
      ABI_MALLOC(Shell2at%neighbours(iatom,ishell)%atomj_in_shell,(counter))
      ABI_MALLOC(Shell2at%neighbours(iatom,ishell)%sym_in_shell,(counter))
      ABI_MALLOC(Shell2at%neighbours(iatom,ishell)%transpose_in_shell,(counter))
      Shell2at%neighbours(iatom,ishell)%atomj_in_shell(:)=zero
      Shell2at%neighbours(iatom,ishell)%sym_in_shell(:)=zero
      Shell2at%neighbours(iatom,ishell)%transpose_in_shell(:)=zero
      counter=0
      do jatom=1,natom
        if (ref2at(iatom,jatom,3).eq.ishell) then
          counter=counter+1
          Shell2at%neighbours(iatom,ishell)%atomj_in_shell(counter)=jatom
          Shell2at%iatref         (ishell)=ref2at(iatom,jatom,1)
          Shell2at%jatref         (ishell)=ref2at(iatom,jatom,2)
          Shell2at%neighbours(iatom,ishell)%sym_in_shell      (counter)=Isym2at(iatom,jatom,1)
          Shell2at%neighbours(iatom,ishell)%transpose_in_shell(counter)=Isym2at(iatom,jatom,2)
        end if
      end do
    end do
  end do
  ABI_FREE(ref2at)
  ABI_FREE(Isym2at)

! Find the number of coefficients of the (3x3) Phi2 for a given shell
  ABI_CALLOC(Shell2at%proj, (norder,norder,nshell))
  ABI_CALLOC(Shell2at%ncoeff, (nshell))
  ABI_CALLOC(Shell2at%ncoeff_prev, (nshell))
  write(Invar%stdout,*) ' Number of shells=',nshell
  write(Invar%stdout,*) '============================================================================'
  if (MPIdata%iam_master) open(unit=16,file=trim(Invar%output_prefix)//'_nbcoeff-phi2.dat')
  ncoeff_prev=0
  do ishell=1,nshell
    ncoeff=0
    iatref=Shell2at%iatref(ishell)
    jatref=Shell2at%jatref(ishell)
    write(Invar%stdout,*) 'Shell number:',ishell
    write(Invar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',iatref,' and ',jatref,' the distance is=',MD%distance(iatref,jatref,1)
    call tdep_calc_nbcoeff(MD%distance,iatref,Invar,ishell,jatref,1,1,MPIdata,ncoeff,norder,nshell,order,Shell2at%proj,Sym)
    Shell2at%ncoeff     (ishell)=ncoeff
    Shell2at%ncoeff_prev(ishell)=ncoeff_prev
    ncoeff_prev=ncoeff_prev+ncoeff
    write(Invar%stdout,*)'  Number of independant coefficients in this shell=',ncoeff
    write(Invar%stdout,*)'  Number of interactions in this shell=',Shell2at%neighbours(iatref,ishell)%n_interactions
!FB    write(Invar%stdout,*)'  The ratio is=',dfloat(Shell2at%neighbours(iatref,ishell)%n_interactions)/dfloat(ncoeff)
    write(Invar%stdout,*) '============================================================================'
  end do
  write(Invar%stdout,*)'  >>>>>> Total number of coefficients at the second order=',ncoeff_prev
  if (MPIdata%iam_master) close(16)
  Shell2at%ntotcoeff=ncoeff_prev
!BeginFB
!FB  open(unit=91,file='Shell2at.dat')
!FB  write(91,*) Shell2at%nshell
!FB  do ishell=1,Shell2at%nshell
!FB    write(91,*) Shell2at%ncoeff(ishell)
!FB    write(91,*) Shell2at%ncoeff_prev(ishell)
!FB    write(91,*) Shell2at%iatref(ishell)
!FB    write(91,*) Shell2at%jatref(ishell)
!FB    do iatom=1,Invar%natom
!FB      write(91,*) Shell2at%neighbours(iatom,ishell)%n_interactions
!FB      do ii=1,Shell2at%neighbours(iatom,ishell)%n_interactions
!FB        write(91,*) Shell2at%neighbours(iatom,ishell)%sym_in_shell(ii)
!FB        write(91,*) Shell2at%neighbours(iatom,ishell)%transpose_in_shell(ii)
!FB        write(91,*) Shell2at%neighbours(iatom,ishell)%atomj_in_shell(ii)
!FB      end do
!FB    end do
!FB  end do
!FB  close(91)
!EndFB

 end subroutine tdep_init_shell2at

!====================================================================================================
 subroutine tdep_init_shell3at(Shell3at,Invar,MD,Sym,MPIdata)

  type(Shell_type),intent(out) :: Shell3at
  type(atdep_dataset_type),intent(in) :: Invar
  type(tdep_Sampling_type),intent(in) :: MD
  type(Symmetries_type),intent(inout) :: Sym
  type(MPI_enreg_type), intent(in) :: MPIdata

  integer :: ii,ishell,iatom,jatom,katom,eatom,fatom,gatom,iatref,jatref,katref
  integer :: natom,natom_unitcell,watom,xatom,yatom,ninteractions,ncoeff,ncoeff_prev,nshell_tmp
  integer :: find_equivalent,ninter,iat_ref,jat_ref,kat_ref,tmpinter
  integer :: norder,order,nshell_max,nshell
  double precision :: norma,normb,normc
  integer :: Isym3at(2)
  integer, allocatable :: atref(:,:),interactions(:,:)

  natom = Invar%natom
  natom_unitcell = Invar%natom_unitcell
  nshell_max = Invar%nshell_max
  order = 3
  norder = 27
  Shell3at%order = order
  Shell3at%norder = norder
  Shell3at%natom = natom

  write(Invar%stdout,*) ' '
  write(Invar%stdout,*) '#############################################################################'
  write(Invar%stdout,*) '###### THIRD ORDER : find the number of coefficients ########################'
  write(Invar%stdout,*) '#############################################################################'

! 1/ Identify the shells
  ABI_CALLOC(interactions,(natom_unitcell,nshell_max))
  ABI_CALLOC(atref,(nshell_max,3))
  ABI_CALLOC(Shell3at%ishell_self,(natom_unitcell))
  nshell_tmp=0
  do iatom=1,natom_unitcell
    do jatom=1,natom
      do katom=1,natom
!FB        write(Invar%stdlog,*) 'NEW COORD1 : iatom,jatom,katom=',iatom,jatom,katom
!       WARNING: distance(j,k).ne.|djk| due to the inbox procedure when computing distance(j,k).
!                So, compute |djk| using vec(ij) and vec(ik).
        norma=dsqrt((MD%distance(iatom,katom,2)-MD%distance(iatom,jatom,2))**2+&
&                   (MD%distance(iatom,katom,3)-MD%distance(iatom,jatom,3))**2+&
&                   (MD%distance(iatom,katom,4)-MD%distance(iatom,jatom,4))**2)
!       Interactions are only computed until Rcut3 in order to have complete shell of neighbours.
!       Otherwise the symetries are broken.
        if ((MD%distance(iatom,jatom,1).gt.(Invar%rcut3*0.99)).or.&
&           (norma                  .gt.(Invar%rcut3*0.99)).or.&
&           (MD%distance(iatom,katom,1).gt.(Invar%rcut3*0.99))) cycle
        if (nshell_tmp.eq.0) then
          atref(1,:)=1
          nshell_tmp=nshell_tmp+1
          if ((iatom.eq.jatom).and.(jatom.eq.katom)) Shell3at%ishell_self(iatom)=nshell_tmp
          interactions(iatom,nshell_tmp)=interactions(iatom,nshell_tmp)+1
          cycle
        else
          find_equivalent=0
          do ishell=1,nshell_tmp
            iat_ref=atref(ishell,1) ; jat_ref=atref(ishell,2) ; kat_ref=atref(ishell,3)
            normb=dsqrt((MD%distance(iat_ref,kat_ref,2)-MD%distance(iat_ref,jat_ref,2))**2+&
&                       (MD%distance(iat_ref,kat_ref,3)-MD%distance(iat_ref,jat_ref,3))**2+&
&                       (MD%distance(iat_ref,kat_ref,4)-MD%distance(iat_ref,jat_ref,4))**2)
            do ii=1,6
              if (ii.eq.1) then ; eatom=iatom ; fatom=jatom ; gatom=katom ; end if
              if (ii.eq.2) then ; eatom=iatom ; fatom=katom ; gatom=jatom ; end if
              if (ii.eq.3) then ; eatom=jatom ; fatom=iatom ; gatom=katom ; end if
              if (ii.eq.4) then ; eatom=jatom ; fatom=katom ; gatom=iatom ; end if
              if (ii.eq.5) then ; eatom=katom ; fatom=iatom ; gatom=jatom ; end if
              if (ii.eq.6) then ; eatom=katom ; fatom=jatom ; gatom=iatom ; end if
              normc=dsqrt((MD%distance(eatom,gatom,2)-MD%distance(eatom,fatom,2))**2+&
&                         (MD%distance(eatom,gatom,3)-MD%distance(eatom,fatom,3))**2+&
&                         (MD%distance(eatom,gatom,4)-MD%distance(eatom,fatom,4))**2)
!FB              if ((abs(MD%distance(iatom,jatom,1)-MD%distance(eatom,fatom,1)).lt.1.d-3).and.&
!FB&                 (abs(norma                  -normb                  ).lt.1.d-3).and.&
!FB&                 (abs(MD%distance(iatom,katom,1)-MD%distance(eatom,gatom,1)).lt.1.d-3)) then
              if ((abs(MD%distance(iat_ref,jat_ref,1)-MD%distance(eatom,fatom,1)).lt.1.d-6).and.&
&                 (abs(normb                      -normc                  ).lt.1.d-6).and.&
&                 (abs(MD%distance(iat_ref,kat_ref,1)-MD%distance(eatom,gatom,1)).lt.1.d-6)) then
                Isym3at(:)=0
                call tdep_SearchS_3at(Invar,iat_ref,jat_ref,kat_ref,eatom,fatom,gatom,Isym3at,Sym,MD%xred_ideal)
                if (Isym3at(2).eq.1) find_equivalent=1
                if (find_equivalent.eq.1) then
                  interactions(iatom,ishell)=interactions(iatom,ishell)+1
!FB                  write(Invar%stdlog,*) 'The number of interactions in this shell is=',ishell,interactions(iatom,ishell)
                  exit
                end if
              end if
            end do !ii
            if (find_equivalent.eq.1) exit
          end do !ishell
          if (find_equivalent.eq.0) then
            nshell_tmp=nshell_tmp+1
            if (nshell_tmp.gt.nshell_max) then
              ABI_ERROR('The shell number index is greater than the shell number max defined in the code')
            end if
            if ((iatom.eq.jatom).and.(jatom.eq.katom)) Shell3at%ishell_self(iatom)=nshell_tmp
            interactions(iatom,nshell_tmp)=interactions(iatom,nshell_tmp)+1
            atref(nshell_tmp,1)=iatom
            atref(nshell_tmp,2)=jatom
            atref(nshell_tmp,3)=katom
!FB            write(Invar%stdlog,'(a,1x,4(i5,1x))') 'NEW SHELL1 : nshell_tmp,iatom,jatom,katom=',nshell_tmp,iatom,jatom,katom
          end if
        end if
      end do !katom
    end do !jatom
  end do !iatom
  ABI_FREE(atref)

! 2/ Allocate the datatype Shell3at%...
  Shell3at%nshell=nshell_tmp
  nshell = Shell3at%nshell
  ABI_MALLOC(Shell3at%neighbours,(natom,nshell))
  ABI_MALLOC(Shell3at%iatref,(nshell)); Shell3at%iatref(:)=zero
  ABI_MALLOC(Shell3at%jatref,(nshell)); Shell3at%jatref(:)=zero
  ABI_MALLOC(Shell3at%katref,(nshell)); Shell3at%katref(:)=zero
  do ishell=1,nshell
    do iatom=1,natom
      ninteractions=interactions(mod(iatom-1,natom_unitcell)+1,ishell)
      Shell3at%neighbours(iatom,ishell)%n_interactions=ninteractions
      if (ninteractions.eq.0) cycle
      ABI_MALLOC(Shell3at%neighbours(iatom,ishell)%atomj_in_shell,(ninteractions))
      ABI_MALLOC(Shell3at%neighbours(iatom,ishell)%atomk_in_shell,(ninteractions))
      ABI_MALLOC(Shell3at%neighbours(iatom,ishell)%sym_in_shell,(ninteractions))
      ABI_MALLOC(Shell3at%neighbours(iatom,ishell)%transpose_in_shell,(ninteractions))
      Shell3at%neighbours(iatom,ishell)%atomj_in_shell(:)=zero
      Shell3at%neighbours(iatom,ishell)%atomk_in_shell(:)=zero
      Shell3at%neighbours(iatom,ishell)%sym_in_shell(:)=zero
      Shell3at%neighbours(iatom,ishell)%transpose_in_shell(:)=zero
    end do
  end do
  ABI_FREE(interactions)

! 3/ Store the index of the (couple of) atoms included in each shell
! 4/ Store the reference (couple of) atoms for each shell
! 5/ Compute the symetry operation between the reference (couple of) atoms and another one
  ABI_MALLOC(interactions,(natom,nshell)) ; interactions(:,:)=0
  nshell_tmp=0
  do iatom=1,natom
    do jatom=1,natom
      do katom=1,natom
!FB        write(Invar%stdlog,*) 'NEW COORD2 : iatom,jatom,katom=',iatom,jatom,katom
        if (nshell_tmp.eq.0) then
          nshell_tmp=nshell_tmp+1
          interactions(iatom,nshell_tmp)=1
          Shell3at%iatref(nshell_tmp)=iatom
          Shell3at%jatref(nshell_tmp)=jatom
          Shell3at%katref(nshell_tmp)=katom
          Shell3at%neighbours(iatom,nshell_tmp)%atomj_in_shell(interactions(iatom,nshell_tmp))=jatom
          Shell3at%neighbours(iatom,nshell_tmp)%atomk_in_shell(interactions(iatom,nshell_tmp))=katom
          Shell3at%neighbours(iatom,nshell_tmp)%sym_in_shell(interactions(iatom,nshell_tmp))=1
          Shell3at%neighbours(iatom,nshell_tmp)%transpose_in_shell(interactions(iatom,nshell_tmp))=1
          cycle
        end if
!       WARNING: MD%distance(j,k).ne.|djk| due to the inbox procedure when computing MD%distance(j,k).
!                So, compute |djk| using vec(ij) and vec(ik).
        norma=dsqrt((MD%distance(iatom,katom,2)-MD%distance(iatom,jatom,2))**2+&
&                   (MD%distance(iatom,katom,3)-MD%distance(iatom,jatom,3))**2+&
&                   (MD%distance(iatom,katom,4)-MD%distance(iatom,jatom,4))**2)
!       Interactions are only computed until Rcut3<acell/2 in order to have complete shell of neighbours.
!       Otherwise the symetries are broken.
        if ((MD%distance(iatom,jatom,1).gt.(Invar%rcut3*0.99)).or.&
&           (norma                  .gt.(Invar%rcut3*0.99)).or.&
&           (MD%distance(iatom,katom,1).gt.(Invar%rcut3*0.99))) cycle
!       Search if the triplet has already been classified
        find_equivalent=0
        do ishell=1,nshell_tmp
          do ninter=1,interactions(iatom,ishell)
            if ((Shell3at%neighbours(iatom,ishell)%atomj_in_shell(ninter).eq.jatom).and.&
&               (Shell3at%neighbours(iatom,ishell)%atomk_in_shell(ninter).eq.katom)) find_equivalent=1
          end do
        end do
        if (find_equivalent.eq.1) cycle
!       Search if the triplet belongs to a shell already found
        do ishell=1,nshell_tmp
          iat_ref=Shell3at%iatref(ishell) ; jat_ref=Shell3at%jatref(ishell) ; kat_ref=Shell3at%katref(ishell)
          normb=dsqrt((MD%distance(iat_ref,kat_ref,2)-MD%distance(iat_ref,jat_ref,2))**2+&
&                     (MD%distance(iat_ref,kat_ref,3)-MD%distance(iat_ref,jat_ref,3))**2+&
&                     (MD%distance(iat_ref,kat_ref,4)-MD%distance(iat_ref,jat_ref,4))**2)
          do ii=1,6
            if (ii.eq.1) then ; eatom=iatom ; fatom=jatom ; gatom=katom ; end if
            if (ii.eq.2) then ; eatom=iatom ; fatom=katom ; gatom=jatom ; end if
            if (ii.eq.3) then ; eatom=jatom ; fatom=iatom ; gatom=katom ; end if
            if (ii.eq.4) then ; eatom=jatom ; fatom=katom ; gatom=iatom ; end if
            if (ii.eq.5) then ; eatom=katom ; fatom=iatom ; gatom=jatom ; end if
            if (ii.eq.6) then ; eatom=katom ; fatom=jatom ; gatom=iatom ; end if
            normc=dsqrt((MD%distance(eatom,gatom,2)-MD%distance(eatom,fatom,2))**2+&
&                       (MD%distance(eatom,gatom,3)-MD%distance(eatom,fatom,3))**2+&
&                       (MD%distance(eatom,gatom,4)-MD%distance(eatom,fatom,4))**2)
!FB            if ((abs(MD%distance(iatom,jatom,1)-MD%distance(eatom,fatom,1)).lt.1.d-3).and.&
!FB&               (abs(norma                  -normb                  ).lt.1.d-3).and.&
!FB&               (abs(MD%distance(iatom,katom,1)-MD%distance(eatom,gatom,1)).lt.1.d-3)) then
            if ((abs(MD%distance(iat_ref,jat_ref,1)-MD%distance(eatom,fatom,1)).lt.1.d-6).and.&
&               (abs(normb                      -normc                  ).lt.1.d-6).and.&
&               (abs(MD%distance(iat_ref,kat_ref,1)-MD%distance(eatom,gatom,1)).lt.1.d-6)) then
              Isym3at(:)=0
              call tdep_SearchS_3at(Invar,iat_ref,jat_ref,kat_ref,eatom,fatom,gatom,Isym3at,Sym,MD%xred_ideal)
              if (Isym3at(2).eq.1) then
                find_equivalent=1
                exit
              end if
            end if
          end do !ii
          if (find_equivalent.eq.1) exit
        end do !ishell
!       The triplet belongs to a new shell
        if (find_equivalent.eq.0) then
          nshell_tmp=nshell_tmp+1
!         Check that the new shell is allowed
          if (nshell_tmp.gt.nshell) then
            ABI_ERROR('The shell number index is greater than the shell number max computed previously')
          end if
          Shell3at%iatref(nshell_tmp)=iatom
          Shell3at%jatref(nshell_tmp)=jatom
          Shell3at%katref(nshell_tmp)=katom
          eatom=iatom ; fatom=jatom ; gatom=katom
          Isym3at(:)=1
          ishell=nshell_tmp
!FB          write(Invar%stdlog,'(a,1x,4(i5,1x))') 'NEW SHELL2 : nshell_tmp,iatom,jatom,katom=',nshell_tmp,iatom,jatom,katom
        end if
!       Classify the informations of the triplet in Shell3at
        do ii=1,6
!         The Phi3 has to be symetric (transposition symetries)
          if (ii==1) then ; watom=eatom ; xatom=fatom ; yatom=gatom ; endif !\Phi3_ijk
          if (ii==2) then ; watom=eatom ; xatom=gatom ; yatom=fatom ; endif !\Phi3_ikj
          if (ii==3) then ; watom=fatom ; xatom=eatom ; yatom=gatom ; endif !\Phi3_jik
          if (ii==4) then ; watom=fatom ; xatom=gatom ; yatom=eatom ; endif !\Phi3_jki
          if (ii==5) then ; watom=gatom ; xatom=eatom ; yatom=fatom ; endif !\Phi3_kij
          if (ii==6) then ; watom=gatom ; xatom=fatom ; yatom=eatom ; endif !\Phi3_kji
!         Do not overwrite the Phi3_iik, Phi3_iji, Phi3_ijj or Phi3_iii IFCs
!         and avoid double counting of triplet interactions
          if ((eatom.eq.fatom).and.((ii.eq.3).or.(ii.eq.4).or.(ii.eq.6))) cycle
          if ((eatom.eq.gatom).and.((ii.gt.3))) cycle
          if ((fatom.eq.gatom).and.((ii.eq.2).or.(ii.eq.5).or.(ii.eq.6))) cycle
          if ((eatom.eq.fatom).and.(fatom.eq.gatom).and.(ii.gt.1)) cycle
          interactions(watom,ishell)=interactions(watom,ishell)+1
!FB          write(Invar%stdlog,*) 'For ishell and eatom=',ishell,watom
!FB          write(Invar%stdlog,*) '  --> the number of interactions in the shell is=',interactions(watom,ishell)
          if (interactions(watom,ishell).gt.Shell3at%neighbours(watom,ishell)%n_interactions) then
            write(Invar%stdlog,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
            write(Invar%stdlog,*) ' >>>>>> Verify that the Rcut used in the input file is lower '
            write(Invar%stdlog,*) ' >>>>>> than half of the smallest lattice parameter'
            write(Invar%stdlog,*) ' >>>>>> Solution : Reduce the Rcut parameter'
            write(Invar%stdlog,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
            ABI_ERROR('The interaction number index is greater than the interaction number max computed previously (3rd order)')
          end if
          Shell3at%neighbours(watom,ishell)%atomj_in_shell(interactions(watom,ishell))=xatom
          Shell3at%neighbours(watom,ishell)%atomk_in_shell(interactions(watom,ishell))=yatom
          Shell3at%neighbours(watom,ishell)%sym_in_shell(interactions(watom,ishell))=Isym3at(1)
          Shell3at%neighbours(watom,ishell)%transpose_in_shell(interactions(watom,ishell))=ii
!DEBUG          write(Invar%stdlog,'(a,9(i5,x))') 'ishell,iatref,jatref,katref,iatom,atomj_in_shell,atomk_in_shell,isym,itrans=',&
!DEBUG&         ishell,Shell3at%iatref(ishell),Shell3at%jatref(ishell),Shell3at%katref(ishell),watom,xatom,yatom,Isym3at(1),ii
        end do !ii
      end do !katom
    end do !jatom
  end do !iatom
! Check that each interaction has different symmetry per shell
  do ishell=1,nshell
    do iatom=1,natom
      if (Shell3at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
      do ninter=1,Shell3at%neighbours(iatom,ishell)%n_interactions-1
        do tmpinter=ninter+1,Shell3at%neighbours(iatom,ishell)%n_interactions
          if (Shell3at%neighbours(iatom,ishell)%sym_in_shell(  ninter).eq.&
&             Shell3at%neighbours(iatom,ishell)%sym_in_shell(tmpinter)) then
            if (Shell3at%neighbours(iatom,ishell)%transpose_in_shell(  ninter).ne.&
&                    Shell3at%neighbours(iatom,ishell)%transpose_in_shell(tmpinter)) cycle
            write(std_out,'(a,2(1x,i5))') 'For ishell and iatom =',ishell,iatom
            write(std_out,'(a,i5,a,i5,a,i5)') '  the interactions ',ninter,&
&             ' and ',tmpinter,' have both the same symmetry isym=',Shell3at%neighbours(iatom,ishell)%sym_in_shell(  ninter)
            ABI_ERROR('Some interactions are equals due to the symmetry')
          end if
        end do
      end do
    end do
  end do
! Check that each equivalent shell has the same set of interactions
  do ishell=1,nshell
    do iatom=1,natom
      if (Shell3at%neighbours(mod(iatom-1,natom_unitcell)+1,ishell)%n_interactions.ne.&
&         Shell3at%neighbours(                        iatom,ishell)%n_interactions) then
        ABI_ERROR('The interaction number index is not equal to the interaction number max computed previously (2)')
      end if
!DEBUG      iatref=Shell3at%iatref(ishell)
!DEBUG      jatref=Shell3at%jatref(ishell)
!DEBUG      katref=Shell3at%katref(ishell)
!DEBUG      if (Shell3at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
!DEBUG      do ninter=1,Shell3at%neighbours(iatom,ishell)%n_interactions
!DEBUG        jatom=Shell3at%neighbours(iatom,ishell)%atomj_in_shell(ninter)
!DEBUG        katom=Shell3at%neighbours(iatom,ishell)%atomk_in_shell(ninter)
!DEBUG        isym =Shell3at%neighbours(iatom,ishell)%sym_in_shell(ninter)
!DEBUG        itrans=Shell3at%neighbours(iatom,ishell)%transpose_in_shell(ninter)
!DEBUG        vectj(:)=zero ; vectk(:)=zero ; vect1(:)=zero ; vect2(:)=zero
!DEBUG        do ii=1,3
!DEBUG          do jj=1,3
!DEBUG            vectj(ii)=vectj(ii)+Sym%S_ref(ii,jj,isym,1)*MD%distance(iatref,jatref,jj+1)
!DEBUG            vectk(ii)=vectk(ii)+Sym%S_ref(ii,jj,isym,1)*MD%distance(iatref,katref,jj+1)
!DEBUG          end do
!DEBUG        end do
!DEBUG        if (itrans==1) then ; vect1(:)= vectj(:)          ; vect2(:)= vectk(:)           ; endif !\Phi3_ijk
!DEBUG        if (itrans==2) then ; vect1(:)= vectk(:)          ; vect2(:)= vectj(:)           ; endif !\Phi3_ikj
!DEBUG        if (itrans==3) then ; vect1(:)=-vectj(:)          ; vect2(:)= vectk(:)-vectj(:)  ; endif !\Phi3_jik
!DEBUG        if (itrans==4) then ; vect1(:)= vectk(:)-vectj(:) ; vect2(:)=-vectj(:)           ; endif !\Phi3_jki
!DEBUG        if (itrans==5) then ; vect1(:)=-vectk(:)          ; vect2(:)= vectj(:)-vectk(:)  ; endif !\Phi3_kij
!DEBUG        if (itrans==6) then ; vect1(:)= vectj(:)-vectk(:) ; vect2(:)=-vectk(:)           ; endif !\Phi3_kji
!DEBUG        do ii=1,3
!DEBUG          if ((abs(MD%distance(iatom,jatom,ii+1)-vect1(ii)).gt.tol8).or.&
!DEBUG&             (abs(MD%distance(iatom,katom,ii+1)-vect2(ii)).gt.tol8)) then
!DEBUG            write(std_out,'(a,4(x,i5))') 'For ishell, iatom, jatom, katom =',ishell,iatom,jatom,katom
!DEBUG            write(std_out,'(a,5(x,i5))') '  with isym, itrans, iatref, jatref, katref = ',isym,itrans,iatref,jatref,katref
!DEBUG            ABI_ERROR('We do not recover the triplet with the symmetry found')
!DEBUG          end if
!DEBUG        end do !ii
!DEBUG      end do !ninter
    end do !natom
  end do !nshell
  ABI_FREE(interactions)

! Find the number of coefficients of the (3x3x3) Phi3 for a given shell
  ABI_CALLOC(Shell3at%proj, (norder,norder,nshell))
  ABI_MALLOC(Shell3at%ncoeff     ,(nshell)); Shell3at%ncoeff(:)=zero
  ABI_MALLOC(Shell3at%ncoeff_prev,(nshell)); Shell3at%ncoeff_prev(:)=zero
  write(Invar%stdout,*) 'Number of shells=',nshell
  write(Invar%stdout,*) '============================================================================'
  if (MPIdata%iam_master) open(unit=16,file=trim(Invar%output_prefix)//'_nbcoeff-phi3.dat')
  ncoeff_prev=0
  do ishell=1,nshell
    ncoeff=0
    iatref=Shell3at%iatref(ishell)
    jatref=Shell3at%jatref(ishell)
    katref=Shell3at%katref(ishell)
    write(Invar%stdout,*) 'Shell number:',ishell
    write(Invar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',iatref,' and ',jatref,' the distance is=',MD%distance(iatref,jatref,1)
    write(Invar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',jatref,' and ',katref,' the distance is=',MD%distance(jatref,katref,1)
    write(Invar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',katref,' and ',iatref,' the distance is=',MD%distance(katref,iatref,1)
    call tdep_calc_nbcoeff(MD%distance,iatref,Invar,ishell,jatref,katref,1,MPIdata,ncoeff,norder,nshell,order,Shell3at%proj,Sym)
    Shell3at%ncoeff     (ishell)=ncoeff
    Shell3at%ncoeff_prev(ishell)=ncoeff_prev
    ncoeff_prev=ncoeff_prev+ncoeff
    write(Invar%stdout,*)'  Number of independant coefficients in this shell=',ncoeff
    write(Invar%stdout,*)'  Number of interactions in this shell=',Shell3at%neighbours(iatref,ishell)%n_interactions
!FB    write(Invar%stdout,*)'  The ratio is=',dfloat(Shell3at%neighbours(iatref,ishell)%n_interactions)/dfloat(ncoeff)
    write(Invar%stdout,*) '============================================================================'
  end do
  write(Invar%stdout,*)'  >>>>>> Total number of coefficients at the third order=',ncoeff_prev
  if (MPIdata%iam_master) close(16)
  Shell3at%ntotcoeff=ncoeff_prev
!BeginFB
!FB  open(unit=91,file='Shell3at.dat')
!FB  write(91,*) Shell3at%nshell
!FB  do ishell=1,Shell3at%nshell
!FB    write(91,*) Shell3at%ncoeff(ishell)
!FB    write(91,*) Shell3at%ncoeff_prev(ishell)
!FB    write(91,*) Shell3at%iatref(ishell)
!FB    write(91,*) Shell3at%jatref(ishell)
!FB    write(91,*) Shell3at%katref(ishell)
!FB    do iatom=1,Invar%natom
!FB      write(91,*) Shell3at%neighbours(iatom,ishell)%n_interactions
!FB      do ii=1,Shell3at%neighbours(iatom,ishell)%n_interactions
!FB        write(91,*) Shell3at%neighbours(iatom,ishell)%sym_in_shell(ii)
!FB        write(91,*) Shell3at%neighbours(iatom,ishell)%transpose_in_shell(ii)
!FB        write(91,*) Shell3at%neighbours(iatom,ishell)%atomj_in_shell(ii)
!FB        write(91,*) Shell3at%neighbours(iatom,ishell)%atomk_in_shell(ii)
!FB      end do
!FB    end do
!FB  end do
!FB  close(91)
!EndFB

 end subroutine tdep_init_shell3at

!====================================================================================================

 subroutine tdep_init_shell4at(Shell4at,Invar,MD,Sym,MPIdata)

  type(Shell_type),intent(out) :: Shell4at
  type(atdep_dataset_type),intent(in) :: Invar
  type(tdep_Sampling_type),intent(in) :: MD
  type(Symmetries_type),intent(inout) :: Sym
  type(MPI_enreg_type), intent(in) :: MPIdata

  integer :: ii,ishell,iatom,jatom,katom,latom,eatom,fatom,gatom,hatom,iatref,jatref,katref,latref
  integer :: natom,natom_unitcell,watom,xatom,yatom,zatom,ninteractions,ncoeff,ncoeff_prev,nshell_tmp
  integer :: already_found,find_equivalent,ninter,iat_ref,jat_ref,kat_ref,lat_ref,tmpinter
  integer :: norder,order,nshell_max,nshell
  double precision :: norma1,norma2,norma3
  double precision :: normb1,normb2,normb3
  double precision :: normc1,normc2,normc3
  integer :: Isym4at(2)
  integer, allocatable :: atref(:,:),interactions(:,:)

  natom = Invar%natom
  natom_unitcell = Invar%natom_unitcell
  nshell_max = Invar%nshell_max
  order = 4
  norder = 81
  Shell4at%order = order
  Shell4at%norder = norder
  Shell4at%natom = natom

  write(Invar%stdout,*) ' '
  write(Invar%stdout,*) '#############################################################################'
  write(Invar%stdout,*) '###### FOURTH ORDER : find the number of coefficients ########################'
  write(Invar%stdout,*) '#############################################################################'

! 1/ Identify the shells
  ABI_CALLOC(interactions,(natom_unitcell,nshell_max))
  ABI_CALLOC(atref,(nshell_max,4))
  ABI_CALLOC(Shell4at%ishell_self,(natom_unitcell))
  nshell_tmp=0
  do iatom=1,natom_unitcell
    do jatom=1,natom
!     Interactions are only computed until Rcut4 in order to have complete shell of neighbours.
!     Otherwise the symetries are broken.
      if (MD%distance(iatom,jatom,1).gt.(Invar%rcut4*0.99)) cycle
      do katom=1,natom
        if (MD%distance(iatom,katom,1).gt.(Invar%rcut4*0.99)) cycle
!         WARNING: distance(j,k).ne.|djk| due to the inbox procedure when computing distance(j,k).
!                  So, compute |djk| using vec(ij) and vec(ik).
          norma1=dsqrt((MD%distance(iatom,katom,2)-MD%distance(iatom,jatom,2))**2+&
&                      (MD%distance(iatom,katom,3)-MD%distance(iatom,jatom,3))**2+&
&                      (MD%distance(iatom,katom,4)-MD%distance(iatom,jatom,4))**2)
        if (norma1                 .gt.(Invar%rcut4*0.99)) cycle
        do latom=1,natom
!FB          write(Invar%stdlog,*) 'NEW COORD1 : iatom,jatom,katom=',iatom,jatom,katom,latom
          if (MD%distance(iatom,latom,1).gt.(Invar%rcut4*0.99)) cycle
          norma2=dsqrt((MD%distance(iatom,latom,2)-MD%distance(iatom,jatom,2))**2+&
&                      (MD%distance(iatom,latom,3)-MD%distance(iatom,jatom,3))**2+&
&                      (MD%distance(iatom,latom,4)-MD%distance(iatom,jatom,4))**2)
          if (norma2                 .gt.(Invar%rcut4*0.99)) cycle
          norma3=dsqrt((MD%distance(iatom,latom,2)-MD%distance(iatom,katom,2))**2+&
&                      (MD%distance(iatom,latom,3)-MD%distance(iatom,katom,3))**2+&
&                      (MD%distance(iatom,latom,4)-MD%distance(iatom,katom,4))**2)
          if (norma3                 .gt.(Invar%rcut4*0.99)) cycle

          if (nshell_tmp.eq.0) then
            atref(1,:)=1
            nshell_tmp=nshell_tmp+1
            if ((iatom.eq.jatom).and.(jatom.eq.katom).and.(katom.eq.latom)) Shell4at%ishell_self(iatom)=nshell_tmp
            interactions(iatom,nshell_tmp)=interactions(iatom,nshell_tmp)+1
            cycle
          end if
          find_equivalent=0
          do ishell=1,nshell_tmp
            iat_ref=atref(ishell,1) ; jat_ref=atref(ishell,2) ; kat_ref=atref(ishell,3) ; lat_ref=atref(ishell,4)
            normb1=dsqrt((MD%distance(iat_ref,kat_ref,2)-MD%distance(iat_ref,jat_ref,2))**2+&
&                        (MD%distance(iat_ref,kat_ref,3)-MD%distance(iat_ref,jat_ref,3))**2+&
&                        (MD%distance(iat_ref,kat_ref,4)-MD%distance(iat_ref,jat_ref,4))**2)
            normb2=dsqrt((MD%distance(iat_ref,lat_ref,2)-MD%distance(iat_ref,jat_ref,2))**2+&
&                        (MD%distance(iat_ref,lat_ref,3)-MD%distance(iat_ref,jat_ref,3))**2+&
&                        (MD%distance(iat_ref,lat_ref,4)-MD%distance(iat_ref,jat_ref,4))**2)
            normb3=dsqrt((MD%distance(iat_ref,lat_ref,2)-MD%distance(iat_ref,kat_ref,2))**2+&
&                        (MD%distance(iat_ref,lat_ref,3)-MD%distance(iat_ref,kat_ref,3))**2+&
&                        (MD%distance(iat_ref,lat_ref,4)-MD%distance(iat_ref,kat_ref,4))**2)
            do ii=1,24
              if (ii.eq.1 ) then ; eatom=iatom ; fatom=jatom ; gatom=katom ; hatom=latom ; end if  !ijkl
              if (ii.eq.2 ) then ; eatom=iatom ; fatom=katom ; gatom=jatom ; hatom=latom ; end if  !ikjl
              if (ii.eq.3 ) then ; eatom=jatom ; fatom=iatom ; gatom=katom ; hatom=latom ; end if  !jikl
              if (ii.eq.4 ) then ; eatom=jatom ; fatom=katom ; gatom=iatom ; hatom=latom ; end if  !jkil
              if (ii.eq.5 ) then ; eatom=katom ; fatom=iatom ; gatom=jatom ; hatom=latom ; end if  !kijl
              if (ii.eq.6 ) then ; eatom=katom ; fatom=jatom ; gatom=iatom ; hatom=latom ; end if  !kjil

              if (ii.eq.7 ) then ; eatom=iatom ; fatom=jatom ; gatom=latom ; hatom=katom ; end if  !ijlk
              if (ii.eq.8 ) then ; eatom=iatom ; fatom=katom ; gatom=latom ; hatom=jatom ; end if  !iklj
              if (ii.eq.9 ) then ; eatom=jatom ; fatom=iatom ; gatom=latom ; hatom=katom ; end if  !jilk
              if (ii.eq.10) then ; eatom=jatom ; fatom=katom ; gatom=latom ; hatom=iatom ; end if  !jkli
              if (ii.eq.11) then ; eatom=katom ; fatom=iatom ; gatom=latom ; hatom=jatom ; end if  !kilj
              if (ii.eq.12) then ; eatom=katom ; fatom=jatom ; gatom=latom ; hatom=iatom ; end if  !kjli

              if (ii.eq.13) then ; eatom=iatom ; fatom=latom ; gatom=jatom ; hatom=katom ; end if  !iljk
              if (ii.eq.14) then ; eatom=iatom ; fatom=latom ; gatom=katom ; hatom=jatom ; end if  !ilkj
              if (ii.eq.15) then ; eatom=jatom ; fatom=latom ; gatom=iatom ; hatom=katom ; end if  !jlik
              if (ii.eq.16) then ; eatom=jatom ; fatom=latom ; gatom=katom ; hatom=iatom ; end if  !jlki
              if (ii.eq.17) then ; eatom=katom ; fatom=latom ; gatom=iatom ; hatom=jatom ; end if  !klij
              if (ii.eq.18) then ; eatom=katom ; fatom=latom ; gatom=jatom ; hatom=iatom ; end if  !klji

              if (ii.eq.19) then ; eatom=latom ; fatom=iatom ; gatom=jatom ; hatom=katom ; end if  !lijk
              if (ii.eq.20) then ; eatom=latom ; fatom=iatom ; gatom=katom ; hatom=jatom ; end if  !likj
              if (ii.eq.21) then ; eatom=latom ; fatom=jatom ; gatom=iatom ; hatom=katom ; end if  !ljik
              if (ii.eq.22) then ; eatom=latom ; fatom=jatom ; gatom=katom ; hatom=iatom ; end if  !ljki
              if (ii.eq.23) then ; eatom=latom ; fatom=katom ; gatom=iatom ; hatom=jatom ; end if  !lkij
              if (ii.eq.24) then ; eatom=latom ; fatom=katom ; gatom=jatom ; hatom=iatom ; end if  !lkji

              normc1=dsqrt((MD%distance(eatom,gatom,2)-MD%distance(eatom,fatom,2))**2+&
&                          (MD%distance(eatom,gatom,3)-MD%distance(eatom,fatom,3))**2+&
&                          (MD%distance(eatom,gatom,4)-MD%distance(eatom,fatom,4))**2)
              normc2=dsqrt((MD%distance(eatom,hatom,2)-MD%distance(eatom,fatom,2))**2+&
&                          (MD%distance(eatom,hatom,3)-MD%distance(eatom,fatom,3))**2+&
&                          (MD%distance(eatom,hatom,4)-MD%distance(eatom,fatom,4))**2)
              normc3=dsqrt((MD%distance(eatom,hatom,2)-MD%distance(eatom,gatom,2))**2+&
&                          (MD%distance(eatom,hatom,3)-MD%distance(eatom,gatom,3))**2+&
&                          (MD%distance(eatom,hatom,4)-MD%distance(eatom,gatom,4))**2)
              if ((abs(MD%distance(iat_ref,jat_ref,1)-MD%distance(eatom,fatom,1)).lt.1.d-6).and.&
&                 (abs(normb1                     -normc1                 ).lt.1.d-6).and.&
&                 (abs(normb2                     -normc2                 ).lt.1.d-6).and.&
&                 (abs(normb3                     -normc3                 ).lt.1.d-6).and.&
&                 (abs(MD%distance(iat_ref,kat_ref,1)-MD%distance(eatom,gatom,1)).lt.1.d-6).and.&
&                 (abs(MD%distance(iat_ref,lat_ref,1)-MD%distance(eatom,hatom,1)).lt.1.d-6)) then
                Isym4at(:)=0
                call tdep_SearchS_4at(Invar,iat_ref,jat_ref,kat_ref,lat_ref,eatom,fatom,gatom,hatom,Isym4at,Sym,MD%xred_ideal)
                if (Isym4at(2).eq.1) find_equivalent=1
                if (find_equivalent.eq.1) then
                  interactions(iatom,ishell)=interactions(iatom,ishell)+1
!FB                  write(Invar%stdlog,*) 'The number of interactions in this shell is=',ishell,interactions(iatom,ishell)
                  exit
                end if
              end if
            end do !ii
            if (find_equivalent.eq.1) exit
          end do !ishell
          if (find_equivalent.eq.0) then
            nshell_tmp=nshell_tmp+1
            if (nshell_tmp.gt.nshell_max) then
              ABI_ERROR('The shell number index is greater than the shell number max defined in the code')
            end if
            if ((iatom.eq.jatom).and.(jatom.eq.katom).and.(katom.eq.latom)) Shell4at%ishell_self(iatom)=nshell_tmp
            interactions(iatom,nshell_tmp)=interactions(iatom,nshell_tmp)+1
            atref(nshell_tmp,1)=iatom
            atref(nshell_tmp,2)=jatom
            atref(nshell_tmp,3)=katom
            atref(nshell_tmp,4)=latom
!FB            write(Invar%stdlog,'(a,1x,5(i5,1x))') 'NEW SHELL1 : nshell_tmp,iatom,jatom,katom,latom=',nshell_tmp,iatom,jatom,katom,latom
          end if
        end do !latom
      end do !katom
    end do !jatom
  end do !iatom
  ABI_FREE(atref)

! 2/ Allocate the datatype Shell4at%...
  Shell4at%nshell = nshell_tmp
  nshell = Shell4at%nshell
  ABI_MALLOC(Shell4at%neighbours,(natom,nshell))
  ABI_CALLOC(Shell4at%iatref,(nshell))
  ABI_CALLOC(Shell4at%jatref,(nshell))
  ABI_CALLOC(Shell4at%katref,(nshell))
  ABI_CALLOC(Shell4at%latref,(nshell))
  do ishell=1,nshell
    do iatom=1,natom
      ninteractions=interactions(mod(iatom-1,natom_unitcell)+1,ishell)
      Shell4at%neighbours(iatom,ishell)%n_interactions=ninteractions
      if (ninteractions.eq.0) cycle
      ABI_MALLOC(Shell4at%neighbours(iatom,ishell)%atomj_in_shell,(ninteractions))
      ABI_MALLOC(Shell4at%neighbours(iatom,ishell)%atomk_in_shell,(ninteractions))
      ABI_MALLOC(Shell4at%neighbours(iatom,ishell)%atoml_in_shell,(ninteractions))
      ABI_MALLOC(Shell4at%neighbours(iatom,ishell)%sym_in_shell,(ninteractions))
      ABI_MALLOC(Shell4at%neighbours(iatom,ishell)%transpose_in_shell,(ninteractions))
      Shell4at%neighbours(iatom,ishell)%atomj_in_shell(:)=zero
      Shell4at%neighbours(iatom,ishell)%atomk_in_shell(:)=zero
      Shell4at%neighbours(iatom,ishell)%atoml_in_shell(:)=zero
      Shell4at%neighbours(iatom,ishell)%sym_in_shell(:)=zero
      Shell4at%neighbours(iatom,ishell)%transpose_in_shell(:)=zero
    end do
  end do
  ABI_FREE(interactions)

! 3/ Store the index of the (couple of) atoms included in each shell
! 4/ Store the reference (couple of) atoms for each shell
! 5/ Compute the symetry operation between the reference (couple of) atoms and another one
  ABI_MALLOC(interactions,(natom,nshell)) ; interactions(:,:)=0
  nshell_tmp=0
  do iatom=1,natom
    do jatom=1,natom
      if (MD%distance(iatom,jatom,1).gt.(Invar%rcut4*0.99)) cycle
      do katom=1,natom
        if (MD%distance(iatom,katom,1).gt.(Invar%rcut4*0.99)) cycle
!         WARNING: distance(j,k).ne.|djk| due to the inbox procedure when computing distance(j,k).
!                  So, compute |djk| using vec(ij) and vec(ik).
        norma1=dsqrt((MD%distance(iatom,katom,2)-MD%distance(iatom,jatom,2))**2+&
&                    (MD%distance(iatom,katom,3)-MD%distance(iatom,jatom,3))**2+&
&                    (MD%distance(iatom,katom,4)-MD%distance(iatom,jatom,4))**2)
        if (norma1                 .gt.(Invar%rcut4*0.99)) cycle
        do latom=1,natom
!FB          write(Invar%stdlog,*) 'NEW COORD2 : iatom,jatom,katom=',iatom,jatom,katom,latom
          if (MD%distance(iatom,latom,1).gt.(Invar%rcut4*0.99)) cycle
          norma2=dsqrt((MD%distance(iatom,latom,2)-MD%distance(iatom,jatom,2))**2+&
&                      (MD%distance(iatom,latom,3)-MD%distance(iatom,jatom,3))**2+&
&                      (MD%distance(iatom,latom,4)-MD%distance(iatom,jatom,4))**2)
          if (norma2                 .gt.(Invar%rcut4*0.99)) cycle
          norma3=dsqrt((MD%distance(iatom,latom,2)-MD%distance(iatom,katom,2))**2+&
&                      (MD%distance(iatom,latom,3)-MD%distance(iatom,katom,3))**2+&
&                      (MD%distance(iatom,latom,4)-MD%distance(iatom,katom,4))**2)
          if (norma3                 .gt.(Invar%rcut4*0.99)) cycle
          if (nshell_tmp.eq.0) then
            nshell_tmp=nshell_tmp+1
            interactions(iatom,nshell_tmp)=1
            Shell4at%iatref(nshell_tmp)=iatom
            Shell4at%jatref(nshell_tmp)=jatom
            Shell4at%katref(nshell_tmp)=katom
            Shell4at%latref(nshell_tmp)=latom
            Shell4at%neighbours(iatom,nshell_tmp)%atomj_in_shell(interactions(iatom,nshell_tmp))=jatom
            Shell4at%neighbours(iatom,nshell_tmp)%atomk_in_shell(interactions(iatom,nshell_tmp))=katom
            Shell4at%neighbours(iatom,nshell_tmp)%atoml_in_shell(interactions(iatom,nshell_tmp))=latom
            Shell4at%neighbours(iatom,nshell_tmp)%sym_in_shell(interactions(iatom,nshell_tmp))=1
            Shell4at%neighbours(iatom,nshell_tmp)%transpose_in_shell(interactions(iatom,nshell_tmp))=1
            cycle
          end if

!         Search if the quadruplet has already been classified
          find_equivalent=0
          do ishell=1,nshell_tmp
            do ninter=1,interactions(iatom,ishell)
              if ((Shell4at%neighbours(iatom,ishell)%atomj_in_shell(ninter).eq.jatom).and.&
&                 (Shell4at%neighbours(iatom,ishell)%atomk_in_shell(ninter).eq.katom).and.&
&                 (Shell4at%neighbours(iatom,ishell)%atoml_in_shell(ninter).eq.latom)) find_equivalent=1
            end do
          end do
          if (find_equivalent.eq.1) cycle
!         Search if the quadruplet belongs to a shell already found
          do ishell=1,nshell_tmp
            iat_ref=Shell4at%iatref(ishell) ; jat_ref=Shell4at%jatref(ishell)
            kat_ref=Shell4at%katref(ishell) ; lat_ref=Shell4at%latref(ishell)
            normb1=dsqrt((MD%distance(iat_ref,kat_ref,2)-MD%distance(iat_ref,jat_ref,2))**2+&
&                        (MD%distance(iat_ref,kat_ref,3)-MD%distance(iat_ref,jat_ref,3))**2+&
&                        (MD%distance(iat_ref,kat_ref,4)-MD%distance(iat_ref,jat_ref,4))**2)
            normb2=dsqrt((MD%distance(iat_ref,lat_ref,2)-MD%distance(iat_ref,jat_ref,2))**2+&
&                        (MD%distance(iat_ref,lat_ref,3)-MD%distance(iat_ref,jat_ref,3))**2+&
&                        (MD%distance(iat_ref,lat_ref,4)-MD%distance(iat_ref,jat_ref,4))**2)
            normb3=dsqrt((MD%distance(iat_ref,lat_ref,2)-MD%distance(iat_ref,kat_ref,2))**2+&
&                        (MD%distance(iat_ref,lat_ref,3)-MD%distance(iat_ref,kat_ref,3))**2+&
&                        (MD%distance(iat_ref,lat_ref,4)-MD%distance(iat_ref,kat_ref,4))**2)
            do ii=1,24
              if (ii.eq.1 ) then ; eatom=iatom ; fatom=jatom ; gatom=katom ; hatom=latom ; end if  !ijkl
              if (ii.eq.2 ) then ; eatom=iatom ; fatom=katom ; gatom=jatom ; hatom=latom ; end if  !ikjl
              if (ii.eq.3 ) then ; eatom=jatom ; fatom=iatom ; gatom=katom ; hatom=latom ; end if  !jikl
              if (ii.eq.4 ) then ; eatom=jatom ; fatom=katom ; gatom=iatom ; hatom=latom ; end if  !jkil
              if (ii.eq.5 ) then ; eatom=katom ; fatom=iatom ; gatom=jatom ; hatom=latom ; end if  !kijl
              if (ii.eq.6 ) then ; eatom=katom ; fatom=jatom ; gatom=iatom ; hatom=latom ; end if  !kjil

              if (ii.eq.7 ) then ; eatom=iatom ; fatom=jatom ; gatom=latom ; hatom=katom ; end if  !ijlk
              if (ii.eq.8 ) then ; eatom=iatom ; fatom=katom ; gatom=latom ; hatom=jatom ; end if  !iklj
              if (ii.eq.9 ) then ; eatom=jatom ; fatom=iatom ; gatom=latom ; hatom=katom ; end if  !jilk
              if (ii.eq.10) then ; eatom=jatom ; fatom=katom ; gatom=latom ; hatom=iatom ; end if  !jkli
              if (ii.eq.11) then ; eatom=katom ; fatom=iatom ; gatom=latom ; hatom=jatom ; end if  !kilj
              if (ii.eq.12) then ; eatom=katom ; fatom=jatom ; gatom=latom ; hatom=iatom ; end if  !kjli

              if (ii.eq.13) then ; eatom=iatom ; fatom=latom ; gatom=jatom ; hatom=katom ; end if  !iljk
              if (ii.eq.14) then ; eatom=iatom ; fatom=latom ; gatom=katom ; hatom=jatom ; end if  !ilkj
              if (ii.eq.15) then ; eatom=jatom ; fatom=latom ; gatom=iatom ; hatom=katom ; end if  !jlik
              if (ii.eq.16) then ; eatom=jatom ; fatom=latom ; gatom=katom ; hatom=iatom ; end if  !jlki
              if (ii.eq.17) then ; eatom=katom ; fatom=latom ; gatom=iatom ; hatom=jatom ; end if  !klij
              if (ii.eq.18) then ; eatom=katom ; fatom=latom ; gatom=jatom ; hatom=iatom ; end if  !klji

              if (ii.eq.19) then ; eatom=latom ; fatom=iatom ; gatom=jatom ; hatom=katom ; end if  !lijk
              if (ii.eq.20) then ; eatom=latom ; fatom=iatom ; gatom=katom ; hatom=jatom ; end if  !likj
              if (ii.eq.21) then ; eatom=latom ; fatom=jatom ; gatom=iatom ; hatom=katom ; end if  !ljik
              if (ii.eq.22) then ; eatom=latom ; fatom=jatom ; gatom=katom ; hatom=iatom ; end if  !ljki
              if (ii.eq.23) then ; eatom=latom ; fatom=katom ; gatom=iatom ; hatom=jatom ; end if  !lkij
              if (ii.eq.24) then ; eatom=latom ; fatom=katom ; gatom=jatom ; hatom=iatom ; end if  !lkji

              normc1=dsqrt((MD%distance(eatom,gatom,2)-MD%distance(eatom,fatom,2))**2+&
&                          (MD%distance(eatom,gatom,3)-MD%distance(eatom,fatom,3))**2+&
&                          (MD%distance(eatom,gatom,4)-MD%distance(eatom,fatom,4))**2)
              normc2=dsqrt((MD%distance(eatom,hatom,2)-MD%distance(eatom,fatom,2))**2+&
&                          (MD%distance(eatom,hatom,3)-MD%distance(eatom,fatom,3))**2+&
&                          (MD%distance(eatom,hatom,4)-MD%distance(eatom,fatom,4))**2)
              normc3=dsqrt((MD%distance(eatom,hatom,2)-MD%distance(eatom,gatom,2))**2+&
&                          (MD%distance(eatom,hatom,3)-MD%distance(eatom,gatom,3))**2+&
&                          (MD%distance(eatom,hatom,4)-MD%distance(eatom,gatom,4))**2)
              if ((abs(MD%distance(iat_ref,jat_ref,1)-MD%distance(eatom,fatom,1)).lt.1.d-6).and.&
&                 (abs(normb1                     -normc1                 ).lt.1.d-6).and.&
&                 (abs(normb2                     -normc2                 ).lt.1.d-6).and.&
&                 (abs(normb3                     -normc3                 ).lt.1.d-6).and.&
&                 (abs(MD%distance(iat_ref,kat_ref,1)-MD%distance(eatom,gatom,1)).lt.1.d-6).and.&
&                 (abs(MD%distance(iat_ref,lat_ref,1)-MD%distance(eatom,hatom,1)).lt.1.d-6)) then
                Isym4at(:)=0
                call tdep_SearchS_4at(Invar,iat_ref,jat_ref,kat_ref,lat_ref,eatom,fatom,gatom,hatom,Isym4at,Sym,MD%xred_ideal)
                if (Isym4at(2).eq.1) then
                  find_equivalent=1
                  exit
                end if
              end if
            end do !ii
            if (find_equivalent.eq.1) exit
          end do !ishell
!         The quadruplet belongs to a new shell
          if (find_equivalent.eq.0) then
            nshell_tmp=nshell_tmp+1
!           Check that the new shell is allowed
            if (nshell_tmp.gt.nshell) then
              ABI_ERROR('The shell number index is greater than the shell number max computed previously')
            end if
            Shell4at%iatref(nshell_tmp)=iatom
            Shell4at%jatref(nshell_tmp)=jatom
            Shell4at%katref(nshell_tmp)=katom
            Shell4at%latref(nshell_tmp)=latom
            eatom=iatom ; fatom=jatom ; gatom=katom ; hatom=latom
            Isym4at(:)=1
            ishell=nshell_tmp
!FB            write(Invar%stdlog,'(a,1x,5(i5,1x))') 'NEW SHELL2 : nshell_tmp,iatom,jatom,katom=',nshell_tmp,iatom,jatom,katom,latom
          end if
!         Classify the informations of the quadruplet in Shell4at
          do ii=1,24
!           The Phi4 has to be symetric (transposition symetries)
            if (ii.eq.1 ) then ; watom=eatom ; xatom=fatom ; yatom=gatom ; zatom=hatom ; end if  !ijkl
            if (ii.eq.2 ) then ; watom=eatom ; xatom=gatom ; yatom=fatom ; zatom=hatom ; end if  !ikjl
            if (ii.eq.3 ) then ; watom=fatom ; xatom=eatom ; yatom=gatom ; zatom=hatom ; end if  !jikl
            if (ii.eq.4 ) then ; watom=fatom ; xatom=gatom ; yatom=eatom ; zatom=hatom ; end if  !jkil
            if (ii.eq.5 ) then ; watom=gatom ; xatom=eatom ; yatom=fatom ; zatom=hatom ; end if  !kijl
            if (ii.eq.6 ) then ; watom=gatom ; xatom=fatom ; yatom=eatom ; zatom=hatom ; end if  !kjil

            if (ii.eq.7 ) then ; watom=eatom ; xatom=fatom ; yatom=hatom ; zatom=gatom ; end if  !ijlk
            if (ii.eq.8 ) then ; watom=eatom ; xatom=gatom ; yatom=hatom ; zatom=fatom ; end if  !iklj
            if (ii.eq.9 ) then ; watom=fatom ; xatom=eatom ; yatom=hatom ; zatom=gatom ; end if  !jilk
            if (ii.eq.10) then ; watom=fatom ; xatom=gatom ; yatom=hatom ; zatom=eatom ; end if  !jkli
            if (ii.eq.11) then ; watom=gatom ; xatom=eatom ; yatom=hatom ; zatom=fatom ; end if  !kilj
            if (ii.eq.12) then ; watom=gatom ; xatom=fatom ; yatom=hatom ; zatom=eatom ; end if  !kjli

            if (ii.eq.13) then ; watom=eatom ; xatom=hatom ; yatom=fatom ; zatom=gatom ; end if  !iljk
            if (ii.eq.14) then ; watom=eatom ; xatom=hatom ; yatom=gatom ; zatom=fatom ; end if  !ilkj
            if (ii.eq.15) then ; watom=fatom ; xatom=hatom ; yatom=eatom ; zatom=gatom ; end if  !jlik
            if (ii.eq.16) then ; watom=fatom ; xatom=hatom ; yatom=gatom ; zatom=eatom ; end if  !jlki
            if (ii.eq.17) then ; watom=gatom ; xatom=hatom ; yatom=eatom ; zatom=fatom ; end if  !klij
            if (ii.eq.18) then ; watom=gatom ; xatom=hatom ; yatom=fatom ; zatom=eatom ; end if  !klji

            if (ii.eq.19) then ; watom=hatom ; xatom=eatom ; yatom=fatom ; zatom=gatom ; end if  !lijk
            if (ii.eq.20) then ; watom=hatom ; xatom=eatom ; yatom=gatom ; zatom=fatom ; end if  !likj
            if (ii.eq.21) then ; watom=hatom ; xatom=fatom ; yatom=eatom ; zatom=gatom ; end if  !ljik
            if (ii.eq.22) then ; watom=hatom ; xatom=fatom ; yatom=gatom ; zatom=eatom ; end if  !ljki
            if (ii.eq.23) then ; watom=hatom ; xatom=gatom ; yatom=eatom ; zatom=fatom ; end if  !lkij
            if (ii.eq.24) then ; watom=hatom ; xatom=gatom ; yatom=fatom ; zatom=eatom ; end if  !lkji
!           Do not overwrite the Phi4_iikl, Phi4_ijil, Phi4_ijjl, Phi4_iiil... IFCs
!           and avoid double counting of quadruplet interactions
            already_found=0
            do ninter=1,interactions(watom,ishell)
              if ((Shell4at%neighbours(watom,ishell)%atomj_in_shell(ninter).eq.xatom).and.&
&                 (Shell4at%neighbours(watom,ishell)%atomk_in_shell(ninter).eq.yatom).and.&
&                 (Shell4at%neighbours(watom,ishell)%atoml_in_shell(ninter).eq.zatom)) then
!FB                write(*,'(a,4(1x,i5))') 'FOR efgh =',eatom,fatom,gatom,hatom
!FB                write(*,'(a,3(1x,i5))') '  --> ishell,ninter,ninter_tot =',ishell,ninter,interactions(watom,ishell)
!FB                write(*,'(a,4(1x,i5))') '  --> ALREADY FOUND =',watom,xatom,yatom,zatom
                already_found=1
                exit
              end if
            end do
            if (already_found==1) cycle
            interactions(watom,ishell)=interactions(watom,ishell)+1
!FB            write(Invar%stdlog,*) 'For ishell and eatom=',ishell,watom
!FB            write(Invar%stdlog,*) '  --> the number of interactions in the shell is=',interactions(watom,ishell)
            if (interactions(watom,ishell).gt.Shell4at%neighbours(watom,ishell)%n_interactions) then
              write(Invar%stdlog,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
              write(Invar%stdlog,*) ' >>>>>> Verify that the Rcut used in the input file is lower '
              write(Invar%stdlog,*) ' >>>>>> than half of the smallest lattice parameter'
              write(Invar%stdlog,*) ' >>>>>> Solution : Reduce the Rcut parameter'
              write(Invar%stdlog,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
              ABI_ERROR('The interaction number index is greater than the interaction number max computed previously (4th order)')
            end if
            Shell4at%neighbours(watom,ishell)%atomj_in_shell(interactions(watom,ishell))=xatom
            Shell4at%neighbours(watom,ishell)%atomk_in_shell(interactions(watom,ishell))=yatom
            Shell4at%neighbours(watom,ishell)%atoml_in_shell(interactions(watom,ishell))=zatom
            Shell4at%neighbours(watom,ishell)%sym_in_shell(interactions(watom,ishell))=Isym4at(1)
            Shell4at%neighbours(watom,ishell)%transpose_in_shell(interactions(watom,ishell))=ii
!DEBUG            write(Invar%stdlog,'(a,9(i5,x))') 'ishell,iatref,jatref,katref,iatom,atomj_in_shell,atomk_in_shell,isym,itrans=',&
!DEBUG&           ishell,Shell4at%iatref(ishell),Shell4at%jatref(ishell),Shell4at%katref(ishell),watom,xatom,yatom,Isym4at(1),ii
          end do !ii
        end do !latom
      end do !katom
    end do !jatom
  end do !iatom
! Check that each interaction has different symmetry per shell
  do ishell=1,nshell
    do iatom=1,natom
      if (Shell4at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
      do ninter=1,Shell4at%neighbours(iatom,ishell)%n_interactions-1
        do tmpinter=ninter+1,Shell4at%neighbours(iatom,ishell)%n_interactions
          if (Shell4at%neighbours(iatom,ishell)%sym_in_shell(  ninter).eq.&
&             Shell4at%neighbours(iatom,ishell)%sym_in_shell(tmpinter)) then
            if (Shell4at%neighbours(iatom,ishell)%transpose_in_shell(  ninter).ne.&
&                    Shell4at%neighbours(iatom,ishell)%transpose_in_shell(tmpinter)) cycle
            write(std_out,'(a,2(1x,i5))') 'For ishell and iatom =',ishell,iatom
            write(std_out,'(a,i5,a,i5,a,i5)') '  the interactions ',ninter,&
&             ' and ',tmpinter,' have both the same symmetry isym=',Shell4at%neighbours(iatom,ishell)%sym_in_shell(  ninter)
            ABI_ERROR('Some interactions are equals due to the symmetry')
          end if
        end do
      end do
    end do
  end do
! Check that each equivalent shell has the same set of interactions
  do ishell=1,nshell
    do iatom=1,natom
      if (Shell4at%neighbours(mod(iatom-1,natom_unitcell)+1,ishell)%n_interactions.ne.&
&         Shell4at%neighbours(                        iatom,ishell)%n_interactions) then
        ABI_ERROR('The interaction number index is not equal to the interaction number max computed previously (2)')
      end if
!DEBUG      iatref=Shell4at%iatref(ishell)
!DEBUG      jatref=Shell4at%jatref(ishell)
!DEBUG      katref=Shell4at%katref(ishell)
!DEBUG      if (Shell4at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
!DEBUG      do ninter=1,Shell4at%neighbours(iatom,ishell)%n_interactions
!DEBUG        jatom=Shell4at%neighbours(iatom,ishell)%atomj_in_shell(ninter)
!DEBUG        katom=Shell4at%neighbours(iatom,ishell)%atomk_in_shell(ninter)
!DEBUG        isym =Shell4at%neighbours(iatom,ishell)%sym_in_shell(ninter)
!DEBUG        itrans=Shell4at%neighbours(iatom,ishell)%transpose_in_shell(ninter)
!DEBUG        vectj(:)=zero ; vectk(:)=zero ; vect1(:)=zero ; vect2(:)=zero
!DEBUG        do ii=1,3
!DEBUG          do jj=1,3
!DEBUG            vectj(ii)=vectj(ii)+Sym%S_ref(ii,jj,isym,1)*MD%distance(iatref,jatref,jj+1)
!DEBUG            vectk(ii)=vectk(ii)+Sym%S_ref(ii,jj,isym,1)*MD%distance(iatref,katref,jj+1)
!DEBUG          end do
!DEBUG        end do
!DEBUG        if (itrans==1) then ; vect1(:)= vectj(:)          ; vect2(:)= vectk(:)           ; endif !\Phi3_ijk
!DEBUG        if (itrans==2) then ; vect1(:)= vectk(:)          ; vect2(:)= vectj(:)           ; endif !\Phi3_ikj
!DEBUG        if (itrans==3) then ; vect1(:)=-vectj(:)          ; vect2(:)= vectk(:)-vectj(:)  ; endif !\Phi3_jik
!DEBUG        if (itrans==4) then ; vect1(:)= vectk(:)-vectj(:) ; vect2(:)=-vectj(:)           ; endif !\Phi3_jki
!DEBUG        if (itrans==5) then ; vect1(:)=-vectk(:)          ; vect2(:)= vectj(:)-vectk(:)  ; endif !\Phi3_kij
!DEBUG        if (itrans==6) then ; vect1(:)= vectj(:)-vectk(:) ; vect2(:)=-vectk(:)           ; endif !\Phi3_kji
!DEBUG        do ii=1,3
!DEBUG          if ((abs(MD%distance(iatom,jatom,ii+1)-vect1(ii)).gt.tol8).or.&
!DEBUG&             (abs(MD%distance(iatom,katom,ii+1)-vect2(ii)).gt.tol8)) then
!DEBUG            write(std_out,'(a,4(x,i5))') 'For ishell, iatom, jatom, katom =',ishell,iatom,jatom,katom
!DEBUG            write(std_out,'(a,5(x,i5))') '  with isym, itrans, iatref, jatref, katref = ',isym,itrans,iatref,jatref,katref
!DEBUG            ABI_ERROR('We do not recover the quadruplet with the symmetry found')
!DEBUG          end if
!DEBUG        end do !ii
!DEBUG      end do !ninter
    end do !natom
  end do !nshell
  ABI_FREE(interactions)

! Find the number of coefficients of the (3x3x3x3) Phi4 for a given shell
  ABI_CALLOC(Shell4at%proj, (norder,norder,nshell))
  ABI_CALLOC(Shell4at%ncoeff     ,(nshell))
  ABI_CALLOC(Shell4at%ncoeff_prev,(nshell))
  write(Invar%stdout,*) 'Number of shells=',nshell
  write(Invar%stdout,*) '============================================================================'
  if (MPIdata%iam_master) open(unit=16,file=trim(Invar%output_prefix)//'_nbcoeff-phi4.dat')
  ncoeff_prev=0
  do ishell=1,nshell
    ncoeff=0
    iatref=Shell4at%iatref(ishell)
    jatref=Shell4at%jatref(ishell)
    katref=Shell4at%katref(ishell)
    latref=Shell4at%latref(ishell)
    write(Invar%stdout,*) 'Shell number:',ishell
    write(Invar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',iatref,' and ',jatref,' the distance is=',MD%distance(iatref,jatref,1)
    write(Invar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',jatref,' and ',katref,' the distance is=',MD%distance(jatref,katref,1)
    write(Invar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',katref,' and ',latref,' the distance is=',MD%distance(katref,latref,1)
    write(Invar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',latref,' and ',iatref,' the distance is=',MD%distance(latref,iatref,1)
    call tdep_calc_nbcoeff(MD%distance,iatref,Invar,ishell,jatref,katref,latref,MPIdata,ncoeff,norder,nshell,order,Shell4at%proj,Sym)
    Shell4at%ncoeff     (ishell)=ncoeff
    Shell4at%ncoeff_prev(ishell)=ncoeff_prev
    ncoeff_prev=ncoeff_prev+ncoeff
    write(Invar%stdout,*)'  Number of independant coefficients in this shell=',ncoeff
    write(Invar%stdout,*)'  Number of interactions in this shell=',Shell4at%neighbours(iatref,ishell)%n_interactions
!FB    write(Invar%stdout,*)'  The ratio is=',dfloat(Shell4at%neighbours(iatref,ishell)%n_interactions)/dfloat(ncoeff)
    write(Invar%stdout,*) '============================================================================'
  end do
  write(Invar%stdout,*)'  >>>>>> Total number of coefficients at the fourth order=',ncoeff_prev
  if (MPIdata%iam_master) close(16)
  Shell4at%ntotcoeff=ncoeff_prev
!BeginFB
!FB  open(unit=91,file='Shell4at.dat')
!FB  write(91,*) Shell4at%nshell
!FB  do ishell=1,Shell4at%nshell
!FB    write(91,*) Shell4at%ncoeff(ishell)
!FB    write(91,*) Shell4at%ncoeff_prev(ishell)
!FB    write(91,*) Shell4at%iatref(ishell)
!FB    write(91,*) Shell4at%jatref(ishell)
!FB    write(91,*) Shell4at%katref(ishell)
!FB    do iatom=1,Invar%natom
!FB      write(91,*) Shell4at%neighbours(iatom,ishell)%n_interactions
!FB      do ii=1,Shell4at%neighbours(iatom,ishell)%n_interactions
!FB        write(91,*) Shell4at%neighbours(iatom,ishell)%sym_in_shell(ii)
!FB        write(91,*) Shell4at%neighbours(iatom,ishell)%transpose_in_shell(ii)
!FB        write(91,*) Shell4at%neighbours(iatom,ishell)%atomj_in_shell(ii)
!FB        write(91,*) Shell4at%neighbours(iatom,ishell)%atomk_in_shell(ii)
!FB      end do
!FB    end do
!FB  end do
!FB  close(91)
!EndFB

 end subroutine tdep_init_shell4at

!====================================================================================================
 subroutine tdep_destroy_shell(Shell)

  type(Shell_type),intent(inout) :: Shell

  integer :: iatom,ishell,natref

  ABI_FREE(Shell%ncoeff)
  ABI_FREE(Shell%ncoeff_prev)
  ABI_FREE(Shell%iatref)
  if (Shell%order.gt.1) then
    ABI_FREE(Shell%jatref)
    ABI_FREE(Shell%ishell_self)
  end if
  if (Shell%order.gt.2) then
    ABI_FREE(Shell%katref)
  end if
  if (Shell%order.gt.3) then
    ABI_FREE(Shell%latref)
  end if
  if (Shell%order.eq.1) then
    natref=1
  else
    natref=Shell%natom
  end if
  do iatom=1,natref
    do ishell=1,Shell%nshell
      if (Shell%neighbours(iatom,ishell)%n_interactions.ne.0) then
        ABI_FREE(Shell%neighbours(iatom,ishell)%atomj_in_shell)
        ABI_FREE(Shell%neighbours(iatom,ishell)%sym_in_shell)
        if (Shell%order.gt.1) then
          ABI_FREE(Shell%neighbours(iatom,ishell)%transpose_in_shell)
        end if
        if (Shell%order.gt.2) then
          ABI_FREE(Shell%neighbours(iatom,ishell)%atomk_in_shell)
        end if
        if (Shell%order.gt.3) then
          ABI_FREE(Shell%neighbours(iatom,ishell)%atoml_in_shell)
        end if
      end if
    end do
  end do
  ABI_FREE(Shell%neighbours)
  ABI_FREE(Shell%proj)

 end subroutine tdep_destroy_shell

!====================================================================================================

subroutine tdep_calc_nbcoeff(distance,iatcell,Invar,ishell,jatom,katom,latom,MPIdata,&
&                            ncoeff,norder,nshell,order,proj,Sym)

  integer,intent(in) :: iatcell,ishell,jatom,katom,latom,nshell,order,norder
  integer,intent(inout) :: ncoeff
  type(atdep_dataset_type),intent(in) :: Invar
  type(Symmetries_type),intent(in) :: Sym
  type(MPI_enreg_type), intent(in) :: MPIdata
  double precision,intent(in) :: distance(Invar%natom,Invar%natom,4)
  double precision,intent(out) :: proj(norder,norder,nshell)

  integer :: ii,jj,kk,ll,isym,LWORK,INFO,const_tot,itemp,nconst_perm,nconst_loc
  integer :: ncount,icoeff,jatcell,katcell,latcell,mu,nu,xi,zeta
  integer :: inv,watom,xatom,yatom,zatom,isyminv,nsyminv,facorder,iseed
  integer, allocatable :: iconst(:)
  double precision :: prod_scal,drandom
  double precision :: eigvec(3,3)
  double precision :: vect_trial(3),vect_trial1(3),vect_trial2(3),vect_trial3(3)
  double precision :: vect_trial4(3),vect_trial5(3),vect_trial6(3)
  double precision :: WR(3),WI(3),VL(3,3),VR(3,3)
  double precision, allocatable :: WORK(:)
  double complex :: eigenvectors(3,3),eigenvalues(3)
  double complex :: pp(3,3),ppp(3,3,3),pppp(3,3,3,3),lambda
  double complex, allocatable :: tab_vec(:,:),temp(:,:),alphaij(:,:,:),constraints(:,:,:)
  logical :: ok
  logical, allocatable :: unchanged(:)
  character(len=500) :: message

  if (iatcell==1.and.order==1) return
  if (jatom==iatcell.and.order==2) return
!FB  if (katom==iatcell.and.jatom==iatcell.and.order==3) return

  if (order==1) then
    facorder=1
  else if (order==2) then
    facorder=2
  else if (order==3) then
    facorder=6
  else if (order==4) then
    facorder=24
  end if

! If we want to remove the constraints coming from the symetries
!FB  if (order==3) then
!FB    do ii=1,norder
!FB      proj(ii,ii,ishell)=1.d0
!FB    end do
!FB    ncoeff=norder
!FB    return
!FB  end if

  nconst_loc=0
  const_tot=0
  nsyminv=Sym%nsym*facorder
  ABI_MALLOC(alphaij,(nsyminv,norder,norder)); alphaij(:,:,:)=czero
  ABI_MALLOC(iconst,(nsyminv))               ; iconst(:)=0
  ABI_MALLOC(unchanged,(nsyminv))            ; unchanged(:)=.false.

! ================================================================================================
! ================ Big loop over symmetries and invariance (nsym*facorder) =======================
! ================================================================================================
  if (MPIdata%iam_master) write(16,'(a)') ' '
  if (MPIdata%iam_master) write(16,'(a,i4)') 'For shell number=',ishell
  do isyminv=1,nsyminv
    isym=(isyminv-1)/facorder+1
    inv=isyminv-(isym-1)*facorder
    if (isym==1) cycle

!   For the 1st order: Search if the atom is let invariant
    if (order==1) then
      if (Sym%indsym(4,isym,iatcell)==iatcell) then
        if (MPIdata%iam_master) then
          write(16,'(a,1x,i3)')'===========The atom is kept invariant for isym=',isym
        end if
      else
        cycle
      end if
    end if

!   For the 2nd order: Search if the bond is kept invariant or reversed
    if (order==2) then
      vect_trial(:)=zero
      if (inv==1) then ; watom=iatcell ; xatom=jatom   ; endif !\Phi_ij
      if (inv==2) then ; watom=jatom   ; xatom=iatcell ; endif !\Phi_ji
      do ii=1,3
        do jj=1,3
          vect_trial(ii)=vect_trial(ii)+Sym%S_ref(ii,jj,isym,1)*distance(watom,xatom,jj+1)
        end do
      end do
      jatcell=mod(jatom-1,Invar%natom_unitcell)+1
      if ((sum(abs(vect_trial(:)-distance(iatcell,jatom,2:4))).lt.tol8).and.&
&         (Sym%indsym(4,isym,watom)==iatcell).and.&
&         (Sym%indsym(4,isym,xatom)==jatcell)) then
        if (MPIdata%iam_master) then
          if (inv==1) write(16,'(a,1x,i3)')'===========The bond is kept invariant for isym=',isym
          if (inv==2) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j) --> (j,i) for isym=',isym
        end if
      else
        cycle
      end if
    end if

!   For the 3rd order : 6 permutations at all
    if (order==3) then
      vect_trial1(:)=zero
      vect_trial2(:)=zero
      vect_trial3(:)=zero
      if (inv==1) then ; watom=iatcell ; xatom=jatom   ; yatom=katom   ; endif !\Phi3_ijk
      if (inv==2) then ; watom=iatcell ; xatom=katom   ; yatom=jatom   ; endif !\Phi3_ikj
      if (inv==3) then ; watom=jatom   ; xatom=iatcell ; yatom=katom   ; endif !\Phi3_jik
      if (inv==4) then ; watom=jatom   ; xatom=katom   ; yatom=iatcell ; endif !\Phi3_jki
      if (inv==5) then ; watom=katom   ; xatom=iatcell ; yatom=jatom   ; endif !\Phi3_kij
      if (inv==6) then ; watom=katom   ; xatom=jatom   ; yatom=iatcell ; endif !\Phi3_kji
      do ii=1,3
        do jj=1,3
          vect_trial1(ii)=vect_trial1(ii)+Sym%S_ref(ii,jj,isym,1)*distance(watom,xatom,jj+1)
          vect_trial2(ii)=vect_trial2(ii)+Sym%S_ref(ii,jj,isym,1)*distance(xatom,yatom,jj+1)
          vect_trial3(ii)=vect_trial3(ii)+Sym%S_ref(ii,jj,isym,1)*distance(yatom,watom,jj+1)
        end do
      end do
      jatcell=mod(jatom-1,Invar%natom_unitcell)+1
      katcell=mod(katom-1,Invar%natom_unitcell)+1
      if ((sum(abs(vect_trial1(:)-distance(iatcell,jatom  ,2:4))).lt.tol8).and.&
&         (sum(abs(vect_trial2(:)-distance(jatom  ,katom  ,2:4))).lt.tol8).and.&
&         (sum(abs(vect_trial3(:)-distance(katom  ,iatcell,2:4))).lt.tol8).and.&
&         (Sym%indsym(4,isym,watom)==iatcell).and.&
&         (Sym%indsym(4,isym,xatom)==jatcell).and.&
&         (Sym%indsym(4,isym,yatom)==katcell)) then
        if (MPIdata%iam_master) then
          if (inv==1) write(16,'(a,1x,i3)')'===========The bond is kept invariant for isym=',isym
          if (inv==2) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k) --> (i,k,j) for isym=',isym
          if (inv==3) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k) --> (j,i,k) for isym=',isym
          if (inv==4) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k) --> (j,k,i) for isym=',isym
          if (inv==5) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k) --> (k,i,j) for isym=',isym
          if (inv==6) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k) --> (k,j,i) for isym=',isym
        end if
      else
        cycle
      end if
    end if

!   For the 4th order : 24 permutations at all
    if (order==4) then
      vect_trial1(:)=zero
      vect_trial2(:)=zero
      vect_trial3(:)=zero
      vect_trial4(:)=zero
      vect_trial5(:)=zero
      vect_trial6(:)=zero
      if (inv==1) then ; watom=iatcell ; xatom=jatom   ; yatom=katom   ; zatom=latom   ; endif !\Phi4_ijkl
      if (inv==2) then ; watom=iatcell ; xatom=katom   ; yatom=jatom   ; zatom=latom   ; endif !\Phi4_ikjl
      if (inv==3) then ; watom=jatom   ; xatom=iatcell ; yatom=katom   ; zatom=latom   ; endif !\Phi4_jikl
      if (inv==4) then ; watom=jatom   ; xatom=katom   ; yatom=iatcell ; zatom=latom   ; endif !\Phi4_jkil
      if (inv==5) then ; watom=katom   ; xatom=iatcell ; yatom=jatom   ; zatom=latom   ; endif !\Phi4_kijl
      if (inv==6) then ; watom=katom   ; xatom=jatom   ; yatom=iatcell ; zatom=latom   ; endif !\Phi4_kjil

      if (inv==7 ) then ; watom=iatcell ; xatom=jatom   ; yatom=latom   ; zatom=katom   ; endif !\Phi4_ijlk
      if (inv==8 ) then ; watom=iatcell ; xatom=katom   ; yatom=latom   ; zatom=jatom   ; endif !\Phi4_iklj
      if (inv==9 ) then ; watom=jatom   ; xatom=iatcell ; yatom=latom   ; zatom=katom   ; endif !\Phi4_jilk
      if (inv==10) then ; watom=jatom   ; xatom=katom   ; yatom=latom   ; zatom=iatcell ; endif !\Phi4_jkli
      if (inv==11) then ; watom=katom   ; xatom=iatcell ; yatom=latom   ; zatom=jatom   ; endif !\Phi4_kilj
      if (inv==12) then ; watom=katom   ; xatom=jatom   ; yatom=latom   ; zatom=iatcell ; endif !\Phi4_kjli

      if (inv==13) then ; watom=iatcell ; xatom=latom   ; yatom=jatom   ; zatom=katom   ; endif !\Phi4_iljk
      if (inv==14) then ; watom=iatcell ; xatom=latom   ; yatom=katom   ; zatom=jatom   ; endif !\Phi4_ilkj
      if (inv==15) then ; watom=jatom   ; xatom=latom   ; yatom=iatcell ; zatom=katom   ; endif !\Phi4_jlik
      if (inv==16) then ; watom=jatom   ; xatom=latom   ; yatom=katom   ; zatom=iatcell ; endif !\Phi4_jlki
      if (inv==17) then ; watom=katom   ; xatom=latom   ; yatom=iatcell ; zatom=jatom   ; endif !\Phi4_klij
      if (inv==18) then ; watom=katom   ; xatom=latom   ; yatom=jatom   ; zatom=iatcell ; endif !\Phi4_klji

      if (inv==19) then ; watom=latom   ; xatom=iatcell ; yatom=jatom   ; zatom=katom   ; endif !\Phi4_lijk
      if (inv==20) then ; watom=latom   ; xatom=iatcell ; yatom=katom   ; zatom=jatom   ; endif !\Phi4_likj
      if (inv==21) then ; watom=latom   ; xatom=jatom   ; yatom=iatcell ; zatom=katom   ; endif !\Phi4_ljik
      if (inv==22) then ; watom=latom   ; xatom=jatom   ; yatom=katom   ; zatom=iatcell ; endif !\Phi4_ljki
      if (inv==23) then ; watom=latom   ; xatom=katom   ; yatom=iatcell ; zatom=jatom   ; endif !\Phi4_lkij
      if (inv==24) then ; watom=latom   ; xatom=katom   ; yatom=jatom   ; zatom=iatcell ; endif !\Phi4_lkji

      do ii=1,3
        do jj=1,3
          vect_trial1(ii)=vect_trial1(ii)+Sym%S_ref(ii,jj,isym,1)*distance(watom,xatom,jj+1)
          vect_trial2(ii)=vect_trial2(ii)+Sym%S_ref(ii,jj,isym,1)*distance(watom,yatom,jj+1)
          vect_trial3(ii)=vect_trial3(ii)+Sym%S_ref(ii,jj,isym,1)*distance(watom,zatom,jj+1)
          vect_trial4(ii)=vect_trial4(ii)+Sym%S_ref(ii,jj,isym,1)*distance(xatom,yatom,jj+1)
          vect_trial5(ii)=vect_trial5(ii)+Sym%S_ref(ii,jj,isym,1)*distance(xatom,zatom,jj+1)
          vect_trial6(ii)=vect_trial6(ii)+Sym%S_ref(ii,jj,isym,1)*distance(yatom,zatom,jj+1)
        end do
      end do
      jatcell=mod(jatom-1,Invar%natom_unitcell)+1
      katcell=mod(katom-1,Invar%natom_unitcell)+1
      latcell=mod(latom-1,Invar%natom_unitcell)+1
      if ((sum(abs(vect_trial1(:)-distance(iatcell,jatom,2:4))).lt.tol8).and.&
&         (sum(abs(vect_trial2(:)-distance(iatcell,katom,2:4))).lt.tol8).and.&
&         (sum(abs(vect_trial3(:)-distance(iatcell,latom,2:4))).lt.tol8).and.&
&         (sum(abs(vect_trial4(:)-distance(jatom  ,katom,2:4))).lt.tol8).and.&
&         (sum(abs(vect_trial5(:)-distance(jatom  ,latom,2:4))).lt.tol8).and.&
&         (sum(abs(vect_trial6(:)-distance(katom  ,latom,2:4))).lt.tol8).and.&
&         (Sym%indsym(4,isym,watom)==iatcell).and.&
&         (Sym%indsym(4,isym,xatom)==jatcell).and.&
&         (Sym%indsym(4,isym,yatom)==katcell).and.&
&         (Sym%indsym(4,isym,zatom)==latcell)) then
        if (MPIdata%iam_master) then
          if (inv==1 ) write(16,'(a,1x,i3)')'===========The bond is kept invariant for isym=',isym
          if (inv==2 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (i,k,j,l) for isym=',isym !\Phi4_ikjl
          if (inv==3 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (j,i,k,l) for isym=',isym !\Phi4_jikl
          if (inv==4 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (j,k,i,l) for isym=',isym !\Phi4_jkil
          if (inv==5 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (k,i,j,l) for isym=',isym !\Phi4_kijl
          if (inv==6 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (k,j,i,l) for isym=',isym !\Phi4_kjil

          if (inv==7 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (i,j,l,k) for isym=',isym !\Phi4_ijlk
          if (inv==8 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (i,k,l,j) for isym=',isym !\Phi4_iklj
          if (inv==9 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (j,i,l,k) for isym=',isym !\Phi4_jilk
          if (inv==10) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (j,k,l,i) for isym=',isym !\Phi4_jkli
          if (inv==11) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (k,i,l,j) for isym=',isym !\Phi4_kilj
          if (inv==12) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (k,j,l,i) for isym=',isym !\Phi4_kjli

          if (inv==13) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (i,l,j,k) for isym=',isym !\Phi4_iljk
          if (inv==14) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (i,l,k,j) for isym=',isym !\Phi4_ilkj
          if (inv==15) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (j,l,i,k) for isym=',isym !\Phi4_jlik
          if (inv==16) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (j,l,k,i) for isym=',isym !\Phi4_jlki
          if (inv==17) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (k,l,i,j) for isym=',isym !\Phi4_klij
          if (inv==18) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (k,l,j,i) for isym=',isym !\Phi4_klji

          if (inv==19) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (l,i,j,k) for isym=',isym !\Phi4_lijk
          if (inv==20) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (l,i,k,j) for isym=',isym !\Phi4_likj
          if (inv==21) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (l,j,i,k) for isym=',isym !\Phi4_ljik
          if (inv==22) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (l,j,k,i) for isym=',isym !\Phi4_ljki
          if (inv==23) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (l,k,i,j) for isym=',isym !\Phi4_lkij
          if (inv==24) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (l,k,j,i) for isym=',isym !\Phi4_lkji
        end if
      else
        cycle
      end if
    end if

!   Write the S_ref matrix
!FB    write(16,'(3(f16.12,1x))') Sym%S_ref(1,1,isym,1),Sym%S_ref(1,2,isym,1),Sym%S_ref(1,3,isym,1)
!FB    write(16,'(3(f16.12,1x))') Sym%S_ref(2,1,isym,1),Sym%S_ref(2,2,isym,1),Sym%S_ref(2,3,isym,1)
!FB    write(16,'(3(f16.12,1x))') Sym%S_ref(3,1,isym,1),Sym%S_ref(3,2,isym,1),Sym%S_ref(3,3,isym,1)

!   Diagonalize the S_ref matrix
    do ii=1,3
      do jj=1,3
        eigvec(ii,jj)=Sym%S_ref(jj,ii,isym,1)
      end do
    end do
    LWORK=4*3
    ABI_MALLOC(WORK,(LWORK)); WORK(:)=zero
!   This one is real and could be non-symmetric
    call dgeev( 'N', 'V', 3, eigvec, 3, WR, WI, VL, 3, VR, 3, WORK, LWORK, INFO)
    ABI_FREE(WORK)

!   Build the real and imaginary parts of the eigenvectors and eigenvalues
    jj=0
    do ii=1,3
      eigenvalues(ii)=dcmplx(WR(ii),-WI(ii))
      if (WI(ii).ne.zero.and.jj==0) then
        do kk=1,3
          eigenvectors(kk,ii)=dcmplx(VR(kk,ii),VR(kk,ii+1))
        end do
        jj=jj+1
      else if (WI(ii).ne.zero.and.jj==1) then
        do kk=1,3
          eigenvectors(kk,ii)=dcmplx(VR(kk,ii-1),-VR(kk,ii))
        end do
        jj=jj+1
      else
        do kk=1,3
          eigenvectors(kk,ii)=dcmplx(VR(kk,ii),zero)
        end do
      end if
    end do

!   Write the eigenvalues and eigenvectors
    ok=.true.
    do ii=1,3
      if ((aimag(eigenvalues(1)).ne.0).or.(aimag(eigenvalues(2)).ne.0).or.(aimag(eigenvalues(3)).ne.0)) then
        ok=.false.
      end if
    end do
    if (.not.ok.and.MPIdata%iam_master) write(16,'(a)') '            WARNING: THERE IS COMPLEX EIGENVALUES'

!   If the transformation matrix keeps the bond invariant:
!       Phi_{\alpha\beta}=\sum_{\mu\nu} S_{\alpha\mu}.S_{\beta\nu}.Phi_{\mu\nu}
!       If lambda and p are the eigenvectors and eigenvalues of the S matrix, then:
!       \sum_{\alpha\beta} p_{\alpha}^l.p_{\beta}^k Phi_{\alpha\beta}
!     = \sum_{\mu\nu,\alpha\beta} p_{\alpha}^l.p_{\beta}^k.S_{\alpha\mu}.S_{\beta\nu}.Phi_{\mu\nu}
!     = lambda^{*l}.lambda^{*k} \sum_{\mu\nu} p_{\mu}^l.p_{\nu}^k.Phi_{\mu\nu}
!   So, if lambda^{*l}.lambda^{*k} = -1, we must have:
!      \sum_{\alpha\beta} p_{\alpha}^l.p_{\beta}^k.Phi_{\alpha\beta}= 0
!
!   In the case of the reversed bond, one obtains the following constraint:
!      \sum_{\alpha\beta} (lambda^{*l}.lambda^{*k}.p_{\alpha}^l.p_{\beta}^k-p_{\beta}^l.p_{\alpha}^k).Phi_{\alpha\beta}= 0
!   which applies whether lambda^{*l}.lambda^{*k} = \pm 1
!
!   We obtain n vectors with norder coefficients (defined in the R^norder space).
!   The space of the independent solutions are in the R^(norder-n) space, orthogonal
!   to the space spanned by the starting n vectors.
    if (order==1) then
      do ii=1,3
        lambda=eigenvalues(ii)
        if ((abs(real(lambda)-1.d0).lt.tol6).and.(abs(aimag(lambda)).lt.tol6)) cycle
        unchanged(isyminv)=.true.
        iconst(isyminv)=iconst(isyminv)+1
!FB        const_tot=const_tot+1
!FB        write(16,*)'  The eigenvalue',ii
!FB        write(16,*)'  is equal to ',lambda
        do mu=1,3
          alphaij(isyminv,mu,iconst(isyminv))=eigenvectors(mu,ii)
        end do
!FB        write(16,*)'  Real & imaginary parts of the eigenvectors product:'
!FB        write(16,'(3(f16.12,1x))')  real(alphaij(isyminv,:,iconst(isyminv)))
!FB        write(16,'(3(f16.12,1x))') aimag(alphaij(isyminv,:,iconst(isyminv)))
      end do !ii
    else if (order==2) then
      do ii=1,3
        do jj=1,3
          lambda=eigenvalues(ii)*eigenvalues(jj)
          if (((abs(real(lambda)-1.d0).lt.tol6).and.(abs(aimag(lambda)).lt.tol6).and.(inv==1)).or.&
&             ((abs(real(lambda)-1.d0).lt.tol6).and.(abs(aimag(lambda)).lt.tol6).and.(inv==2).and.(ii==jj))) cycle
          unchanged(isyminv)=.true.
          iconst(isyminv)=iconst(isyminv)+1
!FB          const_tot=const_tot+1
!FB          write(16,*)'  The product of eigenvalues',ii,jj
!FB          write(16,*)'  is equal to ',lambda
          do mu=1,3
            do nu=1,3
              pp(mu,nu)=eigenvectors(mu,ii)*eigenvectors(nu,jj)
            end do
          end do
          do mu=1,3
            do nu=1,3
              if (inv==1) then
                alphaij(isyminv,(mu-1)*3+nu,iconst(isyminv))=pp(mu,nu)
              else if (inv==2) then
                alphaij(isyminv,(mu-1)*3+nu,iconst(isyminv))=lambda*pp(mu,nu)-pp(nu,mu)
              else
                ABI_BUG('This symetry is neither Keptinvariant nor Reversed')
              end if
            end do
          end do
!FB          write(16,*)'  Real & imaginary parts of the eigenvectors product:'
!FB          write(16,'(9(f16.12,1x))')  real(alphaij(isyminv,:,iconst(isyminv)))
!FB          write(16,'(9(f16.12,1x))') aimag(alphaij(isyminv,:,iconst(isyminv)))
        end do !jj
      end do !ii
    else if (order==3) then
      do ii=1,3
        do jj=1,3
          do kk=1,3
            lambda=eigenvalues(ii)*eigenvalues(jj)*eigenvalues(kk)
            if (((abs(real(lambda)-1.d0).lt.tol6).and.(abs(aimag(lambda)).lt.tol6).and.(inv==1)).or.&
&               ((abs(real(lambda)-1.d0).lt.tol6).and.(abs(aimag(lambda)).lt.tol6).and.(inv==2).and.(jj==kk)).or.&
&               ((abs(real(lambda)-1.d0).lt.tol6).and.(abs(aimag(lambda)).lt.tol6).and.(inv==3).and.(ii==jj)).or.&
&               ((abs(real(lambda)-1.d0).lt.tol6).and.(abs(aimag(lambda)).lt.tol6).and.(inv==4).and.(ii==jj).and.(jj==kk)).or.&
&               ((abs(real(lambda)-1.d0).lt.tol6).and.(abs(aimag(lambda)).lt.tol6).and.(inv==5).and.(ii==jj).and.(jj==kk)).or.&
&               ((abs(real(lambda)-1.d0).lt.tol6).and.(abs(aimag(lambda)).lt.tol6).and.(inv==6).and.(ii==kk))) cycle
            unchanged(isyminv)=.true.
            iconst(isyminv)=iconst(isyminv)+1
!FB            const_tot=const_tot+1
!FB            write(16,*)'  The product of eigenvalues',ii,jj
!FB            write(16,*)'  is equal to ',lambda
            do mu=1,3
              do nu=1,3
                do xi=1,3
                  ppp(mu,nu,xi)=eigenvectors(mu,ii)*eigenvectors(nu,jj)*eigenvectors(xi,kk)
                end do !xi
              end do !nu
            end do !mu
            do mu=1,3
              do nu=1,3
                do xi=1,3
                  if (inv==1) then
                    alphaij(isyminv,(mu-1)*9+(nu-1)*3+xi,iconst(isyminv))=ppp(mu,nu,xi)
                  else if (inv==2) then
                    alphaij(isyminv,(mu-1)*9+(nu-1)*3+xi,iconst(isyminv))=lambda*ppp(mu,nu,xi)-ppp(mu,xi,nu)
                  else if (inv==3) then
                    alphaij(isyminv,(mu-1)*9+(nu-1)*3+xi,iconst(isyminv))=lambda*ppp(mu,nu,xi)-ppp(nu,mu,xi)
                  else if (inv==4) then
                    alphaij(isyminv,(mu-1)*9+(nu-1)*3+xi,iconst(isyminv))=lambda*ppp(mu,nu,xi)-ppp(nu,xi,mu)
                  else if (inv==5) then
                    alphaij(isyminv,(mu-1)*9+(nu-1)*3+xi,iconst(isyminv))=lambda*ppp(mu,nu,xi)-ppp(xi,mu,nu)
                  else if (inv==6) then
                    alphaij(isyminv,(mu-1)*9+(nu-1)*3+xi,iconst(isyminv))=lambda*ppp(mu,nu,xi)-ppp(xi,nu,mu)
                  else
                    ABI_BUG('This symetry is neither Keptinvariant nor Reversed')
                  end if
                end do !xi
              end do !nu
            end do !mu
!FB            write(16,*)'  Real & imaginary parts of the eigenvectors product:'
!FB            write(16,'(27(f16.12,1x))')  real(alphaij(isyminv,:,iconst(isyminv)))
!FB            write(16,'(27(f16.12,1x))') aimag(alphaij(isyminv,:,iconst(isyminv)))
          end do !kk
        end do !jj
      end do !ii
    else if (order==4) then
      do ii=1,3
        do jj=1,3
          do kk=1,3
            do ll=1,3
              lambda=eigenvalues(ii)*eigenvalues(jj)*eigenvalues(kk)*eigenvalues(ll)
              if ((abs(real(lambda)-1.d0).lt.tol6).and.(abs(aimag(lambda)).lt.tol6)) then
                if ((inv==1 )                                       .or.& !\Phi4_ijkl
&                 ((inv==2 ).and.(jj==kk))                          .or.& !\Phi4_ikjl
&                 ((inv==3 ).and.(ii==jj))                          .or.& !\Phi4_jikl
&                 ((inv==4 ).and.(ii==jj).and.(jj==kk))             .or.& !\Phi4_jkil
&                 ((inv==5 ).and.(ii==jj).and.(jj==kk))             .or.& !\Phi4_kijl
&                 ((inv==6 ).and.(ii==kk))                          .or.& !\Phi4_kjil

&                 ((inv==7 ).and.(kk==ll))                          .or.& !\Phi4_ijlk
&                 ((inv==8 ).and.(jj==kk).and.(kk==ll))             .or.& !\Phi4_iklj
&                 ((inv==9 ).and.(ii==jj).and.(kk==ll))             .or.& !\Phi4_jilk
&                 ((inv==10).and.(ii==jj).and.(jj==kk).and.(kk==ll)).or.& !\Phi4_jkli
&                 ((inv==11).and.(ii==jj).and.(jj==kk).and.(kk==ll)).or.& !\Phi4_kilj
&                 ((inv==12).and.(ii==kk).and.(kk==ll))             .or.& !\Phi4_kjli

&                 ((inv==13).and.(jj==kk).and.(kk==ll))             .or.& !\Phi4_iljk
&                 ((inv==14).and.(jj==ll))                          .or.& !\Phi4_ilkj
&                 ((inv==15).and.(ii==jj).and.(jj==kk).and.(kk==ll)).or.& !\Phi4_jlik
&                 ((inv==16).and.(ii==jj).and.(jj==ll))             .or.& !\Phi4_jlki
&                 ((inv==17).and.(ii==kk).and.(jj==ll))             .or.& !\Phi4_klij
&                 ((inv==18).and.(ii==jj).and.(jj==kk).and.(kk==ll)).or.& !\Phi4_klji

&                 ((inv==19).and.(ii==jj).and.(jj==kk).and.(kk==ll)).or.& !\Phi4_lijk
&                 ((inv==20).and.(ii==jj).and.(jj==ll))             .or.& !\Phi4_likj
&                 ((inv==21).and.(ii==kk).and.(kk==ll))             .or.& !\Phi4_ljik
&                 ((inv==22).and.(ii==ll))                          .or.& !\Phi4_ljki
&                 ((inv==23).and.(ii==jj).and.(jj==kk).and.(kk==ll)).or.& !\Phi4_lkij
&                 ((inv==24).and.(ii==ll).and.(jj==kk))) cycle            !\Phi4_lkji
              end if
              unchanged(isyminv)=.true.
              iconst(isyminv)=iconst(isyminv)+1
!FB              const_tot=const_tot+1
!FB              write(16,*)'  The product of eigenvalues',ii,jj
!FB              write(16,*)'  is equal to ',lambda
              do mu=1,3
                do nu=1,3
                  do xi=1,3
                    do zeta=1,3
                      pppp(mu,nu,xi,zeta)=eigenvectors(mu,ii)*eigenvectors(nu,jj)*eigenvectors(xi,kk)*eigenvectors(zeta,ll)
                    end do !zeta
                  end do !xi
                end do !nu
              end do !mu
              do mu=1,3
                do nu=1,3
                  do xi=1,3
                    do zeta=1,3
                      itemp=(mu-1)*27+(nu-1)*9+(xi-1)*3+zeta
                      if (inv==1)       then ; alphaij(isyminv,itemp,iconst(isyminv))=pppp(mu,nu,xi,zeta)
                      else if (inv==2 ) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(mu,xi,nu,zeta)
                      else if (inv==3 ) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(nu,mu,xi,zeta)
                      else if (inv==4 ) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(nu,xi,mu,zeta)
                      else if (inv==5 ) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(xi,mu,nu,zeta)
                      else if (inv==6 ) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(xi,nu,mu,zeta)

                      else if (inv==7 ) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(mu,nu,zeta,xi)
                      else if (inv==8 ) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(mu,xi,zeta,nu)
                      else if (inv==9 ) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(nu,mu,zeta,xi)
                      else if (inv==10) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(nu,xi,zeta,mu)
                      else if (inv==11) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(xi,mu,zeta,nu)
                      else if (inv==12) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(xi,nu,zeta,mu)

                      else if (inv==13) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(mu,zeta,nu,xi)
                      else if (inv==14) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(mu,zeta,xi,nu)
                      else if (inv==15) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(nu,zeta,mu,xi)
                      else if (inv==16) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(nu,zeta,xi,mu)
                      else if (inv==17) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(xi,zeta,mu,nu)
                      else if (inv==18) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(xi,zeta,nu,mu)

                      else if (inv==19) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(zeta,mu,nu,xi)
                      else if (inv==20) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(zeta,mu,xi,nu)
                      else if (inv==21) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(zeta,nu,mu,xi)
                      else if (inv==22) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(zeta,nu,xi,mu)
                      else if (inv==23) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(zeta,xi,mu,nu)
                      else if (inv==24) then ; alphaij(isyminv,itemp,iconst(isyminv))=lambda*pppp(mu,nu,xi,zeta)-pppp(zeta,xi,nu,mu)
                      else ; ABI_BUG('This symetry is neither Keptinvariant nor Reversed')
                      end if
                    end do !zeta
                  end do !xi
                end do !nu
              end do !mu
!FB              write(16,*)'  Real & imaginary parts of the eigenvectors product:'
!FB              write(16,'(81(f16.12,1x))')  real(alphaij(isyminv,:,iconst(isyminv)))
!FB              write(16,'(81(f16.12,1x))') aimag(alphaij(isyminv,:,iconst(isyminv)))
            end do !ll
          end do !kk
        end do !jj
      end do !ii
    else
      ABI_BUG('Only the first, second, third and fourth order are allowed')
    end if


!FB=================================================================
!FB======== TO CLEAN ===============================================
!FB=================================================================
  nconst_loc=const_tot+iconst(isyminv)
  ii=0
  ABI_MALLOC(tab_vec,(norder,nconst_loc)); tab_vec(:,:)=czero
  do itemp=1,isyminv
    if (unchanged(itemp)) then
      do jj=1,iconst(itemp)
        ii=ii+1
        tab_vec(:,ii)=alphaij(itemp,:,jj)
      end do
    end if
  end do

  do kk=2,nconst_loc
    do jj=1,kk-1
      prod_scal=sum( real(tab_vec(:,jj))* real(tab_vec(:,jj))+aimag(tab_vec(:,jj))*aimag(tab_vec(:,jj)))
      if (abs(prod_scal).gt.tol8) then
        tab_vec(:,kk)=tab_vec(:,kk)-sum(tab_vec(:,kk)*conjg(tab_vec(:,jj)))/dcmplx(prod_scal,zero)*tab_vec(:,jj)
        do ii=1,norder
          if (abs( real(tab_vec(ii,kk))).lt.tol8) tab_vec(ii,kk)=dcmplx(zero,aimag(tab_vec(ii,kk)))
          if (abs(aimag(tab_vec(ii,kk))).lt.tol8) tab_vec(ii,kk)=dcmplx( real(tab_vec(ii,kk)),zero)
        end do
      end if
    end do
  end do

! On stocke les vecteurs non-nuls
  ABI_MALLOC(temp   ,(norder,nconst_loc)); temp(:,:)   =czero
  ii=0
  do kk=1,nconst_loc
    prod_scal=sum( real(tab_vec(:,kk))* real(tab_vec(:,kk))+aimag(tab_vec(:,kk))*aimag(tab_vec(:,kk)))
    if (abs(prod_scal).gt.tol8) then
      ii=ii+1
      temp(:,ii)=tab_vec(:,kk)/dsqrt(prod_scal)
    end if
  end do
  ABI_FREE(tab_vec)
  iconst(isyminv)=ii-const_tot
  const_tot=const_tot+iconst(isyminv)

  ii=0
  alphaij(:,:,:)=czero
  do itemp=1,isyminv
    if (unchanged(itemp)) then
      do jj=1,iconst(itemp)
        ii=ii+1
        alphaij(itemp,:,jj)=temp(:,ii)
      end do
    end if
  end do
  ABI_FREE(temp)
!FB=================================================================
!FB======== TO CLEAN ===============================================
!FB=================================================================






!   WARNING: There are some minimum and maximum of constraints
    if (order==1.and.(iconst(isyminv).eq.3)) then
      ncoeff=0
      proj(:,:,ishell)=zero
      ABI_FREE(unchanged)
      ABI_FREE(alphaij)
      ABI_FREE(iconst)
      return
    else if (order==1.and.(iconst(isyminv).gt.3)) then
      ABI_BUG(' First order : There are more than 3 constraints')
    end if
    if (order==2.and.(iconst(isyminv).gt.8)) then
      ABI_BUG(' Second order : There are more than 8 constraints')
    end if
    if (order==3.and.(iconst(isyminv).gt.27)) then
      ABI_BUG(' Third order : There are more than 27 constraints')
    end if
    if (order==4.and.(iconst(isyminv).gt.81)) then
      ABI_BUG(' Fourth order : There are more than 81 constraints')
    end if
  end do !isyminv
! ================================================================================================
! =========== End big loop over symetries and facorder ===========================================
! ================================================================================================
  nconst_perm=0
! The (iik, iji, ijj and iii) third order IFCs are symmetric with respect to some permutations.
! Some constraints have to be added :
  if (order.eq.3) then
    if ((iatcell.eq.jatom).or.(iatcell.eq.katom).or.(jatom.eq.katom)) then
      nconst_perm=5
      if (MPIdata%iam_master) write(16,'(a)')'=========== The IFCs are symmetric'
      const_tot=const_tot+nconst_perm*norder
      ABI_MALLOC(constraints,(nconst_perm,norder,norder)) ; constraints(:,:,:)=czero
      ii=0
      do mu=1,3
        do nu=1,3
          do xi=1,3
            ii=ii+1
            if (iatcell.eq.jatom) then
              if (mu.eq.nu) cycle
              constraints(1,(mu-1)*9+(nu-1)*3+xi,ii)= cone
              constraints(1,(nu-1)*9+(mu-1)*3+xi,ii)=-cone
            end if
            if (iatcell.eq.katom) then
              if (mu.eq.xi) cycle
              constraints(2,(mu-1)*9+(nu-1)*3+xi,ii)= cone
              constraints(2,(xi-1)*9+(nu-1)*3+mu,ii)=-cone
            end if
            if (jatom.eq.katom) then
              if (nu.eq.xi) cycle
              constraints(3,(mu-1)*9+(nu-1)*3+xi,ii)= cone
              constraints(3,(mu-1)*9+(xi-1)*3+nu,ii)=-cone
            end if
            if ((iatcell.eq.jatom).and.(jatom.eq.katom)) then
              if ((nu.eq.xi).and.(nu.eq.mu)) cycle
              constraints(4,(mu-1)*9+(nu-1)*3+xi,ii)= cone
              constraints(4,(xi-1)*9+(mu-1)*3+nu,ii)=-cone
            end if
            if ((iatcell.eq.jatom).and.(jatom.eq.katom)) then
              if ((nu.eq.xi).and.(nu.eq.mu)) cycle
              constraints(5,(mu-1)*9+(nu-1)*3+xi,ii)= cone
              constraints(5,(nu-1)*9+(xi-1)*3+mu,ii)=-cone
            end if
          end do
        end do
      end do
    end if
  end if

! The (iikl, ijil, ijki, ijjl, ijkj, ijkk, iiil, iiki, ijii, ijjj, iiii)
! fourth order IFCs are symmetric with respect to some permutations.
! Some constraints have to be added :
  if (order.eq.4) then
    if ((iatcell.eq.jatom).or.(iatcell.eq.katom).or.(iatcell.eq.latom)&
&                         .or.(jatom.eq.katom).or.(jatom.eq.latom).or.(katom.eq.latom)) then
      nconst_perm=17
      if (MPIdata%iam_master) write(16,'(a)')'=========== The IFCs are symmetric'
      const_tot=const_tot+nconst_perm*norder
      ABI_MALLOC(constraints,(nconst_perm,norder,norder)) ; constraints(:,:,:)=czero
      ii=0
      do mu=1,3
        do nu=1,3
          do xi=1,3
            do zeta=1,3
              ii=ii+1
              if (iatcell.eq.jatom) then
                if (mu.eq.nu) cycle
                constraints(1,(mu-1)*27+(nu-1)*9+(xi-1)*3+zeta,ii)= cone
                constraints(1,(nu-1)*27+(mu-1)*9+(xi-1)*3+zeta,ii)=-cone
              end if
              if (iatcell.eq.katom) then
                if (mu.eq.xi) cycle
                constraints(2,(mu-1)*27+(nu-1)*9+(xi-1)*3+zeta,ii)= cone
                constraints(2,(xi-1)*27+(nu-1)*9+(mu-1)*3+zeta,ii)=-cone
              end if
              if (iatcell.eq.latom) then
                if (mu.eq.zeta) cycle
                constraints(3,(mu  -1)*27+(nu-1)*9+(xi-1)*3+zeta,ii)= cone
                constraints(3,(zeta-1)*27+(nu-1)*9+(xi-1)*3+mu  ,ii)=-cone
              end if
              if (jatom.eq.katom) then
                if (nu.eq.xi) cycle
                constraints(4,(mu-1)*27+(nu-1)*9+(xi-1)*3+zeta,ii)= cone
                constraints(4,(mu-1)*27+(xi-1)*9+(nu-1)*3+zeta,ii)=-cone
              end if
              if (jatom.eq.latom) then
                if (nu.eq.zeta) cycle
                constraints(5,(mu-1)*27+(nu  -1)*9+(xi-1)*3+zeta,ii)= cone
                constraints(5,(mu-1)*27+(zeta-1)*9+(xi-1)*3+nu  ,ii)=-cone
              end if
              if (katom.eq.latom) then
                if (xi.eq.zeta) cycle
                constraints(6,(mu-1)*27+(nu-1)*9+(xi  -1)*3+zeta,ii)= cone
                constraints(6,(mu-1)*27+(nu-1)*9+(zeta-1)*3+xi  ,ii)=-cone
              end if

              if ((iatcell.eq.jatom).and.(jatom.eq.katom)) then
                if ((mu.eq.nu).and.(nu.eq.xi)) cycle
                constraints(7,(mu-1)*27+(nu-1)*9+(xi-1)*3+zeta,ii)= cone
                constraints(7,(xi-1)*27+(mu-1)*9+(nu-1)*3+zeta,ii)=-cone
              end if
              if ((iatcell.eq.jatom).and.(jatom.eq.katom)) then
                if ((mu.eq.nu).and.(nu.eq.xi)) cycle
                constraints(8,(mu-1)*27+(nu-1)*9+(xi-1)*3+zeta,ii)= cone
                constraints(8,(nu-1)*27+(xi-1)*9+(mu-1)*3+zeta,ii)=-cone
              end if

              if ((iatcell.eq.jatom).and.(jatom.eq.latom)) then
                if ((mu.eq.nu).and.(nu.eq.zeta)) cycle
                constraints(9,(mu-1)*27+(nu  -1)*9+(xi-1)*3+zeta,ii)= cone
                constraints(9,(nu-1)*27+(zeta-1)*9+(xi-1)*3+mu  ,ii)=-cone
              end if
              if ((iatcell.eq.jatom).and.(jatom.eq.latom)) then
                if ((mu.eq.nu).and.(nu.eq.zeta)) cycle
                constraints(10,(mu  -1)*27+(nu-1)*9+(xi-1)*3+zeta,ii)= cone
                constraints(10,(zeta-1)*27+(mu-1)*9+(xi-1)*3+nu  ,ii)=-cone
              end if

              if ((iatcell.eq.katom).and.(katom.eq.latom)) then
                if ((mu.eq.xi).and.(xi.eq.zeta)) cycle
                constraints(11,(mu-1)*27+(nu-1)*9+(xi  -1)*3+zeta,ii)= cone
                constraints(11,(xi-1)*27+(nu-1)*9+(zeta-1)*3+mu  ,ii)=-cone
              end if
              if ((iatcell.eq.katom).and.(katom.eq.latom)) then
                if ((mu.eq.xi).and.(xi.eq.zeta)) cycle
                constraints(12,(mu  -1)*27+(nu-1)*9+(xi-1)*3+zeta,ii)= cone
                constraints(12,(zeta-1)*27+(nu-1)*9+(mu-1)*3+xi  ,ii)=-cone
              end if

              if ((jatom.eq.katom).and.(katom.eq.latom)) then
                if ((nu.eq.xi).and.(xi.eq.zeta)) cycle
                constraints(13,(mu-1)*27+(nu-1)*9+(xi  -1)*3+zeta,ii)= cone
                constraints(13,(mu-1)*27+(xi-1)*9+(zeta-1)*3+nu  ,ii)=-cone
              end if
              if ((jatom.eq.katom).and.(katom.eq.latom)) then
                if ((nu.eq.xi).and.(xi.eq.zeta)) cycle
                constraints(14,(mu-1)*27+(nu  -1)*9+(xi-1)*3+zeta,ii)= cone
                constraints(14,(mu-1)*27+(zeta-1)*9+(nu-1)*3+xi  ,ii)=-cone
              end if

              if ((iatcell.eq.jatom).and.(jatom.eq.katom).and.(katom.eq.latom)) then
                if ((mu.eq.nu).and.(nu.eq.xi).and.(xi.eq.zeta)) cycle
                constraints(15,(mu-1)*27+(nu-1)*9+(xi  -1)*3+zeta,ii)= cone
                constraints(15,(nu-1)*27+(xi-1)*9+(zeta-1)*3+mu  ,ii)=-cone
              end if
              if ((iatcell.eq.jatom).and.(jatom.eq.katom).and.(katom.eq.latom)) then
                if ((mu.eq.nu).and.(nu.eq.xi).and.(xi.eq.zeta)) cycle
                constraints(16,(mu-1)*27+(nu  -1)*9+(xi-1)*3+zeta,ii)= cone
                constraints(16,(xi-1)*27+(zeta-1)*9+(mu-1)*3+nu  ,ii)=-cone
              end if
              if ((iatcell.eq.jatom).and.(jatom.eq.katom).and.(katom.eq.latom)) then
                if ((mu.eq.nu).and.(nu.eq.xi).and.(xi.eq.zeta)) cycle
                constraints(17,(mu  -1)*27+(nu-1)*9+(xi-1)*3+zeta,ii)= cone
                constraints(17,(zeta-1)*27+(mu-1)*9+(nu-1)*3+xi  ,ii)=-cone
              end if
            end do
          end do
        end do
      end do
    end if
  end if

! In the case where the matrix has norder**2 inequivalent and non-zero elements
  if (const_tot==0) then
    write(message,'(a,1x,i3,1x,a)') 'For shell number=',ishell,'there is no symetry operation reducing the number of coefficients'
    ABI_WARNING(message)
    proj(:,:,ishell)=0.d0
    do ii=1,norder
      proj(ii,ii,ishell)=1.d0
    end do
    ncoeff=norder
    ABI_FREE(unchanged)
    ABI_FREE(alphaij)
    ABI_FREE(iconst)
    return
  end if

! When some constraints have been found
  ncount=const_tot
  if (MPIdata%iam_master) then
    write(16,'(a,1x,i7,1x,a)') 'There is a total of ',ncount,' non-independant constraints for this shell'
  end if
  ii=0
  ABI_MALLOC(tab_vec,(norder,ncount)); tab_vec(:,:)=czero
  ABI_MALLOC(temp   ,(norder,ncount)); temp(:,:)   =czero
  do isyminv=1,nsyminv
    if (unchanged(isyminv)) then
      do jj=1,iconst(isyminv)
        ii=ii+1
        tab_vec(:,ii)=alphaij(isyminv,:,jj)
      end do
    end if
  end do
  ABI_FREE(unchanged)
  ABI_FREE(alphaij)
  ABI_FREE(iconst)
! Add the constraints coming from the symmetry of the IFCs (at the 3rd order)
  if (nconst_perm.gt.0) then
    do jj=1,norder
      do kk=1,nconst_perm
        ii=ii+1
        tab_vec(:,ii)=constraints(kk,:,jj)
      end do
    end do
    ABI_FREE(constraints)
  end if
  if (ii.ne.ncount) then
    write(message,'(i7,1x,a,1x,i7)') ii,' non equal to ',ncount
    ABI_BUG(message)
  end if
  do ii=1,norder
    do jj=1,ncount
      if (abs( real(tab_vec(ii,jj))).lt.tol8) tab_vec(ii,jj)=dcmplx(zero,aimag(tab_vec(ii,jj)))
      if (abs(aimag(tab_vec(ii,jj))).lt.tol8) tab_vec(ii,jj)=dcmplx( real(tab_vec(ii,jj)),zero)
    end do
  end do

! On stocke les vecteurs non-nuls
  ii=0
  do kk=1,ncount
    prod_scal=sum( real(tab_vec(:,kk))* real(tab_vec(:,kk))+aimag(tab_vec(:,kk))*aimag(tab_vec(:,kk)))
    if (abs(prod_scal).gt.tol8) then
      ii=ii+1
      temp(:,ii)=tab_vec(:,kk)/dsqrt(prod_scal)
    end if
  end do
  ncount=ii
  ABI_FREE(tab_vec)
  ABI_MALLOC(tab_vec,(norder,ncount)); tab_vec(:,1:ncount)=temp(:,1:ncount)
  ABI_FREE(temp)
  ABI_MALLOC(temp   ,(norder,ncount)); temp(:,:)   =czero

! L'ensemble des vecteurs reduisants l'espace de R^norder a R^n ne forment pas une base
! independante. Il faut donc trouver les vecteurs independants.
! --> Orthogonalisation de Gram-Schmidt
  do kk=2,ncount
    do jj=1,kk-1
      prod_scal=sum( real(tab_vec(:,jj))* real(tab_vec(:,jj))+aimag(tab_vec(:,jj))*aimag(tab_vec(:,jj)))
      if (abs(prod_scal).gt.tol8) then
        tab_vec(:,kk)=tab_vec(:,kk)-sum(tab_vec(:,kk)*conjg(tab_vec(:,jj)))/dcmplx(prod_scal,zero)*tab_vec(:,jj)
        do ii=1,norder
          if (abs( real(tab_vec(ii,kk))).lt.tol8) tab_vec(ii,kk)=dcmplx(zero,aimag(tab_vec(ii,kk)))
          if (abs(aimag(tab_vec(ii,kk))).lt.tol8) tab_vec(ii,kk)=dcmplx( real(tab_vec(ii,kk)),zero)
        end do
!FB      else
!FB        write(Invar%stdout,*)'One prod_scal equals zero'
      end if
    end do
  end do

! On stocke les vecteurs non-nuls
  ii=0
  do kk=1,ncount
    prod_scal=sum( real(tab_vec(:,kk))* real(tab_vec(:,kk))+aimag(tab_vec(:,kk))*aimag(tab_vec(:,kk)))
    if (abs(prod_scal).gt.tol8) then
      ii=ii+1
      temp(:,ii)=tab_vec(:,kk)/dsqrt(prod_scal)
    end if
  end do
  ncount=ii
  ABI_FREE(tab_vec)

! On ecrit les vecteurs non-nuls
!FB  write(16,*) ' '
!FB  write(16,*) '  ========The final set of vectors is:'
!FB  do kk=1,ncount
!FB    write(16,'(81(f16.12,1x))')  real(temp(:,kk))
!FB    write(16,'(81(f16.12,1x))') aimag(temp(:,kk))
!FB  end do
  if (MPIdata%iam_master) then
    write(16,'(a,1x,i7,1x,a)') '  ======= Finally, there are ',ncount,' independent vectors'
  end if
  if (ncount.gt.8.and.order==2) then
    ABI_ERROR(' Order 2 : There are too many independent vectors')
  end if
  if (ncount.gt.27.and.order==3) then
    ABI_ERROR(' Order 3 : There are too many independent vectors')
  end if
  if (ncount.gt.81.and.order==4) then
    ABI_ERROR(' Order 4 : There are too many independent vectors')
  end if

! On cherche les (norder-ncount) vecteurs orthogonaux aux vecteurs non-nuls
! --> Orthogonalisation de Gram-Schmidt
  ABI_MALLOC(tab_vec,(norder,norder)); tab_vec(:,:)=czero
  iseed=-5
  do kk=1,norder
    if (kk.le.ncount) then
      tab_vec(:,kk)=temp(:,kk)
    else
      do jj=1,norder
        drandom=uniformrandom(iseed)
        tab_vec(jj,kk)=dcmplx(drandom,zero)
      end do
      do jj=1,kk-1
        prod_scal=sum( real(tab_vec(:,jj))* real(tab_vec(:,jj))+aimag(tab_vec(:,jj))*aimag(tab_vec(:,jj)))
        if (abs(prod_scal).gt.tol8) then
          tab_vec(:,kk)=tab_vec(:,kk)-sum(tab_vec(:,kk)*conjg(tab_vec(:,jj)))/prod_scal*tab_vec(:,jj)
          do ii=1,norder
            if (abs( real(tab_vec(ii,kk))).lt.tol8) tab_vec(ii,kk)=dcmplx(zero,aimag(tab_vec(ii,kk)))
            if (abs(aimag(tab_vec(ii,kk))).lt.tol8) tab_vec(ii,kk)=dcmplx( real(tab_vec(ii,kk)),zero)
          end do
        end if
        prod_scal=sum( real(tab_vec(:,kk))* real(tab_vec(:,kk))+aimag(tab_vec(:,kk))*aimag(tab_vec(:,kk)))
        tab_vec(:,kk)=tab_vec(:,kk)/dsqrt(prod_scal)
      end do
    end if
  end do
  ABI_FREE(temp)

! On ecrit les vecteurs non-nuls
!FB  write(16,*) ' '
!FB  write(16,*) '  ========The orthogonal set of vectors is:'
  do kk=ncount+1,norder
!FB    write(16,'(81(f16.12,1x))')  real(tab_vec(:,kk))
!FB    write(16,'(81(f16.12,1x))') aimag(tab_vec(:,kk))
    if ((abs(aimag(tab_vec(1,kk))).gt.tol6).or.&
&       (abs(aimag(tab_vec(1,kk))).gt.tol6).or.&
&       (abs(aimag(tab_vec(1,kk))).gt.tol6)) then
      ABI_ERROR('the constraint has an imaginary part')
    end if
  end do
  ncoeff=norder-ncount
  if (MPIdata%iam_master) then
    write(16,'(a,1x,i7,1x,a)') '  ======= Finally, there are ',ncoeff,' coefficients'
  end if

! On copie tab_vec dans proj
  do icoeff=1,ncoeff
    proj(:,icoeff,ishell)=tab_vec(:,ncount+icoeff)
  end do
  ABI_FREE(tab_vec)

end subroutine tdep_calc_nbcoeff

!====================================================================================================

end module m_tdep_shell
