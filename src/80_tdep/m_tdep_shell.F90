
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_shell

  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_tdep_readwrite,   only : Input_Variables_type, MPI_enreg_type
  use m_tdep_latt,        only : Lattice_Variables_type
  use m_tdep_sym,         only : Symetries_Variables_type, tdep_SearchS_2at, tdep_SearchS_3at, tdep_SearchS_4at
  use m_tdep_utils,       only : tdep_calc_nbcoeff
  use m_io_tools

  implicit none

  type List_of_neighbours
    integer :: n_interactions
    integer, allocatable :: atomj_in_shell(:)
    integer, allocatable :: atomk_in_shell(:)
    integer, allocatable :: atoml_in_shell(:)
    integer, allocatable :: sym_in_shell(:)
    integer, allocatable :: transpose_in_shell(:)
  end type List_of_neighbours

  type Shell_Variables_type
    integer :: nshell
    integer, allocatable :: ncoeff(:)
    integer, allocatable :: ncoeff_prev(:)
    integer, allocatable :: ishell_self(:)
    integer, allocatable :: iatref(:)
    integer, allocatable :: jatref(:)
    integer, allocatable :: katref(:)
    integer, allocatable :: latref(:)
    type(List_of_neighbours),allocatable :: neighbours(:,:)
  end type Shell_Variables_type

  public :: tdep_init_shell1at
  public :: tdep_init_shell2at
  public :: tdep_init_shell3at
  public :: tdep_init_shell4at
  public :: tdep_destroy_shell

contains

!====================================================================================================
 subroutine tdep_init_shell1at(distance,InVar,MPIdata,norder,nshell_max,ntotcoeff,order,proj,Shell1at,Sym)

  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(out) :: Shell1at
  type(Symetries_Variables_type),intent(inout) :: Sym
  type(MPI_enreg_type), intent(in) :: MPIdata
  integer,intent(in) :: norder,order,nshell_max
  integer,intent(out) :: ntotcoeff
  double precision, intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision, intent(out) :: proj(norder,norder,nshell_max)
  integer :: ishell,iatcell,iatom,eatom,iatref
  integer :: natom,natom_unitcell,counter,ncoeff,ncoeff_prev,isym
  integer, allocatable :: ref1at(:,:),Isym1at(:,:)

  natom         =InVar%natom
  natom_unitcell=InVar%natom_unitcell

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '####### FIRST ORDER : find the number of coefficients #######################'
  write(InVar%stdout,*) '#############################################################################'

! - Identify the shells
! - Store the index of the atoms included in each shell
! - Store the reference atoms for each shell
! - Compute the symetry operation between the reference atom and another one
  write(InVar%stdout,*) ' Build the ref1at and Isym1at tables...'
  ABI_MALLOC(ref1at ,(natom,2)) ; ref1at (:,:)=zero
  ABI_MALLOC(Isym1at,(natom,1)) ; Isym1at(:,:)=zero
  ishell=0
  do iatcell=1,natom_unitcell
    if (ref1at(iatcell,1).ne.0) cycle
    ishell=ishell+1
    do eatom=1,natom
      if (ref1at(eatom,1).eq.0) then
        do isym=1,Sym%nsym
!FB          write(6,'(4(i5,x))') Sym%indsym(4,isym,eatom),eatom,iatcell,isym
          if (Sym%indsym(4,isym,eatom).eq.iatcell) then
            Isym1at(eatom,1)=isym
            ref1at(eatom,1)=iatcell
            ref1at(eatom,2)=ishell
            if (InVar%debug) write(InVar%stdout,'(a,1x,2(i4,1x),a,i4)') &
&             'For:',iatcell,eatom,' direct transformation with isym=',Isym1at(eatom,1)
            exit
          end if  
        end do !isym  
      end if !already treated
    end do !eatom
  end do !iatcell 
  Shell1at%nshell=ishell
  if (Shell1at%nshell.gt.nshell_max) then
    write(InVar%stdout,*) '  STOP : The maximum number of shells allowed by the code is:',nshell_max
    write(InVar%stdout,*) '         In the present calculation, the number of shells is:',Shell1at%nshell
    write(InVar%stdout,*) '         Action: increase nshell_max'
    MSG_ERROR('The maximum number of shells allowed by the code is reached')
  end if  


! Store all the previous quantities in a better way than in ref1at (without using too memory).
  write(InVar%stdout,*) ' Build the Shell1at datatype...'
  ABI_MALLOC(Shell1at%neighbours,(1,Shell1at%nshell))
  ABI_MALLOC(Shell1at%iatref       ,(Shell1at%nshell)); Shell1at%iatref(:)=zero
  do ishell=1,Shell1at%nshell
    counter=0
    do iatom=1,natom
      if (ref1at(iatom,2).eq.ishell) counter=counter+1 
    end do
    Shell1at%neighbours(1,ishell)%n_interactions=counter
    if (counter.eq.0) cycle
    ABI_MALLOC(Shell1at%neighbours(1,ishell)%atomj_in_shell,(counter))
    ABI_MALLOC(Shell1at%neighbours(1,ishell)%sym_in_shell,(counter))
    Shell1at%neighbours(1,ishell)%atomj_in_shell(:)=zero
    Shell1at%neighbours(1,ishell)%sym_in_shell(:)=zero
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
  ABI_MALLOC(Shell1at%ncoeff     ,(Shell1at%nshell)); Shell1at%ncoeff(:)=zero
  ABI_MALLOC(Shell1at%ncoeff_prev,(Shell1at%nshell)); Shell1at%ncoeff_prev(:)=zero
  write(InVar%stdout,*) ' Number of shells=',Shell1at%nshell
  write(InVar%stdout,*) '============================================================================'
  if (MPIdata%iam_master) open(unit=16,file=trim(InVar%output_prefix)//'nbcoeff-phi1.dat')
  ncoeff_prev=0
  do ishell=1,Shell1at%nshell
    ncoeff=0
    iatref=Shell1at%iatref(ishell)
    write(InVar%stdout,*) 'Shell number:',ishell 
    write(InVar%stdout,'(a,i5,a)') '  For atom',iatref,':'
    call tdep_calc_nbcoeff(distance,iatref,InVar,ishell,1,1,1,MPIdata,ncoeff,norder,Shell1at%nshell,order,proj,Sym)
    if (ncoeff.eq.0) then 
      Shell1at%neighbours(1,ishell)%n_interactions=0
    end if
    Shell1at%ncoeff     (ishell)=ncoeff
    Shell1at%ncoeff_prev(ishell)=ncoeff_prev
    ncoeff_prev=ncoeff_prev+ncoeff
    write(InVar%stdout,*)'  Number of independant coefficients in this shell=',ncoeff
    write(InVar%stdout,*)'  Number of interactions in this shell=',Shell1at%neighbours(1,ishell)%n_interactions
!FB    write(InVar%stdout,*)'  The ratio is=',dfloat(Shell1at%neighbours(iatref,ishell)%n_interactions)/dfloat(ncoeff)
    write(InVar%stdout,*) '============================================================================'
  end do  
  write(InVar%stdout,*)'  >>>>>> Total number of coefficients at the first order=',ncoeff_prev
  if (MPIdata%iam_master) close(16)
  ntotcoeff=ncoeff_prev

 end subroutine tdep_init_shell1at

!====================================================================================================
 subroutine tdep_init_shell2at(distance,InVar,MPIdata,norder,nshell_max,ntotcoeff,order,proj,Shell2at,Sym)

  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(out) :: Shell2at
  type(Symetries_Variables_type),intent(inout) :: Sym
  type(MPI_enreg_type), intent(in) :: MPIdata
  integer,intent(in) :: norder,order,nshell_max
  integer,intent(out) :: ntotcoeff
  double precision, intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision, intent(out) :: proj(norder,norder,nshell_max)

  integer :: ishell,iatcell,iatom,jatom,eatom,fatom,iatref,jatref !,ii
  integer :: natom,natom_unitcell,counter,ncoeff,ncoeff_prev
  integer, allocatable :: ref2at(:,:,:),Isym2at(:,:,:)

  natom         =InVar%natom
  natom_unitcell=InVar%natom_unitcell

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '###### SECOND ORDER : find the number of coefficients #######################'
  write(InVar%stdout,*) '#############################################################################'

! - Identify the shells
! - Store the index of the atoms included in each shell
! - Store the reference atoms for each shell
! - Compute the symetry operation between the reference atom and another one
  write(InVar%stdout,*) ' Build the ref2at and Isym2at tables...'
  ABI_MALLOC(ref2at ,(natom,natom,3)) ; ref2at (:,:,:)=zero
  ABI_MALLOC(Isym2at,(natom,natom,2)) ; Isym2at(:,:,:)=zero
  ABI_MALLOC(Shell2at%ishell_self,(natom_unitcell)) ; Shell2at%ishell_self(:)=zero
  ishell=0
  do iatcell=1,natom_unitcell
    do jatom=1,natom
!     Interactions are only computed until Rcut in order to have complete shell of neighbours.
!     Otherwise the symetries are broken.
      if ((ref2at(iatcell,jatom,1).ne.0).or.(distance(iatcell,jatom,1).gt.(InVar%Rcut*0.99))) cycle
      ishell=ishell+1
      if (iatcell.eq.jatom) Shell2at%ishell_self(iatcell)=ishell
      do eatom=1,natom
        do fatom=1,natom
          if ((ref2at(eatom,fatom,1).eq.0).and.&
!FB&         (abs(distance(iatcell,jatom,1)-distance(eatom,fatom,1)).lt.1.d-3)) then
&         (abs(distance(iatcell,jatom,1)-distance(eatom,fatom,1)).lt.1.d-6)) then
            call tdep_SearchS_2at(InVar,iatcell,jatom,eatom,fatom,Isym2at,Sym,InVar%xred_ideal)
            if (Isym2at(eatom,fatom,2)==1) then
              if (InVar%debug) write(InVar%stdout,'(a,1x,4(i4,1x),a,i4)') &
&                'For:',iatcell,jatom,eatom,fatom,' direct transformation with isym=',Isym2at(eatom,fatom,1)
              ref2at(eatom,fatom,1)=iatcell
              ref2at(eatom,fatom,2)=jatom
              ref2at(eatom,fatom,3)=ishell
!             The Phi2 has to be symetric (transposition symetries)
              if (InVar%debug) write(InVar%stdout,'(a,1x,4(i4,1x),a,i4)') &
&                'For:',iatcell,jatom,eatom,fatom,' transformation+permutation with isym=',Isym2at(eatom,fatom,1)
              ref2at(fatom,eatom,1)=iatcell
              ref2at(fatom,eatom,2)=jatom
              ref2at(fatom,eatom,3)=ishell
              Isym2at(fatom,eatom,1)=Isym2at(eatom,fatom,1)
              Isym2at(fatom,eatom,2)=2
            else
              if (InVar%debug) write(InVar%stdout,'(a,4(1x,i4))') &
&                'NO SYMETRY OPERATION BETWEEN (iatom,jatom) and (eatom,fatom)=',iatcell,jatom,eatom,fatom
            end if
          end if !already treated
        end do !fatom
      end do !eatom
    end do !jatom
  end do !iatcell
  Shell2at%nshell=ishell
  if (Shell2at%nshell.gt.nshell_max) then
    write(InVar%stdout,*) '  STOP : The maximum number of shells allowed by the code is:',nshell_max
    write(InVar%stdout,*) '         In the present calculation, the number of shells is:',Shell2at%nshell
    write(InVar%stdout,*) '         Action: increase nshell_max'
    MSG_ERROR('The maximum number of shells allowed by the code is reached')
  end if


! Store all the previous quantities in a better way than in ref2at (without using too memory).
  write(InVar%stdout,*) ' Build the Shell2at datatype...'
  ABI_MALLOC(Shell2at%neighbours,(natom,Shell2at%nshell))
  ABI_MALLOC(Shell2at%iatref          ,(Shell2at%nshell)); Shell2at%iatref     (:)=zero
  ABI_MALLOC(Shell2at%jatref          ,(Shell2at%nshell)); Shell2at%jatref     (:)=zero
  do ishell=1,Shell2at%nshell
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
  ABI_MALLOC(Shell2at%ncoeff     ,(Shell2at%nshell)); Shell2at%ncoeff(:)=zero
  ABI_MALLOC(Shell2at%ncoeff_prev,(Shell2at%nshell)); Shell2at%ncoeff_prev(:)=zero
  write(InVar%stdout,*) ' Number of shells=',Shell2at%nshell
  write(InVar%stdout,*) '============================================================================'
  if (MPIdata%iam_master) open(unit=16,file=trim(InVar%output_prefix)//'nbcoeff-phi2.dat')
  ncoeff_prev=0
  do ishell=1,Shell2at%nshell
    ncoeff=0
    iatref=Shell2at%iatref(ishell)
    jatref=Shell2at%jatref(ishell)
    write(InVar%stdout,*) 'Shell number:',ishell
    write(InVar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',iatref,' and ',jatref,' the distance is=',distance(iatref,jatref,1)
    call tdep_calc_nbcoeff(distance,iatref,InVar,ishell,jatref,1,1,MPIdata,ncoeff,norder,Shell2at%nshell,order,proj,Sym)
    Shell2at%ncoeff     (ishell)=ncoeff
    Shell2at%ncoeff_prev(ishell)=ncoeff_prev
    ncoeff_prev=ncoeff_prev+ncoeff
    write(InVar%stdout,*)'  Number of independant coefficients in this shell=',ncoeff
    write(InVar%stdout,*)'  Number of interactions in this shell=',Shell2at%neighbours(iatref,ishell)%n_interactions
!FB    write(InVar%stdout,*)'  The ratio is=',dfloat(Shell2at%neighbours(iatref,ishell)%n_interactions)/dfloat(ncoeff)
    write(InVar%stdout,*) '============================================================================'
  end do  
  write(InVar%stdout,*)'  >>>>>> Total number of coefficients at the second order=',ncoeff_prev
  if (MPIdata%iam_master) close(16)
  ntotcoeff=ncoeff_prev
!BeginFB 
!FB  open(unit=91,file='Shell2at.dat')
!FB  write(91,*) Shell2at%nshell
!FB  do ishell=1,Shell2at%nshell
!FB    write(91,*) Shell2at%ncoeff(ishell) 
!FB    write(91,*) Shell2at%ncoeff_prev(ishell) 
!FB    write(91,*) Shell2at%iatref(ishell) 
!FB    write(91,*) Shell2at%jatref(ishell) 
!FB    do iatom=1,InVar%natom
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
 subroutine tdep_init_shell3at(distance,InVar,MPIdata,norder,nshell_max,ntotcoeff,order,proj,Shell3at,Sym)

  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(out) :: Shell3at
  type(Symetries_Variables_type),intent(inout) :: Sym
  type(MPI_enreg_type), intent(in) :: MPIdata
  integer, intent(in) :: norder,order,nshell_max
  integer, intent(out) :: ntotcoeff
  double precision, intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision, intent(out) :: proj(norder,norder,nshell_max)

  integer :: ii,ishell,iatom,jatom,katom,eatom,fatom,gatom,iatref,jatref,katref
  integer :: natom,natom_unitcell,watom,xatom,yatom,ninteractions,ncoeff,ncoeff_prev,nshell_tmp
  integer :: find_equivalent,jj,ninter,iat_ref,jat_ref,kat_ref,ok,tmpinter,isym,itrans
  double precision :: norma,normb,normc
  double precision :: vectj(3),vectk(3),vect1(3),vect2(3)
  integer :: Isym3at(2)
  integer, allocatable :: atref(:,:),interactions(:,:)

  natom         =InVar%natom
  natom_unitcell=InVar%natom_unitcell

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '###### THIRD ORDER : find the number of coefficients ########################'
  write(InVar%stdout,*) '#############################################################################'

! 1/ Identify the shells
  ABI_MALLOC(interactions,(natom_unitcell,nshell_max)) ; interactions(:,:)=0
  ABI_MALLOC(atref,(nshell_max,3)) ; atref(:,:)=0
  ABI_MALLOC(Shell3at%ishell_self,(natom_unitcell)) ; Shell3at%ishell_self(:)=zero
  nshell_tmp=0
  do iatom=1,natom_unitcell
    do jatom=1,natom
      do katom=1,natom
!FB        write(6,*) 'NEW COORD1 : iatom,jatom,katom=',iatom,jatom,katom
!       WARNING: distance(j,k).ne.|djk| due to the inbox procedure when computing distance(j,k). 
!                So, compute |djk| using vec(ij) and vec(ik).      
        norma=dsqrt((distance(iatom,katom,2)-distance(iatom,jatom,2))**2+&
&                   (distance(iatom,katom,3)-distance(iatom,jatom,3))**2+&
&                   (distance(iatom,katom,4)-distance(iatom,jatom,4))**2)
!       Interactions are only computed until Rcut3 in order to have complete shell of neighbours.
!       Otherwise the symetries are broken.
        if ((distance(iatom,jatom,1).gt.(InVar%Rcut3*0.99)).or.&
&           (norma                  .gt.(InVar%Rcut3*0.99)).or.&
&           (distance(iatom,katom,1).gt.(InVar%Rcut3*0.99))) cycle
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
            normb=dsqrt((distance(iat_ref,kat_ref,2)-distance(iat_ref,jat_ref,2))**2+&
&                       (distance(iat_ref,kat_ref,3)-distance(iat_ref,jat_ref,3))**2+&
&                       (distance(iat_ref,kat_ref,4)-distance(iat_ref,jat_ref,4))**2)
            do ii=1,6
              if (ii.eq.1) then ; eatom=iatom ; fatom=jatom ; gatom=katom ; end if
              if (ii.eq.2) then ; eatom=iatom ; fatom=katom ; gatom=jatom ; end if
              if (ii.eq.3) then ; eatom=jatom ; fatom=iatom ; gatom=katom ; end if
              if (ii.eq.4) then ; eatom=jatom ; fatom=katom ; gatom=iatom ; end if
              if (ii.eq.5) then ; eatom=katom ; fatom=iatom ; gatom=jatom ; end if
              if (ii.eq.6) then ; eatom=katom ; fatom=jatom ; gatom=iatom ; end if
              normc=dsqrt((distance(eatom,gatom,2)-distance(eatom,fatom,2))**2+&
&                         (distance(eatom,gatom,3)-distance(eatom,fatom,3))**2+&
&                         (distance(eatom,gatom,4)-distance(eatom,fatom,4))**2)
!FB              if ((abs(distance(iatom,jatom,1)-distance(eatom,fatom,1)).lt.1.d-3).and.&
!FB&                 (abs(norma                  -normb                  ).lt.1.d-3).and.&
!FB&                 (abs(distance(iatom,katom,1)-distance(eatom,gatom,1)).lt.1.d-3)) then
              if ((abs(distance(iat_ref,jat_ref,1)-distance(eatom,fatom,1)).lt.1.d-6).and.&
&                 (abs(normb                      -normc                  ).lt.1.d-6).and.&
&                 (abs(distance(iat_ref,kat_ref,1)-distance(eatom,gatom,1)).lt.1.d-6)) then
                Isym3at(:)=0
                call tdep_SearchS_3at(InVar,iat_ref,jat_ref,kat_ref,eatom,fatom,gatom,Isym3at,Sym,InVar%xred_ideal)
                if (Isym3at(2).eq.1) find_equivalent=1
                if (find_equivalent.eq.1) then 
                  interactions(iatom,ishell)=interactions(iatom,ishell)+1
!FB                  write(6,*) 'The number of interactions in this shell is=',ishell,interactions(iatom,ishell)
                  exit
                end if  
              end if
            end do !ii
            if (find_equivalent.eq.1) exit
          end do !ishell
          if (find_equivalent.eq.0) then
            nshell_tmp=nshell_tmp+1
            if (nshell_tmp.gt.nshell_max) then
              MSG_ERROR('The shell number index is greater than the shell number max defined in the code')
            end if  
            if ((iatom.eq.jatom).and.(jatom.eq.katom)) Shell3at%ishell_self(iatom)=nshell_tmp
            interactions(iatom,nshell_tmp)=interactions(iatom,nshell_tmp)+1
            atref(nshell_tmp,1)=iatom
            atref(nshell_tmp,2)=jatom
            atref(nshell_tmp,3)=katom
!FB            write(6,'(a,1x,4(i5,1x))') 'NEW SHELL1 : nshell_tmp,iatom,jatom,katom=',nshell_tmp,iatom,jatom,katom
          end if  
        end if  
      end do !katom
    end do !jatom
  end do !iatom
  ABI_FREE(atref)

! 2/ Allocate the datatype Shell3at%...
  Shell3at%nshell=nshell_tmp
  ABI_MALLOC(Shell3at%neighbours,(natom,Shell3at%nshell))
  ABI_MALLOC(Shell3at%iatref,(Shell3at%nshell)); Shell3at%iatref(:)=zero
  ABI_MALLOC(Shell3at%jatref,(Shell3at%nshell)); Shell3at%jatref(:)=zero
  ABI_MALLOC(Shell3at%katref,(Shell3at%nshell)); Shell3at%katref(:)=zero
  do ishell=1,Shell3at%nshell
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
  ABI_MALLOC(interactions,(natom,Shell3at%nshell)) ; interactions(:,:)=0
  nshell_tmp=0
  do iatom=1,natom
    do jatom=1,natom
      do katom=1,natom
!FB        write(6,*) 'NEW COORD2 : iatom,jatom,katom=',iatom,jatom,katom
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
!       WARNING: distance(j,k).ne.|djk| due to the inbox procedure when computing distance(j,k). 
!                So, compute |djk| using vec(ij) and vec(ik).      
        norma=dsqrt((distance(iatom,katom,2)-distance(iatom,jatom,2))**2+&
&                   (distance(iatom,katom,3)-distance(iatom,jatom,3))**2+&
&                   (distance(iatom,katom,4)-distance(iatom,jatom,4))**2)
!       Interactions are only computed until Rcut3<acell/2 in order to have complete shell of neighbours.
!       Otherwise the symetries are broken.
        if ((distance(iatom,jatom,1).gt.(InVar%Rcut3*0.99)).or.&
&           (norma                  .gt.(InVar%Rcut3*0.99)).or.&
&           (distance(iatom,katom,1).gt.(InVar%Rcut3*0.99))) cycle
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
          normb=dsqrt((distance(iat_ref,kat_ref,2)-distance(iat_ref,jat_ref,2))**2+&
&                     (distance(iat_ref,kat_ref,3)-distance(iat_ref,jat_ref,3))**2+&
&                     (distance(iat_ref,kat_ref,4)-distance(iat_ref,jat_ref,4))**2)
          do ii=1,6
            if (ii.eq.1) then ; eatom=iatom ; fatom=jatom ; gatom=katom ; end if
            if (ii.eq.2) then ; eatom=iatom ; fatom=katom ; gatom=jatom ; end if
            if (ii.eq.3) then ; eatom=jatom ; fatom=iatom ; gatom=katom ; end if
            if (ii.eq.4) then ; eatom=jatom ; fatom=katom ; gatom=iatom ; end if
            if (ii.eq.5) then ; eatom=katom ; fatom=iatom ; gatom=jatom ; end if
            if (ii.eq.6) then ; eatom=katom ; fatom=jatom ; gatom=iatom ; end if
            normc=dsqrt((distance(eatom,gatom,2)-distance(eatom,fatom,2))**2+&
&                       (distance(eatom,gatom,3)-distance(eatom,fatom,3))**2+&
&                       (distance(eatom,gatom,4)-distance(eatom,fatom,4))**2)
!FB            if ((abs(distance(iatom,jatom,1)-distance(eatom,fatom,1)).lt.1.d-3).and.&
!FB&               (abs(norma                  -normb                  ).lt.1.d-3).and.&
!FB&               (abs(distance(iatom,katom,1)-distance(eatom,gatom,1)).lt.1.d-3)) then
            if ((abs(distance(iat_ref,jat_ref,1)-distance(eatom,fatom,1)).lt.1.d-6).and.&
&               (abs(normb                      -normc                  ).lt.1.d-6).and.&
&               (abs(distance(iat_ref,kat_ref,1)-distance(eatom,gatom,1)).lt.1.d-6)) then
              Isym3at(:)=0
              call tdep_SearchS_3at(InVar,iat_ref,jat_ref,kat_ref,eatom,fatom,gatom,Isym3at,Sym,InVar%xred_ideal)
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
          if (nshell_tmp.gt.Shell3at%nshell) then
            MSG_ERROR('The shell number index is greater than the shell number max computed previously')
          end if  
          Shell3at%iatref(nshell_tmp)=iatom
          Shell3at%jatref(nshell_tmp)=jatom
          Shell3at%katref(nshell_tmp)=katom
          eatom=iatom ; fatom=jatom ; gatom=katom
          Isym3at(:)=1
          ishell=nshell_tmp
!FB          write(6,'(a,1x,4(i5,1x))') 'NEW SHELL2 : nshell_tmp,iatom,jatom,katom=',nshell_tmp,iatom,jatom,katom
        end if  
!       Classify the informations of the triplet in Shell3at
        do ii=1,6
!         The Phi3 has to be symetric (transposition symetries)
          if (ii==1) then ; watom=eatom ; xatom=fatom ; yatom=gatom ; endif !\Psi_ijk
          if (ii==2) then ; watom=eatom ; xatom=gatom ; yatom=fatom ; endif !\Psi_ikj
          if (ii==3) then ; watom=fatom ; xatom=eatom ; yatom=gatom ; endif !\Psi_jik
          if (ii==4) then ; watom=fatom ; xatom=gatom ; yatom=eatom ; endif !\Psi_jki
          if (ii==5) then ; watom=gatom ; xatom=eatom ; yatom=fatom ; endif !\Psi_kij
          if (ii==6) then ; watom=gatom ; xatom=fatom ; yatom=eatom ; endif !\Psi_kji
!         Do not overwrite the Psi_iik, Psi_iji, Psi_ijj or Psi_iii IFCs
!         and avoid double counting of triplet interactions
          if ((eatom.eq.fatom).and.((ii.eq.3).or.(ii.eq.4).or.(ii.eq.6))) cycle
          if ((eatom.eq.gatom).and.((ii.gt.3))) cycle
          if ((fatom.eq.gatom).and.((ii.eq.2).or.(ii.eq.5).or.(ii.eq.6))) cycle
          if ((eatom.eq.fatom).and.(fatom.eq.gatom).and.(ii.gt.1)) cycle
          interactions(watom,ishell)=interactions(watom,ishell)+1
!FB          write(6,*) 'For ishell and eatom=',ishell,watom
!FB          write(6,*) '  --> the number of interactions in the shell is=',interactions(watom,ishell)
          if (interactions(watom,ishell).gt.Shell3at%neighbours(watom,ishell)%n_interactions) then
                write(6,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
                write(6,*) ' >>>>>> Verify that the Rcut used in the input file is lower '
                write(6,*) ' >>>>>> than half of the smallest lattice parameter'
                write(6,*) ' >>>>>> Solution : Reduce the Rcut parameter'
                write(6,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
            MSG_ERROR('The interaction number index is greater than the interaction number max computed previously (3rd order)')
          end if
          Shell3at%neighbours(watom,ishell)%atomj_in_shell(interactions(watom,ishell))=xatom
          Shell3at%neighbours(watom,ishell)%atomk_in_shell(interactions(watom,ishell))=yatom
          Shell3at%neighbours(watom,ishell)%sym_in_shell(interactions(watom,ishell))=Isym3at(1)
          Shell3at%neighbours(watom,ishell)%transpose_in_shell(interactions(watom,ishell))=ii
!DEBUG          write(6,'(a,9(i5,x))') 'ishell,iatref,jatref,katref,iatom,atomj_in_shell,atomk_in_shell,isym,itrans=',&
!DEBUG&         ishell,Shell3at%iatref(ishell),Shell3at%jatref(ishell),Shell3at%katref(ishell),watom,xatom,yatom,Isym3at(1),ii
        end do !ii
      end do !katom
    end do !jatom
  end do !iatom
! Check that each interaction has different symmetry per shell
  do ishell=1,Shell3at%nshell 
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
            MSG_ERROR('Some interactions are equals due to the symmetry')
          end if
        end do
      end do
    end do
  end do  
! Check that each equivalent shell has the same set of interactions
  do ishell=1,Shell3at%nshell 
    do iatom=1,natom
      if (Shell3at%neighbours(mod(iatom-1,natom_unitcell)+1,ishell)%n_interactions.ne.&
&         Shell3at%neighbours(                        iatom,ishell)%n_interactions) then
        MSG_ERROR('The interaction number index is not equal to the interaction number max computed previously (2)')
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
!DEBUG            vectj(ii)=vectj(ii)+Sym%S_ref(ii,jj,isym,1)*distance(iatref,jatref,jj+1)
!DEBUG            vectk(ii)=vectk(ii)+Sym%S_ref(ii,jj,isym,1)*distance(iatref,katref,jj+1)
!DEBUG          end do  
!DEBUG        end do
!DEBUG        if (itrans==1) then ; vect1(:)= vectj(:)          ; vect2(:)= vectk(:)           ; endif !\Psi_ijk
!DEBUG        if (itrans==2) then ; vect1(:)= vectk(:)          ; vect2(:)= vectj(:)           ; endif !\Psi_ikj
!DEBUG        if (itrans==3) then ; vect1(:)=-vectj(:)          ; vect2(:)= vectk(:)-vectj(:)  ; endif !\Psi_jik
!DEBUG        if (itrans==4) then ; vect1(:)= vectk(:)-vectj(:) ; vect2(:)=-vectj(:)           ; endif !\Psi_jki
!DEBUG        if (itrans==5) then ; vect1(:)=-vectk(:)          ; vect2(:)= vectj(:)-vectk(:)  ; endif !\Psi_kij
!DEBUG        if (itrans==6) then ; vect1(:)= vectj(:)-vectk(:) ; vect2(:)=-vectk(:)           ; endif !\Psi_kji
!DEBUG        do ii=1,3
!DEBUG          if ((abs(distance(iatom,jatom,ii+1)-vect1(ii)).gt.tol8).or.&
!DEBUG&             (abs(distance(iatom,katom,ii+1)-vect2(ii)).gt.tol8)) then
!DEBUG            write(std_out,'(a,4(x,i5))') 'For ishell, iatom, jatom, katom =',ishell,iatom,jatom,katom
!DEBUG            write(std_out,'(a,5(x,i5))') '  with isym, itrans, iatref, jatref, katref = ',isym,itrans,iatref,jatref,katref
!DEBUG            MSG_ERROR('We do not recover the triplet with the symmetry found')
!DEBUG          end if
!DEBUG        end do !ii 
!DEBUG      end do !ninter
    end do !natom
  end do !nshell
  ABI_FREE(interactions)

! Find the number of coefficients of the (3x3x3) Phi3 for a given shell
  ABI_MALLOC(Shell3at%ncoeff     ,(Shell3at%nshell)); Shell3at%ncoeff(:)=zero
  ABI_MALLOC(Shell3at%ncoeff_prev,(Shell3at%nshell)); Shell3at%ncoeff_prev(:)=zero
  write(InVar%stdout,*) 'Number of shells=',Shell3at%nshell
  write(InVar%stdout,*) '============================================================================'
  if (MPIdata%iam_master) open(unit=16,file=trim(InVar%output_prefix)//'nbcoeff-phi3.dat')
  ncoeff_prev=0
  do ishell=1,Shell3at%nshell
    ncoeff=0
    iatref=Shell3at%iatref(ishell)
    jatref=Shell3at%jatref(ishell)
    katref=Shell3at%katref(ishell)
    write(InVar%stdout,*) 'Shell number:',ishell
    write(InVar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',iatref,' and ',jatref,' the distance is=',distance(iatref,jatref,1)
    write(InVar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',jatref,' and ',katref,' the distance is=',distance(jatref,katref,1)
    write(InVar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',katref,' and ',iatref,' the distance is=',distance(katref,iatref,1)
    call tdep_calc_nbcoeff(distance,iatref,InVar,ishell,jatref,katref,1,MPIdata,ncoeff,norder,Shell3at%nshell,order,proj,Sym)
    Shell3at%ncoeff     (ishell)=ncoeff
    Shell3at%ncoeff_prev(ishell)=ncoeff_prev
    ncoeff_prev=ncoeff_prev+ncoeff
    write(InVar%stdout,*)'  Number of independant coefficients in this shell=',ncoeff
    write(InVar%stdout,*)'  Number of interactions in this shell=',Shell3at%neighbours(iatref,ishell)%n_interactions
!FB    write(InVar%stdout,*)'  The ratio is=',dfloat(Shell3at%neighbours(iatref,ishell)%n_interactions)/dfloat(ncoeff)
    write(InVar%stdout,*) '============================================================================'
  end do  
  write(InVar%stdout,*)'  >>>>>> Total number of coefficients at the third order=',ncoeff_prev
  if (MPIdata%iam_master) close(16)
  ntotcoeff=ncoeff_prev
!BeginFB 
!FB  open(unit=91,file='Shell3at.dat')
!FB  write(91,*) Shell3at%nshell
!FB  do ishell=1,Shell3at%nshell
!FB    write(91,*) Shell3at%ncoeff(ishell) 
!FB    write(91,*) Shell3at%ncoeff_prev(ishell) 
!FB    write(91,*) Shell3at%iatref(ishell) 
!FB    write(91,*) Shell3at%jatref(ishell) 
!FB    write(91,*) Shell3at%katref(ishell) 
!FB    do iatom=1,InVar%natom
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
 subroutine tdep_init_shell4at(distance,InVar,MPIdata,norder,nshell_max,ntotcoeff,order,proj,Shell4at,Sym)

  implicit none
  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(out) :: Shell4at
  type(Symetries_Variables_type),intent(inout) :: Sym
  type(MPI_enreg_type), intent(in) :: MPIdata
  integer, intent(in) :: norder,order,nshell_max
  integer, intent(out) :: ntotcoeff
  double precision, intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision, intent(out) :: proj(norder,norder,nshell_max)

  integer :: ii,ishell,iatcell,iatom,jatom,katom,latom,eatom,fatom,gatom,hatom,iatref,jatref,katref,latref
  integer :: natom,natom_unitcell,watom,xatom,yatom,zatom,ninteractions,ncoeff,ncoeff_prev,nshell_tmp
  integer :: already_found,find_equivalent,jj,ninter,iat_ref,jat_ref,kat_ref,lat_ref,ok,tmpinter,isym,itrans
  double precision :: norma1,norma2,norma3
  double precision :: normb1,normb2,normb3
  double precision :: normc1,normc2,normc3
  double precision :: vectj(3),vectk(3),vect1(3),vect2(3)
  integer :: Isym4at(2)
  integer, allocatable :: atref(:,:),interactions(:,:)

  natom         =InVar%natom
  natom_unitcell=InVar%natom_unitcell

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '###### FOURTH ORDER : find the number of coefficients ########################'
  write(InVar%stdout,*) '#############################################################################'

! 1/ Identify the shells
  ABI_MALLOC(interactions,(natom_unitcell,nshell_max)) ; interactions(:,:)=0
  ABI_MALLOC(atref,(nshell_max,4)) ; atref(:,:)=0
  ABI_MALLOC(Shell4at%ishell_self,(natom_unitcell)) ; Shell4at%ishell_self(:)=zero
  nshell_tmp=0
  do iatom=1,natom_unitcell
    do jatom=1,natom
!     Interactions are only computed until Rcut4 in order to have complete shell of neighbours.
!     Otherwise the symetries are broken.
      if (distance(iatom,jatom,1).gt.(InVar%Rcut4*0.99)) cycle
      do katom=1,natom
        if (distance(iatom,katom,1).gt.(InVar%Rcut4*0.99)) cycle
!         WARNING: distance(j,k).ne.|djk| due to the inbox procedure when computing distance(j,k). 
!                  So, compute |djk| using vec(ij) and vec(ik).      
          norma1=dsqrt((distance(iatom,katom,2)-distance(iatom,jatom,2))**2+&
&                      (distance(iatom,katom,3)-distance(iatom,jatom,3))**2+&
&                      (distance(iatom,katom,4)-distance(iatom,jatom,4))**2)
        if (norma1                 .gt.(InVar%Rcut4*0.99)) cycle
        do latom=1,natom
!FB          write(6,*) 'NEW COORD1 : iatom,jatom,katom=',iatom,jatom,katom,latom
          if (distance(iatom,latom,1).gt.(InVar%Rcut4*0.99)) cycle
          norma2=dsqrt((distance(iatom,latom,2)-distance(iatom,jatom,2))**2+&
&                      (distance(iatom,latom,3)-distance(iatom,jatom,3))**2+&
&                      (distance(iatom,latom,4)-distance(iatom,jatom,4))**2)
          if (norma2                 .gt.(InVar%Rcut4*0.99)) cycle
          norma3=dsqrt((distance(iatom,latom,2)-distance(iatom,katom,2))**2+&
&                      (distance(iatom,latom,3)-distance(iatom,katom,3))**2+&
&                      (distance(iatom,latom,4)-distance(iatom,katom,4))**2)
          if (norma3                 .gt.(InVar%Rcut4*0.99)) cycle
              
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
            normb1=dsqrt((distance(iat_ref,kat_ref,2)-distance(iat_ref,jat_ref,2))**2+&
&                        (distance(iat_ref,kat_ref,3)-distance(iat_ref,jat_ref,3))**2+&
&                        (distance(iat_ref,kat_ref,4)-distance(iat_ref,jat_ref,4))**2)
            normb2=dsqrt((distance(iat_ref,lat_ref,2)-distance(iat_ref,jat_ref,2))**2+&
&                        (distance(iat_ref,lat_ref,3)-distance(iat_ref,jat_ref,3))**2+&
&                        (distance(iat_ref,lat_ref,4)-distance(iat_ref,jat_ref,4))**2)
            normb3=dsqrt((distance(iat_ref,lat_ref,2)-distance(iat_ref,kat_ref,2))**2+&
&                        (distance(iat_ref,lat_ref,3)-distance(iat_ref,kat_ref,3))**2+&
&                        (distance(iat_ref,lat_ref,4)-distance(iat_ref,kat_ref,4))**2)
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

              normc1=dsqrt((distance(eatom,gatom,2)-distance(eatom,fatom,2))**2+&
&                          (distance(eatom,gatom,3)-distance(eatom,fatom,3))**2+&
&                          (distance(eatom,gatom,4)-distance(eatom,fatom,4))**2)
              normc2=dsqrt((distance(eatom,hatom,2)-distance(eatom,fatom,2))**2+&
&                          (distance(eatom,hatom,3)-distance(eatom,fatom,3))**2+&
&                          (distance(eatom,hatom,4)-distance(eatom,fatom,4))**2)
              normc3=dsqrt((distance(eatom,hatom,2)-distance(eatom,gatom,2))**2+&
&                          (distance(eatom,hatom,3)-distance(eatom,gatom,3))**2+&
&                          (distance(eatom,hatom,4)-distance(eatom,gatom,4))**2)
              if ((abs(distance(iat_ref,jat_ref,1)-distance(eatom,fatom,1)).lt.1.d-6).and.&
&                 (abs(normb1                     -normc1                 ).lt.1.d-6).and.&
&                 (abs(normb2                     -normc2                 ).lt.1.d-6).and.&
&                 (abs(normb3                     -normc3                 ).lt.1.d-6).and.&
&                 (abs(distance(iat_ref,kat_ref,1)-distance(eatom,gatom,1)).lt.1.d-6).and.&
&                 (abs(distance(iat_ref,lat_ref,1)-distance(eatom,hatom,1)).lt.1.d-6)) then
                Isym4at(:)=0
                call tdep_SearchS_4at(InVar,iat_ref,jat_ref,kat_ref,lat_ref,eatom,fatom,gatom,hatom,Isym4at,Sym,InVar%xred_ideal)
                if (Isym4at(2).eq.1) find_equivalent=1
                if (find_equivalent.eq.1) then 
                  interactions(iatom,ishell)=interactions(iatom,ishell)+1
!FB                  write(6,*) 'The number of interactions in this shell is=',ishell,interactions(iatom,ishell)
                  exit
                end if  
              end if
            end do !ii
            if (find_equivalent.eq.1) exit
          end do !ishell
          if (find_equivalent.eq.0) then
            nshell_tmp=nshell_tmp+1
            if (nshell_tmp.gt.nshell_max) then
              MSG_ERROR('The shell number index is greater than the shell number max defined in the code')
            end if  
            if ((iatom.eq.jatom).and.(jatom.eq.katom).and.(katom.eq.latom)) Shell4at%ishell_self(iatom)=nshell_tmp
            interactions(iatom,nshell_tmp)=interactions(iatom,nshell_tmp)+1
            atref(nshell_tmp,1)=iatom
            atref(nshell_tmp,2)=jatom
            atref(nshell_tmp,3)=katom
            atref(nshell_tmp,4)=latom
!FB            write(6,'(a,1x,5(i5,1x))') 'NEW SHELL1 : nshell_tmp,iatom,jatom,katom,latom=',nshell_tmp,iatom,jatom,katom,latom
          end if  
        end do !latom  
      end do !katom
    end do !jatom
  end do !iatom
  ABI_FREE(atref)

! 2/ Allocate the datatype Shell4at%...
  Shell4at%nshell=nshell_tmp
  ABI_MALLOC(Shell4at%neighbours,(natom,Shell4at%nshell))
  ABI_MALLOC(Shell4at%iatref,(Shell4at%nshell)); Shell4at%iatref(:)=zero
  ABI_MALLOC(Shell4at%jatref,(Shell4at%nshell)); Shell4at%jatref(:)=zero
  ABI_MALLOC(Shell4at%katref,(Shell4at%nshell)); Shell4at%katref(:)=zero
  ABI_MALLOC(Shell4at%latref,(Shell4at%nshell)); Shell4at%latref(:)=zero
  do ishell=1,Shell4at%nshell
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
  ABI_MALLOC(interactions,(natom,Shell4at%nshell)) ; interactions(:,:)=0
  nshell_tmp=0
  do iatom=1,natom
    do jatom=1,natom
      if (distance(iatom,jatom,1).gt.(InVar%Rcut4*0.99)) cycle
      do katom=1,natom
        if (distance(iatom,katom,1).gt.(InVar%Rcut4*0.99)) cycle
!         WARNING: distance(j,k).ne.|djk| due to the inbox procedure when computing distance(j,k). 
!                  So, compute |djk| using vec(ij) and vec(ik).      
        norma1=dsqrt((distance(iatom,katom,2)-distance(iatom,jatom,2))**2+&
&                    (distance(iatom,katom,3)-distance(iatom,jatom,3))**2+&
&                    (distance(iatom,katom,4)-distance(iatom,jatom,4))**2)
        if (norma1                 .gt.(InVar%Rcut4*0.99)) cycle
        do latom=1,natom
!FB          write(6,*) 'NEW COORD2 : iatom,jatom,katom=',iatom,jatom,katom,latom
          if (distance(iatom,latom,1).gt.(InVar%Rcut4*0.99)) cycle
          norma2=dsqrt((distance(iatom,latom,2)-distance(iatom,jatom,2))**2+&
&                      (distance(iatom,latom,3)-distance(iatom,jatom,3))**2+&
&                      (distance(iatom,latom,4)-distance(iatom,jatom,4))**2)
          if (norma2                 .gt.(InVar%Rcut4*0.99)) cycle
          norma3=dsqrt((distance(iatom,latom,2)-distance(iatom,katom,2))**2+&
&                      (distance(iatom,latom,3)-distance(iatom,katom,3))**2+&
&                      (distance(iatom,latom,4)-distance(iatom,katom,4))**2)
          if (norma3                 .gt.(InVar%Rcut4*0.99)) cycle
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
            iat_ref=Shell4at%iatref(ishell) ; jat_ref=Shell4at%jatref(ishell) ; kat_ref=Shell4at%katref(ishell) ; lat_ref=Shell4at%latref(ishell)
            normb1=dsqrt((distance(iat_ref,kat_ref,2)-distance(iat_ref,jat_ref,2))**2+&
&                        (distance(iat_ref,kat_ref,3)-distance(iat_ref,jat_ref,3))**2+&
&                        (distance(iat_ref,kat_ref,4)-distance(iat_ref,jat_ref,4))**2)
            normb2=dsqrt((distance(iat_ref,lat_ref,2)-distance(iat_ref,jat_ref,2))**2+&
&                        (distance(iat_ref,lat_ref,3)-distance(iat_ref,jat_ref,3))**2+&
&                        (distance(iat_ref,lat_ref,4)-distance(iat_ref,jat_ref,4))**2)
            normb3=dsqrt((distance(iat_ref,lat_ref,2)-distance(iat_ref,kat_ref,2))**2+&
&                        (distance(iat_ref,lat_ref,3)-distance(iat_ref,kat_ref,3))**2+&
&                        (distance(iat_ref,lat_ref,4)-distance(iat_ref,kat_ref,4))**2)
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

              normc1=dsqrt((distance(eatom,gatom,2)-distance(eatom,fatom,2))**2+&
&                          (distance(eatom,gatom,3)-distance(eatom,fatom,3))**2+&
&                          (distance(eatom,gatom,4)-distance(eatom,fatom,4))**2)
              normc2=dsqrt((distance(eatom,hatom,2)-distance(eatom,fatom,2))**2+&
&                          (distance(eatom,hatom,3)-distance(eatom,fatom,3))**2+&
&                          (distance(eatom,hatom,4)-distance(eatom,fatom,4))**2)
              normc3=dsqrt((distance(eatom,hatom,2)-distance(eatom,gatom,2))**2+&
&                          (distance(eatom,hatom,3)-distance(eatom,gatom,3))**2+&
&                          (distance(eatom,hatom,4)-distance(eatom,gatom,4))**2)
              if ((abs(distance(iat_ref,jat_ref,1)-distance(eatom,fatom,1)).lt.1.d-6).and.&
&                 (abs(normb1                     -normc1                 ).lt.1.d-6).and.&
&                 (abs(normb2                     -normc2                 ).lt.1.d-6).and.&
&                 (abs(normb3                     -normc3                 ).lt.1.d-6).and.&
&                 (abs(distance(iat_ref,kat_ref,1)-distance(eatom,gatom,1)).lt.1.d-6).and.&
&                 (abs(distance(iat_ref,lat_ref,1)-distance(eatom,hatom,1)).lt.1.d-6)) then
                Isym4at(:)=0
                call tdep_SearchS_4at(InVar,iat_ref,jat_ref,kat_ref,lat_ref,eatom,fatom,gatom,hatom,Isym4at,Sym,InVar%xred_ideal)
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
            if (nshell_tmp.gt.Shell4at%nshell) then
              MSG_ERROR('The shell number index is greater than the shell number max computed previously')
            end if  
            Shell4at%iatref(nshell_tmp)=iatom
            Shell4at%jatref(nshell_tmp)=jatom
            Shell4at%katref(nshell_tmp)=katom
            Shell4at%latref(nshell_tmp)=latom
            eatom=iatom ; fatom=jatom ; gatom=katom ; hatom=latom
            Isym4at(:)=1
            ishell=nshell_tmp
!FB            write(6,'(a,1x,5(i5,1x))') 'NEW SHELL2 : nshell_tmp,iatom,jatom,katom=',nshell_tmp,iatom,jatom,katom,latom
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
!FB            write(6,*) 'For ishell and eatom=',ishell,watom
!FB            write(6,*) '  --> the number of interactions in the shell is=',interactions(watom,ishell)
            if (interactions(watom,ishell).gt.Shell4at%neighbours(watom,ishell)%n_interactions) then
                  write(6,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
                  write(6,*) ' >>>>>> Verify that the Rcut used in the input file is lower '
                  write(6,*) ' >>>>>> than half of the smallest lattice parameter'
                  write(6,*) ' >>>>>> Solution : Reduce the Rcut parameter'
                  write(6,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
              MSG_ERROR('The interaction number index is greater than the interaction number max computed previously (4th order)')
            end if
            Shell4at%neighbours(watom,ishell)%atomj_in_shell(interactions(watom,ishell))=xatom
            Shell4at%neighbours(watom,ishell)%atomk_in_shell(interactions(watom,ishell))=yatom
            Shell4at%neighbours(watom,ishell)%atoml_in_shell(interactions(watom,ishell))=zatom
            Shell4at%neighbours(watom,ishell)%sym_in_shell(interactions(watom,ishell))=Isym4at(1)
            Shell4at%neighbours(watom,ishell)%transpose_in_shell(interactions(watom,ishell))=ii
!DEBUG            write(6,'(a,9(i5,x))') 'ishell,iatref,jatref,katref,iatom,atomj_in_shell,atomk_in_shell,isym,itrans=',&
!DEBUG&           ishell,Shell4at%iatref(ishell),Shell4at%jatref(ishell),Shell4at%katref(ishell),watom,xatom,yatom,Isym4at(1),ii
          end do !ii
        end do !latom
      end do !katom
    end do !jatom
  end do !iatom
! Check that each interaction has different symmetry per shell
  do ishell=1,Shell4at%nshell 
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
            MSG_ERROR('Some interactions are equals due to the symmetry')
          end if
        end do
      end do
    end do
  end do  
! Check that each equivalent shell has the same set of interactions
  do ishell=1,Shell4at%nshell 
    do iatom=1,natom
      if (Shell4at%neighbours(mod(iatom-1,natom_unitcell)+1,ishell)%n_interactions.ne.&
&         Shell4at%neighbours(                        iatom,ishell)%n_interactions) then
        MSG_ERROR('The interaction number index is not equal to the interaction number max computed previously (2)')
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
!DEBUG            vectj(ii)=vectj(ii)+Sym%S_ref(ii,jj,isym,1)*distance(iatref,jatref,jj+1)
!DEBUG            vectk(ii)=vectk(ii)+Sym%S_ref(ii,jj,isym,1)*distance(iatref,katref,jj+1)
!DEBUG          end do  
!DEBUG        end do
!DEBUG        if (itrans==1) then ; vect1(:)= vectj(:)          ; vect2(:)= vectk(:)           ; endif !\Psi_ijk
!DEBUG        if (itrans==2) then ; vect1(:)= vectk(:)          ; vect2(:)= vectj(:)           ; endif !\Psi_ikj
!DEBUG        if (itrans==3) then ; vect1(:)=-vectj(:)          ; vect2(:)= vectk(:)-vectj(:)  ; endif !\Psi_jik
!DEBUG        if (itrans==4) then ; vect1(:)= vectk(:)-vectj(:) ; vect2(:)=-vectj(:)           ; endif !\Psi_jki
!DEBUG        if (itrans==5) then ; vect1(:)=-vectk(:)          ; vect2(:)= vectj(:)-vectk(:)  ; endif !\Psi_kij
!DEBUG        if (itrans==6) then ; vect1(:)= vectj(:)-vectk(:) ; vect2(:)=-vectk(:)           ; endif !\Psi_kji
!DEBUG        do ii=1,3
!DEBUG          if ((abs(distance(iatom,jatom,ii+1)-vect1(ii)).gt.tol8).or.&
!DEBUG&             (abs(distance(iatom,katom,ii+1)-vect2(ii)).gt.tol8)) then
!DEBUG            write(std_out,'(a,4(x,i5))') 'For ishell, iatom, jatom, katom =',ishell,iatom,jatom,katom
!DEBUG            write(std_out,'(a,5(x,i5))') '  with isym, itrans, iatref, jatref, katref = ',isym,itrans,iatref,jatref,katref
!DEBUG            MSG_ERROR('We do not recover the quadruplet with the symmetry found')
!DEBUG          end if
!DEBUG        end do !ii 
!DEBUG      end do !ninter
    end do !natom
  end do !nshell
  ABI_FREE(interactions)

! Find the number of coefficients of the (3x3x3x3) Phi4 for a given shell
  ABI_MALLOC(Shell4at%ncoeff     ,(Shell4at%nshell)); Shell4at%ncoeff(:)=zero
  ABI_MALLOC(Shell4at%ncoeff_prev,(Shell4at%nshell)); Shell4at%ncoeff_prev(:)=zero
  write(InVar%stdout,*) 'Number of shells=',Shell4at%nshell
  write(InVar%stdout,*) '============================================================================'
  if (MPIdata%iam_master) open(unit=16,file=trim(InVar%output_prefix)//'nbcoeff-phi4.dat')
  ncoeff_prev=0
  do ishell=1,Shell4at%nshell
    ncoeff=0
    iatref=Shell4at%iatref(ishell)
    jatref=Shell4at%jatref(ishell)
    katref=Shell4at%katref(ishell)
    latref=Shell4at%latref(ishell)
    write(InVar%stdout,*) 'Shell number:',ishell
    write(InVar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',iatref,' and ',jatref,' the distance is=',distance(iatref,jatref,1)
    write(InVar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',jatref,' and ',katref,' the distance is=',distance(jatref,katref,1)
    write(InVar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',katref,' and ',latref,' the distance is=',distance(katref,latref,1)
    write(InVar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',latref,' and ',iatref,' the distance is=',distance(latref,iatref,1)
    call tdep_calc_nbcoeff(distance,iatref,InVar,ishell,jatref,katref,latref,MPIdata,ncoeff,norder,Shell4at%nshell,order,proj,Sym)
    Shell4at%ncoeff     (ishell)=ncoeff
    Shell4at%ncoeff_prev(ishell)=ncoeff_prev
    ncoeff_prev=ncoeff_prev+ncoeff
    write(InVar%stdout,*)'  Number of independant coefficients in this shell=',ncoeff
    write(InVar%stdout,*)'  Number of interactions in this shell=',Shell4at%neighbours(iatref,ishell)%n_interactions
!FB    write(InVar%stdout,*)'  The ratio is=',dfloat(Shell4at%neighbours(iatref,ishell)%n_interactions)/dfloat(ncoeff)
    write(InVar%stdout,*) '============================================================================'
  end do  
  write(InVar%stdout,*)'  >>>>>> Total number of coefficients at the fourth order=',ncoeff_prev
  if (MPIdata%iam_master) close(16)
  ntotcoeff=ncoeff_prev
!BeginFB 
!FB  open(unit=91,file='Shell4at.dat')
!FB  write(91,*) Shell4at%nshell
!FB  do ishell=1,Shell4at%nshell
!FB    write(91,*) Shell4at%ncoeff(ishell) 
!FB    write(91,*) Shell4at%ncoeff_prev(ishell) 
!FB    write(91,*) Shell4at%iatref(ishell) 
!FB    write(91,*) Shell4at%jatref(ishell) 
!FB    write(91,*) Shell4at%katref(ishell) 
!FB    do iatom=1,InVar%natom
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
 subroutine tdep_destroy_shell(natom,order,Shell)

  integer, intent(in) :: natom,order
  type(Shell_Variables_type),intent(inout) :: Shell

  integer :: iatom,ishell

  ABI_FREE(Shell%ncoeff)
  ABI_FREE(Shell%ncoeff_prev)
  ABI_FREE(Shell%iatref)
  ABI_FREE(Shell%jatref)
  ABI_FREE(Shell%ishell_self)
  if (order.gt.2) then
    ABI_FREE(Shell%katref)
  end if
  if (order.gt.3) then
    ABI_FREE(Shell%latref)
  end if
  do iatom=1,natom
    do ishell=1,Shell%nshell
      if (Shell%neighbours(iatom,ishell)%n_interactions.ne.0) then
        ABI_FREE(Shell%neighbours(iatom,ishell)%atomj_in_shell)
        ABI_FREE(Shell%neighbours(iatom,ishell)%sym_in_shell)
        ABI_FREE(Shell%neighbours(iatom,ishell)%transpose_in_shell)
        if (order.gt.2) then
          ABI_FREE(Shell%neighbours(iatom,ishell)%atomk_in_shell)
        end if
        if (order.gt.3) then
          ABI_FREE(Shell%neighbours(iatom,ishell)%atoml_in_shell)
        end if
      end if
    end do
  end do
  ABI_FREE(Shell%neighbours)

 end subroutine tdep_destroy_shell

!====================================================================================================
end module m_tdep_shell
