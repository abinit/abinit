
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_shell

  use defs_basis
  use m_errors
  use m_abicore
  use m_tdep_readwrite,   only : Input_Variables_type
  use m_tdep_latt,        only : Lattice_Variables_type
  use m_tdep_sym,         only : Symetries_Variables_type, tdep_SearchS_2at, tdep_SearchS_3at
  use m_tdep_utils,       only : tdep_calc_nbcoeff
  use m_io_tools

  implicit none

  type List_of_neighbours
    integer :: n_interactions
    integer, allocatable :: atomj_in_shell(:)
    integer, allocatable :: atomk_in_shell(:)
    integer, allocatable :: sym_in_shell(:)
    integer, allocatable :: transpose_in_shell(:)
  end type List_of_neighbours

  type Shell_Variables_type
    integer :: nshell
    integer, allocatable :: ncoeff(:)
    integer, allocatable :: ncoeff_prev(:)
    integer, allocatable :: iatref(:)
    integer, allocatable :: jatref(:)
    integer, allocatable :: katref(:)
    type(List_of_neighbours),allocatable :: neighbours(:,:)
  end type Shell_Variables_type

  public :: tdep_init_shell1at
  public :: tdep_init_shell2at
  public :: tdep_init_shell3at
  public :: tdep_destroy_shell

contains

!====================================================================================================
 subroutine tdep_init_shell1at(distance,InVar,norder,nshell_max,ntotcoeff,order,proj,Shell1at,Sym)

  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(out) :: Shell1at
  type(Symetries_Variables_type),intent(inout) :: Sym
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
    ABI_ERROR('The maximum number of shells allowed by the code is reached')
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

! Find the number of coefficients of the (3x3) Phij for a given shell
  ABI_MALLOC(Shell1at%ncoeff     ,(Shell1at%nshell)); Shell1at%ncoeff(:)=zero
  ABI_MALLOC(Shell1at%ncoeff_prev,(Shell1at%nshell)); Shell1at%ncoeff_prev(:)=zero
  write(InVar%stdout,*) ' Number of shells=',Shell1at%nshell
  write(InVar%stdout,*) '============================================================================'
  open(unit=16,file=trim(InVar%output_prefix)//'nbcoeff-pij.dat')
  ncoeff_prev=0
  do ishell=1,Shell1at%nshell
    ncoeff=0
    iatref=Shell1at%iatref(ishell)
    write(InVar%stdout,*) 'Shell number:',ishell 
    write(InVar%stdout,'(a,i5,a)') '  For atom',iatref,':'
    call tdep_calc_nbcoeff(distance,iatref,InVar,ishell,1,1,ncoeff,norder,Shell1at%nshell,order,proj,Sym)
    if (ncoeff.eq.0) then 
      Shell1at%neighbours(1,ishell)%n_interactions=0
    end if
    Shell1at%ncoeff     (ishell)=ncoeff
    Shell1at%ncoeff_prev(ishell)=ncoeff_prev
    ncoeff_prev=ncoeff_prev+ncoeff
    write(InVar%stdout,*)'  Number of independant coefficients in this shell=',ncoeff
    write(InVar%stdout,*) '============================================================================'
  end do  
  write(InVar%stdout,*)'  >>>>>> Total number of coefficients at the first order=',ncoeff_prev
  close(16)
  ntotcoeff=ncoeff_prev

 end subroutine tdep_init_shell1at

!====================================================================================================
 subroutine tdep_init_shell2at(distance,InVar,norder,nshell_max,ntotcoeff,order,proj,Shell2at,Sym)

  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(out) :: Shell2at
  type(Symetries_Variables_type),intent(inout) :: Sym
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
  call flush_unit(InVar%stdout)
  ABI_MALLOC(ref2at ,(natom,natom,3)) ; ref2at (:,:,:)=zero
  ABI_MALLOC(Isym2at,(natom,natom,2)) ; Isym2at(:,:,:)=zero
  ishell=0
  do iatcell=1,natom_unitcell
    do jatom=1,natom
!     Interactions are only computed until Rcut<acell/2 in order to have complete shell of neighbours.
!     Otherwise the symetries are broken.
      if ((ref2at(iatcell,jatom,1).ne.0).or.(distance(iatcell,jatom,1).gt.(InVar%Rcut*0.99))) cycle
      ishell=ishell+1
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
!             The Phij_NN has to be symetric (transposition symetries)
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
    ABI_ERROR('The maximum number of shells allowed by the code is reached')
  end if


! Store all the previous quantities in a better way than in ref2at (without using too memory).
  write(InVar%stdout,*) ' Build the Shell2at datatype...'
  call flush_unit(InVar%stdout)
  ABI_MALLOC(Shell2at%neighbours,(natom,Shell2at%nshell))
  ABI_MALLOC(Shell2at%iatref          ,(Shell2at%nshell)); Shell2at%iatref(:)=zero
  ABI_MALLOC(Shell2at%jatref          ,(Shell2at%nshell)); Shell2at%jatref(:)=zero
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
  write(InVar%stdout,*) "done"
  call flush_unit(InVar%stdout)

! Find the number of coefficients of the (3x3) Phij for a given shell
  ABI_MALLOC(Shell2at%ncoeff     ,(Shell2at%nshell)); Shell2at%ncoeff(:)=zero
  ABI_MALLOC(Shell2at%ncoeff_prev,(Shell2at%nshell)); Shell2at%ncoeff_prev(:)=zero
  write(InVar%stdout,*) ' Number of shells=',Shell2at%nshell
  write(InVar%stdout,*) '============================================================================'
  open(unit=16,file=trim(InVar%output_prefix)//'nbcoeff-phij.dat')
  ncoeff_prev=0
  do ishell=1,Shell2at%nshell
    ncoeff=0
    iatref=Shell2at%iatref(ishell)
    jatref=Shell2at%jatref(ishell)
    write(InVar%stdout,*) 'Shell number:',ishell
    write(InVar%stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',iatref,' and ',jatref,' the distance is=',distance(iatref,jatref,1)
    call tdep_calc_nbcoeff(distance,iatref,InVar,ishell,jatref,1,ncoeff,norder,Shell2at%nshell,order,proj,Sym)
    Shell2at%ncoeff     (ishell)=ncoeff
    Shell2at%ncoeff_prev(ishell)=ncoeff_prev
    ncoeff_prev=ncoeff_prev+ncoeff
    write(InVar%stdout,*)'  Number of independant coefficients in this shell=',ncoeff
    write(InVar%stdout,*) '============================================================================'
  end do  
  write(InVar%stdout,*)'  >>>>>> Total number of coefficients at the second order=',ncoeff_prev
  close(16)
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
 subroutine tdep_init_shell3at(distance,InVar,norder,nshell_max,ntotcoeff,order,proj,Shell3at,Sym)

  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(out) :: Shell3at
  type(Symetries_Variables_type),intent(inout) :: Sym
  integer, intent(in) :: norder,order,nshell_max
  integer, intent(out) :: ntotcoeff
  double precision, intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision, intent(out) :: proj(norder,norder,nshell_max)

  integer :: ii,ishell,iatom,jatom,katom,eatom,fatom,gatom,iatref,jatref,katref
  integer :: natom,natom_unitcell,watom,xatom,yatom,ninteractions,ncoeff,ncoeff_prev,nshell_tmp
  integer :: find_equivalent,ninter,iat_ref,jat_ref,kat_ref,tmpinter !jj, ok,
  double precision :: norm1,norm2,norm3
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
  nshell_tmp=0
  do iatom=1,natom_unitcell
    do jatom=1,natom
      do katom=1,natom
!FB        write(6,*) 'NEW COORD1 : iatom,jatom,katom=',iatom,jatom,katom
!       WARNING: distance(j,k).ne.|djk| due to the inbox procedure when computing distance(j,k). 
!                So, compute |djk| using vec(ij) and vec(ik).      
        norm1=dsqrt((distance(iatom,katom,2)-distance(iatom,jatom,2))**2+&
&                   (distance(iatom,katom,3)-distance(iatom,jatom,3))**2+&
&                   (distance(iatom,katom,4)-distance(iatom,jatom,4))**2)
!       Interactions are only computed until Rcut3<acell/2 in order to have complete shell of neighbours.
!       Otherwise the symetries are broken.
        if ((distance(iatom,jatom,1).gt.(InVar%Rcut3*0.99)).or.&
&           (norm1                  .gt.(InVar%Rcut3*0.99)).or.&
&           (distance(iatom,katom,1).gt.(InVar%Rcut3*0.99))) cycle
        if (nshell_tmp.eq.0) then
          atref(1,:)=1
          nshell_tmp=nshell_tmp+1
          interactions(iatom,nshell_tmp)=interactions(iatom,nshell_tmp)+1
          cycle
        else   
          find_equivalent=0
          do ishell=1,nshell_tmp
            iat_ref=atref(ishell,1) ; jat_ref=atref(ishell,2) ; kat_ref=atref(ishell,3)
            norm2=dsqrt((distance(iat_ref,kat_ref,2)-distance(iat_ref,jat_ref,2))**2+&
&                       (distance(iat_ref,kat_ref,3)-distance(iat_ref,jat_ref,3))**2+&
&                       (distance(iat_ref,kat_ref,4)-distance(iat_ref,jat_ref,4))**2)
            do ii=1,6
              if (ii.eq.1) then ; eatom=iatom ; fatom=jatom ; gatom=katom ; end if
              if (ii.eq.2) then ; eatom=iatom ; fatom=katom ; gatom=jatom ; end if
              if (ii.eq.3) then ; eatom=jatom ; fatom=iatom ; gatom=katom ; end if
              if (ii.eq.4) then ; eatom=jatom ; fatom=katom ; gatom=iatom ; end if
              if (ii.eq.5) then ; eatom=katom ; fatom=iatom ; gatom=jatom ; end if
              if (ii.eq.6) then ; eatom=katom ; fatom=jatom ; gatom=iatom ; end if
              norm3=dsqrt((distance(eatom,gatom,2)-distance(eatom,fatom,2))**2+&
&                         (distance(eatom,gatom,3)-distance(eatom,fatom,3))**2+&
&                         (distance(eatom,gatom,4)-distance(eatom,fatom,4))**2)
!FB              if ((abs(distance(iatom,jatom,1)-distance(eatom,fatom,1)).lt.1.d-3).and.&
!FB&                 (abs(norm1                  -norm2                  ).lt.1.d-3).and.&
!FB&                 (abs(distance(iatom,katom,1)-distance(eatom,gatom,1)).lt.1.d-3)) then
              if ((abs(distance(iat_ref,jat_ref,1)-distance(eatom,fatom,1)).lt.1.d-6).and.&
&                 (abs(norm2                      -norm3                  ).lt.1.d-6).and.&
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
        norm1=dsqrt((distance(iatom,katom,2)-distance(iatom,jatom,2))**2+&
&                   (distance(iatom,katom,3)-distance(iatom,jatom,3))**2+&
&                   (distance(iatom,katom,4)-distance(iatom,jatom,4))**2)
!       Interactions are only computed until Rcut3<acell/2 in order to have complete shell of neighbours.
!       Otherwise the symetries are broken.
        if ((distance(iatom,jatom,1).gt.(InVar%Rcut3*0.99)).or.&
&           (norm1                  .gt.(InVar%Rcut3*0.99)).or.&
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
          norm2=dsqrt((distance(iat_ref,kat_ref,2)-distance(iat_ref,jat_ref,2))**2+&
&                     (distance(iat_ref,kat_ref,3)-distance(iat_ref,jat_ref,3))**2+&
&                     (distance(iat_ref,kat_ref,4)-distance(iat_ref,jat_ref,4))**2)
          do ii=1,6
            if (ii.eq.1) then ; eatom=iatom ; fatom=jatom ; gatom=katom ; end if
            if (ii.eq.2) then ; eatom=iatom ; fatom=katom ; gatom=jatom ; end if
            if (ii.eq.3) then ; eatom=jatom ; fatom=iatom ; gatom=katom ; end if
            if (ii.eq.4) then ; eatom=jatom ; fatom=katom ; gatom=iatom ; end if
            if (ii.eq.5) then ; eatom=katom ; fatom=iatom ; gatom=jatom ; end if
            if (ii.eq.6) then ; eatom=katom ; fatom=jatom ; gatom=iatom ; end if
            norm3=dsqrt((distance(eatom,gatom,2)-distance(eatom,fatom,2))**2+&
&                       (distance(eatom,gatom,3)-distance(eatom,fatom,3))**2+&
&                       (distance(eatom,gatom,4)-distance(eatom,fatom,4))**2)
!FB            if ((abs(distance(iatom,jatom,1)-distance(eatom,fatom,1)).lt.1.d-3).and.&
!FB&               (abs(norm1                  -norm2                  ).lt.1.d-3).and.&
!FB&               (abs(distance(iatom,katom,1)-distance(eatom,gatom,1)).lt.1.d-3)) then
            if ((abs(distance(iat_ref,jat_ref,1)-distance(eatom,fatom,1)).lt.1.d-6).and.&
&               (abs(norm2                      -norm3                  ).lt.1.d-6).and.&
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
            ABI_ERROR('The shell number index is greater than the shell number max computed previously')
          end if  
          Shell3at%iatref(nshell_tmp)=iatom
          Shell3at%jatref(nshell_tmp)=jatom
          Shell3at%katref(nshell_tmp)=katom
          eatom=iatom ; fatom=jatom ; gatom=katom
          Isym3at(:)=1
          ishell=nshell_tmp
!FB          write(6,'(a,1x,4(i5,1x))') 'NEW SHELL2 : nshell_tmp,iatom,jatom,katom=',nshell_tmp,iatom,jatom,katom
        end if  
!       Classify the information of the triplet in Shell3at
        do ii=1,6
!         The Phij_NN has to be symetric (transposition symetries)
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
            ABI_ERROR('The interaction number index is greater than the interaction number max computed previously')
          end if
          Shell3at%neighbours(watom,ishell)%atomj_in_shell(interactions(watom,ishell))=xatom
          Shell3at%neighbours(watom,ishell)%atomk_in_shell(interactions(watom,ishell))=yatom
          Shell3at%neighbours(watom,ishell)%sym_in_shell(interactions(watom,ishell))=Isym3at(1)
          Shell3at%neighbours(watom,ishell)%transpose_in_shell(interactions(watom,ishell))=ii
!DEBUG          write(6,'(a,9(i5,x))') 'ishell,iatref,jatref,katref,iatom,atomj_in_shell,atomk_in_shell,isym,trans=',&
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
&               Shell3at%neighbours(iatom,ishell)%transpose_in_shell(tmpinter)) cycle
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
  do ishell=1,Shell3at%nshell 
    do iatom=1,natom
      if (Shell3at%neighbours(mod(iatom-1,natom_unitcell)+1,ishell)%n_interactions.ne.&
&         Shell3at%neighbours(                        iatom,ishell)%n_interactions) then
        ABI_ERROR('The interaction number index is not equal to the interaction number max computed previously')
      end if
!FB      do ninter=1,Shell3at%neighbours(iatom,ishell)%n_interactions
!FB        jatom=Shell3at%neighbours(iatom,ishell)%atomj_in_shell(ninter)
!FB	katom=Shell3at%neighbours(iatom,ishell)%atomk_in_shell(ninter)
!FB	isym =Shell3at%neighbours(iatom,ishell)%sym_in_shell(ninter)
!FB	trans=Shell3at%neighbours(iatom,ishell)%transpose_in_shell(ninter)
!FB        do ii=1,3
!FB          do jj=1,3
!FB            vectj(ii)=vectj(ii)+Sym%S_ref(ii,jj,isym,1)*distance(iatom,jatom,jj+1)
!FB            vectk(ii)=vectk(ii)+Sym%S_ref(ii,jj,isym,1)*distance(iatom,katom,jj+1)
!FB          end do  
!FB        end do
!FB        if (trans==1) then ; watom=iatcell ; xatom=jatom   ; yatom=katom   ; endif !\Psi_ijk
!FB        if (trans==2) then ; watom=iatcell ; xatom=katom   ; yatom=jatom   ; endif !\Psi_ikj
!FB        if (trans==3) then ; watom=jatom   ; xatom=iatcell ; yatom=katom   ; endif !\Psi_jik
!FB        if (trans==4) then ; watom=jatom   ; xatom=katom   ; yatom=iatcell ; endif !\Psi_jki
!FB        if (trans==5) then ; watom=katom   ; xatom=iatcell ; yatom=jatom   ; endif !\Psi_kij
!FB        if (trans==6) then ; watom=katom   ; xatom=jatom   ; yatom=iatcell ; endif !\Psi_kji
!FB	if 
!FB      end do
!FB
!FB      if (Shell3at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
!FB      do ninter=1,Shell3at%neighbours(iatom,ishell)%n_interactions
!FB        ok=0
!FB        do tmpinter=1,Shell3at%neighbours(iatom,ishell)%n_interactions
!FB          if ((Shell3at%neighbours(mod(iatom-1,natom_unitcell)+1,ishell)%sym_in_shell(  ninter)).eq.&
!FB&             (Shell3at%neighbours(                      iatom,ishell)%sym_in_shell(tmpinter))) ok=1
!FB	end do  
!FB	if (ok.eq.0) then
!FB          write(std_out,'(a,2(x,i5))') 'For ishell and iatom =',ishell,iatom
!FB          write(std_out,'(a,i5,a)') '  the symmetry ',&
!FB&	    Shell3at%neighbours(mod(iatom-1,natom_unitcell)+1,ishell)%sym_in_shell(  ninter),' is not found'
!FB          do tmpinter=1,Shell3at%neighbours(iatom,ishell)%n_interactions
!FB	    write(6,*) 'INTER=',tmpinter
!FB            write(6,*) Shell3at%neighbours(mod(iatom-1,natom_unitcell)+1,ishell)%sym_in_shell(tmpinter),Shell3at%neighbours(iatom,ishell)%sym_in_shell(tmpinter)
!FB            write(6,*) Shell3at%neighbours(mod(iatom-1,natom_unitcell)+1,ishell)%transpose_in_shell(tmpinter),Shell3at%neighbours(iatom,ishell)%transpose_in_shell(tmpinter)
!FB            write(6,*) Shell3at%neighbours(mod(iatom-1,natom_unitcell)+1,ishell)%atomj_in_shell(tmpinter),Shell3at%neighbours(iatom,ishell)%atomj_in_shell(tmpinter)
!FB            write(6,*) Shell3at%neighbours(mod(iatom-1,natom_unitcell)+1,ishell)%atomk_in_shell(tmpinter),Shell3at%neighbours(iatom,ishell)%atomk_in_shell(tmpinter)
!FB	  end do
!FB          ABI_ERROR('Some symmetries are not found')
!FB	end if
!FB      end do
    end do
  end do
  ABI_FREE(interactions)

! Find the number of coefficients of the (3x3x3) Psij for a given shell
  ABI_MALLOC(Shell3at%ncoeff     ,(Shell3at%nshell)); Shell3at%ncoeff(:)=zero
  ABI_MALLOC(Shell3at%ncoeff_prev,(Shell3at%nshell)); Shell3at%ncoeff_prev(:)=zero
  write(InVar%stdout,*) 'Number of shells=',Shell3at%nshell
  write(InVar%stdout,*) '============================================================================'
  open(unit=16,file=trim(InVar%output_prefix)//'nbcoeff-psij.dat')
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
    call tdep_calc_nbcoeff(distance,iatref,InVar,ishell,jatref,katref,ncoeff,norder,Shell3at%nshell,order,proj,Sym)
    Shell3at%ncoeff     (ishell)=ncoeff
    Shell3at%ncoeff_prev(ishell)=ncoeff_prev
    ncoeff_prev=ncoeff_prev+ncoeff
    write(InVar%stdout,*)'  Number of independant coefficients in this shell=',ncoeff
    write(InVar%stdout,*) '============================================================================'
  end do  
  write(InVar%stdout,*)'  >>>>>> Total number of coefficients at the third order=',ncoeff_prev
  close(16)
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
 subroutine tdep_destroy_shell(natom,order,Shell)

  integer, intent(in) :: natom,order
  type(Shell_Variables_type),intent(inout) :: Shell

  integer :: iatom,ishell

  ABI_FREE(Shell%ncoeff)
  ABI_FREE(Shell%ncoeff_prev)
  ABI_FREE(Shell%iatref)
  ABI_FREE(Shell%jatref)
  if (order.gt.2) then
    ABI_FREE(Shell%katref)
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
      end if
    end do
  end do
  ABI_FREE(Shell%neighbours)

 end subroutine tdep_destroy_shell

!====================================================================================================
end module m_tdep_shell
