
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_utils
  
  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_wffile
  use m_tdep_latt,        only : Lattice_Variables_type, tdep_make_inbox
  use m_tdep_readwrite,   only : Input_Variables_type, MPI_enreg_type
  use m_tdep_sym,         only : Symetries_Variables_type, tdep_SearchS_1at
  use m_io_tools

  implicit none

  type Coeff_Moore_type
    integer :: ntotcoeff
    integer :: ntotconst
    integer :: ncoeff1st
    integer :: ncoeff2nd
    integer :: ncoeff3rd
    integer :: ncoeff4th
    integer :: nconst_1st
    integer :: nconst_2nd
    integer :: nconst_3rd
    integer :: nconst_4th
    integer :: nconst_rot2nd
    integer :: nconst_huang
    integer :: nconst_dynmat
    integer :: nconst_rot3rd
    integer :: nconst_asr3rd
    integer :: nconst_rot4th
    integer :: nconst_asr4th
    double precision, allocatable :: fcoeff(:,:)
    double precision, allocatable :: const(:,:)
  end type Coeff_Moore_type

  type S_product
    double precision, allocatable :: SS  (:,:,:)
    double precision, allocatable :: SSS (:,:,:,:)
    double precision, allocatable :: SSSS(:,:,:,:,:)
  end type S_product

  type Asr_Rot
    double precision, allocatable :: ABG (:,:,:)
    double precision, allocatable :: ABGD(:,:,:,:)
    double precision, allocatable :: ABGDE(:,:,:,:,:)
  end type Asr_Rot

  type,public :: Constraints_Variables_type
    type(S_product),allocatable :: Sprod(:,:)
    type(Asr_Rot),allocatable :: AsrRot3(:,:,:)
    type(Asr_Rot),allocatable :: AsrRot4(:,:,:,:)
  end type Constraints_Variables_type


 public :: tdep_calc_MoorePenrose
 public :: tdep_MatchIdeal2Average
 public :: tdep_calc_model
 public :: tdep_calc_nbcoeff

contains

!=====================================================================================================
 subroutine tdep_calc_MoorePenrose(CoeffMoore,Forces,order,InVar,IFC_coeff,MPIdata)

  implicit none 

  type(Input_Variables_type),intent(in) :: InVar
  type(Coeff_Moore_type), intent(in) :: CoeffMoore
  double precision, intent(in)  :: Forces(3*InVar%natom*InVar%my_nstep)
  double precision, intent(out)  :: IFC_coeff(CoeffMoore%ntotcoeff,1)
  type(MPI_enreg_type), intent(in) :: MPIdata
  integer, intent(in) :: order
  
  integer :: ii,LWORK,INFO,ntotcoeff,ntotconst
  integer :: natnstep,nconcoef,ierr,ncoeff_prev,nconst_prev,iconst,icoeff
  integer, allocatable :: IWORK(:),IPIV(:)
  double precision, allocatable :: WORK(:)
  double precision, allocatable :: ffcoeff_tmp(:,:),fforces_tmp(:),b_const(:)
  double precision, allocatable :: A_tot(:,:),A_inv(:,:),b_tot(:),x_tot(:)

  write(InVar%stdout,*) '################### And compute the pseudo-inverse ##########################'
  write(InVar%stdout,*) '#############################################################################'
  natnstep=3*InVar%natom*InVar%my_nstep

  if (order.eq.0) then
!   Simultaneously (InVar%together=0) 
    ncoeff_prev=0
    nconst_prev=0
    ntotcoeff=CoeffMoore%ntotcoeff
    ntotconst=CoeffMoore%ntotconst
  else if (order.eq.1) then  
!   Successively (InVar%together=1 and InVar%Order=2)
    ncoeff_prev=0
    nconst_prev=0
    ntotcoeff=CoeffMoore%ncoeff1st +CoeffMoore%ncoeff2nd
    ntotconst=CoeffMoore%nconst_1st+CoeffMoore%nconst_2nd
  else if (order.eq.2) then  
!   Successively (InVar%together=1 and InVar%Order=3)
    ncoeff_prev=CoeffMoore%ncoeff1st +CoeffMoore%ncoeff2nd
    nconst_prev=CoeffMoore%nconst_1st+CoeffMoore%nconst_2nd
    ntotcoeff=CoeffMoore%ncoeff3rd
    ntotconst=CoeffMoore%nconst_3rd
  else if (order.eq.3) then  
!   Successively (InVar%together=1 and InVar%Order=4)
    ncoeff_prev=CoeffMoore%ncoeff1st +CoeffMoore%ncoeff2nd +CoeffMoore%ncoeff3rd
    nconst_prev=CoeffMoore%nconst_1st+CoeffMoore%nconst_2nd+CoeffMoore%nconst_3rd
    ntotcoeff=CoeffMoore%ncoeff4th
    ntotconst=CoeffMoore%nconst_4th
  end if
  nconcoef=ntotcoeff+ntotconst
!  write(InVar%stdout,*) 'ncoeff_prev=',ncoeff_prev
!  write(InVar%stdout,*) 'nconst_prev=',nconst_prev
!  write(InVar%stdout,*) 'ntotconst=',ntotconst
!  write(InVar%stdout,*) 'ntotcoeff=',ntotcoeff
!  write(InVar%stdout,*) 'nconcoef=',nconcoef
  if ((ntotconst.gt.0).and.(order.ge.2)) then 
    ABI_MALLOC(b_const,(ntotconst)) ; b_const(:)=0.d0
    do iconst=1,ntotconst
      do icoeff=1,ncoeff_prev
        b_const(iconst)=b_const(iconst)+&
&         CoeffMoore%const(nconst_prev+iconst,icoeff)*IFC_coeff(icoeff,1)
      end do
    end do
  end if  

  ABI_MALLOC(ffcoeff_tmp,(ntotcoeff,ntotcoeff)) ; ffcoeff_tmp(:,:)=0.d0
  ABI_MALLOC(fforces_tmp,(ntotcoeff)) ; fforces_tmp(:)=0.d0
  ABI_MALLOC(A_tot,(nconcoef,nconcoef)) ; A_tot(:,:)=0.d0
  ABI_MALLOC(A_inv,(nconcoef,nconcoef)) ; A_inv(:,:)=0.d0
  ABI_MALLOC(b_tot,(nconcoef)) ; b_tot(:)=0.d0
  ABI_MALLOC(x_tot,(nconcoef)) ; x_tot(:)=0.d0
  call DGEMM('T','N',ntotcoeff,ntotcoeff,natnstep,2.d0,&
&            CoeffMoore%fcoeff(:,ncoeff_prev+1:ncoeff_prev+ntotcoeff),natnstep,&
&            CoeffMoore%fcoeff(:,ncoeff_prev+1:ncoeff_prev+ntotcoeff),natnstep,&
&            0.d0,ffcoeff_tmp,ntotcoeff)
! NOTE, we have to solve F_ij = -\sum_j \Phi_ij u_j, so we add a minus sign
  call DGEMV('T',natnstep,ntotcoeff,-2.d0,&
&            CoeffMoore%fcoeff(:,ncoeff_prev+1:ncoeff_prev+ntotcoeff),natnstep,&
&            Forces,1,0.d0,fforces_tmp,1)
  call xmpi_sum(ffcoeff_tmp,MPIdata%comm_step,ierr)
  call xmpi_sum(fforces_tmp,MPIdata%comm_step,ierr)

  A_tot(1:ntotcoeff,1:ntotcoeff)=ffcoeff_tmp(1:ntotcoeff,1:ntotcoeff)
  ABI_FREE(ffcoeff_tmp)
  if (ntotconst.gt.0) then
    A_tot(ntotcoeff+1:nconcoef,1:ntotcoeff)=CoeffMoore%const(nconst_prev+1:nconst_prev+ntotconst,ncoeff_prev+1:ncoeff_prev+ntotcoeff)
!FB    do iconst=1,ntotconst
!FB      write(InVar%stdout,*) 'A_tot=',A_tot(ntotcoeff+iconst,1:ntotcoeff)
!FB      write(InVar%stdout,*) 'const=',CoeffMoore%const(nconst_prev+iconst,:)
!FB    end do  
    A_tot(1:ntotcoeff,ntotcoeff+1:nconcoef)=transpose(CoeffMoore%const(nconst_prev+1:nconst_prev+ntotconst,ncoeff_prev+1:ncoeff_prev+ntotcoeff))
!FB    ABI_FREE(CoeffMoore%const)
  end if  
  b_tot(1:ntotcoeff)=fforces_tmp(:)
  if ((ntotconst.gt.0).and.(order.ge.2)) then
    b_tot(ntotcoeff+1:nconcoef)=-b_const(1:ntotconst)
    ABI_FREE(b_const)
  end if  
  ABI_FREE(fforces_tmp)

  ABI_MALLOC(IPIV,(nconcoef)); IPIV(:)=0
  ABI_MALLOC(WORK,(nconcoef)); WORK(:)=0.d0
  A_inv(:,:)=A_tot(:,:)
  call DGETRF(nconcoef,nconcoef,A_inv,nconcoef,IPIV,INFO)
  call DGETRI(nconcoef,A_inv,nconcoef,IPIV,WORK,nconcoef,INFO)
  if (INFO.ne.0) write(InVar%stdout,*) 'INFO (dgetri)=',INFO
  ABI_FREE(IPIV)
  ABI_FREE(WORK)

  call DGEMV('N',nconcoef,nconcoef,1.d0,A_inv,nconcoef,b_tot,1,0.d0,x_tot,1)
  write(InVar%stdout,*) ' The solutions are:'
  do ii=1,nconcoef
    write(InVar%stdout,'(1x,i4,1x,f15.10)') ii,x_tot(ii)
  end do
  write(InVar%stdout,'(a,1x,f15.10)')'  condition number=',maxval(x_tot(:))/minval(x_tot(:))

  IFC_coeff(ncoeff_prev+1:ncoeff_prev+ntotcoeff,1)=x_tot(1:ntotcoeff)
  ABI_FREE(A_tot)
  ABI_FREE(A_inv)
  ABI_FREE(b_tot)
  ABI_FREE(x_tot)

 end subroutine tdep_calc_MoorePenrose

!====================================================================================================
 subroutine tdep_MatchIdeal2Average(distance,Forces_MD,InVar,Lattice,MPIdata,&
&                              Rlatt_cart,Rlatt4dos,Sym,ucart)

  implicit none 

  type(Input_Variables_type),intent(inout) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Symetries_Variables_type),intent(inout) :: Sym
  type(MPI_enreg_type),intent(in) :: MPIdata
  double precision, intent(out)  :: distance(InVar%natom,InVar%natom,4)
  double precision, intent(out)  :: Forces_MD(3*InVar%natom*InVar%my_nstep)
  double precision, intent(out)  :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)
  double precision, intent(out)  :: Rlatt4dos (3,InVar%natom_unitcell,InVar%natom)
  double precision, intent(out)  :: ucart(3,InVar%natom,InVar%my_nstep)
  
  integer :: ii,jj,kk,max_ijk,iatcell,jatcell,iatom,jatom,eatom,fatom,istep
  integer :: foo,foo2,atom_ref,ierr
  double precision :: tmp(3),tmp1(3),tmp2(3),Rlatt(3),xred_tmp(3),rprimd_MD_tmp(3,3)
  double precision, allocatable :: dist_unitcell(:,:,:),xcart_average(:,:)
  double precision, allocatable :: fcart_tmp(:,:,:),ucart_tmp(:,:,:)
  double precision, allocatable  :: xred_average(:,:)
  double precision, allocatable  :: xred_center(:,:)
  double precision, allocatable  :: Rlatt_red (:,:,:)
  double precision, allocatable  :: xred_ideal(:,:)
! double precision, allocatable  :: distance_average(:,:,:)
  integer, allocatable  :: FromIdeal2Average(:)
  double precision, allocatable  :: xcart(:,:,:)
  double precision, allocatable  :: xcart_ideal(:,:)
  logical :: ok,ok1

  ierr = 0;

  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '###### Find the matching between ideal and average positions  ###############'
  write(InVar%stdout,*) '#############################################################################'
!==========================================================================================
!======== 1/ Determine ideal positions and distances ======================================
!==========================================================================================
  write(InVar%stdout,*)' Determine ideal positions and distances...'
!Check that atoms (defined in the input.in file) are set correctly in the unitcell
  do ii=1,3
    do iatcell=1,InVar%natom_unitcell
      if ((InVar%xred_unitcell(ii,iatcell).le.(-0.5)).or.(InVar%xred_unitcell(ii,iatcell).gt.(0.5))) then
        do while (InVar%xred_unitcell(ii,iatcell).le.(-0.5))
          InVar%xred_unitcell(ii,iatcell)=InVar%xred_unitcell(ii,iatcell)+1.d0
        end do  
        do while (InVar%xred_unitcell(ii,iatcell).gt.(0.5))
          InVar%xred_unitcell(ii,iatcell)=InVar%xred_unitcell(ii,iatcell)-1.d0
        end do  
!FB        write(InVar%stdout,*) 'xred_unitcell='
!FB        write(InVar%stdout,*)  InVar%xred_unitcell(:,1:InVar%natom_unitcell)
!FB        write(InVar%stdout,*) 'Please put the atoms in the ]-0.5;0.5] range'
!FB        stop -1
      endif
    end do
  end do

! Define the bigbox with ideal positions
  ABI_MALLOC(Rlatt_red ,(3,InVar%natom_unitcell,InVar%natom)); Rlatt_red (:,:,:)=0.d0
  ABI_MALLOC(xred_ideal,(3,InVar%natom))                     ; xred_ideal(:,:)=0.d0
  max_ijk=20
  iatom=1
  do ii=-max_ijk,max_ijk
    do jj=-max_ijk,max_ijk
      do kk=-max_ijk,max_ijk
        do iatcell=1,InVar%natom_unitcell
          if (iatcell==1) ok=.false.
          Rlatt(1)=real(ii-1)
          Rlatt(2)=real(jj-1)
          Rlatt(3)=real(kk-1)
!         Then compute the reduced positions
          tmp(:)=Rlatt(:)+InVar%xred_unitcell(:,iatcell)
          call DGEMV('T',3,3,1.d0,Lattice%multiplicitym1(:,:),3,tmp(:),1,0.d0,xred_tmp(:),1)

!         If the first atom of the pattern is in the [0;1[ range then keep all the
!         atoms of the pattern (even if the others are outside the box). Elsewhere, 
!         none are taken.
          if (iatcell==1) then
            if (minval(xred_tmp(:)).lt.0.d0.or.maxval(xred_tmp(:)).ge.(1.d0-1.d-12)) then
              cycle
            else
              ok=.true.
            end if  
          else 
            if (.not.ok) cycle
          end if
          if (iatom.gt.(InVar%natom+1)) then
            MSG_ERROR('The number of atoms found in the bigbox exceeds natom' )
          end if  
          xred_ideal(:,iatom)=xred_tmp(:)
          call DGEMV('T',3,3,1.d0,Lattice%multiplicitym1(:,:),3,Rlatt(:),1,0.d0,Rlatt_red(:,1,iatom),1)
          iatom=iatom+1
        end do   
      end do
    end do
  end do

  if (iatom.lt.InVar%natom) then
    MSG_ERROR('The number of atoms found in the bigbox is lower than natom')
  end if  

! Compute the distances between ideal positions in the SUPERcell
  do eatom=1,InVar%natom
    do fatom=1,InVar%natom
      tmp(:)=xred_ideal(:,fatom)-xred_ideal(:,eatom)
      call tdep_make_inbox(tmp,1,1d-3)
      call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,tmp(:),1,0.d0,distance(eatom,fatom,2:4),1)
      do ii=1,3
!       Remove the rounding errors before writing (for non regression testing purposes)
        if (abs(distance(eatom,fatom,ii+1)).lt.tol8) distance(eatom,fatom,ii+1)=zero
        distance(eatom,fatom,1)=distance(eatom,fatom,1)+(distance(eatom,fatom,ii+1))**2
      end do
      distance(eatom,fatom,1)=distance(eatom,fatom,1)**0.5
    end do  
  end do  

! Compute the distances between ideal positions in the UNITcell
  ABI_MALLOC(dist_unitcell,(InVar%natom_unitcell,InVar%natom_unitcell,3)); dist_unitcell(:,:,:)=zero
  do iatcell=1,InVar%natom_unitcell
    do jatcell=1,InVar%natom_unitcell
      tmp(:)=xred_ideal(:,jatcell)-xred_ideal(:,iatcell)
      call tdep_make_inbox(tmp,1,tol8)
      dist_unitcell(iatcell,jatcell,:)=tmp(:)
    end do
  end do  

!==========================================================================================
!======== 2/ Find the matching between the ideal and average ==============================
!========   (from the MD simulations) positions. ==========================================
!======== NOTE: - xred_center is used to find the matching with the ideal positions ======= 
!========       - xred_average is used to compute the displacements (from MD trajectories) 
!==========================================================================================
  write(InVar%stdout,*)' Compute average positions...'
  ABI_MALLOC(xred_average,(3,InVar%natom))             ; xred_average(:,:)=0.d0
  ABI_MALLOC(xred_center,(3,InVar%natom))              ; xred_center(:,:)=0.d0
! Average positions from MD (on nstep steps)
  do istep=1,InVar%my_nstep
    do iatom=1,InVar%natom
      xred_average(:,iatom)=xred_average(:,iatom)+InVar%xred(:,iatom,istep)
    end do
  end do
  call xmpi_sum(xred_average,MPIdata%comm_step,ierr)
  xred_average(:,:)=xred_average(:,:)/real(InVar%nstep_tot)

  write(InVar%stdout,*)' Search the unitcell basis of atoms in the MD trajectory...'
! Search the basis of atoms in the supercell
  ok=.true.
  xred_center(:,:)=xred_average(:,:)
  do iatom=1,InVar%natom
    foo2=0
    do jatom=1,InVar%natom
      tmp(:)=xred_center(:,jatom)-xred_center(:,iatom)
      call tdep_make_inbox(tmp,1,1d-3)
      iatcell=1
      do jatcell=1,InVar%natom_unitcell
        foo=0
        do ii=1,3
          if ((abs(tmp(ii)-dist_unitcell(iatcell,jatcell,ii)).le.InVar%tolmotif).and.&
&               InVar%typat(iatom).eq.InVar%typat_unitcell(iatcell).and.&
&               InVar%typat(jatom).eq.InVar%typat_unitcell(jatcell)) then
            foo=foo+1
          end if
        end do
        if (foo==3) then
          foo2=foo2+1
          exit
        end if
      end do
    end do
    if (foo2.eq.InVar%natom_unitcell) then
      atom_ref=iatom
!FB      write(6,*) 'natom_unitcell (ok)=',InVar%natom_unitcell
!FB      write(6,*) 'foo2 (ok)=',foo2
!FB      write(6,*) 'ATOM REF (ok)=',iatom
!FB      write(InVar%stdout,*) 'ATOM REF=',atom_ref
      ok=.false.
      exit
    else if (foo2.gt.InVar%natom_unitcell) then
!FB      write(6,*) 'natom_unitcell (bug)=',InVar%natom_unitcell
!FB      write(6,*) 'foo2 (bug)=',foo2
!FB      write(6,*) 'ATOM REF (bug)=',iatom
      MSG_BUG(' Something wrong: WTF')
    endif
  end do  
  if (ok) then
    if (MPIdata%iam_master) then 
      open(unit=31,file=trim(InVar%output_prefix)//'xred_average.xyz')
      do iatom=1,InVar%natom
        write(31,'(a,1x,3(f10.6,1x))') 'C',xred_center(:,iatom)
        write(31,'(a,1x,3(f10.6,1x))') 'I',xred_ideal (:,iatom)
      end do
      close(31)
    end if
    MSG_ERROR_NOSTOP('The basis of atoms written in input.in file does not appear in the MD trajectory',ierr)
    MSG_ERROR('Perhaps, you can adjust the tolerance (tolmotif)')
  end if  
  ABI_FREE(dist_unitcell)

! Modification of xred and Rlatt tabs
! For averaged quantities --> kk=1: xred_center, xred_average et xred 
! For ideal quantities    --> kk=2: Rlatt_red et xred_ideal 
  write(InVar%stdout,*)' Compare ideal and average positions using PBC...'
  do kk=1,2
!   1/ The "atom_ref" (kk=1) or iatom=1 (kk=2) atom is put in (0.0;0.0;0.0)
    if (kk==1) then
      tmp(:)=xred_center(:,atom_ref)
    else if (kk==2) then
      tmp1(:)=xred_ideal(:,1)
      tmp2(:)=Rlatt_red(:,1,1)
    end if  
    do jatom=1,InVar%natom
      if (kk==1) xred_center(:,jatom)=xred_center(:,jatom)-tmp(:)
      if (kk==2) then
        xred_ideal(:,jatom)=  xred_ideal(:,jatom)  -tmp1(:)
        Rlatt_red (:,1,jatom)=Rlatt_red (:,1,jatom)-tmp2(:)
      end if        
    end do  
!   2/ All the atoms are put in the range [-0.5;0.5[ (use of PBC)
    do jatom=1,InVar%natom
      if (kk==1) then
        tmp(:)=xred_center(:,jatom)
        call tdep_make_inbox(tmp,1,InVar%tolinbox,xred_center (:,jatom))
        call tdep_make_inbox(tmp,1,InVar%tolinbox,xred_average(:,jatom))
        do istep=1,InVar%my_nstep
          call tdep_make_inbox(tmp,1,InVar%tolinbox,InVar%xred(:,jatom,istep))
        end do  
      else if (kk==2) then
        tmp(:)=xred_ideal(:,jatom)
        call tdep_make_inbox(tmp,1,tol8,xred_ideal(:,jatom))
        call tdep_make_inbox(tmp,1,tol8,Rlatt_red(:,1,jatom))
!FB        call tdep_make_inbox(Rlatt_red(:,1,jatom),1,tol8)
      end if  
    end do  
  end do  

! When the multiplicity equals 1 along one direction, there is some trouble
! To clean!!!!!!!
  do ii=1,3
    if ((InVar%multiplicity(ii,ii).eq.1).and.(InVar%multiplicity(ii,mod(ii  ,3)+1).eq.0)&
&                                 .and.(InVar%multiplicity(ii,mod(ii+1,3)+1).eq.0)) then
      Rlatt_red(ii,1,:)=0.d0
      write(InVar%stdout,*) 'WARNING: multiplicity=1 for ii=',ii
    end if
  end do

! Define Rlatt for all the atoms in the basis (Rlatt_red varies as a function of iatcell)
  if (InVar%natom_unitcell.gt.1) then
    do iatcell=2,InVar%natom_unitcell
      Rlatt_red(:,iatcell,:)=Rlatt_red(:,1,:)
    end do
  end if  
  do iatom=1,InVar%natom
    do iatcell=1,InVar%natom_unitcell
      tmp(:)=xred_ideal(:,iatom)-xred_ideal(:,iatcell)
      call tdep_make_inbox(tmp,1,tol8,Rlatt_red(:,iatcell,iatom))
    end do
  end do
  if (InVar%debug) then
    do iatcell=1,InVar%natom_unitcell
      write(InVar%stdout,*) 'For iatcell=',iatcell
      do jatom=1,InVar%natom
        write(InVar%stdout,'(a,i4,a,3(f16.10,1x))') 'For jatom=',jatom,', Rlatt=',Rlatt_red(1:3,iatcell,jatom)
      end do  
    end do  
  end if

! Matching between Ideal and Average positions: xred_ideal and xred_center
! Then, write them in the xred_average.xyz file.
  write(InVar%stdout,*)' Write the xred_average.xyz file with ideal and average positions...'
  ABI_MALLOC(FromIdeal2Average,(InVar%natom))             ; FromIdeal2Average(:)=0
  do iatom=1,InVar%natom
    do jatom=1,InVar%natom
      ok =.true.
      ok1=.true.
      foo=0
      do ii=1,3
        if ((abs(xred_center(ii,iatom)-xred_ideal(ii,jatom)     ).le.InVar%tolmatch.and.&
&         (InVar%typat(iatom).eq.InVar%typat_unitcell(mod(jatom-1,InVar%natom_unitcell)+1)))) then
          foo=foo+1
        else if ((abs(xred_center(ii,iatom)-xred_ideal(ii,jatom)-1.d0).le.InVar%tolmatch.and.&
&         (InVar%typat(iatom).eq.InVar%typat_unitcell(mod(jatom-1,InVar%natom_unitcell)+1))).or.&
&           (abs(xred_center(ii,iatom)-xred_ideal(ii,jatom)+1.d0).le.InVar%tolmatch.and.&
&           (InVar%typat(iatom).eq.InVar%typat_unitcell(mod(jatom-1,InVar%natom_unitcell)+1)))) then
          foo=foo+1
          ok1=.false.
        endif
      end do
      if (foo==3.and.ok1) then
        FromIdeal2Average(jatom)=iatom
        ok=.false.
        exit
      else if (foo==3.and..not.ok1) then
!FB        write(InVar%stdout,*) '  THE CODE STOPS'
!FB        write(InVar%stdout,*) '  Some positions are outside the [-0.5;0.5[ range:'
!FB        write(InVar%stdout,*) '  xred_center(:,',iatom,')='
!FB        write(InVar%stdout,*) xred_center(:,iatom)
!FB        write(InVar%stdout,*) '  xred_ideal(:,',jatom,')='
!FB        write(InVar%stdout,*) xred_ideal(:,jatom)
!FB        write(InVar%stdout,*) 'Perhaps, you can adjust the tolerance (tolinbox)'
!FB        stop -1
        do ii=1,3
          if (abs(xred_center(ii,iatom)-xred_ideal(ii,jatom)-1.d0).le.InVar%tolmatch.and.&
&           (InVar%typat(iatom).eq.InVar%typat_unitcell(mod(jatom-1,InVar%natom_unitcell)+1))) then
!JB            write(*,*) iatom,jatom
!JB            write(*,*) xred_center(:,iatom)
!JB            write(*,*) xred_ideal(:,jatom)
            xred_center(ii,iatom)=xred_center(ii,iatom)-1d0
            do istep=1,InVar%my_nstep
              InVar%xred(:,iatom,istep)=InVar%xred(:,iatom,istep)-1d0
            end do
!JB            write(*,*) xred_center(:,iatom)
!JB            write(*,*) xred_ideal(:,jatom)
            FromIdeal2Average(jatom)=iatom
          else if (abs(xred_center(ii,iatom)-xred_ideal(ii,jatom)+1.d0).le.InVar%tolmatch.and.&
&           (InVar%typat(iatom).eq.InVar%typat_unitcell(mod(jatom-1,InVar%natom_unitcell)+1))) then
!jB            write(*,*) iatom,jatom
!jB            write(*,*) xred_center(:,iatom)
!jB            write(*,*) xred_ideal(:,jatom)
            xred_center(ii,iatom)=xred_center(ii,iatom)+1d0
            do istep=1,InVar%my_nstep
              InVar%xred(:,iatom,istep)=InVar%xred(:,iatom,istep)+1d0
            end do
!JB            write(*,*) xred_center(:,iatom)
!JB            write(*,*) xred_ideal(:,jatom)
            FromIdeal2Average(jatom)=iatom
          end if
        end do
        ok=.false.
        exit
      end if  
    end do  
    if (ok) then
      write(InVar%stdout,*) 'Problem to find the average position for iatom=',iatom
      write(InVar%stdout,*) '  Reasons:'
      write(InVar%stdout,*) '    1/ One atom jump to another equilibrium position'
      write(InVar%stdout,*) '    2/ The system is no more solid'
      write(InVar%stdout,*) '    3/ Perhaps, you can adjust the tolerance (tolmatch)'
      write(InVar%stdout,*) '  xred_center=',(xred_center(ii,iatom),ii=1,3)
      do eatom=1,InVar%natom
        write(InVar%stdout,'(a,1x,3(f10.6,1x))') 'I',xred_ideal (:,eatom)
        write(InVar%stdout,'(a,1x,3(f10.6,1x))') 'C',xred_center(:,eatom)
      end do
      MSG_ERROR('Problem to find the average position')
    end if  
  end do

! WARNING: VERY IMPORTANT: The positions are displayed/sorted (and used in the following) 
! according to ideal positions xred_ideal. 
  if (MPIdata%iam_master) then
    open(unit=31,file=trim(InVar%output_prefix)//'xred_average.xyz')
    write(31,'(i4)') InVar%natom*2
    write(31,'(i4)') 1
!   --> In reduced coordinates
    do iatom=1,InVar%natom
      write(31,'(a,1x,3(f10.6,1x))') 'Ired',xred_ideal (:,iatom)
      write(31,'(a,1x,3(f10.6,1x))') 'Cred',xred_center(:,FromIdeal2Average(iatom))
    end do  
!   --> In cartesian coordinates  
    do iatom=1,InVar%natom
      tmp(:)=zero
      call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,xred_ideal (:,iatom),1,0.d0,tmp(:),1)
      write(31,'(a,1x,3(f10.6,1x))') 'Icart',tmp(:)
      tmp(:)=zero
      call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,xred_center(:,FromIdeal2Average(iatom)),1,0.d0,tmp(:),1)
      write(31,'(a,1x,3(f10.6,1x))') 'Ccart',tmp(:)
    end do  
    close(31)
  end if
  ABI_FREE(xred_center)

!FB! Average distances between atoms --> distance_average
!FB  ABI_MALLOC(distance_average,(InVar%natom,InVar%natom,4))      ; distance_average(:,:,:)=0.d0
!FB  do eatom=1,InVar%natom
!FB    do fatom=1,InVar%natom
!FB      tmp(:)=xred_center(:,FromIdeal2Average(fatom))-xred_center(:,FromIdeal2Average(eatom))
!FB      call tdep_make_inbox(tmp,1,1d-3)
!FB      call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,tmp(:),1,0.d0,distance_average(eatom,fatom,2:4),1)
!FB      do ii=1,3
!FB!       Remove the rounding errors before writing (for non regression testing purposes)
!FB        if (abs(distance_average(eatom,fatom,ii+1)).lt.tol8) distance_average(eatom,fatom,ii+1)=zero
!FB        distance_average(eatom,fatom,1)=distance_average(eatom,fatom,1)+(distance_average(eatom,fatom,ii+1))**2
!FB      end do
!FB      distance_average(eatom,fatom,1)=distance_average(eatom,fatom,1)**0.5
!FB    end do  
!FB  end do  
!FB  ABI_FREE(xred_center)
!FB  ABI_FREE(distance_average)

!====================================================================================
!====================== END OF REDUCED COORDINATES ==================================
!====================================================================================
! a/ Get cartesian coordinates from reduced ones
! b/ Compute ucart and fcart tabs
! c/ The atoms are sorted according the IDEAL arrangement
!    The correspondance function is contained in: FromIdeal2Average
!    WARNING : Consequently the arrangement of the xcart* tabs is not modified. 
  write(InVar%stdout,*)' Compute cartesian coordinates and forces...'
  ABI_MALLOC(xcart        ,(3,InVar%natom,InVar%my_nstep)); xcart(:,:,:)=0.d0
  ABI_MALLOC(xcart_ideal  ,(3,InVar%natom))               ; xcart_ideal(:,:)=0.d0
  ABI_MALLOC(xcart_average,(3,InVar%natom))               ; xcart_average(:,:)=0.d0
  ABI_MALLOC(ucart_tmp    ,(3,InVar%natom,InVar%my_nstep)); ucart_tmp(:,:,:)=0.d0
  do iatom=1,InVar%natom
    call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,xred_ideal  (:,iatom),1,0.d0,xcart_ideal  (:,iatom),1)
    call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,xred_average(:,iatom),1,0.d0,xcart_average(:,iatom),1)
    do iatcell=1,InVar%natom_unitcell
      call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,Rlatt_red(:,iatcell,iatom),1,0.d0,Rlatt_cart(:,iatcell,iatom),1)
    end do  
  end do
  do istep=1,InVar%my_nstep
    do iatom=1,InVar%natom
      call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,InVar%xred(:,FromIdeal2Average(iatom),istep),&
&       1,0.d0,xcart(:,FromIdeal2Average(iatom),istep),1)
      if (InVar%Use_ideal_positions.eq.0) then
        ucart_tmp(:,iatom,istep)=xcart(:,FromIdeal2Average(iatom),istep)-xcart_average(:,FromIdeal2Average(iatom))
      else
        ucart_tmp(:,iatom,istep)=xcart(:,FromIdeal2Average(iatom),istep)-xcart_ideal  (:,iatom)
      end if  
    end do
  end do
  ABI_FREE(xred_average)
  ABI_FREE(xcart)
  ABI_FREE(xcart_ideal)
  ABI_FREE(xcart_average)

! Rearrangement of the fcart tabs in column --> Forces_MD
  ABI_MALLOC(fcart_tmp,(3,InVar%natom,InVar%my_nstep)); fcart_tmp(:,:,:)=0.d0
  do istep=1,InVar%my_nstep
    do iatom=1,InVar%natom
      fcart_tmp(:,iatom,istep)=InVar%fcart(:,FromIdeal2Average(iatom),istep)
    end do  
  end do  
  do istep=1,InVar%my_nstep
    do jatom=1,InVar%natom
      do ii=1,3 
        Forces_MD(ii+3*(jatom-1)+3*InVar%natom*(istep-1))=fcart_tmp(ii,jatom,istep)
        ucart(ii,jatom,istep)=ucart_tmp(ii,jatom,istep)
      enddo  
    enddo
  enddo  
  ABI_FREE(FromIdeal2Average)
  ABI_FREE(ucart_tmp)
  ABI_FREE(fcart_tmp)

! Define Rlatt4dos, fulfilling the definition of mkphdos (ABINIT routine)
  do ii=1,3
    rprimd_MD_tmp(ii,:)=Lattice%rprimd_MD(ii,:)/Lattice%acell_unitcell(ii)
  end do  
  do iatom=1,InVar%natom
    do iatcell=1,InVar%natom_unitcell
      call DGEMV('T',3,3,1.d0,rprimd_MD_tmp(:,:),3,Rlatt_red(:,iatcell,iatom),1,0.d0,Rlatt4dos(:,iatcell,iatom),1)
    end do  
  end do

! Find the symetry operation between 2 atoms 
  call tdep_SearchS_1at(InVar,Lattice,MPIdata,Sym,xred_ideal)
  ABI_MALLOC(InVar%xred_ideal,(3,InVar%natom)) ; InVar%xred_ideal(:,:)=0.d0
  InVar%xred_ideal(:,:)=xred_ideal(:,:)
  ABI_FREE(xred_ideal)
  ABI_FREE(Rlatt_red)

 end subroutine tdep_MatchIdeal2Average

!====================================================================================================
 subroutine tdep_calc_model(Forces_MD,Forces_TDEP,InVar,MPIdata,Phi1Ui,Phi2UiUj,&
&                           Phi3UiUjUk,Phi4UiUjUkUl,U0) 

  implicit none 

  type(Input_Variables_type),intent(in) :: InVar
  type(MPI_enreg_type), intent(in) :: MPIdata
  double precision, intent(in)  :: Forces_MD(3*InVar%natom*InVar%my_nstep)
  double precision, intent(in)  :: Forces_TDEP(3*InVar%natom*InVar%my_nstep)
  double precision, intent(out) :: U0
  double precision, intent(in)  :: Phi1Ui(InVar%my_nstep)
  double precision, intent(in)  :: Phi2UiUj(InVar%my_nstep)
  double precision, intent(in)  :: Phi3UiUjUk(InVar%my_nstep)
  double precision, intent(in)  :: Phi4UiUjUkUl(InVar%my_nstep)
  
  integer :: ii,jj,kk,istep,iatom,jatom,katom,my_istep
  double precision :: Delta_F2,Delta_U,Delta_U2
  double precision :: sigma,U_1,U_2,U_3,U_4,UMD
  double precision, allocatable :: tmp(:),Phi_tot(:)
  double precision, allocatable :: U_MD(:),U_TDEP(:)
  integer :: ierr

  ierr = 0

  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '######################### Energies, errors,...  #############################'
  write(InVar%stdout,*) '#############################################################################'
  
! Compute U0, U_TDEP, Delta_U and write them in the data.out file
  write(InVar%stdout,'(a)') ' Thermodynamic quantities and convergence parameters of THE MODEL,' 
  write(InVar%stdout,'(a)') '      as a function of the step number (energies in eV/atom and forces in Ha/bohr) :'
  if (InVar%Order.eq.4) then
    write(InVar%stdout,'(a)') ' <U_TDEP> = U_0 + U_1 + U_2 + U_3 + U_4'
    write(InVar%stdout,'(a)') '       with U_0 = < U_MD - sum_i Phi1 ui - 1/2 sum_ij Phi2 ui uj - 1/6 sum_ijk Phi3 ui uj uk - 1/24 sum_ijkl Phi4 ui uj uk ul >'
    write(InVar%stdout,'(a)') '        and U_1 = <      sum_i    Phi1 ui >'
    write(InVar%stdout,'(a)') '        and U_2 = < 1/2  sum_ij   Phi2 ui uj >'
    write(InVar%stdout,'(a)') '        and U_3 = < 1/6  sum_ijk  Phi3 ui uj uk >'
    write(InVar%stdout,'(a)') '        and U_4 = < 1/24 sum_ijkl Phi4 ui uj uk ul >'
  else if (InVar%Order.eq.3) then
    write(InVar%stdout,'(a)') ' <U_TDEP> = U_0 + U_1 + U_2 + U_3'
    write(InVar%stdout,'(a)') '       with U_0 = < U_MD - sum_i Phi1 ui - 1/2 sum_ij Phi2 ui uj - 1/6 sum_ijk Phi3 ui uj uk >'
    write(InVar%stdout,'(a)') '        and U_1 = <      sum_i    Phi1 ui >'
    write(InVar%stdout,'(a)') '        and U_2 = < 1/2  sum_ij   Phi2 ui uj >'
    write(InVar%stdout,'(a)') '        and U_3 = < 1/6  sum_ijk  Phi3 ui uj uk >'
  else  
    write(InVar%stdout,'(a)') ' <U_TDEP> = U_0 + U_1 + U_2'
    write(InVar%stdout,'(a)') '       with U_0 = < U_MD - sum_i Phi1 ui - 1/2 sum_ij Phi2 ui uj >'
    write(InVar%stdout,'(a)') '        and U_1 = <      sum_i    Phi1 ui >'
    write(InVar%stdout,'(a)') '        and U_2 = < 1/2  sum_ij   Phi2 ui uj >'
  end if  
  write(InVar%stdout,'(a)') '  Delta_U =   < U_MD - U_TDEP > '
  write(InVar%stdout,'(a)') '  Delta_U2= (< (U_MD - U_TDEP)^2 >)**0.5 '
  write(InVar%stdout,'(a)') '  Delta_F2= (< (F_MD - F_TDEP)^2 >)**0.5 '
  write(InVar%stdout,'(a)') '  Sigma   = (< (F_MD - F_TDEP)^2 >/<F_MD**2>)**0.5 '
  if (InVar%Order.eq.4) then
    write(InVar%stdout,'(2a)') '     <U_MD>            U_0              U_1              U_2  ',&
&     '            U_3              U_4            Delta_U          Delta_U2          Delta_F2          Sigma'
  else if (InVar%Order.eq.3) then
    write(InVar%stdout,'(2a)') '     <U_MD>            U_0              U_1              U_2  ',&
&     '            U_3            Delta_U          Delta_U2          Delta_F2          Sigma'
  else
    write(InVar%stdout,'(2a)') '     <U_MD>            U_0              U_1              U_2  ',&
&     '          Delta_U          Delta_U2          Delta_F2          Sigma'
  end if  

! Compute eucledian distance for forces    
  ABI_MALLOC(tmp,(11))                  ; tmp(:)   =0.d0
  do istep=1,InVar%my_nstep
    do iatom=1,InVar%natom
      do ii=1,3
        tmp(4)=tmp(4)+(Forces_MD(ii+3*(iatom-1)+3*InVar%natom*(istep-1))&
&               -Forces_TDEP(ii+3*(iatom-1)+3*InVar%natom*(istep-1)))**2
        tmp(5)=tmp(5)+Forces_MD(ii+3*(iatom-1)+3*InVar%natom*(istep-1))**2
      end do
    end do  
  end do
! Compute energies
  ABI_MALLOC(U_TDEP,(InVar%nstep_tot))  ; U_TDEP(:)=0.d0
  ABI_MALLOC(U_MD,  (InVar%nstep_tot))  ; U_MD(:)  =0.d0
  ABI_MALLOC(Phi_tot,(MPIdata%my_nstep)); Phi_tot(:)=0.d0
  do istep=1,InVar%my_nstep
    tmp(7) =tmp(7) +InVar%etot(istep)
    tmp(10)=tmp(10)+Phi1Ui(istep)
    tmp(6) =tmp(6) +Phi2UiUj(istep)
    tmp(8) =tmp(8) +Phi3UiUjUk(istep)
    tmp(11)=tmp(11)+Phi4UiUjUkUl(istep)
  end do  
  call xmpi_sum(tmp,MPIdata%comm_step,ierr)
  tmp(1) = tmp(7)-tmp(10)-tmp(6)-tmp(8)-tmp(11)
  Phi_tot(:)=tmp(1)/real(InVar%nstep_tot)+Phi1Ui(:)+Phi2UiUj(:)+Phi3UiUjUk(:)+Phi4UiUjUkUl(:)
  call xmpi_gatherv(Phi_tot,InVar%my_nstep,U_TDEP,MPIdata%nstep_all,MPIdata%shft_step,&
&                   MPIdata%master,MPIdata%comm_step,ierr)
  call xmpi_gatherv(InVar%etot,InVar%my_nstep,U_MD,MPIdata%nstep_all,MPIdata%shft_step,&
&                   MPIdata%master,MPIdata%comm_step,ierr)
  do istep=1,InVar%nstep_tot
    tmp(2) =tmp(2) + (U_MD(istep)-U_TDEP(istep))
    tmp(9) =tmp(9) + (U_MD(istep)-U_TDEP(istep))**2
  end do
  U0       =tmp(1) /real(InVar%nstep_tot*InVar%natom)
  UMD      =tmp(7) /real(InVar%nstep_tot*InVar%natom)
  U_1      =tmp(10)/real(InVar%nstep_tot*InVar%natom)
  U_2      =tmp(6) /real(InVar%nstep_tot*InVar%natom)
  U_3      =tmp(8) /real(InVar%nstep_tot*InVar%natom)
  U_4      =tmp(11)/real(InVar%nstep_tot*InVar%natom)
  Delta_U  =tmp(2) /real(InVar%nstep_tot*InVar%natom)
  Delta_U2 =tmp(9) /real(InVar%nstep_tot*InVar%natom)
  Delta_F2 =tmp(4) /real(InVar%nstep_tot*InVar%natom*3)
  sigma    =dsqrt(tmp(4)/tmp(5))
  if (InVar%Order.eq.4) then
    write(InVar%stdout,'(10(f12.5,5x))') UMD*Ha_eV,U0*Ha_eV,U_1*Ha_eV,U_2*Ha_eV,U_3*Ha_eV,U_4*Ha_eV,&
&     Delta_U*Ha_eV,Delta_U2**0.5*Ha_eV,Delta_F2**0.5,sigma
  else if (InVar%Order.eq.3) then
    write(InVar%stdout,'(9(f12.5,5x))') UMD*Ha_eV,U0*Ha_eV,U_1*Ha_eV,U_2*Ha_eV,U_3*Ha_eV,&
&     Delta_U*Ha_eV,Delta_U2**0.5*Ha_eV,Delta_F2**0.5,sigma
  else
    write(InVar%stdout,'(8(f12.5,5x))') UMD*Ha_eV,U0*Ha_eV,U_1*Ha_eV,U_2*Ha_eV,&
&     Delta_U*Ha_eV,Delta_U2**0.5*Ha_eV,Delta_F2**0.5,sigma
  endif 
  ABI_FREE(tmp)
  write(InVar%stdout,'(a)') ' TODO : write all the data in etotMDvsTDEP.dat & fcartMDvsTDEP.dat '
  write(InVar%stdout,'(a,1x,f12.5)') ' NOTE : in the harmonic and classical limit (T>>T_Debye), U_2=3/2*kB*T=',&
&   3.d0/2.d0*kb_HaK*Ha_eV*InVar%temperature

! Write : i) (U_TDEP vs U_MD) in etotMDvsTDEP.dat
!        ii) (Forces_TDEP vs Forces_MD) in fcartMDvsTDEP.dat
  write(InVar%stdout,'(a)') ' '
  if (MPIdata%iam_master) then
    open(unit=32,file=trim(InVar%output_prefix)//'etotMDvsTDEP.dat')
    open(unit=33,file=trim(InVar%output_prefix)//'fcartMDvsTDEP.dat')
    write(32,'(a)') '#   Istep      U_MD(Ha)         U_TDEP(Ha)'
    write(33,'(a)') '# Forces_MD(Ha/bohr) Forces_TDEP(Ha/bohr)'
    do istep=1,InVar%nstep_tot
      write(32,'(i6,1x,2(f17.6,1x))') istep,U_MD(istep),U_TDEP(istep)
    end do  
    do istep=1,InVar%my_nstep
      do iatom=1,InVar%natom
        do ii=1,3
          write(33,'(2(f17.10,1x))') Forces_MD  (ii+3*(iatom-1)+3*InVar%natom*(istep-1)),&
&                                    Forces_TDEP(ii+3*(iatom-1)+3*InVar%natom*(istep-1))
        end do
      end do  
    end do  
    close(32)
    close(33)
  end if  
  ABI_FREE(U_MD)
  ABI_FREE(U_TDEP)
  ABI_FREE(Phi_tot)

 end subroutine tdep_calc_model

!====================================================================================================
subroutine tdep_calc_nbcoeff(distance,iatcell,InVar,ishell,jatom,katom,latom,MPIdata,&
&                            ncoeff,norder,nshell,order,proj,Sym)

  implicit none


  integer,intent(in) :: iatcell,ishell,jatom,katom,latom,nshell,order,norder
  integer,intent(out) :: ncoeff
  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(MPI_enreg_type), intent(in) :: MPIdata
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(out) :: proj(norder,norder,nshell)

  integer :: ii,jj,kk,ll,isym,LWORK,INFO,const_tot,itemp,nconst_perm,nconst_loc
  integer :: ncount,icoeff,jatcell,katcell,latcell,mu,nu,xi,zeta
  integer :: inv,watom,xatom,yatom,zatom,isyminv,nsyminv,facorder
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
      jatcell=mod(jatom-1,InVar%natom_unitcell)+1
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
      if (inv==1) then ; watom=iatcell ; xatom=jatom   ; yatom=katom   ; endif !\Psi_ijk
      if (inv==2) then ; watom=iatcell ; xatom=katom   ; yatom=jatom   ; endif !\Psi_ikj
      if (inv==3) then ; watom=jatom   ; xatom=iatcell ; yatom=katom   ; endif !\Psi_jik
      if (inv==4) then ; watom=jatom   ; xatom=katom   ; yatom=iatcell ; endif !\Psi_jki
      if (inv==5) then ; watom=katom   ; xatom=iatcell ; yatom=jatom   ; endif !\Psi_kij
      if (inv==6) then ; watom=katom   ; xatom=jatom   ; yatom=iatcell ; endif !\Psi_kji
      do ii=1,3
        do jj=1,3
          vect_trial1(ii)=vect_trial1(ii)+Sym%S_ref(ii,jj,isym,1)*distance(watom,xatom,jj+1)
          vect_trial2(ii)=vect_trial2(ii)+Sym%S_ref(ii,jj,isym,1)*distance(xatom,yatom,jj+1)
          vect_trial3(ii)=vect_trial3(ii)+Sym%S_ref(ii,jj,isym,1)*distance(yatom,watom,jj+1)
        end do  
      end do
      jatcell=mod(jatom-1,InVar%natom_unitcell)+1
      katcell=mod(katom-1,InVar%natom_unitcell)+1
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
      if (inv==1) then ; watom=iatcell ; xatom=jatom   ; yatom=katom   ; zatom=latom   ; endif !\Psi_ijkl
      if (inv==2) then ; watom=iatcell ; xatom=katom   ; yatom=jatom   ; zatom=latom   ; endif !\Psi_ikjl
      if (inv==3) then ; watom=jatom   ; xatom=iatcell ; yatom=katom   ; zatom=latom   ; endif !\Psi_jikl
      if (inv==4) then ; watom=jatom   ; xatom=katom   ; yatom=iatcell ; zatom=latom   ; endif !\Psi_jkil
      if (inv==5) then ; watom=katom   ; xatom=iatcell ; yatom=jatom   ; zatom=latom   ; endif !\Psi_kijl
      if (inv==6) then ; watom=katom   ; xatom=jatom   ; yatom=iatcell ; zatom=latom   ; endif !\Psi_kjil

      if (inv==7 ) then ; watom=iatcell ; xatom=jatom   ; yatom=latom   ; zatom=katom   ; endif !\Psi_ijlk
      if (inv==8 ) then ; watom=iatcell ; xatom=katom   ; yatom=latom   ; zatom=jatom   ; endif !\Psi_iklj
      if (inv==9 ) then ; watom=jatom   ; xatom=iatcell ; yatom=latom   ; zatom=katom   ; endif !\Psi_jilk
      if (inv==10) then ; watom=jatom   ; xatom=katom   ; yatom=latom   ; zatom=iatcell ; endif !\Psi_jkli
      if (inv==11) then ; watom=katom   ; xatom=iatcell ; yatom=latom   ; zatom=jatom   ; endif !\Psi_kilj
      if (inv==12) then ; watom=katom   ; xatom=jatom   ; yatom=latom   ; zatom=iatcell ; endif !\Psi_kjli

      if (inv==13) then ; watom=iatcell ; xatom=latom   ; yatom=jatom   ; zatom=katom   ; endif !\Psi_iljk
      if (inv==14) then ; watom=iatcell ; xatom=latom   ; yatom=katom   ; zatom=jatom   ; endif !\Psi_ilkj
      if (inv==15) then ; watom=jatom   ; xatom=latom   ; yatom=iatcell ; zatom=katom   ; endif !\Psi_jlik
      if (inv==16) then ; watom=jatom   ; xatom=latom   ; yatom=katom   ; zatom=iatcell ; endif !\Psi_jlki
      if (inv==17) then ; watom=katom   ; xatom=latom   ; yatom=iatcell ; zatom=jatom   ; endif !\Psi_klij
      if (inv==18) then ; watom=katom   ; xatom=latom   ; yatom=jatom   ; zatom=iatcell ; endif !\Psi_klji

      if (inv==19) then ; watom=latom   ; xatom=iatcell ; yatom=jatom   ; zatom=katom   ; endif !\Psi_lijk
      if (inv==20) then ; watom=latom   ; xatom=iatcell ; yatom=katom   ; zatom=jatom   ; endif !\Psi_likj
      if (inv==21) then ; watom=latom   ; xatom=jatom   ; yatom=iatcell ; zatom=katom   ; endif !\Psi_ljik
      if (inv==22) then ; watom=latom   ; xatom=jatom   ; yatom=katom   ; zatom=iatcell ; endif !\Psi_ljki
      if (inv==23) then ; watom=latom   ; xatom=katom   ; yatom=iatcell ; zatom=jatom   ; endif !\Psi_lkij
      if (inv==24) then ; watom=latom   ; xatom=katom   ; yatom=jatom   ; zatom=iatcell ; endif !\Psi_lkji

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
      jatcell=mod(jatom-1,InVar%natom_unitcell)+1
      katcell=mod(katom-1,InVar%natom_unitcell)+1
      latcell=mod(latom-1,InVar%natom_unitcell)+1
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
          if (inv==2 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (i,k,j,l) for isym=',isym !\Psi_ikjl
          if (inv==3 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (j,i,k,l) for isym=',isym !\Psi_jikl
          if (inv==4 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (j,k,i,l) for isym=',isym !\Psi_jkil
          if (inv==5 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (k,i,j,l) for isym=',isym !\Psi_kijl
          if (inv==6 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (k,j,i,l) for isym=',isym !\Psi_kjil
  
          if (inv==7 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (i,j,l,k) for isym=',isym !\Psi_ijlk
          if (inv==8 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (i,k,l,j) for isym=',isym !\Psi_iklj
          if (inv==9 ) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (j,i,l,k) for isym=',isym !\Psi_jilk
          if (inv==10) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (j,k,l,i) for isym=',isym !\Psi_jkli
          if (inv==11) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (k,i,l,j) for isym=',isym !\Psi_kilj
          if (inv==12) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (k,j,l,i) for isym=',isym !\Psi_kjli

          if (inv==13) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (i,l,j,k) for isym=',isym !\Psi_iljk
          if (inv==14) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (i,l,k,j) for isym=',isym !\Psi_ilkj
          if (inv==15) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (j,l,i,k) for isym=',isym !\Psi_jlik
          if (inv==16) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (j,l,k,i) for isym=',isym !\Psi_jlki
          if (inv==17) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (k,l,i,j) for isym=',isym !\Psi_klij
          if (inv==18) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (k,l,j,i) for isym=',isym !\Psi_klji

          if (inv==19) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (l,i,j,k) for isym=',isym !\Psi_lijk
          if (inv==20) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (l,i,k,j) for isym=',isym !\Psi_likj
          if (inv==21) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (l,j,i,k) for isym=',isym !\Psi_ljik
          if (inv==22) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (l,j,k,i) for isym=',isym !\Psi_ljki
          if (inv==23) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (l,k,i,j) for isym=',isym !\Psi_lkij
          if (inv==24) write(16,'(a,1x,i3)')'===========The bond is reversed with (i,j,k,l) --> (l,k,j,i) for isym=',isym !\Psi_lkji
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
!   These ones could be complex!!!!
    ok=.false.
    do ii=1,3
!FB      write(16,*)'  For eigenvalue number',ii
!FB      write(16,*)'     The eigenvalue is:',eigenvalues(ii)
!FB      write(16,*)'     The eigenvector is:'
!FB      write(16,'(2(f16.12,1x))') eigenvectors(1,ii)
!FB      write(16,'(2(f16.12,1x))') eigenvectors(2,ii)
!FB      write(16,'(2(f16.12,1x))') eigenvectors(3,ii)
      if ((aimag(eigenvalues(1)).ne.0).or.(aimag(eigenvalues(2)).ne.0).or.(aimag(eigenvalues(3)).ne.0)) then
	ok=.true.
      end if  
    end do  
    if (ok.and.MPIdata%iam_master) write(16,'(a)') '            WARNING: THERE IS COMPLEX EIGENVALUES'

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
        if ((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6)) cycle
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
          if (((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6).and.(inv==1)).or.&
&             ((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6).and.(inv==2).and.(ii==jj))) cycle
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
                MSG_BUG('This symetry is neither Keptinvariant nor Reversed')
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
            if (((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6).and.(inv==1)).or.&
&               ((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6).and.(inv==2).and.(jj==kk)).or.&
&               ((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6).and.(inv==3).and.(ii==jj)).or.&
&               ((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6).and.(inv==4).and.(ii==jj).and.(jj==kk)).or.&
&               ((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6).and.(inv==5).and.(ii==jj).and.(jj==kk)).or.&
&               ((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6).and.(inv==6).and.(ii==kk))) cycle
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
                    MSG_BUG('This symetry is neither Keptinvariant nor Reversed')
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
              if ((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6)) then
                if ((inv==1 )                                       .or.& !\Psi_ijkl
&                 ((inv==2 ).and.(jj==kk))                          .or.& !\Psi_ikjl
&                 ((inv==3 ).and.(ii==jj))                          .or.& !\Psi_jikl
&                 ((inv==4 ).and.(ii==jj).and.(jj==kk))             .or.& !\Psi_jkil
&                 ((inv==5 ).and.(ii==jj).and.(jj==kk))             .or.& !\Psi_kijl
&                 ((inv==6 ).and.(ii==kk))                          .or.& !\Psi_kjil

&                 ((inv==7 ).and.(kk==ll))                          .or.& !\Psi_ijlk
&                 ((inv==8 ).and.(jj==kk).and.(kk==ll))             .or.& !\Psi_iklj
&                 ((inv==9 ).and.(ii==jj).and.(kk==ll))             .or.& !\Psi_jilk
&                 ((inv==10).and.(ii==jj).and.(jj==kk).and.(kk==ll)).or.& !\Psi_jkli
&                 ((inv==11).and.(ii==jj).and.(jj==kk).and.(kk==ll)).or.& !\Psi_kilj
&                 ((inv==12).and.(ii==kk).and.(kk==ll))             .or.& !\Psi_kjli

&                 ((inv==13).and.(jj==kk).and.(kk==ll))             .or.& !\Psi_iljk
&                 ((inv==14).and.(jj==ll))                          .or.& !\Psi_ilkj
&                 ((inv==15).and.(ii==jj).and.(jj==kk).and.(kk==ll)).or.& !\Psi_jlik
&                 ((inv==16).and.(ii==jj).and.(jj==ll))             .or.& !\Psi_jlki
&                 ((inv==17).and.(ii==kk).and.(jj==ll))             .or.& !\Psi_klij
&                 ((inv==18).and.(ii==jj).and.(jj==kk).and.(kk==ll)).or.& !\Psi_klji

&                 ((inv==19).and.(ii==jj).and.(jj==kk).and.(kk==ll)).or.& !\Psi_lijk
&                 ((inv==20).and.(ii==jj).and.(jj==ll))             .or.& !\Psi_likj
&                 ((inv==21).and.(ii==kk).and.(kk==ll))             .or.& !\Psi_ljik
&                 ((inv==22).and.(ii==ll))                          .or.& !\Psi_ljki
&                 ((inv==23).and.(ii==jj).and.(jj==kk).and.(kk==ll)).or.& !\Psi_lkij
&                 ((inv==24).and.(ii==ll).and.(jj==kk))) cycle            !\Psi_lkji
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
                      else ; MSG_BUG('This symetry is neither Keptinvariant nor Reversed')
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
      MSG_BUG('Only the first, second, third and fourth order are allowed')
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
          if (abs( real(tab_vec(ii,kk))).lt.1.d-8) tab_vec(ii,kk)=dcmplx(zero,aimag(tab_vec(ii,kk)))
          if (abs(aimag(tab_vec(ii,kk))).lt.1.d-8) tab_vec(ii,kk)=dcmplx( real(tab_vec(ii,kk)),zero)
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
      return
    else if (order==1.and.(iconst(isyminv).gt.3)) then
      MSG_BUG(' First order : There are more than 3 constraints')
    end if
    if (order==2.and.(iconst(isyminv).gt.8)) then
      MSG_BUG(' Second order : There are more than 8 constraints')
    end if
    if (order==3.and.(iconst(isyminv).gt.27)) then
      MSG_BUG(' Third order : There are more than 27 constraints')
    end if
    if (order==4.and.(iconst(isyminv).gt.81)) then
      MSG_BUG(' Fourth order : There are more than 81 constraints')
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
    if ((iatcell.eq.jatom).or.(iatcell.eq.katom).or.(iatcell.eq.latom).or.(jatom.eq.katom).or.(jatom.eq.latom).or.(katom.eq.latom)) then
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
    MSG_WARNING(message)
    proj(:,:,ishell)=0.d0
    do ii=1,norder
      proj(ii,ii,ishell)=1.d0
    end do
    ncoeff=norder
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
    MSG_BUG(message)
  end if  
  do ii=1,norder
    do jj=1,ncount
      if (abs( real(tab_vec(ii,jj))).lt.1.d-8) tab_vec(ii,jj)=dcmplx(zero,aimag(tab_vec(ii,jj)))
      if (abs(aimag(tab_vec(ii,jj))).lt.1.d-8) tab_vec(ii,jj)=dcmplx( real(tab_vec(ii,jj)),zero)
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
          if (abs( real(tab_vec(ii,kk))).lt.1.d-8) tab_vec(ii,kk)=dcmplx(zero,aimag(tab_vec(ii,kk)))
          if (abs(aimag(tab_vec(ii,kk))).lt.1.d-8) tab_vec(ii,kk)=dcmplx( real(tab_vec(ii,kk)),zero)
        end do
!FB      else
!FB        write(InVar%stdout,*)'One prod_scal equals zero'
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
    MSG_ERROR(' Order 2 : There are too many independent vectors')
  end if
  if (ncount.gt.27.and.order==3) then
    MSG_ERROR(' Order 3 : There are too many independent vectors')
  end if
  if (ncount.gt.81.and.order==4) then
    MSG_ERROR(' Order 4 : There are too many independent vectors')
  end if

! On cherche les (norder-ncount) vecteurs orthogonaux aux vecteurs non-nuls
! --> Orthogonalisation de Gram-Schmidt
  ABI_MALLOC(tab_vec,(norder,norder)); tab_vec(:,:)=czero
  do kk=1,norder
    if (kk.le.ncount) then
      tab_vec(:,kk)=temp(:,kk)
    else  
      do jj=1,norder
        call random_number(drandom)
        tab_vec(jj,kk)=dcmplx(drandom,zero)
      end do  
      do jj=1,kk-1
        prod_scal=sum( real(tab_vec(:,jj))* real(tab_vec(:,jj))+aimag(tab_vec(:,jj))*aimag(tab_vec(:,jj)))
        if (abs(prod_scal).gt.tol8) then
          tab_vec(:,kk)=tab_vec(:,kk)-sum(tab_vec(:,kk)*conjg(tab_vec(:,jj)))/prod_scal*tab_vec(:,jj)
          do ii=1,norder
            if (abs( real(tab_vec(ii,kk))).lt.1.d-8) tab_vec(ii,kk)=dcmplx(zero,aimag(tab_vec(ii,kk)))
            if (abs(aimag(tab_vec(ii,kk))).lt.1.d-8) tab_vec(ii,kk)=dcmplx( real(tab_vec(ii,kk)),zero)
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
    if ((abs(aimag(tab_vec(1,kk))).gt.1.d-6).or.&
&       (abs(aimag(tab_vec(1,kk))).gt.1.d-6).or.&
&       (abs(aimag(tab_vec(1,kk))).gt.1.d-6)) then
      MSG_ERROR('the constraint has an imaginary part')
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
!=====================================================================================================

end module m_tdep_utils
