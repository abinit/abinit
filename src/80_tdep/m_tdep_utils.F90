
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_utils
  
  use defs_basis
  use m_errors
  use m_abicore
  use m_tdep_latt,        only : Lattice_Variables_type, tdep_make_inbox
  use m_tdep_readwrite,   only : Input_Variables_type
  use m_tdep_sym,         only : Symetries_Variables_type, tdep_SearchS_1at

  implicit none

  type Coeff_Moore_type
    integer :: ntotcoeff
    integer :: ntotconst
    integer :: ncoeff1st
    integer :: ncoeff2nd
    integer :: ncoeff3rd
    integer :: nconst_1st
    integer :: nconst_2nd
    integer :: nconst_3rd
    integer :: nconst_rot2nd
    integer :: nconst_huang
    integer :: nconst_dynmat
    integer :: nconst_rot3rd
    integer :: nconst_asr3rd
    double precision, allocatable :: fcoeff(:,:)
  end type Coeff_Moore_type

 public :: tdep_calc_MoorePenrose
 public :: tdep_MatchIdeal2Average
 public :: tdep_calc_model
 public :: tdep_calc_nbcoeff

contains

!=====================================================================================================
 subroutine tdep_calc_MoorePenrose(Forces_MD,CoeffMoore,InVar,IFC_coeff)

  implicit none 

  type(Input_Variables_type),intent(in) :: InVar
  type(Coeff_Moore_type), intent(inout) :: CoeffMoore
  double precision, intent(in)  :: Forces_MD(3*InVar%natom*InVar%nstep)
  double precision, intent(out)  :: IFC_coeff(CoeffMoore%ntotcoeff,1)
  
  integer :: ii,LWORK,INFO,ntotcoeff,ntotconst,firstdim
  integer, allocatable :: IWORK(:)
  double precision, allocatable :: WORK(:),pseudo_inverse(:,:),sigma(:),matU(:,:)
  double precision, allocatable :: pseudo_sigma(:,:),transmatV(:,:),tmp1(:,:),fcart_tmp(:,:)

  write(InVar%stdout,*) '################### And compute the pseudo-inverse ##########################'
  write(InVar%stdout,*) '#############################################################################'

  ntotcoeff=CoeffMoore%ntotcoeff
  ntotconst=CoeffMoore%ntotconst
  firstdim=3*InVar%natom*InVar%nstep+ntotconst

  write(InVar%stdout,*) ' Singular value decomposition...'
  ABI_MALLOC(sigma,(ntotcoeff)) ; sigma(:)=0.d0 
  ABI_MALLOC(matU,(firstdim,ntotcoeff)) ; matU(:,:)=0.d0
  ABI_MALLOC(transmatV,(ntotcoeff,ntotcoeff)) ; transmatV(:,:)=0.d0
  LWORK=3*(ntotcoeff)**2+max(firstdim,4*(ntotcoeff)**2+4*(ntotcoeff))
  ABI_MALLOC(WORK,(LWORK)) ; WORK(:)=0.d0
  ABI_MALLOC(IWORK,(8*(ntotcoeff))) ; IWORK(:)=0
  call DGESDD('S',firstdim,ntotcoeff,CoeffMoore%fcoeff(1:firstdim,:),firstdim,&
&   sigma,matU,firstdim,transmatV,ntotcoeff,WORK,LWORK,IWORK,INFO) 
  ABI_FREE(CoeffMoore%fcoeff)
  ABI_FREE(WORK)
  ABI_FREE(IWORK)

  ABI_MALLOC(pseudo_sigma,(ntotcoeff,ntotcoeff)) ; pseudo_sigma(:,:)=0.d0
  sigma(:)=1.d0/sigma(:)
  write(InVar%stdout,*) ' The eigenvalues are:'
  do ii=1,ntotcoeff
    if (sigma(ii).lt.1.d8) then
      pseudo_sigma(ii,ii)=sigma(ii)
    end if  
    write(InVar%stdout,'(1x,i4,1x,f15.10)') ii,sigma(ii)
  end do
  write(InVar%stdout,'(a,1x,f15.10)')'  condition number=',maxval(sigma(:))/minval(sigma(:))
  ABI_FREE(sigma)
  
  write(InVar%stdout,*) ' Calculation of the pseudo-inverse...'
  ABI_MALLOC(tmp1,(ntotcoeff,firstdim)) ; tmp1(:,:)=0.d0
  call DGEMM('N','T',ntotcoeff,firstdim,ntotcoeff,1.d0,pseudo_sigma,ntotcoeff,matU,&
&   firstdim,1.d0,tmp1,ntotcoeff)
  ABI_FREE(matU)
  ABI_FREE(pseudo_sigma)
  ABI_MALLOC(pseudo_inverse,(ntotcoeff,firstdim)) ; pseudo_inverse(:,:)=0.d0
  call DGEMM('T','N',ntotcoeff,firstdim,ntotcoeff,1.d0,transmatV,ntotcoeff,&
&   tmp1,ntotcoeff,1.d0,pseudo_inverse,ntotcoeff)
  ABI_FREE(tmp1)
  ABI_FREE(transmatV)
  ABI_MALLOC(fcart_tmp,(firstdim,1)); fcart_tmp(:,:)=zero  
  fcart_tmp(1:3*InVar%natom*InVar%nstep,1)=Forces_MD(:)
! NOTE, we have to solve F_ij = -\sum_j \Phi_ij u_j, so we add a minus sign to the pseudo_inverse  
  call DGEMM('N','N',ntotcoeff,1,firstdim,1.d0,-pseudo_inverse,ntotcoeff,fcart_tmp,&
&   firstdim,1.d0,IFC_coeff,ntotcoeff)
  ABI_FREE(fcart_tmp)
  ABI_FREE(pseudo_inverse)

 end subroutine tdep_calc_MoorePenrose

!====================================================================================================
 subroutine tdep_MatchIdeal2Average(distance,Forces_MD,InVar,Lattice,&
&                              Rlatt_cart,Rlatt4dos,Sym,ucart)

  implicit none 

  type(Input_Variables_type),intent(inout) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Symetries_Variables_type),intent(inout) :: Sym
  double precision, intent(out)  :: distance(InVar%natom,InVar%natom,4)
  double precision, intent(out)  :: Forces_MD(3*InVar%natom*InVar%nstep)
  double precision, intent(out)  :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)
  double precision, intent(out)  :: Rlatt4dos (3,InVar%natom_unitcell,InVar%natom)
  double precision, intent(out)  :: ucart(3,InVar%natom,InVar%nstep)
  
  integer :: ii,jj,kk,max_ijk,iatcell,jatcell,iatom,jatom,eatom,fatom,istep,jstep
  integer :: foo,foo2,atom_ref
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
  integer :: ierr

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
            ABI_ERROR('The number of atoms found in the bigbox exceeds natom' )
          end if  
          xred_ideal(:,iatom)=xred_tmp(:)
          call DGEMV('T',3,3,1.d0,Lattice%multiplicitym1(:,:),3,Rlatt(:),1,0.d0,Rlatt_red(:,1,iatom),1)
          iatom=iatom+1
        end do   
      end do
    end do
  end do

  if (iatom.lt.InVar%natom) then
    ABI_ERROR('The number of atoms found in the bigbox is lower than natom')
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
  do istep=1,InVar%nstep
    do iatom=1,InVar%natom
      xred_average(:,iatom)=xred_average(:,iatom)+InVar%xred(:,iatom,istep)
    end do
  end do
  xred_average(:,:)=xred_average(:,:)/real(InVar%nstep)

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
!FB      write(InVar%stdout,*) 'ATOM REF=',atom_ref
      ok=.false.
      exit
    else if (foo2.gt.InVar%natom_unitcell) then
      ABI_BUG(' Something wrong: WTF')
    endif
  end do  
  if (ok) then
    open(unit=31,file=trim(InVar%output_prefix)//'xred_average.xyz')
    do iatom=1,InVar%natom
      write(31,'(a,1x,3(f10.6,1x))') 'C',xred_center(:,iatom)
      write(31,'(a,1x,3(f10.6,1x))') 'I',xred_ideal (:,iatom)
    end do
    close(31)
    ABI_ERROR_NOSTOP('The basis of atoms written in input.in file does not appear in the MD trajectory',ierr)
    ABI_ERROR('Perhaps, you can adjust the tolerance (tolmotif)')
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
        do istep=1,InVar%nstep
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
  open(unit=31,file=trim(InVar%output_prefix)//'xred_average.xyz')
  write(31,'(i4)') InVar%natom*2
  write(31,'(i4)') 1
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
            do istep=1,InVar%nstep
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
            do istep=1,InVar%nstep
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
      close(31)
      do eatom=1,InVar%natom
        write(31,'(a,1x,3(f10.6,1x))') 'I',xred_ideal (:,eatom)
        write(31,'(a,1x,3(f10.6,1x))') 'C',xred_center(:,eatom)
      end do
      ABI_ERROR('Problem to find the average position')
    end if  
  end do

! WARNING: VERY IMPORTANT: The positions are displayed/sorted (and used in the following) 
! according to ideal positions xred_ideal. 
! --> In reduced coordinates
  do iatom=1,InVar%natom
    write(31,'(a,1x,3(f10.6,1x))') 'Ired',xred_ideal (:,iatom)
    write(31,'(a,1x,3(f10.6,1x))') 'Cred',xred_center(:,FromIdeal2Average(iatom))
  end do  
! --> In cartesian coordinates  
  do iatom=1,InVar%natom
    tmp(:)=zero
    call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,xred_ideal (:,iatom),1,0.d0,tmp(:),1)
    write(31,'(a,1x,3(f10.6,1x))') 'Icart',tmp(:)
    tmp(:)=zero
    call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,xred_center(:,FromIdeal2Average(iatom)),1,0.d0,tmp(:),1)
    write(31,'(a,1x,3(f10.6,1x))') 'Ccart',tmp(:)
  end do  
  close(31)
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
  ABI_MALLOC(xcart        ,(3,InVar%natom,InVar%nstep)); xcart(:,:,:)=0.d0
  ABI_MALLOC(xcart_ideal  ,(3,InVar%natom))            ; xcart_ideal(:,:)=0.d0
  ABI_MALLOC(xcart_average,(3,InVar%natom))            ; xcart_average(:,:)=0.d0
  ABI_MALLOC(ucart_tmp    ,(3,InVar%natom,InVar%nstep)); ucart_tmp(:,:,:)=0.d0
  do iatom=1,InVar%natom
    call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,xred_ideal  (:,iatom),1,0.d0,xcart_ideal  (:,iatom),1)
    call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,xred_average(:,iatom),1,0.d0,xcart_average(:,iatom),1)
    do iatcell=1,InVar%natom_unitcell
      call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,Rlatt_red(:,iatcell,iatom),1,0.d0,Rlatt_cart(:,iatcell,iatom),1)
    end do  
  end do
  do istep=1,InVar%nstep
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
  ABI_MALLOC(fcart_tmp,(3,InVar%natom,InVar%nstep)); fcart_tmp(:,:,:)=0.d0
  do istep=1,InVar%nstep
    do iatom=1,InVar%natom
      fcart_tmp(:,iatom,istep)=InVar%fcart(:,FromIdeal2Average(iatom),istep)
    end do  
  end do  
  jstep=0
  do istep=1,InVar%nstep
    jstep=jstep+1
    do jatom=1,InVar%natom
      do ii=1,3 
        Forces_MD(ii+3*(jatom-1)+3*InVar%natom*(jstep-1))=fcart_tmp(ii,jatom,istep)
        ucart(ii,jatom,jstep)=ucart_tmp(ii,jatom,istep)
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
  call tdep_SearchS_1at(InVar,Lattice,Sym,xred_ideal)
  ABI_MALLOC(InVar%xred_ideal,(3,InVar%natom)) ; InVar%xred_ideal(:,:)=0.d0
  InVar%xred_ideal(:,:)=xred_ideal(:,:)
  ABI_FREE(xred_ideal)
  ABI_FREE(Rlatt_red)

 end subroutine tdep_MatchIdeal2Average

!====================================================================================================
 subroutine tdep_calc_model(Free_Anh,Forces_MD,Forces_TDEP,InVar,Phij_NN,Pij_N,ucart,U0,&
&                 ftot3) !optional 

  implicit none 

  type(Input_Variables_type),intent(in) :: InVar
  double precision, intent(in)  :: Forces_MD(3*InVar%natom*InVar%nstep)
  double precision, intent(out) :: Forces_TDEP(3*InVar%natom*InVar%nstep)
  double precision, intent(in)  :: Phij_NN(3*InVar%natom,3*InVar%natom)
  double precision, intent(in)  :: Pij_N(3*InVar%natom)
  double precision, intent(in)  :: ucart(3,InVar%natom,InVar%nstep)
  double precision, intent(out) :: U0,Free_Anh
  double precision, intent(in),optional  :: ftot3(3*InVar%natom,InVar%nstep)
  
  integer :: ii,jj,istep,iatom,jatom,islice,nslice,istepmin,istepmax
  !integer :: nu,delta !alpha,beta,gama,lambda,
  double precision :: Delta_F2,Delta_U,Delta_U2
  double precision :: tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,sigma,U_1,U_2,U_3,UMD
  double precision, allocatable :: Fmean(:)
  double precision, allocatable :: U_MD(:),U_TDEP(:),PijUi(:),PhijUiUj(:),PsijUiUjUk(:),residualF(:)
  double precision, allocatable :: ftot2(:),ucart_blas(:)
  integer :: ierr

  ierr = 0;

  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '######################### Energies, errors,...  #############################'
  write(InVar%stdout,*) '#############################################################################'
  ABI_MALLOC(U_MD,(InVar%nstep)); U_MD(:)=0.d0
  do istep=1,InVar%nstep
    U_MD(istep)=InVar%etot(istep)
  end do  
  
! Compute Forces of the model TDEP
  ABI_MALLOC(PijUi     ,(InVar%nstep)); PijUi   (:)=0.d0 
  ABI_MALLOC(PhijUiUj  ,(InVar%nstep)); PhijUiUj(:)=0.d0 
  if (present(ftot3)) then 
    ABI_MALLOC(PsijUiUjUk,(InVar%nstep)); PsijUiUjUk(:)=0.d0 
  end if  
  ABI_MALLOC(ucart_blas,(3*InVar%natom))                      ; ucart_blas(:)=0.d0
  ABI_MALLOC(ftot2     ,(3*InVar%natom))                      ; ftot2     (:)=0.d0
  do istep=1,InVar%nstep
    ftot2(:)=0.d0 ; ucart_blas(:)=0.d0 
    do jatom=1,InVar%natom
      do jj=1,3
        ucart_blas(3*(jatom-1)+jj)=ucart(jj,jatom,istep)
      end do
    end do
    call DGEMM('N','N',3*InVar%natom,1,3*InVar%natom,1.d0,Phij_NN,3*InVar%natom,ucart_blas,3*InVar%natom,0.d0,ftot2,3*InVar%natom)
    call DGEMM('T','N',1,1,3*InVar%natom,1./2.d0,ftot2,3*InVar%natom,ucart_blas,3*InVar%natom,0.d0,PhijUiUj(istep),3*InVar%natom)
    PijUi(istep)=sum(Pij_N(:)*ucart_blas(:))
    Forces_TDEP(3*InVar%natom*(istep-1)+1:3*InVar%natom*istep)=-Pij_N(:)-ftot2(:)
    if (present(ftot3)) then
      call DGEMM('T','N',1,1,3*InVar%natom,1./6.d0,ftot3(:,istep),3*InVar%natom,ucart_blas,&
&       3*InVar%natom,0.d0,PsijUiUjUk(istep),3*InVar%natom)
      Forces_TDEP(3*InVar%natom*(istep-1)+1:3*InVar%natom*istep)=-Pij_N(:)-ftot2(:)-ftot3(:,istep)/2.d0
    end if !Psij
  end do !istep  
  ABI_FREE(ftot2)
  ABI_FREE(ucart_blas)

! Compute U0, U_TDEP, Delta_U and write them in the data.out file
  write(InVar%stdout,'(a)') ' Thermodynamic quantities and convergence parameters of THE MODEL,' 
  write(InVar%stdout,'(a)') '      as a function of the step number (energies in eV/atom and forces in Ha/bohr) :'
  if (present(ftot3)) then
    write(InVar%stdout,'(a)') ' <U_TDEP> = U_0 + U_1 + U_2 + U_3'
    write(InVar%stdout,'(a)') '       with U_0 = < U_MD - sum_i Pii ui - 1/2 sum_ij Phij ui uj - 1/6 sum_ij Psijk ui uj uk >'
    write(InVar%stdout,'(a)') '        and U_1 = < sum_i Pii ui >'
    write(InVar%stdout,'(a)') '        and U_2 = < 1/2 sum_ij Phij  ui uj >'
    write(InVar%stdout,'(a)') '        and U_3 = < 1/6 sum_ij Psijk ui uj uk >'
  else  
    write(InVar%stdout,'(a)') ' <U_TDEP> = U_0 +  U_1 +U_2'
    write(InVar%stdout,'(a)') '       with U_0 = < U_MD - sum_i Pii ui - 1/2 sum_ij Phij ui uj >'
    write(InVar%stdout,'(a)') '        and U_1 = < sum_i Pii ui >'
    write(InVar%stdout,'(a)') '        and U_2 = < 1/2 sum_ij Phij  ui uj >'
  end if  
  write(InVar%stdout,'(a)') '  Delta_U = < U_MD-F_TDEP > '
  write(InVar%stdout,'(a)') '  Delta_U2= (< (U_MD-U_TDEP)^2 >)**0.5 '
  write(InVar%stdout,'(a)') '  Delta_F2= (< (F_MD-F_TDEP)^2 >)**0.5 '
  write(InVar%stdout,'(a)') '  Sigma   = (< (F_MD-F_TDEP)^2 >/<F_MD**2>)**0.5 '
  if (present(ftot3)) then
    write(InVar%stdout,'(2a)') ' Istep         <U_MD>            U_0              U_1              U_2  ',&
&     '            U_3            Delta_U          Delta_U2          Delta_F2          Sigma'
  else
    write(InVar%stdout,'(2a)') ' Istep         <U_MD>            U_0              U_1              U_2  ',&
&     '          Delta_U          Delta_U2          Delta_F2          Sigma'
  end if  
  ABI_MALLOC(U_TDEP,(InVar%nstep)); U_TDEP(:)=0.d0
  tmp0=zero; tmp1=zero; tmp2=zero; tmp3=zero; tmp4=zero; tmp5=zero; tmp6=zero; tmp7=zero; tmp8=zero; tmp9=zero
  nslice=10
  do islice=1,nslice
    istepmin=int(InVar%nstep*(islice-1)/nslice+1)
    istepmax=int(InVar%nstep*islice/nslice)
    do istep=istepmin,istepmax
      tmp0=tmp0+(U_MD(istep)-PijUi(istep)-PhijUiUj(istep))
      if (present(ftot3)) tmp0=tmp0-PsijUiUjUk(istep)
      tmp6=tmp6+U_MD(istep)
      tmp9=tmp9+PijUi(istep)
      tmp5=tmp5+PhijUiUj(istep)
      if (present(ftot3)) tmp7=tmp7+PsijUiUjUk(istep)
    end do  
    U_TDEP(:)=tmp0/real(istepmax)+PijUi(:)+PhijUiUj(:)
    if (present(ftot3)) U_TDEP(:) = U_TDEP(:)+PsijUiUjUk(:)
    do istep=1,istepmax
      tmp1 =tmp1 + (U_MD(istep)-U_TDEP(istep))
      tmp8 =tmp8 + (U_MD(istep)-U_TDEP(istep))**2
      tmp2 =tmp2 -((U_MD(istep)-U_TDEP(istep)) - tmp1/real(istepmax))**2
    end do
!   Compute eucledian distance for forces    
    do istep=istepmin,istepmax
      do iatom=1,InVar%natom
        do ii=1,3
          tmp3=tmp3+(Forces_MD(ii+3*(iatom-1)+3*InVar%natom*(istep-1))&
&                 -Forces_TDEP(ii+3*(iatom-1)+3*InVar%natom*(istep-1)))**2
          tmp4=tmp4+Forces_MD(ii+3*(iatom-1)+3*InVar%natom*(istep-1))**2
        end do
      end do  
    end do
    U0       =tmp0/real(istepmax*InVar%natom)
    UMD      =tmp6/real(istepmax*InVar%natom)
    U_1      =tmp9/real(istepmax*InVar%natom)
    U_2      =tmp5/real(istepmax*InVar%natom)
    U_3      =tmp7/real(istepmax*InVar%natom)
    Delta_U  =tmp1/real(istepmax*InVar%natom)
    Delta_U2 =tmp8/real(istepmax*InVar%natom)
    Delta_F2 =tmp3/real(istepmax*InVar%natom*3)
    sigma    =dsqrt(tmp3/tmp4)
    Free_Anh =tmp2/real(istepmax*InVar%natom)/(2.d0*kb_HaK*InVar%temperature)
    if (present(ftot3)) then
      write(InVar%stdout,'(i5,9(5x,f12.5))') istepmax,UMD*Ha_eV,U0*Ha_eV,U_1*Ha_eV,U_2*Ha_eV,U_3*Ha_eV,&
&       Delta_U*Ha_eV,Delta_U2**0.5*Ha_eV,Delta_F2**0.5,sigma
    else
      write(InVar%stdout,'(i5,8(5x,f12.5))') istepmax,UMD*Ha_eV,U0*Ha_eV,U_1*Ha_eV,U_2*Ha_eV,&
&       Delta_U*Ha_eV,Delta_U2**0.5*Ha_eV,Delta_F2**0.5,sigma
    endif 
  end do  
  write(InVar%stdout,'(a,1x,f12.5)') ' NOTE : in the harmonic and classical limit (T>>T_Debye), U_2=3/2*kB*T=',&
&   3.d0/2.d0*kb_HaK*Ha_eV*InVar%temperature

! Write : i) (U_TDEP vs U_MD) in etotMDvsTDEP.dat
!        ii) (Forces_TDEP vs Forces_MD) in fcartMDvsTDEP.dat
  write(InVar%stdout,'(a)') ' '
  if (present(ftot3)) then
    open(unit=32,file=trim(InVar%output_prefix)//'etotMDvsTDEP3.dat')
    write(32,'(a)') '#   Istep      U_MD(Ha)         U_TDEP(Ha)       PiiUi(Ha)       PhijUiUj/2(Ha)  PsijUiUjUk/6(Ha)'
    open(unit=33,file=trim(InVar%output_prefix)//'fcartMDvsTDEP3.dat')
    write(33,'(a)') '# Forces_MD(Ha/bohr) Forces_TDEP(Ha/bohr)'
  else  
    open(unit=32,file=trim(InVar%output_prefix)//'etotMDvsTDEP2.dat')
    write(32,'(a)') '#   Istep      U_MD(Ha)         U_TDEP(Ha)       PiiUi(Ha)       PhijUiUj/2(Ha)'
    open(unit=33,file=trim(InVar%output_prefix)//'fcartMDvsTDEP2.dat')
    write(33,'(a)') '# Forces_MD(Ha/bohr) Forces_TDEP(Ha/bohr)'
  end if  
  do istep=1,InVar%nstep
    if (present(ftot3)) then
      write(32,'(i6,1x,5(f17.6,1x))') istep,U_MD(istep),U_TDEP(istep),PijUi(istep),PhijUiUj(istep),PsijUiUjUk(istep)
    else  
      write(32,'(i6,1x,4(f17.6,1x))') istep,U_MD(istep),U_TDEP(istep),PijUi(istep),PhijUiUj(istep)
    end if

    do iatom=1,InVar%natom
      do ii=1,3
        write(33,'(2(f17.10,1x))') Forces_MD  (ii+3*(iatom-1)+3*InVar%natom*(istep-1)),&
&                                  Forces_TDEP(ii+3*(iatom-1)+3*InVar%natom*(istep-1))
      end do
    end do  
  end do  
  close(32)
  close(33)
  ABI_FREE(U_MD)
  ABI_FREE(U_TDEP)
  ABI_FREE(PijUi)
  ABI_FREE(PhijUiUj)
  if (present(ftot3)) then
    ABI_FREE(PsijUiUjUk)
  end if

! Compute average on forces on each atom and each direction
  ABI_MALLOC(residualF,(3*InVar%natom)); residualF(:)=zero
  ABI_MALLOC(Fmean    ,(3*InVar%natom)); Fmean(:)    =zero
  if (present(ftot3)) then
    write(InVar%stdout,'(a)') ' See the etotMDvsTDEP3.dat, fcartMDvsTDEP3.dat and residualF3.dat files'
    open(unit=34,file=trim(InVar%output_prefix)//'residualF3.dat')
  else  
    write(InVar%stdout,'(a)') ' See the etotMDvsTDEP2.dat, fcartMDvsTDEP2.dat, residualF1.dat and residualF2.dat files'
    open(unit=33,file=trim(InVar%output_prefix)//'residualF1.dat')
    do istep=1,InVar%nstep
      do iatom=1,InVar%natom
        do ii=1,3
          residualF(ii+3*(iatom-1))=residualF(ii+3*(iatom-1))+    Forces_MD(ii+3*(iatom-1)+3*InVar%natom*(istep-1))
          Fmean    (ii+3*(iatom-1))=Fmean    (ii+3*(iatom-1))+abs(Forces_MD(ii+3*(iatom-1)+3*InVar%natom*(istep-1)))
        end do
      end do
    end do
    Fmean    (:)=Fmean    (:)/real(InVar%nstep)
    residualF(:)=residualF(:)/real(InVar%nstep)
    do iatom=1,InVar%natom
      write(33,'(i6,1x,6(f17.10,1x))') iatom,(residualF(ii+3*(iatom-1)),ii=1,3),(Fmean(jj+3*(iatom-1)),jj=1,3)
    end do  
    close(33)
    open(unit=34,file=trim(InVar%output_prefix)//'residualF2.dat')
  end if  

! Compute residual on forces, Chi^2 and Sigma
  residualF(:)=0.d0
  do istep=1,InVar%nstep
    do iatom=1,InVar%natom
      do ii=1,3
        residualF(ii+3*(iatom-1))=residualF(ii+3*(iatom-1))+&
&          abs(Forces_MD  (ii+3*(iatom-1)+3*InVar%natom*(istep-1))&
&          -Forces_TDEP(ii+3*(iatom-1)+3*InVar%natom*(istep-1)))
      end do  
    end do  
  end do  
  residualF(:)=residualF(:)/real(InVar%nstep)
  do iatom=1,InVar%natom
    write(34,'(i6,1x,3(f17.10,1x))') iatom,(residualF(ii+3*(iatom-1)),ii=1,3)
  end do  
  close(34)
  ABI_FREE(residualF)
  ABI_FREE(Fmean)
  write(InVar%stdout,'(a)') ' NOTE: Large differences on residual forces can be observed between average and '
  write(InVar%stdout,'(a)') '          ideal positions (check with the "Use_Ideal_Positions" input parameter).'

 end subroutine tdep_calc_model

!====================================================================================================
subroutine tdep_calc_nbcoeff(distance,iatcell,InVar,ishell,jatom,katom,ncoeff,norder,nshell,order,proj,Sym)

  implicit none


  integer,intent(in) :: iatcell,ishell,jatom,nshell,katom,order,norder
  integer,intent(out) :: ncoeff
  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(out) :: proj(norder,norder,nshell)

  integer :: ii,jj,isym,LWORK,INFO,const_tot
  integer :: kk,ncount,icoeff,jatcell,katcell,mu,nu,xi
  integer :: inv,watom,xatom,yatom,isyminv,nsyminv,facorder
  integer, allocatable :: iconst(:)
  double precision :: prod_scal,drandom
  double precision :: eigvec(3,3)
  double precision :: vect_trial(3),vect_trial1(3),vect_trial2(3),vect_trial3(3)
  double precision :: WR(3),WI(3),VL(3,3),VR(3,3)
  double precision, allocatable :: WORK(:)
  double complex :: eigenvectors(3,3),eigenvalues(3)
  double complex :: pp(3,3),ppp(3,3,3),lambda
  double complex, allocatable :: tab_vec(:,:),temp(:,:),alphaij(:,:,:),constraints(:,:)
  logical, allocatable :: unchanged(:)

  if (iatcell==1.and.order==1) then
    ncoeff=0
    return
  end if    
  if (jatom==iatcell.and.order==2) then
    ncoeff=0
    return
  end if    
  if (katom==iatcell.and.jatom==iatcell.and.order==3) then
    ncoeff=0
    return
  end if    

  if (order==1) then
    facorder=1
  else if (order==2) then
    facorder=2
  else if (order==3) then
    facorder=6
  end if  

! If we want to remove the constraints coming from the symetries
!FB  if (order==3) then
!FB    do ii=1,norder
!FB      proj(ii,ii,ishell)=1.d0
!FB    end do
!FB    ncoeff=norder
!FB    return
!FB  end if  

  const_tot=0
  nsyminv=Sym%nsym*facorder
  ABI_MALLOC(alphaij,(nsyminv,norder,norder)); alphaij(:,:,:)=czero
  ABI_MALLOC(iconst,(nsyminv))               ; iconst(:)=0
  ABI_MALLOC(unchanged,(nsyminv))            ; unchanged(:)=.false.

! ================================================================================================
! ================ Big loop over symmetries and invariance (nsym*facorder) =======================
! ================================================================================================
  do isyminv=1,nsyminv
    isym=(isyminv-1)/facorder+1
    inv=isyminv-(isym-1)*facorder
    if (isym==1) cycle 

!   For the 1st order: Search if the atom is let invariant
    if (order==1) then
      if (Sym%indsym(4,isym,iatcell)==iatcell) then
        write(16,*) ' '
        write(16,*) 'For shell number=',ishell
        write(16,*)'===========The atom is let invariant for isym=',isym
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
        write(16,*) ' '
        write(16,*) 'For shell number=',ishell
        if (inv==1) write(16,*)'===========The bond is kept invariant for isym=',isym
        if (inv==2) write(16,*)'===========The bond is reversed with (i,j) --> (j,i) for isym=',isym
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
        write(16,*) ' '
        write(16,*) 'For shell number=',ishell
        if (inv==1) write(16,*)'===========The bond is kept invariant for isym=',isym
        if (inv==2) write(16,*)'===========The bond is reversed with (i,j,k) --> (i,k,j) for isym=',isym
        if (inv==3) write(16,*)'===========The bond is reversed with (i,j,k) --> (j,i,k) for isym=',isym
        if (inv==4) write(16,*)'===========The bond is reversed with (i,j,k) --> (j,k,i) for isym=',isym
        if (inv==5) write(16,*)'===========The bond is reversed with (i,j,k) --> (k,i,j) for isym=',isym
        if (inv==6) write(16,*)'===========The bond is reversed with (i,j,k) --> (k,j,i) for isym=',isym
      else
        cycle
      end if  
    end if  

!   Write the S_ref matrix
    write(16,'(3(f16.12,1x))') Sym%S_ref(1,1,isym,1),Sym%S_ref(1,2,isym,1),Sym%S_ref(1,3,isym,1)
    write(16,'(3(f16.12,1x))') Sym%S_ref(2,1,isym,1),Sym%S_ref(2,2,isym,1),Sym%S_ref(2,3,isym,1)
    write(16,'(3(f16.12,1x))') Sym%S_ref(3,1,isym,1),Sym%S_ref(3,2,isym,1),Sym%S_ref(3,3,isym,1)

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
    do ii=1,3
      write(16,*)'  For eigenvalue number',ii
      write(16,*)'     The eigenvalue is:',eigenvalues(ii)
      write(16,*)'     The eigenvector is:'
      write(16,'(2(f16.12,1x))') eigenvectors(1,ii)
      write(16,'(2(f16.12,1x))') eigenvectors(2,ii)
      write(16,'(2(f16.12,1x))') eigenvectors(3,ii)
      if ((aimag(eigenvalues(1)).ne.0).or.(aimag(eigenvalues(2)).ne.0).or.(aimag(eigenvalues(3)).ne.0)) then
        write(16,*) '  WARNING: THERE IS COMPLEX EIGENVALUES:'
      end if  
    end do  
    write(16,*) ' '

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
        const_tot=const_tot+1
        write(16,*)'  The eigenvalue',ii
        write(16,*)'  is equal to ',lambda  
        do mu=1,3
          alphaij(isyminv,mu,iconst(isyminv))=eigenvectors(mu,ii)
        end do  
        write(16,*)'  Real & imaginary parts of the eigenvectors product:'
        write(16,'(9(f16.12,1x))')  real(alphaij(isyminv,:,iconst(isyminv)))
        write(16,'(9(f16.12,1x))') aimag(alphaij(isyminv,:,iconst(isyminv)))
      end do !ii
    else if (order==2) then
      do ii=1,3
        do jj=1,3
          lambda=eigenvalues(ii)*eigenvalues(jj)
          if (((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6).and.(inv==1)).or.&
&             ((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6).and.(inv==2).and.(ii==jj))) cycle
          unchanged(isyminv)=.true.
          iconst(isyminv)=iconst(isyminv)+1
          const_tot=const_tot+1
          write(16,*)'  The product of eigenvalues',ii,jj
          write(16,*)'  is equal to ',lambda  
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
          write(16,*)'  Real & imaginary parts of the eigenvectors product:'
          write(16,'(9(f16.12,1x))')  real(alphaij(isyminv,:,iconst(isyminv)))
          write(16,'(9(f16.12,1x))') aimag(alphaij(isyminv,:,iconst(isyminv)))
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
            const_tot=const_tot+1
            write(16,*)'  The product of eigenvalues',ii,jj
            write(16,*)'  is equal to ',lambda  
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
            write(16,*)'  Real & imaginary parts of the eigenvectors product:'
            write(16,'(9(f16.12,1x))')  real(alphaij(isyminv,:,iconst(isyminv)))
            write(16,'(9(f16.12,1x))') aimag(alphaij(isyminv,:,iconst(isyminv)))
          end do !kk
        end do !jj
      end do !ii
    else
      ABI_BUG('Only the first, second and third order is allowed')
    end if  

!   WARNING: There are some minimum and maximum of constraints
    if (order==1.and.(iconst(isyminv).eq.3)) then
      ncoeff=0
      proj(:,:,ishell)=zero   
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
  end do !isyminv
! ================================================================================================
! =========== End big loop over symetries and facorder ===========================================
! ================================================================================================

! If the third order IFCs are symmetric with respect to the atoms (iik, iji, ijj), 
! some constraints have to be added :
  if (order.eq.3) then
    if ((iatcell.eq.jatom).or.(iatcell.eq.katom).or.(jatom.eq.katom)) then
    write(16,*) ' '
    write(16,*) 'For shell number=',ishell
    write(16,*)'=========== The IFCs are symmetric'
      const_tot=const_tot+norder
      ABI_MALLOC(constraints,(norder,norder)) ; constraints(:,:)=czero
      ii=0
      do mu=1,3
        do nu=1,3
          do xi=1,3
            ii=ii+1
            if (iatcell.eq.jatom) then 
              if (mu.eq.nu) cycle
              constraints((mu-1)*9+(nu-1)*3+xi,ii)= cone
              constraints((nu-1)*9+(mu-1)*3+xi,ii)=-cone
            end if 
            if (iatcell.eq.katom) then 
              if (mu.eq.xi) cycle
              constraints((mu-1)*9+(nu-1)*3+xi,ii)= cone
              constraints((xi-1)*9+(nu-1)*3+mu,ii)=-cone
            end if 
            if (jatom.eq.katom) then 
              if (nu.eq.xi) cycle
              constraints((mu-1)*9+(nu-1)*3+xi,ii)= cone
              constraints((mu-1)*9+(xi-1)*3+nu,ii)=-cone
            end if 
          end do  
        end do  
      end do  
    end if
  end if

! In the case where the matrix has norder elements  
  if (const_tot==0) then 
    write(16,*) ' '
    write(16,*) 'For shell number=',ishell
    write(16,*)'  WARNING: There is no symetry operation leaving the bond unchanged'
    write(16,*) isym,ishell,jatom,iatcell
    proj(:,:,ishell)=0.d0
    do ii=1,norder
      proj(ii,ii,ishell)=1.d0
    end do
    ncoeff=norder
    return
  end if

! When some constraints have been found
  ncount=const_tot
  write(16,*) ' '
  write(16,*) 'There is a total of ',ncount,' non-independant constraints for this shell'
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
  if (order.eq.3) then
    if ((iatcell.eq.jatom).or.(iatcell.eq.katom).or.(jatom.eq.katom)) then
      do jj=1,norder
        ii=ii+1
        tab_vec(:,ii)=constraints(:,jj)
      end do
      ABI_FREE(constraints)
    end if
  end if  
  if (ii.ne.ncount) then
    write(16,*) ii,' non equal to ',ncount
    ABI_BUG('Count error')
  end if  
  do ii=1,norder
    do jj=1,ncount
      if (abs( real(tab_vec(ii,jj))).lt.1.d-8) tab_vec(ii,jj)=dcmplx(zero,aimag(tab_vec(ii,jj)))
      if (abs(aimag(tab_vec(ii,jj))).lt.1.d-8) tab_vec(ii,jj)=dcmplx( real(tab_vec(ii,jj)),zero)
    end do
  end do

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
  write(16,*) ' '
  write(16,*) '  ========The final set of vectors is:'
  do kk=1,ncount
    write(16,'(9(f16.12,1x))')  real(temp(:,kk))
    write(16,'(9(f16.12,1x))') aimag(temp(:,kk))
  end do  
  write(16,*) '  =======Au total, il y a ',ncount,' vecteurs independants'
  if (ncount.gt.8.and.order==2) then
    ABI_ERROR(' Order 2 : There are too many independant vectors')
  end if
  if (ncount.gt.26.and.order==3) then
    ABI_ERROR(' Order 3 : There are too many independant vectors')
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
  write(16,*) ' '
  write(16,*) '  ========The orthogonal set of vectors is:'
  do kk=ncount+1,norder
    write(16,'(9(f16.12,1x))')  real(tab_vec(:,kk))
    write(16,'(9(f16.12,1x))') aimag(tab_vec(:,kk))
    if ((abs(aimag(tab_vec(1,kk))).gt.1.d-6).or.&
&       (abs(aimag(tab_vec(1,kk))).gt.1.d-6).or.&
&       (abs(aimag(tab_vec(1,kk))).gt.1.d-6)) then
      ABI_ERROR('the constraint has an imaginary part')
    end if
  end do  
  ncoeff=norder-ncount
  write(16,*) '  =======Au total, il y a ',ncoeff,' coefficients'

! On copie tab_vec dans proj
  do icoeff=1,ncoeff
    proj(:,icoeff,ishell)=tab_vec(:,ncount+icoeff)
  end do
  ABI_FREE(tab_vec)

end subroutine tdep_calc_nbcoeff
!=====================================================================================================

end module m_tdep_utils
