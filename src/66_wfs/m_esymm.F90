!!****m* ABINIT/m_esymm
!! NAME
!! m_esymm
!!
!! FUNCTION
!! This module defines structures and provides procedures used to find
!! the irreducible representations associated to electronic eigenstates.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_esymm

 use defs_basis
 use m_abicore
 use m_errors

 use m_io_tools,       only : file_exists
 use m_symtk,          only : matr3inv, chkgrp, symrelrot, littlegroup_q
 use m_symfind,        only : symbrav
 use m_fstrings,       only : int2char10, itoa, sjoin
 use m_numeric_tools,  only : print_arr, set2unit, get_trace
 use m_hide_lapack,    only : xgeev, xginv
 use m_crystal,        only : crystal_t
 use m_defs_ptgroups,  only : point_group_t, irrep_t
 use m_ptgroups,       only : get_classes, point_group_init, irrep_free,&
&                             copy_irrep, init_irrep, mult_table, sum_irreps

 implicit none

 private
!!***

 ! Error codes.
 integer,private,parameter :: ESYM_NOERROR             = 0
 integer,private,parameter :: ESYM_ACCDEG_ERROR        = 10
 integer,private,parameter :: ESYM_CLASSIFICATION_ERROR= 11
 integer,private,parameter :: ESYM_ORTHO_ERROR         = 12
 integer,private,parameter :: ESYM_UNITARY_ERROR       = 13
 integer,private,parameter :: ESYM_PTG_WRONG_MAPPING   = 20
 integer,private,parameter :: ESYM_HERRING_WRONG_TEST  = 30
 integer,private,parameter :: ESYM_HEUR_WRONG_NCLASSES = 40
 integer,private,parameter :: ESYM_HEUR_WRONG_DIMS     = 41

!----------------------------------------------------------------------

!!****t* m_esymm/esymm_t
!! NAME
!!  esymm_t
!!
!! FUNCTION
!!  Dataype gathering data and tables needed to analize the symmetries
!!  of electronic states at a given k-point via Group Theory.
!!
!! SOURCE

 type,public :: esymm_t

  integer :: nspinor
  ! Number of spinorial components.

  integer :: first_ib
  ! Index of the first treated band.

  integer :: nbnds
  ! Number of bands for this k-point and spin.

  integer :: nclass
  ! The number of classes in the group of k.

  integer :: nsym_gk
  ! Number of symmetries in the group of k. Namely that the set of symmetries such that Sk = k +G0.

  integer :: nsym_trgk
  ! Number of symmetries in the extended group of k. Namely that the set of symmetries such that -Sk = k + G0.

  integer :: err_status = ESYM_NOERROR
  ! Flag signaling if the classification algorithm succeed or not.

  real(dp) :: tol_deg
  ! Energy tolerance below which two states are considered degenerate.

  logical :: can_use_tr
  ! .TRUE. if time-reversal can be used

  logical :: only_trace
  ! if .TRUE. only the trace of a single matrix per class is calculated
  ! this is the standard way used to analyze bands symmetries. If .FALSE.
  ! the full matrices of the irreducible representations are calculated and stored

  logical :: has_spatial_inv
  ! .TRUE. if the inversion belongs to the space group

  logical :: nonsymmorphic_at_zoneborder
  ! if .TRUE. analysis cannot be performed since kpt is
  ! at border zone and non-zero fractional translations are present in the space group

  logical :: has_chtabs
  ! True if Ref_irreps and character tables are available (tables are initialized either
  ! from point group irreps or from an external database downloaded from the Bilbao server)

  real(dp) :: kpt(3)
  ! The crystalline momentum of the wavefunctions in reduced coordinates.

  character(len=500) :: err_msg="None"

  integer,allocatable :: g0(:,:)
  ! g0(3,nsym_gk)
  ! The umklapp g0 vector associated to each little group operation.

  integer,allocatable :: tr_g0(:,:)
  ! tr_g0(3,nsym_trgk)
  ! The umklapp g0 vector associated to each little group operation.

  integer :: ndegs
  ! Number of degenerate states.

  integer,allocatable :: nelements(:)
  ! nelements(nclass)
  ! Number of symmetry operations in each class.

  integer,allocatable :: sgk2symrec(:)
  ! sgk2symrec(nsym_gk)
  ! Mapping between the symmetries of the group of k and the symrec(l) array.
  ! The symmetries of the little group are always packed in classes to facilitate
  ! the calculation of the character of the irrep. Abinit symmetries are randomly ordered.

  integer,allocatable :: tr_sgk2symrec(:)
  ! trsgk2symrec(nsym_trgk)
  ! Mapping between the symmetries of the group of k and the symrec(l) array.
  ! The symmetries of the little group are always packed in classes to facilitate
  ! the calculation of the character of the irrep. Abinit symmetries are randomly ordered.

  integer,allocatable :: herring_test(:)
  ! herring_test(nclass)
  ! The result of Herring test for each irreducible representantion of the group of k.
  ! Possible values are: +1, 0, -1

  integer,allocatable :: b2irrep(:)
  ! b2irrep(nbnds)
  ! For each band, it gives the index of the irreducible representation in Ref_irreps.

  type(coeffi1_type),allocatable :: irrep2b(:)
  ! irrep2b(0:nclass)%value(:)
  ! Ragged arrays with the mapping between the set of irreducible representation and the band indices.
  ! irrep2b(irp)%value(:) gives the indeces of the states belonging to irrep irp, irp=1,nclass
  ! irrep2b(0)%value(:) stores the indeces of the states that have not been classified due to
  !   the presence of an accidental degeneracy.

  integer,allocatable :: degs_bounds(:,:)
  ! degs_bounds(2,ndegs)
  !   degs_bounds(1,idg)= first band index of the degenerate set idg=1,ndegs
  !   degs_bounds(2,idg)= final band index of the degenerate set idg=1,ndegs

  integer,allocatable :: degs_dim(:)
  ! degs_dim(ndegs)
  ! Number of states in each degenerate subspace. Cannot be larger that nclass provided
  ! that no accidental degeneracy occurs.

  !% integer,allocatable :: class_ids(:,:)
  ! class_ids(2,nclass)
  ! (1,icl) = index of the first symmetry of class icl
  ! (2,icl) = index of the last symmetry of class icl
  ! Note that symmetries in sym are packed in classes.

  type(irrep_t),allocatable :: Calc_irreps(:)
  ! Calc_irreps(ndegs)
  !  The representations of the little group of k calculated from the wavefunctions. <\phi_nk|R_t|\phi_mk>
  !  where R_t belong to the little group of k.
  !  They represent an unitary irreducible representation provided that no accidental degeneracy occurs.

  type(irrep_t),allocatable :: trCalc_irreps(:)
  ! trCalc_irreps(ndegs)
  !  The representations of the little group of k calculated from the wavefunctions. <\phi_nk|R_t|\phi_mk>
  !  where R_t belong to the little group of k.
  !  They represent an unitary irreducible representation provided that no accidental degeneracy occurs.

  type(irrep_t),allocatable :: Ref_irreps(:)
  ! Irreps(nclass)
  !   Reference irreducible representations of the group of k derived from the point group
  !   or from the external database downloaded from the Bilbao web site.

 end type esymm_t

 public :: esymm_init             ! Initialize the object
 public :: esymm_print            ! Print info
 public :: esymm_free             ! Free memory
 public :: esymm_finalize         ! Finalize the object
 public :: esymm_symmetrize_mels  ! Symmetrize given matrix elements
 public :: esymm_failed           ! True if symmetry analysis failed.
!!***

!----------------------------------------------------------------------

 !public :: polish_irreps ! TODO method of Irreps_t, therefore should be moved to m_ptgroups.
                          ! but first one has to solve the dependency on m_abilasi and scalapack
 interface esymm_free
   module procedure esymm_free_0D
   module procedure esymm_free_2D
 end interface esymm_free

CONTAINS  !==============================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_esymm/esymm_init
!! NAME
!! esymm_init
!!
!! FUNCTION
!!  Initialize a esymm_t datatype containing data and parameters
!!  needed to analyze the irreducible representations at a particular k-point....
!!
!! INPUTS
!!  kpt_in(3)=The k-point where the classification of bands is required.
!!  Cryst<crystal_t>=Datatype describing the unit cell and its symmetries.
!!  nspinor=number of spinorial components
!!  nsppol=number of independent polarizations
!!  first_ib=Index of the first band.
!!  nbnds=Number of bands for this k-point.
!!  ene_k(nbnds)=energies for this k-point. ene_k(1) corresponds to band first_ib.
!!  EDIFF_TOL=tolerance below which two states are considered to belong to the same irreducible representation
!!
!! OUTPUT
!!  esymm<esymm_t>= Initialized data type gathering information of the small group
!!     of the k-point as well as the irreducible representations.
!!
!! NOTES
!!   The present implementation does NOT work at zone border if the little group of
!!   kpt_in is non-symmorphic namely thers is at lest a symmetry operation with non-zero tnons.
!!
!! PARENTS
!!      classify_bands
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine esymm_init(esymm,kpt_in,Cryst,only_trace,nspinor,first_ib,nbnds,EDIFF_TOL,ene_k,tolsym)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbnds,nspinor,first_ib
 real(dp),intent(in) :: EDIFF_TOL,tolsym
 logical,intent(in) :: only_trace
 type(crystal_t),intent(in) :: Cryst
 type(esymm_t),intent(out) :: esymm
!arrays
 real(dp),intent(in) :: ene_k(nbnds),kpt_in(3)

!Local variables-------------------------------
!scalars
 integer :: dim_degs,iband,idg,irp,nacc_deg,isym_gk,grp_ierr
 integer :: nsym_fm,idx_fm,idx_gk,idx_trgk,isym,jsym,dummy_timrev !,iholohedry
 integer :: iel,icls,msym,iord !isym1,!iprod,dim_irrep,icls2, isym2,isym_tr,
 integer :: spgroup,chkprim !,ptgroupma
 real(dp) :: mkt
 !complex(dpc) :: phase_k
 character(len=5) :: ptgroup,ptgroup_name
 character(len=10) :: spgroup_str
 character(len=1000) :: msg
 character(len=fnlen) :: lgroup_fname
!arrays
 integer :: inversion(3,3)
 integer,allocatable :: degs_bounds(:,:),dim_irreps(:)
 integer :: bravais(11),sym_axis(3)
 real(dp) :: pmat1(3,3),pmat2(3,3),pmat3(3,3),pmat4(3,3),pmat5(3,3),pmat6(3,3)
 !real(dp) :: genafm(3)
 !integer :: rot2(3,3)
 !integer,allocatable :: mtab(:,:)
 integer,allocatable :: elements_idx(:,:),tmp_nelements(:)
 integer,allocatable :: found(:),symrec_fm(:,:,:),fm2symrec(:)
 integer,allocatable :: ksym_table(:,:,:),sgk(:,:,:),tr_sgk(:,:,:),dum_symafm(:)
 integer,allocatable :: new_idx(:),new_g0(:,:),tmp_symrec(:,:,:),conv_symrec(:,:,:) !,tr_conv_symrec(:,:,:)
 integer,allocatable :: dummy_symafm(:)
 real(dp) :: conv_gprimd(3,3),axes(3,3) !,tau2(3)
 !complex(dpc),allocatable :: her_test(:) !,mat_test(:,:)
 complex(dpc),allocatable :: phase_mkt(:)
 type(point_group_t) :: Ptg

! *************************************************************************

 DBG_ENTER("COLL")

 !@esymm_t
 esymm%err_status= ESYM_NOERROR
 inversion=RESHAPE((/-1,0,0,0,-1,0,0,0,-1/),(/3,3/))
 !
 ! ====================================
 ! ==== Initialize basic variables ====
 ! ====================================
 esymm%nspinor        = nspinor
 esymm%first_ib       = first_ib
 esymm%nbnds          = nbnds
 esymm%only_trace     = only_trace
 esymm%tol_deg        = EDIFF_TOL
 esymm%has_spatial_inv= (cryst%idx_spatial_inversion() /= 0)
 esymm%can_use_tr     = .TRUE. !TODO this should be input
 esymm%has_chtabs     = .FALSE.
 esymm%kpt            = kpt_in(:)
 esymm%nonsymmorphic_at_zoneborder=.FALSE.
 !
 ! ===============================
 ! === Locate degenerate_bands ===
 ! ===============================
 esymm%ndegs=1

 ABI_MALLOC(degs_bounds,(2,nbnds))
 degs_bounds=0; degs_bounds(1,1)=1

 do iband=2,nbnds
   if (ABS(ene_k(iband)-ene_k(iband-1))>EDIFF_TOL) then
     degs_bounds(2,esymm%ndegs) = iband-1 + (first_ib-1)
     esymm%ndegs=esymm%ndegs+1
     degs_bounds(1,esymm%ndegs) = iband + (first_ib-1)
   end if
 end do
 degs_bounds(2,esymm%ndegs)=nbnds + (first_ib-1)

 ABI_MALLOC(esymm%degs_bounds,(2,esymm%ndegs))
 esymm%degs_bounds = degs_bounds(:,1:esymm%ndegs)
 ABI_FREE(degs_bounds)
 !
 ! Each band is initialized as "Unknown".
 ABI_MALLOC(esymm%b2irrep,(esymm%nbnds))
 esymm%b2irrep = 0
 !
 ! ==================================
 ! ==== Find the group of kpt_in ====
 ! ==================================
 ! The small point group is the subset of symrec such that $ S q = q + g0 $
 ! Symmetries are packed in classes.
 ! For the time being, AFM symmetries are not treated.

 write(msg,'(a,3(1x,f7.4))')" Finding the little group of k-point: ",esymm%kpt
 call wrtout(std_out,msg,"COLL")

 ! Only FM symmetries are used.
 nsym_fm = COUNT(Cryst%symafm==1)

 if (nsym_fm /= Cryst%nsym) then
   write(msg,'(4a)')ch10,&
&    "Band classification in terms of magnetic space groups not coded! ",ch10,&
&    "Only the ferromagnetic subgroup will be used "
   MSG_COMMENT(msg)
 end if

 ABI_MALLOC(symrec_fm,(3,3,nsym_fm))
 ABI_MALLOC(fm2symrec,(nsym_fm))

 idx_fm = 0
 do isym=1,Cryst%nsym
   if (Cryst%symafm(isym) == 1) then
     idx_fm = idx_fm+1
     symrec_fm(:,:,idx_fm) = Cryst%symrec(:,:,isym)
     fm2symrec(idx_fm) = isym
   end if
 end do

 ! Find symmetries that preserve k.
 ABI_MALLOC(ksym_table,(4,2,nsym_fm))
 ABI_MALLOC(dummy_symafm,(nsym_fm))

 dummy_symafm = 1
 call littlegroup_q(nsym_fm,esymm%kpt,ksym_table,symrec_fm,dummy_symafm,dummy_timrev,prtvol=0)

 esymm%nsym_gk  =COUNT(ksym_table(4,1,:)==1)  ! # S such that  S k = k +G0

 esymm%nsym_trgk=0
 if (esymm%can_use_tr) esymm%nsym_trgk=COUNT(ksym_table(4,2,:)==1)  ! # S such that -S k = k +G0

 ! Allocate workspace arrays.
 ABI_MALLOC(sgk,(3,3,esymm%nsym_gk))
 ABI_MALLOC(tr_sgk,(3,3,esymm%nsym_trgk))

 ! Allocate mapping little-group --> symrec and table for umklapps.
 ABI_MALLOC(esymm%sgk2symrec,(esymm%nsym_gk))
 ABI_MALLOC(esymm%g0,(3,esymm%nsym_gk))
 ABI_MALLOC(esymm%tr_sgk2symrec,(esymm%nsym_trgk))
 ABI_MALLOC(esymm%tr_g0,(3,esymm%nsym_trgk))

 ! Important NOTE:
 ! If nonsymmorphic_at_zoneborder symmetry analysis cannot be performed unless
 ! an external database retrieved from the bilbao server (REPRES) is found.
 idx_gk=0; idx_trgk=0
 esymm%sgk2symrec=-999; esymm%tr_sgk2symrec=-999
 do isym=1,nsym_fm

   if (ksym_table(4,1,isym)==1) then ! S k = k +G0
     idx_gk=idx_gk+1
     sgk(:,:,idx_gk)=symrec_fm(:,:,isym)
     esymm%g0(:,idx_gk)=ksym_table(1:3,1,isym)
     esymm%sgk2symrec(idx_gk)=fm2symrec(isym)
     if (ANY(ksym_table(1:3,1,isym)/=0).and.(ANY(ABS(Cryst%tnons(:,fm2symrec(isym)))>tol6))) then
        esymm%nonsymmorphic_at_zoneborder=.TRUE.
     end if
   end if

   if (esymm%can_use_tr.and.ksym_table(4,2,isym)==1) then ! -S k = k +G0
     idx_trgk=idx_trgk+1
     tr_sgk(:,:,idx_trgk)=symrec_fm(:,:,isym)
     esymm%tr_g0(:,idx_trgk)=ksym_table(1:3,2,isym)
     esymm%tr_sgk2symrec(idx_trgk)=fm2symrec(isym)
   end if
 end do

 ABI_FREE(ksym_table)
 ABI_FREE(symrec_fm)
 ABI_FREE(fm2symrec)

! ==========================================
! ==== Divide the operations in classes ====
! ==========================================
 ABI_MALLOC(dum_symafm,(esymm%nsym_gk))
 dum_symafm=1

 call chkgrp(esymm%nsym_gk,dum_symafm,sgk,grp_ierr)

 ABI_CHECK(grp_ierr==0,"chkgrp failed")
 ABI_FREE(dum_symafm)

 ABI_MALLOC(tmp_nelements,(esymm%nsym_gk))
 ABI_MALLOC(elements_idx,(esymm%nsym_gk,esymm%nsym_gk))

 call get_classes(esymm%nsym_gk,sgk,esymm%nclass,tmp_nelements,elements_idx)

 ABI_MALLOC(esymm%nelements,(esymm%nclass))
 esymm%nelements = tmp_nelements(1:esymm%nclass)
 ABI_FREE(tmp_nelements)

 ! From the list of symmetry operations and the lattice vectors, determine the
 ! Bravais information including the holohedry, the centering, the coordinate of
 ! the primitive vectors in the conventional vectors, as well as the point group,
 msym=192; if (allocated(Cryst%symrec)) msym=size(Cryst%symrec,3)
 ABI_MALLOC(tmp_symrec,(3,3,msym))
 tmp_symrec(:,:,1:esymm%nsym_gk)=sgk

 call symbrav(bravais,msym,esymm%nsym_gk,ptgroup,Cryst%gprimd,tmp_symrec,tolsym,axis=sym_axis)

 ABI_FREE(tmp_symrec)

 write(std_out,'(a)')" symptgroup returned point group: "//TRIM(ptgroup)
 write(std_out,'(a,i2)')" iholohedry ",bravais(1)
 write(std_out,'(a,i2)')" center     ",bravais(2)
 write(std_out,'(a,9i3)')" gprimd in the axes of the conventional bravais lattice (*2 if center/=0)",bravais(3:11)
 write(std_out,'(a,3i3)')" sym_axis ",sym_axis

 ! Branching:
 ! 1) If the little group is not symmorphic_at_zoneborder we can
 !    classify the states using the irreducible representation of the point group.
 !
 ! 2) If the little group is symmorphic_at_zoneborder, we have to rely on
 !    an external database retrieved from the Bilbao server in order to classify the states.
 !    If the file is not available, we only know the number of classes but neither their
 !    character nor the dimension of the irreducible representation.
 !
 if (esymm%nonsymmorphic_at_zoneborder) then

   spgroup=0
   chkprim=1 ! Cell must be primitive.
   !call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
   !call symspgr(bravais,Cryst%nsym,spgroup,Cryst%symrel,Cryst%tnons,tolsym)

   !call symanal(bravais,chkprim,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tolsym)

   call int2char10(spgroup,spgroup_str)
   lgroup_fname = "lgroup_"//TRIM(spgroup_str)

   if (file_exists(lgroup_fname)) then
     MSG_ERROR("Not coded")

     ! Read little groups from the external database.
     !% call init_groupk_from_file(Lgrp,spgroup,lgroup_fname,ierr)

     ! Save the irreducible representations in esymm.
     ! Reorder symmetries such that they correspond to the Bilbao database.
     !% allocate(esymm%Ref_irreps(esymm%nclass))
     !% call copy_irrep(Irreps, esymm%Ref_irreps)

   else
     write(msg,'(7a)')&
&      "Non-symmorphic small group and zone border. ",ch10,&
&      "External file: ",TRIM(lgroup_fname)," containing Bilbao tables not found ",ch10,&
&      "Character analysis cannot be performed. Accidental degeneracies cannot be detected. "
     MSG_WARNING(msg)

     esymm%has_chtabs = .FALSE.

     ! Reorder indeces such that symmetries are packed in classes.
     ABI_MALLOC(new_idx,(esymm%nsym_gk))
     ABI_MALLOC(new_g0,(3,esymm%nsym_gk))
     new_g0=0; iord = 0
     do icls=1,esymm%nclass
       do iel=1,esymm%nelements(icls)
         iord = iord+1
         jsym = elements_idx(iel,icls)
         new_idx(iord)  = esymm%sgk2symrec(jsym)
         new_g0(:,iord) = esymm%g0(:,jsym)
       end do
     end do

     esymm%sgk2symrec = new_idx
     esymm%g0 = new_g0

     ABI_FREE(new_idx)
     ABI_FREE(new_g0)
   end if ! file exists

 else
   !
   ! **** This part is still under development. It might not work for particular ****
   ! **** orientations of the unit cell or particular lattices.                  ****
   !
   ! The symmetries in the Bilbao database refer to the conventional unit cells.
   ! Therefore we have to map the abinit symmetries (in reduced coordinates)
   ! onto the Bilbao dataset. Bilbao standard settings are:
   !
   ! * unique axis b (cell choice 1) for space groups withing the monoclinic system
   ! * obverse triple hexagonal unit cell R space groups.
   ! * origin choice two - inversion center at (0, 0, 0) - for the centrosymmetric
   !   space groups for which there are two origins choices, within the
   !   orthorombic, tetragonal and cubic system.

   ! 1) Retrieve the rotation matrices and the irreducible representations (Bilbao setting).
   call point_group_init(Ptg,ptgroup)

   esymm%has_chtabs = .TRUE.
   ABI_CHECK(esymm%nclass==Ptg%nclass,"esymm%nclass/=Ptg%nclass!")

   do icls=1,esymm%nclass ! FIXME this is awful, should be done in a cleaner way.
     esymm%nelements(icls)=Ptg%class_ids(2,icls) - Ptg%class_ids(1,icls) + 1
   end do

   ! 2) Generate the symmetry operations in the conventional vector coordinates.
   conv_gprimd(:,1)=bravais(3:5)
   conv_gprimd(:,2)=bravais(6:8)
   conv_gprimd(:,3)=bravais(9:11)

   axes = conv_gprimd
   call matr3inv(conv_gprimd,axes) !; axes=TRANSPOSE(axes)

   conv_gprimd=MATMUL(Cryst%gprimd,TRANSPOSE(axes))
   !conv_gprimd=MATMUL(axes,Cryst%gprimd)
   !conv_gprimd=MATMUL(TRANSPOSE(axes),Cryst%gprimd)
   !write(std_out,*)"conv_gprimd:", conv_gprimd

   ptgroup_name = ADJUSTL(ptgroup)

   select case (ptgroup_name)

   case ("3m","-3m")
     call wrtout(std_out," Changing the conventional cell: rhombohedral --> triple hexagonal","COLL")
     ! Transformation matrices: primitive rhombohedral --> triple hexagonal cell obverse setting. Table 5.1.3.1 ITA page 81.
     pmat1 = RESHAPE( (/ 1,-1, 0, 0, 1,-1, 1, 1, 1/), (/3,3/) ) ! R1
     pmat2 = RESHAPE( (/ 0, 1,-1,-1, 0, 1, 1, 1, 1/), (/3,3/) ) ! R2
     pmat3 = RESHAPE( (/-1, 0, 1, 1,-1, 0, 1, 1, 1/), (/3,3/) ) ! R3
     pmat4 = RESHAPE( (/-1, 1, 0, 0,-1, 1, 1, 1, 1/), (/3,3/) ) ! R1 reverse setting.
     pmat5 = RESHAPE( (/ 0,-1, 1, 1, 0,-1, 1, 1, 1/), (/3,3/) ) ! R2 reverse setting.
     pmat6 = RESHAPE( (/ 1, 0,-1,-1, 1, 0, 1, 1, 1/), (/3,3/) ) ! R3 reverse setting.
     conv_gprimd = MATMUL(conv_gprimd,pmat1)
     !conv_gprimd = MATMUL(conv_gprimd,pmat2)
     !conv_gprimd = MATMUL(conv_gprimd,pmat3)
     !conv_gprimd = MATMUL(conv_gprimd,pmat4)
     !conv_gprimd = MATMUL(conv_gprimd,pmat5)
     !conv_gprimd = MATMUL(conv_gprimd,pmat6)
     !write(std_out,*)" New conv_gprimd:", conv_gprimd

   case ("mm2")
     call wrtout(std_out," Changing the conventional cell: unconventional orthorhombic setting --> conventional","COLL")
     ! Transformation matrices: unconvential orthorhombic --> conventional orthorhombic. Table 5.1.3.1 ITA page 81.
     pmat1 = RESHAPE( (/ 0, 1, 0, 1, 0, 0, 0, 0,-1/), (/3,3/) )  ! ( b, a,-c) --> (a,b,c)
     pmat2 = RESHAPE( (/ 0, 1, 0, 0, 0, 1, 1, 0, 0/), (/3,3/) )  ! ( c, a, b) --> (a,b,c)
     pmat3 = RESHAPE( (/ 0, 0, 1, 0, 1, 0,-1, 0, 0/), (/3,3/) )  ! (-c, b, a) --> (a,b,c)
     pmat4 = RESHAPE( (/ 0, 0, 1, 1, 0, 0, 0, 1, 0/), (/3,3/) )  ! ( b, c, a) --> (a,b,c)
     pmat5 = RESHAPE( (/ 1, 0, 0, 0, 0, 1, 0,-1, 0/), (/3,3/) )  ! ( a,-c, b) --> (a,b,c)
     conv_gprimd = MATMUL(conv_gprimd,pmat2)
     !write(std_out,*)" New conv_gprimd:", conv_gprimd
   case default
     continue
   end select

   ABI_MALLOC(conv_symrec,(3,3,esymm%nsym_gk))
   conv_symrec = sgk

   !axes=zero; axes(1,1)=one ; axes(2,2)=one ; axes(3,3)=one
   !call symrelrot(esymm%nsym_gk,conv_gprimd,axes,conv_symrec,tolsym)
   call symrelrot(esymm%nsym_gk,Cryst%gprimd,conv_gprimd,conv_symrec,tolsym)

   ! 3) Reorder indeces such that symmetries are packed in classes.
   ABI_MALLOC(found,(esymm%nsym_gk))
   ABI_MALLOC(new_idx,(esymm%nsym_gk))
   ABI_MALLOC(new_g0,(3,esymm%nsym_gk))
   new_g0=0; found=0

   do isym=1,esymm%nsym_gk
     do jsym=1,esymm%nsym_gk
       if (ALL(Ptg%sym(:,:,isym) == conv_symrec(:,:,jsym) ))  then
         found(isym)    = found(isym) + 1
         new_idx(isym)  = esymm%sgk2symrec(jsym)
         new_g0(:,isym) = esymm%g0(:,jsym)
         !EXIT
       end if
     end do
   end do
   !
   ! DEBUGGING SECTION
   !do isym=1,esymm%nsym_gk
   !  jsym=esymm%sgk2symrec(isym)
   !  call print_symmetries(1,Cryst%symrec(:,:,jsym),Cryst%tnons(:,jsym),Cryst%symafm(jsym))
   !  write(std_out,*)esymm%g0(:,isym)
   !end do

   if ( Ptg%nsym/=esymm%nsym_gk .or. ANY(found/=1) ) then
     !write(std_out,*)Ptg%nsym, esymm%nsym_gk
     !write(std_out,'(a,(i2))')" found = ",found
     write(std_out,*)" Ptg%sym list, conv_symrec list,  found Ptg% "
     do isym=1,Ptg%nsym
       write(std_out,'(a,i2,a,9i2,4x,a,9i2)')" found ",found(isym)," Ptg ",Ptg%sym(:,:,isym),"conv_symrec ",conv_symrec(:,:,isym)
     end do
     msg = " sgk and esymm%Ptg are inconsistent. Check tables or source"
     MSG_WARNING(msg)
     esymm%err_msg = msg(1:500)
     esymm%err_status = ESYM_PTG_WRONG_MAPPING
     esymm%has_chtabs = .FALSE.

   else ! Reorder symmetries.
     esymm%sgk2symrec = new_idx
     esymm%g0 = new_g0
   end if

   ABI_FREE(new_idx)
   ABI_FREE(new_g0)
   ABI_FREE(found)
   ABI_FREE(conv_symrec)

   if (esymm%has_chtabs) then
     ! Multiply the point group irreps by e^{-ik.\tau} to have the irreps of the little group.
     ! Store the results in esymm%Ref_irreps so that one can classify the states afterwards.
     ABI_DT_MALLOC(esymm%Ref_irreps,(esymm%nclass))
     ABI_MALLOC(phase_mkt,(esymm%nsym_gk))

     do isym_gk=1,esymm%nsym_gk
       isym =  esymm%sgk2symrec(isym_gk)
       mkt = -two_pi * DOT_PRODUCT(esymm%kpt, Cryst%tnons(:,isym))
       phase_mkt(isym_gk) = CMPLX(DCOS(mkt), DSIN(mkt))
     end do

     call copy_irrep(Ptg%Irreps,esymm%Ref_irreps,phase_mkt)
     ABI_FREE(phase_mkt)
   end if

#if 0
   ! Herring test requires the evaluation of the expression:
   !
   !   sum_{S,\tau} \chi^{k,\alpha} ({S|\tau}^2)
   !
   ! where Sk = -k + g0, and \chi is the trace of the \alpha-th
   ! irreducible representation of the little group of k.
   ! \chi^{k,\alpha} = e^{-ik.\tau} \chi(\alpha) provided that
   ! we are not at zone border with a non-symmorphic operation.
   ! The expression is always real and it can only be equal to \pm Ptg%nsym or zero.
   ! FIXME this part has to be rewritten from scratch.
   !if (esymm%err_status/=esymm_NOERROR) then
   !  write(std_out,*)" Skipping Herring test"
   !  goto 110
   !end if

   if (esymm%can_use_tr) then
     ABI_MALLOC(her_test,(esymm%nclass))

     ABI_MALLOC(tr_conv_symrec,(3,3,esymm%nsym_trgk))
     do isym_tr=1,esymm%nsym_trgk
       isym = esymm%tr_sgk2symrec(isym_tr)
       tr_conv_symrec(:,:,isym_tr)=Cryst%symrec(:,:,isym)
     end do

     call symrelrot(esymm%nsym_trgk,Cryst%gprimd,conv_gprimd,tr_conv_symrec_tr,tolsym)

     do isym_tr=1,esymm%nsym_trgk
       isym = esymm%tr_sgk2symrec(isym_tr)
       !rot2 = MATMUL(tr_sgk(:,:,isym),tr_sgk(:,:,isym))
       !tau2 = MATMUL(tr_sgk(:,:,isym),Cryst%tnons(:,isym)) + Cryst%tnons(:,isym)

       rot2 = MATMUL(tr_conv_symrec(:,:,isym_tr),tr_conv_symrec(:,:,isym_tr))
       tau2 = MATMUL(tr_conv_symrec(:,:,isym_tr),Cryst%tnons(:,isym)) + Cryst%tnons(:,isym)

       phase_k = EXP(-j_dpc*two_pi*DOT_PRODUCT(kpoint,tau2))
       call locate_sym(Ptg,rot2,isym2,icls2)

       do irp=1,esymm%nclass
         her_test(irp) = her_test(irp) + phase_k * Ptg%Irreps(irp)%trace(icls2)
       end do
     end do

     ABI_FREE(tr_conv_symrec)

     ! FIXME
     ABI_MALLOC(esymm%herring_test,(esymm%nclass))

     do irp=1,esymm%nclass
       if ( ABS(her_test(irp) - Ptg%nsym) < tol6 ) then
         esymm%herring_test(irp) = +1
       else if ( ABS(her_test(irp)) < tol6 ) then
         esymm%herring_test(irp) =  0
       else if ( ABS(her_test(irp) + Ptg%nsym) < tol6 ) then
         esymm%herring_test(irp) = -1
       else
         write(msg,'(a,i2,2a,i0,a,i2)')&
&          "Herring test for the irreducible representation number ",irp,ch10,&
&          "gave ",esymm%herring_test(irp),", while it should be 0 or +- ",Ptg%nsym
          MSG_WARNING(msg)
          esymm%err_msg   =msg
          esymm%err_status=esymm_HERRING_WRONG_TEST
       end if
     end do

     ABI_FREE(her_test)
   end if ! can_use_tr
#endif
   !
   ! Final check
   !allocate(mtab(esymm%nsym_gk,esymm%nsym_gk))
   !call mult_table(esymm%nsym_gk,Ptg%sym,mtab)

   !do isym=1,esymm%nsym_gk
   !  isym1 = esymm%sgk2symrec(isym)
   !  do jsym=1,esymm%nsym_gk
   !    isym2 = esymm%sgk2symrec(jsym)
   !    rot2 = MATMUL(Cryst%symrec(:,:,isym1),Cryst%symrec(:,:,isym2))

   !    iprod = mtab(isym,jsym)

   !    do irp=1,esymm%nclass
   !       dim_irrep = Ptg%Irreps(irp)%dim
   !       allocate(mat_test(dim_irrep,dim_irrep))
   !       mat_test = Ptg%Irreps(irp)%mat(:,:,isym) * Ptg%Irreps(irp)%mat(:,:,jsym)
   !       !call locate_sym(Ptg,rot2,isym2,icls2)
   !       write(std_out,*)mat_test - Ptg%Irreps(irp)%mat(:,:,iprod)
   !       deallocate(mat_test)
   !    end do
   !
   !  end do
   !end do
   !
   !deallocate(mtab)
 end if

 ABI_FREE(sgk)
 ABI_FREE(tr_sgk)
 ABI_FREE(elements_idx)

 !% allocate(esymm%irrep2b(0:esymm%nclass))
 !% call nullify_coeff(esymm%irrep2b)
 !
 ! 1) Allocate space for the irreducible representations.

 ! 2) Try to determine if we are in presence of an accidental degeneracy. Sufficient condition:
 !    There exists a set of degenerate states whose dimension is greater than the dimension
 !    of the irreducible representations of the point group. The check can be done only
 !    if Character tables are available.

 if (esymm%has_chtabs) then
   ABI_MALLOC(dim_irreps,(esymm%nclass))
   dim_irreps = (/(esymm%Ref_irreps(irp)%dim, irp=1,esymm%nclass)/)
 end if

 nacc_deg=0
 ABI_MALLOC(esymm%degs_dim,(esymm%ndegs))
 ABI_DT_MALLOC(esymm%Calc_irreps,(esymm%ndegs))

 if (esymm%can_use_tr)  then
   ABI_DT_MALLOC(esymm%trCalc_irreps,(esymm%ndegs))
 end if

 do idg=1,esymm%ndegs
   dim_degs=esymm%degs_bounds(2,idg)-esymm%degs_bounds(1,idg)+1

   if (esymm%has_chtabs) then
     if (ALL(dim_degs /= dim_irreps)) then ! An accidental degeneracy is present.
       nacc_deg=nacc_deg+1
     end if
   end if

   esymm%degs_dim(idg) = dim_degs

   call init_irrep(esymm%Calc_irreps(idg),esymm%nsym_gk,dim_degs)
   if (esymm%can_use_tr) then
     call init_irrep(esymm%trCalc_irreps(idg),esymm%nsym_trgk,dim_degs)
   end if
 end do ! idg

 if (esymm%has_chtabs) then
   ABI_FREE(dim_irreps)
   if (nacc_deg/=0) then
     write(msg,'(a,i0,a)')" Detected ",nacc_deg," accidental degeneracies."
     MSG_WARNING(msg)
     esymm%err_status=ESYM_ACCDEG_ERROR
     ! TODO this should signal to the caller that we have to decompose the calculated representation.
     esymm%err_msg =msg(1:500)
   end if
 end if

 DBG_EXIT("COLL")

end subroutine esymm_init
!!***

!----------------------------------------------------------------------

!!****f* m_esymm/esymm_print
!! NAME
!! esymm_print
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      classify_bands
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine esymm_print(esymm,unit,mode_paral,prtvol)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=4),optional,intent(in) :: mode_paral
 type(esymm_t),intent(in) :: esymm

!Local variables-------------------------------
!scalars
 integer :: icl,idg,my_unt,my_prtvol
 integer :: irr_idx,nstates,nunknown,istart,istop,ii
 character(len=4) :: my_mode
 character(len=1000) :: fmt,msg,msg0

! *********************************************************************

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 write(fmt,*)'(2a,3f8.4,3a,i4,2a,i3,2a,i2,2a,i2,a,',esymm%nclass,'i2,a)'
 write(msg,fmt)ch10,&
&  ' ===== Character of bands at k-point: ',esymm%kpt,' ===== ',ch10,&
&  '   Total number of bands analyzed .................. ',esymm%nbnds,ch10,&
&  '   Number of degenerate sets detected .............. ',esymm%ndegs,ch10,&
&  '   Number of operations in the little group of k ... ',esymm%nsym_gk,ch10,&
&  '   Number of classes (irreps) in the group of k .... ',esymm%nclass,' (',(esymm%nelements(icl),icl=1,esymm%nclass),' )'
 call wrtout(my_unt,msg,my_mode)

 if (esymm%nonsymmorphic_at_zoneborder) then
   call wrtout(my_unt," Non-symmorphic small group at zone border. Character analysis not available ",my_mode)
 end if

 if (esymm_failed(esymm)) then
   write(std_out,'(3a)')"Band classification algorithm failed with the error:",ch10,TRIM(esymm%err_msg)
   write(msg,'(3a)')"Band classification algorithm failed with the error:",ch10,TRIM(esymm%err_msg)
   call wrtout(my_unt,msg,my_mode)
 end if

 !nunknown=0
 !do iband=1,esymm%nbnds
 !  irr_idx = esymm%b2irrep(iband)
 !  if (irr_idx /= 0) then
 !    if (     esymm%has_chtabs) irr_name = esymm%Ref_Irreps(irr_idx)%name
 !    if (.not.esymm%has_chtabs) write(irr_name,'(i0)')irr_idx ! use the index instead of the name.
 !  else
 !    irr_name = "???"
 !    nunknown = nunknown +1
 !  end if
 !  write(msg,'(a,i3,2a)')' Band ',iband,' belongs to irrep ',TRIM(irr_name)
 !  call wrtout(my_unt,msg,my_mode)
 !end do

 do irr_idx=1,esymm%nclass
   nstates = size(esymm%irrep2b(irr_idx)%value)
   if (esymm%has_chtabs) then
     write(msg0,'(a,i0,3a)')"  Found ",nstates," states with character ",TRIM(esymm%Ref_irreps(irr_idx)%name),": "
   else
     write(msg0,'(2(a,i0),a)')"   Found ",nstates," states with character index ",irr_idx,": "
   end if
   do istart=1,nstates,20
     istop=istart+11; if (istop>nstates) istop=nstates
     write(msg,'(20(1x,i0))')(esymm%irrep2b(irr_idx)%value(ii), ii=istart,istop)
     if (istart==1) msg = TRIM(msg0)//TRIM(msg)
     if (istart/=1) msg = "   "//TRIM(msg)
     call wrtout(my_unt,msg,my_mode)
   end do
 end do

 nunknown = size(esymm%irrep2b(0)%value)
 if (nunknown > 0) then
   write(msg0,'(a,i0,a)')" WARNING: ",nunknown," states have not been classified:"
   do istart=1,nunknown,20
     istop=istart+11; if (istop>nunknown) istop=nunknown
     write(msg,'(20(1x,i0))')(esymm%irrep2b(0)%value(ii), ii=istart,istop)
     if (istart==1) msg = TRIM(msg0)//TRIM(msg)
     if (istart/=1) msg = "   "//TRIM(msg)
     call wrtout(my_unt,msg,my_mode)
   end do
 end if

 if (my_prtvol>0 .or. nunknown>0 .or. .not.esymm%has_chtabs) then ! print the calculated character table.
   call wrtout(my_unt,ch10//" Calculated character table ",my_mode)
   !write(fmt,*)'(i2,a,i2,1x,',esymm%nclass,'(a,2f6.3),a)'
   write(fmt,*)'(i2,a,i2,1x,',esymm%nclass,'(a,2f5.2),a)'
   do idg=1,esymm%ndegs
     write(msg,fmt)&
&      esymm%degs_bounds(1,idg),'-',esymm%degs_bounds(2,idg),&
&      ('|',esymm%Calc_irreps(idg)%trace(esymm%nelements(icl)), icl=1,esymm%nclass),'|'
     call wrtout(my_unt,msg,my_mode)
   end do
 end if

end subroutine esymm_print
!!***

!----------------------------------------------------------------------

!!****f* m_esymm/esymm_free_0D
!! NAME
!! esymm_free_0D
!!
!! FUNCTION
!!  Deallocate the memory allocated in the esymm_t datatype (scalar version)
!!
!! PARENTS
!!      m_esymm
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine esymm_free_0D(esymm)

!Arguments ------------------------------------
!scalars
 type(esymm_t),intent(inout) :: esymm

!Local variables ------------------------------
!scalars
 integer :: ii
! *************************************************************************

 !@esymm_t
 ABI_SFREE(esymm%g0)
 ABI_SFREE(esymm%tr_g0)
 ABI_SFREE(esymm%nelements)
 ABI_SFREE(esymm%sgk2symrec)
 ABI_SFREE(esymm%tr_sgk2symrec)
 ABI_SFREE(esymm%herring_test)
 ABI_SFREE(esymm%b2irrep)
 ABI_SFREE(esymm%degs_bounds)
 ABI_SFREE(esymm%degs_dim)

 if (allocated(esymm%irrep2b)) then
   do ii=LBOUND(esymm%irrep2b,DIM=1),UBOUND(esymm%irrep2b,DIM=1)
     ABI_FREE(esymm%irrep2b(ii)%value)
   end do
   ABI_DT_FREE(esymm%irrep2b)
 end if

 if (allocated(esymm%Calc_irreps)) call irrep_free(esymm%Calc_irreps)
 if (allocated(esymm%trCalc_irreps)) call irrep_free(esymm%trCalc_irreps)
 if (allocated(esymm%Ref_irreps)) call irrep_free(esymm%Ref_irreps)

end subroutine esymm_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_esymm/esymm_free_2D
!! NAME
!! esymm_free_2D
!!
!! FUNCTION
!!  Deallocate the memory allocated in the esymm_t datatype (2D version)
!!
!! PARENTS
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine esymm_free_2D(esymm)

!Arguments ------------------------------------
!scalars
 type(esymm_t),intent(inout) :: esymm(:,:)

!Local variables ------------------------------
 integer :: id1,id2
! *************************************************************************

 do id2=1,SIZE(esymm,DIM=2)
   do id1=1,SIZE(esymm,DIM=1)
     call esymm_free_0D(esymm(id1,id2))
   end do
 end do

end subroutine esymm_free_2D
!!***

!----------------------------------------------------------------------

!!****f* m_esymm/esymm_finalize
!! NAME
!! esymm_finalize
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      classify_bands
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine esymm_finalize(esymm,prtvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: prtvol
 type(esymm_t),target,intent(inout) :: esymm

!Local variables-------------------------------
 integer :: idg,ib1,ib2,idx,nunknown,dg_dim
 integer :: try,irep,nitems,nseen,isn
 integer :: isym,idg1,idg2,dim_mat,irr_idx2,irr_idx1
 real(dp),parameter :: TOL_TRACE=0.1_dp,TOL_ORTHO=0.1_dp,TOL_UNITARY=0.1_dp ! Large tolerance is needed to avoid problems.
 !real(dp),parameter :: TOL_TRACE=0.01_dp,TOL_ORTHO=0.01_dp,TOL_UNITARY=0.01_dp ! Large tolerance is needed to avoid problems.
 !real(dp),parameter :: TOL_TRACE=tol3,TOL_ORTHO=tol3,TOL_UNITARY=tol3 ! Large tolerance is needed to avoid problems.
 real(dp) :: uerr,max_err
 complex(dpc) :: ctest
 logical :: isnew
 character(len=500) :: msg
!arrays
 integer,allocatable :: dims_seen(:)
 complex(dpc),allocatable :: traces_seen(:,:)
 complex(dpc),pointer :: trace(:)
 complex(dpc),pointer :: calc_mat(:,:),trace1(:),trace2(:)
 complex(dpc),allocatable :: cidentity(:,:)

! *************************************************************************

 !@esymm_t

 ! Each band is initialized as "Unknown".
 esymm%b2irrep = 0

 ! Force the matrices to be unitary.
 call polish_irreps(esymm%Calc_irreps)

 if (.not.esymm%has_chtabs) then

   write(msg,'(5a)')&
&    "Reference character table not available. ",ch10,&
&    "Symmetry analysis not available. Using heuristic method to classify the states.",ch10,&
&    "It might not work, especially if accidental degeneracies are present."
   MSG_WARNING(msg)
   !
   ! The simplest thing we can do here is using the calculated matrices to get the
   ! character and comparing the results hoping everything is OK.
   ABI_MALLOC(traces_seen,(esymm%nsym_gk,esymm%ndegs))
   ABI_MALLOC(dims_seen,(esymm%ndegs))

   traces_seen=czero; nseen=1
   traces_seen(:,1) = esymm%Calc_irreps(1)%trace
   dims_seen(1)     = esymm%Calc_irreps(1)%dim

   do idg=2,esymm%ndegs
     dg_dim = esymm%Calc_irreps(idg)%dim
     trace => esymm%Calc_irreps(idg)%trace
     isnew=.TRUE.
     do isn=1,nseen
       if (ALL (ABS(trace - traces_seen(:,isn)) < TOL_TRACE) ) then
         isnew=.FALSE.; EXIT
       end if
     end do

     if (isnew) then
       nseen = nseen+1
       traces_seen(:,nseen) = trace
       dims_seen(nseen) = dg_dim
     end if
   end do

   if (nseen>esymm%nclass) then
     write(msg,'(3a)')&
&      "The number of different calculated traces is found to be greater than nclasses!",ch10,&
&      "Heuristic method clearly failed. Symmetry analysis cannot be performed."
     MSG_WARNING(msg)
     esymm%err_status = ESYM_HEUR_WRONG_NCLASSES
     esymm%err_msg    = msg

     do isn=1,nseen
       write(msg,'(a,i0)')" Representation: ",isn
       call wrtout(std_out,msg,"COLL")
       call print_arr(traces_seen(:,isn),max_r=esymm%nsym_gk,unit=std_out,mode_paral="COLL")
     end do

   else  ! It seems that the Heuristic method succeeded.
     do idg=1,esymm%ndegs
       ib1=esymm%degs_bounds(1,idg)
       ib2=esymm%degs_bounds(2,idg)
       trace => esymm%Calc_irreps(idg)%trace
       do isn=1,nseen
         if (ALL (ABS(trace - traces_seen(:,isn)) < TOL_TRACE) ) then
           esymm%b2irrep(ib1:ib2)=isn
           if (esymm%Calc_irreps(idg)%dim /= dims_seen(isn)) then
             write(msg,'(3a)')&
&              "Found two set of degenerate states with same character but different dimension!",ch10,&
&              "heuristic method clearly failed. Symmetry analysis cannot be performed."
             MSG_ERROR(msg)
             esymm%err_status = ESYM_HEUR_WRONG_DIMS
             esymm%err_msg    = msg
           end if
           EXIT
         end if
       end do
     end do
   end if

   ABI_FREE(traces_seen)
   ABI_FREE(dims_seen)

 else
   !
   ! * Search in the lookup table definining the irreducible representation
   nunknown = 0
   do idg=1,esymm%ndegs

     ib1=esymm%degs_bounds(1,idg)
     ib2=esymm%degs_bounds(2,idg)
     trace => esymm%Calc_irreps(idg)%trace

     try = which_irrep(esymm, trace, tol3)
     if (try==0) try = which_irrep(esymm, trace, 0.1_dp) ! try again with increased tolerance.
     if (try/=0) then
       esymm%b2irrep(ib1:ib2)=try
     else
       esymm%b2irrep(ib1:ib2)=0
       nunknown = nunknown + (ib2-ib1+1)
     end if
   end do
 end if
 !
 ! %irrep2b(0)) gives the indeces of the states that have not been classified.
 ABI_DT_MALLOC(esymm%irrep2b,(0:esymm%nclass))

 !write(std_out,*)"b2irrep",esymm%b2irrep

 do irep=0,esymm%nclass
   nitems = COUNT(esymm%b2irrep==irep)
   ABI_MALLOC(esymm%irrep2b(irep)%value,(nitems))
   idx=0
   do ib1=1,esymm%nbnds
     if (esymm%b2irrep(ib1) == irep) then
       idx = idx + 1
       esymm%irrep2b(irep)%value(idx) = ib1
     end if
   end do
 end do

 if (size(esymm%irrep2b(0)%value) /= 0) then
   write(msg,'(a,i0,a)')" Band classification algorithm was not able to classify ",size(esymm%irrep2b(0)%value)," states."
   MSG_WARNING(msg)
   esymm%err_status = ESYM_CLASSIFICATION_ERROR
   esymm%err_msg    = msg
 end if
 !
 ! ==============================================================
 ! ==== Test basic properties of irreducible representations ====
 ! ==============================================================

 if (.not.esymm_failed(esymm)) then
   !
   ! 1) \sum_R \chi^*_a(R)\chi_b(R)= N_R \delta_{ab}
   !
   !call wrtout(std_out," \sum_R \chi^*_a(R)\chi_b(R) = N_R \delta_{ab} ","COLL")
   max_err=zero
   do idg2=1,esymm%ndegs
     trace2 => esymm%Calc_irreps(idg2)%trace(1:esymm%nsym_gk)
     ib2 = esymm%degs_bounds(1,idg2)
     irr_idx2 = esymm%b2irrep(ib2)
     if (irr_idx2 == 0) CYCLE

     do idg1=1,idg2
       trace1 => esymm%Calc_irreps(idg1)%trace(1:esymm%nsym_gk)
       ib1 = esymm%degs_bounds(1,idg1)
       irr_idx1 = esymm%b2irrep(ib1)
       if (irr_idx1 == 0) CYCLE
       ctest=DOT_PRODUCT(trace1,trace2)/esymm%nsym_gk
       if (irr_idx1==irr_idx2) ctest=ctest-one
       max_err = MAX(max_err,ABS(ctest))
       if (.FALSE..and.ABS(ctest)>tol3) then
         write(msg,'(a,4i3,2es16.8)')&
&          ' WARNING: should be delta_ij: cx1 cx2, irr1, irr2, ctest: ',idg1,idg2,irr_idx1,irr_idx2,ctest
         call wrtout(std_out,msg,"COLL")
       end if
     end do
   end do

   if (max_err>TOL_ORTHO) then
     write(msg,'(a,es10.2)')" Too large maximum error on \sum_R \chi^*_a(R)\chi_b(R) = N_R \delta_{ab}: ",max_err
     MSG_WARNING(msg)
     esymm%err_status =  ESYM_ORTHO_ERROR
     esymm%err_msg    =  msg
   else
     write(msg,'(a,es10.2)')" maximum error on \sum_R \chi^*_a(R)\chi_b(R) = N_R \delta_{ab}: ",max_err
     call wrtout(std_out,msg,"COLL")
   end if

   if (.not.esymm%only_trace) then
     !call wrtout(std_out," **** Testing the unitary of the calculated irreps ****",my_mode)
     max_err=zero
     do idg1=1,esymm%ndegs
       ib1 = esymm%degs_bounds(1,idg1)
       irr_idx1 = esymm%b2irrep(ib1)
       if (irr_idx1 == 0) CYCLE

       do isym=1,esymm%nsym_gk
         calc_mat => esymm%Calc_irreps(idg1)%mat(:,:,isym)
         dim_mat  =  esymm%Calc_irreps(idg1)%dim
         ABI_MALLOC(cidentity,(dim_mat,dim_mat))
         call set2unit(cidentity)
         uerr = MAXVAL( ABS(MATMUL(calc_mat,TRANSPOSE(DCONJG(calc_mat))) - cidentity) )
         max_err = MAX(max_err,uerr)
         ABI_FREE(cidentity)
         if (.FALSE..and.prtvol>=10) then
           write(std_out,'(a,i3,a,i2,a,es16.8,a)')&
&          " === idg: ",idg1,", isym: ",isym,", Error on U^* U = 1: ",uerr," ==="
           call print_arr(calc_mat,dim_mat,dim_mat,unit=std_out,mode_paral="COLL")
         end if
       end do
     end do

     if (max_err>TOL_UNITARY) then
       write(msg,'(a,es10.2)')" Too large maximum error on the unitary of representions matrices: ",max_err
       MSG_WARNING(msg)
       esymm%err_msg    = msg
       esymm%err_status = ESYM_UNITARY_ERROR
     else
       write(msg,'(a,es10.2)')" maximum error on the unitary of representions matrices: ",max_err
       call wrtout(std_out,msg,"COLL")
     end if

   end if

 end if

end subroutine esymm_finalize
!!***

!----------------------------------------------------------------------

!!****f* m_esymm/which_irrep
!! NAME
!!  m_esymm
!!
!! FUNCTION
!!  Return the index of the irreducible representation with character charact. 0 if not found.
!!
!! INPUTS
!!  esymm<esymm_t>
!!  trace(%nsym_gk)=The trace of the representation to be compared with the internal database (if present).
!!  tolerr=Absolute error on the character.
!!
!! PARENTS
!!
!! SOURCE

function which_irrep(esymm,trace,tolerr)

!Arguments ------------------------------------
!scalars
 integer :: which_irrep
 real(dp),intent(in) :: tolerr
 type(esymm_t),intent(in) :: esymm
!arrays
 complex(dpc),intent(in) :: trace(esymm%nsym_gk)

!Local variables-------------------------------
!scalars
 integer :: irp

! *********************************************************************

 which_irrep = 0
 if (esymm%has_chtabs) then ! Symmetry analysis can be performed.
   do irp=1,esymm%nclass
     if ( ALL( ABS(esymm%Ref_irreps(irp)%trace(:) - trace(:)) < tolerr)) then
       which_irrep = irp; EXIT
     end if
   end do
 end if

end function which_irrep
!!***

!----------------------------------------------------------------------


!!****f* m_esymm/esymm_symmetrize_mels
!! NAME
!!  esymm_symmetrize_mels
!!
!! FUNCTION
!!
!! INPUTS
!!  esymm<esymm_t>
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cohsex_me
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine esymm_symmetrize_mels(esymm,lbnd,ubnd,in_me,out_me)

!Arguments ------------------------------------
!scalars
 integer :: lbnd,ubnd
 type(esymm_t),target,intent(in) :: esymm
!arrays
 complex(dpc),intent(in) :: in_me(2,lbnd:ubnd,lbnd:ubnd)
 complex(dpc),intent(out) :: out_me(lbnd:ubnd,lbnd:ubnd)

!Local variables-------------------------------
!scalars
 integer :: idg1,b1_start,b1_stop,irp1
 integer :: idg2,b2_start,b2_stop,irp2
 integer :: ii,jj,ib,jb,kk,kb,lb,ll
 complex(dpc) :: tr_ofd,ofd,dsd,tr_dsd
 type(irrep_t),pointer :: Irrep1,Irrep2
 type(irrep_t),pointer :: tr_Irrep1,tr_Irrep2

! *********************************************************************

 if (esymm_failed(esymm)) then
   MSG_ERROR("Symmetrization cannot be performed. You should not be here!")
 end if

 do idg1=1,esymm%ndegs  ! First loop over set of degenerate states.
   b1_start = esymm%degs_bounds(1,idg1)
   b1_stop  = esymm%degs_bounds(2,idg1)

   !if (b1_stop<lbnd .or. b2_start >ubnd) then
   !  MSG_ERROR("Wrong band indeces, check esymm initialization")
   !end if

   Irrep1 => esymm%Calc_irreps(idg1)
   if (esymm%can_use_tr) tr_Irrep1 => esymm%trCalc_irreps(idg1)
   irp1 = esymm%b2irrep(b1_start)

   do idg2=1,esymm%ndegs ! Second loop over set of degenerate states.
     !write(std_out,*)" ==> Symmetrizing degenerate set ",idg1,idg2
     b2_start = esymm%degs_bounds(1,idg2)
     b2_stop  = esymm%degs_bounds(2,idg2)
     irp2 = esymm%b2irrep(b2_start)

     if (irp1/=irp2 .or. idg1==idg2) CYCLE  ! Skip diago elements or elements belonging to different irreps.

     Irrep2 => esymm%Calc_irreps(idg2)
     if (esymm%can_use_tr) tr_Irrep2 => esymm%trCalc_irreps(idg2)
     !
     ! Symmetrize the off-diagonal matrix elements.
     ! summing over kk and ll. ii and jj are the indeces of the bands that are symmetrized
     do ii=1,b1_stop-b1_start+1
       ib= ii+b1_start-1
       do jj=1,b2_stop-b2_start+1
         jb= jj+b2_start-1
         !write(std_out,*)" ====> Symmetrizing ",ib,jb

         ofd= czero; tr_ofd=czero
         do kk=1,b1_stop-b1_start+1
           kb= kk+b1_start-1
           do ll=1,b2_stop-b2_start+1
             lb= ll+b2_start-1
             dsd = sum_irreps(Irrep1,Irrep2,kk,ii,ll,jj)
             ofd = ofd + dsd * in_me(1,kb,lb)
             if (esymm%can_use_tr) then
               tr_dsd = sum_irreps(tr_Irrep1,tr_Irrep2,kk,jj,ll,ii) ! Exchange of band indeces.
               tr_ofd = tr_ofd + tr_dsd * in_me(2,kb,lb)            ! Contribution obtained from TR.
             end if
           end do
         end do

         out_me(ib,jb)= ofd/esymm%nsym_gk
         if (esymm%can_use_tr .and. esymm%nsym_trgk>0) out_me(ib,jb)= out_me(ib,jb) + tr_ofd/esymm%nsym_trgk
       end do
     end do
   end do
 end do

end subroutine esymm_symmetrize_mels
!!***

!----------------------------------------------------------------------

!!****f* m_esymm/esymm_failed
!! NAME
!!  esymm_failed
!!
!! FUNCTION
!!
!! INPUTS
!!  esymm<esymm_t>
!!
!! PARENTS
!!
!! SOURCE

function esymm_failed(esymm)

!Arguments ------------------------------------
!scalars
 logical :: esymm_failed
 type(esymm_t),intent(in) :: esymm

! *********************************************************************

 esymm_failed = (esymm%err_status /= ESYM_NOERROR)

end function esymm_failed
!!***

!----------------------------------------------------------------------

!!****f* m_esymm/polish_irreps
!! NAME
!!  polish_irreps
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_esymm
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine polish_irreps(Irreps)

!Arguments ------------------------------------
!scalars
 type(irrep_t),intent(inout) :: Irreps(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: ldvl1=1,ldvr1=1
 integer :: irp,sym,dim,ldvr,ii,ivec,jvec,info
 !character(len=500) :: msg
!arrays
 complex(dpc),allocatable :: vl(:,:),vr(:,:),vrm1(:,:),overlap(:,:)
 complex(dpc),allocatable :: cmat(:,:),eigval(:)

! *********************************************************************

 ! Eigen decomposition: A = V D V^{-1}.
 do irp=1,SIZE(Irreps)
   dim = Irreps(irp)%dim
   ABI_MALLOC(cmat,(dim,dim))
   ABI_MALLOC(eigval,(dim))
   ldvr=dim
   ABI_MALLOC(vl,(ldvl1,dim))
   ABI_MALLOC(vr,(ldvr,dim))
   ABI_MALLOC(vrm1,(dim,dim))
   ABI_MALLOC(overlap,(dim,dim))
   do sym=1,Irreps(irp)%nsym
     cmat = Irreps(irp)%mat(:,:,sym)
     call xgeev("No vectors","Vectors",dim,cmat,dim,eigval,vl,ldvl1,vr,ldvr)
     !
     ! Orthogonalize the eigenvectors using Cholesky orthogonalization.
     do jvec=1,dim
       do ivec=1,jvec
         overlap(ivec,jvec) = DOT_PRODUCT(vr(:,ivec),vr(:,jvec))
       end do
     end do
     !
     ! 2) Cholesky factorization: overlap = U^H U with U upper triangle matrix.
     call ZPOTRF('U',dim,overlap,dim,info)
     ABI_CHECK(info == 0, sjoin('ZPOTRF returned info=', itoa(info)))

     ! 3) Solve X U = Vr, on exit the Vr treated by this node is orthonormalized.
     call ZTRSM('R','U','N','N',dim,dim,cone,overlap,dim,vr,dim)

     !write(std_out,*)"After ortho",MATMUL(TRANSPOSE(CONJG(vr)),vr)

     vrm1 = vr
     call xginv(vrm1,dim)
     do ii=1,dim
       eigval(ii) = eigval(ii)/ABS(eigval(ii)) ! Rescale the eigevalues.
       vrm1(ii,:) =  eigval(ii) * vrm1(ii,:)
     end do
     Irreps(irp)%mat(:,:,sym) = MATMUL(vr,vrm1)
     Irreps(irp)%trace(sym) = get_trace(Irreps(irp)%mat(:,:,sym))
   end do
   ABI_FREE(cmat)
   ABI_FREE(eigval)
   ABI_FREE(vl)
   ABI_FREE(vr)
   ABI_FREE(vrm1)
   ABI_FREE(overlap)
 end do

end subroutine polish_irreps
!!***

!----------------------------------------------------------------------

END MODULE m_esymm

