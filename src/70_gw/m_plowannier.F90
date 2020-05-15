!!****m* ABINIT/m_plowannier
!! NAME
!!  m_plowannier
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2020 ABINIT group (BAmadon,AGerossier,ROuterovitch)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_plowannier


#ifndef HAVE_CRPA_OPTIM
#ifdef FC_INTEL
#if  __INTEL_COMPILER<=1700
!DEC$ NOOPTIMIZE
#endif
#endif
#endif

 use defs_basis
 use m_errors
 use m_abicore
 use m_dtset
 use m_dtfil
 use defs_wvltypes
 use m_xmpi

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes, only : MPI_type
 use m_io_tools,  only : open_file
 use m_mpinfo,    only : proc_distrb_cycle
 use m_crystal, only : crystal_t
 use m_pawtab, only : pawtab_type
 use m_pawcprj, only : pawcprj_type,pawcprj_alloc,pawcprj_get,pawcprj_free
 use m_pawrad, only : pawrad_type, simp_gen


 implicit none

 private

 public :: init_plowannier
 public :: copy_orbital
 public :: compute_coeff_plowannier
 public :: destroy_plowannier
 public :: print_plowannier
 public :: get_plowannier
 public :: fullbz_plowannier
 public :: initialize_operwan
 public :: destroy_operwan
 public :: zero_operwan
 public :: compute_oper_ks2wan
 public :: normalization_plowannier
 public :: print_operwan
 public :: init_operwan_realspace
 public :: reduce_operwan_realspace
 public :: destroy_operwan_realspace
 public :: zero_operwan_realspace
 public :: compute_oper_wank2realspace
!!***


!!****t* m_plowannier/latom_wan_type
!! NAME
!!  latom_wan_type
!!
!! FUNCTION
!!
!!
!! SOURCE

type, public :: latom_wan_type

  integer, allocatable :: lcalc(:)
  ! array of the l we want to compute the psichi with

end type latom_wan_type
!!***


!!****t* m_plowannier/projector_wan_type
!! NAME
!!  projector_wan_type
!!
!! FUNCTION
!!
!!
!! SOURCE

type, public :: projector_wan_type

  integer, allocatable :: lproj(:)
  ! gives the list of the projector chosen

end type projector_wan_type
!!***


!!****t* m_plowannier/position_wan_type
!! NAME
!!  position_wan_type
!!
!! FUNCTION
!!
!!
!! SOURCE

type, public :: position_wan_type

  integer, allocatable :: pos(:,:)
  ! size (number of position,3)

end type position_wan_type
!!***

!!****t* m_plowannier/lorbital_type
!! NAME
!!  lorbital_type
!!
!! FUNCTION
!!
!!
!! SOURCE

type, public :: lorbital_type

  complex(dpc), allocatable :: matl(:,:,:)
  !details for different m

  real(dp), allocatable :: ph0phiint(:)
  ! stocks the values for each projector of the l considered
  ! size(total number of projectors for this l)

end type lorbital_type
!!***

!!****t* m_plowannier/orbital_type
!! NAME
!!  orbital_type
!!
!! FUNCTION
!!
!!
!! SOURCE

type, public :: orbital_type

  type(lorbital_type), allocatable :: atom(:)
  ! details of the psichi coefficients for each atom
  ! size of number of l chosen

end type orbital_type
!!***


!!****t* m_plowannier/lorbital2_type
!! NAME
!!  lorbital2_type
!!
!! FUNCTION
!!
!!
!! SOURCE

type, public :: lorbital2_type

  complex(dpc), allocatable :: matl(:,:,:,:,:)
  ! size (2l1+1,2l2+1,nspppol,nspinor,nspinor)

  real(dp), allocatable :: ph0phiint(:)
  ! size (nproj), stocks the value of ph0phiint we may want


end type lorbital2_type
!!***

!!****t* m_plowannier/operwan_type
!! NAME
!!  operwan_type
!!
!! FUNCTION
!!
!!
!! SOURCE

type, public :: operwan_type

  type(lorbital2_type), allocatable :: atom(:,:)
  ! l chosen for each on of both atoms

end type operwan_type
!!***

!!****t* m_plowannier/atom_index_type
!! NAME
!!  atom_index_type_type
!!
!! FUNCTION
!!
!!
!! SOURCE

type, public :: atom_index_type

  type(operwan_type), allocatable :: position(:,:)
  ! size (number of positions chosen for atom1, number of positions chosen for atom2)

end type atom_index_type
!!***

!!****t* m_plowannier/operwan_realspace_type
!! NAME
!!  operwan_realspace_type
!!
!! FUNCTION
!!
!!
!! SOURCE

type, public :: operwan_realspace_type

  type(atom_index_type), allocatable :: atom_index(:,:)
  ! size (number of atom, number of atom)

end type operwan_realspace_type
!!***

!!****t* m_plowannier/plowannier_type
!! NAME
!!  plowannier_type
!!
!! FUNCTION
!!
!!
!! SOURCE

type, public :: plowannier_type

  integer :: nkpt
  ! number of k points in Brillouin zone

  integer :: bandi_wan
  ! energy band minimum considered

  integer :: bandf_wan
  ! energy band maximum considered

  integer :: natom_wan
  ! number of atoms (used to compute Wannier functions)

  integer :: size_wan
  ! sum of all the m possible for every atom considered

  integer, allocatable :: iatom_wan(:)
  ! array of each atom (used to compute Wannier functions)
 
  integer, allocatable :: nbl_atom_wan(:)
  ! array of the number of l considered for each atom

  type(latom_wan_type), allocatable :: latom_wan(:)
  ! for each atom, it contains an array of the l we are interested in

  integer, allocatable :: nbproj_atom_wan(:)
  ! array of the number of projectors considered for each atom

  type(projector_wan_type), allocatable :: projector_wan(:)
  ! for each atom, it contains an array of the projectors we are interested in

  type(position_wan_type), allocatable :: nposition(:)
  ! array of the number of position considered for each atom

  integer :: nsppol
  ! number of polarization

  integer :: nspinor
  ! number of spinorial components

  type(orbital_type), allocatable :: psichi(:,:,:)
  ! arrays of psichi

  integer, allocatable :: position(:,:)
  ! size natom,3, gives the position of the cell for this atom (rprim coordinates)

  real(dp),allocatable :: kpt(:,:)
  ! gives the coordinates in the BZ of the kpoint
  ! size (3,nkpt)

  real(dp),allocatable :: wtk(:)
  !weight of each kpoint

  real(dp),allocatable :: acell(:)
  !size of the cell

end type plowannier_type
!!***

CONTAINS  !========================================================================================*
!!***


!!***f* m_plowannier/init_plowannier
!! NAME
!!  init_plowannier
!!
!! FUNCTION
!!  initialize the variables useful for the computation
!!
!! INPUTS
!! INPUTS
!! plowan_bandf  = max index of band for Wannier construction
!! plowan_bandi  = min index of band for Wannier construction
!! plowan_compute = keyword to activate Wannier calculation
!! plowan_iatom(plowan_natom) = index of atoms to use for Wannier
!! plowan_it(plowan_nt)= index of atoms for real space calculation
!! plowan_lcalc(sum_plowan_natom Plowan_nbl()) = index of l value for Wannier construction
!! plowan_natom = nb of atoms for Wannier
!! plowan_nbl(plowan_natom) = nb of l values for Wannier for each atoms.
!! plowan_nt = nb of atoms for real space calculation
!! plowan_projcalc(sum_plowan_natom Plowan_nbl()) = index of projectors for Wannier construction
!! acell_orig(3,nimage) = cell parameters
!! kpt(3,nkpt)  = k-points
!! nkpt = nb of k-points
!! nimage
!! nspinor = nb of spinors
!! nsppol = nb of polarization of wfc. 
!! wtk = weight of k-points 
!!
!! OUTPUT
!!  wan : plowannier type
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!! outscfcv,screening_driver,sigma_driver,fullbz_plowannier
!!
!! CHILDREN
!!
!! SOURCE


subroutine init_plowannier(plowan_bandf,plowan_bandi,plowan_compute,plowan_iatom,plowan_it,&
&plowan_lcalc,plowan_natom,plowan_nbl,plowan_nt,plowan_projcalc,acell_orig,kpt,nimage,nkpt,&
&nspinor,nsppol,wtk,wan)

!Arguments ----------------------------------
!scalars
! type(dataset_type), intent(in) :: dtset
 integer,intent(in) ::plowan_bandi,plowan_bandf,plowan_natom,plowan_nt,plowan_compute
 integer,intent(in) ::nkpt,nsppol,nspinor,nimage
 integer,intent(in) ::plowan_iatom(plowan_natom)
 integer,intent(in) ::plowan_nbl(plowan_natom)
 integer,intent(in) ::plowan_lcalc(sum(plowan_nbl(:)))
 integer,intent(in) ::plowan_projcalc(sum(plowan_nbl(:)))
 integer,intent(in) ::plowan_it(plowan_nt*3)
 real(dp),intent(in) :: kpt(3,nkpt)
 real(dp),intent(in) :: wtk(nkpt)
 real(dp),intent(in) :: acell_orig(3,nimage)
 type(plowannier_type), intent(inout) :: wan

!Local --------------------------------------
 integer :: iatom,ikpt,ib,iband,il,iltot,it,ittot,ltemp,nn,norbtot
 character(len=500) :: message
!************************************************************************

 !! generally
 wan%nkpt = nkpt
 wan%bandi_wan = plowan_bandi
 wan%bandf_wan = plowan_bandf
 wan%nsppol = nsppol
 wan%nspinor = nspinor

 !! for this case
 wan%natom_wan = plowan_natom

 !! generally
 ABI_ALLOCATE(wan%kpt,(3,size(kpt,2)))
  wan%kpt = kpt
 ABI_ALLOCATE(wan%iatom_wan,(wan%natom_wan))
 ABI_ALLOCATE(wan%nbl_atom_wan,(wan%natom_wan))
 wan%nbl_atom_wan = 0
 ABI_DATATYPE_ALLOCATE(wan%latom_wan,(wan%natom_wan))
 ABI_ALLOCATE(wan%nbproj_atom_wan,(wan%natom_wan))
 wan%nbproj_atom_wan = 0
 ABI_DATATYPE_ALLOCATE(wan%projector_wan,(wan%natom_wan))
 ABI_ALLOCATE(wan%position,(wan%natom_wan,3))
 wan%position = 0
 ABI_ALLOCATE(wan%wtk,(size(wtk,1)))
 wan%wtk(:) = wtk(:)
 ABI_ALLOCATE(wan%acell,(3))
 wan%acell(1) = acell_orig(1,1)
 wan%acell(2) = acell_orig(2,1)
 wan%acell(3) = acell_orig(3,1)

 ! If we want to study twice the same atom (but at different positions), use the same iatom and modify the positions below.
 ! In this case, the Wannier functions will be orthonormalized for one atom.
 ! For this particular reason, if we use the study twice the same atom at different position, each of them should have exactly the same projectors (it could be improved though by rewriting the normalization routine).



 ABI_DATATYPE_ALLOCATE(wan%nposition,(wan%natom_wan))
     !write(std_out,*)  "plowan_it", dtset%plowan_it

 iltot=0
 do iatom=1,wan%natom_wan

   wan%iatom_wan(iatom)       = plowan_iatom(iatom)
   wan%nbl_atom_wan(iatom)    = plowan_nbl  (iatom)
   wan%nbproj_atom_wan(iatom) = plowan_nbl  (iatom)

  ! Now we define for each atom the selected orbital moments.
   ABI_ALLOCATE(wan%latom_wan(iatom)%lcalc,(wan%nbl_atom_wan(iatom)))
   ABI_ALLOCATE(wan%projector_wan(iatom)%lproj,(wan%nbproj_atom_wan(iatom)))
   do il=1,wan%nbl_atom_wan(iatom)
     iltot=iltot+1
     wan%latom_wan(iatom)%lcalc(il)=plowan_lcalc(iltot)
     wan%projector_wan(iatom)%lproj(il)=plowan_projcalc(iltot)
   enddo

  !For each iatom , pos is an array of two dimensions. The first one is
  !the number of lattice translation and the second one is ist
  !coordinates.
   ABI_ALLOCATE(wan%nposition(iatom)%pos,(plowan_nt,3))
   ittot=0
   do it=1,plowan_nt
     wan%nposition(iatom)%pos(it,1) = plowan_it(ittot+1)
     wan%nposition(iatom)%pos(it,2) = plowan_it(ittot+2)
     wan%nposition(iatom)%pos(it,3) = plowan_it(ittot+3)
     !write(std_out,*)  "position",wan%nposition(iatom)%pos(it,:)
     ittot=ittot+3
   enddo
 enddo

 !!generally
 ABI_DATATYPE_ALLOCATE(wan%psichi,(wan%nkpt,wan%bandf_wan-wan%bandi_wan+1,wan%natom_wan))
 do ikpt = 1,wan%nkpt
   do iband = wan%bandi_wan,wan%bandf_wan
     ib=iband-wan%bandi_wan+1
     do iatom = 1,wan%natom_wan
       ABI_DATATYPE_ALLOCATE(wan%psichi(ikpt,ib,iatom)%atom,(wan%nbl_atom_wan(iatom)))
       do il = 1,wan%nbl_atom_wan(iatom)
         nn=(2*wan%latom_wan(iatom)%lcalc(il)+1)
         ABI_ALLOCATE(wan%psichi(ikpt,ib,iatom)%atom(il)%matl,(nn,wan%nsppol,wan%nspinor))
         wan%psichi(ikpt,ib,iatom)%atom(il)%matl = zero
       end do
     end do
   end do
 end do
 do iatom = 1,wan%natom_wan
   do il = 1,wan%nbl_atom_wan(iatom)
     ABI_ALLOCATE(wan%psichi(1,1,iatom)%atom(il)%ph0phiint,(10)) ! max number of proj for l =10..
   end do
 end do


!sum of all the m possible
 wan%size_wan = 0
 do iatom = 1,wan%natom_wan
   do il = 1,wan%nbl_atom_wan(iatom)
     ltemp = wan%latom_wan(iatom)%lcalc(il)
     wan%size_wan = wan%size_wan + 2*ltemp + 1
   end do
 end do

   write(message,'(2a,i5,i5)') ch10,&
&   ' == Lower and upper values of the selected bands',wan%bandi_wan,wan%bandf_wan
   call wrtout(std_out,message,'COLL') ; call wrtout(ab_out,message,'COLL')
   write(message,'(a,i10)')  ' == Number of atoms                             ',wan%natom_wan
   call wrtout(std_out,message,'COLL') ; call wrtout(ab_out,message,'COLL')
   write(message,'(a,9i2)')  ' == Atoms selected                               ',(wan%iatom_wan(ltemp),ltemp=1,wan%natom_wan)
   call wrtout(std_out,message,'COLL') ; call wrtout(ab_out,message,'COLL')
   write(message,'(a,9i2)')  ' == Nb of angular momenta used for each atom     ',(wan%nbl_atom_wan(ltemp),ltemp=1,wan%natom_wan)
   call wrtout(std_out,message,'COLL') ; call wrtout(ab_out,message,'COLL')
   norbtot=0
   do iatom=1,wan%natom_wan
     write(message,'(a,i2,a,9i2)')  ' == Value of the angular momenta for atom',iatom,' is : ',&
  &   (wan%latom_wan(iatom)%lcalc(ltemp),ltemp=1,wan%nbl_atom_wan(iatom))
     call wrtout(std_out,message,'COLL') ; call wrtout(ab_out,message,'COLL')
     do ltemp=1,wan%nbl_atom_wan(iatom)
       norbtot=norbtot+2*(wan%latom_wan(iatom)%lcalc(ltemp))+1
     enddo
     write(message,'(a,i2,a,9i2)') ' == Value of the projectors      for atom',iatom,' is : ', &
  &    (wan%projector_wan(iatom)%lproj(ltemp),ltemp=1,wan%nbl_atom_wan(iatom))
     call wrtout(std_out,message,'COLL') ; call wrtout(ab_out,message,'COLL')
   enddo
   if(norbtot>wan%bandf_wan-wan%bandi_wan+1) then
     write(message,'(3a,2i6)') "  Number of wannier functions is larger than" ,&
     &" number of Kohn Sham bands used for Wannier functions: decrease the number of Wannier functions", &
     &" or increase the number of bands ",norbtot,wan%bandf_wan-wan%bandi_wan+1
     !MSG_ERROR(message)
   endif
   if(plowan_compute==2) then
     write(message,'(3a)')  ch10,' == plowan_compute=2 => off diag blocks in the k-space Wannier Hamiltonian matrix',&
    &                          'is put to zero before diagonalisation'
     call wrtout(std_out,message,'COLL') ; call wrtout(ab_out,message,'COLL')
   endif

end subroutine init_plowannier
!!***


!!****f* m_plowannier/copy_orbital
!! NAME
!!  copy_orbital
!!
!! FUNCTION
!!  Copy an array of orbital_type
!!
!! INPUTS
!!  lorbital1
!!
!! OUTPUT
!!  lorbital2
!!
!! PARENTS
!!      m_plowannier
!!
!! CHILDREN
!!
!! SOURCE



subroutine copy_orbital(orbital1,orbital2,n1,n2,n3)

 !Arguments----------------
 integer,intent(in) :: n1,n2,n3
 type(orbital_type), intent(in) :: orbital1(n1,n2,n3)
 type(orbital_type),intent(inout) :: orbital2(n1,n2,n3)

 !Local variable-----------
 integer :: n4,n5,n6,n7
 integer :: i,j,k,l,m,p,q

 do i = 1,n1
   do j = 1,n2
     do k = 1,n3
       n4 = size(orbital1(i,j,k)%atom,1)
       do l = 1,n4
         n5 = size(orbital1(i,j,k)%atom(l)%matl,1)
         n6 = size(orbital1(i,j,k)%atom(l)%matl,2)
         n7 = size(orbital1(i,j,k)%atom(l)%matl,3)
         do m = 1,n5
           do p = 1,n6
             do q = 1,n7
               orbital2(i,j,k)%atom(l)%matl(m,p,q) = orbital1(i,j,k)%atom(l)%matl(m,p,q)
             end do
           end do
         end do
       end do
     end do
   end do
 end do

end subroutine copy_orbital
!!***

!!****f* m_plowannier/allocate_orbital
!! NAME
!!  allocate_orbital
!!
!! FUNCTION
!!  allocate an array of orbital_type
!!
!! INPUTS
!!  lorbital1
!!
!! OUTPUT
!!  lorbital2
!!
!! PARENTS
!!      m_plowannier
!!
!! CHILDREN
!!
!! SOURCE



subroutine allocate_orbital(orbital1,orbital2,n1,n2,n3)

 !Arguments----------------
 integer,intent(in) :: n1,n2,n3
 type(orbital_type), intent(in) :: orbital1(n1,n2,n3)
 type(orbital_type),intent(inout) :: orbital2(n1,n2,n3)

 !Local variable-----------
 integer :: n4,n5,n6,n7
 integer :: i,j,k,l

 do i = 1,n1
   do j = 1,n2
     do k = 1,n3
       n4 = size(orbital1(i,j,k)%atom,1)
       ABI_DATATYPE_ALLOCATE(orbital2(i,j,k)%atom,(n4))
       do l = 1,n4
         n5 = size(orbital1(i,j,k)%atom(l)%matl,1)
         n6 = size(orbital1(i,j,k)%atom(l)%matl,2)
         n7 = size(orbital1(i,j,k)%atom(l)%matl,3)
         ABI_ALLOCATE(orbital2(i,j,k)%atom(l)%matl,(n5,n6,n7))
       end do
     end do
   end do
 end do

end subroutine allocate_orbital
!!***


!!****f* m_plowannier/destroy_orbital
!! NAME
!!  destroy_orbital
!!
!! FUNCTION
!!  destroy an array of orbital_type
!!
!! INPUTS
!!  lorbital1
!!
!! OUTPUT
!!  lorbital2
!!
!! PARENTS
!!      m_plowannier
!!
!! CHILDREN
!!
!! SOURCE



subroutine destroy_orbital(orbital2,n1,n2,n3)

 !Arguments----------------
 integer,intent(in) :: n1,n2,n3
 type(orbital_type),intent(inout) :: orbital2(n1,n2,n3)

 !Local variable-----------
 integer :: n4
 integer :: i,j,k,l

 do i = 1,n1
   do j = 1,n2
     do k = 1,n3
       n4 = size(orbital2(i,j,k)%atom,1)
       do l = 1,n4
         ABI_DEALLOCATE(orbital2(i,j,k)%atom(l)%matl)
       end do
       ABI_DATATYPE_DEALLOCATE(orbital2(i,j,k)%atom)
     end do
   end do
 end do

end subroutine destroy_orbital
!!***


!!***f* m_plowannier/compute_coeff_plowannier
!! NAME
!!  compute_coeff_plowannier
!!
!! FUNCTION
!!  Compute the coefficient
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!        -gprimd(3,3)=dimensional reciprocal space primitive translations
!!        -indsym(4,nsym,natom)=indirect indexing array for atom labels
!!        -symrec(3,3,nsym)=symmetry operations in reciprocal space
!!        -nsym= number of symetry operations
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  dimcprj(natom) = dimension for cprj
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  fermie= Fermi energy
!!  mband=maximum number of bands
!!  mbandcprj=
!!  mkmem =number of k points treated by this node
!!  mpi_enreg=information about MPI parallelization
!!  nkpt=number of k points.
!!  my_nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol) = occupancies of KS states.
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  usecprj=
!!  unpaw=file number for cprj
!!  nbandkss
!!  dtfil
!!
!! OUTPUT
!!  wan%psichi: projections <Psi|chi>
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!
!! SOURCE


subroutine compute_coeff_plowannier(cryst_struc,cprj,dimcprj,dtset,eigen,fermie,& 
& mpi_enreg,occ,wan,pawtab,psps,usecprj,unpaw,pawrad,dtfil)


 use m_hide_lapack

!Arguments ------------------------------------
!scalars

 type(plowannier_type),intent(inout) :: wan
 integer,intent(in) :: unpaw,usecprj
 real(dp),intent(in) :: fermie
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(crystal_t),intent(in) :: cryst_struc
!arrays
 integer, intent(in) :: dimcprj(cryst_struc%natom)
 real(dp),intent(in) :: eigen(dtset%mband*wan%nkpt*wan%nsppol)
 real(dp),intent(in) :: occ(dtset%mband*wan%nkpt*wan%nsppol)
 type(pawcprj_type), intent(in) :: cprj(cryst_struc%natom,wan%nspinor*dtset%mband*dtset%mkmem*wan%nsppol*usecprj)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(datafiles_type),intent(in) :: dtfil

!Local variables-------------------------------
!scalars
 integer :: band_index,dimpsichi,facpara
 integer :: iatom,iatom1,iatom2,iband,ibandc,ibg,ierr,ikpt
 integer :: iband1,owrunt,opt
 integer :: ilmn,iorder_cprj,ispinor,isppol,itypat,ilmn2
 integer :: lmn_size
 integer :: m1,maxnproju,me,natom,nband_k,nband_k_cprj
 integer :: nnn,nprocband,spaceComm
 integer :: plowan_greendos,plowan_hybrid,plowan_inter,plowan_computegreen
 real(dp) :: ph0phiint_used
 character(len=500) :: message
 character(len=50) :: mat_writing2,mat_writing2_out
 character(len=5000) :: mat_writing,mat_writing_out
 integer :: l1,count,mesh_size,il,count_total,l,proj
 integer :: il1,il2,im1,im2,index_c,index_l,ispinor1,ispinor2,sizem,pos1,pos2
 real(dp) :: int_current,sum,sum2,sum3

 complex(dpc) :: wbase,wcurrent
 real(dp) :: resolution, wincrease,wmax,wmin
 integer :: iw,dos,shift,unt,unt2,dos_unt,dos_unt2
 integer :: number_of_frequencies,band_struct,prtocc,prtint
 real(dp) :: convert
 complex(dpc):: xsum
 character(len=fnlen) :: owrfile

!arrays
 real(dp) :: chinorm
 complex(dpc), allocatable :: Fff(:)
 complex(dpc), allocatable :: buffer1(:)
 logical :: lprojchi
 type(pawcprj_type),allocatable :: cwaveprj(:,:)
 type(operwan_type), allocatable :: operwan(:,:,:)
 type(operwan_realspace_type) :: operwan_realspace
 type(operwan_realspace_type) :: operocc
 complex(dpc), allocatable :: eigenks(:,:,:,:)
 complex(dpc), allocatable :: operks(:,:,:,:)
 complex(dpc), allocatable :: identityks(:,:,:,:)
 real(dp), allocatable :: ff(:)
 complex(dpc), allocatable :: operwansquare(:,:,:,:)
 complex(dpc), allocatable :: operwansquarereal(:,:,:)
 complex(dpc), allocatable :: matrix_to_diag(:,:)
 complex(dpc), allocatable :: energies(:,:)
 complex(dpc), allocatable :: Ffftable(:,:)
 character(len = 5) :: i2s,x1

!To diagonalize eigenvalues
 real(dp), allocatable :: eig(:), rwork(:)
 complex(dpc), allocatable :: zwork(:)
 integer :: lwork,info,whole_diag
 !complex(dpc), allocatable :: densmat(:,:)
!************************************************************************

! Drive the normalization of the psichis

if (dtset%plowan_projcalc(1)==-2)then
  if (dtset%ucrpa >= 1 .or. dtset%dmft_kspectralfunc==1) then 
    opt = 0
  else
    opt=1
  endif
else
  opt=0
end if

if (opt==0) then
  write(message,*)ch10,"Normalization of plowannier k-point by k-point"
else
  write(message,*)ch10,"Normalization of plowannier on the sum of the k-points"
endif

MSG_COMMENT(message)
!opt=1
        ! 0 : normalization k-point by k-point (normal use of plowan)
        ! 1 : normalization of the sum over k-points (use with crpa old keywords)


! Internal variables (could be put one day as input variables of ABINIT).
 plowan_computegreen   = 0  !
              ! 0 : do nothing do not compute hybri or dos
              ! 1 : Compute hybridization or dos (depends on following keywords)
              ! 2 : not tested, probably with bug included.: compute hybridization in Wannier basis

 ! If plowan_computegreen>0, the following  data is useful
 plowan_greendos = 1  ! For the first atom, plowan_greendos is the index
              ! of the angular momentum in array
              ! wan%latom_wan(iatom)%lcalc: it is thus betwween 0 and  wan%nbl_atom_wan(iatom)
              ! It is not the value of the angular momentum but its
              ! index
 plowan_hybrid   = 0  !  Same convention as for greendos
 plowan_inter    = 1  ! compute all interaction between all atoms all orbitals all neighbours, requires plowan_realspace=1

 !owrfile = trim(dtfil%filnam_ds(4))//"_operwan_realspace"
 owrfile = "__operwan_realspace__"

 dos = 0
 if(plowan_greendos>0) dos=plowan_greendos
 if(plowan_hybrid>0)   dos=-plowan_hybrid
 if(plowan_hybrid>0.and.plowan_greendos>0) then
   write(message,*) " plowan_hybrid and plowan_greendos cannot be both >0"
   MSG_ERROR(message)
 endif

 ! GREEN STUDY PARAMETERS (FREQUENCIES)
 ! ===================================
 !!To choose the frequencies in the Green study
 wmin=-2.d0 ! eV
 wmax= 2.d0 ! eV
 resolution= 0.02  ! eV
 wbase = cmplx(wmin/27.2107,0.001,kind=dp) ! most negative value of frequency
 wincrease = resolution/27.2107  ! step
 number_of_frequencies = int((wmax-wmin)/resolution)



 ! Select the way of diagonalisation
 !===================================
 whole_diag = 1! 1 for diagonalization of the whole matrix, 0 for each orbital
 if(dtset%plowan_compute==2)  whole_diag = 0 ! off diagonal blocks are suppressed in the hamiltonian matrix before diagonalisation
 band_struct = 1 ! 1 for plotting band struct (Wannier bands)

 ! Select the real space calculation of Wannier function: Interpolation
 ! versus Analysis
 !===================================
 prtocc = 0 ! occupations have no meaning for a k-point path so the default is 0
 if(dtset%kptopt>0.and.dtset%plowan_realspace>=1) prtocc = 1  !1 to print the occupation in real space


 ! Select if computation of interactions is done
 !===================================
 prtint = plowan_inter !1 to print sqrt(sum of interaction squared) between orbitals, we do not use the input file to do this one


 ! Plot KS band structure
 !===================================
! !data for printing KS bands
! do ikpt = 1,wan%nkpt
!   print* ,'bandstruct', real(eigen(1+(ikpt-1)*30:30+(ikpt-1)*30))!, real(eigen(7471+(ikpt-1)*30:7500+(ikpt-1)*30))
! end do



!DBG_ENTER("COLL")
!Fake test to keep fermie as argument. REMOVE IT AS SOON AS POSSIBLE ...
 if(fermie>huge(zero))chinorm=zero

 facpara=1 !mpi_enreg%nproc
 if(abs(dtset%pawprtvol)>=3) then
   write(message,*) ch10, " number of k-points used is nkpt = ", dtset%nkpt
   call wrtout(std_out,  message,'COLL')
   write(message,*) " warning: parallelised version        ", dtset%nkpt
   call wrtout(std_out,  message,'COLL')
   write(message,*) " weights k-points used is wtk = wtk"
   call wrtout(std_out,  message,'COLL')
 end if

 if(usecprj==0) then
   write(message,*) "  usecprj=0 : BUG in init_plowannier",usecprj
   MSG_BUG(message)
 end if

 if(wan%nspinor/=dtset%nspinor) then
   write(message,*) "  wan%nspinor=/dtset%nspinor, init_plowannier is not working in this case",&
&   wan%nspinor,dtset%nspinor
   MSG_ERROR(message)
 end if



!----------------------------------- MPI-------------------------------------

!Init parallelism
 spaceComm=mpi_enreg%comm_cell
 if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_kpt
 me=mpi_enreg%me_kpt

!----------------------------------- MPI-------------------------------------


 lprojchi=.false.
 lprojchi=.true.
 natom=cryst_struc%natom


 write(message,'(2a)') ch10,&
& '  == Prepare data for projected local orbital wannier function calculation  '
 call wrtout(std_out,message,'COLL')
 if(abs(dtset%pawprtvol)>=3) then
   write(message, '(a,a)' ) ch10,&
&   '---------------------------------------------------------------'
!  call wrtout(ab_out,message,'COLL');call wrtout(std_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,a,a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&   '  Print useful data (as a check)',ch10,&
&   '  - Overlap of KS wfc with atomic orbital inside sphere',ch10,&
&   '  - Eigenvalues',ch10,&
&   '  - Weights of k-points',ch10,&
&   '  - Number of spins ',ch10,&
&   '  - Number of states'
!  call wrtout(ab_out,message,'COLL');call wrtout(std_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,a)' ) ch10,&
&   '---------------------------------------------------------------'
 end if
 if(dtset%nstep==0) then
   message = 'nstep should be greater than 1'
   MSG_BUG(message)
 end if


!********************* Max Values for U terms.
!maxlpawu=0
 maxnproju=0
 do iatom=1,natom
   if(pawtab(dtset%typat(iatom))%lpawu.ne.-1 .and. pawtab(dtset%typat(iatom))%nproju.gt.maxnproju)&
&   maxnproju=pawtab(dtset%typat(iatom))%nproju
 end do
!*****************   in forlb.eig
 if(me.eq.0.and.abs(dtset%pawprtvol)>=3) then
   if (open_file('forlb.eig',message,newunit=unt,form='formatted',status='unknown') /= 0) then
     MSG_ERROR(message)
   end if
   rewind(unt)
   write(unt,*) " Number of bands,   spins, and  k-point; and spin-orbit flag"
   write(unt,*) dtset%mband,wan%nsppol,wan%nkpt,wan%nspinor,wan%bandi_wan,wan%bandf_wan
   write(unt,*) " For each k-point, eigenvalues for each band"
   write(unt,*) (dtset%wtk(ikpt)*facpara,ikpt=1,wan%nkpt)
   band_index=0
   do isppol=1,wan%nsppol
     write(unt,*) " For spin"
     write(unt,*)  isppol
     do ikpt=1,wan%nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*wan%nkpt)
       write(unt,*) " For k-point"
       write(unt,*)  ikpt
       do iband=wan%bandi_wan,wan%bandf_wan
         write(unt, '(2i6,4x,f20.15)' ) iband-wan%bandi_wan+1,ikpt,eigen(iband+band_index)*2.d0
       end do
       band_index=band_index+nband_k
     end do
   end do
   close(unt)
 end if ! proc=me



 !-----------------------------------------------------------------
 ! Allocates, computes and stocks the ph0phiint
 !-----------------------------------------------------------------
 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 do iatom = 1,wan%natom_wan
   do il = 1,wan%nbl_atom_wan(iatom)
     count = 0
     count_total = 0
     l = wan%latom_wan(iatom)%lcalc(il)
     proj = wan%projector_wan(iatom)%lproj(il)
     itypat = dtset%typat(wan%iatom_wan(iatom))
     lmn_size = pawtab(itypat)%lmn_size
     ! modif proj if proj == -2 ---> usemdft has been used
     if (proj==-2) then
       do ilmn = 1,lmn_size 
         if ( psps%indlmn(1,ilmn,itypat) .eq. l .and. psps%indlmn(2,ilmn,itypat) .eq. 0 .and. proj .eq. -2 ) then
           proj=psps%indlmn(5,ilmn,itypat)
         end if
       end do
     end if
     ! check if the choice of proj is coherent with the value of l
     do ilmn = 1,lmn_size
       if (psps%indlmn(1,ilmn,itypat).eq. l .and. psps%indlmn(2,ilmn,itypat) .eq. 0) then
         count_total = count_total+1 !!counts the number total of projector for this l
          if (psps%indlmn(5,ilmn,itypat) .eq. proj) then
           count = count+1 !!the projector chosen is in the right l
         end if
       end if
     end do
     if (count .eq. 0) then
       write(message,'(a)') " The projector choice is not consistent with the orbital l"
       MSG_ERROR(message)
     else !good choice of projector
       wan%psichi(1,1,iatom)%atom(il)%ph0phiint = zero
       do ilmn2 = 1,lmn_size
         if (psps%indlmn(1,ilmn2,itypat) .eq. l .and. psps%indlmn(2,ilmn2,itypat) .eq. 0) then
           mesh_size = pawtab(itypat)%mesh_size
           ABI_ALLOCATE(ff,(mesh_size))
           ff(1:mesh_size) = pawtab(itypat)%phi(1:mesh_size,proj)*pawtab(itypat)%phi(1:mesh_size,psps%indlmn(5,ilmn2,itypat))
!           ff(1:mesh_size) = pawtab(itypat)%tphi(1:mesh_size,proj)*pawtab(itypat)%tphi(1:mesh_size,psps%indlmn(5,ilmn2,itypat))
           int_current = 0
           call simp_gen(int_current,ff,pawrad(itypat)) !call the subroutine which does the computation
           wan%psichi(1,1,iatom)%atom(il)%ph0phiint(psps%indlmn(3,ilmn2,itypat)) = int_current !we put the values for ikpt = 1 and iband = 1
           ABI_DEALLOCATE(ff)
         end if
       end do
     end if
   end do
 end do
 !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


!==========================================================================
!***************** Compute  <Psi|Chi>=\sum_{proja} <Psi|P_a><phi_a|Chi>
!==========================================================================

!Allocate temporary cwaveprj storage
 ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,wan%nspinor))

 call pawcprj_alloc(cwaveprj,0,dimcprj)

 nprocband=(dtset%mband/dtset%mband)
 ibg=0
 do isppol=1,wan%nsppol
   do ikpt=1,wan%nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*wan%nkpt)
     nband_k_cprj=nband_k/nprocband
     !nband_k is mband for each k, so it is mband most of the time
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle
     do iband=wan%bandi_wan,wan%bandf_wan !loop only over the bands we are interested in
       ibandc=iband-wan%bandi_wan+1
!      Parallelization: treat only some bands
       if (dtset%paral_kgb==1) then
         if (mod((iband-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band)/=mpi_enreg%me_band) cycle
       else
         if (mpi_enreg%proc_distrb(ikpt,ibandc,isppol)/=me) cycle
       end if
       do ispinor=1,wan%nspinor
         do iatom = 1,wan%natom_wan !loop over the atom chosen
           itypat = dtset%typat(wan%iatom_wan(iatom))
           lmn_size = pawtab(itypat)%lmn_size !retrieve the number of different lmn
           call pawcprj_get(cryst_struc%atindx1,cwaveprj,cprj,natom,iband,ibg,ikpt,&
&          iorder_cprj,isppol,dtset%mband,dtset%mkmem,dtset%natom,1,nband_k_cprj,&
&          wan%nspinor,wan%nsppol,unpaw,mpicomm=mpi_enreg%comm_kpt,&
&          proc_distrb=mpi_enreg%proc_distrb)
           chinorm=1.d0
           do l1 = 1,wan%nbl_atom_wan(iatom) !l1 is the index of the orbital
             l = wan%latom_wan(iatom)%lcalc(l1) !l is the value of the orbital (linked to index l1)
             do ilmn = 1,lmn_size
               if (psps%indlmn(1,ilmn,itypat) .eq. l) then!
                 ph0phiint_used = wan%psichi(1,1,iatom)%atom(l1)%ph0phiint(psps%indlmn(3,ilmn,itypat))
                 m1 = psps%indlmn(2,ilmn,itypat)+l+1
                 wan%psichi(ikpt,ibandc,iatom)%atom(l1)%matl(m1,isppol,ispinor)=&
                   wan%psichi(ikpt,ibandc,iatom)%atom(l1)%matl(m1,isppol,ispinor)+&
                   cmplx(cwaveprj(wan%iatom_wan(iatom),ispinor)%cp(1,ilmn)*ph0phiint_used,cwaveprj(&
                   wan%iatom_wan(iatom),ispinor)%cp(2,ilmn)*ph0phiint_used,kind=dp)
               end if
             end do
           end do
         end do ! iatom
       end do ! ispinor
     end do !iband
     ibg=ibg+nband_k_cprj*wan%nspinor !useful to select the right ikpt in pawcprj_get
   end do !ikpt
 end do ! isppol


!===========================================================
!************************ new gather info for MPI
!===========================================================

 dimpsichi=0
 do iatom = 1,wan%natom_wan
   do l1 = 1,wan%nbl_atom_wan(iatom)
     dimpsichi = dimpsichi + wan%nkpt*(wan%bandf_wan-wan%bandi_wan+1)*(2*wan%latom_wan(iatom)%lcalc(l1)+1)*wan%nsppol*wan%nspinor
   end do
 end do
 dimpsichi = 2*dimpsichi !for complex
 ABI_ALLOCATE(buffer1,(dimpsichi))
 buffer1 = zero
 nnn = 0
 do ikpt = 1,wan%nkpt
   do ibandc = 1,wan%bandf_wan-wan%bandi_wan+1
     do iatom=1,wan%natom_wan
       do l1 = 1,wan%nbl_atom_wan(iatom)
         do m1 = 1,2*wan%latom_wan(iatom)%lcalc(l1)+1
           do isppol = 1,wan%nsppol
             do ispinor = 1,wan%nspinor
               nnn=nnn+1
               buffer1(nnn)=wan%psichi(ikpt,ibandc,iatom)%atom(l1)%matl(m1,isppol,ispinor)
             end do
           end do
         end do
       end do
     end do
   end do
 end do
 call xmpi_barrier(spaceComm)
 call xmpi_sum(buffer1,spaceComm,ierr)
 if (dtset%paral_kgb==1 .and. nprocband > 1) then
   call xmpi_sum(buffer1,mpi_enreg%comm_band,ierr) !build sum over band processors
 end if
 call xmpi_barrier(spaceComm)
 nnn = 0
 do ikpt = 1,wan%nkpt
   do ibandc = 1,wan%bandf_wan-wan%bandi_wan+1
     do iatom = 1,wan%natom_wan
       do l1 = 1,wan%nbl_atom_wan(iatom)
         do m1 = 1,2*wan%latom_wan(iatom)%lcalc(l1)+1
           do isppol = 1,wan%nsppol
             do ispinor = 1,wan%nspinor
               nnn=nnn+1
               wan%psichi(ikpt,ibandc,iatom)%atom(l1)%matl(m1,isppol,ispinor)=buffer1(nnn)
             end do
           end do
         end do
      end do
     end do
   end do
 end do
 ABI_DEALLOCATE(buffer1)

 call xmpi_barrier(spaceComm)



 !! -------------------------------------------------------------
 !! COMPUTATION OF THE OCCUPATION MATRIX BEFORE NORMALIZATION
 !! -------------------------------------------------------------
 !! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

 if (dtset%prtvol >= 5) then
   !Inialize an empty Wannier operator
   ABI_DATATYPE_ALLOCATE(operwan,(wan%nkpt,wan%natom_wan,wan%natom_wan))
   call initialize_operwan(wan,operwan)

   !Creation of the KS occupation operator
   ABI_ALLOCATE(eigenks,(wan%nkpt,wan%bandf_wan-wan%bandi_wan+1,wan%bandf_wan-wan%bandi_wan+1,wan%nsppol))
   ABI_ALLOCATE(identityks,(wan%nkpt,wan%bandf_wan-wan%bandi_wan+1,wan%bandf_wan-wan%bandi_wan+1,wan%nsppol))
   eigenks = czero
   identityks=czero
   do isppol = 1,wan%nsppol
     do iband1 = 1,wan%bandf_wan-wan%bandi_wan+1
      ibandc = iband1 + wan%bandi_wan - 1
       do ikpt = 1,wan%nkpt
        eigenks(ikpt,iband1,iband1,isppol) = occ(((ikpt-1)*dtset%mband+ibandc+(isppol-1)*wan%nkpt*dtset%mband))
         !write(6,*) 'eigenks', ikpt,iband1,isppol,((ikpt-1)*dtset%mband+ibandc+(isppol-1)*wan%nkpt*dtset%mband),occ(((ikpt-1)*dtset%mband+ibandc+(isppol-1)*wan%nkpt*dtset%mband))
       end do
     end do
   end do
   
   
   !compute the occupation in wannier basis and print it
  write(message,*)char(10),&
&" Print the occupation levels (not normalized) for 1 atom, 1 orbitals"
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   write(message,*)" Atom =",wan%iatom_wan(1),"orbital =",wan%latom_wan(1)%lcalc(1)
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   do ikpt = 1,wan%nkpt
     call compute_oper_ks2wan(wan,eigenks,operwan,ikpt)
   enddo
   call init_operwan_realspace(wan,operocc)
   write(message,*)char(10)," The occupation matrix before normalization is :"
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   call compute_oper_wank2realspace(wan,operwan,operocc)
    if (wan%nsppol ==1)then
     write(message,*)char(10)," For one spin :"
     call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')     
     do im1=1,2*wan%latom_wan(1)%lcalc(1)+1
       write(mat_writing,'(7f20.10)')real(operocc%atom_index(1,1)%position(1,1)%atom(1,1)%matl(im1,:,1,1,1))
       call wrtout(std_out,mat_writing,'COLL'); call wrtout(ab_out,mat_writing,'COLL')
     enddo
   else 
     write(message,*)char(10)," For spin up :"
     call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
     do im1=1,2*wan%latom_wan(1)%lcalc(1)+1
       write(mat_writing,'(7f20.10)')real(operocc%atom_index(1,1)%position(1,1)%atom(1,1)%matl(im1,:,1,1,1))
       call wrtout(std_out,mat_writing,'COLL'); call wrtout(ab_out,mat_writing,'COLL')
     enddo
     write(message,*)char(10)," For spin down : "
     call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
     do im1=1,2*wan%latom_wan(1)%lcalc(1)+1       
       write(mat_writing,'(7f20.10)')real(operocc%atom_index(1,1)%position(1,1)%atom(1,1)%matl(im1,:,2,1,1))
       call wrtout(std_out,mat_writing,'COLL'); call wrtout(ab_out,mat_writing,'COLL')
     enddo
   endif
 
 endif

 !!-------------------------------------------------------------
 !!NORMALIZATION
 !!--------------------------------------------------------------

 if(dtset%prtvol>=5) then 
   do isppol = 1,wan%nsppol
     do iband1 = 1,wan%bandf_wan-wan%bandi_wan+1
       ibandc = iband1 + wan%bandi_wan - 1
       do ikpt = 1,wan%nkpt
         identityks(ikpt,iband1,iband1,isppol) = one
!write(6,*) "idks", identityks(ikpt,iband1,iband1,isppol) 
       end do
     end do
   end do
   call zero_operwan(wan,operwan)
   do ikpt = 1,wan%nkpt
     call compute_oper_ks2wan(wan,identityks,operwan,ikpt)
   enddo
   do ikpt = 1,wan%nkpt
    ! if (ikpt<=5)then
       write(message,*)char(10)," For ikpt=",ikpt,"the normalization matrix is before normalization :"
       call wrtout(std_out,message,'COLL'); !call wrtout(ab_out,message,'COLL')
       do im1=1,2*wan%latom_wan(1)%lcalc(1)+1
         write(mat_writing,'(7f20.5)')real(operwan(ikpt,1,1)%atom(1,1)%matl(im1,:,1,1,1))
         call wrtout(std_out,mat_writing,'COLL'); !call wrtout(ab_out,mat_writing,'COLL')
       enddo
     !endif
   end do
 endif

 
 call normalization_plowannier(wan,opt)


 if (dtset%prtvol>=5) then
   call zero_operwan(wan,operwan)
   do ikpt = 1,wan%nkpt
     call compute_oper_ks2wan(wan,identityks,operwan,ikpt)
   enddo
   do ikpt = 1,wan%nkpt
     !if (ikpt<=5)then
       write(message,*)char(10)," For ikpt=",ikpt,"the normalization matrix is after normalization :"
       call wrtout(std_out,message,'COLL'); !call wrtout(ab_out,message,'COLL')
       do im1=1,2*wan%latom_wan(1)%lcalc(1)+1
         write(mat_writing,'(7f20.5)')real(operwan(ikpt,1,1)%atom(1,1)%matl(im1,:,1,1,1))
         call wrtout(std_out,mat_writing,'COLL'); !call wrtout(ab_out,mat_writing,'COLL')
       enddo
     !endif
   end do
 endif
 !! -------------------------------------------------------------
 !! COMPUTATION OF THE OCCUPATION MATRIX AFTER NORMALIZATION
 !! -------------------------------------------------------------
 !! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 !compute the occupation in wannier basis and print it

 if (dtset%prtvol >= 5) then  
   write(message,*)char(10),&
     &" Print the occupation levels (normalized) for 1 atom, 1 orbital"
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   write(message,*)"Atom =",wan%iatom_wan(1),"orbital =",wan%latom_wan(1)%lcalc(1)
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   call zero_operwan(wan,operwan)
   do ikpt = 1,wan%nkpt
     call compute_oper_ks2wan(wan,eigenks,operwan,ikpt)
   enddo
   call zero_operwan_realspace(wan,operocc)
   write(message,*)char(10)," The occupation matrix after normalization is :"
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   call compute_oper_wank2realspace(wan,operwan,operocc)
   if (wan%nsppol ==1)then
     write(message,*)char(10)," For one spin :"
     call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')     
     do im1=1,2*wan%latom_wan(1)%lcalc(1)+1
       write(mat_writing,'(7f20.10)')real(operocc%atom_index(1,1)%position(1,1)%atom(1,1)%matl(im1,:,1,1,1))
       call wrtout(std_out,mat_writing,'COLL'); call wrtout(ab_out,mat_writing,'COLL')
     enddo
   else 
     write(message,*)char(10)," For spin up :"
     call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
     do im1=1,2*wan%latom_wan(1)%lcalc(1)+1
       write(mat_writing,'(7f20.10)')real(operocc%atom_index(1,1)%position(1,1)%atom(1,1)%matl(im1,:,1,1,1))
       call wrtout(std_out,mat_writing,'COLL'); call wrtout(ab_out,mat_writing,'COLL')
     enddo
     write(message,*)char(10)," For spin down : "
     call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
     do im1=1,2*wan%latom_wan(1)%lcalc(1)+1       
       write(mat_writing,'(7f20.10)')real(operocc%atom_index(1,1)%position(1,1)%atom(1,1)%matl(im1,:,2,1,1))
       call wrtout(std_out,mat_writing,'COLL'); call wrtout(ab_out,mat_writing,'COLL')
     enddo
   endif


   mat_writing = ""
   
   if (me.eq.0) then
!print operwan in the real space, in a file
     mat_writing = trim(dtfil%filnam_ds(4))//"_wannierocc"
     convert = 1
     call print_operwan(wan,operwan,trim(mat_writing),convert)
   end if
! destroy operators and the occupation matrix
   ABI_DEALLOCATE(eigenks)
   ABI_DEALLOCATE(identityks)
   call destroy_operwan(wan,operwan)
   ABI_DATATYPE_DEALLOCATE(operwan)
   call destroy_operwan_realspace(wan,operocc)!!Destroy the occupation matrix 
 endif


!! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^








 !! -------------------------------------------------------------
 !! TO COMPUTE THE ENERGY MATRIX
 !! -------------------------------------------------------------
 !! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

 ! Initialize an empty Wannier operator
 ABI_DATATYPE_ALLOCATE(operwan,(wan%nkpt,wan%natom_wan,wan%natom_wan))
 call initialize_operwan(wan,operwan)

 ! Creation of the KS occupation operator
 ABI_ALLOCATE(operks,(wan%nkpt,wan%bandf_wan-wan%bandi_wan+1,wan%bandf_wan-wan%bandi_wan+1,wan%nsppol))
 operks = czero
 do isppol = 1,wan%nsppol
   do iband1 = 1,wan%bandf_wan-wan%bandi_wan+1
     ibandc = iband1 + wan%bandi_wan - 1
     do ikpt = 1,wan%nkpt
       operks(ikpt,iband1,iband1,isppol) = eigen(((ikpt-1)*dtset%mband+ibandc+(isppol-1)*wan%nkpt*dtset%mband))-fermie
     end do
   end do
 end do

 ! Compute the energy in wannier basis
 do ikpt = 1,wan%nkpt
   call compute_oper_ks2wan(wan,operks,operwan,ikpt)
 end do
!!In operwan, energies level are stored (shifted with fermi level) in Hartree
 ABI_DEALLOCATE(operks)


 ! check that the eigenvalues are real
 do ikpt = 1,wan%nkpt
   do iatom1 = 1,wan%natom_wan
     do il1 = 1,wan%nbl_atom_wan(iatom1)
       do m1= 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
         do isppol = 1,wan%nsppol
           if (aimag(operwan(ikpt,iatom1,iatom1)%atom(il1,il1)%matl(m1,m1,isppol,1,1)) > 1d-8) then
             write(mat_writing,'(a)') " An eigenvalue has an imaginary part: ikpt, atom, l, m, isppol, value"
             write(mat_writing2,'(i0,i0,i0,i0,i0,E15.6)') ikpt, iatom1, il1, im1, isppol, &
  &                   operwan(ikpt,iatom1,iatom1)%atom(il1,il1)%matl(m1,m1,isppol,1,1)
             MSG_ERROR(message)
             
           end if
         end do
       end do
     end do
   end do
 end do





 !! -------------------------------------------------------------
 !! Transform the Wannier operator in real space (in eV)
 !! 1) allocate operwan_realspace
 !! -------------------------------------------------------------
if (dtset%plowan_realspace >= 1) then
  call init_operwan_realspace(wan,operwan_realspace)
endif

! ABI_DATATYPE_ALLOCATE(operwan_realspace,(wan%natom_wan,wan%natom_wan))
! do iatom1 = 1,wan%natom_wan
!   do iatom2 = 1,wan%natom_wan
!     n1=size(wan%nposition(iatom1)%pos,1)
!     n2=size(wan%nposition(iatom2)%pos,1)
!     ABI_DATATYPE_ALLOCATE(operwan_realspace(iatom1,iatom2)%position,(n1,n2))
!     do pos1 = 1,size(wan%nposition(iatom1)%pos,1)
!       do pos2 = 1,size(wan%nposition(iatom2)%pos,1)
!         n1=wan%nbl_atom_wan(iatom1)
!         n2=wan%nbl_atom_wan(iatom2)
!         ABI_DATATYPE_ALLOCATE(operwan_realspace(iatom1,iatom2)%position(pos1,pos2)%atom,(n1,n2))
!         do il1 = 1,wan%nbl_atom_wan(iatom1)
!           do il2 = 1,wan%nbl_atom_wan(iatom2)
!             n1=2*wan%latom_wan(iatom1)%lcalc(il1)+1
!             n2=2*wan%latom_wan(iatom2)%lcalc(il2)+1
! ABI_ALLOCATE(operwan_realspace(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl,(n1,n2,wan%nsppol,1,1))
!             operwan_realspace(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl = zero
!           end do
!         end do
!       end do
!     end do
!   end do
! end do

 !! -------------------------------------------------------------
 !! Transform the Wannier operator in real space (in eV)
 !! 2) compute the value in real space (only if kptopt>0 ie BZ correctly sampled)
 !! -------------------------------------------------------------
 if (dtset%plowan_realspace >= 1 .and. dtset%kptopt > 0 ) then ! interpolation and kptopt >0 : compute Wannier functions in real space.
   call compute_oper_wank2realspace(wan,operwan,operwan_realspace)
!   do isppol = 1,wan%nsppol
!     do iatom1 = 1,wan%natom_wan
!       do pos1 = 1,size(wan%nposition(iatom1)%pos,1)
!         do il1 = 1,wan%nbl_atom_wan(iatom1)
!           do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
!             do iatom2 = 1,wan%natom_wan
!               do pos2 = 1,size(wan%nposition(iatom2)%pos,1)
!                 do il2 = 1,wan%nbl_atom_wan(iatom2)
!                   do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
!                     !sum over ikpt
!                     do ikpt = 1,wan%nkpt
!                       operwan_realspace%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl(im1,im2,isppol,1,1) =&
!                         operwan_realspace%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl(im1,im2,isppol,1,1)&
!                         + real(operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,1,1)&
!                         * wan%wtk(ikpt) * exp( cmplx(0.0,1.0) * two_pi * ( &
!                         wan%kpt(1,ikpt) * ( wan%nposition(iatom1)%pos(pos1,1) - wan%nposition(iatom2)%pos(pos2,1) )+&
!                         wan%kpt(2,ikpt) * ( wan%nposition(iatom1)%pos(pos1,2) - wan%nposition(iatom2)%pos(pos2,2) )+&
!                         wan%kpt(3,ikpt) * ( wan%nposition(iatom1)%pos(pos1,3) - wan%nposition(iatom2)%pos(pos2,3)))))
!                     end do
!                     !end of the sum
!                   end do
!                 end do
!               end do
!             end do
!           end do
!         end do
!       end do
!     end do
!   end do



 !! -------------------------------------------------------------
 !! Transform the Wannier operator in real space (in eV)
 !! 3) write operwan_realspace in a file unit owrunt
 !! -------------------------------------------------------------


   write(message,'(4a)') ch10,&
&   '  == Write hamiltonian in real space Wannier function to file ',trim(owrfile),' =='
   call wrtout(std_out,message,'COLL')
   if (me.eq.0) then
     if (open_file(owrfile, message, newunit=owrunt, form="unformatted", status="unknown", action="write") /= 0) then
       MSG_ERROR(message)
     end if
     rewind(owrunt)
     do isppol = 1,wan%nsppol
       do iatom1 = 1,wan%natom_wan
         do pos1 = 1,size(wan%nposition(iatom1)%pos,1)
           do il1 = 1,wan%nbl_atom_wan(iatom1)
             do iatom2 = 1,wan%natom_wan
               do pos2 = 1,size(wan%nposition(iatom2)%pos,1)
                 do il2 = 1,wan%nbl_atom_wan(iatom2)
                   write(owrunt) operwan_realspace%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl(:,:,isppol,1,1)
                 end do
               end do
             end do
           end do
         end do
       end do
     end do
     close(owrunt)
   end if
 end if ! plotwan_realspace>0 and kptopt>0

 !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 ! To change until green study (use operwan_realspace instead of the print subroutine)
!!  BA?


 if (dtset%plowan_realspace == 1) then ! We print the matrix of energy in eV
   if (me.eq.0) then
     !print operwan in the real space, in a file
     mat_writing = trim(dtfil%filnam_ds(4))//"_wanniereigen"
     convert = 27.2107
     call print_operwan(wan,operwan,trim(mat_writing),convert)
   end if
 end if





 !!========================================================================================
 !! Computation of the interaction  ( sum(l and l') t_ll' )
 !! This interaction is meaningful in the real space, a transformation is made in this loop
 !!========================================================================================

 if (3==4.and.prtint .eq. 1 .and. dtset%plowan_realspace == 1 ) then
   if (open_file(trim(dtfil%filnam_ds(4))//'_inter',message, newunit=unt) /= 0) then
     MSG_ERROR(message)
   end if
   write(unt,'(a,i0,a,F7.3)') "# nsppol = ",wan%nsppol," and acell = ",dtset%acell_orig(1,1)
   write(unt,'(a)') "# Interaction between an orbital and another (in Hartree) : isppol iatom1 pos1 iproj1 iatom2 pos2 iproj2 value"

!!!to compute all interactions
!   do isppol = 1,wan%nsppol
!     do iatom1 = 1,wan%natom_wan
!       do pos1 = 1,size(wan%nposition(iatom1)%pos,1)
!         do il1 = 1,wan%nbl_atom_wan(iatom1)
!           do iatom2 = 1,wan%natom_wan
!             do pos2 = 1,size(wan%nposition(iatom2)%pos,1)
!               do il2 = 1,wan%nbl_atom_wan(iatom2)
!                 if (iatom1 .ne. iatom2 .or. il1 .ne. il2 .or. pos1 .ne. pos2) then !not the same orbital on the same atom
!                   if (iatom1 .lt. iatom2 .or. iatom1 .eq. iatom2 .and. (pos1 .lt. pos2 .or. (pos1 .eq. pos2 .and. il1 .lt. il2))) then ! to print only once each interac tion
!                     sum = 0
!                     do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
!                       do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
!                         sum2 = 0
!                         do ikpt = 1,wan%nkpt
!                           sum2 = sum2 + abs(real(operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,1,1)*wan%wtk(ikpt)*exp(cmplx(0.0,1.0)*two_pi*(wan%kpt(1,ikpt)*(wan%nposition(iatom1)%pos(pos1,1)-wan%nposition(iatom2)%pos(pos2,1))+wan%kpt(2,ikpt)*(wan%nposition(iatom1)%pos(pos1,2)-wan%nposition(iatom2)%pos(pos2,2))+wan%kpt(3,ikpt)*(wan%nposition(iatom1)%pos(pos1,3)-wan%nposition(iatom2)%pos(pos2,3))))))
!                         end do
!                         sum = sum + sum2
!                       end do
!                     end do
!                     write(unt,'(i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,E15.6)') isppol," ",iatom1," ",pos1," ",wan%projector_wan(iatom1)%lproj(il1)," ",iatom2," ",pos2," ",wan%projector_wan(iatom2)%lproj(il2) ,sum
!                   end if
!                 end if
!               end do
!             end do
!           end do
!         end do
!       end do
!     end do
!   end do



   do isppol = 1,wan%nsppol
     do iatom2 = 1,wan%natom_wan
       do pos2 = 1,size(wan%nposition(iatom2)%pos,1)
         do il2 = 1,wan%nbl_atom_wan(iatom2)
           sum3 = 0
           do im1 = 1,7 ! for f orbitals
             do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
               sum = 0
               sum2 = 0
               do ikpt = 1,wan%nkpt
                 sum2 = sum2 + real(operwan(ikpt,1,iatom2)%atom(1,il2)%matl(im1,im2,isppol,1,1)*wan%wtk(ikpt)*&
&                              exp(cmplx(0.0,1.0)*two_pi*( &
&             wan%kpt(1,ikpt)*(wan%nposition(1)%pos(1,1)-wan%nposition(iatom2)%pos(pos2,1))+ &
&             wan%kpt(2,ikpt)*(wan%nposition(1)%pos(1,2)-wan%nposition(iatom2)%pos(pos2,2))+ &
&             wan%kpt(3,ikpt)*(wan%nposition(1)%pos(1,3)-wan%nposition(iatom2)%pos(pos2,3))  )))
                 sum = sum + real(operwan(ikpt,iatom2,iatom2)%atom(il2,il2)%matl(im2,im2,isppol,1,1)*wan%wtk(ikpt))
               end do
               sum2 = sum2**2
               sum3 = sum3 + sum2/sum
             end do
           end do
           write(unt,'(i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,E15.6)') isppol,&
          &         " ",1," ",1," ",7," ",iatom2," ",pos2," ",wan%projector_wan(iatom2)%lproj(il2) ,sum3
         end do
       end do
     end do
   end do
   close(unt)
 end if


 !! transformation back in reciprocal space with a limited number of neighbors
 !!==========================================================================
 if (dtset%plowan_realspace == 2) then
   !read from the file
   write(message,'(4a)') ch10,&
&   '  == Read hamiltonian in real space Wannier function on file ',trim(owrfile),' =='
   call wrtout(std_out,message,'COLL')

   if (open_file(owrfile, message, newunit=owrunt, form="unformatted", status="old", action="read") /= 0) then
     MSG_ERROR(message)
   end if
   rewind(owrunt)
   do isppol = 1,wan%nsppol
     do iatom1 = 1,wan%natom_wan
       do pos1 = 1,size(wan%nposition(iatom1)%pos,1)
         do il1 = 1,wan%nbl_atom_wan(iatom1)
           do iatom2 = 1,wan%natom_wan
             do pos2 = 1,size(wan%nposition(iatom2)%pos,1)
               do il2 = 1,wan%nbl_atom_wan(iatom2)
                 read(owrunt) operwan_realspace%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl(:,:,isppol,1,1)
               end do
             end do
           end do
         end do
       end do
     end do
   end do
   close(owrunt)

   !-----------------------------------------------------------------
   !Set the xx' interaction to 0
   ! In order to do a Wannier interpolation without a given number of
   ! real space terms in the wannier hamiltonian.
   !-----------------------------------------------------------------
   if(3==4) then
     do isppol = 1,wan%nsppol
       do pos2 = 2,size(wan%nposition(2)%pos,1)
         operwan_realspace%atom_index(1,2)%position(1,pos2)%atom(1,1)%matl(:,:,isppol,1,1) = zero
         operwan_realspace%atom_index(2,1)%position(pos2,1)%atom(1,1)%matl(:,:,isppol,1,1) = zero
        !             atoms 1 and 2 are selected for removal
        !             position(1,pos2): select cell 1 and all other cells  pos2.
        !             atom(1,1): selected index 1 of atom1 and index 1 of atom2
        !             mat1(:,:,isppol,1,1): remove all ml terms.
         operwan_realspace%atom_index(1,1)%position(1,pos2)%atom(1,2)%matl(:,:,isppol,1,1) = zero
       end do
     end do
     write(message,'(2a)') '  == Block suppressed in the real space Wannier hamiltonian'
     call wrtout(std_out,message,'COLL')
   endif



   !-----------------------------------------------------------------
   ! Perform Wannier transform from real space to reciprocal space
   !-----------------------------------------------------------------
   write(message,'(2a)') '  == Perform Wannier transform from real space to reciprocal space =='
   call wrtout(std_out,message,'COLL')
   do isppol = 1,wan%nsppol
     do ikpt = 1,wan%nkpt
       do iatom1 = 1,wan%natom_wan
         do il1 = 1,wan%nbl_atom_wan(iatom1)
           do iatom2 = 1,wan%natom_wan
             do il2 = 1,wan%nbl_atom_wan(iatom2)
               operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(:,:,isppol,1,1) = zero
               do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
                 do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                   !sum over neigbours
                   pos1 = 1
                     do pos2 = 1,size(wan%nposition(iatom2)%pos,1)
                       operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,1,1) =&
                         operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,1,1)&
                         +real(operwan_realspace%atom_index(iatom1,iatom2)%&
                         &position(pos1,pos2)%atom(il1,il2)%matl(im1,im2,isppol,1,1)&
                         * exp( - cmplx(0.0,1.0) * two_pi * ( &
                         wan%kpt(1,ikpt) * ( wan%nposition(iatom1)%pos(pos1,1) - wan%nposition(iatom2)%pos(pos2,1) )+&
                         wan%kpt(2,ikpt) * ( wan%nposition(iatom1)%pos(pos1,2) - wan%nposition(iatom2)%pos(pos2,2) )+&
                         wan%kpt(3,ikpt) * ( wan%nposition(iatom1)%pos(pos1,3) - wan%nposition(iatom2)%pos(pos2,3)))))
                     end do
                   !end do
                   !end of the sum
                 end do
               end do
             end do
           end do
         end do
       end do
     end do
   end do

 ! plowan_realspace==2 (Wannier interpolation)

 !! -------------------------------------------------------------
 !! End Transform the Wannier operator in real space (in eV)
 !! n) deallocate operwan_realspace
 !! ------------------------------------------------------------
!   do iatom1 = 1,wan%natom_wan
!     do iatom2 = 1,wan%natom_wan
!       do pos1 = 1,size(wan%nposition(iatom1)%pos,1)
!         do pos2 = 1,size(wan%nposition(iatom2)%pos,1)
!           do il1 = 1,wan%nbl_atom_wan(iatom1)
!             do il2 = 1,wan%nbl_atom_wan(iatom2)
!               ABI_DATATYPE_DEALLOCATE(operwan_realspace%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl)
!             end do
!           end do
!           ABI_DATATYPE_DEALLOCATE(operwan_realspace%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom)
!         end do
!       end do
!       ABI_DATATYPE_DEALLOCATE(operwan_realspace%atom_index(iatom1,iatom2)%position)
!     end do
!   end do
!   ABI_DATATYPE_DEALLOCATE(operwan_realspace%atom_index)
 endif
 if (dtset%plowan_realspace >= 1) then
   call destroy_operwan_realspace(wan,operwan_realspace)
 endif
 
 ! ----------------------------------------------------------------------------------------
 ! Here each block of the hamiltonian matrix in Wannier basis is diagonalized separately
 ! ----------------------------------------------------------------------------------------

! off diagonal blocks are suppressed in the hamiltonian matrix before diagonalisation
! whole_diag=-1
 if (whole_diag .eq. 0) then
   write(message,'(a,i5,a)')'  plowan_compute =',dtset%plowan_compute,&
 &  ' Off diag blocks are suppressed in the Wannier hamiltonian before diagonalisation'
   call wrtout(std_out,message,'COLL')
   ! !To diagonalize the block matrix for each orbital
   do iatom1 = 1,wan%natom_wan
     do il1 = 1,wan%nbl_atom_wan(iatom1)
       do isppol = 1,wan%nsppol
         count = 2*wan%latom_wan(iatom1)%lcalc(il1)+1
         ABI_ALLOCATE(matrix_to_diag,(count,count))
         ABI_ALLOCATE(eig,(count))
         ABI_ALLOCATE(rwork,(3*count-2))
         lwork = 65*count !Value to optimize the diagonalization
         ABI_ALLOCATE(zwork,(lwork))
         do ikpt = 1,wan%nkpt
           matrix_to_diag(:,:) = operwan(ikpt,iatom1,iatom1)%atom(il1,il1)%matl(:,:,isppol,1,1)
           call zheev('v','u',count,matrix_to_diag,count,eig,zwork,lwork,rwork,info)
           if (info .eq. 0) then !!Correct diagonalization
             matrix_to_diag = zero
             do im1 = 1,count
               matrix_to_diag(im1,im1) = eig(im1)
             end do
             operwan(ikpt,iatom1,iatom1)%atom(il1,il1)%matl(:,:,isppol,1,1) = matrix_to_diag(:,:)
           else
             write(message,'(a)') "Error in the normalization of the Wannier eigenvalues" ! BA?
             MSG_ERROR(message)
           end if
         end do
         ABI_DEALLOCATE(matrix_to_diag)
         ABI_DEALLOCATE(eig)
         ABI_DEALLOCATE(rwork)
         ABI_DEALLOCATE(zwork)
       end do
     end do
   end do
 end if



 ! ----------------------------------------------------------------------------------------
 ! Here the hamiltonian matrix in Wannier basis is diagonalized completely
 ! ----------------------------------------------------------------------------------------

 if (whole_diag .eq. 1) then
   ! To diagonalize the whole matrix
   count = wan%nspinor*wan%size_wan
   ABI_ALLOCATE(matrix_to_diag,(count,count))
   ABI_ALLOCATE(eig,(count))
   ABI_ALLOCATE(rwork,(3*count-2))
   lwork = 65*count ! Value to optimize speed of the diagonalization
   ABI_ALLOCATE(zwork,(lwork))
   !First, write operwan matrix in an inversible matrix
   do isppol = 1,wan%nsppol
     do ikpt = 1,wan%nkpt
       matrix_to_diag = czero
       index_l = 0
       do iatom1 = 1,wan%natom_wan
         do il1 = 1,wan%nbl_atom_wan(iatom1)
           do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
             index_l = index_l + 1 ! the line changes
             index_c = 1 ! index_c is set to one each time the line changes
             do iatom2 = 1,wan%natom_wan
               do il2 = 1,wan%nbl_atom_wan(iatom2)
                 do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                   matrix_to_diag(index_l,index_c) = operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,1,1)
                   index_c = index_c + 1
                 end do !im2
               end do !il2
             end do ! iatom2 (the line changes)
           end do ! im1
         end do ! il1
       end do !iatom1

       !Then, invert the matrix
      call zheev('v','u',count,matrix_to_diag,count,eig,zwork,lwork,rwork,info)
       if (info .eq. 0) then ! Correct diagonalization
         matrix_to_diag = czero
         do im1 = 1,count
           matrix_to_diag(im1,im1) = eig(im1)
         end do
       else
         write(message,'(a)') "Error in the normalization of the Wannier eigenvalues" ! BA?
         MSG_ERROR(message)
       end if
       !Finally, we write the value diagonalized back into operwan
       index_l = 0
       do iatom1 = 1,wan%natom_wan
         do il1 = 1,wan%nbl_atom_wan(iatom1)
           do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
             index_l = index_l + 1
             index_c = 1
             do iatom2 = 1,wan%natom_wan
               do il2 = 1,wan%nbl_atom_wan(iatom2)
                 do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                   operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,1,1) = matrix_to_diag(index_l,index_c)
                   index_c = index_c + 1
                 end do
               end do
             end do
           end do
         end do
       end do
       !ikpt/isppol
     end do
   end do
   ABI_DEALLOCATE(matrix_to_diag)
   ABI_DEALLOCATE(eig)
   ABI_DEALLOCATE(rwork)
   ABI_DEALLOCATE(zwork)
 end if

 !------------------------------------------------------------------------
 !  we write the band structure of the atoms in the Wannier basis (in eV)
 ! Warning:  In the case of the diagonalisation of the whole hamltonien
 ! the band structure is separated in as many files as atoms.
 !------------------------------------------------------------------------
 if (band_struct .eq. 1 ) then

   write(message,'(3a)') ch10,&
&   ' == For each k-point of the path, gives the eigenvalues (in eV) of the Hamiltonian in the Wannier basis'
   call wrtout(std_out,message,'COLL') ; call wrtout(ab_out,message,'COLL')
   write(message,'(2a,f13.4,a)') ch10,&
&   '   (The band structure is shifted by fermie =',fermie*27.211,' eV )'
   call wrtout(std_out,message,'COLL') ; call wrtout(ab_out,message,'COLL')

   write(message,'(a,i10)')  ' == Number of atoms                             ',wan%natom_wan
   do iatom1 = 1,wan%natom_wan
     i2s = '(I0)'               ! trick to add the atom number
     write(x1,i2s) iatom1       ! at the end of the filename
     if (wan%nsppol .eq. 1 .and. me.eq.0 ) then
       if (open_file(trim(dtfil%filnam_ds(4))//"_BANDSTRUCT"//trim(x1),message,newunit=unt) /= 0) then
         MSG_ERROR(message)
       end if
       write(unt,'(a,i0)') "#Wannier band structure for the atom ",iatom1
       write(message,'(2a,i2)') ch10," Wannier band structure for atom ",iatom1
       call wrtout(std_out,message,'COLL') ; call wrtout(ab_out,message,'COLL')
       do ikpt = 1,wan%nkpt
         mat_writing = "#ikpt ="
         write(mat_writing2,'(i0)') ikpt
         mat_writing = trim(mat_writing)//" "//trim(mat_writing2)
         write(unt,*) trim(mat_writing)
         write(mat_writing,'(i0)') ikpt
         write(mat_writing_out,'(i0)') ikpt
         do l1 = 1,wan%nbl_atom_wan(iatom1)
           do m1 = 1,2*wan%latom_wan(iatom1)%lcalc(l1)+1
             write (mat_writing2,'(F25.7)') 27.2107*real(operwan(ikpt,iatom1,iatom1)%atom(l1,l1)%matl(m1,m1,1,1,1))
             write (mat_writing2_out,'(F12.3)') 27.2107*real(operwan(ikpt,iatom1,iatom1)%atom(l1,l1)%matl(m1,m1,1,1,1))
             mat_writing = trim(mat_writing)//trim(mat_writing2)
             mat_writing_out = trim(mat_writing_out)//trim(mat_writing2_out)
           end do
         end do
         write(unt,*) trim(mat_writing)
         write(ab_out,*) trim(mat_writing_out)
         write(std_out,*) trim(mat_writing_out)
       end do
       close(unt)
     else if (wan%nsppol .eq. 2) then
       if (open_file(trim(dtfil%filnam_ds(4))//"_BANDSTRUCTUP"//trim(x1),message,newunit=unt) /= 0) then
         MSG_ERROR(message)
       end if
       if (open_file(trim(dtfil%filnam_ds(4))//"_BANDSTRUCTDN"//trim(x1),message,newunit=unt2) /= 0) then
         MSG_ERROR(message)
       end if
       write(unt,'(a,i0,a)') "#Wannier band structure for the atom ",iatom1, " polarization up"
       write(unt2,'(a,i0,a)') "#Wannier band structure for the atom ",iatom1, " polarization down"
       do isppol = 1,wan%nsppol
         do ikpt = 1,wan%nkpt
           mat_writing = "#ikpt ="
           write(mat_writing2,'(i0)') ikpt
           mat_writing = trim(mat_writing)//" "//trim(mat_writing2)
           if (isppol .eq. 1) then
             write(unt,*) trim(mat_writing)
           else
             write(unt2,*) trim(mat_writing)
           end if
           write(mat_writing,'(i0)') ikpt
           do l1 = 1,wan%nbl_atom_wan(iatom1)
             do m1 = 1,2*wan%latom_wan(iatom1)%lcalc(l1)+1
               write (mat_writing2,'(F12.7)') 27.2107*real(operwan(ikpt,iatom1,iatom1)%atom(l1,l1)%matl(m1,m1,isppol,1,1))
               mat_writing = trim(mat_writing)//trim(mat_writing2)
             end do
           end do
           if (isppol .eq. 1) then
             write(unt,*) trim(mat_writing)
           else
             write(unt2,*) trim(mat_writing)
           end if
           write(ab_out,*) trim(mat_writing)
           write(std_out,*) trim(mat_writing)
         end do
       end do
       close(unt)
       close(unt2)
     end if
   end do
 end if


 !! ----------------------------------------------
 !! GREEN STUDY
 !! ----------------------------------------------

 !!Here we only study the first atom

 if (plowan_computegreen .eq. 1 .and. me.eq.0 ) then
   if (dos .ge. 1) then ! compute partial DOS for l=dos
     if (open_file(trim(dtfil%filnam_ds(4))//"_dosfromgreen",message,newunit=dos_unt) /= 0) then
       MSG_ERROR(message)
     end if
     write(dos_unt,'(a)') "#DOS function for the first atom computed with the Green function"
     write(dos_unt,'(a,i0,a,i0)') "# l = ",wan%latom_wan(1)%lcalc(dos)," bands ; nsppol = ",wan%nsppol
     write(dos_unt,'(a)') "# frequency, DOS, isppol"

     if(wan%nsppol>=2) then
       if (open_file(trim(dtfil%filnam_ds(4))//"_dosfromgreen_b",message,newunit=dos_unt2) /= 0) then
         MSG_ERROR(message)
       end if
       write(dos_unt2,'(a)') "#DOS function for the first atom computed with the Green function"
       write(dos_unt2,'(a,i0,a,i0)') "# l = ",wan%latom_wan(1)%lcalc(dos)," bands ; nsppol = ",wan%nsppol
       write(dos_unt2,'(a)') "# frequency, DOS, isppol"
     endif
   end if

   if (dos .le. -1) then ! compute Hybri for l=|dos|
     if (open_file(trim(dtfil%filnam_ds(4))//"_hybridization",message,newunit=dos_unt) /= 0) then
       MSG_ERROR(message)
     end if
     write(dos_unt,'(a)') "#Hybridization obtained from the green function for the bands selected"
     write(dos_unt,'(a,i0,a,i0)') "# l = ",wan%latom_wan(1)%lcalc(abs(dos))," bands ; nsppol = ",wan%nsppol
     write(dos_unt,'(a)') "#isppol, frequency, F(m=0), F(m=1), F(m=2) ..."
     if(wan%nsppol>=2)   then
       if (open_file(trim(dtfil%filnam_ds(4))//"_hybridization_b",message, newunit=dos_unt2) /= 0) then
         MSG_ERROR(message)
       end if
       write(dos_unt2,'(a)') "#Hybridization obtained from the green function for the bands selected"
       write(dos_unt2,'(a,i0,a,i0)') "# l = ",wan%latom_wan(1)%lcalc(abs(dos))," bands ; nsppol = ",wan%nsppol
       write(dos_unt2,'(a)') "#isppol, frequency, F(m=0), F(m=1), F(m=2) ..."
     endif
   end if

   ! Method 0
   sizem = 2*wan%latom_wan(1)%lcalc(abs(dos))+1 !number of m for the l orbital we want
   ABI_ALLOCATE(energies,(sizem,wan%nsppol))
   energies = czero
   ABI_ALLOCATE(Fff,(2))
   Fff = czero

   !We put the eigenenergies of the first atom in an array to use them later
   do isppol = 1,wan%nsppol
     do ikpt = 1,wan%nkpt
       do m1 = 1,sizem
         energies(m1,isppol) = energies(m1,isppol) +&
 &                real(operwan(ikpt,1,1)%atom(abs(dos),abs(dos))%matl(m1,m1,isppol,1,1))*wan%wtk(ikpt)
       end do
     end do
   end do


   shift = 0
   do il1 = 1,abs(dos)-1
     shift = shift + 2*wan%latom_wan(1)%lcalc(il1)+1!shift for the right choice of indices for the lorbital chosen
   end do

   !We destroy the operwan which was used to compute energies before

   call destroy_operwan(wan,operwan)
   ABI_DATATYPE_DEALLOCATE(operwan)

   !-----------------------------------------------------------
   ! Loop over the frequencies to compute DOS or Hybridization
   !-----------------------------------------------------------
   do iw = 1,number_of_frequencies
     ABI_DATATYPE_ALLOCATE(operwan,(wan%nkpt,wan%natom_wan,wan%natom_wan))
     call initialize_operwan(wan,operwan)
     !!creation of the Green operator
     wcurrent = wbase + (iw-1)*wincrease
    ! if (allocated(operks)) then
    !   ABI_DEALLOCATE(operks)
    ! endif
     ABI_ALLOCATE(operks,(wan%nkpt,wan%bandf_wan-wan%bandi_wan+1,wan%bandf_wan-wan%bandi_wan+1,wan%nsppol))
     operks = czero

     ! Fill diagonal elements to have DFT Green's function.
     !------------------------------------------------------
     do isppol = 1,wan%nsppol
       do iband1 = 1,wan%bandf_wan-wan%bandi_wan+1
         ibandc = iband1 + wan%bandi_wan - 1
         do ikpt = 1,wan%nkpt
           operks(ikpt,iband1,iband1,isppol) = &
                  & 1d0/(wcurrent-eigen(((ikpt-1)*dtset%mband+ibandc+(isppol-1)*wan%nkpt*dtset%mband))+fermie) ! 1/(w-E(kv))
         end do
       end do
     end do

     ! Compute Green's function in wannier basis in recip space.
     !----------------------------------------------------------
     do ikpt = 1,wan%nkpt
       call compute_oper_ks2wan(wan,operks,operwan,ikpt) !in reciprocal space
     end do
     ABI_DEALLOCATE(operks)


     ! Transform the operwan into a better shape for inversion
     !----------------------------------------------------------
     ABI_ALLOCATE(operwansquare,(wan%nkpt,wan%nsppol,wan%nspinor*wan%size_wan,wan%nspinor*wan%size_wan))
     operwansquare = czero
     do ikpt = 1,wan%nkpt
       do isppol = 1,wan%nsppol
         do ispinor1 = 1,wan%nspinor
           do ispinor2 = 1,wan%nspinor
             index_l = 0 !index_l is set to 0 at the beginning
             do iatom1 = 1,wan%natom_wan
               do il1 = 1,wan%nbl_atom_wan(iatom1)
                 do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
                   index_l = index_l + 1 !the line changes
                   index_c = 1 !counter_c is set to one each time the line changes
                   do iatom2 = 1,wan%natom_wan
                     do il2 = 1,wan%nbl_atom_wan(iatom2)
                       do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                         operwansquare(ikpt,isppol,index_l+wan%size_wan*(ispinor1-1),index_c+wan%size_wan*(ispinor2-1)) = &
                         &        operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2)
                         index_c = index_c + 1
                       end do !im2
                     end do !il2
                   end do !iatom2 (the line changes)
                 end do !im1
               end do !il1
             end do !iatom1
           end do
         end do
       end do
     end do

     ABI_ALLOCATE(operwansquarereal,(wan%nsppol,size(operwansquare,3),size(operwansquare,4)))
     operwansquarereal = czero

     ! Transformation in the real space (T=T'=0) : compute the local quantities
     !-------------------------------------------------------------------------
     do isppol = 1,wan%nsppol
       do index_l = 1,size(operwansquare,3)
         do index_c = 1,size(operwansquare,4)
           do ikpt = 1,wan%nkpt
             operwansquarereal(isppol,index_l,index_c) = operwansquarereal(isppol,index_l,index_c) + &
&                     operwansquare(ikpt,isppol,index_l,index_c)*wan%wtk(ikpt)
           end do
         end do
       end do
     end do
     ABI_DEALLOCATE(operwansquare)

     if (dos .ge. 1) then ! either we compute the DOS (-imaginary part/Pi of the green function in the Wannier basis)
       ! Compute the dos
       !-------------------------------------------------------------------------
       sum = 0
       do isppol = 1,wan%nsppol
         do m1 = 1,2*wan%latom_wan(1)%lcalc(dos)+1
           sum = sum - aimag(operwansquarereal(isppol,shift+m1,shift+m1))
         end do
         if(isppol==1) write(dos_unt,'(F8.3,E15.6)') real(27.2101*wcurrent),sum/(27.2107*3.14159)
         if(isppol==2) write(dos_unt2,'(F8.3,E15.6)') real(27.2101*wcurrent),sum/(27.2107*3.14159)
       end do

     else ! either we compute the F part (the residual part) of the invert of the green functions
       ! Compute the hybridization
       !-------------------------------------------------------------------------

       do isppol = 1,wan%nsppol
         ABI_ALLOCATE(matrix_to_diag,(sizem,sizem))
         matrix_to_diag = operwansquarereal(isppol,shift+1:shift+sizem,shift+1:shift+sizem)
     ! attention a isppol ci dessus
         call xginv(matrix_to_diag,sizem)
         operwansquarereal(isppol,shift+1:shift+sizem,shift+1:shift+sizem) = matrix_to_diag
         ABI_DEALLOCATE(matrix_to_diag)
         mat_writing = ""
         do m1 = 1,sizem
           Fff(isppol) = wcurrent - operwansquarereal(isppol,shift+m1,shift+m1) - energies(m1,isppol)
           write(mat_writing2,'(E15.6)') -aimag(27.2107*Fff(isppol))
           mat_writing = trim(mat_writing)//" "//trim(mat_writing2)
         end do
         if(isppol==1) write(dos_unt,'(F10.3,a)') real(27.2107*wcurrent),trim(mat_writing)
         if(isppol==2) write(dos_unt2,'(F10.3,a)') real(27.2107*wcurrent),trim(mat_writing)
       end do
     end if
     call destroy_operwan(wan,operwan)
     ABI_DATATYPE_DEALLOCATE(operwan)
     ABI_DEALLOCATE(operwansquarereal)
   end do
   close(dos_unt)
   close(dos_unt2)
   ABI_DEALLOCATE(energies)
   ABI_DEALLOCATE(Fff)

 end if ! choice of the 1 plowan_computegreen




 if (plowan_computegreen .eq. 2 .and. me.eq.0 ) then !! Not working ! not tested, not up to date with the code
   !Method 1
   ABI_ALLOCATE(energies,(7,wan%nsppol))
   ABI_ALLOCATE(Ffftable,(7,wan%nsppol))
   Ffftable = czero
   energies = czero
   ! Keep energies for later use
   !----------------------------
   do ikpt = 1,wan%nkpt
     do im1 = 1,7
       do isppol = 1,wan%nsppol
         energies(im1,isppol) = energies(im1,isppol) + real(operwan(ikpt,1,1)%atom(1,1)%matl(im1,im1,isppol,1,1))*wan%wtk(ikpt)
       end do
     end do
   end do
   write(std_out,*) "energies", energies*27.211

   ! Loop over frequency
   !----------------------
   do iw = 1,number_of_frequencies

     ABI_ALLOCATE(operwansquare,(wan%nkpt,wan%nsppol,wan%nspinor*wan%size_wan,wan%nspinor*wan%size_wan))
     operwansquare = czero
     wcurrent = wbase + (iw-1)*wincrease
     Ffftable = czero

     ! create operwansquare
     !----------------------------
     do ikpt = 1,wan%nkpt
       do isppol = 1,wan%nsppol
         do ispinor1 = 1,wan%nspinor
           do ispinor2 = 1,wan%nspinor
             index_l = 0 ! index_l is set to 0 at the beginning
             do iatom1 = 1,wan%natom_wan
               do il1 = 1,wan%nbl_atom_wan(iatom1)
                 do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
                   index_l = index_l + 1 ! the line changes
                   index_c = 1 ! counter_c is set to one each time the line changes
                   do iatom2 = 1,wan%natom_wan
                     do il2 = 1,wan%nbl_atom_wan(iatom2)
                       do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                         operwansquare(ikpt,isppol,index_l+wan%size_wan*(ispinor1-1),index_c+wan%size_wan*(ispinor2-1)) &
          &                       = operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2)
                         index_c = index_c + 1
                       end do !im2
                     end do !il2
                   end do ! iatom2 (the line changes)
                 end do ! im1
               end do ! il1
             end do !iatom1
           end do
         end do
       end do

     ! Create inverse of Green's function
     !-----------------------------------
       do isppol = 1,wan%nsppol
         do im1 = 1,size(operwansquare,3)
           do im2 = 1,size(operwansquare,4)
             if (im1 .eq. im2) then
               operwansquare(ikpt,isppol,im1,im2) = wcurrent-operwansquare(ikpt,isppol,im1,im2)
             else
               operwansquare(ikpt,isppol,im1,im2) = -operwansquare(ikpt,isppol,im1,im2)
             end if
           end do
         end do
       end do


     ! Create Green's function
     !-----------------------------------
       ABI_ALLOCATE(matrix_to_diag,(size(operwansquare,3),size(operwansquare,3)))
       do isppol = 1,wan%nsppol
         matrix_to_diag = czero
         matrix_to_diag = operwansquare(ikpt,isppol,:,:)
         call xginv(matrix_to_diag,size(matrix_to_diag,1))
         operwansquare(ikpt,isppol,:,:) = matrix_to_diag
       end do
       ABI_DEALLOCATE(matrix_to_diag)


!!     if (dos .eq. 0) then
!!       !select f bands
!!       do isppol = 1,wan%nsppol
!!         if (allocated(matrix_to_diag)) ABI_DEALLOCATE(matrix_to_diag)
!!         ABI_ALLOCATE(matrix_to_diag,(7,7))
!!         matrix_to_diag = czero
!!         matrix_to_diag = operwansquare(ikpt,isppol,1:7,1:7)
!!         call xginv(matrix_to_diag,size(matrix_to_diag,1))

!!         operwansquare(ikpt,isppol,1:7,1:7) = matrix_to_diag
!!       end do
!!     end if

!!     do im1 = 1,7
!!       do isppol = 1,wan%nsppol
!!         if (dos .eq. 0) then
!!          Ffftable(im1,isppol) = Ffftable(im1,isppol) + (wcurrent - &
!!  &             operwansquare(ikpt,isppol,im1,im1) - energies(im1,isppol))*wan%wtk(ikpt)
!!         end if
!!         if (dos .eq. 1) then
!!           Ffftable(im1,isppol) = Ffftable(im1,isppol) + operwansquare(ikpt,isppol,im1,im1)*wan%wtk(ikpt)
!!         end if
!!       end do
!!     end do

     end do !!loop ikpt

     ABI_ALLOCATE(operwansquarereal,(wan%nsppol,size(operwansquare,3),size(operwansquare,4)))
     operwansquarereal = czero

     ! Compute local Green's function
     !--------------------------------------------
     do isppol = 1,wan%nsppol
       do index_l = 1,size(operwansquare,3)
         do index_c = 1,size(operwansquare,4)
           do ikpt = 1,wan%nkpt
             operwansquarereal(isppol,index_l,index_c) = operwansquarereal(isppol,index_l,index_c)&
       &              + operwansquare(ikpt,isppol,index_l,index_c)*wan%wtk(ikpt)
           end do
         end do
       end do
     end do

     ! Inverse Local Correlated Green's function
     !--------------------------------------------
     if (dos < 0) then
       !select f bands
       do isppol = 1,wan%nsppol
         ABI_ALLOCATE(matrix_to_diag,(7,7))
         matrix_to_diag = operwansquarereal(isppol,1:7,1:7)
         call xginv(matrix_to_diag,size(matrix_to_diag,1))
         operwansquarereal(isppol,1:7,1:7) = matrix_to_diag
         ABI_DEALLOCATE(matrix_to_diag)
       end do
       write(268,*) 27.2107*real(wcurrent),27.2107*real(operwansquarereal(1,1,1)),27.2107*aimag(operwansquarereal(1,1,1))

     end if

     do im1 = 1,7
       do isppol = 1,wan%nsppol
         ! Compute hybridization
         !--------------------------------------------
         if (dos < 0) then
           Ffftable(im1,isppol) = Ffftable(im1,isppol) + wcurrent - operwansquarereal(isppol,im1,im1) - energies(im1,isppol)
         end if
         if (dos > 0) then
         ! Compute Dos
         !--------------------------------------------
           Ffftable(im1,isppol) = Ffftable(im1,isppol) + operwansquarereal(isppol,im1,im1)
         end if
       end do
     end do
     write(269,*) 27.2107*real(wcurrent),27.2107*real(Ffftable(1,1)),27.2107*aimag(Ffftable(1,1))
     write(2699,*) 27.2107*real(wcurrent),27.2107*real(wcurrent - operwansquarereal(1,1,1) - energies(1,1)),&
&      27.211*real(energies(1,1)),27.211*real(operwansquarereal(1,1,1))

     if (dos < 0) then
       xsum=czero
       do im1 = 1,7
         if (wan%nsppol .eq. 1) then
           write(std_out,*)'hybri',im1,27.2107*wcurrent,27.2107*Ffftable(im1,1)
         else
           do isppol = 1,wan%nsppol
             write(std_out,*)'hybri',im1,isppol,27.2107*wcurrent,27.2107*Ffftable(im1,isppol)
             xsum=xsum+Ffftable(im1,isppol)
           end do
         end if
       end do
       write(270,*)27.2107*real(wcurrent),27.2107*real(xsum),27.2107*aimag(xsum)
     end if

     if (dos > 0) then
       xsum=czero
       do isppol = 1,wan%nsppol
         Ffftable(1,isppol) = Ffftable(1,isppol)+Ffftable(2,isppol)+Ffftable(3,isppol)+Ffftable(4,isppol)&
           &      +Ffftable(5,isppol)+Ffftable(6,isppol)+Ffftable(7,isppol)
         if (wan%nsppol .eq. 2) then
           write(std_out,*)'green',isppol,27.2107*wcurrent,Ffftable(1,isppol)/(27.2107*3.14159)
         else
           write(std_out,*)'green',27.2107*wcurrent,Ffftable(1,isppol)/(27.2107*3.14159)
             xsum=xsum+Ffftable(im1,isppol)
         end if
       write(271,*)27.2107*real(wcurrent),real(xsum)/(27.2107*3.14159),aimag(xsum)/(27.2107*3.14159)
       end do
     end if
     ABI_DEALLOCATE(operwansquare)
     ABI_DEALLOCATE(operwansquarereal)
   end do !loop frequencies w
   ABI_DEALLOCATE(energies)
   ABI_DEALLOCATE(Fff)

   call destroy_operwan(wan,operwan)
   ABI_DATATYPE_DEALLOCATE(operwan)
 end if !! choice of the plowan_computegreen

 if(plowan_computegreen==0) then
   call destroy_operwan(wan,operwan)
   ABI_DATATYPE_DEALLOCATE(operwan)
 end if !! choice of the plowan_computegreen

 !deallocate temporary cwaveprj/cprj storage
 call pawcprj_free(cwaveprj)
 ABI_DATATYPE_DEALLOCATE(cwaveprj)


end subroutine compute_coeff_plowannier
!!***

!!****f* m_plowannier/print_plowannier
!! NAME
!!  print_plowannier
!!
!! FUNCTION
!!  print the wannier weight (psichi) on a forlb.ovlp file
!!
!! INPUTS
!! dtset%typat,wan
!!
!! OUTPUT
!!
!! PARENTS
!! outscfcv
!!
!! CHILDREN
!! open_file, wrtout
!! SOURCE



 subroutine print_plowannier(wan)

 use m_abicore
 use m_io_tools,  only : open_file
 use m_specialmsg, only : wrtout
 implicit none

 !Arguments-------------------------
 type(plowannier_type),intent(in) :: wan
 !Local variables-------------------
 character(len=500) :: msg
 integer :: unt,iatom,spin,ikpt,iband,ibandc,il,ispinor,im

 !Creation of the data.plowann file
 if (open_file('data.plowann',msg,newunit=unt,form='formatted',status='replace') /= 0) then
  MSG_ERROR(msg)
 end if
 rewind(unt)

 write(msg,'(2a)') ch10,' Print the psichi coefficients in data.plowann'
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 
 !Header of the file data.plowann
 write(unt,'(a22,i2)')"Total number of atom =", wan%natom_wan
 write(unt,*)"List of atoms", wan%iatom_wan(:)
 write(unt,'(a7,2i4)')"Bands =",wan%bandi_wan,wan%bandf_wan
 write(unt,'(a26,i2)')"Total number of orbitals =",sum(wan%nbl_atom_wan(:))
 do iatom=1,wan%natom_wan
   write(unt,'(a17,i2,a3,4i2)')"Orbitals for atom",wan%iatom_wan(iatom)," = ",wan%latom_wan(iatom)%lcalc(:)
 enddo
 write(unt,'(a16,i2)')"Number of spin =",wan%nsppol
 write(unt,'(a19,i4)')"Number of k-points=",wan%nkpt
 do ikpt=1,wan%nkpt
   write(unt,'(a,2x,i4)')"ikpt =",ikpt
    do spin=1,wan%nsppol
      do ispinor=1,wan%nspinor
        do iband=wan%bandi_wan,wan%bandf_wan
          ibandc=iband-wan%bandi_wan+1
          write(unt,'(2x,a,2x,i2,2x,i2)')"iband =",iband
          do iatom=1,wan%natom_wan
            do il=1,wan%nbl_atom_wan(iatom)
              do im=1,2*wan%latom_wan(iatom)%lcalc(il)+1
                write(unt,'(8x,3i3,2x,2f23.15)')iatom,wan%latom_wan(iatom)%lcalc(il),im,&
                &real(wan%psichi(ikpt,ibandc,iatom)%atom(il)%matl(im,spin,ispinor))&
                &,aimag(wan%psichi(ikpt,ibandc,iatom)%atom(il)%matl(im,spin,ispinor))
              enddo!m
            enddo!l
          enddo!atom
        enddo!band
      enddo!spinor
    enddo!spin
  enddo!k-point
 close(unt)
 end subroutine print_plowannier
!!***


!!****f* m_plowannier/get_plowannier
!! NAME
!!  get_plowannier
!!
!! FUNCTION
!!  get the psichies (Wannier weights) from a data.plowann file
!!
!! INPUTS
!! wan
!!
!! OUTPUT
!! wan 
!! 
!! PARENTS
!! outscfcv,cchi0, cchi0q0, prep_calc_ucrpa
!!
!! CHILDREN
!! open_file,
!! SOURCE

 subroutine get_plowannier(wan_in,wan_out,dtset)

 use m_abicore
 use defs_abitypes
 use m_io_tools,  only : open_file
 use m_specialmsg, only : wrtout
 implicit none

 !Arguments-------------------------
 type(plowannier_type),intent(inout) :: wan_in
 type(plowannier_type),intent(inout) :: wan_out
 type(dataset_type),intent(in) :: dtset
 !Local variables-------------------
 character(len=500) :: msg
 integer :: unt,iatom,spin,ikpt,iband,ibandc,il,ispinor,im,dummy,natom,bandi,bandf,nbl,nspin,nkpt
 real(dp) ::xx,yy

 !Opening of the data.plowann file
 if (open_file('data.plowann',msg,newunit=unt,form='formatted',status='old') /= 0) then
  MSG_ERROR(msg)
 end if
 rewind(unt)
 

 !Reading of the header of data.plowann
 read(unt,'(a22,i2)') msg, natom
 read(unt,*)
 read(unt,'(a7,2i4)') msg, bandi,bandf
 read(unt,'(a26,i2)') msg, nbl
 do iatom=1,wan_in%natom_wan
   read(unt,*)
 enddo
 read(unt,'(a16,i2)') msg, nspin
 read(unt,'(a19,i4)') msg, nkpt

 !Testing the header
 if (natom /= wan_in%natom_wan .OR.&
& nbl/= sum(wan_in%nbl_atom_wan(:)) .OR. nspin /= wan_in%nsppol .OR. nkpt/=wan_in%nkpt ) then
   write(msg,'(a,3i3)')"Not the same atoms or bands in both datasets",natom,bandi,bandf
   MSG_ERROR(msg)
 endif

 call init_plowannier(bandf,bandi,dtset%plowan_compute,&
     &dtset%plowan_iatom,dtset%plowan_it,dtset%plowan_lcalc,dtset%plowan_natom,&
     &dtset%plowan_nbl,dtset%plowan_nt,dtset%plowan_projcalc,dtset%acell_orig,&
     &dtset%kptns,dtset%nimage,dtset%nkpt,dtset%nspinor,dtset%nsppol,dtset%wtk,wan_out)

 call destroy_plowannier(wan_in)

 write(msg,'(a)')"Reading of the Wannier weights from data.plowann"
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 !Reading of the psichis 
 do ikpt=1,wan_out%nkpt
   read(unt,*)
    do spin=1,wan_out%nsppol
      do ispinor=1,wan_out%nspinor
        do iband=wan_out%bandi_wan,wan_out%bandf_wan
          ibandc=iband-wan_out%bandi_wan+1
          read(unt,*)
          do iatom=1,wan_out%natom_wan
            do il=1,wan_out%nbl_atom_wan(iatom)
              do im=1,2*wan_out%latom_wan(iatom)%lcalc(il)+1
                read(unt,'(8x,3i3,2x,2f23.15)')dummy,dummy,dummy,xx,yy
                wan_out%psichi(ikpt,ibandc,iatom)%atom(il)%matl(im,spin,ispinor)=cmplx(xx,yy)
              enddo!m
            enddo!l
          enddo!atom
        enddo!band
      enddo!spinor
    enddo!spin
  enddo!k-point
 close(unt)
end subroutine get_plowannier
!!***


!!****f* m_plowannier/fullbz_plowannier
!! NAME
!!  fullbz_plowannier
!!
!! FUNCTION
!!  Reconstruct the pischis on the full BZ 
!!
!! INPUTS
!! dtset,kmesh,cryst,wanibz
!!
!! OUTPUT
!! wanbz 
!! 
!! PARENTS
!! screening_driver
!!
!! CHILDREN
!! 
!! SOURCE

 subroutine fullbz_plowannier(dtset,kmesh,cryst,pawang,wanibz,wanbz)
   
   use m_abicore
   use m_specialmsg, only : wrtout
   use defs_abitypes
   use m_bz_mesh, only : kmesh_t, get_BZ_item
   use m_crystal, only : crystal_t
   use m_pawang, only  : pawang_type
   implicit none

!Arguments-------------------------
   type(plowannier_type),intent(inout) :: wanibz
   type(plowannier_type),intent(out) :: wanbz
   type(dataset_type),intent(in) :: dtset
   type(kmesh_t),intent(in) :: kmesh
   type(crystal_t),intent(in) :: cryst
   type(pawang_type),intent(in) :: pawang
!Local variables----------------------
   character(len=500) :: msg
   integer :: sym,iatom,spin,ik_bz,iband,ibandc,il,ispinor,im,ik_ibz,isym,itim
   integer ::  at_indx,indx, iat,m1,m2,l
   real(dp) :: kbz(3),wtk(kmesh%nbz)
!*****************************************************************************************
  
   wtk=one
   sym=kmesh%nbz/kmesh%nibz
   call init_plowannier(wanibz%bandf_wan,wanibz%bandi_wan,dtset%plowan_compute,&
     &dtset%plowan_iatom,dtset%plowan_it,dtset%plowan_lcalc,dtset%plowan_natom,&
     &dtset%plowan_nbl,dtset%plowan_nt,dtset%plowan_projcalc,dtset%acell_orig,&
     &kmesh%bz,dtset%nimage,kmesh%nbz,dtset%nspinor,dtset%nsppol,wtk,wanbz)
  
   write(msg,'(a)')" Reconstruction of the full Brillouin Zone using data.plowann in the IBZ"
   call wrtout(std_out,msg,'COLL');call wrtout(ab_out,msg,'COLL')
   if (cryst%nsym==1) then
     do ik_bz=1,kmesh%nbz
       do iband=wanbz%bandi_wan,wanbz%bandf_wan
         ibandc=iband-wanbz%bandi_wan+1
         do iatom=1,wanbz%natom_wan
           do il=1,wanbz%nbl_atom_wan(iatom)
             do im=1,2*wanbz%latom_wan(iatom)%lcalc(il)+1
               do spin=1,wanbz%nsppol
                 do ispinor=1,wanbz%nspinor
                   if (kmesh%tabi(ik_bz)==1) then 
                     wanbz%psichi(ik_bz,ibandc,iatom)%atom(il)%matl(im,spin,ispinor)=&
                       &wanibz%psichi(kmesh%tab(ik_bz),ibandc,iatom)%atom(il)%matl(im,spin,ispinor)
                   else if (kmesh%tabi(ik_bz)==-1) then
                     wanbz%psichi(ik_bz,ibandc,iatom)%atom(il)%matl(im,spin,ispinor)=&
                       &conjg(wanibz%psichi(kmesh%tab(ik_bz),ibandc,iatom)%atom(il)%matl(im,spin,ispinor))
                   endif
                 enddo
               enddo
             enddo
           enddo
         enddo
       enddo
     enddo
   else if (cryst%nsym>1) then
    do ik_bz=1,kmesh%nbz
      call get_BZ_item(kmesh,ik_bz,kbz,ik_ibz,isym,itim)
      do iatom=1,wanbz%natom_wan
        indx=cryst%indsym(4,isym,wanibz%iatom_wan(iatom))
!Link beetween full list and wan list of atom
        do iat=1,wanbz%natom_wan
          if (indx==wanbz%iatom_wan(iat))then
            at_indx=iat
          end if
        end do
!
        do spin=1,wanibz%nsppol
          do ispinor=1,wanibz%nspinor
            do il=1,wanibz%nbl_atom_wan(iatom)
              l=wanibz%latom_wan(iatom)%lcalc(il)
              do m1=1,2*l+1
                do m2=1,2*l+1
                  do iband=wanibz%bandi_wan,wanibz%bandf_wan
                    ibandc=iband-wanibz%bandi_wan+1
                    wanbz%psichi(ik_bz,ibandc,iatom)%atom(il)%matl(m1,spin,ispinor)=&
                      &wanbz%psichi(ik_bz,ibandc,iatom)%atom(il)%matl(m1,spin,ispinor)+&
                      &wanibz%psichi(ik_ibz,ibandc,at_indx)%atom(il)%matl(m2,spin,ispinor)&
                      &*pawang%zarot(m2,m1,l+1,isym)
                  end do
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  end if
  call destroy_plowannier(wanibz)
end subroutine fullbz_plowannier
!!***

!!****f* m_plowannier/destroy_plowannier
!! NAME
!!  destroy_plowannier
!!
!! FUNCTION
!!  deallocate variables
!!
!! INPUTS
!!  wan
!!
!! OUTPUT
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!
!! SOURCE


 subroutine destroy_plowannier(wan)

!Arguments-------------------------------------
 type(plowannier_type), intent(inout) :: wan
!Local variables-------------------------------
 integer :: iatom,ikpt,iband,il

 do iatom=1,wan%natom_wan
   ABI_DEALLOCATE(wan%latom_wan(iatom)%lcalc)
   ABI_DEALLOCATE(wan%projector_wan(iatom)%lproj)
   ABI_DEALLOCATE(wan%nposition(iatom)%pos)
 enddo
 do iatom = 1,wan%natom_wan
   do il = 1,wan%nbl_atom_wan(iatom)
     ABI_DEALLOCATE(wan%psichi(1,1,iatom)%atom(il)%ph0phiint)
   end do
 end do
 do ikpt = 1,wan%nkpt
   do iband = wan%bandi_wan,wan%bandf_wan
     do iatom = 1,wan%natom_wan
       do il = 1,wan%nbl_atom_wan(iatom)
        ABI_DEALLOCATE(wan%psichi(ikpt,iband-wan%bandi_wan+1,iatom)%atom(il)%matl)
       end do
       ABI_DATATYPE_DEALLOCATE(wan%psichi(ikpt,iband-wan%bandi_wan+1,iatom)%atom)
     end do
   end do
 end do

 if (allocated(wan%kpt)) then
   ABI_DEALLOCATE(wan%kpt)
 end if
 if (allocated(wan%iatom_wan)) then
   ABI_DEALLOCATE(wan%iatom_wan)
 end if
 if (allocated(wan%nbl_atom_wan)) then
   ABI_DEALLOCATE(wan%nbl_atom_wan)
 end if
 if (allocated(wan%latom_wan)) then
   ABI_DATATYPE_DEALLOCATE(wan%latom_wan)
 end if
 if (allocated(wan%nbproj_atom_wan)) then
   ABI_DEALLOCATE(wan%nbproj_atom_wan)
 end if
 if (allocated(wan%projector_wan)) then
   ABI_DATATYPE_DEALLOCATE(wan%projector_wan)
 end if
 if (allocated(wan%position)) then
   ABI_DEALLOCATE(wan%position)
 end if
 if (allocated(wan%wtk)) then
   ABI_DEALLOCATE(wan%wtk)
 end if
 if (allocated(wan%acell)) then
   ABI_DEALLOCATE(wan%acell)
 end if


 if (allocated(wan%nposition)) then
   ABI_DATATYPE_DEALLOCATE(wan%nposition)
 end if
 if (allocated(wan%psichi)) then
   ABI_DATATYPE_DEALLOCATE(wan%psichi)
 end if


 end subroutine destroy_plowannier
!!***



!!****f* m_plowannier/initialize_operwan
!! NAME
!!  initialize_operwan
!!
!! FUNCTION
!!  initialize operwan
!!
!! INPUTS
!!  wan
!!
!! OUTPUT
!!  operwan
!!
!! PARENTS
!!      m_plowannier
!!
!! CHILDREN
!!
!! SOURCE

 subroutine initialize_operwan(wan,operwan)

   !Arguments----------------------------------
   type(plowannier_type), intent(in) :: wan
   type(operwan_type), intent(inout) :: operwan(wan%nkpt,wan%natom_wan,wan%natom_wan)

   !Local variables----------------------------
   integer :: ikpt,iatom1,iatom2,il1,il2,n1,n2

   do ikpt = 1,wan%nkpt
     do iatom1 = 1,wan%natom_wan
       do iatom2 = 1,wan%natom_wan
         ABI_DATATYPE_ALLOCATE(operwan(ikpt,iatom1,iatom2)%atom,(wan%nbl_atom_wan(iatom1),wan%nbl_atom_wan(iatom2)))
         do il1 = 1,wan%nbl_atom_wan(iatom1)
           do il2 = 1,wan%nbl_atom_wan(iatom2)
             n1=2*wan%latom_wan(iatom1)%lcalc(il1)+1
             n2=2*wan%latom_wan(iatom2)%lcalc(il2)+1
   ABI_ALLOCATE(operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl,(n1,n2,wan%nsppol,wan%nspinor,wan%nspinor))
             operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl = zero
           end do
         end do
       end do
     end do
   end do

 end subroutine initialize_operwan
!!***


!!****f* m_plowannier/destroy_operwan
!! NAME
!!  destroy_operwan
!!
!! FUNCTION
!!  destroy operwan
!!
!! INPUTS
!!  wan
!!
!! OUTPUT
!!  operwan
!!
!! PARENTS
!!      m_plowannier
!!
!! CHILDREN
!!
!! SOURCE

 subroutine destroy_operwan(wan,operwan)

   !Arguments----------------------------------
   type(plowannier_type), intent(in) :: wan
   type(operwan_type), intent(inout) :: operwan(wan%nkpt,wan%natom_wan,wan%natom_wan)

   !Local variables----------------------------
   integer :: ikpt,iatom1,iatom2,il1,il2


   do ikpt = 1,wan%nkpt
     do iatom1 = 1,wan%natom_wan
       do iatom2 = 1,wan%natom_wan
         do il1 = 1,wan%nbl_atom_wan(iatom1)
           do il2 = 1,wan%nbl_atom_wan(iatom2)
             ABI_DEALLOCATE(operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl)
           end do
         end do
         ABI_DATATYPE_DEALLOCATE(operwan(ikpt,iatom1,iatom2)%atom)
       end do
     end do
   end do
 end subroutine destroy_operwan
!!***

!!****f* m_plowannier/zero_operwan
!! NAME
!!  zero_operwan
!!
!! FUNCTION
!!  zero operwan
!!
!! INPUTS
!!  wan
!!
!! OUTPUT
!!  operwan
!!
!! PARENTS
!!      m_plowannier
!!
!! CHILDREN
!!
!! SOURCE

 subroutine zero_operwan(wan,operwan)

   implicit none

   !Arguments----------------------------------
   type(plowannier_type), intent(in) :: wan
   type(operwan_type), intent(inout) :: operwan(wan%nkpt,wan%natom_wan,wan%natom_wan)

   !Local variables----------------------------
   integer :: ikpt,  iatom1, iatom2, il1, il2, isppol, ispinor1, ispinor2,  im1, im2


   do ikpt = 1,wan%nkpt
     do isppol = 1,wan%nsppol
       do iatom1 = 1,wan%natom_wan
         do iatom2 = 1,wan%natom_wan
           do il1 = 1,wan%nbl_atom_wan(iatom1)
             do il2 = 1,wan%nbl_atom_wan(iatom2)
               do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
                 do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                   do ispinor1 = 1,wan%nspinor
                     do ispinor2 = 1,wan%nspinor
                       operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2)=czero
                     end do
                   end do
                 end do
               end do
             end do
           end do
         end do
       end do
     end do
   end do

 end subroutine zero_operwan
!!***

!!****f* m_plowannier/compute_oper_ks2wan
!! NAME
!!  compute_oper_ks2wan
!!
!! FUNCTION
!!  transform ks operator into wan one
!!
!! INPUTS
!!  wan,operks,option
!!
!! OUTPUT
!!  if option = ikpt, gives the wan operator in reciprocal space (for each k)
!!
!! PARENTS
!!      m_plowannier
!!
!! CHILDREN
!!
!! SOURCE

 subroutine compute_oper_ks2wan(wan,operks,operwan,option)

   !Arguments--------------------------
   type(plowannier_type), intent(in) :: wan
   type(operwan_type), intent(inout) :: operwan(:,:,:)
   complex(dpc), intent(in) :: operks(:,:,:,:)
   integer, intent(in) :: option

   !Local variables--------------------
   integer :: iatom1, iatom2, il1, il2, isppol, ispinor1, ispinor2, iband1, iband2, im1, im2

   ! ----------------------------------
   !Transformation KS2WAN
   ! ----------------------------------


   !!operation on reciprocal space
   do iatom1 = 1,wan%natom_wan
     do iatom2 = 1,wan%natom_wan
       do isppol = 1,wan%nsppol
         do il1 = 1,wan%nbl_atom_wan(iatom1)
           do il2 = 1,wan%nbl_atom_wan(iatom2)
             do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
               do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                 do ispinor1 = 1,wan%nspinor
                   do ispinor2 = 1,wan%nspinor
                     !!sum over the bands
                     do iband1 = 1,wan%bandf_wan-wan%bandi_wan+1
                       do iband2 = 1,wan%bandf_wan-wan%bandi_wan+1
                         operwan(option,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2)=&
                           operwan(option,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2)&
                          + conjg(wan%psichi(option,iband2,iatom2)%atom(il2)%matl(im2,isppol,ispinor2))&
                   *operks(option,iband1,iband2,isppol)*wan%psichi(option,iband1,iatom1)%atom(il1)%matl(im1,isppol,ispinor1)
                       end do
                     end do
                    ! write(6,*) "operwan",im1,im2,operwan(option,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2)
                   end do
                 end do
               end do
             end do
           end do
         end do
       end do
     end do
   end do

end subroutine compute_oper_ks2wan
!!***



!!****f* m_plowannier/normalization_plowannier
!! NAME
!!  normalization_plowannier
!!
!! FUNCTION
!!  Use compute_oper_ks2wan to calculate overlap and do the normalization for the wan%psichi coefficients
!!
!! INPUTS
!!  wan, opt (=0 normalize k-point by k-point; =1 normalize the sum over k)
!!
!! OUTPUT
!!  wan itself is modified
!!
!! PARENTS
!!      m_plowannier
!!
!! CHILDREN
!!
!! SOURCE


subroutine normalization_plowannier(wan,opt)

  use m_matrix, only : invsqrt_matrix

!Arguments------------------
  type(plowannier_type), intent(inout) :: wan
  integer, intent(in) :: opt
!Local----------------------
  complex(dpc), allocatable :: operks(:,:,:,:)
  type(operwan_type), allocatable :: operwan(:,:,:)
  complex(dpc), allocatable :: operwansquare(:,:,:,:)
  complex(dpc), allocatable :: tmp_operwansquare(:,:)
  integer :: ikpt, iband, iband1, iband2, isppol,  ispinor1, ispinor2, iatom1,nb_zeros_tot
  integer :: iatom2, il1, il2, im1, im2, index_c, index_l, n1,n2,n3, nkpt,nb_of_zeros
  type(orbital_type), allocatable :: psichinormalized(:,:,:)
  !character(len = 50) :: mat_writing2
  !character(len = 5000) :: mat_writing
  character(len = 500) :: message

  !Initialize nkpt (wan%nkpt if opt=0, 1 if opt=1)
  if (opt==1) then 
    nkpt=1
  else 
    nkpt=wan%nkpt
  end if

  !First, creation of the ks identity operator
  ABI_ALLOCATE(operks,(wan%nkpt,wan%bandf_wan-wan%bandi_wan+1,wan%bandf_wan-wan%bandi_wan+1,wan%nsppol))
  operks = czero
  do iband1 = 1,wan%bandf_wan-wan%bandi_wan+1
    do iband2 = 1,wan%bandf_wan-wan%bandi_wan+1
      if (iband1.eq.iband2) then
        do ikpt = 1,wan%nkpt
          do isppol= 1,wan%nsppol
            operks(ikpt,iband1,iband2,isppol) = cone
          end do
        end do
      end if
    end do
  end do


  !Allocation of operwan
  ABI_DATATYPE_ALLOCATE(operwan,(wan%nkpt,wan%natom_wan,wan%natom_wan))
  call initialize_operwan(wan,operwan)



  !Computation of the overlap
  do ikpt = 1,wan%nkpt
    call compute_oper_ks2wan(wan,operks,operwan,ikpt)
  end do

  

  !transform the operwan into an inversible matrix
  !!operwansquare is the overlap square matrix (wan%size_wan * wan%size_wan)
  ABI_ALLOCATE(operwansquare,(wan%nkpt,wan%nsppol,wan%nspinor*wan%size_wan,wan%nspinor*wan%size_wan))

  operwansquare = czero

  n1=size(wan%psichi,1)
  n2=size(wan%psichi,2)
  n3=size(wan%psichi,3)
  ABI_DATATYPE_ALLOCATE(psichinormalized,(n1,n2,n3))
  call allocate_orbital(wan%psichi,psichinormalized,n1,n2,n3)
  call copy_orbital(wan%psichi,psichinormalized,n1,n2,n3)

  do isppol = 1,wan%nsppol
    do ispinor1 = 1,wan%nspinor
      do ispinor2 = 1,wan%nspinor
        index_l = 0 ! index_l is set to 0 at the beginning
        do iatom1 = 1,wan%natom_wan
          do il1 = 1,wan%nbl_atom_wan(iatom1)
            do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
              index_l = index_l + 1 ! the line changes
              index_c = 1 ! counter_c is set to one each time the line changes
              do iatom2 = 1,wan%natom_wan
                do il2 = 1,wan%nbl_atom_wan(iatom2)
                  do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                    do ikpt = 1,wan%nkpt
                      if (opt==1) then ! operwansquare is a sum of operwan over k points     
                        operwansquare(1,isppol,index_l+wan%size_wan*(ispinor1-1),index_c+wan%size_wan*(ispinor2-1))  &
&                        =operwansquare(1,isppol,index_l+wan%size_wan*(ispinor1-1),index_c+wan%size_wan*(ispinor2-1))+ &
&                        (operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2)*wan%wtk(ikpt))
                      else !operwansquare is a matrix with a k-point dimension
                        operwansquare(ikpt,isppol,index_l+wan%size_wan*(ispinor1-1),index_c+wan%size_wan*(ispinor2-1))  &
                          &= operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2)
                      end if
                    end do
                    index_c = index_c + 1
                  end do !im2
                end do !il2
              end do ! iatom2 (the line changes)
            end do ! im1
          end do ! il1
        end do !iatom1
      end do
    end do
  end do




  !!Write the overlap matrix for ikpt = 1 in a nice shape
!    do isppol = 1,wan%nsppol
!      do il1 = 1,size(operwansquare,3) !! dummy variable without any meaning
!        write(mat_writing,'(a,i0,i0)') 'Overlap matrix before orthonormalization 1 ',isppol,il1
!        do il2 = 1,size(operwansquare,4)
!         write(mat_writing2,'(F10.6)') real(operwansquare(1,isppol,il1,il2))
!          mat_writing = trim(mat_writing)//trim(mat_writing2)
!        end do
!        print*,trim(mat_writing)
!      end do
!    end do



  !take the square root inverse of operwansquare for normalization purposes
  nb_zeros_tot=0
  ABI_ALLOCATE(tmp_operwansquare,(wan%nspinor*wan%size_wan,wan%nspinor*wan%size_wan))
  do isppol = 1,wan%nsppol
    do ikpt = 1,nkpt
      write(std_out,*)"ikpt = ", ikpt
      tmp_operwansquare(:,:)=operwansquare(ikpt,isppol,:,:)
      call invsqrt_matrix(tmp_operwansquare,wan%nspinor*wan%size_wan,nb_of_zeros)
      operwansquare(ikpt,isppol,:,:)=tmp_operwansquare(:,:)
      nb_zeros_tot=nb_zeros_tot+nb_of_zeros
    end do
  end do
  ABI_DEALLOCATE(tmp_operwansquare)

  do ikpt = 1,wan%nkpt
    do iband = 1,wan%bandf_wan-wan%bandi_wan+1
      do iatom1 = 1,wan%natom_wan
        do il1 = 1,wan%nbl_atom_wan(iatom1)
          psichinormalized(ikpt,iband,iatom1)%atom(il1)%matl = czero
        end do
      end do
    end do
  end do



  ! compute the new psichi normalized
  do isppol = 1,wan%nsppol
    do ispinor1 = 1,wan%nspinor
      do iband = 1,wan%bandf_wan-wan%bandi_wan+1
        index_l = 0
        do iatom1 = 1,wan%natom_wan
          do il1 = 1,wan%nbl_atom_wan(iatom1)
            do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
              ! sum
              do ispinor2 = 1,wan%nspinor
                index_l = index_l + 1 ! the line changes
                index_c = 1 ! when the line changes, index_c is set to 1
                do iatom2 = 1,wan%natom_wan
                  do il2 = 1,wan%nbl_atom_wan(iatom2)
                    do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                      do ikpt = 1,wan%nkpt
                        if (opt==1)then ! all the psichi are normalized with the same operwansquare
                          psichinormalized(ikpt,iband,iatom1)%atom(il1)%matl(im1,isppol,ispinor1) =&
&                           psichinormalized(ikpt,iband,iatom1)%atom(il1)%matl(im1,isppol,ispinor1) +&
&                           wan%psichi(ikpt,iband,iatom2)%atom(il2)%matl(im2,isppol,ispinor2)*&
&                           operwansquare(1,isppol,index_l+wan%size_wan*(ispinor1-1),index_c+&
&                           wan%size_wan*(ispinor2-1))                          
                        else ! each psichi is normalized with his own operwansquare
                          psichinormalized(ikpt,iband,iatom1)%atom(il1)%matl(im1,isppol,ispinor1) =&
&                           psichinormalized(ikpt,iband,iatom1)%atom(il1)%matl(im1,isppol,ispinor1) +&
&                           wan%psichi(ikpt,iband,iatom2)%atom(il2)%matl(im2,isppol,ispinor2)*&
&                           operwansquare(ikpt,isppol,index_l+wan%size_wan*(ispinor1-1),index_c+&
&                           wan%size_wan*(ispinor2-1))
                        end if
                      end do
                      index_c = index_c + 1
                    end do !im2
                  end do !il2
                end do !iatom2
              end do ! ispinor2
            end do !im1
          end do ! il1
        end do ! iatom1
      end do ! iband
    end do ! ispinor1
  end do !isppol
  ! copy the new psichi normalized
  call copy_orbital(psichinormalized,wan%psichi,n1,n2,n3)
  call destroy_orbital(psichinormalized,n1,n2,n3)
  ABI_DATATYPE_DEALLOCATE(psichinormalized)

  call destroy_operwan(wan,operwan)
  ABI_DATATYPE_DEALLOCATE(operwan)








!!
!!  !-------------------------------------------------------------
!!  !check if the new norm is one

  ABI_DATATYPE_ALLOCATE(operwan,(wan%nkpt,wan%natom_wan,wan%natom_wan))
  call initialize_operwan(wan,operwan)
  do ikpt = 1,wan%nkpt
    call compute_oper_ks2wan(wan,operks,operwan,ikpt)
  end do


  do isppol = 1,wan%nsppol
    do ikpt = 1,wan%nkpt
      do iatom1 = 1,wan%natom_wan
        do iatom2 = 1,wan%natom_wan
          do il1 = 1,wan%nbl_atom_wan(iatom1)
            do il2 = 1,wan%nbl_atom_wan(iatom2)
              do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
                do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                  do ispinor1 = 1,wan%nspinor
                    do ispinor2 = 1,wan%nspinor
                      if (opt==0 .and. nb_zeros_tot==0) then
                        if (iatom1.eq.iatom2 .and. il1.eq.il2 .and. im1.eq.im2 .and. ispinor1.eq.ispinor2) then
                          if (abs(cmplx(1.0,0.0,dpc)-&
                            &operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%&
                            &matl(im1,im2,isppol,ispinor1,ispinor2)) > 1d-8) then
                            write(message,'(a,i0,a,F18.11)') 'Normalization error for ikpt =',ikpt,&
                              &' on diag, value = ',&
                              &abs(operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2))
                            MSG_ERROR(message)
                          end if
                        else
                          if (abs(operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2)) > 1d-8) then
                            write(message,'(a,i0,a,F10.3)') 'Normalization error for ikpt =',ikpt,&
                              &' not on diag, value = ',&
                              &abs(operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2))
                            MSG_ERROR(message)
                          end if
                        end if
                      end if
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
  if (opt==0 .and. nb_zeros_tot/=0) then
    write(message,'(a,i2,a)')"The matrix inversion detects ",nb_zeros_tot,&
    " zero(s) on the diagonals. Take results with caution or modify nkpt and/or bands for plowan"
    MSG_COMMENT(message)
  end if

  !!Uncomment to print the overlap matrix in the log file (for ikpt = 1)
!  do isppol = 1,wan%nsppol
!    do ispinor1 = 1,wan%nspinor
!      do ispinor2 = 1,wan%nspinor
!        index_l = 0 ! index_l is set to 0 at the beginning
!        do iatom1 = 1,wan%natom_wan
!          do il1 = 1,wan%nbl_atom_wan(iatom1)
!            do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
!              index_l = index_l + 1 ! the line changes
!              index_c = 1 ! counter_c is set to one each time the line changes
!              do iatom2 = 1,wan%natom_wan
!                do il2 = 1,wan%nbl_atom_wan(iatom2)
!                  do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
!                    do ikpt = 1,wan%nkpt
!                      operwansquare(ikpt,isppol,index_l+wan%size_wan*(ispinor1-1),&
!            &index_c+wan%size_wan*(ispinor2-1)) = operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2)
!                    end do
!                    index_c = index_c + 1
!                  end do !im2
!                end do !il2
!              end do ! iatom2 (the line changes)
!            end do ! im1
!          end do ! il1
!        end do !iatom1
!      end do
!    end do
!  end do

  !!Print the overlap matrix in a nice shape
!  do isppol = 1,wan%nsppol
!    do il1 = 1,size(operwansquare,3) !! dummy variable without any meaning
!      write(mat_writing,'(a,i0,i0)') 'Overlap matrix after orthonormalization ',isppol,il1
!      do il2 = 1,size(operwansquare,4)
!        write(mat_writing2,'(F10.6)') real(operwansquare(9,isppol,il1,il2))
!        mat_writing = trim(mat_writing)//trim(mat_writing2)
!      end do
!      print*,trim(mat_writing)
!    end do
!  end do




!!  !----------------------------------------------------------------
  ABI_DEALLOCATE(operwansquare)
  ABI_DEALLOCATE(operks)
  call destroy_operwan(wan,operwan)
  ABI_DATATYPE_DEALLOCATE(operwan)


end subroutine normalization_plowannier

!!***



!!****f* m_plowannier/print_operwan
!! NAME
!!  print_operwan
!!
!! FUNCTION
!!  Print the Wannier operator (real space) in a latex file
!!
!! INPUTS
!!  wan, operwan, name
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      m_plowannier
!!
!! CHILDREN
!!
!! SOURCE


subroutine print_operwan(wan,operwan,name,convert)

!Arguments----------------------------------
  type(operwan_type),intent(in) :: operwan(:,:,:)
  type(plowannier_type), intent(in) :: wan
  character(len=*), intent(in) :: name
  real(dp), intent(in) :: convert

!Local variables----------------------------
  integer :: iatom1,iatom2,pos1,pos2,il1,il2,im1,im2,isppol,ikpt,unt
  real(dp) :: sum
  character(len = 500) :: str1,str2,msg

  if (open_file(name, msg, newunit=unt) /= 0) then
    MSG_ERROR(msg)
  end if

  write(unt,'(a)') '\documentclass[11pt,a4paper,landscape]{article}'
 write(unt,'(a)') '\usepackage[T1]{fontenc}'
  write(unt,'(a)') '\usepackage{geometry,tabularx,graphicx}'
  write(unt,'(a)') '\geometry{left=0.5cm,right=0.5cm}'

  write(unt,'(a)') '\begin{document}'
  write(unt,'(a)') '\noindent'

!  write(unt,'(a,i0,a,F7.3,a)') "% ",wan%natom_wan," atom in a ",wan%acell(1)," cell"
!  write(unt,'(a)') "% atom isppol proj"


do isppol = 1,wan%nsppol
  write(unt,'(a)') '\begin{figure}'
  write(unt,'(a)') '\resizebox{\linewidth}{!}{%'
  write(unt,'(a)') '$ \left('


  write(str1,'(a)') '\begin{array}{'
  do iatom1 = 1,wan%natom_wan
    do pos1 = 1,size(wan%nposition(iatom1)%pos,1)
      do il1 = 1,wan%nbl_atom_wan(iatom1)
        do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
          str1 = trim(str1)//"c"
        end do
        if (iatom1 .ne. wan%natom_wan .or. il1 .ne. wan%nbl_atom_wan(iatom1) .or. pos1 .ne. size(wan%nposition(iatom1)%pos,1)) then
          str1 = trim(str1)//'|'
        end if
      end do
    end do
  end do
  str1 = trim(str1)//'}'
  write(unt,'(a)') str1


  do iatom1 = 1,wan%natom_wan
    do pos1 = 1,size(wan%nposition(iatom1)%pos,1)
      do il1 = 1,wan%nbl_atom_wan(iatom1)
        do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
          write(str1,'(a)') ""
          do iatom2 = 1,wan%natom_wan
            do pos2 = 1,size(wan%nposition(iatom2)%pos,1)
              do il2 = 1,wan%nbl_atom_wan(iatom2)
                do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                  sum = 0
                  do ikpt = 1,wan%nkpt
                    sum = sum + convert*real(operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,1,1)*wan%wtk(ikpt)*&
        &      exp(cmplx(0.0,1.0)*two_pi*(wan%kpt(1,ikpt)*(  &
        &           wan%nposition(iatom1)%pos(pos1,1)-wan%nposition(iatom2)%pos(pos2,1))+ &
        &           wan%kpt(2,ikpt)*(wan%nposition(iatom1)%pos(pos1,2)-wan%nposition(iatom2)%pos(pos2,2))+ &
        &           wan%kpt(3,ikpt)*(wan%nposition(iatom1)%pos(pos1,3)-wan%nposition(iatom2)%pos(pos2,3)))))
                  end do
                  write(str2,'(F10.6)') real(sum)
                  if ( len_trim(str1) .ge. 2) then
                    str1 = trim(str1)//"&"//trim(str2)
                  else
                    str1 = trim(str2)
                  end if
                  if (iatom2 .eq. wan%natom_wan .and. il2 .eq. wan%nbl_atom_wan(iatom2) .and. im2 &
   &                       .eq. 2*wan%latom_wan(iatom2)%lcalc(il2)+1 .and. pos2 .eq. size(wan%nposition(iatom2)%pos,1)) then
                    str1 = trim(str1)//'\\'
                  end if
                end do
              end do
            end do
          end do
          write(unt,'(a)') trim(str1)
        end do
        if (iatom1 .ne. wan%natom_wan .or. il1 .ne. wan%nbl_atom_wan(iatom1) .or. pos1 .ne. size(wan%nposition(iatom1)%pos,1)) then
          write(unt,'(a)') '\hline'
        end if
      end do
    end do
  end do



  write(unt,'(a)') '\end{array} \right) $ }'

  if (name(len_trim(name)-2:len(trim(name))) == 'gen') then
    write(unt,'(a,i0,a,F7.3,a)')     "\caption{Energy matrix in real space for isppol = ", &
&    isppol," in a ",wan%acell(1), " a.u. cell}"
  end if

  if (name(len_trim(name)-2:len(trim(name))) == 'occ') then
    write(unt,'(a,i0,a,F7.3,a)')     "\caption{Occupation matrix in real space for isppol = ",&
&    isppol," in a ",wan%acell(1), " a.u. cell}"
  end if



  write(unt,'(a)') '\end{figure}'
  write(unt,'(a)') '\end{document}'

end do
  close(unt)

end subroutine print_operwan
!!***


!!****f* m_plowannier/init_operwan_realspace
!! NAME
!!  init_operwan_realspace
!!
!! FUNCTION
!!  Initialize an operwan_realspace type variable
!!
!! INPUTS
!!  wan, operwan_realspace
!!
!! OUTPUT
!! operwan_realspace
!!
!! PARENTS
!! m_plowannier
!!
!! CHILDREN
!!
!! SOURCE
subroutine init_operwan_realspace(wan,oprs)

!Arguments----------------------------------
  type(operwan_realspace_type),intent(inout) :: oprs
  type(plowannier_type), intent(in) :: wan

!Local variables----------------------------
  integer :: i1,i2,n1,n2,p1,p2,l1,l2,sp,pi
  
 !variable names is shorten to achieve not too long line lenght
  sp=wan%nsppol
  pi=wan%nspinor
  ABI_DATATYPE_ALLOCATE(oprs%atom_index,(wan%natom_wan,wan%natom_wan))
  do i1 = 1,wan%natom_wan
    do i2 = 1,wan%natom_wan
      n1=size(wan%nposition(i1)%pos,1)
      n2=size(wan%nposition(i2)%pos,1)
      ABI_DATATYPE_ALLOCATE(oprs%atom_index(i1,i2)%position,(n1,n2))
      do p1 = 1,size(wan%nposition(i1)%pos,1)
        do p2 = 1,size(wan%nposition(i2)%pos,1)
          n1=wan%nbl_atom_wan(i1)
          n2=wan%nbl_atom_wan(i2)
          ABI_DATATYPE_ALLOCATE(oprs%atom_index(i1,i2)%position(p1,p2)%atom,(n1,n2))
          do l1 = 1,wan%nbl_atom_wan(i1)
            do l2 = 1,wan%nbl_atom_wan(i2)
              n1=2*wan%latom_wan(i1)%lcalc(l1)+1
              n2=2*wan%latom_wan(i2)%lcalc(l2)+1
              ABI_ALLOCATE(oprs%atom_index(i1,i2)%position(p1,p2)%atom(l1,l2)%matl,(n1,n2,sp,pi,pi))
              oprs%atom_index(i1,i2)%position(p1,p2)%atom(l1,l2)%matl = czero
            end do
         end do
       end do
     end do
   end do
 end do

end subroutine init_operwan_realspace
!!***

!!****f* m_plowannier/reduce_operwan_realspace
!! NAME
!!  reduce_operwan_realspace
!!
!! FUNCTION
!!  reduce a table of operwan_realspace type
!!
!! INPUTS
!!  wan,rhot1,npwx,nibz,comm,nbz,nsppol
!!
!! OUTPUT
!! rhot1
!!
!! PARENTS
!! m_sigma_driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine reduce_operwan_realspace(wan,rhot1,npwx,nibz,comm,nbz,nsppol)


  use m_xmpi, only : xmpi_barrier,xmpi_sum  
!Arguments---------------------------------------
  type(plowannier_type),intent(in) :: wan
  integer, intent(in) :: npwx,nibz,comm,nbz,nsppol
  type(operwan_realspace_type),target,intent(inout) :: rhot1(npwx,nibz)
!Local variables----------------------------------
  complex(dpc),allocatable ::  buffer(:)
  integer :: dim,pwx,ibz, spin, ispinor1, ispinor2, iatom1, iatom2, pos1, pos2
  integer :: il1, il2, im1, im2, nnn, ierr
  complex(dpc),pointer :: oper_ptr(:,:,:,:,:)

  
   dim=0
     do pwx=1,npwx
     do ibz=1,nibz
       do spin=1,wan%nsppol
       do ispinor1=1,wan%nspinor
       do ispinor2=1,wan%nspinor
         do iatom1=1,wan%natom_wan
         do iatom2=1,wan%natom_wan
           do pos1=1,size(wan%nposition(iatom1)%pos,1)
           do pos2=1,size(wan%nposition(iatom2)%pos,1)
             do il1=1,wan%nbl_atom_wan(iatom1)
             do il2=1,wan%nbl_atom_wan(iatom2)
               do im1=1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
               do im2=1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
     dim=dim+1
               enddo!im2
               enddo!im1
             enddo!il2
             enddo!il1
           enddo!pos2
           enddo!pos1
         enddo!iatom2
         enddo!iatom1
       enddo!ispinor2
       enddo!ispinor1
       enddo!spin
     enddo!ibz
     enddo!pwx
     ABI_ALLOCATE(buffer,(dim))
     nnn=0
     do pwx=1,npwx
     do ibz=1,nibz         
       do iatom1=1,wan%natom_wan
       do iatom2=1,wan%natom_wan
         do pos1=1,size(wan%nposition(iatom1)%pos,1)
         do pos2=1,size(wan%nposition(iatom2)%pos,1)
           do il1=1,wan%nbl_atom_wan(iatom1)
           do il2=1,wan%nbl_atom_wan(iatom2)
             oper_ptr=>rhot1(pwx,ibz)%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl
             do im1=1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
             do im2=1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
               do spin=1,wan%nsppol
                 do ispinor1=1,wan%nspinor
                 do ispinor2=1,wan%nspinor
     nnn=nnn+1
     buffer(nnn)=oper_ptr(im1,im2,spin,ispinor1,ispinor2)
                 enddo!ispinor2
                 enddo!ispinor1
               enddo!spin
             enddo!im2
             enddo!im1
           enddo!il2
           enddo!il1
         enddo!pos2
         enddo!pos1
       enddo!iatom2
       enddo!iatom1
     enddo!ibz  
     enddo!pwx
     call xmpi_barrier(comm)
     call xmpi_sum(buffer,comm,ierr)
     call xmpi_barrier(comm)
     buffer=buffer/nbz/nsppol
     nnn=0
     do pwx=1,npwx
     do ibz=1,nibz
         do iatom1=1,wan%natom_wan
         do iatom2=1,wan%natom_wan
           do pos1=1,size(wan%nposition(iatom1)%pos,1)
           do pos2=1,size(wan%nposition(iatom2)%pos,1)
             do il1=1,wan%nbl_atom_wan(iatom1)
             do il2=1,wan%nbl_atom_wan(iatom2)
               oper_ptr=>rhot1(pwx,ibz)%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl
               do im1=1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
               do im2=1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                 do spin=1,wan%nsppol
                   do ispinor1=1,wan%nspinor
                   do ispinor2=1,wan%nspinor
      nnn=nnn+1
      oper_ptr(im1,im2,spin,ispinor1,ispinor2)=buffer(nnn)
               enddo!im2
               enddo!im1
             enddo!il2
             enddo!il1
           enddo!pos2
           enddo!pos1
         enddo!iatom2
         enddo!iatom1
       enddo!ispinor2
       enddo!ispinor1
       enddo!spin
     enddo!ibz
     enddo!pwx
     ABI_DEALLOCATE(buffer)


end subroutine reduce_operwan_realspace
!!***

!!****f* m_plowannier/destroy_operwan_realspace
!! NAME
!!  destroy_operwan_realspace
!!
!! FUNCTION
!!  Destroy an operwan_realspace type variable
!!
!! INPUTS
!!  wan, operwan_realspace
!!
!! OUTPUT
!! operwan_realspace
!!
!! PARENTS
!! m_plowannier
!!
!! CHILDREN
!!
!! SOURCE
subroutine destroy_operwan_realspace(wan,operwan_realspace)

!Arguments----------------------------------
  type(operwan_realspace_type),intent(inout) :: operwan_realspace
  type(plowannier_type), intent(in) :: wan

!Local variables----------------------------
  integer :: iatom1,iatom2,pos1,pos2,il1,il2


 do iatom1 = 1,wan%natom_wan
   do iatom2 = 1,wan%natom_wan
     do pos1 = 1,size(wan%nposition(iatom1)%pos,1)
       do pos2 = 1,size(wan%nposition(iatom2)%pos,1)
         do il1 = 1,wan%nbl_atom_wan(iatom1)
           do il2 = 1,wan%nbl_atom_wan(iatom2)
            ABI_DEALLOCATE(operwan_realspace%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl)
           end do
         end do
         ABI_DATATYPE_DEALLOCATE(operwan_realspace%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom)
       end do
     end do
     ABI_DATATYPE_DEALLOCATE(operwan_realspace%atom_index(iatom1,iatom2)%position)
   end do
 end do
 ABI_DATATYPE_DEALLOCATE(operwan_realspace%atom_index)

end subroutine destroy_operwan_realspace
!!***


!!****f* m_plowannier/zero_operwan_realspace
!! NAME
!!  zero_operwan_realspace
!!
!! FUNCTION
!!  Set an operwan_realspace to zero
!!
!! INPUTS
!!  wan, operwan_realspace
!!
!! OUTPUT
!! operwan_realspace
!!
!! PARENTS
!! m_plowannier
!!
!! CHILDREN
!!
!! SOURCE
subroutine zero_operwan_realspace(wan,operwan_realspace)

!Arguments----------------------------------
  type(operwan_realspace_type),intent(inout) :: operwan_realspace
  type(plowannier_type), intent(in) :: wan

!Local variables----------------------------
  integer :: isppol,iatom1,iatom2,pos1,pos2,il1,il2,im1,im2,ispinor1,ispinor2


  do isppol = 1,wan%nsppol
    do ispinor1=1,wan%nspinor
      do ispinor2=1,wan%nspinor
        do iatom1 = 1,wan%natom_wan
          do pos1 = 1,size(wan%nposition(iatom1)%pos,1)
            do il1 = 1,wan%nbl_atom_wan(iatom1)
              do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
                do iatom2 = 1,wan%natom_wan
                  do pos2 = 1,size(wan%nposition(iatom2)%pos,1)
                    do il2 = 1,wan%nbl_atom_wan(iatom2)
                      do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                        operwan_realspace%atom_index(iatom1,iatom2)%position(pos1,pos2)%&
                          &atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2)=czero
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
end subroutine zero_operwan_realspace
!!***




!!****f* m_plowannier/compute_oper_wank2realspace
!! NAME
!!  compute_operwan_wanK2realspace
!!
!! FUNCTION
!! Compute an operator from WannierK space to real space
!!
!! INPUTS
!!  wan,operwan, operwan_realspace
!!
!! OUTPUT
!! operwan_realspace
!!
!! PARENTS
!! m_plowannier
!!
!! CHILDREN
!!
!! SOURCE
subroutine compute_oper_wank2realspace(wan,operwan,operwan_realspace)

!Arguments----------------------------------
  type(operwan_realspace_type),intent(inout) :: operwan_realspace
  type(operwan_type),intent(in) :: operwan(:,:,:) 
  type(plowannier_type), intent(in) :: wan
  

!Local variables----------------------------
  integer :: isppol,iatom1,pos1,il1,im1,iatom2,pos2,il2,im2,ikpt,ispinor1,ispinor2
  

  do isppol = 1,wan%nsppol
    do ispinor1=1,wan%nspinor
      do ispinor2=1,wan%nspinor
        do iatom1 = 1,wan%natom_wan
          do pos1 = 1,size(wan%nposition(iatom1)%pos,1)
            do il1 = 1,wan%nbl_atom_wan(iatom1)
              do im1 = 1,2*wan%latom_wan(iatom1)%lcalc(il1)+1
                do iatom2 = 1,wan%natom_wan
                  do pos2 = 1,size(wan%nposition(iatom2)%pos,1)
                    do il2 = 1,wan%nbl_atom_wan(iatom2)
                      do im2 = 1,2*wan%latom_wan(iatom2)%lcalc(il2)+1
                       !sum over ikpt
                        do ikpt = 1,wan%nkpt
                          operwan_realspace%atom_index(iatom1,iatom2)%position(pos1,pos2)%&
                            &atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2) =&
                            operwan_realspace%atom_index(iatom1,iatom2)%position(pos1,pos2)%&
                            &atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2)&
                            + real(operwan(ikpt,iatom1,iatom2)%atom(il1,il2)%matl(im1,im2,isppol,ispinor1,ispinor2)&
                            * wan%wtk(ikpt) * exp( cmplx(0.0,1.0) * two_pi * ( &
                            wan%kpt(1,ikpt) * ( wan%nposition(iatom1)%pos(pos1,1) - wan%nposition(iatom2)%pos(pos2,1) )+&
                            wan%kpt(2,ikpt) * ( wan%nposition(iatom1)%pos(pos1,2) - wan%nposition(iatom2)%pos(pos2,2) )+&
                            wan%kpt(3,ikpt) * ( wan%nposition(iatom1)%pos(pos1,3) - wan%nposition(iatom2)%pos(pos2,3)))))
                        end do
                       !end of the sum
                      end do
                    enddo
                  enddo
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end subroutine compute_oper_wank2realspace
!!***
END MODULE m_plowannier
!!***
