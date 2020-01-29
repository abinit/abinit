!!****m* ABINIT/m_ifc
!! NAME
!!  m_ifc
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  used to handle interatomic force constant sets
!!
!! COPYRIGHT
!! Copyright (C) 2011-2020 ABINIT group (XG,MJV,EB,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_ifc

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_sort
 use m_cgtools
 use m_lebedev
 use m_nctk
 use m_ddb
 use m_ddb_hdr
 use m_symkpt
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_io_tools,    only : open_file
 use m_fstrings,    only : ktoa, int2char4, sjoin, itoa, ltoa, ftoa
 use m_symtk,       only : matr3inv
 use m_special_funcs,  only : abi_derfc
 use m_time,        only : cwtime, cwtime_report, timab
 use m_copy,        only : alloc_copy
 use m_pptools,     only : printbxsf
 use m_ewald,       only : ewald9
 use m_crystal,     only : crystal_t
 use m_geometry,    only : phdispl_cart2red, normv, mkrdim
 use m_kpts,        only : kpts_ibz_from_kptrlatt, smpbz
 use m_dynmat,      only : canct9, dist9 , ifclo9, axial9, q0dy3_apply, q0dy3_calc, asrif9, dynmat_dq, &
                           make_bigbox, canat9, chkrp9, ftifc_q2r, wght9, symdm9, nanal9, gtdyn9, dymfz9, &
                           massmult_and_breaksym, dfpt_phfrq, dfpt_prtph

 implicit none

 private
!!***

!!****t* m_ifc/ifc_type
!! NAME
!! ifc_type
!!
!! FUNCTION
!!  Contains the necessary data to interpolate the
!!  phonon bandstructure and eigenvectors in reciprocal space (ie.
!!  interatomic force constants and corresponding real space grid info).
!!
!! SOURCE

 type,public :: ifc_type

   integer :: natom
     ! Number of atoms in the unit cell.

   integer :: mpert
     ! Maximum number of ipert.

   integer :: asr
     ! Option for the treatment of the Acoustic Sum Rule.

   integer :: brav
     ! Option for the sampling of the BZ (anaddb input variable)

   integer :: dipdip
     ! dipole dipole interaction flag.

   integer :: dipquad
     ! dipole quadrupole interaction flag.

   integer :: quadquad
     ! dipole quadrupole interaction flag.

   integer :: symdynmat
     ! If equal to 1, the dynamical matrix is symmetrized in dfpt_phfrq before the diagonalization.

   integer :: nqshft
     ! Number of shifts in the q-mesh (usually 1 since the mesh is gamma-centered!)

   integer :: nqibz
     ! Number of points in the IBZ

   integer :: nqbz
     ! Number of points in the full BZ

   integer :: nrpt
     ! Number of real space points used to integrate IFC (for interpolation of dynamical matrices)

   integer :: ngqpt(3)
    ! Number of division in the Q mesh.

   integer :: ewald_option
    ! Option for the ewald sum

   real(dp) :: rprim(3,3),gprim(3,3),acell(3)
     ! These values are used to call anaddb routines that don't use rprimd, gprimd

   real(dp) :: dielt(3,3)
     ! Dielectric tensor (Cartesian coordinates)

   real(dp) :: omega_minmax(2)
     ! Min and max frequency obtained on the initial ab-initio q-mesh (-+ 30 cmm1)
     ! Used to generate frequency meshes for DOSes.

   real(dp) :: r_inscribed_sphere
     ! radius of biggest sphere inscribed in the WS supercell

   real(dp),allocatable :: amu(:)
     ! amu(ntypat)
     ! mass of the atoms (atomic mass unit)

   real(dp),allocatable :: atmfrc(:,:,:,:,:)
     ! atmfrc(3,natom,3,natom,nrpt)
     ! Inter atomic forces in real space

   integer,allocatable :: cell(:,:)
     ! cell(nrpt,3)
     ! Give the index of the the cell and irpt

   real(dp),allocatable :: ewald_atmfrc(:,:,:,:,:)
     ! Ewald_atmfrc(3,natom,3,natom,nrpt)
     ! Ewald Inter atomic forces in real space

   real(dp),allocatable :: short_atmfrc(:,:,:,:,:)
     ! short_atmfrc(3,natom,3,natom,nrpt)
     ! Short range part of Inter atomic forces in real space

   real(dp),allocatable :: qshft(:,:)
    ! qshft(3,nqshft)
    ! The shifts of the q-mesh

   real(dp), allocatable :: rpt(:,:)
     ! rpt(3,nrpt)
     ! Real space points in canonical type coordinates.

   real(dp),allocatable :: wghatm(:,:,:)
     ! wghatm(natom,natom,nrpt)
     ! Weights for each point and atom in the Wigner Seitz supercell in real space.

   real(dp),allocatable :: rcan(:,:)
     ! rcan(3,natom)
     ! Atomic position in canonical coordinates.

   real(dp),allocatable :: trans(:,:)
     ! trans(3,natom)
     ! Atomic translations: xred = rcan + trans

   real(dp),allocatable :: dyewq0(:,:,:)
     ! dyewq0(3,3,natom)
     ! Atomic self-interaction correction to the dynamical matrix (only when dipdip=1).

   real(dp),allocatable :: zeff(:,:,:)
     ! zeff(3,3,natom)
     ! Effective charge on each atom, versus electric field and atomic displacement.
     ! Cartesian coordinates

   real(dp),allocatable :: qdrp_cart(:,:,:,:)
     ! qdrp_cart(3,3,3,natom)
     ! Quadrupole tensor on each atom
     ! Cartesian coordinates

   real(dp),allocatable :: qibz(:,:)
     ! qibz(3,nqibz))
     ! List of q-points in the IBZ

   real(dp),allocatable :: wtq(:)
     ! wtq(nqibz))
     ! q-point Weights.

   real(dp),allocatable :: qbz(:,:)
     ! qbz(3,nqbz))
     ! List of q-points in the full BZ

   real(dp),allocatable :: dynmat(:,:,:,:,:,:)
     ! dynmat(2,3,natom,3,natom,nqbz))
     ! dynamical matrices relative to the q points of the B.Z. sampling
     ! Note that the long-range dip-dip part has been removed if dipdip=1
     ! Moreover the array is multiplied by a phase shift in mkifc9.

   !real(dp),allocatable :: dynmat_lr(:,:,:,:,:,:)
    ! dynmat_lr(2,3,natom,3,natom,nqbz))
    ! Long-range part of dynmat in q-space

 contains

    procedure :: free => ifc_free
    ! Release memory

    procedure :: print => ifc_print
     ! Print info on the object.

    procedure :: fourq => ifc_fourq
     ! Use Fourier interpolation to compute interpolated frequencies w(q) and eigenvectors e(q)

    procedure :: speedofsound => ifc_speedofsound
     ! Compute the speed of sound by averaging phonon group velocities.

    procedure :: write => ifc_write
     ! Print the ifc (output, netcdf and text file)

    procedure :: outphbtrap => ifc_outphbtrap
     ! Print out phonon frequencies on regular grid for BoltzTrap code.

    procedure :: printbxsf => ifc_printbxsf
     ! Output phonon isosurface in Xcrysden format.

    procedure :: calcnwrite_nana_terms => ifc_calcnwrite_nana_terms
     ! Compute phonons for q--> 0 with LO-TO

 end type ifc_type

 public :: ifc_init          ! Constructor from DDB datatype
 public :: ifc_init_fromFile ! Constructor from filename
!!***

!----------------------------------------------------------------------

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_free
!! NAME
!! ifc_free
!!
!! FUNCTION
!!  Deallocate memory for the ifc_type structure
!!
!! PARENTS
!!      anaddb,compute_anharmonics,eph,m_anharmonics_terms
!!      m_effective_potential,m_effective_potential_file,m_gruneisen
!!      m_harmonics_terms,m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_free(ifc)

!Arguments ------------------------------------
 class(ifc_type),intent(inout) :: ifc

! ************************************************************************

 ABI_SFREE(ifc%amu)
 ABI_SFREE(ifc%atmfrc)
 ABI_SFREE(ifc%cell)
 ABI_SFREE(ifc%ewald_atmfrc)
 ABI_SFREE(ifc%short_atmfrc)
 ABI_SFREE(ifc%qshft)
 ABI_SFREE(ifc%rpt)
 ABI_SFREE(ifc%wghatm)
 ABI_SFREE(ifc%rcan)
 ABI_SFREE(ifc%trans)
 ABI_SFREE(ifc%dyewq0)
 ABI_SFREE(ifc%qibz)
 ABI_SFREE(ifc%wtq)
 ABI_SFREE(ifc%qbz)
 ABI_SFREE(ifc%zeff)
 ABI_SFREE(ifc%qdrp_cart)
 ABI_SFREE(ifc%dynmat)
 !ABI_SFREE(ifc%dynmat_lr)

end subroutine ifc_free
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_init
!! NAME
!!  ifc_init
!!
!! FUNCTION
!!  Initialize the dynamical matrix as well as the IFCs.
!!  taking into account the dipole-dipole, dipole-quadrupole and quadrupole-quadrupole 
!!  interaction.
!!
!! INPUTS
!! crystal<type(crystal_t)> = Information on the crystalline structure.
!! ddb<type(ddb_type)> = Database with derivatives.
!! brav=bravais lattice (1 or -1=simple lattice,2=face centered lattice, 3=centered lattice,4=hexagonal lattice)
!! asr= Option for the imposition of the ASR
!!   0 => no ASR,
!!   1 => modify "asymmetrically" the diagonal element
!!   2 => modify "symmetrically" the diagonal element
!! symdynmat=if 1, (re)symmetrize the dynamical matrix, except if Gamma wavevector with electric field added.
!! dipdip=
!!   if 0, no dipole-dipole interaction was subtracted in atmfrc
!!   if 1, atmfrc has been build without dipole-dipole part
!! rfmeth =
!!   1 if non-stationary block
!!   2 if stationary block
!!   3 if third order derivatives
!! dielt(3,3)=dielectric tensor.
!! zeff(3,3,natom)=effective charge on each atom, versus electric field and atomic displacement
!! prtsrlr: TODO: TO BE REMOVED
!! enunit: TODO: TO BE REMOVED
!! dielt(3,3)=dielectric tensor
!! ngqpt_in = input values of ngqpt
!! nqshft=Number of shifths in q-grid.
!! q1shft(3,nqshft)=Shifts for q-grid
!! nsphere=number of atoms to be included in the cut-off sphere for interatomic force constant.
!!    0: maximum extent allowed by the grid.
!!  > 0: Apply cutoff
!!   -1: Analyze the effect of different nsphere values on the phonon spectrum, in particular the
!!       frequencies around gamma.
!! rifcsph=radius for cutoff of IFC.
!! comm=MPI communicator.
!! [Ifc_coarse]=Optional.
!! [dipquad] = if 1, atmfrc has been build without dipole-quadrupole part
!! [quadquad] = if 1, atmfrc has been build without quadrupole-quadrupole part
!!
!! OUTPUT
!! Ifc<ifc_type>=Object containing the dynamical matrix and the IFCs.
!!
!! PARENTS
!!      anaddb,eph,m_effective_potential_file,m_gruneisen,m_tdep_abitypes
!!
!! CHILDREN
!!
!! SOURCE

#ifdef MR_DEV
subroutine ifc_init(ifc,crystal,ddb,brav,asr,symdynmat,dipdip,&
  rfmeth,ngqpt_in,nqshft,q1shft,dielt,zeff,qdrp_cart,nsphere,rifcsph,&
  prtsrlr,enunit, & ! TODO: TO BE REMOVED
  comm, &
  Ifc_coarse,dipquad,quadquad) ! Optional
#else
subroutine ifc_init(ifc,crystal,ddb,brav,asr,symdynmat,dipdip,&
  rfmeth,ngqpt_in,nqshft,q1shft,dielt,zeff,qdrp_cart,nsphere,rifcsph,&
  prtsrlr,enunit, & ! TODO: TO BE REMOVED
  comm, &
  Ifc_coarse) ! Optional
#endif

!Arguments ------------------------------------
 integer,intent(in) :: asr,brav,dipdip,symdynmat,nqshft,rfmeth,nsphere,comm
 real(dp),intent(in) :: rifcsph
 type(ifc_type),intent(inout) :: Ifc
 type(crystal_t),intent(in) :: Crystal
 type(ddb_type),intent(in) :: ddb
 type(ifc_type),optional,intent(in) :: Ifc_coarse
#ifdef MR_DEV
 integer,optional,intent(in) :: dipquad, quadquad
#endif

!arrays
 integer,intent(in) :: ngqpt_in(3)
 real(dp),intent(in) :: q1shft(3,nqshft)
 real(dp),intent(in) :: dielt(3,3),zeff(3,3,Crystal%natom)
 real(dp),intent(in) :: qdrp_cart(3,3,3,Crystal%natom)
!anaddb variables (TO BE REMOVED)
 integer,intent(in) :: prtsrlr,enunit
!end anaddb variables

!Local variables -------------------------
!scalars
 integer,parameter :: timrev1=1,iout0=0,chksymbreak0=0
 integer :: mpert,iout,iqpt,mqpt,nsym,ntypat,iq_ibz,iq_bz,ii,natom
 integer :: nqbz,option,plus,sumg0,irpt,irpt_new
 integer :: nprocs,my_rank,my_ierr,ierr
 real(dp),parameter :: qphnrm=one
 real(dp) :: xval,cpu,wall,gflops,rcut_min
 real(dp) :: r_inscribed_sphere,toldist
 character(len=500) :: msg
 type(ifc_type) :: ifc_tmp
!arrays
 integer :: ngqpt(9),qptrlatt(3,3)
 integer,allocatable :: qmissing(:),ibz2bz(:),bz2ibz_smap(:,:)
 real(dp) :: gprim(3,3),rprim(3,3),qpt(3),rprimd(3,3)
 real(dp):: rcan(3,Crystal%natom),trans(3,Crystal%natom),dyewq0(3,3,Crystal%natom)
 real(dp) :: displ_cart(2*3*Crystal%natom*3*Crystal%natom)
 real(dp) :: phfrq(3*Crystal%natom)
 real(dp) :: eigvec(2,3,Crystal%natom,3,Crystal%natom)
 real(dp),allocatable :: dyew(:,:,:,:,:),out_d2cart(:,:,:,:,:)
 real(dp),allocatable :: dynmatfull(:,:,:,:,:,:),dynmat_sr(:,:,:,:,:,:),dynmat_lr(:,:,:,:,:,:) ! for OmegaSRLR
 real(dp),allocatable :: wtq(:),wtq_folded(:),qbz(:,:)

!******************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 call cwtime(cpu, wall, gflops, "start")

 ! TODO: This dimension should be encapsulated somewhere. We don't want to
 ! change the entire code if someone adds a new kind of perturbation.
#ifdef MR_DEV
 mpert = Crystal%natom + MPERT_MAX; iout = ab_out
#else
 mpert = Crystal%natom + 6; iout = ab_out
#endif

 rprim = ddb%rprim; gprim = ddb%gprim

 nsym = Crystal%nsym
 natom = Crystal%natom
 ntypat = Crystal%ntypat
 rprimd = Crystal%rprimd

 ngqpt=0; ngqpt(1:3)=ngqpt_in(1:3)

! Copy important parameters in Ifc
 Ifc%natom = natom
 Ifc%mpert = mpert
 Ifc%asr = asr
 Ifc%brav = brav
 Ifc%dipdip = abs(dipdip)
#ifdef MR_DEV
 if (present(dipquad)) Ifc%dipquad = dipquad
 if (present(quadquad)) Ifc%quadquad = quadquad
#endif
 Ifc%symdynmat = symdynmat
 Ifc%nqshft = nqshft
 call alloc_copy(q1shft(:,1:Ifc%nqshft),Ifc%qshft)
 Ifc%ngqpt = ngqpt_in(1:3)
 Ifc%rprim = ddb%rprim
 Ifc%gprim = ddb%gprim
 Ifc%acell = ddb%acell
 Ifc%ewald_option = 0; if (dipdip<0) Ifc%ewald_option = 1 !HM TODO: expose this in the init?

 ! Check if the rprim are coherent with the choice used in the interatomic forces generation
 call chkrp9(Ifc%brav,rprim)

 dyewq0 = zero
 if (Ifc%dipdip==1 .and. (Ifc%asr==1.or.Ifc%asr==2)) then
   ! Calculation of the non-analytical part for q=0
   sumg0=0
   qpt(:)=zero
   ABI_MALLOC(dyew,(2,3,natom,3,natom))
#ifdef MR_DEV
   if (present(dipquad).or.present(quadquad)) then
   call ewald9(ddb%acell,dielt,dyew,Crystal%gmet,gprim,natom,qpt,Crystal%rmet,rprim,sumg0,Crystal%ucvol,&
               Crystal%xred,zeff,qdrp_cart,ifc%ewald_option,dipquad=dipquad,quadquad=quadquad)
   else
   call ewald9(ddb%acell,dielt,dyew,Crystal%gmet,gprim,natom,qpt,Crystal%rmet,rprim,sumg0,Crystal%ucvol,&
               Crystal%xred,zeff,qdrp_cart,ifc%ewald_option)
   end if
#else
   call ewald9(ddb%acell,dielt,dyew,Crystal%gmet,gprim,natom,qpt,Crystal%rmet,rprim,sumg0,Crystal%ucvol,&
               Crystal%xred,zeff,qdrp_cart,ifc%ewald_option)
#endif

   call q0dy3_calc(natom,dyewq0,dyew,Ifc%asr)
   ABI_FREE(dyew)
 end if

 ! Sample the Brillouin zone
 option=1
 qptrlatt = 0; qptrlatt(1,1)=ngqpt(1); qptrlatt(2,2)=ngqpt(2); qptrlatt(3,3)=ngqpt(3)
 mqpt=ngqpt(1)*ngqpt(2)*ngqpt(3)*nqshft
 if (Ifc%brav==2) mqpt=mqpt/2
 if (Ifc%brav==3) mqpt=mqpt/4

 ABI_MALLOC(qbz,(3,mqpt))
 call smpbz(Ifc%brav,ab_out,qptrlatt,mqpt,nqbz,nqshft,option,q1shft,qbz)

 ! Find the irreducible zone (qibz)
 ABI_MALLOC(ibz2bz, (nqbz))
 ABI_MALLOC(wtq_folded, (nqbz))
 ABI_MALLOC(wtq, (nqbz))
 wtq = one / nqbz ! Weights sum up to one
 ABI_MALLOC(bz2ibz_smap, (6, nqbz))

 call symkpt(chksymbreak0,crystal%gmet,ibz2bz,iout0,qbz,nqbz,ifc%nqibz,crystal%nsym,&
   crystal%symrec,timrev1,wtq,wtq_folded, bz2ibz_smap, xmpi_comm_self)

 ABI_FREE(bz2ibz_smap)

 ABI_MALLOC(ifc%qibz, (3,ifc%nqibz))
 ABI_MALLOC(ifc%wtq, (ifc%nqibz))
 do iq_ibz=1,ifc%nqibz
   ifc%qibz(:,iq_ibz) = qbz(:, ibz2bz(iq_ibz))
   ifc%wtq(iq_ibz) = wtq_folded(ibz2bz(iq_ibz))
 end do
 ABI_FREE(ibz2bz)
 ABI_FREE(wtq_folded)
 ABI_FREE(wtq)

 ABI_MALLOC(Ifc%dynmat,(2,3,natom,3,natom,nqbz))

! Find symmetrical dynamical matrices
 if (.not.present(Ifc_coarse)) then

   ! Each q-point in the BZ mush be the symmetrical of one of the qpts in the ddb file.
   call symdm9(ddb%flg,ddb%nrm,ddb%qpt,ddb%typ,ddb%val,&
&    Ifc%dynmat,gprim,Crystal%indsym,mpert,natom,ddb%nblok,nqbz,nsym,rfmeth,rprim,qbz,&
&    Crystal%symrec, Crystal%symrel, comm)

 else
   ! Symmetrize the qpts in the BZ using the q-points in the ddb.
   ! Then use Ifc_coarse to fill the missing entries with Fourier interpolated matrices.
   !
   ! TODO: The previous version of refineblk was hacking the DDB database to add the q-points in the **IBZ**
   ! Then D(q) was symmetrized in symdm9. This version avoids the symmetrization: the q-points
   ! in the BZ that are not in the coarse q-mesh are obtained by an explicit FT.
   ! This means that the final D(q) may break some symmetry in q-space if the FT does not preserve it.
   ! The most elegant approach would be to get D(q_ibz) via FT if q_ibz is not in the coarse mesh and then
   ! call symdm9 to get D(q) for each q point in the star of q_ibz.
   call wrtout(std_out,"Will fill missing qpoints in the full BZ using the coarse q-mesh","COLL")

   call symdm9(ddb%flg,ddb%nrm,ddb%qpt,ddb%typ,ddb%val,&
     Ifc%dynmat,gprim,Crystal%indsym,mpert,natom,ddb%nblok,nqbz,nsym,rfmeth,rprim,qbz,&
     Crystal%symrec,Crystal%symrel,comm, qmissing=qmissing)

   ! Compute dynamical matrix with Fourier interpolation on the coarse q-mesh.
   write(msg,"(a,i0,a)")"Will use Fourier interpolation to construct D(q) for ",size(qmissing)," q-points"
   call wrtout(std_out,msg)

   ABI_MALLOC(out_d2cart, (2,3,natom,3,natom))
   do ii=1,size(qmissing)
     iq_bz = qmissing(ii)
     qpt = qbz(:,iq_bz)
     ! TODO: check dipdip option and phase, but I think this is correct!
     call ifc_fourq(Ifc_coarse,Crystal,qpt,phfrq,displ_cart,out_d2cart=out_d2cart)
     Ifc%dynmat(:,:,:,:,:,iq_bz) = out_d2cart
   end do

   ABI_FREE(qmissing)
   ABI_FREE(out_d2cart)
 end if

 ! OmegaSRLR: Store full dynamical matrix for decomposition into short- and long-range parts
 ABI_MALLOC(dynmatfull,(2,3,natom,3,natom,nqbz))
 dynmatfull=Ifc%dynmat

 if (Ifc%dipdip==1) then
   ! Take off the dipole-dipole part of the dynamical matrix
   call wrtout(std_out, " Will extract the dipole-dipole part for every wavevector")
   ABI_MALLOC(dyew,(2,3,natom,3,natom))

   do iqpt=1,nqbz
     if (mod(iqpt, nprocs) /= my_rank) then ! mpi-parallelism
       ifc%dynmat(:,:,:,:,:,iqpt) = zero; cycle
     end if
     qpt(:)=qbz(:,iqpt)
     sumg0=0
#ifdef MR_DEV
     if (present(dipquad).or.present(quadquad)) then
       call ewald9(ddb%acell,dielt,dyew,Crystal%gmet,gprim,natom,qpt,Crystal%rmet,rprim,sumg0,Crystal%ucvol,&
                 Crystal%xred,zeff,qdrp_cart,ifc%ewald_option,dipquad=dipquad,quadquad=quadquad)
     else
       call ewald9(ddb%acell,dielt,dyew,Crystal%gmet,gprim,natom,qpt,Crystal%rmet,rprim,sumg0,Crystal%ucvol,&
                 Crystal%xred,zeff,qdrp_cart,ifc%ewald_option)
     end if
#else
     call ewald9(ddb%acell,dielt,dyew,Crystal%gmet,gprim,natom,qpt,Crystal%rmet,rprim,sumg0,Crystal%ucvol,&
                 Crystal%xred,zeff,qdrp_cart,ifc%ewald_option)
#endif
     call q0dy3_apply(natom,dyewq0,dyew)
     plus=0
     call nanal9(dyew,Ifc%dynmat,iqpt,natom,nqbz,plus)
   end do

   call xmpi_sum(ifc%dynmat, comm, ierr)
   ABI_FREE(dyew)
 end if

 ! OmegaSRLR: Store the short-range dynmat and compute long-range as difference
 ABI_MALLOC(dynmat_sr,(2,3,natom,3,natom,nqbz))
 ABI_MALLOC(dynmat_lr,(2,3,natom,3,natom,nqbz))
 dynmat_sr=Ifc%dynmat
 dynmat_lr=dynmatfull-dynmat_sr

 ! Now, take care of the remaining part of the dynamical matrix
 ! Move to canonical normalized coordinates
 call canat9(Ifc%brav,natom,rcan,rprim,trans,Crystal%xred)

 ! Multiply the dynamical matrix by a phase shift
 option=1
 call dymfz9(Ifc%dynmat,natom,nqbz,gprim,option,qbz,trans)

 ! Create the Big Box of R vectors in real space and compute the number of points (cells) in real space
 call make_bigbox(Ifc%brav,ifc_tmp%cell,ngqpt,nqshft,rprim,ifc_tmp%nrpt,ifc_tmp%rpt)

 ! Weights associated to these R points and to atomic pairs
 ABI_MALLOC(ifc_tmp%wghatm, (natom, natom, ifc_tmp%nrpt))

 ! HM: this tolerance is highly dependent on the compilation/architecture
 !     numeric errors in the DDB text file. Try a few tolerances and check whether all the weights are found.
 toldist = tol8
 do while (toldist <= tol6)
   ! Note ngqpt(9) with intent(inout)!
   call wght9(Ifc%brav,gprim,natom,ngqpt,nqbz,nqshft,ifc_tmp%nrpt,q1shft,rcan,&
              ifc_tmp%rpt,rprimd,toldist,r_inscribed_sphere,ifc_tmp%wghatm,my_ierr)
   call xmpi_max(my_ierr, ierr, comm, ii)
   if (ierr > 0) toldist = toldist * 10
   if (ierr == 0) exit
 end do

 if (ierr > 0) then
   write(msg, '(3a,es14.4,2a,i0, 14a)' ) &
    'The sum of the weight is not equal to nqpt.',ch10,&
    'The sum of the weights is: ',sum(ifc_tmp%wghatm),ch10,&
    'The number of q points is: ',nqbz, ch10, &
    'This might have several sources.',ch10,&
    'If toldist is larger than 1.0e-8, the atom positions might be loose.',ch10,&
    'and the q point weights not computed properly.',ch10,&
    'Action: make input atomic positions more symmetric.',ch10,&
    'Otherwise, you might increase "buffer" in m_dynmat.F90 see bigbx9 subroutine and recompile.',ch10,&
    'Actually, this can also happen when ngqpt is 0 0 0,',ch10,&
    'if abs(brav) /= 1, in which case you should change brav to 1.'
   MSG_ERROR(msg)
 end if

 ! Fourier transform of the dynamical matrices (q-->R)
 ABI_MALLOC(ifc_tmp%atmfrc, (3,natom,3,natom,ifc_tmp%nrpt))
 call ftifc_q2r(ifc_tmp%atmfrc,Ifc%dynmat,gprim,natom,nqbz,ifc_tmp%nrpt,ifc_tmp%rpt,qbz, comm)

 ! Eventually impose Acoustic Sum Rule on the interatomic forces
 if (Ifc%asr > 0) call asrif9(Ifc%asr,ifc_tmp%atmfrc,natom,ifc_tmp%nrpt,ifc_tmp%rpt,ifc_tmp%wghatm)

 ! The interatomic forces have been calculated
 write(msg, '(2a)')ch10,' The interatomic forces have been obtained '
 call wrtout([std_out, ab_out], msg,'COLL')
 call cwtime_report(" ifc_init1", cpu, wall, gflops)

 ! Apply cutoff on ifc if needed
 if (nsphere > 0 .or. abs(rifcsph) > tol10) then
   call wrtout(std_out, ' Apply cutoff on IFCs.')
   call wrtout(std_out, sjoin(" nsphere:", itoa(nsphere), ", rifcsph:", ftoa(rifcsph)))
   call corsifc9(ddb%acell,gprim,natom,ifc_tmp%nrpt,nsphere,rifcsph,rcan,rprim,ifc_tmp%rpt,rcut_min,ifc_tmp%wghatm)
   if (Ifc%asr > 0) then
     call wrtout(std_out, ' Enforcing ASR on cutoffed IFCs.')
     call asrif9(Ifc%asr,ifc_tmp%atmfrc,natom,ifc_tmp%nrpt,ifc_tmp%rpt,ifc_tmp%wghatm)
   end if
 end if

 ! Only conserve the necessary points in rpt: in the FT algorithm the order of the points is unimportant
 ! In the case of effective potential, we need to keep all the points
 Ifc%nrpt = 0
 do irpt=1,ifc_tmp%nrpt
   if (sum(ifc_tmp%wghatm(:,:,irpt)) /= 0) Ifc%nrpt = Ifc%nrpt+1
 end do

 ABI_CALLOC(Ifc%atmfrc,(3,natom,3,natom,Ifc%nrpt))
 ABI_CALLOC(Ifc%rpt,(3,Ifc%nrpt))
 ABI_CALLOC(Ifc%cell,(3,Ifc%nrpt))
 ABI_CALLOC(Ifc%wghatm,(natom,natom,Ifc%nrpt))
 ABI_CALLOC(Ifc%short_atmfrc,(3,natom,3,natom,Ifc%nrpt))
 ABI_CALLOC(Ifc%ewald_atmfrc,(3,natom,3,natom,Ifc%nrpt))

 irpt_new = 1
 do irpt = 1, ifc_tmp%nrpt
   if (sum(ifc_tmp%wghatm(:,:,irpt)) /= 0) then
     Ifc%atmfrc(:,:,:,:,irpt_new) = ifc_tmp%atmfrc(:,:,:,:,irpt)
     Ifc%rpt(:,irpt_new) = ifc_tmp%rpt(:,irpt)
     Ifc%wghatm(:,:,irpt_new) = ifc_tmp%wghatm(:,:,irpt)
     Ifc%cell(:,irpt_new) = ifc_tmp%cell(:,irpt)
     Ifc%r_inscribed_sphere = r_inscribed_sphere
     irpt_new = irpt_new + 1
   end if
 end do

 !write(std_out,*)"nrpt before filter:", ifc_tmp%nrpt, ", after: ", ifc%nrpt
 !do irpt=1,ifc%nrpt
 !  write(std_out,*)ifc%rpt(:,irpt), (ifc%wghatm(ii,ii,irpt), ii=1,natom)
 !end do

 ! Copy other useful arrays.
 Ifc%dielt = dielt
 Ifc%nqbz = nqbz

 call alloc_copy(rcan, Ifc%rcan)
 call alloc_copy(trans, Ifc%trans)
 call alloc_copy(dyewq0, Ifc%dyewq0)
 call alloc_copy(qbz(:,1:nqbz), Ifc%qbz)
 call alloc_copy(zeff, Ifc%zeff)
 call alloc_copy(qdrp_cart, Ifc%qdrp_cart)
 call alloc_copy(ddb%amu, Ifc%amu)

 call ifc_free(ifc_tmp)

 ! Compute min/max ph frequency with ab-initio q-mesh.
 ifc%omega_minmax(1) = huge(one); ifc%omega_minmax(2) = -huge(one)
 do iq_ibz=1,ifc%nqibz
   if (mod(iq_ibz, nprocs) /= my_rank) cycle ! mpi-parallelism
   call ifc_fourq(ifc, crystal, ifc%qibz(:,iq_ibz), phfrq, displ_cart)
   ifc%omega_minmax(1) = min(ifc%omega_minmax(1), minval(phfrq))
   ifc%omega_minmax(2) = max(ifc%omega_minmax(2), maxval(phfrq))
 end do
 xval = ifc%omega_minmax(1); call xmpi_min(xval, ifc%omega_minmax(1), comm, ierr)
 xval = ifc%omega_minmax(2); call xmpi_max(xval, ifc%omega_minmax(2), comm, ierr)
 ! Enlarge boundaries by 30 cm-1
 ifc%omega_minmax(1) = ifc%omega_minmax(1) - 30.0_dp/Ha_cmm1
 ifc%omega_minmax(2) = ifc%omega_minmax(2) + 30.0_dp/Ha_cmm1

 ! TODO (This is to be suppressed in a future version)
 if (prtsrlr == 1) then
   ! Check that the starting values are well reproduced.
   write(msg, '(2a)' )' mkifc9 : now check that the starting values ',&
     ' are reproduced after the use of interatomic forces '
   call wrtout(std_out, msg)
   do iqpt=1,nqbz
     qpt(:)=Ifc%qbz(:,iqpt)
     call ifc_fourq(Ifc,Crystal,qpt,phfrq,displ_cart,out_eigvec=eigvec)

     ! OmegaSRLR: Perform decomposition of dynamical matrix
     ! MG: FIXME I don't think the implementation is correct when q !=0
     if (prtsrlr==1) then
       call omega_decomp(ddb%amu,natom,ntypat,Crystal%typat,dynmatfull,dynmat_sr,dynmat_lr,iqpt,nqbz,eigvec)
     end if

     ! Write the phonon frequencies (this is for checking purposes).
     ! Note: these phonon frequencies are not written on unit iout, only on unit std_out.
     call dfpt_prtph(displ_cart,0,enunit,-1,natom,phfrq,qphnrm,qpt)
   end do
 end if

 ! OmegaSRLR: deallocate memory used by dynmat decomposition
 ABI_FREE(dynmatfull)
 ABI_FREE(dynmat_sr)
 ABI_FREE(dynmat_lr)
 ABI_FREE(qbz)

 if (nsphere == -1) call ifc_autocutoff(ifc, crystal, comm)

 call cwtime_report(" ifc_init2", cpu, wall, gflops)

end subroutine ifc_init
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_init_fromFile
!! NAME
!!  ifc_init_fromFile
!!
!! FUNCTION
!!  Need to be updated
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      anaddb,eph,m_effective_potential_file,m_gruneisen,m_tdep_abitypes
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_init_fromFile(dielt,filename,Ifc,natom,ngqpt,nqshift,qshift,ucell_ddb,zeff,qdrp_cart,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqshift,comm
 integer,intent(inout) :: natom
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: qshift(3,nqshift)
 character(len=*),intent(in) :: filename
 type(ifc_type),intent(out) :: Ifc
 real(dp),intent(inout) :: dielt(3,3)
 real(dp),allocatable,intent(inout) :: zeff(:,:,:)
 real(dp),allocatable,intent(inout) :: qdrp_cart(:,:,:,:)
 type(crystal_t),intent(out) :: ucell_ddb

!Local variables -------------------------
!scalars
 integer :: dipdip,i,iblok,iblok_tmp
 logical :: file_exists
 character(len=500) :: msg
 type(ddb_type) :: ddb
 type(ddb_hdr_type) :: ddb_hdr
!arrays
 integer,allocatable :: atifc(:)

!******************************************************************

 !check if ddb file exists
 inquire(file=filename, exist=file_exists)

 if (file_exists .eqv. .true.)then
   !Reading the ddb
   call ddb_hdr_open_read(ddb_hdr,filename,2,DDB_VERSION,comm,dimonly=1)

   natom = ddb_hdr%natom
   ABI_ALLOCATE(atifc,(ddb_hdr%natom))
   do i=1,ddb_hdr%natom
     atifc(i)=i
   end do

   call ddb_from_file(ddb,filename,1,ddb_hdr%natom,ddb_hdr%natom,atifc,ucell_ddb,comm)

 else
   MSG_ERROR(sjoin("File:", filename, "is not present in the directory"))
 end if

 ! Get Dielectric Tensor and Effective Charges
 ABI_ALLOCATE(zeff,(3,3,natom))
 ABI_ALLOCATE(qdrp_cart,(3,3,3,natom))
 iblok = ddb%get_dielt_zeff(ucell_ddb,1,1,0,dielt,zeff)
 iblok = ddb%get_quadrupoles(1,3,qdrp_cart)

 ! Try to get dielt, in case just the DDE are present
 if (iblok == 0) then
   iblok_tmp = ddb%get_dielt(1,dielt)
 end if

 ! ifc to be calculated for interpolation
 write(msg, '(a,a,(80a),a,a,a,a)' ) ch10,('=',i=1,80),ch10,ch10,' Calculation of the interatomic forces ',ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 if ((maxval(abs(zeff)) .lt. tol10) .OR. (maxval(dielt) .gt. 100000.0)) then
   dipdip=0
 else
   dipdip=1
 end if
 call ifc_init(Ifc,ucell_ddb,ddb,1,1,1,dipdip,1,ngqpt,nqshift,qshift,dielt,zeff,qdrp_cart,0,0.0_dp,0,1,comm)

 ! Free them all
 ABI_DEALLOCATE(atifc)
 call ddb%free()
 call ddb_hdr_free(ddb_hdr)

 end subroutine ifc_init_fromFile
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_print
!! NAME
!!  ifc_print
!!
!! FUNCTION
!!  Print info on the object
!!
!! INPUTS
!!  [unit]=Unit number for output. Defaults to std_out
!!  [prtvol]=Verbosity level.
!!  [header]=String to be printed as header for additional info.
!!
!! OUTPUT
!!  Only printing
!!
!! PARENTS
!!      anaddb,eph,m_tdep_abitypes
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_print(ifc, header, unit, prtvol)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=*),optional,intent(in) :: header
 class(ifc_type),intent(in) :: ifc

!Local variables-------------------------------
 integer :: unt,my_prtvol,iatom,ii,idir
 character(len=500) :: msg
! *********************************************************************

 unt = std_out; if (present(unit)) unt = unit
 my_prtvol = 0; if (present(prtvol)) my_prtvol = prtvol

 msg = ' ==== Info on the ifc% object ==== '
 if (present(header)) msg = ' ==== '//trim(adjustl(header))//' ==== '
 call wrtout(unt, msg)

 call wrtout(unt,' Real(R)+Recip(G) space primitive vectors, cartesian coordinates (Bohr,Bohr^-1):')
 do ii=1,3
   write(msg,'(1x,a,i1,a,3f11.7,2x,a,i1,a,3f11.7)')&
    'R(',ii,')=',ifc%rprim(:,ii),'G(',ii,')=',ifc%gprim(:,ii)
   call wrtout(unt,msg)
 end do
 call wrtout(unt, sjoin(" acell:", ltoa(ifc%acell)))

 call wrtout(unt, sjoin(" Acoustic Sum Rule option (asr):", itoa(ifc%asr)))
 call wrtout(unt, sjoin(" Option for the sampling of the BZ (brav):", itoa(ifc%brav)))
 call wrtout(unt, sjoin(" Symmetrization flag (symdynmat):", itoa(ifc%symdynmat)))
 call wrtout(unt, sjoin(" Dipole-dipole interaction flag (dipdip):", itoa(ifc%dipdip)))
 call wrtout(unt, sjoin(" Dielectric tensor: ", ch10, ltoa(reshape(ifc%dielt, [9]), fmt="f10.2")))
 call wrtout(unt, " Effective charges:")
 do iatom=1,ifc%natom
   call wrtout(unt, ltoa(reshape(ifc%zeff(:,:,iatom), [3*3]), fmt="f10.2"))
 end do
 call wrtout(unt, " Quadrupolar terms:")
 do iatom=1,ifc%natom
   do idir=1,3
     call wrtout(unt, ltoa(reshape(ifc%qdrp_cart(:,:,idir,iatom), [3*3]), fmt="f10.2"))
   end do
 end do

 call wrtout(unt, sjoin(" Mass of the atoms (atomic mass unit): ", ltoa(ifc%amu)))
 call wrtout(unt, sjoin(" Number of real-space points for IFC(R): ", itoa(ifc%nrpt)))
 call wrtout(unt, " ")

 call wrtout(unt, " Q-mesh:")
 call wrtout(unt, sjoin(" ngqpt:", ltoa(ifc%ngqpt),", nqshft:", itoa(ifc%nqshft)))
 do ii=1,ifc%nqshft
   call wrtout(unt, sjoin("  ", ktoa(ifc%qshft(:,ii))))
 end do

end subroutine ifc_print
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_fourq
!! NAME
!!  ifc_fourq
!!
!! FUNCTION
!!  Compute the phonon frequencies and the group velocities at the specified q-point by performing
!!  a Fourier transform on the IFCs matrix in real space.
!!
!! INPUTS
!!  Ifc<type(ifc_type)>=Object containing the dynamical matrix and the IFCs.
!!  Crystal<type(crystal_t)> = Information on the crystalline structure.
!!  qpt(3)=q-point in reduced coordinates (unless nanaqdir is specified)
!!  [nanaqdir]=If present, the qpt will be treated as a vector specifying the
!!    direction in q-space along which the non-analytic behaviour of the dynamical
!!    matrix will be treated. Possible values:
!!       "cart" if qpt defines a direction in Cartesian coordinates
!!       "reduced" if qpt defines a direction in reduced coordinates
!!  [comm]: MPI communicator
!! [dipquad] = if 1, atmfrc has been build without dipole-quadrupole part
!! [quadquad] = if 1, atmfrc has been build without quadrupole-quadrupole part
!!
!! OUTPUT
!!  phfrq(3*natom) = Phonon frequencies in Hartree
!!  displ_cart(2,3,natom,3*natom) = Phonon displacement in Cartesian coordinates
!!  [out_d2cart(2,3,3*natom,3,3*natom)] = The (interpolated) dynamical matrix for this q-point
!!  [out_eigvec(2*3*natom*3*natom) = The (interpolated) eigenvectors of the dynamical matrix in Cartesian coords..
!!  [out_displ_red(2*3*natom*3*natom) = The (interpolated) displacement in reduced coordinates.
!!  [dwdq(3,3*natom)] = Group velocities i.e. d(omega(q))/dq in Cartesian coordinates.
!!
!! PARENTS
!!      get_nv_fs_en,get_tau_k,harmonic_thermo,interpolate_gkk,m_gruneisen
!!      m_ifc,m_phgamma,m_phonons,m_phpi,m_sigmaph,m_tdep_phdos,mka2f,mka2f_tr
!!      mka2f_tr_lova,mkph_linwid,read_gkk
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_fourq(ifc, crystal, qpt, phfrq, displ_cart, &
                     nanaqdir, comm, &                              ! Optional [in]
                     out_d2cart, out_eigvec, out_displ_red, dwdq)   ! Optional [out]

!Arguments ------------------------------------
!scalars
 character(len=*),optional,intent(in) :: nanaqdir
 class(ifc_type),intent(in) :: Ifc
 type(crystal_t),intent(in) :: Crystal
 integer,optional,intent(in) :: comm
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: displ_cart(2,3,Crystal%natom,3*Crystal%natom)
 real(dp),intent(out) :: phfrq(3*Crystal%natom)
 real(dp),optional,intent(out) :: out_d2cart(2,3,Crystal%natom,3,Crystal%natom)
 real(dp),optional,intent(out) :: out_eigvec(2,3,Crystal%natom,3*Crystal%natom)
 real(dp),optional,intent(out) :: out_displ_red(2,3,Crystal%natom,3*Crystal%natom)
 real(dp),optional,intent(out) :: dwdq(3,3*crystal%natom)

!Local variables-------------------------------
!scalars
 integer :: natom, comm_
 real(dp) :: qphnrm
!arrays
 real(dp) :: my_qpt(3),eigvec(2,3,Crystal%natom,3*Crystal%natom),eigval(3*Crystal%natom)
 real(dp) :: d2cart(2,3,Ifc%mpert,3,Ifc%mpert),tsec(2)

! ************************************************************************

 ! Keep track of total time spent.
 call timab(1748, 1, tsec)

 natom = Crystal%natom
 ! TODO: Rewrite and Parallelize ewald9 in gtdyn9
 comm_ = xmpi_comm_self; if (present(comm)) comm_ = comm

 ! Use my_qpt because dfpt_phfrq can change the q-point (very bad design)
 qphnrm = one; my_qpt = qpt

 if (present(nanaqdir)) then
   ! This will break backward compatibility because qpt is **always** in reduced coordinates.
   ! while dfpt_phfrq assume cartesian coordinates !!!!!!!!!!!
   ! It does not make sense to change API just to treat this particular case
   ! We should **alwayse use q-points in reduced coordinates.
   qphnrm = zero
   select case (nanaqdir)
   case ("reduced")
     ! Convert to Cartesian.
     my_qpt = matmul(Crystal%gprimd, qpt)
   case ("cart")
     continue
   case default
     MSG_ERROR(sjoin("Wrong value for nanaqdir:", nanaqdir))
   end select
 end if

 ! The dynamical matrix d2cart is calculated here:
#ifdef MR_DEV
 call gtdyn9(Ifc%acell,Ifc%atmfrc,Ifc%dielt,Ifc%dipdip,Ifc%dyewq0,d2cart,Crystal%gmet,Ifc%gprim,Ifc%mpert,natom,&
   Ifc%nrpt,qphnrm,my_qpt,Crystal%rmet,Ifc%rprim,Ifc%rpt,Ifc%trans,Crystal%ucvol,Ifc%wghatm,Crystal%xred,Ifc%zeff,&
   Ifc%qdrp_cart,Ifc%ewald_option,comm_,dipquad=Ifc%dipquad,quadquad=Ifc%quadquad)
#else
 call gtdyn9(Ifc%acell,Ifc%atmfrc,Ifc%dielt,Ifc%dipdip,Ifc%dyewq0,d2cart,Crystal%gmet,Ifc%gprim,Ifc%mpert,natom,&
   Ifc%nrpt,qphnrm,my_qpt,Crystal%rmet,Ifc%rprim,Ifc%rpt,Ifc%trans,Crystal%ucvol,Ifc%wghatm,Crystal%xred,Ifc%zeff,&
   Ifc%qdrp_cart,Ifc%ewald_option,comm_)
#endif

 ! Calculate the eigenvectors and eigenvalues of the dynamical matrix
 call dfpt_phfrq(Ifc%amu,displ_cart,d2cart,eigval,eigvec,Crystal%indsym,&
   Ifc%mpert,Crystal%nsym,natom,Crystal%nsym,Crystal%ntypat,phfrq,qphnrm,my_qpt,&
   Crystal%rprimd,Ifc%symdynmat,Crystal%symrel,Crystal%symafm,Crystal%typat,Crystal%ucvol)

 ! OmegaSRLR: Perform decomposition of dynamical matrix
 !if (srlr==1) call omega_decomp(amu,natom,ntypat,typat,dynmatfull,dynmatsr,dynmatlr,iqpt,nqpt,eigvec)

 ! Return the interpolated dynamical matrix and the eigenvector for this q-point
 if (present(out_d2cart)) out_d2cart = d2cart(:,:,:natom,:,:natom)
 if (present(out_eigvec)) out_eigvec = eigvec

 ! Return phonon displacement in reduced coordinates.
 if (present(out_displ_red)) call phdispl_cart2red(natom, crystal%gprimd, displ_cart, out_displ_red)

 ! Option to get vectors in reduced coordinates?
 !call phdispl_cart2red(natom, crystal%gprimd, out_eigvec, out_eigvec_red)

 ! Compute group velocities.
 if (present(dwdq)) call ifc_get_dwdq(ifc, crystal, my_qpt, phfrq, eigvec, dwdq, comm_)

 call timab(1748, 2, tsec)

end subroutine ifc_fourq
!!***

!!****f* m_ifc/ifc_get_dwdq
!! NAME
!!  ifc_get_dwdq
!!
!! FUNCTION
!!  Compute phonon group velocities at an arbitrary q-point.
!!
!! INPUTS
!!  ifc<ifc_type>=Object containing the dynamical matrix and the IFCs.
!!  crystal<crystal_t> = Information on the crystalline structure.
!!  qpt(3)=q-point in reduced coordinates.
!!  eigvec(2*3*natom*3*natom) = The eigenvectors of the dynamical matrix.
!!  comm: MPI communicator
!!
!! OUTPUT
!!  dwdq(3,3*natom) = Group velocities e.g. d(omega(q))/dq in Cartesian coordinates.
!!
!! NOTES
!!  Using:
!!
!!    D(q) u(q,nu) = w(q, nu)**2 and <u(q,nu) | u(q,nu')> = \delta_{nu, nu'}
!!
!!  one can show, using the Hellman-Feynman theorem, that:
!!
!!    \nabla_q w(q, nu) = 1/(2 w(q, nu))  <u(q, nu)| \nabla_q D(q) | u(q, nu)>
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_get_dwdq(ifc, cryst, qpt, phfrq, eigvec, dwdq, comm)

!Arguments ------------------------------------
!scalars
 type(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: comm
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(in) :: phfrq(3*cryst%natom)
 real(dp),intent(in) :: eigvec(2,3*cryst%natom,3*cryst%natom)
 real(dp),intent(out) :: dwdq(3,3*cryst%natom)

!Local variables-------------------------------
!scalars
 !integer,save :: enough=0
 integer,parameter :: nqpt1=1,option2=2,sumg0=0
 integer :: ii,nu,natom3,jj
 real(dp) :: hh
!arrays
 real(dp) :: dddq(2,3*cryst%natom,3*cryst%natom,3),dot(2),qfd(3)
 real(dp) :: omat(2,3*cryst%natom,3*cryst%natom)
 real(dp) :: dyew(2,3*cryst%natom,3*cryst%natom)

! ************************************************************************

 ABI_UNUSED((/comm/))

 natom3 = cryst%natom * 3

 ! Generate the analytical part from the interatomic forces
 call dynmat_dq(qpt, cryst%natom, ifc%gprim, ifc%nrpt, ifc%rpt, ifc%atmfrc, ifc%wghatm, dddq)

 ! The analytical dynamical matrix dq has been generated
 ! in the normalized canonical coordinate system. Now, the
 ! phase is modified, in order to recover the usual (xred) coordinate of atoms.
 do ii=1,3
   call dymfz9(dddq(:,:,:,ii), cryst%natom, nqpt1, ifc%gprim, option2, qpt, ifc%trans)
   dddq(:,:,:,ii) = dddq(:,:,:,ii) * ifc%acell(ii)
 end do

 if (ifc%dipdip == 1) then
   ! Add the gradient of the non-analytical part.
   ! Note that we dddq is in cartesian cordinates.
   ! For the time being, the gradient is computed with finite difference and step hh.
   ! TODO: should generalize ewald9 to compute dq.
   !enough = enough + 1
   !if (enough <= 5)  MSG_WARNING("phonon velocities with dipdip==1 not yet tested.")
   hh = 0.01_dp
   do ii=1,3
     do jj=-1,1,2
       ! qcart --> qred
       qfd = zero; qfd(ii) = jj
       qfd = matmul(cryst%rprimd, qfd); qfd = qfd / normv(qfd, cryst%gmet, "G")
       !write(std_out,*)"normv:",normv(qfd, cryst%gmet, "G")
       qfd = qpt + hh * qfd

       call ewald9(ifc%acell,ifc%dielt,dyew,cryst%gmet,ifc%gprim,cryst%natom,qfd,&
          cryst%rmet,ifc%rprim,sumg0,cryst%ucvol,cryst%xred,ifc%zeff,ifc%qdrp_cart,&
          ifc%ewald_option,dipquad=ifc%dipquad,quadquad=ifc%quadquad)
       call q0dy3_apply(cryst%natom,ifc%dyewq0,dyew)
       dddq(:,:,:,ii) = dddq(:,:,:,ii) + (jj * half / hh) * dyew
     end do
   end do
 end if

 do ii=1,3
   call massmult_and_breaksym(cryst%natom, cryst%ntypat, cryst%typat, ifc%amu, dddq(:,:,:,ii))
 end do

 ! Compute 1/(2w(q)) <u(q)|dD(q)/dq|u(q)>
 do ii=1,3
   call zgemm('N','N',natom3,natom3,natom3,cone,dddq(:,:,:,ii),natom3,eigvec,natom3,czero,omat,natom3)
   do nu=1,natom3
     if (abs(phfrq(nu)) > tol12) then
       dot = cg_zdotc(natom3, eigvec(1,1,nu), omat(1,1,nu))
       ! abs(w) is needed to get the correct derivative if we have a purely imaginary solution.
       dwdq(ii, nu) = dot(1) / (two * abs(phfrq(nu)))
     else
       dwdq(ii, nu) = zero
     end if
   end do
 end do

end subroutine ifc_get_dwdq
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_speedofsound
!!
!! NAME
!! ifc_speedofsound
!!
!! FUNCTION
!!  Calculate the speed of sound by averaging the phonon group velocities of the
!!  three acoustic modes on a small sphere of radius qrad centered around Gamma.
!!  Perform spherical integration with Lebedev-Laikov grids
!!
!! INPUTS
!! ifc<ifc_type>=Object containing the dynamical matrix and the IFCs.
!! crystal<crystal_t> = Information on the crystalline structure.
!! qrad_tolkms(2):
!!   qrad=Radius of the sphere in reciprocal space
!!   atols_kms=Absolute tolerance in kilometer/second. The code generates spherical meshes
!!     until the results are converged twice within atols_kms.
!! ncid=the id of the open NetCDF file. Use nctk_noid to disable netcdf output.
!! comm=MPI communicator.
!!
!! OUTPUT
!!
!! PARENTS
!!      anaddb,m_gruneisen
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_speedofsound(ifc, crystal, qrad_tolkms, ncid, comm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: comm,ncid
 class(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: crystal
!arrays
 real(dp),intent(in) :: qrad_tolkms(2)

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: ii,nu,igrid,my_rank,nprocs,ierr,converged,npts,num_negw,vs_ierr,ncerr
 integer :: iatom,iatref,num_acoustic,isacoustic
 real(dp) :: min_negw,cpu,wall,gflops
 real(dp) :: qrad,tolkms,diff
 character(len=500) :: msg
 type(lebedev_t) :: lgrid
!arrays
 integer :: asnu(3)
 real(dp) :: qred(3),qvers_cart(3),qvers_red(3),quad(3),prev_quad(3),vs(7,3)
 real(dp) :: phfrqs(3*crystal%natom),dwdq(3,3*crystal%natom)
 real(dp) :: displ_cart(2,3*crystal%natom,3*crystal%natom),eigvec(2,3*crystal%natom,3*crystal%natom)

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 if (ifc%asr == 0) MSG_WARNING("Computing speed of sound with asr == 0! Use asr > 0")
 qrad = qrad_tolkms(1); tolkms = qrad_tolkms(2)
 ABI_CHECK(qrad > zero, "vs_qrad <= 0")
 ABI_CHECK(tolkms > zero, "vs_tolkms <= 0")

 call cwtime(cpu, wall, gflops, "start")

 ! Find the index of the first acoustic modes (needed to handle systems with unstable branches at Gamma
 ! In this case, indeed, we end up computing derivatives for the wrong branches, the integration becomes
 ! unstable and difficult to converge.
 qred = zero
 call ifc_fourq(ifc, crystal, qred, phfrqs, displ_cart,out_eigvec=eigvec)
 !write(std_out,*)"omega(q==Gamma): ",phfrqs

 num_acoustic = 0
 do nu = 1, 3*crystal%natom
   ! Check if this mode is acoustic like: scalar product of all displacement vectors are collinear
   isacoustic = 1
   ! Find reference atom with non-zero displacement
   do iatom=1,crystal%natom
     if(sum(displ_cart(:,(iatom-1)*3+1:(iatom-1)*3+3,nu)**2) >tol16)iatref=iatom
   end do
   ! Now compute scalar product, and check they are all positive
   do iatom = 1, crystal%natom
     if (sum(eigvec(:,(iatom-1)*3+1:(iatom-1)*3+3, nu)*eigvec(:,(iatref-1)*3+1:(iatref-1)*3+3, nu)) < tol16 ) isacoustic = 0
   end do
   if (isacoustic == 1) then
     num_acoustic=num_acoustic+1
     asnu(num_acoustic)=nu
     if (num_acoustic==3) exit
   end if
 end do

 ABI_CHECK(num_acoustic == 3, sjoin("Wrong number of acoustic modes:", itoa(num_acoustic)))

 write(std_out,"(a,3i2,a)") "The bands with indices ",asnu(:)," will be used to calculate the sound velocities"

 ! Speed of sound along reduced directions.
 do ii=1,6
   qred = zero; qred(MOD(ii-1,3)+1) = one
   if (ii >= 4 .and. ii <= 6) qred = matmul(crystal%rprimd, qred) ! Cartesian directions.
   !if (ii >= 7 .and. ii <= 9) qred = matmul(crystal%rprimd, qred) ! Cartesian directions.
   qvers_red = (qred / normv(qred, crystal%gmet, "G"))
   qred = qrad * qvers_red
   !write(std_out,*)"dir",normv(qred, crystal%gmet, "G"), qrad
   call ifc_fourq(ifc, crystal, qred, phfrqs, displ_cart, dwdq=dwdq)

   do nu=1,3
     vs(ii, nu) = sqrt(sum(dwdq(1:3,asnu(nu)) ** 2)) * Bohr_meter * 0.001_dp / Time_Sec
   end do
   write(std_out,"(a,3es12.4,a)")" ||vs(nu)||:",vs(ii,:), " [km/s]"

   qvers_cart = matmul(crystal%gprimd, qvers_red) * two_pi
   do nu=1,3
     vs(ii, nu) = dot_product(dwdq(1:3,asnu(nu)), qvers_cart) * Bohr_meter * 0.001_dp / Time_Sec
   end do
   write(std_out,"(a,3es12.4,a)")" <q|vs(nu)>:",vs(ii,:), " [km/s]"

   !do nu=1,3
   !  write(std_out,"(a,3es12.4,a)")" vs(nu)_vect_red:",&
   !     matmul(crystal%gprimd, dwdq(1:3,asnu(nu))) * Bohr_meter * 0.001_dp / Time_Sec, " [km/s]"
   !end do
 end do

 ! Spherical average with Lebedev-Laikov grids.
 converged = 0
 do igrid=1,lebedev_ngrids
   lgrid = lebedev_new(igrid)
   npts = lgrid%npts; quad = zero; num_negw = 0; min_negw = zero
   do ii=1,npts
     if (mod(ii, nprocs) /= my_rank) cycle ! mpi-parallelism

     ! Build q-point on sphere of radius qrad. qcart --> qred
     qred = matmul(crystal%rprimd, lgrid%versors(:, ii))
     qred = qrad * (qred / normv(qred, crystal%gmet, "G"))
     !write(std_out,*)"lebe",normv(qred, crystal%gmet, "G"), qrad
     call ifc_fourq(ifc, crystal, qred, phfrqs, displ_cart, dwdq=dwdq)
     if (any(phfrqs(asnu) < -tol8)) then
       num_negw = num_negw + 1; min_negw = min(min_negw, minval(phfrqs(asnu)))
     end if

     do nu=1,3
       quad(nu) = quad(nu) + lgrid%weights(ii) * sqrt(sum(dwdq(1:3,asnu(nu)) ** 2))
       !quad(nu) = quad(nu) + lgrid%weights(ii) * abs(dot_product(lgrid%versors(:,ii), dwdq(:,asnu(nu))))
     end do
   end do

   ! Will use km/sec unit for echo purposes
   quad = quad * Bohr_meter * 0.001_dp / Time_Sec
   call xmpi_sum(quad, comm, ierr)
   call xmpi_sum(num_negw, comm, ierr)
   call lebedev_free(lgrid)

   write(std_out,'(2(a,i6),a,3es12.4,a,es12.4)') &
     " Lebedev-Laikov grid: ",igrid,", npts: ", npts, " vs_sphavg(ac_modes): ",quad, " <vs>: ",sum(quad)/3

   if (igrid > 1) then
     diff = zero
     do nu=1,3
       diff = diff + abs(quad(nu) - prev_quad(nu)) / 3
     end do
     !if (abs(sum(quad - prev_quad)/3) < tolkms) then
     if (diff < tolkms) then
        converged = converged + 1
     else
        converged = 0
     end if
   end if
   prev_quad = quad
   vs(7, :) = quad
   if (converged == 2) exit
 end do ! igrid

 if (my_rank == master) then
   ! vs_err: 1 if not converged, < 0 if negative freqs, == 0 if success.
   vs_ierr = 0
   do ii=1,3
     write(ab_out,"(a,3es12.4,a,i1)")" Speed of sound:",vs(ii,:)," [km/s] along reduced direction: ",ii
   end do
   write(ab_out,'(2(a,es12.4),a,i0)') &
     " Lebedev-Laikov integration with qradius: ", qrad, " tolkms: ",tolkms, " [km/s], npts: ", npts
   write(ab_out,"(a,3es12.4,a,es12.4)")" Spherical average:",vs(7,:)," [km/s], ",sum(vs(7,:))/3
   if (converged /= 2) then
     vs_ierr = 1
     write(msg,'(a,es12.4,a)')" WARNING: Results are not converged within: ",tolkms, " [km/s]"
     call wrtout(ab_out, msg)
     MSG_WARNING(msg)
   end if
   if (num_negw > 0) then
     vs_ierr = -num_negw
     write(msg,'(a,i0,a,es12.4,3a)') &
       " WARNING: Detected ",num_negw, " negative frequencies. Minimum was: ",min_negw * Ha_meV, "[meV]",ch10,&
       " Speed of sound could be wrong"
     call wrtout(ab_out, msg)
     MSG_WARNING(msg)
   end if

   ! Dump results to netcdf file.
   if (ncid /= nctk_noid) then
#ifdef HAVE_NETCDF
     ncerr = nctk_def_arrays(ncid, [nctkarr_t("vsound", "dp", "seven, three")], defmode=.True.)
     NCF_CHECK(ncerr)
     ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "vsound_ierr"])
     NCF_CHECK(ncerr)
     ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "vsound_qrad", "vsound_tolkms"])
     NCF_CHECK(ncerr)
     ncerr = nctk_def_arrays(ncid, [nctkarr_t("asnu", "i", "three")], defmode=.True.)
     NCF_CHECK(ncerr)
     ! Write data.
     NCF_CHECK(nctk_set_datamode(ncid))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vsound_ierr"), vs_ierr))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vsound_qrad"), qrad))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vsound_tolkms"), tolkms))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vsound"), vs))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "asnu"), asnu))
#endif
   end if
 end if

 call cwtime_report(" ifc_speedofsound", cpu, wall, gflops)

end subroutine ifc_speedofsound
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_autocutoff
!! NAME
!!  ifc_autocutoff
!!
!! FUNCTION
!! Find the value of nsphere that gives non-negative frequencies around Gamma
!! in a small sphere of radius qrad.
!! Use bisection to reduce the number of attemps although there's no guarantee
!! that the number of negative frequencies is monotonic.
!!
!! INPUTS
!!  crystal<crystal_t> = Information on the crystalline structure.
!!  comm=MPI communicator
!!
!! SIDE EFFECTS
!!  ifc%wghatm(natom,natom,nrpt) = Weights associated to a pair of atoms and to a R vector
!!    with the last cutoff found by the bisection algorithm applied.
!!  ifc%atmfrc(2,3,natom,3,natom,nrpt)= ASR-imposed Interatomic Forces
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_autocutoff(ifc, crystal, comm)

!Arguments ------------------------------------
!scalars
 type(ifc_type),intent(inout) :: ifc
 type(crystal_t),intent(in) :: crystal
 integer,intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iq_ibz,ierr,my_rank,nprocs,ii,nsphere,num_negw,jl,ju,jm,natom,nrpt
 real(dp),parameter :: rifcsph0=zero
 real(dp) :: adiff,qrad,min_negw,xval,rcut_min
 type(lebedev_t) :: lgrid
!arrays
 real(dp) :: displ_cart(2*3*ifc%natom*3*ifc%natom)
 real(dp) :: qred(3),qred_vers(3),phfrqs(3*ifc%natom) !,dwdq(3,3*ifc%natom)
 real(dp),allocatable :: ref_phfrq(:,:),cut_phfrq(:,:)
 real(dp),allocatable :: save_wghatm(:,:,:),save_atmfrc(:,:,:,:,:)

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 natom = ifc%natom; nrpt = ifc%nrpt

 ! Compute frequencies on the ab-initio q-mesh without cutoff.
 ABI_CALLOC(ref_phfrq, (3*natom, ifc%nqibz))
 do iq_ibz=1,ifc%nqibz
   if (mod(iq_ibz, nprocs) /= my_rank) cycle ! mpi-parallelism
   call ifc_fourq(ifc, crystal, ifc%qibz(:,iq_ibz), ref_phfrq(:,iq_ibz), displ_cart)
 end do
 call xmpi_sum(ref_phfrq, comm, ierr)

 ABI_MALLOC(save_wghatm, (natom,natom,nrpt))
 ABI_MALLOC(save_atmfrc, (3,natom,3,natom,ifc%nrpt))
 save_wghatm = ifc%wghatm; save_atmfrc = ifc%atmfrc

 ABI_MALLOC(cut_phfrq, (3*natom, ifc%nqibz))
 qrad = 0.01; lgrid = lebedev_new(16)

 if (my_rank == master) then
   write(ab_out, "(a)")" Apply cutoff on IFCs. Using bisection algorithm to find initial guess for nsphere."
   write(ab_out, "(a,i0)")" Maximum nuber of atom-centered spheres: ",natom * nrpt
   write(ab_out, "(a,i0,a,f5.3)")" Using Lebedev-Laikov grid with npts: ",lgrid%npts, ", qrad: ",qrad
   write(ab_out, "(/,a)")" <adiff>: Average difference between ab-initio frequencies and frequencies with cutoff."
   write(ab_out, "(a)")" num_negw: Number of negative freqs detected in small sphere around Gamma."
   write(ab_out, "(a)")" min_negw: Min negative frequency on the small sphere."
   write(ab_out, "(a,/,/)")" rifcsph: Effective cutoff radius corresponding to nsphere."
   write(ab_out, "(a)")" nsphere   <adiff>[meV]   num_negw   min_negw[meV]   rifcsph"
 end if

 jl = 0; ju = natom * nrpt + 1 ! Initialize lower and upper limits.
 do
   if (ju - jl <= 1) then
     exit
   end if
   jm = (ju + jl) / 2  ! Compute a midpoint
   nsphere = jm

   ifc%wghatm = save_wghatm; ifc%atmfrc = save_atmfrc
   call corsifc9(ifc%acell,ifc%gprim, natom, nrpt,nsphere,rifcsph0,ifc%rcan,ifc%rprim,ifc%rpt,rcut_min,ifc%wghatm)
   if (ifc%asr > 0) call asrif9(ifc%asr,ifc%atmfrc,ifc%natom,ifc%nrpt,ifc%rpt,ifc%wghatm)

   cut_phfrq = zero
   do iq_ibz=1,ifc%nqibz
     if (mod(iq_ibz, nprocs) /= my_rank) cycle ! mpi-parallelism
     call ifc_fourq(ifc, crystal, ifc%qibz(:,iq_ibz), cut_phfrq(:,iq_ibz), displ_cart)
     !write(std_out,*)cut_phfrq(1,iq_ibz),ref_phfrq(1,iq_ibz)
   end do
   call xmpi_sum(cut_phfrq, comm, ierr)

   ! Test wether there are negative frequencies around gamma, including reciprocal lattice vectors.
   num_negw = 0; min_negw = zero
   do ii=1,lgrid%npts+3
     if (mod(ii, nprocs) /= my_rank) cycle ! mpi-parallelism
     if (ii <= 3) then
       qred = zero; qred(ii) = one
     else
       qred = matmul(crystal%rprimd, lgrid%versors(:, ii-3))
     end if
     qred_vers = (qred / normv(qred, crystal%gmet, "G"))
     qred = qrad * qred_vers
     call ifc_fourq(ifc, crystal, qred, phfrqs, displ_cart) !, dwdq=dwdq)
     if (any(phfrqs < +tol8)) then
       num_negw = num_negw + 1; min_negw = min(min_negw, minval(phfrqs))
     end if
     !do jj=1,3
     !  xval = dot_product(dwdq(:,jj), matmul(crystal%gprimd, qred_vers))
     !  if (xval < zero) num_negw = num_negw + 1
     !end do
   end do
   call xmpi_sum(num_negw, comm, ierr)
   xval = min_negw; call xmpi_min(xval, min_negw, comm, ierr)

   adiff = sum(abs(cut_phfrq - ref_phfrq)) / (ifc%nqibz * 3 * natom)
   if (my_rank == master) then
     write(ab_out,"(a,i7,1x,es13.4,4x,i8,1x,es13.4,2x,es13.4)") &
       "-",nsphere, adiff * Ha_meV, num_negw, min_negw * Ha_meV, rcut_min
   end if

   if (num_negw == 0) then
     jl = jm ! Replace lower limit
   else
     ju = jm ! Replace upper limit
   end if
 end do

 ABI_FREE(ref_phfrq)
 ABI_FREE(cut_phfrq)
 ABI_FREE(save_wghatm)
 ABI_FREE(save_atmfrc)
 call lebedev_free(lgrid)

end subroutine ifc_autocutoff
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/corsifc9
!! NAME
!! corsifc9
!!
!! FUNCTION
!! Applies a cutoff on the ifc in real space
!!
!! INPUTS
!! acell(3)=length scales by which rprim is to be multiplied
!! gprim(3,3)=dimensionless primitive translations in reciprocal space
!! natom=number of atoms in unit cell
!! nrpt= Number of R points in the Big Box
!! rcan(3,natom)=canonical coordinates of atoms
!! rprim(3,3)=dimensionless primitive translations in real space
!! rpt(3,nrpt)=canonical coordinates of the points in the BigBox.
!! nsphere=number of atoms to be included in the cut-off sphere for interatomic
!!  force constant; if = 0 : maximum extent allowed by the grid.
!! rifcsph=radius for cutoff of IFC
!! wghatm(natom,natom,nrpt) = Weights associated to a pair of atoms and to a R vector.
!!
!! OUTPUT
!! wghatm(natom,natom,nrpt) = Weights associated to a pair of atoms and to a R vector
!!  with the required cutoff applied.
!! rcut_min=Effective cutoff. Defined by the minimum cutoff radius over the natom sites.
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine corsifc9(acell,gprim,natom,nrpt,nsphere,rifcsph,rcan,rprim,rpt,rcut_min,wghatm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nrpt,nsphere
 real(dp),intent(in) :: rifcsph
 real(dp),intent(out) :: rcut_min
!arrays
 real(dp),intent(in) :: acell(3)
 real(dp),intent(in) :: gprim(3,3),rcan(3,natom)
 real(dp),intent(in) :: rprim(3,3),rpt(3,nrpt)
 real(dp),intent(inout) :: wghatm(natom,natom,nrpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,ii,index,irpt
 real(dp) :: rmax,rsigma,r0
!arrays
 integer,allocatable :: list(:)
 real(dp),allocatable :: dist(:,:,:),wkdist(:)

! *********************************************************************

 ! Compute the distances between atoms
 ! dist(ia,ib,irpt) contains the distance from atom ia to atom ib in unit cell irpt.
 ABI_MALLOC(dist,(natom,natom,nrpt))
 call dist9(acell,dist,gprim,natom,nrpt,rcan,rprim,rpt)

 ABI_MALLOC(list,(natom*nrpt))
 ABI_MALLOC(wkdist,(natom*nrpt))

 ! loop on all generic atoms.
 rcut_min = huge(one)
 do ia=1,natom

   wkdist = reshape(dist(ia,:,:), [natom*nrpt])
   do ii=1,natom*nrpt
     list(ii)=ii
   end do
   ! This sorting algorithm is slow ...
   call sort_dp(natom*nrpt,wkdist,list,tol14)
   rmax = wkdist(natom*nrpt)

   ! zero the outside IFCs: act on wghatm

   ! fix number of spheres
   if(nsphere/=0.and.nsphere<natom*nrpt)then
     rcut_min = min(rcut_min, wkdist(nsphere+1))
     do ii=nsphere+1,natom*nrpt
       index=list(ii)
       irpt=(index-1)/natom+1
       ib=index-natom*(irpt-1)
       wghatm(ia,ib,irpt)=zero
     end do
   end if

   ! or fix radius of maximum ifc
   if(rifcsph>tol10)then
     do ii=nsphere+1,natom*nrpt
       index=list(ii)
       ! preserve weights for atoms inside sphere of radius rifcsph
       if (wkdist(ii) < rifcsph) cycle
       rcut_min = min(rcut_min, wkdist(ii))
       irpt=(index-1)/natom+1
       ib=index-natom*(irpt-1)
       wghatm(ia,ib,irpt)=zero
     end do
   end if

   ! filter smoothly to 0 at edge of WS cell
   if (rifcsph < -tol10) then
     ! Use different filter
     r0 = abs(rifcsph) * rmax; rsigma = half*(rmax-r0) !one
     rcut_min = r0 ! Set it to r0
     do ii=nsphere+1,natom*nrpt
       index=list(ii)
       irpt=(index-1)/natom+1
       ib=index-natom*(irpt-1)
       wghatm(ia,ib,irpt) = wghatm(ia,ib,irpt) * half * abi_derfc((wkdist(ii) - r0) / rsigma)
     end do
   end if

 end do

 ABI_FREE(dist)
 ABI_FREE(list)
 ABI_FREE(wkdist)

end subroutine corsifc9
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_write
!!
!! NAME
!! ifc_write
!!
!! FUNCTION
!! Adds the real-space interatomic force constants to:
!!  the output file,
!!  a NetCDF file which is already open on ncid
!!  if prt_ifc==1, to the ifcinfo.out file
!!  to a TDEP file named outfile.forceconstants_ABINIT
!!
!! INPUTS
!! Ifc<type(ifc_type)>=Object containing the dynamical matrix and the IFCs.
!! ifcana= 0 => no analysis of ifc ; 1 => full analysis
!! atifc(natom) =  atifc(ia) equals 1 if the analysis of ifc has to be done for atom ia; otherwise 0.
!! ifcout= Number of interatomic force constants written in the output file
!! prt_ifc = flag to print out ifc information for dynamical matrix (AI2PS)
!! ncid=the id of the open NetCDF file. Set to nctk_noid if netcdf output is not wanted.
!!
!! OUTPUT
!!   written in the output file and in the NetCDF file
!!
!! NOTES
!! This routine should be executed by one processor only
!!
!! TODO:
!!  1) ifc_write should not have side effects
!!
!!  2) the code is unreadable and horrible - 3/4 different file formats for the
!!  same stuff. We should make different subroutines, even if it duplicates some code
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_write(Ifc,ifcana,atifc,ifcout,prt_ifc,ncid)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ifcout,ifcana,prt_ifc,ncid
 class(ifc_type),intent(inout) :: Ifc
!arrays
 integer,intent(in) :: atifc(Ifc%natom)

!Local variables -------------------------
!scalars
 integer :: ia,ib,ii,ncerr,iatifc,ifcout1,mu,nu,iout, irpt
! unit number to print out ifc information for dynamical matrix (AI2PS)
 integer :: unit_ifc, unit_tdep
 real(dp) :: detdlt
 real(dp) :: maxdist_tdep
 character(len=500) :: message
 character(len=4) :: str1, str2
!arrays
 integer,allocatable :: list(:),indngb(:)
 real(dp) :: invdlt(3,3),ra(3),xred(3),dielt(3,3)
 real(dp),allocatable :: dist(:,:,:),wkdist(:),rsiaf(:,:,:),sriaf(:,:,:),vect(:,:,:)
 real(dp),allocatable :: posngb(:,:)
 real(dp) :: gprimd(3,3),rprimd(3,3)

! *********************************************************************

 iout = ab_out
 dielt = ifc%dielt

 ! Compute the distances between atoms
 ABI_MALLOC(dist,(Ifc%natom,Ifc%natom,Ifc%nrpt))
 call dist9(Ifc%acell,dist,Ifc%gprim,Ifc%natom,Ifc%nrpt,Ifc%rcan,Ifc%rprim,Ifc%rpt)
 ! Now dist(ia,ib,irpt) contains the distance from atom ia to atom ib in unit cell irpt.

 ABI_MALLOC(list,(Ifc%natom*Ifc%nrpt))
 ABI_MALLOC(wkdist,(Ifc%natom*Ifc%nrpt))

 ! Calculating the inverse (transpose) of the dielectric tensor
 call matr3inv(dielt,invdlt)

 ! Calculating the determinant of the dielectric tensor
 detdlt=dielt(1,1)*dielt(2,2)*dielt(3,3)+dielt(1,3)*dielt(2,1)*&
& dielt(3,2)+dielt(1,2)*dielt(2,3)*dielt(3,1)-dielt(1,3)*&
& dielt(2,2)*dielt(3,1)-dielt(1,1)*dielt(2,3)*dielt(3,2)-&
& dielt(1,2)*dielt(2,1)*dielt(3,3)

! echo to log file
 write(std_out,'(a)' )' ifc_write: analysis of interatomic force constants '
 call mkrdim(Ifc%acell,Ifc%rprim,rprimd)
 call matr3inv(rprimd,gprimd)

 if (iout > 0) then
   write(iout, '(/,a,/)' )' Analysis of interatomic force constants '
#ifdef MR_DEV
   if(Ifc%dipdip==1.and.Ifc%dipquad==0.and.Ifc%quadquad==0)then
#else
   if(Ifc%dipdip==1)then
#endif 
     write(iout, '(a)' )' Are given : column(1-3), the total force constant'
     write(iout, '(a)' )'       then  column(4-6), the Ewald part'
     write(iout, '(a)' )'       then  column(7-9), the short-range part'
     write(iout, '(a)' )' Column 1, 4 and 7 are related to the displacement'
     write(iout, '(a)' )'       of the generic atom along x,               '
     write(iout, '(a)' )' column 2, 5 and 8 are related to the displacement'
     write(iout, '(a)' )'       of the generic atom along y,               '
     write(iout, '(a)' )' column 3, 6 and 9 are related to the displacement'
     write(iout, '(a)')'       of the generic atom along z.               '
#ifdef MR_DEV
   else if(Ifc%dipquad==1.or.Ifc%quadquad==1)then
     write(iout, '(a)' )' Are given : column(1-3), ONLY the short-range part!!!!'
     write(iout, '(a)' )' column 1 is related to the displacement'
     write(iout, '(a)' )'        of the generic atom along x,    '
     write(iout, '(a)' )' column 2 is related to the displacement'
     write(iout, '(a)' )'        of the generic atom along y,    '
     write(iout, '(a)' )' column 3 is related to the displacement'
     write(iout, '(a)' )'        of the generic atom along z,    '
#endif
   else if(Ifc%dipdip==0)then
     write(iout, '(a)' )' column 1 is related to the displacement'
     write(iout, '(a)' )'        of the generic atom along x,    '
     write(iout, '(a)' )' column 2 is related to the displacement'
     write(iout, '(a)' )'        of the generic atom along y,    '
     write(iout, '(a)' )' column 3 is related to the displacement'
     write(iout, '(a)' )'        of the generic atom along z,    '
   end if
 end if

 if (ifcout>Ifc%natom*Ifc%nrpt .or. ifcout == -1) then
   ifcout1=Ifc%natom*Ifc%nrpt
   write(message, '(3a,i0,a)' )&
&   'The value of ifcout exceeds the number of atoms in the big box.', ch10, &
&   'Output limited to ',Ifc%natom*Ifc%nrpt,' atoms.'
   MSG_WARNING(message)
 else
   ifcout1=ifcout
 end if

 ! set up file for real space ifc output, if required
 if (prt_ifc == 1) then
   if (open_file('ifcinfo.out', message, newunit=unit_ifc, status="replace") /= 0) then
     MSG_ERROR(message)
   end if
   write(iout, '(a,a)' )ch10,&
&   '  NOTE: Open file ifcinfo.out, for the output of interatomic force constants. This is because prt_ifc==1. '

   if (open_file('outfile.forceconstants_ABINIT', message, newunit=unit_tdep, status="replace") /= 0) then
     MSG_ERROR(message)
   end if
   write(iout, '(a,a,a)' )ch10,&
&   '  NOTE: Open file outfile.forceconstants_ABINIT, for the output of interatomic force',&
&   ' constants in TDEP format. This is because prt_ifc==1. '
   ! Print necessary stuff for TDEP
   write(unit_tdep,"(1X,I10,15X,'How many atoms per unit cell')") Ifc%natom

   ! look at all pairs, find furthest one with weight 1
!   do ia
!Ifc%wghatm(ia,ib,irpt)
   maxdist_tdep = Ifc%r_inscribed_sphere !maxval(dist)*0.8_dp
   write(unit_tdep,"(1X,F20.15,5X,'Realspace cutoff (A)')") maxdist_tdep*Bohr_Ang
 end if

#ifdef HAVE_NETCDF
 if (ncid /= nctk_noid) then
   ! initialize netcdf variables
   ncerr = nctk_def_dims(ncid, [nctkdim_t("natifc", SUM(atifc)), nctkdim_t("number_of_r_points_big_box", Ifc%nrpt), &
     nctkdim_t("number_of_atoms_big_box", Ifc%natom*Ifc%nrpt), nctkdim_t("ifcout", ifcout1)], defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t('ifc_atoms_indices', "i", "natifc"),&
     nctkarr_t('ifc_neighbours_indices', "i", "ifcout, natifc"),&
     nctkarr_t('ifc_distances', "dp", "ifcout, natifc "),&
     nctkarr_t('ifc_matrix_cart_coord', "dp", "number_of_cartesian_directions,number_of_cartesian_directions, ifcout, natifc")])
   NCF_CHECK(ncerr)

   if (Ifc%dipdip==1) then
     ncerr = nctk_def_arrays(ncid, [&
       nctkarr_t('ifc_matrix_cart_coord_short_range', "dp", &
       "number_of_cartesian_directions, number_of_cartesian_directions, ifcout, natifc")])
     NCF_CHECK(ncerr)
   end if

   if (ifcana==1) then
     ncerr = nctk_def_arrays(ncid, [&
       nctkarr_t('ifc_local_vectors', "dp", "number_of_cartesian_directions, number_of_cartesian_directions, ifcout, natifc")])
     NCF_CHECK(ncerr)
   end if

   NCF_CHECK(nctk_set_datamode(ncid))
 end if
#endif

 ABI_MALLOC(rsiaf,(3,3,ifcout1))
 ABI_MALLOC(sriaf,(3,3,ifcout1))
 ABI_MALLOC(vect,(3,3,ifcout1))
 ABI_MALLOC(indngb,(ifcout1))
 ABI_MALLOC(posngb,(3,ifcout1))

 iatifc=0

 ! BIG loop on all generic atoms
 do ia=1,Ifc%natom
   if(atifc(ia)==1)then

     iatifc=iatifc+1

     ! First transform canonical coordinates to reduced coordinates
     do ii=1,3
       xred(ii)=Ifc%gprim(1,ii)*Ifc%rcan(1,ia)+Ifc%gprim(2,ii)*Ifc%rcan(2,ia)+Ifc%gprim(3,ii)*Ifc%rcan(3,ia)
     end do

     ! Then to cartesian coordinates
     ra(:)=xred(1)*Ifc%acell(1)*Ifc%rprim(:,1)+ xred(2)*Ifc%acell(2)*Ifc%rprim(:,2)+ xred(3)*Ifc%acell(3)*Ifc%rprim(:,3)

     ! This sorting algorithm is slow ...
     wkdist(:)=reshape(dist(ia,:,:),(/Ifc%natom*Ifc%nrpt/))
     do ii=1,Ifc%natom*Ifc%nrpt
       list(ii)=ii
     end do
     call sort_dp(Ifc%natom*Ifc%nrpt,wkdist,list,tol14)

     if (iout > 0) then
       write(iout, '(a)' )
       write(std_out,'(a,i4)' )' generic atom number',ia
       write(iout, '(a,i4)' )' generic atom number',ia
       write(std_out,'(a,3es16.8)' ) ' with cartesian coordinates',ra(1:3)
       write(iout,'(a,3es16.8)' ) ' with cartesian coordinates',ra(1:3)
       write(iout, '(a)' )
     end if

     call ifc_getiaf(Ifc,ifcana,ifcout1,iout,ifc%zeff,ia,ra,list,dist,invdlt,&
&                    detdlt,rsiaf,sriaf,vect,indngb,posngb)

     if (prt_ifc == 1) then
       do ii=1,ifcout1
         if (wkdist(ii) > maxdist_tdep) exit
       end do
       ii = ii - 1
       write(unit_tdep,"(1X,I10,15X,'How many neighbours does atom ',I3,' have')") ii, ia

       do ii=1,ifcout1
         ib = indngb(ii)
         irpt = (list(ii)-1)/Ifc%natom+1
         ! limit printing to maximum distance for tdep
         if (wkdist(ii) > maxdist_tdep) cycle

         !TDEP
         call int2char4(ii, str1)
         call int2char4(ia, str2)
         write(unit_tdep,"(1X,I10,15X,a,a,a,a)") ib, &
&            'In the unit cell, what is the index of neighbour ', &
&            trim(str1), " of atom ", trim(str2)
         ! The lattice vector needs to be in reduced coordinates?
         ! TODO: check whether this is correct for TDEP: might need just lattice
         ! vector part and not full vector, and could be in integers instead of
         ! cartesian vector...
         write (unit_tdep,'(3es28.16)') matmul(Ifc%rpt(1:3,irpt),Ifc%gprim)

         !AI2PS
         write(unit_ifc,'(i6,i6)') ia,ii
         write(unit_ifc,'(3es28.16)') posngb(1:3,ii)
         do nu=1,3
           !TDEp
           ! And the actual short ranged forceconstant: TODO: check if
           ! a transpose is needed or a swap between the nu and the mu
           !write(unit_tdep,'(3f28.16)') (sriaf(nu,mu,ii)*Ha_eV/amu_emass, mu=1, 3)
           write(unit_tdep,'(3f28.16)') (Ifc%short_atmfrc(mu,ia,nu,ib,irpt)*Ha_eV/Bohr_Ang**2, mu=1, 3)

           !AI2PS
           write(unit_ifc,'(3f28.16)')(rsiaf(nu,mu,ii),mu=1,3)
         end do
       end do

#ifdef HAVE_NETCDF
       if (ncid /= nctk_noid) then
         NCF_CHECK(nf90_put_var(ncid, vid("ifc_atoms_indices"), ia, start=[iatifc]))
         NCF_CHECK(nf90_put_var(ncid, vid("ifc_neighbours_indices"), indngb, start=[1,iatifc], count=[ifcout1,1]))
         NCF_CHECK(nf90_put_var(ncid, vid("ifc_distances"), wkdist(:ifcout1), start=[1,iatifc],count=[ifcout1,1]))
         ncerr = nf90_put_var(ncid, vid("ifc_matrix_cart_coord"), rsiaf, start=[1,1,1,iatifc], count=[3,3,ifcout1,1])
         NCF_CHECK(ncerr)
         if (Ifc%dipdip==1) then
           ncerr = nf90_put_var(ncid, vid("ifc_matrix_cart_coord_short_range"), sriaf, &
             start=[1,1,1,iatifc], count=[3,3,ifcout1,1])
           NCF_CHECK(ncerr)
         end if
         if (ifcana==1) then
           ncerr = nf90_put_var(ncid, vid("ifc_local_vectors"), vect, start=[1,1,1,iatifc], count=[3,3,ifcout1,1])
           NCF_CHECK(ncerr)
         end if
       end if
#endif
     end if
   end if ! End the condition on atifc
 end do ! End Big loop on atoms in the unit cell, and corresponding test


! NB for future use: in TDEP the following can also be provided.
!        ! Print some auxiliary information, if it is there. Such as norm of
!        ! forceconstant per shell, which shells there are and so on.
!        if ( fc%npairshells .gt. 0 .and. allocated(fc%pairshell) ) then
!            write(u,"(1X,I10,15X,'Number of irreducible coordination shells')") fc%npairshells
!            do i=1,fc%npairshells
!                write(u,"(1X,I10,1X,F16.10,1X,F16.10,15X,'number atoms in shell, radius, norm of forceconstant',I0)") fc%pairshell(i)%n,fc%pairshell(i)%rad,fc%pairshell(i)%norm,i
!            enddo
!            do i=1,fc%npairshells
!                do j=1,fc%pairshell(i)%n
!                    write(u,"(1X,3(1X,F18.12),2(1X,I0))") lo_chop(matmul(p%inv_latticevectors,fc%pairshell(i)%vec(:,j)),lo_sqtol),fc%pairshell(i)%atind(j),fc%pairshell(i)%pairind(j)
!                enddo
!            enddo
!        endif

 if (prt_ifc == 1) then
   close(unit_ifc)
   close(unit_tdep)

   if (open_file('infile.lotosplitting_ABINIT', message, newunit=unit_tdep, status="replace") /= 0) then
     MSG_ERROR(message)
   end if
   write(unit_tdep,'(3es28.16)') dielt(:,1)
   write(unit_tdep,'(3es28.16)') dielt(:,2)
   write(unit_tdep,'(3es28.16)') dielt(:,3)
   do ia = 1, Ifc%natom
     do ii = 1, 3
       write(unit_tdep,'(3es28.16)') ifc%zeff(:,ii,ia)
     end do
   end do
   close(unit_tdep)
 end if

 ABI_FREE(rsiaf)
 ABI_FREE(sriaf)
 ABI_FREE(vect)
 ABI_FREE(indngb)
 ABI_FREE(posngb)
 ABI_FREE(dist)
 ABI_FREE(list)
 ABI_FREE(wkdist)

#ifdef HAVE_NETCDF
contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid
#endif

end subroutine ifc_write
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_getiaf
!!
!! NAME
!! ifc_getiaf
!!
!! FUNCTION
!! Extracts the IFCs needed for the output for one atom in the
!! unit cell. Accumulates the results for writing in the NetCDF file.
!! Prints to the output file
!!
!! INPUTS
!! Ifc<type(ifc_type)>=Object containing the dynamical matrix and the IFCs.
!! ifcana= 0 => no analysis of ifc ; 1 => full analysis
!! ifcout= Number of interatomic force constants written in the output file
!! iout=unit number for nice output
!! zeff(3,3,natom)=effective charge on each atom, versus electric field and atomic displacement
!! ia=index of the atom in the unit cell for which the IFCs are being analyzed
!! ra(3)=position of atom ia in cartesian coordinates
!! list(ifcout)=index of permutation for distances from atom ia in ascending order
!! dist(natom,natom,nrpt)=distance from atom ia to atom ib in unit cell irpt.
!! invdlt(3,3)=inverse (transpose) of the dielectric tensor
!! detdlt=determinant of the dielectric tensor
!!
!! OUTPUT
!! rsiaf(3,3,ifcout)=list of real space IFCs
!! sriaf(3,3,ifcout)=list of the short range part of the real space IFCs
!! vect(3,3,ifcout)=base vectors for local coordinates (longitudinal/transverse
!! indngb(ifcout)=indices in the unit cell of the neighbouring atoms
!! posngb(3,ifcout)=position of the neighbouring atoms in cartesian coordinates
!! output file
!!
!! NOTES
!! This routine should be executed by one processor only
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_getiaf(Ifc,ifcana,ifcout,iout,zeff,ia,ra,list,&
& dist,invdlt,detdlt,rsiaf,sriaf,vect,indngb,posngb)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ia,ifcana,ifcout,iout
 real(dp), intent(in) :: detdlt
 type(ifc_type),intent(inout) :: Ifc
!arrays
 real(dp),intent(in) :: invdlt(3,3),ra(3)
 real(dp),intent(in) :: dist(Ifc%natom,Ifc%natom,Ifc%nrpt)
 real(dp),intent(in) :: zeff(3,3,Ifc%natom)
 integer,intent(in) :: list(Ifc%natom*Ifc%nrpt)
 integer,intent(out) :: indngb(ifcout)
 real(dp),intent(out) :: rsiaf(3,3,ifcout),sriaf(3,3,ifcout),vect(3,3,ifcout),posngb(3,ifcout)

!Local variables -------------------------
!scalars
 integer :: flag,ib,ii,index,jj,kk,mu,nu,irpt
 real(dp) :: ew1,rsq,scprod,trace1,trace2,trace3
 real(dp) :: yy,dist1
 character(len=500) :: message
!arrays
 real(dp) :: ewiaf0(3,3),ewiaf1(3,3),ewloc(3,3),ifcloc(3,3)
 real(dp) :: rcart(3),rdiff(3),rsloc(3,3)
 real(dp) :: srloc(3,3),vect1(3),vect2(3),vect3(3),work(3),xx(3)

! *********************************************************************

 if(ifcana==1)then
   ! Generate the local coordinate system for the atom ia
   index=list(2)
   write(std_out,*)index
   call canct9(Ifc%acell,Ifc%gprim,ib,index,irpt,Ifc%natom,Ifc%nrpt,Ifc%rcan,rcart,Ifc%rprim,Ifc%rpt)
   vect2(1)=rcart(1)-ra(1)
   vect2(2)=rcart(2)-ra(2)
   vect2(3)=rcart(3)-ra(3)
   flag=0
   do ii=3,Ifc%natom*Ifc%nrpt
     index=list(ii)
     call canct9(Ifc%acell,Ifc%gprim,ib,index,irpt,Ifc%natom,Ifc%nrpt,Ifc%rcan,rcart,Ifc%rprim,Ifc%rpt)
     vect1(1)=(rcart(1)-ra(1))-vect2(1)
     vect1(2)=(rcart(2)-ra(2))-vect2(2)
     vect1(3)=(rcart(3)-ra(3))-vect2(3)
     scprod=0.0_dp
     do jj=1,3
       scprod=scprod+vect1(jj)**2
     end do
     do jj=1,3
       vect1(jj)=vect1(jj)/scprod**0.5
     end do
     scprod=0.0_dp
     do jj=1,3
       scprod=scprod+vect2(jj)*vect1(jj)
     end do
     do jj=1,3
       work(jj)=vect2(jj)-vect1(jj)*scprod
     end do
     scprod=0.0_dp
     do jj=1,3
       scprod=scprod+work(jj)**2
     end do
     if(scprod>1.0d-10)then
       flag=1
     end if
     if(flag==1)exit
   end do
   if(flag==0)then
     write(message, '(3a)' )&
&     'Unable to find a third atom not aligned',ch10,&
&     'with the two selected ones.'
     MSG_BUG(message)
   end if
   vect2(1)=work(1)/scprod**0.5
   vect2(2)=work(2)/scprod**0.5
   vect2(3)=work(3)/scprod**0.5
   vect3(1)=vect1(2)*vect2(3)-vect1(3)*vect2(2)
   vect3(2)=vect1(3)*vect2(1)-vect1(1)*vect2(3)
   vect3(3)=vect1(1)*vect2(2)-vect1(2)*vect2(1)
   if (iout > 0) then
     write(iout, '(a)' )' Third atom defining local coordinates : '
     write(iout, '(a,i4,a,i4)' )'     ib = ',ib,'   irpt = ',irpt
   end if
 end if

 ! Analysis and output of force constants, ordered with respect to the distance from atom ia
 do ii=1,ifcout
   index=list(ii)
   call canct9(Ifc%acell,Ifc%gprim,ib,index,irpt,Ifc%natom,Ifc%nrpt,Ifc%rcan,posngb(:,ii),Ifc%rprim,Ifc%rpt)
   indngb(ii)=ib
   dist1=dist(ia,ib,irpt)
   if (iout > 0) then
     write(iout, '(a)' )
     write(iout, '(i4,a,i6,a,i8)' )ii,' interaction with atom',ib,' cell',irpt
     write(iout, '(a,3es16.6)' )' with coordinates ',posngb(1:3,ii)*(one+tol8)
     write(iout, '(a,es16.6)' )' and distance ',dist1
   end if

   if(ifcana==1.and.ii/=1)then
     vect1(1)=(posngb(1,ii)-ra(1))/dist1
     vect1(2)=(posngb(2,ii)-ra(2))/dist1
     vect1(3)=(posngb(3,ii)-ra(3))/dist1
   end if

#ifdef MR_DEV
   if(Ifc%dipdip==0.or.Ifc%dipquad==1.or.Ifc%quadquad==1)then
#else
   if(Ifc%dipdip==0)then
#endif
     ! Get the "total" force constants (=real space FC)
     ! without taking into account the dipole-dipole interaction
     do mu=1,3
       do nu=1,3
         rsiaf(mu,nu,ii)=Ifc%atmfrc(mu,ia,nu,ib,irpt) * Ifc%wghatm(ia,ib,irpt)
       end do
     end do
     ! Output of the ifcs in cartesian coordinates
     if (iout > 0) then
       do nu=1,3
         write(iout, '(1x,3f9.5)' )(rsiaf(mu,nu,ii)+tol10,mu=1,3)
!       transfer short range and long range
         do mu=1,3
           Ifc%short_atmfrc(mu,ia,nu,ib,irpt) = rsiaf(mu,nu,ii) + tol10
         end do

       end do
     end if

     if(ifcana==1)then
       ! Further analysis
       trace1=rsiaf(1,1,ii)+rsiaf(2,2,ii)+rsiaf(3,3,ii)
       if (iout > 0) then
         write(iout, '(a,f9.5)' ) '  Trace         ',trace1+tol10
       end if
       if(ii/=1)then
         call axial9(rsiaf(:,:,ii),vect1,vect2,vect3)
       end if
       if (iout > 0) then
         write(iout, '(a)' )' Transformation to local coordinates '
         write(iout, '(a,3f16.6)' ) ' First  local vector :',vect1
         write(iout, '(a,3f16.6)' ) ' Second local vector :',vect2
         write(iout, '(a,3f16.6)' ) ' Third  local vector :',vect3
       end if
       call ifclo9(rsiaf(:,:,ii),ifcloc,vect1,vect2,vect3)
       if (iout > 0) then
         do nu=1,3
           write(iout, '(1x,3f9.5)' )(ifcloc(mu,nu)+tol10,mu=1,3)
         end do
       end if

       vect(:,1,ii) = vect1
       vect(:,2,ii) = vect2
       vect(:,3,ii) = vect3

     end if ! Further analysis finished

   else if(Ifc%dipdip==1)then

!DEBUG
!    write(iout,'(a)')
!    write(iout,'(a)')' Enter dipdip section, for debugging'
!    write(iout,'(a)')
!ENDDEBUG
     ! Get the Coulomb part
     do jj=1,3
       rdiff(jj)=ra(jj)-posngb(jj,ii)
     end do
     rsq=0.0_dp
     xx(1:3)=0.0_dp
     do jj=1,3
       do kk=1,3
         ewiaf0(jj,kk)=0.0_dp
         rsq=rsq+rdiff(jj)*invdlt(kk,jj)*rdiff(kk)
         xx(kk)=xx(kk)+invdlt(kk,jj)*rdiff(jj)
       end do
     end do
     yy=sqrt(rsq)
     !  Avoid zero denominators in term:
     if (sqrt(rsq)>=tol12) then
       do mu=1,3
         do nu=1,3
           ewiaf0(mu,nu)=(-3*xx(nu)*xx(mu)+invdlt(nu,mu)*yy**2)/yy**5/sqrt(detdlt)
         end do
       end do
     else
       if (ia/=ib)then
         write(message, '(a,a,a,a,a,i5,a,i5,a)' )&
&         'The distance between two atoms vanishes.',ch10,&
&         'This is not allowed.',ch10,&
&         'Action: check the input for the atoms number',ia,' and',ib,'.'
         MSG_ERROR(message)
       end if
     end if

     ! Take into account the effective charge tensor
     do mu=1,3
       do nu=1,3
         ew1=zero
         if(ii==1)then
           ew1=-Ifc%dyewq0(mu,nu,ia)
         end if
         do jj=1,3
           do kk=1,3
             ew1=ew1+zeff(jj,mu,ia)*(zeff(kk,nu,ib)* ewiaf0(jj,kk))
           end do
         end do
         ewiaf1(mu,nu)=ew1
       end do
     end do
     ! Get the short-range force constants and the
     ! "total" force constants (=real space FC)
     do mu=1,3
       do nu=1,3
         sriaf(mu,nu,ii)=Ifc%atmfrc(mu,ia,nu,ib,irpt)* Ifc%wghatm(ia,ib,irpt)
         rsiaf(mu,nu,ii)=ewiaf1(mu,nu)+sriaf(mu,nu,ii)
       end do
     end do

     ! Output of the results
     if (iout > 0) then
       do nu=1,3
         write(iout, '(1x,3(3f9.5,1x))' )&
&         (rsiaf(mu,nu,ii) +tol10,mu=1,3),&
&         (ewiaf1(mu,nu)+tol10,mu=1,3),&
&         (sriaf(mu,nu,ii) +tol10,mu=1,3)

!       transfer short range and long range
         do mu=1,3
           Ifc%short_atmfrc(mu,ia,nu,ib,irpt) = sriaf(mu,nu,ii) + tol10
           Ifc%ewald_atmfrc(mu,ia,nu,ib,irpt) = ewiaf1(mu,nu) + tol10
         end do
       end do
     end if

     if(ifcana==1)then
       ! Further analysis
       if (iout > 0) then
         write(iout, '(a)' )' Traces (and ratios) :'
       end if
       trace1=rsiaf(1,1,ii)+rsiaf(2,2,ii)+rsiaf(3,3,ii)
       trace2=ewiaf1(1,1)+ewiaf1(2,2)+ewiaf1(3,3)
       trace3=sriaf(1,1,ii)+sriaf(2,2,ii)+sriaf(3,3,ii)
       if (iout > 0) then
         write(iout,'(3(f9.5,17x))')trace1+tol10,trace2+tol10,trace3+tol10
         write(iout,'(3(f9.5,17x))')1.0,trace2/trace1+tol10,trace3/trace1+tol10 !
       end if

       if(ii/=1)then
         call axial9(rsiaf(:,:,ii),vect1,vect2,vect3)
       end if
       if (iout > 0) then
         write(iout, '(a)' )' Transformation to local coordinates '
         write(iout, '(a,3f16.6)' )' First  local vector :',vect1
         write(iout, '(a,3f16.6)' )' Second local vector :',vect2
         write(iout, '(a,3f16.6)' )' Third  local vector :',vect3
       end if
       call ifclo9(rsiaf(:,:,ii),rsloc,vect1,vect2,vect3)
       call ifclo9(ewiaf1,ewloc,vect1,vect2,vect3)
       call ifclo9(sriaf(:,:,ii),srloc,vect1,vect2,vect3)
       if (iout > 0) then
         do nu=1,3
           write(iout, '(1x,3(3f9.5,1x))' )&
&           (rsloc(mu,nu)+tol10,mu=1,3),&
&           (ewloc(mu,nu)+tol10,mu=1,3),&
&           (srloc(mu,nu)+tol10,mu=1,3)
         end do
         if(ii/=1)then
           write(iout, '(a)' )' Ratio with respect to the longitudinal ifc'
         else
           write(iout, '(a)' )' Ratio with respect to the (1,1) element'
         end if
         do nu=1,3
           write(iout, '(1x,3(3f9.5,1x))' )&
&           (rsloc(mu,nu)/rsloc(1,1)+tol10,mu=1,3),&
&           (ewloc(mu,nu)/rsloc(1,1)+tol10,mu=1,3),&
&           (srloc(mu,nu)/rsloc(1,1)+tol10,mu=1,3)
         end do
       end if

       vect(:,1,ii) = vect1
       vect(:,2,ii) = vect2
       vect(:,3,ii) = vect3

     end if ! Further analysis finished
   end if ! End the condition on dipdip
 end do ! End loop over all atoms in BigBox:

end subroutine ifc_getiaf
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/omega_decomp
!!
!! NAME
!!  omega_decomp
!!
!! FUNCTION
!! Compute and return the eigenvalues (frequencies) of the short-range and
!! long-range part of the dynamical matrix  See Europhys. Lett. 33 p.713 (1996) for details.
!! (included by U. Aschauer and EB)
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine omega_decomp(amu,natom,ntypat,typat,dynmatfl,dynmatsr,dynmatlr,iqpt,nqpt,eigenvec)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 integer,intent(in) :: iqpt,nqpt
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: amu(ntypat)
 real(dp),intent(inout) :: dynmatfl(2,3,natom,3,natom,nqpt)
 real(dp),intent(inout) :: dynmatsr(2,3,natom,3,natom,nqpt)
 real(dp),intent(inout) :: dynmatlr(2,3,natom,3,natom,nqpt)
 real(dp),intent(in)    :: eigenvec(2*3*natom*3*natom)

!Local variables -------------------------
!scalars
 integer :: i1,i2,idir1,idir2,imode,ipert1,ipert2,index1,index2
 real(dp),parameter :: break_symm=1.0d-12
 real(dp) :: fac
!arrays
 real(dp) :: nearidentity(3,3)
 real(dp) :: omegafl, omegasr, omegalr
 real(dp) :: sumfl,sumlr,sumsr,asr
! *********************************************************************

!write(ab_out,*)''
!write(std_out,*) 'SR/LR decomposition: enter for wavevector number :',iqpt

!apply asr (note the lr part is already asred by construction in mkifc9)
 do ipert1=1,natom
   do idir1=1,3
     do idir2=1,3
       asr=0.0d0
       do ipert2=1,natom
         asr=asr+dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)
       end do
       dynmatfl(1,idir1,ipert1,idir2,ipert1,iqpt)=dynmatfl(1,idir1,ipert1,idir2,ipert1,iqpt)-asr
       dynmatsr(1,idir1,ipert1,idir2,ipert1,iqpt)=dynmatsr(1,idir1,ipert1,idir2,ipert1,iqpt)-asr
     end do
   end do
 end do


!This slight breaking of the symmetry allows the
!results to be more portable between machines
 nearidentity(:,:)=1.0
 nearidentity(1,1)=1.0+break_symm
 nearidentity(3,3)=1.0-break_symm


!Include Mass
 do ipert1=1,natom
   do ipert2=1,natom

     fac=1.0d0/sqrt(amu(typat(ipert1))*amu(typat(ipert2)))/amu_emass

     do idir1=1,3
       do idir2=1,3

         dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&         dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)*&
&         fac*nearidentity(idir1,idir2)

         dynmatsr(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&         dynmatsr(1,idir1,ipert1,idir2,ipert2,iqpt)*&
&         fac*nearidentity(idir1,idir2)

         dynmatlr(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&         dynmatlr(1,idir1,ipert1,idir2,ipert2,iqpt)*&
&         fac*nearidentity(idir1,idir2)

         ! This is to break slightly the translation invariance, and make
         ! the automatic tests more portable
         if(ipert1==ipert2 .and. idir1==idir2)then
           dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&           dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)+&
&           break_symm*natom/amu_emass/idir1*0.01d0

           dynmatsr(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&           dynmatsr(1,idir1,ipert1,idir2,ipert2,iqpt)+&
&           break_symm*natom/amu_emass/idir1*0.01d0

           dynmatlr(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&           dynmatlr(1,idir1,ipert1,idir2,ipert2,iqpt)+&
&           break_symm*natom/amu_emass/idir1*0.01d0
         end if

       end do
     end do
   end do
 end do

!Calculation of <eigvec|Dyn_tot,Dyn_SR,Dyn_LR|eigenvec>=omega**2

!write(ab_out,*)''
!write(ab_out,*)'==============================================================================='
 write(ab_out,*)''
 write(ab_out,*) 'Long-Range/Short-Range decomposed phonon freq. (cm-1)**2'
 write(ab_out,*) 'at wavevector number:',iqpt
 write(ab_out,*)''
 write(ab_out,'(a13,1x,a16,2x,a16,2x,a16)') ' Mode number.','tot**2','SR**2','LR**2'
 write(std_out,'(a13,1x,a16,2x,a16,2x,a16)') ' Mode number.','tot**2','SR**2','LR**2'
!write(ab_out,'(a12,2x,a10,2x,a10,2x,a10,2x,a16,2x,a16,2x,a16)') 'Mode number.','tot','SR','LR','tot**2','SR**2','LR**2'
!write(std_out,'(a12,2x,a10,2x,a10,2x,a10,2x,a16,2x,a16,2x,a16)') 'Mode number.','tot','SR','LR','tot**2','SR**2','LR**2'

 do imode=1,3*natom
   sumfl=zero; sumlr=zero; sumsr=zero

   do ipert1=1,natom
     do ipert2=1,natom
       do i1=1,3
         do i2=1,3

           index1=i1+(ipert1-1)*3+3*natom*(imode-1)
           index2=i2+(ipert2-1)*3+3*natom*(imode-1)
           ! MG FIXME: I don't think these expressions are correct when q != 0
           ! We should also include the imaginary part

           sumfl = sumfl + eigenvec(2*index1-1) * dynmatfl(1,i1,ipert1,i2,ipert2,iqpt) * eigenvec(2*index2-1)
           sumlr = sumlr + eigenvec(2*index1-1) * dynmatlr(1,i1,ipert1,i2,ipert2,iqpt) * eigenvec(2*index2-1)
           sumsr = sumsr + eigenvec(2*index1-1) * dynmatsr(1,i1,ipert1,i2,ipert2,iqpt) * eigenvec(2*index2-1)
         end do
       end do
     end do
   end do

   sumfl = sumfl * Ha_cmm1 * Ha_cmm1
   sumsr = sumsr * Ha_cmm1 * Ha_cmm1
   sumlr = sumlr * Ha_cmm1 * Ha_cmm1

!  Compute omega=sqrt(omega**2)
   if(sumfl>=1.0d-16)then
     omegafl=sqrt(sumfl)
   else if(sumfl>=-1.0d-16)then
     omegafl=zero
   else
     omegafl=-sqrt(-sumfl)
   end if

   if(sumsr>=1.0d-16)then
     omegasr=sqrt(sumsr)
   else if(sumsr>=-1.0d-16)then
     omegasr=zero
   else
     omegasr=-sqrt(-sumsr)
   end if

   if(sumlr>=1.0d-16)then
     omegalr=sqrt(sumlr)
   else if(sumlr>=-1.0d-16)then
     omegalr=zero
   else
     omegalr=-sqrt(-sumlr)
   end if

!  Output
   write(ab_out,'(i4,10x,s,f16.4,2x,f16.4,2x,f16.4)') imode,sumfl,sumsr,sumlr  !vz_d
   write(std_out,'(i4,10x,s,f16.4,2x,f16.4,2x,f16.4)') imode,sumfl,sumsr,sumlr  !vz_d
!  write(ab_out,'(i4,8x,f10.4,2x,f10.4,2x,f10.4,2x,s,f16.6,2x,f16.6,2x,f16.6)') imode,omegafl,omegasr,omegalr,sumfl,sumsr,sumlr
!  write(std_out,'(i4,8x,f10.4,2x,f10.4,2x,f10.4,2x,s,f16.6,2x,f16.6,2x,f16.6)') imode,omegafl,omegasr,omegalr,sumfl,sumsr,sumlr
 end do

end subroutine omega_decomp
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_outphbtrap
!! NAME
!! ifc_outphbtrap
!!
!! FUNCTION
!!  Print out phonon frequencies on regular grid for BoltzTrap
!!  Flag in input file is outboltztrap=1
!!
!! INPUTS
!!  ifc<ifc_type>=Stores data related to interatomic force constants.
!!  Crystal<crystal_t>=Info on the crystal structure
!!  basename = file name for output to disk
!!  ngqpt(3)=Divisions of the q-mesh
!!  nqshft=Number of shifts
!!  qshft(3,nqshft)=Shifts of the q-mesh.
!!
!! OUTPUT
!!  only write to file. This routine should be called by a single processor.
!!
!! PARENTS
!!      anaddb,eph
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_outphbtrap(ifc, cryst, ngqpt, nqshft, qshft, basename)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nqshft
 character(len=*),intent(in) :: basename
 class(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: qshft(3,nqshft)

!Local variables -------------------------
!scalars
 integer,parameter :: qptopt1=1
 integer :: natom,imode,iq_ibz,nqbz,nqibz, nreals,unit_btrap,iatom,idir
 character(len=500) :: msg,format_nreals,format_line_btrap
 character(len=fnlen) :: outfile
!arrays
 integer :: qptrlatt(3,3)
 real(dp) :: d2cart(2,3,cryst%natom,3,cryst%natom),displ(2*3*cryst%natom*3*cryst%natom)
 real(dp) :: phfrq(3*cryst%natom),qphon(3)
 real(dp),allocatable :: qbz(:,:),qibz(:,:),wtq(:)

! *********************************************************************

 DBG_ENTER("COLL")

 natom = cryst%natom

 ! Setup IBZ, weights and BZ. Always use q --> -q symmetry for phonons even in systems wo inversion
 qptrlatt = 0; qptrlatt(1,1) = ngqpt(1); qptrlatt(2,2) = ngqpt(2); qptrlatt(3,3) = ngqpt(3)
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt1, nqshft, qshft, nqibz, qibz, wtq, nqbz, qbz)

 outfile = trim(basename) // '_BTRAP'
 write(msg, '(3a)')ch10,' Will write phonon FREQS in BoltzTrap format to file ',trim(outfile)
 call wrtout([std_out, ab_out], msg)

 if (open_file(outfile,msg,newunit=unit_btrap,status="replace") /= 0) then
   MSG_ERROR(msg)
 end if

 write (unit_btrap,'(a)') '#'
 write (unit_btrap,'(a)') '# ABINIT package : Boltztrap phonon file. With old BT versions remove this header before feeding to BT'
 write (unit_btrap,'(a)') '#    for compatibility with PHON output the freq are in Ry (before the square)'
 write (unit_btrap,'(a)') '#'
 write (unit_btrap,'(a)') '#    nq, nband  '
 write (unit_btrap,'(a)') '#  qx, qy, qz   '
 write (unit_btrap,'(a)') '#  qpt weight   '
 write (unit_btrap,'(a)') '#  freq_1^2, dynmat column for mode 1 '
 write (unit_btrap,'(a)') '#  etc for mode 2,3,4... qpt 2,3,4... '
 write (unit_btrap,'(2I6)') nqibz, 3*natom

! Loop over irreducible q-points
 do iq_ibz=1,nqibz
   qphon(:)=qibz(:,iq_ibz)

   call ifc_fourq(ifc, cryst, qphon, phfrq, displ, out_d2cart=d2cart)

   write (unit_btrap,'(3E20.10)') qphon
   write (unit_btrap,'(E20.10)') wtq(iq_ibz)
   nreals=1+2*3*natom
   call appdig(nreals,'(',format_nreals)
   format_line_btrap=trim(format_nreals)//'E20.10)'
   do iatom = 1, natom
     do idir = 1, 3
       imode = idir + 3*(iatom-1)
!      factor two for Ry output - this may change in definitive BT and abinit formats
       write (unit_btrap,trim(format_line_btrap))phfrq(imode)*two,d2cart(1:2,1:3,1:natom,idir,iatom)
     end do
   end do

 end do !irred q-points
 close (unit_btrap)

 ABI_DEALLOCATE(qibz)
 ABI_DEALLOCATE(qbz)
 ABI_DEALLOCATE(wtq)

 DBG_EXIT("COLL")

end subroutine ifc_outphbtrap
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_printbxsf
!! NAME
!! ifc_printbxsf
!!
!! FUNCTION
!!  Output phonon isosurface in Xcrysden format.
!!
!! INPUTS
!!  ifc<ifc_type>=Stores data related to interatomic force constants.
!!  crystal<crystal_t>=Info on the crystal structure
!!  ngqpt(3)=Divisions of the q-mesh
!!  nqshft=Number of shifts
!!  qshft(3,nqshft)=Shifts of the q-mesh.
!!  path=File name for output to disk
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Only write to file
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_printbxsf(ifc, cryst, ngqpt, nqshft, qshft, path, comm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nqshft,comm
 character(len=*),intent(in) :: path
 class(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: qshft(3,nqshft)

!Local variables -------------------------
!scalars
 integer,parameter :: nsppol1=1,master=0,qptopt1=1
 integer :: my_rank,nprocs,iq_ibz,nqibz,nqbz,ierr
 character(len=500) :: msg
!arrays
 integer :: qptrlatt(3,3),dummy_symafm(cryst%nsym)
 real(dp) :: displ_cart(2,3*cryst%natom,3*cryst%natom)
 real(dp),allocatable :: qibz(:,:),wtq(:),qbz(:,:),freqs_qibz(:,:)

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 ! Setup IBZ, weights and BZ. Always use q --> -q symmetry for phonons even in systems wo inversion
 qptrlatt = 0; qptrlatt(1,1) = ngqpt(1); qptrlatt(2,2) = ngqpt(2); qptrlatt(3,3) = ngqpt(3)
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt1, nqshft, qshft, nqibz, qibz, wtq, nqbz, qbz)
 ABI_FREE(qbz)
 ABI_FREE(wtq)

 ! Compute phonon frequencies in the irreducible wedge.
 ABI_CALLOC(freqs_qibz, (3*cryst%natom, nqibz))

 do iq_ibz=1,nqibz
   if (mod(iq_ibz, nprocs) /= my_rank) cycle ! mpi parallelism.
   call ifc_fourq(ifc, cryst, qibz(:,iq_ibz), freqs_qibz(:,iq_ibz), displ_cart)
 end do
 call xmpi_sum(freqs_qibz, comm, ierr)

 ! Output phonon isosurface.
 if (my_rank == master) then
   dummy_symafm = 1
   call printbxsf(freqs_qibz, zero, zero, cryst%gprimd, qptrlatt, 3*cryst%natom,&
     nqibz, qibz, cryst%nsym, .False., cryst%symrec, dummy_symafm, .True., nsppol1, qshft, nqshft, path, ierr)
   if (ierr /=0) then
     msg = "Cannot produce BXSF file with phonon isosurface, see log file for more info"
     MSG_WARNING(msg)
     call wrtout(ab_out, msg, 'COLL')
   end if
 end if

 ABI_FREE(freqs_qibz)
 ABI_FREE(qibz)

end subroutine ifc_printbxsf
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_calcnwrite_nana_terms
!! NAME
!!  ifc_calcnwrite_nana_terms
!!
!! FUNCTION
!!  Compute frequencies and phonon displacement for q-->0 in the presence of non-analytical behaviour.
!!
!! INPUTS
!!  nph2l=Number of qpoints.
!!  qph2l(3,nph2l)=List of phonon wavevector directions along which the non-analytical correction
!!    to the Gamma-point phonon frequencies will be calculated
!!    The direction is in CARTESIAN COORDINATES
!!  qnrml2(nph2l)=Normalization factor.
!!
!! OUTPUT
!!  (Optional)
!!  phfrq2l(3*crystal%natom,nph2l)=List of phonon frequencies
!!  polarity2l(3,3*crystal%natom,nph2l)=List of mode-polarities
!!     (see Eq.(41) of Veithen et al, PRB71, 125107 (2005) [[cite:Veithen2005]])
!!
!! NOTES:
!!  This routine should be called by master node and when ifcflag == 1.
!!
!! PARENTS
!!      m_gruneisen,m_phonons
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_calcnwrite_nana_terms(ifc, crystal, nph2l, qph2l, &
&  qnrml2, ncid, phfrq2l, polarity2l) ! optional arguments

!Arguments ------------------------------------
 integer,intent(in) :: nph2l
 integer,optional,intent(in) :: ncid
 class(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: crystal
!arrays
 real(dp),intent(in) :: qph2l(3, nph2l)
 real(dp),optional,intent(in) :: qnrml2(nph2l)
 real(dp),optional,intent(out) :: phfrq2l(3*crystal%natom,nph2l), polarity2l(3,3*crystal%natom,nph2l)

!Local variables-------------------------------
!scalars
 integer :: iatom,idir,imode,iphl2
#ifdef HAVE_NETCDF
 integer :: ncerr
#endif
!arrays
 real(dp) :: qphnrm(3),qphon(3,3)
 real(dp),allocatable :: displ_cart(:,:,:),phfrq(:),d2cart(:,:,:),eigvec(:,:,:),eigval(:)

! ************************************************************************

 if (nph2l == 0) return

 !Now treat the second list of vectors (only at the Gamma point, but can include non-analyticities)
 ABI_MALLOC(phfrq, (3*crystal%natom))
 ABI_MALLOC(displ_cart, (2, 3*crystal%natom, 3*crystal%natom))
 ABI_MALLOC(d2cart, (2, 3*ifc%mpert, 3*ifc%mpert))
 ABI_MALLOC(eigvec, (2, 3*crystal%natom, 3*crystal%natom))
 ABI_MALLOC(eigval, (3*crystal%natom))

 ! Before examining every direction or the dielectric tensor, generates the dynamical matrix at gamma
 qphon(:,1)=zero; qphnrm(1)=zero

 ! Generation of the dynamical matrix in cartesian coordinates
 ! Get d2cart using the interatomic forces and the
 ! long-range coulomb interaction through Ewald summation
 call gtdyn9(ifc%acell,ifc%atmfrc,ifc%dielt,ifc%dipdip, &
   ifc%dyewq0,d2cart,crystal%gmet,ifc%gprim,ifc%mpert,crystal%natom, &
   ifc%nrpt,qphnrm(1),qphon,crystal%rmet,ifc%rprim,ifc%rpt, &
   ifc%trans,crystal%ucvol,ifc%wghatm,crystal%xred,ifc%zeff,ifc%qdrp_cart,ifc%ewald_option,xmpi_comm_self)

#ifdef HAVE_NETCDF
 if (present(ncid)) then
   iphl2 = 0
   call nctk_defwrite_nonana_terms(ncid, iphl2, nph2l, qph2l, crystal%natom, phfrq, displ_cart, "define")
   ! Add epsinf, Born effective charges and some useful metadata.
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t('emacro_cart', "dp", 'number_of_cartesian_directions, number_of_cartesian_directions'), &
     nctkarr_t('becs_cart', "dp", "number_of_cartesian_directions, number_of_cartesian_directions, number_of_atoms")], &
     defmode=.True.)
   NCF_CHECK(ncerr)
   ! TODO chneut is missing
   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
       "asr", "dipdip", "symdynmat"])
   NCF_CHECK(ncerr)
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'emacro_cart'), ifc%dielt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'becs_cart'), ifc%zeff))
   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
     "asr", "dipdip", "symdynmat"], &
     [ifc%asr, ifc%dipdip, ifc%symdynmat])
   NCF_CHECK(ncerr)
 end if
#endif

 ! Examine every wavevector of this list
 do iphl2=1,nph2l
   ! Initialisation of the phonon wavevector
   qphon(:,1) = qph2l(:,iphl2)
   qphnrm(1)=zero
   if(present(qnrml2)) qphnrm(1) = qnrml2(iphl2)

   ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
   call dfpt_phfrq(ifc%amu,displ_cart,d2cart,eigval,eigvec,crystal%indsym, &
      ifc%mpert,crystal%nsym,crystal%natom,crystal%nsym,crystal%ntypat,phfrq,qphnrm(1),qphon, &
      crystal%rprimd,ifc%symdynmat,crystal%symrel,crystal%symafm,crystal%typat,crystal%ucvol)

   ! Write the phonon frequencies
   !call dfpt_prtph(displ_cart,inp%eivec,inp%enunit,ab_out,natom,phfrq,qphnrm(1),qphon)

   if(present(phfrq2l)) phfrq2l(:,iphl2)=phfrq(:)

   if(present(polarity2l))then
     polarity2l(:,:,iphl2)=zero
     do imode=1,3*crystal%natom
       do iatom=1,crystal%natom
         do idir=1,3
           polarity2l(:,imode,iphl2)=polarity2l(:,imode,iphl2)+ifc%zeff(:,idir,iatom)*displ_cart(1,idir+(iatom-1)*3,imode)
         enddo
       enddo
     enddo
   endif

#ifdef HAVE_NETCDF
   if(present(ncid))then
     ! Loop is not MPI-parallelized --> no need for MPI-IO API.
     call nctk_defwrite_nonana_terms(ncid, iphl2, nph2l, qph2l, crystal%natom, phfrq, displ_cart, "write")
   endif
#endif
 end do ! iphl2

 ABI_FREE(phfrq)
 ABI_FREE(displ_cart)
 ABI_FREE(d2cart)
 ABI_FREE(eigvec)
 ABI_FREE(eigval)

end subroutine ifc_calcnwrite_nana_terms
!!***

!----------------------------------------------------------------------

end module m_ifc
