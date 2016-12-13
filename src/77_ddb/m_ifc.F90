!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ifc
!! NAME
!!  m_ifc
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  used to handle interatomic force constant sets
!!
!! COPYRIGHT
!! Copyright (C) 2011-2016 ABINIT group (XG,MJV,EB,MG)
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
 use m_profiling_abi
 use m_xmpi
 use m_sort
 use m_bspline
 use m_skw
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_io_tools,    only : open_file
 use m_fstrings,    only : ktoa
 use m_time,        only : cwtime
 use m_numeric_tools, only : wrap2_zero_one, wrap2_pmhalf
 use m_copy,        only : alloc_copy
 use m_ewald,       only : ewald9
 use m_crystal,     only : crystal_t
 use m_geometry,    only : phdispl_cart2red
 use m_kpts,        only : kpts_ibz_from_kptrlatt
 use m_dynmat,      only : canct9, dist9 , ifclo9, axial9, q0dy3_apply, q0dy3_calc, asrif9, &
&                          make_bigbox, canat9, chkrp9, ftifc_q2r, wght9, symdm9, nanal9, gtdyn9, dymfz9, dfpt_phfrq
 use m_ddb,         only : ddb_type

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

   !integer :: ifcflag
   ! 0 if Fourier interpolation should not be used, 1 otherwise

   integer :: asr
     ! Option for the treatment of the Acoustic Sum Rule.

   integer :: brav
     ! Option for the sampling of the BZ (anaddb input variable)

   integer :: dipdip
     ! dipole dipole interaction flag.

   integer :: symdynmat
     ! If equal to 1, the dynamical matrix is symmetrized in dfpt_phfrq before the diagonalization.

   integer :: nqshft
     ! Number of shifts in the q-mesh (usually 1 since the mesh is gamma-centered!)

   integer :: nqbz
     ! Number of points in the full BZ

   integer :: nrpt
     ! Number of real space points used to integrate IFC (for interpolation of dynamical matrices)

   integer :: ngqpt(3)
    ! Number of division in the Q mesh.

   real(dp) :: rprim(3,3),gprim(3,3),acell(3)
     ! These values are used to call anaddb routines the don't use rprimd, gprimd

   ! These values will be stored in ddb but then we have to pass ddb to ifc_fourq
    real(dp) :: dielt(3,3)
     ! Dielectric tensor

   real(dp),allocatable :: amu(:)
     ! amu(ntypat)
     ! mass of the atoms (atomic mass unit)

   real(dp),allocatable :: atmfrc(:,:,:,:,:,:)
     ! atmfrc(2,3,natom,3,natom,nrpt)
     ! Inter atomic forces in real space

   integer,allocatable :: cell(:,:)
     ! cell(nrpt,3)
     ! Give the index of the the cell and irpt

   real(dp),allocatable :: ewald_atmfrc(:,:,:,:,:,:)
     ! Ewald_atmfrc(2,3,natom,3,natom,nrpt)
     ! Ewald Inter atomic forces in real space

   real(dp),allocatable :: short_atmfrc(:,:,:,:,:,:)
     ! short_atmfrc(2,3,natom,3,natom,nrpt)
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

   real(dp),allocatable :: qbz(:,:)
     ! qbz(3,nqpt))
     ! List of q-points in the full BZ

   real(dp),allocatable :: dynmat(:,:,:,:,:,:)
     ! dynmat(2,3,natom,3,natom,nqpt))
     ! dynamical matrices relative to the q points of the B.Z. sampling
     ! Note that the long-range dip-dip part has been removed if dipdip=1
     ! Moreover the array is multiplied by a phase shift in mkifc9.

   !real(dp),allocatable :: dynmat_lr(:,:,:,:,:,:)
    ! dynmat_lr(2,3,natom,3,natom,nqpt))
    ! Long-range part of dynmat in q-space
 end type ifc_type

 public :: ifc_init        ! Constructor
 public :: ifc_free        ! Release memory
 public :: ifc_fourq       ! Use Fourier interpolation to compute interpolated frequencies w(q) and eigenvectors e(q)
 public :: ifc_print       ! Print the ifc (output, netcdf and text file)
 public :: ifc_outphbtrap  ! Print out phonon frequencies on regular grid for BoltzTrap code.
 public :: ifc_printbxsf   ! Output phonon isosurface in Xcrysden format.
!!***

!----------------------------------------------------------------------

!!****t* m_ifc/phbspl_t
!! NAME
!! phbspl_t
!!
!! FUNCTION
!!
!! SOURCE

 type :: bcoeff_t
   real(dp),allocatable :: vals(:,:,:)
 end type bcoeff_t

 type,public :: phbspl_t

   integer :: natom3

   integer :: nqx,nqy,nqz
   ! Number of input data points

   integer :: qxord,qyord,qzord
   ! Order of the spline.

   !real(dp),allocatable :: xvec(:),yvec(:),zvec(:)
   real(dp),allocatable :: xknot(:),yknot(:),zknot(:)
   ! Array of length ndata+korder containing the knot

   type(bcoeff_t),allocatable :: coeff(:)
   ! coeff(natom3)

 end type phbspl_t

 public :: ifc_build_phbspl     ! Build B-spline object.
 public :: phbspl_evalq         ! Interpolate frequencies at arbitrary q-point.
 public :: phbspl_free          ! Free memory.
 public :: ifc_test_phinterp

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
!!      anaddb,eph,m_ifc
!!
!! CHILDREN
!!      appdig,ifc_fourq,smpbz,symkpt,wrtout
!!
!! SOURCE

subroutine ifc_free(ifc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(ifc_type),intent(inout) :: ifc

! ************************************************************************
 ! real
 if (allocated(ifc%amu)) then
   ABI_FREE(ifc%amu)
 end if

 if (allocated(ifc%atmfrc)) then
   ABI_FREE(ifc%atmfrc)
 end if

 if (allocated(ifc%cell)) then
   ABI_FREE(ifc%cell)
 end if

 if (allocated(ifc%ewald_atmfrc)) then
   ABI_FREE(ifc%ewald_atmfrc)
 end if

 if (allocated(ifc%short_atmfrc)) then
   ABI_FREE(ifc%short_atmfrc)
 end if

 if (allocated(ifc%qshft)) then
   ABI_FREE(ifc%qshft)
 end if

 if (allocated(ifc%rpt)) then
   ABI_FREE(ifc%rpt)
 end if

 if (allocated(ifc%wghatm)) then
   ABI_FREE(ifc%wghatm)
 end if

 if (allocated(ifc%rcan)) then
   ABI_FREE(ifc%rcan)
 end if

 if (allocated(ifc%trans)) then
   ABI_FREE(ifc%trans)
 end if

 if (allocated(ifc%dyewq0)) then
   ABI_FREE(ifc%dyewq0)
 end if

 if (allocated(ifc%qbz)) then
   ABI_FREE(ifc%qbz)
 end if

 if (allocated(ifc%zeff))then
   ABI_FREE(ifc%zeff)
 end if

 if (allocated(ifc%dynmat)) then
   ABI_FREE(ifc%dynmat)
 end if

 !if (allocated(ifc%dynmat_lr)) then
 !  ABI_FREE(ifc%dynmat_lr)
 !end if

end subroutine ifc_free
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_init
!! NAME
!!  ifc_init
!!
!! FUNCTION
!!  Initialize the dynamical matrix as well as the IFCs.
!!  taking into account the dipole-dipole interaction.
!!
!! INPUTS
!! crystal<type(crystal_t)> = Information on the crystalline structure.
!! ddb<type(ddb_type)> = Database with derivatives.
!! brav=bravais lattice (1=simple lattice,2=face centered lattice, 3=centered lattice,4=hexagonal lattice)
!! nqshft=Number of shifths in q-grid.
!! q1shft(3,nqshft)=Shifts for q-grid
!! asr= Option for the imposition of the ASR
!!   0 => no ASR,
!!   1 => modify "asymmetrically" the diagonal element
!!   2 => modify "symmetrically" the diagonal element
!! symdynmat=if 1, (re)symmetrize the dynamical matrix, except if Gamma wavevector with electric field added.
!! dipdip= if 0, no dipole-dipole interaction was subtracted in atmfrc
!!   if 1, atmfrc has been build without dipole-dipole part
!! rfmeth = 1 if non-stationary block
!!  2 if stationary block
!!  3 if third order derivatives
!! dielt(3,3)=dielectric tensor TODO: Should be computed from DDB
!! zeff(3,3,natom)=effective charge on each atom, versus electric field and atomic displacement
!! anaddb_dtset variables TODO: TO BE REMOVED
!!    anaddb_dtset%nsphere
!!    anaddb_dtset%prt_ifc
!!    anaddb_dtset%ifcana
!!    anaddb_dtset%rifcsph
!!    atifc(natom)
!!    anaddb_dtset%ifcout
!!    anaddb_dtset%prtsrlr
!!    anaddb_dtset%enunit
!! ddb = storage object for ddb information read in from DDB file
!! dielt(3,3)=dielectric tensor
!! ngqpt_in = input values of ngqpt
!! nsphere=number of atoms to be included in the cut-off sphere for interatomic
!!  force constant; if = 0 : maximum extent allowed by the grid.
!! rifcsph=radius for cutoff of IFC
!! [Ifc_coarse]=Optional.
!! [prtfreq]=True if phonon frequencies should be printed (default: false)
!!
!! OUTPUT
!! Ifc<type(ifc_type)>=Object containing the dynamical matrix and the IFCs.
!!   dyewq0(3,3,natom)=atomic self-interaction correction to the dynamical matrix. (only when dipdip=1)
!!   rcan(3,natom) = Atomic position in canonical coordinates
!!   trans(3,natom) = Atomic translations : xred = rcan + trans
!!   Ifc%nrpt = number of vectors of the lattice in real space
!!   Ifc%atmfrc(2,3,natom,3,natom,nrpt)= Interatomic Forces in real space
!!   Ifc%short_atmfrc(2,3,natom,3,natom,nrpt)= Interatomic Forces in real space
!!   Ifc%ewlad_atmfrc(2,3,natom,3,natom,nrpt)= Interatomic Forces in real space
!!   Ifc%rpt(3,nprt) = Canonical coordinates of the R points in the unit cell
!!   Ifc%wghatm(natom,natom,nrpt)= Weight associated to the couple of atoms and the R vector
!!
!! PARENTS
!!      anaddb,eph
!!
!! CHILDREN
!!      appdig,ifc_fourq,smpbz,symkpt,wrtout
!!
!! SOURCE

subroutine ifc_init(ifc,crystal,ddb,brav,asr,symdynmat,dipdip,&
  rfmeth,ngqpt_in,nqshft,q1shft,dielt,zeff,nsphere,rifcsph,&
  prtsrlr,enunit,& ! TO BE REMOVED
  prtfreq,Ifc_coarse) ! Optional


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_init'
 use interfaces_14_hidewrite
 use interfaces_56_recipspace
 use interfaces_72_response
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: asr,brav,dipdip,symdynmat,nqshft,rfmeth,nsphere
 real(dp),intent(in) :: rifcsph
 type(ifc_type),intent(inout) :: Ifc
 type(crystal_t),intent(in) :: Crystal
 type(ddb_type),intent(in) :: ddb
 type(ifc_type),optional,intent(in) :: Ifc_coarse

!arrays
 integer,intent(in) :: ngqpt_in(3)
 real(dp),intent(in) :: q1shft(3,nqshft)
 real(dp),intent(in) :: dielt(3,3),zeff(3,3,Crystal%natom)

!anaddb variables (TO BE REMOVED)
 integer,intent(in) :: prtsrlr,enunit
 logical,optional,intent(in) :: prtfreq
!end anaddb variables

!Local variables -------------------------
!scalars
 integer :: mpert,iout,iqpt,mqpt,nsym,ntypat,iq_bz,ii,natom
 integer :: nqbz,option,plus,sumg0,irpt,irpt_new
 real(dp),parameter :: qphnrm=one
 character(len=500) :: message
 logical :: do_prtfreq
 type(ifc_type) :: ifc_tmp
!arrays
 integer :: ngqpt(9),qptrlatt(3,3)
 integer,allocatable :: qmissing(:)
 real(dp) :: gprim(3,3),rprim(3,3),qpt(3),rprimd(3,3)
 real(dp):: rcan(3,Crystal%natom),trans(3,Crystal%natom),dyewq0(3,3,Crystal%natom)
 real(dp) :: displ_cart(2*3*Crystal%natom*3*Crystal%natom)
 real(dp) :: phfrq(3*Crystal%natom) !eigval(3,Crystal%natom),
 real(dp) :: eigvec(2,3,Crystal%natom,3,Crystal%natom)
 real(dp),allocatable :: dyew(:,:,:,:,:),spqpt(:,:)
 real(dp),allocatable :: dynmatfull(:,:,:,:,:,:),dynmat_sr(:,:,:,:,:,:),dynmat_lr(:,:,:,:,:,:) ! for OmegaSRLR
 real(dp),allocatable :: out_d2cart(:,:,:,:,:)

!******************************************************************

 ! This dimension should be encapsulated somewhere. We don't want to
 ! change the entire code if someone adds a new kind of perturbation.
 mpert = Crystal%natom + 6; iout = ab_out

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
 Ifc%dipdip = dipdip
 Ifc%symdynmat = symdynmat
 Ifc%nqshft = nqshft
 call alloc_copy(q1shft(:,1:Ifc%nqshft),Ifc%qshft)
 Ifc%ngqpt = ngqpt_in(1:3)
 Ifc%rprim = ddb%rprim
 Ifc%gprim = ddb%gprim
 Ifc%acell = ddb%acell

 ! Check if the rprim are coherent with the choice used in the interatomic forces generation
 call chkrp9(Ifc%brav,rprim)

 dyewq0 = zero
 if (Ifc%dipdip==1 .and. (Ifc%asr==1.or.Ifc%asr==2)) then
   ! Calculation of the non-analytical part for q=0
   sumg0=0
   qpt(:)=zero
   ABI_MALLOC(dyew,(2,3,natom,3,natom))
   call ewald9(ddb%acell,dielt,dyew,Crystal%gmet,gprim,natom,qpt,Crystal%rmet,rprim,sumg0,Crystal%ucvol,Crystal%xred,zeff)
   option=Ifc%asr
   call q0dy3_calc(natom,dyewq0,dyew,option)
   ABI_FREE(dyew)
 end if

 ! Sample the Brillouin zone
 option=1
 qptrlatt(:,:)=0
 qptrlatt(1,1)=ngqpt(1)
 qptrlatt(2,2)=ngqpt(2)
 qptrlatt(3,3)=ngqpt(3)
 mqpt=ngqpt(1)*ngqpt(2)*ngqpt(3)*nqshft
 if (Ifc%brav==2) mqpt=mqpt/2
 if (Ifc%brav==3) mqpt=mqpt/4

 ABI_MALLOC(spqpt,(3,mqpt))
 call smpbz(Ifc%brav,ab_out,qptrlatt,mqpt,nqbz,nqshft,option,q1shft,spqpt)

 ABI_MALLOC(Ifc%dynmat,(2,3,natom,3,natom,nqbz))

! Find symmetrical dynamical matrices
 if (.not.present(Ifc_coarse)) then

   ! Each q-point in the BZ mush be the symmetrical of one of the qpts in the ddb file.
   call symdm9(ddb%flg,ddb%nrm,ddb%qpt,ddb%typ,ddb%val,&
&    Ifc%dynmat,gprim,Crystal%indsym,mpert,natom,ddb%nblok,nqbz,nsym,rfmeth,rprim,spqpt,&
&    Crystal%symrec,Crystal%symrel)

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
&    Ifc%dynmat,gprim,Crystal%indsym,mpert,natom,ddb%nblok,nqbz,nsym,rfmeth,rprim,spqpt,&
&    Crystal%symrec,Crystal%symrel,qmissing=qmissing)

   ! Compute dynamical matrix with Fourier interpolation on the coarse q-mesh.
   write(message,"(a,i0,a)")"Will use Fourier interpolation to construct D(q) for ",size(qmissing)," q-points"
   call wrtout(std_out,message,"COLL")

   ABI_MALLOC(out_d2cart, (2,3,natom,3,natom))
   do ii=1,size(qmissing)
     iq_bz = qmissing(ii)
     qpt = spqpt(:,iq_bz)
     ! TODO: check dipdip option and phase, but I think this is correct!
     call ifc_fourq(Ifc_coarse,Crystal,qpt,phfrq,displ_cart,out_d2cart=out_d2cart)
     Ifc%dynmat(:,:,:,:,:,iq_bz) = out_d2cart

   end do

   ABI_FREE(qmissing)
   ABI_FREE(out_d2cart)
 end if

!OmegaSRLR: Store full dynamical matrix for decomposition into short- and long-range parts
 ABI_MALLOC(dynmatfull,(2,3,natom,3,natom,nqbz))
 dynmatfull=Ifc%dynmat

 if (Ifc%dipdip==1) then
   ! Take off the dipole-dipole part of the dynamical matrix
   write(message, '(3a)' )' mkifc9 : will extract the dipole-dipole part,',ch10,&
&                         ' using ewald9, q0dy3 and nanal9 for every wavevector.'
   call wrtout(std_out,message,"COLL")

   ABI_MALLOC(dyew,(2,3,natom,3,natom))

   do iqpt=1,nqbz
     qpt(:)=spqpt(:,iqpt)
     sumg0=0
     call ewald9(ddb%acell,dielt,dyew,Crystal%gmet,gprim,natom,qpt,Crystal%rmet,rprim,sumg0,Crystal%ucvol,Crystal%xred,zeff)
     call q0dy3_apply(natom,dyewq0,dyew)
     plus=0
     call nanal9(dyew,Ifc%dynmat,iqpt,natom,nqbz,plus)
   end do

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
 call dymfz9(Ifc%dynmat,natom,nqbz,gprim,option,spqpt,trans)

! Create the Big Box of R vectors in real space and compute the number of points (cells) in real space
 call make_bigbox(Ifc%brav,ifc_tmp%cell,ngqpt,nqshft,rprim,ifc_tmp%nrpt,ifc_tmp%rpt)

! Weights associated to these R points and to atomic pairs
! MG FIXME: Why ngqpt is intent(inout)?
 ABI_MALLOC(ifc_tmp%wghatm,(natom,natom,ifc_tmp%nrpt))
 call wght9(Ifc%brav,gprim,natom,ngqpt,nqbz,nqshft,ifc_tmp%nrpt,q1shft,rcan,ifc_tmp%rpt,ifc_tmp%wghatm)

! Fourier transformation of the dynamical matrices (q-->R)
 ABI_MALLOC(ifc_tmp%atmfrc,(2,3,natom,3,natom,ifc_tmp%nrpt))
 call ftifc_q2r(ifc_tmp%atmfrc,Ifc%dynmat,gprim,natom,nqbz,ifc_tmp%nrpt,ifc_tmp%rpt,spqpt)

! Eventually impose Acoustic Sum Rule to the interatomic forces
 if (Ifc%asr>0) then
   call asrif9(Ifc%asr,ifc_tmp%atmfrc,natom,ifc_tmp%nrpt,ifc_tmp%rpt,ifc_tmp%wghatm)
 end if



!*** The interatomic forces have been calculated ! ***
 write(message, '(2a)')ch10,' The interatomic forces have been obtained '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')


 ! Apply cutoff on ifc if needed
 if (nsphere/=0 .or. rifcsph>tol10) then
   call wrtout(std_out, 'ifc_init: apply cutoff on IFCs', "COLL")

   call corsifc9(ddb%acell,gprim,natom,ifc_tmp%nrpt,nsphere,rifcsph,&
&   rcan,rprim,ifc_tmp%rpt,ifc_tmp%wghatm)
 end if

 !ABI_FREE(dynmat)

!Eventually impose Acoustic Sum Rule to the interatomic forces
!(Note : here, after the analysis, in which the range
!of the short-range IFCs may have been changed
!That is why it is asked that asr be negative)
!FIXME: asr < 0 is not tested in abinit suite and does not appear to work
!(I get frequencies of 10^105 Hartree...) Modifying this 12/6/2011
!if(asr<0)then
!asr=-asr

 ! Be careful here: if I move this call to anaddb, several tests will fail
 ! due to the different number of points in the big box (see code below.)
 if (Ifc%asr > 0 .and. (nsphere/=0 .or. rifcsph>tol10)) then
   call asrif9(Ifc%asr,ifc_tmp%atmfrc,natom,ifc_tmp%nrpt,ifc_tmp%rpt,ifc_tmp%wghatm)
 end if

 ! Only conserve the necessary points in rpt: in the FT algorithm the order of the points is unimportant
 ! In the case of effective potential, we need to keep all the points
 Ifc%nrpt = 0
 do irpt=1,ifc_tmp%nrpt
   if (sum(ifc_tmp%wghatm(:,:,irpt)) /= 0) then
     Ifc%nrpt = Ifc%nrpt+1
   end if
 end do
 write(message,"(2(a,i0))")"ifc%nrpt: ",ifc%nrpt,", nqbz: ",nqbz
 call wrtout(std_out,message,"COLL")

 ABI_MALLOC(Ifc%atmfrc,(2,3,natom,3,natom,Ifc%nrpt))
 ABI_MALLOC(Ifc%rpt,(3,Ifc%nrpt))
 ABI_MALLOC(Ifc%cell,(3,Ifc%nrpt))
 ABI_MALLOC(Ifc%wghatm,(natom,natom,Ifc%nrpt))
 ABI_MALLOC(Ifc%short_atmfrc,(2,3,natom,3,natom,Ifc%nrpt))
 ABI_MALLOC(Ifc%ewald_atmfrc,(2,3,natom,3,natom,Ifc%nrpt))

 Ifc%short_atmfrc = zero
 Ifc%ewald_atmfrc = zero
 Ifc%short_atmfrc = zero
 Ifc%ewald_atmfrc = zero
 Ifc%cell = zero

 irpt_new = 1
 do irpt = 1, ifc_tmp%nrpt
   if (sum(ifc_tmp%wghatm(:,:,irpt)) /= 0) then
     Ifc%atmfrc(:,:,:,:,:,irpt_new) = ifc_tmp%atmfrc(:,:,:,:,:,irpt)
     Ifc%rpt(:,irpt_new) = ifc_tmp%rpt(:,irpt)
     Ifc%wghatm(:,:,irpt_new) = ifc_tmp%wghatm(:,:,irpt)
     Ifc%cell(:,irpt_new) = ifc_tmp%cell(:,irpt)
     irpt_new = irpt_new + 1
   end if
 end do

 ! Copy other useful arrays.
 Ifc%dielt = dielt
 Ifc%nqbz = nqbz

 call alloc_copy(rcan, Ifc%rcan)
 call alloc_copy(trans, Ifc%trans)
 call alloc_copy(dyewq0, Ifc%dyewq0)
 call alloc_copy(spqpt(:,1:nqbz), Ifc%qbz)
 call alloc_copy(zeff, Ifc%zeff)
 call alloc_copy(ddb%amu, Ifc%amu)

 call ifc_free(ifc_tmp)


 ! TODO (This is to be suppressed in a future version)
 do_prtfreq = .False.; if (present(prtfreq)) do_prtfreq = prtfreq
 if (do_prtfreq .or. prtsrlr == 1) then
   ! Check that the starting values are well reproduced.
   write(message, '(2a)' )' mkifc9 : now check that the starting values ',&
     ' are reproduced after the use of interatomic forces '
   call wrtout(std_out,message,"COLL")
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

 ABI_FREE(spqpt)

!OmegaSRLR: deallocate memory used by dynmat decomposition
 ABI_FREE(dynmatfull)
 ABI_FREE(dynmat_sr)
 ABI_FREE(dynmat_lr)

end subroutine ifc_init
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_fourq
!! NAME
!!  ifc_fourq
!!
!! FUNCTION
!!  Compute the phonon frequencies at the specified q-point by performing
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
!!
!! OUTPUT
!!  phfrq(3*natom) = Phonon frequencies in Hartree
!!  displ_cart(2,3,natom,3*natom) = Phonon displacement in Cartesian coordinates
!!  [out_d2cart(2,3,3*natom,3,3*natom)] = The (interpolated) dynamical matrix for this q-point
!!  [out_eigvec(2*3*natom*3*natom) = The (interpolated) eigenvectors of the dynamical matrix.
!!  [out_displ_red(2*3*natom*3*natom) = The (interpolated) displacement in reduced coordinates.
!!
!! PARENTS
!!      get_nv_fs_en,get_tau_k,harmonic_thermo,interpolate_gkk,m_ifc,m_phgamma
!!      m_phonons,mka2f,mka2f_tr,mka2f_tr_lova,mkph_linwid,read_gkk
!!
!! CHILDREN
!!      appdig,ifc_fourq,smpbz,symkpt,wrtout
!!
!! SOURCE


subroutine ifc_fourq(ifc, crystal, qpt, phfrq, displ_cart, &
                     nanaqdir, out_d2cart, out_eigvec, out_displ_red)   ! Optional [out]


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_fourq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),optional,intent(in) :: nanaqdir
 type(ifc_type),intent(in) :: Ifc
 type(crystal_t),intent(in) :: Crystal
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: displ_cart(2,3,Crystal%natom,3*Crystal%natom)
 real(dp),intent(out) :: phfrq(3*Crystal%natom)
 real(dp),optional,intent(out) :: out_d2cart(2,3,Crystal%natom,3,Crystal%natom)
 real(dp),optional,intent(out) :: out_eigvec(2,3,Crystal%natom,3*Crystal%natom)
 real(dp),optional,intent(out) :: out_displ_red(2,3,Crystal%natom,3*Crystal%natom)

!Local variables-------------------------------
!scalars
 integer :: natom
 real(dp) :: qphnrm
!arrays
 real(dp) :: my_qpt(3),eigvec(2,3,Crystal%natom,3*Crystal%natom),eigval(3*Crystal%natom)
 real(dp) :: d2cart(2,3,Ifc%mpert,3,Ifc%mpert)

! ************************************************************************

 natom = Crystal%natom

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
     MSG_ERROR("Wrong value for nanaqdir: "//trim(nanaqdir))
   end select
 end if

 ! The dynamical matrix d2cart is calculated here:
 call gtdyn9(Ifc%acell,Ifc%atmfrc,Ifc%dielt,Ifc%dipdip,Ifc%dyewq0,d2cart,Crystal%gmet,Ifc%gprim,Ifc%mpert,natom,&
&  Ifc%nrpt,qphnrm,my_qpt,Crystal%rmet,Ifc%rprim,Ifc%rpt,Ifc%trans,Crystal%ucvol,Ifc%wghatm,Crystal%xred,Ifc%zeff)

 ! Calculate the eigenvectors and eigenvalues of the dynamical matrix
 call dfpt_phfrq(Ifc%amu,displ_cart,d2cart,eigval,eigvec,Crystal%indsym,&
&  Ifc%mpert,Crystal%nsym,natom,Crystal%nsym,Crystal%ntypat,phfrq,qphnrm,my_qpt,&
&  Crystal%rprimd,Ifc%symdynmat,Crystal%symrel,Crystal%symafm,Crystal%typat,Crystal%ucvol)

 ! OmegaSRLR: Perform decomposition of dynamical matrix
 !if (srlr==1)then
 !  call omega_decomp(amu,natom,ntypat,typat,dynmatfull,dynmatsr,dynmatlr,iqpt,nqpt,eigvec)
 !end if

 ! Return the interpolated dynamical matrix and the eigenvector for this q-point
 if (present(out_d2cart)) out_d2cart = d2cart(:,:,:natom,:,:natom)
 if (present(out_eigvec)) out_eigvec = eigvec

 ! Return phonon displacement in reduced coordinates.
 if (present(out_displ_red)) call phdispl_cart2red(natom, crystal%gprimd, displ_cart, out_displ_red)

 ! Option to get vectors in reduced coordinates?
 !call phdispl_cart2red(natom, crystal%gprimd, out_eigvec, out_eigvec_red)

 !if (present(dwdq) call ifc_get_dwdq(ifc, crystal, my_qpt, phfrq, eigvec, dwdq)

end subroutine ifc_fourq
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/corsifc9
!!
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
!!
!! NOTES
!! This routine should be executed by one processor only
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!      appdig,ifc_fourq,smpbz,symkpt,wrtout
!!
!! SOURCE

subroutine corsifc9(acell,gprim,natom,nrpt,nsphere,rifcsph,rcan,rprim,rpt,wghatm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'corsifc9'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nrpt,nsphere
 real(dp),intent(in) :: rifcsph
!arrays
 real(dp),intent(in) :: acell(3)
 real(dp),intent(in) :: gprim(3,3),rcan(3,natom)
 real(dp),intent(in) :: rprim(3,3),rpt(3,nrpt)
 real(dp),intent(inout) :: wghatm(natom,natom,nrpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,ii,index,irpt
!arrays
 integer,allocatable :: list(:)
 real(dp),allocatable :: dist(:,:,:),wkdist(:)

! *********************************************************************

 ! Compute the distances between atoms
 ABI_MALLOC(dist,(natom,natom,nrpt))
 call dist9(acell,dist,gprim,natom,nrpt,rcan,rprim,rpt)
 ! Now dist(ia,ib,irpt) contains the distance from atom ia to atom ib in unit cell irpt.

 ABI_MALLOC(list,(natom*nrpt))
 ABI_MALLOC(wkdist,(natom*nrpt))

 ! loop on all generic atoms.
 do ia=1,natom

   wkdist(:)=reshape(dist(ia,:,:),(/natom*nrpt/))
   do ii=1,natom*nrpt
     list(ii)=ii
   end do
   ! This sorting algorithm is slow ...
   call sort_dp(natom*nrpt,wkdist,list,tol14)

   ! zero the outside IFCs : act on wghatm
   if(nsphere/=0.and.nsphere<natom*nrpt)then
     do ii=nsphere+1,natom*nrpt
       index=list(ii)
       irpt=(index-1)/natom+1
       ib=index-natom*(irpt-1)
       wghatm(ia,ib,irpt)=0._dp
     end do
   end if

   if(rifcsph>tol10)then
     do ii=nsphere+1,natom*nrpt
       index=list(ii)
       ! preserve weights for atoms inside sphere of radius rifcsph
       if (wkdist(ii) < rifcsph) cycle
       irpt=(index-1)/natom+1
       ib=index-natom*(irpt-1)
       wghatm(ia,ib,irpt)=0._dp
     end do
   end if

 end do

 ABI_FREE(dist)
 ABI_FREE(list)
 ABI_FREE(wkdist)

end subroutine corsifc9
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_print
!!
!! NAME
!! ifc_print
!!
!! FUNCTION
!! Adds the real-space interatomic force constants to the output file,
!! the NetCDF file and, if prt_ifc==1, to the ifcinfo.out file.
!!
!! INPUTS
!! Ifc<type(ifc_type)>=Object containing the dynamical matrix and the IFCs.
!! dielt(3,3)=dielectric tensor
!! zeff(3,3,natom)=effective charge on each atom, versus electric field and atomic displacement
!! ifcana= 0 => no analysis of ifc ; 1 => full analysis
!! atifc(natom) =  atifc(ia) equals 1 if the analysis of ifc
!!  has to be done for atom ia; otherwise 0.
!! ifcout= Number of interatomic force constants written in the output file
!! prt_ifc = flag to print out ifc information for dynamical matrix (AI2PS)
!! ncid=the unit of the open NetCDF file.
!!
!! OUTPUT
!!   written in the output file and in the NetCDF file
!!
!! NOTES
!! This routine should be executed by one processor only
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      appdig,ifc_fourq,matr3inv,mkrdim,smpbz,symkpt,wrtout
!!
!! SOURCE

subroutine ifc_print(Ifc,dielt,zeff,ifcana,atifc,ifcout,prt_ifc,ncid)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_print'
 use interfaces_32_util
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ifcout,ifcana
 integer,intent(in) :: prt_ifc
 integer,optional,intent(in) :: ncid
 type(ifc_type),intent(inout) :: Ifc
!arrays
 integer,intent(in) :: atifc(Ifc%natom)
 real(dp),intent(in) :: dielt(3,3)
 real(dp),intent(in) :: zeff(3,3,Ifc%natom)
!Local variables -------------------------
!scalars
 integer :: ia,ii,ncerr,iatifc,ifcout1,mu,nu,iout
! unit number to print out ifc information for dynamical matrix (AI2PS)
 integer :: unit_ifc
 real(dp) :: detdlt
 character(len=500) :: message
!arrays
 integer,allocatable :: list(:)
 real(dp) :: invdlt(3,3),ra(3),xred(3)
 real(dp),allocatable :: dist(:,:,:),wkdist(:),rsiaf(:,:,:),sriaf(:,:,:),vect(:,:,:)
 real(dp) :: gprimd(3,3),rprimd(3,3)
 real(dp),allocatable :: indngb(:),posngb(:,:)

! *********************************************************************

 iout = ab_out

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

 write(std_out,'(a)' )' ifc_print : analysis of interatomic force constants '
 call mkrdim(Ifc%acell,Ifc%rprim,rprimd)
 call matr3inv(rprimd,gprimd)

 if (iout > 0) then
   write(iout, '(/,a,/)' )' Analysis of interatomic force constants '
   if(Ifc%dipdip==1)then
     write(iout, '(a)' )' Are given : column(1-3), the total force constant'
     write(iout, '(a)' )'       then  column(4-6), the Ewald part'
     write(iout, '(a)' )'       then  column(7-9), the short-range part'
     write(iout, '(a)' )' Column 1, 4 and 7 are related to the displacement'
     write(iout, '(a)' )'       of the generic atom along x,               '
     write(iout, '(a)' )' column 2, 5 and 8 are related to the displacement'
     write(iout, '(a)' )'       of the generic atom along y,               '
     write(iout, '(a)' )' column 3, 6 and 9 are related to the displacement'
     write(iout, '(a)')'       of the generic atom along z.               '
   else if(Ifc%dipdip==0)then
     write(iout, '(a)' )' column 1 is related to the displacement'
     write(iout, '(a)' )'        of the generic atom along x,    '
     write(iout, '(a)' )' column 2 is related to the displacement'
     write(iout, '(a)' )'        of the generic atom along y,    '
     write(iout, '(a)' )' column 3 is related to the displacement'
     write(iout, '(a)' )'        of the generic atom along z,    '
   end if
 end if

 if (ifcout>Ifc%natom*Ifc%nrpt) then
   ifcout1=Ifc%natom*Ifc%nrpt
   write(message, '(3a,i0,a)' )&
&   'The value of ifcout exceed the number of atoms in the big box.', ch10, &
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

#ifdef HAVE_NETCDF
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
#endif

 end if
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

     call ifc_getiaf(Ifc,ifcana,ifcout1,iout,zeff,ia,ra,list,dist,invdlt,&
&                    detdlt,rsiaf,sriaf,vect,indngb,posngb)


     if (prt_ifc == 1) then
       do ii=1,ifcout1
         write(unit_ifc,'(i6,i6)') ia,ii
         write(unit_ifc,'(3es28.16)') posngb(1:3,ii)
         do nu=1,3
           write(unit_ifc,   '(3f28.16)'  )(rsiaf(nu,mu,ii),mu=1,3)
         end do
       end do

#ifdef HAVE_NETCDF
       NCF_CHECK(nf90_put_var(ncid, vid("ifc_atoms_indices"), ia, start=(/iatifc/)))
       NCF_CHECK(nf90_put_var(ncid, vid("ifc_neighbours_indices"), indngb, start=(/1,iatifc/), count=(/ifcout1,1/)))
       NCF_CHECK(nf90_put_var(ncid, vid("ifc_distances"), wkdist(:ifcout1), start=(/1,iatifc/),count=(/ifcout1,1/)))
       ncerr = nf90_put_var(ncid, vid("ifc_matrix_cart_coord"), rsiaf, &
         start=(/1,1,1,iatifc/), count=(/3,3,ifcout1,1/))
       NCF_CHECK(ncerr)
       if (Ifc%dipdip==1) then
         ncerr = nf90_put_var(ncid, vid("ifc_matrix_cart_coord_short_range"), sriaf, &
           start=(/1,1,1,iatifc/), count=(/3,3,ifcout1,1/))
         NCF_CHECK(ncerr)
       end if
       if (ifcana==1) then
         ncerr = nf90_put_var(ncid, vid("ifc_local_vectors"), vect, &
           start=(/1,1,1,iatifc/), count=(/3,3,ifcout1,1/))
         NCF_CHECK(ncerr)
       end if
#endif
     end if
   end if ! End the condition on atifc

 end do ! End Big loop on atoms in the unit cell, and corresponding test .

 ABI_FREE(rsiaf)
 ABI_FREE(sriaf)
 ABI_FREE(vect)
 ABI_FREE(indngb)
 ABI_FREE(posngb)

 if (prt_ifc == 1) close(unit_ifc)

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

end subroutine ifc_print

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
!!      appdig,ifc_fourq,smpbz,symkpt,wrtout
!!
!! SOURCE

subroutine ifc_getiaf(Ifc,ifcana,ifcout,iout,zeff,ia,ra,list,&
& dist,invdlt,detdlt,rsiaf,sriaf,vect,indngb,posngb)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_getiaf'
!End of the abilint section

implicit none

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
 real(dp),intent(out) :: rsiaf(3,3,ifcout),sriaf(3,3,ifcout),vect(3,3,ifcout),indngb(ifcout),posngb(3,ifcout)

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

   if(Ifc%dipdip==0)then
     ! Get the "total" force constants (=real space FC)
     ! without taking into account the dipole-dipole interaction
     do mu=1,3
       do nu=1,3
         rsiaf(mu,nu,ii)=Ifc%atmfrc(1,mu,ia,nu,ib,irpt) * Ifc%wghatm(ia,ib,irpt)
       end do
     end do
     ! Output of the ifcs in cartesian coordinates
     if (iout > 0) then
       do nu=1,3
         write(iout, '(1x,3f9.5)' )(rsiaf(mu,nu,ii)+tol10,mu=1,3)
!       transfer short range and long range
         do mu=1,3
           Ifc%short_atmfrc(1,mu,ia,nu,ib,irpt) = rsiaf(mu,nu,ii) + tol10
           Ifc%ewald_atmfrc(1,mu,ia,nu,ib,irpt) = zero
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
         sriaf(mu,nu,ii)=Ifc%atmfrc(1,mu,ia,nu,ib,irpt)* Ifc%wghatm(ia,ib,irpt)
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
           Ifc%short_atmfrc(1,mu,ia,nu,ib,irpt) = sriaf(mu,nu,ii) + tol10
           Ifc%ewald_atmfrc(1,mu,ia,nu,ib,irpt) = ewiaf1(mu,nu) + tol10
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
!!      appdig,ifc_fourq,smpbz,symkpt,wrtout
!!
!! SOURCE

subroutine omega_decomp(amu,natom,ntypat,typat,dynmatfl,dynmatsr,dynmatlr,iqpt,nqpt,eigenvec)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'omega_decomp'
!End of the abilint section

 implicit none

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

   sumfl=0.0d0
   sumlr=0.0d0
   sumsr=0.0d0

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
!!      appdig,ifc_fourq,smpbz,symkpt,wrtout
!!
!! SOURCE

subroutine ifc_outphbtrap(ifc, cryst, ngqpt, nqshft, qshft, basename)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_outphbtrap'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nqshft
 character(len=*),intent(in) :: basename
 type(ifc_type),intent(in) :: ifc
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
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 if (open_file(outfile,msg,newunit=unit_btrap,status="replace") /= 0) then
   MSG_ERROR(msg)
 end if

 write (unit_btrap,'(a)') '#'
 write (unit_btrap,'(a)') '# ABINIT package : Boltztrap phonon file. Remove this header before feeding to BT'
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
 close (unit=unit_btrap)

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_printbxsf'
 use interfaces_14_hidewrite
 use interfaces_61_occeig
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nqshft,comm
 character(len=*),intent(in) :: path
 type(ifc_type),intent(in) :: ifc
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
 real(dp) :: phfrq(3*cryst%natom),displ_cart(2,3*cryst%natom,3*cryst%natom)
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
   if (mod(iq_ibz, nprocs) /= my_rank) cycle
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

!!****f* m_ifc/ifc_build_phbspl
!! NAME
!! ifc_build_phbspl
!!
!! FUNCTION
!! build the phbspl_t object to interpolate the phonon band structure.
!!
!! INPUTS
!!  ifc<type(ifc_type)>=Object containing the dynamical matrix and the IFCs.
!!  crystal<type(crystal_t)> = Information on the crystalline structure.
!!  ngqpt(3)=Divisions of the q-mesh used to produce the B-spline.
!!  nshiftq=Number of shifts in Q-mesh
!!  shiftq(3,nshiftq)=Shifts of the q-mesh.
!!  ords(3)=order of the spline for the three directions. ord(1) must be in [0, nqx] where
!!    nqx is the number of points along the x-axis.
!!  comm=MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

type(phbspl_t) function ifc_build_phbspl(ifc, cryst, ngqpt, nshiftq, shiftq, ords, comm) result(new)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_build_phbspl'
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,nshiftq
 type(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: ngqpt(3),ords(3)
 real(dp),intent(in) :: shiftq(3,nshiftq)

!local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1=1,timrev1=1,qptopt1=1
 integer :: qxord,qyord,qzord,nxknot,nyknot,nzknot,nqibz,nqbz,nqfull,iqf,ierr
 integer :: nu,iq_ibz,ix,iy,iz,nqx,nqy,nqz,my_rank,nprocs,cnt
 real(dp) :: dksqmax
 character(len=500) :: msg
!arrays
 integer :: qptrlatt(3,3)
 integer,allocatable :: bz2ibz(:,:)
 real(dp) :: phfrq(3*cryst%natom),qpt(3),displ_cart(2,3,cryst%natom,3*cryst%natom)
 real(dp),allocatable :: xvec(:),yvec(:),zvec(:),xyzdata(:,:,:)
 real(dp),allocatable :: ibz_freqs(:,:),ibzdata_qnu(:,:)
 real(dp),allocatable :: wtq(:),qbz(:,:),qfull(:,:),qibz(:,:)

! *********************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Check input parameters
 ierr = 0
 if (nshiftq /= 1) then
   MSG_WARNING('Multiple shifts not allowed')
   ierr = ierr + 1
 end if
 if (ierr /= 0) then
   MSG_WARNING("bspline interpolation cannot be performed. See warnings above. Returning")
   return
 end if

 ! Setup IBZ, weights and BZ. Always use q --> -q symmetry for phonons even in systems wo inversion
 qptrlatt = 0; qptrlatt(1,1) = ngqpt(1); qptrlatt(2,2) = ngqpt(2); qptrlatt(3,3) = ngqpt(3)
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt1, nshiftq, shiftq, nqibz, qibz, wtq, nqbz, qbz)

 ABI_FREE(wtq)
 ABI_FREE(qbz)

 ! Build BZ mesh Note that:
 ! 1) q-point coordinates are in [0, 1]
 ! 2) The mesh is closed e.g. (0,0,0) and (1,1,1) are included
 nqx = ngqpt(1)+1; nqy = ngqpt(2)+1; nqz = ngqpt(3)+1

 ABI_MALLOC(xvec, (nqx))
 ABI_MALLOC(yvec, (nqy))
 ABI_MALLOC(zvec, (nqz))

 ! Multiple shifts are not supported here.
 do ix=1,nqx
   xvec(ix) = (ix-one+shiftq(1,1)) / ngqpt(1)
 end do
 do iy=1,nqy
   yvec(iy) = (iy-one+shiftq(2,1)) / ngqpt(2)
 end do
 do iz=1,nqz
   zvec(iz) = (iz-one+shiftq(3,1)) / ngqpt(3)
 end do

 ! Build list of q-points in full BZ (ordered as required by B-spline routines)
 nqfull = nqx*nqy*nqz
 ABI_MALLOC(qfull, (3,nqfull))
 iqf = 0
 do iz=1,nqz
   do iy=1,nqy
     do ix=1,nqx
       iqf = iqf + 1
       qfull(:,iqf) = [xvec(ix), yvec(iy), zvec(iz)]
     end do
   end do
 end do

 ! Build mapping qfull --> IBZ (q --> -q symmetry is always used)
 ABI_MALLOC(bz2ibz, (nqfull*sppoldbl1,6))

 call listkk(dksqmax,cryst%gmet,bz2ibz,qibz,qfull,nqibz,nqfull,cryst%nsym,&
   sppoldbl1,cryst%symafm,cryst%symrec,timrev1,use_symrec=.True.)
 ABI_FREE(qfull)

 if (dksqmax > tol12) then
   write(msg, '(3a,es16.6,4a)' )&
   'At least one of the q-points could not be generated from a symmetrical one.',ch10,&
   'dksqmax=',dksqmax,ch10,&
   'Action: check q-point input variables',ch10,&
   '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
   MSG_ERROR(msg)
 end if

 ! Generate knots (ords is input)
 qxord = ords(1); qyord = ords(2); qzord = ords(3)
 nxknot = nqx + qxord
 nyknot = nqy + qyord
 nzknot = nqz + qzord

 new%nqx = nqx; new%qxord = qxord
 new%nqy = nqy; new%qyord = qyord
 new%nqz = nqz; new%qzord = qzord

 ABI_MALLOC(new%xknot,(nxknot))
 ABI_MALLOC(new%yknot,(nyknot))
 ABI_MALLOC(new%zknot,(nzknot))

 call dbsnak(nqx, xvec, qxord, new%xknot)
 call dbsnak(nqy, yvec, qyord, new%yknot)
 call dbsnak(nqz, zvec, qzord, new%zknot)

 new%natom3 = 3 * cryst%natom

 ! Get phonon frequencies in IBZ
 ABI_CALLOC(ibz_freqs, (new%natom3, nqibz))
 cnt = 0
 do iq_ibz=1,nqibz
   cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! Mpi parallelism.
   call ifc_fourq(ifc, cryst, qibz(:,iq_ibz), ibz_freqs(:,iq_ibz), displ_cart)
 end do
 call xmpi_sum(ibz_freqs, comm, ierr)

 ABI_MALLOC(ibzdata_qnu, (nqibz, new%natom3))
 ibzdata_qnu = transpose(ibz_freqs)
 ABI_FREE(ibz_freqs)

 ABI_MALLOC(xyzdata,(nqx,nqy,nqz))
 ABI_DT_MALLOC(new%coeff, (new%natom3))

 do nu=1,new%natom3

   ABI_MALLOC(new%coeff(nu)%vals, (nqx,nqy,nqz))

   ! Build array in full bz to prepare call to dbs3in.
   iqf = 0
   do iz=1,nqz
     do iy=1,nqy
       do ix=1,nqx
         iqf = iqf + 1
         iq_ibz = bz2ibz(iqf,1)
         xyzdata(ix,iy,iz) = ibzdata_qnu(iq_ibz, nu)
       end do
     end do
   end do

   ! Construct 3D tensor for B-spline. Results in coeff(nu)%vals
   call dbs3in(nqx,xvec,nqy,yvec,nqz,zvec,xyzdata,nqx,nqy,qxord,qyord,qzord,new%xknot,new%yknot,new%zknot,&
      new%coeff(nu)%vals)
 end do

 ABI_FREE(xvec)
 ABI_FREE(yvec)
 ABI_FREE(zvec)
 ABI_FREE(bz2ibz)
 ABI_FREE(xyzdata)
 ABI_FREE(ibzdata_qnu)
 ABI_FREE(qibz)

end function ifc_build_phbspl
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/phbspl_evalq
!! NAME
!! phbspl_evalq
!!
!! FUNCTION
!!  Interpolate phonon frequencies at an arbitrary q-point.
!!
!! INPUTS
!!  qpt(3)=Q-point in reduced coordinate (will be wrapped in the interval [0,1[
!!
!! OUTPUT
!!  ofreqs(%natom3)=Interpolated phonon frequencies.
!!    Note that ofreqs is not necessarily sorted in ascending order.
!!    The routine does not reorder the interpolated frequencies
!!    to be consistent with the interpolation of the derivatives.
!!  [oder1(3,%natom3)]=First order derivatives.
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine phbspl_evalq(phbspl, qpt, ofreqs, oder1)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phbspl_evalq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(phbspl_t),intent(in) :: phbspl
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: ofreqs(phbspl%natom3)
 real(dp),optional,intent(out) :: oder1(3,phbspl%natom3)

!Local variables-------------------------------
!scalars
 integer :: nu,ii
!arrays
 integer :: iders(3)!, iperm(phbspl%natom3)
 real(dp) :: qred(3),shift(3)

! *********************************************************************

 ! Wrap k-point in the interval [0,1[ where 1 is not included (tol12)
 ! This is required because the spline has been constructed in this region.
 call wrap2_zero_one(qpt, qred, shift)

 do nu=1,phbspl%natom3
   ! B-spline interpolation.
   ofreqs(nu) = dbs3vl(qred(1), qred(2), qred(3), phbspl%qxord, phbspl%qyord, phbspl%qzord,&
                       phbspl%xknot, phbspl%yknot, phbspl%zknot, phbspl%nqx, phbspl%nqy, phbspl%nqz,&
                       phbspl%coeff(nu)%vals)
 end do

 ! Sort frequencies.
 !iperm = [(nu, nu=1, phbspl%natom3)]; call sort_dp(phbspl%natom3, ofreqs, iperm, tol14)

 if (present(oder1)) then
   ! Compute first-order derivatives.
   do nu=1,phbspl%natom3
     do ii=1,3
       iders = 0; iders(ii) = 1
       oder1(ii,nu) = dbs3dr(iders(1), iders(2), iders(3), &
                             qred(1), qred(2), qred(3), phbspl%qxord, phbspl%qyord, phbspl%qzord,&
                             phbspl%xknot, phbspl%yknot, phbspl%zknot, phbspl%nqx, phbspl%nqy, phbspl%nqz,&
                             phbspl%coeff(nu)%vals)
     end do
   end do

 end if

end subroutine phbspl_evalq
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/phbspl_free
!! NAME
!! phbspl_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine phbspl_free(phbspl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phbspl_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(phbspl_t),intent(inout) :: phbspl

!Local variables-------------------------------
!scalars
 integer :: ii

! *********************************************************************

 if (allocated(phbspl%xknot)) then
   ABI_FREE(phbspl%xknot)
 end if
 if (allocated(phbspl%yknot)) then
   ABI_FREE(phbspl%yknot)
 end if
 if (allocated(phbspl%zknot)) then
   ABI_FREE(phbspl%zknot)
 end if

 ! Free B-spline coefficients.
 if (allocated(phbspl%coeff)) then
   do ii=1,size(phbspl%coeff, dim=1)
     if (allocated(phbspl%coeff(ii)%vals)) then
       ABI_FREE(phbspl%coeff(ii)%vals)
     end if
   end do
   ABI_DT_FREE(phbspl%coeff)
 end if

end subroutine phbspl_free
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_build_skw
!! NAME
!! ifc_build_skw
!!
!! FUNCTION
!!
!! INPUTS
!!  ifc<type(ifc_type)>=Object containing the dynamical matrix and the IFCs.
!!  crystal<type(crystal_t)> = Information on the crystalline structure.
!!  ngqpt(3)=Divisions of the q-mesh used to produce the B-spline.
!!  nshiftq=Number of shifts in Q-mesh
!!  shiftq(3,nshiftq)=Shifts of the q-mesh.
!!  ords(3)=order of the spline for the three directions. ord(1) must be in [0, nqx] where
!!    nqx is the number of points along the x-axis.
!!  comm=MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(skw_t) function ifc_build_skw(ifc, cryst, ngqpt, nshiftq, shiftq, comm) result(new)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_build_skw'
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,nshiftq
 type(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: shiftq(3,nshiftq)

!local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1=1,timrev1=1,master=0,qptopt1=1
 integer :: nqibz,nqbz,iq_ibz,iq_bz,natom3,ierr,nu
 integer :: my_rank,nprocs,cnt
 real(dp) :: dksqmax
 character(len=500) :: msg
!arrays
 integer :: qptrlatt(3,3)
 integer,allocatable :: bz2ibz(:,:)
 real(dp) :: phfrq(3*cryst%natom),displ_cart(2,3,cryst%natom,3*cryst%natom)
 real(dp),allocatable :: ibz_freqs(:,:),wtq(:),qbz(:,:),qibz(:,:)

! *********************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 natom3 = 3 * cryst%natom

 ! Setup IBZ, weights and BZ. Always use q --> -q symmetry for phonons even in systems wo inversion
 qptrlatt = 0; qptrlatt(1,1) = ngqpt(1); qptrlatt(2,2) = ngqpt(2); qptrlatt(3,3) = ngqpt(3)
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt1, nshiftq, shiftq, nqibz, qibz, wtq, nqbz, qbz)

 ! Get phonon frequencies in IBZ
 ABI_CALLOC(ibz_freqs, (natom3, nqibz))
 cnt = 0
 do iq_ibz=1,nqibz
   cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! Mpi parallelism.
   call ifc_fourq(ifc, cryst, qibz(:,iq_ibz), ibz_freqs(:,iq_ibz), displ_cart)
 end do
 call xmpi_sum(ibz_freqs, comm, ierr)

 new = skw_new(cryst, 1, natom3, nqibz, 1, qibz, ibz_freqs, [0,0], [0,0], comm)

 if (.False. .and. my_rank == master) then
   ! Test whether SKW preserves symmetries.
   ! Build mapping qbz --> IBZ (q --> -q symmetry is always used)
   ABI_MALLOC(bz2ibz, (nqbz*sppoldbl1,6))

   call listkk(dksqmax,cryst%gmet,bz2ibz,qibz,qbz,nqibz,nqbz,cryst%nsym,&
     sppoldbl1,cryst%symafm,cryst%symrec,timrev1,use_symrec=.True.)

   if (dksqmax > tol12) then
     write(msg, '(3a,es16.6,4a)' )&
     'At least one of the q-points could not be generated from a symmetrical one.',ch10,&
     'dksqmax=',dksqmax,ch10,&
     'Action: check q-point input variables',ch10,&
     '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
     MSG_ERROR(msg)
   end if

   do iq_bz=1,nqbz
     iq_ibz = bz2ibz(iq_bz,1)
     do nu=1,natom3
       call skw_eval_bks(new, cryst, nu, qbz(:,iq_bz), 1, phfrq(nu))
     end do
     write(std_out,*)"BZ-IBZ:", maxval(abs(phfrq - ibz_freqs(:, iq_ibz)))
   end do
   ABI_FREE(bz2ibz)
 end if

 ABI_FREE(qibz)
 ABI_FREE(wtq)
 ABI_FREE(qbz)
 ABI_FREE(ibz_freqs)

end function ifc_build_skw
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_test_phinterp
!! NAME
!! ifc_test_phinterp
!!
!! INPUTS
!!  crystal<type(crystal_t)> = Information on the crystalline structure.
!!  ngqpt(3)=Divisions of the q-mesh used to produce the B-spline.
!!  nshiftq=Number of shifts in Q-mesh
!!  shiftq(3,nshiftq)=Shifts of the q-mesh.
!!  ords(3)=order of the spline for the three directions. ord(1) must be in [0, nqx] where
!!    nqx is the number of points along the x-axis.
!!  comm=MPI communicator
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_test_phinterp(ifc, cryst, ngqpt, nshiftq, shiftq, ords, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_test_phinterp'
!End of the abilint section

 implicit none

!arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,nshiftq
 type(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: ngqpt(3)
 integer,intent(in) :: ords(3)
 real(dp),intent(in) :: shiftq(3,nshiftq)

!local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iq,nq,natom3,my_rank,nprocs,ierr,mu
 real(dp) :: mare_bspl,mae_bspl,mare_skw,mae_skw
 real(dp) :: cpu,wall,gflops,cpu_fourq,wall_fourq,gflops_fourq
 real(dp) :: cpu_bspl,wall_bspl,gflops_bspl,cpu_skw,wall_skw,gflops_skw
 type(phbspl_t) :: phbspl
 type(skw_t) :: skw
!arrays
 real(dp) :: phfrq(3*cryst%natom),ofreqs(3*cryst%natom),qpt(3)
 real(dp) :: adiff_mev(3*cryst%natom), rel_err(3*cryst%natom)
 real(dp) :: displ_cart(2,3,cryst%natom,3*cryst%natom)
 real(dp) :: qred(3),shift(3),vals4(4)

! *********************************************************************

 !call make_path(Kpath%nbounds,bounds,Kpath%gmet,"G",ndiv_small,Kpath%ndivs,Kpath%npts,pts,unit=dev_null)
 !call kpath_init(Kpath,bounds,gprimd,ndiv_small)

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 natom3 = 3 * cryst%natom

 ! Build interpolator objects.
 phbspl = ifc_build_phbspl(ifc, cryst, ngqpt, nshiftq, shiftq, ords, comm)
 skw = ifc_build_skw(ifc, cryst, ngqpt, nshiftq, shiftq, comm)

 cpu_fourq = zero; wall_fourq = zero; gflops_fourq = zero
 cpu_bspl = zero; wall_bspl = zero; gflops_bspl = zero
 cpu_skw = zero; wall_skw = zero; gflops_skw = zero
 mae_bspl = zero; mare_bspl = zero
 mae_skw = zero; mare_skw = zero

 nq = 1000
 !call random_seed(put=my_rank*100)
 do iq=1,nq
   if (mod(iq, nprocs) /= my_rank) cycle ! mpi parallelism
   !call random_number(qpt)
   call random_number(qred)
   qred = qred + (0.1_dp * my_rank / nprocs)
   !call wrap2_zero_one(qred, qpt, shift)
   call wrap2_pmhalf(qred, qpt, shift)

   ! Fourier interpolation
   call cwtime(cpu, wall, gflops, "start")
   call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart)
   call cwtime(cpu, wall, gflops, "stop")
   cpu_fourq = cpu_fourq + cpu; wall_fourq = wall_fourq + wall

   ! B-spline interpolation
   call cwtime(cpu, wall, gflops, "start")
   call phbspl_evalq(phbspl, qpt, ofreqs)
   call cwtime(cpu, wall, gflops, "stop")
   cpu_bspl = cpu_bspl + cpu; wall_bspl = wall_bspl + wall

   adiff_meV = abs(phfrq - ofreqs); rel_err = zero
   where (abs(phfrq) > tol16)
     rel_err = adiff_meV / abs(phfrq)
   end where
   rel_err = 100 * rel_err; adiff_meV = adiff_meV * Ha_meV
   mae_bspl = mae_bspl + sum(adiff_meV) / natom3
   mare_bspl = mare_bspl + sum(rel_err) / natom3
   !if (my_rank == master)
   !write(std_out,*)maxval(adiff_meV), maxval(rel_err), trim(ktoa(qpt))
   !write(456, *)phfrq; write(457, *)ofreqs
   !end if

   ! SKW interpolation
   call cwtime(cpu, wall, gflops, "start")
   do mu=1,natom3
     call skw_eval_bks(skw, cryst, mu, qpt, 1, ofreqs(mu))
   end do
   call cwtime(cpu, wall, gflops, "stop")
   cpu_skw = cpu_skw + cpu; wall_skw = wall_skw + wall

   adiff_meV = abs(phfrq - ofreqs); rel_err = zero
   where (abs(phfrq) > tol16)
     rel_err = adiff_meV / abs(phfrq)
   end where
   rel_err = 100 * rel_err; adiff_meV = adiff_meV * Ha_meV
   mae_skw = mae_skw + sum(adiff_meV) / natom3
   mare_skw = mare_skw + sum(rel_err) / natom3
   !if (my_rank == master) then
   !write(std_out,*)maxval(adiff_meV), maxval(rel_err), trim(ktoa(qpt))
   !write(456, *)phfrq; write(457, *)ofreqs
   !endif
 end do

 vals4 = [mae_bspl, mare_bspl, mae_skw, mare_skw]
 call xmpi_sum(vals4, comm, ierr)
 mae_bspl = vals4(1); mare_bspl = vals4(2); mae_skw = vals4(3); mare_skw = vals4(4)

 mae_bspl = mae_bspl / nq; mare_bspl = mare_bspl / nq
 mae_skw = mae_skw / nq; mare_skw = mare_skw / nq

 if (my_rank == master) then
   write(std_out,"(2(a,f6.2),a,i0)")"B-spline MAE: ",mae_bspl," [meV], MARE: ",mare_bspl, "% with nqpt: ",nq
   write(std_out,"(2(a,f6.2),a,i0)")"SKW      MAE: ",mae_skw, " [meV], MARE: ",mare_skw, "% with nqpt: ",nq
   write(std_out,"(2(a,f6.2))")"fourq: cpu: ",cpu_fourq,", wall: ",wall_fourq
   write(std_out,"(2(a,f6.2))")"bspl:  cpu: ",cpu_bspl,", wall: ",wall_bspl
   write(std_out,"(2(a,f6.2))")"skw:   cpu: ",cpu_skw,", wall: ",wall_skw
 end if

 call phbspl_free(phbspl)
 call skw_free(skw)

end subroutine ifc_test_phinterp
!!***

!----------------------------------------------------------------------

end module m_ifc
