!!****m* ABINIT/defs_datatypes
!! NAME
!! defs_datatypes
!!
!! FUNCTION
!! This module contains definitions important structured datatypes for the ABINIT package.
!! If you want to add one new datatype, please, examine first whether
!! another datatype might meet your need (e.g. adding some records to it).
!! Then, if you are sure your new structured datatype is needed,
!! write it here, and DOCUMENT it properly (not all datastructure here are
!! well documented, it is a shame ...).
!! Do not forget: you will likely be the major winner if you document properly.
!! Proper documentation of a structured datatype means:
!!
!!  (1) Mention it in the list just below
!!  (2) Describe it in the NOTES section
!!  (3) Put it in alphabetical order in the the main section of this module
!!  (4) Document each of its records, except if they are described elsewhere
!!      (this exception is typically the case of the dataset associated with
!!      input variables, for which there is a help file)
!!
!! List of datatypes:
!! * ebands_t: different information about the band structure
!! * pseudopotential_type: for norm-conserving pseudopotential, all the information
!! * pspheader_type: for norm-conserving pseudopotentials, the header of the file
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module defs_datatypes

 use defs_basis

 implicit none
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/ebands_t
!! NAME
!! ebands_t
!!
!! FUNCTION
!! It contains different information about the band structure: eigenenergies, residuals, derivative of occupation
!! numbers vs energy in case of metallic occupations and Brillouin zone according to the context: k points,
!! occupation numbers, storage mode of wavefunctions, weights ...
!! For example, the initial Brillouin zone, set up in the dataset, will be treated in the response function part of
!! the code, to give a reduced Brillouin zone different from the original one, due to the breaking of the symmetries
!! related to the existence of a wavevector, or the lack of time-reversal invariance.
!!
!! SOURCE

 type ebands_t

  integer :: bantot                ! total number of bands (sum(nband(:))
  integer :: mband                 ! Max number of bands i.e MAXVAL(nband) (to dimension arrays)
  integer :: nkpt                  ! number of k points
  integer :: nspinor               ! 1 for collinear, 2 for noncollinear.
  integer :: nsppol                ! number of spin-polarizations
  integer :: ntemp                 ! Number of temperatures
  integer :: occopt                ! Occupation option, see input variable.

  real(dp) :: entropy              ! Entropy associated with the smearing (adimensional)
  real(dp) :: fermie               ! Fermi energy
  real(dp) :: nelect               ! Number of electrons.
  real(dp) :: tphysel              ! Physical temperature of electrons.
  real(dp) :: tsmear               ! Temperature of smearing.

  !real(dp) :: spinmagntarget
  ! TODO This should be set via dtset%spinmagntarget to simplify the API.

  integer,allocatable :: istwfk(:)
  ! istwfk(nkpt)
  ! Storage mode at each k point.

  integer,allocatable :: nband(:)
  ! nband(nkpt*nsppol)
  ! Number of bands at each k point and spin-polarisation.

  integer,allocatable :: npwarr(:)
  ! npwarr(nkpt)
  ! Number of plane waves at each k point.

  real(dp),allocatable :: kptns(:,:)
  ! kptns(3,nkpt)
  ! k-point vectors.

  real(dp),allocatable :: eig(:,:,:)
  ! eig(mband,nkpt,nsppol)
  ! Eigenvalues of each band.

  real(dp),allocatable :: linewidth(:,:,:,:)
  ! linewidth(itemp,mband,nkpt,nsppol)
  ! Linewidth of each band
  ! MG: TODO: This array should be removed (I think Yannick introduced it, see also Ktmesh)

  !real(dp),allocatable :: kTmesh(:)
  ! kTmesh(ntemp)
  ! List of temperatures (KT units).

  !real(dp),allocatable :: velocity(:,:,:,:)
  ! velocity(3,mband,nkpt,nsppol)
  ! Group velocity of each band
  ! MG: TODO: This array should be removed (I think HM introduced it)

  real(dp),allocatable :: occ(:,:,:)
  ! occ(mband, nkpt, nsppol)
  ! occupation of each band.

  real(dp),allocatable :: doccde(:,:,:)
  ! doccde(mband, nkpt, nsppol)
  ! derivative of the occupation of each band wrt energy (needed for RF).

  real(dp),allocatable :: wtk(:)
  ! wtk(nkpt)
  ! weight of each k point, normalized to one.

  integer :: kptopt
  ! Option used for k-point generation.

  integer :: nshiftk_orig, nshiftk
  ! original number of shifts given in input and the actual value (changed in inkpts)

  real(dp) :: charge
  ! nelect = zion - charge
  ! Extra charge added to the unit cell when performing GS calculations
  ! To treat a system missing one electron per unit cell, charge is set to +1.
  ! When reading the band structure from an external file,
  ! charge is mainly used as metadata describing the GS calculation that procuded the ebands_t object.
  ! To simulate doping in a post-processing tool, use ebands_set_extrael that defines the value of %extra_el.
  ! and changes %nelect, accordingly.

  real(dp) :: extrael = zero
  ! Extra number of electrons.
  ! This variable is mainly used to simulate doping in the rigid band approximation.
  ! Set by ebands_set_extrael method.

  integer :: kptrlatt_orig(3,3), kptrlatt(3,3)
  ! Original value of kptrlatt and value after the call to inkpts

  real(dp),allocatable :: shiftk_orig(:,:)
  ! shiftk_orig(3, nshiftk_orig)
  ! original shifts given in input (changed in inkpts).

  real(dp),allocatable :: shiftk(:,:)
  ! shiftk(3, nshiftk)

 end type ebands_t
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pseudopotential_gth_type
!! NAME
!! pseudopotential_gth_type
!!
!! FUNCTION
!! This structure is a sub-structure of pseudopotential_type used to
!! store parameters from the GTH pseudo-potentials. All arrays have
!! indices running on 1:npsp for each read pseudo-file. The 'set' array
!! is a check array, since several different pseudo can be used in a simulation
!! it set a flag for each npsp if params have been set or not. This is
!! redundant with psps%pspcod in the way that when psps%pspcod(i) is 2,
!! then gth_params%set(i) is .true.. GTH pseudo previous to wavelets introduction
!! doesn't have geometric information. These have been added on the last line.
!! It is three radius information, the %hasGeometry flag is there to know
!! which kind of pseudo has been read.
!!
!! SOURCE

 type pseudopotential_gth_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  real(dp), allocatable :: psppar(:, :, :)
   ! These are {rloc, C(1...4)} coefficients for psppar(0, :, :) indices,
   ! Followed by the h coefficients for psppar(1:2, :, :) indices.
   !  size (0:2, 0:4, npsp)

  real(dp), allocatable :: radii_cf(:, :)
   ! Cut-off radii for core part and long-range part.
   ! radii_cf(:, 1) is for the long-range cut-off and
   ! radii_cf(:, 2) is for the core cut-off. size (npsp, 2)

  real(dp), allocatable :: psp_k_par(:, :, :)
   ! Spin orbit coefficients in HGH/GTH formats: k11p etc... see psp3ini.F90
   !   dimension = num l channels, 3 coeffs, num psp = (1:lmax+1,1:3,npsp)

  logical, allocatable :: hasGeometry(:)
   ! Flag for geometric information in the pseudo. size (npsp)

  logical, allocatable :: set(:)
   ! Consistency array, used for checking size (npsp)

 end type pseudopotential_gth_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/nctab_t
!! NAME
!! nctab_type
!!
!! FUNCTION
!!  This object contains TABulated data for NC pseudopotentials
!!  It follows the conventions used in pawtab so that we can reuse
!!  the PAW routines for the treatment of model core change as well
!!  as the code in the atm2fft routines used to build approximated densities
!!  from atomic quantities. Methods are defined in m_psps.
!!
!! SOURCE

 type,public :: nctab_t

   integer :: mqgrid_vl = 0
   ! Number of points in the reciprocal space grid on which
   ! the radial functions are specified (same grid as the one used for the local part).

   logical :: has_tvale = .False.
    ! True if the norm-conserving pseudopotential provides the atomic pseudized valence density.
    ! If alchemy, has_tvale is True only if all the mixed pseudos
    ! have the valence charge in the pseudopotential file.

   logical :: has_tcore = .False.
    ! True if the norm-conserving pseudopotential has the model core-charge for NLCC.
    ! If alchemy, has_tcore is True if at least one of the mixed pseudos has NLCC.
    ! See also tcorespl

   real(dp) :: dncdq0 = zero
    ! Gives 1/q d(tNcore(q))/dq for q=0
    ! (tNcore(q) = FT of pseudo core density)

   real(dp) :: d2ncdq0 = zero
    ! Gives contribution of d2(tNcore(q))/d2q for q=0
    ! \int{(16/15)*pi^5*n(r)*r^6* dr}
    ! (tNcore(q) = FT of pseudo core density)

   real(dp) :: dnvdq0 = zero
    ! Gives 1/q d(tNvale(q))/dq for q=0
    ! (tNvale(q) = FT of pseudo valence density)

   real(dp), allocatable :: tvalespl(:,:)
    ! tvalespl(mqgrid_vl,2)
    ! Gives the pseudo valence density in reciprocal space on a regular grid
    ! Initialized only if has_nctvale(itypat)

   real(dp), allocatable :: tcorespl(:,:)
    ! tcorespl(mqgrid_vl,2)
    ! Gives the pseudo core density in reciprocal space on a regular grid.
    ! tcorespl is **always** allocated and initialized with zeros if not has_tcore
    ! A similar approach is used in PAW.

 end type nctab_t
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pseudopotential_type
!! NAME
!! pseudopotential_type
!!
!! FUNCTION
!! This structured datatype contains all the information about one
!! norm-conserving pseudopotential, including the description of the local
!! and non-local parts, the different projectors, the non-linear core
!! correction ...
!!
!! SOURCE

 type pseudopotential_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.


! Integer scalars
  integer :: dimekb
   ! Dimension of Ekb
   ! ->Norm conserving : Max. number of Kleinman-Bylander energies
   !                     for each atom type
   !                     dimekb=lnmax (lnmax: see this file)
   ! ->PAW : Max. number of Dij coefficients connecting projectors
   !                     for each atom type
   !                     dimekb=lmnmax*(lmnmax+1)/2 (lmnmax: see this file)

  integer :: lmnmax
   !  If useylm=0, max number of (l,n)   comp. over all type of psps (lnproj)
   !  If useylm=1, max number of (l,m,n) comp. over all type of psps (lmnproj)
   !  If mpspso is 2, lmnmax takes into account the spin-orbit projectors,
   !  so, it is equal to the max of lmnprojso or lnprojso, see pspheader_type

  integer :: lnmax
   !  Max. number of (l,n) components over all type of psps
   !  If mpspso is 2, lmnmax takes into account the spin-orbit projectors,
   !  so, it is equal to the max of lnprojso, see pspheader_type

  integer :: mproj
   ! Maximum number of non-local projectors over all angular momenta and type of psps
   ! 0 only if all psps are local

  integer :: mpsang
   ! Highest angular momentum of non-local projectors over all type of psps.
   ! shifted by 1 : for all local psps, mpsang=0; for largest s, mpsang=1,
   ! for largest p, mpsang=2; for largest d, mpsang=3; for largest f, mpsang=4
   ! This gives also the number of non-local "channels"

  integer :: mpspso
   ! mpspso is set to 1 if none of the psps is used with a spin-orbit part (that
   !  is, if the user input variable so_psp is not equal
   !  to 1 in at least one case
   ! otherwise, it is set to 2

  integer :: mpssoang
   ! Maximum number of channels, including those for treating the spin-orbit coupling
   ! when mpspso=1, mpssoang=mpsang
   ! when mpspso=2, mpssoang=2*mpsang-1

  integer :: mqgrid_ff
   ! Number of points in the reciprocal space grid on which
   ! the radial functions ffspl are specified

  integer :: mqgrid_vl
   ! Number of points in the reciprocal space grid on which
   ! the radial functions vlspl are specified

  integer :: mtypalch
   ! Maximum number of alchemical pseudo atoms. If non-zero,
   ! the mechanism to generate mixing of pseudopotentials is activated

  integer :: npsp
   ! Number of types of pseudopotentials

  integer :: npspalch
   ! Number of types of pseudopotentials used for alchemical purposes

  integer :: ntypat
   ! Number of types of atoms (might be alchemy wrt pseudopotentials)

  integer :: ntypalch
   ! Number of types of alchemical pseudoatoms

  integer :: ntyppure
   ! Number of types of pure pseudoatoms

  integer :: n1xccc
   ! Number of radial points for the description of the pseudo-core charge
   ! (in the framework of the non-linear XC core correction)

  integer :: optnlxccc
   ! Option for the choice of non-linear XC core correction treatment (see the input variable)
   ! Used only for FHI pseudos.

  integer :: positron
   ! Option for the choice of type of GS calculation (electron or positron)

  integer :: usepaw
   ! if usepaw=0 , use norm-conserving psps part of the code
   ! is usepaw=1 , use paw part of the code

  integer :: usewvl
   ! if usewvl=0 ,  plane waves
   ! is usewvl=1 ,  wavelets

  integer :: useylm
   ! governs the way the nonlocal operator is to be applied:
   !   1=using Ylm, 0=using Legendre polynomials

! Logical scalars
  logical :: vlspl_recipSpace
   ! governs if vlspl is compute in reciprocal space or in real
   ! space (when available).

! Integer arrays
  integer, allocatable :: algalch(:)
   ! algalch(ntypalch)
   ! For each type of pseudo atom, the algorithm to mix the pseudopotentials

  integer, allocatable :: indlmn(:,:,:)
   ! indlmn(6,lmnmax,ntypat)
   ! For each type of psp,
   ! array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
   !                                or i=lmn (if useylm=1)
   ! NB: spin is used for NC pseudos with SOC term: 1 if scalar term (spin diagonal), 2 if SOC term.

  integer, allocatable :: pspdat(:)
   ! pspdat(ntypat)
   ! For each type of psp, the date of psp generation, as given by the psp file

  integer, allocatable :: pspcod(:)
   ! pspcod(npsp)
   ! For each type of psp, the format -or code- of psp generation,
   !  as given by the psp file

  integer, allocatable :: pspso(:)
   ! pspso(npsps)
   ! For each type of psp,
   ! 1 if no spin-orbit component is taken into account
   ! 2 if a spin-orbit component is used

  integer, allocatable :: pspxc(:)
   ! pspxc(npsps)
   ! For each type of psp, the XC functional that was used to generate it,
   ! as given by the psp file

! Real (real(dp)) arrays

  real(dp), allocatable :: ekb(:,:)
   ! ekb(dimekb,ntypat*(1-usepaw))
   !  ->NORM-CONSERVING PSPS ONLY:
   !    (Real) Kleinman-Bylander energies (hartree)
   !           for number of basis functions (l,n) (lnmax)
   !           and number of atom types (ntypat)
   ! NOTE (MT) : ekb (norm-conserving) is now diagonal (one dimension
   !             lnmax); it would be easy to give it a second
   !             (symmetric) dimension by putting
   !             dimekb=lnmax*(lnmax+1)/2
   !             in the place of dimekb=lmnmax.

  real(dp), allocatable :: ffspl(:,:,:,:)
   ! ffspl(mqgrid_ff,2,lnmax,ntypat)
   ! Gives, on the radial grid, the different non-local projectors,
   ! in both the norm-conserving case, and the PAW case

  real(dp), allocatable :: mixalch(:,:)
   ! mixalch(npspalch,ntypalch)
   ! Mixing coefficients to generate alchemical pseudo atoms

  real(dp), allocatable :: qgrid_ff(:)
   ! qgrid_ff(mqgrid_ff)
   ! The coordinates of all the points of the radial grid for the nl form factors

  real(dp), allocatable :: qgrid_vl(:)
   ! qgrid_vl(mqgrid_vl)
   ! The coordinates of all the points of the radial grid for the local part of psp

  real(dp), allocatable :: vlspl(:,:,:)
   ! vlspl(mqgrid_vl,2,ntypat)
   ! Gives, on the radial grid, the local part of each type of psp.

  real(dp), allocatable :: dvlspl(:,:,:)
   ! dvlspl(mqgrid_vl,2,ntypat)
   ! Gives, on the radial grid, the first derivative of the local
   ! part of each type of psp (computed when the flag 'vlspl_recipSpace' is true).

  real(dp), allocatable :: xcccrc(:)
   ! xcccrc(ntypat)
   ! Gives the maximum radius of the pseudo-core charge, for each type of psp.

  real(dp), allocatable :: xccc1d(:,:,:)
   ! xccc1d(n1xccc*(1-usepaw),6,ntypat)
   ! Norm-conserving psps only
   ! The component xccc1d(n1xccc,1,ntypat) is the pseudo-core charge
   ! for each type of atom, on the radial grid. The components
   ! xccc1d(n1xccc,ideriv,ntypat) give the ideriv-th derivative of the
   ! pseudo-core charge with respect to the radial distance.

  real(dp), allocatable :: zionpsp(:)
   ! zionpsp(npsp)
   ! For each pseudopotential, the ionic pseudo-charge
   ! (giving raise to a long-range coulomb potential)

  real(dp), allocatable :: ziontypat(:)
   ! ziontypat(ntypat)
   ! For each type of atom (might be alchemy wrt psps), the ionic pseudo-charge
   ! (giving raise to a long-range coulomb potential)

  real(dp), allocatable :: znuclpsp(:)
   ! znuclpsp(npsp)
   ! The atomic number of each pseudopotential

  real(dp), allocatable :: znucltypat(:)
   ! znucltypat(ntypat)
   ! The atomic number of each type of atom (might be alchemy wrt psps)

! Character arrays
  character(len=fnlen), allocatable :: filpsp(:)
   ! filpsp(npsps)
   ! The filename of the pseudopotential

  character(len=fnlen), allocatable :: title(:)
   ! title(npsp)
   ! The content of first line read from the psp file

  character(len=md5_slen), allocatable :: md5_pseudos(:)
  ! md5pseudos(npsp)
  ! md5 checksums associated to pseudos (read from file)

  type(pseudopotential_gth_type) :: gth_params
   ! Types for pseudo-potentials that are based on parameters. Currently, only
   ! GTH are supported (see pseudopotential_gth_type). To add one, one should
   ! create an initialisation method and a destruction method in 02psp (see
   ! psp2params.F90). These methods are called in driver().

   type(nctab_t),allocatable :: nctab(:)
   ! nctab(ntypat)
   ! Tables storing data for NC pseudopotentials.

   integer :: nc_xccc_gspace = 0
   ! NC pseudos only. Set to 1 if the non-linear core correction should
   ! be treated in G-space similarly to the approach used for PAW.

 end type pseudopotential_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pspheader_paw_type
!! NAME
!! pspheader_paw_type
!!
!! FUNCTION
!! The pspheader_paw_type structured datatype gather additional information
!! about a PAW pseudopotential file, from its header.
!!
!! SOURCE

 type pspheader_paw_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! WARNING : Also pay attention to subroutine pspheads_comm, which broadcasts this datatype.

  integer :: basis_size    ! Number of elements of the wf basis ((l,n) quantum numbers)
  integer :: l_size        ! Maximum value of l+1 leading to a non zero Gaunt coefficient
  integer :: lmn_size      ! Number of elements of the paw basis
  integer :: mesh_size     ! Dimension of (main) radial mesh
  integer :: pawver        ! Version number of paw psp format
  integer :: shape_type    ! Type of shape function
  real(dp) :: rpaw         ! Radius for paw spheres
  real(dp) :: rshp         ! Cut-off radius of shape function

 end type pspheader_paw_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pspheader_type
!! NAME
!! pspheader_type
!!
!! FUNCTION
!! The pspheader_type structured datatype gather different information
!! about a pseudopotential file, from its header.
!!
!! SOURCE

 type pspheader_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.
! WARNING : Also pay attention to subroutine pspheads_comm, which broadcasts this datatype.

  integer :: nproj(0:3) ! number of scalar projectors for each angular momentum

  integer :: nprojso(3) ! number of spin-orbit projectors for each angular momentum

  integer :: lmax       ! maximum l quantum number (-1 if only local)
                        ! Example : s only       -> lmax=0
                        !           s and p      -> lmax=1
                        !           s,p,d        -> lmax=2

  integer :: pspcod
    ! code number of the pseudopotential

  integer :: pspdat
    ! date of generation of the pseudopotential

  integer :: pspxc
    ! exchange-correlation functional

  integer :: pspso
    ! spin-orbit characteristics
    ! 0 if pseudo does not support SOC, 1 or 2 if SOC terms are available in the pp file.
    ! Note that one could have a pseudo with SOC terms but ignore the SOC contribution
    ! For example, it's possible to use nspinor=2 and set so_psp to 0 in the input file
    ! or perform a run with nspinor=1 with pseudos containing SOC terms.

  integer :: usewvl
   ! if usewvl=0 ,  plane waves
   ! is usewvl=1 ,  wavelets

  integer :: xccc
    ! =0 if no XC core correction, non-zero if XC core correction

  real(dp) :: zionpsp
    ! charge of the ion made of core electrons only

  real(dp) :: znuclpsp
    ! atomic number of the nuclei

  real(dp) :: GTHradii(0:4)
    ! Radii values for GTH (and HGH) family potentials

  character(len=fnlen) :: filpsp
    ! name of the psp file

  character(len=fnlen) :: title
    ! content of first line read from the psp file

  character(len=md5_slen) :: md5_checksum = md5_none
    ! md5 checksum read from file.

  type(pspheader_paw_type) :: pawheader
    ! only for PAW psps ; see m_pawpsp.

 end type pspheader_type
!!***

end module defs_datatypes
!!***

