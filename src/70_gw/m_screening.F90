!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_screening
!! NAME
!!  m_screening
!!
!! FUNCTION
!!  This module contains the definition of the object used to deal
!!  with the inverse dielectric matrix as well as related methods.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_screening

 use defs_basis
 use defs_abitypes
 use m_hide_blas
 use m_linalg_interfaces
 use m_xmpi
 use m_errors
 use m_copy
 use m_splines
 use m_abicore
 use m_lebedev
 use m_spectra
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_gwdefs,          only : GW_TOLQ0, czero_gw, GW_Q0_DEFAULT
 use m_fstrings,        only : toupper, endswith, sjoin, itoa
 use m_io_tools,        only : open_file
 use m_numeric_tools,   only : print_arr, hermitianize
 use m_special_funcs,   only : k_fermi, k_thfermi
 use m_geometry,        only : normv, vdotw, metric
 use m_hide_lapack,     only : xginv
 use m_crystal,         only : crystal_t
 use m_bz_mesh,         only : kmesh_t, get_BZ_item, box_len
 use m_fft_mesh,        only : g2ifft
 use m_fftcore,         only : kgindex
 use m_fft,             only : fourdp
 use m_gsphere,         only : gsphere_t
 use m_vcoul,           only : vcoul_t
 use m_io_screening,    only : hscr_free, hscr_io, hscr_print, hscr_from_file, read_screening, write_screening, &
&                              hscr_copy, HSCR_LATEST_HEADFORM, hscr_t, ncname_from_id, em1_ncname
 use m_paw_sphharm,     only : ylmc
 use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_screening/epsilonm1_results
!! NAME
!! epsilonm1_results
!!
!! FUNCTION
!! For the GW part of ABINIT, the epsilonm1_results structured datatype
!! gather the results of screening: the inverse dielectric matrix, and the omega matrices.
!!
!! SOURCE

 type,public :: Epsilonm1_results

  integer :: id
  ! Matrix identifier: O if not yet defined, 1 for chi0,
  ! 2 for chi, 3 for epsilon, 4 for espilon^{-1}, 5 for W.

  integer :: ikxc
  ! Kxc kernel used, 0 for None (RPA), >0 for static TDDFT (=ixc), <0 for TDDFT

  integer :: fform
  ! File format: 1002 for SCR|SUSC files.

  integer :: mqmem
  ! =0 for out-of-core solution, =nqibz if entire matrix is stored in memory.

  integer :: nI,nJ
  ! Number of components (rows,columns) in chi|eps^-1. (1,1) if collinear.

  integer :: nqibz
  ! Number of q-points in the IBZ used.

  integer :: nqlwl
  ! Number of point used for the treatment of the long wave-length limit.

  integer :: nomega
  ! Total number of frequencies.

  integer :: nomega_i
  ! Number of purely imaginary frequencies used.

  integer :: nomega_r
  ! Number of real frequencies used.

  integer :: npwe
  ! Number of G vectors.

  integer :: test_type
  ! 0 for None, 1 for TEST-PARTICLE, 2 for TEST-ELECTRON (only for TDDFT)

  integer :: tordering
  ! 0 if not defined, 1 for Time-Ordered, 2 for Advanced, 3 for Retarded.

  character(len=fnlen) :: fname
  ! Name of the file from which epsm1 is read.

  integer,allocatable :: gvec(:,:)
  ! gvec(3,npwe)
  ! G-vectors used to describe the two-point function (r.l.u.).

  real(dp),allocatable :: qibz(:,:)
  ! qibz(3,nqibz)
  ! q-points in reduced coordinates

  real(dp),allocatable :: qlwl(:,:)
  ! qlwl(3,nqlwl)
  ! q-points used for the long wave-length limit treatment.

  !real(gwp),allocatable :: epsm1_pole(:,:,:,:)
  ! epsm1(npwe,npwe,nomega,nqibz)
  ! Contains the two-point function $\epsilon_{G,Gp}(q,omega)$ in frequency and reciprocal space.

  complex(gwpc),allocatable :: epsm1(:,:,:,:)
  ! epsm1(npwe,npwe,nomega,nqibz)
  ! Contains the two-point function $\epsilon_{G,Gp}(q,omega)$ in frequency and reciprocal space.

  complex(dpc),allocatable :: omega(:)
  ! omega(nomega)
  ! Frequencies used both along the real and the imaginary axis.

  type(hscr_t) :: Hscr
  ! The header reported in the _SCR of _SUSC file.
  ! This object contains information on the susceptibility or the inverse dielectric matrix
  ! as stored in the external file. These quantities do *NOT* correspond to the quantities
  ! used during the GW calculation since some parameters might differ, actually they might be smaller.
  ! For example, the number of G-vectors used can be smaller than the number of G"s stored on file.

 end type Epsilonm1_results

 public :: em1results_free                ! Free memory
 public :: em1results_print               ! Print basic info
 public :: Epsm1_symmetrizer              ! Symmetrize two-point function at a q-point in the BZ.
 public :: Epsm1_symmetrizer_inplace      ! In-place version of the above
 public :: init_Er_from_file              ! Initialize the object from file
 public :: mkdump_Er                      ! Dump the object to a file.
 public :: get_epsm1
 public :: decompose_epsm1
 public :: make_epsm1_driver              !  Calculate the inverse symmetrical dielectric matrix starting from chi0
 public :: mkem1_q0                       ! construct the microscopic dieletric matrix for q-->0
 public :: screen_mdielf                  ! Calculates W_{G,G'}(q,w) for a given q-point in the BZ using a model dielectric function.
 public :: rpa_symepsm1
!!***


!!****t* m_screening/chi_t
!! NAME
!!  chi_t
!!
!! FUNCTION
!!  This object contains the head and the wings of the polarizability
!!  These quantities are used to treat the q-->0 limit
!!
!! SOURCE

 type,public :: chi_t

   integer :: npwe
   ! Number of G vectors.

   integer :: nomega
   ! Number of frequencies

   complex(gwpc),allocatable :: mat(:,:,:)
   ! mat(npwe, npwe, nomega)

   complex(dpc),allocatable :: head(:,:,:)
   ! head(3,3,nomega)

   complex(dpc),allocatable :: lwing(:,:,:)
   ! lwing(3,npwe,nomega)
   ! Lower wings

   complex(dpc),allocatable :: uwing(:,:,:)
   ! uwing(3,npwe,nomega)
   ! Upper wings.

 end type chi_t

 public :: chi_new    ! Create new object (allocate memory)
 public :: chi_free   ! Free memory.
!!***

!!****t* m_screening/lwl_t
!! NAME
!!  lwl_t
!!
!! FUNCTION
!!
!! SOURCE

 type,public :: lwl_t

  integer :: npwe
  ! Number of G vectors.

  integer :: nomega
  ! Number of frequencies

  integer :: method
  ! 1 = Only head
  ! 2 = head + wings
  ! 3 = head + wings + body corrections.

  character(len=fnlen) :: fname
   ! Name of the file from which epsm1 is read.

   complex(dpc),allocatable :: head(:,:,:)
   ! head(3,3,nomega)

   complex(dpc),allocatable :: lwing(:,:,:)
   ! lwing(3,npwe,nomega)
   ! Lower wings

   complex(dpc),allocatable :: uwing(:,:,:)
   ! uwing(3,npwe,nomega)
   ! Upper wings.

   complex(dpc),allocatable :: body(:,:,:)
   ! uwing(npwe,npwe,nomega)
   ! Body terms

 end type lwl_t

 public :: lwl_write
 public :: lwl_init
 !public :: lwl_from_file
 public :: lwl_free
!!***


CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_screening/em1results_free
!! NAME
!! em1results_free
!!
!! FUNCTION
!! Deallocate all the pointers in Er that result to be associated.
!! Perform also a cleaning of the Header.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_screening,mrgscr,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine em1results_free(Er)

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_results),intent(inout) :: Er
! *************************************************************************

 !@Epsilonm1_results
 !integer
 if (allocated(Er%gvec)) then
   ABI_FREE(Er%gvec)
 end if

 !real
 if (allocated(Er%qibz)) then
   ABI_FREE(Er%qibz)
 end if
 if (allocated(Er%qlwl)) then
   ABI_FREE(Er%qlwl)
 end if

 !complex
 if (allocated(Er%epsm1)) then
   ABI_FREE(Er%epsm1)
 end if
 if (allocated(Er%omega)) then
   ABI_FREE(Er%omega)
 end if

 !datatypes
 call hscr_free(Er%Hscr)

end subroutine em1results_free
!!***

!----------------------------------------------------------------------

!!****f* m_screening/em1results_print
!! NAME
!!  em1results_print
!!
!! FUNCTION
!! Print the basic dimensions and the most important
!! quantities reported in the Epsilonm1_results data type.
!!
!! INPUTS
!!  Er<Epsilonm1_results>=The data type.
!!  unit[optional]=the unit number for output.
!!  prtvol[optional]=verbosity level.
!!  mode_paral[optional]=either COLL or PERS.
!!
!! OUTPUT
!!  Only printing.
!!
!! PARENTS
!!      m_screening,mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine em1results_print(Er,unit,prtvol,mode_paral)

!Arguments ------------------------------------
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 type(Epsilonm1_results),intent(in) :: Er

!Local variables-------------------------------
 integer :: iw,iqibz,iqlwl,unt,my_prtvol
 character(len=50) :: rfname,rforder,rfapprox,rftest,kxcname
 character(len=500) :: msg
 character(len=4) :: mode
! *************************************************************************

 unt = std_out; if (present(unit)) unt = unit
 my_prtvol = 0; if (present(prtvol)) my_prtvol = prtvol
 mode   ='COLL'; if (present(mode_paral)) mode = mode_paral

 ! === chi0 or \epsilon^{-1} ? ===
 SELECT CASE (Er%ID)
 CASE (0)
   rfname='Undefined'
 CASE (1)
   rfname='Irreducible Polarizability'
 CASE (2)
   rfname='Polarizability'
 CASE (3)
   rfname='Symmetrical Dielectric Matrix'
 CASE (4)
   rfname='Symmetrical Inverse Dielectric Matrix'
 CASE DEFAULT
   MSG_BUG(sjoin('Wrong Er%ID:',itoa(Er%ID)))
 END SELECT

 ! For chi, \espilon or \epsilon^{-1}, define the approximation.
 rfapprox='None'
 if (Er%ID>=2.or.Er%ID<=4) then
   if (Er%ikxc==0) then
     rfapprox='RPA'
   else if (Er%ikxc>0) then
     rfapprox='Static TDDFT'
   else
     rfapprox='TDDFT'
   end if
 end if

 ! === If TDDFT and \epsilon^{-1}, define the type ===
 rftest='None'
! if (Er%ID==0) then
!  if (Er%test_type==0) then
!   rftest='TEST-PARTICLE'
!  else if (Er%test_type==1) then
!   rftest='TEST-ELECTRON'
!  else
!   write(msg,'(4a,i3)')ch10,&
!&   ' em1results_print : BUG - ',ch10,&
!&   ' Wrong value of Er%test_type = ',Er%test_type
!   MSG_ERROR(msg)
!  end if
! end if

 ! === Define time-ordering ===
 rforder='Undefined'
 if (Er%Tordering==1) then
   rforder='Time-Ordered'
 else if (Er%Tordering==2) then
   rforder='Advanced'
 else if (Er%Tordering==3) then
   rforder='Retarded'
 else
   MSG_BUG(sjoin('Wrong er%tordering= ',itoa(Er%Tordering)))
 end if

 kxcname='None'
 if (Er%ikxc/=0) then
   !TODO Add function to retrieve kxc name
   MSG_ERROR('Add function to retrieve kxc name')
   kxcname='XXXXX'
 end if

 write(msg,'(6a,5(3a))')ch10,&
&  ' ==== Info on the Response Function ==== ',ch10,&
&  '  Associated File ................  ',TRIM(Er%fname),ch10,&
&  '  Response Function Type .......... ',TRIM(rfname),ch10,&
&  '  Type of Approximation ........... ',TRIM(rfapprox),ch10,&
&  '  XC kernel used .................. ',TRIM(kxcname),ch10,&
&  '  Type of probing particle ........ ',TRIM(rftest),ch10,&
&  '  Time-Ordering ................... ',TRIM(rforder),ch10
 call wrtout(unt,msg,mode)
 write(msg,'(a,2i4,a,3(a,i4,a),a,3i4,2a,i4,a)')&
&  '  Number of components ............ ',Er%nI,Er%nJ,ch10,&
&  '  Number of q-points in the IBZ ... ',Er%nqibz,ch10,&
&  '  Number of q-points for q-->0 .... ',Er%nqlwl,ch10,&
&  '  Number of G-vectors ............. ',Er%npwe,ch10,&
&  '  Number of frequencies ........... ',Er%nomega,Er%nomega_r,Er%nomega_i,ch10,&
&  '  Value of mqmem .................. ',Er%mqmem,ch10
 call wrtout(unt,msg,mode)

 if (Er%nqlwl/=0) then
   write(msg,'(a,i3)')' q-points for long wavelength limit: ',Er%nqlwl
   call wrtout(unt,msg,mode)
   do iqlwl=1,Er%nqlwl
     write(msg,'(1x,i5,a,3es16.8)')iqlwl,') ',Er%qlwl(:,iqlwl)
     call wrtout(unt,msg,mode)
   end do
 end if

 if (my_prtvol>0) then ! Print out head and wings in the long-wavelength limit.
   ! TODO add additional stuff.
   write(msg,'(a,i4)')' Calculated Frequencies: ',Er%nomega
   call wrtout(unt,msg,mode)
   do iw=1,Er%nomega
     write(msg,'(i4,es14.6)')iw,Er%omega(iw)*Ha_eV
     call wrtout(unt,msg,mode)
   end do

   write(msg,'(a,i4)')' Calculated q-points: ',Er%nqibz
   call wrtout(unt,msg,mode)
   do iqibz=1,Er%nqibz
     write(msg,'(1x,i4,a,3es16.8)')iqibz,') ',Er%qibz(:,iqibz)
     call wrtout(unt,msg,mode)
   end do
 end if ! my_prtvol>0

end subroutine em1results_print
!!***

!----------------------------------------------------------------------

!!****f* m_screening/Epsm1_symmetrizer
!! NAME
!!  Epsm1_symmetrizer
!!
!! FUNCTION
!!  Symmetrize the inverse dielectric matrix, namely calculate epsilon^{-1} at a generic
!!  q-point in the BZ starting from the knowledge of the matrix at a q-point in the IBZ.
!!  The procedure is quite generic and can be used for every two-point function which has
!!  the same symmetry as the crystal.
!!
!! INPUTS
!!  nomega=Number of frequencies required. All frequencies from 1 up to nomega are symmetrized.
!!  npwc=Number of G vectors in symmetrized matrix, has to be smaller than Er%npwe.
!!  remove_exchange=If .TRUE., return e^{-1}-1 namely remove the exchange part.
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix.
!!  Gsph<gsphere_t>=data related to the G-sphere
!!    %grottb
!!    %phmSGt
!!  Qmesh<kmesh_t>=Structure defining the q-mesh used for Er.
!!    %nbz=Number of q-points in the BZ
!!    %tab(nbz)=Index of the symmetric q-point in the IBZ, for each point in the BZ
!!    %tabo(nbz)=The operation that rotates q_ibz onto \pm q_bz (depending on tabi)
!!    %tabi(nbz)=-1 if time-reversal has to be considered, 1 otherwise
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required.
!!
!! OUTPUT
!!  epsm1_qbz(npwc,npwc,nomega)=The inverse dielectric matrix at the q-point defined by iq_bz.
!!   Exchange part can be subtracted out.
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required
!!  to reconstruct the BZ.
!!
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!    If q_bz=Sq_ibz+G0:
!!
!!    $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then :
!!    $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!!
!! TODO
!!  Symmetrization can be skipped if iq_bz correspond to a point in the IBZ
!!
!! PARENTS
!!      calc_sigc_me,cohsex_me
!!
!! CHILDREN
!!
!! SOURCE

subroutine Epsm1_symmetrizer(iq_bz,nomega,npwc,Er,Gsph,Qmesh,remove_exchange,epsm1_qbz)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,nomega,npwc
 logical,intent(in) :: remove_exchange
 type(Epsilonm1_results),intent(in) :: Er
 type(gsphere_t),target,intent(in) :: Gsph
 type(kmesh_t),intent(in) :: Qmesh
!arrays
 complex(gwpc),intent(out) :: epsm1_qbz(npwc,npwc,nomega)

!Local variables-------------------------------
!scalars
 integer :: iw,ii,jj,iq_ibz,itim_q,isym_q,iq_loc,sg1,sg2
 complex(gwpc) :: phmsg1t,phmsg2t_star
!arrays
 real(dp) :: qbz(3)

! *********************************************************************

 ABI_CHECK(Er%nomega>=nomega,'Too many frequencies required')
 ABI_CHECK(Er%npwe  >=npwc , 'Too many G-vectors required')

 ! * Get iq_ibz, and symmetries from iq_ibz.
 call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)

 ! If out-of-memory, only Er%espm1(:,:,:,1) has been allocated and filled.
 iq_loc=iq_ibz; if (Er%mqmem==0) iq_loc=1

 ! MG: grottb is a 1-1 mapping, hence we can collapse the loops (false sharing is not an issue here).
 !grottb => Gsph%rottb (1:npwc,itim_q,isym_q)
 !phmSgt => Gsph%phmSGt(1:npwc,isym_q)

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(sg2,sg1,phmsg1t,phmsg2t_star)
 do iw=1,nomega
   do jj=1,npwc
     sg2 = Gsph%rottb(jj,itim_q,isym_q)
     phmsg2t_star = CONJG(Gsph%phmSGt(jj,isym_q))
     do ii=1,npwc
       sg1 = Gsph%rottb(ii,itim_q,isym_q)
       phmsg1t = Gsph%phmSGt(ii,isym_q)
       epsm1_qbz(sg1,sg2,iw) = Er%epsm1(ii,jj,iw,iq_loc) * phmsg1t * phmsg2t_star
       !epsm1_qbz(sg1,sg2,iw) = Er%epsm1(ii,jj,iw,iq_loc) * phmSgt(ii) * CONJG(phmSgt(jj))
     end do
   end do
 end do
 !
 ! === Account for time-reversal ===
 !epsm1_qbz(:,:,iw)=TRANSPOSE(epsm1_qbz(:,:,iw))
 if (itim_q==2) then
!$OMP PARALLEL DO IF (nomega > 1)
   do iw=1,nomega
     call sqmat_itranspose(npwc,epsm1_qbz(:,:,iw))
   end do
 end if

 if (remove_exchange) then
   ! === Subtract the exchange contribution ===
   ! If it's a pole screening, the exchange contribution is already removed
!$OMP PARALLEL DO IF (nomega > 1)
   do iw=1,nomega
     do ii=1,npwc
       epsm1_qbz(ii,ii,iw)=epsm1_qbz(ii,ii,iw)-CMPLX(1.0_gwp,0.0_gwp)
     end do
   end do
 endif

end subroutine Epsm1_symmetrizer
!!***

!----------------------------------------------------------------------

!!****f* m_screening/Epsm1_symmetrizer_inplace
!! NAME
!!  Epsm1_symmetrizer_inplace
!!
!! FUNCTION
!!  Same function as Epsm1_symmetrizer, ecept now the array Ep%epsm1 is modified inplace
!!  thorugh an auxiliary work array of dimension (npwc,npwc)
!!
!! INPUTS
!!  nomega=Number of frequencies required. All frequencies from 1 up to nomega are symmetrized.
!!  npwc=Number of G vectors in symmetrized matrix, has to be smaller than Er%npwe.
!!  remove_exchange=If .TRUE., return e^{-1}-1 namely remove the exchange part.
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix.
!!  Gsph<gsphere_t>=data related to the G-sphere
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix.
!!  Gsph<gsphere_t>=data related to the G-sphere
!!    %grottb
!!    %phmSGt
!!  Qmesh<kmesh_t>=Structure defining the q-mesh used for Er.
!!    %nbz=Number of q-points in the BZ
!!    %tab(nbz)=Index of the symmetric q-point in the IBZ, for each point in the BZ
!!    %tabo(nbz)=The operation that rotates q_ibz onto \pm q_bz (depending on tabi)
!!    %tabi(nbz)=-1 if time-reversal has to be considered, 1 otherwise
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required.
!!
!! OUTPUT
!!  Er%epsm1(npwc,npwc,nomega,iq_loc) symmetrised
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required
!!  to reconstruct the BZ.
!!
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!    If q_bz=Sq_ibz+G0:
!!
!!    $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then :
!!    $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!!
!! TODO
!!  Symmetrization can be skipped if iq_bz correspond to a point in the IBZ
!!
!! PARENTS
!!      calc_sigc_me
!!
!! CHILDREN
!!
!! SOURCE

subroutine Epsm1_symmetrizer_inplace(iq_bz,nomega,npwc,Er,Gsph,Qmesh,remove_exchange)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,nomega,npwc
 logical,intent(in) :: remove_exchange
 type(Epsilonm1_results) :: Er
 type(gsphere_t),target,intent(in) :: Gsph
 type(kmesh_t),intent(in) :: Qmesh

!Local variables-------------------------------
!scalars
 integer :: iw,ii,jj,iq_ibz,itim_q,isym_q,iq_loc,sg1,sg2
!arrays
 real(dp) :: qbz(3)
 complex(gwpc) :: phmsg1t,phmsg2t_star
 complex(gwpc),allocatable :: work(:,:)

! *********************************************************************

 ABI_CHECK(Er%nomega>=nomega,'Too many frequencies required')
 ABI_CHECK(Er%npwe  >=npwc , 'Too many G-vectors required')

 ABI_MALLOC(work,(npwc,npwc))

 ! * Get iq_ibz, and symmetries from iq_ibz.
 call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)

 ! If out-of-memory, only Er%espm1(:,:,:,1) has been allocated and filled.
 iq_loc=iq_ibz; if (Er%mqmem==0) iq_loc=1

!$OMP PARALLEL DO PRIVATE(sg2,sg1,phmsg1t,phmsg2t_star) IF (nomega > 1)
 do iw=1,nomega
   do jj=1,npwc
     sg2 = Gsph%rottb(jj,itim_q,isym_q)
     phmsg2t_star = CONJG(Gsph%phmSGt(jj,isym_q))
     do ii=1,npwc
       sg1 = Gsph%rottb(ii,itim_q,isym_q)
       phmsg1t = Gsph%phmSGt(ii,isym_q)
       work(sg1,sg2) = Er%epsm1(ii,jj,iw,iq_loc) * phmsg1t * phmsg2t_star
       !work(grottb(ii),grottb(jj))=Er%epsm1(ii,jj,iw,iq_loc)*phmSgt(ii)*CONJG(phmSgt(jj))
     end do
   end do
   Er%epsm1(:,:,iw,iq_loc) = work(:,:)
 end do
 !
 ! === Account for time-reversal ===
 !Er%epsm1(:,:,iw,iq_loc)=TRANSPOSE(Er%epsm1(:,:,iw,iq_loc))
 if (itim_q==2) then
!$OMP PARALLEL DO IF (nomega > 1)
   do iw=1,nomega
     call sqmat_itranspose(npwc,Er%epsm1(:,:,iw,iq_loc))
   end do
 end if

 ! === Subtract the exchange contribution ===
 if (remove_exchange) then
!$OMP PARALLEL DO IF (nomega > 1)
   do iw=1,nomega
     do ii=1,npwc
       Er%epsm1(ii,ii,iw,iq_loc)=Er%epsm1(ii,ii,iw,iq_loc)-1.0_gwp
     end do
   end do
 endif

 ABI_FREE(work)

end subroutine Epsm1_symmetrizer_inplace
!!***

!----------------------------------------------------------------------

!!****f* m_screening/init_Er_from_file
!! NAME
!!  init_Er_from_file
!!
!! FUNCTION
!!  Initialize basic dimensions and the important (small) arrays in an Epsilonm1_results data type
!!  starting from a file containing either epsilon^{-1} (_SCR) or chi0 (_SUSC).
!!
!! INPUTS
!!  fname=The name of the external file used to read the matrix.
!!  mqmem=0 for out-of-core solution, /=0 if entire matrix has to be stored in memory.
!!  npwe_asked=Number of G-vector to be used in the calculation, if <=0 use Max allowed number.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Er<Epsilonm1_results>=The structure initialized with basic dimensions and arrays.
!!
!! PARENTS
!!      m_screening,mrgscr,setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_Er_from_file(Er,fname,mqmem,npwe_asked,comm)

!Arguments ------------------------------------
 integer,intent(in) :: mqmem,npwe_asked,comm
 character(len=*),intent(in) :: fname
 type(Epsilonm1_results),intent(inout) :: Er

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iw,fform,my_rank,unclassified
 real(dp) :: re, im, tol
 character(len=500) :: msg

! *********************************************************************

 DBG_ENTER("COLL")

 !@Epsilonm1_results
 my_rank = xmpi_comm_rank(comm)

 ! Read header from file.
 call wrtout(std_out,sjoin('init_Er_from_file- testing file: ',fname),'COLL')
 call hscr_from_file(Er%hscr, fname, fform, comm)

 ! Master echoes the header.
 if (my_rank==master) call hscr_print(er%hscr)

 ! Generic Info
 Er%ID         =0       ! Not yet initialized as epsm1 is calculated in mkdump_Er.F90
 Er%fname      =fname
 Er%fform      =fform
 Er%Tordering=Er%Hscr%Tordering

!TODO these quantitities should be checked and initiliazed in mkdump_Er
!BEGIN HARCODED
 Er%nI       = 1
 Er%nJ       = 1
 Er%ikxc     = 0
 Er%test_type=-1

 Er%Hscr%headform=HSCR_LATEST_HEADFORM   ! XG20090912
!END HARDCODED

 Er%nqibz=Er%Hscr%nqibz
 Er%mqmem=mqmem ; if (mqmem/=0) Er%mqmem=Er%nqibz
 ABI_MALLOC(Er%qibz,(3,Er%nqibz))
 Er%qibz(:,:)=Er%Hscr%qibz(:,:)

 Er%nqlwl=Er%Hscr%nqlwl
 ABI_MALLOC(Er%qlwl,(3,Er%nqlwl))
 Er%qlwl(:,:)=Er%Hscr%qlwl(:,:)

 Er%nomega=Er%Hscr%nomega
 ABI_MALLOC(Er%omega,(Er%nomega))
 Er%omega(:)=Er%Hscr%omega(:)

 ! Count number of real, imaginary, and complex frequencies.
 Er%nomega_r = 1
 Er%nomega_i = 0
 if (Er%nomega == 2) then
   Er%nomega_i = 1
 else
   unclassified = 0
   tol = tol6*Ha_eV
   do iw = 2, Er%nomega
     re =  REAL(Er%omega(iw))
     im = AIMAG(Er%omega(iw))
     if (re > tol .and. im < tol) then
       Er%nomega_r = iw ! Real freqs are packed in the first locations.
     else if (re < tol .and. im > tol) then
       Er%nomega_i = Er%nomega_i + 1
     else
       unclassified = unclassified + 1
     end if
   end do
   if (unclassified > 0) then
     write(msg,'(3a,i6)')&
&      'Some complex frequencies are too small to qualify as real or imaginary.',ch10,&
&      'Number of unidentified frequencies = ', unclassified
     MSG_WARNING(msg)
   end if
 end if

 ! Get G-vectors.
 Er%npwe=Er%Hscr%npwe
 if (npwe_asked>0) then
   if (npwe_asked>Er%Hscr%npwe) then
     write(msg,'(a,i8,2a,i8)')&
&     'Number of G-vectors saved on file is less than the value required = ',npwe_asked,ch10,&
&     'Calculation will proceed with Max available npwe = ',Er%Hscr%npwe
     MSG_WARNING(msg)
   else  ! Redefine the no. of G"s for W.
     Er%npwe=npwe_asked
   end if
 end if

 ! pointer to Er%Hscr%gvec ?
 ABI_MALLOC(Er%gvec,(3,Er%npwe))
 Er%gvec=Er%Hscr%gvec(:,1:Er%npwe)

 DBG_EXIT("COLL")

end subroutine init_Er_from_file
!!***

!----------------------------------------------------------------------

!!****f* m_screening/mkdump_Er
!! NAME
!!  mkdump_Er
!!
!! FUNCTION
!!  Dump the content of an Epsilonm1_results data type on file.
!!
!! INPUTS
!!  id_required=Identifier of the matrix to be calculated
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!  ngfft(18)=Info on the FFT mesh.
!!  nfftot=Total number of point on the FFT mesh.
!!  gvec(3,npwe)=Reduced coordinates of plane waves for the response functions
!!  npwe=Number of plane waves.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!
!! PARENTS
!!      mrgscr,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkdump_Er(Er,Vcp,npwe,gvec,nkxc,kxcg,id_required,approx_type,&
&                    ikxc_required,option_test,fname_dump,iomode,&
&                    nfftot,ngfft,comm,fxc_ADA)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: id_required,approx_type,option_test,ikxc_required,nkxc
 integer,intent(in) :: iomode,nfftot,npwe,comm
 type(Epsilonm1_results),intent(inout) :: Er
 type(vcoul_t),intent(in) :: Vcp
 character(len=*),intent(in) :: fname_dump
!arrays
 integer,intent(in) :: ngfft(18),gvec(3,npwe)
 complex(gwpc),intent(in) :: kxcg(nfftot,nkxc)
 complex(gwpc),intent(in), optional :: fxc_ADA(npwe*Er%nI,npwe*Er%nJ,Er%nqibz)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: dim_wing,iqibz,is_qeq0,mqmem_,npwe_asked
 integer :: unt_dump,fform,rdwr,ierr
 integer :: my_rank,comm_self
 real(dp) :: ucvol
 character(len=500) :: msg
 character(len=fnlen) :: ofname
 character(len=nctk_slen) :: in_varname,out_varname
 type(hscr_t) :: Hscr_cp
 type(spectra_t) :: spectra
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 complex(gwpc),allocatable :: epsm1(:,:,:)
 complex(dpc),allocatable :: dummy_lwing(:,:,:),dummy_uwing(:,:,:),dummy_head(:,:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(id_required==4,'Value of id_required not coded')
 ABI_CHECK(npwe==Er%npwe,"mismatch in npwe")

 my_rank = xmpi_comm_rank(comm); comm_self = xmpi_comm_self
 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)

 ! if (Er%ID/=0) call reset_Epsilonm1(Er)
 Er%ID=id_required

 ofname = fname_dump
 in_varname = ncname_from_id(er%hscr%id)
 out_varname = ncname_from_id(id_required)

 write(std_out,*)'Er%ID: ',Er%ID,', Er%Hscr%ID: ',Er%Hscr%ID

 if (Er%ID==Er%Hscr%ID) then
   ! === The two-point function we are asking for is already stored on file ===
   ! * According to mqmem either read and store the entire matrix in memory or do nothing.

   if (Er%mqmem>0) then
     ! In-core solution.
     write(msg,'(a,f12.1,a)')' Memory needed for Er%epsm1 = ',two*gwpc*npwe**2*Er%nomega*Er%nqibz*b2Mb,' [Mb] <<< MEM'
     call wrtout(std_out,msg)
     ABI_MALLOC_OR_DIE(Er%epsm1,(npwe,npwe,Er%nomega,Er%nqibz), ierr)

     if (iomode == IO_MODE_MPI) then
       !call wrtout(std_out, "read_screening with MPI_IO")
       MSG_WARNING("SUSC files is buggy. Using Fortran IO")
       call read_screening(in_varname,Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,IO_MODE_FORTRAN,comm)
     else
       call read_screening(in_varname,Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,iomode,comm)
     end if
   else
     ! Out-of-core solution ===
     MSG_COMMENT("mqmem==0 => allocating a single q-slice of (W|chi0) (slower but less memory).")
     continue
   end if

   return

 else
   ! === The matrix stored on file do not correspond to the quantity required ===
   ! * Presently only the transformation chi0 => e^-1 is coded
   ! * According to Er%mqmem either calculate e^-1 dumping the result to a file
   !   for a subsequent use or calculate e^-1 keeping everything in memory.

   if (Er%mqmem==0) then
     ! === Open file and write the header for the SCR file ===
     ! * For the moment only master works.

     if (my_rank==master) then
       if (iomode == IO_MODE_ETSF) then
#ifdef HAVE_NETCDF
          ofname = nctk_ncify(ofname)
          NCF_CHECK(nctk_open_create(unt_dump, ofname, xmpi_comm_self))
#endif
       else
         if (open_file(ofname,msg,newunit=unt_dump,form="unformatted",status="unknown",action="write") /= 0) then
           MSG_ERROR(msg)
         end if
       end if
       call wrtout(std_out,sjoin('mkdump_Er: calculating and writing epsilon^-1 matrix on file: ',ofname),'COLL')

       ! Update the entries in the header that have been modified.
       ! TODO, write function to return title, just for info
       call hscr_copy(Er%Hscr,Hscr_cp)
       Hscr_cp%ID = id_required
       Hscr_cp%ikxc = ikxc_required
       Hscr_cp%test_type = option_test
       Hscr_cp%titles(1)  = 'SCR file: epsilon^-1'
       Hscr_cp%titles(2)  = 'TESTPARTICLE'
       ! Treat the case in which a smaller matrix is used.
       Hscr_cp%npwe = npwe

       rdwr=2; fform=Hscr_cp%fform
       call hscr_io(hscr_cp,fform,rdwr,unt_dump,comm_self,master,iomode)
       call hscr_free(Hscr_cp)

       ABI_MALLOC_OR_DIE(epsm1,(npwe,npwe,Er%nomega), ierr)

       do iqibz=1,Er%nqibz
         is_qeq0=0
         if (normv(Er%qibz(:,iqibz),gmet,'G')<GW_TOLQ0) is_qeq0=1
         ! FIXME there's a problem with SUSC files and MPI-IO
         !if (iomode == IO_MODE_MPI) then
         !  MSG_WARNING("SUSC files is buggy. Using Fortran IO")
         !  call read_screening(in_varname,Er%fname,npwe,1,Er%nomega,epsm1,IO_MODE_FORTRAN,comm_self,iqiA=iqibz)
         !else
         call read_screening(in_varname,Er%fname,npwe,1,Er%nomega,epsm1,iomode,comm_self,iqiA=iqibz)
         !end if

         dim_wing=0; if (is_qeq0==1) dim_wing=3
         ABI_MALLOC(dummy_lwing,(npwe*Er%nI,Er%nomega,dim_wing))
         ABI_MALLOC(dummy_uwing,(npwe*Er%nJ,Er%nomega,dim_wing))
         ABI_MALLOC(dummy_head,(dim_wing,dim_wing,Er%nomega))

         if (approx_type<2 .or. approx_type>3) then ! bootstrap
           MSG_WARNING('Entering out-of core RPA or Kxc branch')
           call make_epsm1_driver(iqibz,dim_wing,npwe,Er%nI,Er%nJ,Er%nomega,Er%omega,&
&                    approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,dummy_head,&
&                    dummy_lwing,dummy_uwing,epsm1,spectra,comm_self)
         else
           MSG_WARNING('Entering out-of core fxc_ADA branch')
           call make_epsm1_driver(iqibz,dim_wing,npwe,Er%nI,Er%nJ,Er%nomega,Er%omega,&
&                    approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,dummy_head,&
&                    dummy_lwing,dummy_uwing,epsm1,spectra,comm_self,fxc_ADA(:,:,iqibz))
         end if

         ABI_FREE(dummy_head)
         ABI_FREE(dummy_uwing)
         ABI_FREE(dummy_lwing)

         if (is_qeq0==1) then
           call spectra%repr(msg)
           call wrtout(std_out,msg,'COLL')
           call wrtout(ab_out,msg,'COLL')
         end if
         call spectra%free()

         call write_screening(out_varname,unt_dump,iomode,npwe,Er%nomega,iqibz,epsm1)
       end do

       if (iomode == IO_MODE_ETSF) then
#ifdef HAVE_NETCDF
         NCF_CHECK(nf90_close(unt_dump))
#endif
       else
         close(unt_dump)
       endif

       ABI_FREE(epsm1)
     end if !master

     ! Master broadcasts ofname.
     ! NOTE: A synchronization is required here, else the other procs start to read the
     ! SCR file before it is written by the master. xmpi_bcast will synch the procs.
     call xmpi_bcast(ofname,  master, comm, ierr)

     ! Now Er% "belongs" to the file "ofname", thus
     ! each proc has to destroy and re-initialize the object.
     call em1results_free(Er)

     mqmem_=Er%mqmem; npwe_asked=npwe
     call init_Er_from_file(Er,ofname,mqmem_,npwe_asked,comm)

     !Now Er% has been reinitialized and ready-to-use.
     Er%id = id_required
     call em1results_print(Er)
   else
     ! ========================
     ! === In-core solution ===
     ! ========================
     ABI_MALLOC_OR_DIE(Er%epsm1,(npwe,npwe,Er%nomega,Er%nqibz), ierr)

     ! FIXME there's a problem with SUSC files and MPI-IO
     !if (iomode == IO_MODE_MPI) then
     !  !call wrtout(std_out, "read_screening with MPI_IO")
     !  MSG_WARNING("SUSC files is buggy. Using Fortran IO")
     !  call read_screening(in_varname,Er%fname,npwe,Er%nqibz,Er%nomega,Er%epsm1,IO_MODE_FORTRAN,comm)
     !else
     call read_screening(in_varname,Er%fname,npwe,Er%nqibz,Er%nomega,Er%epsm1,iomode,comm)
     !end if

     do iqibz=1,Er%nqibz
       is_qeq0=0; if (normv(Er%qibz(:,iqibz),gmet,'G')<GW_TOLQ0) is_qeq0=1

       dim_wing=0; if (is_qeq0==1) dim_wing=3 ! FIXME
       ABI_MALLOC(dummy_lwing,(npwe*Er%nI,Er%nomega,dim_wing))
       ABI_MALLOC(dummy_uwing,(npwe*Er%nJ,Er%nomega,dim_wing))
       ABI_MALLOC(dummy_head,(dim_wing,dim_wing,Er%nomega))

       if (approx_type<2 .or. approx_type>3) then
         MSG_WARNING('Entering in-core RPA and Kxc branch')
         call make_epsm1_driver(iqibz,dim_wing,npwe,Er%nI,Er%nJ,Er%nomega,Er%omega,&
&                  approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,dummy_head,&
&                  dummy_lwing,dummy_uwing,Er%epsm1(:,:,:,iqibz),spectra,comm)
       else
         MSG_WARNING('Entering in-core fxc_ADA branch')
         call make_epsm1_driver(iqibz,dim_wing,npwe,Er%nI,Er%nJ,Er%nomega,Er%omega,&
&                  approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,dummy_head,&
&                  dummy_lwing,dummy_uwing,Er%epsm1(:,:,:,iqibz),spectra,comm,fxc_ADA=fxc_ADA(:,:,iqibz))
       end if

       ABI_FREE(dummy_lwing)
       ABI_FREE(dummy_uwing)
       ABI_FREE(dummy_head)

       if (is_qeq0==1) then
         call spectra%repr(msg)
         call wrtout(std_out,msg,'COLL')
         call wrtout(ab_out,msg,'COLL')
       end if

       call spectra%free()
     end do

     Er%id = id_required
     call em1results_print(Er)
   end if
 end if

 DBG_EXIT("COLL")

end subroutine mkdump_Er
!!***

!----------------------------------------------------------------------

!!****f* m_screening/get_epsm1
!! NAME
!!  get_epsm1
!!
!! FUNCTION
!!  Work in progress but the main is idea is as follows:
!!
!!  Return the symmetrized inverse dielectric matrix.
!!  This method implements both in-core and the out-of-core solution
!!  In the later, epsilon^-1 or chi0 are read from file.
!!  It is possible to specify options to retrieve (RPA |TDDDT, [TESTCHARGE|TESTPARTICLE]).
!!  All dimensions are already initialized in the Er% object, this method
!!  should act as a wrapper around rdscr and make_epsm1_driver. A better
!!  implementation will be done in the following once the coding of file handlers is completed.
!!
!! INPUTS
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!  iqibzA[optional]=Index of the q-point to be read from file (only for out-of-memory solutions)
!!  iomode=option definig the file format.
!!  option_test
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Er%epsm1
!!
!! TODO
!!  Remove this routine. Now everything should be done with mkdump_Er
!!
!! PARENTS
!!      calc_sigc_me,cohsex_me
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_epsm1(Er,Vcp,approx_type,option_test,iomode,comm,iqibzA)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,option_test,approx_type,comm
 integer,optional,intent(in) :: iqibzA
 type(vcoul_t),intent(in) :: Vcp
 type(Epsilonm1_results),intent(inout) :: Er

!Local variables-------------------------------
!scalars
 integer :: my_approx_type,my_option_test,ng,ierr
! *********************************************************************

 DBG_ENTER("COLL")

 my_approx_type = approx_type; my_option_test = option_test

 ! Vcp not yet used.
 ng = Vcp%ng

 select case (Er%mqmem)
 case (0)
   !  Out-of-core solution
   if (allocated(Er%epsm1))  then
     ABI_FREE(Er%epsm1)
   end if
   ABI_MALLOC_OR_DIE(Er%epsm1,(Er%npwe,Er%npwe,Er%nomega,1), ierr)

   ! FIXME there's a problem with SUSC files and MPI-IO
   !if (iomode == IO_MODE_MPI) then
   !  !write(std_out,*)"read_screening with iomode",iomode,"file: ",trim(er%fname)
   !  MSG_WARNING("SUSC files is buggy. Using Fortran IO")
   !  call read_screening(em1_ncname,Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,IO_MODE_FORTRAN,comm,iqiA=iqibzA)
   !else
   call read_screening(em1_ncname,Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,iomode,comm,iqiA=iqibzA)
   !end if

   if (Er%id == 4) then
     ! If q-slice of epsilon^-1 has been read then return
     !call em1results_print(Er)
     return
   else
     MSG_ERROR(sjoin('Wrong Er%ID', itoa(er%id)))
   end if

 case default
   ! In-core solution.
   MSG_ERROR("you should not be here")
 end select

 DBG_EXIT("COLL")

end subroutine get_epsm1
!!***

!----------------------------------------------------------------------

!!****f* m_screening/decompose_epsm1
!! NAME
!! decompose_epsm1
!!
!! FUNCTION
!! Decompose the complex symmetrized dielectric
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine decompose_epsm1(Er,iqibz,eigs)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz
 type(Epsilonm1_results),intent(in) :: Er
!arrays
 complex(dpc),intent(out) :: eigs(Er%npwe,Er%nomega)

!Local variables-------------------------------
!scalars
 integer :: info,lwork,iw,negw,ig1,ig2,idx,sdim,npwe,ierr
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: ww(:),rwork(:)
 complex(dpc),allocatable :: work(:),Adpp(:),eigvec(:,:),Afull(:,:),vs(:,:),wwc(:)
 logical,allocatable :: bwork(:)
 logical :: sortcplx !BUG in abilint

! *********************************************************************

 ABI_CHECK(Er%mqmem/=0,'mqmem==0 not implemented')

 npwe = Er%npwe

 do iw=1,Er%nomega

   if (ABS(REAL(Er%omega(iw)))>0.00001) then
     ! Eigenvalues for a generic complex matrix
     lwork=4*2*npwe
     ABI_MALLOC(wwc,(npwe))
     ABI_MALLOC(work,(lwork))
     ABI_MALLOC(rwork,(npwe))
     ABI_MALLOC(bwork,(npwe))
     ABI_MALLOC(vs,(npwe,npwe))
     ABI_MALLOC(Afull,(npwe,npwe))

     Afull=Er%epsm1(:,:,iw,iqibz)

     !for the moment no sort, maybe here I should sort using the real part?
     call ZGEES('V','N',sortcplx,npwe,Afull,npwe,sdim,wwc,vs,npwe,work,lwork,rwork,bwork,info)
     if (info/=0) then
       MSG_ERROR(sjoin("ZGEES returned info:",itoa(info)))
     end if

     eigs(:,iw)=wwc(:)

     ABI_FREE(wwc)
     ABI_FREE(work)
     ABI_FREE(rwork)
     ABI_FREE(bwork)
     ABI_FREE(vs)
     ABI_FREE(Afull)

   else
     ! Hermitian version.
     lwork=2*npwe-1
     ABI_MALLOC(ww,(npwe))
     ABI_MALLOC(work,(lwork))
     ABI_MALLOC(rwork,(3*npwe-2))
     ABI_MALLOC(eigvec,(npwe,npwe))
     ABI_MALLOC_OR_DIE(Adpp,(npwe*(npwe+1)/2), ierr)

     idx=0 ! Pack the matrix
     do ig2=1,npwe
       do ig1=1,ig2
         idx=idx+1
         Adpp(idx)=Er%epsm1(ig1,ig2,iw,iqibz)
       end do
     end do

     ! For the moment we require also the eigenvectors.
     call ZHPEV('V','U',npwe,Adpp,ww,eigvec,npwe,work,rwork,info)
     if (info/=0) then
       MSG_ERROR(sjoin('ZHPEV returned info=', itoa(info)))
     end if

     negw=(COUNT((REAL(ww)<tol6)))
     if (negw/=0) then
       write(msg,'(a,i5,a,i3,a,f8.4)')&
        'Found negative eigenvalues. No. ',negw,' at iqibz= ',iqibz,' minval= ',MINVAL(REAL(ww))
       MSG_WARNING(msg)
     end if

     eigs(:,iw)=ww(:)

     ABI_FREE(ww)
     ABI_FREE(work)
     ABI_FREE(rwork)
     ABI_FREE(eigvec)
     ABI_FREE(Adpp)
   end if
 end do !iw

! contains
! function sortcplx(carg) result(res)
!  implicit none
!  complex(dpc),intent(in) :: carg
!  logical :: res
!  res=.TRUE.
! end function sortcplx

end subroutine decompose_epsm1
!!***

!----------------------------------------------------------------------

!!****f* m_screening/make_epsm1_driver
!! NAME
!! make_epsm1_driver
!!
!! FUNCTION
!!  Driver routine to calculate the inverse symmetrical dielectric matrix starting
!!  from the irreducible polarizability. The routine considers a single q-point, and
!!  performs the following tasks:
!!
!!  1) Calculate $\tilde\epsilon^{-1}$ using different approximations:
!!      * RPA
!!      * ALDA within TDDFT
!!
!!  2) Use a special treatment of non-Analytic behavior of heads and wings in reciprocal space
!!     calculating these quantities for different small q-directions specified by the user
!!     (Not yet operative)
!!
!!  3) Output the electron energy loss function and the macroscopic dielectric function with and
!!     without local field effects (only if non-zero real frequencies are available)
!!
!! INPUTS
!!  iqibz=index of the q-point in the array Vcp%qibz where epsilon^-1 has to be calculated
!!  dim_wing=Dimension of the wings (0 or 3 if q-->0)
!!  npwe=Number of G-vectors in chi0.
!!  nI,nJ=Number of rows/columns in chi0_ij (1,1 if collinear case)
!!  nomega=Number of frequencies.
!!  omega(nomega)=Frequencines in Hartree
!!  approx_type=Integer flag defining the type of approximation
!!   == 0 for RPA   ==
!!   == 1 for TDDFT ==
!!  option_test=Only for TDDFT:
!!   == 0 for TESTPARTICLE ==
!!   == 1 for TESTELECTRON ==
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!  nfftot=Total number of points in the FFT mesh.
!!  ngfft(18)=Info on the FFT mesh.
!!  nkxc=Dimension of the kernel in reciprocal space. 0 if kernel is not needed
!!  kxcg(nfftot,nkxc)=TDDFT kernel in reciprocal space on the FFT mesh. Used only if approx_type==1
!!  gvec(3,npwe)=G-vectors
!!  comm=MPI communicator.
!!  chi0_lwing(npwe*nI,nomega,dim_wing)=Lower wings of chi0 (only for q-->0)
!!  chi0_uwing(npwe*nJ,nomega,dim_wing)=Upper wings of chi0 (only for q-->0)
!!  chi0_head(dim_wing,dim_wing,nomega)=Head of of chi0 (only for q-->0)
!!
!! OUTPUT
!!  spectra<spectra_t>Object containing e_macro(w) and EELS(w)
!!
!! SIDE EFFECTS
!!  chi0(npwe*nI,npwe*nJ,nomega): in input the irreducible polarizability, in output
!!   the symmetrized inverse dielectric matrix.
!!
!! PARENTS
!!      m_screening,screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_epsm1_driver(iqibz,dim_wing,npwe,nI,nJ,nomega,omega,&
  approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,chi0_head,&
  chi0_lwing,chi0_uwing,chi0,spectra,comm,&
  fxc_ADA) ! optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,nI,nJ,npwe,nomega,dim_wing,approx_type,option_test,nkxc,nfftot,comm
 type(vcoul_t),target,intent(in) :: Vcp
 type(spectra_t),intent(out) :: Spectra
!arrays
 integer,intent(in) :: ngfft(18),gvec(3,npwe)
 complex(gwpc),intent(in) :: kxcg(nfftot,nkxc)
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(inout) :: chi0_lwing(:,:,:)   !(npwe*nI,nomega,dim_wing)
 complex(dpc),intent(inout) :: chi0_uwing(:,:,:)   !(npwe*nJ,nomega,dim_wing)
 complex(dpc),intent(inout) :: chi0_head(:,:,:)   !(dim_wing,dim_wing,nomega)
 complex(gwpc),intent(inout) :: chi0(npwe*nI,npwe*nJ,nomega)
 complex(gwpc),intent(in),optional :: fxc_ADA(npwe*nI,npwe*nJ)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: i1,i2,ig1,ig2,io,ierr,irank,my_nqlwl !iqlwl
 integer :: nor,my_rank,nprocs,comm_self,g1mg2_idx
 real(dp) :: ucvol
 logical :: is_qeq0,use_MPI
 character(len=500) :: msg
!arrays
 integer :: omega_distrb(nomega)
 integer,allocatable :: istart(:),istop(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: eelf(:,:),tmp_eelf(:)
 complex(dpc),allocatable :: epsm_lf(:,:),epsm_nlf(:,:),tmp_lf(:),tmp_nlf(:)
 complex(dpc),allocatable :: buffer_lwing(:,:),buffer_uwing(:,:)
 complex(gwpc),allocatable :: kxcg_mat(:,:)

!bootstrap @WC
 integer :: istep,nstep
 logical :: converged
 real(dp) :: conv_err
 real(gwpc) :: chi00_head, fxc_head
 !real(gwpc) :: chi00_head, chi00rpa_head, fxc_head
 complex(gwpc),allocatable :: vfxc_boot(:,:), chi0_tmp(:,:), chi0_save(:,:,:)
 !complex(gwpc),allocatable :: fxc_lrc(:,:), vfxc_boot(:,:), chi0_tmp(:,:), chi0_save(:,:,:)
 complex(gwpc), ABI_CONTIGUOUS pointer :: vc_sqrt(:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (nI/=1.or.nJ/=1) then
   MSG_ERROR("nI or nJ=/1 not yet implemented")
 end if

 nprocs  = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! MG TODO We use comm_self for the inversion as the single precision version is not yet available
 comm_self = xmpi_comm_self

 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)

 is_qeq0 = (normv(Vcp%qibz(:,iqibz),gmet,'G')<GW_TOLQ0)

 omega_distrb = my_rank
 use_MPI = .FALSE.
 use_MPI = (nprocs>=nomega)  ! Parallelism is not used

 if (use_MPI) then
   ! * Initialize distribution table for frequencies.
   ABI_MALLOC(istart,(nprocs))
   ABI_MALLOC(istop,(nprocs))
   call xmpi_split_work2_i4b(nomega,nprocs,istart,istop)
   omega_distrb(:)=xmpi_undefined_rank
   do irank=0,nprocs-1
     i1 = istart(irank+1)
     i2 = istop (irank+1)
     if (i1<=i2) omega_distrb(i1:i2) = irank
   end do
   ABI_FREE(istart)
   ABI_FREE(istop)
 end if

 ! Initialize container for spectral results
 do nor=1,nomega
   if (ABS(AIMAG(omega(nor)))>1.e-3) EXIT
 end do
 nor=nor-1; if (nor==0) nor = 1 ! only imag !?

 if (dim_wing==3) then
   call wrtout(std_out,' Analyzing long wavelength limit for several q','COLL')
   call spectra_init(Spectra,nor,REAL(omega(1:nor)),Vcp%nqlwl,Vcp%qlwl)
   my_nqlwl = 1
   !my_nqlwl = dim_wing ! TODO
   !ABI_CHECK(dim_wing==SIZE(Vcp%vcqlwl_sqrt,DIM=2),"WRONG DIMS")
 else
   call spectra_init(Spectra,nor,REAL(omega(1:nor)),1,Vcp%qibz(:,iqibz))
   my_nqlwl = 1
 end if
 !
 ! NOTE: all processors have to perform this operation in order to have the
 !       epsm1 matrix when performing a sigma calculation starting with the file _SUS
 !
 ! Temporary arrays to store spectra.
 ABI_MALLOC(epsm_lf,(nomega,my_nqlwl))
 ABI_MALLOC(epsm_nlf,(nomega,my_nqlwl))
 ABI_MALLOC(eelf,(nomega,my_nqlwl))
 epsm_lf=czero; epsm_nlf=czero; eelf=zero

 ! Temporary arrays used to store output results.
 ABI_MALLOC(tmp_lf, (my_nqlwl))
 ABI_MALLOC(tmp_nlf, (my_nqlwl))
 ABI_MALLOC(tmp_eelf, (my_nqlwl))

 SELECT CASE (approx_type)

 CASE (0)
   ! * RPA: \tepsilon=1 - Vc^{1/2} chi0 Vc^{1/2}
   ! * vc_sqrt contains vc^{1/2}(q,G), complex-valued to allow for a possible cutoff.
   do io=1,nomega
     if (omega_distrb(io) == my_rank) then
       !write(std_out,*)"dim_wing",dim_wing
       call rpa_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,io),my_nqlwl,dim_wing,&
&        chi0_head(:,:,io),chi0_lwing(:,io,:),chi0_uwing(:,io,:),&
&        tmp_lf,tmp_nlf,tmp_eelf,comm_self)

         ! Store results.
         epsm_lf(io,:) = tmp_lf
         epsm_nlf(io,:) = tmp_nlf
         eelf(io,:) = tmp_eelf
     end if
   end do ! nomega

 CASE (1)
   ! Vertex correction from Adiabatic TDDFT. chi_{G1,G2} = [\delta -\chi0 (vc+kxc)]^{-1}_{G1,G3} \chi0_{G3,G2}
   ABI_CHECK(Vcp%nqlwl==1,"nqlwl/=1 not coded")
   ABI_CHECK(nkxc==1,"nkxc/=1 not coded")

   ! Make kxcg_mat(G1,G2) = kxcg(G1-G2) from kxcg defined on the FFT mesh.
   ABI_MALLOC_OR_DIE(kxcg_mat,(npwe,npwe), ierr)

   ierr=0
   do ig2=1,npwe
     do ig1=1,npwe
       g1mg2_idx = g2ifft(gvec(:,ig1)-gvec(:,ig2),ngfft)
       if (g1mg2_idx>0) then
         kxcg_mat(ig1,ig2) = kxcg(g1mg2_idx,1)
       else
         ierr=ierr+1
         kxcg_mat(ig1,ig2) = czero
       end if
     end do
   end do

   if (ierr/=0) then
     write(msg,'(a,i4,3a)')&
&     'Found ',ierr,' G1-G2 vectors falling outside the FFT box. ',ch10,&
&     'Enlarge the FFT mesh to get rid of this problem. '
     MSG_WARNING(msg)
   end if

   !FIXME "recheck TDDFT code and parallel"
   ABI_CHECK(nkxc==1,"nkxc/=1 not coded")
   do io=1,nomega
     if (omega_distrb(io) == my_rank) then
       call atddft_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,io),kxcg_mat,option_test,my_nqlwl,dim_wing,omega(io),&
&        chi0_head(:,:,io),chi0_lwing(:,io,:),chi0_uwing(:,io,:),tmp_lf,tmp_nlf,tmp_eelf,comm)

       ! Store results.
       epsm_lf(io,:) = tmp_lf
       epsm_nlf(io,:) = tmp_nlf
       eelf(io,:) = tmp_eelf
     end if
   end do

   ABI_FREE(kxcg_mat)
   do io=1,nomega
     write(msg,'(a,i4,a,2f9.4,a)')' Symmetrical epsilon^-1(G,G'') at the ',io,' th omega',omega(io)*Ha_eV,' [eV]'
     call wrtout(std_out,msg,'COLL')
     call print_arr(chi0(:,:,io),mode_paral='PERS')
   end do

 CASE (2)
   ! ADA nonlocal vertex correction contained in fxc_ADA
   MSG_WARNING('Entered fxc_ADA branch: EXPERIMENTAL!')
   ! Test that argument was passed
   if (.not.present(fxc_ADA)) then
     MSG_ERROR('make_epsm1_driver was not called with optional argument fxc_ADA')
   end if
   ABI_CHECK(Vcp%nqlwl==1,"nqlwl/=1 not coded")

   do io=1,nomega
     if (omega_distrb(io) == my_rank) then
       call atddft_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,io),fxc_ADA,option_test,my_nqlwl,dim_wing,omega(io),&
&        chi0_head(:,:,io),chi0_lwing(:,io,:),chi0_uwing(:,io,:),tmp_lf,tmp_nlf,tmp_eelf,comm)

       ! Store results.
       epsm_lf(io,:) = tmp_lf
       epsm_nlf(io,:) = tmp_nlf
       eelf(io,:) = tmp_eelf
     end if
   end do

   do io=1,nomega
     write(msg,'(a,i4,a,2f9.4,a)')' Symmetrical epsilon^-1(G,G'') at the ',io,' th omega',omega(io)*Ha_eV,' [eV]'
     call wrtout(std_out,msg,'COLL')
     call print_arr(chi0(:,:,io),mode_paral='PERS')
   end do

 CASE (4)
   !@WC bootstrap vertex correction, Sharma et al. PRL 107, 196401 (2011) [[cite:Sharma2011]]
   ABI_MALLOC_OR_DIE(vfxc_boot,(npwe*nI,npwe*nJ), ierr)
   ABI_MALLOC_OR_DIE(chi0_tmp,(npwe*nI,npwe*nJ), ierr)
   ABI_MALLOC_OR_DIE(chi0_save,(npwe*nI,npwe*nJ,nomega), ierr)

   if (iqibz==1) then
     vc_sqrt => Vcp%vcqlwl_sqrt(:,1)  ! Use Coulomb term for q-->0
   else
     vc_sqrt => Vcp%vc_sqrt(:,iqibz)
   end if

   chi0_save = chi0 ! a copy of chi0 (ks)
   nstep = 50 ! max iteration steps
   chi00_head = chi0(1,1,1)*vc_sqrt(1)**2
   fxc_head = czero; vfxc_boot = czero; chi0_tmp = czero
   epsm_lf = czero; epsm_nlf = czero; eelf = zero
   write(msg,'(a,2f10.6)') ' -> chi0_dft(head): ', chi00_head
   call wrtout(std_out,msg,'COLL')

   do istep=1,nstep
     chi0 = chi0_save
     do io=1,1 ! static
       !if (omega_distrb(io) == my_rank) then
       call atddft_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,io),vfxc_boot,0,my_nqlwl,dim_wing,omega(io),&
&       chi0_head(:,:,io),chi0_lwing(:,io,:),chi0_uwing(:,io,:),tmp_lf,tmp_nlf,tmp_eelf,comm_self)
       epsm_lf(io,:) = tmp_lf
       epsm_nlf(io,:) = tmp_nlf
       eelf(io,:) = tmp_eelf
       !end if
     end do

     conv_err = smallest_real
     do ig2=1,npwe*nJ
       do ig1=1,npwe*nI
         conv_err= MAX(conv_err, ABS(chi0(ig1,ig2,1) - chi0_tmp(ig1,ig2)))
       end do
     end do
     converged = (conv_err <= tol4)
     write(msg,'(a,i4,a,f10.6)') ' => bootstrap itr# ', istep, ', eps^-1 max error: ', conv_err
     call wrtout(std_out,msg,'COLL')
     write(msg,'(a,2f10.6)')  '    eps^-1(head):   ', chi0(1,1,1)
     call wrtout(std_out,msg,'COLL')
     write(msg,'(a,2f10.6)')  '    v^-1*fxc(head): ', fxc_head
     call wrtout(std_out,msg,'COLL')

     if (converged) then
       write(msg,'(a,i4,a)') ' => bootstrap fxc converged after ', istep, ' iterations'
       call wrtout(std_out,msg,'COLL')
       chi0 = chi0_save
       do io=1,nomega
         if (omega_distrb(io) == my_rank) then
           call atddft_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,io),vfxc_boot,option_test,my_nqlwl,dim_wing,omega(io),&
&           chi0_head(:,:,io),chi0_lwing(:,io,:),chi0_uwing(:,io,:),tmp_lf,tmp_nlf,tmp_eelf,comm_self)
           epsm_lf(io,:) = tmp_lf
           epsm_nlf(io,:) = tmp_nlf
           eelf(io,:) = tmp_eelf
         end if
       end do
       write(msg, '(a,2f10.6)') ' ->   eps^-1(head): ', chi0(1,1,1)
       call wrtout(std_out,msg,'COLL')
       write(msg,'(a,2f10.6)')  ' -> v^-1*fxc(head): ', fxc_head
       call wrtout(std_out,msg,'COLL')
       exit
     else if (istep < nstep) then
       chi0_tmp = chi0(:,:,1)
       vfxc_boot = chi0(:,:,1)/chi00_head ! full G vectors
       !vfxc_boot = czero; vfxc_boot(1,1) = chi0(1,1,1)/chi00_head ! head only
       fxc_head = vfxc_boot(1,1)
       do ig1=1,npwe
         vfxc_boot(ig1,:) = vc_sqrt(ig1)*vc_sqrt(:)*vfxc_boot(ig1,:)
       end do
     else
       write(msg,'(a,i4,a)') ' -> bootstrap fxc not converged after ', nstep, ' iterations'
       MSG_WARNING(msg)
       ! proceed to calculate the dielectric function even fxc is not converged
       chi0 = chi0_save
       do io=1,nomega
         if (omega_distrb(io) == my_rank) then
           call atddft_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,io),vfxc_boot,option_test,my_nqlwl,dim_wing,omega(io),&
&            chi0_head(:,:,io),chi0_lwing(:,io,:),chi0_uwing(:,io,:),tmp_lf,tmp_nlf,tmp_eelf,comm_self)
           epsm_lf(io,:) = tmp_lf
           epsm_nlf(io,:) = tmp_nlf
           eelf(io,:) = tmp_eelf
         end if
       end do
     end if
   end do

   ABI_FREE(chi0_tmp)
   ABI_FREE(chi0_save)
   ABI_FREE(vfxc_boot)

   do io=1,nomega
     write(msg,'(a,i4,a,2f9.4,a)')' Symmetrical epsilon^-1(G,G'') at the ',io,' th omega',omega(io)*Ha_eV,' [eV]'
     call wrtout(std_out,msg,'COLL')
     call print_arr(chi0(:,:,io),mode_paral='PERS')
   end do

 CASE(5)
   !@WC: one-shot scalar bootstrap approximation
   ABI_MALLOC_OR_DIE(vfxc_boot,(npwe*nI,npwe*nJ), ierr)
   ABI_MALLOC_OR_DIE(chi0_save,(npwe*nI,npwe*nJ,nomega), ierr)

   if (iqibz==1) then
     vc_sqrt => Vcp%vcqlwl_sqrt(:,1)  ! Use Coulomb term for q-->0
   else
     vc_sqrt => Vcp%vc_sqrt(:,iqibz)
   end if

   chi0_save = chi0 ! a copy of chi0
   fxc_head = czero; vfxc_boot = czero;
   epsm_lf = czero; epsm_nlf = czero; eelf = zero
   chi00_head = chi0(1,1,1)*vc_sqrt(1)**2
   write(msg,'(a,2f10.6)') ' -> chi0_dft(head): ',chi00_head
   call wrtout(std_out,msg,'COLL')

   fxc_head = vc_sqrt(1)**2/chi00_head + vc_sqrt(1)**2/chi00_head - vc_sqrt(1)**2
   fxc_head = 0.5*fxc_head + 0.5*sqrt(fxc_head**2 - 4.0*vc_sqrt(1)**4/(chi00_head*chi00_head))
   vfxc_boot(1,1) = fxc_head
   write(msg,'(a,2f10.6)') ' -> v^-1*fxc(head): ',fxc_head/vc_sqrt(1)**2
   call wrtout(std_out,msg,'COLL')

   chi0 = chi0_save
   do io=1,nomega
     if (omega_distrb(io) == my_rank) then
       call atddft_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,io),vfxc_boot,option_test,my_nqlwl,dim_wing,omega(io),&
&         chi0_head(:,:,io),chi0_lwing(:,io,:),chi0_uwing(:,io,:),tmp_lf,tmp_nlf,tmp_eelf,comm_self)
       epsm_lf(io,:) = tmp_lf
       epsm_nlf(io,:) = tmp_nlf
       eelf(io,:) = tmp_eelf
     end if
   end do
   write(msg,'(a,2f10.6)')  '    eps^-1(head):   ',chi0(1,1,1)
   call wrtout(std_out,msg,'COLL')

   ABI_FREE(chi0_save)
   ABI_FREE(vfxc_boot)

   do io=1,nomega
     write(msg,'(a,i4,a,2f9.4,a)')' Symmetrical epsilon^-1(G,G'') at the ',io,' th omega',omega(io)*Ha_eV,' [eV]'
     call wrtout(std_out,msg,'COLL')
     call print_arr(chi0(:,:,io),mode_paral='PERS')
   end do

CASE(6)
   !@WC: RPA bootstrap by Rigamonti et al. (PRL 114, 146402) [[cite:Rigamonti2015]]
   !@WC: and Berger (PRL 115, 137402) [[cite:Berger2015]]
   ABI_MALLOC_OR_DIE(vfxc_boot,(npwe*nI,npwe*nJ), ierr)
   ABI_MALLOC_OR_DIE(chi0_save,(npwe*nI,npwe*nJ,nomega), ierr)

   if (iqibz==1) then
     vc_sqrt => Vcp%vcqlwl_sqrt(:,1)  ! Use Coulomb term for q-->0
   else
     vc_sqrt => Vcp%vc_sqrt(:,iqibz)
   end if

   chi0_save = chi0 ! a copy of chi0
   fxc_head = czero; vfxc_boot = czero;
   epsm_lf = czero; epsm_nlf = czero; eelf = zero
   chi00_head = chi0(1,1,1)*vc_sqrt(1)**2
   write(msg,'(a,2f10.6)') ' -> chi0_dft(head): ',chi00_head
   call wrtout(std_out,msg,'COLL')

   ! static
   io = 1
   call atddft_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,io),vfxc_boot,0,my_nqlwl,dim_wing,omega(io),&
&    chi0_head(:,:,io),chi0_lwing(:,io,:),chi0_uwing(:,io,:),tmp_lf,tmp_nlf,tmp_eelf,comm_self)
   epsm_lf(1,:) = tmp_lf

   vfxc_boot(1,1) = 1.0/(chi00_head * epsm_lf(1,1))
   fxc_head = vfxc_boot(1,1)
   do ig1=1,npwe
     vfxc_boot(ig1,:) = vc_sqrt(ig1)*vc_sqrt(:)*vfxc_boot(ig1,:)
   end do
   write(msg,'(a,2f10.6)') ' -> v^-1*fxc(head): ',fxc_head
   call wrtout(std_out,msg,'COLL')

   chi0 = chi0_save
   do io=1,nomega
     if (omega_distrb(io) == my_rank) then
       call atddft_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,io),vfxc_boot,option_test,my_nqlwl,dim_wing,omega(io),&
&        chi0_head(:,:,io),chi0_lwing(:,io,:),chi0_uwing(:,io,:),tmp_lf,tmp_nlf,tmp_eelf,comm_self)
       epsm_lf(io,:) = tmp_lf
       epsm_nlf(io,:) = tmp_nlf
       eelf(io,:) = tmp_eelf
     end if
   end do
   write(msg,'(a,2f10.6)')  '    eps^-1(head):   ',chi0(1,1,1)
   call wrtout(std_out,msg,'COLL')

   ABI_FREE(chi0_save)
   ABI_FREE(vfxc_boot)

   do io=1,nomega
     write(msg,'(a,i4,a,2f9.4,a)')' Symmetrical epsilon^-1(G,G'') at the ',io,' th omega',omega(io)*Ha_eV,' [eV]'
     call wrtout(std_out,msg,'COLL')
     call print_arr(chi0(:,:,io),mode_paral='PERS')
   end do

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong approx_type:',itoa(approx_type)))
 END SELECT

 if (use_MPI) then
   ! Collect results on each node.
   ABI_MALLOC(buffer_lwing, (size(chi0_lwing,dim=1), size(chi0_lwing, dim=3)))
   ABI_MALLOC(buffer_uwing, (size(chi0_uwing,dim=1), size(chi0_uwing, dim=3)))

   do io=1,nomega
     if (omega_distrb(io)/=my_rank) then
       ! Zero arrays.
       chi0(:,:,io) = czero_gw
       if(dim_wing>0) then
          chi0_lwing(:,io,:) = zero
          chi0_uwing(:,io,:) = zero
          chi0_head(:,:,io)  = czero
       endif
       epsm_lf(io,:) = czero
       epsm_nlf(io,:) = czero
       eelf(io,:) = zero
     end if

     call xmpi_sum(chi0(:,:,io), comm,ierr)

     if (dim_wing>0) then
        ! Build contiguous arrays
        buffer_lwing = chi0_lwing(:,io,:)
        buffer_uwing = chi0_uwing(:,io,:)
        call xmpi_sum(buffer_lwing,comm,ierr)
        call xmpi_sum(buffer_uwing,comm,ierr)
        chi0_lwing(:,io,:) = buffer_lwing
        chi0_uwing(:,io,:) = buffer_uwing
        if (size(chi0_head(:,:,io))/= zero) then
          call xmpi_sum(chi0_head(:,:,io),comm,ierr)
        end if
     end if

   end do ! iw

   call xmpi_sum(epsm_lf, comm,ierr )
   call xmpi_sum(epsm_nlf,comm,ierr)
   call xmpi_sum(eelf,    comm,ierr)
   ABI_FREE(buffer_lwing)
   ABI_FREE(buffer_uwing)
 end if

 ! Save results in Spectra%, mind the slicing.
 Spectra%emacro_nlf(:,:) = epsm_nlf(1:nor,:)
 Spectra%emacro_lf (:,:) = epsm_lf (1:nor,:)
 Spectra%eelf      (:,:) = eelf    (1:nor,:)

 ABI_FREE(epsm_lf)
 ABI_FREE(epsm_nlf)
 ABI_FREE(eelf)
 ABI_FREE(tmp_lf)
 ABI_FREE(tmp_nlf)
 ABI_FREE(tmp_eelf)

 DBG_EXIT("COLL")

end subroutine make_epsm1_driver
!!***

!----------------------------------------------------------------------

!!****f* m_screening/rpa_symepsm1
!! NAME
!! rpa_symepsm1
!!
!! FUNCTION
!!  Calculate RPA $\tilde\epsilon^{-1}$
!!
!!  Use a special treatment of non-Analytic behavior of heads and wings in reciprocal space
!!  calculating these quantities for different small q-directions specified by the user
!!  (Not yet operative)
!!
!! INPUTS
!!  iqibz=index of the q-point in the array Vcp%qibz where epsilon^-1 has to be calculated
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!  npwe=Number of G-vectors in chi0.
!!  nI,nJ=Number of rows/columns in chi0_ij (1,1 in collinear case)
!!  dim_wing=Dimension of the wings (0 or 3 if q-->0)
!!  chi0_head(dim_wing,dim_wing)=Head of of chi0 (only for q-->0)
!!  chi0_lwing(npwe*nI,dim_wing)=Lower wings of chi0 (only for q-->0)
!!  chi0_uwing(npwe*nJ,dim_wing)=Upper wings of chi0 (only for q-->0)
!!  comm=MPI communicator.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chi0(npwe*nI,npwe*nJ): in input the irreducible polarizability, in output
!!   the symmetrized inverse dielectric matrix.
!!
!! PARENTS
!!      m_screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine rpa_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0,my_nqlwl,dim_wing,chi0_head,chi0_lwing,chi0_uwing,epsm_lf,epsm_nlf,eelf,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,nI,nJ,npwe,dim_wing,my_nqlwl,comm
 type(vcoul_t),target,intent(in) :: Vcp
!arrays
 complex(gwpc),intent(inout) :: chi0(npwe*nI,npwe*nJ)
 complex(dpc),intent(inout) :: chi0_lwing(:,:) !(npwe*nI,dim_wing)
 complex(dpc),intent(inout) :: chi0_uwing(:,:) !(npwe*nJ,dim_wing)
 complex(dpc),intent(inout) :: chi0_head(:,:) !(dim_wing,dim_wing)
 real(dp),intent(out) :: eelf(my_nqlwl)
 complex(dpc),intent(out) :: epsm_lf(my_nqlwl),epsm_nlf(my_nqlwl)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,prtvol=0
 integer :: ig1,ig2,iqlwl,my_rank,nprocs
 real(dp) :: ucvol
 logical :: is_qeq0
 !character(len=500) :: msg
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 complex(gwpc), ABI_CONTIGUOUS pointer :: vc_sqrt(:)
 complex(gwpc),allocatable :: chi0_save(:,:)

! *************************************************************************

 ABI_UNUSED(chi0_head(1,1))

 if (nI/=1.or.nJ/=1) then
   MSG_ERROR("nI or nJ=/1 not yet implemented")
 end if

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)

 is_qeq0 = (normv(Vcp%qibz(:,iqibz),gmet,'G')<GW_TOLQ0)
 if (is_qeq0) then
   ABI_CHECK(iqibz==1,"q is 0 but iq_ibz /= 1")
 end if
 !
 if (my_nqlwl>1) then
   ABI_MALLOC(chi0_save,(npwe*nI,npwe*nJ))
   chi0_save = chi0
 end if
 !
 ! Symmetrized RPA epsilon = 1 - Vc^{1/2} chi0 Vc^{1/2}
 !   * vc_sqrt contains vc^{1/2}(q,G), complex-valued to allow for a possible cutoff.
 !
 ! * Loop over small q"s (if any) to treat the nonanalytical behavior.
 do iqlwl=my_nqlwl,1,-1
   !
   if (my_nqlwl>1) then
     chi0(:,:) = chi0_save           ! restore pristine polarizability
     chi0(:,1) = chi0_lwing(:,iqlwl) ! change the wings
     chi0(1,:) = chi0_uwing(:,iqlwl)
   end if

   if (iqibz==1) then
     vc_sqrt => Vcp%vcqlwl_sqrt(:,iqlwl)  ! Use Coulomb term for q-->0
   else
     vc_sqrt => Vcp%vc_sqrt(:,iqibz)
   end if

   do ig2=1,npwe*nJ
     do ig1=1,npwe*nI
       chi0(ig1,ig2)=-vc_sqrt(ig1)*chi0(ig1,ig2)*vc_sqrt(ig2)
     end do
     chi0(ig2,ig2)=one+chi0(ig2,ig2)
   end do

   epsm_nlf(iqlwl)=chi0(1,1) ! * chi0, now contains \tepsilon.

   if (prtvol > 0) then
     call wrtout(std_out,' Symmetrical epsilon(G,G'') ','COLL')
     call print_arr(chi0, unit=std_out)
   end if
   !
   ! === Invert tepsilon and calculate macroscopic dielectric constant ===
   ! * epsm_lf(w)=1/epsm1(G=0,Gp=0,w).
   ! * Since G=Gp=0 there is no difference btw symmetrical and not symmetrical.
   !
   call xginv(chi0,npwe,comm=comm)

   epsm_lf(iqlwl) = one/chi0(1,1)
   eelf(iqlwl) = -AIMAG(chi0(1,1))

   if (prtvol > 0) then
     call wrtout(std_out," Symmetrical epsilon^-1(G,G'')",'COLL')
     call print_arr(chi0, unit=std_out)
   end if
   !
   ! Save wings of e^-1 overwriting input values.
   if (dim_wing>0.and..FALSE.) then
     chi0_lwing(:,iqlwl) = chi0(:,1)
     chi0_uwing(:,iqlwl) = chi0(1,:)
   end if

 end do !iqlwl

 if (allocated(chi0_save))  then
   ABI_FREE(chi0_save)
 end if

end subroutine rpa_symepsm1
!!***

!----------------------------------------------------------------------

!!****f* m_screening/atddft_symepsm1
!! NAME
!! atddft_symepsm1
!!
!! FUNCTION
!!  Calculate $\tilde\epsilon^{-1}$ using ALDA within TDDFT
!!
!!  2) Use a special treatment of non-Analytic behavior of heads and wings in reciprocal space
!!     calculating these quantities for different small q-directions specified by the user
!!     (Not yet operative)
!!
!!  Output the electron energy loss function and the macroscopic dielectric function with and
!!  without local field effects (only if non-zero real frequencies are available)
!!
!! INPUTS
!!  iqibz=index of the q-point in the array Vcp%qibz where epsilon^-1 has to be calculated
!!  npwe=Number of G-vectors in chi0.
!!  nI,nJ=Number of rows/columns in chi0_ij (1,1 in collinear case)
!!  dim_wing=Dimension of the wings (0 or 3 if q-->0)
!!  option_test=Only for TDDFT:
!!   == 0 for TESTPARTICLE ==
!!   == 1 for TESTELECTRON ==
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!   %nqibz=Number of q-points.
!!   %qibz(3,nqibz)=q-points in the IBZ.
!!  comm=MPI communicator.
!!  chi0_lwing(npwe*nI,dim_wing)=Lower wings of chi0 (only for q-->0)
!!  chi0_uwing(npwe*nJ,dim_wing)=Upper wings of chi0 (only for q-->0)
!!  chi0_head(dim_wing,dim_wing)=Head of of chi0 (only for q-->0)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chi0(npwe*nI,npwe*nJ): in input the irreducible polarizability, in output
!!   the symmetrized inverse dielectric matrix.
!!
!! PARENTS
!!      m_screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine atddft_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0,kxcg_mat,option_test,my_nqlwl,dim_wing,omega,&
& chi0_head,chi0_lwing,chi0_uwing,epsm_lf,epsm_nlf,eelf,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,nI,nJ,npwe,dim_wing,my_nqlwl
 integer,intent(in) :: option_test,comm
 type(vcoul_t),target,intent(in) :: Vcp
!arrays
 complex(gwpc),intent(in) :: kxcg_mat(npwe,npwe)
 complex(dpc),intent(in) :: omega
 complex(dpc),intent(inout) :: chi0_lwing(npwe*nI,dim_wing)
 complex(dpc),intent(inout) :: chi0_uwing(npwe*nJ,dim_wing)
 complex(dpc),intent(inout) :: chi0_head(dim_wing,dim_wing)
 complex(gwpc),intent(inout) :: chi0(npwe*nI,npwe*nJ)
 real(dp),intent(out) :: eelf(my_nqlwl)
 complex(dpc),intent(out) :: epsm_lf(my_nqlwl),epsm_nlf(my_nqlwl)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,prtvol=0
 integer :: ig1,ig2,my_rank,nprocs,ierr
 real(dp) :: ucvol
 logical :: is_qeq0
 character(len=500) :: msg
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 complex(gwpc),allocatable :: chitmp(:,:)
 complex(gwpc), ABI_CONTIGUOUS pointer :: vc_sqrt(:)

! *************************************************************************

 ABI_UNUSED(chi0_head(1,1))
 ABI_UNUSED(chi0_lwing(1,1))
 ABI_UNUSED(chi0_uwing(1,1))

 if (nI/=1.or.nJ/=1) then
   MSG_ERROR("nI or nJ=/1 not yet implemented")
 end if

 ABI_CHECK(Vcp%nqlwl==1,"nqlwl/=1 not coded")
 ABI_CHECK(my_nqlwl==1,"my_nqlwl/=1 not coded")

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)

 is_qeq0 = (normv(Vcp%qibz(:,iqibz),gmet,'G')<GW_TOLQ0)

 if (iqibz==1) then
   !%vc_sqrt => Vcp%vcqlwl_sqrt(:,iqlwl)  ! Use Coulomb term for q-->0
   vc_sqrt => Vcp%vcqlwl_sqrt(:,1)  ! TODO add treatment of non-Analytic behavior
 else
   vc_sqrt => Vcp%vc_sqrt(:,iqibz)
 end if

 write(msg,'(a,f8.2,a)')" chitmp requires: ",npwe**2*gwpc*b2Mb," Mb"
 ABI_MALLOC_OR_DIE(chitmp,(npwe,npwe), ierr)
 !
 ! * Calculate chi0*fxc.
 chitmp = MATMUL(chi0,kxcg_mat)
 ! * First calculate the NLF contribution
 do ig1=1,npwe
   do ig2=1,npwe
     chitmp(ig1,ig2)=-chitmp(ig1,ig2)
   end do
   chitmp(ig1,ig1)=chitmp(ig1,ig1)+one
 end do

 call xginv(chitmp,npwe,comm=comm)

 chitmp = MATMUL(chitmp,chi0)
 !if (.not. ABS(REAL(omega))> tol3) call hermitianize(chitmp,"All")
 chitmp(1,1)=-vc_sqrt(1)*chitmp(1,1)*vc_sqrt(1)
 chitmp(1,1)=chitmp(1,1)+one

 epsm_nlf(1)=chitmp(1,1)

 chitmp = MATMUL(chi0,kxcg_mat)
 ! * Calculate (1-chi0*Vc-chi0*Kxc) and put it in chitmp.
 do ig1=1,npwe
   do ig2=1,npwe
     chitmp(ig1,ig2)=-chitmp(ig1,ig2)-chi0(ig1,ig2)*vc_sqrt(ig2)**2
   end do
   chitmp(ig1,ig1)=chitmp(ig1,ig1)+one
 end do

 ! * Invert (1-chi0*Vc-chi0*Kxc) and Multiply by chi0.
 call xginv(chitmp,npwe,comm=comm)
 chitmp=MATMUL(chitmp,chi0)

 ! * Save result, now chi0 contains chi.
 chi0=chitmp

 select case (option_test)
 case (0)
   ! Symmetrized TESTPARTICLE epsilon^-1
   call wrtout(std_out,' Calculating TESTPARTICLE epsilon^-1(G,G") = 1 + Vc*chi','COLL')
   do ig1=1,npwe
     chi0(ig1,:)=(vc_sqrt(ig1)*vc_sqrt(:))*chi0(ig1,:)
     chi0(ig1,ig1)=one+chi0(ig1,ig1)
   end do

 case (1)
   ! Symmetrized TESTELECTRON epsilon^-1
   call wrtout(std_out,' Calculating TESTELECTRON epsilon^-1(G,G") = 1 + (Vc + fxc)*chi',"COLL")
   chitmp=MATMUL(kxcg_mat,chi0)

   ! Perform hermitianization, only valid along the imaginary axis.
   if (.not. ABS(REAL(omega))> tol3) call hermitianize(chitmp,"All")

   do ig1=1,npwe
     chi0(ig1,:)=(vc_sqrt(ig1)*vc_sqrt(:))*chi0(ig1,:)+chitmp(ig1,:)
     chi0(ig1,ig1)=one+chi0(ig1,ig1)
   end do

 case default
   MSG_BUG(sjoin('Wrong option_test:',itoa(option_test)))
 end select

 ABI_FREE(chitmp)
 !
 ! === chi0 now contains symmetrical epsm1 ===
 ! * Calculate macroscopic dielectric constant epsm_lf(w)=1/epsm1(G=0,Gp=0,w).
 epsm_lf(1) =  one/chi0(1,1)
 eelf   (1) = -AIMAG(chi0(1,1))

 if (prtvol > 0) then
   write(msg,'(a,2f9.4,a)')' Symmetrical epsilon^-1(G,G'') at omega',omega*Ha_eV,' [eV]'
   call wrtout(std_out,msg,'COLL')
   call print_arr(chi0,unit=std_out)
 end if

end subroutine atddft_symepsm1
!!***

!----------------------------------------------------------------------

!!****f* m_screening/mkem1_q0
!! NAME
!! mkem1_q0
!!
!! FUNCTION
!!   This routine construct the microscopic dieletric matrix for q-->0 starting from the heads, wings and the body
!!   of the irreducible polarizability. Afterwards it calculates the symmetrized inverse dieletric matrix
!!   via a block wise inversion thus obtaining the heads and the wings of e^{-1} that can be
!!   used to describe the non-analytic behavior for q-->0.
!!
!! INPUTS
!! npwe=Number of Gs used to describe chi0
!! nomega=Number of frequencies in chi0.
!! n1,n2=Factors used to define the same of the chi0 matrix (1,1 if collinear, the typical case)
!! Cryst<crystal_t>=Structure describing the crystal structure.
!! Vcp<vcoul_t>=datatypes gathering info on the Coulomb term
!! gvec(3,npwe)=G-vector for chi0 in reduced coordinates.
!! comm=MPI communicator
!!
!! OUTPUT
!! eps_head(3,3,nomega)=The macroscopic dieletric tensor in reduced coordinates.
!!   The dieletric matrix along versor \hat q can be obtained with
!!     e(\hat q) = \hat q.eps_head \hat q if all quantities are given in Cartesian coordinates.
!!
!! SIDE EFFECTS
!! chi0(npwe*n1,npwe*n2,nomega)= Input: polarizability. output: inverse dieletric matrix (only the body is used)
!! chi0_lwing(npwe*n1,nomega,3)
!! chi0_uwing(npwe*n2,nomega,3)  Input:  the lower and upper wings of the polarizability
!!                               Output: the "lower" and "upper" wings of the inverse dieletric matrix. See notes below.
!! chi0_head(3,3,nomega)= Input: the polarizability tensor in Cartesian coordinates.
!!                        Ouput: The "head" of the inverse dieletric matrix. See notes below.
!!
!! NOTES
!!  Matrix inversion in block form.
!!
!!         1  n-1
!!  M =  | c  u^t| 1     ==>   M^{-1} =  |  1/k          -u^t A^{-1}/k                    |
!!       | v  A  | n-1                   | -A^{-1} v/k    A^{-1} + (A^{-1}v u^t A^{-1})/k |
!!
!!                             where k = c - u^t A^{-1} v
!!
!!  Let q be a versor in reciprocal space, the symmetrized dielectric matrix with bare coulomb interaction
!!  can be written as
!!
!!  \tilde\epsilon = | q.Eq      q.Wl(G2) |  where   E_ij = \delta_ij -4\pi chi0_head_ij
!!                   | q.Wl(G1)  B(G1,G2  |          Wl(G1) = -4\pi chi0_lwing(G1)
!!                                                   Wu(G2) = -4\pi chi0_uwing(G1)
!!  therefore, in Cartesian coordinates, we have:
!!
!!  1) e^{-1}_{00}(q) = [ q_i q_j (E_{ij} - \sum_{GG'} Wu_i(G)a B_{GG'}^{-1} Wl_j(G')) ]^{-1} = 1/(q.Lq)
!!
!!  2) e^{-1}_{0G'}(q) = -e^{-1}_{00}(q) [ \sum_{iG} q_i Wu_i(G)a B_{GG'}^{-1} ] = (q.Su) /(q.Lq)
!!
!!  3) e^{-1}_{G0}(q)  = -e^{-1}_{00}(q) [ \sum_{iG'} q_i B_{GG'}^{-1} Wl_i(G') ] = (q.Sl) /(q.Lq)
!!
!!  4) e^{-1}_{GG'}(q) = B_{GG'}^{-1} +
!!     [ \sum_{ij} q_i q_j ( \sum_T B^{-1}_{GT}^{-1} Wl_i(T)) (\sum_T' Wu_j(T') B^{-1}_{T'G'}^{-1} ] / (q.Lq)
!!
!!  where Su(G,3) and Sl(G,3) are the "upper" and "lower" wings of the inverse dielectric matrix and
!!  L is the inverse dielectric tensor. Similar equations hold even if vectors and tensors are given in terms
!!  of the reciprocal lattice vectors provided that the metric tensor is taken into account.
!!  The main difference is in the expression for the tensor as only one metric tensor can be
!!  absorbed in the scalar product, the second metric multiplies one of the wings.
!!
!!  *) The present implementation assumes that no cutoff technique is used in the Coulomb term.
!!
!!  *) Once the tensor in know it is possible to average the quadratic form on the sphere exactly.
!!  In Cartesian coordinates one obtains.
!!
!!    \dfrac{1}{4\pi} \int v.Tv d\Omega = Trace(T)/3
!!
!!  For the inverse dielectric matrix we have to resort to a numerical integration
!!
!! PARENTS
!!      m_screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkem1_q0(npwe,n1,n2,nomega,Cryst,Vcp,gvec,chi0_head,chi0_lwing,chi0_uwing,chi0,eps_head,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwe,nomega,n1,n2,comm
 type(crystal_t),intent(in) :: Cryst
 type(vcoul_t),intent(in) :: Vcp
!arrays
 integer,intent(in) :: gvec(3,npwe)
 complex(gwpc),intent(in) :: chi0(npwe*n1,npwe*n2,nomega)
 complex(dpc),intent(inout) :: chi0_lwing(npwe*n1,nomega,3)
 complex(dpc),intent(inout) :: chi0_uwing(npwe*n2,nomega,3)
 complex(dpc),intent(inout) :: chi0_head(3,3,nomega)
 complex(dpc),intent(out) :: eps_head(3,3,nomega)

!Local variables ------------------------------
!scalars
 integer :: iw,ig,ig1,ig2,idir,jdir
!arrays
 real(dp),allocatable :: modg_inv(:)
 complex(dpc),allocatable :: eps_lwing(:,:),eps_uwing(:,:),eps_body(:,:),cvec(:)

!************************************************************************

 ABI_CHECK(npwe/=1,"npwe must be >1")
 ABI_UNUSED(comm)

 ! Precompute 1/|G|.
 ABI_MALLOC(modg_inv,(npwe-1))
 do ig=1,npwe-1
   modg_inv(ig) = one/normv(gvec(:,ig+1),Cryst%gmet,'G')
 end do

 ABI_MALLOC(eps_uwing,((npwe-1)*n1,3))
 ABI_MALLOC(eps_lwing,((npwe-1)*n2,3))
 ABI_MALLOC(eps_body,(npwe-1,npwe-1))
 ABI_MALLOC(cvec,(npwe-1))

 do iw=1,nomega
   !
   ! Head and wings of the symmetrized epsilon.
   eps_head(:,:,iw) = -four_pi*chi0_head(:,:,iw)
   do idir=1,3
     eps_head(idir,idir,iw) = one + eps_head(idir,idir,iw)
     eps_lwing(:,idir) = -four_pi * modg_inv * chi0_lwing(2:,iw,idir)
     eps_uwing(:,idir) = -four_pi * modg_inv * chi0_uwing(2:,iw,idir)
     !eps_lwing(:,idir) = -chi0_lwing(2:,iw,idir) * SQRT(four_pi) * Vcp%vcqlwl_sqrt(2:npwe,1)
     !eps_uwing(:,idir) = -chi0_uwing(2:,iw,idir) * SQRT(four_pi) * Vcp%vcqlwl_sqrt(2:npwe,1)
   end do

   write(std_out,*)" espilon head"
   call print_arr(eps_head(:,:,iw))
   !
   ! Construct the body of the symmetrized epsilon then invert it.
   do ig2=1,npwe-1
     do ig1=1,npwe-1
       eps_body(ig1,ig2) = -four_pi * modg_inv(ig1)*chi0(ig1+1,ig2+1,iw )*modg_inv(ig2)
       !eps_body(ig1,ig2) = -Vcp%vcqlwl_sqrt(ig1+1,1)*chi0(ig1+1,ig2+1,iw)* Vcp%vcqlwl_sqrt(ig2+1,1)
     end do
     eps_body(ig2,ig2) = one + eps_body(ig2,ig2)
   end do

   call xginv(eps_body,npwe-1)
   !
   ! Overwrite chi0_head and chi0_wings with the head and the wings of the inverse dielectric matrix.
   do jdir=1,3
     !
     ! Head.
     cvec=czero
     do idir=1,3
       cvec = cvec + two_pi*Cryst%gmet(jdir,idir)*MATMUL(eps_body,eps_lwing(:,idir)) ! as we work in reciprocal coords.
     end do
     !cvec = MATMUL(eps_body,eps_lwing(:,jdir))
     do idir=1,3
       chi0_head(idir,jdir,iw) = eps_head(idir,jdir,iw) - xdotu(npwe-1,eps_uwing(:,idir),1,cvec,1)
     end do
     !
     ! Now the wings.
     chi0_uwing(2:,iw,jdir) = -MATMUL(eps_uwing(:,jdir),eps_body)
     chi0_lwing(2:,iw,jdir) = -MATMUL(eps_body,eps_lwing(:,jdir))
     !
   end do !jdir

   call wrtout(std_out, "espilon^1 head after block inversion", "COLL")
   call print_arr(chi0_head(:,:,iw))
   !
   ! Change the body but do not add the corrections due to the head and the wings.
   ! since they can be obtained on the fly from eps_body and the wings of eps^{-1}.
   !%chi0(2:,2:,iw) = eps_body
 end do !iw

 ABI_FREE(modg_inv)
 ABI_FREE(cvec)
 ABI_FREE(eps_lwing)
 ABI_FREE(eps_uwing)
 ABI_FREE(eps_body)

 RETURN
 ABI_UNUSED(Vcp%ng)

end subroutine mkem1_q0
!!***

!----------------------------------------------------------------------

!!****f* m_screening/lebedev_laikov_int
!! NAME
!!  lebedev_laikov_int
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine lebedev_laikov_int()

!Arguments ------------------------------------

!Local variables-------------------------------
!scalars
 integer :: on,npts,ii,ll,mm,lmax,leb_idx !ierr,
 real(dp) :: accuracy
 complex(dpc) :: ang_int
!arrays
 real(dp) :: cart_vpt(3) !,real_pars(0)
 real(dp),allocatable :: vx(:),vy(:),vz(:),ww(:)
 complex(dpc) :: tensor(3,3),cplx_pars(9)
 complex(dpc),allocatable :: ref_func(:),expd_func(:) !tmp_momenta(:)

! *************************************************************************

 MSG_ERROR("lebedev_laikov_int is still under development")

 !tensor=RESHAPE((/4.0,2.0,4.0,0.5,2.1,0.0,5.4,2.1,5.0/),(/3,3/))
 tensor=RESHAPE((/4.0,0.0,0.0,0.0,4.0,0.0,0.0,0.0,5.0/),(/3,3/))
 !tensor=RESHAPE((/1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0/),(/3,3/))

 npts=26
 ABI_MALLOC(vx,(npts))
 ABI_MALLOC(vy,(npts))
 ABI_MALLOC(vz,(npts))
 ABI_MALLOC(ww,(npts))

 !call LD0026(vx,vy,vz,ww,on)

 ang_int=czero
 do ii=1,npts
   cart_vpt = [vx(ii),vy(ii),vz(ii)]
   ang_int = ang_int + ww(ii)*DOT_PRODUCT(cart_vpt,MATMUL(tensor,cart_vpt))
 end do

 !write(std_out,*)"quadratic form associated to tensor=",tensor
 write(std_out,*)"on ang_int",on,ang_int

 ABI_FREE(vx)
 ABI_FREE(vy)
 ABI_FREE(vz)
 ABI_FREE(ww)

 !call init_lebedev_gridset()
 cplx_pars = RESHAPE(tensor,(/9/)); accuracy=tol10

 ! This is the function to be expanded evaluated on the lebedev_laikov grid of index leb_idx
 leb_idx=3; npts=lebedev_npts(leb_idx)
 ABI_MALLOC(ref_func,(npts))
 do ii=1,npts
   !cart_vpt = Lgridset(leb_idx)%versor(:,ii)
   ref_func(ii) = one/DOT_PRODUCT(cart_vpt,MATMUL(tensor,cart_vpt))
 end do

 ! Calculate the expansion in angular momenta of 1/{q.Tq}.
 ! Only even l-components contribute thanks to the parity of the integrand.
 ! tol6 seems to be an acceptable error, convergence wrt lmax is very slow even for simple tensors.
 ABI_MALLOC(expd_func,(npts))
 expd_func=czero
 lmax=10
 do ll=0,lmax,2
   !allocate(tmp_momenta(-ll:ll))
   do mm=-ll,ll
     ! MG: Commented becase it causes problems with the new version of abilint
     !call lebedev_quadrature(ylmstar_over_qTq,(/ll,mm/),real_pars,cplx_pars,ang_int,ierr,accuracy)
     write(std_out,*)ll,mm,ang_int
     !tmp_momenta(mm) = ang_int
     do ii=1,npts
       !cart_vpt = Lgridset(leb_idx)%versor(:,ii)
       expd_func(ii) = expd_func(ii) + four_pi*ang_int*ylmc(ll,mm,cart_vpt)
     end do
   end do
   !deallocate(tmp_momenta)
   write(std_out,*)"Error in angular expansion at l=",ll," is ",MAXVAL(ABS(expd_func-ref_func))
 end do

!BEGINDEBUG
! do ii=1,npts
!   write(777,*)ref_func(ii)
!   write(778,*)expd_func(ii)
! end do
!ENDDEBUG

 ABI_FREE(expd_func)
 ABI_FREE(ref_func)

 MSG_ERROR("Exiting from lebedev_laikov_int")

end subroutine lebedev_laikov_int
!!***

!----------------------------------------------------------------------

!!****f* m_screening/ylmstar_over_qTq
!! NAME
!!  ylmstar_over_qTq
!!
!! FUNCTION
!!  Return Ylm(q)^*/(q.Tq) where q is a versor in Cartesian coordinates.
!!  and Ylm is a complex spherical Harmonics whose index (l,m) are
!!  passed via int_pars(1:2). T is a tensore in Cartesian coordinates
!!  passed via cplx_pars(1:9).
!!
!! INPUTS
!!  cart_vers(3)=Cartesian components of the versor
!!  int_pars(1:2)=(l,m) indeces in Ylm. l>=0 and  m \in [-l,l]
!!  cplx_pars(1:9)=Tensor T in Cartesian coordinates.
!!  real_pars=Not used.
!!
!! OUTPUT
!!  Value of Ylm(q)^*/(q.Tq)
!!
!! PARENTS
!!
!! SOURCE

function ylmstar_over_qTq(cart_vers,int_pars,real_pars,cplx_pars)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: cart_vers(3)
 integer,intent(in) :: int_pars(:)
 real(dp),intent(in) :: real_pars(:)
 complex(dpc),intent(in) :: cplx_pars(:)
 complex(dpc) :: ylmstar_over_qTq
!arrays

!Local variables-------------------------------
!scalars
 integer :: ll,mm
!arrays
 complex(dpc) :: tensor(3,3)

! *************************************************************************

 tensor = RESHAPE(cplx_pars(1:9),(/3,3/))
 ll = int_pars(1) ! ll starts from zero.
 mm = int_pars(2) ! m \in [-l,l]

 ylmstar_over_qTq = CONJG(ylmc(ll,mm,cart_vers))/DOT_PRODUCT(cart_vers,MATMUL(tensor,cart_vers))

 RETURN
 ABI_UNUSED(real_pars(1))

end function ylmstar_over_qTq
!!***

!----------------------------------------------------------------------

!!****f* m_screening/ylmstar_wtq_over_qTq
!! NAME
!!  ylmstar_wtq_over_qTq
!!
!! FUNCTION
!!  Return Ylm(q)^* weight(q)/(q.Tq) where q is a versor in Cartesian coordinates.
!!  Ylm is a complex spherical Harmonics whose index (l,m) are
!!  passed via int_pars(1:2). T is a tensor in Cartesian coordinates
!!  passed via cplx_pars(1:9). weight(q) is the weighting function giving
!!  the length of the vector parallel to versor q that connects the origin
!!  of the lattice to one of the boundaries of the small cell centered at Gamma
!!
!! INPUTS
!!  cart_vers(3)=Cartesian components of the versor
!!  int_pars(1:2)=(l,m) indeces in Ylm. l>=0 and  m \in [-l,l]
!!  cplx_pars(1:9)=Tensor T in Cartesian coordinates.
!!  real_pars(1:9)=The Cartesian vectors defining the small box centered around gamma point
!!    when referred to this vectors the points in the box are given by {(x,y,z) | x,y,z \in [-1,1]}.
!!
!! OUTPUT
!!  Value of Ylm(q)^* weigh(q)/(q.Tq)
!!
!! PARENTS
!!
!! SOURCE

function ylmstar_wtq_over_qTq(cart_vers,int_pars,real_pars,cplx_pars)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: cart_vers(3)
 integer,intent(in) :: int_pars(:)
 real(dp),intent(in) :: real_pars(:)
 complex(dpc),intent(in) :: cplx_pars(:)
 complex(dpc) :: ylmstar_wtq_over_qTq
!arrays

!Local variables-------------------------------
!scalars
 integer :: ll,mm
 real(dp) :: wtq
!arrays
 real(dp) :: gprimd(3,3),rprimd(3,3),red_vers(3)
 complex(dpc) :: tensor(3,3)

! *************************************************************************

 MSG_ERROR("Work in progress")
 ! box_len has to be tested

 gprimd = RESHAPE(real_pars(1:9),(/3,3/))
 red_vers = MATMUL(rprimd,cart_vers)
 wtq = box_len(red_vers,gprimd)

 tensor = RESHAPE(cplx_pars(1:9),(/3,3/))
 ll = int_pars(1) ! true ll i.e. not shifted
 mm = int_pars(2)

 ylmstar_wtq_over_qTq = CONJG(ylmc(ll,mm,cart_vers))*wtq/DOT_PRODUCT(cart_vers,MATMUL(tensor,cart_vers))

end function ylmstar_wtq_over_qTq
!!***

!----------------------------------------------------------------------

!!****f* m_screening/mdielf_bechstedt
!! NAME
!!  mdielf_bechstedt
!!
!! FUNCTION
!!  Calculates the model dielectric function for the homogeneous system
!!  as proposed by F. Bechstedt, in Solid State Commun. 84, 765 1992.
!!
!! INPUTS
!!  eps_inf=Dielectric constant of the material
!!  qnrm=The modulus of the q-point.
!!  rhor=The local value of the density
!!
!! PARENTS
!!
!! SOURCE

elemental function mdielf_bechstedt(eps_inf,qnrm,rhor) result(mdielf)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: eps_inf,qnrm,rhor
 real(dp) :: mdielf

! *************************************************************************

 mdielf = one + &
&  one / ( one/(eps_inf-one) + (qnrm/k_thfermi(rhor))**2 + (three*qnrm**4)/(four*k_fermi(rhor)**2 * k_thfermi(rhor)**2) )

end function mdielf_bechstedt
!!***

!----------------------------------------------------------------------

!!****f* m_screening/screen_mdielf
!! NAME
!!  screen_mdielf
!!
!! FUNCTION
!!  Calculates W_{G,G'}(q,w) for a given q-point in the BZ using a model dielectric function.
!!
!! INPUTS
!!  iq_bz=The index of the q-point in the BZ where W(q) is calculated.
!!  npw=Number of plane waves for W
!!  nomega=Number of frequency points.
!!  model_type=Flag defining the model.
!!  eps_inf=Dielectric constant of the material.
!!  Cryst<crystal_t>=Info on the unit cell
!!  Qmesh<kmesh_t>=Info on the set of q-points.
!!  Vcp<vcoul_t datatype>= containing information on the cutoff technique
!!  Gsph<Gsphere>=The G-sphere for W.
!!  nspden=Number of spin density components of the density.
!!  nfft=Number of FFT points on the dense FFT mesh
!!  ngfft(18)=contain all needed information about 3D FFT.
!!  rhor(nfft,nspden)=Electron density in real space (The PAW AE term is included)
!!  which= Set it to "EM1" if the symmetrized inverse dielectric matrix is wanted.
!!   By default the routines returns W.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  w_qbz(npw,npw,nomega)
!!
!! NOTES
!!   W_{G1,G2} =  1/2 {
!!     v(q+G1) \int em1(|q+G1|,r) e^{-i(G1-G2).r} dr  +
!!     v(q+G2) \int em1(|q+G2|,r) e^{-i(G1-G2).r} dr } / \Omega
!!
!! PARENTS
!!      m_screen
!!
!! CHILDREN
!!
!! SOURCE

subroutine screen_mdielf(iq_bz,npw,nomega,model_type,eps_inf,Cryst,Qmesh,Vcp,Gsph,nspden,nfft,ngfft,rhor,which,w_qbz,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nomega,nfft,nspden,iq_bz,comm,model_type
 real(dp),intent(in) :: eps_inf
 character(len=*),intent(in) :: which
 type(kmesh_t),intent(in) :: Qmesh
 type(crystal_t),intent(in) :: Cryst
 type(vcoul_t),target,intent(in) :: Vcp
 type(gsphere_t),intent(in) :: Gsph
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,nspden)
 complex(gwpc),intent(out) :: w_qbz(npw,npw,nomega)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourdp0=0,paral_kgb0=0,cplex1=1
 integer :: my_gstart,my_gstop,iq_ibz,ig,itim_q,isym_q
 integer :: ig1,ig2,g1mg2_fft,iw,ii,ierr,nprocs,isg,ifft !,row ,col
 real(dp) :: qpg2_nrm
 complex(dpc) :: ph_mqbzt
 logical :: is_qeq0,isirred
 !character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: umklp(3)
 integer,allocatable :: igfft(:),g1mg2(:,:)
 real(dp) :: qpg2(3),qpt_bz(3)
 real(dp),allocatable :: em1_qpg2r(:),fofg(:,:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: vc_sqrt_ibz(:)
 complex(gwpc),allocatable :: vc_qbz(:),ctmp(:,:)
 logical,allocatable :: mask(:)

! *************************************************************************

 ABI_CHECK(nomega==1,"screen_mdielf does not support nomega>1")

!Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft(2),ngfft(3),'all')

 nprocs = xmpi_comm_size(comm)
 call xmpi_split_work(npw,comm,my_gstart,my_gstop)

 call get_bz_item(Qmesh,iq_bz,qpt_bz,iq_ibz,isym_q,itim_q,ph_mqbzt,umklp,isirred)

 !if (itim_q/=1.or.isym_q/=1.or.ANY(umklp/=0) ) then
 !  MSG_ERROR("Bug in mdielf_bechstedt")
 !end if
 !
 ! Symmetrize Vc in the full BZ.
 is_qeq0 = (normv(qpt_bz,Cryst%gmet,'G')<GW_TOLQ0) ! Check if q==0
 if (is_qeq0) then
   vc_sqrt_ibz => Vcp%vcqlwl_sqrt(:,1)  ! Use Coulomb term for q-->0, only first Q is used, shall we average if nqlwl>1?
 else
   vc_sqrt_ibz => Vcp%vc_sqrt(:,iq_ibz)
 end if

 ABI_MALLOC(vc_qbz,(npw))
 do ig=1,npw
   isg = Gsph%rottb(ig,itim_q,isym_q)
   vc_qbz(isg) = vc_sqrt_ibz(ig)**2
 end do

 ABI_MALLOC(igfft,(npw))
 ABI_MALLOC(g1mg2,(3,npw))
 ABI_MALLOC(fofg,(2,nfft))
 ABI_MALLOC(em1_qpg2r,(nfft))
 ABI_MALLOC(mask,(npw))

 w_qbz=czero
 do ig2=my_gstart,my_gstop
   !
   ! Compute the index of G-G2 wave in the FFT grid.
   do ii=1,npw
     g1mg2(:,ii) = Gsph%gvec(:,ii) - Gsph%gvec(:,ig2)
   end do
   call kgindex(igfft,g1mg2,mask,MPI_enreg_seq,ngfft,npw)

   ! TODO can use zero-padding FFT to speed up the transform.
   !call sphereboundary(gbound,istwfk1,g1mg2,mgfft,npw)

   ! Evaluate em1_qpg2r = \int em1(|q+G2|,r) e^{-i(G1-G2).r} dr }.
   qpg2 = qpt_bz + Gsph%gvec(:,ig2)
   qpg2_nrm = normv(qpg2,Cryst%gmet,"G")

   do iw=1,nomega
     !
     select case (model_type)
     case (1)
       do ifft=1,nfft
         em1_qpg2r(ifft) = one / mdielf_bechstedt(eps_inf,qpg2_nrm,rhor(ifft,1))
       end do
     case default
       MSG_ERROR(sjoin("Unknown model_type:",itoa(model_type)))
     end select

     call fourdp(cplex1,fofg,em1_qpg2r,-1,MPI_enreg_seq,nfft,1,ngfft,tim_fourdp0)
     !
     ! Here, unlike the other parts of the code, the unsymmetrized e^{-1} is used.
     do ig1=1,npw
       g1mg2_fft = igfft(ig1)
       w_qbz(ig1,ig2,iw) = DCMPLX(fofg(1,g1mg2_fft), fofg(2,g1mg2_fft)) * vc_qbz(ig2) !/ Cryst%ucvol
     end do
   end do ! iw
 end do ! ig2

 ABI_FREE(em1_qpg2r)
 ABI_FREE(fofg)
 ABI_FREE(igfft)
 ABI_FREE(g1mg2)
 ABI_FREE(mask)
 !
 ! W = 1/2 * (A + A^H)
 ! The MPI sum is done inside the loop to avoid problems with the size of the packet.
 ABI_MALLOC_OR_DIE(ctmp,(npw,npw), ierr)

 do iw=1,nomega
   !ctmp = TRANSPOSE(CONJG(w_qbz(:,:,iw)))
   ctmp = GWPC_CONJG(w_qbz(:,:,iw))
   call sqmat_itranspose(npw,ctmp)
   w_qbz(:,:,iw) = half * (ctmp + w_qbz(:,:,iw))
   call xmpi_sum(w_qbz(:,:,iw),comm,ierr)
 end do
 !
 ! Calculate the symmetrized Em1. W = vc(G1)^{1/2} \tilde Em1 vc(G2)^{1/2} -------------------------
 if (toupper(which)=="EM1") then
   do ig=1,npw
     isg = Gsph%rottb(ig,itim_q,isym_q)
     vc_qbz(isg) = vc_sqrt_ibz(ig)  ! Workspace storing vc*{1/2}(q_BZ,G).
   end do

   do ig2=1,npw
     do ig1=1,npw
       ctmp(ig1,ig2) =  one / (vc_qbz(ig1) * vc_qbz(ig2))
     end do
   end do

   do iw=1,nomega
     w_qbz(:,:,iw) = w_qbz(:,:,iw) * ctmp(:,:)
   end do
 end if

 call destroy_mpi_enreg(MPI_enreg_seq)

 ABI_FREE(vc_qbz)
 if (allocated(ctmp)) then
   ABI_FREE(ctmp)
 end if

end subroutine screen_mdielf
!!***

!----------------------------------------------------------------------

!!****f* m_screening/chi_new
!! NAME
!! chi_new
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(chi_t) function chi_new(npwe, nomega) result(chi)

!Arguments ------------------------------------
 integer,intent(in) :: npwe,nomega

! *************************************************************************

 chi%nomega = nomega; chi%npwe = npwe

 ABI_MALLOC(chi%mat, (npwe,npwe,nomega))

 ABI_MALLOC(chi%head, (3,3,nomega))
 ABI_MALLOC(chi%lwing, (npwe, nomega,3))
 ABI_MALLOC(chi%uwing, (npwe, nomega,3))

end function chi_new
!!***

!----------------------------------------------------------------------

!!****f* m_screening/chi_free
!! NAME
!! chi_free
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine chi_free(chi)

!Arguments ------------------------------------
!scalars
 type(chi_t),intent(inout) :: chi

! *************************************************************************

 if (allocated(chi%mat)) then
   ABI_FREE(chi%mat)
 end if
 if (allocated(chi%head)) then
   ABI_FREE(chi%head)
 end if
 if (allocated(chi%lwing)) then
   ABI_FREE(chi%lwing)
 end if
 if (allocated(chi%uwing)) then
   ABI_FREE(chi%uwing)
 end if

end subroutine chi_free
!!***

!----------------------------------------------------------------------

!!****f* m_screening/lwl_write
!! NAME
!! lwl_write
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine lwl_write(path, cryst, vcp, npwe, nomega, gvec, chi0, chi0_head, chi0_lwing, chi0_uwing, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwe,nomega,comm
 character(len=*),intent(in) :: path
 type(crystal_t),intent(in) :: cryst
 type(vcoul_t),intent(in) :: Vcp
!arrays
 integer,intent(in) :: gvec(3,npwe)
 complex(gwpc),intent(in) :: chi0(npwe,npwe,nomega)
 complex(dpc),intent(inout)  :: chi0_head(3,3,nomega),chi0_lwing(npwe,nomega,3),chi0_uwing(npwe,nomega,3)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,prtvol=0
 integer :: iw,ii,iomode,unt,my_rank
 character(len=500) :: msg
 real(dp) :: length
 complex(dpc) :: wng(3),em1_00
! type(hscr_t),intent(out) :: hscr
!arrays
 complex(dpc),allocatable :: wtest(:),eps_head(:,:,:)

! *************************************************************************

 !if (xmpi_comm_rank(comm) /= master) goto 100
 my_rank = xmpi_comm_rank(comm)

 ABI_MALLOC(wtest,(npwe))

 if (prtvol > 0 .and. my_rank == master) then
   call wrtout(std_out, "head of chi0")
   do iw=1,nomega
     call print_arr(chi0_head(:,:,iw),max_r=3,max_c=3)
   end do

   do iw=1,nomega
     call wrtout(std_out, "symmetrized e_00 via tensor")
     wng = MATMUL(chi0_head(:,:,iw),GW_Q0_DEFAULT)
     write(std_out,*) one - vdotw(GW_Q0_DEFAULT,wng,cryst%gmet,"G") * Vcp%vcqlwl_sqrt(1,1)*Vcp%vcqlwl_sqrt(1,1)

     call wrtout(std_out, "symmetrized e_0G via tensor")
     do ii=1,npwe
       wng = chi0_uwing(ii,iw,:)
       wtest(ii) = - vdotw(GW_Q0_DEFAULT,wng,cryst%gmet,"G") * Vcp%vcqlwl_sqrt(1,1) * Vcp%vcqlwl_sqrt(ii,1)
     end do
     call print_arr(wtest,max_r=9,unit=std_out)

     call wrtout(std_out, "symmetrized e_G0 via tensor")
     do ii=1,npwe
       wng = chi0_lwing(ii,iw,:)
       wtest(ii) = - vdotw(GW_Q0_DEFAULT,wng,cryst%gmet,"G") * Vcp%vcqlwl_sqrt(1,1) * Vcp%vcqlwl_sqrt(ii,1)
     end do
     call print_arr(wtest,max_r=9,unit=std_out)
   end do
 end if

 ! Write chi0 data
 iomode = IO_MODE_FORTRAN; if (endswith(path, ".nc")) iomode = IO_MODE_ETSF

 if (my_rank == master) then
   if (iomode == IO_MODE_FORTRAN) then
     if (open_file(path,msg,newunit=unt,form="unformatted", action="write") /= 0) then
       MSG_ERROR(msg)
     end if
     !call hscr_io(er%hscr,fform,rdwr,unt,comm,master,iomode)
     do iw=1,nomega
       write(unt)chi0_head(:,:,iw)
     end do
     !do iw=1,nomega
     !  write(unt)chi0_lwing(:,iw,:)
     !end do
     !do iw=1,nomega
     !  write(unt)chi0_uwing(:,iw,:)
     !end do

   else
     MSG_ERROR(sjoin("iomode", itoa(iomode), "is not supported"))
   end if
 end if

 ABI_MALLOC(eps_head,(3,3,nomega))
 call mkem1_q0(npwe,1,1,nomega,cryst,Vcp,gvec,chi0_head,chi0_lwing,chi0_uwing,chi0,eps_head,comm)

 if (my_rank == master) then
   if (iomode == IO_MODE_FORTRAN) then
     do iw=1,nomega
       write(unt)eps_head(:,:,iw)
     end do
     !do iw=1,nomega
     !  write(unt)chi0_lwing(:,iw,:)
     !end do
     !do iw=1,nomega
     !  write(unt)chi0_uwing(:,iw,:)
     !end do
   else
     MSG_ERROR(sjoin("iomode:", itoa(iomode), "is not supported"))
   end if
 end if

 ABI_FREE(eps_head)

 if (prtvol > 0 .and. my_rank == master) then
   length = normv(GW_Q0_DEFAULT,cryst%gmet,"G")

   do iw=1,nomega
     em1_00 = one / vdotw(GW_Q0_DEFAULT/length, MATMUL(chi0_head(:,:,iw),GW_Q0_DEFAULT/length),cryst%gmet,"G")
     call wrtout(std_out, "e^1_{00} from tensor")
     write(std_out,*) em1_00

     call wrtout(std_out, "symmetrized e^-1_0G via tensor")
     do ii=1,npwe
       wng = chi0_uwing(ii,iw,:)
       wtest(ii) = em1_00*vdotw(GW_Q0_DEFAULT/length,wng,cryst%gmet,"G")
     end do
     wtest(1) = em1_00
     call print_arr(wtest,max_r=9,unit=std_out)

     call wrtout(std_out, "symmetrized e^-1_G0 via tensor")
     do ii=1,npwe
       wng = chi0_lwing(ii,iw,:)
       wtest(ii) = em1_00*vdotw(GW_Q0_DEFAULT/length,wng,cryst%gmet,"G")
     end do
     wtest(1) = em1_00
     call print_arr(wtest,max_r=9,unit=std_out)
   end do !iw
 end if

 ABI_FREE(wtest)

 if (my_rank == master) then
   if (iomode == IO_MODE_FORTRAN) then
     close(unt)
   else
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_close(unt))
#endif
   end if
 end if

!100 call xmpi_barrier(comm)

end subroutine lwl_write
!!***

!----------------------------------------------------------------------

!!****f* m_screening/lwl_init
!! NAME
!! lwl_init
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine lwl_init(lwl, path, method, cryst, vcp, npwe, gvec, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,method,npwe
 character(len=*),intent(in) :: path
 type(crystal_t),intent(in) :: cryst
 type(vcoul_t),intent(in) :: vcp
 type(lwl_t),intent(out) :: lwl
!arrays
 integer,intent(in) :: gvec(3,npwe)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iomode,my_rank,nproc,unt
 character(len=500) :: msg
!arrays

! *************************************************************************

 ABI_UNUSED((/cryst%natom, gvec(1,1), vcp%ng/))

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 lwl%fname = path
 lwl%method = method
 ABI_CHECK(any(method == [1,2,3]), sjoin("Wrong method:", itoa(method)))

 iomode = IO_MODE_FORTRAN; if (endswith(path, ".nc")) iomode = IO_MODE_ETSF

 ! Only master reads.
 if (my_rank == master) then

   select case (iomode)
   case (IO_MODE_FORTRAN)
     if (open_file(path, msg, newunit=unt, action="read", form="unformatted", status="old") /= 0) then
       MSG_ERROR(msg)
     end if

     close(unt)

   case default
     MSG_ERROR(sjoin("iomode:", itoa(iomode), "is not coded"))
   end select
 end if

 ! Broad cast data
 if (nproc > 1) then
 end if

end subroutine lwl_init
!!***

!----------------------------------------------------------------------

!!****f* m_screening/lwl_free
!! NAME
!! lwl_free
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine lwl_free(lwl)

!Arguments ------------------------------------
!scalars
 type(lwl_t),intent(inout) :: lwl

! *************************************************************************

 if (allocated(lwl%head)) then
   ABI_FREE(lwl%head)
 end if
 if (allocated(lwl%lwing)) then
   ABI_FREE(lwl%lwing)
 end if
 if (allocated(lwl%uwing)) then
   ABI_FREE(lwl%uwing)
 end if
 if (allocated(lwl%body)) then
   ABI_FREE(lwl%body)
 end if

end subroutine lwl_free
!!***

END MODULE m_screening
!!***
