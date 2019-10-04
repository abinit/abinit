!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pawpwij
!! NAME
!!  m_pawpwij
!!
!! FUNCTION
!!  This module defines methods to calculate the onsite contribution of a plane wave in the PAW method:
!!
!!    - pawpwff_t: Form factors used to calculate the onsite contributions of a plane wave.
!!    - pawpwij_t: Onsite matrix elements of a plane wave for a given atom type.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (MG,GKA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_pawpwij

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_abicore
 use m_fft

 use m_numeric_tools,  only : arth
 use m_geometry,       only : metric
 use m_paw_numeric,    only : paw_jbessel_4spline, paw_spline
 use m_splines,        only : splfit
 use m_pawang,         only : pawang_type
 use m_paw_sphharm,     only : realgaunt
 use m_pawrad,         only : pawrad_type, pawrad_init, pawrad_free, pawrad_copy, simp_gen
 use m_pawtab,         only : pawtab_type
 use m_pawcprj,        only : pawcprj_type
 use m_mpinfo,         only : destroy_mpi_enreg, initmpi_seq
 use m_initylmg,       only : initylmg

 implicit none

 private
!!***

!!****t* m_pawpwij/pawpwff_t
!! NAME
!! pawpwff_t
!!
!! FUNCTION
!! For PAW, form factors used to evaluate $<phi|e^{-i(q+G).r}|phj> - <tphi|e^{-i(q+G).r}|tphj>$
!!
!! SOURCE

 type,public :: pawpwff_t

  integer :: method
  ! 1 For Arnaud-Alouani"s exact expression.
  ! 2 For Shishkin-Kresse"s approximated expression.

  integer :: dim1
  integer :: dim2
  ! Dimensions of pwff_spl, depending on method.

  integer :: nq_spl
   ! Number of points in the reciprocal space grid on which
   ! the radial integrals are evaluated.

  real(dp) :: gmet(3,3)
  ! Reciprocal space metric tensor in Bohr**-2

  real(dp),allocatable :: qgrid_spl(:)
   ! qgrid_spl(nq_spl)
   ! The coordinates of the points of the radial grid for the integrals used in the spline.

  real(dp),allocatable :: pwff_spl(:,:,:,:)
  ! pwff_spl(nq_spl,2,0:dim1,dim2)
  ! The different integrals on the radial |q| grid, for a given atom type.

 end type pawpwff_t

 public :: pawpwff_init    ! Initialize form factors for spline.
 public :: pawpwff_free    ! Deallocate dynamic memory.
!!***

!----------------------------------------------------------------------

!!****t* m_pawpwij/pawpwij_t
!! NAME
!! pawpwij_t
!!
!! FUNCTION
!! For PAW, object storing $<phi|e^{-i(q+G).r}|phj> - <tphi|e^{-i(q+G).r}|tphj>$
!! for a given q-point, for a particular TYPE of atom. Therefore the phase factor
!! e^{i(q+G).R_at} has to be considered to have the onsite contribution of a particular atom.
!!
!! SOURCE

 type,public :: pawpwij_t

  integer :: istpw
  ! Storage mode (similar to istwfk), not used at present

  integer :: npw
   ! The number of plane waves

  integer :: lmn_size
   ! Number of (l,m,n) elements for the paw basis

  integer :: lmn2_size
   ! lmn2_size=lmn_size*(lmn_size+1)/2
   ! where lmn_size is the number of (l,m,n) elements for the paw basis

  real(dp) :: qpt(3)
  ! The q-point in e^{-i(q+G)}.r}

  real(dp),allocatable :: mqpgij(:,:,:)
  ! pwij(2,npw,lmn2_size)
  ! $<phi|e^{-i(q+G).r}|phj> - <tphi|e^{-i(q+G).r}|tphj>$

 end type pawpwij_t

 public :: pawpwij_init       ! Calculate onsite matrix elements of a set of plane waves.
 public :: pawpwij_free        ! Deallocate dynamic memory in the structure.
 !public :: paw_pwij_bcast
 public :: paw_rho_tw_g        ! Calculate the PAW contribution to the oscillator matrix element.
 public :: paw_cross_rho_tw_g  ! Calculate the PAW cross term contribution to the oscillator matrix element.
!!***

 interface pawpwij_free
   module procedure pawpwij_free_d1
   module procedure pawpwij_free_d2
 end interface pawpwij_free

!----------------------------------------------------------------------

 integer,parameter :: PWIJ_ARNAUD   = 1   ! Arnaud-Alouani exact expression. PRB 62. 4464 [[cite:Arnaud2000]]
 integer,parameter :: PWIJ_SHISHKIN = 2   ! Shishkin-Kresse approximated expression. PRB 74. 035101 [[cite:Shishkin2006]]

CONTAINS  !========================================================================================
!!***

!!****f* m_pawpwij/pawpwff_init
!! NAME
!!  pawpwff_init
!!
!! FUNCTION
!!  Initialize the structure containing integrals used to evaluate the onsite
!!  matrix elements of a plane wave by means of a spline fit technique.
!!
!! INPUTS
!!  method=1 for Arnaud-Alouani, 2 for Shishkin-Kresse.
!!  nq_spl(%ntypat)=Number of points in the mesh used for the spline.
!!  qmax(%ntypat)=Max |q| for the mesh
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  Pawrad(%ntypat) <type(pawrad_type)>=paw radial mesh and related data....
!!  Pawtab(%ntypat) <type(pawtab_type)>=paw tabulated starting data.
!!    %lsize=1+maximum value of l leading to non zero Gaunt coeffs.
!!    %ij_size=Number of (i,j) elements for the symetric paw basis
!!    %lmn2_size=lmn_size*(lmn_size+1)/2
!!  Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!    %ntypat=Number of type of atoms.
!!
!! OUTPUT
!!  Paw_pwff(%ntypat) <pawpwff_t>=Object storing the form factors
!!                                    for the spline used in pawpwij_init.
!! PARENTS
!!      bethe_salpeter,m_shirley,screening,sigma
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many,fftpad
!!
!! SOURCE

subroutine pawpwff_init(Paw_pwff,method,nq_spl,qmax,gmet,Pawrad,Pawtab,Psps)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: method
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: nq_spl(Psps%ntypat)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: qmax(Psps%ntypat)
 type(Pawrad_type),intent(in) :: Pawrad(Psps%ntypat)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(pawpwff_t),intent(out) :: Paw_pwff(Psps%ntypat)

!Local variables-------------------------------
!scalars
 integer :: dim1,dim2,itypat,nq
 real(dp) :: dq

!************************************************************************

 !@pawpwff_t

! === Evaluate form factors for the radial part of phi.phj-tphi.tphj ===
 do itypat=1,Psps%ntypat
   Paw_pwff(itypat)%method = method

   select case (method)
   case (PWIJ_ARNAUD)
     dim1 = Pawtab(itypat)%l_size-1
     dim2 = Pawtab(itypat)%ij_size
   case (PWIJ_SHISHKIN)
     dim1 = Pawtab(itypat)%l_size**2
     dim2 = Pawtab(itypat)%lmn2_size
   case default
     MSG_BUG("Wrong method")
   end select

   Paw_pwff(itypat)%dim1 = dim1
   Paw_pwff(itypat)%dim2 = dim2

   Paw_pwff(itypat)%gmet = gmet

   ! === Setup of the q-mesh for spline ===
   ! * It can be type-dependent.
   nq = nq_spl(itypat)
   dq = qmax(itypat)/(one*(nq-1))
   !write(std_out,*)"nq,dq",nq,dq

   Paw_pwff(itypat)%nq_spl = nq
   ABI_MALLOC(Paw_pwff(itypat)%qgrid_spl,(nq))
   Paw_pwff(itypat)%qgrid_spl = arth(zero,dq,nq)
   !
   ! === Calculate form factors depending on method ===
   ABI_MALLOC(Paw_pwff(itypat)%pwff_spl,(nq,2,0:dim1,dim2))

   call paw_mkrhox_spl(itypat,Psps%ntypat,method,dim1,dim2,nq,&
&    Paw_pwff(itypat)%qgrid_spl,Pawrad,Pawtab,Paw_pwff(itypat)%pwff_spl)
 end do !itypat

end subroutine pawpwff_init
!!***

!----------------------------------------------------------------------

!!****f* m_pawpwij/pawpwff_free
!! NAME
!!  pawpwff_free
!!
!! FUNCTION
!!  Free the memory allocated in a structure of type pawpwff_t
!!
!! SIDE EFFECTS
!!  Paw_pwff(:)=<pawpwff_t>=Object storing form factors for the spline of wf into PAW spheres
!!
!! PARENTS
!!      bethe_salpeter,m_shirley,screening,sigma
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many,fftpad
!!
!! SOURCE

subroutine pawpwff_free(Paw_pwff)

 implicit none

!Arguments ------------------------------------
!scalars
 type(pawpwff_t),intent(inout) :: Paw_pwff(:)

!Local variables-------------------------------
!scalars
 integer :: ii

!************************************************************************

 !@pawpwff_t
 do ii=1,SIZE(Paw_pwff)
   if (allocated(Paw_pwff(ii)%qgrid_spl)) then
     ABI_FREE(Paw_pwff(ii)%qgrid_spl)
   end if
   if (allocated(Paw_pwff(ii)%pwff_spl )) then
     ABI_FREE(Paw_pwff(ii)%pwff_spl)
   end if
 end do

end subroutine pawpwff_free
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhox/pawpwij_init
!! NAME
!!  pawpwij_init
!!
!! FUNCTION
!!  Creation method for pawpwij_t. Calculates the onsite matrix elements
!!   $ <phj|e^{-i(q+G)}|phi> - <tphj|e^{-i(q+G)}|tphi> $
!!  for a given q and a set of G"s for a given __TYPE__ of atom.
!!  Phase factors arising from atom positions are therefore not included.
!!
!! INPUTS
!!  npw=Number of plane waves
!!  Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!    %ntypat=Number of type of atoms,
!   gvec(3,npw)=Plane wave reduced components.
!!  qpt_in(3)=The reduced components of the q-point.
!!  rprim(3,3)=dimensionless real space primitive translations
!!  Pawtab(%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  Paw_pwff(%ntypat) <pawpwff_t>=Object storing the form factors for the spline used in pawpwij_init.
!!  Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!    %ntypat
!!
!! OUTPUT
!!  Pwij(%ntypat)<pawpwij_t>=Structure containing the onsite matrix elements of e^{-i(q+G).r}.
!!   Completely initialized in output.
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cchi0,cchi0q0,cchi0q0_intraband,cohsex_me
!!      exc_build_block,exc_build_ham,m_shirley,prep_calc_ucrpa
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many,fftpad
!!
!! SOURCE

subroutine pawpwij_init(Pwij,npw,qpt_in,gvec,rprimd,Psps,Pawtab,Paw_pwff)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: gvec(3,npw)
 real(dp),intent(in) :: qpt_in(3),rprimd(3,3)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat)
 type(pawpwij_t),intent(out) :: Pwij(Psps%ntypat)

!Local variables-------------------------------
!scalars
 integer,parameter :: unkg0=0,unylm0=0
 integer :: dim1,dim2,method
 integer :: my_mqmem,my_nqpt,optder,two_lmaxp1,itypat
 integer :: dummy_nsppol,lmn_size,lmn2_size,nq_spl,ierr
 real(dp) :: ucvol
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer,allocatable :: npwarr(:),dummy_nband(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: my_qtmp(:,:)
 real(dp),allocatable :: ylm_q(:,:),ylmgr_q(:,:,:)

! *********************************************************************

 !@pawpwij_t

 ! ===============================================
 ! === Get real spherical harmonics in G space ===
 ! ===============================================

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ! * Set up of REAL Ylm(q+G) up to 2*l_max for this q-point.
 my_mqmem=1; two_lmaxp1=2*Psps%mpsang-1; optder=0

 ABI_MALLOC(ylm_q  ,(npw*my_mqmem,two_lmaxp1**2))
 ABI_MALLOC(ylmgr_q,(npw*my_mqmem,3+6*(optder/2),two_lmaxp1**2))

 my_nqpt=1
 ABI_MALLOC(my_qtmp,(3,my_nqpt))
 my_qtmp(:,1)=qpt_in(:)
 !
 ! dummy_nband and dummy_nsppol are not used in sequential mode.
 dummy_nsppol=1
 ABI_MALLOC(dummy_nband,(my_nqpt*dummy_nsppol))
 dummy_nband=0
 ABI_MALLOC(npwarr,(my_nqpt))
 npwarr(:)=npw

 ! Fake MPI_type for sequential part.
 call initmpi_seq(MPI_enreg_seq)

 call initylmg(gprimd,gvec,my_qtmp,my_mqmem,MPI_enreg_seq,two_lmaxp1,npw,dummy_nband,my_nqpt,&
   npwarr,dummy_nsppol,optder,rprimd,ylm_q,ylmgr_q)

 call destroy_mpi_enreg(MPI_enreg_seq)

 ABI_FREE(my_qtmp)
 ABI_FREE(dummy_nband)
 ABI_FREE(npwarr)
 ABI_FREE(ylmgr_q)

 ! Construct the Pwij structure
 do itypat=1,Psps%ntypat

   Pwij(itypat)%istpw = 1
   Pwij(itypat)%npw   = npw

   lmn_size  = Pawtab(itypat)%lmn_size
   lmn2_size = lmn_size*(lmn_size+1)/2
   Pwij(itypat)%lmn_size  = lmn_size
   Pwij(itypat)%lmn2_size = lmn2_size

   Pwij(itypat)%qpt(:) = qpt_in(:)

   ! Prepare the call to paw_mkrhox
   method = Paw_pwff(itypat)%method
   dim1   = Paw_pwff(itypat)%dim1
   dim2   = Paw_pwff(itypat)%dim2
   nq_spl = Paw_pwff(itypat)%nq_spl
   gmet   = Paw_pwff(itypat)%gmet

   ABI_MALLOC_OR_DIE(Pwij(itypat)%mqpgij,(2,npw,lmn2_size), ierr)

   ! Evaluate oscillator matrix elements mqpgij
   call paw_mkrhox(itypat,lmn2_size,method,dim1,dim2,nq_spl,Paw_pwff(itypat)%qgrid_spl,Paw_pwff(itypat)%pwff_spl,&
&    gmet,qpt_in,npw,gvec,ylm_q,Psps,Pawtab,Pwij(itypat)%mqpgij)
 end do !itypat

 ABI_FREE(ylm_q)

end subroutine pawpwij_init
!!***

!----------------------------------------------------------------------

!!****f* m_pawpwij/pawpwij_free_d1
!! NAME
!!  pawpwij_free_d1
!!
!! FUNCTION
!!  Free all memory allocated in a structure of type pawpwij_t
!!
!! SIDE EFFECTS
!!  Paw_pwij(:)=<pawpwij_t>=Structure containing the onsite matrix elements of e^{-i(q+G).r}
!!
!! PARENTS
!!      m_pawpwij
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many,fftpad
!!
!! SOURCE

subroutine pawpwij_free_d1(Pwij)

 implicit none

!Arguments ------------------------------------
!scalars
 type(pawpwij_t),intent(inout) :: Pwij(:)

!Local variables-------------------------------
!scalars
 integer :: ii

!************************************************************************

 !@pawpwij_t
 do ii=1,SIZE(Pwij)
   if (allocated(Pwij(ii)%mqpgij))  then
     ABI_FREE(Pwij(ii)%mqpgij)
   end if
 end do

end subroutine pawpwij_free_d1
!!***

!----------------------------------------------------------------------

!!****f* m_pawpwij/pawpwij_free_d2
!! NAME
!!  pawpwij_free_d2
!!
!! FUNCTION
!!  Free all memory allocated in a structure of type pawpwij_t
!!
!! SIDE EFFECTS
!!  Paw_pwij(:)=<pawpwij_t>=Structure containing the onsite matrix elements of e^{-i(q+G).r}
!!
!! PARENTS
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many,fftpad
!!
!! SOURCE

subroutine pawpwij_free_d2(Pwij)

 implicit none

!Arguments ------------------------------------
!scalars
 type(pawpwij_t),intent(inout) :: Pwij(:,:)

!Local variables-------------------------------
!scalars
 integer :: jj

!************************************************************************

 do jj=1,SIZE(Pwij,DIM=2)
   call pawpwij_free_d1(Pwij(:,jj))
 end do

end subroutine pawpwij_free_d2
!!***

!----------------------------------------------------------------------

!!****f* m_pawpwij/paw_mkrhox_spl
!! NAME
!! paw_mkrhox_spl
!!
!! FUNCTION
!!  Evaluate PAW form factor ff^{aL}_{ij}(q) for each angular momentum L,
!!  each type of atom, a, and each [(in,il),(jn,jl)] channel. These quantities
!!  are used in paw_mkrhox to evaluate $<phi|e^{-i(q+G)}|phj>-<tphi|e^{-i(q+G)}|tphj>$
!!  for an arbitrary q+G vector.
!!
!! INPUTS
!!  dim1
!!   = 2*(Pawtab(itypat)%l_size-1) if method=1
!!   = MAXVAL(Pawtab(:)%l_size)**2 if method=2
!!  dim2
!!   = Pawtab(itypat)%ij_size      if method=1
!!   = MAXVAL(Pawtab(:)%lmn2_size) if method=2
!!  method=integer flag defining the approach used:
!!   1 --> Expression based on the expansion on a plane wave in terms of Bessel functions
!!        and spherical harmonics (Arnaud-Alouani's methos, see PRB 62, 4464 [[cite:Arnaud2000]]
!!   2 --> Approximate expression with correct description of the multipoles. Eq. 9 in PRB 74, 035101 [[cite:Shishkin2006]]
!!  nq_spl=number of grid points in the q-mesh
!!  qgrid_spl(nq_spl)=values where form factors are returned
!!  ntypat=number of type of atoms
!!  Pawrad<type(Pawrad_type)>=datatype containing radial grid information
!!  Pawtab(ntypat)<type(pawtab_type)>=PAW tabulated starting data
!!
!! OUTPUT
!!  pwff_spl(nq_spl,2,0:dim1,dim1_rhox2,ntypat)
!!   form factors ff^{aL}_{ij}(q) and second derivative in packed storage mode
!!  === if method=1 ===
!!    $ff_^{aL}_{ij}(q) =
!!      \int_0^{r_a} j_L(2\pi qr) [phi_{n_i,l_i}.phi_{n_j l_j}(r) - tphi_{n_i l_i}.tph_{n_j l_j}(r)] dr$
!!  === if method=2 ===
!!    $ff_^{aL}_{ij}(q) = q_ij^{LM} \int_0^{r_a} j_L(2\pi qr) g_L(r) r^2 dr$
!!
!! NOTES
!!  * $j_L(2\pi q)$ is a spherical Bessel function
!!  * Output matrix elements are stored in packed storage mode
!!  * Inspired by pawpsp_nl
!!
!! TODO
!!  One might save CPU time taking into account Gaunt selection rules!
!!
!! PARENTS
!!      m_pawpwij
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many,fftpad
!!
!! SOURCE

subroutine paw_mkrhox_spl(itypat,ntypat,method,dim1,dim2,nq_spl,qgrid_spl,Pawrad,Pawtab,pwff_spl)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itypat,ntypat,method,dim2,dim1,nq_spl
!arrays
 real(dp),intent(in) :: qgrid_spl(nq_spl)
 real(dp),intent(out) :: pwff_spl(nq_spl,2,0:dim1,dim2)
 type(Pawrad_type),intent(in) :: Pawrad(ntypat)
 type(Pawtab_type),intent(in) :: Pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: mm,nlmn,jlmn,ilmn,klmn,ij_size,l_size,ider
 integer :: iq,ir,ll,meshsz,mmax,iln,jln,nln,k0ln,kln,qlm
 real(dp),parameter :: EPS=tol14**4,TOLJ=0.001_dp
 real(dp) :: arg,argn,bes,besp,qr,yp1,ypn
 character(len=500) :: msg
 type(Pawrad_type) :: Tmpmesh
!arrays
 real(dp),allocatable :: rrshape_l(:),shape_l(:),ff(:),gg(:),rr(:),rrdphi_ij(:)
 real(dp),allocatable :: dphi_ij(:),tmp_spl(:,:,:,:),tmp_jgl(:,:,:)

!*************************************************************************

 DBG_ENTER("COLL")

 pwff_spl=zero

 SELECT CASE (method)

 CASE (PWIJ_ARNAUD)
   ! === Arnaud-Alouani exact expression PRB 62. 4464 [[cite:Arnaud2000]] ===
   ! * $ff_^{aL}_{ij}(q) =
   !    \int_0^{r_a} j_L(2\pi qr) [phi_{n_i,l_i}.phi_{n_j l_j}(r)-tphi_{n_i l_i}.tph_{n_j l_j}(r)]dr$
   ! * It does not descrive correctly the multipoles of the AE charge density if low cutoff on G
   write(msg,'(a,i3)')' paw_mkrhox_spl: Using Arnaud-Alouani expression for atom type: ',itypat
   call wrtout(std_out,msg,'COLL')

   nln     = Pawtab(itypat)%basis_size
   ij_size = Pawtab(itypat)%ij_size
   l_size  = Pawtab(itypat)%l_size

   ABI_MALLOC(tmp_spl,(nq_spl,2,0:l_size-1,ij_size))

   ! Is mesh beginning with r=0 ?
   if (ABS(Pawrad(itypat)%rad(1))>tol10) then
     MSG_ERROR("Radial mesh starts with r/=0")
   end if
   !
   ! === Initialize temporary arrays and variables ===
   meshsz = Pawtab(itypat)%mesh_size ; mmax=meshsz

   ABI_MALLOC(dphi_ij,(meshsz))
   ABI_MALLOC(rrdphi_ij,(meshsz))
   ABI_MALLOC(ff,(meshsz))
   ABI_CALLOC(gg,(meshsz))
   ABI_MALLOC(rr,(meshsz))
   rr(1:meshsz) = Pawrad(itypat)%rad(1:meshsz)
   !
   ! === Loop on (jln,iln) channels for this type. Packed form ===
   tmp_spl=zero

   do jln=1,nln
     k0ln=jln*(jln-1)/2
     do iln=1,jln
       kln=k0ln+iln

       dphi_ij(1:meshsz) = Pawtab(itypat)%phiphj(1:meshsz,kln)-Pawtab(itypat)%tphitphj(1:meshsz,kln)
       rrdphi_ij(1:meshsz) = rr(1:meshsz)*dphi_ij(1:meshsz) ! will be used for first derivative

       ir=meshsz
       do while (ABS(dphi_ij(ir))<EPS)
         ir=ir-1
       end do
       mmax=MIN(ir+1,meshsz)
       ! ir is equal to meshsz if no point are below EPS
       if (mmax/=Pawrad(itypat)%int_meshsz) then  ! mmax=meshsz
         call pawrad_init(Tmpmesh,mesh_size=Pawtab(itypat)%mesh_size,mesh_type=Pawrad(itypat)%mesh_type, &
&                         rstep=Pawrad(itypat)%rstep,lstep=Pawrad(itypat)%lstep,r_for_intg=rr(mmax))
       else
         call pawrad_copy(Pawrad(itypat),Tmpmesh)
       end if
       !
       ! === Loop on l for Bessel function. Note the starting point ===
       ! TODO Here I should loop only the moments allowed by Gaunt, I should use indklm!
       ! and only on lmax for this atom
       do ll=0,l_size-1
         !
         ! === Compute f_l(q=0) only if l=0, and first derivative fp_l(q=0) (nonzero only if ll==1) ===
         tmp_spl(1,1,ll,kln)=zero; yp1=zero
         if (ll==0) then
           call simp_gen(tmp_spl(1,1,ll,kln),dphi_ij,Tmpmesh)
         end if
         if (ll==1) then
           call simp_gen(yp1,rrdphi_ij,Tmpmesh)
           yp1=yp1*two_pi*third
         end if
         !
         ! === Compute f_l(0<q<qmax) ===
         if (nq_spl>2) then
           do iq=2,nq_spl-1
             arg=two_pi*qgrid_spl(iq)
             do ir=1,mmax
               qr=arg*rr(ir)
               call paw_jbessel_4spline(bes,besp,ll,0,qr,TOLJ)
               ff(ir)=bes*dphi_ij(ir)
             end do
             call simp_gen(tmp_spl(iq,1,ll,kln),ff,Tmpmesh)
           end do
         end if
         !
         ! === Compute f_l(q=qmax) and first derivative ===
         if (nq_spl>1) then
           argn=two_pi*qgrid_spl(nq_spl)
           do ir=1,mmax
             qr=argn*rr(ir)
             call paw_jbessel_4spline(bes,besp,ll,1,qr,TOLJ)
             ff(ir)=bes *  dphi_ij(ir)
             gg(ir)=besp*rrdphi_ij(ir)
           end do
           call simp_gen(tmp_spl(nq_spl,1,ll,kln),ff,Tmpmesh)
           gg(:)=two_pi*gg(:) ! two_pi comes from 2\pi|q| in the Bessel function
           call simp_gen(ypn,gg,Tmpmesh)
         else
           ypn=yp1
         end if
         !
         ! === Compute second derivative of ff^{al}_{ij)(q) ===
         !yp1=zero; ypn=zero
         call paw_spline(qgrid_spl,tmp_spl(:,1,ll,kln),nq_spl,yp1,ypn,tmp_spl(:,2,ll,kln))
       end do !ll

       call pawrad_free(Tmpmesh)

     end do !iln
   end do !jln
   !
   ! === Save values for this atom type, each ll and kln channel ===
   pwff_spl = tmp_spl

   ABI_FREE(dphi_ij)
   ABI_FREE(rrdphi_ij)
   ABI_FREE(ff)
   ABI_FREE(gg)
   ABI_FREE(rr)
   ABI_FREE(tmp_spl)

   if (.FALSE.) then ! write form factors on file for plotting purpose.
     ll=0
     do iq=1,nq_spl
       write(777+itypat,'(50(es16.8))')qgrid_spl(iq),((pwff_spl(iq,ider,ll,iln),ider=1,2),iln=1,dim2)
     end do
   end if

 CASE (PWIJ_SHISHKIN)
   ! ==== Shishkin-Kresse approximated expression ====
   ! $ff_^{aL}_{ij}(q) = q_ij^{LM} \int_0^{r_a} j_L(2\pi qr) g_L(r) r^2 dr$
   ! * Better description of multipoles of AE charge,
   ! * Better results for energy degeneracies in GW band structure
   write(msg,'(a,i3)')' paw_mkrhox_spl: Using Shishkin-Kresse expression for atom type ',itypat
   call wrtout(std_out,msg,'COLL')
   l_size = Pawtab(itypat)%l_size
   nlmn   = Pawtab(itypat)%lmn_size

   !allocate(tmp_jgl(nq_spl,2,0:2*(Psps%mpsang-1)))
   ABI_MALLOC(tmp_jgl,(nq_spl,2,0:l_size-1))

   ! Is mesh beginning with r=0 ?
   if (ABS(Pawrad(itypat)%rad(1))>tol10) then
     MSG_ERROR("Radial mesh starts with r/=0")
   end if
   !
   ! === Initialize temporary arrays and variables ===

   call pawrad_copy(Pawrad(itypat),Tmpmesh)
   meshsz=Pawtab(itypat)%mesh_size ; mmax=meshsz
   ABI_MALLOC(ff,(meshsz))
   ABI_CALLOC(gg,(meshsz))
   ABI_MALLOC(rr,(meshsz))
   ABI_MALLOC(shape_l,(meshsz))
   ABI_MALLOC(rrshape_l,(meshsz))
   rr(1:meshsz)=Tmpmesh%rad(1:meshsz)

   tmp_jgl(:,:,:)=zero
   rrshape_l(1)=zero
     shape_l(1)=zero
   !
   ! TODO Here I should loop only the moments allowed by Gaunt, I should use indklm!
   ! and only lmax for this atom
   !do ll=0,2*(Psps%mpsang-1)
   do ll=0,l_size-1
       shape_l(2:meshsz)=Pawtab(itypat)%shapefunc(2:meshsz,ll+1)*rr(2:meshsz)**2
     rrshape_l(2:meshsz)=Pawtab(itypat)%shapefunc(2:meshsz,ll+1)*rr(2:meshsz)**3
     !
     ! === Compute f_l(q=0) and first derivative fp_l(q=0) (only if ll==1) ===
     tmp_jgl(1,1,ll)=zero ; yp1=zero
     if (ll==0) then
       call simp_gen(tmp_jgl(1,1,ll),shape_l,Tmpmesh)
     end if
     if (ll==1) then
       call simp_gen(yp1,rrshape_l,Tmpmesh) !rr comes from d/dq
       yp1=yp1*two_pi*third
     end if
     !
     ! === Compute f_l(0<q<qmax) ===
     if (nq_spl>2) then
       do iq=2,nq_spl-1
         arg=two_pi*qgrid_spl(iq)
         do ir=1,mmax
           qr=arg*rr(ir)
           call paw_jbessel_4spline(bes,besp,ll,0,qr,TOLJ)
           ff(ir)=bes*shape_l(ir)
         end do
         call simp_gen(tmp_jgl(iq,1,ll),ff,Tmpmesh)
       end do
     end if
     !
     ! === Compute f_l(q=qmax) and first derivative ===
     if (nq_spl>1) then
       argn=two_pi*qgrid_spl(nq_spl)
       do ir=1,mmax
         qr=argn*rr(ir)
         call paw_jbessel_4spline(bes,besp,ll,1,qr,TOLJ)
         ff(ir)=bes *  shape_l(ir)
         gg(ir)=besp*rrshape_l(ir)
       end do
       call simp_gen(tmp_jgl(nq_spl,1,ll),ff,Tmpmesh)
       gg(:)=two_pi*gg(:) !two_pi comes from 2\pi|q|
       call simp_gen(ypn,gg,Tmpmesh)
     else
       ypn=yp1
     end if
     !
     ! === Compute second derivative of ff_^{al}_{ij)(q) ===
     call paw_spline(qgrid_spl,tmp_jgl(:,1,ll),nq_spl,yp1,ypn,tmp_jgl(:,2,ll))
   end do !ll
   !
   ! === Save values for this type, each ll and ilmn,jlm channels ===
   ! * Here we assembly q_{ij}^{lm} \int_0^{r_a} j_l(2\pi(q+G)r) g_l(r)r^2 dr
   ! * Some of the contributions from qijl are zero due to Gaunt selection rules
   do jlmn=1,nlmn
     do ilmn=1,jlmn
       klmn=ilmn+(jlmn-1)*jlmn/2
       !do ll=0,2*(Psps%mpsang-1)
       do ll=0,l_size-1
         do mm=-ll,ll
           qlm=1+ll**2+ll+mm
           pwff_spl(:,:,qlm-1,klmn)=tmp_jgl(:,:,ll)*Pawtab(itypat)%qijl(qlm,klmn)
         end do
       end do
     end do
   end do

   call pawrad_free(Tmpmesh)
   ABI_FREE(shape_l)
   ABI_FREE(rrshape_l)
   ABI_FREE(ff)
   ABI_FREE(gg)
   ABI_FREE(rr)
   ABI_FREE(tmp_jgl)

 CASE DEFAULT
   write(msg,'(a,i3)')' Called with wrong value for method ',method
   MSG_BUG(msg)
 END SELECT

 DBG_EXIT("COLL")

end subroutine paw_mkrhox_spl
!!***

!----------------------------------------------------------------------

!!****f* m_pawpwij/paw_mkrhox
!! NAME
!! paw_mkrhox
!!
!! FUNCTION
!!  Evaluate $<phj|e^{-i(q+G)}|phi>-<tphj|e^{-i(q+G)}|tphi>$
!!  for a fixed q-point and npw G vectors. Matrix elements are stored in packed storage mode.
!!
!! INPUTS
!!  gmet(3,3)=reciprocal lattice metric tensor ($\textrm{Bohr}^{-2}$)
!!  gvec(3,npw)=G vectors in reduced coordinates
!!  npw=numper of G vectors
!!  Psps<pseudopotential_type>:
!!     %lnmax=Max. number of (l,n) components over all type of PAW datasets
!!  nq_spl=Number of points in the reciprocal space grid on which the radial functions pwff_spl are specified
!!  qgrid_spl(nq_spl)=values at which form factors have been evaluated
!!     %mpsang=1+maximum angular momentum
!!  qpt(3)= q-point in reduced coordinates
!!  ylm_q(npw,(2*Psps%mpsang-1)**2)=real spherical harmonics Ylm(q+G) for q-point qpt up to l=2*l_max
!!  pwff_spl(nq_spl,2,0:2*(Psps%mpsang-1),Psps%lnmax*(Psps%lnmax+1)/2))
!!  Pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!     %indklmn(8,lmn2_size)=array giving klm, kln, abs(il-jl) and (il+jl), ilm and jlm, ilmn and jlmn for each klmn=(ilmn,jlmn)
!!
!! OUTPUT
!!  paw_rhox(2,npw,lmn2_size): $<phj|e^{-i(q+G).r}|phi>-<tphj|e^{-i(q+G).r}|tphi>$ in packed form for
!!    a type itypat (phase factor arising from atom position is not included)
!!
!! PARENTS
!!      m_pawpwij
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many,fftpad
!!
!! SOURCE

subroutine paw_mkrhox(itypat,lmn2_size,method,dim1,dim2,nq_spl,qgrid_spl,pwff_spl,&
&  gmet,qpt,npw,gvec,ylm_q,Psps,Pawtab,paw_rhox)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itypat,dim1,dim2,method,npw,nq_spl,lmn2_size
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: gvec(3,npw)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: pwff_spl(nq_spl,2,0:dim1,dim2)
 real(dp),intent(in) :: qpt(3),ylm_q(npw,(2*Psps%mpsang-1)**2)
 real(dp),intent(out) :: paw_rhox(2,npw,lmn2_size)
 real(dp),intent(in) :: qgrid_spl(nq_spl)
 type(Pawtab_type),target,intent(in) :: Pawtab(Psps%ntypat)

!Local variables-------------------------------
!scalars
 integer :: ider,ig,ignt,il,ilm,ilm_G,ilmn,iln,im,ipow,jl,jlm
 integer :: jlmn,jln,jm,k0lm,k0lmn,k0ln,klm,klmn,kln,ll_G,mm_G,mpsang,ngnt
 real(dp) :: rgnt,dummy
 character(len=500) :: msg
!arrays
 integer,allocatable :: gntselect(:,:)
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp) :: mi_l(2,0:3),qpg(3)
 real(dp),allocatable :: derfun(:),newfun(:),qpg_norm(:),realgnt(:),wk_ffnl(:,:)

! *************************************************************************

 !write(std_out,*)itypat,dim1,dim2,method,npw,nq_spl,lmn2_size,Psps%mpsang

 mpsang=Psps%mpsang
 indlmn => Pawtab(itypat)%indlmn
 !
 ! === Pre-calculate (-i)^l ===
 mi_l(1,0)=one  ; mi_l(2,0)=zero
 mi_l(1,1)=zero ; mi_l(2,1)=-one
 mi_l(1,2)=-one ; mi_l(2,2)=zero
 mi_l(1,3)=zero ; mi_l(2,3)=one
 !
 ! === Calculate |q+G| ===
 ! * 2\pi is not included to be consistent with the spline.
 ABI_MALLOC(qpg_norm,(npw))
 do ig=1,npw
   qpg = qpt + gvec(:,ig)
   qpg_norm(ig)=SQRT(DOT_PRODUCT(qpg,MATMUL(gmet,qpg)))
 end do

 ! Check q-grid as %qgrid_spl must be large enoung.
 if (MAXVAL(qpg_norm)>MAXVAL(qgrid_spl)) then
   write(msg,'(3a,f8.4,a,f8.4,2a)')&
&    ' Function values are being requested outside range of data. ',ch10,&
&    ' Max qpg_norm = ',MAXVAL(qpg_norm),' Max qgrid_spl = ',MAXVAL(qgrid_spl),ch10,&
&    ' Increase ecut(wfn), check qrid_ff and gsqcut '
   MSG_ERROR(msg)
 end if

 ABI_MALLOC(wk_ffnl,(nq_spl,2))
 ABI_MALLOC(newfun,(npw))
 ABI_MALLOC(derfun,(npw))

 SELECT CASE (method)

 CASE (PWIJ_ARNAUD)
   ! === Arnaud-Alouani exact expression ===
   ! * It does not describe the multipoles of the AE charge density
   ! * $ 4\pi \sum_{LM} (-i)^l Y_M^L(q+G) G_{\li\mi\lj\mj}^{\LM} ff^{aL}_{ij}(|q+G|) $
   !   where f has been calculated in paw_mkrhox_spl
   !
   ! === Re-evaluate Gaunt coefficients, just to be on the safe side ===
   ! * Note that gntselect is in packed form, thanks to invariance under permutation.
   ! * Could use Pawang% but size of gntselect depends on pawxcdev!

   ABI_MALLOC(  realgnt,((2*mpsang-1)**2*(mpsang)**4))
   ABI_MALLOC(gntselect,((2*mpsang-1)**2, mpsang**2*(mpsang**2+1)/2))
   call realgaunt(mpsang,ngnt,gntselect,realgnt)
   paw_rhox=zero
   !
   ! === Loop on (jl,jm,jn) channels for this atom ===
   do jlmn=1,Pawtab(itypat)%lmn_size
     jl =indlmn(1,jlmn)
     jm =indlmn(2,jlmn)
     jlm=indlmn(4,jlmn)
     jln=indlmn(5,jlmn)

     k0lmn=jlmn*(jlmn-1)/2
     k0lm =jlm *(jlm -1)/2
     k0ln =jln *(jln -1)/2
     !
     ! === Loop on (il,im,in) channels; klmn is index for packed form ===
     do ilmn=1,jlmn
       il =indlmn(1,ilmn)
       im =indlmn(2,ilmn)
       ilm=indlmn(4,ilmn)
       iln=indlmn(5,ilmn)

       klmn=k0lmn+ilmn
       klm =k0lm +ilm
       kln =k0ln +iln
       !
       ! === Summing over allowed (L,M), taking into account Gaunt selection rules ===
       do ll_G=ABS(jl-il),jl+il,2
         ipow=MOD(ll_G,4)
         ider=0
         wk_ffnl(:,:)=pwff_spl(:,:,ll_G,kln)
         call splfit(qgrid_spl,derfun,wk_ffnl,ider,qpg_norm,newfun,nq_spl,npw)
         do mm_G=-ll_G,ll_G
           ilm_G=1+ll_G**2+ll_G+mm_G
           ignt=gntselect(ilm_G,klm)
           if (ignt==0) CYCLE
           rgnt=realgnt(ignt)
           !
           ! === Evaluate matrix elements for each plane wave ===
           do ig=1,npw
             dummy = newfun(ig) * ylm_q(ig,ilm_G) * rgnt
             paw_rhox(1,ig,klmn) = paw_rhox(1,ig,klmn) &
&               + ( dummy * mi_l(1,ipow) )

             paw_rhox(2,ig,klmn) = paw_rhox(2,ig,klmn) &
&               + ( dummy * mi_l(2,ipow) )
           end do
         end do !mm_G
       end do !ll_G

     end do !ilmn
   end do !jlmn

   !
   ! * Multiply by 4\pi arising from the expansion of the plane wave
   paw_rhox = four_pi*paw_rhox
   ABI_FREE(realgnt)
   ABI_FREE(gntselect)

 CASE (PWIJ_SHISHKIN)
   ! === Shishkin-Kresse approximated expression ====
   ! * Better description of multipoles of AE charge,
   ! * Better results for energy degeneracies in GW band structure
   ! * $4\pi \sum_{LM} q_ij^{LM} Y_M^L(q+G) f^{aL}_{ij}(q+G)$ where f has been calculated in paw_mkrhox_spl
   !
   paw_rhox=zero
   !
   ! === Loop on (jl,jm,jn) channels for this atom ===
   !itypat=Cryst%typat(iatm)
   do jlmn=1,Pawtab(itypat)%lmn_size
     jl =indlmn(1,jlmn)
     jm =indlmn(2,jlmn)
     jlm=indlmn(4,jlmn)
     jln=indlmn(5,jlmn)

     k0lmn=jlmn*(jlmn-1)/2
     k0lm =jlm *(jlm -1)/2
     k0ln =jln *(jln -1)/2
     !
     ! === Loop on (il,im,in) channels; klmn is index for packed form ===
     do ilmn=1,jlmn
       il =indlmn(1,ilmn)
       im =indlmn(2,ilmn)
       ilm=indlmn(4,ilmn)
       iln=indlmn(5,ilmn)

       klmn=k0lmn+ilmn
       klm =k0lm +ilm
       kln =k0ln +iln
       !
       ! === Summing over allowed (l,m), taking into account Gaunt selection rules ===
       do ll_G=ABS(jl-il),jl+il,2
         ipow=MOD(ll_G,4)
         do mm_G=-ll_G,ll_G
           ! here I can move splfit before the loop over mm_G but I have to change paw_rhox_spl
           ilm_G=1+ll_G**2+ll_G+mm_G
           ider=0
           wk_ffnl(:,:)=pwff_spl(:,:,ilm_G-1,klmn)  ! Note klmn and ilm_G-1
           call splfit(qgrid_spl,derfun,wk_ffnl,ider,qpg_norm,newfun,nq_spl,npw)
           !
           ! === Evaluate matrix elements for each plane wave ===
           do ig=1,npw
             dummy = newfun(ig) * ylm_q(ig,ilm_G)
             paw_rhox(1,ig,klmn) = paw_rhox(1,ig,klmn) &
&              + dummy * mi_l(1,ipow) !(ph3d(1,ig)*mi_l(1,ipow)-ph3d(2,ig)*mi_l(2,ipow))

             paw_rhox(2,ig,klmn) = paw_rhox(2,ig,klmn) &
&              + dummy * mi_l(2,ipow) !(ph3d(1,ig)*mi_l(2,ipow)+ph3d(2,ig)*mi_l(1,ipow))
           end do
         end do !mm_G
       end do !ll_G

     end do !ilmn
   end do !jlmn
   !
   ! * Multiply by 4\pi arising from the expansion of the plane wave
   paw_rhox=four_pi*paw_rhox

 CASE DEFAULT
   write(msg,'(a,i3)')' Wrong value for method= ',method
   MSG_BUG(msg)
 END SELECT

 ABI_FREE(wk_ffnl)
 ABI_FREE(newfun)
 ABI_FREE(derfun)
 ABI_FREE(qpg_norm)

end subroutine paw_mkrhox
!!***

!----------------------------------------------------------------------

!!****f* m_pawpwij/paw_rho_tw_g
!! NAME
!! paw_rho_tw_g
!!
!! FUNCTION
!!  Evaluates the PAW onsite contribution to the oscillator strengths:
!!  sum_{i,j} <\tpsi_{k-q,b1}|\cprj_i> <\cprj_j|\tpsi_{k,b2}>*
!!   \[ <\phi_i|e^{-i(q+G).r}|\phi_j> - <\tilde\phi_i|e^{-i(q+G).r}|\tilde\phi_j> \].
!!
!! INPUTS
!! dim_rtwg=Define the size of the array rhotwg
!!   === for nspinor==1 ===
!!    dim_rtwg=1
!!   === for nspinor==2 ===
!!    dim_rtwg=2 if only <up|up>, <dwn|dwn> matrix elements are required
!!    dim_rtwg=4 for <up|up>, <dwn|dwn>, <up|dwn> and <dwn|up>.
!!  nspinor=number of spinorial components.
!!  npw=number of plane waves for oscillator matrix elements
!!  natom=number of atoms
!!  Cprj_kmqb1(natom,nspinor),Cprj_kb2(natom,nspinor) <type(pawcprj_type)>=
!!   projected input wave functions <Proj_i|Cnk> with all NL projectors corresponding to
!!   wavefunctions (k-q,b1,s) and (k,b2,s), respectively.
!!
!! SIDE EFFECTS
!!  rhotwg(npw*dim_rtwg)=Updated oscillator strengths with the on-site PAW contributions added.
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cchi0,cchi0q0,cchi0q0_intraband
!!      check_completeness,cohsex_me,exc_build_block,exc_build_ham,m_shirley
!!      prep_calc_ucrpa
!!
!! SOURCE


pure subroutine paw_rho_tw_g(npw,dim_rtwg,nspinor,natom,ntypat,typat,xred,gvec,Cprj_kmqb1,Cprj_kb2,Pwij,rhotwg)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat,npw,nspinor,dim_rtwg
!arrays
 integer,intent(in) :: gvec(3,npw),typat(natom)
 real(dp),intent(in) :: xred(3,natom)
 complex(gwpc),intent(inout) :: rhotwg(npw*dim_rtwg)
 type(pawcprj_type),intent(in) :: Cprj_kmqb1(natom,nspinor),Cprj_kb2(natom,nspinor)
 type(pawpwij_t),intent(in) :: Pwij(ntypat)

!Local variables-------------------------------
!scalars
 integer :: ig,iat,nlmn,ilmn,jlmn,k0lmn,klmn,iab,isp1,isp2,spad,itypat
 real(dp) :: fij,re_psp,im_psp,re_pw,im_pw,arg
!arrays
 integer,parameter :: spinor_idxs(2,4)=RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
 real(dp) :: tmp(2),qpg(3),x0(3),ph3d(2)

! *************************************************************************

 ! === Loop over the four spinorial combinations ===
 do iab=1,dim_rtwg
   isp1=spinor_idxs(1,iab)
   isp2=spinor_idxs(2,iab)
   spad=npw*(iab-1)

   do ig=1,npw
     tmp(:)=zero
     do iat=1,natom
       itypat = typat(iat)
       nlmn   = Pwij(itypat)%lmn_size
       x0(:) = xred(:,iat)

       ! === Structure factor e^{-i(q+G)*xred} ===
       qpg(:)= Pwij(itypat)%qpt(:) + gvec(:,ig)
       arg=-two_pi*DOT_PRODUCT(qpg(:),x0)
       ph3d(1)=COS(arg)
       ph3d(2)=SIN(arg)

       ! === Loop on [(jl,jm,jn),(il,im,in)] channels. packed storage mode ===
       do jlmn=1,nlmn
         k0lmn=jlmn*(jlmn-1)/2
         do ilmn=1,jlmn
           re_psp =  Cprj_kmqb1(iat,isp1)%cp(1,ilmn) * Cprj_kb2(iat,isp2)%cp(1,jlmn) &
&                   +Cprj_kmqb1(iat,isp1)%cp(2,ilmn) * Cprj_kb2(iat,isp2)%cp(2,jlmn) &
&                   +Cprj_kmqb1(iat,isp1)%cp(1,jlmn) * Cprj_kb2(iat,isp2)%cp(1,ilmn) &
&                   +Cprj_kmqb1(iat,isp1)%cp(2,jlmn) * Cprj_kb2(iat,isp2)%cp(2,ilmn)

           im_psp =  Cprj_kmqb1(iat,isp1)%cp(1,ilmn) * Cprj_kb2(iat,isp2)%cp(2,jlmn) &
&                   -Cprj_kmqb1(iat,isp1)%cp(2,ilmn) * Cprj_kb2(iat,isp2)%cp(1,jlmn) &
&                   +Cprj_kmqb1(iat,isp1)%cp(1,jlmn) * Cprj_kb2(iat,isp2)%cp(2,ilmn) &
&                   -Cprj_kmqb1(iat,isp1)%cp(2,jlmn) * Cprj_kb2(iat,isp2)%cp(1,ilmn)

           klmn=k0lmn+ilmn; fij=one; if (jlmn==ilmn) fij=half

           ! Multiply by the phase due to the atom position.
           re_pw =  Pwij(itypat)%mqpgij(1,ig,klmn) * ph3d(1) &
                   -Pwij(itypat)%mqpgij(2,ig,klmn) * ph3d(2)

           im_pw =  Pwij(itypat)%mqpgij(1,ig,klmn) * ph3d(2) &
                   +Pwij(itypat)%mqpgij(2,ig,klmn) * ph3d(1)

           tmp(1)=tmp(1)+ fij * (re_pw*re_psp - im_pw*im_psp)
           tmp(2)=tmp(2)+ fij * (re_pw*im_psp + im_pw*re_psp)
         end do !ilmn
       end do !jlmn
     end do !iat
     !
     !if(ig==1) write(std_out,*) " TOTAL PW     osc str = ",rhotwg(ig+spad)
     !if(ig==1) write(std_out,*) " TOTAL PAW    osc str = ",tmp(1),tmp(2)
     ! Update input data using the appropriate index.
     rhotwg(ig+spad) = rhotwg(ig+spad) + CMPLX(tmp(1),tmp(2),kind=gwpc)
     !if(ig==1) write(std_out,*) " TOTAL PW+PAW osc str = ",rhotwg(ig+spad)
   end do !ig

 end do !dim_rtwg

end subroutine paw_rho_tw_g
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/paw_cross_rho_tw_g
!! NAME
!!  paw_cross_rho_tw_g
!!
!! FUNCTION
!!  Compute the cross term between the PAW onsite part and plane-wave part
!!  in rho_tw
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cchi0,cchi0q0,prep_calc_ucrpa
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many,fftpad
!!
!! SOURCE

subroutine paw_cross_rho_tw_g(nspinor,npwvec,nr,ngfft,map2sphere,use_padfft,igfftg0,gbound,&
& ur_ae1,ur_ae_onsite1,ur_ps_onsite1,i1,ktabr1,ktabp1,spinrot1,&
& ur_ae2,ur_ae_onsite2,ur_ps_onsite2,i2,ktabr2,ktabp2,spinrot2,&
& dim_rtwg,rhotwg)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: i1,i2,npwvec,nr,nspinor,dim_rtwg,map2sphere,use_padfft
 complex(dpc),intent(in) :: ktabp1,ktabp2
!arrays
 integer,intent(in) :: gbound(:,:)
 integer,intent(in) :: igfftg0(npwvec*map2sphere),ngfft(18)
 integer,intent(in) :: ktabr1(nr),ktabr2(nr)
 real(dp),intent(in) :: spinrot1(4),spinrot2(4)
 complex(gwpc),intent(in) :: ur_ae1(nr),ur_ae2(nr)
 complex(gwpc),intent(in) :: ur_ae_onsite1(nr),ur_ae_onsite2(nr)
 complex(gwpc),intent(in) :: ur_ps_onsite1(nr),ur_ps_onsite2(nr)
 complex(gwpc),intent(inout) :: rhotwg(npwvec*dim_rtwg)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1=1
 integer :: ig,igfft,nx,ny,nz,ldx,ldy,ldz,mgfft,isprot1,isprot2
 type(fftbox_plan3_t) :: plan
!arrays
 complex(dpc),allocatable :: usk(:),uu(:),rho(:)

! *************************************************************************

 SELECT CASE (nspinor)

 CASE (1) ! Collinear case.

   ABI_MALLOC(uu,(nr))
   ABI_MALLOC(usk,(nr))
   ABI_MALLOC(rho,(nr))

   uu  = ur_ae1(ktabr1)*ktabp1        - ur_ae_onsite1(ktabr1)*ktabp1; if (i1==1) uu  = CONJG(uu)
   usk = ur_ae_onsite2(ktabr2)*ktabp2 - ur_ps_onsite2(ktabr2)*ktabp2; if (i2==2) usk = CONJG(usk)
   rho = uu * usk

   uu  = ur_ae_onsite1(ktabr1)*ktabp1 - ur_ps_onsite1(ktabr1)*ktabp1; if (i1==1) uu  = CONJG(uu)
   usk = ur_ae2(ktabr2)*ktabp2        - ur_ae_onsite2(ktabr2)*ktabp2; if (i2==2) usk = CONJG(usk)
   rho = rho + uu * usk

   SELECT CASE (map2sphere)

   CASE (0) ! Need results on the full FFT box thus cannot use zero-padded FFT.

     call fftbox_plan3_many(plan,ndat1,ngfft(1:3),ngfft(1:3),ngfft(7),-1)
     call fftbox_execute(plan,rho)

     rhotwg=rhotwg + rho

   CASE (1) ! Need results on the G-sphere. Call zero-padded FFT routines if required.

     if (use_padfft==1) then
       nx =ngfft(1); ny =ngfft(2); nz =ngfft(3); mgfft = MAXVAL(ngfft(1:3))
       ldx=nx      ; ldy=ny      ; ldz=nz
       call fftpad(rho,ngfft,nx,ny,nz,ldx,ldy,ldz,ndat1,mgfft,-1,gbound)
     else

       call fftbox_plan3_many(plan,ndat1,ngfft(1:3),ngfft(1:3),ngfft(7),-1)
       call fftbox_execute(plan,rho)
     end if

     do ig=1,npwvec       ! Have to map FFT to G-sphere.
       igfft=igfftg0(ig)
       if (igfft/=0) then ! G-G0 belong to the FFT mesh.
         rhotwg(ig)=rhotwg(ig)+rho(igfft)
       end if
     end do

   CASE DEFAULT
     MSG_BUG("Wrong map2sphere")
   END SELECT

   RETURN

 CASE (2) ! Spinorial case.

   isprot1=spinrot1(1); isprot2=spinrot2(1) ! This is to bypass abirule
   MSG_ERROR("Spinorial case not implemented yet")

   SELECT CASE (map2sphere)

   CASE (0) ! Need results on the full FFT box thus cannot use zero-padded FFT.
   CASE (1) ! Need results on the G-sphere. Call zero-padded FFT routines if required.
   CASE DEFAULT
     MSG_BUG("Wrong map2sphere")
   END SELECT

   RETURN

 CASE DEFAULT
   MSG_BUG('Wrong nspinor')
 END SELECT

end subroutine paw_cross_rho_tw_g
!!***

!----------------------------------------------------------------------

END MODULE m_pawpwij
!!***
