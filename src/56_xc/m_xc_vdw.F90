!!****m* ABINIT/m_xc_vdw
!! NAME
!!  m_xc_vdw
!!
!! FUNCTION
!!  Calculates van der Waals corrections to exchange-correlation.
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2020 ABINIT group (Yann Pouillon, Camilo Espejo)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  DRSLL04 = doi:10.1103/PhysRevLett.92.246401 [[cite:Dion2004]]
!!  RS09 = doi:10.1103/PhysRevLett.103.096102 [[cite:Romanperez2009]]
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
#include "abi_xc_vdw.h"

module m_xc_vdw

 use defs_basis
 use iso_c_binding
 use m_abicore
 use m_errors
 use libxc_functionals
 use m_splines
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_io_tools,      only : flush_unit, open_file
 use m_numeric_tools, only : simpson_int, cspint
 use m_integrals,     only : radsintr

 implicit none

#if defined DEV_YP_VDWXC

  public :: &
&   xc_vdw_type, &
&   xc_vdw_aggregate, &
&   xc_vdw_energy, &
&   xc_vdw_done, &
&   xc_vdw_get_params, &
&   xc_vdw_init, &
&   xc_vdw_libxc_init, &
&   xc_vdw_memcheck, &
&   xc_vdw_read, &
&   xc_vdw_set_functional, &
&   xc_vdw_show, &
&   xc_vdw_status, &
&   xc_vdw_trigger, &
&   xc_vdw_write

  private
!!***

!!****t* ABINIT/m_xc_vdw/xc_vdw_type
!! NAME
!!  xc_vdw_type
!!
!! FUNCTION
!!  Tunable parameters for calculations relying upon van der Waals XC.
!!
!! SOURCE

  type xc_vdw_type

    integer :: functional = 0
    ! Functional to use for the van der Waals correction
    ! (see the description of the vdw_xc input variable for possible values)

    real(dp) :: zab = -0.8491_dp
    ! Parameter of the vdW functional, introduced by DRSLL04

    integer :: ndpts = 20
    ! Number of points for the d-mesh

    real(dp) :: dcut = 30.0_dp
    ! Cut-off for the d-mesh

    real(dp) :: dratio = 20.0_dp
    ! Ratio between highest and lowest d

    real(dp) :: dsoft = 1.0_dp
    ! Distance within the d-mesh below which the kernel will be softened

    real(dp) :: phisoft = -1.0_dp
    ! Value of the softened kernel for d=0
    ! Will be set automatically if negative (default behaviour)

    integer :: nqpts = 30
    ! Number of points of the q-mesh

    real(dp) :: qcut = 5.0_dp
    ! Cut-off for the q-mesh

    real(dp) :: qratio = 20.0_dp
    ! Ratio between highest and lowest q

    integer :: nrpts = 2048
    ! Number of points in real space for the kernel

    real(dp) :: rcut = 100.0_dp
    ! Real-space cut-off for the kernel

    real(dp) :: rsoft = 0.0_dp
    ! Radius below which the kernel will be softened

    integer :: ngpts = -1
    ! Number of wavevectors for the kernel

    real(dp) :: gcut = 5.0_dp
    ! Wavevector cut-off for the kernel (=sqrt(ecut))

    real(dp) :: acutmin = 10.0_dp
    ! Minimum cut-off for the angular mesh

    real(dp) :: aratio = 30.0_dp
    ! Ratio between highest and lowest a !DEBUG this definition has to
    ! be corrected since its use in the code is different.

    real(dp) :: damax = 0.5_dp
    ! Maximum variation in the angular mesh

    real(dp) :: damin = 1.0e-2_dp
    ! Minimum variation in the angular mesh

    integer :: nsmooth = 12
    ! Saturation level to smoothen q0 near qcut

    real(dp) :: tolerance = 1.0e-13_dp
    ! Global numerical tolerance for the boundary values of the kernel

    integer :: tweaks = 0
    ! Tweaks of the implementation (modify with extreme care)

  end type xc_vdw_type
!!***

!!****t* ABINIT/m_xc_vdw/vdw_df_tweaks_type
!! NAME
!!  vdw_df_tweaks_type
!!
!! FUNCTION
!!  Tweaks for van der Waals XC calculations (use at your own risks).
!!
!! SOURCE

  type vdw_df_tweaks_type

    integer :: amesh_type
    ! A-mesh type

    integer :: deriv_type
    ! Derivation mode

    integer :: dmesh_type
    ! D-mesh type

    integer :: interp_type
    ! Interpolation mode

    integer :: qmesh_type
    ! Q-mesh type

    integer :: run_type
    ! Run mode

  end type vdw_df_tweaks_type
!!***

  ! Internal vdW-DF parameters, protected from outside
  type(xc_vdw_type),save :: my_vdw_params

  ! Internal vdW-DF tweaks, protected from outside
  type(vdw_df_tweaks_type) :: my_vdw_tweaks

  ! D-mesh
  real(dp),allocatable,save :: dmesh(:)

  ! Q-mesh
  real(dp),allocatable,save :: qmesh(:)

  ! Polynomial basis to inpterpolate on the Q-mesh
  real(dp),allocatable,save :: qpoly_basis(:,:,:)

  ! Phi, softened and unsoftened
  ! Actually contains 4 arrays:
  !   * phi
  !   * dphi/dd1
  !   * dphi/dd2
  !   * d2phi/dd1dd2
  real(dp),allocatable,save :: phi(:,:,:)
  real(dp),allocatable,save :: phi_u(:,:,:)

  ! Bicubic spline representation of Phi (for interpolation)
  real(dp),allocatable,save :: phi_bicubic(:,:,:,:)

  ! Bicubic spline representation of Phi_u (for interpolation)
  real(dp),allocatable,save :: phi_u_bicubic(:,:,:,:)

  ! Filtered representation of Phi
  real(dp),allocatable,save :: phir(:,:,:),d2phidr2(:,:,:)
  real(dp),allocatable,save :: phig(:,:,:),d2phidg2(:,:,:)

  ! Unsoftened representation of Phi
  real(dp),allocatable,save :: phir_u(:,:,:)

  type(libxc_functional_type) :: vdw_funcs(2)

! ----------------------------------------------------------------------

! vdW-DF file format
!   * 1.0: original version
!   * 1.1: including unsoftened kernel
  character(len=64),parameter :: vdw_format_string = 'XC vdW-DF 1.1 ABINIT'

! vdW-DF main switch
  logical :: vdw_switch = .false.

! ----------------------------------------------------------------------

contains
!!***

! ----------------------------------------------------------------------
! Public routines
! ----------------------------------------------------------------------

!!****f* ABINIT/m_xc_vdw/xc_vdw_aggregate
!! NAME
!!  xc_vdw_aggregate
!!
!! FUNCTION
!!  Aggregates the various quantities calculated by the other vdW routines and
!!  produces energy, potential and stress tensor.
!!
!! INPUTS
!!  volume= the volume of the cell
!!  gprimd= unit vectors in reciprocal space
!!  npts_rho= number of density points to treat
!!  ngr2= number of components of grho2
!!  nspden= number of spin components
!!  nr1= 1st dimension of real-space mesh
!!  nr2= 2nd dimension of real-space mesh
!!  nr3= 3rd dimension of real-space mesh
!!  rho= electronic density
!!  grho2= gradient of the density
!!  lda_ex= exchange energy computed at the LDA level
!!  lda_ec= correlation energy computed at the LDA level
!!  lda_vx= exchange potential computed at the LDA level
!!  lda_vc= correlation potential computed at the LDA level
!!
!! OUTPUTS
!!  dexc= vdW correction to the exchange-correlation energy
!!  dexcg= vdW correction to the exchange-correlation potential
!!  theta= see RS09
!!
!! NOTES
!!  exc_vdw includes deltae_vdw.
!!
!! PARENTS
!!      rhotoxc
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine xc_vdw_aggregate(volume,gprimd,npts_rho,nspden,ngrad,nr1,nr2,nr3, &
& rho_grho,deltae_vdw,exc_vdw,decdrho_vdw,decdgrho_vdw,stress_vdw)

#if defined HAVE_FFTW3
  include "fftw3.f03"
#endif

!Arguments ------------------------------------
  integer,intent(in) :: npts_rho,nspden,ngrad,nr1,nr2,nr3
  real(dp),intent(in) :: volume,gprimd(3,3)
  real(dp),intent(in) :: rho_grho(npts_rho,nspden,ngrad)
  real(dp),intent(out) :: exc_vdw,deltae_vdw
  real(dp),intent(out) :: decdrho_vdw(nspden),decdgrho_vdw(3,nspden)
  real(dp),intent(out),optional :: stress_vdw(3,3)

!Local variables ------------------------------
  character(len=512) :: msg
  integer :: ig,ip1,ip2,iq1,iq2,ir1,ir2,ir3,is1,is2,ngpts,nqpts
  integer(kind=SIZEOF_PTRDIFF_T) :: fftw3_plan
  real(dp) :: a1,a2,a3,b1,b2,b3,gtmp,gcut,sg
  real(dp) :: ex,ec,vx,vc
  real(dp) :: exc_nl,eps_vdw,deltae_uns,dexc,dexcg(3)
  real(dp) :: dg,dvol,rho_tmp,gtol
  real(dp) :: exc_tmp,decdrho_tmp(nspden),decdgrho_tmp(3,nspden)
  real(dp),allocatable :: exc_lda(:,:),vxc_lda(:,:),vxcg_lda(:,:,:)
  real(dp),allocatable :: gvec(:,:),theta(:,:,:)
  real(dp),allocatable :: t3dr(:,:,:,:)
  complex(dp),allocatable :: t3dg(:,:,:,:)
  real(dp),allocatable :: ptmp(:,:,:),ttmp(:,:),utmp(:),wtmp(:)

! *************************************************************************

  DBG_ENTER("COLL")

  ! Init
  ngpts = my_vdw_params%ngpts
  nqpts = my_vdw_params%nqpts
  gcut = my_vdw_params%gcut
  gtol = my_vdw_params%tolerance
  exc_vdw = zero
  deltae_uns = zero
  deltae_vdw = zero
  if ( present(stress_vdw) ) stress_vdw(:,:) = zero
  decdgrho_tmp(:,:) = zero
  decdrho_tmp(:) = zero
  dvol = volume / npts_rho

  ! Check the status of the switch
  if ( .not. vdw_switch ) return

  ABI_ALLOCATE(exc_lda,(2,npts_rho))
  ABI_ALLOCATE(vxc_lda,(2,npts_rho))
  ABI_ALLOCATE(vxcg_lda,(2,3,npts_rho))
  ABI_ALLOCATE(gvec,(3,npts_rho))
  ABI_ALLOCATE(theta,(nqpts,nspden,5))
  ABI_ALLOCATE(t3dr,(nr1,nr2,nr3,nqpts))
  ABI_ALLOCATE(t3dg,(nr1,nr2,nr3,nqpts))
  ABI_ALLOCATE(ttmp,(nqpts,nr1*nr2*nr3))
  ABI_ALLOCATE(ptmp,(nqpts,nqpts,2))
  ABI_ALLOCATE(utmp,(nqpts))
  ABI_ALLOCATE(wtmp,(nqpts))

  ! Calculate XC energy density from LDA / GGA
  call vdw_df_ldaxc(npts_rho,nspden,ngrad,rho_grho,exc_lda,vxc_lda,vxcg_lda)

  ! Build theta in 3D and g-vectors
  if ( npts_rho /= nr1*nr2*nr3 ) then
    MSG_WARNING('The 3D reconstruction of the density might be wrong (npts /= nr1*nr2*nr3)')
  end if
#if defined DEBUG_VERBOSE
  write(msg,'(a,1x,i8,1x,a)') "Will now call xc_vdw_energy", &
&   nr1 * nr2 * nr3,"times"
  MSG_COMMENT(msg)
#endif

  ip1 = 1
  ip2 = 1
  do ir3=1,nr3
    do ir2=1,nr2
      do ir1=1,nr1
        gvec(:,ip1) = gprimd(:,1) * ir1 + gprimd(:,2) * ir2 + gprimd(:,3) * ir3
        ex = exc_lda(1,ip1)
        ec = exc_lda(2,ip1)
        vx = vxc_lda(1,ip1)
        vc = vxc_lda(2,ip1)
        theta(:,:,:) = zero
        call xc_vdw_energy(nspden,rho_grho(ip1,1:nspden,1), &
&         rho_grho(ip1,1:nspden,2:ngrad), &
&         ex,ec,vx,vc,theta,eps_vdw)
        t3dr(ir1,ir2,ir3,1:nqpts) = theta(1:nqpts,1,1)
        rho_tmp = sum(rho_grho(ip1,1:nspden,1))
        deltae_uns = deltae_uns + rho_tmp * eps_vdw * dvol
        ip1 = ip1 + 1
      end do
    end do
#if defined DEBUG_VERBOSE
    if ( (ip1 - ip2) > 100 ) then
      write(msg,'(1x,a,1x,i3,"% complete")') &
&       '[vdW-DF Energy]',int(ip1*100.0/(nr1*nr2*nr3))
      call wrtout(std_out,msg,'COLL')
      ip2 = ip1
    end if
#endif
  end do

  ! Fourier-transform theta
#if defined HAVE_FFTW3
  do iq1=1,nqpts
    call dfftw_plan_dft_r2c_3d(fftw3_plan,nr1,nr2,nr3, &
&     t3dr(:,:,:,iq1),t3dg(:,:,:,iq1),FFTW_ESTIMATE)
    call dfftw_execute(fftw3_plan)
    call dfftw_destroy_plan(fftw3_plan)
    t3dg(:,:,:,iq1) = t3dg(:,:,:,iq1) / (nr1 * nr2 * nr3)
  end do
#else
  MSG_ERROR('vdW-DF calculations require FFTW3 support')
#endif

  ! Repack theta
  ip1 = 1
  do ir3=1,nr3
    do ir2=1,nr2
      do ir1=1,nr1
        ttmp(:,ip1) = t3dg(ir1,ir2,ir3,:)
        ip1 = ip1 + 1
      end do
    end do
  end do

  ! Reset 3D counters, since theta will be reconstructed in 3D on the fly
  ir1 = 1
  ir2 = 1
  ir3 = 1

  ! Go through reciprocal vectors
  ! FIXME: get values of G-vectors
  do ip1=1,npts_rho
    if ( (ir3 > nr3) .and. (ir1 /= 1) .and. (ir2 /= 1) ) then
      MSG_WARNING('buffer overflow reconstructing theta in 3D')
    end if

    gtmp = sqrt(sum(gvec(:,ip1)**2))

    if ( gtmp < gcut ) then

      ! Interpolate phi in reciprocal space:
      !   * ptmp(:,:,1) = phi(g)
      !   * ptmp(:,:,2) = dphi(g)/dg
      ! Note: do this on the fly, to go as fast as possible (this is equivalent
      !       to a call of the 'splint' routine)
      ptmp(:,:,:) = zero
      dg = pi / my_vdw_params%rcut
      ig = int(gtmp * ngpts / gcut)
      a1 = ((ig + 1) * dg - gtmp) / dg
      b1 = one - a1
      a2 = (3 * a1**2 - one) * dg / six
      b2 = (3 * b1**2 - one) * dg / six
      a3 = (a1**3 - a1) * dg**2 / six
      b3 = (b1**3 - b1) * dg**2 / six
      do iq2 = 1,nqpts
        do iq1 = 1,iq2
          ptmp(iq1,iq2,1) = a1 * phig(ig,iq1,iq2) + b1 * phig(ig+1,iq1,iq2) &
&           + a3 * d2phidg2(ig,iq1,iq2) + b3 * d2phidg2(ig+1,iq1,iq2)
          ptmp(iq1,iq2,2) = (phig(ig+1,iq1,iq2) - phig(ig,iq1,iq2)) / dg &
&           - a2 * d2phidg2(ig,iq1,iq2) + b2 * d2phidg2(ig+1,iq1,iq2)
        end do
      end do

      do iq2 = 1,nqpts
        do iq1 = 1,iq2-1
          ptmp(iq2,iq1,:) = ptmp(iq1,iq2,:)
        end do
      end do

      ! Calculate contributions to integral in Fourier space:
      ! FIXME: find back # of integral in paper
      utmp(:) = matmul(ttmp(:,ip1),ptmp(:,:,1))

      ! Calculate contribution to stress in reciprocal space
      ! Note: contribution of g=0 is zero
      if ( present(stress_vdw) ) then
        if ( gtmp > gtol ) then
          wtmp(:) = matmul(ttmp(:,ip1),ptmp(:,:,2))
          sg = sum(wtmp(:)*ttmp(:,ip1)) * volume / gtmp
          do is2=1,3
            do is1=1,3
              stress_vdw(is1,is2) = stress_vdw(is1,is2) - &
&               sg * gvec(is1,ip1) * gvec(is2,ip1)
            end do
          end do
        end if ! gtmp > gtol
      end if ! present(stress_vdw)

    else

      utmp(:) = zero

    end if ! gtmp < gcut

    ! Reconstruct the integrand in 3D
    t3dg(ir1,ir2,ir3,:) = utmp(:)

    ir1 = ir1 + 1
    if ( ir1 > nr1 ) then
      ir2 = ir2 + 1
      ir1 = 1
    end if
    if ( ir2 > nr2 ) then
      ir3 = ir3 + 1
      ir2 = 1
    end if
  end do !ip1=1,npts_rho

  ! Fourier-transform back the integrand
#if defined HAVE_FFTW3
  !write(msg,'(a)') ch10
  !call wrtout(std_out,msg,'COLL')
  do iq1=1,nqpts
    !write(msg,'(1x,a,i4.4)') "xc_vdw_aggregate: backward FFT #",iq1
    !call wrtout(std_out,msg,'COLL')

    call dfftw_plan_dft_c2r_3d(fftw3_plan,nr1,nr2,nr3, &
&     t3dg(:,:,:,iq1),t3dr(:,:,:,iq1),FFTW_ESTIMATE)
    call dfftw_execute(fftw3_plan)
    call dfftw_destroy_plan(fftw3_plan)
  end do
#else
  MSG_ERROR('vdW-DF calculations require FFTW3 support')
#endif

  ! Repack the integrand
  ip1 = 1
!$OMP  PARALLEL DO COLLAPSE(3) &
!$OMP& DEFAULT(SHARED) PRIVATE(ir1,ir2,ir3)
  do ir3=1,nr3
    do ir2=1,nr2
      do ir1=1,nr1
        ttmp(:,ip1) = t3dr(ir1,ir2,ir3,1:nqpts)
        ip1 = ip1 + 1
      end do
    end do
  end do
!$OMP END PARALLEL DO

#if defined DEBUG_VERBOSE
  write(msg,'(a,1x,i8.8,1x,a)') "Will now call xc_vdw_energy",npts_rho,"times"
  MSG_COMMENT(msg)
#endif

  ! Calculate and integrate vdW corrections at each point
  do ip1=1,npts_rho

    ! Get local contributions
    ex = exc_lda(1,ip1)
    ec = exc_lda(2,ip1)
    vx = vxc_lda(1,ip1)
    vc = vxc_lda(2,ip1)
    theta(:,:,:) = zero
    call xc_vdw_energy(nspden,rho_grho(ip1,1:nspden,1), &
&     rho_grho(ip1,1:nspden,2:ngrad), &
&     ex,ec,vx,vc,theta,exc_tmp,decdrho_tmp,decdgrho_tmp)

    ! Get nonlocal contributons
    ! Note: is2 represents cartesian coordinates here.
    rho_tmp = sum(rho_grho(ip1,1:nspden,1))
    deltae_vdw = deltae_vdw + rho_tmp * exc_tmp * dvol
    exc_vdw = exc_vdw + rho_tmp * &
&     ( half * ( sum(ttmp(:,ip1) * theta(:,1,1)) / &
&       (rho_tmp + tiny(rho_tmp)) ) * dvol ) * dvol
    !Correctness of multiplication by dvol above depends on how fftw3 deals
    !with the volume element in direct or inverse FFT. Here it was assumed
    !that fftw3 does not multiply by the volume upon FFT^-1
    do is1=1,nspden
      decdrho_vdw(is1) = decdrho_vdw(is1) + decdrho_tmp(is1) + &
&       sum(ttmp(:,ip1) * theta(:,is1,2))
      do is2=1,3
        decdgrho_vdw(is2,is1) = decdgrho_vdw(is2,is1) + decdgrho_tmp(is2,is1) + &
&         sum(ttmp(:,ip1) * theta(:,is1,is2+2))
      end do
    end do

#if defined DEBUG_VERBOSE
    if ( modulo(ip1,int(npts_rho/100)) == 0 ) then
      write(msg,'(1x,a,1x,i3,"% complete")') &
&       "[vdW-DF Energy]",int(ip1*100.0/npts_rho)
      call wrtout(std_out,msg,'COLL')
    end if
#endif

  end do ! Loop on npts_rho


  deltae_vdw = deltae_uns + deltae_vdw

#if defined DEBUG_VERBOSE
  write(msg,'(1x,a)') "[vdW-DF Enrgy] 100% complete"
#endif

  ! Display results
  write(msg,'(a,1x,3a,2(3x,a,1x,e12.5,a),3x,a,1(1x,e12.5),a, &
&   3x,a,3(1x,e12.5),a)') &
&   ch10,"[vdW-DF] xc_vdw_aggregate: reporting vdW-DF contributions", &
&   ch10,ch10, &
&   "DeltaE_vdw = ",deltae_vdw,ch10, &
&   "Exc_vdw = ",exc_vdw,ch10, &
&   "dExc_vdw/drho = ",decdrho_vdw(:),ch10, &
&   "dExc_vdw/dgrho = ",decdgrho_vdw(:,:),ch10
  call wrtout(std_out,msg,'COLL')

  ! Final adjustments of stress
  if ( present(stress_vdw) ) then
    stress_vdw(:,:) = stress_vdw(:,:) / volume

    write(msg,'(3x,a,a,3(5x,3(1x,e12.5),a))') &
&   "Stress_vdw = ",ch10, &
&      stress_vdw(1,:),ch10, &
&      stress_vdw(2,:),ch10, &
&      stress_vdw(3,:),ch10
    call wrtout(std_out,msg,'COLL')
  end if

  ! Clean-up the mess
  ABI_DEALLOCATE(exc_lda)
  ABI_DEALLOCATE(vxc_lda)
  ABI_DEALLOCATE(vxcg_lda)
  ABI_DEALLOCATE(gvec)
  ABI_DEALLOCATE(theta)
  ABI_DEALLOCATE(t3dr)
  ABI_DEALLOCATE(t3dg)
  ABI_DEALLOCATE(ptmp)
  ABI_DEALLOCATE(ttmp)
  ABI_DEALLOCATE(utmp)
  ABI_DEALLOCATE(wtmp)

  DBG_EXIT("COLL")

end subroutine xc_vdw_aggregate
!!***

!!****f* ABINIT/m_xc_vdw/xc_vdw_energy
!! NAME
!!  xc_vdw_energy
!!
!! FUNCTION
!!  Calculates the exchange-correlation energy correction due to
!!  van der Waals interactions at one point.
!!
!! INPUTS
!!  nspden= number of spin components
!!  rho= electronic density
!!  grho= gradient of the density
!!  lda_ex= exchange energy computed at the LDA level
!!  lda_ec= correlation energy computed at the LDA level
!!  lda_vx= exchange potential computed at the LDA level
!!  lda_vc= correlation potential computed at the LDA level
!!
!! OUTPUTS
!!  theta= theta and its derivatives (see RS09)
!!  dexc= vdW correction to the exchange-correlation energy
!!  dexcg= vdW correction to the exchange-correlation potential
!!
!! NOTES
!!  The derivatives of theta are calculated only if the optional arguments
!!  are  present.
!!  For performance reasons, this routine should not allocate any variable.
!!
!! PARENTS
!!      m_xc_vdw
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine xc_vdw_energy(nspden,rho,grho,ex_lda,ec_lda,vx_lda,vc_lda, &
& theta,eps,dexc,dexcg)

!Arguments ------------------------------------
  integer,intent(in) :: nspden
  real(dp),intent(in) :: ex_lda,ec_lda,vx_lda,vc_lda
  real(dp),intent(in) :: rho(nspden),grho(nspden,3)
  real(dp),intent(inout) :: theta(my_vdw_params%nqpts,nspden,5)
  real(dp),intent(out),optional :: eps,dexc(2),dexcg(3,nspden)

!Local variables ------------------------------
  character(len=512) :: msg
  logical :: calc_corrections
  integer :: iq0,iq1,iq2,ir,is,ix,nqpts,nrpts,ns
  real(dp) :: decdrho,dexdrho,dvcdrho,dvxdrho,ec,ex,vc,vx
  real(dp) :: dr,kappa,rho_tmp,grho_tmp(3),grho2_tmp,rtol,qcut,zab
  real(dp) :: dq(3),ptmp(2),q0(5),qtmp(2)
  real(dp) :: ztmp(my_vdw_params%nqpts,my_vdw_params%nqpts)
  real(dp) :: qpoly(my_vdw_params%nqpts,3)

! *************************************************************************

  ! Init
  nqpts = my_vdw_params%nqpts
  nrpts = my_vdw_params%nrpts
  ns = my_vdw_params%nsmooth
  qcut = my_vdw_params%qcut
  rtol = my_vdw_params%tolerance
  zab = my_vdw_params%zab

  calc_corrections = .false.
  if ( present(eps) .and. present(dexc) .and. present(dexcg) ) then
    calc_corrections = .true.
  end if

  ! Sum density over spin
  rho_tmp = sum(rho(1:nspden))
  kappa = (three * pi**2 * rho_tmp)**third
  forall(ix=1:3) grho_tmp(ix) = sum(grho(1:nspden,ix))
  grho2_tmp = sum(grho_tmp**2)

  ! Calculate local wavevector q0 of eqs. 11-12 from DRSLL04
  !   * q0(1)   = q0
  !   * q0(2)   = dq0 / drho
  !   * q0(3:5) = dq0 / dgrho2
  ! Notes:
  !   * treating rho->0 separately (divide by 0)
  !   * treating ex_lda->0 separately (divide by 0)
  q0(:) = zero
  if ( (rho_tmp < rtol) .or. (ex_lda < rtol) ) then
    q0(1) = qcut
  else
    q0(1) = kappa * (one + ec_lda / ex_lda - zab / nine * grho2_tmp / &
&     (two * kappa * rho_tmp)**2)

    if ( calc_corrections ) then
     q0(2) = ( (vc_lda - ec_lda) / (rho_tmp * ex_lda) - ec_lda * &
&      (vx_lda - ex_lda) / &
&      (rho_tmp * ex_lda**2) + two * zab / nine * grho2_tmp / &
&      (two * kappa * rho_tmp)**3 * &
&      eight * kappa / three ) * kappa + q0(1) / three * rho_tmp
     q0(3:5) = -two * (zab / nine) / (two * kappa * rho_tmp)**2 * &
&      grho_tmp(:)
    end if
  end if

  if ( q0(1) > qcut ) then
    q0(1) = qcut
    q0(2:5)= zero
  end if

  ! Smoothen q0 near qcut exponentially  eq. 5 RS09. (Saturation)
  !   * qtmp(1) = qs = qcut * (1 - exp(-Sum_i=1,ns 1/i (x/xc)**i))
  !   * qtmp(2) = dqs / dq |q=q0
  qtmp(1) = (q0(1) / qcut) / ns
  qtmp(2) = one / qcut
  do is=ns-1,1,-1
    qtmp(1) = (qtmp(1) + one / is) * q0(1) / qcut
    qtmp(2) = (qtmp(2) * q0(1) + one) / qcut
  end do
  qtmp(2) = qtmp(2) * qcut * exp(-qtmp(1))

  q0(1) = qcut * (one - exp(-qtmp(1)))
  q0(2:5) = q0(2:5) * qtmp(2)

  ! Calculate polynomial coefficients for cubic-spline interpolation at q0
  !   * qpoly(:,1) = coefficients
  !   * qpoly(:,2) = first derivatives
  qpoly(:,:) = zero
  iq0 = vdw_df_indexof(qmesh(:),nqpts,q0(1))
  dq(1) = qmesh(iq0+1) - qmesh(iq0)
  dq(2) = (qmesh(iq0+1) - q0(1)) / dq(1)
  dq(3) = (q0(1) - qmesh(iq0)) / dq(1)
  do iq1=1,nqpts
   qpoly(iq1,1) = dq(2) * qpoly_basis(iq0,iq1,1) + dq(3) * &
&    qpoly_basis(iq0+1,iq1,1) + &
&    ((dq(2)**3 - dq(2)) * qpoly_basis(iq0,iq1,2) + &
&    (dq(3)**3 - dq(3)) * qpoly_basis(iq0+1,iq1,2)) * dq(1)**2 / six

   if ( calc_corrections ) then
    qpoly(iq1,2) = -(qpoly_basis(iq0,iq1,1) - &
&     qpoly_basis(iq0+1,iq1,1)) / dq(1) - &
&     ((three * dq(2)**2 - one) * qpoly_basis(iq0,iq1,2) - &
&     (three * dq(3)**2 - one) * qpoly_basis(iq0+1,iq1,2)) * dq(1) / six
   end if
  end do

  ! Pre-compute integrand for vdW energy correction
  ! by cubic-spline interpolation and store it in
  ! qpoly(:,3). This is done twice: one for the
  ! unsoftened kernel and other one for the
  ! softened kernel
  dr = my_vdw_params%rcut / nrpts
  ztmp(:,:) = zero

  if ( calc_corrections ) then
   do ir=1,nrpts
     do iq1=1,nqpts
       do iq2=1,iq1
         ztmp(iq1,iq2) = ztmp(iq1,iq2) - two * pi * phir(ir,iq1,iq2) * dr
       end do
     end do
   end do
  else
   do ir=1,nrpts
     do iq1=1,nqpts
       do iq2=1,iq1
         ztmp(iq1,iq2) = ztmp(iq1,iq2) + two * pi * phir_u(ir,iq1,iq2) * dr
       end do
     end do
   end do
  end if ! calc_corrections


  do iq1=1,nqpts
    do iq2=1,iq1-1
      ztmp(iq2,iq1) = ztmp(iq1,iq2)
    end do
  end do
  qpoly(:,3) = matmul(qpoly(:,1),ztmp(:,:))

  ! Calculate theta and its derivatives (see RS09)
  !   * theta(:,:,1)   = theta
  !   * theta(:,:,2)   = dtheta / drho
  !   * theta(:,:,3:5) = dtheta / dgrho2
  ! Note: trick to go from Abinit to LibXC conventions

  do is=1,nspden
    theta(:,is,1) = qpoly(:,1) * (2 - nspden) * rho(is)
  end do

  ! Calculate theta derivatives
  if ( calc_corrections ) then
    do is=1,nspden
      theta(:,is,2) = qpoly(:,1) + qpoly(:,2) * q0(2) * rho(is)
      do ix=3,5
        theta(:,is,ix) = qpoly(:,2) * q0(ix) * (2 - nspden) * rho(is)
      end do
    end do
  end if

  ! Calculate energy corrections
    eps = zero
    ptmp(1) = sum(qpoly(:,3) * qpoly(:,1))
    eps = ptmp(1) * rho_tmp
  if ( calc_corrections ) then
     dexc(:) = zero
     dexcg(:,:) = zero
     ptmp(2) = two * sum(qpoly(:,2) * qpoly(:,3))
     dexc(:) = two * ptmp(1) * rho_tmp + ptmp(2) * qtmp(2) * rho_tmp**2
    do ix=3,5
      dexcg(ix-2,:) = ptmp(2) * q0(ix) * rho_tmp**2
    end do
  end if ! calc_corrections

end subroutine xc_vdw_energy
!!***

!!****f* ABINIT/m_xc_vdw/xc_vdw_done
!! NAME
!!  xc_vdw_done
!!
!! FUNCTION
!!  Cleans-up the mess once van der Waals corrections to the energy are
!!  not needed anymore.
!!
!! INPUTS
!!  vdw_params= van der Waals parameters
!!
!! PARENTS
!!      driver,vdw_kernelgen
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine xc_vdw_done(vdw_params)

!Arguments ------------------------------------
  type(xc_vdw_type), intent(in) :: vdw_params

!Local variables ------------------------------
  character(len=512) :: msg
  integer :: i

! *************************************************************************

  DBG_ENTER("COLL")

  call libxc_functionals_end(xc_functionals=vdw_funcs)

  ABI_DEALLOCATE(dmesh)
  ABI_DEALLOCATE(qmesh)
  ABI_DEALLOCATE(qpoly_basis)
  ABI_DEALLOCATE(phi)
  ABI_DEALLOCATE(phi_bicubic)
  ABI_DEALLOCATE(phi_u_bicubic)
  ABI_DEALLOCATE(phir)
  ABI_DEALLOCATE(phir_u)
  ABI_DEALLOCATE(d2phidr2)
  ABI_DEALLOCATE(phig)
  ABI_DEALLOCATE(d2phidg2)

  DBG_EXIT("COLL")

end subroutine xc_vdw_done
!!***

!!****f* ABINIT/m_xc_vdw/xc_vdw_get_params
!! NAME
!!  xc_vdw_get_params
!!
!! FUNCTION
!!  Exports internal vdW-DF parameters.
!!
!! OUTPUT
!!  vdw_params= van der Waals parameters
!!
!! PARENTS
!!      vdw_kernelgen
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine xc_vdw_get_params(vdw_params)

!Arguments ------------------------------------
  type(xc_vdw_type), intent(out)  :: vdw_params

!Local variables ------------------------------
  character(len=512) :: msg

! *************************************************************************

  vdw_params%functional = my_vdw_params%functional
  vdw_params%zab        = my_vdw_params%zab
  vdw_params%ndpts      = my_vdw_params%ndpts
  vdw_params%dcut       = my_vdw_params%dcut
  vdw_params%dratio     = my_vdw_params%dratio
  vdw_params%dsoft      = my_vdw_params%dsoft
  vdw_params%phisoft    = my_vdw_params%phisoft
  vdw_params%nqpts      = my_vdw_params%nqpts
  vdw_params%qcut       = my_vdw_params%qcut
  vdw_params%qratio     = my_vdw_params%qratio
  vdw_params%nrpts      = my_vdw_params%nrpts
  vdw_params%rcut       = my_vdw_params%rcut
  vdw_params%rsoft      = my_vdw_params%rsoft
  vdw_params%ngpts      = my_vdw_params%ngpts
  vdw_params%gcut       = my_vdw_params%gcut
  vdw_params%acutmin    = my_vdw_params%acutmin
  vdw_params%aratio     = my_vdw_params%aratio
  vdw_params%damax      = my_vdw_params%damax
  vdw_params%damin      = my_vdw_params%damin
  vdw_params%nsmooth    = my_vdw_params%nsmooth
  vdw_params%tolerance  = my_vdw_params%tolerance
  vdw_params%tweaks     = my_vdw_params%tweaks

end subroutine xc_vdw_get_params
!!***

!!****f* ABINIT/m_xc_vdw/xc_vdw_init
!! NAME
!!  xc_vdw_init
!!
!! FUNCTION
!!  Calculates the van der Waals kernel.
!!
!! INPUTS
!!  vdw_params= parameters for the van der Waals calculations
!!
!! PARENTS
!!      driver,vdw_kernelgen
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine xc_vdw_init(vdw_params)

!Arguments ------------------------------------
  type(xc_vdw_type), intent(in)  :: vdw_params

!Local variables ------------------------------
  character(len=*),parameter :: dmesh_file = "vdw_df_dmesh.pts"
  character(len=*),parameter :: qmesh_file = "vdw_df_qmesh.pts"

  character(len=512) :: msg
  integer :: id1,id2,iq1,iq2,ir,ir0,ixc,jd1,jd2,ndpts,ngpts,nqpts,nrpts,ntest,sofswt
  real(dp) :: d1,d2,dcut,dd,delta_test,dsoft,phisoft,acutmin,aratio,damax
  real(dp) :: phimm,phimp,phipm,phipp,phipn,phi_tmp,ptmp(2)
  real(dp),allocatable :: deltd(:),kern_spline_der(:,:,:)
  real(dp),allocatable :: btmp(:,:),utmp(:),xde(:)

! *************************************************************************

  DBG_ENTER("COLL")

  ! Global init
  call xc_vdw_set_params(vdw_params)
  call vdw_df_set_tweaks(my_vdw_params%tweaks,my_vdw_tweaks)

  ! Local init
  ndpts   = my_vdw_params%ndpts
  nqpts   = my_vdw_params%nqpts
  nrpts   = my_vdw_params%nrpts
  acutmin = my_vdw_params%acutmin
  aratio  = my_vdw_params%aratio
  damax   = my_vdw_params%damax
  dsoft   = my_vdw_params%dsoft
  phisoft = my_vdw_params%phisoft

  ! Create d-mesh
  ! Note: avoid zero and make sure the last point is exactly dcut
  ABI_ALLOCATE(dmesh,(ndpts))
#if defined DEBUG_VERBOSE
  write(msg,'(1x,a)') "[vdW-DF Build] Generating D-mesh"
  call wrtout(std_out,msg,'COLL')
#endif
  call vdw_df_create_mesh(dmesh,ndpts,my_vdw_tweaks%dmesh_type, &
&   my_vdw_params%dcut,mesh_ratio=my_vdw_params%dratio, &
&   mesh_tolerance=my_vdw_params%tolerance,mesh_file=dmesh_file, &
&   avoid_zero=.true.,exact_max=.true.)

  ! Create q-mesh
  ABI_ALLOCATE(qmesh,(nqpts))
#if defined DEBUG_VERBOSE
  write(msg,'(1x,a)') "[vdW-DF Build] Generating Q-mesh"
  call wrtout(std_out,msg,'COLL')
#endif
  call vdw_df_create_mesh(qmesh,nqpts,my_vdw_tweaks%qmesh_type, &
&   my_vdw_params%qcut,mesh_ratio=my_vdw_params%qratio, &
&   mesh_tolerance=my_vdw_params%tolerance, mesh_file=qmesh_file)

  ! Build polynomial basis for cubic-spline interpolation
  !   * qpoly_basis(:,1) = polynomial basis
  !   * qpoly_basis(:,2) = second derivative
  ABI_ALLOCATE(qpoly_basis,(nqpts,nqpts,2))
  ABI_ALLOCATE(btmp,(nqpts,nqpts))
  ABI_ALLOCATE(utmp,(nqpts))
#if defined DEBUG_VERBOSE
  write(msg,'(1x,a)') "[vdW-DF Build] Generating polynomial basis"
  call wrtout(std_out,msg,'COLL')
#endif
  btmp(:,:) = zero
  forall(iq1=1:nqpts) btmp(iq1,iq1) = one
  utmp(:) = zero
  qpoly_basis(:,:,:) = zero
  do iq1=1,nqpts
    qpoly_basis(iq1,iq1,1) = one
    qpoly_basis(iq1,1,2) = zero
    utmp(1) = zero
    do iq2=2,nqpts-1
      ptmp(1) = (qmesh(iq2) - qmesh(iq2-1)) / &
&     (qmesh(iq2+1) - qmesh(iq2-1))
      ptmp(2) = ptmp(1) * qpoly_basis(iq1,iq2-1,2) + two
      qpoly_basis(iq1,iq2,2) = (ptmp(1) - one) / ptmp(2)
      utmp(iq2) = ( six * ((btmp(iq1,iq2+1) - btmp(iq1,iq2)) / &
&     (qmesh(iq2+1) - qmesh(iq2)) - (btmp(iq1,iq2) - &
&     btmp(iq1,iq2-1)) / (qmesh(iq2) - qmesh(iq2-1))) / &
&     (qmesh(iq2+1) - qmesh(iq2-1)) - ptmp(1) * utmp(iq2-1)) / &
&     ptmp(2)
    end do
    utmp(nqpts) = zero
    qpoly_basis(iq1,nqpts,2) = zero
    do iq2=nqpts-1,1,-1
      qpoly_basis(iq1,iq2,2) = qpoly_basis(iq1,iq2,2) * &
&       qpoly_basis(iq1,iq2+1,2) + utmp(iq2)
    end do
  end do
  ABI_DEALLOCATE(btmp)
  ABI_DEALLOCATE(utmp)

  ! Create kernel and its derivatives
  ! Note: using 4 neighbours for derivatives
  ABI_ALLOCATE(phi,(ndpts,ndpts,4))
#if defined DEBUG_VERBOSE
  write(msg,'(1x,a)') "[vdW-DF Build] Building kernel and its derivatives"
  call wrtout(std_out,msg,'COLL')
#endif

  if ( phisoft < zero ) then
    phisoft = three * half * exp(-dsoft * two / sqrt(two))
    my_vdw_params%phisoft = phisoft
  end if

  ! Building of the softened kernel

  sofswt = 1
  do id1=1,ndpts
    d1 = dmesh(id1)
    phi(id1,id1,1) = vdw_df_kernel(d1,d1,dsoft,phisoft,acutmin,aratio,damax)

    ! Delta(d) should be at least 10^-6
    ! WARNING: changed tol3 to tol6 and changed max to min
    dd = min(my_vdw_params%dratio*d1,tol6)

      phimm = vdw_df_kernel(d1-dd,d1-dd,dsoft,phisoft,acutmin,aratio,damax)
      phipm = vdw_df_kernel(d1+dd,d1-dd,dsoft,phisoft,acutmin,aratio,damax)
      phipp = vdw_df_kernel(d1+dd,d1+dd,dsoft,phisoft,acutmin,aratio,damax)
      phimp = vdw_df_kernel(d1-dd,d1+dd,dsoft,phisoft,acutmin,aratio,damax)

      ! Using five point stencil formula for crossed derivative phi(id1,id2,4)
      phi(id1,id1,2) = (phipp + phipm - phimp - phimm) / (four * dd)
      phi(id1,id1,3) = (phipp - phipm + phimp - phimm) / (four * dd)
      phi(id1,id1,4) = (phipp - phipm - phimp + phimm) / ((two * dd)**2)

    do id2=1,id1-1
      d2 = dmesh(id2)
      phi(id1,id2,1) = vdw_df_kernel(d1,d2,dsoft,phisoft,acutmin,aratio,damax)
      phi(id2,id1,1) = phi(id1,id2,1)

      phimm = vdw_df_kernel(d1-dd,d2-dd,dsoft,phisoft,acutmin,aratio,damax)
      phipm = vdw_df_kernel(d1+dd,d2-dd,dsoft,phisoft,acutmin,aratio,damax)
      phipp = vdw_df_kernel(d1+dd,d2+dd,dsoft,phisoft,acutmin,aratio,damax)
      phimp = vdw_df_kernel(d1-dd,d2+dd,dsoft,phisoft,acutmin,aratio,damax)

      ! Using five point stencil formula for crossed derivative phi(id1,id2,4)
      phi(id1,id2,2) = (phipp + phipm - phimp - phimm) / (four * dd)
      phi(id1,id2,3) = (phipp - phipm + phimp - phimm) / (four * dd)
      phi(id1,id2,4) = (phipp - phipm - phimp + phimm) / ((two * dd)**2)

      phi(id2,id1,2) = phi(id1,id2,2)
      phi(id2,id1,3) = phi(id1,id2,3)
      phi(id2,id1,4) = phi(id1,id2,4)
    end do

    write(msg,'(1x,a,1x,i3,"% complete")') &
&     '[vdW-DF Build]',int(id1*100.0/ndpts)
    call wrtout(std_out,msg,'COLL')
#if defined DEBUG_VERBOSE
    call flush_unit(std_out)
#endif
  end do
  write(msg,'(a)') " "
  call wrtout(std_out,msg,'COLL')

  ! Set boundary conditions
  ! Note: will have to check that the borders are smooth enough
  ! These boundary conditions produce large discontinuities, now testing its removal
  ! phi(ndpts,:,:) = zero
  ! phi(:,ndpts,:) = zero
  ! phi(1,:,2:4)   = zero
  ! phi(:,1,2:4)   = zero
  ! Ar2 results show better convergence of E_vdW than when boundary conditions were used.
  ! Kernel plot show that there are not dicontinuities as well.

#if defined DEBUG_VERBOSE
  write(msg,'(1x,a)') "[vdW-DF Build] Building filtered kernel"
  call wrtout(std_out,msg,'COLL')
#endif

  ! Calculate coefficients for bicubic interpolation
  ABI_ALLOCATE(phi_bicubic,(4,4,ndpts,ndpts))
  call spline_bicubic(ndpts,ndpts,dmesh,dmesh,phi(:,:,1),phi(:,:,2),&
&   phi(:,:,3),phi(:,:,4),phi_bicubic)

  ! Build filtered kernel
  ABI_ALLOCATE(phir,(nrpts,nqpts,nqpts))
  ABI_ALLOCATE(d2phidr2,(nrpts,nqpts,nqpts))
  ABI_ALLOCATE(phig,(nrpts,nqpts,nqpts))
  ABI_ALLOCATE(d2phidg2,(nrpts,nqpts,nqpts))
  call vdw_df_filter(nqpts,nrpts,my_vdw_params%rcut,my_vdw_params%gcut,ngpts,sofswt)
  my_vdw_params%ngpts = ngpts

  ! Find closest indices in dmesh
  ! FIXME: something is missing or should be removed here
  ABI_ALLOCATE(kern_spline_der,(ndpts,ndpts,3))
  ABI_ALLOCATE(deltd,(2))
  ABI_ALLOCATE(xde,(2))

  kern_spline_der(:,:,:) = zero

  do jd1=1,ndpts
    do jd2=1,ndpts

      d1 = dmesh(jd1)
      d2 = dmesh(jd2)

      if ( (d1 < dcut) .or. (d2 < dcut) ) then
        id1 = vdw_df_indexof(dmesh,ndpts,d1)
        id2 = vdw_df_indexof(dmesh,ndpts,d2)

        deltd(1) = dmesh(id1+1) - dmesh(id1)
        deltd(2) = dmesh(id2+1) - dmesh(id2)

        xde(1) = (d1 - dmesh(id1)) / deltd(1)
        xde(2) = (d2 - dmesh(id2)) / deltd(2)

        kern_spline_der(jd1,jd2,1) = phi_bicubic(2,1,id1,id2) / deltd(1)
        kern_spline_der(jd1,jd2,2) = phi_bicubic(1,2,id1,id2) / deltd(2)
        kern_spline_der(jd1,jd2,3) = &
&         phi_bicubic(2,2,id1,id2) / ( deltd(1)*deltd(2) )
      end if

    end do
  end do

  ABI_DEALLOCATE(kern_spline_der)
  ABI_DEALLOCATE(deltd)
  ABI_DEALLOCATE(xde)

  ! Building of the unsoftened kernel

  ! Create unsoftened kernel and its derivatives
  ! Note: using 4 neighbours for derivatives
  ABI_ALLOCATE(phi_u,(ndpts,ndpts,4))
#if defined DEBUG_VERBOSE
  write(msg,'(1x,a)') "[vdW-DF Build] Building unsoftened kernel and its derivatives"
  call wrtout(std_out,msg,'COLL')
#endif

  do id1=1,ndpts
    d1 = dmesh(id1)
    phi_u(id1,id1,1) = vdw_df_kernel_value(d1,d1,acutmin,aratio,damax)

    ! Delta(d) should be at least 10^-6
    ! WARNING: changed tol3 to tol6 and changed max to min
    dd = min(my_vdw_params%dratio*d1,tol6)

      phimm = vdw_df_kernel_value(d1-dd,d1-dd,acutmin,aratio,damax)
      phipm = vdw_df_kernel_value(d1+dd,d1-dd,acutmin,aratio,damax)
      phipp = vdw_df_kernel_value(d1+dd,d1+dd,acutmin,aratio,damax)
      phimp = vdw_df_kernel_value(d1-dd,d1+dd,acutmin,aratio,damax)

      ! Using five point stencil formula for crossed derivative phi(id1,id2,4)
      phi_u(id1,id1,2) = (phipp + phipm - phimp - phimm) / (four * dd)
      phi_u(id1,id1,3) = (phipp - phipm + phimp - phimm) / (four * dd)
      phi_u(id1,id1,4) = (phipp - phipm - phimp + phimm) / ((two * dd)**2)

    do id2=1,id1-1
      d2 = dmesh(id2)
      phi_u(id1,id2,1) = vdw_df_kernel_value(d1,d2,acutmin,aratio,damax)
      phi_u(id2,id1,1) = phi_u(id1,id2,1)

      phimm = vdw_df_kernel_value(d1-dd,d2-dd,acutmin,aratio,damax)
      phipm = vdw_df_kernel_value(d1+dd,d2-dd,acutmin,aratio,damax)
      phipp = vdw_df_kernel_value(d1+dd,d2+dd,acutmin,aratio,damax)
      phimp = vdw_df_kernel_value(d1-dd,d2+dd,acutmin,aratio,damax)

      ! Using five point stencil formula for crossed derivative phi(id1,id2,4)
      phi_u(id1,id2,2) = (phipp + phipm - phimp - phimm) / (four * dd)
      phi_u(id1,id2,3) = (phipp - phipm + phimp - phimm) / (four * dd)
      phi_u(id1,id2,4) = (phipp - phipm - phimp + phimm) / ((two * dd)**2)

      phi_u(id2,id1,2) = phi_u(id1,id2,2)
      phi_u(id2,id1,3) = phi_u(id1,id2,3)
      phi_u(id2,id1,4) = phi_u(id1,id2,4)
    end do

    write(msg,'(1x,a,1x,i3,"% complete")') &
&     '[vdW-DF-Uns Build]',int(id1*100.0/ndpts)
    call wrtout(std_out,msg,'COLL')
#if defined DEBUG_VERBOSE
    call flush_unit(std_out)
#endif
  end do
  write(msg,'(a)') " "
  call wrtout(std_out,msg,'COLL')

  ! Set boundary conditions
  ! Note: will have to check that the borders are smooth enough
  ! These boundary conditions produce large discontinuities, now testing its removal
  ! phi(ndpts,:,:) = zero
  ! phi(:,ndpts,:) = zero
  ! phi(1,:,2:4)   = zero
  ! phi(:,1,2:4)   = zero
  ! Ar2 results show better convergence of E_vdW than when boundary conditions were used.
  ! Kernel plot show that there are not dicontinuities as well.

#if defined DEBUG_VERBOSE
  write(msg,'(1x,a)') "[vdW-DF Build] Building filtered kernel"
  call wrtout(std_out,msg,'COLL')
#endif

  ! Calculate coefficients for bicubic interpolation
  ABI_ALLOCATE(phi_u_bicubic,(4,4,ndpts,ndpts))
  call spline_bicubic(ndpts,ndpts,dmesh,dmesh,phi_u(:,:,1),phi_u(:,:,2),&
&   phi_u(:,:,3),phi_u(:,:,4),phi_u_bicubic)

  ! Build filtered kernel
  ABI_ALLOCATE(phir_u,(nrpts,nqpts,nqpts))
  ! For the local correction we just need the kernel not its
  ! derivatives.
  sofswt = 0
  !ABI_ALLOCATE(d2phidr2,(nrpts,nqpts,nqpts))
  !ABI_ALLOCATE(phig,(nrpts,nqpts,nqpts))
  !ABI_ALLOCATE(d2phidg2,(nrpts,nqpts,nqpts))
  call vdw_df_filter(nqpts,nrpts,my_vdw_params%rcut,my_vdw_params%gcut,ngpts,sofswt)
  !  my_vdw_params%ngpts = ngpts already defined above

  ! The following is not needed for the local correction:
  ! Find closest indices in dmesh
  ! FIXME: something is missing or should be removed here
  ! ABI_ALLOCATE(kern_spline_der,(ndpts,ndpts,3))
  ! ABI_ALLOCATE(deltd,(2))
  ! ABI_ALLOCATE(xde,(2))

  ! kern_spline_der(:,:,:) = zero

  ! do jd1=1,ndpts
  !   do jd2=1,ndpts

  !     d1 = dmesh(jd1)
  !     d2 = dmesh(jd2)

  !     if ( (d1 < dcut) .or. (d2 < dcut) ) then
  !       id1 = vdw_df_indexof(dmesh,ndpts,d1)
  !       id2 = vdw_df_indexof(dmesh,ndpts,d2)

  !       deltd(1) = dmesh(id1+1) - dmesh(id1)
  !       deltd(2) = dmesh(id2+1) - dmesh(id2)

  !       xde(1) = (d1 - dmesh(id1)) / deltd(1)
  !       xde(2) = (d2 - dmesh(id2)) / deltd(2)

  !       kern_spline_der(jd1,jd2,1) = phi_bicubic(2,1,id1,id2) / deltd(1)
  !       kern_spline_der(jd1,jd2,2) = phi_bicubic(1,2,id1,id2) / deltd(2)
  !       kern_spline_der(jd1,jd2,3) = &
  !&         phi_bicubic(2,2,id1,id2) / ( deltd(1)*deltd(2) )
  !     end if

  !   end do
  ! end do

  ! ABI_DEALLOCATE(kern_spline_der)
  ! ABI_DEALLOCATE(deltd)
  ! ABI_DEALLOCATE(xde)

  call vdw_df_internal_checks(1)

  vdw_switch = .true.

  DBG_EXIT("COLL")

end subroutine xc_vdw_init
!!***

!!****f* ABINIT/m_xc_vdw/xc_vdw_libxc_init
!! NAME
!!  xc_vdw_libxc_init
!!
!! FUNCTION
!!  Initializes LibXC parameters for the LDA-based part of vdW-DF
!!  calculations.
!!
!! INPUTS
!!  ixc_vdw= vdW-DF functional to apply
!!
!! OUTPUT
!!  (only writing)
!!
!! SIDE EFFECTS
!!  Internal variable vdw_funcs is set according to specified ixc_vdw.
!!  Internal variable my_vdw_params receives the selected functional.
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine xc_vdw_libxc_init(ixc_vdw)

!Arguments ------------------------------------
  integer,intent(in) :: ixc_vdw

!Local variables ------------------------------
  character(len=*),parameter :: c_vdw1 = "GGA_C_PBE"
  character(len=*),parameter :: x_vdw1 = "GGA_X_PBE_R"
  character(len=*),parameter :: x_vdw2 = "GGA_X_PW86"
  character(len=*),parameter :: x_vdw3 = "GGA_X_C09X"
  integer :: i,ii,ic,ix,ixc
  character(len=500) :: msg

! *************************************************************************

  DBG_ENTER("COLL")

  ! Select LibXC functionals
  select case (ixc_vdw)
    case (1)
      ix = libxc_functionals_getid(x_vdw1)
      ic = libxc_functionals_getid(c_vdw1)
    case (2)
      ix = libxc_functionals_getid(x_vdw2)
      ic = libxc_functionals_getid(c_vdw1)
    case (3)
      ix = libxc_functionals_getid(x_vdw3)
      ic = libxc_functionals_getid(c_vdw1)
      MSG_ERROR('[vdW-DF] C09 not available for now')
    case default
      MSG_ERROR('[vdW-DF] invalid setting of vdw_xc')
  end select
  if ( (ix == -1) .or. (ic == -1) ) then
    MSG_ERROR('[vdW-DF] unable to set LibXC parameters')
  end if
  ixc = -(ix * 1000 + ic)

  ! Propagate to internal parameters
  my_vdw_params%functional = ixc_vdw

  ! XC functional init
  call libxc_functionals_init(ixc,1,xc_functionals=vdw_funcs)

  DBG_EXIT("COLL")

end subroutine xc_vdw_libxc_init
!!***

!!****f* ABINIT/m_xc_vdw/xc_vdw_memcheck
!! NAME
!!  xc_vdw_memcheck
!!
!! FUNCTION
!!  Estimates the memory to be used by the vdW-DF method.
!!
!! INPUTS
!!  unt= unit to write the data to
!!  vp= van der Waals parameters
!!
!! PARENTS
!!      driver,vdw_kernelgen
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine xc_vdw_memcheck(unt,vp)

!Arguments ------------------------------------
  integer,intent(in) :: unt
  type(xc_vdw_type), intent(in), optional :: vp

!Local variables ------------------------------
  character(len=1536) :: msg
  integer :: napts,ndpts,nqpts,nrpts,id1,id2
  real(dp) :: dmax,dtmp,mem_perm,mem_temp
  type(xc_vdw_type) :: my_vp

! *************************************************************************

  if ( present(vp) ) then
    my_vp = vp
  else
    my_vp = my_vdw_params
  end if

  ndpts = my_vp%ndpts
  nqpts = my_vp%nqpts
  nrpts = my_vp%nrpts

  if ( .not. allocated(dmesh) ) return

  dmax = zero
  do id1=1,ndpts
    do id2=1,id1
      dtmp = sqrt(dmesh(id1)**2+dmesh(id2)**2)
      if ( dtmp > dmax ) dmax = dtmp
    end do
  end do
  napts = nint(max(my_vp%acutmin,my_vp%aratio*dmax))

  mem_perm = ( &
&   4 * ndpts * (1 + 5 * ndpts ) + &
&   2 * nqpts + (2 + 4 * nrpts) * nqpts * nqpts &
&   ) * sizeof(one) / 1048576.0_dp

  mem_temp = ( &
&   7 * napts + &
&   4 * nqpts + 3* nqpts * nqpts + &
&   5* nrpts + 2 &
&   ) * sizeof(one) / 1048576.0_dp

  write(msg,'(a,1x,3(a),12(3x,a,1x,i8,1x,"elts",a), &
&   a,10x,a,1x,f10.3,1x,"Mb",a,a,13(3x,a,1x,i8,1x,"elts",a), &
&   a,10x,a,1x,f10.3,1x,"Mb",a,a,14x,a,1x,f10.3,1x,"Mb",2(a),1x,2(a))') &
&   ch10, &
&   '==== vdW-DF memory estimation ====',ch10,ch10, &
&   'd2phidg2 .........',nrpts*nqpts*nqpts,ch10, &
&   'd2phidr2 .........',nrpts*nqpts*nqpts,ch10, &
&   'dmesh ............',ndpts * 4,ch10, &
&   'phi ..............',ndpts*ndpts*4,ch10, &
&   'phi_u ............',ndpts*ndpts*4,ch10, &
&   'phi_bicubic ......',4*4*ndpts*ndpts,ch10, &
&   'phi_u_bicubic ....',4*4*ndpts*ndpts,ch10, &
&   'phig .............',nrpts*nqpts*nqpts,ch10, &
&   'phir .............',nrpts*nqpts*nqpts,ch10, &
&   'phir_u ...........',nrpts*nqpts*nqpts,ch10, &
&   'qmesh ............',nqpts * 4,ch10, &
&   'qpoly_basis.......',nqpts * nqpts * 2,ch10,ch10, &
&   'Permanent =',mem_perm,ch10,ch10, &
&   'amesh ............',napts,ch10, &
&   'amesh_cos ........',napts,ch10, &
&   'amesh_sin ........',napts,ch10, &
&   'btmp .............',nqpts*nqpts,ch10, &
&   'dphida ...........',napts,ch10, &
&   'dphidb ...........',napts,ch10, &
&   'ftmp .............',2*(2*nrpts+1),ch10, &
&   'nu1 ..............',napts,ch10, &
&   'nu2 ..............',napts,ch10, &
&   'qpoly ............',nqpts*3,ch10, &
&   'utmp(q) ..........',nqpts,ch10, &
&   'utmp(r) ..........',nrpts,ch10, &
&   'wtmp .............',nqpts*nqpts*2,ch10,ch10, &
&   'Temporary =',mem_temp,ch10,ch10, &
&   'TOTAL =',mem_perm+mem_temp,ch10, &
&   ch10, &
&   '==================================',ch10

  call wrtout(unt,msg,'COLL')

end subroutine xc_vdw_memcheck
!!***

!!****f* ABINIT/m_xc_vdw/xc_vdw_read
!! NAME
!!  xc_vdw_read
!!
!! FUNCTION
!!  Reads vdW-DF variables from disk.
!!
!! INPUTS
!!  filename= file to read data from
!!
!! TODO
!!  design an extension for ETSF_IO
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine xc_vdw_read(filename)

!Arguments ------------------------------------
  character(len=*), intent(in)  :: filename

!Local variables ------------------------------
  character(len=512) :: msg
  character(len=64) :: format_string
  integer :: ncid,dimids(4),varids(13)
  integer :: nderivs_id,npoly_id,ndpts_id,ngpts_id,nqpts_id,nrpts_id
  integer :: nderivs,npoly,ndpts,ngpts,nqpts,nrpts

! *************************************************************************

  DBG_ENTER("COLL")

#if defined HAVE_NETCDF
  write(msg,'("Reading vdW-DF data from",1x,a)') trim(filename)
  call wrtout(std_out,msg,'COLL')

  ! Open file for reading
  NETCDF_VDWXC_CHECK(nf90_open(trim(filename),NF90_NOWRITE,ncid=ncid))

  ! Get file format and generator
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'file_format',format_string))
  write(msg,'(3x,"File format:",1x,a)') trim(format_string)
  call wrtout(std_out,msg,'COLL')
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'generator',msg))
  write(msg,'(3x,"Generator:",1x,a)') trim(msg)
  call wrtout(std_out,msg,'COLL')

  ! Get attributes
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'acutmin',my_vdw_params%acutmin))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'aratio',my_vdw_params%aratio))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'damax',my_vdw_params%damax))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'damin',my_vdw_params%damin))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'dcut',my_vdw_params%dcut))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'dratio',my_vdw_params%dratio))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'dsoft',my_vdw_params%dsoft))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'functional',my_vdw_params%functional))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'gcut',my_vdw_params%gcut))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'nsmooth',my_vdw_params%nsmooth))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'phisoft',my_vdw_params%phisoft))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'qcut',my_vdw_params%qcut))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'qratio',my_vdw_params%qratio))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'rcut',my_vdw_params%rcut))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'rsoft',my_vdw_params%rsoft))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'tolerance',my_vdw_params%tolerance))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'tweaks',my_vdw_params%tweaks))
  NETCDF_VDWXC_CHECK(nf90_get_att(ncid,NF90_GLOBAL,'zab',my_vdw_params%zab))

  ! Get dimensions
  NETCDF_VDWXC_CHECK(nf90_inq_dimid(ncid,'nderivs',nderivs_id))
  NETCDF_VDWXC_CHECK(nf90_inq_dimid(ncid,'npoly',npoly_id))
  NETCDF_VDWXC_CHECK(nf90_inq_dimid(ncid,'ndpts',ndpts_id))
  NETCDF_VDWXC_CHECK(nf90_inq_dimid(ncid,'ngpts',ngpts_id))
  NETCDF_VDWXC_CHECK(nf90_inq_dimid(ncid,'nqpts',nqpts_id))
  NETCDF_VDWXC_CHECK(nf90_inq_dimid(ncid,'nrpts',nrpts_id))

  ! Get varids
  NETCDF_VDWXC_CHECK(nf90_inq_varid(ncid,'qmesh',varids(2)))
  NETCDF_VDWXC_CHECK(nf90_inq_varid(ncid,'dmesh',varids(3)))
  NETCDF_VDWXC_CHECK(nf90_inq_varid(ncid,'phi',varids(4)))
  NETCDF_VDWXC_CHECK(nf90_inq_varid(ncid,'phi_u',varids(5)))
  NETCDF_VDWXC_CHECK(nf90_inq_varid(ncid,'phi_bicubic',varids(6)))
  NETCDF_VDWXC_CHECK(nf90_inq_varid(ncid,'phi_u_bicubic',varids(7)))
  NETCDF_VDWXC_CHECK(nf90_inq_varid(ncid,'phir',varids(8)))
  NETCDF_VDWXC_CHECK(nf90_inq_varid(ncid,'phir_u',varids(9)))
  NETCDF_VDWXC_CHECK(nf90_inq_varid(ncid,'d2phidr2',varids(10)))
  NETCDF_VDWXC_CHECK(nf90_inq_varid(ncid,'phig',varids(11)))
  NETCDF_VDWXC_CHECK(nf90_inq_varid(ncid,'d2phidg2',varids(12)))
  NETCDF_VDWXC_CHECK(nf90_inq_varid(ncid,'qpoly_basis',varids(13)))

  ! Get dimensions
  NETCDF_VDWXC_CHECK(nf90_inquire_dimension(ncid,nderivs_id,len=nderivs))
  NETCDF_VDWXC_CHECK(nf90_inquire_dimension(ncid,npoly_id,len=npoly))
  NETCDF_VDWXC_CHECK(nf90_inquire_dimension(ncid,ndpts_id,len=ndpts))
  NETCDF_VDWXC_CHECK(nf90_inquire_dimension(ncid,ngpts_id,len=ngpts))
  NETCDF_VDWXC_CHECK(nf90_inquire_dimension(ncid,nqpts_id,len=nqpts))
  NETCDF_VDWXC_CHECK(nf90_inquire_dimension(ncid,nrpts_id,len=nrpts))

  ! Get qmesh
  ABI_ALLOCATE(qmesh,(nqpts))
  NETCDF_VDWXC_CHECK(nf90_get_var(ncid,varids(2),qmesh))

  ! Get dmesh
  ABI_ALLOCATE(dmesh,(ndpts))
  NETCDF_VDWXC_CHECK(nf90_get_var(ncid,varids(3),dmesh))

  ! Get phi
  ABI_ALLOCATE(phi,(ndpts,ndpts,nderivs))
  NETCDF_VDWXC_CHECK(nf90_get_var(ncid,varids(4),phi))

  ! Get phi_u
  ABI_ALLOCATE(phi_u,(ndpts,ndpts,nderivs))
  NETCDF_VDWXC_CHECK(nf90_get_var(ncid,varids(5),phi_u))

  ! Get phi_bicubic
  ABI_ALLOCATE(phi_bicubic,(nderivs,nderivs,ndpts,ndpts))
  NETCDF_VDWXC_CHECK(nf90_get_var(ncid,varids(6),phi_bicubic))

  ! Get phi_u_bicubic
  ABI_ALLOCATE(phi_u_bicubic,(nderivs,nderivs,ndpts,ndpts))
  NETCDF_VDWXC_CHECK(nf90_get_var(ncid,varids(7),phi_u_bicubic))

  ! Get phir
  ABI_ALLOCATE(phir,(nrpts,nqpts,nqpts))
  NETCDF_VDWXC_CHECK(nf90_get_var(ncid,varids(8),phir))

  ! Get phir_u
  ABI_ALLOCATE(phir_u,(nrpts,nqpts,nqpts))
  NETCDF_VDWXC_CHECK(nf90_get_var(ncid,varids(9),phir_u))

  ! Get d2phidr2
  ABI_ALLOCATE(d2phidr2,(nrpts,nqpts,nqpts))
  NETCDF_VDWXC_CHECK(nf90_get_var(ncid,varids(10),d2phidr2))

  ! Get phig
  ABI_ALLOCATE(phig,(nrpts,nqpts,nqpts))
  NETCDF_VDWXC_CHECK(nf90_get_var(ncid,varids(11),phig))

  ! Get d2phidg2
  ABI_ALLOCATE(d2phidg2,(nrpts,nqpts,nqpts))
  NETCDF_VDWXC_CHECK(nf90_get_var(ncid,varids(12),d2phidg2))

  ! Get qpoly_basis
  ABI_ALLOCATE(qpoly_basis,(nqpts,nqpts,npoly))
  NETCDF_VDWXC_CHECK(nf90_get_var(ncid,varids(13),qpoly_basis))

  ! Close file
  NETCDF_VDWXC_CHECK(nf90_close(ncid))

  ! Update my_vdw_params
  my_vdw_params%ndpts = ndpts
  my_vdw_params%ngpts = ngpts
  my_vdw_params%nqpts = nqpts
  my_vdw_params%nrpts = nrpts
#else
  MSG_ERROR('reading vdW-DF variables requires NetCDF')
#endif

  DBG_EXIT("COLL")

end subroutine xc_vdw_read
!!***

!!****f* ABINIT/m_xc_vdw/xc_vdw_set_functional
!! NAME
!!  xc_vdw_set_functional
!!
!! FUNCTION
!!  Sets the vdW-DF parameters directly related to the functional.
!!
!! INPUTS
!!  vdw_func= van der Waals functional
!!  vdw_zab= Zab parameter to use for the calculations
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine xc_vdw_set_functional(vdw_func,vdw_zab)

!Arguments ------------------------------------
  integer, intent(in)  :: vdw_func
  real(dp), optional, intent(in) :: vdw_zab

!Local variables ------------------------------
  character(len=512) :: msg

! *************************************************************************

  my_vdw_params%functional = vdw_func
  if ( present(vdw_zab) ) then
    my_vdw_params%zab = vdw_zab
  end if

end subroutine xc_vdw_set_functional
!!***

!!****f* ABINIT/m_xc_vdw/xc_vdw_set_params
!! NAME
!!  xc_vdw_set_params
!!
!! FUNCTION
!!  Imports external vdW-DF parameters.
!!
!! INPUT
!!  vdw_params= van der Waals parameters
!!
!! PARENTS
!!      m_xc_vdw
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine xc_vdw_set_params(vdw_params)

!Arguments ------------------------------------
  type(xc_vdw_type), intent(in)  :: vdw_params

!Local variables ------------------------------
  character(len=512) :: msg

! *************************************************************************

  my_vdw_params%functional = vdw_params%functional
  my_vdw_params%zab        = vdw_params%zab
  my_vdw_params%ndpts      = vdw_params%ndpts
  my_vdw_params%dcut       = vdw_params%dcut
  my_vdw_params%dratio     = vdw_params%dratio
  my_vdw_params%dsoft      = vdw_params%dsoft
  my_vdw_params%phisoft    = vdw_params%phisoft
  my_vdw_params%nqpts      = vdw_params%nqpts
  my_vdw_params%qcut       = vdw_params%qcut
  my_vdw_params%qratio     = vdw_params%qratio
  my_vdw_params%nrpts      = vdw_params%nrpts
  my_vdw_params%rcut       = vdw_params%rcut
  my_vdw_params%rsoft      = vdw_params%rsoft
  my_vdw_params%ngpts      = vdw_params%ngpts
  my_vdw_params%gcut       = vdw_params%gcut
  my_vdw_params%acutmin    = vdw_params%acutmin
  my_vdw_params%aratio     = vdw_params%aratio
  my_vdw_params%damax      = vdw_params%damax
  my_vdw_params%damin      = vdw_params%damin
  my_vdw_params%nsmooth    = vdw_params%nsmooth
  my_vdw_params%tolerance  = vdw_params%tolerance
  my_vdw_params%tweaks     = vdw_params%tweaks

end subroutine xc_vdw_set_params
!!***

!!****f* ABINIT/m_xc_vdw/xc_vdw_show
!! NAME
!!  xc_vdw_show
!!
!! FUNCTION
!!  Displays the parameters in use for the vdW-DF corrections.
!!
!! INPUTS
!!  unt= unit to write the data to
!!  vp= van der Waals parameters
!!
!! PARENTS
!!      driver,vdw_kernelgen
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine xc_vdw_show(unt,vp)

!Arguments ------------------------------------
  integer,intent(in) :: unt
  type(xc_vdw_type), intent(in), optional :: vp

!Local variables ------------------------------
  character(len=1536) :: msg
  type(xc_vdw_type) :: my_vp

! *************************************************************************

  if ( present(vp) ) then
    my_vp = vp
  else
    my_vp = my_vdw_params
  end if

  write(msg,'(a,1x,3(a),3x,a,1x,i9,a,3x,a,1x,f9.3,a, &
&   5(3x,a,1x,i9,a),14(3x,a,1x,f9.3,a),1(3x,a,1x,i9,a),a,1x,a,a)') &
&   ch10, &
&   '==== vdW-DF internal parameters ====',ch10,ch10, &
&   'functional .............',my_vp%functional,ch10, &
&   'zab ....................',my_vp%zab,ch10, &
&   'd-mesh points ..........',my_vp%ndpts,ch10, &
&   'q-mesh points ..........',my_vp%nqpts,ch10, &
&   'r-mesh points ..........',my_vp%nrpts,ch10, &
&   'g-mesh points ..........',my_vp%ngpts,ch10, &
&   'smoothing iterations ...',my_vp%nsmooth,ch10, &
&   'd-mesh cut-off .........',my_vp%dcut,ch10, &
&   'd-mesh ratio ...........',my_vp%dratio,ch10, &
&   'd-mesh softening .......',my_vp%dsoft,ch10, &
&   'phi softening ..........',my_vp%phisoft,ch10, &
&   'q-mesh cut-off .........',my_vp%qcut,ch10, &
&   'q-mesh ratio ...........',my_vp%qratio,ch10, &
&   'r-mesh cut-off .........',my_vp%rcut,ch10, &
&   'r-mesh softening .......',my_vp%rsoft,ch10, &
&   'g-mesh cut-off .........',my_vp%gcut,ch10, &
&   'a-mesh min cut-off .....',my_vp%acutmin,ch10, &
&   'a-mesh ratio ...........',my_vp%aratio,ch10, &
&   'a-mesh delta max .......',my_vp%damax,ch10, &
&   'a-mesh delta min .......',my_vp%damin,ch10, &
&   'global tolerance .......',my_vp%tolerance,ch10, &
&   'tweaks .................',my_vp%tweaks,ch10, &
&   ch10, &
&   '====================================',ch10

  call wrtout(unt,msg,'COLL')

end subroutine xc_vdw_show
!!***

!!****f* ABINIT/m_xc_vdw/xc_vdw_status
!! NAME
!!  xc_vdw_status
!!
!! FUNCTION
!!  Returns the status of the main vdW-DF switch.
!!
!! INPUTS
!!  None
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

function xc_vdw_status()

!Arguments ------------------------------------

!Local variables ------------------------------
  logical :: xc_vdw_status

! *************************************************************************

  xc_vdw_status = vdw_switch

end function xc_vdw_status
!!***

!!****f* ABINIT/m_xc_vdw/xc_vdw_trigger
!! NAME
!!  xc_vdw_trigger
!!
!! FUNCTION
!!  Switches on and off the calculation of vdW interactions.
!!
!! INPUTS
!!  condition= boolean condition to trigger the calculations
!!
!! PARENTS
!!      driver,scprqt
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine xc_vdw_trigger(condition)

!Arguments ------------------------------------
  logical, intent(in) :: condition

!Local variables ------------------------------
  character(len=512) :: msg

! *************************************************************************

  if ( my_vdw_params%functional /= 0 ) then
    if ( vdw_switch ) then
      write (msg,'(1x,a,a)') &
&       "[vdW-DF] xc_vdw_trigger: keeping xc_vdw_aggregate enabled",ch10
    else
      if ( condition ) then
        write (msg,'(1x,a,a)') &
&         "[vdW-DF] xc_vdw_trigger: enabling xc_vdw_aggregate",ch10
      else
        write (msg,'(1x,a,a)') &
&         "[vdW-DF] xc_vdw_trigger: disabling xc_vdw_aggregate",ch10
      end if

      vdw_switch = condition
    end if

    call wrtout(std_out,msg,'COLL')
  end if

end subroutine xc_vdw_trigger
!!***

!!****f* ABINIT/m_xc_vdw/xc_vdw_write
!! NAME
!!  xc_vdw_write
!!
!! FUNCTION
!!  Writes vdW-DF variables to disk.
!!
!! INPUTS
!!  filename= file to write data to
!!
!! TODO
!!  FIXME: design an extension for ETSF_IO
!!
!! PARENTS
!!      driver,vdw_kernelgen
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine xc_vdw_write(filename)

!Arguments ------------------------------------
  character(len=*), intent(in)  :: filename

!Local variables ------------------------------
  character(len=512) :: msg
  integer :: ncid,dimids(4),varids(13)
  integer :: nderivs_id,npoly_id,ndpts_id,ngpts_id,nqpts_id,nrpts_id

! *************************************************************************

  DBG_ENTER("COLL")

#if defined HAVE_NETCDF
  write(msg,'(a,1x,"Writing vdW-DF data to",1x,a)') ch10,trim(filename)
  call wrtout(std_out,msg,'COLL')

  ! Create file (overwriting any existing data)
  NETCDF_VDWXC_CHECK(nf90_create(trim(filename),NF90_CLOBBER,ncid))

  ! Write attributes
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'file_format',vdw_format_string))
  write(msg,'(a,1x,a,1x,"for",1x,a,1x,"(build",1x,a,")")') &
&   PACKAGE_NAME,ABINIT_VERSION,ABINIT_TARGET,ABINIT_VERSION_BUILD
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'generator',msg))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'acutmin',my_vdw_params%acutmin))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'aratio',my_vdw_params%aratio))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'damax',my_vdw_params%damax))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'damin',my_vdw_params%damin))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'dcut',my_vdw_params%dcut))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'dratio',my_vdw_params%dratio))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'dsoft',my_vdw_params%dsoft))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'functional',my_vdw_params%functional))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'gcut',my_vdw_params%gcut))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'nsmooth',my_vdw_params%nsmooth))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'phisoft',my_vdw_params%phisoft))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'qcut',my_vdw_params%qcut))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'qratio',my_vdw_params%qratio))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'rcut',my_vdw_params%rcut))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'rsoft',my_vdw_params%rsoft))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'tolerance',my_vdw_params%tolerance))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'tweaks',my_vdw_params%tweaks))
  NETCDF_VDWXC_CHECK(nf90_put_att(ncid,NF90_GLOBAL,'zab',my_vdw_params%zab))

  ! Define dimensions
  NETCDF_VDWXC_CHECK(nf90_def_dim(ncid,'nderivs',4,nderivs_id))
  NETCDF_VDWXC_CHECK(nf90_def_dim(ncid,'npoly',2,npoly_id))
  NETCDF_VDWXC_CHECK(nf90_def_dim(ncid,'ndpts',my_vdw_params%ndpts,ndpts_id))
  NETCDF_VDWXC_CHECK(nf90_def_dim(ncid,'ngpts',my_vdw_params%ngpts,ngpts_id))
  NETCDF_VDWXC_CHECK(nf90_def_dim(ncid,'nqpts',my_vdw_params%nqpts,nqpts_id))
  NETCDF_VDWXC_CHECK(nf90_def_dim(ncid,'nrpts',my_vdw_params%nrpts,nrpts_id))

  ! Define qmesh
  dimids(1) = nqpts_id
  NETCDF_VDWXC_CHECK(nf90_def_var(ncid,'qmesh',NF90_DOUBLE,dimids(1),varids(2)))

  ! Define dmesh
  dimids(1) = ndpts_id
  NETCDF_VDWXC_CHECK(nf90_def_var(ncid,'dmesh',NF90_DOUBLE,dimids(1),varids(3)))

  ! Define phi
  dimids(1) = ndpts_id
  dimids(2) = ndpts_id
  dimids(3) = nderivs_id
  NETCDF_VDWXC_CHECK(nf90_def_var(ncid,'phi',NF90_DOUBLE,dimids(1:3),varids(4)))

  ! Define phi_u
  dimids(1) = ndpts_id
  dimids(2) = ndpts_id
  dimids(3) = nderivs_id
  NETCDF_VDWXC_CHECK(nf90_def_var(ncid,'phi_u',NF90_DOUBLE,dimids(1:3),varids(5)))

  ! Define phi_bicubic
  dimids(1) = nderivs_id
  dimids(2) = nderivs_id
  dimids(3) = ndpts_id
  dimids(4) = ndpts_id
  NETCDF_VDWXC_CHECK(nf90_def_var(ncid,'phi_bicubic',NF90_DOUBLE,dimids(1:4),varids(6)))

  ! Define phi_u_bicubic
  dimids(1) = nderivs_id
  dimids(2) = nderivs_id
  dimids(3) = ndpts_id
  dimids(4) = ndpts_id
  NETCDF_VDWXC_CHECK(nf90_def_var(ncid,'phi_u_bicubic',NF90_DOUBLE,dimids(1:4),varids(7)))

  ! Define phir
  dimids(1) = nrpts_id
  dimids(2) = nqpts_id
  dimids(3) = nqpts_id
  NETCDF_VDWXC_CHECK(nf90_def_var(ncid,'phir',NF90_DOUBLE,dimids(1:3),varids(8)))

  ! Define phir_u
  dimids(1) = nrpts_id
  dimids(2) = nqpts_id
  dimids(3) = nqpts_id
  NETCDF_VDWXC_CHECK(nf90_def_var(ncid,'phir_u',NF90_DOUBLE,dimids(1:3),varids(9)))

  ! Define d2phidr2
  dimids(1) = nrpts_id
  dimids(2) = nqpts_id
  dimids(3) = nqpts_id
  NETCDF_VDWXC_CHECK(nf90_def_var(ncid,'d2phidr2',NF90_DOUBLE,dimids(1:3),varids(10)))

  ! Define phig
  dimids(1) = nrpts_id
  dimids(2) = nqpts_id
  dimids(3) = nqpts_id
  NETCDF_VDWXC_CHECK(nf90_def_var(ncid,'phig',NF90_DOUBLE,dimids(1:3),varids(11)))

  ! Define d2phidg2
  dimids(1) = nrpts_id
  dimids(2) = nqpts_id
  dimids(3) = nqpts_id
  NETCDF_VDWXC_CHECK(nf90_def_var(ncid,'d2phidg2',NF90_DOUBLE,dimids(1:3),varids(12)))

  ! Define qpoly_basis
  dimids(1) = nqpts_id
  dimids(2) = nqpts_id
  dimids(3) = npoly_id
  NETCDF_VDWXC_CHECK(nf90_def_var(ncid,'qpoly_basis',NF90_DOUBLE,dimids(1:3),varids(13)))

  ! Stop definitions
  NETCDF_VDWXC_CHECK(nf90_enddef(ncid))

  ! Write variables
  NETCDF_VDWXC_CHECK(nf90_put_var(ncid,varids(2),qmesh))
  NETCDF_VDWXC_CHECK(nf90_put_var(ncid,varids(3),dmesh))
  NETCDF_VDWXC_CHECK(nf90_put_var(ncid,varids(4),phi))
  NETCDF_VDWXC_CHECK(nf90_put_var(ncid,varids(5),phi_u))
  NETCDF_VDWXC_CHECK(nf90_put_var(ncid,varids(6),phi_bicubic))
  NETCDF_VDWXC_CHECK(nf90_put_var(ncid,varids(7),phi_u_bicubic))
  NETCDF_VDWXC_CHECK(nf90_put_var(ncid,varids(8),phir))
  NETCDF_VDWXC_CHECK(nf90_put_var(ncid,varids(9),phir_u))
  NETCDF_VDWXC_CHECK(nf90_put_var(ncid,varids(10),d2phidr2))
  NETCDF_VDWXC_CHECK(nf90_put_var(ncid,varids(11),phig))
  NETCDF_VDWXC_CHECK(nf90_put_var(ncid,varids(12),d2phidg2))
  NETCDF_VDWXC_CHECK(nf90_put_var(ncid,varids(13),qpoly_basis))

  ! Close file
  NETCDF_VDWXC_CHECK(nf90_close(ncid))
#else
  MSG_WARNING('writing vdW-DF variables requires NetCDF - skipping')
#endif

  DBG_EXIT("COLL")

end subroutine xc_vdw_write
!!***

! ----------------------------------------------------------------------
! Private routines
! ----------------------------------------------------------------------

!!****f* ABINIT/m_xc_vdw/vdw_df_filter
!! NAME
!!  vdw_df_filter
!!
!! FUNCTION
!!  Softens the kernel by applying a filter in reciprocal space.
!!
!! INPUTS
!!  nqpts= number of q-mesh points
!!  nrpts= number of mesh points in real space
!!  rcut= cut-off in real space
!!  gcut= cut-off in reciprocal space
!!  sofswt= Driver of the softening in real space
!!
!! OUTPUTS
!!  ngpts= number of mesh points in reciprocal space
!!
!! SIDE EFFECTS
!!  phir= the softened kernel in real space
!!  phir_u= the unsoftened kernel in real space
!!  phig= the softened kernel in reciprocal space
!!
!! PARENTS
!!      m_xc_vdw
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine vdw_df_filter(nqpts,nrpts,rcut,gcut,ngpts,sofswt)

!Arguments ------------------------------------
  integer,intent(in) :: nqpts,nrpts,sofswt
  integer,intent(out) :: ngpts
  real(dp),intent(in) :: rcut,gcut

!Local variables ------------------------------
  integer :: ig,iq1,iq2,ir
  real(dp) :: dg,dr,lstep,ptmp,q1,q2,xr,yq1,yqn,yr1,yrn
  real(dp),allocatable :: gmesh(:),rmesh(:),utmp(:)

! *************************************************************************

  DBG_ENTER("COLL")

  ! Init


  dr = rcut / nrpts
  dg = pi / rcut
  ngpts = 1 + int(gcut/dg)

  ABI_ALLOCATE(utmp,(nrpts))

  ! Create radial mesh for radial FT: src/32_util/radsintr.F90
  ABI_ALLOCATE(rmesh,(nrpts))
  forall(ir=1:nrpts) rmesh(ir) = dr * dble(ir-1)

  if ( sofswt == 1 ) then
  ! Create reciprocal radial mesh
   ABI_ALLOCATE(gmesh,(ngpts))
   forall(ig=1:ngpts) gmesh(ig) = dg * dble(ig-1)
  end if


  ! Build filtered kernel for each (q1,q2) pair
  do iq1=1,nqpts
    do iq2=1,iq1
      q1 = qmesh(iq1)
      q2 = qmesh(iq2)

      if ( sofswt == 1 ) then
      ! Build kernel in real space
      ! Note: smoothly going to zero when approaching rcut
      ! radial mesh is indeed a mesh of |\vec{r1}-\vec{r2}| values.
       do ir=1,nrpts
         xr=rmesh(ir)
         phir(ir,iq1,iq2) = vdw_df_interpolate(q1*xr,q2*xr,sofswt) * &
&          (one - ((ir - 1) / nrpts)**8)**4
       end do

      ! Obtain kernel in reciprocal space
       call radsintr(phir(:,iq1,iq2),phig(:,iq1,iq2),ngpts,nrpts,gmesh,rmesh,yq1,yqn)

      ! Filter in reciprocal space
       do ig=1,ngpts
         phig(ig,iq1,iq2) = phig(ig,iq1,iq2) * &
&          (one - ((ig - 1) / ngpts)**8)**4
       end do
       phig(ngpts+1:nrpts,iq1,iq2) = zero

      ! Go back to real space
       call radsintr(phig(:,iq1,iq2),phir(:,iq1,iq2),nrpts,ngpts,rmesh,gmesh,yr1,yrn)

      ! Calculate second derivative in real space
       d2phidr2(1,iq1,iq2) = zero
       d2phidr2(nrpts,iq1,iq2) = zero
       utmp(1) = zero
       do ir=2,nrpts-1
         ptmp = half * d2phidr2(ir-1,iq1,iq2) + two
         d2phidr2(ir,iq1,iq2) = (half - one) / ptmp
         utmp(ir) = (three * (phir(ir+1,iq1,iq2) + phir(ir-1,iq1,iq2) - &
&          two*phir(ir,iq1,iq2)) / (dr**2) - half * utmp(ir-1)) / ptmp
       end do
       do ir=nrpts-1,1,-1
         d2phidr2(ir,iq1,iq2) = d2phidr2(ir,iq1,iq2) * &
&          d2phidr2(ir+1,iq1,iq2) + utmp(ir)
       end do

      ! Calculate second derivative in reciprocal space
       d2phidg2(1,iq1,iq2) = zero
       d2phidg2(ngpts,iq1,iq2) = zero
       utmp(1) = zero
       do ig=2,ngpts-1
         ptmp = half * d2phidg2(ig-1,iq1,iq2) + two
         d2phidg2(ig,iq1,iq2) = (half - one) / ptmp
         utmp(ig) = (three * (phig(ig+1,iq1,iq2) + phig(ig-1,iq1,iq2) - &
&          two*phig(ig,iq1,iq2)) / (dr**2) - half * utmp(ig-1)) / ptmp
       end do
       do ig=ngpts-1,1,-1
         d2phidg2(ig,iq1,iq2) = d2phidg2(ig,iq1,iq2) * d2phidg2(ig+1,iq1,iq2) + &
&          utmp(ig)
       end do

      ! Symmetrize kernels & derivatives
       phir(:,iq2,iq1) = phir(:,iq1,iq2)
       phig(:,iq2,iq1) = phig(:,iq1,iq2)
       d2phidr2(:,iq2,iq1) = d2phidr2(:,iq1,iq2)
       d2phidg2(:,iq2,iq1) = d2phidg2(:,iq1,iq2)
      end if ! sofswt == 1

      if ( sofswt == 0) then
      ! Build unsoftened kernel in real space
      ! Note: smoothly going to zero when approaching rcut
      ! radial mesh is indeed a mesh of |\vec{r1}-\vec{r2}| values.
       do ir=1,nrpts
         xr=rmesh(ir)
         phir_u(ir,iq1,iq2) = vdw_df_interpolate(q1*xr,q2*xr,sofswt) * &
&          (one - ((ir - 1) / nrpts)**8)**4
       end do
      ! Symmetrize unsoftened kernel
       phir_u(:,iq2,iq1) = phir_u(:,iq1,iq2)
      end if ! sofswt == 0

     end do
   end do

  ABI_DEALLOCATE(utmp)

  DBG_EXIT("COLL")

end subroutine vdw_df_filter
!!***

!!****f* ABINIT/m_xc_vdw/vdw_df_kernel
!! NAME
!!  vdw_df_kernel
!!
!! FUNCTION
!!  Calculates the van der Waals kernel for specified d-coordinates.
!!  Decides whether to use direct integration of Eq.(14) of
!!  Dion et al., PRL 92, 246401 (2004) [[cite:Dion2004]], or to return a 4th-order
!!  polynomial for small distances.
!!
!! INPUTS
!!  d1= first coordinate
!!  d2= second coordinate
!!  dsoft= distance below which the kernel will be softened
!!  acutmin= minimum angular cut-off
!!  aratio= ratio between highest and lowest angular delta
!!  damax= maximum angular delta
!!
!! SIDE EFFECTS
!!  phisoft= value of the softened kernel at d=0
!!           will be automatically set if negative
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

function vdw_df_kernel(d1,d2,dsoft,phisoft,acutmin,aratio,damax)

!Arguments ------------------------------------
  real(dp),intent(in) :: d1,d2,dsoft,phisoft,acutmin,aratio,damax

!Local variables ------------------------------
  real(dp) :: vdw_df_kernel
  character(len=512) :: msg
  integer :: napts
  real(dp) :: deltad,dtol,dtmp,d1m,d1p,d2m,d2p,phid,phim,phip,phix

! *************************************************************************

#if ( defined DEBUG_VERBOSE ) && ( DEBUG_VERBOSE > 1 )
  DBG_ENTER("COLL")
#endif

  ! Init
  dtol = my_vdw_params%tolerance
  dtmp = sqrt(d1**2 + d2**2)

  if ( dtmp < dsoft ) then

    ! Calculate a softened version of phi when (d1,d2) -> (0,0)
    ! using a polynomial expansion for nonzero distances and an
    ! exponential fit on dsoft when d -> 0 too.
    ! Note 1: the exponential fit is inspired from the definition
    !         of phi0_soft in RS09.
    ! Note 2: contrary to RS09, deltad is proportional to dsoft.

    if ( phisoft < 0 ) then
      MSG_ERROR('the softened kernel must be positive')
    end if

    vdw_df_kernel = phisoft

    if ( dtmp >= dtol ) then
      deltad = dsoft / 100.0_dp
      d1m = (dsoft - deltad) * d1 / dtmp
      d1p = (dsoft + deltad) * d1 / dtmp
      d2m = (dsoft - deltad) * d2 / dtmp
      d2p = (dsoft + deltad) * d2 / dtmp

      phim = vdw_df_kernel_value(d1m,d2m,acutmin,aratio,damax)
      phip = vdw_df_kernel_value(d1p,d2p,acutmin,aratio,damax)
      phix = (phim + phip) / two
      phid = (phip - phim) / (two * deltad)

      vdw_df_kernel = phisoft + &
&       (four * (phix - phisoft) - phid * dsoft) &
&       * dtmp**2 / (two * dsoft**2) + &
&       (two  * (phisoft - phix) + phid * dsoft) &
&       * dtmp**4 / (two * dsoft**4)
    end if

  else

    vdw_df_kernel = vdw_df_kernel_value(d1,d2,acutmin,aratio,damax)

  end if ! dtmp < dsoft

#if ( defined DEBUG_VERBOSE ) && ( DEBUG_VERBOSE > 1 )
  DBG_EXIT("COLL")
#endif

end function vdw_df_kernel
!!***

!!****f* ABINIT/m_xc_vdw/vdw_df_kernel_value
!! NAME
!!  vdw_df_kernel_value
!!
!! FUNCTION
!!  Calculates the van der Waals kernel for specified d-coordinates. Uses
!!  direct integration of Eq.(14) of Dion et al., PRL 92, 246401 (2004) [[cite:Dion2004]].
!!
!! INPUTS
!!  d1= first coordinate
!!  d2= second coordinate
!!  acutmin= minimum angular cut-off
!!  aratio= ratio between highest and lowest angular delta
!!  damax= maximum angular delta
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

function vdw_df_kernel_value(d1,d2,acutmin,aratio,damax)

!Arguments ------------------------------------
  real(dp),intent(in) :: d1,d2,acutmin,aratio,damax

!Local variables ------------------------------
  real(dp) :: vdw_df_kernel_value
  character(len=512) :: msg
  integer :: ia,ib,napts,metyp
  real(dp) :: aa,bb,atol,atmp,acut,btmp,dain,damin,tau,ttmp,wtmp
  real(dp),allocatable :: amesh(:),amesh_cos(:),amesh_sin(:),nu1(:),nu2(:)
  real(dp),allocatable :: dphida(:),dphidb(:),dphidd(:),ycoef(:,:),work(:)
  real(dp),allocatable :: eint(:)
! *************************************************************************

#if ( defined DEBUG_VERBOSE ) && ( DEBUG_VERBOSE > 1 )
  DBG_ENTER("COLL")
#endif

  ! Create angular mesh
  !
  ! Note: Contrary to RS09, we use a linear mesh, but leave
  !       the door open to other kinds of meshes. Indeed, the use of a
  !       logarithmic mesh is interesting only when both d1 and d2 are
  !       small, in which case the kernel is usually softened.

  ! Init

!!  Testing the same way in which amesh is obtained in siesta:
!!  napts = nint(max(acutmin,aratio*sqrt(d1**2+d2**2))) This is
!!  very different to the way the number of a-mesh points are determined
!!  in Siesta.  Actually this almost corresponds to their definition of acut

  acut = max(acutmin,aratio*sqrt(d1**2+d2**2))

!! Obtaining the number of amesh points, it is introduced damin
!! and dain:
  damin = 0.01d0
  dain = min(damax,0.1*sqrt(d1**2+d2**2))
  dain = max(dain,damin)

  if ( dain == damax ) then !Linear mesh
   metyp = 2
   napts = nint( acut / damax )
  else
   bb = acut * dain / (damax-dain)
   aa = log( 1 + dain/bb )
   napts = 1 + nint( log(1+acut/bb) / aa )
   metyp = 1
  end if
!! getting napts end

  atol = my_vdw_params%tolerance

#if ( defined DEBUG_VERBOSE ) && ( DEBUG_VERBOSE > 1 )
  write(msg,'(1x,"vdw_df_kernel_value:",1x,"napts = ",i8)') napts
  call wrtout(std_out,msg,'COLL')
#endif

  ABI_ALLOCATE(amesh,(napts))
  ABI_ALLOCATE(amesh_cos,(napts))
  ABI_ALLOCATE(amesh_sin,(napts))
  ABI_ALLOCATE(dphida,(napts))
  ABI_ALLOCATE(dphidb,(napts))
  ABI_ALLOCATE(nu1,(napts))
  ABI_ALLOCATE(nu2,(napts))
  ABI_ALLOCATE(ycoef,(3,napts))
  ABI_ALLOCATE(work,(napts))
  ABI_ALLOCATE(eint,(napts))
  ! Init amesh
  !FIXME: allow for other mesh types
  !Now intruducing metyp, setting mesh_cutoff to acut, removing
  !mesh_inc=damax and including mesh_ratio=
  call vdw_df_create_mesh(amesh,napts,metyp,acut,mesh_ratio=damax/dain)

  ! Init cos(a), sin(a), nu1(a), and nu2(a)
  tau = 4 * pi / 9
!$OMP  PARALLEL DO &
!$OMP& IF ( napts > 511 ) &
!$OMP& DEFAULT(SHARED) PRIVATE(ia)
  do ia = 1,napts
    atmp = amesh(ia)
    amesh_cos(ia) = cos(atmp)
    amesh_sin(ia) = sin(atmp)
    if ( atmp == zero ) then  !it was changed the original condition: atmp < atol
      nu1(ia) = d1**2 / 2 / tau
      nu2(ia) = d2**2 / 2 / tau
    else
      if ( d1 == zero ) then  !it was changed the original condition: d1 < atol
        nu1(ia) = atmp*atmp / 2
      else
        nu1(ia) = atmp*atmp / 2 / (1 - exp(-tau*(atmp/d1)**2))
      end if
      if ( d2 == zero ) then  !it was changed the original condition: d2 < atol
        nu2(ia) = atmp*atmp / 2
      else
        nu2(ia) = atmp*atmp / 2 / (1 - exp(-tau*(atmp/d2)**2))
      end if
    end if
  end do
!$OMP END PARALLEL DO

  ! Integrate
  dphida(1) = 0
!$OMP PARALLEL DO &
!$OMP& IF ( napts > 255 ) &
!$OMP& DEFAULT(SHARED) PRIVATE(ia,ib,atmp,btmp,ttmp,wtmp)
  do ia = 2,napts
    atmp = amesh(ia)
    dphidb(1) = 0

    do ib = 2,napts
      btmp = amesh(ib)

      ! eq.(15) DRSLL04
      ttmp = half * ( 1 / (nu1(ia) + nu1(ib)) + 1 / (nu2(ia) + nu2(ib)) ) * &
&       ( 1 / (nu1(ia) + nu2(ia)) / (nu1(ib) + nu2(ib)) + &
&         1 / (nu1(ia) + nu2(ib)) / (nu2(ia) + nu1(ib)) )

      ! eq.(16) DRSLL04
      wtmp = two * ( (three - atmp*atmp) * btmp * amesh_cos(ib) * &
&       amesh_sin(ia) + (three - btmp*btmp) * atmp * amesh_cos(ia) * &
&       amesh_sin(ib) + (atmp*atmp + btmp*btmp - three) * &
&       amesh_sin(ia) * amesh_sin(ib) - &
&       three * atmp * btmp * amesh_cos(ia) * amesh_cos(ib) ) / &
&       (atmp * btmp)**3

      dphidb(ib) = atmp*atmp * btmp*btmp * wtmp * ttmp

    end do ! ib = 2,napts

!    call spline_integrate(dphida(ia),napts,damax,dphidb)
     call cspint(dphidb,amesh,napts,amesh(1),amesh(napts),ycoef,eint,work,dphida(ia))
  end do ! ia = 2,napts
!$OMP END PARALLEL DO

  ! Final value of the kernel
!  call spline_integrate(vdw_df_kernel_value,napts,damax,dphida)
     call cspint(dphida,amesh,napts,amesh(1),amesh(napts),ycoef,eint,work,vdw_df_kernel_value)
  vdw_df_kernel_value = vdw_df_kernel_value * 2 / pi**2

  ! Clean-up the mess
  ABI_DEALLOCATE(amesh)
  ABI_DEALLOCATE(amesh_cos)
  ABI_DEALLOCATE(amesh_sin)
  ABI_DEALLOCATE(nu1)
  ABI_DEALLOCATE(nu2)
  ABI_DEALLOCATE(dphida)
  ABI_DEALLOCATE(dphidb)
  ABI_DEALLOCATE(ycoef)
  ABI_DEALLOCATE(work)
  ABI_DEALLOCATE(eint)

#if ( defined DEBUG_VERBOSE ) && ( DEBUG_VERBOSE > 1 )
  DBG_EXIT("COLL")
#endif

end function vdw_df_kernel_value
!!***

!!****f* ABINIT/m_xc_vdw/vdw_df_ldaxc
!! NAME
!!  vdw_df_ldaxc
!!
!! FUNCTION
!!  Calculates LDA-based XC energy density for other vdW routines.
!!
!! INPUTS
!!  ixc= functional to apply
!!  nspden= number of spin components
!!  npts_rho= number of density points
!!  rho= density
!!  grho2= density gradient
!!
!! OUTPUT
!!  ex= exchange energy
!!  ec= correlation energy
!!  vx= exchange potential
!!  vc= correlation potential
!!
!! PARENTS
!!      m_xc_vdw
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine vdw_df_ldaxc(npts_rho,nspden,ngrad,rho_grho, &
&                       exc_lda,vxc_lda,vxcg_lda)

!Arguments ------------------------------------
  integer,intent(in) :: nspden,npts_rho,ngrad
  real(dp),intent(in) :: rho_grho(npts_rho,nspden,ngrad)
  real(dp),intent(out) :: exc_lda(2,npts_rho),vxc_lda(2,npts_rho,nspden)
  real(dp),intent(out) :: vxcg_lda(2,3,npts_rho)

!Local variables ------------------------------
  integer :: ii
  logical :: is_gga
  real(dp),allocatable :: rho_tmp(:,:),grho2(:,:)

! *************************************************************************

  DBG_ENTER("COLL")

  is_gga=libxc_functionals_isgga(vdw_funcs)

  ABI_ALLOCATE(rho_tmp,(npts_rho,nspden))
  if (is_gga) then
    ABI_ALLOCATE(grho2,(npts_rho,2*min(nspden,2)-1))
  end if

  ! Convert the quantities provided by ABINIT to the ones needed by LibXC
  !
  ! Notes:
  !
  !   * Abinit passes rho_up in the spin-unpolarized case, while LibXC
  !     expects the total density, but as we use rhonow, we don't need
  !     to convert anything.
  !
  !   * In the case off gradients:
  !
  !       - Abinit passes |grho_up|^2 while LibXC needs |grho_tot|^2;
  !       - Abinit passes |grho_up|^2, |grho_dn|^2, and |grho_tot|^2,
  !         while LibXC needs |grho_up|^2, grho_up.grho_dn,
  !         and |grho_dn|^2;
  !
  !     but again the use of rho_grho lets us build these quantities
  !     directly.
  !
  ! See ~abinit/src/56_xc/rhotoxc.F90 for details.
  if (nspden==1) then
    do ii=1,npts_rho
      rho_tmp(ii,1)=half*rho_grho(ii,1,1)
      if (is_gga) grho2(ii,1)=quarter*(rho_grho(ii,1,2)**2+rho_grho(ii,1,3)**2 &
&                                     +rho_grho(ii,1,4)**2)
    end do
  else
    do ii=1,npts_rho
      rho_tmp(ii,1)=rho_grho(ii,2,1)
      rho_tmp(ii,2)=rho_grho(ii,1,1)-rho_grho(ii,2,1)
      if (is_gga) then
        grho2(ii,1)=rho_grho(ii,2,2)**2+rho_grho(ii,2,3)**2+rho_grho(ii,2,4)**2
        grho2(ii,2)=(rho_grho(ii,1,2)-rho_grho(ii,2,2))**2 &
&                  +(rho_grho(ii,1,3)-rho_grho(ii,2,3))**2 &
&                  +(rho_grho(ii,1,4)-rho_grho(ii,2,4))**2
        grho2(ii,3)=rho_grho(ii,1,2)**2+rho_grho(ii,1,3)**2+rho_grho(ii,1,4)**2
      end if
    end do
  end if

  if (is_gga) then
    call libxc_functionals_getvxc(1,1,npts_rho,nspden,1,rho_tmp,exc_lda,vxc_lda,&
&                                 grho2=grho2,vxcgr=vxcg_lda,xc_functionals=vdw_funcs)
  else
    call libxc_functionals_getvxc(1,1,npts_rho,nspden,1,rho_tmp,exc_lda,vxc_lda,&
&                                 xc_functionals=vdw_funcs)
  end if

  ABI_DEALLOCATE(rho_tmp)
  if (is_gga) then
    ABI_DEALLOCATE(grho2)
  end if

  DBG_EXIT("COLL")

end subroutine vdw_df_ldaxc
!!***

! ----------------------------------------------------------------------
! Private utility routines
! ----------------------------------------------------------------------

!!****f* ABINIT/m_xc_vdw/vdw_df_create_mesh
!! NAME
!!  vdw_df_create_mesh
!!
!! FUNCTION
!!  Creates a 1D mesh following user specifications.
!!
!! INPUTS
!!  npts= length of the mesh
!!  mesh_type= integer representing the type of the mesh
!!  mesh_cutoff= where to cut the mesh
!!  mesh_inc(optional)= mesh increment, if linear
!!  avoid_zero(optional)= whether to avoid zero in the mesh
!!  exact_max(optional)= whether to force the last mesh element to be
!!                       strictly equal to the cutoff
!!  mesh_ratio(optional)= ratio between the first and last segment (or inverse),
!!                        when applies
!!  mesh_tolerance(optional)= tolerance to apply for the mesh
!!  mesh_file(optional)= formatted file (e20.8) where to read the mesh,
!!                       when applies
!!
!! OUTPUT
!!  mesh= the desired mesh
!!
!! PARENTS
!!      m_xc_vdw
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine vdw_df_create_mesh(mesh,npts,mesh_type,mesh_cutoff, &
& mesh_inc,mesh_ratio,mesh_tolerance,mesh_file,avoid_zero,exact_max)

!Arguments ------------------------------------
  integer,intent(in) :: npts,mesh_type
  real(dp),intent(in) :: mesh_cutoff
  real(dp),intent(out) :: mesh(npts)
  logical,optional,intent(in) :: avoid_zero,exact_max
  real(dp),optional,intent(in) :: mesh_inc,mesh_ratio,mesh_tolerance
  character(len=*),optional,intent(in) :: mesh_file

!Local variables ------------------------------
  logical :: my_avoid_zero,my_exact_max
  integer :: im,unt
  real(dp) :: exp_inc,exp_ofs,exp_tol,lin_inc,log_inc,log_ofs,log_tol
  character(len=500) :: msg

! *************************************************************************

#if ( defined DEBUG_VERBOSE ) && ( DEBUG_VERBOSE > 1 )
  DBG_ENTER("COLL")
#endif

  if ( present(avoid_zero) ) then
    my_avoid_zero = avoid_zero
  else
    my_avoid_zero = .false.
  end if

  if ( present(exact_max) ) then
    my_exact_max = exact_max
  else
    my_exact_max = .false.
  end if

#if ( defined DEBUG_VERBOSE ) && ( DEBUG_VERBOSE > 1 )
      write(msg,'(1x,a,a,1x,a,1x,i6,a,1x,a,i6,a,1x,a,1x,e10.3)') &
&       "> Mesh parameters",ch10, &
&       ">   npts    :",npts,ch10, &
&       ">   type    :",mesh_type,ch10, &
&       ">   cut-off :",mesh_cutoff
      call wrtout(std_out,msg,'COLL')
#endif

  select case(mesh_type)

    case(0)

      if ( present(mesh_file) ) then
#if ( defined DEBUG_VERBOSE ) && ( DEBUG_VERBOSE > 1 )
        write(msg,'(1x,a,1x,a)') ">   file    :",mesh_file
        call wrtout(std_out,msg,'COLL')
#endif
        if (open_file(mesh_file,msg, newunit=unt,status="old",action="read") /= 0) then
          MSG_ERROR(msg)
        end if
        read(unt,'(e20.8)') mesh
        close(unit=unt)
      else
        MSG_ERROR("  You must specify a mesh file for this mesh type!")
      end if

    case(1)

      if ( present(mesh_ratio) ) then
        if ( present(mesh_tolerance) ) then
          log_tol = mesh_tolerance
        else
          log_tol = sqrt(my_vdw_params%tolerance)
        end if
!DEBUG
!        log_inc = mesh_ratio / (npts - 2)
        log_inc = log(mesh_ratio) / (npts - 1)
!ENDDEBUG
        log_ofs = mesh_cutoff / (exp(log_inc*(npts - 1)) - 1)

!!DEBUG
!! Previous definition of log_ofs implies that the mesh
!! points are to be calculated as:
!! x(i)=log_ofs( exp(log_inc*(npts - 1)) - 1 )
!! which in turn means that x(0)=0 and is not consistent with the
!! mesh(im) points given at line 3056.
!! If the latter is adopted mesh_ratio above should be computed as:
!! mesh_ratio=log ( (x(n)-x(n-1))/(x(2)-x(1)) ) and not merely as
!! mesh_ratio= (x(n)-x(n-1))/(x(2)-x(1))
!! in order to log_inc definition makes sense.
!! The latter relation between log_inc and mesh_ratio is valid even if
!! mesh points do not start at 0:
!! x(i)=log_ofs( exp(log_inc*(npts - 1)) - 1 ) + x_0
!!ENDDEBUG

#if ( defined DEBUG_VERBOSE ) && ( DEBUG_VERBOSE > 1 )
        write(msg,'(1x,a,1x,e10.3,3(a,1x,a,1x,e10.3))') &
&         ">   ratio   :",mesh_ratio,ch10, &
&         ">   log_inc :",log_inc,ch10, &
&         ">   log_ofs :",log_ofs,ch10, &
&         ">   log_tol :",log_tol
        call wrtout(std_out,msg,'COLL')
#endif
        if ( log_inc < log_tol ) then
          MSG_WARNING("  The logarithmic increment is very low!")
        end if

        do im = 1,npts
!DEBUG
!          mesh(im) = log_ofs * exp(log_inc * (im-1)) this is not consistent
!       with definition of log_ofs at line 3025.
          mesh(im) = log_ofs * (exp(log_inc * (im-1)) - 1)
!ENDDEBUG
        end do

        !FIXME: check the correctness of 0.2
        if ( my_avoid_zero ) then
          mesh(1) = exp(log_inc * 0.2) * mesh_cutoff / &
&           (exp(log_inc*(npts - 1)) - one)
        end if

        if ( my_exact_max ) then
          mesh(npts) = mesh_cutoff
        end if
      else
        MSG_ERROR("  You must specify a mesh ratio for this mesh type!")
      end if

      ! Make sure the last point is qcut
      mesh(npts) = mesh_cutoff

    case(2)

      if ( present(mesh_inc) ) then
        lin_inc = mesh_inc
      else
        lin_inc = mesh_cutoff / npts
      end if

#if ( defined DEBUG_VERBOSE ) && ( DEBUG_VERBOSE > 1 )
      write(msg,'(1x,a,1x,e10.3)') ">   lin_inc :",lin_inc
      call wrtout(std_out,msg,'COLL')
#endif
      do im = 1,npts
        mesh(im) = lin_inc * (im - 1)
      end do

    case(3)

      if ( present(mesh_ratio) ) then
        if ( present(mesh_tolerance) ) then
          exp_tol = mesh_tolerance
        else
          exp_tol = sqrt(my_vdw_params%tolerance)
        end if
        exp_inc = mesh_ratio / (npts - 2)
        exp_ofs = mesh_cutoff * (log(exp_inc*(npts - 1)) + 1)
#if ( defined DEBUG_VERBOSE ) && ( DEBUG_VERBOSE > 1 )
        write(msg,'(1x,a,1x,e10.3,3(a,1x,a,1x,e10.3))') &
&         ">   ratio   :",mesh_ratio,ch10, &
&         ">   exp_inc :",exp_inc,ch10, &
&         ">   exp_ofs :",exp_ofs,ch10, &
&         ">   exp_tol :",exp_tol
        call wrtout(std_out,msg,'COLL')
#endif
        if ( log_inc < log_tol ) then
          MSG_WARNING("  The logarithmic increment is very low!")
        end if

        do im = 1,npts
          mesh(im) = exp_ofs / log(exp_inc * (im-1))
        end do
      else
        MSG_ERROR("  You must specify a mesh ratio for this mesh type!")
      end if

    case default

      write(msg,'(2x,a,1x,i2)') "Unknown mesh type",mesh_type
      MSG_ERROR(msg)

  end select

#if ( defined DEBUG_VERBOSE ) && ( DEBUG_VERBOSE > 1 )
  DBG_EXIT("COLL")
#endif

end subroutine vdw_df_create_mesh
!!***

!!****f* ABINIT/m_xc_vdw/vdw_df_indexof
!! NAME
!!  vdw_df_indexof
!!
!! FUNCTION
!!  Returns the index of a specified value in a sorted array.
!!
!! INPUTS
!!  list_1d= sorted array
!!  npts= length of the sorted array
!!  value= value the index of which has to be found
!!
!! NOTES
!!  The returned index is such that list_1d(index) < value < list_1d(index+1)
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

function vdw_df_indexof(list_1d,npts,value)

!Arguments ------------------------------------
  integer,intent(in) :: npts
  real(dp),intent(in) :: value,list_1d(npts)

!Local variables ------------------------------
  integer :: imin,imax,imid
  integer :: vdw_df_indexof

! *************************************************************************

  imin = 1
  imax = npts

  do
    imid = (imin + imax) / 2
    if ( value < list_1d(imid) ) then
      imax = imid
    else
      imin = imid
    end if
    if ( imin + 1 >= imax ) exit
  end do

  vdw_df_indexof = imin

end function vdw_df_indexof
!!***

!!****f* ABINIT/m_xc_vdw/vdw_df_interpolate
!! NAME
!!  vdw_df_interpolate
!!
!! FUNCTION
!!  Calculates values of Phi(d1,d2) by interpolating the kernel.
!!
!! INPUTS
!!  d1= first coordinate
!!  d2= second coordinate
!!  sofswt= switch for the kernel softening
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

function vdw_df_interpolate(d1,d2,sofswt)

  use m_errors

!Arguments ------------------------------------
  real(dp),intent(in) :: d1,d2
  integer,intent(in) :: sofswt

!Local variables ------------------------------
  real(dp) :: vdw_df_interpolate
  integer :: id1,id2,ix1,ix2,ndpts
  real(dp) :: dcut,xd(4,4)

! *************************************************************************

  vdw_df_interpolate = zero

  ndpts = size(dmesh(:))
  dcut = dmesh(ndpts)

  if ( (d1 >= dcut) .or. (d2 >= dcut) ) then
    return
  end if

  ! Find closest indices in dmesh
  id1 = vdw_df_indexof(dmesh,ndpts,d1)
  id2 = vdw_df_indexof(dmesh,ndpts,d2)

  ! Find intervals
  xd(1,:) = one
  xd(2,1) = (d1 - dmesh(id1)) / (dmesh(id1+1) - dmesh(id1))
  xd(2,2) = (d2 - dmesh(id2)) / (dmesh(id2+1) - dmesh(id2))
  xd(3,:) = xd(2,:)**2
  xd(4,:) = xd(2,:)**3


  select case ( sofswt )

    case (1)
      do ix1=1,4
        do ix2=1,4
          vdw_df_interpolate = vdw_df_interpolate + &
&           phi_bicubic(ix1,ix2,id1,id2) * xd(ix1,1) * xd(ix2,2)
        end do
      end do

    case (0)
      do ix1=1,4
        do ix2=1,4
          vdw_df_interpolate = vdw_df_interpolate + &
&           phi_u_bicubic(ix1,ix2,id1,id2) * xd(ix1,1) * xd(ix2,2)
        end do
      end do

    case default
      MSG_ERROR("  vdW-DF soft switch must be 0 or 1")

  end select

end function vdw_df_interpolate
!!***

!!****f* ABINIT/m_xc_vdw/vdw_df_internal_checks
!! NAME
!!  vdw_df_internal_checks
!!
!! FUNCTION
!!  Performs consistency checks of the vdW-DF kernel.
!!
!! INPUTS
!!  test_mode= bitfield to enable/disable tests
!!
!! PARENTS
!!      m_xc_vdw
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine vdw_df_internal_checks(test_mode)

!Arguments ------------------------------------
  integer,intent(in) :: test_mode

!Local variables ------------------------------
  character(len=512) :: msg
  integer :: ndpts,ntest,id1
  real(dp) :: acutmin,aratio,damax,dsoft,phisoft
  real(dp) :: delta_test,d1,d2
  real(dp),allocatable :: kern_test_c(:),kern_test_i(:)
  real(dp),allocatable :: intg_calc(:),intg_int(:)

! *************************************************************************

  DBG_ENTER("COLL")

  ! Here starts a test which will evaluate the integral over D for
  ! delta=0 (see Fig. 1 of DRSLL04).
  ! The integration will be performed over a finer and linear d mesh,
  ! not dmesh, for which the kernel values will be obtained by both
  ! interpolation and direct calculation. The range spanned by this mesh
  ! is half of dmesh range. Introduced variables: delta_test, kern_test_c,
  ! kern_test_i, ntest, intg_calc, intg_int.
  if ( iand(test_mode,1) /= 0 ) then

    write(msg,'(1x,a)') "[vdW-DF Check] Test #01: Kernel integral over D"
    call wrtout(std_out,msg,'COLL')

    acutmin = my_vdw_params%acutmin
    aratio = my_vdw_params%aratio
    damax = my_vdw_params%damax
    dsoft = my_vdw_params%dsoft
    ndpts = my_vdw_params%ndpts
    phisoft = my_vdw_params%phisoft

    ntest = 600
    delta_test = dmesh(ndpts) / (2 * ntest)

    ABI_ALLOCATE(kern_test_c,(ntest))
    ABI_ALLOCATE(kern_test_i,(ntest))
    ABI_ALLOCATE(intg_calc,(ntest))
    ABI_ALLOCATE(intg_int,(ntest))

#if defined DEBUG_VERBOSE
    write(msg,'(a,4x,3(1x,a10),a)') ch10, "D","Calcul.","Interp.",ch10
    call wrtout(std_out,msg,'COLL')
#endif

    do id1 = 1,ntest
      d1=id1*delta_test
      kern_test_c(id1) = (four_pi * d1**2) * &
&     vdw_df_kernel_value(d1,d1,acutmin,aratio,damax)
      kern_test_i(id1) = (four_pi * d1**2) * vdw_df_interpolate(d1,d1,0)
#if defined DEBUG_VERBOSE
      write(msg,'(4x,3(1x,e10.3))') id1*delta_test, &
&       kern_test_c(id1), kern_test_i(id1)
      call wrtout(std_out,msg,'COLL')
#endif
    end do

    call simpson_int(ntest,delta_test,kern_test_c,intg_calc)
    call simpson_int(ntest,delta_test,kern_test_i,intg_int)

    write(msg,'(a,4x,a,2(a,4x,a,e10.3))') ch10, &
&     "> Integrals over D", ch10, ">   calculated   :",intg_calc(ntest),ch10, &
&     ">   interpolated :",intg_int(ntest)
    call wrtout(std_out,msg,'COLL')
    if ( abs(intg_int(ntest) - intg_calc(ntest)) < tol3 ) then
      write(msg,'(a,1x,a)') ch10,"[vdW-DF Check] Test #01 PASSED"
    else
      write(msg,'(a,1x,a)') ch10,"[vdW-DF Check] Test #01 FAILED"
    end if
    call wrtout(std_out,msg,'COLL')

    ABI_DEALLOCATE(kern_test_c)
    ABI_DEALLOCATE(kern_test_i)
    ABI_DEALLOCATE(intg_calc)
    ABI_DEALLOCATE(intg_int)

  end if

  DBG_EXIT("COLL")

end subroutine vdw_df_internal_checks
!!***

!!****f* m_xc_vdw/vdw_df_netcdf_ioerr
!! NAME
!!  vdw_df_netcdf_ioerr
!!
!! FUNCTION
!!  Reports errors occurring when allocating memory.
!!
!! INPUTS
!!  array_name= name of the array not properly allocated.
!!  array_size= size of the array.
!!  file_name= name of the file where the allocation was performed.
!!  file_line= line number in the file where the allocation was performed.
!!
!! NOTES
!!  This routine is usually interfaced with the macros defined in abi_common.h
!!  and uses this information to define a line offset.
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine vdw_df_netcdf_ioerr(ncerr,file_name,file_line)

!Arguments ------------------------------------
 integer,intent(in) :: ncerr
 character(len=*),optional,intent(in) :: file_name
 integer,optional,intent(in) :: file_line

!Local variables-------------------------------
 character(len=1024) :: msg
 character(len=fnlen) :: my_file_name = "N/D"
 integer :: my_file_line = -1

! *************************************************************************

  if ( present(file_name) ) then
    my_file_name = trim(file_name)
  end if
  if ( present(file_line) ) then
    my_file_line = file_line
  end if

#if defined HAVE_NETCDF
  write(msg,'(a,a,3x,a)') &
&  'NetCDF returned the error:',ch10,trim(nf90_strerror(ncerr))
  if ( ncerr /= NF90_NOERR ) then
    call die(msg,trim(my_file_name),my_file_line)
  end if
#endif

end subroutine vdw_df_netcdf_ioerr
!!***

!!****f* ABINIT/m_xc_vdw/vdw_df_set_tweaks
!! NAME
!!  vdw_df_set_tweaks
!!
!! FUNCTION
!!  Expands user-specified tweaks for internal use.
!!
!! INPUTS
!!  input_tweaks= condensed tweaks
!!
!! OUTPUT
!!  output_tweaks= expanded tweaks
!!
!! PARENTS
!!      m_xc_vdw
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine vdw_df_set_tweaks(input_tweaks,output_tweaks)

!Arguments ------------------------------------
  integer,intent(in) :: input_tweaks
  type(vdw_df_tweaks_type),intent(out) :: output_tweaks

!Local variables-------------------------------

! *************************************************************************

  output_tweaks%amesh_type = iand(ishft(input_tweaks,-10),3)
  output_tweaks%deriv_type = iand(ishft(input_tweaks,-2),3)
  output_tweaks%dmesh_type = iand(ishft(input_tweaks,-6),3)
  output_tweaks%interp_type = iand(ishft(input_tweaks,-4),3)
  output_tweaks%qmesh_type = iand(ishft(input_tweaks,-8),3)
  output_tweaks%run_type = iand(input_tweaks,3)

end subroutine vdw_df_set_tweaks
!!***

!!****f* ABINIT/m_xc_vdw/vdw_df_write_func
!! NAME
!!  vdw_df_write_func
!!
!! FUNCTION
!!  Writes function names for debugging purposes.
!!
!! INPUTS
!!  func_name= name of the function to print
!!  mode= write mode (see wrtout)
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine vdw_df_write_func(func_name,mode)

!Arguments ------------------------------------
  character(len=*),intent(in) :: func_name,mode

!Local variables-------------------------------
  character(len=512) :: msg

! *************************************************************************

  write(msg,'(a,1x,a,1x,a,1x,a)') ch10,"=== subprogram :",func_name,"==="
  call wrtout(std_out,msg,mode)

end subroutine vdw_df_write_func
!!***

#endif

end module m_xc_vdw
!!***
