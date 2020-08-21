!!****m* ABINIT/m_skw
!! NAME
!!  m_skw
!!
!! FUNCTION
!!  Shankland-Koelling-Wood Fourier interpolation scheme.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG)
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

module m_skw

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_crystal
 use m_sort
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings,       only : itoa, sjoin, ktoa, yesno, ftoa
 use m_special_funcs,  only : abi_derfc
 use m_time,           only : cwtime, cwtime_report
 use m_numeric_tools,  only : imax_loc, vdiff_t, vdiff_eval, vdiff_print
 use m_bz_mesh,        only : isamek
 use m_gsphere,        only : get_irredg

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_skw/skw_t
!! NAME
!! skw_t
!!
!! FUNCTION
!!  Object implementing the Shankland-Koelling-Wood Fourier interpolation scheme.
!!  It can be used to interpolate functions in k-space with the periodicity of the
!!  reciprocal lattice and satisfying F(k) = F(Sk) for each rotation S
!!  belonging to the point group of the crystal. For readability reason,
!!  the names of the variables are chosen assuming we are interpolating electronic eigenvalues
!!  but the same object can be used to interpolate phonons as well. Just use nsppol=1 and nband = 3 * natom
!!
!! SOURCE

 type,public :: skw_t

  integer :: cplex
   ! 1 if time-reversal symmetry can be used, 2 otherwise.

  integer :: nr
   ! Number of star functions.

  integer :: nkpt
   ! Number of ab-initio k-points.

  integer :: ptg_nsym
   ! Number of operations in the point group.

  logical :: has_inversion
   ! True if the point group contains spatial inversion.

  integer :: band_block(2)
   ! Initial and final band index.

  integer :: bcount
   ! Number of bands

  integer :: nsppol
   ! Number of independent spin polarizations.

  integer,allocatable :: rpts(:,:)
   ! rpts(3, nr)
   ! Real-space lattice points (in reduced coordinates) ordered with non-decreasing length.

  integer,allocatable :: ptg_symrel(:,:,:)
    ! ptg_symrel(3,3,ptg_nsym)
    ! operations of the point group (real space).

  integer,allocatable :: ptg_symrec(:,:,:)
    ! ptg_symrec(3,3,ptg_nsym)
    ! operations of the point group (reciprocal space).

  complex(dpc),allocatable :: coefs(:,:,:)
   ! coefs(nr, bcount, nsppol).

  complex(dpc),allocatable :: cached_srk(:)
   ! cached_srk(%nr)
   ! The star function for cached_kpt (used in skw_eval_bks).
  real(dp) :: cached_kpt(3)

  complex(dpc),allocatable :: cached_srk_dk1(:,:)
   ! cached_srk_dk1(%nr, 3)
   ! The 1d derivative wrt k of the star function for cached_kpt_dk1 (used in skw_eval_bks).
  real(dp) :: cached_kpt_dk1(3)

  complex(dpc),allocatable :: cached_srk_dk2(:,:,:)
   ! cached_srk_dk2(%nr,3,3)
   ! The 2d derivatives wrt k of the star function for cached_kpt_dk2 (used in skw_eval_bks).
  real(dp) :: cached_kpt_dk2(3)

 contains

   procedure :: print => skw_print
   ! Print info about object.

   procedure :: ncwrite => skw_ncwrite
   ! Write the object in netcdf format

   procedure :: eval_bks => skw_eval_bks
   ! Interpolate eigenvalues, 1st, 2nd derivates wrt k, at an arbitrary k-point.

   procedure :: free => skw_free
   ! Free memory.

 end type skw_t
!!***

 public :: skw_new          ! Create new object.

CONTAINS  !=====================================================================================
!!***

!!****f* m_skw/skw_new
!! NAME
!!  skw_new
!!
!! FUNCTION
!!  Initialize the object.
!!
!! INPUTS
!!  cryst<crystal_t>=Crystalline structure.
!!  params(:)
!!     params(1): Ratio between star functions and ab-initio k-points.
!!     params(2:3): Activate Fourier filtering (Eq 9 of PhysRevB.61.1639 [[cite:Uehara2000]]) if params(2) > tol6
!!       params(2)=rcut, params(3) = rsigma
!!  cplex=1 if time reversal can be used, 2 otherwise.
!!  nband=Total Number of bands in the eig array.
!!  nkpt=Number of ab-initio k-points.
!!  nsppol=Number of independent spin polarizations.
!!  kpts(3,nkpt)=ab-initio k-points in reduced coordinates.
!!  eig(nband,nkpt,nsppol)=ab-initio eigenvalues.
!!  band_block(2)=Initial and final band index to interpolate. If [0,0], all bands are used
!!    This is a global variable i.e. all MPI procs MUST call the routine with the same value.
!!  comm=MPI communicator
!!
!! PARENTS
!!
!! CHILDREN
!!      sort_dp
!!
!! SOURCE

type(skw_t) function skw_new(cryst, params, cplex, nband, nkpt, nsppol, kpts, eig, band_block, comm) result(new)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nband,nkpt,nsppol,comm
 real(dp),intent(in) :: params(:)
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: band_block(2)
 real(dp),intent(in) :: kpts(3,nkpt)
 real(dp),intent(in) :: eig(nband,nkpt,nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,prtvol=1
 integer :: my_rank,nprocs,cnt,bstop,bstart,bcount,lwork
 integer :: ir,ik,ib,ii,jj,nr,band,spin,ierr,lpratio,nrwant
 real(dp),parameter :: c1=0.25_dp,c2=0.25_dp
 real(dp) :: r2,r2min,mare,mae_meV,adiff_meV,rel_err,rcut,rsigma
 real(dp) :: cpu_tot,wall_tot,gflops_tot,cpu,wall,gflops,rval
 character(len=500) :: fmt,msg
!arrays
 integer :: rmax(3)
 integer,allocatable :: ipiv(:)
 real(dp) :: list2(2)
 real(dp),allocatable :: r2vals(:),inv_rhor(:),oeig(:)
 complex(dpc),allocatable :: srk(:,:),hmat(:,:),lambda(:,:,:),work(:)

! *********************************************************************

 ABI_CHECK(nkpt > 1, sjoin("nkpt must be > 1 but got:", itoa(nkpt)))

 call cwtime(cpu_tot, wall_tot, gflops_tot, "start")

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Get slice of bands to be treated.
 new%band_block = band_block; if (all(band_block == 0)) new%band_block = [1, nband]
 bstart = new%band_block(1); bstop = new%band_block(2); bcount = bstop - bstart + 1
 new%cplex = cplex; new%nkpt = nkpt; new%nsppol = nsppol; new%bcount = bcount

 ! Get point group operations.
 call cryst%get_point_group(new%ptg_nsym, new%ptg_symrel, new%ptg_symrec, new%has_inversion, include_timrev=cplex==1)

 ! -----------------------
 ! Find nrwant star points
 ! -----------------------
 lpratio = abs(params(1))
 ABI_CHECK(lpratio > 0, "lpratio must be > 0")
 rmax = nint((one + (lpratio * new%nkpt * new%ptg_nsym) / two) ** third)
 if (new%has_inversion) then
   rmax = nint((one + (lpratio * new%nkpt * new%ptg_nsym / 2) / two) ** third)
 end if
 nrwant = lpratio * new%nkpt

 call cwtime(cpu, wall, gflops, "start")
 do
   call find_rstar_gen(new, cryst, nrwant, rmax, r2vals, comm)
   if (new%nr >= nrwant) then
     !write(std_out,*)"Entered with rmax", rmax," abs(skw%rpts(last)): ", abs(new%rpts(:,new%nr))
     exit
   end if
   write(std_out,*)"rmax: ", rmax," was not large enough to find ", nrwant," R-star points."
   rmax = 2 * rmax
   write(std_out,*)"Will try again with enlarged rmax: ",rmax
   ABI_FREE(r2vals)
 end do
 nr = new%nr
 call cwtime_report(" find_rstar_gen", cpu, wall, gflops)

 if (my_rank == master) call new%print(std_out)

 ! Compute (inverse) roughness function.
 r2min = r2vals(2)
 ABI_MALLOC(inv_rhor, (nr))
 do ir=1,nr
   r2 = r2vals(ir)
   inv_rhor(ir) = one / ((one - c1 * r2/r2min)**2 + c2 * (r2 / r2min)**3)
   ! TODO: Test the two versions.
   !if (params(1) < zero) inv_rhor(ir) = one / (c1 * r2 + c2 * r2**2)
 end do

 ! Construct star functions for the ab-initio k-points.
 ABI_MALLOC(srk, (nr, nkpt))
 do ik=1,nkpt
   call mkstar(new, kpts(:,ik), srk(:,ik))
 end do

 ! Build H(k,k') matrix (Hermitian)
 ABI_CALLOC(hmat, (nkpt-1, nkpt-1))
 cnt = 0
 do jj=1,nkpt-1
   do ii=1,jj
   !do ii=1,nkpt-1
     cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! mpi parallelism.
     do ir=2,nr
       hmat(ii, jj) = hmat(ii, jj) + &
         (srk(ir, ii) - srk(ir, nkpt)) * conjg(srk(ir, jj) - srk(ir, nkpt)) * inv_rhor(ir)
     end do
   end do
 end do
 call xmpi_sum(hmat, comm, ierr)

 ABI_MALLOC(lambda, (nkpt-1, bcount, nsppol))
 do spin=1,nsppol
   do ib=1,bcount
     band = ib + bstart - 1
     lambda(:,ib,spin) = eig(band,1:nkpt-1,spin) - eig(band,nkpt,spin)
   end do
 end do

 ! Solve all bands and spins at once [[cite:Pickett1988]]
 call wrtout(std_out, " Solving system of linear equations to get lambda coeffients (eq. 10 of PRB 38 2721)...", do_flush=.True.)
 call cwtime(cpu, wall, gflops, "start")
 ABI_MALLOC(ipiv, (nkpt-1))

 if (.False.) then
   ! General complex.
   call zgesv(nkpt-1, bcount*nsppol, hmat, nkpt-1, ipiv, lambda, nkpt-1, ierr)
   ABI_CHECK(ierr == 0, sjoin("ZGESV returned:", itoa(ierr)))
 else
   ! Hermitian version
   do ii=1,nkpt-1
     hmat(ii, ii) = real(hmat(ii, ii))
   end do
   lwork = -1
   ABI_MALLOC(work, (1))
   call zhesv("U", nkpt-1, bcount*nsppol, hmat, nkpt-1, ipiv, lambda, nkpt-1, work, lwork, ierr)
   lwork = nint(real(work(1)))
   ABI_FREE(work)
   ABI_MALLOC(work, (lwork))
   call zhesv("U", nkpt-1, bcount*nsppol, hmat, nkpt-1, ipiv, lambda, nkpt-1, work, lwork, ierr)
   ABI_CHECK(ierr == 0, sjoin("ZHESV returned:", itoa(ierr)))
   ABI_FREE(work)
 end if
 call cwtime_report(" ZHESV", cpu, wall, gflops)

 ! Compute coefficients
 ABI_MALLOC(new%coefs, (nr,bcount,nsppol))

 do spin=1,nsppol
   do ib=1,bcount
     band = ib + bstart - 1
     do ir=2,nr
       new%coefs(ir,ib,spin) = inv_rhor(ir) * dot_product(srk(ir,:nkpt-1) - srk(ir,nkpt), lambda(:nkpt-1, ib, spin))
       !new%coefs(ir,ib,spin) = inv_rhor(ir) * dot_product(lambda(:nkpt-1, ib, spin), conjg(srk(ir,:) - srk(ir,nkpt)))
       !new%coefs(ir,ib,spin) = inv_rhor(ir) * dot_product(lambda(:nkpt-1, ib, spin), conjg(srk(ir,:) - srk(ir,1)))
     end do
     new%coefs(1,ib,spin) = eig(band,nkpt,spin) - dot_product(conjg(new%coefs(2:nr, ib,spin)), srk(2:nr, nkpt))
   end do
 end do

 ! Filter high-frequency.
 if (params(2) > tol6) then
   rcut = params(2) * sqrt(r2vals(new%nr))
   rsigma = params(3); if (rsigma <= zero) rsigma = five
   call wrtout(std_out," Applying filter (Eq 9 of PhysRevB.61.1639)") ! [[cite:Uehara2000]]
   do ir=2,nr
     new%coefs(ir,:,:) = new%coefs(ir,:,:) * half * abi_derfc((sqrt(r2vals(ir)) - rcut) / rsigma)
   end do
 end if

 ! Prepare workspace arrays for star functions.
 new%cached_kpt = huge(one)
 ABI_MALLOC(new%cached_srk, (new%nr))
 new%cached_kpt_dk1 = huge(one)
 ABI_MALLOC(new%cached_srk_dk1, (new%nr, 3))
 new%cached_kpt_dk2 = huge(one)
 ABI_MALLOC(new%cached_srk_dk2, (new%nr, 3, 3))

 ABI_FREE(r2vals)
 ABI_FREE(srk)
 ABI_FREE(inv_rhor)
 ABI_FREE(hmat)
 ABI_FREE(lambda)
 ABI_FREE(ipiv)

 ! Compare ab-initio data with interpolated results.
 ABI_MALLOC(oeig, (bcount))
 fmt = sjoin("(a,", itoa(bcount), "(es12.4))")
 bstop = bstart + bcount - 1
 mare = zero; mae_meV = zero; cnt = 0
 call wrtout(std_out, ch10//" Comparing ab-initio energies with SKW interpolated results...")
 do spin=1,nsppol
   do ik=1,nkpt
     cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! mpi parallelism.

     do ib=1,bcount
       band = ib + new%band_block(1) - 1
       call new%eval_bks(band, kpts(:,ik), spin, oeig(ib))

       adiff_meV = abs(eig(band,ik,spin) - oeig(ib)); rel_err = zero
       if (abs(eig(band,ik,spin)) > tol16) rel_err = adiff_meV / abs(eig(band,ik,spin))
       rel_err = 100 * rel_err; adiff_meV = adiff_meV * Ha_meV
       mae_meV = mae_meV + adiff_meV; mare = mare + rel_err
     end do

     if (prtvol > 0) then
       ib = imax_loc(eig(bstart:bstop,ik,spin) - oeig)
       rval = (eig(bstart+ib-1,ik,spin) - oeig(ib)) * Ha_meV
       write(std_out,"(a,es12.4,2a)") &
         " SKW maxerr: ", rval, &
         " (meV), kpt: ", sjoin(ktoa(kpts(:,ik)), "band:",itoa(bstart+ib-1),", spin: ", itoa(spin))
       !write(std_out,fmt)"-- ref ", eig(bstart:bstop,ik,spin) * Ha_meV
       !write(std_out,fmt)"-- int ", oeig * Ha_meV
       !call vdiff_print(vdiff_eval(1, bcount, eig(bstart:bstop,ik,spin), oeig, one))
     end if
   end do
 end do
 ABI_FREE(oeig)

 ! Issue warning if error too large.
 list2 = [mare, mae_meV]; call xmpi_sum(list2, comm, ierr); mare = list2(1); mae_meV = list2(2)
 cnt = bcount * nkpt * nsppol; mare = mare / cnt; mae_meV = mae_meV / cnt
 write(std_out,"(2(a,es12.4),a,/)")" MARE: ",mare, ", MAE: ", mae_meV, " (meV)"
 if (mae_meV > ten) then
   write(msg,"(2a,2(a,es12.4),a)") &
     "Large error in SKW interpolation!",ch10," MARE: ",mare, ", MAE: ", mae_meV, " (meV)"
   call wrtout(ab_out, msg)
   MSG_WARNING(msg)
 end if

 call cwtime_report(" skw_new", cpu_tot, wall_tot, gflops_tot, end_str=ch10)

end function skw_new
!!***

!----------------------------------------------------------------------

!!****f* m_skw/skw_print
!! NAME
!!  skw_print
!!
!! FUNCTION
!!  Print info on object
!!
!! INPUTS
!!  unt=Fortran unit number.
!!
!! OUTPUT
!!  only writing
!!
!! PARENTS
!!
!! CHILDREN
!!      get_irredg,sort_dp,xmpi_allgatherv,xmpi_split_work2_i4b,xmpi_sum
!!
!! SOURCE

subroutine skw_print(skw, unt)

!Arguments ------------------------------------
!scalars
 class(skw_t),intent(in) :: skw
 integer,intent(in) :: unt

! *********************************************************************

 write(unt,"(a)")" === Shankland-Koelling-Wood Fourier interpolation scheme ==="
 write(unt,"(a)")sjoin(" nsppol", itoa(skw%nsppol), ", cplex:", itoa(skw%cplex))
 write(unt,"(a)")sjoin(" Number of ab-initio k-points:", itoa(skw%nkpt))
 write(unt,"(a)")sjoin(" Number of star functions:", itoa(skw%nr))
 write(unt,"(a)")sjoin(" Stars/Nk ratio:", ftoa(skw%nr * one / skw%nkpt))
 write(unt,"(a)")sjoin(" Has spatial inversion:", yesno(skw%has_inversion))

end subroutine skw_print
!!***

!----------------------------------------------------------------------

!!****f* m_skw/skw_ncwrite
!! NAME
!! skw_ncwrite
!!
!! FUNCTION
!!   Write the object in netcdf format
!!
!! INPUTS
!!  ncid=NC file handle.
!!  [prefix]=String prepended to netcdf dimensions/variables (HDF5 poor-man groups)
!!   "skw" if not specified.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function skw_ncwrite(self, ncid, prefix) result(ncerr)

!Arguments ------------------------------------
!scalars
 class(skw_t),intent(in) :: self
 integer,intent(in) :: ncid
 character(len=*),optional,intent(in) :: prefix

#ifdef HAVE_NETCDF
!Local variables-------------------------------
!scalars
 character(len=500) :: prefix_
!arrays
 real(dp),allocatable :: real_coefs(:,:,:,:)

! *************************************************************************

 prefix_ = "skw"; if (present(prefix)) prefix_ = trim(prefix)

 ! Define dimensions.
 ncerr = nctk_def_dims(ncid, [ &
   nctkdim_t("nr", self%nr), nctkdim_t("nkpt", self%nkpt), nctkdim_t("bcount", self%bcount), &
   nctkdim_t("nsppol", self%nsppol)], &
   defmode=.True., prefix=prefix_)
 NCF_CHECK(ncerr)

 ncerr = nctk_def_arrays(ncid, [ &
  ! Atomic structure and symmetry operations
  nctkarr_t("rpts", "dp", "three, number_of_cartesian_directions, number_of_vectors"), &
  nctkarr_t("kpts", "dp", "three, nkpt"), &
  nctkarr_t("coefs", "dp", "two, nr, bcount, nsppol") &
 ], prefix=prefix_)
 NCF_CHECK(ncerr)

 ! Write data.
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("rpts")), self%rpts))
 !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("kpts")), self%kpts))
 ABI_MALLOC(real_coefs, (2, self%nr, self%bcount, self%nsppol))
 real_coefs(1,:,:,:) = real(self%coefs); real_coefs(2,:,:,:) = aimag(self%coefs)
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("coefs")), real_coefs))
 ABI_FREE(real_coefs)

contains
  pure function pre(istr) result(ostr)
    character(len=*),intent(in) :: istr
    character(len=len_trim(prefix_) + len_trim(istr)+1) :: ostr
    ostr = trim(prefix_) // trim(istr)
  end function pre

#endif

end function skw_ncwrite
!!***

!!****f* m_skw/skw_eval_bks
!! NAME
!!  skw_eval_bks
!!
!! FUNCTION
!!  Interpolate the energies for an arbitrary k-point and spin with slow FT.
!!
!! INPUTS
!!  band=Band index (global index associated to the input eigenvalues, i.e. independent of band_block)
!!  kpt(3)=K-point in reduced coordinates.
!!  spin=Spin index.
!!
!! OUTPUT
!!  oeig=interpolated eigenvalues
!!    Note that oeig is not necessarily sorted in ascending order.
!!    The routine does not reorder the interpolated eigenvalues
!!    to be consistent with the interpolation of the derivatives.
!!  [oder1(3)]=First-order derivatives wrt k in reduced coordinates.
!!  [oder2(3,3)]=Second-order derivatives wrt k in reduced coordinates.
!!
!! PARENTS
!!
!! CHILDREN
!!      get_irredg,sort_dp,xmpi_allgatherv,xmpi_split_work2_i4b,xmpi_sum
!!
!! SOURCE

subroutine skw_eval_bks(skw, band, kpt, spin, oeig, oder1, oder2)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,spin
 class(skw_t),intent(inout) :: skw
!arrays
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(out) :: oeig
 real(dp),optional,intent(out) :: oder1(3),oder2(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,ib
! *********************************************************************

 ib = band - skw%band_block(1) + 1
 ABI_CHECK(ib >= 1 .and. ib <= skw%bcount, sjoin("out of range band:", itoa(band)))

 ! Compute star function for this k-point (if not already in memory)
 if (any(kpt /= skw%cached_kpt)) then
   call mkstar(skw, kpt, skw%cached_srk)
   skw%cached_kpt = kpt
 end if

 oeig = dot_product(conjg(skw%coefs(:,ib,spin)), skw%cached_srk)

 ! TODO: Test Derivatives
 if (present(oder1)) then
   ! Compute first-order derivatives.
   if (any(kpt /= skw%cached_kpt_dk1)) then
     call mkstar_dk1(skw, kpt, skw%cached_srk_dk1)
     skw%cached_kpt_dk1 = kpt
   end if

   do ii=1,3
     oder1(ii) = dot_product(conjg(skw%coefs(:,ib,spin)), skw%cached_srk_dk1(:,ii)) * two_pi
   end do
 end if

 if (present(oder2)) then
   ! Compute second-order derivatives.
   if (any(kpt /= skw%cached_kpt_dk2)) then
     call mkstar_dk2(skw, kpt, skw%cached_srk_dk2)
     skw%cached_kpt_dk2 = kpt
   end if

   oder2 = zero
   do jj=1,3
     do ii=1,jj
       oder2(ii, jj) = dot_product(conjg(skw%coefs(:,ib,spin)), skw%cached_srk_dk2(:,ii,jj)) * two_pi**2
       if (ii /= jj) oder2(jj, ii) = oder2(ii, jj)
     end do
   end do
 end if

end subroutine skw_eval_bks
!!***

!----------------------------------------------------------------------

!!****f* m_skw/skw_eval_fft
!! NAME
!!  skw_eval_fft
!!
!! FUNCTION
!!  Interpolate the energies for an arbitrary k-point and spin with slow FT.
!!
!! INPUTS
!!  cryst<crystal_t>=Crystalline structure.
!!  nfft=Number of points in FFT mesh.
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  band=Band index.
!!  spin=Spin index.
!!
!! OUTPUT
!!  oeig_mesh(nfft)=interpolated eigenvalues
!!    Note that oeig is not necessarily sorted in ascending order.
!!    The routine does not reorder the interpolated eigenvalues
!!    to be consistent with the interpolation of the derivatives.
!!  [oder1_mesh(3,nfft))]=First-order derivatives wrt k in reduced coordinates.
!!  [oder2_mesh(3,3,nfft)]=Second-order derivatives wrt k in reduced coordinates.
!!
!! PARENTS
!!
!! CHILDREN
!!      get_irredg,sort_dp,xmpi_allgatherv,xmpi_split_work2_i4b,xmpi_sum
!!
!! SOURCE

subroutine skw_eval_fft(skw, ngfft, nfft, band, spin, oeig_mesh, oder1_mesh, oder2_mesh)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,band,spin
 type(skw_t),intent(in) :: skw
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(out) :: oeig_mesh(nfft)
 real(dp),optional,intent(out) :: oder1_mesh(3,nfft)
 real(dp),optional,intent(out) :: oder2_mesh(3,3,nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourdp0=0,paral_kgb0=0
 integer :: cplex,ix,iy,iz,nx,ny,nz,ldx,ldy,ldz,ifft,ir,ii,jj
!arrays
 real(dp),allocatable :: fofg(:,:),fofr(:),workg(:,:)
! *********************************************************************

 ! Internal MPI_type needed for calling fourdp!
 !call initmpi_seq(skw%mpi_enreg)
 ! Initialize tables to call fourdp in sequential
 !call init_distribfft_seq(skw%mpi_enreg%distribfft,'c',ngfft(2),ngfft(3),'all')
 !call init_distribfft_seq(skw%mpi_enreg%distribfft,'f',ngfft(2),ngfft(3),'all')

 ! Transfer data from the G-sphere to the FFT box.
 ! Use the following indexing (N means ngfft of the adequate direction)
 ! 0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= gc
 ! 1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index ig
 nx = ngfft(1); ny = ngfft(2); nz = ngfft(3)
 ldx = ngfft(4); ldy = ngfft(5); ldz = ngfft(6)

 cplex = 2
 ABI_MALLOC(fofr, (cplex*nfft))
 ABI_MALLOC(fofg, (2,nfft))

 ! TODO: Complete Derivatives
 ! Decide between fourdp and fftbox API
 fofg = zero
 do ir=1,skw%nr
   ix = skw%rpts(1,ir); if (ix<0) ix=ix+nx; ix=ix+1
   iy = skw%rpts(2,ir); if (iy<0) iy=iy+ny; iy=iy+1
   iz = skw%rpts(3,ir); if (iz<0) iz=iz+nz; iz=iz+1
   ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
   !band = ib + bstart - 1
   fofg(1,ifft) = real(skw%coefs(ir, band, spin))
   fofg(2,ifft) = aimag(skw%coefs(ir, band, spin))
 end do

 !call fourdp(cplex, fofg, fofr, +1, skw%mpi_enreg, nfft, ngfft, paral_kgb0, tim_fourdp0)
 !fofr = fofr / skw%ptg_nsym
 if (cplex == 1) oeig_mesh = fofr
 if (cplex == 2) oeig_mesh = fofr(1::2)

 if (present(oder1_mesh)) then
   ! Compute first-order derivatives.
   ABI_MALLOC(workg, (2,nfft))
   do ii=1,3
      workg = zero
      do ir=1,skw%nr
        !ifft
        workg(1,ifft) = -fofg(2,ifft) * skw%rpts(ii,ir)
        workg(2,ifft) =  fofg(1,ifft) * skw%rpts(ii,ir)
      end do
      !call fourdp(cplex, workg, fofr, +1, skw%mpi_enreg, nfft, ngfft, paral_kgb0, tim_fourdp0)
      if (cplex == 1) oder1_mesh(ii,:) = fofr
      if (cplex == 2) oder1_mesh(ii,:) = fofr(1::2)
   end do
   ABI_FREE(workg)
 end if

 if (present(oder2_mesh)) then
   ! Compute second-order derivatives.
   ABI_MALLOC(workg, (2,nfft))
   do jj=1,3
     do ii=1,jj
       workg = zero
       do ir=1,skw%nr
         !ifft
         workg(:,ifft) = -fofg(:,ifft) * skw%rpts(ii,ir) * skw%rpts(jj,ir)
       end do
       !call fourdp(cplex, workg, fofr, +1, skw%mpi_enreg, nfft, ngfft, paral_kgb0, tim_fourdp0)
       if (cplex == 1) oder2_mesh(ii,jj,:) = fofr
       if (cplex == 2) oder2_mesh(ii,jj,:) = fofr(1::2)
       if (ii /= jj) oder2_mesh(jj, ii,:) = oder2_mesh(ii, jj,:)
     end do
   end do
   ABI_FREE(workg)
 end if

 ABI_FREE(fofg)
 ABI_FREE(fofr)

end subroutine skw_eval_fft
!!***

!----------------------------------------------------------------------

!!****f* m_skw/skw_free
!! NAME
!!  skw_free
!!
!! FUNCTION
!!  Free memory
!!
!! PARENTS
!!
!! CHILDREN
!!      get_irredg,sort_dp,xmpi_allgatherv,xmpi_split_work2_i4b,xmpi_sum
!!
!! SOURCE

subroutine skw_free(skw)

!Arguments ------------------------------------
!scalars
 class(skw_t),intent(inout) :: skw

! *********************************************************************

 ABI_SFREE(skw%rpts)
 ABI_SFREE(skw%ptg_symrel)
 ABI_SFREE(skw%ptg_symrec)
 ABI_SFREE(skw%coefs)

 ABI_SFREE(skw%cached_srk)
 skw%cached_kpt = huge(one)
 ABI_SFREE(skw%cached_srk_dk1)
 skw%cached_kpt_dk1 = huge(one)
 ABI_SFREE(skw%cached_srk_dk2)
 skw%cached_kpt_dk2 = huge(one)

end subroutine skw_free
!!***

!----------------------------------------------------------------------

!!****f* m_skw/mkstar
!! NAME
!!  mkstar
!!
!! FUNCTION
!!  Compute the star function for k-point kpt
!!
!! INPUTS
!!  kpt(3)=K-point in reduced coordinates.
!!
!! OUTPUT
!!  srk(%nr)=Star function for this k-point.
!!
!! PARENTS
!!      m_skw
!!
!! CHILDREN
!!      get_irredg,sort_dp,xmpi_allgatherv,xmpi_split_work2_i4b,xmpi_sum
!!
!! SOURCE

subroutine mkstar(skw, kpt, srk)

!Arguments ------------------------------------
!scalars
 type(skw_t),intent(in) :: skw
!arrays
 real(dp),intent(in) :: kpt(3)
 complex(dpc),intent(out) :: srk(skw%nr)

!Local variables-------------------------------
!scalars
 integer :: ir,isym
!arrays
 real(dp) :: sk(3)

! *********************************************************************

 srk = zero
 do isym=1,skw%ptg_nsym
   sk = two_pi * matmul(transpose(skw%ptg_symrel(:,:,isym)), kpt)
   do ir=1,skw%nr
     srk(ir) = srk(ir) + exp(j_dpc * dot_product(sk, skw%rpts(:,ir)))
   end do
 end do
 srk = srk / skw%ptg_nsym

end subroutine mkstar
!!***

!----------------------------------------------------------------------

!!****f* m_skw/mkstar_dk1
!! NAME
!!  mkstar_dk1
!!
!! FUNCTION
!!  Compute the 1st derivative of the star function wrt k
!!
!! INPUTS
!!  kpt(3)=K-point in reduced coordinates.
!!
!! OUTPUT
!!  srk_dk1(%nr,3)=Derivative of the star function wrt k in reduced coordinates.
!!
!! PARENTS
!!      m_skw
!!
!! CHILDREN
!!      get_irredg,sort_dp,xmpi_allgatherv,xmpi_split_work2_i4b,xmpi_sum
!!
!! SOURCE

subroutine mkstar_dk1(skw, kpt, srk_dk1)

!Arguments ------------------------------------
!scalars
 type(skw_t),intent(in) :: skw
!arrays
 real(dp),intent(in) :: kpt(3)
 complex(dpc),intent(out) :: srk_dk1(skw%nr,3)

!Local variables-------------------------------
!scalars
 integer :: ir,isym
!arrays
 real(dp) :: sk(3)
 complex(dpc) :: work(3,skw%nr)

! *********************************************************************

 work = zero
 do isym=1,skw%ptg_nsym
   sk = two_pi * matmul(transpose(skw%ptg_symrel(:,:,isym)), kpt)
   do ir=1,skw%nr
     work(:,ir) = work(:,ir) + exp(j_dpc * dot_product(sk, skw%rpts(:,ir))) * &
        matmul(skw%ptg_symrel(:,:,isym), skw%rpts(:,ir))
   end do
 end do
 work = j_dpc * work / skw%ptg_nsym
 srk_dk1 = transpose(work)

end subroutine mkstar_dk1
!!***

!----------------------------------------------------------------------

!!****f* m_skw/mkstar_dk2
!! NAME
!!  mkstar_dk2
!!
!! FUNCTION
!!  Compute the 2st derivatives of the star function wrt k
!!
!! INPUTS
!!  kpt(3)=K-point in reduced coordinates.
!!
!! OUTPUT
!!  srk_dk2(%nr,3,3)=2nd derivatives of the star function wrt k in reduced coordinates.
!!
!! PARENTS
!!      m_skw
!!
!! CHILDREN
!!      get_irredg,sort_dp,xmpi_allgatherv,xmpi_split_work2_i4b,xmpi_sum
!!
!! SOURCE

subroutine mkstar_dk2(skw, kpt, srk_dk2)

!Arguments ------------------------------------
!scalars
 type(skw_t),intent(in) :: skw
!arrays
 real(dp),intent(in) :: kpt(3)
 complex(dpc),intent(out) :: srk_dk2(skw%nr,3,3)

!Local variables-------------------------------
!scalars
 integer :: ir,isym,ii,jj
 complex(dpc) :: eiskr
!arrays
 integer :: sr(3)
 real(dp) :: sk(3)
 complex(dpc) :: work(3,3,skw%nr)

! *********************************************************************

 work = zero
 do isym=1,skw%ptg_nsym
   sk = two_pi * matmul(transpose(skw%ptg_symrel(:,:,isym)), kpt)
   do ir=1,skw%nr
     sr = matmul(skw%ptg_symrel(:,:,isym), skw%rpts(:,ir))
     eiskr = exp(j_dpc * dot_product(sk, skw%rpts(:,ir)))
     do jj=1,3
       do ii=1,jj
         work(ii,jj,ir) = work(ii,jj,ir) + eiskr * sr(ii) * sr(jj)
       end do
     end do
   end do
 end do
 work = - work / skw%ptg_nsym

 do jj=1,3
   do ii=1,jj
     srk_dk2(:, ii, jj) = work(ii, jj, :)
     if (ii /= jj) srk_dk2(:,jj,ii) = work(:,ii,jj)
   end do
 end do

end subroutine mkstar_dk2
!!***

!----------------------------------------------------------------------

!!****f* m_skw/find_rstar_gen
!! NAME
!!  find_rstar_gen
!!
!! FUNCTION
!!  Find the R-space points generating the stars.
!!  Set skw%nr and skw%rpts.
!!
!! INPUTS
!!  cryst<crystal_t>=Crystalline structure.
!!  nrwant=Number of R-space points wanted
!!  rmax(3)=Max reduced components of supercell.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  or2vals(skw%nr)=||R||**2
!!
!! PARENTS
!!      m_skw
!!
!! CHILDREN
!!      get_irredg,sort_dp,xmpi_allgatherv,xmpi_split_work2_i4b,xmpi_sum
!!
!! SOURCE

subroutine find_rstar_gen(skw, cryst, nrwant, rmax, or2vals, comm)

!Arguments ------------------------------------
!scalars
 type(skw_t),intent(inout) :: skw
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: nrwant,comm
!arrays
 integer,intent(in) :: rmax(3)
 real(dp),allocatable,intent(out) :: or2vals(:)

!Local variables-------------------------------
!scalars
 integer :: cnt,nstars,i1,i2,i3,msize,ir,nsh,ish,ss,ee,nst,ierr,nprocs,my_rank,ii
 real(dp) :: r2_prev
 !character(len=500) :: msg
!arrays
 integer,allocatable :: iperm(:),rtmp(:,:),rgen(:,:),r2sh(:),shlim(:),sh_start(:),sh_stop(:)
 integer,allocatable :: recvcounts(:),displs(:),recvbuf(:,:)
 real(dp),allocatable :: r2tmp(:),cnorm(:)

! *********************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 msize = product(2*rmax + 1)
 ABI_MALLOC(rtmp, (3, msize))
 ABI_MALLOC(r2tmp, (msize))

 cnt = 0
 do i3=-rmax(3),rmax(3)
   do i2=-rmax(2),rmax(2)
     do i1=-rmax(1),rmax(1)
       cnt = cnt + 1
       rtmp(:, cnt) = [i1,i2,i3]
       r2tmp(cnt) = dot_product(rtmp(:,cnt), matmul(cryst%rmet, rtmp(:,cnt)))
     end do
   end do
 end do

 ! Sort r2tmp
 ABI_MALLOC(iperm, (msize))
 iperm = [(i1, i1=1,msize)]
 call sort_dp(msize, r2tmp, iperm, tol12)

 ! Find R-points generating the stars.
 ABI_MALLOC(rgen, (3, msize))
 do ir=1,msize
   rgen(:,ir) = rtmp(:,iperm(ir))
 end do
 rtmp = rgen
 ABI_FREE(iperm)

 ABI_MALLOC(r2sh, (msize))     ! Correspondence between R and the shell index.
 ABI_MALLOC(shlim, (msize+1))  ! For each shell, the index of the initial G-vector.
 nsh = 1; r2sh(1) = 1; shlim(1) = 1; r2_prev = zero
 do ir=2,msize
   if (abs(r2tmp(ir) - r2_prev) > r2tmp(ir) * tol8) then
     r2_prev = r2tmp(ir); nsh = nsh + 1; shlim(nsh) = ir
     !write(std_out,*)"nsh: ",shlim(nsh) - shlim(nsh-1)
   end if
   r2sh(ir) = nsh
 end do
 shlim(nsh+1) = msize + 1
 ABI_FREE(r2tmp)
 ABI_FREE(r2sh)

 !call get_irredg(msize, skw%ptg_nsym, +1, cryst%rprimd, skw%ptg_symrel, rtmp, nstars, rgen, cnorm)
 !write(66,*)nstars; do ish=1,nstars; write(66,*)rgen(:,ish); end do

 ! Distribute shells among processor so that we can parallelize the search algorithm.
 ! Each proc works on a contigous block of shells, then we have to gather the results.
 ABI_MALLOC(sh_start, (0:nprocs-1))
 ABI_MALLOC(sh_stop, (0:nprocs-1))
 call xmpi_split_work2_i4b(nsh, nprocs, sh_start, sh_stop)

 ABI_MALLOC(cnorm, (msize))
 nstars = 0
 do ish=sh_start(my_rank),sh_stop(my_rank)
   ss = shlim(ish); ee = shlim(ish+1) - 1; msize = ee - ss + 1
   call get_irredg(msize, skw%ptg_nsym, + 1, cryst%rprimd, skw%ptg_symrel, rtmp(:,ss:), &
     nst, rgen(:,nstars+1:), cnorm(nstars+1:))
   nstars = nstars + nst
 end do

 ABI_FREE(cnorm)
 ABI_FREE(sh_start)
 ABI_FREE(sh_stop)
 ABI_FREE(rtmp)
 ABI_FREE(shlim)

 if (nprocs > 1) then
   ! Collect star functions.
   ABI_MALLOC(recvcounts, (nprocs))
   recvcounts = 0; recvcounts(my_rank+1) = 3 * nstars
   call xmpi_sum(recvcounts, comm, ierr)
   ABI_MALLOC(displs, (nprocs))
   displs(1) = 0
   do ii=2,nprocs
     displs(ii) = sum(recvcounts(:ii-1))
   end do
   call xmpi_sum(nstars, nst, comm, ierr)   ! Now nst is the total number of star functions.
   ABI_MALLOC(recvbuf, (3, nst))
   call xmpi_allgatherv(rgen, 3*nstars, recvbuf, recvcounts, displs, comm, ierr)
   ABI_FREE(recvcounts)
   ABI_FREE(displs)
   nstars = nst
   rgen(:,1:nstars) = recvbuf
   ABI_FREE(recvbuf)
 end if
 !if (my_rank == 0) then
 !  write(67,*)"nstars",nstars,"nsh",nsh; do ish=1,nstars; write(67,*)rgen(:,ish); end do
 !end if

 ! Store rpts and compute ||R||**2.
 skw%nr = min(nstars, nrwant)
 if (allocated(skw%rpts)) then
   ABI_FREE(skw%rpts)
 end if
 ABI_MALLOC(skw%rpts, (3, skw%nr))
 skw%rpts = rgen(:,1:skw%nr)
 ABI_MALLOC(or2vals, (skw%nr))
 do ir=1,skw%nr
   or2vals(ir) = dot_product(skw%rpts(:,ir), matmul(cryst%rmet, skw%rpts(:,ir)))
 end do

 ABI_FREE(rgen)

end subroutine find_rstar_gen
!!***

end module m_skw
!!***
