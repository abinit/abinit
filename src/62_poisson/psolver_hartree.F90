!{\src2tex{textfont=tt}}
!!****f* ABINIT/Psolver_hartree
!! NAME
!! Psolver_hartree
!!
!! FUNCTION
!! Given rho(r), compute Hartree potential considering the system as
!! an isolated one. This potential is obtained from the convolution
!! of 1/r and rho(r), treated in Fourier space. This method is a wrapper around
!! Psolver() developped for BigDFT.
!! It does not compute the xc energy nor potential. See psolver_rhohxc() to do it.
!! WARNING : the XC energy and potential computation capability has been
!! for spin-polarized case, as everything is done as if nspden=1
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mpi_enreg=MPI-parallelisation information.
!!  rhor(nfft,nspden)=electron density in real space in electrons/bohr**3
!!
!! OUTPUT
!!  enhartr=returned Hartree energy (hartree).
!!  vhartr(nfft)=Hartree potential.
!!
!! NOTE
!!  In PSolver, with nspden == 2, rhor(:,1) = density up and
!!                                rhor(:,2) = density down.
!!  But in ABINIT (dtset%usewvl != 1) rhor(:,1) = total density and
!!                                    rhor(:,2) = density up .
!!  In ABINIT (dtset%usewvl != 1), the same convention is used as in PSolver.
!!
!! PARENTS
!!      mklocl_realspace,nres2vres
!!
!! CHILDREN
!!      h_potential,psolver_kernel,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine psolver_hartree(enhartr, hgrid, icoulomb, me, mpi_comm, nfft, ngfft, nproc, &
     & nscforder, nspden, rhor, vhartr, usewvl)

 use defs_basis
 use m_errors
 use m_profiling_abi
#if defined HAVE_DFT_BIGDFT
 use BigDFT_API,     only : coulomb_operator
 use poisson_solver, only : H_potential
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psolver_hartree'
 use interfaces_14_hidewrite
 use interfaces_62_poisson, except_this_one => psolver_hartree
!End of the abilint section

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer, intent(in)           :: nfft, nspden, icoulomb, usewvl, mpi_comm, me, nproc, nscforder
  real(dp), intent(out)         :: enhartr
  !arrays
  integer, intent(in)    :: ngfft(3)
  real(dp),intent(in)    :: hgrid(3)
  real(dp),intent(in)    :: rhor(nfft,nspden)
  real(dp),intent(out)   :: vhartr(nfft)

  !Local variables-------------------------------
#if defined HAVE_DFT_BIGDFT
  !scalars
  character(len=500) :: message
  character(len = 1) :: datacode, bndcode
  !arrays
  real(dp), dimension(1) :: pot_ion_dummy
  type(coulomb_operator):: kernel
#endif

! *********************************************************************

#if defined HAVE_DFT_BIGDFT

 if (icoulomb == 0) then
!  The kernel is built with 'P'eriodic boundary counditions.
   bndcode = 'P'
 else if (icoulomb == 1) then
!  The kernel is built with 'F'ree boundary counditions.
   bndcode = 'F'
 else if (icoulomb == 2) then
!  The kernel is built with 'S'urface boundary counditions.
   bndcode = 'S'
 end if

 if(nspden > 2 .and. usewvl/=0 )then
   write(message, '(a,a,a,i0)' )&
&   'Only non-spin-polarised or collinear spin is allowed for wavelets,',ch10,&
&   'while the argument nspden = ', nspden
   MSG_BUG(message)
 end if

!We do the computation.
 write(message, "(A,A,A,3I6)") "Psolver_hartree(): compute potential (Vhartree)...", ch10, &
& " | dimension:", ngfft(1:3)
 call wrtout(std_out, message,'COLL')

 if (usewvl == 0) then
   vhartr(:)  = rhor(:, 1)

   datacode = 'G'
!  This may not work with MPI in the planewave code...
 else
   if(nspden==1)vhartr(:)  = rhor(:, 1)
   if(nspden==2)vhartr(:)  = rhor(:, 1) + rhor(:, 2)
!  The data are 'D'istributed in the wavelet case or 'G'lobal otherwise.
   if (nproc > 1) then
     datacode = 'D'
   else
     datacode = 'G'
   end if
 end if

!We get the kernel.
 call psolver_kernel( hgrid, 2, icoulomb, me, kernel, mpi_comm, ngfft, nproc, nscforder)


!We attack PSolver with the total density contained in vhartr.
!This is also valid for spin-polarized (collinear and non-collinear)
!systems. Thus we enter nspden (last arg of PSolver) as being 1.
!Warning : enxc and evxc are meaningless.
! call psolver(bndcode, datacode, me, nproc, ngfft(1), ngfft(2), ngfft(3),&
!& 0, hgrid(1), hgrid(2), hgrid(3), vhartr, kernel%co%kernel, pot_ion_dummy, &
!& enhartr, enxc, evxc, 0.d0, .false., 1)

 call H_potential(datacode,kernel,vhartr,pot_ion_dummy,&
& enhartr,zero,.false.)


#else
 BIGDFT_NOTENABLED_ERROR() 
 if (.false.) write(std_out,*)  nfft,nspden,icoulomb,usewvl,mpi_comm,me,nproc,nscforder,enhartr,&
& ngfft(1),hgrid(1),rhor(1,1),vhartr(1)
#endif

end subroutine psolver_hartree
!!***
