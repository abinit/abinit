!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_psitohpsi
!! NAME
!! wvl_psitohpsi
!!
!! FUNCTION
!! Compute new trial potential and calculate the hamiltionian application into hpsi.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/dpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mpi_enreg=information about MPI parallelization
!!
!! OUTPUT
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree)
!!  vtrial(nfft,nspden)=new potential
!!
!! NOTES
!!
!! PARENTS
!!      afterscfloop,rhotov,setvtr,vtorho
!!
!! CHILDREN
!!      psitohpsi,total_energies,wrtout,wvl_vtrial_abi2big,wvl_vxc_abi2big
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wvl_psitohpsi(alphamix,eexctX, eexcu, ehart, ekin_sum, epot_sum, eproj_sum, eSIC_DC, &
     & itrp, iter, iscf, me, natom, nfft, nproc, nspden, rpnrm, scf, &
     & vexcu, wvl, wvlbigdft, xcart, xcstr,vtrial,vxc)

 use defs_basis
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
 use m_abi2big, only: wvl_vxc_abi2big,wvl_vtrial_abi2big
#if defined HAVE_BIGDFT
 use BigDFT_API, only: psitohpsi, KS_POTENTIAL, total_energies
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_psitohpsi'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments-------------------------------
!scalars
 integer, intent(in) :: me, nproc, itrp, iter, iscf, natom, nfft, nspden
 real(dp), intent(in) :: alphamix
 real(dp), intent(out) :: rpnrm
 logical, intent(in) :: scf
 logical, intent(in) :: wvlbigdft
 type(wvl_data), intent(inout) :: wvl
 real(dp), intent(inout) :: eexctX,eSIC_DC,ehart,eexcu,vexcu, ekin_sum, epot_sum, eproj_sum
 real(dp), dimension(6), intent(out) :: xcstr
 real(dp), intent(inout) :: xcart(3, natom)
!arrays
 real(dp),intent(out), optional :: vxc(nfft,nspden)
 real(dp),intent(out), optional :: vtrial(nfft,nspden)

!Local variables-------------------------------
!scalars
#if defined HAVE_BIGDFT
 character(len=500) :: message
 integer :: linflag = 0
 character(len=3), parameter :: unblock_comms = "OFF"
#endif

! *************************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT

 if(wvl%descr%atoms%npspcode(1)==7) then
   call psitohpsi(me,nproc,wvl%descr%atoms,scf,wvl%den%denspot, &
&   itrp, iter, iscf, alphamix,&
&   wvl%projectors%nlpsp,xcart,linflag,unblock_comms, &
&   wvl%wfs%GPU,wvl%wfs%ks,wvl%e%energs,rpnrm,xcstr,&
&   wvl%projectors%G,wvl%descr%paw)
 else
   call psitohpsi(me,nproc,wvl%descr%atoms,scf,wvl%den%denspot, &
&   itrp, iter, iscf, alphamix,&
&   wvl%projectors%nlpsp,xcart,linflag,unblock_comms, &
&   wvl%wfs%GPU,wvl%wfs%ks,wvl%e%energs,rpnrm,xcstr)
 end if

 if(scf) then
   ehart     = wvl%e%energs%eh
   eexcu     = wvl%e%energs%exc
   vexcu     = wvl%e%energs%evxc
 end if
 eexctX    = wvl%e%energs%eexctX
 eSIC_DC   = wvl%e%energs%evsic
 ekin_sum  = wvl%e%energs%ekin
 eproj_sum = wvl%e%energs%eproj
 epot_sum  = wvl%e%energs%epot

!Correct local potential, since in BigDFT 
!this variable contains more terms
!Do the following only if sumpion==.true. in psolver_rhohxc.
!For the moment it is set to false.

 epot_sum=epot_sum-real(2,dp)*wvl%e%energs%eh
 epot_sum=epot_sum-wvl%e%energs%evxc

 if(wvlbigdft) then
   call total_energies(wvl%e%energs, iter, me)
 end if

!Note: if evxc is not rested here,
!we have to rest this from etotal in prtene, afterscfcv and etotfor.
!check ABINIT-6.15.1.

 if(scf) then
   if (present(vxc)) then
     write(message, '(a,a,a,a)' ) ch10, ' wvl_psitohpsi : but why are you copying vxc :..o('
     call wrtout(std_out,message,'COLL')
     call wvl_vxc_abi2big(2,vxc,wvl%den)
   end if
   if (wvl%den%denspot%rhov_is == KS_POTENTIAL .and. present(vtrial)) then
     write(message, '(a,a,a,a)' ) ch10, ' wvl_psitohpsi : but why are you copying vtrial :..o('
     call wrtout(std_out,message,'COLL')
     call wvl_vtrial_abi2big(2,vtrial,wvl%den)
   end if
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) me,nproc,itrp,iter,iscf,natom,nfft,nspden,alphamix,rpnrm,scf,&
& wvlbigdft,wvl%wfs%ks,eexctX,eSIC_DC,ehart,eexcu,vexcu,ekin_sum,&
& epot_sum,eproj_sum,xcstr(1),xcart(1,1),vxc(1,1),vtrial(1,1)
#endif

 DBG_EXIT("COLL")

end subroutine wvl_psitohpsi
!!***
