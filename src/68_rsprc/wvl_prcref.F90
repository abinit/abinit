!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_prcref
!! NAME
!!  wvl_prcref
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2017 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      newrho,newvtr
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_prcref(dielar,iprcel,my_natom,nfftprc,npawmix,nspden,pawrhoij,&
& rhoijrespc,usepaw,vresid,vrespc)
    
 use defs_basis
 use m_profiling_abi
 use defs_wvltypes
 use m_pawrhoij, only : pawrhoij_type
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_prcref'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: iprcel,nfftprc,my_natom,npawmix,nspden,usepaw
 real(dp), intent(in)  :: dielar(7)
 real(dp), intent(in)  :: vresid(nfftprc,nspden)
 real(dp),intent(out) :: rhoijrespc(npawmix)
 real(dp),intent(out) :: vrespc(nfftprc,nspden)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*usepaw)

!Local variables-------------------------------
 integer :: iatom,index,ispden,klmn,kmix
 real(dp):: diemix,diemixmag,mixfac_eff
 character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************
 
 DBG_ENTER("COLL")

!PENDING:
!optres==1, or 0, for density or potential mixing.
!for potential mixing, we have to average over spins.
!check prcref.F90 and moddiel.F90


#if defined HAVE_BIGDFT
#endif

 if(iprcel .ne. 0) then
   write(message, '(a,i3,a,a,a,a)' )&
&   '  From the calling routine, iprcel=',iprcel,ch10,&
&   '  For wavelets, the only allowed value is 0.',ch10,&
&   '  Action : correct your input file.'
   MSG_ERROR(message)
 end if

!call wvl_moddiel  !PENDING
 diemix=dielar(4)
!dielng=dielar(2) ; diemac=dielar(3) ; diemix=dielar(4) ; 
 diemixmag=abs(dielar(7))
 vrespc(:,1)=diemix*vresid(:,1)
 if (nspden/=1) vrespc(:,2:nspden)=diemixmag*vresid(:,2:nspden)

!3) PAW only : precondition the rhoij quantities (augmentation
!occupancies) residuals. Use a simple preconditionning
!with the same mixing factor as the model dielectric function.

 if (usepaw==1.and.my_natom>0) then
!  mixfac=dielar(4);mixfacmag=abs(dielar(7))
   if (pawrhoij(1)%cplex==1) then
     index=0
     do iatom=1,my_natom
       do ispden=1,pawrhoij(iatom)%nspden
         mixfac_eff=diemix;if (ispden>1) mixfac_eff=diemixmag
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           index=index+1;klmn=pawrhoij(iatom)%kpawmix(kmix)
           rhoijrespc(index)=mixfac_eff*pawrhoij(iatom)%rhoijres(klmn,ispden)
         end do
       end do
     end do
   else
     index=-1
     do iatom=1,my_natom
       do ispden=1,pawrhoij(iatom)%nspden
         mixfac_eff=diemix;if (ispden>1) mixfac_eff=diemixmag
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           index=index+2;klmn=2*pawrhoij(iatom)%kpawmix(kmix)-1
           rhoijrespc(index:index+1)=mixfac_eff*pawrhoij(iatom)%rhoijres(klmn:klmn+1,ispden)
         end do
       end do
     end do
   end if
 end if



 DBG_EXIT("COLL")

end subroutine wvl_prcref
!!***
