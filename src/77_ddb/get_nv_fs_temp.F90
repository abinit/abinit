!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_nv_fs_temp
!! NAME
!!  get_nv_fs_temp
!!
!! FUNCTION
!! This routine calculates the fermi energy, FD smeared DOS(Ef) and
!! Veloc_sq(Ef) at looped temperatures.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2018 ABINIT group (BXU)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  elph_ds
!!    elph_ds%nband = number of bands in ABINIT
!!    elph_ds%k_fine%nkptirr = Number of irreducible points for which there exist at least one band that crosses the Fermi level.
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%k_fine%nkpt = number of k points for fine k-grid
!!    elph_ds%k_phon%nkpt = number of k points for coarse k-grid
!!    elph_ds%tempermin = minimum temperature at which resistivity etc are calculated (in K)
!!    elph_ds%temperinc = interval temperature grid on which resistivity etc are calculated (in K)
!!    elph_ds%ep_b_min= first band taken into account in FS integration (if telphint==2)
!!    elph_ds%ep_b_max= last band taken into account in FS integration (if telphint==2)
!!    elph_ds%telphint = flag for integration over the FS with 0=tetrahedra 1=gaussians
!!    elph_ds%elphsmear = smearing width for gaussian integration
!!           or buffer in energy for calculations with tetrahedra (telphint=0)
!!
!!  eigenGS = Ground State eigenvalues
!!  gprimd = reciprocal lattice vectors (dimensionful)
!!  kptrlatt_fine = k-point grid vectors (if divided by determinant of present matrix)
!!  max_occ = maximal occupancy for a band
!!
!! OUTPUT
!!  elph_ds%fermie=Fermi level at input temperature
!!  elph_tr_ds%dos_n0=DOS(Ef) at looped temperatures
!!  elph_tr_ds%veloc_sq0=FS averaged velocity at Ef at looped temperatures
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      ebands_update_occ,ep_fs_weights,get_veloc_tr,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine get_nv_fs_temp(elph_ds,BSt,eigenGS,gprimd,max_occ,elph_tr_ds)
    

 use defs_basis
 use defs_datatypes
 use defs_elphon
 use m_profiling_abi
 use m_io_tools

 use m_ebands, only : ebands_update_occ

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_nv_fs_temp'
 use interfaces_14_hidewrite
 use interfaces_77_ddb, except_this_one => get_nv_fs_temp
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!data_type
 type(elph_type),intent(inout) :: elph_ds
 type(ebands_t),intent(inout)   :: BSt
 type(elph_tr_type),intent(inout) :: elph_tr_ds

!Scalars
 real(dp), intent(in) :: max_occ

! arrays
 real(dp), intent(in) :: gprimd(3,3)
 real(dp), intent(in) :: eigenGS(elph_ds%nband,elph_ds%k_fine%nkptirr,elph_ds%nsppol)

!Local variables-------------------------------

 integer :: isppol!, ie1
 integer :: itemp, tmp_nenergy  

 character(len=500) :: message

 real(dp) :: Temp, tmp_elphsmear, tmp_delta_e
! real(dp) :: xtr, e1
! real(dp),allocatable :: tmp_wtk(:,:)
 
! *************************************************************************

 ABI_ALLOCATE(elph_tr_ds%dos_n0,(elph_ds%ntemper,elph_ds%nsppol))
 ABI_ALLOCATE(elph_tr_ds%veloc_sq0,(elph_ds%ntemper,3,elph_ds%nsppol))
!if (elph_ds%use_k_fine == 1) then
!ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt))
!else
!ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_phon%nkpt))
!end if
 
 elph_tr_ds%dos_n0 = zero
 elph_tr_ds%veloc_sq0 = zero
 
 tmp_nenergy = 8
 do itemp=1,elph_ds%ntemper  ! runs over temperature in K
   Temp=elph_ds%tempermin + elph_ds%temperinc*dble(itemp)
   tmp_delta_e = kb_HaK*Temp
   Bst%occopt = 3
   Bst%tsmear = Temp*kb_HaK
   tmp_elphsmear = Temp*kb_HaK
   call ebands_update_occ(Bst,-99.99_dp)
   write(message,'(a,f12.6,a,E20.12)')'At T=',Temp,' Fermi level is:',Bst%fermie
   call wrtout(std_out,message,'COLL')
   if (abs(elph_ds%fermie) < tol10) then
     elph_ds%fermie = BSt%fermie
   end if

!  FD smeared DOS and veloc

   call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, tmp_elphsmear, &
&   elph_ds%fermie, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine,&
&   max_occ, elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&   elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)

   do isppol=1,elph_ds%nsppol
     elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt
     write(message,'(a,f12.6,a,f12.6)')'At T=',Temp,' The DOS at Ef is:', elph_ds%n0(isppol)
     call wrtout(std_out,message,'COLL')

!    For the non-LOVA case, N(Ef) is not that important (canceled out eventually).
!    Should not be important for metal, comment out for now
!    tmp_wtk = zero
!    do ie1=-tmp_nenergy,tmp_nenergy ! use ie1 here, hope there is no confusion
!    e1=Bst%fermie+ie1*tmp_delta_e
!    xtr=(e1-Bst%fermie)/(2.0_dp*kb_HaK*Temp)
!    
!    call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
!    &       e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, &
!    &       max_occ, elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
!    &       elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)
!    
!    tmp_wtk(:,:) = tmp_wtk(:,:) + elph_ds%k_fine%wtk(:,:,isppol)* &
!    &       tmp_delta_e/(4.0d0*kb_HaK*Temp)/(COSH(xtr)**2.0d0)
!    end do ! ie1

!    elph_ds%k_fine%wtk(:,:,isppol) = tmp_wtk(:,:)
     elph_tr_ds%dos_n0(itemp,isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt
!    elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr
!    write(message,'(a,f12.6,a,f12.6)')'At T=',Temp,' The eff. DOS at Ef is:', elph_tr_ds%dos_n0(itemp,isppol)
!    call wrtout(std_out,message,'COLL')
   end do ! isppol
   call get_veloc_tr(elph_ds,elph_tr_ds)
   elph_tr_ds%veloc_sq0(itemp,:,:) = elph_tr_ds%FSelecveloc_sq(:,:)

 end do ! temperature

end subroutine get_nv_fs_temp
!!***
