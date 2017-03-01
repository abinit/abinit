!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_nv_fs_en
!! NAME
!!  get_nv_fs_en
!!
!! FUNCTION
!! This routine finds the energy grids for the integration on epsilon
!! and epsilon prime. It then calculates the DOS and FS averaged velocity_sq at
!! these energies. Metals and semiconductors are treated differently, to deal
!! correctly with the gap.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2017 ABINIT group (BXu)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! crystal<crystal_t>=data type gathering info on the crystalline structure.
!! Ifc<ifc_type>=Object containing the interatomic force constants.
!!  elph_ds
!!    elph_ds%nband = number of bands in ABINIT
!!    elph_ds%k_fine%nkptirr = Number of irreducible points for which there exist at least one band that crosses the Fermi level.
!!    elph_ds%nbranch = number of phonon branches = 3*natom
!!    elph_ds%k_phon%nkpt = number of k points
!!    elph_ds%k_fine%irredtoGS = mapping of elph k-points to ground state grid
!!    elph_ds%minFSband = lowest band included in the FS integration
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%fermie = fermi energy
!!    elph_ds%tempermin = minimum temperature at which resistivity etc are calculated (in K)
!!    elph_ds%temperinc = interval temperature grid on which resistivity etc are calculated (in K)
!!    elph_ds%ep_b_min= first band taken into account in FS integration (if telphint==2)
!!    elph_ds%ep_b_max= last band taken into account in FS integration (if telphint==2)
!!    elph_ds%telphint = flag for integration over the FS with 0=tetrahedra 1=gaussians
!!    elph_ds%elphsmear = smearing width for gaussian integration
!!           or buffer in energy for calculations with tetrahedra (telphint=0)
!!
!!  elph_tr_ds
!!    elph_tr_ds%el_veloc = electronic velocities from the fine k-grid
!!
!!  eigenGS = Ground State eigenvalues
!!  kptrlatt_fine = k-point grid vectors (if divided by determinant of present matrix)
!!  max_occ = maximal occupancy for a band
!!
!! OUTPUT
!!  elph_ds%nenergy = number of energy points for integration on epsilon
!!  elph_tr_ds%en_all = energy points
!!  elph_tr_ds%de_all = differences between energy points
!!  elph_tr_ds%dos_n = DOS at selected energy points
!!  elph_tr_ds%veloc_sq = FS averaged velocity square at selected energy points
!!  elph_tr_ds%tmp_gkk_intweight = integration weights at coarse k grid
!!  elph_tr_ds%tmp_velocwtk = velocity times integration weights at coarse k grid
!!  elph_tr_ds%tmp_vvelocwtk = velocity square times integration weights at coarse k grid
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      d2c_weights,ep_el_weights,ep_fs_weights,get_veloc_tr,ifc_fourq,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine get_nv_fs_en(crystal,ifc,elph_ds,eigenGS,max_occ,elph_tr_ds,omega_max)
    
 use defs_basis
 use defs_elphon
 use m_io_tools
 use m_errors
 use m_profiling_abi
 use m_ifc

 use m_crystal,    only : crystal_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_nv_fs_en'
 use interfaces_14_hidewrite
 use interfaces_77_ddb, except_this_one => get_nv_fs_en
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!Scalars
 real(dp), intent(in)  :: max_occ
 real(dp), intent(out) :: omega_max
 type(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: crystal
 type(elph_type),intent(inout) :: elph_ds
 type(elph_tr_type),intent(inout) :: elph_tr_ds
!Arrays 

 real(dp), intent(in)  :: eigenGS(elph_ds%nband,elph_ds%k_fine%nkptirr,elph_ds%nsppol)

!Local variables-------------------------------
!scalars
 integer ::  iFSqpt,isppol,ie1,ierr
 integer ::  i_metal,low_T
 integer ::  in_nenergy, out_nenergy
 integer ::  n_edge1, n_edge2, edge
 integer ::  ie_all, ne_all
 integer ::  sz1, sz2, sz3, sz4
  real(dp) :: e_vb_max, e_cb_min,ucvol
 real(dp) :: e1,max_e,fine_range
 real(dp) :: enemin,enemax
 real(dp) :: Temp,e_tiny,de0
 real(dp) :: eff_mass1, eff_mass2, tmp_dos
 character(len=500) :: message
!arrays
 real(dp) :: gprimd(3,3)
 real(dp) :: kpt_2nd(3), e_cb_2nd(2), en1(2)
 real(dp),allocatable :: dos_e1(:,:),tmp_wtk(:,:,:,:)
 real(dp),allocatable :: phfrq(:,:)
 real(dp),allocatable :: displ(:,:,:,:)

! *************************************************************************

 gprimd = crystal%gprimd
 ucvol = crystal%ucvol

 Temp             = elph_ds%tempermin+elph_ds%temperinc
 elph_ds%delta_e  = kb_HaK*Temp ! about 1000 cm^-1/100, no need to be omega_max 
 max_e            = elph_ds%nenergy*kb_HaK*Temp
 e_tiny           = kb_HaK*0.00001_dp ! this is the min. delta_e
 de0              = kb_HaK*Temp ! Kb*T

 in_nenergy = elph_ds%nenergy
 
 ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%nsppol,4))
 ABI_ALLOCATE(dos_e1,(elph_ds%nsppol,3))

 ABI_ALLOCATE(phfrq,(elph_ds%nbranch, elph_ds%k_phon%nkpt))
 ABI_ALLOCATE(displ,(2, elph_ds%nbranch, elph_ds%nbranch, elph_ds%k_phon%nkpt))

 do iFSqpt=1,elph_ds%k_phon%nkpt
   call ifc_fourq(ifc,crystal,elph_ds%k_phon%kpt(:,iFSqpt),phfrq(:,iFSqpt),displ(:,:,:,iFSqpt))
 end do

 omega_max = maxval(phfrq)*1.1_dp
 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(displ)

 write(message,'(a,E20.12)')' The max phonon energy is  ', omega_max
 call wrtout(std_out,message,'COLL')

 enemin = elph_ds%fermie - max_e*2
 enemax = elph_ds%fermie + max_e
 call ep_el_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
& enemin, enemax, 4, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
& elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
& elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine, tmp_wtk)

 do isppol=1,elph_ds%nsppol
   dos_e1(isppol,1) = sum(tmp_wtk(:,:,isppol,2))/elph_ds%k_fine%nkpt
   dos_e1(isppol,2) = sum(tmp_wtk(:,:,isppol,3))/elph_ds%k_fine%nkpt
   dos_e1(isppol,3) = sum(tmp_wtk(:,:,isppol,4))/elph_ds%k_fine%nkpt

!  ! BXU, only treat metallic case at this moment, as variational method may not
!  ! apply to insulators
!  i_metal = -1
   i_metal = 1
!  if (dos_e1(isppol,1) .gt. 0.1_dp .and. dos_e1(isppol,2) .gt. 0.1_dp .and. &
!  &   dos_e1(isppol,3) .gt. 0.1_dp) then ! metal
!  i_metal = 1
   if (i_metal == 1) then
     write(message,'(a)')' This is a metal.'
     call wrtout(std_out,message,'COLL')

     fine_range = 1.5_dp
     e1 = elph_ds%fermie + omega_max*fine_range
     out_nenergy = 0
     low_T = 1
     if (omega_max*fine_range .lt. max_e) then
       low_T = 0
       de0 = omega_max*fine_range/in_nenergy ! energy spacing within Ef +/- omega_max 
       do while ((e1-elph_ds%fermie) .lt. max_e)
         e1 = e1 + elph_ds%delta_e
         out_nenergy = out_nenergy + 1
       end do
     end if

     if (low_T == 0) max_e = e1 - elph_ds%fermie
     elph_ds%nenergy = in_nenergy*2 + 1 + out_nenergy*2

   else ! semiconductor/insulator, need careful consideration later
     i_metal = 0
!    between CB min and the next k point, use free electron to replace
!    The weights will be proportional to the DOS, relative to the weights
!    calculated with ep_fs_weights, tetrahedron method prefered

!    output VB and CB edges for semiconductor/insulator
     e_vb_max = maxval(eigenGS(elph_ds%minFSband+elph_ds%nFSband/2-1,:,isppol))
     e_cb_min = minval(eigenGS(elph_ds%minFSband+elph_ds%nFSband/2,:,isppol))
     e_cb_2nd(1) = eigenGS(elph_ds%minFSband+elph_ds%nFSband/2,2,isppol)
     e_cb_2nd(2) = eigenGS(elph_ds%minFSband+elph_ds%nFSband/2+1,2,isppol)
     write(message,'(a,E20.12,2x,E20.12)')' elphon : top of VB, bottom of CB = ',&
&     e_vb_max, e_cb_min
     call wrtout(std_out,message,'COLL')
     write(message,'(a,E20.12)')' elphon : energy at the neighbor kpt = ',e_cb_2nd(1)
     call wrtout(std_out,message,'COLL')

     n_edge1 = 4 ! at the very edge
     n_edge2 = 8  ! sparse to the end of free-electron part

     kpt_2nd(:) = gprimd(:,1)*elph_ds%k_fine%kptirr(1,2) + &
&     gprimd(:,2)*elph_ds%k_fine%kptirr(2,2) + &
&     gprimd(:,3)*elph_ds%k_fine%kptirr(3,2)
     write(message,'(a,3E20.12)')' The neighbor k point is:  ', elph_ds%k_fine%kptirr(:,2)
     call wrtout(std_out,message,'COLL')

     if (dabs(elph_ds%fermie-e_cb_min) .lt. dabs(elph_ds%fermie-e_vb_max)) then
       e1 = e_cb_2nd(1)
     else
       e1 = e_vb_max
     end if
     call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&     e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&     elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&     elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)

     elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt

     eff_mass1 = (kpt_2nd(1)*kpt_2nd(1) + kpt_2nd(2)*kpt_2nd(2) + kpt_2nd(3)*kpt_2nd(3)) / &
&     (2.0_dp*(e_cb_2nd(1)-e_cb_min))
     write(message,'(a,E20.12)')' The eff. mass from band1 is: ', eff_mass1
     call wrtout(std_out,message,'COLL')
     eff_mass2 = (kpt_2nd(1)*kpt_2nd(1) + kpt_2nd(2)*kpt_2nd(2) + kpt_2nd(3)*kpt_2nd(3)) / &
&     (2.0_dp*(e_cb_2nd(2)-e_cb_min))
     write(message,'(a,E20.12)')' The eff. mass from band2 is: ', eff_mass2
     call wrtout(std_out,message,'COLL')

!    bxu, but the eff. mass estimated in this way is too small
!    The following is obtained by roughly fitting to the DOS of 48x48x48
     eff_mass1 = 0.91036
     write(message,'(a,E20.12)')' The eff. mass we are using is: ', eff_mass1
     call wrtout(std_out,message,'COLL')

     tmp_dos = (ucvol/2.0_dp/pi**2.0_dp)*(2.0_dp*eff_mass1)**1.5_dp*(e1-e_cb_min)**0.5_dp + &
&     2.0_dp*(ucvol/2.0_dp/pi**2.0_dp)*(2.0_dp*eff_mass2)**1.5_dp*(e1-e_cb_min)**0.5_dp
     write(message,'(a,E20.12)')' The fake DOS at kpt1 =   ', tmp_dos
     call wrtout(std_out,message,'COLL')
     write(message,'(a,E20.12)')' The calculated DOS at kpt1 =   ', elph_ds%n0(isppol)
     call wrtout(std_out,message,'COLL')


     e1 = elph_ds%fermie - max_e
     ie_all = 1
     ne_all = 0
     edge = 0

     call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&     e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&     elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&     elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)
     
     elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt
     do while ((e1-elph_ds%fermie) .lt. max_e)
       if (e1 .lt. e_cb_min .and. elph_ds%n0(isppol) .lt. tol9) then
         e1 = e_cb_2nd(1)
         edge = 1
         e1 = e1 + de0
       end if

       if (e1 .lt. e_cb_2nd(1)) then
         e1 = e_cb_2nd(1)
         edge = 1
         e1 = e1 + de0
       end if

       if (e1 .gt. e_cb_2nd(1)) then
         call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&         e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&         elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&         elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)
         
         elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt

         e1 = e1 + de0
         ie_all = ie_all + 1
       end if
     end do ! e_all
     ne_all = ie_all - 1 + (n_edge1 + n_edge2 - 1)*edge ! energy levels in the free-electron range
     write(message,'(a,i3,a,i3,a)')' For spin', isppol, '  there are ', &
&     ne_all, '  energy levels considered '
     call wrtout(std_out,message,'COLL')

     elph_ds%nenergy = ne_all
   end if ! metal or insulator
 end do ! isppol

 ABI_DEALLOCATE(tmp_wtk)

 if (elph_ds%nenergy .lt. 2) then 
   MSG_ERROR('There are too few energy levels for non-LOVA')
 end if

 sz1=elph_ds%ngkkband;sz2=elph_ds%k_phon%nkpt
 sz3=elph_ds%nsppol;sz4=elph_ds%nenergy+1
 ABI_ALLOCATE(elph_tr_ds%dos_n,(sz4,sz3))
 ABI_ALLOCATE(elph_tr_ds%veloc_sq,(3,sz3,sz4))
 ABI_ALLOCATE(elph_tr_ds%en_all,(sz3,sz4))
 ABI_ALLOCATE(elph_tr_ds%de_all,(sz3,sz4+1))
 ABI_ALLOCATE(elph_tr_ds%tmp_gkk_intweight,(sz1,sz2,sz3,sz4))
 ABI_ALLOCATE(elph_tr_ds%tmp_velocwtk,(sz1,sz2,3,sz3,sz4))
 ABI_ALLOCATE(elph_tr_ds%tmp_vvelocwtk,(sz1,sz2,3,3,sz3,sz4))

 elph_tr_ds%dos_n = zero
 elph_tr_ds%veloc_sq = zero
 elph_tr_ds%tmp_gkk_intweight = zero
 elph_tr_ds%tmp_velocwtk = zero
 elph_tr_ds%tmp_vvelocwtk = zero

 ABI_STAT_ALLOCATE(elph_ds%k_phon%velocwtk,(elph_ds%nFSband,elph_ds%k_phon%nkpt,3,elph_ds%nsppol), ierr)
 ABI_CHECK(ierr==0, 'allocating elph_ds%k_phon%velocwtk')

 ABI_STAT_ALLOCATE(elph_ds%k_phon%vvelocwtk,(elph_ds%nFSband,elph_ds%k_phon%nkpt,3,3,elph_ds%nsppol), ierr)
 ABI_CHECK(ierr==0, 'allocating elph_ds%k_phon%vvelocwtk')

 elph_ds%k_phon%velocwtk = zero
 elph_ds%k_phon%vvelocwtk = zero

!metal
 if (i_metal .eq. 1) then
   e1 = elph_ds%fermie - max_e
   en1(:) = elph_ds%fermie - max_e
   if (low_T .eq. 1) then
     enemin = elph_ds%fermie - max_e - elph_ds%delta_e
     enemax = elph_ds%fermie + max_e

     ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%nsppol,elph_ds%nenergy+1))
     call ep_el_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&     enemin, enemax, elph_ds%nenergy+1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&     elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&     elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine, tmp_wtk)
     
     do isppol=1,elph_ds%nsppol
       do ie1 = 1, elph_ds%nenergy
         elph_tr_ds%en_all(isppol,ie1) = en1(isppol)
         elph_tr_ds%de_all(isppol,ie1) = elph_ds%delta_e

         elph_ds%k_fine%wtk(:,:,isppol) = tmp_wtk(:,:,isppol,ie1+1)
         elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr
         elph_tr_ds%dos_n(ie1,isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt

         call get_veloc_tr(elph_ds,elph_tr_ds)
         elph_tr_ds%veloc_sq(:,isppol,ie1)=elph_tr_ds%FSelecveloc_sq(:,isppol)

         call d2c_weights(elph_ds,elph_tr_ds)

         elph_tr_ds%tmp_gkk_intweight(:,:,isppol,ie1) = elph_ds%k_phon%wtk(:,:,isppol)
         elph_tr_ds%tmp_velocwtk(:,:,:,isppol,ie1) = elph_ds%k_phon%velocwtk(:,:,:,isppol)
         elph_tr_ds%tmp_vvelocwtk(:,:,:,:,isppol,ie1) = elph_ds%k_phon%vvelocwtk(:,:,:,:,isppol)
         en1(isppol) = en1(isppol) + elph_ds%delta_e
       end do
     end do
     ABI_DEALLOCATE(tmp_wtk)

   else ! low_T = 0
     enemin = e1 - elph_ds%delta_e
     enemax = e1 + (out_nenergy-1)*elph_ds%delta_e

     ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%nsppol,out_nenergy+1))
     call ep_el_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&     enemin, enemax, out_nenergy+1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&     elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&     elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine, tmp_wtk)
     do isppol=1,elph_ds%nsppol
       do ie1 = 1, out_nenergy
         elph_tr_ds%en_all(isppol,ie1) = en1(isppol)
         elph_tr_ds%de_all(isppol,ie1) = elph_ds%delta_e
         
         elph_ds%k_fine%wtk(:,:,isppol) = tmp_wtk(:,:,isppol,ie1+1)
         elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr
         elph_tr_ds%dos_n(ie1,isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt

         call get_veloc_tr(elph_ds,elph_tr_ds)
         elph_tr_ds%veloc_sq(:,isppol,ie1)=elph_tr_ds%FSelecveloc_sq(:,isppol)

         call d2c_weights(elph_ds,elph_tr_ds)

         elph_tr_ds%tmp_gkk_intweight(:,:,isppol,ie1) = elph_ds%k_phon%wtk(:,:,isppol)
         elph_tr_ds%tmp_velocwtk(:,:,:,isppol,ie1) = elph_ds%k_phon%velocwtk(:,:,:,isppol)
         elph_tr_ds%tmp_vvelocwtk(:,:,:,:,isppol,ie1) = elph_ds%k_phon%vvelocwtk(:,:,:,:,isppol)

         en1(isppol) = en1(isppol) + elph_ds%delta_e
       end do
     end do
     ABI_DEALLOCATE(tmp_wtk)

     e1 = en1(1)
     enemin = e1 - de0
     enemax = e1 + in_nenergy*2*de0

     ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%nsppol,in_nenergy*2+2))
     call ep_el_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&     enemin, enemax, in_nenergy*2+2, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&     elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&     elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine, tmp_wtk)

     do isppol=1,elph_ds%nsppol
       do ie1 = out_nenergy+1, out_nenergy+in_nenergy*2+1
         elph_tr_ds%en_all(isppol,ie1) = en1(isppol)
         elph_tr_ds%de_all(isppol,ie1) = de0
         
         elph_ds%k_fine%wtk(:,:,isppol) = tmp_wtk(:,:,isppol,ie1-out_nenergy+1)
         elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr
         elph_tr_ds%dos_n(ie1,isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt

         call get_veloc_tr(elph_ds,elph_tr_ds)
         elph_tr_ds%veloc_sq(:,isppol,ie1)=elph_tr_ds%FSelecveloc_sq(:,isppol)

         call d2c_weights(elph_ds,elph_tr_ds)

         elph_tr_ds%tmp_gkk_intweight(:,:,isppol,ie1) = elph_ds%k_phon%wtk(:,:,isppol)
         elph_tr_ds%tmp_velocwtk(:,:,:,isppol,ie1) = elph_ds%k_phon%velocwtk(:,:,:,isppol)
         elph_tr_ds%tmp_vvelocwtk(:,:,:,:,isppol,ie1) = elph_ds%k_phon%vvelocwtk(:,:,:,:,isppol)

         en1(isppol) = en1(isppol) + de0
       end do
     end do
     ABI_DEALLOCATE(tmp_wtk)

     e1 = en1(1)
     enemin = e1 - elph_ds%delta_e
     enemax = e1 + (out_nenergy-1)*elph_ds%delta_e

     ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%nsppol,out_nenergy+1))
     call ep_el_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&     enemin, enemax, out_nenergy+1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&     elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&     elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine, tmp_wtk)

     en1(:) = en1(:) - de0 + elph_ds%delta_e ! adjust to make the points symmetric around Ef
     do isppol=1,elph_ds%nsppol
       do ie1 = out_nenergy+in_nenergy*2+2, in_nenergy*2+1+out_nenergy*2
         elph_tr_ds%en_all(isppol,ie1) = en1(isppol)
         elph_tr_ds%de_all(isppol,ie1) = elph_ds%delta_e
         
         elph_ds%k_fine%wtk(:,:,isppol) = tmp_wtk(:,:,isppol,ie1-out_nenergy-in_nenergy*2)
         elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr
         elph_tr_ds%dos_n(ie1,isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt

         call get_veloc_tr(elph_ds,elph_tr_ds)
         elph_tr_ds%veloc_sq(:,isppol,ie1)=elph_tr_ds%FSelecveloc_sq(:,isppol)

         call d2c_weights(elph_ds,elph_tr_ds)

         elph_tr_ds%tmp_gkk_intweight(:,:,isppol,ie1) = elph_ds%k_phon%wtk(:,:,isppol)
         elph_tr_ds%tmp_velocwtk(:,:,:,isppol,ie1) = elph_ds%k_phon%velocwtk(:,:,:,isppol)
         elph_tr_ds%tmp_vvelocwtk(:,:,:,:,isppol,ie1) = elph_ds%k_phon%vvelocwtk(:,:,:,:,isppol)

         en1(isppol) = en1(isppol) + elph_ds%delta_e
       end do
     end do
     ABI_DEALLOCATE(tmp_wtk)
   end if

!semiconductor         
 else if (i_metal .eq. 0) then
   e1 = elph_ds%fermie - max_e
   ie_all = 1

   call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&   e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&   elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&   elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)
   
   elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt
   do while ((e1-elph_ds%fermie) .lt. max_e)
     if (e1 .lt. e_cb_min .and. elph_ds%n0(isppol) .lt. tol9) then
       e1 = e_cb_min
     end if

     if (ie_all .ge. n_edge1+n_edge2) then
       if (ie_all .eq. n_edge1+n_edge2) e1 = e1 + de0
       call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&       e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&       elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&       elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)
       
       elph_tr_ds%tmp_gkk_intweight(:,:,isppol,ie_all) = elph_ds%k_fine%wtk(:,:,isppol)
       elph_tr_ds%dos_n(ie_all,isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt
       elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr

       elph_tr_ds%en_all(isppol,ie_all) = e1
       call get_veloc_tr(elph_ds,elph_tr_ds)
       elph_tr_ds%veloc_sq(:,isppol,ie_all)=elph_tr_ds%FSelecveloc_sq(:,isppol)
!      bxu
!      veloc_sq(1,isppol,ie_all) is "1" good and general??

       elph_tr_ds%de_all(isppol,ie_all) = de0
       e1 = e1 + elph_tr_ds%de_all(isppol,ie_all)
       ie_all = ie_all + 1
     else ! divided according to the 1/DOS (evenly)
       if (ie_all .lt. n_edge1) then
         elph_tr_ds%en_all(isppol,ie_all) = e_cb_min + &
&         (e_tiny**(-0.5_dp) - ie_all*(e_tiny**(-0.5_dp)-(e_cb_2nd(1)-e_cb_min)**(-0.5_dp))/ &
&         dble(n_edge1))**(-2.0_dp)
         if (ie_all .gt. 1) then
           elph_tr_ds%de_all(isppol,ie_all) = elph_tr_ds%en_all(isppol,ie_all) - elph_tr_ds%en_all(isppol,ie_all-1)
         else
           elph_tr_ds%de_all(isppol,ie_all) = elph_tr_ds%en_all(isppol,ie_all) - e_cb_min - e_tiny
         end if
         e1 = elph_tr_ds%en_all(isppol,ie_all)
       else
         elph_tr_ds%en_all(isppol,ie_all) = e_cb_min + &
&         ((ie_all-n_edge1+1)/dble(n_edge2))**2.0_dp*(e_cb_2nd(1)-e_cb_min)
         if (ie_all .gt. 1) then
           elph_tr_ds%de_all(isppol,ie_all) = elph_tr_ds%en_all(isppol,ie_all) - elph_tr_ds%en_all(isppol,ie_all-1)
         else
           elph_tr_ds%de_all(isppol,ie_all) = (e_cb_2nd(1)-e_cb_min)/(dble(n_edge2)**2.0_dp)
         end if
         e1 = elph_tr_ds%en_all(isppol,ie_all)
       end if

       call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&       e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&       elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&       elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)
       
       elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr

       tmp_dos = (ucvol/2.0_dp/pi**2.0_dp)*(2.0_dp*eff_mass1)**1.5_dp*(e1-e_cb_min)**0.5_dp + &
&       2.0_dp*(ucvol/2.0_dp/pi**2.0_dp)*(2.0_dp*eff_mass2)**1.5_dp*(e1-e_cb_min)**0.5_dp
       elph_tr_ds%dos_n(ie_all,isppol) = tmp_dos
       elph_tr_ds%tmp_gkk_intweight(:,:,isppol,ie_all) = elph_ds%k_fine%wtk(:,:,isppol)*tmp_dos/elph_ds%n0(isppol)

       call get_veloc_tr(elph_ds,elph_tr_ds)
       elph_tr_ds%veloc_sq(:,isppol,ie_all)=elph_tr_ds%FSelecveloc_sq(:,isppol)

       if (ie_all .eq. (n_edge1+n_edge2)) e1 = e_cb_2nd(1) + de0
       ie_all = ie_all + 1
     end if
   end do ! ie_all
 else
   MSG_BUG('check i_metal!')
 end if ! metal or insulator

 ABI_DEALLOCATE(dos_e1)

end subroutine get_nv_fs_en
!!***
