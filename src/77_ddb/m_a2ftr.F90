!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_a2ftr
!! NAME
!! m_a2ftr
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!   Copyright (C) 2004-2019 ABINIT group (JPC, MJV, BXU)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_a2ftr

 use defs_basis
 use defs_elphon
 use m_errors
 use m_abicore
 use m_xmpi
 use m_splines
 use m_ebands

 use defs_datatypes,    only : ebands_t
 use m_io_tools,        only : open_file
 use m_numeric_tools,   only : simpson_int
 use m_hide_lapack,     only : matrginv
 use m_geometry,        only : phdispl_cart2red
 use m_crystal,         only : crystal_t
 use m_ifc,             only : ifc_type
 use m_dynmat,          only : ftgam_init, ftgam
 use m_epweights,       only : d2c_wtq, ep_ph_weights, ep_el_weights, ep_ph_weights
 use m_fstab,           only : mkqptequiv

 implicit none

 private
!!***

 public :: mka2f_tr
 public :: mka2f_tr_lova
 public :: get_tau_k
!!***

contains
!!***

!!****f* ABINIT/mka2f_tr
!!
!! NAME
!! mka2f_tr
!!
!! FUNCTION
!!  calculates the FS averaged Transport alpha^2F_tr function
!!  calculates and outputs the associated electrical conductivity, relaxation time, and Seebeck coefficient
!!  and thermal conductivities
!!  for the first task : copied from mka2F
!!
!! INPUTS
!! crystal<crystal_t>=data type gathering info on the crystalline structure.
!! Ifc<ifc_type>=Object containing the interatomic force constants.
!!  elph_ds
!!    elph_ds%gkk2 = gkk2 matrix elements on full FS grid for each phonon mode
!!    elph_ds%nbranch = number of phonon branches = 3*natom
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%k_phon%nkpt = number of kpts included in the FS integration
!!    elph_ds%k_phon%wtk = integration weights on the FS
!!    delph_ds%n0 = DOS at the Fermi level calculated from the k_phon integration weights
!!    elph_ds%k_phon%kpt = coordinates of all FS kpoints
!!  mustar = coulomb pseudopotential parameter eventually for 2 spin channels
!!  natom = number of atoms
!!  ntemper = number of temperature points to calculate, from tempermin to tempermin+ntemper*temperinc
!!  tempermin = minimum temperature at which resistivity etc are calculated (in K)
!!  temperinc = interval for temperature grid on which resistivity etc are calculated (in K)
!!  elph_tr_ds%dos_n0 = DOS at the Fermi level calculated from the k_phon integration
!!           weights, but has a temperature dependence
!!  elph_tr_ds%dos_n  = DOS at varied energy levels around Fermi level
!!  elph_tr_ds%veloc_sq0 = Fermi velocity square with T dependence
!!
!! OUTPUT
!!  elph_ds
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      d2c_wtq,dgemm,ep_ph_weights,ftgam,ftgam_init,gam_mult_displ,ifc_fourq
!!      matrginv,simpson_int,wrtout,xmpi_sum,zgemm
!!
!! NOTES
!!   copied from ftiaf9.f
!!
!! SOURCE

subroutine mka2f_tr(crystal,ifc,elph_ds,ntemper,tempermin,temperinc,pair2red,elph_tr_ds)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ntemper
 real(dp),intent(in) :: tempermin,temperinc
 type(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: crystal
 type(elph_tr_type),intent(inout) :: elph_tr_ds
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: pair2red(elph_ds%nenergy,elph_ds%nenergy)

!Local variables -------------------------
!x =w/(2kbT)
!scalars
 integer :: iFSqpt,ibranch,iomega,isppol,jbranch,nerr
 integer :: unit_a2f_tr,natom,ii,jj
 integer :: idir, iatom, k1, kdir
 integer :: unit_lor,unit_rho,unit_tau,unit_sbk, unit_therm
 integer :: itemp, tmp_nenergy
 integer :: itrtensor, icomp, jcomp!, kcomp
 integer :: ie, ie_1, ie2, ie_2, ie1, ie_tmp, ssp, s1(4), s2(4)
 integer :: ie2_left, ie2_right
 integer :: ik_this_proc, ierr,nrpt
 logical,parameter :: debug=.False.
 real(dp) :: Temp,chgu,chtu,chsu,chwu,diagerr,ucvol
 real(dp) :: a2fprefactor, gtemp
 real(dp) :: lambda_tr,lor0,lorentz,maxerr,omega
 real(dp) :: rho,tau,wtherm,xtr
 real(dp) :: lambda_tr_trace
 real(dp) :: domega, omega_min, omega_max
 real(dp) :: gaussval, gaussprefactor, gaussfactor, gaussmaxarg, xx
 real(dp) :: qnorm2, tmp_fermie
 real(dp) :: e1, e2, diff, xe
 real(dp) :: occ_omega, occ_e1, occ_e2
 real(dp) :: nv1, nv2, sigma1, sigma2
 real(dp) :: dos_n_e2, veloc_sq_icomp, veloc_sq_jcomp
 real(dp) :: tointegq00_1, tointegq00_2, tointegq01_1, tointegq01_2,tointegq11_1, tointegq11_2
 real(dp) :: j00, j01, j11
 real(dp) :: tointegq00,tointegq01,tointegq11
 real(dp) :: pref_s, pref_w, tmp_veloc_sq0, tmp_veloc_sq1, tmp_veloc_sq2
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 real(dp),parameter :: c0(2)=(/0.d0,0.d0/),c1(2)=(/1.d0,0.d0/)
 real(dp) ::  gprimd(3,3)
 real(dp) :: eigval(elph_ds%nbranch)
 real(dp) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: gam_now(2,elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: tmpa2f(elph_ds%na2f)
 real(dp) :: tmpgam1(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmpgam2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: q11_inv(3,3)
!real(dp) ::  fullq(6,6)
 real(dp),allocatable :: phfrq(:,:)
 real(dp),allocatable :: tmp_a2f_1d_tr(:,:,:,:,:)
 real(dp),allocatable :: displ(:,:,:,:)
 real(dp),allocatable :: pheigvec(:,:)
 real(dp),allocatable :: tmp_wtq(:,:,:)
 real(dp),allocatable :: integrho(:), tointegrho(:)
 real(dp),allocatable :: integrand_q00(:),integrand_q01(:),integrand_q11(:)
 real(dp),allocatable :: q00(:,:,:,:), q01(:,:,:,:),q11(:,:,:,:)
 real(dp),allocatable :: seebeck(:,:,:,:)!, rho_nm(:,:,:,:)
 real(dp),allocatable :: rho_T(:)
 real(dp),allocatable :: coskr(:,:), sinkr(:,:)

! *********************************************************************
!calculate a2f_tr for frequencies between 0 and omega_max


 write(std_out,*) 'mka2f_tr : enter '

 ucvol = crystal%ucvol
 natom = crystal%natom
 gprimd = crystal%gprimd

 ! number of real-space points for FT interpolation
 nrpt = ifc%nrpt
!
!MG: the step should be calculated locally using nomega and the extrema of the spectrum.
!One should not rely on previous calls for the setup of elph_ds%domega
!I will remove elph_ds%domega since mka2f.F90 will become a method of gamma_t
 domega =elph_ds%domega
 omega_min       = elph_ds%omega_min
 omega_max       = elph_ds%omega_max

 if (elph_ds%ep_lova .eq. 1) then
   tmp_nenergy = 1
 else if (elph_ds%ep_lova .eq. 0) then
   tmp_nenergy = elph_ds%nenergy
 end if

!! defaults for number of temperature steps and max T (all in Kelvin...)
!ntemper=1000
!tempermin=zero
!temperinc=one

 ABI_ALLOCATE(rho_T,(ntemper))


 gaussprefactor = sqrt(piinv) / elph_ds%a2fsmear
 gaussfactor = one / elph_ds%a2fsmear
 gaussmaxarg = sqrt(-log(1.d-90))
!lor0=(pi*kb_HaK)**2/3.
 lor0=pi**2/3.0_dp

!maximum value of frequency (a grid has to be chosen for the representation of alpha^2 F)
!WARNING! supposes this value has been set in mkelph_linwid.

!ENDMG

 maxerr=0.
 nerr=0

 ABI_ALLOCATE(tmp_wtq,(elph_ds%nbranch, elph_ds%k_fine%nkpt, elph_ds%na2f+1))
 ABI_ALLOCATE(elph_ds%k_fine%wtq,(elph_ds%nbranch, elph_ds%k_fine%nkpt, elph_ds%na2f))
 ABI_ALLOCATE(elph_ds%k_phon%wtq,(elph_ds%nbranch, elph_ds%k_phon%nkpt, elph_ds%na2f))

 ABI_ALLOCATE(phfrq,(elph_ds%nbranch, elph_ds%k_fine%nkpt))
 ABI_ALLOCATE(displ,(2, elph_ds%nbranch, elph_ds%nbranch, elph_ds%k_fine%nkpt))
 ABI_ALLOCATE(pheigvec,(2*elph_ds%nbranch*elph_ds%nbranch, elph_ds%k_fine%nkpt))

 do iFSqpt=1,elph_ds%k_fine%nkpt
   call ifc%fourq(crystal,elph_ds%k_fine%kpt(:,iFSqpt),phfrq(:,iFSqpt),displ(:,:,:,iFSqpt),out_eigvec=pheigvec(:,iFSqpt))
 end do
 omega_min = omega_min - domega

!bxu, obtain wtq for the q_fine, then condense to q_phon
 call ep_ph_weights(phfrq,elph_ds%a2fsmear,omega_min,omega_max,elph_ds%na2f+1,gprimd,elph_ds%kptrlatt_fine, &
& elph_ds%nbranch,elph_ds%telphint,elph_ds%k_fine,tmp_wtq)
!call ep_ph_weights(phfrq,elph_ds%a2fsmear,omega_min,omega_max,elph_ds%na2f+1,gprimd,elph_ds%kptrlatt_fine, &
!& elph_ds%nbranch,1,elph_ds%k_fine,tmp_wtq)
 omega_min = omega_min + domega

 do iomega = 1, elph_ds%na2f
   elph_ds%k_fine%wtq(:,:,iomega) = tmp_wtq(:,:,iomega+1)
 end do
 ABI_DEALLOCATE(tmp_wtq)

 if (elph_ds%use_k_fine == 1) then
   call d2c_wtq(elph_ds)
 end if

 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(displ)
 ABI_DEALLOCATE(pheigvec)

!reduce the dimension from fine to phon for phfrq and pheigvec
 ABI_ALLOCATE(phfrq,(elph_ds%nbranch, elph_ds%k_phon%nkpt))
 ABI_ALLOCATE(displ,(2, elph_ds%nbranch, elph_ds%nbranch, elph_ds%k_phon%nkpt))
 ABI_ALLOCATE(pheigvec,(2*elph_ds%nbranch*elph_ds%nbranch, elph_ds%k_phon%nkpt))

 do iFSqpt=1,elph_ds%k_phon%nkpt
   call ifc%fourq(crystal,elph_ds%k_phon%kpt(:,iFSqpt),phfrq(:,iFSqpt),displ(:,:,:,iFSqpt),out_eigvec=pheigvec(:,iFSqpt))
 end do

 ABI_ALLOCATE(elph_tr_ds%a2f_1d_tr,(elph_ds%na2f,9,elph_ds%nsppol,4,tmp_nenergy**2,ntemper))
 ABI_ALLOCATE(tmp_a2f_1d_tr,(elph_ds%na2f,9,elph_ds%nsppol,4,tmp_nenergy**2))

! prepare phase factors
 ABI_ALLOCATE(coskr, (elph_ds%k_phon%nkpt, nrpt))
 ABI_ALLOCATE(sinkr, (elph_ds%k_phon%nkpt, nrpt))
 call ftgam_init(ifc%gprim, elph_ds%k_phon%nkpt, nrpt, elph_ds%k_phon%kpt, ifc%rpt, coskr, sinkr)

 elph_tr_ds%a2f_1d_tr = zero
 tmp_a2f_1d_tr = zero


 do ie = 1, elph_ds%n_pair
   do ssp = 1,4
     do isppol = 1, elph_ds%nsppol

!      loop over qpoint in full kpt grid (presumably dense)
       do ik_this_proc =1,elph_ds%k_phon%my_nkpt
         iFSqpt = elph_ds%k_phon%my_ikpt(ik_this_proc)

         qnorm2 = sum(elph_ds%k_phon%kpt(:,iFSqpt)**2)
!        if (flag_to_exclude_soft_modes = .false.) qnorm2 = zero
         do itrtensor=1,9

           if (elph_ds%ep_int_gkk == 1) then
             gam_now(:,:) = elph_tr_ds%gamma_qpt_tr(:,itrtensor,:,isppol,iFSqpt)
           else
!            Do FT from real-space gamma grid to 1 qpt.
             call ftgam(ifc%wghatm,gam_now,elph_tr_ds%gamma_rpt_tr(:,itrtensor,:,isppol,:,ssp,ie),natom,1,nrpt,0,&
&             coskr(iFSqpt,:), sinkr(iFSqpt,:))
           end if

!          Diagonalize gamma matrix at this qpoint (complex matrix).

!          if ep_scalprod==0 we have to dot in the displacement vectors here
           if (elph_ds%ep_scalprod==0) then

             displ_red(:,:,:) = zero
             do jbranch=1,elph_ds%nbranch
               do iatom=1,natom
                 do idir=1,3
                   ibranch=idir+3*(iatom-1)
                   do kdir=1,3
                     k1 = kdir+3*(iatom-1)
                     displ_red(1,ibranch,jbranch) = displ_red(1,ibranch,jbranch) + &
&                     gprimd(kdir,idir)*displ(1,k1,jbranch,iFSqpt)
                     displ_red(2,ibranch,jbranch) = displ_red(2,ibranch,jbranch) + &
&                     gprimd(kdir,idir)*displ(2,k1,jbranch,iFSqpt)
                   end do
                 end do
               end do
             end do

             tmpgam2 = reshape (gam_now, (/2,elph_ds%nbranch,elph_ds%nbranch/))
             call gam_mult_displ(elph_ds%nbranch, displ_red, tmpgam2, tmpgam1)
             do jbranch=1,elph_ds%nbranch
               eigval(jbranch) = tmpgam1(1, jbranch, jbranch)
             end do

           else if (elph_ds%ep_scalprod == 1) then


!            NOTE: in these calls gam_now and pheigvec do not have the right rank, but blas usually does not care

             call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, c1, gam_now, 3*natom,&
&             pheigvec(:,iFSqpt), 3*natom, c0, tmpgam1, 3*natom)
             call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, c1, pheigvec(:,iFSqpt), 3*natom,&
&             tmpgam1, 3*natom, c0, tmpgam2, 3*natom)
             diagerr = zero
             do ibranch=1,elph_ds%nbranch
               eigval(ibranch) = tmpgam2(1,ibranch,ibranch)
               do jbranch=1,ibranch-1
                 diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
               end do
               do jbranch=ibranch+1,elph_ds%nbranch
                 diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
               end do
             end do
             if (diagerr > tol12) then
               nerr=nerr+1
               maxerr=max(diagerr, maxerr)
             end if
           end if  ! end ep_scalprod if

!          Add all contributions from the phonon modes at this qpoint to a2f and the phonon dos.
           do ibranch=1,elph_ds%nbranch
!            if (abs(phfrq(ibranch,iFSqpt)) < tol10) then
             if ( abs(phfrq(ibranch,iFSqpt)) < tol7 .or. &
&             (phfrq(ibranch,iFSqpt) < tol4 .and. qnorm2 > 0.03 )) then
!              note: this should depend on the velocity of sound, to accept acoustic modes!
               a2fprefactor = zero
             else
!              a2fprefactor  = eigval (ibranch)/(two_pi*abs(phfrq(ibranch,iFSqpt))*elph_ds%n0(isppol))
!              Use the dos_n0 at the lowest input temperature, assuming to be low
               a2fprefactor  = eigval (ibranch)/(two_pi*abs(phfrq(ibranch,iFSqpt)))
             end if

             omega = omega_min
             tmpa2f(:) = zero
             do iomega=1,elph_ds%na2f
               xx = (omega-phfrq(ibranch,iFSqpt))*gaussfactor
               omega = omega+domega
               if (abs(xx) > gaussmaxarg) cycle

               gaussval = gaussprefactor*exp(-xx*xx)
               gtemp = gaussval*a2fprefactor

               if (dabs(gtemp) < 1.0d-50) gtemp = zero
               tmpa2f(iomega) = tmpa2f(iomega) + gtemp
             end do

!            tmpa2f(:) = zero
!            do iomega=1,elph_ds%na2f
!            gtemp = a2fprefactor*elph_ds%k_phon%wtq(ibranch,iFSqpt,iomega)
!            if (dabs(gtemp) < 1.0d-50) gtemp = zero
!            tmpa2f(iomega) = tmpa2f(iomega) + gtemp
!            end do

             tmp_a2f_1d_tr (:,itrtensor,isppol,ssp,ie) = tmp_a2f_1d_tr (:,itrtensor,isppol,ssp,ie) + tmpa2f(:)

           end do ! end ibranch
         end do ! end itrtensor
       end do ! end iFSqpt  - loop done in parallel
     end do ! end isppol
   end do ! ss'
 end do ! n_pair

 ! MG: FIXME: Why xmpi_world? besides only one CPU should perform IO (see below)
 ! Likely this routine is never executed in parallel
 call xmpi_sum (tmp_a2f_1d_tr, xmpi_world, ierr)

 ABI_DEALLOCATE(coskr)
 ABI_DEALLOCATE(sinkr)

 do itemp=1,ntemper  ! runs over termperature in K
   do isppol=1,elph_ds%nsppol
     do jj=1,tmp_nenergy**2
       do ii=1,4
         elph_tr_ds%a2f_1d_tr(:,:,isppol,ii,jj,itemp) = tmp_a2f_1d_tr(:,:,isppol,ii,jj)/elph_tr_ds%dos_n0(itemp,isppol)
       end do
     end do
   end do
 end do

 ABI_DEALLOCATE(tmp_a2f_1d_tr)

!second 1 / elph_ds%k_phon%nkpt factor for the integration weights
 elph_tr_ds%a2f_1d_tr  = elph_tr_ds%a2f_1d_tr  / elph_ds%k_phon%nkpt

 if (elph_ds%ep_scalprod == 1) then
   write(std_out,*) 'mka2f_tr: errors in diagonalization of gamma_tr with phon eigenvectors: ', nerr,maxerr
 end if

!output the elph_tr_ds%a2f_1d_tr
 fname = trim(elph_ds%elph_base_name) // '_A2F_TR'
 if (open_file(fname,message,newunit=unit_a2f_tr,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

 write (unit_a2f_tr,'(a)')       '#'
 write (unit_a2f_tr,'(a)')       '# ABINIT package : a2f_tr file'
 write (unit_a2f_tr,'(a)')       '#'
 write (unit_a2f_tr,'(a)')       '# a2f_tr function integrated over the FS. omega in a.u.'
 write (unit_a2f_tr,'(a,I10)')   '#     number of kpoints integrated over : ', elph_ds%k_phon%nkpt
 write (unit_a2f_tr,'(a,I10)')   '#     number of energy points : ',elph_ds%na2f
 write (unit_a2f_tr,'(a,E16.6,a,E16.6,a)') '#       between omega_min = ', omega_min,' Ha and omega_max = ', omega_max, ' Ha'
 write (unit_a2f_tr,'(a,E16.6)') '#   and the smearing width for gaussians is ', elph_ds%a2fsmear
 write (unit_a2f_tr,'(a)')       '#'

!done with header
 do isppol=1,elph_ds%nsppol
   write (unit_a2f_tr,'(a,E16.6)') '# The DOS at Fermi level is ', elph_tr_ds%dos_n0(1,isppol)
   omega = omega_min
   do iomega=1,elph_ds%na2f
!    bxu, at which eps and eps' should I save it
!    better to save them all, but could be too many
     write (unit_a2f_tr,   '(10D16.6)') omega, elph_tr_ds%a2f_1d_tr(iomega,:,isppol,1,INT(elph_ds%n_pair/2)+1,1)
     omega=omega+domega
   end do
   write (unit_a2f_tr,*)
 end do !isppol

 close (unit=unit_a2f_tr)

!calculation of transport properties
 ABI_ALLOCATE(integrho,(elph_ds%na2f))
 ABI_ALLOCATE(tointegrho,(elph_ds%na2f))
 ABI_ALLOCATE(integrand_q00,(elph_ds%na2f))
 ABI_ALLOCATE(integrand_q01,(elph_ds%na2f))
 ABI_ALLOCATE(integrand_q11,(elph_ds%na2f))
 ABI_ALLOCATE(q00,(ntemper,3,3,elph_ds%nsppol))
 ABI_ALLOCATE(q01,(ntemper,3,3,elph_ds%nsppol))
 ABI_ALLOCATE(q11,(ntemper,3,3,elph_ds%nsppol))
 ABI_ALLOCATE(seebeck,(elph_ds%nsppol,ntemper,3,3))
!ABI_ALLOCATE(rho_nm,(elph_ds%nsppol,ntemper,3,3))

 fname = trim(elph_ds%elph_base_name) // '_RHO'
 if (open_file(fname,message,newunit=unit_rho,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if
!print header to resistivity file
 write (unit_rho,*) '# Resistivity as a function of temperature.'
 write (unit_rho,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_rho,*) '#  '
 write (unit_rho,*) '#  Columns are: '
 write (unit_rho,*) '#  temperature[K]   rho[au]   rho [SI]        rho/temp [au]'
 write (unit_rho,*) '#  '

 fname = trim(elph_ds%elph_base_name) // '_TAU'
 if (open_file(fname,message,newunit=unit_tau,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if
!print header to relaxation time file
 write (unit_tau,*) '# Relaxation time as a function of temperature.'
 write (unit_tau,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_tau,*) '#  '
 write (unit_tau,*) '#  Columns are: '
 write (unit_tau,*) '#  temperature[K]   tau[au]   tau [SI]     '
 write (unit_tau,*) '#  '

 fname = trim(elph_ds%elph_base_name) // '_SBK'
 if (open_file(fname,message,newunit=unit_sbk,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if
!print header to relaxation time file
 write (unit_sbk,*) '# Seebeck Coefficint as a function of temperature.'
 write (unit_sbk,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_sbk,*) '#  '
 write (unit_sbk,*) '#  Columns are: '
 write (unit_sbk,*) '#  temperature[K]   S [au]   S [SI]     '
 write (unit_sbk,*) '#  '

 fname = trim(elph_ds%elph_base_name) // '_WTH'
 if (open_file(fname,message,newunit=unit_therm,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to thermal conductivity file
 write (unit_therm,'(a)') '# Thermal conductivity/resistivity as a function of temperature.'
 write (unit_therm,'(a)') '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_therm,'(a)') '#  '
 write (unit_therm,'(a)') '#  Columns are: '
 write (unit_therm,'(a)') '#  temperature[K]   thermal rho[au]   thermal cond [au]   thermal rho [SI]   thermal cond [SI]'
 write (unit_therm,'(a)') '#  '

 fname = trim(elph_ds%elph_base_name) // '_LOR'
 if (open_file(fname,message,newunit=unit_lor,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to lorentz file
 write (unit_lor,*) '# Lorentz number as a function of temperature.'
 write (unit_lor,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_lor,*) '#  '
 write (unit_lor,*) '#  Columns are: '
 write (unit_lor,*) '#  temperature[K]   Lorentz number[au]   Lorentz quantum = (pi*kb_HaK)**2/3'
 write (unit_lor,*) '#  '

 do isppol=1,elph_ds%nsppol
   lambda_tr_trace = zero
   do itrtensor=1,9
     omega = omega_min
     tointegrho = zero
     do iomega=1,elph_ds%na2f
       if(omega<=0) then
         omega=omega+domega
         cycle
       end if
!      bxu, agian, which eps and eps' to use?
!      sometimes Ef is in the gap
       tointegrho(iomega)=two*elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,1,1,1)/omega
       omega=omega+domega
     end do

     integrho = zero
     call simpson_int(elph_ds%na2f,domega,tointegrho,integrho)
     lambda_tr=integrho(elph_ds%na2f)
     write (message, '(a,2i3,a,es16.6)' )&
&     ' mka2f_tr: TRANSPORT lambda for isppol itrtensor', isppol, itrtensor, ' =  ', lambda_tr
     call wrtout(std_out,message,'COLL')
     if (itrtensor == 1 .or. itrtensor == 5 .or. itrtensor == 9) lambda_tr_trace = lambda_tr_trace + lambda_tr
   end do !end itrtensor do

   lambda_tr_trace = lambda_tr_trace / three
   write (message, '(a,i3,a,es16.6)' )&
&   ' mka2f_tr: 1/3 trace of TRANSPORT lambda for isppol ', isppol, ' =  ', lambda_tr_trace
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 end do !end isppol do

!constant to change units of rho from au to SI
 chgu=2.173969*1.0d-7
 chtu=2.4188843265*1.0d-17
 chsu=8.617343101*1.0d-5 ! au to J/C/K
 chwu=9.270955772*1.0d-5 ! au to mK/W

!change the fermi level to zero, as required for q01 to vanish.
 tmp_fermie = elph_ds%fermie
!Get Q00, Q01, Q11, and derive rho, tau
 q00 = zero
 q01 = zero
 q11 = zero
! prepare s1 and s2 arrays
 s1 = (/1, 1, -1, -1/)
 s2 = (/1, -1, 1, -1/)

 do isppol=1,elph_ds%nsppol
   do icomp=1, 3
     do jcomp=1, 3
       itrtensor=(icomp-1)*3+jcomp

       write(unit_rho,*) '# Rho for isppol, itrten = ', isppol, itrtensor
       write(unit_tau,*) '# Tau for isppol, itrten = ', isppol, itrtensor

       do itemp=1,ntemper  ! runs over termperature in K
         Temp=tempermin+temperinc*dble(itemp)
         tmp_veloc_sq0 = sqrt(elph_tr_ds%veloc_sq0(itemp,icomp,isppol)*elph_tr_ds%veloc_sq0(itemp,jcomp,isppol))

         integrand_q00 = zero
         integrand_q01 = zero
         integrand_q11 = zero

         omega = omega_min
         do iomega=1,elph_ds%na2f
           if(omega .le. 0) then
             omega=omega+domega
             cycle
           end if
           xtr=omega/(kb_HaK*Temp)
           occ_omega=1.0_dp/(exp(xtr)-1.0_dp)

           tmp_veloc_sq1 = zero
           tmp_veloc_sq2 = zero
           do ie1=1,elph_ds%nenergy
             e1 = elph_tr_ds%en_all(isppol,ie1)

!! BXU, the tolerance here needs to be used with caution
!! which depends on the dimensions of the system
!! E.g. in 2D system, DOS can be much smaller
             if (elph_tr_ds%dos_n(ie1,isppol)/natom .lt. 0.1d0) cycle ! energy in the gap forbidden

             xtr=(e1-tmp_fermie)/(kb_HaK*Temp)
             occ_e1=1.0_dp/(exp(xtr)+1.0_dp)

             e2 = e1 + omega
             xtr=(e2-tmp_fermie)/(kb_HaK*Temp)
             occ_e2=1.0_dp/(exp(xtr)+1.0_dp)
!            Do we need to change the fermie to the one with T dependence?
!            find which ie2 give the closest energy
             if (e2 .gt. elph_tr_ds%en_all(isppol,elph_ds%nenergy)) then
               ie2 = 0
               cycle
             else
               ie_tmp = 1
               diff = dabs(e2-elph_tr_ds%en_all(isppol,1))
               do ie2 = 2, elph_ds%nenergy
                 if (dabs(e2-elph_tr_ds%en_all(isppol,ie2)) .lt. diff) then
                   diff = dabs(e2-elph_tr_ds%en_all(isppol,ie2))
                   ie_tmp = ie2
                 end if
               end do
               ie2 = ie_tmp

               if (e2 < elph_tr_ds%en_all(isppol,ie2)) then
                 ie2_right = ie2
                 ie2_left  = ie2-1
               else
                 ie2_right = ie2+1
                 ie2_left  = ie2
               end if

             end if

             if (elph_tr_ds%dos_n(ie2,isppol)/natom .lt. 0.1d0) cycle

             tointegq00 = zero
             tointegq01 = zero
             tointegq11 = zero

! BXU linear interpolation of volec_sq and dos_n
             xe=(e2-elph_tr_ds%en_all(isppol,ie2_left))/ &
&             (elph_tr_ds%en_all(isppol,ie2_right)-elph_tr_ds%en_all(isppol,ie2_left))
             veloc_sq_icomp = elph_tr_ds%veloc_sq(icomp,isppol,ie2_left)*(1.0d0-xe) + &
&             elph_tr_ds%veloc_sq(icomp,isppol,ie2_right)*xe
             veloc_sq_jcomp = elph_tr_ds%veloc_sq(jcomp,isppol,ie2_left)*(1.0d0-xe) + &
&             elph_tr_ds%veloc_sq(jcomp,isppol,ie2_right)*xe
             dos_n_e2 = elph_tr_ds%dos_n(ie2_left,isppol)*(1.0d0-xe) + &
&             elph_tr_ds%dos_n(ie2_right,isppol)*xe

             tmp_veloc_sq1 = sqrt(elph_tr_ds%veloc_sq(icomp,isppol,ie1)*elph_tr_ds%veloc_sq(jcomp,isppol,ie1))
!             tmp_veloc_sq2 = sqrt(elph_tr_ds%veloc_sq(icomp,isppol,ie2)*elph_tr_ds%veloc_sq(jcomp,isppol,ie2))
             tmp_veloc_sq2 = sqrt(veloc_sq_icomp*veloc_sq_jcomp)

!            ie_1 = (ie1-1)*elph_ds%nenergy + ie2
!            ie_2 = (ie2-1)*elph_ds%nenergy + ie1
             ie_1 = pair2red(ie1,ie2)
             ie_2 = pair2red(ie2,ie1)
             if (ie_1 .eq. 0 .or. ie_2 .eq. 0) then
               MSG_BUG('CHECK pair2red!')
             end if

             do ssp=1,4  ! (s,s'=+/-1, condense the indices)

               nv1 = 1.0_dp/(elph_tr_ds%dos_n(ie1,isppol)*sqrt(tmp_veloc_sq1))
               sigma1 = sqrt(3.0_dp)*(e1-tmp_fermie)/(pi*Temp*kb_HaK)
!DEBUG
               if (elph_ds%ep_lova .eq. 1) then
                 nv1 = 1.0_dp/(elph_tr_ds%dos_n0(itemp,isppol)*sqrt(tmp_veloc_sq0))
                 sigma1 = sqrt(3.0_dp)*(e1-tmp_fermie)/(pi*Temp*kb_HaK)
               end if
!ENDDEBUG

               tointegq00_1 = zero
               tointegq01_1 = zero
               tointegq11_1 = zero

!DEBUG
               if (elph_ds%ep_lova .eq. 1) then
                 nv2 = 1.0_dp/(elph_tr_ds%dos_n0(itemp,isppol)*sqrt(tmp_veloc_sq0))
                 sigma2 = sqrt(3.0_dp)*(e2-tmp_fermie)/(pi*Temp*kb_HaK)
                 j00 = (nv1 + s1(ssp)*nv2)*(nv1 + s2(ssp)*nv2)/4.0_dp
                 j01 = (nv1 + s1(ssp)*nv2)*(nv1*sigma1 + s2(ssp)*nv2*sigma2)/4.0_dp
                 j11 = (nv1*sigma1 + s1(ssp)*nv2*sigma2)*(nv1*sigma1 + s2(ssp)*nv2*sigma2)/4.0_dp
                 tointegq00_1 = elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,1,itemp)* &
&                 occ_e1*(1.0_dp-occ_e2)*j00*occ_omega
                 tointegq01_1 = elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,1,itemp)* &
&                 occ_e1*(1.0_dp-occ_e2)*j01*occ_omega
                 tointegq11_1 = elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,1,itemp)* &
&                 occ_e1*(1.0_dp-occ_e2)*j11*occ_omega
!END DEBUG
               else if (elph_ds%ep_lova .eq. 0) then
                 nv2 = 1.0_dp/(dos_n_e2*sqrt(tmp_veloc_sq2))
                 sigma2 = sqrt(3.0_dp)*(e2-tmp_fermie)/(pi*Temp*kb_HaK)
                 j00 = (nv1 + s1(ssp)*nv2)*(nv1 + s2(ssp)*nv2)/4.0_dp
                 j01 = (nv1 + s1(ssp)*nv2)*(nv1*sigma1 + s2(ssp)*nv2*sigma2)/4.0_dp
                 j11 = (nv1*sigma1 + s1(ssp)*nv2*sigma2)*(nv1*sigma1 + s2(ssp)*nv2*sigma2)/4.0_dp
!                bxu TEST
                 if (debug) then
                   if (ssp .eq. 1 .and. itrtensor .eq. 1) then
                     write(21,"(3i5,4E20.12)") iomega, ie1, ie2, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp), j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp)*j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp)*j01*occ_e1*(1.0_dp-occ_e2)*occ_omega
                   end if
                   if (ssp .eq. 2 .and. itrtensor .eq. 1) then
                     write(22,"(3i5,4E20.12)") iomega, ie1, ie2, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp), j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp)*j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp)*j01*occ_e1*(1.0_dp-occ_e2)*occ_omega
                   end if
                   if (ssp .eq. 3 .and. itrtensor .eq. 1) then
                     write(23,"(3i5,4E20.12)") iomega, ie1, ie2, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp), j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp)*j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp)*j01*occ_e1*(1.0_dp-occ_e2)*occ_omega
                   end if
                   if (ssp .eq. 4 .and. itrtensor .eq. 1) then
                     write(24,"(3i5,4E20.12)") iomega, ie1, ie2, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp), j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp)*j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp)*j01*occ_e1*(1.0_dp-occ_e2)*occ_omega
                   end if
                 end if
                 tointegq00_1 = elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp)* &
&                 occ_e1*(1.0_dp-occ_e2)*j00*occ_omega
                 tointegq01_1 = elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp)* &
&                 occ_e1*(1.0_dp-occ_e2)*j01*occ_omega
                 tointegq11_1 = elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_1,itemp)* &
&                 occ_e1*(1.0_dp-occ_e2)*j11*occ_omega
               end if

               tointegq00_2 = zero
               tointegq01_2 = zero
               tointegq11_2 = zero

!DEBUG
               if (elph_ds%ep_lova .eq. 1) then
                 nv2 = 1.0_dp/(elph_tr_ds%dos_n0(itemp,isppol)*sqrt(tmp_veloc_sq0))
                 sigma2 = sqrt(3.0_dp)*(e2-tmp_fermie)/(pi*Temp*kb_HaK)
                 j00 = (nv2 + s1(ssp)*nv1)*(nv2 + s2(ssp)*nv1)/4.0_dp
                 j01 = (nv2 + s1(ssp)*nv1)*(nv2*sigma2 + s2(ssp)*nv1*sigma1)/4.0_dp
                 j11 = (nv2*sigma2 + s1(ssp)*nv1*sigma1)*(nv2*sigma2 + s2(ssp)*nv1*sigma1)/4.0_dp
                 tointegq00_2 = elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,1,itemp)* &
&                 occ_e1*(1.0_dp-occ_e2)*j00*(occ_omega+1)
                 tointegq01_2 = elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,1,itemp)* &
&                 occ_e1*(1.0_dp-occ_e2)*j01*(occ_omega+1)
                 tointegq11_2 = elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,1,itemp)* &
&                 occ_e1*(1.0_dp-occ_e2)*j11*(occ_omega+1)
!END DEBUG
               else if (elph_ds%ep_lova .eq. 0) then
                 nv2 = 1.0_dp/(dos_n_e2*sqrt(tmp_veloc_sq2))
                 sigma2 = sqrt(3.0_dp)*(e2-tmp_fermie)/(pi*Temp*kb_HaK)
                 j00 = (nv2 + s1(ssp)*nv1)*(nv2 + s2(ssp)*nv1)/4.0_dp
                 j01 = (nv2 + s1(ssp)*nv1)*(nv2*sigma2 + s2(ssp)*nv1*sigma1)/4.0_dp
                 j11 = (nv2*sigma2 + s1(ssp)*nv1*sigma1)*(nv2*sigma2 + s2(ssp)*nv1*sigma1)/4.0_dp
!DEBUG           bxu TEST
                 if (debug) then
                   if (ssp .eq. 1 .and. itrtensor .eq. 1) then
                     write(31,"(3i5,4E20.12)") iomega, ie2, ie1, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp), j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp)*j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp)*j01*occ_e2*(1.0_dp-occ_e1)*(occ_omega+1)
                   end if
                   if (ssp .eq. 2 .and. itrtensor .eq. 1) then
                     write(32,"(3i5,4E20.12)") iomega, ie2, ie1, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp), j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp)*j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp)*j01*occ_e2*(1.0_dp-occ_e1)*(occ_omega+1)
                   end if
                   if (ssp .eq. 3 .and. itrtensor .eq. 1) then
                     write(33,"(3i5,4E20.12)") iomega, ie2, ie1, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp), j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp)*j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp)*j01*occ_e2*(1.0_dp-occ_e1)*(occ_omega+1)
                   end if
                   if (ssp .eq. 4 .and. itrtensor .eq. 1) then
                     write(34,"(3i5,4E20.12)") iomega, ie2, ie1, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp), j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp)*j01, &
&                     elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp)*j01*occ_e2*(1.0_dp-occ_e1)*(occ_omega+1)
                   end if
                 end if
!ENDDEBUG
                 tointegq00_2 = elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp)* &
&                 occ_e2*(1.0_dp-occ_e1)*j00*(occ_omega+1)
                 tointegq01_2 = elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp)* &
&                 occ_e2*(1.0_dp-occ_e1)*j01*(occ_omega+1)
                 tointegq11_2 = elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,ssp,ie_2,itemp)* &
&                 occ_e2*(1.0_dp-occ_e1)*j11*(occ_omega+1)
               end if ! elph_ds%ep_lova

               tointegq00 = tointegq00 + tointegq00_1 + tointegq00_2
               tointegq01 = tointegq01 + tointegq01_1 + tointegq01_2
               tointegq11 = tointegq11 + tointegq11_1 + tointegq11_2

             end do ! ss' = 4
             integrand_q00(iomega) = integrand_q00(iomega) + elph_tr_ds%de_all(isppol,ie1)*tointegq00
             integrand_q01(iomega) = integrand_q01(iomega) + elph_tr_ds%de_all(isppol,ie1)*tointegq01
             integrand_q11(iomega) = integrand_q11(iomega) + elph_tr_ds%de_all(isppol,ie1)*tointegq11
           end do ! ie1 ~ 20
           omega=omega+domega
           q00(itemp,icomp,jcomp,isppol) = q00(itemp,icomp,jcomp,isppol) +&
&           domega*integrand_q00(iomega)
           q01(itemp,icomp,jcomp,isppol) = q01(itemp,icomp,jcomp,isppol) +&
&           domega*integrand_q01(iomega)
           q11(itemp,icomp,jcomp,isppol) = q11(itemp,icomp,jcomp,isppol) +&
&           domega*integrand_q11(iomega)
         end do ! omega ~ 400

         q00(itemp,icomp,jcomp,isppol)=q00(itemp,icomp,jcomp,isppol)* &
&         ucvol*two_pi*elph_tr_ds%dos_n0(itemp,isppol)/(kb_HaK*Temp)
         q01(itemp,icomp,jcomp,isppol)=q01(itemp,icomp,jcomp,isppol)* &
&         ucvol*two_pi*elph_tr_ds%dos_n0(itemp,isppol)/(kb_HaK*Temp)
         q11(itemp,icomp,jcomp,isppol)=q11(itemp,icomp,jcomp,isppol)* &
&         ucvol*two_pi*elph_tr_ds%dos_n0(itemp,isppol)/(kb_HaK*Temp)

         rho = 0.5_dp*q00(itemp,icomp,jcomp,isppol)
!        is tau energy dependent?
         tau = 2.0d0*ucvol/(q00(itemp,icomp,jcomp,isppol)*elph_tr_ds%dos_n0(itemp,isppol)*tmp_veloc_sq0)
         write(unit_rho,'(4D20.10)')temp,rho,rho*chgu,rho/temp
         write(unit_tau,'(3D20.10)')temp,tau,tau*chtu
         rho_T(itemp)=rho
       end do ! temperature = 1?
       write(unit_rho,*)
       write(unit_tau,*)

     end do ! jcomp = 3
   end do ! icomp = 3
 end do ! isppol = 2

!-----------------------------

 seebeck = zero
!rho_nm  = zero

!do isppol=1,elph_ds%nsppol
!do itemp=1,ntemper
!q11_inv(:,:)=q11(itemp,:,:,isppol)
!call matrginv(q11_inv,3,3)
!do icomp=1,3
!do jcomp=1,3
!do kcomp=1,3
!seebeck(isppol,itemp,icomp,jcomp) = seebeck(isppol,itemp,icomp,jcomp) + &
!&                             q01(itemp,icomp,kcomp,isppol)*q11_inv(kcomp,jcomp)
!end do
!end do
!end do
!end do
!end do

 do isppol=1,elph_ds%nsppol
   do itemp=1,ntemper
     q11_inv(:,:)=q11(itemp,:,:,isppol)
     call matrginv(q11_inv,3,3)
     call DGEMM('N','N',3,3,3,one,q01(itemp,:,:,isppol),3,q11_inv,&
&     3,zero,seebeck(isppol,itemp,:,:),3)
!    call DGEMM('N','N',3,3,3,one,seebeck(isppol,itemp,:,:),3,q01(itemp,:,:,isppol),&
!    &     3,zero,rho_nm(isppol,itemp,:,:),3)
   end do
 end do
 pref_s = pi/sqrt(3.0_dp)
 seebeck=pref_s*seebeck

!fullq = zero
!do icomp=1,3
!do jcomp=1,3
!fullq(icomp,jcomp) = q00(1,icomp,jcomp,1)
!end do
!end do
!do icomp=1,3
!do jcomp=4,6
!fullq(icomp,jcomp) = q01(1,icomp,jcomp-3,1)
!end do
!end do
!do icomp=4,6
!do jcomp=1,3
!fullq(icomp,jcomp) = q01(1,icomp-3,jcomp,1)
!end do
!end do
!do icomp=4,6
!do jcomp=4,6
!fullq(icomp,jcomp) = q11(1,icomp-3,jcomp-3,1)
!end do
!end do
!write(*,*)' fullq'
!write(*,"(6E20.12)") (fullq(1,jcomp),jcomp=1,6)
!write(*,"(6E20.12)") (fullq(2,jcomp),jcomp=1,6)
!write(*,"(6E20.12)") (fullq(3,jcomp),jcomp=1,6)
!write(*,"(6E20.12)") (fullq(4,jcomp),jcomp=1,6)
!write(*,"(6E20.12)") (fullq(5,jcomp),jcomp=1,6)
!write(*,"(6E20.12)") (fullq(6,jcomp),jcomp=1,6)
 write(message,'(a)') 'q00:'
 call wrtout(std_out,message,'COLL')
 write(message,'(3E20.12)') (q00(1,1,jcomp,1),jcomp=1,3)
 call wrtout(std_out,message,'COLL')
 write(message,'(3E20.12)') (q00(1,2,jcomp,1),jcomp=1,3)
 call wrtout(std_out,message,'COLL')
 write(message,'(3E20.12)') (q00(1,3,jcomp,1),jcomp=1,3)
 call wrtout(std_out,message,'COLL')
 write(message,'(a)') 'q01:'
 call wrtout(std_out,message,'COLL')
 write(message,'(3E20.12)') (q01(1,1,jcomp,1),jcomp=1,3)
 call wrtout(std_out,message,'COLL')
 write(message,'(3E20.12)') (q01(1,2,jcomp,1),jcomp=1,3)
 call wrtout(std_out,message,'COLL')
 write(message,'(3E20.12)') (q01(1,3,jcomp,1),jcomp=1,3)
 call wrtout(std_out,message,'COLL')
 write(message,'(a)') 'q11, q11_inv:'
 call wrtout(std_out,message,'COLL')
 write(message,'(6E20.12)') (q11(1,1,jcomp,1),jcomp=1,3),(q11_inv(1,jcomp),jcomp=1,3)
 call wrtout(std_out,message,'COLL')
 write(message,'(6E20.12)') (q11(1,2,jcomp,1),jcomp=1,3),(q11_inv(2,jcomp),jcomp=1,3)
 call wrtout(std_out,message,'COLL')
 write(message,'(6E20.12)') (q11(1,3,jcomp,1),jcomp=1,3),(q11_inv(3,jcomp),jcomp=1,3)
 call wrtout(std_out,message,'COLL')
!q11_inv = zero
!do icomp = 1, 3
!q11_inv(icomp,icomp) = 2.0_dp
!end do

!call matrginv(fullq,6,6)

!do isppol=1,elph_ds%nsppol
!do itemp=1,ntemper
!rho_nm(isppol,itemp,:,:) = q00(itemp,:,:,isppol) - rho_nm(isppol,itemp,:,:)
!end do
!end do
!rho_nm = 0.5_dp*rho_nm

!Output of Seebeck coefficient
 do isppol=1,elph_ds%nsppol
   do icomp=1,3
     do jcomp=1,3
       itrtensor=(icomp-1)*3+jcomp
       write(unit_sbk,*) '# Seebeck for isppol, itrten = ', isppol, itrtensor
!      write(88,*) '# Rho for isppol, itrten = ', isppol, itrtensor
!      write(89,*) '# Rho for isppol, itrten = ', isppol, itrtensor
       do itemp=1,ntemper
         Temp=tempermin+temperinc*dble(itemp)
         write(unit_sbk,'(3D20.10)')temp, seebeck(isppol,itemp,icomp,jcomp), seebeck(isppol,itemp,icomp,jcomp)*chsu
!        write(88,'(3D20.10)')temp, rho_nm(isppol,itemp,icomp,jcomp), rho_nm(isppol,itemp,icomp,jcomp)*chgu
!        write(89,'(3D20.10)')temp, 0.5_dp/fullq(1,1), 0.5_dp*chgu/fullq(1,1)
       end do ! temperature
       write(unit_sbk,*)
!      write(88,*)
!      write(89,*)
     end do ! jcomp
   end do ! icomp
 end do ! isppol

!Get thermal resistivity, based on eqn. (52) in Allen's PRB 17, 3725 (1978) [[cite:Allen1978]]
!WARNING: before 6.13.1 the thermal resistivity and Lorentz number were not in
!atomic units, BUT the SI units are good.
 pref_w = 3.0_dp/(2.0_dp*pi**2.0d0)
 do isppol=1,elph_ds%nsppol
   do icomp=1, 3
     do jcomp=1, 3
       itrtensor=(icomp-1)*3+jcomp

       write(unit_therm,*) '# Thermal resistivity for isppol, itrten= ', isppol
       write(unit_lor,*) '# Lorentz coefficient for isppol, itrten= ', isppol

       do itemp=1,ntemper

         Temp=tempermin + temperinc*dble(itemp)

         wtherm = pref_w*q11(itemp,icomp,jcomp,isppol)/(kb_HaK*Temp)

!        write(unit_therm,'(5D20.10)')temp,wtherm,1./wtherm,wtherm/3.4057d9,1./(wtherm) *3.4057d9
         write(unit_therm,'(5D20.10)')temp,wtherm,1.0_dp/wtherm,wtherm*chwu,1.0_dp/(wtherm*chwu)

         lorentz=rho_T(itemp)/(wtherm*kb_HaK*Temp)
         write(unit_lor,*)temp,lorentz,lor0

       end do
       write(unit_therm,*)
       write(unit_lor,*)
     end do ! jcomp
   end do ! icomp
 end do ! isppol


 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(displ)
 ABI_DEALLOCATE(pheigvec)
 ABI_DEALLOCATE(integrand_q00)
 ABI_DEALLOCATE(integrand_q01)
 ABI_DEALLOCATE(integrand_q11)
 ABI_DEALLOCATE(q00)
 ABI_DEALLOCATE(q01)
 ABI_DEALLOCATE(q11)
 ABI_DEALLOCATE(seebeck)
 ABI_DEALLOCATE(rho_T)
 ABI_DEALLOCATE(integrho)
 ABI_DEALLOCATE(tointegrho)

 close (unit=unit_lor)
 close (unit=unit_rho)
 close (unit=unit_tau)
 close (unit=unit_sbk)
 close (unit=unit_therm)

 ABI_DEALLOCATE(elph_ds%k_fine%wtq)
 ABI_DEALLOCATE(elph_ds%k_phon%wtq)

 ABI_DEALLOCATE(elph_tr_ds%a2f_1d_tr)

 ABI_DEALLOCATE(elph_tr_ds%gamma_qpt_tr)
 ABI_DEALLOCATE(elph_tr_ds%gamma_rpt_tr)
 write(std_out,*) ' mka2f_tr : end '

end subroutine mka2f_tr
!!***


!!****f* ABINIT/mka2f_tr_lova
!!
!! NAME
!! mka2f_tr_lova
!!
!! FUNCTION
!!  calculates the FS averaged Transport alpha^2F_tr alpha^2F_trout alpha^2F_trin functions
!!  calculates and outputs the associated electrical and thermal conductivities
!!  for the first task: copied from mka2F
!!
!! INPUTS
!! crystal<crystal_t>=data type gathering info on the crystalline structure.
!! Ifc<ifc_type>=Object containing the interatomic force constants.
!!  elph_ds
!!    elph_ds%gkk2 = gkk2 matrix elements on full FS grid for each phonon mode
!!    elph_ds%nbranch = number of phonon branches = 3*natom
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%k_fine%nkpt = number of kpts included in the FS integration
!!    elph_ds%k_fine%wtk = integration weights on the FS
!!    delph_ds%n0 = DOS at the Fermi level calculated from the k_fine integration weights
!!    elph_ds%k_fine%kpt = coordinates of all FS kpoints
!!  mustar = coulomb pseudopotential parameter
!!       eventually for 2 spin channels
!!  ntemper = number of temperature points to calculate, from tempermin to tempermin+ntemper*temperinc
!!  tempermin = minimum temperature at which resistivity etc are calculated (in K)
!!  temperinc = interval for temperature grid on which resistivity etc are calculated (in K)
!!
!! OUTPUT
!!  elph_ds
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      ftgam,ftgam_init,gam_mult_displ,ifc_fourq,simpson_int,wrtout,zgemm
!!
!! NOTES
!!   copied from ftiaf9.f
!!
!! SOURCE

subroutine mka2f_tr_lova(crystal,ifc,elph_ds,ntemper,tempermin,temperinc,elph_tr_ds)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ntemper
 real(dp),intent(in) :: tempermin,temperinc
 type(crystal_t),intent(in) :: crystal
 type(ifc_type),intent(in) :: ifc
 type(elph_tr_type),intent(inout) :: elph_tr_ds
 type(elph_type),intent(inout) :: elph_ds

!Local variables -------------------------
!x =w/(2kbT)
!scalars
 integer :: iFSqpt,ibranch,iomega,isppol,jbranch,nerr
 integer :: unit_a2f_tr, unit_a2f_trout, unit_a2f_trin, natom
 integer :: idir, iatom, k1, kdir,unit_lor,unit_rho,unit_tau,unit_therm
 integer :: itemp,nrpt,itrtensor, icomp, jcomp
 real(dp) :: Temp,chgu,chtu,femto,diagerr,firh,firhT,gaussfactor,domega
 real(dp) :: firh_tau,firhT_tau ! added by BX to get Tau
 real(dp) :: a2fprefactor_in, temp_in
 real(dp) :: a2fprefactor_out, temp_out
 real(dp) :: gaussprefactor,gaussval,lambda_tr,lor0,lorentz,maxerr,maxx,omega
 real(dp) :: rho,tau,tolexp,wtherm,xtr,xx
 real(dp) :: lambda_tr_trace,omega_min, omega_max,qnorm2,spinfact
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 real(dp),parameter :: c0(2)=(/0.d0,0.d0/),c1(2)=(/1.d0,0.d0/)
 real(dp) :: gprimd(3,3)
 real(dp) :: eigval_in(elph_ds%nbranch)
 real(dp) :: eigval_out(elph_ds%nbranch)
 real(dp) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: gam_now_in (2,elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: gam_now_out(2,elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: tmpa2f_in (elph_ds%na2f)
 real(dp) :: tmpa2f_out(elph_ds%na2f)
 real(dp) :: tmpgam1(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmpgam2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),allocatable :: phfrq(:,:)
 real(dp),allocatable :: displ(:,:,:,:)
 real(dp),allocatable :: pheigvec(:,:)
 real(dp),allocatable :: integrho(:),integtau(:),tointegrho(:),tointega2f(:),tointegtau(:)
 real(dp),allocatable :: rho_T(:),tau_T(:)
 real(dp),allocatable :: coskr(:,:)
 real(dp),allocatable :: sinkr(:,:)

! *********************************************************************

!calculate a2f_tr for frequencies between 0 and omega_max
 write(std_out,*) 'mka2f_tr_lova : enter '
!
!MG: the step should be calculated locally using nomega and the extrema of the spectrum.
!One should not rely on previous calls for the setup of elph_ds%domega
!I will remove elph_ds%domega since mka2f.F90 will become a method of gamma_t
 domega =elph_ds%domega

 ! Number of points for FFT interpolation
 nrpt = ifc%nrpt
 natom = crystal%natom
 gprimd = crystal%gprimd

 ABI_ALLOCATE(elph_tr_ds%a2f_1d_tr,(elph_ds%na2f,9,elph_ds%nsppol,1,1,1))
 ABI_ALLOCATE(elph_tr_ds%a2f_1d_trin,(elph_ds%na2f,9,elph_ds%nsppol))
 ABI_ALLOCATE(elph_tr_ds%a2f_1d_trout,(elph_ds%na2f,9,elph_ds%nsppol))

!! defaults for number of temperature steps and max T (all in Kelvin...)
!ntemper=1000
!tempermin=zero
!temperinc=one
 ABI_ALLOCATE(rho_T,(ntemper))
 ABI_ALLOCATE(tau_T,(ntemper))


!tolerance on gaussian being = 0
 tolexp = 1.d-100
 maxx = sqrt(-log(tolexp))
 lor0=(pi*kb_HaK)**2/3.

!maximum value of frequency (a grid has to be chosen for the representation of alpha^2 F)
!WARNING! supposes this value has been set in mkelph_linwid.

 gaussprefactor = sqrt(piinv) / elph_ds%a2fsmear
 gaussfactor = one / elph_ds%a2fsmear

!spinfact should be 1 for a normal non sppol calculation without spinorbit
!for spinors it should also be 1 as bands are twice as numerous but n0 has been divided by 2
!for sppol 2 it should be 0.5 as we have 2 spin channels to sum
 spinfact = one / elph_ds%nsppol !/ elph_ds%nspinor

!ENDMG

 elph_tr_ds%a2f_1d_tr = zero
 elph_tr_ds%a2f_1d_trin = zero
 elph_tr_ds%a2f_1d_trout = zero

 maxerr=0.
 nerr=0

 ABI_ALLOCATE(phfrq,(elph_ds%nbranch, elph_ds%k_fine%nkpt))
 ABI_ALLOCATE(displ,(2, elph_ds%nbranch, elph_ds%nbranch, elph_ds%k_fine%nkpt))
 ABI_ALLOCATE(pheigvec,(2*elph_ds%nbranch*elph_ds%nbranch, elph_ds%k_fine%nkpt))

 do iFSqpt=1,elph_ds%k_fine%nkpt
   call ifc%fourq(crystal,elph_ds%k_fine%kpt(:,iFSqpt),phfrq(:,iFSqpt),displ(:,:,:,iFSqpt),out_eigvec=pheigvec(:,iFSqpt))
 end do

 omega_min = minval(phfrq)
 omega_max = maxval(phfrq)

 ABI_ALLOCATE(coskr, (elph_ds%k_fine%nkpt,nrpt))
 ABI_ALLOCATE(sinkr, (elph_ds%k_fine%nkpt,nrpt))
 call ftgam_init(Ifc%gprim, elph_ds%k_fine%nkpt, nrpt, elph_ds%k_fine%kpt, Ifc%rpt, coskr, sinkr)

 do isppol=1,elph_ds%nsppol

!  loop over qpoint in full kpt grid (presumably dense)
   do iFSqpt=1,elph_ds%k_fine%nkpt
     qnorm2 = sum(elph_ds%k_fine%kpt(:,iFSqpt)**2)
!    if (flag_to_exclude_soft_modes = .false.) qnorm2 = zero
     do itrtensor=1,9
!      Do FT from real-space gamma grid to 1 qpt.

       if (elph_ds%ep_int_gkk == 1) then
         gam_now_in(:,:) = elph_tr_ds%gamma_qpt_trin(:,itrtensor,:,isppol,iFSqpt)
         gam_now_out(:,:) = elph_tr_ds%gamma_qpt_trout(:,itrtensor,:,isppol,iFSqpt)
       else
         call ftgam(Ifc%wghatm,gam_now_in, elph_tr_ds%gamma_rpt_trin(:,itrtensor,:,isppol,:),natom,1,nrpt,0,&
&         coskr(iFSqpt,:), sinkr(iFSqpt,:))
         call ftgam(Ifc%wghatm,gam_now_out,elph_tr_ds%gamma_rpt_trout(:,itrtensor,:,isppol,:),natom,1,nrpt,0,&
&         coskr(iFSqpt,:), sinkr(iFSqpt,:))
       end if

!      Diagonalize gamma matrix at this qpoint (complex matrix).

!      if ep_scalprod==0 we have to dot in the displacement vectors here
       if (elph_ds%ep_scalprod==0) then

         displ_red(:,:,:) = zero
         do jbranch=1,elph_ds%nbranch
           do iatom=1,natom
             do idir=1,3
               ibranch=idir+3*(iatom-1)
               do kdir=1,3
                 k1 = kdir+3*(iatom-1)
                 displ_red(1,ibranch,jbranch) = displ_red(1,ibranch,jbranch) + &
&                 gprimd(kdir,idir)*displ(1,k1,jbranch,iFSqpt)
                 displ_red(2,ibranch,jbranch) = displ_red(2,ibranch,jbranch) + &
&                 gprimd(kdir,idir)*displ(2,k1,jbranch,iFSqpt)
               end do
             end do
           end do
         end do

         tmpgam2 = reshape (gam_now_in, (/2,elph_ds%nbranch,elph_ds%nbranch/))
         call gam_mult_displ(elph_ds%nbranch, displ_red, tmpgam2, tmpgam1)
         do jbranch=1,elph_ds%nbranch
           eigval_in(jbranch)   = tmpgam1(1, jbranch, jbranch)
         end do

         tmpgam2 = reshape (gam_now_out, (/2,elph_ds%nbranch,elph_ds%nbranch/))
         call gam_mult_displ(elph_ds%nbranch, displ_red, tmpgam2, tmpgam1)
         do jbranch=1,elph_ds%nbranch
           eigval_out(jbranch)   = tmpgam1(1, jbranch, jbranch)
         end do

       else if (elph_ds%ep_scalprod == 1) then

!
!        NOTE: in these calls gam_now and pheigvec do not have the right rank, but blas usually does not care
         call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, c1, gam_now_in, 3*natom,&
&         pheigvec(:,iFSqpt), 3*natom, c0, tmpgam1, 3*natom)
         call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, c1, pheigvec(:,iFSqpt), 3*natom,&
&         tmpgam1, 3*natom, c0, tmpgam2, 3*natom)
         diagerr = zero

         do ibranch=1,elph_ds%nbranch
           eigval_in(ibranch) = tmpgam2(1,ibranch,ibranch)
           do jbranch=1,ibranch-1
             diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
           end do
           do jbranch=ibranch+1,elph_ds%nbranch
             diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
           end do
         end do
         if (diagerr > tol12) then
           nerr=nerr+1
           maxerr=max(diagerr, maxerr)
         end if

         call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, c1, gam_now_out, 3*natom,&
&         pheigvec(:,iFSqpt), 3*natom, c0, tmpgam1, 3*natom)
         call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, c1, pheigvec(:,iFSqpt), 3*natom,&
&         tmpgam1, 3*natom, c0, tmpgam2, 3*natom)
         diagerr = zero

         do ibranch=1,elph_ds%nbranch
           eigval_out(ibranch) = tmpgam2(1,ibranch,ibranch)
           do jbranch=1,ibranch-1
             diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
           end do
           do jbranch=ibranch+1,elph_ds%nbranch
             diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
           end do
         end do
         if (diagerr > tol12) then
           nerr=nerr+1
           maxerr=max(diagerr, maxerr)
         end if
       end if
!      end ep_scalprod if

!      Add all contributions from the phonon modes at this qpoint to
!      a2f and the phonon dos.
       do ibranch=1,elph_ds%nbranch
!        if (abs(phfrq(ibranch,iFSqpt)) < tol10) then
         if ( abs(phfrq(ibranch,iFSqpt)) < tol7 .or. &
&         (phfrq(ibranch,iFSqpt) < tol4 .and. qnorm2 > 0.03 )) then !
!          note: this should depend on the velocity of sound, to accept acoustic
!          modes!
           a2fprefactor_in = zero
           a2fprefactor_out= zero
         else
           a2fprefactor_in  = eigval_in (ibranch)/(two_pi*abs(phfrq(ibranch,iFSqpt))*elph_ds%n0(isppol))
           a2fprefactor_out = eigval_out(ibranch)/(two_pi*abs(phfrq(ibranch,iFSqpt))*elph_ds%n0(isppol))
         end if

         omega = omega_min
         tmpa2f_in (:) = zero
         tmpa2f_out(:) = zero
         do iomega=1,elph_ds%na2f
           xx = (omega-phfrq(ibranch,iFSqpt))*gaussfactor
           gaussval = gaussprefactor*exp(-xx*xx)

           temp_in = gaussval*a2fprefactor_in
           temp_out = gaussval*a2fprefactor_out

           if (dabs(temp_in) < 1.0d-50) temp_in = zero
           if (dabs(temp_out) < 1.0d-50) temp_out = zero
           tmpa2f_in (iomega) = tmpa2f_in (iomega) + temp_in
           tmpa2f_out(iomega) = tmpa2f_out(iomega) + temp_out
           omega = omega+domega
         end do

         elph_tr_ds%a2f_1d_trin (:,itrtensor,isppol) = elph_tr_ds%a2f_1d_trin (:,itrtensor,isppol) + tmpa2f_in(:)
         elph_tr_ds%a2f_1d_trout(:,itrtensor,isppol) = elph_tr_ds%a2f_1d_trout(:,itrtensor,isppol) + tmpa2f_out(:)

       end do ! end ibranch do
     end do ! end itrtensor do
   end do ! end iFSqpt do
 end do ! end isppol

 ABI_DEALLOCATE(coskr)
 ABI_DEALLOCATE(sinkr)

!second 1 / elph_ds%k_fine%nkpt factor for the integration weights
 elph_tr_ds%a2f_1d_trin  = elph_tr_ds%a2f_1d_trin  / elph_ds%k_fine%nkpt
 elph_tr_ds%a2f_1d_trout = elph_tr_ds%a2f_1d_trout / elph_ds%k_fine%nkpt

 if (elph_ds%ep_scalprod == 1) then
   write(std_out,*) 'mka2f_tr_lova: errors in diagonalization of gamma_tr with phon eigenvectors: ', nerr,maxerr
 end if

 elph_tr_ds%a2f_1d_tr(:,:,:,1,1,1) = elph_tr_ds%a2f_1d_trout(:,:,:) - elph_tr_ds%a2f_1d_trin(:,:,:)

!output the elph_tr_ds%a2f_1d_tr
 fname = trim(elph_ds%elph_base_name) // '_A2F_TR'
 if (open_file (fname,message,newunit=unit_a2f_tr,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

 fname = trim(elph_ds%elph_base_name) // '_A2F_TRIN'
 if (open_file(fname,message,newunit=unit_a2f_trin,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

 fname = trim(elph_ds%elph_base_name) // '_A2F_TROUT'
 if (open_file (fname,message,newunit=unit_a2f_trout,status='unknown') /=0) then
   MSG_ERROR(message)
 end if

 write (unit_a2f_tr,'(a)')       '#'
 write (unit_a2f_tr,'(a)')       '# ABINIT package : a2f_tr file'
 write (unit_a2f_tr,'(a)')       '#'
 write (unit_a2f_tr,'(a)')       '# a2f_tr function integrated over the FS. omega in a.u.'
 write (unit_a2f_tr,'(a,I10)')   '#     number of kpoints integrated over : ', elph_ds%k_fine%nkpt
 write (unit_a2f_tr,'(a,I10)')   '#     number of energy points : ',elph_ds%na2f
 write (unit_a2f_tr,'(a,E16.6,a,E16.6,a)') '#       between omega_min = ', omega_min,' Ha and omega_max = ', omega_max, ' Ha'
 write (unit_a2f_tr,'(a,E16.6)') '#   and the smearing width for gaussians is ', elph_ds%a2fsmear
 write (unit_a2f_tr,'(a)')       '#'

 write (unit_a2f_trin,'(a)')       '#'
 write (unit_a2f_trin,'(a)')       '# ABINIT package : a2f_trin file'
 write (unit_a2f_trin,'(a)')       '#'
 write (unit_a2f_trin,'(a)')       '# a2f_trin function integrated over the FS. omega in a.u.'
 write (unit_a2f_trin,'(a,I10)')   '#     number of kpoints integrated over : ', elph_ds%k_fine%nkpt
 write (unit_a2f_trin,'(a,I10)')   '#     number of energy points : ',elph_ds%na2f
 write (unit_a2f_trin,'(a,E16.6,a,E16.6,a)') '#       between omega_min = ', omega_min,' Ha and omega_max = ', omega_max, ' Ha'
 write (unit_a2f_trin,'(a,E16.6)') '#   and the smearing width for gaussians is ', elph_ds%a2fsmear
 write (unit_a2f_trin,'(a)')       '#'

 write (unit_a2f_trout,'(a)')       '#'
 write (unit_a2f_trout,'(a)')       '# ABINIT package : a2f_trout file'
 write (unit_a2f_trout,'(a)')       '#'
 write (unit_a2f_trout,'(a)')       '# a2f_trout function integrated over the FS. omega in a.u.'
 write (unit_a2f_trout,'(a,I10)')   '#     number of kpoints integrated over : ', elph_ds%k_fine%nkpt
 write (unit_a2f_trout,'(a,I10)')   '#     number of energy points : ',elph_ds%na2f
 write (unit_a2f_trout,'(a,E16.6,a,E16.6,a)') '#       between omega_min = ', omega_min,' Ha and omega_max = ', omega_max, ' Ha'
 write (unit_a2f_trout,'(a,E16.6)') '#   and the smearing width for gaussians is ', elph_ds%a2fsmear
 write (unit_a2f_trout,'(a)')       '#'

!done with header
 do isppol=1,elph_ds%nsppol
   write (unit_a2f_tr,'(a,E16.6)') '# The DOS at Fermi level is ', elph_ds%n0(isppol)
   write (unit_a2f_trin,'(a,E16.6)') '# The DOS at Fermi level is ', elph_ds%n0(isppol)
   write (unit_a2f_trout,'(a,E16.6)') '# The DOS at Fermi level is ', elph_ds%n0(isppol)
!  omega = zero
   omega = omega_min
   do iomega=1,elph_ds%na2f
     write (unit_a2f_tr,   '(10D16.6)') omega, elph_tr_ds%a2f_1d_tr   (iomega,:,isppol,1,1,1)
     write (unit_a2f_trin, '(10D16.6)') omega, elph_tr_ds%a2f_1d_trin (iomega,:,isppol)
     write (unit_a2f_trout,'(10D16.6)') omega, elph_tr_ds%a2f_1d_trout(iomega,:,isppol)
     omega=omega+domega
   end do
   write (unit_a2f_tr,*)
   write (unit_a2f_trin,*)
   write (unit_a2f_trout,*)
 end do !isppol

 close (unit=unit_a2f_tr)
 close (unit=unit_a2f_trin)
 close (unit=unit_a2f_trout)

!calculation of transport properties
 ABI_ALLOCATE(integrho,(elph_ds%na2f))
 ABI_ALLOCATE(tointegrho,(elph_ds%na2f))
 ABI_ALLOCATE(tointega2f,(elph_ds%na2f))
 ABI_ALLOCATE(integtau,(elph_ds%na2f))
 ABI_ALLOCATE(tointegtau,(elph_ds%na2f))

 fname = trim(elph_ds%elph_base_name) // '_RHO'
 if (open_file(fname,message,newunit=unit_rho,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to resistivity file
 write (unit_rho,*) '# Resistivity as a function of temperature.'
 write (unit_rho,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_rho,*) '#  '
 write (unit_rho,*) '#  Columns are: '
 write (unit_rho,*) '#  temperature[K]   rho[au]   rho [SI]        rho/temp [au]'
 write (unit_rho,*) '#  '

 fname = trim(elph_ds%elph_base_name) // '_TAU'
 if (open_file(fname,message,newunit=unit_tau,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to relaxation time file
 write (unit_tau,*) '# Relaxation time as a function of temperature.'
 write (unit_tau,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_tau,*) '#  '
 write (unit_tau,*) '#  Columns are: '
 write (unit_tau,*) '#  temperature[K]   tau[au]   tau [femtosecond]     '
 write (unit_tau,*) '#  '

 fname = trim(elph_ds%elph_base_name) // '_WTH'
 if (open_file(fname,message,newunit=unit_therm,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to thermal conductivity file
 write (unit_therm,'(a)') '# Thermal conductivity/resistivity as a function of temperature.'
 write (unit_therm,'(a)') '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_therm,'(a)') '#  '
 write (unit_therm,'(a)') '#  Columns are: '
 write (unit_therm,'(a)') '#  temperature[K]   thermal rho[au]   thermal cond [au]   thermal rho [SI]   thermal cond [SI]'
 write (unit_therm,'(a)') '#  '

 fname = trim(elph_ds%elph_base_name) // '_LOR'
 if (open_file(fname,message,newunit=unit_lor,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to lorentz file
 write (unit_lor,*) '# Lorentz number as a function of temperature.'
 write (unit_lor,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_lor,*) '#  '
 write (unit_lor,*) '#  Columns are: '
 write (unit_lor,*) '#  temperature[K]   Lorentz number[au]   Lorentz quantum = (pi*kb_HaK)**2/3'
 write (unit_lor,*) '#  '

 do isppol=1,elph_ds%nsppol
   lambda_tr_trace = zero
   do itrtensor=1,9
     omega = omega_min
     tointega2f = zero
     do iomega=1,elph_ds%na2f
       if(omega<=0) then
         omega=omega+domega
         cycle
       end if
       tointega2f(iomega)=elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,1,1,1)/omega
       omega=omega+domega
     end do

     integrho = zero
     call simpson_int(elph_ds%na2f,domega,tointega2f,integrho)
     lambda_tr = two * spinfact * integrho(elph_ds%na2f)
     write (message, '(a,2i3,a,es16.6)' )&
&     ' mka2f_tr_lova : TRANSPORT lambda for isppol itrtensor', isppol, itrtensor, ' =  ', lambda_tr
     call wrtout(std_out,message,'COLL')
     if (itrtensor == 1 .or. itrtensor == 5 .or. itrtensor == 9) lambda_tr_trace = lambda_tr_trace + lambda_tr
   end do !end itrtensor do

   lambda_tr_trace = lambda_tr_trace / three
   write (message, '(a,i3,a,es16.6)' )&
&   ' mka2f_tr_lova: 1/3 trace of TRANSPORT lambda for isppol ', isppol, ' =  ', lambda_tr_trace
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 end do !end isppol do

!constant to change units of rho from au to SI
 chgu=2.173969d-7
 chtu=2.4188843265d-17
 femto=1.0d-15

 do isppol=1,elph_ds%nsppol
   do icomp=1, 3
     do jcomp=1, 3
       itrtensor=(icomp-1)*3+jcomp

!      prefactor for resistivity integral
!      firh=6.d0*pi*crystal%ucvol*kb_HaK/(elph_ds%n0(isppol)*elph_tr_ds%FSelecveloc_sq(isppol))
!      FIXME: check factor of 2 which is different from Savrasov paper. 6 below for thermal conductivity is correct.
       firh=2.d0*pi*crystal%ucvol*kb_HaK/elph_ds%n0(isppol)/&
&       sqrt(elph_tr_ds%FSelecveloc_sq(icomp,isppol)*elph_tr_ds%FSelecveloc_sq(jcomp,isppol))

!      Add by BX to get Tau_elph
       firh_tau = 2.0d0*pi*kb_HaK
!      End Adding

       write(unit_rho,*) '# Rho for isppol, itrten = ', isppol, itrtensor
       write(unit_tau,*) '# Tau for isppol, itrten = ', isppol, itrtensor

! jmb
       tointegtau(:)=0.
       tointegrho(:)=0.
       do itemp=1,ntemper  ! runs over termperature in K
         Temp=tempermin+temperinc*dble(itemp)
         firhT=firh*Temp
         firhT_tau=firh_tau*Temp
         omega = omega_min
         do iomega=1,elph_ds%na2f
           if(omega<=0) then
             omega=omega+domega
             cycle
           end if
           xtr=omega/(2*kb_HaK*Temp)
           if(xtr < log(huge(zero)*tol16)/2)then
             tointegrho(iomega)=spinfact*firhT*omega*elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,1,1,1)  &
&             /(((2*Temp*kb_HaK)**2)*((exp(xtr)-exp(-xtr))/2)**2)
!            Add by BX to get Tau
             tointegtau(iomega)=spinfact*firhT_tau*omega*elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,1,1,1)  &
&             /(((2*Temp*kb_HaK)**2)*((exp(xtr)-exp(-xtr))/2)**2)
           else
             tointegrho(iomega)=zero
             tointegtau(iomega)=zero
           end if
           omega=omega+domega
         end do

         call simpson_int(elph_ds%na2f,domega,tointegrho,integrho)
         call simpson_int(elph_ds%na2f,domega,tointegtau,integtau)
         rho=integrho(elph_ds%na2f)
         tau=1.0d99
         if(dabs(integtau(elph_ds%na2f)) < tol7) then
           write(message,'(a)') ' Cannot get a physical relaxation time '
           MSG_WARNING(message)
         else
           tau=1.0d0/integtau(elph_ds%na2f)
         end if
!         if(elph_ds%na2f < 350.0) then
!           tau=1.0d0/integtau(elph_ds%na2f)
!         end if
         write(unit_rho,'(4D20.10)')temp,rho,rho*chgu,rho/temp
         write(unit_tau,'(3D20.10)')temp,tau,tau*chtu/femto
         rho_T(itemp)=rho
         tau_T(itemp)=tau
       end do ! temperature
       write(unit_rho,*)
       write(unit_tau,*)

     end do ! jcomp
   end do ! icomp
 end do ! isppol

!-----------------------------


 do isppol=1,elph_ds%nsppol
   do icomp=1, 3
     do jcomp=1, 3
       itrtensor=(icomp-1)*3+jcomp
!      prefactor for integral of thermal conductivity
!      firh=(18.*crystal%ucvol)/(pi*kb_HaK*elph_ds%n0(isppol)*elph_tr_ds%FSelecveloc_sq(isppol))
       firh=(6.d0*crystal%ucvol)/(pi*kb_HaK*elph_ds%n0(isppol))/ &
&       sqrt(elph_tr_ds%FSelecveloc_sq(icomp,isppol)*elph_tr_ds%FSelecveloc_sq(jcomp,isppol))


       write(unit_therm,*) '# Thermal resistivity for isppol, itrten= ', isppol
       write(unit_lor,*) '# Lorentz coefficient for isppol, itrten= ', isppol

       tointegrho(:)=0.
       do itemp=1,ntemper

         Temp=tempermin + temperinc*dble(itemp)
         omega = omega_min
         do iomega=1,elph_ds%na2f
           if(omega<=0) then
             omega=omega+domega
             cycle
           end if
           xtr=omega/(2*kb_HaK*Temp)
           if(xtr < log(huge(zero)*tol16)/2)then
             tointegrho(iomega) = spinfact*xtr**2/omega*&
&             ( elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol,1,1,1)+&
&             4*xtr**2*elph_tr_ds%a2f_1d_trout(iomega,itrtensor,isppol)/pi**2+   &
&             2*xtr**2*elph_tr_ds%a2f_1d_trin(iomega,itrtensor,isppol)/pi**2)  &
&             /(((exp(xtr)-exp(-xtr))/2)**2)
           else
             tointegrho(iomega) = zero
           end if
           omega=omega+domega
         end do

         call simpson_int(elph_ds%na2f,domega,tointegrho,integrho)
         wtherm=integrho(elph_ds%na2f)*firh

         if(abs(wtherm) > tol12)then
           write(unit_therm,'(5D20.10)')temp,wtherm,1./wtherm,wtherm/3.4057d9,1./(wtherm) *3.4057d9

           lorentz=rho_T(itemp)/(wtherm*temp)
           write(unit_lor,*)temp,lorentz,lor0
         else
           write(unit_therm,'(5D20.10)')temp,zero,huge(one),zero,huge(one)
           write(unit_lor,*)temp,huge(one),lor0
         end if

       end do
       write(unit_therm,*)
       write(unit_lor,*)
     end do ! jcomp
   end do ! icomp
 end do !end isppol do


 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(displ)
 ABI_DEALLOCATE(pheigvec)
 ABI_DEALLOCATE(rho_T)
 ABI_DEALLOCATE(tau_T)

 close (unit=unit_lor)
 close (unit=unit_rho)
 close (unit=unit_tau)
 close (unit=unit_therm)

 ABI_DEALLOCATE(integrho)
 ABI_DEALLOCATE(integtau)
 ABI_DEALLOCATE(tointega2f)
 ABI_DEALLOCATE(tointegrho)
 ABI_DEALLOCATE(tointegtau)
 ABI_DEALLOCATE(elph_tr_ds%a2f_1d_tr)
 ABI_DEALLOCATE(elph_tr_ds%a2f_1d_trin)
 ABI_DEALLOCATE(elph_tr_ds%a2f_1d_trout)

 ABI_DEALLOCATE(elph_tr_ds%gamma_qpt_trin)
 ABI_DEALLOCATE(elph_tr_ds%gamma_qpt_trout)
 ABI_DEALLOCATE(elph_tr_ds%gamma_rpt_trin)
 ABI_DEALLOCATE(elph_tr_ds%gamma_rpt_trout)

!DEBUG
 write(std_out,*) ' mka2f_tr_lova : end '
!ENDDEBUG

end subroutine mka2f_tr_lova
!!***

!!****f* ABINIT/get_tau_k
!! NAME
!!  get_tau_k
!!
!! FUNCTION
!!  Calculate the k-dependent relaxation time due to EPC. Impelementation based
!!  on derivation from Grmvall's book or
!!  OD Restrepo's paper (PRB 94 212103 (2009) [[cite:Restrepo2009]])
!!
!! INPUTS
!!  Cryst<crystal_t>=Info on the unit cell and on its symmetries.
!!  Ifc<ifc_type>=Object containing the interatomic force constants.
!!  elph_ds = elphon datastructure with data and dimensions
!!  eigenGS = Ground State eigenvalues
!!  max_occ = maximal occupancy for a band
!!
!! OUTPUT
!!  tau_k(nsppol,nkptirr,nband)=mode relaxation time due to electron phonono coupling
!!  rate_e(nene)= scattering rate due to electron phonono coupling vs. energy
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      dgemm,ebands_prtbltztrp_tau_out,ebands_update_occ,ep_el_weights
!!      ep_ph_weights,ftgam,ftgam_init,gam_mult_displ,ifc_fourq,matrginv
!!      mkqptequiv,phdispl_cart2red,spline,splint,wrtout,xmpi_sum,zgemm
!!
!! SOURCE

subroutine get_tau_k(Cryst,ifc,Bst,elph_ds,elph_tr_ds,eigenGS,max_occ)

!Arguments ------------------------------------
 type(crystal_t),intent(in) :: Cryst
 type(ifc_type),intent(in) :: ifc
 type(ebands_t),intent(inout)   :: Bst
 type(elph_type),intent(inout) :: elph_ds
 type(elph_tr_type), intent(inout) :: elph_tr_ds
 real(dp),intent(in) :: max_occ
 real(dp),intent(in) :: eigenGS(elph_ds%nband,elph_ds%k_phon%nkpt,elph_ds%nsppol)

!Local variables-------------------------------
!scalars
 character(len=500) :: message
 character(len=fnlen) :: fname
 integer :: ntemper,nsppol,nbranch,nband,natom
 integer :: nkpt,nqpt,nkptirr,nqptirr,new_nkptirr
 integer :: isppol,iFSkpt,iFSqpt,iqpt,iqpt_fullbz,imqpt_fullbz,ikpt_kpq,ikpt_kmq
 integer :: iband,jband,jpband,jbeff,ibranch,jbranch,itemp
 integer :: irec,ierr,nrpt,ik_this_proc
 integer :: unit_tau,unit_invtau
 integer :: nene,nene_all,iene,iene_fine,unit_taue,unit_mfp
 integer :: icomp,jcomp,itensor
 integer :: ikpt_irr,iomega,unit_cond,unit_therm,unit_sbk
 integer :: nskip,nspline
 real(dp) :: occ_omega,occ_e
 real(dp) :: xx,Temp,therm_factor
 real(dp) :: factor,dfermide
 real(dp) :: e_k,chu_tau,rate_e,mfp_e
 real(dp) :: ene,enemin,enemax,deltaene
 real(dp) :: omega,omega_min,omega_max,domega
 real(dp) :: diagerr
 real(dp) :: chu_mfp,chu_cond,chu_cth,chu_sbk,femto
 real(dp) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: eigval(elph_ds%nbranch),eigval2(elph_ds%nbranch)
 real(dp) :: imeigval(elph_ds%nbranch)
 real(dp) :: tmp_wtkpq, tmp_wtkmq, tol_wtk
 real(dp) :: yp1,ypn
!arrays
 integer,allocatable :: FSfullpktofull(:,:),mqtofull(:)
 integer,allocatable :: kpttokpt(:,:,:)
 real(dp) :: cond_inv(3,3)
 real(dp),allocatable :: fermie(:)
 real(dp),allocatable :: tmp_eigenGS(:,:,:)
 real(dp),allocatable :: tmp_gkk_qpt(:,:,:),tmp_gkk_rpt(:,:,:),tmp_gkk_kpt(:,:)
 real(dp),allocatable :: tmp_gkk_kpt2(:,:,:), gkk_kpt(:,:,:)
 real(dp),allocatable :: tau_k(:,:,:,:),inv_tau_k(:,:,:,:),tmp_tau_k(:,:,:,:)
 real(dp),allocatable :: phfrq(:,:),pheigvec(:,:)
 real(dp),allocatable :: displ(:,:,:,:)
 real(dp),allocatable :: a2f_2d(:),a2f_2d2(:)
 real(dp),allocatable :: tmp_wtk(:,:,:,:),tmp2_wtk(:),tmp_wtk1(:),tmp_wtk2(:)
 real(dp),allocatable :: ene_pt(:),ene_ptfine(:),ff2(:)
 real(dp),allocatable :: wtq(:,:,:),tmp_wtq(:,:,:),tmp2_wtq(:,:)
 real(dp),allocatable :: dos_e(:,:)
 real(dp),allocatable :: coskr1(:,:),sinkr1(:,:)
 real(dp),allocatable :: coskr2(:,:),sinkr2(:,:)
 real(dp),allocatable :: cond_e(:,:,:,:),cond(:,:,:,:),sbk(:,:,:,:),seebeck(:,:,:,:),cth(:,:,:,:)

! *************************************************************************

 write(std_out,*) 'get_tau_k : enter '

 nrpt = ifc%nrpt
 natom = cryst%natom

 nsppol   = elph_ds%nsppol
 nbranch  = elph_ds%nbranch
 nband    = elph_ds%ngkkband
 nkpt     = elph_ds%k_phon%nkpt
 nqpt     = elph_ds%nqpt_full
 nkptirr  = elph_ds%k_phon%nkptirr
 new_nkptirr  = elph_ds%k_phon%new_nkptirr
 nqptirr  = elph_ds%nqptirred
 ntemper  = elph_ds%ntemper
 nene = 2*elph_ds%na2f-1 ! only need e_k +- omega_max range, take deltaene=delta_oemga

 chu_tau  = 2.4188843265*1.0d-17
 chu_mfp  = 5.291772*1.0d-11
 chu_cond = 4.59988159904764*1.0d6
 chu_cth  = 1.078637439971599*1.0d4
 chu_sbk  = 8.617343101*1.0d-5
 femto    = 1.0d-15

 tol_wtk = tol7/nkptirr/nband

 ABI_ALLOCATE(fermie ,(ntemper))
 ABI_ALLOCATE(tmp_gkk_qpt ,(2,nbranch**2,nqpt))
 ABI_ALLOCATE(tmp_gkk_rpt ,(2,nbranch**2,nrpt))
 ABI_ALLOCATE(tmp_gkk_kpt ,(2,nbranch**2))
 ABI_ALLOCATE(tmp_gkk_kpt2 ,(2,nbranch,nbranch))
 ABI_ALLOCATE(gkk_kpt ,(2,nbranch,nbranch))
 ABI_ALLOCATE(a2f_2d, (nene))
 ABI_ALLOCATE(a2f_2d2, (nene))
 ABI_ALLOCATE(inv_tau_k, (ntemper,nsppol,nkpt,nband))
 ABI_ALLOCATE(tau_k, (ntemper,nsppol,nkpt,nband))
 ABI_ALLOCATE(tmp_tau_k ,(ntemper,nsppol,new_nkptirr,nband))

 if (elph_ds%gkqwrite == 0) then
   call wrtout(std_out,' get_tau_k : keeping gkq matrices in memory','COLL')
 else if (elph_ds%gkqwrite == 1) then
   fname=trim(elph_ds%elph_base_name) // '_GKKQ'
   write (message,'(2a)')' get_tau_k : reading gkq matrices from file ',trim(fname)
   call wrtout(std_out,message,'COLL')
 else
   write (message,'(a,i0)')' Wrong value for gkqwrite = ',elph_ds%gkqwrite
   MSG_BUG(message)
 end if

!=========================================================
!Get equivalence between a kpt_phon pair and a qpt in qpt_full
!only works if the qpt grid is complete (identical to
!the kpt one, with a basic shift of (0,0,0)
!=========================================================

!mapping of k + q onto k' for k and k' in full BZ
!for dense k grid
 ABI_ALLOCATE(FSfullpktofull,(nkpt,nkpt))
 ABI_ALLOCATE(mqtofull,(nkpt))

!kpttokpt(itim,isym,iqpt) = kpoint index which transforms to ikpt under isym and with time reversal itim.
 ABI_ALLOCATE(kpttokpt,(2,Cryst%nsym,nkpt))

 call wrtout(std_out,'get_tau_k: calling mkqptequiv to set up the FS kpoint set',"COLL")

 call mkqptequiv (FSfullpktofull,Cryst,elph_ds%k_phon%kpt,nkpt,nkpt,kpttokpt,elph_ds%k_phon%kpt,mqtofull)

!=========================================================
!=========================================================

 omega_max       = elph_ds%omega_max
 omega_min       = elph_ds%omega_min
 domega          = elph_ds%domega
 enemax = maxval(eigenGS(elph_ds%maxFSband,:,:))
 enemin = minval(eigenGS(elph_ds%minFSband,:,:))

 if (enemin < (elph_ds%fermie-0.2)) then
   enemin = elph_ds%fermie-0.2
 end if
 if (enemax > (elph_ds%fermie+0.2)) then
   enemax = elph_ds%fermie+0.2
 end if

 nspline = elph_ds%ep_nspline
 nene_all = INT((enemax-enemin+domega)/(nspline*domega)) + 1
 deltaene = domega
 write(std_out,*) 'E_min= ',enemin, 'E_max= ',enemax
 write(std_out,*) 'Number of energy points= ',nene_all
 write(std_out,'(a,I8)') 'scale factor for spline interpolation in RTA = ', elph_ds%ep_nspline
 write(std_out,*) 'delta_ene before spline interpolation= ',deltaene*nspline
 write(std_out,*) 'delta_ene after spline interpolation= ',deltaene
 write(std_out,*) 'Omega_min= ',omega_min, 'Omega_max= ',omega_max
 write(std_out,*) 'Number of phonon points= ',elph_ds%na2f
 write(std_out,*) 'delta_omega= ',domega
 write(std_out,*) 'number of bands= ', elph_ds%nband, nband

 ABI_ALLOCATE(tmp_wtk,(nband,nkpt,nsppol,nene_all))
 ABI_ALLOCATE(tmp2_wtk,(nene_all))
 ABI_ALLOCATE(ff2,(nene_all))
 ABI_ALLOCATE(ene_pt,(nene_all))
 ABI_ALLOCATE(ene_ptfine,(nene_all*nspline))
 ABI_ALLOCATE(tmp_wtk1,(nene_all*nspline))
 ABI_ALLOCATE(tmp_wtk2,(nene_all*nspline))
 ABI_ALLOCATE(dos_e,(nsppol,nene_all))

!Get energy points for spline interpolation
 do iene = 1, nene_all
   ene_pt(iene) = enemin + (iene-1)*nspline*deltaene
 end do

 do iene = 1, nene_all*nspline
   ene_ptfine(iene) = enemin + (iene-1)*deltaene
 end do

 ABI_ALLOCATE(tmp_wtq,(elph_ds%nbranch, elph_ds%k_phon%nkpt, elph_ds%na2f+1))
 ABI_ALLOCATE(wtq,(elph_ds%nbranch, elph_ds%k_phon%nkpt, elph_ds%na2f))
 ABI_ALLOCATE(tmp2_wtq,(elph_ds%nbranch, elph_ds%na2f))

!phonon
 ABI_ALLOCATE(phfrq,(nbranch, nkptirr))
 ABI_ALLOCATE(displ,(2, nbranch, nbranch, nkptirr))
 ABI_ALLOCATE(pheigvec,(2*nbranch*nbranch, nkptirr))

 do iFSqpt = 1, nkptirr
   call ifc%fourq(cryst,elph_ds%k_phon%kptirr(:,iFSqpt),phfrq(:,iFSqpt),displ(:,:,:,iFSqpt),out_eigvec=pheigvec(:,iFSqpt))
 end do

 omega_min = omega_min - domega

!bxu, obtain wtq for the q_fine, then condense to q_phon
 call ep_ph_weights(phfrq,elph_ds%a2fsmear,omega_min,omega_max,elph_ds%na2f+1,Cryst%gprimd,elph_ds%kptrlatt, &
& elph_ds%nbranch,elph_ds%telphint,elph_ds%k_phon,tmp_wtq)
 omega_min = omega_min + domega

 do iomega = 1, elph_ds%na2f
   wtq(:,:,iomega) = tmp_wtq(:,:,iomega+1)
   !write(1005,*) omega_min+(iomega-1)*domega, sum(tmp_wtq(:,:,iomega+1))/nkpt
 end do
 ABI_DEALLOCATE(tmp_wtq)

! electron
 tmp_wtk =zero
 dos_e = zero
 call ep_el_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS(elph_ds%minFSband:elph_ds%minFSband+nband-1,:,:), &
& elph_ds%elphsmear, &
& enemin, enemax, nene_all, Cryst%gprimd, elph_ds%k_phon%irredtoGS, elph_ds%kptrlatt, max_occ, &
& 1, nband, elph_ds%nFSband, nsppol, elph_ds%telphint, elph_ds%k_phon, tmp_wtk)
!& elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, nsppol, elph_ds%telphint, elph_ds%k_phon, tmp_wtk)

 do isppol = 1, nsppol
   do iene = 1, nene_all
     dos_e(isppol,iene) = sum(tmp_wtk(:,:,isppol,iene))/nkpt
   end do
 end do

 ABI_ALLOCATE(coskr1, (nqpt,nrpt))
 ABI_ALLOCATE(sinkr1, (nqpt,nrpt))
 call ftgam_init(ifc%gprim, nqpt, nrpt, elph_ds%k_phon%kpt, Ifc%rpt, coskr1, sinkr1)
 ABI_ALLOCATE(coskr2, (nkptirr,nrpt))
 ABI_ALLOCATE(sinkr2, (nkptirr,nrpt))
 call ftgam_init(ifc%gprim, nkptirr, nrpt, elph_ds%k_phon%kpt, Ifc%rpt, coskr2, sinkr2)

!get fermie for itemp
 fermie = elph_ds%fermie
 do itemp=1,ntemper  ! runs over termperature in K
   Temp=elph_ds%tempermin+elph_ds%temperinc*dble(itemp)

   Bst%occopt = 3
   Bst%tsmear = Temp*kb_HaK
   call ebands_update_occ(Bst,-99.99_dp)
   write(message,'(a,f12.6,a,E20.12)')'At T=',Temp,' Fermi level is:',Bst%fermie
   call wrtout(std_out,message,'COLL')

   if (abs(elph_ds%fermie) < tol10) then
     fermie(itemp) = Bst%fermie
   end if
 end do

 inv_tau_k = zero
!get a2f_2d = \sum_{q,nbranch,jband'} |gkk|^2*\delta(\epsilon_{k'j'}-\epsilon')*\delta(\omega_q-\omega)
 do isppol=1,nsppol
   write (std_out,*) '##############################################'
   write (std_out,*) 'get_tau_k : Treating spin polarization ', isppol
   write (std_out,*) '##############################################'

!   do iFSkpt =1,nkpt
   do ik_this_proc =1,elph_ds%k_phon%my_nkpt
     iFSkpt = elph_ds%k_phon%my_ikpt(ik_this_proc)
     write (std_out,*) 'get_tau_k : working on kpt # ', iFSkpt, '/', nkpt
     do jband = 1, nband
!          write(*,*)'i am here 1 ', isppol,iFSkpt,jband
       a2f_2d = zero
       a2f_2d2 = zero

!sum from here
       nskip = 0
       do jpband = 1, nband
         jbeff = jpband+(jband-1)*nband

         if (elph_ds%gkqwrite == 0) then
           tmp_gkk_qpt(:,:,:) = elph_ds%gkk_qpt(:,jbeff,:,ik_this_proc,isppol,:)
         else if (elph_ds%gkqwrite == 1) then
           irec = (ik_this_proc-1)*elph_ds%k_phon%my_nkpt + iqpt
           if (iFSkpt == 1) then
             write (std_out,*) ' get_tau_k  read record ', irec
           end if
           read (elph_ds%unitgkq,REC=irec) tmp_gkk_qpt(:,:,iqpt_fullbz)
         end if

!FT to real space
         call ftgam(Ifc%wghatm,tmp_gkk_qpt,tmp_gkk_rpt,natom,nqpt,nrpt,1,coskr1,sinkr1)

!sum over irred q over k_phon, with corresponding weights
         do iFSqpt = 1, nkptirr
           iqpt_fullbz = elph_ds%k_phon%irredtoGS(iFSqpt)
           ikpt_kpq = FSfullpktofull(iFSkpt,iqpt_fullbz)

           imqpt_fullbz = mqtofull(iqpt_fullbz)
           ikpt_kmq = FSfullpktofull(iFSkpt,imqpt_fullbz)

!Do FT from real-space gamma grid to 1 kpt in k_phon%new_kptirr
           call ftgam(Ifc%wghatm,tmp_gkk_kpt,tmp_gkk_rpt,natom,1,nrpt,0,coskr2(iqpt_fullbz,:),sinkr2(iqpt_fullbz,:))
!tmp_gkk_kpt(:,:)=tmp_gkk_qpt(:,:,iFSqpt)

!if ep_scalprod==0 we have to dot in the displacement vectors here
           if (elph_ds%ep_scalprod==0) then

             call phdispl_cart2red(natom,Cryst%gprimd,displ(:,:,:,iFSqpt),displ_red)

             tmp_gkk_kpt2 = reshape (tmp_gkk_kpt(:,:), (/2,nbranch,nbranch/))
             call gam_mult_displ(nbranch, displ_red, tmp_gkk_kpt2, gkk_kpt)

             do jbranch=1,nbranch
               eigval(jbranch) = gkk_kpt(1, jbranch, jbranch)
               imeigval(jbranch) = gkk_kpt(2, jbranch, jbranch)

               if (abs(imeigval(jbranch)) > tol10) then
                 write (message,'(a,i0,a,es16.8)')" real values  branch = ",jbranch,' eigval = ',eigval(jbranch)
                 MSG_WARNING(message)
                 write (message,'(a,i0,a,es16.8)')" imaginary values  branch = ",jbranch,' imeigval = ',imeigval(jbranch)
                 MSG_WARNING(message)
               end if

             end do

!            if ep_scalprod==1 we have to diagonalize the matrix we interpolated.
           else if (elph_ds%ep_scalprod == 1) then

!            MJV NOTE : gam_now is being recast as a (3*natom)**2 matrix here
             call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, cone, tmp_gkk_kpt, 3*natom,&
&             pheigvec(:,iFSqpt), 3*natom, czero, tmp_gkk_kpt2, 3*natom)

             call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, cone, pheigvec(:,iFSqpt), 3*natom,&
&             tmp_gkk_kpt2, 3*natom, czero, gkk_kpt, 3*natom)

             diagerr = zero
             do ibranch=1,nbranch
               eigval(ibranch) = gkk_kpt(1,ibranch,ibranch)
               do jbranch=1,ibranch-1
                 diagerr = diagerr + abs(gkk_kpt(1,jbranch,ibranch))
               end do
               do jbranch=ibranch+1,nbranch
                 diagerr = diagerr + abs(gkk_kpt(1,jbranch,ibranch))
               end do
             end do

             if (diagerr > tol12) then
               write(message,'(a,es15.8)') 'get_tau_k: residual in diagonalization of gamma with phon eigenvectors: ', diagerr
               MSG_WARNING(message)
             end if

           else
             write (message,'(a,i0)')' Wrong value for ep_scalprod = ',elph_ds%ep_scalprod
             MSG_BUG(message)
           end if ! end ep_scalprod if

!For k'=k-q
!Do FT from real-space gamma grid to 1 kpt in k_phon%new_kptirr
           call ftgam(Ifc%wghatm,tmp_gkk_kpt,tmp_gkk_rpt,natom,1,nrpt,0,coskr2(imqpt_fullbz,:),sinkr2(imqpt_fullbz,:))
!tmp_gkk_kpt(:,:)=tmp_gkk_qpt(:,:,iFSqpt)

!if ep_scalprod==0 we have to dot in the displacement vectors here
           if (elph_ds%ep_scalprod==0) then

             call phdispl_cart2red(natom,Cryst%gprimd,displ(:,:,:,iFSqpt),displ_red)

             tmp_gkk_kpt2 = reshape (tmp_gkk_kpt(:,:), (/2,nbranch,nbranch/))
             call gam_mult_displ(nbranch, displ_red, tmp_gkk_kpt2, gkk_kpt)

             do jbranch=1,nbranch
               eigval2(jbranch) = gkk_kpt(1, jbranch, jbranch)
               imeigval(jbranch) = gkk_kpt(2, jbranch, jbranch)

               if (abs(imeigval(jbranch)) > tol10) then
                 write (message,'(a,i0,a,es16.8)')" real values  branch = ",jbranch,' eigval = ',eigval2(jbranch)
                 MSG_WARNING(message)
                 write (message,'(a,i0,a,es16.8)')" imaginary values  branch = ",jbranch,' imeigval = ',imeigval(jbranch)
                 MSG_WARNING(message)
               end if

             end do

!            if ep_scalprod==1 we have to diagonalize the matrix we interpolated.
           else if (elph_ds%ep_scalprod == 1) then

!            MJV NOTE : gam_now is being recast as a (3*natom)**2 matrix here
             call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, cone, tmp_gkk_kpt, 3*natom,&
&             pheigvec(:,iFSqpt), 3*natom, czero, tmp_gkk_kpt2, 3*natom)

             call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, cone, pheigvec(:,iFSqpt), 3*natom,&
&             tmp_gkk_kpt2, 3*natom, czero, gkk_kpt, 3*natom)

             diagerr = zero
             do ibranch=1,nbranch
               eigval2(ibranch) = gkk_kpt(1,ibranch,ibranch)
               do jbranch=1,ibranch-1
                 diagerr = diagerr + abs(gkk_kpt(1,jbranch,ibranch))
               end do
               do jbranch=ibranch+1,nbranch
                 diagerr = diagerr + abs(gkk_kpt(1,jbranch,ibranch))
               end do
             end do

             if (diagerr > tol12) then
               write(message,'(a,es15.8)') 'get_tau_k: residual in diagonalization of gamma with phon eigenvectors: ', diagerr
               MSG_WARNING(message)
             end if

           else
             write (message,'(a,i0)')' Wrong value for ep_scalprod = ',elph_ds%ep_scalprod
             MSG_BUG(message)
           end if ! end ep_scalprod if

           tmp2_wtk(:) = tmp_wtk(jpband,ikpt_kpq,isppol,:)
           yp1 = (tmp2_wtk(2)-tmp2_wtk(1))/nspline/deltaene
           ypn = (tmp2_wtk(nene_all)-tmp2_wtk(nene_all-1))/nspline/deltaene
           call spline(ene_pt,tmp2_wtk,nene_all,yp1,ypn,ff2)
           call splint(nene_all,ene_pt,tmp2_wtk,ff2,nene_all*nspline,ene_ptfine,tmp_wtk1)

           tmp2_wtk(:) = tmp_wtk(jpband,ikpt_kmq,isppol,:)
           yp1 = (tmp2_wtk(2)-tmp2_wtk(1))/nspline/deltaene
           ypn = (tmp2_wtk(nene_all)-tmp2_wtk(nene_all-1))/nspline/deltaene
           call spline(ene_pt,tmp2_wtk,nene_all,yp1,ypn,ff2)
           call splint(nene_all,ene_pt,tmp2_wtk,ff2,nene_all*nspline,ene_ptfine,tmp_wtk2)

           tmp2_wtq(:,:) = wtq(:,iFSqpt,:)
           do iene=1,nene
             e_k = eigenGS(elph_ds%minFSband+jband-1,iFSkpt,isppol)
             ene = e_k - omega_max + (iene-1)*deltaene
             if (ene<enemin .or. ene>enemax) cycle
             iene_fine = NINT((ene-enemin+deltaene)/deltaene)
             tmp_wtkpq = tmp_wtk1(iene_fine) * elph_ds%k_phon%wtkirr(iFSqpt)
             tmp_wtkmq = tmp_wtk2(iene_fine) * elph_ds%k_phon%wtkirr(iFSqpt)

             if (tmp_wtkpq+tmp_wtkmq < tol_wtk ) then
               nskip = nskip +1
               cycle
             end if

             do ibranch = 1, nbranch
               if (abs(phfrq(ibranch,iFSqpt)) < tol7) cycle

               if (ene > e_k) then
                 omega = ene - e_k
                 if (abs(omega) < tol7 .or. abs(omega) > omega_max) cycle
                 iomega = NINT((omega-omega_min+domega)/domega)

                 a2f_2d(iene) = a2f_2d(iene) +&
&                 eigval(ibranch)/phfrq(ibranch,iFSqpt)*&
&                 tmp_wtkpq * tmp2_wtq(ibranch,iomega)
               end if

               if (ene < e_k) then
                 omega = e_k - ene
                 if (abs(omega) < tol7 .or. abs(omega) > omega_max) cycle
                 iomega = NINT((omega-omega_min+domega)/domega)

                 a2f_2d2(iene) = a2f_2d2(iene) +&
&                 eigval(ibranch)/phfrq(ibranch,iFSqpt)*&
&                 tmp_wtkmq * tmp2_wtq(ibranch,iomega)
               end if

             end do ! ibranch 3
           end do ! nene  800
         end do ! kptirr 216
       end do ! j' band 3
!      print *, ' skipped ',  nskip, ' energy points out of ', nene*nband*nkptirr

! get inv_tau_k
       do itemp=1,ntemper  ! runs over termperature in K
         Temp=elph_ds%tempermin+elph_ds%temperinc*dble(itemp)
         do iene=1,nene
           e_k = eigenGS(elph_ds%minFSband+jband-1,iFSkpt,isppol)
           ene = e_k - omega_max + (iene-1)*deltaene
           if (ene<enemin .or. ene>enemax) cycle

           xx=(ene-fermie(itemp))/(kb_HaK*Temp)
           occ_e=1.0_dp/(exp(xx)+1.0_dp)
           if (ene > e_k .and. (ene-e_k) .le. omega_max) then
             omega = ene - e_k
             if (abs(omega) < tol7) cycle
             xx = omega/(kb_HaK*Temp)
             occ_omega=1.0_dp/(exp(xx)-1.0_dp)

             therm_factor = occ_e + occ_omega

             inv_tau_k(itemp,isppol,iFSkpt,jband) = inv_tau_k(itemp,isppol,iFSkpt,jband) +&
             a2f_2d(iene)*therm_factor*deltaene
           end if
           if (ene < e_k .and. (e_k-ene) .le. omega_max) then
             omega = e_k - ene
             if (abs(omega) < tol7) cycle
             xx = omega/(kb_HaK*Temp)
             occ_omega=1.0_dp/(exp(xx)-1.0_dp)

             therm_factor = 1 - occ_e + occ_omega

             inv_tau_k(itemp,isppol,iFSkpt,jband) = inv_tau_k(itemp,isppol,iFSkpt,jband) +&
             a2f_2d2(iene)*therm_factor*deltaene
           end if

         end do ! nene
       end do ! Temp
!          write(*,*)'i am here 2 ', isppol,iFSkpt,jband
     end do ! jband
   end do ! kpt
 end do ! nsppol

!write (300+mpi_enreg%me,*) inv_tau_k
 call xmpi_sum (inv_tau_k, xmpi_world, ierr)

 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(displ)
 ABI_DEALLOCATE(pheigvec)
 ABI_DEALLOCATE(tmp2_wtk)
 ABI_DEALLOCATE(ff2)
 ABI_DEALLOCATE(ene_pt)
 ABI_DEALLOCATE(ene_ptfine)
 ABI_DEALLOCATE(tmp_wtk1)
 ABI_DEALLOCATE(tmp_wtk2)
 ABI_DEALLOCATE(tmp2_wtq)
 ABI_DEALLOCATE(wtq)
 ABI_DEALLOCATE(coskr1)
 ABI_DEALLOCATE(sinkr1)
 ABI_DEALLOCATE(coskr2)
 ABI_DEALLOCATE(sinkr2)
 ABI_DEALLOCATE(kpttokpt)
 ABI_DEALLOCATE(FSfullpktofull)
 ABI_DEALLOCATE(mqtofull)
 ABI_DEALLOCATE(tmp_gkk_qpt)
 ABI_DEALLOCATE(tmp_gkk_rpt)
 ABI_DEALLOCATE(tmp_gkk_kpt)
 ABI_DEALLOCATE(tmp_gkk_kpt2)
 ABI_DEALLOCATE(gkk_kpt)
 ABI_DEALLOCATE(a2f_2d)
 ABI_DEALLOCATE(a2f_2d2)

!output inv_tau_k and tau_k
 fname = trim(elph_ds%elph_base_name) // '_INVTAUK'
 if (open_file(fname,message,newunit=unit_invtau,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to relaxation time file
 write (unit_invtau,*) '# k-dep inverse of the relaxation time as a function of temperature.'
 write (unit_invtau,*) '# '
 write (unit_invtau,*) '# nkptirr= ', nkptirr, 'nband= ', nband
 write (unit_invtau,*) '# number of temperatures=  ', ntemper
 write (unit_invtau,*) '# tau [femtosecond^-1]     '

 fname = trim(elph_ds%elph_base_name) // '_TAUK'
 if (open_file(fname,message,newunit=unit_tau,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to relaxation time file
 write (unit_tau,*) '# k-dep relaxation time as a function of temperature.'
 write (unit_tau,*) '# '
 write (unit_tau,*) '# nkptirr= ', nkptirr, 'nband= ', nband
 write (unit_tau,*) '# number of temperatures=  ', ntemper
 write (unit_tau,*) '# tau [femtosecond]     '

 tau_k = zero
 do itemp=1,ntemper  ! runs over termperature in K
   Temp=elph_ds%tempermin+elph_ds%temperinc*dble(itemp)
   write(unit_invtau,'(a,f16.8)') '# Temperature = ', Temp
   write(unit_tau,'(a,f16.8)') '# Temperature = ', Temp
   do isppol=1,nsppol
     write(unit_invtau,'(a,i6)') '# For isppol = ', isppol
     write(unit_tau,'(a,i6)') '# For isppol = ', isppol
     do iFSkpt = 1,nkpt
!FIXME: check when tau_k is too small, whether there should be a phonon
!scattering or not, and should tau_k be zero or not.
       do jband = 1,nband
         if (abs(inv_tau_k(itemp,isppol,iFSkpt,jband)) < tol9) then
           inv_tau_k(itemp,isppol,iFSkpt,jband) = zero
           tau_k(itemp,isppol,iFSkpt,jband) = zero
         else
!no need to *nkpt due to wtkirr, as we need /nkpt for the sum
!no need to *two_pi due to the missing prefactor in gkk (see mka2f_tr_lova)
           inv_tau_k(itemp,isppol,iFSkpt,jband) = inv_tau_k(itemp,isppol,iFSkpt,jband)*elph_ds%occ_factor
           tau_k(itemp,isppol,iFSkpt,jband) = one/inv_tau_k(itemp,isppol,iFSkpt,jband)
         end if
       end do ! nband
       write(unit_invtau,'(a,i8,a,3f12.6)') '# kpt# ', iFSkpt, '   kpt=', elph_ds%k_phon%kptirr(:,iFSkpt)
       write(unit_invtau,'(100D16.8)') (inv_tau_k(itemp,isppol,iFSkpt,iband)*femto/chu_tau,iband=1,nband)
       write(unit_tau,'(a,i8,a,3f12.6)') '# kpt# ', iFSkpt, '   kpt=', elph_ds%k_phon%kptirr(:,iFSkpt)
       write(unit_tau,'(100D16.8)') (tau_k(itemp,isppol,iFSkpt,iband)*chu_tau/femto,iband=1,nband)
     end do ! nkptirr
     write(unit_invtau,*) ' '
     write(unit_tau,*) ' '
   end do ! nsppol
   write(unit_invtau,*) ' '
   write(unit_invtau,*) ' '
   write(unit_tau,*) ' '
   write(unit_tau,*) ' '
 end do ! ntemper

! Only use the irred k for eigenGS and tau_k
 ABI_ALLOCATE(tmp_eigenGS,(elph_ds%nband,elph_ds%k_phon%new_nkptirr,elph_ds%nsppol))

 do ikpt_irr = 1, new_nkptirr
   tmp_eigenGS(:,ikpt_irr,:) = eigenGS(:,elph_ds%k_phon%new_irredtoGS(ikpt_irr),:)
   tmp_tau_k(:,:,ikpt_irr,:) = tau_k(:,:,elph_ds%k_phon%new_irredtoGS(ikpt_irr),:)*chu_tau
 end do

!BoltzTraP output files in SIESTA format
 if (elph_ds%prtbltztrp == 1) then
   call ebands_prtbltztrp_tau_out (tmp_eigenGS(elph_ds%minFSband:elph_ds%maxFSband,:,:),&
&   elph_ds%tempermin,elph_ds%temperinc,ntemper,fermie, &
&   elph_ds%elph_base_name,elph_ds%k_phon%new_kptirr,nband,elph_ds%nelect,new_nkptirr, &
&   elph_ds%nspinor,nsppol,Cryst%nsym,Cryst%rprimd,Cryst%symrel,tmp_tau_k)
 end if !prtbltztrp
 ABI_DEALLOCATE(tmp_eigenGS)
 ABI_DEALLOCATE(tmp_tau_k)

!Get the energy dependence of tau.
!Eq. (6) in  Restrepo et al. Appl. Phys. Lett. 94, 212103 (2009) [[cite:Restrepo2009]]

 fname = trim(elph_ds%elph_base_name) // '_TAUE'
 if (open_file(fname,message,newunit=unit_taue,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to relaxation time file
 write (unit_taue,*) '# Energy-dep relaxation time as a function of temperature.'
 write (unit_taue,*) '# '
 write (unit_taue,*) '# number of temperatures=  ', ntemper
 write (unit_taue,*) '# ene[Ha] tau [femtosecond] DOS[au]    '

 fname = trim(elph_ds%elph_base_name) // '_MFP'
 if (open_file(fname,message,newunit=unit_mfp,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

 write (unit_mfp,*) '# Energy-dep mean free path as a function of temperature.'
 write (unit_mfp,*) '# '
 write (unit_mfp,*) '# number of temperatures=  ', ntemper
 write (unit_mfp,*) '# ene[Ha] mfp [femtometer]   '

 do itemp=1,ntemper  ! runs over termperature in K
   Temp=elph_ds%tempermin+elph_ds%temperinc*dble(itemp)
   write(unit_taue,'(a,f16.8)') '# Temperature = ', Temp
   do isppol = 1, nsppol
     write(unit_taue,*) '# Tau_e for isppol = ',isppol
     do iene = 1, nene_all
       rate_e = zero
       do iFSkpt = 1, nkpt
         do jband = 1, nband
           rate_e = rate_e + inv_tau_k(itemp,isppol,iFSkpt,jband)* &
&           tmp_wtk(jband,iFSkpt,isppol,iene)
         end do ! jband
       end do ! kpt
       if (dabs(dos_e(isppol,iene)) < tol7) then
         rate_e = zero
       else
         rate_e = rate_e/nkpt/dos_e(isppol,iene)
       end if
       write(unit_taue,"(3D16.8)") enemin+(iene-1)*deltaene*nspline, rate_e*femto/chu_tau, dos_e(isppol,iene)
     end do ! number of energies
     write(unit_taue,*) ' '
   end do ! nsppol
   write(unit_taue,*) ' '
 end do ! ntemperature

! calculate and output mean free path
 do itemp=1,ntemper  ! runs over termperature in K
   Temp=elph_ds%tempermin+elph_ds%temperinc*dble(itemp)
   write(unit_mfp,'(a,f16.8)') '# Temperature = ', Temp
   do isppol = 1, nsppol
     do icomp = 1, 3
       write(unit_mfp,*) '# Mean free path for isppol, icomp= ',isppol,icomp
       do iene = 1, nene_all
         mfp_e = zero
         do iFSkpt = 1, nkptirr
           do jband = 1, nband
             mfp_e = mfp_e + tau_k(itemp,isppol,iFSkpt,jband)* &
&             elph_tr_ds%el_veloc(iFSkpt,elph_ds%minFSband+jband-1,icomp,isppol)* &
&             tmp_wtk(jband,iFSkpt,isppol,iene)
!&                          elph_ds%k_phon%new_wtkirr(iFSqpt)
           end do ! jband
         end do ! kpt
         if (dabs(dos_e(isppol,iene)) < tol7) then
           mfp_e = zero
         else
           mfp_e = mfp_e/nkptirr/dos_e(isppol,iene)
         end if
         write(unit_mfp,"(2D16.8)") enemin+(iene-1)*deltaene*nspline, mfp_e*chu_mfp/femto
       end do ! number of energies
       write(unit_mfp,*) ' '
     end do ! icomp
     write(unit_mfp,*) ' '
   end do ! nsppol
   write(unit_mfp,*) ' '
 end do ! ntemperature

 ABI_ALLOCATE(cond_e ,(ntemper,nsppol,nene_all,9))

!get cond_e
 cond_e = zero
 do itemp=1,ntemper  ! runs over termperature in K
   do isppol = 1, nsppol
     do iene = 1, nene_all
!       do iFSkpt =1,nkpt
       do ik_this_proc =1,elph_ds%k_phon%my_nkpt
         iFSkpt = elph_ds%k_phon%my_ikpt(ik_this_proc)
         do jband = 1, nband
           do icomp = 1, 3
             do jcomp = 1, 3
               itensor = (icomp-1)*3+jcomp
               cond_e(itemp,isppol,iene,itensor) = cond_e(itemp,isppol,iene,itensor) + &
&               tau_k(itemp,isppol,iFSkpt,jband)* &
&               elph_tr_ds%el_veloc(iFSkpt,elph_ds%minFSband+jband-1,icomp,isppol)* &
&               elph_tr_ds%el_veloc(iFSkpt,elph_ds%minFSband+jband-1,jcomp,isppol)* &
&               tmp_wtk(jband,iFSkpt,isppol,iene)
             end do
           end do
         end do ! jband
       end do ! kpt
     end do ! number of energies
   end do ! nsppol
 end do ! ntemperature

 ! MG FIXME: Why xmpi_world, besides only master should perform IO in the section below.
 call xmpi_sum (cond_e, xmpi_world, ierr)

 cond_e = cond_e/nkpt

!get transport coefficients

 fname = trim(elph_ds%elph_base_name) // '_COND'
 if (open_file(fname,message,newunit=unit_cond,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to conductivity file
 write (unit_cond,*) '#  Conductivity as a function of temperature.'
 write (unit_cond,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_cond,*) '#  '
 write (unit_cond,*) '#  Columns are: '
 write (unit_cond,*) '#  temperature[K]   cond[au]   cond [SI]    '
 write (unit_cond,*) '#  '

 fname = trim(elph_ds%elph_base_name) // '_CTH'
 if (open_file(fname,message,newunit=unit_therm,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to thermal conductivity file
 write (unit_therm,'(a)') '# Thermal conductivity as a function of temperature.'
 write (unit_therm,'(a)') '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_therm,'(a)') '#  '
 write (unit_therm,'(a)') '#  Columns are: '
 write (unit_therm,'(a)') '#  temperature[K]   thermal cond [au]   thermal cond [SI]'
 write (unit_therm,'(a)') '#  '

 fname = trim(elph_ds%elph_base_name) // '_SBK'
 if (open_file(fname,message,newunit=unit_sbk,status='unknown') /=0) then
   MSG_ERROR(message)
 end if

!print header to relaxation time file
 write (unit_sbk,*) '# Seebeck Coefficint as a function of temperature.'
 write (unit_sbk,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_sbk,*) '#  '
 write (unit_sbk,*) '#  Columns are: '
 write (unit_sbk,*) '#  temperature[K]   S [au]   S [SI]     '
 write (unit_sbk,*) '#  '

 ABI_ALLOCATE(cond ,(ntemper,nsppol,3,3))
 ABI_ALLOCATE(cth ,(ntemper,nsppol,3,3))
 ABI_ALLOCATE(sbk ,(ntemper,nsppol,3,3))
 ABI_ALLOCATE(seebeck ,(ntemper,nsppol,3,3))

 cond = zero
 cth = zero
 sbk = zero
 seebeck = zero
 do isppol=1,nsppol
   do icomp=1, 3
     do jcomp=1, 3
       itensor=(icomp-1)*3+jcomp
       do itemp=1,ntemper
         Temp=elph_ds%tempermin + elph_ds%temperinc*dble(itemp)
         do iene = 1, nene_all
           factor = (enemin+(iene-1)*deltaene*nspline - fermie(itemp))/(kb_HaK*Temp)
           if (factor < -40.0d0) then
             dfermide = zero
           else if (factor > 40.0d0) then
             dfermide = zero
           else
             dfermide = EXP(factor)/(kb_HaK*Temp*(EXP(factor)+one)**2.0d0)
           end if
           cond(itemp,isppol,icomp,jcomp) = cond(itemp,isppol,icomp,jcomp) + &
&           cond_e(itemp,isppol,iene,itensor)*dfermide*deltaene*nspline
           cth(itemp,isppol,icomp,jcomp) = cth(itemp,isppol,icomp,jcomp) + cond_e(itemp,isppol,iene,itensor)* &
&           (enemin+(iene-1)*deltaene*nspline - fermie(itemp))**2.0d0*dfermide*deltaene*nspline
           sbk(itemp,isppol,icomp,jcomp) = sbk(itemp,isppol,icomp,jcomp) + cond_e(itemp,isppol,iene,itensor)* &
&           (enemin+(iene-1)*deltaene*nspline - fermie(itemp))*dfermide*deltaene*nspline
         end do
       end do ! temperature
     end do ! jcomp
   end do ! icomp
 end do !end isppol

 do isppol=1,nsppol
   do itemp=1,ntemper
     cond_inv(:,:)=cond(itemp,isppol,:,:)
     call matrginv(cond_inv,3,3)
     call DGEMM('N','N',3,3,3,one,sbk(itemp,isppol,:,:),3,cond_inv,&
&     3,zero,seebeck(itemp,isppol,:,:),3)
   end do
 end do

 do isppol=1,nsppol
   do icomp=1, 3
     do jcomp=1, 3
       itensor=(icomp-1)*3+jcomp
       write(unit_cond,*) '# Conductivity for isppol, itrten= ',isppol,itensor
       write(unit_therm,*) '# Thermal conductivity for isppol, itrten= ',isppol,itensor
       write(unit_sbk,*) '# Seebeck coefficient for isppol, itrten= ',isppol,itensor
       do itemp=1,ntemper
         Temp=elph_ds%tempermin + elph_ds%temperinc*dble(itemp)

         seebeck(itemp,isppol,icomp,jcomp) = -1.0d0*seebeck(itemp,isppol,icomp,jcomp)/(kb_HaK*Temp)
         cond(itemp,isppol,icomp,jcomp) = cond(itemp,isppol,icomp,jcomp)/cryst%ucvol
         cth(itemp,isppol,icomp,jcomp) = cth(itemp,isppol,icomp,jcomp)/(kb_HaK*Temp)/cryst%ucvol
         write(unit_cond,'(3D20.10)')Temp,cond(itemp,isppol,icomp,jcomp),cond(itemp,isppol,icomp,jcomp)*chu_cond
         write(unit_therm,'(3D20.10)')Temp,cth(itemp,isppol,icomp,jcomp),cth(itemp,isppol,icomp,jcomp)*chu_cth
         write(unit_sbk,'(3D20.10)')Temp,seebeck(itemp,isppol,icomp,jcomp),seebeck(itemp,isppol,icomp,jcomp)*chu_sbk
       end do ! temperature
       write(unit_cond,*)
       write(unit_therm,*)
       write(unit_sbk,*)
     end do ! jcomp
   end do ! icomp
 end do !end isppol


 ABI_DEALLOCATE(inv_tau_k)
 ABI_DEALLOCATE(tau_k)
 ABI_DEALLOCATE(tmp_wtk)
 ABI_DEALLOCATE(dos_e)
 ABI_DEALLOCATE(cond_e)
 ABI_DEALLOCATE(fermie)
 ABI_DEALLOCATE(cond)
 ABI_DEALLOCATE(sbk)
 ABI_DEALLOCATE(cth)
 ABI_DEALLOCATE(seebeck)

 close (unit=unit_tau)
 close (unit=unit_taue)
 close (unit=unit_mfp)
 close (unit=unit_invtau)
 close (unit=unit_cond)
 close (unit=unit_therm)
 close (unit=unit_sbk)

end subroutine get_tau_k
!!***

end module m_a2ftr
!!***
