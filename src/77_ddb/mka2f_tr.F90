!{\src2tex{textfont=tt}}
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
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (JPC, MJV, BXU)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYINGS=
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mka2f_tr(crystal,ifc,elph_ds,ntemper,tempermin,temperinc,pair2red,elph_tr_ds)

 use defs_basis
 use defs_elphon
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_io_tools,        only : open_file
 use m_numeric_tools,   only : simpson_int
 use m_crystal,         only : crystal_t
 use m_ifc,             only : ifc_type, ifc_fourq
 use m_dynmat,          only : ftgam_init, ftgam

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mka2f_tr'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_77_ddb, except_this_one => mka2f_tr
!End of the abilint section

 implicit none

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
   call ifc_fourq(ifc,crystal,elph_ds%k_fine%kpt(:,iFSqpt),phfrq(:,iFSqpt),displ(:,:,:,iFSqpt),out_eigvec=pheigvec(:,iFSqpt))
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
   call ifc_fourq(ifc,crystal,elph_ds%k_phon%kpt(:,iFSqpt),phfrq(:,iFSqpt),displ(:,:,:,iFSqpt),out_eigvec=pheigvec(:,iFSqpt))
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

!Get thermal resistivity, based on eqn. (52) in Allen's PRB 17, 3725 (1978)
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
