!{\src2tex{textfont=tt}}
!!****f* ABINIT/mka2f
!!
!! NAME
!! mka2f
!!
!! FUNCTION
!!  calculate the FS averaged alpha^2F function
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! Cryst<crystal_t>=data type gathering info on the crystalline structure.
!! Ifc<ifc_type>=Object containing the interatomic force constants.
!!  elph_ds
!!    elph_ds%gkk2 = gkk2 matrix elements on full FS grid for each phonon mode
!!    elph_ds%nbranch = number of phonon branches = 3*natom
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%k_phon%nkpt = number of kpts included in the FS integration
!!    elph_ds%k_phon%kpt = coordinates of all FS kpoints
!!    elph_ds%k_phon%wtk = integration weights on the FS
!!    elph_ds%n0 = DOS at the Fermi level calculated from the k_phon integration weights (event. 2 spin pol)
!!  mustar = coulomb pseudopotential parameter
!!  natom = number of atoms
!!
!! OUTPUT
!!  a2f_1d = 1D alpha
!!  dos_phon = density of states for phonons
!!  elph_ds
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      d2c_wtq,ep_ph_weights,ftgam,ftgam_init,gam_mult_displ,ifc_fourq
!!      phdispl_cart2red,simpson_int,wrtout,zgemm
!!
!! NOTES
!!   copied from ftiaf9.f
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mka2f(Cryst,ifc,a2f_1d,dos_phon,elph_ds,kptrlatt,mustar)

 use defs_basis
 use defs_elphon
 use m_errors
 use m_profiling_abi

 use m_io_tools,        only : open_file
 use m_numeric_tools,   only : simpson_int, simpson
 use m_geometry,        only : phdispl_cart2red
 use m_crystal,         only : crystal_t
 use m_ifc,             only : ifc_type, ifc_fourq
 use m_dynmat,          only : ftgam_init, ftgam

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mka2f'
 use interfaces_14_hidewrite
 use interfaces_77_ddb, except_this_one => mka2f
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: mustar
 type(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: Cryst
 type(elph_type),target,intent(inout) :: elph_ds
!arrays
 integer, intent(in) :: kptrlatt(3,3)
 real(dp),intent(out) :: a2f_1d(elph_ds%na2f),dos_phon(elph_ds%na2f)

!Local variables -------------------------
!scalars
 integer :: natom,iFSqpt,ibranch,iomega,nbranch,na2f,nsppol,nkpt,nrpt
 integer :: isppol,jbranch,unit_a2f,unit_phdos,ep_scalprod
 real(dp) :: a2fprefactor,avgelphg,avglambda,avgomlog,diagerr
 real(dp) :: lambda_2,lambda_3,lambda_4,lambda_5
 real(dp) :: spinfact
 real(dp) :: lambda_iso(elph_ds%nsppol)
 real(dp) :: lqn,omega
 real(dp) :: omegalog(elph_ds%nsppol)
 real(dp) :: omlog_qn
 real(dp) :: tc_macmill,a2fsmear,domega,omega_min,omega_max
 real(dp) :: gaussval, gaussprefactor, gaussfactor, gaussmaxval, xx
 character(len=500) :: msg
 character(len=fnlen) :: fname,base_name
!arrays
 real(dp) :: displ_cart(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: eigval(elph_ds%nbranch)
 real(dp) :: gam_now(2,elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: imeigval(elph_ds%nbranch)
! real(dp) :: pheigvec(2*elph_ds%nbranch*elph_ds%nbranch),phfrq(elph_ds%nbranch)
 real(dp) :: tmp_a2f(elph_ds%na2f)
 real(dp) :: tmp_gam1(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_gam2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_phondos(elph_ds%na2f),n0(elph_ds%nsppol)
 real(dp),pointer :: kpt(:,:)
 real(dp),allocatable :: phfrq(:,:)
 real(dp),allocatable :: pheigvec(:,:)
 real(dp),allocatable :: tmp_wtq(:,:,:)
 real(dp),allocatable :: a2f1mom(:),a2f2mom(:),a2f3mom(:),a2f4mom(:)
 real(dp),allocatable :: a2f_1mom(:),a2flogmom(:)
 real(dp),allocatable :: a2flogmom_int(:)
 real(dp),allocatable :: coskr(:,:)
 real(dp),allocatable :: sinkr(:,:)

! *********************************************************************
!calculate a2f for frequencies between 0 and elph_ds%omega_max

 DBG_ENTER("COLL")

!might need kptrlatt for finer interpolation later
 ABI_UNUSED(kptrlatt(1,1))

 ! nrpt = number of real-space points for FT interpolation
 nrpt = Ifc%nrpt
 natom = Cryst%natom

 nbranch   =  elph_ds%nbranch
 na2f      =  elph_ds%na2f
 nsppol    =  elph_ds%nsppol
 base_name =  elph_ds%elph_base_name
 a2fsmear  =  elph_ds%a2fsmear
 nkpt      =  elph_ds%k_phon%nkpt
 kpt       => elph_ds%k_phon%kpt

 ep_scalprod = elph_ds%ep_scalprod
 n0        = elph_ds%n0

!spinfact should be 1 for a normal non sppol calculation without spinorbit
!for spinors it should also be 1 as bands are twice as numerous but n0 has been divided by 2
!for sppol 2 it should be 0.5 as we have 2 spin channels to sum
 spinfact = one/elph_ds%nsppol !/elph_ds%nspinor

!maximum value of frequency (a grid has to be chosen for the representation of alpha^2 F)
!WARNING! supposes this value has been set in mkelph_linwid.
 domega = (elph_ds%omega_max-elph_ds%omega_min)/(na2f-one)
 elph_ds%domega  = domega  ! MG Why do we need to store domega in elph_ds?
 omega_min       = elph_ds%omega_min
 omega_max       = elph_ds%omega_max

 gaussprefactor = sqrt(piinv) / a2fsmear
 gaussfactor = one / a2fsmear
 gaussmaxval = sqrt(-log(1.d-100))

 ! only open the file for the first sppol
 fname = trim(base_name) // '_A2F'
 if (open_file(fname,msg,newunit=unit_a2f,status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if

 !write (std_out,*) ' a2f function integrated over the FS'

!output the a2f_1d header
 write (unit_a2f,'(a)')                 '#'
 write (unit_a2f,'(a)')                 '# ABINIT package : a2f file'
 write (unit_a2f,'(a)')                 '#'
 write (unit_a2f,'(a)')                 '# a2f function integrated over the FS. omega in a.u.'
 write (unit_a2f,'(a,I10)')             '#  number of kpoints integrated over : ',nkpt
 write (unit_a2f,'(a,I10)')             '#  number of energy points : ',na2f
 write (unit_a2f,'(a,E16.6,a,E16.6,a)') '#  between omega_min = ',omega_min,' Ha and omega_max = ',omega_max,' Ha'
 write (unit_a2f,'(a,E16.6)')           '#  and the smearing width for gaussians is ',a2fsmear

 ! Open file for PH DOS
 fname = trim(base_name) // '_PDS'
 if (open_file(fname,msg,newunit=unit_phdos,status="replace") /= 0) then
   MSG_ERROR(msg)
 end if

 ! output the phonon DOS header
 write (unit_phdos,'(a)')                '#'
 write (unit_phdos,'(a)')                '# ABINIT package : phonon DOS file'
 write (unit_phdos,'(a)')                '#'
 write (unit_phdos,'(a)')                '# Phonon DOS integrated over the FS. omega in a.u. EXPERIMENTAL!!!'
 write (unit_phdos,'(a,I10)')            '# number of kpoints integrated over : ',nkpt
 write (unit_phdos,'(a,I10)')            '# number of energy points : ',na2f
 write (unit_phdos,'(a,E16.6,a,E16.6,a)')'# between omega_min = ',omega_min,' Ha and omega_max = ',omega_max,' Ha'
 write (unit_phdos,'(a,i4,a,E16.6)')     '# The DOS at Fermi level for spin ', 1, ' is ', n0(1)
 if (nsppol==2) then
   write (unit_phdos,'(a,i4,a,E16.6)')   '# The DOS at Fermi level for spin ', 2, ' is ', n0(2)
 end if
 write (unit_phdos,'(a,E16.6)')          '# and the smearing width for gaussians is ',a2fsmear
 write (unit_phdos,'(a)') '#'

!Get the integration weights, using tetrahedron method or gaussian
 ABI_ALLOCATE(tmp_wtq,(nbranch,elph_ds%k_fine%nkpt,na2f+1))
 ABI_ALLOCATE(elph_ds%k_fine%wtq,(nbranch,elph_ds%k_fine%nkpt,na2f))
 ABI_ALLOCATE(elph_ds%k_phon%wtq,(nbranch,nkpt,na2f))

 ABI_ALLOCATE(phfrq,(nbranch,elph_ds%k_fine%nkpt))
 ABI_ALLOCATE(pheigvec,(2*nbranch*nbranch,elph_ds%k_fine%nkpt))

 do iFSqpt=1,elph_ds%k_fine%nkpt
   call ifc_fourq(ifc,cryst,elph_ds%k_fine%kpt(:,iFSqpt),phfrq(:,iFSqpt),displ_cart,out_eigvec=pheigvec(:,iFSqpt))
 end do

 omega_min = omega_min - domega

 call ep_ph_weights(phfrq,elph_ds%a2fsmear,omega_min,omega_max,na2f+1,Cryst%gprimd,elph_ds%kptrlatt_fine, &
& elph_ds%nbranch,elph_ds%telphint,elph_ds%k_fine,tmp_wtq)
!call ep_ph_weights(phfrq,elph_ds%a2fsmear,omega_min,omega_max,na2f+1,Cryst%gprimd,elph_ds%kptrlatt_fine, &
!& elph_ds%nbranch,1,elph_ds%k_fine,tmp_wtq)
 omega_min = omega_min + domega

 do iomega = 1, na2f
   elph_ds%k_fine%wtq(:,:,iomega) = tmp_wtq(:,:,iomega+1)
 end do
 ABI_DEALLOCATE(tmp_wtq)

 if (elph_ds%use_k_fine == 1) then
   call d2c_wtq(elph_ds)
 end if

 ABI_ALLOCATE(coskr, (nkpt,nrpt))
 ABI_ALLOCATE(sinkr, (nkpt,nrpt))
 call ftgam_init(Ifc%gprim, nkpt, nrpt, kpt, Ifc%rpt, coskr, sinkr)

 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(pheigvec)

 do isppol=1,nsppol
   write (std_out,*) '##############################################'
   write (std_out,*) 'mka2f : Treating spin polarization ', isppol
   write (std_out,*) '##############################################'

!  Average of electron phonon coupling over the whole BZ
   avgelphg = zero
!  MG20060607 Do the same for lambda and omega_log
   avglambda = zero
   avgomlog = zero

   a2f_1d(:) = zero
   dos_phon(:) = zero

!  reduce the dimenstion from fine to phon for phfrq and pheigvec
   ABI_ALLOCATE(phfrq,(nbranch,elph_ds%k_phon%nkpt))
   ABI_ALLOCATE(pheigvec,(2*nbranch*nbranch,elph_ds%k_phon%nkpt))

!  loop over qpoint in full kpt grid (presumably dense)
!  MG TODO : This loop can be performed using the IBZ and appropriated weights.
   do iFSqpt=1,nkpt
!    
!    This reduced version of ftgkk supposes the kpoints have been integrated
!    in integrate_gamma. Do FT from real-space gamma grid to 1 qpt.

     if (elph_ds%ep_int_gkk == 1) then
       gam_now(:,:) = elph_ds%gamma_qpt(:,:,isppol,iFSqpt)
     else
       call ftgam(Ifc%wghatm,gam_now,elph_ds%gamma_rpt(:,:,isppol,:),natom,1,nrpt,0, &
&       coskr(iFSqpt,:), sinkr(iFSqpt,:))
     end if

     call ifc_fourq(ifc,cryst,kpt(:,iFSqpt),phfrq(:,iFSqpt),displ_cart,out_eigvec=pheigvec)

!    Diagonalize gamma matrix at qpoint (complex matrix).

!    if ep_scalprod==0 we have to dot in the displacement vectors here
     if (ep_scalprod==0) then

       call phdispl_cart2red(natom,Cryst%gprimd,displ_cart,displ_red)

       tmp_gam2 = reshape (gam_now, (/2,nbranch,nbranch/))
       call gam_mult_displ(nbranch, displ_red, tmp_gam2, tmp_gam1)

       do jbranch=1,nbranch
         eigval(jbranch) = tmp_gam1(1, jbranch, jbranch)
         imeigval(jbranch) = tmp_gam1(2, jbranch, jbranch)

         if (abs(imeigval(jbranch)) > tol8) then
           write (msg,'(a,i0,a,es16.8)')" imaginary values  branch = ",jbranch,' imeigval = ',imeigval(jbranch)
           MSG_WARNING(msg)
         end if

       end do

!      if ep_scalprod==1 we have to diagonalize the matrix we interpolated.
     else if (ep_scalprod == 1) then

!      MJV NOTE : gam_now is being recast as a (3*natom)**2 matrix here
       call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, cone, gam_now, 3*natom,&
&       pheigvec, 3*natom, czero, tmp_gam1, 3*natom)

       call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, cone, pheigvec, 3*natom,&
&       tmp_gam1, 3*natom, czero, tmp_gam2, 3*natom)

       diagerr = zero
       do ibranch=1,nbranch
         eigval(ibranch) = tmp_gam2(1,ibranch,ibranch)
         do jbranch=1,ibranch-1
           diagerr = diagerr + abs(tmp_gam2(1,jbranch,ibranch))
         end do
         do jbranch=ibranch+1,nbranch
           diagerr = diagerr + abs(tmp_gam2(1,jbranch,ibranch))
         end do
       end do

       if (diagerr > tol12) then
         write(msg,'(a,es15.8)') 'mka2f: residual in diagonalization of gamma with phon eigenvectors: ', diagerr
         MSG_WARNING(msg)
       end if

     else
       write (msg,'(a,i0)')' Wrong value for ep_scalprod = ',ep_scalprod
       MSG_BUG(msg)
     end if

!    MG20060603MG
!    there was a bug in the calculation of the phonon DOS
!    since frequencies with small e-ph interaction were skipped inside the loop
!    In this new version all the frequencies (both positive and negative) are taken into account.
!    IDEA: it could be useful to calculate the PH-dos and the a2f
!    using several smearing values to perform a convergence study
!    Now the case ep_scalprod=1 is treated in the right way although it is not default anymore
!    FIXME to be checked
!    ENDMG

!    Add all contributions from the phonon modes at this qpoint to a2f and the phonon dos.
     do ibranch=1,nbranch

!      if (abs(phfrq(ibranch,iFSqpt)) < tol10) then
       if (abs(phfrq(ibranch,iFSqpt)) < tol7) then
         a2fprefactor= zero
         lqn         = zero
         omlog_qn    = zero
       else
         a2fprefactor = eigval(ibranch)/(two_pi*abs(phfrq(ibranch,iFSqpt))*n0(isppol))
         lqn          = eigval(ibranch)/(pi*phfrq(ibranch,iFSqpt)**2*n0(isppol))
         omlog_qn     = lqn*log(abs(phfrq(ibranch,iFSqpt)))
       end if

!      Add contribution to average elphon coupling
!      MANY ISSUES WITH FINITE T SUMS. THIS IS DEFINITELY
!      NOT A CORRECT FORMULATION YET.

!      Added avglambda and avgomglog to calculate lamda and omega_log using the sum over the kpt-grid.
!      If the k-grid is dense enough, these values should be better than the corresponding quantities
!      evaluated through the integration over omega that depends on the a2fsmear

       avgelphg = avgelphg + eigval(ibranch)
       avglambda = avglambda + lqn
       avgomlog= avgomlog + omlog_qn
!      ENDMG

       omega = omega_min
       tmp_a2f(:) = zero
       tmp_phondos(:) = zero
       do iomega=1,na2f
         xx = (omega-phfrq(ibranch,iFSqpt))*gaussfactor
         omega = omega + domega
         if (abs(xx) > gaussmaxval) cycle

         gaussval = gaussprefactor*exp(-xx*xx)
         tmp_a2f(iomega) = tmp_a2f(iomega) + gaussval*a2fprefactor
         tmp_phondos(iomega) = tmp_phondos(iomega) + gaussval
       end do

!      tmp_a2f(:) = zero
!      tmp_phondos(:) = zero
!      do iomega=1,na2f
!      tmp_a2f(iomega) = tmp_a2f(iomega) + a2fprefactor*elph_ds%k_phon%wtq(ibranch,iFSqpt,iomega)
!      tmp_phondos(iomega) = tmp_phondos(iomega) + elph_ds%k_phon%wtq(ibranch,iFSqpt,iomega)
!      end do

       a2f_1d(:) = a2f_1d(:) + tmp_a2f(:)
       dos_phon(:) = dos_phon(:) + tmp_phondos(:)

     end do ! ibranch
   end do  ! iFSqpt do


!  second 1 / nkpt factor for the integration weights
   a2f_1d(:) = a2f_1d(:) / nkpt
   dos_phon(:) = dos_phon(:) / nkpt

!  MG
   avglambda = avglambda/nkpt
   avgomlog= avgomlog/nkpt
   avgomlog = exp (avgomlog/avglambda)
   write(std_out,*) ' from mka2f: for spin ', isppol
   write(std_out,*) ' w/o interpolation lambda = ',avglambda,' omega_log= ',avgomlog
!  ENDMG

   write (std_out,'(a,I4,a,E16.6)') '# The DOS at Fermi level for spin ',isppol,' is ',n0(isppol)

   write (unit_a2f,'(a,I4,a,E16.6)') '# The DOS at Fermi level for spin ',isppol,' is ',n0(isppol)
   write (unit_a2f,'(a)') '#'

   omega = omega_min
   do iomega=1,na2f
     write (unit_a2f,*) omega, a2f_1d(iomega)
     omega=omega + domega
   end do
   write (unit_a2f,*)
!  
!  output the phonon DOS, but only for the first sppol case
   if (isppol == 1) then
     omega = omega_min
     do iomega=1,na2f
       write (unit_phdos,*) omega, dos_phon(iomega)
       omega=omega + domega
     end do
   end if
!  
!  Do isotropic calculation of lambda and output lambda, Tc(MacMillan)
!  
   ABI_ALLOCATE(a2f_1mom,(na2f))
   ABI_ALLOCATE(a2f1mom,(na2f))
   ABI_ALLOCATE(a2f2mom,(na2f))
   ABI_ALLOCATE(a2f3mom,(na2f))
   ABI_ALLOCATE(a2f4mom,(na2f))

   a2f_1mom=zero
   a2f1mom=zero;  a2f2mom=zero
   a2f3mom=zero;  a2f4mom=zero
   
   omega = omega_min
   do iomega=1,na2f
     if (abs(omega) > tol10) then
       a2f_1mom(iomega) =    two*spinfact*a2f_1d(iomega)/abs(omega)   ! first inverse moment of alpha2F
       a2f1mom(iomega)  =    two*spinfact*a2f_1d(iomega)*abs(omega)   ! first positive moment of alpha2F
       a2f2mom(iomega)  =     a2f1mom(iomega)*abs(omega)  ! second positive moment of alpha2F
       a2f3mom(iomega)  =     a2f2mom(iomega)*abs(omega)  ! third positive moment of alpha2F
       a2f4mom(iomega)  =     a2f3mom(iomega)*abs(omega)  ! fourth positive moment of alpha2F
     end if
     omega=omega + domega
   end do
!  
!  From Allen PRL 59 1460
!  \lambda <\omega^n> = 2 \int_0^{\infty} d\omega [\alpha^2F / \omega] \omega^n
!  
   lambda_iso(isppol) = simpson(domega,a2f_1mom)
   lambda_2 = simpson(domega,a2f1mom)
   lambda_3 = simpson(domega,a2f2mom)
   lambda_4 = simpson(domega,a2f3mom)
   lambda_5 = simpson(domega,a2f4mom)

   ABI_DEALLOCATE(phfrq)
   ABI_DEALLOCATE(pheigvec)
   ABI_DEALLOCATE(a2f_1mom)
   ABI_DEALLOCATE(a2f1mom)
   ABI_DEALLOCATE(a2f2mom)
   ABI_DEALLOCATE(a2f3mom)
   ABI_DEALLOCATE(a2f4mom)

   write (std_out,*) 'mka2f: elphon coupling lambdas for spin = ', isppol
   write (std_out,*) 'mka2f: isotropic lambda', lambda_iso(isppol)
   write (std_out,*) 'mka2f: positive moments of alpha2F:'
   write (std_out,*) 'lambda <omega^2> = ', lambda_2
   write (std_out,*) 'lambda <omega^3> = ', lambda_3
   write (std_out,*) 'lambda <omega^4> = ', lambda_4
   write (std_out,*) 'lambda <omega^5> = ', lambda_5
!  
!  Get log moment of alpha^2F
   ABI_ALLOCATE(a2flogmom,(na2f))
   ABI_ALLOCATE(a2flogmom_int,(na2f))
   omega = omega_min
   a2flogmom(:) = zero
   do iomega=1,na2f
     if (abs(omega) > tol10) then
       a2flogmom(iomega) = a2f_1d(iomega)*log(abs(omega))/abs(omega)
     end if
     omega=omega + domega
   end do
   call simpson_int(na2f,domega,a2flogmom,a2flogmom_int)

!  NOTE: omegalog actually stores the log moment of a2F, which is the quantity to sum over spins, instead of 
!  exp(moment/lambda) which is an actual frequency
   omegalog(isppol) = two*spinfact*a2flogmom_int(na2f)

   ABI_DEALLOCATE(a2flogmom)
   ABI_DEALLOCATE(a2flogmom_int)
   
   if (nsppol > 1) then
     write (msg, '(3a)' ) ch10,&
&     ' Warning : some of the following quantities should be integrated over spin', ch10
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end if

   write (msg, '(3a)' ) ch10,&
&   ' Superconductivity : isotropic evaluation of parameters from electron-phonon coupling.',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   if (elph_ds%nsppol > 1) then
     write (msg, '(a,i6,a,es16.6)' )' mka2f: isotropic lambda for spin ', isppol, ' = ', lambda_iso(isppol)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end if

   write (msg, '(a,es16.6)' )' mka2f: lambda <omega^2> = ', lambda_2
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write (msg, '(a,es16.6)' )' mka2f: lambda <omega^3> = ', lambda_3
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write (msg, '(a,es16.6)' )' mka2f: lambda <omega^4> = ', lambda_4
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write (msg, '(a,es16.6)' )' mka2f: lambda <omega^5> = ', lambda_5
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   if (elph_ds%nsppol > 1) then
     write (msg, '(a,i6,a,es16.6,a,es16.6,a)' )' mka2f: omegalog for spin ', isppol, ' = ',&
&     exp(omegalog(isppol)/lambda_iso(isppol)), ' (Ha) ', exp(omegalog(isppol)/lambda_iso(isppol))/kb_HaK, ' (Kelvin) '
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end if

 end do ! isppol

!also print out spin-summed quantities
 lambda_2 = sum(lambda_iso(1:elph_ds%nsppol))
 write (msg, '(a,es16.6)' )' mka2f: isotropic lambda = ', lambda_2
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 omega = exp( sum(omegalog(1:elph_ds%nsppol))/lambda_2 )
 write (msg, '(a,es16.6,a,es16.6,a)' )' mka2f: omegalog  = ', omega, ' (Ha) ', omega/kb_HaK, ' (Kelvin) '
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 write (msg, '(a,es16.6)' )' mka2f: input mustar = ', mustar
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 tc_macmill = omega/1.2_dp * exp((-1.04_dp*(one+lambda_2)) / (lambda_2-mustar*(one+0.62_dp*lambda_2)))
 write ( msg, '(a,es16.6,a,es16.6,a)')'-mka2f: MacMillan Tc = ', tc_macmill, ' (Ha) ', tc_macmill/kb_HaK, ' (Kelvin) '
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 close(unit=unit_a2f)
 close(unit=unit_phdos)

 ABI_DEALLOCATE(elph_ds%k_fine%wtq)
 ABI_DEALLOCATE(elph_ds%k_phon%wtq)

 ABI_DEALLOCATE(coskr)
 ABI_DEALLOCATE(sinkr)

 DBG_EXIT("COLL")

end subroutine mka2f
!!***
