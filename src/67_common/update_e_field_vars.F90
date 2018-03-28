!{\src2tex{textfont=tt}}
!!****f* ABINIT/update_e_field_vars
!! NAME
!! update_e_field_vars
!!
!! FUNCTION
!! This routine updates E field variables 
!!
!! COPYRIGHT
!! Copyright (C) 2003-2018 ABINIT  group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! atindx(natom)=index table for atoms, inverse of atindx (see gstate.f)
!! atindx1(natom)=index table for atoms (see gstate.f)
!! cg(2,mcg)=planewave coefficients of wavefunctions
!! dimcprj(usepaw*natom)=lmn_size for each atom
!! dtfil <type(datafiles_type)>=variables related to files
!! gmet(3,3)=metric in reciprocal space
!! gprimd(3,3)=reciprocal space dimensional primitive translations
!! idir = determines directions for derivatives computed in ctocprj (0 for all)
!! kg(3,mpw*mkmem)=reduced planewave coordinates
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mkmem=number of k points treated by this node.
!! mpw=maximum dimensioned size of npw
!! my_natom=number of atoms treated by current processor
!! natom=number of atoms in cell
!! nattyp(ntypat)=number of atoms of each type
!! ngfft(18)=contain all needed information about 3D FFT, see ~ABINIT/Infos/vargs.htm#ngfft
!! nkpt=number of k-points
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! ntypat=number of types of atoms in unit cell
!! pawrhoij(natom*usepaw) <type(pawrhoij_type)> atomic occupancies
!! pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!! pwind_alloc = first dimension of pwind
!! pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!! rmet(3,3)=metric in real space
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! scfcv_level= 0 if calling before scf loop, 1 if during
!! scfcv_quit=signals whether calling during scf quit (see scfcv.F90)
!! scfcv_step=istep value of loop counter from scfcv.F90
!! ucvol=unit cell volume in bohr**3.
!! unit_out= unit for output of the results (usually the .out file of ABINIT)
!!   The option unit_out = 0 is allowed. In this case, no information is written
!!   to the output file but only to the log file.
!! usepaw= 1: use paw framework. 0:do not use paw.
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for
!!     each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real
!!     spherical harmonics
!!
!! OUTPUT
!! efield_old_cart(3)=updating cartesian values of efield (used in berryopt
!!                    6,16,17)
!! pel_cg(3)=electronic polarization
!! pelev(3)=leading order PAW contribution in pel_cg (for reporting purposes
!!          only)
!! pion(3)=ionic part of polarization
!! ptot(3)=total polarization
!! red_efield2=updating efield used in berryopt 16,17
!! red_efield2_old=updating efield used in berryopt 16.17
!! red_ptot=updating efield used in berryopt 16.17
!!
!! SIDE EFFECTS
!! Input/Output
!! dtset <type(dataset_type)>=all input variables in this dataset
!! dtefield <type(efield_type)> = efield variables
!! hdr <type(hdr_type)>=the header of wf, den and pot files
!! mpi_enreg=information about MPI parallelization
!! ptot_cart(3)=total polarization in cartesian coordinates
!! xred(3,natom)=reduced atomic coordinates
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      berryphase_new,ctocprj,getph,pawcprj_alloc,pawcprj_free,prtefield
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine update_e_field_vars(atindx,atindx1,cg,dimcprj,dtefield,dtfil,dtset,&
&  efield_old_cart,gmet,gprimd,hdr,idir,kg,mcg,&
&  mkmem,mpi_enreg,mpw,my_natom,natom,nattyp,ngfft,nkpt,npwarr,ntypat,&
&  pawrhoij,pawtab,pel_cg,pelev,pion,psps,ptot,ptot_cart,pwind,&
&  pwind_alloc,pwnsfac,red_efield2,red_efield2_old,red_ptot,rmet,rprimd,&
&  scfcv_level,scfcv_quit,scfcv_step,ucvol,unit_out,&
&  usepaw,xred,ylm,ylmgr)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_errors
 use m_efield
 use m_profiling_abi

 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_kg,       only : getph

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'update_e_field_vars'
 use interfaces_14_hidewrite
 use interfaces_56_recipspace
 use interfaces_66_nonlocal
 use interfaces_67_common, except_this_one => update_e_field_vars
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: idir,mcg,mkmem,mpw,my_natom,natom,nkpt,ntypat
 integer, intent(in) :: pwind_alloc,scfcv_level,scfcv_quit,scfcv_step,unit_out,usepaw
 real(dp), intent(in) :: ucvol
 type(datafiles_type), intent(in) :: dtfil
 type(pseudopotential_type),intent(in) :: psps
 type(dataset_type), intent(inout) :: dtset
 type(efield_type), intent(inout) :: dtefield
 type(hdr_type), intent(inout) :: hdr
 type(MPI_type), intent(inout) :: mpi_enreg
!arrays
 integer, intent(in) :: atindx(natom),atindx1(natom),dimcprj(usepaw*natom)
 integer, intent(in) :: kg(3,mpw*mkmem),nattyp(ntypat)
 integer, intent(in) :: ngfft(18),npwarr(nkpt),pwind(pwind_alloc,2,3)
 real(dp), intent(in) :: cg(2,mcg),gmet(3,3),gprimd(3,3)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
 real(dp), intent(in) :: rmet(3,3),rprimd(3,3)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(inout) :: ptot_cart(3),xred(3,natom),efield_old_cart(3) !vz_i
 real(dp), intent(out) :: pel_cg(3),pelev(3),pion(3) !vz_i
 real(dp), intent(inout) :: red_efield2(3),red_efield2_old(3) !vz_i
 real(dp), intent(out) :: ptot(3),red_ptot(3) !vz_i
 type(pawrhoij_type), intent(in) :: pawrhoij(my_natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables -------------------------
!scalars
 character(len=500) :: message
 integer :: ctocprj_choice,iatom,ii,iorder_cprj,mcprj,my_nspinor,ncpgr
 integer :: optberry,usecprj
 logical :: calc_epaw3_force, calc_epaw3_stress, efield
!arrays
 real(dp) :: efield_test_cart(3),red_efield1(3)
 real(dp),allocatable :: ph1d(:,:)
 type(pawcprj_type),allocatable :: cprj(:,:)

! *************************************************************************

 efield = .false.

 if ( dtset%berryopt == 4 .or. &
& dtset%berryopt == 6 .or. &
& dtset%berryopt == 7 .or. &
& dtset%berryopt ==14 .or. &
& dtset%berryopt ==16 .or. &
& dtset%berryopt ==17 ) efield = .true.
 calc_epaw3_force = ( efield .and. dtset%optforces /= 0 .and. usepaw == 1 )
 calc_epaw3_stress = ( efield .and. dtset%optstress /= 0  .and. usepaw == 1 )

 usecprj=1; if (psps%usepaw==0)  usecprj = 0
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 mcprj=my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol

 ncpgr = 0
 ctocprj_choice = 1 ! no derivs
 if ( efield .and. psps%usepaw == 1) then
   ABI_DATATYPE_ALLOCATE(cprj,(dtset%natom,mcprj))
!  finite electric field may need gradients for forces, stress
   if (calc_epaw3_force .and. .not. calc_epaw3_stress) then
     ncpgr = 3; ctocprj_choice = 2 ! derivs w.r.t. position
   else if (.not. calc_epaw3_force .and. calc_epaw3_stress) then
     ncpgr = 6; ctocprj_choice = 3 ! derivs w.r.t strain
   else if (calc_epaw3_force .and. calc_epaw3_stress) then
     ncpgr = 9; ctocprj_choice = 23 ! derivs w.r.t. position and strain
   end if
   call pawcprj_alloc(cprj,ncpgr,dimcprj)
   iatom=0 ; iorder_cprj=1 ! retain ordering of input list
!  all arguments to ctocprj are defined already except ph1d, do that here
   ABI_ALLOCATE(ph1d,(2,3*(2*dtset%mgfft+1)*dtset%natom))
   call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)
   call ctocprj(atindx,cg,ctocprj_choice,cprj,gmet,gprimd,iatom,idir,iorder_cprj,&
&   dtset%istwfk,kg,dtset%kptns,mcg,mcprj,dtset%mgfft,dtset%mkmem,&
&   mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,nattyp,dtset%nband,&
&   dtset%natom,ngfft,dtset%nkpt,dtset%nloalg,npwarr,dtset%nspinor,&
&   dtset%nsppol,dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,&
&   dtset%typat,ucvol,dtfil%unpaw,xred,ylm,ylmgr)
   ABI_DEALLOCATE(ph1d)
 else 
   ABI_DATATYPE_ALLOCATE(cprj,(0,0))
 end if ! end update of cprj

 if ( efield ) then ! compute polarization and if necessary store cprj in efield
   optberry=1 
   pel_cg(:) = zero;pelev=zero
   call berryphase_new(atindx1,cg,cprj,dtefield,dtfil,dtset,psps,gprimd,hdr,psps%indlmn,kg,&
&   psps%lmnmax,dtset%mband,mcg,mcprj,dtset%mkmem,mpi_enreg,dtset%mpw,my_natom,&
&   dtset%natom,npwarr,dtset%nsppol,psps%ntypat,dtset%nkpt,optberry,pawrhoij,pawtab,&
&   pel_cg,pelev,pion,ptot,red_ptot,pwind,&
&   pwind_alloc,pwnsfac,rprimd,dtset%typat,ucvol,&
&   unit_out,usecprj,psps%usepaw,xred,psps%ziontypat)

   dtefield%red_ptot1(:)=red_ptot(:)

 end if ! end compute polarization and store cprj for efield

 if (efield .and. (scfcv_level == 0) ) then ! do this before scfcv loop

   efield_old_cart(:)=dtset%efield(:)   !!HONG
   
!  save this value in order to print the final value of real electric field, comparing with the desired red_fieldbar
   dtefield%efield2(:)=dtset%efield(:)

   if ( dtset%berryopt ==16 .or. dtset%berryopt ==17) then   !!HONG
     do ii=1,3
       red_efield2(ii)=zero
       red_efield2_old(ii)  =(ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,ii))
     end do
   end if
   
   if (dtset%berryopt == 14 .and. scfcv_quit /=1) then
!    ! Convert polarization to cartesian coords

     ptot_cart(:)=zero
     do ii = 1,3
       ptot_cart(ii)=rprimd(ii,1)*red_ptot(1) + rprimd(ii,2)*red_ptot(2) + &
&       rprimd(ii,3)*red_ptot(3)
     end do
     ptot_cart(:)=ptot_cart(:)/ucvol

     do ii=1,3
       dtefield%efield_dot(ii) = dot_product(dtset%efield(:),rprimd(:,ii))
     end do

!    !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
     write(message,'(a,a)')   ch10, 'scfcv: Constant reduced ebar-field:'

     call wrtout(std_out,message,'COLL')
     call prtefield(dtset,dtefield,std_out,rprimd)

     if(dtset%prtvol>=10)then
       call wrtout(ab_out,message,'COLL')
       call prtefield(dtset,dtefield,ab_out,rprimd)
     end if

!    updating E field
     do ii =1,3   ! desired E field
       efield_test_cart(ii)=gprimd(ii,1)*dtset%red_efieldbar(1) + &
&       gprimd(ii,2)*dtset%red_efieldbar(2)+gprimd(ii,3)*dtset%red_efieldbar(3)
     end do

!    if not convergence well, need to add some code here to make sure efield_test_cart(:) not change much
     dtset%efield(:) = efield_test_cart(:)

   end if  ! berryopt ==14

 end if ! end efield .and. scfcv_level 0 tasks

!!! 
!!! Various printing and update steps for the different efield options
!!!

 if (efield .and. (scfcv_level == 1) ) then ! do this each scf step

   if (dtset%prtvol >= 10)then
     write(message,'(6(a),3(e16.9,2x),a,a,3(e16.9,2x))')ch10,&
&     ' scfcv: New value of the polarization:',ch10,&
&     ' (reduced coordinates, a. u.)',ch10,&
&     '     Electronic berry phase:       ', (pel_cg(ii), ii = 1, 3)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     if(psps%usepaw==1) then
       write(message,'(a,3(e16.9,2x))')&
&       '     ...includes PAW on-site term: ', (pelev(ii), ii = 1, 3)
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if
     write(message,'(a,3(e16.9,2x),a,a,3(e16.9,2x))')&
&     '     Ionic:                        ', (pion(ii), ii = 1, 3), ch10, &
&     '     Total:                        ', (red_ptot(ii), ii = 1, 3) !!REC
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if ! end prtvol >= 10 output

   ptot_cart(:)=zero
   do ii = 1,3
     ptot_cart(ii)=rprimd(ii,1)*red_ptot(1) + rprimd(ii,2)*red_ptot(2) + &
&     rprimd(ii,3)*red_ptot(3)
   end do
   ptot_cart(:)=ptot_cart(:)/ucvol

!  !===================================================================================================
!  !                                       OUTPUT  for fixed E
!  !===================================================================================================

   if (dtset%berryopt == 4) then

!    !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
     write(message,'(a,a)')   ch10, 'scfcv: Constant unreduced E-field:'
     call wrtout(std_out,message,'COLL')
     call prtefield(dtset,dtefield,std_out,rprimd)
     if(dtset%prtvol>=10)then
       call wrtout(ab_out,message,'COLL')
       call prtefield(dtset,dtefield,ab_out,rprimd)
     end if
   end if ! end berryopt 4 output

!  =====================================================================================
!  !                                      fixed D calculation
!  !====================================================================================
   if (dtset%berryopt == 6) then
     if (scfcv_step > 1) then

!      ! update efield taking damping into account dfield is in cartesian in dtset structure (contains input value)
!      ! same goes for efield - update the dtset%efield value
       efield_test_cart(:)=dtset%ddamp*(dtset%dfield(:)-4.0d0*pi*ptot_cart(:))+&
&       (1.0d0-dtset%ddamp)*efield_old_cart(:)

!      ! test whether change in efield in any direction exceed maxestep, if so, set the
!      ! change to maxestep instead   ! need optimized !
       do ii = 1,3

         if (dabs(efield_test_cart(ii)-efield_old_cart(ii)) > dabs(dtset%maxestep)) then

           write(std_out,'(a,a,i5)') "JH - ","  E-field component:",ii
           write(std_out,'(a,es13.5,a,es13.5,a,es13.5,a,es13.5)') " E(n)=",efield_test_cart(ii), &
&           ",    E(n-1)=",efield_old_cart(ii), ",    E(n)-E(n-1)=", efield_test_cart(ii)-efield_old_cart(ii), &
&           ",    maxestep=",dtset%maxestep


           if (efield_test_cart(ii) > efield_old_cart(ii)) then
             efield_test_cart(ii) = efield_old_cart(ii) + dabs(dtset%maxestep)
           else
             efield_test_cart(ii) = efield_old_cart(ii) - dabs(dtset%maxestep)
           end if
         end if
       end do

       dtset%efield(:) = efield_test_cart(:)

!      !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
       write(message,'(a,a)')   ch10, 'scfcv: Constant unreduced D-field  - updating E-field:'
       call wrtout(std_out,message,'COLL')
       call prtefield(dtset,dtefield,std_out,rprimd)
       if(dtset%prtvol>=10)then
         call wrtout(ab_out,message,'COLL')
         call prtefield(dtset,dtefield,ab_out,rprimd)
       end if

!      ! need to update dtset%efield_dot(:) with new value
       dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
       dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
       dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

     else

       write(message,'(a,a)')   ch10, 'scfcv: Constant unreduced D-field  - Pre E-field:'
       call wrtout(std_out,message,'COLL')
       call prtefield(dtset,dtefield,std_out,rprimd)
       if(dtset%prtvol>=10)then
         call wrtout(ab_out,message,'COLL')
         call prtefield(dtset,dtefield,ab_out,rprimd)
       end if

     end if  ! scfcv_step >1

     efield_old_cart(:)=dtset%efield(:)
   end if  ! berryopt ==6
!  !===================================================================================================
!  !                                      fixed reduced d calculation
!  !===================================================================================================
   if (dtset%berryopt == 16) then

     if (scfcv_step > 1) then
!      ! update efield taking damping into account reduced red_dfield
!      red_efield2 is reduced electric field, defined by Eq.(25) of Nat. Phys. suppl. (2009)

       red_efield2(:)=dtset%ddamp*(dtset%red_dfield(:)-red_ptot(:))+ (1.0d0-dtset%ddamp)*red_efield2_old(:)

!      to calculate unreduced E
       efield_test_cart(:)=(4*pi/ucvol)*(rprimd(:,1)*red_efield2(1)+rprimd(:,2)*red_efield2(2)+rprimd(:,3)*red_efield2(3))

!      ! test whether change in efield in any direction exceed maxestep, if so, set the
!      ! change to maxestep instead   ! need optimized !
       do ii = 1,3

         if (dabs(efield_test_cart(ii)-efield_old_cart(ii)) > dabs(dtset%maxestep)) then

           write(std_out,'(a,a,i5)') "JH - ","  E-field component:",ii
           write(std_out,'(a,es13.5,a,es13.5,a,es13.5,a,es13.5)') " E(n)=",efield_test_cart(ii), &
&           ",    E(n-1)=",efield_old_cart(ii), ",    E(n)-E(n-1)=", efield_test_cart(ii)-efield_old_cart(ii), &
&           ",    maxestep=",dtset%maxestep

           if (efield_test_cart(ii) > efield_old_cart(ii)) then
             efield_test_cart(ii) = efield_old_cart(ii) + dabs(dtset%maxestep)
           else
             efield_test_cart(ii) = efield_old_cart(ii) - dabs(dtset%maxestep)
           end if
         end if
       end do

       dtset%efield(:) = efield_test_cart(:)

!      !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
       write(message,'(a,a)')   ch10, 'scfcv: Constant reduced d-field  - updating E-field:'
       call wrtout(std_out,message,'COLL')
       call prtefield(dtset,dtefield,std_out,rprimd)
       if(dtset%prtvol>=10)then
         call wrtout(ab_out,message,'COLL')
         call prtefield(dtset,dtefield,ab_out,rprimd)
       end if

!      ! need to update dtset%efield_dot(:) with new value
!      ! This needs to be deleted  when efield_dot is deleted
       dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
       dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
       dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

     else

       write(message,'(a,a)')   ch10, 'scfcv: Constant reduced d-field  - Pre E-field:'
       call wrtout(std_out,message,'COLL')
       call prtefield(dtset,dtefield,std_out,rprimd)
       if(dtset%prtvol>=10)then
         call wrtout(ab_out,message,'COLL')
         call prtefield(dtset,dtefield,ab_out,rprimd)
       end if

     end if  ! scfcv_step > 1

     efield_old_cart(:)=dtset%efield(:)
     red_efield2_old(:)=red_efield2(:)
   end if  ! berryopt ==16


!  !===================================================================================================
!  !                                      fixed reduced d and ebar calculation (mixed BC)
!  !===================================================================================================
   if (dtset%berryopt == 17) then

     if (scfcv_step > 1) then
!      ! update efield taking damping into account reduced red_dfield
!      red_efield1 and red_efield2 is reduced electric field, defined by Eq.(25) of Nat. Phys. suppl. (2009)
!      red_efield1 for fixed ebar, red_efield2 for fixed d calculation

!      save this value in order to print the final value of real electric field, comparing with the desired red_fieldbar
       dtefield%efield2(:)=dtset%efield(:)

!      write(*,'(a,3i4)') "jfielddir=", (dtset%jfielddir(ii),ii=1,3)

       do ii=1,3
         if (dtset%jfielddir(ii) ==2 ) then    ! direction under fixed d
           dtset%red_efieldbar(ii) = dot_product(dtset%efield(:),rprimd(:,ii)) !  update ebar which is not fixed
           dtefield%efield_dot(ii) = dot_product(dtset%efield(:),rprimd(:,ii))
           red_efield2(ii)=dtset%ddamp*(dtset%red_dfield(ii) - red_ptot(ii)) +  &
&           (1.0d0-dtset%ddamp)*red_efield2_old(ii)         ! d(ii) is fixed, update e(ii)  may need ddamping here

!          write(message,'(a,a,i5,a,i5)')   ch10, 'direction  ', ii,'   for fixed d, value is (2)  ', dtset%jfielddir(ii)
!          call wrtout(ab_out,message,'COLL')
!          call wrtout(std_out,message,'COLL')

         else if (dtset%jfielddir(ii) ==1 ) then   ! direction under fixed ebar
           red_efield2(ii)= (ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,ii)) !  update e which is not fixed
           dtset%red_dfield(ii)=red_ptot(ii) +  (ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,ii))  ! update d

!          write(message,'(a,a,i5,a,i5)')   ch10, 'direction  ', ii,'   for fixed ebar, value is (1)  ', dtset%jfielddir(ii)
!          call wrtout(ab_out,message,'COLL')
!          call wrtout(std_out,message,'COLL')

         end if
       end do

       do ii=1,3
         red_efield1(ii)  =(ucvol/(4*pi))*dot_product(dtset%red_efieldbar(:),gmet(:,ii))
       end do


       dtset%red_efield(:)=(red_efield1(:) + red_efield2(:))/2.0d0 ! average reduced efield, 
!      one is from fixed ebar part, 
!      the other is from fixed d part. 
!      This may need to be optimized !!

       write(message,'(a,a,a,a,3(es16.9,2x),a)')   ch10, 'Reduced efield from fixed ebar:', ch10, &
&       '       e:  ', (red_efield1(ii),ii=1,3), ch10

!      call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')

       write(message,'(a,a,a,a,3(es16.9,2x),a)')   ch10, 'Reduced efield from fixed d:', ch10, &
&       '       e:  ', (red_efield2(ii),ii=1,3), ch10

!      call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')

       write(message,'(a,a,a,a,3(es16.9,2x),a)')   ch10, 'Average reduced efield:', ch10, &
&       '       e:  ', (dtset%red_efield(ii),ii=1,3), ch10

!      call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')

!      to calculate unreduced E
       do ii=1,3
         efield_test_cart(ii)  = (4*pi/ucvol)* dot_product(dtset%red_efield(:),rprimd(:,ii))
       end do

!      ! test whether change in efield in any direction exceed maxestep, if so, set the
!      ! change to maxestep instead   ! need optimized !
       do ii = 1,3
         if (dabs(efield_test_cart(ii)-efield_old_cart(ii)) > dabs(dtset%maxestep)) then

           write(std_out,'(a,a,i5)') "JH - ","  E-field component:",ii
           write(std_out,'(a,es13.5,a,es13.5,a,es13.5,a,es13.5)') " E(n)=",efield_test_cart(ii), &
&           ",    E(n-1)=",efield_old_cart(ii), ",    E(n)-E(n-1)=", efield_test_cart(ii)-efield_old_cart(ii), &
&           ",    maxestep=",dtset%maxestep

           if (efield_test_cart(ii) > efield_old_cart(ii)) then
             efield_test_cart(ii) = efield_old_cart(ii) + dabs(dtset%maxestep)
           else
             efield_test_cart(ii) = efield_old_cart(ii) - dabs(dtset%maxestep)
           end if
         end if
       end do

       dtset%efield(:) = efield_test_cart(:)

!      !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
       write(message,'(a,a)')   ch10, 'scfcv: Constant reduced ebar and d-field  - updating E-field:'
       call wrtout(std_out,message,'COLL')
       call prtefield(dtset,dtefield,std_out,rprimd)
       if(dtset%prtvol>=10)then
         call wrtout(ab_out,message,'COLL')
         call prtefield(dtset,dtefield,ab_out,rprimd)
       end if


!      ! need to update dtset%efield_dot(:) with new value
!      ! This needs to be deleted  when efield_dot is deleted
       dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
       dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
       dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

     else

       write(message,'(a,a)')   ch10, 'scfcv: Constant reduced ebar and d-field  - Pre E-field:'
       call wrtout(std_out,message,'COLL')
       call prtefield(dtset,dtefield,std_out,rprimd)
       if(dtset%prtvol>=10)then
         call wrtout(ab_out,message,'COLL')
         call prtefield(dtset,dtefield,ab_out,rprimd)
       end if

     end if  ! scfcv_step > 1

     efield_old_cart(:)=dtset%efield(:)
     red_efield2_old(:)=red_efield2(:)

   end if  ! berryopt ==17

 end if ! end efield .and. scfcv_level 1 tasks

!deallocate cprj
 if ( efield .and. psps%usepaw == 1) then
   call pawcprj_free(cprj)
 end if
 ABI_DATATYPE_DEALLOCATE(cprj)


!DEBUG
!write(std_out,*)'update_e_field_vars exit'
!END_DEBUG
end subroutine update_e_field_vars
!!***
