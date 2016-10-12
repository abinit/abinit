!{\src2tex{textfont=tt}}
!!****f* ABINIT/pspini
!! NAME
!! pspini
!!
!! FUNCTION
!! Looping over atom types 1 ... ntypat,
!! read pseudopotential data filename, then call pspatm for each psp.
!! Might combine the psps to generate pseudoatoms, thanks to alchemy.
!! Also compute ecore=[Sum(i) zion(i)] * [Sum(i) epsatm(i)] by calling pspcor.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | iscf=parameter controlling scf or non-scf calculations
!!   | ixc=exchange-correlation choice as input to main routine
!!   | natom=number of atoms in unit cell
!!   | pawxcdev=choice of XC development in PAW formalism
!!   | prtvol= control output volume
!!   | typat(natom)=type (integer) for each atom
!!   |              main routine, for each type of atom
!!  gsqcut=cutoff for G^2 based on ecut for basis sphere (bohr^-2)
!!  gsqcutdg=PAW only - cutoff for G^2 based on ecutdg (fine grid) for basis sphere (bohr^-2)
!!  level= level of the calling routine
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!   used to estimate real space mesh (if necessary)
!!
!! OUTPUT
!!  ecore=total psp core correction energy*ucvol (hartree*bohr^3)
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  gencond=general condition for new computation of pseudopotentials
!!          (if gencond=1, new psps have been re-computed)
!!
!! SIDE EFFECTS
!!  psps <type(pseudopotential_type)>=at output, psps is completely initialized
!!   At the input, it is already partially or completely initialized.
!!
!! NOTES
!! The interplay with the multi-dataset mode is interesting :
!! the pseudopotentials
!! are independent of the dataset, but the largest q vector, the
!! spin-orbit characteristics, the use of Ylm as well as ixc
!! play a role in the set up of pseudopotentials (ixc plays a very minor
!! role, however). So, the pseudopotential data ought not be recomputed
!! when gsqcut, gsqcutdg, mqgrid_ff, mqgrid_vl, npspso, ixc, dimekb and useylm do not change.
!! In many cases, this routine is also called just to write the psp line
!! of the header, without reading again the psp. This psp line
!! is constant throughout run.
!!
!! PARENTS
!!      bethe_salpeter,eph,gstate,nonlinear,respfn,screening,sigma,wfk_analyze
!!
!! CHILDREN
!!      nctab_free,nctab_init,nctab_mixalch,pawtab_set_flags,pspatm,pspcor
!!      psps_ncwrite,psps_print,timab,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pspini(dtset,dtfil,ecore,gencond,gsqcut,gsqcutdg,level,pawrad,pawtab,psps,rprimd,comm_mpi)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_psps,          only : psps_print, psps_ncwrite, nctab_init, nctab_free, nctab_mixalch
 use m_pawrad,        only : pawrad_type
 use m_pawtab,        only : pawtab_type, pawtab_set_flags

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pspini'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_64_psp, except_this_one => pspini
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: level
 integer, optional,intent(in) :: comm_mpi
 integer,intent(out) :: gencond
 real(dp),intent(in) :: gsqcut,gsqcutdg
 real(dp),intent(out) :: ecore
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
!arrays
 real(dp),intent(in) :: rprimd(3,3)
!no_abirules
 type(pseudopotential_type), target,intent(inout) :: psps
 type(pawrad_type), intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type), intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: npspmax=50
 integer,save :: dimekb_old=0,ifirst=1,ixc_old=-1,lmnmax_old=0,lnmax_old=0
 integer,save :: mpssoang_old=0,mqgridff_old=0,mqgridvl_old=0,optnlxccc_old=-1
 integer,save :: paw_size_old=-1,pawxcdev_old=-1,positron_old=-2,usepaw_old=-1
 integer,save :: usexcnhat_old=-1,usewvl_old=-1,useylm_old=-1
 integer :: comm_mpi_,ierr,ii,ilang,ilmn,ilmn0,iproj,ipsp,ipspalch
 integer :: ispin,itypalch,itypat,mtypalch,npsp,npspalch,ntypalch
 integer :: ntypat,ntyppure,paw_size
 logical :: has_kij,has_tproj,has_tvale,has_nabla,has_shapefncg,has_wvl
 real(dp),save :: ecore_old=zero,gsqcut_old=zero,gsqcutdg_old=zero
 real(dp) :: dq,epsatm_psp,qmax,rmax,xcccrc
 character(len=500) :: message
 type(pawrad_type) :: pawrad_dum
 type(pawtab_type) :: pawtab_dum
 type(nctab_t) :: nctab_dum
 type(nctab_t),pointer :: nctab_ptr
!arrays
 integer :: paw_options(8)
 integer,save :: paw_options_old(8)=(/-1,-1,-1,-1,-1,-1,-1,-1/)
 integer,save :: pspso_old(npspmax),pspso_zero(npspmax)
 integer,allocatable :: indlmn_alch(:,:,:),new_pspso(:)
 integer,pointer :: indlmn(:,:)
 real(dp),save :: epsatm(npspmax)
 real(dp) :: tsec(2)
 real(dp),allocatable :: dvlspl(:,:),dvlspl_alch(:,:,:),ekb(:),ekb_alch(:,:)
 real(dp),allocatable :: epsatm_alch(:),ffspl(:,:,:),ffspl_alch(:,:,:,:)
 real(dp),allocatable :: vlspl(:,:),vlspl_alch(:,:,:),xccc1d(:,:)
 real(dp),allocatable :: xccc1d_alch(:,:,:),xcccrc_alch(:)
 type(nctab_t),target,allocatable :: nctab_alch(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Keep track of time spent in this subroutine
 call timab(15,1,tsec)

!-------------------------------------------------------------
! Some initializations
!-------------------------------------------------------------

!Useful sizes
 ntypat=psps%ntypat
 mtypalch=psps%mtypalch
 npsp=psps%npsp
 if (npsp>npspmax) then
   MSG_BUG("npsp>npspmax in pspini !")
 end if

!Size of grids for atomic data represented in reciprocal space

!Set up q grids, make qmax 20% larger than largest expected:
 qmax=1.2d0 * sqrt(gsqcut)
!ffnl is always computed in reciprocal space
 dq=qmax/(one*(psps%mqgrid_ff-1))
 do ii=1,psps%mqgrid_ff
   psps%qgrid_ff(ii)=(ii-1)*dq
 end do
 if (psps%usepaw==1) qmax=1.2d0 * sqrt(gsqcutdg)
!If vlspl is computed in real space, qgrid contains a real space mesh
!the max is taken as the biggest distance in the box.
 if (psps%vlspl_recipSpace) then
   dq=qmax/(one*(psps%mqgrid_vl-1))
 else
   rmax = (rprimd(1, 1) + rprimd(1, 2) + rprimd(1, 3)) ** 2
   rmax = rmax + (rprimd(2, 1) + rprimd(2, 2) + rprimd(2, 3)) ** 2
   rmax = rmax + (rprimd(3, 1) + rprimd(3, 2) + rprimd(3, 3)) ** 2
   rmax = sqrt(rmax)
   dq = rmax /(one*(psps%mqgrid_vl-1))
 end if
 do ii=1,psps%mqgrid_vl
   psps%qgrid_vl(ii)=(ii-1)*dq
 end do

!Determine whether new optional data requests have changed
 paw_options=0;paw_size=0
 if (psps%usepaw==1) then
   paw_size=size(pawtab)
   has_kij=(dtset%positron/=0)
   has_tvale=.true. ! Will be modified later (depending on PAW dataset format)
   has_nabla=.false.
   has_shapefncg=(dtset%optdriver==RUNL_GSTATE.and.((dtset%iprcel>=20.and.dtset%iprcel<70).or.dtset%iprcel>=80))
   has_wvl=(dtset%usewvl==1.or.dtset%icoulomb/=0)
   has_tproj=(dtset%usewvl==1) ! projectors will be free at the end of the psp reading
   if (has_kij)      paw_options(1)=1
   if (has_tvale)    paw_options(2)=1
   if (has_nabla)    paw_options(5)=1
   if (has_shapefncg)paw_options(6)=1
   if (has_wvl)      paw_options(7)=1
   if (has_tproj)    paw_options(8)=1
 end if

!Determine whether the spin-orbit characteristic has changed
!Do not forget that the SO is not consistent with alchemy presently
 ABI_ALLOCATE(new_pspso,(npsp))
 if (ifirst==1) pspso_old(:)=-1
 if (ifirst==1) pspso_zero(:)=-1
 do ipsp=1,npsp
   new_pspso(ipsp)=1
!  No new characteristics if it is equal to the old one,
!  or, if it is one, the old one is equal to the intrinsic characteristic one.
   if (psps%pspso(ipsp)==pspso_old(ipsp).or. &
&   (psps%pspso(ipsp)==1.and.pspso_old(ipsp)==pspso_zero(ipsp))) then
     new_pspso(ipsp)=0
   end if
!  Prepare the saving of the intrinsic pseudopotential characteristics
   if(psps%pspso(ipsp)==1) pspso_zero(ipsp)=0
 end do

!Compute the general condition for new computation of pseudopotentials
 gencond=0
 if(   ixc_old /= dtset%ixc                &
& .or. mqgridff_old /= psps%mqgrid_ff      &
& .or. mqgridvl_old /= psps%mqgrid_vl      &
& .or. mpssoang_old /= psps%mpssoang       &
& .or. abs(gsqcut_old-gsqcut)>1.0d-10      &
& .or. (psps%usepaw==1.and.abs(gsqcutdg_old-gsqcutdg)>1.0d-10) &
& .or. dimekb_old /= psps%dimekb           &
& .or. lmnmax_old /= psps%lmnmax           &
& .or. lnmax_old  /= psps%lnmax            &
& .or. optnlxccc_old /= psps%optnlxccc     &
& .or. usepaw_old /= psps%usepaw           &
& .or. useylm_old /= psps%useylm           &
& .or. pawxcdev_old /= dtset%pawxcdev      &
& .or. positron_old /= dtset%positron      &
& .or. usewvl_old /= dtset%usewvl          &
& .or. paw_size_old /= paw_size            &
& .or. usexcnhat_old/=dtset%usexcnhat_orig &
& .or. any(paw_options_old(:)/=paw_options(:)) &
& .or. sum(new_pspso(:))/=0                &
& .or. mtypalch>0                          &
& .or. (dtset%usewvl==1.and.psps%usepaw==1)&
& ) gencond=1

 if (present(comm_mpi).and.psps%usepaw==1) then
   if(xmpi_comm_size(comm_mpi)>1)then
     call xmpi_sum(gencond,comm_mpi,ierr)
   end if
   if (gencond/=0) gencond=1
 end if
 ABI_DEALLOCATE(new_pspso)

!-------------------------------------------------------------
! Following section is only reached when new computation
! of pseudopotentials is needed
!-------------------------------------------------------------

 if (gencond==1) then

   write(message, '(a,a)' ) ch10,&
&   '--- Pseudopotential description ------------------------------------------------'
   call wrtout(ab_out,message,'COLL')

   ABI_ALLOCATE(ekb,(psps%dimekb*(1-psps%usepaw)))
   ABI_ALLOCATE(xccc1d,(psps%n1xccc*(1-psps%usepaw),6))
   ABI_ALLOCATE(ffspl,(psps%mqgrid_ff,2,psps%lnmax))
   ABI_ALLOCATE(vlspl,(psps%mqgrid_vl,2))
   if (.not.psps%vlspl_recipSpace) then
     ABI_ALLOCATE(dvlspl,(psps%mqgrid_vl,2))
   else
     ABI_ALLOCATE(dvlspl,(0,0))
   end if

!  PAW: reset flags for optional data
   if (psps%usepaw==1) then
     call pawtab_set_flags(pawtab,has_kij=paw_options(1),has_tvale=paw_options(2),&
&     has_vhnzc=paw_options(3),has_vhtnzc=paw_options(4),&
&     has_nabla=paw_options(5),has_shapefncg=paw_options(6),&
&     has_wvl=paw_options(7),has_tproj=paw_options(8))
   end if

!  Read atomic pseudopotential data and get transforms
!  for each atom type : two cases, alchemy or not.

!  No alchemical pseudoatom, in all datasets, npsp=ntypat
   if(mtypalch==0)then

     do ipsp=1,npsp

       xcccrc=zero
       ekb(:)=zero;ffspl(:,:,:)=zero;vlspl(:,:)=zero
       if (.not.psps%vlspl_recipSpace) dvlspl(:, :)=zero
       if (psps%usepaw==0) xccc1d(:,:)=zero
       indlmn=>psps%indlmn(:,:,ipsp)
       indlmn(:,:)=0

       write(message, '(a,i4,a,t38,a)' ) &
&       '- pspini: atom type',ipsp,'  psp file is',trim(psps%filpsp(ipsp))
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')

!      Read atomic psp V(r) and wf(r) to get local and nonlocal psp:
!      Cannot use the same call in case of bound checking, because of pawrad/pawtab
       if(psps%usepaw==0)then
         call pspatm(dq,dtset,dtfil,ekb,epsatm(ipsp),ffspl,indlmn,ipsp,&
&         pawrad_dum,pawtab_dum,psps,vlspl,dvlspl,xcccrc,xccc1d,psps%nctab(ipsp))
         psps%ekb(:,ipsp)=ekb(:)
         psps%xccc1d(:,:,ipsp)=xccc1d(:,:)
       else
         comm_mpi_=xmpi_comm_self;if (present(comm_mpi)) comm_mpi_=comm_mpi
         call pspatm(dq,dtset,dtfil,ekb,epsatm(ipsp),ffspl,indlmn,ipsp,&
&         pawrad(ipsp),pawtab(ipsp),psps,vlspl,dvlspl,xcccrc,xccc1d,nctab_dum,&
&         comm_mpi=comm_mpi_)
       end if

       ! Copy data to psps datastructure.
       psps%xcccrc(ipsp)=xcccrc
       psps%znucltypat(ipsp)=psps%znuclpsp(ipsp)
       psps%ffspl(:,:,:,ipsp)=ffspl(:,:,:)
       psps%vlspl(:,:,ipsp)=vlspl(:,:)
       if (.not.psps%vlspl_recipSpace) psps%dvlspl(:, :, ipsp) = dvlspl(:, :)
     end do ! ipsp

   else ! if mtypalch/=0

     npspalch=psps%npspalch
     ntyppure=npsp-npspalch
     ntypalch=psps%ntypalch
     ABI_ALLOCATE(epsatm_alch,(npspalch))
     ABI_ALLOCATE(ekb_alch,(psps%dimekb,npspalch*(1-psps%usepaw)))
     ABI_ALLOCATE(ffspl_alch,(psps%mqgrid_ff,2,psps%lnmax,npspalch))
     ABI_ALLOCATE(xccc1d_alch,(psps%n1xccc*(1-psps%usepaw),6,npspalch))
     ABI_ALLOCATE(xcccrc_alch,(npspalch))
     ABI_ALLOCATE(vlspl_alch,(psps%mqgrid_vl,2,npspalch))
     if (.not.psps%vlspl_recipSpace) then
       ABI_ALLOCATE(dvlspl_alch,(psps%mqgrid_vl,2,npspalch))
     end if
     ABI_ALLOCATE(indlmn,(6,psps%lmnmax))
     ABI_ALLOCATE(indlmn_alch,(6,psps%lmnmax,npspalch))

     ! Allocate NC tables used for mixing.
     if (psps%usepaw == 0) then
       ABI_DT_MALLOC(nctab_alch, (npspalch))
       do ipspalch=1,npspalch
         call nctab_init(nctab_alch(ipspalch), psps%mqgrid_vl, .False., .False.)
       end do
     end if

     do ipsp=1,npsp
       write(message, '(a,i4,a,t38,a)' ) &
&       '- pspini: atom type',ipsp,'  psp file is',trim(psps%filpsp(ipsp))
       call wrtout(ab_out,message,'COLL')

       xcccrc=zero
       ekb(:)=zero;ffspl(:,:,:)=zero;vlspl(:,:)=zero
       if (.not.psps%vlspl_recipSpace) dvlspl(:, :)=zero
       if (psps%usepaw==0) xccc1d(:,:)=zero
       indlmn(:,:)=0

!      Read atomic psp V(r) and wf(r) to get local and nonlocal psp:
       if (psps%usepaw==0) then
         if (ipsp <= ntyppure) then
           ! Store data in nctab if pure atom.
           nctab_ptr => psps%nctab(ipsp)
         else
           ! Store data in nctab_alch (to be mixed afterwards).
           nctab_ptr => nctab_alch(ipsp-ntyppure)
         end if

         call pspatm(dq,dtset,dtfil,ekb,epsatm_psp,ffspl,indlmn,ipsp,&
&         pawrad_dum,pawtab_dum,psps,vlspl,dvlspl,xcccrc,xccc1d,nctab_ptr)

       else if (psps%usepaw==1) then
         comm_mpi_=xmpi_comm_self;if (present(comm_mpi)) comm_mpi_=comm_mpi
         call pspatm(dq,dtset,dtfil,ekb,epsatm_psp,ffspl,indlmn,ipsp,&
&         pawrad(ipsp),pawtab(ipsp),psps,vlspl,dvlspl,xcccrc,xccc1d,nctab_dum,&
&         comm_mpi=comm_mpi_)
       end if

       if (ipsp<=ntyppure) then
!        Pure pseudopotentials, leading to pure pseudoatoms
         epsatm(ipsp)=epsatm_psp
         psps%znucltypat(ipsp)=psps%znuclpsp(ipsp)
         if (psps%usepaw==0) psps%ekb(:,ipsp)=ekb(:)
         psps%ffspl(:,:,:,ipsp)=ffspl(:,:,:)
         psps%vlspl(:,:,ipsp)=vlspl(:,:)
         if (.not.psps%vlspl_recipSpace) psps%dvlspl(:, :, ipsp)=dvlspl(:, :)
         if (psps%usepaw==0) psps%xccc1d(:,:,ipsp)=xccc1d(:,:)
         psps%xcccrc(ipsp)=xcccrc
         psps%indlmn(:,:,ipsp)=indlmn(:,:)

       else
!        Pseudopotentials for alchemical generation
         ipspalch=ipsp-ntyppure
         epsatm_alch(ipspalch)=epsatm_psp
         ffspl_alch(:,:,:,ipspalch)=ffspl(:,:,:)
         vlspl_alch(:,:,ipspalch)=vlspl(:,:)
         if (.not.psps%vlspl_recipSpace) dvlspl_alch(:,:,ipspalch)=dvlspl(:,:)
         if (psps%usepaw==0) then
           ekb_alch(:,ipspalch)=ekb(:)
           xccc1d_alch(:,:,ipspalch)=xccc1d(:,:)
         end if
         xcccrc_alch(ipspalch)=xcccrc
         indlmn_alch(:,:,ipspalch)=indlmn(:,:)
!        write(std_out,'(a,6i4)' )' pspini : indlmn_alch(:,1,ipspalch)=',indlmn_alch(:,1,ipspalch)
!        write(std_out,'(a,6i4)' )' pspini : indlmn_alch(:,2,ipspalch)=',indlmn_alch(:,2,ipspalch)
       end if

     end do ! ipsp

     ! Generate data for alchemical pseudos.
     do itypalch=1,ntypalch
       itypat=itypalch+ntyppure
       psps%znucltypat(itypat)=200.0+itypalch    ! Convention for alchemical pseudoatoms
       vlspl(:,:)=zero
       if (.not.psps%vlspl_recipSpace) dvlspl(:, :) = zero
       epsatm(itypat)=zero
       xcccrc=zero
       if (psps%usepaw==0) xccc1d(:,:)=zero

!      Here, linear combination of the quantities
!      MG: FIXME I think that the mixing of xcccrc is wrong when the xxccrc are different!
!      but this is minor bug since alchemical pseudos should not have XCCC (?)
       do ipspalch=1,npspalch
         epsatm(itypat) = epsatm(itypat) + epsatm_alch(ipspalch) * psps%mixalch(ipspalch,itypalch)
         vlspl(:,:) = vlspl(:,:) + vlspl_alch(:,:,ipspalch) * psps%mixalch(ipspalch,itypalch)
         if (.not.psps%vlspl_recipSpace) then
           dvlspl(:,:) = dvlspl(:,:) + dvlspl_alch(:,:,ipspalch) * psps%mixalch(ipspalch,itypalch)
         end if
         xcccrc = xcccrc + xcccrc_alch(ipspalch) * psps%mixalch(ipspalch,itypalch)
         if (psps%usepaw==0) then
           xccc1d(:,:) = xccc1d(:,:) + xccc1d_alch(:,:,ipspalch) * psps%mixalch(ipspalch,itypalch)
         end if
       end do ! ipspalch

       psps%vlspl(:,:,itypat)=vlspl(:,:)
       if (.not.psps%vlspl_recipSpace) psps%dvlspl(:, :, itypat) = dvlspl(:, :)
       if (psps%usepaw==0) then 
         psps%xccc1d(:,:,itypat)=xccc1d(:,:)
       end if
       psps%xcccrc(itypat)=xcccrc

       if (abs(xcccrc) > tol6) then
         write(std_out, *)"xcccrc", xcccrc
         MSG_WARNING("Alchemical pseudopotential with nlcc!")
       end if

!      Combine the different non-local projectors : for the scalar part then
!      the spin-orbit part, treat the different angular momenta
!      WARNING : this coding does not work for PAW
       ilmn=0
       psps%indlmn(:,:,itypat)=0
       do ispin=1,2
         do ilang=0,3
           if(ispin==2 .and. ilang==0)cycle
           iproj=0
           do ipspalch=1,npspalch
             if(abs(psps%mixalch(ipspalch,itypalch))>tol10)then
               do ilmn0=1,psps%lmnmax
                 if(indlmn_alch(5,ilmn0,ipspalch)/=0)then
                   if(indlmn_alch(6,ilmn0,ipspalch)==ispin)then
                     if(indlmn_alch(1,ilmn0,ipspalch)==ilang)then
                       ilmn=ilmn+1         ! increment the counter
                       iproj=iproj+1       ! increment the counter, this does not work for PAW
                       if(ilmn>psps%lmnmax)then
                         MSG_BUG('Problem with the alchemical pseudopotentials : ilmn>lmnmax.')
                       end if
                       psps%indlmn(1,ilmn,itypat)=ilang
                       psps%indlmn(2,ilmn,itypat)=indlmn_alch(2,ilmn0,ipspalch)
                       psps%indlmn(3,ilmn,itypat)=iproj                       ! This does not work for PAW
                       psps%indlmn(4,ilmn,itypat)=ilmn                        ! This does not work for PAW
                       psps%indlmn(5,ilmn,itypat)=ilmn
                       psps%indlmn(6,ilmn,itypat)=ispin
                       ! The two lines below do not work for PAW
                       if (psps%usepaw==0) then
                         psps%ekb(ilmn,itypat)=psps%mixalch(ipspalch,itypalch) *ekb_alch(ilmn0,ipspalch) 
                       end if
                       psps%ffspl(:,:,ilmn,itypat)=ffspl_alch(:,:,ilmn0,ipspalch)
                     end if ! ilang is OK
                   end if ! ispin is OK
                 end if ! ilmn0 exist
               end do ! ilmn0
             end if ! mixalch>tol10
           end do ! ipspalch
         end do ! ilang
       end do ! ispin

     end do ! itypalch

     ABI_DEALLOCATE(epsatm_alch)
     ABI_DEALLOCATE(ekb_alch)
     ABI_DEALLOCATE(ffspl_alch)
     ABI_DEALLOCATE(xccc1d_alch)
     ABI_DEALLOCATE(xcccrc_alch)
     ABI_DEALLOCATE(vlspl_alch)
     if (.not.psps%vlspl_recipSpace) then
       ABI_DEALLOCATE(dvlspl_alch)
     end if
     ABI_DEALLOCATE(indlmn_alch)
     ABI_DEALLOCATE(indlmn)

     ! Mix NC tables.
     if (psps%usepaw == 0) then
       call nctab_mixalch(nctab_alch, npspalch, ntypalch, psps%algalch, psps%mixalch, psps%nctab(ntyppure+1:))
       do ipspalch=1,npspalch
         call nctab_free(nctab_alch(ipspalch))
       end do
       ABI_DT_FREE(nctab_alch)
     end if
   end if ! mtypalch

   ABI_DEALLOCATE(ekb)
   ABI_DEALLOCATE(ffspl)
   ABI_DEALLOCATE(vlspl)
   ABI_DEALLOCATE(xccc1d)

   if (.not.psps%vlspl_recipSpace) then
     ABI_DEALLOCATE(dvlspl)
   end if

 end if !  End condition of new computation needed

!-------------------------------------------------------------
! Following section is always executed
!-------------------------------------------------------------
!One should move this section of code outside of pspini,
!but epsatm is needed, so should be in the psp datastructure.
!Compute pseudo correction energy. Will differ from an already
!computed one if the number of atom differ ...
 call pspcor(ecore,epsatm,dtset%natom,ntypat,dtset%typat,psps%ziontypat)
 if(abs(ecore_old-ecore)>tol8*abs(ecore_old+ecore))then
   write(message, '(2x,es15.8,t50,a)' ) ecore,'ecore*ucvol(ha*bohr**3)'
!  ecore is useless if iscf<=0, but at least it has been initialized
   if(dtset%iscf>=0)then
     call wrtout(ab_out,message,'COLL')
   end if
   call wrtout(std_out,message,'COLL')
 end if

!End of pseudopotential output section
 write(message, '(2a)' )&
& '--------------------------------------------------------------------------------',ch10
 call wrtout(ab_out,message,'COLL')

!-------------------------------------------------------------
! Keep track of this call to the routine
!-------------------------------------------------------------

 if (ifirst==1) ifirst=0

 mqgridff_old=psps%mqgrid_ff
 mqgridvl_old=psps%mqgrid_vl
 mpssoang_old=psps%mpssoang
 ixc_old=dtset%ixc
 gsqcut_old=gsqcut;if (psps%usepaw==1) gsqcutdg_old=gsqcutdg
 lmnmax_old=psps%lmnmax
 lnmax_old=psps%lnmax
 optnlxccc_old=psps%optnlxccc
 usepaw_old=psps%usepaw
 dimekb_old=psps%dimekb
 useylm_old=psps%useylm
 pawxcdev_old=dtset%pawxcdev
 positron_old=dtset%positron
 usewvl_old = dtset%usewvl
 usexcnhat_old=dtset%usexcnhat_orig
 paw_size_old=paw_size
 ecore_old=ecore
 paw_options_old(:)=paw_options(:)

 do ipsp=1,npsp
   pspso_old(ipsp)=psps%pspso(ipsp)
   if(pspso_zero(ipsp)==0)pspso_zero(ipsp)=psps%pspso(ipsp)
 end do
 psps%mproj = maxval(psps%indlmn(3,:,:))

 if (gencond == 1) call psps_print(psps,std_out,dtset%prtvol)

 ! Write the PSPS.nc file and exit here if requested by the user.
 if (abs(dtset%prtpsps) == 1) then
   if (xmpi_comm_rank(xmpi_world) == 0) call psps_ncwrite(psps, trim(dtfil%filnam_ds(4))//"_PSPS.nc")
   if (dtset%prtpsps == -1) then
     MSG_ERROR_NODUMP("prtpsps == -1 ==> aborting now")
   end if
 end if

 call timab(15,2,tsec)

 DBG_EXIT("COLL")

end subroutine pspini
!!***
