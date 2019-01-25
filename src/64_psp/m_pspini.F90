!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pspini
!! NAME
!!  m_pspini
!!
!! FUNCTION
!!  Initialize pseudopotential datastructures from files.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, MT, FrD, AF, DRH)
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

module m_pspini

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_abicore
 use m_xmpi
 use m_psxml2ab
 !use m_psps

 use m_time,      only : timab
 use m_io_tools,  only : open_file
 use m_pawrad,    only : pawrad_type
 use m_pawtab,    only : pawtab_type, pawtab_set_flags
 use m_psps,      only : psps_print, psps_ncwrite, nctab_init, nctab_free, nctab_mixalch, test_xml_xmlpaw_upf, &
                         nctab_eval_tcorespl
 use m_pawpsp,    only : pawpsp_bcast, pawpsp_read_pawheader, pawpsp_read_header_xml,&
                         pawpsp_header_type, pawpsp_wvl, pawpsp_7in, pawpsp_17in
 use m_pawxmlps,  only : paw_setup_free,paw_setuploc
 use m_pspheads,  only : pawpsxml2ab
#if defined HAVE_BIGDFT
 use BigDFT_API, only : dictionary, atomic_info, dict_init, dict_free, UNINITIALIZED
#endif

 use m_psp1,       only : psp1in
 use m_psp5,       only : psp5in
 use m_psp6,       only : psp6in
 use m_psp8,       only : psp8in
 use m_psp9,       only : psp9in
 use m_upf2abinit, only : upf2abinit
 use m_psp_hgh,    only : psp2in, psp3in, psp10in
 use m_wvl_descr_psp,  only : wvl_descr_psp_fill

 implicit none

 private
!!***

 public :: pspini
!!***

contains
!!***

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
!! Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, MT)
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

subroutine pspini(dtset,dtfil,ecore,gencond,gsqcut,gsqcutdg,pawrad,pawtab,psps,rprimd,comm_mpi)

!Arguments ------------------------------------
!scalars
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
 logical :: has_kij,has_tproj,has_tvale,has_nabla,has_shapefncg,has_vminushalf,has_wvl
 real(dp),save :: ecore_old=zero,gsqcut_old=zero,gsqcutdg_old=zero
 real(dp) :: dq,epsatm_psp,qmax,rmax,xcccrc
 character(len=500) :: message
 type(pawrad_type) :: pawrad_dum
 type(pawtab_type) :: pawtab_dum
 type(nctab_t) :: nctab_dum
 type(nctab_t),pointer :: nctab_ptr
!arrays
 integer :: paw_options(9)
 integer,save :: paw_options_old(9)=(/-1,-1,-1,-1,-1,-1,-1,-1,-1/)
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
   has_vminushalf=(maxval(dtset%ldaminushalf)==1)
   if (has_kij)       paw_options(1)=1
   if (has_tvale)     paw_options(2)=1
   if (has_nabla)     paw_options(5)=1
   if (has_shapefncg) paw_options(6)=1
   if (has_wvl)       paw_options(7)=1
   if (has_tproj)     paw_options(8)=1
   if (has_vminushalf)paw_options(9)=1
   !if (dtset%prtvclmb /= 0) then
   paw_options(3) = 1
   paw_options(4) = 1
   !end if
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
!  No new characteristics if PAW
   if (psps%usepaw==1) new_pspso(ipsp)=0
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
! the following have to be included in pawtab_set_flags
     do ipsp=1,psps%ntypat
       pawtab(ipsp)%has_vminushalf=dtset%ldaminushalf(ipsp)
     end do
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

!!****f* ABINIT/pspcor
!! NAME
!! pspcor
!!
!! FUNCTION
!! Compute ecore pseudoion-pseudoion correction energy from epsatm for
!! different types of atoms in unit cell.
!!
!! INPUTS
!!  natom=number of atoms in cell
!!  ntypat=number of types of atoms
!!  typat(natom)=integer label of 'typat' for each atom in cell
!!  epsatm(ntypat)=pseudoatom energy for each type of atom
!!  zion(ntypat)=valence charge on each type of atom in cell
!!
!! OUTPUT
!!  ecore=resulting psion-psion energy in Hartrees
!!
!! PARENTS
!!      pspini
!!
!! CHILDREN
!!
!! SOURCE

subroutine pspcor(ecore,epsatm,natom,ntypat,typat,zion)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 real(dp),intent(out) :: ecore
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: epsatm(ntypat),zion(ntypat)

!Local variables-------------------------------
!scalars
 integer :: ia
 real(dp) :: charge,esum

! *************************************************************************

 charge = 0.d0
 esum = 0.d0
 do ia=1,natom
!  compute pseudocharge:
   charge=charge+zion(typat(ia))
!  add pseudocore energies together:
   esum = esum + epsatm(typat(ia))
 end do

 ecore=charge*esum

end subroutine pspcor
!!***

!!****f* ABINIT/pspatm
!! NAME
!! pspatm
!!
!! FUNCTION
!! Open atomic pseudopotential data file for a given atom,
!! read the three first lines, make some checks, then
!! call appropriate subroutine for the reading of
!! V(r) and wf R(r) data for each angular momentum, and subsequent
!! Fourier and Bessel function transforms for local and nonlocal potentials.
!! Close psp file at end.
!!
!! Handles pseudopotential files produced by (pspcod=1 or 4) Teter code,
!! or from the Goedecker-Teter-Hutter paper (pspcod=2),
!! or from the Hartwigsen-Goedecker-Hutter paper (pspcod=3 or 10)
!! or "Phoney pseudopotentials" (Hamman grid in real space) (pspcod=5)
!! or "Troullier-Martins pseudopotentials" from the FHI (pspcod=6)
!! or "XML format" (pspcod=9)
!! or "UPF PWSCF format" (pspcod=11)
!!
!! INPUTS
!!  dq= spacing of the q-grid
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | ixc=exchange-correlation choice from main routine data file
!!   | pawxcdev=choice of XC development in PAW formalism
!!   | usexcnhat_orig=choice for use of nhat in Vxc in PAW formalism
!!   | xclevel= XC functional level
!!  ipsp=id in the array of the currently read pseudo.
!!
!! OUTPUT
!!  ekb(dimekb)=
!!    ->NORM-CONSERVING PSPS ONLY (pspcod/=7):
!!      (Real) Kleinman-Bylander energies (hartree)
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for number of basis functions (l,n) (dimekb=lnmax)
!!             If any, spin-orbit components begin at l=mpsang+1
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$(hartree)
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  vlspl(mqgrid_vl,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  ffspl(mqgrid_ff,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector; if any, spin-orbit components begin at l=mpsang+1
!!  xcccrc=XC core correction cutoff radius (bohr) from psp file
!!  xccc1d(n1xccc*(1-usepaw),6)=1D core charge function and five derivatives, from psp file (used in NC only)
!!  nctab=<nctab_t>
!!    has_tvale=True if the pseudo provides the valence density (used in NC only)
!!    tvalespl(mqgrid_vl(1-usepaw),2)=the pseudo valence density and 2nd derivative in reciprocal space on a regular grid
!!                                     (used in NC only)
!!
!! SIDE EFFECTS
!! Input/Output :
!!  psps <type(pseudopotential_type)>=at output, values depending on the read
!!                                    pseudo are set.
!!   | dimekb(IN)=dimension of ekb (see module defs_psp.f)
!!   | filpsp(IN)=name of formatted external file containing atomic psp data.
!!   | lmnmax(IN)=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!   |           =if useylm=0, max number of (l,n)   comp. over all type of psps
!!   | lnmax(IN)=max. number of (l,n) components over all type of psps
!!   |           angular momentum of nonlocal pseudopotential
!!   | mpsang(IN)= 1+maximum angular momentum for nonlocal pseudopotentials
!!   | mpssoang(IN)= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!   | mqgrid_ff(IN)=dimension of q (or G) grid for nl form factors (array ffspl)
!!   | mqgrid_vl(IN)=dimension of q (or G) grid for Vloc (array vlspl)
!!   | n1xccc(IN)=dimension of xccc1d ; 0 if no XC core correction is used
!!   | optnlxccc(IN)=option for nl XC core correction
!!   | positron(IN)=0 if electron GS calculation
!!   |              1 if positron GS calculation
!!   |              2 if electron GS calculation in presence of the positron
!!   | pspso(INOUT)=spin-orbit characteristics, govern the content of ffspl and ekb
!!   |          if =0 : this input requires NO spin-orbit characteristics of the psp
!!   |          if =2 : this input requires HGH characteristics of the psp
!!   |          if =3 : this input requires HFN characteristics of the psp
!!   |          if =1 : this input will be changed at output to 1, 2, 3, according
!!   |                  to the intrinsic characteristics of the psp file
!!   | qgrid_ff(mqgrid_ff)(IN)=values of q on grid from 0 to qmax (bohr^-1) for nl form factors
!!   | qgrid_vl(mqgrid_vl)(IN)=values of q on grid from 0 to qmax (bohr^-1) for Vloc
!!   | usepaw(IN)= 0 for non paw calculation; =1 for paw calculation
!!   | useylm(IN)=governs the way the nonlocal operator is to be applied:
!!   |            1=using Ylm, 0=using Legendre polynomials
!!   | vlspl_recipSpace(IN)=.true. if pseudo are expressed in reciprocal space.
!!   | znuclpsp(IN)=atomic number of atom as specified in input file to main routine
!!
!! NOTES
!!  Format expected for the three first lines of pseudopotentials
!!  (1) title (character) line
!!  (2) znucl,zion,pspdat
!!  (3) pspcod,pspxc,lmax,lloc,mmax,r2well
!!
!!  Dimensions of form factors and Vloc q grids must be the same in Norm-Conserving case
!!
!! PARENTS
!!      pspini
!!
!! CHILDREN
!!      nctab_eval_tcorespl,pawpsp_17in,pawpsp_7in,pawpsp_bcast
!!      pawpsp_read_header_xml,pawpsp_read_pawheader,pawpsp_wvl,psp10in,psp1in
!!      psp2in,psp3in,psp5in,psp6in,psp8in,psp9in,psp_dump_outputs
!!      psxml2abheader,test_xml_xmlpaw_upf,timab,upf2abinit,wrtout
!!      wvl_descr_psp_fill
!!
!! SOURCE

subroutine pspatm(dq,dtset,dtfil,ekb,epsatm,ffspl,indlmn,ipsp,pawrad,pawtab,&
&  psps,vlspl,dvlspl,xcccrc,xccc1d,nctab,comm_mpi)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipsp
 integer, optional,intent(in) :: comm_mpi
 real(dp),intent(in) :: dq
 real(dp),intent(out) :: epsatm,xcccrc
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(pawrad_type),intent(inout) :: pawrad
 type(pawtab_type),intent(inout) :: pawtab
 type(nctab_t),intent(inout) :: nctab
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,intent(out) :: indlmn(6,psps%lmnmax)
 real(dp),intent(out) :: dvlspl(psps%mqgrid_vl,2)
 real(dp),intent(inout) :: ekb(psps%dimekb*(1-psps%usepaw))
 real(dp),intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
 real(dp),intent(out) :: vlspl(psps%mqgrid_vl,2)
 real(dp),intent(inout) :: xccc1d(psps%n1xccc*(1-psps%usepaw),6)

!Local variables ---------------------------------------
!scalars
 integer :: ii,il,ilmn,iln,iln0,lloc,lmax,me,mmax
 integer :: paral_mode,pspcod,pspdat,pspxc,useupf,usexml,xmlpaw,unt
 real(dp) :: maxrad,qchrg,r2well,zion,znucl
 character(len=500) :: message,errmsg
 character(len=fnlen) :: title
 character(len=fnlen) :: filnam
 type(pawpsp_header_type):: pawpsp_header
 type(pspheader_type) :: pspheads_tmp
!arrays
 integer,allocatable :: nproj(:)
 real(dp) :: tsec(2),ecut_tmp(3,2)
 real(dp),allocatable :: e990(:),e999(:),ekb1(:),ekb2(:),epspsp(:),rcpsp(:)
 real(dp),allocatable :: rms(:)
#if defined HAVE_PSML
!!  usexml= 0 for non xml ps format ; =1 for xml ps format
 character(len=3) :: atmsymb
 character(len=30) :: creator
 type(pspheader_type) :: psphead
#endif

! ******************************************************************************

!paral_mode defines how we access to the psp file
!  paral_mode=0: all processes access to the file (sequentially)
!  paral_mode=1: only proc. 0 access to the file and then broadcast
 paral_mode=0
 if (present(comm_mpi)) then
   if (psps%usepaw==1.and.xmpi_comm_size(comm_mpi)>1) paral_mode=1
 end if
 me=0;if (paral_mode==1) me=xmpi_comm_rank(comm_mpi)

 if (paral_mode == 1) then
   ABI_CHECK(psps%usepaw==1, "paral_mode==1 is only compatible with PAW, see call to pawpsp_bcast below")
 end if

 nctab%has_tvale = .False.; nctab%has_tcore = .False.

 if (me==0) then
!  Dimensions of form factors and Vloc q grids must be the same in Norm-Conserving case
   if (psps%usepaw==0 .and. psps%mqgrid_ff/=psps%mqgrid_vl) then
     write(message, '(a,a,a,a,a)' )&
&     'Dimension of q-grid for nl form factors (mqgrid_ff)',ch10,&
&     'is different from dimension of q-grid for Vloc (mqgrid_vl) !',ch10,&
&     'This is not allowed for norm-conserving psp.'
     MSG_ERROR(message)
   end if

   write(message, '(a,t38,a)' )'- pspatm: opening atomic psp file',trim(psps%filpsp(ipsp))
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   !write(message, "(2a)")"- md5: ",trim(psps%md5_pseudos(ipsps))

   !  Check if the file pseudopotential file is written in (XML| XML-PAW | UPF)
   call test_xml_xmlpaw_upf(psps%filpsp(ipsp), usexml, xmlpaw, useupf)

!  ----------------------------------------------------------------------------
!  allocate nproj here: can be read in now for UPF
   ABI_ALLOCATE(nproj,(psps%mpssoang))
   nproj(:)=0

   if (usexml /= 1 .and. useupf /= 1) then

!    Open the atomic data file, and read the three first lines
!    These three first lines have a similar format in all allowed psp files

!    Open atomic data file (note: formatted input file)
     if (open_file(psps%filpsp(ipsp), message, unit=tmp_unit, form='formatted', status='old') /= 0) then
       MSG_ERROR(message)
     end if
     rewind (unit=tmp_unit,err=10,iomsg=errmsg)

!    Read and write some description of file from first line (character data)
     read (tmp_unit,'(a)',err=10,iomsg=errmsg) title
     write(message, '(a,a)' ) '- ',trim(title)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

!    Read and write more data describing psp parameters
     read (tmp_unit,*,err=10,iomsg=errmsg) znucl,zion,pspdat
     write(message,'(a,f9.5,f10.5,2x,i8,t47,a)')'-',znucl,zion,pspdat,'znucl, zion, pspdat'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

     read (tmp_unit,*,err=10,iomsg=errmsg) pspcod,pspxc,lmax,lloc,mmax,r2well
     if(pspxc<0) then
       write(message, '(i5,i8,2i5,i10,f10.5,t47,a)' ) &
&       pspcod,pspxc,lmax,lloc,mmax,r2well,'pspcod,pspxc,lmax,lloc,mmax,r2well'
     else
       write(message, '(4i5,i10,f10.5,t47,a)' ) &
&       pspcod,pspxc,lmax,lloc,mmax,r2well,'pspcod,pspxc,lmax,lloc,mmax,r2well'
     end if
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

   else if (usexml == 1 .and. xmlpaw == 0) then

! the following is probably useless - already read in everything in inpspheads
#if defined HAVE_PSML
!     write(message,'(a,a)') &
!&     '- pspatm: Reading pseudopotential header in XML form from ', trim(psps%filpsp(ipsp))
!     call wrtout(ab_out,message,'COLL')
!     call wrtout(std_out,  message,'COLL')

     call psxml2abheader( psps%filpsp(ipsp), psphead, atmsymb, creator, 0 )
     znucl = psphead%znuclpsp
     zion = psphead%zionpsp
     pspdat = psphead%pspdat
     pspcod = psphead%pspcod
     pspxc =  psphead%pspxc
     lmax = psphead%lmax
     !lloc   = 0 ! does this mean s? in psml case the local potential can be different from any l channel
     lloc = -1
     mmax = -1
     r2well = 0

     write(message,'(a,1x,a3,3x,a)') "-",atmsymb,trim(creator)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     write(message,'(a,f9.5,f10.5,2x,i8,t47,a)')'-',znucl,zion,pspdat,'znucl, zion, pspdat'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     if(pspxc<0) then
       write(message, '(i5,i8,2i5,i10,f10.5,t47,a)' ) &
&       pspcod,pspxc,lmax,lloc,mmax,r2well,'pspcod,pspxc,lmax,lloc,mmax,r2well'
     else
       write(message, '(4i5,i10,f10.5,t47,a)' ) &
&       pspcod,pspxc,lmax,lloc,mmax,r2well,'pspcod,pspxc,lmax,lloc,mmax,r2well'
     end if
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
#else
     write(message,'(a,a)')  &
&     'ABINIT is not compiled with XML support for reading this type of pseudopotential ', &
&     trim(psps%filpsp(ipsp))
     MSG_BUG(message)
#endif
! END useless
   else if (usexml == 1 .and. xmlpaw == 1) then
     write(message,'(a,a)')  &
&     '- pspatm : Reading pseudopotential header in XML form from ', trim(psps%filpsp(ipsp))
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

!    Return header informations
     call pawpsxml2ab(psps%filpsp(ipsp),ecut_tmp, pspheads_tmp,0)
     lmax=pspheads_tmp%lmax
     pspxc=pspheads_tmp%pspxc
     znucl=pspheads_tmp%znuclpsp
     pawpsp_header%basis_size=pspheads_tmp%pawheader%basis_size
     pawpsp_header%l_size=pspheads_tmp%pawheader%l_size
     pawpsp_header%lmn_size=pspheads_tmp%pawheader%lmn_size
     pawpsp_header%mesh_size=pspheads_tmp%pawheader%mesh_size
     pawpsp_header%pawver=pspheads_tmp%pawheader%pawver
     pawpsp_header%shape_type=pspheads_tmp%pawheader%shape_type
     pawpsp_header%rpaw=pspheads_tmp%pawheader%rpaw
     pawpsp_header%rshp=pspheads_tmp%pawheader%rshp
     lloc=0; pspcod=17

   else if (useupf == 1) then
     if (psps%usepaw /= 0) then
       MSG_ERROR("UPF format not allowed with PAW (USPP part not read yet)")
     end if

     pspcod = 11
     r2well = 0

!    should initialize znucl,zion,pspxc,lmax,lloc,mmax
     call upf2abinit (psps%filpsp(ipsp), znucl, zion, pspxc, lmax, lloc, mmax, &
&     psps, epsatm, xcccrc, indlmn, ekb, ffspl, nproj, vlspl, xccc1d)

   else
     MSG_ERROR("You should not be here! erroneous type or pseudopotential input")
   end if

!  ------------------------------------------------------------------------------
!  Check data for consistency against main routine input

!  Does required spin-orbit characteristics agree with format
!  TODO: in case of pspcod 5 (phoney) and 8 (oncvpsp) this is not specific enough.
!  they can be non-SOC as well.
!  HGH is ok - can always turn SOC on or off.
!  PAW is ok - can be used with or without SOC
!  write(std_out,*) pspso
   if((pspcod/=3).and.(pspcod/=5).and.(pspcod/=8).and.(pspcod/=10).and. &
&     (pspcod/=7).and.(pspcod/=17))then
!    If pspso requires internal characteristics, set it to 1 for non-HGH psps
     if(psps%pspso(ipsp)==1) psps%pspso(ipsp)=0
     if(psps%pspso(ipsp)/=0)then
       write(message, '(3a,i0,3a)' )&
&       'Pseudopotential file cannot give spin-orbit characteristics,',ch10,&
&       'while pspso(itypat)= ',psps%pspso(ipsp),'.',ch10,&
&       'Action: check your pseudopotential and input files for consistency.'
       MSG_ERROR(message)
     end if
   end if

!  Does nuclear charge znuclpsp agree with psp input znucl
   !write(std_out,*)znucl,ipsp,psps%znuclpsp(ipsp)
   !MGNAG: v5[66] gives NAG in %znuclpsp if -nan
   if (abs(psps%znuclpsp(ipsp)-znucl)>tol8) then
     write(message, '(a,f10.5,2a,f10.5,5a)' )&
&     'Pseudopotential file znucl=',znucl,ch10,&
&     'does not equal input znuclpsp=',psps%znuclpsp(ipsp),' better than 1e-08 .',ch10,&
&     'znucl is read from the psp file in pspatm, while',ch10,&
&     'znuclpsp is read in iofn2.'
     MSG_BUG(message)
   end if

!  Is the highest angular momentum within limits?
!  Recall mpsang is 1+highest l for nonlocal correction.
!  Nonlocal corrections for s, p, d, and f are supported.
   if (lmax+1>psps%mpsang) then
     write(message, '(a,i0,a,i0,a,a)' )&
&     'input lmax+1= ',lmax+1,' exceeds mpsang= ',psps%mpsang,ch10,&
&     'indicates input lmax too large for dimensions.'
     MSG_BUG(message)
   end if

!  Check several choices for ixc against pspxc
!  ixc is from ABINIT code; pspxc is from atomic psp file
   if (dtset%ixc==0) then
     MSG_WARNING('Note that input ixc=0 => no xc is being used.')
   else if(dtset%ixc/=pspxc) then
     write(message, '(a,i0,3a,i0,9a)' )&
&     'Pseudopotential file pspxc= ',pspxc,',',ch10,&
&     'not equal to input ixc= ',dtset%ixc,'.',ch10,&
&     'These parameters must agree to get the same xc ',ch10,&
&     'in ABINIT code as in psp construction.',ch10,&
&     'Action: check psp design or input file.',ch10,&
&     'Assume experienced user. Execution will continue.'
     MSG_WARNING(message)
   end if

   if (lloc>lmax .and. pspcod/=4 .and. pspcod/=8 .and. pspcod/=10) then
     write(message, '(a,2i12,a,a,a,a)' )&
&     'lloc,lmax=',lloc,lmax,ch10,&
&     'chosen l of local psp exceeds range from input data.',ch10,&
&     'Action: check pseudopotential input file.'
     MSG_ERROR(message)
   end if

!  Does the pspcod agree with type of calculation (paw or not)?
   if (((pspcod/=7.and.pspcod/=17).and.psps%usepaw==1).or.((pspcod==7.or.pspcod==17).and.psps%usepaw==0)) then
     write(message, '(a,i2,a,a,i0,a)' )&
&     'In reading atomic psp file, finds pspcod= ',pspcod,ch10,&
&     'This is not an allowed value with usepaw= ',psps%usepaw,'.'
     MSG_BUG(message)
   end if

   if (.not.psps%vlspl_recipSpace .and. &
&   (pspcod /= 2 .and. pspcod /= 3 .and. pspcod /= 10 .and. pspcod /= 7)) then
!    The following "if" statement can substitute the one just before once libBigDFT
!    has been upgraded to include pspcod 10
!    if (.not.psps%vlspl_recipSpace .and. (pspcod /= 2 .and. pspcod /= 3 .and. pspcod /= 10)) then
     write(message, '(a,i2,a,a)' )&
&     'In reading atomic psp file, finds pspcod=',pspcod,ch10,&
&     'This is not an allowed value with real space computation.'
     MSG_BUG(message)
   end if

!  MJV 16/6/2009 added pspcod 11 for upf format
   if( pspcod<1 .or. (pspcod>11.and.pspcod/=17) ) then
     write(message, '(a,i0,4a)' )&
&     'In reading atomic psp file, finds pspcod= ',pspcod,ch10,&
&     'This is not an allowed value. Allowed values are 1 to 11 .',ch10,&
&     'Action: check pseudopotential input file.'
     MSG_ERROR(message)
   end if

!  -----------------------------------------------------------------------
!  Set various terms to 0 in case not defined below
   ABI_ALLOCATE(e990,(psps%mpssoang))
   ABI_ALLOCATE(e999,(psps%mpssoang))
   ABI_ALLOCATE(rcpsp,(psps%mpssoang))
   ABI_ALLOCATE(rms,(psps%mpssoang))
   ABI_ALLOCATE(epspsp,(psps%mpssoang))
   ABI_ALLOCATE(ekb1,(psps%mpssoang))
   ABI_ALLOCATE(ekb2,(psps%mpssoang))
   e990(:)=zero ;e999(:)=zero
   rcpsp(:)=zero;rms(:)=zero
   ekb1(:)=zero ;ekb2(:)=zero
   epspsp(:)=zero
   qchrg=0

!  ----------------------------------------------------------------------
   if(pspcod==1 .or. pspcod==4)then

!    Teter pseudopotential (pspcod=1 or 4)
     call psp1in(dq,ekb,ekb1,ekb2,epsatm,epspsp,&
&     e990,e999,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,&
&     mmax,psps%mpsang,psps%mqgrid_ff,nproj,psps%n1xccc,pspcod,qchrg,psps%qgrid_ff,&
&     rcpsp,rms,psps%useylm,vlspl,xcccrc,xccc1d,zion,psps%znuclpsp(ipsp))

   else if (pspcod==2)then

!    GTH pseudopotential
     call psp2in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps,vlspl,dvlspl,zion)
     xccc1d(:,:)=0.0d0 ; qchrg=0.0d0 ; xcccrc=0.0d0

   else if (pspcod==3)then

!    HGH pseudopotential
     call psp3in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps, psps%pspso(ipsp), &
&     vlspl,zion)
     xccc1d(:,:)=0.0d0 ; qchrg=0.0d0 ; xcccrc=0.0d0

   else if (pspcod==5)then

!    Old phoney pseudopotentials
     call psp5in(ekb,ekb1,ekb2,epsatm,epspsp,&
&     e990,e999,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,&
&     mmax,psps%mpsang,psps%mpssoang,psps%mqgrid_ff,nproj,psps%n1xccc,psps%pspso(ipsp),qchrg,psps%qgrid_ff,&
&     rcpsp,rms,psps%useylm,vlspl,xcccrc,xccc1d,zion,psps%znuclpsp(ipsp))

   else if (pspcod==6)then

!    FHI pseudopotentials
     call psp6in(ekb,epsatm,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,mmax,&
&     psps%mpsang,psps%mqgrid_ff,nproj,psps%n1xccc,psps%optnlxccc,psps%positron,qchrg,psps%qgrid_ff,psps%useylm,vlspl,&
&     xcccrc,xccc1d,zion,psps%znuclpsp(ipsp))

   else if (pspcod==7)then
!    PAW "pseudopotentials"
     call pawpsp_7in(epsatm,ffspl,dtset%icoulomb,dtset%ixc,&
&     lmax,psps%lnmax,mmax,psps%mqgrid_ff,psps%mqgrid_vl,&
&     pawrad,pawtab,dtset%pawxcdev,psps%qgrid_ff,psps%qgrid_vl,&
&     dtset%usewvl,dtset%usexcnhat_orig,vlspl,xcccrc,dtset%xclevel,&
&     dtset%xc_denpos,zion,psps%znuclpsp(ipsp))

   else if (pspcod==8)then

!    DRH pseudopotentials
     call psp8in(ekb,epsatm,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,mmax,&
&     psps%mpsang,psps%mpssoang,psps%mqgrid_ff,psps%mqgrid_vl,nproj,psps%n1xccc,psps%pspso(ipsp),&
&     qchrg,psps%qgrid_ff,psps%qgrid_vl,psps%useylm,vlspl,xcccrc,xccc1d,zion,psps%znuclpsp(ipsp),nctab,maxrad)

#if defined DEV_YP_DEBUG_PSP
     call psp_dump_outputs("DBG",pspcod,psps%lmnmax,psps%lnmax,psps%mpssoang, &
&     psps%mqgrid_ff,psps%n1xccc,mmax,maxrad,epsatm,qchrg,xcccrc,nctab, &
&     indlmn,nproj,ekb,ffspl,vlspl,xccc1d)
#endif

   else if (pspcod==9)then

#if defined HAVE_PSML
     call psp9in(psps%filpsp(ipsp),ekb,epsatm,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,mmax,&
&     psps%mpsang,psps%mpssoang,psps%mqgrid_ff,psps%mqgrid_vl,nproj,psps%n1xccc, &
&     psps%pspso(ipsp),qchrg,psps%qgrid_ff,psps%qgrid_vl,psps%useylm,vlspl,&
&     xcccrc,xccc1d,zion,psps%znuclpsp(ipsp),nctab,maxrad)

#if defined DEV_YP_DEBUG_PSP
     call psp_dump_outputs("DBG",pspcod,psps%lmnmax,psps%lnmax,psps%mpssoang, &
&     psps%mqgrid_ff,psps%n1xccc,mmax,maxrad,epsatm,qchrg,xcccrc,nctab, &
&     indlmn,nproj,ekb,ffspl,vlspl,xccc1d)
#endif
#else
     write(message,'(2a)')  &
&     'ABINIT is not compiled with XML support for reading this type of pseudopotential ', &
&     trim(psps%filpsp(ipsp))
     MSG_BUG(message)
#endif

   else if (pspcod==10)then

!    HGH pseudopotential, full h/k matrix read
     call psp10in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps, psps%pspso(ipsp), &
&     vlspl,zion)
     xccc1d(:,:)=0.0d0 ; qchrg=0.0d0 ; xcccrc=0.0d0

!    NB for pspcod 11 the reading has already been done above.
   else if (pspcod==17)then
!    PAW XML "pseudopotentials"
     call pawpsp_17in(epsatm,ffspl,dtset%icoulomb,ipsp,dtset%ixc,lmax,&
&     psps%lnmax,mmax,psps%mqgrid_ff,psps%mqgrid_vl,pawpsp_header,pawrad,pawtab,&
&     dtset%pawxcdev,psps%qgrid_ff,psps%qgrid_vl,dtset%usewvl,&
&     dtset%usexcnhat_orig,vlspl,xcccrc,&
&     dtset%xclevel,dtset%xc_denpos,pspheads_tmp%zionpsp,psps%znuclpsp(ipsp))
     call paw_setup_free(paw_setuploc)
   end if

   close (unit=tmp_unit)

!  ----------------------------------------------------------------------
   if (pspcod==2 .or. pspcod==3 .or. pspcod==10)then
     write(message, '(a,a,a,a,a,a,a,a,a,a)' )ch10,&
&     ' pspatm : COMMENT -',ch10,&
&     '  the projectors are not normalized,',ch10,&
&     '  so that the KB energies are not consistent with ',ch10,&
&     '  definition in PRB44, 8503 (1991). ',ch10,& ! [[cite:Gonze1991]]
&     '  However, this does not influence the results obtained hereafter.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
!    The following lines are added to keep backward compatibilty
     maxrad=zero
#if defined HAVE_BIGDFT
     do ii=1,size(psps%gth_params%psppar,1)-1 ! psppar first dim begins at 0
       if (psps%gth_params%psppar(ii,0,ipsp)/=zero) maxrad=max(maxrad,psps%gth_params%psppar(ii,0,ipsp))
     end do
     if (abs(maxrad)<=tol12) then
       psps%gth_params%radii_cf(ipsp,3)=zero
     else
       psps%gth_params%radii_cf(ipsp,3)=max( &
&       min(dtset%wvl_crmult*psps%gth_params%radii_cf(ipsp,1),15._dp*maxrad)/dtset%wvl_frmult, &
&       psps%gth_params%radii_cf(ipsp,2))
     end if
#endif
   end if

   if (pspcod/=7.and.pspcod/=17) then
     write(message, '(a,f14.8,a,a)' ) '  pspatm : epsatm=',epsatm,ch10,&
&     '         --- l  ekb(1:nproj) -->'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     iln0=0
     do ilmn=1,psps%lmnmax
       iln=indlmn(5,ilmn)
       if (iln>iln0) then
         il=indlmn(1,ilmn)
         if (indlmn(6,ilmn)==1) then
           iln0=iln0+nproj(il+1)
           !if (dtset%optdriver == RUNL_SIGMA) then
           !  do ii=0,nproj(il+1)-1
           !    ekb(iln+ii) = zero
           !  end do
           !end if
           write(message, '(13x,i1,4f12.6)' ) il,(ekb(iln+ii),ii=0,nproj(il+1)-1)
         else
           iln0=iln0+nproj(il+psps%mpsang)
           write(message, '(2x,a,i1,4f12.6)' ) 'spin-orbit ',il,(ekb(iln+ii),ii=0,nproj(il+psps%mpsang)-1)
         end if
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
       end if
     end do
   end if

   ! NC: Evalute spline-fit of the model core charge in reciprocal space.
   ! TODO: Be careful, because we will be using the PAW part in which tcore is always avaiable!
   ! Should add a test with 2 NC pseudos: one with NLCC and the other without!
   if (psps%usepaw == 0) then
     call nctab_eval_tcorespl(nctab, psps%n1xccc, xcccrc, xccc1d, psps%mqgrid_vl, psps%qgrid_vl)
   end if

   write(message,'(3a)') ' pspatm: atomic psp has been read ',' and splines computed',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   ABI_DEALLOCATE(e990)
   ABI_DEALLOCATE(e999)
   ABI_DEALLOCATE(rcpsp)
   ABI_DEALLOCATE(rms)
   ABI_DEALLOCATE(ekb1)
   ABI_DEALLOCATE(ekb2)
   ABI_DEALLOCATE(epspsp)
   ABI_DEALLOCATE(nproj)

   if (dtset%prtvol > 9 .and. psps%usepaw==0 .and. psps%lmnmax>3) then

     write (filnam, '(a,i0,a)') trim(dtfil%fnameabo_pspdata), ipsp, ".dat"
     if (open_file(filnam, message, newunit=unt) /= 0) then
       MSG_ERROR(message)
     end if
     write (unt,*) '# Pseudopotential data in reciprocal space as used by ABINIT'
     write (unt,'(a)', ADVANCE='NO') '# index       vlocal   '
     if (psps%lnmax > 0) &
&     write (unt,'(a,I3)', ADVANCE='NO')   '           1st proj(l=', indlmn(1,1)
     if (psps%lnmax > 1) &
&     write (unt,'(a,I3)', ADVANCE='NO')   ')            2nd(l=', indlmn(1,2)
     if (psps%lnmax > 2) &
&     write (unt,'(a,I3,a)', ADVANCE='NO') ')            3rd(l=', indlmn(1,3), ')'
     write (unt,*)

     do ii = 1, psps%mqgrid_vl
       write(unt, '(I5,E24.16)', ADVANCE='NO') ii, vlspl(ii,1)
       if (psps%lnmax > 0) write(unt, '(E24.16)', ADVANCE='NO') ffspl(ii,1,1)
       if (psps%lnmax > 1) write(unt, '(E24.16)', ADVANCE='NO') ffspl(ii,1,2)
       if (psps%lnmax > 2) write(unt, '(E24.16)', ADVANCE='NO') ffspl(ii,1,3)
       write(unt, *)
     end do
     close(unt)

     write (filnam, '(a,i0,a)') trim(dtfil%fnameabo_nlcc_derivs), ipsp, ".dat"
     if (open_file(filnam, message, newunit=unt) /= 0) then
       MSG_ERROR(message)
     end if
     write (unt,*) '# Non-linear core corrections'
     write (unt,*) '#  r, pseudocharge, 1st, 2nd, 3rd, 4th, 5th derivatives'
     do ii = 1, psps%n1xccc
       write (unt,*) xcccrc*(ii-1)/(psps%n1xccc-1), xccc1d(ii,1), xccc1d(ii,2), &
&       xccc1d(ii,3), xccc1d(ii,4), &
&       xccc1d(ii,5), xccc1d(ii,6)
     end do
     close(unt)
   end if

 end if !me=0

 if (paral_mode==1) then
   call timab(48,1,tsec)
   call pawpsp_bcast(comm_mpi,epsatm,ffspl,pawrad,pawtab,vlspl,xcccrc)
   call timab(48,2,tsec)
 end if

 if (psps%usepaw==1) then
   indlmn(:,:)=0
   indlmn(1:6,1:pawtab%lmn_size)=pawtab%indlmn(1:6,1:pawtab%lmn_size)
 end if

!--------------------------------------------------------------------
!WVL+PAW:
 if (dtset%usepaw==1 .and. (dtset%icoulomb /= 0 .or. dtset%usewvl==1)) then
#if defined HAVE_BIGDFT
   psps%gth_params%psppar(:,:,ipsp) = UNINITIALIZED(1._dp)
   psps%gth_params%radii_cf(ipsp,:) = UNINITIALIZED(1._dp)
   call wvl_descr_psp_fill(psps%gth_params, ipsp, psps%pspxc(1), int(psps%zionpsp(ipsp)), int(psps%znuclpsp(ipsp)), 0)
#endif

!  The following lines are added to keep backward compatibilty
   maxrad=zero
#if defined HAVE_BIGDFT
   do ii=1,size(psps%gth_params%psppar,1)-1 ! psppar first dim begins at 0
     if (psps%gth_params%psppar(ii,0,ipsp)/=zero) maxrad=max(maxrad,psps%gth_params%psppar(ii,0,ipsp))
   end do
   if (abs(maxrad)<=tol12) then
!== MT COMMENT
!    Damien wants to activate this (in order to directly compare to bigDFT):
     psps%gth_params%radii_cf(ipsp,3)= psps%gth_params%radii_cf(ipsp,2)
!    But, this changes strongly file references.
!    So, I keep this, waiting for Tonatiuh s validation
     psps%gth_params%radii_cf(ipsp,3)= &
&     (psps%gth_params%radii_cf(ipsp,1)+psps%gth_params%radii_cf(ipsp,2))*half
!== MT COMMENT
   else
     psps%gth_params%radii_cf(ipsp,3)=max( &
&     min(dtset%wvl_crmult*psps%gth_params%radii_cf(ipsp,1),15._dp*maxrad)/dtset%wvl_frmult, &
&     psps%gth_params%radii_cf(ipsp,2))
   end if
   if(present(comm_mpi)) then
     call pawpsp_wvl(psps%filpsp(ipsp),pawrad,pawtab,dtset%usewvl,dtset%wvl_ngauss,comm_mpi)
   else
     call pawpsp_wvl(psps%filpsp(ipsp),pawrad,pawtab,dtset%usewvl,dtset%wvl_ngauss)
   end if
#endif
 end if

!end of WVL+PAW section
!----------------------------------------------------

 return

 ! Handle IO error
 10 continue
 MSG_ERROR(errmsg)

end subroutine pspatm
!!***

!!****f* ABINIT/psp_dump_outputs
!! NAME
!! psp_dump_outputs
!!
!! FUNCTION
!! (To be described ...)
!!
!! COPYRIGHT
!! Copyright (C) 2017-2019 ABINIT group (YP)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt
!!
!! INPUTS
!! (to be filled)
!!
!! OUTPUT
!! (to be filled)
!!
!! SIDE EFFECTS
!! (to be filled)
!!
!! PARENTS
!!      pspatm
!!
!! CHILDREN
!!
!! SOURCE

subroutine psp_dump_outputs(pfx,pspcod,lmnmax,lnmax,mpssoang, &
&      mqgrid,n1xccc,mmax,maxrad,epsatm,qchrg,xcccrc,nctab, &
&      indlmn,nproj,ekb,ffspl,vlspl,xccc1d)

 use defs_basis
 use m_errors

 use defs_datatypes, only : nctab_t

!Arguments ------------------------------------
!scalars
 character(len=*), intent(in) :: pfx
 integer,intent(in) :: pspcod,lmnmax,lnmax,mpssoang,mqgrid,n1xccc
 integer,intent(in) :: mmax
 real(dp),intent(in) :: maxrad,epsatm,qchrg,xcccrc
 type(nctab_t),intent(in) :: nctab
!arrays
 integer,intent(in) :: indlmn(6,lmnmax),nproj(mpssoang)
 real(dp),intent(in) :: ekb(lnmax),ffspl(mqgrid,2,lnmax),vlspl(mqgrid,2)
 real(dp),intent(in) :: xccc1d(n1xccc,6)

!Local variables ------------------------------
!scalars
 integer, parameter :: dump = 64
 integer :: ierr, i, j ,k
 character(len=500) :: msg

 ! *********************************************************************

 open(unit=dump, file=trim(pfx)//"_psp_info.yaml", status='REPLACE', err=10, iostat=ierr)

 write(dump,'(3a)') "%YAML 1.2", ch10, "---"

 write(dump, '(2a)') ch10, "# Pseudopotential info"
 write(dump, '(a,1x,i8)') "pspcod:", pspcod

 write(dump, '(2a)') ch10, "# Array dimensions"
 write(dump, '(a)') "dims:"
 write(dump, '(4x,a,1x,i8)') "lmnmax:", lmnmax
 write(dump, '(4x,a,1x,i8)') "lnmax:", lnmax
 write(dump, '(4x,a,1x,i8)') "mpssoang:", mpssoang
 write(dump, '(4x,a,1x,i8)') "mqgrid:", mqgrid
 write(dump, '(4x,a,1x,i8)') "n1xccc:", n1xccc
 write(dump, '(4x,a,1x,i8)') "mmax:", mmax

 write(dump, '(2a)') ch10, "# Quantities"
 write(dump, '(a,1x,e12.5)') "maxrad:", maxrad
 write(dump, '(a,1x,e12.5)') "epsatm:", epsatm
 write(dump, '(a,1x,e12.5)') "qchrg:", qchrg
 write(dump, '(a,1x,e12.5)') "xcccrc:", xcccrc

 write(dump, '(2a)') ch10, "# Structure: nctab"
 write(dump, '(a)') "nctab:"
 write(dump,'(4x,a,":",1x,i4)') "mqgrid_vl", nctab%mqgrid_vl
 write(dump,'(4x,a,":",1x,l4)') "has_tvale", nctab%has_tvale
 write(dump,'(4x,a,":",1x,l4)') "has_tcore", nctab%has_tcore
 write(dump,'(4x,a,":",1x,e12.5)') "dncdq0", nctab%dncdq0
 write(dump,'(4x,a,":",1x,e12.5)') "d2ncdq0", nctab%d2ncdq0
 write(dump,'(4x,a,":",1x,e12.5)') "dnvdq0", nctab%dnvdq0

 if ( nctab%has_tvale ) then
   write(dump, '(2a)') ch10, "# Array: nctab_tvalespl(mqgrid_vl,2)"
   write(dump, '(a)') "nctab_tvalespl:"
   do j=1,2
     do i=1,nctab%mqgrid_vl
       if ( i == 1 ) then
         write(dump,'(4x,a,1x,e12.5)') "- -", nctab%tvalespl(i,j)
       else
         write(dump,'(4x,a,1x,e12.5)') "  -", nctab%tvalespl(i,j)
       end if
     end do
   end do
 end if

 if ( nctab%has_tcore ) then
   write(dump, '(2a)') ch10, "# Array: nctab_tcorespl(mqgrid_vl,2)"
   write(dump, '(a)') "nctab_tcorespl:"
   do j=1,2
     do i=1,nctab%mqgrid_vl
       if ( i == 1 ) then
         write(dump,'(4x,a,1x,e12.5)') "- -", nctab%tcorespl(i,j)
       else
         write(dump,'(4x,a,1x,e12.5)') "  -", nctab%tcorespl(i,j)
       end if
     end do
   end do
 end if

 write(dump, '(2a)') ch10, "# Array: integer indlmn(6,lmnmax)"
 write(dump, '(a)') "indlmn:"
 do i=1,lmnmax
   write(dump,'(4x,a,i4,5(",",i4),a)') "- [", indlmn(:,i), "]"
 end do

 write(dump, '(2a)') ch10, "# Array: integer nproj(mpssoang)"
 write(dump, '(a)') "nproj:"
 do i=1,mpssoang
   write(dump,'(4x,"-",1x,i4)') nproj(i)
 end do

 write(dump, '(2a)') ch10, "# Array: double ekb(lnmax)"
 write(dump, '(a)') "ekb:"
 do i=1,lnmax
   write(dump,'(4x,"-",1x,e12.5)') ekb(i)
 end do

 write(dump, '(2a)') ch10, "# Array: ffspl(mqgrid,2,lnmax)"
 write(dump, '(a)') "ffspl:"
 do k=1,lnmax
   do j=1,2
     do i=1,mqgrid
       if ( (i == 1) .and. (j == 1) ) then
         write(dump,'(4x,a,1x,e12.5)') "- - -", ffspl(i,j,k)
       else if ( i == 1 ) then
         write(dump,'(4x,a,1x,e12.5)') "  - -", ffspl(i,j,k)
       else
         write(dump,'(4x,a,1x,e12.5)') "    -", ffspl(i,j,k)
       end if
     end do
   end do
 end do

 write(dump, '(2a)') ch10, "# Array: vlspl(mqgrid,2)"
 write(dump, '(a)') "vlspl:"
 do j=1,2
   do i=1,mqgrid
     if ( i == 1 ) then
       write(dump,'(4x,a,1x,e12.5)') "- -", vlspl(i,j)
     else
       write(dump,'(4x,a,1x,e12.5)') "  -", vlspl(i,j)
     end if
   end do
 end do

 write(dump, '(2a)') ch10, "# Array: xccc1d(n1xccc,6)"
 write(dump, '(a)') "xccc1d:"
 do j=1,6
   do i=1,mqgrid
     if ( i == 1 ) then
       write(dump,'(4x,a,1x,e12.5)') "- -", xccc1d(i,j)
     else
       write(dump,'(4x,a,1x,e12.5)') "  -", xccc1d(i,j)
     end if
   end do
 end do

 write (dump,'(2a)') ch10, "..."

 close(dump)

 return
 10 continue

 if ( ierr /= 0 ) then
   write(msg,'(a,a,a,i8)') "Error writing pseudopotential information", ch10, "IOSTAT=", ierr
   MSG_WARNING(msg)
 end if

end subroutine psp_dump_outputs
!!***

end module m_pspini
!!***
