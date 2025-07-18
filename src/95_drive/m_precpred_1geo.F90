!!****m* ABINIT/m_precpred_1geo
!! NAME
!!  m_precpred_1geo
!!
!! FUNCTION
!! Single geometry: apply force and stress preconditioner followed by geometry predictor.
!! Choose among the whole set of geometry predictors defined by iommov.
!!
!! COPYRIGHT
!!  Copyright (C) 2018-2025 ABINIT group (DCA, XG, GMR, SE)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_precpred_1geo

 use defs_basis
 use m_errors
 use m_abicore
 use m_abimover
 use m_abihist
 use m_xmpi
 use m_nctk
 use m_pimd
 use netcdf
#if defined HAVE_LOTF
 use lotfpath
 use m_pred_lotf
#endif

 use m_fstrings,           only : strcat, sjoin, itoa
 use m_geometry,           only : chkdilatmx
 use m_crystal,            only : crystal_t
 use m_pred_bfgs,          only : pred_bfgs, pred_lbfgs
 use m_pred_delocint,      only : pred_delocint
 use m_pred_fire,          only : pred_fire
 use m_pred_isokinetic,    only : pred_isokinetic
 use m_pred_diisrelax,     only : pred_diisrelax
 use m_pred_nose,          only : pred_nose
 use m_pred_srkhna14,      only : pred_srkna14
 use m_pred_isothermal,    only : pred_isothermal
 use m_pred_verlet,        only : pred_verlet
 use m_pred_velverlet,     only : pred_velverlet
 use m_pred_moldyn,        only : pred_moldyn
 use m_pred_langevin,      only : pred_langevin
 use m_pred_langevin_pimd, only : pred_langevin_pimd
 use m_pred_steepdesc,     only : pred_steepdesc
 use m_pred_simple,        only : pred_simple, prec_simple
 use m_pred_hmc,           only : pred_hmc
 use m_ipi,                only : ipi_pred
!use m_generate_training_set, only : generate_training_set

 implicit none

 private
!!***

 public :: precpred_1geo
!!***

contains
!!***

!!****f* ABINIT/precpred_1geo
!! NAME
!! mover
!!
!! FUNCTION
!! Single geometry: apply force and stress preconditioner followed by geometry predictor.
!! Choose among the whole set of geometry predictors defined by ionmov.
!!
!! INPUTS
!!  (TO BE DESCRIBED)
!!
!! OUTPUT
!!  (TO BE DESCRIBED)
!!
!! SIDE EFFECTS
!!  hist is the main quantity updated, contains xred, rprimd and acell.
!!  (TO BE DESCRIBED)
!!
!! NOTES
!!
!! SOURCE

subroutine precpred_1geo(ab_mover,ab_xfh,amu_curr,deloc,dt_chkdilatmx,comm_cell,dilatmx,filnam_ds4,hist,hmctt,&
& icycle,iexit,itime,mttk_vars,nctime,ncycle,nerr_dilatmx,npsp,ntime,pimd_param,rprimd_orig,skipcycle,usewvl)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: comm_cell,dt_chkdilatmx,hmctt,icycle,itime,nctime,npsp,ntime
 integer, intent(inout) :: iexit,ncycle,nerr_dilatmx
 integer, intent(in) :: usewvl
 real(dp), intent(in) :: dilatmx
 logical, intent(inout) :: skipcycle
 character(len=fnlen), intent(in) :: filnam_ds4
 type(ab_xfh_type),intent(inout) :: ab_xfh
 type(abihist), intent(inout) :: hist
 type(abimover), intent(in) :: ab_mover
 type(delocint), intent(inout) :: deloc
 type(mttk_type), intent(inout) :: mttk_vars
 type(pimd_type), intent(in) :: pimd_param
!arrays
 real(dp), intent(in) :: amu_curr(ab_mover%ntypat)
 real(dp), intent(in) :: rprimd_orig(3,3)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: ii,me,nloop
!integer :: iatom
 logical,parameter :: DEBUG=.FALSE.
!character(len=500) :: message
 character(len=500) :: dilatmx_errmsg
 character(len=fnlen) :: filename
 type(abiforstr) :: preconforstr ! Preconditioned forces and stress
 type(crystal_t) :: crystal
!arrays
 real(dp) :: acell(3),rprimd(3,3)
 real(dp), allocatable :: xred(:,:)

! ***************************************************************

!DEBUG
!write(std_out,'(a,i4)')' m_precpred_1geo, enter : ab_mover%ionmov=',ab_mover%ionmov
!ABI_MALLOC(xred,(3,ab_mover%natom))
!call hist2var(acell,hist,ab_mover%natom,rprimd,xred,DEBUG)
!do iatom=1,ab_mover%natom
!  write(std_out,'(i4,3es14.6)') iatom,xred(1:3,iatom)
!enddo
!ABI_FREE(xred)
!ENDDEBUG

 me=xmpi_comm_rank(comm_cell)

!Precondition forces, stress and energy
 call abiforstr_ini(preconforstr,ab_mover%natom)

 if (ab_mover%goprecon>0 .or. iexit==1)then
   call prec_simple(ab_mover,preconforstr,hist,icycle,itime,iexit)
 end if

!Call to each predictor
!MT->GAF: dirty trick to predict vel(t)
!do a double loop: 1- compute vel, 2- exit
 nloop=1

 if (nctime>0.and.iexit==1) then
   iexit=0;nloop=2
 end if

 do ii=1,nloop
   if (ii==2) iexit=1

   select case (ab_mover%ionmov)
   case (1)
     call pred_moldyn(ab_mover,hist,icycle,itime,ncycle,ntime,DEBUG,iexit)
   case (2,3)
     call pred_bfgs(ab_mover,ab_xfh,preconforstr,hist,ab_mover%ionmov,itime,DEBUG,iexit)
   case (4,5)
     call pred_simple(ab_mover,hist,iexit)
   case (6,7)
     call pred_verlet(ab_mover,hist,ab_mover%ionmov,itime,ntime,DEBUG,iexit)
   case (8)
     call pred_nose(ab_mover,hist,itime,ntime,DEBUG,iexit)
   case (9)
     call pred_langevin(ab_mover,hist,icycle,itime,ncycle,ntime,DEBUG,iexit,skipcycle)
   case (10,11)
     call pred_delocint(ab_mover,ab_xfh,deloc,preconforstr,hist,ab_mover%ionmov,itime,DEBUG,iexit)
   case (12)
     call pred_isokinetic(ab_mover,hist,itime,ntime,DEBUG,iexit)
   case (13)
     call pred_isothermal(ab_mover,hist,itime,mttk_vars,ntime,DEBUG,iexit)
   case (14)
     call pred_srkna14(ab_mover,hist,icycle,DEBUG,iexit,skipcycle)
   case (15)
     call pred_fire(ab_mover, ab_xfh,preconforstr,hist,ab_mover%ionmov,itime,DEBUG,iexit)
   case (16)
     call pred_langevin_pimd(ab_mover,hist,itime,DEBUG,pimd_param)
   case (20)
     call pred_diisrelax(ab_mover,hist,itime,ntime,DEBUG,iexit)
   case (21)
     call pred_steepdesc(ab_mover,preconforstr,hist,itime,DEBUG,iexit)
   case (22)
     call pred_lbfgs(ab_mover,ab_xfh,preconforstr,hist,ab_mover%ionmov,itime,DEBUG,iexit)
#if defined HAVE_LOTF
   case (23)
     call pred_lotf(ab_mover,hist,itime,icycle,DEBUG,iexit)
#endif
   case (24)
     call pred_velverlet(ab_mover,hist,itime,ntime,DEBUG,iexit)
   case (25)
     call pred_hmc(ab_mover,hist,itime,icycle,ntime,hmctt,mttk_vars,DEBUG,iexit)
   case (27)
     !In case of ionmov 27, all the atomic configurations have been computed at the
     !beginning of the routine in generate_training_set, thus we just need to increase the indexes
     !in the hist
     hist%ihist = abihist_findIndex(hist,+1)
   case (28)
     call ipi_pred(ab_mover, hist, itime, ntime, DEBUG, iexit, comm_cell)
   case default
     ABI_ERROR(sjoin("Wrong value of ionmov:", itoa(ab_mover%ionmov)))
   end select

 end do

 ABI_MALLOC(xred,(3,ab_mover%natom))
 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,DEBUG)

 ! check dilatmx here and correct if necessary
 if (usewvl == 0) then
   call chkdilatmx(dt_chkdilatmx,dilatmx,rprimd,rprimd_orig,dilatmx_errmsg)
   if (LEN_TRIM(dilatmx_errmsg) /= 0) then
     ABI_WARNING(dilatmx_errmsg)
     nerr_dilatmx = nerr_dilatmx+1
     if (nerr_dilatmx > 3) then
       ! Write last structure before aborting, so that we can restart from it.
       ! zion is not available, but it's not useful here.
       if (me == master) then
         ! Init crystal
         hist%ihist = abihist_findIndex(hist,-1)
         call hist2var(acell,hist,ab_mover%natom,rprimd,xred,DEBUG)
         call crystal%init(amu_curr,0,ab_mover%natom,&
           npsp,ab_mover%ntypat,ab_mover%nsym,rprimd,ab_mover%typat,xred,&
           [(-one, ii=1,ab_mover%ntypat)],ab_mover%znucl,2,.False.,.False.,"dilatmx_structure",&
           symrel=ab_mover%symrel,tnons=ab_mover%tnons,symafm=ab_mover%symafm)

         ! Write netcdf file
         filename = strcat(filnam_ds4, "_DILATMX_STRUCT.nc")
         NCF_CHECK(crystal%ncwrite_path(filename))
         call crystal%free()
       end if
       call xmpi_barrier(comm_cell)
       write (dilatmx_errmsg, '(a,i0,9a)') &
        'Dilatmx has been exceeded too many times (', nerr_dilatmx, ')',ch10, &
        'See the description of dilatmx and chkdilatmx input variables.',ch10, &
        'Action: either first do a calculation with chkdilatmx=0, or ',ch10,&
        'restart your calculation with a larger dilatmx, or larger lattice vectors.',ch10,&
        'Warning: With chkdilatmx = 0 the final computation of lattice parameters might be inaccurate.'
       ABI_ERROR_CLASS(dilatmx_errmsg, "DilatmxError")
     end if
   else
     nerr_dilatmx=0
   end if
 end if

!DEBUG
!write(std_out,'(a,i4)')' m_precpred_1geo, exit '
!do iatom=1,ab_mover%natom
!  write(std_out,'(i4,3es14.6)') iatom,xred(1:3,iatom)
!enddo
!ENDDEBUG

 call abiforstr_fin(preconforstr)
 ABI_FREE(xred)

end subroutine precpred_1geo
!!***

end module m_precpred_1geo
!!***
