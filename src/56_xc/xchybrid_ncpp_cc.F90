!{\src2tex{textfont=tt}}
!!****f* ABINIT/xchybrid_ncpp_cc
!! NAME
!! xchybrid_ncpp_cc
!!
!! FUNCTION
!! XC Hybrid Norm-Conserving PseudoPotential Core Correction:
!! Relevant only for Norm-Conserving PseudoPotentials (NCPP) and for hybrid functionals.
!! Compute the correction to the XC energy/potential due to the lack of core wave-functions:
!! As Fock exchange cannot be computed for core-core and core-valence interactions, these
!! contribution have to be also substracted from the GGA exchange-correlation.
!!
!! COPYRIGHT
!!  Copyright (C) 2015-2017 ABINIT group (FA,MT,FJ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>= all input variables in this dataset
!!  mpi_enreg= information about MPI parallelization
!!  nfft= number of fft grid points.
!!  ngfft(1:3)= integer fft box dimensions, see getng for ngfft(4:8).
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  optstr= calculate corrected vxc if optstr=1 
!!  rhor(nfft,nspden)= electron density in real space in electrons/bohr**3
!!  rprimd(3,3)= dimensional primitive translations for real space in Bohr.
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc*(1-usepaw),6,ntypat)=1D core charge function and five derivatives,
!!                          for each type of atom, from psp (used in Norm-conserving only)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!! OUTPUT
!!
!! SIDE EFFECTS
!!  enxc= exchange correlation energy
!!  grxc= correction yo the forces 
!!  strsxc(6)= exchange correlation contribution to stress tensor
!!  vxc= exchange correlation potential
!!  vxcavg= unit cell average of Vxc
!!
!! NOTES
!!  The final expression of the XC potential (that should be added to alpha*VFock[Psi_val]) is:
!!   Vxc=Vx[rho_core+rho_val] - alpha*Vx[rho_val] + Vc[rho_core+rho_val]
!!  To accomodate libXC convention, Vxc is computed as follows:
!!   Vxc=Vxc_libXC[rho_val] + Vxc_gga[rho_core+rho_val] - Vxc_gga[rho_val]
!!
!! PARENTS
!!      forces,forstr,rhotov,setvtr
!!
!! CHILDREN
!!      dtset_copy,dtset_free,libxc_functionals_end,libxc_functionals_init
!!      metric,mkcore,rhotoxc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine xchybrid_ncpp_cc(dtset,enxc,mpi_enreg,nfft,ngfft,n3xccc,rhor,rprimd,strsxc,vxcavg,xccc3d,vxc,grxc,xcccrc,xccc1d,&
&                           xred,n1xccc,optstr)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xcdata
 use libxc_functionals

 use defs_abitypes, only : MPI_type, dataset_type
 use m_dtset,       only : dtset_copy, dtset_free

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xchybrid_ncpp_cc'
 use interfaces_41_geometry
 use interfaces_56_xc, except_this_one => xchybrid_ncpp_cc
!End of the abilint section

 implicit none

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: nfft,n3xccc
 integer,optional,intent(in) :: n1xccc,optstr
 real(dp),intent(out) :: enxc,vxcavg
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,dtset%nspden),rprimd(3,3),xccc3d(n3xccc)
 real(dp),intent(out) :: strsxc(6)
 real(dp),optional,intent(in) :: xcccrc(dtset%ntypat),xred(3,dtset%natom),xccc1d(:,:,:)
 real(dp),optional,intent(out) :: grxc(3,dtset%natom),vxc(nfft,dtset%nspden)

!Local variables -------------------------------------------------------
!scalars
 integer :: ixc_gga,ndim,nkxc,n3xccc_null,option,optstr_loc,usexcnhat
 real(dp) :: enxc_corr,ucvol,vxcavg_corr
 character(len=500) :: msg
 type(xcdata_type) :: xcdata_gga,xcdata_hybrid 
 logical :: calcgrxc
!arrays
 integer :: gga_id(2)
 real(dp) :: nhat(1,0),nhatgr(1,1,0),strsxc_corr(6),gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: kxc_dum(:,:),vxc_corr(:,:),xccc3d_null(:),dyfrx2_dum(:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!Not relevant for PAW
 if (dtset%usepaw==1) return
 if(n3xccc==0) return
 calcgrxc=(present(grxc).and.present(n1xccc).and.present(xcccrc).and.present(xred).and.present(xccc1d)) 
 optstr_loc=0
 if(present(optstr)) optstr_loc=1
!Not applicable for electron-positron
 if (dtset%positron>0) then
   msg='NCPP+Hybrid functionals not applicable for electron-positron calculations!'
   MSG_ERROR(msg)
 end if

!Select the GGA on which the hybrid functional is based on
!or return if not applicable
 if (dtset%ixc==41.or.dtset%ixc==42) then
   ixc_gga = 11
 else if (dtset%ixc<0) then
   if (libxc_functionals_gga_from_hybrid(gga_id=gga_id)) then
     ixc_gga=-gga_id(1)*1000-gga_id(2)
   else
     return
   end if
 else
   return
 end if

!Define xcdata_hybrid as well as xcdata_gga
 call xcdata_init(dtset%intxc,dtset%ixc,&
&    dtset%nelect,dtset%tphysel,dtset%usekden,dtset%vdw_xc,dtset%xc_tb09_c,dtset%xc_denpos,xcdata_hybrid)
 call xcdata_init(dtset%intxc,ixc_gga,&
&    dtset%nelect,dtset%tphysel,dtset%usekden,dtset%vdw_xc,dtset%xc_tb09_c,dtset%xc_denpos,xcdata_gga)

 nkxc=0;ndim=0;usexcnhat=0;n3xccc_null=0
 ABI_ALLOCATE(kxc_dum,(nfft,nkxc))
 ABI_ALLOCATE(vxc_corr,(nfft,dtset%nspden))

 if (present(vxc).and.optstr_loc==0) then
!Initialize args for rhotoxc
   option=0 ! XC only
   ABI_ALLOCATE(xccc3d_null,(n3xccc_null))
!  Compute Vxc^Hybrid(rho_val)
   call rhotoxc(enxc,kxc_dum,mpi_enreg,nfft,ngfft,nhat,ndim,nhatgr,ndim,nkxc,nkxc,&
&   dtset%nspden,n3xccc_null,option,dtset%paral_kgb,rhor,rprimd,strsxc,usexcnhat,vxc,vxcavg,xccc3d_null,xcdata_hybrid)

!  Initialize GGA functional
   if (ixc_gga<0) then
     call libxc_functionals_end()
     call libxc_functionals_init(ixc_gga,dtset%nspden)
   end if

!Add Vxc^GGA(rho_core+rho_val)
   call rhotoxc(enxc_corr,kxc_dum,mpi_enreg,nfft,ngfft,nhat,ndim,nhatgr,ndim,nkxc,nkxc,&
&   dtset%nspden,n3xccc,option,dtset%paral_kgb,rhor,rprimd,strsxc_corr,usexcnhat,vxc_corr,vxcavg_corr,xccc3d,xcdata_gga)
   enxc=enxc+enxc_corr
   vxc(:,:)=vxc(:,:)+vxc_corr(:,:)
   vxcavg=vxcavg+vxcavg_corr
   strsxc(:)=strsxc(:)+strsxc_corr(:)

!Substract Vxc^GGA(rho_val)
   call rhotoxc(enxc_corr,kxc_dum,mpi_enreg,nfft,ngfft,nhat,ndim,nhatgr,ndim,nkxc,nkxc,&
&   dtset%nspden,n3xccc_null,option,dtset%paral_kgb,rhor,rprimd,strsxc_corr,usexcnhat,&
&   vxc_corr,vxcavg_corr,xccc3d_null,xcdata_gga)
   enxc=enxc-enxc_corr
   vxc(:,:)=vxc(:,:)-vxc_corr(:,:)
   vxcavg=vxcavg-vxcavg_corr
   strsxc(:)=strsxc(:)-strsxc_corr(:)
   
!Release memory
   ABI_DEALLOCATE(xccc3d_null)
 end if

 if (calcgrxc) then

   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
!  Initialize GGA functional
   if (ixc_gga<0) then
     call libxc_functionals_end()
     call libxc_functionals_init(ixc_gga,dtset%nspden)
   end if
   ABI_ALLOCATE(dyfrx2_dum,(3,3,dtset%natom))
   ABI_ALLOCATE(xccc3d_null,(n3xccc))
   option=1
! calculate xccc3d in this case
   call mkcore(strsxc_corr,dyfrx2_dum,grxc,mpi_enreg,dtset%natom,nfft,dtset%nspden,dtset%ntypat,ngfft(1),n1xccc,ngfft(2),&
&   ngfft(3),option,rprimd,dtset%typat,ucvol,vxc_corr,xcccrc,xccc1d,xccc3d_null,xred)
!Add Vxc^GGA(rho_core+rho_val)
   call rhotoxc(enxc_corr,kxc_dum,mpi_enreg,nfft,ngfft,nhat,ndim,nhatgr,ndim,nkxc,nkxc,&
&   dtset%nspden,n3xccc,option,dtset%paral_kgb,&
&   rhor,rprimd,strsxc_corr,usexcnhat,vxc_corr,vxcavg_corr,xccc3d_null,xcdata_gga)
   option=2
   call mkcore(strsxc_corr,dyfrx2_dum,grxc,mpi_enreg,dtset%natom,nfft,dtset%nspden,dtset%ntypat,ngfft(1),n1xccc,ngfft(2),&
&   ngfft(3),option,rprimd,dtset%typat,ucvol,vxc_corr,xcccrc,xccc1d,xccc3d_null,xred)
   ABI_DEALLOCATE(dyfrx2_dum)
   ABI_DEALLOCATE(xccc3d_null)
 end if

 if(optstr_loc==1) then
!  Initialize GGA functional
   if (ixc_gga<0) then
     call libxc_functionals_end()
     call libxc_functionals_init(ixc_gga,dtset%nspden)
   end if
!calculate Vxc^GGA(rho_core+rho_val)
   option=0
   call rhotoxc(enxc_corr,kxc_dum,mpi_enreg,nfft,ngfft,nhat,ndim,nhatgr,ndim,nkxc,nkxc,&
&   dtset%nspden,n3xccc,option,dtset%paral_kgb,rhor,rprimd,strsxc_corr,usexcnhat,vxc,vxcavg_corr,xccc3d,xcdata_gga)
 end if 

! Revert libxc to original settings
 if (ixc_gga<0) then
   call libxc_functionals_end()
   call libxc_functionals_init(dtset%ixc,dtset%nspden)
 end if
 ABI_DEALLOCATE(vxc_corr)
 ABI_DEALLOCATE(kxc_dum)

 DBG_EXIT("COLL")

end subroutine xchybrid_ncpp_cc
!!***
