!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_denspot_set
!! NAME
!!  wvl_denspot_set
!!
!! FUNCTION
!!  Fill in denspot datatype with information related
!!  to density and potential data.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2017 ABINIT group (DC)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      allocaterhopot,density_descriptors,dpbox_set
!!      initialize_dft_local_fields,wrtout,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_denspot_set(den,gth_params,ixc,natom,nsppol,rprimd,wvl,&
&                          wvl_crmult,wvl_frmult,wvl_mpi_comm,xred)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

#if defined HAVE_BIGDFT
 use BigDFT_API,only: initialize_DFT_local_fields,allocateRhoPot, &
&                     input_variables,dpbox_set,density_descriptors
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_denspot_set'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
  integer,intent(in):: ixc,natom,nsppol,wvl_mpi_comm
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(in) :: wvl_frmult,wvl_crmult
  real(dp), intent(inout)  :: xred(3,natom)
  type(wvl_denspot_type), intent(out) :: den
  type(wvl_internal_type),intent(in)  :: wvl
  type(pseudopotential_gth_type),intent(in)::gth_params

!Local variables-------------------------------
#if defined HAVE_BIGDFT
  integer :: groupsize,me,nproc
  real(dp), allocatable :: xcart(:,:)
  character(len=3),parameter :: rho_commun='DBL'
  character(len=500) :: message
  character(len=4) :: SICapproach
  type(local_zone_descriptors) :: Lzd
#endif

  ! *************************************************************************
 
!DEBUG
!write (std_out,*) ' wvl_denspot_set : enter'
!ENDDEBUG
 
#if defined HAVE_BIGDFT

 write(message, '(a,a)' ) ch10,&
& ' wvl_denspot_set: Create wavelet type denspot.'
 call wrtout(std_out,message,'COLL')

 nproc=xmpi_comm_size(wvl_mpi_comm)
 me=xmpi_comm_rank(wvl_mpi_comm)
 groupsize=0

!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, natom))
 call xred2xcart(natom, rprimd, xcart, xred)

 call initialize_DFT_local_fields(den%denspot, ixc, nsppol)

!number of planes for the density
!dpbox%nscatterarr(jproc, 1) = ngfft3_density
!number of planes for the potential
!dpbox%nscatterarr(jproc, 2) = ngfft3_potential
!starting offset for the potential
!dpbox%nscatterarr(jproc, 3) = density_start + potential_shift - 1
!GGA XC shift between density and potential
!dpbox%nscatterarr(jproc, 4) = potential_shift

 SICapproach="NONE"
 Lzd%hgrids(1:3)=wvl%h(1:3)
 Lzd%Glr%d%n1i=wvl%Glr%d%n1i
 Lzd%Glr%d%n2i=wvl%Glr%d%n2i
 Lzd%Glr%d%n3i=wvl%Glr%d%n3i
 call dpbox_set(den%denspot%dpbox,Lzd,den%denspot%xc,me,nproc,wvl_mpi_comm,groupsize,&
& SICapproach,wvl%atoms%astruct%geocode,nsppol)

!here dpbox can be put as input
 call density_descriptors(me,nproc,den%denspot%xc,nsppol,wvl_crmult,wvl_frmult,wvl%atoms,&
 den%denspot%dpbox,rho_commun,xcart,gth_params%radii_cf,den%denspot%rhod)

!Note: change allocateRhoPot
 call allocateRhoPot(wvl%Glr,nsppol,wvl%atoms,xcart,den%denspot)

!Aditional informations.
 den%symObj = wvl%atoms%astruct%sym%symObj

 ABI_DEALLOCATE(xcart)

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) ixc,natom,nsppol,wvl_mpi_comm,rprimd(1,1),wvl_frmult,wvl_crmult,&
& xred(1,1),den%symObj,wvl%h(1),gth_params%psppar
#endif

!DEBUG
!write (std_out,*) ' wvl_denspot_set : exit'
!ENDDEBUG

end subroutine wvl_denspot_set
!!***
