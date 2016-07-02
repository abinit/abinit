!{\src2tex{textfont=tt}}
!!****f* ABINIT/newfermie1
!! NAME
!! newfermie1
!!
!! FUNCTION
!! This routine computes the derivative of the fermi energy wrt
!! the active perturbation for use in evaluating the edocc term
!! and active subspace contribution to the first-order wavefunctions
!! in the case of metals.  This is presently used only for the
!! strain perturbation, and only for Q = 0.
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  fe1fixed=fixed contribution to the first-order Fermi energy
!!  ipert=index of perturbation
!!  istep=index of the number of steps in the routine scfcv
!!  ixc= choice of exchange-correlation scheme
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nfftot= total number of FFT grid points
!!  nhatfermi(nfft,nspden)=fermi-level compensation charge density (PAW only)
!!  nspden=number of spin-density components
!!  ntypat=number of atom types
!!  occopt=option for occupancies
!!  paw_an(natom) <type(paw_an_type)>=paw arrays for 0th-order quantities given on angular mesh
!!  paw_an1(natom) <type(paw_an_type)>=paw arrays for 1st-order quantities given on angular mesh
!!  paw_ij1(natom) <type(paw_ij_type)>=(1st-order) paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawnzlm=-- PAW only -- option for the computation of non-zero
!!          lm moments of the on-sites densities
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij1(natom) <type(pawrhoij_type)>= paw rhoij 1st-order occupancies
!!  pawrhoijfermi(natom) <type(pawrhoij_type)>=paw rhoij occupancies at Fermi level
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  xclevel= XC functional level
!!  prtvol=control print volume and debugging output
!!  rhorfermi(nfft,nspden)=fermi-level electronic density
!!  ucvol=unit cell volume in bohr**3
!!  usepaw=1 if PAW is activated
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vtrial1(cplex*nfft,nspden)=1-st order potential
!!  vxc1(cplex*nfft,nspden)=1-st order XC potential
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  fermie1=derivative of fermi energy wrt perturbation
!!   at input  : old value
!!   at output : updated value
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      dotprod_vn,free_my_atmtab,get_my_atmtab,pawdfptenergy,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine newfermie1(cplex,fermie1,fe1fixed,ipert,istep,ixc,my_natom,natom,nfft,nfftot,&
&                     nhatfermi,nspden,ntypat,occopt,paw_an,paw_an1,paw_ij1,pawang,pawnzlm,pawrad,&
&                     pawrhoij1,pawrhoijfermi,pawtab,pawxcdev,prtvol,rhorfermi,&
&                     ucvol,usepaw,usexcnhat,vtrial1,vxc1,xclevel,&
&                     mpi_atmtab,comm_atom) ! optional arguments (parallelism)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi, only : xmpi_comm_self

 use m_pawang,     only : pawang_type
 use m_pawrad,     only : pawrad_type
 use m_pawtab,     only : pawtab_type
 use m_paw_an,     only : paw_an_type
 use m_paw_ij,     only : paw_ij_type
 use m_pawrhoij,   only : pawrhoij_type
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'newfermie1'
 use interfaces_14_hidewrite
 use interfaces_53_spacepar
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,ipert,istep,ixc,my_natom,natom,nfft,nfftot,nspden,ntypat
 integer,intent(in) :: occopt,pawnzlm,pawxcdev,prtvol,usepaw,usexcnhat,xclevel
 integer,optional,intent(in) :: comm_atom
 real(dp),intent(in) :: fe1fixed,ucvol
 real(dp),intent(inout) :: fermie1
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: rhorfermi(nfft,nspden),vtrial1(cplex*nfft,nspden)
 real(dp),intent(in) :: nhatfermi(:,:),vxc1(:,:)
 type(paw_an_type),intent(in) :: paw_an(my_natom*usepaw)
 type(paw_an_type),intent(inout) :: paw_an1(my_natom*usepaw)
 type(paw_ij_type),intent(inout) :: paw_ij1(my_natom*usepaw)
 type(pawrad_type),intent(in) :: pawrad(ntypat*usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij1(my_natom*usepaw),pawrhoijfermi(my_natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: ipert0,my_comm_atom,nzlmopt,nzlmopt_fermi,option,pawprtvol
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: doti,fe1_scf,fe1_tmp,fermie1_new,fermie1rs
 character(len=500) :: msg
!arrays
 integer, pointer :: my_atmtab(:)
 real(dp) :: fe1_paw(2)
 real(dp), allocatable :: rhor_nonhat(:,:),vtrial1_novxc(:,:)

! *********************************************************************

!Tests
 if (cplex==2) then
   msg='Not compatible with cplex=2!'
   MSG_BUG(msg)
 end if
 if (usepaw==1.and.usexcnhat==0.and.(size(nhatfermi)<=0.or.size(vxc1)<=0)) then
   msg='Should have nhatfermi and vxc1 allocated with usexcnhat=0!'
   MSG_BUG(msg)
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 if(occopt>=3 .and. occopt <=8) then

!  The product of the current trial potential and the so-called Fermi level
!  density is integrated to give the local potential contributions to the
!  first-order Fermi level.
   option=1
   if (usepaw==1.and.usexcnhat==0) then
     ABI_ALLOCATE(rhor_nonhat,(nfft,nspden))
     ABI_ALLOCATE(vtrial1_novxc,(nfft,nspden))
     rhor_nonhat(1:nfft,1:nspden)=rhorfermi(1:nfft,1:nspden)-nhatfermi(1:nfft,1:nspden)
     vtrial1_novxc(1:nfft,1:nspden)=vtrial1(1:nfft,1:nspden)-vxc1(1:nfft,1:nspden)
     call dotprod_vn(cplex,rhor_nonhat,fe1_scf,doti,nfft,nfftot,&
&     nspden,option,vtrial1,ucvol)
     call dotprod_vn(cplex,nhatfermi,fe1_tmp,doti,nfft,nfftot,&
&     nspden,option,vtrial1_novxc,ucvol)
     fe1_scf=fe1_scf+fe1_tmp
     ABI_DEALLOCATE(rhor_nonhat)
     ABI_DEALLOCATE(vtrial1_novxc)
   else
     call dotprod_vn(cplex,rhorfermi,fe1_scf,doti,nfft,nfftot,&
&     nspden,option,vtrial1,ucvol)
   end if

   fe1_paw(:)=zero
!  PAW on-site contribution (use Fermi level occupation matrix)
   if (usepaw==1) then
     ipert0=0;pawprtvol=0
     nzlmopt=0;if (istep>1) nzlmopt=pawnzlm
     if (istep==1.and.pawnzlm>0) nzlmopt=-1
     nzlmopt_fermi=0;if (pawnzlm>0) nzlmopt_fermi=-1
     call pawdfptenergy(fe1_paw,ipert,ipert0,ixc,my_natom,natom,ntypat,nzlmopt,&
&     nzlmopt_fermi,paw_an,paw_an1,paw_ij1,pawang,pawprtvol,pawrad,&
&     pawrhoij1,pawrhoijfermi,pawtab,pawxcdev,xclevel,&
&     mpi_atmtab=my_atmtab, comm_atom=my_comm_atom)
   end if

!  The fixed contributions consisting of non-local potential and kinetic terms
!  are added
   fermie1_new=fe1fixed+fe1_scf+fe1_paw(1)
   fermie1rs=(fermie1-fermie1_new)**2
   fermie1=fermie1_new

   if(prtvol>=10)then
     write(msg, '(a,i5,2es18.8)' ) ' fermie1, residual squared',istep,fermie1,fermie1rs
     call wrtout(std_out,msg,'COLL')
   end if

 else
   fermie1=zero
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine newfermie1
!!***
