!{\src2tex{textfont=tt}}
!!****f* ABINIT/mklocl
!! NAME
!! mklocl
!!
!! FUNCTION
!! This method is a wrapper for mklocl_recipspace and mklocl_realspace.
!! It does some consistency checks before calling one of the two methods.
!!
!! Optionally compute :
!!  option=1 : local ionic potential throughout unit cell
!!  option=2 : contribution of local ionic potential to E gradient wrt xred
!!  option=3 : contribution of local ionic potential to
!!                stress tensor (only with reciprocal space computations)
!!  option=4 : contribution of local ionic potential to
!!                second derivative of E wrt xred  (only with reciprocal space computations)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  if(option==3) eei=local pseudopotential part of total energy (hartree)
!!  gmet(3,3)=reciprocal space metric ($\textrm{Bohr}^{-2}$).
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff on $|G|^2$: see setup1 for definition (doubled sphere).
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in unit cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms.
!!  option= (see above)
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qprtrb(3)= integer wavevector of possible perturbing potential
!!   in basis of reciprocal lattice translations
!!  rhog(2,nfft)=electron density rho(G) (electrons/$\textrm{Bohr}^3$)
!!    (needed if option==2 or if option==4)
!!  rhor(nfft,nspden)=electron density in electrons/bohr**3.
!!    (needed if option==2 or if option==4)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol=unit cell volume ($\textrm{Bohr}^3$).
!!  vprtrb(2)=complex amplitude of possible perturbing potential; if nonzero,
!!   perturbing potential is added of the form
!!   $V(G)=(vprtrb(1)+I*vprtrb(2))/2$ at the values G=qprtrb and
!!   $(vprtrb(1)-I*vprtrb(2))/2$ at $G=-qprtrb$ (integers)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (if option==1) vpsp(nfft)=local crystal pseudopotential in real space.
!!  (if option==2) grtn(3,natom)=grads of Etot wrt tn.
!!  (if option==3) lpsstr(6)=components of local psp part of stress tensor
!!   (Cartesian coordinates, symmetric tensor) in hartree/$\textrm{bohr}^3$
!!   Store 6 unique components in order 11, 22, 33, 32, 31, 21
!!  (if option==4) dyfrlo(3,3,natom)=d2 Eei/d tn(i)/d tn(j).  (Hartrees)
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Note that the present routine is tightly connected to the dfpt_vlocal.f routine,
!! that compute the derivative of the local ionic potential
!! with respect to one atomic displacement. The argument list
!! and the internal loops to be considered were sufficiently different
!! as to make the two routine different.
!!
!! PARENTS
!!      forces,prcref,prcref_PMA,respfn,setvtr
!!
!! CHILDREN
!!      mklocl_realspace,mklocl_recipspace,mklocl_wavelets,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mklocl(dtset, dyfrlo,eei,gmet,gprimd,grtn,gsqcut,lpsstr,mgfft,&
&  mpi_enreg,natom,nattyp,nfft,ngfft,nspden,ntypat,option,pawtab,ph1d,psps,qprtrb,&
&  rhog,rhor,rprimd,ucvol,vprtrb,vpsp,wvl,wvl_den,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_errors

 use m_pawtab, only : pawtab_type

#if defined HAVE_BIGDFT
 use BigDFT_API, only : ELECTRONIC_DENSITY
 use m_abi2big, only : wvl_rho_abi2big
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mklocl'
 use interfaces_41_geometry
 use interfaces_67_common, except_this_one => mklocl
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,natom,nfft,nspden,ntypat,option
 real(dp),intent(in) :: eei,gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
 type(wvl_denspot_type), intent(inout) :: wvl_den
 type(pawtab_type),intent(in)  :: pawtab(ntypat*psps%usepaw)
!arrays
 integer,intent(in) :: nattyp(ntypat),ngfft(18),qprtrb(3)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rhog(2,nfft),rprimd(3,3)
 real(dp),intent(in) :: vprtrb(2),xred(3,natom)
 real(dp),intent(in),target :: rhor(nfft,nspden)
 real(dp),intent(out) :: dyfrlo(3,3,natom),grtn(3,natom),lpsstr(6)
 real(dp),intent(inout) :: vpsp(nfft)

!Local variables-------------------------------
!scalars
 character(len=500) :: message
!arrays
 real(dp),allocatable :: xcart(:,:)
 real(dp),pointer :: rhor_ptr(:,:)

! *************************************************************************

 if (option < 1 .or. option > 4) then
   write(message,'(a,i0,a,a)')&
&   'From the calling routine, option=',option,ch10,&
&   'The only allowed values are between 1 and 4.'
   MSG_ERROR(message)
 end if
 if (option > 2 .and. .not.psps%vlspl_recipSpace) then
   write(message,'(a,i0,a,a,a,a)')&
&   'From the calling routine, option=',option,ch10,&
&   'but the local part of the pseudo-potential is in real space.',ch10,&
&   'Action: set icoulomb = 0 to turn-off real space computations.'
   MSG_ERROR(message)
 end if
 if (option > 2 .and. dtset%usewvl == 1) then
   write(message,'(a,i0,a,a)')&
&   'From the calling routine, option=',option,ch10,&
&   'but this is not implemented yet from wavelets.'
   MSG_ERROR(message)
 end if

 if (dtset%usewvl == 0) then
!  Plane wave case
   if (psps%vlspl_recipSpace) then
     call mklocl_recipspace(dyfrlo,eei,gmet,gprimd,grtn,gsqcut,lpsstr,mgfft, &
&     mpi_enreg,psps%mqgrid_vl,natom,nattyp,nfft,ngfft, &
&     ntypat,option,dtset%paral_kgb,ph1d,psps%qgrid_vl,qprtrb,rhog,ucvol, &
&     psps%vlspl,vprtrb,vpsp)
   else
     call mklocl_realspace(grtn,dtset%icoulomb,mpi_enreg,natom,nattyp,nfft, &
&     ngfft,dtset%nscforder,nspden,ntypat,option,pawtab,psps,rhog,rhor, &
&     rprimd,dtset%typat,ucvol,dtset%usewvl,vpsp,xred)
   end if
 else
!  Store xcart for each atom
   ABI_ALLOCATE(xcart,(3, dtset%natom))
   call xred2xcart(dtset%natom, rprimd, xcart, xred)
!  Eventually retrieve density
#if defined HAVE_BIGDFT
   if (option>1.and.wvl_den%denspot%rhov_is/=ELECTRONIC_DENSITY) then
     rhor_ptr => rhor ! Just to bypass intent(inout)
     call wvl_rho_abi2big(1,rhor_ptr,wvl_den)
   end if
#endif
!  Wavelets case
   call mklocl_wavelets(dtset%efield, grtn, mpi_enreg, dtset%natom, &
&   nfft, nspden, option, rprimd, vpsp, &
&   wvl_den, wvl, xcart)
   ABI_DEALLOCATE(xcart)
 end if

end subroutine mklocl
!!***
