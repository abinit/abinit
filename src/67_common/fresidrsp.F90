!{\src2tex{textfont=tt}}
!!****f* ABINIT/fresidrsp
!!
!! NAME
!! fresidrsp
!!
!! FUNCTION
!! Compute the forces due to the residual of the potential (or density)
!! in RECIPROCAL SPACE, using
!!  - the atomic density read in psp file (PAW or NC with nctval_spl e.g. psp8 format)
!!  - a gaussian atomic density (norm-conserving psps if nctval_spl is not available)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! atindx1(natom)=index table for atoms, inverse of atindx
!! dtset <type(dataset_type)>=all input variables in this dataset
!!  | densty(ntypat,4)=parameters for initialisation of the density of each atom type
!!  | icoulomb=0 periodic treatment of Hartree potential, 1 use of Poisson solver
!!  | ixc= choice of exchange-correlation scheme
!!  | natom=number of atoms in cell.
!!  | nspden=number of spin-density components
!!  | typat(natom)=integer type for each atom in cell
!! gmet(3,3)=reciprocal space metric
!! gprimd(3,3)=reciprocal space dimensional primitive translations
!! gsqcut=cutoff value on G**2 for sphere inside fft box
!! mgfft=maximum size of 1D FFTs
!! mpi_enreg=information about MPI parallelization
!! mqgrid=number of grid pts in q array for atomic density spline n^AT(q)
!! nattyp(ntypat)=number of atoms of each type in cell
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! ntypat=number of types of atoms in cell.
!! psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase information for given atom coordinates.
!! qgrid(mqgrid)=q grid for spline atomic valence density n^AT(q) from 0 to qmax
!! ucvol=unit cell volume (bohr**3).
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! vresid(nfft,nspden)=potential residual - (non-collinear magn. : only V11 and V22 are used)
!! zion(ntypat)=charge on each type of atom (real number)
!! znucl(ntypat)=atomic number, for each type of atom
!!
!! OUTPUT
!! gresid(3,natom)=forces due to the residual of the potential
!!
!! PARENTS
!!      forces
!!
!! CHILDREN
!!      atm2fft,fourdp,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fresidrsp(atindx1,dtset,gmet,gprimd,gresid,gsqcut,mgfft,mpi_enreg,mqgrid,nattyp,nfft,&
&          ngfft,ntypat,psps,pawtab,ph1d,qgrid,ucvol,usepaw,vresid,zion,znucl)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi

 use defs_datatypes, only : pseudopotential_type
 use m_atomdata,     only : atom_length
 use m_pawtab,       only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fresidrsp'
 use interfaces_14_hidewrite
 use interfaces_53_ffts
 use interfaces_64_psp
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mqgrid,nfft,ntypat,usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(pseudopotential_type),intent(in) :: psps
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: atindx1(dtset%natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(in) :: qgrid(mqgrid),vresid(nfft,dtset%nspden),zion(ntypat)
 real(dp),intent(in) :: znucl(ntypat)
 real(dp),intent(out) :: gresid(3,dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: itypat,optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv
 logical :: usegauss
!arrays
 integer :: dummy3(3)
 real(dp) :: dummy2(2)
 real(dp) :: dummy_in1(0),dummy_in2(0) 
 real(dp) :: dummy_out1(0),dummy_out2(0),dummy_out3(0),dummy_out4(0),dummy_out5(0),dummy_out6(0) 
 real(dp) :: strn_dummy6(6),strv_dummy6(6)
 real(dp),allocatable :: gauss(:,:),vresg(:,:),work(:)

! *************************************************************************

!Inits
 optatm=0;optdyfr=0;opteltfr=0;optgr=1;optstr=0;optv=0;optn=1
 ABI_ALLOCATE(vresg,(2,nfft))

!Transfer potential residual to reciprocal space
!Use only Vres=Vres11+Vres22=Vres_up+Vres_dn
 ABI_ALLOCATE(work,(nfft))
 work(:)=vresid(:,1)
 if (dtset%nspden>=2) work(:)=work(:)+vresid(:,2)
 call fourdp(1,vresg,work,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
 ABI_DEALLOCATE(work)

!Determine wether a gaussan atomic density has to be used or not
 usegauss=.true.
 if (usepaw==0) usegauss = any(.not.psps%nctab(1:ntypat)%has_tvale)
 if (usepaw==1) usegauss=(minval(pawtab(1:ntypat)%has_tvale)==0)
 if (usegauss) then
   optn2=3
   ABI_ALLOCATE(gauss,(2,ntypat))
   do itypat=1,ntypat
     gauss(1,itypat)=zion(itypat)
     gauss(2,itypat) = atom_length(dtset%densty(itypat,1),zion(itypat),znucl(itypat))
   end do
   call wrtout(std_out," Computing residual forces using gaussian functions as atomic densities", "COLL")
 else
   optn2=2
   ABI_ALLOCATE(gauss,(2,0))
   call wrtout(std_out," Computing residual forces using atomic densities taken from pseudos", "COLL")
 end if

!Compute forces due to residual
 call atm2fft(atindx1,dummy_out1,dummy_out2,dummy_out3,dummy_out4,&
& dummy_out5,gauss,gmet,gprimd,gresid,dummy_out6,gsqcut,mgfft,&
& mqgrid,dtset%natom,nattyp,nfft,ngfft,ntypat,optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,&
& psps,pawtab,ph1d,qgrid,dummy3,dummy_in1,strn_dummy6,strv_dummy6,ucvol,usepaw,vresg,vresg,vresg,dummy2,dummy_in2,&
& comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
& paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)

!In case of nspden>=2, has to apply 1/2 factor
 if (dtset%nspden>=2) gresid=gresid*half

 ABI_DEALLOCATE(gauss)
 ABI_DEALLOCATE(vresg)

end subroutine fresidrsp
!!***
