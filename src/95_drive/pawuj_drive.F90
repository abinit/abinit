!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawuj_drive
!! NAME
!! pawuj_drive
!!
!! FUNCTION
!!  Drive for automatic determination of U
!!  Relevant only in PAW+U context
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DJA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cpus= cpu time limit in seconds
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs for the "coarse" grid (see NOTES below)
!!   | mkmem =number of k points treated by this node.
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in cell.
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |      for the "coarse" grid (see NOTES below)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  ecore=core psp energy (part of total energy) (hartree)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)= # atoms of each type.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points and spins
!!
!! SIDE EFFECTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=updated wavefunctions.
!!  dtefield <type(efield_type)> = variables related to Berry phase calculations (see initberry.f)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  initialized= if 0 the initialization of the gstate run is not yet finished
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!       for the "fine" grid (see NOTES below)
!!  occ(mband*nkpt*nsppol)=occupation number for each band (often 2) at each k point
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!   (should be made a pure output quantity)
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in el./bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  wffnew,wffnow=struct info for wf disk files.
!!  wvl <type(wvl_data)>=all wavelets data.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)= at input, previous reduced dimensionless atomic coordinates
!!                     at output, current xred is transferred to xred_old
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      pawuj_det,pawuj_free,pawuj_ini,scfcv_run
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawuj_drive(scfcv, &
& dtset,electronpositron,rhog,rhor,rprimd,&
& xred,xred_old)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_electronpositron, only : electronpositron_type
 use m_scfcv

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawuj_drive'
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(scfcv_t), intent(inout) :: scfcv
 type(dataset_type),intent(inout) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 !type(wffile_type),intent(inout) :: wffnew,wffnow
!arrays
 real(dp), intent(inout) :: rprimd(3,3)
 real(dp), pointer :: rhog(:,:),rhor(:,:)
 real(dp), intent(inout) :: xred(3,dtset%natom),xred_old(3,dtset%natom)

!Local variables -------------------------
!scalars
 integer,target :: ndtpawuj=4
 integer :: iuj,conv_retcode
 real(dp) :: ures
 character(len=500) :: message
!arrays
 real(dp),allocatable :: cgstart(:,:)
 type(macro_uj_type),allocatable,target :: dtpawuj(:)
! *********************************************************************

 DBG_ENTER("COLL")

 if (dtset%macro_uj==0) then
   message = 'Macro_uj must be set !'
   MSG_BUG(message)
 end if

 ABI_DATATYPE_ALLOCATE(dtpawuj,(0:ndtpawuj))
 ABI_ALLOCATE(cgstart,(2,scfcv%mcg))

!DEBUG
!write(std_out,*)'pawuj_drive: before ini dtpawuj(:)%iuj ', dtpawuj(:)%iuj
!END DEBUG
 call pawuj_ini(dtpawuj,ndtpawuj)

 cgstart=scfcv%cg
 do iuj=1,ndtpawuj
!  allocate(dtpawuj(iuj)%rprimd(3,3)) ! this has already been done in pawuj_ini
   dtpawuj(iuj)%macro_uj=dtset%macro_uj
   dtpawuj(iuj)%pawprtvol=dtset%pawprtvol
   dtpawuj(iuj)%diemix=dtset%diemix
   dtpawuj(iuj)%pawujat=dtset%pawujat
   dtpawuj(iuj)%nspden=dtset%nspden
   dtpawuj(iuj)%rprimd=dtset%rprimd_orig(1:3,1:3,1)
 end do

!allocate(dtpawuj(0)%vsh(0,0),dtpawuj(0)%occ(0,0))

 do iuj=1,2
   if (iuj>1) scfcv%cg(:,:)=cgstart(:,:)

!  DEBUG
!  write(std_out,*)'drive_pawuj before count dtpawuj(:)%iuj ', dtpawuj(:)%iuj
!  END DEBUG

   dtpawuj(iuj*2-1)%iuj=iuj*2-1

   scfcv%ndtpawuj=>ndtpawuj
   scfcv%dtpawuj=>dtpawuj

   !call scfcv_new(ab_scfcv_in,ab_scfcv_inout,dtset,electronpositron,&
!&   paw_dmft,rhog,rhor,rprimd,wffnew,wffnow,xred,xred_old,conv_retcode)
   call scfcv_run(scfcv,electronpositron,rhog,rhor,rprimd,xred,xred_old,conv_retcode)

   scfcv%fatvshift=scfcv%fatvshift*(-one)

 end do

!Calculate Hubbard U (or J)
 call pawuj_det(dtpawuj,ndtpawuj,trim(scfcv%dtfil%filnam_ds(4))//"_UJDET.nc",ures)
 dtset%upawu(dtset%typat(dtset%pawujat),1)=ures/Ha_eV

!Deallocations
 do iuj=0,ndtpawuj
   call pawuj_free(dtpawuj(iuj))
 end do

 ABI_DATATYPE_DEALLOCATE(dtpawuj)
 ABI_DEALLOCATE(cgstart)

 DBG_EXIT("COLL")

end subroutine pawuj_drive
!!***
