!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawmkrho
!! NAME
!! pawmkrho
!!
!! FUNCTION
!! PAW only:
!! Build total pseudo (compensated) density (\tild_rho + \hat_rho)
!! Build compensation charge density (\hat_rho)
!! Build occupation matrix (packed storage)
!!
!! COPYRIGHT
!! Copyright (C) 2010-2018 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!         1 for GS calculations
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  ipert=index of perturbation if pawrhoij is a pertubed rhoij
!!        no meaning for ground-state calculations (should be 0)
!!  idir=direction of atomic displacement (in case of atomic displ. perturb.)
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell
!!  nspden=number of spin-density components
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of types of atoms in unit cell.
!!  paral_kgb=option for (kpt,g vectors,bands) parallelism
!!  pawang <type(pawang_type)>=angular mesh discretization and related data
!!  pawang_sym <type(pawang_type)>=angular data used for symmetrization
!!                                 optional parameter only needed for RF calculations
!!  pawfgr <type(paw_fgr_type)>=fine rectangular grid parameters
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrhoij0(natom) <type(pawrhoij_type)>= GS paw rhoij occupancies and related data (used only if ipert>0)
!!                                          optional parameter only needed for RF calculations
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  qphon(3)=wavevector of the phonon (RF only)
!!  rhopsg(2,pawfgr%nfftc)= pseudo density given on the coarse grid in reciprocal space
!!  rhopsr(pawfgr%nfftc,nspden)= pseudo density given on the coarse grid in real space
!!  rprimd(3,3)=real space primitive translations.
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!                   reciprocal space primitive translations
!!  typat(natom)=type for each atom
!!  ucvol=volume of the unit cell
!!  xred(3,natom)= reduced atomic coordinates
!!
!! OUTPUT
!!  compch_fft=compensation charge inside spheres integrated over fine fft grid
!!  pawnhat(pawfgr%nfft,nspden)=compensation charge density on fine rectangular grid (optional argument)
!!  rhog(2,pawfgr%nfft)= compensated pseudo density given on the fine grid in reciprocal space
!!                       This output is optional
!!  rhor(pawfgr%nfft,nspden)= compensated pseudo density given on the fine grid in real space
!!
!! SIDE EFFECTS
!!  pawrhoij(my_natom)= PAW occupancies
!!                   At input : values at previous step  in packed storage (pawrhoij()%rhoijp)
!!                   At output: values (symmetrized)     in packed storage (pawrhoij()%rhoijp)
!!  pawrhoij_unsym(:)= unsymmetrized PAW occupancies
!!                   At input : values (unsymmetrized) in unpacked storage (pawrhoij()%rhoij_)
!!                   At output: values in unpacked storage (pawrhoij()%rhoij_) are destroyed
!!
!! NOTES
!!  pawrhoij and pawrhoij_unsym can be identical (refer to the same pawrhoij datastructure).
!!  They should be different only if pawrhoij is distributed over atomic sites
!!  (in that case pawrhoij_unsym should not be distributed over atomic sites).
!!
!! PARENTS
!!      afterscfloop,dfpt_nstpaw,dfpt_rhofermi,dfpt_vtorho,scfcv,vtorho
!!
!! CHILDREN
!!      fourdp,pawmknhat,pawrhoij_copy,pawrhoij_free,pawrhoij_free_unpacked
!!      pawrhoij_nullify,symrhoij,timab,transgrid
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawmkrho(compch_fft,cplex,gprimd,idir,indsym,ipert,mpi_enreg,&
&          my_natom,natom,nspden,nsym,ntypat,paral_kgb,pawang,pawfgr,pawfgrtab,pawprtvol,&
&          pawrhoij,pawrhoij_unsym,&
&          pawtab,qphon,rhopsg,rhopsr,rhor,rprimd,symafm,symrec,typat,ucvol,usewvl,xred,&
&          pawang_sym,pawnhat,pawrhoij0,rhog) ! optional arguments

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_pawang,   only : pawang_type
 use m_pawtab,   only : pawtab_type
 use m_pawfgrtab,only : pawfgrtab_type
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_copy, pawrhoij_free_unpacked, &
&                       pawrhoij_nullify, pawrhoij_free, symrhoij
 use m_pawfgr,   only : pawfgr_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawmkrho'
 use interfaces_18_timing
 use interfaces_53_ffts
 use interfaces_65_paw, except_this_one => pawmkrho
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,my_natom,natom,nspden,nsym,ntypat,paral_kgb,pawprtvol
 integer,intent(in) :: usewvl
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: compch_fft
 type(MPI_type),intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(pawang_type),intent(in),optional :: pawang_sym
 type(pawfgr_type),intent(in) :: pawfgr
!arrays
 integer,intent(in) :: indsym(4,nsym,natom)
 integer,intent(in) :: symafm(nsym),symrec(3,3,nsym),typat(natom)
 real(dp),intent(in) :: gprimd(3,3),qphon(3),rprimd(3,3),xred(3,natom)
 real(dp),intent(inout),target,optional :: pawnhat(cplex*pawfgr%nfft,nspden) !vz_i
 real(dp),intent(inout) :: rhor(cplex*pawfgr%nfft,nspden)
 real(dp),intent(out),optional :: rhog(2,pawfgr%nfft)
 real(dp),intent(inout) :: rhopsg(2,pawfgr%nfftc),rhopsr(cplex*pawfgr%nfftc,nspden)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawrhoij_type),intent(inout),target :: pawrhoij(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij_unsym(:)
 type(pawrhoij_type),intent(in),target,optional :: pawrhoij0(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: choice,ider,izero,option
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp) :: rhodum(0,0,0)
 real(dp),pointer :: pawnhat_ptr(:,:)
 type(pawrhoij_type),pointer :: pawrhoij_ptr(:),pawrhoij0_ptr(:)

! ***********************************************************************

 DBG_ENTER("COLL")

 call timab(556,1,tsec)

!Compatibility tests
 if (size(pawrhoij_unsym)>0) then
   if (pawrhoij_unsym(1)%use_rhoij_==0) then
     msg='  rhoij_ field must be allocated in pawrhoij_unsym !'
     MSG_BUG(msg)
   end if
 end if
 if (ipert>0.and.(.not.present(pawrhoij0))) then
   msg='  pawrhoij0 must be present when ipert>0 !'
   MSG_BUG(msg)
 end if

!Symetrize PAW occupation matrix and store it in packed storage
 call timab(557,1,tsec)
 option=1;choice=1
 if (present(pawang_sym)) then
   call symrhoij(pawrhoij,pawrhoij_unsym,choice,gprimd,indsym,ipert,natom,nsym,ntypat,&
&   option,pawang_sym,pawprtvol,pawtab,rprimd,symafm,symrec,typat,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&   qphon=qphon)
 else
   call symrhoij(pawrhoij,pawrhoij_unsym,choice,gprimd,indsym,ipert,natom,nsym,ntypat,&
&   option,pawang,pawprtvol,pawtab,rprimd,symafm,symrec,typat,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&   qphon=qphon)
 end if
 call pawrhoij_free_unpacked(pawrhoij_unsym)
 call timab(557,2,tsec)

!In somes cases (parallelism), has to distribute the PAW occupation matrix
 if (size(pawrhoij)==natom.and.(my_natom/=natom)) then
   ABI_DATATYPE_ALLOCATE(pawrhoij_ptr,(my_natom))
   call pawrhoij_nullify(pawrhoij_ptr)
   call pawrhoij_copy(pawrhoij,pawrhoij_ptr,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom, &
&   keep_cplex=.false.,keep_itypat=.false.,keep_nspden=.false.)       
 else
   pawrhoij_ptr=>pawrhoij
 end if

!Compute compensation charge density
 ider=0;izero=0
 if (present(pawnhat)) then
   pawnhat_ptr => pawnhat
 else
   ABI_ALLOCATE(pawnhat_ptr,(pawfgr%nfft,nspden))
 end if
 if (present(pawrhoij0)) then
   pawrhoij0_ptr => pawrhoij0
 else
   pawrhoij0_ptr => pawrhoij_ptr
 end if

 call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,my_natom,natom,&
& pawfgr%nfft,pawfgr%ngfft,ider,nspden,ntypat,pawang,pawfgrtab,&
& rhodum,pawnhat_ptr,pawrhoij_ptr,pawrhoij0_ptr,pawtab,qphon,rprimd,ucvol,usewvl,xred,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
& comm_fft=mpi_enreg%comm_fft,paral_kgb=paral_kgb,me_g0=mpi_enreg%me_g0,&
& distribfft=mpi_enreg%distribfft,mpi_comm_wvl=mpi_enreg%comm_wvl) 

!Transfer pseudo density from coarse grid to fine grid
 if(usewvl==0) then
   call transgrid(cplex,mpi_enreg,nspden,+1,1,0,paral_kgb,pawfgr,rhopsg,rhodum,rhopsr,rhor)
 end if

!Add pseudo density and compensation charge density (on fine grid)
 rhor(:,:)=rhor(:,:)+pawnhat_ptr(:,:)

!Free temporary memory spaces
 if (.not.present(pawnhat)) then
   ABI_DEALLOCATE(pawnhat_ptr)
 end if
 if (size(pawrhoij)==natom.and.(my_natom/=natom)) then
   call pawrhoij_free(pawrhoij_ptr)
   ABI_DATATYPE_DEALLOCATE(pawrhoij_ptr)
 end if
 nullify(pawnhat_ptr)
 nullify(pawrhoij_ptr)

!Compute compensated pseudo density in reciprocal space
 if (present(rhog)) then
   call fourdp(cplex,rhog,rhor(:,1),-1,mpi_enreg,pawfgr%nfft,pawfgr%ngfft,paral_kgb,0)
 end if

 call timab(556,2,tsec)

 DBG_EXIT("COLL")

end subroutine pawmkrho
!!***
