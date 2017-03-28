!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_dyfro
!! NAME
!! dfpt_dyfro
!!
!! FUNCTION
!! Compute the different parts of the frozen-wavefunction part of
!! the dynamical matrix, except the non-local one, computed previously.
!! Also (when installed) symmetrize the different part and their sum.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2017 ABINIT group (XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dyfr_cplex=1 if dyfrnl is real, 2 if it is complex
!!  dyfr_nondiag=1 if dyfrnl and dyfrwf are non diagonal with respect to atoms; 0 otherwise
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!     (bohr**-1)
!!  gsqcut=cutoff on G^2 based on ecut
!!  indsym(4,nsym,natom)=index showing transformation of atom labels
!!   under symmetry operations (computed in symatm)
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  mqgrid=dimensioned number of q grid points for local psp spline
!!  natom=number of atoms in unit cell
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  qgrid(mqgrid)=q point array for local psp spline fits
!!  qphon(3)=wavevector of the phonon
!!  rhog(2,nfft)=electron density in G space
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  symq(4,2,nsym)=1 if symmetry preserves present qpoint. From littlegroup_q
!!  symrec(3,3,nsym)=symmetries in reciprocal space
!!  typat(natom)=integer type for each atom in cell
!!  ucvol=unit cell volume (bohr**3).
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vlspl(mqgrid,2,ntypat)=q^2 v(q) spline for each type of atom.
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree) in real
!!   space--only used when n1xccc/=0
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and five derivatives,
!!   for each type of atom, from psp
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced coordinates for atoms in unit cell
!!
!! OUTPUT
!!  dyfrlo(3,3,natom)=frozen wavefunctions part of the dynamical matrix
!!                    (local only)
!!  dyfrwf(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)=
!!                    frozen wavefunctions part of the dynamical matrix
!!                    (local + non-local)
!!                    If NCPP, it depends on one atom
!!                    If PAW,  it depends on two atoms
!!  dyfrxc(3,3,natom)=frozen wavefunctions part of the dynamical matrix
!!                    (non-linear xc core correction)
!!
!! SIDE EFFECTS
!! Input/Output
!!  dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)=
!!                    frozen wavefunctions part of the dynamical matrix
!!                    (non-local only)
!!                    If NCPP, it depends on one atom
!!                    If PAW,  it depends on two atoms
!!
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      atm2fft,dfpt_sydy,fourdp,mkcore,mklocl_recipspace,timab,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_dyfro(atindx1,dyfrnl,dyfrlo,dyfrwf,dyfrxc,dyfr_cplex,dyfr_nondiag,&
&  gmet,gprimd,gsqcut,indsym,mgfft,mpi_enreg,mqgrid,natom,nattyp,&
&  nfft,ngfft,nspden,nsym,ntypat,n1xccc,n3xccc,paral_kgb,psps,pawtab,ph1d,qgrid,&
&  qphon,rhog,rprimd,symq,symrec,typat,ucvol,usepaw,vlspl,vxc,&
&  xcccrc,xccc1d,xccc3d,xred)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use defs_datatypes,only : pseudopotential_type
 use m_pawtab,      only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_dyfro'
 use interfaces_18_timing
 use interfaces_53_ffts
 use interfaces_56_xc
 use interfaces_64_psp
 use interfaces_67_common
 use interfaces_72_response, except_this_one => dfpt_dyfro
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dyfr_cplex,dyfr_nondiag,mgfft,mqgrid,n1xccc,n3xccc,natom,nfft,nspden
 integer,intent(in) :: nsym,ntypat,paral_kgb,usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(pseudopotential_type),intent(in) :: psps
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),indsym(4,nsym,natom),nattyp(ntypat)
 integer,intent(in) :: ngfft(18),symq(4,2,nsym),symrec(3,3,nsym),typat(natom)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),qphon(3),rhog(2,nfft),rprimd(3,3)
 real(dp),intent(in) :: vlspl(mqgrid,2,ntypat),vxc(nfft,nspden)
 real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat),xcccrc(ntypat),xred(3,natom)
 real(dp),intent(inout) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag),xccc3d(n3xccc)
 real(dp),intent(out) :: dyfrlo(3,3,natom),dyfrwf(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
 real(dp),intent(inout) :: dyfrxc(3,3,natom) !vz_i
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 logical, parameter :: do_final_sym=.true.
 integer :: iatom,jatom,n1,n2,n3,optatm,optdyfr,opteltfr,optgr,option
 integer :: optn,optn2,optstr,optv
 real(dp) :: eei
!arrays
 integer  :: qprtrb(3)
 real(dp) :: dummy6(6),dum_strn(6),dum_strv(6)
 real(dp) :: tsec(2),vprtrb(2)
 real(dp) :: dum_atmrho(0),dum_atmvloc(0),dum_gauss(0),dum_grn(0),dum_grv(0),dum_eltfrxc(0)
 real(dp),allocatable :: dyfrlo_tmp1(:,:,:),dyfrlo_tmp2(:,:,:,:,:),dyfrsym_tmp(:,:,:,:,:)
 real(dp),allocatable :: gr_dum(:,:),v_dum(:),vxctotg(:,:)

! *************************************************************************

 if(nspden==4)then
   MSG_WARNING('dfpt_dyfro : DFPT with nspden=4 works at the moment just for insulators and norm-conserving psp!')
 end if

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)

 if (usepaw==1 .or. psps%nc_xccc_gspace==1) then

!  PAW or NC with nc_xccc_gspace: compute local psp and core charge contribs together
!  in reciprocal space
!  -----------------------------------------------------------------------
   call timab(563,1,tsec)
   if (n3xccc>0) then
     ABI_ALLOCATE(v_dum,(nfft))
     ABI_ALLOCATE(vxctotg,(2,nfft))
     v_dum(:)=vxc(:,1);if (nspden>=2) v_dum(:)=0.5_dp*(v_dum(:)+vxc(:,2))
     call fourdp(1,vxctotg,v_dum,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     call zerosym(vxctotg,2,ngfft(1),ngfft(2),ngfft(3),&
&     comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     ABI_DEALLOCATE(v_dum)
   else
     ABI_ALLOCATE(vxctotg,(0,0))
   end if
   optatm=0;optdyfr=1;optgr=0;optstr=0;optv=1;optn=n3xccc/nfft;optn2=1;opteltfr=0
   call atm2fft(atindx1,dum_atmrho,dum_atmvloc,dyfrxc,dyfrlo,dum_eltfrxc,&
&   dum_gauss,gmet,gprimd,dum_grn,dum_grv,gsqcut,mgfft,mqgrid,natom,nattyp,nfft,ngfft,ntypat,&
&   optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,psps,pawtab,ph1d,qgrid,qprtrb,rhog,&
&   dum_strn,dum_strv,ucvol,usepaw,vxctotg,vxctotg,vxctotg,vprtrb,vlspl)
   ABI_DEALLOCATE(vxctotg)
   if (n3xccc==0) dyfrxc=zero
 else

!  Norm-conserving: compute local psp contribution in reciprocal space
!  and core charge contribution in real space
!  -----------------------------------------------------------------------
   option=4
   ABI_ALLOCATE(dyfrlo_tmp1,(3,3,natom))
   ABI_ALLOCATE(gr_dum,(3,natom))
   ABI_ALLOCATE(v_dum,(nfft))
   call mklocl_recipspace(dyfrlo_tmp1,eei,gmet,gprimd,&
&   gr_dum,gsqcut,dummy6,mgfft,mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,&
&   ntypat,option,paral_kgb,ph1d,qgrid,qprtrb,rhog,ucvol,vlspl,vprtrb,v_dum)
   do iatom=1,natom
!    Reestablish correct order of atoms
     dyfrlo(1:3,1:3,atindx1(iatom))=dyfrlo_tmp1(1:3,1:3,iatom)
   end do
   ABI_DEALLOCATE(dyfrlo_tmp1)
   ABI_DEALLOCATE(v_dum)
   if(n1xccc/=0)then
     call mkcore(dummy6,dyfrxc,gr_dum,mpi_enreg,natom,nfft,nspden,ntypat,&
&     n1,n1xccc,n2,n3,option,rprimd,typat,ucvol,vxc,xcccrc,xccc1d,xccc3d,xred)
   end if
   ABI_DEALLOCATE(gr_dum)
 end if

!Symmetrize dynamical matrix explicitly for given space group:

!Symmetrize local part of the dynamical matrix dyfrlo:
 ABI_ALLOCATE(dyfrsym_tmp,(1,3,3,natom,1))
 ABI_ALLOCATE(dyfrlo_tmp2,(1,3,3,natom,1))
 dyfrsym_tmp(1,:,:,:,1)=dyfrlo(:,:,:)
 call dfpt_sydy(1,dyfrsym_tmp,indsym,natom,0,nsym,qphon,dyfrlo_tmp2,symq,symrec)
 dyfrlo(:,:,:)=dyfrlo_tmp2(1,:,:,:,1)
 if (do_final_sym) then
   dyfrsym_tmp(1,:,:,:,1)=dyfrxc(:,:,:)
   call dfpt_sydy(1,dyfrsym_tmp,indsym,natom,0,nsym,qphon,dyfrlo_tmp2,symq,symrec)
   dyfrxc(:,:,:)=dyfrlo_tmp2(1,:,:,:,1)
 end if
 ABI_DEALLOCATE(dyfrsym_tmp)
 ABI_DEALLOCATE(dyfrlo_tmp2)

!Symmetrize nonlocal part of the dynamical matrix dyfrnl:
!atindx1 is used to reestablish the correct order of atoms
 if (dyfr_nondiag==0) then
   ABI_ALLOCATE(dyfrsym_tmp,(dyfr_cplex,3,3,natom,1))
   do iatom=1,natom
     dyfrsym_tmp(:,:,:,atindx1(iatom),1)=dyfrnl(:,:,:,iatom,1)
   end do
 else
   ABI_ALLOCATE(dyfrsym_tmp,(dyfr_cplex,3,3,natom,natom))
   do jatom=1,natom
     do iatom=1,natom
       dyfrsym_tmp(:,:,:,atindx1(iatom),atindx1(jatom))=dyfrnl(:,:,:,iatom,jatom)
     end do
   end do
 end if
 call dfpt_sydy(dyfr_cplex,dyfrsym_tmp,indsym,natom,dyfr_nondiag,nsym,qphon,dyfrnl,symq,symrec)
 ABI_DEALLOCATE(dyfrsym_tmp)

!Collect local, nl xc core, and non-local part
!of the frozen wf dynamical matrix.
 dyfrwf(:,:,:,:,:)=dyfrnl(:,:,:,:,:)
 if (dyfr_nondiag==0) then
   dyfrwf(1,:,:,:,1)=dyfrwf(1,:,:,:,1)+dyfrlo(:,:,:)+dyfrxc(:,:,:)
 else
   do iatom=1,natom
     dyfrwf(1,:,:,iatom,iatom)=dyfrwf(1,:,:,iatom,iatom)+dyfrlo(:,:,iatom)+dyfrxc(:,:,iatom)
   end do
 end if

end subroutine dfpt_dyfro
!!***
