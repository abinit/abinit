!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_longwave
!! NAME
!!  m_longwave
!!
!! FUNCTION
!!  DFPT long-wave calculation of spatial dispersion properties
!!
!! COPYRIGHT
!!  Copyright (C) 2019 ABINIT group (MR, MS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

module m_longwave
    
 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xcdata
 use m_hdr
 use m_ebands
 use m_wffile

 use m_pspini,      only : pspini
 use m_dfpt_lw,     only : dfpt_qdrpole, dfpt_flexo
 use m_common,      only : setup1
 use m_pawfgr,      only : pawfgr_type, pawfgr_init
 use m_pawrhoij,    only : pawrhoij_type
 use m_paw_dmft,    only : paw_dmft_type
 use m_pawang,      only : pawang_type
 use m_pawrad,      only : pawrad_type
 use m_pawtab,      only : pawtab_type
 use m_drivexc,     only : check_kxc
 use m_rhotoxc,     only : rhotoxc
 use m_ioarr,       only : read_rhor
 use m_symtk,       only : matr3inv,symmetrize_xred 
 use m_kg,          only : kpgio
 use m_inwffil,     only : inwffil
 use m_spacepar,    only : setsym
 use m_mkrho,       only : mkrho

 implicit none

 private
!!***

 public :: longwave
!!***

! *************************************************************************

contains 
!!***

!!****f* ABINIT/longwave
!! NAME
!!  longwave
!!
!! FUNCTION
!! Primary routine for conducting DFPT calculations of spatial dispersion properties
!!
!! INPUTS
!!  codvsn = code version
!!  dtfil <type(datafiles_type)> = variables related to files
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  etotal = new total energy (no meaning at output)
!!  iexit= exit flag
!!  mpi_enreg=informations about MPI pnarallelization
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!  npwtot(nkpt) = total number of plane waves at each k point
!!
!! SIDE EFFECTS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!
!! NOTES
!!
!! PARENTS
!!  driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine longwave(codvsn,dtfil,dtset,etotal,iexit,mpi_enreg,npwtot,occ,&
&                   pawang,pawrad,pawtab,psps,xred)
    

 implicit none

!Arguments ------------------------------------
 !scalars
 integer,intent(in) :: iexit
 real(dp),intent(inout) :: etotal
 character(len=6),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
 !arrays
 integer,intent(out) :: npwtot(dtset%nkpt)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol),xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
 !scalars
 integer,parameter :: cplex1=1,formeig=0,response=1
 integer :: ask_accurate,bantot,gscase,iatom,ierr,indx,ireadwf0,iscf_eff,itypat
 integer :: mcg,mgfftf,natom,nfftf,nfftot,nfftotf,nhatdim,nhatgrdim,mk1mem
 integer :: nkxc,nk3xc,ntypat,n3xccc
 integer :: option,optorth,psp_gencond,rdwrpaw,spaceworld,tim_mkrho,usexcnhat
 real(dp) :: ecore,ecutdg_eff,ecut_eff,enxc,etot,fermie,gsqcut_eff,gsqcutc_eff,residm
 real(dp) :: ucvol,vxcavg
 logical :: non_magnetic_xc
 character(len=500) :: msg
 type(ebands_t) :: bstruct
 type(paw_dmft_type) :: paw_dmft
 type(pawfgr_type) :: pawfgr
 type(hdr_type) :: hdr,hdr_den
 type(xcdata_type) :: xcdata
 type(wvl_data) :: wvl
 type(wffile_type) :: wffgs,wfftgs
 !arrays
 integer :: ngfft(18),ngfftf(18)
 real(dp) :: gmet(3,3),gmet_for_kg(3,3),gprimd(3,3),gprimd_for_kg(3,3)
 real(dp) :: rmet(3,3),rprimd(3,3),rprimd_for_kg(3,3)
 real(dp) :: strsxc(6)
 integer,allocatable :: atindx(:),atindx1(:),indsym(:,:,:),irrzon(:,:,:),kg(:,:)
 integer,allocatable :: nattyp(:),npwarr(:),symrec(:,:,:) 
 real(dp),allocatable :: cg(:,:)
 real(dp),allocatable :: doccde(:),eigen0(:),kxc(:,:),vxc(:,:),nhat(:,:),nhatgr(:,:,:)
 real(dp),allocatable :: phnons(:,:,:),rhog(:,:),rhor(:,:)
 real(dp),allocatable :: xccc3d(:)
 type(pawrhoij_type),allocatable :: pawrhoij(:),pawrhoij_read(:)
! *************************************************************************

 DBG_ENTER("COLL")

!Not valid for PAW
 if (psps%usepaw==1) then
   msg='This routine cannot be used for PAW (use pawnst3 instead) !'
   MSG_BUG(msg)
 end if

!Not valid for finite wave-vector perturbations
 if (sqrt(sum(dtset%qptn**2))/=0_dp) then
   msg='This routine cannot be used for q=/0.d0'
   MSG_BUG(msg)
 end if

!Only usable with spherical harmonics
 if (dtset%useylm/=1) then
   msg='This routine cannot be used for uselim/=1'
   MSG_BUG(msg)
 end if

!Not valid for spin-dependent calculations
 if (dtset%nspinor/=1.or.dtset%nsppol/=1.or.dtset%nspden/=1) then
   msg='This routine cannot be used for spin-dependent calculations'
   MSG_BUG(msg)
 end if

!Not usable with core electron density corrections
 if (psps%n1xccc/=0) then
   msg='This routine cannot be used for n1xccc/=0'
   MSG_BUG(msg)
 end if

!Define some data 
 ntypat=psps%ntypat
 natom=dtset%natom
! mk1mem=mkmems(3)

!Init spaceworld
 spaceworld=mpi_enreg%comm_cell

!Define FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfft,ngfftf)
 nfftot=product(ngfft(1:3))
 nfftotf=product(ngfftf(1:3))

!Set up for iterations
 call setup1(dtset%acell_orig(1:3,1),bantot,dtset,&
& ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,&
& natom,ngfftf,ngfft,dtset%nkpt,dtset%nsppol,&
& response,rmet,dtset%rprim_orig(1:3,1:3,1),rprimd,ucvol,psps%usepaw)

!In some cases (e.g. getcell/=0), the plane wave vectors have
! to be generated from the original simulation cell
 rprimd_for_kg=rprimd
 if (dtset%getcell/=0.and.dtset%usewvl==0) rprimd_for_kg=dtset%rprimd_orig(:,:,1)
 call matr3inv(rprimd_for_kg,gprimd_for_kg)
 gmet_for_kg=matmul(transpose(gprimd_for_kg),gprimd_for_kg)

!Set up the basis sphere of planewaves
 ABI_ALLOCATE(kg,(3,dtset%mpw*dtset%mkmem))
 ABI_ALLOCATE(npwarr,(dtset%nkpt))
 call kpgio(ecut_eff,dtset%exchn2n3d,gmet_for_kg,dtset%istwfk,kg,&
& dtset%kptns,dtset%mkmem,dtset%nband,dtset%nkpt,'PERS',mpi_enreg,dtset%mpw,npwarr,npwtot,&
& dtset%nsppol)

!Open and read pseudopotential files
ecore=zero
 call pspini(dtset,dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,pawrad,pawtab,&
& psps,rprimd,comm_mpi=mpi_enreg%comm_cell)

!Initialize band structure datatype
 bstruct = ebands_from_dtset(dtset, npwarr)

!Initialize PAW atomic occupancies to zero
 ABI_DATATYPE_ALLOCATE(pawrhoij,(0))

!Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,hdr,pawtab,gscase,psps,wvl%descr, &
& comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab)

!Update header, with evolving variables, when available
!Here, rprimd, xred and occ are available
 etot=hdr%etot ; fermie=hdr%fermie ; residm=hdr%residm
!If parallelism over atom, hdr is distributed
 call hdr_update(hdr,bantot,etot,fermie,&
& residm,rprimd,occ,pawrhoij,xred,dtset%amu_orig(:,1), &
& comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab)

!Clean band structure datatype (should use it more in the future !)
 call ebands_free(bstruct)

!Initialize wavefunction files and wavefunctions.
 ireadwf0=1

 mcg=dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
 ABI_STAT_ALLOCATE(cg,(2,mcg), ierr)
 ABI_CHECK(ierr==0, "out-of-memory in cg")

 ABI_ALLOCATE(eigen0,(dtset%mband*dtset%nkpt*dtset%nsppol))
 eigen0(:)=zero ; ask_accurate=1
 optorth=0

 hdr%rprimd=rprimd_for_kg ! We need the rprimd that was used to generate de G vectors
 call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen0,dtset%exchn2n3d,&
& formeig,hdr,ireadwf0,dtset%istwfk,kg,dtset%kptns,&
& dtset%localrdwf,dtset%mband,mcg,dtset%mkmem,mpi_enreg,dtset%mpw,&
& dtset%nband,ngfft,dtset%nkpt,npwarr,dtset%nsppol,dtset%nsym,&
& occ,optorth,dtset%symafm,dtset%symrel,dtset%tnons,&
& dtfil%unkg,wffgs,wfftgs,dtfil%unwffgs,dtfil%fnamewffk,wvl)
 hdr%rprimd=rprimd

!Close wffgs, if it was ever opened (in inwffil)
 if (ireadwf0==1) then
   call WffClose(wffgs,ierr)
 end if

!Do symmetry stuff
 ABI_ALLOCATE(irrzon,(nfftot**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(phnons,(2,nfftot**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(indsym,(4,dtset%nsym,natom))
 ABI_ALLOCATE(symrec,(3,3,dtset%nsym))
 irrzon=0;indsym=0;symrec=0;phnons=zero
!If the density is to be computed by mkrho, need irrzon and phnons
 iscf_eff=0;if(dtset%getden==0)iscf_eff=1
 call setsym(indsym,irrzon,iscf_eff,natom,&
& nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
& phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

!Symmetrize atomic coordinates over space group elements:
 call symmetrize_xred(indsym,natom,dtset%nsym,dtset%symrel,dtset%tnons,xred)

!Generate an index table of atoms, in order for them to be used
!type after type.
 ABI_ALLOCATE(atindx,(natom))
 ABI_ALLOCATE(atindx1,(natom))
 ABI_ALLOCATE(nattyp,(ntypat))
 indx=1
 do itypat=1,ntypat
   nattyp(itypat)=0
   do iatom=1,natom
     if(dtset%typat(iatom)==itypat)then
       atindx(iatom)=indx
       atindx1(indx)=iatom
       indx=indx+1
       nattyp(itypat)=nattyp(itypat)+1
     end if
   end do
 end do

!Derivative of occupations is always zero for non metallic systems
 ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
 doccde(:)=zero

!Read ground-state charge density from diskfile in case getden /= 0
!or compute it from wfs that were read previously : rhor 

 ABI_ALLOCATE(rhog,(2,nfftf))
 ABI_ALLOCATE(rhor,(nfftf,dtset%nspden))

 if (dtset%getden /= 0 .or. dtset%irdden /= 0) then
   ! Read rho1(r) from a disk file and broadcast data.
   ! This part is not compatible with MPI-FFT (note single_proc=.True. below)

   rdwrpaw=psps%usepaw
   ABI_DATATYPE_ALLOCATE(pawrhoij_read,(0))

!  MT july 2013: Should we read rhoij from the density file ?
   call read_rhor(dtfil%fildensin, cplex1, dtset%nspden, nfftf, ngfftf, rdwrpaw, mpi_enreg, rhor, &
   hdr_den, pawrhoij_read, spaceworld, check_hdr=hdr)
   etotal = hdr_den%etot; call hdr_free(hdr_den)

   ABI_DATATYPE_DEALLOCATE(pawrhoij_read)
 else
!  Obtain the charge density from read wfs
!  Be careful: in PAW, compensation density has to be added !
   tim_mkrho=4
   paw_dmft%use_sc_dmft=0 ! respfn with dmft not implemented
   paw_dmft%use_dmft=0 ! respfn with dmft not implemented

     call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&     mpi_enreg,npwarr,occ,paw_dmft,phnons,rhog,rhor,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs)
 end if ! getden

!Set up xc potential. Compute kxc here.
 option=2 ; nk3xc=1
 nkxc=2*min(dtset%nspden,2)-1;if(dtset%xclevel==2)nkxc=12*min(dtset%nspden,2)-5
 call check_kxc(dtset%ixc,dtset%optdriver)
 ABI_ALLOCATE(kxc,(nfftf,nkxc))
 ABI_ALLOCATE(vxc,(nfftf,dtset%nspden))

 call xcdata_init(xcdata,dtset=dtset)

 nhatgrdim=0;nhatdim=0
 ABI_ALLOCATE(nhat,(0,0))
 ABI_ALLOCATE(nhatgr,(0,0,0))
 n3xccc=0
 ABI_ALLOCATE(xccc3d,(n3xccc))
 non_magnetic_xc=.false.
 
 enxc=zero; usexcnhat=0
 call rhotoxc(enxc,kxc,mpi_enreg,nfftf,ngfftf,&
& nhat,nhatdim,nhatgr,nhatgrdim,nkxc,nk3xc,non_magnetic_xc,n3xccc,option,dtset%paral_kgb,rhor,&
& rprimd,strsxc,usexcnhat,vxc,vxcavg,xccc3d,xcdata)

!Deallocations
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(atindx1)
 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(irrzon)
 ABI_DEALLOCATE(nattyp)
 ABI_DEALLOCATE(kg)
 ABI_DEALLOCATE(kxc)
 ABI_DEALLOCATE(npwarr)
 ABI_DEALLOCATE(phnons)
 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(rhor)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(vxc)
 

 DBG_EXIT("COLL")

end subroutine longwave
!!***

end module m_longwave
!!***
