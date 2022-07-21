!!****m* ABINIT/m_longwave
!! NAME
!!  m_longwave
!!
!! FUNCTION
!!  DFPT long-wave calculation of spatial dispersion properties
!!
!! COPYRIGHT
!!  Copyright (C) 2019-2022 ABINIT group (MR, MS)
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
 use m_errors
 use m_xmpi
 use defs_datatypes
 use defs_abitypes, only : MPI_type
 use defs_wvltypes
 use m_dtset
 use m_dtfil
 use m_xcdata
 use m_hdr
 use m_ebands
 use m_wffile

 use m_pspini,      only : pspini
 use m_common,      only : setup1
 use m_pawfgr,      only : pawfgr_type, pawfgr_init
 use m_pawrhoij,    only : pawrhoij_type
 use m_paw_dmft,    only : paw_dmft_type
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
 use m_fft,         only : fourdp
 use m_ddb,         only : DDB_VERSION,dfpt_lw_doutput,lwcart
 use m_ddb_hdr,     only : ddb_hdr_type, ddb_hdr_init
 use m_mkcore,      only : mkcore
 use m_dfptlw_loop, only : dfptlw_loop
 use m_dfptlw_nv,   only : dfptlw_nv
 use m_dfptlw_pert, only : preca_ffnl
 use m_initylmg,    only : initylmg
 use m_dynmat,      only : d3lwsym, sylwtens
 use m_geometry,    only : symredcart

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
!!  mpi_enreg=information about MPI pnarallelization
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!  npwtot(nkpt) = total number of plane waves at each k point
!!
!! SIDE EFFECTS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!
!! NOTES
!!
!! PARENTS
!!      m_driver
!!
!! CHILDREN
!!      check_kxc,ddb_hdr%free,ddb_hdr%open_write,ddb_hdr_init,
!!      dfpt_lw_doutput,ebands_free
!!      fourdp,hdr%free,hdr%update,hdr_init,inwffil,kpgio,matr3inv,mkcore,mkrho
!!      pawfgr_init,pspini,read_rhor,rhotoxc,setsym,setup1,symmetrize_xred
!!      wffclose,xcdata_init
!!
!! SOURCE

subroutine longwave(codvsn,dtfil,dtset,etotal,mpi_enreg,npwtot,occ,&
&                   pawrad,pawtab,psps,xred)

#ifdef FC_INTEL
!DEC$ NOOPTIMIZE
#endif


 implicit none

!Arguments ------------------------------------
 !scalars
 real(dp),intent(inout) :: etotal
 character(len=8),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
 !arrays
 integer,intent(out) :: npwtot(dtset%nkpt)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol),xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
 !scalars
 integer,parameter :: cplex1=1,formeig=0,response=1
 integer :: ask_accurate,bantot,coredens_method,dimffnl,gscase,iatom,ierr,indx,ireadwf0,iscf_eff,itypat
 integer :: ider,idir0
 integer :: i1dir,i1pert,i2dir,ii,i2pert,i3dir,i3pert
 integer :: isym,mcg,mgfftf,natom,nfftf,nfftot,nfftotf,nhatdim,nhatgrdim
 integer :: mpert,my_natom,n1,nkxc,nk3xc,ntypat,n3xccc,nylmgr
 integer :: option,optorth,psp_gencond,rdwrpaw,spaceworld,timrev,tim_mkrho
 integer :: usexcnhat,useylmgr
! integer :: idir,ipert,
 real(dp) :: ecore,ecutdg_eff,ecut_eff,enxc,etot,fermie,fermih,gsqcut_eff,gsqcutc_eff,residm ! CP added fermih
 real(dp) :: ucvol,vxcavg
 logical :: has_strain,non_magnetic_xc
 character(len=fnlen) :: dscrpt
 character(len=500) :: msg
 type(ebands_t) :: bstruct
 type(ddb_hdr_type) :: ddb_hdr
 type(paw_dmft_type) :: paw_dmft
 type(pawfgr_type) :: pawfgr
 type(hdr_type) :: hdr,hdr_den
 type(xcdata_type) :: xcdata
 type(wvl_data) :: wvl
 type(wffile_type) :: wffgs,wfftgs
 !arrays
 integer :: ngfft(18),ngfftf(18),perm(6)
 real(dp) :: dummy6(6),gmet(3,3),gmet_for_kg(3,3),gprimd(3,3),gprimd_for_kg(3,3)
 real(dp) :: rmet(3,3),rprimd(3,3),rprimd_for_kg(3,3)
 real(dp) :: strsxc(6)
 integer,allocatable :: atindx(:),atindx1(:)
 integer,allocatable :: blkflg(:,:,:,:,:,:),blkflg_car(:,:,:,:,:,:)
 integer,allocatable :: d3e_pert1(:),d3e_pert2(:),d3e_pert3(:)
 integer,allocatable :: indsym(:,:,:),irrzon(:,:,:),kg(:,:)
 integer,allocatable :: nattyp(:),npwarr(:),pertsy(:,:),symrec(:,:,:)
 integer,allocatable :: rfpert(:,:,:,:,:,:)
 real(dp),allocatable :: cg(:,:)
 real(dp),allocatable :: d3etot(:,:,:,:,:,:,:),d3etot_car(:,:,:,:,:,:,:)
 real(dp),allocatable :: d3etot_nv(:,:,:,:,:,:,:),doccde(:)
 real(dp),allocatable :: eigen0(:),ffnl(:,:,:,:,:)
 real(dp),allocatable :: grxc(:,:),kxc(:,:),vxc(:,:),nhat(:,:),nhatgr(:,:,:)
 real(dp),allocatable :: phnons(:,:,:),rhog(:,:),rhor(:,:),dummy_dyfrx2(:,:,:)
 real(dp),allocatable :: symrel_cart(:,:,:),work(:),xccc3d(:)
 real(dp),allocatable :: ylm(:,:),ylmgr(:,:,:)
 type(pawrhoij_type),allocatable :: pawrhoij(:),pawrhoij_read(:)
! *************************************************************************

 DBG_ENTER("COLL")

!Not valid for PAW
 if (psps%usepaw==1) then
   msg='This routine cannot be used for PAW!'
   ABI_BUG(msg)
 end if

!Not valid for finite wave-vector perturbations
 if (sqrt(sum(dtset%qptn**2))/=0_dp) then
   msg='This routine cannot be used for q /= 0 '
   ABI_BUG(msg)
 end if

!Only usable with spherical harmonics
 if (dtset%useylm/=1) then
   msg='This routine cannot be used for useylm/=1'
   ABI_BUG(msg)
 end if

!Not valid for spin-dependent calculations
 if (dtset%nspinor/=1.or.dtset%nsppol/=1.or.dtset%nspden/=1) then
   msg='This routine cannot be used for spin-dependent calculations'
   ABI_BUG(msg)
 end if

!Not usable with core electron density corrections
 if (psps%n1xccc/=0) then
   msg='This routine cannot be used for n1xccc/=0'
   ABI_BUG(msg)
 end if

!Define some data
 ntypat=psps%ntypat
 natom=dtset%natom
 timrev=1

!Init spaceworld
 spaceworld=mpi_enreg%comm_cell
 my_natom=mpi_enreg%my_natom

!Define FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfft,ngfftf)
 nfftot=product(ngfft(1:3))
 nfftotf=product(ngfftf(1:3))

!Set up for iterations
 call setup1(dtset%acell_orig(1:3,1),bantot,dtset,&
& ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,&
& ngfftf,ngfft,dtset%nkpt,dtset%nsppol,&
& response,rmet,dtset%rprim_orig(1:3,1:3,1),rprimd,ucvol,psps%usepaw)

!Define the set of admitted perturbations taking into account
!the possible permutations
!  -> natom+8 refers to ddq perturbation
 mpert=natom+8
 ABI_MALLOC(blkflg,(3,mpert,3,mpert,3,mpert))
 ABI_MALLOC(d3etot,(2,3,mpert,3,mpert,3,mpert))
 ABI_MALLOC(d3etot_nv,(2,3,mpert,3,mpert,3,mpert))
 ABI_MALLOC(rfpert,(3,mpert,3,mpert,3,mpert))
 ABI_MALLOC(d3e_pert1,(mpert))
 ABI_MALLOC(d3e_pert2,(mpert))
 ABI_MALLOC(d3e_pert3,(mpert))
 blkflg(:,:,:,:,:,:) = 0
 d3etot(:,:,:,:,:,:,:) = zero
 d3etot_nv(:,:,:,:,:,:,:) = zero
 rfpert(:,:,:,:,:,:) = 0
 d3e_pert1(:) = 0 ; d3e_pert2(:) = 0 ; d3e_pert3(:) = 0

 d3e_pert3(natom+8)=1
 if (dtset%lw_qdrpl==1) then
   d3e_pert1(natom+2)=1
   d3e_pert2(1:natom)=1
 end if

 if (dtset%lw_flexo==2.or.dtset%lw_flexo==1) then
   d3e_pert1(natom+2)=1
   d3e_pert2(natom+3:natom+4)=1
 end if

 if (dtset%lw_flexo==3.or.dtset%lw_flexo==1) then
   d3e_pert1(natom+2)=1 ; d3e_pert1(1:natom)=1
   d3e_pert2(1:natom)=1
 end if

 if (dtset%lw_flexo==4.or.dtset%lw_flexo==1) then
   d3e_pert1(1:natom)=1
   d3e_pert2(natom+3:natom+4)=1
 end if

 if (dtset%lw_natopt==1) then
   d3e_pert1(natom+2)=1
   !AZ*******************
   d3e_pert2(1:natom)=1  
   !AZ*******************
 end if

 perm(:)=0
 do i1pert = 1, mpert
   do i2pert = 1, mpert
     do i3pert = 1, mpert
       perm(1)=d3e_pert1(i1pert)*d3e_pert2(i2pert)*d3e_pert3(i3pert)
!       perm(2)=d3e_pert1(i1pert)*d3e_pert2(i3pert)*d3e_pert3(i2pert)
!       perm(3)=d3e_pert1(i2pert)*d3e_pert2(i1pert)*d3e_pert3(i3pert)
!       perm(4)=d3e_pert1(i2pert)*d3e_pert2(i3pert)*d3e_pert3(i1pert)
!       perm(5)=d3e_pert1(i3pert)*d3e_pert2(i2pert)*d3e_pert3(i1pert)
!       perm(6)=d3e_pert1(i3pert)*d3e_pert2(i1pert)*d3e_pert3(i2pert)
       if ( sum(perm(:)) > 0 ) rfpert(:,i1pert,:,i2pert,:,i3pert)=1
     end do
   end do
 end do 

!Do symmetry stuff
 ABI_MALLOC(irrzon,(nfftot**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_MALLOC(phnons,(2,nfftot**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_MALLOC(indsym,(4,dtset%nsym,natom))
 ABI_MALLOC(symrec,(3,3,dtset%nsym))
 irrzon=0;indsym=0;symrec=0;phnons=zero
!If the density is to be computed by mkrho, need irrzon and phnons
 iscf_eff=0;if(dtset%getden==0)iscf_eff=1
 call setsym(indsym,irrzon,iscf_eff,natom,&
& nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
& phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

!Symmetrize atomic coordinates over space group elements:
 call symmetrize_xred(natom,dtset%nsym,dtset%symrel,dtset%tnons,xred,indsym=indsym)

! Get symmetries in cartesian coordinates
 ABI_MALLOC(symrel_cart, (3, 3, dtset%nsym))
 do isym =1,dtset%nsym
   call symredcart(rprimd, gprimd, symrel_cart(:,:,isym), dtset%symrel(:,:,isym))
   ! purify operations in cartesian coordinates.
   where (abs(symrel_cart(:,:,isym)) < tol14)
     symrel_cart(:,:,isym) = zero
   end where
 end do

 call sylwtens(indsym,mpert,natom,dtset%nsym,rfpert,symrec,dtset%symrel,symrel_cart)

 write(msg,'(a,a,a,a,a)') ch10, &
& ' The list of irreducible elements of the spatial-dispersion third-order energy derivatives is: ', ch10,& 
& ' (in reduced coordinates except for strain pert.) ', ch10
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 write(msg,'(12x,a)') 'i1dir   i1pert  i2dir   i2pert  i3dir  i3pert'
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')
 n1 = 0
 do i3pert = 1, mpert
   do i3dir = 1, 3
     do i2pert = 1, mpert
       do i2dir = 1,3
         do i1pert = 1, mpert
           do i1dir = 1, 3
             if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
               n1 = n1 + 1
               write(msg,'(2x,i4,a,6(5x,i3))') n1,')', &
             & i1dir,i1pert,i2dir,i2pert,i3dir,i3pert
               call wrtout(ab_out,msg,'COLL')
               call wrtout(std_out,msg,'COLL')
             else if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==-2) then
               blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1
               if (dtset%prtvol>=10) then
                 n1 = n1 + 1
                 write(msg,'(2x,i4,a,6(5x,i3),a)') n1,')', &
  &               i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,' => must be zero, not computed'
                 call wrtout(ab_out,msg,'COLL')
                 call wrtout(std_out,msg,'COLL')
               end if
             else if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==-1) then
               if (dtset%prtvol>=10) then
                 n1 = n1 + 1
                 write(msg,'(2x,i4,a,6(5x,i3),a)') n1,')', &
  &               i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,' => symmetric of another element, not computed'
                 call wrtout(ab_out,msg,'COLL')
                 call wrtout(std_out,msg,'COLL')
               end if
             end if
           end do
         end do
       end do
     end do
   end do
 end do
 write(msg,'(a,a)') ch10,ch10
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

!In some cases (e.g. getcell/=0), the plane wave vectors have
!to be generated from the original simulation cell
 rprimd_for_kg=rprimd
 if (dtset%getcell/=0.and.dtset%usewvl==0) rprimd_for_kg=dtset%rprimd_orig(:,:,1)
 call matr3inv(rprimd_for_kg,gprimd_for_kg)
 gmet_for_kg=matmul(transpose(gprimd_for_kg),gprimd_for_kg)

!Set up the basis sphere of planewaves
 ABI_MALLOC(kg,(3,dtset%mpw*dtset%mkmem))
 ABI_MALLOC(npwarr,(dtset%nkpt))
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
 ABI_MALLOC(pawrhoij,(0))

!Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,hdr,pawtab,gscase,psps,wvl%descr, &
& comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab)

!Update header, with evolving variables, when available
!Here, rprimd, xred and occ are available
 etot=hdr%etot ; fermie=hdr%fermie ; fermih=hdr%fermih ; residm=hdr%residm ! CP added fermih
!If parallelism over atom, hdr is distributed
! CP modified
 call hdr%update(bantot,etot,fermie,fermih,&
& residm,rprimd,occ,pawrhoij,xred,dtset%amu_orig(:,1), &
& comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab)
! call hdr%update(bantot,etot,fermie,&
!& residm,rprimd,occ,pawrhoij,xred,dtset%amu_orig(:,1), &
!& comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab)
 ! End CP modified

!Clean band structure datatype (should use it more in the future !)
 call ebands_free(bstruct)

!Initialize wavefunction files and wavefunctions.
 ireadwf0=1

 mcg=dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
 ABI_STAT_MALLOC(cg,(2,mcg), ierr)
 ABI_CHECK(ierr==0, "out-of-memory in cg")

 ABI_MALLOC(eigen0,(dtset%mband*dtset%nkpt*dtset%nsppol))
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


!Generate an index table of atoms, in order for them to be used
!type after type.
 ABI_MALLOC(atindx,(natom))
 ABI_MALLOC(atindx1,(natom))
 ABI_MALLOC(nattyp,(ntypat))
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
 ABI_MALLOC(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
 doccde(:)=zero

!Read ground-state charge density from diskfile in case getden /= 0
!or compute it from wfs that were read previously : rhor

 ABI_MALLOC(rhog,(2,nfftf))
 ABI_MALLOC(rhor,(nfftf,dtset%nspden))

 if (dtset%getden /= 0 .or. dtset%irdden /= 0) then
   ! Read rho1(r) from a disk file and broadcast data.
   ! This part is not compatible with MPI-FFT (note single_proc=.True. below)

   rdwrpaw=psps%usepaw
   ABI_MALLOC(pawrhoij_read,(0))

!  MT july 2013: Should we read rhoij from the density file ?
   call read_rhor(dtfil%fildensin, cplex1, dtset%nspden, nfftf, ngfftf, rdwrpaw, mpi_enreg, rhor, &
   hdr_den, pawrhoij_read, spaceworld, check_hdr=hdr)
   etotal = hdr_den%etot; call hdr_den%free()

   ABI_FREE(pawrhoij_read)

!  Compute up+down rho(G) by fft
   ABI_MALLOC(work,(nfftf))
   work(:)=rhor(:,1)
   call fourdp(1,rhog,work,-1,mpi_enreg,nfftf,1,ngfftf,0)
   ABI_FREE(work)
 else
!  Obtain the charge density from read wfs
!  Be careful: in PAW, compensation density has to be added !
   tim_mkrho=4
   paw_dmft%use_sc_dmft=0 ! respfn with dmft not implemented
   paw_dmft%use_dmft=0 ! respfn with dmft not implemented

     call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&     mpi_enreg,npwarr,occ,paw_dmft,phnons,rhog,rhor,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs)
 end if ! getden
! ABI_FREE(cg)

!Pseudo core electron density by method 2
!TODO: The code is not still adapted to consider n3xccc in the long-wave
!driver.
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=nfftf
 ABI_MALLOC(xccc3d,(n3xccc))
 coredens_method=2
 if (coredens_method==2.and.psps%n1xccc/=0) then
   option=1
   ABI_MALLOC(dummy_dyfrx2,(3,3,natom)) ! dummy
   ABI_MALLOC(vxc,(0,0)) ! dummy
   ABI_MALLOC(grxc,(3,natom))
   call mkcore(dummy6,dummy_dyfrx2,grxc,mpi_enreg,natom,nfftf,dtset%nspden,ntypat,&
&   ngfftf(1),psps%n1xccc,ngfftf(2),ngfftf(3),option,rprimd,dtset%typat,ucvol,vxc,&
&   psps%xcccrc,psps%xccc1d,xccc3d,xred)
   ABI_FREE(dummy_dyfrx2) ! dummy
   ABI_FREE(vxc) ! dummy
   ABI_FREE(grxc) ! dummy
 end if

!Set up xc potential. Compute kxc here.
!TODO: Iclude nonlinear core corrections (see m_respfn_driver.F90)
 option=2 ; nk3xc=1
 nkxc=2*min(dtset%nspden,2)-1;if(dtset%xclevel==2)nkxc=12*min(dtset%nspden,2)-5
 call check_kxc(dtset%ixc,dtset%optdriver)
 ABI_MALLOC(kxc,(nfftf,nkxc))
 ABI_MALLOC(vxc,(nfftf,dtset%nspden))

 nhatgrdim=0;nhatdim=0
 ABI_MALLOC(nhat,(0,0))
 ABI_MALLOC(nhatgr,(0,0,0))
! n3xccc=0
! ABI_MALLOC(xccc3d,(n3xccc))
 non_magnetic_xc=.false.

 enxc=zero; usexcnhat=0

 call xcdata_init(xcdata,dtset=dtset)
 call rhotoxc(enxc,kxc,mpi_enreg,nfftf,ngfftf,&
& nhat,nhatdim,nhatgr,nhatgrdim,nkxc,nk3xc,non_magnetic_xc,n3xccc,option,rhor,&
& rprimd,strsxc,usexcnhat,vxc,vxcavg,xccc3d,xcdata)

!Set up the spherical harmonics (Ylm) and gradients at each k point 
 useylmgr=1; option=2 ; nylmgr=9
 ABI_MALLOC(ylm,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))               
 ABI_MALLOC(ylmgr,(dtset%mpw*dtset%mkmem,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr))
 call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,&
& psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,npwarr,dtset%nsppol,option,&
& rprimd,ylm,ylmgr)                                   

!Compute nonlocal form factors ffnl1, for all atoms and all k-points.
 if (dtset%ffnl_lw == 0) then 
   if (dtset%lw_qdrpl==1.or.dtset%lw_flexo==3) ider=1; idir0=4; dimffnl=4
   if (dtset%lw_flexo==1.or.dtset%lw_flexo==2.or.dtset%lw_flexo==4) then
     ider=2; idir0=4; dimffnl=10
   end if
   ABI_MALLOC(ffnl,(dtset%mkmem,dtset%mpw,dimffnl,psps%lmnmax,psps%ntypat))
   call preca_ffnl(dimffnl,ffnl,gmet,gprimd,ider,idir0,kg, &
 & dtset%kptns,dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw, &
 & dtset%nkpt,npwarr,nylmgr,psps,rmet,useylmgr,ylm,ylmgr)
   useylmgr=0
   ABI_FREE(ylmgr)
   ABI_MALLOC(ylmgr,(dtset%mpw*dtset%mkmem,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr))
 else if (dtset%ffnl_lw == 1) then 
   dimffnl=0
   ABI_MALLOC(ffnl,(dtset%mkmem,dtset%mpw,dimffnl,psps%lmnmax,psps%ntypat))
 end if

!TODO: This part of the implementation does not work properly to select specific directions
!      for each perturbation. This development is temporarily frozen.
!Initialize the list of perturbations rfpert
! mpert=natom+11
! ABI_MALLOC(rfpert,(mpert))
! rfpert(:)=0
! rfpert(natom+1)=1
! if (dtset%lw_qdrpl==1.or.dtset%lw_flexo==1.or.dtset%lw_flexo==3.or.dtset%lw_flexo==4 &
!&.or.dtset%d3e_pert1_phon==1.or.dtset%d3e_pert2_phon==1) then
!   if (dtset%d3e_pert1_phon==1) rfpert(dtset%d3e_pert1_atpol(1):dtset%d3e_pert1_atpol(2))=1
!   if (dtset%d3e_pert2_phon==1) rfpert(dtset%d3e_pert2_atpol(1):dtset%d3e_pert2_atpol(2))=1
! end if
! if (dtset%lw_qdrpl==1.or.dtset%lw_flexo==1.or.dtset%lw_flexo==2.or.dtset%lw_flexo==3.or.&
!& dtset%d3e_pert1_elfd==1) then
!   rfpert(natom+2)=1
!   rfpert(natom+10)=1
!   rfpert(natom+11)=1
! end if
! if (dtset%lw_flexo==1.or.dtset%lw_flexo==2.or.dtset%lw_flexo==4.or.dtset%d3e_pert2_strs/=0) then
!   if (dtset%d3e_pert2_strs==1.or.dtset%d3e_pert2_strs==3) rfpert(natom+3)=1
!   if (dtset%d3e_pert2_strs==2.or.dtset%d3e_pert2_strs==3) rfpert(natom+4)=1
! endif
!
!!Determine which directions treat for each type of perturbation
! ABI_MALLOC(pertsy,(3,natom+6))
! pertsy(:,:)=0
! !atomic displacement
! do ipert=1,natom
!   if (rfpert(ipert)==1.and.dtset%d3e_pert1_phon==1) then
!     do idir=1,3
!       if (dtset%d3e_pert1_dir(idir)==1) pertsy(idir,ipert)=1
!     end do
!   endif
!   if (rfpert(ipert)==1.and.dtset%d3e_pert2_phon==1) then
!     do idir=1,3
!       if (dtset%d3e_pert2_dir(idir)==1) pertsy(idir,ipert)=1
!     end do
!   end if
! end do
! !ddk
! do idir=1,3
!   if (dtset%d3e_pert3_dir(idir)==1) pertsy(idir,natom+1)=1
! end do
! !electric field
! if (rfpert(natom+2)==1) then
!   do idir=1,3
!     if (dtset%d3e_pert1_dir(idir)==1) pertsy(idir,natom+2)=1
!   end do
! end if
! !strain
! if (rfpert(natom+3)==1) pertsy(:,natom+3)=1
! if (rfpert(natom+4)==1) pertsy(:,natom+4)=1

!TODO:Add perturbation symmetries. See m_respfn_driver.F90.
!........

! All perturbations and directions are temporarily activated
 ABI_MALLOC(pertsy,(3,natom+6))
 pertsy(:,:)=1


!#############  SPATIAL-DISPERSION PROPERTIES CALCULATION  ###########################

!Anounce start of spatial-dispersion calculation
 write(msg, '(a,80a,a,a,a)' ) ch10,('=',ii=1,80),ch10,&
&   ' ==> Compute spatial-dispersion 3rd-order energy derivatives <== ',ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 if (dtset%prtvol>=10) then
   write(msg,'(5a)') ' CAUTION: Individual contributions to the 3rd-order energy derivatives ',ch10, &
                   & ' are not written in a unified form. Mixed cartesian/reduced coordinates ',ch10, &
                   & ' and/or type-I/type-II forms are used.'
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if 
 
!Calculate the nonvariational Ewald terms
 if (dtset%lw_flexo==1.or.dtset%lw_flexo==3.or.dtset%lw_flexo==4) then
   call dfptlw_nv(d3etot_nv,dtset,gmet,gprimd,mpert,my_natom,rfpert,rmet,rprimd,ucvol,xred,psps%ziontypat, & 
  & mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
 end if

!Main loop over the perturbations to calculate the stationary part
 call dfptlw_loop(atindx,blkflg,cg,d3e_pert1,d3e_pert2,d3etot,dimffnl,dtfil,dtset,&
& eigen0,ffnl,gmet,gprimd,&
& hdr,kg,kxc,dtset%mband,dtset%mgfft,mgfftf,&
& dtset%mkmem,dtset%mk1mem,mpert,mpi_enreg,dtset%mpw,natom,nattyp,ngfftf,nfftf,nhat,&
& dtset%nkpt,nkxc,dtset%nspinor,dtset%nsppol,npwarr,nylmgr,occ,&
& pawfgr,pawrad,pawrhoij,pawtab,&
& psps,rfpert,rhog,rhor,rmet,rprimd,ucvol,useylmgr,vxc,xred,ylm,ylmgr)

!Merge stationay and nonvariational contributions
 d3etot(:,:,:,:,:,:,:)=d3etot(:,:,:,:,:,:,:) + d3etot_nv(:,:,:,:,:,:,:)

!Complete missing elements using symmetry operations
 has_strain=.false.
 if (dtset%lw_flexo==1.or.dtset%lw_flexo==2.or.dtset%lw_flexo==4) has_strain=.true.
 call d3lwsym(blkflg,d3etot,has_strain,indsym,mpert,natom,dtset%nsym,symrec,dtset%symrel,symrel_cart)

!Deallocate global proc_distrib
 if(xmpi_paral==1) then
   ABI_FREE(mpi_enreg%proc_distrb)
 end if

!Open the formatted derivative database file, and write the
!preliminary information
 if (mpi_enreg%me == 0) then
   dscrpt=' Note : temporary (transfer) database '

   call ddb_hdr_init(ddb_hdr,dtset,psps,pawtab,DDB_VERSION,dscrpt,1,xred=xred,occ=occ)

   call ddb_hdr%open_write(dtfil%fnameabo_ddb, dtfil%unddb)

   call ddb_hdr%free()

!  Call main output routine
   call dfpt_lw_doutput(blkflg,d3etot,mpert,dtset%natom,dtset%ntypat,dtfil%unddb)

!  Close DDB
   close(dtfil%unddb)

   !Calculate spatial-dispersion quantities in Cartesian coordinates and write
   !them in abi_out
   ABI_MALLOC(blkflg_car,(3,mpert,3,mpert,3,mpert))
   ABI_MALLOC(d3etot_car,(2,3,mpert,3,mpert,3,mpert))
   call lwcart(blkflg,blkflg_car,d3etot,d3etot_car,gprimd,mpert,natom,rprimd)
   call dfptlw_out(blkflg_car,d3etot_car,dtset%lw_flexo,dtset%lw_qdrpl,dtset%lw_natopt,mpert,natom,ucvol)
 end if

!Deallocate arrays
 ABI_FREE(atindx)
 ABI_FREE(atindx1)
 ABI_FREE(blkflg)
 ABI_FREE(doccde)
 ABI_FREE(eigen0)
 ABI_FREE(ffnl)
 ABI_FREE(indsym)
 ABI_FREE(irrzon)
 ABI_FREE(nattyp)
 ABI_FREE(kg)
 ABI_FREE(kxc)
 ABI_FREE(npwarr)
 ABI_FREE(phnons)
 ABI_FREE(rhog)
 ABI_FREE(rhor)
 ABI_FREE(symrec)
 ABI_FREE(symrel_cart)
 ABI_FREE(vxc)
 ABI_FREE(d3etot)
 ABI_FREE(pertsy)
 ABI_FREE(rfpert)
 ABI_FREE(d3e_pert1)
 ABI_FREE(d3e_pert2)
 ABI_FREE(d3e_pert3)
 ABI_FREE(ylm)
 ABI_FREE(ylmgr)

 ! Clean the header
 call hdr%free()

 DBG_EXIT("COLL")

end subroutine longwave
!!***

!!****f* ABINIT/m_dfptlw_loop/dfptlw_out
!! NAME
!!  dfptlw_out
!!
!! FUNCTION
!!  Write the relevant spatial-dispersion quantities in Cartesian coordinates
!!
!! COPYRIGHT
!!  Copyright (C) 2022 ABINIT group (MR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  blkflg_car(3,mpert,3,mpert,3,mpert) =flags for each element of the 3DTE
!!  d3etot_car(2,3,mpert,3,mpert,3,mpert) =array with the cartesian thir-order derivatives
!!  lw_qdrpl= flag that activates quadrupoles calculation
!!  lw_flexo= flag that activates flexoelectric tensor calculation
!!  mpert =maximum number of ipert
!!  natom = number of atoms in unit cell
!!
!! OUTPUT
!!
!! SIDE EFFECTS
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


subroutine dfptlw_out(blkflg_car,d3etot_car,lw_flexo,lw_qdrpl,lw_natopt,mpert,natom,ucvol)

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lw_flexo,lw_qdrpl,lw_natopt,mpert,natom
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: blkflg_car(3,mpert,3,mpert,3,mpert) 
 real(dp),intent(out) :: d3etot_car(2,3,mpert,3,mpert,3,mpert)

!Local variables-------------------------------
!scalar 
 integer :: beta,delta,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,istr
!arrays
 integer,save :: idx(18)=(/1,1,2,2,3,3,3,2,3,1,2,1,2,3,1,3,1,2/)
 real(dp),allocatable :: qdrp(:,:,:,:,:,:,:)
 real(dp) :: piezoci(2),piezofr(2),celastci(2)

! *************************************************************************

 DBG_ENTER("COLL")

 i3pert=natom+8
 if (lw_qdrpl==1.or.lw_flexo==3.or.lw_flexo==1) then
   write(ab_out,'(a)')' First real-space moment of the polarization response '
   write(ab_out,'(a)')' to an atomic displacementatom, in cartesian coordinates,'
   write(ab_out,'(a)')' (1/ucvol factor not included),'
   write(ab_out,'(a)')' efidir   atom   atdir    qgrdir          real part        imaginary part'
   i1pert=natom+2
   do i3dir=1,3
     do i1dir=1,3
       do i2pert=1,natom
         do i2dir=1,3
           if (blkflg_car(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
             write(ab_out,'(4(i5,3x),2(1x,f20.10))') i1dir,i2pert,i2dir,i3dir, &
           & d3etot_car(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert), &
           & -d3etot_car(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
           end if
         end do
       end do
     end do
     write(ab_out,*)' '
   end do

   !Calculate cuadrupoles (symmetrize i1dir/i3dir)
   ABI_MALLOC(qdrp,(2,3,mpert,3,mpert,3,mpert))
   i1pert=natom+2
   do i2pert=1,natom
     do i2dir=1,3
       do i1dir=1,3
         do i3dir=1,i1dir-1
           if (blkflg_car(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
             !real part
             qdrp(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)=&
           & d3etot_car(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) + &
           & d3etot_car(2,i3dir,i1pert,i2dir,i2pert,i1dir,i3pert) 

             qdrp(1,i3dir,i1pert,i2dir,i2pert,i1dir,i3pert)=&
           & qdrp(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)

             !imaginary part
             qdrp(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)=&
           & -(d3etot_car(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) + &
           &   d3etot_car(1,i3dir,i1pert,i2dir,i2pert,i1dir,i3pert) ) 

             qdrp(2,i3dir,i1pert,i2dir,i2pert,i1dir,i3pert)=&
           & qdrp(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
           end if
         end do
         if (blkflg_car(i1dir,i1pert,i2dir,i2pert,i1dir,i3pert)==1) then
           !real part
           qdrp(1,i1dir,i1pert,i2dir,i2pert,i1dir,i3pert)=&
         & two*d3etot_car(2,i1dir,i1pert,i2dir,i2pert,i1dir,i3pert)

           !imaginary part
           qdrp(2,i1dir,i1pert,i2dir,i2pert,i1dir,i3pert)=&
         &-two*d3etot_car(1,i1dir,i1pert,i2dir,i2pert,i1dir,i3pert)
         end if
       end do
     end do
   end do

   write(ab_out,'(a)')' Quadrupole tensor, in cartesian coordinates,'
   write(ab_out,'(a)')' efidir   atom   atdir    qgrdir          real part        imaginary part'
   do i3dir=1,3
     do i1dir=1,3
       do i2pert=1,natom
         do i2dir=1,3
           if (blkflg_car(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
             write(ab_out,'(4(i5,3x),2(1x,f20.10))') i1dir,i2pert,i2dir,i3dir, &
           & qdrp(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert), &
           & qdrp(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
           end if
         end do
       end do
     end do
     write(ab_out,*)' '
   end do
   ABI_FREE(qdrp)

   write(ab_out,'(a)')' Electronic (clamped-ion) contribution to the piezoelectric tensor,'
   write(ab_out,'(a)')' in cartesian coordinates, (from sum rule of dynamic quadrupoles or P^1 tensor)'
   write(ab_out,'(a)')' efidir   atdir    qgrdir        real part           imaginary part'
   do i3dir=1,3
     do i1dir=1,3
       do i2dir=1,3
         piezoci=zero
         do i2pert=1,natom
           if (blkflg_car(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
             piezoci(1)=piezoci(1)+d3etot_car(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
             piezoci(2)=piezoci(2)-d3etot_car(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
           end if
         end do
         piezoci(1)=-piezoci(1)/ucvol
         piezoci(2)=-piezoci(2)/ucvol
         write(ab_out,'(3(i5,3x),2(1x,f20.10))') i1dir,i2dir,i3dir,piezoci(1),piezoci(2)
       end do
     end do
     write(ab_out,*)' '
   end do
 end if

 if (lw_flexo==2.or.lw_flexo==1) then
   write(ab_out,'(a)')' Clamped-ion flexoelectric tensor (type-II), in cartesian coordinates,'
   write(ab_out,'(a)')' efidir  qgrdir  strdir1  strdir2         real part          imaginary part'
   i1pert=natom+2
   do i3dir=1,3
     do i2pert=natom+3,natom+4
       do i2dir=1,3
         istr=(i2pert-natom-3)*3+i2dir
         beta=idx(2*istr-1); delta=idx(2*istr)
         do i1dir=1,3
           if (blkflg_car(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
             write(ab_out,'(4(i5,3x),2(1x,f20.10))') i1dir,i3dir,beta,delta, &
           & d3etot_car(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)/ucvol, &
           & d3etot_car(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)/ucvol
           end if
         end do
       end do
       write(ab_out,*)' '
     end do
   end do
 end if

 if (lw_flexo==3.or.lw_flexo==1) then
   write(ab_out,'(a)')' 1st real-space moment of IFCs, in cartesian coordinates,'
   write(ab_out,'(a)')' iatdir   iatom    jatdir   jatom    qgrdir           real part          imaginary part'
   do i3dir=1,3
     do i1pert=1,natom
       do i1dir=1,3
         do i2pert=1,natom
           do i2dir=1,3
             if (blkflg_car(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
               write(ab_out,'(5(i5,4x),2(1x,f20.10))') i1dir,i1pert,i2dir,i2pert,i3dir, &
             & -d3etot_car(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert),&
             &  d3etot_car(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
             end if
           end do
         end do
       end do
     end do
     write(ab_out,*)' '
   end do

   write(ab_out,'(a)')' Piezoelectric force-response tensor, in cartesian coordinates '
   write(ab_out,'(a)')' (from sum rule of 1st moment of IFCs),'
   write(ab_out,'(a)')' (for non-vanishing forces in the cell it lacks an improper contribution),'
   write(ab_out,'(a)')' iatom   iatddir  jatddir   qgrdir           real part          imaginary part'
   do i3dir=1,3
     do i1pert=1,natom
       do i1dir=1,3
         do i2dir=1,3
           piezofr=zero
           do i2pert=1,natom
             if (blkflg_car(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
               piezofr(1)=piezofr(1)-d3etot_car(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
               piezofr(2)=piezofr(2)+d3etot_car(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
             end if
           end do
           write(ab_out,'(4(i5,4x),2(1x,f20.10))') i1pert,i1dir,i2dir,i3dir, &
         & piezofr(1), piezofr(2)
         end do
       end do
     end do
     write(ab_out,*)' '
   end do
 end if

 if (lw_flexo==4.or.lw_flexo==1) then
   write(ab_out,'(a)')' Clamped-ion flexoelectric force-response tensor (type-II),  in cartesian coordinates,'
   write(ab_out,'(a)')'  atom   atdir   qgrdir  strdir1 strdir2          real part          imaginary part'
   do i3dir=1,3
     do i1pert=1,natom
       do i1dir=1,3
         do i2pert=natom+3, natom+4
           do i2dir=1,3
             istr=(i2pert-natom-3)*3+i2dir
             beta=idx(2*istr-1); delta=idx(2*istr)
             if (blkflg_car(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
               write(ab_out,'(5(i5,3x),2(1x,f20.10))') i1pert,i1dir,i3dir,beta,delta, &
             & d3etot_car(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert),&
             & d3etot_car(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
             end if
           end do
         end do
         write(ab_out,*)' '
       end do
     end do
   end do


   write(ab_out,'(a)')' Clamped-ion elastic tensor, in cartesian coordinates '
   write(ab_out,'(a)')' (from sum rule of clamped-ion flexoelectric force-response tensor),'
   write(ab_out,'(a)')' (for stressed cells it lacks an improper contribution),'
   write(ab_out,'(a)')' atdir   qgrdir  strdir1  strdir2         real part          imaginary part'
   do i1dir=1,3
     do i3dir=1,i1dir
       do i2pert=natom+3, natom+4
         do i2dir=1,3
           istr=(i2pert-natom-3)*3+i2dir
           beta=idx(2*istr-1); delta=idx(2*istr)
           celastci=zero
           do i1pert=1,natom
             if (blkflg_car(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
               celastci(1)=celastci(1)+d3etot_car(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
               celastci(2)=celastci(2)+d3etot_car(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
             end if
           end do
           write(ab_out,'(4(i5,3x),2(1x,f20.10))') i1dir,i3dir,beta,delta, &
         & celastci(1)/ucvol,celastci(2)/ucvol
         end do
       end do
       write(ab_out,*)' '
     end do
   end do
 end if

 if (lw_natopt==1) then
   write(ab_out,'(a)')' Natural optical activity tensor, in cartesian coordinates,'
   write(ab_out,'(a)')' (1/ucvol factor not included),'
   write(ab_out,'(a)')' efidir1   efidir2   qgrdir          real part          imaginary part'
   i1pert=natom+2
   i2pert=natom+2
   do i3dir=1,3
     do i1dir=1,3
       do i2dir=1,3
         if (blkflg_car(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
           write(ab_out,'(3(i5,3x),2(1x,f20.10))') i1dir,i2dir,i3dir, &
         & four*pi*d3etot_car(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert), &
         & -four*pi*d3etot_car(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
         end if
       end do
     end do
     write(ab_out,*)' '
   end do
 end if

 DBG_EXIT("COLL")

end subroutine dfptlw_out
!!***

end module m_longwave
!!***
