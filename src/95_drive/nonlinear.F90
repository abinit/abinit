!{\src2tex{textfont=tt}}
!!****f* ABINIT/nonlinear
!! NAME
!! nonlinear
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations of
!! non linear response functions.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2017 ABINIT group (MVeithen, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  codvsn = code version
!!  dtfil <type(datafiles_type)> = variables related to files
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  etotal = new total energy (no meaning at output)
!!  iexit= exit flag
!!  mband = maximum number of bands
!!  mgfft = maximum single fft dimension
!!  mkmem = maximum number of k points which can fit in core memory
!!  mpi_enreg=informations about MPI pnarallelization
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  natom = number of atoms in unit cell
!!  nfft  = (effective) number of FFT grid points (for this processor)
!!  nkpt  = number of k points
!!  nspden = number of spin-density components
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  nsym   = number of symmetry elements in space group
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!
!!  npwtot(nkpt) = total number of plane waves at each k point
!!
!! SIDE EFFECTS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      d3sym,ddb_hdr_free,ddb_hdr_init,ddb_hdr_open_write,dfptnl_doutput
!!      dfptnl_loop,ebands_free,fourdp,getcut,getkgrid,getshell,hdr_free
!!      hdr_init,hdr_update,initmv,inwffil,kpgio,mkcore,nlopt,pspini,read_rhor
!!      rhotoxc,setsym,setup1,status,symmetrize_xred,sytens,timab,wffclose
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine nonlinear(codvsn,dtfil,dtset,etotal,iexit,&
&  mband,mgfft,mkmem,mpi_enreg,mpw,natom,nfft,nkpt,npwtot,nspden,&
&  nspinor,nsppol,nsym,occ,pawrad,pawtab,psps,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_wffile
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_hdr
 use m_ebands
 use m_xcdata

 use m_dynmat,   only : d3sym, sytens
 use m_ddb,      only : nlopt, DDB_VERSION
 use m_ddb_hdr,  only : ddb_hdr_type, ddb_hdr_init, ddb_hdr_free, ddb_hdr_open_write
 use m_ioarr,    only : read_rhor
 use m_pawrad,   only : pawrad_type
 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nonlinear'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_56_xc
 use interfaces_64_psp
 use interfaces_67_common
 use interfaces_72_response
 use interfaces_79_seqpar_mpi
 use interfaces_95_drive, except_this_one => nonlinear
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iexit,mband,mgfft,mkmem,mpw,nfft
 integer,intent(in) :: natom,nkpt,nspden,nspinor,nsppol,nsym
 logical :: non_magnetic_xc
 real(dp),intent(inout) :: etotal
 character(len=6),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,intent(out) :: npwtot(nkpt)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol),xred(3,natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat,psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat,psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=50,response=1,cplex1=1
 integer :: ask_accurate,bantot,dum_nshiftk,flag
 integer :: formeig,gencond,gscase,i1dir,i1pert,i2dir,i2pert,i3dir
 integer :: i3pert,ierr,ireadwf,mcg,mkmem_max,mpert,n1,n2,n3,n3xccc
 integer :: nkpt3,nkxc,nk3xc,nneigh,option,optorth,rdwrpaw,comm_cell
 real(dp) :: boxcut,ecore,ecut_eff,enxc,fermie,gsqcut,gsqcut_eff,gsqcut_eff2,gsqcutdg_eff
 real(dp) :: rdum,residm,ucvol,vxcavg
 character(len=500) :: message
 character(len=fnlen) :: dscrpt
 type(ebands_t) :: bstruct
 type(hdr_type) :: hdr,hdr_den
 type(ddb_hdr_type) :: ddb_hdr
 type(wffile_type) :: wffgs,wfftgs
 type(wvl_data) :: wvl
 type(xcdata_type) :: xcdata
!arrays
 integer :: dum_kptrlatt(3,3),dum_vacuum(3),perm(6)
 integer,allocatable :: blkflg(:,:,:,:,:,:),carflg(:,:,:,:,:,:),cgindex(:,:)
 integer,allocatable :: d3e_pert1(:),d3e_pert2(:),d3e_pert3(:)
 integer,allocatable :: indsym(:,:,:),irrzon(:,:,:),kg(:,:),kneigh(:,:),kg_neigh(:,:,:)
 integer,allocatable :: kptindex(:,:),npwarr(:),pwind(:,:,:),rfpert(:,:,:,:,:,:),symrec(:,:,:)
 real(dp) :: dum_shiftk(3,210),dummy2(6),gmet(3,3),gprimd(3,3),k0(3)
 real(dp) :: rmet(3,3),rprimd(3,3),strsxc(6),tsec(2)
 real(dp),allocatable :: amass(:),cg(:,:),d3cart(:,:,:,:,:,:,:)
 real(dp),allocatable :: d3lo(:,:,:,:,:,:,:),dum_kptns(:,:)
 real(dp),allocatable :: dum_wtk(:),dyfrx2(:,:,:),eigen(:),grxc(:,:),k3xc(:,:)
 real(dp),allocatable :: kpt3(:,:),kxc(:,:),mvwtk(:,:),phnons(:,:,:),rhog(:,:)
 real(dp),allocatable :: rhor(:,:),vhartr(:),vxc(:,:),work(:),xccc3d(:)
 type(pawrhoij_type),allocatable :: pawrhoij(:)

! ***********************************************************************

 call timab(501,1,tsec)
 call status(0,dtfil%filstat,iexit,level,'enter         ')

 comm_cell = mpi_enreg%comm_cell

! Initialise non_magnetic_xc for rhohxc
 non_magnetic_xc=(dtset%usepawu==4).or.(dtset%usepawu==14)

!Check if the perturbations asked in the input file can be computed

 if (((dtset%d3e_pert1_phon == 1).and.(dtset%d3e_pert2_phon == 1)).or. &
& ((dtset%d3e_pert1_phon == 1).and.(dtset%d3e_pert3_phon == 1)).or. &
& ((dtset%d3e_pert2_phon == 1).and.(dtset%d3e_pert3_phon == 1))) then
   write(message,'(7a)')&
&   'You have asked for a third-order derivative with respect to',ch10,&
&   '2 or more atomic displacements.',ch10,&
&   'This is not allowed yet.',ch10,&
&   'Action : change d3e_pert1_phon, d3e_pert2_phon or d3e_pert3_phon in your input file.'
   MSG_ERROR(message)
 end if

!Define the set of admitted perturbations taking into account
!the possible permutations
 mpert=natom+6
 ABI_ALLOCATE(blkflg,(3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(carflg,(3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(rfpert,(3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3e_pert1,(mpert))
 ABI_ALLOCATE(d3e_pert2,(mpert))
 ABI_ALLOCATE(d3e_pert3,(mpert))
 ABI_ALLOCATE(d3lo,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3cart,(2,3,mpert,3,mpert,3,mpert))
 blkflg(:,:,:,:,:,:) = 0
 d3lo(:,:,:,:,:,:,:) = 0_dp
 rfpert(:,:,:,:,:,:) = 0
 d3e_pert1(:) = 0 ; d3e_pert2(:) = 0 ; d3e_pert3(:) = 0

 if (dtset%d3e_pert1_phon==1) d3e_pert1(dtset%d3e_pert1_atpol(1):dtset%d3e_pert1_atpol(2))=1
 if (dtset%d3e_pert2_phon==1) d3e_pert2(dtset%d3e_pert2_atpol(1):dtset%d3e_pert2_atpol(2))=1
 if (dtset%d3e_pert3_phon==1) d3e_pert3(dtset%d3e_pert3_atpol(1):dtset%d3e_pert3_atpol(2))=1
 if (dtset%d3e_pert1_elfd/=0) d3e_pert1(natom+2)=1
 if (dtset%d3e_pert2_elfd/=0) d3e_pert2(natom+2)=1
 if (dtset%d3e_pert3_elfd/=0) d3e_pert3(natom+2)=1

 do i1pert = 1, mpert
   do i1dir = 1, 3
     do i2pert = 1, mpert
       do i2dir = 1, 3
         do i3pert = 1, mpert
           do i3dir = 1, 3
             perm(1) = &
&             d3e_pert1(i1pert)*dtset%d3e_pert1_dir(i1dir) &
&             *d3e_pert2(i2pert)*dtset%d3e_pert2_dir(i2dir) &
&             *d3e_pert3(i3pert)*dtset%d3e_pert3_dir(i3dir)
             perm(2) = &
&             d3e_pert1(i1pert)*dtset%d3e_pert1_dir(i1dir) &
&             *d3e_pert2(i3pert)*dtset%d3e_pert2_dir(i3dir) &
&             *d3e_pert3(i2pert)*dtset%d3e_pert3_dir(i2dir)
             perm(3) = &
&             d3e_pert1(i2pert)*dtset%d3e_pert1_dir(i2dir) &
&             *d3e_pert2(i1pert)*dtset%d3e_pert2_dir(i1dir) &
&             *d3e_pert3(i3pert)*dtset%d3e_pert3_dir(i3dir)
             perm(4) = &
&             d3e_pert1(i2pert)*dtset%d3e_pert1_dir(i2dir) &
&             *d3e_pert2(i3pert)*dtset%d3e_pert2_dir(i3dir) &
&             *d3e_pert3(i1pert)*dtset%d3e_pert3_dir(i1dir)
             perm(5) = &
&             d3e_pert1(i3pert)*dtset%d3e_pert1_dir(i3dir) &
&             *d3e_pert2(i2pert)*dtset%d3e_pert2_dir(i2dir) &
&             *d3e_pert3(i1pert)*dtset%d3e_pert3_dir(i1dir)
             perm(6) = &
&             d3e_pert1(i3pert)*dtset%d3e_pert1_dir(i3dir) &
&             *d3e_pert2(i1pert)*dtset%d3e_pert2_dir(i1dir) &
&             *d3e_pert3(i2pert)*dtset%d3e_pert3_dir(i2dir)
             if (sum(perm(:)) > 0) rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1
           end do
         end do
       end do
     end do
   end do
 end do

!Determine the symmetrical perturbations
 ABI_ALLOCATE(irrzon,(dtset%nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4)))
 ABI_ALLOCATE(phnons,(2,dtset%nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4)))
 ABI_ALLOCATE(indsym,(4,nsym,natom))
 ABI_ALLOCATE(symrec,(3,3,nsym))
 call status(0,dtfil%filstat,iexit,level,'call setsym   ')
 call setsym(indsym,irrzon,dtset%iscf,natom,&
& nfft,dtset%ngfft,nspden,nsppol,nsym,&
& phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

 call status(0,dtfil%filstat,iexit,level,'call sym_xred ')
 call symmetrize_xred(indsym,natom,nsym,dtset%symrel,dtset%tnons,xred)
 call status(0,dtfil%filstat,iexit,level,'call sytens   ')
 call sytens(indsym,mpert,natom,nsym,rfpert,symrec,dtset%symrel)

 write(message, '(a,a,a,a,a)' ) ch10, &
& ' The list of irreducible elements of the Raman and non-linear',&
& ch10,' optical susceptibility tensors is:',ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message,'(12x,a)')&
& 'i1pert  i1dir   i2pert  i2dir   i3pert  i3dir'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 n1 = 0
 do i1pert = 1, natom + 2
   do i1dir = 1, 3
     do i2pert = 1, natom + 2
       do i2dir = 1, 3
         do i3pert = 1, natom + 2
           do i3dir = 1,3
             if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
               n1 = n1 + 1
               write(message,'(2x,i4,a,6(5x,i3))') n1,')', &
&               i1pert,i1dir,i2pert,i2dir,i3pert,i3dir
               call wrtout(ab_out,message,'COLL')
               call wrtout(std_out,message,'COLL')
             else if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==-2) then
               blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1
             end if
           end do
         end do
       end do
     end do
   end do
 end do
 write(message,'(a,a)') ch10,ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 !if (dtset%paral_rf == -1) then
 write(std_out,'(a)')"--- !IrredPerts"
 write(std_out,'(a)')'# List of irreducible perturbations for nonlinear'
 write(std_out,'(a)')'irred_perts:'

 n1 = 0
 do i1pert = 1, natom + 2
   do i1dir = 1, 3
     do i2pert = 1, natom + 2
       do i2dir = 1, 3
         do i3pert = 1, natom + 2
           do i3dir = 1,3
             if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
               n1 = n1 + 1
               write(std_out,'(a,i0)')"   - i1pert: ",i1pert
               write(std_out,'(a,i0)')"     i1dir: ",i1dir
               write(std_out,'(a,i0)')"     i2pert: ",i2pert
               write(std_out,'(a,i0)')"     i2dir: ",i2dir
               write(std_out,'(a,i0)')"     i3pert: ",i3pert
               write(std_out,'(a,i0)')"     i3dir: ",i3dir
             end if
           end do
         end do
       end do
     end do
   end do
 end do
 write(std_out,'(a)')"..."
   !MSG_ERROR_NODUMP("aborting now")
 !end if

!Set up for iterations
 ecut_eff= (dtset%ecut) * (dtset%dilatmx) **2
 ABI_ALLOCATE(amass,(natom))
 call status(0,dtfil%filstat,iexit,level,'call setup1   ')
 call setup1(dtset%acell_orig(1:3,1),amass,dtset%amu_orig(:,1),bantot,dtset,&
& ecut_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcut_eff2,&
& natom,dtset%ngfft,dtset%ngfft,nkpt,nsppol,&
& response,rmet,dtset%rprim_orig(1:3,1:3,1),rprimd,ucvol,psps%usepaw)

!Set up the basis sphere of planewaves
 ABI_ALLOCATE(kg,(3,mpw*dtset%mk1mem))
 ABI_ALLOCATE(npwarr,(nkpt))
 call status(0,dtfil%filstat,iexit,level,'call kpgio    ')
 call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg,&
& dtset%kptns,mkmem,dtset%nband,nkpt,'PERS',mpi_enreg,&
& mpw,npwarr,npwtot,nsppol)

!Recompute first large sphere cut-off gsqcut,
!without taking into account dilatmx
 k0(:)=0.0_dp
 call status(0,dtfil%filstat,iexit,level,'call getcut   ')
 call getcut(boxcut,dtset%ecut,gmet,gsqcut,dtset%iboxcut,std_out,k0,dtset%ngfft)

!Open and read pseudopotential files
 ecore = 0_dp
 call status(0,dtfil%filstat,iexit,level,'call pspini   ')
 call pspini(dtset,dtfil,ecore,gencond,gsqcut_eff,gsqcutdg_eff,level,&
& pawrad,pawtab,psps,rprimd,comm_mpi=mpi_enreg%comm_cell)

!Initialize band structure datatype
 bstruct = ebands_from_dtset(dtset, npwarr)

!Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,hdr,pawtab,gscase,psps,wvl%descr)

!Update header, with evolving variables, when available
!Here, rprimd, xred and occ are available
 residm=hdr%residm ; fermie=hdr%fermie
 call hdr_update(hdr,bantot,etotal,fermie,residm,rprimd,occ,pawrhoij,xred,dtset%amu_orig(:,1))

!Clean band structure datatype (should use it more in the future !)
 call ebands_free(bstruct)

!Read ground-state wavefunctions
 mcg=dtset%mpw*dtset%nspinor*mband*dtset%mkmem*dtset%nsppol
 ABI_ALLOCATE(cg,(2,mcg))
 ABI_ALLOCATE(eigen,(mband*dtset%nkpt*dtset%nsppol))
 optorth=1;if (psps%usepaw==1) optorth=0
 ireadwf=1 ; formeig=0
 eigen(:)=0_dp ; ask_accurate=1
 call status(0,dtfil%filstat,iexit,level,'call inwffil  ')
 call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen,dtset%exchn2n3d,&
& formeig,hdr,ireadwf,dtset%istwfk,kg,dtset%kptns,&
& dtset%localrdwf,mband,mcg,dtset%mkmem,mpi_enreg,mpw,&
& dtset%nband,dtset%ngfft,dtset%nkpt,&
& npwarr,dtset%nsppol,dtset%nsym,&
& occ,optorth,dtset%symafm,dtset%symrel,dtset%tnons,&
& dtfil%unkg,wffgs,wfftgs,dtfil%unwffgs,&
& dtfil%fnamewffk,wvl)

 if (ireadwf==1) then
   call WffClose(wffgs,ierr)
 end if

 ABI_DEALLOCATE(eigen)

 ABI_ALLOCATE(rhog,(2,nfft))
 ABI_ALLOCATE(rhor,(nfft,nspden))
!Get the ground state charge density

 if (dtset%getden /= 0 .or. dtset%irdden /= 0) then
   rdwrpaw = 0
   call status(0,dtfil%filstat,iexit,level,'call ioarr    ')

   call read_rhor(dtfil%fildensin, cplex1, dtset%nspden, nfft, dtset%ngfft, rdwrpaw, &
   mpi_enreg, rhor, hdr_den, pawrhoij, comm_cell, check_hdr=hdr)
   call hdr_free(hdr_den)

!  Compute up+down rho(G) by fft
   ABI_ALLOCATE(work,(nfft))
   work(:)=rhor(:,1)
   call status(0,dtfil%filstat,iexit,level,'call fourdp   ')
   call fourdp(1,rhog,work,-1,mpi_enreg,nfft,dtset%ngfft,dtset%paral_kgb,0)
   ABI_DEALLOCATE(work)
 end if

!Compute core electron density xccc3d
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 ABI_ALLOCATE(grxc,(3,natom))
 ABI_ALLOCATE(vxc,(nfft,nspden))
 ABI_ALLOCATE(vhartr,(nfft))
 n3xccc=0
 if (psps%n1xccc/=0) n3xccc=nfft
 ABI_ALLOCATE(xccc3d,(n3xccc))
 if (psps%n1xccc/=0) then
   option=1
   ABI_ALLOCATE(dyfrx2,(3,3,natom))
   call status(0,dtfil%filstat,iexit,level,'call mkcore   ')
   call mkcore(dummy2,dyfrx2,grxc,mpi_enreg,natom,nfft,nspden,psps%ntypat,&
&   n1,psps%n1xccc,n2,n3,option,rprimd,dtset%typat,ucvol,vxc,&
&   psps%xcccrc,psps%xccc1d,xccc3d,xred)
   ABI_DEALLOCATE(dyfrx2)
 end if

 call hartre(1,gsqcut,psps%usepaw,mpi_enreg,nfft,dtset%ngfft,dtset%paral_kgb,rhog,rprimd,vhartr)

!Compute kxc (second- and third-order exchange-correlation kernel)
 option=3
 nkxc=2*nspden-1 ! LDA
 if(dtset%xclevel==2.and.nspden==1) nkxc=7  ! non-polarized GGA
 if(dtset%xclevel==2.and.nspden==2) nkxc=19 ! polarized GGA
 nk3xc=3*nspden-2
 ABI_ALLOCATE(kxc,(nfft,nkxc))
 ABI_ALLOCATE(k3xc,(nfft,nk3xc))

 call status(0,dtfil%filstat,iexit,level,'call rhotoxc   ')
 ABI_ALLOCATE(work,(0))
 call xcdata_init(dtset%auxc_ixc,dtset%intxc,dtset%ixc,&
&    dtset%nelect,dtset%tphysel,dtset%usekden,dtset%vdw_xc,dtset%xc_tb09_c,dtset%xc_denpos,xcdata)
 call rhotoxc(enxc,kxc,mpi_enreg,nfft,dtset%ngfft,&
& work,0,work,0,nkxc,nk3xc,non_magnetic_xc,nspden,n3xccc,option,dtset%paral_kgb,rhor,rprimd,strsxc,1,&
& vxc,vxcavg,xccc3d,xcdata,k3xc=k3xc,vhartr=vhartr)
 ABI_DEALLOCATE(work)

 ABI_DEALLOCATE(vhartr)
 ABI_DEALLOCATE(vxc)
 ABI_DEALLOCATE(xccc3d)

!Initialize finite difference calculation of the ddk

 call status(0,dtfil%filstat,iexit,level,'call getshell ')
 nkpt3 = 0

!Prepare first call to getkgrid (obtain number of k points in FBZ)
 dum_kptrlatt(:,:) = dtset%kptrlatt(:,:)
 dum_nshiftk = dtset%nshiftk
 ABI_CHECK(dum_nshiftk<=210,"dum_nshiftk must be <= 210!")
 dum_shiftk(:,:) = zero
 dum_shiftk(:,1:dtset%nshiftk) = dtset%shiftk(:,1:dtset%nshiftk)
 dum_vacuum(:) = 0

 ABI_ALLOCATE(dum_kptns,(3,0))
 ABI_ALLOCATE(dum_wtk,(0))
 call getkgrid(0,0,dtset%iscf,dum_kptns,3,dum_kptrlatt,&
& rdum,dtset%nsym,0,nkpt3,dum_nshiftk,dtset%nsym,&
& rprimd,dum_shiftk,dtset%symafm,dtset%symrel,&
& dum_vacuum,dum_wtk)
 ABI_DEALLOCATE(dum_kptns)
 ABI_DEALLOCATE(dum_wtk)

 write(std_out,*) 'nonlinear : nkpt, nkpt3 = ',nkpt,nkpt3
!call flush(6)
!jmb : malloc() problem with gcc461_openmpi under max2 : change order of allocations works ?!?
!allocate(kneigh(30,nkpt),kg_neigh(30,nkpt,3),mvwtk(30,nkpt))
 ABI_ALLOCATE(kg_neigh,(30,nkpt,3))
 ABI_ALLOCATE(mvwtk,(30,nkpt))
 ABI_ALLOCATE(kneigh,(30,nkpt))

 ABI_ALLOCATE(kptindex,(2,nkpt3))
 ABI_ALLOCATE(kpt3,(3,nkpt3))

 call getshell(gmet,kneigh,kg_neigh,kptindex,dtset%kptopt,&
& dtset%kptrlatt,dtset%kptns,kpt3,mkmem,mkmem_max,mpi_enreg,mvwtk,&
& nkpt,nkpt3,nneigh,dtset%nshiftk,rmet,rprimd,dtset%shiftk,dtset%wtk)

 ABI_ALLOCATE(pwind,(mpw,nneigh,mkmem))
 ABI_ALLOCATE(cgindex,(nkpt,nsppol))
 ABI_ALLOCATE(mpi_enreg%kpt_loc2ibz_sp,(0:mpi_enreg%nproc-1,1:mkmem_max, 1:2))
 ABI_ALLOCATE(mpi_enreg%mkmem,(0:mpi_enreg%nproc-1))

 call status(0,dtfil%filstat,iexit,level,'call initmv   ')
 call initmv(cgindex,dtset,gmet,kg,kneigh,kg_neigh,kptindex,&
& kpt3,mband,mkmem,mpi_enreg,mpw,dtset%nband,nkpt,&
& nkpt3,nneigh,npwarr,nsppol,occ,pwind)

 call status(0,dtfil%filstat,iexit,level,'call dfptnl_loop ')
 call dfptnl_loop(blkflg,cg,cgindex,dtfil,dtset,d3lo,gmet,gprimd,gsqcut,&
& hdr,kg,kneigh,kg_neigh,kptindex,kpt3,kxc,k3xc,mband,mgfft,&
& mkmem,mkmem_max,dtset%mk1mem,mpert,mpi_enreg,mpw,mvwtk,natom,nfft,&
& nkpt,nkpt3,nkxc,nk3xc,nneigh,nspinor,nsppol,npwarr,occ,psps,pwind,&
& rfpert,rprimd,ucvol,xred)

 write(message,'(a,a,a)')ch10,&
& ' --- Third order energy calculation completed --- ',ch10
 call wrtout(ab_out,message,'COLL')

!Complete missing elements using symmetry operations
 call status(0,dtfil%filstat,iexit,level,'call d3sym    ')
 call d3sym(blkflg,d3lo,indsym,mpert,natom,nsym,&
& symrec,dtset%symrel)

!Open the formatted derivative database file, and write the
!preliminary information
 if (mpi_enreg%me == 0) then
   call status(0,dtfil%filstat,iexit,level,'call ioddb8_ou')

   dscrpt=' Note : temporary (transfer) database '

   call ddb_hdr_init(ddb_hdr,dtset,psps,pawtab,DDB_VERSION,dscrpt,&
&   1,xred=xred,occ=occ)

   call ddb_hdr_open_write(ddb_hdr, dtfil%fnameabo_ddb, dtfil%unddb)

   call ddb_hdr_free(ddb_hdr)

!  Call main output routine
   call dfptnl_doutput(blkflg,d3lo,dtset%mband,mpert,dtset%nkpt,dtset%natom,dtset%ntypat,dtfil%unddb)

!  Close DDB
   close(dtfil%unddb)

!  Compute tensors related to third-order derivatives
   call nlopt(blkflg,carflg,d3lo,d3cart,gprimd,mpert,natom,rprimd,ucvol)

   if ((d3e_pert1(natom+2)==1).and.(d3e_pert2(natom+2)==1).and. &
&   (d3e_pert3(natom+2)==1)) then

     flag = 1
     i1pert = natom+2

     d3cart(:,:,i1pert,:,i1pert,:,i1pert) = &
&     d3cart(:,:,i1pert,:,i1pert,:,i1pert)*16*(pi**2)*(Bohr_Ang**2)*1.0d-8*eps0/e_Cb

     write(ab_out,*)ch10
     write(ab_out,*)' Non-linear optical susceptibility tensor d (pm/V)'
     write(ab_out,*)' in cartesian coordinates'
     write(ab_out,*)'  i1dir  i2dir  i3dir             d'

     do i1dir = 1, 3
       do i2dir = 1, 3
         do i3dir = 1, 3
           write(ab_out,'(3(5x,i2),5x,f16.9)') i1dir,i2dir,i3dir,&
&           d3cart(1,i1dir,i1pert,i2dir,i1pert,i3dir,i1pert)
           if ((blkflg(i1dir,i1pert,i2dir,i1pert,i3dir,i1pert)/=1).or.&
&           (carflg(i1dir,i1pert,i2dir,i1pert,i3dir,i1pert)/=1)) flag = 0
         end do
       end do
     end do

     if (flag == 0) then
       write(message,'(a,a,a,a,a,a)')ch10,&
&       ' dfptnl_doutput: WARNING -',ch10,&
&       '  matrix of third-order energies incomplete,',ch10,&
&       '  non-linear optical coefficients may be wrong.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if

   end if  ! d3e_pert1,d3e_pert2,d3e_pert3

   if (((maxval(d3e_pert1(1:natom))/=0).and.(d3e_pert2(natom+2)/=0).and. &
&   (d3e_pert3(natom+2)/=0)).or.&
   ((maxval(d3e_pert2(1:natom))/=0).and.(d3e_pert1(natom+2)/=0).and. &
&   (d3e_pert3(natom+2)/=0)).or.&
   ((maxval(d3e_pert3(1:natom))/=0).and.(d3e_pert2(natom+2)/=0).and. &
&   (d3e_pert1(natom+2)/=0))) then
!    Perform a check if all relevant elements are available

     flag = 1
     do i1pert = 1, natom
       do i1dir = 1, 3
         do i2dir = 1, 3
           do i3dir = 1, 3
             if ((blkflg(i1dir,i1pert,i2dir,natom+2,i3dir,natom+2) /= 1).or.&
             (blkflg(i1dir,natom+2,i2dir,i1pert,i3dir,natom+2) /= 1).or.&
             (blkflg(i1dir,natom+2,i2dir,natom+2,i3dir,i1pert) /= 1)) flag = 0
             if ((carflg(i1dir,i1pert,i2dir,natom+2,i3dir,natom+2) /= 1).or.&
             (carflg(i1dir,natom+2,i2dir,i1pert,i3dir,natom+2) /= 1).or.&
             (carflg(i1dir,natom+2,i2dir,natom+2,i3dir,i1pert) /= 1)) flag = 0
           end do
         end do
       end do
     end do

     write(ab_out,*)ch10
     write(ab_out,*)' First-order change in the electronic dielectric '
     write(ab_out,*)' susceptibility tensor (Bohr^-1)'
     write(ab_out,*)' induced by an atomic displacement'
     write(ab_out,*)'  atom  displacement'

     do i1pert = 1,natom
       do i1dir = 1,3
         write(ab_out,'(1x,i4,9x,i2,3(3x,f16.9))')i1pert,i1dir,&
&         d3cart(1,i1dir,i1pert,1,natom+2,:,natom+2)
         write(ab_out,'(16x,3(3x,f16.9))')&
&         d3cart(1,i1dir,i1pert,2,natom+2,:,natom+2)
         write(ab_out,'(16x,3(3x,f16.9))')&
&         d3cart(1,i1dir,i1pert,3,natom+2,:,natom+2)
       end do
       write(ab_out,*)
     end do

     if (flag == 0) then
       write(message,'(a,a,a,a,a,a)')ch10,&
&       ' dfptnl_doutput: WARNING -',ch10,&
&       '  matrix of third-order energies incomplete,',ch10,&
&       '  changes in the dielectric susceptibility may be wrong.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if

   end if  ! d3e_pert1,d3e_pert2,d3e_pert3
 end if   ! mpi_enreg%me

 ABI_DEALLOCATE(blkflg)
 ABI_DEALLOCATE(carflg)
 ABI_DEALLOCATE(cg)
 ABI_DEALLOCATE(d3lo)
 ABI_DEALLOCATE(d3cart)
 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(rhor)
 ABI_DEALLOCATE(kneigh)
 ABI_DEALLOCATE(kg_neigh)
 ABI_DEALLOCATE(kptindex)
 ABI_DEALLOCATE(kpt3)
 ABI_DEALLOCATE(mvwtk)
 ABI_DEALLOCATE(pwind)
 ABI_DEALLOCATE(cgindex)
 ABI_DEALLOCATE(d3e_pert1)
 ABI_DEALLOCATE(d3e_pert2)
 ABI_DEALLOCATE(d3e_pert3)
 ABI_DEALLOCATE(rfpert)
 ABI_DEALLOCATE(amass)
 ABI_DEALLOCATE(grxc)
 ABI_DEALLOCATE(kg)
 ABI_DEALLOCATE(kxc)
 ABI_DEALLOCATE(k3xc)
 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(npwarr)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(irrzon)
 ABI_DEALLOCATE(phnons)
 ABI_DEALLOCATE(mpi_enreg%kpt_loc2ibz_sp)
 ABI_DEALLOCATE(mpi_enreg%mkmem)

!Clean the header
 call hdr_free(hdr)

!As the etotal energy has no meaning here, we set it to zero
!(to avoid meaningless side-effects when comparing ouputs...)
 etotal = zero

 call status(0,dtfil%filstat,iexit,level,' exit         ')
 call timab(501,2,tsec)

end subroutine nonlinear
!!***
