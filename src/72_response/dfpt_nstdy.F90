!!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_nstdy
!! NAME
!! dfpt_nstdy
!!
!! FUNCTION
!! This routine compute the non-stationary expression for the
!! second derivative of the total energy, for a whole row of
!! mixed derivatives.
!! Only for norm-conserving pseudopotentials (no PAW)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG, DCA, GMR, MM, AR, MV, MB, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions at k
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF wavefunctions at k,q.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  eigen1(2*mband*mband*nkpt_rbz*nsppol)=array for holding eigenvalues
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  idir=direction of the perturbation
!!  indkpt1(nkpt_rbz)=non-symmetrized indices of the k-points
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  ipert=type of the perturbation
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points in the reduced BZ
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!  mkmem =number of k points treated by this node (GS data)
!!  mk1mem =number of k points treated by this node (RF data)
!!  mpert =maximum number of ipert
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  nattyp(ntypat)= # atoms of each type.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  nfft=(effective) number of FFT grid points (for this proc)
!!  ngfft(1:18)=integer array with FFT box dimensions and other
!!  nkpt=number of k points in the full BZ
!!  nkpt_rbz=number of k points in the reduced BZ for this perturbation
!!  nkxc=second dimension of the kxc array. If /=0, the XC kernel must be computed.
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym1=number of symmetry elements in space group consistent with i perturbation
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band
!!   and k in the reduced Brillouin zone (usually =2)
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhor1(cplex*nfft,nspden)=RF electron density in electrons/bohr**3.
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  ucvol=unit cell volume in bohr**3.
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point in the reduced BZ
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k+q point
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert)=flags for each element of the 2DTE (=1 if computed)
!!  d2bbb(2,3,3,mpert,mband,mband*prtbbb)=band by band decomposition of some
!!                                        second order derivatives
!!  d2lo(2,3,mpert,3,mpert)=local contributions to the 2DTEs
!!  d2nl(2,3,mpert,3,mpert)=non-local contributions to the 2DTEs
!!
!! NOTES
!! Note that the ddk perturbation should not be treated here.
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      appdig,destroy_hamiltonian,dfpt_mkcore,dfpt_mkvxc,dfpt_mkvxc_noncoll
!!      dfpt_nstwf,dfpt_sygra,dfpt_vlocal,dotprod_vn,hdr_skip,init_hamiltonian
!!      load_spin_hamiltonian,mati3inv,timab,wffclose,wffopen,wffreadnpwrec
!!      wffreadskipk,wffreadskiprec,wfk_close,wfk_open_read,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_nstdy(atindx,blkflg,cg,cg1,cplex,dtfil,dtset,d2bbb,d2lo,d2nl,eigen0,eigen1,&
&          gmet,gsqcut,idir,indkpt1,indsy1,ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mkmem,mk1mem,&
&          mpert,mpi_enreg,mpw,mpw1,nattyp,nband_rbz,nfft,ngfft,nkpt,nkpt_rbz,nkxc,&
&          npwarr,npwar1,nspden,nsppol,nsym1,occ_rbz,ph1d,psps,rhor1,rmet,rprimd,&
&          symrc1,ucvol,wtk_rbz,xred,ylm,ylm1,rhor)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_wffile
 use m_wfk
 use m_nctk
 use m_hamiltonian

 use m_io_tools,  only : file_exists
 use m_hdr,       only : hdr_skip
 use m_pawtab,    only : pawtab_type
 use m_paw_ij,    only : paw_ij_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_nstdy'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_53_spacepar
 use interfaces_56_xc
 use interfaces_62_iowfdenpot
 use interfaces_72_response, except_this_one => dfpt_nstdy
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mk1mem,mkmem,mpert,mpw,mpw1,nfft,nkpt,nkpt_rbz,nkxc,nspden,nsppol,nsym1
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx(dtset%natom),indkpt1(nkpt_rbz),indsy1(4,nsym1,dtset%natom)
 integer,intent(in) :: istwfk_rbz(nkpt_rbz),kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem)
 integer,intent(in) :: nattyp(dtset%ntypat),nband_rbz(nkpt_rbz*nsppol),ngfft(18)
 integer,intent(in) :: npwar1(nkpt_rbz),npwarr(nkpt_rbz),symrc1(3,3,nsym1)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert) !vz_i
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*nsppol)
 real(dp),intent(in) :: eigen0(dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: gmet(3,3),kpt_rbz(3,nkpt_rbz)
 real(dp),intent(in) :: kxc(nfft,nkxc),occ_rbz(dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
 real(dp),intent(in) :: rhor1(cplex*nfft,nspden),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xred(3,dtset%natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*dtset%prtbbb)!vz_i
 real(dp),intent(inout) :: d2lo(2,3,mpert,3,mpert),d2nl(2,3,mpert,3,mpert) !vz_i
! optional
 real(dp),optional,intent(in) :: rhor(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: formeig1=1
 integer :: ban2tot,bantot,bdtot_index,ddkcase,iband,icg,icg1,idir1
 integer :: ierr,ifft,ii,ikg,ikg1,ikpt,ikpt_dum,ilm,ipert1,ispden,isppol
 integer :: istwf_k,isym,jj,master,me,n1,n2,n3,n3xccc,n4,n5,n6,nband_dum
 integer :: nband_k,nfftot,npw1_k,npw_k,nskip,nspinor_,option,spaceworld
 integer :: optnc,optxc
 real(dp) :: doti,dotr,wtk_k
 logical :: t_exist
 character(len=500) :: msg
 character(len=fnlen) :: fiwfddk
 type(gs_hamiltonian_type) :: gs_hamkq
!arrays
 integer :: ddkfil(3),ikpt_fbz(3),ikpt_fbz_previous(3),skipddk(3)
 integer,allocatable :: kg1_k(:,:),kg_k(:,:),symrl1(:,:,:)
 real(dp) :: d2nl_elfd(2,3),d2nl_mgfd(2,3),kpoint(3),kpq(3),sumelfd(2),summgfd(2),tsec(2)
 real(dp),allocatable :: buffer1(:),buffer2(:),d2bbb_k(:,:,:,:),d2nl_k(:,:,:)
 real(dp),allocatable :: eig1_k(:),eig_k(:),occ_k(:)
 real(dp) :: rhodummy(0,0)
 real(dp),allocatable :: vpsp1(:),vxc1(:,:),work1(:,:,:),xccc3d1(:),ylm1_k(:,:),ylm_k(:,:)
 type(paw_ij_type) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawtab_type) :: pawtab(dtset%ntypat*psps%usepaw)
 type(wffile_type) :: wffddk(3)
 type(wfk_t) :: ddks(3)

! *********************************************************************

 ABI_UNUSED((/nkpt, ii, ikpt_dum, nband_dum/))

 DBG_ENTER("COLL")

!Not valid for PAW
 if (psps%usepaw==1) then
   msg='This routine cannot be used for PAW (use pawnst3 instead) !'
   MSG_BUG(msg)
 end if

!Keep track of total time spent in dfpt_nstdy
 call timab(101,1,tsec)

!Init parallelism
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt

 master =0

!Zero only portion of nonlocal matrix to be computed here
 d2nl(:,:,1:dtset%natom+2,idir,ipert)=zero

 ABI_ALLOCATE(d2bbb_k,(2,3,dtset%mband,dtset%mband*dtset%prtbbb))
 ABI_ALLOCATE(d2nl_k,(2,3,mpert))
 ABI_ALLOCATE(eig_k,(nsppol*dtset%mband))
 ABI_ALLOCATE(eig1_k,(2*nsppol*dtset%mband**2))
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(kg1_k,(3,mpw1))

!Do not try to open electric field file
 ddkfil(:)=0
!The treatment of homogeneous electric field potential need the existence of d/dk files.
 do idir1=1,3
   ddkcase=idir1+dtset%natom*3
   call appdig(ddkcase,dtfil%fnamewffddk,fiwfddk)

!  Check that ddk file exists
   t_exist = file_exists(fiwfddk)
   if (.not. t_exist) then
     ! Try netcdf file.
     t_exist = file_exists(nctk_ncify(fiwfddk))
     if (t_exist) then
       fiwfddk = nctk_ncify(fiwfddk)
       write(msg,"(3a)")"- File: ",trim(fiwfddk)," does not exist but found netcdf file with similar name."
       call wrtout(std_out,msg,'COLL')
     end if
   end if

   if (t_exist) then
!    Note the use of unit numbers 21, 22 and 23
     ddkfil(idir1)=20+idir1
     write(msg, '(a,a)') '-open ddk wf file :',trim(fiwfddk)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
#ifndef DEV_MG_WFK
     call WffOpen(dtset%iomode,spaceworld,fiwfddk,ierr,wffddk(idir1),master,me,ddkfil(idir1))
#else
     call wfk_open_read(ddks(idir1),fiwfddk,formeig1,dtset%iomode,ddkfil(idir1),spaceworld)
#endif
   end if
 end do

!Update list of computed matrix elements
 if (ipert /= dtset%natom + 1) then
   do ipert1=1,mpert
     do idir1=1,3
       if(ipert1 <= dtset%natom .or. ipert1==dtset%natom+2 &
&       .and. ddkfil(idir1)/=0) then
         blkflg(idir1,ipert1,idir,ipert)=1
       end if
     end do
   end do
 else
   ipert1 = dtset%natom + 1
   do idir1=1,3
!    If was already computed in another run or dataset, or if is to be computed in the present one
     if ((ddkfil(idir1) /= 0).or. (dtset%rfdir(idir1)/=0.and. idir1<=idir) ) then
!      if ((ddkfil(idir1) /= 0).or. (idir1==idir) ) then
       blkflg(idir1,ipert1,idir,ipert)=1
     end if
   end do
 end if

 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)
 nspinor_=dtset%nspinor

!Initialisation of the ddk files
#ifndef DEV_MG_WFK
 do idir1=1,3
   if (ddkfil(idir1)/=0)then
     call hdr_skip(wffddk(idir1),ierr)
   end if
 end do
#endif

 bantot = 0
 ban2tot = 0
 skipddk(:) = 0

!==== Initialize most of the Hamiltonian ====
!1) Allocate all arrays and initialize quantities that do not depend on k and spin.
!2) Perform the setup needed for the non-local factors:
!3) Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 call init_hamiltonian(gs_hamkq,psps,pawtab,dtset%nspinor,nspden,dtset%natom,&
& dtset%typat,xred,nfft,dtset%mgfft,ngfft,rprimd,dtset%nloalg,ph1d=ph1d,&
& use_gpu_cuda=dtset%use_gpu_cuda)

!LOOP OVER SPINS
 bdtot_index=0
 icg=0;icg1=0
 do isppol=1,nsppol

!  In case isppol = 2, skip the records that correspond to isppol = 1 and that have not been read
#ifndef DEV_MG_WFK
   if (isppol == 2) then
     do idir1 = 1, 3
       if ((ddkfil(idir1)/=0).and.(skipddk(idir1) < nkpt)) then
         do ikpt = 1, (nkpt - skipddk(idir1))
           call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw_k,nspinor_,wffddk(idir1))
           call WffReadSkipRec(ierr,1,wffddk(idir1))
           do iband = 1, nband_k
             call WffReadSkipRec(ierr,2,wffddk(idir1))
           end do
         end do
       end if
     end do
   end if
#endif

   ikg=0;ikg1=0

   ikpt_fbz(1:3)=0

!  Continue to initialize the Hamiltonian
   call load_spin_hamiltonian(gs_hamkq,isppol,paw_ij=paw_ij)

!  BIG FAT k POINT LOOP
   do ikpt=1,nkpt_rbz

     nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
     istwf_k=istwfk_rbz(ikpt)
     npw_k=npwarr(ikpt)
     npw1_k=npwar1(ikpt)

     eig_k(1:nband_k) = eigen0(1+bantot:nband_k+bantot)
     eig1_k(1:2*nband_k**2) = eigen1(1+ban2tot:2*nband_k**2+ban2tot)
     bantot = bantot + nband_k
     ban2tot = ban2tot + 2*nband_k**2

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k
!      The wavefunction blocks for ddk file is skipped elsewhere in the loop
!      Skip the rest of the k-point loop
       cycle
     end if

     ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
     ABI_ALLOCATE(ylm1_k,(npw1_k,psps%mpsang*psps%mpsang*psps%useylm))

!    In case of electric field pert1, read ddk wfs file
!    Note that the symmetries are not used for ddk, so read each k point
!    Also take into account implicitely the parallelism over k points

     do idir1=1,3
       if (ddkfil(idir1)/=0) then
!        Must select the corresponding k point in the full set of k points
!        used in the ddk file : compute the number of k points to skip
         ikpt_fbz_previous(idir1)=ikpt_fbz(idir1)
         ikpt_fbz(idir1)=indkpt1(ikpt)

         nskip=ikpt_fbz(idir1)-ikpt_fbz_previous(idir1)-1
         skipddk(idir1) = skipddk(idir1) + 1 + nskip
#ifndef DEV_MG_WFK
         if(nskip/=0)then
           do ikpt_dum=1+ikpt_fbz_previous(idir1),ikpt_fbz(idir1)-1
             nband_dum=dtset%nband(ikpt_dum+(isppol-1)*nkpt)
!            Skip the records whose information is not needed (in case of parallelism)
             call WffReadSkipK(1,0,ikpt_dum,isppol,mpi_enreg,wffddk(idir1))
           end do
         end if
#else
         ii = wfk_findk(ddks(idir1), kpt_rbz(:, ikpt))
         ABI_CHECK(ii == indkpt1(ikpt),  "ii !=  indkpt1")
#endif
       end if
     end do

     ABI_ALLOCATE(occ_k,(nband_k))
     occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)
     kpoint(:)=kpt_rbz(:,ikpt)
     kpq(:)=kpoint(:)+dtset%qptn(:)
     wtk_k=wtk_rbz(ikpt)
     d2nl_k(:,:,:)=zero
     if(dtset%prtbbb==1)d2bbb_k(:,:,:,:)=zero

!    Get plane-wave vectors and related data at k
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if

!    Get plane-wave vectors and related data at k+q
     kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm1_k(1:npw1_k,ilm)=ylm1(1+ikg1:npw1_k+ikg1,ilm)
       end do
     end if

!    Compute the eigenvalues, wavefunction,
!    contributions to kinetic energy, nonlocal energy, forces,
!    and update of rhor1 to this k-point and this spin polarization.
!    Note that dfpt_nstwf is called with kpoint, while kpt is used inside dfpt_vtowfk
     call dfpt_nstwf(cg,cg1,ddkfil,dtset,d2bbb_k,d2nl_k,eig_k,eig1_k,gs_hamkq,&
&     icg,icg1,idir,ikpt,ipert,isppol,istwf_k,kg_k,kg1_k,kpoint,kpq,mkmem,mk1mem,mpert,&
&     mpi_enreg,mpw,mpw1,nband_k,npw_k,npw1_k,nsppol,&
&     occ_k,psps,rmet,wffddk,ddks,wtk_k,ylm_k,ylm1_k)

     d2nl(:,:,:,idir,ipert)=d2nl(:,:,:,idir,ipert)+d2nl_k(:,:,:)
     if(dtset%prtbbb==1)d2bbb(:,:,idir,ipert,:,:)=d2bbb(:,:,idir,ipert,:,:)+d2bbb_k(:,:,:,:)

!    Keep track of total number of bands
     bdtot_index=bdtot_index+nband_k

!    Shift arrays memory
     if (mkmem/=0) then
       icg=icg+npw_k*dtset%nspinor*nband_k
       ikg=ikg+npw_k
     end if
     if (mk1mem/=0) then
       icg1=icg1+npw1_k*dtset%nspinor*nband_k
       ikg1=ikg1+npw1_k
     end if

     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylm1_k)

!    End big k point loop
   end do

!  End loop over spins
 end do

!if(xmpi_paral==1)then
!call timab(161,1,tsec)
!call wrtout(std_out,' dfpt_nstdy: loop on k-points and spins done in parallel','COLL')
!call xmpi_barrier(spaceworld)
!call timab(161,2,tsec)
!end if

 call destroy_hamiltonian(gs_hamkq)

!Treat fixed occupation numbers (as in vtorho)
 if(xmpi_paral==1)then
   ABI_ALLOCATE(buffer1,(2*3*mpert))
   ABI_ALLOCATE(buffer2,(2*3*mpert))
!  Pack d2nl
   buffer1(1:2*3*mpert)=reshape(d2nl(:,:,:,idir,ipert),(/2*3*mpert/))
!  Build sum of everything
   call timab(48,1,tsec)
   call xmpi_sum(buffer1,buffer2,2*3*mpert,spaceworld,ierr)
   call timab(48,2,tsec)
!  Unpack the final result
   d2nl(:,:,:,idir,ipert)=reshape(buffer2(:),(/2,3,mpert/))
   ABI_DEALLOCATE(buffer1)
   ABI_DEALLOCATE(buffer2)
   if(dtset%prtbbb==1)then
     ABI_ALLOCATE(buffer1,(2*3*dtset%mband*dtset%mband))
     ABI_ALLOCATE(buffer2,(2*3*dtset%mband*dtset%mband))
!    Pack d2bbb
     buffer1(1:2*3*dtset%mband*dtset%mband)=reshape(d2bbb(:,:,idir,ipert,:,:),(/2*3*dtset%mband*dtset%mband/))
!    Build sum of everything
     call timab(48,1,tsec)
     call xmpi_sum(buffer1,buffer2,2*3*dtset%mband*dtset%mband,spaceworld,ierr)
     call timab(48,2,tsec)
!    Unpack the final result
     d2bbb(:,:,idir,ipert,:,:)=reshape(buffer2(:),(/2,3,dtset%mband,dtset%mband/))
     ABI_DEALLOCATE(buffer1)
     ABI_DEALLOCATE(buffer2)
   end if
 end if ! xmpi_paral==1

!In the case of the strain perturbation time-reversal symmetry will always
!be true so imaginary part of d2nl will be must be set to zero here since
!the symmetry-reduced kpt set will leave a non-zero imaginary part.
 if(ipert==dtset%natom+3 .or. ipert==dtset%natom+4) d2nl(2,:,:,idir,ipert)=zero

!In case of electric field ipert1, close the ddk wf files
 do idir1=1,3
   if (ddkfil(idir1)/=0)then
#ifndef DEV_MG_WFK
     call WffClose(wffddk(idir1),ierr)
#else
     call wfk_close(ddks(idir1))
#endif
   end if
 end do

!Symmetrize the non-local contributions,
!as was needed for the forces in a ground-state calculation
!However, here the quantity is complex, and there are phases !

!Do the transform
 ABI_ALLOCATE(work1,(2,3,dtset%natom))
 do ipert1=1,dtset%natom
   do idir1=1,3
     work1(1,idir1,ipert1)=d2nl(1,idir1,ipert1,idir,ipert)
     work1(2,idir1,ipert1)=d2nl(2,idir1,ipert1,idir,ipert)
   end do
 end do
 call dfpt_sygra(dtset%natom,d2nl(:,:,:,idir,ipert),work1,indsy1,ipert,nsym1,dtset%qptn,symrc1)
 ABI_DEALLOCATE(work1)

!Must also symmetrize the electric/magnetic field perturbation response !
!(XG 000803 This was not implemented until now)
 if(sum(ddkfil(:))/=0)then
!  Get the symmetry matrices in terms of real space basis
   ABI_ALLOCATE(symrl1,(3,3,nsym1))
   do isym=1,nsym1
     call mati3inv(symrc1(:,:,isym),symrl1(:,:,isym))
   end do
!  There should not be any imaginary part, but stay general (for debugging)
   d2nl_elfd(:,:)=d2nl(:,:,dtset%natom+2,idir,ipert)
   do ii=1,3
     sumelfd(:)=zero
     summgfd(:)=zero
     do isym=1,nsym1
       do jj=1,3
         if(symrl1(ii,jj,isym)/=0)then
           if(ddkfil(jj)==0)then
             blkflg(ii,dtset%natom+2,idir,ipert)=0
           end if
         end if
       end do
       sumelfd(:)=sumelfd(:)+dble(symrl1(ii,1,isym))*d2nl_elfd(:,1)+&
&       dble(symrl1(ii,2,isym))*d2nl_elfd(:,2)+&
&       dble(symrl1(ii,3,isym))*d2nl_elfd(:,3)
       summgfd(:)=summgfd(:)+dble(symrl1(ii,1,isym))*d2nl_mgfd(:,1)+&
&       dble(symrl1(ii,2,isym))*d2nl_mgfd(:,2)+&
&       dble(symrl1(ii,3,isym))*d2nl_mgfd(:,3)
     end do
     d2nl(:,ii,dtset%natom+2,idir,ipert)=sumelfd(:)/dble(nsym1)
   end do

   if ((dtset%prtbbb==1).and.(ipert<=dtset%natom)) then
     do iband = 1,dtset%mband
       d2nl_elfd(:,:)=d2bbb(:,:,idir,ipert,iband,iband)
       do ii=1,3
         sumelfd(:)=zero
         do isym=1,nsym1
           sumelfd(:)=sumelfd(:)+dble(symrl1(ii,1,isym))*d2nl_elfd(:,1)+&
&           dble(symrl1(ii,2,isym))*d2nl_elfd(:,2)+&
&           dble(symrl1(ii,3,isym))*d2nl_elfd(:,3)
         end do
         d2bbb(:,ii,idir,ipert,iband,iband)=sumelfd(:)/dble(nsym1)
       end do
     end do  !iband
   end if

   ABI_DEALLOCATE(symrl1)
 end if

!----------------------------------------------------------------------------
!Now, treat the local contribution

 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 ABI_ALLOCATE(vpsp1,(cplex*nfft))
 if (ipert /= dtset%natom + 1) then
   n3xccc=0;if(psps%n1xccc/=0) n3xccc=nfft
   ABI_ALLOCATE(xccc3d1,(cplex*n3xccc))
   ABI_ALLOCATE(vxc1,(cplex*nfft,nspden))

   do ipert1=1,mpert
     do idir1=1,3
       if(ipert1 <= dtset%natom)then

!        Get first-order local potential and first-order pseudo core density
         call dfpt_vlocal(atindx,cplex,gmet,gsqcut,idir1,ipert1,mpi_enreg,psps%mqgrid_ff,dtset%natom,&
&         nattyp,nfft,ngfft,dtset%ntypat,n1,n2,n3,dtset%paral_kgb,ph1d,psps%qgrid_ff,&
&         dtset%qptn,ucvol,psps%vlspl,vpsp1,xred)
         if(psps%n1xccc/=0)then
           call dfpt_mkcore(cplex,idir1,ipert1,dtset%natom,dtset%ntypat,n1,psps%n1xccc,&
&           n2,n3,dtset%qptn,rprimd,dtset%typat,ucvol,psps%xcccrc,psps%xccc1d,xccc3d1,xred)
         end if

!        Get first-order exchange-correlation potential (core-correction contribution only !)
         if(psps%n1xccc/=0)then
           option=0
!FR EB non-collinear magnetism
! the second nkxc should be nkxc_cur (see 67_common/nres2vres.F90)
           if (nspden==4.and.present(rhor)) then
             optnc=1
             optxc=1
             call dfpt_mkvxc_noncoll(cplex,dtset%ixc,kxc,mpi_enreg,nfft,ngfft,rhodummy,0,rhodummy,0,&
&             nkxc,nkxc,nspden,n3xccc,optnc,option,optxc,dtset%paral_kgb,dtset%qptn,rhodummy,rhodummy,&
&             rprimd,0,vxc1,xccc3d1)
           else
             call dfpt_mkvxc(cplex,dtset%ixc,kxc,mpi_enreg,nfft,ngfft,rhodummy,0,rhodummy,0,&
&             nkxc,nspden,n3xccc,option,dtset%paral_kgb,dtset%qptn,rhodummy,&
&             rprimd,0,vxc1,xccc3d1)
           end if
         else
           vxc1(:,:)=zero
         end if

!        Norm-conserving pseudpopotential case:
!        Combines density j2 with local potential j1 (vpsp1 and vxc1)
!        XG030514 : this is a first possible coding, however, each dotprod contains
!        a parallel section (reduction), so it is better to use only one dotprod ...
!        call dotprod_vn(cplex,rhor1,dr_psp1,di_psp1,mpi_enreg,nfft,nfftot,1,2,vpsp1,ucvol)
!        call dotprod_vn(cplex,rhor1,dr_xc1,di_xc1,mpi_enreg,nfft,nfftot,nspden,2,vxc1,ucvol)
!        dotr=dr_psp1+dr_xc1;doti=di_psp1+di_xc1... but then, one needs to overload vxc1
         do ispden=1,min(nspden,2)
           do ifft=1,cplex*nfft
             vxc1(ifft,ispden)=vxc1(ifft,ispden)+vpsp1(ifft)
           end do
         end do
         call dotprod_vn(cplex,rhor1,dotr,doti,nfft,nfftot,nspden,2,vxc1,ucvol)

!        MVeithen 021212 : in case ipert = 2, these lines compute the local part
!        of the Born effective charges from phonon and electric
!        field type perturbations, see eq. 43 of
!        X. Gonze and C. Lee, PRB 55, 10355 (1997)
!        The minus sign is due to the fact that the effective charges
!        are minus the second derivatives of the energy
         if (ipert == dtset%natom+2) then
           d2lo(1,idir1,ipert1,idir,ipert)=-dotr
           d2lo(2,idir1,ipert1,idir,ipert)=-doti
         else
           d2lo(1,idir1,ipert1,idir,ipert)=dotr
           d2lo(2,idir1,ipert1,idir,ipert)=doti
         end if
!        Endif ipert1<=natom
       end if
     end do
   end do

   ABI_DEALLOCATE(vxc1)
   ABI_DEALLOCATE(xccc3d1)

 end if ! ipert /= natom +1

 ABI_DEALLOCATE(d2bbb_k)
 ABI_DEALLOCATE(d2nl_k)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(kg1_k)
 ABI_DEALLOCATE(vpsp1)
 ABI_DEALLOCATE(eig_k)
 ABI_DEALLOCATE(eig1_k)

 call timab(101,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_nstdy
!!***
