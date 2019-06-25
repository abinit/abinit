!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_phpi
!! NAME
!!
!! FUNCTION
!!  Tools for the computation of phonon self-energy.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (GKA)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_phpi

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_abicore
 use m_xmpi
 use m_errors
 use m_ifc
 use m_ebands
 use iso_c_binding
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_wfk
 use m_ddb
 use m_dvdb
 use m_fft
 use m_hamiltonian
 use m_pawcprj

 use m_time,            only : cwtime
 use m_fstrings,        only : sjoin, itoa, ftoa, ktoa, ltoa, strcat
 use m_io_tools,        only : iomode_from_fname
 use m_cgtools,         only : dotprod_g
 use m_kg,              only : getph
 use m_fftcore,         only : get_kg
 use m_crystal,         only : crystal_t
 use m_bz_mesh,         only : findqg0
 use m_wfd,             only : wfd_init, wfd_t
 use m_pawang,          only : pawang_type
 use m_pawrad,          only : pawrad_type
 use m_pawtab,          only : pawtab_type
 use m_pawfgr,          only : pawfgr_type
 use m_getgh1c,         only : getgh1c, rf_transgrid_and_pack, getgh1c_setup

 implicit none

 private
!!***

 public :: eph_phpi


contains  !=================================================================================
!!***

!!****f* m_phpi/eph_phpi
!! NAME
!!  eph_phpi
!!
!! FUNCTION
!!  Compute phonon self-energy.
!!
!! INPUTS
!! wk0_path=String with the path to the GS unperturbed WFK file.
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!! ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!! pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!! pawang<pawang_type)>=PAW angular mesh and related data.
!! pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! comm=MPI communicator.
!!
!! OUTPUT
!!
!! PARENTS
!!      eph
!!
!! NOTES
!!
!! CHILDREN
!!
!! SOURCE

subroutine eph_phpi(wfk0_path,wfq_path,dtfil,ngfft,ngfftf,dtset,cryst,ebands_k,ebands_kq,dvdb,ifc,&
                       pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path, wfq_path
 integer,intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands_k, ebands_kq
 type(dvdb_t),intent(inout) :: dvdb
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(ifc_type),intent(in) :: ifc
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_getgh1c=1,berryopt0=0
 integer,parameter :: useylmgr1=0,master=0
 integer :: my_rank,nproc,iomode,mband,mband_kq,my_minb,my_maxb,nsppol,nkpt,nkpt_kq,idir,ipert
 integer :: cplex,db_iqpt,natom,natom3,ipc,nspinor,onpw,imode
 integer :: ib1,ib2
 integer :: ik,ikq
 integer :: spin,istwf_k,istwf_kq,npw_k,npw_kq
 integer :: mpw,my_mpw,ierr,my_kstart,my_kstop,cnt
 integer :: n1,n2,n3,n4,n5,n6,nspden
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg,nkpg1
 real(dp) :: cpu,wall,gflops
 real(dp) :: ecut,eshift,eig0nk,eig0mkq,dotr,doti
 real(dp) :: eta,f_nk,f_mkq,omega,wtk,gkk2,term1,term2
 logical :: i_am_master,gen_eigenpb
 type(wfd_t) :: wfd_k,wfd_kq
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 character(len=500) :: msg
!arrays
 integer :: g0_k(3)
 integer,allocatable :: kg_k(:,:),kg_kq(:,:),gtmp(:,:),nband(:,:),nband_kq(:,:),blkflg(:,:), wfd_istwfk(:)
 real(dp) :: kk(3),kq(3),qpt(3),phfrq(3*cryst%natom)
 real(dp) :: displ_cart(2,3,cryst%natom,3*cryst%natom),displ_red(2,3,cryst%natom,3*cryst%natom)
 real(dp) :: Pi_ph(3*cryst%natom)
 real(dp),allocatable :: grad_berry(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:),gkk(:,:,:,:,:), gkk_m(:,:,:)
 real(dp),allocatable :: bras_kq(:,:,:),kets_k(:,:,:),h1kets_kq(:,:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:),gvnlx1(:,:)
 real(dp),allocatable ::  gs1c(:,:)
 logical,allocatable :: bks_mask(:,:,:),bks_mask_kq(:,:,:),keep_ur(:,:,:),keep_ur_kq(:,:,:)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:) !natom,nspinor*usecprj)

!************************************************************************

 write(msg, '(3a)') ch10, "Computation of the real part of the phonon self-energy", ch10
 call wrtout(ab_out, msg, "COLL", do_flush=.True.)
 call wrtout(std_out, msg, "COLL", do_flush=.True.)

 if (psps%usepaw == 1) then
   MSG_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 i_am_master = (my_rank == master)

 ! Copy important dimensions
 natom = cryst%natom
 natom3 = 3 * natom
 nsppol = ebands_k%nsppol
 nspinor = ebands_k%nspinor
 nspden = dtset%nspden
 nkpt = ebands_k%nkpt
 mband = ebands_k%mband
 nkpt_kq = ebands_kq%nkpt
 mband_kq = ebands_kq%mband
 ecut = dtset%ecut

! GKA TODO: Make sure there is a single q-point present.
 qpt = dtset%qptn(:)

 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)

 ! Initialize the wave function descriptors.
 ! For the time being, no memory distribution, each node has the full set of states.
 my_minb = 1; my_maxb = mband

 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(bks_mask,(mband, nkpt, nsppol))
 ABI_MALLOC(keep_ur,(mband, nkpt ,nsppol))
 nband=mband; bks_mask=.False.; keep_ur=.False.

 ABI_MALLOC(nband_kq, (nkpt_kq, nsppol))
 ABI_MALLOC(bks_mask_kq,(mband_kq, nkpt_kq, nsppol))
 ABI_MALLOC(keep_ur_kq,(mband_kq, nkpt_kq ,nsppol))
 nband_kq=mband_kq; bks_mask_kq=.False.; keep_ur_kq=.False.

 ! Distribute the k-points over the processors
 call xmpi_split_work(nkpt,comm,my_kstart,my_kstop)
 do ik=1,nkpt
 if (.not. ((ik .ge. my_kstart) .and. (ik .le. my_kstop))) cycle

   kk = ebands_k%kptns(:,ik)
   kq = kk + qpt
   call findqg0(ikq,g0_k,kq,nkpt_kq,ebands_kq%kptns(:,:),(/1,1,1/))  ! Find the index of the k+q point

   bks_mask(:,ik,:) = .True.
   bks_mask_kq(:,ikq,:) = .True.
 end do

 ! Initialize the wavefunction descriptors

 ! Impose istwfk=1 for all k points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkpt))
 wfd_istwfk = 1

 call wfd_init(wfd_k,cryst,pawtab,psps,keep_ur,mband,nband,nkpt,nsppol,bks_mask,&
   nspden,nspinor,ecut,dtset%ecutsm,dtset%dilatmx,wfd_istwfk,ebands_k%kptns,ngfft,&
   dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm)
 ABI_FREE(wfd_istwfk)

 call wfd_k%print(header="Wavefunctions on the k-points grid",mode_paral='PERS')

 ABI_MALLOC(wfd_istwfk, (nkpt_kq))
 wfd_istwfk = 1

 call wfd_init(wfd_kq,cryst,pawtab,psps,keep_ur_kq,mband_kq,nband_kq,nkpt_kq,nsppol,bks_mask_kq,&
   nspden,nspinor,ecut,dtset%ecutsm,dtset%dilatmx,wfd_istwfk,ebands_kq%kptns,ngfft,&
   dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm)

 ABI_FREE(wfd_istwfk)

 call wfd_kq%print(header="Wavefunctions on the q-shifted k-points grid",mode_paral='PERS')

 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(nband_kq)
 ABI_FREE(bks_mask_kq)
 ABI_FREE(keep_ur_kq)

 ! Read wafefunctions on the k-points grid and q-shifted k-points grid.
 iomode = iomode_from_fname(wfk0_path)
 call wfd_k%read_wfk(wfk0_path,iomode)
 if (.False.) call wfd_k%test_ortho(cryst,pawtab,unit=std_out,mode_paral="PERS")

 iomode = iomode_from_fname(wfq_path)
 call wfd_kq%read_wfk(wfq_path,iomode)
 if (.False.) call wfd_kq%test_ortho(cryst,pawtab,unit=std_out,mode_paral="PERS")

 ! ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2,3*(2*mgfft+1)*natom))
 call getph(cryst%atindx,natom,n1,n2,n3,ph1d,cryst%xred)

 ! Find the appropriate value of mpw
 mpw = 0; cnt=0
 do spin=1,nsppol
   do ik=1,nkpt
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
     kk = ebands_k%kptns(:,ik)
     call get_kg(kk,1,ecut,cryst%gmet,onpw,gtmp)
     ABI_FREE(gtmp)
     mpw = max(mpw, onpw)
   end do
 end do
 cnt=0
 do spin=1,nsppol
   do ikq=1,nkpt_kq
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
     kq = ebands_kq%kptns(:,ikq)
     call get_kg(kq,1,ecut,cryst%gmet,onpw,gtmp)
     ABI_FREE(gtmp)
     mpw = max(mpw, onpw)
   end do
 end do
 my_mpw = mpw; call xmpi_max(my_mpw, mpw, comm, ierr)

 ! Allow PW-arrays dimensioned with mpw
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(kg_kq, (3, mpw))

 ! Spherical Harmonics for useylm==1.
 ABI_MALLOC(ylm_k,(mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylm_kq,(mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylmgr_kq,(mpw, 3, psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))

 ! TODO FOR PAW
 usecprj = 0
 ABI_DT_MALLOC(cwaveprj0, (natom, nspinor*usecprj))

 ! Prepare call to getgh1c
 usevnl = 0
 optlocal = 1  ! local part of H^(1) is computed in gh1c=<G|H^(1)|C>
 optnl = 2     ! non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
 opt_gvnlx1 = 0 ! gvnlx1 is output
 ABI_MALLOC(gvnlx1, (2,usevnl))
 ABI_MALLOC(grad_berry, (2,nspinor*(berryopt0/4)))

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 !1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 !2) Perform the setup needed for the non-local factors:
 !* Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 !* PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_hamkq,psps,pawtab,nspinor,NSPPOL,nspden,natom,&
&  dtset%typat,cryst%xred,nfft,mgfft,ngfft,cryst%rprimd,dtset%nloalg,&
&  comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
&  usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

 ! Allocate vlocal. Note nvloc
 ! I set vlocal to huge to trigger possible bugs (DFPT routines should not access the data)
 ABI_MALLOC(vlocal,(n4,n5,n6,gs_hamkq%nvloc))
 vlocal = huge(one)

 ! Allocate work space arrays.
 ABI_MALLOC(blkflg, (natom3,natom3))
 ABI_CALLOC(dummy_vtrial, (nfftf,nspden))

 call cwtime(cpu,wall,gflops,"start")

 ! Find the index of the q-point in the DVDB.
 db_iqpt = dvdb%findq(qpt)

 if (db_iqpt /= -1) then
   if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Found: ",ktoa(qpt)," in DVDB with index ",itoa(db_iqpt)))
   ! Read or reconstruct the dvscf potentials for all 3*natom perturbations.
   ! This call allocates v1scf(cplex, nfftf, nspden, 3*natom))
   call dvdb%readsym_allv1(db_iqpt, cplex, nfftf, ngfftf, v1scf, comm)
 else
   MSG_ERROR(sjoin("Could not find symmetric of q-point:", ktoa(qpt), "in DVDB"))
 end if

 ! Allocate vlocal1 with correct cplex. Note nvloc
 ABI_MALLOC_OR_DIE(vlocal1,(cplex*n4,n5,n6,gs_hamkq%nvloc,natom3), ierr)

 ! Allocate el-ph coupling matrix elements
 ABI_MALLOC(gkk, (2, mband_kq, mband, natom, 3))
 ABI_MALLOC(gkk_m, (2, mband_kq, mband))

 ! Compute displacement vectors and phonon frequencies
 call ifc_fourq(ifc,cryst,qpt,phfrq,displ_cart, out_displ_red=displ_red)

 ! Broadening parameter
 if (dtset%elph2_imagden .gt. tol12) then
   eta = dtset%elph2_imagden
 else
   eta = 0.0001_dp
 end if

 ! Kpoints weights (not using symmetries at the moment)
 wtk = 1.0 / nkpt

 ! Initialize phonon self-energy
 Pi_ph = zero


 ! Examine the symmetries of the q wavevector
 ! call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=dtset%prtvol)

 ! ----------------------------------------------------------------------------------------------- !
 ! Begin loop over states
 ! ----------------------------------------------------------------------------------------------- !
 do spin=1,nsppol

   ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
   do ipc=1,natom3
     call rf_transgrid_and_pack(spin,nspden,psps%usepaw,cplex,nfftf,nfft,ngfft,gs_hamkq%nvloc,&
             pawfgr,mpi_enreg,dummy_vtrial,v1scf(:,:,:,ipc),vlocal,vlocal1(:,:,:,:,ipc))
   end do

   ! Continue to initialize the Hamiltonian
   call load_spin_hamiltonian(gs_hamkq,spin,vlocal=vlocal,with_nonlocal=.true.)

   do ik=1,nkpt

     ! Only do a subset a k-points
     if (.not. ((ik .ge. my_kstart) .and. (ik .le. my_kstop))) cycle

     ! Allocate workspace for wavefunctions. Make npw larger than expected.
     ABI_MALLOC(bras_kq, (2, mpw*nspinor, mband))
     ABI_MALLOC(kets_k, (2, mpw*nspinor, mband))
     ABI_MALLOC(h1kets_kq, (2, mpw*nspinor, mband))

     kk = ebands_k%kptns(:,ik)

     kq = kk + qpt
     call findqg0(ikq,g0_k,kq,nkpt_kq,ebands_kq%kptns(:,:),(/1,1,1/))  ! Find the index of the k+q point


     ! Copy u_k(G)
     istwf_k = wfd_k%istwfk(ik); npw_k = wfd_k%npwarr(ik)
     ABI_CHECK(mpw >= npw_k, "mpw < npw_k")
     kg_k(:,1:npw_k) = wfd_k%kdata(ik)%kg_k
     do ib2=1,mband
       call wfd_k%copy_cg(ib2, ik, spin, kets_k(1,1,ib2))
     end do

     ! Copy u_kq(G)
     istwf_kq = wfd_kq%istwfk(ikq); npw_kq = wfd_kq%npwarr(ikq)
     ABI_CHECK(mpw >= npw_kq, "mpw < npw_kq")
     kg_kq(:,1:npw_kq) = wfd_kq%kdata(ikq)%kg_k
     do ib1=1,mband_kq
       call wfd_kq%copy_cg(ib1, ikq, spin, bras_kq(1,1,ib1))
     end do

     ! if PAW, one has to solve a generalized eigenproblem
     ! BE careful here because I will need sij_opt==-1
     gen_eigenpb = (psps%usepaw==1)
     sij_opt = 0; if (gen_eigenpb) sij_opt = 1
     ABI_MALLOC(gs1c, (2,npw_kq*nspinor*((sij_opt+1)/2)))

     gkk = zero

     ! Loop over all 3*natom perturbations.
     do ipc=1,natom3
       idir = mod(ipc-1, 3) + 1
       ipert = (ipc - idir) / 3 + 1

       !write(msg, '(a,2i4)')  "Treating ipert, idir = ", ipert, idir
       !call wrtout(std_out, msg, "COLL", do_flush=.True.)

       ! Prepare application of the NL part.
       call init_rf_hamiltonian(cplex,gs_hamkq,ipert,rf_hamkq,has_e1kbsc=.true.)
           !&paw_ij1=paw_ij1,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
           !&mpi_spintab=mpi_enreg%my_isppoltab)
       call load_spin_rf_hamiltonian(rf_hamkq,spin,vlocal1=vlocal1(:,:,:,:,ipc),with_nonlocal=.true.)

       ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
       call getgh1c_setup(gs_hamkq,rf_hamkq,dtset,psps,kk,kq,idir,ipert,&                   ! In
         cryst%natom,cryst%rmet,cryst%gprimd,cryst%gmet,istwf_k,&                           ! In
         npw_k,npw_kq,useylmgr1,kg_k,ylm_k,kg_kq,ylm_kq,ylmgr_kq,&                          ! In
         dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1)                      ! Out

       ! Calculate dvscf * psi_k, results stored in h1kets_kq on the k+q sphere.
       ! Compute H(1) applied to GS wavefunction Psi(0)
       do ib2=1,mband
         eig0nk = ebands_k%eig(ib2,ik,spin)
         ! Use scissor shift on 0-order eigenvalue
         eshift = eig0nk - dtset%dfpt_sciss

         call getgh1c(berryopt0,kets_k(:,:,ib2),cwaveprj0,h1kets_kq(:,:,ib2),&
&                     grad_berry,gs1c,gs_hamkq,gvnlx1,idir,ipert,eshift,mpi_enreg,optlocal,&
&                     optnl,opt_gvnlx1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)
       end do

       ABI_FREE(kinpw1)
       ABI_FREE(kpg1_k)
       ABI_FREE(kpg_k)
       ABI_FREE(dkinpw)
       ABI_FREE(ffnlk)
       ABI_FREE(ffnl1)
       ABI_FREE(ph3d)

       if (allocated(gs1c)) then
         ABI_FREE(gs1c)
       end if

       if (allocated(ph3d1)) then
         ABI_FREE(ph3d1)
       end if

       ! Calculate elphmat(j,i) = <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation.
       !The array eig1_k contains:
       !
       ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)>                           (NC psps)
       ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(band,k)^(0)> (PAW)
       do ib2=1,mband
         do ib1=1,mband_kq
           call dotprod_g(dotr,doti,istwf_kq,npw_kq*nspinor,2,bras_kq(1,1,ib1),h1kets_kq(1,1,ib2),&
             mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

           gkk(1,ib1,ib2,ipert,idir) = dotr
           gkk(2,ib1,ib2,ipert,idir) =  doti
         end do
       end do

     end do  ! ipc

     ! Loop over 3*natom phonon branches.
     do imode=1,natom3

       omega = phfrq(imode)

       ! Do not compute Pi for negative or too small frequencies
       if (omega .lt. tol6) cycle

       gkk_m = zero

       ! Transform the gkk from atom,cart basis to mode basis
       do idir=1,3
         do ipert=1,natom
            gkk_m(1,:,:) = gkk_m(1,:,:) &
&                          + gkk(1,:,:,ipert,idir) * displ_red(1,idir,ipert,imode) &
&                          - gkk(2,:,:,ipert,idir) * displ_red(2,idir,ipert,imode)
            gkk_m(2,:,:) = gkk_m(2,:,:) &
&                          + gkk(1,:,:,ipert,idir) * displ_red(2,idir,ipert,imode) &
&                          + gkk(2,:,:,ipert,idir) * displ_red(1,idir,ipert,imode)
         end do
       end do

       gkk_m = gkk_m / sqrt(two * omega)

       ! sum contribution to phonon self-energy
       do ib2=1,mband
         do ib1=1,mband_kq

           f_nk = ebands_k%occ(ib2,ik,spin)
           f_mkq = ebands_kq%occ(ib1,ikq,spin)

           if (abs(f_mkq - f_nk) .le. tol12) cycle

           eig0nk = ebands_k%eig(ib2,ik,spin)
           eig0mkq = ebands_kq%eig(ib1,ikq,spin)

           gkk2 = gkk_m(1,ib1,ib2) ** 2 + gkk_m(2,ib1,ib2) ** 2

           term1 = (f_mkq - f_nk) * (eig0mkq - eig0nk - omega) / ((eig0mkq - eig0nk - omega) ** 2 + eta ** 2)
           term2 = (f_mkq - f_nk) * (eig0mkq - eig0nk        ) / ((eig0mkq - eig0nk        ) ** 2 + eta ** 2)

           Pi_ph(imode) = Pi_ph(imode) + wtk * gkk2 * (term1 - term2)

         end do
       end do

     end do  ! imode

     ABI_FREE(bras_kq)
     ABI_FREE(kets_k)
     ABI_FREE(h1kets_kq)
   end do ! ikfs

   call destroy_rf_hamiltonian(rf_hamkq)
 end do ! spin

 ! Gather the k-points computed by all processes
 call xmpi_sum_master(Pi_ph,master,comm,ierr)

 ! Output the results
 if (i_am_master) then
   call out_phpi(ab_out, Pi_ph, phfrq, qpt, natom3)
   call out_phpi(std_out, Pi_ph, phfrq, qpt, natom3)
 end if

#ifdef HAVE_NETCDF
 if (i_am_master) call out_phpi_nc(dtfil, cryst, Pi_ph, phfrq, qpt, natom3)
#endif

 ! Free memory
 call cwtime(cpu,wall,gflops,"stop")

 write(msg, '(3a)') "Computation of the real part of the phonon self-energy completed", ch10, &
&                   "--------------------------------------------------------------------------------"
 call wrtout(ab_out, msg, "COLL", do_flush=.True.)
 call wrtout(std_out, msg, "COLL", do_flush=.True.)

 ! Free memory
 ABI_FREE(gkk)
 ABI_FREE(gkk_m)
 ABI_FREE(v1scf)
 ABI_FREE(vlocal1)
 ABI_FREE(gvnlx1)
 ABI_FREE(grad_berry)
 ABI_FREE(dummy_vtrial)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(kg_k)
 ABI_FREE(kg_kq)
 ABI_FREE(ylm_k)
 ABI_FREE(ylm_kq)
 ABI_FREE(ylmgr_kq)
 ABI_FREE(blkflg)

 call destroy_hamiltonian(gs_hamkq)
 call wfd_k%free()
 call wfd_kq%free()

 call pawcprj_free(cwaveprj0)
 ABI_DT_FREE(cwaveprj0)

end subroutine eph_phpi
!!***

!----------------------------------------------------------------------

!!****f* m_phpi/out_phpi
!! NAME
!!  out_phpi
!!
!! FUNCTION
!!  Output the phonon self-energy.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_phpi
!!
!! NOTES
!!
!! CHILDREN
!!
!! SOURCE

subroutine out_phpi(iout, Pi_ph, phfrq, qpt, natom3)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout
 integer,intent(in) :: natom3
!arrays
 real(dp),intent(in) :: Pi_ph(natom3),phfrq(natom3),qpt(3)

!Local variables ------------------------------
!scalars
 integer :: imode

 write(iout,'(a)')' '
 !write(iout,'(a)')' ----------------------------------------'
 !write(iout,'(a)')' '
 write(iout,'(a)')' Phonon self-energy (Hartree)'
 write(iout,'(a)')' '
 write(iout,'(a,3f14.8)')' qpt =',qpt
 write(iout,'(a)')' '
 write(iout,'(1x,a,10x,a)')'omega','Pi(omega)'

 do imode=1,natom3
    write(iout,'(1x,f12.8,1x,es14.6)') phfrq(imode), Pi_ph(imode)
 end do

 write(iout,'(a)')' '
 !write(iout,'(a)')' ----------------------------------------'

end subroutine out_phpi
!!***

!----------------------------------------------------------------------

!!****f* m_phpi/out_phpi_nc
!! NAME
!!  out_phpi_nc
!!
!! FUNCTION
!!  Output the phonon self-energy in netCDF format.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_phpi
!!
!! NOTES
!!
!! CHILDREN
!!
!! SOURCE

subroutine out_phpi_nc(dtfil, cryst, Pi_ph, phfrq, qpt, natom3)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom3
 type(datafiles_type), intent(in) :: dtfil
 type(crystal_t),intent(in) :: cryst
!arrays
 real(dp),intent(in) :: Pi_ph(natom3),phfrq(natom3),qpt(3)

!Local variables ------------------------------
!scalars
 integer :: natom,one_dim,cplex,cart_dir
 integer :: ncid, ncerr
 character(len=fnlen) :: fname

#ifdef HAVE_NETCDF

 ! Initialize NetCDF file.
 fname = strcat(dtfil%filnam_ds(4),"_Pi.nc")
 NCF_CHECK(nctk_open_create(ncid, fname, xmpi_comm_self))

 ! Write information of the crystal
 NCF_CHECK(cryst%ncwrite(ncid))

 ! Write the dimensions specified by ETSF
 one_dim = 1
 cplex = 2
 cart_dir = 3

 natom = natom3 / 3

 ncerr = nctk_def_dims(ncid, [&
   nctkdim_t('current_one_dim', one_dim), &
   nctkdim_t('number_of_atoms', natom), &
   nctkdim_t('number_of_cartesian_directions', cart_dir), &
   nctkdim_t('number_of_perturbations', natom3), &
   nctkdim_t('cplex',cplex)], defmode=.True.)
 NCF_CHECK(ncerr)

 ! Create the arrays
 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t('q_point_reduced_coord', "dp", 'number_of_cartesian_directions'),&
   nctkarr_t('phonon_frequencies', "dp", 'number_of_perturbations'), &
   nctkarr_t('phonon_self_energy_realpart', "dp", 'number_of_perturbations')])
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_atomic_units(ncid, 'phonon_frequencies'))
 NCF_CHECK(nctk_set_atomic_units(ncid, 'phonon_self_energy_realpart'))

! Write data
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'q_point_reduced_coord'), qpt))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'phonon_frequencies'), phfrq))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'phonon_self_energy_realpart'), Pi_ph))

 ! Close file
 NCF_CHECK(nf90_close(ncid))

#else
 MSG_ERROR("NETCDF support required to write Pi.nc file.")
#endif

end subroutine out_phpi_nc
!!***

!----------------------------------------------------------------------

end module m_phpi
!!***
