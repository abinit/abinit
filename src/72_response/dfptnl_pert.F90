!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfptnl_pert
!! NAME
!! dfptnl_pert
!!
!! FUNCTION
!! Compute the linear response part to the 3dte
!!
!! COPYRIGHT
!! Copyright (C) 2016-2016 ABINIT group (LB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol) = array for planewave
!!                                          coefficients of wavefunctions
!!  cg1 = first-order wavefunction relative to the perturbations i1pert
!!  cg3 = first-order wavefunction relative to the perturbations i3pert
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!          if 2, COMPLEX
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  i1dir,i2dir,i3dir=directions of the corresponding perturbations
!!  i1pert,i2pert,i3pert = type of perturbation that has to be computed
!!  kg(3,mpw*mkmem)=reduced planewave coordinates
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  mband = maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem = maximum number of k points which can fit in core memory
!!  mk1mem = maximum number of k points for first-order WF
!!           which can fit in core memory
!!  mpert =maximum number of ipert
!!  mpi_enreg=MPI-parallelisation information
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  natom = number of atoms in unit cell
!!  nfft  = (effective) number of FFT grid points (for this processor)
!!  nkpt  = number of k points
!!  nspden = number of spin-density components
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  npwarr(nkpt) = array holding npw for each k point
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  rprimd(3,3) = dimensional primitive translations (bohr)
!!  vtrial1(cplex*nfft,nspden)=firs-order local potential
!!  xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!  d3etot(2,3,mpert,3,mpert,3,mpert) = matrix of the 3DTEs
!!
!! PARENTS
!!      dfptnl_loop
!!
!! CHILDREN
!!      destroy_hamiltonian,dotprod_g,fftpac,fourwf,init_hamiltonian
!!      load_k_hamiltonian,mkffnl,mkkpg,nonlop,status,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfptnl_pert(cg,cg1,cg3,cplex,dtfil,dtset,d3etot,eigen2,gs_hamkq,i1dir,i2dir,i3dir,&
& i1pert,i2pert,i3pert,kg,mband,mgfft,mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,natom,nfftf,nkpt,&
& nspden,nspinor,nsppol,npwarr,occ,pawfgr,ph1d,psps,rf_hamkq,rprimd,vtrial,vtrial1,wffddk,ddk_f,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_wffile
 use m_wfk
 use m_xmpi

 use m_cgtools,    only : dotprod_g
 use m_hamiltonian,only : gs_hamiltonian_type,load_k_hamiltonian,load_spin_hamiltonian,&
                          load_spin_rf_hamiltonian,load_spin_rf_hamiltonian,rf_hamiltonian_type
 use m_pawtab,     only : pawtab_type
 use m_pawcprj,    only : pawcprj_type
 use m_pawfgr,     only : pawfgr_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptnl_pert'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_56_recipspace
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,mband,mgfft
 integer,intent(in) :: mk1mem,mkmem,mpert,mpsang,mpw,natom,nfftf,nkpt,nspden
 integer,intent(in) :: nspinor,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(pawfgr_type),intent(in) :: pawfgr
 type(rf_hamiltonian_type),intent(inout) :: rf_hamkq
 type(wffile_type),intent(inout) :: wffddk(2)
 type(wfk_t),intent(inout) :: ddk_f(2)

!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: cg3(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: eigen2(2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom),rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom),vtrial(cplex*nfftf,nspden)
 real(dp),intent(inout) :: vtrial1(cplex*nfftf,nspden),d3etot(2,3,mpert,3,mpert,3,mpert)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=52
 integer :: bantot,choice,counter,cpopt,dimffnl,iband,icg0,ider,ierr,iexit
 integer :: ii,igs,ikg,ikg1,ikpt,ilm,ipw,isppol,ispinor,istwf_k,jband,jj
 integer :: me,n1,n2,n3,n4,n5,n6,nband_k,nkpg,nkpg1,nnlout,npw_k,npw1_k
 integer :: offset_cg,offset_eigen,option,paw_opt,esigns,size_wf,spaceComm,tim_fourwf,tim_nonlop,useylmgr1
 real(dp) :: dot1i,dot1r,dot2i,dot2r,doti,dotr,lagi,lagr,sumi,sumr,tol_test,weight
!arrays
 integer,allocatable :: kg_k(:,:),kg1_k(:,:)
 real(dp) :: buffer(2),enlout(3),kpq(3),kpt(3)
 real(dp) :: dum_svectout(1,1),dum(1),rmet(3,3),dummy_ylmgr(1,1,1)
 real(dp),allocatable :: cwave0(:,:),cwavef3(:,:),dkinpw(:),eig1_k_tmp(:)
 real(dp),allocatable :: ffnl1(:,:,:,:),ffnlk(:,:,:,:)
 real(dp),allocatable :: dudk(:,:),dudkde(:,:),kinpw1(:),dummy_array(:)
 real(dp),allocatable :: gh0(:,:),gh1(:,:),gvnl(:,:),kpg_k(:,:),kpg1_k(:,:)
 real(dp),allocatable :: ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: vlocal(:,:,:,:),vlocal1(:,:,:,:),wfraug(:,:,:,:),work1(:,:)
 real(dp),allocatable :: ylm(:,:),ylm1(:,:),ylmgr1(:,:,:)
 real(dp),allocatable :: ylm_k(:,:),ylm1_k(:,:),ylmgr1_k(:,:,:)
 type(pawcprj_type) :: cprj_dum(1,1)
 type(pawtab_type) :: pawtab_dum(0)
 type(rf_hamiltonian_type) :: rf_ham_dum
!***********************************************************************

 me = mpi_enreg%me
 spaceComm=mpi_enreg%comm_cell

 call status(0,dtfil%filstat,iexit,level,'enter         ')

 bantot = 0
 icg0 = 0
 option = 2
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)

 ABI_ALLOCATE(vlocal,(n4,n5,n6,gs_hamkq%nvloc))
 ABI_ALLOCATE(vlocal1,(cplex*n4,n5,n6,gs_hamkq%nvloc))

 ABI_ALLOCATE(wfraug,(2,n4,n5,n6))

!!Initialize Hamiltonian (k-independent terms) - NCPP only
! call init_hamiltonian(gs_hamk,psps,pawtab_dum,nspinor,nspden,natom,&
!& dtset%typat,xred,nfft,mgfft,dtset%ngfft,rprimd,dtset%nloalg,ph1d=ph1d)
 rmet = MATMUL(TRANSPOSE(rprimd),rprimd)

 sumr = zero ; sumi = zero

!Set up the Ylm for each k point
 ABI_ALLOCATE(ylm,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))
 if (psps%useylm==1) then
   call status(0,dtfil%filstat,iexit,level,'call initylmg ')
   option=2
   call initylmg(gs_hamkq%gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%nband,&
   dtset%nkpt,npwarr,dtset%nsppol,option,rprimd,ylm,dummy_ylmgr)
 end if
 
!Set up the spherical harmonics (Ylm) at k+q
 useylmgr1=0; option=0
 if (psps%useylm==1.and. &
& (i2pert==natom+1.or.i2pert==natom+3.or.i2pert==natom+4.or.(psps%usepaw==1.and.i2pert==natom+2))) then
   useylmgr1=1; option=1
 end if
 ABI_ALLOCATE(ylm1,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))
 ABI_ALLOCATE(ylmgr1,(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))
!To change the following when q/=0
 if (psps%useylm==1) then
   call initylmg(gs_hamkq%gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%nband,&
   dtset%nkpt,npwarr,dtset%nsppol,option,rprimd,ylm1,ylmgr1)
 end if

 size_wf = dtset%nspinor*npw_k

!Loop over spins

 do isppol = 1, nsppol

!   call status(0,dtfil%filstat,iexit,level,'call fftpac   ')
!   call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,vtrial1,vlocal1,option)

!  Set up local potential vlocal1 with proper dimensioning, from vtrial1
!  Same thing for vlocal from vtrial Also take into account the spin.
   call rf_transgrid_and_pack(isppol,nspden,psps%usepaw,cplex,nfftf,dtset%nfft,dtset%ngfft,&
&   gs_hamkq%nvloc,pawfgr,mpi_enreg,vtrial,vtrial1,vlocal,vlocal1)

!  Continue to initialize the Hamiltonian
   call load_spin_hamiltonian(gs_hamkq,isppol,vlocal=vlocal, &
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
!   call load_spin_hamiltonian(gs_hamkq,isppol,paw_ij=paw_ij,vlocal=vlocal, &
!&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

   call load_spin_rf_hamiltonian(rf_hamkq,gs_hamkq,isppol,vlocal1=vlocal1, &
   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
!   call load_spin_rf_hamiltonian(rf_hamkq,gs_hamkq,isppol,paw_ij1=paw_ij1,vlocal1=vlocal1, &
!   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

!  Loop over k-points

   ikg = 0
   do ikpt = 1, nkpt

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,-1,mpi_enreg%me))cycle

     counter = 100*ikpt

     nband_k = dtset%nband(ikpt+(isppol-1)*nkpt)
     npw_k = npwarr(ikpt)
     npw1_k = npw_k ! To change for q/=0
     istwf_k = dtset%istwfk(ikpt)

     kpt(:) = dtset%kptns(:,ikpt)
     kpq(:) = dtset%kptns(:,ikpt) ! In case of non zero q, kpt = kpt + q

     ABI_ALLOCATE(cwave0,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(cwavef3,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(gh0,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(gvnl,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(gh1,(2,npw_k*dtset%nspinor))

     ABI_ALLOCATE(kg_k,(3,npw_k))
     ABI_ALLOCATE(kg1_k,(3,npw1_k))
     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     ABI_ALLOCATE(ylm1_k,(npw1_k,mpsang*mpsang*psps%useylm))
     ABI_ALLOCATE(ylmgr1_k,(npw1_k,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))

!    Get (k+G) wave vectors and associated spherical harmonics
     kg_k(:,1:npw_k) = kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,mpsang*mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if

!    Get (k+q+G) wave vectors and associated spherical harmonics
     kg1_k(:,1:npw1_k)=kg(:,1+ikg1:npw1_k+ikg1) ! To change for q/=0
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm1_k(1:npw1_k,ilm)=ylm1(1+ikg1:npw1_k+ikg1,ilm)
       end do
       if (useylmgr1==1) then
         do ilm=1,psps%mpsang*psps%mpsang
           do ii=1,3
             ylmgr1_k(1:npw1_k,ii,ilm)=ylmgr1(1+ikg1:npw1_k+ikg1,ii,ilm)
           end do
         end do
       end if
     end if

!    Set up the ground-state Hamiltonian, and some parts of the 1st-order Hamiltonian
     call getgh1c_setup(gs_hamkq,rf_hamkq,dtset,psps,&                              ! In
     kpt,kpq,i2dir,i2pert,natom,rmet,gs_hamkq%gprimd,gs_hamkq%gmet,istwf_k,&        ! In
     npw_k,npw1_k,useylmgr1,kg_k,ylm_k,kg1_k,ylm1_k,ylmgr1_k,&                      ! In
     dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1,&                 ! Out
     dummy_array,dummy_array,rf_ham_dum)                                            ! Out

     ABI_STAT_ALLOCATE(dudk,  (2,nband_k*size_wf), ierr)
     ABI_STAT_ALLOCATE(dudkde,(2,nband_k*size_wf), ierr)
     ABI_STAT_ALLOCATE(eig1_k_tmp,(2*nband_k), ierr)
     ABI_STAT_ALLOCATE(work1,(2,size_wf), ierr)
     do iband = 1,nband_k
!      TO READ dudkde
#ifndef DEV_MG_WFK
       call WffReadDataRec(eig1_k_tmp,ierr,2*nband_k,wffddk(file_index(kdir1)))
       call WffReadDataRec(work1,ierr,2,size_wf,wffddk(file_index(kdir1)))
#else
       call wfk_read_bks(ddk_f(1), iband, ikpt, isppol, xmpio_single, cg_bks=work1,eig1_bks=eig1_k_tmp)
#endif
!      Filter the wavefunctions for large modified kinetic energy
!      The GS wavefunctions should already be non-zero
       do ispinor=1,gs_hamkq%nspinor
         igs=(ispinor-1)*npw_k
         do ipw=1+igs,npw_k+igs
           if(gs_hamkq%kinpw_kp(ipw-igs)>huge(zero)*1.d-11)then
             work1(1,ipw)=zero
             work1(2,ipw)=zero
           end if
         end do
       end do
!      Copy work1 in "dudk"
       dudk(:,1+(iband-1)*size_wf:size_wf+(iband-1)*size_wf)=work1(:,:)
!!      Copy eig1_k_tmp in "eig1_k_stored"
!       shift_band1=(iband-1)*2*nband_k
!       shift_dir1_lambda=2*(kdir1-1)*nband_k**2
!       eig1_k_stored(1+shift_band1+shift_dir1_lambda:2*nband_k+shift_band1+shift_dir1_lambda)=eig1_k_tmp(:)
     end do
     ABI_DEALLOCATE(eig1_k_tmp)

!    Loop over bands
     do jband = 1,nband_k

       offset_cg    = size_wf      *(jband-1+nband_k*    (ikpt-1+nkpt*(isppol-1)))
       offset_eigen = 2*dtset%mband*(jband-1+dtset%mband*(ikpt-1+nkpt*(isppol-1)))

!!      Test if < u^(0) | ( H^(1) - eps^(0) S^(1) ) | u^(0) > = eig^(1)
!       tol_test = tol8
!       ABI_ALLOCATE(cgj,(2,size_wf))
!       ABI_ALLOCATE(work2,(2,size_wf))
!       ABI_ALLOCATE(iddk,(2,size_wf))
!       cgj(:,:) = cg(:,1+offset_cg:size_wf+offset_cg)
!       iddk(:,:) = zero
!       sij_opt = 0
!       optlocal = 1
!       opt_gvnl1 = 1
!       usevnl = 1
!       call getgh1c(berryopt,0,cgj,cwaveprj,work1,dum1,work2,gs_hamkq,iddk,i2dir,i2pert,zero,&
!                    mpi_enreg,optlocal,optnl,opt_gvnl1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)
!       do iband=1,nband_k
!         offset_cg    = size_wf      *(iband-1+nband_k*    (ikpt-1+nkpt*(isppol-1)))
!         cgj(:,:) = cg(:,1+offset_cg:size_wf+offset_cg)
!         call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cgj,work1,mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
!         dotr = dotr - 0 ! TO CHANGE
!         doti = doti - 0 ! TO CHANGE
!         dotr = sqrt(dotr**2+doti**2)
!        if (dotr > tol_test) then
!           write(message,'(4(a,i2),a,es22.13E3)') 'RF2 TEST GETGH1 : ipert=',i2pert-natom,' idir=',i2dir,&
!                                              ' jband=',jband,' iband=',iband,' NOT PASSED dotr = ',dotr
!           call wrtout(std_out,message)
!        end if
!       end do ! end iband
!       ABI_DEALLOCATE(cgj)
!       ABI_DEALLOCATE(work1)
!       ABI_DEALLOCATE(work2)
!       ABI_DEALLOCATE(iddk)

!       cwave0(:,:)=cg(:,1+(iband - 1)*npw_k*dtset%nspinor+icg0:&
!&       iband*npw_k*dtset%nspinor+icg0)
!       cwavef3(:,:)=cg3(:,1+(iband-1)*npw_k*dtset%nspinor+icg0:&
!&       iband*npw_k*dtset%nspinor+icg0)

!!      Compute vtrial1 | cwafef3 >
!       tim_fourwf = 0 ; weight = one
!       call status(counter,dtfil%filstat,iexit,level,'call fourwf  ')
!       call fourwf(cplex,vlocal1,cwavef3,gh1,wfraug,gs_hamkq%gbound_k,gs_hamkq%gbound_k,&
!&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,dtset%ngfft,npw_k,npw_k,n4,n5,n6,option,&
!&       dtset%paral_kgb,tim_fourwf,weight,weight,&
!&       use_gpu_cuda=dtset%use_gpu_cuda)

!!      In case i2pert = phonon-type perturbation
!!      add first-order change in the nonlocal potential
!       if (i2pert<natom+1) then
!         signs=2 ; choice=2 ; nnlout=3 ; tim_nonlop = 0 ; paw_opt=0 ; cpopt=-1
!         call status(counter,dtfil%filstat,iexit,level,'call nonlop  ')
!         call nonlop(choice,cpopt,cprj_dum,enlout,gs_hamkq,i2dir,dum,mpi_enreg,1,nnlout,paw_opt,&
!&         signs,dum_svectout,tim_nonlop,cwavef3,gvnl,iatom_only=i2pert)
!         gh1(:,:) = gh1(:,:) + gvnl(:,:)
!       end if

!       ii = (iband-1)*npw_k*dtset%nspinor + icg0
!       call dotprod_g(dotr,doti,istwf_k,npw_k,2,cg1(:,ii+1:ii+npw_k),gh1,mpi_enreg%me_g0,xmpi_comm_self)

!!      Compute vtrial1 | cwave0 >
!       tim_fourwf = 0 ; weight = one
!       call status(counter,dtfil%filstat,iexit,level,'call fourwf  ')
!       call fourwf(cplex,vlocal1,cwave0,gh0,wfraug,gs_hamkq%gbound_k,gs_hamkq%gbound_k,&
!&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,dtset%ngfft,npw_k,npw_k,n4,n5,n6,option,&
!&       dtset%paral_kgb,tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)

!!      In case i2pert = phonon-type perturbation
!!      add first-order change in the nonlocal potential
!       if (i2pert<natom+1) then
!         signs=2 ; choice=2 ; nnlout=3 ; tim_nonlop = 0 ; paw_opt=0 ; cpopt=-1
!         call status(counter,dtfil%filstat,iexit,level,'call nonlop  ')
!         call nonlop(choice,cpopt,cprj_dum,enlout,gs_hamkq,i2dir,dum,mpi_enreg,1,nnlout,paw_opt,&
!&         signs,dum_svectout,tim_nonlop,cwave0,gvnl,iatom_only=i2pert)
!         gh0(:,:) = gh0(:,:) + gvnl(:,:)
!       end if

!!      Compute the dft contribution to the Lagrange multiplier
!!      cwavef3 and cwave0 have been transferred to gh1 and gh0
!!      these vectors will be used to store the wavefunctions of band iband
!!      cg1 and gh0 contain the wavefunctions of band jband

!       lagr = zero ; lagi = zero
!       do jband = 1, nband_k

!         ii = (jband - 1)*npw_k*dtset%nspinor + icg0
!         jj = (iband - 1)*npw_k*dtset%nspinor + icg0

!!        dot1r and dot1i contain < u_mk | v^(1) | u_nk >
!!        dot2r and dot2i contain < u_nk^(1) | u_mk^(1) >
!!        m -> jband and n -> iband

!         dot1r = zero ; dot1i = zero
!         dot2r = zero ; dot2i = zero
!         do ipw = 1, npw_k
!           ii = ii + 1 ; jj = jj + 1
!           dot1r = dot1r + cg(1,ii)*gh0(1,ipw) + cg(2,ii)*gh0(2,ipw)
!           dot1i = dot1i + cg(1,ii)*gh0(2,ipw) - cg(2,ii)*gh0(1,ipw)
!           dot2r = dot2r + cg1(1,jj)*cg3(1,ii) + &
!&           cg1(2,jj)*cg3(2,ii)
!           dot2i = dot2i + cg1(1,jj)*cg3(2,ii) - &
!&           cg1(2,jj)*cg3(1,ii)
!         end do  !  ipw

!         lagr = lagr + dot1r*dot2r - dot1i*dot2i
!         lagi = lagi + dot1r*dot2i + dot1i*dot2r

!       end do    ! jband

       sumr = sumr + dtset%wtk(ikpt)*occ(bantot+iband)*(dotr-lagr)
       sumi = sumi + dtset%wtk(ikpt)*occ(bantot+iband)*(doti-lagi)

     end do   ! end loop over bands

     bantot = bantot + nband_k
     icg0 = icg0 + npw_k*dtset%nspinor*nband_k
     ikg = ikg + npw_k

     ABI_DEALLOCATE(cwave0)
     ABI_DEALLOCATE(cwavef3)
     ABI_DEALLOCATE(gh0)
     ABI_DEALLOCATE(gh1)
     ABI_DEALLOCATE(gvnl)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(kg1_k)
     ABI_DEALLOCATE(dudk)
     ABI_DEALLOCATE(dudkde)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylm1_k)
     ABI_DEALLOCATE(ylmgr1_k)
     ABI_DEALLOCATE(ffnlk)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(kpg1_k)

   end do   ! end loop over k-points

 end do   ! end loop over spins

 if (xmpi_paral == 1) then
   buffer(1) = sumr ; buffer(2) = sumi
   call xmpi_sum(buffer,spaceComm,ierr)
   sumr = buffer(1) ; sumi = buffer(2)
 end if


 d3etot(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sumr
!d3etot(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sumi

!In some cases, the imaginary part is /= 0 because of the
!use of time reversal symmetry
 d3etot(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = zero

 ABI_DEALLOCATE(ylm)
 ABI_DEALLOCATE(ylm1)
 ABI_DEALLOCATE(ylmgr1)
 ABI_DEALLOCATE(vlocal1)
 ABI_DEALLOCATE(wfraug)

 call status(0,dtfil%filstat,iexit,level,'exit          ')

end subroutine dfptnl_pert
!!***
