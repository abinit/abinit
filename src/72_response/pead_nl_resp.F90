!{\src2tex{textfont=tt}}
!!****f* ABINIT/pead_nl_resp
!! NAME
!! pead_nl_resp
!!
!! FUNCTION
!! Compute the linear response part to the 3dte
!!
!! COPYRIGHT
!! Copyright (C) 2002-2017 ABINIT group (MVeithen)
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
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= spherical harmonics for
!!       each G and k point
!!
!! OUTPUT
!!  d3lo(2,3,mpert,3,mpert,3,mpert) = matrix of the 3DTEs
!!
!! PARENTS
!!      pead_nl_loop
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


subroutine pead_nl_resp(cg,cg1,cg3,cplex,dtfil,dtset,d3lo,&
& i1dir,i2dir,i3dir,i1pert,i2pert,i3pert,&
& kg,mband,mgfft,mkmem,mk1mem,&
& mpert,mpi_enreg,mpsang,mpw,natom,nfft,nkpt,nspden,nspinor,nsppol,&
& npwarr,occ,ph1d,psps,rprimd,vtrial1,xred,ylm)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi

 use m_cgtools,    only : dotprod_g
 use m_pawtab,     only : pawtab_type
 use m_pawcprj,    only : pawcprj_type
 use m_hamiltonian,only : init_hamiltonian,destroy_hamiltonian,&
&                         load_k_hamiltonian,gs_hamiltonian_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pead_nl_resp'
 use interfaces_32_util
 use interfaces_53_ffts
 use interfaces_66_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,mband,mgfft
 integer,intent(in) :: mk1mem,mkmem,mpert,mpsang,mpw,natom,nfft,nkpt,nspden
 integer,intent(in) :: nspinor,nsppol
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: cg3(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom),rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom),ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(inout) :: vtrial1(cplex*nfft,nspden)
 real(dp),intent(inout) :: d3lo(2,3,mpert,3,mpert,3,mpert)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=52
 integer :: bantot,choice,counter,cpopt,dimffnl,iband,icg0,ider,ierr,iexit
 integer :: ii,ikg,ikpt,ilm,ipw,isppol,istwf_k,jband,jj
 integer :: me,n1,n2,n3,n4,n5,n6,nband_k,nkpg,nnlout,npw_k
 integer :: option,paw_opt,signs,spaceComm,tim_fourwf,tim_nonlop
 real(dp) :: dot1i,dot1r,dot2i,dot2r,doti,dotr,lagi,lagr,sumi,sumr,weight
 type(gs_hamiltonian_type) :: gs_hamk
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: buffer(2),enlout(3),kpq(3),kpt(3)
 real(dp) :: dum_svectout(1,1),dum(1),rmet(3,3),ylmgr_dum(1,1,1)
 real(dp),allocatable :: cwave0(:,:),cwavef3(:,:),ffnlk(:,:,:,:)
 real(dp),allocatable :: gh0(:,:),gh1(:,:),gvnl(:,:),kpg_k(:,:)
 real(dp),allocatable :: vlocal1(:,:,:),wfraug(:,:,:,:),ylm_k(:,:)
 type(pawcprj_type) :: cprj_dum(1,1)
 type(pawtab_type) :: pawtab_dum(0)

!***********************************************************************

 me = mpi_enreg%me
 spaceComm=mpi_enreg%comm_cell

 call status(0,dtfil%filstat,iexit,level,'enter         ')

 bantot = 0
 icg0 = 0
 option = 2
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)

 ABI_ALLOCATE(vlocal1,(cplex*n4,n5,n6))
 ABI_ALLOCATE(wfraug,(2,n4,n5,n6))

!Initialize Hamiltonian (k-independent terms) - NCPP only
 call init_hamiltonian(gs_hamk,psps,pawtab_dum,nspinor,nsppol,nspden,natom,&
& dtset%typat,xred,nfft,mgfft,dtset%ngfft,rprimd,dtset%nloalg,ph1d=ph1d,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab)
!& paw_ij=paw_ij)
 rmet = MATMUL(TRANSPOSE(rprimd),rprimd)

 sumr = zero ; sumi = zero

!Loop over spins

 do isppol = 1, nsppol

   call status(0,dtfil%filstat,iexit,level,'call fftpac   ')
   call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,vtrial1,vlocal1,option)

!  Loop over k-points

   ikg = 0
   do ikpt = 1, nkpt

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,-1,mpi_enreg%me))cycle

     counter = 100*ikpt

     nband_k = dtset%nband(ikpt+(isppol-1)*nkpt)
     npw_k = npwarr(ikpt)
     istwf_k = dtset%istwfk(ikpt)

     kpt(:) = dtset%kptns(:,ikpt)
     kpq(:) = dtset%kptns(:,ikpt) ! In case of non zero q, kpt = kpt + q

     ABI_ALLOCATE(cwave0,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(cwavef3,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(gh0,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(gvnl,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(gh1,(2,npw_k*dtset%nspinor))

     ABI_ALLOCATE(kg_k,(3,npw_k))
     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     kg_k(:,1:npw_k) = kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,mpsang*mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if

!    Compute (k+G) and (k+q+G) vectors (only if useylm=1)
     nkpg=0;if (i2pert<natom+1) nkpg=3*dtset%nloalg(3)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpt,nkpg,npw_k)
     end if

!    Compute nonlocal form factors ffnl at (k+G), for all atoms
     dimffnl=1
     ABI_ALLOCATE(ffnlk,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
     if (i2pert<natom+1) then
       ider=0
       call status(counter,dtfil%filstat,iexit,level,'call mkffnl  ')
       call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnlk,psps%ffspl,gs_hamk%gmet,gs_hamk%gprimd,&
&       ider,ider,psps%indlmn,kg_k,kpg_k,kpt,psps%lmnmax,psps%lnmax,psps%mpsang,&
&       psps%mqgrid_ff,nkpg,npw_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&       psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
     end if

!    Load k-dependent part in the Hamiltonian datastructure
     call load_k_hamiltonian(gs_hamk,kpt_k=kpt,npw_k=npw_k,istwf_k=istwf_k,&
&     kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnlk,compute_gbound=.true.)
!    Load k+q-dependent part in the Hamiltonian datastructure
!    call load_kprime_hamiltonian...  !! To be activated when q/=0

!    Loop over bands

     do iband = 1,nband_k

       cwave0(:,:)=cg(:,1+(iband - 1)*npw_k*dtset%nspinor+icg0:&
&       iband*npw_k*dtset%nspinor+icg0)
       cwavef3(:,:)=cg3(:,1+(iband-1)*npw_k*dtset%nspinor+icg0:&
&       iband*npw_k*dtset%nspinor+icg0)

!      Compute vtrial1 | cwafef3 >
       tim_fourwf = 0 ; weight = one
       call status(counter,dtfil%filstat,iexit,level,'call fourwf  ')
       call fourwf(cplex,vlocal1,cwavef3,gh1,wfraug,gs_hamk%gbound_k,gs_hamk%gbound_k,&
&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,dtset%ngfft,npw_k,npw_k,n4,n5,n6,option,&
&       dtset%paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=dtset%use_gpu_cuda)

!      In case i2pert = phonon-type perturbation
!      add first-order change in the nonlocal potential
       if (i2pert<natom+1) then
         signs=2 ; choice=2 ; nnlout=3 ; tim_nonlop = 0 ; paw_opt=0 ; cpopt=-1
         call status(counter,dtfil%filstat,iexit,level,'call nonlop  ')
         call nonlop(choice,cpopt,cprj_dum,enlout,gs_hamk,i2dir,dum,mpi_enreg,1,nnlout,paw_opt,&
&         signs,dum_svectout,tim_nonlop,cwavef3,gvnl,iatom_only=i2pert)
         gh1(:,:) = gh1(:,:) + gvnl(:,:)
       end if

       ii = (iband-1)*npw_k*dtset%nspinor + icg0
       call dotprod_g(dotr,doti,istwf_k,npw_k,2,cg1(:,ii+1:ii+npw_k),gh1,mpi_enreg%me_g0,xmpi_comm_self)

!      Compute vtrial1 | cwave0 >
       tim_fourwf = 0 ; weight = one
       call status(counter,dtfil%filstat,iexit,level,'call fourwf  ')
       call fourwf(cplex,vlocal1,cwave0,gh0,wfraug,gs_hamk%gbound_k,gs_hamk%gbound_k,&
&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,dtset%ngfft,npw_k,npw_k,n4,n5,n6,option,&
&       dtset%paral_kgb,tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)

!      In case i2pert = phonon-type perturbation
!      add first-order change in the nonlocal potential
       if (i2pert<natom+1) then
         signs=2 ; choice=2 ; nnlout=3 ; tim_nonlop = 0 ; paw_opt=0 ; cpopt=-1
         call status(counter,dtfil%filstat,iexit,level,'call nonlop  ')
         call nonlop(choice,cpopt,cprj_dum,enlout,gs_hamk,i2dir,dum,mpi_enreg,1,nnlout,paw_opt,&
&         signs,dum_svectout,tim_nonlop,cwave0,gvnl,iatom_only=i2pert)
         gh0(:,:) = gh0(:,:) + gvnl(:,:)
       end if

!      Compute the dft contribution to the Lagrange multiplier
!      cwavef3 and cwave0 have been transferred to gh1 and gh0
!      these vectors will be used to store the wavefunctions of band iband
!      cg1 and gh0 contain the wavefunctions of band jband

       lagr = zero ; lagi = zero
       do jband = 1, nband_k

         ii = (jband - 1)*npw_k*dtset%nspinor + icg0
         jj = (iband - 1)*npw_k*dtset%nspinor + icg0

!        dot1r and dot1i contain < u_mk | v^(1) | u_nk >
!        dot2r and dot2i contain < u_nk^(1) | u_mk^(1) >
!        m -> jband and n -> iband

         dot1r = zero ; dot1i = zero
         dot2r = zero ; dot2i = zero
         do ipw = 1, npw_k
           ii = ii + 1 ; jj = jj + 1
           dot1r = dot1r + cg(1,ii)*gh0(1,ipw) + cg(2,ii)*gh0(2,ipw)
           dot1i = dot1i + cg(1,ii)*gh0(2,ipw) - cg(2,ii)*gh0(1,ipw)
           dot2r = dot2r + cg1(1,jj)*cg3(1,ii) + &
&           cg1(2,jj)*cg3(2,ii)
           dot2i = dot2i + cg1(1,jj)*cg3(2,ii) - &
&           cg1(2,jj)*cg3(1,ii)
         end do  !  ipw

         lagr = lagr + dot1r*dot2r - dot1i*dot2i
         lagi = lagi + dot1r*dot2i + dot1i*dot2r

       end do    ! jband

       sumr = sumr + &
&       dtset%wtk(ikpt)*occ(bantot+iband)*(dotr-lagr)
       sumi = sumi + &
&       dtset%wtk(ikpt)*occ(bantot+iband)*(doti-lagi)

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
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ffnlk)
     ABI_DEALLOCATE(kpg_k)

   end do   ! end loop over k-points

 end do   ! end loop over spins

 if (xmpi_paral == 1) then
   buffer(1) = sumr ; buffer(2) = sumi
   call xmpi_sum(buffer,spaceComm,ierr)
   sumr = buffer(1) ; sumi = buffer(2)
 end if


 d3lo(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sumr
!d3lo(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sumi

!In some cases, the imaginary part is /= 0 because of the
!use of time reversal symmetry
 d3lo(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = zero

 call destroy_hamiltonian(gs_hamk)

 ABI_DEALLOCATE(vlocal1)
 ABI_DEALLOCATE(wfraug)

 call status(0,dtfil%filstat,iexit,level,'exit          ')

end subroutine pead_nl_resp
!!***
