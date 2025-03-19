!!****m* ABINIT/m_datafordmft
!! NAME
!!  m_datafordmft
!!
!! FUNCTION
!! This module produces inputs for the DMFT calculation
!!
!! COPYRIGHT
!! Copyright (C) 2006-2025 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_datafordmft

 use defs_abitypes
 use defs_basis
 use m_abi_linalg, only : abi_xgemm
 use m_crystal, only : crystal_t
 use m_dtset
 use m_errors
 use m_fstrings, only : int2char4
 use m_io_tools, only : open_file
 use m_matlu, only : add_matlu,checkdiag_matlu,copy_matlu,destroy_matlu,diff_matlu, &
                   & init_matlu,matlu_type,print_matlu,sym_matlu,xmpi_matlu
 use m_matrix, only : invsqrt_matrix
 use m_mpinfo, only : proc_distrb_cycle
 use m_oper, only : copy_oper,destroy_oper,diff_oper,downfold_oper,identity_oper,init_oper,oper_type,prod_oper
 use m_paw_dmft, only : paw_dmft_type
 use m_paw_ij, only : paw_ij_type
 use m_pawcprj, only : pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_type
 use m_pawrad, only : pawrad_type,simp_gen
 use m_pawtab, only : pawtab_type
 use m_xmpi

 implicit none

 private

 public :: datafordmft
 public :: chipsi_print
 public :: compute_levels
 public :: chipsi_renormalization
 public :: hybridization_asymptotic_coefficient
 public :: compute_wannier
 public :: print_wannier
!!***

contains

!!****f* ABINIT/datafordmft
!! NAME
!! datafordmft
!!
!! FUNCTION
!!  Compute chipsi (and print some data for check)
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!        -gprimd(3,3)=dimensional reciprocal space primitive translations
!!        -indsym(4,nsym,natom)=indirect indexing array for atom labels
!!        -symrec(3,3,nsym)=symmetry operations in reciprocal space
!!        - nsym= number of symetry operations
!!  dft_occup <type(oper_type)> = DFT occupations of the correlated orbitals
!!  dimcprj(natom) = dimension for cprj
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding DFT eigenvalues (hartree)
!!  mband_cprj=number of bands on each process of the band communicator
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=information about MPI parallelization
!!  my_nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  occ(mband*nkpt*nsppol) = occupancies of KS states.
!!  paw_dmft <type(paw_dmft_type)>= paw+dmft related data
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  usecprj=1 if cprj datastructure is stored in memory
!!  nbandkss=number of bands in the KSS file
!!
!! OUTPUT
!!  paw_dmft%chipsi((2*maxlpawu+1)*nspinor,mbandc,nkpt,nsppol,natom): projections <Chi|Psi>
!!  paw_dmft%eigen(paw_dmft%mbandc,paw_dmft%nkpt,paw_dmft%nsppol)
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! SOURCE

subroutine datafordmft(cg,cprj,cryst_struc,dft_occup,dimcprj,dtset,eigen,mband_cprj,mcg,&
                     & mpi_enreg,my_nspinor,occ,paw_dmft,paw_ij,pawtab,usecprj,nbandkss)

!Arguments ------------------------------------
 integer, intent(in) :: mband_cprj,mcg,my_nspinor,usecprj
 integer, optional, intent(in) :: nbandkss
 type(MPI_type), intent(in) :: mpi_enreg
 type(dataset_type), intent(in) :: dtset
 type(oper_type), intent(inout) :: dft_occup !vz_i
 type(crystal_t), intent(in) :: cryst_struc
 type(paw_dmft_type), intent(inout) :: paw_dmft
 integer, intent(in) :: dimcprj(paw_dmft%natom)
 real(dp), intent(in) :: occ(paw_dmft%mband*paw_dmft%nkpt*paw_dmft%nsppol)
 real(dp), target, intent(in) :: eigen(paw_dmft%mband*paw_dmft%nkpt*paw_dmft%nsppol)
 real(dp), intent(in) :: cg(2,mcg)
 type(paw_ij_type), intent(in) :: paw_ij(paw_dmft%natom)
! type(pawcprj_type) :: cprj(cryst_struc%natom,my_nspinor*mband*mkmem*nsppol)
 type(pawcprj_type), intent(in) :: cprj(paw_dmft%natom,my_nspinor*mband_cprj*dtset%mkmem*paw_dmft%nsppol*usecprj)
 type(pawtab_type), intent(in) :: pawtab(paw_dmft%ntypat)
!Local variables-------------------------------
 integer :: band_index,comm_band,comm_kpt,iatom,ib,iband,ibandc,ibuf_chipsi,ibuf_psi,icg,icgb
 integer :: icprj,idijeff,ierr,ik,ikpt,ilmn,im,iorder_cprj,iproj,ir,irank,ispinor,ispinor1
 integer :: isppol,itypat,lmn_size,lpawu,lpawu1,maxlpawu,maxmeshsize,maxnproju,mband,mbandc
 integer :: me_band,me_kpt,mkmem,natom,nband_k,nband_k_cprj,nbandf,nbandi,ndim,nkpt
 integer :: nproc_band,nproc_spkpt,nproju,npw,nspinor,nsploop,nsppol,nsppol_mem,opt_renorm
 integer :: option,paral_kgb,pawprtvol,siz_buf,siz_buf_psi,siz_paw,siz_proj,siz_wan,unt
 logical :: prt_wan,t2g,use_full_chipsi,verif_band,x2my2d
 real(dp) :: rint
 character(len=500) :: message
 type(oper_type) :: loc_norm_check
 integer, allocatable :: displs(:),recvcounts(:)
 complex(dpc), allocatable :: buf_chipsi(:),buf_chipsi_tot(:),chipsi_tmp(:),cwprj(:,:)
 type(pawcprj_type), allocatable :: cwaveprj(:,:)
 type(matlu_type), allocatable :: matlu_temp(:)
 integer, parameter :: spinor_idxs(2,4) = RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
!************************************************************************

!DBG_ENTER("COLL")

 mband   = paw_dmft%mband
 mbandc  = paw_dmft%mbandc
 mkmem   = paw_dmft%mkmem
 natom   = paw_dmft%natom
 nkpt    = paw_dmft%nkpt
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 nsppol_mem = sum(mpi_enreg%my_isppoltab(1:nsppol))
 pawprtvol  = dtset%pawprtvol
 prt_wan    = (paw_dmft%dmft_prtwan == 1)

 if (abs(pawprtvol) >= 3) then
   write(message,*) " number of k-points used is nkpt=nkpt ",nkpt
   call wrtout(std_out,message,'COLL')
   write(message,*) " warning: parallelized version        ",nkpt
   call wrtout(std_out,message,'COLL')
   write(message,*) " weights k-points used is wtk=wtk"
   call wrtout(std_out,message,'COLL')
 end if ! abs(pawprtvol)>=3

 if (usecprj == 0) then
   write(message,*) "  usecprj=0 : BUG in datafordmft",usecprj
   ABI_BUG(message)
 end if

 if (my_nspinor /= nspinor) then
   write(message,*) "  my_nspinor=/dtset%nspinor, datafordmft not working in this case",my_nspinor,nspinor
   ABI_ERROR(message)
 end if


!do ib=1,my_nspinor*mband_cprj*mkmem*nsppol*usecprj
!write(std_out,'(a,i6,3e16.7)') "cprj",ib,cprj(1,ib)%cp(1,19),cprj(1,ib)%cp(2,19),cprj(1,ib)%cp(1,19)**2+cprj(1,ib)%cp(2,19)**2
!enddo

!----------------------------------- MPI-------------------------------------

! Init parallelism
 paral_kgb = mpi_enreg%paral_kgb
 comm_kpt  = mpi_enreg%comm_cell
 if (paral_kgb == 1) comm_kpt = mpi_enreg%comm_kpt
 comm_band = mpi_enreg%comm_band
 me_kpt  = mpi_enreg%me_kpt
 me_band = mpi_enreg%me_band
 nproc_band  = mpi_enreg%nproc_band
 nproc_spkpt = mpi_enreg%nproc_spkpt

 if (nproc_band /= (mband/mband_cprj)) then
   message = "Inconsistency in datafordmft: nproc_band should be equal to mband/mband_cprj"
   ABI_BUG(message)
 end if

 iorder_cprj = 0
 ABI_CHECK(dtset%mkmem/=0,"mkmem=0 not supported anymore!")
!todo_ab: extract cprj from file unpaw in the following..
!call abi_abort('COLL')

!----------------------------------- MPI-------------------------------------

 nbandi = paw_dmft%dmftbandi
 nbandf = paw_dmft%dmftbandf
 t2g    = (paw_dmft%dmft_t2g == 1)
 x2my2d = (paw_dmft%dmft_x2my2d == 1)
 use_full_chipsi = (paw_dmft%dmft_use_full_chipsi == 1)

 if (use_full_chipsi .and. mpi_enreg%nproc_fft > 1) then
   message = "datafordmft not working when nproc_fft > 1 and use_full_chipsi=1"
   ABI_ERROR(message)
 end if

!if(mpi_enreg%me==0) write(7886,*) "in datafordmft", mpi_enreg%me, mpi_enreg%nproc
!if(mpi_enreg%me==1) write(7887,*) "in datafordmft", mpi_enreg%me, mpi_enreg%nproc
!if(mpi_enreg%me==2) write(7888,*) "in datafordmft", mpi_enreg%me, mpi_enreg%nproc
 write(message,'(2a)') ch10,' == Prepare data for DMFT calculation'
 call wrtout(std_out,message,'COLL')
 if (abs(pawprtvol) >= 3) then
   write(message,'(2a)') ch10,'---------------------------------------------------------------'
!  call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message,'(12a)') ch10,'  Print useful data (as a check)',ch10,&
    & '  - Overlap of KS wfc with atomic orbital inside sphere',ch10,&
    & '  - Eigenvalues',ch10,&
    & '  - Weights of k-points',ch10,&
    & '  - Number of spins ',ch10,&
    & '  - Number of states'
!  call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message,'(2a)') ch10,'---------------------------------------------------------------'
   call wrtout(std_out,message,'COLL')
 end if ! abs(pawprtvol)>=3

 if (dtset%nstep == 0 .and. dtset%nbandkss == 0) then
   message = 'nstep should be greater than 1'
   ABI_BUG(message)
 end if

!********************* Max Values for U terms.
 maxlpawu    = paw_dmft%maxlpawu
 maxmeshsize = paw_dmft%maxmeshsize
 maxnproju   = paw_dmft%maxnproju

!*****************   in forlb.eig
 if (paw_dmft%myproc == 0 .and. abs(pawprtvol) >= 3) then
   if (open_file('forlb.eig',message,newunit=unt,form='formatted',status='unknown') /= 0) ABI_ERROR(message)
   rewind(unt)
   write(unt,*) "Number of bands,   spins, and  k-point; and spin-orbit flag"
   write(unt,*) mband,nsppol,nkpt,my_nspinor,nbandi,nbandf
   write(unt,*) " For each k-point, eigenvalues for each band"
   write(unt,*) (dtset%wtk(ikpt),ikpt=1,nkpt)
   band_index = 0
   do isppol=1,nsppol
     write(unt,*) " For spin"
     write(unt,*) isppol
     do ikpt=1,nkpt
       nband_k = dtset%nband(ikpt+(isppol-1)*nkpt)
       ibandc = 0
       write(unt,*) " For k-point"
       write(unt,*) ikpt
       do iband=1,mband
         if (paw_dmft%band_in(iband)) then
           ibandc = ibandc + 1
           write(unt,'(2i6,4x,f20.15)') ibandc,ikpt,eigen(iband+band_index)*two
         end if
       end do ! iband
       band_index = band_index + nband_k
     end do ! ikpt
   end do ! isppol
   close(unt)
 end if ! proc=me

!==   put eigen into eigen_dft
 paw_dmft%eigen => eigen(:)
 band_index = 0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k = dtset%nband(ikpt+(isppol-1)*nkpt)
     ibandc = 0
     do iband=1,mband
       if (paw_dmft%band_in(iband)) then
         ibandc = ibandc + 1
         paw_dmft%eigen_dft(ibandc,ikpt,isppol) = eigen(iband+band_index) ! in Ha
        ! paw_dmft%eigen_dft(isppol,ikpt,ibandc)=fermie
       end if
     end do ! iband
     band_index = band_index + nband_k
   end do ! ikpt
 end do ! isppol

 if (abs(pawprtvol) >= 3) then
   write(message,'(2a)') ch10,' datafordmft : eigenvalues written on file'
   call wrtout(std_out,message,'COLL')
 end if
!==========================================================================
!***************** Compute  <Chi|Psi>= <Chi|Psi_tilde> + \sum_{proja} <P_a|Psi><Chi|phi_a-phi_tilde_a>
!==========================================================================
!write(std_out,*) "size(cprj,dim=1)",size(cprj,dim=1),size(cprj,dim=2),dtset%mband,dtset%mkmem,dtset%nkpt

!Allocate temporary cwaveprj storage
 ABI_MALLOC(cwaveprj,(natom,my_nspinor))
!write(std_out,*) "before alloc cprj"
!write(std_out,*) size(cwaveprj,dim=1),size(cwaveprj,dim=2),size(dimcprj,dim=1)

 call pawcprj_alloc(cwaveprj(:,:),0,dimcprj(:))
!write(std_out,*) "after alloc cprj"

 ABI_MALLOC(cwprj,(maxnproju,2*maxlpawu+1))
 ABI_MALLOC(chipsi_tmp,(maxmeshsize))

 siz_buf = 0
 siz_buf_psi = 0
 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   itypat = paw_dmft%typat(iatom)
   ndim = 2*lpawu + 1
   siz_buf = siz_buf + ndim
   siz_buf_psi = siz_buf_psi + ndim*paw_dmft%radgrid(itypat)%mesh_size
 end do ! iatom
 siz_buf = siz_buf * mbandc * nspinor * mkmem * nsppol_mem
 siz_buf_psi = siz_buf_psi * mbandc * nspinor * mkmem * nsppol_mem

 ABI_MALLOC(recvcounts,(nproc_spkpt))
 ABI_MALLOC(displs,(nproc_spkpt))
 call xmpi_allgather(siz_buf,recvcounts(:),comm_kpt,ierr)

 displs(1) = 0
 do irank=2,nproc_spkpt
   displs(irank) = displs(irank-1) + recvcounts(irank-1)
 end do ! irank

 ABI_MALLOC(buf_chipsi,(siz_buf))
 ABI_MALLOC(buf_chipsi_tot,(recvcounts(nproc_spkpt)+displs(nproc_spkpt)))
 if (prt_wan) then
   ABI_MALLOC(paw_dmft%buf_psi,(siz_buf_psi))
 end if

 icprj = 0
 icg = 0
 ibuf_psi = 0
 ibuf_chipsi = 0

 buf_chipsi(:) = czero

 do isppol=1,nsppol

   if (mpi_enreg%my_isppoltab(isppol) == 0) cycle
   ik = 0

   do ikpt=1,nkpt

     nband_k = dtset%nband(ikpt+(isppol-1)*nkpt)
     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_kpt)) cycle

     nband_k_cprj = nband_k / nproc_band
     npw = paw_dmft%npwarr(ikpt)
     icgb = icg
     ik = ik + 1
     ibandc = 0
!    LOOP OVER BANDS
     ib = 0

     do iband=1,nband_k

       ! Parallelization: treat only some bands
       verif_band = .true.
       if (paral_kgb == 1) then
         if (mod((iband-1)/mpi_enreg%bandpp,nproc_band) /= me_band) verif_band = .false.
       else
         if (mpi_enreg%proc_distrb(ikpt,iband,isppol) /= me_kpt) verif_band = .false.
       end if

       if (verif_band) ib = ib + 1

       if (paw_dmft%band_in(iband)) then
         ibandc = ibandc + 1
       else
         icgb = icgb + npw*nspinor
         cycle
       end if

       if (verif_band) then
         call pawcprj_get(cryst_struc%atindx1(:),cwaveprj(:,:),cprj(:,:),natom,ib,icprj,ikpt,&
                        & iorder_cprj,isppol,mband_cprj,dtset%mkmem,natom,1,nband_k_cprj,&
                        & my_nspinor,nsppol,paw_dmft%unpaw,mpicomm=mpi_enreg%comm_kpt,&
                        & proc_distrb=mpi_enreg%proc_distrb(:,:,:))
       end if

       do ispinor=1,my_nspinor
         do iatom=1,natom
           lpawu = paw_dmft%lpawu(iatom)
           if (lpawu == -1) cycle
           lpawu1 = lpawu
           if (t2g .or. x2my2d) lpawu1 = 2
           itypat = paw_dmft%typat(iatom)
           ndim = 2*lpawu + 1
           nproju = pawtab(itypat)%nproju
           siz_proj = paw_dmft%siz_proj(itypat)
           siz_wan = paw_dmft%radgrid(itypat)%mesh_size
           siz_paw = min(siz_wan,paw_dmft%int_meshsz(itypat))
           rint = paw_dmft%radgrid(itypat)%rad(siz_proj)

           if (verif_band) then
             lmn_size = pawtab(itypat)%lmn_size
             do ilmn=1,lmn_size
               ! ------------ Select l=lpawu.
               if (pawtab(itypat)%indlmn(1,ilmn) /= lpawu1) cycle

               im = pawtab(itypat)%indlmn(2,ilmn) + lpawu1 + 1
               if (x2my2d) then
                 if (im /= 5) cycle
                 im = 1
               else if (t2g) then
                 if (im == 3 .or. im == 5) cycle
                 if (im == 4) im = 3
               end if
               iproj = pawtab(itypat)%indlmn(3,ilmn)
               cwprj(iproj,im) = cmplx(cwaveprj(iatom,ispinor)%cp(1,ilmn), &
                                     & cwaveprj(iatom,ispinor)%cp(2,ilmn),kind=dp)
             end do ! ilmn
           end if ! verif

           do im=1,ndim

             ibuf_chipsi = ibuf_chipsi + 1

             if (use_full_chipsi) then

               buf_chipsi(ibuf_chipsi) = sum(cmplx(cg(1,icgb+1:icgb+npw),cg(2,icgb+1:icgb+npw),kind=dp)* &
                                           & paw_dmft%dpro(1:npw,iatom,ik)*paw_dmft%ylm(1:npw,im,lpawu+1,ik)* &
                                           & paw_dmft%bessel_int(1:npw,itypat,ik))

               if (prt_wan) then

                 do ir=1,siz_wan

                   ! Compute <Ylm|Psi_tilde>(r) = sum_g c_g * <Ylm|exp(j*(k+G)*(r+Rat))> / sqrt(ucvol)
                   ! using exp(j*(k+G)*r) = 4*pi*sum_{lm} j**l * jl(|k+G|*r) * ylm(k+G) * ylm(theta,phi)
                   ! (spherical harmonics expansion of planewave)

                   paw_dmft%buf_psi(ibuf_psi+ir) = sum(paw_dmft%dpro(1:npw,iatom,ik) * paw_dmft%bessel(1:npw,ir,itypat,ik) * &
                                                     & cmplx(cg(1,icgb+1:icgb+npw),cg(2,icgb+1:icgb+npw),kind=dp) * &
                                                     & paw_dmft%ylm(1:npw,im,lpawu+1,ik))

                 end do ! ir

               end if

               if (verif_band) then
                 do iproj=1,nproju
                   buf_chipsi(ibuf_chipsi) = buf_chipsi(ibuf_chipsi) + cwprj(iproj,im)*paw_dmft%phimtphi_int(iproj,itypat)
                   if (prt_wan) paw_dmft%buf_psi(ibuf_psi+1:ibuf_psi+siz_paw) = paw_dmft%buf_psi(ibuf_psi+1:ibuf_psi+siz_paw) + &
                                                          & cwprj(iproj,im)*paw_dmft%phimtphi(1:siz_paw,iproj,itypat)
                 end do ! iproj
               end if ! verif

             else

               ! In that case, simply assume |Psi> = \sum_{proja} <P_a|Psi><Chi|phi_a> (only true inside the PAW sphere
               ! and if your PAW basis is complete)
               ! Do not use DOT_PRODUCT
               if (verif_band) buf_chipsi(ibuf_chipsi) = buf_chipsi(ibuf_chipsi) + sum(cwprj(1:nproju,im)*paw_dmft%phi_int(1:nproju,itypat))

             end if ! use_full_chipsi

             ibuf_psi = ibuf_psi + siz_wan

           end do ! im
         end do ! iatom
         icgb = icgb + npw
       end do ! ispinor
     end do ! iband

     icprj = icprj + nband_k_cprj*nspinor
     icg = icg + nband_k*npw*nspinor
   end do ! ikpt
 end do ! isppol

!do isppol=1,nsppol
!do ikpt=1,nkpt
!do ispinor=1,my_nspinor
!write(std_out,*) "psichi integers",isppol,ikpt,ispinor
!write(std_out,*) "psichi IB3 iAT1 IM1",&
!&             real(paw_dmft%psichi(isppol,ikpt,3,ispinor,1,1)), imag(paw_dmft%psichi(isppol,ikpt,3,ispinor,1,1))
!
!enddo
!enddo
!enddo
!call abi_abort('COLL')
 !if (abs(pawprtvol) >= 3) then
 !  write(message,*) "chinorm used here =",chinorm
 !  call wrtout(std_out,message,'COLL')
 !end if

!deallocate temporary cwaveprj/cprj storage
 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

 ABI_FREE(cwprj)
 ABI_FREE(chipsi_tmp)

!==========================================================================
!********************* Gather information for MPI before printing
!==========================================================================

 !call xmpi_barrier(comm_kpt)
 if (paral_kgb == 1 .and. nproc_band > 1) then
   ! Build sum over band processors
   call xmpi_sum(buf_chipsi(:),comm_band,ierr)
 end if
 call xmpi_allgatherv(buf_chipsi(:),siz_buf,buf_chipsi_tot(:),recvcounts(:),displs(:),comm_kpt,ierr)

 ABI_FREE(displs)
 ABI_FREE(recvcounts)
 ABI_FREE(buf_chipsi)

 ! Reorder the chipsi since the kpts can have any arbitrary distribution over the different MPI processes
 ibuf_chipsi = 0
 do irank=0,nproc_spkpt-1
   do isppol=1,nsppol
     do ikpt=1,nkpt

       nband_k = dtset%nband(ikpt+(isppol-1)*nkpt)
       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,irank)) cycle

       do ibandc=1,mbandc
         do ispinor=1,nspinor
           do iatom=1,natom
             lpawu = paw_dmft%lpawu(iatom)
             if (lpawu == -1) cycle
             ndim = 2*lpawu + 1
             paw_dmft%chipsi(1+(ispinor-1)*ndim:ispinor*ndim,ibandc,ikpt,isppol,iatom) = buf_chipsi_tot(ibuf_chipsi+1:ibuf_chipsi+ndim)
             ibuf_chipsi = ibuf_chipsi + ndim
           end do ! iatom
         end do ! ispinor
       end do ! ibandc

     end do ! ikpt
   end do ! isppol
 end do ! irank

 !call xmpi_barrier(comm_kpt)

 ABI_FREE(buf_chipsi_tot)

!do isppol=1,nsppol
!do ikpt=1,nkpt
!do ibandc=1,paw_dmft%mbandc
!do ispinor=1,my_nspinor
!write(std_out,*) "psichigather",isppol,ikpt,ibandc,&
!&             real(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,1,1))**2+&
!&             imag(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,1,1))**2
!
!enddo
!enddo
!enddo
!enddo

!if(mpi_enreg%me.eq.0) write(177,*) "end",psichi
!if(mpi_enreg%me.eq.1) write(178,*) "end",psichi
!if(mpi_enreg%me.eq.2) write(179,*) "end",psichi

!==========================================================================
!********* WRITE unnormalized chipsi in file for reference
!==========================================================================
 if (paw_dmft%myproc == 0) then
   call chipsi_print(paw_dmft,pawtab(:))
 end if ! proc=0

!********************* Check normalization and occupations ***************
! Only if the complete BZ is sampled (ie paw_dmft%kspectralfunc=0)
!==========================================================================
 if (paw_dmft%dmft_kspectralfunc == 0) then

   call init_oper(paw_dmft,loc_norm_check,opt_ksloc=2)
   call chipsi_check(paw_dmft,dft_occup,loc_norm_check)
!==========================================================================
!***************  write checks  *******************************************
!==========================================================================
   !if (abs(pawprtvol) >= 3) then
   !  write(message,*) "normalization computed"
   !  call wrtout(std_out,message,'COLL')
   !end if

   write(message,'(2a)') ch10," == The DMFT orbitals are now projected on the correlated bands"
   call wrtout(std_out,message,'COLL')

   write(message,'(2a,i4)') ch10," == Check: Downfolded Occupations and Norm of unnormalized projected orbitals"
   call wrtout(std_out,message,'COLL')

   if (paw_dmft%dmftcheck >= 1) then
     ! print occupations
     write(message,'(2a,i4)') ch10,'  ------ Unsymmetrized Occupations'
     call wrtout(std_out,message,'COLL')

     call print_matlu(dft_occup%matlu(:),natom,pawprtvol)

     ! print norms
     write(message,'(2a,i4)') ch10,'  ------ Unsymmetrized Norm'
     call wrtout(std_out,message,'COLL')

     call print_matlu(loc_norm_check%matlu(:),natom,pawprtvol)
   end if ! dmftcheck>=1

   ! symmetrize and print occupations
   call sym_matlu(dft_occup%matlu(:),paw_dmft)

   write(message,'(2a,i4)') ch10,'  ------ Symmetrized Occupations'
   call wrtout(std_out,message,'COLL')

   call print_matlu(dft_occup%matlu(:),natom,pawprtvol)

   ! symmetrize and print norms
   call sym_matlu(loc_norm_check%matlu(:),paw_dmft)

   write(message,'(2a,i4)') ch10,'  ------ Symmetrized Norm'
   call wrtout(std_out,message,'COLL')

   call print_matlu(loc_norm_check%matlu(:),natom,pawprtvol)

   ! Tests density matrix DFT+U and density matrix computed here.
   if (paw_dmft%dmftcheck == 2 .or. paw_dmft%dmftbandi == 1) then
     ABI_MALLOC(matlu_temp,(natom))
     call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu_temp(:))
     isppol   = 1
     ispinor  = 1
     ispinor1 = 1
     nsploop = max(nsppol,nspinor**2)
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       ndim = 2*lpawu + 1
       do idijeff=1,nsploop
         if (nsploop <= 2) then
           isppol = idijeff
         else if (nsploop == 4) then
           ispinor  = spinor_idxs(1,idijeff)
           ispinor1 = spinor_idxs(2,idijeff)
         else
           ABI_BUG(" BUG in datafordmft: nsploop should be equal to 2 or 4")
         end if ! nsploop
         if (my_nspinor == 2) matlu_temp(iatom)%mat(1+(ispinor-1)*ndim:ispinor*ndim,1+(ispinor1-1)*ndim:ispinor1*ndim,isppol) = &
               & cmplx(paw_ij(iatom)%noccmmp(1,1:ndim,1:ndim,idijeff),paw_ij(iatom)%noccmmp(2,1:ndim,1:ndim,idijeff),kind=dp)
         if (my_nspinor == 1) matlu_temp(iatom)%mat(1+(ispinor-1)*ndim:ispinor*ndim,1+(ispinor1-1)*ndim:ispinor1*ndim,isppol) = &
               & cmplx(paw_ij(iatom)%noccmmp(1,1:ndim,1:ndim,idijeff),zero,kind=dp)
       end do ! idijeff
     end do ! iatom
     if (paw_dmft%dmftcheck == 2) option = 1
     if (paw_dmft%dmftcheck <= 1) option = 0
     call diff_matlu("DFT+U density matrix from INPUT wfk",&
       & "Direct calculation of density matrix with chipsi from DIAGONALIZED wfk",&
       & matlu_temp(:),dft_occup%matlu(:),natom,option,tol3,ierr) !tol1 tol2 tol3
     if (ierr == -1) then
       write(message,'(10a)') ch10,&
        & '    -> These two quantities should agree if three conditions are fullfilled',ch10,&
        & '         -  input wavefunctions come from the same Hamiltonian (e.g LDA/GGA)',ch10,&
        & '         -  dmatpuopt is equal to 1',ch10,&
        & '         -  all valence states are in the valence',ch10,&
        & '    (for experts users: it is not compulsory that these conditions are fullfilled)'
       call wrtout(std_out,message,'COLL')
     end if
!    write(message,'(2a)') ch10,&
!    &   '  ***** => Calculations of density matrices with projections and in DFT+U are coherent****'
!    call wrtout(std_out,message,'COLL')

     call destroy_matlu(matlu_temp(:),natom)
     ABI_FREE(matlu_temp)
   else
     write(message,'(2a)') ch10,&
       & '  Warning: Consistency of density matrices computed from projection has not been checked: use dmftcheck>=2 '
     call wrtout(std_out,message,'COLL')
   end if

   call destroy_oper(loc_norm_check)
 end if ! dmft_kspectralfunc=0

 if (present(nbandkss)) then
   if ((me_kpt == 0 .and. nbandkss /= 0) .or. (paw_dmft%dmft_kspectralfunc == 1)) then
!     opt_renorm=1 ! if ucrpa==1, no need for individual orthonormalization
     opt_renorm = 3
     if (dtset%ucrpa >= 1 .or. paw_dmft%dmft_kspectralfunc == 1) opt_renorm = 2
     call chipsi_renormalization(paw_dmft,opt=opt_renorm)
     if (paw_dmft%myproc == 0) then
       call chipsi_print(paw_dmft,pawtab(:))
     end if
   end if ! proc=me
 end if

 CONTAINS

!!****f* m_datafordmft/chipsi_check
!! NAME
!!  chipsi_check
!!
!! FUNCTION
!!  Check chipsi: compute norm and occupations
!!
!! INPUTS
!!  paw_dmft <type(paw_dmft)>=paw data for the self-consistency
!!
!!  OUTPUTS:
!!  xocc_check: density matrix
!!  xnorm_check: matrix of norms
!!
!! SIDE EFFECTS
!!
!! SOURCE

subroutine chipsi_check(paw_dmft,xocc_check,xnorm_check)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(oper_type), intent(inout) :: xnorm_check,xocc_check
!Local variables ------------------------------------
 real(dp), allocatable :: occ_dft(:,:,:)
! *********************************************************************

 ABI_MALLOC(occ_dft,(mbandc,nkpt,nsppol))

 band_index = 0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k = paw_dmft%nband(ikpt+(isppol-1)*nkpt)
     ibandc = 0
     do iband=1,nband_k
       if (paw_dmft%band_in(iband)) then
         ibandc = ibandc + 1
         occ_dft(ibandc,ikpt,isppol) = occ(iband+band_index)
       end if
     end do ! iband
     band_index = band_index + nband_k
   end do ! ikpt
 end do ! isppol

 if (nsppol == 1 .and. my_nspinor == 1) occ_dft(:,:,:) = occ_dft(:,:,:) * half

 call downfold_oper(xnorm_check,paw_dmft,procb=paw_dmft%distrib%procb(:),iproc=paw_dmft%distrib%me_kpt,option=2)
 call xmpi_matlu(xnorm_check%matlu(:),natom,paw_dmft%distrib%comm_kpt)
 call downfold_oper(xocc_check,paw_dmft,procb=paw_dmft%distrib%procb(:), &
                  & iproc=paw_dmft%distrib%me_kpt,option=3,op_ks_diag=occ_dft(:,:,:))
 call xmpi_matlu(xocc_check%matlu(:),natom,paw_dmft%distrib%comm_kpt)

 ABI_FREE(occ_dft)

 end subroutine chipsi_check
!DBG_EXIT("COLL")
!!***
end subroutine datafordmft
!!***

!!****f* m_datafordmft/chipsi_print
!! NAME
!!  chipsi_print
!!
!! FUNCTION
!!  Print chipsi for reference
!!
!! INPUTS
!!  paw_dmft <type(paw_dmft)>=paw data for the self-consistency
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! SIDE EFFECTS
!!  print chipsi in forlb.ovlp
!!
!! SOURCE

subroutine chipsi_print(paw_dmft,pawtab)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(pawtab_type), intent(in) :: pawtab(paw_dmft%ntypat)
!Local variables ------------------------------------
 integer :: iatom,iband,ibandc,ikpt,im,im1,ispinor,isppol
 integer :: itypat,lpawu,nband_k,ndim,unt
 logical :: t2g,x2my2d
 character(len=500) :: msg
 integer, parameter :: mt2g(3) = (/1,2,4/)
! *********************************************************************

   t2g = (paw_dmft%dmft_t2g == 1)
   x2my2d = (paw_dmft%dmft_x2my2d == 1)

   if (open_file('forlb.ovlp',msg,newunit=unt,form='formatted',status='unknown') /= 0) ABI_ERROR(msg)
   rewind(unt)

!  Header for calc_uCRPA.F90
   write(unt,*) "# isppol   nspinor   natom   m    Re(<chi|psi>)   Im(<chi|psi>)"
   if (count(pawtab(:)%lpawu /= -1) == 1) then
     do itypat=1,paw_dmft%ntypat
       lpawu = pawtab(itypat)%lpawu
       if (lpawu == -1) cycle
       if (t2g) then
         write(unt,*) "l= ",1,itypat
       else if (x2my2d) then
         write(unt,*) "l= ",0,itypat
       else
         write(unt,*) "l= ",lpawu,itypat
       end if
     end do ! itypat
   else
     write(unt,*) "More than one correlated species"
   end if

   write(unt,*) "Bands ",paw_dmft%dmftbandi,paw_dmft%dmftbandf

   do isppol=1,paw_dmft%nsppol
     do ikpt=1,paw_dmft%nkpt
!      rewind(1023)
       write(unt,'(a6,2x,i6)') "ikpt =",ikpt
       nband_k = paw_dmft%nband(ikpt+(isppol-1)*paw_dmft%nkpt)
       ibandc = 0
       do iband=1,nband_k
         if (paw_dmft%band_in(iband)) then
           ibandc = ibandc + 1
           write(unt,'(a8,2x,i6)') " iband =",iband
         else
           cycle
         end if
         do ispinor=1,paw_dmft%nspinor
           do iatom=1,paw_dmft%natom
             lpawu = paw_dmft%lpawu(iatom)
             if (lpawu == -1) cycle
             ndim = 2*lpawu + 1
             do im=1,ndim
               if (t2g) then
                 im1 = mt2g(im)
               else if (x2my2d) then
                 im1 = 5
               else
                 im1 = im
               end if

               write(unt,'(4i6,3x,2f23.15)') isppol,ispinor,iatom,im1,&
                        & dble(paw_dmft%chipsi(im+(ispinor-1)*ndim,ibandc,ikpt,isppol,iatom)),&
                        & aimag(paw_dmft%chipsi(im+(ispinor-1)*ndim,ibandc,ikpt,isppol,iatom))
             end do !im
           end do ! iatom
         end do ! ispinor
       end do !iband
     end do !ikpt
   end do ! isppol
!   write(unt,*) "Fermi level (in Ryd)="
!   write(unt,*) fermie*two
   close(unt)
 end subroutine chipsi_print
!!***

!!****f* m_datafordmft/compute_levels
!! NAME
!! compute_levels
!!
!! FUNCTION
!! Compute correlated electronic levels for ctqmc
!!
!! INPUTS
!!  hdc= double counting
!!  paw_dmft <type(paw_dmft)>=paw data for the self-consistency
!!
!! OUTPUT
!!  energy_level= local electronic levels
!!  nondiag= true if the levels are not diagonal
!!
!! NOTES
!!
!! SOURCE

 subroutine compute_levels(energy_level,hdc,paw_dmft,nondiag)

!Arguments ------------------------------------
 type(oper_type), intent(in) :: hdc
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(oper_type), intent(inout) :: energy_level
 logical, optional, intent(out) :: nondiag
!Local variables ------------------------------
 integer :: iatom,im,isppol,lpawu,natom,ndim
 character(len=13) :: tag
 character(len=500) :: message
!************************************************************************

 natom = paw_dmft%natom
 if (present(nondiag)) nondiag = .false.

!======================================================================
!Compute atomic levels from projection of \epsilon_{nks} and symmetrize
!======================================================================

 call downfold_oper(energy_level,paw_dmft,option=3,procb=paw_dmft%distrib%procb(:), &
                  & iproc=paw_dmft%distrib%me_kpt,op_ks_diag=paw_dmft%eigen_dft(:,:,:))
 call xmpi_matlu(energy_level%matlu(:),natom,paw_dmft%distrib%comm_kpt)
! write(message,'(a,2x,a,f13.5)') ch10," == Print Energy levels before sym and only DFT"
! call wrtout(std_out,message,'COLL')
! call print_matlu(energy_level%matlu,natom,1)
 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   ndim = paw_dmft%nspinor * (2*lpawu+1)
   do isppol=1,paw_dmft%nsppol
     do im=1,ndim
       energy_level%matlu(iatom)%mat(im,im,isppol) = energy_level%matlu(iatom)%mat(im,im,isppol) &
          & - paw_dmft%fermie
     end do ! im
   end do ! isppol
   energy_level%matlu(iatom)%mat(:,:,:) = energy_level%matlu(iatom)%mat(:,:,:) - hdc%matlu(iatom)%mat(:,:,:)
!    write(std_out,*) "DC,fermie",hdc%matlu(iatom)%mat(1,1,1,1,1),paw_dmft%fermie
 end do ! iatom

 call sym_matlu(energy_level%matlu(:),paw_dmft)
 if (present(nondiag)) then
   call checkdiag_matlu(energy_level%matlu(:),natom,tol7,nondiag)
 end if

 write(tag,'(f13.5)') paw_dmft%fermie
 write(message,'(a,2x,2a)') ch10," == Print Energy levels in cubic basis for Fermi Level= ",adjustl(tag)
 call wrtout(std_out,message,'COLL')
!call print_oper(energy_level,1,paw_dmft,1)
 call print_matlu(energy_level%matlu(:),natom,1)

 end subroutine compute_levels
!!***

!!****f* m_datafordmft/chipsi_renormalization
!! NAME
!! chipsi_renormalization
!!
!! FUNCTION
!! Orthonormalize chipsi.
!!
!! INPUTS
!!  paw_dmft =  data for DFT+DMFT calculations.
!!  opt = 2 : orthonormalize all the atoms at the same time,
!!            and for each individual kpt
!!      = 3 (default) : orthonormalize the sum over all kpt,
!!          and for each individual atom
!!
!! OUTPUT
!!  paw_dmft%chipsi((2*maxlpawu+1)*nspinor,mbandc,nkpt,nsppol,natom):
!!                    orthonormalized projections <Chi|Psi>
!!
!! NOTES
!!
!! SOURCE

subroutine chipsi_renormalization(paw_dmft,opt)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(inout) :: paw_dmft
 integer, optional, intent(in) :: opt
!Local variables ------------------------------
 integer :: jkpt,natom,nkpt,option
 real(dp) :: pawprtvol
 type(oper_type) :: norm,oper_temp
 character(len=500) :: message
 real(dp), allocatable :: wtk_tmp(:)
! real(dp),allocatable :: e0pde(:,:,:),omegame0i(:)
!************************************************************************

 DBG_ENTER("COLL")

 option = 3
 if (present(opt)) then
   if (opt == 2 .or. opt == 3) option = opt
 end if
 pawprtvol = 2

 natom = paw_dmft%natom
 nkpt  = paw_dmft%nkpt

!== Normalize psichi
 !if (option == 1) then
!  ====================================
!  == simply renormalize psichi =======
!  ====================================
 !  write(message,'(2a)') ch10," Psichi are renormalized  "
 !  call wrtout(std_out,  message,'COLL')
 !  do isppol=1,nsppol
 !    do ikpt=1,nkpt
 !      do ib=1,mbandc
 !        do iatom=1,natom
 !          if(paw_dmft%lpawu(iatom).ne.-1) then
 !            ndim=2*paw_dmft%lpawu(iatom)+1
 !            do im=1,ndim
 !              do ispinor=1,nspinor
!                write(std_out,*) "psichi1",paw_dmft%psichi(isppol,ikpt,ib,ispinor,iatom,im)
 !                paw_dmft%psichi(isppol,ikpt,ib,ispinor,iatom,im)=     &
!&                 paw_dmft%psichi(isppol,ikpt,ib,ispinor,iatom,im)/    &
!&                 sqrt(real(norm%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)))
 !              end do ! ispinor
 !            end do ! im
 !          end if
 !        end do ! iatom
 !      end do ! ib
 !    end do ! ikpt
 !  end do ! isppol
!  todo_ab introduce correct orthonormalization in the general case.

 write(message,'(6a)') ch10, &
   &      '  =================================================== ',&
   & ch10,'  == The DMFT orbitals will now be orthonormalized == ',&
   & ch10,'  =================================================== '
 call wrtout(std_out,message,'COLL')

 if (option == 2) then ! option==2
!  ====================================
!  == renormalize k-point after k-point
!  ====================================

   write(message,'(2a,i5)') ch10,' Nb of K-point',nkpt
   call wrtout(std_out,message,'COLL')
   do jkpt=1,nkpt  ! jkpt
     write(message,'(2a,i5)') ch10,'  == Renormalization for K-point',jkpt
     call wrtout(std_out,message,'COLL')
     if (paw_dmft%distrib%procb(jkpt) /= paw_dmft%distrib%me_kpt) cycle
     call normalizechipsi(1,paw_dmft,jkpt=jkpt)
   end do ! jkpt
   write(message,'(2a)') ch10,'  ===== K-points all renormalized'
   call wrtout(std_out,message,'COLL')

 else if (option == 3) then  ! option==3
!  ====================================
!  == renormalize the sum over k-points
!  ====================================
   write(message,'(6a)') ch10, &
     &      '  ====================================== ',&
     & ch10,'  == Renormalization for all K-points == ',&
     & ch10,'  ====================================== '
   call wrtout(std_out,message,'COLL')
   call normalizechipsi(nkpt,paw_dmft)

 end if ! option

 ! Gather contribution from each CPU
 call chipsi_gather(paw_dmft)

!== Change back repr for norm

!===============================================
!==  Compute norm with new chipsi
!===============================================

 write(message,'(2a)') ch10,'  ===== Compute new norm after renormalization'
 call wrtout(std_out,message,'COLL')
 call init_oper(paw_dmft,oper_temp,opt_ksloc=2)
 call identity_oper(oper_temp,2)

 if (paw_dmft%dmft_kspectralfunc == 1) then
   ABI_MALLOC(wtk_tmp,(nkpt))
   wtk_tmp(:) = one
   call init_oper(paw_dmft,norm,nkpt=1,wtk=wtk_tmp(:),opt_ksloc=2)
   do jkpt=1,nkpt ! jkpt
     norm%shiftk = jkpt - 1
     call downfold_oper(norm,paw_dmft,option=2)
     write(message,'(2a,i5)') &
       & ch10,"  == Check: Overlap after renormalization for k-point",jkpt
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm%matlu(:),natom,prtopt=1)
!== Check that norm is now the identity
     call diff_matlu('Overlap after renormalization','Identity',&
       & norm%matlu(:),oper_temp%matlu(:),natom,1,tol6,zero_or_one=1)
   end do ! jkpt
   ABI_FREE(wtk_tmp)
 else !dmft_kspectralfunc
   write(message,'(2a)') ch10,'  ===== Starting downfold'
   call wrtout(std_out,message,'COLL')
   call init_oper(paw_dmft,norm,opt_ksloc=2)
   call downfold_oper(norm,paw_dmft,procb=paw_dmft%distrib%procb(:),iproc=paw_dmft%distrib%me_kpt,option=2)
   call xmpi_matlu(norm%matlu(:),natom,paw_dmft%distrib%comm_kpt)
   write(message,'(2a)') ch10,'  ===== Finished downfold'
   call wrtout(std_out,message,'COLL')

!== Print unsymmetrized norm%matlu with new chipsi
   if (pawprtvol > 2) then
     write(message,'(4a,2a)') &
      &  ch10,"  == Check: Overlap with renormalized chipsi without symmetrization is == "
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm%matlu(:),natom,prtopt=1)
   end if ! pawprtvol>2

!== Symmetrize norm%matlu with new chipsi
   call sym_matlu(norm%matlu(:),paw_dmft)

!== Print symmetrized norm%matlu with new chipsi
   if (pawprtvol > 2) then
     write(message,'(4a,2a)') &
      & ch10,"  == Check: Overlap with renormalized chipsi and symmetrization is =="
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm%matlu(:),natom,prtopt=1,opt_diag=-1)
   end if

!== Check that norm is now the identity
   call diff_oper('Overlap after renormalization','Identity',norm,oper_temp,0,tol6)

 end if ! dmft_kspectralfunc

 call destroy_oper(norm)
 call destroy_oper(oper_temp)

 paw_dmft%lchipsiortho = 1

 DBG_EXIT("COLL")

 CONTAINS
!===========================================================
!!***

!!****f* chipsi_renormalization/normalizechipsi
!! NAME
!!  normalizechipsi
!!
!! FUNCTION
!!  Orthonormalize chipsi
!!
!! INPUTS
!!  nkpt = number of kpt
!!  paw_dmft =  data for DFT+DMFT calculations.
!!  jkpt = if present, index of the kpt to be orthonormalized
!!
!! SIDE EFFECTS
!!
!! SOURCE

subroutine normalizechipsi(nkpt,paw_dmft,jkpt)

!Arguments ------------------------------------
 integer, intent(in) :: nkpt
 integer, optional, intent(in) :: jkpt
 type(paw_dmft_type), intent(inout) :: paw_dmft
!Local variables ------------------------------
 integer :: dimoverlap,dum,iatom,ib,ikpt,isppol,itot,itot1,lpawu,mbandc
 integer :: natom,ndim,ndim_max,nspinor,nsppol,pawprtvol
 type(oper_type) :: norm1,norm2,norm3
 complex(dpc), allocatable :: chipsivect(:,:),largeoverlap(:,:),mat_tmp(:,:)
 character(len=500) :: message
! real(dp),allocatable :: e0pde(:,:,:),omegame0i(:)
 !complex(dpc), allocatable :: wan(:,:,:),sqrtmatinv(:,:),wanall(:)
 !type(coeff2c_type), allocatable :: overlap(:)
!************************************************************************

 mbandc    = paw_dmft%mbandc
 natom     = paw_dmft%natom
 nspinor   = paw_dmft%nspinor
 nsppol    = paw_dmft%nsppol
 ndim_max  = (2*paw_dmft%maxlpawu+1)*nspinor
 pawprtvol = 3

 if (nkpt /= 1 .and. present(jkpt)) ABI_BUG('BUG in chipsi_normalization')

  ! iortho=1
  ! write(6,*) "nkpt, iortho",nkpt,iortho
   !if (natomcor>1) iortho=2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  First case: usual case (this numerically guarantees downfold(upfold)=Id only for one atom and nkpt=1)
 if (.not. present(jkpt)) then ! .and.iortho==1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    *********************************************************************

   call init_oper(paw_dmft,norm1,opt_ksloc=2)

   call downfold_oper(norm1,paw_dmft,procb=paw_dmft%distrib%procb(:),iproc=paw_dmft%distrib%me_kpt,option=2)
   call xmpi_matlu(norm1%matlu(:),natom,paw_dmft%distrib%comm_kpt)
   if (nkpt > 1) then
     call sym_matlu(norm1%matlu(:),paw_dmft)
   end if

   if (pawprtvol > 2) then
     write(message,'(2a)') ch10,'  - Print current norm and overlap (before orthonormalization)'
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm1%matlu(:),natom,prtopt=1,opt_exp=1)
   end if

!    build large overlap matrix
   !write(message,'(2a)') ch10,'  - Overlap (before orthonormalization) -'
   !call wrtout(std_out,message,'COLL')

!  ==-------------------------------------
!  == Start loop over atoms

     !ABI_MALLOC(overlap,(natom))
     !do iatom=1,natom
     !  if(paw_dmft%lpawu(iatom).ne.-1) then
     !    ndim=2*paw_dmft%lpawu(iatom)+1
     !    tndim=nsppol*nspinor*ndim
     !    ABI_MALLOC(overlap(iatom)%value,(tndim,tndim))
     !    overlap(iatom)%value=czero
     !  end if
     !end do
!    ==-------------------------------------

!    built large overlap matrix
     !write(message,'(2a)') ch10,'  - Overlap (before orthonormalization) -'
     !call wrtout(std_out,message,'COLL')
     !call gather_matlu(norm1%matlu,overlap,cryst_struc%natom,option=1,prtopt=1)
     !call destroy_oper(norm1)

   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     ndim = nspinor * (2*lpawu+1)
     ABI_MALLOC(mat_tmp,(ndim,mbandc))

!        == Compute Inverse Square root of overlap : O^{-0.5}
         !do im=1,tndim
         !  do im1=1,tndim
         !    !write(message,'(a,1x,a,e21.14,a,e21.14,a)') "overlap", &
         !    !"(",real(overlap(1)%value(im,im1)),",",aimag(overlap(1)%value(im,im1)),")"
         !    write(6,*) "overlap",overlap(iatom)%value(im,im1)
         !  enddo
         !enddo
         !stop
         !call wrtout(std_out,message,'COLL')
         !if(diag==0) then
         !call invsqrt_matrix(overlap(iatom)%value,tndim,dum)
         !sqrtmatinv=overlap(iatom)%value
     do isppol=1,nsppol
       ! if(diag==0) then
       call invsqrt_matrix(norm1%matlu(iatom)%mat(:,:,isppol),ndim,dum)
       !else
       !  do im1=1,tndim
       !    do im=1,tndim
       !      if (im==im1) then
       !        norm1%matlu(iatom)%mat(im,im1,isppol)=cone/sqrt(norm1%matlu(iatom)%mat(im,im,isppol)
       !      else
       !        norm1%matlu(iatom)%mat(im,im1,isppol)=czero
       !      end if
       !    end do
       ! end do
       !end if

!        == Apply O^{-0.5} on chipsi
       do ikpt=1,nkpt
         if (paw_dmft%distrib%procb(ikpt) /= paw_dmft%distrib%me_kpt) cycle
         call abi_xgemm("n","n",ndim,mbandc,ndim,cone,norm1%matlu(iatom)%mat(:,:,isppol),ndim, &
                      & paw_dmft%chipsi(:,:,ikpt,isppol,iatom),ndim_max,czero,mat_tmp(:,:),ndim)
         paw_dmft%chipsi(1:ndim,:,ikpt,isppol,iatom) = mat_tmp(:,:)
       end do ! ikpt
     end do ! isppol


       !       ABI_MALLOC(wan,(nsppol,nspinor,ndim))
!        write(std_out,*) mbandc,nsppol,nspinor,ndim
!        write(std_out,*)  paw_dmft%psichi(1,1,1,1,1,1)
        ! do ikpt=1,nkpt
        !   do ib=1,mbandc
        !     if(present(jkpt)) then
        !       ikpt1=jkpt
        !     else
        !       ikpt1=ikpt
        !     end if
        !     jc=0
        !     wan=czero
        !     do isppol=1,nsppol
        !       do ispinor=1,nspinor
        !         do im=1,ndim
!                  write(std_out,*) "psichi", paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)
         !          jc=jc+1
         !          jc1=0
         !          do isppol1=1,nsppol
         !            do ispinor1=1,nspinor
         !              do im1=1,ndim
         !                jc1=jc1+1
         !                wan(isppol,ispinor,im)= wan(isppol,ispinor,im) &
!&                         + paw_dmft%psichi(isppol1,ikpt1,ib,ispinor1,iatom,im1)*sqrtmatinv(jc,jc1)
          !             end do ! ispinor1
          !           end do ! isppol1
          !         end do ! im1
          !       end do ! im
          !     end do ! ispinor
          !   end do !  isppol
          !   do isppol=1,nsppol
          !     do ispinor=1,nspinor
          !       do im=1,ndim
          !         paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)=wan(isppol,ispinor,im)
!                  write(std_out,*) "psichi2", paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)
          !       end do ! ispinor
          !     end do ! isppol
          !   end do ! im
          ! end do ! ib
        ! end do ! ikpt
        ! ABI_FREE(wan)
        ! ABI_FREE(sqrtmatinv)

     ABI_FREE(mat_tmp)

   end do ! iatom

!  == End loop over atoms
!  ==-------------------------------------

!  ======================================================================
!  == Check norm with new chipsi.
!  ======================================================================

   call downfold_oper(norm1,paw_dmft,procb=paw_dmft%distrib%procb(:),iproc=paw_dmft%distrib%me_kpt,option=2)
   call xmpi_matlu(norm1%matlu(:),natom,paw_dmft%distrib%comm_kpt)

   if (nkpt > 1) then
     call sym_matlu(norm1%matlu(:),paw_dmft)
   end if

   if (pawprtvol > 2) then
     write(message,'(2a)') ch10,'  - Print new norm after orthonormalization'
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm1%matlu(:),natom,prtopt=1)
   end if

!  ======================================================================
!  == Check that norm-identity is zero
!  ======================================================================
   call init_oper(paw_dmft,norm2,opt_ksloc=2)
   call init_oper(paw_dmft,norm3,opt_ksloc=2)
   call identity_oper(norm2,2)
   call add_matlu(norm1%matlu(:),norm2%matlu(:),norm3%matlu(:),natom,-1)
   call destroy_oper(norm2)
   if (pawprtvol > 2) then
     write(message,'(2a)') ch10,'  - Print new norm minus Identity '
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm3%matlu(:),natom,prtopt=1,opt_exp=1)
   end if
   call destroy_oper(norm3)

   call destroy_oper(norm1)
!    call flush(std_out)           ! debug debug  debug   debug
!    ABI_ERROR("Stop for debugging")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  New implementation, several atoms, general case.
 else if (present(jkpt)) then !.or.iortho==2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   dimoverlap = 0
   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     dimoverlap = dimoverlap + 2*lpawu + 1
      ! write(6,*) "atom, dimoverlap",iatom,dimoverlap,natomcor
   end do ! iatom

   dimoverlap = dimoverlap * nspinor

   ABI_MALLOC(largeoverlap,(dimoverlap,dimoverlap))
   ABI_MALLOC(chipsivect,(dimoverlap,mbandc))
   ABI_MALLOC(mat_tmp,(dimoverlap,mbandc))

!  Big loop over isppol
   do isppol=1,nsppol
     do ib=1,mbandc
       itot = 0
       do iatom=1,natom
         lpawu = paw_dmft%lpawu(iatom)
         if (lpawu == -1) cycle
         ndim = nspinor * (2*lpawu+1)
         chipsivect(itot+1:itot+ndim,ib) = paw_dmft%chipsi(1:ndim,ib,jkpt,isppol,iatom)
         itot = itot + ndim
         !if(itot>dimoverlap) write(std_out,*) "itot>ndim",itot,ndim
         ! write(6,*) "ib,iatom,im,ispinor",ib,iatom,im,ispinor,jkpt
       end do ! iatom
     end do ! ib


!    Calculation of overlap
     call abi_xgemm("n","c",dimoverlap,dimoverlap,mbandc,cone,chipsivect(:,:),dimoverlap,&
                  & chipsivect(:,:),dimoverlap,czero,largeoverlap(:,:),dimoverlap)

    !   largeoverlap=czero
    !   do ib=1,mbandc
    !     do itot=1,dimoverlap
    !       do itot1=1,dimoverlap
    !          largeoverlap(itot,itot1)=largeoverlap(itot,itot1)+ &
!&              psichivect(ib,itot)*conjg(psichivect(ib,itot1))
    !       enddo ! itot1
    !     enddo ! itot
    !   enddo ! ib

!    Math: orthogonalization of overlap
     write(std_out,*) "jkpt=",jkpt
     do itot=1,dimoverlap
       write(std_out,'(100f7.3)') (largeoverlap(itot,itot1),itot1=1,dimoverlap)
     end do
     call invsqrt_matrix(largeoverlap(:,:),dimoverlap,dum)
     write(std_out,*) "jkpt=",jkpt
     do itot=1,dimoverlap
       write(std_out,'(100f7.3)') (largeoverlap(itot,itot1),itot1=1,dimoverlap)
     end do
     write(std_out,*) "jkpt=",jkpt
     do itot=1,dimoverlap
       write(std_out,'(100e9.3)') (largeoverlap(itot,itot1),itot1=1,dimoverlap)
     end do

     call abi_xgemm("n","n",dimoverlap,mbandc,dimoverlap,cone,largeoverlap(:,:),dimoverlap,&
                  & chipsivect(:,:),dimoverlap,czero,mat_tmp(:,:),dimoverlap)

     !  do ib=1,mbandc
     !    wanall=czero
     !    do itot=1,dimoverlap
     !      do itot1=1,dimoverlap
     !         wanall(itot)= wanall(itot)+psichivect(ib,itot1)*sqrtmatinv(itot,itot1)
      !     enddo ! itot1
     !      write(std_out,'(3i3,2x,i3,2x,2e15.5,2x,2e15.5)') jkpt,isppol,ib,itot,psichivect(ib,itot),wanall(itot)
      !   enddo ! itot
      !   iatomcor=0
      !   do itot=1,dimoverlap
      !     psichivect(ib,itot)=wanall(itot)
      !   enddo
        ! do iatom=1,natom
        !   if(paw_dmft%lpawu(iatom).ne.-1) then
        !     ndim=2*paw_dmft%lpawu(iatom)+1
        !     iatomcor=iatomcor+1
        !     do im=1,ndim
        !       do ispinor=1,nspinor
        !         paw_dmft%psichi(isppol,jkpt,ib,ispinor,iatom,im)=wanall(iatomcor,isppol,ispinor,im)
        !       end do ! ispinor
        !     end do ! im
        !   endif
        ! enddo ! iatom
      ! enddo ! ib


!    Calculation of overlap (check)
     call abi_xgemm("n","c",dimoverlap,dimoverlap,mbandc,cone,mat_tmp(:,:),dimoverlap,&
                  & mat_tmp(:,:),dimoverlap,czero,largeoverlap(:,:),dimoverlap)

      ! largeoverlap=czero
      ! do ib=1,mbandc
      !   do itot=1,dimoverlap
      !     do itot1=1,dimoverlap
      !        largeoverlap(itot,itot1)=largeoverlap(itot,itot1)+ &
!&              psichivect(ib,itot)*conjg(psichivect(ib,itot1))
      !     enddo ! itot1
      !   enddo ! itot
      ! enddo ! ib

     write(std_out,*) "jkpt=",jkpt
     do itot=1,dimoverlap
       write(std_out,'(100f7.3)') (largeoverlap(itot,itot1),itot1=1,dimoverlap)
     end do

!    chipsivect -> chipsi
     do ib=1,mbandc
       itot = 0
       do iatom=1,natom
         lpawu = paw_dmft%lpawu(iatom)
         if (lpawu == -1) cycle
         ndim = nspinor * (2*lpawu+1)
         paw_dmft%chipsi(1:ndim,ib,jkpt,isppol,iatom) = mat_tmp(itot+1:itot+ndim,ib)
         itot = itot + ndim
       end do ! iatom
     end do ! ib

!  End big loop over isppol
   end do !isppol

   ABI_FREE(chipsivect)
   ABI_FREE(largeoverlap)
   ABI_FREE(mat_tmp)

 end if ! option

 end subroutine normalizechipsi
!!***

!!****f* chipsi_renormalization/chipsi_gather
!! NAME
!! chipsi_gather
!!
!! FUNCTION
!! Gather chipsi from every CPU (parallelization over kpts).
!!
!! INPUTS
!!  paw_dmft =  data for DFT+DMFT calculations.
!!
!! OUTPUT
!!
!! NOTES
!!
!! SOURCE

subroutine chipsi_gather(paw_dmft)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(inout) :: paw_dmft
!Local variables ------------------------------
 integer :: iatom,ib,ibuf,ierr,ikpt,irank,isppol,lpawu,mbandc,me_kpt
 integer :: mkmem,natom,ndim,nproc,nspinor,nsppol,shift,siz_buf
 integer, allocatable :: displs(:),recvcounts(:)
 complex(dpc), allocatable :: buffer(:),buffer_tot(:)
!************************************************************************

 me_kpt  = paw_dmft%distrib%me_kpt
 mkmem   = paw_dmft%distrib%nkpt_mem(me_kpt+1)
 mbandc  = paw_dmft%mbandc
 natom   = paw_dmft%natom
 nproc   = paw_dmft%nproc
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 shift   = paw_dmft%distrib%shiftk

 siz_buf = 0
 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   siz_buf = siz_buf + 2*lpawu + 1
 end do ! iatom
 siz_buf = siz_buf * nspinor * mbandc * nsppol

 ABI_MALLOC(recvcounts,(nproc))
 ABI_MALLOC(displs,(nproc))

 recvcounts(:) = siz_buf * paw_dmft%distrib%nkpt_mem(:)
 displs(1) = 0
 do irank=2,nproc
   displs(irank) = displs(irank-1) + recvcounts(irank-1)
 end do ! irank

 ABI_MALLOC(buffer,(recvcounts(me_kpt+1)))
 ABI_MALLOC(buffer_tot,(recvcounts(nproc)+displs(nproc)))

 ibuf = 0
 do ikpt=1,mkmem
   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     ndim = nspinor * (2*lpawu+1)
     do isppol=1,nsppol
       do ib=1,mbandc
         buffer(ibuf+1:ibuf+ndim) = paw_dmft%chipsi(1:ndim,ib,ikpt+shift,isppol,iatom)
         ibuf = ibuf + ndim
       end do ! ib
     end do ! isppol
   end do ! iatom
 end do ! ikpt

 call xmpi_allgatherv(buffer(:),recvcounts(me_kpt+1),buffer_tot(:),recvcounts(:),displs(:),paw_dmft%distrib%comm_kpt,ierr)

 ibuf = 0
 do ikpt=1,paw_dmft%nkpt
   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     ndim = nspinor * (2*lpawu+1)
     do isppol=1,nsppol
       do ib=1,mbandc
         paw_dmft%chipsi(1:ndim,ib,ikpt,isppol,iatom) = buffer_tot(ibuf+1:ibuf+ndim)
         ibuf = ibuf + ndim
       end do ! ib
     end do ! isppol
   end do ! iatom
 end do ! ikpt

 ABI_FREE(recvcounts)
 ABI_FREE(displs)
 ABI_FREE(buffer)
 ABI_FREE(buffer_tot)

 end subroutine chipsi_gather
!!***

end subroutine chipsi_renormalization
!!***

!!****f* m_datafordmft/hybridization_asymptotic_coefficient
!! NAME
!! hybridization_asymptotic_coefficient
!!
!! FUNCTION
!! Compute some components for the limit of hybridization
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  dft_occup
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!  pawtab <type(pawtab)>
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!
!! NOTES
!!
!! SOURCE

subroutine hybridization_asymptotic_coefficient(cryst_struc,paw_dmft,hybri_coeff)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(paw_dmft_type), intent(in) :: paw_dmft
 !type(pawang_type), intent(in) :: pawang
 type(matlu_type), intent(inout) :: hybri_coeff(paw_dmft%natom)
!Local variables ------------------------------
 type(oper_type)  :: ham_a
 type(oper_type)  :: ham_b
 type(oper_type)  :: ham_squarelocal
 type(oper_type)  :: ham_squareks
 integer :: iband1,iband2,ikpt,isppol
!************************************************************************

! call init_oper(paw_dmft,self_minus_hdc)
 call init_oper(paw_dmft,ham_a)
 call init_oper(paw_dmft,ham_b)
 call init_oper(paw_dmft,ham_squareks)
 call init_oper(paw_dmft,ham_squarelocal)

! Create self_minus_hdc%matlu = Sigma-Hdc in local orbitals
! call add_matlu(self%oper(paw_dmft%dmft_nwlo)%matlu,self%hdc%matlu,&
!&             self_minus_hdc%matlu,cryst_struc%natom,-1)

!   write(message,'(a,2x,a)') ch10,        "  == self_minus_hdc (1)"
!   call wrtout(std_out,message,'COLL')
!   call print_matlu(self_minus_hdc%matlu,paw_dmft%natom,1,opt_exp=1)

!! Create self_minus_hdc%ks = Upfold Sigma-Hdc
! call upfold_oper(self_minus_hdc,paw_dmft,1)
! call loc_oper(self_minus_hdc,paw_dmft,1)

!   write(message,'(a,2x,a)') ch10,        "  == self_minus_hdc (2)"
!   call wrtout(std_out,message,'COLL')
!   call print_matlu(self_minus_hdc%matlu,paw_dmft%natom,1,opt_exp=1)

! Create ham_a%ks = H_ks  in KS basis
!----------------------------------------------------
 do iband1=1,paw_dmft%mbandc
   do iband2=1,paw_dmft%mbandc
     do ikpt=1,paw_dmft%nkpt
       do isppol=1,paw_dmft%nsppol
         if(iband1==iband2) then
           ham_a%ks(iband1,iband2,ikpt,isppol) = cmplx(paw_dmft%eigen_dft(iband1,ikpt,isppol),0.d0,kind=dp)
         else
           ham_a%ks(iband1,iband2,ikpt,isppol) = czero
         end if
       end do
     end do
   end do
 end do

! Create ham_a%matlu = H_ks in local orbitals
!---------------------------------------------
 call downfold_oper(ham_a,paw_dmft)

! Symetrise the local quantity (energy levels)
!---------------------------------------------
 call sym_matlu(ham_a%matlu,paw_dmft)

! Create ham_b%ks : Duplicate both ks and local part of ham_a into ham_b
!-----------------------------------------------------------------------
 call copy_oper(ham_a,ham_b)

! Compute ham_squareks%ks   : Make a product of the two KS version
!------------------------------------------------------------------
 call prod_oper(ham_a,ham_b,ham_squareks,1)

! Compute ham_squareks%matlu
!---------------------------
 call downfold_oper(ham_squareks,paw_dmft)

! Symetrise ham_squareks%matlu
!------------------------------
 call sym_matlu(ham_squareks%matlu(:),paw_dmft)

!   write(message,'(a,2x,a)') ch10,        "  == squareks"
!   call wrtout(std_out,message,'COLL')
!   call print_matlu(ham_squareks%matlu,paw_dmft%natom,1,opt_exp=1)


! Compute ham_squarelocal%matlu
!-------------------------------
 call prod_oper(ham_a,ham_b,ham_squarelocal,2)

! Compute the product in local orbitals
!--------------------------------------
 call sym_matlu(ham_squarelocal%matlu(:),paw_dmft)

!   write(message,'(a,2x,a)') ch10,        "  == squarelocal"
!   call wrtout(std_out,message,'COLL')
!   call print_matlu(ham_squarelocal%matlu,paw_dmft%natom,1,opt_exp=1)

! The difference of ham_squareks and ham_squarelocal
! gives the coefficient that we are looking for ( such that F_ij(iw_n) = -C_ij/(iw_n) ).
!----------------------------------------------------------------------------------------
 call add_matlu(ham_squareks%matlu(:),ham_squarelocal%matlu(:),hybri_coeff,cryst_struc%natom,-1)

 !  write(message,'(a,2x,a)') ch10,        "  == Coeff C_ij before sym"
 !  call wrtout(std_out,message,'COLL')
 !  call print_matlu(hybri_coeff,paw_dmft%natom,1,opt_exp=1)

! Symetrise the local quantity
!------------------------------
 call sym_matlu(hybri_coeff,paw_dmft)

 !  write(message,'(a,2x,a)') ch10,        "  == Coeff C_ij after sym"
 !  call wrtout(std_out,message,'COLL')
 !  call print_matlu(hybri_coeff,paw_dmft%natom,1,opt_exp=1)

 call destroy_oper(ham_squarelocal)
! call destroy_oper(self_minus_hdc)
 call destroy_oper(ham_a)
 call destroy_oper(ham_b)
 call destroy_oper(ham_squareks)


end subroutine hybridization_asymptotic_coefficient
!!***

!!****f* m_datafordmft/compute_wannier
!! NAME
!! compute_wannier
!!
!! FUNCTION
!!  Compute the projected Wannier function in real space.
!!
!! INPUTS
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!  mpi_enreg=information about MPI parallelization
!!
!! OUTPUT
!!
!! NOTES
!!
!! SOURCE

subroutine compute_wannier(paw_dmft,mpi_enreg)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(MPI_type), intent(in) :: mpi_enreg
!Local variables ------------------------------
 integer :: iatom,ibuf_psi,ibandc,ierr,iflavor,ik,ikpt,im,ispinor,isppol
 integer :: itypat,lpawu,natom,nband_k,ndim,nkpt,nspinor,nsppol,siz_wan
!************************************************************************

 natom   = paw_dmft%natom
 nkpt    = paw_dmft%nkpt
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol

 ABI_MALLOC(paw_dmft%wannier,(paw_dmft%maxmeshsize,nspinor*(2*paw_dmft%maxlpawu+1)*nsppol,natom))
 paw_dmft%wannier(:,:,:) = czero
 ibuf_psi = 0

 do isppol=1,nsppol

   if (mpi_enreg%my_isppoltab(isppol) == 0) cycle
   ik = 0

   do ikpt=1,nkpt

     nband_k = paw_dmft%nband(ikpt+(isppol-1)*nkpt)
     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,mpi_enreg%me_kpt)) cycle

     ik = ik + 1

     do ibandc=1,paw_dmft%mbandc
       do ispinor=1,nspinor
         do iatom=1,natom

           lpawu = paw_dmft%lpawu(iatom)
           if (lpawu == -1) cycle
           itypat = paw_dmft%typat(iatom)
           ndim = 2*lpawu + 1
           siz_wan = paw_dmft%radgrid(itypat)%mesh_size

           do im=1,ndim

             iflavor = im + (isppol*ispinor-1)*ndim

             paw_dmft%wannier(1:siz_wan,iflavor,iatom) = paw_dmft%wannier(1:siz_wan,iflavor,iatom) + &
                   & conjg(paw_dmft%chipsi(im+(ispinor-1)*ndim,ibandc,ikpt,isppol,iatom))* &
                   & paw_dmft%buf_psi(ibuf_psi+1:ibuf_psi+siz_wan)*paw_dmft%wtk(ikpt)

             ibuf_psi = ibuf_psi + siz_wan

           end do ! im

         end do ! iatom
       end do ! ispinor
     end do ! ibandc

   end do ! ikpt

 end do ! isppol

 ! No need to broadcast on every CPU
 if (mpi_enreg%paral_kgb == 1 .and. mpi_enreg%nproc_band > 1) then
   call xmpi_sum_master(paw_dmft%wannier(:,:,:),0,mpi_enreg%comm_band,ierr)
 end if
 call xmpi_sum_master(paw_dmft%wannier(:,:,:),0,mpi_enreg%comm_kpt,ierr)

 ABI_FREE(paw_dmft%buf_psi)

end subroutine compute_wannier
!!***

!!****f* m_datafordmft/print_wannier
!! NAME
!! print_wannier
!!
!! FUNCTION
!!  Write projected Wannier functions on file.
!!
!! INPUTS
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!  istep = iteration step
!!
!! OUTPUT
!!
!! NOTES
!!
!! SOURCE

subroutine print_wannier(paw_dmft,istep)

 use m_io_tools, only : get_unit

!Arguments ------------------------------------
 integer, intent(in) :: istep
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables ------------------------------
 integer :: iatom,iflavor,ir,itypat,lpawu,nflavor,unt
 character(len=3) :: tag_iter
 character(len=10) :: tag_at
 character(len=500) :: message
 character(len=fnlen) :: tmpfil
!************************************************************************

 if (istep < 10) then
   write(tag_iter,'("00",i1)') istep
 else if (istep >= 10 .and. istep < 100) then
   write(tag_iter,'("0",i2)') istep
 else if (istep >= 100 .and. istep < 1000) then
   write(tag_iter,'(i3)') istep
 else
   tag_iter="xxx"
 end if ! istep

 do iatom=1,paw_dmft%natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   itypat = paw_dmft%typat(iatom)
   nflavor = (2*lpawu+1) * paw_dmft%nspinor * paw_dmft%nsppol
   call int2char4(iatom,tag_at)
   ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')
   tmpfil = trim(paw_dmft%filapp)//'Wannier_functions_iatom_'//trim(tag_at)//'_'//tag_iter
   if (open_file(tmpfil,message,newunit=unt) /= 0) ABI_ERROR(message)
   do ir=1,paw_dmft%radgrid(itypat)%mesh_size
     write(unt,*) paw_dmft%radgrid(itypat)%rad(ir),(dble(paw_dmft%wannier(ir,iflavor,iatom)),iflavor=1,nflavor)
   end do ! ir
   close(unt)
 end do ! iatom

end subroutine print_wannier
!!***

END MODULE m_datafordmft
!!***
