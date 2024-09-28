!!****m* ABINIT/m_datafordmft
!! NAME
!!  m_datafordmft
!!
!! FUNCTION
!! This module produces inputs for the DMFT calculation
!!
!! COPYRIGHT
!! Copyright (C) 2006-2024 ABINIT group (BAmadon)
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
!!  Compute psichi (and print some data for check)
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!        -gprimd(3,3)=dimensional reciprocal space primitive translations
!!        -indsym(4,nsym,natom)=indirect indexing array for atom labels
!!        -symrec(3,3,nsym)=symmetry operations in reciprocal space
!!        - nsym= number of symetry operations
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  dimcprj(natom) = dimension for cprj
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  fermie= Fermi energy
!!  dft_occup <type(oper_type)> = occupations in the correlated orbitals in DFT
!!  mband=maximum number of bands
!!  mkmem =number of k points treated by this node
!!  mpi_enreg=information about MPI parallelization
!!  nkpt=number of k points.
!!  my_nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol) = occupancies of KS states.
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  unpaw = file number for cprj
!!
!! OUTPUT
!!  paw_dmft%psichi(nsppol,nkpt,mband,my_nspinor,dtset%natom,(2*maxlpawu+1))): projections <Psi|chi>
!!  paw_dmft%eigen(paw_dmft%nsppol,paw_dmft%nkpt,paw_dmft%mband)
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! SOURCE

subroutine datafordmft(cg,cprj,cryst_struc,dft_occup,dimcprj,dtset,eigen,mband_cprj,mpi_enreg,&
  & my_nspinor,occ,paw_dmft,paw_ij,pawtab,usecprj,nbandkss)
  
!Arguments ------------------------------------
!scalars
 integer, intent(in) :: mband_cprj,my_nspinor,usecprj
 integer, optional, intent(in) :: nbandkss
 !real(dp),intent(in) :: fermie
 type(MPI_type), intent(in) :: mpi_enreg
 type(dataset_type), intent(in) :: dtset
 type(oper_type), intent(inout) :: dft_occup !vz_i
 type(crystal_t), intent(in) :: cryst_struc
!arrays
 type(paw_dmft_type), intent(inout) :: paw_dmft
 integer, intent(in) :: dimcprj(paw_dmft%natom)
 real(dp), intent(in) :: occ(paw_dmft%mband*paw_dmft%nkpt*paw_dmft%nsppol)
 real(dp), target, intent(in) :: eigen(paw_dmft%mband*paw_dmft%nkpt*paw_dmft%nsppol) 
 real(dp), contiguous, intent(in) :: cg(:,:)
 type(paw_ij_type), intent(in) :: paw_ij(:)
! type(pawcprj_type) :: cprj(cryst_struc%natom,my_nspinor*mband*mkmem*nsppol)
 type(pawcprj_type), intent(in) :: cprj(paw_dmft%natom,my_nspinor*mband_cprj*dtset%mkmem*paw_dmft%nsppol*usecprj)
 type(pawtab_type), intent(in) :: pawtab(paw_dmft%ntypat)
!Local variables-------------------------------
 integer :: band_index,comm_band,comm_kpt,iatom,ib,iband,ibandc,ibuf_chipsi,ibuf_psi
 integer :: icg,icgb,icprj,idijeff,ierr,ik,ikpt,ilmn,im,im1,iproj,ir,irank,ispinor,ispinor1,isppol,itypat,lmn_size
 integer :: lpawu,lpawu1,maxlpawu,maxmeshsize,maxnproju,mband,mbandc
 integer :: me_band,me_kpt,mkmem,natom,nband_k,nband_k_cprj,nbandf,nbandi,ndim,nkpt
 integer :: nproc_band,nproc_spkpt,nproju,npw,nspinor,nsploop,nsppol,nsppol_mem,opt_renorm,option,paral_kgb,pawprtvol
 integer :: siz_buf,siz_buf_psi,siz_paw,siz_proj,siz_wan,unt
 logical :: t2g,prt_wan,use_full_chipsi,verif,x2my2d
 real(dp) :: ima,re,rint
 complex(dpc) :: tmp
 character(len=500) :: message
 type(oper_type) :: loc_norm_check
 integer, allocatable :: displs(:),recvcounts(:)
 complex(dpc), allocatable :: buf_chipsi(:),buf_chipsi_tot(:),buf_psi(:),cwprj(:,:),psi_tmp(:)
 type(pawcprj_type), allocatable :: cwaveprj(:,:)
 type(matlu_type), allocatable :: matlu_temp(:)
 integer,parameter :: spinor_idxs(2,4) = RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
!************************************************************************

!DBG_ENTER("COLL")
!Fake test to keep fermie as argument. REMOVE IT AS SOON AS POSSIBLE ...
 !if(fermie>huge(zero))chinorm=zero
 
 !facpara=1 !mpi_enreg%nproc
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
   write(message,*) " warning: parallelised version        ",nkpt
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

!Init parallelism
 !spaceComm=mpi_enreg%comm_cell
 !if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_kpt
 !me=mpi_enreg%me_kpt
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
   ABI_ERROR(message)
 end if 
 
 !iorder_cprj=0
 ABI_CHECK(dtset%mkmem/=0,"mkmem==0 not supported anymore!")
!todo_ab: extract cprj from file unpaw in the following..
!call abi_abort('COLL')

!----------------------------------- MPI-------------------------------------

 nbandi = paw_dmft%dmftbandi 
 nbandf = paw_dmft%dmftbandf
 !lprojchi=.false.
 !lprojchi = .true.
 t2g    = (paw_dmft%dmft_t2g == 1)
 x2my2d = (paw_dmft%dmft_x2my2d == 1)
 use_full_chipsi = (paw_dmft%dmft_use_full_chipsi == 1)

 if ((use_full_chipsi) .and. mpi_enreg%nproc_fft > 1) then
   message = "datafordmft not working when nproc_fft > 1 and either use_full_chipsi or prt_wan"
   ABI_ERROR(message)
 end if
 !natom=cryst_struc%natom

!if(mpi_enreg%me==0) write(7886,*) "in datafordmft", mpi_enreg%me, mpi_enreg%nproc
!if(mpi_enreg%me==1) write(7887,*) "in datafordmft", mpi_enreg%me, mpi_enreg%nproc
!if(mpi_enreg%me==2) write(7888,*) "in datafordmft", mpi_enreg%me, mpi_enreg%nproc
 write(message,'(2a)') ch10,'  == Prepare data for DMFT calculation  '
 call wrtout(std_out,message,'COLL')
 if (abs(pawprtvol) >= 3) then
   write(message,'(a,a)') ch10,'---------------------------------------------------------------'
!  call wrtout(ab_out,message,'COLL')
!  call wrtout(std_out,  message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message,'(a,a,a,a,a,a,a,a,a,a,a,a)') ch10,'  Print useful data (as a check)',ch10,&
    & '  - Overlap of KS wfc with atomic orbital inside sphere',ch10,&
    & '  - Eigenvalues',ch10,&
    & '  - Weights of k-points',ch10,&
    & '  - Number of spins ',ch10,&
    & '  - Number of states'
!  call wrtout(ab_out,message,'COLL')
!  call wrtout(std_out,  message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message,'(a,a)') ch10,'---------------------------------------------------------------'
   call wrtout(std_out,message,'COLL')
 end if ! abs(pawprtvol)>=3
 
 if (dtset%nstep == 0 .and. dtset%nbandkss == 0) then
   message = 'nstep should be greater than 1'
   ABI_BUG(message)
 end if

!********************* Max Values for U terms.
!maxlpawu=0
 !maxnproju=0
 !do iatom=1,natom
 !  if(pawtab(dtset%typat(iatom))%lpawu.ne.-1 .and. pawtab(dtset%typat(iatom))%nproju.gt.maxnproju)&
!&   maxnproju=pawtab(dtset%typat(iatom))%nproju
! end do
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

 if(abs(pawprtvol) >= 3) then
   write(message,'(a,a)') ch10,'   datafordmft :  eigenvalues written'
   call wrtout(std_out,message,'COLL')
 end if
!==========================================================================
!***************** Compute  <Psi|Chi>=\sum_{proja} <Psi|P_a><phi_a|Chi>
!==========================================================================
!write(std_out,*) "size(cprj,dim=1)",size(cprj,dim=1),size(cprj,dim=2),dtset%mband,dtset%mkmem,dtset%nkpt

!Allocate temporary cwaveprj storage
 ABI_MALLOC(cwaveprj,(natom,my_nspinor))
!write(std_out,*) "before alloc cprj"
!write(std_out,*) size(cwaveprj,dim=1),size(cwaveprj,dim=2),size(dimcprj,dim=1)

 call pawcprj_alloc(cwaveprj(:,:),0,dimcprj(:))
!write(std_out,*) "after alloc cprj"
 
 ABI_MALLOC(cwprj,(maxnproju,2*maxlpawu+1))
 ABI_MALLOC(psi_tmp,(maxmeshsize))
 
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
 siz_buf = siz_buf * mbandc * nspinor  
 siz_buf = siz_buf * mkmem * nsppol_mem
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
 if (prt_wan) ABI_MALLOC(paw_dmft%buf_psi,(siz_buf_psi))

 !nprocband=(mband/mband_cprj)
 icprj = 0 
 icg = 0
 ibuf_psi = 0
 ibuf_chipsi = 0
 
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
     ib = 0
!    write(2011,*) ikpt
!    ibsp=ibg
     !ibandc = 0
!    LOOP OVER BANDS
     !ib = 0
     
     do iband=1,nband_k

!      Parallelization: treat only some bands
       verif = .true.
       if (paral_kgb == 1) then
         if (mod((iband-1)/mpi_enreg%bandpp,nproc_band) /= me_band) verif = .false.
       else
         if (mpi_enreg%proc_distrb(ikpt,iband,isppol) /= me_kpt) verif = .false.
       end if
       
       if (verif) ib = ib + 1
       
       if (paw_dmft%band_in(iband)) then
         ibandc = ibandc + 1
       else
         icgb = icgb + npw*nspinor
         cycle
       end if       
              
       if (verif) call pawcprj_get(cryst_struc%atindx1(:),cwaveprj(:,:),cprj(:,:),natom,ib,icprj,ikpt,&
                                 & 0,isppol,mband_cprj,dtset%mkmem,natom,1,nband_k_cprj,&
                                 & nspinor,nsppol,paw_dmft%unpaw,mpicomm=mpi_enreg%comm_kpt,&
                                 & proc_distrb=mpi_enreg%proc_distrb(:,:,:))

       do ispinor=1,nspinor
!        ibsp =~ (ikpt-1)*nband*my_nspinor+iband
!        ibsp=ibsp+1
         !icat=1
!        write(std_out,*) isppol,ikpt,iband,ispinor
         !iat=0 ! to declare
         !do itypat=1,dtset%ntypat
          ! lmn_size=pawtab(itypat)%lmn_size
!          write(std_out,*) isppol,ikpt,iband,ispinor
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
          
           if (verif) then
             lmn_size = pawtab(itypat)%lmn_size
             do ilmn=1,lmn_size
               if (pawtab(itypat)%indlmn(1,ilmn) /= lpawu1) cycle
               im = pawtab(itypat)%indlmn(2,ilmn) + lpawu1 + 1
               if (x2my2d) then
                 if (im /= 5) cycle
                 im = 1
               end if 
               if (t2g) then
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
             
               do ir=1,siz_wan
               
                 tmp = sum(paw_dmft%dpro(1:npw,iatom,ik) * paw_dmft%bessel(1:npw,ir,itypat,ik) * &
                         & cmplx(cg(1,icgb+1:icgb+npw),cg(2,icgb+1:icgb+npw),kind=dp) * &
                         & paw_dmft%ylm(1:npw,im,lpawu+1,ik))
               
                 if (ir <= siz_proj) psi_tmp(ir) = tmp * pawtab(itypat)%proj(ir)
                 if (prt_wan) paw_dmft%buf_psi(ibuf_psi+ir) = tmp
               
               end do ! ir
             
               call simp_gen(re,dble(psi_tmp(1:siz_proj)),paw_dmft%radgrid(itypat),r_for_intg=rint)
               call simp_gen(ima,aimag(psi_tmp(1:siz_proj)),paw_dmft%radgrid(itypat),r_for_intg=rint)
             
               buf_chipsi(ibuf_chipsi) = cmplx(re,ima,kind=dp)
             
               if (verif) then 
                 do iproj=1,nproju
                   buf_chipsi(ibuf_chipsi) = buf_chipsi(ibuf_chipsi) + &
                                           & cwprj(iproj,im) * paw_dmft%phimtphi_int(iproj,itypat)
                   if (prt_wan) paw_dmft%buf_psi(ibuf_psi+1:ibuf_psi+siz_paw) = paw_dmft%buf_psi(ibuf_psi+1:ibuf_psi+siz_paw) + &
                     & cwprj(iproj,im) * paw_dmft%phimtphi(1:siz_paw,iproj,itypat)
                 end do ! iproj
               end if ! verif
               
             else
             
               if (verif) then
                 do iproj=1,nproju
                   buf_chipsi(ibuf_chipsi) = buf_chipsi(ibuf_chipsi) + &
                                           & cwprj(iproj,im) * paw_dmft%phi_int(iproj,itypat)
                 end do ! iproj
               else
                 buf_chipsi(ibuf_chipsi) = czero
               end if ! verif
             
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

!          ----------   t2g case
           !if (t2g) m1_t2g = 0
            ! if (lpawu == 2) then
!                lpawu==2 must be chosen in input and thus in
!                pawtab. On the contrary, paw_dmft now has
!                lpawu=1
              ! m1_t2g = 0 ! index for psichi which has a dimension 3
            ! else
            !   write(message,'(a,a,i4,i4,2a)')  ch10,&
            !    & '  Wrong use of dmftqmc_t2g',paw_dmft%dmftqmc_t2g,lpawu,ch10,&
            !    & ' Action: desactivate qmftqmc_t2g or use lpawu=2'
            !     ABI_ERROR(message)
            ! end if ! lpawu=2
          ! end if ! t2g
!          ----------   t2g case
!          ----------   x2my2d case
           !if (x2my2d) m1_x2my2d = 0
            ! if (lpawu==2) then
!                lpawu==2 must be chosen in input and thus in
!                pawtab. On the contrary, paw_dmft now has
!                lpawu=1
            !     m1_x2my2d=0 ! index for psichi which has a dimension 1
            !   else
            !     write(message,'(a,a,i4,i4,2a)')  ch10,&
!&                 '  Wrong use of dmftqmc_x2my2d',paw_dmft%dmftqmc_x2my2d,lpawu,ch10,&
!&                 ' Action: desactivate dmftqmc_x2my2d or use lpawu=2'
           !      ABI_ERROR(message)
           !    end if
           !  end if
!            ----------   x2my2d case
!            if(isppol==2) write(std_out,*) "ee",size(cprj(iatom,ibsp)%cp(:,:))

            ! iat=iat+1
            ! jj1=0
            ! if(lpawu.ne.-1) then

               !call pawcprj_get(cryst_struc%atindx1,cwaveprj,cprj,natom,ib,ibg,ikpt,&
!&               iorder_cprj,isppol,mband_cprj,mkmem,natom,1,nband_k_cprj,&
!&               my_nspinor,nsppol,unpaw,mpicomm=mpi_enreg%comm_kpt,&
!&               proc_distrb=mpi_enreg%proc_distrb)
!              write(std_out,*) "M-2",isppol,ikpt,iband,iatom,&
!              &             (cwaveprj(iatom,ispinor)%cp(1,13)**2+cwaveprj(iatom,ispinor)%cp(2,13)**2) ,ibsp

!              chinorm=(pawtab(itypat)%phiphjint(1))
               !chinorm=1.d0
!              write(std_out,*) isppol,ikpt,iband,ispinor,iat
              ! do ilmn=1,lmn_size
!                write(std_out,*) ilmn
!                ------------ Select l=lpawu.
              !   if (psps%indlmn(1,ilmn,itypat)==lpawu) then
!                  ------------ Check that the band is choosen (within nbandi and nbandf)
               !    if(paw_dmft%band_in(iband)) then
!                    if(ilmn==13) then
!                    write(std_out,*) "M-2c",isppol,ikpt,iband,iatom
!                    write(std_out,*) "M-2b",isppol,ikpt,ibandc,&
!                    &             (cwaveprj(iatom,ispinor)%cp(1,13)**2+cwaveprj(iatom,ispinor)%cp(2,13)**2)
!                    endif
!                    write(std_out,*) "inside paw_dmft%band_in",iband
                 !    jj1=jj1+1
                 !    if(jj1>pawtab(itypat)%nproju*(2*lpawu+1)) then
                  !     write(message,'(a,a,a,a)')  ch10,&
!&                 !      ' jj1 is not correct in datafordmft',ch10,&
!&                       ' Action: CONTACT Abinit group'
                  !     ABI_BUG(message)
                   !  end if ! jj1
                  !   icount_proj_ilmn=psps%indlmn(3,ilmn,itypat)  ! n: nb of projector
!                    write(std_out,*) "icount_proj_ilmn",icount_proj_ilmn
                 !    m1=psps%indlmn(2,ilmn,itypat)+psps%indlmn(1,ilmn,itypat)+1
!                    ---- if lprochi=true do the sum over every projectors
!                    ---- if lprochi=false: . do the sum over only ground state projector
!                    ----                   . and take ph0phiint(icount_proj_ilmn)=1

!                    call pawcprj_get(cryst_struc%atindx1,cwaveprj,cprj,natom,&
!                    &             ib,ibg,ikpt,iorder_cprj,isppol,mband_cprj,mkmem,&
!                    &             natom,1,nband_k_cprj,my_nspinor,nsppol,unpaw,&
!                    &             mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)

                  !   if(lprojchi.or.icount_proj_ilmn==1) then
                   !    if(.not.lprojchi) ph0phiint_used=one
                    !   if(lprojchi) ph0phiint_used=pawtab(itypat)%ph0phiint(icount_proj_ilmn)
                     !  if(paw_dmft%dmftqmc_t2g==1) then ! t2g case

                      !   if(m1==1.or.m1==2.or.m1==4) then ! t2g case
                      !     m1_t2g=m1_t2g+1  ! t2g case1
                      !     m1_t2g_mod=mod(m1_t2g-1,3)+1
!                          write(std_out,*) "M0",isppol,ikpt,iband,ilmn,cprj(iatom,ibsp)%cp(1,ilmn)
!                          write(std_out,*) "M1",m1,m1_t2g,iband,ilmn,icount_proj_ilmn
                       !    paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_t2g_mod)=&  ! t2g case
!&                      !    paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_t2g_mod)+&  ! t2g case
!                          &                           cmplx(cprj(iatom,ibsp)%cp(1,ilmn)*ph0phiint_used,&  ! t2g case
!                          &                           cprj(iatom,ibsp)%cp(2,ilmn)*ph0phiint_used,kind=dp)  ! t2g case
!&                          cmplx(cwaveprj(iatom,ispinor)%cp(1,ilmn)*ph0phiint_used,&  ! t2g case
!&                          cwaveprj(iatom,ispinor)%cp(2,ilmn)*ph0phiint_used,kind=dp)  ! t2g case
                   !      end if  !t2g case
                    !   else if(paw_dmft%dmftqmc_x2my2d==1) then ! x2my2d case
                    !     if(m1==5) then ! x2my2d case
                    !       m1_x2my2d=1  ! x2my2d case1
                   !        paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_x2my2d)=&      ! x2my2d case
!&                          paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_x2my2d)+&      ! x2my2d case
!&                          cmplx(cwaveprj(iatom,ispinor)%cp(1,ilmn)*ph0phiint_used,&       ! x2my2d case
!&                          cwaveprj(iatom,ispinor)%cp(2,ilmn)*ph0phiint_used,kind=dp)      ! x2my2d case
              !           end if  !x2my2d case
               !        else
      !                   paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)=&
!&                         paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)+&
!                        &                         cmplx(cprj(iatom,ibsp)%cp(1,ilmn)*ph0phiint_used,&
!                        &                         cprj(iatom,ibsp)%cp(2,ilmn)*ph0phiint_used,kind=dp)
!&                         cmplx(cwaveprj(iatom,ispinor)%cp(1,ilmn)*ph0phiint_used,&
!&                         cwaveprj(iatom,ispinor)%cp(2,ilmn)*ph0phiint_used,kind=dp)

!                        if(ibandc==3.and.iat==1.and.m1==1) then
!                        write(std_out,'(a,3i5)') "psichi integers",isppol,ikpt,ispinor
!                        write(std_out,'(a,2i5,2e16.7)') "psichi IB3 iAT1 IM1",ilmn,icount_proj_ilmn,&
!                        &             real(paw_dmft%psichi(isppol,ikpt,3,ispinor,1,1)), imag(paw_dmft%psichi(isppol,ikpt,3,ispinor,1,1))
!                        write(std_out,'(a,2i5,2e16.7)') "cwaveprj IB3 iAT1 IM1",ilmn,icount_proj_ilmn,cwaveprj(iatom,ispinor)%cp(1,ilmn) &
!                        &                    , cwaveprj(iatom,ispinor)%cp(2,ilmn)
!                        endif
              !         end if
              !       end if ! lprojchi=.true. (always)
              !     end if ! nband belongs to paw_dmft%band_in
              !   end if ! L=lpawu
              ! end do !ilmn : loop over L,M,N
   !          end if ! If lpawu/=-1
   !        end do ! iatom
    !       icat=icat+cryst_struc%nattyp(itypat)
    !     end do ! itypat
    !   end do ! ispinor
    ! end do !iband
    ! ibg=ibg+nband_k_cprj*my_nspinor
!    bdtot_index=bdtot_index+nband_k ! useless  (only for occ)
 !  end do !ikpt
 !end do ! isppol
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
 ABI_FREE(psi_tmp)

!==========================================================================
!********************* Gather information for MPI before printing
!==========================================================================

 !dimpsichi=2*nsppol*nkpt*mband*my_nspinor*natom*(2*paw_dmft%maxlpawu+1)
 !ABI_MALLOC(buffer1,(dimpsichi))
 !buffer1 = zero
 !nnn=0
!write(176,*) "beg",psichi
 !do isppol=1,nsppol
 !  do ikpt=1,nkpt
 !    do ibandc=1,paw_dmft%mbandc
 !      do ispinor=1,my_nspinor
 !        do iat=1,natom
 !          do m1=1,2*paw_dmft%maxlpawu+1
!            do m=1,2
 !            nnn=nnn+1
 !            buffer1(nnn)=paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)
!            enddo
 !          end do
 !        end do
 !      end do
 !    end do
 !  end do
 !end do
 call xmpi_barrier(comm_kpt)
 !call xmpi_sum(buffer1,spaceComm ,ierr)
 if (paral_kgb == 1 .and. nproc_band > 1) call xmpi_sum(buf_chipsi(:),comm_band,ierr) !Build sum over band processors
 call xmpi_allgatherv(buf_chipsi(:),siz_buf,buf_chipsi_tot(:),recvcounts(:),displs(:),comm_kpt,ierr)
 
 ABI_FREE(displs)
 ABI_FREE(recvcounts)
 ABI_FREE(buf_chipsi)
 
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
             do im=1,ndim
               ibuf_chipsi = ibuf_chipsi + 1
               paw_dmft%chipsi(im+(ispinor-1)*ndim,ibandc,ikpt,isppol,iatom) = buf_chipsi_tot(ibuf_chipsi)
             end do ! im
           end do ! iatom
         end do ! ispinor
       end do ! ibandc
       
     end do ! ikpt
   end do ! isppol
 end do ! irank
 
 call xmpi_barrier(comm_kpt)
 
 ABI_FREE(buf_chipsi_tot)
 
 !nnn=0
 !do isppol=1,nsppol
 !  do ikpt=1,nkpt
 !    do ibandc=1,paw_dmft%mbandc
 !      do ispinor=1,my_nspinor
 !        do iat=1,natom
 !          do m1=1,2*paw_dmft%maxlpawu+1
!            do m=1,2
 !            nnn=nnn+1
 !            paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)=buffer1(nnn)
!            enddo
             ! if(ibandc==1) then
             !   write(6,*) "m1",m1
             !   write(6,*) "psichi datafordmft",paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)
 !            ! endif
 !          end do
 !        end do
 !      end do
 !    end do
!   end do
! end do
! ABI_FREE(buffer1)
!write(177,*) "end",psichi

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

 !call xmpi_barrier(spaceComm)
!if(mpi_enreg%me.eq.0) write(177,*) "end",psichi
!if(mpi_enreg%me.eq.1) write(178,*) "end",psichi
!if(mpi_enreg%me.eq.2) write(179,*) "end",psichi
 !call xmpi_barrier(spaceComm)
!==========================================================================
!********* WRITE psichi in file for reference
!==========================================================================
 !if (me_kpt == 0) call chipsi_print(paw_dmft,pawtab(:),t2g,x2my2d)
 !end if ! proc=me

!********************* Check normalization  and occupations ***************
! Only if the complete BZ is sampled (ie paw_dmft%kspectralfunc=0)
!==========================================================================
 if (paw_dmft%dmft_kspectralfunc == 0) then
 
   !ABI_MALLOC(loc_occ_check,(natom))
   !ABI_MALLOC(loc_norm_check,(natom))
   !call init_matlu(natom,my_nspinor,nsppol,paw_dmft%lpawu(:),loc_occ_check(:))
   !call init_matlu(natom,my_nspinor,nsppol,paw_dmft%lpawu(:),loc_norm_check(:))
   call init_oper(paw_dmft,loc_norm_check,opt_ksloc=2)
   call chipsi_check(paw_dmft,dft_occup,loc_norm_check)
!==========================================================================
!***************  write checks  *******************************************
!==========================================================================
   if (abs(pawprtvol) >= 3) then
     write(message,*) "normalization computed"
     call wrtout(std_out,message,'COLL')
   end if

   !ABI_MALLOC(loc_occ_check,(natom))
   !ABI_MALLOC(loc_norm_check,(natom))
   !call init_matlu(natom,my_nspinor,nsppol,paw_dmft%lpawu(:),loc_occ_check(:))
   !call init_matlu(natom,my_nspinor,nsppol,paw_dmft%lpawu(:),loc_norm_check(:))
   !call copy_matlu(xocc_check(:),loc_occ_check(:),natom)
   !call copy_matlu(xnorm_check(:),loc_norm_check(:),natom)

   write(message,'(2a,i4)') ch10," == Check: Occupations and Norm from chipsi are"
   call wrtout(std_out,message,'COLL')

   if (paw_dmft%dmftcheck >= 1) then
!    print occupations
     write(message,'(2a,i4)') ch10,'  ------ Unsymmetrized Occupations'
     call wrtout(std_out,message,'COLL')

     call print_matlu(dft_occup%matlu(:),natom,pawprtvol)

!    print norms
     write(message,'(2a,i4)') ch10,'  ------ Unsymmetrized Norm'
     call wrtout(std_out,message,'COLL')

     call print_matlu(loc_norm_check%matlu(:),natom,pawprtvol)
   end if ! dmftcheck>=1

!  symetrise and print occupations
   call sym_matlu(dft_occup%matlu(:),paw_dmft)

   write(message,'(2a,i4)') ch10,'  ------ Symmetrized Occupations'
   call wrtout(std_out,message,'COLL')

   call print_matlu(dft_occup%matlu(:),natom,pawprtvol)

!  symetrise and print norms
   call sym_matlu(loc_norm_check%matlu(:),paw_dmft)

   write(message,'(2a,i4)') ch10,'  ------ Symmetrized Norm'
   call wrtout(std_out,message,'COLL')

   call print_matlu(loc_norm_check%matlu(:),natom,pawprtvol)

!  deallocations
   !do iatom=1,natom
   !  dft_occup%matlu(iatom)%mat=loc_occ_check(iatom)%mat
   !end do

!  Tests density matrix DFT+U and density matrix computed here.
   if (paw_dmft%dmftcheck == 2 .or. (paw_dmft%dmftbandi == 1)) then
     ABI_MALLOC(matlu_temp,(natom))
     call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu_temp(:))
     isppol   = 1
     ispinor  = 1
     ispinor1 = 1
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       !if(paw_dmft%lpawu(iatom).ne.-1) then
       ndim = 2*lpawu + 1
       nsploop = max(nsppol,nspinor**2)
       do idijeff=1,nsploop
         !ispinor = 0
         !ispinor1 = 0
         if (nsploop <= 2) then
           isppol = idijeff ! spinor_idxs(1,idijeff)
           !ispinor = 1
           !ispinor1 = 1
         else if (nsploop == 4) then
           !isppol = 1
           ispinor  = spinor_idxs(1,idijeff)
           ispinor1 = spinor_idxs(2,idijeff)
         !else if (nsploop == 1) then
         !  isppol = 1
         !  ispinor = 1
         !  ispinor1 = 1
         else
           write(message,'(2a)') " BUG in datafordmft: nsploop should be equal to 2 or 4"
           call wrtout(std_out,message,'COLL')
         end if ! nsploop
         do im1=1,ndim
           do im=1,ndim
             if (my_nspinor == 2) matlu_temp(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol) = &
               & cmplx(paw_ij(iatom)%noccmmp(1,im,im1,idijeff),paw_ij(iatom)%noccmmp(2,im,im1,idijeff),kind=dp)
             if (my_nspinor == 1) matlu_temp(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol) = &
               & cmplx(paw_ij(iatom)%noccmmp(1,im,im1,idijeff),zero,kind=dp)
           end do ! im1
         end do ! im
       end do ! idijeff
     end do ! iatom
     if (paw_dmft%dmftcheck == 2) option = 1
     if (paw_dmft%dmftcheck <= 1) option = 0
     call diff_matlu("DFT+U density matrix from INPUT wfk",&
       & "Direct calculation of density matrix with chipsi from DIAGONALIZED wfk",&
       & matlu_temp(:),dft_occup%matlu(:),natom,option,tol3,ierr) !tol1 tol2 tol3
     if (ierr == -1) then
       write(message,'(10a)') ch10,&
        & '    -> These two quantities should agree if three conditions are fulfilled',ch10,&
        & '         -  input wavefunctions come from the same Hamiltonien (e.g LDA/GGA)',ch10,&
        & '         -  dmatpuopt is equal to 1',ch10,&
        & '         -  all valence states are in the valence',ch10,&
        & '    (for experts users: it is not compulsary that these conditions are fulfilled)'
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
   !call destroy_matlu(loc_norm_check(:),natom)
   !ABI_FREE(loc_norm_check)
   !call destroy_matlu(loc_occ_check(:),natom)
   !ABI_FREE(loc_occ_check)

   !call destroy_matlu(xocc_check,natom)
   !call destroy_matlu(xnorm_check,natom)
   !ABI_FREE(xocc_check)
   !ABI_FREE(xnorm_check)
 end if ! dmft_kspectralfunc=0

 if (present(nbandkss)) then
   if ((me_kpt == 0 .and. nbandkss /= 0) .or. (paw_dmft%dmft_kspectralfunc == 1)) then
!     opt_renorm=1 ! if ucrpa==1, no need for individual orthonormalization
     opt_renorm = 3
     if (dtset%ucrpa >= 1 .or. paw_dmft%dmft_kspectralfunc == 1) opt_renorm = 2
     call chipsi_renormalization(paw_dmft,opt=opt_renorm)
     call chipsi_print(paw_dmft,pawtab(:))
   end if ! proc=me
 end if
!!***

 CONTAINS

!!****f* m_datafordmft/psichi_check
!! NAME
!!  psichi_check
!!
!! FUNCTION
!!  Check psichi: compute occupations
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  maxnproju = maximum number of projector for DFT+U species
!!  nattyp(ntypat)= # atoms of each type
!!  mband= number of bands
!!  nkpt=number of k points
!!  my_nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat= number of species
!!  paw_dmft <type(paw_dmft)>=paw data for the self-consistency
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!!  OUTPUTS:
!!  xocc_check(nsppol,my_nspinor,my_nspinor,natom,2*maxlpawu+1,2*maxlpawu+1): density matrix
!!  xnorm_check(nsppol,my_nspinor,my_nspinor,natom,2*maxlpawu+1,2*maxlpawu+1): matrix of norms
!!
!! SIDE EFFECTS
!!  check psichi: compute norm and occupations
!!
!! SOURCE

subroutine chipsi_check(paw_dmft,xocc_check,xnorm_check)

!Arguments ------------------------------------
!scalars
 !integer,intent(in) :: nkpt,my_nspinor,nsppol,ntypat
!arrays
 !integer, intent(in) :: nattyp(ntypat)
 !type(dataset_type),intent(in) :: dtset
 !type(pseudopotential_type),intent(in) :: psps
 !integer, intent(in) :: my_nspinor
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(oper_type), intent(inout) :: xnorm_check,xocc_check !vz_i
 !type(pawtab_type), intent(in) :: pawtab(paw_dmft%ntypat)
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
 
  ! facnsppol=1
  ! if(my_nspinor==1.and.nsppol==1) then
  !   facnsppol=2
  ! end if

  ! ibg=0
  ! band_index=0
  ! do isppol=1,nsppol
   !  do ikpt=1,nkpt
   !    nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
   !    ibandc=0
   !    do iband=1,nband_k
   !      if(paw_dmft%band_in(iband)) ibandc=ibandc+1
   !      do ispinor=1,my_nspinor
   !        icat=1
   !        iat=0
   !        do itypat=1,dtset%ntypat
   !          lmn_size=pawtab(itypat)%lmn_size
   !          do iatom=icat,icat+nattyp(itypat)-1
   !            iat=iat+1
!              ------------ Select correlated atoms
    !           if(paw_dmft%lpawu(iatom).ne.-1) then
    !             chinorm=1.d0
    !             do ilmn=1,lmn_size
!                  ------------ Select l=lpawu.
    !               if (psps%indlmn(1,ilmn,itypat)==paw_dmft%lpawu(iatom).and.&
!&                   psps%indlmn(3,ilmn,itypat)==1) then
  !                   do ilmn1=1,lmn_size
!                      ------------ Select l=lpawu and do not sum over projectors
!                      (this is already done in paw_dmft%psichi)
   !                    if (psps%indlmn(1,ilmn1,itypat)==paw_dmft%lpawu(iatom).and.&
!&                       psps%indlmn(3,ilmn1,itypat)==1) then
!                        ------------ Check that the band is choosen (within nbandi and nbandf)
    !                     if(paw_dmft%band_in(iband)) then
    !                       m=psps%indlmn(2,ilmn,itypat)+psps%indlmn(1,ilmn,itypat)+1
    !                       m1=psps%indlmn(2,ilmn1,itypat)+psps%indlmn(1,ilmn,itypat)+1
    !                       if(psps%indlmn(3,ilmn,itypat)==1) then
    !                         do ispinor1=1,my_nspinor
    !                           psichic=paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m)
    !                           psichic1=paw_dmft%psichi(isppol,ikpt,ibandc,ispinor1,iat,m1)
!                              ------------ Compute occupation matrix
    !                           xocc_check(iatom)%mat(m,m1,isppol,ispinor,ispinor1)=&
!&                               xocc_check(iatom)%mat(m,m1,isppol,ispinor,ispinor1)&
!                              &               +conjg(psichic)*psichic1*dtset%wtk(ikpt)*facpara*occ(iband+band_index)
!&                               +conjg(psichic1)*psichic*dtset%wtk(ikpt)*facpara*occ(iband+band_index)/facnsppol
!                              ------------ Compute norm (should be equal to noccmmp
!                              (dmatpuopt=1) if all bands are taken into account)
 !                              xnorm_check(iatom)%mat(m,m1,isppol,ispinor,ispinor1)=&
!&                               xnorm_check(iatom)%mat(m,m1,isppol,ispinor,ispinor1)&
!                              &               +conjg(psichic)*psichic1*dtset%wtk(ikpt)*facpara
!&                               +conjg(psichic1)*psichic*dtset%wtk(ikpt)*facpara
 !                            end do ! ispinor1
 !                          end if
 !                        end if ! paw_dmft%band_in
 !                      end if
 !                    end do !ilmn1
 !                  end if
 !                end do !ilmn
 !              end if ! lpawu.ne.-1
 !            end do ! iatom
     !        icat=icat+nattyp(itypat)
 !          end do ! itypat
 !        end do ! ispinor
 !      end do !iband
 !      band_index=band_index+nband_k
 !      ibg=ibg+nband_k*my_nspinor
 !    end do !ikpt
!   end do ! isppol

 end subroutine chipsi_check
!DBG_EXIT("COLL")
end subroutine datafordmft
!!***

!!****f* m_datafordmft/psichi_print
!! NAME
!!  psichi_print
!!
!! FUNCTION
!!  Print psichi for reference
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  maxnproju = maximum number of projector for DFT+U species
!!  nattyp(ntypat)= # atoms of each type
!!  mband= number of bands
!!  nkpt=number of k points
!!  my_nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  paw_dmft <type(paw_dmft)>=paw data for the self-consistency
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psichi(2,nsppol,nkpt,mband,my_nspinor,dtset%natom,(2*paw_dmft%maxlpawu+1))) projections for DMFT
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!!
!! SIDE EFFECTS
!!  print psichi in forlb.ovlp
!!
!! SOURCE

subroutine chipsi_print(paw_dmft,pawtab)

!Arguments ------------------------------------
!scalars
 !integer,intent(in) :: nkpt,my_nspinor,nsppol,ntypat
!arrays
 !integer, intent(in) :: nattyp(ntypat)
 !type(dataset_type),intent(in) :: dtset
 !type(pseudopotential_type),intent(in) :: psps
 !logical, intent(in) :: t2g,x2my2d
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(pawtab_type), intent(in) :: pawtab(paw_dmft%ntypat)
!Local variables ------------------------------------
 integer :: iatom,iband,ibandc,ikpt,im,im1,ispinor,isppol
 integer :: itypat,lpawu,nband_k,ndim,unt
 logical :: t2g,x2my2d
 character(len=500) :: msg
 integer, parameter :: mt2g(3) = (/1,2,4/)
! *********************************************************************
   !ll = 1
   t2g = paw_dmft%dmft_t2g == 1
   x2my2d = paw_dmft%dmft_x2my2d == 1

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

   !ibg = 0
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
             !itypat = paw_dmft%typat(iatom)
             !lpawu1 = lpawu
             !if (t2g .or. x2my2d) lpawu1 = 2
           !icat = 1
!          write(std_out,*) isppol,ikpt,iband,ispinor
           !iat = 0 ! to declare
           !do itypat=1,dtset%ntypat
             !lmn_size = pawtab(itypat)%lmn_size
!            write(std_out,*) isppol,ikpt,iband,ispinor
             !do iatom=icat,icat+nattyp(itypat)-1
              ! iat = iat + 1
               !jj1 = 0
               !if(pawtab(itypat)%lpawu.ne.-1) then
!                chinorm=(pawtab(itypat)%phiphjint(1))
               !  chinorm=1.d0
!                write(std_out,*) isppol,ikpt,iband,ispinor,iat
                ! m1_t2g=0
                ! m1_x2my2d=0
             do im=1,ndim
!                  write(std_out,*) ilmn
!                  ------------ Select l=lpawu.  ---------------------------------------
               !if (pawtab(itypat)%indlmn(1,ilmn) /= lpawu1) cycle
                  ! if (psps%indlmn(1,ilmn,itypat)==pawtab(itypat)%lpawu.and.psps%indlmn(3,ilmn,itypat)==1) then
!                    ------------ Check that the band is choosen (within nbandi and nbandf)
                   !  if(paw_dmft%band_in(iband)) then
                    !   jj1=jj1+1
                     !  if(jj1>(2*pawtab(itypat)%lpawu+1)) then
                      !   write(message,'(a,a,i4,i5,i4)') ch10," Error 2 in datafordmft",jj1,pawtab(itypat)%lpawu
                       !  call wrtout(std_out,  message,'COLL')
                      !   ABI_ERROR("Aborting now")
                      ! end if ! jj1
!                      if(jj1>pawtab(dtset%typat(iatom))%nproju*(2*lpawu+1)) then
!                      write(message,'(a,a,a,a)')  ch10,&
!                      &                         'BUG: jj1 is not correct in datafordmft psichi_print',ch10,&
!                      &                         'Action: CONTACT Abinit group'
!                      call wrtout(std_out,  message,'COLL')
!                      call abi_abort('COLL')
!                      end if ! jj1
               !m1 = pawtab(itypat)%indlmn(2,ilmn) + pawtab(itypat)%indlmn(1,ilmn) + 1
!                      ----- Print only when the sum over projectors is done
!                      write(std_out,*) ilmn,m1
               if (t2g) then
                 im1 = mt2g(im)
               else if (x2my2d) then
                 im1 = 5
               else
                 im1 = im
               end if
                         !m1_t2g=m1_t2g+1
               write(unt,'(4i6,3x,2f23.15)') isppol,ispinor,iatom,im1,&
                        & dble(paw_dmft%chipsi(im+(ispinor-1)*ndim,ibandc,ikpt,isppol,iatom)),&
                        & aimag(paw_dmft%chipsi(im+(ispinor-1)*ndim,ibandc,ikpt,isppol,iatom))
                      ! else if(x2my2d) then
                      !   if(m1==5) then
                      !     m1_x2my2d=1
                      !     write(unt,'(4i6,3x,2f23.15)') isppol, ispinor, iat, m1,&
!&                           real(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_x2my2d))/chinorm,&
!&                           aimag(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_x2my2d))/chinorm
               !          end if
                !       else
                !         write(unt,'(4i6,3x,2f23.15)') isppol, ispinor, iat, m1,&
!&                         real(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1))/chinorm,&
!&                         aimag(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1))/chinorm
              !         end if
!                      if natom=1 then jj1 maximal value should be 2*lpawu+1
               !      end if ! paw_dmft%band_in
                 !  end if
             end do !im
             !  end if ! lpawu.ne.-1
           end do ! iatom
            ! icat=icat+nattyp(itypat)
           !end do ! itypat
         end do ! ispinor
       end do !iband
       !ibg=ibg+nband_k*my_nspinor
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
!! Compute levels for ctqmc
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! SOURCE

 subroutine compute_levels(energy_level,hdc,paw_dmft,nondiag,loc_levels)

!Arguments ------------------------------------
!scalars
 !type(crystal_t),intent(in) :: cryst_struc
 !type(pawang_type), intent(in) :: pawang
 type(oper_type), intent(in) :: hdc
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(oper_type), intent(inout) :: energy_level !vz_i
 logical, optional, intent(out) :: nondiag
 type(oper_type), optional, intent(in) :: loc_levels
!Local variables ------------------------------
! scalars
 integer :: iatom,im,isppol,lpawu,natom,ndim
 character(len=500) :: message
! arrays
!************************************************************************

 !mbandc=paw_dmft%mbandc
 !nkpt=paw_dmft%nkpt
 !nsppol=paw_dmft%nsppol
 !nspinor=paw_dmft%nspinor
 natom = paw_dmft%natom
 if (present(nondiag)) nondiag = .false.
 
!========================
!Get KS eigenvalues
!========================
 !do iband=1,mbandc
 !  do ikpt=1,nkpt
 !    do isppol=1,nsppol
!      Take \epsilon_{nks}
!      ========================
  !     energy_level%ks(isppol,ikpt,iband,iband)=paw_dmft%eigen_dft(isppol,ikpt,iband)
 !    end do
  ! end do
 !end do


!======================================================================
!Compute atomic levels from projection of \epsilon_{nks} and symetrize
!======================================================================
 !call loc_oper(energy_level,paw_dmft,1)
 if (present(loc_levels)) then
   call copy_matlu(loc_levels%matlu(:),energy_level%matlu(:),natom)
 else
   call downfold_oper(energy_level,paw_dmft,option=3,op_ks_diag=paw_dmft%eigen_dft(:,:,:))
 end if 
! write(message,'(a,2x,a,f13.5)') ch10," == Print Energy levels before sym and only DFT"
! call wrtout(std_out,message,'COLL')
! call print_matlu(energy_level%matlu,natom,1)
 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   ndim = paw_dmft%nspinor * (2*lpawu+1)
   !if(lpawu/=-1) then
   do isppol=1,paw_dmft%nsppol
     do im=1,ndim
         !do im1=1,2*lpawu+1
       energy_level%matlu(iatom)%mat(im,im,isppol) = energy_level%matlu(iatom)%mat(im,im,isppol) &
          & - paw_dmft%fermie
        ! end do
     end do ! im
   end do ! isppol
   energy_level%matlu(iatom)%mat(:,:,:) = energy_level%matlu(iatom)%mat(:,:,:) - hdc%matlu(iatom)%mat(:,:,:)
!    write(std_out,*) "DC,fermie",hdc%matlu(iatom)%mat(1,1,1,1,1),paw_dmft%fermie
   !end if
 end do ! iatom
 
 call sym_matlu(energy_level%matlu(:),paw_dmft)
 if (present(nondiag)) call checkdiag_matlu(energy_level%matlu(:),natom,tol7,nondiag)

 write(message,'(a,2x,a,f13.5)') ch10," == Print Energy levels for Fermi Level=",paw_dmft%fermie
 call wrtout(std_out,message,'COLL')
!call print_oper(energy_level,1,paw_dmft,1)
 call print_matlu(energy_level%matlu(:),natom,1)

 end subroutine compute_levels
!!***


!!****f* m_datafordmft/psichi_renormalization
!! NAME
!! psichi_renormalization
!!
!! FUNCTION
!! Renormalize psichi.
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>= crystal structure data.
!!  paw_dmft =  data for DFT+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!
!! OUTPUT
!!  paw_dmft%psichi(nsppol,nkpt,mband,nspinor,dtset%natom,(2*maxlpawu+1))): projections <Psi|chi> are orthonormalized.
!!
!! NOTES
!!
!! SOURCE

subroutine chipsi_renormalization(paw_dmft,opt)

!Arguments ------------------------------------
!scalars
 !type(crystal_t),intent(in) :: cryst_struc
 type(paw_dmft_type), intent(inout) :: paw_dmft
 !type(pawang_type), intent(in) :: pawang
 integer, optional, intent(in) :: opt
!Local variables ------------------------------
!scalars
 integer :: jkpt,natom,nkpt,option
 !real(dp), pointer :: temp_wtk(:) => null()
 real(dp) :: pawprtvol
 type(oper_type) :: norm,oper_temp
 character(len=500) :: message
 real(dp), allocatable :: wtk_tmp(:)
!arrays
! real(dp),allocatable :: e0pde(:,:,:),omegame0i(:)
!************************************************************************

 DBG_ENTER("COLL")

 option = 3
 if (present(opt)) then
   if (opt == 2 .or. opt == 3) option = opt
  ! if(opt==3) option=3
 end if
 pawprtvol = 2

 !nsppol  = paw_dmft%nsppol
 natom = paw_dmft%natom
 nkpt  = paw_dmft%nkpt
 !mbandc  = paw_dmft%mbandc
 !nspinor = paw_dmft%nspinor


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

 if (option == 2) then ! option==2
!  ====================================
!  == renormalize k-point after k-point
!  ====================================
   !ABI_MALLOC(temp_wtk,(1))

   write(message,'(2a,i5)') ch10,' Nb of K-point',nkpt
   call wrtout(std_out,message,'COLL')
   do jkpt=1,nkpt  ! jkpt
     write(message,'(2a,i5)') ch10,'  == Renormalization for K-point',jkpt
     call wrtout(std_out,message,'COLL')
    ! temp_wtk(1)=one
     if (paw_dmft%distrib%procb(jkpt) /= paw_dmft%distrib%me_kpt) cycle
     call normalizechipsi(1,paw_dmft,jkpt=jkpt)
   end do ! jkpt
   write(message,'(2a)') ch10,'  ===== K-points all renormalized'
   call wrtout(std_out,message,'COLL')
   !ABI_FREE(temp_wtk)

 else if (option == 3) then  ! option==3
!  ====================================
!  == renormalize the sum over k-points
!  ====================================
!  allocate(temp_wtk(nkpt))
   !temp_wtk=>paw_dmft%wtk
   write(message,'(6a)') ch10,'  ====================================== ',&
     & ch10,'  == Renormalization for all K-points == ',&
     & ch10,'  ======================================='
   call wrtout(std_out,message,'COLL')
   call normalizechipsi(nkpt,paw_dmft)
!  deallocate(temp_wtk)

 end if ! option

 call chipsi_gather(paw_dmft)

!== Change back repr for norm


!===============================================
!==  Compute norm with new psichi
!===============================================
 write(message,'(2a)') ch10,'  ===== Compute norm with new psichi'
 call wrtout(std_out,message,'COLL')
 call init_oper(paw_dmft,oper_temp,opt_ksloc=2)
 call identity_oper(oper_temp,2)

!== Build identity for norm%ks (option=1)
 !call identity_oper(norm,1)
 ABI_MALLOC(wtk_tmp,(nkpt))
!== Deduce norm%matlu from norm%ks with new psichi
 if (paw_dmft%dmft_kspectralfunc == 1) then
   !call identity_oper(oper_temp,2)
   wtk_tmp(:) = one
   call init_oper(paw_dmft,norm,nkpt=1,wtk=wtk_tmp(:),opt_ksloc=2)
   do jkpt=1,nkpt  ! jkpt
     norm%shiftk = jkpt - 1
     call downfold_oper(norm,paw_dmft,option=2)
     write(message,'(2a,i5)') &
       & ch10,"  == Check: Overlap with renormalized psichi for k-point",jkpt
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm%matlu(:),natom,prtopt=1)
!== Check that norm is now the identity
     !call identity_oper(oper_temp,2)
     call diff_matlu('Overlap after renormalization','Identity',&
       & norm%matlu(:),oper_temp%matlu(:),natom,1,tol6,zero_or_one=1)
   end do ! jkpt
   !call destroy_oper(oper_temp)
 else !dmft_kspectralfunc
   write(message,'(2a)') ch10,'  ===== Calling loc_oper'
   call wrtout(std_out,message,'COLL')
   call init_oper(paw_dmft,norm,opt_ksloc=2)
   call downfold_oper(norm,paw_dmft,procb=paw_dmft%distrib%procb(:),iproc=paw_dmft%distrib%me_kpt,option=2)
   call xmpi_matlu(norm%matlu(:),natom,paw_dmft%distrib%comm_kpt)
   write(message,'(2a)') ch10,'  ===== Finished loc_oper'
   call wrtout(std_out,message,'COLL')

!== Print norm%matlu unsymetrized with new psichi
   if (pawprtvol > 2) then
     write(message,'(4a,2a)') &
      &  ch10,"  == Check: Overlap with renormalized psichi without symmetrization is == "
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm%matlu(:),natom,prtopt=1)
   end if ! pawprtvol>2


!== Symetrise norm%matlu with new psichi
   call sym_matlu(norm%matlu(:),paw_dmft)

!== Print norm%matlu symetrized with new psichi
   if (pawprtvol > 2) then
     write(message,'(4a,2a)') &
      & ch10,"  == Check: Overlap with renormalized psichi and symetrization is =="
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm%matlu(:),natom,prtopt=1,opt_diag=-1)
   end if

!== Check that norm is now the identity
  ! call init_oper(paw_dmft,oper_temp)
   !call identity_oper(oper_temp,2)
   call diff_oper('Overlap after renormalization','Identity',norm,oper_temp,0,tol6)

 end if ! dmft_kspectralfunc

 call destroy_oper(norm)
 call destroy_oper(oper_temp)

 ABI_FREE(wtk_tmp)

 paw_dmft%lchipsiortho = 1

 DBG_EXIT("COLL")

 CONTAINS
!===========================================================
!!***

!!****f* psichi_renormalization/normalizepsichi
!! NAME
!!  normalizepsichi
!!
!! FUNCTION
!!  normalizepsichi psichi from norm
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  change psichi: normalizepsichi it
!!
!! SOURCE

subroutine normalizechipsi(nkpt,paw_dmft,jkpt)

!Arguments ------------------------------------
 integer, intent(in) :: nkpt
 integer, optional, intent(in) :: jkpt
 !real(dp),pointer :: temp_wtk(:)
!scalars
 !type(crystal_t),intent(in) :: cryst_struc
 type(paw_dmft_type), intent(inout) :: paw_dmft
 !type(pawang_type), intent(in) :: pawang
!Local variables ------------------------------
 integer :: dimoverlap,dum,iatom,ib,ikpt,im,isppol,itot,itot1,lpawu,mbandc
 integer :: natom,ndim,ndim_max,nspinor,nsppol,pawprtvol
 type(oper_type) :: norm1,norm2,norm3
 complex(dpc), allocatable :: chipsivect(:,:),largeoverlap(:,:),mat_tmp(:,:)
 character(len=500) :: message
!scalars
!arrays
! real(dp),allocatable :: e0pde(:,:,:),omegame0i(:)
!************************************************************************
   !nsppol  = paw_dmft%nsppol
   !mbandc  = paw_dmft%mbandc
   !natom   = cryst_struc%natom
   !nspinor = paw_dmft%nspinor
 mbandc    = paw_dmft%mbandc
 natom     = paw_dmft%natom
 nspinor   = paw_dmft%nspinor
 nsppol    = paw_dmft%nsppol
 ndim_max  = (2*paw_dmft%maxlpawu+1)*nspinor
 pawprtvol = 3
   !diag = 0
   !natomcor=0
  ! dimoverlap=0
  ! do iatom=1,natom
  !   if(paw_dmft%lpawu(iatom).ne.-1) then
  !     natomcor=natomcor+1
  !     ndim=2*paw_dmft%lpawu(iatom)+1
  !     tndim=nspinor*ndim
  !     dimoverlap=dimoverlap+tndim
      ! write(6,*) "atom, dimoverlap",iatom,dimoverlap,natomcor
  !   end if
  ! end do

 if (nkpt /= 1 .and. present(jkpt)) then
   message = 'BUG in chipsi_normalization'
   ABI_ERROR(message)
 end if

  ! iortho=1
  ! write(6,*) "nkpt, iortho",nkpt,iortho
   !if (natomcor>1) iortho=2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  First case: usual case (should in fact be used only for one atom and nkpt=1)
 if (.not. present(jkpt)) then ! .and.iortho==1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    *********************************************************************

   call init_oper(paw_dmft,norm1,opt_ksloc=2)

!    == Build identity for norm1%ks (option=1)
     !call identity_oper(norm1,1)

     !if(nkpt==1.and.present(jkpt)) then
     !  call loc_oper(norm1,paw_dmft,1,jkpt=jkpt)
     !end if
     !if(.not.present(jkpt)) then
   call downfold_oper(norm1,paw_dmft,procb=paw_dmft%distrib%procb(:),iproc=paw_dmft%distrib%me_kpt,option=2)
   call xmpi_matlu(norm1%matlu(:),natom,paw_dmft%distrib%comm_kpt)
     !end if
   if (nkpt > 1) call sym_matlu(norm1%matlu(:),paw_dmft)

   if (pawprtvol > 2) then
     write(message,'(2a)') ch10,'  - Print norm with current psichi '
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm1%matlu(:),natom,prtopt=1,opt_exp=1)
   end if
!    ==-------------------------------------
!    == Start loop on atoms
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
   write(message,'(2a)') ch10,'  - Overlap (before orthonormalization) -'
   call wrtout(std_out,message,'COLL')
   !call gather_matlu(norm1%matlu,overlap,cryst_struc%natom,option=1,prtopt=1)
   !call destroy_oper(norm1)

   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
       !if(paw_dmft%lpawu(iatom).ne.-1) then
     ndim = nspinor * (2*lpawu+1)
     ABI_MALLOC(mat_tmp,(ndim,mbandc))
     !tndim = nsppol*nspinor*ndim
         !ABI_MALLOC(sqrtmatinv,(tndim,tndim))

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
     do isppol=1,nsppol
       call invsqrt_matrix(norm1%matlu(iatom)%mat(:,:,isppol),ndim,dum)
       do ikpt=1,nkpt
         if (paw_dmft%distrib%procb(ikpt) /= paw_dmft%distrib%me_kpt) cycle
         call abi_xgemm("n","n",ndim,mbandc,ndim,cone,norm1%matlu(iatom)%mat(:,:,isppol), &
                  & ndim,paw_dmft%chipsi(:,:,ikpt,isppol,iatom),ndim_max,czero,mat_tmp(:,:),ndim)
         paw_dmft%chipsi(1:ndim,:,ikpt,isppol,iatom) = mat_tmp(:,:)
       end do ! ikpt
     end do ! isppol

     ABI_FREE(mat_tmp)

   end do ! iatom

           !sqrtmatinv=overlap(iatom)%value
         !else
         !  sqrtmatinv(:,:)=czero
         !  do ib=1,tndim
         !    sqrtmatinv(ib,ib)=cone/(sqrt(overlap(iatom)%value(ib,ib)))
         !  end do
         !end if

!        == Apply O^{-0.5} on psichi
       !  ABI_MALLOC(wan,(nsppol,nspinor,ndim))
!        write(std_out,*) mbandc,nsppol,nspinor,ndim
!        write(std_out,*)  paw_dmft%psichi(1,1,1,1,1,1)
       !  do ikpt=1,nkpt
        !   do ib=1,mbandc
         !    if(present(jkpt)) then
         !      ikpt1=jkpt
         !    else
         !      ikpt1=ikpt
         !    end if
         !    jc=0
          !   wan=czero
          !   do isppol=1,nsppol
          !     do ispinor=1,nspinor
          !       do im=1,ndim
!                  write(std_out,*) "psichi", paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)
          !         jc=jc+1
           !        jc1=0
            !       do isppol1=1,nsppol
           !          do ispinor1=1,nspinor
           !            do im1=1,ndim
           !              jc1=jc1+1
           !              wan(isppol,ispinor,im)= wan(isppol,ispinor,im) &
!&                         + paw_dmft%psichi(isppol1,ikpt1,ib,ispinor1,iatom,im1)*sqrtmatinv(jc,jc1)
          !             end do ! ispinor1
           !          end do ! isppol1
           !        end do ! im1
          !       end do ! im
         !      end do ! ispinor
          !   end do !  isppol
         !    do isppol=1,nsppol
          !     do ispinor=1,nspinor
          !       do im=1,ndim
          !         paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)=wan(isppol,ispinor,im)
!                  write(std_out,*) "psichi2", paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)
         !        end do ! ispinor
         !      end do ! isppol
        !     end do ! im
        !   end do ! ib
       !  end do ! ikpt
       !  ABI_FREE(wan)
        ! ABI_FREE(sqrtmatinv)
!        write(std_out,*)  paw_dmft%psichi(1,1,1,1,1,1)

!        ==-------------------------------------
     !  end if ! lpawu.ne.-1
    ! end do ! iatom
!    == End loop on atoms
!    ==-------------------------------------
    ! do iatom=1,natom
     !  if(paw_dmft%lpawu(iatom).ne.-1) then
    !     ABI_FREE(overlap(iatom)%value)
    !   end if
   !  end do
   !  ABI_FREE(overlap)

!    ======================================================================
!    == Check norm with new psichi.
!    ======================================================================

     !call init_oper(paw_dmft,norm1,nkpt=nkpt,wtk=temp_wtk)

     !call identity_oper(norm1,1)

     !if(nkpt==1.and.present(jkpt)) then
     !  call loc_oper(norm1,paw_dmft,1,jkpt=jkpt)
     !end if
     !if(.not.present(jkpt)) then
   call downfold_oper(norm1,paw_dmft,procb=paw_dmft%distrib%procb(:),iproc=paw_dmft%distrib%me_kpt,option=2)
   call xmpi_matlu(norm1%matlu(:),natom,paw_dmft%distrib%comm_kpt)

     !end if

   if (nkpt > 1) call sym_matlu(norm1%matlu(:),paw_dmft)
   !  end if

   if (pawprtvol > 2) then
     write(message,'(2a)') ch10,'  - Print norm with new psichi '
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm1%matlu(:),natom,prtopt=1)
   end if

!    ======================================================================
!    == Check that norm-identity is zero
!    ======================================================================
   call init_oper(paw_dmft,norm2,opt_ksloc=2)
   call init_oper(paw_dmft,norm3,opt_ksloc=2)
   call identity_oper(norm2,2)
   call add_matlu(norm1%matlu(:),norm2%matlu(:),norm3%matlu(:),natom,-1)
   call destroy_oper(norm2)
   if (pawprtvol > 2) then
     write(message,'(2a)') ch10,'  - Print norm with new psichi minus Identity '
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
     !if(paw_dmft%lpawu(iatom).ne.-1) then
       !natomcor=natomcor+1
       !ndim=2*paw_dmft%lpawu(iatom)+1
       !tndim=nspinor*ndim
     dimoverlap = dimoverlap + 2*lpawu + 1
      ! write(6,*) "atom, dimoverlap",iatom,dimoverlap,natomcor
     !end if
   end do ! iatom

   dimoverlap = dimoverlap * nspinor

   ABI_MALLOC(largeoverlap,(dimoverlap,dimoverlap))
   ABI_MALLOC(chipsivect,(dimoverlap,mbandc))
     !ABI_MALLOC(sqrtmatinv,(dimoverlap,dimoverlap))
     !ABI_MALLOC(wanall,(dimoverlap))
   ABI_MALLOC(mat_tmp,(dimoverlap,mbandc))

!    Big loop over isppol
   do isppol=1,nsppol
     do ib=1,mbandc
       itot = 0
       do iatom=1,natom
         lpawu = paw_dmft%lpawu(iatom)
         if (lpawu == -1) cycle
           !if(paw_dmft%lpawu(iatom).ne.-1) then
         ndim = nspinor * (2*lpawu+1)
         do im=1,ndim
               !do ispinor=1,nspinor
           itot = itot + 1
                 !if(itot>dimoverlap) write(std_out,*) "itot>ndim",itot,ndim
                ! write(6,*) "ib,iatom,im,ispinor",ib,iatom,im,ispinor,jkpt
           chipsivect(itot,ib) = paw_dmft%chipsi(im,ib,jkpt,isppol,iatom)
               !enddo ! ispinor
         end do ! im
           !endif
       end do ! iatom
     end do ! ib


!     calculation of overlap
     call abi_xgemm("n","c",dimoverlap,dimoverlap,mbandc,cone,chipsivect(:,:),dimoverlap,&
                  & chipsivect(:,:),dimoverlap,czero,largeoverlap(:,:),dimoverlap)
     !largeoverlap(:,:) = czero
     !  do ib=1,mbandc
     !    do itot=1,dimoverlap
     !      do itot1=1,dimoverlap
     !         largeoverlap(itot,itot1)=largeoverlap(itot,itot1)+ &
!&              psichivect(ib,itot)*conjg(psichivect(ib,itot1))
 !          enddo ! itot1
 !        enddo ! itot
 !      enddo ! ib

!     Math: orthogonalisation of overlap
     write(std_out,*)"jkpt=",jkpt
     do itot=1,dimoverlap
       write(std_out,'(100f7.3)') (largeoverlap(itot,itot1),itot1=1,dimoverlap)
     end do
     call invsqrt_matrix(largeoverlap(:,:),dimoverlap,dum)
       !sqrtmatinv=largeoverlap
     write(std_out,*)"jkpt=",jkpt
     do itot=1,dimoverlap
       write(std_out,'(100f7.3)') (largeoverlap(itot,itot1),itot1=1,dimoverlap)
     end do
     write(std_out,*)"jkpt=",jkpt
     do itot=1,dimoverlap
       write(std_out,'(100e9.3)') (largeoverlap(itot,itot1),itot1=1,dimoverlap)
     end do
     
     call abi_xgemm("n","n",dimoverlap,mbandc,dimoverlap,cone,largeoverlap(:,:),dimoverlap,&
                  & chipsivect(:,:),dimoverlap,czero,mat_tmp(:,:),dimoverlap)
     

      ! do ib=1,mbandc
      !   wanall=czero
      !   do itot=1,dimoverlap
      !     do itot1=1,dimoverlap
      !        wanall(itot)= wanall(itot)+psichivect(ib,itot1)*sqrtmatinv(itot,itot1)
       !    enddo ! itot1
     !      write(std_out,'(3i3,2x,i3,2x,2e15.5,2x,2e15.5)') jkpt,isppol,ib,itot,psichivect(ib,itot),wanall(itot)
       !  enddo ! itot
       !  iatomcor=0
       !  do itot=1,dimoverlap
        !   psichivect(ib,itot)=wanall(itot)
       !  enddo
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
       !enddo ! ib
        
     call abi_xgemm("n","c",dimoverlap,dimoverlap,mbandc,cone,mat_tmp(:,:),dimoverlap,&
                  & mat_tmp(:,:),dimoverlap,czero,largeoverlap(:,:),dimoverlap)
!     calculation of overlap (check)
      ! largeoverlap=czero
      ! do ib=1,mbandc
      !   do itot=1,dimoverlap
      !     do itot1=1,dimoverlap
      !        largeoverlap(itot,itot1)=largeoverlap(itot,itot1)+ &
!&              psichivect(ib,itot)*conjg(psichivect(ib,itot1))
      !     enddo ! itot1
      !   enddo ! itot
      ! enddo ! ib

     write(std_out,*)"jkpt=",jkpt
     do itot=1,dimoverlap
       write(std_out,'(100f7.3)') (largeoverlap(itot,itot1),itot1=1,dimoverlap)
     end do

!      psichivect -> psichi
     do ib=1,mbandc
       itot = 0
       do iatom=1,natom
         lpawu = paw_dmft%lpawu(iatom)
         if (lpawu == -1) cycle
           !if(paw_dmft%lpawu(iatom).ne.-1) then
         ndim = nspinor * (2*lpawu+1)
             !iatomcor=iatomcor+1
         do im=1,ndim
               !do ispinor=1,nspinor
           itot = itot + 1
           paw_dmft%chipsi(im,ib,jkpt,isppol,iatom) = mat_tmp(itot,ib)
               !end do ! ispinor
         end do ! im
           !endif
       end do ! iatom
     end do ! ib

!   End big loop over isppol
   end do !isppol

   ABI_FREE(chipsivect)
     !ABI_FREE(sqrtmatinv)
     !ABI_FREE(wanall)
   ABI_FREE(largeoverlap)
   ABI_FREE(mat_tmp)

 end if ! option

 end subroutine normalizechipsi

 !!****f* m_datafordmft/psichi_renormalization
!! NAME
!! psichi_renormalization
!!
!! FUNCTION
!! Renormalize psichi.
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>= crystal structure data.
!!  paw_dmft =  data for DFT+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!
!! OUTPUT
!!  paw_dmft%psichi(nsppol,nkpt,mband,nspinor,dtset%natom,(2*maxlpawu+1))): projections <Psi|chi> are orthonormalized.
!!
!! NOTES
!!
!! SOURCE

subroutine chipsi_gather(paw_dmft)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(inout) :: paw_dmft
!Local variables ------------------------------
 integer :: iatom,ib,ibuf,ierr,ikpt,im,irank,isppol,lpawu,mbandc,mkmem
 integer :: myproc,natom,ndim,nproc,nspinor,nsppol,shift,siz_buf
 integer, allocatable :: displs(:),recvcounts(:)
 complex(dpc), allocatable :: buffer(:),buffer_tot(:)
!************************************************************************

 myproc  = paw_dmft%myproc
 mkmem   = paw_dmft%distrib%nkpt_mem(myproc)
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

 ABI_MALLOC(buffer,(recvcounts(myproc+1)))
 ABI_MALLOC(buffer_tot,(recvcounts(nproc)+displs(nproc)))

 ibuf = 0
 do ikpt=1,mkmem
   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     ndim = nspinor * (2*lpawu+1)
     do isppol=1,nsppol
       do ib=1,mbandc
         do im=1,ndim
           ibuf = ibuf + 1
           buffer(ibuf) = paw_dmft%chipsi(im,ib,ikpt+shift,isppol,iatom)
         end do ! im
       end do ! ib
     end do ! isppol
   end do ! iatom
 end do ! ikpt

 call xmpi_allgatherv(buffer(:),recvcounts(myproc+1),buffer_tot(:),recvcounts(:),displs(:),paw_dmft%distrib%comm_kpt,ierr)

 ibuf = 0
 do ikpt=1,paw_dmft%nkpt
   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     ndim = nspinor * (2*lpawu+1)
     do isppol=1,nsppol
       do ib=1,mbandc
         do im=1,ndim
           ibuf = ibuf + 1
           paw_dmft%chipsi(im,ib,ikpt,isppol,iatom) = buffer_tot(ibuf)
         end do ! im
       end do ! ib
     end do ! isppol
   end do ! iatom
 end do ! ikpt

 ABI_FREE(recvcounts)
 ABI_FREE(displs)
 ABI_FREE(buffer)
 ABI_FREE(buffer_tot)

 end subroutine chipsi_gather

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
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab <type(pawtab)>
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!
!! NOTES
!!
!! SOURCE

subroutine hybridization_asymptotic_coefficient(cryst_struc,paw_dmft,pawang,hybri_coeff)

 use m_pawang, only : pawang_type

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(pawang_type), intent(in) :: pawang
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

!!****f* m_datafordmft/psichi_renormalization
!! NAME
!! psichi_renormalization
!!
!! FUNCTION
!!
!! INPUTS
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
 integer :: iatom,ibuf_psi,ibandc,ierr,iflavor,ik,ikpt,im,ispinor
 integer :: isppol,itypat,lpawu,natom,nband_k,ndim,nkpt,nspinor,nsppol,siz_wan
!************************************************************************

 natom   = paw_dmft%natom
 nkpt    = paw_dmft%nkpt
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol

 ABI_MALLOC(paw_dmft%wannier,(paw_dmft%maxmeshsize,nspinor*(2*paw_dmft%maxlpawu+1)*nsppol,natom))
 paw_dmft%wannier(:,:,:) = zero
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

 if (mpi_enreg%paral_kgb == 1 .and. mpi_enreg%nproc_band > 1) call xmpi_sum_master(paw_dmft%wannier(:,:,:),0,mpi_enreg%comm_band,ierr)
 if (mpi_enreg%nproc_spkpt > 1) call xmpi_sum_master(paw_dmft%wannier(:,:,:),0,mpi_enreg%comm_kpt,ierr)

 ABI_FREE(paw_dmft%buf_psi)

end subroutine compute_wannier

!!****f* m_datafordmft/print_wannier
!! NAME
!! print_wannier
!!
!! FUNCTION
!!
!! INPUTS
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
   tmpfil = trim(paw_dmft%filapp)//'Wannier_functions_iatom'//trim(tag_at)//'_'//tag_iter
   unt = get_unit()
   open(unit=unt,file=tmpfil,status='unknown',form='formatted')
   rewind(unt)
   do ir=1,paw_dmft%radgrid(itypat)%mesh_size
     write(unt,*) paw_dmft%radgrid(itypat)%rad(ir),(paw_dmft%wannier(ir,iflavor,iatom),iflavor=1,nflavor)
   end do ! ir
   close(unt)
 end do ! iatom

end subroutine print_wannier
!!***

END MODULE m_datafordmft
!!***

