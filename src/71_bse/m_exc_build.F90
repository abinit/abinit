!!****m* ABINIT/m_exc_build
!! NAME
!!  m_exc_build
!!
!! FUNCTION
!!  Build BSE Hamiltonian in e-h reprensentation.
!!
!! COPYRIGHT
!!  Copyright (C) 1992-2009 EXC group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida)
!!  Copyright (C) 2009-2020 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_exc_build

 use defs_basis
 use m_abicore
 use m_bs_defs
 use m_bse_io
 use m_xmpi
 use m_errors
 use m_screen
#if defined HAVE_MPI2
 use mpi
#endif
 use m_hdr

 use defs_datatypes, only : pseudopotential_type
 use m_gwdefs,       only : czero_gw, cone_gw, GW_TOLQ0
 use m_time,         only : cwtime, timab
 use m_io_tools,     only : get_unit, open_file
 use m_hide_blas,    only : xdotc, xgemv
 use m_geometry,     only : normv
 use m_crystal,      only : crystal_t
 use m_gsphere,      only : gsphere_t, gsph_fft_tabs
 use m_vcoul,        only : vcoul_t
 use m_bz_mesh,      only : kmesh_t, get_BZ_item, get_BZ_diff, has_BZ_item, isamek, findqg0
 use m_pawpwij,      only : pawpwff_t, pawpwij_t, pawpwij_init, pawpwij_free, paw_rho_tw_g
 use m_pawang,       only : pawang_type
 use m_pawtab,       only : pawtab_type
 use m_pawcprj,      only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_paw_sym,      only : paw_symcprj_op
 use m_wfd,          only : wfd_t
 use m_oscillators,  only : rho_tw_g, sym_rhotwgq0

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

 private
!!***

 public :: exc_build_ham  !  Calculate and write the excitonic Hamiltonian to file.
!!***

contains
!!***

!!****f* m_exc_build/exc_build_block
!! NAME
!!  exc_build_block
!!
!! FUNCTION
!!  Calculate and write the excitonic Hamiltonian on an external binary file (Fortran file open
!!  in random mode) for subsequent treatment in the Bethe-Salpeter code.
!!
!! INPUTS
!!  BSp<excparam>=The parameters for the Bethe-Salpeter calculation.
!!  Cryst<crystal_t>=Info on the crystalline structure.
!!  Kmesh<kmesh_t>=The list of k-points in the BZ, IBZ and symmetry tables.
!!  Qmesh<kmesh_t>=The list of q-points for epsilon^{-1} and related symmetry tables.
!!  ktabr(nfftot_osc,BSp%nkbz)=The FFT index of $(R^{-1}(r-\tau))$ where R is symmetry needed to obtains
!!    the k-points from the irreducible image.  Used to symmetrize u_Sk where S = \transpose R^{-1}
!!  Gsph_x<gsphere_t>=Info on the G-sphere used to describe wavefunctions and W (the largest one is actually stored).
!!  Gsph_c<gsphere_t>=Info on the G-sphere used to describe the correlation part.
!!  Vcp<vcoul_t>=The Coulomb interaction in reciprocal space. A cutoff can be used
!!  W<screen_t>=Data type gathering info and data for W.
!!  nfftot_osc=Total Number of FFT points used for the oscillator matrix elements.
!!  ngfft_osc(18)=Info on the FFT algorithm used to calculate the oscillator matrix elements.
!!  Psps<Pseudopotential_type>=Variables related to pseudopotentials
!!  Pawtab(Psps%ntypat)<pawtab_type>=PAW tabulated starting data.
!!  Pawang<pawang_type>=PAW angular mesh and related data.
!!  Paw_pwff(Cryst%ntypat*Wfd%usepaw)<pawpwff_t>=Form factor used to calculate the onsite matrix
!!    elements of a plane wave.
!!  Wfd<wfd_t>=Handler for the wavefunctions.
!!    prtvol=Verbosity level.
!!  rhxtwg_q0
!!  is_resonant
!!  fname
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  The excitonic Hamiltonian is saved on an external binary file (see below).
!!
!! NOTES
!!  *) Version for K_V = K_C (q=0), thus KP_V = KP_C
!!  *) No exchange limit: use LDA energies in case.
!!  *) Symmetry of H(-k-k') = H*(k k') not used.
!!  *) Coulomb term can be approssimateed as diagonal in G.
!!  *) Valence bands treated from lomo on.
!!  *) Symmetries of the sub-blocks are used to reduce the number of elements to calculate.
!!
!!            ____________
!!           |_(cv)__(vc)_|
!!   H_exc = |  R      C  |
!!           | -C*    -R* |
!!
!!   where C is symmetric and R is Hermitian provided that the QP energies are real.
!!
!!  For nsppol=1 ==> R = diag-W+2v; C = -W+2v
!!  since the Hamiltonian can be diagonalized in the spin-singlet basis set thanks to
!!  the fact that spin triplet does not contribute to the optical limit of epsilon.
!!
!!  For nsppol=2 ==> R = diag-W+v; C = -W+v
!!  Now the matrix elements depend on the spin of the transitions but only those
!!  transitions in which the spin of the electron and of the hole are equal contribute
!!  to the macroscopic dielectric function. Moreover only the exchange term can connect transitions of different spin.
!!  When nsppol==2 the transitions are ordered using | (cv up) | (cv dwn) | (vc up) | (vc down) |
!!
!!  The resonant block is given by:
!!
!!      |  (v'c' up)       | (v'c' dwn)   |
!!      -----------------------------------           where v_{-+} = v_{+-}^H when the momentum of the photon is neglected.
!!      | [diag-W+v]++     |      v+-     | (vc up)   Note that v_{+-} is not Hermitian due to the presence of different spins.
!!  R = -----------------------------------           Actually it reduces to a Hermitian matrix when the system is not spin polarized.
!!      |     v-+          | [diag-W+v]-- | (vc dwn)  but in this case one should use nsppol=1.
!!      -----------------------------------           As a consequence the entire matrix is calculated and stored on file.
!!
!!  The coupling block is given by:
!!
!!      |  (c'v' up)   |    (c'v dwn)     |
!!      -----------------------------------           where v_{-+} = v_{+-}^t when the momentum of the photon is neglected.
!!      | [-W+v]++     |      v+-         | (vc up)   Also in this case the entire matrix v_{+-} has to be calculated
!!  C = -----------------------------------           and stored on file.
!!      |     v-+      |    [-W+v]--      | (vc dwn)
!!      -----------------------------------
!!
!! PARENTS
!!      exc_build_ham
!!
!! CHILDREN
!!
!! SOURCE

subroutine exc_build_block(BSp,Cryst,Kmesh,Qmesh,ktabr,Gsph_x,Gsph_c,Vcp,Wfd,W,Hdr_bse,&
&  nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff,rhxtwg_q0,is_resonant,fname)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot_osc
 character(len=*),intent(in) :: fname
 logical,intent(in) :: is_resonant
 type(excparam),intent(in) :: BSp
 type(screen_t),intent(inout) :: W
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(crystal_t),intent(in) :: Cryst
 type(vcoul_t),intent(in) :: Vcp
 type(gsphere_t),intent(in) :: Gsph_x,Gsph_c
 type(Pseudopotential_type),intent(in) :: Psps
 type(Hdr_type),intent(inout) :: Hdr_bse
 type(pawang_type),intent(in) :: Pawang
 type(wfd_t),target,intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfft_osc(18)
 integer,intent(in) :: ktabr(nfftot_osc,Kmesh%nbz)
 complex(gwpc),intent(in) :: rhxtwg_q0(BSp%npweps,BSp%lomo_min:BSp%humo_max,BSp%lomo_min:BSp%humo_max,Wfd%nkibz,Wfd%nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Wfd%usepaw)
 type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: map2sphere=1,ndat1=1,master=0
 integer(i8b) :: bsize_my_block
 integer :: nspinor,nsppol,ISg,mpi_err,tmp_size,ngx
 integer :: ik_bz,ikp_bz,col_glob,itpk_min,itpk_max
 integer :: dim_rtwg,bsh_unt,ncol,dump_unt,npweps
#ifdef HAVE_MPI_IO
 integer :: amode,mpi_fh,hmat_type,offset_err,old_type
 integer(XMPI_OFFSET_KIND) :: ehdr_offset,my_offset
 logical,parameter :: is_fortran_file=.TRUE.
#endif
 integer :: neh1,neh2,ig,nblocks
 integer :: ik_ibz,itim_k,ikp_ibz,itim_kp,isym_k,isym_kp
 integer :: iq_bz,iq_ibz,isym_q,itim_q,iqbz0,rank
 integer :: iv,ivp,ic,icp,jj,nrows,sender,my_ncols
 integer :: use_padfft,prev_nrows,spin1,spin2,block
 integer :: ierr,nproc,my_rank,mgfft_osc,fftalga_osc,comm
 integer(i8b) :: tot_nels,prev_nels,prev_ncols,nels,ir,it,itp,ist,iend,my_hsize
 real(dp) :: faq,kx_fact,cputime,walltime,gflops
 complex(spc) :: http,ctemp
 complex(dpc) :: ph_mkpt,ph_mkt,ene_t,ene_tp
 logical,parameter :: with_umklp=.FALSE.
 logical :: use_mpiio,do_coulomb_term,do_exchange_term,w_is_diagonal,isirred
 logical :: is_qeq0
 character(len=500) :: msg
!arrays
 integer :: bidx(2,4),g0(3),spin_ids(2,3)
 integer(i8b) :: nels_block(3)
 integer :: my_cols(2),my_rows(2),proc_end(2),proc_start(2)
 integer :: my_extrema(2,2),sender_extrema(2,2),my_starts(2),my_ends(2)
 integer,allocatable :: igfftg0(:),ktabr_k(:),ktabr_kp(:),id_tab(:)
 integer,allocatable :: ncols_of(:)
 integer(i8b),allocatable :: t_start(:),t_stop(:),hsize_of(:)
 integer,allocatable :: col_start(:),col_stop(:)
 integer,allocatable :: gbound(:,:)
 real(dp) :: kbz(3),kpbz(3),qbz(3),spinrot_k(4),spinrot_kp(4),kmkp(3),tsec(2)
 complex(dpc),allocatable :: my_bsham(:),buffer(:),buffer_2d(:,:),my_kxssp(:,:),prev_col(:)
!DBYG
 complex(dpc),allocatable :: acoeffs(:),bcoeffs(:),ccoeffs(:) ! Coeff of W = a/q^2 + b/q + c
 integer :: a_unt, b_unt, c_unt
 complex(dpc) :: aatmp, bbtmp, cctmp
 complex(gwpc),allocatable :: aa_vpv(:),aa_cpc(:),aa_ctccp(:)
 complex(gwpc),allocatable :: bb_vpv1(:),bb_cpc1(:),bb_ctccp1(:)
 complex(gwpc),allocatable :: bb_vpv2(:),bb_cpc2(:),bb_ctccp2(:)
 complex(gwpc),allocatable :: cc_vpv(:),cc_cpc(:),cc_ctccp(:)
 complex(dpc),allocatable :: abuffer(:),aprev_col(:)
 complex(dpc),allocatable :: bbuffer(:),bprev_col(:)
 complex(dpc),allocatable :: cbuffer(:),cprev_col(:)
 character(len=fnlen) :: tmpfname
 integer :: ii
!END DBYG
 complex(gwpc),allocatable :: vc_sqrt_qbz(:)
 complex(gwpc),allocatable :: rhotwg1(:),rhotwg2(:),rhxtwg_vpv(:),rhxtwg_cpc(:),ctccp(:)
 complex(gwpc),target,allocatable :: ur_ckp(:),ur_vkp(:),ur_vk(:),ur_ck(:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: ptur_ckp(:),ptur_vkp(:),ptur_vk(:),ptur_ck(:)
 type(pawcprj_type),target,allocatable :: Cp_tmp1(:,:),Cp_tmp2(:,:)
 type(pawcprj_type),target,allocatable :: Cp_tmp3(:,:),Cp_tmp4(:,:)
 type(pawcprj_type),allocatable :: Cp_ckp(:,:),Cp_vkp(:,:)
 type(pawcprj_type),allocatable :: Cp_vk(:,:),Cp_ck(:,:)
 type(pawcprj_type),pointer :: ptcp_ckp(:,:),ptcp_vkp(:,:),ptcp_vk(:,:),ptcp_ck(:,:)
 type(pawpwij_t),allocatable :: Pwij_q(:)
#ifdef HAVE_MPI_IO
 integer(XMPI_OFFSET_KIND) :: tmp_off,my_offpad
 integer(XMPI_OFFSET_KIND),allocatable :: bsize_frecord(:),offset_of_block(:)
#endif
#ifdef DEV_MG_DEBUG_MODE
 integer,allocatable :: ttp_check(:,:)
#endif

!************************************************************************

 call timab(680,1,tsec)
 call timab(681,1,tsec)

 DBG_ENTER("COLL")

 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")
 ABI_CHECK(nfftot_osc==PRODUCT(ngfft_osc(1:3)),"mismatch in FFT size")

 if (Wfd%nsppol==2) then
   MSG_WARNING("nsppol==2 is still under testing")
 end if
 !
 ! MPI variables.
 comm    = Wfd%comm
 nproc   = Wfd%nproc
 my_rank = Wfd%my_rank

 !
 ! Basic constants.
 nspinor = Wfd%nspinor
 nsppol  = Wfd%nsppol
 dim_rtwg=1; faq = one/(Cryst%ucvol*Kmesh%nbz)
 npweps = Bsp%npweps
 !
 ! Prepare the FFT tables to have u(r) on the ngfft_osc mesh.
 mgfft_osc = MAXVAL(ngfft_osc(1:3))
 fftalga_osc = ngfft_osc(7)/100
 if ( ANY(ngfft_osc(1:3) /= Wfd%ngfft(1:3)) ) then
   call wfd%change_ngfft(Cryst,Psps,ngfft_osc)
 end if

 ABI_MALLOC(igfftg0,(npweps))
 ABI_MALLOC(ktabr_k,(nfftot_osc))
 ABI_MALLOC(ktabr_kp,(nfftot_osc))
 ABI_MALLOC(id_tab,(nfftot_osc))
 id_tab = (/(ic, ic=1,nfftot_osc)/)
 !
 ! Workspace arrays for wavefunctions and oscillator matrix elements.
 ABI_MALLOC(rhxtwg_vpv,(npweps))
 ABI_MALLOC(rhxtwg_cpc,(npweps))

 if (BSp%prep_interp) then
   call wrtout(std_out,"Preparing BSE interpolation","COLL")
   ABI_MALLOC(aa_vpv,(npweps))
   ABI_MALLOC(bb_vpv1,(npweps))
   ABI_MALLOC(bb_vpv2,(npweps))
   ABI_MALLOC(cc_vpv,(npweps))

   ABI_MALLOC(aa_cpc,(npweps))
   ABI_MALLOC(bb_cpc1,(npweps))
   ABI_MALLOC(bb_cpc2,(npweps))
   ABI_MALLOC(cc_cpc,(npweps))
 end if

 ABI_MALLOC(ur_ckp,(nspinor*nfftot_osc))
 ABI_MALLOC(ur_vkp,(nspinor*nfftot_osc))
 ABI_MALLOC(ur_ck ,(nspinor*nfftot_osc))
 ABI_MALLOC(ur_vk ,(nspinor*nfftot_osc))

 if (Wfd%usepaw==1) then
   ABI_DT_MALLOC(Cp_vk,(Wfd%natom,nspinor))
   call pawcprj_alloc(Cp_vk,0,Wfd%nlmn_atm)
   ABI_DT_MALLOC(Cp_ck,(Wfd%natom,nspinor))
   call pawcprj_alloc(Cp_ck,0,Wfd%nlmn_atm)
   ABI_DT_MALLOC(Cp_ckp,(Wfd%natom,nspinor))
   call pawcprj_alloc(Cp_ckp,0,Wfd%nlmn_atm)
   ABI_DT_MALLOC(Cp_vkp,(Wfd%natom,nspinor))
   call pawcprj_alloc(Cp_vkp,0,Wfd%nlmn_atm)

   ABI_DT_MALLOC(Cp_tmp1,(Wfd%natom,nspinor))
   call pawcprj_alloc(Cp_tmp1,0,Wfd%nlmn_atm)
   ABI_DT_MALLOC(Cp_tmp2,(Wfd%natom,nspinor))
   call pawcprj_alloc(Cp_tmp2,0,Wfd%nlmn_atm)
   ABI_DT_MALLOC(Cp_tmp3,(Wfd%natom,nspinor))
   call pawcprj_alloc(Cp_tmp3,0,Wfd%nlmn_atm)
   ABI_DT_MALLOC(Cp_tmp4,(Wfd%natom,nspinor))
   call pawcprj_alloc(Cp_tmp4,0,Wfd%nlmn_atm)
 end if
 !
 ! Identify the index of q==0
 iqbz0=0
 do iq_bz=1,Qmesh%nbz
   if (ALL(ABS(Qmesh%bz(:,iq_bz))<tol3)) iqbz0 = iq_bz
 end do
 ABI_CHECK(iqbz0/=0,"q=0 not found")
 !
 ! Treat the spin polarization.
 spin_ids(:,1) = (/1,1/)
 spin_ids(:,2) = (/2,2/)
 spin_ids(:,3) = (/1,2/)

 nblocks=1
 kx_fact=two
 nels_block(:)=0
 nels_block(1)=BSp%nreh(1)*(BSp%nreh(1)+1_i8b)/2
 tot_nels=nels_block(1)

 if (nsppol==2) then
   nblocks=3
   kx_fact=one
   nels_block(1) = BSp%nreh(1)*(BSp%nreh(1)+1_i8b)/2   ! Only the upper triangle for block 1 and 2
   nels_block(2) = BSp%nreh(2)*(BSp%nreh(2)+1_i8b)/2
   nels_block(3) = BSp%nreh(1)*BSp%nreh(2)*1_i8b       ! Note: Block 3 does not have symmetries.
   tot_nels= SUM(nels_block)
 end if
 !
 ! Distribute the calculation of the matrix elements among the nodes.
 ! * tstart and t_stop give the initial and final transition index treated by each node.
 ! * my_hsize is the number of transitions treated by this processor
 ! * my_cols(1:2) gives the initial and final column treated by this node.
 !
 use_mpiio=.FALSE.
#ifdef HAVE_MPI_IO
 use_mpiio = (nproc>1)
#endif
 use_mpiio=.FALSE.
 !use_mpiio=.TRUE.

 if (is_resonant) then
   if (use_mpiio) then
     write(msg,'(2a,f6.2,a)')&
&      ". Writing resonant excitonic Hamiltonian on file "//TRIM(fname)," via MPI-IO; file size= ",two*tot_nels*dpc*b2Gb," [Gb]."
   else
     write(msg,'(2a,f6.2,a)')&
&      ". Writing resonant excitonic Hamiltonian on file "//TRIM(fname),"; file size= ",two*dpc*tot_nels*b2Gb," [Gb]."
   end if
 else
   if (use_mpiio) then
     write(msg,'(2a,f6.2,a)')&
&      ". Writing coupling excitonic Hamiltonian on file "//TRIM(fname)," via MPI-IO; file size= ",tot_nels*2*dpc*b2Gb," [Gb]."
   else
     write(msg,'(2a,f6.2,a)')&
&      ". Writing coupling excitonic Hamiltonian on file "//TRIM(fname),"; file size= ",two*dpc*tot_nels*b2Gb," [Gb]."
   end if
 end if
 call wrtout(std_out,msg,"COLL",do_flush=.True.)
 call wrtout(ab_out,msg,"COLL",do_flush=.True.)
 !
 ! Master writes the BSE header with Fortran IO.
 if (my_rank==master) then

   if (open_file(fname,msg,newunit=bsh_unt,form="unformatted",action="write") /= 0) then
      MSG_ERROR(msg)
   end if
   call exc_write_bshdr(bsh_unt,Bsp,Hdr_bse)
   ! To force the writing (needed for MPI-IO).
   close(bsh_unt)

   if (.not.use_mpiio) then ! Reopen the file and skip the header.
     if (open_file(fname,msg,newunit=bsh_unt,form="unformatted",action="readwrite") /= 0) then
        MSG_ERROR(msg)
     end if
     call exc_skip_bshdr(bsh_unt,ierr)
   end if

   if (BSp%prep_interp) then
     tmpfname = fname
     ii = LEN_TRIM(fname)
     tmpfname(ii-2:ii+1) = 'ABSR'
     if (open_file(tmpfname,msg,newunit=a_unt,form='unformatted',action="write") /= 0) then
       MSG_ERROR(msg)
     end if
     tmpfname(ii-2:ii+1) = 'BBSR'
     if (open_file(tmpfname,msg,newunit=b_unt,form='unformatted',action="write") /= 0) then
       MSG_ERROR(msg)
     end if
     tmpfname(ii-2:ii+1) = 'CBSR'
     if (open_file(tmpfname,msg,newunit=c_unt,form='unformatted',action="write") /= 0) then
       MSG_ERROR(msg)
     end if
     call exc_write_bshdr(a_unt,Bsp,Hdr_bse)
     call exc_write_bshdr(b_unt,Bsp,Hdr_bse)
     call exc_write_bshdr(c_unt,Bsp,Hdr_bse)
     close(a_unt)
     close(b_unt)
     close(c_unt)
     if (.not.use_mpiio) then ! Reopen the file and skip the header.
       tmpfname(ii-2:ii+1) = 'ABSR'
       if (open_file(tmpfname,msg,newunit=a_unt,form='unformatted',action="readwrite") /= 0) then
          MSG_ERROR(msg)
       end if
       call exc_skip_bshdr(a_unt,ierr)
       tmpfname(ii-2:ii+1) = 'BBSR'
       if (open_file(tmpfname,msg,newunit=b_unt,form='unformatted',action="readwrite") /= 0) then
          MSG_ERROR(msg)
       end if
       call exc_skip_bshdr(b_unt,ierr)
       tmpfname(ii-2:ii+1) = 'CBSR'
       if (open_file(tmpfname,msg,newunit=c_unt,form='unformatted',action="readwrite") /= 0) then
          MSG_ERROR(msg)
       end if
       call exc_skip_bshdr(c_unt,ierr)
     end if
   end if
 end if

 call xmpi_barrier(comm)

 if (use_mpiio) then
#ifdef HAVE_MPI_IO
   ! Open the file with MPI-IO
   amode = MPI_MODE_RDWR
   !amode = MPI_MODE_CREATE + MPI_MODE_RDWR,

   call MPI_FILE_OPEN(comm, fname, amode, MPI_INFO_NULL, mpi_fh, mpi_err)
   ABI_CHECK_MPI(mpi_err,"opening: "//TRIM(fname))

   ! Skip the header.
   call exc_skip_bshdr_mpio(mpi_fh,xmpio_collective,ehdr_offset)

   ! Precompute the offset of the each block including the Fortran markers.
   ABI_MALLOC(offset_of_block,(nblocks))
   offset_of_block(1) = ehdr_offset
   do block=2,nblocks
     tmp_off = offset_of_block(block-1) + nels_block(block-1)*xmpi_bsize_dpc
     tmp_off = tmp_off + Bsp%nreh(block-1)*2*xmpio_bsize_frm  ! markers.
     offset_of_block(block) = tmp_off
   end do
#endif
 end if

 call timab(681,2,tsec)

 do block=1,nsppol
   !
   ! Indices used to loop over bands.
   ! bidx contains the starting and final indices used to loop over bands.
   !
   !      (b3,b4)
   !         |... ...|
   ! (b1,b2) |... ...|
   !
   ! Resonant matrix is given by
   !      (v',c')
   !       |... ...|
   ! (v,c) |... ...|
   !
   ! Coupling matrix is given by
   !       (c',v')
   !       |... ...|
   ! (v,c) |... ...|

   if (is_resonant) then
     bidx(:,1) = [BSp%lomo_spin(block),BSp%homo_spin(block)] ! range for b1
     bidx(:,2) = [BSp%lumo_spin(block),BSp%humo_spin(block)] ! range for b2
     bidx(:,3) = [BSp%lomo_spin(block),BSp%homo_spin(block)] ! range for b3
     bidx(:,4) = [BSp%lumo_spin(block),BSp%humo_spin(block)] ! range for b4
   else
     bidx(:,1) = [BSp%lomo_spin(block),BSp%homo_spin(block)] ! range for b1
     bidx(:,2) = [BSp%lumo_spin(block),BSp%humo_spin(block)] ! range for b2
     bidx(:,3) = [BSp%lumo_spin(block),BSp%humo_spin(block)] ! range for b3
     bidx(:,4) = [BSp%lomo_spin(block),BSp%homo_spin(block)] ! range for b4
   end if

   spin1 = spin_ids(1,block)
   spin2 = spin_ids(2,block)

   do_coulomb_term  = (Bsp%use_coulomb_term .and. (spin1==spin2))
   do_exchange_term = (Bsp%exchange_term>0)
   w_is_diagonal    = BSp%use_diagonal_Wgg
   !
   ! Distribution of the matrix elements among the nodes.
   ! Note that rank0 will get the first transitions.
   nels=nels_block(block)
   ABI_MALLOC(t_start,(0:nproc-1))
   ABI_MALLOC(t_stop,(0:nproc-1))
   call xmpi_split_work2_i8b(nels,nproc,t_start,t_stop)

   ABI_MALLOC(hsize_of,(0:nproc-1))
   hsize_of=0
   do rank=0,nproc-1
     if (t_stop(rank)>=t_start(rank)) hsize_of(rank) = t_stop(rank)-t_start(rank)+1
     !write(std_out,*)"nels",nels,hsize_of(rank)
   end do

   my_hsize = hsize_of(my_rank)
   if (my_hsize<=0) then
     write(msg,'(a,i0)')"Wrong number of transitions: my_hsize= ",my_hsize
     MSG_ERROR(msg)
   end if
   if (my_hsize /= INT(my_hsize,KIND=i4b)) then
     write(msg,'(a,i0)')"Size of local block too large for a default integer, Increase the number of CPUs: my_hsize= ",my_hsize
     MSG_ERROR(msg)
   end if

   my_cols=0
   do itp=1,Bsp%nreh(block)
     do it=1,itp
       ir = it + itp*(itp-1_i8b)/2
       if (ir==t_start(my_rank)) then
         my_rows(1) = it
         my_cols(1) = itp
       end if
       if (ir==t_stop(my_rank)) then
         my_rows(2) = it
         my_cols(2) = itp
       end if
     end do
   end do

   my_starts = [my_rows(1),my_cols(1)]
   my_ends   = [my_rows(2),my_cols(2)]
   !
   ! Announce the treatment of submatrix treated by each node.
   bsize_my_block = 2*dpc*my_hsize
   write(msg,'(4(a,i0))')' Treating ',my_hsize,'/',nels,' matrix elements, from column ',my_cols(1),' up to column ',my_cols(2)
   call wrtout(std_out,msg,'PERS')

   if (is_resonant) then
     write(msg,'(a,f8.1,a)')' Calculating resonant blocks. Memory required: ',bsize_my_block*b2Mb,' [Mb] <<< MEM'
   else
     write(msg,'(a,f8.1,a)')' Calculating coupling blocks. Memory required: ',bsize_my_block*b2Mb,' [Mb] <<< MEM'
   end if
   call wrtout(std_out,msg,"COLL")

   ! Allocate big (scalable) buffer to store the BS matrix on this node.
   ABI_MALLOC_OR_DIE(my_bsham,(t_start(my_rank):t_stop(my_rank)), ierr)

   if (BSp%prep_interp) then
     ! Allocate big (scalable) buffers to store a,b,c coeffients
     ABI_MALLOC_OR_DIE(acoeffs,(t_start (my_rank):t_stop(my_rank)), ierr)

     ABI_MALLOC_OR_DIE(bcoeffs,(t_start(my_rank):t_stop(my_rank)), ierr)

     ABI_MALLOC_OR_DIE(ccoeffs,(t_start(my_rank):t_stop(my_rank)), ierr)
   end if

   if (do_coulomb_term) then ! Construct Coulomb term.

     call timab(682,1,tsec) ! exc_build_ham(Coulomb)

     write(msg,'(a,2i2,a)')" Calculating direct Coulomb term for (spin1,spin2) ",spin1,spin2," using full W_{GG'} ..."
     if (w_is_diagonal) then
        write(msg,'(a,2i2,a)')&
&        " Calculating direct Coulomb term for (spin1, spin2) ",spin1,spin2," using diagonal approximation for W_{GG'} ..."
     end if
     call wrtout(std_out,msg,"COLL")

     ABI_MALLOC(ctccp,(npweps))

     if (BSp%prep_interp) then
       ABI_MALLOC(aa_ctccp,(npweps))
       ABI_MALLOC(bb_ctccp1,(npweps))
       ABI_MALLOC(bb_ctccp2,(npweps))
       ABI_MALLOC(cc_ctccp,(npweps))
     end if

     ABI_MALLOC(vc_sqrt_qbz,(npweps))

#ifdef DEV_MG_DEBUG_MODE
     ABI_MALLOC(ttp_check,(BSp%nreh(block),BSp%nreh(block)))
     ttp_check=0
#endif

     do ikp_bz=1,BSp%nkbz ! Loop over kp
       ! NOTE: this way of looping is good for bulk but it's not optimal in the
       !       case of systems sampled only at Gamma e.g. isolated systems in which
       !       one should take advantage of Hermiticity by looping over c-v !!!!

       ! Check whether (vp,cp,ikp_bz,spin2) belongs to the set of columns treated by me for some vp,cp
       ! Be careful since vcks2t contains zeros corresponding to transitions that should be skipped.
       itpk_min = MINVAL(Bsp%vcks2t(:,:,ikp_bz,spin2), MASK=(Bsp%vcks2t(:,:,ikp_bz,spin2)>0) )
       itpk_max = MAXVAL(Bsp%vcks2t(:,:,ikp_bz,spin2))
       if (my_cols(2)<itpk_min .or. my_cols(1)>itpk_max) CYCLE

       write(msg,'(3(a,i0))')" status: ",ikp_bz,"/",BSp%nkbz," done by node ",my_rank
       call wrtout(std_out,msg,"PERS",do_flush=.True.)

       ! * Get ikp_ibz, non-symmorphic phase, ph_mkpt, and symmetries from ikp_bz.
       call get_BZ_item(Kmesh,ikp_bz,kpbz,ikp_ibz,isym_kp,itim_kp,ph_mkpt,isirred=isirred)
       !ABI_CHECK(isirred,"not irred!")
       !ABI_CHECK(ph_mkpt == cone, "Wrong phase!")

       ktabr_kp(:) = ktabr(:,ikp_bz)
       spinrot_kp(:)=Cryst%spinrot(:,isym_kp)
       !ABI_CHECK(ALL(ktabr_kp == id_tab), "wrong tab")

       do ik_bz=1,ikp_bz ! Loop over k
         !
         ! * Get ik_ibz, non-symmorphic phase, ph_mkt, and symmetries from ik_bz
         call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k,ph_mkt,isirred=isirred)
         !ABI_CHECK(isirred,"not irred!")
         !ABI_CHECK(ph_mkt == cone, "Wrong phase!")

         ktabr_k(:) = ktabr(:,ik_bz)
         spinrot_k(:)=Cryst%spinrot(:,isym_k)
         !ABI_CHECK(ALL(ktabr_k == id_tab), "wrong tab")
         !if(itim_k==2) CYCLE ! time-reversal or not
         !
         ! * Find q = K-KP-G0 in the full BZ.
         kmkp = Kmesh%bz(:,ik_bz) - Kmesh%bz(:,ikp_bz)
         call findqg0(iq_bz,g0,kmkp,Qmesh%nbz,Qmesh%bz,BSp%mG0)

         ! Evaluate the tables needed for the padded FFT performed in rhotwg. Note that we have
         ! to pass G-G0 to sphereboundary instead of G as we need FFT results on the shifted G-sphere,
         ! If Gamma is not inside G-G0 one has to disable FFT padding as sphereboundary will give wrong tables.
         ! * Get the G-G0 shift for the FFT of the oscillators.
         !
         ABI_MALLOC(gbound,(2*mgfft_osc+8,2))
         call gsph_fft_tabs(Gsph_c,g0,mgfft_osc,ngfft_osc,use_padfft,gbound,igfftg0)
#ifdef FC_IBM
 ! XLF does not deserve this optimization (problem with [v67mbpt][t03])
 use_padfft = 0
#endif
         if ( ANY(fftalga_osc == (/2,4/)) ) use_padfft=0 ! Pad-FFT is not coded in rho_tw_g
         if (use_padfft==0) then
           ABI_FREE(gbound)
           ABI_MALLOC(gbound,(2*mgfft_osc+8,2*use_padfft))
         end if
         !
         ! * Get iq_ibz, and symmetries from iq_bz
         call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)
         is_qeq0 = (normv(qbz,Cryst%gmet,'G')<GW_TOLQ0)

         ! Symmetrize em1(omega=0)
         call screen_symmetrizer(W,iq_bz,Cryst,Gsph_c,Qmesh,Vcp)
         !
         ! * Set up table of |q_BZ+G|
         if (iq_ibz==1) then
           do ig=1,npweps
             isg = Gsph_c%rottb(ig,itim_q,isym_q)
             vc_sqrt_qbz(isg)=Vcp%vcqlwl_sqrt(ig,1)
           end do
         else
           do ig=1,npweps
             isg = Gsph_c%rottb(ig,itim_q,isym_q)
             vc_sqrt_qbz(isg) = Vcp%vc_sqrt(ig,iq_ibz)
           end do
         end if

         ! === Evaluate oscillator matrix elements ===
         ! * $ <phj/r|e^{-i(q+G)}|phi/r> - <tphj/r|e^{-i(q+G)}|tphi/r> $ in packed form.
         if (Wfd%usepaw==1.and.ik_bz/=ikp_bz) then
           ABI_DT_MALLOC(Pwij_q,(Cryst%ntypat))
           call pawpwij_init(Pwij_q,npweps,Qmesh%bz(:,iq_bz),Gsph_c%gvec,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
         end if

         ! =======================================
         ! === Loop over the four band indeces ===
         ! =======================================
         do ic=bidx(1,2),bidx(2,2) !do ic=BSp%lumo,BSp%nbnds

           if (wfd%ihave_ur(ic,ik_ibz,spin1,how="Stored")) then
             ptur_ck => Wfd%Wave(ic,ik_ibz,spin1)%ur
           else
             call wfd%get_ur(ic,ik_ibz,spin1,ur_ck)
             ptur_ck => ur_ck
           end if
           !
           ! Get cprj for this (c,kbz,s1) in the BZ.
           ! phase due to the umklapp G0 in k-q is already included.
           if (Wfd%usepaw==1) then
             if (wfd%ihave_cprj(ic,ik_ibz,spin1,how="Stored")) then
               ptcp_ck => Wfd%Wave(ic,ik_ibz,spin1)%cprj
             else
               call wfd%get_cprj(ic,ik_ibz,spin1,Cryst,Cp_tmp1,sorted=.FALSE.)
               ptcp_ck => Cp_tmp1
             end if
             call paw_symcprj_op(ik_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,ptcp_ck,Cp_ck)
           end if

           do icp=bidx(1,4),bidx(2,4)  !do icp=BSp%lumo,BSp%nbnds
             ! Calculate matrix-elements rhxtwg_cpc
             !
             if (ik_bz==ikp_bz) then ! Already in memory.
               rhxtwg_cpc(:) = sym_rhotwgq0(itim_k,isym_k,dim_rtwg,npweps,rhxtwg_q0(:,icp,ic,ik_ibz,spin1),Gsph_c)

             else
               ! Calculate matrix element from wfr.
               ! TODO: change the order of the loops.

               if (wfd%ihave_ur(icp,ikp_ibz,spin2,how="Stored")) then
                 ptur_ckp => Wfd%Wave(icp,ikp_ibz,spin2)%ur
               else
                 call wfd%get_ur(icp,ikp_ibz,spin2,ur_ckp)
                 ptur_ckp => ur_ckp
               end if

               ! Load cprj for this (c,k,s2) in the BZ.
               ! Do not care about umklapp G0 in k-q as the phase is already included.
               if (Wfd%usepaw==1) then
                 if (wfd%ihave_cprj(icp,ikp_ibz,spin2,how="Stored")) then
                   ptcp_ckp =>  Wfd%Wave(icp,ikp_ibz,spin2)%cprj
                 else
                   call wfd%get_cprj(icp,ikp_ibz,spin2,Cryst,Cp_tmp2,sorted=.FALSE.)
                   ptcp_ckp =>  Cp_tmp2
                 end if
                 call paw_symcprj_op(ikp_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,ptcp_ckp,Cp_ckp)
               end if

               call rho_tw_g(nspinor,npweps,nfftot_osc,ndat1,ngfft_osc,map2sphere,use_padfft,igfftg0,gbound,&
&                ptur_ckp,itim_kp,ktabr_kp,ph_mkpt,spinrot_kp,&
&                ptur_ck ,itim_k ,ktabr_k ,ph_mkt ,spinrot_k ,&
&                dim_rtwg,rhxtwg_cpc)

               if (Wfd%usepaw==1) then ! Add PAW onsite contribution.
                 call paw_rho_tw_g(npweps,dim_rtwg,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_c%gvec,&
&                  Cp_ckp,Cp_ck,Pwij_q,rhxtwg_cpc)
               end if
             end if

             if (BSp%prep_interp) then
               aa_cpc = rhxtwg_cpc
               aa_cpc(2:) = czero
               bb_cpc1 = vc_sqrt_qbz*rhxtwg_cpc
               bb_cpc1(1) = czero
               bb_cpc2 = rhxtwg_cpc
               bb_cpc2(2:) = czero

               if(ik_bz == ikp_bz) then
                  ! Enforce orthogonality of the wavefunctions.
                  if(icp == ic) then
                    aa_cpc(1) = cone
                    bb_cpc2(1) = cone
                  else
                    aa_cpc(1) = czero
                    bb_cpc2(1) = czero
                  end if
               end if

               ! MG TODO: a does not require a call to w0gemv
               call screen_w0gemv(W,"C",npweps,nspinor,w_is_diagonal,cone_gw,czero_gw,aa_cpc,aa_ctccp)
               call screen_w0gemv(W,"C",npweps,nspinor,w_is_diagonal,cone_gw,czero_gw,bb_cpc1,bb_ctccp1)
               call screen_w0gemv(W,"C",npweps,nspinor,w_is_diagonal,cone_gw,czero_gw,bb_cpc2,bb_ctccp2)

               cc_cpc = vc_sqrt_qbz*rhxtwg_cpc
               cc_cpc(1) = czero

               call screen_w0gemv(W,"C",npweps,nspinor,w_is_diagonal,cone_gw,czero_gw,cc_cpc,cc_ctccp)
             end if

             ! Prepare sum_GG' rho_c'c*(G) W_qbz(G,G') rho_v'v(G')
             ! First sum on G: sum_G rho_c'c(G) W_qbz*(G,G') (W_qbz conjugated)
             rhxtwg_cpc = rhxtwg_cpc * vc_sqrt_qbz
             call screen_w0gemv(W,"C",npweps,nspinor,w_is_diagonal,cone_gw,czero_gw,rhxtwg_cpc,ctccp)

             do iv=bidx(1,1),bidx(2,1)    !do iv=BSp%lomo,BSp%homo
               it = BSp%vcks2t(iv,ic,ik_bz,spin1); if (it==0) CYCLE ! ir-uv-cutoff
               ene_t = BSp%Trans(it,spin1)%en

               ! TODO: use this but change the order of the loops.
               if (wfd%ihave_ur(iv,ik_ibz,spin1,how="Stored")) then
                 ptur_vk => Wfd%Wave(iv,ik_ibz,spin1)%ur
               else
                 call wfd%get_ur(iv,ik_ibz,spin1,ur_vk)
                 ptur_vk => ur_vk
               end if
               !
               ! Load cprj for this (v,k,s1) in the BZ.
               ! * Do not care about umklapp G0 in k-q as the phase is already included.
               if (Wfd%usepaw==1) then
                 if (wfd%ihave_cprj(iv,ik_ibz,spin1,how="Stored")) then
                   ptcp_vk => Wfd%Wave(iv,ik_ibz,spin1)%cprj
                 else
                   call wfd%get_cprj(iv,ik_ibz,spin1,Cryst,Cp_tmp3,sorted=.FALSE.)
                   ptcp_vk => Cp_tmp3
                 end if
                 call paw_symcprj_op(ik_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,ptcp_vk,Cp_vk)
               end if

               do ivp=bidx(1,3),bidx(2,3) !do ivp=BSp%lomo,BSp%homo

                 if (is_resonant) then
                   itp = BSp%vcks2t(ivp,icp,ikp_bz,spin2)
                 else ! have to exchange band indeces
                   itp = BSp%vcks2t(icp,ivp,ikp_bz,spin2)
                 end if

                 if (itp==0) CYCLE ! ir-uv-cutoff

                 ! FIXME Temporary work around, when ikp_bz == ik it might happen that itp<it
                 ! should rewrite the loops using contracted k-dependent indeces for bands
                 if (itp<it) CYCLE

                 ir = it + itp*(itp-1)/2
                 if (ir<t_start(my_rank).or.ir>t_stop(my_rank)) CYCLE

                 ene_tp = BSp%Trans(itp,spin2)%en

                 ! ============================================
                 ! === Calculate matrix elements rhxtwg_vpv ===
                 ! ============================================
                 if (ik_bz==ikp_bz) then
                   ! Already in memory.
                   rhxtwg_vpv(:) = sym_rhotwgq0(itim_k,isym_k,dim_rtwg,npweps,rhxtwg_q0(:,ivp,iv,ik_ibz,spin1),Gsph_c)

                 else
                   ! Calculate matrix element from wfr.
                   if (wfd%ihave_ur(ivp,ikp_ibz,spin2,how="Stored")) then
                     ptur_vkp => Wfd%Wave(ivp,ikp_ibz,spin2)%ur
                   else
                     call wfd%get_ur(ivp,ikp_ibz,spin2,ur_vkp)
                     ptur_vkp => ur_vkp
                   end if
                   !
                   ! Load cprj for this (vp,kp,s2) in the BZ.
                   ! Do not care about umklapp G0 in k-q as the phase is already included.
                   if (Wfd%usepaw==1) then
                     if (wfd%ihave_cprj(ivp,ikp_ibz,spin2,how="Stored")) then
                       ptcp_vkp =>  Wfd%Wave(ivp,ikp_ibz,spin2)%cprj
                     else
                       call wfd%get_cprj(ivp,ikp_ibz,spin2,Cryst,Cp_tmp4,sorted=.FALSE.)
                       ptcp_vkp => Cp_tmp4
                     end if
                     call paw_symcprj_op(ikp_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,ptcp_vkp,Cp_vkp)
                   end if

                   call rho_tw_g(nspinor,npweps,nfftot_osc,ndat1,ngfft_osc,map2sphere,use_padfft,igfftg0,gbound,&
&                     ptur_vkp,itim_kp,ktabr_kp,ph_mkpt,spinrot_kp,&
&                     ptur_vk ,itim_k ,ktabr_k ,ph_mkt ,spinrot_k ,&
&                     dim_rtwg,rhxtwg_vpv)

                   if (Wfd%usepaw==1) then ! Add PAW onsite contribution.
                     call paw_rho_tw_g(npweps,dim_rtwg,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,&
&                      Gsph_c%gvec,Cp_vkp,Cp_vk,Pwij_q,rhxtwg_vpv)
                   end if
                 end if

                 ! Index in the global Hamiltonian matrix.
                 ir = it + itp*(itp-1_i8b)/2

                 if (ir<t_start(my_rank).or.ir>t_stop(my_rank)) then
                   write(msg,'(a,3(1x,i0))')" Gonna SIGFAULT, ir, t_start, t_stop ",ir,t_start(my_rank),t_stop(my_rank)
                   MSG_ERROR(msg)
                 end if
                 !ABI_CHECK(itp >= it,"itp < it")

                 if (BSp%prep_interp) then
                   ! Save a,b, c coefficients.
                   aa_vpv = rhxtwg_vpv
                   aa_vpv(2:) = czero
                   bb_vpv1 = rhxtwg_vpv
                   bb_vpv1(2:) = czero
                   bb_vpv2 = vc_sqrt_qbz*rhxtwg_vpv
                   bb_vpv2(1) = czero

                   if (ik_bz == ikp_bz) then
                      ! Enforce orthogonality of the wavefunctions.
                      if (ivp == iv) then
                        aa_vpv(1) = cone
                        bb_vpv1(1) = cone
                      else
                        aa_vpv(1) = czero
                        bb_vpv1(1) = czero
                      end if
                   end if

                   cc_vpv = vc_sqrt_qbz*rhxtwg_vpv
                   cc_vpv(1) = czero

                   aatmp = -faq * xdotc(npweps,aa_ctccp,1,aa_vpv,1)
                   bbtmp = -faq * xdotc(npweps,bb_ctccp1,1,bb_vpv1,1)-faq*xdotc(npweps,bb_ctccp2,1,bb_vpv2,1)
                   cctmp = -faq * xdotc(npweps,cc_ctccp,1,cc_vpv,1)

                   acoeffs(ir) = aatmp
                   bcoeffs(ir) = bbtmp
                   ccoeffs(ir) = cctmp
                 end if

                 ! sum_G2 rho_c'c(G) W_qbz(G,G') rho_v'v(G')
                 rhxtwg_vpv = vc_sqrt_qbz * rhxtwg_vpv
                 http = - faq * xdotc(npweps,ctccp,1,rhxtwg_vpv,1)

                 ! Save result taking into account the symmetry of the matrix.
                 ! Note that the diagonal of the resonant block is not forced to be real
                 my_bsham(ir) = http

#ifdef DEV_MG_DEBUG_MODE
                 ttp_check(it,itp) = ttp_check(it,itp)+1
#endif
               end do !ivp
             end do !iv
           end do !icp
         end do !ic

         ABI_FREE(gbound)

         if (Wfd%usepaw==1.and.ik_bz/=ikp_bz) then ! Free the onsite contribution for this q.
           call pawpwij_free(Pwij_q)
           ABI_DT_FREE(Pwij_q)
         end if

       end do ! ik_bz
     end do ! Fat loop over ikp_bz

#ifdef DEV_MG_DEBUG_MODE
     do itp=1,BSp%nreh(block)
       do it=1,BSp%nreh(block)
        ir = it + itp*(itp-1_i8b)/2
         if (itp>=it .and. ttp_check(it,itp) /= 1) then
           if (ir>=t_start(my_rank).and.ir<=t_stop(my_rank)) then
             write(std_out,*)"WARN: upper triangle is not 1 ",it,itp,ttp_check(it,itp)
             write(std_out,*)TRIM(repr_trans(Bsp%Trans(it ,spin1)))
             write(std_out,*)TRIM(repr_trans(Bsp%Trans(itp,spin2)))
           end if
         end if
         if (itp< it .and. ttp_check(it,itp) /= 0) then
           write(std_out,*)"WARN: then lower triangle is not 0 ",it,itp,ttp_check(it,itp)
           write(std_out,*)TRIM(repr_trans(Bsp%Trans(it ,spin1)))
           write(std_out,*)TRIM(repr_trans(Bsp%Trans(itp,spin2)))
         end if
       end do
     end do
     ierr = SUM(SUM(ttp_check,DIM=2),DIM=1)
     if (ierr/=my_hsize) then
       write(msg,'(a,2i0)')"ierr/=my_hsize",ierr,my_hsize
       MSG_ERROR(msg)
     end if
     ABI_FREE(ttp_check)
#endif

     ABI_FREE(ctccp)
     if(Bsp%prep_interp) then
       ABI_FREE(aa_ctccp)
       ABI_FREE(bb_ctccp1)
       ABI_FREE(bb_ctccp2)
       ABI_FREE(cc_ctccp)
     end if

     ABI_FREE(vc_sqrt_qbz)
     call wrtout(std_out,' Coulomb term completed',"COLL")

     call timab(682,2,tsec) ! exc_build_ham(Coulomb)
   end if ! do_coulomb_term
   !
   ! =====================
   ! === Exchange term ===
   ! =====================
   ! TODO might add treatment of <psi|q+G|psi> for q+G -> 0
   ! TODO might used enlarged G-sphere for better convergence.
   if (do_exchange_term) then

     !call exc_build_v(spin1,spin2,nsppol,npweps,Bsp,Cryst,Kmesh,Qmesh,Gsph_x,Gsph_c,Vcp,&
     ! &  is_resonant,rhxtwg_q0,nproc,my_rank,t_start,t_stop,my_bsham,comm)

     call timab(683,1,tsec) ! exc_build_ham(exchange)

     write(msg,'(a,2i2,a)')" Calculating exchange term for (spin1,spin2) ",spin1,spin2," ..."
     call wrtout(std_out,msg,"COLL")

     ABI_MALLOC(rhotwg1,(npweps))
     ABI_MALLOC(rhotwg2,(npweps))

     ngx = Gsph_x%ng
     ABI_MALLOC(vc_sqrt_qbz,(ngx))

     ! * Get iq_ibz, and symmetries from iq_bz.
     iq_bz = iqbz0 ! q = 0 -> iqbz0
     call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)

     ! * Set up table of |q(BZ)+G|
     if (iq_ibz==1) then
       do ig=1,ngx
         ISg = Gsph_x%rottb(ig,itim_q,isym_q)
         vc_sqrt_qbz(ISg)=Vcp%vcqlwl_sqrt(ig,1)
       end do
     else
        MSG_ERROR("iq_ibz should be 1")
     end if

     do itp=1,BSp%nreh(block) ! Loop over transition tp = (kp,vp,cp,spin2)

       if (itp<my_cols(1) .or. itp>my_cols(2)) CYCLE ! I dont have this column.
       ene_tp = Bsp%Trans(itp,spin2)%en
       ikp_bz = Bsp%Trans(itp,spin2)%k
       ivp    = Bsp%Trans(itp,spin2)%v
       icp    = Bsp%Trans(itp,spin2)%c

       ikp_ibz = Kmesh%tab (ikp_bz)
       isym_kp = Kmesh%tabo(ikp_bz)
       itim_kp = (3-Kmesh%tabi(ikp_bz))/2

       if (is_resonant) then
         rhotwg2(:) = sym_rhotwgq0(itim_kp,isym_kp,dim_rtwg,npweps,rhxtwg_q0(:,ivp,icp,ikp_ibz,spin2),Gsph_c)
       else ! Code for coupling block.
         rhotwg2(:) = sym_rhotwgq0(itim_kp,isym_kp,dim_rtwg,npweps,rhxtwg_q0(:,icp,ivp,ikp_ibz,spin2),Gsph_c)
       end if
       !
       ! Multiply by the Coulomb term.
        do ig=2,npweps
          rhotwg2(ig) = rhotwg2(ig) * vc_sqrt_qbz(ig) * vc_sqrt_qbz(ig)
        end do

       do it=1,itp ! Loop over transition t = (k,v,c,spin1)
         ir = it + itp*(itp-1_i8b)/2
         if (ir<t_start(my_rank) .or. ir>t_stop(my_rank)) CYCLE

         ene_t = Bsp%Trans(it,spin1)%en
         ik_bz = Bsp%Trans(it,spin1)%k
         iv    = Bsp%Trans(it,spin1)%v
         ic    = Bsp%Trans(it,spin1)%c

         ik_ibz = Kmesh%tab(ik_bz)
         isym_k = Kmesh%tabo(ik_bz)
         itim_k = (3-Kmesh%tabi(ik_bz))/2
         !if (itim_k==2) CYCLE ! time-reversal or not

         rhotwg1(:) = sym_rhotwgq0(itim_k,isym_k,dim_rtwg,npweps,rhxtwg_q0(:,iv,ic,ik_ibz,spin1),Gsph_c)
         !
         ! sum over G/=0
         ctemp = xdotc(npweps-1,rhotwg1(2:),1,rhotwg2(2:),1)
         ctemp = faq * kx_fact * ctemp

         ! exchange term is non divergent !
         if (BSp%prep_interp) then
           ccoeffs(ir) = ccoeffs(ir) + ctemp
         end if

         my_bsham(ir) = my_bsham(ir) + ctemp
       end do !it
     end do !itp

     ABI_FREE(rhotwg1)
     ABI_FREE(rhotwg2)
     ABI_FREE(vc_sqrt_qbz)

     call timab(683,2,tsec) ! exc_build_ham(exchange)
   end if ! do_exchange_term
   !
   ! =====================
   ! === Diagonal term ===
   ! =====================
   if (is_resonant .and. spin1==spin2) then
     write(msg,'(a,2i2,a)')" Adding diagonal term for (spin1,spin2) ",spin1,spin2," ..."
     call wrtout(std_out,msg,"COLL")
     do it=1,BSp%nreh(block)
       ir = it + it*(it-1_i8b)/2
       if (ir>=t_start(my_rank) .and. ir<=t_stop(my_rank)) my_bsham(ir) = my_bsham(ir) + Bsp%Trans(it,spin1)%en
     end do
   end if

   if (.FALSE.) then
     dump_unt = get_unit()
     msg=' Coupling Hamiltonian matrix elements: '
     if (is_resonant) msg=' Reasonant Hamiltonian matrix elements: '
     call wrtout(dump_unt,msg,"PERS")
     call wrtout(dump_unt,'    k  v  c  s      k" v" c" s"       H',"PERS")
     do itp=1,BSp%nreh(block)
       ikp_bz = Bsp%Trans(itp,spin2)%k
       ivp    = Bsp%Trans(itp,spin2)%v
       icp    = Bsp%Trans(itp,spin2)%c
       do it=1,itp
         ik_bz = Bsp%Trans(it,spin1)%k
         iv    = Bsp%Trans(it,spin1)%v
         ic    = Bsp%Trans(it,spin1)%c
         ir = it + itp*(itp-1_i8b)/2
         if (ir>=t_start(my_rank).and.ir<=t_stop(my_rank)) then
           http = my_bsham(ir)
           !if (ABS(http) > tol3) then
           write(msg,'(2(i0,1x),2(i5,3i3,3x),2f7.3)')it,itp, ik_bz,iv,ic,spin1, ikp_bz,ivp,icp,spin2, http
           call wrtout(dump_unt,msg,"PERS")
           !end if
         end if
       end do
     end do
   end if

!DBYG
   if (.False.) then
     dump_unt = get_unit()
     dump_unt = 999
     msg=' Coupling Hamiltonian matrix elements: '
     if (is_resonant) msg=' Resonant Hamiltonian matrix elements: '
     call wrtout(dump_unt,msg,"PERS")
     call wrtout(dump_unt,'    k v  c  s      k" v" c" s"       H',"PERS")
     do itp=1,BSp%nreh(block)
       ikp_bz = Bsp%Trans(itp,spin2)%k
       ivp    = Bsp%Trans(itp,spin2)%v
       icp    = Bsp%Trans(itp,spin2)%c
       do it=1,BSp%nreh(block)
         ik_bz = Bsp%Trans(it,spin1)%k
         iv    = Bsp%Trans(it,spin1)%v
         ic    = Bsp%Trans(it,spin1)%c
         if(it > itp) then
           ir = itp+it*(it-1_i8b)/2
         else
           ir = it + itp*(itp-1_i8b)/2
         end if
         if (ir>=t_start(my_rank).and.ir<=t_stop(my_rank)) then
           if(it > itp) then
             http = CONJG(my_bsham(ir))
             if (BSp%prep_interp) then
               aatmp = CONJG(acoeffs(ir))
               bbtmp = CONJG(bcoeffs(ir))
               cctmp = CONJG(ccoeffs(ir))
             end if
           else
             http = my_bsham(ir)
             if (BSp%prep_interp) then
               aatmp = acoeffs(ir)
               bbtmp = bcoeffs(ir)
               cctmp = ccoeffs(ir)
             end if
           end if
           if (it == itp) http = http - Bsp%Trans(it,spin1)%en
           !if (ABS(http) > tol3) then
           if (BSp%prep_interp) then
             write(msg,'(2(i0,1x),2(i5,3i3,3x),2f24.20,2f24.20,2f24.20,2f24.20)')it,itp, ik_bz,iv,ic,spin1, ikp_bz,ivp,icp,&
&   spin2, http, aatmp, bbtmp, cctmp
           else
             write(msg,'(2(i0,1x),2(i5,3i3,3x),2f24.20)')it,itp, ik_bz,iv,ic,spin1, ikp_bz,ivp,icp,spin2, http
           end if
           call wrtout(dump_unt,msg,"PERS")
           !end if
         end if
       end do
     end do
   end if

   call timab(684,1,tsec) ! exc_build_ham(synchro)
   call xmpi_barrier(comm)
   call timab(684,2,tsec) ! exc_build_ham(synchro)
   !
   ! =================================
   ! === Write Hamiltonian on disk ===
   ! =================================
   call timab(685,1,tsec) ! exc_build_ham(write_ham)
   if (use_mpiio) then
#ifdef HAVE_MPI_IO
     ! Write the Hamiltonian with collective MPI-IO.
     if (BSp%prep_interp) then
       MSG_ERROR("Preparation of interpolation technique not yet coded with MPI-IO")
     end if
     ABI_CHECK(nsppol==1,"nsppol==2 not coded, offset is wrong")
     !
     old_type = MPI_DOUBLE_COMPLEX
     call xmpio_create_fherm_packed(my_starts,my_ends,is_fortran_file,my_offset,old_type,hmat_type,offset_err)

     if (offset_err/=0) then
       write(msg,"(3a)")&
&        "Global position index cannot be stored in a standard Fortran integer. ",ch10,&
&        "BSE matrix cannot be written with a single MPI-IO call. "
       MSG_ERROR(msg)
     end if
     !
     ! Each node uses a different offset to skip the header and the blocks written by the other CPUs.
     my_offset = offset_of_block(block) + my_offset

     call MPI_FILE_SET_VIEW(mpi_fh, my_offset, MPI_BYTE, hmat_type, 'native', MPI_INFO_NULL, mpi_err)
     ABI_CHECK_MPI(mpi_err,"SET_VIEW")

     call MPI_TYPE_FREE(hmat_type,mpi_err)
     ABI_CHECK_MPI(mpi_err,"MPI_TYPE_FREE")

     if (hsize_of(my_rank) /= INT(hsize_of(my_rank),kind=i4b) ) then
       MSG_ERROR("Wraparound error")
     end if

     tmp_size = INT(hsize_of(my_rank))
     call MPI_FILE_WRITE_ALL(mpi_fh, my_bsham, tmp_size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpi_err)
     ABI_CHECK_MPI(mpi_err,"FILE_WRITE")

     ! It seems that personal calls in make the code stuck
     !if (is_fortran_file .and. my_rank==master) then ! Master writes the Fortran record markers.
     ! Write the Fortran record markers.
     neh2=BSp%nreh(block)
     ABI_MALLOC(bsize_frecord,(neh2))
     bsize_frecord = (/(col_glob * xmpi_bsize_dpc, col_glob=1,neh2)/)
     ! ehdr_offset points to the end of the header.
     !call xmpio_write_frmarkers(mpi_fh,ehdr_offset,xmpio_collective,neh2,bsize_frecord,mpi_err)
     my_offset = offset_of_block(block)
     call xmpio_write_frmarkers(mpi_fh,my_offset,xmpio_collective,neh2,bsize_frecord,ierr)
     ABI_CHECK(ierr==0,"Error while writing Fortran markers")
     ABI_FREE(bsize_frecord)
#else
     MSG_BUG("You should not be here!")
#endif
   else
     ! Use FORTRAN IO with sequential access mode.
     ! * Each node sends its data to master node.
     ! * Blocks are distributed according to the rank of the node.
     ! * Matrix is written by columns hence make sure that the last column is completely written.
     call cwtime(cputime,walltime,gflops,"start")

     if (my_rank==master) then
       prev_nrows=0; if (my_cols(2) /= my_rows(2)) prev_nrows = my_rows(2)
       ncol = my_cols(2)-my_cols(1)+1
       ist=1
       do jj=1,ncol
         col_glob = my_starts(2) + jj - 1
         nrows = col_glob; if (jj==ncol) nrows=my_rows(2)
         iend = ist + nrows -1
         write(bsh_unt) my_bsham(ist:iend)
         if (BSp%prep_interp) then
           write(a_unt) acoeffs(ist:iend)
           write(b_unt) bcoeffs(ist:iend)
           write(c_unt) ccoeffs(ist:iend)
         end if
         ist=iend+1
       end do
       write(msg,'(2(a,i0))')" Wraparound error: iend=",iend," my_hsize=",hsize_of(my_rank)
       ABI_CHECK(iend==hsize_of(my_rank),msg)
       ABI_FREE(my_bsham)
       if (BSp%prep_interp) then
         ABI_FREE(acoeffs)
         ABI_FREE(bcoeffs)
         ABI_FREE(ccoeffs)
       end if
     end if

     call xmpi_barrier(comm)
     !
     ! Collect data from the other nodes.
     do sender=1,nproc-1
       ! If I'm not involved, jump to the end of the loop and wait there (sequential IO? Of course!)
       if (all(my_rank /= [sender, master])) goto 100

       if (my_rank==master)  then
         ABI_MALLOC(buffer,(hsize_of(sender)))
         if (BSp%prep_interp) then
           ABI_MALLOC(abuffer,(hsize_of(sender)))
           ABI_MALLOC(bbuffer,(hsize_of(sender)))
           ABI_MALLOC(cbuffer,(hsize_of(sender)))
         end if
       end if
       tmp_size = INT(hsize_of(sender),kind=i4b)
       call xmpi_exch(my_bsham,tmp_size,sender,buffer,master,comm,mpi_err)
       if (BSp%prep_interp) then
         call xmpi_exch(acoeffs,tmp_size,sender,abuffer,master,comm,mpi_err)
         call xmpi_exch(bcoeffs,tmp_size,sender,bbuffer,master,comm,mpi_err)
         call xmpi_exch(ccoeffs,tmp_size,sender,cbuffer,master,comm,mpi_err)
       end if

       ! TODO Be careful with the MPI TAG here, add optional Arguments in xmpi_exch so that the TAG can be specified!
       proc_start = (/my_rows(1),my_cols(1)/)
       proc_end   = (/my_rows(2),my_cols(2)/)
       my_extrema(:,1) = proc_start
       my_extrema(:,2) = proc_end

       sender_extrema = my_extrema ! just to avoid NAN on sender. xechh_mpi is not well designed
       call xmpi_exch(my_extrema,4,sender,sender_extrema,master,comm,mpi_err)

       if (my_rank==master) then
          proc_start = sender_extrema(:,1)
          proc_end   = sender_extrema(:,2)
          !write(std_out,*)"proc_start, proc_end",proc_start,proc_end

         if (prev_nrows>0) then ! backspace the file if the last record written was not complete.
           !write(std_out,*)" master node had to call backspace"
           backspace(bsh_unt)
           ABI_MALLOC(prev_col,(prev_nrows))
           read(bsh_unt) prev_col
           backspace(bsh_unt)

           if (BSp%prep_interp) then
             backspace(a_unt)
             ABI_MALLOC(aprev_col,(prev_nrows))
             read(a_unt) aprev_col
             backspace(a_unt)

             backspace(b_unt)
             ABI_MALLOC(bprev_col,(prev_nrows))
             read(b_unt) bprev_col
             backspace(b_unt)

             backspace(c_unt)
             ABI_MALLOC(cprev_col,(prev_nrows))
             read(c_unt) cprev_col
             backspace(c_unt)
           end if
         end if
         !
         ! Write the columns owned by sender.
         ncol = proc_end(2)-proc_start(2)+1
         ist=1
         do jj=1,ncol
           col_glob = proc_start(2) + jj-1
           nrows = col_glob
           if (jj==1   )  nrows=col_glob - proc_start(1) + 1
           if (jj==ncol) then
             nrows=proc_end(1)
             if (ncol==1)  nrows=proc_end(1) - proc_start(1) + 1
           end if
           iend = ist + nrows -1
           !write(std_out,*)"Using nrows, ist, iend=",nrows,ist,iend
           if (jj==1 .and. prev_nrows>0) then ! join prev_col and this subcolumn.
             write(bsh_unt) CMPLX(prev_col,kind=dpc),CMPLX(buffer(ist:iend),kind=dpc)
             if (BSp%prep_interp) then
               write(a_unt) CMPLX(aprev_col,kind=dpc),CMPLX(abuffer(ist:iend),kind=dpc)
               write(b_unt) CMPLX(bprev_col,kind=dpc),CMPLX(bbuffer(ist:iend),kind=dpc)
               write(c_unt) CMPLX(cprev_col,kind=dpc),CMPLX(cbuffer(ist:iend),kind=dpc)
             end if
             prev_nrows = prev_nrows + iend-ist+1
           else
             write(bsh_unt) CMPLX(buffer(ist:iend),kind=dpc)
             if (BSp%prep_interp) then
               write(a_unt) CMPLX(abuffer(ist:iend),kind=dpc)
               write(b_unt) CMPLX(bbuffer(ist:iend),kind=dpc)
               write(c_unt) CMPLX(cbuffer(ist:iend),kind=dpc)
             end if
             prev_nrows=0
           end if
           ist=iend+1
         end do
         if (ncol>1) then ! Reset prev_nrows if a new column has begun.
           prev_nrows = proc_end(1)
           if (proc_end(1) == proc_end(2)) prev_nrows = 0
         end if
         if (iend/=hsize_of(sender)) then
           write(msg,'(2(a,i0))')" Wraparound error: iend=",iend," my_hsize=",hsize_of(sender)
           MSG_ERROR(msg)
         end if
         ABI_SFREE(prev_col)
         if (BSp%prep_interp) then
           ABI_SFREE(aprev_col)
           ABI_SFREE(bprev_col)
           ABI_SFREE(cprev_col)
         end if
         ABI_FREE(buffer)
         if (BSp%prep_interp) then
           ABI_FREE(abuffer)
           ABI_FREE(bbuffer)
           ABI_FREE(cbuffer)
         end if
       end if ! master
       !
100    call xmpi_barrier(comm)
     end do ! sender

     call cwtime(cputime,walltime,gflops,"stop")
     write(msg,'(2(a,f9.1),a)')" Fortran-IO completed. cpu_time: ",cputime,"[s], walltime: ",walltime," [s]"
     call wrtout(std_out,msg,"COLL",do_flush=.True.)
   end if ! use_mpiio
   call timab(685,2,tsec) ! exc_build_ham(write_ham)
   !
   ABI_SFREE(my_bsham)
   if (BSp%prep_interp) then
     ABI_SFREE(acoeffs)
     ABI_SFREE(bcoeffs)
     ABI_SFREE(ccoeffs)
   end if
   ABI_FREE(t_start)
   ABI_FREE(t_stop)
   ABI_FREE(hsize_of)
 end do ! block
 !
 ! ===========================================
 ! === Exchange term for spin_up spin_down ===
 ! ===========================================

 if (nsppol==2) then
   call timab(686,2,tsec) ! exc_build_ham(exch.spin)
   block=3
   neh1=BSp%nreh(1)
   neh2=BSp%nreh(2)
   !
   ! The oscillators at q=0 are available on each node for both spin.
   ! Here the calculation of the block is parallelized over columns.
   ABI_MALLOC(col_start,(0:nproc-1))
   ABI_MALLOC(col_stop,(0:nproc-1))
   call xmpi_split_work2_i4b(neh2,nproc,col_start,col_stop)

   my_cols(1) = col_start(my_rank)
   my_cols(2) = col_stop (my_rank)
   if (my_cols(2)-my_cols(1)<=0) then
     MSG_ERROR("One of the processors has zero columns!")
   end if

   ABI_MALLOC(ncols_of,(0:nproc-1))
   ncols_of=0
   do rank=0,nproc-1
     if (col_stop(rank)>=col_start(rank)) ncols_of(rank) = col_stop(rank)-col_start(rank)+1
   end do

   ABI_FREE(col_start)
   ABI_FREE(col_stop)
   !
   ! TODO might add treatment of <psi|q+G|psi> for q+G -> 0
   ! TODO might used enlarged G-sphere for better convergence.
   ! Note that my_kxssp is always written on file when nsppol=2, even when
   ! non-local field effects are neglected.
   ABI_MALLOC(my_kxssp,(neh1,my_cols(1):my_cols(2)))
   my_kxssp=czero

   if (do_exchange_term) then
     spin1=1; spin2=2
     write(msg,'(a,2i2,a)')" Calculating exchange term for (spin1,spin2) ",spin1,spin2," ..."
     call wrtout(std_out,msg,"COLL")

     ABI_MALLOC(rhotwg1,(npweps))
     ABI_MALLOC(rhotwg2,(npweps))

     ngx = Gsph_x%ng
     ABI_MALLOC(vc_sqrt_qbz,(ngx))
     !
     ! * Get iq_ibz, and symmetries from iq_bz.
     iq_bz = iqbz0 ! q = 0 -> iqbz0
     call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)
     !
     ! * Set up table of |q(BZ)+G|
     if (iq_ibz==1) then
       do ig=1,ngx
         ISg = Gsph_x%rottb(ig,itim_q,isym_q)
         vc_sqrt_qbz(ISg)=Vcp%vcqlwl_sqrt(ig,1)
       end do
     else
        MSG_ERROR("iq_ibz should be 1")
     end if

     do itp=1,neh2 ! Loop over transition tp = (kp,vp,cp,spin2)

       if (itp<my_cols(1) .or. itp>my_cols(2)) CYCLE ! I dont have this column.
       ene_tp = Bsp%Trans(itp,spin2)%en
       ikp_bz = Bsp%Trans(itp,spin2)%k
       ivp    = Bsp%Trans(itp,spin2)%v
       icp    = Bsp%Trans(itp,spin2)%c

       ikp_ibz = Kmesh%tab (ikp_bz)
       isym_kp = Kmesh%tabo(ikp_bz)
       itim_kp = (3-Kmesh%tabi(ikp_bz))/2

       if (is_resonant) then
         rhotwg2(:) = sym_rhotwgq0(itim_kp,isym_kp,dim_rtwg,npweps,rhxtwg_q0(:,ivp,icp,ikp_ibz,spin2),Gsph_c)
       else ! Code for coupling block.
         rhotwg2(:) = sym_rhotwgq0(itim_kp,isym_kp,dim_rtwg,npweps,rhxtwg_q0(:,icp,ivp,ikp_ibz,spin2),Gsph_c)
       end if
       !
       ! Multiply by the Coulomb term.
        do ig=2,npweps
          rhotwg2(ig) = rhotwg2(ig) * vc_sqrt_qbz(ig) * vc_sqrt_qbz(ig)
        end do

       do it=1,neh1 ! Loop over transition t = (k,v,c,spin1) FULL matrix.

         ene_t = Bsp%Trans(it,spin1)%en
         ik_bz = Bsp%Trans(it,spin1)%k
         iv    = Bsp%Trans(it,spin1)%v
         ic    = Bsp%Trans(it,spin1)%c

         ik_ibz = Kmesh%tab(ik_bz)
         isym_k = Kmesh%tabo(ik_bz)
         itim_k = (3-Kmesh%tabi(ik_bz))/2
         !if (itim_k==2) CYCLE ! time-reversal or not

         rhotwg1(:) = sym_rhotwgq0(itim_k,isym_k,dim_rtwg,npweps,rhxtwg_q0(:,iv,ic,ik_ibz,spin1),Gsph_c)
         !
         ! sum over G/=0
         ctemp = XDOTC(npweps-1,rhotwg1(2:),1,rhotwg2(2:),1)
         ctemp = faq * kx_fact * ctemp

         my_kxssp(it,itp) = ctemp
       end do !it
     end do !itp

     ABI_FREE(rhotwg1)
     ABI_FREE(rhotwg2)
     ABI_FREE(vc_sqrt_qbz)
   end if ! do_exchange_term
   call timab(686,2,tsec) ! exc_build_ham(exch.spin)
   !
   ! =====================================
   ! === Write the Hamiltonian on disk ===
   ! =====================================
   call timab(685,1,tsec) ! exc_build_ham(write_ham)

   if (use_mpiio) then
#ifdef HAVE_MPI_IO
     my_ncols=ncols_of(my_rank); old_type=MPI_DOUBLE_COMPLEX
     call xmpio_create_fsubarray_2D((/neh1,my_ncols/),(/neh1,my_ncols/),(/1,1/),old_type,hmat_type,my_offpad,mpi_err)
     ABI_CHECK_MPI(mpi_err,"fsubarray_2D")
     !
     ! Each node uses a different offset to skip the header and the blocks written by the other CPUs.
     prev_nels=0
     prev_ncols=0
     if (my_rank>0) then
       prev_ncols = SUM(ncols_of(0:my_rank-1))
       prev_nels = neh1*prev_ncols
     end if
     tmp_off = prev_nels*xmpi_bsize_dpc + prev_ncols*2*xmpio_bsize_frm

     my_offset = offset_of_block(block) + tmp_off + my_offpad

     call MPI_FILE_SET_VIEW(mpi_fh, my_offset, MPI_BYTE, hmat_type, 'native', MPI_INFO_NULL, mpi_err)
     ABI_CHECK_MPI(mpi_err,"SET_VIEW")

     call MPI_TYPE_FREE(hmat_type,mpi_err)
     ABI_CHECK_MPI(mpi_err,"MPI_TYPE_FREE")

     tmp_size = INT(neh1*my_ncols)
     call MPI_FILE_WRITE_ALL(mpi_fh, my_kxssp,tmp_size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpi_err)
     ABI_CHECK_MPI(mpi_err,"FILE_WRITE")

     ! It seems that personal calls in make the code stuck
     ! Master writes the Fortran record markers.
     ABI_MALLOC(bsize_frecord,(neh2))
     bsize_frecord = neh1 * xmpi_bsize_dpc
     ! ehdr_offset points to the end of the header.
     !call xmpio_write_frmarkers(mpi_fh,ehdr_offset,xmpio_collective,neh2,bsize_frecord,mpi_err)
     my_offset = offset_of_block(block)
     call xmpio_write_frmarkers(mpi_fh,my_offset,xmpio_collective,neh2,bsize_frecord,ierr)
     ABI_CHECK(ierr==0,"Error while writing Fortran markers")
     ABI_FREE(bsize_frecord)
#else
     MSG_BUG("You should not be here")
#endif
   else
     ! Use FORTRAN IO with sequential access mode.
     ! * Each node sends its data to master node.
     ! * Columns are distributed according to the rank of the node.
     if (my_rank==master) then
       do jj=my_cols(1),my_cols(2)
         write(bsh_unt) my_kxssp(:,jj)
       end do
       ABI_FREE(my_kxssp)
     end if

     call xmpi_barrier(comm)
     !
     ! Collect data from the other nodes.
     do sender=1,nproc-1
       ! If I'm not involved, jump to the end of the loop and wait there (sequential IO? Of course!)
       if (all(my_rank /= [sender, master])) goto 200

       if (my_rank==master)  then
         ABI_MALLOC(buffer_2d,(neh1,ncols_of(sender)))
       end if
       call xmpi_exch(my_kxssp,neh1*ncols_of(sender),sender,buffer_2d,master,comm,mpi_err)
       !
       if (my_rank==master) then ! Write the columns owned by sender.
         do jj=1,ncols_of(sender)
           write(bsh_unt) buffer_2d(:,jj)
         end do
         ABI_FREE(buffer_2d)
       end if ! master
       !
200    call xmpi_barrier(comm)
     end do ! sender
   end if
   call timab(685,2,tsec) ! exc_build_ham(write_ham)

   ABI_FREE(ncols_of)
   ABI_SFREE(my_kxssp)
 end if

 ! Close the file.
 if (use_mpiio) then
#ifdef HAVE_MPI_IO
   call MPI_FILE_CLOSE(mpi_fh, mpi_err)
   ABI_CHECK_MPI(mpi_err,"FILE_CLOSE")
   ABI_FREE(offset_of_block)
#endif
 end if

 ! master closes the Fortran files.
 if (my_rank==master) then
   close(bsh_unt)
   if (BSp%prep_interp) then
     close(a_unt)
     close(b_unt)
     close(c_unt)
   end if
 end if

 ! Free memory.
 ABI_FREE(igfftg0)
 ABI_FREE(ktabr_k)
 ABI_FREE(id_tab)
 ABI_FREE(ktabr_kp)
 ABI_FREE(rhxtwg_vpv)
 ABI_FREE(rhxtwg_cpc)
 if (BSp%prep_interp) then
   ABI_FREE(aa_vpv)
   ABI_FREE(bb_vpv1)
   ABI_FREE(bb_vpv2)
   ABI_FREE(cc_vpv)
   ABI_FREE(aa_cpc)
   ABI_FREE(bb_cpc1)
   ABI_FREE(bb_cpc2)
   ABI_FREE(cc_cpc)
 end if
 ABI_FREE(ur_ckp)
 ABI_FREE(ur_vkp)
 ABI_FREE(ur_vk)
 ABI_FREE(ur_ck)

 ! Deallocation for PAW.
 if (Wfd%usepaw==1) then
   call pawcprj_free(Cp_vk)
   ABI_DT_FREE(Cp_vk)
   call pawcprj_free(Cp_ck)
   ABI_DT_FREE(Cp_ck)
   call pawcprj_free(Cp_ckp)
   ABI_DT_FREE(Cp_ckp)
   call pawcprj_free(Cp_vkp)
   ABI_DT_FREE(Cp_vkp)
   call pawcprj_free(Cp_tmp1)
   ABI_DT_FREE(Cp_tmp1)
   call pawcprj_free(Cp_tmp2)
   ABI_DT_FREE(Cp_tmp2)
   call pawcprj_free(Cp_tmp3)
   ABI_DT_FREE(Cp_tmp3)
   call pawcprj_free(Cp_tmp4)
   ABI_DT_FREE(Cp_tmp4)
 end if

 call xmpi_barrier(comm)

 DBG_EXIT("COLL")

 call timab(680,2,tsec)

end subroutine exc_build_block
!!***

!!****f* m_exc_build/exc_build_v
!! NAME
!!  exc_build_v
!!
!! FUNCTION
!!  Calculate and write the excitonic Hamiltonian on an external binary file (Fortran file open
!!  in random mode) for subsequent treatment in the Bethe-Salpeter code.
!!
!! INPUTS
!!  BSp<excparam>=The parameters for the Bethe-Salpeter calculation.
!!  Cryst<crystal_t>=Info on the crystalline structure.
!!  Kmesh<kmesh_t>=The list of k-points in the BZ, IBZ and symmetry tables.
!!  Qmesh<kmesh_t>=The list of q-points for epsilon^{-1} and related symmetry tables.
!!  Gsph_x<gsphere_t>=Info on the G-sphere used to describe wavefunctions and W (the largest one is actually stored).
!!  Gsph_c<gsphere_t>=Info on the G-sphere used to describe the correlation part.
!!  Vcp<vcoul_t>=The Coulomb interaction in reciprocal space. A cutoff can be used
!!  rhxtwg_q0
!!  is_resonant
!!  comm=MPI communicator.
!!
!! OUTPUT
!!
!! NOTES
!!  *) Version for K_V = K_C (q=0), thus KP_V = KP_C
!!  *) No exchange limit: use LDA energies in case.
!!  *) Symmetry of H(-k-k') = H*(k k') not used.
!!  *) Coulomb term can be approssimateed as diagonal in G.
!!  *) Valence bands treated from lomo on.
!!  *) Symmetries of the sub-blocks are used to reduce the number of elements to calculate.
!!
!!            ____________
!!           |_(cv)__(vc)_|
!!   H_exc = |  R      C  |
!!           | -C*    -R* |
!!
!!   where C is symmetric and R is Hermitian provided that the QP energies are real.
!!
!!  For nsppol=1 ==> R = diag-W+2v; C = -W+2v
!!  since the Hamiltonian can be diagonalized in the spin-singlet basis set thanks to
!!  the fact that spin triplet does not contribute to the optical limit of epsilon.
!!
!!  For nsppol=2 ==> R = diag-W+v; C = -W+v
!!  Now the matrix elements depend on the spin of the transitions but only those
!!  transitions in which the spin of the electron and of the hole are equal contribute
!!  to the macroscopic dielectric function. Moreover only the exchange term can connect transitions of different spin.
!!  When nsppol==2 the transitions are ordered using | (cv up) | (cv dwn) | (vc up) | (vc down) |
!!
!!  The resonant block is given by:
!!
!!      |  (v'c' up)       | (v'c' dwn)   |
!!      -----------------------------------           where v_{-+} = v_{+-}^H when the momentum of the photon is neglected.
!!      | [diag-W+v]++     |      v+-     | (vc up)   Note that v_{+-} is not Hermitian due to the presence of different spins.
!!  R = -----------------------------------           Actually it reduces to a Hermitian matrix when the system is not spin polarized.
!!      |     v-+          | [diag-W+v]-- | (vc dwn)  but in this case one should use nsppol=1.
!!      -----------------------------------           As a consequence the entire matrix is calculated and stored on file.
!!
!!  The coupling block is given by:
!!
!!      |  (c'v' up)   |    (c'v dwn)     |
!!      -----------------------------------           where v_{-+} = v_{+-}^t when the momentum of the photon is neglected.
!!      | [-W+v]++     |      v+-         | (vc up)   Also in this case the entire matrix v_{+-} has to be calculated
!!  C = -----------------------------------           and stored on file.
!!      |     v-+      |    [-W+v]--      | (vc dwn)
!!      -----------------------------------
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine exc_build_v(spin1,spin2,nsppol,npweps,Bsp,Cryst,Kmesh,Qmesh,Gsph_x,Gsph_c,Vcp,&
&  is_resonant,rhxtwg_q0,nproc,my_rank,t_start,t_stop,my_bsham)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin1,spin2,nsppol,npweps,nproc,my_rank
 logical,intent(in) :: is_resonant
 type(excparam),intent(in) :: BSp
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(crystal_t),intent(in) :: Cryst
 type(vcoul_t),intent(in) :: Vcp
 type(gsphere_t),intent(in) :: Gsph_x,Gsph_c
!arrays
 integer(i8b),intent(in) :: t_start(0:nproc-1),t_stop(0:nproc-1)
 complex(gwpc),intent(in) :: rhxtwg_q0(npweps,BSp%lomo_min:BSp%humo_max,BSp%lomo_min:BSp%humo_max,Kmesh%nibz,nsppol)
 complex(dpc),intent(inout) :: my_bsham(t_start(my_rank):t_stop(my_rank))

!Local variables ------------------------------
!scalars
 integer :: ISg,ngx,ik_bz,ikp_bz,dim_rtwg
 integer :: neh1,neh2,ig,nblocks
 integer :: ik_ibz,itim_k,ikp_ibz,itim_kp,isym_k,isym_kp
 integer :: iq_bz,iq_ibz,isym_q,itim_q,iqbz0,rank
 integer :: iv,ivp,ic,icp
 integer :: block
 integer(i8b) :: tot_nels,ir,it,itp
 real(dp) :: faq,kx_fact
 complex(spc) :: ctemp
 character(len=500) :: msg
!arrays
 integer :: bidx(2,4),spin_ids(2,3)
 integer(i8b) :: nels_block(3)
 integer :: my_cols(2),my_rows(2) !,proc_end(2),proc_start(2)
 integer,allocatable :: ncols_of(:)
 integer,allocatable :: col_start(:),col_stop(:)
 real(dp) :: qbz(3),tsec(2) !kbz(3),kpbz(3),
 complex(dpc),allocatable :: my_kxssp(:,:)
 complex(gwpc),allocatable :: vc_sqrt_qbz(:),rhotwg1(:),rhotwg2(:)

!************************************************************************

 DBG_ENTER("COLL")

 write(msg,'(a,2i2,a)')" Calculating exchange term for (spin1,spin2) ",spin1,spin2," ..."
 call wrtout(std_out,msg,"COLL")

 ! Basic constants.
 dim_rtwg=1; faq = one/(Cryst%ucvol*Kmesh%nbz)

 ! Identify the index of q==0
 iqbz0=0
 do iq_bz=1,Qmesh%nbz
   if (ALL(ABS(Qmesh%bz(:,iq_bz))<tol3)) iqbz0 = iq_bz
 end do
 ABI_CHECK(iqbz0/=0,"q=0 not found")
 !
 ! Treat the spin polarization.
 spin_ids(:,1) = (/1,1/)
 spin_ids(:,2) = (/2,2/)
 spin_ids(:,3) = (/1,2/)

 nblocks=1
 kx_fact=two
 nels_block(:)=0
 nels_block(1)=BSp%nreh(1)*(BSp%nreh(1)+1_i8b)/2
 tot_nels=nels_block(1)

 if (nsppol==2) then
   nblocks=3
   kx_fact=one
   nels_block(1) = BSp%nreh(1)*(BSp%nreh(1)+1_i8b)/2   ! Only the upper triangle for block 1 and 2
   nels_block(2) = BSp%nreh(2)*(BSp%nreh(2)+1_i8b)/2
   nels_block(3) = BSp%nreh(1)*BSp%nreh(2)*1_i8b       ! Note: Block 3 does not have symmetries.
   tot_nels= SUM(nels_block)
 end if
 !
 ! Distribute the calculation of the matrix elements among the nodes.
 ! * tstart and t_stop give the initial and final transition index treated by each node.
 ! * my_hsize is the number of transitions treated by this processor
 ! * my_cols(1:2) gives the initial and final column treated by this node.
 !
 do block=1,nsppol
   !
   ! Indices used to loop over bands.
   ! bidx contains the starting and final indices used to loop over bands.
   !
   !      (b3,b4)
   !         |... ...|
   ! (b1,b2) |... ...|
   !
   ! Resonant matrix is given by
   !      (v',c')
   !       |... ...|
   ! (v,c) |... ...|
   !
   ! Coupling matrix is given by
   !       (c',v')
   !       |... ...|
   ! (v,c) |... ...|

   if (is_resonant) then
     bidx(:,1) = [BSp%lomo_spin(block),BSp%homo_spin(block)] ! range for b1
     bidx(:,2) = [BSp%lumo_spin(block),BSp%humo_spin(block)] ! range for b2
     bidx(:,3) = [BSp%lomo_spin(block),BSp%homo_spin(block)] ! range for b3
     bidx(:,4) = [BSp%lumo_spin(block),BSp%humo_spin(block)] ! range for b4
   else
     bidx(:,1) = [BSp%lomo_spin(block),BSp%homo_spin(block)] ! range for b1
     bidx(:,2) = [BSp%lumo_spin(block),BSp%humo_spin(block)] ! range for b2
     bidx(:,3) = [BSp%lumo_spin(block),BSp%humo_spin(block)] ! range for b3
     bidx(:,4) = [BSp%lomo_spin(block),BSp%homo_spin(block)] ! range for b4
   end if

   !spin1 = spin_ids(1,block)
   !spin2 = spin_ids(2,block)

   my_cols=0
   do itp=1,Bsp%nreh(block)
     do it=1,itp
       ir = it + itp*(itp-1_i8b)/2
       if (ir==t_start(my_rank)) then
         my_rows(1) = it
         my_cols(1) = itp
       end if
       if (ir==t_stop(my_rank)) then
         my_rows(2) = it
         my_cols(2) = itp
       end if
     end do
   end do

   ! Allocate big (scalable) buffer to store the BS matrix on this node.
   !ABI_MALLOC(my_bsham,(t_start(my_rank):t_stop(my_rank)))
   !
   ! =====================
   ! === Exchange term ===
   ! =====================
   ! TODO might add treatment of <psi|q+G|psi> for q+G -> 0
   ! TODO might used enlarged G-sphere for better convergence.
!if (do_exchange_term) then
   call timab(683,1,tsec) ! exc_build_ham(exchange)

   ABI_MALLOC(rhotwg1,(npweps))
   ABI_MALLOC(rhotwg2,(npweps))

   ngx = Gsph_x%ng
   ABI_MALLOC(vc_sqrt_qbz,(ngx))

   ! * Get iq_ibz, and symmetries from iq_bz.
   iq_bz = iqbz0 ! q = 0 -> iqbz0
   call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)

   ! * Set up table of |q(BZ)+G|
   if (iq_ibz==1) then
     do ig=1,ngx
       ISg = Gsph_x%rottb(ig,itim_q,isym_q)
       vc_sqrt_qbz(ISg)=Vcp%vcqlwl_sqrt(ig,1)
     end do
   else
      MSG_ERROR("iq_ibz should be 1")
   end if

   do itp=1,BSp%nreh(block) ! Loop over transition tp = (kp,vp,cp,spin2)

     if (itp<my_cols(1) .or. itp>my_cols(2)) CYCLE ! I dont have this column.
     ikp_bz = Bsp%Trans(itp,spin2)%k
     ivp    = Bsp%Trans(itp,spin2)%v
     icp    = Bsp%Trans(itp,spin2)%c

     ikp_ibz = Kmesh%tab (ikp_bz)
     isym_kp = Kmesh%tabo(ikp_bz)
     itim_kp = (3-Kmesh%tabi(ikp_bz))/2

     if (is_resonant) then
       rhotwg2(:) = sym_rhotwgq0(itim_kp,isym_kp,dim_rtwg,npweps,rhxtwg_q0(:,ivp,icp,ikp_ibz,spin2),Gsph_c)
     else ! Code for coupling block.
       rhotwg2(:) = sym_rhotwgq0(itim_kp,isym_kp,dim_rtwg,npweps,rhxtwg_q0(:,icp,ivp,ikp_ibz,spin2),Gsph_c)
     end if
     !
     ! Multiply by the Coulomb term.
      do ig=2,npweps
        rhotwg2(ig) = rhotwg2(ig) * vc_sqrt_qbz(ig) * vc_sqrt_qbz(ig)
      end do

     do it=1,itp ! Loop over transition t = (k,v,c,spin1)
       ir = it + itp*(itp-1_i8b)/2
       if (ir<t_start(my_rank) .or. ir>t_stop(my_rank)) CYCLE

       ik_bz   = Bsp%Trans(it,spin1)%k
       iv      = Bsp%Trans(it,spin1)%v
       ic      = Bsp%Trans(it,spin1)%c

       ik_ibz = Kmesh%tab(ik_bz)
       isym_k = Kmesh%tabo(ik_bz)
       itim_k = (3-Kmesh%tabi(ik_bz))/2
       !if (itim_k==2) CYCLE ! time-reversal or not

       rhotwg1(:) = sym_rhotwgq0(itim_k,isym_k,dim_rtwg,npweps,rhxtwg_q0(:,iv,ic,ik_ibz,spin1),Gsph_c)
       !
       ! sum over G/=0
       ctemp = xdotc(npweps-1,rhotwg1(2:),1,rhotwg2(2:),1)
       ctemp = faq * kx_fact * ctemp

       ! exchange term is non divergent !
       !if (BSp%prep_interp) then
       !  ccoeffs(ir) = ccoeffs(ir) + ctemp
       !end if

       my_bsham(ir) = my_bsham(ir) + ctemp
     end do !it
   end do !itp

   ABI_FREE(rhotwg1)
   ABI_FREE(rhotwg2)
   ABI_FREE(vc_sqrt_qbz)

   call timab(683,2,tsec) ! exc_build_ham(exchange)
!end if ! do_exchange_term
 end do ! block

 !
 ! ===========================================
 ! === Exchange term for spin_up spin_down ===
 ! ===========================================

if (nsppol==2) then
 call timab(686,2,tsec) ! exc_build_ham(exch.spin)
 block=3
 neh1=BSp%nreh(1)
 neh2=BSp%nreh(2)
 !
 ! The oscillators at q=0 are available on each node for both spin.
 ! Here the calculation of the block is parallelized over columns.
 ABI_MALLOC(col_start,(0:nproc-1))
 ABI_MALLOC(col_stop,(0:nproc-1))
 call xmpi_split_work2_i4b(neh2,nproc,col_start,col_stop) !check this but it should be OK.

 my_cols(1) = col_start(my_rank)
 my_cols(2) = col_stop (my_rank)
 if (my_cols(2)-my_cols(1)<=0) then
   MSG_ERROR("One of the processors has zero columns!")
 end if

 ABI_MALLOC(ncols_of,(0:nproc-1))
 ncols_of=0
 do rank=0,nproc-1
   if (col_stop(rank)>=col_start(rank)) ncols_of(rank) = col_stop(rank)-col_start(rank)+1
 end do

 ABI_FREE(col_start)
 ABI_FREE(col_stop)
 !
 ! TODO might add treatment of <psi|q+G|psi> for q+G -> 0
 ! TODO might used enlarged G-sphere for better convergence.
 ! Note that my_kxssp is always written on file when nsppol=2, even when
 ! non-local field effects are neglected.
 ABI_MALLOC(my_kxssp,(neh1,my_cols(1):my_cols(2)))
 my_kxssp=czero

 !if (do_exchange_term) then
   !spin1=1; spin2=2
   ABI_MALLOC(rhotwg1,(npweps))
   ABI_MALLOC(rhotwg2,(npweps))

   ngx = Gsph_x%ng
   ABI_MALLOC(vc_sqrt_qbz,(ngx))
   !
   ! * Get iq_ibz, and symmetries from iq_bz.
   iq_bz = iqbz0 ! q = 0 -> iqbz0
   call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)
   !
   ! * Set up table of |q(BZ)+G|
   if (iq_ibz==1) then
     do ig=1,ngx
       ISg = Gsph_x%rottb(ig,itim_q,isym_q)
       vc_sqrt_qbz(ISg)=Vcp%vcqlwl_sqrt(ig,1)
     end do
   else
      MSG_ERROR("iq_ibz should be 1")
   end if

   do itp=1,neh2 ! Loop over transition tp = (kp,vp,cp,spin2)

     if (itp<my_cols(1) .or. itp>my_cols(2)) CYCLE ! I dont have this column.
     ikp_bz = Bsp%Trans(itp,spin2)%k
     ivp    = Bsp%Trans(itp,spin2)%v
     icp    = Bsp%Trans(itp,spin2)%c

     ikp_ibz = Kmesh%tab (ikp_bz)
     isym_kp = Kmesh%tabo(ikp_bz)
     itim_kp = (3-Kmesh%tabi(ikp_bz))/2

     if (is_resonant) then
       rhotwg2(:) = sym_rhotwgq0(itim_kp,isym_kp,dim_rtwg,npweps,rhxtwg_q0(:,ivp,icp,ikp_ibz,spin2),Gsph_c)
     else ! Code for coupling block.
       rhotwg2(:) = sym_rhotwgq0(itim_kp,isym_kp,dim_rtwg,npweps,rhxtwg_q0(:,icp,ivp,ikp_ibz,spin2),Gsph_c)
     end if
     !
     ! Multiply by the Coulomb term.
      do ig=2,npweps
        rhotwg2(ig) = rhotwg2(ig) * vc_sqrt_qbz(ig) * vc_sqrt_qbz(ig)
      end do

     do it=1,neh1 ! Loop over transition t = (k,v,c,spin1) FULL matrix.
       ik_bz = Bsp%Trans(it,spin1)%k
       iv    = Bsp%Trans(it,spin1)%v
       ic    = Bsp%Trans(it,spin1)%c

       ik_ibz = Kmesh%tab(ik_bz)
       isym_k = Kmesh%tabo(ik_bz)
       itim_k = (3-Kmesh%tabi(ik_bz))/2
       !if (itim_k==2) CYCLE ! time-reversal or not

       rhotwg1(:) = sym_rhotwgq0(itim_k,isym_k,dim_rtwg,npweps,rhxtwg_q0(:,iv,ic,ik_ibz,spin1),Gsph_c)
       !
       ! sum over G/=0
       ctemp = XDOTC(npweps-1,rhotwg1(2:),1,rhotwg2(2:),1)
       ctemp = faq * kx_fact * ctemp

       my_kxssp(it,itp) = ctemp
     end do !it
   end do !itp

   ABI_FREE(rhotwg1)
   ABI_FREE(rhotwg2)
   ABI_FREE(vc_sqrt_qbz)
 !end if ! do_exchange_term
 call timab(686,2,tsec) ! exc_build_ham(exch.spin)

 ABI_FREE(ncols_of)
 ABI_SFREE(my_kxssp)
 end if

 DBG_EXIT("COLL")

end subroutine exc_build_v
!!***

!!****f* m_exc_build/exc_build_ham
!! NAME
!!  exc_build_ham
!!
!! FUNCTION
!!  Calculate and write the excitonic Hamiltonian on an external binary file (Fortran file open
!!  in random mode) for subsequent treatment in the Bethe-Salpeter code.
!!
!! INPUTS
!!  BSp<excparam>=The parameters for the Bethe-Salpeter calculation.
!!  BS_files<excfiles>=File names internally used in the BS code.
!!  Cryst<crystal_t>=Info on the crystalline structure.
!!  Kmesh<kmesh_t>=The list of k-points in the BZ, IBZ and symmetry tables.
!!  Qmesh<kmesh_t>=The list of q-points for epsilon^{-1} and related symmetry tables.
!!  ktabr(nfftot_osc,BSp%nkbz)=The FFT index of $(R^{-1}(r-\tau))$ where R is symmetry needed to obtains
!!    the k-points from the irreducible image.  Used to symmetrize u_Sk where S = \transpose R^{-1}
!!  Gsph_x<gsphere_t>=Info on the G-sphere used to describe wavefunctions and W (the largest one is actually stored).
!!  Gsph_c<gsphere_t>=Info on the G-sphere used to describe the correlation part.
!!  Vcp<vcoul_t>=The Coulomb interaction in reciprocal space. A cutoff can be used
!!  W<screen_t>=Data type gathering info and data for W.
!!  nfftot_osc=Total Number of FFT points used for the oscillator matrix elements.
!!  ngfft_osc(18)=Info on the FFT algorithm used to calculate the oscillator matrix elements.
!!  Psps<Pseudopotential_type>=Variables related to pseudopotentials
!!  Pawtab(Psps%ntypat)<pawtab_type>=PAW tabulated starting data.
!!  Pawang<pawang_type>=PAW angular mesh and related data.
!!  Paw_pwff(Cryst%ntypat*Wfd%usepaw)<pawpwff_t>=Form factor used to calculate the onsite mat. elements of a plane wave.
!!  Wfd<wfd_t>=Handler for the wavefunctions.
!!
!! OUTPUT
!!  The excitonic Hamiltonian is saved on an external binary file (see below).
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!
!! SOURCE

subroutine exc_build_ham(BSp,BS_files,Cryst,Kmesh,Qmesh,ktabr,Gsph_x,Gsph_c,Vcp,&
& Wfd,W,Hdr_bse,nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot_osc
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(screen_t),intent(inout) :: W
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(crystal_t),intent(in) :: Cryst
 type(vcoul_t),intent(in) :: Vcp
 type(gsphere_t),intent(in) :: Gsph_x,Gsph_c
 type(Pseudopotential_type),intent(in) :: Psps
 type(Hdr_type),intent(inout) :: Hdr_bse
 type(pawang_type),intent(in) :: Pawang
 type(wfd_t),target,intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfft_osc(18)
 integer,intent(in) :: ktabr(nfftot_osc,Kmesh%nbz)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Wfd%usepaw)
 type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 logical :: do_resonant,do_coupling
 !character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 complex(gwpc),allocatable :: all_mgq0(:,:,:,:,:)

!************************************************************************

 call timab(670,1,tsec)

 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")
 ABI_CHECK(nfftot_osc==PRODUCT(ngfft_osc(1:3)),"mismatch in FFT size")

 if (BSp%have_complex_ene) then
   MSG_ERROR("Complex energies are not supported yet")
 end if

 ! Do we have to compute some block?
 do_resonant = (BS_files%in_hreso == BSE_NOFILE)
 do_coupling = (BS_files%in_hcoup == BSE_NOFILE)

 if (BSp%use_coupling == 0) then
   if (.not.do_resonant) then
     call wrtout(std_out,"Will skip the calculation of resonant block (will use BSR file)","COLL")
     goto 100
   end if
 else
   if (.not. do_resonant .and. .not. do_coupling) then
     call wrtout(std_out,"Will skip the calculation of both resonant and coupling block (will use BSR and BSC files)","COLL")
     goto 100
   end if
 end if

 ! Compute M_{k,q=0}^{b,b}(G) for all k-points in the IBZ and each pair b, b'
 ! used for the exchange part and part of the Coulomb term.
 call wrtout(std_out," Calculating all matrix elements for q=0 to save CPU time","COLL")

 call wfd_all_mgq0(Wfd,Cryst,Qmesh,Gsph_x,Vcp,Psps,Pawtab,Paw_pwff,&
&  Bsp%lomo_spin,Bsp%homo_spin,Bsp%humo_spin,nfftot_osc,ngfft_osc,Bsp%npweps,all_mgq0)

 ! ========================
 ! ==== Resonant Block ====
 ! ========================
 if (do_resonant) then
   call timab(672,1,tsec)
   call exc_build_block(BSp,Cryst,Kmesh,Qmesh,ktabr,Gsph_x,Gsph_c,Vcp,&
&    Wfd,W,Hdr_bse,nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff,all_mgq0,.TRUE.,BS_files%out_hreso)
   call timab(672,2,tsec)
 end if

 ! ========================
 ! ==== Coupling Block ====
 ! ========================
 if (do_coupling.and.BSp%use_coupling>0) then
   call timab(673,1,tsec)
   call exc_build_block(BSp,Cryst,Kmesh,Qmesh,ktabr,Gsph_x,Gsph_c,Vcp,&
&    Wfd,W,Hdr_bse,nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff,all_mgq0,.FALSE.,BS_files%out_hcoup)
   call timab(673,2,tsec)
 end if

 ! Free memory.
 ABI_FREE(all_mgq0)

100 call timab(670,2,tsec)

end subroutine exc_build_ham
!!***

!!****f* m_exc_build/wfd_all_mgq0
!! NAME
!!  wfd_all_mgq0
!!
!! FUNCTION
!!
!! INPUTS
!!  Wfd<wfd_t>=Handler for the wavefunctions.
!!  Cryst<crystal_t>=Info on the crystalline structure.
!!  Qmesh<kmesh_t>=The list of q-points for epsilon^{-1} and related symmetry tables.
!!  Gsph_x<gsphere_t>=G-sphere with the G-vectors in mgq0.
!!  Vcp<vcoul_t>=The Coulomb interaction in reciprocal space. A cutoff can be used
!!  Psps<Pseudopotential_type>=Variables related to pseudopotentials
!!  Pawtab(Psps%ntypat)<pawtab_type>=PAW tabulated starting data.
!!  Paw_pwff(Cryst%ntypat*Wfd%usepaw)<pawpwff_t>=Form factor used to calculate the onsite mat. elements of a plane wave.
!!  lomo_spin(Wfd%nsppol)=Lowest occupied band for each spin
!!  homo_spin(Wfd%nsppol)=Highest occupied band for each spin
!!  humo_spin(Wfd%nsppol)=Highest unoccupied band for each spin
!!  nfftot_osc=Total Number of FFT points used for the oscillator matrix elements.
!!  ngfft_osc(18)=Info on the FFT algorithm used to calculate the oscillator matrix elements.
!!  npweps=Number of G-vectors in mgq0.
!!
!! OUTPUT
!!   mgq0(npweps,lomo_min:humo_max,lomo_min:humo_max,Wfd%nkibz,Wfd%nsppol)
!!     Allocated here and filled with the matrix elements on each node.
!!
!! PARENTS
!!      exc_build_ham
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_all_mgq0(Wfd,Cryst,Qmesh,Gsph_x,Vcp,&
& Psps,Pawtab,Paw_pwff,lomo_spin,homo_spin,humo_spin,nfftot_osc,ngfft_osc,npweps,mgq0)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot_osc,npweps
 type(kmesh_t),intent(in) :: Qmesh
 type(crystal_t),intent(in) :: Cryst
 type(vcoul_t),intent(in) :: Vcp
 type(gsphere_t),intent(in) :: Gsph_x
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfd_t),target,intent(inout) :: Wfd
!arrays
 integer,intent(in) :: lomo_spin(Wfd%nsppol),homo_spin(Wfd%nsppol),humo_spin(Wfd%nsppol)
 integer,intent(in) :: ngfft_osc(18)
 complex(gwpc),allocatable,intent(out) :: mgq0(:,:,:,:,:)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: map2sphere1=1,dim_rtwg1=1,ndat1=1
 integer :: use_padfft,mgfft_osc,fftalga_osc,ii
 integer :: ik_ibz,itim_k,isym_k,iq_bz,iq_ibz,isym_q,itim_q,iqbz0
 integer :: ierr,iv,ic,spin,lomo_min,humo_max !,inv_ipw,ipw
 real(dp) :: cpu,wall,gflops !q0vol,fcc_const
 complex(dpc) :: ph_mkt
 character(len=500) :: msg
!arrays
 integer,allocatable :: igfftg0(:),task_distrib(:,:,:,:)
 integer,allocatable :: gbound(:,:),id_tab(:)
 real(dp) :: qbz(3),spinrot_k(4),tsec(2)
 complex(gwpc),allocatable :: rhotwg1(:)
 complex(gwpc),target,allocatable :: ur1(:),ur2(:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: ptr_ur1(:),ptr_ur2(:)
 type(pawcprj_type),allocatable :: Cp1(:,:),Cp2(:,:)
 type(pawpwij_t),allocatable :: Pwij_q0(:)

!************************************************************************

 call timab(671,1,tsec)

 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")
 ABI_CHECK(nfftot_osc==PRODUCT(ngfft_osc(1:3)),"mismatch in FFT size")

 lomo_min = MINVAL(lomo_spin); humo_max = MAXVAL(humo_spin)

 if ( ANY(ngfft_osc(1:3) /= Wfd%ngfft(1:3)) ) then
   call wfd%change_ngfft(Cryst,Psps,ngfft_osc)
 end if

 mgfft_osc   = MAXVAL(ngfft_osc(1:3))
 fftalga_osc = ngfft_osc(7)/100 !; fftalgc_osc=MOD(ngfft_osc(7),10)

 ! (temporary) Table used for the wavefunction in the IBZ.
 ABI_MALLOC(id_tab, (Wfd%nfft))
 id_tab = (/(ii, ii=1,Wfd%nfft)/)

 ! Analytic integration of 4pi/q^2 over the volume element:
 ! $4pi/V \int_V d^3q 1/q^2 =4pi bz_geometric_factor V^(-2/3)$
 ! i_sz=4*pi*bz_geometry_factor*q0_vol**(-two_thirds) where q0_vol= V_BZ/N_k
 ! bz_geometry_factor: sphere=7.79, fcc=7.44, sc=6.188, bcc=6.946, wz=5.255
 ! (see gwa.pdf, appendix A.4)

 ! If q=0 and C=V then set up rho-twiddle(G=0) to reflect an
 ! analytic integration of q**-2 over the volume element:
 ! <q**-2> = 7.44 V**(-2/3)   (for fcc cell)

 ! q0vol = (8.0*pi**3) / (Cryst%ucvol*Kmesh%nbz)
 ! fcc_const = SQRT(7.44*q0vol**(-2.0/3.0))
 ! rtw = (6.0*pi**2/(Cryst%ucvol*Kmesh%nkbz))**(1./3.)
 ! Average of (q+q')**-2 integration for head of Coulomb matrix
 ! INTRTW(QL) = (2*pi*rtw + pi*(rtw**2/QL-QL)*LOG((QL+rtw)/(QL-rtw)))
 ! &              * (Cryst%ucvol*Kmesh%nbz)/(2*pi)**3. * QL*QL

 if (Wfd%usepaw==1) then
   ABI_DT_MALLOC(Cp1,(Wfd%natom,Wfd%nspinor))
   call pawcprj_alloc(Cp1,0,Wfd%nlmn_atm)
   ABI_DT_MALLOC(Cp2,(Wfd%natom,Wfd%nspinor))
   call pawcprj_alloc(Cp2,0,Wfd%nlmn_atm)
 end if

 ABI_MALLOC(ur1,(nfftot_osc*Wfd%nspinor))
 ABI_MALLOC(ur2,(nfftot_osc*Wfd%nspinor))

 ! Identify q==0
 iqbz0=0
 do iq_bz=1,Qmesh%nbz
   if (ALL(ABS(Qmesh%bz(:,iq_bz))<tol3)) iqbz0 = iq_bz
 end do
 ABI_CHECK(iqbz0/=0,"q=0 not found in q-point list!")

 ! * Get iq_ibz, and symmetries from iqbz0.
 call get_BZ_item(Qmesh,iqbz0,qbz,iq_ibz,isym_q,itim_q)

 if (Wfd%usepaw==1) then ! Prepare onsite contributions at q==0
   ABI_DT_MALLOC(Pwij_q0,(Cryst%ntypat))
   call pawpwij_init(Pwij_q0,npweps,Qmesh%bz(:,iqbz0),Gsph_x%gvec,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
 end if
 !
 ! Tables for the FFT of the oscillators.
 !  a) FFT index of the G sphere (only vertical transitions, unlike cchi0, no need to shift the sphere).
 !  b) gbound table for the zero-padded FFT performed in rhotwg.
 ABI_MALLOC(igfftg0,(Gsph_x%ng))
 ABI_MALLOC(gbound,(2*mgfft_osc+8,2))
 call gsph_fft_tabs(Gsph_x,(/0,0,0/),mgfft_osc,ngfft_osc,use_padfft,gbound,igfftg0)
 if ( ANY(fftalga_osc == (/2,4/)) ) use_padfft=0 ! Pad-FFT is not coded in rho_tw_g
#ifdef FC_IBM
 ! XLF does not deserve this optimization (problem with [v67mbpt][t03])
 use_padfft = 0
#endif
 if (use_padfft==0) then
   ABI_FREE(gbound)
   ABI_MALLOC(gbound,(2*mgfft_osc+8,2*use_padfft))
 end if

 ABI_MALLOC(rhotwg1,(npweps))

 ABI_MALLOC_OR_DIE(mgq0, (npweps,lomo_min:humo_max,lomo_min:humo_max,Wfd%nkibz,Wfd%nsppol), ierr)
 mgq0 = czero

 call cwtime(cpu,wall,gflops,"start")

 do spin=1,Wfd%nsppol
   ! Distribute the calculation of the matrix elements.
   ! processors have the entire set of wavefunctions hence we divide the workload
   ! without checking if the pair of states is available. Last dimension is fake.
   ABI_MALLOC(task_distrib,(lomo_spin(spin):humo_spin(spin),lomo_spin(spin):humo_spin(spin),Wfd%nkibz,1))
   call xmpi_distab(Wfd%nproc,task_distrib)

   ! loop over the k-points in IBZ
   do ik_ibz=1,Wfd%nkibz
     if ( ALL(task_distrib(:,:,ik_ibz,1)/= Wfd%my_rank) ) CYCLE

     ! Don't need to symmetrize the wavefunctions.
     itim_k=1; isym_k=1; ph_mkt=cone; spinrot_k=Cryst%spinrot(:,isym_k)

     do iv=lomo_spin(spin),humo_spin(spin) ! Loop over band V
       if ( ALL(task_distrib(:,iv,ik_ibz,1)/=Wfd%my_rank) ) CYCLE

       if (wfd%ihave_ur(iv,ik_ibz,spin,how="Stored")) then
         ptr_ur1 =>  Wfd%Wave(iv,ik_ibz,spin)%ur
       else
         call wfd%get_ur(iv,ik_ibz,spin,ur1)
         ptr_ur1 =>  ur1
       end if

       if (Wfd%usepaw==1) then
         call wfd%get_cprj(iv,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)
       end if

       ! Loop over band C
       do ic=lomo_spin(spin),humo_spin(spin)
         if ( task_distrib(ic,iv,ik_ibz,1)/=Wfd%my_rank ) CYCLE

         if (wfd%ihave_ur(ic,ik_ibz,spin,how="Stored")) then
           ptr_ur2 =>  Wfd%Wave(ic,ik_ibz,spin)%ur
         else
           call wfd%get_ur(ic,ik_ibz,spin,ur2)
           ptr_ur2 =>  ur2
         end if

         if (Wfd%usepaw==1) then
           call wfd%get_cprj(ic,ik_ibz,spin,Cryst,Cp2,sorted=.FALSE.)
         end if

         call rho_tw_g(Wfd%nspinor,npweps,nfftot_osc,ndat1,ngfft_osc,map2sphere1,use_padfft,igfftg0,gbound,&
&          ptr_ur1,1,id_tab,ph_mkt,spinrot_k,&
&          ptr_ur2,1,id_tab,ph_mkt,spinrot_k,&
&          dim_rtwg1,rhotwg1)

         if (Wfd%usepaw==1) then ! Add PAW onsite contribution.
           call paw_rho_tw_g(npweps,dim_rtwg1,Wfd%nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_x%gvec,&
&            Cp1,Cp2,Pwij_q0,rhotwg1)
         end if

         ! If q=0 treat Exchange and Coulomb-term independently
         if (iv <= homo_spin(spin) .and. ic <= homo_spin(spin) .or. &
&            iv >  homo_spin(spin) .and. ic >  homo_spin(spin)) then

           if (iv/=ic) then !COULOMB term: C/=V: ignore them
             rhotwg1(1) = czero_gw
           else
             ! If q=0 and C=V then set up rho-twiddle(G=0) to reflect an
             ! analytic integration of q**-2 over the volume element:
             ! <q**-2> = 7.44 V**(-2/3)   (for fcc cell)
             !rhotwg1(1) = fcc_const * qpg(1,iqbz0)
             rhotwg1(1) = SQRT(GWPC_CMPLX(Vcp%i_sz,zero)) / Vcp%vcqlwl_sqrt(1,1)
             !if (vcut) rhotwg1(1) = 1.0
           end if

         else
           ! At present this term is set to zero
           ! EXCHANGE term: limit value.
           ! Set up rho-twiddle(G=0) using small vector q instead of zero and k.p perturbation theory (see notes)
           rhotwg1(1) = czero_gw
         end if

         mgq0(:,iv,ic,ik_ibz,spin) = rhotwg1(:)
       end do !ic
     end do !iv
   end do !ik_ibz

   ABI_FREE(task_distrib)
 end do !spin

 ! TODO: One can speedup the calculation by computing the upper triangle of the
 ! matrix in (b,b') space and then take advantage of the symmetry property:
 !
 !   M_{k,0}{{bb'}(G)^* = M{k,0}{b'b'}(-G)

#if 0
 !!!! $OMP PARALLEL DO COLLAPSE(3) PRIVATE(inv_ipw)
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     do iv=lomo_spin(spin),humo_spin(spin)
       do ic=1,iv-1
         do ipw=1,npweps
           inv_ipw = gsph_x%g2mg(ipw)
           mgq0(inv_ipw,ic,iv,ik_ibz,spin) = mgq0(ipw,iv,ic,ik_ibz,spin)
         end do
       end do
     end do
   end do
 end do
#endif
 !
 ! Gather matrix elements on each node.
 call xmpi_sum(mgq0,Wfd%comm,ierr)

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f9.6))')"cpu_time = ",cpu,", wall_time = ",wall
 call wrtout(std_out,msg,"PERS")

 ABI_FREE(rhotwg1)
 ABI_FREE(igfftg0)
 ABI_FREE(gbound)
 ABI_FREE(ur1)
 ABI_FREE(ur2)
 ABI_FREE(id_tab)

 if (Wfd%usepaw==1) then
   ! Deallocation for PAW.
   call pawpwij_free(Pwij_q0)
   ABI_DT_FREE(Pwij_q0)
   call pawcprj_free(Cp1)
   ABI_DT_FREE(Cp1)
   call pawcprj_free(Cp2)
   ABI_DT_FREE(Cp2)
 end if

 call timab(671,2,tsec)

end subroutine wfd_all_mgq0
!!***

end module m_exc_build
!!***
