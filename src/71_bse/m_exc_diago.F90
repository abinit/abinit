!!****m* ABINIT/m_exc_diago
!! NAME
!! m_exc_diago
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT and EXC groups (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
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

MODULE m_exc_diago

 use defs_basis
 use m_slk
 use m_bs_defs
 use m_abicore
 use m_errors
 use m_xmpi
#if defined HAVE_MPI2
 use mpi
#endif
 use m_hdr
 use m_sort

 use defs_datatypes,    only : pseudopotential_type, ebands_t
 use m_io_tools,        only : open_file
 use m_fstrings,        only : int2char4
 use m_numeric_tools,   only : print_arr, hermitianize
 use m_crystal,         only : crystal_t
 use m_kpts,            only : listkk
 use m_bz_mesh,         only : kmesh_t
 use m_ebands,          only : ebands_report_gap
 use m_eprenorms,       only : eprenorms_t
 use m_wfd,             only : wfd_t
 use m_paw_hr,          only : pawhur_t
 use m_pawtab,          only : pawtab_type
 use m_exc_itdiago,     only : exc_iterative_diago
 use m_hide_lapack,     only : xheev, xheevx, xgeev, xhegvx, xginv, xhdp_invert, xhegv
 use m_hide_blas,       only : xdotc, xgemm
 use m_bse_io,          only : exc_fullh_from_blocks, offset_in_file, rrs_of_glob, ccs_of_glob, &
&                              exc_read_bshdr, exc_skip_bshdr_mpio, exc_read_rblock_fio
 use m_exc_spectra,     only : build_spectra

 implicit none

 private

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!#define DEV_MG_DEBUG_THIS

 public ::  exc_diago_driver ! Driver routine for the direct diagonalization of the BSE Hamiltonian (main entry point)
!!***

contains

!!****f* m_exc_diago/exc_diago_driver
!! NAME
!!  exc_diago_driver
!!
!! FUNCTION
!!  Driver routine for the direct diagonalization of the Hermitian excitonic Hamiltonian.
!!
!! INPUTS
!!  neh=Rank of the resonant block of the Hamiltonian.
!!  BS_files<excfiles>=Datatype storing names and files used in the Bethe-Salpeter code.
!!    %exh=Name of the file storing the excitonic resonant part.
!!
!! OUTPUT
!!  Eigenvalues and eigenvectors are written on file.
!!
!! PARENTS
!!      m_bethe_salpeter
!!
!! CHILDREN
!!
!! SOURCE

subroutine exc_diago_driver(Wfd,Bsp,BS_files,KS_BSt,QP_BSt,Cryst,Kmesh,Psps,&
&  Pawtab,Hur,Hdr_bse,drude_plsmf,Epren)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: drude_plsmf
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) ::  BS_files
 type(Hdr_type),intent(in) :: Hdr_bse
 type(crystal_t),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 type(kmesh_t),intent(in) :: Kmesh
 type(ebands_t),intent(in) :: KS_BSt,QP_BSt
 type(wfd_t),intent(inout) :: Wfd
 type(eprenorms_t),intent(in) :: Epren
!arrays
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(pawhur_t),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer :: my_rank,master,comm,prtvol
 complex(dpc) :: exc_gap,gw_gap
 logical :: eval_eigenstates
 character(len=500) :: msg
!arrays
 real(dp) :: gaps(3,QP_BSt%nsppol)

!************************************************************************

 DBG_ENTER("COLL")

 comm    = Wfd%comm
 my_rank = Wfd%my_rank
 master  = Wfd%master
 prtvol  = Wfd%prtvol

 if (BSp%have_complex_ene) then
   ABI_ERROR("Complex energies are not supported yet")
 end if
 !
 ! This trick is needed to restart a CG run, use DDIAGO to calculate the spectra reusing an old BSEIG file.
 eval_eigenstates = (BS_files%in_eig == BSE_NOFILE) .or. (Bsp%algorithm == BSE_ALGO_CG)

 if (eval_eigenstates) then
   !
   select case (BSp%algorithm)
   case (BSE_ALGO_DDIAGO)
     if (BSp%use_coupling==0) then
       call exc_diago_resonant(BSp,BS_files,Hdr_bse,prtvol,comm)
       if(Bsp%do_ep_renorm) then
         call exc_diago_resonant(BSp,BS_files,Hdr_bse,prtvol,comm,Epren=Epren,Kmesh=Kmesh,Cryst=Cryst,elph_lifetime=.TRUE.)
       end if
     else
       if (Bsp%have_complex_ene) then
         ! Solve Hv = ev with generic complex matrix.
         call exc_diago_coupling(BSp,BS_files,Hdr_bse,prtvol,comm)
       else
         ! Solve generalized eigenvalue problem F Hbar with Hbar Hermitian definitive positive matrix.
         call exc_diago_coupling_hegv(BSp,BS_files,Hdr_bse,prtvol,comm)
       end if
     end if

   case (BSE_ALGO_CG)
     if (BSp%use_coupling==0) then
       call exc_iterative_diago(Bsp,BS_files,Hdr_bse,prtvol,comm)
     else
       ABI_ERROR("CG + coupling not coded")
     end if

   case default
     write(msg,'(a,i0)')" Wrong value for Bsp%algorithm: ",Bsp%algorithm
     ABI_ERROR(msg)
   end select
   !
   if (my_rank==master) then
     call ebands_report_gap(QP_BSt,header="QP bands",unit=std_out,gaps=gaps)
     gw_gap = MINVAL(gaps(2,:))
     call exc_print_eig(BSp,BS_files%out_eig,gw_gap,exc_gap)
   end if
   call xmpi_barrier(comm)
   !
 end if

 call build_spectra(BSp,BS_files,Cryst,Kmesh,KS_BSt,QP_BSt,Psps,Pawtab,Wfd,Hur,drude_plsmf,comm)

 ! Electron-phonon renormalization !
 if (BSp%algorithm == BSE_ALGO_DDIAGO .and. BSp%use_coupling == 0 .and. BSp%do_ep_renorm) then
   call build_spectra(BSp,BS_files,Cryst,Kmesh,KS_BSt,QP_BSt,Psps,Pawtab,Wfd,Hur,drude_plsmf,comm,Epren=Epren)
 end if

 DBG_EXIT("COLL")

end subroutine exc_diago_driver
!!***

!----------------------------------------------------------------------

!!****f* m_exc_diago/exc_diago_resonant
!! NAME
!!  exc_diago_resonant
!!
!! FUNCTION
!!  Calculates eigenvalues and eigenvectors of the Hermitian excitonic Hamiltonian (coupling is neglected).
!!
!! INPUTS
!!  Bsp
!!  BS_files<excfiles>=Datatype storing names and files used in the Bethe-Salpeter code.
!!  comm=MPI communicator.
!!  bseig_fname=The name of the output file
!!  prtvol=Verbosity level.
!!
!! OUTPUT
!!  Eigenvalues and eigenvectors are written on file bseig_fname
!!
!! PARENTS
!!      m_exc_diago
!!
!! CHILDREN
!!
!! SOURCE

subroutine exc_diago_resonant(Bsp,BS_files,Hdr_bse,prtvol,comm,Epren,Kmesh,Cryst,elph_lifetime)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,prtvol
 logical,optional,intent(in) :: elph_lifetime
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(Hdr_type),intent(in) :: Hdr_bse
 type(eprenorms_t),optional,intent(in) :: Epren
 type(kmesh_t),optional,intent(in) :: Kmesh
 type(crystal_t),optional,intent(in) :: Cryst

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: ii,it,mi,hreso_unt,eig_unt,exc_size,neh1,neh2,j
 integer :: nsppol,il,iu,mene_found,nstates
 integer :: nprocs,my_rank,fform,nene_printed,ierr
 real(dp) :: exc_gap,exc_maxene,abstol
 real(dp) :: vl,vu
 logical :: use_scalapack,do_full_diago,diagonal_is_real
 character(len=500) :: msg
 character(len=fnlen) :: hreso_fname,bseig_fname
!arrays
 real(dp),allocatable :: exc_ene(:)
 complex(dpc),allocatable :: exc_mat(:,:),exc_vec(:,:)
#if defined HAVE_LINALG_SCALAPACK && defined HAVE_MPI_IO
 integer :: amode,mpi_fh,istwf_k,tbloc,tmp_unt
 integer :: itloc,jj,jtloc,itglob,jtglob
 integer(XMPI_OFFSET_KIND) :: ehdr_offset,fmarker
 integer :: block_sizes(2,3),array_of_sizes(2),gsub(2,2)
 logical,parameter :: is_fortran_file=.TRUE.
 real(dp),external :: PDLAMCH
 type(matrix_scalapack)    :: Slk_mat,Slk_vec
 type(processor_scalapack) :: Slk_processor
#endif

 integer :: ik, ic, iv, isppol, ireh, ep_ik, itemp
 complex(dpc) :: en

 real(dp) :: dksqmax
 integer,allocatable :: bs2eph(:,:)
 integer :: sppoldbl, timrev
 logical :: do_ep_renorm, do_ep_lifetime
 integer :: ntemp
 character(len=4) :: ts
 complex(dpc),allocatable :: exc_vl(:,:),exc_ene_c(:)
 complex(dpc) :: ctemp
!! complex(dpc),allocatable :: ovlp(:,:)

!************************************************************************

 DBG_ENTER("PERS")

 nprocs  = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 if (BSp%have_complex_ene) then ! QP lifetimes are not included
   ABI_ERROR("complex energies not coded yet")
 end if

 if (ANY(Bsp%nreh/=Bsp%nreh(1))) then
   write(std_out,*)" Bsp%nreh: ",Bsp%nreh
   write(msg,'(a)')" BSE code does not support different number of transitions for the two spin channels"
   ABI_WARNING(msg)
 end if

 nsppol   = Hdr_bse%nsppol
 exc_size = SUM(BSp%nreh)
 nstates  = BSp%nstates; do_full_diago=(Bsp%nstates==exc_size)

 neh1 = Bsp%nreh(1); neh2 = neh1
 if (Hdr_bse%nsppol==2) neh2 = Bsp%nreh(2)

 ! Scalapack is disabled due to portability issues in slk_read
 ! This part should be rewritten  with hdf5 + mpi-io

 use_scalapack = .FALSE.
!#if defined HAVE_LINALG_SCALAPACK
! use_scalapack = (nprocs > 1)
!#endif
 if (use_scalapack .and. nsppol == 2) then
   use_scalapack = .False.
   msg = "Scalapack with nsppol==2 not yet available. Using sequential version"
   ABI_WARNING(msg)
 end if

 if (.not.use_scalapack .and. my_rank/=master) GOTO 10 ! Inversion is done by master only.

 nene_printed = MIN(32*nsppol,nstates); if (prtvol>10) nene_printed = nstates

 if (BS_files%in_hreso /= BSE_NOFILE) then
   hreso_fname = BS_files%in_hreso
 else
   hreso_fname = BS_files%out_hreso
 end if

 bseig_fname = BS_files%out_eig
 if (BS_files%in_eig /= BSE_NOFILE) then
   ABI_ERROR("BS_files%in_eig is defined!")
 end if

 write(msg,'(a,i0)')' Direct diagonalization of the resonant excitonic Hamiltonian, Matrix size= ',exc_size
 call wrtout([std_out, ab_out], msg)

 ABI_MALLOC_OR_DIE(exc_ene,(exc_size), ierr)

 ABI_MALLOC_OR_DIE(exc_ene_c,(exc_size), ierr)

 do_ep_renorm = .FALSE.
 ntemp = 1
 do_ep_lifetime = .FALSE.
 if(BSp%do_ep_renorm .and. present(Epren)) then
   do_ep_renorm = .TRUE.
   ntemp = Epren%ntemp
   if(present(elph_lifetime)) then
     do_ep_lifetime = elph_lifetime
   end if
 end if

 if (do_ep_renorm) then
   ABI_CHECK(nsppol == 1, "Nsppol == 2 not supported with elphon renormalizations")
 end if

 SELECT CASE (use_scalapack)
 CASE (.FALSE.)

   write(msg,'(a)')". Using LAPACK sequential version. "
   call wrtout([std_out, ab_out], msg)
   write(msg,'(a,f8.1,a)')' Allocating excitonic eigenvalues. Memory required: ', exc_size*dp*b2Mb,' Mb. '
   call wrtout(std_out, msg)

   write(msg,'(a,f8.1,a)')' Allocating excitonic hamiltonian.  Memory required: ',exc_size**2*dpc*b2Mb,' Mb.'
   call wrtout(std_out, msg, do_flush=.True.)

   ABI_MALLOC_OR_DIE(exc_mat,(exc_size,exc_size), ierr)

   if (do_ep_renorm) then
     ABI_MALLOC_OR_DIE(exc_vl,(exc_size,exc_size),ierr)
   end if
   !exc_mat = HUGE(zero)
   !
   ! Read data from file.
   if (open_file(hreso_fname,msg,newunit=hreso_unt,form="unformatted",status="old",action="read") /= 0) then
     ABI_ERROR(msg)
   end if
   !
   ! Read the header and perform consistency checks.
   call exc_read_bshdr(hreso_unt,Bsp,fform,ierr)
   ABI_CHECK(ierr==0,"Fatal error, cannot continue")
   !
   ! Construct full resonant block using Hermiticity.
   diagonal_is_real = .not.Bsp%have_complex_ene
   call exc_read_rblock_fio(hreso_unt,diagonal_is_real,nsppol,Bsp%nreh,exc_size,exc_mat,ierr)
   ABI_CHECK(ierr==0,"Fatal error, cannot continue")

   close(hreso_unt)

   if (do_ep_renorm) then
     write(std_out,'(a)') "Mapping kpts from bse to eph"
     sppoldbl = 1 !; if (any(Cryst%symafm == -1) .and. Epren%nsppol == 1) nsppoldbl=2
     ABI_MALLOC(bs2eph, (Kmesh%nbz*sppoldbl, 6))
     timrev = 1
     call listkk(dksqmax, Cryst%gmet, bs2eph, Epren%kpts, Kmesh%bz, Epren%nkpt, Kmesh%nbz, Cryst%nsym, &
        sppoldbl, Cryst%symafm, Cryst%symrel, timrev, xmpi_comm_self, use_symrec=.False.)
   end if

   do itemp = 1, ntemp

     !TODO should find a way not to read again and again !
     !  but without storing it twice !!!
     !exc_mat = HUGE(zero)
     !
     ! Read data from file.
     if (open_file(hreso_fname,msg,newunit=hreso_unt,form="unformatted",status="old",action="read") /= 0) then
       ABI_ERROR(msg)
     end if
     !
     ! Read the header and perform consistency checks.
     call exc_read_bshdr(hreso_unt,Bsp,fform,ierr)
     ABI_CHECK(ierr==0,"Fatal error, cannot continue")
     !
     ! Construct full resonant block using Hermiticity.
     diagonal_is_real = .not.Bsp%have_complex_ene
     call exc_read_rblock_fio(hreso_unt,diagonal_is_real,nsppol,Bsp%nreh,exc_size,exc_mat,ierr)
     ABI_CHECK(ierr==0,"Fatal error, cannot continue")

     close(hreso_unt)

     bseig_fname = BS_files%out_eig

     if (do_ep_renorm) then
       write(std_out,'(a,i4)') "Will perform elphon renormalization for itemp = ",itemp

       call int2char4(itemp,ts)

       bseig_fname = TRIM(BS_files%out_eig) // TRIM("_T") // ts

       ! Should patch the diagonal of exc_mat

       do isppol = 1, BSp%nsppol
         do ireh = 1, BSp%nreh(isppol)
           ic = BSp%Trans(ireh,isppol)%c
           iv = BSp%Trans(ireh,isppol)%v
           ik = BSp%Trans(ireh,isppol)%k ! In the full bz
           en = BSp%Trans(ireh,isppol)%en

           ep_ik = bs2eph(ik,1)

           !TODO support multiple spins !
           if(ABS(en - (Epren%eigens(ic,ep_ik,isppol)-Epren%eigens(iv,ep_ik,isppol)+BSp%mbpt_sciss)) > tol3) then
             ABI_ERROR("Eigen from the transition does not correspond to the EP file !")
           end if
           exc_mat(ireh,ireh) = exc_mat(ireh,ireh) + (Epren%renorms(1,ic,ik,isppol,itemp) - Epren%renorms(1,iv,ik,isppol,itemp))

           ! Add lifetime
           if(do_ep_lifetime) then
             exc_mat(ireh,ireh) = exc_mat(ireh,ireh) - j_dpc*(Epren%linewidth(1,ic,ik,isppol,itemp) &
&                + Epren%linewidth(1,iv,ik,isppol,itemp))
           end if

         end do
       end do

     end if

     if (do_full_diago) then
       if(do_ep_renorm) then
         call wrtout(std_out," Full diagonalization with XGEEV... ")
         ABI_MALLOC(exc_vec,(exc_size,exc_size))
         call xgeev('V','V',exc_size,exc_mat,exc_size,exc_ene_c,exc_vl,exc_size,exc_vec,exc_size)
         exc_mat(:,1:nstates) = exc_vec
         ABI_FREE(exc_vec)
       else
         call wrtout(std_out," Full diagonalization with XHEEV... ")
         call xheev("Vectors","Upper",exc_size,exc_mat,exc_ene)
         exc_ene_c(:) = exc_ene(:)
       end if
     else
       call wrtout(std_out," Partial diagonalization with XHEEVX... ")
       abstol=zero; il=1; iu=nstates
       ABI_MALLOC_OR_DIE(exc_vec,(exc_size,nstates),ierr)
       call xheevx("Vectors","Index","Upper",exc_size,exc_mat,vl,vu,il,iu,abstol,mene_found,exc_ene,exc_vec,exc_size)
       exc_mat(:,1:nstates) = exc_vec
       exc_ene_c(:) = exc_ene(:)
       ABI_FREE(exc_vec)
     end if
     !
     ! ==============================================
     ! === Now exc_mat contains the eigenvectors ====
     ! ==============================================

     ! * Write the final results.
     call wrtout(std_out,' Writing eigenvalues and eigenvectors to file: '//TRIM(bseig_fname))

     if (open_file(bseig_fname,msg,newunit=eig_unt,form="unformatted",action="write") /= 0) then
       ABI_ERROR(msg)
     end if

     !!! !DBYG
     !!! !Compute overlap matrix
     !!! ABI_MALLOC(ovlp,(exc_size,exc_size))
     !!! do mi=1,nstates
     !!!   do ireh=1,nstates
     !!!     ovlp(mi,ireh) = xdotc(exc_size,exc_vl(:,mi),1,exc_mat(:,ireh),1)
     !!!     if(mi==ireh) then
     !!!       !if(ABS(ovlp(mi,ireh)) < 0.999) then
     !!!       !  write(*,*) "it,itp = ",mi,ireh,"ovlp = ",ovlp(mi,ireh)
     !!!       !end if
     !!!     else
     !!!       if(ABS(ovlp(mi,ireh)) > 0.001) then
     !!!         write(*,*) "it,itp = ",mi,ireh,"ovlp = ",ovlp(mi,ireh)
     !!!       end if
     !!!     end if
     !!!   end do
     !!! end do
     !!! !call xgemm("C","N",exc_size,nstates,nstates,cone,exc_vl,exc_size,exc_mat,exc_size,czero,ovlp,nstates)

     !!! write(777,*) ovlp
     !!! ABI_FREE(ovlp)
     !!! !ENDDBYG

     !% fform = 1002 ! FIXME
     !% call hdr_io_int(fform,Hdr_bse,2,eig_unt)

     write(eig_unt) do_ep_lifetime
     write(eig_unt) exc_size, nstates
     write(eig_unt) exc_ene_c(1:nstates)
     do mi=1,nstates
       write(eig_unt) exc_mat(1:exc_size,mi)
       if(do_ep_lifetime) then
         write(eig_unt) exc_vl(1:exc_size,mi)
       end if
     end do

     close(eig_unt)

   end do ! itemp

   ABI_FREE(exc_mat)
   if (do_ep_renorm) then
     ABI_FREE(exc_vl)
     ABI_FREE(bs2eph)
   end if

 CASE (.TRUE.)

#if defined HAVE_LINALG_SCALAPACK && defined HAVE_MPI_IO
   if (nsppol==2) then
     ABI_WARNING("nsppol==2 + scalapack not coded yet")
   end if

   istwf_k=1; tbloc=50
   write(msg,'(2(a,i0))')". Using scaLAPACK version with nprocs= ",nprocs,"; block size= ",tbloc
   call wrtout([std_out, ab_out], msg)

   write(msg,'(a,f8.1,a)')' Allocating excitonic eigenvalues. Memory required: ',exc_size*dp*b2Mb,' Mb. '
   call wrtout(std_out, msg)
   !
   ! Init scaLAPACK environment.
   call init_scalapack(Slk_processor,comm)
   !
   ! Init scaLAPACK matrices
   call init_matrix_scalapack(Slk_mat,exc_size,exc_size,Slk_processor,istwf_k,tbloc=tbloc)

   call init_matrix_scalapack(Slk_vec,exc_size,exc_size,Slk_processor,istwf_k,tbloc=tbloc)
   !
   ! Open the file with MPI-IO and skip the record.
   amode=MPI_MODE_RDONLY

   call MPI_FILE_OPEN(comm, hreso_fname, amode, MPI_INFO_NULL, mpi_fh, ierr)
   ABI_CHECK_MPI(ierr,"MPI_IO error opening file: "//TRIM(hreso_fname))

   ! Skip the header and find the offset for reading the matrix.
   call exc_skip_bshdr_mpio(mpi_fh,xmpio_collective,ehdr_offset)
   !
   ! Read scaLAPACK matrix from the file.
   if (nsppol==1) then
     call slk_read(Slk_mat,"Upper","Hermitian",is_fortran_file,mpi_fh=mpi_fh,offset=ehdr_offset)
   else
     array_of_sizes = (/exc_size,exc_size/)
     block_sizes(:,1) = (/neh1,neh1/)
     block_sizes(:,2) = (/neh2,neh2/)
     block_sizes(:,3) = (/neh1,neh2/)
     ABI_ERROR("Not tested")
     !call slk_read_from_blocks(Slk_mat,array_of_sizes,block_sizes,is_fortran_file,mpi_fh=mpi_fh,offset=ehdr_offset)
   end if

   call MPI_FILE_CLOSE(mpi_fh, ierr)
   ABI_CHECK_MPI(ierr,"FILE_CLOSE")

   if (do_full_diago) then
     call wrtout(std_out," Performing full diagonalization with scaLAPACK...")

     call slk_pzheev("Vectors","Upper",Slk_mat,Slk_vec,exc_ene)
   else
     call wrtout(std_out," Performing partial diagonalization with scaLAPACK...")
     il=1; iu=nstates; abstol=zero !ABSTOL = PDLAMCH(comm,'U')
     call slk_pzheevx("Vectors","Index","Upper",Slk_mat,vl,vu,il,iu,abstol,Slk_vec,mene_found,exc_ene)
   end if

   exc_ene_c(:) = exc_ene(:)

   call Slk_mat%free()

   call wrtout(std_out,' Writing eigenvalues/vectors to file: '//TRIM(bseig_fname), do_flush=.True.)

   ! Write distributed matrix on file bseig_fname with Fortran records.
   if (my_rank==master) then ! Write exc eigenvalues. Vectors will be appended in slk_write.
     if (open_file(bseig_fname,msg,newunit=eig_unt,form="unformatted",action="write") /= 0) then
       ABI_ERROR(msg)
     end if
     write(eig_unt) exc_size, nstates
     write(eig_unt) exc_ene_c(1:nstates)
     close(eig_unt)
   end if

   call xmpi_barrier(comm)
   !
   ! Open the file with MPI-IO and skip the record.
   amode=MPI_MODE_RDWR

   call MPI_FILE_OPEN(comm, bseig_fname, amode, MPI_INFO_NULL, mpi_fh, ierr)
   ABI_CHECK_MPI(ierr,"MPI_IO error opening file: "//TRIM(hreso_fname))

   !call MPI_FILE_SYNC(mpi_fh,ierr)

   ehdr_offset = 0
   call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_collective,fmarker,ierr)
   call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_collective,fmarker,ierr)

   write(std_out,*)"Writing nstates ",nstates
   gsub(:,1) = (/1,1/)
   gsub(:,2) = (/exc_size,nstates/)
   call slk_write(Slk_vec,"All",is_fortran_file,mpi_fh=mpi_fh,offset=ehdr_offset,glob_subarray=gsub)

   call MPI_FILE_CLOSE(mpi_fh, ierr)
   ABI_CHECK_MPI(ierr,"FILE_CLOSE")

   call Slk_vec%free()
   call end_scalapack(Slk_processor)
   call xmpi_barrier(comm)
#else
   ABI_BUG("You should not be here!")
#endif

 END SELECT

 ! Order the eigenvalues
 do ii=nstates,2,-1
   do j=1,ii-1
     if (DBLE(exc_ene_c(j)) > DBLE(exc_ene_c(j+1))) then
       ctemp = exc_ene_c(j)
       exc_ene_c(j) = exc_ene_c(j+1)
       exc_ene_c(j+1) = ctemp
     end if
   end do
 end do


 write(msg,'(a,i4)')' Excitonic eigenvalues in eV up to n= ',nene_printed
 call wrtout([std_out, ab_out], msg)

 do it=0,(nene_printed-1)/8
   write(msg,'(8f10.5)') ( DBLE(exc_ene_c(ii))*Ha_eV, ii=1+it*8,MIN(it*8+8,nene_printed) )
   call wrtout([std_out, ab_out], msg)
 end do

 exc_gap    = MINVAL(DBLE(exc_ene_c(1:nstates)))
 exc_maxene = MAXVAL(DBLE(exc_ene_c(1:nstates)))

 write(msg,'(a,2(a,f7.2,2a),a)')ch10,&
  " First excitonic eigenvalue= ",exc_gap*Ha_eV,   " [eV]",ch10,&
  " Last  excitonic eigenvalue= ",exc_maxene*Ha_eV," [eV]",ch10,ch10
 call wrtout([std_out, ab_out], msg, do_flush=.True.)

 ABI_FREE(exc_ene_c)
 ABI_FREE(exc_ene)

10 call xmpi_barrier(comm)

 DBG_EXIT("PERS")

end subroutine exc_diago_resonant
!!***

!----------------------------------------------------------------------

!!****f* m_exc_diago/exc_print_eig
!! NAME
!!  exc_print_eig
!!
!! FUNCTION
!!  Print excitonic eigenvalues on std_out and ab_out.
!!
!! INPUTS
!!  gw_gap=GW direct gap.
!!  bseig_fname=The name of file containing eigenvalues and eigenvectors
!!
!! OUTPUT
!!  exc_gap=Excitonic direct gap.
!!  Additional info on the Excitonic spectrum are reported on standard output.
!!
!! PARENTS
!!      m_exc_diago
!!
!! CHILDREN
!!
!! SOURCE

subroutine exc_print_eig(BSp,bseig_fname,gw_gap,exc_gap)

!Arguments ------------------------------------
!scalars
 complex(dpc),intent(in) :: gw_gap
 complex(dpc),intent(out) :: exc_gap
 character(len=*),intent(in) :: bseig_fname
 type(excparam),intent(in) :: BSp

!Local variables ------------------------------
!scalars
 integer :: nstates_read,ii,j,k,eig_unt,ieig,hsize_exp
 integer :: hsize_read !,nstates
 complex(dpc) :: bind_energy,ctemp
 character(len=500) :: msg
 !type(Hdr_type) :: tmp_Hdr
!arrays
 integer,allocatable :: iperm(:)
 real(dp),allocatable :: exc_rene(:)
 complex(dpc),allocatable :: exc_cene(:)

!************************************************************************

 exc_gap = czero

 if (open_file(bseig_fname,msg,newunit=eig_unt,form="unformatted",status="old",action="read") /= 0) then
   ABI_ERROR(msg)
 end if

 read(eig_unt) ! do_ep_lifetime
 read(eig_unt) hsize_read, nstates_read

 if (BSp%use_coupling==0) hsize_exp =   SUM(Bsp%nreh)
 if (BSp%use_coupling>0)  hsize_exp = 2*SUM(Bsp%nreh)

 if (hsize_exp /= hsize_read) then
   write(msg,'(2(a,i0))')" Wrong dimension: read: ",hsize_read," expected= ",hsize_exp
   ABI_ERROR(msg)
 end if

 ABI_MALLOC(exc_cene,(nstates_read))
 read(eig_unt) exc_cene(:)

 ABI_MALLOC(exc_rene,(nstates_read))
 exc_rene = DBLE(exc_cene)

 ABI_MALLOC(iperm,(nstates_read))
 iperm = (/(ii, ii=1,nstates_read)/)

 call sort_dp(nstates_read,exc_rene,iperm,tol6)

 ABI_FREE(exc_rene)
 ABI_FREE(iperm)

 ! put in ascending order
 do ii=nstates_read,2,-1
   do j=1,ii-1
     if (DBLE(exc_cene(j)) > DBLE(exc_cene(j+1))) then
       ctemp = exc_cene(j)
       exc_cene(j) = exc_cene(j+1)
       exc_cene(j+1) = ctemp
     end if
   end do
 end do

 exc_gap = DCMPLX(ABS(DBLE(exc_cene(1))),AIMAG(exc_cene(1)))

 do ii=1,nstates_read
   if (ABS(DBLE(exc_cene(ii))) < DBLE(exc_gap)) then
     exc_gap = DCMPLX(ABS(DBLE(exc_cene(ii))),AIMAG(exc_cene(ii)))
   end if
 end do

 bind_energy = gw_gap - exc_gap

 write(msg,"(3(a,2f6.2,2a))")&
&  " GW  direct gap     ",gw_gap*Ha_eV,     " [eV] ",ch10,&
&  " EXC direct gap     ",exc_gap*Ha_eV,    " [eV] ",ch10,&
&  " EXC binding energy ",bind_energy*Ha_eV," [eV] ",ch10
 call wrtout([std_out, ab_out], msg)

 msg=' Excitonic eigenvalues up to the GW energy gap [eV]'
 call wrtout([std_out, ab_out], msg)

 do ii=1,nstates_read
   if (DBLE(exc_cene(ii)) > zero) EXIT
 end do

 do j=ii,nstates_read
   if (DBLE(exc_cene(j)) > DBLE(gw_gap)) EXIT
 end do
 j=j-1

 do ieig=ii,j
   write(msg,'(i3,a,2f6.2,a)')ieig," (",exc_cene(ieig)*Ha_eV,")"
   call wrtout([std_out, ab_out], msg)
 end do

 ii=ii-1
 do j=ii,1,-1
   if (ABS(DBLE(exc_cene(j))) > DBLE(gw_gap)) EXIT
 end do
 j=j+1

 ! This coding is not portable, write to ab_out has been disabled.
 if (ii>0) then
   do k=ii,j,-1
     write(msg,'(i3,a,2f6.2,a)')k," (",exc_cene(k)*Ha_eV,")"
     call wrtout(std_out, msg)
   end do
 end if

 ABI_FREE(exc_cene)

 close(eig_unt)

end subroutine exc_print_eig
!!***

!----------------------------------------------------------------------

!!****f* m_exc_diago/exc_diago_coupling
!! NAME
!!  exc_diago_coupling
!!
!! FUNCTION
!!  Calculate excitonic eigenvalues and eigenvectors by performing a direct diagonalization.
!!  of the non Hermitian excitonic Hamiltoninan (resonant + coupling).
!!
!! INPUTS
!!  bseig_fname=The name of the output file.
!!  Bsp
!!    neh=Rank of the resonant block of the Hamiltoninan (equal to the rank of the coupling part)
!!  comm=MPI communicator.
!!  BS_files<excfiles>=Datatype storing names and files used in the Bethe-Salpeter code.
!!
!! OUTPUT
!!  Excitonic eigenvectors and eigenvalues are written on file BS_files%out_eig.
!!
!! PARENTS
!!      m_exc_diago
!!
!! CHILDREN
!!
!! SOURCE

subroutine exc_diago_coupling(Bsp,BS_files,Hdr_bse,prtvol,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,prtvol
 type(excfiles),intent(in) :: BS_files
 type(excparam),intent(in) :: BSp
 type(Hdr_type),intent(in) :: Hdr_bse

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0,ldvl=1
 integer(i8b) :: bsize_ham
 integer :: ii,exc_size,hreso_unt,hcoup_unt,eig_unt,nsppol,nstates
 integer :: bsz,block,bs1,bs2,jj
 integer :: fform,row_sign
 integer :: mi,it,nprocs,my_rank !itp
 integer :: nene_printed,ierr
 real(dp) :: exc_gap,exc_maxene,temp
 logical :: diago_is_real,do_full_diago
 logical :: do_ep_lifetime
 character(len=500) :: msg
 character(len=fnlen) :: hreso_fname,hcoup_fname,bseig_fname
!arrays
 complex(dpc),allocatable :: exc_ham(:,:),exc_rvect(:,:),exc_ene(:),ovlp(:,:)
 complex(dpc),allocatable :: cbuff(:,:)
 complex(dpc) :: vl_dpc(ldvl,1)

!************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs  = xmpi_comm_size(comm)

 nsppol = Hdr_bse%nsppol
 if (nsppol==2) then
   ABI_WARNING("nsppol==2 with coupling is still under development")
 end if

 if (nprocs > 1) then
   ABI_WARNING("Scalapack does not provide ZGEEV, diagonalization is done in sequential!")
 end if

 exc_size = 2*SUM(BSp%nreh)
 nstates  = BSp%nstates
 do_full_diago = (exc_size==nstates)
 ABI_CHECK(do_full_diago,"Partial diago not coded yet")

 bseig_fname = BS_files%out_eig
 if (BS_files%in_eig /= BSE_NOFILE) then
   ABI_ERROR("BS_files%in_eig is defined!")
 end if
 !
 ! Only master performs the diagonalization since ScaLAPACK does not provide the parallel version of ZGEEV.
 if (my_rank/=master) GOTO 10

 write(msg,'(a,i0)')' Direct diagonalization of the full excitonic Hamiltonian, Matrix size= ',exc_size
 call wrtout([std_out, ab_out], msg)

 bsize_ham = 2*dpc*exc_size**2
 write(msg,'(a,f9.2,a)')' Allocating full excitonic Hamiltonian. Memory requested: ',bsize_ham*b2Gb,' Gb. '
 call wrtout(std_out, msg)

 ABI_MALLOC_OR_DIE(exc_ham,(exc_size,exc_size), ierr)

 write(msg,'(3a,f8.1,3a,f8.1,a)')&
  ' Allocating excitonic eigenvalues and eigenvectors. ',ch10,&
  ' Memory-space requested: ',2*dpc*exc_size*b2Gb,' Gb. ',ch10,&
  ' Memory-space requested: ',bsize_ham*b2Gb,' Gb. '
 call wrtout(std_out, msg)

 ABI_MALLOC_OR_DIE(exc_ene,(exc_size), ierr)

 if (BS_files%in_hreso /= BSE_NOFILE) then
   hreso_fname = BS_files%in_hreso
 else
   hreso_fname = BS_files%out_hreso
 end if

 call wrtout(std_out,' Reading resonant excitonic Hamiltonian from '//TRIM(hreso_fname))

 if (open_file(hreso_fname,msg,newunit=hreso_unt,form="unformatted",status="old",action="read") /= 0) then
   ABI_ERROR(msg)
 end if
 !
 ! Read the header and perform consistency checks.
 call exc_read_bshdr(hreso_unt,Bsp,fform,ierr)
 ABI_CHECK(ierr==0,"Wrong header")
 !
 ! Construct resonant and anti-resonant part of the excitonic Hamiltonian using Hermiticity. File is always in double precision.
 ! Fill exc_ham with ( R  0 )
 !                   ( 0 -R*)
!BEGINDEBUG
 exc_ham = HUGE(one)
!ENDDEBUG

 row_sign=-1; diago_is_real=(.not.BSp%have_complex_ene)
 call exc_fullh_from_blocks(hreso_unt,"Resonant",nsppol,row_sign,diago_is_real,BSp%nreh,exc_size,exc_ham)
 close(hreso_unt)

 if (BS_files%in_hcoup /= BSE_NOFILE) then
   hcoup_fname =  BS_files%in_hcoup
 else
   hcoup_fname =  BS_files%out_hcoup
 end if

 call wrtout(std_out,' Reading coupling excitonic Hamiltonian from '//TRIM(hcoup_fname))
 if (open_file(hcoup_fname,msg,newunit=hcoup_unt,form="unformatted",status="old",action="read") /= 0) then
   ABI_ERROR(msg)
 end if
 !
 ! Read the header and perform consistency checks.
 call exc_read_bshdr(hcoup_unt,Bsp,fform,ierr)
 ABI_CHECK(ierr==0,"Wrong header")
 !
 ! Fill exc_ham with ( 0  C) to have ( R   C )
 !                   (-C* 0)         (-C* -R*)
 row_sign=-1; diago_is_real=(.not.BSp%have_complex_ene) ! not used here
 call exc_fullh_from_blocks(hcoup_unt,"Coupling",nsppol,row_sign,diago_is_real,BSp%nreh,exc_size,exc_ham)

!BEGINDEBUG
 if (ANY(exc_ham==HUGE(one))) then
   write(msg,'(a,2(1x,i0))')"There is a bug in exc_fullh_from_blocks",COUNT(exc_ham==HUGE(one)),exc_size**2
   ABI_WARNING(msg)
   bsz = Bsp%nreh(1)
   ABI_MALLOC(cbuff,(bsz,bsz))
   block=0
   do jj=1,2*nsppol
     do ii=1,2*nsppol
       block=block+1
       bs1 = (ii-1)*bsz+1
       bs2 = (jj-1)*bsz+1
       cbuff = exc_ham(bs1:bs1+bsz-1,bs2:bs2+bsz-1)
       if (ANY(cbuff==HUGE(one))) then
         write(std_out,*)" for block ",ii,jj," found ",COUNT(cbuff==HUGE(one))," wrong entries"
       end if
     end do
   end do

   ABI_FREE(cbuff)
   ABI_ERROR("Cannot continue")
 end if
!ENDDEBUG

 close(hcoup_unt)
 !
 ! ======================================================
 ! ==== Calculate right eigenvectors and eigenvalues ====
 ! ======================================================
 ABI_MALLOC_OR_DIE(exc_rvect,(exc_size,exc_size), ierr)

 if (do_full_diago) then
   call wrtout(std_out,"Complete direct diagonalization with xgeev...")
   call xgeev("No_left_eigen","Vectors",exc_size,exc_ham,exc_size,exc_ene,vl_dpc,ldvl,exc_rvect,exc_size)
 else
   ABI_ERROR("Not implemented error")
 end if

 ABI_FREE(exc_ham)

 exc_gap    = MINVAL(ABS(DBLE (exc_ene(1:nstates))))
 exc_maxene = MAXVAL(ABS(DBLE (exc_ene(1:nstates))))
 temp       = MAXVAL(ABS(AIMAG(exc_ene(1:nstates))))

 write(msg,'(2(a,f7.2,2a),a,es9.2,2a)')&
  " First excitonic eigenvalue: ",exc_gap*Ha_eV,   " [eV].",ch10,&
  " Last  excitonic eigenvalue: ",exc_maxene*Ha_eV," [eV].",ch10,&
  " Largest imaginary part:     ",temp*Ha_eV,      " [eV] ",ch10
 call wrtout([std_out, ab_out], msg)

 nene_printed = MIN(32*nsppol,nstates); if (prtvol>10) nene_printed = nstates

 ! This is not portable as the the eigenvalues calculated by ZGEEV are not sorted.
 ! Even two subsequent calculations with the same input on the same machine
 ! might produce different orderings. Might sort the eigenvalues though, just for printing.

 write(msg,'(a,i0)')' Complex excitonic eigenvalues in eV up to n= ',nene_printed
 call wrtout(std_out,msg)

 do it=0,(nene_printed-1)/4
   write(msg,'(8f10.5)') ( exc_ene(ii)*Ha_eV, ii=1+it*4,MIN(it*4+4,nene_printed) )
   call wrtout(std_out,msg)
 end do

 call wrtout(std_out,ch10//" Writing eigenvalues and eigenvectors on file "//TRIM(bseig_fname))

 if (open_file(bseig_fname,msg,newunit=eig_unt,form="unformatted",action="write") /= 0) then
   ABI_ERROR(msg)
 end if

!YG : new version with lifetime
 do_ep_lifetime = .FALSE.
 write(eig_unt) do_ep_lifetime

 write(eig_unt)exc_size,nstates
 write(eig_unt)CMPLX(exc_ene(1:nstates),kind=dpc)
 do mi=1,nstates
   write(eig_unt) exc_rvect(:,mi)
 end do

 ABI_FREE(exc_ene)

 ABI_MALLOC_OR_DIE(ovlp,(nstates,nstates), ierr)

 call wrtout(std_out,' Calculating overlap matrix... ')

 !do itp=1,nstates
 !  do it=1,nstates
 !    ovlp(it,itp) = xdotc(exc_size,exc_rvect(:,it),1,exc_rvect(:,itp),1)
 !  end do
 !end do
 call xgemm("C","N",exc_size,nstates,nstates,cone,exc_rvect,exc_size,exc_rvect,exc_size,czero,ovlp,nstates)
 ABI_FREE(exc_rvect)

 call wrtout(std_out," Inverting overlap matrix... ")

 ! Version for generic complex matrix.
 !call xginv(ovlp,exc_size)

 ! The overlap is Hermitian definite positive.
 call xhdp_invert("Upper",ovlp,nstates)
 call hermitianize(ovlp,"Upper")

 call wrtout(std_out,' Writing overlap matrix S^-1 on file: '//TRIM(bseig_fname))

 do it=1,nstates
   write(eig_unt) CMPLX(ovlp(:,it),kind=dpc)
 end do

 ABI_FREE(ovlp)

 close(eig_unt)

10 call xmpi_barrier(comm)

end subroutine exc_diago_coupling
!!***

!----------------------------------------------------------------------

!!****f* m_exc_diago/exc_diago_coupling_hegv
!! NAME
!!  exc_diago_coupling_hegv
!!
!! FUNCTION
!!  Calculate excitonic eigenvalues and eigenvectors by performing a direct diagonalization.
!!  of the non Hermitian excitonic Hamiltonian (resonant + coupling).
!!
!! INPUTS
!!  bseig_fname=The name of the output file.
!!  Bsp
!!    neh=Rank of the resonant block of the Hamiltoninan (equal to the rank of the coupling part)
!!  comm=MPI communicator.
!!  BS_files<excfiles>=Datatype storing names and files used in the Bethe-Salpeter code.
!!
!! OUTPUT
!!  Excitonic eigenvectors and eigenvalues are written to file BS_files%out_eig.
!!
!! PARENTS
!!      m_exc_diago
!!
!! CHILDREN
!!
!! SOURCE

subroutine exc_diago_coupling_hegv(Bsp,BS_files,Hdr_bse,prtvol,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,prtvol
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(Hdr_type),intent(in) :: Hdr_bse

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer(i8b) :: bsize_ham
 integer :: itype,il,iu,spin,row1,row2,pad_r1,pad_r2,neh1,neh2
 integer :: ii,exc_size,hreso_unt,hcoup_unt,eig_unt
 integer :: fform,neig_found,nstates
 integer :: mi,it,nprocs,my_rank
 integer :: nene_printed,nsppol,row_sign,ierr
 real(dp) :: exc_gap,exc_maxene,abstol,vl,vu
 character(len=500) :: msg
 character(len=fnlen) :: reso_fname,coup_fname,bseig_fname
 logical :: use_scalapack,do_full_diago,diago_is_real
 logical :: do_ep_lifetime
!arrays
 real(dp),allocatable :: exc_ene(:) !,test_ene(:)
 complex(dpc),allocatable :: exc_ham(:,:),exc_rvect(:,:),fmat(:,:),ovlp(:,:)
#if defined HAVE_LINALG_SCALAPACK && defined HAVE_MPI_IO
 integer,parameter :: istwfk1=1
 integer :: amode,mpi_fh,tbloc,tmp_unt,mene_found,mpi_err,my_nel,nsblocks
 integer :: iloc,jj,jloc,iglob,jglob,etype,slk_mask_type,offset_err,el,rrs_kind,ccs_kind
 integer :: max_r,max_c
 integer(XMPI_OFFSET_KIND) :: ehdr_offset,fmarker,my_offset
 integer :: gsub(2,2)
 logical,parameter :: is_fortran_file=.TRUE.
 complex(dpc) :: ctmp
 integer,allocatable :: sub_block(:,:,:)
 integer,pointer :: myel2loc(:,:)
 complex(dpc),allocatable :: tmp_cbuffer(:)
 character(50) :: uplo
 real(dp),external :: PDLAMCH
 type(matrix_scalapack)    :: Slk_F,Slk_Hbar,Slk_vec,Slk_ovlp,Slk_tmp
 type(processor_scalapack) :: Slk_processor
#endif

!************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs  = xmpi_comm_size(comm)

 nsppol = Hdr_bse%nsppol
 if (nsppol==2) then
   ABI_WARNING("nsppol==2 is still under development!")
 end if

 neh1 = BSp%nreh(1); neh2=neh1
 if (nsppol==2) neh2 = BSp%nreh(2)

 exc_size = 2*SUM(Bsp%nreh)
 nstates  = Bsp%nstates
 do_full_diago=(nstates==exc_size)

 write(msg,'(a,i0)')'. Direct diagonalization of the full excitonic Hamiltonian, Matrix size= ',exc_size
 call wrtout([std_out, ab_out], msg)

 bseig_fname = BS_files%out_eig
 if (BS_files%in_eig /= BSE_NOFILE) then
   ABI_ERROR("BS_files%in_eig is defined!")
 end if

 if (BS_files%in_hreso /= BSE_NOFILE) then
   reso_fname = BS_files%in_hreso
 else
   reso_fname = BS_files%out_hreso
 end if
 call wrtout(std_out,' Reading resonant excitonic Hamiltonian from '//TRIM(reso_fname))

 if (BS_files%in_hcoup /= BSE_NOFILE) then
   coup_fname =  BS_files%in_hcoup
 else
   coup_fname =  BS_files%out_hcoup
 end if
 call wrtout(std_out,' Reading coupling excitonic Hamiltonian from '//TRIM(coup_fname))

 ! TODO: Reintegrate SCALAPACK: use new format
! --- !ERROR
! src_file: m_exc_diago.F90
! src_line: 1398
! mpi_rank: 0
! message: |
!     SET_VIEW
!     Other I/O error , error stack:
!     ADIO_Set_view(48):  **iobadoverlap displacements of filetype must be in a monotonically nondecreasing order
! ...

 use_scalapack = .FALSE.
!#ifdef HAVE_LINALG_SCALAPACK
! use_scalapack = (nprocs > 1)
!#endif
 !use_scalapack = .FALSE.
 !use_scalapack = .TRUE.

 if (.not.use_scalapack .and. my_rank/=master) GOTO 10

 ABI_MALLOC_OR_DIE(exc_ene,(exc_size), ierr)

 SELECT CASE (use_scalapack)

 CASE (.FALSE.)
   write(msg,'(a)')". Using LAPACK sequential version to solve FHv = ev with H positive definite. "
   call wrtout([std_out, ab_out], msg)

   bsize_ham = 2*dpc*exc_size**2
   write(msg,'(a,f9.2,a)')' Allocating full excitonic Hamiltonian. Memory requested: ',2*bsize_ham*b2Gb,' Gb. '
   call wrtout(std_out, msg)

   ABI_MALLOC_OR_DIE(exc_ham,(exc_size,exc_size), ierr)

   ABI_MALLOC_OR_DIE(fmat,(exc_size,exc_size), ierr)

   write(msg,'(3a,f8.1,3a,f8.1,a)')&
    ' Allocating excitonic eigenvalues and eigenvectors. ',ch10,&
    ' Memory-space requested: ',2*dpc*exc_size*b2Gb,' Gb. ',ch10,&
    ' Memory-space requested: ',bsize_ham*b2Gb,' Gb. '
   call wrtout(std_out, msg)

   if (open_file(reso_fname,msg,newunit=hreso_unt,form="unformatted",status="old",action="read") /= 0) then
     ABI_ERROR(msg)
   end if
   !
   ! Read the header and perform consistency checks.
   call exc_read_bshdr(hreso_unt,Bsp,fform,ierr)
   ABI_CHECK(ierr==0,"Wrong header")
   !
   ! Construct Hbar = ( R   C )
   !                  ( C*  R*)
   !
   row_sign=+1; diago_is_real=(.not.BSp%have_complex_ene)
   call exc_fullh_from_blocks(hreso_unt,"Resonant",nsppol,row_sign,diago_is_real,Bsp%nreh,exc_size,exc_ham)
   close(hreso_unt)

   if (open_file(coup_fname,msg,newunit=hcoup_unt,form="unformatted",status="old",action="read") /= 0) then
     ABI_ERROR(msg)
   end if
   !
   ! Read the header and perform consistency checks.
   call exc_read_bshdr(hcoup_unt,Bsp,fform,ierr)
   ABI_CHECK(ierr==0,"Wrong header")

   row_sign=+1; diago_is_real=(.not.BSp%have_complex_ene) ! not used here.
   call exc_fullh_from_blocks(hcoup_unt,"Coupling",nsppol,row_sign,diago_is_real,Bsp%nreh,exc_size,exc_ham)
   close(hcoup_unt)

#ifdef DEV_MG_DEBUG_THIS
write(666)exc_ham
#endif
   !
   ! Fill fmat = (1  0)
   !             (0 -1)
   fmat = czero
   do spin=1,nsppol
     pad_r1 = (spin-1)*Bsp%nreh(1)
     pad_r2 = SUM(Bsp%nreh)
     if (spin==2) pad_r2 = pad_r2 + Bsp%nreh(1)
     do it=1,Bsp%nreh(spin)
       row1 = it + pad_r1
       row2 = it + pad_r2
       fmat(row1,row1) =  cone
       fmat(row2,row2) = -cone
     end do
   end do
   !
   ! ==================================================
   ! ==== Solve generalized EV problem F H u = e u ====
   ! ==================================================
   ! The eigenvectors Z are normalized as follows: if ITYPE = 1 or 2, Z**T*B*Z = I; if ITYPE = 3, Z**T*inv(B)*Z = I.
   !
   itype=2
   if (do_full_diago) then
     call wrtout(std_out," Full diagonalization with XHEGV... ")
     call xhegv(itype,"Vectors","Upper",exc_size,fmat,exc_ham,exc_ene)
   else
     call wrtout(std_out," Partial diagonalization with XHEGVX... ")
     ABI_MALLOC_OR_DIE(exc_rvect,(exc_size,nstates), ierr)
     il=1; iu=1; abstol=zero
     call xhegvx(itype,"Vectors","All","Upper",exc_size,fmat,exc_ham,vl,vu,il,iu,abstol,neig_found,exc_ene,exc_rvect,exc_size)
   end if

   ABI_FREE(exc_ham)

   if (do_full_diago) then
     ABI_MALLOC(exc_rvect,(exc_size,nstates))
     exc_rvect = fmat(:,1:nstates)
   end if

   ABI_FREE(fmat)

   call wrtout(std_out," Writing eigenvalues and eigenvectors on file: "//TRIM(bseig_fname))

   if (open_file(bseig_fname,msg,newunit=eig_unt,form="unformatted",action="write") /= 0) then
     ABI_ERROR(msg)
   end if

   do_ep_lifetime = .FALSE.
   write(eig_unt) do_ep_lifetime
   write(eig_unt) exc_size, nstates
   write(eig_unt) CMPLX(exc_ene(1:nstates),kind=dpc)
   do mi=1,nstates
     write(eig_unt) CMPLX(exc_rvect(:,mi),kind=dpc)
   end do

#ifdef DEV_MG_DEBUG_THIS
   write(888)exc_rvect
   write(888)exc_ene
#endif

   ABI_MALLOC_OR_DIE(ovlp,(nstates,nstates), ierr)

   call wrtout(std_out,' Calculating overlap matrix...')

   call xgemm("C","N",exc_size,nstates,nstates,cone,exc_rvect,exc_size,exc_rvect,exc_size,czero,ovlp,nstates)
   ABI_FREE(exc_rvect)

#ifdef DEV_MG_DEBUG_THIS
write(667)ovlp
#endif

   call wrtout(std_out," Inverting overlap matrix... ")
   !
   ! The overlap is Hermitian definite positive.
   call xhdp_invert("Upper",ovlp,nstates)
   call hermitianize(ovlp,"Upper")

   ! Version for generic complex matrix.
   !call xginv(ovlp,nstates)

#ifdef DEV_MG_DEBUG_THIS
write(668,*)ovlp
#endif

   call wrtout(std_out,' Writing overlap matrix O^-1 on file: '//TRIM(bseig_fname))
   do it=1,nstates
     write(eig_unt) ovlp(:,it)
   end do

   ABI_FREE(ovlp)
   close(eig_unt)

 CASE (.TRUE.)

#if defined HAVE_LINALG_SCALAPACK && defined HAVE_MPI_IO
   !
   ! Init scaLAPACK matrix Hbar = ( R   C )
   !                              ( C*  R*)
   ! Battle plan:
   !   Here the reading is complicated by the fact that R and C are stored on two different files
   !   and moreover the matrices are in packed storage mode.
   !   For initializing the local part of the resonant and anti-resonant block we have to allocate
   !   a temporary buffer. Then we read the buffer from file and the corresponding elements of the
   !   scaLAPACK matrix are initialized taking into account the symmetries of R (Hermitian)
   !   The same procedure is used to read the coupling and the anti-coupling part (Symmetric).
   !
   tbloc=50
   write(msg,'(2(a,i0))')". Using MPI-IO + scaLAPACK version with nprocs= ",nprocs,"; block size= ",tbloc
   call wrtout([std_out, ab_out], msg, do_flush=.True.)
   !
   ! Init scaLAPACK environment.
   call init_scalapack(Slk_processor,comm)
   !
   ! Open the Resonant file with MPI-IO and skip the record.
   amode=MPI_MODE_RDONLY

   call MPI_FILE_OPEN(comm, reso_fname, amode, MPI_INFO_NULL, mpi_fh, mpi_err)
   msg = " MPI_IO error opening file: "//TRIM(reso_fname)
   ABI_CHECK_MPI(mpi_err,msg)
   !
   ! Skip the header and find the offset for reading the matrix.
   call exc_skip_bshdr_mpio(mpi_fh,xmpio_collective,ehdr_offset)
   !
   ! Read  = ( R  - )
   !         ( -  R*)
   call init_matrix_scalapack(Slk_Hbar,exc_size,exc_size,Slk_processor,istwfk1,tbloc=tbloc)

   nullify(myel2loc)
   nsblocks=nsppol
   ABI_MALLOC(sub_block,(2,2,nsblocks))
   ABI_CHECK(nsppol==1,"nsppol==2 not coded yet")

   call slk_single_fview_read_mask(Slk_Hbar,rrs_of_glob,offset_in_file,nsblocks,sub_block,my_nel,myel2loc,etype,slk_mask_type,&
&    offset_err,is_fortran_file)

   if (offset_err/=0) then
     write(msg,"(3a)")&
&      " Global position index cannot be stored in a standard Fortran integer ",ch10,&
&      " Excitonic matrix cannot be read with a single MPI-IO call."
     ABI_ERROR(msg)
   end if

   ! Shift the offset because the view starts at the fist matrix element!
   ! TODO should rationalize the treatment of the offset
   my_offset = ehdr_offset + xmpio_bsize_frm
   call MPI_FILE_SET_VIEW(mpi_fh, my_offset, etype, slk_mask_type, 'native', MPI_INFO_NULL, mpi_err)
   ABI_CHECK_MPI(mpi_err,"SET_VIEW")

   call MPI_TYPE_FREE(slk_mask_type,mpi_err)
   ABI_CHECK_MPI(mpi_err,"MPI_TYPE_FREE")
   !
   ! Read my portion of the R,-R* sublocks and store the values in a temporary buffer.
   ABI_MALLOC_OR_DIE(tmp_cbuffer,(my_nel), ierr)

   call xmpi_barrier(comm)

   call MPI_FILE_READ_ALL(mpi_fh, tmp_cbuffer, my_nel, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpi_err)
   ABI_CHECK_MPI(mpi_err,"READ_ALL")
   !
   ! Symmetrize my Resonant part.
   do el=1,my_nel
     iloc = myel2loc(1,el)
     jloc = myel2loc(2,el)
     call idx_glob(Slk_Hbar,iloc,jloc,iglob,jglob)
     ctmp = tmp_cbuffer(el)
     if (iglob==jglob.and..not.Bsp%have_complex_ene) ctmp = DBLE(ctmp) ! Force the diagonal to be real.
     rrs_kind = rrs_of_glob(iglob,jglob,Slk_Hbar%sizeb_global)
     if (rrs_kind==1.and.jglob<iglob) then ! Lower resonant
       ctmp = DCONJG(ctmp)
     else if (rrs_kind==-1.and.jglob>=iglob) then  ! Lower Anti-resonant (Diagonal is included).
       ctmp = DCONJG(ctmp)
     end if
     Slk_Hbar%buffer_cplx(iloc,jloc) = ctmp
   end do

   ABI_FREE(tmp_cbuffer)
   ABI_FREE(myel2loc)

   call MPI_FILE_CLOSE(mpi_fh, mpi_err)
   ABI_CHECK_MPI(mpi_err,"FILE_CLOSE")
   !
   ! Read  = ( -  C)
   !         (-C* -)
   !
   call MPI_FILE_OPEN(comm, coup_fname, amode, MPI_INFO_NULL, mpi_fh, mpi_err)
   msg = " MPI_IO error opening file: "//TRIM(coup_fname)
   ABI_CHECK_MPI(mpi_err,msg)
   !
   ! Skip the header and find the offset for reading the matrix.
   call exc_skip_bshdr_mpio(mpi_fh,xmpio_collective,ehdr_offset)

   nullify(myel2loc)
   call slk_single_fview_read_mask(Slk_Hbar,ccs_of_glob,offset_in_file,nsblocks,sub_block,my_nel,myel2loc,etype,slk_mask_type,&
&    offset_err,is_fortran_file)

   ABI_FREE(sub_block)

   if (offset_err/=0) then
     write(msg,"(3a)")&
&      " Global position index cannot be stored in a standard Fortran integer ",ch10,&
&      " Excitonic matrix cannot be read with a single MPI-IO call."
     ABI_ERROR(msg)
   end if
   !
   ! Shift the offset because the view starts at the fist matrix element!
   ! TODO should rationalize the treatment of the offset so that the client code
   ! will automatically receive my_offset.
   my_offset = ehdr_offset + xmpio_bsize_frm
   call MPI_FILE_SET_VIEW(mpi_fh, my_offset, etype, slk_mask_type, 'native', MPI_INFO_NULL, mpi_err)
   ABI_CHECK_MPI(mpi_err,"SET_VIEW")

   call MPI_TYPE_FREE(slk_mask_type,mpi_err)
   ABI_CHECK_MPI(mpi_err,"MPI_TYPE_FREE")
   !
   ! Read my portion of the C-C* blocks and store the values in a temporary buffer.
   ABI_MALLOC_OR_DIE(tmp_cbuffer,(my_nel), ierr)

   call MPI_FILE_READ_ALL(mpi_fh, tmp_cbuffer, my_nel, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpi_err)
   ABI_CHECK_MPI(mpi_err,"READ_ALL")
   !
   ! Symmetrize my coupling part.
   ! Coupling block is symmetric => No symmetrization of the lower triangle.
   do el=1,my_nel
     iloc = myel2loc(1,el)
     jloc = myel2loc(2,el)
     call idx_glob(Slk_Hbar,iloc,jloc,iglob,jglob)
     ccs_kind = ccs_of_glob(iglob,jglob,Slk_Hbar%sizeb_global)
     ctmp = tmp_cbuffer(el)
     if (ccs_kind==-1) ctmp = DCONJG(ctmp) ! Anti-coupling (Diagonal is included).
     Slk_Hbar%buffer_cplx(iloc,jloc) = ctmp
   end do

   ABI_FREE(tmp_cbuffer)
   ABI_FREE(myel2loc)

   !max_r=20; max_c=10
   !call print_arr(Slk_Hbar%buffer_cplx,max_r=max_r,max_c=max_c,unit=std_out)

#ifdef DEV_MG_DEBUG_THIS
   ABI_MALLOC(exc_ham,(exc_size,exc_size))
   read(666)exc_ham

   write(std_out,*)"Error Hbar: ",MAXVAL(ABS(exc_ham-Slk_Hbar%buffer_cplx))
   ABI_FREE(exc_ham)
#endif

   call MPI_FILE_CLOSE(mpi_fh, mpi_err)
   ABI_CHECK_MPI(mpi_err,"FILE_CLOSE")
   !
   ! Init scaLAPACK matrix F
   call init_matrix_scalapack(Slk_F,exc_size,exc_size,Slk_processor,istwfk1,tbloc=tbloc)
   !
   ! Global F = (1  0)
   !            (0 -1)
   do jloc=1,Slk_F%sizeb_local(2)
     do iloc=1,Slk_F%sizeb_local(1)
       call idx_glob(Slk_F,iloc,jloc,iglob,jglob)
       if (iglob==jglob) then
         if (iglob<=SUM(Bsp%nreh)) then
           Slk_F%buffer_cplx(iloc,jloc) =  cone
         else
           Slk_F%buffer_cplx(iloc,jloc) = -cone
         end if
       else
         Slk_F%buffer_cplx(iloc,jloc) =  czero
       end if
     end do
   end do
   !
   ! ===========================================================
   ! ==== Solve generalized EV problem H u = F Hbar u = e u ====
   ! ===========================================================
   call init_matrix_scalapack(Slk_vec,exc_size,exc_size,Slk_processor,istwfk1,tbloc=tbloc)
   !
   itype=2; vl=1; vu=1; il=1; iu=nstates
   abstol=zero !ABSTOL = PDLAMCH(comm,'U')

#if 1
   if (do_full_diago) then
     call slk_pzhegvx(itype,"Vectors","All","Upper",Slk_F,Slk_Hbar,vl,vu,il,iu,abstol,Slk_vec,mene_found,exc_ene)
   else
     ABI_WARNING("Partial diago is still under testing")
     call slk_pzhegvx(itype,"Vectors","Index","Upper",Slk_F,Slk_Hbar,vl,vu,il,iu,abstol,Slk_vec,mene_found,exc_ene)
   end if
#else
   call xhegv(itype,"Vectors","Upper",exc_size,Slk_F%buffer_cplx,Slk_Hbar%buffer_cplx,exc_ene)
   Slk_vec%buffer_cplx = Slk_F%buffer_cplx
#endif

#ifdef DEV_MG_DEBUG_THIS
   if (PRODUCT(Slk_Hbar%sizeb_local) /= exc_size**2) then
     ABI_ERROR("Wrong size")
   end if

   ABI_MALLOC(exc_ham,(exc_size,exc_size))
   read(888)exc_ham

   write(std_out,*)"Error rvec: ",MAXVAL(ABS(exc_ham-Slk_vec%buffer_cplx))
   ABI_FREE(exc_ham)

   ABI_MALLOC(test_ene,(exc_size))
   read(888)test_ene
   write(std_out,*)"Error ene: ",MAXVAL(ABS(exc_ene-test_ene))
   ABI_FREE(test_ene)
#endif

   call Slk_F%free()
   call Slk_Hbar%free()

   call wrtout(std_out,ch10//" Writing eigenvalues and eigenvectors on file: "//TRIM(bseig_fname))
   !
   ! Open the file with Fortran-IO to write the Header.
   if (my_rank==master) then
     if (open_file(bseig_fname,msg,newunit=eig_unt,form="unformatted",action="write") /= 0) then
       ABI_ERROR(msg)
     end if

     write(eig_unt) exc_size,nstates
     write(eig_unt) CMPLX(exc_ene(1:nstates),kind=dpc)
     !do mi=1,exc_size
     !  write(eig_unt) CMPLX(exc_rvect(:,mi),kind=dpc)
     !end do
     close(eig_unt)
   end if
   !
   ! Open the file with MPI-IO and write the distributed eigevectors.
   call xmpi_barrier(comm)
   amode=MPI_MODE_RDWR
   call MPI_FILE_OPEN(comm, bseig_fname, amode, MPI_INFO_NULL, mpi_fh, mpi_err)
   ABI_CHECK_MPI(mpi_err,"FILE_OPEN: "//TRIM(bseig_fname))
   !
   ! Skip the header and find the offset for writing the matrix.
   ehdr_offset=0
   !call hdr_mpio_skip(mpi_fh,fform,ehdr_offset)

   call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_collective,fmarker,mpi_err)
   write(std_out,*)" fmarker1 = ",fmarker
   call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_collective,fmarker,mpi_err)
   write(std_out,*)" fmarker2 = ",fmarker

   write(std_out,*)" Writing nstates ",nstates
   gsub(:,1) = (/1,1/)
   gsub(:,2) = (/exc_size,nstates/)
   call slk_write(Slk_vec,"All",is_fortran_file,mpi_fh=mpi_fh,offset=ehdr_offset,glob_subarray=gsub)

   call wrtout(std_out,' Calculating overlap matrix... ')
   if (.not.do_full_diago) then
     ABI_ERROR(" Init of Slk_ovlp is wrong")
   end if

   call init_matrix_scalapack(Slk_ovlp,exc_size,exc_size,Slk_processor,istwfk1,tbloc=tbloc)

   ! Calculate the overlap matrix.
   ! FIXME
   ! The ESLL manual says that "matrices matrix1 and matrix2 must have no common elements; otherwise, results are unpredictable."
   ! However the official scaLAPACK documentation does not report this (severe) limitation.

   !call init_matrix_scalapack(Slk_tmp,exc_size,exc_size,Slk_processor,istwfk1,tbloc=tbloc)
   !Slk_tmp%buffer_cplx = Slk_vec%buffer_cplx
   !call slk_pzgemm("C","N",Slk_tmp,cone,Slk_vec,czero,Slk_ovlp)
   !call Slk_tmp%free()

   call slk_pzgemm("C","N",Slk_vec,cone,Slk_vec,czero,Slk_ovlp)

#ifdef DEV_MG_DEBUG_THIS
   ABI_MALLOC(exc_ham,(exc_size,exc_size))
   read(667)exc_ham

   write(std_out,*)"Error Ovlp: ",MAXVAL(ABS(exc_ham-Slk_ovlp%buffer_cplx))
   !Slk_ovlp%buffer_cplx = exc_ham
#endif

   !max_r=20; max_c=10
   !call print_arr(Slk_ovlp%buffer_cplx,max_r=max_r,max_c=max_c,unit=std_out)

   call Slk_vec%free()

   call wrtout(std_out," Inverting overlap matrix... ")
   uplo="Upper"

#if 0
!DEBUG
   call xhdp_invert(uplo,Slk_ovlp%buffer_cplx,exc_size)

   !call slk_symmetrize(Slk_ovlp,uplo,"Hermitian")
   call hermitianize(Slk_ovlp%buffer_cplx,uplo)

   exc_ham = MATMUL(exc_ham,Slk_ovlp%buffer_cplx)
   do it=1,exc_size
     exc_ham(it,it) = exc_ham(it,it) - cone
   end do

   write(std_out,*)"Error Inversion: ",MAXVAL(ABS(exc_ham))
   ABI_FREE(exc_ham)
!END DEBUG

#else
   ! call slk_zdhp_invert(Slk_ovlp,uplo)
   ! call hermitianize(Slk_ovlp%buffer_cplx,uplo)
   ! !call slk_symmetrize(Slk_ovlp,uplo,"Hermitian")

   call slk_zinvert(Slk_ovlp)  ! Version for generic complex matrix.
#endif

   if (allocated(exc_ham))  then
     ABI_FREE(exc_ham)
   end if

#ifdef DEV_MG_DEBUG_THIS
   ABI_MALLOC(exc_ham,(exc_size,exc_size))
   read(668)exc_ham
   write(std_out,*)"Error in Inv Ovlp: ",MAXVAL(ABS(exc_ham-Slk_ovlp%buffer_cplx))

   !exc_ham = exc_ham-Slk_ovlp%buffer_cplx
   !do it=1,exc_size
   !  if ( MAXVAL(ABS(exc_ham(:,it))) > 0.1 ) write(std_out,*)"it: ",it,exc_ham(:,it)
   !end do

   !Slk_ovlp%buffer_cplx = exc_ham
   ABI_FREE(exc_ham)

   !write(std_out,*)"MAX ERR",MAXVAL(ABS(Slk_ovlp%buffer_cplx - TRANSPOSE(DCONJG(Slk_ovlp%buffer_cplx))))
#endif

   call wrtout(std_out,' Writing overlap matrix S^-1 on file: '//TRIM(bseig_fname))

   call slk_write(Slk_ovlp,"All",is_fortran_file,mpi_fh=mpi_fh,offset=ehdr_offset)

   call MPI_FILE_CLOSE(mpi_fh, mpi_err)
   ABI_CHECK_MPI(mpi_err,"FILE_CLOSE")

   call Slk_ovlp%free()
   call end_scalapack(Slk_processor)
#else
   ABI_BUG("You should not be here!")
#endif

 END SELECT

 exc_gap    = MINVAL(ABS(exc_ene(1:nstates)))
 exc_maxene = MAXVAL(ABS(exc_ene(1:nstates)))

 write(msg,'(2(a,f7.2,2a))')&
  " First excitonic eigenvalue: ",exc_gap*Ha_eV,   " [eV].",ch10,&
  " Last  excitonic eigenvalue: ",exc_maxene*Ha_eV," [eV].",ch10
 call wrtout([std_out, ab_out], msg)

 nene_printed = MIN(32*nsppol,nstates); if (prtvol>10) nene_printed = nstates
 write(msg,'(a,i0)')' Complex excitonic eigenvalues in eV up to n= ',nene_printed
 call wrtout(std_out, msg)

 do it=0,(nene_printed-1)/4
   write(msg,'(4f10.5)') ( exc_ene(ii)*Ha_eV, ii=1+it*4,MIN(it*4+4,nene_printed) )
   call wrtout(std_out, msg)
 end do

 ABI_FREE(exc_ene)

10 call xmpi_barrier(comm)

end subroutine exc_diago_coupling_hegv
!!***

!----------------------------------------------------------------------

end module m_exc_diago
!!***
