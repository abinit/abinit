!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_hexc
!! NAME
!! m_hexc
!!
!! FUNCTION
!! module for excitonic hamiltonian for Haydock
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (M.Giantomassi, Y. Gillet)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_hexc

 use defs_basis
 use m_abicore
 use m_bs_defs
 use m_xmpi
 use m_errors
 use m_nctk
 use m_haydock_io
 use m_linalg_interfaces
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_time,              only : timab
 use m_fstrings,          only : indent, strcat, sjoin, itoa
 use defs_datatypes,      only : ebands_t, pseudopotential_type
 use m_hide_blas,         only : xdotc, xgemv
 use m_numeric_tools,     only : print_arr, symmetrize, hermitianize, wrap2_pmhalf
 use m_crystal,           only : crystal_t
 use m_bz_mesh,           only : kmesh_t, findqg0, get_bz_item
 use m_double_grid,       only : double_grid_t, get_kpt_from_indices_coarse, compute_corresp
 use m_wfd,               only : wfd_t
 use m_bse_io,            only : exc_read_rcblock, exc_write_optme, exc_ham_ncwrite
 use m_pawtab,            only : pawtab_type
 use m_vcoul,             only : vcoul_t
 use m_bseinterp,         only : interpolator_t, interpolator_init, interpolator_normalize, &
&                                interpolator_free, int_alloc_work, int_free_work

 implicit none

 private
!!***

!!****t* m_haydock/hexc_t
!! NAME
!! hexc_t
!!
!! FUNCTION
!!  Store the excitonic hamiltonian and other related information
!!
!! SOURCE

 type,public :: hexc_t

 !scalars
    integer :: comm
    ! MPI communicator

    integer :: hsize_coarse
    ! Size of the coarse hamiltonian

    integer :: hsize
    ! Size of the hamiltonian and the kets
    ! (= hsize_coarse without interpolation, = hsize_dense with interpolation)

    integer :: nbz
    ! Number of kpoints for the full problem
    ! (= nbz_coarse without interpolation, = nbz_dense with interpolation)

    integer :: my_t1
    ! Lower limit of MPI paral

    integer :: my_t2
    ! Upper limit of MPI paral

    integer :: my_nt
    ! Number of transitions treat by node
    ! = my_t2 - my_t1 + 1

    integer :: nbnd_coarse
    ! Product of number of bands conduction X valence

    ! Pointers to data that are already in memory
    type(excparam),pointer :: bsp => null()
    ! parameters for BS

    type(excfiles),pointer :: bs_files => null()
    ! files for BSE

    type(crystal_t),pointer :: crystal => null()
    ! crystal info

    type(kmesh_t),pointer :: kmesh_coarse => null()
    ! kmesh of the coarse mesh

    type(kmesh_t),pointer :: kmesh => null()
    ! kmesh of the full problem

    type(ebands_t),pointer :: ks_bst => null()
    type(ebands_t),pointer :: qp_bst => null()
    ! band structures of the full problem

    type(wfd_t),pointer :: wfd_coarse => null()
    ! Wfd of the coarse problem

    type(wfd_t),pointer :: wfd => null()
    ! wfd of the full problem

 !arrays
    complex(dpc),allocatable :: hreso(:,:)
    ! Resonant part of the hamiltonian

    complex(dpc),allocatable :: hcoup(:,:)
    ! Coupling part of the hamiltonian

    complex(dpc),allocatable :: diag_coarse(:)
    ! Diagonal part of the hamiltonian with transition energies

 end type hexc_t
!!***

!----------------------------------------------------------------------

!!****t* m_hexc/hexc_interp_t
!! NAME
!! hexc_interp_t
!!
!! FUNCTION
!!  Store information about interpolation of excitonic hamiltonian
!!
!! SOURCE
 type,public :: hexc_interp_t

 !scalars
    integer :: hsize_dense
    ! Size of the dense hamiltonian

    real(dp) :: m3_width
    ! Width of the region where M3 is applied instead of M1

    type(interpolator_t) :: interpolator
    ! Interpolator containing overlaps and interpolation info

    ! Pointers to datatypes that are already in memory
    type(kmesh_t),pointer :: kmesh_dense => null()
    ! kmesh of the dense mesh

    type(vcoul_t),pointer :: vcp_dense => null()
    ! coulomb interaction on the dense mesh

 !arrays
    integer,allocatable :: kdense2div(:)
    ! kdense2div(nbz_dense)
    ! Index of kpoint -> Index of division

    integer,allocatable :: div2kdense(:,:)
    ! div2kdense(nbz_coarse,ndiv)
    ! Index of kpoint coarse + Index of division -> Index of kdense

    complex(dpc),allocatable :: diag_dense(:)
    ! diag_dense(hsize_dense)
    ! Diagonal part of the dense hamiltonian

    complex(dpc),allocatable :: hinterp(:,:)
    ! hinterp(hsize_dense,hsize_dense)
    ! Interpolated hamiltonian

    complex(dpc),allocatable :: all_hmat(:,:)
    ! all_hmat,(hsize,hsize))
    ! Coarse excitonic matrix in a format suitable for interpolation in k-space

    complex(dpc),allocatable :: all_acoeffs(:,:)
    ! all_acoeffs(hsize,hsize))
    ! a coefficients in a format suitable for interpolation in k-space

    complex(dpc),allocatable :: all_bcoeffs(:,:)
    ! all_bcoeffs(hsize,hsize))
    ! b coefficients in a format suitable for interpolation in k-space

    complex(dpc),allocatable :: all_ccoeffs(:,:)
    ! all_ccoeffs(hsize,hsize))
    ! c coefficients in a format suitable for interpolation in k-space

 end type hexc_interp_t

!!***

 public :: hexc_init           ! Construct the object
 public :: hexc_interp_init    ! Construct the object for interpolated ham
 public :: hexc_free           ! Free memory
 public :: hexc_interp_free    ! Free memory for interpolated ham
 public :: hexc_build_hinterp  ! Interpolate the Hamiltonian and store it in memory
 public :: hexc_matmul_tda     ! Matrix-vector multiplication (TDA)
 public :: hexc_matmul_full    ! Matrix-vector multiplication (TDA + Coupling)
 public :: hexc_matmul_elphon  ! Matrix-vector multiplication (TDA + elphon)

!----------------------------------------------------------------------

CONTAINS  !=======================================================================
!!***

!!****f* m_hexc/hexc_init
!! NAME
!! hexc_init
!!
!! FUNCTION
!! Construct the hexc object
!!
!! INPUTS
!! BSp<excparam>=Parameters of BS
!! BS_files<excparam>=Files for BS
!! Cryst<crystal_t>=Info on the crystalline structure
!! Kmesh_coarse<kmesh_t>=Kmesh info
!! Wfd_coarse<wfd_t>=Wavefunction descriptor
!! KS_BSt<ebands_t>=Kohn-Sham band structure
!! QP_BSt<ebands_t>=Quasi-Particle band structure
!! comm=communicator
!!
!! OUTPUT
!! hexc<hexc_t>=Excitonic Hamiltonian
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

subroutine hexc_init(hexc, BSp, BS_files, Cryst, Kmesh_coarse, Wfd_coarse, KS_BSt, QP_BSt, comm)

!Arguments ---------------------------
!scalars
 integer,intent(in) :: comm
 type(hexc_t),intent(inout) :: hexc
 type(excparam),intent(in),target :: BSp
 type(excfiles),intent(in),target :: BS_files
 type(crystal_t),intent(in),target :: Cryst
 type(kmesh_t),intent(in),target :: Kmesh_coarse
 type(wfd_t),intent(in),target :: Wfd_coarse
 type(ebands_t),intent(in),target :: KS_BSt, QP_BSt
!arrays

!Local variables ---------------------
!scalars
 integer :: ierr,ncid,ncerr
 integer :: max_r, max_c
 integer :: hsize
 integer :: spin, spad, itt ! For diagonal !
 logical :: is_resonant, diago_is_real, use_mpio=.FALSE.
 character(len=fnlen) :: hreso_fname, hcoup_fname
 !character(len=500) :: msg
!arrays
 complex(dpc),allocatable :: test(:,:)

!*****************************************************************************

 hexc%bsp => BSp
 hexc%bs_files => BS_files
 hexc%crystal => Cryst
 hexc%kmesh_coarse => Kmesh_coarse


 hexc%comm = comm
 hsize = SUM(BSp%nreh)

 hexc%ks_bst => KS_BSt
 hexc%qp_bst => QP_BSt
 hexc%wfd => Wfd_coarse
 hexc%wfd_coarse => Wfd_coarse
 hexc%kmesh => Kmesh_coarse
 hexc%hsize_coarse = hsize
 hexc%hsize = hsize
 hexc%nbz = Kmesh_coarse%nbz

 hexc%nbnd_coarse = BSp%maxnbndv*BSp%maxnbndc

 ! Divide the columns of the Hamiltonian among the nodes.
 call xmpi_split_work(hsize,comm,hexc%my_t1,hexc%my_t2)

 hexc%my_nt = hexc%my_t2 - hexc%my_t1 + 1
 ABI_CHECK(hexc%my_nt>0,"found processor with 0 rows")

 ABI_MALLOC_OR_DIE(hexc%hreso,(hsize,hexc%my_t1:hexc%my_t2), ierr)

 ! Read the resonant block from file.
 if (BS_files%in_hreso /= BSE_NOFILE) then
   hreso_fname = BS_files%in_hreso
 else
   hreso_fname = BS_files%out_hreso
 end if

 is_resonant=.TRUE.; diago_is_real=(.not.BSp%have_complex_ene)
 call exc_read_rcblock(hreso_fname,Bsp,is_resonant,diago_is_real,BSp%nsppol,BSp%nreh,hsize,&
&   hexc%my_t1,hexc%my_t2,hexc%hreso,use_mpio,comm)

 !BEGIN DEBUG
 if (use_mpio) then
   MSG_WARNING("Testing MPI-IO routines")
   ABI_MALLOC_OR_DIE(test,(hsize,hexc%my_t1:hexc%my_t2), ierr)
   diago_is_real=(.not.BSp%have_complex_ene)
   call exc_read_rcblock(hreso_fname,Bsp,is_resonant,diago_is_real,Bsp%nsppol,Bsp%nreh,hsize,&
&     hexc%my_t1,hexc%my_t2,test,.FALSE.,comm)
   test = test-hexc%hreso
   write(std_out,*)"DEBUG: Diff MPI-IO - Fortran ",MAXVAL(ABS(test))
   max_r=20; max_c=10
   write(std_out,*)" **** Testing resonant block **** "
   call print_arr(test,max_r=max_r,max_c=max_c,unit=std_out)
   if (BSp%nsppol==2) then
     write(std_out,*)" **** D down down ****"
     call print_arr(test(hsize/2+1:,hsize/2+1:),max_r=max_r,max_c=max_c,unit=std_out)
     write(std_out,*)" **** V up down ****"
     call print_arr(test(1:hsize/2,hsize/2+1:),max_r=max_r,max_c=max_c,unit=std_out)
     write(std_out,*)" **** V down up ****"
     call print_arr(test(hsize/2+1:,1:hsize/2),max_r=max_r,max_c=max_c,unit=std_out)
   end if
   ABI_FREE(test)
 end if
 !END DEBUG

 !
 ! Read coupling block.
 if (BSp%use_coupling>0) then
   ABI_CHECK(.not. Bsp%use_interp,"interpolation with coupling not coded!")
   if (BS_files%in_hcoup /= BSE_NOFILE) then
     hcoup_fname = BS_files%in_hcoup
   else
     hcoup_fname = BS_files%out_hcoup
   end if

   ABI_MALLOC_OR_DIE(hexc%hcoup,(hsize,hexc%my_t1:hexc%my_t2), ierr)
   is_resonant=.FALSE.; diago_is_real=.FALSE.
   call exc_read_rcblock(hcoup_fname,Bsp,is_resonant,diago_is_real,BSp%nsppol,BSp%nreh,hsize,&
&     hexc%my_t1,hexc%my_t2,hexc%hcoup,use_mpio,comm)
   !call symmetrize(hcoup,"ALL")

   if (use_mpio) then
     MSG_WARNING("Testing MPI-IO routines")
     ABI_MALLOC_OR_DIE(test,(hsize,hexc%my_t1:hexc%my_t2), ierr)
     diago_is_real=.FALSE.
     call exc_read_rcblock(hcoup_fname,Bsp,is_resonant,diago_is_real,BSp%nsppol,Bsp%nreh,hsize,&
&       hexc%my_t1,hexc%my_t2,test,.FALSE.,comm)
     test = test-hexc%hcoup
     write(std_out,*)"DEBUG: Diff MPI-IO - Fortran ",MAXVAL(ABS(test))
     max_r=20; max_c=10
     write(std_out,*)" **** Testing coupling block **** "
     call print_arr(test,max_r=max_r,max_c=max_c,unit=std_out)
     if (BSp%nsppol==2) then
       write(std_out,*)" **** D down down ****"
       call print_arr(test(hsize/2+1:,hsize/2+1:),max_r=max_r,max_c=max_c,unit=std_out)
       write(std_out,*)" **** V up down ****"
       call print_arr(test(1:hsize/2,hsize/2+1:),max_r=max_r,max_c=max_c,unit=std_out)
       write(std_out,*)" **** V down up ****"
       call print_arr(test(hsize/2+1:,1:hsize/2),max_r=max_r,max_c=max_c,unit=std_out)
     end if
     ABI_FREE(test)
   end if
 end if

 if(BSp%prt_ncham .or. BSp%use_interp) then
   ! I want to store the diagonal part for future use (printing or interpolation) !
   ABI_MALLOC(hexc%diag_coarse,(hexc%hsize_coarse))
   spad=0
   do spin=1,BSp%nsppol
     if(spin==2) spad=BSp%nreh(1)
     do itt=1,BSp%nreh(spin) ! 1 is for spin 1
       hexc%diag_coarse(spad+itt) = Bsp%Trans(itt,spin)%en
     end do
   end do

   if (BSp%prt_ncham) then
#ifdef HAVE_NETCDF
     ncerr = nctk_open_create(ncid, trim(hexc%BS_files%out_basename)//"_HEXC.nc", xmpi_comm_self)
     NCF_CHECK_MSG(ncerr, "Creating HEXC file")
     call exc_ham_ncwrite(ncid, hexc%Kmesh_coarse, hexc%BSp, hexc%hsize_coarse, hexc%BSp%nreh, &
&         hexc%BSp%vcks2t,hexc%hreso,hexc%diag_coarse)
     NCF_CHECK(nf90_close(ncid))
#else
     ABI_UNUSED(ncid)
#endif
   end if
 end if

end subroutine hexc_init
!!***

!-------------------------------------------------------------------

!!****f* m_hexc/hexc_interp_init
!! NAME
!! hexc_interp_init
!!
!! FUNCTION
!! Construct the hexc_interp object
!!
!! INPUTS
!! hexc<hexc_t>=Excitonic hamiltonian
!! Kmesh_dense<kmesh_t>=Kmesh info
!! Vcp_dense<vcoul_t>=Dense mesh info about coulomb
!! double_grid<double_grid_t>=Link between dense and coarse mesh
!! Wfd_dense<wfd_t>=Wavefunction descriptor
!! KS_BSt_dense<ebands_t>=Kohn-Sham band structure
!! QP_BSt_dense<ebands_t>=Quasi-Particle band structure
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials.
!! Pawtab(Cryst%ntypat*usepaw)<pawtab_type>=PAW tabulated starting data.
!! comm=communicator
!!
!! OUTPUT
!! hexc_i<hexc_interp_t>=Interpolated excitonic hamiltonian
!!
!! SIDE EFFECTS
!! hexc
!!   Will be modified so that the size of the problem is full interpolated hamiltonian
!! Wfd, Wfd_dense
!!   The memory might be modified by computing wavefunctions
!!
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

subroutine hexc_interp_init(hexc_i, hexc, m3_width, method, Kmesh_dense, Vcp_dense, &
&    double_grid, Wfd_dense, KS_BSt_dense, QP_BSt_dense, Psps, Pawtab)

!Arguments ---------------------------
!scalars
 integer,intent(in) :: method
 real(dp),intent(in) :: m3_width
 type(hexc_t),intent(inout) :: hexc
 type(hexc_interp_t),intent(inout) :: hexc_i
 type(double_grid_t),intent(in),target :: double_grid
 type(wfd_t),intent(inout),target :: Wfd_dense !, Wfd
 type(kmesh_t),intent(in),target :: Kmesh_dense
 type(pseudopotential_type),intent(in) :: Psps
 type(vcoul_t),intent(in),target :: Vcp_dense
 type(ebands_t),intent(in),target :: KS_BSt_dense, QP_BSt_dense
!arrays
 type(pawtab_type),intent(in) :: Pawtab(hexc%crystal%ntypat*hexc%Wfd_coarse%usepaw)

!Local variables ---------------------
!scalars
 integer,parameter :: spin1 = 1
 integer :: nsppol,ierr,ii,itt,nproc, my_rank,hsize
 real(dp),parameter :: threshold = 0.1_dp
 type(excparam) :: BSp
 logical :: is_resonant, diago_is_real, use_mpio
!arrays
 character(len=fnlen) :: tmpfname, hreso_fname

!*****************************************************************************

 BSp = hexc%bsp
 ABI_CHECK(BSp%nsppol == 1,"nsppol > 1 not implemented yet")

 hsize = hexc%hsize_coarse

 nproc  = xmpi_comm_size(hexc%comm); my_rank= xmpi_comm_rank(hexc%comm)
 nsppol = hexc%Bsp%nsppol

 ABI_CHECK(nproc == 1,"Parallelization not available in interpolation")

 hexc_i%m3_width = m3_width

 hexc_i%kmesh_dense => Kmesh_dense
 hexc_i%vcp_dense => Vcp_dense
 hexc_i%hsize_dense = SUM(BSp%nreh_interp)
 hexc%hsize = hexc_i%hsize_dense
 hexc%nbz = Kmesh_dense%nbz

 hexc%ks_bst => KS_BSt_dense
 hexc%qp_bst => QP_BSt_dense
 hexc%wfd => Wfd_dense
 hexc%kmesh => Kmesh_dense

 ! No parallelization !
 hexc%my_t1 = 1
 hexc%my_t2 = hexc_i%hsize_dense

 ! Initialize the interpolator
 call interpolator_init(hexc_i%interpolator, double_grid, Wfd_dense, hexc%Wfd_coarse, Kmesh_dense, &
&    hexc%Kmesh_coarse, hexc%BSp, hexc%crystal, Psps, Pawtab, method)

 if(BSp%sum_overlaps) then
   call interpolator_normalize(hexc_i%interpolator)
 end if

 ABI_MALLOC(hexc_i%kdense2div,(double_grid%nbz_dense))
 ABI_MALLOC(hexc_i%div2kdense,(double_grid%nbz_coarse,double_grid%ndiv))

 call compute_corresp(double_grid,hexc_i%div2kdense,hexc_i%kdense2div)

 if (any(BSp%interp_mode == [2,3])) then
   ! Read a, b, c coefficient matrices from file.
   ! For the time being, we read the full matrix in a temporary array, and
   ! then we store the data in a form suitable for the interpolation.
   is_resonant=.TRUE.; diago_is_real=(.not.BSp%have_complex_ene); use_mpio=.FALSE.

   if (hexc%bs_files%in_hreso /= BSE_NOFILE) then
     hreso_fname = hexc%bs_files%in_hreso
   else
     hreso_fname = hexc%bs_files%out_hreso
   end if

   tmpfname = hreso_fname; ii = LEN_TRIM(hreso_fname)

   ! TODO: Write new IO routines to read MPI-distributed data in a format suitable for the interpolation
   tmpfname(ii-2:ii+1) = 'ABSR'
   ABI_MALLOC_OR_DIE(hexc_i%all_acoeffs,(hsize,hsize), ierr)
   call exc_read_rcblock(tmpfname,Bsp,is_resonant,diago_is_real,nsppol,BSp%nreh,hsize,1,hsize,&
&     hexc_i%all_acoeffs,use_mpio,hexc%comm)

   tmpfname(ii-2:ii+1) = 'BBSR'
   ABI_MALLOC_OR_DIE(hexc_i%all_bcoeffs,(hsize,hsize), ierr)
   call exc_read_rcblock(tmpfname,Bsp,is_resonant,diago_is_real,nsppol,BSp%nreh,hsize,1,hsize,&
&     hexc_i%all_bcoeffs,use_mpio,hexc%comm)

   tmpfname(ii-2:ii+1) = 'CBSR'
   ABI_MALLOC_OR_DIE(hexc_i%all_ccoeffs,(hsize,hsize), ierr)
   call exc_read_rcblock(tmpfname,Bsp,is_resonant,diago_is_real,nsppol,BSp%nreh,hsize,1,hsize,&
&     hexc_i%all_ccoeffs,use_mpio,hexc%comm)
 end if

 ! Compute overlaps & compute all hmat
 ABI_MALLOC_OR_DIE(hexc_i%all_hmat,(hsize,hsize), ierr)

 hexc_i%all_hmat(:,:) = hexc%hreso(:,:)

 do itt=1,hsize
   hexc_i%all_hmat(itt,itt) = hexc_i%all_hmat(itt,itt) - hexc%diag_coarse(itt)
 end do

 ! I don't need the diag_coarse any more
 ABI_FREE(hexc%diag_coarse)

 ! Compute diagonal part of the dense Ham
 ABI_MALLOC(hexc_i%diag_dense,(hexc_i%hsize_dense))
 do itt=1,BSp%nreh_interp(spin1) ! 1 is for spin 1
   hexc_i%diag_dense(itt) = Bsp%Trans_interp(itt,spin1)%en
 end do

end subroutine hexc_interp_init
!!***

!-------------------------------------------------------------------

!!****f* m_hexc/hexc_build_hinterp
!! NAME
!! hexc_build_hinterp
!!
!! FUNCTION
!! Pre-compute interpolated hamiltonian and store it in memory
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!! hexc, hexc_i
!!   Pre-compute info to save CPU time when computing matmul
!!
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

subroutine hexc_build_hinterp(hexc,hexc_i)

!Arguments ---------------------------
 type(hexc_t),intent(inout) :: hexc
 type(hexc_interp_t),intent(inout) :: hexc_i

!Local variables ---------------------
 integer :: ierr,ncerr,ncid
 character(len=500) :: msg

!*****************************************************************************

 write(msg,"(a,f8.1,a)")"Memory needed for hinterp = ",one*(hexc_i%hsize_dense**2)*2*dpc*b2Mb," Mb"
 call wrtout(std_out,msg,"COLL")

 ABI_MALLOC_OR_DIE(hexc_i%hinterp,(hexc_i%hsize_dense,hexc_i%hsize_dense), ierr)

 call hexc_compute_hinterp(hexc%BSp, hexc%hsize_coarse, hexc_i%hsize_dense, hexc_i%all_hmat, &
&  hexc_i%interpolator%double_grid,hexc%nbnd_coarse, hexc_i%interpolator, &
&  hexc_i%kdense2div, hexc_i%all_acoeffs,hexc_i%all_bcoeffs, hexc_i%all_ccoeffs, &
&  hexc_i%Kmesh_dense, hexc_i%Vcp_dense, hexc%crystal%gmet, hexc_i%hinterp, hexc_i%m3_width)

 if( allocated(hexc_i%all_acoeffs) ) then
   ABI_FREE(hexc_i%all_acoeffs)
 end if
 if( allocated(hexc_i%all_bcoeffs) ) then
   ABI_FREE(hexc_i%all_bcoeffs)
 end if
 if( allocated(hexc_i%all_ccoeffs) ) then
   ABI_FREE(hexc_i%all_ccoeffs)
 end if


 if (hexc%BSp%prt_ncham) then
#ifdef HAVE_NETCDF
   MSG_COMMENT("Printing HEXC_I.nc file")
   ncerr = nctk_open_create(ncid, trim(hexc%BS_files%out_basename)//"_HEXC_I.nc", xmpi_comm_self)
   NCF_CHECK_MSG(ncerr, "Creating HEXC_I file")
   call exc_ham_ncwrite(ncid, hexc_i%Kmesh_dense, hexc%BSp, hexc_i%hsize_dense, hexc%BSp%nreh_interp, &
&       hexc%BSp%vcks2t_interp, hexc_i%hinterp, hexc_i%diag_dense)
   NCF_CHECK(nf90_close(ncid))
#else
   ABI_UNUSED(ncid)
#endif
 end if

end subroutine hexc_build_hinterp
!!***

!-------------------------------------------------------------------

!!****f* m_hexc/hexc_compute_subhinterp
!! NAME
!! hexc_compute_subhinterp
!!
!! FUNCTION
!! Compute the interpolation for work_coeffs in dense mesh
!!
!! INPUTS
!! BSp<excparam>=Parameters for BS run
!! grid<double_grid_t>=Double grid info
!! nbnd_coarse=Number of bands (= nbndv * nbndc)
!! interpolator<interpolator_t>=Interpolation info
!! kdense2div=Mapping between dense point and coarse point
!! work_coeffs=Coefficients to be interpolated
!! ikp_dense=Current kpoint
!! overlaps=Wavefunction overlaps
!!
!! OUTPUT
!! Cmat(nbnd_coarse) = Interpolated coefficients
!!
!! PARENTS
!!      m_hexc
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

subroutine hexc_compute_subhinterp(BSp,grid,nbnd_coarse,&
&  interpolator,kdense2div,work_coeffs,Cmat,ikp_dense,overlaps)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbnd_coarse
 integer,intent(in) :: ikp_dense
 type(excparam),intent(in) :: BSp
 type(double_grid_t),intent(in) :: grid
 type(interpolator_t),target,intent(inout) :: interpolator
!arrays
 integer,intent(in) :: kdense2div(grid%nbz_dense)
 complex(gwpc),intent(in) :: overlaps(interpolator%mband_coarse,interpolator%mband_dense,interpolator%nvert)
 complex(dpc),intent(in) :: work_coeffs(nbnd_coarse,interpolator%nvert)
 complex(dpc),intent(out) :: Cmat(nbnd_coarse)

!Local variables ------------------------------
!scalars
 integer,parameter :: spin1=1,spin2=1
 integer :: iv1,ic1
 integer :: icp,ivp,idivp,ibndp_coarse,ibndp_coarse1,ineighbourp
 integer :: indwithnb
 integer :: lumo2,lomo2,humo2,homo2
 complex(dpc) :: tmp_val, tmp2, tmp4
!arrays
 complex(dpc),ABI_CONTIGUOUS pointer :: btemp(:),ctemp(:)

!*********************************************************************

 btemp => interpolator%btemp
 ctemp => interpolator%ctemp

 btemp = czero
 ctemp = czero

 lumo2 = BSp%lumo_spin(spin2)
 lomo2 = BSp%lomo_spin(spin2)
 humo2 = BSp%humo_spin(spin2)
 homo2 = BSp%homo_spin(spin2)

 Cmat = czero

 idivp = kdense2div(ikp_dense)

 do ineighbourp = 1,interpolator%nvert

   btemp(((ineighbourp-1)*nbnd_coarse+1):(ineighbourp*nbnd_coarse)) = &
&       interpolator%interp_factors(ineighbourp,idivp)*work_coeffs(:,ineighbourp)

 end do !ineighbourp

 ! Loop over the (c', v') part of the right transition
 do ivp = lomo2,homo2
   do icp = lumo2,humo2

     ibndp_coarse = (ivp-lomo2)*BSp%maxnbndc+(icp-lumo2+1)
     ! Now we now it_dense, and itp_dense

     do ineighbourp = 1,interpolator%nvert

       do iv1 = lomo2, homo2
         ! BSp%lumo_spin(spin2),BSp%humo_spin(spin2)

         tmp4 = overlaps(iv1,ivp,ineighbourp)

         do ic1 = lumo2,humo2
           tmp2 = GWPC_CONJG(overlaps(ic1,icp,ineighbourp))
           ! BSp%lomo_spin(spin2),BSp%homo_spin(spin2)

           ibndp_coarse1 = (iv1-lomo2)*BSp%maxnbndc+(ic1-lumo2+1)
           indwithnb = (ineighbourp-1)*nbnd_coarse+ibndp_coarse1

           ctemp(indwithnb) = &
&             tmp4 &
&            *tmp2

         end do ! iv1
       end do !ic1

     end do !ineighbourp

     tmp_val = xdotc(interpolator%nvert*nbnd_coarse,ctemp,1,btemp,1)
     !tmp_val = DOT_PRODUCT(ctemp,btemp)

     Cmat(ibndp_coarse) = tmp_val
   end do !ivp
 end do !icp

 nullify(btemp)
 nullify(ctemp)

end subroutine hexc_compute_subhinterp
!!***

!----------------------------------------------------------------------

!!****f* m_hexc/hexc_compute_hinterp
!! NAME
!! hexc_compute_hinterp
!!
!! FUNCTION
!! Compute interpolated matrix elements for methods 2 and 3
!!
!! INPUTS
!! BSp<type(excparam)=The parameter for the Bethe-Salpeter run.
!! hsize_coarse=Size of the coarse Hamiltonian
!! hsize_dense=Size of the dense Hamiltonian
!! hmat(hsize_coarse,hsize_coarse,8)=Excitonic matrix
!! grid<double_grid_t> = Correspondence between coarse and dense k-mesh.
!! nbnd_coarse = Total number of bands
!! interpolator<interpolator_t> = Interpolator
!! kdense2div = Mapping kdense2div
!! acoeffs, bcoeffs, ccoeffs = decomposition "W = a/q^2 + b/q + c"
!! Kmesh_dense<type(kmesh_t)>=The list of k-points in the BZ, IBZ and symmetry tables.
!! Vcp_dense<vcoul_t>=Coulomb interation in G-space on the dense Q-mesh
!! gmet(3,3)=Metric tensor in G-space
!!
!! OUTPUT
!!   hinterp = Interpolated hamiltonian
!!
!! PARENTS
!!      m_hexc
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

subroutine hexc_compute_hinterp(BSp,hsize_coarse,hsize_dense,hmat,grid,nbnd_coarse,&
&  interpolator,kdense2div,acoeffs,bcoeffs,ccoeffs,Kmesh_dense,Vcp_dense,gmet,hinterp,&
&  m3_width)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: hsize_coarse,hsize_dense,nbnd_coarse !,ntrans
 real(dp),intent(in) :: m3_width
 type(excparam),intent(in) :: BSp
 type(double_grid_t),intent(in) :: grid
 type(vcoul_t),intent(in) :: Vcp_dense
 type(kmesh_t),intent(in) :: Kmesh_dense
 type(interpolator_t),target,intent(inout) :: interpolator
!arrays
 integer,intent(in) :: kdense2div(grid%nbz_dense)
 real(dp),intent(in) :: gmet(3,3)
 complex(dpc),intent(in) :: hmat(hsize_coarse,hsize_coarse)
 complex(dpc),intent(in) :: acoeffs(hsize_coarse,hsize_coarse)
 complex(dpc),intent(in) :: bcoeffs(hsize_coarse,hsize_coarse)
 complex(dpc),intent(in) :: ccoeffs(hsize_coarse,hsize_coarse)
 complex(dpc),intent(out) :: hinterp(hsize_dense,hsize_dense)

!Local variables ------------------------------
!scalars
 integer,parameter :: spin1=1, spin2=1
 integer :: ic,iv,iv1,ic1,ik_dense,ik_coarse,it_coarse,it_dense,idiv,ibnd_coarse,ibnd_coarse1,ineighbour
 integer :: icp,ivp,ikp_dense,ikp_coarse,itp_coarse,itp_dense,idivp,ibndp_coarse,ibndp_coarse1,ineighbourp,itp_coarse1
 integer :: itc,it_dense1,indwithnb, corresp_ind
 integer :: limitnbz, inb,ierr
 real(dp) :: factor,vc_sqrt_qbz,qnorm
 complex(dpc) :: term
 logical :: newway, use_herm
!arrays
 real(dp) :: kmkp(3),q2(3),shift(3),qinred(3),tsec(2)
 complex(dpc),allocatable :: Cmat(:,:,:) !Temp matrices for optimized version
 complex(dpc),allocatable :: tmp_Cmat(:)
 complex(dpc),allocatable :: work_coeffs(:,:)
 integer,allocatable :: band2it(:)
 complex(dpc),ABI_CONTIGUOUS pointer :: btemp(:),ctemp(:)

!************************************************************************

 call timab(696,1,tsec)

 newway = .True.
 use_herm = .True.

 if (any(BSp%interp_mode == [2,3])) then
   if(Vcp_dense%mode /= 'CRYSTAL' .and. Vcp_dense%mode /= 'AUXILIARY_FUNCTION') then
     MSG_BUG('Vcp_dense%mode not implemented yet !')
   end if
 end if

 if(BSp%nsppol > 1) then
   MSG_BUG("nsppol > 1 not yet implemented")
 end if

 factor = one/grid%ndiv

 hinterp = czero; term = czero

 ABI_MALLOC_OR_DIE(Cmat,(nbnd_coarse,nbnd_coarse,interpolator%nvert), ierr)
 Cmat = czero

 ABI_MALLOC(band2it,(nbnd_coarse))
 call int_alloc_work(interpolator,nbnd_coarse*interpolator%nvert)

 btemp => interpolator%btemp
 ctemp => interpolator%ctemp

 if(newway) then
   ABI_MALLOC(tmp_Cmat,(nbnd_coarse))
   ABI_MALLOC(work_coeffs,(nbnd_coarse,interpolator%nvert))
 end if

 do ik_dense = 1,grid%nbz_dense
   write(std_out,*) "Kdense = ",ik_dense,"/",grid%nbz_dense
   ik_coarse = grid%dense_to_coarse(ik_dense)
   if(use_herm) then
     limitnbz = ik_dense
   else
     limitnbz = grid%nbz_dense
   end if

   do ikp_dense = 1,limitnbz
     ikp_coarse = grid%dense_to_coarse(ikp_dense)

     do iv1 = BSp%lomo_spin(spin2),BSp%homo_spin(spin2)
       do ic1 = BSp%lumo_spin(spin2),BSp%humo_spin(spin2)
         itp_coarse1 = BSp%vcks2t(iv1,ic1,ikp_coarse,spin2)
         ibndp_coarse1 = (iv1-BSp%lomo_spin(spin2))*BSp%maxnbndc+(ic1-BSp%lumo_spin(spin2)+1)

         band2it(ibndp_coarse1) = itp_coarse1
       end do
     end do


     if (any(BSp%interp_mode == [2,3])) then
       ! Check if we are along the diagonal
       kmkp = Kmesh_dense%bz(:,ik_dense) - Kmesh_dense%bz(:,ikp_dense)

       call wrap2_pmhalf(kmkp(:),q2(:),shift(:))
       qinred = MATMUL(grid%kptrlatt_coarse,q2)

       ! We are outside the diagonal
       if (BSp%interp_mode==3 .and. ANY((ABS(qinred)-tol7) > m3_width)) cycle

       qnorm = two_pi*SQRT(DOT_PRODUCT(q2,MATMUL(gmet,q2)))

       if(ALL(ABS(q2(:)) < 1.e-3)) then
         vc_sqrt_qbz = SQRT(Vcp_dense%i_sz)
       else
         vc_sqrt_qbz = SQRT(four_pi/qnorm**2)
       end if

       !!DEBUG CHK !
       !!COMPUTE Qpoint
       !call findqg0(iq_bz,g0,kmkp,Qmesh_dense%nbz,Qmesh_dense%bz,BSp%mG0)

       !! * Get iq_ibz, and symmetries from iq_bz
       !call get_BZ_item(Qmesh_dense,iq_bz,qbz,iq_ibz,isym_q,itim_q)

       !if(iq_ibz > 1 .and. ABS(vc_sqrt_qbz - Vcp_dense%vc_sqrt(1,iq_ibz)) > 1.e-3) then
       !   write(*,*) "vc_sqrt_qbz = ",vc_sqrt_qbz
       !   write(*,*) "Vcp_dense%vc_sqrt(1,iq_ibz) = ",Vcp_dense%vc_sqrt(1,iq_ibz)
       !   MSG_ERROR("vcp are not the same !")
       !else if(iq_ibz == 1 .and. ABS(vc_sqrt_qbz - SQRT(Vcp_dense%i_sz)) > 1.e-3) then
       !   write(*,*) "vc_sqrt_qbz = ",vc_sqrt_qbz
       !   write(*,*) "SQRT(Vcp_dense%i_sz) = ",SQRT(Vcp_dense%i_sz)
       !   MSG_ERROR("vcp are not the same !")
       !end if
       !!END DEBUG CHK !
     end if

     if(newway) then

       Cmat = czero

       work_coeffs = czero

       do ineighbour = 1,interpolator%nvert

         ! Loop over the (c, v) part of the left transition
         do iv = BSp%lomo_spin(spin1),BSp%homo_spin(spin1)
           do ic = BSp%lumo_spin(spin1),BSp%humo_spin(spin1)

             it_dense = BSp%vcks2t_interp(iv,ic,ik_dense,spin1)
             it_coarse = BSp%vcks2t(iv,ic,ik_coarse,spin1)
             ibnd_coarse = (iv-BSp%lomo_spin(spin1))*BSp%maxnbndc+(ic-BSp%lumo_spin(spin1)+1)

             itc = interpolator%corresp(it_coarse,ineighbour,spin1)

             if (any(BSp%interp_mode == [1,3,4])) then

               !work_coeffs(:,:) = hmat(itc,band2it(:),:)
               do inb = 1,interpolator%nvert
                 work_coeffs(:,inb) = hmat(itc,interpolator%corresp(band2it(:),inb,spin2))
               end do

               call hexc_compute_subhinterp(BSp,grid,nbnd_coarse,&
&        interpolator,kdense2div,&
&        work_coeffs,tmp_Cmat,ikp_dense,&
&        interpolator%overlaps(:,:,:,ikp_dense,spin2))

               if(any(BSp%interp_mode == [1,4])) then
                 Cmat(ibnd_coarse,:,ineighbour) = tmp_Cmat
               else if (BSp%interp_mode == 3) then
                 Cmat(ibnd_coarse,:,ineighbour) = -tmp_Cmat
               end if
             end if


             if (any(BSp%interp_mode == [2,3])) then
               !work_coeffs(:,:) = acoeffs(itc,band2it(:),:)
               do inb = 1,interpolator%nvert
                 work_coeffs(:,inb) = acoeffs(itc,interpolator%corresp(band2it(:),inb,spin2))
               end do

               call hexc_compute_subhinterp(BSp,grid,nbnd_coarse,&
&        interpolator,kdense2div,&
&        work_coeffs,tmp_Cmat,ikp_dense,&
&        interpolator%overlaps(:,:,:,ikp_dense,spin2))

               tmp_Cmat = tmp_Cmat * (vc_sqrt_qbz**2)
               Cmat(ibnd_coarse,:,ineighbour) = Cmat(ibnd_coarse,:,ineighbour) + tmp_Cmat
             end if


             if (any(BSp%interp_mode == [2,3])) then
               !work_coeffs(:,:) = bcoeffs(itc,band2it(:),:)
               do inb = 1,interpolator%nvert
                 work_coeffs(:,inb) = bcoeffs(itc,interpolator%corresp(band2it(:),inb,spin2))
               end do

               call hexc_compute_subhinterp(BSp,grid,nbnd_coarse,&
&        interpolator,kdense2div,&
&        work_coeffs,tmp_Cmat,ikp_dense,&
&        interpolator%overlaps(:,:,:,ikp_dense,spin2))

               tmp_Cmat = tmp_Cmat * (vc_sqrt_qbz)
               Cmat(ibnd_coarse,:,ineighbour) = Cmat(ibnd_coarse,:,ineighbour) + tmp_Cmat
             end if


             if (any(BSp%interp_mode == [2,3])) then
               !work_coeffs(:,:) = ccoeffs(itc,band2it(:),:)
               do inb = 1,interpolator%nvert
                 work_coeffs(:,inb) = ccoeffs(itc,interpolator%corresp(band2it(:),inb,spin2))
               end do

               call hexc_compute_subhinterp(BSp,grid,nbnd_coarse,&
&        interpolator,kdense2div,&
&        work_coeffs,tmp_Cmat,ikp_dense,&
&        interpolator%overlaps(:,:,:,ikp_dense,spin2))

               Cmat(ibnd_coarse,:,ineighbour) = Cmat(ibnd_coarse,:,ineighbour) + tmp_Cmat
             end if
           end do ! ic
         end do ! iv
       end do ! ineighbour

     else
       ! Loop over the (c, v) part of the left transition
       do iv = BSp%lomo_spin(spin1),BSp%homo_spin(spin1)
         do ic = BSp%lumo_spin(spin1),BSp%humo_spin(spin1)

           it_dense = BSp%vcks2t_interp(iv,ic,ik_dense,spin1)
           it_coarse = BSp%vcks2t(iv,ic,ik_coarse,spin1)
           ibnd_coarse = (iv-BSp%lomo_spin(spin1))*BSp%maxnbndc+(ic-BSp%lumo_spin(spin1)+1)

           ! Loop over the (c', v') part of the right transition
           do ivp = BSp%lomo_spin(spin2),BSp%homo_spin(spin2)
             do icp = BSp%lumo_spin(spin2),BSp%humo_spin(spin2)

               itp_dense = BSp%vcks2t_interp(ivp,icp,ikp_dense,spin2)
               itp_coarse = BSp%vcks2t(ivp,icp,ikp_coarse,spin2)
               ibndp_coarse = (ivp-Bsp%lomo_spin(spin2))*BSp%maxnbndc+(icp-BSp%lumo_spin(spin2)+1)
               ! Now we now it_dense, and itp_dense

               idivp = kdense2div(ikp_dense)

               btemp = czero; ctemp = czero

               ! MG TODO: This way of looping is not optimal
               do ineighbour = 1,interpolator%nvert
                 itc = interpolator%corresp(it_coarse,ineighbour,spin1)

                 do ineighbourp = 1,interpolator%nvert

                   do iv1 = BSp%lomo_spin(spin2),BSp%homo_spin(spin2)
                     do ic1 = BSp%lumo_spin(spin2),BSp%humo_spin(spin2)

                       ibndp_coarse1 = (iv1-BSp%lomo_spin(spin2))*BSp%maxnbndc+(ic1-BSp%lumo_spin(spin2)+1)
                       indwithnb = (ineighbourp-1)*nbnd_coarse+ibndp_coarse1
                       itp_coarse1 = BSp%vcks2t(iv1,ic1,ikp_coarse,spin2)
                       corresp_ind = interpolator%corresp(itp_coarse1,ineighbourp,spin2)

                       select case (BSp%interp_mode)
                       case (1,4)
                         interpolator%btemp(indwithnb) = hmat(itc,corresp_ind)
                       case (2)
                         interpolator%btemp(indwithnb) = acoeffs(itc,corresp_ind)*(vc_sqrt_qbz**2) &
&                                         + bcoeffs(itc,corresp_ind)*(vc_sqrt_qbz) &
&                                         + ccoeffs(itc,corresp_ind)
                       case (3)
                         ! Diff between divergence and hmat
                         interpolator%btemp(indwithnb) = acoeffs(itc,corresp_ind)*(vc_sqrt_qbz**2) &
&                                         + bcoeffs(itc,corresp_ind)*(vc_sqrt_qbz) &
&                                         + ccoeffs(itc,corresp_ind) &
&                                         - hmat(itc,corresp_ind)
                       case default
                         MSG_ERROR("Wrong Bsp%interp_mode")
                       end select

                       ctemp(indwithnb) = &
&                        interpolator%overlaps(iv1,ivp,ineighbourp,ikp_dense,spin2) &
&                        * GWPC_CONJG(interpolator%overlaps(ic1,icp,ineighbourp,ikp_dense,spin2)) &
&                        *interpolator%interp_factors(ineighbourp,idivp)
                     end do ! ic1
                   end do !iv1

                 end do !ineighbourp
                 Cmat(ibnd_coarse,ibndp_coarse,ineighbour) = xdotc(interpolator%nvert*nbnd_coarse,&
&                                   ctemp,1,btemp,1)
               end do !ineighbour

             end do !icp
           end do !ivp

         end do !ic
       end do !iv

     end if

     do iv = BSp%lomo_spin(spin1),BSp%homo_spin(spin1)
       do ic = BSp%lumo_spin(spin1),BSp%humo_spin(spin1)
         it_dense = BSp%vcks2t_interp(iv,ic,ik_dense,spin1)
         it_coarse = BSp%vcks2t(iv,ic,ik_coarse,spin1)
         ibnd_coarse = (iv-BSp%lomo_spin(spin1))*BSp%maxnbndc+(ic-BSp%lumo_spin(spin1)+1)

         idiv = kdense2div(ik_dense)

         do ivp = BSp%lomo_spin(spin2),BSp%homo_spin(spin2)
           do icp = BSp%lumo_spin(spin2),BSp%humo_spin(spin2)
             itp_dense = BSp%vcks2t_interp(ivp,icp,ikp_dense,spin2)
             itp_coarse = BSp%vcks2t(ivp,icp,ikp_coarse,spin2)
             ibndp_coarse = (ivp-Bsp%lomo_spin(spin2))*BSp%maxnbndc+(icp-BSp%lumo_spin(spin2)+1)

             ! Hermicity
             if(use_herm .and. it_dense < itp_dense) then
               continue
             end if

             !btemp = czero; ctemp = czero

             do ineighbour = 1,interpolator%nvert
               do iv1 = BSp%lomo_spin(spin1),BSp%homo_spin(spin1)
                 do ic1 = BSp%lumo_spin(spin1),BSp%humo_spin(spin1)
                   ibnd_coarse1 = (iv1-BSp%lomo_spin(spin1))*BSp%maxnbndc+(ic1-BSp%lumo_spin(spin1)+1)
                   it_dense1 = BSp%vcks2t_interp(iv1,ic1,ik_dense,spin1)
                   indwithnb = (ineighbour-1)*nbnd_coarse+ibnd_coarse1

                   btemp(indwithnb) = Cmat(ibnd_coarse1,ibndp_coarse,ineighbour)

                   ctemp(indwithnb) = GWPC_CONJG(interpolator%overlaps(iv1,iv,ineighbour,ik_dense,spin1)) &
&                                    *interpolator%overlaps(ic1,ic,ineighbour,ik_dense,spin1) &
&                                    *interpolator%interp_factors(ineighbour,idiv)
                 end do !ic1
               end do !iv1
             end do !ineighbour

             ! Save interpolated value.
             hinterp(it_dense,itp_dense) = xdotc(interpolator%nvert*nbnd_coarse,ctemp,1,btemp,1)
             !DOT_PRODUCT(ctemp,btemp)

           end do !icp
         end do !ivp

       end do !ic
     end do !iv

   end do !ikp
 end do !ik

 ! Enforce hermiticity
 if(use_herm) then
   do itp_dense = 1,BSp%nreh_interp(spin2)
     do it_dense = itp_dense,BSp%nreh_interp(spin1)
       if(it_dense == itp_dense) then
         hinterp(itp_dense,it_dense) = DBLE(hinterp(it_dense,itp_dense))
       else
         hinterp(itp_dense,it_dense) = CONJG(hinterp(it_dense,itp_dense))
       end if
     end do
   end do
 end if

 ABI_FREE(Cmat)
 ABI_FREE(band2it)

 if(newway) then
   ABI_FREE(tmp_Cmat)
   ABI_FREE(work_coeffs)
 end if

 nullify(ctemp)
 nullify(btemp)
 call int_free_work(interpolator)

 hinterp = hinterp*factor

 call timab(696,2,tsec)

end subroutine hexc_compute_hinterp
!!***

!----------------------------------------------------------------------

!!****f* m_hexc/hexc_free
!! NAME
!! hexc_free
!!
!! FUNCTION
!! Destroy the interpolator object in memory
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

subroutine hexc_free(hexc)

!Arguments ---------------------------
 type(hexc_t),intent(inout) :: hexc

!*****************************************************************************

 if( associated(hexc%bsp) ) then
   nullify(hexc%bsp)
 end if

 if( associated(hexc%crystal) ) then
   nullify(hexc%crystal)
 end if

 if( associated(hexc%kmesh_coarse) ) then
   nullify(hexc%kmesh_coarse)
 end if

 if( allocated(hexc%hreso) ) then
   ABI_FREE(hexc%hreso)
 end if

 if( allocated(hexc%hcoup) ) then
   ABI_FREE(hexc%hcoup)
 end if

 if( allocated(hexc%diag_coarse) ) then
   ABI_FREE(hexc%diag_coarse)
 end if

end subroutine hexc_free
!!***

!-------------------------------------------------------------------

!!****f* m_hexc/hexc_interp_free
!! NAME
!! hexc_interp_free
!!
!! FUNCTION
!!  Destroy the interpolator object in memory
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

subroutine hexc_interp_free(hexc_i)

!Arguments ---------------------------
 type(hexc_interp_t),intent(inout) :: hexc_i

!*****************************************************************************

 if( allocated(hexc_i%kdense2div) ) then
   ABI_FREE(hexc_i%kdense2div)
 end if

 if( allocated(hexc_i%div2kdense) ) then
   ABI_FREE(hexc_i%div2kdense)
 end if

 if( allocated(hexc_i%diag_dense) ) then
   ABI_FREE(hexc_i%diag_dense)
 end if

 if( allocated(hexc_i%hinterp) ) then
   ABI_FREE(hexc_i%hinterp)
 end if

 if( allocated(hexc_i%all_hmat) ) then
   ABI_FREE(hexc_i%all_hmat)
 end if

 if( allocated(hexc_i%all_acoeffs) ) then
   ABI_FREE(hexc_i%all_acoeffs)
 end if

 if( allocated(hexc_i%all_bcoeffs) ) then
   ABI_FREE(hexc_i%all_bcoeffs)
 end if

 if( allocated(hexc_i%all_ccoeffs) ) then
   ABI_FREE(hexc_i%all_ccoeffs)
 end if

 if( associated(hexc_i%kmesh_dense) ) then
   nullify(hexc_i%kmesh_dense)
 end if

 if( associated(hexc_i%vcp_dense) ) then
   nullify(hexc_i%vcp_dense)
 end if

 call interpolator_free(hexc_i%interpolator)

end subroutine hexc_interp_free
!!***

!----------------------------------------------------------------------

!!****f* m_hexc/hexc_interp_matmul
!! NAME
!! hexc_interp_matmul
!!
!! FUNCTION
!! Compute matrix-vector product Hmat * phi by interpolating coarse Hmat
!!
!! INPUTS
!!  BSp<type(excparam)>=Parameters defining the BS calculation
!!  hsize_coarse = Size of the coarse hamiltonian
!!  hsize_dense = Size of the dense hamiltonian
!!  hmat(hsize_coarse,hsize_coarse,8) = coarse hamiltonian
!!  phi(hsize_dense) = ket on which apply the matrix
!!  grid <double_grid_t> = Correspondence between coarse and dense k-mesh.
!!  nbnd_coarse = Total number of bands
!!  interpolator<interpolator_t> = Interpolator
!!  div2kdense = Mapping from coarse and division -> dense mesh
!!  kdense2div = Mapping from dense mesh -> division
!!
!! OUTPUT
!!  hphi(hsize_dense) = Interp(hmat)*phi
!!
!! PARENTS
!!      m_hexc
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

subroutine hexc_interp_matmul(BSp,hsize_coarse,hsize_dense,hmat,phi,hphi,grid,&
&   nbnd_coarse,interpolator,div2kdense,kdense2div)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: hsize_coarse,hsize_dense,nbnd_coarse !,ntrans
 type(excparam),intent(in) :: BSp
 type(double_grid_t),intent(in) :: grid
 type(interpolator_t),intent(in) :: interpolator
!arrays
 integer,intent(in) :: div2kdense(grid%nbz_coarse,grid%ndiv), kdense2div(grid%nbz_dense)
 complex(dpc),intent(in) :: phi(hsize_dense)
 complex(dpc),intent(in) :: hmat(hsize_coarse,hsize_coarse)
 complex(dpc),intent(inout) :: hphi(hsize_dense)

!Local variables ------------------------------
!scalars
 integer :: itt,ik_dense,ik_coarse,it_coarse
 integer :: ic,iv,iv1,ic1, ibnd_coarse
 integer :: ibnd_coarse1
 integer :: ineighbour,idense,ikpt
 integer :: my_k1,my_k2,ind_with_nb,is, is1
 real(dp) :: factor
 complex(dpc) :: tmp
 logical,parameter :: use_blas=.True.
!arrays
 integer :: allindices(nbnd_coarse)
 complex(dpc) :: allp(hsize_coarse,interpolator%nvert), test(hsize_coarse)
 complex(dpc) :: ophi(grid%nbz_dense,interpolator%nvert,nbnd_coarse)
 complex(dpc),allocatable :: b(:), c(:),A(:,:)
 complex(dpc),allocatable :: tmp_array(:), tmp_array2(:,:)

!************************************************************************

 factor = one/grid%ndiv

 !hphi = czero

 ! Outer index : k point in the dense zone
 ! Sum over vc
 ! Index of result : k point in the dense zone, v2,c2,neighbour

 ! Parallelization on nbz in the coarse mesh !
 my_k1 = 1
 my_k2 = grid%nbz_coarse

 ABI_MALLOC(A,(interpolator%nvert*nbnd_coarse,nbnd_coarse))
 ABI_MALLOC(b,(nbnd_coarse))
 ABI_MALLOC(c,(interpolator%nvert*nbnd_coarse))

 c = czero; ophi = czero

 do ik_dense = 1,grid%nbz_dense
   ! if( ik_dense is not in my set of k-points)
   !   ! continue
   !
   do is1 = 1, BSp%nsppol
     do iv1 = BSp%lomo_spin(is1),Bsp%homo_spin(is1)
       do ic1 = BSp%lumo_spin(is1),Bsp%humo_spin(is1)
         ibnd_coarse = (iv1-BSp%lomo_spin(is1))*BSp%maxnbndc+(ic1-BSp%lumo_spin(is1)+1)
         itt = BSp%vcks2t_interp(iv1,ic1,ik_dense,is1)
         allindices(ibnd_coarse) = itt
       end do !ic1
     end do !iv1
   end do !is1

   b(:) = phi(allindices(:))

   do is = 1, BSp%nsppol
     do iv = BSp%lomo_spin(is),Bsp%homo_spin(is)
       do ic = BSp%lumo_spin(is),Bsp%humo_spin(is)
         ibnd_coarse = (iv-BSp%lomo_spin(is))*BSp%maxnbndc+(ic-BSp%lumo_spin(is)+1)
         idense = Bsp%vcks2t_interp(iv,ic,ik_dense,is)

         do ineighbour = 1,interpolator%nvert
           ind_with_nb = (ineighbour-1)*(nbnd_coarse)+ibnd_coarse

           !A(ind_with_nb,:) = overlaps(allindices(:),ibnd_coarse,ineighbour)

           ! Should be optimized !!!
           do iv1 = BSp%lomo_spin(is),Bsp%homo_spin(is)
             do ic1 = BSp%lumo_spin(is),Bsp%humo_spin(is)
               ibnd_coarse1 = (iv1-BSp%lomo_spin(is))*BSp%maxnbndc+(ic1-BSp%lumo_spin(is)+1)
               A(ind_with_nb,ibnd_coarse1) = GWPC_CONJG(interpolator%overlaps(iv,iv1,ineighbour,ik_dense,is)) &
&                                          *interpolator%overlaps(ic,ic1,ineighbour,ik_dense,is)
             end do !ic1
           end do !iv1
         end do !ineighbour
       end do !ic
     end do !iv
   end do !is

   if(use_blas) then
     call xgemv('N',interpolator%nvert*nbnd_coarse,nbnd_coarse,cone,A,interpolator%nvert*nbnd_coarse,b,1,czero,c,1)
   else
     c = MATMUL(A,b)
   end if

   do is = 1, BSp%nsppol
     do iv = BSp%lomo_spin(is),BSp%homo_spin(is)
       do ic = BSp%lumo_spin(is),BSp%humo_spin(is)
         ibnd_coarse = (iv-BSp%lomo_spin(is))*BSp%maxnbndc+(ic-BSp%lumo_spin(is)+1)
         do ineighbour = 1,interpolator%nvert
           ind_with_nb = (ineighbour-1)*(nbnd_coarse)+ibnd_coarse
           ophi(ik_dense,ineighbour,ibnd_coarse) = c(ind_with_nb)
         end do !ineighbour
       end do !ic
     end do !iv
   end do !is

 end do !ik_dense

 ABI_FREE(A)
 ABI_FREE(b)
 ABI_FREE(c)

 !call xmpi_sum_(ophi,comm,ierr)

 ! Outer index : k,v,c in the coarse zone, ineighbour
 ! Sum over all k-dense relative to one coarse point
 ! Index of result : k,v,c in the coarse zone, ineighbour

 ABI_MALLOC(b,(grid%ndiv))
 ABI_MALLOC(c,(grid%ndiv))

 allp = czero

 do is = 1, BSp%nsppol
   do ineighbour = 1,interpolator%nvert

     do it_coarse = 1, BSp%nreh(is)
       ibnd_coarse = (Bsp%trans(it_coarse,is)%v-BSp%lomo_spin(is))*BSp%maxnbndc+&
&            (BSp%Trans(it_coarse,is)%c-BSp%lumo_spin(is)+1)
       ik_coarse = BSp%trans(it_coarse,is)%k
       !b(:) = interp_factors(it_coarse,ineighbour,:)
       b(:) = interpolator%interp_factors(ineighbour,:)
       !c(:) = ophi(indices(it_coarse,:),ineighbour,ibnd_coarse)
       c(:) = ophi(div2kdense(ik_coarse,:),ineighbour,ibnd_coarse)
       tmp = DOT_PRODUCT(b,c)
       allp(it_coarse,ineighbour) = tmp
     end do

   end do
 end do

 !call xmpi_sum_(allp,comm,ierr)

 ABI_FREE(b)
 ABI_FREE(c)

 ABI_MALLOC(tmp_array,(hsize_coarse))
 ABI_MALLOC(tmp_array2,(hsize_coarse,hsize_coarse))
 tmp_array(:) = czero
 tmp_array2(:,:) = czero

 test = czero

 ! Second step : Multiplication by hmat
 do ineighbour = 1,interpolator%nvert
   if(use_blas) then
     !call xgemv('N',hsize_coarse,hsize_coarse,cone,factor*(hmat(:,:,ineighbour)),hsize_coarse,allp(:,ineighbour),1,czero,tmp_array,1)
     !tmp_array2 = hmat(:,:,ineighbour)
     tmp_array2 = hmat(:,interpolator%corresp(:,ineighbour,1)) ! 1 is for spin
     tmp_array2 = factor*tmp_array2
     call xgemv('N',hsize_coarse,hsize_coarse,cone,tmp_array2,hsize_coarse,allp(:,ineighbour),1,czero,tmp_array,1)
     test = test + tmp_array
   else
     test = test+MATMUL(factor*(hmat(:,interpolator%corresp(:,ineighbour,1))),allp(:,ineighbour))
   end if
 end do

 ABI_FREE(tmp_array)
 ABI_FREE(tmp_array2)

 ! Outer index : ineighbour
 ! Sum over all v c
 ! Index of result : ineighbour, k_dense, v,c
 ABI_MALLOC(A,(nbnd_coarse,nbnd_coarse))
 ABI_MALLOC(b,(nbnd_coarse))
 ABI_MALLOC(c,(nbnd_coarse))
 c = czero

 do ineighbour = 1,interpolator%nvert
   do ik_dense = 1,grid%nbz_dense

     do is1 = 1, Bsp%nsppol
       do iv1 = Bsp%lomo_spin(is1),Bsp%homo_spin(is1)
         do ic1 = BSp%lumo_spin(is1), Bsp%humo_spin(is1)
           ibnd_coarse = (iv1-BSp%lomo_spin(is1))*BSp%maxnbndc+(ic1-BSp%lumo_spin(is1)+1)

           ik_coarse = grid%dense_to_coarse(ik_dense)
           itt = BSp%vcks2t(iv1,ic1,ik_coarse,is1)
           b(ibnd_coarse) = test(interpolator%corresp(itt,ineighbour,is1))
         end do ! ic1
       end do ! iv1
     end do ! is1

     do is = 1, BSp%nsppol
       do iv = BSp%lomo_spin(is),Bsp%homo_spin(is)
         do ic = BSp%lumo_spin(is),BSp%humo_spin(is)
           ibnd_coarse = (iv-BSp%lomo_spin(is))*Bsp%maxnbndc+(ic-BSp%lumo_spin(is)+1)
           idense = BSp%vcks2t_interp(iv,ic,ik_dense,is)

           !A(ibnd_coarse,:) = CONJG(overlaps(idense,:,ineighbour))

           ! Should be optimized !!!
           do iv1 = BSp%lomo_spin(is),Bsp%homo_spin(is)
             do ic1 = BSp%lumo_spin(is),Bsp%humo_spin(is)
               ibnd_coarse1 = (iv1-BSp%lomo_spin(is))*BSp%maxnbndc+(ic1-BSp%lumo_spin(is)+1)
               A(ibnd_coarse,ibnd_coarse1) = (interpolator%overlaps(iv1,iv,ineighbour,ik_dense,is)) &
&                                          *GWPC_CONJG(interpolator%overlaps(ic1,ic,ineighbour,ik_dense,is))
             end do !ic1
           end do !iv1
         end do ! ic
       end do !iv
     end do !is

     if(use_blas) then
       call xgemv('N',nbnd_coarse,nbnd_coarse,cone,A,nbnd_coarse,b,1,czero,c,1)
     else
       c = MATMUL(A,b)
     end if

     do is = 1, BSp%nsppol
       do iv = BSp%lomo_spin(is),Bsp%homo_spin(is)
         do ic = BSp%lumo_spin(is),BSp%humo_spin(is)
           ibnd_coarse = (iv-BSp%lomo_spin(is))*BSp%maxnbndc+(ic-BSp%lumo_spin(is)+1)
           idense = Bsp%vcks2t_interp(iv,ic,ik_dense,is)
           !ophi(ik_dense,ineighbour,ibnd_coarse) = c(idense)
           ophi(ik_dense,ineighbour,ibnd_coarse) = c(ibnd_coarse)
         end do
       end do
     end do

   end do ! ik_dense
 end do ! ineighbour

 !call xmpi_sum_(ophi,comm,ierr)

 ABI_FREE(A)
 ABI_FREE(b)
 ABI_FREE(c)

 ! Outer indices : it_dense
 ! Sum over neighbours
 ! Index of result : it_dense (ik,ic,iv dense)

 ABI_MALLOC(b,(interpolator%nvert))
 ABI_MALLOC(c,(interpolator%nvert))

 do is = 1, BSp%nsppol
   do itt = 1,BSp%nreh_interp(is)
    ! From itt -> ik_ibz,ic,iv
    ik_dense = BSp%Trans_interp(itt,is)%k
    ic = BSp%Trans_interp(itt,is)%c
    iv = BSp%Trans_interp(itt,is)%v

    ! From ik_ibz in the dense mesh -> indices_dense
    ik_coarse = grid%dense_to_coarse(ik_dense)
    it_coarse = BSp%vcks2t(iv,ic,ik_coarse,is)

    ibnd_coarse = (iv-BSp%lomo_spin(is))*BSp%maxnbndc+(ic-BSp%lumo_spin(is)+1)

    ikpt = kdense2div(ik_dense)
    !ikpt = -1
    !do ix = 1,grid%ndiv
    !  if (indices(it_coarse,ix) == ik_dense) then
    !    ikpt = ix
    !    exit
    !  end if
    !end do
    !ABI_CHECK(ikpt/=-1,"Cannot find ik_dense")

    !b = interp_factors(it_coarse,:,ikpt)
    b = interpolator%interp_factors(:,ikpt)
    c =  ophi(ik_dense,:,ibnd_coarse)

    hphi(itt) = hphi(itt) + xdotc(interpolator%nvert, b, 1, c, 1)
   end do
 end do

 ABI_FREE(b)
 ABI_FREE(c)

end subroutine hexc_interp_matmul
!!***

!-------------------------------------------------------------------

!!****f* m_hexc/hexc_matmul_tda
!! NAME
!! hexc_matmul_tda
!!
!! FUNCTION
!! Compute H | \psi >
!!
!! INPUTS
!! hexc<hexc_t> = Excitonic hamiltonian
!! hexc_i<hexc_interp_t> = Interpolated excitonic hamiltonian
!! phi = Input ket
!!
!! OUTPUT
!! hphi = hreso * phi
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

subroutine hexc_matmul_tda(hexc, hexc_i, phi, hphi)

!Arguments ---------------------------
 type(hexc_t),intent(in) :: hexc
 type(hexc_interp_t),intent(in) :: hexc_i
 complex(dpc),intent(in) :: phi(hexc%hsize)
 complex(dpc),intent(out) :: hphi(hexc%hsize)

!Local variables ---------------------
 integer :: ierr
 real(dp) :: tsec(2)

!*****************************************************************************

 call timab(697,1,tsec)

 if(hexc%BSp%use_interp) then
   hphi = hexc_i%diag_dense * phi

   if (any(hexc%BSp%interp_mode == [2,3,4])) then
     ! hphi = hphi + MATMUL(hinterp,phi)
     call xgemv('N',hexc_i%hsize_dense,hexc_i%hsize_dense,cone,hexc_i%hinterp,hexc_i%hsize_dense,phi,1,cone,hphi,1)
     if (any(hexc%BSp%interp_mode == [2,4])) then
       call timab(697,2,tsec)
       return ! We are done
     end if
   end if

   call hexc_interp_matmul(hexc%bsp, hexc%hsize_coarse, hexc_i%hsize_dense, hexc_i%all_hmat, phi, hphi, &
&     hexc_i%interpolator%double_grid,hexc%nbnd_coarse, hexc_i%interpolator, hexc_i%div2kdense, hexc_i%kdense2div)

 else ! No interpolation
   call xgemv('N',hexc%hsize,hexc%my_nt,cone,hexc%hreso,hexc%hsize,phi,1,czero,hphi,1)
   call xmpi_sum(hphi,hexc%comm,ierr)
 end if

 call timab(697,2,tsec)

end subroutine hexc_matmul_tda
!!***

!-------------------------------------------------------------------

!!****f* m_hexc/hexc_matmul_elphon
!! NAME
!! hexc_matmul_elphon
!!
!! FUNCTION
!! Compute H | \psi > + E_{elphon} | \psi >
!!
!! INPUTS
!! hexc<hexc_t> = Excitonic hamiltonian
!! phi = Input ket
!! ep_renorm = vector with electron-phonon renorms
!! op = 'N' for H | psi >, 'C' for H^\dagger | psi >
!!
!! OUTPUT
!! hphi = hreso * phi + ep_renorm * phi
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

subroutine hexc_matmul_elphon(hexc, phi, hphi, op, ep_renorm)

!Arguments ---------------------------
 type(hexc_t),intent(in) :: hexc
 character,intent(in) :: op
 complex(dpc),intent(in) :: phi(hexc%my_nt)
 complex(dpc),intent(out) :: hphi(hexc%hsize)
 complex(dpc),intent(in) :: ep_renorm(hexc%hsize)

!Local variables ---------------------
 integer :: ierr
 real(dp) :: tsec(2)

!*****************************************************************************

 call timab(697,1,tsec)

 if(hexc%BSp%use_interp) then
   MSG_ERROR('Not yet implemented with interpolation !')
 else ! No interpolation
   ! As our matrix is hermitian (hreso), we should always use 'N' here (it is stored column-wise !)
   call xgemv('N',hexc%hsize,hexc%my_nt,cone,hexc%hreso,hexc%hsize,phi,1,czero,hphi,1)

   !!! ep_renorm is stored on each cpu
   if (op == 'N') then
     hphi(hexc%my_t1:hexc%my_t2) = hphi(hexc%my_t1:hexc%my_t2) + ep_renorm(hexc%my_t1:hexc%my_t2) * phi
   else if(op == 'C') then
     hphi(hexc%my_t1:hexc%my_t2) = hphi(hexc%my_t1:hexc%my_t2) + CONJG(ep_renorm(hexc%my_t1:hexc%my_t2)) * phi
   end if
   call xmpi_sum(hphi,hexc%comm,ierr)
 end if

 call timab(697,2,tsec)

end subroutine hexc_matmul_elphon
!!***

!-------------------------------------------------------------------

!!****f* m_hexc/hexc_matmul_full
!! NAME
!! hexc_matmul_full
!!
!! FUNCTION
!! Compute H | \psi >
!!
!! INPUTS
!! hexc<hexc_t> = Excitonic hamiltonian
!! hexc_i<hexc_interp_t> = Interpolated hamiltonian
!! phi = Input ket
!! parity = -1 or +1 parameter
!!
!! OUTPUT
!! hphi = hreso * phi + parity * hcoup * CONJ(phi)
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

subroutine hexc_matmul_full(hexc, hexc_i, phi, hphi, parity)

!Arguments ---------------------------
 integer,intent(in) :: parity
 type(hexc_t),intent(in) :: hexc
 type(hexc_interp_t),intent(in) :: hexc_i
 complex(dpc),intent(in) :: phi(hexc%hsize)
 complex(dpc),intent(out) :: hphi(hexc%hsize)

!Local variables ---------------------
 real(dp) :: tsec(2)

!*****************************************************************************

 call timab(697,1,tsec)

 ABI_UNUSED(hexc_i%hsize_dense)

 if(hexc%BSp%use_interp) then
   MSG_ERROR("Coupling is not yet implemented with interpolation")
 else ! No interpolation
   hphi = MATMUL(hexc%hreso,phi) + parity * MATMUL(hexc%hcoup,CONJG(phi))
 end if

 call timab(697,2,tsec)

end subroutine hexc_matmul_full
!!***

!-------------------------------------------------------------------

end module m_hexc
!!***
