!!****m* ABINIT/m_mlwfovlp
!! NAME
!!  m_mlwfovlp
!!
!! FUNCTION
!!  Interface with Wannier90
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2020 ABINIT group (BAmadon, CEspejo, FJollet, TRangel, DRH)
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

module m_mlwfovlp

 use defs_basis
 use defs_wannier90
 use m_abicore
 use m_errors
 use m_atomdata
 use m_xmpi
 use m_sort
#ifdef FC_NAG
 use f90_unix_dir
#endif
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk
 use m_hdr
 use m_dtset
 use m_dtfil

 use defs_datatypes, only : pseudopotential_type, ebands_t
 use defs_abitypes, only : MPI_type
 use m_io_tools, only : delete_file, get_unit, open_file
 use m_hide_lapack,     only : matrginv
 use m_fstrings,      only : strcat, sjoin
 use m_numeric_tools, only : uniformrandom, simpson_int, c2r, l2int
 use m_special_funcs,   only : besjm
 use m_geometry,  only : xred2xcart, rotmat, wigner_seitz
 use m_fftcore,  only : sphereboundary
 use m_crystal,  only : crystal_t
 use m_ebands,   only : ebands_ncwrite
 use m_pawang,   only : pawang_type
 use m_pawrad,   only : pawrad_type, simp_gen
 use m_pawtab,   only : pawtab_type
 use m_pawcprj,  only : pawcprj_type
 use m_paw_sphharm, only : ylm_cmplx, initylmr
 use m_paw_overlap, only : smatrix_pawinit
 use m_evdw_wannier, only : evdw_wannier
 use m_fft,            only : fourwf

 implicit none

 private
!!***

 public :: mlwfovlp
!!***

contains
!!***

!!****f* m_mlwfovlp/mlwfovlp
!! NAME
!! mlwfovlp
!!
!! FUNCTION
!! Routine which computes overlap M_{mn}(k,b) and projection A_{mn}(k)
!! for Wannier code (www.wannier.org f90 version).
!! Various file are written (wannier90.*) which can be used to run a
!! separate wannier calculation with the wannier90 code.
!!
!! INPUTS
!!  crystal<crystal_t>=Info on the crystalline structure.
!!  ebands<ebands_t>=The object describing the band structure.
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cg(2,mcg)=planewave coefficients of wavefunctions.
!!  cprj(natom,mcprj)= <p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  dtfil <type(datafiles_type)>=variables related to files
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mgfft=maximum size of 1D FFTs
!!  mgfftc=maximum size of 1D FFTs (coarse grid)
!!  mkmem =number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nfft=(effective) number of FFT grid points (for this processor) (see NOTES at beginning of scfcv)
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  occ(mband*nkpt*nsppol) Occupation number for each band (often 2) for each k point.
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  ucvol=unit cell volume (bohr**3)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      m_outscfcv
!!
!! CHILDREN
!!      initylmr,matrginv,rotmat
!!
!! SOURCE

 subroutine mlwfovlp(crystal, ebands, hdr, atindx1,cg,cprj,dtset,dtfil,eigen,gprimd,kg,&
& mband,mcg,mcprj,mgfftc,mkmem,mpi_enreg,mpw,natom,&
& nattyp,nfft,ngfft,nkpt,npwarr,nsppol,ntypat,occ,&
& pawang,pawrad,pawtab,prtvol,psps,rprimd,ucvol,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mcprj,mgfftc,mkmem,mpw,natom,nfft,nkpt
 integer,intent(in) :: nsppol,ntypat,prtvol
 real(dp),intent(in) :: ucvol
 type(crystal_t),intent(in) :: crystal
 type(ebands_t),intent(in) :: ebands
 type(hdr_type),intent(in) :: hdr
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx1(natom)
 integer :: kg(3,mpw*mkmem),nattyp(ntypat),ngfft(18),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),gprimd(3,3),rprimd(3,3)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: xred(3,natom)
 type(pawcprj_type) :: cprj(natom,mcprj)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: band_index,cplex,i,iatom,iband,iband1,iband2,icgtemp
 integer :: ig,ii,ikg,ierr
 integer :: ikpt,ikpt1,ikpt2,ilmn,intot,isppol,itypat
 integer :: iun(nsppol),iun_plot,iwan,jband,jband1,jband2,jj,jj1,jj2,jj3
 integer :: lmn_size,lproj,lwanniersetup,mwan,mgfft,n1
#if defined HAVE_WANNIER90
 integer :: kk
#ifdef HAVE_NETCDF
 integer :: ncid, ncerr, nrpts
 character(len=fnlen) :: abiwan_fname
 integer :: have_disentangled_spin(nsppol)
 integer,allocatable :: irvec(:,:),ndegen(:)
#endif
#endif
 integer :: n1tmp,n2,n2tmp,n3,n3tmp,n4,n5,n6,nband_k
 integer :: nntot,npw_k,num_nnmax,spacing
 integer :: tim_fourwf
 integer :: master,max_num_bands,nprocs,spaceComm,rank
 integer  :: nwan(nsppol),nband_inc(nsppol),num_bands(nsppol)
 real(dp) :: weight
#if defined HAVE_WANNIER90
 real(dp) :: corrvdw
 complex(dpc) :: caux,caux2,caux3
#endif
 logical :: gamma_only,leig,lmmn,lwannierrun,spinors !,have_disentangled
 character(len=fnlen) :: wfnname
 character(len=1000) :: message
 character(len=fnlen) :: seed_name(nsppol)
 character(len=fnlen) :: fname,filew90_win(nsppol),filew90_wout(nsppol),filew90_amn(nsppol),filew90_ramn(nsppol)
 character(len=fnlen) :: filew90_mmn(nsppol),filew90_eig(nsppol)
!arrays
 integer :: g1temp(3),ngkpt(3)
 integer,allocatable :: g1(:,:,:),gbound(:,:),icg(:,:)
 integer,allocatable:: iwav(:,:,:),kg_k(:,:),ovikp(:,:)
 integer,allocatable :: proj_l(:,:),proj_m(:,:),proj_radial(:,:)
 integer,allocatable :: proj_s_loc(:)
 real(dp) :: real_lattice(3,3)
 real(dp) :: recip_lattice(3,3)
 real(dp),allocatable :: cm1(:,:,:,:,:,:),cm2_paw(:,:,:),cwavef(:,:)
 real(dp),allocatable :: denpot(:,:,:)
 real(dp),allocatable :: eigenvalues_w(:,:,:),fofgout(:,:),fofr(:,:,:,:)
 real(dp),allocatable :: proj_site(:,:,:),proj_x(:,:,:),proj_z(:,:,:),proj_zona(:,:)
 real(dp),allocatable :: wann_centres(:,:,:),wann_spreads(:,:),xcart(:,:)
 real(dp),allocatable :: proj_s_qaxis_loc(:,:)
 complex(dpc),allocatable :: A_paw(:,:,:,:)
 complex(dpc),allocatable :: M_matrix(:,:,:,:,:),U_matrix(:,:,:,:)
 complex(dpc),allocatable :: U_matrix_opt(:,:,:,:)
 complex(dpc),pointer :: A_matrix(:,:,:,:)
 logical,allocatable :: band_in(:,:),lwindow(:,:,:)
 character(len=3),allocatable :: atom_symbols(:)
 logical,allocatable::just_augmentation(:,:)
#if defined HAVE_WANNIER90
 real(dp) :: spreadw(3,nsppol)
 real(dp),allocatable :: csix(:,:,:,:)
 real(dpc),allocatable :: occ_arr(:,:,:),occ_wan(:,:,:)
 real(dp),allocatable :: tdocc_wan(:,:)
#endif

!************************************************************************

 ABI_UNUSED((/crystal%natom, ebands%nkpt, hdr%nkpt/))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!1) Initialize variables and allocations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Some initialization and checks
!
 lwanniersetup=1 ! 1 is mandatory ( 0 is for debug)
!to use lwanniersetup=0, one would need
!to define which bands to exclude.
 lwannierrun=.true.   ! .false. and .true. are possible
 lmmn=.true.          ! .false. and .true. are possible
 leig=.true.          ! .false. and .true. are possible
!
 gamma_only=.false. !not yet implemented
 spinors=.false. !not yet implemented
!
!mpi initialization
!
 spaceComm=MPI_enreg%comm_cell
 nprocs=xmpi_comm_size(spaceComm)
 rank=MPI_enreg%me_kpt
 master=0

!write(std_out,'("master ",i3," rank ",i3," nprocs ",i3)') master,rank,nprocs
!
!Generate seed names for wannier90 files, and file names
!
 call mlwfovlp_seedname(dtfil%fnameabo_w90,filew90_win,filew90_wout,filew90_amn,&
& filew90_ramn,filew90_mmn,filew90_eig,nsppol,seed_name)
!
!Check the validity of input variables
!FIXME: this is not a check, and prints a warning even if the input is fine!
!must be changed to not print anything if kptopt 3 and istwfk 1 (the latter is easier to check)
!
 if(rank==master) then
   write(message, '(a,a,a,a)' ) ch10,&
&   '   mlwfovlp:  you should give k-point in the full brillouin zone ',ch10,&
&   '   with explicit k-points (or kptopt=3) and istwfk 1'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if
!
 if(MPI_enreg%paral_spinor==1) then
   message = ' Parallelization over spinorial components not yet available !'
   ABI_ERROR(message)
 end if

 if(nsppol==2) then
   write(message, '(3a)' ) ch10,&
&   '   mlwfovlp:  Calculating matrices for both spin polarization  ',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if
!
!get lattice parameters in wannier90 format
!
 do i=1, 3
   real_lattice(:,i)=Bohr_Ang*rprimd(i,:)
   recip_lattice(:,i)=two_pi*gprimd(i,:)/Bohr_Ang
 end do
!
 if(psps%npsp/=psps%ntypat) then
   ABI_ERROR("prb npsp")
 end if
!
!Allocations.
!
 num_nnmax=12 !limit fixed for compact structure in wannier_setup.
 ABI_ALLOCATE(g1,(3,nkpt,num_nnmax))
 ABI_ALLOCATE(ovikp,(nkpt,num_nnmax))
 ABI_ALLOCATE(atom_symbols,(natom))
 ABI_ALLOCATE(xcart,(3,natom))
 ABI_ALLOCATE(band_in,(mband,nsppol))
 ABI_ALLOCATE(proj_site,(3,mband,nsppol))
 ABI_ALLOCATE(proj_l,(mband,nsppol))
 ABI_ALLOCATE(proj_m,(mband,nsppol))
 ABI_ALLOCATE(proj_radial,(mband,nsppol))
 ABI_ALLOCATE(proj_x,(3,mband,nsppol))
 ABI_ALLOCATE(proj_s_loc,(mband))
 ABI_ALLOCATE(proj_s_qaxis_loc,(3,mband))
 ABI_ALLOCATE(proj_z,(3,mband,nsppol))
 ABI_ALLOCATE(proj_zona,(mband,nsppol))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!2) Call to  Wannier setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 nullify(A_matrix)

 !
 call mlwfovlp_setup(atom_symbols,band_in,dtset,filew90_win,gamma_only,&
&  g1,lwanniersetup,mband,natom,nband_inc,nkpt,&
&  nntot,num_bands,num_nnmax,nsppol,nwan,ovikp,&
&  proj_l,proj_m,proj_radial,proj_site,proj_s_loc, proj_s_qaxis_loc, proj_x,proj_z,proj_zona,&
&  real_lattice,recip_lattice,rprimd,seed_name,spinors,xcart,xred)

 do isppol=1, nsppol
   write(message, '(6a)' ) ch10,&
&   '   mlwfovlp :  mlwfovlp_setup done -',ch10,&
&   '-  see ',trim(filew90_wout(isppol)),' for details.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end do

!
!some allocations after wannier90 setup
!
 max_num_bands=maxval(num_bands(:))
 mwan=maxval(nwan(:))
 ABI_ALLOCATE(eigenvalues_w,(max_num_bands,nkpt,nsppol))
 ABI_ALLOCATE(M_matrix,(max_num_bands,max_num_bands,nntot,nkpt,nsppol))
 ABI_ALLOCATE(A_matrix,(max_num_bands,mwan,nkpt,nsppol))
 ABI_ALLOCATE(iwav,(mband,nkpt,nsppol))

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!3) Write Eigenvalues (file seed_name.eig)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 if(leig) then
!
!  Assign file unit numbers
   if(rank==master) then
     iun(1)=444
     if(nsppol==2) iun(2)=455
     do isppol=1,nsppol
       open(unit=iun(isppol),file=filew90_eig(isppol),form='formatted',status='unknown')
     end do
   end if !rank==master
!  Loop to write eigenvalues
   band_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
       jband=0
       do iband=1,mband
         if(band_in(iband,isppol)) then
           jband=jband+1
!          Writing data
           if(rank==master) write(iun(isppol), '(2i6,4x,f10.5)' ) jband,ikpt,Ha_eV*eigen(iband+band_index)
!          Finish writing, now save eigenvalues
           eigenvalues_w(jband,ikpt,isppol)=Ha_eV*eigen(iband+band_index)
         end if
       end do !iband
       band_index=band_index+nband_k
     end do !ikpt
   end do  !nsppol
   if(rank==master) then
     do isppol=1,nsppol
       close(iun(isppol))
     end do
     write(message, '(a,a)' ) ch10,&
&     '   mlwfovlp :  eigenvalues written'
     call wrtout(std_out,  message,'COLL')
   end if !master
 end if !leig
!else if( leig . and. lwannierun) then
!read .eig file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!4) Calculate overlaps (file seed_name.mmn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!First calculate indices and shift
!
!write(std_out,*) "Computes shift for cg"
 write(message, '(a,a)' ) ch10,&
& '   mlwfovlp : compute shifts for g-points '
 call wrtout(std_out,  message,'COLL')
!----------------------------------------------------------------------
!Compute shifts for g points (icg,iwav)
!(here mband is not used, because shifts are internal variables of abinit)
!----------------------------------------------------------------------
!write(std_out,*) mpw*dtset%nspinor*mband*mkmem*nsppol
 ABI_ALLOCATE(icg,(nsppol,nkpt))
 icg=0
 icgtemp=0
 iwav(:,:,:)=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
!
!    MPI:cycle over k-points not treated by this node
!

     if (nprocs>1 ) then !sometimes we can have just one processor
       if ( ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-rank)  /=0) CYCLE
     end if

!    write(std_out,*)'rank',rank,'ikpt',ikpt,'isppol',isppol
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
!    write(std_out,*) ikpt+(isppol-1)*nkpt,nkpt
     npw_k=npwarr(ikpt)
     do iband=1,nband_k
       if(iband.gt.mband) then
         write(message,'(a,3i0)')" mband",iband,mband,nband_k
         ABI_ERROR(message)
       end if
       iwav(iband,ikpt,isppol)= &
&       (iband-1)*npw_k*dtset%nspinor+icgtemp
     end do ! iband
     icgtemp=icgtemp+ npw_k*dtset%nspinor*nband_k
!    icg(isppol,ikpt)=icgtemp
!    write(std_out,*) "icg", isppol,ikpt,icg(isppol,ikpt)
   end do  ! ikpt
 end do   ! isppol
!write(std_out,*) "shift for cg computed"
 ABI_DEALLOCATE(icg)
!
!Shifts computed.
!
 if( lmmn) then
!
!  In case of parallelization write out cg for all k-points
!
   if (nprocs > 1) then
!
     if(prtvol>0) then
       write(message, '(3a)' ) ch10,&
&       '   mlwfovlp :  Creating temporary files with cg and cprj (PAW)',ch10
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end if
!
     do isppol=1,nsppol
       do ikpt=1,nkpt
!
!        MPI:cycle over k-points not treated by this node
!
         if (nprocs>1 ) then !sometimes we can have just one processor
           if ( ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-rank)  /=0) CYCLE
         end if

!        write(std_out,*)'writing kpt ',ikpt,'isppol',isppol,' by node ', rank
         write(wfnname,'(a,I5.5,".",I1)') trim(dtfil%fnametmp_cg),ikpt,isppol
         iun_plot=1000+ikpt+ikpt*(isppol-1)

         open (unit=iun_plot, file=wfnname,form='unformatted')
         npw_k=npwarr(ikpt)
         do iband=1,mband
           do ig=1,npw_k*dtset%nspinor
             write(iun_plot) (cg(i,ig+iwav(iband,ikpt,isppol)),i=1,2)
           end do
         end do
         close(iun_plot)
       end do !ikpt
     end do !isppol
!
!    In the PAW case we also need to write out cprj into files
!
     if(psps%usepaw==1) then
!
!      big loop on atoms, kpts, bands and lmn
!
       ikpt2=0
       do isppol=1,nsppol
         do ikpt=1,nkpt
!
!          MPI:cycle over k-points not treated by this node
!
           if (nprocs>1 ) then !sometimes we can have just one processor
             if ( ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-MPI_enreg%me)  /=0) CYCLE
           end if

           ikpt2=ikpt2+1 !sums just on the k-points treated by this node
!
           write(wfnname,'(a,I5.5,".",I1)') trim(dtfil%fnametmp_cprj),ikpt,isppol
           iun_plot=1000+ikpt
           open (unit=iun_plot, file=wfnname,form='unformatted')
!
           do iband=1,mband*dtset%nspinor
             ig=iband+(ikpt2-1)*mband*dtset%nspinor +(isppol-1)*nkpt*mband*dtset%nspinor !index for cprj(:,ig)
!
             do iatom=1,natom
               itypat=dtset%typat(iatom)
               lmn_size=pawtab(itypat)%lmn_size
!
               do ilmn=1,lmn_size
                 write(iun_plot) (( cprj(iatom,ig)%cp(i,ilmn)),i=1,2)
               end do !ilmn
             end do !iatom
           end do !iband

           close(iun_plot)
         end do !ikpt
       end do !isppol
     end if !usepaw==1

!
!
   end if !MPI nprocs>1
!
!  End of MPI preliminarities
!  Calculate PW contribution of overlaps
!
   ABI_ALLOCATE(cm1,(2,mband,mband,nntot,nkpt,nsppol))
   ! this loops over spin internally
   call mlwfovlp_pw(cg,cm1,g1,iwav,kg,mband,&
&   mkmem,mpi_enreg,mpw,nfft,ngfft,nkpt,nntot,&
&   npwarr,dtset%nspinor,nsppol,ovikp,dtfil%fnametmp_cg)
   write(message, '(a,a)' ) ch10,&
&   '   mlwfovlp : PW part of overlap computed   '
   call wrtout(std_out,  message,'COLL')
!
!  compute PAW Contribution and add it to PW contribution
!
   if(psps%usepaw==1) then
     write(message, '(a,a)' ) ch10,&
&     '** smatrix_pawinit : PAW part of overlap  '
     call wrtout(std_out,  message,'COLL')
     ABI_ALLOCATE(cm2_paw,(2,mband,mband))
     do isppol=1,nsppol
       do ikpt1=1,nkpt
!
!        MPI:cycle over k-points not treated by this node
!
         if (nprocs>1 ) then !sometimes we can have just one processor
           if ( ABS(MPI_enreg%proc_distrb(ikpt1,1,isppol)-rank)  /=0) CYCLE
         end if

         write(message, '(a,i6,a,2i6)' ) &
&         '   processor',rank,' computes PAW part for kpt and spin',ikpt1,isppol
         call wrtout(std_out,  message,'COLL')

         do intot=1,nntot
           ikpt2= ovikp(ikpt1,intot)
           g1temp(:)=g1(:,ikpt1,intot)
           call smatrix_pawinit(atindx1,cm2_paw,cprj,ikpt1,ikpt2,isppol,&
&           g1temp,gprimd,dtset%kpt,mband,mband,mkmem,mpi_enreg,&
&           natom,dtset%nband,nkpt,dtset%nspinor,nsppol,dtset%ntypat,pawang,pawrad,pawtab,rprimd,&
&           dtfil%fnametmp_cprj,dtset%typat,xred)
!          cm1(:,:,:,intot,ikpt1,isppol)=four_pi*cm2_paw(:,:,:)
!           write(6,*) "ikpt1=",ikpt1
!           do iband=1,mband
!             write(6,*) "iband=",iband
!             write(6,*) "Wannier PW       overlap",cm1(:,iband,iband,intot,ikpt1,isppol)
!             write(6,*) "Wannier PAW      overlap",four_pi*cm2_paw(:,iband,iband)
!             write(6,*) "Wannier PW+PAW   overlap",cm1(:,iband,iband,intot,ikpt1,isppol)+four_pi*cm2_paw(:,iband,iband)
!           enddo
           cm1(:,:,:,intot,ikpt1,isppol)=cm1(:,:,:,intot,ikpt1,isppol)+four_pi*cm2_paw(:,:,:)
         end do ! intot
       end do ! ikpt1
     end do ! isppol
     ABI_DEALLOCATE(cm2_paw)
     write(message, '(a,a)' ) ch10,&
&     '   mlwfovlp : PAW part of overlap computed '
     call wrtout(std_out,  message,'COLL')
   end if ! usepaw
!
   call xmpi_barrier(spaceComm)
   call xmpi_sum(cm1,spaceComm,ierr)
!
!  write overlap for separate calculation of wannier functions
!
   if(rank==master) then
     do isppol=1,nsppol !we write separate output files for each isppol
       iun(isppol)=220+isppol
       open(unit=iun(isppol),file=filew90_mmn(isppol),form='formatted',status='unknown')
       write(iun(isppol),*) "nnkp version 90"
       write(iun(isppol),*) num_bands(isppol),nkpt,nntot
     end do
   end if ! rank==master

   do isppol=1,nsppol
     do ikpt1=1,nkpt
       do intot=1,nntot
         if( rank==master) write(iun(isppol),'(2i6,3x,3x,3i5)') ikpt1,ovikp(ikpt1,intot),(g1(jj,ikpt1,intot),jj=1,3)
         jband2=0
         do iband2=1,mband ! the first index is faster
           if(band_in(iband2,isppol)) then
             jband2=jband2+1
             jband1=0
             do iband1=1,mband
               if(band_in(iband1,isppol)) then
                 jband1=jband1+1
                 if(rank==master) write(iun(isppol),*) &
&                 cm1(1,iband1,iband2,intot,ikpt1,isppol),cm1(2,iband1,iband2,intot,ikpt1,isppol)
                 M_matrix(jband1,jband2,intot,ikpt1,isppol)=&
&                 cmplx(cm1(1,iband1,iband2,intot,ikpt1,isppol),cm1(2,iband1,iband2,intot,ikpt1,isppol))
!                write(2211,*) ikpt1,intot,iband1,iband2
!                write(2211,*) cm1(1,iband1,iband2,intot,ikpt1,isppol),cm1(2,iband1,iband2,intot,ikpt1,isppol)
               end if ! band_in(iband1)
             end do ! iband1
           end if ! band_in(iband2)
         end do ! iband2
       end do !intot
     end do !ikpt
     if( rank==master ) then
       close(iun(isppol))
       write(message, '(3a)' )  '   ',trim(filew90_mmn(isppol)),' written'
       call wrtout(std_out,  message,'COLL')
     end if !rank==master
   end do !isppol
!
   ABI_DEALLOCATE(cm1)
!
!  Write down part of the matrix to the output file
!  This is for the automatic tests
!
   if(rank==master) then
     write(message, '(4a)' ) ch10,&
&     '   Writing top of the overlap matrix: M_mn(ikb,ik)',ch10,&
&     '   m=n=1:3, ikb=1, ik=1'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
!
!    just write down the first 3 elements
!
     do isppol=1,nsppol
       write(message, '( " " )')
       if (nsppol>1 ) then
         if (isppol==1) write(message,'(2a)')trim(message),'   spin up:'
         if (isppol==2) write(message,'(2a)')trim(message),'   spin down:'
       end if
       do ii=1,3
         if(ii>num_bands(isppol)) cycle
         write(message,'(3a)') trim(message),ch10,';   ( '
         do jj=1,3
           if(jj>num_bands(isppol))cycle
           write(message, '(a,2f11.6,a)') trim(message),&
&           M_matrix(ii,jj,1,1,isppol),' , '
         end do
         write(message,'(2a)') trim(message),'    ) '
       end do
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end do
!
!    Now write down bottom of the matrix
!
     write(message, '(4a)' ) ch10,&
&     '   Writing bottom of the overlap matrix: M_mn(ikb,ik)',ch10,&
&     '   m=n=num_bands-2:num_bands, ikb=nntot, ik=nkpt'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
!
     do isppol=1,nsppol
       write(message, '( " " )')
       if (nsppol>1 ) then
         if (isppol==1) write(message,'(2a)')trim(message),'   spin up:'
         if (isppol==2) write(message,'(2a)')trim(message),'   spin down:'
       end if
       do ii=num_bands(isppol)-2,num_bands(isppol)
         if(ii<1) cycle
         write(message,'(3a)') trim(message),ch10,';   ( '
         do jj=num_bands(isppol)-2,num_bands(isppol)
           if(jj<1)cycle
           write(message, '(a,2f11.6,a)') trim(message),&
&           M_matrix(ii,jj,nntot,nkpt,isppol),' , '
         end do !j
         write(message,'(2a)') trim(message),'    ) '
       end do !ii
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end do !isppol
   end if !rank==master
!
!  erase temporary files created for parallel runs
!
   if (nprocs > 1) then
!
     if(prtvol>0) then
       write(message, '(3a)' ) ch10,&
&       '   mlwfovlp :  Removing temporary files with cg and cprj (PAW)',ch10
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end if
!
!    Just master  node will remove the files
!
     if(rank==master) then
       do isppol=1,nsppol
         do ikpt=1,nkpt
           write(wfnname,'(a,I5.5,".",I1)') trim(dtfil%fnametmp_cg),ikpt,isppol
           call delete_file(wfnname,ierr)
           if(psps%usepaw==1) then
             write(wfnname,'(a,I5.5,".",I1)') trim(dtfil%fnametmp_cprj),ikpt,isppol
             call delete_file(wfnname,ierr)
           end if
         end do !ikpt
       end do !isppol
     end if
   end if !MPI nprocs>1
!
 end if !lmmn
!if ( lmmn== .false. .and. lwannierun ) the
!read .mmn file
!
!Deallocate arrays  no longer used
!
 ABI_DEALLOCATE(ovikp)
 ABI_DEALLOCATE(g1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!5) Calculate initial projections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if(dtset%w90iniprj/=0 )  then


!
!  Set value for lproj (type of projections to be computed)
!  In PAW, options 5 and 6 are not in use.
!  5 means that there will be a contribution from inside the spheres and another from the PW part
!  6 means that we take into account just the inside-spheres contribution
!  2 means that PW part will be calculated

!
   lproj=dtset%w90iniprj
   if(dtset%w90iniprj == 5 ) lproj=2 ! Necessary to calculate PW contribution
!
   ABI_ALLOCATE(just_augmentation,(mwan,nsppol))
   just_augmentation(:,:)=.false.

   if( psps%usepaw==1 .and. (dtset%w90iniprj==2 .or. dtset%w90iniprj>4)) then
     if (dtset%w90iniprj==6) just_augmentation(:,:)=.true.
     if (dtset%w90iniprj==5) then
       do isppol=1,nsppol
         do iwan=1,nwan(isppol)
!
!          Trick to skip the planewave contribution for some Wannier functions
!          (Not in production).
!
           if(proj_radial(iwan,isppol) > 4) then
             just_augmentation(iwan,isppol)=.true.
             proj_radial(iwan,isppol)=proj_radial(iwan,isppol)-3
             write(message, '(2a,2i4)' ) &
&             '   ','Skiping planewave contribution for iwan, ispin=',iwan,isppol
             call wrtout(std_out,  message,'COLL')
           end if !proj_radial>4
         end do !iwan
       end do !isppol
     end if !w90iniprj == 5
   end if !paw
!
!  Call mlwfovlp_proj (plane waves part of projections)
!
   if (dtset%w90iniprj/=6) then ! option 6 not yet in use
     call mlwfovlp_proj(A_matrix,band_in,cg,cprj,dtset,gprimd,just_augmentation,kg,&
&     lproj,max_num_bands,mband,mkmem,mpi_enreg,mpw,mwan,natom,&
&     nattyp,nkpt,npwarr,&
&     dtset%nspinor,nsppol,ntypat,num_bands,nwan,pawtab,proj_l,proj_m,&
&     proj_radial,proj_site,proj_x,proj_z,proj_zona,psps,ucvol)
     write(message, '(a,a,a,a)' ) ch10,&
&     '   mlwfovlp:  mlwfovlp_proj done -',ch10,&
&     '   Projectors computed.'
     call wrtout(std_out,  message,'COLL')
   end if !w90proj/=6
!
!  Calculate inside-sphere part of projections (PAW)
!
   if (psps%usepaw ==1 .and. ( dtset%w90iniprj>4)) then
     ABI_ALLOCATE(A_paw,(max_num_bands,mwan,nkpt,nsppol))
     call mlwfovlp_projpaw(A_paw,band_in,cprj,just_augmentation,max_num_bands,mband,mkmem,&
&     mwan,natom,dtset%nband,nkpt,&
&     dtset%nspinor,nsppol,dtset%ntypat,nwan,pawrad,pawtab,&
&     proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,psps,&
&     rprimd,dtset%typat,xred)
!
     write(message, '(a,a,a,a)' ) ch10,&
&     '   mlwfovlp:  mlwfovlp_proj_paw done -',ch10,&
&     '   Inside-spheres part of projectors computed.'
     call wrtout(std_out,  message,'COLL')
!
!    Add in-sphere contribution to A_matrix
!
!
!    w90iniprj==5. Plane waves + augmentation contributions
!
     if(dtset%w90iniprj==5) A_matrix(:,:,:,:)=A_matrix(:,:,:,:)+A_paw(:,:,:,:)
!
!    w90iniprj==6. Just augmentation contribution
!
     if(dtset%w90iniprj==6) A_matrix(:,:,:,:)=A_paw(:,:,:,:)
!
!    deallocations
!
     ABI_DEALLOCATE(A_paw)
   end if !usepaw==1

   ABI_DEALLOCATE(just_augmentation)
!
   call xmpi_barrier(spaceComm)
   call xmpi_sum(A_matrix,spaceComm,ierr)

!
!  write      projections  to a file
!
   if(rank==master) then
     if(dtset%w90iniprj==1) then
       do isppol=1,nsppol
         iun(isppol)=219+isppol
         open(unit=iun(isppol),file=trim(filew90_ramn(isppol)),form='formatted',status='unknown')
         write(iun(isppol),*) 'Projections from Abinit : mband,nkpt,nwan. indices: iband1,iwan,ikpt'
         write(iun(isppol),*) num_bands(isppol),nkpt,nwan(isppol)
       end do
     else
       do isppol=1,nsppol
         iun(isppol)=220+isppol
         open(unit=iun(isppol),file=trim(filew90_amn(isppol)),form='formatted',status='unknown')
         write(iun(isppol),*) 'Projections from Abinit : mband,nkpt,nwan. indices: iband1,iwan,ikpt'
         write(iun(isppol),*) num_bands(isppol),nkpt,nwan(isppol)
       end do
     end if
!
     do isppol=1,nsppol
       do ikpt=1,nkpt
         do iwan=1,nwan(isppol)
           jband=0
           do iband=1,mband
             if(band_in(iband,isppol)) then
               jband=jband+1
               write(iun(isppol),'(3i6,13x,3x,2f18.14)')jband,iwan,ikpt,A_matrix(jband,iwan,ikpt,isppol)
             end if !band_in
           end do !iband
         end do !iwan
       end do !ikpt
     end do !isppol
!
     if(dtset%w90iniprj==1) then
       do isppol=1,nsppol
         close(iun(isppol))
         write(message, '(3a)' ) &
&         '   ',trim(filew90_ramn(isppol)),' written'
         call wrtout(std_out,  message,'COLL')
       end do
     else
       do isppol=1,nsppol
         close(iun(isppol))
         write(message, '(3a)' ) &
&         '   ',trim(filew90_amn(isppol)),' written'
         call wrtout(std_out,  message,'COLL')
       end do
     end if
   end if !rank==master

!
!
!  Write down part of the matrix to the output file
!  This is for the automatic tests
!
   if(rank==master) then
     write(message, '(4a)' ) ch10,&
&     '   Writing top of the initial projections matrix: A_mn(ik)',ch10,&
&     '   m=1:3, n=1:3, ik=1'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
!
!    just write down the first 3 elements
!
     do isppol=1,nsppol
       write(message, '( " " )')
       if (nsppol>1 ) then
         if (isppol==1) write(message,'(2a)')trim(message),'   spin up:'
         if (isppol==2) write(message,'(2a)')trim(message),'   spin down:'
       end if
       do ii=1,3
         if(ii>num_bands(isppol)) cycle
         write(message,'(3a)') trim(message),ch10,';   ( '
         do jj=1,3
           if(jj>nwan(isppol))cycle
           write(message, '(a,2f11.6,a)') trim(message),&
&           A_matrix(ii,jj,1,isppol),' , '
         end do
         write(message,'(2a)') trim(message),'    ) '
       end do
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end do
!
!    Now write down bottom of the matrix
!
     write(message, '(4a)' ) ch10,&
&     '   Writing bottom of the initial projections matrix: A_mn(ik)',ch10,&
&     '   m=num_bands-2:num_bands, n=nwan-2:nwan, ik=nkpt'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
!
     do isppol=1,nsppol
       write(message, '( " " )')
       if (nsppol>1 ) then
         if (isppol==1) write(message,'(2a)')trim(message),'   spin up:'
         if (isppol==2) write(message,'(2a)')trim(message),'   spin down:'
       end if
       do ii=num_bands(isppol)-2,num_bands(isppol)
         if(ii<1) cycle
         write(message,'(3a)') trim(message),ch10,';   ( '
         do jj=nwan(isppol)-2,nwan(isppol)
           if(jj<1)cycle
           write(message, '(a,2f11.6,a)') trim(message),&
&           A_matrix(ii,jj,nkpt,isppol),' , '
         end do
         write(message,'(2a)') trim(message),'    ) '
       end do
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end do !isppol
   end if !rank==master
 end if !dtset%w90iniprj/=0
!
!Deallocations
!
 ABI_DEALLOCATE(proj_site)
 ABI_DEALLOCATE(proj_l)
 ABI_DEALLOCATE(proj_m)
 ABI_DEALLOCATE(proj_radial)
 ABI_DEALLOCATE(proj_x)
 ABI_DEALLOCATE(proj_z)
 ABI_DEALLOCATE(proj_zona)
 ABI_DEALLOCATE(proj_s_loc)
 ABI_DEALLOCATE(proj_s_qaxis_loc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!6) write files for wannier function plot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if( dtset%w90prtunk>0) then
   if(psps%usepaw==1) then
     write(message, '( a,a,a,a,a,a,a,a,a)')ch10,&
&     "   WARNING: The UNK matrices will not contain the correct wavefunctions ",ch10,&
&     "   since we are just writing the plane wave contribution.",ch10,&
&     "   The contribution from inside the spheres is missing. ",ch10,&
&     "   However, these files can be used for plotting purposes",ch10
     call wrtout(std_out,  message,'COLL')
   end if
!
   spacing = dtset%w90prtunk
   write(message, '( 8a,i3,2a)')ch10,&
&   "   UNK files will be written.",ch10,&
&   "   According to the chosen value of w90prtunk",ch10,&
&   "   the wavefunctions are to be written ",ch10, &
&   "   at every ", spacing," records.",ch10
   call wrtout(std_out,  message,'COLL')
!
   ABI_ALLOCATE(kg_k,(3,mpw))
   n1=ngfft(1)
   n2=ngfft(2)
   n3=ngfft(3)
   n4=ngfft(4)
   n5=ngfft(5)
   n6=ngfft(6)
   cplex=1
   mgfft=mgfftc ! error
   do isppol=1,nsppol
     ikg=0
     do ikpt=1,nkpt
!
!      MPI:cycle over k-points not treated by this node
!
       if (nprocs>1 ) then !sometimes we can have just one processor
         if ( ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-rank)  /=0) CYCLE
       end if
!
       npw_k=npwarr(ikpt)
       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       ABI_ALLOCATE(denpot,(cplex*n4,n5,n6))
       ABI_ALLOCATE(cwavef,(2,npw_k))
       ABI_ALLOCATE(fofr,(2,n4,n5,n6))
       ABI_ALLOCATE(gbound,(2*mgfft+8,2))
       ABI_ALLOCATE(fofgout,(2,npw_k))
       iun_plot=1000+ikpt+ikpt*(isppol-1)
       write(wfnname,'("UNK",I5.5,".",I1)') ikpt, isppol
!      open (unit=iun_plot, file=wfnname,form='formatted')
       open(unit=iun_plot, file=wfnname,form='unformatted')
!      optimizing grid for UNK files
       n1tmp = n1/spacing
       n2tmp = n2/spacing
       n3tmp = n3/spacing
       if( mod(n1,spacing) /= 0) then
         n1tmp = n1tmp + 1
       end if
       if( mod(n2,spacing) /= 0) then
         n2tmp = n2tmp + 1
       end if
       if( mod(n3,spacing) /= 0) then
         n3tmp = n3tmp + 1
       end if
!      write(iun_plot,*) n1tmp,n2tmp,n3tmp,ikpt,nband_inc
       write(iun_plot) n1tmp,n2tmp,n3tmp,ikpt,nband_inc(isppol)
!      gbound=zero
       call sphereboundary(gbound,dtset%istwfk(ikpt),kg_k,mgfft,npw_k)
       write(std_out,*) "  writes UNK file for ikpt, spin=",ikpt,isppol
       denpot(:,:,:)=zero
       weight = one
       do iband=1,mband
         if(band_in(iband,isppol)) then
           do ig=1,npw_k*dtset%nspinor
             cwavef(1,ig)=cg(1,ig+iwav(iband,ikpt,isppol))
             cwavef(2,ig)=cg(2,ig+iwav(iband,ikpt,isppol))
           end do
           tim_fourwf=0
           call fourwf(cplex,denpot,cwavef,fofgout,fofr,&
&           gbound,gbound,dtset%istwfk(ikpt),kg_k,kg_k,mgfft,&
&           mpi_enreg,1,ngfft,npw_k,npw_k,n4,n5,n6,0,&
&           tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)
!          do jj3=1,n3,spacing
!          do jj2=1,n2,spacing
!          do jj1=1,n1,spacing
!          write(iun_plot,*) fofr(1,jj1,jj2,jj3),&
!          & fofr(2,jj1,jj2,jj3)
!          end do !jj1
!          end do !jj2
!          end do !jj3
!          unformatted (must be one record)
           write(iun_plot) (((fofr(1,jj1,jj2,jj3),fofr(2,jj1,jj2,jj3),&
&           jj1=1,n1,spacing),jj2=1,n2,spacing),jj3=1,n3,spacing)
         end if !iband
       end do ! iband
       ABI_DEALLOCATE(cwavef)
       ABI_DEALLOCATE(fofr)
       ABI_DEALLOCATE(gbound)
       ABI_DEALLOCATE(denpot)
       ABI_DEALLOCATE(fofgout)
       ikg=ikg+npw_k
       close(iun_plot)
     end do  ! ikpt
   end do  ! nsppol
   ABI_DEALLOCATE(kg_k)
!
   write(message, '(4a)' )ch10, &
&   '   ','UNK files written',ch10
   call wrtout(std_out,  message,'COLL')
 end if !dtset%w90prtunk
!
 ABI_DEALLOCATE(iwav)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!7) Call to  Wannier90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if(lwannierrun) then
   if(lwanniersetup.ne.1) ABI_ERROR("lwanniersetup.ne.1")
   ABI_ALLOCATE(U_matrix,(mwan,mwan,nkpt,nsppol))
   ABI_ALLOCATE(U_matrix_opt,(max_num_bands,mwan,nkpt,nsppol))
   ABI_ALLOCATE(lwindow,(max_num_bands,nkpt,nsppol))
   ABI_ALLOCATE(wann_centres,(3,mwan,nsppol))
   ABI_ALLOCATE(wann_spreads,(mwan,nsppol))
!  Initialize
   U_matrix(:,:,:,:)=czero
   U_matrix_opt(:,:,:,:)=czero
   lwindow(:,:,:)=.false.
   wann_centres(:,:,:)=zero
   wann_spreads(:,:)=zero
!
!  write(std_out,*) seed_name
!  write(std_out,*) ngkpt
   ngkpt(1)=dtset%kptrlatt(1,1)
   ngkpt(2)=dtset%kptrlatt(2,2) ! ajouter test de verif que kptrlatt est bien diagonal
   ngkpt(3)=dtset%kptrlatt(3,3)
!  write(std_out,*) nkpt
!  write(std_out,*) rprimd*Bohr_Ang
!  write(std_out,*) two_pi*gprimd/Bohr_Ang
!  write(std_out,*) mband
!  write(std_out,*) "nwan",nwan
!  write(std_out,*) nntot
!  write(std_out,*) natom
!  write(std_out,*) atom_symbols
!  write(std_out,*) xcart
!  write(std_out,*) num_bands,num_bands,nntot,nkpt
!  write(std_out,*) wann_spreads
!  wann_spreads=2
!  do i=1, nkpt
!  do j=1, nntot
!  write(std_out,*) i,j
!  do k=1, num_bands
!  do l=1, num_bands
!  write(std_out,*) "m",M_matrix(l,k,j,i,1)
!  enddo
!  enddo
!  enddo
!  enddo



#if defined HAVE_WANNIER90
   do isppol=1,nsppol
!    when nsppol>1, master runs isppol 1 and rank==1 runs isppol 2
     if(nprocs>1 .and. isppol==1.and.rank.ne.master) cycle
     if(nprocs>1 .and. isppol==2.and.rank.ne.1) cycle

     write(message, '(8a)' ) ch10,&
&     '** mlwfovlp :   call wannier90 library subroutine wannier_run ',ch10,&
&     '   Calculation is running         ',ch10,&
&     '-  see ',trim(filew90_wout(isppol)),' for details.'
     call wrtout(std_out,  message,'COLL')
!
     call wannier_run(trim(seed_name(isppol)),ngkpt,nkpt,&            !input
&    real_lattice,recip_lattice,dtset%kpt,num_bands(isppol),& !input
&    nwan(isppol),nntot,natom,atom_symbols,&                  !input
&    xcart*Bohr_Ang,gamma_only,M_matrix(:,:,:,:,isppol),A_matrix(:,:,:,isppol),eigenvalues_w(:,:,isppol),& !input
&    U_matrix(1:nwan(isppol),1:nwan(isppol),:,isppol),& !output
&    U_matrix_opt(1:num_bands(isppol),1:nwan(isppol),:,isppol),& !output
&    lwindow_loc=lwindow(1:num_bands(isppol),:,isppol),& !output
&    wann_centres_loc=wann_centres(:,1:nwan(isppol),isppol),&     !output
&    wann_spreads_loc=wann_spreads(1:nwan(isppol),isppol),spread_loc=spreadw(:,isppol))                            !output

!    ----------------------------------------------------------------------------------------------

     write(message, '(7a)' ) ch10,&
&     '   mlwfovlp :  mlwfovlp_run completed -',ch10,&
&     '-  see ',trim(filew90_wout(isppol)),' for details.',ch10
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

   end do !isppol

!  collect output of  wannier90 from different processors
   call xmpi_barrier(spaceComm)

   call xmpi_sum(U_matrix,spaceComm,ierr)
   call xmpi_sum(U_matrix_opt,spaceComm,ierr)
   call xmpi_lor(lwindow,spaceComm)
   call xmpi_sum(wann_centres,spaceComm,ierr)
   call xmpi_sum(wann_spreads,spaceComm,ierr)

   ! Output ABIWAN.nc file
#ifdef HAVE_NETCDF
   if (dtset%kptopt == 0) then
     ABI_WARNING("Output of ABIWAN.nc requires kptopt /= 0. ABIWAN.nc file won't be produced!")
     ! Need kptrlatt in wigner_seitz and client code need to know the k-grid.
   end if
   if (rank == master .and. dtset%kptopt /= 0) then
     abiwan_fname = strcat(dtfil%filnam_ds(4), "_ABIWAN.nc")
     call wrtout(std_out, sjoin("Saving wannier90 ouput results in:", abiwan_fname))
     call wigner_seitz([zero, zero, zero], [2, 2, 2], dtset%kptrlatt, crystal%rmet, nrpts, irvec, ndegen)
     ! We know if disentanglement has been done by looking at the output values of lwindow
     ! Not elegant but it is the only way to avoid the parsing of the wannier input.
     ! In wannier_run lwindow is set to True if not disentanglement
     have_disentangled_spin = 0
     do isppol=1,nsppol
       !if nwan(isppol) < num_bands(isppol)
       if (.not. all(lwindow(:,:,isppol))) have_disentangled_spin(isppol) = 1
     end do

     NCF_CHECK(nctk_open_create(ncid, abiwan_fname, xmpi_comm_self))
     NCF_CHECK(hdr%ncwrite(ncid, fform_from_ext("ABIWAN"), nc_define=.True.))
     NCF_CHECK(crystal%ncwrite(ncid))
     NCF_CHECK(ebands_ncwrite(ebands, ncid))

     ncerr = nctk_def_dims(ncid, [ &
       nctkdim_t("mwan", mwan), &
       nctkdim_t("max_num_bands", max_num_bands), &
       nctkdim_t("nrpts", nrpts) &
     ], defmode=.True.)
     NCF_CHECK(ncerr)

     ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "nntot"])
     NCF_CHECK(ncerr)
     !ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "fermi_energy", "smearing_width"])
     !NCF_CHECK(ncerr)

     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t("nwan", "int", "number_of_spins"), &
       nctkarr_t("num_bands", "int", "number_of_spins"), &
       nctkarr_t("band_in_int", "int", "max_number_of_states, number_of_spins"), &
       nctkarr_t("lwindow_int", "int", "max_num_bands, number_of_kpoints, number_of_spins"), &
       !nctkarr_t("exclude_bands", "int", "max_number_of_states, number_of_spins"), &
       !nctkarr_t("eigenvalues_w", "int", "max_num_bands, number_of_kpoints, number_of_spins"), &
       nctkarr_t("spread", "dp", "three, number_of_spins"), &
       !nctkarr_t("A_matrix", "dp", "two, max_num_bands, mwan, number_of_kpoints, number_of_spins"), &
       nctkarr_t("irvec", "int", "three, nrpts"), &
       nctkarr_t("ndegen", "int", "nrpts"), &
       nctkarr_t("have_disentangled_spin", "int", "number_of_spins"), &
       nctkarr_t("U_matrix", "dp", "two, mwan, mwan, number_of_kpoints, number_of_spins"), &
       nctkarr_t("U_matrix_opt", "dp", "two, max_num_bands, mwan, number_of_kpoints, number_of_spins"), &
       nctkarr_t("wann_centres", "dp", "three, mwan, number_of_spins"), &
       nctkarr_t("wann_spreads", "dp", "mwan, number_of_spins") &
       ])
     NCF_CHECK(ncerr)

     ! Write data.
     NCF_CHECK(nctk_set_datamode(ncid))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nntot"), nntot))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nwan"), nwan))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "num_bands"), num_bands))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "band_in_int"), l2int(band_in)))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "lwindow_int"), l2int(lwindow)))
     !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "exclude_bands"), exclude_bands))
     !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eigenvalues_w"), eigenvalues_w))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "spread"), spreadw))
     !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "A_matrix"), c2r(A_matrix)))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "irvec"), irvec))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ndegen"), ndegen))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "have_disentangled_spin"), have_disentangled_spin))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "U_matrix"), c2r(U_matrix)))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "U_matrix_opt"), c2r(U_matrix_opt)))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "wann_centres"), wann_centres))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "wann_spreads"), wann_spreads))
     NCF_CHECK(nf90_close(ncid))

     ABI_FREE(irvec)
     ABI_FREE(ndegen)
   end if
#endif

!  CALL SILVESTRELLI'S APPROACH TO EVALUATE vdW INTERACTION ENERGY USING MLWF!!
!  ----------------------------------------------------------------------------------------------
   if (dtset%vdw_xc==10.or.dtset%vdw_xc==11.or.dtset%vdw_xc==12.or.dtset%vdw_xc==14.and.rank==master) then
!    vdw_xc==10,11,12,14 starts the
!    vdW interaction using MLWFs
     write(std_out,*) 'nwan(nsppol)=',ch10
     do ii=1,nsppol
       write(std_out,*) 'nsppol=',ii, 'nwan(nsppol)=',nwan(ii),ch10
     end do
     write(std_out,*) 'mwan=', mwan, ch10

     ABI_ALLOCATE(occ_arr,(mband,nkpt,isppol))
     ABI_ALLOCATE(occ_wan,(mwan,nkpt,nsppol))
     ABI_ALLOCATE(tdocc_wan,(mwan,nsppol))

     occ_arr(:,:,:)=zero
     occ_wan(:,:,:)=zero
     tdocc_wan(:,:)=zero
     jj = 0
     do isppol=1,nsppol
       do ikpt=1,nkpt
         do iband=1,num_bands(isppol)
           jj = jj + 1
           occ_arr(iband,ikpt,isppol) = occ(jj)
         end do
       end do
     end do

     do isppol=1,nsppol
       do ikpt=1,nkpt
         do iwan=1,nwan(isppol)
           caux=czero
           caux2=czero
           caux3=czero
           do iband=1,num_bands(isppol) !nband_inc(isppol) !nwan(isppol)
             do ii=1,nwan(isppol)
               caux=U_matrix(ii,iwan,ikpt,isppol)*U_matrix_opt(iband,ii,ikpt,isppol)
!              DEBUG
!              if(ISNAN(dble(caux))) then
!              write(std_out,*) 'NaN: caux(ikpt,iwan,iband,ii):',ikpt,iwan,iband,ii,ch10
!              end if
!              END DEBUG
               do kk=1,nwan(isppol)
                 caux2=conjg(U_matrix(kk,iwan,ikpt,isppol))*conjg(U_matrix_opt(iband,kk,ikpt,isppol))
                 caux3= caux3+caux*caux2*occ_arr(iband,ikpt,isppol) !take care here as exclude_bands case is not well
!                DEBUG
!                if(ISNAN(dble(caux2))) then
!                write(std_out,*) 'NaN: caux2(ikpt,iwan,iband,kk):',ikpt,iwan,iband,kk,ch10
!                end if
!                if(ISNAN(dble(caux3))) then
!                write(std_out,*) 'NaN: caux3(ikpt,iwan,iband,kk,jj):',ikpt,iwan,iband,kk,jj
!                end if
!                END DEBUG
               end do
             end do
           end do
           occ_wan(iwan,ikpt,isppol) = dble(caux3)
!          DEBUG
!          write(std_out,*) occ_wan(iwan,ikpt,isppol)
!          END DEBUG
!          end do
         end do
       end do
     end do

     write(std_out,*) ch10,'MLWFs Occupation Matrix diagonal terms:',ch10

     do jj=1,nsppol
       forall(iwan=1:nwan(jj)) tdocc_wan(iwan,jj) = sum(occ_wan(iwan,1:nkpt,jj)) / real(nkpt,dp)
       write(std_out,*) 'tdocc_wan(iwan),isppol:',ch10
       write(std_out,*) (tdocc_wan(iwan,jj),iwan=1,nwan(jj)),jj
     end do

     ABI_ALLOCATE(csix,(mwan,mwan,nsppol,nsppol))

     call evdw_wannier(csix,corrvdw,mwan,natom,nsppol,nwan,tdocc_wan,dtset%vdw_nfrag,&
&     dtset%vdw_supercell,dtset%vdw_typfrag,dtset%vdw_xc,rprimd,wann_centres,wann_spreads,xcart)

     ABI_DEALLOCATE(csix)
     ABI_DEALLOCATE(occ_arr)
     ABI_DEALLOCATE(occ_wan)
     ABI_DEALLOCATE(tdocc_wan)
   end if
#else
   ABI_UNUSED(occ)
#endif
!  FIXME: looks like there is no automatic test which goes through here: g95 bot did not catch
!  the missing deallocations
   ABI_DEALLOCATE(wann_centres)
   ABI_DEALLOCATE(wann_spreads)
   ABI_DEALLOCATE(U_matrix)
   ABI_DEALLOCATE(U_matrix_opt)
   ABI_DEALLOCATE(lwindow)

 end if !lwannierrun
!
!deallocation
!
 ABI_DEALLOCATE(band_in)
 ABI_DEALLOCATE(atom_symbols)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(eigenvalues_w)
 ABI_DEALLOCATE(M_matrix)
 ABI_DEALLOCATE(A_matrix)

contains
!!***


!!****f* mlwfovlp/read_chkunit
!! NAME
!! read_chkunit
!!
!! FUNCTION
!! Function which reads the .chk file produced by Wannier90
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      initylmr,matrginv,rotmat
!!
!! SOURCE

 subroutine read_chkunit(seed_name,nkpt,ndimwin,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
 character(len=*),intent(in) :: seed_name
 integer,intent(out) :: ierr
!arrays
 integer,intent(out) :: ndimwin(nkpt)

!Local variables-------------------------------
!scalars
 integer :: chk_unit,ios,ikpt
 logical :: have_disentangled

!************************************************************************

   chk_unit=get_unit()
   fname=TRIM(seed_name)//'.chk'
   open(unit=chk_unit,file=fname,form='unformatted',status='old',iostat=ios)

   ierr=0
   read(chk_unit) ! header                                   ! Date and time
   read(chk_unit) ! ((real_lattice(i,j),i=1,3),j=1,3)        ! Real lattice
   read(chk_unit) ! ((recip_lattice(i,j),i=1,3),j=1,3)       ! Reciprocal lattice
   read(chk_unit) ! num_kpts
   read(chk_unit) ! ((kpt_latt(i,nkp),i=1,3),nkp=1,num_kpts) ! K-points
   read(chk_unit) ! nntot                  ! Number of nearest k-point neighbours
   read(chk_unit) ! num_wann               ! Number of wannier functions
   read(chk_unit) ! chkpt1                 ! Position of checkpoint
   read(chk_unit) have_disentangled        ! Whether a disentanglement has been performed
   if (have_disentangled) then
!    read(chk_unit) ! omega_invariant     ! Omega invariant
!    read(chk_unit) ((lwindow(i,nkp),i=1,num_bands),nkp=1,num_kpts)
     read(chk_unit) (ndimwin(ikpt),ikpt=1,nkpt)
!    read(chk_unit) (((u_matrix_opt(i,j,nkp),i=1,num_bands),j=1,num_wann),nkp=1,num_kpts)
   else
!    this is not expected. we should have disentanglement. Report the error.
     ierr=-1
   end if
!  read(chk_unit)  (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),k=1,num_kpts)               ! U_matrix
!  read(chk_unit)  ((((m_matrix(i,j,k,l,1),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts) ! M_matrix
!  read(chk_unit)  ((wannier_centres(i,j),i=1,3),j=1,num_wann)
   close(chk_unit)

end subroutine read_chkunit
!!***

end subroutine mlwfovlp
!!***

!!****f* m_mlwfovlp/mlwfovlp_seedname
!! NAME
!! mlwfovlp_seedname
!!
!! FUNCTION
!! Get seed name and file names of all wannier90 related files
!!
!! INPUTS
!! fname_w90=root name of file appended with _w90
!!
!! OUTPUT
!! filew90_win= main input file for Wannier90
!! filew90_wout= main output file for Wannier90
!! filew90_amn= file containing Amn matrix
!! filew90_ramn= file containing Amn matrix (random initial projections)
!! filew90_mmn= file containing Mmn matrix
!! filew90_eig= file containing eigenvalues
!! nsppol= number of spin polarizations
!! seed_name= common seed name for all wannier90 related files
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_mlwfovlp
!!
!! CHILDREN
!!      initylmr,matrginv,rotmat
!!
!! SOURCE


subroutine mlwfovlp_seedname(fname_w90,filew90_win,filew90_wout,filew90_amn,&
& filew90_ramn,filew90_mmn,filew90_eig,nsppol,seed_name)

!Arguments ------------------------------------
 integer,intent(in) :: nsppol
 character(len=fnlen),intent(out) :: filew90_win(nsppol),filew90_wout(nsppol),filew90_amn(nsppol),filew90_ramn(nsppol)
 character(len=fnlen),intent(out) :: filew90_mmn(nsppol),filew90_eig(nsppol),seed_name(nsppol)
 character(len=fnlen),intent(in) :: fname_w90

!Local variables-------------------------------
 integer::isppol
 character(len=fnlen) :: test_win1,test_win2,test_win3
 logical :: lfile
 character(len=2000) :: message
 character(len=10)::postfix
! *************************************************************************

 seed_name(:)=trim(fname_w90)
 do isppol=1,nsppol
   if(nsppol==1)postfix='.win'
   if(nsppol==2 .and. isppol==1)postfix='_up.win'
   if(nsppol==2 .and. isppol==2)postfix='_down.win'
!
   filew90_win(isppol)=trim(seed_name(isppol))//trim(postfix)
   test_win1=filew90_win(isppol)
   inquire(file=filew90_win(isppol),exist=lfile)

   if(.not.lfile) then
     seed_name(isppol)='wannier90'
     filew90_win(isppol)=trim(seed_name(isppol))//trim(postfix)
     test_win2=filew90_win(isppol)
     inquire(file=filew90_win(isppol),exist=lfile)
   end if

   if(.not.lfile) then
     seed_name(isppol)='w90'
     filew90_win=trim(seed_name(isppol))//trim(postfix)
     test_win3=filew90_win(isppol)
     inquire(file=filew90_win(isppol),exist=lfile)
   end if

   if(.not.lfile) then
     write(message,'(17a)')ch10,&
&     ' mlwfovlp_seedname : ERROR - ',ch10,&
&     ' wannier90 interface needs one of the following files:',ch10,&
&     '      ',trim(test_win1),ch10,&
&     '      ',trim(test_win2),ch10,&
&     '      ',trim(test_win3),ch10,&
&     ' Action: read wannier90 tutorial and/or user manual',ch10,&
&     '  and supply proper *.win file'
     ABI_ERROR(message)
   end if
 end do !isppol


!Files having different names for
!different spin polarizations
 if(nsppol==1) then
   filew90_win(1) =trim(seed_name(1))//'.win'
   filew90_wout(1)=trim(seed_name(1))//'.wout'
   filew90_ramn(1)=trim(seed_name(1))//'random.amn'
   filew90_amn(1) =trim(seed_name(1))//'.amn'
   filew90_mmn(1) =trim(seed_name(1))//'.mmn'
   filew90_eig(1) =trim(seed_name(1))//'.eig'
 elseif(nsppol==2) then
   filew90_win(1) =trim(seed_name(1))//'_up.win'
   filew90_win(2) =trim(seed_name(2))//'_down.win'
!
   filew90_wout(1)=trim(seed_name(1))//'_up.wout'
   filew90_wout(2)=trim(seed_name(2))//'_down.wout'
!
   filew90_ramn(1)=trim(seed_name(1))//'random_up.amn'
   filew90_ramn(2)=trim(seed_name(2))//'random_down.amn'
!
   filew90_amn(1)=trim(seed_name(1))//'_up.amn'
   filew90_amn(2)=trim(seed_name(2))//'_down.amn'
!
   filew90_mmn(1)=trim(seed_name(1))//'_up.mmn'
   filew90_mmn(2)=trim(seed_name(2))//'_down.mmn'
!
   filew90_eig(1)=trim(seed_name(1))//'_up.eig'
   filew90_eig(2)=trim(seed_name(2))//'_down.eig'
 end if
!change also seed_name for nsppol=2
 if(nsppol==2) then
   seed_name(1)=trim(seed_name(1))//'_up'
   seed_name(2)=trim(seed_name(2))//'_down'
 end if
!End file-name section

 write(message, '(a,a)' ) ch10,&
& '---------------------------------------------------------------'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')
 write(message, '(5a)' ) ch10,&
& '  Calculation of overlap and call to wannier90 library ',ch10,&
& '  to obtain maximally localized wannier functions ',ch10

 call wrtout(std_out,  message,'COLL')
 call wrtout(ab_out,message,'COLL')

 if(nsppol==1) then
   write(message, '(23a)' ) &
&   '  - ',trim(filew90_win(1)),' is a mandatory secondary input',ch10,&
&   '  - ',trim(filew90_wout(1)),' is the output for the library',ch10,&
&   '  - ',trim(filew90_ramn(1)),' contains random projections',ch10,&
&   '  - ',trim(filew90_amn(1)),' contains projections',ch10,&
&   '  - ',trim(filew90_mmn(1)),' contains the overlap',ch10,&
&   '  - ',trim(filew90_eig(1)),' contains the eigenvalues'
 elseif(nsppol==2) then
   write(message, '(41a)' ) &
&   '  - ',trim(filew90_win(1)),&
&   ' and ',trim(filew90_win(2)),ch10,'are mandatory secondary input',ch10,&
&   '  - ',trim(filew90_wout(1)),&
&   ' and ',trim(filew90_wout(2)),ch10,' are the output for the library',ch10,&
&   '  - ',trim(filew90_ramn(1)),&
&   ' and ',trim(filew90_ramn(2)),ch10,' contain random projections',ch10,&
&   '  - ',trim(filew90_amn(1)),&
&   ' and ',trim(filew90_amn(2)),ch10,' contain projections',ch10,&
&   '  - ',trim(filew90_mmn(1)),&
&   ' and ',trim(filew90_mmn(2)),ch10,' contain the overlap',ch10,&
&   '  - ',trim(filew90_eig(1)),&
&   ' and ',trim(filew90_eig(2)),ch10,' contain the eigenvalues'
 end if
 call wrtout(std_out,  message,'COLL')
 call wrtout(ab_out,message,'COLL')

 write(message, '(a,a)' ) ch10,&
& '---------------------------------------------------------------'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

end subroutine mlwfovlp_seedname
!!***

!!****f* m_mlwfovlp/mlwfovlp_setup
!! NAME
!! mlwfovlp_setup
!!
!! FUNCTION
!! Routine which creates table g1 and ovikp  necessary to compute
!! overlap for Wannier code (www.wannier.org f90 version).
!!
!! INPUTS
!!  atom_symbols(natom)= table of symbol for each atom
!!                                          and each |p_lmn> non-local projector
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  filew90_win(nsppol) secondary input files for w90
!!  lwanniersetup= flag: only 1 is fully working.
!!  natom              =number of atoms in cell.
!!  mband=maximum number of bands
!!  natom=number of atoms in cell.
!!  nkpt=number of k points.
!!  num_bands(isppol)=number of bands actually used to construct the wannier function
!!  nwan(isppol)= number of wannier fonctions (read in wannier90.win).
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  real_lattice(3,3)=dimensional primitive translations for real space
!!                 in format required by wannier90
!!  recip_lattice(3,3)=dimensional primitive translations for reciprocal space
!!                 in format required by wannier90
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  seed_name=character string for generating wannier90 filenames
!!  xcart(3,natom)=atomic coordinates in bohr
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  band_in(mband,nsppol)   = band to take into account for wannier calculation
!!  g1(3,nkpt,nntot) = G vector shift which is necessary to obtain k1+b
!!                     from k2 in the case where k1+b does not belong to the 1st BZ.
!!  nband_inc(nsppol) = # of included bands
!!  nntot            = number of k-point neighbour
!!  ovikp(nkpt,nntot)= gives  nntot value of k2 (in the BZ) for each k1  (k2=k1+b mod(G))
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      m_mlwfovlp
!!
!! CHILDREN
!!      initylmr,matrginv,rotmat
!!
!! SOURCE

 subroutine mlwfovlp_setup(atom_symbols,band_in,dtset,filew90_win,gamma_only,&
& g1,lwanniersetup,mband,natom,nband_inc,nkpt,&
& nntot,num_bands,num_nnmax,nsppol,nwan,ovikp,&
& proj_l,proj_m,proj_radial,proj_site,proj_s_loc,proj_s_qaxis_loc,proj_x,proj_z,proj_zona,&
& real_lattice,recip_lattice,rprimd,seed_name,spinors,xcart,xred)

!Arguments---------------------------
! scalars
!scalars
 integer,intent(in) :: lwanniersetup,mband,natom,nkpt,nsppol
 integer,intent(in) :: num_nnmax
 integer,intent(out) :: nband_inc(nsppol),nntot,num_bands(nsppol),nwan(nsppol)
 logical,intent(in) :: gamma_only,spinors
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(out) :: g1(3,nkpt,num_nnmax),ovikp(nkpt,num_nnmax)
 integer,intent(out) :: proj_l(mband,nsppol),proj_m(mband,nsppol),proj_radial(mband,nsppol)
 real(dp),intent(in) :: real_lattice(3,3)
 real(dp),intent(in) :: recip_lattice(3,3),rprimd(3,3),xred(3,natom)
 real(dp),intent(out) :: proj_site(3,mband,nsppol),proj_x(3,mband,nsppol),proj_z(3,mband,nsppol)
 real(dp),intent(out) :: proj_zona(mband,nsppol),xcart(3,natom)
 logical,intent(out) :: band_in(mband,nsppol)
 character(len=3),intent(out) :: atom_symbols(natom)
 character(len=fnlen),intent(in) :: seed_name(nsppol),filew90_win(nsppol)

 integer, optional, intent(out) :: proj_s_loc(mband)
 real(dp), optional, intent(out) :: proj_s_qaxis_loc(3,mband)

!Local variables---------------------------
!scalars
 integer :: iatom,icb,ikpt,ikpt1,intot,isppol,itypat,jj,mband_,unt
 real(dp) :: znucl1
 character(len=2) :: symbol
 character(len=500) :: message
 character(len=fnlen) :: filew90_nnkp
 type(atomdata_t) :: atom
!arrays
 integer :: exclude_bands(mband,nsppol),ngkpt(3)

! *************************************************************************

!^^^^^^^^^^^^^^^^read wannier90.nnkp^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 if(lwanniersetup==0) then  !this part is not coded for nsppol>1
   isppol=1
   filew90_nnkp=trim(seed_name(isppol))//'.nnkp'
   if (open_file(filew90_nnkp,message,newunit=unt,form='formatted',status='old') /= 0) then
     ABI_ERROR(message)
   end if
   read(unt,*)
   read(unt,*) nntot , mband_, nwan(1)
   write(message, '(a,a,i6,i6,i6)' )ch10,&
&   ' mlwfovlp_setup nntot,mband,nwan ', nntot,mband_,nwan(1)
   call wrtout(std_out,message,'COLL')
   if(mband_.ne.mband) then
     write(message, '(4a)' )&
&     'mband_ is not equal to mband ',ch10,&
&     'Action: check ',trim(filew90_nnkp)
     ABI_ERROR(message)
   end if
   if(nwan(1)>mband) then
     write(message, '(4a)' )&
&     'nwan > mband ',ch10,&
&     'Action: check ',trim(filew90_nnkp)
     ABI_ERROR(message)
   end if
   if(nwan(1)==0) then
     write(message, '(4a)' )&
&     'nwan = 0 ',ch10,&
&     'Action: check ',trim(filew90_nnkp)
     ABI_ERROR(message)
   end if
   do ikpt=1,nkpt
     do intot=1,nntot
!      ikpt1: k point  (ikpt=ikpt1)
!      ovikp(intot,ikpt): neighbour number intot for ikpt
!      g1(1:3,intot,ikpt): non reciprocal space vector between the 2 k-points
       read(unt,*)  &
&       ikpt1,ovikp(ikpt,intot),(g1(jj,ikpt,intot),jj=1,3)
       if(ikpt1.ne.ikpt) write(std_out,*) "warning: ikpt1 .ne ikpt : ?"
     end do
   end do
   close(unt)
   write(message, '(3a)' )ch10,&
&   trim(filew90_nnkp),'wannier90.nnkp has been read !'
   call wrtout(std_out,message,'COLL')

   message = ' exclude bands is not given in this case (not implemented) '
   ABI_ERROR(message)

!  ^^^^^^^^^^^^^^^^^^^^^^^ call wannier_setup begin^^^^^^^^^^^^^^^^^^^^^^^^
 else if (lwanniersetup==1) then
   num_bands(:)=mband
!  num_nnmax=12 !limit fixed for compact structure in wannier_setup.
   ovikp=0.d0
!  "When nshiftk=1, kptrlatt is initialized as a diagonal (3x3) matrix, whose diagonal
!  elements are the three values ngkpt(1:3)"
   ngkpt(1)=dtset%kptrlatt(1,1)
   ngkpt(2)=dtset%kptrlatt(2,2) !  have to verif kptrlatt is diagonal
   ngkpt(3)=dtset%kptrlatt(3,3)
   do iatom=1,natom
     itypat=dtset%typat(iatom)
     znucl1=dtset%znucl(itypat)
     call atomdata_from_znucl(atom, znucl1)
     symbol=trim(adjustl(atom%symbol))
!    write(309,*) symbol
     atom_symbols(iatom)=symbol
     xcart(:,iatom)=rprimd(:,1)*xred(1,iatom)+&
&     rprimd(:,2)*xred(2,iatom)+&
&     rprimd(:,3)*xred(3,iatom)
   end do ! iatom
!  write(std_out,*) xcart
!  write(std_out,*) Bohr_Ang
!  write(std_out,*) rprimd*Bohr_Ang
!  write(std_out,*) seed_name
!  write(std_out,*) ngkpt
!  write(std_out,*) nkpt
!  write(std_out,*) mband
!  write(std_out,*) natom
!  write(std_out,*) atom_symbols
   write(message, '(a,a)' )ch10,&
&   '** mlwfovlp_setup:  call wannier90 library subroutine wannier_setup'
   call wrtout(std_out,message,'COLL')
#if defined HAVE_WANNIER90
   nwan(:)=0
   num_bands(:)=0
   do isppol=1,nsppol
#ifdef HAVE_WANNIER90_V1
       call wannier_setup(seed_name(isppol),ngkpt,nkpt&            !input
&      ,real_lattice,recip_lattice,dtset%kpt&                      !input
&      ,mband,natom,atom_symbols,xcart*Bohr_Ang&                   !input
&      ,gamma_only,spinors&                                        !input
&      ,nntot,ovikp,g1,num_bands(isppol),nwan(isppol)&             !output
&      ,proj_site(:,:,isppol),proj_l(:,isppol)&                    !output
&      ,proj_m(:,isppol),proj_radial(:,isppol)&                    !output
&      ,proj_z(:,:,isppol),proj_x(:,:,isppol)&                     !output
&      ,proj_zona(:,isppol),exclude_bands(:,isppol))               !output
#else
!WANNIER90_V2 has the 2 optional arguments
     if (present(proj_s_loc)) then
       call wannier_setup(seed_name(isppol),ngkpt,nkpt&            !input
&      ,real_lattice,recip_lattice,dtset%kpt&                      !input
&      ,mband,natom,atom_symbols,xcart*Bohr_Ang&                   !input
&      ,gamma_only,spinors&                                        !input
&      ,nntot,ovikp,g1,num_bands(isppol),nwan(isppol)&             !output
&      ,proj_site(:,:,isppol),proj_l(:,isppol)&                    !output
&      ,proj_m(:,isppol),proj_radial(:,isppol)&                    !output
&      ,proj_z(:,:,isppol),proj_x(:,:,isppol)&                     !output
&      ,proj_zona(:,isppol),exclude_bands(:,isppol)&               !output
&      ,proj_s_loc,proj_s_qaxis_loc)                               !output
     else
!no proj_s_loc provided
       call wannier_setup(seed_name(isppol),ngkpt,nkpt&            !input
&      ,real_lattice,recip_lattice,dtset%kpt&                      !input
&      ,mband,natom,atom_symbols,xcart*Bohr_Ang&                   !input
&      ,gamma_only,spinors&                                        !input
&      ,nntot,ovikp,g1,num_bands(isppol),nwan(isppol)&             !output
&      ,proj_site(:,:,isppol),proj_l(:,isppol)&                    !output
&      ,proj_m(:,isppol),proj_radial(:,isppol)&                    !output
&      ,proj_z(:,:,isppol),proj_x(:,:,isppol)&                     !output
&      ,proj_zona(:,isppol),exclude_bands(:,isppol))               !output
     end if
#endif
   end do !isppol
! if we do not have w90, avoid complaints about unused input variables
#else
   ABI_UNUSED(gamma_only)
   ABI_UNUSED(real_lattice)
   ABI_UNUSED(recip_lattice)
   ABI_UNUSED(spinors)
#endif
!  do isppol=1,nsppol
!  write(std_out,*)  "1", nntot,nwan(isppol)
!  write(std_out,*)  "2",num_bands(isppol)  ! states on which wannier functions are computed
!  write(std_out,*)  "3", proj_site(:,1:nwan(isppol),isppol)
!  write(std_out,*)  "4",proj_l(1:nwan(isppol),isppol)
!  write(std_out,*)  "5",proj_m(1:nwan(isppol),isppol)
!  write(std_out,*)  "6", proj_radial(1:nwan(isppol),isppol)
!  write(std_out,*)  "7", proj_z(:,1:nwan(isppol),isppol)
!  write(std_out,*)  "8", proj_x(:,1:nwan(isppol),isppol)
!  write(std_out,*)  "9",proj_zona(1:nwan(isppol),isppol)
!  write(std_out,*)  "10",exclude_bands(:,isppol)
!  end do!isppol
!  testdebug:  ovikp(1,1)=1
!  ^^^^^^^^^^^^^^^^^^^^^^^ end ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 end if  ! lwanniersetup
 do isppol=1,nsppol
   band_in(:,isppol)=.true.
   do icb=1,mband
     if(exclude_bands(icb,isppol).ne.0)  band_in(exclude_bands(icb,isppol),isppol)=.false.
   end do
   nband_inc(isppol)=0
   do icb=1, mband
     if (band_in(icb,isppol)) then
       nband_inc(isppol)=nband_inc(isppol)+1
     end if
   end do
 end do !isppol
 if(any(mband.gt.num_bands(:))) then
   write(message, '(a,a)' )ch10,&
&   '   Following bands are excluded from the calculation of wannier functions:'
   call wrtout(std_out,message,'COLL')

   do isppol=1,nsppol
     if(nsppol==2) then
       write(message,'("For spin",i4)')isppol
!      write(message,'(a,i)')'For spin=',isppol
       call wrtout(std_out,message,'COLL')
     end if !nsppol
     do jj=1,mband-num_bands(isppol),10
       write(message,'(10i7)') exclude_bands(jj:min(jj+9,mband-num_bands(isppol)),isppol)
       call wrtout(std_out,message,'COLL')
     end do
   end do !isppol
 end if

 do isppol=1,nsppol
   if(nsppol==2) then
     write(message,'("For spin",i4)')isppol
     call wrtout(std_out,message,'COLL')
   end if !nsppol
   write(message, '(a,i6,3a)' )ch10,&
&   nwan(isppol),' wannier functions will be computed (see ',trim(filew90_win(isppol)),')'
   call wrtout(std_out,message,'COLL')
!  write(std_out,*) exclude_bands(icb),band_in(icb)
!  ^^^^^^^^^^^^^^^END OF READING
   write(message, '(a,i6,a)' )ch10,&
&   num_bands(isppol),' bands will be used to extract wannier functions'
   call wrtout(std_out,message,'COLL')
   if(num_bands(isppol).lt.nwan(isppol)) then
     write(message, '(4a)' )&
&     ' number of bands is lower than the number of wannier functions',ch10,&
&     ' Action : check input file and ',trim(filew90_win(isppol))
     ABI_ERROR(message)
   else if (num_bands(isppol)==nwan(isppol)) then
     write(message, '(a,a,a,a)' )ch10,&
&     '   Number of bands is equal to the number of wannier functions',ch10,&
&     '   Disentanglement will not be necessary'
     call wrtout(std_out,message,'COLL')
   else if  (num_bands(isppol).gt.nwan(isppol)) then
     write(message, '(a,a,a,a)' )ch10,&
&     '   Number of bands is larger than the number of wannier functions',ch10,&
&     '   Disentanglement will be necessary'
     call wrtout(std_out,message,'COLL')
   end if
   write(message, '(2x,a,a,i3,1x,a)' )ch10,&
&   '   Each k-point has', nntot,'neighbours'
   call wrtout(std_out,message,'COLL')

 end do !isppol

end subroutine mlwfovlp_setup
!!***

!!****f* m_mlwfovlp/mlwfovlp_pw
!! NAME
!! mlwfovlp_pw
!!
!! FUNCTION
!! Routine which computes PW part of overlap M_{mn}(k,b)
!! for Wannier code (www.wannier.org f90 version).
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!  g1(3,nkpt,nntot) = G vector shift which is necessary to obtain k1+b
!!  iwav(mband,nkpt,nsppol): shift for pw components in cg.
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw.
!!  nfft=(effective) number of FFT grid points (for this processor) (see NOTES at beginning of scfcv)
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ovikp(nkpt,nntot)= gives  nntot value of k2 (in the BZ) for each k1  (k2=k1+b mod(G))
!!  seed_name= seed_name of files containing cg for all k-points to be used with MPI
!!
!! OUTPUT
!!  cm1(2,mband,mband,nntot,nkpt,nsppol): overlap <u_(nk1)|u_(mk1+b)>.

!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      m_mlwfovlp
!!
!! CHILDREN
!!      initylmr,matrginv,rotmat
!!
!! SOURCE

subroutine mlwfovlp_pw(cg,cm1,g1,iwav,kg,mband,mkmem,mpi_enreg,mpw,nfft,ngfft,nkpt,nntot,&
&  npwarr,nspinor,nsppol,ovikp,seed_name)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,mpw,nfft,nkpt,nntot
 integer,intent(in) :: nspinor,nsppol
 character(len=fnlen) ::  seed_name  !seed names of files containing cg info used in case of MPI
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: g1(3,nkpt,nntot),kg(3,mpw*mkmem),ngfft(18),npwarr(nkpt)
 integer,intent(in) :: iwav(mband,nkpt,nsppol)
 integer,intent(in) :: ovikp(nkpt,nntot)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(out) :: cm1(2,mband,mband,nntot,nkpt,nsppol)

!Local variables-------------------------------
!scalars
 integer :: iband1,iband2,ierr,ig,ig1,ig1b,ig2,ig2b
 integer :: ig3,ig3b,igk1,igk2,igks1,igks2,ii,ikg,ikpt,ikpt1,ikpt2,imntot,index,intot,ios,ipw
 integer :: ispinor,isppol,iunit,me,n1,n2,n3,npoint,npoint2,npw_k,npw_k2
 integer :: nprocs,spaceComm
 integer,allocatable::indpwk(:,:),kg_k(:,:)
 integer,allocatable :: invpwk(:,:)
 character(len=500) :: message
 character(len=fnlen) ::  cg_file  !file containing cg info used in case of MPI
 logical::lfile
 real(dp),allocatable :: cg_read(:,:) !to be used in case of MPI

!************************************************************************

 write(message, '(a,a)' ) ch10,&
& '** mlwfovlp_pw : compute pw part of overlap'
 call wrtout(std_out,  message,'COLL')

!initialize flags
 lfile=.false.
!mpi initialization
 spaceComm=MPI_enreg%comm_cell
 nprocs=xmpi_comm_size(spaceComm)
 me=MPI_enreg%me_kpt

 if(nprocs>1) then
   ABI_ALLOCATE(cg_read,(2,nspinor*mpw*mband))
 end if


!****************compute intermediate quantities  (index, shifts) ******
!------------compute index for g points--------------------------------
!ig is a plane waves which belongs to the sphere ecut for ikpt (they
!are npwarr(ikpt))
!npoint is the position in the grid of planes waves
!(they are nfft)
!indpwk is a application ig-> npoint
!invpwk is not an application (some npoint have no ig corresponding)
!cg are ordered with npw_k !
!----------------------------------------------------------------------
!------------compute index for g points--------------------------------
!----------------------------------------------------------------------
 write(message, '(a,a)' ) ch10,&
& '   first compute index for g-points'
 call wrtout(std_out,  message,'COLL')
!
!Allocations
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(indpwk,(nkpt,mpw))
 ABI_ALLOCATE(invpwk,(nkpt,nfft))
!
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 invpwk=0
 indpwk=0
 kg_k=0
 do isppol=1,1  !invpwk is not spin dependent
!  so we just do it once
   ikg=0
   do ikpt=1,nkpt
!
!    MPI:cycle over k-points not treated by this node
!
     if ( ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-me)  /=0) CYCLE

!
!    write(std_out,*)'me',me,'ikpt',ikpt,'isppol',isppol
     do npoint=1,nfft
       if(invpwk(ikpt,npoint)/=0 )then
         write(std_out,*) "error0 , invpwk is overwritten"
         write(std_out,*) ikpt,npoint
         ABI_ERROR("Aborting now")
       end if
     end do
     npw_k=npwarr(ikpt)
!    write(std_out,*) ikpt,npw_k,nfft
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     do ig=1,npw_k
       if(ig.gt.mpw) then
         write(std_out,*)"error ig",ig,"greater than mpw ",mpw
         ABI_ERROR("Aborting now")
       end if
       if(indpwk(ikpt,ig)/=0) then
         write(std_out,*) "error, indpwk is overwritten"
         write(std_out,*) ikpt,ig,indpwk(ikpt,ig)
         ABI_ERROR("Aborting now")
       end if
       ig1=modulo(kg_k(1,ig),n1)
       ig2=modulo(kg_k(2,ig),n2)
       ig3=modulo(kg_k(3,ig),n3)
       indpwk(ikpt,ig)=ig1+1+n1*(ig2+n2*ig3)
       npoint=indpwk(ikpt,ig)
       if(npoint.gt.nfft) then
         ABI_ERROR("error npoint")
       end if
!      write(std_out,*) ikpt,ig,npoint,invpwk(ikpt,npoint)
       if(invpwk(ikpt,npoint)/=0) then
         write(std_out,*) "error, invpwk is overwritten"
         write(std_out,*) ikpt,ig,npoint,invpwk(ikpt,npoint)
         ABI_ERROR("Aborting now")
       end if
       invpwk(ikpt,npoint)=ig
!      write(std_out,*)'ikpt,npoint,invpwk',ikpt,npoint,invpwk(ikpt,npoint)
!      if(ikpt.eq.1) write(std_out,*) "ig npoint",ig, npoint
!      write(std_out,*) "ikpt ig npoint",ikpt,ig, npoint
     end do
     ikg=ikg+npw_k

   end do !ikpt
 end do !isppol
!write(std_out,*) "index for g points has been computed"

 call xmpi_barrier(spaceComm)
 call xmpi_sum(invpwk,spaceComm,ierr)

!----------------------------------------------------------------------
!------------test invpwk-----------------------------------------------
!----------------------------------------------------------------------
!write(std_out,*) "TEST INVPWK"
!ikpt=3
!isppol=1
!do ig=1,npwarr(ikpt)
!npoint=indpwk(ikpt,ig)
!write(std_out,*) "ig npoint    ",ig, npoint
!write(std_out,*) "ig npoint inv",invpwk(ikpt,npoint),npoint
!end do
!do ig3=1,n3
!do ig2=1,n2
!do ig1=1,n1
!npoint=ig1+(ig2-1)*n1+(ig3-1)*n2*n1
!ig=invpwk(ikpt,npoint)
!!   if(ig/=0)  write(std_out,*) "ig npoint",ig, npoint
!end do
!end do
!end do

!Deallocate unused variables
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(indpwk)

!***********************************************************************
!**calculate overlap M_{mn}(k,b)=<\Psi_{k,m}|e^{-ibr}|\Psi_{k+b,n}>*****
!***********************************************************************
 write(message, '(a,a)' ) ch10,&
& '   mlwfovlp_pw : compute overlaps '
 call wrtout(std_out,  message,'COLL')
 write(message, '(a,a)' ) ch10,&
& "     nkpt  nntot  mband "
 call wrtout(std_out,  message,'COLL')
 write(message, '(i6,2x,i6,2x,i6,2x,i6)' ) &
& nkpt,nntot,mband
 call wrtout(std_out,  message,'COLL')
 cm1=zero
 write(message, '(a)' )  '  '
 call wrtout(std_out,  message,'COLL')
 do isppol=1,nsppol
   imntot=0
   do ikpt1=1,nkpt
!
!    MPI:cycle over k-points not treated by this node
!
     if ( ABS(MPI_enreg%proc_distrb(ikpt1,1,isppol)-me)  /=0) CYCLE
!
     write(message, '(a,i6,a,i6,a,i6)' ) &
&     '     Processor',me,' computes k-point',ikpt1,' and spin=',isppol
     call wrtout(std_out,  message,'COLL')
!    write(std_out,*)trim(message)

     do intot=1,nntot
       lfile=.false. !flag to know if this kpt will be read from a file, see below
       imntot=imntot+1
       ikpt2= ovikp(ikpt1,intot)
!      write(std_out,*)'me',me,'ikpt1',ikpt1,'ikpt2',ikpt2,'intot',intot,'isppol',isppol

!
!      MPI: if ikpt2 not found in this processor then
!      read info from an unformatted file
!
       if ( ABS(MPI_enreg%proc_distrb(ikpt2,1,isppol)-me)  /=0) then
         lfile=.true.
         write(cg_file,'(a,I5.5,".",I1)') trim(seed_name),ikpt2,isppol
         iunit=1000+ikpt2+ikpt2*(isppol-1)
         npw_k2=npwarr(ikpt2)
         open (unit=iunit, file=cg_file,form='unformatted',status='old',iostat=ios)
         if(ios /= 0) then
           write(message,*) " mlwfovlp_pw: file",trim(cg_file), "not found"
           ABI_ERROR(message)
         end if
!
         do iband2=1,mband
           do ipw=1,npw_k2*nspinor
             index=ipw+(iband2-1)*npw_k2*nspinor
             read(iunit) (cg_read(ii,index),ii=1,2)
!            if(me==0 .and. ikpt2==4)write(300,*)'ipw,iband2,index',ipw,iband2,index,cg_read(:,index)
!            if(me==1 .and. ikpt2==4)write(301,*)'ipw,iband2,index',ipw,iband2,index,cg_read(:,index)
           end do
         end do
         close(iunit)
       end if
!
       npw_k=npwarr(ikpt1)
       npw_k2=npwarr(ikpt2)
       do ig3=1,n3
         do ig2=1,n2
           do ig1=1,n1
!            write(std_out,*) isppol,ikpt1,iband1,iband2,intot
             npoint=ig1+(ig2-1)*n1+(ig3-1)*n2*n1
             if(npoint.gt.nfft) then
               write(std_out,*) "error npoint  b"
               ABI_ERROR("Aborting now")
             end if
             ig1b=ig1+g1(1,ikpt1,intot)
             ig2b=ig2+g1(2,ikpt1,intot)
             ig3b=ig3+g1(3,ikpt1,intot)
!            write(std_out,*) ig1,ig2,ig3
!            write(std_out,*) ig1b,ig2b,ig3b
             if(ig1b.lt.1) ig1b=ig1b+n1
             if(ig2b.lt.1) ig2b=ig2b+n2
             if(ig3b.lt.1) ig3b=ig3b+n3
             if(ig1b.gt.n1) ig1b=ig1b-n1
             if(ig2b.gt.n2) ig2b=ig2b-n2
             if(ig3b.gt.n3) ig3b=ig3b-n3
             npoint2=ig1b+(ig2b-1)*n1+(ig3b-1)*n2*n1
             if(npoint2.gt.nfft) then
               write(std_out,*)"error npoint  c"
               ABI_ERROR("Aborting now")
             end if
             igk1=invpwk(ikpt1,npoint)
             igk2=invpwk(ikpt2,npoint2)

!            if(intot==10) write(std_out,*)'Before igk1 and igk2',ikpt1,ikpt2,isppol

             if(igk1/=0.and.igk2/=0) then
               do iband2=1,mband
                 do iband1=1,mband
                   do ispinor=1,nspinor
!                    igks1= (igk1*nspinor)-(nspinor-ispinor)
!                    igks2= (igk2*nspinor)-(nspinor-ispinor)
                     igks1= igk1+ (ispinor-1)*npw_k
                     igks2= igk2+ (ispinor-1)*npw_k2

!                    Here the igks is to include the spinor component missing in igk
                     if(lfile) index=igks2+npw_k2*nspinor*(iband2-1) !In case of MPI, see below
!
!                    If MPI sometimes the info was read from an unformatted file
!                    If that is the case lfile==.true.
!
! TODO: this filter should be outside, not inside 1000 loops!!!
                     if(lfile) then
                       cm1(1,iband1,iband2,intot,ikpt1,isppol)=cm1(1,iband1,iband2,intot,ikpt1,isppol)+ &
&                       cg(1,igks1+iwav(iband1,ikpt1,isppol))*cg_read(1,index)&
&                       + cg(2,igks1+iwav(iband1,ikpt1,isppol))*cg_read(2,index)
                       cm1(2,iband1,iband2,intot,ikpt1,isppol)=cm1(2,iband1,iband2,intot,ikpt1,isppol)+ &
&                       cg(1,igks1+iwav(iband1,ikpt1,isppol))*cg_read(2,index)&
&                       - cg(2,igks1+iwav(iband1,ikpt1,isppol))*cg_read(1,index)
!
                     else
                       cm1(1,iband1,iband2,intot,ikpt1,isppol)=cm1(1,iband1,iband2,intot,ikpt1,isppol)+ &
&                       cg(1,igks1+iwav(iband1,ikpt1,isppol))*cg(1,igks2+iwav(iband2,ikpt2,isppol))&
&                       + cg(2,igks1+iwav(iband1,ikpt1,isppol))*cg(2,igks2+iwav(iband2,ikpt2,isppol))
                       cm1(2,iband1,iband2,intot,ikpt1,isppol)=cm1(2,iband1,iband2,intot,ikpt1,isppol)+ &
&                       cg(1,igks1+iwav(iband1,ikpt1,isppol))*cg(2,igks2+iwav(iband2,ikpt2,isppol))&
&                       - cg(2,igks1+iwav(iband1,ikpt1,isppol))*cg(1,igks2+iwav(iband2,ikpt2,isppol))
                     end if
                   end do !ispinor
                 end do ! iband1
               end do ! iband2
             end if
           end do ! ig1
         end do ! ig2
       end do ! ig3
     end do ! intot
   end do ! ikpt1
 end do ! isppol
!
!Deallocations
!
 ABI_DEALLOCATE(invpwk)
 if(nprocs>1)  then
   ABI_DEALLOCATE(cg_read)
 end if

 end subroutine mlwfovlp_pw
!!***

!!****f* m_mlwfovlp/mlwfovlp_proj
!! NAME
!! mlwfovlp_proj
!!
!! FUNCTION
!! Routine which computes projection A_{mn}(k) for Wannier code (www.wannier.org f90 version).
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  filew90_win = secondary input file for wannier90   (WAS NOT USED IN v6.7.1 - so has been temporarily removed)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  lproj= flag 0: no projections, 1: random projections,
!!              2: projections on atomic orbitals
!!              3: projections on projectors
!!  mband=maximum number of bands
!!  mkmem =number of k points treated by this node.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nkpt=number of k points.
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  num_bands=number of bands actually used to construct the wannier function
!!  nwan= number of wannier fonctions (read in wannier90.win).
!!  proj_l(mband)= angular part of the projection function (quantum number l)
!!  proj_m(mband)= angular part of the projection function (quantum number m)
!!  proj_radial(mband)= radial part of the projection.
!!  proj_site(3,mband)= site of the projection.
!!  proj_x(3,mband)= x axis for the projection.
!!  proj_z(3,mband)= z axis for the projection.
!!  proj_zona(mband)= extension of the radial part.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!! OUTPUT
!!  A_matrix(num_bands,nwan,nkpt,nsppol)= Matrix of projections needed by wannier_run
!!  ( also wannier90random.amn is written)
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      m_mlwfovlp
!!
!! CHILDREN
!!      initylmr,matrginv,rotmat
!!
!! SOURCE

 subroutine mlwfovlp_proj(A_matrix,band_in,cg,cprj,dtset,gprimd,just_augmentation,kg,&
&lproj,max_num_bands,mband,mkmem,mpi_enreg,mpw,mwan,natom,nattyp,&
&nkpt,npwarr,nspinor,&
&nsppol,ntypat,num_bands,nwan,pawtab,proj_l,proj_m,proj_radial,&
&proj_site,proj_x,proj_z,proj_zona,psps,ucvol)

!Arguments ------------------------------------
!scalars
 complex(dpc),parameter :: c1=(1._dp,0._dp)
 integer,intent(in) :: lproj,max_num_bands,mband,mkmem,mpw,mwan,natom,nkpt,nspinor,nsppol
 integer,intent(in) :: ntypat
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer ::nattyp(ntypat)
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt),num_bands(nsppol),nwan(nsppol),proj_l(mband,nsppol)
 integer,intent(in) :: proj_m(mband,nsppol)
 integer,intent(inout)::proj_radial(mband,nsppol)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: gprimd(3,3),proj_site(3,mband,nsppol)
 real(dp),intent(in) :: proj_x(3,mband,nsppol),proj_z(3,mband,nsppol),proj_zona(mband,nsppol)
 complex(dpc),intent(out) :: A_matrix(max_num_bands,mwan,nkpt,nsppol)
!character(len=fnlen),intent(in) :: filew90_win(nsppol)
 logical,intent(in) :: band_in(mband,nsppol)
 logical,intent(in)::just_augmentation(mwan,nsppol)
 type(pawcprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: iatom,iatprjn,iband,iband1,iband2,ibg,icat,icg,icg_shift
 integer :: idum,idx,ikg,ikpt,ilmn,ipw,iproj
 integer :: ispinor,isppol,itypat,iwan,jband,jj1,libprjn
 integer :: lmn_size,natprjn,nband_k,nbprjn,npw_k
 integer :: sumtmp
 integer :: max_lmax,max_lmax2,mproj,nprocs,spaceComm,rank
 real(dp),parameter :: qtol=2.0d-8
 real(dp) :: arg,norm_error,norm_error_bar
 real(dp) :: ucvol,x1,x2,xnorm,xnormb,xx,yy,zz
 complex(dpc) :: amn_tmp(nspinor)
 complex(dpc) :: cstr_fact
 character(len=500) :: message
!arrays
 integer :: kg_k(3,mpw),lmax(nsppol),lmax2(nsppol),nproj(nsppol)
 integer,allocatable :: lprjn(:),npprjn(:)
 real(dp) :: kpg(3),kpt(3)
 real(dp),allocatable :: amn(:,:,:,:,:),amn2(:,:,:,:,:,:,:)
 real(dp),allocatable :: gsum2(:),kpg2(:),radial(:)
 complex(dpc),allocatable :: gf(:,:),gft_lm(:)
 complex(dpc),allocatable :: ylmc_fac(:,:,:),ylmcp(:)

!no_abirules
!Tables 3.1 & 3.2, User guide
 integer,save :: orb_l_defs(-5:3)=(/2,2,1,1,1,0,1,2,3/)
! integer,parameter :: mtransfo(0:3,7)=&
!&  reshape((/1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,-2,-1,2,1,0,0,0,-1,1,2,-2,-3,3/),(/4,7/))

!************************************************************************

!mpi initialization
 spaceComm=MPI_enreg%comm_cell
 nprocs=xmpi_comm_size(spaceComm)
 rank=MPI_enreg%me_kpt

!Check input variables
 if ((lproj/=1).and.(lproj/=2).and.(lproj/=5)) then
   write(message, '(3a)' )&
&   ' Value of lproj no allowed ',ch10,&
&   ' Action : change lproj.'
   ABI_ERROR(message)
 end if

 write(message, '(a,a)' )ch10,&
& '** mlwfovlp_proj:  compute A_matrix of initial guess for wannier functions'
 call wrtout(std_out,message,'COLL')

!Initialize to 0.d0
 A_matrix(:,:,:,:)=cmplx(0.d0,0.d0)

!End of preliminarities

!********************* Write Random projectors
 if(lproj==1) then
   idum=123456
!  Compute random projections
   ABI_ALLOCATE(amn,(2,mband,mwan,nkpt,nsppol))
   amn=zero
   do isppol=1,nsppol
     do ikpt=1,nkpt
!
!      MPI: cycle over kpts not treated by this node
!
       if (ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-rank)/=0) CYCLE
!      write(std_out,'("kpt loop2: ikpt",i3," rank ",i3)') ikpt,rank

!
       do iband1=1,mband
         xnormb=0.d0
         do iband2=1,nwan(isppol)
           x1=uniformrandom(idum)
           x2=uniformrandom(idum)
           xnorm=sqrt(x1**2+x2**2)
           xnormb=xnormb+xnorm
           amn(1,iband1,iband2,ikpt,isppol)=x1
           amn(2,iband1,iband2,ikpt,isppol)=x2
         end do
         do iband2=1,nwan(isppol)
           amn(1,iband1,iband2,ikpt,isppol)=amn(1,iband1,iband2,ikpt,isppol)/xnormb
           amn(2,iband1,iband2,ikpt,isppol)=amn(2,iband1,iband2,ikpt,isppol)/xnormb
         end do !iband2
       end do !iband1
     end do !ikpt
   end do !isppol
   do isppol=1,nsppol
     do ikpt=1,nkpt
!
!      MPI: cycle over kpts not treated by this node
!
       if (ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-rank)/=0) CYCLE
!
       do iband2=1,nwan(isppol)
         jband=0
         do iband1=1,mband
           if(band_in(iband1,isppol)) then
             jband=jband+1
             if(jband.gt.num_bands(isppol)) then
               write(message, '(3a)' )&
&               'Value of jband is above num_bands ',ch10,&
&               'Action: contact Abinit group'
               ABI_ERROR(message)
             end if
             A_matrix(jband,iband2,ikpt,isppol)=cmplx(amn(1,iband1,iband2,ikpt,isppol),amn(2,iband1,iband2,ikpt,isppol))
           end if
         end do !iband1
       end do !iband2
     end do !ikpt
   end do !isppol
   ABI_DEALLOCATE(amn)
 end if

!********************* Projection on atomic orbitals based on .win file
 if( lproj==2) then !based on .win file
   nproj(:)=nwan(:)/nspinor !if spinors, then the number of projections are
   mproj=maxval(nproj(:))
!  half the total of wannier functions
!
!  obtain lmax and lmax2
   lmax(:)=0
   lmax2(:)=0
!
   do isppol=1,nsppol
     do iproj=1,nproj(isppol)
       lmax(isppol)=max(lmax(isppol),orb_l_defs(proj_l(iproj,isppol)))
     end do !iproj
     lmax2(isppol)=(lmax(isppol)+1)**2
   end do !isppol
   max_lmax=maxval(lmax(:))
   max_lmax2=maxval(lmax2(:))
!  Allocate arrays
   ABI_ALLOCATE(ylmc_fac,(max_lmax2,mproj,nsppol))
!
!  get ylmfac, factor used for rotations and hybrid orbitals
   do isppol=1,nsppol
     call mlwfovlp_ylmfac(ylmc_fac(1:lmax2(isppol),1:nproj(isppol),isppol),lmax(isppol),lmax2(isppol),&
&     mband,nproj(isppol),proj_l(:,isppol),proj_m(:,isppol),proj_x(:,:,isppol),proj_z(:,:,isppol))
   end do
!
   norm_error=zero
   norm_error_bar=zero
   icg=0
!
   do isppol=1,nsppol
!    Allocate arrays
!      this has to be done this way because the variable icg changes at the end of the
!      cycle. We cannot just skip the hole cycle.
     ABI_ALLOCATE(gf,(mpw,nproj(isppol)))
     ABI_ALLOCATE(gft_lm,(lmax2(isppol)))
     ABI_ALLOCATE(gsum2,(nproj(isppol)))
     ABI_ALLOCATE(kpg2,(mpw))
     ABI_ALLOCATE(radial,(lmax2(isppol)))
     ABI_ALLOCATE(ylmcp,(lmax2(isppol)))
!
     ikg=0
     do ikpt=1, nkpt
!
!      MPI: cycle over kpts not treated by this node
!
       if (ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-rank)/=0) CYCLE
!
       write(message, '(a,i6,a,2i6)' ) &
&       '   processor',rank,' will compute k-point,spin=',ikpt,isppol
       call wrtout(std_out,  message,'COLL')

!      Initialize variables
       npw_k=npwarr(ikpt)
       gsum2(:)=0.d0
       gf(:,:) = (0.d0,0.d0)
       kpt(:)=dtset%kpt(:,ikpt)
       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)

       do ipw=1, npw_k
         kpg(1)= (kpt(1) + real(kg_k(1,ipw),dp))     !k+G
         kpg(2)= (kpt(2) + real(kg_k(2,ipw),dp))
         kpg(3)= (kpt(3) + real(kg_k(3,ipw),dp))
!
!        Calculate modulus of k+G
         xx=gprimd(1,1)*kpg(1)+gprimd(1,2)*kpg(2)+gprimd(1,3)*kpg(3)
         yy=gprimd(2,1)*kpg(1)+gprimd(2,2)*kpg(2)+gprimd(2,3)*kpg(3)
         zz=gprimd(3,1)*kpg(1)+gprimd(3,2)*kpg(2)+gprimd(3,3)*kpg(3)
         kpg2(ipw)= two_pi*sqrt(xx**2+yy**2+zz**2)
!
!        Complex Y_lm for k+G
         if(lmax(isppol)==0) then
           ylmcp(1)=c1/sqrt(four_pi)
         else
           call ylm_cmplx(lmax(isppol),ylmcp,xx,yy,zz)
         end if
!
!        !
         do iproj=1,nproj(isppol)
!
!          In PAW, we can use proj_radial > 4 to indicate that we just
!          want the in-sphere contribution
!
           if( psps%usepaw==1) then
             if( just_augmentation(iproj,isppol)) cycle
           end if
!
!          obtain radial part
           call mlwfovlp_radial(proj_zona(iproj,isppol),lmax(isppol),lmax2(isppol)&
&           ,radial,proj_radial(iproj,isppol),kpg2(ipw))
!
!          scale complex representation of projector orbital with radial functions
!          of appropriate l
           gft_lm(:)=radial(:)*ylmc_fac(1:lmax2(isppol),iproj,isppol)
!
!          complex structure factor for projector orbital position
           arg = ( kpg(1)*proj_site(1,iproj,isppol) + &
&           kpg(2)*proj_site(2,iproj,isppol) + &
&           kpg(3)*proj_site(3,iproj,isppol) ) * 2*pi
           cstr_fact = cmplx(cos(arg), -sin(arg) )
!
!          obtain guiding functions
           gf(ipw,iproj)=cstr_fact*dot_product(ylmcp,gft_lm)
!
           gsum2(iproj)=gsum2(iproj)+real(gf(ipw,iproj))**2+aimag(gf(ipw,iproj))**2
         end do !iproj
       end do !ipw
!
       do iproj=1,nproj(isppol)
!
!        In PAW, we can use proj_radial > 4 to indicate that we just
!        want the in-sphere contribution
!
         if(psps%usepaw==1 ) then
           if (just_augmentation(iproj,isppol)) cycle
         end if
!
         gsum2(iproj)=16._dp*pi**2*gsum2(iproj)/ucvol
         gf(:,iproj)=gf(:,iproj)/sqrt(gsum2(iproj))
         norm_error=max(abs(gsum2(iproj)-one),norm_error)
         norm_error_bar=norm_error_bar+(gsum2(iproj)-one)**2
       end do !iproj
!
!      Guiding functions are computed.
!      compute overlaps of gaussian projectors and wave functions
       do iproj=1,nproj(isppol)
!
!        In PAW, we can use proj_radial > 4 to indicate that we just
!        want the in-sphere contribution
!
         if(psps%usepaw==1 ) then
           if ( just_augmentation(iproj,isppol)) cycle
         end if
!
         jband=0
         do iband=1,mband
           if(band_in(iband,isppol)) then
             icg_shift=npw_k*nspinor*(iband-1)+icg
             jband=jband+1
             amn_tmp(:)=cmplx(0.d0,0.d0)
             do ispinor=1,nspinor
               do ipw=1,npw_k
!
!                The case of spinors is tricky, we have nproj =  nwan/2
!                so we project to spin up and spin down separately, to have at
!                the end an amn matrix with nwan projections.
!
!                idx=ipw*nspinor - (nspinor-ispinor)
                 idx=ipw+(ispinor-1)*npw_k
                 amn_tmp(ispinor)=amn_tmp(ispinor)+gf(ipw,iproj)*cmplx(cg(1,idx+icg_shift),-cg(2,idx+icg_shift))
               end do !ipw
             end do !ispinor
             do ispinor=1,nspinor
               iwan=(iproj*nspinor)- (nspinor-ispinor)
               A_matrix(jband,iwan,ikpt,isppol)=amn_tmp(ispinor)
             end do
           end if !band_in
         end do !iband
       end do !iproj
       icg=icg+npw_k*nspinor*mband
       ikg=ikg+npw_k
     end do !ikpt
!    Deallocations
     ABI_DEALLOCATE(gf)
     ABI_DEALLOCATE(gft_lm)
     ABI_DEALLOCATE(gsum2)
     ABI_DEALLOCATE(kpg2)
     ABI_DEALLOCATE(radial)
     ABI_DEALLOCATE(ylmcp)
   end do !isppol
!
!  if(isppol==1) then
!    norm_error_bar=sqrt(norm_error_bar/real(nkpt*(nwan(1)),dp))
!  else
!    norm_error_bar=sqrt(norm_error_bar/real(nkpt*(nwan(1)+nwan(2)),dp))
!  end if
!  if(norm_error>0.05_dp) then
!  write(message, '(6a,f6.3,a,f6.3,12a)' )ch10,&
!  &     ' mlwfovlp_proj : WARNING',ch10,&
!  &     '  normalization error for wannier projectors',ch10,&
!  &     '  is',norm_error_bar,' (average) and',norm_error,' (max).',ch10,&
!  &     '  this may indicate more cell-to-cell overlap of the radial functions',ch10,&
!  &     '  than you want.',ch10,&
!  &     '  Action : modify zona (inverse range of radial functions)',ch10,&
!  '  under "begin projectors" in ',trim(filew90_win),' file',ch10
!  call wrtout(std_out,message,'COLL')
!  end if
!
!  !Deallocate
!  deallocate(ylmc_fac)
!
   ABI_DEALLOCATE(ylmc_fac)
 end if !lproj==2


!*************** computes projection  from PROJECTORS ********************
 if(lproj==3) then  !! if LPROJPRJ
!  ----- set values for projections --------------------- ! INPUT
!  nbprjn:number of  different l-values for projectors
!  lprjn: value of l for each projectors par ordre croissant
!  npprjn: number of projectors for each lprjn
   natprjn=1  ! atoms with wannier functions are first
   if(natprjn/=1) then ! in this case lprjn should depend on iatprjn
     ABI_ERROR("natprjn/=1")
   end if
   nbprjn=2
   ABI_ALLOCATE(lprjn,(nbprjn))
   lprjn(1)=0
   lprjn(2)=1
   ABI_ALLOCATE(npprjn,(0:lprjn(nbprjn)))
   npprjn(0)=1
   npprjn(1)=1
!  --- test coherence of nbprjn and nwan
   sumtmp=0
   do iatprjn=1,natprjn
     do libprjn=0,lprjn(nbprjn)
       sumtmp=sumtmp+(2*libprjn+1)*npprjn(libprjn)
     end do
   end do
   if(sumtmp/=nwan(1)) then
     write(std_out,*) "Number of Wannier orbitals is not equal to number of projections"
     write(std_out,*) "Action: check values of lprjn,npprjn % nwan"
     write(std_out,*) "nwan, sumtmp=",nwan,sumtmp
     ABI_ERROR("Aborting now")
   end if
!  --- end test of coherence
   ABI_ALLOCATE(amn2,(2,natom,nsppol,nkpt,mband,nspinor,nwan(1)))
   if(psps%usepaw==1) then
     amn2=zero
     ibg=0
     do isppol=1,nsppol
       do ikpt=1,nkpt
         nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
         do iband=1,nband_k
!          write(std_out,*)"amn2",iband,ibg,ikpt
           do ispinor=1,nspinor
             icat=1
             do itypat=1,dtset%ntypat
               lmn_size=pawtab(itypat)%lmn_size
               do iatom=icat,icat+nattyp(itypat)-1
                 jj1=0
                 do ilmn=1,lmn_size
                   if(iatom.le.natprjn) then
!                    do iwan=1,nwan
                     do libprjn=0,lprjn(nbprjn)
!                      if (psps%indlmn(1,ilmn,itypat)==proj_l(iwan)) then
!                      if (psps%indlmn(2,ilmn,itypat)==mtransfo(proj_l(iwan),proj_m(iwan))) then
                       if (psps%indlmn(1,ilmn,itypat)==libprjn) then
                         if (psps%indlmn(3,ilmn,itypat)<=npprjn(libprjn)) then
                           if(band_in(iband,isppol)) then
                             jj1=jj1+1
                             if(jj1>nwan(isppol)) then
                               write(std_out,*) "number of wannier orbitals is lower than lmn_size"
                               write(std_out,*) jj1,nwan(isppol)
                               ABI_ERROR("Aborting now")
                             end if
                             amn2(1,iatom,isppol,ikpt,iband,ispinor,jj1)=cprj(iatom,iband+ibg)%cp(1,ilmn)
                             amn2(2,iatom,isppol,ikpt,iband,ispinor,jj1)=cprj(iatom,iband+ibg)%cp(2,ilmn)
                           end if
                         end if
                       end if
                     end do ! libprjn
!                    endif
!                    endif
!                    enddo ! iwan
                   end if ! natprjn
                 end do !ilmn
               end do ! iatom
               icat=icat+nattyp(itypat)
             end do ! itypat
           end do ! ispinor
         end do !iband
         ibg=ibg+nband_k*nspinor
!        write(std_out,*)'amn2b',iband,ibg,ikpt
       end do !ikpt
     end do ! isppol

!    -----------------------  Save Amn   --------------------
     do isppol=1,nsppol
       do ikpt=1,nkpt
         do iband2=1,nwan(isppol)
           jband=0
           do iband1=1,mband
             if(band_in(iband1,isppol)) then
               jband=jband+1
               A_matrix(jband,iband2,ikpt,isppol)=&
&               cmplx(amn2(1,1,1,ikpt,iband1,1,iband2),amn2(2,1,1,ikpt,iband1,1,iband2))
             end if
           end do
         end do
       end do
     end do
   end if !usepaw
   ABI_DEALLOCATE(amn2)
   ABI_DEALLOCATE(npprjn)
   ABI_DEALLOCATE(lprjn)

 end if ! lproj==3

end subroutine mlwfovlp_proj
!!***

!!****f* m_mlwfovlp/mlwfovlp_projpaw
!! NAME
!! mlwfovlp_projpaw
!!
!! FUNCTION
!! Calculates the functions that are given to Wannier90 as an starting guess.
!! Here we project them inside the PAW spheres
!!
!! INPUTS
!!  band_in(mband)= logical array which indicates the bands to be excluded from the calculation
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  just_augmentation= flag used to indicate that we are just going
!!                     to compute augmentation part of the matrix
!!                     and we are excluding the plane wave part.
!!  mband= maximum number of bands
!!  mkmem= number of k points which can fit in memory; set to 0 if use disk
!!  natom= number of atoms in cell.
!!  nband(nkpt*nsppol)= array cointaining number of bands at each k-point and isppol
!!  nkpt=number of k points.
!!  num_bands=number of bands actually used to construct the wannier function (NOT USED IN 6.7.1 SO WAS TEMPORARILY REMOVED)
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  nwan= number of wannier fonctions (read in wannier90.win).
!!  pawrad(ntypat)= type(pawrad_type) radial information of paw objects
!!  pawtab(ntypat)= For PAW, TABulated data initialized at start
!!  proj_l(mband)= angular part of the projection function (quantum number l)
!!  proj_m(mband)= angular part of the projection function (quantum number m)
!!  proj_radial(mband)= radial part of the projection.
!!  proj_site(3,mband)= site of the projection.
!!  proj_x(3,mband)= x axis for the projection.
!!  proj_z(3,mband)= z axis for the projection.
!!  proj_zona(mband)= extension of the radial part.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)= Direct lattice vectors, Bohr units.
!!  typat(natom)= atom type
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  A_paw(max_num_bands,nwan,nkpt) = A matrix containing initial guess for MLWFs
!!                          (augmentation part of the matrix)
!!
!! SIDE EFFECTS
!!
!! NOTES
!! This routine is still under developement
!!
!! PARENTS
!!      m_mlwfovlp
!!
!! CHILDREN
!!      initylmr,matrginv,rotmat
!!
!! SOURCE

subroutine mlwfovlp_projpaw(A_paw,band_in,cprj,just_augmentation,max_num_bands,mband,mkmem,&
&mwan,natom,nband,nkpt,&
&nspinor,nsppol,ntypat,nwan,pawrad,pawtab,&
&proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,psps,&
&rprimd,typat,xred)

!Arguments ------------------------------------
 integer,intent(in) :: max_num_bands,mband,mkmem,mwan,natom,nkpt
 integer,intent(in) :: nspinor,nsppol,ntypat
 !arrays
 integer,intent(in) :: nband(nsppol*nkpt),nwan(nsppol)
 integer,intent(in) :: proj_l(mband,nsppol),proj_m(mband,nsppol),proj_radial(mband,nsppol)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in):: proj_site(3,mband,nsppol)
 real(dp),intent(in) :: proj_x(3,mband,nsppol),proj_z(3,mband,nsppol),proj_zona(mband,nsppol)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 complex(dpc),intent(out) :: A_paw(max_num_bands,mwan,nkpt,nsppol)
 logical,intent(in) :: band_in(mband,nsppol)
 logical,intent(in)::just_augmentation(mwan,nsppol)
 type(pawcprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)
 type(pseudopotential_type),intent(in) :: psps

!Local variables-------------------------------
 !local variables
 integer :: basis_size,iatom,iband,ii
 integer :: ikpt,ir,isppol,itypat,iwan,jband
 integer :: ll,lm,ln,mm,ilmn
 integer :: lmn_size,max_lmax2, mesh_size,nn
 integer :: lmax(nsppol),lmax2(nsppol)
 real(dp):: aa,int_rad2,prod_real,prod_imag
 real(dp),parameter :: dx=0.015d0,rmax=10.d0,xmin=0.d0
 real(dp):: sum,wan_lm_fac,x
 complex(dpc)::prod
 character(len=500) :: message
 !arrays
 integer :: index(mband,nkpt,nsppol)
 real(dp):: dist,norm(mwan,nsppol)
 real(dp):: proj_cart(3,mwan,nsppol),proj_site_unit(3,mwan,nsppol)
 real(dp):: xcart_unit(3,natom),xred_unit(3,natom)
 real(dp),allocatable ::aux(:),ff(:),r(:),int_rad(:),rad_int(:)
 real(dp),allocatable::ylmr_fac(:,:,:)

!no_abirules
!Tables 3.1 & 3.2, User guide
 integer,save :: orb_l_defs(-5:3)=(/2,2,1,1,1,0,1,2,3/)
!real(dp),allocatable :: ylm(:,:)


! *************************************************************************

 write(message, '(a,a)' )ch10,'** mlwfovlp_proj:  compute in-sphere part of A_matrix'
 call wrtout(std_out,message,'COLL')

!Check input variables
 do isppol=1,nsppol
   do iwan=1,nwan(nsppol)
     if(proj_radial(iwan,isppol)<1 .or. proj_radial(iwan,isppol)>4)then
       write(message,'(a,a,a,i6)')&
&       '  proj_radial should be between 1 and 4,',ch10,&
&       '  however, proj_radial=',proj_radial(iwan,isppol)
       ABI_BUG(message)
     end if
   end do
 end do

!Initialize
 A_paw(:,:,:,:)=cmplx(0.d0,0.d0)

!Get index for cprj
 ii=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     do iband=1,nband(ikpt)
       ii=ii+1
       index(iband,ikpt,isppol)=ii
     end do
   end do
 end do

!obtain lmax and lmax2
 lmax(:)=0
 lmax2(:)=0
 do isppol=1,nsppol
   do iwan=1,nwan(isppol)
     lmax(isppol)=max(lmax(isppol),orb_l_defs(proj_l(iwan,isppol)))
   end do !iwan
   lmax2(isppol)=(lmax(isppol)+1)**2
 end do
 max_lmax2=maxval(lmax2(:))
!
!get ylmfac, factor used for rotations and hybrid orbitals
!
 ABI_ALLOCATE(ylmr_fac,(max_lmax2,mwan,nsppol))


 do isppol=1,nsppol
   call mlwfovlp_ylmfar(ylmr_fac(1:lmax2(isppol),1:nwan(isppol),isppol),&
&   lmax(isppol),lmax2(isppol),mband,nwan(isppol),proj_l(:,isppol),proj_m(:,isppol),&
&   proj_x(:,:,isppol),proj_z(:,:,isppol))
!
!  Shift projection centers and atom centers to the primitive cell
!  This will be useful after, when we check if the Wannier function
!  lies on one specific atom
!
   proj_site_unit(:,:,:)=0.d0
   do iwan=1,nwan(isppol)
     do ii=1,3
       proj_site_unit(ii,iwan,isppol)=ABS(proj_site(ii,iwan,isppol)-AINT(proj_site(ii,iwan,isppol)) )
     end do
   end do
   do iatom=1,natom
     do ii=1,3
       xred_unit(ii,iatom)=ABS(xred(ii,iatom)-AINT(xred(ii,iatom)) )
     end do
   end do
   call xred2xcart(natom,rprimd,xcart_unit,xred_unit)
   call xred2xcart(mwan,rprimd,proj_cart(:,:,isppol),proj_site_unit(:,:,isppol))
!
!  Normalize the Wannier functions
!
!  Radial part
   mesh_size= nint((rmax - xmin ) / dx + 1)
   ABI_ALLOCATE( ff,(mesh_size))
   ABI_ALLOCATE(r,(mesh_size))
   ABI_ALLOCATE(rad_int,(mesh_size))
   ABI_ALLOCATE(aux,(mesh_size))
   do ir=1, mesh_size
     x=xmin+DBLE(ir-1)*dx
     r(ir)=x
   end do   !ir
   do iwan=1,nwan(isppol)
!    write(std_out,*)'iwan',iwan
!    radial functions shown in table 3.3 of wannier90 manual
     if(proj_radial(iwan,isppol)==1) ff(:) = 2.d0 * proj_zona(iwan,isppol)**(1.5d0) * exp(-proj_zona(iwan,isppol)*r(:))
     if(proj_radial(iwan,isppol)==2) ff(:) = 1.d0/(2.d0*sqrt(2.d0))*proj_zona(iwan,isppol)**(1.5d0) *&
&     (2.d0 - proj_zona(iwan,isppol)*r(:))*exp(-proj_zona(iwan,isppol)*r(:)/2.d0)
     if(proj_radial(iwan,isppol)==3) ff(:) = sqrt(4.d0/27.d0)*proj_zona(iwan,isppol)**(1.5d0)&
&     * (1.d0 - 2.d0*proj_zona(iwan,isppol)*r(:)/3.d0 + 2.d0*proj_zona(iwan,isppol)**2*r(:)**2/27.d0)&
&     * exp(-proj_zona(iwan,isppol) * r(:)/3.d0)

     if(proj_radial(iwan,isppol)/=4) then
       aux(:)=ff(:)**2*r(:)**2
       call simpson_int(mesh_size,dx,aux,rad_int)
       sum=0.d0
       do ir=1,mesh_size
         sum=sum+rad_int(ir)
       end do
       int_rad2=sum/real(mesh_size,dp)
!
!      do ir=1,mesh_size
!      if(iwan==1) write(400,*)r(ir),aux(ir),rad_int(ir)
!      end do
     else
!
!      ==4: gaussian function
!      f(x)=\exp(-1/4(x/aa)**2)
!      \int f(x)f(x) dx = \int \exp(-1/2(x/aa)**2) = aa*sqrt(2pi)
!
       int_rad2=sqrt(2.d0*pi)*proj_zona(iwan,isppol)
     end if

!
!    Now angular part
!
     prod_real=0.d0
     do lm=1,lmax2(isppol)
       wan_lm_fac=ylmr_fac(lm,iwan,isppol)
!      write(std_out,*)'wan_lm_fac',wan_lm_fac
!      write(std_out,*)'int_rad2',int_rad2
       prod_real= prod_real + wan_lm_fac**2 * int_rad2
     end do
     norm(iwan,isppol)=sqrt(prod_real)
   end do !iwan
   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(r)
   ABI_DEALLOCATE(rad_int)
   ABI_DEALLOCATE(aux)
!
!  Now that we found our guiding functions
!  We proceed with the internal product of
!  our guiding functions and the wave function
!  Amn=<G_m|\Psi_n> inside the sphere.
!  The term <G_m|\Psi_n> inside the sphere is:
!  = \sum_i <G_n | \phi_i - \tphi_i> <p_im|\Psi_m>
!
!
!  G_n \phi and \tphi can be decomposed in
!  a radial function times an angular function.
!
!
!  Big loop on iwan and iatom
!
   do iwan=1,nwan(isppol)
     do iatom=1,natom
!
!      check if center of wannier function coincides
!      with the center of the atom
!
       dist=((proj_cart(1,iwan,isppol)-xcart_unit(1,iatom))**2 + &
&       (proj_cart(2,iwan,isppol)-xcart_unit(2,iatom))**2 + &
&       (proj_cart(3,iwan,isppol)-xcart_unit(3,iatom))**2)**0.5
!
!      if the distance between the centers is major than 0.1 angstroms skip
!
       if( dist > 0.188972613) cycle
!
       write(message, '(2a,i4,a,i4,2a)')ch10, '   Wannier function center',iwan,' is on top of atom',&
&       iatom,ch10,'      Calculating in-sphere contribution'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
!
!      Get useful quantities
!
       itypat=typat(iatom)
       lmn_size=pawtab(itypat)%lmn_size
       basis_size=pawtab(itypat)%basis_size
       mesh_size=pawtab(itypat)%mesh_size
       ABI_ALLOCATE(int_rad,(basis_size))
       ABI_ALLOCATE(ff,(mesh_size))
       ABI_ALLOCATE(aux,(mesh_size))

!
!      Integrate first the radial part
!      and save it into an array
!
!
!      radial functions shown in table 3.3 of wannier90 manual
!
       if(proj_radial(iwan,isppol)==1) aux(1:mesh_size) = 2.d0 * proj_zona(iwan,isppol)**(1.5d0) *&
&       exp(-proj_zona(iwan,isppol)*pawrad(itypat)%rad(1:mesh_size))
       if(proj_radial(iwan,isppol)==2) aux(1:mesh_size) = 1.d0/(2.d0*sqrt(2.d0))*proj_zona(iwan,isppol)**(1.5d0) *&
&       (2.d0 - proj_zona(iwan,isppol)*pawrad(itypat)%rad(1:mesh_size)) &
&       * exp(-proj_zona(iwan,isppol)*pawrad(itypat)%rad(1:mesh_size)/2.d0)
       if(proj_radial(iwan,isppol)==3) aux(1:mesh_size) = sqrt(4.d0/27.d0)*proj_zona(iwan,isppol)**(1.5d0)&
&       * (1.d0 - 2.d0*proj_zona(iwan,isppol)*pawrad(itypat)%rad(1:mesh_size)/3.d0 &
&       + 2.d0*proj_zona(iwan,isppol)**2 *pawrad(itypat)%rad(1:mesh_size)**2/27.d0)&
&       * exp(-proj_zona(iwan,isppol) * pawrad(itypat)%rad(1:mesh_size)/3.d0)
!
!      ==4: gaussian function
!      f(x)=\exp(-1/4(x/aa)**2)
!
       if(proj_radial(iwan,isppol)==4) then
         aa=1.d0/proj_zona(iwan,isppol)
         aux(1:mesh_size)= exp(-0.25d0*(pawrad(itypat)%rad(1:mesh_size)*aa)**2)
       end if
!
!      Normalize aux
       aux(:)=aux(:)/norm(iwan,isppol)
!
       do ln=1,basis_size
         if(just_augmentation(iwan,isppol)) then
!
!          just augmentation region contribution
!          In this case there is no need to use \tphi
!          ff= \int R_wan(r) (R_phi(ln;r)/r ) r^2 dr
!
           ff(1:mesh_size)= aux(1:mesh_size) * pawtab(itypat)%phi(1:mesh_size,ln) &
&           * pawrad(itypat)%rad(1:mesh_size)
         else
!          Inside sphere contribution = \phi - \tphi
!          ff= \int R_wan(r) (R_phi(ln;r)/r - R_tphi(ln;r)/r) r^2 dr
           ff(1:mesh_size)= aux(1:mesh_size) * (pawtab(itypat)%phi(1:mesh_size,ln)-pawtab(itypat)%tphi(1:mesh_size,ln)) &
&           * pawrad(itypat)%rad(1:mesh_size)
         end if
!
!        Integration with simpson routine
!
         call simp_gen(int_rad(ln),ff,pawrad(itypat))
!        do ii=1,mesh_size
!        unit_ln=400+ln
!        if( iwan==1 ) write(unit_ln,*)pawrad(itypat)%rad(ii),ff(ii),int_rad(ln)
!        end do
       end do !ln
       ABI_DEALLOCATE(ff)
       ABI_DEALLOCATE(aux)
!
!      Now integrate the angular part
!      Cycle on i indices
!
!      prod_real=0.d0
       do ilmn=1, lmn_size
         ll=Psps%indlmn(1,ilmn,itypat)
         mm=Psps%indlmn(2,ilmn,itypat)
         nn=Psps%indlmn(3,ilmn,itypat)
         lm=Psps%indlmn(4,ilmn,itypat)
         ln=Psps%indlmn(5,ilmn,itypat)
!        write(std_out,*)'ll ',ll,' mm ',mm,'nn',nn,"lm",lm,"ln",ln
!
!        Get wannier factor for that lm component
         if(lm <=lmax2(isppol)) then
           wan_lm_fac=ylmr_fac(lm,iwan,isppol)
!          Make delta product
!          Here we integrate the angular part
!          Since the integral of the product of two spherical harmonics
!          is a delta function
           if( abs(wan_lm_fac) > 0.0d0) then
!            write(std_out,*) 'll',ll,'mm',mm,'lm',lm,'ln',ln,'factor',wan_lm_fac !lm index for wannier function
!
!            Calculate Amn_paw, now that the radial and angular integrations are done
!
             prod=cmplx(0.d0,0.d0)
             do ikpt=1,nkpt
               jband=0
               do iband=1,nband(ikpt)
                 if(band_in(iband,isppol)) then
                   jband=jband+1

                   prod_real= cprj(iatom,index(iband,ikpt,isppol))%cp(1,ilmn) * int_rad(ln) * wan_lm_fac
                   prod_imag= cprj(iatom,index(iband,ikpt,isppol))%cp(2,ilmn) * int_rad(ln) * wan_lm_fac
                   prod=cmplx(prod_real,prod_imag)

                   A_paw(jband,iwan,ikpt,isppol)=A_paw(jband,iwan,ikpt,isppol)+prod
                 end if !band_in
               end do !iband
             end do !ikpt
!
           end if !lm<=lmax2
         end if  ! abs(wan_lm_fac) > 0.0d0
       end do !ilmn=1, lmn_size
       ABI_DEALLOCATE(int_rad)
     end do !iatom
   end do !iwan
 end do !isppol

!Deallocate quantities
 ABI_DEALLOCATE(ylmr_fac)

end subroutine mlwfovlp_projpaw
!!***

!!****f* m_mlwfovlp/mlwfovlp_radial
!! NAME
!! mlwfovlp_radial
!!
!! FUNCTION
!! Calculates the radial part of the initial functions given to Wannier90
!! as an starting point for the minimization.
!! The trial functions are a set of solutions to the radial part of the hydrogenic
!! Schrodinger equation as it is explained in Table 3.3 of the Wannier90 user guide.
!!
!! INPUTS
!!  alpha= Z/a = zona
!!  lmax= maximum value of l
!!  rvalue= integer defining the choice for radial functions R(r).
!!   It can take values from 1-3.
!!   It is associted to the radial part of the hydrogenic Schrodinger equation for l=0,
!!   See the manual of Wannier90 for more information. (www.wannier.org)
!!  xx= scalar number used to calculate the spherical bessel function. J_il(xx)
!!
!! OUTPUT
!!  mlwfovlp_radial= radial part for initial projections used to construct MLWF
!!
!! SIDE EFFECTS
!!  None
!!
!! NOTES
!!  Calculates the radial part of the initial functions given as an initial
!!  guess by the user to construct the MLWF.
!!
!! PARENTS
!!      m_mlwfovlp
!!
!! CHILDREN
!!      initylmr,matrginv,rotmat
!!
!! SOURCE

subroutine mlwfovlp_radial(alpha,lmax,lmax2,radial,rvalue,xx)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmax,lmax2,rvalue
 real(dp),intent(in) :: alpha,xx
!arrays
 real(dp),intent(out) :: radial(lmax2)

!Local variables
!scalars
 integer :: ir,ll,lm,mesh,mm
 real(dp),parameter :: dx=0.015d0,rmax=10.d0,xmin=0.d0
 real(dp) :: aa,ftmp,gauss,rtmp,x
 character(len=500) :: message
!arrays
 real(dp),save :: dblefact(4)=(/1_dp,3_dp,15_dp,105_dp/)
 real(dp),allocatable :: aux(:),bes(:),cosr(:),func_r(:),r(:),rad_int(:)
 real(dp),allocatable :: sinr(:)

! *************************************************************************

!Radial functions in the form of hydrogenic orbitals as defined in the
!wannier90 manual.
 if(( rvalue > 0 ).and.(rvalue < 4)) then

!  mesh
   mesh= nint((rmax - xmin ) / dx + 1)
   ABI_ALLOCATE( bes,(mesh))
   ABI_ALLOCATE(func_r,(mesh))
   ABI_ALLOCATE(r,(mesh))
   ABI_ALLOCATE(rad_int,(mesh))
   ABI_ALLOCATE( aux,(mesh))
   ABI_ALLOCATE(cosr,(mesh))
   ABI_ALLOCATE(sinr,(mesh))
   do ir=1, mesh
     x=xmin+DBLE(ir-1)*dx
     r(ir)=x
   end do   !ir

!  radial functions shown in table 3.3 of wannier90 manual
   if (rvalue==1) func_r(:) = 2.d0 * alpha**(3.d0/2.d0) * exp(-alpha*r(:))
   if (rvalue==2) func_r(:) = 1.d0/(2.d0*sqrt(2.d0))*alpha**(3.d0/2.d0) *&
&   (2.d0 - alpha*r(:))*exp(-alpha*r(:)/2.d0)
   if (rvalue==3) func_r(:) = sqrt(4.d0/27.d0)*alpha**(3.d0/2.d0)&
&   * (1.d0 - 2.d0*alpha*r(:)/3.d0 + 2.d0*alpha**2*r(:)**2/27.d0)&
&   * exp(-alpha * r(:)/3.d0)

!  compute spherical bessel functions
   cosr(:)=cos(xx*r(:))
   sinr(:)=sin(xx*r(:))
   lm=0
   do ll=0,lmax
     call besjm(xx,bes,cosr,ll,mesh,sinr,r)
     aux(:)=bes(:)*func_r(:)*r(:)
!    do ir=1,mesh
!    write(310,*) r(ir),bes(ir)
!    end do
     call simpson_int(mesh,dx,aux,rad_int)
     rtmp=rad_int(mesh)/mesh
     do mm=-ll,ll
       lm=lm+1
       radial(lm)=rtmp
     end do !mm
   end do !ll
   ABI_DEALLOCATE(bes)
   ABI_DEALLOCATE(func_r)
   ABI_DEALLOCATE(r)
   ABI_DEALLOCATE(aux)
   ABI_DEALLOCATE(rad_int)
   ABI_DEALLOCATE(cosr)
   ABI_DEALLOCATE(sinr)

!  Radial part in the form of Gaussian functions of a given width
!  Taken by code of made by drh.
 elseif ( rvalue == 4) then
   aa=1._dp/alpha
   gauss=exp(-0.25_dp*(aa*xx)**2)
   lm=0
   do ll=0,lmax
     ftmp=(0.5_dp*pi)**(0.25_dp)*aa*sqrt(aa/dblefact(ll+1))*(aa*xx)**ll*gauss
     do mm=-ll,ll
       lm=lm+1
       radial(lm)=ftmp
     end do
   end do
 else ! rvalue < 0 of rvalue > 4
   write(message,'(a,i6,5a)')&
&   '  Radial function r=',rvalue,ch10,&
&   '  is not defined',ch10,&
&   '  Modify .win file',ch10
   ABI_BUG(message)
 end if !rvalue

end subroutine mlwfovlp_radial
!!***

!!****f* m_mlwfovlp/mlwfovlp_ylmfac
!! NAME
!! mlwfovlp_ylmfac
!!
!! FUNCTION
!! Routine that produces a factor by which the initial
!! guess of functions will be multiplied for the Wannier90 interface.
!! It is just used if there are rotations, or if the functions required
!! are linear combinations of the ylm real functions.
!!
!! Example,
!! For a function G(r)= 1/2 s + 1/3 px - 1/2 pz
!!   it would produce a matrix of the following form:
!!   [1/2,-1/2,1/3,0,0...0]
!!
!! The real spherical harmonics are given as factors of complex spherical harmonics
!! The real spherical harmonics are given in table 3.1 of Wannier90 user guide.
!!
!! INPUTS
!!  lmax= maximum l value for spherical harmonics
!!  lmax2=number of ylm functions
!!  mband=maximum number of bands
!!  nwan = number of wannier functions
!!  proj_l(mband)= angular part of the projection function (quantum number l)
!!  proj_m(mband)= angular part of the projection function (quantum number m)
!!  proj_x(3,mband)= x axis for the projection.
!!  proj_z(3,mband)= z axis for the projection.
!!
!! OUTPUT
!!  ylmc_fac(lmax2,nwan)=matrix containig a factor for ylm hybrid orbitals
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      m_mlwfovlp
!!
!! CHILDREN
!!      initylmr,matrginv,rotmat
!!
!! SOURCE


subroutine mlwfovlp_ylmfac(ylmc_fac,lmax,lmax2,mband,nwan,proj_l,proj_m,proj_x,proj_z)

!Arguments ------------------------------------
 integer, intent(in):: lmax,lmax2,nwan,mband
! arrays
 integer,intent(in) :: proj_l(mband),proj_m(mband)
 real(dp),intent(in) :: proj_x(3,mband),proj_z(3,mband)
 complex(dp),intent(out)::ylmc_fac(lmax2,nwan)
!
!Local variables-------------------------------
!
 integer :: orb_idx(16)=(/1,3,4,2,7,8,6,9,5,13,14,12,15,11,16,10/) !Tab3.1 Wannier90 user guide
 integer :: idum,ii,info,inversion_flag
 integer :: ir,iwan,jj,ll,lm,lmc,mm,mr
 real(dp):: onem,test
! arrays
 integer:: ipiv(lmax2)
 real(dp)::r(3,lmax2),rp(3,lmax2)
 real(dp)::rs2,rs3,rs6,rs12,umat(3,3)
 complex(dp)::crot(lmax2,lmax2),ctor(lmax2,lmax2),orb_lm(lmax2,-5:3,7)
 complex(dp):: ylmcp(lmax2)
 complex(dp):: ylmc_rr(lmax2,lmax2),ylmc_rr_save(lmax2,lmax2)
 complex(dp):: ylmc_rrinv(lmax2,lmax2),ylmc_rp(lmax2,lmax2)
 complex(dp),parameter :: c0=(0._dp,0._dp),c1=(1._dp,0._dp),ci=(0._dp,1._dp)
 character(len=500) :: message                   ! to be uncommented, if needed

! *************************************************************************


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!DEBUG
!write(std_out,*)'lmax ',lmax,'lmax2 ',lmax2
!write(std_out,*)'mband ',mband,'nwan ',nwan
!
!do iwan=1,nwan
!write(std_out,*)'iwan,proj_l, proj_m',proj_l(iwan),proj_m(iwan)
!write(std_out,*)'iwan,proj_x, proj_z',iwan,proj_x(:,iwan),proj_z(:,iwan)
!end do
!!END DEBUG

!constants for linear combinations of ylm's
 rs2=1._dp/sqrt(2._dp)
 rs3=1._dp/sqrt(3._dp)
 rs6=1._dp/sqrt(6._dp)
 rs12=1._dp/sqrt(12._dp)

!complex lm coefficients for real spherical harmonics in conventional order
!s, py,pz,px, dxy,dyz,dz2,dxz,dx2-y2, fy(3x2-y2),fxyz,fyz2,fz3,fxz2,
!fz(x2-y2),fx(x2-3y2)
 ctor(:,:)=c0
 do ll=0,lmax
   mm=0
   lm= ll**2+ll+mm+1
   ctor(lm,lm)=c1
   if(ll>0) then
     onem=one
     do mm=1,ll
       onem=-onem !(-1^mm)
       lm= ll**2+ll+mm+1
       lmc=ll**2+ll-mm+1
       ctor(lm ,lm )=rs2*c1
       ctor(lmc,lm )=onem*rs2*c1
       ctor(lm ,lmc)=rs2*ci
       ctor(lmc,lmc)=-onem*rs2*ci
     end do
   end if
 end do

 lm=0
 do ll=0,lmax
   do mm=-ll,ll
     lm=lm+1
     ctor(:,lm)=ctor(:,lm)*conjg(ci)**ll
   end do !mm
 end do !ll


!coefficients for basic wannier orbitals in Table 3.1 order
 orb_lm(:,:,:)=c0
 ii=0
 do ll=0,lmax
   do mr=1,2*ll+1
     ii=ii+1
     orb_lm(:,ll,mr)=ctor(:,orb_idx(ii))
   end do
 end do



!coefficients for linear combinations in table 3.2 order
 if(lmax>=1) then
!  s            px
   orb_lm(:,-1,1)=rs2*ctor(:,1)+rs2*ctor(:,4)
   orb_lm(:,-1,2)=rs2*ctor(:,1)-rs2*ctor(:,4)
!  s            px            py
   orb_lm(:,-2,1)=rs3*ctor(:,1)-rs6*ctor(:,4)+rs2*ctor(:,2)
   orb_lm(:,-2,2)=rs3*ctor(:,1)-rs6*ctor(:,4)-rs2*ctor(:,2)
   orb_lm(:,-2,3)=rs3*ctor(:,1)+2._dp*rs6*ctor(:,4)
!  s        px        py        pz
   orb_lm(:,-3,1)=half*(ctor(:,1)+ctor(:,4)+ctor(:,2)+ctor(:,3))
   orb_lm(:,-3,2)=half*(ctor(:,1)+ctor(:,4)-ctor(:,2)-ctor(:,3))
   orb_lm(:,-3,3)=half*(ctor(:,1)-ctor(:,4)+ctor(:,2)-ctor(:,3))
   orb_lm(:,-3,4)=half*(ctor(:,1)-ctor(:,4)-ctor(:,2)+ctor(:,3))
 end if
 if(lmax>=2) then
!  s            px            py
   orb_lm(:,-4,1)=rs3*ctor(:,1)-rs6*ctor(:,4)+rs2*ctor(:,2)
   orb_lm(:,-4,2)=rs3*ctor(:,1)-rs6*ctor(:,4)-rs2*ctor(:,2)
   orb_lm(:,-4,3)=rs3*ctor(:,1)+2._dp*rs6*ctor(:,4)
!  pz           dz2
   orb_lm(:,-4,4)= rs2*ctor(:,3)+rs2*ctor(:,7)
   orb_lm(:,-4,5)=-rs2*ctor(:,3)+rs2*ctor(:,7)
!  s            px            dz2         dx2-y2
   orb_lm(:,-5,1)=rs6*ctor(:,1)-rs2*ctor(:,4)-rs12*ctor(:,7)+half*ctor(:,9)
   orb_lm(:,-5,2)=rs6*ctor(:,1)+rs2*ctor(:,4)-rs12*ctor(:,7)+half*ctor(:,9)
!  s            py            dz2         dx2-y2
   orb_lm(:,-5,3)=rs6*ctor(:,1)-rs2*ctor(:,2)-rs12*ctor(:,7)-half*ctor(:,9)
   orb_lm(:,-5,4)=rs6*ctor(:,1)+rs2*ctor(:,2)-rs12*ctor(:,7)-half*ctor(:,9)
!  s            pz           dz2
   orb_lm(:,-5,5)=rs6*ctor(:,1)-rs2*ctor(:,3)+rs3*ctor(:,7)
   orb_lm(:,-5,6)=rs6*ctor(:,1)+rs2*ctor(:,3)+rs3*ctor(:,7)
 end if

!stuff complex wannier orbital coefficient array
 do iwan=1,nwan
   ylmc_fac(:,iwan)=orb_lm(:,proj_l(iwan),proj_m(iwan))
 end do


!setup to rotate ylmc_fac to new axes if called for
!skip if only s projectors are used
 if ( lmax>0 ) then
!  generate a set of nr=lmax2 random vectors
!  idum=123456
   do ir=1,lmax2
     do ii=1,3
       r(ii,ir) = uniformrandom(idum)-0.5d0
     end do !ii
     call ylm_cmplx(lmax,ylmcp,r(1,ir),r(2,ir),r(3,ir))
     ylmc_rr(ir,:)=conjg(ylmcp(:))
     ylmc_rr_save(ir,:)=conjg(ylmcp(:))
   end do !ir

   ylmc_rrinv(:,:)=c0
   do ii=1,lmax2
     ylmc_rrinv(ii,ii)=c1
   end do !ii
!  calculate inverse of ylmc(ir,lm) matrix
   call ZGESV(lmax2,lmax2,ylmc_rr,lmax2,ipiv,ylmc_rrinv,lmax2,info)

!  check that r points are independent (ie., that matrix inversion wasn't
!  too close to singular)
   ylmc_rr=matmul(ylmc_rrinv,ylmc_rr_save)
   test=zero
   do ii=1,lmax2
     ylmc_rr(ii,ii)=ylmc_rr(ii,ii)-c1
     do jj=1,lmax2
       test=max(abs(ylmc_rr(ii,jj)),test)
     end do !ii
   end do !jj
   if(test>tol8) then
     write(message, '(5a)' )&
&     '  matrix inversion error for wannier rotations',ch10,&
&     '  random vectors r(j,1:nr) are not all independent !! ',ch10,&
&     '  Action : re-seed uniformrandom or maybe just try again'
     ABI_ERROR(message)
   end if !test>tol8

!  end of the preliminaries, now to the rotations of the wannier orbitals
   do iwan=1,nwan
!    don't bother for s orbitals
     if(proj_l(iwan)==0) cycle
!    check for default axes and cycle if found
     if(proj_z(1,iwan)==zero .and. proj_z(2,iwan)==zero .and.&
&     proj_z(3,iwan)== one .and. proj_x(1,iwan)==one .and.&
&     proj_x(2,iwan)==zero .and. proj_x(3,iwan)==zero) cycle

!    get the u matrix that rotates the reference frame
     call rotmat(proj_x(:,iwan),proj_z(:,iwan),inversion_flag,umat)

!    find rotated r-vectors. Optional inversion
!    operation is an extension of the wannier90 axis-setting options
!    which only allow for proper axis rotations
     if(inversion_flag==1) then
       rp(:,:)= -matmul ( umat(:,:),  r(:,:) )
     else
       rp(:,:) = matmul ( umat(:,:) , r(:,:) )
     end if !inversion_flag

     do ir=1,lmax2
!      get the ylm representation of the rotated vectors
       call ylm_cmplx(lmax,ylmcp,rp(1,ir),rp(2,ir),rp(3,ir))
       ylmc_rp(ir,:)=conjg(ylmcp(:))
     end do !ir
!    the matrix product sum(ir) ylmc_rrinv(lm,ir)*ylmc_rp(ir,lm') gives the
!    the complex lmXlm matrix representation of the coordinate rotation
     crot(:,:)=matmul(ylmc_rrinv(:,:),ylmc_rp(:,:))

!    now rotate the current wannier orbital
     ylmcp(:)=matmul(crot(:,:),ylmc_fac(:,iwan))
     ylmc_fac(:,iwan)=ylmcp(:)

!    write(std_out,*)'ylmc_fac',ylmc_fac(:,iwan)
   end do !iwan
 end if !lmax>0

end subroutine mlwfovlp_ylmfac
!!***

!!****f* m_mlwfovlp/mlwfovlp_ylmfar
!! NAME
!! mlwfovlp_ylmfar
!!
!! FUNCTION
!! Routine that produces a fator by which the initial
!! guess of functions will be multiplied for the Wannier90 interface.
!! It is just used if there are rotations, or if the functions required
!! are linear combinations of the ylm real functions.
!!
!! Example,
!! For a function G(r)= 1/2 s + 1/3 px - 1/2 pz
!!   it would produce a matrix of the following form:
!!   [1/2,-1/2,1/3,0,0...0]
!!
!! This function is similar to mlwfovlp_ylmfac, but the factors it uses
!! real spherical harmonics instead of complex
!! spherical harmonics. Remember that real spherical harmonics
!! are linear combinations of complex
!! spherical harmonics
!!
!! INPUTS
!!  lmax= maximum l value for spherical harmonics
!!  lmax2=number of ylm functions
!!  mband=maximum number of bands
!!  nwan = number of wannier functions
!!  proj_l(mband)= angular part of the projection function (quantum number l)
!!  proj_m(mband)= angular part of the projection function (quantum number m)
!!  proj_x(3,mband)= x axis for the projection.
!!  proj_z(3,mband)= z axis for the projection.
!!
!! OUTPUT
!!  ylmc_fac(lmax2,nwan)=matrix containig a factor for ylm hybrid orbitals
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_mlwfovlp
!!
!! CHILDREN
!!      initylmr,matrginv,rotmat
!!
!! SOURCE

subroutine mlwfovlp_ylmfar(ylmr_fac,lmax,lmax2,mband,nwan,proj_l,proj_m,proj_x,proj_z)

!Arguments ------------------------------------
 integer, intent(in):: lmax,lmax2,nwan,mband
! arrays
 integer,intent(in) :: proj_l(mband),proj_m(mband)
 real(dp),intent(in) :: proj_x(3,mband),proj_z(3,mband)
 real(dp),intent(out)::ylmr_fac(lmax2,nwan)
!
!Local variables-------------------------------
!
 integer :: idum,ii,inversion_flag
 integer :: ir,iwan,jj,ll,lm,mm,mr
 real(dp):: onem,test
! arrays
 real(dp),allocatable::dummy(:,:),nrm(:)
 real(dp)::r(3,lmax2),rp(3,lmax2)
 real(dp)::rs2,rs3,rs6,rs12,umat(3,3)
 real(dp)::rot(lmax2,lmax2),tor(lmax2,lmax2),orb_lm(lmax2,-5:3,7)
 real(dp):: ylmrp(lmax2)
 real(dp):: ylmr_rr(lmax2,lmax2),ylmr_rr_save(lmax2,lmax2)
 real(dp):: ylmr_rrinv(lmax2,lmax2),ylmr_rp(lmax2,lmax2)
 character(len=500) :: message                   ! to be uncommented, if needed
!no_abirules
!integer :: orb_idx(16)=(/1,3,4,2,7,8,6,9,5,13,14,12,15,11,16,10/) !Tab3.1 Wannier90 user guide

! *************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!DEBUG
!write(std_out,*)'lmax ',lmax,'lmax2 ',lmax2
!write(std_out,*)'mband ',mband,'nwan ',nwan
!
!do iwan=1,nwan
!write(std_out,*)'iwan,proj_l, proj_m',proj_l(iwan),proj_m(iwan)
!write(std_out,*)'iwan,proj_x, proj_z',iwan,proj_x(:,iwan),proj_z(:,iwan)
!end do
!!END DEBUG

!constants for linear combinations of ylm's
 rs2=1._dp/sqrt(2._dp)
 rs3=1._dp/sqrt(3._dp)
 rs6=1._dp/sqrt(6._dp)
 rs12=1._dp/sqrt(12._dp)

!
!mapping lm coefficients for real spherical harmonics
!table 3.1 of Wannier90 user guide with real spherical harmonics in routine initylmr
!s, py,pz,px, dxy,dyz,dz2,dxz,dx2-y2, fy(3x2-y2),fxyz,fyz2,fz3,fxz2,
!fz(x2-y2),fx(x2-3y2)
!note: check ordering of f orbitals, it might be wrong

 tor(:,:)=0.d0
 lm=0
 do ll=0,lmax
   do mm=-ll,ll
     onem=(-1.d0)**mm
     lm=lm+1
     if(ll == 0) then
       tor(lm,lm)=1.d0
     else
       tor(lm,lm)=onem*1.d0
     end if
   end do !mm
 end do !ll
!do lm=1,16
!write(std_out,*)'tor lm=',lm,tor(:,lm)
!end do

!coefficients for basic wannier orbitals in Table 3.1 order
 orb_lm(:,:,:)=0.d0
 ii=0
 do ll=0,lmax
   do mr=1,2*ll+1
     ii=ii+1
     orb_lm(:,ll,mr)= tor(:,ii)
!    write(std_out,*)'ii',ii,'orb_lm',orb_lm(:,ll,mr)
   end do
 end do



!coefficients for linear combinations in table 3.2 order
 if(lmax>=1) then
!  s            px
   orb_lm(:,-1,1)=rs2*tor(:,1)+rs2*tor(:,4)
   orb_lm(:,-1,2)=rs2*tor(:,1)-rs2*tor(:,4)
!  s            px            py
   orb_lm(:,-2,1)=rs3*tor(:,1)-rs6*tor(:,4)+rs2*tor(:,2)
   orb_lm(:,-2,2)=rs3*tor(:,1)-rs6*tor(:,4)-rs2*tor(:,2)
   orb_lm(:,-2,3)=rs3*tor(:,1)+2._dp*rs6*tor(:,4)
!  s        px        py        pz
   orb_lm(:,-3,1)=half*(tor(:,1)+tor(:,4)+tor(:,2)+tor(:,3))
   orb_lm(:,-3,2)=half*(tor(:,1)+tor(:,4)-tor(:,2)-tor(:,3))
   orb_lm(:,-3,3)=half*(tor(:,1)-tor(:,4)+tor(:,2)-tor(:,3))
   orb_lm(:,-3,4)=half*(tor(:,1)-tor(:,4)-tor(:,2)+tor(:,3))
 end if
 if(lmax>=2) then
!  s            px            py
   orb_lm(:,-4,1)=rs3*tor(:,1)-rs6*tor(:,4)+rs2*tor(:,2)
   orb_lm(:,-4,2)=rs3*tor(:,1)-rs6*tor(:,4)-rs2*tor(:,2)
   orb_lm(:,-4,3)=rs3*tor(:,1)+2._dp*rs6*tor(:,4)
!  pz           dz2
   orb_lm(:,-4,4)= rs2*tor(:,3)+rs2*tor(:,7)
   orb_lm(:,-4,5)=-rs2*tor(:,3)+rs2*tor(:,7)
!  s            px            dz2         dx2-y2
   orb_lm(:,-5,1)=rs6*tor(:,1)-rs2*tor(:,4)-rs12*tor(:,7)+half*tor(:,9)
   orb_lm(:,-5,2)=rs6*tor(:,1)+rs2*tor(:,4)-rs12*tor(:,7)+half*tor(:,9)
!  s            py            dz2         dx2-y2
   orb_lm(:,-5,3)=rs6*tor(:,1)-rs2*tor(:,2)-rs12*tor(:,7)-half*tor(:,9)
   orb_lm(:,-5,4)=rs6*tor(:,1)+rs2*tor(:,2)-rs12*tor(:,7)-half*tor(:,9)
!  s            pz           dz2
   orb_lm(:,-5,5)=rs6*tor(:,1)-rs2*tor(:,3)+rs3*tor(:,7)
   orb_lm(:,-5,6)=rs6*tor(:,1)+rs2*tor(:,3)+rs3*tor(:,7)
 end if

!real wannier orbital coefficient array
 do iwan=1,nwan
   ylmr_fac(:,iwan)=orb_lm(:,proj_l(iwan),proj_m(iwan))
 end do


!setup to rotate ylmr_fac to new axes if called for
!skip if only s projetors are used
 if ( lmax>0 ) then
!  generate a set of nr=lmax2 random vetors
   idum=123456
   do ir=1,lmax2
     do ii=1,3
       r(ii,ir) = uniformrandom(idum)-0.5d0
     end do !ii
   end do !ir
   ABI_ALLOCATE(nrm,(lmax2))
   nrm(:)=sqrt(r(1,:)**2+r(2,:)**2+r(3,:)**2)**0.5
   call initylmr(lmax+1,1,lmax2,nrm,1,r(:,:),ylmr_rr_save(:,:),dummy)
   ylmr_rr(:,:)=ylmr_rr_save(:,:)
   do ir=1,lmax2
     ylmr_rr_save(ir,:)=ylmr_rr(:,ir)
   end do
   ABI_DEALLOCATE(nrm)

   ylmr_rrinv(:,:)=0.d0
   do ii=1,lmax2
     ylmr_rrinv(ii,ii)=1.d0
   end do !ii
!  calculate inverse of ylmr(ir,lm) matrix
   ylmr_rrinv(:,:)=ylmr_rr_save(:,:)
   call matrginv(ylmr_rrinv,lmax2,lmax2)

!  check that r points are independent (ie., that matrix inversion wasn't
!  too close to singular)
   ylmr_rr=matmul(ylmr_rrinv,ylmr_rr_save)
   test=0.d0
   do ii=1,lmax2
     ylmr_rr(ii,ii)=ylmr_rr(ii,ii)-1.d0
     do jj=1,lmax2
       test=max(abs(ylmr_rr(ii,jj)),test)
     end do !ii
   end do !jj
   if(test>tol8) then
     write(message, '(5a)' )&
&     '  matrix inversion error for wannier rotations',ch10,&
&     '  random vetors r(j,1:nr) are not all independent !! ',ch10,&
&     '  Action : re-seed uniformrandom or maybe just try again'
     ABI_ERROR(message)
   end if !test>tol8

!  end of the preliminaries, now to the rotations of the wannier orbitals
   do iwan=1,nwan
!    don't bother for s orbitals
     if(proj_l(iwan)==0) cycle
!    check for default axes and cycle if found
     if(proj_z(1,iwan)==0.d0 .and. proj_z(2,iwan)==0.d0 .and.&
&     proj_z(3,iwan)== 1.d0 .and. proj_x(1,iwan)==1.d0 .and.&
&     proj_x(2,iwan)==0.d0 .and. proj_x(3,iwan)==0.d0) cycle

!    get the u matrix that rotates the reference frame
     call rotmat(proj_x(:,iwan),proj_z(:,iwan),inversion_flag,umat)
!
!    find rotated r-vetors. Optional inversion
!    operation is an extension of the wannier90 axis-setting options
!    which only allow for proper axis rotations
     if(inversion_flag==1) then
       rp(:,:)= -matmul ( umat(:,:),  r(:,:) )
     else
       rp(:,:) = matmul ( umat(:,:) , r(:,:) )
     end if !inversion_flag

!    get the ylm representation of the rotated vetors
     ABI_ALLOCATE(nrm,(lmax2))
     nrm(:)=sqrt(rp(1,:)**2+rp(2,:)**2+rp(3,:)**2)**0.5
     call initylmr(lmax+1,1,lmax2,nrm,1,rp(:,:),ylmr_rp(:,:),dummy)
     ylmr_rr(:,:)=ylmr_rp(:,:)
     do ir=1,lmax2
       ylmr_rp(ir,:)=ylmr_rr(:,ir)
     end do
     ABI_DEALLOCATE(nrm)
!    the matrix product sum(ir) ylmr_rrinv(lm,ir)*ylmr_rp(ir,lm') gives the
!    the  lmXlm matrix representation of the coordinate rotation

     rot(:,:)=matmul(ylmr_rrinv(:,:),ylmr_rp(:,:))
!
!    now rotate the current wannier orbital
     ylmrp(:)=matmul(rot(:,:),ylmr_fac(:,iwan))
     ylmr_fac(:,iwan)=ylmrp(:)
   end do !iwan
 end if !lmax>0

end subroutine mlwfovlp_ylmfar
!!***

end module m_mlwfovlp
!!***
