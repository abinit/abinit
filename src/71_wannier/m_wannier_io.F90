!!****m* ABINIT/m_wannier_io
!! NAME
!!  m_wannier_io
!!
!! FUNCTION
!!  subroutines for writting Wannier90 related files.
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2022 ABINIT group (BAmadon, CEspejo, FJollet, TRangel, DRH, hexu)
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



module m_wannier_io
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
  use m_dtset, only : dataset_type

  use defs_datatypes, only : pseudopotential_type, ebands_t
  use defs_abitypes, only : MPI_type
  use m_fft, only: fourwf
  use m_io_tools,        only : open_file, get_unit

  use m_fftcore,  only : sphereboundary
  use m_dtfil
  use m_pawtab,   only : pawtab_type
  use m_pawcprj,  only : pawcprj_type
  use m_abstract_wf, only: abstract_wf, cg_cprj, wfd_wf


  implicit none
  private
  public :: write_Amn
  public :: compute_and_write_unk
  public :: write_eigenvalues
  public :: write_Mmn
  public :: read_chkunit
  !!***

contains


  ! Write amn file
  subroutine write_Amn(A_matrix, fname, nsppol, mband, nkpt, num_bands, nwan, band_in)
    ! TODO use the A_matrix sizes instead of the nsppol, mband, nkpt, nwan
    complex(dpc),pointer :: A_matrix(:,:,:,:)
    !type(dataset_type),intent(in) :: dtset
    logical, intent(in) :: band_in(:, :)
    integer, intent(in) :: nsppol, num_bands(nsppol), nwan(nsppol), mband, nkpt
    character(len=fnlen), intent(in) :: fname(nsppol)
    character(len=1000) :: message
    integer :: iun(nsppol)
    integer :: isppol, ikpt, iband, iwan, jband,  ii, jj

    ! below is copied/modified from m_mlwfovlp.F90
    do isppol=1,nsppol
       ! TODO : relpace this.
       if (open_file(trim(fname(isppol)), message, newunit=iun(isppol), &
            & form="formatted", status="unknown", action="write") /= 0) then
          ABI_ERROR(message)
       end if
       write(iun(isppol),*) 'Projections from Abinit : mband,nkpt,nwan. indices: iband1,iwan,ikpt'
       write(iun(isppol),*) num_bands(isppol),nkpt,nwan(isppol)
    end do

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
    do isppol=1,nsppol
       close(iun(isppol))
       write(message, '(3a)' ) &
            &         '   ',trim(fname(isppol)),' written'
       call wrtout(std_out,  message,'COLL')
    end do
    !
    !  Write down part of the matrix to the output file
    !  This is for the automatic tests
    !
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

  end subroutine write_Amn


   !-----------------------------------------------------------------------
   !@brief  Write down the unk matrix
   !-----------------------------------------------------------------------
  subroutine compute_and_write_unk(wfnname, usepaw, w90prtunk, &
       & mpi_enreg, ngfft, nsppol, nspinor,  &
       & nkpt, mband,  mpw, mgfftc, mkmem,  nprocs, rank, npwarr, &
       & band_in,  dtset, kg, mywfc)
    !TODO split the calculation of unk with writting.
    character(len=fnlen), intent(inout) :: wfnname
    integer ,intent(in) :: usepaw
    integer, intent(in) ::w90prtunk
    type(MPI_type),intent(in) :: mpi_enreg
    integer, intent(in) :: ngfft(18), nsppol, nspinor, nkpt, mpw, mgfftc, mband, mkmem
    integer, intent(in) :: nprocs, rank, npwarr(:)
    logical, intent(in) :: band_in(:, :)
    type(dataset_type),intent(in) :: dtset

    !integer, intent(in) :: iwav(:, :,:,:),kg(3,mpw*mkmem)
    integer, intent(in) :: kg(3,mpw*mkmem)
    !real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
    class(abstract_wf), intent(inout) ::  mywfc

    integer :: iun_plot
    integer :: isppol, ikpt, ikg, iband, ig
    integer :: n1, n2, n3, n4, n5, n6, cplex, mgfft, npw_k
    integer :: tim_fourwf
    real(dp) :: weight
    integer :: spacing, nband_inc(nsppol)
    integer,allocatable:: kg_k(:,:)
    real(dp),allocatable :: denpot(:,:,:), cwavef(:,:), fofgout(:,:),fofr(:,:,:,:)
    integer,allocatable :: gbound(:,:)
    integer :: n1tmp, n2tmp, n3tmp, jj1, jj2, jj3, ipw, ispinor
    character(len=1000) :: message

    ABI_UNUSED(nspinor)
    if(usepaw==1) then
       write(message, '( a,a,a,a,a,a,a,a,a)')ch10,&
            &     "   WARNING: The UNK matrices will not contain the correct wavefunctions ",ch10,&
            &     "   since we are just writing the plane wave contribution.",ch10,&
            &     "   The contribution from inside the spheres is missing. ",ch10,&
            &     "   However, these files can be used for plotting purposes",ch10
       call wrtout(std_out,  message,'COLL')
    end if
    !
    spacing = w90prtunk
    write(message, '( 8a,i3,2a)')ch10,&
         &   "   UNK files will be written.",ch10,&
         &   "   According to the chosen value of w90prtunk",ch10,&
         &   "   the wavefunctions are to be written ",ch10, &
         &   "   at every ", spacing," records.",ch10
    call wrtout(std_out,  message,'COLL')
    !
    ABI_MALLOC(kg_k,(3,mpw))
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
          ABI_MALLOC(denpot,(cplex*n4,n5,n6))
          ABI_MALLOC(cwavef,(2,npw_k))
          ABI_MALLOC(fofr,(2,n4,n5,n6))
          ABI_MALLOC(gbound,(2*mgfft+8,2))
          ABI_MALLOC(fofgout,(2,npw_k))


          !iun_plot=1000+ikpt+ikpt*(isppol-1)
          write(wfnname,'("UNK",I5.5,".",I1)') ikpt, isppol
          if (open_file(trim(wfnname), message, newunit=iun_plot, &
               & form="unformatted", status="unknown", action="write") /= 0) then
             ABI_ERROR(message)
          endif


          !      open (unit=iun_plot, file=wfnname,form='formatted')
          !       open(unit=iun_plot, file=wfnname,form='unformatted')
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
          call sphereboundary(gbound,mywfc%hdr%istwfk(ikpt),kg_k,mgfft,npw_k)
          write(std_out,*) "  writes UNK file for ikpt, spin=",ikpt,isppol
          denpot(:,:,:)=zero
          weight = one
          do iband=1,mband
             if(band_in(iband,isppol)) then
                ! TODO: check if this is the right order.
                do ispinor = 1, nsppol
                   do ipw = 1, npw_k
                      !do ig=1,npw_k*dtset%nspinor
                      !cwavef(1,ig)=cg(1,ipw+iwav(ispinor, iband,ikpt,isppol))
                      !cwavef(2,ig)=cg(2,ipw+iwav(ispinor, iband,ikpt,isppol))
                      ig = ipw + (ispinor-1)*npw_k
                      cwavef(1,ig)=mywfc%cg_elem(1,ipw, ispinor, iband, ikpt, isppol)
                      cwavef(2,ig)=mywfc%cg_elem(2,ipw, ispinor, iband, ikpt, isppol)
                   end do
                end do
                tim_fourwf=0
                call fourwf(cplex,denpot,cwavef,fofgout,fofr,&
                     &           gbound,gbound,mywfc%hdr%istwfk(ikpt),kg_k,kg_k,mgfft,&
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
          ABI_FREE(cwavef)
          ABI_FREE(fofr)
          ABI_FREE(gbound)
          ABI_FREE(denpot)
          ABI_FREE(fofgout)
          ikg=ikg+npw_k
          close(iun_plot)
       end do  ! ikpt
    end do  ! nsppol
    ABI_FREE(kg_k)
    !
    write(message, '(4a)' )ch10, &
         &   '   ','UNK files written',ch10
    call wrtout(std_out,  message,'COLL')
  end subroutine compute_and_write_unk


!-----------------------------------------------------------------------
!     This subroutine writes the eigenvalues in the w90 format
!     (see w90io.F90 in w90 package)
!-----------------------------------------------------------------------
  subroutine write_eigenvalues(filew90_eig,eigen, band_in,  eigenvalues_w, &
       & nsppol, nkpt, mband,  dtset, rank, master )
    integer, intent(in) ::  nsppol, nkpt, mband, rank, master
    logical, intent(in) :: band_in(:, :)
    real(dp), intent(in) :: eigen(mband,nkpt,nsppol)
    !real(dp), intent(in) :: eigen(:, :, :)
    real(dp), intent(inout) :: eigenvalues_w(:, :, :)
    type(dataset_type),intent(in) :: dtset
    character(len=fnlen) :: filew90_eig(nsppol)
    integer :: iun(nsppol), isppol, band_index, iband, jband, nband_k, ikpt
    character(len=1000) :: message
    !  Assign file unit numbers
    if(rank==master) then
       do isppol=1,nsppol
          if (open_file(trim(filew90_eig(isppol)), message, newunit=iun(isppol), &
               & form="formatted", status="unknown", action="write") /= 0) then
             ABI_ERROR(message)
          endif
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
                if(rank==master) then
                   write(iun(isppol), '(2i6,4x,f10.5)' ) &
                        & jband,ikpt,Ha_eV*eigen(iband, ikpt, isppol)
                end if
                !eigen(iband+band_index)
                !          Finish writing, now save eigenvalues
                eigenvalues_w(jband,ikpt,isppol)=Ha_eV*eigen(iband, ikpt, isppol)
                !eigen(iband+band_index)
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

  end subroutine write_eigenvalues


!-----------------------------------------------------------------------
!     This subroutine writes the overlap matrix in the w90 format
!     (see w90io.F90 in w90 package)
!-----------------------------------------------------------------------
  subroutine write_Mmn(filew90_mmn, band_in, cm1, ovikp, g1, M_matrix, &
       &  nkpt, nsppol,  nntot, mband, num_bands,  message, iam_master)
    ! input and output vars
    integer, intent(in) :: nsppol,  nntot, nkpt
    integer, intent(in) :: mband, num_bands(:)
    character(len=fnlen), intent(in) :: filew90_mmn(nsppol)
    logical, intent(in) :: band_in(mband, nsppol)
    integer,intent(in):: ovikp(:,:)
    integer, intent(in) :: g1(:, :, :)
    real(dp), intent(in) :: cm1(:,:,:,:,:,:)
    logical, intent(in) :: iam_master
    complex(dpc),intent(inout) :: M_matrix(:,:,:,:,:)


    ! temporary vars
    integer :: isppol, ikpt1, intot,ii, jj,jband1, iband1, jband2, iband2
    integer :: iun(nsppol)
    character(len=1000) :: message

     do isppol=1,nsppol !we write separate output files for each isppol
       iun(isppol)=220+isppol
       if( iam_master) then
          open(unit=iun(isppol),file=filew90_mmn(isppol),form='formatted',status='unknown')
          write(iun(isppol),*) "nnkp version 90"
          write(iun(isppol),*) num_bands(isppol),nkpt,nntot
       end if
     end do

   do isppol=1,nsppol
     do ikpt1=1,nkpt
       do intot=1,nntot
          if( iam_master) then
             write(iun(isppol),'(2i6,3x,3x,3i5)') ikpt1,ovikp(ikpt1,intot),(g1(jj,ikpt1,intot),jj=1,3)
          end if
         jband2=0
         do iband2=1,mband ! the first index is faster
           if(band_in(iband2,isppol)) then
             jband2=jband2+1
             jband1=0
             do iband1=1,mband
               if(band_in(iband1,isppol)) then
                 jband1=jband1+1
                 if(iam_master) write(iun(isppol),*) &
&                 cm1(1,iband1,iband2,intot,ikpt1,isppol),cm1(2,iband1,iband2,intot,ikpt1,isppol)
                 M_matrix(jband1,jband2,intot,ikpt1,isppol)=&
&                 cmplx(cm1(1,iband1,iband2,intot,ikpt1,isppol),cm1(2,iband1,iband2,intot,ikpt1,isppol), kind=dpc )
!                write(2211,*) ikpt1,intot,iband1,iband2
!                write(2211,*) cm1(1,iband1,iband2,intot,ikpt1,isppol),cm1(2,iband1,iband2,intot,ikpt1,isppol)
               end if ! band_in(iband1)
             end do ! iband1
           end if ! band_in(iband2)
         end do ! iband2
       end do !intot
     end do !ikpt
     if( iam_master ) then
       close(iun(isppol))
       write(message, '(3a)' )  '   ',trim(filew90_mmn(isppol)),' written'
       call wrtout(std_out,  message,'COLL')
     end if !rank==master
   end do !isppol


   if(iam_master) then
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
 end subroutine write_Mmn

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
 !string
 character(len=fnlen) :: fname
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


end module m_wannier_io
