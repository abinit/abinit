!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_tddft
!! NAME
!!  m_tddft
!!
!! FUNCTION
!!  Routines for computing excitation energies within TDDFT
!!
!! COPYRIGHT
!! Copyright (C) 1999-2019 ABINIT group (XG, JYR, MB, MBELAND, SHAMEL)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tddft

 use defs_basis
 use defs_abitypes
 use m_abicore
 use m_xmpi
 use m_errors
 use m_wffile
 use m_sort
#if defined HAVE_MPI2
 use mpi
#endif

 use m_io_tools, only : get_unit
 use m_symtk,    only : matr3inv
 use m_time,     only : timab
 use m_fftcore,  only : sphereboundary
 use m_spacepar, only : hartre
 use m_mpinfo,   only : proc_distrb_cycle
 use m_fft,      only : fourwf, fourdp

 implicit none

 private
!!***

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

 public :: tddft
!!***

contains

!!****f* m_tddft/tddft
!! NAME
!! tddft
!!
!! FUNCTION
!! Compute the excitation energies within TDLDA
!! from input wavefunctions, eigenenergies, and band occupations.
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=wf in G space
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  etotal=total energy of the ground-state (Ha)
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kxc(nfft,nkxc)=exchange-correlation kernel
!!  mband=maximum number of bands
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of the dielectric matrix
!!  mkmem=number of k-points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum allowed value for npw
!!  nfft=(effective) number of FFT grid points (for this processor)
!!       WARNING about parallelization: see below
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!                see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  nkxc=second dimension of the array kxc (see rhotoxc for a description)
!!  npwarr(nkpt)=number of planewaves at each k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  ucvol=unit cell volume (Bohr**3)
!!  wffnew=unit number for current wf disk file
!!
!! OUTPUT
!!  (only writing)
!!
!! WARNING:
!! This routine should not be parallelized on space for the time being,
!!    because the already existing parallelisation is not the usual one, found
!!    in the majority of ABINIT routines.
!!
!! NOTES
!! * Only accept nspinor=1, nsppol=1, nkpt=1 (Gamma point), and occopt<3
!!   (insulating occupation numbers).
!!   It is expected to make it work for nsppol=2 in the future.
!!
!! * For the oscillator strengths, see the paper
!!   ''Time-Dependent Density Functional Response Theory of Molecular
!!     systems: Theory, Computational Methods, and Functionals'', by M.E. Casida,
!!   in Recent Developments and Applications of Modern Density Functional
!!   Theory, edited by J.M. Seminario (Elsevier, Amsterdam, 1996).
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      fourdp,fourwf,hartre,matr3inv,mpi_bcast,mpi_gatherv,mpi_reduce
!!      mpi_scatterv,sort_dp,sphereboundary,timab,wrtout,xmpi_barrier
!!      xmpi_bcast,xmpi_exch,xmpi_sum,zhpev
!!
!! SOURCE

 subroutine tddft(cg,dtfil,dtset,eigen,etotal,gmet,gprimd,gsqcut,&
&  kg,kxc,mband,mgfftdiel,mkmem,mpi_enreg,mpw,nfft,ngfftdiel,nkpt,nkxc,&
&  npwarr,nspinor,nsppol,occ,ucvol,wffnew)

!Arguments ------------------------------------
 integer, intent(in) :: mband,mgfftdiel,mkmem,mpw,nfft,nkpt,nkxc,nsppol
 integer, intent(in) :: nspinor
 real(dp), intent(in) :: etotal,gsqcut,ucvol
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in) :: dtset
 type(MPI_type), intent(in) :: mpi_enreg
 type(wffile_type), intent(inout) :: wffnew
 integer, intent(in) :: kg(3,mpw*mkmem),ngfftdiel(18),npwarr(nkpt)
 real(dp), intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),eigen(mband*nkpt*nsppol)
 real(dp), intent(in) :: gmet(3,3),gprimd(3,3),kxc(nfft,nkxc),occ(mband*nkpt*nsppol)

!Local variables-------------------------------
 integer,parameter :: nexcitout=20
 integer :: cplex,i1,i2,i3,iband,idir,ier,ierr,iexcit,iexcit1,iexcit2,ifft
 integer :: old_iexcit,ii,jj,isppol,jsppol,isppol1,isppol2,isppol_l,isppol_n
 integer :: isppol_n1,isppol_n2,iocc_n1,iocc_n2,iunocc_n1,iunocc_n2,temp_unit2
 integer :: ikpt,index,iocc,iocc1,iocc2,iocc_l,iocc_n
 integer :: istwf_k,iunocc,iunocc1,iunocc2,iunocc_l,iunocc_n,istate,jexcit
 integer :: jexcit_cbase,master,mcg_disk,me_loc
 integer :: nband_k(nsppol), nband_occ(nsppol), nband_unocc(nsppol)
 integer :: nstate_k, nstate_occ, nstate_unocc, nexcit_pol(nsppol)
 integer :: nstate_win,ndiel,ndiel1,ndiel2,ndiel3,ndiel4
 integer :: ndiel5,ndiel6,nexcit,nexcit_max,nexcit_win,nfftdiel,nlargest,nnext
 integer :: nnext1,nnext2
 integer :: nproc_loc,npw_k,pole_approx,sing_trip,spaceComm,tim_fourwf
 integer :: tim_rwwf,save_iomode
 integer :: rec,recl,idummy,jdummy
 real(dp) :: buffer,buffer_inv,diffeig,eigunocc,emax_win
 real(dp) :: factor,ff,flargest,fnext,fr_invsquare,fr_power
 real(dp) :: fnext1,fnext2
 real(dp) :: normint,myproduct,saa,sab,sbb
 real(dp) :: sumx
 real(dp) :: sum_kernel(2/nsppol)
 real(dp) :: weight,xx
 logical :: am_master,file_exist
 logical, allocatable :: done_excit(:,:),done_sexc(:),done_sexc2(:)
 character(len=18) :: chain1,chain2
 character(len=500) :: message
 integer,allocatable :: flag_state_win(:),gbound(:,:),indarr(:),index_state(:)
 integer,allocatable :: kg_k(:,:)
 integer,allocatable :: excit_coords(:,:)
 integer :: count_to_do, count, displ, countmax, displmax
 integer :: ijexcit, ijexcit2, sendcount
 real(dp) :: f_sing_trip(2/nsppol),sendbuf(5-nsppol)
 real(dp) :: cauchy(7),poscart(3),rprimd(3,3),tsec(2),dummy(2,1)
 integer :: iomode,action,me,nmaster,sender,source,sread,sskip
 integer :: formeig,icg,ikg,nband_k_
 logical :: mydata, tmaster, swrite
 integer,allocatable ::  kg_disk(:,:)
 real(dp),allocatable :: cg_disk(:,:),cg_tmp(:,:)
 real(dp),allocatable :: cwavef(:,:),eexcit(:)
 real(dp),allocatable :: eexcit2(:)
 real(dp),allocatable :: matr(:)
 real(dp),allocatable :: kxc_for_tddft(:,:,:,:,:,:),omega_tddft_casida(:,:,:,:,:,:,:)
 real(dp),allocatable :: osc_str(:,:),pos(:,:),rhoaug(:,:,:),rhog(:,:)
 real(dp),allocatable :: sexc(:,:),sqrtks(:),vec(:,:,:),vhartr(:),wfprod(:,:,:)
 real(dp) :: omega_tddft_casida_dummy(2/nsppol)
 real(dp),allocatable :: wfraug(:,:,:,:),wfrspa(:,:,:,:),work(:),zhpev1(:,:)
 real(dp),allocatable :: zhpev2(:)
#if defined HAVE_MPI
 integer :: iproc
 integer :: ipwnbd
 integer,allocatable :: counts(:),displs(:),recvcounts(:),tmpbuf(:)
 real(dp), allocatable :: recvbuf(:,:)
#endif

! *************************************************************************

!Init mpi_comm
 spaceComm=mpi_enreg%comm_cell

 am_master=.true.
 master = 0
 nproc_loc = xmpi_comm_size(spaceComm) !Init ntot proc max
 me_loc    = xmpi_comm_rank(spaceComm) !Define who i am

#if defined HAVE_MPI
 if (me_loc/=0) then
   am_master=.FALSE.
 end if
 write(message, '(a,i3,a)' ) ' TDDFT ',nproc_loc,' CPU synchronized'
 call wrtout(std_out,message,'COLL')
 write(message, '(a,3D12.5,a,3D12.5,a,3D12.5)' ) ' gmet ',&
& gmet(1,1),gmet(1,2),gmet(1,3),ch10,&
& gmet(2,1),gmet(2,2),gmet(2,3),ch10,&
& gmet(3,1),gmet(3,2),gmet(3,3)
 call wrtout(std_out,message,'COLL')
#endif


!COMMENT these values should become arguments
!the two first define the energy window

 emax_win=greatest_real*tol6
 if(dtset%td_maxene>tol6)then
   emax_win = dtset%td_maxene
 end if

 call timab(95,1,tsec)

 istwf_k=dtset%istwfk(1)

 if(nkpt/=1 .or. &
& abs(dtset%kptns(1,1))+abs(dtset%kptns(2,1))+abs(dtset%kptns(3,1))>1.0d-6 )then
   write(message, '(a,a,a,a,a,i4,a,3es14.6,a,a,a,a,a)' )&
&   'The computation of excited states using TDDFT is only allowed',ch10,&
&   'with nkpt=1, kpt=(0 0 0), but the following values are input:',ch10,&
&   'nkpt=',nkpt,', kpt=',dtset%kptns(1:3,1),'.',ch10,&
&   'Action: in the input file, set nkpt to 1 and kpt to 0 0 0 ,',ch10,&
&   'or change iscf.'
   MSG_ERROR(message)
 end if

 if(nspinor/=1)then
   write(message, '(a,a,a,a,a,a,a)' )&
&   'The computation of excited states using TDDFT is restricted',ch10,&
&   'for the time being to nspinor=1, while input nspinor=2.',ch10,&
&   'Action: if you want to compute excited states within TDDFT,',ch10,&
&   'set nsppol to 1 in the input file. Otherwise, do not use iscf=-1.'
   MSG_ERROR(message)
 end if


 if(nsppol==2 .and. (dtset%ixc==22 .or. dtset%ixc==20))then
   write(message, '(a,a,a,a,a,a,a,a,a,a,a)' )&
&   'The computation of excited states using TDDFT in the spin',ch10,&
&   'polarized case for the time being cannot be used with ixc=20',ch10,&
&   'or ixc=22',ch10,&
&   'Action: if you want to compute excited states within TDDFT,',ch10,&
&   'set ixc different from 20 or 22. Otherwise, do not use iscf=-1',ch10,&
&   'with nsppol=2.'
   MSG_ERROR(message)
 end if


 if(dtset%occopt>2)then
   write(message, '(a,a,a,i2,a,a,a,a,a)' )&
&   'The computation of excited states using TDDFT is only allowed',ch10,&
&   'with occopt=0, 1, or 2, while input occopt=',dtset%occopt,'.',ch10,&
&   'Action: if you want to compute excited states within TDDFT,',ch10,&
&   'set occopt=0, 1, or 2 in the input file. Otherwise, do not use iscf=-1.'
   MSG_ERROR(message)
 end if

!Examine the occupation numbers, and determine the number of
!occupied and unoccupied states and band.
!States are numerated as usual in Abinit, before all spin up band
!and after all spin down bands.
!Note that if nsppol==1 nstate=nband_k
 do isppol=1,nsppol
   nband_k(isppol)=dtset%nband(isppol)
   nband_occ(isppol)=0
   do iband=1,nband_k(isppol)
     if(abs(occ(iband+(isppol-1)*nband_k(1))-two/nsppol)<tol6)  &
&     nband_occ(isppol)=nband_occ(isppol)+1
   end do
   nband_unocc(isppol)=nband_k(isppol)-nband_occ(isppol)
!  next line make no sense if spin flip is taken into account
   nexcit_pol(isppol)=nband_occ(isppol)*nband_unocc(isppol)
 end do
 nstate_k=nband_k(1)+(nsppol-1)*nband_k(nsppol)
 nstate_occ=nband_occ(1)+(nsppol-1)*nband_occ(nsppol)
 nstate_unocc=nstate_k-nstate_occ
!next line to be changed if spin fli is taken into account
 nexcit=nexcit_pol(1)+(nsppol-1)*nexcit_pol(nsppol)

!number of plane wave (does it work even for nsppol=2 ??)
 npw_k=npwarr(1)

!mux number of excitations that is taken into account
 if(dtset%td_mexcit==0)then
   nexcit_max=nexcit
 else
   nexcit_max =dtset%td_mexcit
 end if

!DEBUG
!write(std_out,*) nband_occ(1),nband_unocc(1)
!write(std_out,*) nband_occ(nsppol),nband_unocc(nsppol)
!END DEBUG


 if(nsppol==1)then
   write(message, '(a,a,a,a,i4,a,i4,a,a,i4,a,a,a,i6,a)' )ch10,&
&   ' *** TDDFT : computation of excited states *** ',ch10,&
&   ' Splitting of',dtset%nband(1),' states in',nband_occ(1),' occupied states,',&
&   ' and',nband_unocc(1),' unoccupied states,',ch10,&
&   ' giving',nexcit,' excitations.'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 else
   write(message, '(a,a,a,a,i4,a,i4,a,a,i4,a,a,a,i6,a,a,a)' )ch10,&
&   ' *** TDDFT : computation of excited states *** ',ch10,&
&   ' Splitting of',nstate_k,' states in',nstate_occ,' occupied states,',&
&   ' and',nstate_unocc,' unoccupied states,',ch10,&
&   ' giving',nexcit,' excitations. Note that spin flip is not possible actually.',ch10,&
&   ' So the number of excitation is the half of the product of the number of state'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 end if

!Allocate the matrices to be diagonalized.
!Use a simple storage mode, to be improved in the future.
 ii=max(nband_occ(1),nband_occ(nsppol))
 jj=max(nband_unocc(1),nband_unocc(nsppol))
 ABI_ALLOCATE(omega_tddft_casida,(ii,jj,nsppol,ii,jj,nsppol,2/nsppol))
 ABI_ALLOCATE(eexcit,(nexcit))
 ABI_ALLOCATE(sqrtks,(nexcit))
 ABI_ALLOCATE(flag_state_win,(nstate_k))
 omega_tddft_casida(:,:,:,:,:,:,:)=zero


!Fill the diagonal elements with square of differences of KS eigenvalues
!(also not very efficient, but OK for the present first coding)
!Also compute the square root of Kohn-Sham eigenvalue differences
 do isppol=1,nsppol
   do iunocc=1,nband_unocc(isppol)
     eigunocc=eigen(iunocc+nband_occ(isppol)+(isppol-1)*nband_k(1))
     do iocc=1,nband_occ(isppol)
       iexcit=iocc+(isppol-1)*nexcit_pol(1)+nband_occ(isppol)*(iunocc-1)
       diffeig=eigunocc-eigen(iocc+(isppol-1)*nband_k(1))
       do sing_trip=1,2/nsppol
         omega_tddft_casida(iocc,iunocc,isppol,iocc,iunocc,isppol,sing_trip)=diffeig**2
       end do
       eexcit(iexcit)=diffeig
       sqrtks(iexcit)=sqrt(diffeig)
     end do
   end do
 end do



!Sort the excitation energies : note that the array eexcit is reordered
 ABI_ALLOCATE(indarr,(nexcit))
 indarr(:)=(/ (ii,ii=1,nexcit) /)
 call sort_dp(nexcit,eexcit,indarr,tol14)

!Determine an energy window for the excitations
!to take into account. This is necessary for large systems

 nexcit_win = 0
 do iexcit = 1, nexcit
   if ((eexcit(iexcit) < emax_win ).and.(nexcit_win < nexcit_max)) then
     nexcit_win = nexcit_win + 1

!    DEBUG
!    write(message,'(a,F12.5,a,a,i2,a,a,i2)') 'excitation energy:', eexcit(indarr(iexcit)),ch10, &
!    &                                        'excitation number:', indarr(iexcit),ch10,         &
!    &                                        'nexcit_win:       ', nexcit_win
!    call wrtout(std_out,message,'COLL')
!    ENDDEBUG

   end if
 end do

!identification of the bands contributing to the
!nexcit_win  excitations within the window


 nstate_win = 0
 flag_state_win(:) = 0
 do iexcit = 1, nexcit_win
   iexcit1 = indarr(iexcit)
   isppol1 = min((iexcit1-1)/nexcit_pol(1) +1,2)
   iunocc1 = (iexcit1-(isppol1-1)*nexcit_pol(1)-1)/nband_occ(isppol1)+1
   iocc1   = iexcit1-(isppol1-1)*nexcit_pol(1)-(iunocc1-1)*nband_occ(isppol1)
   if (flag_state_win(nband_occ(isppol1)+(isppol1-1)*nband_k(1)+iunocc1)==0) &
&   flag_state_win(nband_occ(isppol1)+(isppol1-1)*nband_k(1)+iunocc1) =1
   if (flag_state_win(iocc1+(isppol1-1)*nband_k(1))==0) &
&   flag_state_win(iocc1+(isppol1-1)*nband_k(1)) =1
!  DEBUG
!  write(message,'(a,i2,a,a,i2,i2,a,a,i2,a,a,i2)') 'isppol:', isppol1,ch10, &
!  &                                       'iocc,iunocc:', iocc1,iunocc1,ch10,         &
!  &                                       'flag_state_win:', flag_state_win(iocc1+(isppol1-1)*nband_k(1)),  &
!  &                               ch10,   'flag_state_win:', flag_state_win(nband_occ(isppol1)+(isppol1-1)*nband_k(1)+iunocc1)
!  call wrtout(std_out,message,'COLL')
!  END DEBUG

 end do

 do isppol=1,nsppol
   do iband=1,nband_k(isppol)
     nstate_win=nstate_win+flag_state_win(iband+(isppol-1)*nband_k(1))
   end do
 end do


 write(message,'(a,a,i5)') ch10,'Nr of states to Fourier transform : ',nstate_win
 call wrtout(std_out,message,'COLL')

 ndiel1=ngfftdiel(1) ; ndiel2=ngfftdiel(2) ; ndiel3=ngfftdiel(3)
!ndiel4,ndiel5,ndiel6 are FFT dimensions, modified to avoid cache trashing
 ndiel4=ngfftdiel(4) ; ndiel5=ngfftdiel(5) ; ndiel6=ngfftdiel(6)

!The evaluation of integrals, later, needs the following factor
 normint=one/(ucvol*dble(ndiel1*ndiel2*ndiel3))

!Setup the positions in real space for later integration
 call matr3inv(gprimd,rprimd)

 ABI_ALLOCATE(pos,(max(ndiel1,ndiel2,ndiel3),3))

!Select the reduced position of the point with respect to the box center,
!in the interval ]-0.5,0.5].
 buffer=0.05_dp ; buffer_inv=one/buffer
 do idir=1,3
   if(idir==1)ndiel=ndiel1
   if(idir==2)ndiel=ndiel2
   if(idir==3)ndiel=ndiel3
   do ii=1,ndiel
!    dtset%boxcenter(3)=reduced coordinates of the center of the box,
!    in view of the computation of the oscillator strength
     pos(ii,idir)=(ii-1)/(one*ndiel)-dtset%boxcenter(idir)
     pos(ii,idir)=pos(ii,idir)-nint(pos(ii,idir)-tol12)
!    The linear behaviour is cut-off when one becomes
!    close to the boundaries : the buffer allows to match smoothly
!    one side of the cell to the other. This is important
!    to get rid of small breakings of symmetry, that are
!    confusing in accurate tests
     if(abs(pos(ii,idir))>half-buffer)then
!      xx is always positive, and goes linearly from 1 to 0
!      in the buffer region
       xx=(half-abs(pos(ii,idir)))*buffer_inv
!      The cut-off is applied to pos(:,:)
       pos(ii,idir)=pos(ii,idir)*xx*(two-xx)
!      DEBUG
!      if (idir==1)then
!      write(std_out,'(i2)') ndiel
!      write(std_out,'(a,i2,a,F12.5,F12.5)')'idiel : ',ii,'   x : ',pos(ii,idir),&
!      &    dtset%boxcenter(idir)
!      endif
!      ENDDEBUG
     end if
   end do ! ii
 end do ! idir

!need to run in MPI I/O case
 if (wffnew%iomode == IO_MODE_MPI ) then
   save_iomode=wffnew%iomode
   wffnew%iomode = IO_MODE_FORTRAN
 else
!  Do not store value but set to have save_iomode /= 1
   save_iomode = IO_MODE_FORTRAN
 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!2009 Chunping Hu
!need to collect wavefunctions from processors to master
 me=me_loc

 tim_rwwf =0
 source = master
 sread = master
 tmaster=(master==me)
 swrite=tmaster
 sender=-1

 iomode=wffnew%iomode

 if(am_master)then
#if defined HAVE_MPI
   ABI_ALLOCATE(cg_tmp,(2,mpw*nspinor*mband*nsppol))
#endif
 end if

 ABI_ALLOCATE(kg_disk,(3,mpw))
 mcg_disk=mpw*nspinor*mband
 formeig=0

#if defined HAVE_MPI
 call xmpi_barrier(spaceComm)
 ABI_ALLOCATE(cg_disk,(2,mcg_disk))
#endif

 icg=0
 if(mpi_enreg%paralbd==0) tim_rwwf=6
 if(mpi_enreg%paralbd==1)tim_rwwf=12

 do isppol=1,nsppol
   ikg=0
   do ikpt=1,nkpt
     nband_k_=dtset%nband(ikpt+(isppol-1)*nkpt)
     npw_k=npwarr(ikpt)
#if defined HAVE_MPI
     if (dtset%usewvl == 0) then
       call xmpi_barrier(spaceComm)
!      Must transfer the wavefunctions to the master processor
!      Separate sections for paralbd=1 or other values ; might be merged
       if(mpi_enreg%paralbd==0) then
         nmaster=0
         source=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k_,isppol))
         mydata=.false.
         if(source==me)mydata=.true.
         action=0
!        I am the master node, and I have the data in cg or cg_disk
         if((tmaster).and.(mydata))action=1
!        I am not the master, and I have the data => send to master
         if((.not.tmaster).and.(mydata))action=2
!        I am the master, and I receive the data
         if((tmaster).and.(.not.mydata))action=3
!        I have the data in cg or cg_disk ( MPI_IO case)
         if (iomode==IO_MODE_MPI) then
           action = 0
           sender=-1
           swrite=.false.
           if (mydata)then
             action=1
             swrite=.true.
             sender=me
           end if
         end if
!        I am the master node, and I have the data in cg or cg_disk
!        I have the data in cg or cg_disk ( MPI_IO case)
         if(action==1)then
!          Copy from kg to kg_disk
           kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
!          Copy from cg to cg_disk
           do ipwnbd=1,nband_k_*npw_k*nspinor
             cg_disk(1,ipwnbd)=cg(1,ipwnbd+icg)
             cg_disk(2,ipwnbd)=cg(2,ipwnbd+icg)
           end do
         end if
!        I am not the master, and I have the data => send to master
!        I am the master, and I receive the data
         if ( action==2.or.action==3) then
           call timab(48,1,tsec)
           if(action==2)then
             call xmpi_exch(kg(:,1+ikg:npw_k+ikg),3*npw_k,source,kg_disk,nmaster,spaceComm,ierr)
             call xmpi_exch(cg(:,icg+1:icg+nband_k_*npw_k*nspinor),2*nband_k_*npw_k*nspinor &
&             ,source,cg_disk,nmaster,spaceComm,ierr)
           else
             call xmpi_exch(kg_disk,3*npw_k,source,kg_disk,nmaster,spaceComm,ierr)
             call xmpi_exch(cg_disk,2*nband_k_*npw_k*nspinor,source,cg_disk,nmaster,spaceComm,ierr)
           end if
           call timab(48,2,tsec)
         end if
       else if(mpi_enreg%paralbd==1)then
         nmaster=0
#if defined HAVE_MPI_IO
         sender=-1
         if( iomode ==IO_MODE_MPI ) then
           nmaster=mpi_enreg%proc_distrb(ikpt,1,isppol)
           sender=nmaster
         end if
#endif
!        Note the loop over bands
         do iband=1,nband_k_
!          The message passing related to kg is counted as one band
           action=0
!          I am the master node, and I have the data in cg or cg_disk
           if( mpi_enreg%proc_distrb(ikpt,iband,isppol)==nmaster .and. me==nmaster) then
             action=1
!            I am not the master, and I have the data => send to master
           elseif( mpi_enreg%proc_distrb(ikpt,iband,isppol)==me .and. me/=nmaster ) then
             action = 2
!            I am the master, and I receive the data
           elseif( mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me .and. me==nmaster ) then
             action=3
           end if
           if(action==1) then
!            I am the master node, and I have the data in cg or cg_disk
!            Copy from kg to kg_disk
             if(iband==1)kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
!            Copy from cg to cg_disk
             do ipwnbd=1,npw_k*nspinor
               cg_disk(1,(iband-1)*npw_k*nspinor+ipwnbd)= cg(1,(iband-1)*npw_k*nspinor+ipwnbd+icg)
               cg_disk(2,(iband-1)*npw_k*nspinor+ipwnbd)= cg(2,(iband-1)*npw_k*nspinor+ipwnbd+icg)
             end do
           end if  ! action=1
           if ( action==2.or.action==3) then
!            action=2 :  I am not the master, and I have the data => send to master
!            action=3 :  I am the master, and I receive the data
             call timab(48,1,tsec)
             if ( iband == 1 ) then
               if (action==2) then
                 call xmpi_exch(kg(:,1+ikg:npw_k+ikg),3*npw_k,mpi_enreg%proc_distrb(ikpt,iband,isppol) &
&                 ,kg_disk,nmaster,spaceComm,ierr)
               else
                 call xmpi_exch(kg_disk,3*npw_k,mpi_enreg%proc_distrb(ikpt,iband,isppol)  &
&                 ,kg_disk,nmaster,spaceComm,ierr)
               end if
             end if       ! iband =1
             ipwnbd=(iband-1)*npw_k*nspinor
             if (action==2)then
               call xmpi_exch( cg(:,ipwnbd+icg+1:ipwnbd+icg+npw_k*nspinor),2*npw_k*nspinor &
&               ,mpi_enreg%proc_distrb(ikpt,iband,isppol)                    &
&               ,cg_disk(:,ipwnbd+1:ipwnbd+npw_k*nspinor),nmaster,spaceComm,ierr)
             else
               call xmpi_exch( cg_disk(:,ipwnbd+1:ipwnbd+npw_k*nspinor),2*npw_k*nspinor    &
&               ,mpi_enreg%proc_distrb(ikpt,iband,isppol)                    &
&               ,cg_disk(:,ipwnbd+1:ipwnbd+npw_k*nspinor),nmaster,spaceComm,ierr)
             end if
             call timab(48,2,tsec)
           end if        ! action=2 or action=3
           if(iomode ==IO_MODE_MPI) then
!            I have the data in cg or cg_disk
             swrite=.false.
             if (nmaster == me) then
               swrite=.true.
             end if
           end if
!          End of loop over bands
         end do
!        End of paralbd=1
       end if
     end if
#endif

!    The wavefunctions for the present k point and spin are stored into cg_tmp
     if(am_master)then
#if defined HAVE_MPI
       cg_tmp(:,icg+1:icg+nband_k_*npw_k*nspinor)=cg_disk(:,:)
#endif
     end if

     sskip=1
#if defined HAVE_MPI
     if (dtset%usewvl == 0) then
       sskip=0
       if(.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k_,isppol,me)))sskip=1
     end if
#endif
     if(sskip==1)then
       icg=icg+npw_k*nspinor*nband_k_
       ikg=ikg+npw_k
     end if

   end do ! ikpt
 end do ! isppol
 ABI_DEALLOCATE(kg_disk)
#if defined HAVE_MPI
 ABI_DEALLOCATE(cg_disk)
#endif
!!!!!!!end of collecting wavefunction to master!!!!!!


 if(am_master)then
!  -----------------------------------------------------------
!  The disk access is only done by master...

   ABI_ALLOCATE(gbound,(2*mgfftdiel+8,2))
   ABI_ALLOCATE(kg_k,(3,npw_k))

   ikpt=1
!  Only one k point
!  do isppol=1,nsppol
   kg_k(:,1:npw_k)=kg(:,1:npw_k)
   call sphereboundary(gbound,istwf_k,kg_k,mgfftdiel,npw_k)
!  enddo

 end if ! am_master

!need to run in MPI I/O case
 if ( save_iomode == 1 )  wffnew%iomode = IO_MODE_MPI
!call wrtout(std_out,'After reading the wavefunction','COLL')

!Use a simple implementation for the computation of the kernel elements
 if (am_master) then
   ABI_ALLOCATE(cwavef,(2,mpw))
   ABI_ALLOCATE(rhoaug,(ndiel4,ndiel5,ndiel6))
   ABI_ALLOCATE(wfraug,(2,ndiel4,ndiel5,ndiel6))
 end if
 ABI_ALLOCATE(index_state,(nstate_k))

! all real-space states are kept in memory
 ABI_ALLOCATE(wfrspa,(ndiel4,ndiel5,ndiel6,nstate_win))

!DEBUG
!write(message,'(a)') 'After allocating wfrspa'
!call wrtout(std_out,message,'COLL')
!ENDDEBUG

 weight=zero

!Generate states in real space, only for states contributing to excitations in window
 istate=0

 do isppol=1,nsppol
   do iband=1,nband_k(isppol)

     if(flag_state_win(iband+(isppol-1)*nband_k(1)) == 1) then
       istate=istate+1
       index_state(iband+(isppol-1)*nband_k(1))=istate

       if (am_master) then
#if defined HAVE_MPI
!        Obtain Fourier transform in fft box
         cwavef(:,1:npw_k)=cg_tmp(:,1+(iband-1)*npw_k+(isppol-1)* &
&         (npw_k*nband_k(1)) : iband*npw_k+(isppol-1)*(npw_k*nband_k(1)))
#else
         cwavef(:,1:npw_k)=cg(:,1+(iband-1)*npw_k+(isppol-1)* (npw_k*nband_k(1)) : iband*npw_k+(isppol-1)*(npw_k*nband_k(1)))
#endif

!        write(std_out,*)' iband : ',iband, ' isppol', isppol, '  -> index ', &
!        &            istate,index_state(iband+(isppol-1)*nband_k(1))

         tim_fourwf=14
!        This call should be made by master, and then the results be sent to the other procs

         call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&         istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg,1,ngfftdiel,npw_k,1,ndiel4,ndiel5,ndiel6,&
&         0,tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)

!        write(std_out,'(a,i5)')' After Fourier proc ',me_loc

!        Fix the phase, and checks that the wavefunction is real
!        (should be merged with routine fxphas)
         saa=zero ; sab=zero ; sbb=zero
         do i3=1,ndiel3
           do i2=1,ndiel2
             do i1=1,ndiel1
               saa=saa+wfraug(1,i1,i2,i3)**2
               sbb=sbb+wfraug(2,i1,i2,i3)**2
               sab=sab+wfraug(1,i1,i2,i3)*wfraug(2,i1,i2,i3)
             end do
           end do
         end do

         if(sbb>5.0d-9)then

           write(message, '(a,a,a,es20.10,a,i4,a,i2,a)' )&
&           'The imaginary part of wavefunctions should be practically zero.',ch10,&
&           'This is not the case, since sbb=',sbb,' for iband=',iband,'with sppol=',1,'.'
           MSG_WARNING(message)
           if(sbb>1.0d-7)then
             MSG_ERROR("sbb>1.0d-7")
           end if
         end if

!        Possibility of writing to disk

         wfrspa(:,:,:,istate)=wfraug(1,:,:,:)
       end if !am_master

     end if

!    End loop on iband
   end do
!  End loop on nsppol
 end do

 if (am_master) then
   ABI_DEALLOCATE(gbound)
   ABI_DEALLOCATE(kg_k)
   ABI_DEALLOCATE(cwavef)
   ABI_DEALLOCATE(rhoaug)
   ABI_DEALLOCATE(wfraug)
#if defined HAVE_MPI
   ABI_DEALLOCATE(cg_tmp)
#else
#endif
 end if
 ABI_DEALLOCATE(flag_state_win)

! send wfrspa from master to world
 call xmpi_bcast(wfrspa,master,spaceComm,ierr)
 !call MPI_BCAST(wfrspa,nbuf,MPI_DOUBLE_PRECISION,master,spaceComm,ierr)


!DEBUG
!#if defined HAVE_MPI
!call xmpi_barrier(spaceComm)
!write(message,'(a)')' after iband loop synchronization done...'
!call  wrtout(std_out,message,'COLL')
!#endif
!ENDDEBUG

!Compute the xc kernel, in the form needed for the singlet or triplet
!excitation energy.
!In the case ixc=20, kxc vanishes, but no change is made here, for simplicity.
!(ixc=20 implemented only in the not spin polyrized case)

!DEBUG
!write(std_out,*)' tddft : xc kernel '
!do ifft=1,nkxc,41
!write(std_out,*)ifft,kxc(ifft,1),kxc(ifft,2),kxc(ifft,3)
!enddo
!stop
!ENDDEBUG

 ABI_ALLOCATE(kxc_for_tddft,(ndiel1,ndiel2,ndiel3,nsppol,nsppol,2/nsppol))
 if(dtset%ixc/=22)then
   do isppol=1,nsppol
     do jsppol=1,nsppol
       index=1
       do i3=1,ndiel3
         do i2=1,ndiel2
           do i1=1,ndiel1
             do sing_trip=1,2/nsppol
               kxc_for_tddft(i1,i2,i3,isppol,jsppol,sing_trip)=two/nsppol* &
&               (kxc(index,isppol+jsppol-1)-(sing_trip-1)*kxc(index,2))
             end do
             index=index+1
           end do
         end do
       end do
     end do
   end do
 else
!  This is for the Burke-Petersilka-Gross hybrid, with ixc=22
!  However, the implementation in case of spin-polarized system should not be expected to be the correct one !
   do isppol=1,nsppol
     do jsppol=1,nsppol
       index=1
       do i3=1,ndiel3
         do i2=1,ndiel2
           do i1=1,ndiel1
             do sing_trip=1,2/nsppol
               kxc_for_tddft(i1,i2,i3,isppol,jsppol,sing_trip)=((-1)**(sing_trip+1))*kxc(index,2)
             end do
             index=index+1
           end do
         end do
       end do
     end do
   end do
 end if

 pole_approx=0

 ABI_ALLOCATE(excit_coords,(nexcit_win**2,2))

#if defined HAVE_MPI
 ABI_ALLOCATE(counts,(0:nproc_loc-1))
 ABI_ALLOCATE(displs,(0:nproc_loc-1))
 ABI_ALLOCATE(recvcounts,(0:nproc_loc-1))
 ABI_ALLOCATE(recvbuf,(5-nsppol,nproc_loc-1))
#endif

!DEBUG
!write(std_out,*)'before first loop'
!ENDDEBUG

!0000000000000000000000000000000000000000000000000000000
!check if matrix file fname_tdexcit exists on disk
!if the file is present,calculation is a continuation
 if (am_master) then
   inquire(file=trim(dtfil%fnametmp_tdexcit),exist=file_exist)
!  for direct access to excitation file
   if(nsppol==1)then
     inquire(iolength=recl) omega_tddft_casida(1,1,1,1,1,1,1),omega_tddft_casida(1,1,1,1,1,1,1), iexcit,jexcit
   else
     inquire(iolength=recl) omega_tddft_casida(1,1,1,1,1,1,1), iexcit,jexcit
   end if

   temp_unit2 = get_unit()
   open(temp_unit2, file=trim(dtfil%fnametmp_tdexcit),form='unformatted', recl=recl, access='DIRECT')

   ABI_ALLOCATE(done_excit,(nexcit_win,nexcit_win))

   if(file_exist)then
     write(std_out,*)'TDDFT continues from a previous run'
     rec=0
     do iexcit=1,nexcit_win
       iexcit2 = indarr(iexcit)
       isppol2 = min((iexcit2-1)/nexcit_pol(1) +1,2)
       iunocc2 = (iexcit2-(isppol2-1)*nexcit_pol(1)-1)/nband_occ(isppol2)+1
       iocc2   = iexcit2-(isppol2-1)*nexcit_pol(1)-(iunocc2-1)*nband_occ(isppol2)
       do jexcit=1,nexcit_win
         iexcit1 = indarr(jexcit)
         isppol1 = min((iexcit1-1)/nexcit_pol(1) +1,2)
         iunocc1 = (iexcit1-(isppol1-1)*nexcit_pol(1)-1)/nband_occ(isppol1)+1
         iocc1   = iexcit1-(isppol1-1)*nexcit_pol(1)-(iunocc1-1)*nband_occ(isppol1)

         rec=rec+1 ! record of the entry in the excitation file
         if(nsppol==1)then
           read(temp_unit2,rec=rec) omega_tddft_casida_dummy(1), omega_tddft_casida_dummy(2), idummy, jdummy
         else
           read(temp_unit2,rec=rec) omega_tddft_casida_dummy(1), idummy, jdummy
         end if
         done_excit(jexcit,iexcit)= ( idummy /= -1 .and. jdummy /= -1 ) ! if true, eigsqr_singlet and eigsqr_triplet are ok
!        and a true is marked in the logical array done_excit
         if (done_excit(jexcit,iexcit)) then
           do sing_trip=1,2/nsppol
             omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,sing_trip)= &
&             omega_tddft_casida_dummy(sing_trip)
           end do
         end if
       end do
     end do

   else
     write(std_out,*)'no excitation matrix on disk'
     write(std_out,*)'TDDFT starts from scratch'
     done_excit = .false. ! initialize the logical array to false
     rec=0
     do iexcit=1,nexcit_win
       do jexcit=1,nexcit_win
         rec=rec+1
         if(nsppol==1)then
           write(temp_unit2,rec=rec) zero,zero,-1,-1
         else
           write(temp_unit2,rec=rec) zero,-1,-1
         end if
       end do
     end do
   end if

!  Need to list the elements to compute, taking the symmetry into account: valid only for Gamma point but this is already the case
   count_to_do=0
   do iexcit=1,nexcit_win
     do jexcit=1,iexcit
       if (.not. done_excit(jexcit,iexcit)) then
         count_to_do=count_to_do+1
         excit_coords(count_to_do,1)=iexcit
         excit_coords(count_to_do,2)=jexcit
       end if
     end do
   end do

   ABI_DEALLOCATE(done_excit)



#if defined HAVE_MPI
!  Compute limits for load balancing
   do iproc=0,nproc_loc-1
     displs(iproc)=(iproc*count_to_do)/nproc_loc
     counts(iproc)=min(((iproc+1)*count_to_do)/nproc_loc,count_to_do)-displs(iproc)
   end do
#endif

 end if ! am_master


#if defined HAVE_MPI
 call MPI_BCAST(count_to_do,1,MPI_INTEGER,master,spaceComm,ierr)
#endif

 displ=(me_loc*count_to_do)/nproc_loc
 count=min(((me_loc+1)*count_to_do)/nproc_loc,count_to_do)-displ
 displmax=((nproc_loc-1)*count_to_do)/nproc_loc
 countmax=count_to_do-displmax

 write(message,'(A,I6)') 'Maximum number of matrix elements per processor = ',countmax
 call wrtout(std_out,message,'COLL')


#if defined HAVE_MPI
!Need to dispatch the elements to compute to the different processes
 ABI_ALLOCATE(tmpbuf,(nexcit_win**2))
 tmpbuf=0
 call MPI_Scatterv(excit_coords(1,1),counts,displs,MPI_INTEGER,tmpbuf,count,MPI_INTEGER,0,spaceComm,ierr)
 excit_coords(:,1)=tmpbuf(:)
 tmpbuf=0
 call MPI_Scatterv(excit_coords(1,2),counts,displs,MPI_INTEGER,tmpbuf,count,MPI_INTEGER,0,spaceComm,ierr)
 excit_coords(:,2)=tmpbuf(:)
 ABI_DEALLOCATE(tmpbuf)
#endif

 nfftdiel=ndiel1*ndiel2*ndiel3
 ABI_ALLOCATE(wfprod,(ndiel1,ndiel2,ndiel3))
 ABI_ALLOCATE(work,(nfftdiel))
 ABI_ALLOCATE(sexc,(3,nexcit_win))
 ABI_ALLOCATE(done_sexc,(nexcit_win))
 ABI_ALLOCATE(done_sexc2,(nexcit_win))
 ABI_ALLOCATE(rhog,(2,nfftdiel))
 ABI_ALLOCATE(vhartr,(nfftdiel))

 sexc(:,:)=zero
 done_sexc(:)=.false.

!----------------------------------------------------------
!Main double loop

 old_iexcit=0
 do ijexcit=1,countmax
!  we really loop only through count, but we need to go through countmax
!  to make sure that all processes execute MPI_Gatherv below
   if (ijexcit <= count) then
     iexcit=excit_coords(ijexcit,1)
     jexcit=excit_coords(ijexcit,2)

     iexcit2 = indarr(iexcit)
     isppol2 = min((iexcit2-1)/nexcit_pol(1) +1,2)
     iunocc2 = (iexcit2-(isppol2-1)*nexcit_pol(1)-1)/nband_occ(isppol2)+1
     iocc2   = iexcit2-(isppol2-1)*nexcit_pol(1)-(iunocc2-1)*nband_occ(isppol2)

     iexcit1 = indarr(jexcit)
     isppol1 = min((iexcit1-1)/nexcit_pol(1) +1,2)
     iunocc1 = (iexcit1-(isppol1-1)*nexcit_pol(1)-1)/nband_occ(isppol1)+1
     iocc1   = iexcit1-(isppol1-1)*nexcit_pol(1)-(iunocc1-1)*nband_occ(isppol1)

     if (old_iexcit /= iexcit) then
!      We start a new column of the matrix
!      DEBUG
!      write(message,'(a,i5,a,i3)')'treating  iexcit =  ',iexcit,&
!      &'  with proc ',me_loc
!      call wrtout(std_out,message,'PERS')
!      ENDDEBUG

!      DEBUG
!      write(message,'(a,i3)')'Multiplicating phi1 phi2, on proc ',me_loc
!      call wrtout(std_out,message,'PERS')
!      ENDDEBUG

       ifft=1
       do i3=1,ndiel3
         do i2=1,ndiel2
           do i1=1,ndiel1
             wfprod(i1,i2,i3)=wfrspa(i1,i2,i3,index_state(iocc2+(isppol2-1)*nband_k(1))) &
&             *wfrspa(i1,i2,i3,index_state(iunocc2+nband_occ(isppol2)+(isppol2-1)*nband_k(1)))
             work(ifft)=wfprod(i1,i2,i3)
             ifft=ifft+1
           end do
         end do
       end do
       if (jexcit == 1) then
         do i3=1,ndiel3
           do i2=1,ndiel2
             do i1=1,ndiel1
               do idir=1,3
                 poscart(idir)=rprimd(idir,1)*pos(i1,1)+&
&                 rprimd(idir,2)*pos(i2,2)+&
&                 rprimd(idir,3)*pos(i3,3)
                 sexc(idir,iexcit)=sexc(idir,iexcit)+poscart(idir)*wfprod(i1,i2,i3)
               end do
             end do
           end do
           done_sexc(iexcit)=.true.
         end do
       end if

!      For the singlet correction, must compute the hartre potential created
!      by the product of wavefunctions
       cplex=1

!      DEBUG
!      write(message,'(a,i3)')'Before Fourdp, on proc ',me_loc
!      call wrtout(std_out,message,'PERS')
!      ENDDEBUG

       call fourdp(cplex,rhog,work,-1,mpi_enreg,nfftdiel,1,ngfftdiel,0)

!      DEBUG
!      write(message,'(a,i3)')'Before Hartree, on proc ',me_loc
!      call wrtout(std_out,message,'PERS')
!      write(std_out,*)'CPU ',me_loc,ch10,&
!      &            '   cplex : ',cplex,ch10,&
!      &            '   gmet(3,3)  : ',gmet(3,3),ch10,&
!      &            '   gsqcut : ',gsqcut,ch10,&
!      &            '   rhog(1,1) :,',rhog(1,1),ch10,&
!      &            '   vhartr(1) :,',vhartr(1)
!      ENDDEBUG

       call hartre(cplex,gsqcut,0,mpi_enreg,nfftdiel,ngfftdiel,rhog,rprimd,vhartr)

!      DEBUG
!      write(message,'(a,i3)')'After Hartree, on proc ',me_loc
!      call wrtout(std_out,message,'PERS')
!      ENDDEBUG
     end if
     old_iexcit=iexcit

!    DEBUG
!    write(std_out,*)'  treating  iexcit =  ',jexcit
!    write(std_out,*)'   indarr(iexcit) =',iexcit1,iocc1,iunocc1
!    write(std_out,*)'   index ',index_state(iocc1+(isppol1-1)*nband_k(1)), &
!    &                       index_state(iunocc1+nband_occ(isppol1)+(isppol1-1)*nband_k(1))
!    ENDDEBUG

     if(pole_approx==0 .or. (iunocc1==iunocc2 .and. iocc1==iocc2 .and. isppol1==isppol2))then
       sum_kernel(:)=zero
       f_sing_trip(1)=two/dble(nsppol)
       if(nsppol==1) f_sing_trip(2)=zero
!      For Fermi-Amaldi kxc, the xc contribution is -1/2 the Hartree contribution
!      to the triplet state. The following factors combines both contributions.
       if(dtset%ixc==20 .or. dtset%ixc==22)then
         if(nsppol==1)then
           f_sing_trip(1)= one
           f_sing_trip(2)=-one
         end if
       end if
       ifft=1
       do i3=1,ndiel3
         do i2=1,ndiel2
           do i1=1,ndiel1
             myproduct=wfrspa(i1,i2,i3,index_state(iocc1+(isppol1-1)*nband_k(1))) &
&             *wfrspa(i1,i2,i3,index_state(iunocc1+nband_occ(isppol1)+(isppol1-1)*nband_k(1)))
             do sing_trip=1,2/nsppol
               sum_kernel(sing_trip)=sum_kernel(sing_trip)+&
&               myproduct*(f_sing_trip(sing_trip)*vhartr(ifft)+kxc_for_tddft(i1,i2,i3,isppol1,isppol2,sing_trip) &
               *wfprod(i1,i2,i3))
             end do
             ifft=ifft+1
           end do ! i1
         end do ! i2
       end do ! i3

!      The factor two is coherent with the formulas of Vasiliev et al
       factor=two*sqrtks(iexcit1)*sqrtks(iexcit2)*normint
       do sing_trip=1,2/nsppol
         omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,sing_trip)=   &
&         omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,sing_trip)+  &
&         factor*sum_kernel(sing_trip)
       end do

!      End condition of being diagonal element if pole approximation
     end if

!    Continue writing excitation matrix
     if (am_master) then
!      the master writes its results to disk
       if(nsppol==1)then
         write(temp_unit2, rec=(iexcit-1)*nexcit_win+jexcit ) &
&         omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1),&
&         omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,2),&
&         iexcit, jexcit
       else
         write(temp_unit2, rec=(iexcit-1)*nexcit_win+jexcit ) &
&         omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
&         iexcit, jexcit
       end if

!      DEBUG
!      if(nsppol==1)then
!      write(std_out,*)'singlet: ',omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1),&
!      &             'iexcit: ',iexcit,'jexcit :',jexcit
!      write(std_out,*)'triplet: ',omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,2),&
!      &             'iexcit: ',iexcit,'jexcit :',jexcit
!      else
!      write(std_out,*)'excitation: ',omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1),&
!      &             'iexcit: ',iexcit,'jexcit :',jexcit
!      endif
!      ENDDEBUG

       sendcount=0
     else
       sendcount=5-nsppol

       if(nsppol==1)then
         sendbuf=(/ omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
&         omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,2), &
&         real(iexcit,dp), real(jexcit,dp) /)
       else
         sendbuf=(/ omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
&         real(iexcit,dp), real(jexcit,dp) /)
       end if
     end if ! am_master
   else
!    ijexcit > count

!    done with local work, so send message of zero length
     sendcount=0

   end if ! ijexcit <= count

#if defined HAVE_MPI
   if (am_master) then

!    Compute displacements and counts for the gathering of the results
     displs(0)=0
     recvcounts(0)=0
     do iproc=1,nproc_loc-1
       recvcounts(iproc)=min(((iproc+1)*count_to_do)/nproc_loc,count_to_do)-(iproc*count_to_do)/nproc_loc
       if (recvcounts(iproc) < countmax .and. ijexcit==countmax) then
         recvcounts(iproc)=0
       else
         recvcounts(iproc)=5-nsppol
       end if
       displs(iproc)=displs(iproc-1)+recvcounts(iproc-1)
     end do
   end if

!  ***********************************************
!  ***** I have to ask about that ****************
!  ***********************************************
   call MPI_Gatherv(sendbuf,sendcount,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs, &
&   MPI_DOUBLE_PRECISION,0,spaceComm,ierr)

   if (am_master) then

!    Extract eigsqr_singlet, eigsqr_triplet, iexcit, jexcit from receive buffer and
!    write to file
     do ijexcit2=1,sum(recvcounts)/(5-nsppol)
       iexcit=int(recvbuf(4-nsppol,ijexcit2))
       jexcit=int(recvbuf(5-nsppol,ijexcit2))

       iexcit2 = indarr(iexcit)
       isppol2 = min((iexcit2-1)/nexcit_pol(1) +1,2)
       iunocc2 = (iexcit2-(isppol2-1)*nexcit_pol(1)-1)/nband_occ(isppol2)+1
       iocc2   = iexcit2-(isppol2-1)*nexcit_pol(1)-(iunocc2-1)*nband_occ(isppol2)

       iexcit1 = indarr(jexcit)
       isppol1 = min((iexcit1-1)/nexcit_pol(1) +1,2)
       iunocc1 = (iexcit1-(isppol1-1)*nexcit_pol(1)-1)/nband_occ(isppol1)+1
       iocc1   = iexcit1-(isppol1-1)*nexcit_pol(1)-(iunocc1-1)*nband_occ(isppol1)

       do sing_trip=1,2/nsppol
         omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,sing_trip)= &
&         recvbuf(sing_trip,ijexcit2)
       end do

       if(nsppol==1)then
         write(temp_unit2, rec=(iexcit-1)*nexcit_win+jexcit ) &
&         omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
&         omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,2), &
&         iexcit, jexcit
       else
         write(temp_unit2, rec=(iexcit-1)*nexcit_win+jexcit ) &
         omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
&         iexcit, jexcit
       end if
!      DEBUG
!      if(nsppol==1)then
!      write(std_out,*)'singlet: ',omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
!      &                       'iexcit: ',iexcit,'jexcit :',jexcit
!      write(std_out,*)'triplet: ',omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,2),
!      &                       'iexcit: ',iexcit,'jexcit :',jexcit
!      else
!      write(std_out,*)'excitation: ',omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
!      &                        'iexcit: ',iexcit,'jexcit :',jexcit
!      endif
!      ENDDEBUG

     end do

   end if
#endif

!  End indices loops
 end do ! ijexcit

!End of the main double loop
!--------------------------------------------------------------------


#if defined HAVE_MPI
!sexc needs to be summed here since it used only by master
 call xmpi_barrier(spaceComm)
!call xmpi_sum_master(sexc,master,spaceComm,ierr) ! Does not work on some machines
 call xmpi_sum(sexc,spaceComm,ierr)
 done_sexc2=done_sexc
 call MPI_Reduce(done_sexc2,done_sexc,nexcit_win,MPI_LOGICAL,MPI_LOR,master,spaceComm,ierr)
#endif


 if (am_master) then
!  We compute sexc again if it was not done. Will only be executed if
!  there was a restart from values read from logical unit temp_unit2.

   do iexcit=1,nexcit_win

!    DEBUG
!    write(std_out,*)'do on excitation',iexcit
!    END DEBUG

     if (.not.done_sexc(iexcit)) then
       iexcit2 = indarr(iexcit)
       isppol2 = min((iexcit2-1)/nexcit_pol(1) +1,2)
       iunocc2 = (iexcit1-(isppol2-1)*nexcit_pol(1)-1)/nband_occ(isppol2)+1
       iocc2   = iexcit2-(isppol2-1)*nexcit_pol(1)-(iunocc2-1)*nband_occ(isppol2)
       do i3=1,ndiel3
         do i2=1,ndiel2
           do i1=1,ndiel1
             wfprod(i1,i2,i3)=wfrspa(i1,i2,i3,index_state(iocc2+(isppol2-1)*nband_k(1))) &
&             *wfrspa(i1,i2,i3,index_state(iunocc2+nband_occ(isppol2)+ (isppol2-1)*nband_k(1)))
             do idir=1,3
               poscart(idir)=rprimd(idir,1)*pos(i1,1)+&
&               rprimd(idir,2)*pos(i2,2)+&
&               rprimd(idir,3)*pos(i3,3)
               sexc(idir,iexcit)=sexc(idir,iexcit)+poscart(idir)*wfprod(i1,i2,i3)
             end do
           end do
         end do
       end do
     end if
   end do
 end if

 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(pos)
 ABI_DEALLOCATE(vhartr)
 ABI_DEALLOCATE(kxc_for_tddft)
 ABI_DEALLOCATE(wfprod)
 ABI_DEALLOCATE(index_state)
 ABI_DEALLOCATE(excit_coords)
 ABI_DEALLOCATE(wfrspa)

!Write the first excitation energies
 write(message, '(a,a,es18.8,a,a,a,a,a,a,a,a,a)' )ch10,&
& '  Ground state total energy (Ha) :',etotal,ch10,ch10,&
& '  Kohn-Sham energy differences,',ch10,&
& '  corresponding total energies and oscillator strengths (X,Y,Z and average)-',ch10,&
& '  (oscillator strengths smaller than 1.e-6 are set to zero)',ch10,&
& '  Transition  (Ha)  and   (eV)   Tot. Ene. (Ha)  Aver     XX       YY       ZZ'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

#if defined HAVE_MPI
 ABI_DEALLOCATE(counts)
 ABI_DEALLOCATE(displs)
 ABI_DEALLOCATE(recvbuf)
 ABI_DEALLOCATE(recvcounts)
#endif

 if (am_master) then

   ABI_ALLOCATE(osc_str,(7,nexcit))

   do iexcit=1,nexcit_win
     iexcit2 = indarr(iexcit)
     isppol = min((iexcit2-1)/nexcit_pol(1) +1,2)
     iunocc = (iexcit2-(isppol-1)*nexcit_pol(1)-1)/nband_occ(isppol)+1
     iocc   = iexcit2-(isppol-1)*nexcit_pol(1)-(iunocc-1)*nband_occ(isppol)

     osc_str(1,iexcit)=zero
     do idir=1,3
!      One of the factor of two comes from the spin degeneracy,
!      the other comes from Eq.(40) of Casida
       osc_str(idir+1,iexcit)=&
&       (sexc(idir,iexcit)*sqrtks(iexcit2)*normint*ucvol)**2*two*two/nsppol
       osc_str(1,iexcit)=osc_str(1,iexcit)&
&       +osc_str(idir+1,iexcit)*third
     end do
     do ii=1,4
       if(abs(osc_str(ii,iexcit))<tol6)osc_str(ii,iexcit)=zero
     end do
!    Changed, the whole spectrum is written
!    The array eexcit has been reordered previously, the others also
     if(nsppol==1)then
       write(message, '(i4,a,i3,2es12.5,es13.5,es11.4,3es9.2)' ) &
&       iocc,'->',iunocc+nband_occ(isppol),            &
&       eexcit(iexcit), eexcit(iexcit)*Ha_eV,eexcit(iexcit)+etotal, &
!      XG 020209 : Jean-Yves, I assume that the printout of sexc is for debugging ?!
!      &   osc_str(1:4,iexcit),sexc(1:3,iexcit)
&       osc_str(1:4,iexcit)
     else
       write(message, '(i4,a,i3,a,i1,2es12.5,es13.5,es11.4,3es9.2)' ) &
&       iocc,'->',iunocc+nband_occ(isppol),' s:',isppol,            &
&       eexcit(iexcit), eexcit(iexcit)*Ha_eV,eexcit(iexcit)+etotal, &
&       osc_str(1:4,iexcit)
     end if
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

   end do

!  Check of the Sum rule for Casida eq.47,
!  only exact if complete basis of excitations, as well as local potentials only.
   sumx=zero
   do iexcit=1,nexcit_win
     sumx=sumx+osc_str(1,iexcit)
   end do
   write(message, '(a,es16.6)' )'  Sum of osc. strength : ',sumx
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

!  -Diagonalize the excitation matrices----------------------------

   ABI_ALLOCATE(eexcit2,(nexcit_win))
   ABI_ALLOCATE(vec,(2,nexcit_win,nexcit_win))

   do sing_trip=1,2/nsppol

     if(pole_approx==0)then

       ABI_ALLOCATE(matr,(nexcit_win*(nexcit_win+1)))
       ABI_ALLOCATE(zhpev1,(2,2*nexcit_win-1))
       ABI_ALLOCATE(zhpev2,(3*nexcit_win-2))
       matr(:)=zero
       ier=0
!      DEBUG
!      write(std_out,*)' after allocation matrices     '
!      ENDDEBUG


!      Store the matrix in proper mode before calling zhpev

       index=1
       do iexcit=1,nexcit_win
         iexcit2 = indarr(iexcit)
         isppol2 = min((iexcit2-1)/nexcit_pol(1) +1,2)
         iunocc2 = (iexcit2-(isppol2-1)*nexcit_pol(1)-1)/nband_occ(isppol2)+1
         iocc2   = iexcit2-(isppol2-1)*nexcit_pol(1)-(iunocc2-1)*nband_occ(isppol2)
         do jexcit=1,iexcit
           iexcit1 = indarr(jexcit)
           isppol1 = min((iexcit1-1)/nexcit_pol(1) +1,2)
           iunocc1 = (iexcit1-(isppol1-1)*nexcit_pol(1)-1)/nband_occ(isppol1)+1
           iocc1   = iexcit1-(isppol1-1)*nexcit_pol(1)-(iunocc1-1)*nband_occ(isppol1)

           matr(index)=omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,sing_trip)
           matr(index+1)=zero
           index=index+2
         end do

       end do

!      DEBUG
!      write(std_out,*)' after filling    matrices     '
!      ENDDEBUG

       call ZHPEV ('V','U',nexcit_win,matr,eexcit2,vec,nexcit_win,zhpev1,&
&       zhpev2,ier)

       ABI_DEALLOCATE(matr)
       ABI_DEALLOCATE(zhpev1)
       ABI_DEALLOCATE(zhpev2)
!      DEBUG
!      write(std_out,*)' after deallocating matrices     '
!      ENDDEBUG


     else

       vec(:,:,:)=zero
       do isppol=1, nsppol
         do iunocc=1,nband_k(isppol)
           do iocc=1,nband_k(isppol)
             index=iocc+nband_k(isppol)*(iunocc-1)+(isppol-1)*nexcit_pol(1)
             eexcit2(index)=omega_tddft_casida(iocc,iunocc,isppol,iocc,iunocc,isppol,sing_trip)
             vec(1,index,index)=one
           end do
         end do
       end do

     end if

!    Compute the excitation energies from the square root of eexcit2
!    eexcit(:)=sqrt(eexcit2(:)

     ABI_DEALLOCATE(eexcit)
     ABI_ALLOCATE(eexcit,(nexcit_win))

     eexcit(:)=sqrt(dabs(eexcit2(:))+tol10**2)
!    Write the first excitation energies
     if(sing_trip==1)then
       if(nsppol==1)then
         write(message, '(a,a,a,a,a,a)' )ch10,&
&         '  TDDFT singlet excitation energies (at most 20 of them are printed),',ch10,&
&         '  and corresponding total energies.                ',ch10,&
&         '  Excit#   (Ha)    and    (eV)    total energy (Ha)    major contributions '
       else
         write(message, '(a,a,a,a,a,a)' )ch10,&
&         '  TDDFT   mixed excitation energies (at most 40 of them are printed),',ch10,&
&         '  and corresponding total energies.                ',ch10,&
&         '  Excit#   (Ha)    and    (eV)    total energy (Ha)    major contributions '
       end if
     else
       write(message, '(a,a,a,a,a,a)' )ch10,&
&       '  TDDFT triplet excitation energies (at most 20 of them are printed),',ch10,&
&       '  and corresponding total energies.                ',ch10,&
&       '  Excit#   (Ha)    and    (eV)    total energy (Ha)    major contributions '
     end if
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call wrtout(std_out,' tddft : before iexcit loop',"COLL")

     do iexcit=1,min(nexcit_win,nexcitout)
       write(std_out,*)' tddft : iexcit=',iexcit
!      Select largest and next contributions
       flargest=zero ; fnext=zero
       nlargest=0 ; nnext=0
       if(nsppol==2)then
         fnext1=zero  ;  fnext2=zero
         nnext1=0  ;  nnext2=0
       end if
       do jexcit=1,nexcit_win
         ff=vec(1,jexcit,iexcit)**2+vec(2,jexcit,iexcit)**2
         if(ff>flargest+tol12)then
           if(nsppol==2)then
             nnext2=nnext1  ; fnext2=fnext1
             nnext1=nnext   ; fnext1=fnext
           end if
           nnext=nlargest ; fnext=flargest
           nlargest=indarr(jexcit) ; flargest=ff
         else if(ff>fnext+tol12)then
           if(nsppol==2)then
             nnext2=nnext1 ; fnext2=fnext1
             nnext1=nnext  ; fnext1=fnext
           end if
           nnext=indarr(jexcit) ; fnext=ff
         else if(nsppol==2)then
           if(ff>fnext1+tol12)then
             nnext2=nnext1  ; fnext2=fnext1
             nnext1=indarr(jexcit) ; fnext1=ff
           else if(ff>fnext2+tol12)then
             nnext2=indarr(jexcit) ; fnext2=ff
           end if
         end if

       end do

       isppol_l = min((nlargest-1)/nexcit_pol(1) +1,nsppol)
       iunocc_l = (nlargest-(isppol_l-1)*nexcit_pol(1)-1)/nband_occ(isppol_l)+1
       iocc_l   = nlargest-(isppol_l-1)*nexcit_pol(1)-(iunocc_l-1)*nband_occ(isppol_l)
       isppol_n = min((nnext-1)/nexcit_pol(1) +1,nsppol)
       iunocc_n = (nnext-(isppol_n-1)*nexcit_pol(1)-1)/nband_occ(isppol_n)+1
       iocc_n   = nnext-(isppol_n-1)*nexcit_pol(1)-(iunocc_n-1)*nband_occ(isppol_n)
       if(nsppol==2)then
         isppol_n1 = min((nnext1-1)/nexcit_pol(1) +1,nsppol)
         iunocc_n1 = (nnext1-(isppol_n1-1)*nexcit_pol(1)-1)/nband_occ(isppol_n1)+1
         iocc_n1   = nnext1-(isppol_n1-1)*nexcit_pol(1)-(iunocc_n1-1)*nband_occ(isppol_n1)
         isppol_n2 = min((nnext2-1)/nexcit_pol(1) +1,nsppol)
         iunocc_n2 = (nnext2-(isppol_n2-1)*nexcit_pol(1)-1)/nband_occ(isppol_n2)+1
         iocc_n2   = nnext2-(isppol_n2-1)*nexcit_pol(1)-(iunocc_n2-1)*nband_occ(isppol_n2)
       end if

       if(nsppol==1)then
         write(message,'(i4,es15.5,es14.5,es16.6,f8.2,a,i3,a,i3,a,f6.2,a,i3,a,i3,a)') &
&         iexcit,eexcit(iexcit),&
&         eexcit(iexcit)*Ha_eV,eexcit(iexcit)+etotal,&
&         flargest,'(',iocc_l,'->',iunocc_l+nband_occ(1),')',&
&         fnext,   '(',iocc_n,'->',iunocc_n+nband_occ(1),')'
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
       else
         write(chain1,'(f8.2,a,i3,a,i3,a)')flargest,'(',iocc_l,'->',iunocc_l+nband_occ(isppol_l),')'
         write(chain2,'(f8.2,a,i3,a,i3,a)')fnext,'(',iocc_n,'->',iunocc_n+nband_occ(isppol_n),')'
         if(trim(chain1)==trim(chain2))then
           write(message,'(i4,es15.5,es14.5,es16.6,a,a,a,a)') &
&           iexcit,eexcit(iexcit),&
&           eexcit(iexcit)*Ha_eV,eexcit(iexcit)+etotal,trim(chain1),'(1)',trim(chain2),'(2)'
         else
           write(message,'(i4,es15.5,es14.5,es16.6,a,a,i1,a,a,a,i1,a)') &
&           iexcit,eexcit(iexcit),&
&           eexcit(iexcit)*Ha_eV,eexcit(iexcit)+etotal,trim(chain1),'(',isppol_l,')',&
&           trim(chain2),'(',isppol_n,')'
         end if
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
         write(chain1,'(f8.2,a,i3,a,i3,a)')fnext1,'(',iocc_n1,'->',iunocc_n1+nband_occ(isppol_n1),')'
         write(chain2,'(f8.2,a,i3,a,i3,a)')fnext2,'(',iocc_n2,'->',iunocc_n2+nband_occ(isppol_n2),')'
         if(trim(chain1)==trim(chain2))then
           write(message,'(a,a,a,a,a)' ) &
&           '                                                 ',&
&           chain1,'(1)',chain2,'(2)'
         else
           write(message,'(a,a,a,i1,a,a,a,i1,a)' ) &
&           '                                                 ',&
&           chain1,'(',isppol_n1,')',&
&           chain2,'(',isppol_n2,')'
         end if
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
       end if
     end do

!    For each iexcit excitation, compute the oscillator strength (Casida, eq 47)
     write(message, '(a,a,a,a)' )ch10,&
&     '  Oscillator strengths :  (elements smaller than 1.e-6 are set to zero)',ch10,&
&     '  Excit#   (Ha)   Average    XX        YY        ZZ         XY        XZ        YZ'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

     do iexcit=1,nexcit_win

!      One of the factor of two comes from the spin degeneracy,
!      the other comes from Eq.(40) of Casida
       factor=(normint*ucvol)**2*two*two/nsppol

       osc_str(:,iexcit)=zero

!      factor=(normint*ucvol)**2*two*two/nsppol

       do jexcit=1,nexcit_win
         jexcit_cbase=indarr(jexcit)
         do idir=1,3
           osc_str(idir+1,iexcit)=osc_str(idir+1,iexcit)+ &
&           sexc(idir,jexcit)*sqrtks(jexcit_cbase)*sqrt(factor)* &
&           vec(1,jexcit,iexcit)
         end do ! idir
       end do ! jexcit

!      The "standard" definition of the oscillator strength is the square
!      of the matrix elements.
!      So, instead of the coding
!      do idir=1,3
!      osc_str(1,iexcit)=osc_str(1,iexcit)+osc_str(idir+1,iexcit)**2*third
!      enddo
!      I think that the following is more "standard"
!      Now, osc_str(2:4,iexcit) are the X, Y and Z matrix elements, not
!      yet the oscillator strengths
       osc_str(5,iexcit)=osc_str(2,iexcit)*osc_str(3,iexcit)   ! off diag XY
       osc_str(6,iexcit)=osc_str(2,iexcit)*osc_str(4,iexcit)   ! off diag XZ
       osc_str(7,iexcit)=osc_str(3,iexcit)*osc_str(4,iexcit)   ! off diag ZZ
       do idir=1,3
!        Here the X,Y, and Z matrix elements are combined to give diagonal osc. strengths
         osc_str(idir+1,iexcit)=osc_str(idir+1,iexcit)**2
         osc_str(1,iexcit)=osc_str(1,iexcit)+osc_str(idir+1,iexcit)*third ! compute the trace
       end do
!      At this stage, osc_str(1,iexcit) is exactly the same as from your coding
!      ***End of section to be checked

       do ii=1,7
         if(abs(osc_str(ii,iexcit))<tol6)osc_str(ii,iexcit)=zero
       end do
!      XG 020209 : Jean-Yves, the off-diagonal oscillator strengths
!      can become negative. It is important for automatic
!      checking that the numbers are separated by a blank, even
!      if they are negative. So replace the following format, to have at least one blank.
       write(message, '(i4,es12.5,es10.3,3es10.3,3es10.2)' )iexcit,eexcit(iexcit),osc_str(1:7,iexcit)
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end do

!    Check of the Sum rule for Casida eq.47,
!    only exact if complete basis of excitations, as well as local potentials only.
     sumx=zero
     do iexcit=1,nexcit_win
       sumx=sumx+osc_str(1,iexcit)
     end do
     write(message, '(a,es16.6)' )'  Sum of osc. strength : ',sumx
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

!    If singlet, compute Cauchy coefficients
     if(sing_trip==1.AND.nsppol==1)then
       cauchy(:)=zero
       do iexcit=1,nexcit_win
         fr_invsquare=one/(eexcit(iexcit)**2)
         fr_power=one
         do ii=1,7
           fr_power=fr_power*fr_invsquare
           cauchy(ii)=cauchy(ii)+osc_str(1,iexcit)*fr_power
         end do
       end do
       write(message, '(a,es11.3,a,es11.3,a,es11.3,a,a,es11.3,a,es11.3,a,es11.3,a,es11.3)' ) &
&       '  Cauchy coeffs (au) : ( -2)->',cauchy(1),&
&       ', ( -4)->',cauchy(2),', ( -6)->',cauchy(3),ch10,&
&       '    (-8)->',cauchy(4),', (-10)->',cauchy(5),', (-12)->',cauchy(6),', (-14)->',cauchy(7)
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if

!    End the loop on singlet or triplet
   end do

   ABI_DEALLOCATE(eexcit2)
   ABI_DEALLOCATE(vec)
   ABI_DEALLOCATE(osc_str)

!! The temporary files should be deleted at the end of this routine
   close(temp_unit2,status='delete')
   call timab(95,2,tsec)
 end if  ! end of am_master

 ABI_DEALLOCATE(omega_tddft_casida)
 ABI_DEALLOCATE(eexcit)
 ABI_DEALLOCATE(sqrtks)
 ABI_DEALLOCATE(sexc)
 ABI_DEALLOCATE(done_sexc)
 ABI_DEALLOCATE(indarr)
 ABI_DEALLOCATE(done_sexc2)

end subroutine tddft
!!***

end module m_tddft
!!***
