!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_vtorhorec
!! NAME
!!  m_vtorhorec
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (SLeroux, MMancini).
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

module m_vtorhorec

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use defs_rectypes
 use m_xmpi
 use m_pretty_rec
 use m_errors
 use m_abicore
 use m_per_cond
 use m_dtset

 use m_time,             only : timein, timab
 use m_rec,              only : Calcnrec, init_nlpsprec, cpu_distribution
 use m_rec_tools,        only : reshape_pot, trottersum, get_pt0_pt1
 use m_spacepar,         only : symrhg
 use m_fourier_interpol, only : transgrid
 use m_fft,              only : fourdp

#ifdef HAVE_GPU_CUDA
 use m_gpu_toolbox
 use m_hidecudarec
 use m_xredistribute
#endif

 implicit none

 private
!!***

 public :: vtorhorec
 public :: first_rec
!!***

contains
!!***

!!****f* ABINIT/vtorhorec
!! NAME
!! vtorhorec
!!
!! FUNCTION
!! This routine computes the new density from a fixed potential (vtrial)
!! using a recursion method
!!
!! INPUTS
!!  deltastep= if 0 the iteration step is equal to dtset%nstep
!!  initialized= if 0 the initialization of the gstate run is not yet finished
!!  operator (ground-state symmetries)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in space group
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  vtrial(nfft,nspden)=INPUT Vtrial(r).
!!  rset <type(recursion_type)> all variables for recursion
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  gprimd(3,3)=dimensional primitive translations in reciprocal space
!!
!! OUTPUT
!!  ek=kinetic energy part of total energy.
!!  enlx=nonlocal psp + potential Fock ACE part of total energy.
!!  entropy=entropy due to the occupation number smearing (if metal)
!!  e_eigenvalues=Sum of the eigenvalues - Band energy (Hartree)
!!  fermie=fermi energy (Hartree)
!!  grnl(3*natom)=stores grads of nonlocal energy wrt length scales
!!   (3x3 tensor) and grads wrt atomic coordinates (3*natom)
!!
!! SIDE EFFECTS
!!  rhog(2,nfft)=array for Fourier transform of electron density
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rset%efermi= fermi energy
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      calcnrec,cudarec,destroy_distribfft,entropyrec,fermisolverec
!!      gran_potrec,init_distribfft,nlenergyrec,prtwork,recursion,reshape_pot
!!      symrhg,timab,transgrid,wrtout,xmpi_allgatherv,xmpi_barrier,xmpi_max
!!      xmpi_sum,xredistribute
!!
!! NOTES
!!  at this time :
!!       - grnl in not implemented
!!       - symetrie usage not implemented (irrzon not used and nsym should be 1)
!!       - spin-polarized not implemented (nsppol must be 1, nspden ?)
!!       - phnons used only in symrhg
!!       - need a rectangular box (ngfft(1)=ngfft(2)=ngfft(3))
!!
!! SOURCE

subroutine vtorhorec(dtset,&
&  ek,enlx,entropy,e_eigenvalues,fermie,&
&  grnl,initialized,irrzon,nfftf,phnons,&
&  rhog, rhor, vtrial,rset,deltastep,rprimd,gprimd)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: initialized
 integer,intent(in) :: nfftf,deltastep
 real(dp),intent(out) :: e_eigenvalues,ek,enlx,entropy,fermie
 type(dataset_type),intent(in) :: dtset
 type(recursion_type),intent(inout) :: rset
!arrays
 integer, intent(in) :: irrzon(:,:,:)
 real(dp),intent(in) :: rprimd(3,3),gprimd(3,3)
 real(dp),intent(in) :: phnons(:,:,:)
 real(dp),intent(in) :: vtrial(:,:)
 real(dp),intent(inout) :: rhog(:,:)
 real(dp),intent(out) :: grnl(:)  !vz_i
 real(dp),intent(inout) :: rhor(:,:) !vz_i

!Local variables-------------------------------
!scalars
 integer :: nfftrec, dim_entro,ii1,jj1,kk1
 integer :: ierr,ii,ilmn,ipsp,dim_trott
 integer :: ipoint,ipointlocal,jj,kk,irec
 integer :: n1,n2,n3,n4,n_pt_integ_entropy
 integer :: nrec,iatom,min_pt,max_pt,swt_tm
 integer ::  get_K_S_G
 integer :: tim_fourdp,trotter
 real(dp),parameter :: perc_vmin=one
 real(dp) :: beta,drho,drhomax
 real(dp) :: entropy1,entropy2,entropy3,entropy4
 real(dp) :: entropylocal,entropylocal1,entropylocal2
 real(dp) :: entropylocal3,entropylocal4,gran_pot,gran_pot1,gran_pot2
 real(dp) :: gran_pot3,gran_pot4,gran_pot_local,gran_pot_local1
 real(dp) :: gran_pot_local2,gran_pot_local3,gran_pot_local4
 real(dp) :: inf_ucvol,intrhov,factor
 real(dp) :: nelect,potmin,rtrotter,toldrho,tolrec,tsmear
 real(dp) :: xmax,nlpotmin,ratio1,ratio2,ratio4,ratio8
 type(recparall_type) :: recpar
 character(len=500) :: msg
 !arrays
 integer :: ngfftrec(18),trasl(3)
 integer,pointer :: gcart_loc(:,:)
 integer,allocatable :: bufsize(:), bufdispl(:)
 real(dp) :: tsec(2),tsec2(2)
 real(dp) :: exppot(0:dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)-1)
 real(dp),target :: rholocal(1:rset%par%ntranche)
 real(dp),target :: alocal(0:rset%min_nrec,1:rset%par%ntranche)
 real(dp),target :: b2local(0:rset%min_nrec,1:rset%par%ntranche)
 real(dp),pointer :: rho_wrk(:)
 real(dp),pointer :: a_wrk(:,:),b2_wrk(:,:)
 real(dp),allocatable :: rholocal_f(:), rholoc_2(:)
 real(dp),allocatable :: rhogf(:,:),rhogc(:,:),aloc_copy(:,:),b2loc_copy(:,:)
 real(dp),allocatable :: gran_pot_v_f(:,:),gran_pot_v_c(:,:),gran_pot_v_2(:,:)
 real(dp),allocatable :: entropy_v_f(:,:),entropy_v_c(:,:),entropy_v_2(:,:)
 real(dp),allocatable :: ablocal_1(:,:,:),ablocal_2(:,:,:),ablocal_f(:,:,:)
 real(dp),allocatable :: exppotloc(:)
 real(dp),allocatable :: projec(:,:,:,:,:)
#if defined HAVE_GPU_CUDA
 integer :: max_rec
 integer,allocatable :: vcount_0(:), displs_0(:)
 integer,allocatable :: vcount_1(:), displs_1(:)
 real(dp),allocatable,target :: rho_hyb(:)
 real(dp),allocatable,target :: a_hyb(:,:),b2_hyb(:,:)
 real(cudap),allocatable :: an_dev(:,:)
 real(cudap),allocatable :: bn2_dev(:,:)
#endif

! *********************************************************************
 if(rset%debug)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' vtorhorec : enter '
   call wrtout(std_out,msg,'PERS')
 end if

 call timab(21,1,tsec)
 call timab(600,1,tsec2)
!##################################################################################################
!!--Initalization in the FIRST time in VTORHOREC is made in SCFCV by FIRST_REC routine
!!--Parameters for the recursion method AND  Initialisation

 trotter = dtset%recptrott  !--trotter parameter
 nelect  = dtset%nelect     !--number of electrons

 toldrho = dtset%rectolden  !--tollerance for density
 tolrec  = toldrho*1.d-2    !--tollerance for local density

 tsmear  = dtset%tsmear     !--temperature
 beta    = one/tsmear       !--inverse of temperature

 factor = real(dtset%recgratio**3*100*rset%mpi%nproc,dp)/real(nfftf,dp)

!--Assignation of the rset variable:
 nrec       = rset%min_nrec
 nfftrec    = rset%nfftrec
 ngfftrec   = rset%ngfftrec
 inf_ucvol  = rset%inf%ucvol

 min_pt = rset%par%displs(rset%mpi%me)+1
 max_pt = min_pt+rset%par%vcount(rset%mpi%me)-1

!--In the last self-constistent loop or if density is converged the
!thermodynamics quantities are calculated
 get_K_S_G = 0; if(deltastep==0 .or. rset%quitrec/=0) get_K_S_G = 1;

!--Rewriting the trotter parameter
 rtrotter  = max(half,real(trotter,dp))
 dim_trott = max(0,2*trotter-1)

!--Variables Optimisation
 ratio1 = beta/rtrotter
 ratio2 = ratio1/two
 ratio4 = ratio1/four
 ratio8 = ratio1/eight

!--Integration points for entropy
 n_pt_integ_entropy = max(25,dtset%recnpath)

!-- energies non-local: at day not implemented
 enlx = zero
 grnl = zero
!jmb
 ek = zero

!--only a copy of ngfft(1:3) and nfft  (to purge!!)
 n1 = dtset%ngfft(1) ; n2 = dtset%ngfft(2) ; n3 = dtset%ngfft(3)
 n4 = n3

!--time switch to measure gpu-cpu syncrhonisation
 swt_tm = 0 ; !no gpu initally

 exppot = zero
 nullify(gcart_loc)

 if(dtset%rectesteg==1)then
!  --Free electron gas case
   exppot = one
 else
   if(.not.(rset%nl%nlpsp)) then
!    --Local case
!    --COMPUTATION OF exp( -beta*pot/(4*rtrotter))
     ABI_ALLOCATE(gcart_loc,(0,0))
     gcart_loc = 0
     exppot = exp( -(ratio4*vtrial(:,1)))
   else
!    --Non-Local case
!    --COMPUTATION OF exp(-beta*pot/(8*rtrotter))
     exppot = exp( -(ratio8*vtrial(:,1)))

     ABI_ALLOCATE(gcart_loc,(3,dtset%natom))
     gcart_loc = rset%inf%gcart
     ABI_ALLOCATE(projec,(0:ngfftrec(1)-1,0:ngfftrec(2)-1,0:ngfftrec(3)-1,rset%nl%lmnmax,dtset%natom))
     projec = zero

     if(.not.(rset%tronc)) then
       do iatom =1, dtset%natom
         ipsp = dtset%typat(iatom)
         do ilmn = 1,rset%nl%lmnmax
           projec(:,:,:,ilmn,iatom) = reshape(rset%nl%projec(:,ilmn,ipsp),shape=shape(projec(:,:,:,1,1)))
           do ii=1,3
             projec(:,:,:,ilmn,iatom) = cshift(projec(:,:,:,ilmn,iatom),shift=ngfftrec(ii)/2-gcart_loc(ii,iatom),dim=ii)
           end do
         end do
       end do
     end if
   end if
 end if

!###################################################################################
!MAIN LOOP

 rholocal = zero; alocal = zero; b2local = zero
 ipointlocal = 1

!--Allocation: if hybrid calculation is done then I have to use
!balanced work on devices.
 nullify(rho_wrk,a_wrk,b2_wrk)

 if(rset%load == 1)then
#ifdef HAVE_GPU_CUDA
   ABI_ALLOCATE(rho_hyb,(1:rset%GPU%par%npt))
   ABI_ALLOCATE(a_hyb,(0:nrec,1:rset%GPU%par%npt))
   ABI_ALLOCATE(b2_hyb,(0:nrec,1:rset%GPU%par%npt))
   rho_hyb = zero; a_hyb = zero; b2_hyb = zero
   rho_wrk => rho_hyb
   a_wrk   => a_hyb
   b2_wrk  => b2_hyb
   recpar  = rset%GPU%par
#endif
 else
   rho_wrk => rholocal
   a_wrk   => alocal
   b2_wrk  => b2local
   recpar  = rset%par
 end if

!#if defined HAVE_GPU_CUDA
!if(rset%debug)then
!write (std_out,*) 'rset%recGPU%nptrec ',rset%recGPU%nptrec
!write (std_out,*) 'rset%gpudevice ',rset%gpudevice
!write (std_out,*) 'rset%ngfftrec ',rset%ngfftrec(1:3)
!write (std_out,*) 'rset%min_nrec ',rset%min_nrec
!write (std_out,*) 'rset%par%ntranche ',rset%par%ntranche
!write (std_out,*) 'rset%par%min_pt ',rset%par%min_pt,min_pt
!write (std_out,*) 'rset%par%max_pt ',rset%par%max_pt,max_pt
!write (std_out,*) 'pt0 ',rset%par%pt0%x,rset%par%pt0%y,rset%par%pt0%z
!write (std_out,*) 'pt1 ',rset%par%pt1%x,rset%par%pt1%y,rset%par%pt1%z
!end if
!#endif

 if(rset%gpudevice>=0) then
#if defined HAVE_GPU_CUDA
   swt_tm = 1;
   call timab(607,1,tsec2)
   ABI_ALLOCATE(an_dev,(0:recpar%npt-1,0:nrec))
   ABI_ALLOCATE(bn2_dev,(0:recpar%npt-1,0:nrec))
   an_dev = zero
   bn2_dev = zero; bn2_dev(:,0) = one

   call cudarec( rset, exppot,an_dev,bn2_dev,&
&   beta,trotter,tolrec,dtset%recgratio,dtset%ngfft(:3),max_rec)

   max_rec = min(max_rec,nrec)
   a_wrk(0:max_rec,1:recpar%npt)  = transpose(an_dev(0:,0:max_rec))
   b2_wrk(0:max_rec,1:recpar%npt) = transpose(bn2_dev(0:,0:max_rec))

   ABI_DEALLOCATE(an_dev)
   ABI_DEALLOCATE(bn2_dev)
   call timab(607,2,tsec2)

   ipointlocal = recpar%npt+1

!  !DEBUG CUDA
!  if(rset%debug)then
!  if( rset%mpi%me==0)then
!  do ipoint = 1,rset%par%npt,1
!  kk=ipoint/(product(dtset%ngfft(:2)))
!  jj=ipoint/dtset%ngfft(1)-kk*dtset%ngfft(2)
!  ii=ipoint-jj*dtset%ngfft(1)-kk*dtset%ngfft(2)*dtset%ngfft(1)
!  write(msg,'(a,4i8,2(a,a9,5f12.6))')&
!  & 'pt',ipoint,ii,jj,kk,&
!  & ch10,'an-gpu   ',real(alocal(:4,ipoint)),&
!  & ch10,'b2n-gpu  ',real(b2local(:4,ipoint))
!  call wrtout(std_out,msg,'COLL')
!  end do
!  endif
!  endif
!  !ENDDEBUG CUDA
#endif

 else
   if (.not.(rset%tronc)) then
     graou1 : do kk = recpar%pt0%z,recpar%pt1%z,dtset%recgratio
       do jj = 0,n2-1,dtset%recgratio
         do ii = 0,n1-1,dtset%recgratio
           ipoint = ii+(jj+kk*n2)*n1
!          --Local position of atoms
           if (ipoint<recpar%min_pt) cycle
!          --Computation done by that proc
           tim_fourdp=6
           call recursion(exppot,ii,jj,kk, &
&           a_wrk(:,ipointlocal), &
&           b2_wrk(:,ipointlocal), &
&           rho_wrk(ipointlocal),&
&           nrec, rset%efermi,tsmear,rtrotter,dim_trott, &
&           rset%ZT_p, &
&           tolrec,dtset%typat,rset%nl,&
&           rset%mpi,nfftrec,ngfftrec,rset%inf,&
&           tim_fourdp,dtset%natom,projec,1)
           ipointlocal = ipointlocal + 1
!          write(std_out,*)'ipointlocal',ipoint,ipointlocal,ii,jj,kk
           call prtwork(dtset%recgratio**3*ipointlocal*100*rset%mpi%nproc/nfftrec)
           if(ipoint==recpar%max_pt) exit graou1
         end do
       end do
     end do graou1
   else !--We use a troncation
     ABI_ALLOCATE(exppotloc,(0:nfftrec-1))
     graou2 : do kk = recpar%pt0%z,recpar%pt1%z,dtset%recgratio
       do jj = 0,n2-1,dtset%recgratio
         do ii = 0,n1-1,dtset%recgratio
           ipoint = ii+(jj+kk*n2)*n1
           if (ipoint<recpar%min_pt) cycle
!          computation done by that proc
           exppotloc = zero
!          --Traslation to move position on the ngfftrec grid center
           trasl = -(/ii,jj,kk/)+ngfftrec(:3)/2
           if(rset%nl%nlpsp) then
             do iatom=1,dtset%natom
!              --local position of atoms
               gcart_loc(:,iatom) = rset%inf%gcart(:,iatom)+trasl
               gcart_loc(:,iatom) = modulo(gcart_loc(:,iatom),(/n1,n2,n3/))
!              --Traslation of non-local projectors
               do ilmn = 1,rset%nl%lmnmax
                 projec(:,:,:,ilmn,iatom) = reshape(rset%nl%projec(:,ilmn,dtset%typat(iatom)),shape=shape(projec(:,:,:,1,1)))
                 do ii1=1,3
                   projec(:,:,:,ilmn,iatom) = eoshift(projec(:,:,:,ilmn,iatom),shift=ngfftrec(ii1)/2-gcart_loc(ii1,iatom),dim=ii1)
                 end do
               end do
             end do
           end if

           call reshape_pot(trasl,nfftf,nfftrec,&
&           dtset%ngfft(:3),ngfftrec(:3),&
&           exppot,exppotloc)

           tim_fourdp=6
           call recursion(exppotloc,ngfftrec(1)/2,ngfftrec(2)/2,ngfftrec(3)/2, &
&           a_wrk(:,ipointlocal), &
&           b2_wrk(:,ipointlocal), &
&           rho_wrk(ipointlocal),&
&           nrec, rset%efermi,tsmear,rtrotter,dim_trott, &
&           rset%ZT_p, &
&           tolrec,dtset%typat,rset%nl,&
&           rset%mpi,nfftrec,ngfftrec,rset%inf,&
&           tim_fourdp,dtset%natom,projec,1)
           ipointlocal = ipointlocal + 1
           call prtwork(factor*real(ipointlocal,dp))
           if(ipoint==recpar%max_pt) exit graou2
         end do
       end do
     end do graou2
     ABI_DEALLOCATE(exppotloc)
   end if

 end if
 write(msg,'( a12,i12)')'ipointlocal',ipointlocal
 call wrtout(std_out,msg,'PERS')
 call timab(613+swt_tm,1,tsec2)  !!--start time-counter: sync gpu-cpu
 call xmpi_barrier(rset%mpi%comm_bandfft)
 call timab(613+swt_tm,2,tsec2)  !!--stop time-counter: sync gpu-cpu

!!#############################################################
!!--ASSIGNATION PARAMETERS TO MENAGE PARALLELISME OF TRANSGRID
 if((rset%load==1 .or. dtset%recgratio>1))then
!  --Bufsize contains the values of the number of points calculated
!  by any proc on the coarse grid; bufsize_f on the fine grid
   call timab(604,1,tsec2) !--start time-counter: transgrid
   ABI_ALLOCATE(bufsize,(0:rset%mpi%nproc-1))
   ABI_ALLOCATE(bufdispl,(0:rset%mpi%nproc-1))
   bufsize = 0;
   bufsize(rset%mpi%me) = rset%par%npt
   call xmpi_sum(bufsize,rset%mpi%comm_bandfft,ierr)

   bufdispl(0) = 0;
   if(rset%mpi%nproc>1) bufdispl(1:) = (/(sum(bufsize(:ii)),ii=0,rset%mpi%nproc-1)/)
   call timab(604,2,tsec2) !--stop time-counter: transgrid
 end if
!!####################################################################
!!--REDISTRIBUTION OF LOAD ON PROCS AFTER RECURSION IF CUDA IS USED
#ifdef HAVE_GPU_CUDA
 if(rset%load==1)then
   call timab(604,1,tsec2) !--start time-counter: transgrid
   call xredistribute(rho_hyb,rset%GPU%par%vcount,rset%GPU%par%displs,&
&   rholocal,bufsize,bufdispl,rset%mpi%me,&
&   rset%mpi%nproc,&
&   rset%mpi%comm_bandfft,ierr)


   ABI_ALLOCATE(vcount_0,(0:rset%mpi%nproc-1))
   ABI_ALLOCATE(displs_0,(0:rset%mpi%nproc-1))
   ABI_ALLOCATE(vcount_1,(0:rset%mpi%nproc-1))
   ABI_ALLOCATE(displs_1,(0:rset%mpi%nproc-1))

   vcount_0 = 0
   vcount_0(rset%mpi%me) = rset%par%npt*(nrec+1)
   call xmpi_sum(vcount_0,rset%mpi%comm_bandfft,ierr)
   displs_0 = 0
   if(rset%mpi%nproc>1) displs_0(1:) = (/(sum(vcount_0(:ii)),ii=0,rset%mpi%nproc-1)/)


   vcount_1 = 0
   vcount_1(rset%mpi%me) = rset%GPU%par%npt*(nrec+1)
   call xmpi_sum(vcount_1,rset%mpi%comm_bandfft,ierr)
   displs_1 = 0
   if(rset%mpi%nproc>1) displs_1(1:) = (/(sum(vcount_1(:ii)),ii=0,rset%mpi%nproc-1)/)


   call xredistribute(a_hyb,vcount_1,displs_1,&
&   alocal,vcount_0,displs_0,&
&   rset%mpi%me,rset%mpi%nproc,&
&   rset%mpi%comm_bandfft,ierr)

   call xredistribute(b2_hyb,vcount_1,displs_1,&
&   b2local,vcount_0,displs_0,&
&   rset%mpi%me,rset%mpi%nproc,&
&   rset%mpi%comm_bandfft,ierr)

   nullify(rho_wrk,a_wrk,b2_wrk)
   ABI_DEALLOCATE(rho_hyb)
   ABI_DEALLOCATE(a_hyb)
   ABI_DEALLOCATE(b2_hyb)
   ABI_DEALLOCATE(vcount_0)
   ABI_DEALLOCATE(displs_0)
   ABI_DEALLOCATE(vcount_1)
   ABI_DEALLOCATE(displs_1)
   call timab(604,2,tsec2) !--start time-counter: transgrid
 end if
#endif

!#############################################################
!--TRANSGRID FOR THE DENSITY RHO AND THE COEFFICIENTS AN AND B2N
 if (dtset%recgratio>1) then
!  --variables allocation and initialisation-------
   write (msg,'(a)')' - TRANSGRID USING -----'
   call wrtout(std_out,msg,'COLL')
   call timab(604,1,tsec2) !--start time-counter: transgrid

   ABI_ALLOCATE(rholocal_f,(rset%pawfgr%nfft))
   ABI_ALLOCATE(rhogf,(2,rset%pawfgr%nfft))
   ABI_ALLOCATE(rhogc,(2,rset%pawfgr%nfftc))
   ABI_ALLOCATE(rholoc_2,(1:rset%pawfgr%nfftc))
   ABI_ALLOCATE(ablocal_2,(1:rset%pawfgr%nfftc,0:nrec,2))
   ABI_ALLOCATE(ablocal_f,(1:rset%pawfgr%nfft,0:nrec,2))
   ABI_ALLOCATE(ablocal_1,(1:rset%par%npt,0:nrec,2))

   call destroy_distribfft(rset%mpi%distribfft)
   call init_distribfft(rset%mpi%distribfft,'c',rset%mpi%nproc_fft,rset%pawfgr%ngfftc(2) ,rset%pawfgr%ngfftc(3))
   call init_distribfft(rset%mpi%distribfft,'f',rset%mpi%nproc_fft,rset%pawfgr%ngfft(2) ,rset%pawfgr%ngfft(3))

   rholocal_f = zero; ablocal_f = zero; ablocal_2 = zero

   ablocal_1(:,:,1) = transpose(alocal(:,1:rset%par%npt))
   ablocal_1(:,:,2) = transpose(b2local(:,1:rset%par%npt))

   if(get_K_S_G==1 .and. dtset%recgratio>1 ) then
     ABI_ALLOCATE(aloc_copy,(0:nrec,1:rset%par%npt))
     ABI_ALLOCATE(b2loc_copy,(0:nrec,1:rset%par%npt))
     aloc_copy = alocal(:,1:rset%par%npt)
     b2loc_copy = b2local(:,1:rset%par%npt)
   end if

   if(rset%mpi%nproc ==1) then
!    --SEQUENTIAL CASE--
     rholoc_2 = rholocal(1:rset%par%npt)
     ablocal_2 = ablocal_1
!    --Transigrid: coarse->fine
     rhogf = zero; rhogc = zero
     call transgrid(1,rset%mpi,dtset%nspden,1,0,0,1,rset%pawfgr,rhogc,rhogf,rholoc_2,rholocal_f)
     do jj1 = 1,2
       do ipoint = 0,nrec
         rhogf = zero; rhogc = zero
         call transgrid(1,rset%mpi,dtset%nspden,1,0,0,1,rset%pawfgr,rhogc,rhogf,&
&         ablocal_2(:,ipoint,jj1),ablocal_f(:,ipoint,jj1))
       end do
     end do
!    --Assignation of the interpolated results on the fine grid--
     rholocal = rholocal_f
     alocal = transpose(ablocal_f(:,:,1))
     b2local = transpose(ablocal_f(:,:,2))

   else
!    --PARALLEL CASE--
     rholoc_2 = zero
!    --Send on all procs rho,an,bn--
     call xmpi_allgatherv(rholocal(1:rset%par%npt),bufsize(rset%mpi%me),rholoc_2,&
&     bufsize,bufdispl,rset%mpi%comm_bandfft,ierr)
     do irec = 0,nrec
       call xmpi_allgatherv(ablocal_1(1:rset%par%npt,irec,1),bufsize(rset%mpi%me),&
&       ablocal_2(:,irec,1),bufsize,bufdispl,rset%mpi%comm_bandfft,ierr)
       call xmpi_allgatherv(ablocal_1(1:rset%par%npt,irec,2),bufsize(rset%mpi%me),&
&       ablocal_2(:,irec,2),bufsize,bufdispl,rset%mpi%comm_bandfft,ierr)
     end do


!    --Transigrid: coarse->fine on differents procs (with respect
!    the number of recursion)
     rhogf = zero;  rhogc = zero
     call transgrid(1,rset%mpi,dtset%nspden,1,0,0,1,rset%pawfgr,rhogc,rhogf,rholoc_2,rholocal_f)

     do irec = 0,2*(nrec+1)-1
       ii1 = modulo(irec,rset%mpi%nproc)
       jj1 = 1+modulo(irec,2)
       kk1 = floor(irec/2.)
       if(maxval(abs(ablocal_2(:,kk1,jj1))) > tol10 .and. rset%mpi%me == ii1) then
         rhogf = zero; rhogc = zero
         call transgrid(1,rset%mpi,dtset%nspden,1,0,0,1,rset%pawfgr,rhogc,&
&         rhogf,ablocal_2(:,kk1,jj1),ablocal_f(:,kk1,jj1))
       end if
     end do

!    --Recuperation of all interpolated results
!    from procs to allprocs
     call xmpi_sum(ablocal_f,rset%mpi%comm_bandfft,ierr)

!    --Assignation the interpolated results on the fine grid
!    any procs to obtain the same point as in the standard recursion
     do ii1 = 0, rset%mpi%nproc-1
       jj1 = rset%par%displs(ii1)+1
       if(ii1 == rset%mpi%me) then
         alocal = transpose(ablocal_f(jj1:jj1+rset%par%vcount(ii1)-1,:,1))
         b2local = transpose(ablocal_f(jj1:jj1+rset%par%vcount(ii1)-1,:,2))
         rholocal = rholocal_f(jj1:jj1+rset%par%vcount(ii1)-1)
       end if
     end do
   end if
   ABI_DEALLOCATE(ablocal_f)
   ABI_DEALLOCATE(ablocal_2)
   ABI_DEALLOCATE(ablocal_1)
   ABI_DEALLOCATE(rhogf)
   ABI_DEALLOCATE(rholocal_f)
   ABI_DEALLOCATE(rhogc)
   ABI_DEALLOCATE(rholoc_2)

   call timab(604,2,tsec2) !--stop time-counter: transgrid
 else
   write(msg,'(a)')' - TRANSGRID NOT USED --'
   call wrtout(std_out,msg,'COLL')
 end if
!!--End transgrid
!!###############################################################
!###################################
!--Fermi energy computation
!--find the good mu by imposing the electrons number
 call fermisolverec(rset%efermi,rholocal,alocal,b2local,rset%debug,&
& rset%min_nrec,tsmear,trotter,nelect,tol10,100, &
& rset%par%ntranche,rset%mpi,inf_ucvol,& !.False. .and.&
& (rset%tp==2 .or. rset%tp==3) .and. trotter>1)

!#################################################################
!######### ENTROPY AND GRAN POTENTIAL COMPUTATION  ##################
 entropy = zero
 gran_pot = zero
 noentropie : if(get_K_S_G==1)then
   entropy1 = zero; entropy2 = zero ;entropy3 = zero; entropy4 = zero
   gran_pot1 = zero ; gran_pot2 = zero; gran_pot3 = zero; gran_pot4 = zero

!  --Seek for the min of the path integral
   potmin = zero;  nlpotmin = zero
   if(dtset%rectesteg/=1) potmin = minval(vtrial(:,1))
   if(rset%nl%nlpsp)  nlpotmin = minval(rset%nl%eival(:,:,:))
   xmax = exp(-ratio2*(potmin+nlpotmin-rset%efermi))

   dim_entro = 0;  if(rset%debug) dim_entro = 4;

   if(dtset%recgratio>1) then
!    --Recgratio>1
     ABI_ALLOCATE(rhogf,(2,rset%pawfgr%nfft))
     ABI_ALLOCATE(rhogc,(2,rset%pawfgr%nfftc))
     ABI_ALLOCATE(entropy_v_f,(rset%pawfgr%nfft,0:4))
     ABI_ALLOCATE(entropy_v_c,(rset%pawfgr%nfftc,0:4))
     ABI_ALLOCATE(entropy_v_2,(1:rset%par%npt,0:4))
     ABI_ALLOCATE(gran_pot_v_f,(rset%pawfgr%nfft,0:4))
     ABI_ALLOCATE(gran_pot_v_c,(rset%pawfgr%nfftc,0:4))
     ABI_ALLOCATE(gran_pot_v_2,(1:rset%par%npt,0:4))

     entropy_v_c = zero; entropy_v_f = zero; entropy_v_2 = zero
     gran_pot_v_c = zero; gran_pot_v_f = zero; gran_pot_v_2 = zero


     do ipoint = 1,rset%par%npt
       call entropyrec(exp(rset%efermi*ratio2)*aloc_copy(:,ipoint), &
&       exp(rset%efermi*ratio1)*b2loc_copy(:,ipoint), &
&       nrec,trotter,entropy_v_2(ipoint,0),two,&
&       rset%debug,n_pt_integ_entropy,perc_vmin*xmax,&
&       entropy_v_2(ipoint,1),&
&       entropy_v_2(ipoint,2),&
&       entropy_v_2(ipoint,3),&
&       entropy_v_2(ipoint,4))

       call gran_potrec(exp(rset%efermi*ratio2)*aloc_copy(:,ipoint), &
&       exp(rset%efermi*ratio1)*b2loc_copy(:,ipoint), &
&       nrec,trotter,gran_pot_v_2(ipoint,0),two,&
&       rset%debug,n_pt_integ_entropy,perc_vmin*xmax,&
&       gran_pot_v_2(ipoint,1),&
&       gran_pot_v_2(ipoint,2),&
&       gran_pot_v_2(ipoint,3),&
&       gran_pot_v_2(ipoint,4))
     end do
     ABI_DEALLOCATE(aloc_copy)
     ABI_DEALLOCATE(b2loc_copy)

     call timab(613+swt_tm,1,tsec2)  !!--start time-counter: sync gpu-cpu
     call xmpi_barrier(rset%mpi%comm_bandfft)
     call timab(613+swt_tm,2,tsec2)  !!--stop time-counter: sync gpu-cpu

     call timab(604,1,tsec2) !--start time-counter: transgrid
     if(rset%mpi%nproc==1) then
       entropy_v_c = entropy_v_2
       gran_pot_v_c = gran_pot_v_2
     end if
     do ii1 = 0,dim_entro
       call xmpi_allgatherv(entropy_v_2(:,ii1),bufsize(rset%mpi%me),&
&       entropy_v_c(:,ii1),bufsize,bufdispl,&
&       rset%mpi%comm_bandfft,ierr)
       call xmpi_allgatherv(gran_pot_v_2(:,ii1),bufsize(rset%mpi%me),&
&       gran_pot_v_c(:,ii1),bufsize,bufdispl,&
&       rset%mpi%comm_bandfft,ierr)

       if(maxval(abs(entropy_v_c(:,ii1))) > tol10) then
         rhogf = zero; rhogc = zero
         call transgrid(1,rset%mpi,dtset%nspden,1,0,0,1,&
&         rset%pawfgr,rhogc,rhogf,entropy_v_c(:,ii1),entropy_v_f(:,ii1))
       end if
       if(maxval(abs(gran_pot_v_c(:,ii1))) >tol10) then
         rhogf = zero; rhogc = zero
         call transgrid(1,rset%mpi,dtset%nspden,1,0,0,1,&
&         rset%pawfgr,rhogc,rhogf,gran_pot_v_c(:,ii1),gran_pot_v_f(:,ii1))
       end if
     end do
     call timab(604,2,tsec2) !--stop time-counter: transgrid

     entropy  = sum(entropy_v_f(:,0))
     gran_pot  = sum(gran_pot_v_f(:,0))

     if(rset%debug)then
       entropy1 = sum(entropy_v_f(:,1))
       entropy2 = sum(entropy_v_f(:,2))
       entropy3 = sum(entropy_v_f(:,3))
       entropy4 = sum(entropy_v_f(:,4))
       gran_pot1 = sum(gran_pot_v_f(:,1))
       gran_pot2 = sum(gran_pot_v_f(:,2))
       gran_pot3 = sum(gran_pot_v_f(:,3))
       gran_pot4 = sum(gran_pot_v_f(:,4))
     end if

     ABI_DEALLOCATE(entropy_v_f)
     ABI_DEALLOCATE(entropy_v_c)
     ABI_DEALLOCATE(entropy_v_2)
     ABI_DEALLOCATE(rhogf)
     ABI_DEALLOCATE(rhogc)
     ABI_DEALLOCATE(gran_pot_v_f)
     ABI_DEALLOCATE(gran_pot_v_c)
     ABI_DEALLOCATE(gran_pot_v_2)


   else
!    --Recgratio=1
     do ipoint = 1,rset%par%ntranche
       call entropyrec(exp(rset%efermi*ratio2)*alocal(:,ipoint), &
&       exp(rset%efermi*ratio1)*b2local(:,ipoint), &
&       nrec,trotter,entropylocal,two,&
&       rset%debug,n_pt_integ_entropy,perc_vmin*xmax,&
&       entropylocal1,entropylocal2,&
&       entropylocal3,entropylocal4)
       call gran_potrec(exp(rset%efermi*ratio2)*alocal(:,ipoint), &
&       exp(rset%efermi*ratio1)*b2local(:,ipoint), &
&       nrec,trotter,gran_pot_local,two,& !/ucvol,&
&       rset%debug,n_pt_integ_entropy,perc_vmin*xmax,&
&       gran_pot_local1,gran_pot_local2,&
&       gran_pot_local3,gran_pot_local4)

       entropy = entropy + entropylocal
       gran_pot = gran_pot + gran_pot_local
       if(rset%debug)then
         entropy1 = entropy1 + entropylocal1
         entropy2 = entropy2 + entropylocal2
         entropy3 = entropy3 + entropylocal3
         entropy4 = entropy4 + entropylocal4
         gran_pot1 = gran_pot1 + gran_pot_local1
         gran_pot2 = gran_pot2 + gran_pot_local2
         gran_pot3 = gran_pot3 + gran_pot_local3
         gran_pot4 = gran_pot4 + gran_pot_local4
       end if

     end do

     call xmpi_sum(entropy,rset%mpi%comm_bandfft ,ierr)
     call xmpi_sum(gran_pot,rset%mpi%comm_bandfft ,ierr)
     if(rset%debug)then
       call xmpi_sum(entropy1,rset%mpi%comm_bandfft ,ierr)
       call xmpi_sum(entropy2,rset%mpi%comm_bandfft ,ierr)
       call xmpi_sum(entropy3,rset%mpi%comm_bandfft ,ierr)
       call xmpi_sum(entropy4,rset%mpi%comm_bandfft ,ierr)
       call xmpi_sum(gran_pot1,rset%mpi%comm_bandfft ,ierr)
       call xmpi_sum(gran_pot2,rset%mpi%comm_bandfft ,ierr)
       call xmpi_sum(gran_pot3,rset%mpi%comm_bandfft ,ierr)
       call xmpi_sum(gran_pot4,rset%mpi%comm_bandfft ,ierr)
     end if
   end if

   if(rset%debug)then
     write(msg,'(2(2a,4(2a,es11.4,a)))')&
&     ' --------------------------'        ,ch10, &
&     '  entropy, horiz path=',' ',entropy1,ch10, &
&     '  entropy, xmax  path=',' ',entropy2,ch10, &
&     '  entropy, xmin  path=',' ',entropy3,ch10, &
&     '  entropy, zero  path=',' ',entropy4,ch10, &
&     ' --------------------------'        ,ch10, &
&     ' -omega/T, horiz path=',' ',gran_pot1,ch10, &
&     ' -omega/T, xmax  path=',' ',gran_pot2,ch10, &
&     ' -omega/T, xmin  path=',' ',gran_pot3,ch10, &
&     ' -omega/T, zero  path=',' ',gran_pot4,ch10
     call wrtout(std_out,msg,'COLL')
   end if

   e_eigenvalues=tsmear*(entropy-gran_pot) + rset%efermi*nelect
!  --In reality gran_pot is not the gran potential but the
!  potential omega=-PV (Landau-potential or grand-potential)
!  divided by -T so the internal energy
!  U:=e_eigenvalues= TS+omega+muN = ST-T*sum(ln(1-n))+muN =
!  T(S-gran_pot)+muN


   if(rset%nl%nlpsp) then
     call nlenergyrec(rset,enlx,exppot,dtset%ngfft,dtset%natom,&
&     dtset%typat,tsmear,trotter,tolrec)
   end if
 end if noentropie
!##### END ENTROPY AND GRAN POTENTIAL COMPUTATION  ##################
!#################################################################

!if(associated(projec))
 if(rset%nl%nlpsp)  then
   ABI_DEALLOCATE(projec)
 end if
 if(associated(gcart_loc))  then
   ABI_DEALLOCATE(gcart_loc)
 end if
 if((dtset%recgratio/=1 .or. rset%load==1))  then
   ABI_DEALLOCATE(bufdispl)
   ABI_DEALLOCATE(bufsize)
 end if
!------------------------------------------------------------------
!--Check if the convergence is reached for rho
 drho = maxval(abs(rhor(min_pt:max_pt,1)-rholocal(:)))
 drhomax = drho
 call xmpi_max(drho,drhomax,rset%mpi%comm_bandfft,ierr)

!write(std_out,*)'drhomax,toldrho',drhomax,toldrho
 if(drhomax<toldrho)then
   rset%quitrec = rset%quitrec+1
 else
   rset%quitrec = 0
 end if

!-------------------------------------------------------------------
!--Density on all procs
 rhor(min_pt:max_pt,1) = rholocal(:)
 if(rset%mpi%nproc /= 1)then
   call xmpi_allgatherv(rholocal,rset%par%ntranche,rhor(:,1),&
&   rset%par%vcount,rset%par%displs,&
&   rset%mpi%comm_band,ierr)
 end if

!--------------------------------------------------------------------
!--2nd EKIN CALCULATION: this method is used
 noekin2 : if(get_K_S_G==1)then
   intrhov = (inf_ucvol)*sum(rholocal*vtrial(min_pt:max_pt,1))
   call xmpi_sum(intrhov,rset%mpi%comm_bandfft ,ierr)

   ek = e_eigenvalues-intrhov-enlx


   if(rset%debug) then
     write (msg,'(2a,3f15.10,2a,3f15.10,2a,f15.10,a)') ch10,&
&     ' ek,int(rho*V),ek+int(rho*V) ', ek, intrhov,  ek+ intrhov,ch10, &
&     ' kT*S, kT*sum(ln(...)), diff ', tsmear*entropy, tsmear*gran_pot, tsmear*(entropy-gran_pot),ch10, &
&     ' kT(S-sum(ln(...)))+mu*nelect', tsmear*(entropy-gran_pot)+rset%efermi*nelect,ch10
     call wrtout(std_out,msg,'COLL')
   end if


 end if noekin2
!--------------------------------------------------------------------
 fermie = rset%efermi

!--------------------------------------------------------
!!--At the first step to find the max number of recursion
!!  needed to convergence, then redefine nrec.
 if(initialized==0 .and. dtset%ntime>0) then
   call  Calcnrec(rset,b2local)
 end if

!--------------------------------------------------------
 call destroy_distribfft(rset%mpi%distribfft)
 call init_distribfft(rset%mpi%distribfft,'c',rset%mpi%nproc_fft,rset%ngfftrec(2),rset%ngfftrec(3))
 call init_distribfft(rset%mpi%distribfft,'f',rset%mpi%nproc_fft,dtset%ngfft(2),dtset%ngfft(3))

!--Printing results
 write(msg,'(3a,f15.10)')&
& ' -- Results: --------------------------------------',ch10,&
& ' mu          =',rset%efermi
 call wrtout(std_out,msg,'COLL')
 if(get_K_S_G==1)then
   write(msg,'(a,f15.10,6(2a,f15.10))')&
&   ' potmin      =',potmin,ch10,&
&   ' <V_eff>     =',intrhov,ch10,&
&   ' entropy     =',entropy,ch10,&
&   ' -omega/T    =',gran_pot,ch10,&
&   ' eigenvalues =',e_eigenvalues,ch10,&
&   ' kinetic     =',ek,ch10,&
&   ' non-loc ene =',enlx
   call wrtout(std_out,msg,'COLL')
 end if
 write(msg,'(a,50a)')' ',('-',ii=1,50)
 call wrtout(std_out,msg,'COLL')
!write(std_out,*)'is the pressure ',gran_pot*tsmear/(rset%inf%ucvol*real(nfftrec,dp))

!--Structured debugging : if rset%debug=T, stop here.
 if(.false.)then !(rset%debug)
   call wrtout(std_out,'  rhor ','PERS')
   write(std_out,*)rhor(:,1)
   call wrtout(std_out,' ','COLL')
   write(msg,'(a,2d10.3)')'  temps recursion    ',tsec
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a,l1,a)') ' vtorhorec : rset%debug=-',rset%debug,', debugging mode => stop '
   MSG_ERROR(msg)
 end if

 call timab(600,2,tsec2)
 call timab(21,2,tsec)

 call symrhg(1,gprimd,irrzon,rset%mpi,nfftf,&
& dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%ngfft,dtset%nspden,&
& dtset%nsppol,dtset%nsym,phnons,rhog,rhor,rprimd,dtset%symafm,dtset%symrel)

end subroutine vtorhorec
!!***

!!****f* ABINIT/entropyrec
!! NAME
!! entropyrec
!!
!! FUNCTION
!! This routine computes the local part of the entropy at a point using a path integral,
!! in the recursion method.
!!
!!  an, bn2 : coefficient given by the recursion.
!!  nrec=order of recursion
!!  trotter=trotter parameter
!!  multce=a multiplicator for computing entropy ; 2 for non-spin-polarized system
!!  debug_rec=debug variable
!!  n_pt_integ=number of points of integration for the path integral
!!  xmax =max point of integration on the real axis

!! OUTPUT
!!  ent_out=entropy at the point
!!  ent_out1,ent_out2,ent_out3,ent_out4=debug entropy at the point
!!
!! PARENTS
!!      vtorhorec
!!
!! CHILDREN
!!      timab,wrtout
!!
!! NOTES
!!  at this time :
!!       - multce should be not used
!!       - the routine should be integraly rewrited and use the routine recursion.
!!       - only modified for p /= 0
!!
!! SOURCE

subroutine entropyrec(an,bn2,nrec,trotter,ent_out,multce,debug_rec, &
&                     n_pt_integ,xmax,&
&                     ent_out1,ent_out2,ent_out3,ent_out4)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: n_pt_integ,nrec,trotter
 logical,intent(in) :: debug_rec
 real(dp), intent(in) :: multce,xmax
 real(dp),intent(out) :: ent_out,ent_out1,ent_out2,ent_out3,ent_out4
!arrays
 real(dp),intent(in) :: an(0:nrec),bn2(0:nrec)

!Local variables-------------------------------
!scalars
 integer, parameter :: level = 7
 integer, save :: first_en = 1
 integer :: ii,kk,n_pt_integ_path2,n_pt_integ_path3
 real(dp) :: arg,epsilo,step,twotrotter,xmin,dr_step
 complex(dpc) :: D,Dnew,Dold,N,Nnew,Nold,dz_path,ent_acc,ent_acc1,ent_acc2
 complex(dpc) :: ent_acc3,ent_acc4
 complex(dpc) :: funczero,z_path,zj
 complex(dpc) ::delta_calc
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp) :: iif,factor


! *************************************************************************


 call timab(610,1,tsec)

!structured debugging if debug_rec=T : print detailled result the first time we enter entropyrec

 if(debug_rec .and. first_en==1)then
   write(msg,'(a)')' '
   call wrtout(std_out,msg,'PERS')
   write(msg,'(a)')' entropyrec : enter '
   call wrtout(std_out,msg,'PERS')
   write(msg,'(a,i6)')'n_pt_integ ' , n_pt_integ
   call wrtout(std_out,msg,'COLL')
 end if

 ent_out = zero
 ent_out1 = zero
 ent_out2 = zero
 ent_out3 = zero
 ent_out4 = zero
 ent_acc = czero
 ent_acc1 = czero
 ent_acc2 = czero
 ent_acc3 = czero
 ent_acc4 = czero

!path parameters
 twotrotter = max(two*real(trotter,dp),one)
 if(trotter==0)then
   factor = tol5
   arg =pi*three_quarters
   zj = cmplx(-one,one-sin(arg),dp)
 else
   factor = xmax/ten
   arg = pi/twotrotter
   zj = cmplx( cos(arg) , sin(arg),dp )
 end if

 epsilo = factor*sin( arg )
 xmin = factor*cos( arg )
 step = (xmax-xmin)/real(n_pt_integ,dp)

!####################################################################
![xmax + i*epsilo,xmin + i*epsilo]
 dr_step = one/real(n_pt_integ,dp)
 path1:  do ii = 0,n_pt_integ
   z_path = cmplx(xmin+real(ii,dp)*(xmax-xmin)*dr_step,epsilo,dp)
   dz_path = -cmplx((xmax-xmin)*dr_step,zero,dp)

   Nold = czero
   Dold = cone
   N = cone
   D = z_path - cmplx(an(0),zero,dp)

   do kk=1,nrec
     Nnew = (z_path - cmplx(an(kk),zero,dp))*N - cmplx(bn2(kk),zero,dp)*Nold
     Dnew = (z_path - cmplx(an(kk),zero,dp))*D - cmplx(bn2(kk),zero,dp)*Dold

     Nold = N
     Dold = D
     N = Nnew
     D = Dnew

     if(kk/=nrec)then
       if((bn2(kk+1)<tol14))exit
     end if
   end do

!  <r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   delta_calc = func1_rec(z_path**twotrotter)*(N/D)*dz_path
   if(ii==0.or.ii==n_pt_integ)then
     ent_acc  = ent_acc  + half*delta_calc
     ent_acc1 = ent_acc1 + half*delta_calc
   else
     ent_acc  = ent_acc  + delta_calc
     ent_acc1 = ent_acc1 + delta_calc
   end if
 end do path1


!####################################################################
![1/2zj,0]
 if(epsilo/step>100.d0)then
   n_pt_integ_path2 = int((factor*abs(zj))/step)+1
 else
   n_pt_integ_path2 = 100
 end if

 if(trotter/=0)then
   n_pt_integ_path3 = 0
   dr_step = one/real(n_pt_integ_path2,dp)
   dz_path = -cmplx(xmin,epsilo,dp)*dr_step
   path5:  do ii = 0,n_pt_integ_path2
     z_path = cmplx(real(ii,dp)*xmin,real(ii,dp)*epsilo,dp)*dr_step
     if(abs(z_path)>tol14)then
       Nold = czero
       Dold = cone
       N = cone
       D = z_path - cmplx(an(0),zero,dp)
       do kk=1,nrec
         Nnew = (z_path - cmplx(an(kk),zero,dp))*N - cmplx(bn2(kk),zero,dp)*Nold
         Dnew = (z_path - cmplx(an(kk),zero,dp))*D - cmplx(bn2(kk),zero,dp)*Dold
         Nold = N
         Dold = D
         N = Nnew
         D = Dnew
         if(kk/=nrec)then
           if((bn2(kk+1)<tol14))exit
         end if
       end do

!      <r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
       if(abs(z_path)**twotrotter>tiny(one)) then
         funczero = func1_rec(z_path**twotrotter)
       else
         funczero = czero
       end if
       delta_calc = funczero*N/D*dz_path
       if(ii==0.or.ii==n_pt_integ_path2)then
         ent_acc  = ent_acc  + half*delta_calc
         if(debug_rec) ent_acc3 = ent_acc3 + half*delta_calc
       else
         ent_acc  = ent_acc  + funczero*delta_calc
         if(debug_rec) ent_acc3 = ent_acc3 + funczero*delta_calc
       end if
     end if
   end do path5

 else  ! trotter==0

   n_pt_integ_path3 = max(100,int((epsilo*half*pi)/real(step,dp))+1)
   dr_step = one/real(n_pt_integ_path3,dp)
   path6:  do ii = 0,n_pt_integ_path3
     iif=half*pi*real(ii,dp)*dr_step
     z_path = epsilo*cmplx(-cos(iif),1-sin(iif),dp)
     dz_path = epsilo*cmplx(sin(iif),-cos(iif),dp)*half*pi*dr_step
     if(abs(z_path)**twotrotter>tol14)then
       Nold = czero
       Dold = cone
       N = cone
       D = z_path - cmplx(an(0),zero,dp)
       do kk=1,nrec
         Nnew = (z_path - cmplx(an(kk),zero,dp))*N - cmplx(bn2(kk),zero,dp)*Nold
         Dnew = (z_path - cmplx(an(kk),zero,dp))*D - cmplx(bn2(kk),zero,dp)*Dold
         Nold = N
         Dold = D
         N = Nnew
         D = Dnew
         if(kk/=nrec .and. bn2(kk+1)<tol14) exit !-EXIT
       end do

!      <r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
       delta_calc = func1_rec(z_path**twotrotter) * N/D * dz_path
       if(ii==0.or.ii==n_pt_integ_path3)then
         ent_acc  = ent_acc + half*delta_calc
         if(debug_rec) ent_acc3 = ent_acc3 + half*delta_calc
       else
         ent_acc  = ent_acc + delta_calc    !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
         if(debug_rec) ent_acc3 = ent_acc3 + delta_calc  !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
       end if
     end if
   end do path6

 end if

 if(first_en==1 .and. debug_rec) then
   write(msg,'(a,i5,2a,i5,2a,i5,2a,es11.4,2a,es11.4,2a,es11.4)')&
&   'n_pt_path  =',n_pt_integ,ch10,&
&   'n_pt_path2 =',n_pt_integ_path2,ch10,&
&   'n_pt_path3 =',n_pt_integ_path3,ch10,&
&   'xmin       =',xmin,ch10,&
&   'xmax       =',xmax,ch10,&
&   'epsilon    =',epsilo
   call wrtout(std_out,msg,'COLL')
   first_en = 0
 end if

!####################################################################
![xmax,xmax+i*epsilo]
 dr_step = one/real(n_pt_integ_path2,dp)
 dz_path = cmplx(zero,epsilo*dr_step,dp)
 path4:  do ii = 0,n_pt_integ_path2
   z_path = cmplx(xmax,real(ii,dp)*epsilo*dr_step,dp)

   Nold = czero
   Dold = cone
   N = cone
   D = z_path - cmplx(an(0),zero,dp)

   do kk=1,nrec
     Nnew = (z_path - cmplx(an(kk),zero,dp))*N - cmplx(bn2(kk),zero,dp)*Nold
     Dnew = (z_path - cmplx(an(kk),zero,dp))*D - cmplx(bn2(kk),zero,dp)*Dold

     Nold = N
     Dold = D
     N = Nnew
     D = Dnew

     if(kk/=nrec)then
       if((bn2(kk+1)<tol14))exit
     end if
   end do

!  <r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   delta_calc = func1_rec(z_path**twotrotter)*N/D*dz_path
   if(ii==0.or.ii==n_pt_integ_path2)then

     ent_acc =  ent_acc  + half*delta_calc
     if(debug_rec) ent_acc2 = ent_acc2 + half*delta_calc
   else
     ent_acc  = ent_acc  + delta_calc
     if(debug_rec) ent_acc2 = ent_acc2 + delta_calc
   end if
 end do path4


 ent_out  = multce*real(ent_acc*cmplx(zero,-piinv,dp),dp)
 if(debug_rec) then
   ent_out1 = multce*real(ent_acc1*cmplx(zero,-piinv,dp),dp)
   ent_out2 = multce*real(ent_acc2*cmplx(zero,-piinv,dp),dp)
   ent_out3 = multce*real(ent_acc3*cmplx(zero,-piinv,dp),dp)
   ent_out4 = multce*real(ent_acc4*cmplx(zero,-piinv,dp),dp)
 end if

 call timab(610,2,tsec)

 contains

!function to integrate over the path
!func1_rec(z_path,twotrotter) =  ( z_path**twotrotter/(1+z_path**twotrotter)*log(1+1/z_path**twotrotter)+&    !- f*ln(f)
!&1/(1+z_path**twotrotter)*log(1+z_path**twotrotter))       !- (1-f)*ln(1-f)

!func1_rec(z_path_pow) =   z_path_pow/(cone+z_path_pow)*log(cone+cone/z_path_pow)+&    !- f*ln(f)
!&cone/(cone+z_path_pow)*log(cone+z_path_pow)       !- (1-f)*ln(1-f)

!other expression of func for a path like ro(t)*exp(2*i*pi/(2*p)*(j+1/2))

   function func1_rec(z)

   complex(dpc) :: func1_rec
   complex(dpc),intent(in) :: z

   func1_rec =   z/(cone+z)*log(cone+cone/z)+ cone/(cone+z)*log(cone+z)

 end function func1_rec

end subroutine entropyrec
!!***

!!****f* ABINIT/fermisolverec
!! NAME
!! fermisolverec
!!
!! FUNCTION
!! This routine computes the fermi energy in order to have a given number of
!! valence electrons in the recursion method, using a Ridder s Method
!!
!! INPUTS
!!  debug_rec=debugging variable
!!  nb_rec=order of recursion
!!  nb_point=number of discretization point in one dimension (=n1=n2=n3)
!!  temperature=temperature (Hartree)
!!  trotter=trotter parameter
!!  nelect=number of valence electrons (dtset%nelect)
!!  acc=accuracy for the fermi energy
!!  max_it=maximum number of iteration for the Ridder's Method
!!  long_tranche=number of point computed by thi proc
!!  mpi_enreg=information about MPI parallelization
!!  inf_ucvol=infinitesimal unit cell volume
!!  gputopo=true if topology gpu-cpu= 2 or 3
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  fermie=fermi energy
!!  rho=density, recomputed for the new fermi energy
!!  a, b2 : coefficient given by recursion recomputed for the new fermi energy
!!
!! PARENTS
!!      vtorhorec
!!
!! CHILDREN
!!      alloc_dens_cuda,dealloc_dens_cuda,density_cuda,density_rec,timab,wrtout
!!      xmpi_barrier,xmpi_sum
!!
!! NOTES
!!  at this time :
!!
!! SOURCE

subroutine fermisolverec(fermie,rho,a,b2,debug_rec,nb_rec, &
  &                      temperature,trotter,nelect, &
  &                      acc, max_it, &
  &                      long_tranche,mpi_enreg,&
  &                      inf_ucvol,gputopo)

!Arguments -------------------------------
 !scalars
 integer,intent(in) :: long_tranche,max_it,nb_rec,trotter
 logical,intent(in) :: debug_rec,gputopo
 real(dp),intent(in) :: acc,inf_ucvol,nelect,temperature
 real(dp), intent(inout) :: fermie
 type(MPI_type),intent(in) :: mpi_enreg
 !arrays
 real(dp), intent(inout) :: a(0:nb_rec,long_tranche), b2(0:nb_rec,long_tranche)
 real(dp), intent(inout) :: rho(long_tranche)

!Local variables-------------------------------
 !scalars
 integer  ::  ierr,ii,ipointlocal,nn,dim_trott
 real(dp) :: beta,fermieh,fermiel,fermiem,fermienew,nelecth,nelectl,nelectm
 real(dp) :: nelectnew,res_nelecth,res_nelectl,res_nelectm,res_nelectnew
 real(dp) :: rtrotter,ss,fermitol
 character(len=500) :: msg
 !arrays
 real(dp) :: tsec(2)
 real(dp) :: rhotry(long_tranche)
 !no_abirules
#ifdef HAVE_GPU_CUDA
 integer :: swt_tm,npitch
 real(cudap) :: rhocu(long_tranche)
 real(dp) :: tsec2(2)
#endif

! *************************************************************************

#ifdef HAVE_GPU_CUDA
 swt_tm = 0
#endif

 call timab(609,1,tsec)

 beta = one/temperature
 rtrotter  = max(half,real(trotter,dp))
 dim_trott = max(0,2*trotter-1)

 write(msg,'(a)')' -- fermisolverec ---------------------------------'
 call wrtout(std_out,msg,'COLL')
 if(debug_rec) then
   write (msg,'(a,d10.3)')' nelect= ',nelect
   call wrtout(std_out,msg,'COLL')
 end if
!initialisation of fermiel
 fermiel = fermie
 call timab(609,2,tsec)

!initialisation fermitol
 fermitol = acc
#ifdef HAVE_GPU_CUDA_SP
 if(gputopo)  fermitol = 1.d-3
#endif

 if(gputopo) then
#ifdef HAVE_GPU_CUDA
   swt_tm = 1
!  allocate array an and bn2 on gpu for computation of trotter formula
   call alloc_dens_cuda(long_tranche,nb_rec,dim_trott,npitch,&
&   real(a,cudap),real(b2,cudap))

   call timab(617,1,tsec)
   call density_cuda(npitch,long_tranche,nb_rec,dim_trott,&
&   real(fermiel,cudap),real(temperature,cudap),&
&   real(rtrotter,cudap),real(inf_ucvol,cudap),&
&   real(tol14,cudap),&
&   rhocu)
   rhotry = real(rhocu,dp)
   call timab(617,2,tsec)
#endif
 else
   do ipointlocal = 1,long_tranche
     call density_rec(a(:,ipointlocal),&
&     b2(:,ipointlocal),&
&     rhotry(ipointlocal),&
&     nb_rec,fermiel,temperature,rtrotter,dim_trott, &
&     tol14,inf_ucvol)
   end do
 end if

 call timab(609,1,tsec)
 nelectl = sum(rhotry)
 call xmpi_sum( nelectl,mpi_enreg%comm_bandfft,ierr)
 res_nelectl = inf_ucvol*nelectl - nelect

 if (res_nelectl /= zero) then
!  initialisation of fermih
!  excess of electrons -> smaller fermi
   res_nelecth = zero
   ii = 1
   fermieh = fermie - ten*sign(one,res_nelectl)*temperature
   do while(ii<6 .and. res_nelecth*res_nelectl>=0)
     fermieh = fermieh - ten*sign(one,res_nelectl)*temperature
     call timab(609,2,tsec)

     if(gputopo) then
#ifdef HAVE_GPU_CUDA
       call timab(617,1,tsec)
       call density_cuda(npitch,long_tranche,nb_rec,dim_trott,&
&       real(fermieh,cudap),real(temperature,cudap),&
&       real(rtrotter,cudap),real(inf_ucvol,cudap),&
&       real(tol14,cudap),&
&       rhocu)
       rhotry = real(rhocu,dp)
       call timab(617,2,tsec)
#endif
     else
       do ipointlocal = 1,long_tranche
         call density_rec(a(:,ipointlocal),  &
&         b2(:,ipointlocal), &
&         rhotry(ipointlocal), &
&         nb_rec,fermieh,temperature,rtrotter,dim_trott, &
&         tol14,inf_ucvol)
       end do
     end if
     call timab(609,1,tsec)
     nelecth = sum(rhotry)
     call xmpi_sum( nelecth,mpi_enreg%comm_bandfft ,ierr);
     res_nelecth = inf_ucvol*nelecth - nelect

     if(debug_rec) then
       write (msg,'(a,es11.4e2,a,es11.4e2)') ' Fermi energy interval',fermieh,' ',fermiel
       call wrtout(std_out,msg,'COLL')
     end if
     ii = ii +1
   end do

   if (res_nelecth*res_nelectl>0) then
     write (msg,'(4a)')' fermisolverec : ERROR- ',ch10,&
&     ' initial guess for fermi energy doesnt permit to  find solutions in solver',ch10
     MSG_ERROR(msg)
   end if

!  MAIN LOOP   ------------------------------------------------------
   main : do nn=1,max_it
!    fermiem computation
     fermiem = 0.5d0*(fermiel+fermieh)

!    nelectm = zero
     call timab(609,2,tsec)

     if(gputopo) then
#ifdef HAVE_GPU_CUDA
       call timab(617,1,tsec)
       call density_cuda(npitch,long_tranche,nb_rec,dim_trott,&
&       real(fermiem,cudap),real(temperature,cudap),&
&       real(rtrotter,cudap),real(inf_ucvol,cudap),&
&       real(tol14,cudap),&
&       rhocu)
       rhotry = real(rhocu,dp)
       call timab(617,2,tsec)
#endif
     else
       do ipointlocal = 1,long_tranche
         call density_rec(a(:,ipointlocal),  &
&         b2(:,ipointlocal), &
&         rhotry(ipointlocal), &
&         nb_rec,fermiem,temperature,rtrotter,dim_trott, &
&         tol14,inf_ucvol)
       end do
     end if

     call timab(609,1,tsec)
     nelectm = sum(rhotry)
     call xmpi_sum( nelectm,mpi_enreg%comm_bandfft,ierr)
     res_nelectm = inf_ucvol*nelectm - nelect

!    new guess
     ss = sqrt(res_nelectm**two-res_nelectl*res_nelecth)
     fermienew = fermiem + (fermiem-fermiel)*sign(one, res_nelectl-res_nelecth)*res_nelectm/ss

     call timab(609,2,tsec)
     if(gputopo) then
#ifdef HAVE_GPU_CUDA
       call timab(617,1,tsec)
       call density_cuda(npitch,long_tranche,nb_rec,dim_trott,&
&       real(fermienew,cudap),real(temperature,cudap),&
&       real(rtrotter,cudap),real(inf_ucvol,cudap),&
&       real(tol14,cudap),&
&       rhocu)
       rhotry = real(rhocu,dp)
       call timab(617,2,tsec)
#endif
     else
       do ipointlocal = 1,long_tranche
         call density_rec(a(:,ipointlocal),  &
&         b2(:,ipointlocal), &
&         rhotry(ipointlocal), &
&         nb_rec,fermienew,temperature,rtrotter,dim_trott, &
&         tol14,inf_ucvol)
       end do
     end if

     call timab(609,1,tsec)
     nelectnew = sum(rhotry)
     call xmpi_sum( nelectnew,mpi_enreg%comm_bandfft ,ierr);
     res_nelectnew = inf_ucvol*nelectnew - nelect

!    fermiel et fermieh for new iteration
     if (sign(res_nelectm,res_nelectnew) /= res_nelectm) then
       fermiel = fermiem
       res_nelectl = res_nelectm
       fermieh = fermienew
       res_nelecth = res_nelectnew
     else if (sign(res_nelectl,res_nelectnew) /= res_nelectl) then
       fermieh = fermienew
       res_nelecth = res_nelectnew
     else if (sign(res_nelecth,res_nelectnew) /= res_nelecth) then
       fermiel = fermienew
       res_nelectl = res_nelectnew
     end if

!    are we within the tolerance ?
     if ((abs(res_nelectnew) < fermitol).or.(nn == max_it)) then
       fermie = fermienew
       rho = rhotry
       if(debug_rec) then
         write (msg,'(a,es11.4e2,a,i4)')' err, num_iter ', res_nelectnew, ' ',nn
         call wrtout(std_out,msg,'COLL')
         write(msg,'(a,50a)')' ',('-',ii=1,50)
         call wrtout(std_out,msg,'COLL')
       end if
       exit main
     end if

   end do main

 end if

#ifdef HAVE_GPU_CUDA
!deallocate array on GPU
 if(gputopo) then
   call dealloc_dens_cuda()
 end if
 call timab(613+swt_tm,1,tsec2)  !!--start time-counter: sync gpu-cpu
 call xmpi_barrier(mpi_enreg%comm_bandfft)
 call timab(613+swt_tm,2,tsec2)  !!--stop time-counter: sync gpu-cpu
#endif

 call timab(609,2,tsec)
end subroutine fermisolverec
!!***

!!****f* ABINIT/density_rec
!! NAME
!! density_rec
!!
!! FUNCTION
!! This routine computes the density using  the coefficients corresponding to
!! continued fraction at a point from a fixed potential.
!!
!! INPUTS
!!  coordx, coordy, coordz=coordonnees of the computed point
!!  an, bn2 : coefficient given by density_rec. Input if get_rec_coef=0, output else
!!  nrec=order of density_rec
!!  fermie=fermi energy (Hartree)
!!  tsmear=temperature (Hartree)
!!  rtrotter=real trotter parameter
!!  tol=tolerance criteria for stopping density_rec
!!  inf_ucvol=infinitesimal unit cell volume
!!  dim_trott = max(0,2*trotter-1)
!!
!! OUTPUT
!!  rho_out=result of the continued fraction multiplied by a multiplicator
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      fermisolverec
!!
!! CHILDREN
!!      timab,trottersum
!!
!! NOTES
!!  at this time :
!!       - exppot should be replaced by ?
!!       - coord should be replaced by ?
!!       - need a rectangular box (rmet diagonal matrix)
!!
!! SOURCE

subroutine density_rec(an,bn2,rho_out,nrec, &
&                     fermie,tsmear,rtrotter, &
&                     dim_trott,tol,inf_ucvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nrec
 integer,intent(in) :: dim_trott
 real(dp),intent(in) :: fermie,tol,tsmear,inf_ucvol,rtrotter
 real(dp), intent(out) :: rho_out
!arrays
 real(dp),intent(in) :: an(0:nrec),bn2(0:nrec)
!Local variables-------------------------------
!not used, debugging purpose only
!for debugging purpose, detailled printing only once for density and ekin
!scalars
 integer, parameter :: minrec = 3
 integer  :: irec
 real(dp) :: beta,mult,prod_b2,error,errold
 real(dp) :: pi_on_rtrotter,twortrotter,exp1,exp2
 complex(dpc) :: cinv2rtrotter,coeef_mu,facrec0
! character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 complex(dpc) :: acc_rho(0:nrec)
 complex(dpc) :: D(0:dim_trott),Dold(0:dim_trott)
 complex(dpc) :: N(0:dim_trott),Nold(0:dim_trott)
!**************************************************************************

 call timab(605,1,tsec)

!##############################################################
!--Initialisation of metrics
 mult = two/inf_ucvol   !non-spined system
 beta = one/tsmear

!--Variables for optimisation
 pi_on_rtrotter = pi/rtrotter
 twortrotter = two*rtrotter
 exp1 = exp((beta*fermie)/(rtrotter))
 exp2 = exp(beta*fermie/(twortrotter))
 cinv2rtrotter = cmplx(one/twortrotter,zero,dp)
 coeef_mu = cmplx(one/exp2,zero,dp)

 N = czero;  D = cone
 facrec0 = cone
 Nold = czero; Dold = czero
!--Initialisation of accumulated density
 acc_rho = czero
!--Initialisation of estimated error
 prod_b2 = twortrotter/exp1
 errold = zero


!##############################################################
!--Main loop
 maindo : do irec = 0, nrec

!  ######################################################
!  --Density computation
!  !--using the property that: sum_i(bi*c)^2|(z-ai*c)=1/c*sum_i(bi)^2|(z/c-ai)
!  !and for c =exp(-beta*fermie/(two*rtrotter)

   call trottersum(dim_trott,error,prod_b2,pi_on_rtrotter,&
&   facrec0,coeef_mu,exp1,&
&   an(irec),bn2(irec),&
&   N,D,Nold,Dold)

   if(irec/=nrec .and. irec>=minrec)then
     if((bn2(irec+1)<tol14).or.(mult*error<tol.and.errold<tol)) exit maindo
   end if
   errold = mult*error
 end do maindo
!--Accumulated density
 rho_out = mult*real(cone-sum(N/D)*cinv2rtrotter,dp)

 call timab(605,2,tsec)

 end subroutine density_rec
!!***

!!****f* ABINIT/gran_potrec
!! NAME
!! gran_potrec
!!
!! FUNCTION
!! This routine computes the local part of the grand-potential at a point using a path integral,
!! in the recursion method.
!!
!! INPUTS
!!  an, bn2 : coefficient given by the recursion.
!!  nrec=order of recursion
!!  trotter=trotter parameter
!!  mult=a multiplicator for computing grand-potential (2 for non-spin-polarized system)
!!  debug_rec=debugging variable
!!  n_pt_integ=points for computation of path integral
!!  xmax= maximum point on the x-axis for integration
!!
!! OUTPUT
!!  ene_out=grand-potential at the point
!!  if debug_rec=T then ene_out1,ene_out2,ene_out3,ene_out4 are
!!  the different path branch contriubutions to the grand-potential.
!!  In reality it is not the gren potential but the
!!  grand-potential (omega=-PV) divided by -T
!!
!! PARENTS
!!      vtorhorec
!!
!! CHILDREN
!!      timab,wrtout
!!
!! NOTES
!!  in reality it is not the gren potential but the grand-potential (omega=-PV) divided by -T
!!  at this time :
!!       - mult should be not used
!!       - the routine should be integraly rewrited and use the routine recursion.
!!       - only modified for p /= 0
!!
!! SOURCE

subroutine gran_potrec(an,bn2,nrec,trotter,ene_out, mult, &
&                     debug_rec,n_pt_integ,xmax,&
&                     ene_out1,ene_out2,ene_out3,ene_out4)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: n_pt_integ,nrec,trotter
 logical,intent(in) :: debug_rec
 real(dp), intent(in) :: mult,xmax
 real(dp),intent(inout) :: ene_out,ene_out1,ene_out2,ene_out3,ene_out4 !vz_i
!arrays
 real(dp), intent(in) :: an(0:nrec),bn2(0:nrec)

!Local variables-------------------------------
!scalars
 integer, parameter :: level = 7
 integer, save :: first = 1
 integer :: ii,kk,n_pt_integ_path2
 real(dp) :: epsilon,step,twotrotter,xmin,dr_step
 complex(dpc) :: D,Dnew,Dold,N,Nnew,Nold,dz_path,ene_acc,ene_acc1,ene_acc2
 complex(dpc) :: ene_acc3,ene_acc4
 complex(dpc) :: z_path,delta_calc
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)

! *************************************************************************


 call timab(611,1,tsec)

!structured debugging if debug_rec=T : print detailled result the first time we enter gran_potrec
 if(debug_rec .and. first==1)then
   write(message,'(a)')' '
   call wrtout(std_out,message,'PERS')
   write(message,'(a)')' gran_potrec : enter '
   call wrtout(std_out,message,'PERS')
   write(message,'(a,i8)')'n_pt_integ ' , n_pt_integ
   call wrtout(std_out,message,'COLL')
   first=0
 end if

 ene_out = zero
 ene_acc = czero
 ene_acc1 = czero
 ene_acc2 = czero
 ene_acc3 = czero
 ene_acc4 = czero


!path parameters
!n_pt_integ = 2500
 xmin = -half
 step = (xmax-xmin)/real(n_pt_integ,dp)
 if(trotter==0)then
   twotrotter = one
   epsilon = .5d-1
 else
   twotrotter = two*real(trotter,dp)
   epsilon = half*sin( pi/twotrotter)
 end if

!xmin = -abs(xmin)**(1.d0/twotrotter)

!####################################################################
![xmax + i*epsilon,xmin + i*epsilon]
 dr_step = one/real(n_pt_integ,dp)
 dz_path = -cmplx((xmax-xmin)*dr_step,zero,dp)
 path1:  do ii = 0,n_pt_integ
!  z_path = cmplx(xmin + real(ii,dp)*(xmax-xmin)*dr_step,epsilon,dp)
   z_path = cmplx(xmin,epsilon,dp) - real(ii,dp)*dz_path
   Nold = czero
   Dold = cone
   N = cone
   D = z_path - cmplx(an(0),zero,dp)

   do kk=1,nrec
     Nnew = (z_path - cmplx(an(kk),zero,dp))*N - cmplx(bn2(kk),zero,dp)*Nold
     Dnew = (z_path - cmplx(an(kk),zero,dp))*D - cmplx(bn2(kk),zero,dp)*Dold

     Nold = N
     Dold = D
     N = Nnew
     D = Dnew

     if(kk/=nrec)then
       if((bn2(kk+1)<tol14))exit
     end if

   end do

!  <r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   delta_calc = func_rec(z_path,twotrotter)* N/D *dz_path
   if(ii==0.or.ii==n_pt_integ)then
     ene_acc = ene_acc + half*delta_calc
     if(debug_rec)  ene_acc1 = ene_acc1 + half*delta_calc
   else
     ene_acc = ene_acc + delta_calc
     if(debug_rec)  ene_acc1 = ene_acc1 + delta_calc
   end if
 end do path1

!####################################################################
![xmin + i*epsilon,xmin]
 if(epsilon/step>4.d0)then
   n_pt_integ_path2 = int(epsilon/step)+1
 else
   n_pt_integ_path2 = 5
 end if
 n_pt_integ_path2 = n_pt_integ
 dr_step = one/real(n_pt_integ_path2,dp)
 dz_path = -cmplx(zero,epsilon*dr_step,dp)
 path2:  do ii = 0,n_pt_integ_path2
!  z_path = cmplx(xmin,real(ii,dp)*epsilon*dr_step,dp)
   z_path = cmplx(xmin,zero,dp)-dz_path*real(ii,dp)
   Nold = czero
   Dold = cone
   N = cone
   D = z_path - cmplx(an(0),zero,dp)

   do kk=1,nrec
     Nnew = (z_path - cmplx(an(kk),zero,dp))*N - cmplx(bn2(kk),zero,dp)*Nold
     Dnew = (z_path - cmplx(an(kk),zero,dp))*D - cmplx(bn2(kk),zero,dp)*Dold

     Nold = N
     Dold = D
     N = Nnew
     D = Dnew

     if(kk/=nrec)then
       if((bn2(kk+1)<tol14))exit
     end if

   end do

!  <r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   delta_calc = func_rec(z_path,twotrotter)* N/D *dz_path
   if(ii==0.or.ii==n_pt_integ_path2)then
     ene_acc = ene_acc + half*delta_calc
     if(debug_rec) ene_acc3 = ene_acc3 + half*delta_calc
   else
     ene_acc = ene_acc + delta_calc
     if(debug_rec) ene_acc3 = ene_acc3 + delta_calc
   end if
 end do path2



!####################################################################
![xmin,0]
 if(xmin/=czero)then
   dr_step = one/real(n_pt_integ,dp)
   dz_path = cmplx(xmin*dr_step,zero,dp)
   path3:  do ii = 1,n_pt_integ !the integrand is 0 at 0
!    z_path = cmplx(real(ii,dp)*xmin*dr_step,zero,dp)
     z_path = real(ii,dp)*dz_path

     Nold = czero
     Dold = cone
     N = cone
     D = z_path - cmplx(an(0),zero,dp)

     do kk=1,nrec
       Nnew = (z_path - cmplx(an(kk),zero,dp))*N - cmplx(bn2(kk),zero,dp)*Nold
       Dnew = (z_path - cmplx(an(kk),zero,dp))*D - cmplx(bn2(kk),zero,dp)*Dold

       Nold = N
       Dold = D
       N = Nnew
       D = Dnew

       if(kk/=nrec)then
         if((bn2(kk+1)<tol14))exit
       end if
     end do

!    <r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
     delta_calc = func_rec(z_path,twotrotter) * N/D *dz_path
     if(ii==n_pt_integ)then
       ene_acc = ene_acc +half*delta_calc
       if(debug_rec) ene_acc4 = ene_acc4 + half*delta_calc
     else
       ene_acc = ene_acc + delta_calc
       if(debug_rec) ene_acc4 = ene_acc4 +delta_calc
     end if
   end do path3
 end if

!####################################################################
![xmax,xmax+i*epsilon]
 dr_step = one/real(n_pt_integ_path2,dp)
 dz_path = cmplx(zero,epsilon*dr_step,dp)
 path4:  do ii = 0,n_pt_integ_path2
!  z_path = cmplx(xmax,real(ii,dp)*epsilon*dr_step,dp)
   z_path = cmplx(xmax,0,dp)+real(ii,dp)*dz_path

   Nold = czero
   Dold = cone
   N = cone
   D = z_path - cmplx(an(0),zero,dp)

   do kk=1,nrec
     Nnew = (z_path - cmplx(an(kk),zero,dp))*N - cmplx(bn2(kk),zero,dp)*Nold
     Dnew = (z_path - cmplx(an(kk),zero,dp))*D - cmplx(bn2(kk),zero,dp)*Dold

     Nold = N
     Dold = D
     N = Nnew
     D = Dnew

     if(kk/=nrec)then
       if((bn2(kk+1)<tol14))exit
     end if

   end do

!  <r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   delta_calc = func_rec(z_path,twotrotter) * N/D *dz_path
   if(ii==0.or.ii==n_pt_integ_path2)then
     ene_acc = ene_acc + half*delta_calc
     if(debug_rec) ene_acc2 = ene_acc2 + half*delta_calc
   else
     ene_acc = ene_acc + delta_calc
     if(debug_rec) ene_acc2 = ene_acc2 + delta_calc
   end if
 end do path4

 ene_out = mult*real(ene_acc*cmplx(zero,-piinv,dp),dp)
 if(debug_rec) then
   ene_out1 = mult*real(ene_acc1*cmplx(zero,-piinv,dp),dp)
   ene_out2 = mult*real(ene_acc2*cmplx(zero,-piinv,dp),dp)
   ene_out3 = mult*real(ene_acc3*cmplx(zero,-piinv,dp),dp)
   ene_out4 = mult*real(ene_acc4*cmplx(zero,-piinv,dp),dp)
 end if

 call timab(611,2,tsec)

 contains

!func_rec(z_path,twotrotter) = log(cone+z_path**twotrotter)

   function func_rec(z,x)

   complex(dpc) :: func_rec
   complex(dpc),intent(in) :: z
   real(dp),intent(in) :: x

   func_rec = log(cone+z**x)

 end function func_rec

end subroutine gran_potrec
!!***

!!****f* ABINIT/nlenergyrec
!! NAME
!! nlenergyrec
!!
!! FUNCTION
!! During recursion, it computes the non-local energy
!!
!! INPUTS
!!  rset<recursion_type>=contains all recursion parameters
!!  exppot=exponential of -1/tsmear*vtrial (computed only once in vtorhorec)
!!  tsmear=temperature (Hartree)
!!  trotter=trotter parameter
!!  tol=tolerance criteria for stopping recursion_nl
!!  ngfft=information about FFT(dtset%ngfft a priori different from ngfftrec)
!!  mpi_enreg=information about MPI paralelisation
!!  rset<recursion_type> contains all parameter of recursion
!!  typat(natom)=type of pseudo potential associated to any atom
!!  natom=number of atoms
!!
!! OUTPUT
!!  enlx=non-local energy
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      vtorhorec
!!
!! CHILDREN
!!      recursion_nl,reshape_pot,timab,wrtout,xmpi_sum
!!
!! SOURCE

subroutine nlenergyrec(rset,enlx,exppot,ngfft,natom,typat,tsmear,trotter,tol)

!Arguments ------------------------------------
!Scalar
 integer , intent(in)  :: natom,trotter
 real(dp), intent(in)  :: tsmear,tol
 type(recursion_type),intent(in) :: rset
 real(dp), intent(out) :: enlx
!Arrays
 integer , intent(in)  :: typat(natom),ngfft(18)
 real(dp), intent(in)  :: exppot(0:ngfft(1)*ngfft(2)*ngfft(3)-1)
!Local variables-------------------------------
 integer :: iatom,jatom
 integer :: ii,ipsp,dim_trott
 integer :: ierr,me_count
 integer :: ilmn,jlmn,ilm,jlm,in,jn,il
 character(len=500) :: msg
 logical  :: tronc
 real(dp) :: rho_nl,normali,mult
 type(mpi_type):: mpi_loc
!Arrays
 integer  :: gcart_loc(3,natom)
 integer  :: ngfftrec(3),trasl(3)
 real(dp) :: tsec(2)
 real(dp) :: un0(0:rset%nfftrec)
 real(dp),pointer :: projec(:,:,:,:,:)
 real(dp),allocatable ::  exppotloc(:)
 real(dp) :: proj_arr(0:rset%ngfftrec(1)-1,0:rset%ngfftrec(2)-1,0:rset%ngfftrec(3)-1)

! *************************************************************************


 call timab(612,1,tsec) !!--start time-counter: nlenergyrec

 if(rset%debug)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' nlenergyrec : enter'
   call wrtout(std_out,msg,'PERS')
 end if

 write(msg,'(a)')' -- nlenergyrec -----------------------------------'
 call wrtout(std_out,msg,'COLL')

!--Initialisation variables
 enlx = zero
 mult = two !--is twice for non-spinned systems
 ngfftrec = rset%ngfftrec(:3)
 gcart_loc = rset%inf%gcart
 mpi_loc = rset%mpi
 me_count = 0
 dim_trott = max(0,2*trotter-1)
 nullify(projec)
 ABI_ALLOCATE(projec,(0:rset%ngfftrec(1)-1,0:rset%ngfftrec(2)-1,0:rset%ngfftrec(3)-1,rset%nl%lmnmax,natom))
 projec = zero

 tronc = rset%tronc  !--True if troncation is used
 if(tronc)   then
   ABI_ALLOCATE(exppotloc,(0:rset%nfftrec-1))
 end if


!--LOOP ON ATOMS to create projectors-vector
 atomloop1: do iatom = 1, natom
   ipsp = typat(iatom)
!  --Aquisition,reshape,translation,rotation of the projectors vector
   do ilmn = 1,rset%nl%lmnmax
     in = rset%nl%indlmn(3,ilmn,ipsp)
!    --Projectors vector in 3-composant vector
     projec(:,:,:,ilmn,iatom) = reshape(rset%nl%projec(:,ilmn,ipsp),shape=shape(projec(:,:,:,1,1)))
!    --Moving the projectors vector on the center of the grid
     do ii=1,3
       projec(:,:,:,ilmn,iatom) = cshift(projec(:,:,:,ilmn,iatom),shift=ngfftrec(ii)/2-gcart_loc(ii,iatom),dim=ii)
     end do
   end do

 end do atomloop1


!##################################################################
!--LOOP ON ATOMS (MAIN LOOP)
 atomloop: do iatom = 1, natom
   ipsp = typat(iatom)

!  --If troncation is present, the considered atom has to be in the
!  center of the grid so atoms, potential and projectors have to be translated
   if(tronc) then
     trasl = -rset%inf%gcart(:,iatom)+ngfftrec/2
!    --Translation of atoms
     do jatom=1,natom
       gcart_loc(:,jatom) = rset%inf%gcart(:,jatom)+trasl
       gcart_loc(:,jatom) = modulo(gcart_loc(:,jatom),ngfft(:3))
!      --Translation of non-local projectors
       do ilmn = 1,rset%nl%lmnmax
         projec(:,:,:,ilmn,jatom) = reshape(rset%nl%projec(:,ilmn,typat(jatom)),shape=shape(projec(:,:,:,1,1)))
         do ii=1,3
           projec(:,:,:,ilmn,jatom) = eoshift(projec(:,:,:,ilmn,jatom),shift=ngfftrec(ii)/2-gcart_loc(ii,jatom),dim=ii)
         end do
       end do
     end do

!    --Translation of the potential
     call reshape_pot(trasl,ngfft(1)*ngfft(2)*ngfft(3),rset%nfftrec,ngfft(:3),ngfftrec,exppot,exppotloc)
   end if

!  --Loop on projectors
   projloop: do ilmn = 1,rset%nl%lmnmax
     me_count = iatom+ilmn*natom-2 !--counter of the number of iteration
!    --Only the proc me compute
     if(mpi_loc%me==mod(me_count,mpi_loc%nproc)) then
       ilm = rset%nl%indlmn(4,ilmn,ipsp)
       proj_arr = zero
       do jlmn = 1,rset%nl%lmnmax
         jlm = rset%nl%indlmn(4,jlmn,ipsp)
         if(ilm==jlm) then
           in = rset%nl%indlmn(3,ilmn,ipsp)
           jn = rset%nl%indlmn(3,jlmn,ipsp)
           il = rset%nl%indlmn(1,ilmn,ipsp)+1
           proj_arr(:,:,:) = proj_arr(:,:,:) + rset%nl%eivec(jn,in,il,ipsp)*projec(:,:,:,jlmn,iatom)
!          write(std_out,*)'l,m,lm,n,n',il-1,rset%nl%indlmn(2,ilmn,ipsp),ilm,in,jn
!          write(std_out,*)'eigevectors',rset%nl%eivec(jn,in,il,ipsp)

         end if
       end do

       un0 = pack(proj_arr(:,:,:),mask=.true.)
       normali = sum(un0*un0)*rset%inf%ucvol
       un0 = (one/sqrt(normali))*un0

       if(tronc)then
         call recursion_nl(exppotloc,un0,rho_nl,rset,rset%ngfftrec,&
&         tsmear,trotter,dim_trott,tol,typat,&
&         natom,projec)
       else
         call recursion_nl(exppot,un0,rho_nl,rset,rset%ngfftrec,&
&         tsmear,trotter,dim_trott,tol,typat,&
&         natom,projec)
       end if

       enlx = enlx+mult*rho_nl*rset%nl%eival(in,il,ipsp)*normali
     end if

   end do projloop
 end do atomloop

!--Sum the contribution to the non-local energy computed by any procs
 call xmpi_sum(enlx,mpi_loc%comm_bandfft,ierr)

 if(associated(projec))  then
   ABI_DEALLOCATE(projec)
 end if
 if(tronc)  then
   ABI_DEALLOCATE(exppotloc)
 end if

 if(rset%debug)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' nlenergyrec : exit'
   call wrtout(std_out,msg,'PERS')
 end if

 call timab(612,2,tsec)  !--stop  time-counter: nlenergyrec

end subroutine nlenergyrec
!!***


!{\src2tex{textfont=tt}}
!!****f* ABINIT/first_rec
!! NAME
!! first_rec
!!
!! FUNCTION
!! When recursion method is used, in the first step this routine
!! compute some quantities which are used in the rest of the calculation.
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2019 ABINIT group (MMancini)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset:
!!   | recgratio =fine/coarse grid ratio
!!   | recptrott =trotter parameter
!!   | tsmear    =temperature
!!   | recrcut   =tut radius in recursion (range of iteration)
!!   | ngfft(18) =FFT grid used as real (fine) grid in recursion
!!  psps <type(pseudopotential_type)>=variables related to pseudo-potentials
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  rset <type(recursion_type)>=variables related to recursion method
!!   | debug<logical> = T if debugging is used
!!   | inf <type(metricrec_type)>=information concerning the infinitesimal metrics
!!   | ngfftrec(18) =truncated (or not, if not ngfftrec=ngfft)FFT grid used as real grid in recursion.
!!   | nfftrec =product(ngfftrec(1:3))
!!   | tronc<logical> = T if truncation is effectively used
!!   | ZT_p = fourier transform of the green_kernel calculated on the fine grid
!!
!!
!! NOTES
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      cpu_distribution,cudarec,get_pt0_pt1,green_kernel,init_nlpsprec
!!      random_number,recursion,reshape_pot,timab,timein,wrtout,xmpi_sum
!!
!! SOURCE

subroutine first_rec(dtset,psps,rset)

!Arguments ------------------------------------
! scalars
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(recursion_type),intent(inout) :: rset
!Local variables-------------------------------
!scalars
 integer  :: nfftrec,trotter,ii,dim_trott
 real(dp) :: tsmear,beta,rtrotter
 character(len=500) :: msg
!arrays
 integer  :: ngfftrec(18)
 real(dp) :: tsec(2)
#ifdef HAVE_GPU_CUDA
 integer  :: max_rec,ierr,testpts,swt_tm
 real(dp) :: rho,tm_ratio
 real(dp) :: time_cu,time_f
 type(recursion_type) :: rset_test
 type(recparall_type) :: parold
 integer :: trasl(3)
 real(dp) :: tsec2(2),tsec3(2)
 real(dp) :: aloc(0,1),b2loc(0,1)
 real(dp) :: dm_projec(0,0,0,1,1)
 real(dp) :: exppot(0:dtset%nfft-1)
 real(dp),allocatable :: exppotloc(:)
 real(cudap),allocatable :: aloc_cu(:),b2loc_cu(:)
#endif

! *************************************************************************

 call timab(601,1,tsec)  !!--Start time-counter: initialisation

 MSG_WARNING("RECURSION")
 if(dtset%recgratio>1) then
   write(msg,'(a)')'COARSE GRID IS USED'
   call wrtout(std_out,msg,'COLL')
 end if

!--Initialisation
 trotter = dtset%recptrott  !--Trotter parameter
 tsmear  = dtset%tsmear     !--Temperature
 beta    = one/tsmear       !--Inverse of temperature

!--Rewriting the trotter parameter
 dim_trott = max(0,2*trotter-1)
 rtrotter  = max(half,real(trotter,dp))

 write (msg,'(2a)')ch10,'==== FIRST CYCLE RECURSION ========================='
 call wrtout(std_out,msg,'COLL')


 ngfftrec = rset%ngfftrec
 nfftrec = rset%nfftrec
!------------------------------------------------
!--TRONCATION OF THE BOX: determines new dimensions
!--Now in InitRec
!--------------------------------------------------------
!--DEFINITION PAW VARIABLES COARSE-FINE GRID  TO USE TRANSGRID--INGRID FUNCTIONS
!--Now these variables are defined into gstate by InitRec

!--------------------------------------------------------
!--COMPUTATION OF THE FOURIER TRANSFORM OF THE GREEN KERNEL (only once)
 write (msg,'(a)')' - green kernel calculation -----------------------'
 call wrtout(std_out,msg,'COLL')
 ABI_ALLOCATE(rset%ZT_p,(1:2,0: nfftrec-1))
 call timab(601,2,tsec)
 call green_kernel(rset%ZT_p,rset%inf%rmet,rset%inf%ucvol,rtrotter/beta,rset%mpi,ngfftrec,nfftrec)
 call timab(601,1,tsec)
 write(msg,'(a,50a)')' ',('-',ii=1,50)
 call wrtout(std_out,msg,'COLL')
!!--end computation of the fourier transform of the Green kernel

!!-----------------------------------
!!--ROUTINE FOR THE CALCULATION OF THE NON-LOCAL PSEUDO
!--Now these variables here by  Init_nlpspRec
 call Init_nlpspRec(four*tsmear*rtrotter,psps,rset%nl,rset%inf,rset%ngfftrec,rset%debug)

!!-----------------------------------
!--Load distribution on procs when GPU are present
#if defined HAVE_GPU_CUDA

!--Test timing only if exists GPU and they are not equal to the cpus
 if(rset%tp == 4) then
   parold = rset%par
   ii = 0
   time_f = zero
   time_cu = zero
   call random_number(exppot)  !   exppot = one

   if(rset%gpudevice == -1) then
!    --Test CPUS
     swt_tm = 0
     testpts = min(rset%par%npt, 20)
     call timein(tsec2(1),tsec2(2))
     if(rset%tronc) then
       ABI_ALLOCATE(exppotloc,(0:nfftrec-1))
       do while(ii< testpts)
         trasl = -(/1,2,3/)+ngfftrec(:3)/2
         call reshape_pot(trasl,dtset%nfft,nfftrec,dtset%ngfft(:3),ngfftrec(:3),&
&         exppot,exppotloc)
         call recursion(exppotloc,0,0,0, &
&         aloc, &
&         b2loc, &
&         rho,&
&         0, rset%efermi,tsmear,rtrotter,dim_trott, &
&         rset%ZT_p, &
&         dtset%rectolden,dtset%typat, &
&         rset%nl,&
&         rset%mpi,nfftrec,ngfftrec,rset%inf,&
&         6,dtset%natom,dm_projec,0)
         ii=ii+1
       end do
       ABI_DEALLOCATE(exppotloc)
     else
       do while(ii< testpts)
         call recursion(exppot,0,0,0, &
&         aloc, &
&         b2loc, &
&         rho,&
&         0, rset%efermi,tsmear,rtrotter,dim_trott, &
&         rset%ZT_p, &
&         dtset%rectolden,dtset%typat, &
&         rset%nl,&
&         rset%mpi,nfftrec,ngfftrec,rset%inf,&
&         6,dtset%natom,dm_projec,0)
         ii=ii+1
       end do
     end if
     call timein(tsec3(1),tsec3(2))
     time_f = (tsec3(1)-tsec2(1))/real(testpts,dp)
     time_f = time_f*time_f
   else
!    --Test GPUS
     swt_tm = 1
     rset_test = rset
     rset_test%GPU%par%npt = max(rset%GPU%nptrec,100)
     rset_test%min_nrec = 0
     call get_pt0_pt1(dtset%ngfft(:3),dtset%recgratio,0,&
&     rset_test%GPU%par%npt,rset_test%GPU%par)


     ABI_ALLOCATE(aloc_cu,(rset_test%GPU%par%npt))
     ABI_ALLOCATE(b2loc_cu,(rset_test%GPU%par%npt))
     call timein(tsec2(1),tsec2(2))
     call cudarec(rset_test, exppot,aloc_cu,b2loc_cu,&
&     beta,trotter,dtset%rectolden,dtset%recgratio,dtset%ngfft,max_rec)
     call timein(tsec3(1),tsec3(2))
     ABI_DEALLOCATE(aloc_cu)
     ABI_DEALLOCATE(b2loc_cu)

     time_cu = (tsec3(1)-tsec2(1))/real(rset_test%GPU%par%npt,dp)
     time_cu = time_cu*time_cu
   end if


!  --Get Total Times
   call xmpi_sum(time_f,rset%mpi%comm_bandfft,ierr)
   call xmpi_sum(time_cu,rset%mpi%comm_bandfft,ierr)

!  --Average Total Times
   time_f   = sqrt(time_f/real(rset%mpi%nproc-rset%ngpu,dp))
   time_cu  = sqrt(time_cu/real(rset%ngpu,dp))
   tm_ratio = time_f/time_cu


   write(msg,'(3(a25,f10.5,a))')&
&   ' Time for cpu recursion ',time_f,ch10,&
&   ' Time for gpu recursion ',time_cu,ch10,&
&   ' Time ratio             ',tm_ratio,ch10
   call wrtout(std_out,msg,'COLL')


!  tm_ratio =1.20d2! 0.d0! 1.21d0
   rset%par = parold
!  --Compute the work-load distribution on devices (gpu,cpu)
   if(tm_ratio>1.5d0 .and. time_cu>zero)then
     rset%load = 1
     call cpu_distribution(dtset%recgratio,rset,dtset%ngfft(:3),tm_ratio,1)
   else
     rset%gpudevice = -1
   end if
 end if

#endif


!------------------------------------------------------------
!--DETERMINING WHICH POINT WILL COMPUTE THAT PROC
!--Now these variables are defined into gstate by Init_rec

 write (msg,'(2a)')ch10,'==== END FIRST CYCLE RECURSION ====================='
 call wrtout(std_out,msg,'COLL')
 call timab(601,2,tsec) !!--stop time-counter: initialisation

end subroutine first_rec
!!***


!!****f* ABINIT/green_kernel
!! NAME
!! green_kernel
!!
!! FUNCTION
!! this routine compute the fourrier transform of the Green kernel for the
!! recursion method
!!
!! INPUTS
!!  inf_rmet=define the  infinitesimal metric : rprimd*(transpose(rprimd)) divided
!!    by the number of discretisation point
!!  inf_ucvol=volume of infinitesimal cell
!!  mult=variance of the Gaussian (=rtrotter/beta)
!!  mpi_enreg=information about MPI parallelization
!!  ngfft=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nfft=total number of fft grid points
!!  debug_rec=debugging variable
!!
!! OUTPUT
!!  ZT_p=fourier transforme of the Green kernel
!!
!! PARENTS
!!      first_rec
!!
!! CHILDREN
!!      fourdp,timab,wrtout
!!
!! NOTES
!!  at this time :
!!       - need a rectangular box
!!
!! SOURCE


subroutine green_kernel(ZT_p,inf_rmet,inf_ucvol,mult,mpi_enreg,ngfft,nfft)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nfft
 real(dp),intent(in) :: inf_ucvol,mult
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: inf_rmet(3,3)
 real(dp),intent(out) :: ZT_p(1:2,0:nfft-1)

!Local variables-------------------------------
!scalars
 integer,parameter :: n_green_max=5
 integer :: ii,isign,jj,kk,n_green,xx,yy,zz
 real(dp) :: acc, norme
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: T_p(:)

! *************************************************************************

 call timab(603,1,tsec)

 norme = (mult/pi)**(onehalf)

 ABI_ALLOCATE(T_p,(0:nfft-1))

!n_green should be better chosen for non rectangular cell
 do xx=1, n_green_max
   n_green = xx
   if(exp(-mult*dsq_green(xx*ngfft(1),0,0,inf_rmet))<tol14 &
&   .and. exp(-mult*dsq_green(0,xx*ngfft(2),0,inf_rmet))<tol14 &
&   .and. exp(-mult*dsq_green(0,0,xx*ngfft(3),inf_rmet))<tol14 ) exit
 end do

 acc = zero
 T_p = zero
 do kk = 0,ngfft(3)-1
   do jj = 0,ngfft(2)-1
     do ii = 0,ngfft(1)-1

       do xx=-n_green,n_green-1
         do yy=-n_green,n_green-1
           do zz=-n_green,n_green-1

             T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk) = T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk)+ &
&             exp(-mult*dsq_green(ii+xx*ngfft(1),jj+yy*ngfft(2),kk+zz*ngfft(3),inf_rmet))

           end do
         end do
       end do

       T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk) = norme*T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk)
       acc = acc + inf_ucvol* T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk)

     end do
   end do
 end do

 T_p(:)= (one/acc)*T_p(:)

!if(debug_rec)then
 write(msg,'(a,d12.3,2(2a,i8),2(2a,3d12.3),2a,d16.6)')&
& ' on the boundary    ', exp(-mult*dsq_green(ngfft(1),0,0,inf_rmet)),ch10, &
& ' no zero            ', count(T_p>tol14),ch10, &
& ' n_green            ', n_green,ch10, &
& ' erreur_n_green     ', exp(-mult*dsq_green(n_green*ngfft(1),0,0,inf_rmet)), &
& exp(-mult*dsq_green(0,n_green*ngfft(2),0,inf_rmet)), &
& exp(-mult*dsq_green(0,0,n_green*ngfft(3),inf_rmet)),ch10,&
& ' erreur_troncat     ', T_p(ngfft(1)/2),  &
& T_p(ngfft(1)*(ngfft(2)/2)), &
& T_P(ngfft(1)*ngfft(2)*(ngfft(3)/2)),ch10, &
& ' erreurT_p          ',abs(acc-1.d0)
 call wrtout(std_out,msg,'COLL')
!endif


 isign = -1
 call fourdp(1,ZT_p,T_p,isign,mpi_enreg,nfft,1,ngfft,0)

 ABI_DEALLOCATE(T_p)

 ZT_p(:,:) = real(nfft,dp)*ZT_p


 call timab(603,2,tsec)

 contains

   function dsq_green(ii,jj,kk,inf_rmet)

   real(dp) :: dsq_green
   integer,intent(in) :: ii,jj,kk
   real(dp),intent(in) :: inf_rmet(3,3)
   dsq_green= inf_rmet(1,1)*dble(ii**2)&
&   +inf_rmet(2,2)*dble(jj**2)&
&   +inf_rmet(3,3)*dble(kk**2)&
&   +two*(inf_rmet(1,2)*dble(ii*jj)&
&   +inf_rmet(2,3)*dble(jj*kk)&
&   +inf_rmet(3,1)*dble(kk*ii))
 end function dsq_green

end subroutine green_kernel
!!***


!!****f* ABINIT/recursion
!! NAME
!! recursion
!!
!! FUNCTION
!! This routine computes the recursion coefficients and the corresponding
!! continued fraction to get the density at a point from a fixed potential.
!!
!! INPUTS
!!  exppot=exponential of -1/tsmear*vtrial (computed only once in vtorhorec)
!!  coordx, coordy, coordz=coordonnees of the computed point
!!  nrec=order of recursion
!!  fermie=fermi energy (Hartree)
!!  tsmear=temperature (Hartree)
!!  dim_trott=dimension of the partial fraction decomposition
!!  rtrotter=trotter parameter (real)
!!  ZT_p=fourier transform of the Green krenel (computed only once in vtorhorec)
!!  typat(:)=type of psp associated to any atom
!!  tol=tolerance criteria for stopping recursion
!!  debug=debugging variable
!!  mpi_enreg=information about MPI paralelisation
!!  nfft=number of points in FFT grid
!!  ngfft=information about FFT
!!  metrec<type(metricrec_type)>=information concerning the infinitesimal metrics
!!  inf_ucvol=infinitesimal unit cell volume
!!  tim_fourdp=time counter for fourdp
!!  natom=number of atoms
!!  projec(ngfftrec(1),ngfftrec(2),ngfftrec(3),lmnmax,natom) is the  vector, on the ngfftrec grid containing
!!  the non-lacal projector $Y_{lm}(r-R_A)f_{lk}(r-R_A)
!!  tim= 0 if the time spent in the routine is not taken into account,1 otherwise. For example
!!  when measuring time for loading  balancing, we don't want to add the time spent in this to the
!!  total time calculation
!!
!! OUTPUT
!!  rho_out=result of the continued fraction multiplied by a multiplicator
!!  an, bn2 : coefficient given by recursion.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      first_rec,vtorhorec
!!
!! CHILDREN
!!      fourdp,timab,trottersum,vn_nl_rec
!!
!! NOTES
!!  at this time :
!!       - exppot should be replaced by ?
!!       - coord should be replaced by ?
!!       - need a rectangular box (rmet diagonal matrix)
!!
!! SOURCE

subroutine recursion(exppot,coordx,coordy,coordz,an,bn2,rho_out, &
&                    nrec,fermie,tsmear,rtrotter,dim_trott, &
&                    ZT_p, tol,typat, &
&                    nlrec,mpi_enreg,&
&                    nfft,ngfft,metrec,&
&                    tim_fourdp,natom,projec,tim)


 use m_linalg_interfaces

!Arguments -------------------------------
!scalars
 integer,intent(in) :: coordx,coordy,coordz,nfft,nrec,tim
 integer,intent(in) :: tim_fourdp,natom,dim_trott
 real(dp),intent(in) :: fermie,tol,tsmear,rtrotter
 real(dp), intent(out) :: rho_out
 type(MPI_type),intent(in) :: mpi_enreg
 type(nlpsprec_type),intent(in) :: nlrec
 type(metricrec_type),intent(in) :: metrec
!arrays
 integer, intent(in) :: ngfft(18)
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: ZT_p(1:2, 0:nfft-1)
 real(dp), intent(in) :: exppot(0:nfft-1)
 real(dp), intent(in) :: projec(0:,0:,0:,1:,1:)
 real(dp), intent(out) :: an(0:nrec),bn2(0:nrec)
!Local variables-------------------------------
!not used, debugging purpose only
!for debugging purpose, detailled printing only once for density and ekin
!scalars
 integer, parameter :: level = 7, minrec = 3
 integer  :: irec,isign,timab_id,ii
 real(dp) :: switchimu,switchu
 real(dp) :: bb,beta,mult,prod_b2,error,errold
 real(dp) :: inf_ucvol,pi_on_rtrotter,twortrotter,exp1,exp2
 complex(dpc) :: cinv2rtrotter,coeef_mu,facrec0
! character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp) :: inf_tr(3)
 real(dp) :: Zvtempo(1:2, 0:nfft-1)
 real(dp) :: unold(0:nfft-1),vn(0:nfft-1),un(0:nfft-1)
 complex(dpc) :: acc_rho(0:nrec)
 complex(dpc) :: D(0:dim_trott),Dold(0:dim_trott)
 complex(dpc) :: N(0:dim_trott),Nold(0:dim_trott)

! *************************************************************************

!--If count time or not
 timab_id = 616; if(tim/=0) timab_id = 606;

 call timab(timab_id,1,tsec)

!##############################################################
!--Initialisation of metrics
 inf_ucvol = metrec%ucvol
 inf_tr = metrec%tr
 mult = two/inf_ucvol    !non-spined system

 beta = one/tsmear
!--Variables for optimisation
 pi_on_rtrotter = pi/rtrotter
 twortrotter = two*rtrotter
 exp1 = exp((beta*fermie)/(rtrotter))
 exp2 = exp(beta*fermie/(twortrotter))
 cinv2rtrotter = cmplx(one/twortrotter,zero,dp)
 coeef_mu = cmplx(one/exp2,zero,dp)

!--Initialisation of  an,bn,un....
 N = czero;  D = cone
 facrec0 = cone
 Nold = czero; Dold = czero

 an = zero; bn2 = zero;  bn2(0) = one
 bb = zero; vn  = zero;  unold  = zero
!--u0 is a Dirac function
 un = zero
 un(coordx+ngfft(1)*(coordy+ngfft(2)*coordz)) = one/sqrt(inf_ucvol)

!--Initialisation of accumulated density
 acc_rho = czero
!--Initialisation of estimated error
 prod_b2 = twortrotter/exp1
 errold = zero

!##############################################################
!--Main loop
 maindo : do irec = 0, nrec

!  --Get an and bn2 coef by the lanczos method

!  --Computation of exp(-beta*V/8*p)*un or exp(-beta*V/4*p)*un
!  depending on if nl part has to be calculated or not.
   vn = exppot * un

!  --First Non-local psp contribution: (Id+sum_atom int dr1(E(r,r1))vn(r1))
!  --Computation of exp(-beta*V_NL/4*p)*vn
   if(nlrec%nlpsp) then
     call timab(timab_id,2,tsec)
     call vn_nl_rec(vn,natom,typat,ngfft(:3),inf_ucvol,nlrec,projec)
     call timab(timab_id,1,tsec)

!    --Computation of exp(-beta*V/8*p)*vn in nonlocal case
     vn = exppot * vn
   end if !--End if on nlrec%nlpsp

!  --Convolution with the Green kernel
!  --FFT of vn
   isign = -1
   call fourdp(1,Zvtempo,vn,isign,mpi_enreg,nfft,1,ngfft,tim_fourdp)

!  --F(T)F(vn)
   do ii = 0,nfft-1
     switchu   = Zvtempo(1,ii)
     switchimu = Zvtempo(2,ii)
     Zvtempo(1,ii) = switchu*ZT_p(1,ii) - switchimu*ZT_p(2,ii)
     Zvtempo(2,ii) = switchu*ZT_p(2,ii) + switchimu*ZT_p(1,ii)
   end do

!  --F^-1(F(T)F(vn))
   isign = 1
   call fourdp(1,Zvtempo,vn,isign,mpi_enreg,nfft,1,ngfft,tim_fourdp)

!  --Computation of exp(-beta*V/8*p)*un or exp(-beta*V/4*p)*un
!  depending on if nl part has to be calculated or not.

   vn = inf_ucvol * exppot * vn

!  --Second Non-local psp contribution: (Id+sum_atom E(r,r1))vn
   if(nlrec%nlpsp) then
     call timab(timab_id,2,tsec)
     call vn_nl_rec(vn,natom,typat,ngfft(:3),inf_ucvol,nlrec,projec)
     call timab(timab_id,1,tsec)

!    --Computation of exp(-beta*V/8*p)*vn in nonlocal case
     vn = exppot * vn
   end if !--End if on nlrec%nlpsp

!  --Multiplication of a and b2 coef by exp(beta*fermie/(two*rtrotter)) must be done in the continued fraction computation
!  --Computation of a and b2
   an(irec) = inf_ucvol*ddot(nfft,vn,1,un,1)

!  --an must be positive real
!  --We must compute bn2 and prepare for the next iteration
   if(irec<nrec)then
     do ii = 0,nfft-1
       switchu = un(ii)
       un(ii) = vn(ii)-an(irec)*un(ii)-bb*unold(ii)
       unold(ii) = switchu
       bn2(irec+1) = bn2(irec+1)+inf_ucvol*un(ii)*un(ii)
     end do
     bb = sqrt(bn2(irec+1))
     un = (one/bb)*un
   end if

!  ######################################################
!  --Density computation
!  density computation is done inside the main looping, juste after the calculus of a and b2, in order to make
!  it possible to stop the recursion at the needed accuracy, without doing more recursion loop than needed -
!  further developpement

!  !--using the property that: sum_i(bi*c)^2|(z-ai*c)=1/c*sum_i(bi)^2|(z/c-ai)
!  !and for c =exp(-beta*fermie/(two*rtrotter)


   call trottersum(dim_trott,error,prod_b2,pi_on_rtrotter,&
&   facrec0,coeef_mu,exp1,&
&   an(irec),bn2(irec),&
&   N,D,Nold,Dold)


   if(irec/=nrec .and. irec>=minrec)then
     if((bn2(irec+1)<tol14).or.(mult*error<tol.and.errold<tol)) exit
   end if
   errold = mult*error
 end do maindo
!--Accumulated density
 rho_out = mult*real(cone-sum(N/D)*cinv2rtrotter,dp)


 call timab(timab_id,2,tsec)

 end subroutine recursion
!!***


!!****f* ABINIT/recursion_nl
!! NAME
!! recursion_nl
!!
!! FUNCTION
!! Given a $|un>$ vector on the real-space grid this routine calculates
!! the density in  $|un>$ by recursion method.
!!
!! INPUTS
!!  exppot=exponential of -1/tsmear*vtrial (computed only once in vtorhorec)
!!  trotter=trotter parameter
!!  dim_trott=dimension of the partial fraction decomposition
!!  tsmear=temperature (Hartree)
!!  tol=tolerance criteria for stopping recursion_nl
!!  ngfft=information about FFT(dtset%ngfft a priori different from ngfftrec)
!!  rset<recursion_type> contains all parameter of recursion
!!  typat(natom)=type of pseudo potential associated to any atom
!!  natom=number of atoms
!!  projec(ngfftrec(1),ngfftrec(2),ngfftrec(3),lmnmax,natom) is the  vector, on the ngfftrec grid containing
!!  the non-lacal projector $Y_{lm}(r-R_A)f_{lk}(r-R_A)
!!
!! OUTPUT
!!  rho_out=result of the continued fraction multiplied by a multiplicator
!!
!! SIDE EFFECTS
!!  un(:,:,:)=initial vector on the grid. it is changed in output
!!
!! PARENTS
!!      nlenergyrec
!!
!! CHILDREN
!!      fourdp,timab,trottersum,vn_nl_rec,wrtout
!!
!! NOTES
!!  at this time :
!!       - need a rectangular box (rmet diagonal matrix)
!!
!! SOURCE

subroutine recursion_nl(exppot,un,rho_out,rset,ngfft, &
  &                     tsmear,trotter,dim_trott,tol,typat,&
  &                     natom,projec)


 use m_linalg_interfaces

!Arguments -------------------------------
!scalars
 integer,intent(in) :: trotter,natom,dim_trott
 real(dp),intent(in) :: tol,tsmear
 type(recursion_type),intent(in) :: rset
 real(dp), intent(out) :: rho_out
!arrays
 integer,intent(in) ::  typat(natom),ngfft(18)
 real(dp),intent(in) :: exppot(0:ngfft(1)*ngfft(2)*ngfft(3)-1)
 real(dp),intent(inout) :: un(0:rset%nfftrec-1)
 real(dp),pointer :: projec(:,:,:,:,:)
!Local variables-------------------------------
!scalars
 integer, parameter ::  minrec = 3
 integer  :: irec,isign,ii
 real(dp) :: bb,beta,mult,prod_b2,rtrotter
 real(dp) :: inf_ucvol,pi_on_rtrotter,twortrotter,exp1
 real(dp) :: exp2,error,errold
 real(dp) :: switchu,switchimu
 complex(dpc) :: facrec0,cinv2rtrotter,coeef_mu
 character(len=500) :: msg
 type(mpi_type),pointer:: mpi_loc
!arrays
 real(dp):: tsec(2)
 real(dp):: inf_tr(3)
 real(dp):: an(0:rset%min_nrec),bn2(0:rset%min_nrec)
 real(dp):: vn(0:rset%nfftrec-1)
 real(dp):: unold(0:rset%nfftrec-1)
 real(dp):: Zvtempo(1:2,0:rset%nfftrec-1)
 complex(dpc) :: acc_rho(0:rset%min_nrec)
 complex(dpc) :: D(0:dim_trott),Dold(0:dim_trott)
 complex(dpc) :: N(0:dim_trott),Nold(0:dim_trott)

! *************************************************************************

 call timab(608,1,tsec) !--start time-counter: recursion_nl
 if(rset%debug)then
   msg=' '
   call wrtout(std_out,msg,'COLL')
 end if

!##############################################################
 beta = one/tsmear

!--Rewriting the trotter parameter
 rtrotter  = max(half,real(trotter,dp))

!--Initialisation of mpi
 mpi_loc => rset%mpi

!--Initialisation of metrics
 inf_ucvol = rset%inf%ucvol
 inf_tr = rset%inf%tr
 mult = one   !--In the case of the calculus of the NL-energy

!--Initialisation of  an,bn,un....
 N = czero;  D = cone
 facrec0 = cone
 Nold = czero; Dold = czero

 an = zero; bn2 = zero;  bn2(0) = one
 bb = zero; vn  = zero;  unold  = zero

!--Variables for optimisation
 pi_on_rtrotter = pi/rtrotter
 twortrotter = two*rtrotter
 exp1 = exp((beta*rset%efermi)/(rtrotter))
 exp2 = exp(beta*rset%efermi/(twortrotter))
 cinv2rtrotter = cmplx(one/twortrotter,zero,dp)
 coeef_mu = cmplx(one/exp2,zero,dp)

!--Initialisation of accumulated density
 acc_rho = czero
!--Initialisation of estimated error
 prod_b2 = twortrotter/exp1
 errold = zero

!##############################################################
!--Main loop
 maindo : do irec = 0, rset%min_nrec
!  --Get an and bn2 coef by the lanczos method

!  --Computation of exp(-beta*V/8*p)*un
   vn = exppot * un

!  --First Non-local psp contribution: (Id+sum_atom E(r,r1))vn
   call timab(608,2,tsec)
   call vn_nl_rec(vn,natom,typat,rset%ngfftrec(:3),inf_ucvol,rset%nl,projec)
   call timab(608,1,tsec)

!  --Computation of exp(-beta*V/8*p)*un
   vn = exppot * vn

!  --Convolution with the Green kernel
!  --FFT of vn
   isign = -1
   call fourdp(1,Zvtempo,vn,isign,mpi_loc,rset%nfftrec,1,rset%ngfftrec,6)

!  --F(T)F(vn)
   do ii = 0,rset%nfftrec-1
     switchu   = Zvtempo(1,ii)
     switchimu = Zvtempo(2,ii)
     Zvtempo(1,ii) = switchu*rset%ZT_p(1,ii) - switchimu*rset%ZT_p(2,ii)
     Zvtempo(2,ii) = switchu*rset%ZT_p(2,ii) + switchimu*rset%ZT_p(1,ii)
   end do

!  --F^-1(F(T)F(vn))
   isign = 1
   call fourdp(1,Zvtempo,vn,isign,mpi_loc,rset%nfftrec,1,rset%ngfftrec,6)

!  --Computation of exp(-beta*V/2*p)*vn
   vn = inf_ucvol * exppot * vn

!  --Second Non-local psp contribution: (Id+sum_atom E(r,r1))vn
   call timab(608,2,tsec)
   call vn_nl_rec(vn,natom,typat,rset%ngfftrec(:3),inf_ucvol,rset%nl,projec)
   call timab(608,1,tsec)

!  --Computation of exp(-beta*V/8*p)*vn
   vn = exppot * vn


!  --Multiplication of a and b2 coef by exp(beta*fermie/(2.d0*rtrotter)) must be done in the continued fraction computation
!  --Computation of a and b2
   an(irec) = inf_ucvol*ddot(rset%nfftrec,vn,1,un,1)      !--an must be positive real

!  --We must compute bn2 and prepare for the next iteration
   if(irec<rset%min_nrec)then
     do ii = 0,rset%nfftrec-1
       switchu = un(ii)
       un(ii) = vn(ii)-an(irec)*un(ii)-bb*unold(ii)
       unold(ii) = switchu
       bn2(irec+1) = bn2(irec+1)+inf_ucvol*un(ii)*un(ii)
     end do
     bb = sqrt(bn2(irec+1))
     un = (one/bb)*un
   end if

!  ######################################################
!  --Density computation
!  in order to make it possible to stop the recursion_nl at the
!  needed accuracy, without doing more recursion_nl loop than needed further developpement

   call trottersum(dim_trott,error,&
&   prod_b2,pi_on_rtrotter,&
&   facrec0,coeef_mu,exp1,&
&   an(irec),bn2(irec),&
&   N,D,Nold,Dold)


   if(irec/=rset%min_nrec .and. irec>=minrec)then
     if((bn2(irec+1)<tol14).or.(mult*error<tol.and.errold<tol)) exit
   end if
   errold = mult*error
 end do maindo
!--Accumulated density
 rho_out = mult*real(cone-sum(N/D)*cinv2rtrotter,dp)

 call timab(608,2,tsec) !--stop time-counter: recursion_nl

end subroutine recursion_nl
!!***


!!****f* ABINIT/vn_nl_rec
!! NAME
!! vn_nl_rec
!!
!! FUNCTION
!! this routine computes the contribution to the vector vn, during
!! recursion, due to the non-local psp.
!!
!! INPUTS
!!  vn(:,:,:)=the vector on the real-space grid.
!!  inf_ucvol=volume of infinitesimal cell
!!  natom=number of atoms
!!  typat(natom)=the type of psps associated to the atoms
!!  ngfftrec(3)=first 3 components of ngfftrec (truncated box, if different from ngfft) for the real-space grid
!!  nlrec<type(nlpsprec_type)> in recursion_type containing information concerning psp
!!  projec(ngfftrec(1),ngfftrec(2),ngfftrec(3),lmnmax,natom) is the  vector, on the ngfftrec grid containing
!!  the non-lacal projector $Y_{lm}(r-R_A)f_{lk}(r-R_A)
!!
!! OUTPUT
!! vn_nl(:,:,:)=the non_local contribution to vn
!!
!! PARENTS
!!      recursion,recursion_nl
!!
!! CHILDREN
!!      timab
!!
!! NOTES
!!
!! SOURCE

subroutine vn_nl_rec(vn,natom,typat,ngfftrec,inf_ucvol,nlrec,projec)


 use m_linalg_interfaces

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom
 real(dp),intent(in) :: inf_ucvol
 type(nlpsprec_type),intent(in) :: nlrec
!arrays
 integer,intent(in) :: ngfftrec(3),typat(natom)
 real(dp),intent(in) :: projec(0:,0:,0:,1:,1:)
 real(dp),intent(inout):: vn(0:ngfftrec(1)*ngfftrec(2)*ngfftrec(3)-1)
!Local variables-------------------------------
!scalars
 integer :: iatom,nfftrec
 integer :: jlmn,il,in,jn
 integer :: ipsp,ilmn
 integer :: npsp,lmnmax
 real(dp):: vn_nl_loc
!arrays
 real(dp):: vn_nl(0:ngfftrec(1)-1,0:ngfftrec(2)-1,0:ngfftrec(3)-1)
 real(dp):: vtempo(0:ngfftrec(1)-1,0:ngfftrec(2)-1,0:ngfftrec(3)-1)
 real(dp):: tsec(2)
! *************************************************************************

 call timab(615,1,tsec)
!--Initialisation

 vn_nl = zero
 npsp = nlrec%npsp
 lmnmax = nlrec%lmnmax
 nfftrec = product(ngfftrec)
 vtempo(:,:,:) = reshape(source=vn,shape=ngfftrec(:3))

!--Sum_iatom \int dr1 E(r-r_a,r1-r_a)vn(r1) *infucvol
 do iatom=1,natom !--Loop on atoms
   ipsp = typat(natom)

!  --If psp(typat(iatom)) is local then cycle
   if(all(nlrec%pspinfo(:,ipsp)==0))  cycle


!  write(std_out,*)'lmnmax',nlrec%lmnmax,lmnmax

   do ilmn = 1, lmnmax
     do jlmn = 1,lmnmax
       if(nlrec%indlmn(4,ilmn,ipsp)==nlrec%indlmn(4,jlmn,ipsp)) then
         il = 1+nlrec%indlmn(1,jlmn,ipsp)
         in = nlrec%indlmn(3,ilmn,ipsp)
         jn = nlrec%indlmn(3,jlmn,ipsp)
         vn_nl_loc = ddot(nfftrec,projec(:,:,:,jlmn,iatom),1,vtempo,1)
         vn_nl = vn_nl+projec(:,:,:,ilmn,iatom)*vn_nl_loc*nlrec%mat_exp_psp_nl(in,jn,il,ipsp)
       end if
     end do
   end do
 end do !--End loop on atoms
 vtempo = vtempo + vn_nl*inf_ucvol

 vn = reshape(source=vtempo,shape=(/nfftrec/))

 call timab(615,2,tsec)

end subroutine vn_nl_rec
!!***

end module m_vtorhorec
!!***
