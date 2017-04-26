!!{\src2tex{textfont=tt}}
!!****f* ABINIT/vtorhorec
!! NAME
!! vtorhorec
!!
!! FUNCTION
!! This routine computes the new density from a fixed potential (vtrial)
!! using a recursion method
!!
!! COPYRIGHT
!! Copyright (C) 2008-2017 ABINIT group (SLeroux,MMancini).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!  enl=nonlocal pseudopotential part of total energy.
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine vtorhorec(dtset,&
&  ek,enl,entropy,e_eigenvalues,fermie,&
&  grnl,initialized,irrzon,nfftf,phnons,&
&  rhog, rhor, vtrial,rset,deltastep,rprimd,gprimd)

 use defs_basis
 use defs_abitypes
 use defs_rectypes
 use m_xmpi
 use m_pretty_rec
 use m_errors
 use m_profiling_abi

 use m_rec,            only : Calcnrec
 use m_rec_tools,      only : reshape_pot
#ifdef HAVE_GPU_CUDA
 use m_initcuda,       only : cudap
 use m_hidecudarec
 use m_xredistribute
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vtorhorec'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_65_paw
 use interfaces_67_common
 use interfaces_68_recursion, except_this_one => vtorhorec
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: initialized
 integer,intent(in) :: nfftf,deltastep
 real(dp),intent(out) :: e_eigenvalues,ek,enl,entropy,fermie
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
 enl = zero
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
     call nlenergyrec(rset,enl,exppot,dtset%ngfft,dtset%natom,&
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

   ek = e_eigenvalues-intrhov-enl


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
&   ' non-loc ene =',enl
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
& dtset%nsppol,dtset%nsym,dtset%paral_kgb,&
& phnons,rhog,rhor,rprimd,dtset%symafm,dtset%symrel)

end subroutine vtorhorec
!!***
