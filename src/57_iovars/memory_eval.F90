!{\src2tex{textfont=tt}}
!!****f* ABINIT/memory_eval
!! NAME
!! memory_eval
!!
!! FUNCTION
!! Big loop on the datasets :
!! - for each of the datasets, write one line about the crystallographic data
!! - compute the memory needs for this data set.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group ()
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number of output file
!!  mpi_enregs=information about MPI parallelization
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!      one data set.
!!  npsp=number of pseudopotentials
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!   printing only
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      abi_io_redirect,getdim_nloc,getmpw,getng,libpaw_write_comm_set
!!      littlegroup_q,mati3inv,memorf,memory,metric,mkrdim,prtspgroup,setmqgrid
!!      wvl_memory
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine memory_eval(dtsets,iout,mpi_enregs,ndtset,ndtset_alloc,npsp,pspheads)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_libpaw_tools, only : libpaw_write_comm_set

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memory_eval'
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_52_fft_mpi_noabirule
 use interfaces_56_recipspace
 use interfaces_57_iovars, except_this_one => memory_eval
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,ndtset,ndtset_alloc,npsp
 type(MPI_type),intent(inout) :: mpi_enregs(0:ndtset_alloc)
!arrays
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables -------------------------------
!scalars
 integer :: cplex,exchn2n3d,extrapwf,getcell,idtset,ii,intxc,densfor_pred,iprcel
 integer :: iscf,isym,jdtset,lmnmax,mem_test
 integer :: lmnmax_eff,lmnmaxso,lnmax,lnmax_eff,lnmaxso,mband
 integer :: me_fft,mffmem,mgfftdiel,mgfftf,mkmem,mpsang,mpspso
 integer :: mpssoang,mpw,mqgrid,mqgriddg,mqgrid_ff,mqgrid_vl,n1xccc,natom
 integer :: nfftdiel,nfftf,nkpt,nproc_fft,nptsgvec,npulayit,npwdiel,nspden,nspinor
 integer :: nsppol,nsym,ntypat,occopt,optddk,optforces,optphon,optstress
 integer :: optstrs,paral_fft,pawcpxocc,pawmixdg,pawnhatxc,pawspnorb,pawstgylm,prtvol,ptgroupma,response
 integer :: spgroup,timrev,usepaw,useylm,use_gpu_cuda,xclevel
 real(dp) :: diecut,dilatmx,ecut,ecut_eff,ecutdg_eff,ecutsus,ucvol
!arrays
 integer :: bravais(11),mkmems(3),ngfftdiel(18)
 integer :: ngfftf(18),nloalg(3)
 integer,allocatable :: nband(:),symq(:,:,:),symrec(:,:,:),symrel(:,:,:)
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: genafm(3),gmet(3,3),gprimd(3,3),kpt_diel(3),qphon(3),rmet(3,3),rprimd(3,3)

!*************************************************************************

 do idtset=1,ndtset_alloc
   if(mpi_enregs(idtset)%me<0) cycle
   call abi_io_redirect(new_io_comm=mpi_enregs(idtset)%comm_world)
   call libpaw_write_comm_set(mpi_enregs(idtset)%comm_world)

!  Initialisations
   bravais(:)=dtsets(idtset)%bravais(:)
   exchn2n3d=dtsets(idtset)%exchn2n3d
   extrapwf=dtsets(idtset)%extrapwf
   genafm(:) =dtsets(idtset)%genafm(:)
   getcell=dtsets(idtset)%getcell
   intxc=dtsets(idtset)%intxc
   densfor_pred=dtsets(idtset)%densfor_pred
   iprcel=dtsets(idtset)%iprcel
   iscf=dtsets(idtset)%iscf
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
   me_fft=mpi_enregs(idtset)%me_fft
   mffmem=dtsets(idtset)%mffmem
   mpw=dtsets(idtset)%mpw
   mqgrid=dtsets(idtset)%mqgrid
   mqgriddg=dtsets(idtset)%mqgriddg
   natom=dtsets(idtset)%natom
   nkpt  =dtsets(idtset)%nkpt
   nloalg(:)=dtsets(idtset)%nloalg(:)
   nproc_fft=mpi_enregs(idtset)%nproc_fft
   npulayit=dtsets(idtset)%npulayit
   nspden=dtsets(idtset)%nspden
   nspinor=dtsets(idtset)%nspinor
   nsppol=dtsets(idtset)%nsppol
   nsym     =dtsets(idtset)%nsym
   ntypat=dtsets(idtset)%ntypat
   occopt=dtsets(idtset)%occopt
   optforces=dtsets(idtset)%optforces
   paral_fft=mpi_enregs(idtset)%paral_kgb
   pawcpxocc=dtsets(idtset)%pawcpxocc
   pawmixdg=dtsets(idtset)%pawmixdg
   pawnhatxc=dtsets(idtset)%pawnhatxc
   pawspnorb=dtsets(idtset)%pawspnorb
   pawstgylm=dtsets(idtset)%pawstgylm
   prtvol=dtsets(idtset)%prtvol
   ptgroupma =dtsets(idtset)%ptgroupma
   qphon(:)=dtsets(idtset)%qptn(:)
   spgroup   =dtsets(idtset)%spgroup
   usepaw=dtsets(idtset)%usepaw
   useylm=dtsets(idtset)%useylm
   use_gpu_cuda=dtsets(idtset)%use_gpu_cuda
   xclevel=dtsets(idtset)%xclevel

   ABI_ALLOCATE(symrel,(3,3,nsym))
   symrel(:,:,1:nsym)=dtsets(idtset)%symrel(:,:,1:nsym)

!  Space group output
   call prtspgroup(bravais,genafm,iout,jdtset,ptgroupma,spgroup)

   if (dtsets(idtset)%toldff>tol16.and.optforces==0) optforces=1
   if (dtsets(idtset)%tolrff>tol16.and.optforces==0) optforces=1
   if (dtsets(idtset)%ionmov>tol16.and.optforces==0) optforces=1
   if (dtsets(idtset)%imgmov>tol16.and.optforces==0) optforces=1
   optstress=dtsets(idtset)%optstress
   optddk=0;optphon=0;optstrs=0
   if (dtsets(idtset)%rfddk>0.or.dtsets(idtset)%rf2_dkdk>0.or.dtsets(idtset)%rf2_dkde>0) optddk=1
   if (dtsets(idtset)%rfelfd>0.or.dtsets(idtset)%d3e_pert1_elfd>0.or.&
&   dtsets(idtset)%d3e_pert2_elfd>0.or.dtsets(idtset)%d3e_pert3_elfd>0) optddk=1
   if (dtsets(idtset)%rfphon>0.or.dtsets(idtset)%d3e_pert1_phon>0.or.&
&   dtsets(idtset)%d3e_pert2_phon>0.or.dtsets(idtset)%d3e_pert3_phon>0) optphon=1
   if (dtsets(idtset)%rfstrs>0) optstrs=1

   ABI_ALLOCATE(nband,(nkpt*nsppol))
   nband(1:nkpt*nsppol)=dtsets(idtset)%nband(1:nkpt*nsppol)
   mband=maxval(nband(1:nkpt*nsppol))
   dtsets(idtset)%mband=mband
   
!  mpsang=max(maxval(pspheads(1:npsp)%lmax)+1,1) ! Likely problems with the HP compiler
!  n1xccc=maxval(pspheads(1:npsp)%xccc)
   mpsang=1
   n1xccc=pspheads(1)%xccc
   do ii=1,npsp
     mpsang=max(pspheads(ii)%lmax+1,mpsang)
     n1xccc=max(pspheads(ii)%xccc,n1xccc)
   end do

!  Determine the maximum number of projectors, for the set of pseudo atom
   call getdim_nloc(lmnmax,lmnmaxso,lnmax,lnmaxso,dtsets(idtset)%mixalch_orig,&
&   dtsets(idtset)%nimage,npsp,dtsets(idtset)%npspalch,ntypat,dtsets(idtset)%ntypalch,pspheads)

!  Treatment of the effect of using a spin-orbit part
!  Warning : mpspso is different for each dataset.
   mpspso=1
   do ii=1,npsp
     if(nspinor/=1)then
       if(pspheads(ii)%pspso/=0)then
         if(dtsets(idtset)%so_psp(ii)/=0)then
           mpspso=2
         end if
       end if
     end if
   end do
!  In case of no spin-orbit
   if(mpspso==1)then
     mpssoang=mpsang ; lmnmax_eff =lmnmax; lnmax_eff =lnmax
   else ! spin-orbit will be used
     mpssoang=2*mpsang-1 ; lmnmax_eff =lmnmaxso ; lnmax_eff =lnmaxso
   end if
!  lmnmax is not used if the Ylm are not used
   if (useylm==0) lmnmax_eff =lnmax_eff

   ecut     =dtsets(idtset)%ecut
   dilatmx  =dtsets(idtset)%dilatmx
   ecut_eff=ecut*dilatmx**2
   ecutdg_eff=dtsets(idtset)%pawecutdg*dtsets(idtset)%dilatmx**2

!  Compute mgfft,mpw,nfft for this data set
   call mkrdim(dtsets(idtset)%acell_orig(1:3,1),dtsets(idtset)%rprim_orig(1:3,1:3,1),rprimd)
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

   if (usepaw==0) then
     mgfftf=dtsets(idtset)%mgfft;nfftf=dtsets(idtset)%nfft;ngfftf(:)=dtsets(idtset)%ngfft(:)
   else
     mgfftf=dtsets(idtset)%mgfftdg;nfftf=dtsets(idtset)%nfftdg;ngfftf(:)=dtsets(idtset)%ngfftdg(:)
   end if
   response=0
   if(dtsets(idtset)%rfddk/=0  .or. dtsets(idtset)%rf2_dkdk/=0 .or. dtsets(idtset)%rf2_dkde/=0 .or. &
&   dtsets(idtset)%rfphon/=0 .or. dtsets(idtset)%rfelfd/=0 .or. &
&   dtsets(idtset)%rfstrs/=0 .or. dtsets(idtset)%rfuser/=0 .or. &
&   dtsets(idtset)%rfmagn/=0    ) response=1

!  Compute mgfftdiel,npwdiel,nfftdiel for this data set
   if((modulo(iprcel,100)>=20 .and.modulo(iprcel,100) < 71).or. iscf==-1)then
!    Get diecut, and the fft grid to be used for the susceptibility computation
     diecut=abs(dtsets(idtset)%diecut)
     if( dtsets(idtset)%diecut < zero )then
       ecutsus=ecut
     else
       ecutsus= ( sqrt(ecut) *0.5_dp + sqrt(diecut) *0.25_dp )**2
     end if
!    Beware, for the dielectric matrix fftalg=ngfftdiel(7) is default here
     ngfftdiel(1:3)=0 ; ngfftdiel(7)=101 ; ngfftdiel(8:18)=dtsets(idtset)%ngfft(8:18)
     if(iscf==-1)ngfftdiel(7)=102
     ecut_eff=ecutsus*dilatmx**2
     call getng(dtsets(idtset)%boxcutmin,ecut_eff,gmet,k0,me_fft,mgfftdiel,nfftdiel,&
&     ngfftdiel,nproc_fft,nsym,paral_fft,symrel,&
&     use_gpu_cuda=dtsets(idtset)%use_gpu_cuda)
!    Compute the size of the dielectric matrix : npwdiel
     kpt_diel(1:3)=(/ 0.0_dp, 0.0_dp, 0.0_dp /)
     ecut_eff=diecut*dilatmx**2
     call getmpw(ecut_eff,exchn2n3d,gmet,(/1/),kpt_diel,mpi_enregs(idtset),npwdiel,1)
   else
     npwdiel=1 ; mgfftdiel=1 ; nfftdiel=1 ; ngfftdiel(1:8)=1
   end if

!  Special treatment for the value of mqgrid to be fed in memory.F90

   nptsgvec         = 200 ! At present, this has to be chosen once and for all ... 
   if ( dtsets(idtset)%usewvl == 0) then
     call setmqgrid(mqgrid,mqgriddg,ecut_eff,ecutdg_eff,gprimd,nptsgvec,usepaw)
   else
     call setmqgrid(mqgrid,mqgriddg,one,one,gprimd,nptsgvec,usepaw)
   end if
   mqgrid_ff=mqgrid
   if (usepaw==0) mqgrid_vl=mqgrid
   if (usepaw==1) mqgrid_vl=mqgriddg

!  Compute the memory needs for this data set.
   if(response==0)then

     if (dtsets(idtset)%usewvl == 0) then
       mkmem=dtsets(idtset)%mkmem
       mband=maxval(dtsets(idtset)%nband(1:nkpt*nsppol))

       ! Don't perform memory tests if MBPT.
       mem_test = dtsets(idtset)%mem_test 
       if (any(dtsets(idtset)%optdriver == [RUNL_SIGMA, RUNL_SCREENING, RUNL_BSE])) mem_test = 0

       call memory(n1xccc,extrapwf,getcell,idtset,dtsets(idtset)%icoulomb,&
&       intxc,dtsets(idtset)%ionmov,iout,densfor_pred,&
&       iprcel,iscf,jdtset,lmnmax_eff,lnmax_eff,mband,mffmem,dtsets(idtset)%mgfft,mgfftdiel,mgfftf,mkmem,&
&       mpi_enregs(idtset),mpsang,mpssoang,mpw,mqgrid_ff,mqgrid_vl,natom,nband,dtsets(idtset)%nfft,nfftdiel,nfftf,&
&       dtsets(idtset)%ngfft,ngfftdiel,ngfftf,dtsets(idtset)%nimage,nkpt,nloalg,npsp,npulayit,npwdiel,nspden,nspinor,&
&       nsppol,nsym,ntypat,occopt,optforces,mem_test,optstress,pawcpxocc,pawmixdg,&
&       pawnhatxc,pawspnorb,pawstgylm,prtvol,pspheads,dtsets(idtset)%tfkinfunc,&
&       dtsets(idtset)%typat,ucvol,usepaw,useylm,use_gpu_cuda,xclevel)
     else if( dtsets(idtset)%usepaw==0) then
       if (mpi_enregs(idtset)%me == 0) then
         call wvl_memory(dtsets(idtset), idtset, mpi_enregs(idtset), npsp, 1, pspheads)
       end if
     end if

   else
!    Compute the value of cplex, for which one needs symrec
     ABI_ALLOCATE(symq,(4,2,nsym))
     ABI_ALLOCATE(symrec,(3,3,nsym))
     do isym=1,nsym
       call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
     end do
     call littlegroup_q(nsym,qphon,symq,symrec,dtsets(idtset)%symafm,timrev)
     cplex=2-timrev
     ABI_DEALLOCATE(symq)
     ABI_DEALLOCATE(symrec)
     mkmems(1)=dtsets(idtset)%mkmem
     mkmems(2)=dtsets(idtset)%mkqmem
     mkmems(3)=dtsets(idtset)%mk1mem

     mem_test = dtsets(idtset)%mem_test 

     call memorf(cplex,n1xccc,getcell,idtset,intxc,iout,iprcel,&
&     iscf,jdtset,lmnmax_eff,lnmax_eff,mband,mffmem,dtsets(idtset)%mgfft,&
&     mkmems,mpi_enregs(idtset),mpsang,mpssoang,mpw,mqgrid_ff,natom,nband,dtsets(idtset)%nfft,&
&     dtsets(idtset)%ngfft,nkpt,nloalg,nspden,nspinor,nsppol,nsym,&
&     ntypat,occopt,optddk,optphon,mem_test,optstrs,prtvol,useylm,use_gpu_cuda,xclevel)
   end if

!  Deallocate temporary arrays (when they will really be temporary !)
   ABI_DEALLOCATE(nband)
   ABI_DEALLOCATE(symrel)

 end do ! idtset

end subroutine memory_eval
!!***
