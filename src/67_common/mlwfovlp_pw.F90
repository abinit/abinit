!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp_pw
!! NAME
!! mlwfovlp_pw
!!
!! FUNCTION
!! Routine which computes PW part of overlap M_{mn}(k,b) 
!! for Wannier code (www.wannier.org f90 version).
!!
!! COPYRIGHT
!! Copyright (C) 2005-2018 ABINIT group (BAmadon,FJollet,T Rangel)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!  g1(3,nkpt,nntot) = G vector shift which is necessary to obtain k1+b
!!  iwav(mband,nkpt,nsppol): shift for pw components in cg.
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points treated by this node.
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum dimensioned size of npw.
!!  nfft=(effective) number of FFT grid points (for this processor) (see NOTES at beginning of scfcv)
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ovikp(nkpt,nntot)= gives  nntot value of k2 (in the BZ) for each k1  (k2=k1+b mod(G))
!!  seed_name= seed_name of files containing cg for all k-points to be used with MPI
!!  spin = just used for nsppol>1 ; 0 both, 1 just spin up, 2 just spin down
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
!!      mlwfovlp
!!
!! CHILDREN
!!      wrtout,xmpi_barrier,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mlwfovlp_pw(cg,cm1,g1,iwav,kg,mband,mkmem,mpi_enreg,mpw,nfft,ngfft,nkpt,nntot,&
&  npwarr,nspinor,nsppol,ovikp,seed_name,spin)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mlwfovlp_pw'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,mpw,nfft,nkpt,nntot
 integer,intent(in) :: nspinor,nsppol,spin
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
         MSG_ERROR("Aborting now")
       end if
     end do
     npw_k=npwarr(ikpt)
!    write(std_out,*) ikpt,npw_k,nfft
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     do ig=1,npw_k
       if(ig.gt.mpw) then
         write(std_out,*)"error ig",ig,"greater than mpw ",mpw
         MSG_ERROR("Aborting now")
       end if
       if(indpwk(ikpt,ig)/=0) then
         write(std_out,*) "error, indpwk is overwritten"
         write(std_out,*) ikpt,ig,indpwk(ikpt,ig)
         MSG_ERROR("Aborting now")
       end if
       ig1=modulo(kg_k(1,ig),n1)
       ig2=modulo(kg_k(2,ig),n2)
       ig3=modulo(kg_k(3,ig),n3)
       indpwk(ikpt,ig)=ig1+1+n1*(ig2+n2*ig3)
       npoint=indpwk(ikpt,ig)
       if(npoint.gt.nfft) then
         MSG_ERROR("error npoint")
       end if
!      write(std_out,*) ikpt,ig,npoint,invpwk(ikpt,npoint)
       if(invpwk(ikpt,npoint)/=0) then
         write(std_out,*) "error, invpwk is overwritten"
         write(std_out,*) ikpt,ig,npoint,invpwk(ikpt,npoint)
         MSG_ERROR("Aborting now")
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



 
!
!Deallocate unused variables
!
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
   if(spin.ne.0 .and. spin.ne.isppol) cycle
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
           MSG_ERROR(message)
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
               MSG_ERROR("Aborting now")
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
               MSG_ERROR("Aborting now")
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
           end do
         end do
       end do
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
