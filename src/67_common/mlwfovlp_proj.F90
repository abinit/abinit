!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp_proj
!! NAME
!! mlwfovlp_proj
!!
!! FUNCTION
!! Routine which computes projection A_{mn}(k)
!! for Wannier code (www.wannier.org f90 version).
!!
!! COPYRIGHT
!! Copyright (C) 2005-2016 ABINIT group (BAmadon,FJollet,TRangel,drh)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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
!!  mpi_enreg=informations about MPI parallelization
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
!!  spin = just used for nsppol>1 ; 0 both, 1 just spin up, 2 just spin down
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
!!      mlwfovlp
!!
!! CHILDREN
!!      mlwfovlp_radial,mlwfovlp_ylmfac,wrtout,ylm_cmplx
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine mlwfovlp_proj(A_matrix,band_in,cg,cprj,dtset,gprimd,just_augmentation,kg,&
&lproj,max_num_bands,mband,mkmem,mpi_enreg,mpw,mwan,natom,nattyp,&
&nkpt,npwarr,nspinor,&
&nsppol,ntypat,num_bands,nwan,pawtab,proj_l,proj_m,proj_radial,&
&proj_site,proj_x,proj_z,proj_zona,psps,spin,ucvol)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wannier90
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_pawtab,  only : pawtab_type
 use m_pawcprj, only : pawcprj_type
 use m_paw_sphharm, only : ylm_cmplx

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mlwfovlp_proj'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_67_common, except_this_one => mlwfovlp_proj
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 complex(dpc),parameter :: c1=(1._dp,0._dp)
 integer,intent(in) :: lproj,max_num_bands,mband,mkmem,mpw,mwan,natom,nkpt,nspinor,nsppol
 integer,intent(in) :: ntypat,spin
 type(MPI_type),intent(inout) :: mpi_enreg
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
   MSG_ERROR(message)
 end if

 write(message, '(a,a)' )ch10,&
& '** mlwfovlp_proj:  compute A_matrix of initial guess for wannier functions'
 call wrtout(std_out,message,'COLL')

!Initialize to 0.d0
 A_matrix(:,:,:,:)=cmplx(0.d0,0.d0)


!
!End of preliminarities
!

!********************* Write Random projectors
 if(lproj==1) then
   idum=123456
!  Compute random projections
   ABI_ALLOCATE(amn,(2,mband,mwan,nkpt,nsppol))
   amn=zero
   do isppol=1,nsppol
     if(spin.ne.0 .and. spin.ne.isppol) cycle
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
     if(spin.ne.0 .and. spin.ne.isppol) cycle
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
&               '  Value of jband is above num_bands ',ch10,&
&               '  Action : contact Abinit group'
               MSG_ERROR(message)
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
     if(spin.ne.0 .and. spin.ne.isppol) cycle
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
     if(spin.ne.0 .and. spin.ne.isppol) cycle
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
     if(spin.eq.0 .or. spin.eq.isppol) then
!      this has to be done this way because the variable icg changes at the end of the
!      cycle. We cannot just skip the hole cycle.
       ABI_ALLOCATE(gf,(mpw,nproj(isppol)))
       ABI_ALLOCATE(gft_lm,(lmax2(isppol)))
       ABI_ALLOCATE(gsum2,(nproj(isppol)))
       ABI_ALLOCATE(kpg2,(mpw))
       ABI_ALLOCATE(radial,(lmax2(isppol)))
       ABI_ALLOCATE(ylmcp,(lmax2(isppol)))
     end if
!    
     ikg=0
     do ikpt=1, nkpt
!      
!      MPI: cycle over kpts not treated by this node
!      
       if (ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-rank)/=0) CYCLE
!      
       if(spin.eq.0 .or. spin.eq.isppol) then
         write(message, '(a,i6,a,2i6)' ) &
&         '   processor',rank,' will compute k-point,spin=',ikpt,isppol
         call wrtout(std_out,  message,'COLL')
       end if
!      
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
         if(spin.eq.0 .or. spin.eq.isppol) then
!          !
           do iproj=1,nproj(isppol)
!            
!            In PAW, we can use proj_radial > 4 to indicate that we just 
!            want the in-sphere contribution
!            
             if( psps%usepaw==1) then
               if( just_augmentation(iproj,isppol)) cycle
             end if
!            
!            obtain radial part
             call mlwfovlp_radial(proj_zona(iproj,isppol),lmax(isppol),lmax2(isppol)&
&             ,radial,proj_radial(iproj,isppol),kpg2(ipw))
!            
!            scale complex representation of projector orbital with radial functions
!            of appropriate l
             gft_lm(:)=radial(:)*ylmc_fac(1:lmax2(isppol),iproj,isppol)
!            
!            complex structure factor for projector orbital position
             arg = ( kpg(1)*proj_site(1,iproj,isppol) + &
&             kpg(2)*proj_site(2,iproj,isppol) + &
&             kpg(3)*proj_site(3,iproj,isppol) ) * 2*pi
             cstr_fact = cmplx(cos(arg), -sin(arg) )
!            
!            obtain guiding functions
             gf(ipw,iproj)=cstr_fact*dot_product(ylmcp,gft_lm)
!            
             gsum2(iproj)=gsum2(iproj)+real(gf(ipw,iproj))**2+aimag(gf(ipw,iproj))**2
           end do !iproj
         end if !spin
       end do !ipw
!      
       if(spin.eq.0 .or. spin.eq.isppol) then
         do iproj=1,nproj(isppol)
!          
!          In PAW, we can use proj_radial > 4 to indicate that we just 
!          want the in-sphere contribution
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
!        Guiding functions are computed.
!        compute overlaps of gaussian projectors and wave functions
         do iproj=1,nproj(isppol)
!          
!          In PAW, we can use proj_radial > 4 to indicate that we just 
!          want the in-sphere contribution
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
!                  The case of spinors is tricky, we have nproj =  nwan/2
!                  so we project to spin up and spin down separately, to have at 
!                  the end an amn matrix with nwan projections.
!                  
!                  idx=ipw*nspinor - (nspinor-ispinor)
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
       end if !spin==isppol
       icg=icg+npw_k*nspinor*mband
       ikg=ikg+npw_k
     end do !ikpt
!    Deallocations
     if(spin.eq.0 .or. spin.eq.isppol) then
       ABI_DEALLOCATE(gf)
       ABI_DEALLOCATE(gft_lm)
       ABI_DEALLOCATE(gsum2)
       ABI_DEALLOCATE(kpg2)
       ABI_DEALLOCATE(radial)
       ABI_DEALLOCATE(ylmcp)
     end if
   end do !isppol
!  
!  if(isppol==1) then
!  norm_error_bar=sqrt(norm_error_bar/real(nkpt*(nwan(1)),dp))
!  else
!  if(spin==0)    norm_error_bar=sqrt(norm_error_bar/real(nkpt*(nwan(1)+nwan(2)),dp))
!  if(spin==1)    norm_error_bar=sqrt(norm_error_bar/real(nkpt*nwan(1),dp))
!  if(spin==2)    norm_error_bar=sqrt(norm_error_bar/real(nkpt*nwan(2),dp))
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
     MSG_ERROR("natprjn/=1")
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
     MSG_ERROR("Aborting now")
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
                               MSG_ERROR("Aborting now")
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
