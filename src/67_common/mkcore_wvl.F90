!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkcore_wvl
!! NAME
!!  mkcore_wvl
!!
!! FUNCTION
!! Optionally compute
!!  (1) core electron density throughout unit cell, or
!!  (2) contribution to Exc gradient wrt xred, or
!!  (3) contribution to stress tensor, or (response function code)
!!  (4) contribution to frozen-wavefunction part of
!!      the dynamical matrix (part 2)
!!
!! COPYRIGHT
!!  Copyright (C) 2011-2017 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Based on mkcore.F90
!!
!! PARENTS
!!      setvtr
!!
!! CHILDREN
!!      ext_buffers,ind_positions,metric,mkcore_inner,mkdenpos,pawrad_free
!!      pawrad_init,strconv,timab,wrtout,xcart2xred,xmpi_sum,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

subroutine mkcore_wvl(atindx1,corstr,dyfrx2,geocode,grxc,h,natom,&
& nattyp,nfft,nscatterarr,nspden,ntypat,n1,n1i,n2,n2i,n3,n3pi,&
&n3xccc,option,pawrad,pawtab,psppar,rprimd,ucvol,&
& vxc,xccc3d,xred,mpi_comm_wvl)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_splines

 use m_pawrad,  only : pawrad_type, pawrad_init, pawrad_free
 use m_pawtab,  only : pawtab_type
 use m_xmpi,    only : xmpi_comm_size,xmpi_sum

#if defined HAVE_BIGDFT
  use BigDFT_API, only: ext_buffers,ind_positions
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkcore_wvl'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_41_xc_lowlevel
 use interfaces_67_common, except_this_one => mkcore_wvl
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat,nfft,nspden
 integer,intent(in) ::n1,n2,n3,n1i,n2i,n3pi,n3xccc,option
 integer,intent(in),optional:: mpi_comm_wvl
 real(dp),intent(in) :: h(3),ucvol
 character(1),intent(in)::geocode
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat)
 integer,intent(in) :: nscatterarr(4)
 real(dp),intent(in) :: psppar(0:4,0:6,ntypat),rprimd(3,3)
 real(dp),intent(in)::xred(3,natom)
 real(dp),intent(in)::vxc(nfft,nspden)
 real(dp),intent(out)::xccc3d(n3xccc)
 real(dp),intent(out) :: corstr(6),dyfrx2(3,3,natom),grxc(3,natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)
 type(pawrad_type),intent(in) :: pawrad(ntypat)

!Local variables ------------------------------
#if defined HAVE_BIGDFT
!scalars
!buffer to be added at the end of the last dimension of an array to control bounds_check
 integer :: i1,i2,i3,iat,iatm,iatom,iatom_tot
 integer :: itypat,iwarn
 integer :: nproc_wvl=1
 integer :: ier,iex,iey,iez,isx,isy,isz,ind,i3s
 integer :: j1,j2,j3,msz
 integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3
 integer :: ncmax,nfgd,nfgd_r0,opt_mkdenpos,shift
 real(dp) :: cutoff,factor,grxc1,grxc2,grxc3
 real(dp) :: rloc,rr2,rx,ry,rz
 real(dp) :: rshp,r2shp,rshpm1,strdia,t1,t2,t3
 real(dp) :: xx,yy,zz,ucvol_
 character(len=500) :: message
 logical :: perx,pery,perz,gox,goy,goz
 type(pawrad_type)::core_mesh
!arrays
 real(dp) :: hh(3) !fine grid spacing for wavelets
 real(dp) :: gmet(3,3),gprimd(3,3),rcart(3),rmet(3,3),tsec(2),xcart(3,natom)
 real(dp) :: corfra(3,3)
!allocatable arrays
 integer,allocatable :: ifftsph_tmp(:)!,iindex(:)
! real(dp),allocatable:: raux(:)
! real(dp),pointer ::raux2(:),raux3(:)
 real(dp),pointer::rred(:,:)
 real(dp),allocatable:: rr(:)
#endif

! *************************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT

 if(nspden >1) then
   write(message, '(a)')'mkcore_wvl: this is not yet generalized to npsden>1'
   MSG_ERROR(message)
 end if
 if(option>4 .or. option<1 )then
   write(message,'(a,a,a,a,a,a,i6)') ch10,&
&   ' mkcore_wvl: BUG -',ch10,&
&   '  The argument option should be between 1 and 4,',ch10,&
&   '  however, option=',option
   MSG_BUG(message)
 end if
 if(nfft .ne. n3xccc)then
   write(message,'(a,a,a,a,a,a,2i6)') ch10,&
&   ' mkcore_wvl: BUG -',ch10,&
&   '  nfft and n3xccc should be equal,',ch10,&
&   '  however, nfft and n3xccc=',nfft,n3xccc
   MSG_BUG(message)
 end if

!mpi
 if(present(mpi_comm_wvl)) nproc_wvl=xmpi_comm_size(mpi_comm_wvl)
 i3s  =nscatterarr(3)+1-nscatterarr(4)
 shift=n1i*n2i*nscatterarr(4)

 if (option==1) then
!  Zero out array to permit accumulation over atom types below:
   xccc3d(:)=zero
 else if (option==2) then
!  Zero out gradient of Exc array
   grxc(:,:)=zero
 else if (option==3) then
!  Zero out locally defined stress array
   corfra(:,:)=zero
   strdia=zero
 else if (option==4) then
!  Zero out fr-wf part of the dynamical matrix
   dyfrx2(:,:,:)=zero
 else
   write(message, '(a,a,a,a)' ) ch10,&
&   ' mkcore_wvl: BUG -',ch10,&
&   '  Can''t be here ! (bad option).'
   MSG_BUG(message)
 end if

!ENDDEBUG

 write(message,'(a,a)') ch10,&
& ' mkcore_wvl: Compute core density'
 call wrtout(std_out,message,'COLL')

!Fine grid
 hh(:)=0.5d0*h(:)

!Compute xcart from xred
 call xred2xcart(natom,rprimd,xcart,xred)

!Compute metric tensors and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol_)
!correct ucvol for wavelets case is given as an input:
!ucvol_local = product(wvl%den%denspot%dpbox%hgrids) * real(product(wvl%den%denspot%dpbox%ndims), dp)
!ucvol_local = (half * dtset%wvl_hgrid) ** 3 * ngfft(1)*ngfft(2)*ngfft(3)

!Conditions for periodicity in the three directions
 perx=(geocode /= 'F')
 pery=(geocode == 'P')
 perz=(geocode /= 'F')

!Compute values of external buffers
 call ext_buffers(perx,nbl1,nbr1) 
 call ext_buffers(pery,nbl2,nbr2)
 call ext_buffers(perz,nbl3,nbr3)


 iatm=0
!Big loop on atom types
 do itypat=1,ntypat
!  
   rloc=psppar(0,0,itypat)
   cutoff=10.d0*rloc

!  Set radius size:
   rshp=pawtab(itypat)%rcore
   r2shp=1.0000001_dp*rshp**2
   rshpm1=one/rshp

!  allocate arrays
   if (n3pi > 0) then
!    ncmax=1+int(1.1_dp*nfft*four_pi/(three*ucvol)*rshp**3)
!    ncmax=1+int(1.1_dp*nfft*four_pi/(three*ucvol)*rshp**3)
!    1+int(1.1* factors are included just for cautioness
     ncmax=1+int(1.1d0*((rshp/hh(1))*(rshp/hh(2))*pi))
   else
     ncmax=1
   end if
!  
   ABI_ALLOCATE(ifftsph_tmp,(ncmax))
!  ABI_ALLOCATE(iindex,(ncmax))
   ABI_ALLOCATE(rr,(ncmax))
!  ABI_ALLOCATE(raux,(ncmax))
!  nullify(raux2,raux3)
   nullify(rred)
   if(option>2) then
     ABI_ALLOCATE(rred,(3,ncmax))
!    ABI_ALLOCATE(raux2,(ncmax))
   end if
!  if(option==4) then
!  ABI_ALLOCATE(raux3,(ncmax))
!  end if

!  Create mesh_core object
!  since core_mesh_size can be bigger than pawrad%mesh_size, 
   msz=pawtab(itypat)%core_mesh_size
   call pawrad_init(core_mesh,mesh_size=msz,mesh_type=pawrad(itypat)%mesh_type,&
&   rstep=pawrad(itypat)%rstep,lstep=pawrad(itypat)%lstep)

!  Big loop on atoms  
   do iat=1,nattyp(itypat)
     iatm=iatm+1;iatom=atindx1(iatm)
     iatom_tot=iatom; 
!    
     if(option==2) then
       grxc1=zero
       grxc2=zero
       grxc3=zero
     end if

!    Define a "box" around each atom
     rx=xcart(1,iatom_tot)
     ry=xcart(2,iatom_tot)
     rz=xcart(3,iatom_tot)
!    
     isx=floor((rx-cutoff)/hh(1))
     isy=floor((ry-cutoff)/hh(2))
     isz=floor((rz-cutoff)/hh(3))
     
     iex=ceiling((rx+cutoff)/hh(1))
     iey=ceiling((ry+cutoff)/hh(2))
     iez=ceiling((rz+cutoff)/hh(3))
!    
     do i3=isz,iez
       zz=real(i3,kind=8)*hh(3)-rz
       call ind_positions(perz,i3,n3,j3,goz)
       j3=j3+nbl3+1

!      Initialize counters
       nfgd=0
       nfgd_r0=0
!      
       do i2=isy,iey
         yy=real(i2,kind=8)*hh(2)-ry
         call ind_positions(pery,i2,n2,j2,goy)
!        
         do i1=isx,iex
           xx=real(i1,kind=8)*hh(1)-rx
           call ind_positions(perx,i1,n1,j1,gox)
           rr2=xx**2+yy**2+zz**2
           if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
!            
             if(rr2<=r2shp) then
               if(rr2>tol5) then
                 ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s)*n1i*n2i
                 nfgd=nfgd+1
                 rcart(1)=xx; rcart(2)=yy; rcart(3)=zz
                 rr(nfgd)=(rr2)**0.5
                 ifftsph_tmp(nfgd)=shift+ind
                 if(option>1) then 
                   call xcart2xred(1,rprimd,rcart,rred(:,nfgd))
                 end if
!                DEBUG     
!                write(itmp,'(i10,3(f13.7,x))')ind,x+rx,y+ry,z+rz
!                write(itmp,'(6(f13.7,x))')rcart,rred(:,nfgd)
!                ENDDEBUG
               elseif( option == 4) then
!                We save r=0 vectors only for option==4:
!                for other options this is ignored
                 ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s)*n1i*n2i

!                We reuse the same variable "ifftshp_tmp", 
!                but we start from the higher index
                 nfgd_r0=nfgd_r0+1
                 ifftsph_tmp(ncmax-nfgd_r0+1)=shift+ind

               end if !rr2>tol5
             end if !rr2<r2shp
           end if !j3..
         end do !i1
       end do !i2

!      All of the following  could be done inside or outside the loops (i2,i1,i3)
!      Outside the loops: the memory consuption increases.
!      Inside the inner loop: the time of calculation increases.
!      Here, I chose to do it here, somewhere in the middle.

       if(option .ne.4 ) then
         if(nfgd==0)      cycle
       else
         if(nfgd==0 .and. nfgd_r0==0) cycle
       end if
       call mkcore_inner(corfra,core_mesh,dyfrx2,&
&       grxc1,grxc2,grxc3,ifftsph_tmp,msz,&
&       natom,ncmax,nfft,nfgd,nfgd_r0,nspden,n3xccc,option,pawtab(itypat),&
&       rmet,rprimd,rr,rred,rshpm1,strdia,ucvol,vxc,xccc3d)
     end do !i3

     if(option==2) then
       factor=(ucvol/real(nfft,dp))/rshp
       grxc(1,iatom)=grxc1*factor
       grxc(2,iatom)=grxc2*factor
       grxc(3,iatom)=grxc3*factor

       if( nproc_wvl>1 ) then
         call timab(539,1,tsec)
         call xmpi_sum(grxc1,mpi_comm_wvl,ier)
         call xmpi_sum(grxc2,mpi_comm_wvl,ier)
         call xmpi_sum(grxc3,mpi_comm_wvl,ier)
         call timab(539,2,tsec)
       end if
     end if

   end do !iat

!  Deallocate
   call pawrad_free(core_mesh)
   ABI_DEALLOCATE(ifftsph_tmp)
!  ABI_DEALLOCATE(iindex)
   ABI_DEALLOCATE(rr)
!  ABI_DEALLOCATE(raux)
   if(option>2) then
     if(associated(rred)) then
       ABI_DEALLOCATE(rred)
     end if
!    if(associated(raux2))ABI_DEALLOCATE(raux2)
   end if
!  if(option>4) then
!  if(associated(raux3))   ABI_DEALLOCATE(raux3)
!  end if
 end do !itypat



 if (option==2) then
!  Apply rmet as needed to get reduced coordinate gradients
   do iatom=1,natom
     t1=grxc(1,iatom)
     t2=grxc(2,iatom)
     t3=grxc(3,iatom)
     grxc(:,iatom)=rmet(:,1)*t1+rmet(:,2)*t2+rmet(:,3)*t3
   end do

 elseif( option==3) then

!  Transform stress tensor from full storage mode to symmetric storage mode
   corstr(1)=corfra(1,1)
   corstr(2)=corfra(2,2)
   corstr(3)=corfra(3,3)
   corstr(4)=corfra(3,2)
   corstr(5)=corfra(3,1)
   corstr(6)=corfra(2,1)

!  Transform stress tensor from reduced coordinates to cartesian coordinates
   call strconv(corstr,rprimd,corstr)


!  Compute diagonal contribution to stress tensor (need input xccc3d)
!  strdia = (1/N) Sum(r) [mu_xc_avg(r) * rho_core(r)]
   strdia=zero
   do i3=1,n3pi
     do i2=1,n2i
       do i1=1,n1i
         ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
         strdia=strdia+vxc(ind,1)*xccc3d(ind)
!        write(17,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind),maxdiff
       end do
     end do
   end do
   strdia=strdia/real(nfft,dp)

!  Add diagonal term to stress tensor
   corstr(1)=corstr(1)+strdia
   corstr(2)=corstr(2)+strdia
   corstr(3)=corstr(3)+strdia
 end if

 if(nproc_wvl > 1) then
   call timab(539,1,tsec)
   if(option==3) then
     call xmpi_sum(corstr,mpi_comm_wvl,ier)
   end if
   if(option==2) then
     call xmpi_sum(grxc,mpi_comm_wvl,ier)
   end if
   call timab(539,2,tsec)
 end if

!Make xccc3d positive to avoid numerical instabilities in V_xc
 iwarn=0 ; opt_mkdenpos=0
 call mkdenpos(iwarn, n3xccc, nspden, opt_mkdenpos, xccc3d, tol20 )

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) natom,ntypat,nfft,nspden,n1,n2,n3,n1i,n2i,n3pi,n3xccc,option,&
& mpi_comm_wvl,h(1),ucvol,geocode,atindx1(1),nattyp(1),nscatterarr(1),psppar(1,1,1),rprimd(1,1),&
& xred(1,1),vxc(1,1),xccc3d(1),corstr(1),dyfrx2(1,1,1),grxc(1,1),pawtab(1)%mesh_size,pawrad(1)%mesh_size
#endif

 DBG_EXIT("COLL")

end subroutine mkcore_wvl
!!***
