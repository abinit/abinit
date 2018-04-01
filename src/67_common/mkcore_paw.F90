!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkcore_paw
!! NAME
!!  mkcore_paw
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2018 ABINIT group (TRangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      ind_positions_,mkcore_inner,mkgrid_fft,pawrad_free,pawrad_init
!!      ptabs_fourdp,strconv,timab,wrtout,xcart2xred,xmpi_sum,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkcore_paw(atindx1,corstr,dyfrx2,grxc,icoulomb,natom,mpi_enreg,&
& nattyp,nfft,ngfft,nspden,ntypat,n3xccc,option,pawrad,pawtab,psppar,rprimd,&
& ucvol,vxc,xccc3d,xred)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_time,     only : timab
 use m_geometry, only : xcart2xred, xred2xcart
 use m_pawrad,   only : pawrad_type, pawrad_init, pawrad_free
 use m_pawtab,   only : pawtab_type
 use m_mpinfo,   only : ptabs_fourdp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkcore_paw'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_67_common, except_this_one => mkcore_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icoulomb
 integer,intent(in) :: natom,ntypat,nfft,nspden,n3xccc,option
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat)
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: psppar(0:4,0:6,ntypat)
 real(dp),intent(in) :: rprimd(3,3),vxc(nfft,nspden)
 real(dp),intent(in) :: ucvol
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(inout) :: xccc3d(nfft)
 real(dp),intent(out) :: corstr(6),dyfrx2(3,3,natom),grxc(3,natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)
 type(pawrad_type),intent(in) :: pawrad(ntypat)

!Local variables-------------------------------
 integer :: iat,iatm,iatom,iatom_tot
 integer :: ier,iex,iey,iez,ind
 integer :: isx,isy,isz,itypat
 integer :: i1,i2,i3,i3loc
 integer :: j1,j2,j3,msz
 integer :: ncmax,nfgd,nfgd_r0
 integer :: me,me_fft,nfftot,nproc,nproc_fft,nu
 integer :: n1,n2,n3,n3d
 real(dp) :: cutoff,factor,grxc1,grxc2,grxc3
 real(dp) :: rloc,rr2,rshp,rshpm1,rx,ry,rz
 real(dp) :: r2shp,strdia
 character(len=500) :: message
 character(len=1) :: geocode
 logical :: perx,pery,perz,gox,goy,goz
 type(pawrad_type)::core_mesh
!arrays
 integer,allocatable :: ifftsph_tmp(:)
 real(dp) :: corfra(3,3)
 real(dp) :: hh(3) !fine grid spacing
 real(dp) :: rmet(3,3),tsec(2),xcart(3,natom)
 real(dp),allocatable :: gridcart(:,:),rr(:),rred(:,:)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)

! *************************************************************************

 DBG_ENTER("COLL")

 if(nspden >1) then
   write(message, '(a)')'mkcore_paw: this is not yet generalized to npsden>1'
   MSG_ERROR(message)
 end if

 geocode='P'
 if (icoulomb==1) geocode='F'
 if (icoulomb==2) geocode='S'

!Compute metric tensor in real space rmet
 do nu=1,3
   rmet(:,nu)=rprimd(1,:)*rprimd(1,nu)+rprimd(2,:)*rprimd(2,nu)+&
&   rprimd(3,:)*rprimd(3,nu)
 end do

!MPI
 nproc =xmpi_comm_size(mpi_enreg%comm_fft); nproc_fft= ngfft(10)
 me    =xmpi_comm_rank(mpi_enreg%comm_fft);    me_fft= ngfft(11)

 if(me /= me_fft .or. nproc /= nproc_fft) then
   MSG_BUG("mkcore_paw: comm_size or comm_rank not equal to the corresponding values in ngfft")
 end if

 n1 = ngfft(1)
 n2 = ngfft(2)
 n3 = ngfft(3)
 n3d = ngfft(13)
 if(nproc==1) n3d=n3

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Store xcart for each atom
 call xred2xcart(natom, rprimd, xcart, xred)

!Store cartesian coordinates for each grid points
 ABI_ALLOCATE(gridcart,(3, nfft))
 call mkgrid_fft(ffti3_local,fftn3_distrib,gridcart,nfft,ngfft,rprimd)

!definition of the grid spacings
 hh(1) = rprimd(1,1)/(ngfft(1))
 hh(2) = rprimd(2,2)/(ngfft(2))
 hh(3) = rprimd(3,3)/(ngfft(3))

 if(nfft .ne. n3xccc)then
   write(message,'(a,a,a,a,a,a,2i6)') ch10,&
&   ' mkcore_paw: BUG -',ch10,&
&   '  nfft and n3xccc should be equal,',ch10,&
&   '  however, nfft and n3xccc=',nfft,n3xccc
   MSG_BUG(message)
 end if

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
&   ' mkcore_paw: BUG -',ch10,&
&   '  Can''t be here ! (bad option).'
   MSG_BUG(message)
 end if

 write(message,'(a,a)') ch10,&
& ' mkcore_paw: Compute core density'
 call wrtout(std_out,message,'COLL')

!conditions for periodicity in the three directions
 perx=(geocode /= 'F')
 pery=(geocode == 'P')
 perz=(geocode /= 'F')

 iatm=0
!Big loop on atom types
 do itypat=1,ntypat

   rloc=psppar(0,0,itypat)
   cutoff=10.d0*rloc

!  Set radius size:
   rshp=pawtab(itypat)%rcore
   r2shp=1.0000001_dp*rshp**2
   rshpm1=one/rshp

!  allocate arrays
!  ncmax=1+int(1.1_dp*nfft*four_pi/(three*ucvol)*rshp**3)
!  1+int(1.1* factors are included just for cautioness
   ncmax=1+int(1.1d0*((rshp/hh(1))*(rshp/hh(2))*pi))

   ABI_ALLOCATE(ifftsph_tmp,(ncmax))
   ABI_ALLOCATE(rr,(ncmax))
   if(option>1) then
     ABI_ALLOCATE(rred,(3,ncmax))
   else
     ABI_ALLOCATE(rred,(0,0))
   end if

!  Create mesh_core object
!  since core_mesh_size can be bigger than pawrad%mesh_size,
   msz=pawtab(itypat)%core_mesh_size
   call pawrad_init(core_mesh,mesh_size=msz,mesh_type=pawrad(itypat)%mesh_type,&
&   rstep=pawrad(itypat)%rstep,lstep=pawrad(itypat)%lstep)

!  Big loop on atoms
   do iat=1,nattyp(itypat)
     iatm=iatm+1;iatom=atindx1(iatm)
     iatom_tot=iatom; !if (mpi_enreg%nproc_atom>1) iatom_tot=mpi_enreg%atom_indx(iatom)

     if(option==2) then
       grxc1=zero
       grxc2=zero
       grxc3=zero
     end if

!    Define a "box" around each atom
     rx=xcart(1,iatom)
     ry=xcart(2,iatom)
     rz=xcart(3,iatom)

     isx=floor((rx-cutoff)/hh(1))
     isy=floor((ry-cutoff)/hh(2))
     isz=floor((rz-cutoff)/hh(3))
     iex=ceiling((rx+cutoff)/hh(1))
     iey=ceiling((ry+cutoff)/hh(2))
     iez=ceiling((rz+cutoff)/hh(3))

     do i3=isz,iez
       call ind_positions_(perz,i3,n3,j3,goz)

       if(fftn3_distrib(j3)==me_fft) then
         i3loc=ffti3_local(j3)

!        Initialize counters
         nfgd=0
         nfgd_r0=0

         do i2=isy,iey
           call ind_positions_(pery,i2,n2,j2,goy)
           do i1=isx,iex
             call ind_positions_(perx,i1,n1,j1,gox)
!            r2=x**2+y**2+z**2
             if (goz  .and. goy  .and. gox) then
               ind=j1+(j2-1)*n1+(i3loc-1)*n1*n2
               rr2=(gridcart(1,ind)-rx)**2+(gridcart(2,ind)-ry)**2+(gridcart(3,ind)-rz)**2

               if(rr2<=r2shp) then
                 if(rr2>tol5) then
                   nfgd=nfgd+1
                   rr(nfgd)=sqrt(rr2)
                   ifftsph_tmp(nfgd)=ind
                   if(option>1) then
                     call xcart2xred(1,rprimd,gridcart(:,ind),rred(:,nfgd))
                   end if
                 elseif (option==4) then
!                  We save r=0 vectors only for option==4:
!                  for other options this is ignored

!                  We reuse the same variable "ifftshp_tmp",
!                  but we start from the higher index
                   nfgd_r0=nfgd_r0+1
                   ifftsph_tmp(ncmax-nfgd_r0+1)=ind

                 end if !rr2>tol5
               end if !rr2<r2shp
             end if !gox..
           end do !i1
         end do !i2

!        All of the following  could be done inside or outside the loops (i2,i1,i3)
!        Outside the loops: the memory consuption increases.
!        Inside the inner loop: the time of calculation increases.
!        Here, I choose to do it here, somewhere in the middle.
         if (option/=4.and.nfgd==0) cycle
         if (option==4.and.nfgd==0.and.nfgd_r0==0) cycle
         call mkcore_inner(corfra,core_mesh,dyfrx2,&
&         grxc1,grxc2,grxc3,ifftsph_tmp,msz,&
&         natom,ncmax,nfft,nfgd,nfgd_r0,nspden,n3xccc,option,pawtab(itypat),&
&         rmet,rr,strdia,vxc,xccc3d,rred=rred)

       end if !parallel fftn3
     end do !i3

     if(option==2) then
       nfftot=product(ngfft(1:3))
       factor=(ucvol/real(nfftot,dp)) !/rshp
       grxc(1,iatom)=grxc1*factor
       grxc(2,iatom)=grxc2*factor
       grxc(3,iatom)=grxc3*factor

       if(nproc_fft > 1) then
         call timab(539,1,tsec)
         call xmpi_sum(grxc1,mpi_enreg%comm_fft,ier)
         call xmpi_sum(grxc2,mpi_enreg%comm_fft,ier)
         call xmpi_sum(grxc3,mpi_enreg%comm_fft,ier)
         call timab(539,2,tsec)
       end if

     end if
   end do !iatom

!  Deallocate
   call pawrad_free(core_mesh)
   ABI_DEALLOCATE(ifftsph_tmp)
   ABI_DEALLOCATE(rr)
   ABI_DEALLOCATE(rred)

 end do !itypat

 if (option==2) then
!  Apply rmet as needed to get reduced coordinate gradients
!   do iatom=1,natom
!     t1=grxc(1,iatom)
!     t2=grxc(2,iatom)
!     t3=grxc(3,iatom)
!     grxc(:,iatom)=rmet(:,1)*t1+rmet(:,2)*t2+rmet(:,3)*t3
!!    grxc(:,iatom)=rprimd(1,:)*t1+rprimd(2,:)*t2+rprimd(3,:)*t3
!   end do

 else if (option==3) then

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
   do i3=1,n3d
     do i2=1,n2
       do i1=1,n1
         ind=i1+(i2-1)*n1+(i3-1)*n1*n2
         strdia=strdia+vxc(ind,1)*xccc3d(ind)
       end do
     end do
   end do
   strdia=strdia/real(nfft,dp)
!  Add diagonal term to stress tensor
   corstr(1)=corstr(1)+strdia
   corstr(2)=corstr(2)+strdia
   corstr(3)=corstr(3)+strdia
 end if

 if(nproc_fft > 1) then
   call timab(539,1,tsec)
   if(option==3) then
     call xmpi_sum(corstr,mpi_enreg%comm_fft,ier)
   end if
   if(option==2) then
     call xmpi_sum(grxc,mpi_enreg%comm_fft,ier)
   end if
   call timab(539,2,tsec)
 end if

 DBG_EXIT("COLL")

end subroutine mkcore_paw
!!***
