!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_initro
!! NAME
!!  wvl_initro
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2016 ABINIT group (TRangel)
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
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      ext_buffers,ind_positions,metric,mkdenpos,pawrad_free,pawrad_init
!!      sort_dp,splint,wrtout,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

subroutine wvl_initro(&
& atindx1,geocode,h,me,&
& natom,nattyp,nfft,nspden,ntypat,&
& n1,n1i,n2,n2i,n3,&
& pawrad,pawtab,psppar,&
& rhor,rprimd,spinat,wvl_den,xc_denpos,xred,zion)

 use defs_basis
 use m_profiling_abi
 use m_splines
 use m_errors
 use defs_wvltypes
 use m_pawrad,  only : pawrad_type, pawrad_init, pawrad_free
 use m_pawtab, only : pawtab_type

#if defined HAVE_BIGDFT
  use BigDFT_API, only : ELECTRONIC_DENSITY, ext_buffers, ind_positions
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_initro'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_41_geometry
 use interfaces_41_xc_lowlevel
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: me,natom,ntypat,nfft,nspden
 integer,intent(in)::n1,n2,n1i,n2i,n3
 real(dp),intent(in) :: h(3)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)
 real(dp),intent(in) :: spinat(3,natom),zion(ntypat)
 real(dp),intent(inout) :: rhor(nfft,nspden)
 real(dp),intent(in) :: xc_denpos
 character(1),intent(in)::geocode
 type(wvl_denspot_type), intent(inout) :: wvl_den
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat)
 real(dp),intent(in) :: psppar(0:4,0:6,ntypat),rprimd(3,3)
 real(dp),intent(inout)::xred(3,natom)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer  :: ia1,ia2
 integer  :: iat,iatm,iatom,iatom_tot,iex,iey,iez,ii,ind
 integer  :: ifft,ispden
 integer  :: isx,isy,isz,itypat,i1,i2,i3,iwarn,i3s
 integer  :: j1,j2,j3,msz
 integer  :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3
 integer  :: ncmax,nfgd,nspden_updn,n3pi,shift
 real(dp) :: cutoff,fact,fact0
 real(dp) :: rloc,rr2,rx,ry,rz
 real(dp) :: rshp,r2shp
 real(dp) :: ucvol,xx,yy,zz
 type(pawrad_type)::vale_mesh
!arrays
 logical :: perx,pery,perz,gox,goy,goz
 real(dp) :: hh(3) !fine grid spacing for wavelets
 real(dp) :: gmet(3,3),gprimd(3,3),rcart(3),rmet(3,3),xcart(3,natom)
 character(len=500) :: message                   ! to be uncommented, if needed
!allocatable arrays
 integer,allocatable :: ifftsph_tmp(:),iindex(:)
 real(dp),allocatable:: raux(:),raux2(:)
 real(dp),allocatable:: rr(:)!,rred(:,:)
#endif
 
! ************************************************************************* 

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT

!PENDING: PARALLELIZATION OVER ATOMS

 write(message,'(a,a)') ch10,&
& ' wvl_initro: Initialize valence density from atomic data by splines'
 call wrtout(std_out,message,'COLL')

!initialize
 rhor(:,:)=zero

 if(nspden==4)then
   write(std_out,*)' initro : might work yet for nspden=4 (not checked)'
   write(std_out,*)'spinat',spinat(1:3,1:natom)
!  stop
 end if

!Check whether the values of spinat are acceptable
 if(nspden==2)then
   do itypat=1,ntypat
     do iat=1,nattyp(itypat)
       iatm=iatm+1;iatom=atindx1(iatm)
       iatom_tot=iatom; !if (mpi_enreg%nproc_atom>1) iatom_tot=mpi_enreg%atom_indx(iatom)

       if( sqrt(spinat(1,iatom)**2+spinat(2,iatom)**2+spinat(3,iatom)**2) &
&       > abs(zion(itypat))*(1.0_dp + epsilon(0.0_dp)) ) then
         write(message, '(a,a,a,a,i4,a,a,3es11.4,a,a,a,es11.4)' ) ch10,&
&         ' initro : WARNING - ',ch10,&
&         '  For atom number ',iatom,ch10,&
&         '  input spinat=',spinat(:,iatom),'  is larger, in magnitude,',ch10,&
&         '  than zion(ia)=',zion(itypat)
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,message,'COLL')
       end if
     end do
     ia1=ia2+1
   end do
 end if

!Fine grid
 hh(:)=0.5d0*h(:)

!mpi:
!Obtain n3pi, BigDFT quantity:
 n3pi=wvl_den%denspot%dpbox%n3pi
 i3s=wvl_den%denspot%dpbox%nscatterarr(me,3)+1-wvl_den%denspot%dpbox%nscatterarr(me,4)
 shift=n1i*n2i*wvl_den%denspot%dpbox%nscatterarr(me,4)
 
!Compute xcart from xred
 call xred2xcart(natom,rprimd,xcart,xred)

!Compute metric tensors and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

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

!  Create mesh_core object
!  since tnvale_mesh_size can be bigger than pawrad%mesh_size,
   msz=pawtab(itypat)%tnvale_mesh_size
   call pawrad_init(vale_mesh,mesh_size=msz,mesh_type=pawrad(itypat)%mesh_type,&
&   rstep=pawrad(itypat)%rstep,lstep=pawrad(itypat)%lstep)
!  
!  Set radius size:
   rshp=vale_mesh%rmax
   r2shp=1.0000001_dp*rshp**2

!  allocate arrays
   if (n3pi > 0) then
!    sphere: cycle i1,i2,i3
!    ncmax=1+int(1.1_dp*nfft*four_pi/(three*ucvol)*rshp**3)
!    ncmax=1+int(1.1_dp*nfft*four_pi/(three*ucvol)*rshp**3)
!    1+int(1.1* factors are included just for cautioness
!    circle: cycle only i1 and i2
!    ncmax=1+int(1.1d0*((rshp/hh(1))*(rshp/hh(2))*pi))
!    line:
     ncmax=1+int(1.1_dp*rshp/hh(1)*2.d0)
   else
     ncmax=1
   end if
!  
   ABI_ALLOCATE(ifftsph_tmp,(ncmax))
   ABI_ALLOCATE(iindex,(ncmax))
   ABI_ALLOCATE(rr,(ncmax))
   ABI_ALLOCATE(raux,(ncmax))
   if(nspden==2) then
     ABI_ALLOCATE(raux2,(ncmax))
   end if

!  Big loop on atoms  
   do iat=1,nattyp(itypat)
     iatm=iatm+1;iatom=atindx1(iatm)
     iatom_tot=iatom; !if (mpi_enreg%nproc_atom>1) iatom_tot=mpi_enreg%atom_indx(iatom)

!    Spin
     if(nspden==2) then 
       fact0=half/zion(itypat)
       fact=fact0*(zion(itypat)+spinat(3,iatom))
     end if

!    
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
!      
       do i2=isy,iey
         yy=real(i2,kind=8)*hh(2)-ry
         call ind_positions(pery,i2,n2,j2,goy)
!        
!        Initialize counters
         nfgd=0
!        nfgd_r0=0
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
                 rcart=[xx,yy,zz]
                 rr(nfgd)=(rr2)**0.5
                 ifftsph_tmp(nfgd)=shift+ind
!                DEBUG     
!                write(itmp,'(i10,3(f13.7,x))')ind,xx+rx,yy+ry,zz+rz
!                write(itmp,'(6(f13.7,x))')rcart,rred(:,nfgd)
!                ENDDEBUG
!                else
!                !              We save r=0 vectors 
!                ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s)*n1i*n2i
!                !              We reuse the same variable "ifftshp_tmp", 
!                !              but we start from the higher index
!                nfgd_r0=nfgd_r0+1
!                ifftsph_tmp(ncmax-nfgd_r0+1)=shift+ind
               end if !rr2>tol5
             end if !rr2<r2shp
           end if !j3..
         end do !i1

!        All of the following  could be done inside or outside the loops (i2,i1,i3)
!        Outside the loops: the memory consuption increases.
!        Inside the inner loop: the time of calculation increases.

         if(nfgd==0)      cycle

!        Evaluate spline fit of 1st der of core charge density
!        from tcoredens(:,2) and tcoredens(:,4)
         do ii=1,nfgd
           iindex(ii)=ii
         end do
!        write(600,'(i4,x,9999f14.7)')nfgd, rr(1:nfgd)
         call sort_dp(nfgd,rr(1:nfgd),iindex(1:nfgd),tol16)
         call splint(msz,vale_mesh%rad,&
&         pawtab(itypat)%tvalespl(:,1),pawtab(itypat)%tvalespl(:,2),&
&         nfgd,rr(1:nfgd),raux(1:nfgd))


!        Accumulate contributions to valence density on the entire cell
         rhor(ifftsph_tmp(1:nfgd),1)=rhor(ifftsph_tmp(1:nfgd),1)+raux(iindex(1:nfgd))

         if(nspden==2) then
           raux2(1:nfgd)=raux(iindex(1:nfgd))*fact
           rhor(ifftsph_tmp(1:nfgd),2)=rhor(ifftsph_tmp(1:nfgd),1)+raux2(1:nfgd)
         end if
!        DEBUG
!        do ii=1,msz
!        write(itmp,'(2(f15.7,1x))')vale_mesh%rad(ii),pawtab(itypat)%tvalespl(ii,1)
!        end do
!        do ii=1,nfgd
!        write(itmp,'(2i10)')ii,iindex(ii)
!        write(itmp,'(2(f15.7,1x))')rr(iindex(ii)),rhor(ifftsph_tmp(ii),1)!,raux(iindex(ii))
!        end do
!        END DEBUG

       end do !i2
     end do !i1
   end do !iat

!  Deallocate
   call pawrad_free(vale_mesh)
   ABI_DEALLOCATE(ifftsph_tmp)
   ABI_DEALLOCATE(iindex)
   ABI_DEALLOCATE(rr)
   ABI_DEALLOCATE(raux)
   if(nspden==2) then
     ABI_DEALLOCATE(raux2)
   end if

 end do !itypat

!nspden_updn: 1 for non-polarized, 2 for polarized
 nspden_updn=min(nspden,2)

!Make the density positive everywhere
 call mkdenpos(iwarn,nfft,nspden_updn,1,rhor(:,1:nspden_updn),xc_denpos)

!There seems to be a bug in the intel11 compiler
!rhor = reshape(wvl_den%denspot%rhov, shape(rhor))
 do ispden=1,nspden
   do ifft=1,nfft
     ii=ifft+nfft*(ispden-1)
!    rhor(ifft,ispden)=wvl_den%denspot%rhov(ii)
     wvl_den%denspot%rhov(ii)=rhor(ifft,ispden)
   end do
 end do
 wvl_den%denspot%rhov_is = ELECTRONIC_DENSITY
 write(message, '(a,a,a,a)' ) ch10, ' wvl_initro : but why are you copying me :..o('
 call wrtout(std_out,message,'COLL')


#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) me,natom,ntypat,nfft,nspden,n1,n2,n1i,n2i,n3,h(1),&
& pawrad(1)%mesh_size,pawtab(1)%mesh_size,spinat(1,1),zion(1),rhor(1,1),xc_denpos,&
& geocode,wvl_den%symObj,atindx1(1),nattyp(1),psppar(1,1,1),rprimd(1,1),xred(1,1)
#endif

 DBG_EXIT("COLL")

end subroutine wvl_initro
!!***
