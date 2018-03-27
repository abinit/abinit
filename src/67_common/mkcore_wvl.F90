!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkcore_wvl
!! NAME
!!  mkcore_wvl
!!
!! FUNCTION
!! Optionally compute (in a WVL representation):
!!  (1) pseudo core electron density throughout unit cell
!!  (2) pseudo-core contribution to forces
!!  (3) pseudo-core contribution to stress tensor
!!
!! COPYRIGHT
!!  Copyright (C) 2016-2018 ABINIT group (MT,TRangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  [mpi_comm_wvl]=MPI communicator (optional)
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=dimension of vxc (XC potential)
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of xccc3d (pseudo core charge)
!!  option: 1 for computing core charge density
!!          2 for computing core charge contribution to forces
!!          3 for computing core charge contribution to stress tensor
!!          4 for computing contribution to frozen-wf part of dynamical matrix
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree)
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and 5 derivatives for each atom type
!!  xred(3,natom)=reduced coordinates for atoms in unit cell
!!  wvl_den=density-potential BigDFT object
!!  wvl_descr=wavelet BigDFT object
!!
!! OUTPUT
!!  === if option==1 ===
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  === if option==2 ===
!!  grxc(3,natom)=core charge contribution to forces
!!  === if option==3 ===
!!  corstr(6)=core charge contribution to stress tensor
!!
!! SIDE EFFECTS
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!   (computed and returned when option=1, needed as input when option=3)
!!
!! NOTES
!! Based on mkcore.F90. Adapted to WVL case.
!!
!! PARENTS
!!      forces,setvtr
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

subroutine mkcore_wvl(atindx1,corstr,grxc,natom,nattyp,nfft,nspden,ntypat,n1xccc,n3xccc,option,&
&                     pawrad,pawtab,rprimd,vxc,xccc1d,xccc3d,xcccrc,xred,wvl_den,wvl_descr,&
&                     mpi_comm_wvl) ! optional argument

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_sort, only : sort_dp
 use m_geometry,   only : xcart2xred, xred2xcart, metric
 use m_paw_numeric, only : paw_splint
 use m_pawrad, only : pawrad_type,pawrad_init,pawrad_free
 use m_pawtab, only : pawtab_type

 use defs_wvltypes
#if defined HAVE_BIGDFT
 use BigDFT_API, only : PSPCODE_PAW,ind_positions
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkcore_wvl'
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_41_xc_lowlevel
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat,n1xccc,n3xccc,option
 integer,intent(in),optional :: mpi_comm_wvl
 type(wvl_denspot_type), intent(inout) :: wvl_den
 type(wvl_internal_type), intent(in) :: wvl_descr
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat)
 real(dp),intent(in) :: rprimd(3,3),xccc1d(n1xccc,6,ntypat),xcccrc(ntypat),xred(3,natom)
 real(dp),intent(in),target :: vxc(nfft,nspden)
 real(dp),intent(out) :: corstr(6),grxc(3,natom)
 real(dp),intent(inout) :: xccc3d(n3xccc)
 type(pawrad_type),intent(in) :: pawrad(:)
 type(pawtab_type),intent(in) :: pawtab(:)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: iat,iatm,iatom,iex,iey,iez,ind,ioffset,ipts,isx,isy,isz,itypat
 integer :: iwarn=0,i1,i2,i3,i3s,j1,j2,j3,jj,jpts,me_wvl,nproc_wvl
 integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,npts,npts12
 integer :: ntot,n1,n1i,n2,n2i,n3,n3i,n3pi
 logical :: perx,pery,perz,gox,goy,goz,USE_PAW
 real(dp) :: aa,arg,bb,cc,cutoff,dd,delta=0,deltam1=0,delta2div6=0,diff
 real(dp) :: hxh,hyh,hzh,range,range2,rangem1,rr,rx,ry,rz,r2,strdia
 real(dp) :: term,ucvol,xx,yy,zz
 character(len=500) :: msg
 type(pawrad_type) :: core_mesh
!arrays
 integer,allocatable :: indx(:),ivec(:)
 real(dp) :: corfra(3,3),corgr(3),gmet(3,3),gprimd(3,3),rmet(3,3),tsec(2),tt(3),xcart(3)
 real(dp),allocatable :: dtcore(:),d2tcore(:),rnorm(:),tcore(:),vecx(:),vecy(:)
 real(dp),pointer :: vxc_eff(:)
#endif

!************************************************************************

#if defined HAVE_BIGDFT

 call timab(12,1,tsec)

!Make sure option is acceptable
 if (option<0.or.option>3) then
   write(msg,'(a,i2,a)') 'Option= ',option,' is not allowed!'
   MSG_BUG(MSG)
 end if
 if(nfft/=n3xccc)then
   write(msg,'(a)') 'nfft and n3xccc should be equal!'
   MSG_BUG(msg)
 end if

!MPI
 nproc_wvl=1;if (present(mpi_comm_wvl)) nproc_wvl=xmpi_comm_size(mpi_comm_wvl)
 me_wvl=0;if (present(mpi_comm_wvl)) me_wvl=xmpi_comm_rank(mpi_comm_wvl)

!Zero out only the appropriate array according to option
 if (option==1) then
   xccc3d(:)=zero
 else if (option==2) then
   grxc(:,:)=zero
 else if (option==3) then
   corfra(:,:)=zero
   corstr(:)=zero
   strdia=zero
 end if

!Nothing to do if no xy plane to handle
 n3pi=wvl_den%denspot%dpbox%n3pi
 if (n3pi==0) return

!Show how calculation runs
 if (me_wvl==0) then
   if (option==1) write(std_out,'(a)',advance='no') ' Compute pseudo core density...'
   if (option==2) write(std_out,'(a)',advance='no') ' Compute forces due to core density...'
   if (option==3) write(std_out,'(a)',advance='no') ' Compute stresses due to core density...'
   if (option==4) write(std_out,'(a)',advance='no') ' Compute dyn. matrix due to core density...'
 end if

!PAW or NCPP ?
 USE_PAW=any(wvl_descr%atoms%npspcode==PSPCODE_PAW)

!Conditions for periodicity in the three directions
 perx=(wvl_descr%atoms%astruct%geocode /= 'F')
 pery=(wvl_descr%atoms%astruct%geocode == 'P')
 perz=(wvl_descr%atoms%astruct%geocode /= 'F')
 call ext_buffers(perx,nbl1,nbr1)
 call ext_buffers(pery,nbl2,nbr2)
 call ext_buffers(perz,nbl3,nbr3)

 if (option>=2) then
!  For spin-polarization, replace vxc by (1/2)*(vxc(up)+vxc(down))
!  For non-collinear magnetism, replace vxc by (1/2)*(vxc^{11}+vxc^{22})
   if (nspden>=2) then
     ABI_ALLOCATE(vxc_eff,(nfft))
     do jj=1,nfft
       vxc_eff(jj)=half*(vxc(jj,1)+vxc(jj,2))
     end do
   else
     vxc_eff => vxc(1:nfft,1)
   end if
 end if

!Some constants
 n1=wvl_descr%Glr%d%n1
 n2=wvl_descr%Glr%d%n2
 n3=wvl_descr%Glr%d%n3
 n1i=wvl_den%denspot%dpbox%ndims(1)
 n2i=wvl_den%denspot%dpbox%ndims(2)
 n3i=wvl_den%denspot%dpbox%ndims(3)
 ntot=n1i*n2i*n3i
 ioffset=n1i*n2i*wvl_den%denspot%dpbox%i3xcsh
 i3s=1+wvl_den%denspot%dpbox%nscatterarr(me_wvl,3)-wvl_den%denspot%dpbox%i3xcsh
 hxh=wvl_den%denspot%dpbox%hgrids(1)
 hyh=wvl_den%denspot%dpbox%hgrids(2)
 hzh=wvl_den%denspot%dpbox%hgrids(3)
 call metric(gmet,gprimd,-1,rmet,rprimd,arg)
 ucvol=real(ntot,dp)*(hxh*hyh*hzh)
 if (.not.USE_PAW) then
   delta=one/(n1xccc-1)
   deltam1=n1xccc-1
   delta2div6=delta**2/6.0_dp
 end if

!Loop over atom types
 iatm=0
 do itypat=1,ntypat

!  Set search range (density cuts off perfectly beyond range)
   range=xcccrc(itypat);if (USE_PAW) range=pawtab(itypat)%rcore
   range2=range**2 ; rangem1=one/range

!  Skip loop if this type has no core charge
   if (abs(range)<1.d-16) cycle

!  PAW: create mesh for core density
   if (USE_PAW) then
     call pawrad_init(core_mesh,mesh_size=pawtab(itypat)%core_mesh_size,&
&     mesh_type=pawrad(itypat)%mesh_type,&
&     rstep=pawrad(itypat)%rstep,lstep=pawrad(itypat)%lstep)
   end if

!  Loop over atoms of the type
   do iat=1,nattyp(itypat)
     iatm=iatm+1;iatom=atindx1(iatm)

     if (option==2) corgr(:)=zero

!    Coordinates of the center
     call xred2xcart(1,rprimd,xcart,xred(:,iatom))
     rx=xcart(1) ; ry=xcart(2) ; rz=xcart(3)

!    Range of points to explore
     cutoff=1.1_dp*range
     isx=floor((rx-cutoff)/hxh)
     isy=floor((ry-cutoff)/hyh)
     isz=floor((rz-cutoff)/hzh)
     iex=ceiling((rx+cutoff)/hxh)
     iey=ceiling((ry+cutoff)/hyh)
     iez=ceiling((rz+cutoff)/hzh)

!    Allocate temporary memory
     !npts12=1+int(((range/hh(1))*(range/hh(2))*pi))
     npts12=(iex-isx+1)*(iey-isy+1)
     ABI_ALLOCATE(rnorm,(npts12))
     ABI_ALLOCATE(vecx,(npts12))
     ABI_ALLOCATE(vecy,(npts12))
     ABI_ALLOCATE(ivec,(npts12))
     ABI_ALLOCATE(indx,(npts12))
     if (option==1.or.option==3) then
       ABI_ALLOCATE(tcore,(npts12))
     end if
     if (option>=2) then
       ABI_ALLOCATE(dtcore,(npts12))
     end if
     if (option==4) then
       ABI_ALLOCATE(d2tcore,(npts12))
     end if

!    Explore range of vectors
     do i3=isz,iez
       zz=real(i3,kind=dp)*hzh-rz
       call ind_positions(perz,i3,n3,j3,goz)
       j3=j3+nbl3+1

!      Select the vectors located around the current atom
!        TR: all of the following  could be done inside or
!        outside the loops (i2,i1,i3).
!        Outside: the memory consumption increases.
!        Inside: the time of calculation increases.
!        Here, I choose to do it here, somewhere in the middle.
       npts=0
       do i2=isy,iey
         yy=real(i2,kind=dp)*hyh-ry
         call ind_positions(pery,i2,n2,j2,goy)
         do i1=isx,iex
           xx=real(i1,kind=dp)*hxh-rx
           call ind_positions(perx,i1,n1,j1,gox)
           r2=xx**2+yy**2+zz**2
           if ((j3>=i3s.and.j3<=i3s+n3pi-1) .and. (gox.and.goy)) then
             ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
!            Only accept contributions inside defined range
             if (r2<range2) then
               npts=npts+1 ; indx(npts)=npts
               ivec(npts)=ioffset+ind
               rnorm(npts)=sqrt(r2)
               vecx(npts)=xx;vecy(npts)=yy
             end if
           end if
         end do
       end do
       if (npts==0) cycle
       if (npts>npts12) then
         msg='npts>npts12!'
         MSG_BUG(msg)
       end if

!      Evaluate core density (and derivatives) on the set of selected points
       if (USE_PAW) then
!        PAW: use splint routine
         call sort_dp(npts,rnorm(1:npts),indx(1:npts),tol16)
         if (option==1.or.option==3) then
!          Evaluate fit of core density
           call paw_splint(core_mesh%mesh_size,core_mesh%rad, &
&           pawtab(itypat)%tcoredens(:,1), &
&           pawtab(itypat)%tcoredens(:,3),&
&           npts,rnorm(1:npts),tcore(1:npts))
         end if
         if (option>=2) then
!          Evaluate fit of 1-der of core density
           call paw_splint(core_mesh%mesh_size,core_mesh%rad, &
&           pawtab(itypat)%tcoredens(:,2), &
&           pawtab(itypat)%tcoredens(:,4),&
&           npts,rnorm(1:npts),dtcore(1:npts))
         end if
         if (option==4) then
!          Evaluate fit of 2nd-der of core density
           call paw_splint(core_mesh%mesh_size,core_mesh%rad, &
&           pawtab(itypat)%tcoredens(:,3), &
&           pawtab(itypat)%tcoredens(:,5),&
&           npts,rnorm(1:npts),d2tcore(1:npts))
         end if
       else
!        Norm-conserving PP:
!          Evaluate spline fit with method from Numerical Recipes
         do ipts=1,npts
           rr=rnorm(ipts)*rangem1
           jj=1+int(rr*(n1xccc-1))
           diff=rr-(jj-1)*delta
           bb = diff*deltam1 ; aa = one-bb
           cc = aa*(aa**2-one)*delta2div6
           dd = bb*(bb**2-one)*delta2div6
           if (option==1.or.option==3) then
             tcore(ipts)=aa*xccc1d(jj,1,itypat)+bb*xccc1d(jj+1,1,itypat) +&
&             cc*xccc1d(jj,3,itypat)+dd*xccc1d(jj+1,3,itypat)
           end if
           if (option>=2) then
             dtcore(ipts)=aa*xccc1d(jj,2,itypat)+bb*xccc1d(jj+1,2,itypat) +&
&             cc*xccc1d(jj,4,itypat)+dd*xccc1d(jj+1,4,itypat)
           end if
           if (option==4) then
             d2tcore(ipts)=aa*xccc1d(jj,3,itypat)+bb*xccc1d(jj+1,3,itypat) +&
&             cc*xccc1d(jj,5,itypat)+dd*xccc1d(jj+1,5,itypat)
           end if
         end do
       end if

!      Now, perform the loop over selected grid points
       do ipts=1,npts
         rr=rnorm(ipts)
         xx=vecx(indx(ipts))
         yy=vecy(indx(ipts))
         jpts=ivec(indx(ipts))

!        === Evaluate charge density
         if (option==1) then
           xccc3d(jpts)=xccc3d(jpts)+tcore(ipts)

!        === Accumulate contributions to forces
         else if (option==2) then
           if (rr>tol10) then
             term=vxc_eff(jpts)*dtcore(ipts)/rr
             corgr(1)=corgr(1)+xx*term
             corgr(2)=corgr(2)+yy*term
             corgr(3)=corgr(3)+yy*term
           end if

!        === Accumulate contributions to stress tensor (in red. coordinates)
         else if (option==3) then
           if (rr>tol10) then
             term=vxc_eff(jpts)*dtcore(ipts)*rangem1/rr/real(ntot,dp)
!            Write out the 6 symmetric components
             corfra(1,1)=corfra(1,1)+term*xx*xx
             corfra(2,2)=corfra(2,2)+term*yy*yy
             corfra(3,3)=corfra(3,3)+term*zz*zz
             corfra(3,2)=corfra(3,2)+term*zz*yy
             corfra(3,1)=corfra(3,1)+term*zz*xx
             corfra(2,1)=corfra(2,1)+term*yy*xx
!            (the above still needs to be transformed to cartesian coords)
           end if
!          Also compute a diagonal term
           strdia=strdia+vxc_eff(jpts)*tcore(ipts)

         end if ! Choice of option

       end do ! ipts (i1,i2)

     end do ! i3

!    Release temporary memory
     ABI_DEALLOCATE(rnorm)
     ABI_DEALLOCATE(vecx)
     ABI_DEALLOCATE(vecy)
     ABI_DEALLOCATE(ivec)
     ABI_DEALLOCATE(indx)
     if (allocated(tcore)) then
       ABI_DEALLOCATE(tcore)
     end if
     if (allocated(dtcore)) then
       ABI_DEALLOCATE(dtcore)
     end if
     if (allocated(d2tcore)) then
       ABI_DEALLOCATE(d2tcore)
     end if

     if (option==2) then
       arg=-(ucvol/real(ntot,dp))
       !arg=-(ucvol/real(ntot,dp))/range  !????
       grxc(:,iatom)=corgr(:)*arg
     end if

!  End loop on atoms
   end do

   if (USE_PAW) then
     call pawrad_free(core_mesh)
   end if

!End loop over atom types
 end do

 if(option>=2.and.nspden>=2)  then
   ABI_DEALLOCATE(vxc_eff)
 end if

!Density: make it positive
 if (option==1) then
   call mkdenpos(iwarn,n3xccc,nspden,0,xccc3d,tol20)
 end if

!Forces: translate into reduced coordinates
 if (option==2) then
   do iatom=1,natom
     tt(1:3)=grxc(1:3,iatom)
     grxc(:,iatom)= rprimd(1,:)*tt(1)+rprimd(2,:)*tt(2)+rprimd(3,:)*tt(3)
    !grxc(:,iatom)=rmet(:,1)*tt(1)+rmet(:,2)*tt(2)+rmet(:,3)*tt(3)
   end do
 end if

!Stress tensor: symmetrize, translate into cartesian coord., add diagonal part
 if (option==3) then
   corstr(1)=corfra(1,1) ; corstr(2)=corfra(2,2)
   corstr(3)=corfra(3,3) ; corstr(4)=corfra(3,2)
   corstr(5)=corfra(3,1) ; corstr(6)=corfra(2,1)
   call strconv(corstr,rprimd,corstr)
   corstr(1)=corstr(1)+strdia/real(ntot,dp)
   corstr(2)=corstr(2)+strdia/real(ntot,dp)
   corstr(3)=corstr(3)+strdia/real(ntot,dp)
 end if

!If needed, sum over MPI processes
 if(nproc_wvl>1) then
   call timab(539,1,tsec)
   if (option==2) then
     call xmpi_sum(grxc,mpi_comm_wvl,iex)
   end if
   if (option==3) then
     call xmpi_sum(corstr,mpi_comm_wvl,iex)
   end if
   call timab(539,2,tsec)
 end if

 if (me_wvl==0) write(std_out,'(a)') 'done.'

 call timab(12,1,tsec)

#else
 BIGDFT_NOTENABLED_ERROR()
 ABI_UNUSED(xcccrc)
 if (.false.) write(std_out,*) natom,nfft,nspden,ntypat,n1xccc,n3xccc,option,mpi_comm_wvl,&
& wvl_den%symObj,wvl_descr%h(1),atindx1(1),nattyp(1),rprimd(1,1),vxc(1,1),&
& xred(1,1),xccc1d(1,1,1),corstr(1),grxc(1,1),xccc3d(1),pawrad(1)%mesh_size,pawtab(1)%lmn_size
#endif

end subroutine mkcore_wvl
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/mkcore_wvl_old
!! NAME
!!  mkcore_wvl_old
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
!!  Copyright (C) 2011-2018 ABINIT group (T. Rangel)
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

subroutine mkcore_wvl_old(atindx1,corstr,dyfrx2,geocode,grxc,h,natom,&
& nattyp,nfft,nscatterarr,nspden,ntypat,n1,n1i,n2,n2i,n3,n3pi,&
& n3xccc,option,pawrad,pawtab,psppar,rprimd,ucvol,&
& vxc,xccc3d,xred,mpi_comm_wvl)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_splines

 use m_geometry, only : xred2xcart, xcart2xred
 use m_pawrad,  only : pawrad_type, pawrad_init, pawrad_free
 use m_pawtab,  only : pawtab_type
 use m_xmpi,    only : xmpi_comm_size,xmpi_sum

#if defined HAVE_BIGDFT
  use BigDFT_API, only: ext_buffers,ind_positions
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkcore_wvl_old'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_41_xc_lowlevel
 use interfaces_67_common, except_this_one => mkcore_wvl_old
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
 integer,allocatable :: ifftsph_tmp(:)
 real(dp),allocatable:: rr(:),rred(:,:)
#endif

! *************************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT

 if(nspden >1) then
   write(message, '(a)')'mkcore_wvl_old: this is not yet generalized to npsden>1'
   MSG_ERROR(message)
 end if
 if(option>4 .or. option<1 )then
   write(message,'(a,a,a,a,a,a,i6)') ch10,&
&   ' mkcore_wvl_old: BUG -',ch10,&
&   '  The argument option should be between 1 and 4,',ch10,&
&   '  however, option=',option
   MSG_BUG(message)
 end if
 if(nfft .ne. n3xccc)then
   write(message,'(a,a,a,a,a,a,2i6)') ch10,&
&   ' mkcore_wvl_old: BUG -',ch10,&
&   '  nfft and n3xccc should be equal,',ch10,&
&   '  however, nfft and n3xccc=',nfft,n3xccc
   MSG_BUG(message)
 end if

!MPI
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
&   ' mkcore_wvl_old: BUG -',ch10,&
&   '  Can''t be here ! (bad option).'
   MSG_BUG(message)
 end if

 write(message,'(a,a)') ch10,&
& ' mkcore_wvl_old: Compute core density'
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
     iatom_tot=iatom

     if(option==2) then
       grxc1=zero
       grxc2=zero
       grxc3=zero
     end if

!    Define a "box" around each atom
     rx=xcart(1,iatom_tot)
     ry=xcart(2,iatom_tot)
     rz=xcart(3,iatom_tot)

     isx=floor((rx-cutoff)/hh(1))
     isy=floor((ry-cutoff)/hh(2))
     isz=floor((rz-cutoff)/hh(3))

     iex=ceiling((rx+cutoff)/hh(1))
     iey=ceiling((ry+cutoff)/hh(2))
     iez=ceiling((rz+cutoff)/hh(3))

     do i3=isz,iez
       zz=real(i3,kind=8)*hh(3)-rz
       call ind_positions(perz,i3,n3,j3,goz)
       j3=j3+nbl3+1

!      Initialize counters
       nfgd=0
       nfgd_r0=0

       do i2=isy,iey
         yy=real(i2,kind=8)*hh(2)-ry
         call ind_positions(pery,i2,n2,j2,goy)

         do i1=isx,iex
           xx=real(i1,kind=8)*hh(1)-rx
           call ind_positions(perx,i1,n1,j1,gox)
           rr2=xx**2+yy**2+zz**2
           if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then

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
               else if (option==4) then
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
&       rmet,rr,strdia,vxc,xccc3d,rred=rred)

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
   ABI_DEALLOCATE(rr)
   ABI_DEALLOCATE(rred)

 end do !itypat

 if (option==2) then
!  Apply rmet as needed to get reduced coordinate gradients
   do iatom=1,natom
     t1=grxc(1,iatom)
     t2=grxc(2,iatom)
     t3=grxc(3,iatom)
     grxc(:,iatom)=rmet(:,1)*t1+rmet(:,2)*t2+rmet(:,3)*t3
   end do

 elseif (option==3) then

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

end subroutine mkcore_wvl_old
!!***
