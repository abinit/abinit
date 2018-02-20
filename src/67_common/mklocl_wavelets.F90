!{\src2tex{textfont=tt}}
!!****f* ABINIT/mklocl_wavelets
!!
!! NAME
!! mklocl_wavelets
!!
!! FUNCTION
!! Compute the ionic local potential when the pseudo-potentials are GTH, using
!! the special decomposition of these pseudo. The resulting potential is computed with
!! free boundary conditions. It gives the same result than mklocl_realspace for the
!! GTH pseudo only with a different way to compute the potential.
!!
!! Optionally compute :
!!  option=1 : local ionic potential throughout unit cell
!!  option=2 : contribution of local ionic potential to E gradient wrt xred
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DC,TRangel,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  efield (3)=external electric field
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms
!!  nfft=size of vpsp (local potential)
!!  nspden=number of spin-density components
!!  option=type of calculation (potential, forces, ...)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  xcart(3,natom)=cartesian atomic coordinates.
!!  wvl_den=density-potential BigDFT object
!!  wvl_descr=wavelet BigDFT object
!!
!! OUTPUT
!!  (if option==1) vpsp(nfft)=the potential resulting from the ionic
!!                 density of charge.
!!  (if option==2) grtn(3,natom)=grads of Etot wrt tn. These gradients are in
!!                 reduced coordinates. Multiply them by rprimd to get
!!                 gradients in cartesian coordinates.
!!
!! PARENTS
!!      mklocl,wvl_wfsinp_scratch
!!
!! CHILDREN
!!      calcdvloc_wvl,derf_ab,paw_splint_der
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mklocl_wavelets(efield, grtn, mpi_enreg, natom, nfft, &
     & nspden, option, rprimd, vpsp, wvl_den, wvl_descr, xcart)

 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_abi2big, only : wvl_rhov_abi2big
 use m_profiling_abi
 use m_xmpi
 use m_errors
#if defined HAVE_BIGDFT
 use BigDFT_API, only : ELECTRONIC_DENSITY,createIonicPotential,local_forces
 use poisson_solver, only : H_potential
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mklocl_wavelets'
 use interfaces_14_hidewrite
 use interfaces_67_common, except_this_one => mklocl_wavelets
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option, natom, nfft, nspden
 type(MPI_type),intent(in) :: mpi_enreg
 type(wvl_denspot_type), intent(inout) :: wvl_den
 type(wvl_internal_type), intent(in) :: wvl_descr
!arrays
 real(dp),intent(in) :: rprimd(3,3),efield(3)
 real(dp),intent(inout) :: grtn(3,natom)
 real(dp), intent(inout) :: vpsp(nfft)
 real(dp),intent(inout) :: xcart(3,natom)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: i,i1,i2,i3,ia,ierr,igeo,me,nproc,shift,comm
 real(dp) :: energ
 character(len=500) :: message
!arrays
 real(dp) :: epot(3)
 real(dp) :: elecfield(3)=(/zero,zero,zero/) ! Not used here
 real(dp),allocatable :: gxyz(:,:),vhartr(:),rhov(:,:)
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 elecfield=zero !not used here

!Manage parallelism
 comm=mpi_enreg%comm_wvl
 nproc=xmpi_comm_size(comm)
 me=xmpi_comm_rank(comm)

!----------------------------------------------------------------------
! ----- Option 1: compute local ionic potential                   -----
!----------------------------------------------------------------------
 if (option == 1) then

!  We get the kernel for the Poisson solver (used to go from the ionic
!  charge to the potential).If the kernel is uncomputed, it does it now.
!  call psolver_kernel(wvl_den%denspot%dpbox%hgrids, 2, icoulomb, me, kernel, &
!&                     comm, wvl_den%denspot%dpbox%ndims, nproc, nscforder)
!  if (.not.associated(kernel%co%kernel)) then
!    call psolver_kernel(wvl_den%denspot%dpbox%hgrids, 1, icoulomb, me, kernel, &
!&                       comm, wvl_den%denspot%dpbox%ndims, nproc, nscforder)
!  end if

   message=ch10//' mklocl_wavelets: Create local potential from ions.'
   call wrtout(std_out,message,'COLL')

   shift = 1 + wvl_den%denspot%dpbox%ndims(1) * wvl_den%denspot%dpbox%ndims(2) &
&   * wvl_den%denspot%dpbox%i3xcsh

!  Call the BigDFT routine
   call createIonicPotential(wvl_descr%atoms%astruct%geocode, me, nproc, (me == 0), wvl_descr%atoms, &
&   xcart, wvl_den%denspot%dpbox%hgrids(1), wvl_den%denspot%dpbox%hgrids(2), &
&   wvl_den%denspot%dpbox%hgrids(3), &
&   elecfield, wvl_descr%Glr%d%n1, wvl_descr%Glr%d%n2, wvl_descr%Glr%d%n3, &
&   wvl_den%denspot%dpbox%n3pi, wvl_den%denspot%dpbox%i3s + wvl_den%denspot%dpbox%i3xcsh, &
&   wvl_den%denspot%dpbox%ndims(1), wvl_den%denspot%dpbox%ndims(2), &
&   wvl_den%denspot%dpbox%ndims(3), wvl_den%denspot%pkernel, vpsp(shift), 0.d0,wvl_descr%rholoc)

!  Copy vpsp into bigdft object:
   call wvl_rhov_abi2big(1,vpsp,wvl_den%denspot%v_ext,shift=(shift-1))

!  Eventually add the electric field
   if (maxval(efield) > tol12) then
     message=ch10//'mklocl_wavelets: Add the electric field.'
     call wrtout(std_out,message,'COLL')
!    We add here the electric field since in BigDFT, the field must be on x...
     epot(:) = real(0.5, dp) * efield(:) * wvl_den%denspot%dpbox%hgrids(:)
     do i3 = 1, wvl_den%denspot%dpbox%n3pi, 1
       ia = (i3 - 1) * wvl_den%denspot%dpbox%ndims(1) * wvl_den%denspot%dpbox%ndims(2)
       do i2 = -14, 2 * wvl_descr%Glr%d%n2 + 16, 1
         i = ia + (i2 + 14) * wvl_den%denspot%dpbox%ndims(1)
         do i1 = -14, 2 * wvl_descr%Glr%d%n1 + 16, 1
           i = i + 1
           vpsp(shift + i) = vpsp(shift + i) + &
&           epot(1) * real(i1 - wvl_descr%Glr%d%n1, dp) + &
&           epot(2) * real(i2 - wvl_descr%Glr%d%n2, dp) + &
&           epot(3) * real(i3 - wvl_descr%Glr%d%n3, dp)
         end do
       end do
     end do
   end if

!----------------------------------------------------------------------
! ----- Option 2: compute forces induced by local ionic potential -----
!----------------------------------------------------------------------
 else if (option == 2) then

   message=ch10//' mklocl_wavelets: compute local forces.'
   call wrtout(std_out,message,'COLL')

   if (wvl_den%denspot%rhov_is/=ELECTRONIC_DENSITY) then
     message='denspot bigdft datstructure should contain rhor!'
     MSG_BUG(message)
   end if

!  Extract density rhor from bigDFT datastructure
   ABI_ALLOCATE(rhov,(nfft, nspden))
   ABI_ALLOCATE(vhartr,(nfft))
   shift = wvl_den%denspot%dpbox%ndims(1) * wvl_den%denspot%dpbox%ndims(2) &
&   * wvl_den%denspot%dpbox%i3xcsh
   do i = 1, nfft
     rhov(i, 1) = wvl_den%denspot%rhov(i + shift)
     vhartr(i)  = wvl_den%denspot%rhov(i + shift)
   end do
   if (nspden == 2) then
     shift = shift + wvl_den%denspot%dpbox%ndims(1) * wvl_den%denspot%dpbox%ndims(2) &
&     * wvl_den%denspot%dpbox%n3d
     do i = 1, nfft
       rhov(i, 2) =             wvl_den%denspot%rhov(i + shift)
       vhartr(i)  = vhartr(i) + wvl_den%denspot%rhov(i + shift)
     end do
   end if

!  Compute Hartree potential from rhor
   call H_potential('D',wvl_den%denspot%pkernel,vhartr,vhartr,energ,zero,.false.)

!  Allocate temporary array for forces
   ABI_ALLOCATE(gxyz,(3, natom))

!  Calculate local part of the forces grtn (modified BigDFT routine)
   call local_forces_wvl(me,natom,xcart,&
&   wvl_den%denspot%dpbox%hgrids(1),&
&   wvl_den%denspot%dpbox%hgrids(2),&
&   wvl_den%denspot%dpbox%hgrids(3),&
&   wvl_descr%Glr%d%n1,wvl_descr%Glr%d%n2,wvl_descr%Glr%d%n3,&
&   wvl_den%denspot%dpbox%n3p,&
&   wvl_den%denspot%dpbox%i3s+wvl_den%denspot%dpbox%i3xcsh,&
&   wvl_den%denspot%dpbox%ndims(1),wvl_den%denspot%dpbox%ndims(2),&
&   rhov,vhartr,gxyz,wvl_descr)
!    call local_forces(me, wvl_descr%atoms, xcart, &
! &   wvl_den%denspot%dpbox%hgrids(1), wvl_den%denspot%dpbox%hgrids(2), &
! &   wvl_den%denspot%dpbox%hgrids(3), &
! &   wvl_descr%Glr%d%n1, wvl_descr%Glr%d%n2, wvl_descr%Glr%d%n3, &
! &   wvl_den%denspot%dpbox%n3p, wvl_den%denspot%dpbox%i3s + wvl_den%denspot%dpbox%i3xcsh, &
! &   wvl_den%denspot%dpbox%ndims(1), wvl_den%denspot%dpbox%ndims(2), &
! &   rhov, vhartr,gxyz,locstrten,charge)

!    Pending: floc,locstrten and charge are not used here.
!    Pending: put mpi_enreg%nscatterarr... in object denspot, initialize object, etc.

   if (nproc > 1) then
     call xmpi_sum(gxyz, comm, ierr)
   end if

!  Forces should be in reduced coordinates.
   do ia = 1, natom, 1
     do igeo = 1, 3, 1
       grtn(igeo, ia) = - rprimd(1, igeo) * gxyz(1, ia) &
&       - rprimd(2, igeo) * gxyz(2, ia) &
&       - rprimd(3, igeo) * gxyz(3, ia)
     end do
   end do

!  Deallocate local variables
   ABI_DEALLOCATE(vhartr)
   ABI_DEALLOCATE(rhov)
   ABI_DEALLOCATE(gxyz)

!----------------------------------------------------------------------

 else ! option switch
   message = 'Internal error, option should be 1 or 2!'
   MSG_ERROR(message)
 end if
 
#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) option,natom,nfft,nspden,mpi_enreg%me,&
& wvl_den%symObj,wvl_descr%h(1),rprimd(1,1),efield(1),grtn(1,1),vpsp(1),xcart(1,1)
#endif

end subroutine mklocl_wavelets
!!***

!----------------------------------------------------------------------

!!****f* mklocl_wavelets/local_forces_wvl
!! NAME
!!  local_forces_wvl
!!
!! FUNCTION
!!
!! INPUTS
!!  hxh,hyh,hzh=wavelet grid spacings
!!  iproc=current MPI process number
!!  n1,n2,n3=number of wavelet points in each direction
!!  n1i,n2i=size of intermediate xy wvl grid
!!  n3pi=number of xy wvl planes handled by current MPI process
!!  i3s=starting index of local potential for current MPI process
!!  natom=number of atoms
!!  pot(*)=Hartree ionic potential
!!  rho(*)=electronic density
!!  rxyz(3,natom)=cartesian coordinates of atoms
!!  wvl=wavelet BigDFT object
!!
!! OUTPUT
!!  floc(3,natom)=local ionic potential contribution to forces
!!
!! PARENTS
!!      mklocl_wavelets
!!
!! CHILDREN
!!      calcdvloc_wvl,derf_ab,paw_splint_der
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine local_forces_wvl(iproc,natom,rxyz,hxh,hyh,hzh,n1,n2,n3,n3pi,i3s,n1i,n2i,&
&                           rho,pot,floc,wvl)

 use defs_basis
 use defs_wvltypes
 use m_errors
 use m_paw_numeric, only : paw_splint_der
#if defined HAVE_BIGDFT
 use BigDFT_API, only : PSPCODE_PAW,ind_positions
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'local_forces_wvl'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: i3s,iproc,n1,n1i,n2,n2i,n3,n3pi,natom
 real(dp),intent(in) :: hxh,hyh,hzh
 type(wvl_internal_type),intent(in) :: wvl
!arrays
 real(dp),intent(in) :: rxyz(3,natom)
 real(dp),dimension(*),intent(in) :: rho,pot
 real(dp),intent(out) :: floc(3,natom)

!Local variables -------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: i1,i2,i3,iat,iex,iey,iez,iloc,ind,isx,isy,isz,ityp
 integer :: j1,j2,j3,msz=0
 integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,nloc
 logical :: perx,pery,perz,gox,goy,goz,USE_PAW
 real(dp) :: arg,charge,cutoff,dvhgh,forceloc,fxerf,fyerf,fzerf,fxgau,fygau,fzgau
 real(dp) :: forceleaked,prefactor,r2,rhoel,rloc,rloc2,rx,ry,rz,rzero,tt,vel,x,xp,y,z
!arrays
 real(dp) :: cprime(4),dvpawdr(1),rr(1)
 real(dp),pointer :: psppar(:,:),rad(:),vloc(:),d2vloc(:)
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 if (iproc==0) write(std_out,'(a)',advance='no') ' Calculate local forces...'

!PAW or NCPP ?
 USE_PAW=any(wvl%atoms%npspcode==PSPCODE_PAW)

!Conditions for periodicity in the three directions
 perx=(wvl%atoms%astruct%geocode /= 'F')
 pery=(wvl%atoms%astruct%geocode == 'P')
 perz=(wvl%atoms%astruct%geocode /= 'F')
 call ext_buffers(perx,nbl1,nbr1)
 call ext_buffers(pery,nbl2,nbr2)
 call ext_buffers(perz,nbl3,nbr3)

 forceleaked=zero

 do iat=1,natom
   ityp=wvl%atoms%astruct%iatype(iat)

!  Coordinates of the center
   rx=rxyz(1,iat)
   ry=rxyz(2,iat)
   rz=rxyz(3,iat)

!  Initialization of the forces
!  ion-electron term, error function part
   fxerf=zero
   fyerf=zero
   fzerf=zero
!  ion-electron term, gaussian part
   fxgau=zero
   fygau=zero
   fzgau=zero

!  Building array of coefficients of the derivative of the gaussian part
   psppar => wvl%atoms%psppar(:,:,ityp)
   cprime(1)=2._dp*wvl%atoms%psppar(0,2,ityp)-wvl%atoms%psppar(0,1,ityp)
   cprime(2)=4._dp*wvl%atoms%psppar(0,3,ityp)-wvl%atoms%psppar(0,2,ityp)
   cprime(3)=6._dp*wvl%atoms%psppar(0,4,ityp)-wvl%atoms%psppar(0,3,ityp)
   cprime(4)=-wvl%atoms%psppar(0,4,ityp)

!  Determine number of local terms (HGH pot)
   nloc=0
   do iloc=1,4
     if (wvl%atoms%psppar(0,iloc,ityp).ne.zero) nloc=iloc
   end do

!  Some constants depending on the atom type
   rloc=wvl%atoms%psppar(0,0,ityp) ; rloc2=rloc**2
   charge=real(wvl%atoms%nelpsp(ityp),kind=dp)
   prefactor=charge/(2._dp*pi*sqrt(2._dp*pi)*rloc**5)

!  PAW specifics
   if (USE_PAW) then
     msz=wvl%rholoc%msz(ityp)
     rad    => wvl%rholoc%rad(1:msz,ityp)
     vloc   => wvl%rholoc%d(1:msz,3,ityp)
     d2vloc => wvl%rholoc%d(1:msz,4,ityp)
     rzero=rad(1);if (rzero<=1.d-10) rzero=rad(2)
   end if

!  Maximum extension of the gaussian
   cutoff=10._dp*rloc
   isx=floor((rx-cutoff)/hxh)
   isy=floor((ry-cutoff)/hyh)
   isz=floor((rz-cutoff)/hzh)
   iex=ceiling((rx+cutoff)/hxh)
   iey=ceiling((ry+cutoff)/hyh)
   iez=ceiling((rz+cutoff)/hzh)

!  Calculate the forces near the atom due to the gaussian
!  and error function parts of the potential
   if (n3pi>0) then
     do i3=isz,iez
       z=real(i3,kind=dp)*hzh-rz
       call ind_positions(perz,i3,n3,j3,goz)
       j3=j3+nbl3+1
       do i2=isy,iey
         y=real(i2,kind=dp)*hyh-ry
         call ind_positions(pery,i2,n2,j2,goy)
         do i1=isx,iex
           x=real(i1,kind=dp)*hxh-rx
           call ind_positions(perx,i1,n1,j1,gox)

           r2=x**2+y**2+z**2
           xp=exp(-0.5_dp*r2/rloc2)

           if ((j3>=i3s.and.j3<=i3s+n3pi-1) .and. (gox.and.goy)) then
             ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i

!            Short range part
             rhoel=rho(ind)
!            HGH: V_S^prime=gaussian
             if (.not.USE_PAW) then
               if (nloc/=0) then
                 arg=r2/rloc2
                 tt=cprime(nloc)
                 do iloc=nloc-1,1,-1
                   tt=arg*tt+cprime(iloc)
                 end do
                 forceloc=xp*tt*rhoel
               else
                 forceloc=zero
               end if
!            PAW: V_PAW^prime-V_L^prime
             else
               rr(1)=sqrt(r2)

               if (rr(1)>=rzero) then
                 call paw_splint_der(msz,rad,vloc,d2vloc,1,rr,dvpawdr)
                 call calcdVloc_wvl(dvhgh,rr(1),rloc,charge)
                 forceloc=rhoel*rloc2*(dvpawdr(1)-dvhgh)/rr(1)
               else
                 forceloc=rhoel*rloc2*dvloc_zero_wvl(charge,rloc,msz,rad,vloc,d2vloc)
               end if
             end if

             fxgau=fxgau+forceloc*x
             fygau=fygau+forceloc*y
             fzgau=fzgau+forceloc*z

!            Long range part: error function
             vel=pot(ind)
             fxerf=fxerf+xp*vel*x
             fyerf=fyerf+xp*vel*y
             fzerf=fzerf+xp*vel*z

           else if ((.not.goz).and.(nloc>0)) then
             arg=r2/rloc2
             tt=cprime(nloc)
             do iloc=nloc-1,1,-1
               tt=arg*tt+cprime(iloc)
             end do
             forceleaked=forceleaked+prefactor*xp*tt*rho(1)
           end if

         end do ! i1
       end do ! i2
     end do ! i3
   end if ! n3pi>0

!  Final result of the forces
   floc(1,iat)=(hxh*hyh*hzh*prefactor)*fxerf+(hxh*hyh*hzh/rloc2)*fxgau
   floc(2,iat)=(hxh*hyh*hzh*prefactor)*fyerf+(hxh*hyh*hzh/rloc2)*fygau
   floc(3,iat)=(hxh*hyh*hzh*prefactor)*fzerf+(hxh*hyh*hzh/rloc2)*fzgau

 end do ! iat

 forceleaked=forceleaked*hxh*hyh*hzh
 if (iproc.eq.0) write(std_out,'(a,1pe12.5)') 'done. Leaked force: ',forceleaked

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) i3s,iproc,n1,n1i,n2,n2i,n3,n3pi,natom,hxh,hyh,hzh,&
& rxyz(1,1),floc(1,1),rho(1),pot(1),wvl%h(1)
#endif

 CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* mklocl_wavelets/calcdVloc_wvl
!! NAME
!!  calcdVloc_wvl
!!
!! FUNCTION
!!  Compute 1st-derivative of long-range HGH local ionic potential (derf)
!!
!! COPYRIGHT
!! Copyright (C) 2016-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mklocl_wavelets
!!
!! CHILDREN
!!      calcdvloc_wvl,derf_ab,paw_splint_der
!!
!! SOURCE

subroutine calcdVloc_wvl(yy,xx,rloc,Z)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calcdVloc_wvl'
 use interfaces_43_wvl_wrappers
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 real(dp),intent(in)  :: xx,rloc,Z
 real(dp),intent(out) :: yy

!Local variables-------------------------------
 !scalars
 real(dp):: arg,tt

! *************************************************************************

   arg=xx/(sqrt(2._dp)*rloc)
   call derf_ab(tt,arg)
   yy=(Z/(xx**2))* ( tt - 2._dp/sqrt(pi)*arg*exp(-arg**2) )

 end subroutine calcdVloc_wvl
!!***

!----------------------------------------------------------------------

!!****f* mklocl_wavelets/dvloc_zero_wvl
!! NAME
!!  dvloc_zero_wvl
!!
!! FUNCTION
!!  Use a quadratic interpolation to get limit of (1/x).dVloc(x)/dx at x->0
!!
!! COPYRIGHT
!! Copyright (C) 2013-2016 ABINIT group (TRangel,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

function dvloc_zero_wvl(charge,rloc,msz,rad,vloc,d2vloc)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvloc_zero_wvl'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msz
 real(dp) :: dvloc_zero_wvl
 real(dp),intent(in)  :: charge,rloc
!arrays
 real(dp) :: rad(msz),vloc(msz),d2vloc(msz)

!Local variables-------------------------------
!scalars
 real(dp) :: y1,y2,y3,zz=0._dp
!arrays
 real(dp) :: ll(3),xx(3),yy(3)

! *************************************************************************

!Select 3 points x1,x2,x3 near 0
   if (rad(1)>1.d-10) then
     xx(1:3)=rad(1:3)
   else
     xx(1:3)=rad(2:4)
   end if

!Find the corresponding values of y=(V^PAW(x)-V^HGH(x))/x
   call paw_splint_der(msz,rad,vloc,d2vloc,3,xx,yy)
   call calcdVloc_wvl(y1,xx(1),rloc,charge)
   call calcdVloc_wvl(y2,xx(2),rloc,charge)
   call calcdVloc_wvl(y3,xx(3),rloc,charge)
   yy(1)=(yy(1)-y1)/xx(1)
   yy(2)=(yy(2)-y2)/xx(2)
   yy(3)=(yy(3)-y3)/xx(3)

!Find a polynomial of the form (z=0):
!P(z) = y1.L1(z) + y2.L2(z) + y3.L3(z)

!L1(z) = (z-x2)(z-x3)/((x1-x2)(x1-x3))
   ll(1)=(zz-xx(2))*(zz-xx(3))/((xx(1)-xx(2))*(xx(1)-xx(3)))
!L2(z) = (z-x1)(z-x3)/((x2-x1)(x2-x3))
   ll(2)=(zz-xx(1))*(zz-xx(3))/((xx(2)-xx(1))*(xx(2)-xx(3)))
!L3(z) = (z-x1)(z-x2)/((x3-x1)(x3-x2))
   ll(3)=(zz-xx(1))*(zz-xx(2))/((xx(3)-xx(1))*(xx(3)-xx(2)))

   dvloc_zero_wvl=yy(1)*ll(1)+yy(2)*ll(2)+yy(3)*ll(3)

 end function dvloc_zero_wvl
!!***

end subroutine local_forces_wvl
!!***
