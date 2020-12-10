!!****m* ABINIT/m_dfpt_rhotov
!! NAME
!!  m_dfpt_rhotov
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2020 ABINIT group (XG, DRH, MT, SPr)
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

module m_dfpt_rhotov

 use defs_basis
 use m_abicore
 use m_errors
 use m_cgtools

 use defs_abitypes, only : MPI_type
 use m_time,        only : timab
 use m_spacepar,    only : hartrestr, hartre
 use m_dfpt_mkvxc,    only : dfpt_mkvxc, dfpt_mkvxc_noncoll
 use m_dfpt_mkvxcstr, only : dfpt_mkvxcstr

 implicit none

 private
!!***

 public :: dfpt_rhotov
!!***

contains
!!***

!!****f* ABINIT/dfpt_rhotov
!! NAME
!! dfpt_rhotov
!!
!! FUNCTION
!! This routine is called to compute, from a given 1st-order total density
!!   - the trial (local) 1st-order potential and/or the residual potential,
!!   - some contributions to the 2nd-order energy
!!
!! INPUTS
!!  cplex: if 1, real space 1-order WF on FFT grid are REAL; if 2, COMPLEX
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  idir=direction of atomic displacement (=1,2 or 3 : displacement of atom ipert along the 1st, 2nd or 3rd axis).
!!  ipert=type of the perturbation
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nfft,nkxc)=exchange-correlation kernel
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*nhatdim)= -PAW only- compensation density
!!  nhat1(cplex*nfft,2nspden*usepaw)= -PAW only- 1st-order compensation density
!!  nhat1gr(cplex*nfft,nspden,3*nhat1grdim)= -PAW only- gradients of 1st-order compensation density
!!  nhat1grdim= -PAW only- 1 if nhat1gr array is used ; 0 otherwise
!!  nkxc=second dimension of the array kxc, see rhotoxc.f for a description
!!  non_magnetic_xc= if true, handle density/potential as non-magnetic (even if it is)
!!  nspden=number of spin-density components
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used
!!  optene=0: the contributions to the 2nd order energy are not computed
!!         1: the contributions to the 2nd order energy are computed
!!  optres=0: the trial potential residual is computed ; the input potential value is kept
!!         1: the new value of the trial potential is computed in place of the input value
!!  qphon(3)=reduced coordinates for the phonon wavelength
!!  rhog(2,nfft)=array for Fourier transform of GS electron density
!!  rhog1(2,nfft)=RF electron density in reciprocal space
!!  rhor(nfft,nspden)=array for GS electron density in electrons/bohr**3.
!!  rhor1(cplex*nfft,nspden)=RF electron density in real space (electrons/bohr**3).
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol=unit cell volume in ($\textrm{bohr}^{3}$)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vpsp1(cplex*nfft)=first-order derivative of the ionic potential
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density, see n3xccc
!!
!! OUTPUT
!!  vhartr1(cplex*nfft)=1-order Hartree potential (not output if size=0)
!!  vxc1(cplex*nfft,nspden)= 1st-order XC potential (not output if size=0)
!!  ==== if optene==1
!!    ehart01=inhomogeneous 1st-order Hartree part of 2nd-order total energy
!!    ehart1=1st-order Hartree part of 2nd-order total energy
!!    exc1=1st-order exchange-correlation part of 2nd-order total energy
!!    elpsp1=1st-order local pseudopot. part of 2nd-order total energy.
!!  ==== if optres==0
!!    vresid1(cplex*nfft,nspden)=potential residual
!!    vres2=square of the norm of the residual
!!
!! SIDE EFFECTS
!!  ==== if optres==1
!!    vtrial1(cplex*nfft,nspden)= new value of 1st-order trial potential
!!
!! PARENTS
!!      m_dfpt_lw,m_dfpt_scfcv,m_dfptnl_loop
!!
!! CHILDREN
!!
!! SOURCE

 subroutine dfpt_rhotov(cplex,ehart01,ehart1,elpsp1,exc1,elmag1,gsqcut,idir,ipert,&
&           ixc,kxc,mpi_enreg,natom,nfft,ngfft,nhat,nhat1,nhat1gr,nhat1grdim,nkxc,nspden,n3xccc,&
&           non_magnetic_xc,optene,optres,qphon,rhog,rhog1,rhor,rhor1,rprimd,ucvol,&
&           usepaw,usexcnhat,vhartr1,vpsp1,vresid1,vres2,vtrial1,vxc,vxc1,xccc3d1,ixcrot)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,ixc,n3xccc,natom,nfft,nhat1grdim,nkxc,nspden
 integer,intent(in) :: optene,optres,usepaw,usexcnhat,ixcrot
 logical,intent(in) :: non_magnetic_xc
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(inout) :: ehart01
 real(dp),intent(out) :: vres2
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in) :: kxc(nfft,nkxc)
 real(dp),intent(in) :: vxc(nfft,nspden)
 real(dp),intent(in) :: nhat(nfft,nspden)
 real(dp),intent(in) :: nhat1(cplex*nfft,nspden)  !vz_d
 real(dp),intent(in) :: nhat1gr(cplex*nfft,nspden,3*nhat1grdim)
 real(dp),intent(in) :: qphon(3),rhog(2,nfft)
 real(dp),intent(in) :: rhog1(2,nfft)
 real(dp),target,intent(in) :: rhor(nfft,nspden),rhor1(cplex*nfft,nspden)
 real(dp),intent(in) :: rprimd(3,3),vpsp1(cplex*nfft)
 real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 real(dp),intent(inout) :: vtrial1(cplex*nfft,nspden),elpsp1,ehart1,exc1,elmag1
 real(dp),intent(out) :: vresid1(cplex*nfft,nspden)
 real(dp),target,intent(out) :: vhartr1(:),vxc1(:,:)

!Local variables-------------------------------
!scalars
 integer :: ifft,ispden,nfftot,option
 integer :: optnc,nkxc_cur
 logical :: vhartr1_allocated,vxc1_allocated
 real(dp) :: doti,elpsp10
!arrays
 integer,intent(in)   :: ngfft(18)
 real(dp)             :: tsec(20)
 real(dp),parameter   :: dummyvgeo(3)=zero
 real(dp),allocatable :: rhor1_nohat(:,:),vhartr01(:),vxc1val(:,:)
 real(dp),pointer     :: rhor1_(:,:),vhartr1_(:),vxc1_(:,:),v1zeeman(:,:)

! *********************************************************************

 call timab(157,1,tsec)

 !FR EB SPr
 if (nspden==4) then
   if(usepaw==1) then
     ABI_ERROR('DFPT with nspden=4 works only for norm-conserving psp!')
   end if
 end if

!Get size of FFT grid
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)

!Eventually allocate temporary memory space
 vhartr1_allocated=(size(vhartr1)>0)
 if (vhartr1_allocated) then
   vhartr1_ => vhartr1
 else
   ABI_MALLOC(vhartr1_,(cplex*nfft))
 end if
 vxc1_allocated=(size(vxc1)>0)
 if (vxc1_allocated) then
   vxc1_ => vxc1
 else
   ABI_MALLOC(vxc1_,(cplex*nfft,nspden))
 end if

!If needed, store pseudo density without charge compensation
 if (usepaw==1.and.usexcnhat==0) then
   ABI_MALLOC(rhor1_,(cplex*nfft,nspden))
   rhor1_(:,:)=rhor1(:,:)-nhat1(:,:)
 else
   rhor1_ => rhor1
 end if


 if(ipert==natom+5)then
   ABI_MALLOC(v1zeeman,(cplex*nfft,nspden))
   call dfpt_v1zeeman(nspden,nfft,cplex,idir,v1zeeman)
 end if

!------ Compute 1st-order Hartree potential (and energy) ----------------------

 call hartre(cplex,gsqcut,3,0,mpi_enreg,nfft,ngfft,1,zero,rhog1,rprimd,dummyvgeo,vhartr1_,qpt=qphon)

 if (optene>0) then
   call dotprod_vn(cplex,rhor1,ehart1,doti,nfft,nfftot,1,1,vhartr1_,ucvol)
 end if

 if (optene>0) ehart01=zero
 if(ipert==natom+3 .or. ipert==natom+4) then
   ABI_MALLOC(vhartr01,(cplex*nfft))
   call hartrestr(gsqcut,idir,ipert,mpi_enreg,natom,nfft,ngfft,rhog,rprimd,vhartr01)
   if (optene>0) then
     call dotprod_vn(cplex,rhor1,ehart01,doti,nfft,nfftot,1,1,vhartr01,ucvol)
     ehart01=two*ehart01
     ehart1=ehart1+ehart01
   end if
!  Note that there is a factor 2.0_dp difference with the similar GS formula
   vhartr1_(:)=vhartr1_(:)+vhartr01(:)

   ABI_FREE(vhartr01)
 end if

!------ Compute 1st-order XC potential (and energy) ----------------------
!(including the XC core correction)

!Compute Vxc^(1) (with or without valence contribution according to options)
 option=0;if (optene==0) option=1
 if(ipert==natom+3.or.ipert==natom+4) then
   call dfpt_mkvxcstr(cplex,idir,ipert,kxc,mpi_enreg,natom,nfft,ngfft,nhat,&
&   nhat1,nkxc,non_magnetic_xc,nspden,n3xccc,option,qphon,rhor,rhor1,rprimd,&
&   usepaw,usexcnhat,vxc1_,xccc3d1)
 else
! FR EB non-collinear magnetism
! the second nkxc should be nkxc_cur (see 67_common/nres2vres.F90)
   if (nspden==4) then
     optnc=1
     nkxc_cur=nkxc ! TODO: remove nkxc_cur?

     call dfpt_mkvxc_noncoll(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat,usepaw,nhat1,usepaw,nhat1gr,nhat1grdim,nkxc,&
&     non_magnetic_xc,nspden,n3xccc,optnc,option,qphon,rhor,rhor1,rprimd,usexcnhat,vxc,vxc1_,xccc3d1,ixcrot=ixcrot)

   else
     call dfpt_mkvxc(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1,usepaw,nhat1gr,nhat1grdim,nkxc,&
&     non_magnetic_xc,nspden,n3xccc,option,qphon,rhor1,rprimd,usexcnhat,vxc1_,xccc3d1)
   end if !nspden==4
 end if

!Compute local contribution to 2nd-order energy (includes Vxc and Vpsp and Vmag)
 if (optene>0) then
   if (usepaw==0) then
     call dotprod_vn(cplex,rhor1,elpsp10,doti,nfft,nfftot,nspden,1,vxc1_,ucvol)
     call dotprod_vn(cplex,rhor1,elpsp1 ,doti,nfft,nfftot,1     ,1,vpsp1,ucvol)
     if (ipert==natom+5) then
       call dotprod_vn(cplex,rhor1,elmag1 ,doti,nfft,nfftot,nspden,1,v1zeeman,ucvol)
     end if
   else
     if (usexcnhat/=0) then
       ABI_MALLOC(rhor1_nohat,(cplex*nfft,1))
       rhor1_nohat(:,1)=rhor1(:,1)-nhat1(:,1)
       call dotprod_vn(cplex,rhor1      ,elpsp10,doti,nfft,nfftot,nspden,1,vxc1_,ucvol)
       call dotprod_vn(cplex,rhor1_nohat,elpsp1 ,doti,nfft,nfftot,1     ,1,vpsp1,ucvol)
       ABI_FREE(rhor1_nohat)
     else
       call dotprod_vn(cplex,rhor1_,elpsp10,doti,nfft,nfftot,nspden,1,vxc1_,ucvol)
       call dotprod_vn(cplex,rhor1_,elpsp1 ,doti,nfft,nfftot,1     ,1,vpsp1,ucvol)
     end if
   end if

!  Note that there is a factor 2 difference with the similar GS formula
   elpsp1=two*(elpsp1+elpsp10)
 end if


!Compute XC valence contribution exc1 and complete eventually Vxc^(1)
 if (optene>0) then
   ABI_MALLOC(vxc1val,(cplex*nfft,nspden))
   vxc1val=zero
   option=2
!FR SPr EB non-collinear magnetism
   if (nspden==4) then
     optnc=1
     nkxc_cur=nkxc
     call dfpt_mkvxc_noncoll(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat,usepaw,nhat1,usepaw,nhat1gr,nhat1grdim,nkxc,&
&     non_magnetic_xc,nspden,n3xccc,optnc,option,qphon,rhor,rhor1,rprimd,usexcnhat,vxc,vxc1val,xccc3d1,ixcrot=ixcrot)
   else
     call dfpt_mkvxc(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1,usepaw,nhat1gr,nhat1grdim,nkxc,&
&     non_magnetic_xc,nspden,n3xccc,option,qphon,rhor1,rprimd,usexcnhat,vxc1val,xccc3d1)
   end if !nspden==4
   vxc1_(:,:)=vxc1_(:,:)+vxc1val(:,:)
   call dotprod_vn(cplex,rhor1_,exc1,doti,nfft,nfftot,nspden,1,vxc1val,ucvol)
   ABI_FREE(vxc1val)
 end if

 if (usepaw==1.and.usexcnhat==0) then
   ABI_FREE(rhor1_)
 end if

!DEBUG (do not take away)
!Compute NSC energy ensc1 associated with rhor1 in vtrial1, for debugging purposes
!call dotprod_vn(cplex,rhor1,ensc1,doti,nfft,nfftot,nspden,1,vtrial1,ucvol)
!write(std_out,*)' ek0+eeig0+eloc0=',ek0+eeig0+eloc0
!write(std_out,*)' ensc1=',ensc1
!Compute NSC energy associated with vtrial1, for debugging purposes
!call dotprod_vn(cplex,rhor1,ensc1,doti,mpi_enreg,nfft,nfftot,nspden,1,vtrial1,ucvol)
!ensc1=ensc1+half*enl1
!write(std_out,*)' dfpt_rhotov : check NSC energy, diff=',&
!&  ek0+edocc+eeig0+eloc0+enl0+ensc1
!write(std_out,*)' evarNSC=',ek0+edocc+eeig0+eloc0+enl0
!write(std_out,*)' ensc1,exc1=',ensc1,exc1
!ENDDEBUG

!Here, vhartr1 contains Hartree potential, vpsp1 contains local psp,
!while vxc1 contain xc potential

!------ Produce residual vector and square of norm of it -------------
!(only if requested ; if optres==0)
 if (optres==0) then
!$OMP PARALLEL DO COLLAPSE(2)
   do ispden=1,min(nspden,2)
     do ifft=1,cplex*nfft
       vresid1(ifft,ispden)=vhartr1_(ifft)+vxc1_(ifft,ispden)+vpsp1(ifft)-vtrial1(ifft,ispden)
     end do
   end do
   if(nspden==4)then
!$OMP PARALLEL DO COLLAPSE(2)
     do ispden=3,4
       do ifft=1,cplex*nfft
         vresid1(ifft,ispden)=vxc1_(ifft,ispden)-vtrial1(ifft,ispden)
       end do
     end do
   end if

   if (ipert==natom+5) then
     vresid1 = vresid1 + v1zeeman
   end if
!  Compute square norm vres2 of potential residual vresid
   call sqnorm_v(cplex,nfft,vres2,nspden,optres,vresid1)

 else

!  ------ Produce new value of trial potential-------------
!  (only if requested ; if optres==1)

!$OMP PARALLEL DO COLLAPSE(2)
   do ispden=1,min(nspden,2)
     do ifft=1,cplex*nfft
       vtrial1(ifft,ispden)=vhartr1_(ifft)+vxc1_(ifft,ispden)+vpsp1(ifft)
     end do
   end do
   if(nspden==4)then
!$OMP PARALLEL DO COLLAPSE(2)
     do ispden=3,4
       do ifft=1,cplex*nfft
         vtrial1(ifft,ispden)=vxc1_(ifft,ispden)
       end do
     end do
   end if

   if (ipert==natom+5) then
     vtrial1 = vtrial1 + v1zeeman
   end if

 end if

!Release temporary memory space
 if (.not.vhartr1_allocated) then
   ABI_FREE(vhartr1_)
 end if
 if (.not.vxc1_allocated) then
   ABI_FREE(vxc1_)
 end if

 if (ipert==natom+5) then
   ABI_FREE(v1zeeman)
 end if

 call timab(157,2,tsec)

end subroutine dfpt_rhotov
!!***

!!****f* ABINIT/dfpt_v1zeeman
!! NAME
!!  dfpt_v1zeeman
!!
!! FUNCTION
!!  Calculate 1st order Zeeman potential = -vec{\sigma}.\vec{b}, where
!!  sigma is the vector of Pauli matrices and \vec{b} is the unit
!!  vector indicating the perturbing field direction.
!!
!! INPUTS
!!  nspden = number of density matrix components
!!  nfft   = numbder of fft grid points
!!  cplex  = complex or real density matrix
!!  idir   = direction of the perturbing field in Cartesian frame
!!           1: along x
!!           2: along y
!!           3: along z
!!           4: identity matrix at each fft point is returned (for density-density response)
!!
!! OUTPUT
!!  v1zeeman(nfft*cplex,nspden)= 1st order Zeeman potential, or Identity matrix (electrostatic potential) for idir=4
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  The definition of components of the potential matrix differ depending on cplex
!!  for nspden=4:
!!  For cplex=1, the potential is defined as (V_upup,V_dndn,Re[V_updn],Im[V_updn])
!!  For cplex=2, the definition is (V_upup,V_dndn,V_updn,i.V_updn)
!!
!! PARENTS
!!      m_dfpt_rhotov
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfpt_v1zeeman(nspden,nfft,cplex,idir,v1zeeman)

!Arguments ------------------------------------
 integer , intent(in)    :: idir,nfft,cplex,nspden
 real(dp), intent(inout) :: v1zeeman(cplex*nfft,nspden)

!Local variables-------------------------------
 integer :: ifft
!character(len=500) :: msg

! *************************************************************************

 DBG_ENTER("COLL")

! if (option/=1 .and. option/=2 ) then
!   write(msg,'(3a,i0)')&
!&   'The argument option should be 1 or 2,',ch10,&
!&   'however, option=',option
!   ABI_BUG(msg)
! end if
!
! if (sizein<1) then
!   write(msg,'(3a,i0)')&
!&   'The argument sizein should be a positive number,',ch10,&
!&   'however, sizein=',sizein
!   ABI_ERROR(msg)
! end if

 DBG_EXIT("COLL")

 select case(cplex)
 case(1)
   if (nspden==4) then
     if(idir==3)then       ! Zeeman field along the 3rd axis (z)
       v1zeeman(:,1)=-0.5d0
       v1zeeman(:,2)=+0.5d0
       v1zeeman(:,3)= 0.0d0
       v1zeeman(:,4)= 0.0d0
     else if(idir==2)then  ! Zeeman field along the 2nd axis (y)
       v1zeeman(:,1)= 0.0d0
       v1zeeman(:,2)= 0.0d0
       v1zeeman(:,3)= 0.0d0
       v1zeeman(:,4)=+0.5d0
     else                  ! Zeeman field along the 1st axis (x)
       v1zeeman(:,1)= 0.0d0
       v1zeeman(:,2)= 0.0d0
       v1zeeman(:,3)=-0.5d0
       v1zeeman(:,4)= 0.0d0
     end if
   else if (nspden==2) then
     v1zeeman(:,1)=-0.5e0
     v1zeeman(:,2)= 0.5e0
   else
     v1zeeman(:,1)= 0.0e0
   end if
 case(2)
   if (nspden==2) then
     do ifft=1,nfft
       v1zeeman(2*ifft-1,1)  =-0.5e0
       v1zeeman(2*ifft  ,1)  = 0.0e0
       v1zeeman(2*ifft-1,2)  = 0.5e0
       v1zeeman(2*ifft  ,2)  = 0.0e0
     end do
   else if (nspden==4) then
     select case(idir)
     case(1) !along x, v1=-sigma_x
       do ifft=1,nfft
         v1zeeman(2*ifft-1,1)= 0.0e0 !Re[V^11]
         v1zeeman(2*ifft  ,1)= 0.0e0 !Im[V^11]
         v1zeeman(2*ifft-1,2)= 0.0e0 !Re[V^22]
         v1zeeman(2*ifft  ,2)= 0.0e0 !Im[V^22]
         v1zeeman(2*ifft-1,3)=-0.5e0 !Re[V^12]
         v1zeeman(2*ifft  ,3)= 0.0e0 !Im[V^12]
         v1zeeman(2*ifft-1,4)= 0.0e0 !Re[i.V^21]=Im[V^12]
         v1zeeman(2*ifft  ,4)=-0.5e0 !Im[i.V^21]=Re[V^12]
       end do
     case(2) !along y, v1 = -sigma_y
       do ifft=1,nfft
         v1zeeman(2*ifft-1,1)= 0.0e0 !Re[V^11]
         v1zeeman(2*ifft  ,1)= 0.0e0 !Im[V^11]
         v1zeeman(2*ifft-1,2)= 0.0e0 !Re[V^22]
         v1zeeman(2*ifft  ,2)= 0.0e0 !Im[V^22]
         v1zeeman(2*ifft-1,3)= 0.0e0 !Re[V^12]
         v1zeeman(2*ifft  ,3)=+0.5e0 !Im[V^12]
         v1zeeman(2*ifft-1,4)=+0.5e0 !Re[i.V^21]=Im[V^12]
         v1zeeman(2*ifft  ,4)= 0.0e0 !Im[i.V^21]=Re[V^12]
       end do
     case(3)
       do ifft=1,nfft
         v1zeeman(2*ifft-1,1)=-0.5e0 !Re[V^11]
         v1zeeman(2*ifft  ,1)= 0.0e0 !Im[V^11]
         v1zeeman(2*ifft-1,2)= 0.5e0 !Re[V^22]
         v1zeeman(2*ifft  ,2)= 0.0e0 !Im[V^22]
         v1zeeman(2*ifft-1,3)= 0.0e0 !Re[V^12]
         v1zeeman(2*ifft  ,3)= 0.0e0 !Im[V^12]
         v1zeeman(2*ifft-1,4)= 0.0e0 !Re[i.V^21]
         v1zeeman(2*ifft  ,4)= 0.0e0 !Im[i.V^21]
       end do
     end select
   end if
 end select !cplex

end subroutine dfpt_v1zeeman
!!***

end module m_dfpt_rhotov
!!***
