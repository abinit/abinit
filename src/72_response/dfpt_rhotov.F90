!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_rhotov
!!
!! NAME
!! dfpt_rhotov
!!
!! FUNCTION
!! This routine is called to compute, from a given 1st-order total density
!!   - the trial (local) 1st-order potential and/or the residual potential,
!!   - some contributions to the 2nd-order energy
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (XG, DRH, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex: if 1, real space 1-order WF on FFT grid are REAL; if 2, COMPLEX
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  idir=direction of atomic displacement (=1,2 or 3 : displacement of atom ipert along the 1st, 2nd or 3rd axis).
!!  ipert=type of the perturbation
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nfft,nkxc)=exchange-correlation kernel
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*nhatdim)= -PAW only- compensation density
!!  nhat1(cplex*nfft,2nspden*usepaw)= -PAW only- 1st-order compensation density
!!  nhat1gr(cplex*nfft,nspden,3*nhat1grdim)= -PAW only- gradients of 1st-order compensation density
!!  nhat1grdim= -PAW only- 1 if nhat1gr array is used ; 0 otherwise
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  nspden=number of spin-density components
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used
!!  optene=0: the contributions to the 2nd order energy are not computed
!!         1: the contributions to the 2nd order energy are computed
!!  optres=0: the trial potential residual is computed ; the input potential value is kept
!!         1: the new value of the trial potential is computed in place of the input value
!!  paral_kgb=flag controlling (k,g,bands) parallelization
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
!!      dfpt_scfcv
!!
!! CHILDREN
!!      dfpt_mkvxc,dfpt_mkvxc_noncoll,dfpt_mkvxcstr,dotprod_vn,hartre,hartrestr
!!      sqnorm_v,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


 subroutine dfpt_rhotov(cplex,ehart01,ehart1,elpsp1,exc1,gmet,gprimd,gsqcut,idir,ipert,&
&           ixc,kxc,mpi_enreg,natom,nfft,ngfft,nhat,nhat1,nhat1gr,nhat1grdim,nkxc,nspden,n3xccc,&
&           optene,optres,paral_kgb,qphon,rhog,rhog1,rhor,rhor1,rprimd,ucvol,&
&           usepaw,usexcnhat,vhartr1,vpsp1,vresid1,vres2,vtrial1,vxc1,xccc3d1)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_cgtools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_rhotov'
 use interfaces_18_timing
 use interfaces_53_spacepar
 use interfaces_56_xc
 use interfaces_72_response, except_this_one => dfpt_rhotov
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,ixc,n3xccc,natom,nfft,nhat1grdim,nkxc,nspden
 integer,intent(in) :: optene,optres,paral_kgb,usepaw,usexcnhat
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(inout) :: ehart01 !vz_i
 real(dp),intent(out) :: vres2
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),kxc(nfft,nkxc)
 real(dp),intent(in) :: nhat(nfft,nspden)
 real(dp),intent(in) :: nhat1(cplex*nfft,nspden)  !vz_d
 real(dp),intent(in) :: nhat1gr(cplex*nfft,nspden,3*nhat1grdim)
 real(dp),intent(in) :: qphon(3),rhog(2,nfft)
 real(dp),intent(in) :: rhog1(2,nfft)
 real(dp),target,intent(in) :: rhor(nfft,nspden),rhor1(cplex*nfft,nspden)
 real(dp),intent(in) :: rprimd(3,3),vpsp1(cplex*nfft)
 real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 real(dp),intent(inout) :: vtrial1(cplex*nfft,nspden),elpsp1,ehart1,exc1 !vz_d
 real(dp),intent(out) :: vresid1(cplex*nfft,nspden)
 real(dp),target,intent(out) :: vhartr1(:),vxc1(:,:)

!Local variables-------------------------------
!scalars
 integer :: ifft,ispden,nfftot,option
 integer :: optnc,optxc,nkxc_cur
 logical :: vhartr1_allocated,vxc1_allocated
 real(dp) :: doti,elpsp10
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp) :: tsec(20)
 real(dp),allocatable :: rhor1_nohat(:,:),vhartr01(:),vxc1val(:,:)
 real(dp),pointer :: rhor1_(:,:),vhartr1_(:),vxc1_(:,:)

! *********************************************************************

 call timab(157,1,tsec)

 !FR EB
 if (nspden==4) then
   MSG_WARNING('DFPT under development for nspden=4!')
 end if

!Get size of FFT grid
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)

!Eventually allocate temporary memory space
 vhartr1_allocated=(size(vhartr1)>0)
 if (vhartr1_allocated) then
   vhartr1_ => vhartr1
 else
   ABI_ALLOCATE(vhartr1_,(cplex*nfft))
 end if
 vxc1_allocated=(size(vxc1)>0)
 if (vxc1_allocated) then
   vxc1_ => vxc1
 else
   ABI_ALLOCATE(vxc1_,(cplex*nfft,nspden))
 end if

!If needed, store pseudo density without charge compensation
 if (usepaw==1.and.usexcnhat==0) then
   ABI_ALLOCATE(rhor1_,(cplex*nfft,nspden))
   rhor1_(:,:)=rhor1(:,:)-nhat1(:,:)
 else
   rhor1_ => rhor1
 end if

!------ Compute 1st-order Hartree potential (and energy) ----------------------

 call hartre(cplex,gmet,gsqcut,0,mpi_enreg,nfft,ngfft,paral_kgb,qphon,rhog1,vhartr1_)

 if (optene>0) then
   call dotprod_vn(cplex,rhor1,ehart1,doti,nfft,nfftot,1,1,vhartr1_,ucvol)
 end if

 if (optene>0) ehart01=zero
 if(ipert==natom+3 .or. ipert==natom+4) then
   ABI_ALLOCATE(vhartr01,(cplex*nfft))
   call hartrestr(gmet,gprimd,gsqcut,idir,ipert,mpi_enreg,natom,nfft,ngfft,paral_kgb,rhog,vhartr01)
   if (optene>0) then
     call dotprod_vn(cplex,rhor1,ehart01,doti,nfft,nfftot,1,1,vhartr01,ucvol)
     ehart01=two*ehart01
     ehart1=ehart1+ehart01
   end if
!  Note that there is a factor 2.0_dp difference with the similar GS formula
   vhartr1_(:)=vhartr1_(:)+vhartr01(:)
   ABI_DEALLOCATE(vhartr01)
 end if

!------ Compute 1st-order XC potential (and energy) ----------------------
!(including the XC core correction)

!Compute Vxc^(1) (with or without valence contribution according to options)
 option=0;if (optene==0) option=1
 if(ipert==natom+3.or.ipert==natom+4) then
   call dfpt_mkvxcstr(cplex,idir,ipert,kxc,mpi_enreg,natom,nfft,ngfft,nhat,&
&   nhat1,nkxc,nspden,n3xccc,option,paral_kgb,qphon,rhor,rhor1,rprimd,&
&   usepaw,usexcnhat,vxc1_,xccc3d1)
 else
! FR EB non-collinear magnetism
! the second nkxc should be nkxc_cur (see 67_common/nres2vres.F90)
   if (nspden==4) then
     optnc=1
     optxc=1
     nkxc_cur=nkxc ! TODO: remove nkxc_cur?
!     print*, 'first time call to mkvxc-noncoll: option= ', option
     call dfpt_mkvxc_noncoll(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1,usepaw,nhat1gr,nhat1grdim,nkxc,&
&     nkxc_cur,nspden,n3xccc,optnc,option,optxc,paral_kgb,qphon,rhor,rhor1,rprimd,usexcnhat,vxc1_,xccc3d1)
   else
     call dfpt_mkvxc(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1,usepaw,nhat1gr,nhat1grdim,nkxc,&
&     nspden,n3xccc,option,paral_kgb,qphon,rhor1,rprimd,usexcnhat,vxc1_,xccc3d1)
   end if !nspden==4
 end if

!Compute local contribution to 2nd-order energy (includes Vxc and Vpsp)
 if (optene>0) then
   if (usepaw==0) then
     call dotprod_vn(cplex,rhor1,elpsp10,doti,nfft,nfftot,nspden,1,vxc1_,ucvol)
     call dotprod_vn(cplex,rhor1,elpsp1 ,doti,nfft,nfftot,1     ,1,vpsp1,ucvol)
   else
     if (usexcnhat/=0) then
       ABI_ALLOCATE(rhor1_nohat,(cplex*nfft,1))
       rhor1_nohat(:,1)=rhor1(:,1)-nhat1(:,1)
       call dotprod_vn(cplex,rhor1      ,elpsp10,doti,nfft,nfftot,nspden,1,vxc1_,ucvol)
       call dotprod_vn(cplex,rhor1_nohat,elpsp1 ,doti,nfft,nfftot,1     ,1,vpsp1,ucvol)
       ABI_DEALLOCATE(rhor1_nohat)
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
   ABI_ALLOCATE(vxc1val,(cplex*nfft,nspden))
   option=2
!FR EB non-collinear magnetism
   if (nspden==4) then
     optnc=1
     optxc=1
     nkxc_cur=nkxc
!     print*, 'second time call to mkvxc-noncoll: option= ', option
     call dfpt_mkvxc_noncoll(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1,usepaw,nhat1gr,nhat1grdim,nkxc,&
&     nkxc_cur,nspden,n3xccc,optnc,option,optxc,paral_kgb,qphon,rhor,rhor1,rprimd,usexcnhat,vxc1val,xccc3d1)
   else
     call dfpt_mkvxc(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1,usepaw,nhat1gr,nhat1grdim,nkxc,&
&     nspden,n3xccc,option,paral_kgb,qphon,rhor1,rprimd,usexcnhat,vxc1val,xccc3d1)
   end if !nspden==4
   vxc1_(:,:)=vxc1_(:,:)+vxc1val(:,:)
   call dotprod_vn(cplex,rhor1_,exc1,doti,nfft,nfftot,nspden,1,vxc1val,ucvol)
   ABI_DEALLOCATE(vxc1val)
 end if

 if (usepaw==1.and.usexcnhat==0) then
   ABI_DEALLOCATE(rhor1_)
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
         vresid1(ifft,ispden)=vxc1_(ifft,ispden)
       end do
     end do
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

 end if

!Release temporary memory space
 if (.not.vhartr1_allocated) then
   ABI_DEALLOCATE(vhartr1_)
 end if
 if (.not.vxc1_allocated) then
   ABI_DEALLOCATE(vxc1_)
 end if

 call timab(157,2,tsec)

end subroutine dfpt_rhotov
!!***
