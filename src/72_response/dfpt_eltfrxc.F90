!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_eltfrxc
!! NAME
!! dfpt_eltfrxc
!!
!! FUNCTION
!! Compute the 2nd derivatives of exchange-correlation energy
!! with respect to all pairs of strain and strain-atomic displacement
!! for the frozen wavefunction contribution to the elastic
!! and internal strain tensors
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DRH, DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms ordered by type
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | natom=number of atoms in unit cell
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | nspden=number of spin-density components
!!   | ntypat=number of types of atoms in cell.
!!   | typat(natom)=integer type for each atom in cell
!!  enxc=exchange and correlation energy (hartree)
!!  gsqcut=Fourier cutoff on G^2 for "large sphere" of radius double that of the basis sphere
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngfftf(18)= -PAW ONLY- contain all needed information about 3D FFT for the fine grid
!!              (ngfftf=ngfft for norm-conserving potential runs)
!!  nkxc=2nd dimension of kxc
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of xccc3d (0 if no core charge, nfft otherwise)
!!  nhat(nfft,nspden*nhatdim)= -PAW only- compensation density 
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase (structure factor) information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhor(nfft,nspden)=electron density in r space
!!   (if spin polarized, array contains total density in first half and
!!    spin-up density in second half)
!!   (for non-collinear magnetism, first element: total density,
!!    3 next ones: mx,my,mz)
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!  vxc(nfft,nspden)=xc potential (spin up in first half and spin down in
!!   second half if nspden=2)
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and five derivatives,
!!   for each type of atom, from psp
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced coordinates for atoms in unit cell
!!
!! OUTPUT
!!  eltfrxc(6+3*natom,6) = xc frozen wavefunction contribution to the
!!   elastic tensor
!!
!! SIDE EFFECTS
!!
!! NOTES
!!      Much of the code in versions of this routine prior to 4.4.5
!!      has been transfered to its child eltxccore.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      atm2fft,dfpt_atm2fft,dfpt_mkcore,dfpt_mkvxcstr,dotprod_vn,eltxccore
!!      fourdp,metric,paw_spline,pawpsp_cg,pawrad_free,pawrad_init,pawtab_free
!!      pawtab_nullify,redgr,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_eltfrxc(atindx,dtset,eltfrxc,enxc,gsqcut,kxc,mpi_enreg,mgfft,&
& nattyp,nfft,ngfft,ngfftf,nhat,nkxc,n3xccc,pawtab,ph1d,psps,rhor,rprimd,&
& usexcnhat,vxc,xccc3d,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 
 use m_pawtab,      only : pawtab_type,pawtab_free,pawtab_nullify
 use m_pawrad,      only : pawrad_type,pawrad_init,pawrad_free
 use m_pawpsp,      only : pawpsp_cg
 use m_paw_numeric, only : paw_spline

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_eltfrxc'
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_53_ffts
 use interfaces_53_spacepar
 use interfaces_64_psp
 use interfaces_72_response, except_this_one => dfpt_eltfrxc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
!scalars
 integer,intent(in) :: mgfft,n3xccc,nfft,nkxc,usexcnhat
 real(dp),intent(in) :: enxc,gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,intent(in) :: atindx(dtset%natom),nattyp(dtset%ntypat),ngfft(18)
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(in) :: nhat(nfft,dtset%nspden*psps%usepaw)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom),rprimd(3,3)
 real(dp),intent(in) :: vxc(nfft,dtset%nspden),xccc3d(n3xccc)
 real(dp),intent(in) :: xred(3,dtset%natom)
 real(dp),intent(in),target :: rhor(nfft,dtset%nspden)
 real(dp),intent(inout) :: kxc(nfft,nkxc)
 real(dp),intent(out) :: eltfrxc(6+3*dtset%natom,6)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: mshift=401
 integer :: cplex,fgga,ia,idir,ielt,ieltx,ierr,ifft,ii,ipert,is1,is2,ispden,ispden_c,jj,ka,kb
 integer :: kd,kg,n1,n1xccc,n2,n3,n3xccc_loc,optatm,optdyfr,opteltfr,optgr
 integer :: option,optn,optn2,optstr,optv
 real(dp) :: d2eacc,d2ecdgs2,d2exdgs2,d2gsds1ds2,d2gstds1ds2,decdgs,dexdgs
 real(dp) :: dgsds10,dgsds20,dgstds10,dgstds20,rstep,spnorm,tmp0,tmp0t
 real(dp) :: ucvol,valuei,yp1,ypn
 type(pawrad_type) :: core_mesh
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 real(dp) :: corstr(6),dummy6(0),dummy_in(0,0)
 real(dp) :: dummy_out1(0),dummy_out2(0),dummy_out3(0),dummy_out4(0),dummy_out5(0),dummy_out6(0)
 real(dp) :: eltfrxc_test1(6+3*dtset%natom,6),eltfrxc_test2(6+3*dtset%natom,6)
 real(dp) :: gmet(3,3),gprimd(3,3),qphon(3),rmet(3,3),tsec(2)
 real(dp) :: strn_dummy6(0), strv_dummy6(0)
 real(dp),allocatable :: d2gm(:,:,:,:),dgm(:,:,:),eltfrxc_tmp(:,:)
 real(dp),allocatable :: eltfrxc_tmp2(:,:),elt_work(:,:),rho0_redgr(:,:,:)
 real(dp),allocatable :: vxc10(:,:),vxc10_core(:),vxc10_coreg(:,:)
 real(dp),allocatable :: vxc1is_core(:),vxc1is_coreg(:,:),vxc_core(:)
 real(dp),allocatable :: vxc_coreg(:,:),work(:),workgr(:,:),xccc1d(:,:,:)
 real(dp),allocatable :: xccc3d1(:),xccc3d1_temp(:,:),xcccrc(:)
 real(dp),pointer :: rhor_(:,:)
 type(pawtab_type),allocatable :: pawtab_test(:)

! *************************************************************************

!Initialize variables
 cplex=1
 qphon(:)=zero
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)

 n1xccc = psps%n1xccc
 if(psps%usepaw==0)then
   ABI_ALLOCATE(xcccrc,(dtset%ntypat))
   ABI_ALLOCATE(xccc1d,(n1xccc,6,dtset%ntypat))
   xcccrc = psps%xcccrc
   xccc1d = psps%xccc1d
 end if

 if (usexcnhat==0.and.dtset%usepaw==1) then     
   ABI_ALLOCATE(rhor_,(nfft,dtset%nspden))
   rhor_(:,:) = rhor(:,:)-nhat(:,:)
 else
   rhor_ => rhor
 end if

!HACK - should be fixed globally
 if(n1xccc==0) then
   n3xccc_loc=0
 else
   n3xccc_loc=n3xccc
 end if

 fgga=0 ; if(nkxc==7.or.nkxc==19) fgga=1

 ABI_ALLOCATE(eltfrxc_tmp,(6+3*dtset%natom,6))
 ABI_ALLOCATE(eltfrxc_tmp2,(6+3*dtset%natom,6))
 ABI_ALLOCATE(vxc10,(nfft,dtset%nspden))
 ABI_ALLOCATE(xccc3d1,(cplex*nfft))

 if(n1xccc/=0) then
   ABI_ALLOCATE(vxc_core,(nfft))
   ABI_ALLOCATE(vxc10_core,(nfft))
   ABI_ALLOCATE(vxc1is_core,(nfft))

   if(dtset%nspden==1) then
     vxc_core(:)=vxc(:,1)
   else
     vxc_core(:)=0.5_dp*(vxc(:,1)+vxc(:,2))
   end if
 end if

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 
!For GGA case, prepare quantities needed to evaluate contributions
!arising from the strain dependence of the gradient operator itself

 if(fgga==1) then
   ABI_ALLOCATE(rho0_redgr,(3,nfft,dtset%nspden))
   ABI_ALLOCATE(work,(nfft))
   ABI_ALLOCATE(workgr,(nfft,3))

!  Set up metric tensor derivatives
   ABI_ALLOCATE(dgm,(3,3,6))
   ABI_ALLOCATE(d2gm,(3,3,6,6))
!  Loop over 2nd strain index
   do is2=1,6
     kg=idx(2*is2-1);kd=idx(2*is2)
     do jj = 1,3
       dgm(:,jj,is2)=-(gprimd(kg,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(kg,jj))
     end do

!    Loop over 1st strain index
     do is1=1,6
       ka=idx(2*is1-1);kb=idx(2*is1)
       d2gm(:,:,is1,is2)=0._dp
       do jj = 1,3
         if(ka==kg) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&         +gprimd(kb,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(kb,jj)
         if(ka==kd) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&         +gprimd(kb,:)*gprimd(kg,jj)+gprimd(kg,:)*gprimd(kb,jj)
         if(kb==kg) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&         +gprimd(ka,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(ka,jj)
         if(kb==kd) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&         +gprimd(ka,:)*gprimd(kg,jj)+gprimd(kg,:)*gprimd(ka,jj)
       end do
       d2gm(:,:,is1,is2)=0.5_dp*d2gm(:,:,is1,is2)
     end do
   end do

!  Compute the reduced gradients of the zero-order charge density.
!  Note that in the spin-polarized case, we are computing the reduced
!  gradients of 2 X the spin-up or spin-down charge.  This simplifies
!  subsequent code for the non-spin-polarized case.
   if(dtset%nspden==1) then
     work(:)=rhor_(:,1)
   else
     work(:)=2.0_dp*rhor_(:,2)
   end if
   if(n1xccc/=0) then
     work(:)=work(:)+xccc3d(:)
   end if
   call redgr (work,workgr,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb)
   do ifft=1,nfft
     rho0_redgr(:,ifft,1)=workgr(ifft,:)
   end do
   if(dtset%nspden==2) then
     work(:)=2.0_dp*(rhor_(:,1)-rhor_(:,2))
     if(n1xccc/=0) then
       work(:)=work(:)+xccc3d(:)
     end if
     call redgr(work,workgr,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb)
     do ifft=1,nfft
       rho0_redgr(:,ifft,2)=workgr(ifft,:)
     end do
   end if
   ABI_DEALLOCATE(work)
   ABI_DEALLOCATE(workgr)
 end if !GGA


!Null the elastic tensor accumulator
 eltfrxc(:,:)=zero;eltfrxc_tmp(:,:)=zero;eltfrxc_tmp2(:,:) = zero

!Normalization factor
 if(dtset%nspden==1) then
   spnorm=one
 else
   spnorm=half
 end if

!Big loop over 2nd strain index
 do is2=1,6

!  Translate strain index as needed by dfpt_mkcore below.
   if(is2<=3) then
     ipert=dtset%natom+3
     idir=is2
   else
     ipert=dtset%natom+4
     idir=is2-3
   end if

!  Generate first-order core charge for is2 strain if core charges are present.
   if(n1xccc/=0)then

     if (psps%usepaw==1 .or. psps%nc_xccc_gspace==1) then
!      Calculation in Reciprocal space for paw or NC with nc_xccc_gspace
       ABI_ALLOCATE(xccc3d1_temp,(cplex*nfft,1))
       xccc3d1_temp = zero

       call dfpt_atm2fft(atindx,cplex,gmet,gprimd,gsqcut,is2,ipert,&
&       mgfft,psps%mqgrid_vl,dtset%natom,1,nfft,ngfftf,dtset%ntypat,&
&       ph1d,psps%qgrid_vl,qphon,dtset%typat,ucvol,psps%usepaw,xred,psps,pawtab,&
&       atmrhor1=xccc3d1_temp,optn2_in=1,&
&       comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
&       paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)
       xccc3d1(:) = xccc3d1_temp(:,1)
       ABI_DEALLOCATE(xccc3d1_temp)

     else
!      Calculation in direct space for norm conserving:
       call dfpt_mkcore(cplex,idir,ipert,dtset%natom,dtset%ntypat,n1,n1xccc,&
&       n2,n3,qphon,rprimd,dtset%typat,ucvol,&
&       xcccrc,xccc1d,xccc3d1,xred)
     end if
   else
     xccc3d1(:)=zero
   end if

!  Compute the first-order potentials.
!  Standard first-order potential for LDA and GGA with core charge
   if(fgga==0 .or. (fgga==1 .and. n1xccc/=0)) then
     option=0
     call dfpt_mkvxcstr(cplex,idir,ipert,kxc,mpi_enreg,dtset%natom,nfft,ngfft,nhat,&
&     dummy_in,nkxc,dtset%nspden,n3xccc_loc,option,mpi_enreg%paral_kgb,qphon,rhor,rhor,&
&     rprimd,dtset%usepaw,usexcnhat,vxc10,xccc3d1)
     if(n1xccc/=0)then
       if(dtset%nspden==1) then
         vxc10_core(:)=vxc10(:,1)
         vxc1is_core(:)=vxc10(:,1)
       else
         vxc10_core(:)=0.5_dp*(vxc10(:,1)+vxc10(:,2))
         vxc1is_core(:)=0.5_dp*(vxc10(:,1)+vxc10(:,2))
       end if
     end if
   end if

!  For GGA, first-order potential with doubled gradient operator strain
!  derivative terms needed for elastic tensor but not internal strain.
   if(fgga==1) then
     option=2
     call dfpt_mkvxcstr(cplex,idir,ipert,kxc,mpi_enreg,dtset%natom,nfft,ngfft,nhat,&
&     dummy_in,nkxc,dtset%nspden,n3xccc_loc,option,mpi_enreg%paral_kgb,qphon,rhor,rhor,&
&     rprimd,dtset%usepaw,usexcnhat,vxc10,xccc3d1)
     if(n1xccc/=0)then
       if(dtset%nspden==1) then
         vxc10_core(:)=vxc10(:,1)
       else
         vxc10_core(:)=0.5_dp*(vxc10(:,1)+vxc10(:,2))
       end if
     end if
   end if


!  Additional term for diagonal strains.
   if(is2<=3) then
     vxc10(:,:)=vxc10(:,:)+vxc(:,:)
     if(n1xccc/=0) then
       vxc10_core(:)=vxc10_core(:)+2.0_dp*vxc_core(:)
       vxc1is_core(:)=vxc1is_core(:)+vxc_core(:)
     end if
   end if

!  For GGA, compute the contributions from the strain derivatives acting
!  on the gradient operators.
   if(fgga==1) then

     if (dtset%nspden==1) then
       do ifft=1,nfft
!        Collect the needed derivatives of Exc.  The factors introduced
!        deal with the difference between density as used here and
!        spin density as used with these kxc terms in other contexts.
         dexdgs  =half   *kxc(ifft,2)
         d2exdgs2=quarter*kxc(ifft,4)
!        Loop over 1st strain index
         do is1=1,6
!          The notation here is .gs... for the derivatives of the squared-
!          gradient of (2X) each spin density, and .gst... for the total density.
           dgsds10=zero;dgsds20=zero;d2gsds1ds2=zero
           do jj=1,3
             do ii=1,3
               tmp0=rho0_redgr(ii,ifft,1)*rho0_redgr(jj,ifft,1)
               dgsds10=dgsds10+dgm(ii,jj,is1)*tmp0
               dgsds20=dgsds20+dgm(ii,jj,is2)*tmp0
               d2gsds1ds2=d2gsds1ds2+d2gm(ii,jj,is1,is2)*tmp0
             end do
           end do
!          Volume derivative terms added
           if(is1<=3) d2gsds1ds2=d2gsds1ds2+dgsds20
           if(is2<=3) d2gsds1ds2=d2gsds1ds2+dgsds10
!          Add the gradient derivative terms to eltfrxc.
           eltfrxc(is1,is2)=eltfrxc(is1,is2)+d2exdgs2*dgsds10*dgsds20+dexdgs*d2gsds1ds2
         end do !is1
       end do !ifft

     else ! nspden==2

       do ispden=1,dtset%nspden
         ispden_c=dtset%nspden-ispden+1

         do ifft=1,nfft

!          Collect the needed derivatives of Exc.  The factors introduced
!          deal with the difference between density as used here and
!          spin density as used with these kxc terms in other contexts.
           dexdgs  =quarter       *kxc(ifft,3+ispden)
           d2exdgs2=quarter*eighth*kxc(ifft,7+ispden)
           decdgs  =eighth        *kxc(ifft,10)
           d2ecdgs2=eighth*eighth *kxc(ifft,13)

!          Loop over 1st strain index
           do is1=1,6

!            The notation here is .gs... for the derivatives of the squared-
!            gradient of (2X) each spin density, and .gst... for the total
!            density.  Note the hack that the the total density is given
!            by the same expression for either the non-polarized or spin-
!            polarized case, implemented with the "complementary" index ispden_c
!            in the expression for tmp0t below.
             dgsds10=zero;dgsds20=zero;d2gsds1ds2=zero
             dgstds10=zero;dgstds20=zero;d2gstds1ds2=zero
             do jj=1,3
               do ii=1,3
                 tmp0=rho0_redgr(ii,ifft,ispden)*rho0_redgr(jj,ifft,ispden)
                 tmp0t=(rho0_redgr(ii,ifft,ispden)+rho0_redgr(ii,ifft,ispden_c))&
&                 *(rho0_redgr(jj,ifft,ispden)+rho0_redgr(jj,ifft,ispden_c))
                 dgsds10=dgsds10+dgm(ii,jj,is1)*tmp0
                 dgsds20=dgsds20+dgm(ii,jj,is2)*tmp0
                 dgstds10=dgstds10+dgm(ii,jj,is1)*tmp0t
                 dgstds20=dgstds20+dgm(ii,jj,is2)*tmp0t
                 d2gsds1ds2=d2gsds1ds2+d2gm(ii,jj,is1,is2)*tmp0
                 d2gstds1ds2=d2gstds1ds2+d2gm(ii,jj,is1,is2)*tmp0t
               end do
             end do
!            Volume derivative terms added
             if(is1<=3) then
               d2gsds1ds2=d2gsds1ds2+dgsds20
               d2gstds1ds2=d2gstds1ds2+dgstds20
             end if
             if(is2<=3) then
               d2gsds1ds2=d2gsds1ds2+dgsds10
               d2gstds1ds2=d2gstds1ds2+dgstds10
             end if

!            Add the gradient derivative terms to eltfrxc.
             eltfrxc(is1,is2)=eltfrxc(is1,is2)+spnorm*&
&             (d2exdgs2*(dgsds10 *dgsds20) + dexdgs*d2gsds1ds2&
&             +d2ecdgs2*(dgstds10*dgstds20)+ decdgs*d2gstds1ds2)

           end do !is1
         end do !ifft
       end do !ispden

     end if ! nspden

   end if !GGA

!  Compute valence electron 1st-order charge contributions.  Recall that
!  the diagonal strain derivatives of the valence charge are minus the
!  zero-order density.  The explicit symmetrization avoids the need
!  to store vxc10 for strain indices other than is2.

   call dotprod_vn(1,rhor_,d2eacc,valuei,nfft,nfft,dtset%nspden,1,&
&   vxc10,ucvol)
   do is1=1,3
     eltfrxc_tmp(is1,is2)=eltfrxc_tmp(is1,is2)-0.5_dp*d2eacc
     eltfrxc_tmp(is2,is1)=eltfrxc_tmp(is2,is1)-0.5_dp*d2eacc
   end do

!  Compute additional core contributions from is1 perturbation
!  Internal strain terms calculated here.
   if(n1xccc/=0) then

     if (psps%usepaw==1 .or. psps%nc_xccc_gspace==1) then
!      Calculation in Reciprocal space for paw or NC with nc_xccc_gspace
       optatm=0;optdyfr=0;optgr=0;optstr=0;optv=0;optn=n3xccc/nfft;optn2=1;opteltfr=1
       ABI_ALLOCATE(vxc10_coreg,(2,nfft))
       ABI_ALLOCATE(vxc_coreg,(2,nfft))
       ABI_ALLOCATE(vxc1is_coreg,(2,nfft))

       vxc10_coreg(:,:)=zero;vxc10_coreg(:,:)=zero;vxc1is_coreg(:,:)=zero;

!      Fourier transform of Vxc_core/vxc10_core to use in atm2fft (reciprocal space calculation) 
       call fourdp(1,vxc10_coreg,vxc10_core,-1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,0)
       call fourdp(1,vxc_coreg,vxc_core,-1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,0)
       call fourdp(1,vxc1is_coreg,vxc1is_core,-1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,0)

       call atm2fft(atindx,dummy_out1,dummy_out2,dummy_out3,dummy_out4,eltfrxc_tmp2,dummy_in,gmet,gprimd,&
&       dummy_out5,dummy_out6,gsqcut,mgfft,psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,dtset%ntypat,&
&       optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,psps,pawtab,ph1d,psps%qgrid_vl,dtset%qprtrb,&
&       dummy_in,strn_dummy6,strv_dummy6,ucvol,psps%usepaw,vxc_coreg,vxc10_coreg,vxc1is_coreg,dtset%vprtrb,psps%vlspl,is2_in=is2,&
&       comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
&       paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)

!       The indexing array atindx is used to reestablish the correct order of atoms
       ABI_ALLOCATE(elt_work,(6+3*dtset%natom,6))
       elt_work(1:6,1:6)=eltfrxc_tmp2(1:6,1:6)
       do ia=1,dtset%natom
         ielt=7+3*(ia-1)
         ieltx=7+3*(atindx(ia)-1)
         elt_work(ielt:ielt+2,1:6)=eltfrxc_tmp2(ieltx:ieltx+2,1:6)
       end do
       eltfrxc_tmp2(:,:)=elt_work(:,:)
       ABI_DEALLOCATE(elt_work)


       ABI_DEALLOCATE(vxc10_coreg)
       ABI_DEALLOCATE(vxc_coreg)
       ABI_DEALLOCATE(vxc1is_coreg)
       eltfrxc(:,:)= eltfrxc(:,:) + eltfrxc_tmp2(:,:)

     else

       call eltxccore(eltfrxc,is2,mpi_enreg%my_natom,dtset%natom,nfft,dtset%ntypat,&
&       n1,n1xccc,n2,n3,rprimd,dtset%typat,ucvol,vxc_core,vxc10_core,vxc1is_core,&
&       xcccrc,xccc1d,xred,mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)

!DEBUG
!      TEST ZONE (DO NOT REMOVE) USE TO RECIPROCAL SPACE IN NC CASE
       if (dtset%userid==567) then
         eltfrxc_test1(:,is2)=zero;eltfrxc_test2(:,is2)=zero
         call eltxccore(eltfrxc_test1,is2,mpi_enreg%my_natom,dtset%natom,nfft,dtset%ntypat,&
&         n1,n1xccc,n2,n3,rprimd,dtset%typat,ucvol,vxc_core,vxc10_core,vxc1is_core,&
&         xcccrc,xccc1d,xred,mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
!        if (is2==1) print*,"elt-frxc from eltxccore",is2,eltfrxc_test1(1,1)*ucvol/dble(nfft)
         ABI_DATATYPE_ALLOCATE(pawtab_test,(dtset%ntypat))
         call pawtab_nullify(pawtab_test)
         do jj=1,dtset%ntypat
           pawtab_test(jj)%mqgrid=psps%mqgrid_vl
           ABI_ALLOCATE(pawtab_test(jj)%tcorespl,(pawtab_test(jj)%mqgrid,2))
           rstep=xcccrc(jj)/dble(n1xccc-1)
           call pawrad_init(mesh=core_mesh,mesh_size=n1xccc,mesh_type=1,rstep=rstep)
           call pawpsp_cg(pawtab_test(jj)%dncdq0,pawtab_test(jj)%d2ncdq0,psps%mqgrid_vl,psps%qgrid_vl,&
&           pawtab_test(jj)%tcorespl(:,1),core_mesh,xccc1d(:,1,jj),yp1,ypn)
           call paw_spline(psps%qgrid_vl,pawtab_test(jj)%tcorespl(:,1),psps%mqgrid_vl,yp1,ypn,pawtab_test(jj)%tcorespl(:,2))
!          if (is2==1) then
!            do ii=1,n1xccc;write(100+jj,*) (ii-1)*rstep,xccc1d(ii,1,jj);enddo
!            do ii=1,psps%mqgrid_vl;write(200+jj,*) psps%qgrid_vl(ii),pawtab_test(jj)%tcorespl(ii,1);enddo
!          end if
         end do    
         ABI_ALLOCATE(vxc10_coreg,(2,nfft))
         ABI_ALLOCATE(vxc_coreg,(2,nfft))
         ABI_ALLOCATE(vxc1is_coreg,(2,nfft))
         vxc10_coreg(:,:)=zero;vxc10_coreg(:,:)=zero;vxc1is_coreg(:,:)=zero;
         call fourdp(1,vxc10_coreg,vxc10_core,-1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,0)
         call fourdp(1,vxc_coreg,vxc_core,-1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,0)
         call fourdp(1,vxc1is_coreg,vxc1is_core,-1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,0)   
         optatm=0;optdyfr=0;optgr=0;optstr=0;optv=0;optn=1;optn2=1;opteltfr=1;corstr=zero
         call atm2fft(atindx,dummy_out1,dummy_out2,dummy_out3,dummy_out4,eltfrxc_test2,dummy_in,gmet,gprimd,&
&         dummy_out5,dummy_out6,gsqcut,mgfft,psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,dtset%ntypat,&
&         optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,psps,pawtab_test,ph1d,psps%qgrid_vl,dtset%qprtrb,&
&         dummy_in,corstr,dummy6,ucvol,psps%usepaw,vxc_coreg,vxc10_coreg,vxc1is_coreg,dtset%vprtrb,psps%vlspl,is2_in=is2)
         ABI_DEALLOCATE(vxc10_coreg)
         ABI_DEALLOCATE(vxc_coreg)
         ABI_DEALLOCATE(vxc1is_coreg)
         call pawrad_free(core_mesh)
         call pawtab_free(pawtab_test)
         ABI_DATATYPE_DEALLOCATE(pawtab_test)
         eltfrxc(:,:)= eltfrxc(:,:)+eltfrxc_test2(:,:)
!        if (is2==1) print*,"cor-str from atm2fft",is2,corstr*ucvol
!        if (is2==1) print*,"elt-frxc from atm2fft  ",is2,eltfrxc_test2(1,1)
       end if
!DEBUG

     end if
   end if

!  Additional term for diagonal strains
   if(is2<=3) then
     do is1=1,3
       eltfrxc_tmp(is1,is2)=eltfrxc_tmp(is1,is2)+enxc
     end do
   end if
 end do !is2 outermost strain loop

!Accumulate eltfrxc accross processors
 call timab(48,1,tsec)
 call xmpi_sum(eltfrxc,mpi_enreg%comm_fft,ierr)
 call timab(48,2,tsec)

 !Normalize accumulated 2nd derivatives in NC case
 if(psps%usepaw==1)then
   eltfrxc(:,:)=eltfrxc_tmp(:,:)+eltfrxc
 else
   eltfrxc(:,:)=eltfrxc_tmp(:,:)+eltfrxc*ucvol/dble(nfft)
 end if

 ABI_DEALLOCATE(eltfrxc_tmp)
 ABI_DEALLOCATE(eltfrxc_tmp2)
 ABI_DEALLOCATE(vxc10)
 ABI_DEALLOCATE(xccc3d1)
 if(psps%usepaw==0)then
   ABI_DEALLOCATE(xccc1d)
   ABI_DEALLOCATE(xcccrc)
 end if
 if (usexcnhat==0.and.dtset%usepaw==1) then     
   ABI_DEALLOCATE(rhor_)
 end if

 if(n1xccc/=0) then
   ABI_DEALLOCATE(vxc_core)
   ABI_DEALLOCATE(vxc10_core)
   ABI_DEALLOCATE(vxc1is_core)
 end if

 if(fgga==1) then
   ABI_DEALLOCATE(rho0_redgr)
   ABI_DEALLOCATE(dgm)
   ABI_DEALLOCATE(d2gm)
 end if

end subroutine dfpt_eltfrxc
!!***
